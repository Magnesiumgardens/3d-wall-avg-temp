"""CAD reconstruction checks against known fields and finite CAD geometry."""

import json
from pathlib import Path
import tempfile
import unittest
import csv

try:
    import numpy as np
    from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder
    from OCP.BRep import BRep_Builder
    from OCP.TopoDS import TopoDS_Compound
    from OCP.gp import gp_Pnt
    from OCP.STEPControl import STEPControl_Writer, STEPControl_AsIs
    from paneltemps.cad import load_cad, load_transform, coordinate_hash, sha256, refine_triangles
    from paneltemps.reconstruction import FaceInterpolator, integrate, map_samples, unique_samples
    from paneltemps.reconstruct_cli import panels_for_faces
    from paneltemps.reconstruct_cli import main as reconstruct_main
except ImportError as exc:
    raise unittest.SkipTest('CAD checks require requirements-cad.txt: '+str(exc))

from paneltemps.model import Node


class CADReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.step = cls.root/'cube.step'
        writer = STEPControl_Writer()
        # STEP writer defaults to mm. A 1000 mm cube must import as a 1 m cube.
        writer.Transfer(BRepPrimAPI_MakeBox(1000., 1000., 1000.).Shape(), STEPControl_AsIs)
        writer.Write(str(cls.step))
        cls.faces = load_cad(cls.step)
        grid = np.array([0, .05, .13, .29, .5, .75, 1.0])
        points = set()
        for axis in range(3):
            for value in (0., 1.):
                for a in grid:
                    for b in grid:
                        point = [a, b]
                        point.insert(axis, value)
                        points.add(tuple(point))
        cls.points = np.array(sorted(points))
        cls.temperatures = np.array([300+cls.points @ np.array([10.,20.,30.]), np.full(len(points),300.)])
        cls.mappings, cls.membership, cls.primary, cls.distances = map_samples(cls.points, cls.faces, 1e-6)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_step_unit_conversion_and_exact_area(self):
        self.assertEqual(len(self.faces),6)
        self.assertAlmostEqual(sum(f.area for f in self.faces),6,places=7)
        for f in self.faces:
            self.assertTrue(np.all(f.bbox >= -1e-6))
            self.assertTrue(np.all(f.bbox <= 1+1e-6))

    def test_shared_edges_reuse_nodes_without_losing_faces(self):
        self.assertTrue((self.primary>0).all())
        self.assertLess(self.distances.max(),1e-8)
        self.assertTrue(any(len(m)==3 for m in self.membership))
        self.assertTrue(all(len(self.mappings[f.index])>=len(self.points)//8 for f in self.faces))

    def test_linear_field_surface_mean_std_and_uniform_case(self):
        for face in self.faces:
            with self.subTest(face=face.index):
                interpolation = FaceInterpolator(face,self.mappings[face.index],self.points,self.temperatures)
                result = integrate(face,interpolation,self.temperatures,.16,2.,.001)
                self.assertAlmostEqual(result['area'],1,places=7)
                self.assertAlmostEqual(result['coverage'],1,places=7)
                expected_mean = 300+np.dot(face.center,[10.,20.,30.])
                fixed_axis = np.argmin(face.bbox[1]-face.bbox[0])
                expected_std = np.sqrt(sum(c*c for i,c in enumerate([10.,20.,30.]) if i!=fixed_axis)/12)
                self.assertAlmostEqual(result['cases'][0]['mean'],expected_mean,places=7)
                self.assertAlmostEqual(result['cases'][0]['std'],expected_std,places=7)
                self.assertAlmostEqual(result['cases'][1]['mean'],300,places=7)
                self.assertLess(result['cases'][1]['std'],1e-8)

    def test_nonuniform_source_density_does_not_bias_linear_area_average(self):
        face=self.faces[0]
        interp=FaceInterpolator(face,self.mappings[face.index],self.points,self.temperatures)
        a=integrate(face,interp,self.temperatures,.25,2.,.001)
        b=integrate(face,interp,self.temperatures,.12,2.,.001)
        self.assertAlmostEqual(a['cases'][0]['mean'],b['cases'][0]['mean'],places=7)
        self.assertAlmostEqual(a['cases'][0]['std'],b['cases'][0]['std'],places=7)

    def test_distant_or_other_surface_samples_do_not_supply_a_face(self):
        shifted=self.points+10
        mapping,membership,primary,distances=map_samples(shifted,self.faces,1e-5)
        self.assertTrue((primary==0).all())
        self.assertTrue(all(not v for v in membership))
        self.assertTrue(all(not v for v in mapping.values()))

    def test_equidistant_unrelated_surfaces_are_ambiguous_not_arbitrarily_chosen(self):
        compound=TopoDS_Compound()
        builder=BRep_Builder();builder.MakeCompound(compound)
        builder.Add(compound,BRepPrimAPI_MakeBox(gp_Pnt(0,0,-10),1000,1000,10).Shape())
        builder.Add(compound,BRepPrimAPI_MakeBox(gp_Pnt(0,0,.2),1000,1000,10).Shape())
        step=self.root/'two_skins.step'
        writer=STEPControl_Writer();writer.Transfer(compound,STEPControl_AsIs);writer.Write(str(step))
        faces=load_cad(step)
        samples=np.array([[.5,.5,.0001],[.5,.5,.00001]])
        _,membership,primary,distances=map_samples(samples,faces,.0005)
        self.assertEqual(primary[0],0)
        self.assertFalse(membership[0])
        self.assertTrue(np.isfinite(distances[0]))
        self.assertNotEqual(primary[1],0)
        self.assertEqual(len(membership[1]),1)

    def test_missing_and_insufficient_source_support_are_not_fabricated(self):
        face=self.faces[0]
        with self.assertRaises(ValueError):
            FaceInterpolator(face,self.mappings[face.index][:2],self.points,self.temperatures)
        result=integrate(face,None,self.temperatures,.2,.02,.001)
        self.assertEqual(result['coverage'],0)
        self.assertEqual(result['cases'],[None,None])

    def test_distance_limit_leaves_large_gaps_uncovered(self):
        face=self.faces[0]
        interp=FaceInterpolator(face,self.mappings[face.index],self.points,self.temperatures)
        result=integrate(face,interp,self.temperatures,.2,.001,.0001)
        self.assertLess(result['coverage'],.1)

    def test_duplicate_coordinates_require_agreement(self):
        a=Node(1,0,0,0,1,300)
        b=Node(2,0,0,0,2,300)
        points,temperatures,ids,duplicates=unique_samples({1:a,2:b})
        self.assertEqual(len(points),1)
        self.assertEqual(ids,[[1,2]])
        self.assertEqual(temperatures.tolist(),[300])
        with self.assertRaises(ValueError):
            unique_samples({1:a,2:Node(2,0,0,0,2,301)})

    def test_projection_aliases_use_geometry_and_reject_conflicting_ties(self):
        class Plane:
            index=100
            def evaluate(self,p):return np.array([*p,0.])
            def inside(self,p):return True
        face=Plane()
        points=np.array([[0,0,.0001],[0,0,.0002],[1,0,0],[0,1,0]],dtype=float)
        mapping=[(0,[0,0]),(1,[0,0]),(2,[1,0]),(3,[0,1])]
        temps=np.array([[300,400,310,320]],dtype=float)
        interp=FaceInterpolator(face,mapping,points,temps)
        self.assertIn(0,interp.indices)
        self.assertNotIn(1,interp.indices)
        self.assertEqual(interp.projected_collisions[0]['selected_sample_index'],0)
        points[1,2]=-.0001
        with self.assertRaises(ValueError):FaceInterpolator(face,mapping,points,temps)

    def test_curved_surface_jacobian_integrates_cylinder_side_area(self):
        step=self.root/'cylinder.step'
        writer=STEPControl_Writer()
        writer.Transfer(BRepPrimAPI_MakeCylinder(100.,200.).Shape(),STEPControl_AsIs)
        writer.Write(str(step))
        faces=load_cad(step)
        side=max(faces,key=lambda f:f.area)
        result=integrate(side,None,np.zeros((1,0)),.04,.1,.001)
        self.assertAlmostEqual(result['area'],2*np.pi*.1*.2,places=7)

    def test_transforms_reject_scale_reflection_and_wrong_geometry(self):
        path=self.root/'transform.json'
        base={'rotation_BC_to_CAD':np.eye(3).tolist(),'translation_m':[0,0,0],
              'cad_sha256':sha256(self.step),'source_coordinates_sha256':coordinate_hash(self.points)}
        path.write_text(json.dumps(base))
        points,_=load_transform(path,self.step,self.points)
        np.testing.assert_allclose(points,self.points)
        for rotation in (np.eye(3)*2,np.diag([-1,1,1])):
            path.write_text(json.dumps(dict(base,rotation_BC_to_CAD=rotation.tolist())))
            with self.assertRaises(ValueError):load_transform(path,self.step,self.points)
        path.write_text(json.dumps(base))
        with self.assertRaises(ValueError):load_transform(path,self.step,self.points+1)
        path.write_text(json.dumps(dict(base,cad_sha256='wrong')))
        with self.assertRaises(ValueError):load_transform(path,self.step,self.points)

    def test_face_groups_reject_overlaps_and_silent_omissions(self):
        groups=self.root/'groups.json'
        base={'cad_sha256':sha256(self.step),'panels':[{'id':'a','name':'a','face_indices':[1,2]}]}
        groups.write_text(json.dumps(base))
        with self.assertRaises(ValueError):panels_for_faces(groups,self.faces,sha256(self.step))
        base['excluded_face_indices']=[3,4,5,6]
        groups.write_text(json.dumps(base))
        panels,excluded=panels_for_faces(groups,self.faces,sha256(self.step))
        self.assertEqual(excluded,[3,4,5,6])
        base['panels'].append({'id':'b','name':'b','face_indices':[2]})
        groups.write_text(json.dumps(base))
        with self.assertRaises(ValueError):panels_for_faces(groups,self.faces,sha256(self.step))

    def test_trimmed_hole_is_not_bridged_by_a_linear_stencil(self):
        class TrimmedFace:
            index=99
            def inside(self,p):
                return not (abs(p[0])<.2 and abs(p[1])<.2)
            def evaluate(self,p):
                return np.array([*p,0.])
        face=TrimmedFace()
        uv=np.array([[-1,-1],[-1,1],[1,-1],[1,1]],dtype=float)
        points=np.column_stack([uv,np.zeros(4)])
        interp=FaceInterpolator(face,list(enumerate(uv)),points,np.full((1,4),300.))
        stencil=interp.stencil(np.array([[-.5,-.5]]),np.array([[-.5,-.5,0.]]),3.,.01)
        self.assertFalse(stencil.valid[0])

    def test_skinny_CAD_triangles_refine_without_exploding_in_size(self):
        class Plane:
            def evaluate(self,p):return np.array([*p,0.])
        triangles=refine_triangles(Plane(),np.array([[[0.,0.],[1.,0.],[0.,.0001]]]),.01)
        # A fixed resource budget catches the previous four-way subdivision,
        # while permitting bisection to resolve both long edges of the sliver.
        self.assertLess(len(triangles),4096)
        lengths=np.linalg.norm(triangles-np.roll(triangles,1,axis=1),axis=2)
        self.assertLessEqual(lengths.max(),.01)
        a,b=triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]
        self.assertAlmostEqual(np.abs(a[:,0]*b[:,1]-a[:,1]*b[:,0]).sum()/2,.00005,places=10)

    def test_cli_exports_two_separate_known_cases_and_reports_partial_failure(self):
        case_a,case_b=self.root/'linear.csv',self.root/'fixed.csv'
        header=['nodenumber','x-coordinate','y-coordinate','z-coordinate','pressure','temperature','wall-temperature']
        for path,values in ((case_a,self.temperatures[0]),(case_b,self.temperatures[1])):
            with path.open('w',newline='') as stream:
                writer=csv.writer(stream);writer.writerow(header)
                writer.writerows([[i+1,*p,100,280,t] for i,(p,t) in enumerate(zip(self.points,values))])
        transform=self.root/'identity.json'
        transform.write_text(json.dumps({'rotation_BC_to_CAD':np.eye(3).tolist(),'translation_m':[0,0,0]}))
        args=['--step',str(self.step),'--bc',str(case_a),str(case_b),'--transform',str(transform),
              '--edge-m','.25','--max-support-m','2']
        out=self.root/'cli_complete'
        self.assertEqual(reconstruct_main([*args,'--output-dir',str(out)]),0)
        report=json.loads((out/'reconstruction_report.json').read_text())
        self.assertTrue(report['all_numerical_checks_passed'])
        self.assertFalse(report['cad_selections_verified'])
        self.assertEqual(len(report['cases']),2)
        self.assertEqual((out/'linear.csv'/'original_boundary_conditions.txt').read_bytes(),case_a.read_bytes())
        with (out/'fixed.csv'/'panel_temperatures.csv').open() as stream:
            rows=list(csv.DictReader(stream))
        self.assertEqual(len(rows),6)
        self.assertTrue(all(abs(float(r['temperature_avg_K'])-300)<1e-7 for r in rows))
        partial=self.root/'cli_partial'
        self.assertEqual(reconstruct_main([*args,'--max-support-m','.001','--extrapolation-m','.0001',
                                          '--output-dir',str(partial)]),3)
        with (partial/'linear.csv'/'panel_temperatures.csv').open() as stream:
            rows=list(csv.DictReader(stream))
        self.assertTrue(all(r['status']=='UNRESOLVED' and not r['temperature_avg_K'] for r in rows))


if __name__=='__main__':
    unittest.main()
