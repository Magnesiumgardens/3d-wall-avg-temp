from dataclasses import dataclass, field


@dataclass(frozen=True)
class Node:
    node_id: int
    x: float
    y: float
    z: float
    pressure: float
    temperature: float

    @property
    def xyz(self):
        return self.x, self.y, self.z


@dataclass(frozen=True)
class Triangle:
    element_id: int
    node_ids: tuple[int, int, int]
    panel_id: str | None = None


@dataclass(frozen=True)
class Panel:
    panel_id: str
    name: str
    solidworks_component: str = ""
    solidworks_region: str = ""
    selection_notes: str = ""
    selector: dict = field(default_factory=dict)
