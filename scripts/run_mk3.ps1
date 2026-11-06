param(
    [string]$DataDir = (Join-Path $PSScriptRoot '..\data'),
    [string]$OutputDir = '',
    [switch]$AllowPartial
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
Push-Location $projectRoot
try {
    $step = Join-Path $DataDir 'MK3_watertight.STEP'
    $iso = Join-Path $DataDir 'M08_ALT30_A00_B0_base_Iso_5spe_2_redo'
    $adiab = Join-Path $DataDir 'M08_ALT30_A00_B0_base_Adiab_5spe_9'
    foreach ($file in @($step, $iso, $adiab)) {
        if (-not (Test-Path -LiteralPath $file -PathType Leaf)) {
            throw "Missing input: $file. Put the three original files in data, or use the manual command in docs/RECONSTRUCTION.md for other filenames."
        }
    }
    $python = Join-Path $projectRoot '.venv-cad\Scripts\python.exe'
    if (-not (Test-Path -LiteralPath $python)) {
        if (-not (Get-Command py -ErrorAction SilentlyContinue)) {
            throw 'Install Python 3.12 with the Windows Python launcher, then run this script again.'
        }
        & py -3.12 -m venv .venv-cad
        if ($LASTEXITCODE -ne 0) { throw 'Python 3.12 is required for this CAD environment. Install it and retry.' }
    }
    & $python -m pip install -r requirements-cad.txt
    if ($LASTEXITCODE -ne 0) { throw 'CAD dependency installation failed; review the pip output above.' }
    if (-not $OutputDir) {
        $OutputDir = Join-Path 'outputs' ('mk3_' + (Get-Date -Format 'yyyyMMdd_HHmmss_fff'))
    }
    $arguments = @('-m', 'paneltemps', 'reconstruct', '--step', $step,
        '--bc', $iso, $adiab, '--transform', 'examples/mk3_alignment.json',
        '--coordinate-unit', 'm', '--temperature-unit', 'K',
        '--temperature-kinds', 'prescribed', 'adiabatic', '--output-dir', $OutputDir)
    if ($AllowPartial) { $arguments += '--allow-partial' }
    & $python @arguments
    $result = $LASTEXITCODE
    Write-Host "Output directory: $OutputDir"
    if ($result -eq 3) {
        Write-Host 'Reconstruction finished with unresolved panels or source points. Read reconstruction_report.json; unresolved averages remain blank.'
    }
    elseif ($result -ne 0) { throw "Reconstruction failed with exit code $result." }
    exit $result
}
finally { Pop-Location }
