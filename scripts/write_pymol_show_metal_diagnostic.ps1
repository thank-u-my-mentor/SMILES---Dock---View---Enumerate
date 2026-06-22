param(
    [string]$Out = "E:/TJ/pymol/diagnose_show_metal.pml",
    [string]$Structure = "E:/TJ/md_handoff/gromacs/md_centered_compact.gro"
)

$text = @"
reinitialize
run E:/Codex/scripts/pymol_show_metal_override.py
load $Structure, complex
hide everything
show cartoon, polymer.protein
color gray70, polymer.protein
show sticks, resn LIG
util.cbag resn LIG
show_metal

python
from pymol import cmd
print("[diagnose] show_metal command exists:", "show_metal" in cmd.keyword)
print("[diagnose] ple_show_metal command exists:", "ple_show_metal" in cmd.keyword)
python end
"@

$outPath = [System.IO.Path]::GetFullPath($Out)
New-Item -ItemType Directory -Force -Path ([System.IO.Path]::GetDirectoryName($outPath)) | Out-Null
Set-Content -LiteralPath $outPath -Value $text -Encoding UTF8
Write-Host $outPath
