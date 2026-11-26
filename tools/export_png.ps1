# Export every slide of a pptx to PNG with PowerPoint (COM).
# Usage: pwsh -File export_png.ps1 <pptx> <out folder> [width]
param([string]$Pptx, [string]$OutDir, [int]$Width = 1920)
$ErrorActionPreference = "Stop"
$Pptx = (Resolve-Path $Pptx).Path
New-Item -ItemType Directory -Force $OutDir | Out-Null
$OutDir = (Resolve-Path $OutDir).Path
$app = New-Object -ComObject PowerPoint.Application
try {
    # ReadOnly, Untitled, WithWindow = false
    $pres = $app.Presentations.Open($Pptx, -1, 0, 0)
    $h = [int]($Width * $pres.PageSetup.SlideHeight / $pres.PageSetup.SlideWidth)
    $base = [IO.Path]::GetFileNameWithoutExtension($Pptx)
    foreach ($s in $pres.Slides) {
        $f = Join-Path $OutDir ("{0}_{1}.png" -f $base, $s.SlideIndex)
        $s.Export($f, "PNG", $Width, $h)
        Write-Output $f
    }
    $pres.Close()
}
finally {
    if ($app.Presentations.Count -eq 0) { $app.Quit() }
    [void][Runtime.InteropServices.Marshal]::ReleaseComObject($app)
}
