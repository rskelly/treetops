#Requires -Version 5.1
[CmdletBinding()]
param(
    [string]$BuildDir = "",
    [string]$OutputDir = ""
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Ensure-Directory([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) {
        New-Item -ItemType Directory -Path $Path -Force | Out-Null
    }
}

$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
if (-not $BuildDir) { $BuildDir = Join-Path $Root "build\windows-installer" }
if (-not $OutputDir) { $OutputDir = Join-Path $Root "dist" }
Ensure-Directory $BuildDir
Ensure-Directory $OutputDir

$repo = $Root
$build = $BuildDir
$dist = $OutputDir

Write-Host "==> Building CLI for Windows packaging"
& "$repo\packaging\windows\build-portable.ps1" -BuildDir "$build\portable" -OutputDir "$dist\treetops-portable-win64"

if (-not (Get-Command npm -ErrorAction SilentlyContinue)) {
    throw "npm was not found on PATH. Install Node.js first."
}

Push-Location "$repo\app"
try {
    npm install
    npm run build
} finally {
    Pop-Location
}

$bundleDir = Join-Path $repo "app\src-tauri\target\release\bundle"
if (-not (Test-Path -LiteralPath $bundleDir)) {
    throw "No Tauri bundle output found in $bundleDir"
}

Write-Host "==> Windows packaging workflow completed"
Write-Host "Artifacts are available under $dist"
