#Requires -Version 5.1
<#
.SYNOPSIS
    Build a portable Windows distribution of treetops-cli.

.DESCRIPTION
    Produces a self-contained folder (and zip) that can run without installation.
    Dependencies are resolved via vcpkg manifest mode or an existing OSGeo4W install.

.PARAMETER VcpkgRoot
    Path to a bootstrapped vcpkg checkout. When omitted, the script clones vcpkg
    into build/vcpkg if no OSGeo4W installation is found.

.PARAMETER OsGeo4WRoot
    Optional OSGeo4W root (for example C:\OSGeo4W). When set, CMake uses OSGeo4W
    packages instead of vcpkg.

.PARAMETER BuildDir
    CMake build directory. Defaults to build/windows-portable.

.PARAMETER OutputDir
    Staging directory for the portable package.

.PARAMETER Triplet
    vcpkg triplet. Defaults to x64-windows.
#>
[CmdletBinding()]
param(
    [string]$VcpkgRoot = "",
    [string]$OsGeo4WRoot = "",
    [string]$BuildDir = "",
    [string]$OutputDir = "",
    [string]$Triplet = "x64-windows"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Write-Step([string]$Message) {
    Write-Host "==> $Message"
}

function Ensure-Directory([string]$Path) {
    if (-not (Test-Path -LiteralPath $Path)) {
        New-Item -ItemType Directory -Path $Path | Out-Null
    }
}

function Copy-Tree([string]$Source, [string]$Destination) {
    if (-not (Test-Path -LiteralPath $Source)) {
        return
    }
    Ensure-Directory $Destination
    Copy-Item -LiteralPath $Source -Destination $Destination -Recurse -Force
}

function Copy-VcRuntime([string]$Destination) {
    $system32 = Join-Path $env:WINDIR "System32"
    foreach ($dll in @("msvcp140.dll", "vcruntime140.dll", "vcruntime140_1.dll")) {
        $source = Join-Path $system32 $dll
        if (Test-Path -LiteralPath $source) {
            Copy-Item -LiteralPath $source -Destination $Destination -Force
        }
    }
}

function Get-DumpbinPath {
    $vswhere = Join-Path ${env:ProgramFiles(x86)} "Microsoft Visual Studio\Installer\vswhere.exe"
    if (-not (Test-Path -LiteralPath $vswhere)) {
        return $null
    }

    $installationPath = & $vswhere -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $installationPath) {
        return $null
    }

    $dumpbin = Get-ChildItem -Path (Join-Path $installationPath "VC\Tools\MSVC") -Filter dumpbin.exe -Recurse -ErrorAction SilentlyContinue |
        Sort-Object FullName -Descending |
        Select-Object -First 1

    if ($dumpbin) {
        return $dumpbin.FullName
    }
    return $null
}

function Add-DllDependencies {
    param(
        [string]$BinaryPath,
        [string]$SearchDir,
        [hashtable]$Seen,
        [string]$DumpbinPath
    )

    if (-not (Test-Path -LiteralPath $BinaryPath)) {
        return
    }

    $key = (Resolve-Path -LiteralPath $BinaryPath).Path.ToLowerInvariant()
    if ($Seen.ContainsKey($key)) {
        return
    }
    $Seen[$key] = $true

    $output = & $DumpbinPath /nologo /dependents $BinaryPath
    foreach ($line in $output) {
        $dllName = $line.Trim()
        if ($dllName -notmatch '\.dll$') {
            continue
        }
        if ($dllName -match '^(KERNEL32|USER32|ADVAPI32|OLE32|SHELL32|WS2_32|MSVCRT|VCRUNTIME|UCRTBASE)\.DLL$') {
            continue
        }

        $candidate = Join-Path $SearchDir $dllName
        if (Test-Path -LiteralPath $candidate) {
            Add-DllDependencies -BinaryPath $candidate -SearchDir $SearchDir -Seen $Seen -DumpbinPath $DumpbinPath
        }
    }
}

function Copy-DependenciesFromDumpbin {
    param(
        [string]$ExePath,
        [string]$SearchDir,
        [string]$Destination
    )

    $dumpbin = Get-DumpbinPath
    if (-not $dumpbin) {
        Write-Warning "dumpbin not found; copying all DLLs from $SearchDir"
        Get-ChildItem -Path $SearchDir -Filter *.dll | Copy-Item -Destination $Destination -Force
        return
    }

    $seen = @{}
    Add-DllDependencies -BinaryPath $ExePath -SearchDir $SearchDir -Seen $seen -DumpbinPath $dumpbin

    foreach ($binary in $seen.Keys) {
        Copy-Item -LiteralPath $binary -Destination $Destination -Force
    }
}

$Root = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
if (-not $BuildDir) {
    $BuildDir = Join-Path $Root "build\windows-portable"
}
if (-not $OutputDir) {
    $OutputDir = Join-Path $Root "dist\treetops-portable-win64"
}

$useOsGeo = $false
if ($OsGeo4WRoot) {
    if (-not (Test-Path -LiteralPath $OsGeo4WRoot)) {
        throw "OSGeo4W root not found: $OsGeo4WRoot"
    }
    $useOsGeo = $true
} elseif (Test-Path -LiteralPath "C:\OSGeo4W") {
    $OsGeo4WRoot = "C:\OSGeo4W"
    $useOsGeo = $true
}

$toolchainFile = ""
$depSearchDir = ""
$gdalDataSource = ""
$projDataSource = ""

if ($useOsGeo) {
    Write-Step "Building with OSGeo4W at $OsGeo4WRoot"
    $env:PATH = (Join-Path $OsGeo4WRoot "bin") + ";" + $env:PATH
    $depSearchDir = Join-Path $OsGeo4WRoot "bin"
    $gdalDataSource = Join-Path $OsGeo4WRoot "share\gdal"
    $projDataSource = Join-Path $OsGeo4WRoot "share\proj"
} else {
    if (-not $VcpkgRoot) {
        $VcpkgRoot = Join-Path $Root "build\vcpkg"
    }

    if (-not (Test-Path -LiteralPath (Join-Path $VcpkgRoot "vcpkg.exe"))) {
        Write-Step "Bootstrapping vcpkg in $VcpkgRoot"
        Ensure-Directory (Split-Path -Parent $VcpkgRoot)
        if (-not (Test-Path -LiteralPath $VcpkgRoot)) {
            git clone https://github.com/microsoft/vcpkg.git $VcpkgRoot
        }
        Push-Location $VcpkgRoot
        try {
            & .\bootstrap-vcpkg.bat -disableMetrics
        } finally {
            Pop-Location
        }
    }

    $toolchainFile = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
    $installedRoot = Join-Path $VcpkgRoot "installed\$Triplet"
    $depSearchDir = Join-Path $installedRoot "bin"
    $gdalDataSource = Join-Path $installedRoot "share\gdal"
    $projDataSource = Join-Path $installedRoot "share\proj"
    Write-Step "Building with vcpkg at $VcpkgRoot ($Triplet)"
}

Write-Step "Configuring CMake"
Ensure-Directory $BuildDir

$cmakeArgs = @(
    "-S", $Root,
    "-B", $BuildDir,
    "-DCMAKE_BUILD_TYPE=Release"
)

if ($toolchainFile) {
    $cmakeArgs += @(
        "-DCMAKE_TOOLCHAIN_FILE=$toolchainFile",
        "-DVCPKG_TARGET_TRIPLET=$Triplet"
    )
} else {
    $cmakeArgs += "-DCMAKE_PREFIX_PATH=$OsGeo4WRoot"
}

cmake @cmakeArgs

Write-Step "Compiling treetops-cli"
cmake --build $BuildDir --config Release --target treetops-cli

$builtExe = Join-Path $BuildDir "bin\treetops-cli.exe"
if (-not (Test-Path -LiteralPath $builtExe)) {
    $builtExe = Join-Path $BuildDir "bin\Release\treetops-cli.exe"
}
if (-not (Test-Path -LiteralPath $builtExe)) {
    throw "Built executable not found in $BuildDir\bin"
}

Write-Step "Staging portable package at $OutputDir"
if (Test-Path -LiteralPath $OutputDir) {
    Remove-Item -LiteralPath $OutputDir -Recurse -Force
}
Ensure-Directory $OutputDir

Copy-Item -LiteralPath $builtExe -Destination (Join-Path $OutputDir "treetops-cli.exe") -Force
Copy-Item -LiteralPath (Join-Path $PSScriptRoot "treetops.bat") -Destination (Join-Path $OutputDir "treetops.bat") -Force

$existingDlls = Get-ChildItem -Path (Split-Path -Parent $builtExe) -Filter *.dll -ErrorAction SilentlyContinue
if ($existingDlls) {
    $existingDlls | Copy-Item -Destination $OutputDir -Force
}

Copy-DependenciesFromDumpbin -ExePath (Join-Path $OutputDir "treetops-cli.exe") -SearchDir $depSearchDir -Destination $OutputDir
Copy-VcRuntime -Destination $OutputDir
Copy-Tree -Source $gdalDataSource -Destination (Join-Path $OutputDir "share\gdal")
Copy-Tree -Source $projDataSource -Destination (Join-Path $OutputDir "share\proj")

$settingsExample = Join-Path $Root "_data\settings.json"
if (Test-Path -LiteralPath $settingsExample) {
    Ensure-Directory (Join-Path $OutputDir "examples")
    Copy-Item -LiteralPath $settingsExample -Destination (Join-Path $OutputDir "examples\settings.json") -Force
}

$zipPath = "$OutputDir.zip"
if (Test-Path -LiteralPath $zipPath) {
    Remove-Item -LiteralPath $zipPath -Force
}
Compress-Archive -Path (Join-Path $OutputDir "*") -DestinationPath $zipPath

Write-Step "Portable package ready"
Write-Host "Folder: $OutputDir"
Write-Host "Zip:    $zipPath"
Write-Host ""
Write-Host "Run without installation:"
Write-Host "  $OutputDir\treetops.bat -h"
Write-Host "  $OutputDir\treetops-cli.exe -h"
