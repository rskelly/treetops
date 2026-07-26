# Treetops

Treetops is an improvement upon a treetop-delineation algorithm developed by [Rob Skelly](rob@dijital.ca) at the University of Victoria's Hyperspectral and LiDAR Research Group, led by Dr. K. Olaf Niemann, funded jointly by [Terra Remote Sensing Inc.](https://www.terraremote.com/) and the University of Victoria.

Treetops detects tree tops and delineates tree crowns in canopy height models (CHMs/DSMs) derived from aerial LiDAR or photogrammetry.

[Project wiki](https://github.com/rskelly/treetops/wiki)

## Build, package, and run

### 1. Prerequisites

#### Linux (Debian/Ubuntu)

```bash
sudo apt-get install -y \
  build-essential cmake \
  libgdal-dev libgeos-dev \
  libsqlite3-dev libspatialite-dev \
  nlohmann-json3-dev \
  libwebkit2gtk-4.1-dev curl wget file \
  libxdo-dev libssl-dev libayatana-appindicator3-dev librsvg2-dev
```

#### Windows

Install a recent Visual Studio Build Tools setup with C++ support and a working Node.js installation. The Windows packaging flow expects PowerShell and the standard Windows runtime toolchain.

### 2. Build the command-line application

From the repository root:

```bash
packaging/linux/build.sh
```

This produces the CLI binary at `build/bin/treetops-cli`.

You can also build it manually:

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j"$(nproc)"
```

### 3. Run the command-line application

Show the CLI help:

```bash
./build/bin/treetops-cli -h
```

Show the version:

```bash
./build/bin/treetops-cli --version
```

Run with a configuration file:

```bash
./build/bin/treetops-cli -c _data/settings.json
```

Run with an explicit raster input and output directory:

```bash
./build/bin/treetops-cli -i _data/J5_10cm_CHM.tif -o _data
```

### 4. Run the desktop GUI

The desktop UI is a Tauri application that launches the CLI as a sidecar.

Install the web/frontend dependencies and start it in development mode:

```bash
cd app
npm install
npm run dev
```

Tauri requires Rust 1.88+; the repo pins this in `app/rust-toolchain.toml`. If you see an old-toolchain error, run:

```bash
rustup toolchain install 1.88.0
cd app && rustup override set 1.88.0
```

### 5. Package release artifacts

#### Linux packages

Create a Linux desktop bundle (AppImage, DEB, and RPM):

```bash
./packaging/linux/build-gui.sh
```

Artifacts are written to:

- `app/src-tauri/target/release/bundle/appimage/`
- `app/src-tauri/target/release/bundle/deb/`
- `app/src-tauri/target/release/bundle/rpm/`

#### Windows portable build

Build a portable Windows folder:

```powershell
packaging\windows\build-portable.ps1
```

This writes a portable package to `dist/treetops-portable-win64` and a zip archive under `dist/`.

#### Windows installer

Build a Windows installer (MSI-style packaging flow):

```powershell
packaging\windows\build-installer.ps1
```

The script prepares the Windows build and emits the packaged outputs to `dist/`.

### 6. Install and launch packaged builds

#### Linux

- AppImage: run `./Treetops_0.1.0_amd64.AppImage`
- DEB: install with `sudo apt install ./Treetops_0.1.0_amd64.deb`
- RPM: install with `sudo rpm -i ./Treetops-0.1.0-1.x86_64.rpm`

#### Windows

- Extract the portable folder from the zip archive or run the installer produced by the Windows packaging script.
- The GUI can be launched from the installed Start Menu entry or from the packaged application folder.

### CLI options

```text
  -c, --config FILE     JSON settings file
  -i, --input FILE      Input CHM/DSM raster (overrides config)
  -o, --output-dir DIR  Output directory (derives output paths from input name)
      --no-smooth        Skip Gaussian smoothing
      --no-crowns        Detect treetops only
  -h, --help            Show help
      --version         Show the application version
```
