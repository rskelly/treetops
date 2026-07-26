# Treetops

Treetops is an improvement upon a treetop-delineation algorithm developed by [Rob Skelly](rob@dijital.ca) at the University of Victoria's Hyperspectral and LiDAR Research Group, led by Dr. K. Olaf Niemann, funded jointly by [Terra Remote Sensing Inc.](https://www.terraremote.com/) and the University of Victoria.

Treetops detects tree tops and delineates tree crowns in canopy height models (CHMs/DSMs) derived from aerial LiDAR or photogrammetry.

[Project wiki](https://github.com/rskelly/treetops/wiki)

## Building

### Linux

Install dependencies (Debian/Ubuntu):

```bash
sudo apt-get install -y \
  build-essential cmake \
  libgdal-dev libgeos-dev \
  libsqlite3-dev libspatialite-dev \
  nlohmann-json3-dev
```

Build:

```bash
packaging/linux/build.sh
```

The binary is written to `build/bin/treetops-cli`.

You can also build manually:

```bash
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j"$(nproc)"
```

### Windows (portable package)

On Windows, build a self-contained folder that runs without installation:

```powershell
packaging\windows\build-portable.ps1
```

This produces:

- `dist\treetops-portable-win64\` — portable application folder
- `dist\treetops-portable-win64.zip` — distributable archive

The build script uses **vcpkg** by default (bootstrapped into `build\vcpkg` on first run). If **OSGeo4W** is installed at `C:\OSGeo4W`, it is used automatically. Override with `-OsGeo4WRoot` or `-VcpkgRoot` if needed.

The portable folder includes the executable, required DLLs, VC++ runtime DLLs, GDAL/PROJ data files, and a launcher batch file. No registry changes or system install is required.

Pre-built Windows packages are also available from GitHub Actions workflow artifacts (`treetops-portable-win64`).

### Desktop UI (Tauri)

For a packaged desktop build, use the helper scripts:

```bash
./packaging/linux/build-gui.sh
```

On Windows, build an installer with:

```powershell
packaging\windows\build-installer.ps1
```

These scripts build the CLI, prepare the sidecar, and invoke the desktop packaging workflow for the selected platform.

The Tauri app in `app/` provides settings configuration and processing controls with a status bar.

Linux prerequisites:

```bash
sudo apt-get install -y \
  libwebkit2gtk-4.1-dev build-essential curl wget file \
  libxdo-dev libssl-dev libayatana-appindicator3-dev librsvg2-dev
```

Build the CLI, then run the UI:

```bash
packaging/linux/build.sh
cd app
npm install
npm run dev
```

Tauri requires **Rust 1.88+** for current dependencies. The `app/rust-toolchain.toml` pins `1.88.0`; `npm run dev` runs `check-rust` to install it via rustup if needed. If you see a `dlopen2` / `edition2024` error, your Rust toolchain is too old — update with:

```bash
rustup toolchain install 1.88.0
cd app && rustup override set 1.88.0
```

The UI saves settings to JSON, launches `treetops-cli` as a sidecar, and displays `@status:` progress updates in the status bar.

Release artifacts are written to:
- Linux: `app/src-tauri/target/release/bundle/`
- Windows: `dist/` plus an MSI installer when the Windows packaging script completes.

## Usage

### Linux

```bash
./build/bin/treetops-cli -h
./build/bin/treetops-cli -c _data/settings.json
./build/bin/treetops-cli -i _data/J5_10cm_CHM.tif -o _data
```

### Windows (portable)

Unzip `dist\treetops-portable-win64.zip` anywhere, then run:

```bat
treetops.bat -h
treetops.bat -c examples\settings.json
treetops.bat -i C:\data\chm.tif -o C:\data\output
```

`treetops.bat` sets `PATH`, `GDAL_DATA`, and `PROJ_LIB` relative to the application folder. You can also run `treetops-cli.exe` directly if those environment variables are already configured.

### Options

```
  -c, --config FILE     JSON settings file
  -i, --input FILE      Input CHM/DSM raster (overrides config)
  -o, --output-dir DIR  Output directory (derives output paths from input name)
      --no-smooth        Skip Gaussian smoothing
      --no-crowns        Detect treetops only
  -h, --help            Show help
```
