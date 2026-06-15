#!/usr/bin/env bash
# Build treetops-cli on Linux using system GDAL/GEOS libraries.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD="${BUILD_DIR:-$ROOT/build}"
BUILD_TYPE="${BUILD_TYPE:-Release}"

mkdir -p "$BUILD"
cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
cmake --build "$BUILD" -j"$(nproc)"

echo "Built: $BUILD/bin/treetops-cli"
