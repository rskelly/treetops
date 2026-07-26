#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

bash "$ROOT/packaging/linux/build.sh"
cd "$ROOT/app"
npm install
npm run build

echo "Linux desktop bundle generated in app/src-tauri/target/release/bundle"
