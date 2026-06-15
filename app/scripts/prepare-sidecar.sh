#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$ROOT/.." && pwd)"
SRC="$REPO/build/bin/treetops-cli"
TARGET_DIR="$ROOT/src-tauri/binaries"

mkdir -p "$TARGET_DIR"

if [[ ! -f "$SRC" ]]; then
  echo "Build treetops-cli first:" >&2
  echo "  packaging/linux/build.sh" >&2
  exit 1
fi

TRIPLE="$(rustc -vV | sed -n 's/^host: //p')"
DEST="$TARGET_DIR/treetops-cli-$TRIPLE"

cp "$SRC" "$DEST"
chmod +x "$DEST"
echo "Sidecar ready: $DEST"
