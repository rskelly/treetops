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

if [[ "${OSTYPE:-}${OS:-}" == msys* || "${OSTYPE:-}${OS:-}" == cygwin* || "${OS:-}" == "Windows_NT" ]]; then
  DEST="$TARGET_DIR/treetops-cli.exe"
  ALT_DEST="$TARGET_DIR/treetops-cli-$TRIPLE.exe"
else
  DEST="$TARGET_DIR/treetops-cli"
  ALT_DEST="$TARGET_DIR/treetops-cli-$TRIPLE"
fi

cp "$SRC" "$DEST"
cp "$SRC" "$ALT_DEST"
if [[ "$DEST" != *.exe ]]; then
  chmod +x "$DEST" "$ALT_DEST"
fi

echo "Sidecar ready: $DEST"
echo "Sidecar ready: $ALT_DEST"
