#!/usr/bin/env bash
set -euo pipefail

MIN_VERSION="1.88.0"

version_ge() {
  local current="$1"
  local required="$2"
  [ "$(printf '%s\n%s\n' "$required" "$current" | sort -V | head -n1)" = "$required" ]
}

if ! command -v rustup >/dev/null 2>&1; then
  echo "rustup is required. Install from https://rustup.rs" >&2
  exit 1
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if ! rustup toolchain list | grep -q '^1\.88\.0'; then
  echo "Installing Rust $MIN_VERSION for Tauri..."
  rustup toolchain install "$MIN_VERSION" --profile minimal -c rustc -c cargo -c rust-std
fi

if ! rustup show active-toolchain | grep -q '1\.88\.0'; then
  rustup override set "$MIN_VERSION"
fi

CURRENT="$(rustc --version | sed -E 's/rustc ([0-9]+\.[0-9]+\.[0-9]+).*/\1/')"
if ! version_ge "$CURRENT" "$MIN_VERSION"; then
  echo "Rust $MIN_VERSION or newer is required (found $CURRENT)." >&2
  echo "Run: rustup update stable" >&2
  echo "Or:  rustup toolchain install $MIN_VERSION" >&2
  exit 1
fi

echo "Rust toolchain OK: $(rustc --version) ($(cargo --version))"
