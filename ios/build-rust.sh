#!/usr/bin/env bash
# Build the real FrankenSim laboratory kernel for iPhone, simulator, and Mac
# Catalyst, then package it as one non-embedded static XCFramework.
set -euo pipefail

cd "$(dirname "$0")/.."

APPLE_DEPLOYMENT_TARGET="${APPLE_DEPLOYMENT_TARGET:-17.0}"
# Apple mobile and Catalyst triples require Xcode's SDK and linker, which are
# unavailable on the Linux RCH fleet. Resolve the pinned local Cargo binary
# explicitly so an installed Cargo wrapper cannot silently fall back remotely.
APPLE_RUST_TOOLCHAIN="${APPLE_RUST_TOOLCHAIN:-nightly-2026-08-25-aarch64-apple-darwin}"
APPLE_CARGO="${APPLE_CARGO:-$(rustup which --toolchain "$APPLE_RUST_TOOLCHAIN" cargo)}"
TARGET_BASE="${RCH_TARGET_BASE:-${CARGO_TARGET_DIR:-target}}"
TARGET_ROOT="${FRANKENSIM_APPLE_TARGET_DIR:-${TARGET_BASE}/frankensim-apple}"
MANIFEST="ios/rust/Cargo.toml"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "error: FrankenSim Apple slices require a Darwin host with Xcode" >&2
  exit 2
fi

# The nested Apple workspace owns an intentionally separate lockfile. Refresh
# it only when the caller opts in after changing path dependencies; ordinary
# release builds remain locked and therefore fail closed on dependency drift.
if [[ "${APPLE_REFRESH_LOCKFILE:-0}" == "1" ]]; then
  RUSTUP_TOOLCHAIN="$APPLE_RUST_TOOLCHAIN" \
    RCH_CARGO_WRAPPER_BYPASS=1 \
    "$APPLE_CARGO" generate-lockfile --offline --manifest-path "$MANIFEST"
fi

for target in \
  aarch64-apple-ios \
  aarch64-apple-ios-sim \
  aarch64-apple-ios-macabi \
  x86_64-apple-ios-macabi
do
  rustup target list --toolchain "$APPLE_RUST_TOOLCHAIN" --installed | grep -qx "$target" || \
    rustup target add --toolchain "$APPLE_RUST_TOOLCHAIN" "$target"
  IPHONEOS_DEPLOYMENT_TARGET="$APPLE_DEPLOYMENT_TARGET" \
    RUSTUP_TOOLCHAIN="$APPLE_RUST_TOOLCHAIN" \
    RCH_CARGO_WRAPPER_BYPASS=1 \
    CARGO_TARGET_DIR="$TARGET_ROOT" \
    "$APPLE_CARGO" build --release --locked --manifest-path "$MANIFEST" --target "$target"
done

HEADER_ROOT="$(mktemp -d /tmp/frankensim-apple-headers.XXXXXX)"
cp ios/rust/include/frankensim_apple.h ios/rust/include/module.modulemap "$HEADER_ROOT/"

CATALYST_ROOT="$(mktemp -d /tmp/frankensim-apple-catalyst.XXXXXX)"
lipo -create \
  "$TARGET_ROOT/aarch64-apple-ios-macabi/release/libfrankensim_apple.a" \
  "$TARGET_ROOT/x86_64-apple-ios-macabi/release/libfrankensim_apple.a" \
  -output "$CATALYST_ROOT/libfrankensim_apple.a"

OUTPUT_ROOT="$(mktemp -d /tmp/frankensim-apple-xcframework.XXXXXX)"
STAGED="$OUTPUT_ROOT/FrankenSimCore.xcframework"
xcodebuild -create-xcframework \
  -library "$TARGET_ROOT/aarch64-apple-ios/release/libfrankensim_apple.a" -headers "$HEADER_ROOT" \
  -library "$TARGET_ROOT/aarch64-apple-ios-sim/release/libfrankensim_apple.a" -headers "$HEADER_ROOT" \
  -library "$CATALYST_ROOT/libfrankensim_apple.a" -headers "$HEADER_ROOT" \
  -output "$STAGED"

FRAMEWORK="ios/FrankenSimCore.xcframework"
if [[ -e "$FRAMEWORK" ]]; then
  mv "$FRAMEWORK" "$FRAMEWORK.previous-$(date +%Y%m%d-%H%M%S)"
fi
mv "$STAGED" "$FRAMEWORK"
echo "built $FRAMEWORK"
