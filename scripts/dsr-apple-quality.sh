#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root/ios"

build_root="${FRANKEN_APPLE_BUILD_ROOT:-${DSR_QUALITY_RUN_DIR:-$repo_root/ios/build/dsr-apple-quality}}"
mkdir -p "$build_root"
sbh check --need 20G "$build_root"
command -v xcodegen >/dev/null
xcodegen generate --spec project.yml
git diff --exit-code -- FrankenSim.xcodeproj Sources/Info.plist
git ls-files -z -- '*.swift' | xargs -0 xcrun swiftc -parse
plutil -lint Sources/Info.plist
plutil -lint Sources/PrivacyInfo.xcprivacy
/Users/jemanuel/.local/bin/ensure-simulator-audio-safe prepare
xcodebuild -project FrankenSim.xcodeproj -scheme FrankenSim \
  -destination 'generic/platform=iOS Simulator' \
  -derivedDataPath "$build_root/derived-data" \
  CODE_SIGNING_ALLOWED=NO build
xcodebuild -project FrankenSim.xcodeproj -scheme FrankenSim \
  -destination 'platform=macOS,variant=Mac Catalyst' \
  -derivedDataPath "$build_root/derived-data" \
  CODE_SIGNING_ALLOWED=NO test -only-testing:FrankenSimTests

/Users/jemanuel/.local/bin/ensure-simulator-audio-safe prepare
iphone_udid="${FRANKENSIM_IPHONE_UDID:-$(
  xcrun simctl list devices available --json \
    | jq -r '
        [.devices[][] | select(.name | contains("iPhone"))] as $iphones
        | (($iphones | map(select(.name | test("FrankenSim"; "i"))))
            + ($iphones | map(select(.state == "Booted")))
            + $iphones)
        | .[0].udid // empty
      '
)}"
if [[ -z "$iphone_udid" ]]; then
  echo "No available iPhone Simulator found for the FrankenSim UI proof lane" >&2
  exit 1
fi

/Users/jemanuel/.local/bin/ensure-simulator-audio-safe prepare
xcodebuild -project FrankenSim.xcodeproj -scheme FrankenSim \
  -destination "platform=iOS Simulator,id=$iphone_udid" \
  -derivedDataPath "$build_root/derived-data" \
  -resultBundlePath "$build_root/frankensim-iphone-ui.xcresult" \
  -parallel-testing-enabled NO \
  CODE_SIGNING_ALLOWED=NO test \
  -only-testing:FrankenSimUITests/FrankenSimAppearanceUITests
