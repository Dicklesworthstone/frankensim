#!/usr/bin/env bash
set -eu

corpus_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="${corpus_root}/sources"
repository="https://raw.githubusercontent.com/adafruit/Adafruit_CAD_Parts"
commit="ab3dfc47c32468ba87e7652556cab25efd906eb0"

mkdir -p "${source_root}"

fetch_one() {
    destination_name="$1"
    upstream_path="$2"
    expected_blob="$3"
    destination="${source_root}/${destination_name}"
    encoded_path="${upstream_path// /%20}"
    url="${repository}/${commit}/${encoded_path}"

    if [ -e "${destination}" ]; then
        observed_blob="$(git hash-object "${destination}")"
        if [ "${observed_blob}" != "${expected_blob}" ]; then
            echo "REFUSED existing ${destination_name}: blob ${observed_blob}, expected ${expected_blob}" >&2
            return 1
        fi
        echo "verified existing ${destination_name}"
        return 0
    fi

    remote_blob="$(curl --fail --location --retry 3 --silent --show-error "${url}" | git hash-object --stdin)"
    if [ "${remote_blob}" != "${expected_blob}" ]; then
        echo "REFUSED remote ${destination_name}: blob ${remote_blob}, expected ${expected_blob}" >&2
        return 1
    fi

    curl --fail --location --retry 3 --silent --show-error --output "${destination}" "${url}"
    observed_blob="$(git hash-object "${destination}")"
    if [ "${observed_blob}" != "${expected_blob}" ]; then
        echo "REFUSED retained ${destination_name}: blob ${observed_blob}, expected ${expected_blob}" >&2
        return 1
    fi
    echo "retained ${destination_name}"
}

fetch_http_snapshot() {
    destination_name="$1"
    url="$2"
    expected_sha256="$3"
    destination="${source_root}/${destination_name}"

    if [ -e "${destination}" ]; then
        observed_sha256="$(shasum -a 256 "${destination}" | awk '{print $1}')"
        if [ "${observed_sha256}" != "${expected_sha256}" ]; then
            echo "REFUSED existing ${destination_name}: sha256 ${observed_sha256}, expected ${expected_sha256}" >&2
            return 1
        fi
        echo "verified existing ${destination_name}"
        return 0
    fi

    remote_sha256="$(curl --fail --location --retry 3 --silent --show-error "${url}" | shasum -a 256 | awk '{print $1}')"
    if [ "${remote_sha256}" != "${expected_sha256}" ]; then
        echo "REFUSED remote ${destination_name}: sha256 ${remote_sha256}, expected ${expected_sha256}" >&2
        return 1
    fi

    curl --fail --location --retry 3 --silent --show-error --output "${destination}" "${url}"
    observed_sha256="$(shasum -a 256 "${destination}" | awk '{print $1}')"
    if [ "${observed_sha256}" != "${expected_sha256}" ]; then
        echo "REFUSED retained ${destination_name}: sha256 ${observed_sha256}, expected ${expected_sha256}" >&2
        return 1
    fi
    echo "retained ${destination_name}"
}

fetch_one "adafruit-0258-1200mah-lipo-step.step" \
    "258 1200mAh lipo/258 1200mAh lipo.step" \
    "829f07a8493d53a840334b2cab9fcd3437609ace"
fetch_one "adafruit-0258-1200mah-lipo-stl.stl" \
    "258 1200mAh lipo/258 1200mAh lipo.stl" \
    "823579716c691ab4fa5f83ca990a57d0eefc952b"
fetch_one "adafruit-2278-rgb-matrix-4mm-step.step" \
    "2278 RGB Matrix 4mm/2278 RGB Matrix 4mm.step" \
    "509544807bc32d94f12dff0c87f7c77c1881f7ce"
fetch_one "adafruit-2278-rgb-matrix-4mm-stl.stl" \
    "2278 RGB Matrix 4mm/2278 RGB Matrix 4mm.stl" \
    "f4b398f266b191c0aa1d4b3c0c17248ee8708a44"
fetch_one "adafruit-2719-oled-24in-step.step" \
    "2719 OLED 2.4in/2719 OLED 2.4in.step" \
    "e8096437534dff6e41bb053f03ac42b30f9e7dde"
fetch_one "adafruit-2719-oled-24in-stl.stl" \
    "2719 OLED 2.4in/2719 OLED 2.4in.stl" \
    "76530164662e3ca11c4873e8ac24f53ef7237d62"
fetch_one "adafruit-3258-microusb-panel-mount-step.step" \
    "3258 USB Panel Mout Cable/3258 microUSB panel mount.step" \
    "6d6734c61bc3b40947d544c2f61da4b6de5a8ca6"
fetch_one "adafruit-3258-microusb-panel-mount-stl.stl" \
    "3258 USB Panel Mout Cable/3258 microUSB panel mount.stl" \
    "ddd6f0d54eba1a15235f05c013eb57d6abd4f949"
fetch_one "adafruit-3898-400mah-battery-step.step" \
    "3898 400mah Battery/3898 400mAh Battery.step" \
    "e5e42a71381bbccf5880eb55df94358c9c91ac14"
fetch_one "adafruit-3898-400mah-battery-stl.stl" \
    "3898 400mah Battery/3898 400mAh Battery.stl" \
    "0a389e2402e608840e09c29020c69b5eaa7395c9"
fetch_one "adafruit-3923-mini-oval-speaker-step.step" \
    "3923 Mini Oval Speaker/3923 Mini Oval Speaker.step" \
    "c5a4c79a35a363c05201fd2f66117eb8172716fe"
fetch_one "adafruit-3923-mini-oval-speaker-stl.stl" \
    "3923 Mini Oval Speaker/3923 Mini Oval Speaker.stl" \
    "c7417b5b2c48d97bbf8284aaac54c20139ecb363"
fetch_one "adafruit-4056-usbc-microb-cable-step.step" \
    "4056 USBC-to-microB-cable/4056-USB-C-to-micro-B.step" \
    "e43f9405f562837dc84665db859bd411ae470705"
fetch_one "adafruit-4056-usbc-microb-cable-stl.stl" \
    "4056 USBC-to-microB-cable/4056-USB-C-to-micro-B.stl" \
    "04365fc5e805c6a643e22cb4e97a0afffc697821"
fetch_one "adafruit-413-large-solenoid-step.step" \
    "413 Large Solenoid/413 Large Solenoid.step" \
    "57ce955ca0656982874a751c995dfc570b37bfd8"
fetch_one "adafruit-413-large-solenoid-stl.stl" \
    "413 Large Solenoid/413 Large Solenoid.stl" \
    "18d15347f10ef2f912c14408efe775cb380d8204"
fetch_one "adafruit-5128-macropad-bottom-step.step" \
    "5128 MacroPad RP2040 Kit/5128 MacroPad RP2040-bottom-plate.step" \
    "e235286429efca235ea72dc40553d19b764754ac"
fetch_one "adafruit-5128-macropad-bottom-stl.stl" \
    "5128 MacroPad RP2040 Kit/5128 MacroPad RP2040-bottom-plate.stl" \
    "22a54360ca8424792e6a7249c7b6c5aad8709fcd"
fetch_one "adafruit-5128-macropad-keyplate-step.step" \
    "5128 MacroPad RP2040 Kit/5128 MacroPad RP2040-Keyplate.step" \
    "b396e9e6f47c68bdbde95a0ea78c08988c10ad6b"
fetch_one "adafruit-5128-macropad-keyplate-stl.stl" \
    "5128 MacroPad RP2040 Kit/5128 MacroPad RP2040-Keyplate.stl" \
    "52a4bb16940d1194e8d0fdeaf71607f08af9be37"

fetch_http_snapshot "smithsonian-usnm174698-cuneiform3-right.ply" \
    "https://3d-api.si.edu/content/document/3d_package:c9105426-6818-4c25-b04c-135e79203b20/resources/USNM174698_cuneiform3_right_-300.ply" \
    "30ecab43080af0900a339faa7ae0196e6132966b515334e37f121df39b544286"
