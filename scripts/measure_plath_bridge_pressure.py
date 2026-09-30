#!/usr/bin/env python3
"""Inspect paired piano bridge acceleration and ear pressure without a 4.6-GB download.

Source: Niko Plath, Zenodo record 3274772 (2019), Data.zip. The source
describes an anechoic, lid-off concert grand; its model is not identified.
This is a descriptive paired-signal check, not a bridge-mobility measurement.
Only Python's standard library and NumPy are needed.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
import io
import json
import math
import urllib.request
import zipfile

import numpy as np


URL = "https://zenodo.org/api/records/3274772/files/Data.zip/content"
ARCHIVE_BYTES = 4_611_115_433
SAMPLE_RATE = 50_000
BANDS = ((1_200, 3_000), (3_000, 8_000), (8_000, 16_000))


class RangeFile(io.RawIOBase):
    """Seekable ZIP input backed by checked, bounded HTTP byte ranges."""

    def __init__(self) -> None:
        self.pos = 0
        self.block_bytes = 1 << 20
        self.cache: OrderedDict[int, bytes] = OrderedDict()
        self.fetched_bytes = 0

    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            position = offset
        elif whence == io.SEEK_CUR:
            position = self.pos + offset
        elif whence == io.SEEK_END:
            position = ARCHIVE_BYTES + offset
        else:
            raise ValueError("invalid seek origin")
        if not 0 <= position <= ARCHIVE_BYTES:
            raise ValueError("seek outside pinned archive")
        self.pos = position
        return position

    def _block(self, number: int) -> bytes:
        if number in self.cache:
            self.cache.move_to_end(number)
            return self.cache[number]
        first = number * self.block_bytes
        last = min(first + self.block_bytes, ARCHIVE_BYTES) - 1
        request = urllib.request.Request(
            URL,
            headers={
                "Range": f"bytes={first}-{last}",
                "User-Agent": "OpenAI File Downloader, XaiImageApiFetch/1.0",
            },
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            expected = f"bytes {first}-{last}/{ARCHIVE_BYTES}"
            if response.status != 206 or response.headers.get("Content-Range") != expected:
                raise OSError("Zenodo did not return the requested pinned byte range")
            data = response.read()
        if len(data) != last - first + 1:
            raise OSError("truncated Zenodo range response")
        self.fetched_bytes += len(data)
        self.cache[number] = data
        if len(self.cache) > 32:
            self.cache.popitem(last=False)
        return data

    def read(self, count: int = -1) -> bytes:
        count = min(ARCHIVE_BYTES - self.pos, ARCHIVE_BYTES - self.pos if count < 0 else count)
        if count > 32 * self.block_bytes:
            raise OSError("refusing an unbounded read of the multi-GB source archive")
        result = bytearray()
        while count:
            block = self._block(self.pos // self.block_bytes)
            offset = self.pos % self.block_bytes
            part = block[offset : offset + count]
            result.extend(part)
            self.pos += len(part)
            count -= len(part)
        return bytes(result)


def power_spectrum(values: np.ndarray) -> np.ndarray:
    frame_size = 4096
    hop = frame_size // 2
    if len(values) < frame_size:
        raise ValueError("analysis window is too short")
    window = np.hanning(frame_size)
    powers = [
        np.abs(np.fft.rfft(values[start : start + frame_size] * window)) ** 2
        for start in range(0, len(values) - frame_size + 1, hop)
    ]
    return np.mean(powers, axis=0)


def analyze(archive: zipfile.ZipFile, key: int, state: str, take: int) -> dict:
    # The source ZIP uses two underscores in PROD07 names, one in PROD08.
    candidates = (
        f"D1_{state}_{key:02d}__0000_M{take:04d}.csv",
        f"D1_{state}_{key:02d}_0000_M{take:04d}.csv",
    )
    matches = [name for name in candidates if name in archive.NameToInfo]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one source member for {state} key {key} take {take}")
    member = matches[0]
    info = archive.getinfo(member)
    # ZipFile.read checks this member's CRC. The whole 4.6-GB MD5 is not checked.
    values = np.loadtxt(io.BytesIO(archive.read(member)), delimiter=";", skiprows=1)
    if values.ndim != 2 or values.shape[1] != 5:
        raise ValueError(f"unexpected five-channel format: {member}")
    times = values[:, 0]
    if not np.allclose(np.diff(times), 1 / SAMPLE_RATE, rtol=0, atol=1e-8):
        raise ValueError(f"unexpected timebase: {member}")
    active = (times >= 0.02) & (times < 0.28)
    noise = (times >= -0.28) & (times < -0.02)
    frequencies = np.fft.rfftfreq(4096, 1 / SAMPLE_RATE)
    channels = ("bridge_m_s2", "left_pa", "right_pa")
    active_power = [power_spectrum(values[active, col]) for col in (2, 3, 4)]
    noise_power = [power_spectrum(values[noise, col]) for col in (2, 3, 4)]
    bands = []
    for low, high in BANDS:
        indices = (frequencies >= low) & (frequencies < high)
        signal = [float(np.sum(power[indices])) for power in active_power]
        background = [float(np.sum(power[indices])) for power in noise_power]
        snr = [10 * math.log10(s / n) if s > 0 and n > 0 else None
               for s, n in zip(signal, background)]
        pressure_per_acceleration = [
            math.sqrt(signal[col] / signal[0]) if signal[0] > 0 else None
            for col in (1, 2)
        ]
        bands.append({
            "hz": [low, high],
            "active_to_pretrigger_db": dict(zip(channels, snr)),
            "pressure_per_bridge_acceleration_pa_s2_per_m": pressure_per_acceleration,
            "usable_at_10db_snr": all(level is not None and level >= 10 for level in snr),
        })
    return {
        "member": member,
        "member_crc32": f"{info.CRC:08x}",
        "sample_count": len(times),
        "bands": bands,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keys", default="40,49,64", help="piano key indices: C4=40, A4=49, C6=64")
    parser.add_argument("--state", choices=("PROD07", "PROD08"), default="PROD07")
    parser.add_argument("--takes", default="0", help="comma-separated take indices 0..4")
    args = parser.parse_args()
    try:
        keys = [int(item) for item in args.keys.split(",")]
        takes = [int(item) for item in args.takes.split(",")]
    except ValueError:
        parser.error("--keys and --takes need comma-separated integer indices")
    if not keys or any(not 1 <= key <= 88 for key in keys) or len(set(keys)) != len(keys):
        parser.error("--keys needs unique indices from 1 through 88")
    if not takes or any(not 0 <= take <= 4 for take in takes) or len(set(takes)) != len(takes):
        parser.error("--takes needs unique indices from 0 through 4")
    source = RangeFile()
    with zipfile.ZipFile(source) as archive:
        results = [analyze(archive, key, args.state, take) for take in takes for key in keys]
    print(json.dumps({
        "source": URL,
        "record_doi": "10.5281/zenodo.3274772",
        "archive_bytes": ARCHIVE_BYTES,
        "range_bytes_fetched": source.fetched_bytes,
        "window_seconds": [0.02, 0.28],
        "pretrigger_seconds": [-0.28, -0.02],
        "interpretation": "paired band-energy ratio, not a force-to-velocity mobility or a confirmed Steinway D",
        "whole_archive_md5_verified": False,
        "notes": results,
    }, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
