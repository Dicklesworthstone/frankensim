#!/usr/bin/env python3
"""Build the actual weak-constraint numerical source closure, without bootstrap.

This is deliberately NOT a full-workspace or full fs-ascent package check.
Temporary Cargo manifests select complete production modules; no numerical
functions or tests are rewritten or replaced. All dependency sources are from
this checkout. The normal workspace remains authoritative for integration.
"""
from __future__ import annotations

import json
from pathlib import Path
import re
import subprocess
import tempfile
import tomllib


def prepare(root: Path, out: Path) -> Path:
    package = tomllib.loads((root / "Cargo.toml").read_text())["workspace"]["package"]
    version = package["version"]

    def module(name: str, path: str) -> str:
        source = root / path
        if not source.is_file():
            raise FileNotFoundError(source)
        return f"#[path = {json.dumps(str(source))}] pub mod {name};\n"

    def crate(name: str, source: str | None, deps: tuple[str, ...] = (), *, lib: str | None = None) -> Path:
        folder = out / name
        folder.mkdir(parents=True)
        lib_path = root / lib if lib else folder / "lib.rs"
        if source is not None:
            lib_path.write_text(source)
        elif not lib_path.is_file():
            raise FileNotFoundError(lib_path)
        text = f'[package]\nname = "{name}"\nversion = {json.dumps(version)}\nedition = "2024"\n'
        text += f'[lib]\npath = {json.dumps(str(lib_path))}\n'
        text += '[dependencies]\n' + ''.join(f'{dep} = {{ path = "../{dep}" }}\n' for dep in deps)
        (folder / "Cargo.toml").write_text(text)
        return folder

    crate("fs-math", None, lib="crates/fs-math/src/lib.rs")
    crate("fs-blake3", None, lib="crates/fs-blake3/src/lib.rs")
    crate("fs-ad", module("revolve", "crates/fs-ad/src/revolve.rs"))
    time_source = module("adaptive", "crates/fs-time/src/adaptive.rs") + "pub use adaptive::*;\n"
    crate("fs-time", time_source, ("fs-math", "fs-blake3", "fs-ad"))

    # Import the exact owned callback alias, not a test-only substitute.
    owner = (root / "crates/fs-ascent/src/lib.rs").read_text()
    alias = re.search(r"pub type FnGrad<'a>\s*=.*?;", owner, re.S)
    if alias is None:
        raise RuntimeError("fs-ascent FnGrad alias not found")
    source = ''.join(module(name, f"crates/fs-ascent/src/{name}.rs") for name in ("stop", "wolfe", "lbfgs"))
    source += "pub use stop::{StopReason, StopRule};\npub use lbfgs::{LbfgsError, LbfgsReport, LbfgsState};\n"
    source += alias.group() + "\npub mod transient {\n"
    source += module("variational", "crates/fs-ascent/src/transient/variational.rs") + "}\n"
    ascent = crate("fs-ascent", source, ("fs-time", "fs-math"))
    example = root / "crates/fs-ascent/examples/weak_constraint_heat.rs"
    if not example.is_file():
        raise FileNotFoundError(example)
    with (ascent / "Cargo.toml").open("a") as manifest:
        manifest.write(f'[[example]]\nname = "weak_constraint_heat"\npath = {json.dumps(str(example))}\n')
    (out / "Cargo.toml").write_text('[workspace]\nmembers = ["fs-math", "fs-blake3", "fs-ad", "fs-time", "fs-ascent"]\nresolver = "3"\n')
    return out / "Cargo.toml"


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="frankensim-variational-") as folder:
        manifest = prepare(root, Path(folder))
        for profile in ([], ["--release"]):
            base = ["cargo", "+stable", "test", "--offline", "--manifest-path", str(manifest), "-p", "fs-ascent", *profile]
            subprocess.run([*base, "--lib", "transient::variational::"], check=True, cwd=manifest.parent)
            subprocess.run([*base, "--example", "weak_constraint_heat"], check=True, cwd=manifest.parent)
        subprocess.run(["cargo", "+stable", "run", "--offline", "--release", "--manifest-path", str(manifest), "-p", "fs-ascent", "--example", "weak_constraint_heat"], check=True, cwd=manifest.parent)


if __name__ == "__main__":
    main()
