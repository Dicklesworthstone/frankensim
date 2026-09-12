# Dependency update and release preparation

Date: 2026-09-12 UTC. Release status: **blocked; nothing published**.

## Completed dependency updates

`apps/wright-flyer` was updated one dependency at a time. Each step passed
the existing 21 facade/harness tests; TypeScript, Vite, Three and its type
definitions also passed `tsc --noEmit`. The final full Node test run passed
295 tests, and `npm audit --json` reported zero known vulnerabilities.
The facade tests explicitly did not exercise a live native engine (`WF_PKG`
was unset); these results are not full browser/native simulation validation.

| Dependency | Previous declaration | New declaration |
| --- | --- | --- |
| puppeteer-core | ^25.8.0 | ^25.10.0 |
| typescript | ~5.7.0 | ~7.0.2 |
| vite | ^6.0.0 | ^8.3.0 |
| three | ^0.172.0 | ^0.186.0 |
| @types/three | ^0.185.4 | ^0.186.0 |

Research used the upstream [Puppeteer release notes](https://github.com/puppeteer/puppeteer/releases/tag/puppeteer-core-v25.10.0),
[TypeScript 7 announcement](https://devblogs.microsoft.com/typescript/announcing-typescript-7-0/),
[Vite migration guide](https://vite.dev/guide/migration), and
[Three migration guide](https://github.com/mrdoob/three.js/wiki/Migration-Guide).
Vite resolved the existing ES2022/module-worker configuration successfully.
Three's deprecated `PCFSoftShadowMap` use was replaced with `PCFShadowMap`;
the new implementation also provides soft shadows. Geometry and the GLTF/SMAA
module imports passed a smoke check. Puppeteer connected to Chrome
153.0.8010.36, but browser teardown was slow; this was not a full browser test.

Final local logs:

- `/tmp/frankensim-flyer-final-types-20260912.log`
- `/tmp/frankensim-flyer-all-tests-20260912.log`
- `/tmp/frankensim-flyer-final-audit-20260912.json`

The production build passed through DSR's explicit browser lane:
`DSR_REPOS_FILE=/tmp/frankensim-browser-dsr-20260912.yaml dsr quality --tool frankensim-browser --work-dir /Users/jemanuel/projects/frankensim/apps/wright-flyer`.
Its Vite build passed in 57.776 seconds and produced
`/tmp/frankensim-flyer-dist-20260912`. The first aggregate did not pass:
three parallel timing tests failed under host load, and shared source state
moved during the run. The serial retry passed all 295 tests, the typecheck,
and the production build: DSR 3/3 with a stable source-state receipt at
`/Users/jemanuel/.local/state/dsr/quality-logs/frankensim-browser/20260911T221237-38533/receipt.json`.
No test assertions were weakened. The narrow browser lane does not replace the
ten configured repository checks. Log: `/tmp/frankensim-browser-dsr-20260912.log`.
The earlier required-remote `rch exec -- npm run build` refused with
RCH-E301 (unsupported command classification).

UBS scanned the changed source file: zero critical findings, 25 heuristic
warnings in pre-existing code (typed Three property access, already guarded
button access, and module-lifetime listeners). No warning implicated the
one-line shadow-map migration. UBS reported its ast-grep rule-pack mode was
unavailable; this is not full AST scanner coverage.

## Rust dependencies

The ordinary top-level registry dependencies are current according to the
live crates.io API: serde 1.0.229, serde_json 1.0.151, libm 0.2.16, and
the isolated high-precision oracle's rug 1.30.0. Target-specific scanning also
found wasm-bindgen: root/fs-wasm locks use 0.2.126 and fs-flyer-wasm uses
0.2.127, versus latest 0.2.128. That update remains pending matching CLI/glue
regeneration and WASM qualification; generated artifacts were not relabeled.
The direct getrandom 0.4 dependency resolves to current 0.4.3; the separate
transitive 0.2.17 line was not forced across a breaking version boundary.
No Rust manifest or lockfile changes were made. Sibling path dependencies,
nightly toolchains, and generated WASM package metadata were preserved.
This is not an audit of every transitive dependency owned by sibling projects.

## Reference environment updates

Updated `tools/vvref` one dependency at a time under CPython 3.12.12, with
the environment and uv cache outside the repository. Each upgrade passed
both known-answer checks and solved all four Level-B thermal case decks.
The 126 probe/specification rows remained identical after scikit-fem and
NumPy updates. SciPy changed probe temperatures by at most
8.185452315956354e-12 K; all agreed with the old environment within 1e-9 K.
The final blake3 update preserved those results and every deck/mesh hash.
The final run also treated Python deprecation warnings as errors.
Existing frozen references were not regenerated.
`pip-audit` checked all four packages from the frozen lock export and reported
no known vulnerabilities (`/tmp/frankensim-vvref-audit-20260912.log`).

| Dependency | Previous pin | Updated pin |
| --- | --- | --- |
| scikit-fem | 11.0.0 | 12.0.2 |
| numpy | 2.3.1 | 2.5.3 |
| scipy | 1.16.0 | 1.18.1 |
| blake3 | 1.0.5 | 1.0.9 |

Reviewed upstream [scikit-fem changes](https://github.com/kinnala/scikit-fem/releases),
[NumPy 2.5 migration notes](https://numpy.org/devdocs/release/2.5.0-notes.html),
[SciPy 1.17](https://docs.scipy.org/doc/scipy/release/1.17.0-notes.html)
and [1.18 changes](https://docs.scipy.org/doc/scipy/release/1.18.0-notes.html),
and [blake3 releases](https://github.com/oconnor663/blake3-py/releases).
Candidate outputs: `/tmp/frankensim-vvref-{before,skfem12,numpy25,scipy118,final}-20260912.tsv`.

## Release blockers and remaining work

- Live origin has no tags and the GitHub releases API returns zero releases.
  CASS searches found DSR build references but no prior `dsr release` execution.
  The historical index was stale, so this is not an exhaustive session-history claim.
- Canonical `scripts/ci/checkout_constellation.sh --verify-only` refused
  asupersync: expected `03a0a298d07f565c56b3d80ad85cf2bbc69ad3c2`, observed
  `d47a54b5063173dbaaee3e83e706702355233870`, with 336 dirty paths observed.
  The expected commit is an ancestor of the observed commit. The lock is stale,
  and the live sibling is also dirty; no shared checkout was reset or repinned.
  The same preflight subsequently refused `franken_networkx` (expected
  `ab73033517f47250ac9477148fc055f148e838b0`, actual
  `ccb33b5be7950bb759db3a36700d3a21a452e645`) and `franken_numpy` (expected
  `9b6b5828317dbeda36a6a9e53c8fa754527f0d0c`, actual
  `35e21651211a49162e32d9fdf2d28c60d13b7d3d`). The remaining expensive scan
  was stopped after these decisive refusals (wrapper exit 143); this was not
  a complete seven-sibling preflight. Only asupersync ancestry was classified.
- DSR dry run planned ten checks and returned exit 2 (not a passing gate).
  `/tmp/frankensim-release-quality-dryrun-20260912.log` records the inventory.
  Full DSR tests, native release builds and installer qualification remain open.
- The local DSR registration omitted the `frankensim` CLI. Added it to
  `workspace_binaries` in `~/.config/dsr/repos.d/frankensim.yaml`.
  The configured artifact target remains macOS ARM64; no artifacts were built.
- The shared tree contains extensive staged and unstaged peer work. Establish
  a coherent release commit and clean, qualified constellation before versioning.
- crates.io publication is not configured end to end: the CLI's internal path
  dependencies lack registry versions and the closure includes nonpublishable
  packages. Qualify a publishable dependency closure before attempting uploads.
- No FrankenSim Homebrew formula was found in the local `homebrew-tap`.
  A formula needs a verified published artifact and checksum first.
- Finish the reference-environment upgrades, obtain production browser build
  evidence, pass DSR gates, then commit/version/tag/build/publish and verify the
  release. No GitHub Actions were dispatched or used as verification.
