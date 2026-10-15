# FrankenSim Observatory

A local browser workspace for physical results. The native experiment adapter runs
`fs-couple`'s existing `heated_section_sound` executable; there is no JavaScript
physics surrogate. JSON/CSV visualization also works without a Rust toolchain.

## Start locally

Requires Node.js 22 or newer. No npm install, bundler, CDN, account or cloud service.
From the repository root:

```sh
node apps/observatory/server.mjs --open
```

The server binds **127.0.0.1:7070 only**. `--port 7171` changes the port.
Build the optional native engine using the repository's normal pinned-sibling setup:

```sh
cargo build --release -p fs-couple --example heated_section_sound
```

Click the engine indicator to refresh after building. An alternative binary can be
selected with `--engine /absolute/path/to/heated_section_sound`. The browser cannot
select executables, issue shell commands, or browse the filesystem.

For your own material data, start with both `--core-pack PATH --sleeve-pack PATH`.
These are canonical packs already accepted by the example, not invented presets.
Both need expansion data for fixed-end mode. Resistivity scaling is disabled for
packs; material substitution stays data-only.

## Boundaries

The adapter consumes the example's existing CSV and structured diagnostic fields.
It fails closed when required fields or the 800-sample/16 kHz layout change. Geometry
is the example's fixed 0.65 m length and 0.25/0.5 mm layer radii. A future geometry
change must update this adapter. It never fabricates intermediate thermal frames.

Snapshots are visualization data, **not authenticated Design Ledger certificates**.
Synthetic material cards, supplied packs, imported data and locally executed runs
remain distinguishable. Thermal residual is accounting, not an accuracy certificate.

The HTTP companion admits one bounded native child at a time, kills it on cancel,
disconnect or timeout, and never publishes partial results. It restricts Host and
Origin, requires a per-session token for execution, and serves only an asset allowlist.
This is a single-user loopback tool, not an Internet-facing service.

```sh
node --test apps/observatory/tests/*.test.mjs
```

Process-adapter tests use an explicitly synthetic executable fixture; they are not
native Rust integration or physical validation. Frontend usage and browser checks
are documented alongside the visual workspace implementation.
