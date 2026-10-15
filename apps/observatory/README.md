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

## The workspace

The left inspector edits the **next** native experiment; changing a slider never
alters a previously accepted result. Fixed ends expose engineering strain and its
limit. Solver refusals, cancellations and bad imports leave the accepted run intact.
Only completed native results enter the run library. The native example currently
returns 50 ms of acoustic data: playback is intentionally short, not looped into a
fabricated performance.

The central cutaway supports drag-to-orbit, wheel/button zoom, keyboard rotation,
layer picking, and temperature or signed-stress coloring. Cold and heated are
separate physical states, not interpolation through an invented heating history.
Cross-sectional size is exaggerated and explicitly labeled. Each layer has a
uniform field value; lighting and wire lines are visual aids, not additional
simulation samples. Missing geometry remains a labeled reference diagram.

Traces retain explicit units, separate incompatible-unit channels, and expose
nearest-sample inspection and window zoom. Pixel-bin min/max envelopes keep narrow
peaks visible without changing exported data. Pin a run to overlay its last
matching-unit channel as a dashed baseline. There is no automatic time alignment,
resampling, peak normalization or out-of-range baseline extrapolation.

Uniform audio-rate pressure data also receives a full-record, single-sided Hann
amplitude spectrum. This is signal analysis, not a second physical model. The two
listening buttons share one explicit full-scale pressure in Pa; clipping is counted
and reported, never hidden by per-channel normalization. Playback requires a click.

The run library keeps up to eight runs in memory and as many as fit its 3 MiB browser
storage budget. Storage denial or quota exhaustion does not lose the open result;
export JSON for durable retention. Restored and imported results are always labeled
as such, even if their metadata names a native solver. JSON and numeric CSV export
retain all samples; PNG export captures the current camera and displayed source tag.

Dark/light themes, responsive layouts, reduced-motion support, visible keyboard
focus, accessible native controls and a command palette are included. Use **R** to
run, **I** to import, **Ctrl/Command K** for actions and **?** for the guide. Focus
the model for arrow-key orbit, **+/-** zoom and **0** reset; focus the trace for
adjacent-sample inspection. Escape cancels an active experiment outside dialogs.

## Import format

Numeric CSV uses `time_s` followed by up to twelve uniquely named channels:

```csv
time_s,receiver_pressure_pa,core_temperature_k
0,0.00001,300
0.001,-0.00002,300.1
```

These numbers illustrate the format only. CSV has no geometry or authenticated
provenance. The viewer infers units only from the declared suffixes `_pa`, `_k`,
`_n`, `_j`, `_m`, `_w`, `_hz`; other columns remain separately plotted with
unspecified units. Times must be strictly increasing and values finite. Maximum:
32,768 samples/channel and 8 MiB/file. Quoted headers, BOM and CRLF are accepted.

Snapshot JSON uses `schema: "frankensim.observatory/v1"`, a `title`,
`source: {model, material, notes}`, and `traces: {time_s, series}`. Each series has
`label`, `unit`, and `values`. Optional `inputs`, `section`, and `metrics` carry the
heated-section experiment's reported state. `protocol.mjs` is the shared normative
validator; export a run to inspect the complete shape. Unknown JSON fields are not
preserved. A snapshot is never an authenticated certificate.

## Focused verification

The dependency-free Node suite covers numerical import boundaries, native-output
parsing, HTTP admission, actual subprocess lifecycle using a test double, FFT
amplitude/frequency, extrema-preserving decimation, and CSV round trips:

```sh
node --test apps/observatory/tests/*.test.mjs
```

An optional Python Playwright check exercises the actual browser UI. With the local
server already running and Playwright/Chromium installed:

```sh
python apps/observatory/tests/browser_smoke.py
```

`--browser /path/to/chromium` selects an existing browser; `--screenshots DIR`
retains a screenshot whose input is explicitly labeled as a UI test fixture.
API responses in this smoke check are test doubles, never Rust physics evidence.

In hosts that prohibit browser navigation, `--in-memory` loads the same HTML/CSS/JS
through Playwright, changing only module URLs. This verifies browser rendering and
interactions but **not** served-asset/CSP integration or durable browser storage.
No browser security policy is disabled. The Node tests independently cover actual
HTTP. Native Rust execution, platform listening fidelity and Safari/Firefox remain
separate checks; neither fixture suite establishes them.
