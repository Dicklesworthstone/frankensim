# fs-img — CONTRACT

In-house image plumbing (plan §10.5): PNG and OpenEXR writers/readers, an
à-trous denoiser whose outputs are permanently labeled biased, and
deterministic film/display transforms. Everything is pure Rust from first
principles — no image, compression, or color-management crates (P1).

Ambition tags: PNG/EXR subset writers [S]; denoiser (SVGF-lineage single
frame) [S]; film transforms [S].

## Purpose and layer

Layer **L5** (LUMEN support). Runtime deps: `std`, `fs-math` (deterministic
`pow` for sRGB encoding). Renders ship in EXR (lossless f32/f16 AOVs); PNG
is the preview/report format. The Ledger stores both as artifacts, so the
readers exist to round-trip **our own** outputs, not the world's files.

## Public types and semantics

- `PngColor` (Gray/Rgb/Rgba), `write_png8`, `write_png16` — 8/16-bit PNG
  with sRGB chunk, filter type None on every row, zlib streams built from
  STORED deflate blocks. `read_png` → `DecodedPng` (`bytes`, `samples16()`).
- `Channel { name, ty, data: Vec<f32> }`, `PixelType` (Half/Float),
  `write_exr` — single-part scanline EXR, version 2, NONE compression,
  channels stored in the spec's alphabetical order regardless of argument
  order. `ExrAttribute { name, ty, value }` and
  `write_exr_with_attributes` add alphabetically ordered opaque custom header
  attributes; `DecodedExr.attributes` preserves their names, type names, and
  payload bytes exactly. `SOURCE_ARTIFACT_HASH_ATTRIBUTE` standardizes the
  L5/L6 composition key without teaching this crate ledger semantics. An empty
  attribute slice is byte-identical to `write_exr`. `read_exr` → `DecodedExr`.
  `f32_to_f16_bits` / `f16_bits_to_f32` — IEEE 754 half conversion with
  round-to-nearest-even, including subnormals, ±inf, and NaN (payload
  preserved as a quiet bit).
- `LabeledPlane { width, height, data, provenance }` with mandatory
  `PixelProvenance` tag: `RawEstimate` or `BiasedDenoised { iterations }`.
  `atrous_denoise(noisy, albedo?, params)` — iterated 5×5 B3-spline à-trous
  convolution with edge-stopping weights; the result is PERMANENTLY tagged
  `BiasedDenoised`. `mse` is the improvement metric.
- `film`: `exposure`, `white_balance`, `hable_filmic` (Hable/Uncharted 2
  operator, W = 11.2), `srgb_encode` (via `fs_math::det::pow`), `quantize8`,
  `display_transform` (the full chain, HDR f32 → display u8).

## Invariants

1. **Byte-exact deterministic encodes (P2)**: same pixels → same bytes,
   every run, every ISA. Writers are pure integer/bit code; the only float
   math is f32→f16 conversion, which is exact bit manipulation.
2. **Lossless AOV round-trip**: `read_exr(write_exr(x))` returns exactly
   the input samples for FLOAT channels; HALF channels return exactly the
   RNE-converted value (and exactly the input when it is representable).
   Custom EXR attribute payloads, including NUL and non-UTF-8 bytes, round-trip
   exactly; built-in names cannot be shadowed.
3. **The bias label cannot be dropped**: `atrous_denoise` output is always
   `BiasedDenoised`; there is no API to relabel a plane `RawEstimate`.
4. **Structured rejection**: readers never decode garbage silently — every
   checksum (CRC-32, Adler-32) is verified, every length is bounds-checked,
   truncation at any byte fails.
5. Half round-trip: `f32_to_f16_bits(f16_bits_to_f32(h)) == h` for every
   finite half (tested exhaustively).

## Error model

`ImgError`: `Shape { expected, got, context }` (buffer/shape disagreement),
`Malformed { what }` (structurally invalid bytes — corruption), and
`Unsupported { what }` (valid-looking bytes outside our subset). No panics
on any byte input to the readers (fuzzed); writers panic never — shape
defects return `Err`.

## Determinism class

**D0 (bit-exact)** for both writers and all film transforms (`srgb_encode`
uses `fs_math::det::pow`, not libm). The denoiser accumulates in f64 with a
fixed traversal order and uses `f64::exp`; it is run-to-run deterministic on
a given target and documented as cross-ISA reproducible only to the extent
`f64::exp` is (edge-stopping weights; the *tagged bias* is the honest
qualifier, not the last ulp).

## Cancellation behavior

All entry points are bounded, allocation-up-front, single-pass functions —
no long-running loops that need cancellation tokens (P7 satisfied by
boundedness). Callers cancel between frames.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` emits canonical `fs_obs::EventKind::ConformanceCase`
aggregate verdicts under suite `fs-img/conformance`. Passing cases use `Info`,
failures use `Error`; every reached verdict passes the failure-record lint,
serializes through `to_jsonl`, validates against the fs-obs wire schema, and
prints before its final assertion. Fixed-input cases im-001/002 record seed
zero. The randomized fixtures record their literal root input seeds:
`0x5EED_D401_5E00_0003` for im-003 and `0x5EED_F077_0000_0004` for im-004.
There is no execution seed in this suite.

When `sips` is unavailable, im-002 emits a validated `Warn` `Custom` capability
row under the same suite/case identity and returns without fabricating an
aggregate verdict. Im-003 emits a validated `Info` `Custom` MSE companion under
the distinct `im-003/measurement` scope; finite measurements remain JSON
numbers, non-finite measurements are represented as `null`, and the companion
carries the same root input seed. The scope and Custom name distinguish
supplemental evidence from the aggregate decision without reusing its sequence
identity.

Fixture construction and intermediate `unwrap`/`expect`/assertion operations
remain outside the aggregate boundary. If one aborts before `verdict`, no
aggregate event is fabricated; absence of a verdict means the case did not
complete, never that it passed.

- **im-001** — PNG8/PNG16/EXR encodes are byte-identical across repeated
  calls; PNG round-trips samples exactly; EXR AOV set (FLOAT + on-grid
  HALF) and source-artifact-hash metadata round-trip losslessly; empty
  metadata preserves the legacy EXR bytes exactly.
- **im-002** — external oracle: macOS `sips` (CoreImage) parses our PNG and
  EXR and reports the correct dimensions. Dev-only; **skips with an explicit
  JSON note** when `sips` is absent (Linux CI).
- **im-003** — the denoiser reduces MSE by >2× on a seeded noisy-gradient
  fixture and the output carries `BiasedDenoised { iterations: 3 }`.
- **im-004** — 2000 seeded junk buffers produce 4000 PNG/EXR reader attempts,
  all of which are rejected; a valid PNG truncated at **every** prefix length
  is rejected.

Unit tests additionally pin CRC-32/Adler-32 known-answer vectors, PNG
signature/chunk structure, the exhaustive f16 round-trip, film-transform
known answers, and denoiser partition-of-unity on constant images.

## No-claim boundaries

- **Not general-purpose decoders.** `read_png`/`read_exr` cover exactly the
  subset our writers emit (None-filtered stored-block PNG; single-part
  scanline NONE-compression v2 EXR) and return structured `Unsupported`
  errors beyond it. They are for round-trips and Ledger artifacts.
- **Metadata is opaque at L5.** `fs-img` validates EXR header syntax and
  preserves custom attribute bytes; it does not validate hash algorithms,
  artifact existence, lineage, or whether a claimed source hash matches the
  rendered field. Those checks belong to the L6 composition layer.
- **No compression-ratio claim.** PNG zlib streams use STORED deflate
  blocks: universally decodable, ~0% compression. EXR is NONE compression.
  Compact storage is out of scope for this bead.
- **No color management beyond sRGB.** One transfer function (IEC sRGB via
  deterministic `pow`), one tone map (Hable). No ICC profiles, no wide
  gamuts.
- **The denoiser is biased, and says so in the type system.** Its output
  must never be used as ground truth in a comparison; the Gauntlet compares
  raw estimates.
- **`sips` oracle is dev-only.** External validation runs where macOS is
  available; CI relies on the structural + round-trip suites.
- No SIMD, no threading — planes at preview sizes; performance is not a
  claim here.
