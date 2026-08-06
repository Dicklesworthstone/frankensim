# fs-img — CONTRACT

In-house image plumbing (plan §10.5): PNG and OpenEXR writers/readers, an
à-trous denoiser whose outputs are permanently labeled biased, and
deterministic film/display transforms. Everything is pure Rust from first
principles — no image, compression, or color-management crates (P1).

Ambition tags: PNG/EXR subset writers [S]; denoiser (SVGF-lineage single
frame) [S]; film transforms [S].

## Purpose and layer

Layer **L5** (LUMEN support). Runtime deps: `std`, `fs-blake3` (typed content
identities), and `fs-math` (deterministic `pow` for sRGB encoding). Renders
ship in EXR (lossless f32/f16 AOVs); PNG is the preview/report format. The
Ledger stores both as artifacts, so the readers exist to round-trip **our
own** outputs, not the world's files.

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
  attribute slice is byte-identical to `write_exr`.
  `exr_write_requirements_for_layout` admits dimensions/channel metadata
  without image planes; `exr_write_requirements` additionally validates owned
  plane shapes. Both return the exact NONE-compressed encoded length and the
  logical reference-vector scratch. `write_exr_with_attributes_budgeted`
  refuses either caller ceiling before allocation, uses fallible reservations,
  and preserves the compatibility writer's bytes exactly. `read_exr` →
  `DecodedExr`.
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
  `display_transform` (the legacy full chain, HDR f32 → display u8).
- `CinematicColorConfig` freezes the version-one scene-linear-to-display
  transform: linear-sRGB/D65 input, sRGB/D65 output, Hable-v1 or the explicitly
  named Narkowicz five-coefficient ACES-filmic *fit*, clip or RGB-ratio gamut
  handling,
  counted clamp-to-zero negative policy, exact-power-of-two exposure, positive
  RGB white-balance gains, 8/16-bit output, keyed half-LSB dither, and optional
  `BoxBloomV1`. Its exact 64-byte canonical codec is closed, rejects
  trailing/non-canonical bytes and unknown
  tags, normalizes floating-point zero, and records every semantic parameter.
  `transform_cinematic_preview` consumes immutable planar f32 linear-RGB,
  checks dimensions, plane shapes, numeric values, a hard 8K-UHD pixel ceiling,
  and caller-supplied pixel/working-byte limits before allocation. It returns
  interleaved `CinematicPreviewSamples` compatible with `write_png8` or
  `write_png16`, plus `CinematicColorMetadata` containing canonical replay
  parameters, negative/over-range handling counts, gamut-operation counts,
  linear bright-pass/addition sums, and admitted output-plus-scratch bytes.
  A zero-strength bloom is non-canonical and must be expressed as `Disabled`.
- `sequence`: `FrameSequenceContext` binds the shot, trajectory, render-config,
  scene, producer-build, and image-profile `ContentHash` assertions;
  `ExpectedFrameArtifact` declares each keyed output and its exact byte
  reservation; and `FrameSequenceManifest` registers observations, emits or
  resumes strict canonical snapshots, audits completed byte identities, and
  produces a `FrameSequenceSeal` only after complete finalization. Artifact
  descriptors carry exact time bits, format, dimensions, sorted named channel
  types, sampling statistics, and raw/derived source linkage. The public
  lifecycle is `Incomplete` or immutable `Finalized`; exact registration
  retries return `AlreadyRecorded`, while conflicting retries refuse.

### Frame-sequence artifact contract (h7xu5.6.4)

- The complete expected inventory is admitted before rendering. Its stable key
  is `(frame_index, segment_index, role)`: segments distinguish multiple
  event-delimited renders of one presentation frame, and roles keep raw EXR,
  biased denoised EXR, display preview, and scientific overlay artifacts
  separate. Entries and canonical snapshots are sorted by this key.
- Paths are generated from typed context and descriptors, are relative to an
  unspecified artifact root, and contain no caller-selected path components.
  Moving the root cannot change any path, snapshot byte, or snapshot identity.
  Absolute paths, parent/current-directory components, platform separators,
  and alternate spellings are outside the grammar.
- Frame time is one finite binary64 value stored and compared by its exact
  `to_bits()` representation. Both signed zeros canonicalize to positive zero;
  NaN and infinity refuse. Channel descriptors are nonempty, name-unique, and
  sorted by bytewise channel name with their scalar type retained exactly, so
  caller insertion order cannot affect identity.
- A raw master has no source. Every derived row names an earlier role for the
  same frame and segment and, when complete, repeats that source artifact's
  exact byte hash. This prevents a derived output from being silently attached
  to a different registered master; it is a structural lineage assertion, not
  producer authentication or semantic reconstruction.
- Every expected file reserves its own nonzero maximum byte count before work
  begins. Admission uses the checked sum of those reservations; registration
  checks the file against its own reservation and increments a checked total
  of actual completed bytes. Pending reservations and completed actual bytes
  are distinct quantities. Equality with a limit is admitted; limit plus one
  refuses before state mutation.
- An incomplete canonical snapshot is resumable. Decode revalidates its closed
  grammar, canonical order, limits, paths, source graph, completion totals, and
  the currently available capacity for pending reservations. A finalized
  snapshot contains every expected completion and is immutable. Its
  domain-separated identity is computed from the exact canonical bytes and
  must be pinned and compared by an authority outside the snapshot; a snapshot
  cannot authenticate or bless its own identity.
- Snapshot encoding, audit, and finalization poll at artifact boundaries, with
  identity hashing additionally polling at bounded byte chunks. Cancellation,
  missing/stale observations, resource refusal, or encoding failure publishes
  no seal and leaves the mutable manifest unchanged and resumable. The state
  transition to `Finalized` occurs only after complete re-observation,
  canonical encoding, and identity calculation succeed.

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
6. **Raw-master isolation**: cinematic color accepts shared slices and has no
   mutation surface. Its result is permanently
   `DisplayReferredDerivativeV1`; it cannot be relabeled as a raw estimate.
7. **Visible highlight handling**: NaN, infinity, adjusted overflow, and
   magnitudes above `1e12` refuse in row-major/channel order. Finite negative
   and above-one channels are counted before explicit handling; no such value
   disappears silently.
8. **Budget-before-allocation**: output and optional bloom-scratch sizes use
   checked arithmetic and are compared with the caller's admitted envelope
   before `try_reserve_exact`. Allocation refusal is structured and no partial
   preview is returned.
9. **EXR exact-size admission**: the budgeted EXR writer validates shapes,
   version-2 short names, duplicate/reserved attributes, signed-i32
   scanline/header lengths, logical ordering scratch, and exact encoded output
   bytes before allocating either reference vector or output. It writes the
   channel list directly into the once-reserved output and uses stack storage
   for the fixed data window.
   Scratch accounting is logical reference payload, not allocator bookkeeping;
   caller-owned channel planes and attribute payloads remain caller-budgeted.

## Error model

`ImgError`: `Shape { expected, got, context }` (buffer/shape disagreement),
`Malformed { what }` (structurally invalid bytes — corruption),
`Unsupported { what }` (valid-looking bytes outside our subset),
`ResourceLimit` (exact logical requirement exceeds a caller ceiling),
`AllocationRefused` (fallible reservation failed after admission), and
`SizeOverflow` (checked size arithmetic cannot represent the artifact).
`CinematicColorError` separately reports stable config/canonical field paths,
shape and resource requirements, the first invalid pixel/channel/stage,
checked-arithmetic overflow, and allocation refusal. No panics occur on byte
input to the readers (fuzzed); writers and the cinematic transform return
structured errors for admitted defects.

## Determinism class

**D0 (bit-exact)** for both writers and all film transforms (`srgb_encode`
uses `fs_math::det::pow`, not libm). The denoiser accumulates in f64 with a
fixed traversal order and uses `f64::exp`; it is run-to-run deterministic on
a given target and documented as cross-ISA reproducible only to the extent
`f64::exp` is (edge-stopping weights; the *tagged bias* is the honest
qualifier, not the last ulp).

The cinematic pipeline is also D0: its tone curves use frozen arithmetic,
sRGB uses deterministic `fs-math`, dither is a specified SplitMix64-derived
keyed variate, and bloom traverses fixed row/column orders. A changed seed or
any changed semantic parameter changes the canonical configuration bytes.

## Cancellation behavior

All entry points are bounded, allocation-up-front functions. The cinematic
path is hard-capped at 8K UHD and uses linear-time sliding-window bloom rather
than a radius-squared convolution. Callers cancel between frames; this module
does not claim intra-frame cancellation latency.

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
known answers, and denoiser partition-of-unity on constant images. Cinematic
G0/G3/G5 cases cover the closed canonical codec and every semantic field,
known curve anchors/monotonicity, gray neutrality and exact exposure stops,
non-finite/shape/pixel/memory refusal, visible negative/over-range counts,
raw-input immutability, deterministic seed-sensitive dither with exact
endpoints, local/non-wrapping bloom and constant-field interiors, exact
working-byte accounting, and direct 8/16-bit PNG round trips.

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
- **One explicit display target, not general color management.** Version one
  accepts linear sRGB/D65 and emits sRGB/D65 only. It has no ICC profiles,
  chromatic-adaptation engine, wide-gamut target, HDR transfer function, or
  monitor calibration. `AcesFittedNarkowiczV1` is the published compact fit,
  not an ACES reference transform or an OCIO compatibility claim.
- **Bloom is a labeled display effect.** `BoxBloomV1` thresholds exposed,
  white-balanced scene-linear RGB and applies two normalized zero-boundary box
  passes. The recorded RGB sums make edge loss and added signal visible, but
  they are not radiometric energy, lens-scattering calibration, diffraction,
  flare, or a perceptual-quality certificate.
- **The denoiser is biased, and says so in the type system.** Its output
  must never be used as ground truth in a comparison; the Gauntlet compares
  raw estimates.
- **`sips` oracle is dev-only.** External validation runs where macOS is
  available; CI relies on the structural + round-trip suites.
- No SIMD, no threading — planes at preview sizes; performance is not a
  claim here.
