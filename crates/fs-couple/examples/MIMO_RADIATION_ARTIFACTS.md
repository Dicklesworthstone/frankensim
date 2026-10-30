# Computed and imported mutual radiation

The native flute and material-modal valve formats can select the same original
multiport impedance artifact. Each source goes through
`bernoulli_aperture::radiation::environment::multiport::input`, then the existing
`MimoEnvironment` and `WaveguideNetwork` terminal host. Coupled outlet pressure
returns to the live jet or moving valve. The microphone observes the same
ordered apertures in the artifact's exact exterior gas.

The `mimo` selection is an importer for externally produced coefficients. It does not calculate
a new radiation model when geometry, materials or atmosphere change, or attest
a declared offline bake. Supply a corresponding new artifact when the changed
source requires one. A mismatch refuses construction before audio is produced.
The included operators are **manufactured integration data**, with no measured
instrument or broadband radiation accuracy claim.

## Compute radiation from actual geometry and gas

The `computed-baffled` selection derives the complete pressure/volume-flow
matrix from the performance's actual circular outlet radii, ordered planar
centers and resolved exterior gas. Changing the atmosphere or aperture data
recomputes the source. It requires no radiation artifact, material-name rule,
or hand-edited impedance coefficients. The existing `fs-bem` baffled-piston
array supplies the self and mutual terms; shared `fs-vfit` identifies the
proper rational matrix, and the same `MimoEnvironment` terminal owner feeds
reaction back to the live instrument and supplies microphone flows.

For flute inputs:

```text
radiation computed-baffled 120 65 4 8 12 20 0.02 67108864 2 100 1e-8 64 1048576 257
```

[flue-computed-radiation.performance](flue-computed-radiation.performance)
is a complete 48 kHz example with two different physical openings and regional
interior gases. Render it through the existing `music_render wind` command.
This manufactured geometry is a functional example, not a measured flute.

For material-modal valve inputs, use the same fourteen fields per side:

```text
radiation computed-baffled 1
radiation_group downstream 20 65 4 8 12 20 0.01 67108864 2 100 1e-8 64 1048576 257
```

These are explicit numerical choices, in order:

| Fields | Meaning |
| --- | --- |
| `120` | Lower physical fit frequency in Hz. The upper edge is always the actual mechanical Nyquist frequency. |
| `65` | Logarithmic training frequencies; 64 interleaved frequencies are withheld, and a separate 65-point uniform grid checks source causality. |
| `4 8 12` | Three strictly increasing shared-pole orders. All three are attempted; the smallest passing fit is selected. |
| `20` | Maximum vector-fitting iterations per order. |
| `0.02` | Maximum relative matrix Frobenius error at both training and withheld points. |
| `67108864` | Separate fitting/certification scratch allowance in bytes. |
| `2 100 1e-8` | Causality subtractions, assumed high-frequency tail bound ratio and relative source-data tolerance. |
| `64 1048576 257` | Maximum real states, runtime memory bytes and interval-mean sampled-passivity grid points. |

At most 33 separated circular openings, 513 frequencies, 64 shared poles
per order and 256 real states are admitted. Requested working memory is capped
at 512 MiB, independently of the terminal's retained memory. The complete
three-order scratch requirement is admitted before the source sweep. Failed
quadrature, causality, fit accuracy or passivity returns a construction error;
the runtime does not clip residues, add damping or fall back to independent
openings. A numerically more accurate higher order can still fail passivity.
Each attempted order and its actual refusal are retained in the CLI sidecar.

Large arrays retain every mutual response while the original shared vector
fitter eliminates each response's local numerator with Householder QR before
solving one common denominator. Every array uses that same projected solver,
including its local rank admission. For 17 openings, 65 samples and
orders 4/8/12, the conservative source/fitting scratch estimate is 22,508,704
bytes, down from 1,544,916,032 bytes for the previous dense global solve. The
order-12 realization has 204 real states, so it still needs an explicitly
sufficient runtime-state allowance. Passivity repair receives the remaining working allowance for its model,
certification, constraint and KKT storage and refuses excess storage without dropping
constraints. These capacity changes do not relax accuracy, causality or
passivity requirements.

The physical image is uniform circular apertures in an infinite rigid baffle
and homogeneous inviscid exterior gas. The force/velocity matrix is converted
to pressure/volume-flow using both aperture areas and the phasor convention
exactly once. Existing arbitrary gas mixtures and thermochemical pack-backed
gas records are resolved before building, including their original source
receipts. Reaction and observation share that exact exterior.

The finite rational approximation represents propagation without extracting
separate delays. Its sampled errors do not bound every intermediate frequency,
closed-loop sound error, or instrument accuracy. The causality result is a
consistency test under the explicitly declared tail/data assumptions; a large
tail allowance weakens that test. Continuous passivity keeps the original
certificate class, while interval-mean ZOH admission is explicitly `GridOnly`,
not a discrete storage certificate. This does not complete the plan's delayed
radiation, passivity-repair, or measured-instrument qualification work.

## Select the source

Flue performance input uses the optional existing radiation position:

```text
radiation mimo flue-mutual.mimo-radiation 32 1048576 257
```

The full example is [flue-mimo.performance](flue-mimo.performance), with
[flue-mutual.mimo-radiation](flue-mutual.mimo-radiation). Existing `prescribed`,
`compact-baffled` and omitted selections retain their original behavior.

Material-modal valves place the selection immediately after their `opening`
rows, before `valve`. Each selected acoustic side has one group:

```text
radiation mimo 1
radiation_group downstream mutual.radiation 256 8388608 65
```

[paired-mutual.performance](material_modal_valve/paired-mutual.performance)
and [mutual.radiation](material_modal_valve/mutual.radiation) form a complete
two-port example. Every observed side must have a group. Within that side,
the source matrix follows the performance's opening order, including when
upstream/downstream openings are interleaved globally. A second group on the
same side, an unobserved selected side, or an observed side without a group
refuses. The upstream/downstream groups are independent; this does not model
mutual exterior coupling across the moving valve.

Only each group's explicitly observed memoryless endpoints receive its coupled
load. Every unselected outlet and unobserved acoustic side retains its authored
reflection or reactive load, including its original state. Interior material
walls, supply passage, gas regions, loss subdivisions and fractional travel
also retain their existing owners. MIMO is not added on top of an existing
reactive terminal. Material-modal omission and explicit `radiation compact-baffled`
keep the original compact reaction path.

The four selector fields are a caller-relative artifact path, maximum real
state dimension, maximum terminal memory bytes, and sampled-passivity grid
points. Paths cannot be absolute or contain empty, `.` or `..` components,
backslashes or colons. Native reads are capped at 8 MiB. State caps are 0–256,
memory is 1–67108864 bytes, and grid points are 17–65537. A zero state cap is
useful only for an admitted static matrix; a dynamic model still needs its
actual complete state dimension. The runtime applies the selected tighter
budget without truncation or fallback.

## Original artifact grammar

Use UTF-8, one exact ordered record per LF line, and a final LF. Empty records,
comments, CRLF, extra fields and trailing records refuse. All values must be
finite. This two-port example contains one real pole:

```text
frankensim-mimo-radiation-v1
declaration 1111111111111111111111111111111111111111111111111111111111111111
realization 2222222222222222222222222222222222222222222222222222222222222222
exterior 293.15 101325 0
band 1 1000000
apertures 2
aperture 0 0 0.007
aperture 0.03 0 0.007
reciprocal true
poles 1
pole real -7000
residue 4000000 0 1000000 0 1000000 0 4000000 0
direct 2000 300 300 2000
improper 0 0 0 0
```

`declaration` and `realization` retain the original nonzero 32-byte geometry/
material and rational-bake identities, written as 64 hexadecimal characters.
The importer separately hashes the **complete artifact bytes**, using BLAKE3
domain `org.frankensim.fs-couple.mimo-radiation.v1`. The sidecar names this
`artifact_blake3`; it is distinct from both declared identities. Two byte
encodings of identical numeric coefficients have distinct artifact hashes.
None of these fields alone proves that an operator accurately represents a
physical specimen.

The `exterior` record uses the existing complete gas grammar: `T P RH`, named
`mixture T P mole|mass COUNT SPECIES FRACTION ...`, or the independently bound
thermochemical form below. Pressure is never inherited from the instrument.
The artifact's resolved gas must exactly match the common exterior used for
reaction and observation. Interior regions continue using their original gas.

`band` is the original physical frequency interval in Hz, positive and
increasing, with its upper edge covering the actual mechanics Nyquist limit.
`apertures` is 1–128. Each aperture contains its x/y center and radius in metres
in the common z=0 infinite rigid baffle, in exact matrix order. The runtime
checks these values against the actual degree-one duct endpoints.

`poles` is the number of common pole terms, 0–256. Each term is either
`pole real REAL_LOCATION` or `pole pair REAL_LOCATION POSITIVE_IMAGINARY_LOCATION`,
immediately followed by `residue` and exactly `N*N` complex `RE IM` pairs in
row-major order. A pair stores only the positive-imaginary pole and its matrix;
the conjugate pole/matrix are implicit in the shared rational owner. Real
poles require real residue entries at runtime. The original real state count
is `N * (REAL_TERMS + 2*PAIR_TERMS)`, capped at 256 before matrix allocation.

`direct` and `improper` each contain exactly `N*N` real row-major values.
The source convention is impedance in Pa s/m³ evaluated at `s=+i omega`.
Acoustic `exp(-i omega t)` input must already have been conjugated by its source
producer. This host admits proper, undelayed models: every improper coefficient
must be zero. The existing runtime checks stable poles, actual reciprocal
entries, continuous and sampled passivity, finite dimensions, explicit storage
budget, source clock, geometry and gas. Parsing alone does not imply admission.

## Thermochemical source authority

Thermochemical species names require the existing explicit source registry.
Immediately after `exterior thermochemical ...`, the artifact must also carry
its **original** canonical source receipts and complete gas snapshot:

```text
gas_sources COUNT
gas_source SPECIES_ID SPECIES_PACK_HASH NASA9_PACK_HASH TRANSPORT_PACK_HASH
gas_state T P RHO C MU KAPPA GAMMA R_SPECIFIC CP PRANDTL RHO_C WATER_MOLE_FRACTION
```

Repeat `gas_source` exactly `COUNT` times, in the existing gas owner's canonical
species order; count is bounded at 128. All three hashes are the complete
normalized artifact identities, not filenames. `gas_state` contains exactly
12 SI `f64` values in the shown order. Write enough digits to reproduce the
original floating-point values exactly. The importer uses the shared gas
decoder, then compares every receipt and every snapshot field bit-for-bit.
It retains that original resolved `GasRecord`; the snapshot never overwrites
or reconstructs the state.

Replacing a NASA-9 or transport pack under the same species ID therefore
refuses with unchanged radiation bytes. The snapshot also detects numerical
changes from source selection or water policy when pack identities alone are
insufficient. These extra records are required only for thermochemical gas;
legacy and named-mixture forms remain self-contained under their existing
material models.

## Library and physical scope

Both native loaders expose `from_bytes_with_all_sources`, with separate
material-pack and radiation-byte callbacks plus the optional existing
thermochemical registry. Library admission performs no filesystem access.
Existing constructors continue accepting existing source files; they refuse
a requested MIMO artifact when no radiation loader was supplied.

Signed coupled terminal work remains one external exchange account in the
original work gate. It is not dissipation or a certified rational stored-energy
coordinate. Runtime callbacks preserve the original wave, material and pole
histories and retain their existing atomic cancellation/refusal rules.
The current coplanar infinite-baffle observer and sampled `GridOnly` passivity
report remain diagnostics. This import path does not add arbitrary unbaffled
surfaces, delayed operators, automatic refitting, evolving source geometry or
gas, cross-side coupling, measured flute/reed validation, or M2.5-W3c closure.

Focused source, native-performer and actual CLI regressions retain the original
gas/material histories, exercise off-diagonal feedback and microphone changes,
compare callback partitions and PCM, and refuse stale gas, reordered geometry,
invalid bands/resources, source-pack replacements and cancellation. Native
execution is required before treating these authored tests as passed evidence.
