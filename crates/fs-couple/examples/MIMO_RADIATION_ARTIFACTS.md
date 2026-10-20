# Source-bound mutual radiation artifacts

The native flute and material-modal valve formats can select the same original
multiport impedance artifact. Each source goes through
`bernoulli_aperture::radiation::environment::multiport::input`, then the existing
`MimoEnvironment` and `WaveguideNetwork` terminal host. Coupled outlet pressure
returns to the live jet or moving valve. The microphone observes the same
ordered apertures in the artifact's exact exterior gas.

This is an importer for externally produced coefficients. It does not calculate
a new radiation model when geometry, materials or atmosphere change, or attest
a declared offline bake. Supply a corresponding new artifact when the changed
source requires one. A mismatch refuses construction before audio is produced.
The included operators are **manufactured integration data**, with no measured
instrument or broadband radiation accuracy claim.

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
