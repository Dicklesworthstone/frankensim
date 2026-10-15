# TMM reed performances in named gas mixtures

Both wind performance formats now accept the same material-owned ambient record:

```text
ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3
```

This guide covers `frankensim-reed-performance-v1`: the existing primitive reed
and transfer-matrix duct path. The geometry-derived plate-valve network is
covered by [GAS_MIXTURES.md](GAS_MIXTURES.md). Neither parser contains gas mixing
equations. Both call `fs_material::gas_species::input::read_gas_state`, which
uses the existing composition and constitutive owners.

`mole` and `mass` select fractions, not percentages. The named model admits N2,
O2, Ar, CO2 and H2O, with the aliases and temperature/pressure/saturation limits
in that guide. The old `ambient T P RH` syntax and dry-air limit are unchanged.
Unknown species, malformed fractions and unsupported states refuse. A mixture
is not implicitly eligible for air-only absorption correlations.

## Compare actual physical outputs

The following complete inputs differ in their ambient record ONLY:

| File | Mole fractions | Pressure |
| --- | --- | --- |
| `reed-mixture.performance` | N2 0.7, O2 0.3 | 100000 Pa |
| `reed-mixture-co2.performance` | N2 0.4, O2 0.2, CO2 0.4 | 100000 Pa |
| `reed-mixture-low-pressure.performance` | N2 0.7, O2 0.3 | 50000 Pa |

Temperature is 293.15 K throughout. The shared cylinder radius is explicitly
0.0015 m, so the slower mixture can meet the existing unflanged-radiation band
limit at 48 kHz. The parser does not resize a bore or remove radiation to make a
new gas admissible. These are short synthetic pressure phrases, not measured
instrument reconstructions or calibrated exterior recordings.

```bash
for case in mixture mixture-co2 mixture-low-pressure; do
  cargo run --release -p fs-couple --bin music_render -- \
    wind "crates/fs-couple/examples/reed-${case}.performance" \
    "/tmp/reed-${case}.wav" --block 37
done
```

The admitted gas supplies the same existing reed, bore viscothermal response,
termination and compact-jet observation. It is available through
`ReedPerformance::ambient_gas()`. At fixed temperature and composition, the
ideal-gas sound speed does not change with pressure; density, impedance and
reed/bore loading do. No pitch shift, output normalization or source-gain
correction is inserted. `ensemble --reed` consumes the same finite source;
source PCM full scales remain encoding declarations, not mixing gains.

Focused native commands, from a bootstrapped repository:

```bash
cargo test -p fs-material --lib gas_species::input::tests
cargo test -p fs-couple --lib render::schedule::reed::
cargo test -p fs-couple --lib driving_point::regional::
cargo test -p fs-couple --test reed_gas_mixture
```

The last target exercises the actual executable, native and explicitly
decimated output, complete WAV equality with the existing PCM owner, all five
species on a mass basis, two-mixture ensembles and refusal before file creation.
These tests do not certify experimental timbre or remove the original TMM,
finite-history, reed and radiation approximations.

## Different gases along the same duct

An axial or tone-hole record can append `gas` followed by the shared material
syntax. For example, the complete `reed-regional-mixture.performance` supplies:

```text
segments 2
cylinder 0.0015 0.2
cylinder 0.0015 0.3 gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4
```

The first segment inherits `ambient`; the second supplies its own mixture.
Omitted `gas` always inherits ambient, NOT the preceding segment. The alternative
`gas 303.15 0.5` means moist air at 303.15 K and relative humidity 0.5 at ambient
pressure. `mass` and `mole` work on segment declarations as on ambient. All
connected segments must have the same static pressure as ambient: unequal
pressures would require a mean-flow model, which this renderer does not supply.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/reed-regional-mixture.performance \
  /tmp/reed-regional-mixture.wav --block 37
```

`ReedPerformance::segment_gases()` exposes the actual media in authored order.
The inlet gas supplies Bernoulli density, the characteristic port and the
existing compact-jet proxy. Each constant-medium run uses the original
`fs-duct` cylinder/cone/wall/chimney laws; its downstream load is the following
run's input impedance. The final gas supplies terminal radiation, and physical
travel time sums local length/sound-speed contributions. The same FIR realization,
reed initialization, sample stepper, pressure schedule, PCM encoder and ensemble
path are retained. No averaging, retuning, gain correction or alternate solver
is introduced; adjacent identical gases use the original uniform TMM call.

All open/vented holes must lie in one contiguous constant-medium run so that the
existing mutual-hole corrections remain complete. Open holes across gas-region
boundaries refuse; their cross-region near-field law is not supplied. A hole's
assigned gas describes its compact junction and chimney together, not a separate
exterior gas. Closed holes can lie in different regions. The 64-record cap counts
holes as well as axial segments. Geometry and radiation-band limits remain active.

These are prescribed, stationary gas regions. Neither wind route evolves
composition during playback, solves mean flow/mixing or thermal equilibration,
or replaces a live medium while transferring its stored acoustic state.
The reed observation remains bore pressure plus a compact-jet proxy, not a
calibrated exterior microphone. Regional tests cover direct construction,
uniform bit parity, downstream changes in massive reed motion, local radiation
admission, and native/decimated wind and ensemble PCM through the actual CLI.
