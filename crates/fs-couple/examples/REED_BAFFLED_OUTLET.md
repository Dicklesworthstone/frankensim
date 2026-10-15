# Physical terminal sound from regional gas reed performances

The `frankensim-reed-performance-v1` source can select terminal-only Rayleigh
pressure rather than the historical bore-pressure/compact-jet proxy. Replace
`listener_m` with the explicit physical observer record:

```text
baffled_outlet 0 0 0.5 4 16 1000
termination flanged
```

The fields are X/Y/Z [m] in the actual outlet's baffle frame, radial rings,
angular points, and the declared observation-band upper limit [Hz]. The existing
receiver requires z >= 0.05 m, 1..=128 rings, 8..=512 angular points, and a positive
band no greater than one tenth of the mechanical rate. Its geometry, quadrature
and history limits remain active. This declaration does not add a low-pass filter.
Closed, ideal-open and unflanged terminations refuse instead of being replaced.

The complete `reed-baffled-mixture.performance` example contains two prescribed
gas mixtures. Its geometry, material inputs and pressure schedule follow the
regional reed example. The output scale is explicitly 2 Pa; it is not a gain or
peak normalization. From a bootstrapped repository:

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/reed-baffled-mixture.performance \
  /tmp/reed-baffled-mixture.wav --block 37
```

## Independent internal and exterior atmospheres

An optional `exterior` record immediately after `baffled_outlet` selects the
homogeneous gas used by BOTH the terminal radiation load and the microphone:

```text
baffled_outlet 0 0 0.5 4 16 1000
exterior mixture 300 100000 mass 2 N2 0.7 O2 0.3
termination flanged
```

The same material-owned syntax supports mole fractions or legacy `exterior T P RH`.
No species constants or constitutive equations are duplicated. The exterior must
have the common static pressure of the duct. Composition and temperature may
otherwise differ within the existing material admission domain. Omission uses
the actual final segment's gas, NOT ambient or the inlet; an explicit identical
exterior retains the original calculation. `reed-baffled-exterior.performance`
is a complete comparison whose only difference is the exterior record.

Internal density, viscosity, conductivity, sound speed and source impedance
remain those of the original segment gases. Exterior changes affect radiation
loading, the resulting reed feedback, microphone weighting and exterior travel;
they do not replace bore material or alter geometric length. Unsupported
exposed-hole/two-medium chimney combinations refuse through the existing owner.
An exterior record with the legacy proxy, duplicates, invalid fractions or an
unequal static pressure also refuses before CLI output creation.

## One physical source and its selected observation

`ReedPerformance` hosts the original `ReedBoreVoice`, pressure scheduler and
`RegionalDuct::baffled_outlet`. The latter samples the original TMM terminal-flow
transfer and existing Rayleigh observer. Both feedback and observation use one
bound regional duct, flanged load and admitted FFT grid. The original 5 Pa prime
is supplied exactly once at source sample -1 to each history. The physical
observer REPLACES the bore/jet proxy; the two pressure images are never summed.

The terminal adapter observes every accepted outgoing wave by asking the existing
scheduler for one mechanical sample at a time. It does not implement another
mechanical solver or control schedule. There is no real-time performance claim.
Legacy inputs keep their original block path. Invalid request sizes/horizons
refuse before mutation; late physics failures may leave a partial callback and
poison the complete source. Discard a failed callback rather than resuming only
one of its histories.

`pressure_block` retains typed `ReedPerformanceError::Outlet` failures.
`last_outlet_error()` retains the original observer failure even when using the
older shared `PressureRenderer` string diagnostic envelope. `outlet_receiver()`
reports the actual exterior gas, geometry, Fourier grid and delay support.

Standalone wind and `ensemble --reed` consume this same finite pressure source.
Their sidecars describe the selected observation, ambient and segment gases,
actual receiver medium and physical delays. Explicit `--decimate` retains the
existing output conversion. Axial and exterior travel are NOT filter latency
and are not removed by ensemble alignment. Source PCM scales are not mixing gains.

## Scope and focused checks

These are stationary, prescribed gases, with no composition transport, live
medium replacement, thermal equilibration or mean flow. The observation is
terminal-only: side-hole, jet, structural and room radiation are absent.
Finite-history Fourier sampling, fractional-delay and backward-difference
approximations remain those of the existing owners. No certified broadband
radiation match, measured timbre or instrument calibration is implied.

The existing instrument workflow already includes both test targets:

```bash
cargo test -p fs-couple --lib render::schedule::reed:: -- --test-threads=1
cargo test -p fs-couple --test reed_gas_mixture -- --test-threads=1
```

Regressions cover direct owner equality, massive/massless reed-state parity,
physical delays, callback partitioning, explicit exterior feedback, exact matched
exterior parity, material/interface refusal, native/decimated CLI PCM and ensemble
summation. Authored tests do not themselves establish a passing native run.
