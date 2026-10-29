# Molecular gas data in reed feedback and physical sound

`reed-molecular-synthetic.performance` drives a 15 cm cylindrical bore with a
pressure-controlled reed and observes its flanged outlet at one metre. Its
N2/O2 vibrational reservoirs are **synthetic demonstration coefficients**, not
measured nitrogen/oxygen data or a calibrated instrument reconstruction.

From a bootstrapped repository, use the existing 48 kHz CLI:

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/reed-molecular-synthetic.performance \
  /tmp/reed-molecular-synthetic.wav --block 37
```

The WAV and its `.provenance.json` sidecar must name fresh output paths. The
source declares 9,606 samples at 48 kHz. Its `audio ... 2` full scale is 2 Pa
for PCM encoding; it does not change the physical pressure or normalize peaks.
The targeted library tests use an 8 kHz version with the same duration and
schedule; `music_render wind` does not upsample it.

## Native declarations

Insert these three records after geometry and before `compile_limits`:

```text
molecular <shared fs-material relaxation declaration>
molecular_resolution CONE_SLICES MAX_COMPOSITION_BLOCKS
molecular_wall diffusive POLES OMEGA_MIN OMEGA_MAX MAX_RELATIVE_ERROR
```

The example uses `molecular_resolution 8 4096` and
`molecular_wall diffusive 32 0.0001 160000 0.05`. Frequencies here are **rad/s**.
The band must cover the positive near-DC probe at 0.001 rad/s and every line
sample through Nyquist (`pi * sample_rate_hz`). Diffusion admission, cone
resolution and the converged FIR history remain independent checks.
`molecular_wall none` selects adiabatic slip walls instead.

The shared law syntax is `declared N`, followed by N `SPECIES PROCESS_COUNT`
entries. Each process pairs a reservoir (`vibration MODES THETA_K DEGENERACY ...`
or `rotation DOF`) with a rate (`power P_TAU_REF T_REF EXPONENT T_MIN T_MAX`
or the shared `partners` record). The alternative `iso9613-terrestrial-air`
retains its existing composition and temperature restrictions. Every present
gas species needs a declaration, including an explicit zero process count
when appropriate; the solver does not invent missing molecular properties.

## Changes to make in data

Change `ambient mixture 293.15 101325 mole 2 N2 0.7 O2 0.3` to change temperature,
pressure or composition. For example, use N2/O2 fractions `0.6 0.4`, or pressure
`80000`. The same source supplies gas properties, reed flow, bore dispersion,
loss and default exterior loading. Additional species require their own gas
source and relaxation declaration. Source-domain limits still apply.

Change the first `power 2.0265` to `power 0.20265` to change only the N2
reservoir's pressure-time coefficient by tenfold. Geometry, gas composition,
reed mechanics and pressure gestures stay identical, while the original
constitutive law changes bore feedback and physical outlet pressure.

`molecular` describes gas **inside the bore**. `observer_law none` explicitly
omits atmospheric absorption between the outlet and listener; it preserves
the receiver's propagation delay. Select an appropriate existing observer law
for that exterior path. An optional `exterior` record after `baffled_outlet`
binds both terminal radiation and receiver, at the common static pressure.

This route retains frozen homogeneous linear acoustics and the existing reed,
FIR and baffled receiver approximations. Explicit segment gases, axial or
steady-generated gas fields, and tone holes refuse; their molecular coupling
is not supplied by this harmonic owner. It does not model evolving species,
mean flow or compliant bore walls. Ordinary inputs without `molecular` keep
their existing propagation path.
