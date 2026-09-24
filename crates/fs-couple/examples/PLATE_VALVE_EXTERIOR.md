# Exterior outlet pressure from supplied valve performances

A supplied stationary receiver can observe the **actual outward terminal volume
flow** through a uniform circular outlet in an infinite rigid baffle:

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-exterior.performance \
  /tmp/plate-valve-exterior.wav --decimate --block 37
```

This example preserves the mesh, contact, material, pressure program and tube of
`plate-valve.performance`. It explicitly selects a receiver 0.2 m above the
outlet, 8 equal-area radial rings, 32 azimuthal samples and a 4 kHz resolution
band. Its declared PCM full scale is 0.01 Pa, suitable for examining this small
illustrative excitation; that changes encoding only, never physical pressure.
The prior source file, primitive reed input and internal observations are unchanged.

The new observation record is:

```text
observation baffled-outlet X_M Y_M Z_M RADIAL_RINGS ANGULAR_POINTS MAX_FREQUENCY_HZ
```

Coordinates are relative to the outlet center, with the mouth in z=0 and outward
normal +z. Radius, density and sound speed come from the **actual tube and gas**;
there is no separate tunable acoustic radius, area or volume-flow-to-Pa gain.
Uniform normal velocity is Q/(pi a²), where Q is the existing waveguide's accepted
terminal flow, not its pressure. A perfectly closed termination has zero flow
and no output from this observer even if its internal pressure is nonzero.

The shared `pcm_wav::baffled::BaffledPressure` integrates the time-domain Rayleigh
half-space kernel over positive-area disk samples. Distance, travel time and
interference are retained, including off-axis receivers. This is the same owner
used by the piano's baffled-surface microphone; its existing loaded-mode projection
and numerical stepping order are unchanged. No new structural or flow integrator
is introduced. `ApertureObservation::TubeBaffled` selects it for a library tube;
network-node observations remain internal and cannot masquerade as this outlet.

Admission requires z >= 0.05 m, 1..128 rings and 8..512 azimuthal points. The
positive declared band is at most one tenth of the mechanical sample rate, and
the disk cell's conservative phase-span estimate must be at most pi/4 in that
band. These are **resolution checks, not a continuum error enclosure or filter**.
Nonlinear source content may exceed that band; quantitative work needs temporal,
spatial and source-physics refinement. Acceleration is differenced from consecutive
flow samples, centered half a mechanical sample earlier; fractional propagation
uses linear interpolation with less than one-sample wavefront uncertainty.

The receiver runs at every mechanical sample before the existing explicit output
decimator. Source pressure retains its physical travel delay. Ensemble alignment
aligns only decimator latency; it must not cancel geometric propagation. The
sidecar records receiver geometry, quadrature, declared band and propagation delay
in **mechanical samples**, separately from output filter delay. `ensemble --valve`
accepts the same file and retains this distinction. Supply a compatible common
observation/time origin when composing different sources. Source PCM scales still
never become ensemble gains.

**This is one-way prescribed-flow radiation.** The independently supplied terminal
reflection still determines the mechanical solution. The observer does not add a
matched radiation impedance, subtract energy from mechanics, or make the total
source/radiator acoustically self-consistent. It assumes a uniform outlet profile,
a stationary infinite rigid baffle and homogeneous exterior air, not an unbaffled
bell, reed-jet noise, body/room scattering or empirical microphone calibration.
Moving the receiver leaves valve, contact, material and wave energy unchanged.

Keep the complete `AperturePerformance` through cancellation/resume to retain both
wave and receiver history. A physical or acoustic callback failure poisons the
pressure adapter, rather than resuming mismatched histories. Finite windows do
not synthesize padding or flush propagation tails beyond their declared duration.

## Optional compact radiation feedback

The fixed-reflection example above remains unchanged. To let exterior radiation
react on the valve and returning tube waves, select a **physical terminal load**
instead of the numeric reflection:

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-radiating.performance \
  /tmp/plate-valve-radiating.wav --decimate --block 37
```

The new alternative `tube` record is:

```text
tube LENGTH_M RADIUS_M baffled-low-ka MAX_FREQUENCY_HZ LENGTH_ERROR_M MAX_WAVE_BYTES
```

The committed example uses `tube 0.25 0.007 baffled-low-ka 1500 0.002 1048576`
and a 1.5 kHz receiver band. Everything else—the original mesh, material-memory
law, contact profile, pressure phrase, receiver position and PCM scale—is the
same as the one-way example. The two terminal forms are mutually exclusive;
there is no simultaneous numeric reflectance or separately tuned end correction.
The tube length is the physical interior length, not a radiation-extended length.

`BaffledRadiationLoad` derives the load from the actual radius a, gas density rho
and sound speed c. With Zc=rho*c/(pi*a²), alpha=8/(3*pi), its impedance is

```text
Z(s) = R s/(s+p)
R = 2 alpha² Zc
p = 2 alpha c/a
```

These coefficients match the leading resistance `(ka)²/2` and inertance
`alpha*ka` of the uniform baffled piston. They are not fitted to the generated
waveform. The existing passive relaxation-impedance owner advances its retained
flow state and returns the reflected wave. The existing characteristic network
and nonlinear valve solve provide the rest of the reciprocal dynamics. No new
integrator, delayed force correction, extra damping or second propagation path
is introduced. Near-field energy belongs to the load's stored energy; positive
radiation loss is counted once in the joint work balance, not confused with all
terminal port work. Internal pressure observations can select the same load.

This is a **compact approximation**, not an exact broadband radiation boundary.
Admission requires ka<=0.5 and dt*fmax<=0.1. Before plate reduction, 32 positive
frequencies are compared with the existing full Rayleigh disk oracle. Both the
analogue response and the actual bilinear response must meet fixed 5% relative
complex-impedance AND resistance checks. Resistance has its own check because
small radiation losses can be hidden by much larger reactive pressure. Failed
checks refuse the source; no coefficients or tolerances are retuned. The oracle
has finite spatial quadrature, and the checks are sampled estimates, not a
continuum error bound. The pressure phrase is not automatically band-limited;
nonlinear out-of-band content still needs refinement and physical assessment.

The receiver's declared band must not exceed the load's. It observes the same
accepted terminal volume flow that participates in the load equation, so the
radiation-modified mechanics reach the exterior waveform. The receiver itself
adds no second energy sink. Moving only the receiver retains exactly the same
loaded mechanics. Its numerical Rayleigh transfer and the compact load are
compatible approximations to the same baffled geometry, not a claim of exact
all-frequency source/radiator consistency. Unbaffled bells, viscothermal tube
loss, body scattering, room acoustics and measured calibration remain outside
this path.

The file uses a single existing network section and its passive terminal; this
requires slightly more admitted wave payload than a bare memoryless tube. Source
metadata retains the load coefficients, checked discrepancies and declared band.
`ensemble --valve` accepts the same file: each part keeps its complete material,
load, propagation and receiver history and its source clock. Geometric receiver
delay is not cancelled by filter-latency alignment, and source PCM scales remain
encoding declarations rather than hidden ensemble gains. Existing output-file,
finite-window, cancellation and failed-callback rules are unchanged.
