# Radiation feedback in the primary piano render

`--rayleigh-loading sides,min_hz,max_hz,samples` lets the acoustic pressure
react on the same strings and soundboard that the primary Rayleigh microphones
observe. It changes physical motion, resonance and decay through the existing
passive radiation owner. It is available for a flat geometric board with
`--render`, including the reconstructed Steinway-D preset:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --render piano.wav --duration 2 \
  --rayleigh-loading 2,20,1000,65
```

The four values specify the number of radiating half-spaces, the positive
frequency-band endpoints in Hz, and an odd count of 33 through 257 logarithmic
frequency samples. There is no inferred side count or separate fit-order knob.
Without this option, the existing one-way receiver path is unchanged.

## Physical model

The load uses the same planar, infinite rigid baffle assumption as the primary
pressure observer. `sides=1` loads one half-space. `sides=2` adds two independent,
identical half-spaces on opposite faces, doubling the reaction impedance. It
does not double either microphone's pressure or duplicate its mechanical step.
Two sides cannot accompany `--cavity`, whose explicit air volume already
occupies the underside. One side can compose with that existing cavity owner.
Crowned/curved boards are refused. Cabinet edges, lid scattering and a finite
body require the existing `piano_exterior render-loaded` path and supplied body
geometry; this option does not infer them from a flat panel.

Every supplied area sample and retained board mode participates. The bank's
actual string-mass-loaded coordinate map projects the source shapes once,
exactly as it does for the microphones. The shared BEM Rayleigh kernel forms
the complete reciprocal load, including signed cross-mode terms. Under
`exp(-i omega t)` its added-mass reactance is negative imaginary. Each positive
area is included once in the pressure integral and once in its work-conjugate
force projection. Equivalent-disc self integrals treat each cell's local
singularity; other cells use the existing outgoing Green kernel.

The unchanged shared passive fitter uses even frequency samples for its fixed
dictionary fit and odd samples for withheld validation. Complete complex
peak/RMS errors must meet 5%/2%; resistance errors must meet 10%/5%. No failed
load is replaced by zeros, no signed cross term is removed, and no acceptance
limit is relaxed. The resulting model enters `Instrument::configure_radiation`
before excitation. Its existing collective exchange, acoustic storage, damping
and rollback participate in the mechanical energy balance. Both microphones
then observe the reacted motion using their existing pressure conversion.

Declared board damping should represent intrinsic structural loss when adding
this load. Damping identified from total measured decay may already include
radiation; adding an acoustic load to that total would count the air loss again.

## Resolution and scope

The complete sampled band must include every loaded board reference oscillator.
Those are the split step's diagonal references, not all coupled string/board
poles. The upper endpoint must be below 45% of the output sample rate, and the
shared fitter's poles through eight times that endpoint must satisfy the
existing mechanical-rate guard. At 48 kHz with four substeps, the upper endpoint
must therefore be below 10.8 kHz. Broadband excitation can still reach outside
the sampled band; the fitted load has no validated accuracy there.

The initial capability admits up to 32 complete board modes, 8,192 area samples
and 64 billion cell-mode-frequency visits. It caches the spatial pair distances
and projects each frequency in `O(cells^2*modes + cells*modes^2)` work. No acoustic
port reduction or automatic spatial coarsening is performed. Preparation is
offline work, with no real-time claim.

The default eight-division preset produces 3,918 samples. Even its worst admitted
32-mode, 65-frequency sweep uses about 31.93 billion visits, below the ceiling.
Its largest equivalent-disc cell admits approximately 1,097 Hz at 343 m/s under
the existing six-cells-per-wavelength floor, so the example's 1,000 Hz endpoint
passes this area-resolution check. These geometry/work checks do not establish
that a particular full preset fit or WAV render will succeed. The complete
matrix fit, loaded-reference coverage and runtime admissions still apply.

Equivalent-disc cells per wavelength is an area-based diagnostic, not a bound
on triangle aspect ratio, quadrature error or acoustic accuracy between sampled
frequencies. Acoustic refinement also increases the number of retained cells
and can exhaust the explicit budget. The printed report gives the actual
geometry count, work, sampled band, reference range, pole count and fit errors.

## Focused checks

The numerical-owner tests compare a uniform circular field with the unchanged
piston owner and independent small-wavenumber resistance/added-mass formulas.
They also compare modal resistance with independently integrated outgoing
hemisphere pressure power and exercise signed modal coordinate changes.
The piano tests pass real loaded-bank shapes through the actual passive fitter,
check two-half-space scaling and refusals, and evolve two strings, a two-mode
board and the actual radiation runtime through half a second of free decay.
That last test requires positive radiative loss, reactive acoustic storage and
closure of the combined physical energy ledger.

```sh
cargo test -p fs-bem --lib rayleigh_modal
cargo test --release -p fs-couple --example grand_piano rayleigh_loading
```

The focused native owner tests establish those numerical and physical coupling
identities. They do not constitute measured-instrument calibration, an executed
full preset render, or soundboard mesh convergence evidence.
