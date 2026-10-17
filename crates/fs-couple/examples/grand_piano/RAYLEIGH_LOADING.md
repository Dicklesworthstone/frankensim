# Radiation feedback in the primary piano render

`--rayleigh-loading sides,min_hz,max_hz,samples` lets the acoustic pressure
react on the same strings and soundboard that the primary Rayleigh microphones
observe. It changes physical motion, resonance and decay through the existing
passive radiation owner. It is available for a flat geometric board with
`--render`, including the reconstructed Steinway-D preset:

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --render piano.wav --duration 2 \
  --rayleigh-loading 2,20,1000,129
```

The four values specify the number of radiating half-spaces, the positive
frequency-band endpoints in Hz, and an odd count of 33 through 257 frequency
samples (linear for a complete load, logarithmic with `--radiation-ports`).
There is no inferred side count or separate fit-order knob.
Without this option, the existing one-way receiver path is unchanged.

`--radiation-ports N`, with an explicit ceiling from 1 through 128, connects the
existing projected acoustic runtime to this load. It can attempt a complete
board image of up to 128 coordinates while fitting fewer radiating combinations.
The complete physical impedance is still sampled; every mechanical coordinate
and microphone input remains. Geometry alone selects the orthonormal acoustic
rows before any frequency response or microphone position is consulted. The
lifted fit must pass the existing complex and resistance limits against every
original full matrix, including the held-out frequencies. Insufficient rank
refuses without increasing the requested ceiling. The report distinguishes
board coordinates, acoustic ports, surface-field omission and complete-load
errors. This option requires `--rayleigh-loading`.

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

The shared passive fitter uses even frequency samples for its fit and odd
samples for withheld validation. Complete complex peak/RMS errors must meet
5%/2%; resistance errors must meet 10%/5%. No failed load is replaced by zeros,
no signed cross term is removed, and no acceptance limit is relaxed.

A complete load (no `--radiation-ports`) is sampled LINEARLY in frequency and
fitted with a dictionary derived from the board's own geometry. The Rayleigh
kernel is a superposition of delays up to `tau = D/c`, `D` the largest distance
between two area cells, so the impedance varies on a scale of `1/tau` in
frequency at the top of the band as much as at the bottom. In-band resonators
are spaced `c/(6D)` apart with half-power bandwidth twice that, continued six
spacings past the band top, with lossless storage just above the band (the
frequency-dependent added mass) and the shared sub- and high-band terms. The
training samples must be no further apart than that spacing, or the fit
refuses. The shared fixed 32-term dictionary, log-spaced at 12% damping, cannot
follow that detail over a wide band: on the 2.4 kHz production board it missed
by 21%/17% (complex peak/RMS) even fully converged, and the 32-port projection
of the same load by 13%/11%. The resolved fit of the same board passes at
1.9%/1.5% (resistance 1.9%/1.5%) with 12,912 acoustic oscillators (D = 2.355 m,
24.3 Hz spacing, 257 samples over 90-2400 Hz). An explicit projection keeps the
logarithmic grid and the shared fixed dictionary. The resulting model enters `Instrument::configure_radiation`
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

The direct piano path admits up to 128 complete board modes (every board the
bank admits), 16,384 area samples and 4 trillion cell-mode-frequency visits.
The fitter and runtime model admit 128 complete velocity ports and 16,384
acoustic oscillators. The runtime's collective exchange bounds its largest
rotation per half step, `sigma_max(L) dt/2 <= 0.125`; the 2.4 kHz board's
complete load turns by 0.029. Each oscillator's free transition is independent,
so they are prepared in groups within the modal time owner's own ceiling.
The opt-in projection retains the complete board and attaches a smaller
acoustic load through the existing orthonormal port map. Its full-matrix error
checks include every omitted acoustic direction; a small surface-field error
alone cannot establish an adequate load. All cell, work, clock, resolution and
fit limits remain. The kernel caches the spatial pair distances
and projects each frequency in `O(cells^2*modes + cells*modes^2)` work; the
frequencies are independent and are evaluated in parallel (each matrix is
bit-identical to a serial sweep). No acoustic port reduction or automatic
spatial coarsening is performed. Preparation is offline work, with no
real-time claim: on the 2.4 kHz production board (9,432 cells, 128 modes, 257
frequencies) the sweep takes about 1.5 minutes on 28 cores and the fit about
12 minutes on one.

The default eight-division preset produces 3,918 samples. Its largest
equivalent-disc cell admits approximately 1,097 Hz at 343 m/s under the existing
six-cells-per-wavelength floor, so the example's 1,000 Hz endpoint passes this
area-resolution check. A 2.4 kHz band needs a finer board, for example
`piano_geometry --divisions 12 --rib-normal-step-m 0.03` (9,432 samples, 6.32
cells per wavelength), supplied with `--board-geometry`. These geometry/work checks do not establish
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
