# Full physical preparation in finite-body piano playback

Both `piano_exterior render` and `render-loaded` accept the physical controls
already used by `grand_piano`. The former remains one-way; the latter retains
passive radiation feedback. No new piano engine or acoustic solver is selected.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded settled.fss strings.csv body.obj acoustics.fspe loaded.wav 6 \
  performance.mid --modes 128 --substeps 8 \
  --hammers materials.fsh --hammer-footprints faces.fshp --dampers pads.fspd
```

Each control is optional. Old positional MIDI invocations and the default
4-substep, 24-partial source-shank image remain available unchanged. Without a
score the existing A4 demonstration is used. Outputs remain 48 kHz physical-Pa
PCM, without automatic normalization. The render duration remains 0.05–60 s.

`--modes 1..512` changes the per-string retention ceiling, including unplayed
sympathetic strings and duplex segments. Bass strings are no longer restricted
to the 24-partial demonstration through this front door. Existing partials retain
their geometry/material-derived frequencies. The original 21.6 kHz retention
ceiling still stops inadmissible modes; raising a count is not a full-band or
convergence certificate. It changes the retained endpoint mass and therefore
the loaded board basis: the acoustic solve is prepared from THAT same new bank.

`--substeps 1..16` selects the mechanical rate, 48 kHz times this count. Both
receivers observe every substep on one performance clock; the existing decimator
is prepared with the same factor. Raising the rate does not silently add string
partials or change the frequency band in the acoustic specification. Every
acoustic pole and coupling-strength guard still applies. In particular, the
8-times-band storage pole requires the top fit frequency below
`2700 * substeps` Hz, in addition to the independent acoustic/output guards.
The 10.8 kHz figure in the original feedback description is the DEFAULT 4x case.

`--hammers materials.fsh` uses the complete per-key WoolFelt/Prony input described
in `HAMMERS.md`. `--hammer-footprints faces.fshp` uses `HAMMER_FOOTPRINTS.md` to
select point or finite longitudinal contact faces. Each site keeps independent
felt/relaxation history, sharing its physical hammer, rather than filtering a
point-contact waveform. `--dampers pads.fspd` uses the supplied spatial viscous
pads from `DAMPERS.md`; `--dampers estimated` explicitly selects the existing
approximate spans and drag values. These are not new measured Steinway data.

All supplied cards must cover EVERY admitted scale key, not just the notes in
the score. Missing files, incomplete cards, duplicate options and invalid
budgets refuse before structural/BEM preparation; no source-default fallback
occurs. Constitutive and spatial admission remain with their existing owners.
The output report identifies the chosen controls, actual mechanical rate,
retained string-coordinate count and contact-site count.

`response` and `admittance` also accept `--modes` and `--substeps` after their
output path. Use the same values for a harmonic comparison of a played render.
They reject hammer/damper/score options because those experiments have no
nonlinear contact or key-damper state. The ordinary pressure-only response and
unfitted, fully coupled BEM bridge-force experiment are otherwise unchanged.

More retained partials, contact sites and substeps increase work and memory.
Nothing here certifies real-time performance, spatial convergence, calibrated
materials, full keyboard action, or pressure accuracy outside the sampled band.
