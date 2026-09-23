# Played piano with passive radiation feedback

`piano_exterior render-loaded` connects the supplied body's acoustic pressure
reaction to the nonlinear piano, rather than observing an unaffected piano.
The ordinary `render` command stays one-way for comparisons. `admittance`
remains the unfitted harmonic bridge-force analysis.

```sh
cargo run --release -p fs-couple --example piano_exterior -- \
  render-loaded settled.fss strings.csv acoustic-body.obj acoustic.fspe \
  loaded.wav 2 performance.mid
```

The MIDI argument is optional; without it the existing A4 post-escapement
2 m/s gesture is used. `steinway-d` can replace the CSV to select all 88 source
courses. Reuse the string tensions that produced the supplied board equilibrium.
The board, explicitly labeled closed outward OBJ skin, receiver positions,
medium, physical units and provenance use `EXTERIOR_ACOUSTICS.md`. No missing
cabinet, lid, wood constants or body geometry is manufactured by this command.

Loaded playback adds these explicit limits: retain the **complete** board slice
with at most 32 modes; use an odd 33..257-frequency grid; resolve every prepared
acoustic pole below the mechanical Nyquist guard. The fixed wrapper runs at
48 kHz, four mechanics substeps, and at most 24 partials per string. Its off-band
poles extend to eight times the top fit frequency, so that top frequency must
be **below 10.8 kHz**. The complete per-string/duplex and geometry admission
rules still apply. Larger bases are refused, not truncated to fit the budget.
This is bounded offline rendering, not a real-time performance claim.

## The load is from the same acoustic solve as the receivers

At each supplied frequency one existing Helmholtz BEM batch solves all retained
unit board-velocity fields. The force projection `Z = G^T A P` uses physical
panel areas once, with signed full-vector skin motion in the actual string-
mass-loaded board basis. Rigid scatterers participate in this same solve.
The same fields supply the receiver transfers; pressure per acceleration is
pressure per velocity times `i/omega` under `exp(-i omega t)`.

The acoustic impedance is approximated as

```
Z(s) = sum_j l_j^T l_j s / (s^2 + 2 zeta_j Omega_j s + Omega_j^2).
```

Signed vectors `l_j` retain cross-mode coupling. Every residue is a real
positive-semidefinite matrix, so the complete multiport model is passive by
construction. Checking only the diagonal entries, or fitting every entry as
an unrelated SISO filter, would not establish this property.

This is a **restricted passive matrix fit**, not a general MIMO vector fitter.
There are 31 logarithmic damped-resonator dictionary entries (damping ratio
0.2) and one off-band lossless storage entry. The latter allows a compact
inertive load without inventing a minimum acoustic resistance. Full residue
matrices are fitted with bounded block-coordinate least squares, projecting
fit iterates through the existing symmetric eigenanalysis owner. Raw BEM
matrices are never overwritten or symmetrized to make validation pass. Each
positive residue direction becomes one unit-mass acoustic coordinate; the
worst admitted image has 1024 coordinates, not a silently rank-truncated fit.

Only even grid samples determine coefficients, relative resistive weighting
and convergence stopping. Odd samples remain held out. Both sets must meet
fixed full-matrix bounds: complex peak/RMS error at most 5%/2%, and Hermitian
power-part peak/RMS error at most 10%/5%. Peak normalization is by each set's
peak Frobenius reference norm; RMS is the relative aggregate norm. The power
part is checked separately so a large reactive load cannot hide a bad estimate
of weak radiation loss. Reported errors are the worse of training and held-out
results. These are sampled estimates, not interval or interpolation certificates.

A failed passive fit, receiver fit, modal solve, resolution check or runtime
admission **rejects loaded playback**. There is no automatic one-way fallback,
relaxed error threshold, discarded mode or fabricated damping to obtain a WAV.

## How the pressure reacts during a note

Acoustic coordinates carry storage `sum(p_air^2 + Omega^2 q_air^2)/2`.
The power connection drives their momenta with `l*v_board` and applies the
opposing `-l^T*p_air` force to the board. The cross-power cancels exactly in
the continuous equations; wood, string, hammer, felt and damper losses remain
owned by their original models.

At every mechanics substep a symmetric composition of exact pairwise velocity
rotations and existing exact acoustic free-oscillator flows surrounds the
original nonlinear piano step. This is **second-order operator splitting**,
not an exact full-system propagator. Coupling rotations conserve combined
kinetic energy; damped free flows remove their reported energy. A separate
rate-resolution bound rejects excessive coupling strength. Passivity alone
is not a time-step accuracy certificate: convergence with mechanical rate is
still required, especially during hammer attacks and near acoustic poles.

Acoustic history is checkpointed with strings, board, hammers, felt and the
accounting ledger. A rejected output sample restores all of them. Installation
is cold and once-before-excitation; it cannot reset acoustic memory mid-note.
The unselected path is unchanged, and a zero-coupling model preserves the
original mechanical trajectory.

Both receivers observe every **reacted** mechanical substep on one score clock.
The existing causal anti-alias decimator, propagation delays, receiver filters
and physical-Pa PCM writer are reused. There is no per-channel mechanical
advance, independent normalization, limiter or output equalizer substituting
for pressure reaction. The report includes acoustic stored energy and acoustic
dissipation separately, each counted once in the combined work balance.

## Accuracy boundary

The fit represents the declared BEM load only over the sampled band. The
lossless off-band storage pole and all other extrapolated behavior are numerical
realization choices, not identified physical air eigenmodes. A broadband attack
can excite this unvalidated region. Fitted acoustic loss is the passive model's
energy removal, not an exact broadband far-field energy integral. The pressure
observer has its own fit bounds; its out-of-band energy is not independently
certified against the load's dissipation.

The soundboard remains linearized about the supplied equilibrium, the dynamic
strings use the existing transverse model, and cabinet/lid components remain
rigid. This does not add room acoustics, flexible rim mechanics, axial string
motion, a full keyboard action, new external Steinway mesh bytes or factory
material measurements. No measured-instrument or full-band realism claim follows
from completing this time-domain coupling.
