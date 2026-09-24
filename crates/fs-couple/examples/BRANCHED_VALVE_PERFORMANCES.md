# Branched valve ducts, side chambers and compliant walls

A plate-valve performance can supply a complete connected acoustic graph in
place of its `tube` record. It still uses **one moving valve, one plate/material
history and one reciprocal characteristic network**. Branches are not separately
rendered voices. A side chamber or moving wall changes the returning pressure,
which changes the valve's actual contact and motion.

```bash
cargo run --release -p fs-couple --bin music_render -- \
  wind crates/fs-couple/examples/plate-valve-branched.performance \
  /tmp/branched-valve.wav --decimate --block 37
```

This illustrative input keeps the preceding example's plate, spatial lay,
Maxwell material history and pressure phrase. It supplies a 70/80/100 mm main
path with a wider 8 mm-radius outlet, a 43 mm side branch ending in a 10 cm³
Helmholtz chamber, and a locally reacting wall at the main path's second junction.
All dimensions and wall/contact/material values are **authored estimates, not a
measured instrument**. The declared outlet/load band is 1.5 kHz. PCM full scale
is explicitly 0.01 Pa; it never changes physical pressure.

## Graph records

Immediately after `ambient`, use either the original `tube ...` line unchanged,
or the following block:

```text
network NODE_COUNT SECTION_COUNT MAX_WAVE_BYTES
duct_node inlet
duct_node junction
duct_node wall AREA_M2 SURFACE_MASS_KG_M2 STIFFNESS_PA_M RESISTANCE_PA_S_M
duct_node baffled-low-ka MAX_FREQUENCY_HZ
duct_node cavity VOLUME_M3 NECK_RADIUS_M EFFECTIVE_NECK_LENGTH_M RESISTANCE_PA_S_M3
duct_section NODE_A NODE_B LENGTH_M RADIUS_M MAX_LENGTH_ERROR_M
```

Repeat `duct_node` and `duct_section` to match the declared counts. Node indices
are zero-based **in declaration order**, independent of the later plate-mesh
node indices. The graph has exactly one degree-one inlet and is connected. Its
sections are fixed cylindrical waveguides; unequal radii produce the existing
work-conserving impedance-step scattering. Every section has a positive integer
transit within its own declared length-error allowance. No short branch is
silently deleted or assigned an extra junction delay. Parallel paths and loops
are supported by the existing graph owner.

Supported node records are:

```text
duct_node inlet
duct_node junction
duct_node reflection R
duct_node baffled-low-ka MAX_FREQUENCY_HZ
duct_node cavity VOLUME_M3 NECK_RADIUS_M EFFECTIVE_NECK_LENGTH_M RESISTANCE_PA_S_M3
duct_node wall AREA_M2 SURFACE_MASS_KG_M2 STIFFNESS_PA_M RESISTANCE_PA_S_M
duct_node impedance RESISTANCE_PA_S_M3 INERTANCE_PA_S2_M3 COMPLIANCE_M3_PA_OR_none
duct_node series RESISTANCE_PA_S_M3 INERTANCE_PA_S2_M3 COMPLIANCE_M3_PA_OR_none
duct_node shunt RESISTANCE_PA_S_M3 INERTANCE_PA_S2_M3 COMPLIANCE_M3_PA_OR_none
```

`reflection`, `baffled-low-ka`, `cavity` and `impedance` are degree-one terminals.
`junction`, `wall` and `shunt` join at least two sections. `series` connects
exactly two; positive load flow follows section declaration order, not node
numbering. `none` omits compliance; it does not mean a rigid termination. Negative
or nonfinite passive coefficients refuse. No implicit damping or loss law is
selected from a material name.

The chamber uses the existing `HelmholtzLoadSpec`: its physical volume and
supplied effective neck length determine compliance and inertance in the shared
gas. The effective neck length includes the user's intended end corrections.
**Do not also include that neck in a propagating section.** The example's 43 mm
branch is the connecting duct; its separate 20 mm effective neck belongs only to
the lumped cavity. Uniform cavity pressure and a short neck must be justified
for the intended band; this input does not synthesize a three-dimensional cavity.

The wall uses the existing `WallPatch` law, `sigma*x'' + r*x' + K*x = p`, over
its supplied wetted area. Its flow participates in the same node solve. Wall
storage/loss and chamber storage are included once in the network's energy;
they are not output filters. The wall is linear and locally reacting, with no
axial shell coupling or automatic exterior radiation. Its displacement is
inspectable through `WallPatch::observe` on the retained network. The tube's
reference geometry remains fixed. Distributed gas boundary-layer losses are
not implied by these explicit local elements.

## Choose an actual observation

A graph cannot use ambiguous `observation terminal` or `baffled-outlet` aliases:

```text
observation inlet
observation network-node NODE_INDEX
observation network-baffled NODE_INDEX X_M Y_M Z_M RADIAL_RINGS ANGULAR_POINTS MAX_FREQUENCY_HZ
```

An internal node can be a junction, wall or chamber; its pressure is not labelled
as an exterior microphone. A baffled observation must select a declared outlet,
not an enclosed cavity, wall, series node or inlet. Its radius comes from the
adjacent physical section, **not from the inlet**. If that outlet has a compact
radiation load, the receiver's band cannot exceed the admitted load band. All
original compactness, finite-sample fit and receiver-resolution limits apply.

Multiple radiating terminals are allowed. Each receives its own physical load;
the output observes **only the named outlet**, not their exterior sum. These are
independent baffled terminal impedances, with no mutual exterior radiation,
shared bell/room geometry or calibrated microphone claim. Moving only the
receiver leaves the complete valve, wall, cavity and traveling-wave state intact.

## Output, mixing and continuation

The same file works with `ensemble --valve`. Each source keeps its own mechanical
clock and complete physical duration. Decimation must be explicit; physical
propagation delay is not cancelled by filter alignment. Source PCM scales never
become hidden ensemble gains. The sidecar records each graph section's requested
and represented length, radius and one-way sample transit, every lowered load,
and the observed node. **A sum of all section lengths is not an inlet-to-outlet
path length**, especially for branches and loops. The legacy aggregate length
fields are totals for a graph; the individual section records resolve the paths.

The entire runtime remains in `AperturePerformance` during cancellation and
resume. Neither a callback boundary nor pressure release clears any branch,
wall, cavity or material history. Input/clock/topology refusal occurs before WAV
creation; failed physical callbacks retain the established poison contract.
Existing output files are never overwritten. Limits are 64 acoustic nodes,
128 propagating sections, a declared network payload no greater than 64 MiB,
and the original plate, input-byte and schedule budgets. No omitted nodes,
geometry repairs, solver substitutions or new integration scheme are introduced.

Focused native targets are the existing `plate_aperture` and `music_render_wind`
tests. They cover direct physical-network comparisons and work accounting,
side-chamber/wall effects on actual valve motion, inline load pressure drops,
cancellation with retained load history, native/decimated WAV, and mixed-rate
ensembles. These implementation tests do not establish continuum convergence,
measured instrument fidelity or real-time throughput.
