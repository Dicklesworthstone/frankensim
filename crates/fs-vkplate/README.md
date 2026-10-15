# fs-vkplate: a drum kit from one solver

Every vibrating surface of a standard kit is one object here: a thin plate on a square
finite-difference grid. A drum head is a clamped plate with in-plane tension and the film's own
small bending stiffness; a cymbal is a free-edged plate with bending, a rest shape (a shallow
shell) and full von Karman coupling. Kick, snare, toms, hi-hat, crash and ride are parameter sets of the same
equations, couplings and output stage (`src/kit.rs`).

Why it exists: FrankenSim's modal percussion (`crates/fs-couple/examples/percussion`) keeps at
most 63 head modes, admits audio only for head bands inside 40-1640 Hz, and caps a cymbal shell
at 32 modes (an 18-inch crash stops near 400 Hz). A handful of lightly damped modes is a bell,
and that is what those renders sounded like. A grid solver needs no mode list: an 18-inch crash
on a 7 mm grid carries 3356 nodes and resolves the whole audio band.

## Physics

- **Membranes**: `rho t w_tt = T Lap(w) - D Lap^2(w) + ...`, rim clamped, with the Mylar film's
  bending stiffness `D = E t^3 / 12(1 - nu^2)`, which takes over from tension for waves with
  `D k^2 > T` (above roughly 0.7 kHz on the slack 2-ply kick and floor-tom batters, 5 kHz on the
  snare and rack tom). Berger tension modulation (the mean stretch raises the tension) gives the pitch drop of a
  hard-struck tom.
- **Plates and shells**: bending energy summed over 3x3 curvature stencils with the full
  `(1 - nu)` twist term; free edges are the natural boundary of that energy. Von Karman coupling
  through the Airy stress function with a rest shape `z0(r)`:
  `Lap^2 F = -E t [l(w + z0) - l(z0)]`, `l(u) = u_xx u_yy - u_xy^2`, `F = dF/dn = 0` at the rim.
  The nonlinear force is the exact gradient of the discrete potential.
- **Loss**: `s0 + s1 beta^2` per mode, fixed by two (frequency, T60) points.
- **Contacts** (stick, felt beater, snare wires, hi-hat rims): one-sided power laws with
  Hunt-Crossley loss, `f = K d^alpha (1 + mu d_dot)`, each solved implicitly as a monotone scalar
  equation, so a hickory tip on bronze is stable at audio rates.
- **Air**: a lumped adiabatic cavity couples a drum's two heads; each head's volume motion
  carries the radiation mass of a baffled piston (`8 rho a^3 / 3`, implicit); a monopole
  radiation resistance acts on the drum's NET volume flow (two heads moving the same way in
  space radiate as a weak dipole, which is why a tom's fundamental rings). An optional port
  models a vented front head.
- **Sound**: the Rayleigh integral of surface acceleration at each microphone, per-node delays.
  A cymbal is not in a baffle: its microphones apply the unbaffled-disc limit (a dipole,
  `(4/3pi) ka` times the baffled field, for `ka << 1`), joined to the baffled field for `ka >> 1`
  by a first-order high-pass with `tau = 4a/(3 pi c)`.

## Numerics

Leapfrog in time. The nonlinear potential (Berger or von Karman) enters through a scalar
auxiliary variable `psi = sqrt(2V)` (Bilbao, Ducceschi and Zama): the update stays explicit,
needs one rank-one correction (Sherman-Morrison; with the air mass a 2x2 Woodbury solve), and
conserves a discrete energy, so the nonlinearity cannot destabilise the step. The plain explicit
von Karman update blew up a crash struck at 8 m/s within 7000 steps; this one stays bounded.
`psi` is advanced along the gradient only, so it drifts from `sqrt(2V)`; a shallow shell carries its
curvature stiffness in `V`, so that drift scales the shell's stiffness (on a live grid it wandered
between 0.5 and 1.5 within 2 s and the 8-16 kHz bands grew back after a second). Each step pulls `psi`
back to `sqrt(2V)` at the half step: lowering it only removes energy, and raising it spends at most
the energy the loss terms removed in the step before, so the discrete energy still never grows.
The linear part's CFL limit is checked by power iteration and the internal rate is the smallest
multiple of 48 kHz that uses at most 60 % of it. Output is decimated with a windowed sinc.

The Airy operator is factored once (banded Cholesky, row-major ordering); both triangular sweeps
stream the one stored factor, which dominates a cymbal's cost on a fine grid. Contacts warm-start
their Newton solve from the previous step's force. Contacts on disjoint nodes (the snare wires, the
hi-hat rims) are solved four at a time, interleaved, each with the same operations it would take
alone, so the result is the one-at-a-time result bit for bit. A snare wire's solve stops when a
Newton step lands exactly on its own iterate (an exact fixed point) instead of bisecting away from
it: about 3 iterations per solve instead of 10, and the applied force is within a few ulps of the
root rather than 1e-13.

## Hi-hat

Two shells on one rod. The bottom shell sits on the stand's felt; the top shell hangs from the
clutch (felt washers gripping the annulus round its hole) on the pull rod, which a return spring
holds up and the pedal pulls down. With 80 N on the pedal about 62 N reaches the rims (the spring
keeps the rest). An earlier version also tied the top shell's centre to the ground through felt,
which took ~75 N of the pedal force, so the "closed" hat was barely pressed.

No rim is flat: the rim contact carries a two-lobe waviness of up to 0.15 mm (offline grid), and
with no pedal force the shells rest at most 0.14 mm apart there, so a stroke on the open hat sets
the rims chattering against each other (its sizzle above 4 kHz). The rim contact is lightly damped
(Hunt-Crossley `mu` 0.08 s/m); tangential friction between the rims is still not modelled.

## Validation (`cargo test --release`)

- clamped membrane: fundamental within 1.5 % of `2.405 c / (2 pi R)`;
- free disc: lowest elastic mode within 4 % of `lambda^2 = 5.253` (nu = 0.33);
- crash struck at 8 m/s: bounded, finite output;
- a hi-hat's second stroke retires the stick's earlier contact (each stroke adds one);
- a shallow shell's auxiliary variable keeps tracking `sqrt(2V)` (live-grid ride, 2 s);
- the stencil rows of the sparse products sum bit for bit as their CSR rows (five pieces, signed zeros);
- the four-lane contact solves match solving one contact at a time, force, penetration and every
  update bit for bit, with node-sharing and other-law contacts in the list;
- every piece builds on the live grid, sounds, and the kick and snare stop being stepped once
  silent; a clone, taken settled or mid-ring, renders bit for bit what the original would (`live`).

## Real time (`live`)

`live::Live` is one piece on a grid `scale` times coarser (5 runs several times faster than real
time in WebAssembly), struck anywhere on its playing surface, stepped in 48 kHz blocks, and not
stepped at all once its motion is 70 dB below its peak. On grids too coarse to put four nodes in
the 6.5-20 mm annulus round a cymbal's hole, the felt washers and the hi-hat clutch spread their
total stiffness and damping over the first ring of nodes instead. The live hi-hats keep the rim
contact law from before the retuning (Hunt-Crossley `mu` 3.0 s/m, set in `live.rs`) and a 5 mm
stick tip; the offline hat's lighter rim damping (0.08 s/m) and 4 mm tip were not carried over to
real-time play. A clone of a `Live` steps bit for bit as the original would, so the browser kit
settles each grid once: it plays a clone kept from before the timing stroke.
`apps/instrument-lab/engine-kit` compiles it for the browser.

## Graded against a real kit

The kit's sizes follow the Salamander Drumkit (Alexander Holm, CC BY-SA 3.0): 18x18 kick,
14x6 snare, 12x7 rack tom, 14x14 floor tom, 14-inch hi-hat, Paiste 18-inch crash and 20-inch ride,
recorded with overhead microphones. Head tensions are the drummer's tuning, set so each
fundamental lands on the recording's; loss points come from the recordings' band decay times.
`apps/instrument-lab/tools/kitgrade.py` compares every piece band by band.

## Estimates (not measurements)

Film and bronze properties are handbook values; head tensions are tuned; the stick's effective
tip mass (6-7 g), its Hertz constant (hickory on the target), the felt beater law (a 32 g crown in
point contact), strike points, snare-wire mass, tension and preload, felt washer stiffness and
loss, and the hi-hat rim contact law, rest gap and waviness are declared estimates in
`src/kit.rs`. The lumped cavity has no cavity modes; there is no shell, hoop or rim
vibration, no tangential friction between hi-hat plates, and each piece is simulated alone (a kit
mix adds independent renders: the snare does not buzz when the kick is struck).

## Usage

```sh
cargo run --release --example kit -- crash crash.wav          # PIECE OUT.wav [SPEED_M_S] [SECONDS]
cargo run --release --example perform -- kick hits.csv kick.wav 17   # time_s,speed_m_s per stroke
cargo run --release --example livebench -- 5                   # real-time factor per piece at grid x5
cargo test --release
```
Pieces: `kick snare snare-off tom-12 tom-14 hihat-closed hihat-open crash ride`.
