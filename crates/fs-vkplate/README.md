# fs-vkplate: a drum kit from one solver

Every vibrating surface of a standard kit is one object here: a thin plate on a square
finite-difference grid. A drum head is a clamped plate with in-plane tension and negligible
bending; a cymbal is a free-edged plate with bending, a rest shape (a shallow shell) and full
von Karman coupling. Kick, snare, toms, hi-hat, crash and ride are parameter sets of the same
equations, couplings and output stage (`src/kit.rs`).

Why it exists: FrankenSim's modal percussion (`crates/fs-couple/examples/percussion`) keeps at
most 63 head modes, admits audio only for head bands inside 40-1640 Hz, and caps a cymbal shell
at 32 modes (an 18-inch crash stops near 400 Hz). A handful of lightly damped modes is a bell,
and that is what those renders sounded like. A grid solver needs no mode list: an 18-inch crash
on a 7 mm grid carries 3356 nodes and resolves the whole audio band.

## Physics

- **Membranes**: `rho t w_tt = T Lap(w) + ...`, rim clamped. Berger tension modulation (the mean
  stretch raises the tension) gives the pitch drop of a hard-struck tom.
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
The linear part's CFL limit is checked by power iteration and the internal rate is the smallest
multiple of 48 kHz that uses at most 60 % of it. Output is decimated with a windowed sinc.

The Airy operator is factored once (banded Cholesky, row-major ordering).

## Validation (`cargo test --release`)

- clamped membrane: fundamental within 1.5 % of `2.405 c / (2 pi R)`;
- free disc: lowest elastic mode within 4 % of `lambda^2 = 5.253` (nu = 0.33);
- crash struck at 8 m/s: bounded, finite output.

## Graded against a real kit

The kit's sizes follow the Salamander Drumkit (Alexander Holm, CC BY-SA 3.0): 18x18 kick,
14x6 snare, 12x7 rack tom, 14x14 floor tom, 14-inch hi-hat, Paiste 18-inch crash and 20-inch ride,
recorded with overhead microphones. Head tensions are the drummer's tuning, set so each
fundamental lands on the recording's; loss points come from the recordings' band decay times.
`apps/instrument-lab/tools/kitgrade.py` compares every piece band by band.

## Estimates (not measurements)

Film and bronze properties are handbook values; head tensions are tuned; the stick's effective
mass (12-20 g), its Hertz constant (hickory on the target), the felt beater law, snare-wire mass,
tension and preload, felt washer stiffness and loss, and the hi-hat rim contact law are declared
estimates in `src/kit.rs`. The lumped cavity has no cavity modes; there is no shell, hoop or rim
vibration, no tangential friction between hi-hat plates, and each piece is simulated alone (a kit
mix adds independent renders: the snare does not buzz when the kick is struck).

## Usage

```sh
cargo run --release --example kit -- crash crash.wav          # PIECE OUT.wav [SPEED_M_S] [SECONDS]
cargo test --release
```
Pieces: `kick snare snare-off tom-12 tom-14 hihat-closed hihat-open crash ride`.
