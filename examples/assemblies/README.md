# Assembly Gallery (`examples/assemblies/`)

Tracked `.fsim` assembly files that run directly through `frankensim run` with
measured `fs-matdb` card packs and independent closed-form analytical oracles.
No physics behavior is tuned to the expected answer; every parameter is declared
in the `.fsim` file or queried from a measured card with provenance.

## Assemblies and Expected Emergence

| Assembly | Measured cards | Modes | Expected emergence & independent closed-form oracle |
| --- | --- | --- | --- |
| `forced_oscillator.fsim` | `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (`rho = 7800 kg/m^3`, box `0.1 x 0.05 x 0.05 m` -> `m = 1.95 kg`) | `transient`, `modal`, `harmonic`, `stationary` | 1-DOF damped oscillator (`k = 200 N/m`, `c = 4 N*s/m`, constant `F = 10 N`). Modal: `omega_n = sqrt(k/m) = 10.12739367 rad/s`, damped `f_d = omega_n*sqrt(1 - zeta^2)/(2*pi) = 1.60353769 Hz`, `zeta = c/(2*sqrt(k*m)) = 0.10127394` (`<= 1e-9`). Harmonic (`1, 2, 5 Hz`, spring root-velocity input -> spring extension): `|H(i*w)| = |(c + i*m*w)/(k - m*w^2 + i*c*w)|` (`<= 1e-9`). Stationary: `|x_eq| = F/k = 0.05 m`, `v_eq = 0` (`<= 1e-9`). Transient (`50` steps, `dt = 1 ms`): underdamped step response `x(t), v(t)` (`<= 1e-5 m`, `<= 1e-4 m/s`) and energy residual `<= 1e-12 J`. |
| `two_body_conduction.fsim` | `packs/aluminum-pure-thermal-nbs-c447.fsmcdpk` (`rho = 2700 kg/m^3`, `c_p = 945.584 J/(kg*K)`) + `../system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk` | `transient` | Two pure-Al blocks of distinct size (`0.06 x 0.05 x 0.05 m` -> `C_h = 382.96152 J/K`; `0.04 x 0.05 x 0.05 m` -> `C_c = 255.30768 J/K`) sharing a `0.0025 m^2` contact face at `P = 313 kPa` (`R'' = 3.1606e-4 m^2*K/W`, `G = A/R'' = 7.909890527 W/K`), `500 W` heater on `hot` for `5.0 s` (`dt = 0.05 s`). Both blocks remain inside the Atkins validity window (`289.35..305.85 K`). Exact two-lump exponential `T_hot(t), T_cold(t)` (`<= 1e-4 K`) and total enthalpy conservation `|C_h*(T_h - T_0) + C_c*(T_c - T_0) - Q*t| <= 1e-8 J`. |
| `thermo_mech_slider.fsim` | `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (`m = 1.95 kg`, `c_p = 486 J/(kg*K)` -> `C = 947.7 J/K`) | `transient` | Coupled thermo-mechanical damped slider (`k = 400 N/m`, `c = 12 N*s/m`, `F(t) = 5 + 15*sin(2*pi*3*t + 0.25) N`, `horizon = 0.05 s`, `dt = 1 ms`) where viscous damping heats the body and a thermocouple lags raw temperature (`tau = 5 ms`). Closed-form forced damped oscillator `x(t), v(t)` (`<= 1e-5 m`, `<= 1e-4 m/s`), exact step-mean midpoint-quadrature energy identity `|E_mech + E_th - W_ext| <= 1e-9 J` and `|E_th - Q_damp| <= 1e-9 J`, and discrete lag filter on `slider_temp` (`<= 1e-12 K`). |
| `shear_film_pair.fsim` | `packs/aluminum-pure-nbs-c447-sensible-heat.fsmcdpk` (`rho = 2700 kg/m^3`, `c_p = 945.584 J/(kg*K)`, `sensible-heat` law valid `288.15..298.15 K`) + `packs/dowfrost-pg50-40-60c.fsmcdpk` (`mu = 3.10 mPa*s` at `313.15 K`) | `transient` | Two free pure-Al discs (`r = 0.05 m`, `h = 0.01 m` -> `m = 0.21205750 kg`, `I_z = m*r^2/2 = 2.6507188e-4 kg*m^2`) stacked across a PG50 shear film (`gap = 0.1 mm`, full face `A = pi*r^2`), counter-rotating at `10 rad/s` each about their shared axis as in a parallel-plate rheometer. The film lowers to `c_t = mu*J/gap` with `J = A^2/(2*pi)` (`3.0434179e-4 N*m*s`). `0.7` and `0.2` of its heat go to two insulated pure-Al cells (`20.4246144 J/K`, declared `293.15 K`), and the rest to `unallocated`. Free-pair exponential (the spin axis joins the centres, so the film has no lever): dissipated `= (1/2)*(I_z/2)*(20 rad/s)^2*(1 - exp(-4*c_t*t/I_z))` (`<= 1e-12` relative, `0.026238761 J` at `t = 1 s`). Window ledger, each `<= 1e-12` of dissipated: cells' stored change `0.9x`, `unallocated` `0.1x`, the network's stored change `-1x`, energy residual and numerical loss; environment exactly `0`. Probes: each disc's `angular-velocity-z` `= +-10 rad/s * exp(-2*c_t*t/I_z)` (`<= 1e-12` relative; the spin decays at half the energy's rate, so about `10 %` of it is left at `t = 1 s`), and each cell's `C*(T - 293.15 K)` is its share of the closed-form heat within `2*ulp(T)*C`, two roundings at `T`'s size, from each sample's own `T`. |
| `srm1453_slab.fsim` | `packs/eps-nist-srm1453-2013.fsmcdpk` (NIST SRM 1453 / SP 260-175 `Surface2D` `lambda(Tm, rho)`, `1.5 %` relative half-width at `95 %`, `k = 2`) | `stationary` | 1-D steady conduction across a single SRM 1453 board thickness `0.01343 x 0.25 x 0.25 m` (SP 260-175 Table 8 certification thickness `13.43 mm`, within Table 2 production range `[12.86, 13.88] mm`; `:fidelity (one-d-field 1)`, `T_hot = 307.15 K`, `T_cold = 287.15 K`, `Tm = 297.15 K`, `rho = 41 kg/m^3`, `dT = 20 K`). Closed-form Fourier heat rate `Q = k * (0.25*0.25/0.01343) * 20 = 3.122100707371556 W` (`<= 1e-12` relative), `stationary.residual_norm = 0`. Eight-term stationary QoI budget against the declared design limit `3.35 W` (nominal margin `0.2279 W`, total conservative half-width `0.1394 W` from SRM `1.5 %` conductivity `0.0468 W` plus declared apparatus tolerances, not measured: box-axis geometry corners `0.0720 W`, face-temperature BC corners `0.0156 W`, and measurement `0.0050 W`; remaining margin `0.0885 W` -> determinate `pass`). The SRM 1453 certificate's expanded uncertainty (SP 260-175 Eq. (17), sec. 7.2, `k = 2`; `lambda_exp` budget Eq. (8), incl. the imbalance term `dQ`, A1.3.2) already includes the 1-D guarded-hot-plate apparatus model error for metered-section heat rate under certified conditions (`:covered-by-source`, gated on all certificate conditions including thickness `[12.86, 13.88] mm`). |
| `two_mass_chain.fsim` | `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (`rho = 7800 kg/m^3`: boxes `0.1 x 0.05 x 0.05 m` -> `m1 = 1.95 kg`, `0.06 x 0.05 x 0.05 m` -> `m2 = 1.17 kg`) | `transient`, `modal`, `stationary` | **Added as a file only: nothing in `crates/` changed for it.** Two steel blocks in a chain, the inner joined to the root and the outer to the inner by prismatic joints (`k1 = 200 N/m`, `c1 = 4 N*s/m`; `k2 = 150 N/m`, `c2 = 3 N*s/m`), with a constant `10 N` on the outer block from rest. The damping is stiffness-proportional (`c = 0.02 k`), so the undamped modes decouple: `omega^2` are the roots of `m1*m2*w^4 - (m1*k2 + m2*(k1 + k2))*w^2 + k1*k2 = 0` (`omega_1 = 7.16114874 rad/s`, `omega_2 = 16.01281538 rad/s`) and `zeta_i = 0.02*omega_i/2`. Modal: each mode's `natural_rad_s`, `frequency_hz` (damped) and `damping_ratio`, reported once per conjugate eigenvalue (`<= 1e-9`). Stationary: the joint extensions `-F/k1 = -0.05 m` and `-F/k2 = -0.0667 m`, both velocities zero (`<= 1e-9`). Transient (`1000` steps, `dt = 1 ms`): both positions against the exact modal superposition (`<= 5e-6 m`), and again at `dt = 0.5 ms`: the error falls by `4` (`3.9..4.1`), the integrator's second order; energy residual `<= 1e-12 J`. |

### Material Card Disclosure (`two_body_conduction.fsim`)

The bulk heat capacities use elemental aluminum from NBS Circular C447 Table 67
(`aluminum-pure-nbs-c447`), whereas the vacuum contact-resistance curve is
Atkins (1964) Table II runs 23-28 measured on lathe-cut 6061-T6 aluminum alloy
faces (`6061-T6 aluminum alloy`). Both parts are aluminum so only the alloy
differs between bulk and surface cards; this mismatch is disclosed rather than
masked under a single material label.

**One-way film (`shear_film_pair.fsim`):** `mu` is read once at the declared `313.15 K`
and held; it does not follow the film's temperature, and the receipt marks that key `held`.

## Reproducing the Tracked Card Packs

The material packs under `examples/assemblies/packs/` are built deterministically
from `data/matdb/seed-v1/` via `xtask matdb-pack`:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/steel-music-wire-string-heptonstall-ligo-p1400229/manifest.tsv \
  --material-state "steel-music-wire-ligo-p1400229" "solid" "drawn untreated music wire" 0 \
  --out examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aluminum-pure-thermal-nbs-c447/manifest.tsv \
  --material-state "aluminum-pure-nbs-c447" "solid" "elemental aluminum, NBS C447 Table 67" 0 \
  --out examples/assemblies/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/aluminum-pure-thermal-nbs-c447/manifest.tsv \
  --material-state "aluminum-pure-nbs-c447" "solid" "elemental aluminum, NBS C447 Table 67, sensible-heat law" 0 \
  --models data/matdb/seed-v1/aluminum-pure-sensible-heat-nbs-c447/manifest.tsv \
  --out examples/assemblies/packs/aluminum-pure-nbs-c447-sensible-heat.fsmcdpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/dowfrost-pg50-40-60c/manifest.tsv \
  --material-state "dowfrost-pg50" "liquid" "inhibited DOWFROST aqueous solution, 50 volume percent propylene glycol, Dow 2007 guide" 0 \
  --out examples/assemblies/packs/dowfrost-pg50-40-60c.fsmcdpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/eps-nist-srm1453-2013/manifest.tsv \
  --material-state "eps-nist-srm1453-2013" "solid" "expanded polystyrene thermal insulation board, NIST SRM 1453 (SP 260-175)" 0 \
  --out examples/assemblies/packs/eps-nist-srm1453-2013.fsmcdpk
```

The third pack has the second one's claims plus the `sensible-heat` model, so
its thermal-only solids lower to the owners' thermal timeline; the two-body file
keeps the second pack and its lumped lowering.

The thermal contact interface pack is reused in place from
`examples/system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk`
(built from `data/matdb/seed-v1/atkins-1964-al6061t6-lathe-0p30um-vacuum-interface/`).

## Running the End-to-End Lane

```bash
scripts/e2e/assemblies.sh --artifact-dir /tmp/frankensim-assemblies-e2e \
  --prebuilt-binary target/debug/frankensim
```

Omit `--prebuilt-binary` on a clean git tree to build and run directly from
committed sources.

## Scope & No-Claim Boundary

- **Requirement verdicts & eight-term QoI budgets:** Every gallery `.fsim`
  run—including `shear_film_pair.fsim` (`cell_temp_limit <= 294 K`, declared
  limit (illustrative)) and all four `forced_oscillator.fsim` modes (`pos_limit`,
  `mode1_freq >= 1.5 Hz`, `harm_peak <= 0.25 s`, and `eq_extension >= -0.06 m` on
  `stationary.state.joint.spring_damper.spring.state`, which follows the lowered
  node name; all declared limits (illustrative))—ends with a populated eight-term
  verdict table. Nominal values come from the declared step/analysis run,
  discretization comes from `2 * |stat(dt) - stat(dt/2)|` (bounding the
  closed-form errors: `2.01e-7 >= 1.34e-7 m` on `forced_oscillator`,
  `1.68e-6 >= 1.12e-6 K` on `two_body_conduction`; halved linearisation/probe
  step on modal, harmonic, and `dim > 0` stationary), and roundoff stays
  `NO-DATA` where no verified conversion exists (modal, harmonic, `dim > 0`
  stationary). Built-in state/transfer probe units come from each component's
  declared state/port dimensions; requirements on untyped `"SI"` probes refuse by
  name. Gallery cards other than SRM 1453 carry `UncertaintyModel::Unstated`, so
  those verdicts remain `indeterminate` (`Estimate`-only). `srm1453_slab.fsim` is the
  first gallery assembly with a determinate `pass`: its `3.35 W` threshold is a
  declared design limit (not a measured or certified value), its geometry,
  boundary-condition, and measurement half-widths are declared apparatus
  tolerances (not measured), and its `:covered-by-source` model-form term
  reflects that the SRM 1453 certificate's `k = 2` expanded uncertainty (SP
  260-175 Eq. (17), sec. 7.2; `lambda_exp` budget Eq. (8), incl. imbalance `dQ`,
  A1.3.2) already includes the 1-D ASTM C177 guarded-hot-plate apparatus
  for heat rate through a single board (`[12.86, 13.88] mm`) under certificate
  conditions.
- **Heatsink-as-assembly (gallery item 6):** `examples/heatsink-fan/heatsink-fan.fsim`
  is not yet expressible as a `.fsim` assembly because it requires STL geometry
  ingestion, non-lumped CHT/conduction fidelity, and 3-D field QoIs beyond the
  current lumped assembly lowering path.
- **Authority:** These runs provide **Estimate-only** numerical, conservation,
  determinism, and provenance verification against closed-form analytical
  oracles; they make no experimental or physical-validation claim.
