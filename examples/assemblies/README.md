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
| `gas_piston.fsim` | `packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk` (measures dilute-gas `N2` molar mass & NASA-9 ideal-gas heat capacity, McBride et al. 2002, plus dilute-gas viscosity/thermal conductivity, Lemmon & Jacobsen 2004) + `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (measures solid steel music-wire mass density `7800 kg/m^3`, Heptonstall LIGO P1400229; `0.1 m` cube -> `m = 7.8 kg`) | `modal`, `transient` | Insulated (adiabatic) closed `N2` `GasCavity` (`r = 0.15 m`, `h = 0.20 m`, `T0 = 300 K`, `p0 = 90 kPa`) coupled to a translating steel piston through a card-less `piston-face` (`common_flow_ratios`) with explicit exterior ambient pressure force `-p0*pi*r^2 = -6361.725123519331 N` (no hidden ambient force). Modal: ideal-gas adiabatic Rüchardt frequency `f_ad = sqrt(gamma(T0)*p0*A^2/(m*V0))/(2*pi)` (`<= 1e-6` relative, `damping_ratio <= 1e-6`) plus eight-term verdict table on `ruchardt_freq >= 11.5 Hz` (declared illustrative limit). Transient: rest equilibrium holds (`|v| <= 1e-11 m/s`, `|T - 300 K| <= 1e-11 K`, energy residual `<= 1e-12 J`, eight-term verdict table on `gas_temp_max <= 310 K`); a `+500 N` compression step matches an independent adiabatic RK4 ODE (`<= 5e-6 m`, `<= 1e-4 m/s`, `<= 2e-3 K`, energy residual `<= 1e-9 J`). The finite-`G` ambient-bath thermal-wall cubic check stays in `crates/fs-system/tests/gas_cavity_slice.rs` and `crates/fs-cli/tests/assembly_run.rs` rather than the gallery because no measured gas-to-wall heat-transfer card is in `seed-v1` (Atkins 1964 measures solid-solid 6061-T6 contact resistance in vacuum). |
| `block_on_belt.fsim` | `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (`0.05 x 0.05 x 0.02564102564102564 m` -> `m = 0.5 kg`) + `packs/katsivardi-2025-pla-metal-dry-air-interface.fsintpk` (Katsivardi et al. 2025, *Materials* 18(12):2730, CC BY 4.0, ASTM D1894 sliding-sled Table 3: pure PLA circular pattern on dry metallic plate at `N = 5 N`, `v = 0.0025 m/s`, `mu_s = 0.371`, `mu_k = 0.195`) | `transient` | Non-smooth Coulomb stick-slip oscillator (device design spring `k = 1000 N/m`, `N = 5 N` as printed for the 500 g weight, constant belt speed `v = 0.0025 m/s`, `dt = 1 ms`, `horizon = 3.6 s`; neither side has thermal physics so all friction heat routes to ambient under `default: no thermal side; all heat to ambient`). Closed-form stick-slip period `T_exact = 2*d/v + 2*(pi - atan(omega*d/v))/omega = 0.777085244 s` (`d = (mu_s - mu_k)*N/k = 0.88 mm`, `omega = sqrt(k/m) = 44.72135955 rad/s`): every measured cycle period is within `< dt` (step-boundary breakaway/re-catch by `Slider` design); average dissipated work per cycle matches `W_cyc = mu_k*N*v*T_exact = 1.894145 mJ` within `dt/T_exact` plus closure; peak slip speed matches `v + sqrt(v^2 + (d*omega)^2) = 41.9 mm/s` (`16.8 x v`; peak block speed `39.4 mm/s`); ambient heat equals dissipated work (`<= 1e-12` rel), energy residual `<= 1e-12 J`, and `block_pos_cap <= 0.005 m` (declared illustrative limit) emits an eight-term verdict table (`indeterminate`, Estimate-only). |
| `rubber_mount_creep_288k.fsim`, `rubber_mount_creep_298k.fsim` | `packs/rouleau-2015-self-adhesive-rubber-shear.fsintpk` (Rouleau, Pirk, Pluymers and Desmet 2015, J. Aerosp. Technol. Manag. 7(2) 200-208: a self-adhesive synthetic rubber in DMA double shear, Table 1's generalized Maxwell fit at `Tr = 20 C`, `G0 = 0.15 MPa` and 12 branches `G_k`, `tau_k`, with WLF `C1 = 9.23`, `C2 = 141.2 K`; each stiffness is `k = G * 2A/h = G * (1/11) m` for the source's two `10 x 10 mm` samples at their nominal `2.2 mm`; the specimen is shear-only, so the card states axis `x` and leaves the other five axes unstated, `0` with a receipt each) + `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (boxes `0.1 x 0.05 x 0.05 m` -> `1.95 kg` and `0.05 m` cube -> `0.975 kg`) + `packs/aluminum-pure-nbs-c447-sensible-heat.fsmcdpk` (`0.05 m` chamber cube, `C = 319.1346 J/K`) | `transient` | Creep of the rubber mount between two free steel boxes: `drive_x` holds `0.1 N` from relaxed rest (every branch's and the equilibrium spring's displacement declared `0 m`), and all of the mount's heat goes to the aluminium chamber cell, whose temperature shifts the branch times by `a_T`: `2.1819` at `288.15 K` and `0.49373` at `298.0 K`. The two temperatures are bounded by the aluminium card's `sensible-heat` band `288.15..298.15 K`, not by the source's TTSP range `-30..30 C` (`298.0 K` leaves room for the cell's warming); a colder pair needs a sensible-heat card with a wider band. Two `0.25 s` steps per run (`0.5 s` of creep, against `a_T*tau_1` of `0.71 s` and `0.16 s`), exact at their ends: the owner advances the velocities and the 13 stated `x` memories as one exponential a step, and each unstated axis's memory (`0` here) by its own branch's exact decay. The two runs take `7.9 s` and `6.4 s` in a debug build on one core and `2.1 s` and `1.0 s` in release; before bead `frankensim-ks7e6` took the unstated axes out of the exponential they took `185 s` and `195 s`, and `13.4 s` and `12.4 s`. Oracle: the moving box's `velocity-x` at each step end against `exp(A dt)` of the relative system `[v, e_eq, e_1..e_12, 1]` in 50-digit decimal (`<= 1e-9` of its largest value), and each of three planted mutations (equilibrium spring dropped, `a_T = 1`, `tau_1` and `tau_2` swapped) misses by `>= 100x` that. Window ledger: drive work `F*u`, dissipated `W - E` and the network's net port energy `E` (memory plus relative kinetic energy), each `<= 1e-8` relative; `unallocated` exactly `0`; the cell's `C*(T - T0)` is the dissipated heat within `2*ulp(T)*C`; energy residual `<= 1e-18 J`. The cell's warming moves `a_T` by `<= 1e-9` of itself, reported per run. |
| `servo_mass.fsim` | `packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk` (`0.1 m` cube -> `m = 7.8 kg`) | `transient` | A free steel cube pushed by a constant `-12 N` and held at `0.05 m` by a declared PID position controller (`kp = 72 N/m`, `ki = 36 N/(m*s)`, `kd = 38.4 N*s/m`) sampled at `16 Hz`, its command held over the four `1/64 s` steps of each sample. Oracle: the same discrete loop written out (at sample `k`, `e = r - x`, `S += e`, `u = kp*e + ki*T*S + kd*(e - e_prev)/T`, then four midpoint steps of the free body under `u - 12 N`): position and velocity at all `769` samples (`<= 1e-12`). After `12 s` the integral term holds the push and the body is within `1e-4 m` of the setpoint (without that term it would settle at `r + push/kp = -0.1167 m`); energy residual `<= 1e-12 J`. |

### Material Card Disclosure (`two_body_conduction.fsim`)

The bulk heat capacities use elemental aluminum from NBS Circular C447 Table 67
(`aluminum-pure-nbs-c447`), whereas the vacuum contact-resistance curve is
Atkins (1964) Table II runs 23-28 measured on lathe-cut 6061-T6 aluminum alloy
faces (`6061-T6 aluminum alloy`). Both parts are aluminum so only the alloy
differs between bulk and surface cards; this mismatch is disclosed rather than
masked under a single material label.

**One-way film (`shear_film_pair.fsim`):** `mu` is read once at the declared `313.15 K`
and held; it does not follow the film's temperature, and the receipt marks that key `held`.

**Material card & model disclosures (`block_on_belt.fsim`):** bulk mass is sized to the
500 g dead weight from `steel-music-wire-ligo-p1400229` (`0.05 x 0.05 x 0.02564102564102564 m`
-> `m = 0.5 kg`; the PLA sample's own mass is neglected, and `N = 5 N` follows Section 2.5.1's
printed rounding of `0.5 kg * g = 4.90 N`). Friction values `mu_s = 0.371` and `mu_k = 0.195`
come from Katsivardi et al. (2025) Table 3 (3D circle geometric pattern, not Table 2; `n = 3`,
printed `+/- 0.006` and `+/- 0.003` are undefined in the source -> `UncertaintyModel::Unstated`).
Unstated in Section 2.5.1: counterface metallic-plate alloy/finish and test temperature/humidity
(Section 2.5.2's ball-on-disk conditions are not carried over). Model assumptions: `k = 1000 N/m`
is the device's design choice rather than a tribometer property; `mu_s` is the first peak from rest
(Section 3.4.1) reused at every re-stick; `mu_k` was measured at `2.5 mm/s` and is held constant
(one-way, marked `held` in the receipt) even though relative slip speed peaks near `41.9 mm/s` (`v + sqrt(v^2 + (d*omega)^2) = 16.8 x v`)
with no rate data; and because nothing measured a heat partition between the PLA sample and plate,
neither side declares thermal physics so all friction heat routes to ambient
(`default: no thermal side; all heat to ambient`).

## Reproducing the Tracked Card Packs

The material and interface packs under `examples/assemblies/packs/` are built deterministically
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

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-transport-lemmon2004/manifest.tsv \
  --models data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/manifest.tsv \
  --species data/matdb/seed-v1/nitrogen/manifest.tsv \
  --material-state N2 gas dilute-ideal-gas 0 \
  --out examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/katsivardi-2025-pla-metal-dry-air-interface/manifest.tsv \
  --out examples/assemblies/packs/katsivardi-2025-pla-metal-dry-air-interface.fsintpk

cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/rouleau-2015-self-adhesive-rubber-shear-interface/manifest.tsv \
  --out examples/assemblies/packs/rouleau-2015-self-adhesive-rubber-shear.fsintpk
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
  limit (illustrative)), `gas_piston.fsim` (`ruchardt_freq >= 11.5 Hz` on modal
  `modal.frequency_hz` and `gas_temp_max <= 310 K` on transient `gas_temp`,
  declared limits (illustrative)), `block_on_belt.fsim` (`block_pos_cap <= 0.005 m`,
  declared limit (illustrative)), and all four
  `forced_oscillator.fsim` modes (`pos_limit`, `mode1_freq >= 1.5 Hz`,
  `harm_peak <= 0.25 s` bounding channel `(0, 0)`—the spring's extension per
  unit velocity at its ground-side port—and `eq_extension >= -0.06 m` bounding
  the signed spring state `-0.05 m` under the `+10 N` load, whose sign follows
  the spring's port orientation on `stationary.state.joint.spring_damper.spring.state`,
  which follows the lowered node name; all declared limits (illustrative))—ends
  with a populated eight-term
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
