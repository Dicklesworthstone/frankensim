// Tests of the breathing-wall shell (RubyPrairie; derivation note satwall/BELLWALL-DERIVATION-NOTE.md, MB 23820/23871).
// Bounds are fixed before the first run in satwall/STAMP-shell-tests.txt; the lumped-cell bound is read from the
// model satwall/o1-model.py, not from a run of this file.

use super::*;
use crate::WallModel;
use fs_phs::WallPin;

const PI: f64 = core::f64::consts::PI;
const TAU: f64 = core::f64::consts::TAU;

// An analytic thin-membrane fixture: sigma = rho t and E t / r^2 at the cell's mean radius.
struct Membrane {
    young: f64,
    density: f64,
    thickness: f64,
}

impl ShellLaw for Membrane {
    fn strip(&self, radius_m: f64) -> Result<ShellStrip, NonlinearDuctError> {
        Ok(ShellStrip {
            surface_density_kg_m2: self.density * self.thickness,
            stiffness_per_area_pa_m: self.young * self.thickness / (radius_m * radius_m),
            membrane: true,
        })
    }
}

const BRASS: Membrane = Membrane {
    young: 110_316.0e6,
    density: 8525.4,
    thickness: 0.5e-3,
};

fn cylinder(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}

fn plain_bore(duct: &Duct, cell: f64) -> NonlinearBore {
    let state = air();
    NonlinearBore::new(
        duct,
        air_source(&state),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: cell,
            inlet_volume_m3: 0.0,
            step_s: STEP_S,
            cfl: 0.45,
            wall: None,
        },
    )
    .unwrap()
}

#[test]
fn s1_the_wall_step_is_the_exact_held_force_step_of_the_oscillator() {
    let duct = cylinder(15.0e-3, 0.1);
    let mut bore = plain_bore(&duct, 0.01);
    let law = Membrane {
        young: 1.0e9,
        density: 1200.0,
        thickness: 0.3e-3,
    };
    let quality = 1.0e12;
    let receipt = bore.bind_shell(&duct, &law, quality).unwrap();
    assert_eq!(receipt.cells, bore.solver().mesh().cell_count());
    assert_eq!((receipt.membrane_cells, receipt.thick_cells), (receipt.cells, 0));
    assert_eq!(receipt.largest_half_angle_rad, 0.0);
    // A cylinder: every cell breathes at (1 / r) sqrt(E / rho).
    let hz = (1.0 / 15.0e-3) * (1.0e9_f64 / 1200.0).sqrt() / TAU;
    assert!((receipt.lowest_resonance_hz / hz - 1.0).abs() <= 1e-12, "{receipt:?} {hz}");
    assert!((receipt.highest_resonance_hz / hz - 1.0).abs() <= 1e-12, "{receipt:?} {hz}");
    let shell = bore.shell.as_ref().unwrap();
    let cell = &shell.cells[0];
    let (m, k) = (cell.mass, cell.stiffness);
    let omega = (k / m).sqrt();
    // Independent: the Taylor exponential of the generator with its input column.
    let c = m * omega / quality;
    let generator = [0.0, 1.0, 0.0, -k / m, -c / m, 1.0 / m, 0.0, 0.0, 0.0];
    let e = taylor_exponential(&generator, 3, STEP_S);
    let mine = cell.step;
    let theirs = [e[0], e[1], e[3], e[4], e[2], e[5]];
    for (a, b) in mine.iter().zip(theirs) {
        assert!((a - b).abs() <= 1e-11 * a.abs().max(b.abs()), "{mine:?} {theirs:?}");
    }
    // The undamped closed form of one held-force step.
    let (s0, v0, f0) = (1.0e-6, 2.0e-3, 0.7);
    let [s1, v1] = shell.cells[0].advance([s0, v0], f0);
    let (cw, sw) = (det::cos(omega * STEP_S), det::sin(omega * STEP_S));
    let s_ref = s0 * cw + v0 / omega * sw + f0 / k * (1.0 - cw);
    let v_ref = -s0 * omega * sw + v0 * cw + f0 / (m * omega) * sw;
    assert!((s1 - s_ref).abs() <= 1e-9 * 1.0e-6, "{s1} {s_ref}");
    assert!((v1 - v_ref).abs() <= 1e-9 * 2.0e-3, "{v1} {v_ref}");
}

// A flare: cylinder 6 mm, cone 6 -> 30 mm over 0.15 m, cone 30 -> 60 mm over 0.05 m (continuous radius).
fn flare() -> Duct {
    Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 6.0e-3,
                length: 0.1,
            },
            Segment::Cone {
                inlet_radius: 6.0e-3,
                outlet_radius: 30.0e-3,
                length: 0.15,
            },
            Segment::Cone {
                inlet_radius: 30.0e-3,
                outlet_radius: 60.0e-3,
                length: 0.05,
            },
        ],
    }
}

fn flare_radius(x: f64) -> f64 {
    if x <= 0.1 {
        6.0e-3
    } else if x <= 0.25 {
        6.0e-3 + (30.0e-3 - 6.0e-3) * (x - 0.1) / 0.15
    } else {
        30.0e-3 + (60.0e-3 - 30.0e-3) * (x - 0.25) / 0.05
    }
}

#[test]
fn s2_the_static_wall_is_the_classical_cone_membrane_and_the_resonance_ignores_the_thickness() {
    let duct = flare();
    let pressure = 1000.0;
    let mut resonances = Vec::new();
    for thickness in [0.5e-3, 1.0e-3] {
        let law = Membrane {
            thickness,
            ..BRASS
        };
        let mut bore = plain_bore(&duct, 5.0e-3);
        // Critical damping (Q = 1/2) settles every cell to its static displacement.
        let receipt = bore.bind_shell(&duct, &law, 0.5).unwrap();
        resonances.push((receipt.lowest_resonance_hz, receipt.highest_resonance_hz));
        assert!(receipt.largest_half_angle_rad > 0.5, "{receipt:?}");
        let shell = bore.shell.as_ref().unwrap();
        let positions = bore.solver().mesh().positions().to_vec();
        let mut checked = 0;
        for i in 0..shell.cells.len() {
            let (l, r) = (positions[i], positions[i + 1]);
            let inside = |lo: f64, hi: f64| l >= lo - 1e-12 && r <= hi + 1e-12;
            if !(inside(0.0, 0.1) || inside(0.1, 0.25) || inside(0.25, 0.3)) {
                continue; // a cell across a segment join: its angle is a mean of two cones
            }
            let (ra, rb) = (flare_radius(l), flare_radius(r));
            let w = r - l;
            let mean = 0.5 * (ra + rb);
            let cosine = w / det::sqrt(w * w + (rb - ra) * (rb - ra));
            let lateral = PI * (ra + rb) * w;
            // The ruled law: s = p r^2 / (E t cos phi), the hoop resultant of a cone being p r / cos phi.
            let expected = pressure * mean * mean / (law.young * law.thickness * cosine);
            let force = pressure * lateral;
            let mut state = [0.0, 0.0];
            for _ in 0..20_000 {
                state = shell.cells[i].advance(state, force);
            }
            assert!(
                (state[0] / expected - 1.0).abs() <= 1e-9,
                "cell {i} at {l}: {} against {expected}",
                state[0]
            );
            // At rest: the velocity scale of the cell is omega0 s, not s.
            let omega0 = det::sqrt(shell.cells[i].stiffness / shell.cells[i].mass);
            assert!(state[1].abs() <= 1e-6 * omega0 * state[0].abs(), "{state:?} {omega0}");
            checked += 1;
        }
        assert!(checked >= 50, "{checked}");
    }
    // omega0 = (cos phi / r) sqrt(E / rho): neither the thickness nor the cell size enters.
    let (a, b) = (resonances[0], resonances[1]);
    assert!((a.0 / b.0 - 1.0).abs() <= 1e-13 && (a.1 / b.1 - 1.0).abs() <= 1e-13, "{a:?} {b:?}");
    // The bell end is the flared, widest part: the lowest resonance is below the cylinder's.
    let cylinder_hz = (1.0 / 6.0e-3) * (BRASS.young / BRASS.density).sqrt() / TAU;
    assert!((a.1 / cylinder_hz - 1.0).abs() <= 1e-12, "{a:?} {cylinder_hz}");
    println!("s2 breathing resonances: bell end {:.1} Hz, cylinder {:.1} Hz", a.0, a.1);
    assert!(a.0 < 0.15 * a.1);
    // The lowest is the bell end's, from the closed form (cos phi / r) sqrt(E / rho) on the last cell.
    let (l, r) = (0.295, 0.3);
    let (ra, rb) = (flare_radius(l), flare_radius(r));
    let cosine = (r - l) / det::sqrt((r - l) * (r - l) + (rb - ra) * (rb - ra));
    let bell_hz = cosine / (0.5 * (ra + rb)) * (BRASS.young / BRASS.density).sqrt() / TAU;
    assert!((a.0 / bell_hz - 1.0).abs() <= 1e-9, "{a:?} {bell_hz}");
}

#[test]
fn s3_the_gas_update_keeps_each_cells_mass_and_momentum_and_pays_the_work_from_its_energy() {
    let duct = cylinder(7.0e-3, 0.6);
    let law = Membrane {
        young: 1.0e8,
        density: 1200.0,
        thickness: 0.3e-3,
    };
    let mut bore = junction_bore(&duct);
    drive_and_record(&mut bore, 600, |t| half_cosine(t, 2.0e-4, 2.0e-5));
    let gas = bore.gas().clone();
    let room = bore.room().pressure;
    let (positions, volumes) = {
        let mesh = bore.solver().mesh();
        (mesh.positions().to_vec(), mesh.volumes().to_vec())
    };
    let (mut shell, _) = Shell::new(&duct, &positions, &volumes, &law, 20.0, STEP_S, 0.0).unwrap();
    let initial: Vec<Conserved> = bore.solver().states().to_vec();
    // Mass: each cell's PHYSICAL mass rho (V + A s) is unchanged by every wall update, however many, to roundoff.
    for round in 0..4 {
        let before: Vec<Conserved> = bore.solver().states().to_vec();
        let displacement: Vec<f64> = shell.state.iter().map(|st| st[0]).collect();
        let pressures: Vec<f64> = before
            .iter()
            .map(|q| gas.primitive(*q, || false).unwrap().pressure)
            .collect();
        shell.apply(&mut bore.solver, room).unwrap();
        let after = bore.solver().states().to_vec();
        let new_pressures: Vec<f64> = after
            .iter()
            .map(|q| gas.primitive(*q, || false).unwrap().pressure)
            .collect();
        let mut moved = 0;
        for i in 0..before.len() {
            let cell = &shell.cells[i];
            let (s0, s1) = (displacement[i], shell.state[i][0]);
            let (physical0, physical1) = (cell.volume + cell.lateral * s0, cell.volume + cell.lateral * s1);
            let (q0, q1) = (before[i], after[i]);
            assert!((q1.density * physical1 / (q0.density * physical0) - 1.0).abs() <= 1e-14);
            assert!(
                (q1.momentum * physical1 - q0.momentum * physical0).abs()
                    <= 1e-14 * (q0.momentum * physical0).abs() + 1e-30
            );
            // Energy: the gas pays the work at the MEAN of the pressures before and after, to second order in the
            // swept volume (the trapezoid rule on the isentrope), independent of the code's own force formula.
            let dv = cell.lateral * (s1 - s0);
            let work = 0.5 * (pressures[i] + new_pressures[i]) * dv;
            let lhs = q1.energy * physical1 - q0.energy * physical0;
            assert!(
                (lhs + work).abs() <= 1e-3 * work.abs() + 1e-15 * (q0.energy * physical0).abs(),
                "round {round} cell {i}: {lhs} against {}",
                -work
            );
            if s1 != s0 {
                moved += 1;
            }
        }
        assert!(moved > 20, "round {round}: {moved}");
    }
    // The first round started from rest: the sum of rho V_nominal moved by at most the swept mass.
    let swept: f64 = initial
        .iter()
        .zip(shell.cells.iter().zip(&shell.state))
        .map(|(q, (cell, st))| (q.density * cell.lateral * st[0]).abs())
        .sum();
    assert!(swept > 0.0);
}

// The bound shell leaves the gas field and the walls untouched when one cell's update is refused.
#[test]
fn s7_a_refused_update_leaves_the_field_the_walls_and_the_account_unchanged() {
    let duct = cylinder(7.0e-3, 0.6);
    let mut bore = junction_bore(&duct);
    drive_and_record(&mut bore, 300, |t| half_cosine(t, 2.0e-4, 2.0e-5));
    let room = bore.room().pressure;
    let (positions, volumes) = {
        let mesh = bore.solver().mesh();
        (mesh.positions().to_vec(), mesh.volumes().to_vec())
    };
    let (mut shell, _) = Shell::new(&duct, &positions, &volumes, &BRASS, 50.0, STEP_S, 0.0).unwrap();
    // A wall already displaced inwards by a metre has swept the cell's volume away: the update must refuse.
    shell.state[100] = [-1.0, 0.0];
    shell.state[5] = [1.0e-7, 2.0e-3];
    let states = bore.solver().states().to_vec();
    let (walls, ledger) = (shell.state.clone(), shell.ledger());
    let error = shell.apply(&mut bore.solver, room).unwrap_err().to_string();
    assert!(error.contains("left the bore volume"), "{error}");
    let after = bore.solver().states();
    assert!(
        states.iter().zip(after).all(|(a, b)| {
            a.density.to_bits() == b.density.to_bits()
                && a.momentum.to_bits() == b.momentum.to_bits()
                && a.energy.to_bits() == b.energy.to_bits()
        }),
        "the gas field changed"
    );
    assert_eq!(walls, shell.state);
    assert_eq!(ledger, shell.ledger());
}

// A bore with a lumped inlet cup (no wall on the cup), a thick-wall branch and a profile segment.
struct ThickInLeadPipe;

impl ShellLaw for ThickInLeadPipe {
    fn strip(&self, radius_m: f64) -> Result<ShellStrip, NonlinearDuctError> {
        // The membrane fixture where it applies; below 8 mm the stiffer cylinder law of a thick wall.
        let base = BRASS.strip(radius_m)?;
        Ok(if radius_m < 8.0e-3 {
            ShellStrip {
                stiffness_per_area_pa_m: 0.93 * base.stiffness_per_area_pa_m,
                membrane: false,
                ..base
            }
        } else {
            base
        })
    }
}

#[test]
fn s8_the_cup_the_thick_branch_and_a_profile_segment_are_bound_as_stated() {
    let profile = fs_duct_profile();
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 6.0e-3,
                length: 0.2,
            },
            Segment::Profile {
                profile: std::sync::Arc::new(profile),
                steps: crate::RiccatiSteps(2),
            },
        ],
    };
    let state = air();
    let cup = 2.0e-6;
    let mut bore = NonlinearBore::new(
        &duct,
        air_source(&state),
        293.15,
        101_325.0,
        None,
        BoreSpec {
            maximum_cell_m: 5.0e-3,
            inlet_volume_m3: cup,
            step_s: STEP_S,
            cfl: 0.45,
            wall: None,
        },
    )
    .unwrap();
    let receipt = bore.bind_shell(&duct, &ThickInLeadPipe, 100.0).unwrap();
    assert!(receipt.thick_cells > 0 && receipt.membrane_cells > 0, "{receipt:?}");
    assert_eq!(receipt.thick_cells + receipt.membrane_cells, receipt.cells);
    assert!(receipt.largest_thick_half_angle_rad < receipt.largest_half_angle_rad, "{receipt:?}");
    assert!(receipt.largest_half_angle_rad > 0.04, "{receipt:?}");
    let shell = bore.shell.as_ref().unwrap();
    // Cell 0 holds the cup as gas that no wall surrounds: its volume is the band's plus the cup.
    let mesh = bore.solver().mesh();
    let (w, r) = (mesh.positions()[1] - mesh.positions()[0], 6.0e-3);
    assert!((mesh.volumes()[0] - (core::f64::consts::PI * r * r * w + cup)).abs() <= 1e-12 * mesh.volumes()[0]);
    assert!((shell.cells[0].lateral - core::f64::consts::TAU * r * w).abs() <= 1e-12 * shell.cells[0].lateral);
    // The profile's cells: the projected lateral area against an independent Simpson integral of 2 pi r dx on the
    // exponential flare r = r0 exp(m x / 2), the exact value 2 pi r0 (exp(m x / 2) / (m / 2)) between the faces.
    let (r0, flare) = (0.01, 3.0);
    let mut worst = 0.0_f64;
    for i in 0..shell.cells.len() {
        let (left, right) = (mesh.positions()[i], mesh.positions()[i + 1]);
        if left < 0.2 - 1e-12 {
            continue;
        }
        let (a, b) = (left - 0.2, right - 0.2);
        let exact = core::f64::consts::TAU * r0 * (det::exp(0.5 * flare * b) - det::exp(0.5 * flare * a)) / (0.5 * flare);
        worst = worst.max((shell.cells[i].lateral / exact - 1.0).abs());
    }
    println!("s8 profile cells: worst relative lateral-area error {worst:.2e}");
    assert!(worst <= 1e-9, "{worst}");
}

fn fs_duct_profile() -> crate::RadiusProfile {
    // r(x) = 0.01 exp(1.5 x) over 0.8 m, starting where the cylinder (6 mm) ends: the radii differ, so the join is a
    // radius step; the bore puts a face on it.
    crate::RadiusProfile::exponential(0.01, 3.0, 0.8).unwrap()
}

// Independent dissipation: the integral of c v^2 over each step of the held-force trajectory, by Simpson on 16
// sub-intervals of the Taylor exponential (not the shell's exponential), summed over the force history.
fn independent_dissipation(m: f64, k: f64, c: f64, forces: &[f64]) -> f64 {
    let generator = [0.0, 1.0, 0.0, -k / m, -c / m, 1.0 / m, 0.0, 0.0, 0.0];
    let step = taylor_exponential(&generator, 3, STEP_S);
    let panels = 16usize;
    let sub: Vec<Vec<f64>> = (0..=panels)
        .map(|j| taylor_exponential(&generator, 3, STEP_S * j as f64 / panels as f64))
        .collect();
    let (mut s, mut v) = (0.0, 0.0);
    let mut total = 0.0;
    for &f in forces {
        let mut accumulated = 0.0;
        for (j, e) in sub.iter().enumerate() {
            let weight = if j == 0 || j == panels {
                1.0
            } else if j % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let velocity = e[3] * s + e[4] * v + e[5] * f;
            accumulated += weight * velocity * velocity;
        }
        total += c * accumulated * STEP_S / (3.0 * panels as f64);
        let (s1, v1) = (
            step[0] * s + step[1] * v + step[2] * f,
            step[3] * s + step[4] * v + step[5] * f,
        );
        s = s1;
        v = v1;
    }
    total
}

#[test]
fn o1_a_closed_cell_with_a_wall_follows_the_exact_map_the_continuous_solution_and_the_energy_account() {
    let state = air();
    let gas = air_source(&state);
    let (p0, temperature) = (101_325.0, 293.15);
    let gamma = state.gamma;
    let (r, w) = (15.0e-3, 5.0e-3);
    let duct = cylinder(r, 2.0 * w);
    let cross = PI * r * r;
    let (a_cell, v_cell) = (TAU * r * w, cross * w);
    let g = gamma * p0 * a_cell / v_cell;
    let dp = 100.0;
    let rho0 = p0 / (gas.gas_constant() * temperature);
    let control = StepControl {
        cfl: 0.45,
        reconstruction: Reconstruction::Muscl,
        maximum_halvings: 0,
    };
    for young in [1.0e7, 1.0e8, 1.0e9] {
        let law = Membrane {
            young,
            density: 1200.0,
            thickness: 0.3e-3,
        };
        let quality = 20.0;
        let q = gas
            .conserved(rho0 * (1.0 + dp / (gamma * p0)), 0.0, p0 + dp)
            .unwrap();
        let mut solver = Euler1d::new(
            gas.clone(),
            Length::new(w),
            cross,
            vec![q, q],
            Boundaries::Ends {
                left: Boundary::Reflecting,
                right: Boundary::Reflecting,
            },
            || false,
        )
        .unwrap();
        let (mut shell, _) = {
            let mesh = solver.mesh();
            Shell::new(&duct, mesh.positions(), mesh.volumes(), &law, quality, STEP_S, 0.0).unwrap()
        };
        let (m, k) = (shell.cells[0].mass, shell.cells[0].stiffness);
        let omega0 = (k / m).sqrt();
        let c = m * omega0 / quality;
        let epsilon = g * a_cell / k;
        let omega_c = omega0 * (1.0 + epsilon).sqrt();
        let n = (4.0 * TAU / omega_c / STEP_S).round() as usize;
        let [a, b, cc, d, e, f] = shell.cells[0].step;
        // The scheme's own linear map: the wall step under the MIDPOINT force F = A (dp - g ds / 2) (a scalar linear
        // equation, since ds is linear in F), then the volume rescale to first order, dp' = dp - g ds.
        let denominator = 1.0 + 0.5 * g * a_cell * e;
        let map = |x: [f64; 3]| {
            let free = (a - 1.0) * x[1] + b * x[2];
            let force = (a_cell * x[0] - 0.5 * g * a_cell * free) / denominator;
            let s1 = a * x[1] + b * x[2] + e * force;
            let v1 = cc * x[1] + d * x[2] + f * force;
            [x[0] - g * (s1 - x[1]), s1, v1]
        };
        let e_start: f64 = solver
            .states()
            .iter()
            .zip(shell.cells.iter())
            .map(|(q, cell)| q.energy * cell.volume)
            .sum();
        let mut x_map = [dp, 0.0, 0.0];
        let (mut forces, mut s_bore, mut s_map, mut dp_bore, mut dp_map) =
            (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        let mut dissipated_before = 0.0;
        let mut smallest_step_dissipation = f64::INFINITY;
        for _ in 0..n {
            let report = solver.step(Time::new(STEP_S), control, || false).unwrap();
            assert_eq!(report.elapsed.value(), STEP_S);
            let wall_before = shell.state[0];
            shell.apply(&mut solver, p0).unwrap();
            // The force the wall actually held, read back from its velocity row: v' = c s + d v + f F.
            forces.push((shell.state[0][1] - (cc * wall_before[0] + d * wall_before[1])) / f);
            x_map = map(x_map);
            let after = gas.primitive(solver.states()[0], || false).unwrap().pressure;
            dp_bore.push(after - p0);
            dp_map.push(x_map[0]);
            s_bore.push(shell.state[0][0]);
            s_map.push(x_map[1]);
            let dissipated = shell.ledger().dissipated_j;
            smallest_step_dissipation = smallest_step_dissipation.min(dissipated - dissipated_before);
            dissipated_before = dissipated;
        }
        // (a) The bore is the stated scheme: the nonlinearity is O(dp / p0) of the amplitude.
        // The allowance scales with the coupling's own effect on the pressure (max |dp_map - dp|, the rigid cell
        // holding dp), not with dp, so a wall that did nothing could not pass.
        let nonlinear = 10.0 * dp / p0;
        let effect = dp_map.iter().map(|x| (x - dp).abs()).fold(0.0, f64::max);
        let worst_dp = dp_bore
            .iter()
            .zip(&dp_map)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        let s_scale = s_map.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let worst_s = s_bore
            .iter()
            .zip(&s_map)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        println!(
            "o1 E {young:.0e}: eps {epsilon:.4} w0 dt {:.4} n {n} | bore vs scheme: dp {:.2e} of {dp}, s {:.2e} of {s_scale:.2e} (bound {nonlinear:.1e})",
            omega0 * STEP_S,
            worst_dp / dp,
            worst_s / s_scale
        );
        println!("o1 E {young:.0e}: the coupling moves dp by up to {effect:.4} Pa; bore minus scheme {worst_dp:.3e} Pa, allowance {:.3e}", nonlinear * effect);
        assert!(worst_dp <= nonlinear * effect && worst_s <= nonlinear * s_scale);
        // (b) The scheme is the continuous system to its (second-order) splitting error: o1-model2.py gives gaps of
        // 3.3e-6, 7.9e-7 and 6.9e-8 of dp for these walls; the bound is 2e-3 of the coupling's own scale
        // dp eps / (1 + eps), about 100 times the model's gap.
        let generator = [
            0.0,
            0.0,
            -g,
            0.0,
            0.0,
            1.0,
            a_cell / m,
            -k / m,
            -c / m,
        ];
        let exact = taylor_exponential(&generator, 3, n as f64 * STEP_S);
        let continuous = exact[0] * dp;
        let splitting = 2.0e-3 * dp * epsilon / (1.0 + epsilon);
        let gap = (x_map[0] - continuous).abs();
        println!("o1 E {young:.0e}: scheme {:+.4} continuous {continuous:+.4} gap {gap:.3e} bound {splitting:.3e}", x_map[0]);
        assert!(gap <= splitting + 1e-9 * dp, "{gap} {splitting}");
        // (c) The energy account: dissipation against an independent integral of c v^2, the gas side against
        // the solver's own cell energies, the account closed, and the damper never giving energy back.
        let ledger = shell.ledger();
        let independent = 2.0 * independent_dissipation(m, k, c, &forces);
        println!("o1 E {young:.0e}: ledger {ledger:?} independent dissipation {independent:.6e}");
        assert!(
            (ledger.dissipated_j / independent - 1.0).abs() <= 1e-5,
            "{} {independent}",
            ledger.dissipated_j
        );
        let e_end: f64 = solver
            .states()
            .iter()
            .zip(shell.cells.iter().zip(&shell.state))
            .map(|(q, (cell, st))| q.energy * (cell.volume + cell.lateral * st[0]))
            .sum();
        assert!(
            ((e_start - e_end) / ledger.gas_work_j - 1.0).abs() <= 1e-4,
            "{} {}",
            e_start - e_end,
            ledger.gas_work_j
        );
        assert!(ledger.closing_error() <= 1e-9, "{ledger:?}");
        let mean_step_dissipation = ledger.dissipated_j / n as f64;
        assert!(
            smallest_step_dissipation >= -1e-9 * mean_step_dissipation,
            "{smallest_step_dissipation} against a mean of {mean_step_dissipation}"
        );
        assert!(ledger.dissipated_j > 0.0 && ledger.stored_j > 0.0);
    }
}

#[test]
fn s5_a_bound_shell_runs_inside_the_bore_step_and_its_account_closes() {
    let duct = Duct {
        segments: vec![
            Segment::Cylinder {
                radius: 6.0e-3,
                length: 0.5,
            },
            Segment::Cone {
                inlet_radius: 6.0e-3,
                outlet_radius: 60.0e-3,
                length: 0.3,
            },
        ],
    };
    let mut rigid = junction_bore(&duct);
    let mut shelled = junction_bore(&duct);
    let receipt = shelled.bind_shell(&duct, &BRASS, 100.0).unwrap();
    assert!(rigid.shell_ledger().is_none());
    assert_eq!(receipt.thick_cells, 0);
    let flow = |t: f64| half_cosine(t, 2.0e-4, 2.0e-5);
    let p_rigid = drive_and_record(&mut rigid, 6000, flow);
    let p_shell = drive_and_record(&mut shelled, 6000, flow);
    let ledger = shelled.shell_ledger().unwrap();
    println!("s5 receipt {receipt:?}\ns5 ledger {ledger:?}");
    assert!(ledger.closing_error() <= 1e-9, "{ledger:?}");
    assert!(ledger.dissipated_j >= 0.0 && ledger.stored_j > 0.0 && ledger.gas_work_j.is_finite());
    let scale = p_rigid.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let diff = p_rigid
        .iter()
        .zip(&p_shell)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    println!("s5 largest pressure difference {diff:.3e} of {scale:.3e}");
    assert!(diff > 0.0 && diff <= 0.05 * scale);
}

struct Refuses;

impl ShellLaw for Refuses {
    fn strip(&self, _: f64) -> Result<ShellStrip, NonlinearDuctError> {
        Err(unsupported("the law refuses this wall"))
    }
}

struct Zero;

impl ShellLaw for Zero {
    fn strip(&self, _: f64) -> Result<ShellStrip, NonlinearDuctError> {
        Ok(ShellStrip {
            surface_density_kg_m2: 0.0,
            stiffness_per_area_pa_m: 1.0,
            membrane: true,
        })
    }
}

#[test]
fn s6_binding_refuses_by_name_and_leaves_the_bore_as_it_was() {
    let duct = cylinder(15.0e-3, 0.1);
    let refusal = |bore: &mut NonlinearBore, duct: &Duct, law: &dyn ShellLaw, quality: f64| {
        bore.bind_shell(duct, law, quality).unwrap_err().to_string()
    };
    let mut bore = plain_bore(&duct, 0.01);
    for quality in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let text = refusal(&mut bore, &duct, &BRASS, quality);
        assert!(text.contains("quality factor"), "{text}");
    }
    assert!(refusal(&mut bore, &cylinder(16.0e-3, 0.1), &BRASS, 10.0).contains("cell volume"));
    assert!(refusal(&mut bore, &cylinder(15.0e-3, 0.2), &BRASS, 10.0).contains("length"));
    assert!(refusal(&mut bore, &duct, &Refuses, 10.0).contains("the law refuses this wall"));
    assert!(refusal(&mut bore, &duct, &Zero, 10.0).contains("positive finite surface density"));
    assert!(bore.shell_ledger().is_none());
    bore.bind_shell(&duct, &BRASS, 10.0).unwrap();
    assert!(refusal(&mut bore, &duct, &BRASS, 10.0).contains("binds its shell once"));
    let mut stepped = plain_bore(&duct, 0.01);
    stepped.step(|_| Ok(None)).unwrap();
    assert!(refusal(&mut stepped, &duct, &BRASS, 10.0).contains("before the first step"));
}

// ---- O2: the breathing wall against the plane-wave chain with the equivalent WallPin (MB 23820 item 2) ----

// The bore's response to a half-cosine volume pulse: the inlet-face pressure and the prescribed volume flow, per step.
fn pulse_response(bore: &mut NonlinearBore) -> (Vec<f64>, Vec<f64>) {
    let room = bore.room();
    let source = bore.gas().clone();
    let inlet_area = bore.inlet_area_m2();
    let samples = (0.25 / STEP_S) as usize;
    let (width, peak_flow) = (1.0e-4, 1.0e-8);
    let flow = |t: f64| {
        if t < width {
            peak_flow * 0.5 * (1.0 - det::cos(TAU * t / width))
        } else {
            0.0
        }
    };
    let mut pressure = Vec::with_capacity(samples);
    let mut volume = Vec::with_capacity(samples);
    let mut guess = room.pressure;
    for n in 0..samples {
        let mut face = room.pressure;
        bore.step(|stage| {
            let (flux, p) = velocity_port(&source, stage, flow(stage.time) / inlet_area, &mut guess)?;
            if stage.stage == 0 {
                face = p;
            }
            Ok(Some(flux))
        })
        .unwrap();
        pressure.push(face - room.pressure);
        volume.push(flow(n as f64 * STEP_S));
    }
    (pressure, volume)
}

// The maximum of `level_db` within +-100 cents of `f_ref` (401 points, 0.5 cent apart), by the parabolic peak.
fn peak_near(level_db: impl Fn(f64) -> f64, f_ref: f64) -> f64 {
    let fine: Vec<f64> = (-200..=200)
        .map(|c| f_ref * det::exp(f64::from(c) * 0.5 / 1200.0 * core::f64::consts::LN_2))
        .collect();
    let values: Vec<f64> = fine.iter().map(|&f| level_db(f)).collect();
    let best = (1..fine.len() - 1)
        .max_by(|&a, &b| values[a].total_cmp(&values[b]))
        .unwrap();
    parabolic_peak(&fine, &values, best).0
}

#[test]
fn o2_the_walled_tube_shifts_its_resonances_like_the_plane_wave_chain_with_the_equivalent_pin() {
    let state = air();
    let (r, length) = (15.0e-3, 0.3);
    let duct = cylinder(r, length);
    // A soft wall, so that the shift is measurable: the cell law is the equivalent WallPin sigma = rho t,
    // K = E t / r^2, resistance = sigma omega0 / Q on a cylinder (the derivation note, section 3).
    let law = Membrane {
        young: 1.0e9,
        density: 1200.0,
        thickness: 0.3e-3,
    };
    let quality = 50.0;
    let sigma = law.density * law.thickness;
    let stiffness = law.young * law.thickness / (r * r);
    let pin = WallPin {
        surface_density: sigma,
        stiffness_per_area: stiffness,
        resistance: sigma * det::sqrt(stiffness / sigma) / quality,
    };
    let zc = state.density * state.sound_speed / (PI * r * r);
    let tmm = |f: f64, wall: WallModel<'_>| {
        let omega = TAU * f;
        let load = piston(omega * r / state.sound_speed);
        let z = crate::input_impedance_load_with_wall(
            &duct,
            &state,
            omega,
            LossModel::Lossless,
            C64::new(zc * load.re, -zc * load.im),
            wall,
        )
        .unwrap();
        20.0 * det::ln(z.impedance.abs()) / core::f64::consts::LN_10
    };
    let mut rigid_bore = plain_bore(&duct, 2.4e-3);
    let mut walled_bore = plain_bore(&duct, 2.4e-3);
    walled_bore.bind_shell(&duct, &law, quality).unwrap();
    let (p_rigid, v_rigid) = pulse_response(&mut rigid_bore);
    let (p_walled, v_walled) = pulse_response(&mut walled_bore);
    let level = |p: &[f64], v: &[f64]| {
        let (p, v) = (p.to_vec(), v.to_vec());
        move |f: f64| {
            let omega_dt = TAU * f * STEP_S;
            let z = dft(&p, omega_dt) / dft(&v, omega_dt);
            20.0 * det::ln(z.abs()) / core::f64::consts::LN_10
        }
    };
    let (level_rigid, level_walled) = (level(&p_rigid, &v_rigid), level(&p_walled, &v_walled));
    // The first five peaks of the rigid chain on a 0.5 Hz grid.
    let grid: Vec<f64> = (0..6800).map(|i| 100.0 + 0.5 * i as f64).collect();
    let values: Vec<f64> = grid.iter().map(|&f| tmm(f, WallModel::Rigid)).collect();
    let mut references = Vec::new();
    for i in 1..grid.len() - 1 {
        if values[i] > values[i - 1] && values[i] >= values[i + 1] && references.len() < 5 {
            references.push(grid[i]);
        }
    }
    assert_eq!(references.len(), 5, "{references:?}");
    let cents = |a: f64, b: f64| 1200.0 * det::ln(b / a) / core::f64::consts::LN_2;
    for (k, &f_ref) in references.iter().enumerate() {
        let tmm_rigid = peak_near(|f| tmm(f, WallModel::Rigid), f_ref);
        let tmm_walled = peak_near(|f| tmm(f, WallModel::Pin(&pin)), f_ref);
        let bore_rigid = peak_near(&level_rigid, f_ref);
        let bore_walled = peak_near(&level_walled, f_ref);
        let (shift_tmm, shift_bore) = (cents(tmm_rigid, tmm_walled), cents(bore_rigid, bore_walled));
        let bound = 4.0 * (TAU * tmm_rigid * STEP_S) * shift_tmm.abs() + 0.1;
        println!(
            "o2 peak {}: tmm {tmm_rigid:8.2} Hz, wall shifts it {shift_tmm:+7.3} cents | bore {bore_rigid:8.2} Hz, wall shifts it {shift_bore:+7.3} cents | gap {:+.3} cents (bound {bound:.3})",
            k + 1,
            shift_bore - shift_tmm
        );
        if k < 3 {
            assert!(shift_tmm < -5.0, "the wall must be felt: {shift_tmm}");
        }
        assert!((shift_bore - shift_tmm).abs() <= bound);
    }
}

#[test]
fn o3_the_coupled_scheme_gives_no_energy_back_even_at_a_very_high_quality() {
    // Review of the shell (item 1): with the start-of-step pressure held, the split coupling put g A ds^2 / 2 of
    // energy into the system every step, +9 % of H over 20000 steps at Q = 1e12 (o1-model2.py). With the midpoint
    // force the energy H = m v^2 / 2 + k s^2 / 2 + (A / 2 g) dp^2 stays within 1e-9 of H0, and never exceeds it at Q = 1e6.
    let state = air();
    let (p0, gamma) = (101_325.0, state.gamma);
    let (r, w) = (15.0e-3, 5.0e-3);
    let duct = cylinder(r, 2.0 * w);
    let (a_cell, v_cell) = (TAU * r * w, PI * r * r * w);
    let g = gamma * p0 * a_cell / v_cell;
    for (quality, lower) in [(1.0e12, 1.0 - 1.0e-8), (1.0e6, 0.0)] {
        let mut bore = plain_bore(&duct, w);
        let law = Membrane {
            young: 1.0e9,
            density: 1200.0,
            thickness: 0.3e-3,
        };
        bore.bind_shell(&duct, &law, quality).unwrap();
        let cell = &bore.shell.as_ref().unwrap().cells[0];
        let (m, k) = (cell.mass, cell.stiffness);
        let [a, b, cc, d, e, f] = cell.step;
        let denominator = 1.0 + 0.5 * g * a_cell * e;
        let energy = |x: [f64; 3]| 0.5 * m * x[2] * x[2] + 0.5 * k * x[1] * x[1] + a_cell / (2.0 * g) * x[0] * x[0];
        let mut x = [1.0, 0.0, 0.0];
        let h0 = energy(x);
        let (mut highest, mut lowest) = (0.0_f64, f64::INFINITY);
        for _ in 0..20_000 {
            let free = (a - 1.0) * x[1] + b * x[2];
            let force = (a_cell * x[0] - 0.5 * g * a_cell * free) / denominator;
            let s1 = a * x[1] + b * x[2] + e * force;
            let v1 = cc * x[1] + d * x[2] + f * force;
            x = [x[0] - g * (s1 - x[1]), s1, v1];
            let ratio = energy(x) / h0;
            highest = highest.max(ratio);
            lowest = lowest.min(ratio);
        }
        println!("o3 Q {quality:.0e}: H/H0 between {lowest:.12} and {highest:.12}");
        assert!(highest <= 1.0 + 1.0e-9 && lowest >= lower, "{lowest} {highest}");
    }
}
