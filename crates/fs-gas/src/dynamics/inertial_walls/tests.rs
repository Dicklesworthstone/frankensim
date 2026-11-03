use super::*;
use crate::dynamics::tests::gas;
use fs_qty::Length;

const AREA: f64 = 0.7;
const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn field(boost: f64, nonuniform: bool) -> Euler1d {
    // Original manufactured NASA-9 source, R_specific=1, Cp/R=3.5+0.4*T.
    // This is deliberately not a constant-gamma replacement of the field.
    let source = gas(0.4, 0.0);
    let n = 20;
    let states = (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            let temperature = if nonuniform { 0.8 + 0.4 * x } else { 1.0 };
            source.conserved(1.0, boost, temperature).unwrap()
        })
        .collect();
    Euler1d::new(
        source,
        Length::new(1.0 / n as f64),
        AREA,
        states,
        CLOSED,
        || false,
    )
    .unwrap()
}

fn controls(boost: f64) -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[boost - 1.0, boost + 1.0]; 2],
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 32,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn momentum(flow: &Euler1d, walls: [Piston; 2]) -> f64 {
    flow.totals().momentum
        + walls
            .iter()
            .map(|w| w.mass_kg() * w.velocity_m_s())
            .sum::<f64>()
}
fn energy(flow: &Euler1d, walls: [Piston; 2]) -> f64 {
    flow.totals().energy
        + walls
            .iter()
            .map(|w| 0.5 * w.mass_kg() * w.velocity_m_s().powi(2))
            .sum::<f64>()
}

#[test]
fn balanced_exterior_pressure_preserves_both_resting_walls() {
    let mut flow = field(0.0, false);
    let states = flow.states().to_vec();
    let mesh = flow.mesh().clone();
    let mut walls = [
        Piston::new(0.4, 0.0).unwrap(),
        Piston::new(0.9, 0.0).unwrap(),
    ];
    let original = walls;
    let report = flow
        .step_with_inertial_walls(
            &mut walls,
            Time::new(0.003),
            StepControl::default(),
            [InertialWallLoad {
                exterior_pressure_pa: 1.0,
                external_force_n: 0.0,
            }; 2],
            controls(0.0),
            || false,
        )
        .unwrap();
    assert_eq!(walls, original);
    assert_eq!(flow.mesh(), &mesh);
    assert_eq!(report.displacement_m, [0.0; 2]);
    assert_eq!(report.gas_trials, 1);
    assert_eq!(report.gas.boundary_change.energy, 0.0);
    assert!(report.gas_impulse_n_s[0] < 0.0 && report.gas_impulse_n_s[1] > 0.0);
    for (actual, original) in flow.states().iter().zip(states) {
        close(actual.density, original.density, 3e-14);
        close(actual.momentum, original.momentum, 3e-14);
        close(actual.energy, original.energy, 3e-14);
    }
}

#[test]
fn symmetric_release_recoils_both_walls_and_conserves_closed_system() {
    let mut flow = field(0.0, false);
    let mut walls = [Piston::new(0.4, 0.0).unwrap(); 2];
    let initial_mass = flow.totals().mass;
    let initial_energy = energy(&flow, walls);
    for _ in 0..12 {
        let before = walls;
        let report = flow
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.003),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                controls(0.0),
                || false,
            )
            .unwrap();
        assert!(report.admitted(controls(0.0)));
        assert_eq!(report.gas.boundary_change.mass, 0.0);
        assert!(report.displacement_m[0] < 0.0 && report.displacement_m[1] > 0.0);
        for (i, gas_change) in [
            report.gas.left_boundary_change,
            report.gas.right_boundary_change,
        ]
        .into_iter()
        .enumerate()
        {
            close(
                walls[i].mass_kg() * (walls[i].velocity_m_s() - before[i].velocity_m_s()),
                -gas_change.momentum,
                2e-12,
            );
            close(
                gas_change.energy,
                report.midpoint_velocity_m_s[i] * gas_change.momentum,
                3e-14,
            );
        }
        close(walls[0].velocity_m_s(), -walls[1].velocity_m_s(), 2e-11);
        close(momentum(&flow, walls), 0.0, 3e-11);
        close(energy(&flow, walls), initial_energy, 2e-10);
        close(flow.totals().mass, initial_mass, 2e-13);
    }
    assert!(walls[0].velocity_m_s() < -0.02);
    assert!(walls[1].velocity_m_s() > 0.02);
}

#[test]
fn unequal_walls_and_nonuniform_pressure_close_actual_external_impulse_and_work() {
    let mut flow = field(0.07, true);
    let mut walls = [
        Piston::new(0.25, 0.07).unwrap(),
        Piston::new(1.1, 0.07).unwrap(),
    ];
    let loads = [
        InertialWallLoad {
            exterior_pressure_pa: 0.2,
            external_force_n: -0.03,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.4,
            external_force_n: 0.06,
        },
    ];
    let initial_p = momentum(&flow, walls);
    let initial_e = energy(&flow, walls);
    let mut external_impulse = 0.0;
    let mut external_work = 0.0;
    for _ in 0..15 {
        let old_faces = [
            flow.mesh().positions()[0],
            *flow.mesh().positions().last().unwrap(),
        ];
        let report = flow
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.003),
                StepControl::default(),
                loads,
                controls(0.0),
                || false,
            )
            .unwrap();
        let faces = [
            flow.mesh().positions()[0],
            *flow.mesh().positions().last().unwrap(),
        ];
        for i in 0..2 {
            let force =
                loads[i].external_force_n + [1.0, -1.0][i] * AREA * loads[i].exterior_pressure_pa;
            external_impulse += force * report.gas.elapsed.value();
            external_work += force * (faces[i] - old_faces[i]);
        }
        close(momentum(&flow, walls) - initial_p, external_impulse, 5e-11);
        close(energy(&flow, walls) - initial_e, external_work, 3e-10);
    }
    assert!(walls[0].velocity_m_s() < 0.07);
    assert!(walls[1].velocity_m_s() > 0.07);
    // The resolved gas itself carries momentum; a two-wall reduced-mass
    // oscillator that silently removes this field cannot satisfy this case.
    assert!((flow.totals().momentum - 0.07 * AREA).abs() > 1e-4);
}

#[test]
fn common_boost_retains_resolved_nonuniform_pressure_and_relative_wall_motion() {
    let boost = 0.37;
    let mut ordinary = field(0.0, true);
    let mut boosted = field(boost, true);
    let mut walls = [
        Piston::new(0.25, 0.0).unwrap(),
        Piston::new(1.1, 0.0).unwrap(),
    ];
    let mut other = [
        Piston::new(0.25, boost).unwrap(),
        Piston::new(1.1, boost).unwrap(),
    ];
    // Conserved-variable MUSCL limiting need not be Galilean invariant; this
    // test isolates the original first-order ALE flux and two-wall mechanics.
    let gas_controls = StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    };
    for _ in 0..6 {
        ordinary
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.003),
                gas_controls,
                [InertialWallLoad::default(); 2],
                controls(0.0),
                || false,
            )
            .unwrap();
        boosted
            .step_with_inertial_walls(
                &mut other,
                Time::new(0.003),
                gas_controls,
                [InertialWallLoad::default(); 2],
                controls(boost),
                || false,
            )
            .unwrap();
        close(ordinary.time().value(), boosted.time().value(), 3e-14);
        for i in 0..2 {
            close(
                other[i].velocity_m_s() - boost,
                walls[i].velocity_m_s(),
                2e-9,
            );
        }
        for (&x, &y) in ordinary
            .mesh()
            .positions()
            .iter()
            .zip(boosted.mesh().positions())
        {
            close(y - x, boost * ordinary.time().value(), 2e-11);
        }
        for (q, b) in ordinary.states().iter().zip(boosted.states()) {
            close(b.density, q.density, 2e-9);
            close(b.momentum, q.momentum + boost * q.density, 2e-9);
            close(
                b.energy,
                q.energy + boost * q.momentum + 0.5 * boost * boost * q.density,
                3e-9,
            );
        }
    }
}

#[test]
fn massive_left_wall_approaches_original_fixed_support_piston() {
    let mut original = field(0.0, false);
    let mut two_walls = field(0.0, false);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut walls = [Piston::new(1e10, 0.0).unwrap(), piston];
    let controls = controls(0.0);
    for _ in 0..5 {
        original
            .step_with_piston(
                &mut piston,
                Time::new(0.003),
                StepControl::default(),
                PistonLoad::default(),
                PistonControl {
                    midpoint_velocity_bracket_m_s: [-1.0, 1.0],
                    momentum_tolerance_n_s: controls.momentum_tolerance_n_s,
                    energy_tolerance_j: controls.energy_tolerance_j,
                    maximum_iterations: 64,
                },
                || false,
            )
            .unwrap();
        two_walls
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.003),
                StepControl::default(),
                [InertialWallLoad::default(); 2],
                controls,
                || false,
            )
            .unwrap();
        close(walls[1].velocity_m_s(), piston.velocity_m_s(), 2e-10);
        close(walls[0].velocity_m_s(), 0.0, 2e-11);
    }
}

#[test]
fn cancellation_and_failed_solve_preserve_both_walls_field_geometry_and_clock() {
    let duration = Time::new(0.003);
    let mut baseline = field(0.0, true);
    let original = [
        Piston::new(0.25, 0.0).unwrap(),
        Piston::new(1.1, 0.0).unwrap(),
    ];
    let mut walls = original;
    let mut polls = 0;
    baseline
        .step_with_inertial_walls(
            &mut walls,
            duration,
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut flow = field(0.0, true);
        let states = flow.states().to_vec();
        let mesh = flow.mesh().clone();
        let mut walls = original;
        let mut seen = 0;
        let result = flow.step_with_inertial_walls(
            &mut walls,
            duration,
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || {
                seen += 1;
                seen == cutoff
            },
        );
        assert_eq!(result, Err(PistonError::Cancelled));
        assert_eq!(seen, cutoff);
        assert_eq!(flow.states(), states);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(walls, original);
    }
    let mut flow = field(0.0, true);
    let states = flow.states().to_vec();
    let mesh = flow.mesh().clone();
    let mut walls = original;
    let result = flow.step_with_inertial_walls(
        &mut walls,
        duration,
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        [InertialWallLoad::default(); 2],
        InertialWallsControl {
            midpoint_velocity_bracket_m_s: [[-1e-9, 1e-9]; 2],
            ..controls(0.0)
        },
        || false,
    );
    assert!(
        matches!(result, Err(PistonError::StepRetries { .. })),
        "{result:?}"
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.mesh(), &mesh);
    assert_eq!(flow.time(), Time::new(0.0));
    assert_eq!(walls, original);
    flow.boundaries = Boundaries::Periodic;
    assert!(matches!(
        flow.step_with_inertial_walls(
            &mut walls,
            duration,
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(0.0),
            || false
        ),
        Err(PistonError::Invalid(_))
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(walls, original);
}
