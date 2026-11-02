use super::*;
use crate::dynamics::tests::gas;
use fs_qty::Length;

const AREA: f64 = 0.7;
const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn field(n: usize, cp: f64, formation: f64) -> Euler1d {
    let gas = gas(cp, formation);
    let q = gas.conserved(1.0, 0.0, 1.0).unwrap();
    Euler1d::new(
        gas,
        Length::new(1.0 / n as f64),
        AREA,
        vec![q; n],
        CLOSED,
        || false,
    )
    .unwrap()
}

fn controls() -> PistonControl {
    PistonControl {
        midpoint_velocity_bracket_m_s: [-0.5, 0.5],
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 64,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g0_balanced_pressure_keeps_a_free_piston_and_field_at_rest() {
    let mut field = field(24, 0.4, 0.0);
    let original = field.states().to_vec();
    let mesh = field.mesh().clone();
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    for _ in 0..8 {
        let report = field
            .step_with_piston(
                &mut piston,
                Time::new(0.003),
                StepControl::default(),
                PistonLoad {
                    exterior_pressure_pa: 1.0,
                    external_force_n: 0.0,
                },
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(piston.velocity_m_s(), 0.0);
        assert_eq!(field.mesh(), &mesh);
        assert_eq!(report.displacement_m, 0.0);
        assert_eq!(report.kinetic_energy_change_j, 0.0);
        assert_eq!(report.gas.boundary_change.energy, 0.0);
        assert_eq!(report.gas_trials, 1);
        close(
            report.gas_impulse_n_s,
            AREA * report.gas.elapsed.value(),
            1e-14,
        );
        close(
            report.fixed_wall_impulse_n_s,
            -AREA * report.gas.elapsed.value(),
            1e-14,
        );
        assert!(report.admitted(controls()));
    }
    for (&q, &initial) in field.states().iter().zip(&original) {
        close(q.density, initial.density, 3e-14);
        close(q.momentum, initial.momentum, 3e-14);
        close(q.energy, initial.energy, 3e-14);
    }
}

#[test]
fn g1_pressure_drives_compression_with_gas_wall_load_and_support_balances() {
    let mut ordinary = field(32, 0.4, 0.0);
    let mut shifted = field(32, 0.4, 1e18);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut other = piston;
    let before = ordinary.totals();
    let load = PistonLoad {
        exterior_pressure_pa: 1.8,
        external_force_n: -0.15,
    };
    let mut external_work = 0.0;
    let mut external_impulse = 0.0;
    let mut support_impulse = 0.0;
    let mut steps = 0;
    while ordinary.time().value() < 0.12 {
        let duration = Time::new(0.12 - ordinary.time().value());
        let report = ordinary
            .step_with_piston(
                &mut piston,
                duration,
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        let reference = shifted
            .step_with_piston(
                &mut other,
                duration,
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report, reference);
        assert_eq!(piston, other);
        assert_eq!(ordinary.states(), shifted.states());
        assert_eq!(ordinary.mesh(), shifted.mesh());
        assert_eq!(report.gas.boundary_change.mass, 0.0);
        assert_eq!(report.gas.left_boundary_change.energy, 0.0);
        assert!(report.gas.right_boundary_change.energy > 0.0);
        assert!(report.displacement_m < 0.0);
        assert!(report.gas.maximum_cfl <= 0.40000000000001);
        assert!(report.admitted(controls()));
        close(
            report.midpoint_velocity_m_s,
            f64::midpoint(report.velocity_m_s[0], report.velocity_m_s[1]),
            2e-15,
        );
        close(
            report.gas.right_boundary_change.energy,
            -report.midpoint_velocity_m_s * report.gas_impulse_n_s,
            2e-14,
        );
        external_work += report.external_force_work_j - report.exterior_pressure_work_j;
        external_impulse +=
            report.gas.elapsed.value() * (load.external_force_n - AREA * load.exterior_pressure_pa);
        support_impulse += report.fixed_wall_impulse_n_s;
        steps += 1;
        assert!(steps < 500);
    }
    let after = ordinary.totals();
    let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    close(after.mass, before.mass, 3e-12);
    close(after.energy - before.energy + kinetic, external_work, 3e-10);
    close(
        after.momentum - before.momentum + piston.mass_kg() * piston.velocity_m_s(),
        external_impulse - support_impulse,
        3e-10,
    );
    assert!(piston.velocity_m_s() < -0.05);
    let left = ordinary
        .gas()
        .primitive(ordinary.states()[0], || false)
        .unwrap();
    let right = ordinary
        .gas()
        .primitive(*ordinary.states().last().unwrap(), || false)
        .unwrap();
    assert!(right.pressure > 1.05);
    assert!(
        right.pressure > left.pressure + 0.03,
        "the field must resolve the compression wave"
    );
    assert_eq!(ordinary.mesh().positions()[0], 0.0);
    assert!(ordinary.mesh().positions()[32] < 1.0);
    close(
        ordinary.mesh().volumes().iter().sum(),
        AREA * ordinary.mesh().positions()[32],
        5e-14,
    );
}

#[test]
fn g0_accepted_field_replays_the_final_solved_wall_speed_and_duration() {
    let mut coupled = field(16, 0.4, 0.0);
    let mut prescribed = field(16, 0.4, 0.0);
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let report = coupled
        .step_with_piston(
            &mut piston,
            Time::new(0.004),
            StepControl::default(),
            PistonLoad {
                exterior_pressure_pa: 2.0,
                external_force_n: 0.1,
            },
            controls(),
            || false,
        )
        .unwrap();
    assert!(report.gas_trials > 3);
    let replay = prescribed
        .step_with_moving_walls(
            report.gas.elapsed,
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            [0.0, report.midpoint_velocity_m_s],
            || false,
        )
        .unwrap();
    assert_eq!(report.gas.elapsed, replay.elapsed);
    assert_eq!(coupled.mesh(), prescribed.mesh());
    assert_eq!(coupled.time(), prescribed.time());
    for (&q, &reference) in coupled.states().iter().zip(prescribed.states()) {
        close(q.density, reference.density, 2e-12);
        close(q.momentum, reference.momentum, 2e-12);
        close(q.energy, reference.energy, 2e-12);
    }
    close(
        report.gas.right_boundary_change.momentum,
        replay.right_boundary_change.momentum,
        2e-13,
    );
    close(
        report.gas.right_boundary_change.energy,
        replay.right_boundary_change.energy,
        2e-13,
    );
}

#[test]
fn g3_refinement_approaches_independent_acoustic_piston_impedance_solution() {
    fn error(n: usize) -> f64 {
        let mut flow = field(n, 0.0, 0.0);
        let mass = 0.2;
        let mut piston = Piston::new(mass, 0.0).unwrap();
        let excess = 1e-4;
        let end = 0.2;
        let control = PistonControl {
            midpoint_velocity_bracket_m_s: [-0.02, 0.02],
            momentum_tolerance_n_s: 1e-14,
            energy_tolerance_j: 1e-11,
            ..controls()
        };
        let mut steps = 0;
        while flow.time().value() < end {
            flow.step_with_piston(
                &mut piston,
                Time::new(end - flow.time().value()),
                StepControl::default(),
                PistonLoad {
                    exterior_pressure_pa: 1.0 + excess,
                    external_force_n: 0.0,
                },
                control,
                || false,
            )
            .unwrap();
            steps += 1;
            assert!(steps < 1000);
        }
        // Before a reflection can return from the opposite wall, independent
        // linear acoustics gives p_wall-p0 = -rho*c*v and
        // M*dv/dt = -A*Delta_p - A*rho*c*v. Its step-load solution has no
        // production EOS, wall closure, spatial operator or root solve in it.
        let impedance = 1.4_f64.sqrt();
        let exact = -(excess / impedance) * (1.0 - (-AREA * impedance * end / mass).exp());
        (piston.velocity_m_s() - exact).abs() / exact.abs()
    }
    let errors = [error(16), error(32), error(64)];
    assert!(errors[1] < 0.9 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.9 * errors[1], "{errors:?}");
    assert!(errors[2] < 0.03, "{errors:?}");
}

#[test]
fn g4_cancellation_bracket_source_and_geometry_refusal_roll_back_every_owner() {
    let load = PistonLoad {
        exterior_pressure_pa: 1.8,
        external_force_n: 0.0,
    };
    let duration = Time::new(0.003);
    let mut baseline = field(12, 0.4, 0.0);
    let mut wall = Piston::new(0.4, 0.0).unwrap();
    let mut polls = 0;
    baseline
        .step_with_piston(
            &mut wall,
            duration,
            StepControl::default(),
            load,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, 8, polls / 2, polls] {
        let mut flow = field(12, 0.4, 0.0);
        let before = flow.states().to_vec();
        let mesh = flow.mesh().clone();
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let original = piston;
        let mut seen = 0;
        let result = flow.step_with_piston(
            &mut piston,
            duration,
            StepControl::default(),
            load,
            controls(),
            || {
                seen += 1;
                seen == cutoff
            },
        );
        assert_eq!(result, Err(PistonError::Cancelled));
        assert_eq!(seen, cutoff, "cancellation must not become a retry");
        assert_eq!(flow.states(), before);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(piston, original);
    }
    for bracket in [[-1e-9, 1e-9], [-100.0, 100.0]] {
        let mut flow = field(12, 0.4, 0.0);
        let before = flow.states().to_vec();
        let mesh = flow.mesh().clone();
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let original = piston;
        let result = flow.step_with_piston(
            &mut piston,
            duration,
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            load,
            PistonControl {
                midpoint_velocity_bracket_m_s: bracket,
                ..controls()
            },
            || false,
        );
        assert!(result.is_err(), "{result:?}");
        assert_eq!(flow.states(), before);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(piston, original);
    }
    let mut flow = field(12, 0.4, 0.0);
    flow.boundaries = Boundaries::Periodic;
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    assert!(matches!(
        flow.step_with_piston(
            &mut piston,
            duration,
            StepControl::default(),
            load,
            controls(),
            || false
        ),
        Err(PistonError::Invalid(_))
    ));
    assert_eq!(flow.time(), Time::new(0.0));
    assert_eq!(piston.velocity_m_s(), 0.0);
}
