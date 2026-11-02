use super::*;
use crate::dynamics::tests::gas;
use fs_qty::Length;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn flow(n: usize, cp: f64, formation: f64, velocity: impl Fn(f64) -> f64) -> Euler1d {
    let gas = gas(cp, formation);
    let states = (0..n)
        .map(|i| {
            gas.conserved(1.0, velocity((i as f64 + 0.5) / n as f64), 1.0)
                .unwrap()
        })
        .collect();
    Euler1d::new(
        gas,
        Length::new(1.0 / n as f64),
        0.7,
        states,
        CLOSED,
        || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g0_moving_face_flux_matches_lab_reynolds_transport_and_galilean_transform() {
    let gas = gas(0.4, 0.0);
    let (q, p) = gas.state_at(1.2, 0.3, 1.7).unwrap();
    let fixed = gas.interface(q, p, q, p, &mut || false).unwrap();
    for speed in [-0.7, 0.0, 1.3] {
        let moving = gas
            .moving_interface(q, p, q, p, speed, &mut || false)
            .unwrap();
        for k in 0..3 {
            close(
                moving.array()[k],
                fixed.array()[k] - speed * q.array()[k],
                2e-12,
            );
        }
        close(
            moving.signal_speed,
            (p.velocity - speed).abs() + p.sound_speed,
            2e-14,
        );
        if speed == 0.0 {
            assert_eq!(moving, fixed);
        }
        let boost = 2.0;
        let (boosted, primitive) = gas
            .state_at(p.density, p.velocity + boost, p.pressure)
            .unwrap();
        let translated = gas
            .moving_interface(
                boosted,
                primitive,
                boosted,
                primitive,
                speed + boost,
                &mut || false,
            )
            .unwrap();
        close(translated.mass, moving.mass, 2e-14);
        close(
            translated.momentum,
            moving.momentum + boost * moving.mass,
            2e-13,
        );
        close(
            translated.energy,
            moving.energy + boost * moving.momentum + 0.5 * boost * boost * moving.mass,
            2e-12,
        );
    }
}

#[test]
fn g0_zero_speed_retains_stationary_field_mesh_report_and_clock() {
    let mut original = flow(24, 0.4, 0.0, |x| 0.1 * (6.0 * x).sin());
    let mut moving = flow(24, 0.4, 0.0, |x| 0.1 * (6.0 * x).sin());
    for _ in 0..6 {
        let expected = original
            .step(Time::new(0.003), StepControl::default(), || false)
            .unwrap();
        let actual = moving
            .step_with_moving_walls(Time::new(0.003), StepControl::default(), [0.0; 2], || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(moving.states(), original.states());
        assert_eq!(moving.mesh(), original.mesh());
        assert_eq!(moving.time(), original.time());
        assert!(moving.moving.is_none());
    }
}

#[test]
fn g0_comoving_nonuniform_contact_field_preserves_mass_energy_and_geometry() {
    let gas = gas(0.4, 0.0);
    let n = 24;
    let speed = 0.125;
    let positions: Vec<_> = (0..=n).map(|i| (i as f64 / n as f64).powi(2)).collect();
    let volumes: Vec<_> = positions.windows(2).map(|x| 0.7 * (x[1] - x[0])).collect();
    let mesh = Mesh1d::new(positions.clone(), vec![0.7; n + 1], volumes.clone()).unwrap();
    let states: Vec<_> = (0..n)
        .map(|i| {
            gas.conserved(0.8 + 0.4 * i as f64 / n as f64, speed, 1.0)
                .unwrap()
        })
        .collect();
    let original = states.clone();
    let mut flow = Euler1d::new_on_mesh(gas, mesh, states, CLOSED, || false).unwrap();
    for _ in 0..32 {
        let report = flow
            .step_with_moving_walls(Time::new(0.002), StepControl::default(), [speed; 2], || {
                false
            })
            .unwrap();
        assert!(report.maximum_cfl <= 0.40000000000001);
        close(report.boundary_change.mass, 0.0, 0.0);
        close(report.conservation_residual.energy, 0.0, 2e-13);
        assert_eq!(flow.mesh().volumes(), volumes);
    }
    for (&x, &start) in flow.mesh().positions().iter().zip(&positions) {
        close(x, start + speed * flow.time().value(), 3e-14);
    }
    for (&actual, &expected) in flow.states().iter().zip(&original) {
        for k in 0..3 {
            close(actual.array()[k], expected.array()[k], 3e-11);
        }
    }
}

#[test]
fn g1_resolved_piston_compression_closes_pressure_work_and_formation_gauge() {
    let mut ordinary = flow(48, 0.4, 0.0, |_| 0.0);
    let mut shifted = flow(48, 0.4, 1e18, |_| 0.0);
    let before = ordinary.totals();
    let mut work = 0.0;
    let mut impulse = 0.0;
    let mut steps = 0;
    while ordinary.time().value() < 0.12 {
        let duration = Time::new(0.12 - ordinary.time().value());
        let report = ordinary
            .step_with_moving_walls(duration, StepControl::default(), [0.2, 0.0], || false)
            .unwrap();
        let other = shifted
            .step_with_moving_walls(duration, StepControl::default(), [0.2, 0.0], || false)
            .unwrap();
        assert_eq!(report, other);
        assert_eq!(ordinary.states(), shifted.states());
        assert_eq!(ordinary.mesh(), shifted.mesh());
        assert!(report.left_boundary_change.energy > 0.0);
        assert_eq!(report.right_boundary_change.energy, 0.0);
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.wall_heat, 0.0);
        assert!(report.maximum_cfl <= 0.40000000000001);
        close(report.conservation_residual.mass, 0.0, 2e-13);
        close(report.conservation_residual.energy, 0.0, 2e-12);
        work += report.boundary_change.energy;
        impulse += report.boundary_change.momentum;
        steps += 1;
        assert!(steps < 1000);
    }
    let after = ordinary.totals();
    close(after.mass, before.mass, 3e-12);
    close(after.energy - before.energy, work, 3e-12);
    close(after.momentum - before.momentum, impulse, 3e-12);
    close(
        ordinary.mesh().positions()[0],
        0.2 * ordinary.time().value(),
        3e-14,
    );
    close(ordinary.mesh().positions()[48], 1.0, 0.0);
    close(
        ordinary.mesh().volumes().iter().sum(),
        0.7 * (1.0 - 0.2 * 0.12),
        3e-14,
    );
    assert!(
        ordinary
            .gas()
            .primitive(ordinary.states()[0], || false)
            .unwrap()
            .pressure
            > 1.1
    );
}

#[test]
fn g3_homologous_compression_converges_to_independent_isentropic_euler_solution() {
    fn error(n: usize) -> f64 {
        let right_speed = -0.2;
        let mut field = flow(n, 0.0, 0.0, |x| right_speed * x);
        let end = 0.3;
        let mut steps = 0;
        while field.time().value() < end {
            field
                .step_with_moving_walls(
                    Time::new(end - field.time().value()),
                    StepControl::default(),
                    [0.0, right_speed],
                    || false,
                )
                .unwrap();
            steps += 1;
            assert!(steps < 2000);
        }
        // Exact smooth Euler solution: particle velocities are constant,
        // J=L/L0, rho=rho0/J and p=p0/J^gamma. No production wall/Riemann
        // solve is used to obtain this oracle (gamma=3.5/2.5 for the fixture).
        let jacobian = 1.0 + right_speed * end;
        let density = 1.0 / jacobian;
        let pressure = jacobian.powf(-1.4);
        let mut error = 0.0;
        for (i, &q) in field.states().iter().enumerate() {
            let p = field.gas().primitive(q, || false).unwrap();
            error += field.mesh().width(i)
                * ((p.density - density).abs() / density
                    + (p.pressure - pressure).abs() / pressure);
        }
        close(field.mesh().positions()[n], jacobian, 3e-14);
        error / jacobian
    }
    let errors = [error(16), error(32), error(64)];
    eprintln!("moving homogeneous compression 16/32/64 relative L1: {errors:?}");
    assert!(errors[1] < 0.8 * errors[0]);
    assert!(errors[2] < 0.8 * errors[1]);
}

#[test]
fn g4_moving_geometry_source_and_cancellation_refusals_are_atomic() {
    let fine = Mesh1d::uniform(1000, Length::new(0.001), 0.7).unwrap();
    assert!(
        MovingWorkspace::new(&fine).is_ok(),
        "retain declared uniform geometry"
    );
    let factory = || flow(12, 0.4, 0.0, |_| 0.0);
    let mut counted = factory();
    let mut polls = 0;
    counted
        .step_with_moving_walls(Time::new(0.003), StepControl::default(), [0.2, 0.0], || {
            polls += 1;
            false
        })
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut field = factory();
        let states = field.states().to_vec();
        let mesh = field.mesh().clone();
        let mut seen = 0;
        assert_eq!(
            field.step_with_moving_walls(
                Time::new(0.003),
                StepControl::default(),
                [0.2, 0.0],
                || {
                    seen += 1;
                    seen >= cutoff
                },
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(field.states(), states);
        assert_eq!(field.mesh(), &mesh);
        assert_eq!(field.time().value(), 0.0);
    }
    for velocity in [[f64::NAN, 0.0], [100.0, 0.0]] {
        let mut field = factory();
        let states = field.states().to_vec();
        let mesh = field.mesh().clone();
        assert!(
            field
                .step_with_moving_walls(
                    Time::new(0.003),
                    StepControl {
                        maximum_halvings: 0,
                        ..StepControl::default()
                    },
                    velocity,
                    || false,
                )
                .is_err()
        );
        assert_eq!(field.states(), states);
        assert_eq!(field.mesh(), &mesh);
        assert_eq!(field.time().value(), 0.0);
    }
    let mut open = factory();
    open.set_boundaries(Boundaries::Periodic, || false).unwrap();
    assert!(
        open.step_with_moving_walls(Time::new(0.003), StepControl::default(), [0.2, 0.0], || {
            false
        })
        .is_err()
    );
    let gas = gas(0.0, 0.0);
    let q = gas.conserved(1.0, 0.0, 1.0).unwrap();
    for mesh in [
        Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0, 1.1, 1.2], vec![0.5; 2]).unwrap(),
        Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0; 3], vec![0.4; 2]).unwrap(),
        Mesh1d::new(
            vec![1e16, 1e16 + 2.0, 1e16 + 4.0],
            vec![1.0; 3],
            vec![2.0; 2],
        )
        .unwrap(),
    ] {
        let original = mesh.clone();
        let mut field =
            Euler1d::new_on_mesh(gas.clone(), mesh, vec![q; 2], CLOSED, || false).unwrap();
        assert!(
            field
                .step_with_moving_walls(
                    Time::new(0.001),
                    StepControl {
                        maximum_halvings: 0,
                        ..StepControl::default()
                    },
                    [0.1; 2],
                    || false,
                )
                .is_err()
        );
        assert_eq!(field.mesh(), &original);
        assert_eq!(field.states(), &[q; 2]);
        assert_eq!(field.time().value(), 0.0);
    }
}
