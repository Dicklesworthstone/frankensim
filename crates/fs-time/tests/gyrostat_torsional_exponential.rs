//! Full original spring/drag dynamics, independent heat and native split refinement.

// Exact zero loss, unchanged states and bit-for-bit delegation are requirements.
#![allow(clippy::float_cmp)]

use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_la::factor::lu;
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, MAX_EXPONENTIAL_TORSIONAL_ROTORS, RotorBank, RotorState, TorsionalExponentialError,
};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tol={tolerance:e}"
    );
}
fn zero() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn coaxial<const N: usize>() -> RotorBank<N> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap(); N],
    )
    .unwrap()
}
fn oblique() -> RotorBank<3> {
    RotorBank::from_reduced(
        InertiaTensor3::new(Mat3 {
            m: [2.4, 0.3, -0.2, 0.3, 3.1, 0.25, -0.2, 0.25, 4.0],
        })
        .unwrap(),
        [
            AxialRotor::new(Vec3::new(1.0, 0.2, 0.0), 0.4).unwrap(),
            AxialRotor::new(Vec3::new(0.3, 1.0, 0.4), 0.7).unwrap(),
            AxialRotor::new(Vec3::new(-0.2, 0.1, 1.0), 0.9).unwrap(),
        ],
    )
    .unwrap()
}

// Independent original-coordinate harmonic solution, including critical and
// overdamped limits. This never calls the cached exponential/Gramian owner.
fn oscillator(q0: f64, s0: f64, k: f64, d: f64, drive: f64, rest: f64, t: f64) -> (f64, f64) {
    let equilibrium = rest + drive / k;
    let a = q0 - equilibrium;
    let beta = 0.75 * d;
    let squared = 1.5 * k;
    let discriminant = beta * beta - squared;
    if discriminant.abs() < 1e-13 * squared {
        let b = s0 + beta * a;
        let decay = (-beta * t).exp();
        (
            equilibrium + decay * (a + b * t),
            decay * (b - beta * (a + b * t)),
        )
    } else if discriminant < 0.0 {
        let w = (-discriminant).sqrt();
        let b = (s0 + beta * a) / w;
        let c = (w * t).cos();
        let s = (w * t).sin();
        let decay = (-beta * t).exp();
        (
            equilibrium + decay * (a * c + b * s),
            decay * ((-beta * a + w * b) * c + (-beta * b - w * a) * s),
        )
    } else {
        let root = discriminant.sqrt();
        let r1 = -squared / (beta + root);
        let r2 = -beta - root;
        let c1 = (s0 - r2 * a) / (r1 - r2);
        let c2 = a - c1;
        let e1 = (r1 * t).exp();
        let e2 = (r2 * t).exp();
        (equilibrium + c1 * e1 + c2 * e2, r1 * c1 * e1 + r2 * c2 * e2)
    }
}

#[test]
fn under_critical_overdamped_and_stiff_losses_match_analytic_motion_and_independent_heat() {
    let bank = coaxial::<1>();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.2), [0.8], [0.1])
        .unwrap();
    let (k, drive, rest, time) = (3.0, 0.3, 0.2, 0.6);
    for damping in [0.4, 2.0 * (k / 1.5_f64).sqrt(), 8.0, 400.0] {
        let mut prepared = bank
            .torsional_exponential_stepper([k], [rest], [damping], time)
            .unwrap();
        let step = prepared
            .step(&initial, [drive], zero())
            .unwrap_or_else(|error| panic!("d={damping:e}: {error:?}"));
        let (q, s) = oscillator(0.1, 0.8, k, damping, drive, rest, time);
        close(step.state_after.relative_angle_rad[0], q, 3e-12);
        close(bank.relative_rates(&step.state_after).unwrap()[0], s, 2e-11);
        close(step.drive_work_j[0], drive * (q - 0.1), 1e-12);
        // Positive independent Simpson integral of the original dashpot power.
        // The largest fast relaxation has > 700 quadrature nodes per e-fold.
        let intervals = 262_144_usize;
        let mut integral = 0.0;
        for i in 0..=intervals {
            let t = time * i as f64 / intervals as f64;
            let speed = oscillator(0.1, 0.8, k, damping, drive, rest, t).1;
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * damping * speed * speed;
        }
        integral *= time / (3.0 * intervals as f64);
        close(step.dissipated_energy_j, integral, 2e-11);
        assert!(step.dissipated_energy_j > 0.0);
        close(step.work_energy_residual_j, 0.0, 3e-12);
        close(
            step.elastic_energy_change_j[0],
            0.5 * k * ((q - rest).powi(2) - (0.1 - rest).powi(2)),
            2e-12,
        );
        if damping > 100.0 {
            assert!(
                bank.torsional_stepper([k], [rest], [damping], time)
                    .is_err()
            );
        }
    }
}

#[test]
fn unresolved_midpoint_oscillations_are_retained_and_zero_heat_stays_zero() {
    let bank = coaxial::<1>();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.2), [0.7], [0.3])
        .unwrap();
    let (k, rest, time) = (800.0, 0.1, 1.2);
    assert!(bank.torsional_stepper([k], [rest], [0.0], time).is_err());
    let mut prepared = bank
        .torsional_exponential_stepper([k], [rest], [0.0], time)
        .unwrap();
    let step = prepared.step(&initial, [0.0], zero()).unwrap();
    let (q, s) = oscillator(0.3, 0.7, k, 0.0, 0.0, rest, time);
    close(step.state_after.relative_angle_rad[0], q, 3e-12);
    close(bank.relative_rates(&step.state_after).unwrap()[0], s, 8e-11);
    assert_eq!(step.dissipated_energy_j, 0.0);
    close(
        step.kinetic_energy_change_j + step.elastic_energy_change_j[0],
        0.0,
        4e-11,
    );
    assert!(step.angular_momentum_residual_world.norm() < 2e-13);
}

fn packed(state: &RotorState<3>) -> [f64; 14] {
    [
        state.total_momentum_body.x,
        state.total_momentum_body.y,
        state.total_momentum_body.z,
        state.axial_momentum_nms[0],
        state.axial_momentum_nms[1],
        state.axial_momentum_nms[2],
        state.relative_angle_rad[0],
        state.relative_angle_rad[1],
        state.relative_angle_rad[2],
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}

// Independent dense original-frame Euler equations and positive power integrals.
// Entries 9..14 are physical heat, external work, then the three drive works.
#[allow(clippy::too_many_arguments)] // Keep the original independent physical inputs explicit.
fn reference(
    bank: &RotorBank<3>,
    state: &RotorState<3>,
    stiffness: [f64; 3],
    rest: [f64; 3],
    damping: [f64; 3],
    drive: [f64; 3],
    external: Vec3,
    time: f64,
    count: usize,
) -> [f64; 14] {
    let factor = lu(&bank.reduced_inertia().matrix().m, 3).unwrap();
    let derivative = |y: [f64; 14]| {
        let total = Vec3::new(y[0], y[1], y[2]);
        let mut carrier = total;
        for i in 0..3 {
            carrier = carrier - bank.rotors()[i].axis_body().scale(y[3 + i]);
        }
        let mut omega = [carrier.x, carrier.y, carrier.z];
        factor.solve(&mut omega);
        let omega = Vec3::new(omega[0], omega[1], omega[2]);
        let rate = total.cross(omega) + external;
        let mut out = [0.0; 14];
        out[..3].copy_from_slice(&[rate.x, rate.y, rate.z]);
        for i in 0..3 {
            out[6 + i] = y[3 + i] / bank.rotors()[i].axial_inertia_kg_m2()
                - bank.rotors()[i].axis_body().dot(omega);
            out[3 + i] = drive[i] - stiffness[i] * (y[6 + i] - rest[i]) - damping[i] * out[6 + i];
            out[9] += damping[i] * out[6 + i] * out[6 + i];
            out[11 + i] = drive[i] * out[6 + i];
        }
        out[10] = external.dot(omega);
        out
    };
    let dt = time / count as f64;
    let mut y = packed(state);
    for _ in 0..count {
        let k1 = derivative(y);
        let k2 = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]));
        let k3 = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]));
        let k4 = derivative(core::array::from_fn(|i| y[i] + dt * k3[i]));
        y = core::array::from_fn(|i| y[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
    }
    y
}

#[test]
fn full_coupled_stiff_stages_keep_zero_spring_and_undamped_channels() {
    let bank = oblique();
    // Total L=0 makes body L remain exactly zero during native drift. All
    // coupled p/q dynamics are then the SAME linear shaft flow across stages.
    let initial = RotorState {
        rotation: rotation(Vec3::new(0.2, -0.1, 0.3)),
        total_momentum_body: Vec3::default(),
        axial_momentum_nms: [0.6, -0.4, 0.3],
        relative_angle_rad: [0.1, -0.2, 0.3],
    };
    let stiffness = [2.0, 0.0, 3.0];
    let rest = [0.0, 1e200, 0.1];
    let damping = [300.0, 0.0, 0.2];
    let drive = [0.2, 0.1, -0.1];
    let time = 0.2;
    let oracle = reference(
        &bank,
        &initial,
        stiffness,
        rest,
        damping,
        drive,
        Vec3::default(),
        time,
        32768,
    );
    let coarse = reference(
        &bank,
        &initial,
        stiffness,
        rest,
        damping,
        drive,
        Vec3::default(),
        time,
        16384,
    );
    assert!(
        oracle
            .into_iter()
            .zip(coarse)
            .all(|(a, b)| (a - b).abs() < 8e-10)
    );
    let mut prepared = bank
        .torsional_exponential_stepper(stiffness, rest, damping, time)
        .unwrap();
    let step = prepared.step(&initial, drive, zero()).unwrap();
    for (a, b) in packed(&step.state_after).into_iter().take(9).zip(oracle) {
        close(a, b, 8e-11);
    }
    close(step.dissipated_energy_j, oracle[9], 8e-10);
    for i in 0..3 {
        close(step.drive_work_j[i], oracle[11 + i], 2e-11);
    }
    assert_eq!(step.elastic_energy_change_j[1], 0.0);
    close(step.work_energy_residual_j, 0.0, 4e-11);
    assert!(
        bank.torsional_stepper(stiffness, rest, damping, time)
            .is_err()
    );
    let mut half = bank
        .torsional_exponential_stepper(stiffness, rest, damping, 0.5 * time)
        .unwrap();
    let first = half.step(&initial, drive, zero()).unwrap();
    let second = half.step(&first.state_after, drive, zero()).unwrap();
    for i in 0..3 {
        close(
            second.state_after.axial_momentum_nms[i],
            step.state_after.axial_momentum_nms[i],
            2e-12,
        );
        close(
            second.state_after.relative_angle_rad[i],
            step.state_after.relative_angle_rad[i],
            2e-12,
        );
    }
    close(
        first.dissipated_energy_j + second.dissipated_energy_j,
        step.dissipated_energy_j,
        5e-12,
    );
}

#[test]
fn oblique_gyrostat_trajectory_heat_and_external_work_refine_at_second_order() {
    let bank = oblique();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.2, 0.3, 0.1)),
            Vec3::new(0.6, -0.3, 0.8),
            [1.0, -0.5, 0.7],
            [0.2, -0.1, 0.3],
        )
        .unwrap();
    let stiffness = [0.7, 1.1, 0.5];
    let damping = [0.3, 2.0, 5.0];
    let rest = [0.1, 0.2, -0.1];
    let drive = [0.1, -0.2, 0.4];
    let external = Vec3::new(0.1, 0.2, -0.1);
    let time = 0.7;
    let oracle = reference(
        &bank, &initial, stiffness, rest, damping, drive, external, time, 8192,
    );
    let coarse = reference(
        &bank, &initial, stiffness, rest, damping, drive, external, time, 4096,
    );
    assert!(
        oracle
            .into_iter()
            .zip(coarse)
            .all(|(a, b)| (a - b).abs() < 1e-11)
    );
    let mut previous = None;
    for count in [8_usize, 16, 32] {
        let mut prepared = bank
            .torsional_exponential_stepper(stiffness, rest, damping, time / count as f64)
            .unwrap();
        let mut state = initial;
        let (mut heat, mut external_work, mut drive_work) = (0.0, 0.0, [0.0; 3]);
        for _ in 0..count {
            let step = prepared
                .step(&state, drive, HeldTorque::Body(external))
                .unwrap();
            close(step.work_energy_residual_j, 0.0, 2e-12);
            assert!(step.dissipated_energy_j >= 0.0);
            assert!(step.angular_momentum_residual_world.norm() < 2e-13);
            heat += step.dissipated_energy_j;
            external_work += step.external_work_j;
            for (total, work) in drive_work.iter_mut().zip(step.drive_work_j) {
                *total += work;
            }
            state = step.state_after;
        }
        let mut actual = packed(&state);
        actual[9] = heat;
        actual[10] = external_work;
        actual[11..].copy_from_slice(&drive_work);
        let error = actual
            .into_iter()
            .zip(oracle)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        eprintln!("coupled torsional exponential: count={count}, ODE/work/heat error={error:e}");
        if let Some(old) = previous {
            assert!(old / error > 3.7 && old / error < 4.3, "{old:e}/{error:e}");
        }
        previous = Some(error);
    }
}

#[test]
fn exact_zero_stiffness_delegation_and_world_frame_covariance() {
    let bank = oblique();
    let state = bank
        .state_from_velocities(
            rotation(Vec3::new(0.1, 0.2, -0.3)),
            Vec3::new(0.2, 0.1, -0.3),
            [0.8, -0.2, 0.4],
            [0.1, -0.2, 0.3],
        )
        .unwrap();
    let drive = [0.2, -0.3, 0.1];
    let torque = Vec3::new(0.1, -0.2, 0.3);
    for damping in [[0.0; 3], [2.0, 0.0, 0.4]] {
        let expected = bank
            .viscous_exponential_stepper(damping, 0.05)
            .unwrap()
            .step(&state, drive, HeldTorque::Space(torque))
            .unwrap();
        let actual = bank
            .torsional_exponential_stepper([0.0; 3], [123.0; 3], damping, 0.05)
            .unwrap()
            .step(&state, drive, HeldTorque::Space(torque))
            .unwrap();
        assert_eq!(actual.state_after, expected.state_after);
        assert_eq!(actual.drive_work_j, expected.drive_work_j);
        assert_eq!(actual.dissipated_energy_j, expected.dissipated_energy_j);
        assert_eq!(actual.net_shaft_work_j, expected.net_shaft_work_j);
        assert_eq!(
            actual.work_energy_residual_j,
            expected.work_energy_residual_j
        );
        assert_eq!(actual.elastic_energy_change_j, [0.0; 3]);
    }
    let mut prepared = bank
        .torsional_exponential_stepper([1.0, 2.0, 3.0], [0.0; 3], [0.2, 0.3, 0.4], 0.1)
        .unwrap();
    let original = prepared
        .step(&state, drive, HeldTorque::Space(torque))
        .unwrap();
    let transform = rotation(Vec3::new(-0.5, 0.2, 0.3));
    let transformed = RotorState {
        rotation: transform.compose(state.rotation).unwrap(),
        ..state
    };
    let actual = prepared
        .step(
            &transformed,
            drive,
            HeldTorque::Space(transform.rotate(torque).unwrap()),
        )
        .unwrap();
    for i in 0..3 {
        close(
            actual.state_after.axial_momentum_nms[i],
            original.state_after.axial_momentum_nms[i],
            2e-14,
        );
        close(
            actual.state_after.relative_angle_rad[i],
            original.state_after.relative_angle_rad[i],
            2e-14,
        );
    }
    let wanted = transform.compose(original.state_after.rotation).unwrap();
    assert!(
        wanted
            .inverse()
            .compose(actual.state_after.rotation)
            .unwrap()
            .log()
            .angular
            .norm()
            < 2e-13
    );
    close(
        actual.dissipated_energy_j,
        original.dissipated_energy_j,
        2e-14,
    );
    close(actual.external_work_j, original.external_work_j, 2e-14);
}

#[test]
fn cancellation_capacity_bad_sources_and_unresolved_phases_refuse_atomically() {
    const OVER: usize = MAX_EXPONENTIAL_TORSIONAL_ROTORS + 1;
    let bank = coaxial::<1>();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [0.8], [0.1])
        .unwrap();
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            bank.torsional_exponential_stepper([bad], [0.0], [0.1], 0.1)
                .is_err()
        );
        assert!(
            bank.torsional_exponential_stepper([1.0], [0.0], [bad], 0.1)
                .is_err()
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            bank.torsional_exponential_stepper([1.0], [0.0], [0.1], bad)
                .is_err()
        );
    }
    assert!(matches!(
        bank.torsional_exponential_stepper_with_cancel([1.0], [0.0], [0.1], 0.1, || true),
        Err(TorsionalExponentialError::Cancelled)
    ));
    let mut prepared = bank
        .torsional_exponential_stepper([1.0], [0.0], [0.1], 0.1)
        .unwrap();
    let expected = prepared.step(&initial, [0.2], zero()).unwrap();
    for cancel_at in 1..=4 {
        let mut polls = 0;
        assert!(matches!(
            prepared.step_with_cancel(&initial, [0.2], zero(), || {
                polls += 1;
                polls == cancel_at
            }),
            Err(TorsionalExponentialError::Cancelled)
        ));
        assert_eq!(prepared.step(&initial, [0.2], zero()).unwrap(), expected);
    }
    assert!(prepared.step(&initial, [f64::NAN], zero()).is_err());
    assert_eq!(prepared.step(&initial, [0.2], zero()).unwrap(), expected);
    let coarse = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [1.0], [1e16])
        .unwrap();
    assert!(
        bank.torsional_exponential_stepper([1.0], [1e16], [0.1], 0.1)
            .unwrap()
            .step(&coarse, [0.0], zero())
            .is_err()
    );
    assert_eq!(coarse.relative_angle_rad[0], 1e16);
    let large = coaxial::<OVER>();
    assert!(
        large
            .torsional_exponential_stepper([1.0; OVER], [0.0; OVER], [0.0; OVER], 0.1)
            .is_err()
    );
    assert!(
        large
            .torsional_exponential_stepper([0.0; OVER], [0.0; OVER], [0.0; OVER], 0.1)
            .is_ok()
    );
}

#[test]
fn force_supported_equilibrium_has_no_motion_or_fictitious_heat() {
    let bank = coaxial::<1>();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [0.0], [0.375])
        .unwrap();
    let mut prepared = bank
        .torsional_exponential_stepper([8.0], [0.25], [1e4], 1.0)
        .unwrap();
    let step = prepared.step(&initial, [1.0], zero()).unwrap();
    assert_eq!(step.state_after, initial);
    assert_eq!(step.dissipated_energy_j, 0.0);
    assert_eq!(step.drive_work_j, [0.0]);
    assert_eq!(step.elastic_energy_change_j, [0.0]);
    assert_eq!(step.work_energy_residual_j, 0.0);
}
