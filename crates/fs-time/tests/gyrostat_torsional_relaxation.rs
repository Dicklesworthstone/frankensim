//! G0/G1/G3: causal shaft loss, retained memory, native recoil and positive heat.

#![allow(clippy::float_cmp)]

use fs_ga::{So3, Vec3};
use fs_time::lie::gyrostat::rotors::{
    AxialRotor, RotorBank, TorsionalExponentialError, TorsionalRelaxationBranch,
};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn zero() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}
fn bank() -> RotorBank<1> {
    RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap()],
    )
    .unwrap()
}
fn arm(k: f64, tau: f64) -> TorsionalRelaxationBranch<1> {
    TorsionalRelaxationBranch {
        projection: [1.0],
        stiffness_nm_per_rad: k,
        relaxation_time_s: tau,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}, allowance {tolerance:e}");
}

// Independent solution: x = q-z obeys x'' + x'/tau + (1/J+1/I)*k*x=0,
// s=x'+x/tau and integral(x)=(s0-s)/((1/J+1/I)*k).
fn maxwell(x0: f64, s0: f64, k: f64, tau: f64, t: f64) -> (f64, f64, f64) {
    let w2 = 1.5 * k;
    let a = 0.5 / tau;
    let d = a * a - w2;
    let derivative0 = s0 - x0 / tau;
    let (x, derivative) = if d.abs() < 1e-12 * w2 {
        let b = derivative0 + a * x0;
        let decay = (-a * t).exp();
        (decay * (x0 + b * t), decay * (b - a * (x0 + b * t)))
    } else if d < 0.0 {
        let w = (-d).sqrt();
        let b = (derivative0 + a * x0) / w;
        let (s, c) = (w * t).sin_cos();
        let decay = (-a * t).exp();
        (
            decay * (x0 * c + b * s),
            decay * ((b * w - a * x0) * c - (x0 * w + a * b) * s),
        )
    } else {
        let r1 = -w2 / (a + d.sqrt());
        let r2 = -a - d.sqrt();
        let c1 = (derivative0 - r2 * x0) / (r1 - r2);
        let c2 = x0 - c1;
        (
            c1 * (r1 * t).exp() + c2 * (r2 * t).exp(),
            r1 * c1 * (r1 * t).exp() + r2 * c2 * (r2 * t).exp(),
        )
    };
    let s = derivative + x / tau;
    (x, s, x - x0 + (s0 - s) / (tau * w2))
}

#[test]
fn one_maxwell_arm_matches_analytic_under_critical_overdamped_and_stiff_response() {
    let bank = bank();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.2), [0.8], [0.3])
        .unwrap();
    let (k, time, x0, s0) = (3.0_f64, 0.6, 0.1, 0.8);
    for tau in [1.0, 0.5 / (1.5 * k).sqrt(), 0.04, 0.001] {
        let mut stepper = bank
            .torsional_relaxation_stepper([0.0], [0.1], [0.0], [arm(k, tau)], time)
            .unwrap();
        let before = [0.1];
        let stored_before = stepper.internal_energy_j(&initial, &before).unwrap();
        let result = stepper.step(&initial, &before, [0.0], zero()).unwrap();
        let (x, s, dq) = maxwell(x0, s0, k, tau, time);
        close(result.motion.state_after.relative_angle_rad[0], 0.3 + dq, 1e-10);
        close(bank.relative_rates(&result.motion.state_after).unwrap()[0], s, 2e-10);
        close(result.viscous_angle_after_rad[0], 0.2 + dq - x, 1e-10);
        let stored_after = stepper
            .internal_energy_j(&result.motion.state_after, &result.viscous_angle_after_rad)
            .unwrap();
        close(result.motion.internal_energy_change_j, stored_after - stored_before, 2e-12);
        // Independent positive quadrature of k*x(t)^2/tau, not an energy gap.
        let intervals = 65_536;
        let mut heat = 0.0;
        for i in 0..=intervals {
            let x = maxwell(x0, s0, k, tau, time * i as f64 / intervals as f64).0;
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            heat += weight * k * x * x / tau;
        }
        heat *= time / (3.0 * intervals as f64);
        close(result.motion.dissipated_energy_j, heat, 1e-9);
        assert!(result.motion.dissipated_energy_j > 0.0);
        close(result.motion.work_energy_residual_j, 0.0, 5e-12);
        assert!(result.motion.angular_momentum_residual_world.norm() < 1e-12);
    }
}

#[test]
fn harmonic_response_refines_to_the_original_maxwell_complex_stiffness() {
    let bank = bank();
    let (equilibrium, k, tau, d, force) = (2.0, 3.0, 0.3, 0.15, 0.2);
    for omega in [0.5, 2.0, 8.0] {
        let a = omega * tau;
        let re = equilibrium - (2.0 / 3.0) * omega * omega + k * a * a / (1.0 + a * a);
        let im = omega * d + k * a / (1.0 + a * a);
        let den = re * re + im * im;
        let (qr, qi) = (force * re / den, -force * im / den);
        let zr = (qr + a * qi) / (1.0 + a * a);
        let speed0 = -omega * qi;
        let period = core::f64::consts::TAU / omega;
        let exact_heat = 0.5 * period * (qr * qr + qi * qi)
            * (d * omega * omega + (k / tau) * a * a / (1.0 + a * a));
        let mut errors = Vec::new();
        for count in [128, 256, 512] {
            let h = period / count as f64;
            let mut stepper = bank
                .torsional_relaxation_stepper([equilibrium], [0.0], [d], [arm(k, tau)], h)
                .unwrap();
            let mut state = bank
                .state_from_velocities(
                    So3::identity(), Vec3::new(0.0, 0.0, -speed0 / 3.0), [speed0], [qr],
                )
                .unwrap();
            let mut memory = [zr];
            let mut heat = 0.0;
            for i in 0..count {
                let load = force * (omega * h * (i as f64 + 0.5)).cos();
                let result = stepper.step(&state, &memory, [load], zero()).unwrap();
                close(result.motion.work_energy_residual_j, 0.0, 5e-12);
                heat += result.motion.dissipated_energy_j;
                state = result.motion.state_after;
                memory = result.viscous_angle_after_rad;
            }
            let error = (state.relative_angle_rad[0] - qr)
                .abs()
                .max((bank.relative_rates(&state).unwrap()[0] - speed0).abs() / omega)
                .max((memory[0] - zr).abs());
            errors.push(error);
            if count == 512 {
                close(heat / exact_heat, 1.0, 8e-5);
            }
        }
        assert!(
            errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
            "omega={omega}: held-input refinement {errors:?}"
        );
    }
}

#[test]
fn common_rotation_keeps_relaxed_memory_and_zero_physical_heat() {
    let bank = bank();
    // k=10 and z=0.1 exercise NONZERO normalized memory with a nonbinary
    // sqrt(k). Reconstructing its absolute physical value by divide/multiply
    // can change z by one ULP even when the material state has not evolved.
    // Equilibrium stiffness is zero, so force-relative translation cannot
    // collapse the normalized history to zero and mask this regression.
    for (k, strain) in [(4.0, 0.25), (10.0, 0.1)] {
        let state = bank
            .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.5), [0.0], [strain])
            .unwrap();
        let mut stepper = bank
            .torsional_relaxation_stepper([0.0], [0.0], [0.0], [arm(k, 0.2)], 0.1)
            .unwrap();
        let memory = stepper.relaxed_memory(&state).unwrap();
        assert_eq!(memory, [strain]);
        assert_eq!(stepper.internal_energy_j(&state, &memory).unwrap(), 0.0);
        let result = stepper.step(&state, &memory, [0.0], zero()).unwrap();
        assert_eq!(result.motion.dissipated_energy_j, 0.0);
        assert_eq!(result.motion.internal_energy_change_j, 0.0);
        assert_eq!(result.viscous_angle_after_rad, memory);
        assert_eq!(result.motion.state_after.relative_angle_rad, state.relative_angle_rad);
        assert_eq!(result.motion.state_after.axial_momentum_nms, state.axial_momentum_nms);
    }
}

#[test]
fn coupled_material_ports_keep_both_shaft_reactions_and_all_history() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(3.0, 3.0, 2.0), So3::identity()).unwrap(),
        [
            AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap(),
            AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 0.5).unwrap(),
        ],
    )
    .unwrap();
    let branches = [
        TorsionalRelaxationBranch {
            projection: [1.0, -2.0],
            stiffness_nm_per_rad: 2.0,
            relaxation_time_s: 0.3,
        },
        TorsionalRelaxationBranch {
            projection: [0.5, 1.0],
            stiffness_nm_per_rad: 0.7,
            relaxation_time_s: 2.0,
        },
    ];
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.1), [0.5, -0.4], [0.3, -0.4])
        .unwrap();
    let mut stepper = bank
        .torsional_relaxation_stepper(
            [1.0, 0.5], [0.1, -0.2], [0.2, 0.1], branches, 0.6,
        )
        .unwrap();
    let result = stepper.step(&initial, &[0.1, -0.2], [0.3, -0.1], zero()).unwrap();
    // Independent original-coordinate RK4: q, p, the two physical viscous
    // strains, positive heat and signed applied work. L remains fixed.
    let total = initial.total_momentum_body.z;
    let derivative = |x: [f64; 8]| {
        let omega = (total - x[2] - x[3]) / 2.0;
        let s = [x[2] - omega, 2.0 * x[3] - omega];
        let e = [x[0] - 2.0 * x[1] - x[4], 0.5 * x[0] + x[1] - x[5]];
        [
            s[0], s[1],
            0.3 - x[0] - 0.2 * s[0] - 2.0 * e[0] - 0.35 * e[1],
            -0.1 - 0.5 * x[1] - 0.1 * s[1] + 4.0 * e[0] - 0.7 * e[1],
            e[0] / 0.3, e[1] / 2.0,
            0.2 * s[0] * s[0] + 0.1 * s[1] * s[1] + 2.0 * e[0] * e[0] / 0.3
                + 0.7 * e[1] * e[1] / 2.0,
            0.3 * s[0] - 0.1 * s[1],
        ]
    };
    let mut x = [
        0.2, -0.2, initial.axial_momentum_nms[0], initial.axial_momentum_nms[1],
        0.1, -0.2, 0.0, 0.0,
    ];
    let h = 0.6 / 4096.0;
    for _ in 0..4096 {
        let a = derivative(x);
        let b = derivative(core::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = derivative(core::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = derivative(core::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..8 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    for i in 0..2 {
        close(result.motion.state_after.relative_angle_rad[i], x[i] + [0.1, -0.2][i], 2e-11);
        close(result.motion.state_after.axial_momentum_nms[i], x[2 + i], 2e-11);
        close(result.viscous_angle_after_rad[i], x[4 + i], 2e-11);
    }
    close(result.motion.dissipated_energy_j, x[6], 2e-11);
    close(result.motion.drive_work_j.iter().sum(), x[7], 2e-11);
    close(result.motion.work_energy_residual_j, 0.0, 3e-12);
    assert!(result.motion.angular_momentum_residual_world.norm() < 1e-12);
}

#[test]
fn no_memory_delegates_exactly_and_cancellation_retries_the_same_history() {
    let bank = bank();
    let state = bank
        .state_from_velocities(So3::identity(), Vec3::new(0.0, 0.0, 0.2), [0.8], [0.1])
        .unwrap();
    let mut ordinary = bank
        .torsional_exponential_stepper([3.0], [0.2], [0.4], 0.1)
        .unwrap();
    let expected = ordinary.step(&state, [0.3], zero()).unwrap();
    let mut empty = bank
        .torsional_relaxation_stepper::<0>([3.0], [0.2], [0.4], [], 0.1)
        .unwrap();
    assert_eq!(empty.step(&state, &[], [0.3], zero()).unwrap().motion, expected);
    let mut stepper = bank
        .torsional_relaxation_stepper([3.0], [0.2], [0.4], [arm(2.0, 0.2)], 0.1)
        .unwrap();
    let expected = stepper.step(&state, &[-0.05], [0.3], zero()).unwrap();
    let mut polls = 0;
    stepper
        .step_with_cancel(&state, &[-0.05], [0.3], zero(), || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=polls {
        let mut seen = 0;
        assert!(matches!(
            stepper.step_with_cancel(&state, &[-0.05], [0.3], zero(), || {
                seen += 1;
                seen == stop
            }),
            Err(TorsionalExponentialError::Cancelled)
        ));
        assert_eq!(stepper.step(&state, &[-0.05], [0.3], zero()).unwrap(), expected);
    }
    for bad in [
        arm(-1.0, 0.2), arm(1.0, 0.0), arm(1.0, f64::INFINITY),
        TorsionalRelaxationBranch {
            projection: [f64::NAN],
            ..arm(1.0, 0.2)
        },
    ] {
        assert!(bank.torsional_relaxation_stepper([0.0], [0.0], [0.0], [bad], 0.1).is_err());
    }
    assert!(stepper.step(&state, &[f64::NAN], [0.3], zero()).is_err());
    assert_eq!(stepper.step(&state, &[-0.05], [0.3], zero()).unwrap(), expected);
}
