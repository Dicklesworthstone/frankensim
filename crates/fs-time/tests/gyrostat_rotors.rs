//! Physical rotor inertia, finite shaft work, native precession and relative phase.
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_time::lie::gyrostat::rotors::{AxialRotor, MAX_AXIAL_ROTORS, RotorBank, RotorState};
use fs_time::lie::gyrostat::{HeldTorque, InertiaTensor3};

fn rotation(v: Vec3) -> So3 {
    So3::exp(So3Tangent::new(v)).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, error {}",
        (a - b).abs()
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!((a - b).norm() <= tolerance, "{a:?} != {b:?}");
}
fn attitude(a: So3, b: So3, tolerance: f64) {
    assert!(
        a.body_minus(b).unwrap().angular.norm() <= tolerance,
        "{a:?} != {b:?}"
    );
}
fn zero_torque() -> HeldTorque {
    HeldTorque::Body(Vec3::default())
}

fn locked<const N: usize>(inertia: Mat3, rotors: &[AxialRotor; N]) -> Mat3 {
    let mut m = inertia.m;
    for rotor in rotors {
        let axis = rotor.axis_body();
        let a = [axis.x, axis.y, axis.z];
        for i in 0..3 {
            for j in 0..3 {
                m[3 * i + j] += rotor.axial_inertia_kg_m2() * a[i] * a[j];
            }
        }
    }
    Mat3 { m }
}

// Independent complete locked-tensor kinetic expression, before eliminating
// rotor rates: .5 omega.I_locked.omega + sum(J*s*a.omega + .5 J*s*s).
fn locked_energy<const N: usize>(bank: &RotorBank<N>, state: &RotorState<N>) -> f64 {
    let omega = bank.angular_velocity(state).unwrap();
    let rate = bank.relative_rates(state).unwrap();
    let full = locked(bank.reduced_inertia().matrix(), bank.rotors());
    let mut energy = 0.5 * omega.dot(full.apply(omega));
    for (rotor, s) in bank.rotors().iter().zip(rate) {
        energy += rotor.axial_inertia_kg_m2() * (s * rotor.axis_body().dot(omega) + 0.5 * s * s);
    }
    energy
}

#[test]
fn finite_multi_axis_impulse_has_full_locked_energy_and_coupled_carrier_recoil() {
    let basis = rotation(Vec3::new(0.4, -0.3, 0.2));
    let reduced = InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), basis).unwrap();
    let rotors = [
        AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 0.6).unwrap(),
        AxialRotor::new(Vec3::new(1.0, 1.0, 0.0), 0.8).unwrap(),
        AxialRotor::new(Vec3::new(0.2, -0.3, 1.0), 0.4).unwrap(),
    ];
    let bank = RotorBank::from_locked(locked(reduced.matrix(), &rotors), rotors).unwrap();
    for (a, b) in bank
        .reduced_inertia()
        .matrix()
        .m
        .iter()
        .zip(reduced.matrix().m)
    {
        close(*a, b, 5e-15);
    }
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.2, 0.5, 0.3)),
            Vec3::new(0.8, -0.5, 0.6),
            [2.0, -1.0, 0.4],
            [0.1, -0.3, 0.7],
        )
        .unwrap();
    let impulses = [0.3, -0.2, 0.1];
    let mut stored_delta = Vec3::default();
    for (rotor, impulse) in rotors.iter().zip(impulses) {
        stored_delta = stored_delta + rotor.axis_body().scale(impulse);
    }
    for external in [Vec3::default(), Vec3::new(0.15, 0.05, -0.07)] {
        let kick = bank.apply_impulse(&initial, impulses, external).unwrap();
        let omega0 = bank.angular_velocity(&initial).unwrap();
        let omega1 = bank.angular_velocity(&kick.state_after).unwrap();
        let delta_omega = reduced.velocity(external - stored_delta).unwrap();
        vector(omega1, omega0 + delta_omega, 3e-14);
        vector(
            kick.state_after.total_momentum_body,
            initial.total_momentum_body + external,
            1e-15,
        );
        assert_eq!(kick.state_after.rotation, initial.rotation);
        assert_eq!(
            kick.state_after.relative_angle_rad,
            initial.relative_angle_rad
        );
        let rates0 = bank.relative_rates(&initial).unwrap();
        let rates1 = bank.relative_rates(&kick.state_after).unwrap();
        for i in 0..3 {
            close(
                rates1[i] - rates0[i],
                impulses[i] / rotors[i].axial_inertia_kg_m2()
                    - rotors[i].axis_body().dot(delta_omega),
                3e-14,
            );
            close(
                kick.shaft_work_j[i],
                impulses[i] * 0.5 * (rates0[i] + rates1[i]),
                3e-14,
            );
        }
        let before = locked_energy(&bank, &initial);
        let after = locked_energy(&bank, &kick.state_after);
        close(bank.kinetic_energy(&initial).unwrap(), before, 3e-14);
        close(
            bank.kinetic_energy(&kick.state_after).unwrap(),
            after,
            3e-14,
        );
        close(kick.kinetic_energy_change_j, after - before, 3e-14);
        close(
            kick.external_work_j,
            external.dot((omega0 + omega1).scale(0.5)),
            3e-14,
        );
        close(kick.work_energy_residual_j, 0.0, 3e-14);
    }
}

#[test]
fn coaxial_spin_up_and_braking_have_exact_relative_phase_work_and_counterrotation() {
    let reduced =
        InertiaTensor3::from_principal(Vec3::new(2.4, 2.4, 2.0), So3::identity()).unwrap();
    let rotor = AxialRotor::new(Vec3::new(0.0, 0.0, 1.0), 1.0).unwrap();
    let bank = RotorBank::from_locked(locked(reduced.matrix(), &[rotor]), [rotor]).unwrap();
    let r0 = rotation(Vec3::new(0.2, -0.4, 0.1));
    let initial = bank
        .state_from_velocities(r0, Vec3::new(0.0, 0.0, 0.3), [1.2], [0.4])
        .unwrap();
    let total_world = r0.rotate(initial.total_momentum_body).unwrap();
    let mut state = initial;
    let mut work = 0.0;
    for stage in 0..2 {
        let torque = if stage == 0 { 0.8 } else { -0.8 };
        let start = state;
        let start_rate = bank.relative_rates(&start).unwrap()[0];
        let start_omega = bank.angular_velocity(&start).unwrap().z;
        for n in 1..=20 {
            let previous = state;
            let step = bank.step(&state, [torque], zero_torque(), 0.05).unwrap();
            state = step.state_after;
            let time = 0.05 * f64::from(n);
            let omega = start_omega - torque * time / 2.0;
            let relative = start_rate + torque * time * 1.5;
            close(
                state.axial_momentum_nms[0],
                start.axial_momentum_nms[0] + torque * time,
                4e-14,
            );
            close(bank.angular_velocity(&state).unwrap().z, omega, 4e-14);
            close(bank.relative_rates(&state).unwrap()[0], relative, 4e-14);
            close(
                state.relative_angle_rad[0],
                start.relative_angle_rad[0] + start_rate * time + 0.75 * torque * time * time,
                5e-14,
            );
            attitude(
                state.rotation,
                start
                    .rotation
                    .compose(rotation(Vec3::new(
                        0.0,
                        0.0,
                        start_omega * time - 0.25 * torque * time * time,
                    )))
                    .unwrap(),
                5e-14,
            );
            vector(
                state.rotation.rotate(state.total_momentum_body).unwrap(),
                total_world,
                5e-14,
            );
            close(
                step.shaft_work_j[0],
                torque * (state.relative_angle_rad[0] - previous.relative_angle_rad[0]),
                5e-14,
            );
            close(step.drift_energy_defect_j, 0.0, 5e-14);
            close(step.work_energy_residual_j, 0.0, 5e-14);
            work += step.shaft_work_j[0];
        }
    }
    close(state.relative_angle_rad[0], 4.0, 1e-13);
    close(
        bank.kinetic_energy(&state).unwrap(),
        bank.kinetic_energy(&initial).unwrap(),
        1e-13,
    );
    close(work, 0.0, 1e-13);
}

#[test]
fn finite_spin_change_matches_closed_three_dimensional_gyrostatic_precession() {
    let a = 2.0;
    let c = 3.5;
    let j = 0.7;
    let tau = 0.4;
    let total_time = 1.2;
    let basis = rotation(Vec3::new(0.3, -0.5, 0.2));
    let axis = basis.rotate(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(a, a, c), basis).unwrap(),
        [AxialRotor::new(axis, j).unwrap()],
    )
    .unwrap();
    let w0 = Vec3::new(0.6, -0.4, 0.3);
    let r0 = rotation(Vec3::new(-0.2, 0.1, 0.4));
    let initial = bank
        .state_from_velocities(r0, basis.rotate(w0).unwrap(), [2.0], [0.25])
        .unwrap();
    let p0 = initial.axial_momentum_nms[0];
    let l_world = r0.rotate(initial.total_momentum_body).unwrap();
    let phase = ((c - a) * w0.z + p0) / a * total_time + tau * total_time * total_time / (2.0 * c);
    let (sin, cos) = phase.sin_cos();
    let exact_omega = basis
        .rotate(Vec3::new(
            cos * w0.x - sin * w0.y,
            sin * w0.x + cos * w0.y,
            w0.z - tau * total_time / c,
        ))
        .unwrap();
    let exact_rotation = rotation(l_world.scale(total_time / a))
        .compose(r0)
        .unwrap()
        .compose(rotation(axis.scale(-phase)))
        .unwrap();
    let exact_angle = 0.25
        + (p0 / j - w0.z) * total_time
        + 0.5 * tau * (1.0 / j + 1.0 / c) * total_time * total_time;
    let mut previous_error = None;
    for count in [24, 48, 96] {
        let mut state = initial;
        let h = total_time / f64::from(count);
        for _ in 0..count {
            let step = bank.step(&state, [tau], zero_torque(), h).unwrap();
            vector(step.angular_momentum_residual_world, Vec3::default(), 2e-13);
            close(step.work_energy_residual_j, 0.0, 2e-13);
            state = step.state_after;
        }
        let error = (bank.angular_velocity(&state).unwrap() - exact_omega).norm()
            + state
                .rotation
                .body_minus(exact_rotation)
                .unwrap()
                .angular
                .norm()
            + (state.relative_angle_rad[0] - exact_angle).abs();
        // The aligned native axisymmetric fast path can make this exact to
        // roundoff; a transformed stored vector may instead use native splitting.
        if let Some(old) = previous_error {
            if old > 1e-9 && error > 1e-9 {
                assert!(old / error > 3.7 && old / error < 4.3, "{old} / {error}");
            }
        }
        previous_error = Some(error);
        close(state.axial_momentum_nms[0], p0 + tau * total_time, 5e-13);
        vector(
            state.rotation.rotate(state.total_momentum_body).unwrap(),
            l_world,
            2e-12,
        );
    }
    assert!(previous_error.unwrap() < 2e-5);
}

#[test]
fn external_torque_frames_work_and_reversed_steps_use_the_actual_rotor_trajectory() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(
            Vec3::new(1.5, 2.4, 3.2),
            rotation(Vec3::new(0.2, 0.3, -0.4)),
        )
        .unwrap(),
        [
            AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 0.4).unwrap(),
            AxialRotor::new(Vec3::new(0.2, 1.0, 0.4), 0.6).unwrap(),
        ],
    )
    .unwrap();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(-0.4, 0.1, 0.2)),
            Vec3::new(0.3, -0.6, 0.8),
            [1.0, -0.7],
            [0.2, -0.3],
        )
        .unwrap();
    let torque = Vec3::new(0.3, -0.2, 0.1);
    let shaft = [0.25, -0.4];
    let h = 0.08;
    for external in [HeldTorque::Body(torque), HeldTorque::Space(torque)] {
        let step = bank.step(&initial, shaft, external, h).unwrap();
        let after = step.state_after;
        let impulse = match external {
            HeldTorque::Body(value) => (initial.rotation.rotate(value).unwrap()
                + after.rotation.rotate(value).unwrap())
            .scale(0.5 * h),
            HeldTorque::Space(value) => value.scale(h),
        };
        vector(step.external_angular_impulse_world, impulse, 8e-15);
        vector(step.angular_momentum_residual_world, Vec3::default(), 3e-14);
        for i in 0..2 {
            close(
                step.shaft_work_j[i],
                shaft[i] * (after.relative_angle_rad[i] - initial.relative_angle_rad[i]),
                3e-14,
            );
        }
        close(
            step.kinetic_energy_change_j,
            step.shaft_work_j.iter().sum::<f64>()
                + step.external_work_j
                + step.drift_energy_defect_j,
            3e-14,
        );
        let back = bank.step(&after, shaft, external, -h).unwrap().state_after;
        attitude(back.rotation, initial.rotation, 3e-14);
        vector(back.total_momentum_body, initial.total_momentum_body, 3e-14);
        for i in 0..2 {
            close(
                back.axial_momentum_nms[i],
                initial.axial_momentum_nms[i],
                3e-14,
            );
            close(
                back.relative_angle_rad[i],
                initial.relative_angle_rad[i],
                3e-14,
            );
        }
    }
}

#[test]
fn unactuated_phase_refines_with_native_asymmetric_motion_and_keeps_drift_defect_visible() {
    let bank = RotorBank::from_reduced(
        InertiaTensor3::from_principal(Vec3::new(1.2, 2.0, 3.3), So3::identity()).unwrap(),
        [
            AxialRotor::new(Vec3::new(1.0, 0.2, 0.0), 0.4).unwrap(),
            AxialRotor::new(Vec3::new(-0.1, 0.3, 1.0), 0.7).unwrap(),
        ],
    )
    .unwrap();
    let initial = bank
        .state_from_velocities(
            rotation(Vec3::new(0.2, -0.3, 0.1)),
            Vec3::new(0.7, 0.3, -0.5),
            [1.2, -0.8],
            [0.0, 0.1],
        )
        .unwrap();
    let coarse = bank.step(&initial, [0.0; 2], zero_torque(), 0.2).unwrap();
    assert!(coarse.drift_energy_defect_j.abs() > 1e-7);
    assert_eq!(coarse.shaft_work_j, [0.0; 2]);
    close(
        coarse.kinetic_energy_change_j,
        coarse.drift_energy_defect_j,
        0.0,
    );
    close(coarse.work_energy_residual_j, 0.0, 1e-15);
    let run = |count: usize| {
        let mut state = initial;
        for _ in 0..count {
            state = bank
                .step(&state, [0.0; 2], zero_torque(), 0.8 / count as f64)
                .unwrap()
                .state_after;
        }
        state
    };
    // Independent exact spin-up is tested above. This is a refinement check
    // for the generic native splitting plus the new phase quadrature.
    let fine = run(4096);
    let mut previous = None;
    for count in [16, 32, 64] {
        let state = run(count);
        let error = state
            .rotation
            .body_minus(fine.rotation)
            .unwrap()
            .angular
            .norm()
            + (state.relative_angle_rad[0] - fine.relative_angle_rad[0]).abs()
            + (state.relative_angle_rad[1] - fine.relative_angle_rad[1]).abs();
        if let Some(old) = previous {
            assert!(old / error > 3.8 && old / error < 4.2, "{old} / {error}");
        }
        previous = Some(error);
        assert_eq!(state.axial_momentum_nms, initial.axial_momentum_nms);
    }
    let empty: RotorBank<0> = RotorBank::from_reduced(bank.reduced_inertia(), []).unwrap();
    let body = empty
        .state_from_velocities(initial.rotation, Vec3::new(0.2, 0.4, -0.3), [], [])
        .unwrap();
    let expected = empty
        .reduced_inertia()
        .momentum_step(
            body.rotation,
            body.total_momentum_body,
            Vec3::default(),
            0.05,
        )
        .unwrap();
    let actual = empty
        .step(&body, [], zero_torque(), 0.05)
        .unwrap()
        .state_after;
    assert_eq!((actual.rotation, actual.total_momentum_body), expected);
    let rest = bank
        .state_from_velocities(initial.rotation, Vec3::default(), [4.0, -3.0], [0.2, -0.3])
        .unwrap();
    let advanced = bank
        .step(&rest, [0.0; 2], zero_torque(), 0.1)
        .unwrap()
        .state_after;
    attitude(advanced.rotation, rest.rotation, 2e-14);
    vector(
        bank.angular_velocity(&advanced).unwrap(),
        Vec3::default(),
        2e-14,
    );
    close(advanced.relative_angle_rad[0], 0.6, 2e-14);
    close(advanced.relative_angle_rad[1], -0.6, 2e-14);
}

#[test]
fn aggregate_spin_must_resolve_simultaneous_carrier_recoil() {
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(2.0, 3.0, 4.0), So3::identity()).unwrap();
    let rotor = AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 1.0).unwrap();
    let bank = RotorBank::from_reduced(inertia, [rotor; 2]).unwrap();
    let initial = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [1e16, 0.0], [0.0; 2])
        .unwrap();
    let before = initial;
    // Both kicks are exactly representable in the second rotor. The native
    // aggregate h loses the first entirely and rounds the second to 2: neither
    // may silently lose/overstate the opposite carrier recoil.
    for impulse in [1.0, 1.5] {
        assert!(
            bank.apply_impulse(&initial, [0.0, impulse], Vec3::default())
                .is_err()
        );
    }
    // A simultaneous external impulse does not excuse the partial transfer.
    assert!(
        bank.apply_impulse(&initial, [0.0, 1.0], Vec3::new(2.0, 0.0, 0.0))
            .is_err()
    );
    assert!(
        bank.step(&initial, [0.0, 200.0], zero_torque(), 0.01)
            .is_err()
    );
    assert_eq!(initial, before);

    // The neighboring representable transfer remains usable. Independent
    // I_x=2, J=1 mechanics give omega=-1 and Delta T=1 carrier + 2 rotor = 3 J.
    // The complete locked tensor diag(4,3,4) is physically admissible too.
    let resolved = bank
        .apply_impulse(&initial, [0.0, 2.0], Vec3::default())
        .unwrap();
    vector(
        bank.angular_velocity(&resolved.state_after).unwrap(),
        Vec3::new(-1.0, 0.0, 0.0),
        0.0,
    );
    assert_eq!(
        resolved.state_after.total_momentum_body,
        initial.total_momentum_body
    );
    close(resolved.shaft_work_j[1], 3.0, 0.0);
    close(resolved.kinetic_energy_change_j, 3.0, 0.0);
    close(resolved.work_energy_residual_j, 0.0, 0.0);

    // Explicit lost-component refusal also applies when a huge carrier scale
    // would otherwise make a relative floating-point allowance permissive.
    let moving = bank
        .state_from_velocities(
            So3::identity(),
            Vec3::new(1e16, 0.0, 0.0),
            [0.0, -1e16],
            [0.0; 2],
        )
        .unwrap();
    assert!(
        bank.apply_impulse(&moving, [0.0, 1.0], Vec3::default())
            .is_err()
    );

    // Actual cancelling shaft transfers need no carrier recoil and are valid
    // even beside a huge passive spin. Their independently known work is 1 J.
    let cancelling = RotorBank::from_reduced(inertia, [rotor; 3]).unwrap();
    let rest = cancelling
        .state_from_velocities(So3::identity(), Vec3::default(), [1e16, 0.0, 0.0], [0.0; 3])
        .unwrap();
    let balanced = cancelling
        .apply_impulse(&rest, [0.0, 1.0, -1.0], Vec3::default())
        .unwrap();
    vector(
        cancelling.angular_velocity(&balanced.state_after).unwrap(),
        Vec3::default(),
        0.0,
    );
    close(balanced.shaft_work_j[1], 0.5, 0.0);
    close(balanced.shaft_work_j[2], 0.5, 0.0);
    close(balanced.kinetic_energy_change_j, 1.0, 0.0);
    close(balanced.work_energy_residual_j, 0.0, 0.0);
}

#[test]
fn invalid_inertia_and_unresolved_impulses_or_phase_refuse_without_mutation() {
    assert!(AxialRotor::new(Vec3::default(), 1.0).is_err());
    assert!(AxialRotor::new(Vec3::new(1.0, 0.0, 0.0), 0.0).is_err());
    let rotor = AxialRotor::new(Vec3::new(f64::from_bits(1), 0.0, 0.0), 1.0).unwrap();
    assert_eq!(rotor.axis_body(), Vec3::new(1.0, 0.0, 0.0));
    let inertia =
        InertiaTensor3::from_principal(Vec3::new(1.0, 1.0, 1.0), So3::identity()).unwrap();
    assert!(RotorBank::from_reduced(inertia, [rotor; MAX_AXIAL_ROTORS + 1]).is_err());
    assert!(
        RotorBank::from_locked(
            Mat3 {
                m: [0.5, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
            },
            [rotor]
        )
        .is_err()
    );
    let bank = RotorBank::from_reduced(inertia, [rotor]).unwrap();
    let ordinary = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [1.0], [0.0])
        .unwrap();
    let before = ordinary;
    assert!(bank.step(&ordinary, [0.0], zero_torque(), 0.0).is_err());
    assert!(
        bank.step(&ordinary, [f64::NAN], zero_torque(), 0.1)
            .is_err()
    );
    assert!(
        bank.apply_impulse(&ordinary, [f64::INFINITY], Vec3::default())
            .is_err()
    );
    let huge = bank
        .state_from_velocities(So3::identity(), Vec3::default(), [1e20], [0.0])
        .unwrap();
    assert!(bank.apply_impulse(&huge, [1.0], Vec3::default()).is_err());
    assert!(
        bank.apply_impulse(&huge, [0.0], Vec3::new(1.0, 0.0, 0.0))
            .is_err()
    );
    assert!(
        bank.state_from_velocities(So3::identity(), Vec3::new(1.0, 0.0, 0.0), [1e20], [0.0],)
            .is_err()
    );
    let mut phase = ordinary;
    phase.relative_angle_rad[0] = 1e30;
    assert!(bank.step(&phase, [0.0], zero_torque(), 0.1).is_err());
    assert_eq!(ordinary, before);
    assert_eq!(huge.axial_momentum_nms, [1e20]);
    assert_eq!(phase.relative_angle_rad, [1e30]);
}

#[test]
fn native_boundaries_preserve_original_axes_and_absolute_shaft_momenta() {
    let axis = Vec3::new(1.0 + 1e-13, 0.0, 0.0);
    let admitted = AxialRotor::from_unit_axis(axis, 0.17).unwrap();
    assert_eq!(admitted.axis_body(), axis);
    assert!(AxialRotor::from_unit_axis(Vec3::new(2.0, 0.0, 0.0), 0.17).is_err());
    assert!(AxialRotor::from_unit_axis(Vec3::default(), 0.17).is_err());
    assert!(AxialRotor::from_unit_axis(axis, 0.0).is_err());

    let inertia =
        InertiaTensor3::from_principal(Vec3::new(1.0, 1.5, 2.0), So3::identity()).unwrap();
    let rotor = AxialRotor::from_unit_axis(Vec3::new(1.0, 0.0, 0.0), 0.17).unwrap();
    let bank = RotorBank::from_reduced(inertia, [rotor]).unwrap();
    let omega = Vec3::new(0.31, 0.0, 0.0);
    let state = bank
        .state_from_axial_momenta(So3::identity(), omega, [0.123], [0.37])
        .unwrap();
    assert_eq!(state.axial_momentum_nms, [0.123]);
    assert_eq!(state.relative_angle_rad, [0.37]);
    let speed_roundtrip = 0.17 * ((0.123 / 0.17 - omega.x) + omega.x);
    assert_ne!(speed_roundtrip, state.axial_momentum_nms[0]);
    vector(bank.angular_velocity(&state).unwrap(), omega, 1e-15);
    assert!(
        bank.state_from_axial_momenta(
            So3::identity(),
            Vec3::new(1.0, 0.0, 0.0),
            [1e30],
            [0.0],
        )
        .is_err()
    );

    let full = RotorBank::from_reduced(inertia, [rotor; 64]).unwrap();
    let stationary = full
        .state_from_axial_momenta(So3::identity(), Vec3::default(), [0.0; 64], [0.0; 64])
        .unwrap();
    assert_eq!(
        full.step(&stationary, [0.0; 64], zero_torque(), 0.01)
            .unwrap()
            .state_after,
        stationary
    );
    assert!(RotorBank::from_reduced(inertia, [rotor; 65]).is_err());
}
