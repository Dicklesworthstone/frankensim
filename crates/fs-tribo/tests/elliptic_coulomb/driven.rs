//! G1/G3/G4: original directional law under held external drive.
use super::{close, close2};
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::elliptic::driven::*;
use fs_tribo::tangential::elliptic::relaxation::*;
use fs_tribo::{TangentialLaw, TangentialSolveConfig};

fn controls(h: f64) -> EllipticDrivenControls {
    EllipticDrivenControls {
        maximum_step_s: h,
        minimum_step_s: 1e-14,
        maximum_substeps: 20_000,
        rest_velocity_tolerance_m_s: 1e-12,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
        relaxation: EllipticRelaxationControls {
            time_tolerance_s: 1e-11,
            travel_tolerance_m: 1e-11,
            maximum_evaluations: 100_000,
            maximum_root_iterations: 96,
        },
    }
}
fn input(v: [f64; 2], w: [[f64; 2]; 2], a: [f64; 2], dt: f64) -> EllipticDrivenInput {
    EllipticDrivenInput {
        velocity_m_s: v,
        mobility_per_kg: w,
        driving_acceleration_m_s2: a,
        mode: if v == [0.0; 2] {
            EllipticDrivenMode::Resting
        } else {
            EllipticDrivenMode::Sliding
        },
        normal_n: 10.0,
        dt_s: dt,
        maximum_force_n: 100.0,
        contact_temperature_k: Some(305.0),
        nominal_pressure_pa: Some(2000.0),
    }
}
fn apply(law: &EllipticCoulomb, i: EllipticDrivenInput, h: f64) -> EllipticDrivenReport {
    law.advance_driven(&law.initial_state(), &i, controls(h), &mut || false)
        .unwrap()
}
fn mat(w: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [
        w[0][0] * x[0] + w[0][1] * x[1],
        w[1][0] * x[0] + w[1][1] * x[1],
    ]
}
fn dot(x: [f64; 2], y: [f64; 2]) -> f64 {
    x[0] * y[0] + x[1] * y[1]
}
fn norm(x: [f64; 2]) -> f64 {
    x[0].hypot(x[1])
}

#[test]
fn full_matrix_static_hold_and_breakaway_preserve_the_explicit_rest_mode() {
    let axis = [0.6, 0.8];
    let rotate = |x: [f64; 2]| {
        [
            axis[0] * x[0] - axis[1] * x[1],
            axis[1] * x[0] + axis[0] * x[1],
        ]
    };
    let w = [[1.0, 0.3], [0.3, 2.0]];
    let law = EllipticCoulomb::new([0.6, 0.4], [0.3, 0.2], axis).unwrap();
    let q = rotate([4.0, 1.0]); // Inside static, outside kinetic.
    let mut i = input([1e-15, -2e-15], w, mat(w, q), 0.4);
    i.mode = EllipticDrivenMode::Resting;
    let r = apply(&law, i, 0.01);
    assert_eq!(r.mode, EllipticDrivenMode::Resting);
    assert_eq!(r.velocity_after_m_s, [0.0; 2]);
    assert_eq!(r.slip_increment_m, [0.0; 2]);
    assert_eq!(r.distance_m, 0.0);
    assert_eq!(r.dissipated_j, 0.0);
    assert_eq!(r.drive_work_j, 0.0);
    assert_eq!(r.energy_residual_j, 0.0);
    close2(r.traction_impulse_ns, q.map(|x| x * i.dt_s), 1e-15);
    close2(r.driving_force_n, q, 2e-15);
    assert_eq!(r.substeps, 1);
    assert_eq!(r.first_arrest_time_s, None);
    // The same roundoff-sized velocity with Sliding authority is not silently
    // promoted to static hold in the kinetic/static hysteresis band.
    i.mode = EllipticDrivenMode::Sliding;
    let sliding = law.advance_driven(&law.initial_state(), &i, controls(0.01), &mut || false);
    if let Ok(sliding) = sliding {
        assert_eq!(sliding.mode, EllipticDrivenMode::Sliding);
    }

    let q = rotate([8.0, 3.0]);
    let i = input([0.0; 2], w, mat(w, q), 0.2);
    let r = apply(&law, i, 0.001);
    assert_eq!(r.mode, EllipticDrivenMode::Sliding);
    assert_eq!(r.substeps, 1); // Homogeneous outgoing acceleration is exact.
    let v = r.velocity_after_m_s;
    let principal = [
        axis[0] * v[0] + axis[1] * v[1],
        -axis[1] * v[0] + axis[0] * v[1],
    ];
    let denominator = (3.0 * principal[0]).hypot(2.0 * principal[1]);
    let force = rotate([
        9.0 * principal[0] / denominator,
        4.0 * principal[1] / denominator,
    ]);
    close2(r.traction_impulse_ns, force.map(|x| x * i.dt_s), 3e-12);
    let response = mat(w, force);
    close2(
        v,
        core::array::from_fn(|j| i.dt_s * (i.driving_acceleration_m_s2[j] - response[j])),
        3e-12,
    );
    close2(r.slip_increment_m, v.map(|x| 0.5 * i.dt_s * x), 2e-12);
    close(r.dissipated_j, dot(force, r.slip_increment_m), 3e-12);
    close(r.drive_work_j, dot(q, r.slip_increment_m), 3e-12);
    close(r.energy_residual_j, 0.0, 3e-12);
}

#[test]
fn exact_opposed_motion_stops_then_holds_or_restarts_under_the_original_drive() {
    let law = EllipticCoulomb::new([0.6; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let w = [[0.75, 0.0], [0.0, 1.5]];
    for (q, duration, holds) in [(3.0, 2.0, true), (8.0, 1.0, false)] {
        let stop = 2.0 / (0.75 * (q + 2.0));
        let remaining = duration - stop;
        let acceleration_after = if holds { 0.0 } else { 0.75 * (q - 2.0) };
        let forward = 0.5 * acceleration_after * remaining * remaining;
        let x = -stop + forward;
        let distance = stop + forward;
        let r = apply(
            &law,
            input([-2.0, 0.0], w, [0.75 * q, 0.0], duration),
            0.001,
        );
        close(r.first_arrest_time_s.unwrap(), stop, 1e-15);
        assert!(r.arrest_time_is_exact);
        assert_eq!(
            r.mode,
            if holds {
                EllipticDrivenMode::Resting
            } else {
                EllipticDrivenMode::Sliding
            }
        );
        close2(
            r.velocity_after_m_s,
            [acceleration_after * remaining, 0.0],
            3e-15,
        );
        close2(r.slip_increment_m, [x, 0.0], 3e-15);
        close(r.distance_m, distance, 3e-15);
        close(r.dissipated_j, 2.0 * distance, 3e-15);
        close(r.drive_work_j, q * x, 6e-15);
        close(
            r.kinetic_energy_change_j,
            (r.velocity_after_m_s[0].powi(2) - 4.0) / 1.5,
            6e-15,
        );
        close(r.energy_residual_j, 0.0, 8e-15);
        assert_eq!(r.substeps, 2);
        // Exact integrated momentum checks include the static reaction after stop.
        close(
            r.traction_impulse_ns[0],
            (0.75 * q * duration - 2.0 - r.velocity_after_m_s[0]) / 0.75,
            3e-15,
        );
    }
}

// Independent explicit RK4 integration of the continuous physical differential
// equation, augmented with displacement, unsigned distance and friction work.
// This uses neither the native endpoint projection nor its substep trajectory.
fn reference(w: [[f64; 2]; 2], a: [f64; 2], v: [f64; 2], duration: f64, n: usize) -> [f64; 6] {
    let f = |y: [f64; 6]| {
        let denominator = (2.0 * y[0]).hypot(4.0 * y[1]);
        let force = [4.0 * y[0] / denominator, 16.0 * y[1] / denominator];
        let response = mat(w, force);
        [
            a[0] - response[0],
            a[1] - response[1],
            y[0],
            y[1],
            y[0].hypot(y[1]),
            force[0] * y[0] + force[1] * y[1],
        ]
    };
    let mut y = [v[0], v[1], 0.0, 0.0, 0.0, 0.0];
    let h = duration / n as f64;
    for _ in 0..n {
        let k1 = f(y);
        let k2 = f(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = f(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = f(core::array::from_fn(|i| y[i] + h * k3[i]));
        y = core::array::from_fn(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
    }
    y
}

#[test]
fn full_spd_curved_driven_motion_has_a_three_rung_continuum_ladder_and_exact_discrete_work() {
    let law = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let i = input([1.1, 0.7], [[0.8, 0.25], [0.25, 1.2]], [0.4, -0.3], 0.2);
    let oracle = reference(
        i.mobility_per_kg,
        i.driving_acceleration_m_s2,
        i.velocity_m_s,
        i.dt_s,
        4096,
    );
    let check = reference(
        i.mobility_per_kg,
        i.driving_acceleration_m_s2,
        i.velocity_m_s,
        i.dt_s,
        8192,
    );
    for j in 0..6 {
        close(oracle[j], check[j], 1e-12);
    }
    let mut errors = Vec::new();
    for h in [0.05, 0.025, 0.0125] {
        let r = apply(&law, i, h);
        errors.push([
            norm([
                r.velocity_after_m_s[0] - oracle[0],
                r.velocity_after_m_s[1] - oracle[1],
            ]),
            norm([
                r.slip_increment_m[0] - oracle[2],
                r.slip_increment_m[1] - oracle[3],
            ]),
            (r.distance_m - oracle[4]).abs(),
            (r.dissipated_j - oracle[5]).abs(),
        ]);
        assert!(r.distance_m > norm(r.slip_increment_m));
        close(r.energy_residual_j, 0.0, 2e-11);
        close(
            r.drive_work_j,
            dot(r.driving_force_n, r.slip_increment_m),
            2e-14,
        );
        let momentum = mat(i.mobility_per_kg, r.traction_impulse_ns);
        close2(
            r.velocity_after_m_s,
            core::array::from_fn(|j| {
                i.velocity_m_s[j] + i.dt_s * i.driving_acceleration_m_s2[j] - momentum[j]
            }),
            3e-12,
        );
        assert_eq!(r.first_arrest_time_s, None);
        assert_eq!(r.mode, EllipticDrivenMode::Sliding);
    }
    for j in 0..4 {
        assert!(errors[2][j] < 0.02, "refinement errors {errors:?}");
        for pair in errors.windows(2) {
            let ratio = pair[0][j] / pair[1][j];
            assert!(
                ratio > 1.4 && ratio < 2.6,
                "first-order driven continuum ladder {errors:?}"
            );
        }
    }
}

#[test]
fn transverse_superkinetic_drive_does_not_fabricate_static_capture() {
    let law = EllipticCoulomb::new([0.6; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let i = input([-1.0, 0.1], [[1.0, 0.0], [0.0, 1.0]], [3.0, 0.0], 1.0);
    let r = apply(&law, i, 0.2);
    assert_eq!(r.first_arrest_time_s, None);
    assert_eq!(r.mode, EllipticDrivenMode::Sliding);
    assert!(r.velocity_after_m_s[0] > 0.5);
    assert!(r.substeps > 5); // Coarse false-arrest/assisting-force proposals refine.
    close(r.energy_residual_j, 0.0, 1e-11);
}

#[test]
fn frame_rotation_zero_drive_and_rank_one_limits_reuse_the_original_physics() {
    let law = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let i = input([1.1, 0.7], [[0.8, 0.25], [0.25, 1.2]], [0.4, -0.3], 0.2);
    let base = apply(&law, i, 0.025);
    let axis = [0.6, 0.8];
    let rotate = |x: [f64; 2]| {
        [
            axis[0] * x[0] - axis[1] * x[1],
            axis[1] * x[0] + axis[0] * x[1],
        ]
    };
    let [c, s] = axis;
    let [[a, b], [_, d]] = i.mobility_per_kg;
    let off = c * s * (a - d) + (c * c - s * s) * b;
    let rotated = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], axis).unwrap();
    let r = apply(
        &rotated,
        EllipticDrivenInput {
            velocity_m_s: rotate(i.velocity_m_s),
            driving_acceleration_m_s2: rotate(i.driving_acceleration_m_s2),
            mobility_per_kg: [
                [c * c * a - 2.0 * c * s * b + s * s * d, off],
                [off, s * s * a + 2.0 * c * s * b + c * c * d],
            ],
            ..i
        },
        0.025,
    );
    close2(r.velocity_after_m_s, rotate(base.velocity_after_m_s), 3e-11);
    close2(r.slip_increment_m, rotate(base.slip_increment_m), 3e-11);
    close2(
        r.traction_impulse_ns,
        rotate(base.traction_impulse_ns),
        3e-11,
    );
    close(r.distance_m, base.distance_m, 3e-11);
    close(r.dissipated_j, base.dissipated_j, 3e-11);
    close(r.drive_work_j, base.drive_work_j, 3e-11);

    let rest = EllipticDrivenInput {
        driving_acceleration_m_s2: [0.0; 2],
        ..i
    };
    let r = apply(&law, rest, 0.001);
    let native = law
        .relax_velocity(
            &law.initial_state(),
            &EllipticRelaxationInput {
                velocity_m_s: rest.velocity_m_s,
                mobility_per_kg: rest.mobility_per_kg,
                normal_n: rest.normal_n,
                dt_s: rest.dt_s,
                maximum_force_n: rest.maximum_force_n,
                contact_temperature_k: rest.contact_temperature_k,
                nominal_pressure_pa: rest.nominal_pressure_pa,
            },
            controls(0.001).relaxation,
            &mut || false,
        )
        .unwrap();
    assert_eq!(r.velocity_after_m_s, native.velocity_after_m_s);
    assert_eq!(r.traction_impulse_ns, native.traction_impulse_ns);
    assert_eq!(r.slip_increment_m, native.slip_increment_m);
    assert_eq!(r.distance_m, native.distance_m);
    assert_eq!(r.dissipated_j, native.dissipated_j);
    assert_eq!(r.state, native.state);

    let rank = EllipticCoulomb::new([0.5, 0.0], [0.2, 0.0], [1.0, 0.0]).unwrap();
    let r = apply(
        &rank,
        input([1.0, 2.0], [[1.0, 0.3], [0.3, 1.0]], [0.1, -0.4], 0.2),
        0.05,
    );
    close2(r.velocity_after_m_s, [0.62, 1.8], 2e-15);
    close2(r.slip_increment_m, [0.162, 0.38], 2e-15);
    close2(r.traction_impulse_ns, [0.4, 0.0], 2e-15);
    close(r.dissipated_j, 0.324, 2e-15);
    close(r.energy_residual_j, 0.0, 3e-15);
}

#[test]
fn invalid_mode_mobility_budgets_and_cancellation_refuse_without_history_change() {
    let law = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let initial = law.initial_state();
    let before = initial.clone();
    let original = input([1.1, 0.7], [[0.8, 0.25], [0.25, 1.2]], [0.4, -0.3], 0.2);
    for mode in 0..6 {
        let mut i = original;
        let mut c = controls(0.01);
        match mode {
            0 => i.mode = EllipticDrivenMode::Resting,
            1 => i.velocity_m_s = [0.0; 2],
            2 => i.mobility_per_kg = [[1.0, 1.0], [1.0, 1.0]],
            3 => c.maximum_step_s = f64::NAN,
            4 => c.maximum_substeps = 1,
            _ => i.maximum_force_n = 1.0,
        }
        assert!(law.advance_driven(&initial, &i, c, &mut || false).is_err());
        assert_eq!(initial, before);
    }
    let mut polls = 0;
    let refusal = law.advance_driven(&initial, &original, controls(0.01), &mut || {
        polls += 1;
        polls > 5
    });
    assert!(matches!(
        refusal,
        Err(EllipticDrivenError::Cancelled)
            | Err(EllipticDrivenError::Relaxation(
                EllipticRelaxationError::Cancelled
            ))
    ));
    assert_eq!(initial, before);
    let retry = law
        .advance_driven(&initial, &original, controls(0.01), &mut || false)
        .unwrap();
    assert_ne!(retry.state, initial);
    close(retry.energy_residual_j, 0.0, 1e-11);
}
