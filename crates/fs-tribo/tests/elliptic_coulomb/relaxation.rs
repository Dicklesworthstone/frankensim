//! G1/G3/G4: continuous directional slip, finite arrest and original history.
use super::{close, close2};
use fs_tribo::TangentialLaw;
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::elliptic::relaxation::*;

fn controls() -> EllipticRelaxationControls {
    EllipticRelaxationControls {
        time_tolerance_s: 1e-11,
        travel_tolerance_m: 1e-11,
        maximum_evaluations: 100_000,
        maximum_root_iterations: 64,
    }
}
fn input(v: [f64; 2], w: [[f64; 2]; 2], dt: f64) -> EllipticRelaxationInput {
    EllipticRelaxationInput {
        velocity_m_s: v,
        mobility_per_kg: w,
        dt_s: dt,
        normal_n: 10.0,
        maximum_force_n: 100.0,
        contact_temperature_k: Some(305.0),
        nominal_pressure_pa: Some(20_000.0),
    }
}
fn relax(law: &EllipticCoulomb, i: EllipticRelaxationInput) -> EllipticRelaxationReport {
    law.relax_velocity(&law.initial_state(), &i, controls(), &mut || false)
        .unwrap()
}

#[test]
fn isotropic_flow_stops_at_its_analytic_time_and_continues_without_extra_heat() {
    let law = EllipticCoulomb::new([0.9; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let i = input([3.0, 4.0], [[0.75, 0.0], [0.0, 0.75]], 10.0);
    let full = relax(&law, i);
    close(full.kinetic_phase_end_s.unwrap(), 10.0 / 3.0, 2e-15);
    assert_eq!(full.velocity_after_m_s, [0.0; 2]);
    close2(full.traction_impulse_ns, [4.0, 16.0 / 3.0], 2e-15);
    close2(full.slip_increment_m, [5.0, 20.0 / 3.0], 3e-15);
    close(full.distance_m, 25.0 / 3.0, 3e-15);
    close(full.dissipated_j, 50.0 / 3.0, 4e-15);
    assert_eq!(full.evaluations, 0);
    let first = relax(&law, EllipticRelaxationInput { dt_s: 1.0, ..i });
    close2(first.velocity_after_m_s, [2.1, 2.8], 1e-15);
    let second = law
        .relax_velocity(
            &first.state,
            &EllipticRelaxationInput {
                velocity_m_s: first.velocity_after_m_s,
                dt_s: 9.0,
                ..i
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    close2(second.state.slip_m(), full.state.slip_m(), 4e-15);
    close(
        first.dissipated_j + second.dissipated_j,
        full.dissipated_j,
        4e-15,
    );
    let rest = law
        .relax_velocity(
            &full.state,
            &EllipticRelaxationInput {
                velocity_m_s: full.velocity_after_m_s,
                ..i
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(rest.dissipated_j, 0.0);
    assert_eq!(rest.state, full.state);
    assert_eq!(law.static_coefficients(), [0.9; 2]);
}

fn clock_primitive(z: f64) -> f64 {
    0.5 * (z * (1.0 + z * z).sqrt() + z.asinh())
}
fn x_primitive(z: f64) -> f64 {
    (z * (1.0 + z * z).sqrt() * (1.0 + 2.0 * z * z) - z.asinh()) / 8.0
}
fn y_primitive(z: f64) -> f64 {
    2.0 * (1.0 + z * z).powf(1.5) / 3.0
}
fn path_reference(z: f64) -> f64 {
    // Independent composite Simpson integral of the analytic geometric curve.
    let n = 8192;
    let h = (1.0 - z) / n as f64;
    let f = |x: f64| x * ((x * x + 4.0) * (x * x + 1.0)).sqrt();
    let mut sum = f(z) + f(1.0);
    for k in 1..n {
        sum += if k % 2 == 0 {
            2.0 * f(z + k as f64 * h)
        } else {
            4.0 * f(z + k as f64 * h)
        };
    }
    sum * h / 3.0
}

#[test]
fn curved_flow_clock_impulse_displacement_and_heat_match_closed_continuum_solution() {
    // u=(z^2,z), dt/dz=-sqrt(1+z^2), v=(z^2,2z).
    let law = EllipticCoulomb::new([0.8; 2], [0.2, 0.1], [1.0, 0.0]).unwrap();
    for z in [0.0, 0.5, 0.9] {
        let elapsed = clock_primitive(1.0) - clock_primitive(z);
        let report = relax(
            &law,
            input(
                [1.0, 2.0],
                [[1.0, 0.0], [0.0, 2.0]],
                if z == 0.0 { 4.0 } else { elapsed },
            ),
        );
        close2(report.velocity_after_m_s, [z * z, 2.0 * z], 3e-11);
        close2(report.traction_impulse_ns, [1.0 - z * z, 1.0 - z], 3e-11);
        close2(
            report.slip_increment_m,
            [
                x_primitive(1.0) - x_primitive(z),
                y_primitive(1.0) - y_primitive(z),
            ],
            3e-11,
        );
        close(report.distance_m, path_reference(z), 3e-11);
        close(report.dissipated_j, 1.5 - 0.5 * z.powi(4) - z * z, 3e-11);
        assert!(report.distance_m > report.slip_increment_m[0].hypot(report.slip_increment_m[1]));
        if z == 0.0 {
            close(report.kinetic_phase_end_s.unwrap(), elapsed, 3e-12);
        } else {
            assert!(report.kinetic_phase_end_s.is_none());
        }
    }
}

#[test]
fn rotated_material_and_full_offdiagonal_mobility_preserve_the_curve_and_work() {
    let a = [0.6, 0.8];
    let rotate = |x: [f64; 2]| [a[0] * x[0] - a[1] * x[1], a[1] * x[0] + a[0] * x[1]];
    let w = [
        [a[0] * a[0] + 2.0 * a[1] * a[1], -a[0] * a[1]],
        [-a[0] * a[1], a[1] * a[1] + 2.0 * a[0] * a[0]],
    ];
    let law = EllipticCoulomb::new([0.8; 2], [0.2, 0.1], a).unwrap();
    let report = relax(&law, input(rotate([1.0, 2.0]), w, 0.4));
    let original = EllipticCoulomb::new([0.8; 2], [0.2, 0.1], [1.0, 0.0]).unwrap();
    let reference = relax(&original, input([1.0, 2.0], [[1.0, 0.0], [0.0, 2.0]], 0.4));
    close2(
        report.velocity_after_m_s,
        rotate(reference.velocity_after_m_s),
        3e-11,
    );
    close2(
        report.traction_impulse_ns,
        rotate(reference.traction_impulse_ns),
        3e-11,
    );
    close2(
        report.slip_increment_m,
        rotate(reference.slip_increment_m),
        3e-11,
    );
    close(report.dissipated_j, reference.dissipated_j, 3e-11);
    close(report.distance_m, reference.distance_m, 3e-11);
    let continuation = law
        .relax_velocity(
            &report.state,
            &input(report.velocity_after_m_s, w, 4.0),
            controls(),
            &mut || false,
        )
        .unwrap();
    let full = relax(&law, input(rotate([1.0, 2.0]), w, 4.4));
    close2(continuation.state.slip_m(), full.state.slip_m(), 4e-11);
    close(
        report.distance_m + continuation.distance_m,
        full.distance_m,
        4e-11,
    );
    close(
        report.dissipated_j + continuation.dissipated_j,
        full.dissipated_j,
        4e-11,
    );
}

#[test]
fn one_frictionless_axis_retains_cross_recoil_and_free_motion_after_active_arrest() {
    let law = EllipticCoulomb::new([0.8, 0.4], [0.2, 0.0], [1.0, 0.0]).unwrap();
    let report = relax(&law, input([2.0, 4.0], [[2.0, 1.0], [1.0, 3.0]], 1.0));
    close(report.kinetic_phase_end_s.unwrap(), 0.5, 1e-15);
    close2(report.velocity_after_m_s, [0.0, 3.0], 1e-15);
    close2(report.traction_impulse_ns, [1.0, 0.0], 1e-15);
    close2(report.slip_increment_m, [0.5, 3.25], 1e-15);
    close(report.dissipated_j, 1.0, 1e-15);
    assert!(report.distance_m > 0.5f64.hypot(3.25));
    let free = relax(&law, input([0.0, 3.0], [[2.0, 1.0], [1.0, 3.0]], 1.0));
    assert_eq!(free.dissipated_j, 0.0);
    assert_eq!(free.slip_increment_m, [0.0, 3.0]);
    assert_eq!(free.evaluations, 0);
    let unloaded = relax(
        &law,
        EllipticRelaxationInput {
            normal_n: 0.0,
            ..input([2.0, 4.0], [[2.0, 1.0], [1.0, 3.0]], 2.0)
        },
    );
    assert_eq!(unloaded.velocity_after_m_s, [2.0, 4.0]);
    assert_eq!(unloaded.slip_increment_m, [4.0, 8.0]);
    assert_eq!(unloaded.dissipated_j, 0.0);
}

#[test]
fn cancellation_budgets_invalid_mobility_and_foreign_history_publish_no_candidate() {
    let law = EllipticCoulomb::new([0.8; 2], [0.2, 0.1], [1.0, 0.0]).unwrap();
    let state = law.initial_state();
    let before = state.clone();
    let i = input([1.0, 2.0], [[1.0, 0.0], [0.0, 2.0]], 0.3);
    let mut calls = 0;
    assert!(matches!(
        law.relax_velocity(&state, &i, controls(), &mut || {
            calls += 1;
            calls > 20
        }),
        Err(EllipticRelaxationError::Cancelled)
    ));
    assert!(matches!(
        law.relax_velocity(
            &state,
            &i,
            EllipticRelaxationControls {
                maximum_evaluations: 4,
                ..controls()
            },
            &mut || false
        ),
        Err(EllipticRelaxationError::EvaluationBudget { .. })
    ));
    assert!(matches!(
        law.relax_velocity(
            &state,
            &i,
            EllipticRelaxationControls {
                maximum_root_iterations: 1,
                ..controls()
            },
            &mut || false
        ),
        Err(EllipticRelaxationError::UnresolvedTime)
    ));
    assert!(
        law.relax_velocity(
            &state,
            &EllipticRelaxationInput {
                mobility_per_kg: [[1.0, 1.0], [1.0, 1.0]],
                ..i
            },
            controls(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        law.relax_velocity(
            &state,
            &EllipticRelaxationInput {
                maximum_force_n: 1.0,
                ..i
            },
            controls(),
            &mut || false
        )
        .is_err()
    );
    let foreign = EllipticCoulomb::new([0.9; 2], [0.2, 0.1], [1.0, 0.0])
        .unwrap()
        .initial_state();
    assert!(
        law.relax_velocity(&foreign, &i, controls(), &mut || false)
            .is_err()
    );
    assert_eq!(state, before);
    assert_eq!(
        law.relax_velocity(&state, &i, controls(), &mut || false)
            .unwrap(),
        law.relax_velocity(&state, &i, controls(), &mut || false)
            .unwrap()
    );
}

#[test]
fn short_time_and_stiff_anisotropy_retain_small_impulse_and_physical_work() {
    let law = EllipticCoulomb::new([2.0; 2], [1.0, 0.01], [1.0, 0.0]).unwrap();
    let i = input([1.0, 1.0], [[1.0, 0.0], [0.0, 1.0]], 1e-6);
    let r = law
        .relax_velocity(
            &law.initial_state(),
            &i,
            EllipticRelaxationControls {
                time_tolerance_s: 1e-16,
                travel_tolerance_m: 1e-12,
                maximum_evaluations: 500_000,
                ..controls()
            },
            &mut || false,
        )
        .unwrap();
    let fx = 100.0 / (100.0f64 + 0.01).sqrt();
    let fy = 0.01 / (100.0f64 + 0.01).sqrt();
    close(r.traction_impulse_ns[0], fx * i.dt_s, 2e-12);
    close(r.traction_impulse_ns[1], fy * i.dt_s, 2e-13);
    assert!(r.dissipated_j > 0.0);
    assert!(r.kinetic_phase_end_s.is_none());
    close(
        r.dissipated_j,
        0.5 * r.traction_impulse_ns[0] * (1.0 + r.velocity_after_m_s[0])
            + 0.5 * r.traction_impulse_ns[1] * (1.0 + r.velocity_after_m_s[1]),
        1e-13,
    );
}

#[test]
fn ordinary_anisotropy_resolves_fractional_endpoint_and_finite_stop() {
    // The original source fixture has caps=(4,3), W=diag(1,2), v=(.4,.1).
    // Its eigenvalue ratio is 9/8. The original z chart has an endpoint term
    // z^(1/4), which exhausted depth 48 even with ordinary coefficients.
    // Independent reference: q=z^(1/8), v=(.4*q^8,.1*q^9),
    // dt/dq=2*q^7*sqrt(.16+.005625*q^2). The time primitive is the polynomial
    // 2/B^4*[t^9/9-3*A*t^7/7+3*A^2*t^5/5-A^3*t^3/3],
    // t=sqrt(A+B*q^2), A=.16, B=.005625. Clock inversion used 65-digit
    // decimal arithmetic; independent 16384/32768-panel Simpson displacement
    // and path references differed by <8e-17 m.
    let law = EllipticCoulomb::new([0.6, 0.7], [0.4, 0.3], [1.0, 0.0]).unwrap();
    for (dt, velocity, displacement, distance) in [
        (
            0.1,
            [0.005557644962472804, 0.0008141129425462053],
            [0.0203061852166916, 0.004778781050303254],
            0.02086270525153858,
        ),
        (
            1.0,
            [0.0, 0.0],
            [0.020310066799811902, 0.00477931621688067],
            0.020866623678189287,
        ),
    ] {
        let report = relax(&law, input([0.4, 0.1], [[1.0, 0.0], [0.0, 2.0]], dt));
        close2(report.velocity_after_m_s, velocity, 3e-11);
        close2(report.slip_increment_m, displacement, 3e-11);
        close(report.distance_m, distance, 3e-11);
        close2(
            report.traction_impulse_ns,
            [0.4 - velocity[0], (0.1 - velocity[1]) / 2.0],
            3e-11,
        );
        close(
            report.dissipated_j,
            0.0825 - velocity[0] * velocity[0] / 2.0 - velocity[1] * velocity[1] / 4.0,
            3e-11,
        );
        if dt == 1.0 {
            close(
                report.kinetic_phase_end_s.unwrap(),
                0.10139610258249089,
                3e-12,
            );
        } else {
            assert!(report.kinetic_phase_end_s.is_none());
        }
    }
}
