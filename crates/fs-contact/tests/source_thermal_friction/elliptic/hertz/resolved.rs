//! G1/G3/G4: actual local source speeds, original receipts and all-or-none output.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::hertz::resolved::{
    SourceResolvedHertzError, SourceResolvedHertzReport,
};
use fs_tribo::sliding_patch::elliptic::resolved::{ResolvedHertzError, ResolvedHertzQuadrature};

fn controls(radial: usize, angular: usize) -> ResolvedHertzQuadrature {
    ResolvedHertzQuadrature {
        radial_panels: radial,
        angular_panels: angular,
        maximum_evaluations: 16_384,
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}

fn verify_local_receipts(
    law: &SourceThermalEllipticCoulomb,
    state: &SourceThermalEllipticState,
    report: &SourceResolvedHertzReport,
    input: SourceEllipticHertzInput<'_>,
) {
    assert_eq!(report.samples.len(), report.response.evaluations);
    let mut minimum = f64::INFINITY;
    let mut maximum = 0.0_f64;
    for sample in &report.samples {
        let point = sample.point;
        let expected_velocity = [
            input.slip_velocity_m_s[0] - input.spin_rad_s * point.position_m[1],
            input.slip_velocity_m_s[1] + input.spin_rad_s * point.position_m[0],
        ];
        for i in 0..2 {
            near(point.slip_velocity_m_s[i], expected_velocity[i], 3e-16);
        }
        let speed = point.slip_velocity_m_s[0].hypot(point.slip_velocity_m_s[1]);
        minimum = minimum.min(speed);
        maximum = maximum.max(speed);
        near(
            sample.law.initial_state().principal_axis()[0],
            state.constitutive_state().principal_axis()[0],
            2e-16,
        );
        for (index, source) in sample.sources.iter().enumerate() {
            near(source.state_coordinate("speed").unwrap(), speed, 5e-16);
            assert_eq!(source.state_coordinate("Ta"), Some(290.0));
            assert_eq!(source.state_coordinate("Tb"), Some(310.0));
            assert_eq!(source.state_coordinate("T"), Some(300.0));
            assert_eq!(source.state_coordinate("path"), Some(0.25));
            for property in source.properties() {
                let original = law.principal_sources()[index]
                    .initial_source()
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    original.answer().receipt.selected
                );
                law.principal_sources()[index]
                    .interface()
                    .card()
                    .claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
    }
    assert!(maximum > minimum && minimum >= report.speed_envelope_m_s[0]);
    assert!(maximum <= report.speed_envelope_m_s[1]);
}

#[test]
fn source_rate_spin_and_shear_recover_independent_hertz_second_moments() {
    let beta = 0.2;
    let (law, normal) =
        fixture_with_rate_pair([290.0, 310.0], [0.6, 0.8], true, false, Some([beta; 2]));
    let state = law.initial_state();
    let saved = state.clone();
    let mut input = motion();
    input.indentation_rate_m_s = 0.0;
    for slip in [[0.0; 2], [0.003, -0.002], [0.03, -0.02]] {
        input.slip_velocity_m_s = slip;
        let report = law
            .evaluate_hertz_patch_resolved(
                &state,
                &normal,
                input,
                &conditions(true),
                controls(8, 16),
                &mut || false,
            )
            .unwrap();
        let axes = report.normal.elliptic_patch_axes.unwrap();
        let second = (axes.semi_major_axis_m.powi(2) + axes.semi_minor_axis_m.powi(2)) / 5.0;
        let scale = beta * report.normal.normal_force_n;
        for i in 0..2 {
            near(report.response.wrench.force_n[i], -scale * slip[i], 1e-12);
        }
        near(
            report.response.wrench.spin_moment_nm,
            -scale * input.spin_rad_s * second,
            1e-15,
        );
        near(
            report.response.wrench.dissipated_power_w,
            scale * (slip[0].powi(2) + slip[1].powi(2) + input.spin_rad_s.powi(2) * second),
            1e-14,
        );
        near(report.response.normal_force_residual_n, 0.0, 1e-10);
        near(report.response.wrench.power_residual_w, 0.0, 1e-18);
        verify_local_receipts(&law, &state, &report, input);
        for bound in &report.speed_envelope_samples {
            verify_sources(&law, &bound.state, [290.0, 310.0]);
            assert_eq!(
                bound.state.constitutive_state().slip_m(),
                saved.constitutive_state().slip_m()
            );
        }
        assert_eq!(state, saved);
    }
}

#[test]
fn anisotropic_rate_sources_use_transported_texture_and_rotate_the_complete_wrench() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], true, false);
    let initial = law.initial_state();
    let input = motion();
    let c = controls(4, 8);
    let a = law
        .evaluate_hertz_patch_resolved(&initial, &normal, input, &conditions(true), c, &mut || {
            false
        })
        .unwrap();
    let old = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let new = old.rotated(0.7).unwrap();
    let rotation = TangentialHistoryRotation::between(old, new, 0.0).unwrap();
    let moved = law.transport_state(&initial, rotation).unwrap();
    let rotate = |v: [f64; 2]| {
        let (s, c) = 0.7_f64.sin_cos();
        [c * v[0] + s * v[1], -s * v[0] + c * v[1]]
    };
    let transformed = SourceEllipticHertzInput {
        footprint_axis: rotate(input.footprint_axis),
        slip_velocity_m_s: rotate(input.slip_velocity_m_s),
        ..input
    };
    let b = law
        .evaluate_hertz_patch_resolved(
            &moved,
            &normal,
            transformed,
            &conditions(true),
            c,
            &mut || false,
        )
        .unwrap();
    let force = rotate(a.response.wrench.force_n);
    for i in 0..2 {
        near(b.response.wrench.force_n[i], force[i], 2e-13);
    }
    near(
        b.response.wrench.spin_moment_nm,
        a.response.wrench.spin_moment_nm,
        2e-15,
    );
    near(
        b.response.wrench.dissipated_power_w,
        a.response.wrench.dissipated_power_w,
        2e-14,
    );
    verify_local_receipts(&law, &moved, &b, transformed);
    // Compare both coefficients at actual local speed, not the source template's 99 m/s.
    for sample in &b.samples {
        let speed = sample.sources[0].state_coordinate("speed").unwrap();
        near(
            sample.law.kinetic_coefficients()[0],
            0.2 + 0.2 * speed,
            2e-15,
        );
        near(sample.law.kinetic_coefficients()[1], 0.4, 2e-15);
    }
}

#[test]
fn entire_patch_speed_admission_budget_and_late_cancellation_preserve_accepted_history() {
    let (law, normal) = fixture([290.0, 310.0], [1.0, 0.0], true, false);
    let state = law.initial_state();
    let saved = state.clone();
    let mut input = motion();
    input.slip_velocity_m_s = [1.999, 0.0];
    assert!(
        law.evaluate_hertz_patch_resolved(
            &state,
            &normal,
            input,
            &conditions(true),
            controls(1, 1),
            &mut || false,
        )
        .is_err()
    );
    input.slip_velocity_m_s = [1.0, 0.0];
    input.spin_rad_s = 1e-200;
    assert!(matches!(
        law.evaluate_hertz_patch_resolved(
            &state,
            &normal,
            input,
            &conditions(true),
            controls(1, 1),
            &mut || false,
        ),
        Err(SourceResolvedHertzError::Invalid(
            "local-speed variation is not representable"
        ))
    ));
    input = motion();
    let mut c = controls(1, 2);
    c.maximum_evaluations = 1;
    assert!(matches!(
        law.evaluate_hertz_patch_resolved(
            &state,
            &normal,
            input,
            &conditions(true),
            c,
            &mut || false,
        ),
        Err(SourceResolvedHertzError::Native(
            ResolvedHertzError::EvaluationBudget { .. }
        ))
    ));
    let c = controls(1, 2);
    let mut polls = 0;
    let a = law
        .evaluate_hertz_patch_resolved(&state, &normal, input, &conditions(true), c, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let target = polls - 2;
    polls = 0;
    assert!(matches!(
        law.evaluate_hertz_patch_resolved(
            &state,
            &normal,
            input,
            &conditions(true),
            c,
            &mut || {
                polls += 1;
                polls == target
            },
        ),
        Err(SourceResolvedHertzError::Cancelled)
    ));
    let b = law
        .evaluate_hertz_patch_resolved(&state, &normal, input, &conditions(true), c, &mut || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(state, saved);
    input.indentation_m = 0.0;
    let empty = law
        .evaluate_hertz_patch_resolved(&state, &normal, input, &conditions(true), c, &mut || false)
        .unwrap();
    assert!(empty.patch.is_none() && empty.samples.is_empty());
    assert_eq!(empty.response.wrench.force_n, [0.0; 2]);
    assert_eq!(empty.response.wrench.spin_moment_nm, 0.0);
    assert_eq!(empty.response.wrench.dissipated_power_w, 0.0);
    let (foreign, foreign_normal) = fixture([290.0, 310.0], [1.0, 0.0], true, true);
    assert!(
        law.evaluate_hertz_patch_resolved(
            &state,
            &foreign_normal,
            input,
            &conditions(true),
            c,
            &mut || false,
        )
        .is_err()
    );
    assert!(
        law.evaluate_hertz_patch_resolved(
            &foreign.initial_state(),
            &normal,
            input,
            &conditions(true),
            c,
            &mut || false,
        )
        .is_err()
    );
}
