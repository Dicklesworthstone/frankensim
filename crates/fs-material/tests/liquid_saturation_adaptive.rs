//! Adaptive capillary diffusion against an independent nonlinear transient and
//! finite-volume water balances. The retention parameters are synthetic.

use fs_material::porous::wet::liquid::{
    WetPoreLiquidLaw, WetPoreLiquidSupport,
    transient::{
        LiquidSaturationAdvanceControls as Accuracy, LiquidSaturationAdvanceError as Error,
        LiquidSaturationBoundary as Boundary, LiquidSaturationBoundaryStepInput as Input,
        LiquidSaturationCell as Cell, LiquidSaturationColumn as Column,
        LiquidSaturationImplicitControls as Nonlinear, LiquidSaturationImplicitError,
    },
};

fn cell(thickness_m: f64, porosity: f64, permeability: f64) -> Cell {
    Cell {
        thickness_m,
        porosity,
        intrinsic_permeability_m2: permeability,
        law: WetPoreLiquidLaw::try_new(
            0.1,
            0.1,
            0.5,
            0.5,
            10_000.0,
            WetPoreLiquidSupport {
                liquid_saturation: (0.1, 0.9),
                temperature_k: (280.0, 330.0),
                liquid_pressure_pa: (50_000.0, 200_000.0),
                maximum_superficial_speed_m_s: 0.1,
            },
            "synthetic adaptive capillary fixture; no measured validation",
        )
        .unwrap(),
    }
}
fn nonlinear() -> Nonlinear {
    Nonlinear {
        absolute_saturation_residual_tolerance: 2e-13,
        relative_saturation_residual_tolerance: 0.0,
        volume_balance_tolerance_m: 1e-14,
        maximum_iterations: 128,
        maximum_line_search_trials: 48,
    }
}
fn accuracy(tolerance: f64) -> Accuracy {
    Accuracy {
        minimum_step_s: 1e-8,
        maximum_step_s: 0.1,
        maximum_evaluations: 100_000,
        absolute_saturation_tolerance: tolerance,
        relative_saturation_tolerance: 0.0,
        volume_balance_tolerance_m: 2e-11,
    }
}
fn input(boundaries: [Boundary; 2], duration_s: f64) -> Input {
    Input {
        duration_s,
        temperature_k: 300.0,
        gas_pressure_pa: 150_000.0,
        dynamic_viscosity_pa_s: 1e-3,
        body_force_pressure_gradient_pa_m: 0.0,
        boundaries,
    }
}
fn pressure(saturation: f64) -> f64 {
    let effective = (saturation - 0.1) / 0.8;
    150_000.0 - 10_000.0 * (1.0 / effective.powi(2) - 1.0).sqrt()
}
fn head(saturation: f64) -> Boundary {
    Boundary::Pressure {
        liquid_pressure_pa: pressure(saturation),
        surface_resistance_pa_s_m: 0.0,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
// Independent m=1/2 Darcy/retention formulas and RK4. Does not call the native
// constitutive evaluator, face operator, tangent or time integration.
fn derivative(saturation: f64, imposed: f64) -> f64 {
    let effective = (saturation - 0.1) / 0.8;
    let integral = 1.0 - (1.0 - effective * effective).sqrt();
    let mobility = 2e-10 * effective.sqrt() * integral * integral / 1e-3;
    2.0 * mobility * (imposed - pressure(saturation)) / (0.005 * 0.004)
}
fn reference(duration: f64, imposed: f64) -> f64 {
    let h = duration / 4096.0;
    let mut s = 0.5;
    for _ in 0..4096 {
        let a = derivative(s, imposed);
        let b = derivative(s + 0.5 * h * a, imposed);
        let c = derivative(s + 0.5 * h * b, imposed);
        let d = derivative(s + h * c, imposed);
        s += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    s
}

#[test]
fn controlled_wetting_refines_toward_independent_nonlinear_ode() {
    let original = Column::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let saved = original.clone();
    let forcing = input([head(0.6); 2], 0.1);
    let truth = reference(0.1, pressure(0.6));
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for tolerance in [2e-5, 2e-6, 2e-7] {
        let candidate = original
            .advance_candidate_implicit_with_boundaries(
                forcing,
                nonlinear(),
                accuracy(tolerance),
                || false,
            )
            .unwrap();
        let report = candidate.report();
        assert_eq!(report.interval_s, [0.0, 0.1]);
        assert_eq!(candidate.state().elapsed_s(), 0.1);
        assert!(report.rejected_intervals > 0);
        assert!(report.maximum_error_ratio <= 1.0);
        assert!(report.balance_residual_m.abs() <= accuracy(tolerance).volume_balance_tolerance_m);
        close(
            0.004 * (candidate.state().saturations()[0] - 0.5),
            report.boundary_volume_per_area_m,
            2e-11,
        );
        close(
            report.face_volume_per_area_m[0],
            -report.face_volume_per_area_m[1],
            1e-13,
        );
        errors.push((candidate.state().saturations()[0] - truth).abs());
        counts.push(report.accepted_intervals);
        let replay = original
            .advance_candidate_implicit_with_boundaries(
                forcing,
                nonlinear(),
                accuracy(tolerance),
                || false,
            )
            .unwrap();
        assert_eq!(replay, candidate);
    }
    assert!(errors[2] < errors[1] && errors[1] < errors[0], "{errors:?}");
    assert!(errors[2] < 1.5e-4, "{errors:?}");
    assert!(counts[2] > counts[1] && counts[1] > counts[0], "{counts:?}");
    assert_eq!(original, saved);
}

#[test]
fn closed_heterogeneous_column_preserves_inventory_and_resolves_relaxation() {
    let original = Column::new(
        vec![cell(0.01, 0.4, 2e-10), cell(0.03, 0.6, 8e-11)],
        vec![0.5, 0.65],
    )
    .unwrap();
    let candidate = original
        .advance_candidate_implicit_with_boundaries(
            input([Boundary::Sealed; 2], 1.0),
            nonlinear(),
            accuracy(1e-6),
            || false,
        )
        .unwrap();
    let final_s = candidate.state().saturations();
    assert!(final_s[0] > 0.5 && final_s[1] < 0.65);
    assert!(final_s[0] < final_s[1]);
    let report = candidate.report();
    assert_eq!(report.face_volume_per_area_m[0], 0.0);
    assert_eq!(report.face_volume_per_area_m[2], 0.0);
    assert_eq!(report.boundary_volume_per_area_m, 0.0);
    close(
        0.004 * (final_s[0] - 0.5),
        -report.face_volume_per_area_m[1],
        2e-11,
    );
    close(
        0.018 * (final_s[1] - 0.65),
        report.face_volume_per_area_m[1],
        2e-11,
    );
    close(
        candidate.state().liquid_volume_per_area_m().unwrap(),
        original.liquid_volume_per_area_m().unwrap(),
        2e-11,
    );
}

#[test]
fn prescribed_flux_uses_only_accepted_time_integrated_water() {
    let original = Column::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let q = 1e-5;
    let forcing = input(
        [
            Boundary::Flux {
                superficial_flux_m_s: q,
            },
            Boundary::Sealed,
        ],
        0.73,
    );
    let candidate = original
        .advance_candidate_implicit_with_boundaries(forcing, nonlinear(), accuracy(1e-8), || false)
        .unwrap();
    close(
        candidate.state().saturations()[0],
        0.5 + q * 0.73 / 0.004,
        2e-12,
    );
    close(
        candidate.report().boundary_volume_per_area_m,
        q * 0.73,
        1e-19,
    );
    close(
        candidate.report().face_volume_per_area_m[0],
        q * 0.73,
        1e-19,
    );
    assert_eq!(candidate.report().face_volume_per_area_m[1], 0.0);
    assert_eq!(candidate.state().elapsed_s(), 0.73);
}

#[test]
fn all_refusals_leave_original_column_and_clock_unchanged() {
    let original = Column::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let saved = original.clone();
    let forcing = input([head(0.6); 2], 0.1);
    let mut bounded = accuracy(1e-9);
    bounded.maximum_evaluations = 3;
    assert!(matches!(
        original
            .advance_candidate_implicit_with_boundaries(forcing, nonlinear(), bounded, || false),
        Err(Error::WorkBudget { evaluations: 3 })
    ));
    bounded = accuracy(1e-9);
    bounded.minimum_step_s = bounded.maximum_step_s;
    assert!(matches!(
        original
            .advance_candidate_implicit_with_boundaries(forcing, nonlinear(), bounded, || false),
        Err(Error::Accuracy { .. })
    ));
    let mut invalid_source = forcing;
    invalid_source.temperature_k = 100.0;
    assert!(matches!(
        original.advance_candidate_implicit_with_boundaries(
            invalid_source,
            nonlinear(),
            accuracy(1e-6),
            || false
        ),
        Err(Error::Step {
            error: LiquidSaturationImplicitError::Material(_),
            ..
        })
    ));
    let mut calls = 0;
    assert!(matches!(
        original.advance_candidate_implicit_with_boundaries(
            forcing,
            nonlinear(),
            accuracy(1e-6),
            || {
                calls += 1;
                calls > 40
            }
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(original, saved);
}

#[test]
fn equilibrium_and_exact_seals_preserve_state_without_hidden_mobility() {
    let original = Column::new(
        vec![cell(0.01, 0.4, 0.0), cell(0.03, 0.6, 8e-11)],
        vec![0.6, 0.5],
    )
    .unwrap();
    let candidate = original
        .advance_candidate_implicit_with_boundaries(
            input([Boundary::Sealed; 2], 3.0),
            nonlinear(),
            accuracy(1e-10),
            || false,
        )
        .unwrap();
    assert_eq!(candidate.state().saturations(), original.saturations());
    assert_eq!(candidate.report().face_volume_per_area_m, vec![0.0; 3]);
    assert_eq!(candidate.report().linear_solves, 0);
    assert_eq!(candidate.report().maximum_error_ratio, 0.0);
}

#[test]
fn whole_window_water_budget_is_not_multiplied_by_step_count() {
    let original = Column::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let saved = original.clone();
    let forcing = input(
        [
            Boundary::Flux {
                superficial_flux_m_s: 1e-5,
            },
            Boundary::Sealed,
        ],
        1.0,
    );
    // Each deliberately loose endpoint solve can accept its original state.
    // The complete window must still reject the accumulated missing uptake.
    let mut loose = nonlinear();
    loose.absolute_saturation_residual_tolerance = 1e-3;
    loose.volume_balance_tolerance_m = 1e-5;
    let mut controlled = accuracy(1e-6);
    controlled.volume_balance_tolerance_m = 2e-6;
    let error = original
        .advance_candidate_implicit_with_boundaries(forcing, loose, controlled, || false)
        .unwrap_err();
    match error {
        Error::Conservation {
            residual_m,
            maximum_m,
        } => {
            close(residual_m, -1e-5, 1e-18);
            assert_eq!(maximum_m, 2e-6);
        }
        other => panic!("expected fixed whole-window conservation refusal: {other:?}"),
    }
    assert_eq!(original, saved);
}

#[test]
fn unequal_representable_halves_retain_the_correct_local_error_estimate() {
    let initial = Column::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let original = initial
        .step_candidate_implicit_with_boundaries(
            input([Boundary::Sealed; 2], 1e16),
            nonlinear(),
            || false,
        )
        .unwrap()
        .0
        .into_parts()
        .0;
    // At this clock, midpoint rounding produces four and two second substeps.
    let forcing = input([head(0.5001); 2], 6.0);
    let mut controlled = accuracy(1e-2);
    controlled.minimum_step_s = 6.0;
    controlled.maximum_step_s = 6.0;
    let candidate = original
        .advance_candidate_implicit_with_boundaries(forcing, nonlinear(), controlled, || false)
        .unwrap();
    let coarse = original
        .step_candidate_implicit_with_boundaries(forcing, nonlinear(), || false)
        .unwrap()
        .0;
    let half = original
        .step_candidate_implicit_with_boundaries(
            Input {
                duration_s: 4.0,
                ..forcing
            },
            nonlinear(),
            || false,
        )
        .unwrap()
        .0;
    let fine = half
        .state()
        .step_candidate_implicit_with_boundaries(
            Input {
                duration_s: 2.0,
                ..forcing
            },
            nonlinear(),
            || false,
        )
        .unwrap()
        .0;
    assert_eq!(candidate.state(), fine.state());
    assert_eq!(candidate.report().accepted_intervals, 1);
    assert_eq!(candidate.report().minimum_accepted_substep_s, 2.0);
    assert_eq!(candidate.report().maximum_accepted_substep_s, 4.0);
    close(
        candidate.report().maximum_error_ratio,
        1.25 * (fine.state().saturations()[0] - coarse.state().saturations()[0]).abs() / 1e-2,
        1e-15,
    );
}
