//! Public unsaturated-liquid endpoints against independently evaluated Darcy
//! balances and a nonlinear one-cell ODE; all calibrations here are synthetic.

use fs_material::porous::wet::liquid::{
    LiquidCapillaryPressure, WetPoreLiquidLaw, WetPoreLiquidSupport,
    transient::{
        LiquidSaturationCell, LiquidSaturationColumn, LiquidSaturationImplicitControls,
        LiquidSaturationImplicitError, LiquidSaturationStepInput,
    },
};

fn law_with_support(saturation: (f64, f64), maximum_speed: f64) -> WetPoreLiquidLaw {
    WetPoreLiquidLaw::try_new(
        0.1,
        0.1,
        0.5,
        0.5,
        10_000.0,
        WetPoreLiquidSupport {
            liquid_saturation: saturation,
            temperature_k: (280.0, 330.0),
            liquid_pressure_pa: (50_000.0, 200_000.0),
            maximum_superficial_speed_m_s: maximum_speed,
        },
        "synthetic implicit Richards fixture; no measured calibration claim",
    )
    .unwrap()
}

fn cell(thickness: f64, porosity: f64, permeability: f64) -> LiquidSaturationCell {
    LiquidSaturationCell {
        thickness_m: thickness,
        porosity,
        intrinsic_permeability_m2: permeability,
        law: law_with_support((0.1, 0.9), 0.1),
    }
}

fn input(boundary: [f64; 2], duration: f64) -> LiquidSaturationStepInput {
    LiquidSaturationStepInput {
        duration_s: duration,
        temperature_k: 300.0,
        gas_pressure_pa: 150_000.0,
        dynamic_viscosity_pa_s: 1e-3,
        body_force_pressure_gradient_pa_m: 0.0,
        boundary_liquid_pressure_pa: boundary,
    }
}

fn controls() -> LiquidSaturationImplicitControls {
    LiquidSaturationImplicitControls {
        absolute_saturation_residual_tolerance: 1e-11,
        relative_saturation_residual_tolerance: 0.0,
        volume_balance_tolerance_m: 1e-13,
        maximum_iterations: 128,
        maximum_line_search_trials: 48,
    }
}

// Independent m=1/2 closed form, using neither the native samples nor its
// logarithmic evaluation, derivative, face-resistance or nonlinear solvers.
fn pressure(saturation: f64) -> f64 {
    let effective = (saturation - 0.1) / 0.8;
    150_000.0 - 10_000.0 * (1.0 / (effective * effective) - 1.0).sqrt()
}

fn mobility(saturation: f64, permeability: f64) -> f64 {
    let effective = (saturation - 0.1) / 0.8;
    let integral = 1.0 - (1.0 - effective * effective).sqrt();
    permeability * effective.sqrt() * integral * integral / 1e-3
}

fn native_pressure(state: &LiquidSaturationColumn, index: usize) -> f64 {
    match state.states()[index].capillary_pressure() {
        LiquidCapillaryPressure::Finite(pc) => 150_000.0 - pc,
        LiquidCapillaryPressure::DryLimit => panic!("test requires a mobile cell"),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}",
    );
}

#[test]
fn heterogeneous_implicit_endpoint_matches_manufactured_darcy_balance() {
    let cells = vec![cell(0.01, 0.4, 2e-10), cell(0.03, 0.6, 8e-11)];
    let original = LiquidSaturationColumn::new(cells, vec![0.5, 0.56]).unwrap();
    let saved = original.clone();
    let target = [0.52, 0.54];
    let dt = 2.0;
    let body = 150.0;
    let half = [0.005, 0.015];
    let capacity = [0.004, 0.018];
    let p = target.map(pressure);
    let m = [mobility(target[0], 2e-10), mobility(target[1], 8e-11)];
    let middle_flux =
        (p[0] - p[1] + body * (half[0] + half[1])) / (half[0] / m[0] + half[1] / m[1]);
    let left_flux = middle_flux + capacity[0] * (target[0] - 0.5) / dt;
    let right_flux = middle_flux - capacity[1] * (target[1] - 0.56) / dt;
    let boundary = [
        p[0] + half[0] * (left_flux / m[0] - body),
        p[1] - half[1] * (right_flux / m[1] - body),
    ];
    let mut forcing = input(boundary, dt);
    forcing.body_force_pressure_gradient_pa_m = body;
    let (candidate, report) = original
        .step_candidate_implicit(forcing, controls(), || false)
        .unwrap();
    for i in 0..2 {
        close(candidate.state().saturations()[i], target[i], 2e-10);
        assert_eq!(
            candidate.state().states()[i].law(),
            original.states()[i].law()
        );
        let independently_integrated =
            dt * (candidate.step().face_flux_m_s[i] - candidate.step().face_flux_m_s[i + 1]);
        close(
            capacity[i] * (candidate.state().saturations()[i] - original.saturations()[i]),
            independently_integrated,
            controls().absolute_saturation_residual_tolerance * capacity[i],
        );
    }
    for (actual, expected) in
        candidate
            .step()
            .face_flux_m_s
            .iter()
            .zip([left_flux, middle_flux, right_flux])
    {
        close(*actual, expected, 2e-12);
    }
    assert!(report.saturation_residual_ratio <= 1.0);
    assert!(report.volume_balance_residual_m.abs() <= controls().volume_balance_tolerance_m);
    close(candidate.state().elapsed_s(), dt, 0.0);
    assert_eq!(original, saved);

    // The explicit owner must use the same actual half-cell series resistance,
    // not an unweighted harmonic average appropriate only to equal widths.
    forcing.duration_s = 1e-5;
    let explicit = original.step_candidate(forcing, || false).unwrap();
    let expected = (pressure(0.5) - pressure(0.56) + body * (half[0] + half[1]))
        / (half[0] / mobility(0.5, 2e-10) + half[1] / mobility(0.56, 8e-11));
    close(explicit.step().face_flux_m_s[1], expected, 2e-16);
}

#[test]
fn long_implicit_wetting_survives_explicit_limit_and_checks_full_endpoint_law() {
    let original = LiquidSaturationColumn::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let forcing = input([pressure(0.6); 2], 20.0);
    assert!(original.step_candidate(forcing, || false).is_err());
    let mut budget = controls();
    // A loose local stop must not bypass the independent fixed column budget.
    budget.absolute_saturation_residual_tolerance = 1e-3;
    let (candidate, report) = original
        .step_candidate_implicit(forcing, budget, || false)
        .unwrap();
    let s = candidate.state().saturations()[0];
    assert!(s > 0.5 && s < 0.6);
    assert!(report.frozen_mobility_iterations > 0);
    assert_eq!(report.linear_solves, report.iterations);
    assert!(report.state_evaluations >= report.iterations + 1);
    let q = mobility(s, 2e-10) * (forcing.boundary_liquid_pressure_pa[0] - pressure(s)) / 0.005;
    close(candidate.step().face_flux_m_s[0], q, 2e-16);
    close(candidate.step().face_flux_m_s[1], -q, 2e-16);
    close(0.004 * (s - 0.5), 40.0 * q, 1.1e-13);
    assert!(report.volume_balance_residual_m.abs() <= budget.volume_balance_tolerance_m);
    assert_eq!(original.saturations(), &[0.5]);
    assert_eq!(original.elapsed_s(), 0.0);
}

fn independent_one_cell_ode(saturation: f64, boundary: f64) -> f64 {
    4.0 * mobility(saturation, 2e-10) * (boundary - pressure(saturation)) / (0.4 * 0.01 * 0.01)
}

fn reference_rk4(duration: f64, boundary: f64) -> f64 {
    let h = duration / 4096.0;
    let mut s = 0.5;
    for _ in 0..4096 {
        let a = independent_one_cell_ode(s, boundary);
        let b = independent_one_cell_ode(s + 0.5 * h * a, boundary);
        let c = independent_one_cell_ode(s + 0.5 * h * b, boundary);
        let d = independent_one_cell_ode(s + h * c, boundary);
        s += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    s
}

#[test]
fn backward_euler_refines_toward_independent_nonlinear_transient() {
    let duration = 0.1;
    let boundary = pressure(0.6);
    let reference = reference_rk4(duration, boundary);
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut state =
            LiquidSaturationColumn::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
        for _ in 0..count {
            let (candidate, report) = state
                .step_candidate_implicit(
                    input([boundary; 2], duration / f64::from(count)),
                    controls(),
                    || false,
                )
                .unwrap();
            assert!(report.saturation_residual_ratio <= 1.0);
            state = candidate.into_parts().0;
        }
        close(state.elapsed_s(), duration, 8.0 * f64::EPSILON * duration);
        errors.push((state.saturations()[0] - reference).abs());
    }
    assert!(errors[2] < errors[1] && errors[1] < errors[0], "{errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(ratio > 1.6 && ratio < 2.4, "first-order errors: {errors:?}");
    }
}

#[test]
fn uniform_equilibrium_and_exact_dry_or_intrinsic_seals_need_no_nonlinear_work() {
    let wet = LiquidSaturationColumn::new(
        vec![cell(0.01, 0.4, 2e-10), cell(0.03, 0.6, 8e-11)],
        vec![0.5; 2],
    )
    .unwrap();
    let p = native_pressure(&wet, 0);
    let (candidate, report) = wet
        .step_candidate_implicit(input([p; 2], 1e8), controls(), || false)
        .unwrap();
    assert_eq!(candidate.state().saturations(), wet.saturations());
    assert_eq!(candidate.step().face_flux_m_s, vec![0.0; 3]);
    assert_eq!(report.iterations, 0);
    assert_eq!(report.linear_solves, 0);
    assert_eq!(report.state_evaluations, 1);
    for (s, permeability) in [(0.1, 2e-10), (0.6, 0.0)] {
        let sealed = LiquidSaturationColumn::new(
            vec![cell(0.01, 0.4, permeability), cell(0.02, 0.5, 2e-10)],
            vec![s, 0.5],
        )
        .unwrap();
        let p = native_pressure(&sealed, 1);
        let (candidate, report) = sealed
            .step_candidate_implicit(input([150_000.0, p], 1e8), controls(), || false)
            .unwrap();
        assert_eq!(candidate.state().saturations(), sealed.saturations());
        assert_eq!(candidate.step().face_flux_m_s, vec![0.0; 3]);
        assert_eq!(report.iterations, 0);
        assert_eq!(report.volume_balance_residual_m, 0.0);
        if s == 0.1 {
            assert_eq!(candidate.step().liquid_pressure_pa_after[0], None);
        }
    }

    // Exercise assembly, not only the equilibrium fast path: the admitted
    // sealed cell's dt/capacity overflows, but every connected tangent is zero.
    let thin_seal = LiquidSaturationColumn::new(
        vec![cell(1e-308, 0.4, 0.0), cell(0.01, 0.4, 2e-10)],
        vec![0.6, 0.5],
    )
    .unwrap();
    let (candidate, report) = thin_seal
        .step_candidate_implicit(input([150_000.0, pressure(0.6)], 1.0), controls(), || false)
        .unwrap();
    assert!(report.iterations > 0);
    assert_eq!(candidate.state().saturations()[0], 0.6);
    assert!(candidate.state().saturations()[1] > 0.5);
    assert_eq!(&candidate.step().face_flux_m_s[..2], &[0.0, 0.0]);
    assert!(report.volume_balance_residual_m.abs() <= controls().volume_balance_tolerance_m);
}

#[test]
fn source_refusals_remain_distinct_from_bounded_trial_refusal() {
    let original = LiquidSaturationColumn::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let saved = original.clone();
    let mut forcing = input([pressure(0.6); 2], 1.0);
    forcing.temperature_k = 340.0;
    assert!(matches!(
        original.step_candidate_implicit(forcing, controls(), || false),
        Err(LiquidSaturationImplicitError::Material(_))
    ));
    forcing.temperature_k = 300.0;
    forcing.gas_pressure_pa = 60_000.0; // Actual liquid pressure is below support.
    assert!(matches!(
        original.step_candidate_implicit(forcing, controls(), || false),
        Err(LiquidSaturationImplicitError::Material(_))
    ));

    let mut speed_cell = cell(0.01, 0.4, 2e-10);
    speed_cell.law = law_with_support((0.1, 0.9), 1e-6);
    let speed_limited = LiquidSaturationColumn::new(vec![speed_cell], vec![0.5]).unwrap();
    assert!(matches!(
        speed_limited.step_candidate_implicit(input([pressure(0.6); 2], 1.0), controls(), || false),
        Err(LiquidSaturationImplicitError::Material(_))
    ));

    let mut narrow_cell = cell(0.01, 0.4, 2e-10);
    narrow_cell.law = law_with_support((0.1, 0.51), 0.1);
    let narrow = LiquidSaturationColumn::new(vec![narrow_cell], vec![0.5]).unwrap();
    let mut one_trial = controls();
    one_trial.maximum_line_search_trials = 1;
    let result = narrow.step_candidate_implicit(input([150_000.0; 2], 1.0), one_trial, || false);
    match result {
        Err(LiquidSaturationImplicitError::NoAdmissibleTrial { report, .. }) => {
            assert_eq!(report.state_evaluations, 2);
            assert_eq!(report.iterations, 0);
            assert_eq!(report.linear_solves, 1);
        }
        other => panic!("expected a distinct rejected source trial: {other:?}"),
    }
    assert_eq!(narrow.saturations(), &[0.5]);
    assert_eq!(narrow.elapsed_s(), 0.0);
    assert_eq!(original, saved);
}

#[test]
fn nonlinear_work_and_late_cancellation_never_publish_a_partial_column() {
    let original = LiquidSaturationColumn::new(vec![cell(0.01, 0.4, 2e-10)], vec![0.5]).unwrap();
    let saved = original.clone();
    let forcing = input([pressure(0.6); 2], 20.0);
    let mut short = controls();
    short.maximum_iterations = 1;
    match original.step_candidate_implicit(forcing, short, || false) {
        Err(LiquidSaturationImplicitError::NotConverged { report, .. }) => {
            assert_eq!(report.iterations, 1);
            assert_eq!(report.linear_solves, 1);
            assert!(report.saturation_residual_ratio > 1.0);
        }
        other => panic!("expected bounded nonlinear refusal: {other:?}"),
    }
    let mut polls = 0;
    let cancelled = original.step_candidate_implicit(forcing, controls(), || {
        polls += 1;
        polls > 20
    });
    assert!(matches!(
        cancelled,
        Err(LiquidSaturationImplicitError::Cancelled)
    ));
    assert!(polls > 20);
    assert_eq!(original, saved);
}
