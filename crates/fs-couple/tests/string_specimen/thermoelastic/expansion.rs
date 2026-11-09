//! Actual source-table expansion through the live thermal/string transaction.
//! These manufactured coefficients test the model, not a particular wire.

use super::*;

const SOURCE: [(f64, f64); 4] = [
    (250.0, 1.5e-5),
    (302.0, 2.0e-5),
    (306.0, 2.04e-5),
    (550.0, 4.0e-5),
];

fn table(source: &[(f64, f64)]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: fs_qty::Temperature::DIMS,
        knots: source.to_vec(),
        dims: INVERSE_TEMPERATURE_DIMS,
    }
}

fn with_source(
    cx: &fs_exec::Cx<'_>,
    source: &[(f64, f64)],
    maximum_delta_t: f64,
) -> ThermalMaterialStringRuntime {
    let (mechanical, card, curve) = parts_with_expansion(cx, false, true, table(source));
    ThermalMaterialStringRuntime::try_new_thermoelastic(
        cx, mechanical, card, curve, 0.0, maximum_delta_t, controls(),
    ).unwrap()
}

// Independent analytic integral of the manufactured table, with the same
// declared 300 K stress-free datum. No production expansion evaluator is called.
fn alpha_and_free(source: &[(f64, f64)], temperature: f64) -> (f64, f64) {
    let (from, to, sign) = if temperature >= 300.0 {
        (300.0, temperature, 1.0)
    } else {
        (temperature, 300.0, -1.0)
    };
    let mut free = 0.0;
    let mut alpha = None;
    for pair in source.windows(2) {
        let [(lo, a), (hi, b)] = [pair[0], pair[1]];
        let slope = (b-a)/(hi-lo);
        let lower = from.max(lo);
        let upper = to.min(hi);
        if upper > lower {
            let delta = upper-lower;
            free += (a+slope*(lower-lo))*delta+0.5*slope*delta*delta;
        }
        if lo <= temperature && temperature <= hi {
            alpha = Some(a+slope*(temperature-lo));
        }
    }
    (alpha.unwrap(), sign*free)
}

fn near(actual: f64, expected: f64) {
    let scale = actual.abs().max(expected.abs()).max(1e-10);
    assert!((actual-expected).abs() <= 2e-12*scale, "{actual} != {expected}");
}

fn internal_energy(runtime: &ThermalMaterialStringRuntime, source: &[(f64, f64)]) -> f64 {
    let specimen = runtime.mechanical().specimen();
    let string = specimen.string();
    let temperature = runtime.thermal().temperature_k();
    let (alpha, free) = alpha_and_free(source, temperature);
    let young = 2e9+YOUNG_SLOPE*(temperature-300.0);
    let elastic = (string.length_m-REFERENCE_LENGTH)/REFERENCE_LENGTH-free;
    let volume = specimen.area_m2()*REFERENCE_LENGTH;
    let axial = 0.5*volume*young*elastic*elastic;
    let mut entropy = volume*(young*alpha*elastic-0.5*YOUNG_SLOPE*elastic*elastic);
    for (index, state) in runtime.mechanical().states().iter().enumerate() {
        let k = (index+1) as f64*core::f64::consts::PI/string.length_m;
        let derivative = (YOUNG_SLOPE*specimen.second_moment_m4()*k.powi(4)
            +specimen.area_m2()*(YOUNG_SLOPE*elastic-young*alpha)*k*k)
            / string.lin_density_kg_m;
        entropy -= 0.5*derivative*state.displacement_m_sqrt_kg.powi(2);
    }
    modal_energy(runtime)+axial+temperature*entropy
}

fn same_history(actual: &ThermalMaterialStringRuntime, before: &ThermalMaterialStringRuntime) {
    assert_eq!(actual.thermal(), before.thermal());
    assert_eq!(actual.mechanical().specimen(), before.mechanical().specimen());
    assert_eq!(actual.mechanical().states(), before.mechanical().states());
    assert_eq!(actual.mechanical().epoch(), before.mechanical().epoch());
    assert_eq!(actual.mechanical().accepted_samples(), before.mechanical().accepted_samples());
}

#[test]
fn tabulated_expansion_crosses_its_knot_in_live_pitch_pressure_and_first_law() {
    with_string_cx(|cx, _| {
        let mut runtime = with_source(cx, &SOURCE, 5.0);
        let original = runtime.mechanical().specimen().clone();
        let original_frequency = runtime.mechanical().modes()[0].angular_frequency_rad_s;
        let original_claim = original.thermal_expansion().unwrap().selected_claim();
        let mut crossed_knot = false;
        let mut pressure_difference = 0.0;
        let mut entropy_exchange = 0.0;
        for epoch in 0..64 {
            let before = internal_energy(&runtime, &SOURCE);
            let previous_temperature = runtime.thermal().temperature_k();
            let mut frozen = runtime.mechanical().clone();
            let frozen_frame = frozen.step(cx, &[0.0]).unwrap();
            let frame = runtime.step(cx, epoch, &[0.0], 0.002).unwrap();
            let temperature = runtime.thermal().temperature_k();
            crossed_knot |= previous_temperature < 302.0 && temperature >= 302.0;
            assert_eq!(runtime.mechanical().states(), frozen.states());
            let specimen = runtime.mechanical().specimen();
            let (alpha, free) = alpha_and_free(&SOURCE, temperature);
            let young = 2e9+YOUNG_SLOPE*(temperature-300.0);
            let elastic = (0.502-REFERENCE_LENGTH)/REFERENCE_LENGTH-free;
            near(specimen.string().tension_n, young*specimen.area_m2()*elastic);
            near(specimen.thermal_expansion().unwrap().free_linear_strain(), free);
            near(specimen.thermal_expansion().unwrap().current().resolved()
                .property(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY).unwrap().value_si(), alpha);
            assert_eq!(specimen.thermal_expansion().unwrap().selected_claim(), original_claim);
            assert_eq!(specimen.mass_kg(), original.mass_kg());
            assert_eq!(specimen.area_m2(), original.area_m2());
            assert_eq!(specimen.string().lin_density_kg_m, original.string().lin_density_kg_m);
            near(runtime.mechanical().modes()[0].angular_frequency_rad_s,
                string_mode_omega(&specimen.string(), 1));
            let independent = frame.vibration.acoustic.input_work_j+frame.external_heat_j
                -(internal_energy(&runtime, &SOURCE)-before+frame.enthalpy_change_j
                    +frame.vibration.dissipation.air_loss_j
                    +frame.vibration.dissipation.authored_loss_j
                    +frame.vibration.dissipation.roundoff_residual_j);
            let allowance = frame.energy_roundoff_tolerance_j+frame.thermal_solve_tolerance_j;
            assert!(independent.abs() < allowance+8e-13, "independent first law: {independent}");
            entropy_exchange += frame.mechanical_entropy_energy_change_j.abs();
            pressure_difference += (frame.vibration.acoustic.observer_pressure_pa
                -frozen_frame.acoustic.observer_pressure_pa).abs();
        }
        assert!(crossed_knot, "a real accepted thermal interval must cross the CTE slope change");
        assert!(runtime.mechanical().modes()[0].angular_frequency_rad_s < 0.99*original_frequency);
        assert!(entropy_exchange > 1e-5);
        assert!(pressure_difference > 1e-9);
    });
}

#[test]
fn constant_table_preserves_scalar_expansion_motion_and_heat() {
    with_string_cx(|cx, _| {
        let mut scalar = runtime(cx, 5.0);
        let mut tabulated = with_source(cx, &[(250.0, ALPHA), (300.0, ALPHA), (550.0, ALPHA)], 5.0);
        for epoch in 0..16 {
            let a = scalar.step(cx, epoch, &[0.0], 0.002).unwrap();
            let b = tabulated.step(cx, epoch, &[0.0], 0.002).unwrap();
            near(a.reversible_heat_j, b.reversible_heat_j);
            near(a.vibration.acoustic.observer_pressure_pa, b.vibration.acoustic.observer_pressure_pa);
            near(scalar.thermal().temperature_k(), tabulated.thermal().temperature_k());
            near(scalar.mechanical().specimen().string().tension_n,
                tabulated.mechanical().specimen().string().tension_n);
            for (x, y) in scalar.mechanical().states().iter().zip(tabulated.mechanical().states()) {
                near(x.displacement_m_sqrt_kg, y.displacement_m_sqrt_kg);
                near(x.velocity_m_sqrt_kg_per_s, y.velocity_m_sqrt_kg_per_s);
            }
        }
    });
}

#[test]
fn negative_tabulated_expansion_contracts_the_natural_wire_and_raises_pitch() {
    with_string_cx(|cx, _| {
        let source = [(250.0, -0.5e-5), (302.0, -1.0e-5), (550.0, -3.0e-5)];
        let mut runtime = with_source(cx, &source, 5.0);
        let original_frequency = runtime.mechanical().modes()[0].angular_frequency_rad_s;
        let original_tension = runtime.mechanical().specimen().string().tension_n;
        for epoch in 0..64 {
            let before = internal_energy(&runtime, &source);
            let frame = runtime.step(cx, epoch, &[0.0], 0.002).unwrap();
            let independent = frame.vibration.acoustic.input_work_j+frame.external_heat_j
                -(internal_energy(&runtime, &source)-before+frame.enthalpy_change_j
                    +frame.vibration.dissipation.air_loss_j
                    +frame.vibration.dissipation.authored_loss_j
                    +frame.vibration.dissipation.roundoff_residual_j);
            assert!(independent.abs()
                < frame.energy_roundoff_tolerance_j+frame.thermal_solve_tolerance_j+8e-13);
        }
        assert!(runtime.thermal().temperature_k() > 303.0);
        assert!(runtime.mechanical().specimen().thermal_expansion().unwrap().free_linear_strain() < 0.0);
        assert!(runtime.mechanical().specimen().string().tension_n > original_tension);
        assert!(runtime.mechanical().modes()[0].angular_frequency_rad_s > original_frequency);
    });
}

#[test]
fn table_slope_capacity_and_original_source_domain_are_admitted_before_motion() {
    with_string_cx(|cx, _| {
        // Alpha vanishes at launch. A narrow later interval has negative alpha'
        // and positive prestress, hence a large POSITIVE Psi_TT. Checking only
        // launch alpha, or dropping E*alpha', misses the negative total capacity.
        let sharp = [(250.0, 0.0), (350.0, 0.0), (350.001, -1e-6), (550.0, -1e-6)];
        let (mechanical, card, curve) = parts_with_expansion(cx, false, true, table(&sharp));
        let result = ThermalMaterialStringRuntime::try_new_thermoelastic(
            cx, mechanical, card, curve, 0.0, 5.0, controls());
        assert!(matches!(result, Err(ThermalStringError::Admission(message))
            if message.contains("capacity")));

        // Supported launch/reference points do not authorize extrapolation
        // over the caloric interval used by the coupled capacity admission.
        let narrow = [(250.0, 1e-5), (500.0, 2e-5)];
        let (mechanical, card, curve) = parts_with_expansion(cx, false, true, table(&narrow));
        let result = ThermalMaterialStringRuntime::try_new_thermoelastic(
            cx, mechanical, card, curve, 0.0, 5.0, controls());
        assert!(matches!(result, Err(ThermalStringError::Admission(message))
            if message.contains("continuously cover")));
    });
}

#[test]
fn tabulated_expansion_keeps_electrical_heat_on_the_same_native_transaction() {
    with_string_cx(|cx, _| {
        let mut powered = with_source(cx, &SOURCE, 5.0);
        let mut direct_heat = powered.clone();
        let mut proposals = 0;
        for epoch in 0..48 {
            let frame = powered.step_with_current_and_thermal_transport(cx, epoch, &[0.0], 15.0,
                |input| {
                    proposals += 1;
                    Ok::<_, fs_material::phase::PhaseStateError>(UniformEnthalpyStep {
                        state: input.curve.advance_specific_energy(input.initial,
                            input.internal_heat_j/input.mass_kg)?,
                        external_heat_j: 0.0,
                        energy_residual_tolerance_j: 0.0,
                        report: (),
                    })
                }).unwrap();
            let reference = direct_heat.step(cx, epoch, &[0.0], frame.joule_heat_j).unwrap();
            assert_eq!(frame.coupled, reference);
            assert_eq!(powered.mechanical().states(), direct_heat.mechanical().states());
            assert_eq!(frame.conductor_after.material().query_point(),
                powered.mechanical().specimen().material().query_point());
        }
        assert!(proposals > 48);
        assert!(powered.thermal().temperature_k() > 302.0);
    });
}

#[test]
fn slack_domain_exit_and_late_cancellation_preserve_tabulated_history_and_retry() {
    with_string_cx(|cx, _| {
        let mut runtime = with_source(cx, &SOURCE, 300.0);
        let mut reference = runtime.clone();
        assert!(runtime.step(cx, 0, &[0.0], 8.0).is_err());
        same_history(&runtime, &reference);
        assert!(runtime.step(cx, 0, &[0.0], 100.0).is_err());
        same_history(&runtime, &reference);
        with_string_cx(|cancelled, gate| {
            let mut proposed = false;
            assert!(runtime.step_with_thermal_transport(cancelled, 0, &[0.0], |input| {
                let state = input.curve.advance_specific_energy(input.initial,
                    (input.internal_heat_j+0.002)/input.mass_kg)?;
                proposed = true;
                gate.request();
                Ok::<_, fs_material::phase::PhaseStateError>(UniformEnthalpyStep {
                    state,
                    external_heat_j: 0.002,
                    energy_residual_tolerance_j: 0.0,
                    report: (),
                })
            }).is_err());
            assert!(proposed, "cancellation must follow an original caloric proposal");
        });
        same_history(&runtime, &reference);
        assert_eq!(runtime.step(cx, 0, &[0.0], 0.002).unwrap(),
            reference.step(cx, 0, &[0.0], 0.002).unwrap());
    });
}
