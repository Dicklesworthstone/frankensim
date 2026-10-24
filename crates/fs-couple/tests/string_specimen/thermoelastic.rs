//! Live source-driven thermoelastic pitch and its complete energy account.
//! The deliberately simple coefficients are synthetic, not a measured wire.

use super::*;
use fs_couple::acoustic_realize::thermoelastic::ThermoelasticStringSolve;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStep};

const RATE: u32 = 8_000;
const ALPHA: f64 = 2e-5;
const YOUNG_SLOPE: f64 = -2e6;
const REFERENCE_LENGTH: f64 = 0.5;

fn parts(
    cx: &fs_exec::Cx<'_>,
    nonaffine: bool,
    sourced_expansion: bool,
) -> (LinearMaterialStringRuntime, MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, PropertyValue::Scalar { value: 1_000.0, dims: Density::DIMS }),
        ("young_modulus", Pressure::DIMS, PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(250.0, 2.1e9), (400.0, if nonaffine { 1.7e9 } else { 1.8e9 }), (550.0, 1.5e9)],
            dims: Pressure::DIMS,
        }),
        (KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY, DynViscosity::DIMS,
            PropertyValue::Scalar { value: 1e7, dims: DynViscosity::DIMS }),
        (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS,
            PropertyValue::Scalar { value: ALPHA, dims: INVERSE_TEMPERATURE_DIMS }),
        (fs_material::conductor::ELECTRICAL_RESISTIVITY_PROPERTY,
            fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS, PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(250.0, 0.9e-6), (550.0, 1.5e-6)],
                dims: fs_material::conductor::ELECTRICAL_RESISTIVITY_DIMS,
            }),
    ] {
        let interpolation = if matches!(value, PropertyValue::Curve { .. }) {
            InterpolationPolicy::LinearInside
        } else { InterpolationPolicy::ConstantWithinValidity };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value,
            validity: ValidityDomain::unconstrained().with("T", 250.0, 550.0)
                .with("omega", 1.0, 20_000.0),
            uncertainty: UncertaintyModel::Unstated, interpolation,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic live thermoelastic string".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        if name != LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY {
            requirements.push(ScalarPropertyRequirement::try_new(
                name, dims, ScalarAdmissibility::NonNegative,
            ).unwrap());
        }
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic expanding conducting wire".into(),
        phase: "solid".into(), process: "synthetic".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap().with("omega", 1.0).unwrap();
    let material = resolve_material_state_point(&card, &point, &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let template = PrestressedString { length_m: 0.502, ..loss_template() };
    let geometry = StringGeometryConstraint::FixedRadius(0.0015);
    let specimen = if sourced_expansion {
        let expansion = integrate_isotropic_thermal_expansion(&card, &point, &point,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        with_uniform_circular_thermal_extension(template, geometry, &material, &expansion,
            REFERENCE_LENGTH, 0.01).unwrap()
    } else {
        with_uniform_circular_material_and_constraints(template, geometry, &material,
            StringPrestress::FixedThermalExtension {
                reference_stress_free_length_m: REFERENCE_LENGTH,
                free_thermal_strain: 0.0, linear_strain_limit: 0.01,
            }).unwrap()
    }.with_kelvin_voigt_bending_loss().unwrap();
    let mechanical = LinearMaterialStringRuntime::try_new(cx, specimen,
        Some(Pluck { station_frac: 0.4, height_m: 1e-3 }),
        incremental_ambient(), 1.0, RATE).unwrap();
    // Constant synthetic cp=10 J/(kg K); caloric mass retains loaded geometry.
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), vec![
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: -500.0, temperature_k: 250.0,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1_000.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 2_500.0, temperature_k: 550.0,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1_000.0 },
    ]).unwrap();
    (mechanical, card, curve)
}

fn controls() -> ThermoelasticStringSolve {
    ThermoelasticStringSolve {
        reversible_heat_bracket_j: [-0.004, 0.004],
        energy_tolerance_j: 1e-12,
        maximum_iterations: 96,
    }
}

fn runtime(cx: &fs_exec::Cx<'_>, maximum_delta_t: f64) -> ThermalMaterialStringRuntime {
    let (mechanical, card, curve) = parts(cx, false, true);
    ThermalMaterialStringRuntime::try_new_thermoelastic(cx, mechanical, card, curve,
        0.0, maximum_delta_t, controls()).unwrap()
}

fn modal_energy(runtime: &ThermalMaterialStringRuntime) -> f64 {
    runtime.mechanical().states().iter().zip(runtime.mechanical().modes())
        .map(|(state, mode)| 0.5 * (state.velocity_m_sqrt_kg_per_s.powi(2)
            + (mode.angular_frequency_rad_s * state.displacement_m_sqrt_kg).powi(2)))
        .sum()
}

// Independent free-energy derivatives, including the vibrating string's
// entropy. Merely adding axial stored energy misses the T*S contribution.
fn mechanical_internal_energy(runtime: &ThermalMaterialStringRuntime) -> f64 {
    let specimen = runtime.mechanical().specimen();
    let string = specimen.string();
    let t = runtime.thermal().temperature_k();
    let young = 2e9 + YOUNG_SLOPE * (t - 300.0);
    let strain = (string.length_m - REFERENCE_LENGTH) / REFERENCE_LENGTH
        - ALPHA * (t - 300.0);
    let volume = specimen.area_m2() * REFERENCE_LENGTH;
    let axial = 0.5 * volume * young * strain.powi(2);
    let mut entropy = volume * (young * ALPHA * strain - 0.5 * YOUNG_SLOPE * strain.powi(2));
    for (index, state) in runtime.mechanical().states().iter().enumerate() {
        let k = (index + 1) as f64 * core::f64::consts::PI / string.length_m;
        let stiffness_slope = (YOUNG_SLOPE * specimen.second_moment_m4() * k.powi(4)
            + specimen.area_m2() * (YOUNG_SLOPE * strain - young * ALPHA) * k.powi(2))
            / string.lin_density_kg_m;
        entropy -= 0.5 * stiffness_slope * state.displacement_m_sqrt_kg.powi(2);
    }
    modal_energy(runtime) + axial + t * entropy
}

#[test]
fn heating_changes_live_tension_pitch_and_pressure_without_resetting_motion() {
    with_string_cx(|cx, _| {
        let mut runtime = runtime(cx, 5.0);
        let original_frequency = runtime.mechanical().modes()[0].angular_frequency_rad_s;
        let original_claim = runtime.mechanical().specimen().material()
            .property("young_modulus").unwrap().answer().receipt.selected;
        let mut pressure_difference = 0.0;
        let mut entropy_exchange = 0.0;
        for epoch in 0..64 {
            let before = mechanical_internal_energy(&runtime);
            let mut frozen = runtime.mechanical().clone();
            let frozen_frame = frozen.step(cx, &[0.0]).unwrap();
            let frame = runtime.step(cx, epoch, &[0.0], 0.002).unwrap();
            // The thermal update changes the operator, never q or p.
            assert_eq!(runtime.mechanical().states(), frozen.states());
            assert_eq!(runtime.mechanical().epoch(), epoch + 1);
            assert_eq!(runtime.mechanical().accepted_samples(), epoch + 1);
            let specimen = runtime.mechanical().specimen();
            let t = runtime.thermal().temperature_k();
            let young = 2e9 + YOUNG_SLOPE * (t - 300.0);
            let strain = (0.502 - REFERENCE_LENGTH) / REFERENCE_LENGTH - ALPHA * (t - 300.0);
            close(specimen.string().tension_n, young * specimen.area_m2() * strain);
            close(runtime.mechanical().modes()[0].angular_frequency_rad_s,
                string_mode_omega(&specimen.string(), 1));
            close(frame.vibration.acoustic.total_modal_energy_j, modal_energy(&runtime));
            assert_eq!(specimen.material().property("young_modulus").unwrap()
                .answer().receipt.selected, original_claim);
            assert!(specimen.thermal_expansion().is_some());
            let independent = frame.vibration.acoustic.input_work_j + frame.external_heat_j
                - (mechanical_internal_energy(&runtime) - before + frame.enthalpy_change_j
                    + frame.vibration.dissipation.air_loss_j
                    + frame.vibration.dissipation.authored_loss_j
                    + frame.vibration.dissipation.roundoff_residual_j);
            let allowance = frame.energy_roundoff_tolerance_j + frame.thermal_solve_tolerance_j;
            assert!(frame.energy_balance_residual_j.abs() <= allowance);
            assert!(independent.abs() < allowance + 5e-13, "independent first law: {independent}");
            assert!((frame.reversible_heat_j + frame.modal_parameter_work_j
                + frame.axial_free_energy_change_j + frame.mechanical_entropy_energy_change_j).abs() < 5e-13);
            entropy_exchange += frame.mechanical_entropy_energy_change_j.abs();
            pressure_difference += (frame.vibration.acoustic.observer_pressure_pa
                - frozen_frame.acoustic.observer_pressure_pa).abs();
        }
        assert!(runtime.thermal().temperature_k() > 303.0);
        assert!(runtime.mechanical().modes()[0].angular_frequency_rad_s < 0.99 * original_frequency);
        assert!(entropy_exchange > 1e-5, "the entropy term must be physically exercised");
        assert!(pressure_difference > 1e-9);
    });
}

#[test]
fn unsourced_or_nonaffine_laws_refuse_and_slack_failure_preserves_live_history() {
    with_string_cx(|cx, _| {
        for (nonaffine, sourced) in [(true, true), (false, false)] {
            let (mechanical, card, curve) = parts(cx, nonaffine, sourced);
            assert!(ThermalMaterialStringRuntime::try_new_thermoelastic(cx, mechanical,
                card, curve, 0.0, 300.0, controls()).is_err());
        }
        let mut runtime = runtime(cx, 300.0);
        let mut retry = runtime.clone();
        // Caloric T stays inside the chart, but expansion at T>500 K would
        // remove all tensile prestress. No clamped frequency or new history.
        assert!(runtime.step(cx, 0, &[0.0], 8.0).is_err());
        assert_eq!(runtime.thermal(), retry.thermal());
        assert_eq!(runtime.mechanical().epoch(), retry.mechanical().epoch());
        assert_eq!(runtime.mechanical().specimen(), retry.mechanical().specimen());
        assert_eq!(runtime.mechanical().states(), retry.mechanical().states());
        assert_eq!(runtime.step(cx, 0, &[0.0], 0.002).unwrap(),
            retry.step(cx, 0, &[0.0], 0.002).unwrap());
    });
}

#[test]
fn electrical_transport_drives_the_same_thermoelastic_state_and_retained_sound() {
    with_string_cx(|cx, _| {
        let mut powered = runtime(cx, 5.0);
        let mut direct_heat = powered.clone();
        let original_frequency = powered.mechanical().modes()[0].angular_frequency_rad_s;
        let mut proposals = 0;
        for epoch in 0..64 {
            let frame = powered.step_with_current_and_thermal_transport(cx, epoch, &[0.0], 15.0,
                |input| {
                    proposals += 1;
                    Ok::<_, fs_material::phase::PhaseStateError>(UniformEnthalpyStep {
                        state: input.curve.advance_specific_energy(input.initial,
                            input.internal_heat_j / input.mass_kg)?,
                        external_heat_j: 0.0,
                        energy_residual_tolerance_j: 0.0,
                        report: (),
                    })
                }).unwrap();
            let reference = direct_heat.step(cx, epoch, &[0.0], frame.joule_heat_j).unwrap();
            close(frame.joule_heat_j, 15.0_f64.powi(2)
                * frame.conductor_before.resistance_ohm() / f64::from(RATE));
            assert_eq!(frame.coupled, reference);
            assert_eq!(powered.mechanical().states(), direct_heat.mechanical().states());
            assert!(frame.conductor_after.resistance_ohm() > frame.conductor_before.resistance_ohm());
            assert_eq!(frame.conductor_after.material().query_point(),
                powered.mechanical().specimen().material().query_point());
        }
        assert!(proposals > 64, "reversible heat must close with the real transport owner");
        assert!(powered.thermal().temperature_k() > 303.0);
        assert!(powered.mechanical().modes()[0].angular_frequency_rad_s < original_frequency);
    });
}
