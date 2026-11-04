//! Synthetic calibration fixtures; these are physical/ownership regressions,
//! not empirical qualification of a material specimen.
use core::convert::Infallible;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use fs_material::thermovisco::source::{
    ResolvedThermoMaxwell, SourceError, TemperatureShiftSourceQuery, ThermoMaxwellSourceQuery,
};
use fs_material::thermovisco::thermal::FeedbackConfig;
use fs_material::thermovisco::{HistoryIntegration, Limits, ThermalKnot, ThermoMaxwell};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::{Density, Dims, Pressure, Temperature, Time};

const ONE: Dims = Dims([0; 6]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const HEAT_CAPACITY: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn limits() -> Limits {
    Limits {
        max_strain: 0.04,
        max_strain_rate_s: 10.0,
    }
}

fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 100.0,
        },
        300.0,
        [280.0, 320.0],
    )
    .unwrap()
}

fn law(scales: [f64; 2], tau: f64, integration: HistoryIntegration) -> ThermoMaxwell {
    ThermoMaxwell::new(
        1e8,
        vec![(3e8, tau)],
        vec![
            ThermalKnot {
                temperature_k: 280.0,
                modulus_scale: 1.0,
                relaxation_scale: scales[0],
                free_strain: 0.0,
            },
            ThermalKnot {
                temperature_k: 320.0,
                modulus_scale: 1.0,
                relaxation_scale: scales[1],
                free_strain: 0.0,
            },
        ],
        limits(),
    )
    .unwrap()
    .with_history_integration(integration)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * actual.abs().max(expected.abs()).max(1.0),
        "{actual:e} != {expected:e}",
    );
}

#[test]
fn analytic_off_grid_kinetics_reach_both_original_integrators_and_tangents() {
    for integration in [
        HistoryIntegration::BackwardEuler,
        HistoryIntegration::ExactStrainRamp,
    ] {
        // The original non-unit reference multiplier must not be discarded.
        let reference = law([2.0; 2], 1.0, integration);
        let shifted = reference.clone().with_temperature_shift(shift()).unwrap();
        assert_eq!(shifted.temperature_shift(), Some(&shift()));
        assert_eq!(shifted.history_integration(), integration);
        let factor = shift().evaluate(307.25).unwrap().factor();
        let old = shifted.birth(300.0, 0.0).unwrap();
        let reference_old = reference.birth(300.0, 0.0).unwrap();
        let at_reference = old.advance(300.0, 0.01, 0.2, 1e-6).unwrap();
        let unchanged = reference_old.advance(300.0, 0.01, 0.2, 1e-6).unwrap();
        assert_eq!(
            at_reference.point.branch_strains(),
            unchanged.point.branch_strains()
        );
        assert_eq!(at_reference.stress_pa, unchanged.stress_pa);

        // At one temperature the exact horizontal shift equals reference-time
        // rescaling, including history, physical heat and work-conjugate stress.
        let actual = old.advance(307.25, 0.01, 0.2, 1e-6).unwrap();
        let expected = reference_old
            .advance(307.25, 0.01, 0.2 / factor, 1e-6)
            .unwrap();
        close(
            actual.point.branch_strains()[0],
            expected.point.branch_strains()[0],
        );
        close(actual.stress_pa, expected.stress_pa);
        close(
            actual.work_conjugate_stress_pa,
            expected.work_conjugate_stress_pa,
        );
        close(actual.energy.dissipated, expected.energy.dissipated);
        close(
            actual.energy.numerical_dissipation,
            expected.energy.numerical_dissipation,
        );
        for (a, b) in [
            (
                old.linearize(307.25, 0.2).unwrap(),
                reference_old.linearize(307.25, 0.2 / factor).unwrap(),
            ),
            (
                old.linearize_work(307.25, 0.2).unwrap(),
                reference_old.linearize_work(307.25, 0.2 / factor).unwrap(),
            ),
        ] {
            close(a.tangent_pa, b.tangent_pa);
            close(a.offset_pa, b.offset_pa);
        }

        // An endpoint-resampled WLF curve would give a materially wrong interior
        // response. Two thermal knots must still evaluate the analytic kinetics.
        let sampled = law(
            [
                2.0 * shift().evaluate(280.0).unwrap().factor(),
                2.0 * shift().evaluate(320.0).unwrap().factor(),
            ],
            1.0,
            integration,
        );
        let sampled_step = sampled
            .birth(300.0, 0.0)
            .unwrap()
            .advance(307.25, 0.01, 0.2, 1e-6)
            .unwrap();
        assert!(
            (actual.point.branch_strains()[0] - sampled_step.point.branch_strains()[0]).abs()
                > 1e-4
        );
        assert_eq!(old.branch_strains(), &[0.0]);
    }
}

#[test]
fn heating_and_cooling_preserve_branch_memory_and_only_release_stored_energy() {
    let material = law([2.0; 2], 1.0, HistoryIntegration::ExactStrainRamp)
        .with_temperature_shift(shift())
        .unwrap();
    let mut accepted = material
        .birth(300.0, 0.0)
        .unwrap()
        .advance(300.0, 0.02, 0.05, 1e-6)
        .unwrap()
        .point;
    for (temperature, duration) in [(312.5, 0.2), (289.0, 0.35), (305.0, 0.1)] {
        let old = accepted.clone();
        let q0 = old.branch_strains()[0];
        let tau = 2.0 * shift().evaluate(temperature).unwrap().factor();
        let q1 = q0 * (-duration / tau).exp();
        let step = old.advance(temperature, 0.02, duration, 1e-6).unwrap();
        close(step.point.branch_strains()[0], q1);
        close(step.energy.dissipated, 0.5 * 3e8 * (q0 * q0 - q1 * q1));
        close(step.energy.stored_change, -step.energy.dissipated);
        assert_eq!(step.energy.mechanical_work, 0.0);
        assert_eq!(step.energy.thermal_strain_work, 0.0);
        assert_eq!(step.energy.modulus_work, 0.0);
        assert_eq!(step.energy.numerical_dissipation, 0.0);
        assert!(step.energy.dissipated > 0.0 && step.energy.residual.abs() < 1e-6);
        assert_eq!(accepted, old);
        accepted = step.point;
    }
    let snapshot = accepted.clone();
    assert!(accepted.advance(325.0, 0.02, 0.1, 1e-6).is_err());
    assert_eq!(accepted, snapshot);
}

#[test]
fn attachment_refuses_double_scaling_incompatible_support_and_unrepresentable_times() {
    let integration = HistoryIntegration::ExactStrainRamp;
    assert!(
        law([1.0, 0.5], 1.0, integration)
            .with_temperature_shift(shift())
            .is_err()
    );
    let original = law([1.0; 2], 1.0, integration);
    let shifted = original.clone().with_temperature_shift(shift()).unwrap();
    assert!(shifted.with_temperature_shift(shift()).is_err());
    assert!(original.temperature_shift().is_none());
    for (reference, band) in [(300.0, [290.0, 310.0]), (330.0, [270.0, 340.0])] {
        let incompatible = TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 40000.0,
            },
            reference,
            band,
        )
        .unwrap();
        assert!(
            original
                .clone()
                .with_temperature_shift(incompatible)
                .is_err()
        );
    }
    let strong = TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 10.0,
            c2_k: 100.0,
        },
        300.0,
        [280.0, 320.0],
    )
    .unwrap();
    assert!(
        law([1.0; 2], 1e307, integration)
            .with_temperature_shift(strong)
            .is_err()
    );
    assert!(
        law([1.0; 2], f64::from_bits(1), integration)
            .with_temperature_shift(strong)
            .is_err()
    );
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}

fn add(claims: &mut ClaimSet, name: &str, dims: Dims, knots: Vec<(f64, f64)>) {
    let low = knots[0].0;
    let high = knots.last().unwrap().0;
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", low, high),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic temperature-shift calibration".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
}

fn card(c1_knots: Vec<(f64, f64)>, time_knots: Vec<(f64, f64)>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("E_inf", Pressure::DIMS, 1e8),
        ("E0", Pressure::DIMS, 3e8),
        ("tau0", Time::DIMS, 1.0),
        ("scale", ONE, 1.0),
        ("free", ONE, 0.0),
        ("C2", Temperature::DIMS, 100.0),
        ("activation", MOLAR_ENERGY, 40000.0),
        ("density", Density::DIMS, 1000.0),
        ("cp", HEAT_CAPACITY, 5.0),
    ] {
        add(
            &mut claims,
            name,
            dims,
            vec![(280.0, value), (320.0, value)],
        );
    }
    add(
        &mut claims,
        "T_ref",
        Temperature::DIMS,
        vec![(280.0, 300.0), (293.0, 300.0), (320.0, 300.0)],
    );
    add(&mut claims, "C1", ONE, c1_knots);
    add(&mut claims, "time_scale", ONE, time_knots);
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "original calibrated Maxwell test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn constant_card(c1: f64) -> MaterialCard {
    card(
        vec![(280.0, c1), (307.0, c1), (320.0, c1)],
        vec![(280.0, 1.0), (320.0, 1.0)],
    )
}

fn resolve(card: &MaterialCard) -> ResolvedThermoMaxwell {
    ResolvedThermoMaxwell::try_new(
        card,
        ThermoMaxwellSourceQuery {
            equilibrium_modulus: &PropertyKey::new("E_inf", Pressure::DIMS),
            branches: &[(
                PropertyKey::new("E0", Pressure::DIMS),
                PropertyKey::new("tau0", Time::DIMS),
            )],
            modulus_scale: &PropertyKey::new("scale", ONE),
            relaxation_scale: &PropertyKey::new("time_scale", ONE),
            free_strain: &PropertyKey::new("free", ONE),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        limits(),
        || false,
    )
    .unwrap()
}

fn bind(
    source: ResolvedThermoMaxwell,
    card: &MaterialCard,
) -> Result<ResolvedThermoMaxwell, SourceError> {
    source.with_temperature_shift_from_source(
        card,
        TemperatureShiftSourceQuery::WilliamsLandelFerry {
            reference_temperature: &PropertyKey::new("T_ref", Temperature::DIMS),
            c1: &PropertyKey::new("C1", ONE),
            c2: &PropertyKey::new("C2", Temperature::DIMS),
        },
        SelectionPolicy::SingleClaimOnly,
        || false,
    )
}

#[test]
fn calibration_keeps_original_card_complete_knots_keys_and_mechanical_receipts() {
    let card = constant_card(2.0);
    let original = resolve(&card);
    let original_samples = original.samples().to_vec();
    let original_keys = original.keys().to_vec();
    let source = bind(original, &card)
        .unwrap()
        .with_history_integration(HistoryIntegration::ExactStrainRamp);
    assert_eq!(source.samples(), original_samples.as_slice());
    assert_eq!(source.keys(), original_keys.as_slice());
    let calibration = source.temperature_shift_source().unwrap();
    assert_eq!(calibration.declaration(), &shift());
    assert_eq!(
        source.law().temperature_shift(),
        Some(calibration.declaration())
    );
    let temperatures: Vec<_> = calibration
        .samples()
        .iter()
        .map(|sample| sample.state_coordinate("T").unwrap())
        .collect();
    assert_eq!(temperatures, vec![280.0, 293.0, 300.0, 307.0, 320.0]);
    assert_eq!(
        calibration
            .keys()
            .iter()
            .map(|key| key.name())
            .collect::<Vec<_>>(),
        vec!["T_ref", "C1", "C2"]
    );
    for sample in calibration.samples() {
        assert_eq!(sample.card_identity(), card.content_hash());
        source.admit_state_point(sample).unwrap();
        for key in calibration.keys() {
            let actual = sample.property_by_key(key).unwrap();
            let first = calibration.samples()[0].property_by_key(key).unwrap();
            assert_eq!(
                actual.answer().receipt.selected,
                first.answer().receipt.selected
            );
        }
    }
    let direct = law([1.0; 2], 1.0, HistoryIntegration::ExactStrainRamp)
        .with_temperature_shift(shift())
        .unwrap();
    let step = source
        .law()
        .birth(300.0, 0.0)
        .unwrap()
        .advance(307.25, 0.02, 0.2, 1e-6)
        .unwrap();
    let expected = direct
        .birth(300.0, 0.0)
        .unwrap()
        .advance(307.25, 0.02, 0.2, 1e-6)
        .unwrap();
    assert_eq!(step.point.branch_strains(), expected.point.branch_strains());
    assert_eq!(step.stress_pa, expected.stress_pa);
}

#[test]
fn original_arrhenius_molar_energy_reaches_the_same_history_owner() {
    let card = constant_card(2.0);
    let source = resolve(&card)
        .with_temperature_shift_from_source(
            &card,
            TemperatureShiftSourceQuery::Arrhenius {
                reference_temperature: &PropertyKey::new("T_ref", Temperature::DIMS),
                activation_energy: &PropertyKey::new("activation", MOLAR_ENERGY),
            },
            SelectionPolicy::SingleClaimOnly,
            || false,
        )
        .unwrap();
    let expected = TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 40000.0,
        },
        300.0,
        [280.0, 320.0],
    )
    .unwrap();
    assert_eq!(
        source.temperature_shift_source().unwrap().declaration(),
        &expected
    );
    let at = 311.5;
    let tau = expected.evaluate(at).unwrap().factor();
    let step = source
        .law()
        .birth(300.0, 0.0)
        .unwrap()
        .advance(at, 0.02, 0.2, 1e-6)
        .unwrap();
    close(step.point.branch_strains()[0], 0.02 * tau / (tau + 0.2));
    // Energy per particle cannot silently be interpreted as molar activation.
    assert!(
        resolve(&card)
            .with_temperature_shift_from_source(
                &card,
                TemperatureShiftSourceQuery::Arrhenius {
                    reference_temperature: &PropertyKey::new("T_ref", Temperature::DIMS),
                    activation_energy: &PropertyKey::new("activation", Dims([2, 1, -2, 0, 0, 0])),
                },
                SelectionPolicy::SingleClaimOnly,
                || false,
            )
            .is_err()
    );
}

#[test]
fn calibration_refuses_hidden_variation_support_holes_foreign_cards_and_double_shifts() {
    let constant_time = vec![(280.0, 1.0), (320.0, 1.0)];
    for c1 in [
        vec![(280.0, 2.0), (303.0, 2.1), (320.0, 2.0)],
        vec![(290.0, 2.0), (310.0, 2.0)],
    ] {
        let material = card(c1, constant_time.clone());
        assert!(bind(resolve(&material), &material).is_err());
    }
    let material = card(
        vec![(280.0, 2.0), (320.0, 2.0)],
        vec![(280.0, 1.0), (300.0, 0.8), (320.0, 1.0)],
    );
    assert!(bind(resolve(&material), &material).is_err());
    let original_card = constant_card(2.0);
    let foreign = constant_card(3.0);
    let original = resolve(&original_card);
    assert!(bind(original.clone(), &foreign).is_err());
    assert!(
        bind(
            bind(original.clone(), &original_card).unwrap(),
            &original_card
        )
        .is_err()
    );
    assert!(original.law().temperature_shift().is_none());
    assert!(original.temperature_shift_source().is_none());
    assert!(matches!(
        original.clone().with_temperature_shift_from_source(
            &original_card,
            TemperatureShiftSourceQuery::WilliamsLandelFerry {
                reference_temperature: &PropertyKey::new("T_ref", Temperature::DIMS),
                c1: &PropertyKey::new("C1", ONE),
                c2: &PropertyKey::new("C2", Temperature::DIMS),
            },
            SelectionPolicy::SingleClaimOnly,
            || true,
        ),
        Err(SourceError::Cancelled),
    ));
    assert!(original.law().temperature_shift().is_none());
}

fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("cp", HEAT_CAPACITY),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(280.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -100.0,
        1e-7,
    )
    .unwrap()
}

fn adiabatic(input: UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<()>, Infallible> {
    Ok(UniformEnthalpyStep {
        state: input
            .curve
            .state_at_specific_enthalpy(
                input.initial.specific_enthalpy_j_kg() + input.internal_heat_j / input.mass_kg,
            )
            .unwrap(),
        external_heat_j: 0.0,
        energy_residual_tolerance_j: 1e-10,
        report: (),
    })
}

fn feedback() -> FeedbackConfig {
    FeedbackConfig {
        max_iterations: 100,
        temperature_tolerance_k: 1e-10,
        energy_tolerance_j: 1e-9,
        relaxation: 1.0,
    }
}

#[test]
fn original_calibration_drives_actual_self_heating_with_no_stale_temperature_or_extra_heat() {
    for integration in [
        HistoryIntegration::BackwardEuler,
        HistoryIntegration::ExactStrainRamp,
    ] {
        let mut outputs = Vec::new();
        for c1 in [0.5, 3.0] {
            let card = constant_card(c1);
            let source = bind(resolve(&card), &card)
                .unwrap()
                .with_history_integration(integration);
            let thermal = caloric(&card);
            let old = source
                .caloric_point(&thermal, 300.0, 0.0, 1e-6, 1e-3)
                .unwrap();
            let step = old
                .advance(0.02, 0.2, feedback(), adiabatic, || false)
                .unwrap();
            let temperature = step.next.caloric_state().temperature_k();
            assert!(temperature > 300.0);
            assert!(step.iterations > 1);
            let actual = old
                .mechanical()
                .advance(temperature, 0.02, 0.2, 1e-3)
                .unwrap();
            assert_eq!(step.next.mechanical(), &actual.point);
            close(step.applied_heat_j, 1e-6 * actual.energy.dissipated);
            assert_eq!(step.reversible_heat_j, 0.0);
            assert!(step.combined_residual_j.abs() < 1e-9);
            assert!((temperature - (300.0 + step.applied_heat_j / 0.005)).abs() < 1e-8);
            if integration == HistoryIntegration::ExactStrainRamp {
                assert_eq!(step.mechanical_energy.numerical_dissipation, 0.0);
            } else {
                assert!(step.mechanical_energy.numerical_dissipation > 0.0);
            }
            // Same accepted history underlies every trial and a cancelled trial
            // cannot change the accepted strain, memory, temperature or clock.
            assert_eq!(old.time_s(), 0.0);
            assert_eq!(old.mechanical().branch_strains(), &[0.0]);
            assert!(
                old.advance(0.02, 0.2, feedback(), adiabatic, || true)
                    .is_err()
            );
            assert_eq!(old.time_s(), 0.0);
            outputs.push((temperature, step.next.mechanical().stress_pa()));
        }
        // In this deliberately bounded small-dt experiment faster hot kinetics
        // increases physical heat and lowers the retained Maxwell stress.
        assert!(outputs[1].0 > outputs[0].0);
        assert!(outputs[1].1 < outputs[0].1);
    }
}
