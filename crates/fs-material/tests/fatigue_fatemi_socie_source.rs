//! Original engineering-shear strain-life fits and native J2 source binding.
//! Synthetic coefficients exercise the implementation; they are not measured data.
#![deny(unsafe_code)]

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialAnswer, MaterialCard,
    MaterialStateId, ObservationDataset, ObservationId, PropertyClaim, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::critical_plane::fatemi_socie::FatemiSocieNormalStress;
use fs_material::fatigue::critical_plane::fatemi_socie::source::{
    FatemiSociePowerBasis, FatemiSocieProperty, FatemiSocieSourceRequest,
    resolve_isotropic_fatemi_socie,
};
use fs_material::fatigue::strain_life::source::CmbProperty;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const PHYSICAL: [FatemiSocieProperty; 3] = [
    FatemiSocieProperty::YoungsModulus,
    FatemiSocieProperty::PoissonRatio,
    FatemiSocieProperty::InitialYieldStress,
];

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn request(basis: FatemiSociePowerBasis) -> FatemiSocieSourceRequest {
    FatemiSocieSourceRequest {
        power_basis: basis,
        normal_stress: FatemiSocieNormalStress::Signed,
        temperature_axis: "specimen-temperature".into(),
        shear_stress_ratio_axis: "torsion-calibration-ratio".into(),
        surface_and_geometry: "synthetic smooth torsion coupon".into(),
        environment: "synthetic dry atmosphere".into(),
        load_basis: "engineering shear strain amplitude; fully reversed torsion".into(),
    }
}

fn state(request: &FatemiSocieSourceRequest, temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            &request.temperature_axis,
            absolute_temperature(),
            temperature,
        )
        .unwrap()
}

fn point(request: &FatemiSocieSourceRequest) -> QueryPoint {
    state(request, 300.0)
        .with_quantity(
            &request.shear_stress_ratio_axis,
            QuantitySpec::dimensional(Dims::NONE),
            -1.0,
        )
        .unwrap()
        .with_quantity(
            request.power_basis.life_axis(),
            QuantitySpec::dimensional(Dims::NONE),
            100.0,
        )
        .unwrap()
}

fn temperature_domain(request: &FatemiSocieSourceRequest) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        &request.temperature_axis,
        absolute_temperature(),
        300.0,
        400.0,
    )
}

fn rows(request: &FatemiSocieSourceRequest) -> Vec<PropertyClaim> {
    let factor = if request.power_basis == FatemiSociePowerBasis::FullCycles {
        2.0_f64.powf(-0.5)
    } else {
        1.0
    };
    let life = if request.power_basis == FatemiSociePowerBasis::FullCycles {
        (0.5, 10_000.0)
    } else {
        (1.0, 20_000.0)
    };
    let values = [
        1200.0,
        0.2,
        100.0,
        10.0 * factor,
        -0.5,
        0.03 * factor,
        -0.5,
        0.2,
    ];
    FatemiSocieProperty::ALL
        .into_iter()
        .map(|property| {
            let key = property.property_key(request.normal_stress);
            let mut validity = temperature_domain(request);
            if property.index() >= 3 {
                validity = validity
                    .with_quantity(
                        &request.shear_stress_ratio_axis,
                        QuantitySpec::dimensional(Dims::NONE),
                        -1.0,
                        -1.0,
                    )
                    .with_quantity(
                        request.power_basis.life_axis(),
                        QuantitySpec::dimensional(Dims::NONE),
                        life.0,
                        life.1,
                    );
            }
            PropertyClaim {
                value: PropertyValue::Scalar {
                    value: values[property.index()],
                    dims: key.dims(),
                },
                key,
                validity,
                uncertainty: UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.05,
                    confidence: 0.9,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: match property {
                        FatemiSocieProperty::YoungsModulus => "synthetic elastic modulus",
                        FatemiSocieProperty::PoissonRatio => "synthetic Poisson ratio",
                        FatemiSocieProperty::InitialYieldStress => "synthetic initial yield stress",
                        FatemiSocieProperty::NormalStressSensitivity => {
                            "synthetic separate normal-sensitivity fit"
                        }
                        _ => "synthetic common torsional strain-life fit",
                    }
                    .into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            }
        })
        .collect()
}

fn card_with_observations(
    rows: Vec<PropertyClaim>,
    observations: Vec<ObservationDataset>,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for observation in observations {
        claims.register_observation(observation).unwrap();
    }
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "synthetic coupon process".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    card_with_observations(rows, vec![])
}

fn observation(number: u8) -> ObservationDataset {
    ObservationDataset {
        specimen: format!("synthetic population {number}"),
        method: "synthetic fixture; not a measured test".into(),
        artifact: ContentHash([number; 32]),
        caveats: "test-only".into(),
        provenance: Provenance {
            source: "synthetic observations".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn native_answers(
    material: &MaterialCard,
    request: &FatemiSocieSourceRequest,
    at: &QueryPoint,
) -> [MaterialAnswer; 3] {
    PHYSICAL.map(|property| {
        material
            .claims()
            .query_typed(
                &property.property_key(request.normal_stress),
                at,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap()
    })
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1e-10),
        "{actual:e} vs {expected:e}"
    );
}

#[test]
fn full_cycle_and_reversal_fits_preserve_engineering_shear_and_bounded_life() {
    let cycles_request = request(FatemiSociePowerBasis::FullCycles);
    let reversals_request = request(FatemiSociePowerBasis::Reversals);
    let cycles_card = card(rows(&cycles_request));
    let reversals_card = card(rows(&reversals_request));
    let cycles = resolve_isotropic_fatemi_socie(
        &cycles_card,
        &point(&cycles_request),
        &cycles_request,
        [SINGLE; 8],
    )
    .unwrap();
    let reversals = resolve_isotropic_fatemi_socie(
        &reversals_card,
        &point(&reversals_request),
        &reversals_request,
        [SINGLE; 8],
    )
    .unwrap();
    assert_eq!(cycles.cycle_domain(), (0.5, 10_000.0));
    assert_eq!(cycles.cycle_domain(), reversals.cycle_domain());
    close(cycles.parameters().shear_modulus_pa, 500.0);
    close(
        cycles.parameters().shear_fatigue_strength_coefficient_pa,
        10.0,
    );
    close(
        cycles.parameters().shear_fatigue_ductility_coefficient,
        0.03,
    );
    for n in [0.5, 1.0, 50.0, 1000.0, 10_000.0] {
        // G=500, tau_f'=10, gamma_f'=0.03 and b=c=-0.5 give this
        // independent engineering-shear amplitude, with no tensor-shear factor.
        let expected = 0.05 / (2.0_f64 * n).sqrt();
        let cycles_parameter = cycles.parameter_at_cycles(n).unwrap();
        let reversals_parameter = reversals.parameter_at_cycles(n).unwrap();
        close(cycles_parameter, expected);
        close(reversals_parameter, expected);
        close(
            cycles.cycles_to_failure(cycles_parameter).unwrap().unwrap(),
            n,
        );
        close(
            reversals
                .cycles_to_failure(reversals_parameter)
                .unwrap()
                .unwrap(),
            n,
        );
    }
    close(
        cycles
            .source()
            .unwrap()
            .property(FatemiSocieProperty::ShearFatigueStrengthCoefficient)
            .evidence
            .value
            .value,
        10.0 / 2.0_f64.sqrt(),
    );
    close(
        cycles
            .source()
            .unwrap()
            .property(FatemiSocieProperty::ShearFatigueDuctilityCoefficient)
            .evidence
            .value
            .value,
        0.03 / 2.0_f64.sqrt(),
    );
    assert!(cycles.parameter_at_cycles(0.49).is_err());
    assert!(cycles.parameter_at_cycles(10_001.0).is_err());
    assert!(cycles.cycles_to_failure(0.06).is_err());
}

#[test]
fn all_original_claims_uncertainty_and_physical_source_receipts_survive() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let original = rows(&request);
    let material = card(original.clone());
    let at = point(&request);
    let calibration =
        resolve_isotropic_fatemi_socie(&material, &at, &request, [SINGLE; 8]).unwrap();
    let source = calibration.source().unwrap();
    assert_eq!(source.card_identity(), material.content_hash());
    assert_eq!(source.material_state(), material.id());
    assert_eq!(source.request(), &request);
    assert_eq!(source.query(), &at);
    assert_eq!(source.selections(), &[SINGLE; 8]);
    for property in FatemiSocieProperty::ALL {
        assert_eq!(source.claim(property), &original[property.index()]);
        let answer = material
            .claims()
            .query_typed(
                &property.property_key(request.normal_stress),
                &at,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(source.property(property), &answer);
        material
            .claims()
            .verify_receipt(&source.property(property).receipt)
            .unwrap();
    }
    close(
        calibration.parameters().shear_modulus_pa,
        1200.0 / (2.0 * (1.0 + 0.2)),
    );
    close(calibration.parameters().initial_yield_stress_pa, 100.0);
    assert_eq!(
        calibration.parameters().normal_stress,
        FatemiSocieNormalStress::Signed
    );
    assert_eq!(calibration.scope().material, material.id().chemistry);
    assert_eq!(calibration.scope().process, material.id().process);
    assert_eq!(calibration.scope().temperature_k, 300.0);
}

#[test]
fn source_temperature_interpolation_drives_actual_isotropic_and_fatigue_parameters() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let mut original = rows(&request);
    for (property, a, b) in [
        (FatemiSocieProperty::YoungsModulus, 1200.0, 1000.0),
        (FatemiSocieProperty::PoissonRatio, 0.2, 0.3),
        (FatemiSocieProperty::InitialYieldStress, 100.0, 80.0),
        (
            FatemiSocieProperty::ShearFatigueStrengthCoefficient,
            10.0,
            8.0,
        ),
    ] {
        let row = &mut original[property.index()];
        row.value = PropertyValue::Curve {
            abscissa: request.temperature_axis.clone(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, a), (400.0, b)],
            dims: row.key.dims(),
        };
        row.interpolation = InterpolationPolicy::LinearInside;
    }
    let material = card(original.clone());
    let at = point(&request)
        .with(&request.temperature_axis, 350.0)
        .unwrap();
    let calibration =
        resolve_isotropic_fatemi_socie(&material, &at, &request, [SINGLE; 8]).unwrap();
    close(calibration.parameters().shear_modulus_pa, 440.0);
    close(calibration.parameters().initial_yield_stress_pa, 90.0);
    close(
        calibration
            .parameters()
            .shear_fatigue_strength_coefficient_pa,
        9.0,
    );
    close(
        calibration.parameter_at_cycles(50.0).unwrap(),
        (9.0 / 440.0 + 0.03) / 10.0,
    );
    assert_eq!(calibration.cycle_domain(), (0.5, 10_000.0));
    for property in FatemiSocieProperty::ALL {
        assert_eq!(
            calibration.source().unwrap().claim(property),
            &original[property.index()]
        );
        material
            .claims()
            .verify_receipt(&calibration.source().unwrap().property(property).receipt)
            .unwrap();
    }
    assert!(
        resolve_isotropic_fatemi_socie(
            &material,
            &at.with(&request.temperature_axis, 401.0).unwrap(),
            &request,
            [SINGLE; 8]
        )
        .is_err()
    );
}

#[test]
fn missing_shear_data_are_not_replaced_by_axial_fits_and_ambiguity_needs_pins() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let at = point(&request);
    for property in FatemiSocieProperty::ALL {
        let mut incomplete = rows(&request);
        incomplete.remove(property.index());
        assert!(
            resolve_isotropic_fatemi_socie(&card(incomplete), &at, &request, [SINGLE; 8]).is_err()
        );
        let mut original = rows(&request);
        let keep = ClaimId(original[property.index()].content_hash());
        let mut alternative = original[property.index()].clone();
        alternative.uncertainty = UncertaintyModel::Unstated;
        original.push(alternative);
        let material = card(original);
        assert!(resolve_isotropic_fatemi_socie(&material, &at, &request, [SINGLE; 8]).is_err());
        let mut selections = [SINGLE; 8];
        selections[property.index()] = ClaimSelection::Pinned(keep);
        let calibration =
            resolve_isotropic_fatemi_socie(&material, &at, &request, selections).unwrap();
        assert_eq!(
            calibration
                .source()
                .unwrap()
                .property(property)
                .receipt
                .selected,
            keep
        );
    }
    for (shear, axial) in [
        (
            FatemiSocieProperty::ShearFatigueStrengthCoefficient,
            CmbProperty::FatigueStrengthCoefficient,
        ),
        (
            FatemiSocieProperty::ShearFatigueDuctilityCoefficient,
            CmbProperty::FatigueDuctilityCoefficient,
        ),
    ] {
        let mut original = rows(&request);
        original[shear.index()].key = axial.property_key();
        assert!(
            resolve_isotropic_fatemi_socie(&card(original), &at, &request, [SINGLE; 8]).is_err()
        );
    }
}

#[test]
fn independent_sensitivity_fit_limits_life_and_cannot_switch_normal_convention() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let mut original = rows(&request);
    let observations = vec![observation(1), observation(2)];
    for row in &mut original[3..7] {
        row.observations = vec![ObservationId(observations[0].content_hash())];
    }
    let sensitivity = &mut original[FatemiSocieProperty::NormalStressSensitivity.index()];
    sensitivity.observations = vec![ObservationId(observations[1].content_hash())];
    sensitivity.validity = sensitivity.validity.clone().with_quantity(
        request.power_basis.life_axis(),
        QuantitySpec::dimensional(Dims::NONE),
        20.0,
        2000.0,
    );
    let material = card_with_observations(original.clone(), observations);
    let calibration =
        resolve_isotropic_fatemi_socie(&material, &point(&request), &request, [SINGLE; 8]).unwrap();
    assert_eq!(calibration.cycle_domain(), (10.0, 1000.0));
    assert_eq!(
        calibration
            .source()
            .unwrap()
            .claim(FatemiSocieProperty::NormalStressSensitivity),
        &original[7]
    );
    assert!(calibration.parameter_at_cycles(9.0).is_err());
    assert!(calibration.parameter_at_cycles(1001.0).is_err());
    let mut tension = request.clone();
    tension.normal_stress = FatemiSocieNormalStress::TensionOnly;
    assert_ne!(
        FatemiSocieProperty::NormalStressSensitivity.property_key(request.normal_stress),
        FatemiSocieProperty::NormalStressSensitivity.property_key(tension.normal_stress)
    );
    assert!(
        resolve_isotropic_fatemi_socie(&material, &point(&tension), &tension, [SINGLE; 8]).is_err()
    );
    let tension_card = card(rows(&tension));
    let explicit =
        resolve_isotropic_fatemi_socie(&tension_card, &point(&tension), &tension, [SINGLE; 8])
            .unwrap();
    assert_eq!(
        explicit.parameters().normal_stress,
        FatemiSocieNormalStress::TensionOnly
    );
    assert_eq!(
        explicit.source().unwrap().request().normal_stress,
        FatemiSocieNormalStress::TensionOnly
    );
}

#[test]
fn torsional_coefficients_require_one_fit_population_and_life_interval() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let at = point(&request);
    let mut mixed = rows(&request);
    mixed[6].provenance.source = "different synthetic torsion fit".into();
    assert!(resolve_isotropic_fatemi_socie(&card(mixed), &at, &request, [SINGLE; 8]).is_err());
    let mut mixed = rows(&request);
    let observations = vec![observation(1), observation(2)];
    for row in &mut mixed[3..7] {
        row.observations = vec![ObservationId(observations[0].content_hash())];
    }
    mixed[6].observations = vec![ObservationId(observations[1].content_hash())];
    assert!(
        resolve_isotropic_fatemi_socie(
            &card_with_observations(mixed, observations),
            &at,
            &request,
            [SINGLE; 8]
        )
        .is_err()
    );
    let mut mixed = rows(&request);
    mixed[6].validity = mixed[6].validity.clone().with_quantity(
        request.power_basis.life_axis(),
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
        10_000.0,
    );
    assert!(resolve_isotropic_fatemi_socie(&card(mixed), &at, &request, [SINGLE; 8]).is_err());
}

#[test]
fn intersection_cannot_hide_invalid_original_life_support() {
    for basis in [
        FatemiSociePowerBasis::FullCycles,
        FatemiSociePowerBasis::Reversals,
    ] {
        let request = request(basis);
        for invalid_group in [3..7, 7..8] {
            let mut original = rows(&request);
            for row in &mut original[invalid_group] {
                row.validity = row.validity.clone().with_quantity(
                    request.power_basis.life_axis(),
                    QuantitySpec::dimensional(Dims::NONE),
                    0.1,
                    20_000.0,
                );
            }
            // These finite source intervals are accepted by matdb and contain
            // the query, but invalid original support must not be repaired by
            // intersecting with the other fit's valid half-cycle lower bound.
            let material = card(original);
            assert!(
                resolve_isotropic_fatemi_socie(&material, &point(&request), &request, [SINGLE; 8])
                    .is_err()
            );
        }
    }
}

#[test]
fn typed_temperature_reversed_baseline_and_parameter_dependencies_are_enforced() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let material = card(rows(&request));
    let at = point(&request);
    for bad in [
        at.clone()
            .with(&request.shear_stress_ratio_axis, 0.0)
            .unwrap(),
        at.clone()
            .with_quantity(
                &request.temperature_axis,
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
            )
            .unwrap(),
        at.clone()
            .with_quantity(
                &request.temperature_axis,
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::TemperatureDifference,
                    ValueForm::Static,
                )),
                300.0,
            )
            .unwrap(),
        state(&request, 300.0)
            .with(&request.shear_stress_ratio_axis, -1.0)
            .unwrap()
            .with(request.power_basis.life_axis(), 100.0)
            .unwrap(),
    ] {
        assert!(resolve_isotropic_fatemi_socie(&material, &bad, &request, [SINGLE; 8]).is_err());
    }
    for property in FatemiSocieProperty::ALL {
        let mut original = rows(&request);
        original[property.index()].validity = ValidityDomain::unconstrained();
        assert!(
            resolve_isotropic_fatemi_socie(&card(original), &at, &request, [SINGLE; 8]).is_err()
        );
    }
    let mut unsupported = rows(&request);
    unsupported[3].value = PropertyValue::Curve {
        abscissa: request.power_basis.life_axis().into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(1.0, 10.0), (20_000.0, 8.0)],
        dims: unsupported[3].key.dims(),
    };
    unsupported[3].interpolation = InterpolationPolicy::LinearInside;
    assert!(
        resolve_isotropic_fatemi_socie(&card(unsupported), &at, &request, [SINGLE; 8]).is_err()
    );
    for property in PHYSICAL {
        let mut original = rows(&request);
        original[property.index()].validity =
            original[property.index()].validity.clone().with_quantity(
                request.power_basis.life_axis(),
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                20_000.0,
            );
        assert!(
            resolve_isotropic_fatemi_socie(&card(original), &at, &request, [SINGLE; 8]).is_err()
        );
    }
}

#[test]
fn nonphysical_isotropic_yield_or_fatigue_parameters_refuse() {
    let request = request(FatemiSociePowerBasis::Reversals);
    for (property, invalid) in [
        (FatemiSocieProperty::YoungsModulus, 0.0),
        (FatemiSocieProperty::PoissonRatio, -1.0),
        (FatemiSocieProperty::PoissonRatio, 0.5),
        (FatemiSocieProperty::InitialYieldStress, 0.0),
        (FatemiSocieProperty::ShearFatigueStrengthCoefficient, -1.0),
        (FatemiSocieProperty::ShearFatigueDuctilityCoefficient, -0.1),
        (FatemiSocieProperty::NormalStressSensitivity, -0.1),
    ] {
        let mut original = rows(&request);
        let row = &mut original[property.index()];
        row.value = PropertyValue::Scalar {
            value: invalid,
            dims: row.key.dims(),
        };
        assert!(
            resolve_isotropic_fatemi_socie(
                &card(original),
                &point(&request),
                &request,
                [SINGLE; 8]
            )
            .is_err()
        );
    }
}

#[test]
fn native_binding_requires_the_actual_elastic_poisson_and_initial_yield_sources() {
    let request = request(FatemiSociePowerBasis::Reversals);
    let original = rows(&request);
    let material = card(original.clone());
    let calibration =
        resolve_isotropic_fatemi_socie(&material, &point(&request), &request, [SINGLE; 8]).unwrap();
    let source = calibration.source().unwrap();
    let native = state(&request, 300.0);
    let answers = native_answers(&material, &request, &native);
    assert!(source.matches_material_state(
        material.content_hash(),
        material.id(),
        &native,
        &answers[0],
        &answers[1],
        &answers[2]
    ));
    let extra = native
        .clone()
        .with_quantity(
            "native-only-coordinate",
            QuantitySpec::dimensional(Dims::NONE),
            0.25,
        )
        .unwrap();
    assert!(source.matches_material_state(
        material.content_hash(),
        material.id(),
        &extra,
        &answers[0],
        &answers[1],
        &answers[2]
    ));
    assert!(!source.matches_material_state(
        ContentHash([77; 32]),
        material.id(),
        &native,
        &answers[0],
        &answers[1],
        &answers[2]
    ));
    let mut other_state = material.id().clone();
    other_state.revision += 1;
    assert!(!source.matches_material_state(
        material.content_hash(),
        &other_state,
        &native,
        &answers[0],
        &answers[1],
        &answers[2]
    ));
    for different in [
        state(&request, 301.0),
        QueryPoint::new()
            .with_quantity(
                &request.temperature_axis,
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
            )
            .unwrap(),
    ] {
        assert!(!source.matches_material_state(
            material.content_hash(),
            material.id(),
            &different,
            &answers[0],
            &answers[1],
            &answers[2]
        ));
    }
    for property in PHYSICAL {
        let mut changed = original.clone();
        changed[property.index()].provenance.source =
            "independent physical source with coincident value".into();
        let other = card(changed);
        let mut other_answers = answers.clone();
        other_answers[property.index()] = other
            .claims()
            .query_typed(
                &property.property_key(request.normal_stress),
                &native,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        // Holding the declared card/state/value fixed cannot hide a source switch.
        assert_eq!(
            answers[property.index()].evidence.value,
            other_answers[property.index()].evidence.value
        );
        assert!(!source.matches_material_state(
            material.content_hash(),
            material.id(),
            &native,
            &other_answers[0],
            &other_answers[1],
            &other_answers[2]
        ));
    }
}
