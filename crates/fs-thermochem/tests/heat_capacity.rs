//! G0/G3/G5 public-API checks for the source-owned heat-capacity queries.

use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Temperature};
use fs_thermochem::mixture::{
    FrozenIdealGasMixtureErrorV1, FrozenIdealGasMixtureModelV1, FrozenMixturePropertyV1,
};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, MolarEnergyQuantityV1,
    MolarThermalQuantityV1, NASA9_LAW_ID_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, required_nasa9_parameter_dims,
};
use std::collections::BTreeMap;

fn source(name: &str, mass: f64, regions: &[(f64, f64, [f64; 9])]) -> Nasa9StandardStateModelV1 {
    let regions = regions
        .iter()
        .map(|&(low, high, coefficients)| {
            let mut parameters = BTreeMap::new();
            for ((name, dims), value) in required_nasa9_parameter_dims()
                .into_iter()
                .take(9)
                .zip(coefficients)
            {
                parameters.insert(name.into(), LawParameter { value, dims });
            }
            parameters.insert(
                "reference_pressure".into(),
                LawParameter {
                    value: 100000.0,
                    dims: Pressure::DIMS,
                },
            );
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId(NASA9_LAW_ID_V1.into()),
                law_version: 1,
                parameters,
                state_schema_version: 0,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", low, high),
                sources: Vec::new(),
                provenance: Provenance {
                    source: "analytic heat-capacity fixture".into(),
                    license: "test fixture".into(),
                    artifact: None,
                },
            })
            .unwrap()
        })
        .collect();
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("test-elements").unwrap(),
        )
        .unwrap(),
        regions,
    )
    .unwrap()
}

fn affine(name: &str, mass: f64, cp_over_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_over_r;
    coefficients[3] = slope;
    source(name, mass, &[(200.0, 6000.0, coefficients)])
}

#[test]
fn heat_capacity_keeps_original_polynomial_bits_and_reconciled_source_seams() {
    let low = [80.0, -12.0, 3.5, 1e-3, -1e-8, 3e-12, -1e-16, 100.0, 3.0];
    let mut high = low;
    high[3] = 1.5e-3;
    let raw = source("N2", 0.028, &[(200.0, 1000.0, low), (1000.0, 6000.0, high)]);
    let corrected = raw
        .with_caloric_continuity(
            Temperature::new(300.0),
            MolarEnergyQuantityV1::new(1e7),
            MolarThermalQuantityV1::new(1e5),
        )
        .unwrap();
    for model in [&raw, &corrected] {
        for t in [
            200.0,
            300.0,
            1000.0_f64.next_down(),
            1000.0,
            1000.0_f64.next_up(),
            2500.0,
            6000.0,
        ] {
            let capacity = model.evaluate_heat_capacity(Temperature::new(t)).unwrap();
            let full = model.evaluate(Temperature::new(t)).unwrap();
            assert_eq!(
                capacity.value().to_bits(),
                full.properties().cp().value().to_bits()
            );
            assert_eq!(
                capacity.value().to_bits(),
                raw.evaluate_heat_capacity(Temperature::new(t))
                    .unwrap()
                    .value()
                    .to_bits()
            );
        }
        for t in [
            0.0,
            -1.0,
            199.0,
            6000.0_f64.next_up(),
            f64::INFINITY,
            f64::NAN,
        ] {
            assert_eq!(
                model
                    .evaluate_heat_capacity(Temperature::new(t))
                    .unwrap_err(),
                model.evaluate(Temperature::new(t)).unwrap_err()
            );
        }
    }
    let gap = source("N2", 0.028, &[(200.0, 500.0, low), (700.0, 1000.0, high)]);
    assert!(gap.evaluate_heat_capacity(Temperature::new(600.0)).is_err());
}

#[test]
fn heat_capacity_polynomial_channels_match_the_nasa9_definition() {
    let weights = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 0.0, 0.0];
    for (index, weight) in weights.into_iter().enumerate() {
        let mut coefficients = [0.0; 9];
        coefficients[index] = 1.0;
        let model = source("N2", 0.028, &[(1.0, 10.0, coefficients)]);
        let capacity = model.evaluate_heat_capacity(Temperature::new(2.0)).unwrap();
        assert_eq!(capacity.value().to_bits(), (R * weight).to_bits());
    }
}

#[test]
fn heat_capacity_mixture_tracks_temperature_composition_and_mass_basis() {
    let components = vec![
        affine("O2", 0.032, 3.0, 2e-3),
        affine("N2", 0.028, 3.5, 1e-3),
    ];
    for basis in [
        CompositionBasis::MoleFraction,
        CompositionBasis::MassFraction,
    ] {
        let mixture = FrozenIdealGasMixtureModelV1::new(
            components.clone(),
            Composition::new(basis, vec![0.25, 0.75]).unwrap(),
        )
        .unwrap();
        for t in [200.0, 298.15, 1000.0, 6000.0] {
            let heat = mixture
                .evaluate_heat_capacities(Temperature::new(t))
                .unwrap();
            let x = mixture.mole_composition().fractions();
            let expected = R * (x[0] * (3.5 + 1e-3 * t) + x[1] * (3.0 + 2e-3 * t));
            assert!((heat.cp().value() / expected - 1.0).abs() < 3e-15);
            for p in [0.001, 100000.0, 1e9] {
                let full = mixture
                    .evaluate(Temperature::new(t), Pressure::new(p))
                    .unwrap();
                let full = full.properties();
                assert_eq!(
                    [
                        heat.cp().value(),
                        heat.cv().value(),
                        heat.mass_specific_cp().value(),
                        heat.mass_specific_cv().value()
                    ]
                    .map(f64::to_bits),
                    [
                        full.cp().value(),
                        full.cv().value(),
                        full.mass_specific_cp().value(),
                        full.mass_specific_cv().value()
                    ]
                    .map(f64::to_bits),
                    "T={t}, p={p}, basis={basis:?}",
                );
            }
        }
        let midpoint = mixture
            .evaluate_heat_capacities(Temperature::new(300.0))
            .unwrap();
        let delta_h = mixture
            .specific_enthalpy_change(Temperature::new(299.0), Temperature::new(301.0))
            .unwrap()
            .value();
        assert!((delta_h / (2.0 * midpoint.mass_specific_cp().value()) - 1.0).abs() < 3e-15);
        for t in [0.0, 199.0, 6001.0, f64::INFINITY, f64::NAN] {
            assert_eq!(
                mixture
                    .evaluate_heat_capacities(Temperature::new(t))
                    .unwrap_err(),
                mixture
                    .evaluate_properties(Temperature::new(t), Pressure::new(100000.0))
                    .unwrap_err()
            );
        }
    }
}

#[test]
fn heat_capacity_queries_keep_narrow_arithmetic_refusals() {
    for offset in [7, 8] {
        let mut coefficients = [0.0; 9];
        coefficients[2] = 3.5;
        coefficients[offset] = f64::MAX;
        let model = source("N2", 0.028, &[(200.0, 6000.0, coefficients)]);
        assert!(model.evaluate_properties(Temperature::new(300.0)).is_err());
        assert_eq!(
            model
                .evaluate_heat_capacity(Temperature::new(300.0))
                .unwrap()
                .value(),
            3.5 * R
        );
    }
    let mut overflow = [0.0; 9];
    overflow[6] = f64::MAX;
    let overflow = source("N2", 0.028, &[(200.0, 6000.0, overflow)]);
    assert_eq!(
        overflow
            .evaluate_heat_capacity(Temperature::new(300.0))
            .unwrap_err(),
        overflow.evaluate(Temperature::new(300.0)).unwrap_err()
    );
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![affine("N2", f64::MIN_POSITIVE, 3.5, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    assert!(matches!(
        mixture.evaluate_heat_capacities(Temperature::new(300.0)),
        Err(FrozenIdealGasMixtureErrorV1::NonFiniteEvaluation {
            property: FrozenMixturePropertyV1::MassSpecificCp,
            ..
        })
    ));
    // A local source derivative is not a positive-Cv interval admission.
    let unstable = FrozenIdealGasMixtureModelV1::new(
        vec![affine("N2", 0.028, 0.5, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    assert!(
        unstable
            .evaluate_heat_capacities(Temperature::new(300.0))
            .unwrap()
            .cv()
            .value()
            < 0.0
    );
}
