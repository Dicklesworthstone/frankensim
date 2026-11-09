//! Original material-card fixtures plus an independent two-species WMS oracle.
use super::*;
use crate::thermochemical_gas::{SpeciesTransportMaterial, tests as fixture, typed};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::semantic::QuantityKind;
use fs_qty::Dims;

const MU_A: &[(f64, f64)] = &[(200.0, 1e-5), (400.0, 2e-5), (800.0, 7e-5)];
const MU_B: &[(f64, f64)] = &[(200.0, 6e-5), (650.0, 2e-5), (800.0, 1.5e-5)];
const K_A: &[(f64, f64)] = &[(200.0, 0.02), (350.0, 0.06), (800.0, 0.09)];
const K_B: &[(f64, f64)] = &[(200.0, 0.25), (500.0, 0.13), (800.0, 0.10)];

fn claim(name: &str, dims: Dims, knots: &[(f64, f64)]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: knots.to_vec(),
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                knots[0].0,
                knots[knots.len() - 1].0,
            )
            .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured original conductivity paths".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn transport(name: &str, claims: Vec<PropertyClaim>) -> SpeciesTransportMaterial {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "manufactured".into(),
                revision: 0,
            },
            set,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn bank() -> ThermochemicalGasModel {
    let original = fixture::model(vec![0.6, 0.4], 0.001);
    ThermochemicalGasModel::new(
        original.thermodynamics().clone(),
        vec![
            transport(
                "A",
                vec![
                    claim("dynamic_viscosity", DynViscosity::DIMS, MU_A),
                    claim(
                        THERMAL_CONDUCTIVITY_PROPERTY,
                        THERMAL_CONDUCTIVITY_DIMS,
                        K_A,
                    ),
                ],
            ),
            transport(
                "B",
                vec![
                    claim("dynamic_viscosity", DynViscosity::DIMS, MU_B),
                    claim(
                        THERMAL_CONDUCTIVITY_PROPERTY,
                        THERMAL_CONDUCTIVITY_DIMS,
                        K_B,
                    ),
                ],
            ),
        ],
        None,
    )
    .unwrap()
}

fn controls() -> GasConductivityQuadrature {
    GasConductivityQuadrature {
        absolute_tolerance_w_m_k: 1e-13,
        relative_tolerance: 1e-11,
        maximum_subdivisions: 2048,
        maximum_depth: 24,
    }
}

fn interpolate(knots: &[(f64, f64)], t: f64) -> f64 {
    if t == knots[knots.len() - 1].0 {
        return knots[knots.len() - 1].1;
    }
    let i = knots.partition_point(|&(x, _)| x <= t) - 1;
    let fraction = (t - knots[i].0) / (knots[i + 1].0 - knots[i].0);
    knots[i].1 * (1.0 - fraction) + knots[i + 1].1 * fraction
}

fn oracle(t: f64) -> f64 {
    let mu = [interpolate(MU_A, t), interpolate(MU_B, t)];
    let k = [interpolate(K_A, t), interpolate(K_B, t)];
    let mass: [f64; 2] = [0.028, 0.004];
    let fraction = [0.6, 0.4];
    let phi = |i: usize, j: usize| {
        let numerator = 1.0 + (mu[i] / mu[j]).sqrt() * (mass[j] / mass[i]).sqrt().sqrt();
        numerator * numerator / (8.0 * (1.0 + mass[i] / mass[j])).sqrt()
    };
    fraction[0] * k[0] / (fraction[0] + fraction[1] * phi(0, 1))
        + fraction[1] * k[1] / (fraction[1] + fraction[0] * phi(1, 0))
}

#[test]
fn point_values_keep_original_wilke_wms_sources_context_and_composition() {
    let bank = bank();
    let context = fixture::point(300.0, 101325.0);
    let law = PreparedThermochemicalConductivity::new(&bank, &context, [200.0, 800.0], || Ok(()))
        .unwrap();
    assert_eq!(
        law.temperatures_k(),
        &[200.0, 350.0, 400.0, 500.0, 650.0, 800.0]
    );
    assert_eq!(law.context(), &context);
    assert_eq!(
        law.mole_composition(),
        bank.thermodynamics().mole_composition()
    );
    for (line, original) in law.components().iter().zip(bank.transport_materials()) {
        assert_eq!(line.species, original.species);
        assert_eq!(line.material_identity, original.material.content_hash());
    }
    for t in [200.0, 275.0, 350.0, 400.0, 523.0, 650.0, 799.0, 800.0] {
        let original = bank.evaluate(&fixture::point(t, 101325.0)).unwrap();
        let mut calls = 0;
        let k = law
            .conductivity_at_temperature(t, || {
                calls += 1;
                Ok(())
            })
            .unwrap();
        fixture::close(k, original.state().thermal_conductivity);
        fixture::close(k, oracle(t));
        assert_eq!(calls, 1);
        for (source, material) in law.components().iter().zip(original.components()) {
            assert_eq!(
                source.viscosity.knots[0].1.receipt.selected,
                material
                    .property("dynamic_viscosity")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert_eq!(
                source.conductivity.knots[0].1.receipt.selected,
                material
                    .property(THERMAL_CONDUCTIVITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
}

#[test]
fn nonlinear_mean_crosses_every_original_knot_and_differs_from_midpoint() {
    let law = PreparedThermochemicalConductivity::new(
        &bank(),
        &fixture::point(300.0, 101325.0),
        [200.0, 800.0],
        || Ok(()),
    )
    .unwrap();
    let mut calls = 0;
    let result = law
        .mean_conductivity(200.0, 800.0, controls(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert_eq!(calls, result.point_evaluations);
    assert_eq!(result.source_segments, 5);
    // Independent composite Simpson integration of explicit two-species
    // equations, split at the authored pure-species knots.
    let mut integral = 0.0;
    for pair in law.temperatures_k().windows(2) {
        let step = (pair[1] - pair[0]) / 2048.0;
        let mut sum = oracle(pair[0]) + oracle(pair[1]);
        for i in 1..2048 {
            sum += if i % 2 == 0 { 2.0 } else { 4.0 } * oracle(pair[0] + i as f64 * step);
        }
        integral += step * sum / 3.0;
    }
    assert!((result.mean_w_m_k - integral / 600.0).abs() < 2e-11);
    assert!((result.mean_w_m_k - oracle(500.0)).abs() > 1e-3);
    assert_eq!(
        result,
        law.mean_conductivity(800.0, 200.0, controls(), || Ok(()))
            .unwrap()
    );
    let equal = law
        .mean_conductivity(523.0, 523.0, controls(), || Ok(()))
        .unwrap();
    assert_eq!(
        equal.mean_w_m_k,
        law.conductivity_at_temperature(523.0, || Ok(())).unwrap()
    );
    assert_eq!(equal.estimated_absolute_error_w_m_k, 0.0);
    assert_eq!(equal.point_evaluations, 1);
}

#[test]
fn source_holes_nonpositive_interior_and_wrong_context_refuse() {
    let original = fixture::model(vec![0.5, 0.5], 0.001)
        .with_mole_amounts(&[fs_qty::Amount::new(1.0), fs_qty::Amount::new(0.0)])
        .unwrap();
    let make = |mu: Vec<PropertyClaim>| {
        let mut claims = mu;
        claims.push(claim(
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            &[(200.0, 0.02), (800.0, 0.08)],
        ));
        ThermochemicalGasModel::new(
            original.thermodynamics().clone(),
            vec![transport("A", claims)],
            None,
        )
        .unwrap()
    };
    let gap = make(vec![
        claim(
            "dynamic_viscosity",
            DynViscosity::DIMS,
            &[(200.0, 1e-5), (400.0, 2e-5)],
        ),
        claim(
            "dynamic_viscosity",
            DynViscosity::DIMS,
            &[(600.0, 3e-5), (800.0, 4e-5)],
        ),
    ]);
    let context = fixture::point(300.0, 101325.0);
    assert!(gap.evaluate(&fixture::point(200.0, 101325.0)).is_ok());
    assert!(gap.evaluate(&fixture::point(800.0, 101325.0)).is_ok());
    assert!(
        PreparedThermochemicalConductivity::new(&gap, &context, [200.0, 800.0], || Ok(())).is_err()
    );
    let zero = make(vec![claim(
        "dynamic_viscosity",
        DynViscosity::DIMS,
        &[(200.0, 1e-5), (400.0, 0.0), (800.0, 4e-5)],
    )]);
    assert!(
        PreparedThermochemicalConductivity::new(&zero, &context, [200.0, 800.0], || Ok(()))
            .is_err()
    );
    let untyped = QueryPoint::new()
        .with("temperature", 300.0)
        .unwrap()
        .with("pressure", 101325.0)
        .unwrap();
    assert!(
        PreparedThermochemicalConductivity::new(&bank(), &untyped, [200.0, 800.0], || Ok(()))
            .is_err()
    );
    assert!(
        PreparedThermochemicalConductivity::new(
            &bank(),
            &fixture::point(300.0, 1e8),
            [200.0, 800.0],
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn actual_query_cancellation_work_limits_and_domains_are_explicit() {
    let model = bank();
    let context = fixture::point(300.0, 101325.0);
    assert!(matches!(
        PreparedThermochemicalConductivity::new(&model, &context, [200.0, 800.0], || Err(
            GasConductivityError::Cancelled
        ),),
        Err(GasConductivityError::Cancelled)
    ));
    let law = PreparedThermochemicalConductivity::new(&model, &context, [200.0, 800.0], || Ok(()))
        .unwrap();
    assert!(matches!(
        law.conductivity_at_temperature(300.0, || Err(GasConductivityError::WorkLimit)),
        Err(GasConductivityError::WorkLimit)
    ));
    let mut calls = 0;
    assert!(matches!(
        law.mean_conductivity(200.0, 800.0, controls(), || {
            calls += 1;
            if calls == 4 {
                Err(GasConductivityError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(GasConductivityError::Cancelled)
    ));
    assert_eq!(calls, 4);
    for t in [199.0, 801.0, f64::NAN, f64::INFINITY] {
        assert!(law.conductivity_at_temperature(t, || Ok(())).is_err());
        assert!(
            law.mean_conductivity(300.0, t, controls(), || Ok(()))
                .is_err()
        );
    }
    assert!(
        law.mean_conductivity(
            300.0,
            400.0,
            GasConductivityQuadrature {
                relative_tolerance: f64::NAN,
                ..controls()
            },
            || Ok(())
        )
        .is_err()
    );
}
