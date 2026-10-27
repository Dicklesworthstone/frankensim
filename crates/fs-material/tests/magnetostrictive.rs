//! Synthetic reversible tangent checks, not measured magnetostrictive validation.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::MAGNETIC_PERMEABILITY_DIMS;
use fs_material::magnetostrictive::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, QuantitySpec};

fn card(s: f64, d: f64, mu: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_H_PROPERTY, COMPLIANCE_DIMS, s),
        (PIEZOMAGNETIC_33_PROPERTY, PIEZOMAGNETIC_DIMS, d),
        (PERMEABILITY_33_T_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
        ("density", Density::DIMS, 3.0),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 250.0, 400.0)
                    .with_quantity(
                        "bias_H",
                        QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS),
                        -20.0,
                        20.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic reciprocal magnetic tangent".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic biased solid".into(),
            phase: "solid".into(),
            process: "axis-3 reversible tangent".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(h: f64, t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with_quantity("bias_H", QuantitySpec::dimensional(MAGNETIC_FIELD_DIMS), h)
        .unwrap()
}
fn resolve(c: &MaterialCard) -> LongitudinalMagnetostrictive {
    LongitudinalMagnetostrictive::resolve(
        c,
        &point(10.0, 300.0),
        "bias_H",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g1_direct_converse_and_energy_follow_independent_symmetric_block() {
    for d in [-0.2, 0.0, 0.2] {
        let c = card(0.1, d, 2.0);
        let law = resolve(&c);
        assert_eq!(law.material().card_identity(), c.content_hash());
        assert_eq!(law.bias_field_a_m(), 10.0);
        for stress in [-0.3, 0.0, 0.3] {
            for field in [-0.04, 0.0, 0.02] {
                let strain = 0.1 * stress + d * field;
                let b = d * stress + 2.0 * field;
                let response = law.response(strain, b).unwrap();
                near(response.stress_increment_pa, stress, 3e-15);
                near(response.field_increment_a_m, field, 3e-15);
                let determinant = 0.1 * 2.0 - d * d;
                let energy = (2.0 * strain * strain - 2.0 * d * strain * b + 0.1 * b * b)
                    / (2.0 * determinant);
                near(response.incremental_energy_j_m3, energy, 3e-15);
                near(
                    response.incremental_energy_j_m3,
                    0.5 * (stress * strain + field * b),
                    3e-15,
                );
                let h = 1e-6;
                let de = (law.response(strain + h, b).unwrap().incremental_energy_j_m3
                    - law.response(strain - h, b).unwrap().incremental_energy_j_m3)
                    / (2.0 * h);
                let db = (law.response(strain, b + h).unwrap().incremental_energy_j_m3
                    - law.response(strain, b - h).unwrap().incremental_energy_j_m3)
                    / (2.0 * h);
                near(de, stress, 3e-11);
                near(db, field, 3e-11);
            }
        }
        near(law.modulus_at_fixed_field_pa(), 10.0, 1e-14);
        near(law.modulus_at_fixed_flux_pa(), 2.0 / (0.2 - d * d), 1e-13);
    }
}

#[test]
fn g3_bias_temperature_and_signed_material_data_select_actual_tangents() {
    let original = card(0.1, 0.2, 2.0);
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut source = source.clone();
        if source.key.name() == PIEZOMAGNETIC_33_PROPERTY {
            source.value = PropertyValue::Curve {
                abscissa: "bias_H".into(),
                abscissa_dims: MAGNETIC_FIELD_DIMS,
                knots: vec![(-20.0, -0.3), (20.0, 0.3)],
                dims: PIEZOMAGNETIC_DIMS,
            };
            source.interpolation = InterpolationPolicy::LinearInside;
        }
        if source.key.name() == COMPLIANCE_33_H_PROPERTY {
            source.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(250.0, 0.05), (400.0, 0.2)],
                dims: COMPLIANCE_DIMS,
            };
            source.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(source).unwrap();
    }
    let c = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    for bias in [-10.0, 10.0] {
        for temperature in [300.0, 350.0] {
            let law = LongitudinalMagnetostrictive::resolve(
                &c,
                &point(bias, temperature),
                "bias_H",
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            near(law.piezomagnetic_coefficient_m_a(), 0.015 * bias, 1e-15);
            near(
                law.compliance_at_fixed_field_pa_inv(),
                0.05 + 0.001 * (temperature - 250.0),
                1e-15,
            );
            assert_eq!(law.material().card_identity(), c.content_hash());
            assert_eq!(law.material().query().axes()["bias_H"], bias);
        }
    }
}

#[test]
fn g0_instability_unowned_conditions_and_missing_bias_refuse() {
    for (s, d, mu) in [
        (0.1, 0.5, 2.0),
        (0.125, 0.5, 2.0),
        (-0.1, 0.2, 2.0),
        (0.1, 0.2, -2.0),
    ] {
        assert!(
            LongitudinalMagnetostrictive::resolve(
                &card(s, d, mu),
                &point(10.0, 300.0),
                "bias_H",
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let c = card(0.1, 0.2, 2.0);
    for p in [
        point(0.0, 300.0),
        point(30.0, 300.0),
        point(10.0, 500.0),
        QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("bias_H", 10.0)
            .unwrap(),
    ] {
        assert!(
            LongitudinalMagnetostrictive::resolve(
                &c,
                &p,
                "bias_H",
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let mut claims = ClaimSet::new();
    for (_, source) in c.claims().claims_ordered() {
        let mut source = source.clone();
        source.validity = source.validity.with("strain_rate", 0.0, 1.0);
        claims.insert_claim(source).unwrap();
    }
    let dynamic = MaterialCard::assemble(c.id().clone(), claims, vec![]).unwrap();
    assert!(
        LongitudinalMagnetostrictive::resolve(
            &dynamic,
            &point(10.0, 300.0).with("strain_rate", 0.5).unwrap(),
            "bias_H",
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let law = resolve(&c);
    assert!(law.response(f64::NAN, 0.0).is_err());
    assert!(law.response(f64::MAX, 1.0).is_err());
}
