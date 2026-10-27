//! G1/G3 synthetic constitutive and source-substitution checks; no measured material claim.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use fs_material::piezoelectric::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims};

fn source(s: f64, d: f64, epsilon: f64, density: f64, dynamic: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (COMPLIANCE_33_E_PROPERTY, COMPLIANCE_DIMS, s),
        (STRAIN_CHARGE_33_PROPERTY, STRAIN_CHARGE_DIMS, d),
        (
            PERMITTIVITY_33_T_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            epsilon,
        ),
        ("density", Density::DIMS, density),
    ] {
        let validity = ValidityDomain::unconstrained().with("T", 250.0, 400.0);
        let validity = if dynamic {
            validity.with("electric_field", -1.0, 1.0)
        } else {
            validity
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "authored longitudinal piezoelectric fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal solid".into(),
            phase: "solid".into(),
            process: "authored axis-3 reduction".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point() -> QueryPoint {
    QueryPoint::new().with("T", 300.0).unwrap()
}
fn resolve(card: &MaterialCard) -> LongitudinalPiezoelectric {
    LongitudinalPiezoelectric::resolve(card, &point(), MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance {tolerance}"
    );
}

#[test]
fn g1_direct_converse_response_is_one_reciprocal_energy() {
    for d in [-0.2, 0.0, 0.2] {
        let card = source(0.1, d, 2.0, 3.0, false);
        let law = resolve(&card);
        assert_eq!(law.material().card_identity(), card.content_hash());
        assert_eq!(law.material().query_point(), &[("T".to_owned(), 300.0)]);
        for stress in [-0.3, 0.0, 0.3] {
            for field in [-0.04, 0.0, 0.02] {
                // Independent strain-charge equations in the original source convention.
                let strain = 0.1 * stress + d * field;
                let displacement = d * stress + 2.0 * field;
                let response = law.response(strain, displacement).unwrap();
                near(response.stress_pa, stress, 2e-15);
                near(response.electric_field_v_m, field, 2e-15);
                near(
                    response.stored_energy_j_m3,
                    0.5 * (strain * stress + displacement * field),
                    2e-15,
                );
                let h = 1e-6;
                let ds = (law
                    .response(strain + h, displacement)
                    .unwrap()
                    .stored_energy_j_m3
                    - law
                        .response(strain - h, displacement)
                        .unwrap()
                        .stored_energy_j_m3)
                    / (2.0 * h);
                let dd = (law
                    .response(strain, displacement + h)
                    .unwrap()
                    .stored_energy_j_m3
                    - law
                        .response(strain, displacement - h)
                        .unwrap()
                        .stored_energy_j_m3)
                    / (2.0 * h);
                near(ds, stress, 2e-11);
                near(dd, field, 2e-11);
            }
        }
        near(
            law.axially_clamped_permittivity_f_m(),
            2.0 - d * d / 0.1,
            1e-15,
        );
        near(law.short_circuit_modulus_pa(), 10.0, 1e-15);
        near(
            law.open_circuit_modulus_pa(),
            1.0 / (0.1 - d * d / 2.0),
            1e-13,
        );
    }
}

#[test]
fn g1_si_sized_coefficients_keep_compliance_charge_and_absolute_permittivity_units() {
    // Plausible magnitudes, explicitly synthetic values rather than a named PZT dataset.
    let law = resolve(&source(2e-11, 100e-12, 1e-8, 7500.0, false));
    let response = law.response(2.1e-5, 2e-4).unwrap();
    near(response.stress_pa, 1e6, 1e-8);
    near(response.electric_field_v_m, 1e4, 1e-10);
    near(response.stored_energy_j_m3, 11.5, 1e-12);
    assert_eq!(COMPLIANCE_DIMS, Dims([1, -1, 2, 0, 0, 0]));
    assert_eq!(STRAIN_CHARGE_DIMS, Dims([-1, -1, 3, 0, 1, 0]));
}

#[test]
fn g0_whole_coupled_block_source_domain_and_dynamic_roles_are_enforced() {
    for (s, d, epsilon, rho) in [
        (0.1, 0.5, 2.0, 3.0),
        (0.0, 0.2, 2.0, 3.0),
        (0.1, 0.2, 2.0, 0.0),
    ] {
        assert!(
            LongitudinalPiezoelectric::resolve(
                &source(s, d, epsilon, rho, false),
                &point(),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let card = source(0.1, 0.2, 2.0, 3.0, false);
    assert!(
        LongitudinalPiezoelectric::resolve(
            &card,
            &point().with("T", 401.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let dynamic = source(0.1, 0.2, 2.0, 3.0, true);
    assert!(
        LongitudinalPiezoelectric::resolve(
            &dynamic,
            &point().with("electric_field", 0.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let law = resolve(&card);
    assert!(law.response(f64::NAN, 0.0).is_err());
    assert!(law.response(0.0, f64::MAX).is_err());
}
