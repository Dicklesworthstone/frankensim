//! G1/G3: one original density/bulk/shear bundle drives wave and pressure
//! response. Manufactured temperature curves are implementation fixtures.

use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::acoustic::{
    ElasticAcousticState,
    sourced::{IsotropicBulkShearRequirements, resolve_isotropic_bulk_shear_state},
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), temperature)
        .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
}

fn requirements() -> IsotropicBulkShearRequirements {
    IsotropicBulkShearRequirements::try_new(
        requirement("source-rho", Density::DIMS),
        requirement("source-K", Pressure::DIMS),
        requirement("source-G", Pressure::DIMS),
    )
    .unwrap()
}

fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("source-rho", Density::DIMS, [2000.0, 1800.0]),
        ("source-K", Pressure::DIMS, [30e9, 20e9]),
        ("source-G", Pressure::DIMS, [12e9, 6e9]),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_quantity(),
                    300.0,
                    400.0,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured bulk/shear source, no measurement claim".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "isotropic numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn original_curves_refresh_pressure_and_all_wave_modes_together() {
    let card = card();
    let cold = resolve_isotropic_bulk_shear_state(
        &card,
        &point(300.0),
        &requirements(),
        0.001,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let warm = cold.at_query_point(&card, &point(350.0)).unwrap();
    for (state, rho, bulk, shear) in [(&cold, 2000.0, 30e9, 12e9), (&warm, 1900.0, 25e9, 9e9)] {
        assert_eq!(state.material().card_identity(), card.content_hash());
        assert_eq!(state.material().properties().len(), 3);
        assert_eq!(state.requirements(), &requirements());
        for property in state.material().properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
        near(state.acoustic().bulk_modulus_pa(), bulk);
        let pressure = state.acoustic().hydrostatic_response(1e6).unwrap();
        near(pressure.volumetric_strain, -1e6 / bulk);
        near(pressure.stored_energy_density_j_m3, 1e12 / (2.0 * bulk));
        for direction in [[1.0, 0.0, 0.0], [2.0, -3.0, 1.0]] {
            let modes = state.acoustic().directional_waves(direction).unwrap().modes;
            near(modes[0].phase_speed_m_s, (shear / rho).sqrt());
            near(modes[1].phase_speed_m_s, (shear / rho).sqrt());
            near(
                modes[2].phase_speed_m_s,
                ((bulk + 4.0 * shear / 3.0) / rho).sqrt(),
            );
        }
    }
    assert_ne!(cold.identity(), warm.identity());
    let original = cold.clone();
    assert!(cold.at_query_point(&card, &point(401.0)).is_err());
    assert!(
        cold.at_query_point(&card, &QueryPoint::new().with("other-T", 350.0).unwrap())
            .is_err()
    );
    assert_eq!(cold, original);
    assert!(cold.acoustic().hydrostatic_response(1e9).is_err());
}

#[test]
fn parameter_roles_and_admissibility_are_explicit() {
    assert!(
        IsotropicBulkShearRequirements::try_new(
            requirement("rho", Density::DIMS),
            requirement("same", Pressure::DIMS),
            requirement("same", Pressure::DIMS),
        )
        .is_err()
    );
    assert!(
        IsotropicBulkShearRequirements::try_new(
            requirement("rho", Pressure::DIMS),
            requirement("K", Pressure::DIMS),
            requirement("G", Pressure::DIMS),
        )
        .is_err()
    );
    let card = card();
    assert!(
        resolve_isotropic_bulk_shear_state(
            &card,
            &point(300.0),
            &IsotropicBulkShearRequirements::canonical().unwrap(),
            0.001,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    ); // no silent source-key aliasing
    for (bulk, shear) in [
        (0.0, 1e9),
        (1e9, -1.0),
        (f64::INFINITY, 1e9),
        (1e9, f64::NAN),
    ] {
        assert!(
            ElasticAcousticState::from_bulk_and_shear(
                2000.0,
                bulk,
                shear,
                ContentHash([7; 32]),
                0.001,
            )
            .is_err()
        );
    }
}

#[test]
fn positive_moduli_support_auxetic_media_without_a_poisson_shortcut() {
    // K < 2G/3 gives negative lambda and negative Poisson ratio, while the
    // complete isotropic law is still positive definite.
    let state =
        ElasticAcousticState::from_bulk_and_shear(1000.0, 1e9, 3e9, ContentHash([8; 32]), 0.001)
            .unwrap();
    near(state.bulk_modulus_pa(), 1e9);
    let modes = state.directional_waves([1.0, 1.0, 1.0]).unwrap().modes;
    near(modes[0].directional_modulus_pa, 3e9);
    near(modes[1].directional_modulus_pa, 3e9);
    near(modes[2].directional_modulus_pa, 5e9);
}
