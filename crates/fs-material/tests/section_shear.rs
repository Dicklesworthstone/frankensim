//! G0/G1/G3: bounded canonical torsion and exact per-layer peak/work recovery.
//! Coefficients are synthetic; analytical references establish no material fit.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Pressure};

fn section(bore: f64, layers: &[(f64, f64)], correction: f64) -> ConcentricShearSection {
    let layers = layers
        .iter()
        .map(|&(radius, shear)| {
            let mut claims = ClaimSet::new();
            let mut requirements = Vec::new();
            for (name, value, dims) in [
                ("density", 1000.0, Density::DIMS),
                ("young_modulus", 2.5 * shear, Pressure::DIMS),
                ("shear_modulus", shear, Pressure::DIMS),
            ] {
                claims
                    .insert_claim(PropertyClaim {
                        key: PropertyKey::new(name, dims),
                        value: PropertyValue::Scalar { value, dims },
                        validity: ValidityDomain::unconstrained(),
                        interpolation: InterpolationPolicy::ConstantWithinValidity,
                        uncertainty: UncertaintyModel::Unstated,
                        observations: vec![],
                        provenance: Provenance {
                            source: "synthetic torsion section".into(),
                            license: "test-only".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
                requirements.push(
                    ScalarPropertyRequirement::try_new(
                        name,
                        dims,
                        ScalarAdmissibility::StrictlyPositive,
                    )
                    .unwrap(),
                );
            }
            let card = MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "synthetic".into(),
                    phase: "solid".into(),
                    process: "torsion test".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap();
            ConcentricMaterialLayer {
                outer_radius_m: radius,
                material: resolve_material_state_point(
                    &card,
                    &QueryPoint::new(),
                    &requirements,
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap(),
            }
        })
        .collect();
    ConcentricElasticSection::try_new(bore, layers)
        .unwrap()
        .with_shear_modulus("shear_modulus", correction)
        .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn bounded_point_observations_delegate_to_the_canonical_torsion_law() {
    let section = section(0.0, &[(0.02, 8.0e10)], 0.8);
    let twist = 0.002;
    for point in [[0.0; 2], [0.012, 0.016], [0.0, -0.01]] {
        let bounded = section.torsional_stress_at_twist_bounded(0, point, twist, 0.01).unwrap();
        assert_eq!(bounded, section.torsional_stress_at_twist(0, point, twist).unwrap());
    }
    // A centerline query has zero local shear but still enforces the bound at R.
    assert!(section.torsional_stress_at_twist(0, [0.0; 2], 1.0).is_ok());
    assert!(section.torsional_stress_at_twist_bounded(0, [0.0; 2], 1.0, 0.01).is_err());
}

#[test]
fn g1_bonded_layers_keep_their_own_peak_and_native_torque_energy() {
    let section = section(0.005, &[(0.01, 8e10), (0.02, 2e10)], 0.4);
    let twist = -0.003;
    let mut torque = 0.0;
    let mut energy = 0.0;
    let mut inner: f64 = 0.005;
    for (layer, (outer, shear)) in [(0.01_f64, 8e10), (0.02_f64, 2e10)].into_iter().enumerate() {
        let actual = section.torsional_layer_response(layer, twist, 0.01).unwrap();
        let polar = core::f64::consts::PI * (outer.powi(4) - inner.powi(4)) / 2.0;
        near(actual.maximum_shear_pa, shear * twist.abs() * outer);
        near(actual.torque_n_m, shear * polar * twist);
        near(
            actual.energy_per_length_j_m,
            0.5 * shear * polar * twist * twist,
        );
        torque += actual.torque_n_m;
        energy += actual.energy_per_length_j_m;
        inner = outer;
    }
    near(torque, section.torsional_stiffness_n_m2() * twist);
    near(
        energy,
        0.5 * section.torsional_stiffness_n_m2() * twist * twist,
    );
}

#[test]
fn g3_layer_peak_storage_and_reversal_ignore_transverse_correction() {
    let a = section(0.0, &[(0.02, 3e10)], 0.2);
    let b = section(0.0, &[(0.02, 3e10)], 1.0);
    for twist in [-0.01, 0.0, 0.01] {
        let forward = a.torsional_layer_response(0, twist, 0.1).unwrap();
        let reverse = a.torsional_layer_response(0, -twist, 0.1).unwrap();
        assert_eq!(forward, b.torsional_layer_response(0, twist, 0.1).unwrap());
        assert_eq!(forward.maximum_shear_pa, reverse.maximum_shear_pa);
        assert_eq!(forward.energy_per_length_j_m, reverse.energy_per_length_j_m);
        assert_eq!(forward.torque_n_m, -reverse.torque_n_m);
    }
    let rest_layer = a.torsional_layer_response(0, 0.0, 0.1).unwrap();
    assert_eq!(rest_layer.maximum_shear_pa, 0.0);
    assert_eq!(rest_layer.torque_n_m, 0.0);
    assert_eq!(rest_layer.energy_per_length_j_m, 0.0);
}

#[test]
fn g0_invalid_layers_points_and_hidden_surface_overstrain_refuse() {
    let section = section(0.005, &[(0.01, 8e10), (0.02, 2e10)], 0.8);
    assert!(section.torsional_stress_at_twist_bounded(2, [0.02, 0.0], 0.001, 0.01).is_err());
    assert!(section.torsional_layer_response(2, 0.001, 0.01).is_err());
    for bound in [0.0, -1.0, 1.0, f64::NAN, f64::INFINITY] {
        assert!(section.torsional_layer_response(0, 0.0, bound).is_err());
    }
    for twist in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1.01] {
        assert!(
            section
                .torsional_stress_at_twist_bounded(0, [0.005, 0.0], twist, 0.01)
                .is_err()
        );
    }
    assert!(
        section
            .torsional_stress_at_twist_bounded(0, [0.005, 0.0], 1.0, 0.01)
            .is_ok()
    );
    assert!(section.torsional_layer_response(0, 1e-200, 0.01).is_err());
}
