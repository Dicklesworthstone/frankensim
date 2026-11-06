//! G0/G3: source-backed magnetic tensors, inverse response and parameter work.
use super::*;
use fs_material::{
    field_storage::spatial::{
        PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES, SpatialMagneticLaw, SpatialPermeability,
    },
    state_point::{
        ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};

const AXES: [[f64; 3]; 3] = [
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
    [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
];

fn material(
    properties: &[(&str, f64, Dims)],
    validity: ValidityDomain,
    point: &QueryPoint,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, value, dims) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic principal magnetic response".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic magnetic material".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn principal(values: [f64; 3]) -> [(&'static str, f64, Dims); 3] {
    std::array::from_fn(|i| {
        (
            PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES[i],
            values[i],
            MAGNETIC_PERMEABILITY_DIMS,
        )
    })
}
fn resolve(values: [f64; 3]) -> SpatialPermeability {
    SpatialMagneticLaw::orthotropic(AXES)
        .unwrap()
        .resolve(&material(
            &principal(values),
            ValidityDomain::unconstrained(),
            &QueryPoint::new(),
        ))
        .unwrap()
}
fn projection(vector: [f64; 3]) -> [f64; 3] {
    AXES.map(|axis| (0..3).map(|j| axis[j] * vector[j]).sum())
}
fn world(vector: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|j| (0..3).map(|i| AXES[i][j] * vector[i]).sum())
}

#[test]
fn isotropic_magnetic_response_retains_the_uniform_inductor_source() {
    let mu = 4e-6;
    let state = material(
        &[(
            MAGNETIC_PERMEABILITY_PROPERTY,
            mu,
            MAGNETIC_PERMEABILITY_DIMS,
        )],
        ValidityDomain::unconstrained(),
        &QueryPoint::new(),
    );
    let coefficient = SpatialMagneticLaw::Isotropic.resolve(&state).unwrap();
    assert_eq!(coefficient.isotropic_h_m(), Some(mu));
    assert_eq!(coefficient.maximum_permeability_h_m(), mu);
    assert_eq!(
        coefficient,
        SpatialPermeability::isotropic(linear_magnetic_permeability(&state).unwrap()).unwrap()
    );
    let h = [200.0, -300.0, 500.0];
    let b = coefficient.flux_density_t(h);
    assert_eq!(b, h.map(|value| mu * value));
    assert_eq!(
        coefficient.field_strength_a_m(b).unwrap(),
        b.map(|value| value / mu)
    );
    near(
        coefficient.reluctivity_matrix_m_h().unwrap()[1][1],
        1.0 / mu,
    );
}

#[test]
fn rotated_principal_sources_give_full_flux_energy_and_inverse_response() {
    let values = [2e-6, 5e-6, 13e-6];
    let coefficient = resolve(values);
    let law = SpatialMagneticLaw::orthotropic(AXES).unwrap();
    assert_eq!(law.principal_axes(), Some(AXES));
    assert_eq!(
        law.source_property_names(),
        &PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES
    );
    assert_eq!(coefficient.isotropic_h_m(), None);
    assert_eq!(coefficient.maximum_permeability_h_m(), 13e-6);
    let h = [200.0, -300.0, 500.0];
    let hp = projection(h);
    let expected_b = world(std::array::from_fn(|i| values[i] * hp[i]));
    let b = coefficient.flux_density_t(h);
    for i in 0..3 {
        near(b[i], expected_b[i]);
    }
    let energy = 0.5 * (0..3).map(|i| h[i] * b[i]).sum::<f64>();
    near(
        energy,
        0.5 * (0..3).map(|i| values[i] * hp[i].powi(2)).sum::<f64>(),
    );
    let inverse = coefficient.reluctivity_matrix_m_h().unwrap();
    let arbitrary_b = [0.2, -0.3, 0.5];
    let bp = projection(arbitrary_b);
    let expected_h = world(std::array::from_fn(|i| bp[i] / values[i]));
    let actual_h = coefficient.field_strength_a_m(arbitrary_b).unwrap();
    for i in 0..3 {
        near(actual_h[i], expected_h[i]);
        near(
            (0..3).map(|j| inverse[i][j] * arbitrary_b[j]).sum(),
            expected_h[i],
        );
    }
    let reversed = SpatialMagneticLaw::orthotropic(AXES.map(|row| row.map(|v| -v)))
        .unwrap()
        .resolve(&material(
            &principal(values),
            ValidityDomain::unconstrained(),
            &QueryPoint::new(),
        ))
        .unwrap();
    assert_eq!(reversed.matrix_h_m(), coefficient.matrix_h_m());
    assert_eq!(reversed.reluctivity_matrix_m_h().unwrap(), inverse);
}

#[test]
fn parameter_work_uses_original_principal_changes_and_a_fixed_frame() {
    let before = [2e-6, 5e-6, 13e-6];
    let after = [3e-6, 2e-6, 16e-6];
    let (a, b) = (resolve(before), resolve(after));
    let flux = [0.2, -0.3, 0.5];
    let bp = projection(flux);
    let expected = world(std::array::from_fn(|i| {
        (1.0 / after[i] - 1.0 / before[i]) * bp[i]
    }));
    let actual = a.reluctivity_change_action(b, flux).unwrap();
    for i in 0..3 {
        near(actual[i], expected[i]);
    }
    assert_eq!(a.reluctivity_change_action(a, flux).unwrap(), [0.0; 3]);
    assert_eq!(a.reluctivity_change_action(b, [0.0; 3]).unwrap(), [0.0; 3]);
    assert!(
        a.reluctivity_change_action(SpatialPermeability::isotropic(2e-6).unwrap(), flux)
            .is_err()
    );
    let low = SpatialPermeability::isotropic(1e-300).unwrap();
    let high = SpatialPermeability::isotropic(1e300).unwrap();
    assert_eq!(
        low.reluctivity_change_action(high, [1e-300; 3]).unwrap(),
        [-1.0; 3]
    );
}

#[test]
fn missing_principal_wrong_units_and_unstable_orientations_refuse() {
    let law = SpatialMagneticLaw::orthotropic(AXES).unwrap();
    let mut properties = principal([2e-6, 5e-6, 13e-6]);
    properties[2].0 = MAGNETIC_PERMEABILITY_PROPERTY;
    assert!(
        law.resolve(&material(
            &properties,
            ValidityDomain::unconstrained(),
            &QueryPoint::new()
        ))
        .is_err()
    );
    for (value, dims) in [
        (13e-6, Dims([0; 6])),
        (13e-6, Dims([2, 0, 0, 0, 0, 0])),
        (0.0, MAGNETIC_PERMEABILITY_DIMS),
        (-1.0, MAGNETIC_PERMEABILITY_DIMS),
    ] {
        properties[2] = (PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES[2], value, dims);
        assert!(
            law.resolve(&material(
                &properties,
                ValidityDomain::unconstrained(),
                &QueryPoint::new()
            ))
            .is_err()
        );
    }
    for axes in [
        [[1.0, 0.0, 0.0]; 3],
        [[2.0, 0.0, 0.0], AXES[1], AXES[2]],
        [[f64::NAN, 0.0, 0.0], AXES[1], AXES[2]],
    ] {
        assert!(SpatialMagneticLaw::orthotropic(axes).is_err());
    }
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(SpatialPermeability::isotropic(value).is_err());
    }
    assert!(
        SpatialPermeability::isotropic(f64::from_bits(1))
            .unwrap()
            .reluctivity_matrix_m_h()
            .is_err()
    );
}

#[test]
fn field_and_frequency_bounded_sources_are_not_frozen_into_linear_tensors() {
    let law = SpatialMagneticLaw::orthotropic(AXES).unwrap();
    for (name, quantity) in [
        (
            "frequency",
            static_quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
        ),
        ("B", QuantitySpec::dimensional(Dims([0, 1, -2, 0, -1, 0]))),
    ] {
        let validity = ValidityDomain::unconstrained().with_quantity(name, quantity, 0.0, 1.0);
        let point = QueryPoint::new()
            .with_quantity(name, quantity, 0.5)
            .unwrap();
        assert!(
            law.resolve(&material(&principal([2e-6, 5e-6, 13e-6]), validity, &point))
                .is_err()
        );
    }
}
