//! G0/G3 source admission and independently projected 3D dielectric response.
use super::*;
use fs_material::{
    field_storage::spatial::{
        PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES, SpatialDielectricLaw, SpatialPermittivity,
    },
    state_point::{
        ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};

const AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn material(
    properties: &[(&str, f64, Dims)],
    validity: ValidityDomain,
    query: &QueryPoint,
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
                    source: "manufactured 3D dielectric".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let source = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &source,
        query,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn principal(values: [f64; 3]) -> [(&'static str, f64, Dims); 3] {
    std::array::from_fn(|i| {
        (
            PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[i],
            values[i],
            ELECTRIC_PERMITTIVITY_DIMS,
        )
    })
}

#[test]
fn isotropic_storage_matches_original_uniform_capacitor_source() {
    let source = material(
        &[(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            4e-11,
            ELECTRIC_PERMITTIVITY_DIMS,
        )],
        ValidityDomain::unconstrained(),
        &QueryPoint::new(),
    );
    let coefficient = SpatialDielectricLaw::Isotropic.resolve(&source).unwrap();
    assert_eq!(coefficient, SpatialPermittivity::Scalar(4e-11));
    let uniform = ResolvedCapacitor::try_from_material(source, 0.3, 0.02).unwrap();
    near(
        uniform.capacitance_f(),
        coefficient.matrix_f_m()[0][0] * 0.3 / 0.02,
    );
    for (actual, expected) in coefficient
        .displacement_c_m2([2.0, -3.0, 5.0])
        .into_iter()
        .zip([8e-11, -12e-11, 20e-11])
    {
        near(actual, expected);
    }
}

#[test]
fn three_dimensional_rotation_preserves_principal_response_and_energy() {
    let axes = [
        [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
        [2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0],
        [-2.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0],
    ];
    let values = [2e-11, 5e-11, 13e-11];
    let source = material(
        &principal(values),
        ValidityDomain::unconstrained(),
        &QueryPoint::new(),
    );
    let law = SpatialDielectricLaw::orthotropic(axes).unwrap();
    assert_eq!(law.principal_axes(), Some(axes));
    assert_eq!(
        law.source_property_names(),
        &PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES
    );
    let coefficient = law.resolve(&source).unwrap();
    let matrix = coefficient.matrix_f_m();
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(matrix[i][j], matrix[j][i]);
        }
    }
    // Independent principal-frame projection; no tensor::rotate oracle reuse.
    let field = [2.0, -3.0, 5.0];
    let projection = axes.map(|axis| (0..3).map(|j| axis[j] * field[j]).sum::<f64>());
    let expected: [f64; 3] =
        std::array::from_fn(|j| (0..3).map(|i| values[i] * projection[i] * axes[i][j]).sum());
    let actual = coefficient.displacement_c_m2(field);
    for j in 0..3 {
        near(actual[j], expected[j]);
    }
    let energy = 0.5 * (0..3).map(|j| field[j] * actual[j]).sum::<f64>();
    let expected_energy = 0.5
        * (0..3)
            .map(|i| values[i] * projection[i].powi(2))
            .sum::<f64>();
    near(energy, expected_energy);
    // Reversing a principal direction changes neither the material nor D.
    let reversed = SpatialDielectricLaw::orthotropic(axes.map(|row| row.map(|v| -v))).unwrap();
    assert_eq!(reversed.resolve(&source).unwrap(), coefficient);
}

#[test]
fn orthotropic_law_never_invents_missing_principal_or_accepts_wrong_units() {
    let law = SpatialDielectricLaw::orthotropic(AXES).unwrap();
    let mut properties = principal([2e-11, 5e-11, 13e-11]);
    properties[2] = (
        ELECTRIC_PERMITTIVITY_PROPERTY,
        13e-11,
        ELECTRIC_PERMITTIVITY_DIMS,
    );
    let incomplete = material(
        &properties,
        ValidityDomain::unconstrained(),
        &QueryPoint::new(),
    );
    assert!(law.resolve(&incomplete).is_err());
    for (value, dims) in [
        (13e-11, Dims([0; 6])),
        (0.0, ELECTRIC_PERMITTIVITY_DIMS),
        (-1e-11, ELECTRIC_PERMITTIVITY_DIMS),
    ] {
        properties[2] = (PRINCIPAL_ELECTRIC_PERMITTIVITY_PROPERTIES[2], value, dims);
        let invalid = material(
            &properties,
            ValidityDomain::unconstrained(),
            &QueryPoint::new(),
        );
        assert!(law.resolve(&invalid).is_err());
    }
}

#[test]
fn orientation_and_dynamic_source_applicability_are_not_silently_repaired() {
    for invalid in [
        [[1.0, 0.0, 0.0]; 3],
        [[2.0, 0.0, 0.0], AXES[1], AXES[2]],
        [[f64::NAN, 0.0, 0.0], AXES[1], AXES[2]],
        [[f64::INFINITY, 0.0, 0.0], AXES[1], AXES[2]],
    ] {
        assert!(SpatialDielectricLaw::orthotropic(invalid).is_err());
    }
    let law = SpatialDielectricLaw::orthotropic(AXES).unwrap();
    for (axis, quantity) in [
        (
            "frequency",
            static_quantity(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
        ),
        ("E", QuantitySpec::dimensional(Dims([1, 1, -3, 0, -1, 0]))),
    ] {
        let source = material(
            &principal([2e-11, 5e-11, 13e-11]),
            ValidityDomain::unconstrained().with_quantity(axis, quantity, 0.0, 1.0),
            &QueryPoint::new()
                .with_quantity(axis, quantity, 0.5)
                .unwrap(),
        );
        assert!(law.resolve(&source).is_err());
    }
}
