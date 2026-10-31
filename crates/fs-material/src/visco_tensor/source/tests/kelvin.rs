//! Manufactured material data; no measured viscosity is invented here.
use super::*;
use crate::visco_tensor::TensorProjectionBudget;
use fs_qty::{DynViscosity, Temperature};

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", t)
        .unwrap()
        .with("omega", 100.0)
        .unwrap()
}
fn data(
    matrix: [[f64; 6]; 6],
) -> (
    TensorPronyProperties,
    TensorKelvinViscosityProperties,
    Vec<PropertyClaim>,
) {
    let (properties, mut claims) = fixture();
    claims.push(scalar("density", Density::DIMS, 1200.0));
    let kelvin = TensorKelvinViscosityProperties {
        viscosity_upper: core::array::from_fn(|i| format!("kelvin_{i:02}")),
    };
    let mut entry = 0;
    for i in 0..6 {
        for j in i..6 {
            claims.push(scalar(
                &kelvin.viscosity_upper[entry],
                DynViscosity::DIMS,
                matrix[i][j],
            ));
            entry += 1;
        }
    }
    (properties, kelvin, claims)
}
fn bound(
    card: &MaterialCard,
    properties: TensorPronyProperties,
    kelvin: TensorKelvinViscosityProperties,
) -> ResolvedTensorMaterial {
    ResolvedTensorMaterial::resolve(
        card,
        &point(300.0),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .with_kelvin_viscosity(card, kelvin, MaterialPropertySelection::SingleClaimOnly)
    .unwrap()
}
fn curve(claims: &mut [PropertyClaim], name: &str, knots: Vec<(f64, f64)>) {
    let c = claims.iter_mut().find(|c| c.key.name() == name).unwrap();
    c.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots,
        dims: DynViscosity::DIMS,
    };
    c.interpolation = InterpolationPolicy::LinearInside;
}

#[test]
fn full_original_kelvin_source_survives_frequency_rotation_and_physical_shear_projection() {
    let factor = [0.0, 0.0, 0.0, 2.0, 0.0, -1.0];
    let matrix = core::array::from_fn(|i| core::array::from_fn(|j| factor[i] * factor[j]));
    let (properties, kelvin, mut claims) = data(matrix);
    claims
        .iter_mut()
        .find(|c| c.key.name() == kelvin.viscosity_upper[0])
        .unwrap()
        .validity = ValidityDomain::unconstrained()
        .with("T", 273.0, 373.0)
        .with("omega", 0.0, 700.0);
    let card = card(claims);
    let dry = ResolvedTensorMaterial::resolve(
        &card,
        &point(300.0),
        properties.clone(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let material = bound(&card, properties, kelvin);
    assert_eq!(
        material.kelvin_source().unwrap().state().properties().len(),
        21
    );
    assert_eq!(material.law().kv_pa_s(), Some(&matrix));
    assert_eq!(
        material.law().maximum_strain_norm(),
        dry.law().maximum_strain_norm()
    );
    assert_eq!(material.law().band_rad_s(), (0.0, 700.0));
    let (storage, loss) = material.law().modulus_mandel(100.0).unwrap();
    let (dry_storage, dry_loss) = dry.law().modulus_mandel(100.0).unwrap();
    assert_eq!(storage, dry_storage);
    for i in 0..6 {
        for j in 0..6 {
            assert!((loss[i][j] - dry_loss[i][j] - 100.0 * matrix[i][j]).abs() < 1e-7);
        }
    }
    let a = 0.63_f64;
    let rotation = [
        [a.cos(), -a.sin(), 0.0],
        [a.sin(), a.cos(), 0.0],
        [0.0, 0.0, 1.0],
    ];
    let rotated = material.rotated_law(rotation).unwrap();
    assert_eq!(
        rotated.kv_pa_s(),
        Some(&crate::tensor::rotate_mandel(matrix, rotation))
    );
    let h = 0.001;
    let columns = [
        [0.0, 0.0, 0.0, 0.5 / h, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.5 / h],
    ];
    let projected = material
        .law()
        .project_strain(
            &columns,
            0.02 * h,
            TensorProjectionBudget {
                max_coordinates: 2,
                max_scalar_terms: 2000,
            },
            || false,
        )
        .unwrap();
    let rows = projected.viscosity_gram_rows().unwrap();
    assert_eq!(rows.len(), 6);
    for i in 0..2 {
        for j in 0..2 {
            let cij = rows.iter().map(|r| r[i] * r[j]).sum::<f64>();
            let expected = 0.02 / (2.0 * h) * matrix[[3, 5][i]][[3, 5][j]];
            assert!((cij - expected).abs() < 1e-12);
            assert_eq!(cij, projected.viscosity_matrix()[2 * i + j]);
        }
    }
}

#[test]
fn piecewise_linear_kelvin_refresh_preserves_original_receipts_and_checks_the_whole_interval() {
    let matrix = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 4.0 } else { 0.0 }));
    let (properties, kelvin, mut claims) = data(matrix);
    curve(
        &mut claims,
        &kelvin.viscosity_upper[15],
        vec![(273.0, 3.0), (300.0, 4.0), (333.0, 6.0), (373.0, 8.0)],
    );
    curve(
        &mut claims,
        &kelvin.viscosity_upper[17],
        vec![(273.0, -0.3), (310.0, 0.2), (373.0, 0.5)],
    );
    let card = card(claims);
    let material = bound(&card, properties, kelvin);
    material
        .admit_thermal_relaxation(&card, "T", [273.0, 373.0])
        .unwrap();
    let changed = material.at_query_point(&card, &point(353.0)).unwrap();
    assert_eq!(changed.law().kv_pa_s().unwrap()[3][3], 7.0);
    assert_eq!(material.law().kv_pa_s().unwrap()[3][3], 4.0);
    assert_eq!(changed.law().c_eq(), material.law().c_eq());
    for (a, b) in material
        .kelvin_source()
        .unwrap()
        .state()
        .properties()
        .iter()
        .zip(changed.kelvin_source().unwrap().state().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert_eq!(
        changed.kelvin_source().unwrap().state().query_point(),
        changed.state().query_point()
    );
    assert!(material.at_query_point(&card, &point(374.0)).is_err());
}

#[test]
fn an_indefinite_interior_kelvin_knot_cannot_hide_between_positive_endpoints() {
    let matrix = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    let (properties, kelvin, mut claims) = data(matrix);
    curve(
        &mut claims,
        &kelvin.viscosity_upper[17],
        vec![(273.0, 0.0), (300.0, 0.0), (326.0, 2.0), (373.0, 0.0)],
    );
    let card = card(claims);
    let material = bound(&card, properties, kelvin);
    assert!(material.at_query_point(&card, &point(273.0)).is_ok());
    assert!(material.at_query_point(&card, &point(373.0)).is_ok());
    assert!(
        material
            .admit_thermal_relaxation(&card, "T", [273.0, 373.0])
            .is_err()
    );
    assert!(thermal::AffineTensorThermal::try_new(&material, &card, "T", [273.0, 373.0]).is_err());
}

#[test]
fn missing_kelvin_coefficient_wrong_units_or_foreign_card_refuse() {
    let matrix = core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    let (properties, kelvin, claims) = data(matrix);
    let full = card(claims.clone());
    let material = bound(&full, properties.clone(), kelvin.clone());
    assert!(
        material
            .admit_kelvin_thermal_interval(&full, "omega", [50.0, 200.0])
            .is_err()
    );
    let mut missing = claims.clone();
    missing.retain(|c| c.key.name() != kelvin.viscosity_upper[20]);
    let incomplete = card(missing);
    let dry = ResolvedTensorMaterial::resolve(
        &incomplete,
        &point(300.0),
        properties.clone(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        dry.with_kelvin_viscosity(
            &incomplete,
            kelvin.clone(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let mut wrong = claims;
    let c = wrong
        .iter_mut()
        .find(|c| c.key.name() == kelvin.viscosity_upper[20])
        .unwrap();
    *c = scalar(&kelvin.viscosity_upper[20], Pressure::DIMS, 1.0);
    let bad = card(wrong);
    let dry = ResolvedTensorMaterial::resolve(
        &bad,
        &point(300.0),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        dry.with_kelvin_viscosity(&bad, kelvin, MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    assert!(material.at_query_point(&bad, &point(310.0)).is_err());
    assert_eq!(material.state().state_coordinate("T"), Some(300.0));
}
