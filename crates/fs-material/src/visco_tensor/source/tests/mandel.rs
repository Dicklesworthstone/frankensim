//! Full measured-equilibrium tensors reach the original hereditary/thermal owners.
use super::*;
use crate::tensor;
use crate::visco_tensor::TensorProjectionBudget;
use crate::{Tangent6, Voigt};

// A strictly diagonally dominant anisotropic tensor with signed normal/shear
// and shear/shear couplings. It does not commute with the Maxwell branch.
const EQUILIBRIUM: Tangent6 = [
    [12e9, 2.4e9, -0.6e9, 1.5e9, -0.4e9, 0.8e9],
    [2.4e9, 9e9, 1.2e9, -0.3e9, 0.7e9, -0.2e9],
    [-0.6e9, 1.2e9, 7e9, 0.4e9, -0.2e9, 0.3e9],
    [1.5e9, -0.3e9, 0.4e9, 5e9, 0.4e9, -0.2e9],
    [-0.4e9, 0.7e9, -0.2e9, 0.4e9, 4.5e9, 0.3e9],
    [0.8e9, -0.2e9, 0.3e9, -0.2e9, 0.3e9, 3.5e9],
];
const SLOPE: Tangent6 = [
    [-1e6, 0.2e6, -0.3e6, 2e6, -0.5e6, 0.4e6],
    [0.2e6, 2e6, 0.1e6, -0.4e6, 0.2e6, 0.3e6],
    [-0.3e6, 0.1e6, -3e6, 0.3e6, -0.2e6, 0.1e6],
    [2e6, -0.4e6, 0.3e6, 0.5e6, 0.2e6, -0.3e6],
    [-0.5e6, 0.2e6, -0.2e6, 0.2e6, -0.4e6, 0.1e6],
    [0.4e6, 0.3e6, 0.1e6, -0.3e6, 0.1e6, 0.7e6],
];
const ROTATION: [[f64; 3]; 3] = [[0.8, -0.6, 0.0], [0.48, 0.64, -0.6], [0.36, 0.48, 0.8]];

fn fixture(thermal: bool) -> (TensorPronyProperties, Vec<PropertyClaim>) {
    let (mut properties, mut claims) = super::fixture();
    claims.retain(|claim| !matches!(claim.key.name(), "young_modulus" | "poisson_ratio"));
    claims.push(scalar("density", Density::DIMS, 1200.0));
    let keys = core::array::from_fn(|i| format!("equilibrium_c{i:02}"));
    let mut entry = 0;
    for row in 0..6 {
        for column in row..6 {
            let mut claim = scalar(&keys[entry], Pressure::DIMS, EQUILIBRIUM[row][column]);
            if thermal {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: [273.0, 300.0, 373.0]
                        .map(|t| {
                            (
                                t,
                                SLOPE[row][column].mul_add(t - 300.0, EQUILIBRIUM[row][column]),
                            )
                        })
                        .to_vec(),
                    dims: Pressure::DIMS,
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
            claims.push(claim);
            entry += 1;
        }
    }
    properties.elasticity = TensorElasticity::Mandel21 {
        stiffness_upper: keys,
    };
    (properties, claims)
}

fn point() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 100.0)
        .unwrap()
}

fn source(properties: TensorPronyProperties, claims: Vec<PropertyClaim>) -> ResolvedTensorMaterial {
    ResolvedTensorMaterial::resolve(
        &card(claims),
        &point(),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn apply(c: &Tangent6, strain: &Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| c[i][j] * strain[j]).sum())
}
fn energy(c: &Tangent6, strain: &Voigt) -> f64 {
    let stress = apply(c, strain);
    0.5 * (0..6).map(|i| strain[i] * stress[i]).sum::<f64>()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * actual.abs().max(expected.abs()).max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_full_equilibrium_source_drives_coupled_relaxation_and_joint_frame_rotation() {
    let (properties, mut claims) = fixture(false);
    // The domain of one signed coupling restricts the entire physical law.
    claims
        .iter_mut()
        .find(|c| c.key.name() == "equilibrium_c03")
        .unwrap()
        .validity = ValidityDomain::unconstrained()
        .with("T", 273.0, 373.0)
        .with("omega", 20.0, 600.0);
    let material = source(properties, claims);
    let saved = material.clone();
    assert!(material.state().property("young_modulus").is_none());
    assert!(material.state().property("poisson_ratio").is_none());
    assert_eq!(material.law().equilibrium_mandel(), &EQUILIBRIUM);
    assert_eq!(material.law().band_rad_s(), (20.0, 600.0));
    assert!(material.law().modulus_mandel(601.0).is_err());

    let strain = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mandel = tensor::to_mandel(&strain);
    let rotated = material.rotated_law(ROTATION).unwrap();
    let rotated_strain = tensor::rotate(&strain, &ROTATION);
    let mut history = material.law().history();
    let mut rotated_history = rotated.history();
    history.advance(strain, 0.0, || false).unwrap();
    rotated_history
        .advance(rotated_strain, 0.0, || false)
        .unwrap();
    // A constrained axial strain develops actual shear stress. Dropping the
    // off-diagonal entries would turn this into a different constitutive law.
    assert!(history.stress()[3] > 1e6);
    let before_hold = history.stored_energy();
    let held = history.advance(strain, 0.025, || false).unwrap();
    let rotated_held = rotated_history
        .advance(rotated_strain, 0.025, || false)
        .unwrap();
    let decay = (-0.025_f64 / 0.01).exp();
    let expected_mandel: Voigt = core::array::from_fn(|i| {
        (0..6)
            .map(|j| (EQUILIBRIUM[i][j] + decay * BRANCH[i][j]) * mandel[j])
            .sum()
    });
    for (actual, expected) in held
        .stress
        .into_iter()
        .zip(tensor::from_mandel(&expected_mandel))
    {
        close(actual, expected);
    }
    close(
        held.dissipated,
        energy(&BRANCH, &mandel) * (1.0 - decay * decay),
    );
    close(
        held.stored_energy,
        energy(&EQUILIBRIUM, &mandel) + decay * decay * energy(&BRANCH, &mandel),
    );
    close(held.work, 0.0);
    close(before_hold - held.stored_energy, held.dissipated);
    for (actual, expected) in rotated_held
        .stress
        .into_iter()
        .zip(tensor::rotate(&held.stress, &ROTATION))
    {
        close(actual, expected);
    }
    close(rotated_held.stored_energy, held.stored_energy);
    close(rotated_held.dissipated, held.dissipated);
    assert_eq!(
        material, saved,
        "native evolution cannot mutate source data"
    );
}

#[test]
fn g1_independent_affine_couplings_reach_rotated_free_energy_and_entropy() {
    let (properties, claims) = fixture(true);
    let card = card(claims);
    let material = ResolvedTensorMaterial::resolve(
        &card,
        &point(),
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = material.clone();
    let thermal =
        thermal::TensorThermalLaw::try_new(&material, &card, "T", [280.0, 360.0]).unwrap();
    let affine =
        thermal::AffineTensorThermal::try_new(&material, &card, "T", [280.0, 360.0]).unwrap();
    let derivative = thermal.in_frame(ROTATION).unwrap();
    let affine_derivative = affine.in_frame(ROTATION).unwrap();
    assert!(derivative.is_affine());
    assert_eq!(
        derivative.equilibrium_pa_k,
        affine_derivative.equilibrium_pa_k
    );
    let column: Voigt = [1.0, -0.3, 0.2, 0.7, -0.2, 0.4];
    let displacement = 0.001;
    let strain = tensor::to_mandel(&column.map(|v| v * displacement));
    let volume = 0.02;
    let expected_entropy = -volume * energy(&SLOPE, &strain);
    assert!(
        expected_entropy.abs() > 0.01,
        "thermal normal/shear coupling must be live"
    );
    for temperature in [280.0, 300.0, 360.0] {
        let endpoint = material
            .at_query_point(&card, &point().with("T", temperature).unwrap())
            .unwrap();
        assert_eq!(endpoint.properties(), material.properties());
        let law = endpoint.rotated_law(ROTATION).unwrap();
        let projection = law
            .project_strain(
                &[tensor::rotate(&column, &ROTATION)],
                volume,
                TensorProjectionBudget {
                    max_coordinates: 1,
                    max_scalar_terms: 256,
                },
                || false,
            )
            .unwrap();
        let response = projection
            .thermal_state(
                &derivative,
                &[displacement],
                &[0.0; 6],
                temperature,
                1e-10,
                || false,
            )
            .unwrap();
        let current = core::array::from_fn(|i| {
            core::array::from_fn(|j| SLOPE[i][j].mul_add(temperature - 300.0, EQUILIBRIUM[i][j]))
        });
        let expected_free_energy = volume * energy(&current, &strain);
        close(response.free_energy_j, expected_free_energy);
        close(response.entropy_j_k, expected_entropy);
        close(
            response.internal_energy_j,
            expected_free_energy + temperature * expected_entropy,
        );
        close(response.heat_capacity_j_k, 0.0);
    }
    assert_eq!(material, original);
}

#[test]
fn g0_full_equilibrium_refuses_missing_units_dynamic_samples_and_unstable_tensors() {
    for mutation in ["missing", "units", "frequency", "indefinite"] {
        let (properties, mut claims) = fixture(false);
        let name = "equilibrium_c03";
        if mutation == "missing" {
            claims.retain(|c| c.key.name() != name);
        } else {
            let claim = claims.iter_mut().find(|c| c.key.name() == name).unwrap();
            match mutation {
                "units" => {
                    claim.key = PropertyKey::new(name, Dims::NONE);
                    claim.value = PropertyValue::Scalar {
                        value: 1.5e9,
                        dims: Dims::NONE,
                    };
                }
                "frequency" => {
                    claim.value = PropertyValue::Curve {
                        abscissa: "omega".into(),
                        abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
                        knots: vec![(0.0, 1.5e9), (1000.0, 2e9)],
                        dims: Pressure::DIMS,
                    };
                    claim.interpolation = InterpolationPolicy::LinearInside;
                }
                _ => {
                    claim.value = PropertyValue::Scalar {
                        value: 1e12,
                        dims: Pressure::DIMS,
                    }
                }
            }
        }
        assert!(
            ResolvedTensorMaterial::resolve(
                &card(claims),
                &point(),
                properties,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err(),
            "accepted {mutation}"
        );
    }
}

#[test]
fn g0_thermal_coupling_requires_complete_affine_stable_source_paths() {
    for mutation in ["nonaffine", "unstable_endpoint"] {
        let (properties, mut claims) = fixture(true);
        let claim = claims
            .iter_mut()
            .find(|c| c.key.name() == "equilibrium_c03")
            .unwrap();
        let PropertyValue::Curve { knots, .. } = &mut claim.value else {
            unreachable!()
        };
        if mutation == "nonaffine" {
            knots[2].1 += 1e7;
        } else {
            for (temperature, value) in knots {
                *value = 1e9_f64.mul_add(*temperature - 300.0, 1.5e9);
            }
        }
        let card = card(claims);
        let material = ResolvedTensorMaterial::resolve(
            &card,
            &point(),
            properties,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let saved = material.clone();
        assert!(
            thermal::TensorThermalLaw::try_new(&material, &card, "T", [280.0, 360.0]).is_err(),
            "accepted {mutation}"
        );
        assert_eq!(material, saved);
    }
}
