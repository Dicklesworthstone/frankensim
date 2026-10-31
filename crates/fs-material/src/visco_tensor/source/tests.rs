//! G0/G3: complete material-card data, source domains and original tensor laws.
mod kelvin;
use super::*;
use crate::state_point::{MaterialPropertySelection, resolve_material_state_point};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};

const BRANCH: [[f64; 6]; 6] = [
    [9e8, 2e7, -1e7, 3e7, 4e7, -2e7],
    [2e7, 8e8, 5e7, -3e7, 2e7, 1e7],
    [-1e7, 5e7, 7e8, 1e7, -4e7, 2e7],
    [3e7, -3e7, 1e7, 6e8, 2e7, -1e7],
    [4e7, 2e7, -4e7, 2e7, 5e8, 3e7],
    [-2e7, 1e7, 2e7, -1e7, 3e7, 4e8],
];

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with("T", 273.0, 373.0)
            .with("omega", 0.0, 1000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured tensor source fixture".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn fixture() -> (TensorPronyProperties, Vec<PropertyClaim>) {
    let branch = TensorBranchProperties {
        stiffness_upper: core::array::from_fn(|i| format!("branch_c{i:02}")),
        relaxation_time: "branch_tau".into(),
    };
    let mut claims = vec![
        scalar("young_modulus", Pressure::DIMS, 12e9),
        scalar("poisson_ratio", Dims::NONE, 0.25),
        scalar("strain_limit", Dims::NONE, 0.02),
        scalar("branch_tau", Dims([0, 0, 1, 0, 0, 0]), 0.01),
    ];
    let mut index = 0;
    for row in 0..6 {
        for column in row..6 {
            claims.push(scalar(
                &branch.stiffness_upper[index],
                Pressure::DIMS,
                BRANCH[row][column],
            ));
            index += 1;
        }
    }
    (
        TensorPronyProperties {
            elasticity: TensorElasticity::Isotropic,
            strain_limit: "strain_limit".into(),
            branches: vec![branch],
        },
        claims,
    )
}

fn card(claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured tensor".into(),
            phase: "solid".into(),
            process: "source-binding test".into(),
            revision: 0,
        },
        set,
        vec![],
    )
    .unwrap()
}

fn resolve(
    properties: &TensorPronyProperties,
    claims: Vec<PropertyClaim>,
    permissive: bool,
) -> core::result::Result<ResolvedMaterialStatePoint, MaterialStatePointError> {
    let requirements = properties
        .requirements()
        .unwrap()
        .into_iter()
        .map(|r| {
            if permissive {
                Requirement::try_new(r.name(), r.dims(), Admit::Finite).unwrap()
            } else {
                r
            }
        })
        .collect::<Vec<_>>();
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 100.0)
        .unwrap();
    resolve_material_state_point(
        &card(claims),
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
}

#[test]
fn full_source_mandel_matrix_preserves_shear_cross_terms_and_relaxation() {
    let (properties, claims) = fixture();
    let state = resolve(&properties, claims, false).unwrap();
    let before = state.clone();
    let law = properties.bind(&state).unwrap();
    assert_eq!(state, before, "binding must preserve original evidence");
    assert_eq!(law.branches()[0].c_pa, BRANCH);
    assert_eq!(law.branches()[0].tau_s, 0.01);
    assert_eq!(law.c_eq()[0][0], 14.4e9);
    assert_eq!(law.c_eq()[0][1], 4.8e9);
    assert_eq!(law.c_eq()[3][3], 9.6e9, "Mandel shear stiffness is 2 mu");
    let (storage, loss) = law.modulus_mandel(100.0).unwrap();
    for row in 0..6 {
        for column in 0..6 {
            // omega*tau = 1 gives equal storage/loss increments C_j/2.
            let expected = 0.5 * BRANCH[row][column];
            assert!(
                (storage[row][column] - law.c_eq()[row][column] - expected).abs()
                    < 1e-12 * law.c_eq()[row][row]
            );
            assert!((loss[row][column] - expected).abs() < 1e-12 * BRANCH[row][row]);
        }
    }
    assert!(
        loss[0][2] < 0.0,
        "valid signed coupling must not be clipped"
    );
    let mut history = law.history();
    let step = history
        .advance([0.001, -0.002, 0.001, 0.001, 0.0, -0.001], 0.02, || false)
        .unwrap();
    assert!(step.dissipated > 0.0);
    assert!(step.balance_residual.abs() < 1e-10 * step.work.abs());
    assert!(history.advance([0.03; 6], 0.01, || false).is_err());
}

fn isotropic_pair_fixture(elasticity: TensorElasticity) -> (TensorPronyProperties, Vec<PropertyClaim>) {
    let (mut properties, mut claims) = fixture();
    properties.elasticity = elasticity;
    match elasticity {
        TensorElasticity::Isotropic => {}
        TensorElasticity::IsotropicYoungShear => {
            claims.retain(|c| c.key.name() != "poisson_ratio");
            claims.push(scalar("shear_modulus", Pressure::DIMS, 4.8e9));
        }
        TensorElasticity::IsotropicBulkShear => {
            claims.retain(|c| c.key.name() != "young_modulus" && c.key.name() != "poisson_ratio");
            claims.push(scalar("bulk_modulus", Pressure::DIMS, 8e9));
            claims.push(scalar("shear_modulus", Pressure::DIMS, 4.8e9));
        }
        TensorElasticity::Orthotropic => unreachable!("isotropic fixture"),
    }
    (properties, claims)
}

#[test]
fn independent_equilibrium_pairs_preserve_complete_tensor_histories_and_source_receipts() {
    let close = |a: f64, b: f64| {
        assert!((a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
    };
    let (reference_properties, reference_claims) = fixture();
    let reference = reference_properties
        .bind(&resolve(&reference_properties, reference_claims, false).unwrap())
        .unwrap();
    for elasticity in [
        TensorElasticity::Isotropic,
        TensorElasticity::IsotropicYoungShear,
        TensorElasticity::IsotropicBulkShear,
    ] {
        let (properties, mut claims) = isotropic_pair_fixture(elasticity);
        claims.push(scalar("density", Density::DIMS, 900.0));
        let restricted = if elasticity == TensorElasticity::Isotropic {
            "poisson_ratio"
        } else {
            "shear_modulus"
        };
        claims.iter_mut().find(|c| c.key.name() == restricted).unwrap().validity =
            ValidityDomain::unconstrained().with("T", 273.0, 373.0).with("omega", 20.0, 600.0);
        let material = card(claims);
        let point = QueryPoint::new().with("T", 300.0).unwrap().with("omega", 100.0).unwrap();
        let source = ResolvedTensorMaterial::resolve(
            &material, &point, properties.clone(), MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let original = source.state().clone();
        assert_eq!(source.properties(), &properties);
        assert_eq!(source.density_kg_m3(), 900.0);
        assert_eq!(source.law().band_rad_s(), (20.0, 600.0));
        assert_eq!(source.law().branches()[0].c_pa, BRANCH);
        assert_eq!(source.law().branches()[0].tau_s, 0.01);
        if elasticity != TensorElasticity::Isotropic {
            assert!(source.state().property("poisson_ratio").is_none());
        }
        if elasticity == TensorElasticity::IsotropicBulkShear {
            assert!(source.state().property("young_modulus").is_none());
        }
        for (a, b) in source.law().c_eq().iter().flatten().zip(reference.c_eq().iter().flatten()) {
            close(*a, *b);
        }
        let mut history = source.law().history();
        let mut expected = reference.history();
        for (strain, dt) in [
            ([0.001, -0.002, 0.001, 0.001, 0.0, -0.001], 0.02),
            ([0.001, -0.002, 0.001, 0.001, 0.0, -0.001], 0.04),
            ([0.0002, 0.0001, -0.0003, 0.0, -0.0005, 0.0001], 0.03),
        ] {
            let actual = history.advance(strain, dt, || false).unwrap();
            let expected = expected.advance(strain, dt, || false).unwrap();
            for (a, b) in actual.stress.into_iter().zip(expected.stress) {
                close(a, b);
            }
            close(actual.stored_energy, expected.stored_energy);
            close(actual.work, expected.work);
            close(actual.dissipated, expected.dissipated);
            assert!(actual.dissipated > 0.0);
            assert!(actual.balance_residual.abs() < 1e-10 * actual.work.abs().max(1.0));
        }
        let before_cancel = history.clone();
        assert!(history.advance([0.0; 6], 0.01, || true).is_err());
        assert_eq!(history, before_cancel);
        let refreshed = source.at_query_point(
            &material, &point.clone().with("T", 310.0).unwrap(),
        ).unwrap();
        assert_eq!(refreshed.law(), source.law());
        assert_eq!(refreshed.properties(), source.properties());
        for (before, after) in original.properties().iter().zip(refreshed.state().properties()) {
            assert_eq!(before.requirement(), after.requirement());
            assert_eq!(before.answer().receipt.selected, after.answer().receipt.selected);
            material.claims().verify_receipt(&after.answer().receipt).unwrap();
        }
        assert_eq!(source.state(), &original);
        // Constant stored tensors retain the kinetics-only model and are
        // the zero-derivative limit of the same affine thermodynamic owner.
        source.admit_thermal_relaxation(&material, "T", [280.0, 320.0]).unwrap();
        let thermal = thermal::AffineTensorThermal::try_new(
            &source, &material, "T", [280.0, 320.0],
        ).unwrap();
        let derivative = thermal.in_frame(
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        ).unwrap();
        assert!(derivative.equilibrium_pa_k.iter().flatten().all(|v| *v == 0.0));
        assert!(derivative.branches_pa_k.iter().flatten().flatten().all(|v| *v == 0.0));
    }
}

#[test]
fn independent_tensor_equilibrium_rechecks_pair_stability_and_requires_original_shear() {
    for (elasticity, name, value) in [
        (TensorElasticity::IsotropicYoungShear, "young_modulus", 14.4e9),
        (TensorElasticity::IsotropicYoungShear, "shear_modulus", 0.0),
        (TensorElasticity::IsotropicBulkShear, "bulk_modulus", -1.0),
        (TensorElasticity::IsotropicBulkShear, "shear_modulus", 0.0),
    ] {
        let (properties, mut claims) = isotropic_pair_fixture(elasticity);
        let claim = claims.iter_mut().find(|c| c.key.name() == name).unwrap();
        claim.value = PropertyValue::Scalar { value, dims: Pressure::DIMS };
        // Even a caller's finite-only bundle cannot bypass constitutive stability.
        let state = resolve(&properties, claims, true).unwrap();
        assert!(properties.bind(&state).is_err(), "accepted {elasticity:?}: {name}={value}");
    }
    for elasticity in [TensorElasticity::IsotropicYoungShear, TensorElasticity::IsotropicBulkShear] {
        let (properties, mut claims) = isotropic_pair_fixture(elasticity);
        claims.retain(|c| c.key.name() != "shear_modulus");
        assert!(resolve(&properties, claims, false).is_err());
        let (properties, mut claims) = isotropic_pair_fixture(elasticity);
        let shear = claims.iter_mut().find(|c| c.key.name() == "shear_modulus").unwrap();
        shear.value = PropertyValue::Curve {
            abscissa: "omega".into(),
            abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
            knots: vec![(0.0, 4.8e9), (1000.0, 6e9)],
            dims: Pressure::DIMS,
        };
        shear.interpolation = InterpolationPolicy::LinearInside;
        let state = resolve(&properties, claims, false).unwrap();
        assert!(properties.bind(&state).is_err(), "a sampled shear spectrum is not equilibrium");
    }
}

#[test]
fn orthotropic_source_retains_all_three_elastic_directions() {
    let (mut properties, mut claims) = fixture();
    properties.elasticity = TensorElasticity::Orthotropic;
    claims.retain(|c| c.key.name() != "young_modulus" && c.key.name() != "poisson_ratio");
    let young = [12e9, 3e9, 1.5e9];
    let poisson = [0.2, 0.15, 0.1];
    let shear = [1.2e9, 0.8e9, 0.6e9];
    for i in 0..3 {
        claims.push(scalar(
            ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i],
            Pressure::DIMS,
            young[i],
        ));
        claims.push(scalar(
            ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i],
            Dims::NONE,
            poisson[i],
        ));
        claims.push(scalar(
            ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i],
            Pressure::DIMS,
            shear[i],
        ));
    }
    let law = properties
        .bind(&resolve(&properties, claims.clone(), false).unwrap())
        .unwrap();
    let expected = TensorProny::from_tensor_shear(
        OrthotropicElastic::new(young, poisson, shear, 0.02)
            .unwrap()
            .stiffness(),
        vec![],
        0.02,
    )
    .unwrap();
    assert_eq!(law.c_eq(), expected.c_eq());
    assert_eq!(law.branches()[0].c_pa, BRANCH);
    let third_axis = claims
        .iter_mut()
        .find(|c| c.key.name() == ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[2])
        .unwrap();
    third_axis.value = PropertyValue::Scalar {
        value: 1.8e9,
        dims: Pressure::DIMS,
    };
    let changed = properties
        .bind(&resolve(&properties, claims.clone(), false).unwrap())
        .unwrap();
    assert_ne!(law.c_eq()[2][2], changed.c_eq()[2][2]);
    claims.retain(|c| c.key.name() != ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[2]);
    assert!(resolve(&properties, claims, false).is_err());
}

#[test]
fn repeated_source_properties_are_deduplicated_and_incompatible_roles_refuse() {
    let (mut properties, _) = fixture();
    properties.branches.push(properties.branches[0].clone());
    assert_eq!(properties.requirements().unwrap().len(), 25);
    properties.branches[1].relaxation_time = "young_modulus".into();
    assert!(properties.requirements().is_err());
    properties.branches.pop();
    properties.strain_limit = "young_modulus".into();
    assert!(properties.requirements().is_err());
    properties.strain_limit = "strain_limit".into();
    properties.branches = vec![properties.branches[0].clone(); 64];
    assert_eq!(properties.requirements().unwrap().len(), 25);
    properties.branches.push(properties.branches[0].clone());
    assert!(properties.requirements().is_err());
}

#[test]
fn all_selected_source_domains_intersect_and_the_strain_claim_is_required() {
    let (properties, mut claims) = fixture();
    for claim in &mut claims {
        let band = match claim.key.name() {
            "young_modulus" => (1.0, 900.0),
            "poisson_ratio" => (10.0, 1000.0),
            "strain_limit" => (5.0, 700.0),
            "branch_tau" => (0.0, 600.0),
            _ => (20.0, 800.0),
        };
        claim.validity = ValidityDomain::unconstrained()
            .with("T", 273.0, 373.0)
            .with("omega", band.0, band.1);
    }
    let law = properties
        .bind(&resolve(&properties, claims.clone(), false).unwrap())
        .unwrap();
    assert_eq!(law.band_rad_s(), (20.0, 600.0));
    let original = card(claims.clone());
    let outside = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("omega", 100.0)
        .unwrap();
    assert!(
        resolve_material_state_point(
            &original,
            &outside,
            &properties.requirements().unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let outside = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 650.0)
        .unwrap();
    assert!(
        resolve_material_state_point(
            &original,
            &outside,
            &properties.requirements().unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    claims.retain(|c| c.key.name() != "strain_limit");
    assert!(resolve(&properties, claims, false).is_err());
}

#[test]
fn invalid_source_coefficients_dimensions_and_frequency_sampled_springs_refuse() {
    let (properties, claims) = fixture();
    // A caller may supply a bundle admitted with looser scalar requirements;
    // constitutive admission still has to reject it at this binding boundary.
    for (name, value) in [
        ("young_modulus", -1.0),
        ("poisson_ratio", 0.5),
        ("strain_limit", 0.0),
        ("branch_tau", 0.0),
        ("branch_c00", -1.0),
    ] {
        let mut altered = claims.clone();
        let claim = altered.iter_mut().find(|c| c.key.name() == name).unwrap();
        let PropertyValue::Scalar { value: entry, .. } = &mut claim.value else {
            unreachable!()
        };
        *entry = value;
        let state = resolve(&properties, altered, true).unwrap();
        assert!(properties.bind(&state).is_err(), "admitted {name}={value}");
    }
    let mut wrong_units = claims.clone();
    let claim = wrong_units
        .iter_mut()
        .find(|c| c.key.name() == "branch_c00")
        .unwrap();
    *claim = scalar("branch_c00", Dims::NONE, 9e8);
    assert!(resolve(&properties, wrong_units, false).is_err());
    for value in [f64::NAN, f64::INFINITY] {
        assert!(
            ClaimSet::new()
                .insert_claim(scalar("branch_c00", Pressure::DIMS, value))
                .is_err()
        );
    }
    let mut sampled = claims;
    let claim = sampled
        .iter_mut()
        .find(|c| c.key.name() == "young_modulus")
        .unwrap();
    claim.value = PropertyValue::Curve {
        abscissa: "omega".into(),
        abscissa_dims: Dims([0, 0, -1, 0, 0, 0]),
        knots: vec![(0.0, 12e9), (1000.0, 18e9)],
        dims: Pressure::DIMS,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    let state = resolve(&properties, sampled, false).unwrap();
    assert!(
        properties.bind(&state).is_err(),
        "one frequency sample is not an equilibrium law"
    );
}

#[test]
fn held_source_moisture_reaches_inertia_equilibrium_and_complete_tensor_memory() {
    use crate::state_point::query_at_fixed_environment;
    use fs_qty::semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm};
    let moisture = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass), ValueForm::Static));
    let (properties, mut claims) = fixture();
    claims.push(scalar("density", Density::DIMS, 850.0));
    for claim in &mut claims {
        claim.validity = claim.validity.clone()
            .with_quantity("moisture_content", moisture, 0.0, 0.2);
        let factor = match claim.key.name() {
            "density" => 1.2,
            "young_modulus" => 0.5,
            "branch_tau" => 2.0,
            name if name.starts_with("branch_c") => 0.5,
            _ => continue,
        };
        let PropertyValue::Scalar { value, dims } = claim.value.clone() else { unreachable!() };
        claim.value = PropertyValue::Curve { abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE, knots: vec![(0.0, value), (0.2, factor*value)], dims };
        claim.interpolation = InterpolationPolicy::LinearInside;
    }
    let source = card(claims);
    let at = |value| query_at_fixed_environment(&source, 300.0, 100.0,
        &[("moisture_content".into(), value)]).unwrap();
    let dry = ResolvedTensorMaterial::resolve(&source, &at(0.0), properties,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let wet = dry.at_query_point(&source, &at(0.2)).unwrap();
    assert_eq!(wet.density_kg_m3(), 1.2*dry.density_kg_m3());
    assert_eq!(wet.law().branches()[0].tau_s, 2.0*dry.law().branches()[0].tau_s);
    for i in 0..6 {
        for j in 0..6 {
            assert_eq!(wet.law().c_eq()[i][j], 0.5*dry.law().c_eq()[i][j]);
            assert_eq!(wet.law().branches()[0].c_pa[i][j], 0.5*dry.law().branches()[0].c_pa[i][j]);
        }
    }
    let (dry_storage, dry_loss) = dry.law().modulus_mandel(100.0).unwrap();
    let (wet_storage, wet_loss) = wet.law().modulus_mandel(100.0).unwrap();
    assert_ne!(dry_storage, wet_storage);
    assert_ne!(dry_loss, wet_loss);
    for law in [dry.law(), wet.law()] {
        let mut history = law.history();
        let step = history.advance([0.001, -0.002, 0.001, 0.001, 0.0, -0.001],
            0.02, || false).unwrap();
        assert!(step.dissipated > 0.0);
        assert!(step.balance_residual.abs() < 1e-10*step.work.abs());
    }
    assert_eq!(dry.state().card_identity(), wet.state().card_identity());
    assert_ne!(dry.state().identity(), wet.state().identity());
    assert!(dry.at_query_point(&source, &at(0.3)).is_err());
}

fn affine_equilibrium_pair_fixture(
    elasticity: TensorElasticity,
    coefficients: [(&str, f64, f64); 2],
) -> (MaterialCard, ResolvedTensorMaterial) {
    let (mut properties, mut claims) = fixture();
    properties.elasticity = elasticity;
    claims.retain(|claim| !matches!(claim.key.name(), "young_modulus" | "poisson_ratio"));
    claims.push(scalar("density", Density::DIMS, 850.0));
    for (name, reference, slope) in coefficients {
        let mut claim = scalar(name, Pressure::DIMS, reference);
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: [273.0, 373.0]
                .into_iter()
                .map(|t| (t, slope.mul_add(t - 300.0, reference)))
                .collect(),
            dims: Pressure::DIMS,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        claims.push(claim);
    }
    let card = card(claims);
    let query = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("omega", 100.0)
        .unwrap();
    let material = ResolvedTensorMaterial::resolve(
        &card,
        &query,
        properties,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, material)
}

#[test]
fn independent_affine_bulk_shear_reaches_native_tensor_entropy_and_internal_energy() {
    use crate::visco_tensor::TensorProjectionBudget;
    let (card, material) = affine_equilibrium_pair_fixture(
        TensorElasticity::IsotropicBulkShear,
        [("bulk_modulus", 8e9, -10e6), ("shear_modulus", 4.8e9, 6e6)],
    );
    let before = material.clone();
    let thermal =
        thermal::AffineTensorThermal::try_new(&material, &card, "T", [280.0, 360.0]).unwrap();
    let derivative = thermal
        .in_frame([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        .unwrap();
    let column: [f64; 6] = [1.0, -0.5, 0.25, 0.3, -0.1, 0.2];
    let displacement = 0.001;
    let strain = column.map(|v| v * displacement);
    let volume = 0.01;
    let trace = strain[0] + strain[1] + strain[2];
    let dev_norm = strain[..3]
        .iter()
        .map(|v| (v - trace / 3.0).powi(2))
        .sum::<f64>()
        + 2.0 * strain[3..].iter().map(|v| v * v).sum::<f64>();
    // Independent isotropic continuum potential, including tensor-shear weights.
    let expected_entropy = -volume * (0.5 * -10e6 * trace * trace + 6e6 * dev_norm);
    let close = |actual: f64, expected: f64| {
        assert!(
            (actual - expected).abs() < 2e-11 * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    };
    for temperature in [280.0, 300.0, 360.0] {
        let endpoint = material
            .at_query_point(
                &card,
                &material.state().query().with("T", temperature).unwrap(),
            )
            .unwrap();
        let projection = endpoint
            .law()
            .project_strain(
                &[column],
                volume,
                TensorProjectionBudget {
                    max_coordinates: 1,
                    max_scalar_terms: 256,
                },
                || false,
            )
            .unwrap();
        let state = projection
            .thermal_state(
                &derivative,
                &[displacement],
                &[0.0; 6],
                temperature,
                1e-10,
                || false,
            )
            .unwrap();
        let bulk = (-10e6_f64).mul_add(temperature - 300.0, 8e9);
        let shear = 6e6_f64.mul_add(temperature - 300.0, 4.8e9);
        let expected_free_energy = volume * (0.5 * bulk * trace * trace + shear * dev_norm);
        close(state.free_energy_j, expected_free_energy);
        close(state.entropy_j_k, expected_entropy);
        close(
            state.internal_energy_j,
            expected_free_energy + temperature * expected_entropy,
        );
        assert!(endpoint.state().property("young_modulus").is_none());
        assert!(endpoint.state().property("poisson_ratio").is_none());
        assert_eq!(endpoint.law().branches()[0].c_pa, BRANCH);
    }
    assert_eq!(material, before);
}

#[test]
fn affine_young_shear_accepts_common_scale_and_refuses_changing_ratio() {
    for (shear_slope, admitted) in [(-4.8e6, true), (6e6, false)] {
        let (card, material) = affine_equilibrium_pair_fixture(
            TensorElasticity::IsotropicYoungShear,
            [
                ("young_modulus", 12e9, -12e6),
                ("shear_modulus", 4.8e9, shear_slope),
            ],
        );
        let before = material.clone();
        assert_eq!(
            thermal::AffineTensorThermal::try_new(&material, &card, "T", [280.0, 360.0]).is_ok(),
            admitted
        );
        assert_eq!(material, before);
    }
}

#[test]
fn g1_rational_young_shear_reaches_native_entropy_capacity_and_source_refresh() {
    use crate::visco_tensor::TensorProjectionBudget;
    let (card, material) = affine_equilibrium_pair_fixture(
        TensorElasticity::IsotropicYoungShear,
        [
            ("young_modulus", 12e9, -12e6),
            ("shear_modulus", 4.8e9, 6e6),
        ],
    );
    let before = material.clone();
    assert!(thermal::AffineTensorThermal::try_new(&material, &card, "T", [280.0, 360.0]).is_err());
    let law = thermal::TensorThermalLaw::try_new(&material, &card, "T", [280.0, 360.0]).unwrap();
    assert_eq!(law.source(), &material);
    assert_eq!(law.temperature_axis(), "T");
    assert_eq!(law.temperature_interval_k(), [280.0, 360.0]);
    let reference = law
        .in_frame([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        .unwrap();
    let bounds = reference.equilibrium_bounds.unwrap();
    let column: [f64; 6] = [1.0, -0.5, 0.25, 0.3, -0.1, 0.2];
    let displacement = 0.001;
    let strain = column.map(|v| v * displacement);
    let volume = 0.01;
    let trace = strain[0] + strain[1] + strain[2];
    let dev_norm = strain[..3]
        .iter()
        .map(|v| (v - trace / 3.0).powi(2))
        .sum::<f64>()
        + 2.0 * strain[3..].iter().map(|v| v * v).sum::<f64>();
    let close = |a: f64, b: f64| {
        assert!(
            (a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1.0),
            "{a} != {b}"
        );
    };
    // Quotient differentiation, independently of the implementation's
    // normalized constant-cross-product evaluation.
    let coefficients = |t: f64| {
        let e = 12e9 - 12e6 * (t - 300.0);
        let g = 4.8e9 + 6e6 * (t - 300.0);
        let numerator = e * g;
        let numerator_prime = -12e6 * g + e * 6e6;
        let numerator_second = 2.0 * -12e6 * 6e6;
        let denominator = 3.0 * (3.0 * g - e);
        let denominator_prime: f64 = 3.0 * (3.0 * 6e6 + 12e6);
        let bulk = numerator / denominator;
        let first =
            numerator_prime / denominator - numerator * denominator_prime / denominator.powi(2);
        let second = numerator_second / denominator
            - 2.0 * numerator_prime * denominator_prime / denominator.powi(2)
            + 2.0 * numerator * denominator_prime.powi(2) / denominator.powi(3);
        (bulk, g, first, second)
    };
    let evaluate = |temperature: f64| {
        let endpoint = material
            .at_query_point(
                &card,
                &material.state().query().with("T", temperature).unwrap(),
            )
            .unwrap();
        let point = endpoint
            .law()
            .project_strain(
                &[column],
                volume,
                TensorProjectionBudget {
                    max_coordinates: 1,
                    max_scalar_terms: 256,
                },
                || false,
            )
            .unwrap();
        let derivative = reference.at_temperature(temperature).unwrap();
        assert!(std::sync::Arc::ptr_eq(
            &reference.branches_pa_k,
            &derivative.branches_pa_k
        ));
        assert_eq!(derivative.resolved_temperature_k, Some(temperature));
        if temperature != 300.0 {
            assert!(
                point
                    .thermal_state(
                        &reference,
                        &[displacement],
                        &[0.0; 6],
                        temperature,
                        1e-10,
                        || false
                    )
                    .is_err(),
                "stale nonlinear derivative accepted"
            );
        }
        let value = point
            .thermal_state(
                &derivative,
                &[displacement],
                &[0.0; 6],
                temperature,
                1e-10,
                || false,
            )
            .unwrap();
        let (bulk, shear, bulk_prime, bulk_second) = coefficients(temperature);
        let expected_free = volume * (0.5 * bulk * trace * trace + shear * dev_norm);
        let expected_entropy = -volume * (0.5 * bulk_prime * trace * trace + 6e6 * dev_norm);
        let expected_capacity = -temperature * volume * 0.5 * bulk_second * trace * trace;
        close(value.free_energy_j, expected_free);
        close(value.entropy_j_k, expected_entropy);
        close(
            value.internal_energy_j,
            expected_free + temperature * expected_entropy,
        );
        close(value.heat_capacity_j_k, expected_capacity);
        close(derivative.equilibrium_pa_k[0][0], bulk_prime + 8e6);
        close(derivative.equilibrium_pa_k2[0][1], bulk_second);
        assert!(value.heat_capacity_j_k < 0.0);
        assert!(3.0 * bulk <= bounds.stiffness_norm_pa);
        assert!(2.0 * shear <= bounds.stiffness_norm_pa);
        assert!(3.0 * bulk_prime.abs() <= bounds.slope_norm_pa_k);
        assert!(3.0 * bulk_second <= bounds.curvature_norm_pa_k2);
        assert!(
            point
                .thermal_capacity_bound_j_k2(&derivative, temperature, [0.0; 6])
                .unwrap()
                >= -value.heat_capacity_j_k / temperature
        );
        assert!(endpoint.state().property("poisson_ratio").is_none());
        assert_eq!(endpoint.law().branches()[0].c_pa, BRANCH);
        for property in material.state().properties() {
            assert_eq!(
                property.answer().receipt.selected,
                endpoint
                    .state()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        value
    };
    for temperature in [280.0, 295.0, 300.0, 315.0, 360.0] {
        evaluate(temperature);
    }
    let h = 0.01;
    let difference =
        (evaluate(320.0 + h).internal_energy_j - evaluate(320.0 - h).internal_energy_j) / (2.0 * h);
    let capacity = evaluate(320.0).heat_capacity_j_k;
    assert!((difference - capacity).abs() < 2e-7 * capacity.abs());
    // The same nonlinear source with a constant free-strain derivative and
    // fixed, nonzero viscous strain exercises curvature and expansion jointly.
    let alpha = [2e-6, 2e-6, 2e-6, 0.0, 0.0, 0.0];
    let total_mandel: [f64; 6] = core::array::from_fn(|i| {
        strain[i]
            * if i < 3 {
                1.0
            } else {
                core::f64::consts::SQRT_2
            }
    });
    let branch_bilinear = |a: &[f64; 6], b: &[f64; 6]| {
        let mut value = 0.0;
        for i in 0..6 {
            for j in 0..6 {
                value += a[i] * BRANCH[i][j] * b[j];
            }
        }
        value
    };
    let evaluate_expansion = |temperature: f64| {
        let endpoint = material
            .at_query_point(
                &card,
                &material.state().query().with("T", temperature).unwrap(),
            )
            .unwrap();
        let free = alpha.map(|a| a * (temperature - 300.0));
        let point = endpoint
            .law()
            .project_strain(
                &[column],
                volume,
                TensorProjectionBudget {
                    max_coordinates: 1,
                    max_scalar_terms: 256,
                },
                || false,
            )
            .unwrap()
            .with_eigenstrain_mandel(free)
            .unwrap();
        let free_memory = point.eigenstrain_factor_coordinates().unwrap();
        let memory: Vec<_> = point.branches()[0]
            .rows()
            .iter()
            .zip(free_memory)
            .map(|(row, offset)| 0.6 * displacement * row[0] - offset)
            .collect();
        let spring: [f64; 6] = core::array::from_fn(|i| 0.6 * total_mandel[i] - free[i]);
        let derivative = reference.at_temperature(temperature).unwrap();
        let value = point
            .thermal_state_with_expansion(
                &derivative,
                &[displacement],
                &memory,
                temperature,
                alpha,
                1e-10,
                || false,
            )
            .unwrap();
        let (bulk, shear, bulk_prime, bulk_second) = coefficients(temperature);
        let elastic_trace = trace - 6e-6 * (temperature - 300.0);
        let free_energy = volume
            * (0.5 * bulk * elastic_trace.powi(2)
                + shear * dev_norm
                + 0.5 * branch_bilinear(&spring, &spring));
        let entropy = volume
            * (6e-6 * bulk * elastic_trace
                - 0.5 * bulk_prime * elastic_trace.powi(2)
                - 6e6 * dev_norm
                + branch_bilinear(&alpha, &spring));
        let curvature = volume
            * (0.5 * bulk_second * elastic_trace.powi(2) - 12e-6 * bulk_prime * elastic_trace
                + 36e-12 * bulk
                + branch_bilinear(&alpha, &alpha));
        close(value.free_energy_j, free_energy);
        close(value.entropy_j_k, entropy);
        close(value.internal_energy_j, free_energy + temperature * entropy);
        close(value.heat_capacity_j_k, -temperature * curvature);
        assert!(
            point
                .thermal_capacity_bound_j_k2(&derivative, temperature, alpha)
                .unwrap()
                >= curvature
        );
        value
    };
    for temperature in [280.0, 295.0, 300.0, 315.0, 360.0] {
        evaluate_expansion(temperature);
    }
    let difference = (evaluate_expansion(320.0 + h).internal_energy_j
        - evaluate_expansion(320.0 - h).internal_energy_j)
        / (2.0 * h);
    let capacity = evaluate_expansion(320.0).heat_capacity_j_k;
    assert!((difference - capacity).abs() < 2e-7 * capacity.abs());
    assert!(reference.at_temperature(279.0).is_err());
    assert!(reference.at_temperature(f64::NAN).is_err());
    assert_eq!(material, before);
}

#[test]
fn g3_rational_tensor_law_retains_complete_source_path_admission() {
    for mutation in ["kink", "narrow", "unstable"] {
        let (properties, mut claims) =
            isotropic_pair_fixture(TensorElasticity::IsotropicYoungShear);
        claims.push(scalar("density", Density::DIMS, 850.0));
        for claim in &mut claims {
            let (value, slope): (f64, f64) = match claim.key.name() {
                "young_modulus" => (12e9, -12e6),
                "shear_modulus" => (4.8e9, if mutation == "unstable" { -25e6 } else { 6e6 }),
                _ => continue,
            };
            let mut knots: Vec<_> = [273.0, 300.0, 373.0]
                .into_iter()
                .map(|t| (t, slope.mul_add(t - 300.0, value)))
                .collect();
            if claim.key.name() == "shear_modulus" {
                if mutation == "kink" {
                    knots[1].1 += 1e8;
                }
                if mutation == "narrow" {
                    claim.validity = ValidityDomain::unconstrained()
                        .with("T", 290.0, 330.0)
                        .with("omega", 0.0, 1000.0);
                }
            }
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots,
                dims: Pressure::DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        let card = card(claims);
        let point = QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("omega", 100.0)
            .unwrap();
        let material = ResolvedTensorMaterial::resolve(
            &card,
            &point,
            properties,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let before = material.clone();
        assert!(
            thermal::TensorThermalLaw::try_new(&material, &card, "T", [280.0, 360.0]).is_err(),
            "accepted {mutation}"
        );
        assert_eq!(material, before);
    }
}
