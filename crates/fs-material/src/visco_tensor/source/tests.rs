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
