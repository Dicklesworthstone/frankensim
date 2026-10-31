//! Independent source pairs through the shared thermal/solid lowering.
//! Synthetic source curves; numerical references do not validate real materials.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    BULK_MODULUS_PROPERTY, DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, SHEAR_MODULUS_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{
    Density, Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::linear3::source::{SolidStiffnessModel, ThermalExpansionModel};
use fs_solid::linear3::{
    TetElasticMaterial, TetLinearElasticProblem, TetMaterialField, TetThermalStrainField,
};

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const FIXED: [usize; 11] = [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11];
const MODELS: [SolidStiffnessModel; 3] = [
    SolidStiffnessModel::Isotropic,
    SolidStiffnessModel::IsotropicYoungShear,
    SolidStiffnessModel::IsotropicBulkShear,
];
const FRAME: ContentHash = ContentHash([117; 32]);

fn temperature_kind() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_kind(), t)
        .unwrap()
}
fn card(omit: &[&str]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, start, end, upper) in [
        (DENSITY_PROPERTY, Density::DIMS, 6.0, 6.0, 300.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 60.0, 30.0, 400.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.25, 0.25, 400.0),
        (BULK_MODULUS_PROPERTY, Pressure::DIMS, 40.0, 20.0, 400.0),
        (SHEAR_MODULUS_PROPERTY, Pressure::DIMS, 24.0, 12.0, 400.0),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            1e-5,
            1e-5,
            400.0,
        ),
    ] {
        if omit.contains(&name) {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: if start == end {
                    PropertyValue::Scalar { value: start, dims }
                } else {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, start), (400.0, end)],
                        dims,
                    }
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_kind(),
                    300.0,
                    upper,
                ),
                interpolation: if start == end {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic independent-modulus heating fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "independent-modulus fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 117,
                kernel_id: 117,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn problem(material: &TetElasticMaterial) -> TetLinearElasticProblem<'_> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: &FIXED,
        budget: Default::default(),
    }
}

#[test]
fn g1_explicit_source_pairs_share_native_stress_and_stiffness() {
    let source = card(&[]);
    let reference = TetElasticMaterial::try_new(6.0, 60.0, 0.25, [117; 32]).unwrap();
    let strain = [0.001, -0.002, 0.003, -0.004, 0.005, 0.006];
    let mut identities = Vec::new();
    with_cx(|cx| {
        let expected = problem(&reference).assemble(cx).unwrap();
        for model in MODELS {
            let state = model
                .resolve_reference(
                    &source,
                    &point(300.0),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
            assert_eq!(state.properties().len(), 3, "no fabricated modulus claims");
            let material = model.material(&state, 6.0, FRAME, 0.01).unwrap();
            let tensor = material.stiffness_mandel_pa();
            for i in 0..6 {
                let stress: f64 = tensor[i].iter().zip(strain).map(|(c, e)| c * e).sum();
                let expected_stress = 48.0 * strain[i]
                    + if i < 3 {
                        24.0 * (strain[0] + strain[1] + strain[2])
                    } else {
                        0.0
                    };
                close(stress, expected_stress);
            }
            let actual = problem(&material).assemble(cx).unwrap();
            for (a, b) in actual
                .stiffness
                .to_dense()
                .iter()
                .zip(expected.stiffness.to_dense())
            {
                close(*a, b);
            }
            assert_eq!(actual.mass.to_dense(), expected.mass.to_dense());
            close(actual.total_mass_kg, 1.0);
            if model == SolidStiffnessModel::Isotropic {
                assert_eq!(material.material_state_identity, FRAME.0);
            }
            identities.push(material.material_state_identity);
        }
    });
    assert_ne!(identities[0], identities[1]);
    assert_ne!(identities[1], identities[2]);
}

#[test]
fn g1_g3_heating_rebuilds_stress_and_compliance_without_requerying_reference_mass() {
    // Only K and G are sourced. The hot state has neither density support nor
    // any Young/Poisson measurements that an alternate path could fall back to.
    let source = card(&[YOUNG_MODULUS_PROPERTY, POISSON_RATIO_PROPERTY]);
    let model = SolidStiffnessModel::IsotropicBulkShear;
    let reference = model
        .resolve_reference(
            &source,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    let density = reference.property(DENSITY_PROPERTY).unwrap().value_si();
    assert!(
        model
            .resolve_reference(
                &source,
                &point(400.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
    );
    let mut masses = Vec::new();
    let mut identities = Vec::new();
    with_cx(|cx| {
        for (temperature, stiffness, bulk) in [(300.0, 12.0, 40.0), (400.0, 6.0, 20.0)] {
            let state = model
                .resolve(
                    &source,
                    &point(temperature),
                    MaterialPropertySelection::SingleClaimOnly,
                )
                .unwrap();
            assert_eq!(state.properties().len(), 2);
            assert!(state.property(DENSITY_PROPERTY).is_none());
            let material = model
                .material(&state, density, reference.identity(), 0.01)
                .unwrap();
            let expansion = ThermalExpansionModel::Isotropic
                .resolve(
                    &source,
                    &point(300.0),
                    &point(temperature),
                    MaterialPropertySelection::SingleClaimOnly,
                    &material,
                    FRAME,
                    0.01,
                )
                .unwrap();
            let field = TetThermalStrainField::Uniform(expansion.thermal_strain());
            let problem = problem(&material);
            let matrix = problem.assemble(cx).unwrap();
            close(matrix.stiffness.to_dense()[0], stiffness);
            close(matrix.mass.to_dense()[0], 0.1);
            close(matrix.total_mass_kg, 1.0);
            let forces = [(3, 0.006)];
            let solved = problem
                .solve_thermal_displacement_with_nodal_forces(
                    field,
                    &forces,
                    Default::default(),
                    cx,
                )
                .unwrap();
            let thermal_strain = (temperature - 300.0) * 1e-5;
            let expected_displacement = (0.006 + bulk * thermal_strain / 2.0) / stiffness;
            close(solved.displacement_m()[1][0], expected_displacement);
            let stress = problem
                .recover_thermal_stress_with_nodal_forces(&solved, field, &forces, FRAME, 0.01, cx)
                .unwrap();
            close(stress.elements()[0].stress_mandel_pa()[0], 0.036);
            masses.push(matrix.mass.to_dense());
            identities.push(material.material_state_identity);
        }
    });
    assert_eq!(masses[0], masses[1]);
    assert_ne!(
        identities[0], identities[1],
        "hot stiffness must invalidate downstream operators"
    );
}

#[test]
fn g0_source_pair_is_explicit_and_missing_or_invalid_data_never_falls_back() {
    let source = card(&[SHEAR_MODULUS_PROPERTY]);
    // E/nu is present and usable, but cannot silently replace a selected E/G or K/G law.
    assert!(
        SolidStiffnessModel::Isotropic
            .resolve(
                &source,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_ok()
    );
    for model in [
        SolidStiffnessModel::IsotropicYoungShear,
        SolidStiffnessModel::IsotropicBulkShear,
    ] {
        assert!(
            model
                .resolve(
                    &source,
                    &point(300.0),
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
        );
    }
    let source = card(&[]);
    for model in MODELS {
        let accepted = model
            .resolve(
                &source,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
        let before = accepted.clone();
        assert!(
            model
                .resolve(
                    &source,
                    &point(401.0),
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
        );
        assert!(
            model
                .material(&accepted, 6.0, ContentHash([0; 32]), 0.01)
                .is_err()
        );
        assert!(model.material(&accepted, f64::NAN, FRAME, 0.01).is_err());
        assert_eq!(accepted, before);
    }
}
