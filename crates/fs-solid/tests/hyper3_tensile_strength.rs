//! G1/G3/G4: finite-strain tensile limits on the original source-owned mechanics.
//! Synthetic Neo-Hookean coefficients verify implementation, not failure data.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::cards::ResolvedHyperelastic;
use fs_material::state_point::MaterialPropertySelection;
use fs_material::strength::{IsotropicTensileLimit, ResolvedTensileHyperelastic};
use fs_qty::{Dims, Pressure};
use fs_solid::hyper3::prepared::tensile::HyperTetTensileError;
use fs_solid::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField, TetAssemblyBudget,
};

const LIMIT: &str = "specified_cauchy_tensile_limit";
const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [-1.0, 0.0, 0.0],
];
const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [0, 4, 3, 2]];

fn with_cx(run: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 191,
                    kernel_id: 193,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn material(mu: f64, strength: f64, tag: u8) -> ResolvedTensileHyperelastic {
    let provenance = Provenance {
        source: "synthetic true-stress hyperelastic fixture".into(),
        license: "test-only".into(),
        artifact: Some(ContentHash([tag; 32])),
    };
    let validity = ValidityDomain::unconstrained().with("T", 280.0, 500.0);
    let model = ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            ("mu", mu, Pressure::DIMS),
            ("lambda", 0.0, Pressure::DIMS),
            ("stretch_limit", 4.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: validity.clone(),
        sources: vec![ContentHash([tag; 32])],
        provenance: provenance.clone(),
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(LIMIT, Pressure::DIMS),
            value: PropertyValue::Scalar {
                value: strength,
                dims: Pressure::DIMS,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance,
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic hyperelastic solid".into(),
            phase: "solid".into(),
            process: format!("manufactured specimen {tag}"),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let elastic = ResolvedHyperelastic::try_from_model_card_at(&card.models()[0], &point).unwrap();
    let limit = IsotropicTensileLimit::resolve(
        &card,
        &point,
        LIMIT,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedTensileHyperelastic::from_cauchy_limit(&card, elastic, limit).unwrap()
}

fn problem<'a>(
    tetrahedra: &'a [[usize; 4]],
    material: HyperelasticLawField<'a, ResolvedTensileHyperelastic>,
) -> HyperTetProblem<'a, ResolvedTensileHyperelastic> {
    HyperTetProblem {
        nodes_m: &NODES,
        tetrahedra,
        material,
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}

fn affine(f: [f64; 9]) -> Vec<f64> {
    let mut u = Vec::new();
    for x in NODES {
        for i in 0..3 {
            u.push(
                (0..3)
                    .map(|j| (f[3 * i + j] - if i == j { 1.0 } else { 0.0 }) * x[j])
                    .sum(),
            );
        }
    }
    u
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_hydrostatic_extension_uses_current_area_and_original_source_binding() {
    with_cx(|_, cx| {
        let material = material(100.0, 50.0, 1);
        let p = problem(&TETS[..1], HyperelasticLawField::Uniform(&material));
        let prepared = p.prepare(cx).unwrap();
        let u = affine([2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0]);
        let result = prepared.assess_tensile_strength(&u, 1, cx).unwrap();
        let sample = result.governing_sample().unwrap();
        // Independent NH formula at lambda=0: P=mu(s-1/s) I and
        // sigma=P/s^2, so using the Piola value would overstate demand fourfold.
        near(sample.response().first_piola_pa[0], 150.0);
        near(sample.response().cauchy_pa[0], 37.5);
        near(sample.response().jacobian, 8.0);
        near(sample.assessment().tensile_stress_pa(), 37.5);
        near(sample.assessment().utilization(), 0.75);
        assert!(!result.reaches_limit());
        assert_eq!(sample.material_binding_identity(), material.identity());
        assert_eq!(
            sample.limit_source_identity(),
            material.limit().resolved().identity()
        );
        assert_eq!(
            sample.response(),
            &prepared.element_responses(&u, cx).unwrap()[0]
        );
        assert_eq!(result, p.assess_tensile_strength(&u, 1, cx).unwrap());

        let compressed = affine([0.8, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.8]);
        let result = prepared
            .assess_tensile_strength(&compressed, 1, cx)
            .unwrap();
        assert_eq!(
            result
                .governing_sample()
                .unwrap()
                .assessment()
                .utilization(),
            0.0
        );
    });
}

#[test]
fn g3_native_shear_rotation_and_heterogeneous_limits_keep_local_tensile_demand() {
    with_cx(|_, cx| {
        let materials = [material(100.0, 50.0, 1), material(200.0, 50.0, 2)];
        let p = problem(&TETS, HyperelasticLawField::PerElement(&materials));
        let prepared = p.prepare(cx).unwrap();
        let gamma = 0.5_f64;
        let f = [1.0, gamma, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        // Superposed 90-degree rotation about z: R F.
        let rf = [0.0, -1.0, 0.0, 1.0, gamma, 0.0, 0.0, 0.0, 1.0];
        let original = prepared.assess_tensile_strength(&affine(f), 2, cx).unwrap();
        let rotated = prepared
            .assess_tensile_strength(&affine(rf), 2, cx)
            .unwrap();
        for (i, mu) in [100.0, 200.0].into_iter().enumerate() {
            // sigma = mu(B-I), positive root of [[gamma^2,gamma],[gamma,0]].
            let demand = 0.5 * mu * (gamma * gamma + gamma * (gamma * gamma + 4.0).sqrt());
            near(
                original.samples()[i].assessment().tensile_stress_pa(),
                demand,
            );
            near(
                rotated.samples()[i].assessment().tensile_stress_pa(),
                demand,
            );
            near(
                original.samples()[i].assessment().utilization(),
                demand / 50.0,
            );
            assert_eq!(
                original.samples()[i].material_binding_identity(),
                materials[i].identity()
            );
        }
        assert_eq!(original.governing_sample().unwrap().response().element, 1);
        assert!(original.reaches_limit());

        let homogeneous = problem(&TETS, HyperelasticLawField::Uniform(&materials[0]));
        let tied = homogeneous
            .assess_tensile_strength(&affine(f), 2, cx)
            .unwrap();
        assert_eq!(tied.governing_sample().unwrap().response().element, 0);
    });
}

#[test]
fn g4_budget_inversion_malformed_input_and_cancellation_leave_prepared_owner_reusable() {
    with_cx(|gate, cx| {
        let material = material(100.0, 50.0, 1);
        let p = problem(&TETS[..1], HyperelasticLawField::Uniform(&material));
        let prepared = p.prepare(cx).unwrap();
        let u = affine([1.1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        let original = prepared.assess_tensile_strength(&u, 1, cx).unwrap();
        assert_eq!(
            prepared.assess_tensile_strength(&u, 0, cx),
            Err(HyperTetTensileError::SampleBudget {
                required: 1,
                maximum: 0
            })
        );
        assert!(matches!(
            prepared.assess_tensile_strength(&u[..12], 1, cx),
            Err(HyperTetTensileError::Native(
                HyperTetError::InvalidInput { .. }
            ))
        ));
        let inverted = affine([-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]);
        assert!(matches!(
            prepared.assess_tensile_strength(&inverted, 1, cx),
            Err(HyperTetTensileError::Native(
                HyperTetError::MaterialRefused { .. }
            ))
        ));
        assert_eq!(
            prepared.assess_tensile_strength(&u, 1, cx).unwrap(),
            original
        );
        gate.request();
        assert_eq!(
            prepared.assess_tensile_strength(&u, 1, cx),
            Err(HyperTetTensileError::Native(HyperTetError::Cancelled))
        );
        assert_eq!(original.samples().len(), 1);
    });
}
