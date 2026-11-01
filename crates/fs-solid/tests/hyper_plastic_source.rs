//! Source-bound heterogeneous finite plasticity through actual tetrahedral assembly.
//!
//! The fully prescribed affine patch isolates source selection, reference-volume
//! integration and retained histories. Free-DOF equilibrium is covered by the
//! companion native finite-plastic solid tests; this is not a locking benchmark.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::plastic::finite::source::{FINITE_J2_LAW_ID, ResolvedFiniteJ2};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetSolid,
};
use fs_solid::{HyperTetSettings, TetAssemblyBudget};

const MU: f64 = 80.0e6;
const KAPPA: f64 = 160.0e6;
const H: f64 = 8.0e6;
const YIELDS: [f64; 2] = [2.0e6, 70.0e6];
const VOLUMES: [f64; 2] = [1.0 / 6.0, 1.0 / 3.0];
const NODES: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
];
const TETS: [[usize; 4]; 2] = [[0, 1, 2, 3], [1, 2, 3, 4]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn source_card(index: usize, elastic_limit: f64) -> ConstitutiveModelCard {
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let source = fs_blake3::hash_bytes(
        format!("synthetic heterogeneous finite-plastic solid source {index}").as_bytes(),
    );
    ConstitutiveModelCard {
        law: LawId(FINITE_J2_LAW_ID.into()),
        law_version: 1,
        parameters: [
            ("shear_modulus_pa", MU, pa),
            ("bulk_modulus_pa", KAPPA, pa),
            ("yield_stress_pa", YIELDS[index], pa),
            ("hardening_modulus_pa", H, pa),
            ("minimum_jacobian", 0.25, Dims::NONE),
            ("maximum_elastic_log_strain", elastic_limit, Dims::NONE),
            ("maximum_plastic_log_strain", 2.0, Dims::NONE),
            ("maximum_accumulated_plastic_strain", 3.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity(
            "temperature",
            absolute,
            290.0 + index as f64,
            310.0 + index as f64,
        ),
        sources: vec![source],
        provenance: Provenance {
            source: format!("synthetic finite solid constituent {index}; not measured calibration"),
            license: "test fixture".into(),
            artifact: Some(source),
        },
    }
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap()
}

fn sources(cards: &[ConstitutiveModelCard; 2]) -> [ResolvedFiniteJ2; 2] {
    cards
        .each_ref()
        .map(|card| ResolvedFiniteJ2::try_from_model_card_at(card, &context()).unwrap())
}

fn diagonal_f(volume_log: f64, a: f64) -> [f64; 3] {
    [
        (volume_log + a).exp(),
        (volume_log - 0.5 * a).exp(),
        (volume_log - 0.5 * a).exp(),
    ]
}

fn prescribed(volume_log: f64, a: f64) -> Vec<(usize, f64)> {
    let f = diagonal_f(volume_log, a);
    NODES
        .iter()
        .enumerate()
        .flat_map(|(node, x)| (0..3).map(move |i| (3 * node + i, (f[i] - 1.0) * x[i])))
        .collect()
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1.0e-5,
        // This deliberately coarse finite increment has a large signed
        // endpoint quadrature defect, independently checked below, never heat.
        maximum_endpoint_work_defect_j: 2.0e6,
        yield_tolerance_pa: 1.0e-3,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2.0e-9 * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn assert_sources(
    solid: &FinitePlasticTetSolid<'_, ResolvedFiniteJ2>,
    original: &[ResolvedFiniteJ2; 2],
) {
    let FinitePlasticLawField::PerElement(retained) = solid.material() else {
        panic!("the heterogeneous source field was replaced");
    };
    assert!(std::ptr::eq(retained.as_ptr(), original.as_ptr()));
    for (retained, original) in retained.iter().zip(original) {
        assert_eq!(retained.model_card(), original.model_card());
        assert_eq!(retained.context(), original.context());
    }
    assert_ne!(
        retained[0].model_card().content_hash(),
        retained[1].model_card().content_hash()
    );
}

#[test]
fn connected_source_elements_yield_independently_and_keep_physical_energy_separate() {
    let cards = [source_card(0, 0.5), source_card(1, 0.5)];
    let materials = sources(&cards);
    let (v, a) = (0.03, 0.20);
    let f = diagonal_f(v, a);
    let alpha = YIELDS.map(|yield_stress| (3.0 * MU * a - yield_stress).max(0.0) / (3.0 * MU + H));
    let initial = with_cx(|cx| {
        FinitePlasticTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::PerElement(&materials),
            TetAssemblyBudget::standard(),
            cx,
        )
    })
    .unwrap();
    let step = with_cx(|cx| {
        initial.advance(
            &prescribed(v, a),
            &[],
            [0.0; 3],
            HyperTetSettings {
                load_steps: 1,
                ..HyperTetSettings::default()
            },
            controls(),
            cx,
        )
    })
    .unwrap();

    let mut stored = 0.0;
    let mut dissipation = 0.0;
    let mut internal_work = 0.0;
    for index in 0..2 {
        let returned = a - alpha[index];
        let tau = [
            3.0 * KAPPA * v + 2.0 * MU * returned,
            3.0 * KAPPA * v - MU * returned,
            3.0 * KAPPA * v - MU * returned,
        ];
        close(step.solid.element_states()[index].alpha, alpha[index]);
        for component in 0..3 {
            close(
                step.solid.element_deformation()[index][4 * component],
                f[component],
            );
            close(
                step.solid.element_first_piola_pa()[index][4 * component],
                tau[component] / f[component],
            );
            internal_work += VOLUMES[index] * tau[component] / f[component] * (f[component] - 1.0);
        }
        stored += VOLUMES[index]
            * (1.5 * MU * returned * returned
                + 0.5 * KAPPA * (3.0 * v).powi(2)
                + 0.5 * H * alpha[index] * alpha[index]);
        dissipation += VOLUMES[index] * YIELDS[index] * alpha[index];
    }
    assert!(alpha[0] > 0.18);
    assert_eq!(alpha[1], 0.0);
    assert_eq!(
        step.solid.element_states()[1],
        materials[1].law().initial_state()
    );
    close(step.energy.stored_before_j, 0.0);
    close(step.energy.stored_after_j, stored);
    close(step.solid.stored_energy_j(), stored);
    close(step.energy.plastic_dissipation_j, dissipation);
    close(step.energy.incremental_potential_j, stored + dissipation);
    close(step.equilibrium.strain_energy_j, stored + dissipation);
    close(step.energy.internal_endpoint_work_j, internal_work);
    close(step.energy.support_endpoint_work_j, internal_work);
    close(step.energy.dead_load_endpoint_work_j, 0.0);
    close(
        step.energy.endpoint_work_defect_j,
        internal_work - stored - dissipation,
    );
    assert!(step.energy.endpoint_work_defect_j.abs() > 1.0e6);
    assert!(step.energy.plastic_dissipation_j < 1.0e5);
    assert_sources(&initial, &materials);
    assert_sources(&step.solid, &materials);
    assert!(
        initial
            .element_states()
            .iter()
            .all(|state| state.alpha == 0.0)
    );

    // A modest unloading increment leaves BOTH independently selected histories
    // unchanged, while the previously yielded element retains its plastic metric.
    let unloaded_a = 0.18;
    let unloaded = with_cx(|cx| {
        step.solid.advance(
            &prescribed(v, unloaded_a),
            &[],
            [0.0; 3],
            HyperTetSettings {
                load_steps: 1,
                ..HyperTetSettings::default()
            },
            controls(),
            cx,
        )
    })
    .unwrap();
    assert_eq!(unloaded.solid.element_states(), step.solid.element_states());
    assert_eq!(unloaded.energy.plastic_dissipation_j, 0.0);
    let exact_unloaded_storage: f64 = (0..2)
        .map(|index| {
            let returned = unloaded_a - alpha[index];
            VOLUMES[index]
                * (1.5 * MU * returned * returned
                    + 0.5 * KAPPA * (3.0 * v).powi(2)
                    + 0.5 * H * alpha[index] * alpha[index])
        })
        .sum();
    close(unloaded.energy.stored_before_j, stored);
    close(unloaded.energy.stored_after_j, exact_unloaded_storage);
    close(
        unloaded.energy.incremental_potential_j,
        exact_unloaded_storage,
    );
    assert_sources(&unloaded.solid, &materials);
}

#[test]
fn one_constituent_source_domain_refuses_the_whole_candidate_without_history_loss() {
    let cards = [source_card(0, 0.5), source_card(1, 0.1)];
    let materials = sources(&cards);
    let initial = with_cx(|cx| {
        FinitePlasticTetSolid::try_new(
            &NODES,
            &TETS,
            FinitePlasticLawField::PerElement(&materials),
            TetAssemblyBudget::standard(),
            cx,
        )
    })
    .unwrap();
    let old_states = initial.element_states().to_vec();
    let old_displacement = initial.displacement_m().to_vec();
    let rejected = with_cx(|cx| {
        initial.advance(
            &prescribed(0.03, 0.20),
            &[],
            [0.0; 3],
            HyperTetSettings {
                load_steps: 1,
                ..HyperTetSettings::default()
            },
            controls(),
            cx,
        )
    });
    assert!(rejected.is_err());
    assert_eq!(initial.element_states(), old_states);
    assert_eq!(initial.displacement_m(), old_displacement);
    assert_eq!(initial.stored_energy_j(), 0.0);
    assert_sources(&initial, &materials);
}
