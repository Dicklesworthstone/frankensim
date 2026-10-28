//! G0/G1: prepared geometry retains current-state physics and load admission.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField, TetAssemblyBudget,
};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
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
            ),
        );
    });
}

fn problem(material: &Hyperelastic) -> HyperTetProblem<'_> {
    HyperTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}

fn material() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn prepared_evaluations_follow_each_current_stretch_and_exact_native_assembly() {
    with_cx(|_, cx| {
        let law = material();
        let p = problem(&law);
        let prepared = p.prepare(cx).unwrap();
        // Reusing a reference preparation must not freeze the stress or tangent.
        for stretch in [1.0_f64, 1.4, 0.75, 1.1] {
            let mut u = [0.0; 12];
            u[3] = stretch - 1.0;
            let log_j = stretch.ln();
            let expected_energy =
                (stretch * stretch - 1.0 - 2.0 * log_j + 1.5 * log_j * log_j) / 6.0;
            let expected_force = (2.0 * (stretch - 1.0 / stretch) + 3.0 * log_j / stretch) / 6.0;
            let expected_tangent = (2.0 * (1.0 + 1.0 / (stretch * stretch))
                + 3.0 * (1.0 - log_j) / (stretch * stretch))
                / 6.0;
            let force = prepared.internal_force(&u, cx).unwrap();
            let tangent = prepared.tangent_matrix(&u, cx).unwrap();
            near(prepared.stored_energy_j(&u, cx).unwrap(), expected_energy);
            near(force[3], expected_force);
            near(tangent.get(3, 3), expected_tangent);
            for component in 0..3 {
                near((0..4).map(|node| force[3 * node + component]).sum(), 0.0);
            }
            assert_eq!(force, p.internal_force(&u, cx).unwrap());
            assert_eq!(
                prepared.stored_energy_j(&u, cx).unwrap(),
                p.stored_energy_j(&u, cx).unwrap()
            );
            let native = p.tangent_matrix(&u, cx).unwrap();
            for row in 0..12 {
                assert_eq!(tangent.row(row), native.row(row));
            }
        }
    });
}

#[test]
fn interval_loads_are_revalidated_and_use_original_reference_volume() {
    with_cx(|_, cx| {
        let law = material();
        let prepared = problem(&law).prepare(cx).unwrap();
        let nodal = [(0, 4.0), (0, -1.0), (7, 2.0)];
        let force = prepared
            .external_force_n(&nodal, [6.0, -12.0, 18.0], cx)
            .unwrap();
        for component in 0..3 {
            near(
                (0..4).map(|node| force[3 * node + component]).sum(),
                [4.0, 0.0, 3.0][component],
            );
        }
        assert_eq!(
            prepared.external_force_n(&[], [0.0; 3], cx).unwrap(),
            vec![0.0; 12]
        );
        for nodal in [&[(12, 1.0)][..], &[(3, f64::NAN)][..]] {
            assert!(prepared.external_force_n(nodal, [0.0; 3], cx).is_err());
        }
        assert!(
            prepared
                .external_force_n(&[], [f64::INFINITY, 0.0, 0.0], cx)
                .is_err()
        );
        assert!(
            prepared
                .external_force_n(&[(0, f64::MAX), (0, f64::MAX)], [0.0; 3], cx)
                .is_err()
        );
    });
}

#[test]
fn preparation_preserves_inversion_shape_budget_and_cancellation_refusals() {
    with_cx(|gate, cx| {
        let law = material();
        let p = problem(&law);
        let prepared = p.prepare(cx).unwrap();
        let mut inverted = [0.0; 12];
        inverted[3] = -2.0;
        assert!(matches!(
            prepared.internal_force(&inverted, cx),
            Err(HyperTetError::MaterialRefused { .. })
        ));
        assert!(matches!(
            prepared.stored_energy_j(&inverted, cx),
            Err(HyperTetError::MaterialRefused { .. })
        ));
        assert!(matches!(
            prepared.tangent_matrix(&inverted, cx),
            Err(HyperTetError::MaterialRefused { .. })
        ));
        assert!(prepared.internal_force(&[0.0; 11], cx).is_err());
        assert!(prepared.stored_energy_j(&[f64::NAN; 12], cx).is_err());
        let mut over_budget = p;
        over_budget.budget.maximum_nodes = 3;
        assert!(over_budget.prepare(cx).is_err());
        gate.request();
        assert!(matches!(p.prepare(cx), Err(HyperTetError::Cancelled)));
        assert!(matches!(
            prepared.internal_force(&[0.0; 12], cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            prepared.stored_energy_j(&[0.0; 12], cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            prepared.tangent_matrix(&[0.0; 12], cx),
            Err(HyperTetError::Cancelled)
        ));
        assert!(matches!(
            prepared.external_force_n(&[], [0.0; 3], cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}
