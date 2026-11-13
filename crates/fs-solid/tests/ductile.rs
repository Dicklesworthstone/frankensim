//! Ductility on actual continuum updates, using explicit synthetic coefficients.
//! These are mechanics and history oracles, not calibration or fracture claims.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::{CreepState, plastic::NortonJ2State, thermal::expansion::ThermalCreepState};
use fs_material::ductile::{DuctileDomain, DuctileState, TabulatedDuctileLocus};
use fs_material::plastic::chaboche::ChabocheState;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::reduction::{ReducedState, ReductionConfig};
use fs_material::plastic::J2Plasticity;
use fs_material::{IsotropicElastic, SmallStrainLaw, Voigt};
use fs_solid::ductile::{
    DuctileAdvanceError, DuctileHistory, DuctileLawField, DuctileSolidError,
    DuctileStepConditions, DuctileTemperatures, PlasticStrainHistory,
};
use fs_solid::hyper3::HyperTetSettings;
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetSolid,
};
use fs_solid::{
    Formulation, HexAssemblyBudget, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind,
    SmallStrainHexProblem, SmallStrainLawField, SmallStrainPlaneProblem,
    SmallStrainTetProblem, TetAssemblyBudget,
};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
// Axial displacement at node 1, with free lateral contraction at nodes 2 and 3.
const SUPPORTS: [(usize, f64); 10] = [
    (0, 0.0), (1, 0.0), (2, 0.0), (3, 1.0), (4, 0.0),
    (5, 0.0), (6, 0.0), (8, 0.0), (9, 0.0), (10, 0.0),
];

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 13, kernel_id: 67, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}");
}

fn metal() -> J2Plasticity {
    J2Plasticity::new(IsotropicElastic::new(1000.0, 0.3, 0.05).unwrap(), 1.0, 100.0).unwrap()
}

fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-10.0, 10.0],
        strain_rate_s_inv: [1e-12, 100.0],
        temperature_k: [300.0, 300.0],
    }
}

fn locus(capacity: f64) -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-10.0, capacity), (10.0, capacity)], domain()).unwrap()
}

fn stress_sensitive_locus() -> TabulatedDuctileLocus {
    // Exactly epsilon_f = 1 - 0.02 eta on this explicit synthetic support.
    TabulatedDuctileLocus::new(vec![(-10.0, 1.2), (10.0, 0.8)], domain()).unwrap()
}

fn conditions() -> DuctileStepConditions<'static> {
    DuctileStepConditions { duration_s: 1.0, temperatures_k: DuctileTemperatures::Uniform(300.0) }
}

fn tet<'a>(law: &'a J2Plasticity, path: &'a [f64]) -> SmallStrainTetProblem<'a, J2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: &NODES, tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(law),
        prescribed_m: &SUPPORTS, nodal_forces_n: &[], load_path: path,
        budget: TetAssemblyBudget::standard(), settings: IncrementSettings::default(),
    }
}

fn triaxiality(stress: &Voigt) -> f64 {
    let q2 = 0.5 * ((stress[0] - stress[1]).powi(2)
        + (stress[1] - stress[2]).powi(2) + (stress[2] - stress[0]).powi(2))
        + 3.0 * stress[3..].iter().map(|s| s * s).sum::<f64>();
    (stress[0] + stress[1] + stress[2]) / (3.0 * q2.sqrt())
}

#[test]
fn native_tet_closed_form_initiation_refusal_and_refined_retry() {
    with_cx(|_, cx| {
        let law = metal();
        let calibration = locus(0.004);
        let virgin = DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 1, 1).unwrap();
        let states = [law.initial_state()];
        let first = virgin.advance_tet(&tet(&law, &[0.004]), &[0.0; 12], &states, conditions(), cx).unwrap();
        let alpha = (1000.0 * 0.004 - 1.0) / 1100.0;
        close(first.mechanics.element_states[0].alpha, alpha, 2e-11);
        close(first.mechanics.element_stress[0][0], 1.0 + 100.0 * alpha, 2e-8);
        close(first.ductility.points()[0].utilization(), alpha / 0.004, 2e-8);
        assert_eq!(virgin.points(), &[DuctileState::virgin()]);
        assert_eq!(states, [law.initial_state()]);

        let old_u = first.mechanics.displacement_m.clone();
        let old_states = first.mechanics.element_states.clone();
        let old_damage = first.ductility.points().to_vec();
        let refused = first.ductility.advance_tet(
            &tet(&law, &[0.007]), old_u.as_flattened(), &old_states, conditions(), cx,
        ).unwrap_err();
        match refused {
            DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { point, candidate }) => {
                assert_eq!(point, 0);
                assert!(candidate.state.initiated());
                let fraction = candidate.initiation_fraction.unwrap();
                assert!(fraction > 0.0 && fraction < 1.0);
            }
            other => panic!("expected initiation, got {other:?}"),
        }
        assert_eq!(first.mechanics.displacement_m, old_u);
        assert_eq!(first.mechanics.element_states, old_states);
        assert_eq!(first.ductility.points(), old_damage.as_slice());
        let retry = first.ductility.advance_tet(
            &tet(&law, &[0.005]), old_u.as_flattened(), &old_states, conditions(), cx,
        ).unwrap();
        close(retry.mechanics.element_states[0].alpha, (5.0 - 1.0) / 1100.0, 2e-11);
        close(retry.ductility.points()[0].utilization(), (5.0 - 1.0) / (1100.0 * 0.004), 2e-8);
        assert!(!retry.ductility.points()[0].initiated());
    });
}

#[test]
fn native_unload_reverse_and_restart_keep_irreversible_history_bitwise() {
    with_cx(|_, cx| {
        let law = metal();
        let calibration = locus(0.1);
        let virgin = DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 1, 1).unwrap();
        let first = virgin.advance_tet(
            &tet(&law, &[0.004]), &[0.0; 12], &[law.initial_state()], conditions(), cx,
        ).unwrap();
        // This decrement is elastic; a complete return to zero strain would
        // already cause reverse yielding for these synthetic coefficients.
        let unload = first.ductility.advance_tet(
            &tet(&law, &[0.003]), first.mechanics.displacement_m.as_flattened(),
            &first.mechanics.element_states, conditions(), cx,
        ).unwrap();
        assert_eq!(unload.ductility.points(), first.ductility.points());
        assert_eq!(unload.mechanics.element_states[0].alpha, first.mechanics.element_states[0].alpha);
        assert!(unload.mechanics.element_stress[0][0] < first.mechanics.element_stress[0][0]);

        let restart = unload.ductility.clone();
        let path = [-0.004];
        let problem = tet(&law, &path);
        let reverse = unload.ductility.advance_tet(
            &problem, unload.mechanics.displacement_m.as_flattened(),
            &unload.mechanics.element_states, conditions(), cx,
        ).unwrap();
        let replay = restart.advance_tet(
            &problem, unload.mechanics.displacement_m.as_flattened(),
            &unload.mechanics.element_states, conditions(), cx,
        ).unwrap();
        assert!(reverse.mechanics.element_stress[0][0] < 0.0);
        assert!(reverse.mechanics.element_states[0].alpha > unload.mechanics.element_states[0].alpha);
        assert!(reverse.ductility.points()[0].utilization() > unload.ductility.points()[0].utilization());
        assert_eq!(reverse.mechanics.displacement_m, replay.mechanics.displacement_m);
        assert_eq!(reverse.mechanics.element_states, replay.mechanics.element_states);
        assert_eq!(reverse.mechanics.element_stress, replay.mechanics.element_stress);
        assert_eq!(reverse.ductility.points(), replay.ductility.points());
        assert_eq!(restart.points(), unload.ductility.points());
    });
}

#[test]
fn skipped_history_multiple_targets_and_cancellation_refuse_without_commit() {
    with_cx(|gate, cx| {
        let law = metal();
        let calibration = locus(0.1);
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), 1, 1).unwrap();
        let mut wrong = law.initial_state();
        wrong.alpha = 0.001;
        assert!(matches!(history.advance_tet(
            &tet(&law, &[0.004]), &[0.0; 12], &[wrong], conditions(), cx,
        ), Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
        let states = [law.initial_state()];
        assert!(matches!(history.advance_tet(
            &tet(&law, &[0.002, 0.004]), &[0.0; 12], &states, conditions(), cx,
        ), Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))));
        gate.request();
        assert!(matches!(history.advance_tet(
            &tet(&law, &[0.004]), &[0.0; 12], &states, conditions(), cx,
        ), Err(DuctileAdvanceError::Ductility(DuctileSolidError::Cancelled))));
        assert_eq!(states, [law.initial_state()]);
        assert_eq!(history.points(), &[DuctileState::virgin()]);
    });
}

#[test]
fn plane_reductions_use_full_physical_stress_and_native_point_histories() {
    with_cx(|_, cx| {
        let law = metal();
        let calibration = stress_sensitive_locus();
        for mesh in [Mesh2::quads(1.0, 1.0, 1, 1), Mesh2::triangles(1.0, 1.0, 1, 1)] {
            let prescribed: Vec<_> = mesh.nodes.iter().enumerate()
                .flat_map(|(node, x)| [(2 * node, 0.004 * x[0]), (2 * node + 1, 0.0)])
                .collect();
            // 2x2 Gauss points per quadrilateral; the degree-2 rule's three
            // edge midpoints per triangle (mesh2::quad_points).
            let count = if mesh.elems[0].len() == 4 {
                4 * mesh.elems.len()
            } else {
                3 * mesh.elems.len()
            };
            for plane in [PlaneKind::Stress, PlaneKind::Strain] {
                let states = vec![ReducedState { strain: [0.0; 6], history: law.initial_state() }; count];
                let history = DuctileHistory::virgin(DuctileLawField::Uniform(&calibration), count, count).unwrap();
                let problem = SmallStrainPlaneProblem {
                    mesh: &mesh, laws: SmallStrainLawField::Uniform(&law), plane,
                    formulation: Formulation::Standard, thickness_m: 0.2,
                    reduction: ReductionConfig { max_iterations: 40, stress_tolerance_pa: 1e-12, strain_limit: 0.02 },
                    rotation_limit_rad: 0.02, prescribed_m: &prescribed,
                    nodal_forces_n: &[], load_path: &[1.0],
                    budget: PlaneAssemblyBudget::default(), settings: IncrementSettings::default(),
                };
                let next = history.advance_plane(
                    &problem, &vec![0.0; 2 * mesh.nodes.len()], &states, conditions(), cx,
                ).unwrap();
                assert_eq!(next.ductility.points().len(), count);
                for point in 0..count {
                    let state = &next.mechanics.quadrature_states[point];
                    let stress = next.mechanics.quadrature_stress_pa[point];
                    let eta = triaxiality(&stress);
                    let expected = state.history.alpha / (1.0 - 0.02 * eta);
                    close(next.ductility.points()[point].utilization(), expected, 2e-11);
                    assert_eq!(next.ductility.points()[point].accumulated_plastic_strain(), state.history.alpha);
                    match plane {
                        PlaneKind::Stress => {
                            close(stress[2], 0.0, 2e-9);
                            assert!(state.strain[2] < 0.0);
                        }
                        PlaneKind::Strain => {
                            assert_eq!(state.strain[2], 0.0);
                            assert!(stress[2] > 0.1);
                            let flattened = [stress[0], stress[1], 0.0, stress[3], 0.0, 0.0];
                            let wrong = state.history.alpha / (1.0 - 0.02 * triaxiality(&flattened));
                            assert!((expected - wrong).abs() > 1e-7);
                        }
                    }
                }
                assert!(history.points().iter().all(|p| *p == DuctileState::virgin()));
                assert!(states.iter().all(|s| s.history == law.initial_state()));
            }
        }
    });
}

#[test]
fn bbar_hex_gates_all_eight_points_and_localizes_a_single_weak_calibration() {
    with_cx(|_, cx| {
        let nodes = [
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
        ];
        let cells = [[0, 1, 2, 3, 4, 5, 6, 7]];
        let prescribed: Vec<_> = nodes.iter().enumerate().flat_map(|(n, x)| [
            (3 * n, 0.004 * x[0]), (3 * n + 1, 0.0), (3 * n + 2, 0.0),
        ]).collect();
        let law = metal();
        let problem = SmallStrainHexProblem {
            nodes_m: &nodes, hexahedra: &cells, laws: SmallStrainLawField::Uniform(&law),
            formulation: Formulation::BBar, strain_limit: 0.02, rotation_limit_rad: 0.02,
            prescribed_m: &prescribed, nodal_forces_n: &[], load_path: &[1.0],
            budget: HexAssemblyBudget::default(), settings: IncrementSettings::default(),
        };
        let states = vec![law.initial_state(); 8];
        let temperatures = [300.0];
        let step = DuctileStepConditions {
            duration_s: 1.0, temperatures_k: DuctileTemperatures::PerHexElement(&temperatures),
        };
        let strong = locus(1.0);
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&strong), 8, 8).unwrap();
        let next = history.advance_hex(&problem, &[0.0; 24], &states, step, cx).unwrap();
        assert_eq!(next.mechanics.element_offsets, [0, 8]);
        let shear = 1000.0 / (2.0 * 1.3);
        let alpha = (2.0 * shear * 0.004 - 1.0) / (3.0 * shear + 100.0);
        for point in 0..8 {
            close(next.mechanics.quadrature_states[point].alpha, alpha, 2e-11);
            close(next.ductility.points()[point].utilization(), alpha, 2e-11);
        }
        let mut calibrations = vec![locus(1.0); 8];
        calibrations[5] = locus(1e-4);
        let weak = DuctileHistory::virgin(DuctileLawField::PerPoint(&calibrations), 8, 8).unwrap();
        assert!(matches!(weak.advance_hex(&problem, &[0.0; 24], &states, step, cx),
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { point: 5, .. }))));
        assert!(weak.points().iter().all(|p| *p == DuctileState::virgin()));
        assert!(states.iter().all(|s| *s == law.initial_state()));
    });
}

fn finite_metal() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: 80.0, bulk_modulus_pa: 160.0,
            yield_stress_pa: 1.0, hardening_modulus_pa: 8.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25, maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0, maximum_accumulated_plastic_strain: 3.0,
        },
    ).unwrap()
}

#[test]
fn finite_rotated_stretch_uses_spatial_stress_and_refuses_without_mutation() {
    with_cx(|_, cx| {
        let law = finite_metal();
        let solid = FinitePlasticTetSolid::try_new(
            &NODES, &TETS, FinitePlasticLawField::Uniform(&law), TetAssemblyBudget::standard(), cx,
        ).unwrap();
        let stretch = 1.15_f64;
        let lateral = stretch.sqrt().recip();
        let f = [0.0, -lateral, 0.0, stretch, 0.0, 0.0, 0.0, 0.0, lateral];
        let prescribed: Vec<_> = NODES.iter().enumerate().flat_map(|(node, x)| {
            (0..3).map(move |i| {
                (3 * node + i, (0..3).map(|j| f[3 * i + j] * x[j]).sum::<f64>() - x[i])
            })
        }).collect();
        let controls = FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-8, maximum_endpoint_work_defect_j: 10.0,
            yield_tolerance_pa: 1e-8,
        };
        let strong = stress_sensitive_locus();
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&strong), 1, 1).unwrap();
        let next = history.advance_finite_tet(
            &solid, &prescribed, &[], [0.0; 3], conditions(),
            HyperTetSettings::default(), controls, cx,
        ).unwrap();
        let alpha = (3.0 * 80.0 * stretch.ln() - 1.0) / (3.0 * 80.0 + 8.0);
        close(next.mechanics.solid.element_states()[0].alpha, alpha, 2e-10);
        // Isochoric principal stretches have zero mean Kirchhoff stress before
        // and after rotation. A triaxiality-dependent locus detects raw-Piola use.
        close(next.ductility.points()[0].utilization(), alpha, 2e-10);
        for (actual, expected) in next.mechanics.solid.element_deformation()[0].iter().zip(f) {
            close(*actual, expected, 2e-10);
        }
        let weak = locus(0.01);
        let refusal = DuctileHistory::virgin(DuctileLawField::Uniform(&weak), 1, 1).unwrap();
        assert!(matches!(refusal.advance_finite_tet(
            &solid, &prescribed, &[], [0.0; 3], conditions(),
            HyperTetSettings::default(), controls, cx,
        ), Err(DuctileAdvanceError::Ductility(DuctileSolidError::Initiated { point: 0, .. }))));
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        assert_eq!(solid.element_states(), &[law.initial_state()]);
        assert_eq!(history.points(), &[DuctileState::virgin()]);
        assert_eq!(refusal.points(), &[DuctileState::virgin()]);
    });
}

#[test]
fn plastic_history_accessors_exclude_independent_creep_accumulation() {
    let law = metal();
    let plastic = law.update_state(&[0.004, 0.0, 0.0, 0.0, 0.0, 0.0], &law.initial_state());
    let expected = plastic.alpha;
    assert!(expected > 0.0);
    let combined = NortonJ2State {
        creep: CreepState { creep_strain: [0.0; 6], equivalent_creep_strain: 0.75 },
        plastic: plastic.clone(),
    };
    assert_eq!(combined.plastic_strain_accumulation(), expected);
    let reduced = ReducedState {
        strain: [0.0; 6], history: ThermalCreepState::from_parts(combined, 300.0),
    };
    assert_eq!(reduced.plastic_strain_accumulation(), expected);
    let cyclic = ChabocheState::<2> { plastic, backstress_pa: [[0.0; 6]; 2] };
    assert_eq!(cyclic.plastic_strain_accumulation(), expected);
    let finite_law = finite_metal();
    let stretch = 1.15_f64;
    let lateral = stretch.sqrt().recip();
    let f = [stretch, 0.0, 0.0, 0.0, lateral, 0.0, 0.0, 0.0, lateral];
    let mut finite = finite_law.response(&f, &finite_law.initial_state()).unwrap().state;
    let plastic_alpha = finite.alpha;
    finite.creep_alpha = 0.75;
    assert_eq!(finite.plastic_strain_accumulation(), plastic_alpha);
    assert!(plastic_alpha > 0.0 && plastic_alpha < finite.creep_alpha);
}
