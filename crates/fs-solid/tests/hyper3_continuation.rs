//! Native 3-D finite-strain paths: independent constitutive equilibria,
//! the original force/tangent, physical normalization and transactional refusal.

use std::cell::Cell;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::MaterialError;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::continuation::{self, ArcSettings, PathEvent, PathResidual, PathState};
use fs_solid::hyper3::continuation::{HyperTetContinuationError, HyperTetPathScales};
use fs_solid::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField, SolidError,
    TetAssemblyBudget,
};
use fs_sparse::{Coo, Csr};

#[path = "hyper3_continuation/loads.rs"]
mod loads;

const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const DIAGONAL_SUPPORTS: [(usize, f64); 9] = [
    (0, 0.0),
    (1, 0.0),
    (2, 0.0),
    (4, 0.0),
    (5, 0.0),
    (6, 0.0),
    (8, 0.0),
    (9, 0.0),
    (10, 0.0),
];

fn nodes(length: f64) -> [[f64; 3]; 4] {
    [
        [0.0; 3],
        [length, 0.0, 0.0],
        [0.0, length, 0.0],
        [0.0, 0.0, length],
    ]
}

fn with_cx<R>(run: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 739,
                kernel_id: 318,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&gate, &cx)
    })
}

fn problem<'a, L: HyperelasticMaterial>(
    nodes: &'a [[f64; 3]],
    material: &'a L,
    supports: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
) -> HyperTetProblem<'a, L> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: &TETS,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: supports,
        nodal_forces_n: forces,
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings {
            load_steps: 8,
            relative_tolerance: 1e-11,
            ..HyperTetSettings::default()
        },
    }
}

fn settings() -> ArcSettings {
    ArcSettings {
        ds: 0.04,
        ds_max: 0.08,
        tol: 2e-11,
        ..ArcSettings::default()
    }
}

fn same_state(a: &PathState, b: &PathState) {
    assert_eq!(a.u, b.u);
    assert_eq!(a.lambda, b.lambda);
    assert_eq!(a.tangent, b.tangent);
    assert_eq!(a.ds, b.ds);
    assert_eq!(a.step, b.step);
    assert_eq!(a.events, b.events);
    assert_eq!(a.trace, b.trace);
    assert_eq!(a.last_definite, b.last_definite);
    assert_eq!(a.pending_switch, b.pending_switch);
}

#[test]
fn native_path_matches_uniaxial_neo_hookean_contraction_energy_and_static_owner() {
    with_cx(|_, cx| {
        let x = nodes(1.0);
        let (mu, lame) = (3.2, 8.0);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda: lame }, 5.0).unwrap();
        let reference = mu / 6.0;
        // A load at an eliminated support contributes to reactions, never the path metric.
        let forces = [(3, reference), (0, 3.5 * reference)];
        let p = problem(&x, &material, &DIAGONAL_SUPPORTS, &forces);
        let path = p
            .prepare_continuation(
                HyperTetPathScales {
                    displacement_m: 1.0,
                    force_n: reference,
                },
                cx,
            )
            .unwrap();
        assert_eq!(path.free_dofs(), &[3, 7, 11]);
        let controls = settings();
        let initial = path.start(&controls, cx).unwrap();
        let result = path.advance(&initial, &controls, 18, cx).unwrap();
        let response = path.response(&result, cx).unwrap();
        let s = 1.0 + response.displacement_m[1][0];
        let t = 1.0 + response.displacement_m[2][1];
        assert!(s > 1.15 && t < 1.0);
        assert!((t - 1.0 - response.displacement_m[3][2]).abs() < 1e-11);
        let j = s * t * t;
        assert!((mu * (t * t - 1.0) + lame * j.ln()).abs() < 1e-9);
        let nominal_stress = mu * (s - 1.0 / s) + lame * j.ln() / s;
        assert!((nominal_stress / mu - result.lambda).abs() < 1e-9);
        let energy = (0.5 * mu * (s * s + 2.0 * t * t - 3.0) - mu * j.ln()
            + 0.5 * lame * j.ln().powi(2))
            / 6.0;
        assert!((response.strain_energy_j - energy).abs() < 1e-12);
        let reaction: f64 = response
            .reactions_n
            .iter()
            .filter(|(dof, _)| dof % 3 == 0)
            .map(|(_, f)| f)
            .sum();
        assert!((reaction + result.lambda * 4.5 * reference).abs() < 1e-10);
        let loaded_forces = forces.map(|(dof, force)| (dof, result.lambda * force));
        let native = HyperTetProblem {
            nodal_forces_n: &loaded_forces,
            ..p
        }
        .solve(cx)
        .unwrap();
        for (a, b) in native
            .displacement_m
            .iter()
            .flatten()
            .zip(response.displacement_m.iter().flatten())
        {
            assert!((a - b).abs() < 2e-9);
        }
        let half = path.advance(&initial, &controls, 7, cx).unwrap();
        let resumed = path.advance(&half, &controls, 11, cx).unwrap();
        same_state(&result, &resumed);
        same_state(&initial, &path.start(&controls, cx).unwrap());
        let full: Vec<_> = response.displacement_m.iter().flatten().copied().collect();
        let restarted = path
            .start_from(&full, result.lambda, &controls, cx)
            .unwrap();
        assert_eq!(restarted.u, result.u);
        assert_eq!(restarted.lambda, result.lambda);
    });
}

#[test]
fn normalized_path_and_native_tangent_preserve_force_and_length_scales() {
    with_cx(|_, cx| {
        let mut reference: Option<PathState> = None;
        for (scale, length) in [(1.0, 1.0), (1e-120, 1e-3), (1e120, 1e3)] {
            let x = nodes(length);
            let material = Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 3.2 * scale,
                    lambda: 8.0 * scale,
                },
                5.0,
            )
            .unwrap();
            let force = 3.2 * scale * length * length / 6.0;
            let loads = [(3, force)];
            let p = problem(&x, &material, &DIAGONAL_SUPPORTS, &loads);
            let path = p
                .prepare_continuation(
                    HyperTetPathScales {
                        displacement_m: length,
                        force_n: force,
                    },
                    cx,
                )
                .unwrap();
            let controls = settings();
            let state = path
                .advance(&path.start(&controls, cx).unwrap(), &controls, 12, cx)
                .unwrap();
            if let Some(expected) = &reference {
                assert!((state.lambda - expected.lambda).abs() < 2e-10);
                for (a, b) in state.u.iter().zip(&expected.u) {
                    assert!((a - b).abs() < 2e-10);
                }
            } else {
                reference = Some(state.clone());
            }
            // Independent numerical derivative of the normalized physical residual.
            let k = path.tangent(&state.u, state.lambda, cx).unwrap();
            let direction = [0.7, -0.2, 0.3];
            let h = 2e-6;
            let plus: Vec<_> = state
                .u
                .iter()
                .zip(direction)
                .map(|(q, d)| q + h * d)
                .collect();
            let minus: Vec<_> = state
                .u
                .iter()
                .zip(direction)
                .map(|(q, d)| q - h * d)
                .collect();
            let rp = path.residual(&plus, state.lambda, cx).unwrap();
            let rm = path.residual(&minus, state.lambda, cx).unwrap();
            let mut action = [0.0; 3];
            k.spmv(&direction, &mut action);
            for i in 0..3 {
                assert!(((rp[i] - rm[i]) / (2.0 * h) - action[i]).abs() < 2e-8);
            }
        }
    });
}

#[test]
fn native_constitutive_limit_point_retains_the_descending_equilibrium_branch() {
    // This is an analytic fold of the ORIGINAL logarithmic Neo-Hookean law at
    // large constrained dilation, not a rubber calibration or stability claim.
    // With F=diag(s,1,1), P=mu(s-1/s)+lame*ln(s)/s; its derivative changes sign.
    with_cx(|_, cx| {
        let x = nodes(1.0);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 1.0,
                lambda: 100.0,
            },
            10.0,
        )
        .unwrap();
        let supports: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let force = 100.0 / 6.0;
        let loads = [(3, force)];
        let p = problem(&x, &material, &supports, &loads);
        let path = p
            .prepare_continuation(
                HyperTetPathScales {
                    displacement_m: 1.0,
                    force_n: force,
                },
                cx,
            )
            .unwrap();
        let controls = ArcSettings {
            ds: 0.08,
            ds_max: 0.12,
            tol: 2e-12,
            ..ArcSettings::default()
        };
        let state = path
            .advance(&path.start(&controls, cx).unwrap(), &controls, 38, cx)
            .unwrap();
        let nominal = |s: f64| s - 1.0 / s + 100.0 * s.ln() / s;
        let slope = |s: f64| 1.0 + (1.0 + 100.0 * (1.0 - s.ln())) / (s * s);
        let (mut lo, mut hi) = (2.5, 3.5);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if slope(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let fold = 0.5 * (lo + hi);
        assert!(
            state
                .events
                .iter()
                .any(|event| matches!(event, PathEvent::LimitPoint { .. }))
        );
        assert!(1.0 + state.u[0] > fold + 0.2);
        assert!(state.lambda < nominal(fold) / 100.0);
        assert!(state.tangent.1 < 0.0);
        for &(lambda, q) in &state.trace {
            assert!((lambda - nominal(1.0 + q) / 100.0).abs() < 1e-10);
        }
        let response = path.response(&state, cx).unwrap();
        assert!(response.free_force_residual_n[0].abs() < force * 1e-10);
        assert!(response.min_det_f > 1.0);
    });
}

#[test]
fn native_domain_input_and_checkpoint_refusals_leave_the_complete_path_unchanged() {
    with_cx(|gate, cx| {
        let x = nodes(1.0);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 3.2,
                lambda: 8.0,
            },
            1.01,
        )
        .unwrap();
        let loads = [(3, 1.0)];
        let p = problem(&x, &material, &DIAGONAL_SUPPORTS, &loads);
        let scales = HyperTetPathScales {
            displacement_m: 1.0,
            force_n: 1.0,
        };
        let path = p.prepare_continuation(scales, cx).unwrap();
        let controls = ArcSettings {
            ds: 0.5,
            ds_min: 0.5,
            ds_max: 0.5,
            ..settings()
        };
        let original = path.start(&controls, cx).unwrap();
        let snapshot = original.clone();
        assert!(matches!(
            path.advance(&original, &controls, 1, cx),
            Err(HyperTetContinuationError::Solid(
                HyperTetError::MaterialRefused { element: 0, .. }
            ))
        ));
        same_state(&original, &snapshot);
        assert!(path.residual(&[0.0; 2], 0.0, cx).is_err());
        assert!(matches!(
            path.residual(&[-2.0, 0.0, 0.0], 0.0, cx),
            Err(HyperTetContinuationError::Solid(
                HyperTetError::MaterialRefused { element: 0, .. }
            ))
        ));
        let mut forged = original.clone();
        forged.lambda = 1.0;
        assert!(path.advance(&forged, &controls, 0, cx).is_err());
        assert!(path.start_from(&[0.0; 11], 0.0, &controls, cx).is_err());
        let mut shifted_support = [0.0; 12];
        shifted_support[0] = 0.01;
        assert!(
            path.start_from(&shifted_support, 0.0, &controls, cx)
                .is_err()
        );
        assert!(
            p.prepare_continuation(
                HyperTetPathScales {
                    force_n: 0.0,
                    ..scales
                },
                cx
            )
            .is_err()
        );
        let mut moved = DIAGONAL_SUPPORTS;
        moved[0].1 = 0.01;
        assert!(
            HyperTetProblem {
                prescribed_m: &moved,
                ..p
            }
            .prepare_continuation(scales, cx)
            .is_err()
        );
        assert!(
            HyperTetProblem {
                nodal_forces_n: &[],
                ..p
            }
            .prepare_continuation(scales, cx)
            .is_err()
        );
        let mut limited = p;
        limited.budget.maximum_nodes = 3;
        assert!(limited.prepare_continuation(scales, cx).is_err());
        gate.request();
        assert!(matches!(
            path.advance(&original, &controls, 1, cx),
            Err(HyperTetContinuationError::Solid(HyperTetError::Cancelled))
        ));
        same_state(&original, &snapshot);
    });
}

struct InterruptingMaterial<'a> {
    native: Hyperelastic,
    gate: &'a CancelGate,
    remaining: Cell<usize>,
}
impl HyperelasticMaterial for InterruptingMaterial<'_> {
    fn hyperelastic(&self) -> &Hyperelastic {
        &self.native
    }
    fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        let remaining = self.remaining.get();
        if remaining == 0 {
            self.gate.request();
        } else {
            self.remaining.set(remaining - 1);
        }
        self.native.piola(f)
    }
}

#[test]
fn cancellation_during_native_constitutive_work_publishes_no_path_step() {
    with_cx(|gate, cx| {
        let x = nodes(1.0);
        let material = InterruptingMaterial {
            native: Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 3.2,
                    lambda: 8.0,
                },
                5.0,
            )
            .unwrap(),
            gate,
            remaining: Cell::new(usize::MAX),
        };
        let loads = [(3, 1.0)];
        let path = problem(&x, &material, &DIAGONAL_SUPPORTS, &loads)
            .prepare_continuation(
                HyperTetPathScales {
                    displacement_m: 1.0,
                    force_n: 1.0,
                },
                cx,
            )
            .unwrap();
        let controls = settings();
        let before = path.start(&controls, cx).unwrap();
        let snapshot = before.clone();
        material.remaining.set(2);
        assert!(matches!(
            path.advance(&before, &controls, 4, cx),
            Err(HyperTetContinuationError::Solid(HyperTetError::Cancelled))
        ));
        same_state(&before, &snapshot);
    });
}

struct InterruptingSparsePath {
    polls: Cell<usize>,
}
impl PathResidual for InterruptingSparsePath {
    fn checkpoint(&self) -> Result<(), SolidError> {
        self.polls.set(self.polls.get() + 1);
        if self.polls.get() >= 15 {
            Err(SolidError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn ndof(&self) -> usize {
        1024
    }
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        let mut r = u.to_vec();
        r[0] -= lambda;
        Ok(r)
    }
    fn tangent(&self, _: &[f64], _: f64) -> Result<Csr, SolidError> {
        let mut coo = Coo::new(1024, 1024);
        for i in 0..1024 {
            coo.push(i, i, 1.0);
        }
        Ok(coo.assemble())
    }
    fn load_vector(&self) -> Vec<f64> {
        let mut f = vec![0.0; 1024];
        f[0] = 1.0;
        f
    }
}

#[test]
fn sparse_bordered_krylov_polls_the_existing_problem_checkpoint() {
    let problem = InterruptingSparsePath {
        polls: Cell::new(0),
    };
    let controls = settings();
    let mut state = PathState::start(1024, &controls);
    let snapshot = state.clone();
    assert_eq!(
        continuation::advance(&problem, &mut state, &controls, 1),
        Err(SolidError::Cancelled)
    );
    assert_eq!(problem.polls.get(), 15);
    same_state(&state, &snapshot);
}
