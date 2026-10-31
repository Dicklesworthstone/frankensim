//! Physical material mass through native filtering, curvature, and shared Newton views.
use std::ops::ControlFlow;

use fs_solid::linear3::{TetElasticMaterial, TetMaterialField};
use fs_topopt::elasticity::second_order::{ComplianceError, ComplianceSolveOptions};
use fs_topopt::pipeline::LoadCase;
use fs_topopt::simp_second_order::SimpComplianceGoal;
use fs_topopt::{
    DensityElasticity, DensityFilter, DesignPipeline, EvaluationStop, SimpParams, SolveBudget,
    SolveControl,
};

fn fixture(
    beta: f64,
    sourced: bool,
) -> (
    DesignPipeline,
    DensityElasticity,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
) {
    let (mesh, mut positions) = fs_feec::kuhn_cube(2);
    for p in &mut positions {
        p[0] *= p[0];
    }
    let materials: Vec<_> = (0..mesh.tets.len())
        .map(|i| {
            TetElasticMaterial::try_new(
                2.0 + (i % 7) as f64,
                1.0 + 0.1 * (i % 3) as f64,
                0.29,
                [i as u8 + 1; 32],
            )
            .unwrap()
        })
        .collect();
    let volumes: Vec<_> = fs_feec::element_geometry(&mesh, &positions)
        .vol_signed
        .iter()
        .map(|v| v.abs())
        .collect();
    let masses = volumes
        .iter()
        .zip(&materials)
        .map(|(v, m)| v * m.density_kg_m3)
        .collect();
    let fixed = |p: [f64; 3]| p[0] == 0.0;
    let el = if sourced {
        DensityElasticity::with_material_field(
            &mesh,
            &positions,
            TetMaterialField::PerElement(&materials),
            &fixed,
        )
        .unwrap()
    } else {
        DensityElasticity::new(&mesh, &positions, 1.0, 0.29, &fixed)
    };
    let mut force = vec![0.0; el.n()];
    for (i, p) in positions.iter().enumerate() {
        if p[0] == 1.0 {
            force[3 * i + 2] = -0.1;
        }
    }
    let pipeline = DesignPipeline {
        filter: DensityFilter::new(&mesh, &positions, 0.15),
        params: SimpParams {
            beta,
            penal: 3.0,
            e_min: 0.01,
            eta: 0.5,
        },
    };
    (pipeline, el, force, volumes, masses)
}

fn options() -> ComplianceSolveOptions {
    ComplianceSolveOptions {
        tolerance: 1e-12,
        max_iterations: 2000,
    }
}
fn close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    let error = a
        .iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let scale = b.iter().map(|v| v.abs()).fold(1.0, f64::max);
    assert!(
        error < tolerance * scale,
        "{error:e} / {scale:e} >= {tolerance:e}"
    );
}

#[test]
fn mass_is_unnormalized_projected_material_even_below_the_modal_void_threshold() {
    let (pipeline, el, _, volumes, masses) = fixture(0.0, true);
    close(el.reference_cell_masses_kg().unwrap(), &masses, 1e-13);
    let total: f64 = masses.iter().sum();
    let rho = vec![0.04; el.cells()];
    let (mass, gradient) = pipeline.mass_and_gradient(&el, &rho);
    let (volume, _) = pipeline.volume_and_gradient(&rho, &volumes);
    close(
        &[mass, gradient.iter().sum(), volume],
        &[0.04 * total, total, 0.04],
        1e-10,
    );
    let modal_mass = total * fs_topopt::eigenfreq::mass_interp(0.04);
    assert!(mass > 10.0 * modal_mass);
    let (zero, _) = pipeline.mass_and_gradient(&el, &vec![0.0; el.cells()]);
    assert_eq!(zero, 0.0);
    assert!(el.moduli.iter().all(|&m| m == 1.0));
}

#[test]
fn heterogeneous_mass_has_the_full_filter_projection_gradient_and_cancels_cleanly() {
    let (pipeline, el, _, volumes, masses) = fixture(3.0, true);
    let rho: Vec<_> = (0..el.cells())
        .map(|i| 0.25 + 0.07 * (i % 7) as f64)
        .collect();
    let (mass, gradient) = pipeline.mass_and_gradient(&el, &rho);
    let (_, projected, _) = pipeline.forward(&rho);
    let expected: f64 = masses.iter().zip(&projected).map(|(w, r)| w * r).sum();
    close(&[mass], &[expected], 1e-12);
    let (volume, _) = pipeline.volume_and_gradient(&rho, &volumes);
    assert!((mass - volume * masses.iter().sum::<f64>()).abs() > 1e-5);
    let eps = 1e-4;
    let mut fd = Vec::with_capacity(rho.len());
    for i in 0..rho.len() {
        let mut plus = rho.clone();
        let mut minus = rho.clone();
        plus[i] += eps;
        minus[i] -= eps;
        fd.push(
            (pipeline.mass_and_gradient(&el, &plus).0 - pipeline.mass_and_gradient(&el, &minus).0)
                / (2.0 * eps),
        );
    }
    close(&gradient, &fd, 3e-6);
    let mut stop = |_| ControlFlow::Break(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut stop);
    assert_eq!(
        pipeline.try_mass_and_gradient(&el, &rho, &mut control),
        Err(EvaluationStop::Cancelled)
    );
    assert_eq!(control.work().linear_solves, 0);
    assert_eq!(control.work().linear_iterations, 0);
    assert_eq!(pipeline.mass_and_gradient(&el, &rho), (mass, gradient));
}

#[test]
fn legacy_native_goals_refuse_mass_including_zero_hessian_directions() {
    let (pipeline, el, force, _, _) = fixture(2.0, false);
    assert!(el.reference_cell_masses_kg().is_none());
    let loads = [LoadCase {
        force: &force,
        weight: 1.0,
    }];
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let goal = SimpComplianceGoal::new(
        &pipeline,
        &el,
        &vec![0.45; el.cells()],
        &loads,
        options(),
        &mut control,
    )
    .unwrap();
    let before = control.work();
    assert!(matches!(
        goal.mass_value_gradient(&mut control),
        Err(ComplianceError::InvalidInput(_))
    ));
    assert!(matches!(
        goal.mass_hessian_vector_product(&vec![0.0; el.cells()], &mut control),
        Err(ComplianceError::InvalidInput(_))
    ));
    assert_eq!(control.work(), before);
}

#[cfg(feature = "newton-design")]
mod newton {
    use super::*;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, TilePool};
    use fs_opt::reverse::ReverseLimits;
    use fs_opt::reverse::physics::{PhysicsBinding, PhysicsError, PhysicsHessian, PhysicsModel};
    use fs_opt::{ConstraintKind, Manifold, ProblemBuilder, ReverseProblem, Sense};
    use fs_qty::Dims;
    use fs_topopt::RobustPipeline;
    use fs_topopt::compliance_newton::{COMPLIANCE_DIMS, ComplianceNewtonBudget};
    use fs_topopt::filter::heaviside;
    use fs_topopt::simp_newton::{MASS_DIMS, SimpNewtonModel};

    fn budget() -> ComplianceNewtonBudget {
        ComplianceNewtonBudget {
            solve_reservations: 20_000,
            linear_iterations: 2_000_000,
        }
    }
    fn with_cx(test: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let pool = TilePool::for_host(1, 37);
        let gate = CancelGate::new_clock_free();
        pool.arena_pool().scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 37,
                    kernel_id: 84,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            test(&gate, &cx);
        });
    }

    #[test]
    fn ordinary_and_robust_mass_share_cache_units_curvature_and_failure_semantics() {
        for (beta, offset) in [(0.0, 0.0), (3.0, 0.0), (3.0, 0.1)] {
            let (pipeline, el, force, volumes, masses) = fixture(beta, true);
            let other = fixture(beta, true).0;
            let robust = RobustPipeline {
                filter: other.filter,
                params: other.params,
                eta_offset: offset,
            };
            let loads = [LoadCase {
                force: &force,
                weight: 1.0,
            }];
            let model = if offset == 0.0 {
                SimpNewtonModel::new(&pipeline, &el, &loads, &volumes, options(), budget())
            } else {
                SimpNewtonModel::new_robust(&robust, &el, &loads, &volumes, options(), budget())
            }
            .unwrap();
            let c = model.compliance("compliance").unwrap();
            let v = model.volume("volume").unwrap();
            let m = model.mass("mass").unwrap();
            let x: Vec<_> = (0..el.cells())
                .map(|i| 0.36 + 0.025 * (i % 9) as f64)
                .collect();
            let d: Vec<_> = (0..el.cells())
                .map(|i| 0.03 * ((i % 5) as f64 - 2.0))
                .collect();
            let cs = c.value_gradient(&x, None).unwrap();
            let prepared = model.work();
            let sample = m.value_gradient(&x, None).unwrap();
            let vs = v.value_gradient(&x, None).unwrap();
            assert_eq!(model.work(), prepared);
            assert_eq!(prepared.solve_reservations, loads.len() + 4);
            assert_eq!(prepared.started_solves, loads.len() + 4);
            assert_eq!(MASS_DIMS, Dims([0, 1, 0, 0, 0, 0]));
            assert_eq!(m.signature().residual_dims, MASS_DIMS);
            assert_eq!(m.signature().parameter_dims, Dims::NONE);
            let filtered = pipeline.filter.apply(&x);
            let expected: f64 = masses
                .iter()
                .zip(&filtered)
                .map(|(w, r)| w * heaviside(*r, beta, 0.5 - offset))
                .sum();
            close(&[sample.value], &[expected], 1e-11);
            if offset > 0.0 {
                let eroded: f64 = masses
                    .iter()
                    .zip(&filtered)
                    .map(|(w, r)| w * heaviside(*r, beta, 0.5 + offset))
                    .sum();
                assert!(sample.value > eroded);
            }
            let mut callback = |_| ControlFlow::Continue(());
            let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
            let native = if offset == 0.0 {
                SimpComplianceGoal::new(&pipeline, &el, &x, &loads, options(), &mut control)
            } else {
                SimpComplianceGoal::new_robust(&robust, &el, &x, &loads, options(), &mut control)
            }
            .unwrap();
            let (value, gradient) = native.mass_value_gradient(&mut control).unwrap();
            assert_eq!((sample.value, &sample.gradient), (value, &gradient));
            let hv = m.hessian_vector_product(&x, &d, None).unwrap();
            assert_eq!(
                hv,
                native
                    .mass_hessian_vector_product(&d, &mut control)
                    .unwrap()
            );
            assert_eq!(
                model.work().solve_reservations,
                prepared.solve_reservations + if beta == 0.0 { 0 } else { 2 }
            );
            let spent = model.work();
            assert_eq!(
                m.hessian_vector_product(&x, &vec![0.0; x.len()], None)
                    .unwrap(),
                vec![0.0; x.len()]
            );
            assert_eq!(model.work(), spent);
            let eps = 1e-4;
            let plus: Vec<_> = x.iter().zip(&d).map(|(x, d)| x + eps * d).collect();
            let minus: Vec<_> = x.iter().zip(&d).map(|(x, d)| x - eps * d).collect();
            let gp = m.value_gradient(&plus, None).unwrap().gradient;
            let gm = m.value_gradient(&minus, None).unwrap().gradient;
            let fd: Vec<_> = gp
                .iter()
                .zip(&gm)
                .map(|(p, m)| (p - m) / (2.0 * eps))
                .collect();
            close(&hv, &fd, 3e-6);
            assert_eq!(m.value_gradient(&x, None).unwrap(), sample);
            assert_eq!(c.value_gradient(&x, None).unwrap(), cs);
            assert_eq!(v.value_gradient(&x, None).unwrap(), vs);
            // Real IR: dimensionless minimum mass with a normalized compliance cap.
            let total: f64 = masses.iter().sum();
            let mut b = ProblemBuilder::new();
            let var = b
                .var(
                    "rho",
                    Manifold::Rn {
                        dim: x.len() as u32,
                    },
                    Dims::NONE,
                )
                .unwrap();
            let mn = b.pde_residual("mass", var, true, MASS_DIMS).unwrap();
            let cn = b
                .pde_residual("compliance", var, true, COMPLIANCE_DIMS)
                .unwrap();
            let ms = b.konst(total, MASS_DIMS).unwrap();
            let objective = b.div(mn, ms).unwrap();
            b.objective(objective, Sense::Minimize, 1.0).unwrap();
            let scale = b.konst(cs.value, COMPLIANCE_DIMS).unwrap();
            let normalized = b.div(cn, scale).unwrap();
            let cap = b.konst(1.1, Dims::NONE).unwrap();
            let residual = b.sub(normalized, cap).unwrap();
            b.constraint(residual, ConstraintKind::LeZero, "compliance cap")
                .unwrap();
            let problem = b.finish();
            let oracle = ReverseProblem::new_with_physics(
                &problem,
                ReverseLimits {
                    max_nodes: 128,
                    max_scalar_slots: 8192,
                },
                &[
                    PhysicsBinding {
                        node: mn,
                        model: &m,
                    },
                    PhysicsBinding {
                        node: cn,
                        model: &c,
                    },
                ],
                None,
            )
            .unwrap();
            let before_ir = model.work();
            with_cx(|_, cx| {
                let tape = oracle.evaluate_cancellable(&[x.clone()], cx).unwrap();
                let gradient = tape.lagrangian_gradient_cancellable(&[0.7], cx).unwrap();
                let expected: Vec<_> = sample
                    .gradient
                    .iter()
                    .zip(&cs.gradient)
                    .map(|(m, c)| m / total + 0.7 * c / cs.value)
                    .collect();
                close(&gradient[0], &expected, 1e-10);
                assert_eq!(model.work(), before_ir);
                let ir_hv = tape
                    .lagrangian_hessian_vector_product(&[0.7], &[d.clone()], Some(cx))
                    .unwrap();
                let ch = c.hessian_vector_product(&x, &d, Some(cx)).unwrap();
                let expected: Vec<_> = hv
                    .iter()
                    .zip(&ch)
                    .map(|(m, c)| m / total + 0.7 * c / cs.value)
                    .collect();
                close(&ir_hv[0], &expected, 5e-9);
            });
            let accepted = model.work();
            model.set_budget(ComplianceNewtonBudget {
                solve_reservations: accepted.solve_reservations,
                ..budget()
            });
            assert!(matches!(
                m.value_gradient(&plus, None),
                Err(PhysicsError::Budget { .. })
            ));
            with_cx(|gate, cx| {
                gate.request();
                assert_eq!(
                    m.value_gradient(&plus, Some(cx)),
                    Err(PhysicsError::Cancelled)
                );
                assert_eq!(
                    m.hessian_vector_product(&x, &vec![0.0; x.len()], Some(cx)),
                    Err(PhysicsError::Cancelled)
                );
            });
            assert_eq!(model.work(), accepted);
            assert_eq!(m.value_gradient(&x, None).unwrap(), sample);
            assert!(el.moduli.iter().all(|&m| m == 1.0));
        }
    }

    #[test]
    fn legacy_newton_mass_refusal_spends_no_work() {
        let (pipeline, el, force, volumes, _) = fixture(2.0, false);
        let model = SimpNewtonModel::new(
            &pipeline,
            &el,
            &[LoadCase {
                force: &force,
                weight: 1.0,
            }],
            &volumes,
            options(),
            budget(),
        )
        .unwrap();
        assert!(model.mass("mass").is_err());
        assert_eq!(model.work().solve_reservations, 0);
        assert_eq!(model.work().started_solves, 0);
    }
}
