//! Independent physical re-solves for the complete reflected-radiation FEM dual.
use super::*;

use fs_conduction::adjoint::RobinResponse;
use fs_conduction::radiation::{
    ConvectiveEnclosureSolution, EnclosureRadiationGradient, EnclosureRadiationLinearization,
    pullback_gray_diffuse_enclosure,
    solve_with_convective_gray_diffuse_enclosure,
};
use fs_conduction::ConductivityTable;

#[derive(Clone, Copy)]
struct Controls {
    density: f64,
    h: f64,
    reference: f64,
    epsilon: f64,
    second_epsilon: f64,
    direction: f64,
}

fn nominal() -> Controls {
    Controls {
        density: 100.0, h: 5.0, reference: 300.0,
        epsilon: 0.8, second_epsilon: 0.6, direction: 0.0,
    }
}

struct State {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
    enclosure: GrayDiffuseEnclosure,
    solved: ConvectiveEnclosureSolution,
}

impl State {
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem {
            mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            element_materials: None, source: &self.source,
        }
    }

    fn weights(&self) -> Vec<f64> {
        let mut weights = vec![0.0; self.mesh.vertex_count()];
        *weights.last_mut().unwrap() = 1.0;
        weights
    }

    fn goal(&self) -> f64 { *self.solved.conduction.temperature.last().unwrap() }

    fn gradient(&self, cx: &fs_exec::Cx<'_>, weights: &[f64]) -> EnclosureRadiationGradient {
        pullback_gray_diffuse_enclosure(
            cx, self.problem(), None, linear(), &self.solved.conduction.temperature,
            &["lower", "upper"], &self.enclosure, 1e-7, weights, 8 * self.mesh.vertex_count(),
        ).unwrap()
    }
}

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-9, ..config().linear }
}

fn direction(position: &[f64; 3]) -> f64 {
    1.0 + 2.0 * position[0] - 3.0 * position[1] + 0.25 * position[2]
}

// Two equal-area patches of the same wall use an explicitly manufactured
// algebraic two-patch enclosure. This is not finite-geometry view-factor data.
// A spatially varying source and smooth k(T) make the full FEM tangent
// nonsymmetric and the wall non-isothermal.
fn physical(cx: &fs_exec::Cx<'_>, c: Controls) -> State {
    physical_with_fixed(cx, c, 0, 400.0)
}

fn physical_with_fixed(
    cx: &fs_exec::Cx<'_>, c: Controls, fixed_axis: usize, fixed_temperature: f64,
) -> State {
    let (complex, positions) = box_grid([3, 2, 2], [1.0; 3]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("fixed", |f| on_box_face(f.centroid[fixed_axis], 0.0),
            ThermalBc::dirichlet(fixed_temperature).unwrap()).unwrap()
        .region("lower", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] < 0.5,
            ThermalBc::robin(c.h, c.reference).unwrap()).unwrap()
        .region("upper", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] >= 0.5,
            ThermalBc::robin(7.0, 310.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let lower = RadiationSurface::new(
        &mesh, "lower", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] < 0.5,
        emissivity("lower", c.epsilon, "manufactured-lower"),
    ).unwrap();
    let upper = RadiationSurface::new(
        &mesh, "upper", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] >= 0.5,
        emissivity("upper", c.second_epsilon, "manufactured-upper"),
    ).unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(lower.area_m2()).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![lower, upper], factors).unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 6.0), (500.0, 16.0)],
    ).unwrap());
    let source = ScalarField::Nodal(mesh.positions().iter()
        .map(|p| c.density + 500.0 * p[1] + c.direction * direction(p)).collect());
    let mut cfg = config();
    cfg.stop.residual_rtol = 1e-12;
    cfg.stop.max_iterations = 80;
    let solved = solve_with_convective_gray_diffuse_enclosure(
        cx, ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source,
        }, None, &enclosure, cfg, CoupledRadiationConfig {
            max_iterations: 256, surface_temperature_atol_k: 1e-10,
            surface_temperature_rtol: 0.0, relaxation: 0.5,
        },
    ).unwrap();
    State { mesh, boundary, material, source, enclosure, solved }
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 8e-5 * expected.abs().max(1e-4),
        "complete enclosure adjoint {actual:e} versus independent physical difference {expected:e}");
}

fn uniform_density_bar(state: &State, lambda: &[f64]) -> f64 {
    state.mesh.complex().tets.iter().enumerate().map(|(e, tet)|
        state.mesh.element_volume(e) / 4.0
            * tet.iter().map(|&v| lambda[v as usize]).sum::<f64>()).sum()
}

#[test]
fn g1_enclosure_dual_matches_source_convection_and_both_emissivity_resolves() {
    with_cx(|cx| {
        let base = physical(cx, nominal());
        let before = base.solved.clone();
        let gradient = base.gradient(cx, &base.weights());
        assert_eq!(base.solved, before, "a pullback must leave the accepted physical solve unchanged");
        assert_eq!(gradient.regions, ["lower", "upper"]);
        assert_eq!(gradient.radiation_regions, ["lower", "upper"]);
        let g = &gradient.convection;
        let density = uniform_density_bar(&base, &g.nodal_load);
        let analytic = [
            density, g.log_htc[0] / nominal().h, g.references[0],
            gradient.emissivities[0], gradient.emissivities[1],
        ];
        for (i, step) in [0.25, 0.001, 0.01, 0.0001, 0.0001].into_iter().enumerate() {
            let perturb = |sign: f64| {
                let mut c = nominal();
                let value = match i {
                    0 => &mut c.density, 1 => &mut c.h, 2 => &mut c.reference,
                    3 => &mut c.epsilon, _ => &mut c.second_epsilon,
                };
                *value += sign * step;
                physical(cx, c).goal()
            };
            close(analytic[i], (perturb(1.0) - perturb(-1.0)) / (2.0 * step));
        }
        let means: Vec<_> = base.enclosure.surfaces().iter().map(|surface|
            surface.mean_temperature(&base.mesh, &base.solved.conduction.temperature).unwrap()).collect();
        assert!((means[0] - means[1]).abs() > 0.1, "reflection derivative needs unequal patch temperatures");
        let wall: Vec<_> = base.mesh.boundary().iter()
            .filter(|face| on_box_face(face.centroid[0], 1.0))
            .flat_map(|face| face.vertices.map(|v| base.solved.conduction.temperature[v as usize])).collect();
        assert!(wall.iter().copied().fold(f64::NEG_INFINITY, f64::max)
            - wall.iter().copied().fold(f64::INFINITY, f64::min) > 0.1);
        let frozen = RobinResponse::pullback_mean_htc(
            cx, ConductionProblem { boundary: &base.solved.combined_boundary, ..base.problem() },
            None, linear(), &base.solved.conduction.temperature, &["lower", "upper"],
            &[0.0; 2], &base.weights(), 4 * base.mesh.vertex_count(),
        ).unwrap();
        assert!((g.references[0] - frozen.references[0]).abs() > 1e-4 * g.references[0].abs(),
            "freezing reflected radiation must change the physical derivative");
        assert!(g.relative_residual < linear().tolerance);
        for &v in fs_conduction::DofMap::new(&base.boundary, base.mesh.vertex_count()).unwrap().fixed() {
            assert_eq!(g.nodal_load[v], 0.0);
        }
    });
}

#[test]
fn g1_enclosure_true_transpose_matches_a_signed_goal_and_nonuniform_load_direction() {
    with_cx(|cx| {
        let base = physical(cx, nominal());
        let weights: Vec<_> = base.mesh.positions().iter()
            .map(|p| p[0] - 0.8 * p[1] + 0.3 * p[2] - 0.25).collect();
        let gradient = base.gradient(cx, &weights);
        // Contract the independent consistent tetrahedral source mass with
        // the returned dual. Neither an isothermal patch nor a lumped source
        // can supply this generally signed virtual-work equality.
        let mut adjoint = 0.0;
        for (e, tet) in base.mesh.complex().tets.iter().enumerate() {
            for (a, &v) in tet.iter().enumerate() {
                for (b, &w) in tet.iter().enumerate() {
                    let mass = base.mesh.element_volume(e) / 20.0
                        * if a == b { 2.0 } else { 1.0 };
                    adjoint += gradient.convection.nodal_load[v as usize] * mass
                        * direction(&base.mesh.positions()[w as usize]);
                }
            }
        }
        let step = 0.1;
        let mut lo = nominal(); lo.direction = -step;
        let mut hi = nominal(); hi.direction = step;
        let minus = physical(cx, lo);
        let plus = physical(cx, hi);
        let forward: f64 = weights.iter().zip(&plus.solved.conduction.temperature)
            .zip(&minus.solved.conduction.temperature)
            .map(|((&w, &p), &m)| w * (p - m) / (2.0 * step)).sum();
        close(adjoint, forward);
    });
}

#[test]
fn g0_g4_enclosure_dual_rechecks_physics_and_honors_exact_work_limits() {
    with_cx(|cx| {
        let base = physical(cx, nominal());
        let n = base.mesh.vertex_count();
        let weights = base.weights();
        let run = |cx: &fs_exec::Cx<'_>, temperature: &[f64], cfg, cap, weights: &[f64]| {
            pullback_gray_diffuse_enclosure(
                cx, base.problem(), None, cfg, temperature, &["lower", "upper"],
                &base.enclosure, 1e-7, weights, cap,
            )
        };
        let mut changed = base.solved.conduction.temperature.clone();
        *changed.last_mut().unwrap() += 1.0;
        assert!(run(cx, &changed, linear(), 8 * n, &weights).is_err(),
            "an unchanged-field residual gate must reject a substituted field");
        assert!(run(cx, &base.solved.conduction.temperature, linear(), 8 * n - 1, &weights).is_err());
        assert!(run(cx, &base.solved.conduction.temperature,
            LinearConfig { max_iterations: 1, restart: 1, ..linear() }, 8 * n, &weights).is_err());
        let zero = run(cx, &base.solved.conduction.temperature, linear(), 8 * n, &vec![0.0; n]).unwrap();
        assert_eq!(zero.convection.iterations, 0);
        assert!(zero.convection.nodal_load.iter().all(|&v| v == 0.0));
        assert_eq!(zero.emissivities, [0.0; 2]);
        assert!(pullback_gray_diffuse_enclosure(
            cx, base.problem(), None, linear(), &base.solved.conduction.temperature,
            &["lower", "lower"], &base.enclosure, 1e-7, &weights, 8 * n,
        ).is_err());
        with_cancelled_cx(|cancelled| assert!(matches!(
            run(cancelled, &base.solved.conduction.temperature, linear(), 8 * n, &weights),
            Err(ConductionError::Cancelled { .. })
        )));
    });
}

#[test]
fn enclosure_prescribed_lift_includes_fixed_wall_means_and_direct_goal_terms() {
    with_cx(|cx| {
        // This fixed face intersects the lower radiating patch. Prescribed
        // vertices therefore contribute to both its mean and J_FD radiation.
        let base = physical_with_fixed(cx, nominal(), 1, 400.0);
        let mut weights = base.weights();
        let fixed_goal = base.mesh.positions().iter().position(|p|
            on_box_face(p[0], 1.0) && on_box_face(p[1], 0.0)).unwrap();
        weights[fixed_goal] += 0.7;
        let full = base.gradient(cx, &weights);
        let bound = EnclosureRadiationLinearization::bind(
            cx, base.problem(), linear(), &base.solved.conduction.temperature,
            &["lower", "upper"], &[0.0; 2], &[0.0; 4], &base.enclosure,
            1e-7, 8 * base.mesh.vertex_count(),
        ).unwrap();
        let problem = ConductionProblem { boundary: bound.boundary(), ..base.problem() };
        let lift = |feedback: &[f64]| {
            let bars = RobinResponse::prescribed_temperature_pullback(
                cx, problem, None, &base.solved.conduction.temperature, &weights,
                &full.convection.nodal_load, &["lower", "upper"],
                bound.convective_mean_slopes(), bound.local_reference_slopes(),
                feedback, 65536,
            ).unwrap();
            base.boundary.dirichlet().iter().map(|&(v, _)| bars[v]).sum::<f64>()
        };
        let actual = lift(bound.reference_feedback());
        let step = 0.01;
        let evaluate = |value| {
            let state = physical_with_fixed(cx, nominal(), 1, value);
            weights.iter().zip(&state.solved.conduction.temperature)
                .map(|(&w, &t)| w*t).sum::<f64>()
        };
        close(actual, (evaluate(400.0+step)-evaluate(400.0-step))/(2.0*step));
        assert!((actual-lift(&[0.0; 4])).abs() > 1e-7,
            "fixed vertices inside wall means must carry reflected-radiation feedback");
        assert!(actual > 0.7, "the direct fixed-node objective term must be retained");
    });
}
