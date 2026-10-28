//! G1/G3 mean-flux response with no artificial Robin ports. Fresh physical
//! fixed points exercise the shared steady, capacity and prescribed-lift joins.
mod support;

use fs_conduction::adjoint::{MeanFluxFeedback, RobinResponse};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig, StepSolution};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    InitialGuess, LinearConfig, ScalarField, SolveConfig, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder};
use fs_exec::Cx;
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([2, 1, 1], [1.0, 1.0, 1.0]);
    ConductionMesh::new(complex, positions).unwrap()
}

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-11, max_iterations: 1000, restart: 24 }
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}");
}

fn weights(mesh: &ConductionMesh) -> Vec<f64> {
    mesh.positions().iter().map(|p| 0.1 + 0.3 * p[0] + 0.2 * p[1] + 0.4 * p[2]).collect()
}

fn mean(mesh: &ConductionMesh, boundary: &ThermalBoundary, field: &[f64], name: &str) -> f64 {
    let region = boundary.region_names().iter().position(|found| found == name).unwrap();
    let mut area = 0.0;
    let mut integral = 0.0;
    for (slot, face) in mesh.boundary().iter().enumerate() {
        if boundary.region_for(slot) == Some(region) {
            area += face.area;
            integral += face.area * face.vertices.iter().map(|&v| field[v as usize]).sum::<f64>() / 3.0;
        }
    }
    integral / area
}

const STEADY_SLOPES: [f64; 4] = [6.0, -2.0, -1.0, 4.0];

fn steady_boundary(mesh: &ConductionMesh, flux: [f64; 2], bath: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("bath", |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(bath).unwrap()).unwrap()
        .region("east", |face| on_box_face(face.centroid[0], 1.0),
            ThermalBc::neumann(flux[0]).unwrap()).unwrap()
        .region("north", |face| on_box_face(face.centroid[1], 1.0),
            ThermalBc::neumann(flux[1]).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

// Independently iterate the actual affine outward-flux VALUES through the
// existing forward solver. The new feedback or adjoint code is never called.
fn steady(
    cx: &Cx<'_>, mesh: &ConductionMesh, material: &ConductivityModel,
    control: [f64; 3],
) -> (ThermalBoundary, Vec<f64>) {
    let source = ScalarField::Uniform(20.0 + control[0]);
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(300.0);
    config.linear = linear();
    config.stop.residual_rtol = 1e-12;
    config.stop.step_atol = 0.0;
    let mut flux = [0.0; 2];
    for _ in 0..160 {
        let boundary = steady_boundary(mesh, flux, 300.0 + control[2]);
        let solved = fs_conduction::solve(cx, ConductionProblem {
            mesh, boundary: &boundary, material, element_materials: None, source: &source,
        }, config).unwrap();
        let means = [mean(mesh, &boundary, &solved.temperature, "east") - 300.0,
            mean(mesh, &boundary, &solved.temperature, "north") - 300.0];
        let next = [STEADY_SLOPES[0] * means[0] + STEADY_SLOPES[1] * means[1] + control[1],
            STEADY_SLOPES[2] * means[0] + STEADY_SLOPES[3] * means[1]];
        if next.iter().zip(flux).all(|(new, old)| (new - old).abs() <= 1e-11) {
            return (boundary, solved.temperature);
        }
        for index in 0..2 { flux[index] = 0.5 * (flux[index] + next[index]); }
    }
    panic!("independent steady mean-flux fixed point did not converge");
}

#[test]
fn g1_zero_robin_flux_dual_and_complete_prescribed_lift_match_physical_resolves() {
    with_cx(|cx| {
        let mesh = mesh();
        let material = ConductivityModel::isotropic_declared(8.0).unwrap();
        let (boundary, temperature) = steady(cx, &mesh, &material, [0.0; 3]);
        let source = ScalarField::Uniform(20.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source };
        let n = mesh.vertex_count();
        let feedback = MeanFluxFeedback::new(cx, &mesh, &boundary,
            &["east", "north"], &STEADY_SLOPES, 4 * n).unwrap();
        assert_eq!(feedback.vertex_count(), n);
        assert_eq!(feedback.region_count(), 2);
        assert_eq!(feedback.factor_entries(), 4 * n);
        assert_eq!(feedback.region_names(), &["east".to_string(), "north".to_string()]);
        let weights = weights(&mesh);
        let gradient = RobinResponse::pullback_mean_robin_with_reference_and_flux_feedback(
            cx, problem, None, linear(), &temperature, &[], &[], &[], &[],
            Some(&feedback), &weights, 4 * n,
        ).unwrap();
        assert!(gradient.references.is_empty() && gradient.log_htc.is_empty());
        assert!(gradient.relative_residual < linear().tolerance);
        let source_bar: f64 = mesh.complex().tets.iter().enumerate()
            .map(|(element, tet)| mesh.element_volume(element) / 4.0
                * tet.iter().map(|&v| gradient.nodal_load[v as usize]).sum::<f64>()).sum();
        let flux_bars = feedback.flux_pullback(cx, &gradient.nodal_load).unwrap();
        let lift = RobinResponse::prescribed_temperature_pullback_with_flux_feedback(
            cx, problem, None, &temperature, &weights, &gradient.nodal_load,
            &[], &[], &[], &[], Some(&feedback), 100_000,
        ).unwrap();
        let derivatives = [source_bar, -flux_bars[0], lift.iter().sum()];
        for axis in 0..3 {
            let h = 0.001;
            let mut plus = [0.0; 3]; plus[axis] = h;
            let mut minus = [0.0; 3]; minus[axis] = -h;
            let expected = (dot(&weights, &steady(cx, &mesh, &material, plus).1)
                - dot(&weights, &steady(cx, &mesh, &material, minus).1)) / (2.0 * h);
            close(derivatives[axis], expected, 3e-6);
        }
        // The north patch touches the prescribed bath; dropping its fixed
        // mean weights loses a real boundary response, including direct w_D.
        let missing = RobinResponse::prescribed_temperature_pullback(
            cx, problem, None, &temperature, &weights, &gradient.nodal_load,
            &[], &[], &[], &[], 100_000,
        ).unwrap();
        assert!((lift.iter().sum::<f64>() - missing.iter().sum::<f64>()).abs() > 1e-3);
        for &(vertex, _) in boundary.dirichlet() { assert_eq!(gradient.nodal_load[vertex], 0.0); }
        assert!(MeanFluxFeedback::new(cx, &mesh, &boundary,
            &["east", "north"], &STEADY_SLOPES, 4 * n - 1).is_err());
        assert!(MeanFluxFeedback::new(cx, &mesh, &boundary,
            &["east", "east"], &STEADY_SLOPES, 4 * n).is_err());
        assert!(RobinResponse::pullback_mean_robin_with_reference_and_flux_feedback(
            cx, problem, None, linear(), &temperature, &[], &[], &[], &[],
            Some(&feedback), &vec![0.0; n], 4 * n - 1,
        ).is_err(), "a zero goal must not bypass the combined factor allowance");
        let old = RobinResponse::pullback_mean_robin_with_reference_feedback(
            cx, problem, None, linear(), &temperature, &[], &[], &[], &[], &weights, 0,
        ).unwrap();
        let compatible = RobinResponse::pullback_mean_robin_with_reference_and_flux_feedback(
            cx, problem, None, linear(), &temperature, &[], &[], &[], &[], None, &weights, 0,
        ).unwrap();
        assert_eq!(old, compatible);
        with_cancelled_cx(|cancelled| {
            assert!(MeanFluxFeedback::new(cancelled, &mesh, &boundary,
                &["east", "north"], &STEADY_SLOPES, 4 * n).is_err());
            assert!(feedback.flux_pullback(cancelled, &gradient.nodal_load).is_err());
        });
    });
}

const CAPACITY_SLOPES: [f64; 4] = [4.0, -4.0, -4.0, 4.0];

fn capacity_boundary(mesh: &ConductionMesh, heat: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("west", |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::neumann(heat).unwrap()).unwrap()
        .region("east", |face| on_box_face(face.centroid[0], 1.0),
            ThermalBc::neumann(-heat).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

fn step_config() -> StepConfig {
    StepConfig { linear: linear(), energy_tolerance_j: 1e-8 }
}

fn capacity(
    cx: &Cx<'_>, mesh: &ConductionMesh, material: &ConductivityModel,
    engine: &BackwardEuler<'_>, old: &[f64], source: &ScalarField, flux_shift: f64,
) -> (ThermalBoundary, StepSolution) {
    let mut heat = 0.0;
    for _ in 0..160 {
        let boundary = capacity_boundary(mesh, heat);
        let solved = engine.advance_nonlinear_prescribed(cx, ConductionProblem {
            mesh, boundary: &boundary, material, element_materials: None, source,
        }, None, old, 0.5, step_config(), NonlinearStepConfig::default()).unwrap().step;
        let next = 4.0 * (mean(mesh, &boundary, &solved.temperature, "west")
            - mean(mesh, &boundary, &solved.temperature, "east")) + flux_shift;
        if (next - heat).abs() <= 1e-11 { return (boundary, solved); }
        heat = 0.5 * (heat + next);
    }
    panic!("independent capacity mean-flux fixed point did not converge");
}

#[test]
fn g1_capacity_zero_robin_feedback_retains_k_t_storage_history_and_flux_controls() {
    with_cx(|cx| {
        let mesh = mesh();
        let n = mesh.vertex_count();
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 4.0), (450.0, 8.0)]).unwrap());
        let source = ScalarField::Uniform(20.0);
        let old: Vec<_> = mesh.positions().iter().map(|p| 320.0 + 3.0 * p[1] - 2.0 * p[0]).collect();
        let original = old.clone();
        let direction: Vec<_> = mesh.positions().iter().map(|p| p[0] - 0.4 * p[1] + 0.25).collect();
        let engine = BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(20.0).unwrap()).unwrap();
        let (boundary, accepted) = capacity(cx, &mesh, &material, &engine, &old, &source, 0.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source };
        let step = engine.linearize_prescribed_endpoint(cx, problem, None, &old,
            &accepted.temperature, 0.5, step_config(), 1e-8, &[]).unwrap();
        assert_eq!(step.primal().temperature, accepted.temperature);
        assert_eq!(old, original);
        let feedback = MeanFluxFeedback::new(cx, &mesh, &boundary,
            &["west", "east"], &CAPACITY_SLOPES, 4 * n).unwrap();
        let weights = weights(&mesh);
        let gradient = step.pullback_with_robin_and_flux_feedback(cx, &[], &[], &[],
            Some(&feedback), &weights, 4 * n).unwrap();
        let history = step.previous_temperature_pullback(cx, &gradient.nodal_load).unwrap();
        let derivative = 2.0 * dot(&history, &direction)
            + 0.1 * step.capacity_multiplier_pullback(cx, &gradient.nodal_load).unwrap()
            + 0.7 * step.source_multiplier_pullback(cx, &gradient.nodal_load).unwrap();
        let flux = feedback.flux_pullback(cx, &gradient.nodal_load).unwrap();
        let probe = |delta: f64, flux_shift: f64| {
            let engine = BackwardEuler::uniform(cx, &mesh,
                VolumetricHeatCapacity::declared(20.0 * (1.0 + 0.1 * delta)).unwrap()).unwrap();
            let changed: Vec<_> = old.iter().zip(&direction).map(|(t, d)| t + 2.0 * d * delta).collect();
            let source = ScalarField::Uniform(20.0 * (1.0 + 0.7 * delta));
            dot(&weights, &capacity(cx, &mesh, &material, &engine, &changed, &source, flux_shift).1.temperature)
        };
        let h = 0.001;
        close(derivative, (probe(h, 0.0) - probe(-h, 0.0)) / (2.0 * h), 5e-6);
        close(-(flux[0] - flux[1]), (probe(0.0, h) - probe(0.0, -h)) / (2.0 * h), 5e-6);
        let frozen = step.pullback_with_robin_feedback(cx, &[], &[], &[], &weights, 0).unwrap();
        let frozen_history = step.previous_temperature_pullback(cx, &frozen.nodal_load).unwrap();
        assert!((dot(&history, &direction) - dot(&frozen_history, &direction)).abs() > 1e-4);
        let compatible = step.pullback_with_robin_and_flux_feedback(
            cx, &[], &[], &[], None, &weights, 0,
        ).unwrap();
        assert_eq!(frozen, compatible);
        assert!(step.pullback_with_robin_and_flux_feedback(cx, &[], &[], &[],
            Some(&feedback), &vec![0.0; n], 4 * n - 1).is_err());
    });
}
