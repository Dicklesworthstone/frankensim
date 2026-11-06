//! G0/G1/G4: original DKT elasticity, plastic-section oracle, forming/unloading,
//! consistent tangent, retained history and refusal/cancellation atomicity.

use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::{J2Plasticity, J2State};
use fs_plate::plastic::{BendingConfig, BendingError, BendingSolveControls, PlasticBendingPlate};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateError, PlateMesh, PlateSection, assemble};

const E: f64 = 200e9;
const NU: f64 = 0.3;
const YIELD: f64 = 60e6;
const H: f64 = 0.001;

fn config(cells: usize) -> BendingConfig {
    BendingConfig {
        thickness_m: H,
        thickness_cells: cells,
        reduction: ReductionConfig {
            max_iterations: 32,
            stress_tolerance_pa: 1e-7,
            strain_limit: 0.01,
        },
        maximum_slope: 0.05,
        maximum_dofs: 100,
        maximum_material_points: 10_000,
    }
}
fn controls() -> BendingSolveControls {
    BendingSolveControls {
        maximum_iterations: 32,
        maximum_backtracks: 20,
        force_tolerance_n: 1e-8,
        moment_tolerance_nm: 1e-10,
        maximum_factor_bytes: 1_000_000,
    }
}
fn plate(mesh: &PlateMesh, cells: usize) -> PlasticBendingPlate<J2Plasticity> {
    let material =
        J2Plasticity::new(IsotropicElastic::new(E, NU, 0.01).unwrap(), YIELD, 0.0).unwrap();
    PlasticBendingPlate::new(mesh, material, config(cells), || false).unwrap()
}
fn field(mesh: &PlateMesh, k: [f64; 3]) -> Vec<f64> {
    // DKT's original beta = -grad(w), so this produces engineering curvature k.
    mesh.nodes
        .iter()
        .flat_map(|&(x, y)| {
            [
                -0.5 * (k[0] * x * x + k[1] * y * y + k[2] * x * y),
                -k[0] * x - 0.5 * k[2] * y,
                -k[1] * y - 0.5 * k[2] * x,
            ]
        })
        .collect()
}
fn elastic_matrix(mesh: &PlateMesh) -> fs_sparse::Csr {
    let section = PlateSection::isotropic(E, NU, H, 7800.0).unwrap();
    assemble(
        mesh,
        &section,
        &[],
        &[],
        &AssemblyOptions {
            pretension: 0.0,
            support: EdgeSupport::Clamped,
        },
    )
    .unwrap()
    .k
}
fn peak(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |a, v| a.max(v.abs()))
}
fn yield_curvature() -> f64 {
    2.0 * (YIELD / 3.0_f64.sqrt()) / (E / (1.0 + NU) * H)
}

#[test]
fn g1_elastic_thickness_integration_matches_original_dkt_including_twist() {
    let mesh = PlateMesh::rectangle(0.02, 0.017, 2, 1);
    let k = elastic_matrix(&mesh);
    let reference = k.to_dense();
    let u = field(&mesh, [0.03, -0.02, 0.04]);
    let mut force = vec![0.0; u.len()];
    k.spmv(&u, &mut force);
    for cells in [1, 3, 12] {
        let model = plate(&mesh, cells);
        let response = model
            .evaluate(&model.initial_state(), &u, || false)
            .unwrap();
        let actual = response.tangent.to_dense();
        assert!(
            peak(
                &actual
                    .iter()
                    .zip(&reference)
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>()
            ) < 2e-12 * peak(&reference)
        );
        assert!(
            peak(
                &response
                    .internal_force
                    .iter()
                    .zip(&force)
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>()
            ) < 2e-12 * peak(&force)
        );
        assert!(
            response
                .candidate
                .material_points()
                .iter()
                .all(|p| p.history.alpha == 0.0)
        );
        let integral: f64 = model
            .thickness_points()
            .iter()
            .map(|(z, w)| z * z * w)
            .sum();
        assert!((integral - H.powi(3) / 12.0).abs() < 1e-15 * H.powi(3));
    }
}

#[test]
fn g1_j2_saddle_bending_forms_and_springs_back_under_actual_nodal_loads() {
    let mesh = PlateMesh::rectangle(0.02, 0.02, 1, 1);
    let model = plate(&mesh, 64);
    let ky = yield_curvature();
    let virgin = model.initial_state();
    let below_yield = model
        .evaluate(&virgin, &field(&mesh, [0.99 * ky, -0.99 * ky, 0.0]), || {
            false
        })
        .unwrap();
    assert!(
        below_yield
            .candidate
            .material_points()
            .iter()
            .all(|p| p.history.alpha == 0.0)
    );
    let first_yield = model
        .evaluate(&virgin, &field(&mesh, [1.01 * ky, -1.01 * ky, 0.0]), || {
            false
        })
        .unwrap();
    assert!(
        first_yield.candidate.material_points()[0].history.alpha > 0.0,
        "the actual bottom face yields at the analytical section threshold"
    );
    let maximum_curvature = 3.0 * ky;
    // With k_y=-k_x, sigma_y=-sigma_x and sigma_z=0. Thus q=sqrt(3)|sigma_x|,
    // c=2G, and the elastic-perfectly-plastic section has the independent law
    // M = (Y/sqrt(3)) * (h²/4 - a²/3), a = Y/(sqrt(3)*c*k).
    let c = E / (1.0 + NU);
    let i = H.powi(3) / 12.0;
    let a = YIELD / 3.0_f64.sqrt() / (c * maximum_curvature);
    let moment = YIELD / 3.0_f64.sqrt() * (H * H / 4.0 - a * a / 3.0);
    let residual_curvature = maximum_curvature - moment / (c * i);
    let unit = field(&mesh, [1.0, -1.0, 0.0]);
    let mut load = vec![0.0; unit.len()];
    elastic_matrix(&mesh).spmv(&unit, &mut load);
    for value in &mut load {
        *value *= moment / (c * i);
    }
    let supports = [(0, 0.0), (1, 0.0), (2, 0.0)];
    let mut state = model.initial_state();
    for fraction in [0.25, 0.5, 0.75, 1.0] {
        let applied: Vec<_> = load.iter().map(|f| fraction * f).collect();
        state = model
            .equilibrate(&state, &applied, &supports, controls(), || false)
            .unwrap()
            .response
            .candidate;
    }
    let formed = field(&mesh, [maximum_curvature, -maximum_curvature, 0.0]);
    assert!(
        peak(
            &state
                .displacement()
                .iter()
                .zip(&formed)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>()
        ) < 5e-4 * peak(&formed)
    );
    assert!(
        state
            .material_points()
            .iter()
            .any(|p| p.history.alpha > 0.0)
    );
    let loaded_snapshot = state.clone();
    let unloaded = model
        .equilibrate(
            &state,
            &vec![0.0; unit.len()],
            &supports,
            controls(),
            || false,
        )
        .unwrap();
    let expected = field(&mesh, [residual_curvature, -residual_curvature, 0.0]);
    assert!(
        peak(
            &unloaded
                .response
                .candidate
                .displacement()
                .iter()
                .zip(&expected)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>()
        ) < 1e-3 * peak(&expected)
    );
    assert!(peak(unloaded.response.candidate.displacement()) > 0.4 * peak(state.displacement()));
    assert!(
        unloaded
            .response
            .moments
            .iter()
            .flatten()
            .all(|m| m.abs() < 2e-7)
    );
    assert!(
        unloaded
            .response
            .membrane_resultants
            .iter()
            .flatten()
            .all(|n| n.abs() < 2e-6)
    );
    for (before, after) in state
        .material_points()
        .iter()
        .zip(unloaded.response.candidate.material_points())
    {
        assert!(
            (before.history.alpha - after.history.alpha).abs() < 1e-12,
            "elastic unloading must retain work-hardening/plastic history"
        );
    }
    assert_eq!(state, loaded_snapshot);
}

#[test]
fn g0_algorithmic_dkt_tangent_differentiates_plastic_force_from_frozen_history() {
    let mesh = PlateMesh::rectangle(0.02, 0.02, 1, 1);
    let model = plate(&mesh, 12);
    let old = model.initial_state();
    let ky = yield_curvature();
    let u = field(&mesh, [2.31 * ky, -0.71 * ky, 0.29 * ky]);
    let response = model.evaluate(&old, &u, || false).unwrap();
    assert!(
        response
            .candidate
            .material_points()
            .iter()
            .any(|p| p.history.alpha > 0.0)
    );
    let tangent = response.tangent.to_dense();
    for column in 0..u.len() {
        let step = if column % 3 == 0 { 2e-10 } else { 2e-8 };
        let mut plus = u.clone();
        let mut minus = u.clone();
        plus[column] += step;
        minus[column] -= step;
        let fp = model
            .evaluate(&old, &plus, || false)
            .unwrap()
            .internal_force;
        let fm = model
            .evaluate(&old, &minus, || false)
            .unwrap()
            .internal_force;
        let fd: Vec<_> = fp
            .iter()
            .zip(&fm)
            .map(|(p, m)| (p - m) / (2.0 * step))
            .collect();
        let error: Vec<_> = (0..u.len())
            .map(|r| fd[r] - tangent[r * u.len() + column])
            .collect();
        assert!(
            peak(&error) < 2e-5 * peak(&fd),
            "column {column}: {} / {}",
            peak(&error),
            peak(&fd)
        );
    }
    assert!(old.material_points().iter().all(|p| p.history
        == J2State {
            plastic_strain: [0.0; 6],
            alpha: 0.0
        }));
}

#[test]
fn g4_foreign_history_failed_steps_and_mid_assembly_cancellation_are_atomic() {
    let mesh = PlateMesh::rectangle(0.02, 0.02, 1, 1);
    let model = plate(&mesh, 8);
    let virgin = model.initial_state();
    let u = field(
        &mesh,
        [2.0 * yield_curvature(), -2.0 * yield_curvature(), 0.0],
    );
    let formed = model.evaluate(&virgin, &u, || false).unwrap().candidate;
    let snapshot = formed.clone();
    let foreign = plate(&mesh, 8);
    assert!(matches!(
        foreign.evaluate(&formed, &u, || false),
        Err(BendingError::Invalid(_))
    ));
    let mut checks = 0;
    let cancelled = model.evaluate(&formed, &u, || {
        checks += 1;
        checks == 15
    });
    assert!(matches!(
        cancelled,
        Err(BendingError::Plate(PlateError::Cancelled))
    ));
    let too_large = field(&mesh, [10.0, -10.0, 0.0]);
    assert!(model.evaluate(&formed, &too_large, || false).is_err());
    let mut budget = controls();
    budget.maximum_factor_bytes = 1;
    assert!(
        model
            .equilibrate(
                &formed,
                &vec![0.0; u.len()],
                &[(0, 0.0), (1, 0.0), (2, 0.0)],
                budget,
                || false
            )
            .is_err()
    );
    assert_eq!(formed, snapshot);
    let retried = model.evaluate(&virgin, &u, || false).unwrap().candidate;
    assert_eq!(formed, retried);
}
