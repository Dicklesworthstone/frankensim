//! G0/G1/G3/G4 mechanics checks for real P1/Q1 sheet forming and unloading.
//! Coefficients are synthetic; these tests do not calibrate metal formability.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Material, FiniteJ2Parameters};
use fs_solid::hyper2d::plastic::{
    FinitePlasticSheet, FinitePlasticSheetError, SheetAssemblyBudget,
};
use fs_solid::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use fs_solid::mesh2::Mesh2;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

const MU: f64 = 80.0e6;
const KAPPA: f64 = 160.0e6;
const YIELD: f64 = 2.0e6;
const HARDENING: f64 = 8.0e6;
const E: f64 = 9.0 * KAPPA * MU / (3.0 * KAPPA + MU);
const NU: f64 = (3.0 * KAPPA - 2.0 * MU) / (2.0 * (3.0 * KAPPA + MU));
const THICKNESS: f64 = 0.001;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            StreamKey { seed: 37, kernel_id: 5, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn law(yield_stress_pa: f64) -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters { shear_modulus_pa: MU, bulk_modulus_pa: KAPPA,
            yield_stress_pa, hardening_modulus_pa: HARDENING },
        FiniteJ2Limits { minimum_jacobian: 0.2, maximum_elastic_log_strain: 0.6,
            maximum_plastic_log_strain: 2.0, maximum_accumulated_plastic_strain: 3.0 },
    ).unwrap()
}
fn budget() -> SheetAssemblyBudget {
    SheetAssemblyBudget { maximum_nodes: 100, maximum_elements: 100,
        maximum_free_dofs: 200, minimum_scaled_jacobian: 1e-4 }
}
fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls { equilibrium_work_tolerance_j: 2e-5,
        maximum_endpoint_work_defect_j: 1000.0, yield_tolerance_pa: 2e-4 }
}
fn close(value: f64, expected: f64, tolerance: f64) {
    assert!((value - expected).abs() <= tolerance,
        "actual {value:e}, expected {expected:e}, tolerance {tolerance:e}");
}

// Only left x and bottom y are held. Right x is prescribed during forming,
// then released; all lateral contraction and the interior motion are solved.
fn rollers(mesh: &Mesh2, stretch: Option<f64>) -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
        if x == 0.0 { prescribed.push((2 * node, 0.0)); }
        if y == 0.0 { prescribed.push((2 * node + 1, 0.0)); }
        if x == 1.0 && let Some(stretch) = stretch { prescribed.push((2 * node, stretch - 1.0)); }
    }
    prescribed
}

#[test]
fn tension_and_compression_solve_free_lateral_motion_and_permanent_shape() {
    with_cx(|cx, _| {
        for triangular in [false, true] {
            let mesh = if triangular { Mesh2::triangles(1.0, 1.0, 2, 2) } else { Mesh2::quads(1.0, 1.0, 2, 2) };
            let thickness = vec![THICKNESS; mesh.elems.len()];
            let material = law(YIELD);
            for sign in [-1.0, 1.0] {
                let mut sheet = FinitePlasticSheet::try_new(&mesh, FinitePlasticLawField::Uniform(&material),
                    &thickness, PlaneStressSettings::default(), budget(), cx).unwrap();
                let mut total_heat = 0.0;
                for step in 1..=4 {
                    let stretch = 1.0 + sign * 0.05 * f64::from(step);
                    let endpoint = sheet.advance(&rollers(&mesh, Some(stretch)), &[],
                        IncrementSettings::default(), controls(), cx).unwrap();
                    assert!(endpoint.equilibrium.residual_history.len() > 1);
                    assert!(endpoint.energy.plastic_dissipation_j > 0.0);
                    close(endpoint.energy.equilibrium_work_residual_j, 0.0, controls().equilibrium_work_tolerance_j);
                    total_heat += endpoint.energy.plastic_dissipation_j;
                    sheet = endpoint.sheet;
                }
                // Independent closed form for coaxial Hencky/J2 uniaxial stress.
                let stretch: f64 = 1.0 + sign * 0.2;
                let alpha = (stretch.ln().abs() - YIELD / E) / (1.0 + HARDENING / E);
                let tau = sign * (YIELD + HARDENING * alpha);
                let lateral = (-NU * tau / E - sign * alpha / 2.0).exp();
                for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                    close(sheet.displacement_m()[node][0], (stretch - 1.0) * x, 2e-8);
                    close(sheet.displacement_m()[node][1], (lateral - 1.0) * y, 2e-8);
                }
                for ((state, f), p) in sheet.quadrature_states().iter().zip(sheet.quadrature_deformation()).zip(sheet.quadrature_first_piola_pa()) {
                    close(state.alpha, alpha, 3e-8);
                    close(f[8], lateral, 3e-8);
                    close(p[0], tau / stretch, 0.5);
                    close(p[4], 0.0, 0.5);
                    close(p[8], 0.0, 0.01);
                }
                close(total_heat, THICKNESS * YIELD * alpha, 2e-4);
                let formed = sheet.quadrature_states().to_vec();
                let released = sheet.advance(&rollers(&mesh, None), &[],
                    IncrementSettings::default(), controls(), cx).unwrap();
                let residual_x = (sign * alpha).exp();
                let residual_y = (-sign * alpha / 2.0).exp();
                for (node, &[x, y]) in mesh.nodes.iter().enumerate() {
                    close(released.sheet.displacement_m()[node][0], (residual_x - 1.0) * x, 3e-8);
                    close(released.sheet.displacement_m()[node][1], (residual_y - 1.0) * y, 3e-8);
                }
                for (old, new) in formed.iter().zip(released.sheet.quadrature_states()) {
                    close(new.alpha, old.alpha, 2e-10);
                }
                for thickness in released.sheet.quadrature_thickness_m() { close(thickness, THICKNESS * residual_y, 3e-11); }
                close(released.energy.plastic_dissipation_j, 0.0, 1e-6);
                // Remaining energy is the energetic hardening storage, not heat.
                close(released.sheet.stored_energy_j(), THICKNESS * HARDENING * alpha * alpha / 2.0, 2e-4);
                assert!((residual_x - 1.0).abs() > 0.1);
            }
        }
    });
}

#[test]
fn native_sheet_tangent_matches_force_differences_at_frozen_formed_history() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = law(YIELD);
        let thickness = [THICKNESS];
        let initial = FinitePlasticSheet::try_new(&mesh, FinitePlasticLawField::Uniform(&material),
            &thickness, PlaneStressSettings::default(), budget(), cx).unwrap();
        // Eight DOFs also exercises the shared solver's actual 2-D indexing.
        let sheet = initial.advance(&rollers(&mesh, Some(1.08)), &[], IncrementSettings::default(), controls(), cx).unwrap().sheet;
        let states = sheet.quadrature_states().to_vec();
        let trial: Vec<_> = mesh.nodes.iter().zip(sheet.displacement_m()).map(|(&[x, y], &u)| [u[0] + 0.01*x + 0.004*y, u[1] - 0.003*y]).collect();
        let direction: Vec<_> = mesh.nodes.iter().map(|&[x, y]| [0.02*x + 0.007*y, -0.013*y + 0.005*x]).collect();
        let result = sheet.evaluate_trial(&trial, true, cx).unwrap();
        let tangent = result.tangent_n_m.unwrap();
        let mut action = vec![0.0; 2*mesh.nodes.len()];
        tangent.spmv(direction.as_flattened(), &mut action);
        let h = 1e-5;
        let perturbed = |sign| trial.iter().zip(&direction).map(|(u, d)| [u[0] + sign*h*d[0], u[1] + sign*h*d[1]]).collect::<Vec<_>>();
        let plus = sheet.evaluate_trial(&perturbed(1.0), false, cx).unwrap();
        let minus = sheet.evaluate_trial(&perturbed(-1.0), false, cx).unwrap();
        let mut error = 0.0_f64;
        let mut norm = 0.0_f64;
        for (i, (&plus, &minus)) in plus.internal_force_n.as_flattened().iter().zip(minus.internal_force_n.as_flattened()).enumerate() {
            let numerical = (plus - minus) / (2.0*h);
            error = error.hypot(numerical - action[i]);
            norm = norm.hypot(action[i]);
        }
        assert!(error / norm < 3e-6, "relative tangent error {}", error / norm);
        let energy_derivative = (plus.incremental_potential_j - minus.incremental_potential_j) / (2.0*h);
        let force_work: f64 = result.internal_force_n.as_flattened().iter().zip(direction.as_flattened()).map(|(f,d)| f*d).sum();
        close(energy_derivative, force_work, 2e-5);
        let translated: Vec<_> = trial.iter().map(|u| [u[0] + 10.0, u[1] - 20.0]).collect();
        let shifted = sheet.evaluate_trial(&translated, false, cx).unwrap();
        for (a,b) in result.internal_force_n.as_flattened().iter().zip(shifted.internal_force_n.as_flattened()) { close(*a,*b, 1e-7); }
        assert_eq!(sheet.quadrature_states(), states);
    });
}

#[derive(Debug)]
struct BoundMaterial { name: &'static str, law: FiniteJ2 }
impl FiniteJ2Material for BoundMaterial { fn finite_j2(&self) -> &FiniteJ2 { &self.law } }

#[test]
fn heterogeneous_original_laws_and_reference_thickness_survive_forming() {
    with_cx(|cx, _| {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 1);
        let laws = [BoundMaterial { name: "soft", law: law(YIELD) }, BoundMaterial { name: "strong", law: law(5.0*YIELD) }];
        let thickness = [THICKNESS, 2.0*THICKNESS];
        let sheet = FinitePlasticSheet::try_new(&mesh, FinitePlasticLawField::PerElement(&laws),
            &thickness, PlaneStressSettings::default(), budget(), cx).unwrap();
        close(sheet.reference_volume_m3(0).unwrap(), 0.5*THICKNESS, 1e-16);
        close(sheet.reference_volume_m3(1).unwrap(), THICKNESS, 1e-16);
        let candidate = sheet.advance(&rollers(&mesh, Some(1.05)), &[], IncrementSettings::default(), controls(), cx).unwrap();
        let FinitePlasticLawField::PerElement(retained) = candidate.sheet.laws() else { panic!("original per-element binding lost"); };
        assert!(std::ptr::eq(&retained[0], &laws[0]));
        assert!(std::ptr::eq(&retained[1], &laws[1]));
        assert_eq!(retained[1].name, "strong");
        let states = candidate.sheet.quadrature_states();
        assert_eq!(states.len(), 8);
        assert!(states[..4].iter().map(|s| s.alpha).sum::<f64>() > 0.05);
        assert!(states[4..].iter().all(|s| s.alpha < 1e-9));
        assert!(candidate.energy.plastic_dissipation_j > 0.0);
        assert_eq!(sheet.displacement_m(), &[[0.0;2];6]);
    });
}

#[test]
fn refusals_and_cancellation_leave_shape_and_plastic_history_unpublished() {
    with_cx(|cx, gate| {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = law(YIELD);
        let thickness = [THICKNESS];
        let sheet = FinitePlasticSheet::try_new(&mesh, FinitePlasticLawField::Uniform(&material),
            &thickness, PlaneStressSettings::default(), budget(), cx).unwrap();
        let states = sheet.quadrature_states().to_vec();
        let deformation = sheet.quadrature_deformation().to_vec();
        let mut strict = controls();
        strict.maximum_endpoint_work_defect_j = 1e-12;
        assert!(matches!(sheet.advance(&rollers(&mesh, Some(1.05)), &[], IncrementSettings::default(), strict, cx), Err(FinitePlasticSheetError::Energy { .. })));
        assert!(sheet.advance(&[(99,0.0)], &[], IncrementSettings::default(), controls(), cx).is_err());
        assert!(sheet.advance(&[(0,0.0),(0,1.0)], &[], IncrementSettings::default(), controls(), cx).is_err());
        assert!(sheet.evaluate_trial(&[[0.0,0.0]], false, cx).is_err());
        let inverted: Vec<_> = mesh.nodes.iter().map(|&[x,_]| [-2.0*x,0.0]).collect();
        assert!(sheet.evaluate_trial(&inverted, true, cx).is_err());
        assert_eq!(sheet.quadrature_states(), states);
        assert_eq!(sheet.quadrature_deformation(), deformation);
        let mut folded = mesh.clone();
        folded.elems[0].swap(1,2);
        assert!(FinitePlasticSheet::try_new(&folded, FinitePlasticLawField::Uniform(&material), &thickness,
            PlaneStressSettings::default(), budget(), cx).is_err());
        gate.request();
        assert_eq!(sheet.advance(&rollers(&mesh, Some(1.05)), &[], IncrementSettings::default(), controls(), cx).unwrap_err(),
            FinitePlasticSheetError::Solid(SmallStrainTetError::Cancelled));
        assert_eq!(sheet.quadrature_states(), states);
        assert_eq!(sheet.displacement_m(), &[[0.0;2];4]);
    });
}
