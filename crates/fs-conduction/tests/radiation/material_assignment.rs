//! Regressions for preserving the complete solid model across radiosity iterations.
//! The two-plate factors are the parent fixture's declared algebraic model,
//! not an independently verified view factor for the finite slab geometry.
use super::*;
use fs_conduction::{ConductivityTable, ElementMaterials, MaterialId, MaterialTable};

fn assigned(mesh: &ConductionMesh, models: [ConductivityModel; 2]) -> ElementMaterials {
    let [hot, cold] = models;
    let table = MaterialTable::new([(MaterialId(11), hot), (MaterialId(29), cold)])
        .expect("material table");
    let labels = mesh.complex().tets.iter().map(|tet| {
        if mesh.positions()[tet[0] as usize][0] < 1.0 { MaterialId(11) }
        else { MaterialId(29) }
    }).collect();
    ElementMaterials::new(table, labels).expect("complete element assignments")
}

fn outer_config() -> CoupledRadiationConfig {
    CoupledRadiationConfig {
        surface_temperature_rtol: 0.0,
        surface_temperature_atol_k: 1e-8,
        relaxation: 0.25,
        max_iterations: 256,
    }
}

fn run_assigned(
    mesh: &ConductionMesh, materials: &ElementMaterials, fallback: &ConductivityModel,
) -> fs_conduction::CoupledRadiationSolution {
    let boundary = two_slab_boundary(mesh, false);
    let enclosure = two_slab_enclosure(mesh);
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| solve_with_gray_diffuse_enclosure(
        cx, ConductionProblem { mesh, boundary: &boundary, material: fallback,
            element_materials: Some(materials), source: &source },
        None, &enclosure, config(), outer_config(),
    )).expect("assigned enclosure solve")
}

/// Independent scalar balance: each slab has its own L/k, in series with
/// the gray-plate radiative resistance. No production matrix/iteration is used.
fn reference(k_hot: f64, k_cold: f64) -> (f64, f64, f64) {
    let rh = 0.5 / k_hot;
    let rc = 0.5 / k_cold;
    let mut low = 0.0;
    let mut high = 100.0 / (rh + rc);
    for _ in 0..100 {
        let q = 0.5 * low + 0.5 * high;
        let th = 400.0 - rh * q;
        let tc = 300.0 + rc * q;
        let radiation = STEFAN_BOLTZMANN_W_M2_K4
            * (th - tc) * (th + tc) * (th * th + tc * tc) / 1.5;
        if q < radiation { low = q; } else { high = q; }
    }
    let q = 0.5 * low + 0.5 * high;
    (400.0 - rh * q, 300.0 + rc * q, q)
}

#[test]
fn unequal_materials_match_both_slab_profiles_and_independent_heat_balance() {
    let mesh = two_disconnected_slabs();
    let fallback = ConductivityModel::isotropic_declared(10.0).unwrap();
    for (hot, cold) in [(4.0, 16.0), (16.0, 4.0)] {
        let materials = assigned(&mesh, [
            ConductivityModel::isotropic_declared(hot).unwrap(),
            ConductivityModel::isotropic_declared(cold).unwrap(),
        ]);
        let result = run_assigned(&mesh, &materials, &fallback);
        let (_, _, q) = reference(hot, cold);
        assert!(result.radiation.iterations > 1);
        for (point, &actual) in mesh.positions().iter().zip(&result.conduction.temperature) {
            let expected = if point[0] < 1.0 { 400.0 - q * point[0] / hot }
                else { 300.0 + q * (2.0 - point[0]) / cold };
            assert!((actual - expected).abs() < 2e-6,
                "x={} actual={actual} expected={expected}", point[0]);
        }
        let heats = &result.radiation.radiosity.net_outward_heat_w;
        assert!((heats[0] - q).abs() < 2e-5);
        assert!((heats[1] + q).abs() < 2e-5);
        assert!(result.conduction.report.energy.closure_w.abs() < 2e-5);
        // The old loop replaced both assigned laws with this fallback.
        let (wrong_hot, wrong_cold, _) = reference(10.0, 10.0);
        let (right_hot, right_cold, _) = reference(hot, cold);
        assert!((wrong_hot - right_hot).abs().max((wrong_cold - right_cold).abs()) > 1.0);
    }
}

#[test]
fn assigned_enclosure_is_independent_of_unused_fallback_even_for_nonlinear_laws() {
    let mesh = two_disconnected_slabs();
    for nonlinear in [false, true] {
        let models = if nonlinear {
            [
                ConductivityModel::isotropic(ConductivityTable::declared_curve(
                    vec![(250.0, 7.0), (500.0, 12.0)]).unwrap()),
                ConductivityModel::isotropic(ConductivityTable::declared_curve(
                    vec![(250.0, 17.0), (500.0, 22.0)]).unwrap()),
            ]
        } else {
            [
                ConductivityModel::isotropic_declared(8.0).unwrap(),
                ConductivityModel::isotropic_declared(18.0).unwrap(),
            ]
        };
        let materials = assigned(&mesh, models);
        let original = materials.clone();
        let a = run_assigned(&mesh, &materials, &ConductivityModel::isotropic_declared(6.0).unwrap());
        let b = run_assigned(&mesh, &materials, &ConductivityModel::isotropic_declared(30.0).unwrap());
        assert_eq!(a.conduction.temperature, b.conduction.temperature,
            "unused fallback changed the assigned physical problem");
        assert_eq!(a.radiation, b.radiation);
        assert_eq!(materials, original, "radiation must not mutate the material assignments");
    }
}

#[test]
fn homogeneous_assignment_retains_the_ordinary_unassigned_solution() {
    let mesh = two_disconnected_slabs();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let materials = assigned(&mesh, [material.clone(), material.clone()]);
    let assigned_result = run_assigned(&mesh, &materials, &material);
    let boundary = two_slab_boundary(&mesh, false);
    let enclosure = two_slab_enclosure(&mesh);
    let source = ScalarField::Uniform(0.0);
    let ordinary = with_cx(|cx| solve_with_gray_diffuse_enclosure(
        cx, ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source },
        None, &enclosure, config(), outer_config(),
    )).unwrap();
    assert_eq!(assigned_result.conduction.temperature, ordinary.conduction.temperature);
    assert_eq!(assigned_result.radiation, ordinary.radiation);
}

#[test]
fn published_field_resolves_under_its_original_nonlinear_materials_and_applied_fluxes() {
    let mesh = two_disconnected_slabs();
    let materials = assigned(&mesh, [
        ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 7.0), (500.0, 12.0)]).unwrap()),
        ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(250.0, 17.0), (500.0, 22.0)]).unwrap()),
    ]);
    let fallback = ConductivityModel::isotropic_declared(30.0).unwrap();
    let result = run_assigned(&mesh, &materials, &fallback);
    let enclosure = two_slab_enclosure(&mesh);
    let rows = enclosure.surfaces().iter().zip(&result.radiation.radiosity.net_outward_flux_w_m2)
        .map(|(surface, &flux)| (surface.name().to_owned(), surface.face_slots().to_vec(), flux))
        .collect::<Vec<_>>();
    let boundary = two_slab_boundary(&mesh, false)
        .with_uniform_outward_flux_overlays(&mesh, &rows).unwrap();
    let source = ScalarField::Uniform(0.0);
    let replay = with_cx(|cx| fs_conduction::solve(
        cx, ConductionProblem { mesh: &mesh, boundary: &boundary, material: &fallback,
            element_materials: Some(&materials), source: &source }, config(),
    )).expect("independent solid solve at the exact applied enclosure loads");
    assert_eq!(replay.temperature.len(), result.conduction.temperature.len());
    for (&actual, &expected) in result.conduction.temperature.iter().zip(&replay.temperature) {
        assert!((actual - expected).abs() < 2e-6,
            "published enclosure field is not a solution of its declared material problem");
    }
}
