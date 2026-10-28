//! Genuine volumetric-capacity enclosure steps share only geometry/radiation fixtures.
use super::*;
use fs_conduction::transient::{
    VolumetricHeatCapacity,
    backward_euler::{BackwardEuler, NonlinearStepConfig, StepConfig},
};
const CV: f64 = 100.0;

fn capacity_config() -> StepConfig {
    StepConfig {
        linear: LinearConfig { tolerance: 1e-12, max_iterations: 500, restart: 8 },
        energy_tolerance_j: 1e-8,
    }
}

fn nonlinear_config() -> NonlinearStepConfig {
    NonlinearStepConfig {
        residual_atol_j: 1e-10, residual_rtol: 1e-12,
        ..NonlinearStepConfig::default()
    }
}

fn independent_step(old: [f64; 2], htc: [f64; 2], imposed: [f64; 2]) -> ([f64; 2], f64) {
    let area = 3.0_f64.sqrt();
    let capacity = CV * 2.0_f64.sqrt() / 12.0;
    let diagonal = htc.map(|h| capacity + DT * area * h);
    let free: [f64; 2] = std::array::from_fn(|i|
        (capacity * old[i] + DT * area * (htc[i] * 300.0 - imposed[i])) / diagonal[i]);
    let temperature = |q: f64| [free[0] - DT * q / diagonal[0], free[1] + DT * q / diagonal[1]];
    let equilibrium = (free[0] - free[1]) / (DT / diagonal[0] + DT / diagonal[1]);
    let (mut lo, mut hi) = (equilibrium.min(0.0), equilibrium.max(0.0));
    let factor = area * STEFAN_BOLTZMANN_W_M2_K4 / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
    for _ in 0..90 {
        let q = (lo + hi) / 2.0;
        let t = temperature(q);
        let fourth_difference = (t[0] - t[1]) * (t[0] + t[1]) * (t[0] * t[0] + t[1] * t[1]);
        if q - factor * fourth_difference > 0.0 { hi = q; } else { lo = q; }
    }
    let q = (lo + hi) / 2.0;
    (temperature(q), q)
}

#[test]
fn capacity_enclosure_preserves_finite_storage_original_convection_and_prescribed_flux() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let constant = ConductivityModel::isotropic_declared(5.0).unwrap();
    let old = [440.0, 440.0, 440.0, 440.0, 300.0, 300.0, 300.0, 300.0];
    with_cx(|cx| {
        let engine = BackwardEuler::uniform(cx, &fixture.mesh, VolumetricHeatCapacity::declared(CV).unwrap()).unwrap();
        for (htc, imposed) in [([HTC; 2], [0.0; 2]), ([0.0, HTC], [-25.0, 0.0]), ([0.0; 2], [-25.0, 7.0])] {
            let law = |i: usize| if htc[i] > 0.0 {
                ThermalBc::robin(htc[i], 300.0).unwrap()
            } else { ThermalBc::neumann(imposed[i]).unwrap() };
            let boundary = ThermalBoundaryBuilder::new(&fixture.mesh)
                .region("left", |f| f.vertices.iter().all(|v| *v < 4), law(0)).unwrap()
                .region("right", |f| f.vertices.iter().all(|v| *v >= 4), law(1)).unwrap()
                .finish().unwrap();
            let (expected, heat) = independent_step([440.0, 300.0], htc, imposed);
            for nonlinear in [false, true] {
                let problem = ConductionProblem {
                    boundary: &boundary,
                    material: if nonlinear { &fixture.material } else { &constant },
                    ..fixture.problem(&source)
                };
                let got = engine.advance_prescribed_with_enclosure_radiation(cx, problem, None,
                    &old, DT, capacity_config(), nonlinear.then(nonlinear_config),
                    &fixture.enclosure, radiation_config()).unwrap();
                for (body, expected) in expected.iter().enumerate() {
                    for actual in &got.conduction.temperature[body * 4..body * 4 + 4] {
                        close(*actual, *expected, 2e-10);
                    }
                }
                let area = 3.0_f64.sqrt();
                let neumann = area * (imposed[0] + imposed[1]);
                let convection = area * (htc[0] * (expected[0] - 300.0)
                    + htc[1] * (expected[1] - 300.0));
                close(got.prescribed_neumann_out_w, neumann, 1e-12);
                close(got.convective_out_w, convection, 2e-9);
                close(got.radiation.radiosity.net_outward_heat_w[0], heat, 2e-9);
                close(got.radiation.radiosity.net_outward_heat_w[1], -heat, 2e-9);
                close(got.conduction.stored_energy_change_j, -DT * (neumann + convection), 1e-8);
                assert!(got.physical_residual_norm_j <= got.physical_residual_tolerance_j);
                assert!(got.physical_energy_residual_j.abs() <= capacity_config().energy_tolerance_j);
                assert_eq!(got.nonlinear.is_some(), nonlinear);
                assert_eq!(old, [440.0, 440.0, 440.0, 440.0, 300.0, 300.0, 300.0, 300.0]);
                let audit = engine.evaluate_prescribed_enclosure_endpoint(cx, problem, None,
                    &old, &got.conduction.temperature, DT, capacity_config(),
                    &fixture.enclosure, radiation_config()).unwrap();
                assert_eq!(audit.radiosity, got.radiation.radiosity);
                close(audit.physical.energy_residual_j, got.physical_energy_residual_j, 1e-11);
                close(audit.physical.convective_out_w, convection, 2e-9);
                if htc[0] == 0.0 && htc[1] > 0.0 {
                    assert!((got.conduction.neumann_out_w - neumann).abs() > 1.0,
                        "internal exchange must be removed from authored Neumann heat");
                }
            }
        }
    });
}

#[test]
fn capacity_enclosure_retains_full_temperature_dependent_conductivity() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let old = [350.0, 400.0, 410.0, 420.0, 280.0, 320.0, 330.0, 340.0];
    with_cx(|cx| {
        let engine = BackwardEuler::uniform(cx, &fixture.mesh, VolumetricHeatCapacity::declared(CV).unwrap()).unwrap();
        assert!(engine.advance_prescribed_with_enclosure_radiation(cx, fixture.problem(&source),
            None, &old, DT, capacity_config(), None, &fixture.enclosure, radiation_config()).is_err(),
            "the linear path must not freeze the declared k(T)");
        let got = engine.advance_prescribed_with_enclosure_radiation(cx, fixture.problem(&source),
            None, &old, DT, capacity_config(), Some(nonlinear_config()),
            &fixture.enclosure, radiation_config()).unwrap();
        let audit = engine.evaluate_prescribed_enclosure_endpoint(cx, fixture.problem(&source),
            None, &old, &got.conduction.temperature, DT, capacity_config(),
            &fixture.enclosure, radiation_config()).unwrap();
        assert!(audit.physical.residual_norm_j <= got.physical_residual_tolerance_j);
        assert!(got.radiation.solid_iterations > 0);
        let frozen = ConductivityModel::isotropic_declared(5.0).unwrap();
        let wrong = engine.advance_prescribed_with_enclosure_radiation(cx,
            ConductionProblem { material: &frozen, ..fixture.problem(&source) },
            None, &old, DT, capacity_config(), None, &fixture.enclosure, radiation_config()).unwrap();
        let rejected = engine.evaluate_prescribed_enclosure_endpoint(cx, fixture.problem(&source),
            None, &old, &wrong.conduction.temperature, DT, capacity_config(),
            &fixture.enclosure, radiation_config()).unwrap();
        assert!(rejected.physical.residual_norm_j > 1e-3,
            "rechecking a frozen-k trajectory must expose its nonlinear material defect");
        assert!(got.conduction.temperature.iter().zip(&wrong.conduction.temperature)
            .any(|(a, b)| (a - b).abs() > 1e-4));
    });
}

#[test]
fn capacity_enclosure_bath_reaction_counts_initial_prescribed_storage_once() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(10.0);
    let bath_face = |f: &fs_conduction::mesh::BoundaryFace| {
        (f.vertices.contains(&0) && f.vertices.contains(&1))
            || (f.vertices.contains(&4) && f.vertices.contains(&5))
    };
    let left = |f: &fs_conduction::mesh::BoundaryFace| f.vertices.iter().all(|v| *v < 4);
    let boundary = ThermalBoundaryBuilder::new(&fixture.mesh)
        .region("bath-left", |f| left(f) && bath_face(f), ThermalBc::dirichlet(410.0).unwrap()).unwrap()
        .region("bath-right", |f| !left(f) && bath_face(f), ThermalBc::dirichlet(320.0).unwrap()).unwrap()
        .region("left", |f| left(f) && !bath_face(f), ThermalBc::neumann(-5.0).unwrap()).unwrap()
        .region("right", |f| !left(f) && !bath_face(f), ThermalBc::robin(HTC, 300.0).unwrap()).unwrap()
        .finish().unwrap();
    let a = RadiationSurface::new(&fixture.mesh, "left", |f| left(f) && !bath_face(f),
        emissivity("left", EPS[0])).unwrap();
    let b = RadiationSurface::new(&fixture.mesh, "right", |f| !left(f) && !bath_face(f),
        emissivity("right", EPS[1])).unwrap();
    let area = a.area_m2();
    let factors = ViewFactorMatrix::infinite_parallel_plates(area).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![a, b], factors).unwrap();
    let material = ConductivityModel::isotropic_declared(5.0).unwrap();
    let problem = ConductionProblem { boundary: &boundary, material: &material, ..fixture.problem(&source) };
    let old = [350.0, 350.0, 350.0, 350.0, 300.0, 300.0, 300.0, 300.0];
    with_cx(|cx| {
        let engine = BackwardEuler::uniform(cx, &fixture.mesh, VolumetricHeatCapacity::declared(CV).unwrap()).unwrap();
        let first = engine.advance_prescribed_with_enclosure_radiation(cx, problem, None,
            &old, DT, capacity_config(), None, &enclosure, radiation_config()).unwrap();
        let volume = 2.0_f64.sqrt() / 12.0;
        let storage = CV * volume * 80.0;
        let source_w = 20.0 * volume;
        let external = area * (-5.0 + HTC * 20.0);
        close(first.conduction.stored_energy_change_j, storage, 1e-12);
        close(first.physical_dirichlet_in_w, storage / DT - source_w + external, 1e-12);
        close(first.prescribed_neumann_out_w, -5.0 * area, 1e-12);
        close(first.convective_out_w, HTC * area * 20.0, 1e-12);
        close(first.physical_residual_norm_j, 0.0, 0.0);
        assert!(first.physical_energy_residual_j.abs() <= capacity_config().energy_tolerance_j);
        assert_eq!(first.physical_residual_tolerance_j, 0.0);
        // The accepted all-prescribed endpoint has a genuine zero free
        // residual target. Reverse preparation must preserve it, not inflate
        // it to admit a state with unconverged free equations.
        let audit = engine.evaluate_prescribed_enclosure_endpoint(
            cx, problem, None, &old, &first.conduction.temperature, DT,
            capacity_config(), &enclosure, radiation_config(),
        ).unwrap();
        assert_eq!(audit.physical.free_dofs, 0);
        let response = engine.linearize_prescribed_endpoint(
            cx, ConductionProblem { boundary: &audit.physical.combined_boundary, ..problem },
            None, &old, &first.conduction.temperature, DT, capacity_config(),
            first.physical_residual_tolerance_j, &["right"],
        ).unwrap();
        let gradient = response.pullback_with_robin_feedback(
            cx, &[0.0], &[0.0], &[0.0], &[1.0; 8], 32,
        ).unwrap();
        // Every goal node is prescribed. Source, capacity, convection and
        // previous-state changes cannot change this endpoint temperature goal.
        assert_eq!(gradient.nodal_load, [0.0; 8]);
        assert_eq!(gradient.references, [0.0]);
        assert_eq!(gradient.log_htc, [0.0]);
        assert_eq!(response.previous_temperature_pullback(cx, &gradient.nodal_load).unwrap(), [0.0; 8]);
        assert_eq!(response.source_multiplier_pullback(cx, &gradient.nodal_load).unwrap(), 0.0);
        assert_eq!(response.capacity_multiplier_pullback(cx, &gradient.nodal_load).unwrap(), 0.0);
        assert_eq!(response.primal().krylov_iterations, 0);
        let next = engine.advance_prescribed_with_enclosure_radiation(cx, problem, None,
            &first.conduction.temperature, DT, capacity_config(), None, &enclosure, radiation_config()).unwrap();
        close(next.conduction.stored_energy_change_j, 0.0, 0.0);
        close(next.physical_dirichlet_in_w, -source_w + external, 1e-12);
        assert_eq!(old, [350.0, 350.0, 350.0, 350.0, 300.0, 300.0, 300.0, 300.0]);
    });
}

