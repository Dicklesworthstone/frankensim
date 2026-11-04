//! Reflected patch exchange must use the same live contact endpoint.
use super::*;
use fs_conduction::radiation::{CoupledRadiationConfig, GrayDiffuseEnclosure,
    OpenDiffuseEnclosure, RadiationSurface, ViewFactorEvidence, ViewFactorMatrix,
    ViewFactorTolerance, solve_with_convective_gray_diffuse_enclosure};
use fs_conduction::DofMap;

fn boundary(mesh: &ConductionMesh, prescribed: bool) -> ThermalBoundary {
    let builder = ThermalBoundaryBuilder::new(mesh)
        .region("cold", |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0,
            if prescribed { ThermalBc::neumann(0.0).unwrap() }
            else { ThermalBc::robin(2.0, 297.0).unwrap() }).unwrap()
        .region("hot", |f| f.centroid[2] > 0.0 && f.centroid[1] == 0.0,
            ThermalBc::neumann(0.0).unwrap()).unwrap();
    let builder = if prescribed {
        builder.region("reservoir", |f| f.centroid[2] < 0.0 && f.centroid[0] == 0.0,
            ThermalBc::dirichlet(300.0).unwrap()).unwrap()
    } else { builder };
    builder.adiabatic_remainder().finish().unwrap()
}

fn enclosure(mesh: &ConductionMesh) -> GrayDiffuseEnclosure {
    let emissivity = patch().emissivity().clone();
    let cold = RadiationSurface::new(mesh, "cold",
        |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0, emissivity.clone()).unwrap();
    let hot = RadiationSurface::new(mesh, "hot",
        |f| f.centroid[2] > 0.0 && f.centroid[1] == 0.0, emissivity).unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(cold.area_m2()).unwrap();
    GrayDiffuseEnclosure::new(vec![cold, hot], factors).unwrap()
}

fn controls() -> CoupledRadiationConfig {
    CoupledRadiationConfig { max_iterations: 128, surface_temperature_atol_k: 1e-10,
        surface_temperature_rtol: 0.0, relaxation: 0.5 }
}

fn open_enclosure(cx: &Cx<'_>, mesh: &ConductionMesh) -> OpenDiffuseEnclosure {
    OpenDiffuseEnclosure::new(cx, enclosure(mesh).surfaces().to_vec(),
        vec![vec![0.0, 0.25], vec![0.25, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "synthetic reciprocal open exchange".into() },
        ViewFactorTolerance::default()).unwrap()
}

fn surface_card() -> InterfaceSystemCard {
    let original = card(false);
    let mut claim = original.claims().claims_ordered().next().unwrap().1.clone();
    claim.value = PropertyValue::Surface2D {
        axes: ["T".into(), "normal_pressure".into()],
        axis_dims: [fs_conduction::TEMPERATURE_DIMS, fs_qty::Dims([-1, 1, -2, 0, 0, 0])],
        samples: vec![([280.0, 999_999.0], 0.09), ([380.0, 999_999.0], 0.69),
            ([280.0, 1_000_001.0], 0.11), ([380.0, 1_000_001.0], 0.71)],
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    };
    let mut claims = ClaimSet::new(); claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(original.surface_a().clone(), original.surface_b().clone(),
        original.context().clone(), claims, vec![]).unwrap()
}

#[test]
fn enclosure_newton_and_source_surface_picard_share_fresh_physics_and_count_every_response() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let bc = boundary(&mesh, false); let model = enclosure(&mesh);
        let card = card(false); let surface = surface_card();
        let curve = contacts(cx, &mesh, &bc, &card);
        let surface = contacts(cx, &mesh, &bc, &surface);
        let material = material(true); let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
            source: &source, element_materials: None };
        let initial = [300.0; 8];
        let newton = curve.solve_steady_newton_with_gray_diffuse_enclosure(cx, problem,
            &initial, config(), &model, controls()).unwrap();
        let mut solver = SolveConfig::default();
        solver.stop.step_atol = 0.0; solver.stop.residual_rtol = 1e-12;
        solver.stop.residual_atol = 1e-10; solver.linear.tolerance = 1e-12;
        let picard = surface.solve_steady_with_gray_diffuse_enclosure(cx, problem, &initial,
            solver.clone(), SteadyContactConfig { max_iterations: 64, relaxation: 0.5,
                residual_tolerance_w: 1e-8, energy_tolerance_w: 1e-8 }, &model, controls()).unwrap();
        for (&a, &b) in newton.contact.temperature.iter().zip(&picard.contact.conduction.temperature) {
            near(a, b, 2e-7);
        }
        assert_eq!(newton.work.thermal_responses, newton.radiation.iterations + 1);
        assert!(picard.work.thermal_responses > picard.radiation.iterations + 1);
        assert!(picard.work.nonlinear_iterations >= picard.contact.nonlinear_iterations);
        assert!(picard.work.krylov_iterations >= picard.contact.krylov_iterations);
        assert!(newton.work.residual_evaluations >= newton.contact.residual_evaluations);
        near(newton.convective_out_w, 5.0, 1e-8);
        near(picard.convective_out_w, 5.0, 1e-8);
        let actual = &newton.radiation.equilibrium.nonlinear_radiosity;
        near(actual.enclosure_energy_closure_w, 0.0, 1e-12);
        assert!(actual.net_outward_heat_w[0].abs() > 0.1);
        let state = &newton.contact.physical_contacts.states()[0];
        near(state.resistance.value_m2_k_per_w(),
            0.1 + 0.006 * (state.interface_temperature_k - 280.0), 1e-14);
        let frozen = curve.evaluate(cx, &bc, &initial).unwrap();
        let frozen = solve_with_convective_gray_diffuse_enclosure(cx, problem,
            Some(frozen.frozen_interfaces()), &model, solver, controls()).unwrap();
        assert!(newton.contact.temperature.iter().zip(&frozen.conduction.temperature)
            .any(|(a, b)| (a - b).abs() > 1e-3), "a frozen source is a different physical problem");
        assert_eq!(initial, [300.0; 8]);
    });
}

#[test]
fn enclosure_direct_loads_reach_free_rows_and_prescribed_reservoir_reactions() {
    context(&CancelGate::new_clock_free(), |cx| {
        for prescribed in [false, true] {
            let mesh = mesh(); let bc = boundary(&mesh, prescribed); let model = enclosure(&mesh);
            let card = card(false); let contacts = contacts(cx, &mesh, &bc, &card);
            let material = material(true); let source = source();
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
                source: &source, element_materials: None };
            let accepted = contacts.solve_steady_newton_with_gray_diffuse_enclosure(cx, problem,
                &[300.0; 8], config(), &model, controls()).unwrap();
            let t = &accepted.contact.temperature;
            let fresh = contacts.evaluate(cx, &bc, t).unwrap();
            assert_eq!(fresh.states(), accepted.contact.physical_contacts.states());
            let system = fs_conduction::assemble::assemble_operator_with_interfaces(cx, &mesh, &bc,
                &material, &source, t, fresh.frozen_interfaces()).unwrap();
            let mut residual = fs_conduction::assemble::full_residual(&system, t);
            let dofs = DofMap::new(&bc, 8).unwrap();
            let old_reaction: f64 = dofs.fixed().iter().map(|&v| residual[v]).sum();
            let radiosity = &accepted.radiation.equilibrium.nonlinear_radiosity;
            for (patch, &flux) in model.surfaces().iter().zip(&radiosity.net_outward_flux_w_m2) {
                for &slot in patch.face_slots() {
                    let face = &mesh.boundary()[slot];
                    for &vertex in &face.vertices { residual[vertex as usize] += face.area * (flux / 3.0); }
                }
            }
            let norm = dofs.free().iter().map(|&v| residual[v].powi(2)).sum::<f64>().sqrt();
            near(norm, accepted.contact.residual_norm_w, 1e-11);
            assert!(norm <= config().residual_tolerance_w);
            for &(vertex, value) in &accepted.contact.dirichlet_reactions_w {
                near(value, residual[vertex], 1e-12);
            }
            let reaction: f64 = dofs.fixed().iter().map(|&v| residual[v]).sum();
            near(accepted.contact.energy.dirichlet_in_w, reaction, 1e-12);
            near(5.0 + reaction - accepted.convective_out_w, 0.0, 1e-8);
            if prescribed {
                assert!((old_reaction - reaction).abs() > 0.1,
                    "radiation on fixed trace vertices changes the actual imposed reservoir heat");
            }
        }
    });
}

#[test]
fn closed_heat_cancellation_and_a_loose_temperature_gate_cannot_accept_one_loaded_response() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let bc = boundary(&mesh, false); let model = enclosure(&mesh);
        let card = card(false); let contacts = contacts(cx, &mesh, &bc, &card);
        let material = material(true); let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
            source: &source, element_materials: None };
        let initial = [300.0; 8];
        let refused = contacts.solve_steady_newton_with_gray_diffuse_enclosure(cx, problem,
            &initial, config(), &model, CoupledRadiationConfig { max_iterations: 1,
                surface_temperature_atol_k: 1e6, ..controls() });
        assert!(matches!(refused, Err(ConductionError::Radiation { .. })));
        assert_eq!(initial, [300.0; 8]);
    });
}

#[test]
fn open_contact_newton_and_surface_picard_keep_the_black_bath_out_of_solid_heat() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let bc = boundary(&mesh, false); let model = open_enclosure(cx, &mesh);
        let card = card(false); let surface = surface_card();
        let curve = contacts(cx, &mesh, &bc, &card);
        let surface = contacts(cx, &mesh, &bc, &surface);
        let material = material(true); let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
            source: &source, element_materials: None };
        let initial = [300.0; 8];
        for bath in [290.0_f64, 320.0] {
            let newton = curve.solve_steady_newton_with_open_gray_diffuse_enclosure(cx, problem,
                &initial, config(), &model, bath, controls()).unwrap();
            let mut solver = SolveConfig::default();
            solver.stop.step_atol = 0.0; solver.stop.residual_rtol = 1e-12;
            solver.stop.residual_atol = 1e-10; solver.linear.tolerance = 1e-12;
            let picard = surface.solve_steady_with_open_gray_diffuse_enclosure(cx, problem, &initial,
                solver, SteadyContactConfig { max_iterations: 64, relaxation: 0.5,
                    residual_tolerance_w: 1e-8, energy_tolerance_w: 1e-8 },
                &model, bath, controls()).unwrap();
            for (&a, &b) in newton.contact.temperature.iter().zip(&picard.contact.conduction.temperature) {
                near(a, b, 2e-7);
            }
            assert_eq!(newton.work.thermal_responses, newton.radiation.iterations + 1);
            assert!(picard.work.thermal_responses > picard.radiation.iterations + 1);
            for (contact, convection, radiation) in [(&newton.contact.energy, newton.convective_out_w,
                &newton.radiation), (&picard.contact.energy, picard.convective_out_w, &picard.radiation)]
            {
                let actual = &radiation.equilibrium.nonlinear_radiosity;
                assert_eq!(actual.net_outward_heat_w.len(), 3, "two physical traces and one bath");
                // Independently eliminate two reflected radiosities relative to
                // the bath: dJ_i - (1-epsilon) F_ij dJ_j = epsilon sigma (Ti^4-Tb^4).
                let a = 0.8 * 0.25;
                let b: Vec<_> = actual.surface_temperatures_k[..2].iter()
                    .map(|t| 0.2 * STEFAN_BOLTZMANN_W_M2_K4 * (t.powi(4) - bath.powi(4)))
                    .collect();
                let radiosity = [(b[0] + a*b[1]) / (1.0-a*a), (b[1] + a*b[0]) / (1.0-a*a)];
                for i in 0..2 {
                    near(actual.net_outward_heat_w[i], 0.5*(radiosity[i] - 0.25*radiosity[1-i]), 2e-12);
                }
                let outward = actual.net_outward_heat_w[..2].iter().sum::<f64>();
                assert!(outward.abs() > 0.1);
                assert_eq!(outward.is_sign_positive(), bath == 290.0);
                near(outward, -actual.net_outward_heat_w[2], 1e-12);
                near(convection + outward, 5.0, 1e-8);
                near(contact.robin_out_w, convection + outward, 1e-12);
                near(contact.closure_w, 0.0, 1e-8);
            }
            let fresh = curve.evaluate(cx, &bc, &newton.contact.temperature).unwrap();
            assert_eq!(fresh.states(), newton.contact.physical_contacts.states());
        }
        assert_eq!(initial, [300.0; 8]);
    });
}

#[test]
fn open_contact_radiation_enters_prescribed_reactions_and_actual_free_residual() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let bc = boundary(&mesh, true); let model = open_enclosure(cx, &mesh);
        let card = card(false); let contacts = contacts(cx, &mesh, &bc, &card);
        let material = material(true); let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
            source: &source, element_materials: None };
        let accepted = contacts.solve_steady_newton_with_open_gray_diffuse_enclosure(cx, problem,
            &[300.0; 8], config(), &model, 290.0, controls()).unwrap();
        let t = &accepted.contact.temperature;
        let fresh = contacts.evaluate(cx, &bc, t).unwrap();
        let system = fs_conduction::assemble::assemble_operator_with_interfaces(cx, &mesh, &bc,
            &material, &source, t, fresh.frozen_interfaces()).unwrap();
        let mut residual = fs_conduction::assemble::full_residual(&system, t);
        let dofs = DofMap::new(&bc, 8).unwrap();
        let old_reaction: f64 = dofs.fixed().iter().map(|&v| residual[v]).sum();
        let actual = &accepted.radiation.equilibrium.nonlinear_radiosity;
        for (patch, &flux) in model.surfaces().iter().zip(&actual.net_outward_flux_w_m2) {
            for &slot in patch.face_slots() {
                let face = &mesh.boundary()[slot];
                for &vertex in &face.vertices { residual[vertex as usize] += face.area * (flux / 3.0); }
            }
        }
        let norm = dofs.free().iter().map(|&v| residual[v].powi(2)).sum::<f64>().sqrt();
        near(norm, accepted.contact.residual_norm_w, 1e-11);
        assert!(norm <= config().residual_tolerance_w);
        let reaction: f64 = dofs.fixed().iter().map(|&v| residual[v]).sum();
        for &(vertex, value) in &accepted.contact.dirichlet_reactions_w {
            near(value, residual[vertex], 1e-12);
        }
        assert!((old_reaction - reaction).abs() > 0.1);
        near(accepted.contact.energy.dirichlet_in_w, reaction, 1e-12);
        near(5.0 + reaction - accepted.convective_out_w + actual.net_outward_heat_w[2], 0.0, 1e-8);
    });
}

#[test]
fn open_contact_enclosure_retains_domain_work_and_cancellation_gates() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let bc = boundary(&mesh, false); let model = open_enclosure(cx, &mesh);
        let card = card(false); let contacts = contacts(cx, &mesh, &bc, &card);
        let material = material(true); let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &bc, material: &material,
            source: &source, element_materials: None };
        let initial = [300.0; 8];
        for bath in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(contacts.solve_steady_newton_with_open_gray_diffuse_enclosure(cx, problem,
                &initial, config(), &model, bath, controls()).is_err());
            assert!(contacts.solve_steady_with_open_gray_diffuse_enclosure(cx, problem,
                &initial, SolveConfig::default(), SteadyContactConfig { max_iterations: 64,
                    relaxation: 0.5, residual_tolerance_w: 1e-8, energy_tolerance_w: 1e-8 },
                &model, bath, controls()).is_err());
        }
        let refused = contacts.solve_steady_newton_with_open_gray_diffuse_enclosure(cx, problem,
            &initial, config(), &model, 290.0, CoupledRadiationConfig {
                max_iterations: 1, surface_temperature_atol_k: 1e6, ..controls() });
        assert!(matches!(refused, Err(ConductionError::Radiation { .. })));
        let gate = CancelGate::new_clock_free(); gate.request();
        context(&gate, |cx| assert!(matches!(
            contacts.solve_steady_newton_with_open_gray_diffuse_enclosure(cx, problem,
                &initial, config(), &model, 290.0, controls()), Err(ConductionError::Cancelled { .. }))));
        assert_eq!(initial, [300.0; 8]);
    });
}
