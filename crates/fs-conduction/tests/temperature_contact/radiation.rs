//! Fresh contact and exact radiative heat must balance on the returned field.

use super::*;
#[path = "radiation/enclosure.rs"]
mod enclosure;
use fs_conduction::interface::temperature::steady::{SteadyContactConfig, SteadyContactNewtonConfig};
use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, ConductionSolver, EMISSIVITY_DIMS, InitialGuess,
    STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ThermalBc,
};

fn patch() -> AmbientRadiationPatch {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.2,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 380.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic steady contact/radiation fixture".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = fs_matdb::MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray surface".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    AmbientRadiationPatch::new(
        "exterior",
        SurfaceEmissivity::from_card("exterior", &card, 300.0, SelectionPolicy::SingleClaimOnly)
            .unwrap(),
        300.0,
    )
    .unwrap()
}

fn config() -> SteadyContactNewtonConfig {
    let mut config = SteadyContactNewtonConfig::default();
    config.newton.absolute_tolerance = 1.0e-10;
    config.newton.forcing_minimum = 1.0e-12;
    config.newton.forcing_maximum = 1.0e-10;
    config.newton.linear_restart = 8;
    config.max_face_pairs = 1;
    config
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
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(original.surface_a().clone(), original.surface_b().clone(),
        original.context().clone(), claims, vec![]).unwrap()
}

#[test]
fn steady_radiation_and_live_contact_close_the_same_nodal_and_external_heat_balance() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let source = source();
        let material = material(true);
        let card = card(false);
        for convection in [false, true] {
            let condition = if convection {
                ThermalBc::robin(2.0, 297.0).unwrap()
            } else {
                ThermalBc::neumann(0.0).unwrap()
            };
            let bc = ThermalBoundaryBuilder::new(&mesh)
                .region(
                    "exterior",
                    |face| face.centroid[2] < 0.0 && face.centroid[1] == 0.0,
                    condition,
                )
                .unwrap()
                .adiabatic_remainder()
                .finish()
                .unwrap();
            let contacts = contacts(cx, &mesh, &bc, &card);
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &bc,
                material: &material,
                element_materials: None,
                source: &source,
            };
            let solver = SolveConfig {
                initial: InitialGuess::Uniform(300.0),
                ..Default::default()
            };
            let initial = ConductionSolver::initial_temperature(problem, &solver).unwrap();
            if !convection {
                assert!(matches!(
                    ConductionSolver::new(problem, solver),
                    Err(ConductionError::SingularPureNeumann)
                ));
            }
            let before = initial.clone();
            let solved = contacts
                .solve_steady_newton_with_ambient_radiation(
                    cx,
                    problem,
                    &initial,
                    config(),
                    &[patch()],
                    AmbientRadiationConfig::default(),
                )
                .unwrap();
            assert_eq!(
                initial, before,
                "radiation trials only change their own warm start"
            );
            assert!(solved.radiation.iterations > 1);
            assert!(solved.radiation.solid_iterations >= solved.contact.newton.iterations);
            assert!(solved.residual_evaluations >= solved.contact.residual_evaluations);
            assert!(solved.jacobian_assemblies >= solved.contact.jacobian_assemblies);
            let temperature = &solved.contact.temperature;
            let mean = (temperature[0] + temperature[1] + temperature[3]) / 3.0;
            let radiation =
                0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5 * (mean.powi(4) - 300.0_f64.powi(4));
            let convection_heat = if convection {
                2.0 * 0.5 * (mean - 297.0)
            } else {
                0.0
            };
            near(solved.radiation.nonlinear_radiation_out_w, radiation, 1e-10);
            near(solved.convective_out_w, convection_heat, 1e-10);
            near(
                radiation + convection_heat,
                5.0,
                config().energy_tolerance_w,
            );
            assert_eq!(solved.convective_robin_fluxes.is_empty(), !convection);
            let state = &solved.contact.physical_contacts.states()[0];
            near(
                state.resistance.value_m2_k_per_w(),
                0.1 + 0.006 * (state.interface_temperature_k - 280.0),
                1e-14,
            );
            assert!((state.interface_temperature_k - 300.0).abs() > 0.1);
            let current = contacts
                .evaluate(cx, &solved.physical_boundary, temperature)
                .unwrap();
            assert_eq!(current.states(), solved.contact.physical_contacts.states());
            let system = fs_conduction::assemble::assemble_operator_with_interfaces(
                cx,
                &mesh,
                &solved.physical_boundary,
                &material,
                &source,
                temperature,
                current.frozen_interfaces(),
            )
            .unwrap();
            let residual = fs_conduction::assemble::full_residual(&system, temperature);
            let norm = residual
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
            near(norm, solved.contact.residual_norm_w, 1e-11);
            assert!(norm <= config().residual_tolerance_w);
            assert!(solved.contact.energy.closure_w.abs() <= config().energy_tolerance_w);
            let heat: f64 = current
                .frozen_interfaces()
                .fluxes(temperature)
                .unwrap()
                .iter()
                .map(|row| row.heat_rate_a_to_b_w)
                .sum();
            near(heat.abs(), 5.0, config().energy_tolerance_w);
        }
    });
}

#[test]
fn steady_radiative_contact_refuses_incomplete_outer_work_without_changing_the_input() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "exterior",
                |face| face.centroid[2] < 0.0 && face.centroid[1] == 0.0,
                ThermalBc::neumann(0.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let card = card(false);
        let contacts = contacts(cx, &mesh, &bc, &card);
        let source = source();
        let material = material(true);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &bc,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let initial = [300.0; 8];
        let refused = contacts.solve_steady_newton_with_ambient_radiation(
            cx,
            problem,
            &initial,
            config(),
            &[patch()],
            AmbientRadiationConfig {
                max_iterations: 1,
                ..Default::default()
            },
        );
        assert!(matches!(
            refused,
            Err(ConductionError::AmbientRadiationNotConverged { .. })
        ));
        assert_eq!(initial, [300.0; 8]);
        assert!(
            contacts
                .solve_steady_newton_with_ambient_radiation(
                    cx,
                    problem,
                    &initial,
                    config(),
                    &[],
                    AmbientRadiationConfig::default(),
                )
                .is_err()
        );
    });
}

#[test]
fn steady_source_surfaces_share_the_radiative_endpoint_audit_without_a_source_tangent() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let original = card(false);
        let surface = surface_card();
        let material = material(true);
        let source = source();
        let mut solver = SolveConfig::default();
        solver.stop.step_atol = 0.0;
        solver.stop.residual_rtol = 1e-12;
        solver.stop.residual_atol = 1e-11;
        solver.linear.tolerance = 1e-12;
        let coupling = SteadyContactConfig { max_iterations: 64, relaxation: 0.5,
            residual_tolerance_w: 1e-8, energy_tolerance_w: 1e-8 };
        for convection in [false, true] {
            let bc = ThermalBoundaryBuilder::new(&mesh)
                .region("exterior", |face| face.centroid[2] < 0.0 && face.centroid[1] == 0.0,
                    if convection { ThermalBc::robin(2.0, 297.0).unwrap() }
                    else { ThermalBc::neumann(0.0).unwrap() }).unwrap()
                .adiabatic_remainder().finish().unwrap();
            let contacts = contacts(cx, &mesh, &bc, &surface);
            let reference = super::contacts(cx, &mesh, &bc, &original);
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc,
                material: &material, element_materials: None, source: &source };
            let initial = [300.0; 8];
            let accepted = contacts.solve_steady_with_ambient_radiation(cx, problem, &initial,
                solver.clone(), coupling, &[patch()], AmbientRadiationConfig::default()).unwrap();
            let curve = reference.solve_steady_newton_with_ambient_radiation(cx, problem, &initial,
                config(), &[patch()], AmbientRadiationConfig::default()).unwrap();
            let temperature = &accepted.contact.conduction.temperature;
            for (&a, &b) in temperature.iter().zip(&curve.contact.temperature) { near(a, b, 2e-7); }
            assert!(accepted.thermal_responses > accepted.radiation.iterations);
            assert!(accepted.thermal_responses >= accepted.contact.iterations);
            assert!(accepted.radiation.solid_iterations >= accepted.contact.nonlinear_iterations);
            assert!(accepted.radiation.krylov_iterations >= accepted.contact.krylov_iterations);
            assert_eq!(accepted.convective_robin_fluxes.is_empty(), !convection);
            near(accepted.radiation.nonlinear_radiation_out_w + accepted.convective_out_w,
                5.0, coupling.energy_tolerance_w);
            let current = contacts.evaluate(cx, &accepted.physical_boundary, temperature).unwrap();
            assert_eq!(current.states(), accepted.contact.physical_contacts.states());
            let state = &current.states()[0];
            near(state.resistance.value_m2_k_per_w(),
                0.1 + 0.006 * (state.interface_temperature_k - 280.0), 1e-14);
            let system = fs_conduction::assemble::assemble_operator_with_interfaces(cx, &mesh,
                &accepted.physical_boundary, &material, &source, temperature,
                current.frozen_interfaces()).unwrap();
            let residual = fs_conduction::assemble::full_residual(&system, temperature);
            let norm = residual.iter().map(|v| v * v).sum::<f64>().sqrt();
            near(norm, accepted.contact.residual_norm_w, 1e-11);
            assert!(norm <= coupling.residual_tolerance_w);
            assert!(accepted.contact.energy.closure_w.abs() <= coupling.energy_tolerance_w);
            let refused = contacts.solve_steady_with_ambient_radiation(cx, problem, &initial,
                solver.clone(), coupling, &[patch()], AmbientRadiationConfig {
                    max_iterations: 1, ..Default::default()
                });
            assert!(matches!(refused, Err(ConductionError::AmbientRadiationNotConverged { .. })));
            assert_eq!(initial, [300.0; 8]);
        }
    });
}

fn capacity_boundary(mesh: &ConductionMesh, convection: bool, bath: Option<f64>) -> ThermalBoundary {
    let builder = ThermalBoundaryBuilder::new(mesh)
        .region("exterior", |face| face.centroid[2] < 0.0 && face.centroid[1] == 0.0,
            if convection { ThermalBc::robin(2.0, 297.0).unwrap() }
            else { ThermalBc::neumann(0.0).unwrap() }).unwrap();
    let builder = if let Some(temperature) = bath {
        builder.region("bath", |face| face.centroid[2] > 0.0 && face.centroid[1] == 0.0,
            ThermalBc::dirichlet(temperature).unwrap()).unwrap()
    } else { builder };
    builder.adiabatic_remainder().finish().unwrap()
}

#[allow(clippy::too_many_arguments)]
fn check_capacity_physics(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, old: &[f64], dt: f64,
    solved: &fs_conduction::interface::temperature::capacity::ContactCapacityRadiationStep,
    convection: bool, source_w: f64,
) {
    let t = &solved.conduction.temperature;
    let current = contacts.evaluate(cx, problem.boundary, t).unwrap();
    assert_eq!(current.states(), solved.physical_contacts.states());
    let fresh = engine.evaluate_prescribed_endpoint(cx, problem,
        Some(current.frozen_interfaces()), old, t, dt, step_config(), Some(&[patch()])).unwrap();
    assert_eq!(fresh.residual_norm_j, solved.physical.residual_norm_j);
    assert!(fresh.residual_norm_j <= solved.physical_residual_tolerance_j);
    assert!(fresh.energy_residual_j.abs() <= step_config().energy_tolerance_j);
    let system = fs_conduction::assemble::assemble_operator_with_interfaces(
        cx, problem.mesh, &fresh.combined_boundary, problem.material,
        problem.source, t, current.frozen_interfaces(),
    ).unwrap();
    let spatial = fs_conduction::assemble::full_residual(&system, t);
    // Each original right tetra has volume 1/6. Its four lumped nodal
    // capacities are independently 10/24 or 20/24 J/K.
    let nodal_storage = t.iter().zip(old).enumerate().map(|(i, (&new, &old))|
        if i < 4 { 10.0 * (new - old) / 24.0 } else { 20.0 * (new - old) / 24.0 })
        .collect::<Vec<_>>();
    let storage: f64 = nodal_storage.iter().sum();
    let norm = spatial.iter().zip(&nodal_storage).enumerate()
        .filter(|(i, _)| !problem.boundary.dirichlet().iter().any(|&(v, _)| v == *i))
        .map(|(_, (&load, &stored))| (stored + dt * load).powi(2)).sum::<f64>().sqrt();
    near(norm, fresh.residual_norm_j, 1e-10);
    let reaction: f64 = problem.boundary.dirichlet().iter()
        .map(|&(vertex, _)| spatial[vertex] + nodal_storage[vertex] / dt).sum();
    near(reaction, fresh.dirichlet_in_w, 1e-10);
    let mean = (t[0] + t[1] + t[3]) / 3.0;
    let radiation = 0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5
        * (mean.powi(4) - 300.0_f64.powi(4));
    let convection_heat = if convection { 2.0 * 0.5 * (mean - 297.0) } else { 0.0 };
    near(fresh.source_w, source_w, 1e-12);
    near(fresh.radiation_out_w, radiation, 1e-10);
    near(solved.radiation.nonlinear_radiation_out_w, radiation, 1e-10);
    near(fresh.convective_out_w, convection_heat, 1e-10);
    assert_eq!(fresh.convective_robin_fluxes.is_empty(), !convection);
    near(fresh.stored_energy_change_j, storage, 1e-10);
    near(storage, dt * (source_w + reaction - convection_heat - radiation),
        step_config().energy_tolerance_j);
    let state = &current.states()[0];
    near(state.resistance.value_m2_k_per_w(),
        0.1 + 0.006 * (state.interface_temperature_k - 280.0), 1e-14);
    assert!(state.resistance.receipt().query_point.contains(&("normal_pressure".into(), 1e6)));
    assert_eq!(solved.radiation.solid_iterations, solved.nonlinear_updates);
    assert_eq!(solved.radiation.krylov_iterations, solved.krylov_iterations);
    assert!(solved.krylov_iterations >= solved.conduction.krylov_iterations);
    assert!(solved.krylov_iterations <= solved.thermal_responses * step_config().linear.max_iterations);
    assert!(solved.nonlinear_updates <= solved.thermal_responses * nonlinear().max_iterations);
    assert!(solved.nonlinear_backtracks <= solved.nonlinear_updates * nonlinear().line_search.max_backtracks);
}

#[test]
fn g1_capacity_contact_radiation_matches_fresh_nodal_and_independent_heat_accounts() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let curve = card(false);
        let surface = surface_card();
        let source = source();
        let material = material(true);
        let history = old();
        let before = history.clone();
        for convection in [false, true] {
            let boundary = capacity_boundary(&mesh, convection, None);
            let curve_contacts = contacts(cx, &mesh, &boundary, &curve);
            let surface_contacts = contacts(cx, &mesh, &boundary, &surface);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, source: &source, element_materials: None };
            let newton = curve_contacts.advance_capacity_newton_with_ambient_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
                &[patch()], AmbientRadiationConfig::default(),
            ).unwrap();
            let surface = surface_contacts.advance_capacity_with_ambient_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), Some(nonlinear()), coupling(),
                &[patch()], AmbientRadiationConfig::default(),
            ).unwrap();
            assert_eq!(history, before);
            assert!(newton.radiation.iterations > 1);
            assert_eq!(newton.thermal_responses, newton.radiation.iterations);
            assert!(surface.thermal_responses > surface.radiation.iterations);
            assert!(surface.thermal_responses <= coupling().max_iterations * surface.radiation.iterations);
            for (&a, &b) in newton.conduction.temperature.iter().zip(&surface.conduction.temperature) {
                near(a, b, 2e-7);
            }
            check_capacity_physics(cx, &engine, problem, &curve_contacts, &history, 0.35,
                &newton, convection, 5.0);
            check_capacity_physics(cx, &engine, problem, &surface_contacts, &history, 0.35,
                &surface, convection, 5.0);
        }
    });
}

#[test]
fn g1_capacity_contact_radiation_charges_prescribed_bath_storage_once_in_both_directions() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let card = card(false);
        let material = material(true);
        let source = ScalarField::Uniform(0.0);
        let history = [300.0; 8];
        for bath in [290.0, 350.0] {
            let boundary = capacity_boundary(&mesh, false, Some(bath));
            let contacts = contacts(cx, &mesh, &boundary, &card);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, source: &source, element_materials: None };
            let solved = contacts.advance_capacity_newton_with_ambient_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
                &[patch()], AmbientRadiationConfig::default(),
            ).unwrap();
            assert_eq!(history, [300.0; 8]);
            for &(vertex, prescribed) in boundary.dirichlet() {
                assert_eq!(solved.conduction.temperature[vertex], prescribed);
            }
            assert_eq!(solved.physical.dirichlet_in_w.is_sign_positive(), bath > 300.0);
            check_capacity_physics(cx, &engine, problem, &contacts, &history, 0.35,
                &solved, false, 0.0);
        }
    });
}

#[test]
fn g4_capacity_contact_radiation_refusal_preserves_history_and_exact_retry() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = capacity_boundary(&mesh, false, None);
        let card = card(false);
        let contacts = contacts(cx, &mesh, &boundary, &card);
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let material = material(true);
        let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, source: &source, element_materials: None };
        let history = old();
        let before = history.clone();
        let run = |cx: &Cx<'_>, radiation| contacts.advance_capacity_newton_with_ambient_radiation(
            cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
            &[patch()], radiation,
        );
        assert!(matches!(run(cx, AmbientRadiationConfig { max_iterations: 1, ..Default::default() }),
            Err(ConductionError::AmbientRadiationNotConverged { .. })));
        assert!(run(cx, AmbientRadiationConfig { max_iterations: usize::MAX, ..Default::default() }).is_err());
        let gate = CancelGate::new_clock_free();
        gate.request();
        context(&gate, |stopped| assert!(matches!(run(stopped, AmbientRadiationConfig::default()),
            Err(ConductionError::Cancelled { .. }))));
        assert_eq!(history, before);
        let retry = run(cx, AmbientRadiationConfig::default()).unwrap();
        let fresh = run(cx, AmbientRadiationConfig::default()).unwrap();
        assert_eq!(retry.conduction.temperature, fresh.conduction.temperature);
        assert_eq!(retry.physical_contacts.states(), fresh.physical_contacts.states());
        assert_eq!(retry.radiation, fresh.radiation);
        assert_eq!(retry.thermal_responses, fresh.thermal_responses);
        assert_eq!(history, before);
    });
}

fn closed_capacity_boundary(mesh: &ConductionMesh, convection: bool, bath: Option<f64>)
    -> ThermalBoundary {
    let builder = ThermalBoundaryBuilder::new(mesh)
        .region("lower", |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0,
            if convection { ThermalBc::robin(2.0, 297.0).unwrap() }
            else { ThermalBc::neumann(-0.4).unwrap() }).unwrap()
        .region("upper", |f| f.centroid[2] > 0.0 && f.centroid[1] == 0.0,
            ThermalBc::neumann(0.2).unwrap()).unwrap();
    let builder = if let Some(temperature) = bath {
        builder.region("bath", |f| f.centroid[2] > 0.0 && f.centroid[0] == 0.0,
            ThermalBc::dirichlet(temperature).unwrap()).unwrap()
    } else { builder };
    builder.adiabatic_remainder().finish().unwrap()
}

fn closed_capacity_enclosure(mesh: &ConductionMesh) -> fs_conduction::radiation::GrayDiffuseEnclosure {
    use fs_conduction::radiation::{GrayDiffuseEnclosure, RadiationSurface, ViewFactorMatrix};
    let lower = RadiationSurface::new(mesh, "lower",
        |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0, patch().emissivity().clone()).unwrap();
    let upper = RadiationSurface::new(mesh, "upper",
        |f| f.centroid[2] > 0.0 && f.centroid[1] == 0.0, patch().emissivity().clone()).unwrap();
    // Declared ideal reciprocal two-patch factors exercise the area-mean
    // model; this fixture does not claim geometric view-factor inference.
    GrayDiffuseEnclosure::new(vec![lower, upper],
        ViewFactorMatrix::infinite_parallel_plates(0.5).unwrap()).unwrap()
}

#[allow(clippy::too_many_arguments)]
fn check_closed_capacity_physics(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, old: &[f64],
    enclosure: &fs_conduction::radiation::GrayDiffuseEnclosure,
    solved: &fs_conduction::interface::temperature::capacity::ContactCapacityEnclosureStep,
    convection: bool, source_w: f64,
) {
    let dt = 0.35;
    let t = &solved.conduction.temperature;
    let current = contacts.evaluate(cx, problem.boundary, t).unwrap();
    assert_eq!(current.states(), solved.physical_contacts.states());
    let audit = engine.evaluate_prescribed_enclosure_endpoint(cx, problem,
        Some(current.frozen_interfaces()), old, t, dt, step_config(), enclosure,
        AmbientRadiationConfig::default()).unwrap();
    assert_eq!(audit.radiosity, solved.radiation.radiosity);
    assert_eq!(audit.physical.residual_norm_j, solved.physical.residual_norm_j);
    assert!(audit.physical.residual_norm_j <= coupling().residual_tolerance_j);
    assert!(audit.physical.energy_residual_j.abs() <= step_config().energy_tolerance_j);
    let system = fs_conduction::assemble::assemble_operator_with_interfaces(cx,
        problem.mesh, &audit.physical.combined_boundary, problem.material, problem.source,
        t, current.frozen_interfaces()).unwrap();
    let spatial = fs_conduction::assemble::full_residual(&system, t);
    let stored = t.iter().zip(old).enumerate().map(|(i, (&new, &old))|
        if i < 4 { 10.0 * (new - old) / 24.0 } else { 20.0 * (new - old) / 24.0 })
        .collect::<Vec<_>>();
    let storage: f64 = stored.iter().sum();
    let norm = spatial.iter().zip(&stored).enumerate()
        .filter(|(i, _)| !problem.boundary.dirichlet().iter().any(|&(v, _)| v == *i))
        .map(|(_, (&load, &storage))| (storage + dt * load).powi(2)).sum::<f64>().sqrt();
    near(norm, audit.physical.residual_norm_j, 1e-10);
    let reaction: f64 = problem.boundary.dirichlet().iter()
        .map(|&(v, _)| spatial[v] + stored[v] / dt).sum();
    near(reaction, audit.physical.dirichlet_in_w, 1e-10);
    let means = [(t[0] + t[1] + t[3]) / 3.0, (t[4] + t[5] + t[7]) / 3.0];
    let q = STEFAN_BOLTZMANN_W_M2_K4 * 0.5
        * (means[0].powi(4) - means[1].powi(4)) / (1.0 / 0.2 + 1.0 / 0.2 - 1.0);
    near(audit.radiosity.net_outward_heat_w[0], q, 1e-10);
    near(audit.radiosity.net_outward_heat_w[1], -q, 1e-10);
    near(audit.physical.radiation_out_w, 0.0, 1e-10);
    let convective = if convection { 2.0 * 0.5 * (means[0] - 297.0) } else { 0.0 };
    let neumann = 0.5 * (0.2 + if convection { 0.0 } else { -0.4 });
    near(audit.physical.convective_out_w, convective, 1e-10);
    near(audit.physical.neumann_out_w, neumann, 1e-12);
    near(audit.physical.source_w, source_w, 1e-12);
    assert_eq!(audit.physical.convective_robin_fluxes.is_empty(), !convection);
    near(audit.physical.stored_energy_change_j, storage, 1e-10);
    near(storage, dt * (source_w + reaction - neumann - convective),
        step_config().energy_tolerance_j);
    let state = &current.states()[0];
    near(state.resistance.value_m2_k_per_w(),
        0.1 + 0.006 * (state.interface_temperature_k - 280.0), 1e-14);
    assert!(state.resistance.receipt().query_point.contains(&("normal_pressure".into(), 1e6)));
    assert_eq!(solved.radiation.solid_iterations, solved.nonlinear_updates);
    assert_eq!(solved.radiation.krylov_iterations, solved.krylov_iterations);
    assert!(solved.krylov_iterations >= solved.conduction.krylov_iterations);
    assert!(solved.krylov_iterations <= solved.thermal_responses * step_config().linear.max_iterations);
    assert!(solved.nonlinear_updates <= solved.thermal_responses * nonlinear().max_iterations);
    assert!(solved.nonlinear_backtracks <= solved.nonlinear_updates * nonlinear().line_search.max_backtracks);
}

#[test]
fn g1_capacity_live_contact_enclosure_preserves_fresh_balance_and_internal_heat() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let enclosure = closed_capacity_enclosure(&mesh);
        let curve = card(false);
        let surface = surface_card();
        let source = source();
        let material = material(true);
        let history = old();
        let before = history.clone();
        for convection in [false, true] {
            let boundary = closed_capacity_boundary(&mesh, convection, None);
            let curve_contacts = contacts(cx, &mesh, &boundary, &curve);
            let surface_contacts = contacts(cx, &mesh, &boundary, &surface);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, source: &source, element_materials: None };
            let newton = curve_contacts.advance_capacity_newton_with_enclosure_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
                &enclosure, AmbientRadiationConfig::default(),
            ).unwrap();
            let fixed_point = surface_contacts.advance_capacity_with_enclosure_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), Some(nonlinear()), coupling(),
                &enclosure, AmbientRadiationConfig::default(),
            ).unwrap();
            assert_eq!(history, before);
            assert!(newton.radiation.iterations > 1);
            assert_eq!(newton.thermal_responses, newton.radiation.iterations);
            assert!(fixed_point.thermal_responses > fixed_point.radiation.iterations);
            assert!(fixed_point.thermal_responses <= coupling().max_iterations * fixed_point.radiation.iterations);
            for (&a, &b) in newton.conduction.temperature.iter().zip(&fixed_point.conduction.temperature) {
                near(a, b, 2e-7);
            }
            check_closed_capacity_physics(cx, &engine, problem, &curve_contacts, &history,
                &enclosure, &newton, convection, 5.0);
            check_closed_capacity_physics(cx, &engine, problem, &surface_contacts, &history,
                &enclosure, &fixed_point, convection, 5.0);
        }
    });
}

#[test]
fn g1_capacity_live_contact_enclosure_counts_both_bath_jump_directions_once() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let enclosure = closed_capacity_enclosure(&mesh);
        let card = card(false);
        let material = material(true);
        let source = ScalarField::Uniform(0.0);
        let history = [300.0; 8];
        for bath in [290.0, 350.0] {
            let boundary = closed_capacity_boundary(&mesh, false, Some(bath));
            let contacts = contacts(cx, &mesh, &boundary, &card);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, source: &source, element_materials: None };
            let solved = contacts.advance_capacity_newton_with_enclosure_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
                &enclosure, AmbientRadiationConfig::default(),
            ).unwrap();
            assert_eq!(history, [300.0; 8]);
            for &(vertex, prescribed) in boundary.dirichlet() {
                assert_eq!(solved.conduction.temperature[vertex], prescribed);
            }
            assert_eq!(solved.physical.dirichlet_in_w.is_sign_positive(), bath > 300.0);
            check_closed_capacity_physics(cx, &engine, problem, &contacts, &history,
                &enclosure, &solved, false, 0.0);
        }
    });
}

#[test]
fn g4_capacity_live_contact_enclosure_refusal_preserves_history_and_retry() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = closed_capacity_boundary(&mesh, true, None);
        let enclosure = closed_capacity_enclosure(&mesh);
        let card = card(false);
        let contacts = contacts(cx, &mesh, &boundary, &card);
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let material = material(true);
        let source = source();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, source: &source, element_materials: None };
        let history = old();
        let before = history.clone();
        let run = |cx: &Cx<'_>, radiation| contacts.advance_capacity_newton_with_enclosure_radiation(
            cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
            &enclosure, radiation,
        );
        assert!(matches!(run(cx, AmbientRadiationConfig { max_iterations: 1, ..Default::default() }),
            Err(ConductionError::Radiation { .. })));
        assert!(run(cx, AmbientRadiationConfig { max_iterations: usize::MAX, ..Default::default() }).is_err());
        let gate = CancelGate::new_clock_free();
        gate.request();
        context(&gate, |stopped| assert!(matches!(run(stopped, AmbientRadiationConfig::default()),
            Err(ConductionError::Cancelled { .. }))));
        assert_eq!(history, before);
        let retry = run(cx, AmbientRadiationConfig::default()).unwrap();
        let fresh = run(cx, AmbientRadiationConfig::default()).unwrap();
        assert_eq!(retry.conduction.temperature, fresh.conduction.temperature);
        assert_eq!(retry.physical_contacts.states(), fresh.physical_contacts.states());
        assert_eq!(retry.radiation.radiosity, fresh.radiation.radiosity);
        assert_eq!(retry.radiation.applied_heat_w, fresh.radiation.applied_heat_w);
        assert_eq!(retry.thermal_responses, fresh.thermal_responses);
        assert_eq!(history, before);
    });
}

fn open_capacity_enclosure(cx: &Cx<'_>, mesh: &ConductionMesh)
    -> fs_conduction::radiation::OpenDiffuseEnclosure {
    use fs_conduction::radiation::{OpenDiffuseEnclosure, ViewFactorEvidence, ViewFactorTolerance};
    OpenDiffuseEnclosure::new(cx, closed_capacity_enclosure(mesh).surfaces().to_vec(),
        vec![vec![0.0, 0.25], vec![0.25, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "manufactured reciprocal open exchange".into() },
        ViewFactorTolerance { row_sum_abs: 1e-12, reciprocity_rel: 1e-12 }).unwrap()
}

#[allow(clippy::too_many_arguments)]
fn check_open_capacity_physics(
    cx: &Cx<'_>, engine: &BackwardEuler<'_>, problem: ConductionProblem<'_>,
    contacts: &TemperatureContacts<'_, '_>, old: &[f64],
    model: &fs_conduction::radiation::OpenDiffuseEnclosure, ambient: f64,
    solved: &fs_conduction::interface::temperature::capacity::ContactCapacityEnclosureStep,
    convection: bool, source_w: f64,
) {
    let dt = 0.35;
    let t = &solved.conduction.temperature;
    let current = contacts.evaluate(cx, problem.boundary, t).unwrap();
    assert_eq!(current.states(), solved.physical_contacts.states());
    let audit = engine.evaluate_prescribed_open_enclosure_endpoint(cx, problem,
        Some(current.frozen_interfaces()), old, t, dt, step_config(), model, ambient,
        AmbientRadiationConfig::default()).unwrap();
    assert_eq!(audit.radiosity, solved.radiation.radiosity);
    assert_eq!(audit.physical.residual_norm_j, solved.physical.residual_norm_j);
    assert!(audit.physical.residual_norm_j <= coupling().residual_tolerance_j);
    assert!(audit.physical.energy_residual_j.abs() <= step_config().energy_tolerance_j);
    assert_eq!(solved.radiation.applied_heat_w.len(), 2, "only physical traces receive loads");
    assert_eq!(audit.radiosity.net_outward_heat_w.len(), 3, "the bath keeps its own row");

    // Independent two-by-two gray radiosity solve relative to the black bath.
    // Surface area is 0.5 m2, epsilon=0.2, and mutual view factor is 0.25.
    let means = [(t[0] + t[1] + t[3]) / 3.0, (t[4] + t[5] + t[7]) / 3.0];
    let emission = means.map(|mean| 0.2 * STEFAN_BOLTZMANN_W_M2_K4
        * (mean - ambient) * (mean + ambient) * (mean * mean + ambient * ambient));
    let reflected = 0.8 * 0.25;
    let radiosity = [(emission[0] + reflected * emission[1]) / (1.0 - reflected * reflected),
        (emission[1] + reflected * emission[0]) / (1.0 - reflected * reflected)];
    let heat = [0.5 * (radiosity[0] - 0.25 * radiosity[1]),
        0.5 * (radiosity[1] - 0.25 * radiosity[0])];
    let radiation: f64 = heat.iter().sum();
    for (i, q) in heat.iter().enumerate() { near(audit.radiosity.net_outward_heat_w[i], *q, 1e-10); }
    near(-audit.radiosity.net_outward_heat_w[2], radiation, 1e-10);
    near(audit.physical.radiation_out_w, radiation, 1e-10);

    // Assemble original conduction/contact and add physical surface watts
    // directly, independently of the production shifted-boundary adapter.
    let system = fs_conduction::assemble::assemble_operator_with_interfaces(cx,
        problem.mesh, problem.boundary, problem.material, problem.source,
        t, current.frozen_interfaces()).unwrap();
    let mut spatial = fs_conduction::assemble::full_residual(&system, t);
    for (vertices, q) in [([0, 1, 3], heat[0]), ([4, 5, 7], heat[1])] {
        for vertex in vertices { spatial[vertex] += q / 3.0; }
    }
    let stored = t.iter().zip(old).enumerate().map(|(i, (&new, &old))|
        if i < 4 { 10.0 * (new - old) / 24.0 } else { 20.0 * (new - old) / 24.0 })
        .collect::<Vec<_>>();
    let storage: f64 = stored.iter().sum();
    let residual = spatial.iter().zip(&stored).enumerate()
        .filter(|(i, _)| !problem.boundary.dirichlet().iter().any(|&(v, _)| v == *i))
        .map(|(_, (&load, &storage))| (storage + dt * load).powi(2)).sum::<f64>().sqrt();
    near(residual, audit.physical.residual_norm_j, 1e-10);
    let reaction: f64 = problem.boundary.dirichlet().iter()
        .map(|&(v, _)| spatial[v] + stored[v] / dt).sum();
    near(reaction, audit.physical.dirichlet_in_w, 1e-10);
    let convective = if convection { means[0] - 297.0 } else { 0.0 };
    let neumann = 0.5 * (0.2 + if convection { 0.0 } else { -0.4 });
    near(audit.physical.convective_out_w, convective, 1e-10);
    near(audit.physical.neumann_out_w, neumann, 1e-12);
    near(audit.physical.source_w, source_w, 1e-12);
    near(audit.physical.stored_energy_change_j, storage, 1e-10);
    near(storage, dt * (source_w + reaction - neumann - convective - radiation),
        step_config().energy_tolerance_j);
    assert!(radiation.abs() > 1.0, "a zero external radiation account would fail this balance");
    assert_eq!(solved.radiation.solid_iterations, solved.nonlinear_updates);
    assert_eq!(solved.radiation.krylov_iterations, solved.krylov_iterations);
    assert!(solved.krylov_iterations <= solved.thermal_responses * step_config().linear.max_iterations);
    assert!(solved.nonlinear_updates <= solved.thermal_responses * nonlinear().max_iterations);
}

#[test]
fn g1_open_capacity_contact_keeps_storage_and_black_reservoir_heat_separate() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let model = open_capacity_enclosure(cx, &mesh);
        let curve = card(false); let surface = surface_card();
        let source = source(); let material = material(true); let history = old();
        for convection in [false, true] {
            let boundary = closed_capacity_boundary(&mesh, convection, None);
            let curve_contacts = contacts(cx, &mesh, &boundary, &curve);
            let surface_contacts = contacts(cx, &mesh, &boundary, &surface);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, source: &source, element_materials: None };
            for ambient in [290.0, 370.0] {
                let newton = curve_contacts.advance_capacity_newton_with_open_enclosure_radiation(
                    cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
                    &model, ambient, AmbientRadiationConfig::default()).unwrap();
                let picard = surface_contacts.advance_capacity_with_open_enclosure_radiation(
                    cx, &engine, problem, &history, 0.35, step_config(), Some(nonlinear()), coupling(),
                    &model, ambient, AmbientRadiationConfig::default()).unwrap();
                assert_eq!(newton.thermal_responses, newton.radiation.iterations);
                assert!(picard.thermal_responses > picard.radiation.iterations);
                assert!(picard.thermal_responses <= coupling().max_iterations * picard.radiation.iterations);
                assert_eq!(newton.physical.radiation_out_w.is_sign_positive(), ambient == 290.0);
                for (&a, &b) in newton.conduction.temperature.iter().zip(&picard.conduction.temperature) {
                    near(a, b, 2e-7);
                }
                check_open_capacity_physics(cx, &engine, problem, &curve_contacts, &history,
                    &model, ambient, &newton, convection, 5.0);
                check_open_capacity_physics(cx, &engine, problem, &surface_contacts, &history,
                    &model, ambient, &picard, convection, 5.0);
                if !convection && ambient == 290.0 {
                    let fixed = engine.advance_prescribed_with_open_enclosure_radiation(
                        cx, problem, Some(newton.physical_contacts.frozen_interfaces()), &history,
                        0.35, step_config(), Some(nonlinear()), &model, ambient,
                        AmbientRadiationConfig::default()).unwrap();
                    for (&a, &b) in newton.conduction.temperature.iter().zip(&fixed.conduction.temperature) {
                        near(a, b, 2e-7);
                    }
                }
            }
        }
        assert_eq!(history, old());
    });
}

#[test]
fn g1_open_capacity_contact_keeps_prescribed_node_storage_and_radiation_reactions() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let model = open_capacity_enclosure(cx, &mesh);
        let card = card(false); let material = material(true);
        let source = ScalarField::Uniform(0.0); let history = [300.0; 8];
        for bath in [290.0, 350.0] {
            let boundary = closed_capacity_boundary(&mesh, false, Some(bath));
            let contacts = contacts(cx, &mesh, &boundary, &card);
            let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
                material: &material, source: &source, element_materials: None };
            let solved = contacts.advance_capacity_newton_with_open_enclosure_radiation(
                cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
                &model, 280.0, AmbientRadiationConfig::default()).unwrap();
            for &(vertex, prescribed) in boundary.dirichlet() {
                assert_eq!(solved.conduction.temperature[vertex], prescribed);
            }
            check_open_capacity_physics(cx, &engine, problem, &contacts, &history,
                &model, 280.0, &solved, false, 0.0);
        }
        assert_eq!(history, [300.0; 8]);
    });
}

#[test]
fn g4_open_capacity_contact_refusals_keep_history_and_retry_unchanged() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let boundary = closed_capacity_boundary(&mesh, true, None);
        let model = open_capacity_enclosure(cx, &mesh); let card = card(false);
        let contacts = contacts(cx, &mesh, &boundary, &card);
        let engine = BackwardEuler::per_element(cx, &mesh, &[
            VolumetricHeatCapacity::declared(10.0).unwrap(),
            VolumetricHeatCapacity::declared(20.0).unwrap(),
        ]).unwrap();
        let material = material(true); let source = source(); let history = old();
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, source: &source, element_materials: None };
        let run = |cx: &Cx<'_>, ambient, radiation| contacts.advance_capacity_newton_with_open_enclosure_radiation(
            cx, &engine, problem, &history, 0.35, step_config(), nonlinear(), coupling(), 1,
            &model, ambient, radiation);
        for invalid in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
            assert!(run(cx, invalid, AmbientRadiationConfig::default()).is_err());
        }
        assert!(run(cx, 290.0, AmbientRadiationConfig { max_iterations: 1,
            temperature_tolerance_k: 1e6, ..Default::default() }).is_err());
        assert!(run(cx, 290.0, AmbientRadiationConfig { max_iterations: usize::MAX, ..Default::default() }).is_err());
        let outside = [400.0; 8];
        assert!(contacts.advance_capacity_newton_with_open_enclosure_radiation(
            cx, &engine, problem, &outside, 0.35, step_config(), nonlinear(), coupling(), 1,
            &model, 290.0, AmbientRadiationConfig::default()).is_err());
        assert_eq!(outside, [400.0; 8]);
        let gate = CancelGate::new_clock_free(); gate.request();
        context(&gate, |stopped| assert!(matches!(run(stopped, 290.0, AmbientRadiationConfig::default()),
            Err(ConductionError::Cancelled { .. }))));
        let retry = run(cx, 290.0, AmbientRadiationConfig::default()).unwrap();
        let fresh = run(cx, 290.0, AmbientRadiationConfig::default()).unwrap();
        assert_eq!(retry.conduction.temperature, fresh.conduction.temperature);
        assert_eq!(retry.physical_contacts.states(), fresh.physical_contacts.states());
        assert_eq!(retry.radiation.radiosity, fresh.radiation.radiosity);
        assert_eq!(retry.thermal_responses, fresh.thermal_responses);
        assert_eq!(history, old());
    });
}
