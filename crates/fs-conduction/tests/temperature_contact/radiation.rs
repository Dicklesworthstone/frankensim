//! Fresh contact and exact radiative heat must balance on the returned field.

use super::*;
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
        let surface = InterfaceSystemCard::assemble(original.surface_a().clone(),
            original.surface_b().clone(), original.context().clone(), claims, vec![]).unwrap();
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
