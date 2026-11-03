//! Fresh contact and exact radiative heat must balance on the returned field.

use super::*;
use fs_conduction::interface::temperature::steady::SteadyContactNewtonConfig;
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
