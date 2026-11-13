//! Coupled finite-vapor/swelling references through the public native owners.
//! The coefficients are authored numerical fixtures, not measured wood data.
use super::*;
use crate::atmosphere as air;
use fs_couple::moisture::atmosphere::AtmosphereMoistureError;
use fs_couple::moisture::solid::MoistureSpatialAtmosphereControls;

fn coupling_controls() -> MoistureSpatialAtmosphereControls {
    MoistureSpatialAtmosphereControls {
        atmosphere: air::controls(),
        maximum_geometry_residual_m: 1e-9,
        maximum_transport_evaluations: 100,
        maximum_native_solid_solves: 4096,
    }
}

#[test]
fn g1_finite_vapor_and_swelling_solve_the_same_endpoint_water_equations() {
    with_cx(|cx| {
        // First case isolates bulk 1/lambda_x²; the second also resolves the
        // surface resistance using the actual changing transverse area.
        for (beta, h, dt) in [
            ([0.1, 0.0, 0.0], None, 0.25),
            ([0.1, 0.08, 0.04], Some(0.4), 20.0),
        ] {
            let mut stack = one_cell(beta, cx);
            if let Some(h) = h {
                stack = stack
                    .with_transfer(
                        [MoistureActivityTransfer::ConductanceKgM2S(h); 2],
                        vec![],
                        cx,
                    )
                    .unwrap();
            }
            let solid = bind_stack(&stack, free_supports(), cx);
            let gas = air::atmosphere(0.8, 1.0);
            let saturation = gas.vapor_mass_kg().unwrap() / gas.relative_humidity().unwrap();
            let dry = stack.dry_mass_kg();
            let expected = root(0.1, 0.24, |m| {
                let stretch = beta.map(|b| 1.0 + b * (m - 0.1));
                let area = 0.001 * stretch[1] * stretch[2];
                let half = 2.0 * dry * 0.01 * 0.3 / (0.1 * stretch[0]).powi(2);
                let conductance = h.map_or(half, |h| 1.0 / (1.0 / half + 1.0 / (h * area)));
                let rh = 0.8 - dry * (m - 0.1) / saturation;
                dry * (m - 0.1) - 2.0 * dt * conductance * (rh - m / 0.3)
            });
            let (candidate, state) = solid
                .advance_spatial_diffusion_finite_atmosphere(
                    &gas,
                    &stack,
                    dt,
                    [true; 2],
                    config(),
                    nonlinear(),
                    coupling_controls(),
                    cx,
                )
                .unwrap();
            let m = candidate.moisture().layers()[0].slab().moisture()[0];
            let step = candidate.step();
            let spatial = &step.transport.spatial;
            near(m, expected, 4e-10);
            let expected_rh = 0.8 - dry * (expected - 0.1) / saturation;
            near(
                candidate.atmosphere().relative_humidity().unwrap(),
                expected_rh,
                2e-9,
            );
            near(spatial.axial_stretch[0], 1.0 + beta[0] * (m - 0.1), 2e-10);
            near(
                spatial.current_cross_section_m2,
                0.001 * (1.0 + beta[1] * (m - 0.1)) * (1.0 + beta[2] * (m - 0.1)),
                2e-12,
            );
            near(step.water_to_solid_kg, dry * (m - 0.1), 2e-12);
            near(
                spatial.transport.boundary_water_kg.iter().sum(),
                step.water_to_solid_kg,
                2e-12,
            );
            near(
                candidate.atmosphere().vapor_mass_kg().unwrap()
                    + candidate.moisture().water_mass_kg(),
                gas.vapor_mass_kg().unwrap() + stack.water_mass_kg(),
                2e-12,
            );
            near(
                air::carrier_mass(candidate.atmosphere()),
                air::carrier_mass(&gas),
                1e-12,
            );
            near(
                candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                gas.gas().mass_kg() + stack.wet_mass_kg(),
                2e-12,
            );
            assert_eq!(
                candidate.atmosphere().gas().volume_m3(),
                gas.gas().volume_m3()
            );
            assert!(m < 0.3 * expected_rh); // finite-time uptake, not bulk equilibrium
            near(
                state.wet_mass_kg(),
                candidate.moisture().wet_mass_kg(),
                1e-12,
            );
            assert_eq!(state.geometry().identity(), spatial.geometry_identity);
            assert_eq!(state.time_s(), candidate.moisture().time_s());
            assert!(step.transport.transport_evaluations > 1);
            assert!(step.transport.native_solid_solves > spatial.native_solid_solves);
            assert!(
                step.transport.native_solid_solves
                    <= coupling_controls().maximum_native_solid_solves
            );
            for side in 0..2 {
                let source = step.sorption[side].as_ref().unwrap();
                let activity = spatial.transport.surface_activities[side].unwrap();
                near(
                    source.state_coordinate("relative_humidity").unwrap(),
                    activity,
                    0.0,
                );
                assert_eq!(
                    source.card_identity(),
                    stack.layers()[0].slab().card().content_hash()
                );
                if let Some(h) = h {
                    near(
                        spatial.transport.boundary_water_kg[side],
                        dt * h * spatial.current_cross_section_m2 * (expected_rh - activity),
                        3e-12,
                    );
                }
            }
        }
    });
}

#[test]
fn g1_zero_swelling_reduces_to_original_finite_atmosphere_and_sealed_gas_is_exact() {
    with_cx(|cx| {
        let stack = one_cell([0.0; 3], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let gas = air::atmosphere(0.8, 1.0);
        let (ordinary, _) = solid
            .advance_finite_atmosphere(
                &gas,
                &stack,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                air::controls(),
                cx,
            )
            .unwrap();
        let (spatial, state) = solid
            .advance_spatial_diffusion_finite_atmosphere(
                &gas,
                &stack,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                coupling_controls(),
                cx,
            )
            .unwrap();
        near(
            spatial.moisture().layers()[0].slab().moisture()[0],
            ordinary.moisture().layers()[0].slab().moisture()[0],
            2e-12,
        );
        near(
            spatial.atmosphere().relative_humidity().unwrap(),
            ordinary.atmosphere().relative_humidity().unwrap(),
            2e-12,
        );
        assert_eq!(spatial.step().transport.spatial.axial_stretch, vec![1.0]);
        near(state.stress().elastic_energy_j(), 0.0, 1e-14);

        let sealed_stack = one_cell([0.1, 0.0, 0.0], cx)
            .with_transfer([MoistureActivityTransfer::Sealed; 2], vec![], cx)
            .unwrap();
        let sealed_solid = bind_stack(&sealed_stack, free_supports(), cx);
        let (sealed, _) = sealed_solid
            .advance_spatial_diffusion_finite_atmosphere(
                &gas,
                &sealed_stack,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                coupling_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(sealed.atmosphere(), &gas);
        assert_eq!(sealed.step().water_to_solid_kg, 0.0);
        assert_eq!(sealed.step().transport.transport_evaluations, 1);
        assert_eq!(
            sealed.step().transport.spatial.transport.boundary_water_kg,
            [0.0; 2]
        );
    });
}

#[test]
fn g3_wet_dry_reservoir_cycles_preserve_each_inventory_and_source_reference_shape() {
    with_cx(|cx| {
        let beta = [0.1, 0.0, 0.0];
        let mut stack = one_cell(beta, cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        // Two separately retained finite reservoirs are alternately connected.
        // Only the exposed reservoir exchanges water; neither is reset to RH.
        let mut reservoirs = [air::atmosphere(0.9, 1.0), air::atmosphere(0.05, 1.0)];
        let initial_water = reservoirs
            .iter()
            .map(|gas| gas.vapor_mass_kg().unwrap())
            .sum::<f64>()
            + stack.water_mass_kg();
        let carriers = reservoirs.each_ref().map(air::carrier_mass);
        for cycle in 0..4 {
            let side = cycle % 2;
            let old_m = stack.layers()[0].slab().moisture()[0];
            let (candidate, state) = solid
                .advance_spatial_diffusion_finite_atmosphere(
                    &reservoirs[side],
                    &stack,
                    0.25,
                    [true; 2],
                    config(),
                    nonlinear(),
                    coupling_controls(),
                    cx,
                )
                .unwrap();
            let (gas, next, step) = candidate.into_parts();
            let m = next.layers()[0].slab().moisture()[0];
            assert_eq!(m > old_m, side == 0);
            near(
                step.transport.spatial.axial_stretch[0],
                1.0 + beta[0] * (m - 0.1),
                2e-10,
            );
            near(state.wet_mass_kg(), next.wet_mass_kg(), 2e-12);
            assert!(state.tangent().is_some());
            reservoirs[side] = gas;
            stack = next;
            near(
                reservoirs
                    .iter()
                    .map(|gas| gas.vapor_mass_kg().unwrap())
                    .sum::<f64>()
                    + stack.water_mass_kg(),
                initial_water,
                5e-12,
            );
            for (gas, carrier) in reservoirs.iter().zip(carriers) {
                near(air::carrier_mass(gas), carrier, 1e-12);
            }
            assert_eq!(stack.time_s(), (cycle + 1) as f64 * 0.25);
        }
    });
}

#[test]
fn g3_cumulative_rejected_trial_budgets_and_late_failures_never_advance_inputs() {
    with_cx(|cx| {
        let stack = one_cell([0.1, 0.0, 0.0], cx);
        let solid = bind_stack(&stack, free_supports(), cx);
        let gas = air::atmosphere(0.8, 1.0);
        let original_gas = gas.clone();
        let run = |control| {
            solid.advance_spatial_diffusion_finite_atmosphere(
                &gas,
                &stack,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                control,
                cx,
            )
        };
        let (accepted, _) = run(coupling_controls()).unwrap();
        let total_work = accepted.step().transport.native_solid_solves;
        let mut bounded = coupling_controls();
        bounded.maximum_transport_evaluations = 1;
        assert!(
            matches!(run(bounded), Err(AtmosphereMoistureError::Invalid(message))
            if message.contains("transport-evaluation budget"))
        );
        bounded = coupling_controls();
        bounded.maximum_native_solid_solves = total_work - 1;
        assert!(matches!(run(bounded),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Invalid(message)))
                if message.contains("solid-solve budget")));
        bounded = coupling_controls();
        bounded.atmosphere.maximum_iterations = 1;
        assert!(matches!(
            run(bounded),
            Err(AtmosphereMoistureError::NotConverged { .. })
        ));
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(stack.time_s(), 0.0);
        assert_eq!(gas, original_gas);
        let (retry, _) = run(coupling_controls()).unwrap();
        assert_eq!(retry.moisture().activity(), accepted.moisture().activity());
        assert_eq!(retry.atmosphere(), accepted.atmosphere());

        let warped = one_cell([0.01, 0.08, 0.04], cx);
        let clamped = bind_stack(&warped, (0..12).collect(), cx);
        assert!(
            matches!(clamped.advance_spatial_diffusion_finite_atmosphere(
            &gas, &warped, 0.25, [true; 2], config(), nonlinear(), coupling_controls(), cx,
        ), Err(AtmosphereMoistureError::Moisture(MoistureError::Invalid(message)))
            if message.contains("axis-aligned"))
        );
        assert_eq!(warped.layers()[0].slab().moisture(), &[0.1]);

        let restricted = MoistureLayeredSlab::new(
            vec![layer(
                material(500.0, [0.1, 0.0, 0.0], 0.1001, 0.3),
                0.1,
                vec![0.1],
                cx,
            )],
            16,
            cx,
        )
        .unwrap();
        let limited_solid = bind_stack(&restricted, free_supports(), cx);
        assert!(
            limited_solid
                .advance_spatial_diffusion_finite_atmosphere(
                    &gas,
                    &restricted,
                    0.25,
                    [true; 2],
                    config(),
                    nonlinear(),
                    coupling_controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(restricted.layers()[0].slab().moisture(), &[0.1]);
        assert_eq!(restricted.time_s(), 0.0);
        assert_eq!(gas, original_gas);
    });
}

#[test]
fn g0_cancelled_spatial_atmosphere_keeps_gas_water_and_mechanics_at_their_start() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let stack = one_cell([0.1, 0.0, 0.0], &cx);
        let solid = bind_stack(&stack, free_supports(), &cx);
        let gas = air::atmosphere(0.8, 1.0);
        let before = gas.clone();
        gate.request();
        assert!(matches!(
            solid.advance_spatial_diffusion_finite_atmosphere(
                &gas,
                &stack,
                0.25,
                [true; 2],
                config(),
                nonlinear(),
                coupling_controls(),
                &cx,
            ),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert_eq!(gas, before);
        assert_eq!(stack.time_s(), 0.0);
        assert_eq!(stack.layers()[0].slab().moisture(), &[0.1]);
    });
}
