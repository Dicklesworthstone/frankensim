//! Physical-time refinement of retained gas/wall compression, not a second EOS.
use super::*;
use fs_couple::solid_enthalpy::finite_plastic::porous::cavity::heat::march::{
    AdaptiveCavityHeatSolve, CavityHeatMarchError,
};

fn adaptive() -> AdaptiveCavityHeatSolve {
    AdaptiveCavityHeatSolve {
        initial_step_s: 1.0,
        minimum_step_s: 1e-5,
        maximum_step_s: 1.0,
        maximum_intervals: 512,
        maximum_heat_trials: 32_768,
        absolute_displacement_tolerance_m: 1e-8,
        absolute_velocity_tolerance_mps: 1e-7,
        absolute_temperature_tolerance_k: 0.005,
        absolute_specific_enthalpy_tolerance_j_kg: 0.005 * CP,
        absolute_gas_volume_tolerance_m3: 1e-12,
        absolute_plastic_state_tolerance: 1e-5,
        absolute_energy_tolerance_j: 1e-5,
        relative_tolerance: 1e-3,
        window_energy_tolerance_j: 1e-5,
    }
}

#[test]
fn rigid_heat_refines_toward_the_independent_continuous_two_capacity_solution() {
    context(|cx| {
        let fixture = thermally_constant_fixture();
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let (gas, _) = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.0, 1e24)
            .at_temperature(360.0)
            .unwrap();
        let gas_capacity =
            gas.mass_kg() * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028 * 2.5;
        let wall_capacity: f64 = thermal.cells().values().map(|cell| cell.mass_kg * CP).sum();
        let rate = 0.4;
        let system = ThermalPorousGasSystem::try_new(
            solid,
            &thermal,
            vec![cavity_entry(Arc::new(cavity), gas)],
            links(&thermal, rate),
            cx,
        )
        .unwrap();
        let supports = affine(&nodes, 1.0);
        let observer = Times::default();
        let external = 0.02;
        let coarse = system
            .advance(
                &thermal,
                &observer,
                1.0,
                &supports,
                &[],
                [0.0; 3],
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                &[external],
                &[0.0],
                solve_heat(),
                cx,
            )
            .unwrap();
        let refined = system
            .advance_adaptive(
                &thermal,
                &observer,
                1.0,
                &supports,
                &[],
                [0.0; 3],
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                &[external],
                &[0.0],
                solve_heat(),
                adaptive(),
                cx,
            )
            .unwrap();
        // Closed ODE solution: Cg*Tg' = P-G*(Tg-Tw), Cw*Tw' = G*(Tg-Tw).
        // The supplied whole-stroke external heat is constant power for 1 s.
        let conductance = rate * wall_capacity;
        let decay_rate = conductance * (1.0 / gas_capacity + 1.0 / wall_capacity);
        let decay = (-decay_rate).exp();
        let difference = 60.0 * decay + external / (gas_capacity * decay_rate) * (1.0 - decay);
        let energy = gas_capacity * 360.0 + wall_capacity * 300.0 + external;
        let expected_gas = (energy + wall_capacity * difference) / (gas_capacity + wall_capacity);
        let expected_wall = (energy - gas_capacity * difference) / (gas_capacity + wall_capacity);
        let gas_temperature = refined.endpoint.system.cavities()[0]
            .gas
            .state()
            .state()
            .temperature;
        let coarse_error =
            (coarse.system.cavities()[0].gas.state().state().temperature - expected_gas).abs();
        assert!((gas_temperature - expected_gas).abs() < 0.2 * coarse_error);
        close(gas_temperature, expected_gas, 0.15);
        for cell in refined.endpoint.thermal.cells().values() {
            close(
                refined.endpoint.thermal.temperature_k(cell.id).unwrap(),
                expected_wall,
                0.15,
            );
        }
        assert!(refined.report.accepted_intervals > 1);
        assert!(refined.report.rejected_intervals > 0);
        assert_eq!(
            refined.report.advance_calls,
            3 * (refined.report.accepted_intervals + refined.report.rejected_intervals)
        );
        assert!(refined.report.heat_trials >= refined.report.advance_calls);
        assert!(refined.report.maximum_accepted_error.normalized <= 1.0);
        let e = &refined.report.energy;
        close(
            e.gas_energy_change_j + e.thermal_storage_change_j,
            external,
            1e-6,
        );
        close(e.external_gas_heat_j, external, 1e-14);
        close(e.external_solid_heat_j, 0.0, 1e-6);
        close(
            e.combined_residual_j,
            0.0,
            adaptive().window_energy_tolerance_j,
        );
        assert!(e.absolute_residual_sum_j <= adaptive().window_energy_tolerance_j);
        assert!(e.gas_storage_defect_j <= adaptive().window_energy_tolerance_j);
        assert_eq!(refined.endpoint.system.time_s().to_bits(), 1.0f64.to_bits());
        assert_eq!(
            refined.endpoint.thermal.time_s().to_bits(),
            1.0f64.to_bits()
        );
        assert_eq!(system.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(observer.0.is_empty());
    });
}

#[derive(Clone, Default)]
struct RefuseAfter {
    time_s: f64,
    source_calls: usize,
}
impl DepositionObserver for RefuseAfter {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        self.source_calls += 1;
        if thermal.time_s() >= self.time_s {
            return Err(DepositionError::Cancelled);
        }
        Ok(BTreeMap::new())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn a_late_half_step_refusal_and_global_probe_cap_publish_no_partial_stroke() {
    context(|cx| {
        let fixture = thermally_constant_fixture();
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let gas = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.0, 1e24);
        let initial = gas.clone();
        let system = ThermalPorousGasSystem::try_new(
            solid,
            &thermal,
            vec![cavity_entry(Arc::new(cavity), gas)],
            vec![],
            cx,
        )
        .unwrap();
        let supports = affine(&nodes, 1.0);
        let old_cells = thermal.cells().clone();
        let refusal = RefuseAfter {
            time_s: 0.15,
            source_calls: 0,
        };
        let time = AdaptiveCavityHeatSolve {
            initial_step_s: 0.1,
            minimum_step_s: 0.01,
            maximum_step_s: 0.1,
            ..adaptive()
        };
        let refused = system.advance_adaptive(
            &thermal,
            &refusal,
            0.3,
            &supports,
            &[],
            [0.0; 3],
            None,
            cavity_settings(),
            controls(),
            thermal_config(),
            config(),
            CouplingMethod::RelaxedPicard,
            &[0.0],
            &[0.0],
            solve_heat(),
            time,
            cx,
        );
        match refused {
            Err(CavityHeatMarchError::Native {
                start_s,
                source:
                    CavityHeatError::Feedback(FiniteThermoplasticError::Observer(
                        CavityHeatObserverError::Observer(DepositionError::Cancelled),
                    )),
                ..
            }) => assert!(start_s >= 0.1),
            _ => panic!("expected a late native observer refusal"),
        }
        let native = system
            .advance(
                &thermal,
                &(),
                0.1,
                &supports,
                &[],
                [0.0; 3],
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                &[0.0],
                &[0.0],
                solve_heat(),
                cx,
            )
            .unwrap();
        let budget = system.advance_adaptive(
            &thermal,
            &(),
            0.3,
            &supports,
            &[],
            [0.0; 3],
            None,
            cavity_settings(),
            controls(),
            thermal_config(),
            config(),
            CouplingMethod::RelaxedPicard,
            &[0.0],
            &[0.0],
            solve_heat(),
            AdaptiveCavityHeatSolve {
                maximum_heat_trials: native.trials,
                ..time
            },
            cx,
        );
        assert!(
            matches!(budget, Err(CavityHeatMarchError::Budget { heat_trials, .. })
            if heat_trials == native.trials)
        );
        assert_eq!(system.time_s(), 0.0);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(system.cavities()[0].gas, initial);
        assert_eq!(refusal.source_calls, 0);
    });
}

#[test]
fn adaptive_compression_retains_plastic_volume_velocity_and_each_native_energy_account() {
    context(|cx| {
        let fixture = thermally_constant_fixture();
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx)
            .with_inertia(&vec![[0.0; 3]; nodes.len()], cx)
            .unwrap();
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let (gas, _) = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.001, 1e24)
            .at_temperature(360.0)
            .unwrap();
        let initial_volume = gas.volume_m3();
        let system = ThermalPorousGasSystem::try_new(
            solid,
            &thermal,
            vec![cavity_entry(Arc::new(cavity), gas)],
            links(&thermal, 0.005),
            cx,
        )
        .unwrap();
        let supports = affine(&nodes, 0.94);
        let next = system
            .advance_adaptive(
                &thermal,
                &(),
                1.0,
                &supports,
                &[],
                [0.0; 3],
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                &[0.005],
                &[0.0],
                solve_heat(),
                AdaptiveCavityHeatSolve {
                    initial_step_s: 0.25,
                    maximum_step_s: 0.25,
                    absolute_temperature_tolerance_k: 0.05,
                    absolute_specific_enthalpy_tolerance_j_kg: 0.05 * CP,
                    absolute_plastic_state_tolerance: 1e-4,
                    absolute_energy_tolerance_j: 1e-4,
                    relative_tolerance: 0.01,
                    ..adaptive()
                },
                cx,
            )
            .unwrap();
        let solid = next.endpoint.system.solid().mechanics();
        assert!(solid.element_states().iter().all(|x| x.alpha > 0.0));
        assert!(solid.element_plastic_jacobians().iter().all(|x| *x < 1.0));
        for &(dof, target) in &supports {
            close(solid.displacement_m().as_flattened()[dof], target, 1e-14);
            close(
                solid.velocity_mps().unwrap().as_flattened()[dof],
                target,
                1e-10,
            );
        }
        close(
            next.endpoint.system.cavities()[0].gas.volume_m3(),
            initial_volume * 0.94f64.powi(3),
            1e-18,
        );
        let e = &next.report.energy;
        assert!(e.plastic_dissipation_j > 0.0);
        assert!(e.kinetic_change_j > 0.0);
        assert!(e.kinetic_numerical_loss_j >= 0.0);
        assert!(e.gas_to_solid_heat_j > 0.0);
        close(e.external_gas_heat_j, 0.005, 1e-14);
        // Independent recombination includes numerical work explicitly; it is
        // not used as a physical source for either gas or original wall calorics.
        close(
            e.gas_energy_change_j
                + e.thermal_storage_change_j
                + e.solid_internal_energy_change_j
                + e.kinetic_change_j
                + e.kinetic_numerical_loss_j
                + e.mechanical_work_defect_j
                - e.external_work_j
                - e.external_gas_heat_j
                - e.external_solid_heat_j
                - e.incoming_thermal_energy_j,
            0.0,
            adaptive().window_energy_tolerance_j,
        );
        assert_eq!(system.time_s(), 0.0);
        assert!(
            system
                .solid()
                .mechanics()
                .element_states()
                .iter()
                .all(|x| x.alpha == 0.0)
        );
    });
}
