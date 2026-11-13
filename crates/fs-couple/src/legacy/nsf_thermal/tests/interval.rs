use super::*;

fn interval(steps: usize) -> NsfThermalIntervalControls {
    NsfThermalIntervalControls {
        fluid_steps: steps,
        maximum_fluid_cell_steps: 4 * steps * controls().maximum_endpoint_evaluations,
    }
}

// Independent zero-dimensional energy balance for this uniform resolved grid.
// Each fluid substep uses its original trapezoidal molecular boundary law;
// only the solid endpoint is held over the complete coupling interval.
fn uniform_fluid_endpoint(law: &NsfLaw, capacity: f64, wall_k: f64, steps: usize) -> f64 {
    let mut temperature = 300.0;
    for _ in 0..steps {
        let incoming = scalar_flux(law, temperature, wall_k, 1.2e-4);
        let (mut lo, mut hi) = (temperature, wall_k);
        for _ in 0..70 {
            let endpoint = f64::midpoint(lo, hi);
            let joules = 0.5
                * law.dt_s
                * (2.0 * PLANE_AREA)
                * (incoming + scalar_flux(law, endpoint, wall_k, 1.2e-4));
            if capacity * (endpoint - temperature) > joules {
                hi = endpoint;
            } else {
                lo = endpoint;
            }
        }
        temperature = f64::midpoint(lo, hi);
    }
    temperature
}

#[test]
fn g1_multirate_wall_closes_each_native_transfer_and_two_finite_capacities() {
    with_cx(|cx, _| {
        let law = Law { datum: 1e8 };
        let conductivity = K(0.1);
        let mut gas = fluid();
        let gas_capacity = gas.total_mass_kg() * gas.law().gas_constant / (gas.law().gamma - 1.0);
        let wall_capacity = 2.0 * gas_capacity;
        let mass = wall_capacity / CP;
        let mut wall = solid(&law, &[310.0], mass, false, false);
        let old_h = wall.cells()[&CellId(0)].specific_h_j_kg;
        let mapping = contacts(
            &law,
            NsfSolidLayer::Conducting {
                half_distance_m: 1e-5,
                contact_resistance_m2_k_w: 2e-5,
                conductivity: &conductivity,
            },
            false,
        );
        let steps = 4;
        let (mut low, mut high) = (0.0, 10.0 / (1.0 / gas_capacity + 1.0 / wall_capacity));
        for _ in 0..80 {
            let heat = f64::midpoint(low, high);
            let wall_k = 310.0 - heat / wall_capacity;
            let endpoint = uniform_fluid_endpoint(gas.law(), gas_capacity, wall_k, steps);
            if heat > gas_capacity * (endpoint - 300.0) {
                high = heat;
            } else {
                low = heat;
            }
        }
        let expected_heat = f64::midpoint(low, high);
        let report = advance_nsf_solid_interval(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            interval(steps),
            &[0.0],
            unforced,
        )
        .unwrap();
        near(report.coupled_heat_into_fluid_j, expected_heat, 1e-17, 5e-6);
        near(
            gas.cell([0, 0, 0]).2,
            300.0 + expected_heat / gas_capacity,
            3e-6,
            0.0,
        );
        near(
            wall.temperature_k(CellId(0)).unwrap(),
            310.0 - expected_heat / wall_capacity,
            3e-6,
            0.0,
        );
        near(
            mass * (wall.cells()[&CellId(0)].specific_h_j_kg - old_h),
            -report.coupled_heat_into_fluid_j,
            controls().thermal_tolerance_j,
            0.0,
        );
        near(
            report
                .exchanges
                .iter()
                .map(|face| face.face.heat_into_fluid_j)
                .sum(),
            report.coupled_heat_into_fluid_j,
            1e-24,
            1e-14,
        );
        assert_eq!(report.fluid_steps, steps);
        assert_eq!(
            report.fluid_cell_steps,
            4 * steps * report.endpoint_evaluations
        );
        assert_eq!(report.duration_s, steps as f64 * gas.law().dt_s);
        assert_eq!(gas.steps(), steps);
        assert_eq!(wall.time_s(), steps as f64 * gas.law().dt_s);
        assert_eq!(report.thermal.substeps, 1);
        assert!(report.total_energy_residual_j.abs() <= controls().total_tolerance_j);

        let continued = advance_nsf_solid_interval(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            interval(2),
            &[0.0],
            unforced,
        )
        .unwrap();
        assert_eq!(gas.steps(), 6);
        assert_eq!(wall.time_s(), 6.0 * gas.law().dt_s);
        assert_eq!(continued.duration_s, 2.0 * gas.law().dt_s);
        assert!(continued.coupled_heat_into_fluid_j > 0.0);
    });
}

#[test]
fn g3_one_substep_is_the_original_entry_with_identical_physical_accounts() {
    with_cx(|cx, _| {
        let law = Law { datum: 0.0 };
        let mapping = contacts(&law, NsfSolidLayer::Isothermal, true);
        let mut old_gas = fluid();
        let mut interval_gas = fluid();
        let mut old_wall = solid(&law, &[310.0, 290.0], 1e-11, true, true);
        let mut interval_wall = old_wall.clone();
        let old = advance_nsf_solid(
            cx,
            &mut old_gas,
            &mut old_wall,
            &mapping,
            controls(),
            &[0.0; 2],
            unforced,
        )
        .unwrap();
        let many = advance_nsf_solid_interval(
            cx,
            &mut interval_gas,
            &mut interval_wall,
            &mapping,
            controls(),
            interval(1),
            &[0.0; 2],
            unforced,
        )
        .unwrap();
        assert_eq!(old_gas.state_hash(), interval_gas.state_hash());
        assert_eq!(old_wall.cells(), interval_wall.cells());
        assert_eq!(old_wall.time_s(), interval_wall.time_s());
        assert_eq!(old.exchanges, many.exchanges);
        assert_eq!(
            old.coupled_heat_into_fluid_j,
            many.coupled_heat_into_fluid_j
        );
        assert_eq!(old.solid_external_heat_j, many.solid_external_heat_j);
        assert_eq!(old.total_energy_residual_j, many.total_energy_residual_j);
        assert_eq!(old.endpoint_evaluations, many.endpoint_evaluations);
    });
}

#[test]
fn g1_dynamic_channel_exchanges_heat_with_the_same_finite_wall_over_native_steps() {
    use fs_lbm::nsf3::NsfChannelDynamicThermalSgsConfig;
    use fs_lbm::sgs3::SgsUse;
    use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalarConfig, ChannelDynamicSgsConfig};
    with_cx(|cx, _| {
        // Same dimensions, scaling and exploration bounds as the native
        // dynamic-channel fixture; this checks owner composition, not LES
        // calibration or a physical turbulent boundary-layer benchmark.
        let shape = [5, 4, 5];
        let dynamic = NsfChannelDynamicThermalSgsConfig {
            momentum: ChannelDynamicSgsConfig {
                maximum_coefficient_squared: 16.0,
                maximum_eddy_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
            heat: ChannelDynamicScalarConfig {
                maximum_coefficient: 16.0,
                maximum_diffusivity_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
        };
        let mut gas = NsfGrid3::<D3q27>::new(
            NsfSpec {
                dims: shape,
                periodic: [true, true, false],
                gas: GasSpec::dry_air_ussa1976(),
                temperature_k: 293.15,
                pressure_pa: 101325.0,
                dx_m: 5e-5,
                temperature_window_k: [273.15, 333.15],
                maximum_mach: 0.2,
                maximum_temperature_ratio: 0.2,
                bulk_to_shear_ratio: 0.73,
                bulk_ratio_source: "native dynamic-channel finite wall fixture",
                maximum_cells: 100,
                wall: WallThermal::Adiabatic,
                workers: 1,
                transport_points: 9,
            },
            |_| false,
            |p| {
                let y = (p[2] as f64 + 0.5) / shape[2] as f64;
                let profile = 4.0 * y * (1.0 - y);
                let wave = core::f64::consts::TAU * (p[0] as f64 + 0.5) / shape[0] as f64;
                (
                    1.0,
                    [0.001 * profile, 0.0001 * profile * wave.sin(), 0.0],
                    293.15 + 0.1 * (1.0 - y),
                )
            },
        )
        .unwrap()
        .with_channel_dynamic_thermal_sgs(dynamic, 1 << 20)
        .unwrap();
        let retained_workspace = gas.dynamic_sgs_workspace_bytes();
        let law = Law { datum: 0.0 };
        let mut wall = solid(&law, &[310.0], 1e-11, false, false);
        let old_h = wall.cells()[&CellId(0)].specific_h_j_kg;
        let mapping = [NsfSolidContact {
            fluid_cell: [0, 0, 0],
            axis: 2,
            high: false,
            solid_cell: CellId(0),
            caloric: &law,
            layer: NsfSolidLayer::Isothermal,
        }];
        let mut control = controls();
        control.maximum_fluid_cells = 100;
        control.fluid_tolerance_j = 1e-14;
        control.total_tolerance_j = 2e-14;
        let report = advance_nsf_solid_interval(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            control,
            NsfThermalIntervalControls {
                fluid_steps: 2,
                maximum_fluid_cell_steps: 200 * control.maximum_endpoint_evaluations,
            },
            &[0.0],
            unforced,
        )
        .unwrap();
        assert_eq!(gas.channel_dynamic_thermal_sgs_config(), Some(dynamic));
        assert_eq!(gas.dynamic_sgs_workspace_bytes(), retained_workspace);
        assert_eq!(gas.steps(), 2);
        assert_eq!(report.exchanges.len(), 1);
        assert!(report.coupled_heat_into_fluid_j > 0.0);
        assert!(wall.temperature_k(CellId(0)).unwrap() < 310.0);
        near(
            1e-11 * (wall.cells()[&CellId(0)].specific_h_j_kg - old_h),
            -report.coupled_heat_into_fluid_j,
            control.cell_tolerance_j,
            0.0,
        );
        near(
            report.fluid.wall_heat_j,
            report.coupled_heat_into_fluid_j,
            1e-22,
            1e-12,
        );
        assert_eq!(report.fluid_cell_steps, 200 * report.endpoint_evaluations);
        assert!(report.total_energy_residual_j.abs() <= control.total_tolerance_j);
    });
}

#[test]
fn g4_multirate_shared_work_and_late_source_refusal_preserve_complete_states() {
    with_cx(|cx, gate| {
        let law = Law { datum: 0.0 };
        let mut gas = fluid();
        let mut wall = solid(&law, &[310.0], 1e-11, false, false);
        let mapping = contacts(&law, NsfSolidLayer::Isothermal, false);
        let before = (
            gas.state_hash(),
            wall.cells().clone(),
            wall.time_s(),
            gas.steps(),
        );
        let limited = NsfThermalIntervalControls {
            fluid_steps: 4,
            maximum_fluid_cell_steps: 16,
        };
        let error = advance_nsf_solid_interval(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            limited,
            &[0.0],
            unforced,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            NsfThermalError::FluidWorkBudget {
                performed: 16,
                requested: 16,
                maximum: 16,
            }
        ));
        for bad in [
            NsfThermalIntervalControls {
                fluid_steps: 0,
                ..interval(4)
            },
            NsfThermalIntervalControls {
                fluid_steps: usize::MAX,
                maximum_fluid_cell_steps: usize::MAX,
            },
        ] {
            assert!(
                advance_nsf_solid_interval(
                    cx,
                    &mut gas,
                    &mut wall,
                    &mapping,
                    controls(),
                    bad,
                    &[0.0],
                    unforced
                )
                .is_err()
            );
        }
        assert_eq!(
            before,
            (
                gas.state_hash(),
                wall.cells().clone(),
                wall.time_s(),
                gas.steps()
            )
        );
        let queries = Cell::new(0);
        let error = advance_nsf_solid_interval(
            cx,
            &mut gas,
            &mut wall,
            &mapping,
            controls(),
            interval(4),
            &[0.0],
            |cx, candidate| {
                queries.set(queries.get() + 1);
                if queries.get() == 2 {
                    gate.request();
                }
                unforced(cx, candidate)
            },
        )
        .unwrap_err();
        assert!(matches!(error, NsfThermalError::Cancelled));
        assert_eq!(queries.get(), 2);
        assert_eq!(
            before,
            (
                gas.state_hash(),
                wall.cells().clone(),
                wall.time_s(),
                gas.steps()
            )
        );
    });
}
