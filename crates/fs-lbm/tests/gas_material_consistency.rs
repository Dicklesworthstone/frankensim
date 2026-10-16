//! Frozen material snapshots drive physical flow, heat and acoustic properties
//! together. This is not a compressible, reacting or evolving-EOS LBM claim.
use fs_lbm::material::{MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::{GasSpec, GasState};

#[test]
fn transport_planning_refuses_a_gas_that_other_physics_would_reject() {
    let admitted = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap();
    let mutations: [fn(&mut GasState); 5] = [
        |gas| gas.pressure *= 2.0,
        |gas| gas.sound_speed *= 2.0,
        |gas| gas.prandtl *= 2.0,
        |gas| gas.characteristic_impedance *= 2.0,
        |gas| gas.specific_heat_cp *= 2.0,
    ];
    for mutate in mutations {
        let mut gas = admitted;
        mutate(&mut gas);
        assert!(gas.validate_coherence().is_err());
        assert!(gas.transport_state().is_err());
        assert!(plan_scaling_from_material(&gas, 0.001, 0.00002, 0.1).is_err());
    }
}

#[test]
fn changing_only_gas_data_changes_heating_acceleration_diffusion_and_sound() {
    let dx = 0.001;
    let dt = 0.00002;
    let power_density = 1e7;
    let mut responses = Vec::new();
    // Steam at 450 K avoids pretending that room-temperature, one-atmosphere
    // pure water vapor is a stable phase. No phase equilibrium is inferred.
    for spec in [GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()] {
        let gas = GasState::try_new(&spec, 450.0, 100000.0).unwrap();
        let transport = gas.transport_state().unwrap();
        let plan = plan_scaling_from_material(&gas, dx, dt, 0.1).unwrap();
        let mut slab = plan.thermal_slab3(MaterialThermalSlab3 {
            nx: 3, ny: 2, nz: 8,
            bottom_temperature_k: 460.0, top_temperature_k: 440.0,
            gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 60,
        }).unwrap();
        let n = slab.lattice.gpop.len();
        let mut power = vec![0.0; n];
        power[6..n - 6].fill(power_density);
        let before = slab.temperature_si(slab.lattice.temperature(1, 1, 4)).unwrap();
        let report = slab.step_with_power_density(&power, &vec![[1.0, 0.0, 0.0]; n],
            HeatSourceLimits { max_cells: 512, temperature_range_k: [400.0, 500.0],
                balance_tolerance_j: 1e-9 }).unwrap();
        let rise = slab.temperature_si(slab.lattice.temperature(1, 1, 4)).unwrap() - before;
        let expected = power_density * dt / (gas.density * gas.specific_heat_cp);
        assert!((rise - expected).abs() < 1e-9);
        assert!((report.source_j - 48.0 * power_density * dx.powi(3) * dt).abs() < 1e-12);
        assert!(report.defect_j.abs() < 1e-9);
        let cell = slab.lattice.grid.idx(1, 1, 4);
        let speed = plan.velocity_si(slab.lattice.grid.moments3(cell).1[0]).unwrap();
        assert!(speed > 0.0);
        assert!((gas.sound_speed.powi(2) / (gas.gamma * gas.pressure / gas.density) - 1.0).abs() < 1e-12);
        responses.push((gas, transport, rise, speed, report.source_j));
    }
    let (air, a, air_rise, air_speed, air_energy) = responses[0];
    let (steam, b, steam_rise, steam_speed, steam_energy) = responses[1];
    assert_ne!(a.kinematic_viscosity_m2_s(), b.kinematic_viscosity_m2_s());
    assert_ne!(a.thermal_diffusivity_m2_s(), b.thermal_diffusivity_m2_s());
    assert_ne!(air.sound_speed, steam.sound_speed);
    assert_ne!(air_rise, steam_rise);
    assert!((air_speed / steam_speed - steam.density / air.density).abs() < 1e-7);
    assert!((air_energy - steam_energy).abs() < 1e-12,
        "material changes response, not authored energy");
}
