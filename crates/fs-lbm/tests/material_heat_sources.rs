use fs_lbm::material::{MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};

fn limits() -> HeatSourceLimits {
    HeatSourceLimits { max_cells: 512, temperature_range_k: [280.0, 340.0], balance_tolerance_j: 1e-9 }
}
fn input2() -> MaterialThermalSlab {
    MaterialThermalSlab { nx: 3, ny: 8, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 30 }
}
fn input3() -> MaterialThermalSlab3 {
    MaterialThermalSlab3 { nx: 3, ny: 2, nz: 8, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 60 }
}
fn power(n: usize, plane: usize, value: f64) -> Vec<f64> {
    let mut p = vec![0.0; n]; p[plane..n-plane].fill(value); p
}

#[test]
fn material_heat_capacity_changes_temperature_not_the_authored_energy() {
    let mut rises = Vec::new();
    for cp in [2000.0, 4000.0] {
        let material = FluidTransportState::try_new(1000.0, cp, 1.0, 0.1).unwrap();
        let plan = plan_scaling_from_material(&material, 0.001, 0.001, 0.01).unwrap();
        let mut slab = plan.thermal_slab(input2()).unwrap(); let n = slab.lattice.gpop.len();
        let before = slab.temperature_si(slab.lattice.temperature(1, 4)).unwrap();
        let report = slab.step_with_power_density(&power(n, 3, 4e7), &vec![[1.0, 0.0]; n], 0.003, limits()).unwrap();
        let rise = slab.temperature_si(slab.lattice.temperature(1, 4)).unwrap() - before;
        assert!((rise - 4e7*0.001/(1000.0*cp)).abs() < 1e-10);
        assert!((report.source_j - 4e7*0.001*0.001*0.003*0.001*24.0).abs() < 1e-12);
        assert!(report.defect_j.abs() < 1e-12); rises.push(rise);
    }
    assert!((rises[0]/rises[1]-2.0).abs() < 1e-8);
}

#[test]
fn gas_pressure_data_changes_heating_and_acceleration_together() {
    let spec = GasSpec::dry_air_ussa1976(); let mut rises = Vec::new(); let mut speeds = Vec::new();
    for pressure in [100000.0, 200000.0] {
        let gas = GasState::try_new(&spec, 300.0, pressure).unwrap();
        let transport = gas.transport_state().unwrap();
        let plan = plan_scaling_from_material(&gas, 0.001, 0.00002, 0.1).unwrap();
        let mut slab = plan.thermal_slab3(input3()).unwrap(); let n = slab.lattice.gpop.len();
        let before = slab.temperature_si(slab.lattice.temperature(1, 1, 4)).unwrap();
        let report = slab.step_with_power_density(&power(n, 6, 1e7), &vec![[1.0, 0.0, 0.0]; n], limits()).unwrap();
        let rise = slab.temperature_si(slab.lattice.temperature(1, 1, 4)).unwrap()-before;
        assert!((rise-1e7*0.00002/(transport.density_kg_m3()*transport.specific_heat_j_kg_k())).abs() < 1e-9);
        assert!((report.source_j - 48.0*1e7*1e-9*0.00002).abs() < 1e-12);
        let i = slab.lattice.grid.idx(1, 1, 4);
        speeds.push(plan.velocity_si(slab.lattice.grid.moments3(i).1[0]).unwrap()); rises.push(rise);
    }
    assert!((rises[0]/rises[1]-2.0).abs() < 1e-8);
    assert!((speeds[0]/speeds[1]-2.0).abs() < 1e-7);
}

#[test]
fn explicit_2d_depth_scales_joules_without_changing_evolution() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 1.0, 0.1).unwrap();
    let plan = plan_scaling_from_material(&material, 0.001, 0.001, 0.01).unwrap();
    let mut a = plan.thermal_slab(input2()).unwrap(); let mut b = plan.thermal_slab(input2()).unwrap();
    let n = a.lattice.gpop.len(); let p = power(n, 3, 4e7); let f = vec![[0.0; 2]; n];
    let ra = a.step_with_power_density(&p, &f, 0.002, limits()).unwrap();
    let rb = b.step_with_power_density(&p, &f, 0.004, limits()).unwrap();
    assert_eq!(a.lattice.gpop, b.lattice.gpop); assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
    assert_eq!(rb.source_j, 2.0*ra.source_j); assert_eq!(rb.stored_change_j, 2.0*ra.stored_change_j);
    assert_eq!(rb.bottom_wall_j, 2.0*ra.bottom_wall_j); assert_eq!(rb.top_wall_j, 2.0*ra.top_wall_j);
}

#[test]
fn positive_power_heats_with_reversed_wall_temperature_coordinates() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 1.0, 0.1).unwrap();
    let plan = plan_scaling_from_material(&material, 0.001, 0.001, 0.01).unwrap();
    let mut input = input3(); input.bottom_temperature_k = 300.0; input.top_temperature_k = 310.0;
    let mut slab = plan.thermal_slab3(input).unwrap(); let n = slab.lattice.gpop.len();
    let before = slab.temperature_si(slab.lattice.temperature(1, 1, 4)).unwrap();
    let report = slab.step_with_power_density(&power(n, 6, 4e7), &vec![[0.0; 3]; n], limits()).unwrap();
    assert!((slab.temperature_si(slab.lattice.temperature(1, 1, 4)).unwrap()-before-0.01).abs() < 1e-10);
    assert!(report.source_j > 0.0); assert!(report.stored_change_j > 0.0); assert!(report.defect_j.abs() < 1e-12);
}

#[test]
fn physical_limit_refusals_do_not_publish_partial_heat_or_force_updates() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 1.0, 0.1).unwrap();
    let plan = plan_scaling_from_material(&material, 0.001, 0.001, 0.01).unwrap();
    let mut slab = plan.thermal_slab(input2()).unwrap(); let n = slab.lattice.gpop.len();
    let before = (slab.lattice.grid.f.clone(), slab.lattice.gpop.clone(), slab.lattice.grid.fext.clone());
    let f = vec![[1.0, 0.0]; n]; let p = power(n, 3, 4e11);
    assert!(slab.step_with_power_density(&p, &f, 0.003, limits()).is_err());
    assert!(slab.step_with_power_density(&p, &f, 0.0, limits()).is_err());
    let mut bad = limits(); bad.temperature_range_k[0] = 0.0;
    assert!(slab.step_with_power_density(&p, &f, 0.003, bad).is_err());
    assert_eq!((&slab.lattice.grid.f, &slab.lattice.gpop, &slab.lattice.grid.fext), (&before.0, &before.1, &before.2));
}

#[test]
fn heat_uses_the_retained_snapshot_and_survives_conservative_transport_rebinding() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get()+1); FluidTransportState::try_new(1000.0, 4000.0, 1.0, 0.1)
        }
    }
    let original = Counted(std::cell::Cell::new(0));
    let changed = FluidTransportState::try_new(1000.0, 4000.0, 2.0, 0.2).unwrap();
    let old = plan_scaling_from_material(&original, 0.001, 0.001, 0.01).unwrap();
    let new = plan_scaling_from_material(&changed, 0.001, 0.001, 0.01).unwrap();
    let mut slab = old.thermal_slab3(input3()).unwrap(); let n = slab.lattice.gpop.len();
    let p = power(n, 6, 4e7); let f = vec![[0.0; 3]; n];
    let first = slab.step_with_power_density(&p, &f, limits()).unwrap();
    let heat_before = slab.lattice.gpop.clone();
    new.rebind_thermal_slab3(&mut slab, &old, 0.0, 0.0).unwrap();
    assert_eq!(slab.lattice.gpop, heat_before);
    let next = slab.step_with_power_density(&p, &f, limits()).unwrap();
    assert_eq!(next.source_j, first.source_j); assert_eq!(original.0.get(), 1);
    let off = slab.step_with_power_density(&vec![0.0; n], &f, limits()).unwrap();
    assert_eq!(off.source_j, 0.0);
}
