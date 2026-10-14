//! Data-only material changes must drive both heat and momentum in 2-D and 3-D.
use fs_lbm::material::{
    MaterialLbmError, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};

fn boundary() -> MaterialThermalSlab3 {
    MaterialThermalSlab3 {
        nx: 4, ny: 3, nz: 6,
        bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 9.81, expansion_per_k: 2e-4, max_cells: 96,
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-11 * a.abs().max(b.abs()).max(1e-12), "{a} != {b}");
}

#[test]
fn one_material_preserves_transport_and_rayleigh_across_dimensions() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.0).unwrap();
    let b = boundary();
    let slab3 = scale.thermal_slab3(b).unwrap();
    let slab2 = scale.thermal_slab(MaterialThermalSlab {
        nx: b.nx, ny: b.nz,
        bottom_temperature_k: b.bottom_temperature_k,
        top_temperature_k: b.top_temperature_k,
        gravity_m_s2: b.gravity_m_s2, expansion_per_k: b.expansion_per_k,
        max_cells: b.nx * (b.nz + 2),
    }).unwrap();
    close(slab3.lattice.grid.tau[0], scale.flow().tau);
    close(slab3.lattice.diffusivity(), scale.thermal_diffusivity_lattice());
    close((slab3.lattice.grid.tau[0] - 0.5) / 3.0 / slab3.lattice.diffusivity(), material.prandtl());
    close(slab3.lattice.gbeta, slab2.lattice.gbeta);
    close(slab3.lattice.rayleigh(), material.rayleigh(
        b.gravity_m_s2, b.expansion_per_k, 10.0, b.nz as f64 * scale.cell_size_m(),
    ).unwrap());
    // Reusing the 2-D relaxation in D3Q7 would lose 25% of the diffusivity.
    assert_ne!(slab3.lattice.tau_g, scale.thermal_relaxation_time());
    close(slab3.temperature_si(0.25).unwrap(), slab2.temperature_si(0.25).unwrap());
}

#[test]
fn gas_pressure_changes_both_3d_diffusion_rates_without_solver_changes() {
    let spec = GasSpec::dry_air_ussa1976();
    let low = GasState::try_new(&spec, 300.0, 100000.0).unwrap();
    let high = GasState::try_new(&spec, 300.0, 200000.0).unwrap();
    let a = plan_scaling_from_material(&low, 1e-3, 1e-3, 0.0).unwrap();
    let b = plan_scaling_from_material(&high, 1e-3, 1e-3, 0.0).unwrap();
    let slab_a = a.thermal_slab3(boundary()).unwrap();
    let slab_b = b.thermal_slab3(boundary()).unwrap();
    close(b.density_si(1.0).unwrap(), 2.0 * a.density_si(1.0).unwrap());
    close(slab_b.lattice.grid.tau[0] - 0.5, (slab_a.lattice.grid.tau[0] - 0.5) / 2.0);
    close(slab_b.lattice.diffusivity(), slab_a.lattice.diffusivity() / 2.0);
    close(b.transport().prandtl(), a.transport().prandtl());
}

#[test]
fn actual_3d_solver_preserves_quiescent_conduction_and_mass() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.0).unwrap();
    let mut b = boundary();
    b.gravity_m_s2 = 0.0;
    let mut slab = scale.thermal_slab3(b).unwrap();
    let mass = slab.lattice.grid.total_mass();
    for _ in 0..40 { slab.lattice.step(); }
    close(slab.lattice.grid.total_mass(), mass);
    for z in 1..=b.nz {
        let expected = b.bottom_temperature_k - 10.0 * (z as f64 - 0.5) / b.nz as f64;
        close(slab.temperature_si(slab.lattice.temperature(1, 1, z)).unwrap(), expected);
    }
    close(slab.lattice.nusselt(), 1.0);
}

#[test]
fn signed_temperature_and_expansion_reverse_buoyancy_not_materials() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.0).unwrap();
    let mut b = boundary();
    let hot_bottom = scale.thermal_slab3(b).unwrap();
    std::mem::swap(&mut b.bottom_temperature_k, &mut b.top_temperature_k);
    let hot_top = scale.thermal_slab3(b).unwrap();
    close(hot_top.lattice.gbeta, -hot_bottom.lattice.gbeta);
    close(hot_top.temperature_si(0.0).unwrap(), 310.0);
    close(hot_top.temperature_si(1.0).unwrap(), 300.0);
    b.expansion_per_k = -b.expansion_per_k;
    close(scale.thermal_slab3(b).unwrap().lattice.gbeta, hot_bottom.lattice.gbeta);
}

#[test]
fn dimensions_and_budget_are_refused_before_raw_constructor_panics() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.0).unwrap();
    for (nx, ny, nz, budget) in [
        (0, 3, 6, 96), (4, 0, 6, 96), (4, 3, 1, 96), (4, 3, 6, 95),
        (4, 3, usize::MAX, usize::MAX), (usize::MAX, 3, 6, usize::MAX),
    ] {
        assert!(scale.flow_plates3(nx, ny, nz, budget).is_err());
        let mut b = boundary();
        b.nx = nx; b.ny = ny; b.nz = nz; b.max_cells = budget;
        assert!(scale.thermal_slab3(b).is_err());
    }
    let grid = scale.flow_plates3(4, 3, 6, 96).unwrap();
    assert_eq!(grid.f.len(), 96);
    assert!(grid.tau.iter().all(|&tau| tau == scale.flow().tau));
}

#[test]
fn three_dimensional_boundaries_reject_nonphysical_values() {
    let material = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.0).unwrap();
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY, 300.0] {
        let mut b = boundary(); b.bottom_temperature_k = value;
        assert!(scale.thermal_slab3(b).is_err());
    }
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        let mut b = boundary(); b.gravity_m_s2 = value;
        assert!(scale.thermal_slab3(b).is_err());
    }
}

#[test]
fn all_builders_share_one_snapshot_and_propagate_material_errors() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.0.set(self.0.get() + 1);
            FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001)
        }
    }
    let material = Counted(std::cell::Cell::new(0));
    let scale = plan_scaling_from_material(&material, 1e-4, 1e-3, 0.0).unwrap();
    let _ = scale.flow_grid(4, 4, 16).unwrap();
    let _ = scale.flow_plates3(4, 3, 6, 96).unwrap();
    let _ = scale.thermal_slab3(boundary()).unwrap();
    assert_eq!(material.0.get(), 1);

    struct Refused;
    impl FluidTransport for Refused {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            Err(MaterialError::State { what: "outside liquid calibration".into() })
        }
    }
    assert!(matches!(
        plan_scaling_from_material(&Refused, 1e-4, 1e-3, 0.0),
        Err(MaterialLbmError::Transport(MaterialError::State { .. }))
    ));
}
