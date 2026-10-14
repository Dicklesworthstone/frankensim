//! Material-to-solver regressions use existing numerical kernels. Synthetic
//! coefficients are test inputs, not measured material cards or validation data.
use fs_lbm::material::{MaterialThermalSlab, plan_scaling_from_material};
use fs_lbm::{Cell, CS2};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{ConductivityModel, GasSpec, GasState};
use fs_qty::semantic::{Composition, CompositionBasis};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "actual={actual}, expected={expected}");
}
fn slab_input() -> MaterialThermalSlab {
    MaterialThermalSlab { nx: 4, ny: 8, bottom_temperature_k: 310.0,
        top_temperature_k: 300.0, gravity_m_s2: 9.81, expansion_per_k: 1.0 / 305.0,
        max_cells: 40 }
}

#[test]
fn pressure_changes_both_populations_without_inventing_a_sound_speed_change() {
    let a = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let b = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200_000.0).unwrap();
    let pa = plan_scaling_from_material(&a, 0.001, 0.001, 0.01).unwrap();
    let pb = plan_scaling_from_material(&b, 0.001, 0.001, 0.01).unwrap();
    near(pb.flow().viscosity, pa.flow().viscosity / 2.0, 1e-14);
    near(pb.thermal_diffusivity_lattice(), pa.thermal_diffusivity_lattice() / 2.0, 1e-14);
    near(pb.density_si(1.0).unwrap(), 2.0 * pa.density_si(1.0).unwrap(), 1e-14);
    assert_eq!(a.sound_speed, b.sound_speed);
    assert_eq!(pa.pressure_perturbation_si(1.0).unwrap(), 0.0);
    near(pb.pressure_perturbation_si(1.01).unwrap(), 2.0 * pa.pressure_perturbation_si(1.01).unwrap(), 1e-14);
    assert_eq!(pb.flow_grid(4, 3, 12).unwrap().tau, vec![pb.flow().tau; 12]);
}

#[test]
fn gas_composition_reaches_actual_thermal_evolution() {
    let air_spec = GasSpec::dry_air_ussa1976();
    let light_spec = GasSpec { molar_mass: 0.004, gamma: 5.0 / 3.0,
        sutherland_beta: 1e-6, sutherland_s: 80.0, conductivity: ConductivityModel::Eucken };
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap();
    let air = GasState::try_new(&air_spec, 300.0, 100_000.0).unwrap();
    let mixture = GasState::try_new_mixture(&[air_spec, light_spec], &composition, 300.0, 100_000.0, None).unwrap();
    let mut results = Vec::new();
    for gas in [&air, &mixture] {
        let plan = plan_scaling_from_material(gas, 0.001, 0.001, 0.01).unwrap();
        let mut slab = plan.thermal_slab(MaterialThermalSlab { gravity_m_s2: 0.0, ..slab_input() }).unwrap();
        near(plan.flow().viscosity / plan.thermal_diffusivity_lattice(), gas.prandtl, 2e-14);
        // Perturb real temperature populations, not a diagnostic coefficient.
        let center = slab.lattice.grid.idx(0, 4);
        slab.lattice.gpop[center][0] += 0.1;
        for _ in 0..20 { slab.lattice.step(); }
        let theta = slab.lattice.temperature(0, 4);
        assert!(theta.is_finite());
        assert!((300.0..310.0).contains(&slab.temperature_si(theta).unwrap()));
        results.push(theta);
    }
    assert!((results[0] - results[1]).abs() > 1e-6);
}

#[test]
fn physical_rayleigh_and_temperature_mapping_keep_both_signs() {
    let fluid = FluidTransportState::try_new(1.2, 1000.0, 0.026, 1.8e-5).unwrap();
    let plan = plan_scaling_from_material(&fluid, 0.001, 0.001, 0.01).unwrap();
    for (bottom, top, beta) in [(310.0, 300.0, 0.003), (300.0, 310.0, 0.003),
        (310.0, 300.0, -0.003), (300.0, 310.0, -0.003)] {
        let input = MaterialThermalSlab { bottom_temperature_k: bottom, top_temperature_k: top,
            expansion_per_k: beta, ..slab_input() };
        let slab = plan.thermal_slab(input).unwrap();
        let expected = 9.81 * beta * (bottom - top) * 0.008_f64.powi(3)
            / (fluid.kinematic_viscosity_m2_s() * fluid.thermal_diffusivity_m2_s());
        near(slab.lattice.rayleigh(), expected, 2e-14);
        near(slab.temperature_si(1.0).unwrap(), bottom, 1e-14);
        near(slab.temperature_si(0.0).unwrap(), top, 1e-14);
        near(slab.temperature_si(0.5).unwrap(), 305.0, 1e-14);
    }
}

#[test]
fn changing_viscosity_changes_a_real_channel_at_fixed_physical_geometry_and_force() {
    let mut velocities = Vec::new();
    for viscosity in [0.001, 0.002] {
        // Shared constructor order is rho, cp, conductivity, dynamic viscosity.
        let fluid = FluidTransportState::try_new(1000.0, 1000.0, 1.0, viscosity).unwrap();
        assert_eq!(fluid.dynamic_viscosity_pa_s(), viscosity);
        assert_eq!(fluid.thermal_conductivity_w_m_k(), 1.0);
        let plan = plan_scaling_from_material(&fluid, 1e-4, 1e-3, 1e-3).unwrap();
        let mut grid = plan.flow_grid(4, 10, 40).unwrap();
        grid.periodic_y = false;
        let force = plan.acceleration_lattice(0.001).unwrap();
        for y in 0..10 {
            for x in 0..4 {
                let i = grid.idx(x, y);
                if y == 0 || y == 9 { grid.flags[i] = Cell::Wall; }
                else { grid.fext[i] = [force, 0.0]; }
            }
        }
        let mut scratch = Vec::new();
        for _ in 0..5000 { grid.step(&mut scratch); }
        let velocity = plan.velocity_si(grid.moments(grid.idx(0, 4)).u[0]).unwrap();
        let z = 3.5e-4;
        let exact = 0.001 * z * (8e-4 - z) / (2.0 * fluid.kinematic_viscosity_m2_s());
        near(velocity, exact, 0.025); // finite-grid BGK wall/forcing error
        velocities.push(velocity);
    }
    near(velocities[0] / velocities[1], 2.0, 0.025);
}

#[test]
fn physical_units_round_trip_and_prandtl_survive_resolution_changes() {
    let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    for (dx, dt) in [(0.001, 0.001), (0.0001, 0.0002)] {
        let plan = plan_scaling_from_material(&fluid, dx, dt, 0.01).unwrap();
        near(CS2 * (plan.flow().tau - 0.5) * dx * dx / dt, fluid.kinematic_viscosity_m2_s(), 1e-12);
        near(CS2 * (plan.thermal_relaxation_time() - 0.5) * dx * dx / dt, fluid.thermal_diffusivity_m2_s(), 1e-12);
        near(plan.flow().viscosity / plan.thermal_diffusivity_lattice(), fluid.prandtl(), 1e-14);
        near(plan.velocity_si(plan.velocity_lattice(-0.007).unwrap()).unwrap(), -0.007, 1e-14);
    }
}

#[test]
fn source_is_sampled_once_and_retained_while_multiple_solvers_are_built() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.0.set(self.0.get() + 1);
            FluidTransportState::try_new(1.0, 1000.0, 0.01, 1e-5)
        }
    }
    let source = Counted(std::cell::Cell::new(0));
    let erased: &dyn FluidTransport = &source;
    let plan = plan_scaling_from_material(erased, 0.001, 0.001, 0.01).unwrap();
    plan.flow_grid(4, 4, 16).unwrap();
    plan.thermal_slab(slab_input()).unwrap();
    assert_eq!(source.0.get(), 1);
    assert!(core::ptr::eq(plan.material(), erased));
}

#[test]
fn invalid_domains_and_allocations_refuse_before_solver_construction() {
    let fluid = FluidTransportState::try_new(1.0, 1000.0, 0.01, 1e-5).unwrap();
    for (dx, dt, speed) in [(0.0, 0.001, 0.01), (0.001, -1.0, 0.01),
        (0.001, 0.1, 100.0), (0.001, 0.001, -0.1), (f64::NAN, 0.001, 0.0),
        (1e100, 1e-100, 0.0)] {
        assert!(plan_scaling_from_material(&fluid, dx, dt, speed).is_err());
    }
    let plan = plan_scaling_from_material(&fluid, 0.001, 0.001, 0.0).unwrap();
    assert!(plan.flow_grid(usize::MAX, 2, usize::MAX).is_err());
    assert!(plan.flow_grid(4, 4, 15).is_err());
    assert!(plan.flow_grid(0, 4, 16).is_err());
    assert!(plan.thermal_slab(MaterialThermalSlab { ny: usize::MAX, ..slab_input() }).is_err());
    assert!(plan.thermal_slab(MaterialThermalSlab { max_cells: 39, ..slab_input() }).is_err());
    assert!(plan.thermal_slab(MaterialThermalSlab { bottom_temperature_k: 300.0, ..slab_input() }).is_err());
    assert!(plan.thermal_slab(MaterialThermalSlab { bottom_temperature_k: 0.0, ..slab_input() }).is_err());
    assert!(plan.thermal_slab(MaterialThermalSlab { expansion_per_k: f64::NAN, ..slab_input() }).is_err());
}
