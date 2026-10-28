use super::*;
use crate::conjugate::{FaceContactResistance, PressureSolver, TimeScheme, Voxel,
    simple_unsteady};
use crate::d3q19::Face3;

fn fluid() -> FluidProperties {
    FluidProperties { density_kg_m3: 1.0, specific_heat_j_kg_k: 1000.0,
        conductivity_w_m_k: 0.04, kinematic_viscosity_m2_s: 1e-5 }
}
fn controls() -> WallHeatConfig {
    WallHeatConfig { turbulent_prandtl: 0.9, relative_tolerance: 1e-9, max_refinements: 14 }
}
fn simple() -> SimpleConfig {
    let mut c = SimpleConfig::new([FvBoundary::wall(); 6]);
    c.turbulence = Turbulence::Lvel;
    c.pressure_solver = PressureSolver::IluBicgstab;
    c
}
fn boundary(cell: usize, face: usize, condition: ThermalFace) -> BoundaryWallHeat {
    BoundaryWallHeat { wall: WallHeatFace { fluid_cell: cell, face: Face3::ALL[face] },
        condition, contact_m2_k_w: 0.0 }
}
fn time(dt: f64, steps: usize) -> UnsteadyConfig {
    let mut t = UnsteadyConfig::new(dt, steps);
    t.scheme = TimeScheme::BackwardEuler;
    t.inner_iterations = 400;
    t.inner_tolerance = 1e-8;
    t.average_from_step = 1;
    t
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8 * b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn rest_cooling_matches_native_storage_and_does_not_accumulate_wall_sinks() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let walls = WallHeatMarchConfig { internal_walls: vec![],
        boundary_walls: vec![boundary(0, 0, ThermalFace::Temperature(300.0))], transport: controls() };
    let clock = time(0.02, 8);
    let result = march_conjugate_with_wall_heat(&domain, &fluid(), &[], &setup, &simple(),
        &clock, &walls, None, &[310.0], |_| 1.0, |_| 1.0,
        &EnergyConfig::default(), None, &gate).unwrap();
    let capacity = fluid().volumetric_heat_capacity() * domain.volume(0);
    let g = domain.face_area(0, 0) / (0.005 / fluid().conductivity_w_m_k);
    let expected = 300.0 + 10.0 * (1.0 + clock.time_step_s * g / capacity).powi(-(clock.steps as i32));
    close(result.march.temperature[0], expected);
    close(result.boundary_outflow_j, capacity * (310.0 - expected));
    assert_eq!(result.wall_records.len(), clock.steps);
    for (wall, energy) in result.wall_records.iter().zip(&result.march.energy_records) {
        assert_eq!(wall.time_s, energy.time_s);
        assert_eq!(wall.maximum_y_plus, 0.0);
        assert_eq!(wall.quadrature_evaluations, 0);
        close(energy.stored_energy_change_j + clock.time_step_s * wall.boundary_outflow_w, 0.0);
        close(energy.closure_j, 0.0);
    }
}

#[test]
fn prescribed_flux_is_not_scaled_or_accumulated_with_authored_volume_power() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![7.0];
    let original = setup.clone();
    let walls = WallHeatMarchConfig { internal_walls: vec![],
        boundary_walls: vec![boundary(0, 0, ThermalFace::HeatFlux(12.0))], transport: controls() };
    let clock = time(0.02, 8);
    let result = march_conjugate_with_wall_heat(&domain, &fluid(), &[], &setup, &simple(),
        &clock, &walls, None, &[300.0], |_| 1.0, |_| 0.0,
        &EnergyConfig::default(), None, &gate).unwrap();
    let watts = 12.0 * domain.face_area(0, 0);
    let energy = watts * clock.time_step_s * clock.steps as f64;
    close(result.march.temperature[0], 300.0 + energy / (fluid().volumetric_heat_capacity() * domain.volume(0)));
    close(result.boundary_outflow_j, -energy);
    for (wall, native) in result.wall_records.iter().zip(&result.march.energy_records) {
        close(wall.prescribed_flux_input_w, watts);
        close(native.source_j, watts * clock.time_step_s);
        close(native.stored_energy_change_j, native.source_j);
        close(native.closure_j, 0.0);
    }
    assert_eq!(setup, original);
}

#[test]
fn ramped_native_flow_refreshes_wall_transport_without_changing_momentum() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(6, 4, 2, 0.01).unwrap();
    let mut flow = simple();
    flow.faces[0] = FvBoundary::Inlet { velocity: [2.0, 0.0, 0.0] };
    flow.faces[1] = FvBoundary::Outlet;
    flow.faces[4] = FvBoundary::Symmetry;
    flow.faces[5] = FvBoundary::Symmetry;
    let clock = time(0.01, 4);
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.faces[0] = ThermalFace::Inflow { temperature: 300.0 };
    thermal.faces[1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    let mut selection = Vec::new();
    for z in 0..2 { for x in 0..6 {
        selection.push(boundary(domain.index(x, 0, z), 2, ThermalFace::Temperature(310.0)));
    }}
    let walls = WallHeatMarchConfig { internal_walls: vec![], boundary_walls: selection, transport: controls() };
    let ramp = |t: f64| t / 0.04;
    let result = march_conjugate_with_wall_heat(&domain, &fluid(), &[], &thermal, &flow,
        &clock, &walls, None, &vec![300.0; domain.cell_count()], ramp, |_| 1.0,
        &EnergyConfig::default(), None, &gate).unwrap();
    let reference = simple_unsteady(&domain, &fluid(), &flow, &clock, ramp, &gate).unwrap();
    assert_eq!(result.march.flow, reference, "thermal observations cannot advance or alter momentum");
    assert!(result.wall_records.iter().all(|r| r.maximum_y_plus > 0.0));
    assert!(result.wall_records.windows(2).any(|r|
        (r[1].maximum_heat_transfer_w_m2_k - r[0].maximum_heat_transfer_w_m2_k).abs() > 1e-6));
    assert!(result.final_boundary_faces.iter().all(|face| face.surface_temperature_k == 310.0));
    assert!(result.march.energy_records.iter().all(|r| r.closure_j.abs() < 1e-7));
    assert!(result.march.temperature.iter().any(|t| *t > 300.0));
}

#[test]
fn internal_contact_and_external_cooling_follow_the_same_two_cell_storage_system() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::graded_from_fn([vec![0.006, 0.008], vec![0.02], vec![0.02]],
        |p| if p[0] < 0.006 { Voxel::Fluid } else { Voxel::Solid(0) }).unwrap();
    let mut slab = SolidMaterial::new("synthetic storage fixture", 16.0);
    slab.volumetric_heat_capacity_j_m3_k = Some(1000.0);
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![0.0, 0.004];
    setup.face_contacts.push(FaceContactResistance { cells: [0, 1], resistance_m2_k_w: 0.002 });
    let walls = WallHeatMarchConfig {
        internal_walls: vec![WallHeatFace { fluid_cell: 0, face: Face3::ALL[1] }],
        boundary_walls: vec![boundary(0, 0, ThermalFace::Temperature(300.0))], transport: controls() };
    let clock = time(0.02, 8);
    let result = march_conjugate_with_wall_heat(&domain, &fluid(), &[slab], &setup,
        &simple(), &clock, &walls, None, &[300.0, 300.0], |_| 1.0, |_| 1.0,
        &EnergyConfig::default(), None, &gate).unwrap();
    let area = domain.face_area(0, 0);
    let gf = area / (0.003 / fluid().conductivity_w_m_k);
    let gi = area / (0.003 / fluid().conductivity_w_m_k + 0.002 + 0.004 / 16.0);
    let cf = 1000.0 * domain.volume(0);
    let cs = 1000.0 * domain.volume(1);
    let mut temperature = [300.0, 300.0];
    for _ in 0..clock.steps {
        let a = cf / clock.time_step_s + gf + gi;
        let d = cs / clock.time_step_s + gi;
        let b = cf / clock.time_step_s * temperature[0] + gf * 300.0;
        let e = cs / clock.time_step_s * temperature[1] + 0.004;
        let determinant = a * d - gi * gi;
        temperature = [(b * d + gi * e) / determinant, (a * e + gi * b) / determinant];
    }
    close(result.march.temperature[0], temperature[0]);
    close(result.march.temperature[1], temperature[1]);
    close(cf * (temperature[0] - 300.0) + cs * (temperature[1] - 300.0)
        + result.boundary_outflow_j, 0.004 * clock.time_step_s * clock.steps as f64);
    assert_eq!(result.final_internal_faces[0].physical_contact_m2_k_w, 0.002);
    assert_eq!(result.final_internal_faces[0].correction_m2_k_w, 0.0);
}

#[test]
fn time_refinement_converges_to_independent_exponential_cooling() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let walls = WallHeatMarchConfig { internal_walls: vec![],
        boundary_walls: vec![boundary(0, 0, ThermalFace::Temperature(300.0))], transport: controls() };
    let exact = 300.0 + 10.0 * (-0.8f64).exp();
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let result = march_conjugate_with_wall_heat(&domain, &fluid(), &[], &setup,
            &simple(), &time(1.0 / steps as f64, steps), &walls, None, &[310.0],
            |_| 1.0, |_| 1.0, &EnergyConfig::default(), None, &gate).unwrap();
        errors.push((result.march.temperature[0] - exact).abs());
    }
    assert!(errors.windows(2).all(|e| e[0] / e[1] > 1.8 && e[0] / e[1] < 2.2), "{errors:?}");
}

#[test]
fn transport_conflicts_and_late_cancellation_publish_no_march_or_mutation() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let before = setup.clone();
    let walls = WallHeatMarchConfig { internal_walls: vec![],
        boundary_walls: vec![boundary(0, 0, ThermalFace::Temperature(300.0))], transport: controls() };
    let clock = time(0.02, 4);
    let mut laminar = simple(); laminar.turbulence = Turbulence::Laminar;
    assert!(march_conjugate_with_wall_heat(&domain, &fluid(), &[], &setup, &laminar,
        &clock, &walls, None, &[310.0], |_| 1.0, |_| 1.0,
        &EnergyConfig::default(), None, &gate).is_err());
    let mut conflict = setup.clone(); conflict.eddy_conductivity_w_m_k = vec![0.0];
    assert!(march_conjugate_with_wall_heat(&domain, &fluid(), &[], &conflict, &simple(),
        &clock, &walls, None, &[310.0], |_| 1.0, |_| 1.0,
        &EnergyConfig::default(), None, &gate).is_err());
    let result = march_conjugate_with_wall_heat(&domain, &fluid(), &[], &setup, &simple(),
        &clock, &walls, None, &[310.0], |t| { if t >= 0.04 { gate.request(); } 1.0 },
        |_| 1.0, &EnergyConfig::default(), None, &gate);
    assert!(matches!(result, Err(ChtError::Cancelled)));
    assert_eq!(setup, before);
}
