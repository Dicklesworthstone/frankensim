use super::*;
use crate::conjugate::{FlowFace, FlowField, SimpleReport, ThermalFace, Voxel};
use crate::conjugate::turbulence::law_of_the_wall_ratios;

fn config() -> WallHeatConfig {
    WallHeatConfig { turbulent_prandtl: 0.9, relative_tolerance: 1e-9, max_refinements: 14 }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-8 * b.abs().max(1.0), "{a:e} != {b:e}");
}

// Three graded cells: prescribed solid, free solid, prescribed fluid (or its
// reflection). Uniform tangential advection enters/leaves the fluid at 300 K.
// Only the accepted tangential velocity and the Spalding eddy trace are
// manufactured; the energy solve and solid-to-fluid observer are production.
fn fixture(face: Face3, speed: f64) -> (VoxelDomain, FluidProperties, FvFlow, ThermalSetup, WallHeatFace, usize) {
    let f = face as usize;
    let axis = f / 2;
    let positive = f % 2 == 1;
    let tangent = (axis + 1) % 3;
    let mut widths = [vec![0.02], vec![0.02], vec![0.02]];
    widths[axis] = if positive { vec![0.006, 0.008, 0.004] } else { vec![0.004, 0.008, 0.006] };
    let domain = VoxelDomain::graded_from_fn(widths, |point| {
        let fluid = if positive { point[axis] < 0.006 } else { point[axis] > 0.012 };
        if fluid { Voxel::Fluid } else { Voxel::Solid(0) }
    }).unwrap();
    let cell = (0..domain.cell_count()).find(|&c| domain.is_fluid(c)).unwrap();
    let solid = domain.neighbor(cell, f).unwrap();
    let hot = (0..domain.cell_count()).find(|&c| c != cell && c != solid).unwrap();
    let fluid = FluidProperties { density_kg_m3: 1.0, specific_heat_j_kg_k: 1000.0,
        conductivity_w_m_k: 0.04, kinematic_viscosity_m2_s: 1e-5 };
    let re = speed * (0.5 * domain.widths(cell)[axis]) / fluid.kinematic_viscosity_m2_s;
    let (tangent_eddy, secant) = law_of_the_wall_ratios(re);
    let mut faces = [FlowFace::Wall; 6];
    faces[2 * tangent] = FlowFace::Fixed;
    faces[2 * tangent + 1] = FlowFace::Free;
    let field = FlowField::from_face_velocity(&domain, faces, |a, _| {
        if a == tangent { speed } else { 0.0 }
    }).unwrap();
    let mut velocity = vec![[0.0; 3]; domain.cell_count()];
    velocity[cell][tangent] = speed;
    let mut eddy = vec![0.0; domain.cell_count()];
    eddy[cell] = tangent_eddy * fluid.kinematic_viscosity_m2_s;
    let mut old_heat = vec![0.0; domain.cell_count()];
    old_heat[cell] = (secant - 1.0) * fluid.kinematic_viscosity_m2_s;
    let flow = FvFlow { field, velocity_m_s: velocity, pressure_pa: vec![0.0; domain.cell_count()],
        eddy_viscosity_m2_s: eddy, heat_eddy_viscosity_m2_s: old_heat,
        report: SimpleReport { iterations: 0, mass_residual: 0.0, momentum_residual: 0.0,
            max_divergence_m3_s: 0.0, inflow_m3_s: 0.0, outflow_m3_s: 0.0,
            max_cell_reynolds: 0.0, momentum_krylov_iterations: 0, pressure_krylov_iterations: 0,
            fan: None, internal_fans: Vec::new() } };
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.faces[2 * tangent] = ThermalFace::Inflow { temperature: 300.0 };
    thermal.faces[2 * tangent + 1] = ThermalFace::Outflow { backflow_temperature: 300.0 };
    thermal.fixed_temperature = vec![(hot, 400.0), (cell, 300.0)];
    thermal.face_contacts.push(FaceContactResistance { cells: [solid, cell], resistance_m2_k_w: 0.002 });
    (domain, fluid, flow, thermal, WallHeatFace { fluid_cell: cell, face }, hot)
}

#[test]
fn graded_conjugate_wall_has_correct_joint_temperature_heat_and_contact_on_all_faces() {
    let gate = CancelGate::new_clock_free();
    let solids = [SolidMaterial::new("declared slab", 16.0)];
    for face in Face3::ALL {
        let (domain, fluid, flow, thermal, wall, hot) = fixture(face, 20.0);
        let result = solve_energy_with_wall_heat(&domain, &fluid, &solids, &flow, &thermal,
            &[wall], config(), &EnergyConfig::default(), &gate).unwrap();
        let report = result.faces[0];
        let axis = face as usize / 2;
        let solid = report.solid_cell;
        let area = domain.face_area(solid, axis);
        let solid_link = 0.5 * (domain.widths(hot)[axis] + domain.widths(solid)[axis]) / 16.0;
        let wall_link = 0.5 * domain.widths(solid)[axis] / 16.0
            + report.transfer.resistance_m2_k_w + 0.002;
        let heat = area * 100.0 / (solid_link + wall_link);
        close(result.energy.temperature[solid], 400.0 - heat * solid_link / area);
        close(result.energy.solid_to_fluid_heat_w(&domain), heat);
        close(result.energy.report.balance.fixed_cell_injection_w, 0.0);
        close(result.energy.report.balance.residual_w, 0.0);
        assert_eq!(report.physical_contact_m2_k_w, 0.002);
        assert!(report.correction_m2_k_w > 0.0);
        let key = (solid.min(wall.fluid_cell), solid.max(wall.fluid_cell));
        close(result.energy.face_contact_resistance_m2_k_w[&key], 0.002 + report.correction_m2_k_w);
        close(report.native_fluid_resistance_m2_k_w + report.correction_m2_k_w,
            report.transfer.resistance_m2_k_w);
        assert!(result.energy.temperature[solid] > 300.0 && result.energy.temperature[solid] < 400.0);
    }
}

#[test]
fn bulk_tangent_transport_and_wall_resistance_each_enter_once() {
    let gate = CancelGate::new_clock_free();
    let (domain, fluid, flow, thermal, wall, _) = fixture(Face3::ALL[0], 20.0);
    let original = thermal.clone();
    let prepared = prepare_wall_heat(&domain, &fluid, &flow, &thermal, &[wall], config(), &gate).unwrap();
    let cell = wall.fluid_cell;
    close(prepared.thermal.eddy_conductivity_w_m_k[cell],
        fluid.volumetric_heat_capacity() * flow.eddy_viscosity_m2_s[cell] / config().turbulent_prandtl);
    assert_ne!(prepared.thermal.eddy_conductivity_w_m_k[cell], flow.eddy_conductivity(&fluid)[cell]);
    assert_eq!(thermal, original);
    assert_eq!(prepared.thermal.face_contacts.len(), 1, "augment one physical contact, never duplicate it");
    assert!(prepare_wall_heat(&domain, &fluid, &flow, &prepared.thermal, &[wall], config(), &gate).is_err());
    let solids = [SolidMaterial::new("slab", 16.0)];
    let direct = solve_energy(&domain, &fluid, &solids, &flow.field, &prepared.thermal,
        &EnergyConfig::default(), &gate).unwrap();
    let wrapped = solve_energy_with_wall_heat(&domain, &fluid, &solids, &flow, &thermal,
        &[wall], config(), &EnergyConfig::default(), &gate).unwrap();
    assert_eq!(wrapped.energy, direct);
    let mut normal = flow.clone();
    normal.velocity_m_s[cell][0] = 500.0;
    assert_eq!(prepare_wall_heat(&domain, &fluid, &normal, &thermal, &[wall], config(), &gate).unwrap(), prepared,
        "normal speed must not manufacture wall friction");
}

#[test]
fn rest_wall_is_bit_identical_to_existing_conduction() {
    let gate = CancelGate::new_clock_free();
    let (domain, fluid, flow, thermal, wall, _) = fixture(Face3::ALL[2], 0.0);
    let solids = [SolidMaterial::new("slab", 16.0)];
    let ordinary = solve_energy(&domain, &fluid, &solids, &flow.field, &thermal,
        &EnergyConfig::default(), &gate).unwrap();
    let result = solve_energy_with_wall_heat(&domain, &fluid, &solids, &flow, &thermal,
        &[wall], config(), &EnergyConfig::default(), &gate).unwrap();
    assert_eq!(ordinary, result.energy);
    assert_eq!(result.faces[0].correction_m2_k_w, 0.0);
}

#[test]
fn changing_physical_contact_changes_heat_without_changing_wall_friction() {
    let gate = CancelGate::new_clock_free();
    let (domain, fluid, flow, mut thermal, wall, _) = fixture(Face3::ALL[1], 20.0);
    let solids = [SolidMaterial::new("slab", 16.0)];
    let before = solve_energy_with_wall_heat(&domain, &fluid, &solids, &flow, &thermal,
        &[wall], config(), &EnergyConfig::default(), &gate).unwrap();
    thermal.face_contacts[0].resistance_m2_k_w *= 10.0;
    let after = solve_energy_with_wall_heat(&domain, &fluid, &solids, &flow, &thermal,
        &[wall], config(), &EnergyConfig::default(), &gate).unwrap();
    assert_eq!(before.faces[0].transfer, after.faces[0].transfer);
    assert_eq!(before.faces[0].correction_m2_k_w, after.faces[0].correction_m2_k_w);
    assert!(after.energy.solid_to_fluid_heat_w(&domain) < before.energy.solid_to_fluid_heat_w(&domain));
    assert!(after.energy.temperature[after.faces[0].solid_cell] > before.energy.temperature[before.faces[0].solid_cell]);
}

#[test]
fn incompatible_or_duplicate_inputs_and_cancel_leave_originals_unchanged() {
    let gate = CancelGate::new_clock_free();
    let (domain, fluid, flow, thermal, wall, _) = fixture(Face3::ALL[0], 20.0);
    let original = thermal.clone();
    assert!(prepare_wall_heat(&domain, &fluid, &flow, &thermal, &[wall, wall], config(), &gate).is_err());
    let wrong = WallHeatFace { face: Face3::ALL[1], ..wall };
    assert!(prepare_wall_heat(&domain, &fluid, &flow, &thermal, &[wrong], config(), &gate).is_err());
    let mut incompatible = flow.clone();
    incompatible.eddy_viscosity_m2_s.fill(0.0);
    assert!(prepare_wall_heat(&domain, &fluid, &incompatible, &thermal, &[wall], config(), &gate).is_err(),
        "cannot compensate an incompatible bulk coefficient with negative thermal resistance");
    let mut duplicate = thermal.clone();
    let mut contact = duplicate.face_contacts[0];
    contact.cells.reverse();
    duplicate.face_contacts.push(contact);
    assert!(prepare_wall_heat(&domain, &fluid, &flow, &duplicate, &[wall], config(), &gate).is_err());
    let mut malformed = flow.clone();
    malformed.velocity_m_s.pop();
    assert!(prepare_wall_heat(&domain, &fluid, &malformed, &thermal, &[wall], config(), &gate).is_err());
    let tight = WallHeatConfig { max_refinements: 3, relative_tolerance: 1e-12, ..config() };
    assert!(prepare_wall_heat(&domain, &fluid, &flow, &thermal, &[wall], tight, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(prepare_wall_heat(&domain, &fluid, &flow, &thermal, &[wall], config(), &cancelled), Err(ChtError::Cancelled)));
    assert_eq!(thermal, original);
}
