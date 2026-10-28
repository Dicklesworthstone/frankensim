use super::*;
use crate::conjugate::{FlowFace, FlowField, SimpleReport, Voxel};
use crate::conjugate::turbulence::law_of_the_wall_ratios;
use crate::d3q19::Face3;

fn controls() -> WallHeatConfig {
    WallHeatConfig { turbulent_prandtl: 0.9, relative_tolerance: 1e-9, max_refinements: 14 }
}
fn fluid() -> FluidProperties {
    FluidProperties { density_kg_m3: 1.0, specific_heat_j_kg_k: 1000.0,
        conductivity_w_m_k: 0.04, kinematic_viscosity_m2_s: 1e-5 }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-8 * b.abs().max(1.0), "{a:e} != {b:e}");
}
// Manufactured accepted transport image; all heat assembly/solves are native.
fn flow(domain: &VoxelDomain, face: Face3, speed: f64) -> FvFlow {
    let axis = face as usize / 2;
    let mut velocity = vec![[0.0; 3]; domain.cell_count()];
    let mut eddy = vec![0.0; domain.cell_count()];
    for cell in 0..domain.cell_count() {
        if domain.is_fluid(cell) {
            velocity[cell][(axis + 1) % 3] = speed;
            let re = speed.abs() * 0.5 * domain.widths(cell)[axis] / fluid().kinematic_viscosity_m2_s;
            eddy[cell] = law_of_the_wall_ratios(re).0 * fluid().kinematic_viscosity_m2_s;
        }
    }
    FvFlow { field: FlowField::quiescent(domain), velocity_m_s: velocity,
        pressure_pa: vec![0.0; domain.cell_count()], heat_eddy_viscosity_m2_s: eddy.clone(),
        eddy_viscosity_m2_s: eddy,
        report: SimpleReport { iterations: 0, mass_residual: 0.0, momentum_residual: 0.0,
            max_divergence_m3_s: 0.0, inflow_m3_s: 0.0, outflow_m3_s: 0.0,
            max_cell_reynolds: 0.0, momentum_krylov_iterations: 0, pressure_krylov_iterations: 0,
            fan: None, internal_fans: Vec::new() } }
}
fn wall(face: Face3, condition: ThermalFace) -> BoundaryWallHeat {
    BoundaryWallHeat { wall: WallHeatFace { fluid_cell: 0, face }, condition, contact_m2_k_w: 0.0 }
}

#[test]
fn exterior_temperature_and_film_preserve_series_heat_on_six_graded_faces() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::graded([vec![0.006], vec![0.008], vec![0.012]]).unwrap();
    for face in Face3::ALL {
        let flow = flow(&domain, face, 20.0);
        for condition in [ThermalFace::Temperature(300.0),
            ThermalFace::Convective { h: 25.0, ambient: 300.0 }]
        {
            let boundary = BoundaryWallHeat { contact_m2_k_w: 0.002, ..wall(face, condition) };
            let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
            setup.power_w = vec![0.002];
            let result = solve_energy_with_wall_heat_boundaries(&domain, &fluid(), &[], &flow,
                &setup, &[FvBoundary::wall(); 6], &[], &[boundary], controls(),
                &EnergyConfig::default(), &gate).unwrap();
            let exchange = result.boundary_faces[0];
            let r = exchange.coefficients;
            let total = r.transfer.resistance_m2_k_w + r.exterior_resistance_m2_k_w;
            close(result.energy.temperature[0], 300.0 + 0.002 * total / r.area_m2);
            close(exchange.surface_temperature_k, 300.0 + 0.002 * r.exterior_resistance_m2_k_w / r.area_m2);
            close(exchange.outward_w, 0.002);
            close(result.energy.report.balance.sink_outflow_w, 0.002);
            close(result.energy.report.balance.residual_w, 0.0);
            assert_eq!(result.energy.report.balance.boundary_outflow_w, 0.0,
                "delegated walls are counted once by native linear exchanges");
        }
    }
}

#[test]
fn imposed_flux_sets_watts_and_wall_law_sets_surface_temperature() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let hot = wall(Face3::ALL[0], ThermalFace::HeatFlux(12.0));
    let cold = wall(Face3::ALL[1], ThermalFace::Temperature(300.0));
    for speed in [0.0, 20.0] {
        let flow = flow(&domain, Face3::ALL[0], speed);
        let result = solve_energy_with_wall_heat_boundaries(&domain, &fluid(), &[], &flow,
            &setup, &[FvBoundary::wall(); 6], &[], &[hot, cold], controls(),
            &EnergyConfig::default(), &gate).unwrap();
        let input = result.boundary_faces[0];
        let output = result.boundary_faces[1];
        close(input.outward_w, -0.0012);
        close(output.outward_w, 0.0012);
        close(input.surface_temperature_k,
            result.energy.temperature[0] + 12.0 * input.coefficients.transfer.resistance_m2_k_w);
        assert_eq!(output.surface_temperature_k, 300.0);
        close(result.energy.report.balance.source_w, 0.0012);
        close(result.energy.report.balance.residual_w, 0.0);
    }
}

#[test]
fn wall_relative_speed_is_frame_invariant_and_repeated_queries_are_cached() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let walls = [wall(Face3::ALL[0], ThermalFace::Temperature(300.0)),
        wall(Face3::ALL[1], ThermalFace::Convective { h: 30.0, ambient: 310.0 })];
    let original = flow(&domain, Face3::ALL[0], 20.0);
    let a = prepare_wall_heat_boundaries(&domain, &fluid(), &original, &setup,
        &[FvBoundary::wall(); 6], &[], &walls, controls(), &gate).unwrap();
    assert!(a.quadrature_evaluations > 0);
    assert_eq!(a.quadrature_evaluations, a.boundary_faces[0].transfer.quadrature_evaluations);
    let mut shifted = original.clone();
    shifted.velocity_m_s[0][1] += 7.0;
    shifted.velocity_m_s[0][0] = 500.0; // not a tangential wall-friction input
    let mut mechanical = [FvBoundary::wall(); 6];
    mechanical[0] = FvBoundary::Wall { velocity: [0.0, 7.0, 0.0] };
    mechanical[1] = mechanical[0];
    let b = prepare_wall_heat_boundaries(&domain, &fluid(), &shifted, &setup,
        &mechanical, &[], &walls, controls(), &gate).unwrap();
    assert_eq!(a, b);
    shifted.velocity_m_s[0][1] = 7.0;
    let rest = prepare_wall_heat_boundaries(&domain, &fluid(), &shifted, &setup,
        &mechanical, &[], &walls, controls(), &gate).unwrap();
    assert_eq!(rest.boundary_faces[0].transfer.resistance_m2_k_w, 0.005 / fluid().conductivity_w_m_k);
    assert_eq!(rest.quadrature_evaluations, 0);
}

#[test]
fn internal_solid_and_exterior_wall_share_one_conservative_solve() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::graded_from_fn([vec![0.006, 0.008], vec![0.02], vec![0.02]],
        |p| if p[0] < 0.006 { Voxel::Fluid } else { Voxel::Solid(0) }).unwrap();
    let flow = flow(&domain, Face3::ALL[0], 20.0);
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.power_w = vec![0.0, 0.004];
    let internal = WallHeatFace { fluid_cell: 0, face: Face3::ALL[1] };
    let external = wall(Face3::ALL[0], ThermalFace::Temperature(300.0));
    let result = solve_energy_with_wall_heat_boundaries(&domain, &fluid(),
        &[SolidMaterial::new("slab", 16.0)], &flow, &setup,
        &[FvBoundary::wall(); 6], &[internal], &[external], controls(),
        &EnergyConfig::default(), &gate).unwrap();
    let area = domain.face_area(0, 0);
    let r_out = result.boundary_faces[0].coefficients.transfer.resistance_m2_k_w;
    let r_in = result.internal_faces[0].transfer.resistance_m2_k_w + 0.004 / 16.0;
    close(result.energy.temperature[0], 300.0 + 0.004 * r_out / area);
    close(result.energy.temperature[1], 300.0 + 0.004 * (r_out + r_in) / area);
    close(result.energy.solid_to_fluid_heat_w(&domain), 0.004);
    close(result.boundary_faces[0].outward_w, 0.004);
    close(result.energy.report.balance.residual_w, 0.0);
}

#[test]
fn invalid_boundary_ownership_motion_and_cancellation_do_not_modify_inputs() {
    let gate = CancelGate::new_clock_free();
    let domain = VoxelDomain::new(1, 1, 1, 0.01).unwrap();
    let flow = flow(&domain, Face3::ALL[0], 20.0);
    let setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let before = setup.clone();
    let boundary = wall(Face3::ALL[0], ThermalFace::Temperature(300.0));
    let mechanical = [FvBoundary::wall(); 6];
    let prepare = |f: &FvFlow, s: &ThermalSetup, m: &[FvBoundary; 6], w: &[BoundaryWallHeat], g: &CancelGate|
        prepare_wall_heat_boundaries(&domain, &fluid(), f, s, m, &[], w, controls(), g);
    assert!(prepare(&flow, &setup, &mechanical, &[boundary, boundary], &gate).is_err());
    let mut owned = setup.clone(); owned.faces[0] = ThermalFace::Temperature(300.0);
    assert!(prepare(&flow, &owned, &mechanical, &[boundary], &gate).is_err());
    for condition in [FvBoundary::Symmetry, FvBoundary::Outlet,
        FvBoundary::Wall { velocity: [1.0, 0.0, 0.0] }]
    {
        let mut bad = mechanical; bad[0] = condition;
        assert!(prepare(&flow, &setup, &bad, &[boundary], &gate).is_err());
    }
    let mut crossing = flow.clone();
    crossing.field = FlowField::from_face_velocity(&domain, [FlowFace::Fixed; 6],
        |axis, _| if axis == 0 { 1.0 } else { 0.0 }).unwrap();
    assert!(prepare(&crossing, &setup, &mechanical, &[boundary], &gate).is_err());
    let prepared = prepare(&flow, &setup, &mechanical, &[boundary], &gate).unwrap();
    assert!(prepare(&flow, &prepared.thermal, &mechanical, &[boundary], &gate).is_err());
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(prepare(&flow, &setup, &mechanical, &[boundary], &cancelled), Err(ChtError::Cancelled)));
    assert_eq!(setup, before);
}
