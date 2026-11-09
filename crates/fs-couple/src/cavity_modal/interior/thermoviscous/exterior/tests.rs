use super::*;
use fs_bem::reciprocal_radiation::Budget;
use fs_feec::acoustics::helmholtz::multiport::BoundaryMode;
use fs_feec::acoustics::thermoviscous::moving_wall::MovingWall;
use fs_feec::acoustics::thermoviscous::moving_wall::relaxation::MovingWallRelaxationLimits;
use fs_feec::acoustics::thermoviscous::{
    ThermoviscousBoundary, ThermoviscousFe, ThermoviscousLimits,
    ThermoviscousThermodynamics,
};
use fs_feec::fixtures::kuhn_cube;
use fs_material::gas::GasSpec;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_modal::time::{ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_rep_mesh::TetComplex;

const SIZE: [f64; 3] = [0.02, 0.006, 0.005];
const FREQUENCY: f64 = 620.0;
const BAND: [f64; 2] = [20.0, 2000.0];
const REFERENCES: [f64; 3] = [1.2e7, 2.1e7, 8.0e6];

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn exterior() -> Medium {
    // Independently declared exterior, deliberately unlike the source interior.
    Medium { density: 1.09, sound_speed: 352.0 }
}
fn budget() -> ThermoviscousRadiationBudget {
    ThermoviscousRadiationBudget {
        maximum_mixed_dofs: 1024, maximum_panels: 128,
        maximum_trace_entries: 1024, maximum_projection_terms: 32768,
    }
}
fn route() -> ApertureRoute {
    ApertureRoute::Galerkin {
        budget: Budget {
            maximum_unknowns: 256, maximum_cache_bytes: 32 * 1024 * 1024,
            maximum_work: 200_000_000, angular_refinement: 0,
        },
    }
}
fn on(points: &[[f64; 3]], face: [u32; 3], axis: usize, value: f64) -> bool {
    face.iter().all(|&v| (points[v as usize][axis] - value).abs() < 1e-14)
}
fn fixture() -> (TetComplex, Vec<[f64; 3]>, ThermoviscousBoundary, SpherePanels) {
    let (mesh, unit) = kuhn_cube(2);
    let points: Vec<_> = unit.iter().map(|p| std::array::from_fn(|i| p[i] * SIZE[i])).collect();
    let mut patches = [Vec::new(), Vec::new(), Vec::new()];
    let mut walls = Vec::new();
    let mut triangles = Vec::new();
    for &face in &mesh.faces {
        let boundary = (0..3).find_map(|axis| {
            if on(&points, face, axis, 0.0) { Some((axis, -1.0)) }
            else if on(&points, face, axis, SIZE[axis]) { Some((axis, 1.0)) }
            else { None }
        });
        let Some((axis, sign)) = boundary else { continue };
        if on(&points, face, 0, 0.0) { patches[0].push(face); }
        else if on(&points, face, 0, SIZE[0]) { patches[1].push(face); }
        else if on(&points, face, 2, SIZE[2]) { patches[2].push(face); }
        else { walls.push(face); }
        let mut triangle = face.map(|v| points[v as usize]);
        let a: [f64; 3] = std::array::from_fn(|i| triangle[1][i] - triangle[0][i]);
        let b: [f64; 3] = std::array::from_fn(|i| triangle[2][i] - triangle[0][i]);
        let cross = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
        if cross[axis] * sign < 0.0 { triangle.swap(1, 2); }
        triangle.rotate_left(1); // Orientation is physical, not FE vertex order.
        triangles.push(triangle);
    }
    let boundary = ThermoviscousBoundary {
        no_slip_faces: walls.clone(), isothermal_faces: walls,
        ports: patches.iter().map(|faces| BoundaryMode::uniform(faces)).collect(),
    };
    (mesh, points, boundary, SpherePanels::from_triangles(triangles).unwrap())
}
fn model<'c>(
    mesh: &'c TetComplex, points: &[[f64; 3]], boundary: &ThermoviscousBoundary,
    moving: bool,
) -> ThermoviscousInterior<'c> {
    let gas = gas();
    let medium = super::super::source_medium(&gas, 1.5e-5).unwrap();
    let limits = ThermoviscousLimits::default();
    let fe = if moving {
        let faces: Vec<_> = boundary.no_slip_faces.iter().copied()
            .filter(|&face| on(points, face, 1, 0.0)).collect();
        let mut shape = vec![[0.0; 3]; points.len() + mesh.edges.len()];
        // A conforming manufactured structural mode, zero where it meets any
        // other boundary. P2 edge coefficients are zero hierarchical bubbles.
        for (v, p) in points.iter().enumerate() {
            if p[1] == 0.0 && p[0] > 0.0 && p[0] < SIZE[0] && p[2] > 0.0 && p[2] < SIZE[2] {
                shape[v][1] = -20.0;
            }
        }
        let structure = ModalAcousticTimeModel::try_new(
            48_000,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: core::f64::consts::TAU * 900.0,
                damping_ratio: 0.02, pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        ).unwrap();
        ThermoviscousFe::new_with_moving_wall(
            mesh, points, medium, ThermoviscousThermodynamics::ideal_gas(medium).unwrap(),
            boundary, MovingWall { faces, velocity_shapes: vec![shape], structure },
            limits, &mut || false,
        ).unwrap().with_moving_wall_relaxation_channels(
            vec![vec![1000.0]], vec![0.002],
            MovingWallRelaxationLimits {
                maximum_channels: 4, maximum_workspace_bytes: 4096,
                maximum_projection_terms: 1024,
            }, &mut || false,
        ).unwrap()
    } else {
        ThermoviscousFe::new(mesh, points, medium, boundary, limits, &mut || false).unwrap()
    };
    // Manufactured original reservoir isolates complex caloric loss; this is
    // not a claim of measured air absorption or resolved-mesh convergence.
    let data = "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let relaxing = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    ThermoviscousInterior::from_relaxation(fe, law.evaluate(&relaxing).unwrap(), BAND).unwrap()
}

/// Independently author the constant exterior lifts from physical box planes.
/// This does not read the adapter's geometry map, apertures or native accessor.
fn native_bem(surface: &SpherePanels, selected: &[usize]) -> Vec<C64> {
    let planes = [(0, 0.0, -1.0), (0, SIZE[0], 1.0), (2, SIZE[2], 1.0)];
    let apertures: Vec<_> = selected.iter().map(|&mode| {
        let (axis, coordinate, sign) = planes[mode];
        let mut normal = [0.0; 3]; normal[axis] = sign;
        let mut reference_point = [0.0; 3]; reference_point[axis] = coordinate;
        let shape = surface.triangles().unwrap().iter().map(|triangle| {
            if triangle.iter().all(|p| (p[axis] - coordinate).abs() < 1e-14) { 1.0 } else { 0.0 }
        }).collect();
        Aperture { shape, reference_point, normal }
    }).collect();
    aperture_impedance(surface, &apertures, &[core::f64::consts::TAU * FREQUENCY],
        exterior(), &route(), &CancelGate::new()).unwrap().acoustic(0).unwrap()
        .into_iter().map(C64::conj).collect()
}

/// Independent native traction solve: Q=Y P+Q0, driven Q prescribed, and
/// P_r+Z Q_r=0. No scattering, numerical power references or adapter closure
/// is used. Every basis column and final field comes from the original LNS LU.
fn direct(
    model: &ThermoviscousInterior<'_>, selected: &[usize], q: &[C64], z: &[C64],
    wall_forces: &[C64], force: &[(&ThermoviscousBodyForceSource<'_, '_>, C64)],
    heat: &[(&ThermoviscousHeatSource<'_, '_>, C64)],
) -> ThermoviscousField {
    let native = model.prepare_response(FREQUENCY, &mut || false).unwrap();
    let n = model.inner().port_mode_count();
    let feeds: Vec<_> = (0..n).filter(|i| !selected.contains(i)).collect();
    let zeros = vec![C64::ZERO; n];
    let offset = native.response_with_volume_sources(&zeros, wall_forces, force, heat, &mut || false).unwrap();
    let mut y = vec![C64::ZERO; n*n];
    for column in 0..n {
        let mut drive = zeros.clone(); drive[column] = C64::ONE;
        let basis = native.response(&drive, &mut || false).unwrap();
        for row in 0..n { y[row*n+column] = basis.inward_flow_m3_s[row]; }
    }
    let mut matrix = vec![C64::ZERO; n*n];
    let mut rhs = zeros;
    for (&row, &flow) in feeds.iter().zip(q) {
        // Scale physical feed equations for the tiny dense reference solve.
        let scale = model.gas().density * model.gas().sound_speed / model.inner().port_areas_m2()[row];
        for column in 0..n { matrix[row*n+column] = y[row*n+column].scale(scale); }
        rhs[row] = (flow - offset.inward_flow_m3_s[row]).scale(scale);
    }
    for (r, &row) in selected.iter().enumerate() {
        matrix[row*n+row] = C64::ONE;
        for (s, &port) in selected.iter().enumerate() {
            for column in 0..n {
                matrix[row*n+column] = matrix[row*n+column] + z[r*selected.len()+s]*y[port*n+column];
            }
            rhs[row] = rhs[row] - z[r*selected.len()+s]*offset.inward_flow_m3_s[port];
        }
    }
    lu_complex_with_cancel(&matrix, n, || false).unwrap().solve(&mut rhs);
    native.response_with_volume_sources(&rhs, wall_forces, force, heat, &mut || false).unwrap()
}
fn near(a: &[C64], b: &[C64]) {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().chain(b).map(|z| z.abs()).fold(1e-30_f64, f64::max);
    for (&x, &y) in a.iter().zip(b) {
        assert!((x-y).abs() < 3e-8*scale, "{x:?} != {y:?}, scale={scale:e}");
    }
}
fn fields(actual: &ThermoviscousField, reference: &ThermoviscousField) {
    near(&actual.pressure_drive_pa, &reference.pressure_drive_pa);
    near(&actual.inward_flow_m3_s, &reference.inward_flow_m3_s);
    near(&actual.pressure_pa, &reference.pressure_pa);
    near(&actual.temperature_k, &reference.temperature_k);
    let a: Vec<_> = actual.velocity_m_s.iter().flatten().copied().collect();
    let b: Vec<_> = reference.velocity_m_s.iter().flatten().copied().collect();
    near(&a, &b);
    near(&actual.moving_wall.modal_velocity_m_sqrt_kg_per_s, &reference.moving_wall.modal_velocity_m_sqrt_kg_per_s);
    near(&actual.moving_wall.relaxation_coordinates_sqrt_j, &reference.moving_wall.relaxation_coordinates_sqrt_j);
}
fn balance(response: &ThermoviscousRadiationResponse) {
    let scale = response.total_source_power_w().abs() + response.exterior_power_w.abs()
        + response.interior.power.dissipated_w().abs();
    assert!(response.power_balance_residual_w.abs() < 5e-8*scale.max(1e-25));
    assert!(response.port_relative_residual <= DEFAULT_PORT_RELATIVE_TOLERANCE);
    assert!(response.exterior_power_w >= 0.0);
}

#[test]
fn full_molecular_lns_fields_match_native_traction_bem_closure_and_reference_changes() {
    let (mesh, points, boundary, surface) = fixture();
    let original = model(&mesh, &points, &boundary, false);
    let gate = CancelGate::new();
    let coupled = RadiatingThermoviscousInterior::new(&original, &[1, 2], &surface, exterior(), budget(), &gate).unwrap();
    assert!(core::ptr::eq(coupled.interior(), &original));
    assert_eq!(coupled.exterior_medium(), exterior());
    let prepared = coupled.prepare(FREQUENCY, &REFERENCES, &route(), &gate).unwrap();
    let z = native_bem(&surface, &[1, 2]);
    near(prepared.reaction_impedance(), &z);
    assert!(z[1].abs() > 0.0 && z[2].abs() > 0.0);
    let q = [C64::new(2e-8, -7e-9)];
    let actual = prepared.response(&q, &gate).unwrap();
    let reference = direct(&original, &[1, 2], &q, &z, &[], &[], &[]);
    fields(&actual.interior, &reference);
    balance(&actual);
    for value in [actual.interior.power.viscous_w, actual.interior.power.bulk_w,
        actual.interior.power.conductive_w, actual.interior.power.caloric_w] { assert!(value > 0.0); }
    let physical_work = 0.5*(reference.pressure_drive_pa[0].conj()*q[0]).re;
    assert!((actual.source_power_w-physical_work).abs() < 3e-8*physical_work.abs());
    assert!((actual.interior.power.input_w + actual.exterior_power_w - actual.total_source_power_w()).abs()
        < 3e-8*actual.total_source_power_w().abs());
    let alternate = coupled.prepare(FREQUENCY, &[3.7e6, 5.3e7, 1.9e7], &route(), &gate).unwrap();
    let changed = alternate.response(&q, &gate).unwrap();
    fields(&changed.interior, &reference);
    balance(&changed);
}

#[test]
fn original_wall_memory_body_force_and_heat_radiate_coherently_without_duplicate_work() {
    let (mesh, points, boundary, surface) = fixture();
    let original = model(&mesh, &points, &boundary, true);
    let volume = SIZE.iter().product::<f64>();
    let force = original.volume_force_source(4, &mut |_, x| {
        Ok([(1.0+x[0]/SIZE[0])/volume, -0.3/volume, 0.2/volume])
    }, &mut || false).unwrap();
    let heat = original.volume_heat_source(4, &mut |_, x| {
        Ok((0.5+x[2]/SIZE[2])/volume)
    }, &mut || false).unwrap();
    let gate = CancelGate::new();
    let coupled = RadiatingThermoviscousInterior::new(&original, &[1, 2], &surface, exterior(), budget(), &gate).unwrap();
    let prepared = coupled.prepare(FREQUENCY, &REFERENCES, &route(), &gate).unwrap();
    let walls = [C64::new(2e-5, 1e-5)];
    let forces = [(&force, C64::new(8e-7, -2e-7))];
    let heats = [(&heat, C64::new(3e-5, 8e-6))];
    let q = [C64::ZERO];
    let actual = prepared.response_with_volume_sources(&q, &walls, &forces, &heats, &gate).unwrap();
    let reference = direct(&original, &[1, 2], &q, &native_bem(&surface, &[1, 2]), &walls, &forces, &heats);
    fields(&actual.interior, &reference);
    balance(&actual);
    assert_eq!(actual.source_power_w, 0.0);
    assert!(actual.exterior_power_w > 0.0);
    assert!(actual.interior.power.body_force_work_w.abs() > 0.0);
    assert!(actual.interior.power.heat_source_work_w.abs() > 0.0);
    assert!(actual.interior.power.structural_damping_w > 0.0);
    assert!(actual.interior.power.structural_relaxation_w > 0.0);
    let quiet = prepared.response_with_volume_sources(
        &q, &[C64::ZERO], &[(&force, C64::ONE), (&force, -C64::ONE)],
        &[(&heat, C64::ONE), (&heat, -C64::ONE)], &gate,
    ).unwrap();
    assert!(quiet.interior.pressure_pa.iter().all(|&v| v == C64::ZERO));
    assert_eq!(quiet.total_source_power_w(), 0.0);
    assert_eq!(quiet.exterior_power_w, 0.0);
    let other = model(&mesh, &points, &boundary, false);
    let foreign = other.volume_heat_source(4, &mut |_, _| Ok(1.0/volume), &mut || false).unwrap();
    assert!(prepared.response_with_volume_sources(&q, &walls, &forces, &[(&foreign, C64::ONE)], &gate).is_err());
    let stopped = CancelGate::new(); stopped.request();
    assert!(matches!(prepared.response_with_volume_sources(&q, &walls, &forces, &heats, &stopped),
        Err(ThermoviscousRadiationError::Cancelled)));
    let retry = prepared.response_with_volume_sources(&q, &walls, &forces, &heats, &gate).unwrap();
    assert_eq!(retry.interior.pressure_pa, actual.interior.pressure_pa);
    assert_eq!(retry.interior.velocity_m_s, actual.interior.velocity_m_s);
    assert_eq!(retry.interior.temperature_k, actual.interior.temperature_k);
}

#[test]
fn original_geometry_duality_budgets_band_and_cancellation_are_admitted_before_response() {
    let (mesh, points, boundary, surface) = fixture();
    let original = model(&mesh, &points, &boundary, false);
    let gate = CancelGate::new();
    let triangles = original.inner().plane_wave_port_triangles(1).unwrap();
    let area: f64 = triangles.map(|t| t.area_m2).sum();
    assert!((area/(SIZE[1]*SIZE[2])-1.0).abs() < 1e-14);
    assert!(original.inner().plane_wave_port_triangles(3).is_none());
    assert!(RadiatingThermoviscousInterior::new(&original, &[1, 1], &surface, exterior(), budget(), &gate).is_err());
    let mut small = budget(); small.maximum_mixed_dofs = original.inner().ndof()-1;
    assert!(RadiatingThermoviscousInterior::new(&original, &[1], &surface, exterior(), small, &gate).is_err());
    let shifted = SpherePanels::from_triangles(surface.triangles().unwrap().iter()
        .map(|triangle| triangle.map(|p| [p[0]+0.1, p[1], p[2]])).collect()).unwrap();
    assert!(RadiatingThermoviscousInterior::new(&original, &[1], &shifted, exterior(), budget(), &gate).is_err());
    let mut signed = boundary.clone();
    for (face, values) in &mut signed.ports[1].traces {
        *values = face.map(|v| 2.0*points[v as usize][1]/SIZE[1]-1.0);
    }
    let signed_owner = model(&mesh, &points, &signed, false);
    assert!(signed_owner.inner().plane_wave_port_triangles(1).is_none());
    assert!(RadiatingThermoviscousInterior::new(&signed_owner, &[1], &surface, exterior(), budget(), &gate).is_err());
    let coupled = RadiatingThermoviscousInterior::new(&original, &[1, 2], &surface, exterior(), budget(), &gate).unwrap();
    assert!(coupled.prepare(BAND[1]+1.0, &REFERENCES, &route(), &gate).is_err());
    assert!(coupled.prepare(FREQUENCY, &[0.0, 1.0, 1.0], &route(), &gate).is_err());
    let stopped = CancelGate::new(); stopped.request();
    assert!(matches!(coupled.prepare(FREQUENCY, &REFERENCES, &route(), &stopped), Err(ThermoviscousRadiationError::Cancelled)));
}
