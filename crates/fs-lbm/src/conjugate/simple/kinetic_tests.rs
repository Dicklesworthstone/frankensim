//! Analytic native face-inertia norms, including physical boundary ownership.

use super::*;
use crate::conjugate::{UnsteadyConfig, Voxel, simple_unsteady};

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 2.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

#[test]
fn opposed_staggered_faces_retain_vortex_kinetic_energy() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = SimpleConfig::new([FvBoundary::wall(); 6]);
    let mut solver = Solver::new(&domain, &fluid(), &config);
    // One divergence-free circulation: four native face masses of 2 kg,
    // moving at 2 m/s. Averaging first loses half of its 16 J.
    solver.vel[0][1] = 2.0;
    solver.vel[0][4] = -2.0;
    solver.vel[1][2] = -2.0;
    solver.vel[1][3] = 2.0;
    assert_eq!(solver.field().max_divergence(&domain), 0.0);
    let centred: f64 = (0..domain.cell_count())
        .map(|c| solver.cell_velocity(c).iter().map(|v| v * v).sum::<f64>())
        .sum();
    assert_eq!(centred, 8.0);
    assert_eq!(solver.kinetic_energy_j(&CancelGate::new()).unwrap(), 16.0);
}

#[test]
fn graded_open_faces_own_only_their_physical_half_cells() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![0.5]]).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet {
        velocity: [2.0, 0.0, 0.0],
    };
    faces[1] = FvBoundary::Outlet;
    let config = SimpleConfig::new(faces);
    let mut solver = Solver::new(&domain, &fluid(), &config);
    solver.vel[0].copy_from_slice(&[2.0, 3.0, 5.0]);
    // Face area = 1 m^2. The domain owns masses [1, 4, 3] kg,
    // including half of each exposed cell, never the outlet's ghost volume.
    assert_eq!(solver.kinetic_energy_j(&CancelGate::new()).unwrap(), 57.5);
}

#[test]
fn uniform_motion_on_graded_faces_recovers_total_fluid_mass() {
    let domain =
        VoxelDomain::graded([vec![0.25, 1.0, 2.0], vec![1.0, 3.0], vec![0.5, 2.0]]).unwrap();
    let config = SimpleConfig::new([FvBoundary::Outlet; 6]);
    let mut solver = Solver::new(&domain, &fluid(), &config);
    for (component, speed) in solver.vel.iter_mut().zip([2.0, -3.0, 4.0]) {
        component.fill(speed);
    }
    // rho*V/2 = 32.5 kg, |u|^2 = 29 m^2/s^2.
    assert_eq!(solver.kinetic_energy_j(&CancelGate::new()).unwrap(), 942.5);
}

#[test]
fn stationary_solid_voxels_and_blocked_faces_supply_no_fluid_inertia() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet {
        velocity: [2.0, 0.0, 0.0],
    };
    faces[1] = FvBoundary::Outlet;
    let config = SimpleConfig::new(faces);
    let mut solver = Solver::new(&domain, &fluid(), &config);
    solver.vel[0][3] = 4.0;
    // Two exposed half-cell masses of 1 kg; the solid faces stay at rest.
    assert_eq!(solver.kinetic_energy_j(&CancelGate::new()).unwrap(), 10.0);
}

#[test]
fn ordinary_unsteady_march_records_native_staggered_storage() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [0.01, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let config = SimpleConfig::new(faces);
    let mut time = UnsteadyConfig::new(0.1, 1);
    time.scheme = TimeScheme::BackwardEuler;
    time.inner_iterations = 1000;
    time.inner_tolerance = 1e-11;
    let result = simple_unsteady(
        &domain,
        &fluid(),
        &config,
        &time,
        |_| 1.0,
        &CancelGate::new(),
    )
    .unwrap();
    let centred: f64 = result
        .flow
        .velocity_m_s
        .iter()
        .map(|v| v.iter().map(|x| x * x).sum::<f64>())
        .sum();
    assert!(centred > 1e-12);
    // Each of the four moving interior faces enters two cell averages at
    // half amplitude. Their native kinetic norm is therefore twice this one.
    let native = result.records[0].kinetic_energy_j;
    assert!((native - 2.0 * centred).abs() < 1e-13 * native);
}
