//! G1/G3: physical loads and material swaps reach the actual FEEC flow solver.

use fs_flux::{FluxSystem, MaterialFluxNumerics, MaterialFluxSystem, TriMesh};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_solid::Mesh2;

fn mesh() -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3))
}
fn no_slip(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn zero_force(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn transport(density: f64, viscosity: f64) -> FluidTransportState {
    // Numerical material fixtures, not claimed measurements of named fluids.
    FluidTransportState::try_new(density, 1000.0, 0.5, viscosity).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} versus {expected} (tol {tolerance})"
    );
}

#[test]
fn g1_hydrostatic_pressure_is_force_density_integral_in_pascals() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    for material in [transport(1.2, 0.03), transport(1000.0, 0.7)] {
        let flow =
            MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
        let solution = flow.solve_stokes(&|_| [3.0, -7.0], &no_slip).unwrap();
        let expected: Vec<_> = mesh
            .centroids
            .iter()
            .map(|point| 3.0 * point[0] - 7.0 * point[1])
            .collect();
        let area: f64 = mesh.areas.iter().sum();
        let mean: f64 = expected
            .iter()
            .zip(&mesh.areas)
            .map(|(p, a)| p * a)
            .sum::<f64>()
            / area;
        for (actual, p) in solution.pressure_pa().iter().zip(expected) {
            near(*actual, p - mean, 1e-10);
        }
        assert!(solution.velocity_dofs().iter().all(|u| u.abs() < 1e-10));
        assert!(std::ptr::eq(flow.material(), &material));
    }
}

#[test]
fn g3_swapping_dynamic_viscosity_changes_flow_without_copying_parameters() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let first = transport(2.0, 0.3);
    let second = transport(7.0, 0.6);
    let first_flow =
        MaterialFluxSystem::new(&system, &first, MaterialFluxNumerics::default()).unwrap();
    let second_flow =
        MaterialFluxSystem::new(&system, &second, MaterialFluxNumerics::default()).unwrap();
    let force = |point: [f64; 2]| [-point[1], point[0]];
    let a = first_flow.solve_stokes(&force, &no_slip).unwrap();
    let b = second_flow.solve_stokes(&force, &no_slip).unwrap();
    assert!(a.velocity_dofs().iter().any(|u| u.abs() > 1e-4));
    // Fixed volumetric load: Stokes velocity scales as 1/mu, NOT 1/nu.
    for (ua, ub) in a.velocity_dofs().iter().zip(b.velocity_dofs()) {
        near(*ub, 0.5 * ua, 2e-12);
    }
    for (pa, pb) in a.pressure_pa().iter().zip(b.pressure_pa()) {
        near(*pb, *pa, 2e-12);
    }
    // Energy must change with BOTH inertia rho and the changed velocity.
    near(
        b.kinetic_energy_j_per_m() / a.kinetic_energy_j_per_m(),
        (7.0 / 2.0) / 4.0,
        1e-10,
    );
}

#[test]
fn g1_shared_gas_pressure_changes_inertial_response_in_physical_seconds() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let spec = GasSpec::dry_air_ussa1976();
    let mut velocities = Vec::new();
    for pressure in [101_325.0, 202_650.0] {
        let gas = GasState::try_new(&spec, 293.15, pressure).unwrap();
        let flow = MaterialFluxSystem::new(&system, &gas, MaterialFluxNumerics::default()).unwrap();
        let dt = 0.01;
        let force = 0.25; // same N/m³ in both gases
        let expected_velocity = force * dt / flow.transport().density_kg_m3();
        let previous = vec![0.0; system.n_u];
        let solution = flow
            .bdf1_step(
                &|_| [force, 0.0],
                &|_| [expected_velocity, 0.0],
                &previous,
                dt,
            )
            .unwrap();
        let (error, divergence) =
            system.velocity_error(&solution.discrete().x, &|_| [expected_velocity, 0.0]);
        assert!(error < 1e-12, "velocity error {error}");
        assert!(divergence < 1e-12);
        assert!(solution.pressure_pa().iter().all(|p| p.abs() < 1e-10));
        near(
            solution.kinetic_energy_j_per_m(),
            0.5 * flow.transport().density_kg_m3() * expected_velocity.powi(2),
            1e-14,
        );
        velocities.push(system.velocity_at(&solution.discrete().x, [0.5, 0.5])[0]);
        assert_eq!(previous, vec![0.0; system.n_u]);
    }
    near(velocities[1], 0.5 * velocities[0], 1e-12);
}

#[test]
fn g3_picard_retains_the_same_material_and_creeping_flow_limit() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = transport(1000.0, 1.0);
    let flow =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let force = |point: [f64; 2]| [-1e-6 * point[1], 1e-6 * point[0]];
    let stokes = flow.solve_stokes(&force, &no_slip).unwrap();
    let (navier_stokes, iterations, update) =
        flow.solve_picard(&force, &no_slip, 12, 1e-9).unwrap();
    assert!(iterations > 0);
    assert!(update < 1e-9);
    for (s, n) in stokes
        .velocity_dofs()
        .iter()
        .zip(navier_stokes.velocity_dofs())
    {
        near(*n, *s, 1e-13);
    }
    assert_eq!(navier_stokes.transport(), &material);
}

#[test]
fn invalid_material_loads_and_steps_do_not_publish_a_solution() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = transport(1.0, 0.1);
    let flow =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    assert!(flow.solve_stokes(&|_| [f64::NAN, 0.0], &no_slip).is_err());
    assert!(
        flow.solve_stokes(&zero_force, &|_| [0.0, f64::INFINITY])
            .is_err()
    );
    assert!(flow.bdf1_step(&zero_force, &no_slip, &[], 0.1).is_err());
    let previous = vec![0.0; system.n_u];
    for dt in [0.0, -0.1, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(
            flow.bdf1_step(&zero_force, &no_slip, &previous, dt)
                .is_err()
        );
    }
    assert!(flow.solve_picard(&zero_force, &no_slip, 0, 1e-6).is_err());
    assert!(
        MaterialFluxSystem::new(
            &system,
            &material,
            MaterialFluxNumerics {
                penalty: f64::NAN,
                ..MaterialFluxNumerics::default()
            }
        )
        .is_err()
    );
    let mut gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    gas.density *= 2.0; // A stale public tuple must not evade the shared gas admission.
    assert!(MaterialFluxSystem::new(&system, &gas, MaterialFluxNumerics::default()).is_err());
    assert_eq!(previous, vec![0.0; system.n_u]);
}
