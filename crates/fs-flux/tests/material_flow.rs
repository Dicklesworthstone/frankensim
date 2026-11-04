//! G1/G3: physical loads and material swaps reach the actual FEEC flow solver.

use fs_flux::{
    FluxRheologyError, FluxSystem, MaterialFluxError, MaterialFluxNumerics,
    MaterialFluxRheologyPoint, MaterialFluxSystem, MaterialStokesBoundary, StokesBoundaryError,
    StokesPressureReference, TriMesh,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::gas::{GasSpec, GasState};
use fs_solid::Mesh2;
use std::cell::Cell;
use std::rc::Rc;

#[path = "material_flow/navier_stokes.rs"]
mod navier_stokes;

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
fn viscous_heat_is_objective_and_matches_couette_and_extensional_work() {
    // A = [[a, s-w], [w, -a]] has gamma²=4*a²+s². Translation and
    // rigid spin w do not contribute, even when their kinetic energy is large.
    // For Couette shear the independent wall-power check is tau*U*L=mu*s²*area.
    for n in [2, 4] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(2.0, 1.0, 2 * n, n));
        let system = FluxSystem::new(&mesh);
        for (a, shear, spin, boost) in [
            (0.0, 0.0, 1.7, [2.0, -3.0]),
            (0.0, 1.2, 0.0, [0.0, 0.0]),
            (0.35, -0.8, 2.1, [-3.0, 4.0]),
        ] {
            for density in [2.0, 7.0] {
                let mu = 0.3;
                let material = transport(density, mu);
                let flow =
                    MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default())
                        .unwrap();
                let velocity = |p: [f64; 2]| {
                    [
                        boost[0] + a * p[0] + (shear - spin) * p[1],
                        boost[1] + spin * p[0] - a * p[1],
                    ]
                };
                let solution = flow.solve_stokes(&zero_force, &velocity).unwrap();
                let heating = mu * (4.0 * a * a + shear * shear);
                let dissipation = solution.viscous_dissipation();
                assert_eq!(dissipation.cell_heating_w_m3().len(), mesh.tris.len());
                for &cell_heat in dissipation.cell_heating_w_m3() {
                    assert!(cell_heat >= 0.0);
                    near(cell_heat, heating, 2e-9);
                }
                near(dissipation.total_power_w_per_m(), 2.0 * heating, 2e-9);
                // The same result comes from the resolved volume, with no
                // dependence on the number or areas of the triangles.
                let integral: f64 = dissipation
                    .cell_heating_w_m3()
                    .iter()
                    .zip(&mesh.areas)
                    .map(|(q, area)| q * area)
                    .sum();
                near(dissipation.total_power_w_per_m(), integral, 2e-14);
            }
        }
    }
}

#[test]
fn transient_pressure_impulse_uses_original_gas_density_and_reports_pascals() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let spec = GasSpec::dry_air_ussa1976();
    let gradient = 2.5; // physical pressure drop per metre [Pa/m]
    let dt = 0.04;
    let mut speeds = Vec::new();
    for thermodynamic_pressure in [101_325.0, 202_650.0] {
        let gas = GasState::try_new(&spec, 293.15, thermodynamic_pressure).unwrap();
        let flow = MaterialFluxSystem::new(&system, &gas, MaterialFluxNumerics::default()).unwrap();
        let density = flow.transport().density_kg_m3();
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let pressure = 10.0 - gradient * point[0];
            [-pressure * normal[0], -pressure * normal[1]]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, _)| MaterialStokesBoundary::TotalTraction {
                edge,
                value: &traction,
            })
            .collect();
        let previous = vec![0.0; system.n_u];
        let solved = flow
            .bdf1_stokes_mixed_step(&zero_force, &boundaries, &previous, dt)
            .unwrap();
        let expected_speed = gradient * dt / density;
        let (error, divergence) =
            system.velocity_error(&solved.discrete().x, &|_| [expected_speed, 0.0]);
        assert!(error < 2e-11 && divergence < 2e-11);
        assert_eq!(
            solved.pressure_reference(),
            StokesPressureReference::TractionDatum
        );
        for (cell, point) in mesh.centroids.iter().enumerate() {
            near(solved.pressure_pa()[cell], 10.0 - gradient * point[0], 2e-9);
            assert_eq!(
                solved.pressure_reference_at(cell),
                Some(StokesPressureReference::TractionDatum)
            );
        }
        near(
            solved.kinetic_energy_j_per_m(),
            0.5 * density * expected_speed.powi(2),
            2e-12,
        );
        assert_eq!(solved.transport(), flow.transport());
        assert_eq!(previous, vec![0.0; system.n_u]);
        speeds.push(system.velocity_at(&solved.discrete().x, [0.4, 0.4])[0]);
    }
    near(speeds[1] / speeds[0], 0.5, 2e-9);
}

#[test]
fn transient_mixed_material_shear_uses_dynamic_viscosity_and_inertia_together() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let dt = 0.2;
    // u(t,y) = t*(1.1 + 0.6*y), p = 8 - 2*x. This affine field is in BDM1.
    let velocity = |point: [f64; 2]| [dt * (1.1 + 0.6 * point[1]), 0.0];
    for material in [transport(2.0, 0.3), transport(7.0, 1.7)] {
        let flow =
            MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
        let density = flow.transport().density_kg_m3();
        let viscosity = flow.transport().dynamic_viscosity_pa_s();
        let force = |point: [f64; 2]| [density * (1.1 + 0.6 * point[1]) - 2.0, 0.0];
        let traction = |point: [f64; 2], normal: [f64; 2]| {
            let pressure = 8.0 - 2.0 * point[0];
            let shear = viscosity * dt * 0.6;
            [
                shear * normal[1] - pressure * normal[0],
                shear * normal[0] - pressure * normal[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, geometry)| {
                if geometry.normal[1].abs() > 0.5 {
                    MaterialStokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                } else {
                    MaterialStokesBoundary::TotalTraction {
                        edge,
                        value: &traction,
                    }
                }
            })
            .collect();
        let previous = vec![0.0; system.n_u];
        let solved = flow
            .bdf1_stokes_mixed_step(&force, &boundaries, &previous, dt)
            .unwrap();
        let (error, divergence) = system.velocity_error(&solved.discrete().x, &velocity);
        assert!(error < 3e-11 && divergence < 3e-11);
        for (pressure, point) in solved.pressure_pa().iter().zip(&mesh.centroids) {
            near(*pressure, 8.0 - 2.0 * point[0], 3e-9);
        }
        let expected_energy =
            0.5 * density * dt * dt * (1.1_f64.powi(2) + 1.1 * 0.6 + 0.6_f64.powi(2) / 3.0);
        near(solved.kinetic_energy_j_per_m(), expected_energy, 2e-11);
        let heating = viscosity * (0.6 * dt).powi(2);
        for &cell_heat in solved.viscous_dissipation().cell_heating_w_m3() {
            near(cell_heat, heating, 3e-11);
        }
        near(
            solved.viscous_dissipation().total_power_w_per_m(),
            heating,
            3e-11,
        );
        let queries = Cell::new(0);
        let refused_force = |_: [f64; 2]| {
            queries.set(queries.get() + 1);
            [0.0; 2]
        };
        assert!(matches!(
            flow.bdf1_stokes_mixed_step(&refused_force, &boundaries, &previous, 0.0),
            Err(MaterialFluxError::Boundary(
                StokesBoundaryError::InvalidInput(_)
            ))
        ));
        assert_eq!(queries.get(), 0);
        assert_eq!(previous, vec![0.0; system.n_u]);
    }
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

// G0: a pressure pin cannot replace the closed-domain mass constraint.
#[test]
fn material_flow_refuses_opposite_leaks_in_disconnected_chambers() {
    let mesh = TriMesh::from_mesh2(&Mesh2 {
        nodes: vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [2.0, 1.0],
        ],
        elems: vec![vec![0, 1, 2], vec![3, 4, 5]],
        patches: Vec::new(),
    });
    let system = FluxSystem::new(&mesh);
    let material = transport(1.2, 0.03);
    let flow =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let force = |_| panic!("force must not be evaluated before local flux admission");
    let boundary = |p: [f64; 2]| {
        if p[0] < 1.5 {
            [p[0], 0.0]
        } else {
            [2.0 - p[0], 0.0]
        }
    };
    assert!(matches!(
        flow.assemble_stokes(&force, &boundary),
        Err(MaterialFluxError::Boundary(StokesBoundaryError::IncompatibleFlux { component: 0, net_flux, .. }))
            if (net_flux - 0.5).abs() < 1e-13
    ));
}

// G0: a pressure pin cannot replace the closed-domain mass constraint.
#[test]
fn incompatible_boundary_flux_is_refused_before_force_evaluation() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = transport(1.2, 0.03);
    let flow =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.125; system.n_u];
    let retained = previous.clone();
    let force_calls = Cell::new(0usize);
    let force = |_| {
        force_calls.set(force_calls.get() + 1);
        [0.0; 2]
    };
    for sign in [1.0, -1.0] {
        let boundary = |point: [f64; 2]| [sign * point[0], 0.0];
        for entry in 0..4 {
            let result = match entry {
                0 => flow.assemble_stokes(&force, &boundary).map(|_| ()),
                1 => flow.solve_stokes(&force, &boundary).map(|_| ()),
                2 => flow.solve_picard(&force, &boundary, 12, 1e-9).map(|_| ()),
                _ => flow
                    .bdf1_step(&force, &boundary, &previous, 0.1)
                    .map(|_| ()),
            };
            assert_eq!(
                force_calls.get(),
                0,
                "entry {entry} queried a force before refusal"
            );
            let Err(MaterialFluxError::Boundary(StokesBoundaryError::IncompatibleFlux {
                net_flux,
                tolerance,
                ..
            })) = result
            else {
                panic!("entry {entry} did not refuse incompatible flux: {result:?}");
            };
            near(net_flux, sign, 1e-14);
            assert!(net_flux.abs() > tolerance && tolerance > 0.0);
            assert_eq!(previous, retained);
        }
    }
}

// G1: balance is integral, not a blanket refusal of nonzero boundary flow.
#[test]
fn balanced_translation_and_affine_strain_preserve_incompressible_pressure() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = transport(1.2, 0.03);
    let flow =
        MaterialFluxSystem::new(&system, &material, MaterialFluxNumerics::default()).unwrap();
    for strain in [0.0, 0.1] {
        let boundary = |point: [f64; 2]| [0.3 + strain * point[0], -0.2 - strain * point[1]];
        let solution = flow.solve_stokes(&zero_force, &boundary).unwrap();
        let (error, divergence) = system.velocity_error(&solution.discrete().x, &boundary);
        assert!(error < 1e-11, "affine velocity error {error}");
        assert!(divergence < 1e-11, "affine divergence {divergence}");
        assert!(solution.pressure_pa().iter().all(|p| p.abs() < 1e-10));
    }
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

// Synthetic constitutive owner with explicit query conditions, not empirical
// data for a named fluid. Its additional fields must survive the flow solve.
struct RheologyResponse {
    transport: FluidTransportState,
    point: MaterialFluxRheologyPoint,
    temperature_k: f64,
    reads: Rc<Cell<usize>>,
}
impl FluidTransport for RheologyResponse {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        Ok(self.transport)
    }
}

#[test]
fn g1_material_rheology_closes_temperature_layers_and_physical_shear_traction() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 4, 4));
    let system = FluxSystem::new(&mesh);
    let lower_rate = 0.5 * (5.0_f64.sqrt() - 1.0);
    let upper_rate = 0.5 * (2.0_f64.sqrt() - 1.0);
    let exact = |point: [f64; 2]| {
        [
            lower_rate * point[1].min(0.5) + upper_rate * (point[1] - 0.5).max(0.0),
            0.0,
        ]
    };
    let dt = 0.1;
    for density in [2.0, 7.0] {
        let reference = transport(density, 1.0);
        let flow = MaterialFluxSystem::new(
            &system,
            &reference,
            MaterialFluxNumerics {
                linear_tolerance: 1e-11,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        let previous = vec![0.0; system.n_u];
        let force = |point| [density * exact(point)[0] / dt, 0.0];
        let reads = Rc::new(Cell::new(0));
        let mut calls = 0;
        let mut last_query = vec![None; mesh.tris.len()];
        let step = flow
            .bdf1_rheology_step(
                &force,
                &exact,
                &previous,
                dt,
                60,
                1e-10,
                |point, retained| {
                    assert_eq!(retained, &reference);
                    calls += 1;
                    last_query[point.cell] = Some(point);
                    let (temperature_k, scale) = if point.position_m[1] < 0.5 {
                        (330.0, 1.0)
                    } else {
                        (300.0, 4.0)
                    };
                    // mu(T,gamma)=mu_T*(1 + gamma * 1 s). The independent rates
                    // above solve mu*gamma=1 Pa on BOTH sides of the interface.
                    Ok(RheologyResponse {
                        transport: transport(density, scale * (1.0 + point.shear_rate_per_s)),
                        point,
                        temperature_k,
                        reads: Rc::clone(&reads),
                    })
                },
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&step.discrete().x, &exact);
        assert!(error < 2e-10, "layered nonlinear velocity error {error}");
        assert!(
            divergence < 2e-11,
            "layered nonlinear divergence {divergence}"
        );
        assert!(step.relative_update <= 1e-10 && step.nonlinear_residual <= 1e-10);
        assert!(
            step.iterations > 1,
            "the source law must be solved at the new velocity"
        );
        assert_eq!(step.samples().len(), mesh.tris.len());
        assert_eq!(
            reads.get(),
            calls,
            "each original source response is read once"
        );
        assert!(calls <= 2 * mesh.tris.len() * step.iterations);
        let final_rates = system.shear_rates(&step.discrete().x);
        for (cell, sample) in step.samples().iter().enumerate() {
            let point = sample.point();
            assert_eq!(point.cell, cell);
            assert_eq!(point.position_m, mesh.centroids[cell]);
            assert_eq!(Some(point), last_query[cell]);
            assert_eq!(sample.material().point, point);
            assert_eq!(point.shear_rate_per_s, final_rates[cell]);
            assert_eq!(sample.transport(), &sample.material().transport);
            near(
                sample.transport().dynamic_viscosity_pa_s() * point.shear_rate_per_s,
                1.0,
                2e-9,
            );
            // The independent uniform shear traction is 1 Pa, so local heat
            // is tau*gamma, despite different temperature/rheology on each side.
            near(
                step.viscous_dissipation().cell_heating_w_m3()[cell],
                if point.position_m[1] < 0.5 {
                    lower_rate
                } else {
                    upper_rate
                },
                3e-9,
            );
            near(
                sample.material().temperature_k,
                if point.position_m[1] < 0.5 {
                    330.0
                } else {
                    300.0
                },
                0.0,
            );
        }
        let expected_energy = 0.5
            * density
            * (lower_rate.powi(2) / 6.0
                + lower_rate * upper_rate / 8.0
                + upper_rate.powi(2) / 24.0);
        near(step.kinetic_energy_j_per_m(), expected_energy, 2e-10);
        near(
            step.viscous_dissipation().total_power_w_per_m(),
            0.5 * (lower_rate + upper_rate),
            3e-9,
        );
        assert_eq!(
            reads.get(),
            calls,
            "heat reporting does not requery sources"
        );
        assert!(step.pressure_pa().iter().all(|p| p.abs() < 2e-9));
        // A single reference viscosity cannot reproduce the sourced shear
        // traction balance, even though loads and boundary data are identical.
        let frozen = flow.bdf1_step(&force, &exact, &previous, dt).unwrap();
        let frozen_error = system.velocity_error(&frozen.discrete().x, &exact).0;
        assert!(frozen_error > 1e-3, "frozen-material error {frozen_error}");
        assert_eq!(previous, vec![0.0; system.n_u]);
    }
}

#[test]
fn material_rheology_refuses_source_density_and_nonconvergence_without_publication() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let reference = transport(2.0, 1.0);
    let flow =
        MaterialFluxSystem::new(&system, &reference, MaterialFluxNumerics::default()).unwrap();
    let previous = vec![0.0; system.n_u];
    let boundary = |point: [f64; 2]| [0.3 * point[1], 0.0];
    let refusal = MaterialError::State {
        what: "outside endpoint material support".into(),
    };
    let mut calls = 0;
    let result =
        flow.bdf1_rheology_step(&zero_force, &boundary, &previous, 0.1, 30, 1e-9, |_, _| {
            calls += 1;
            Err::<FluidTransportState, _>(refusal.clone())
        });
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert_eq!(calls, 1, "no further sources are queried after refusal");

    calls = 0;
    let result = flow.bdf1_rheology_step(
        &zero_force,
        &boundary,
        &previous,
        0.1,
        30,
        1e-9,
        |point, _| {
            calls += 1;
            if point.shear_rate_per_s > 1e-6 {
                Err(refusal.clone())
            } else {
                Ok(reference)
            }
        },
    );
    assert!(matches!(result, Err(MaterialFluxError::Material(error)) if error == refusal));
    assert!(
        calls > mesh.tris.len(),
        "refusal must occur after the initial solve"
    );
    assert!(
        calls <= 2 * mesh.tris.len(),
        "refused endpoint must abort promptly"
    );

    let result =
        flow.bdf1_rheology_step(&zero_force, &boundary, &previous, 0.1, 30, 1e-9, |_, _| {
            Ok(transport(3.0, 1.0))
        });
    assert!(matches!(result, Err(MaterialFluxError::Invalid(_))));
    let result =
        flow.bdf1_rheology_step(&zero_force, &boundary, &previous, 0.1, 1, 1e-9, |_, _| {
            Ok(reference)
        });
    assert!(matches!(
        result,
        Err(MaterialFluxError::Rheology(
            FluxRheologyError::NonlinearConvergence { iterations: 1, .. }
        ))
    ));
    assert_eq!(previous, vec![0.0; system.n_u]);
}
