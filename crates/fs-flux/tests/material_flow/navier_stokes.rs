//! Source density scales convection, physical traction and energy together.

use super::*;
use fs_flux::{MaterialNavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy};

#[path = "navier_stokes_rheology.rs"]
mod rheology;
#[path = "reservoir.rs"]
mod reservoir;

#[test]
fn nonlinear_mixed_flow_keeps_shared_gas_density_and_viscous_traction_in_si() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let spec = GasSpec::dry_air_ussa1976();
    let dt = 0.025;
    let target = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let previous = vec![0.0; system.n_u];
    let mut energy_per_density = Vec::new();
    for source_pressure in [101_325.0, 202_650.0] {
        let gas = GasState::try_new(&spec, 293.15, source_pressure).unwrap();
        let owner =
            MaterialFluxSystem::new(&system, &gas, MaterialFluxNumerics::default()).unwrap();
        let density = owner.transport().density_kg_m3();
        let mu = owner.transport().dynamic_viscosity_pa_s();
        let force = |p: [f64; 2]| {
            let u = target(p);
            // grad(u)*u = [.04*x+.11, .04*y-.02].
            [
                density * (u[0] / dt + 0.04 * p[0] + 0.11) + 2.0,
                density * (u[1] / dt + 0.04 * p[1] - 0.02) - 3.0,
            ]
        };
        let traction = |p: [f64; 2], n: [f64; 2]| {
            let pressure = 7.0 + 2.0 * p[0] - 3.0 * p[1];
            [
                (0.4 * mu - pressure) * n[0] + 0.3 * mu * n[1],
                0.3 * mu * n[0] + (-0.4 * mu - pressure) * n[1],
            ]
        };
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e))
            .map(|(edge, e)| {
                if e.normal[0].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &target,
                    }
                }
            })
            .collect();
        let step = owner
            .bdf1_navier_stokes_mixed_step(
                &force,
                &boundaries,
                &previous,
                dt,
                NavierStokesControls {
                    relative_tolerance: 1e-10,
                    ..NavierStokesControls::default()
                },
            )
            .unwrap();
        let (error, divergence) = system.velocity_error(&step.solution.discrete().x, &target);
        assert!(error < 2e-9 && divergence < 2e-9, "{error}, {divergence}");
        for (actual, p) in step.solution.pressure_pa().iter().zip(&mesh.centroids) {
            near(*actual, 7.0 + 2.0 * p[0] - 3.0 * p[1], 3e-8);
        }
        assert_eq!(step.solution.transport(), owner.transport());
        assert_eq!(
            step.solution.pressure_reference(),
            StokesPressureReference::TractionDatum
        );
        assert!(step.iterations > 1 && step.nonlinear_residual < 1e-10);
        assert!(
            step.open_boundary_policies
                .iter()
                .all(|(_, p)| *p == OpenBoundaryPolicy::NaturalTraction)
        );
        energy_per_density.push(step.solution.kinetic_energy_j_per_m() / density);
    }
    // Independent integral on [0,1]² of |u|²/2 for the affine target.
    let x_energy = 0.2_f64.powi(2) / 3.0
        + 0.3_f64.powi(2) / 3.0
        + 0.4_f64.powi(2)
        + 0.2 * 0.3 / 2.0
        + 0.2 * 0.4
        + 0.3 * 0.4;
    let y_energy = 0.2_f64.powi(2) / 3.0 - 0.2 * 0.1 + 0.1_f64.powi(2);
    for value in energy_per_density {
        near(value, 0.5 * (x_energy + y_energy), 2e-9);
    }
    assert!(previous.iter().all(|value| *value == 0.0));
}
