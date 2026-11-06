//! Dense source-rheology regressions reuse the existing original-source fixture.

use super::{DT, DX, card, law, near, solve};
use fs_exec::CancelGate;
use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::material::free_surface::rheology::{
    FreeSurfaceRheologyError, MaterialFreeSurfaceRheology2,
};
use fs_lbm::rheology::{Rheology, update_tau};
use fs_lbm::{Cell, E, Grid};

fn patch_shear(grid: &mut Grid, cell: usize, epsilon: f64) {
    // Diagonal +/- populations preserve density and momentum and give
    // Pi_xy=4*epsilon. Thus sqrt(2 D:D)*tau = 12*abs(epsilon).
    for (q, (x, y)) in E.iter().copied().enumerate() {
        grid.f[cell][q] += epsilon * f64::from(x * y);
    }
}

fn same_native(left: &FreeSurface, right: &FreeSurface) {
    assert_eq!(left.grid.f, right.grid.f);
    assert_eq!(left.grid.flags, right.grid.flags);
    assert_eq!(left.grid.tau, right.grid.tau);
    assert_eq!(left.grid.g, right.grid.g);
    assert_eq!(left.grid.fext, right.grid.fext);
    assert_eq!(left.mass, right.mass);
    assert_eq!(left.carry, right.carry);
    assert_eq!(left.conversions.to_fluid, right.conversions.to_fluid);
    assert_eq!(left.conversions.to_gas, right.conversions.to_gas);
    assert_eq!(left.conversions.gas_to_interface, right.conversions.gas_to_interface);
    assert_eq!(left.conversions.fluid_to_interface, right.conversions.fluid_to_interface);
}

#[test]
fn g1_constant_source_keeps_native_vof_bits_and_sourced_capillarity() {
    let source = card(0.0, 100.0, 0.002);
    let binding = law(&source);
    let gate = CancelGate::new_clock_free();
    for contact in [ContactModel::Neutral, ContactModel::Wetting] {
        let mut rheology =
            MaterialFreeSurfaceRheology2::try_new(&binding, [8, 8], DX, DT, 0.05, solve()).unwrap();
        let sigma = rheology.scaling().surface_tension_lattice();
        near(sigma, 0.002, 1e-16);
        let mut grid = Grid::uniform(8, 8, 0.8);
        grid.periodic_x = false;
        grid.periodic_y = false;
        for y in 0..8 {
            for x in 0..8 {
                let cell = grid.idx(x, y);
                grid.flags[cell] = if x == 0 || x == 7 || y == 0 || y == 7 {
                    Cell::Wall
                } else if y <= 3 {
                    Cell::Fluid
                } else {
                    Cell::Gas
                };
            }
        }
        grid.rebind_gravity([1e-6, -2e-6]).unwrap();
        let mut actual = FreeSurface::new(grid, sigma, contact);
        let mut direct = actual.clone();
        for _ in 0..4 {
            let stats = update_tau(
                &mut direct.grid,
                Rheology::Newtonian {
                    nu: rheology.scaling().flow().viscosity,
                },
            );
            assert_eq!((stats.floored, stats.capped), (0, 0));
            let expected = direct.step_with_momentum(&gate).unwrap();
            let report = rheology.step(&mut actual, &gate).unwrap();
            assert_eq!(report, expected);
            same_native(&actual, &direct);
        }
        let material = rheology.material_at([2, 2]).unwrap();
        assert_eq!(material.card_identity(), source.content_hash());
        near(
            material.property("surface_tension").unwrap().value_si(),
            0.002,
            0.0,
        );
        near(
            rheology.rheology_at([2, 2]).unwrap().collision_tau,
            0.8,
            1e-15,
        );
    }
}

#[test]
fn g1_source_curve_matches_an_independent_quadratic_root_and_changes_collision() {
    let source = card(0.01, 100.0, 0.001);
    let binding = law(&source);
    let mut rheology =
        MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve()).unwrap();
    let epsilon = 0.0004;
    let mut grid = Grid::uniform(5, 4, 0.8);
    for cell in 0..grid.f.len() {
        patch_shear(&mut grid, cell, epsilon);
    }
    let mut liquid = FreeSurface::new(
        grid,
        rheology.scaling().surface_tension_lattice(),
        ContactModel::Neutral,
    );
    let mut frozen = liquid.clone();
    let gate = CancelGate::new_clock_free();
    let initial_mass = liquid.ledger_mass();
    rheology.step(&mut liquid, &gate).unwrap();
    frozen.step_with_momentum(&gate).unwrap();

    // mu = 0.1*(1 + 0.1*gamma_SI), rho=1000. Eliminating gamma
    // yields tau²-(1/2+a)tau-a*0.1*N/dt=0, independently of solve_local.
    let numerator = 12.0 * epsilon;
    let a = 3.0 * 0.1 * DT / (1000.0 * DX * DX);
    let linear = 0.5 + a;
    let expected_tau = 0.5 * (linear + (linear * linear + 4.0 * a * 0.1 * numerator / DT).sqrt());
    let expected_rate = numerator / (expected_tau * DT);
    let sample = rheology.rheology_at([2, 2]).unwrap();
    near(sample.collision_tau, expected_tau, 2e-11);
    near(sample.shear_rate_per_s, expected_rate, 2e-10);
    let material = rheology.material_at([2, 2]).unwrap();
    near(
        material.property("dynamic_viscosity").unwrap().value_si(),
        0.1 + 0.01 * expected_rate,
        2e-11,
    );
    assert_eq!(material.card_identity(), source.content_hash());
    assert!(material.query_point().contains(&("T".into(), 400.0)));
    assert!(
        material
            .query_point()
            .contains(&("pressure".into(), 100000.0))
    );
    assert!(sample.collision_tau > 0.8);
    assert_ne!(liquid.grid.f, frozen.grid.f);
    near(liquid.ledger_mass(), initial_mass, 1e-11);
}

#[test]
fn g0_late_source_refusal_preserves_flow_and_the_previous_collision_receipts() {
    let source = card(0.01, 1.0, 0.001);
    let binding = law(&source);
    let mut rheology =
        MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve()).unwrap();
    let mut liquid = FreeSurface::new(
        Grid::uniform(5, 4, 0.8),
        rheology.scaling().surface_tension_lattice(),
        ContactModel::Neutral,
    );
    let gate = CancelGate::new_clock_free();
    rheology.step(&mut liquid, &gate).unwrap();
    let receipt = rheology.material_at([0, 0]).unwrap().clone();
    let sample = *rheology.rheology_at([0, 0]).unwrap();
    patch_shear(&mut liquid.grid, 0, 1e-5); // admissible, changes proposed tau
    patch_shear(&mut liquid.grid, 19, 1e-3); // late source-domain refusal
    let before = liquid.clone();
    assert!(matches!(
        rheology.step(&mut liquid, &gate),
        Err(FreeSurfaceRheologyError::Material(_)),
    ));
    same_native(&liquid, &before);
    assert_eq!(*rheology.material_at([0, 0]).unwrap(), receipt);
    assert_eq!(*rheology.rheology_at([0, 0]).unwrap(), sample);

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(rheology.step(&mut liquid, &cancelled).is_err());
    same_native(&liquid, &before);
    assert_eq!(*rheology.material_at([0, 0]).unwrap(), receipt);
    assert_eq!(*rheology.rheology_at([0, 0]).unwrap(), sample);
}
