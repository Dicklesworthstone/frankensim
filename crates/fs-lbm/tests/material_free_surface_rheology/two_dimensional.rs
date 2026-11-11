//! Dense source-rheology regressions reuse the existing original-source fixture.

#[path = "two_dimensional/bulk.rs"]
mod bulk;

use super::{DT, DX, card, law, near, solve};
use fs_exec::CancelGate;
use fs_lbm::core2::CollisionModel2;
use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::material::free_surface::rheology::{
    FreeSurfaceRheologyError, MaterialFreeSurfaceRheology2,
};
use fs_lbm::rheology::{Rheology, update_tau};
use fs_lbm::{Cell, Grid};
use fs_material::fluid_transport::viscous::DeformationRate;

// Independent quadrature: the native D2Q9 directions are crate-private.
const E: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

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
    match (left.grid.collision, right.grid.collision) {
        (
            CollisionModel2::CentralMomentBulk { bulk_tau: left },
            CollisionModel2::CentralMomentBulk { bulk_tau: right },
        ) => {
            assert_eq!(left.to_bits(), right.to_bits());
        }
        (left, right) => assert_eq!(left, right),
    }
    let trace_bits = |grid: &Grid| {
        grid.central_moment_bulk_tau
            .as_ref()
            .map(|field| field.iter().map(|tau| tau.to_bits()).collect::<Vec<_>>())
    };
    assert_eq!(trace_bits(&left.grid), trace_bits(&right.grid));
    assert_eq!(left.grid.g, right.grid.g);
    assert_eq!(left.grid.fext, right.grid.fext);
    assert_eq!(left.mass, right.mass);
    assert_eq!(left.carry, right.carry);
    assert_eq!(left.conversions.to_fluid, right.conversions.to_fluid);
    assert_eq!(left.conversions.to_gas, right.conversions.to_gas);
    assert_eq!(
        left.conversions.gas_to_interface,
        right.conversions.gas_to_interface
    );
    assert_eq!(
        left.conversions.fluid_to_interface,
        right.conversions.fluid_to_interface
    );
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
        grid.rebind_gravity([1e-6, -2e-6], &mut Vec::new()).unwrap();
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
fn planar_extension_queries_the_physical_material_tensor_in_the_guo_frame() {
    let source = card(0.01, 100.0, 0.001);
    let binding = law(&source);
    for strain in [[[0.5, 0.0], [0.0, 0.5]], [[0.5, 0.3], [0.3, -0.2]]] {
        let expected_rate =
            DeformationRate::from_tensor([strain[0][0], strain[1][1], 0.0, strain[0][1], 0.0, 0.0])
                .unwrap()
                .shear_rate_per_s();
        let viscosity = 0.1 + 0.01 * expected_rate;
        let tau = 0.5 + 3.0 * viscosity * DT / (1000.0 * DX * DX);
        for force in [[0.0; 2], [0.0003, -0.0002]] {
            let mut rheology =
                MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve())
                    .unwrap();
            let mut grid = Grid::uniform(5, 4, 1.6);
            let cell = grid.idx(2, 2);
            let u = [0.01, -0.015];
            let equilibrium = fs_lbm::equilibrium(1.0, u[0], u[1]);
            // Prescribe original physical moments with the Hermite projection,
            // including both Guo half-force terms; rho=1, D_lattice=D_SI*dt.
            grid.f[cell] = std::array::from_fn(|q| {
                let e = [f64::from(E[q].0), f64::from(E[q].1)];
                let momentum = -0.5 * (e[0] * force[0] + e[1] * force[1]);
                let mut contraction = 0.0;
                for a in 0..2 {
                    for b in 0..2 {
                        let stress = -2.0 * tau * strain[a][b] * DT / 3.0
                            - 0.5 * (u[a] * force[b] + force[a] * u[b]);
                        contraction +=
                            (e[a] * e[b] - if a == b { 1.0 / 3.0 } else { 0.0 }) * stress;
                    }
                }
                equilibrium[q] + W[q] * (3.0 * momentum + 4.5 * contraction)
            });
            grid.fext[cell] = force;
            let mut liquid = FreeSurface::new(
                grid,
                rheology.scaling().surface_tension_lattice(),
                ContactModel::Neutral,
            );
            rheology
                .step(&mut liquid, &CancelGate::new_clock_free())
                .unwrap();
            let sample = rheology.rheology_at([2, 2]).unwrap();
            near(sample.shear_rate_per_s, expected_rate, 2e-10);
            near(sample.collision_tau, tau, 2e-11);
            let material = rheology.material_at([2, 2]).unwrap();
            near(
                material.property("dynamic_viscosity").unwrap().value_si(),
                viscosity,
                2e-11,
            );
            assert_eq!(material.card_identity(), source.content_hash());
            assert!(
                material
                    .query_point()
                    .contains(&("rate".into(), sample.shear_rate_per_s))
            );
        }
    }
}

#[test]
fn central_trace_fields_resolve_the_physical_source_and_reuse_native_populations() {
    let source = card(0.01, 100.0, 0.001);
    let binding = law(&source);
    let gate = CancelGate::new_clock_free();
    for strain in [
        [[0.6, 0.0], [0.0, 0.6]],
        [[0.0, 0.4], [0.4, 0.0]],
        [[1.1, 0.2], [0.2, -0.4]],
    ] {
        let expected_rate =
            DeformationRate::from_tensor([strain[0][0], strain[1][1], 0.0, strain[0][1], 0.0, 0.0])
                .unwrap()
                .shear_rate_per_s();
        let mu = 0.1 + 0.01 * expected_rate;
        let shear_tau = 0.5 + 3.0 * mu * DT / (1000.0 * DX * DX);
        for rho in [0.8, 1.2] {
            for trace_mode in 0..3 {
                for force in [[0.0; 2], [0.0003, -0.0002]] {
                    let mut rheology = MaterialFreeSurfaceRheology2::try_new(
                        &binding,
                        [5, 4],
                        DX,
                        DT,
                        0.05,
                        solve(),
                    )
                    .unwrap();
                    let mut grid = Grid::uniform(5, 4, 1.6);
                    grid.collision = if trace_mode == 0 {
                        CollisionModel2::CentralMoment
                    } else {
                        CollisionModel2::CentralMomentBulk { bulk_tau: 1.7 }
                    };
                    let trace: Vec<_> = (0..20)
                        .map(|cell| match trace_mode {
                            0 => 1.0,
                            1 => 1.7,
                            _ => [0.65, 1.25, 1.8][cell % 3],
                        })
                        .collect();
                    if trace_mode == 2 {
                        grid.set_central_moment_bulk_field(&trace).unwrap();
                    }
                    let u = [0.01, -0.015];
                    let equilibrium = fs_lbm::equilibrium(rho, u[0], u[1]);
                    let mean = 0.5 * (strain[0][0] + strain[1][1]);
                    for cell in 0..20 {
                        // Independent original moments: trace and planar
                        // deviator relax separately, both in the Guo frame.
                        grid.f[cell] = std::array::from_fn(|q| {
                            let e = [f64::from(E[q].0), f64::from(E[q].1)];
                            let momentum = -0.5 * (e[0] * force[0] + e[1] * force[1]);
                            let mut contraction = 0.0;
                            for a in 0..2 {
                                for b in 0..2 {
                                    let isotropic = if a == b { mean } else { 0.0 };
                                    let stress = -2.0 * rho * DT / 3.0
                                        * (shear_tau * (strain[a][b] - isotropic)
                                            + trace[cell] * isotropic)
                                        - 0.5 * (u[a] * force[b] + force[a] * u[b]);
                                    contraction += (e[a] * e[b]
                                        - if a == b { 1.0 / 3.0 } else { 0.0 })
                                        * stress;
                                }
                            }
                            equilibrium[q] + W[q] * (3.0 * momentum + 4.5 * contraction)
                        });
                        grid.fext[cell] = force;
                    }
                    let mut actual = FreeSurface::new(
                        grid,
                        rheology.scaling().surface_tension_lattice(),
                        ContactModel::Neutral,
                    );
                    let mut direct = actual.clone();
                    direct.grid.tau.fill(shear_tau);
                    let expected = direct.step_with_momentum(&gate).unwrap();
                    let report = rheology.step(&mut actual, &gate).unwrap();
                    assert_eq!(report, expected);
                    assert_eq!(actual.grid.flags, direct.grid.flags);
                    assert_eq!(
                        actual.grid.central_moment_bulk_tau,
                        direct.grid.central_moment_bulk_tau
                    );
                    for cell in 0..20 {
                        let position = [cell % 5, cell / 5];
                        let sample = rheology.rheology_at(position).unwrap();
                        near(sample.shear_rate_per_s, expected_rate, 2e-9);
                        near(sample.collision_tau, shear_tau, 2e-11);
                        let material = rheology.material_at(position).unwrap();
                        near(
                            material.property("dynamic_viscosity").unwrap().value_si(),
                            mu,
                            2e-11,
                        );
                        assert_eq!(material.card_identity(), source.content_hash());
                        assert!(
                            material
                                .query_point()
                                .contains(&("rate".into(), sample.shear_rate_per_s),)
                        );
                        assert_eq!(
                            material
                                .property("dynamic_viscosity")
                                .unwrap()
                                .answer()
                                .receipt
                                .selected,
                            binding
                                .material()
                                .property("dynamic_viscosity")
                                .unwrap()
                                .answer()
                                .receipt
                                .selected,
                        );
                        for q in 0..9 {
                            near(actual.grid.f[cell][q], direct.grid.f[cell][q], 2e-12);
                        }
                        near(actual.mass[cell], direct.mass[cell], 2e-12);
                    }
                    near(actual.ledger_mass(), direct.ledger_mass(), 1e-11);
                }
            }
        }
    }
}

#[test]
fn g0_late_source_refusal_preserves_flow_and_the_previous_collision_receipts() {
    let source = card(0.01, 1.0, 0.001);
    let binding = law(&source);
    let mut rheology =
        MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve()).unwrap();
    let mut grid = Grid::uniform(5, 4, 0.8);
    let trace: Vec<_> = (0..20).map(|cell| 0.7 + 0.05 * cell as f64).collect();
    grid.set_central_moment_bulk_field(&trace).unwrap();
    let mut liquid = FreeSurface::new(
        grid,
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

    for bad in [None, Some(0.5), Some(f64::NAN)] {
        let mut invalid = before.clone();
        let field = invalid.grid.central_moment_bulk_tau.as_mut().unwrap();
        if let Some(value) = bad {
            field[19] = value;
        } else {
            field.pop();
        }
        let saved = invalid.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            rheology.step(&mut invalid, &gate)
        }));
        assert!(matches!(
            result.expect("bad trace fields must refuse without panicking"),
            Err(FreeSurfaceRheologyError::Invalid(_)),
        ));
        same_native(&invalid, &saved);
        assert_eq!(*rheology.material_at([0, 0]).unwrap(), receipt);
        assert_eq!(*rheology.rheology_at([0, 0]).unwrap(), sample);
    }

    for bulk_tau in [0.5, f64::NAN] {
        let mut invalid = before.clone();
        invalid.grid.collision = CollisionModel2::CentralMomentBulk { bulk_tau };
        let saved = invalid.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            rheology.step(&mut invalid, &gate)
        }));
        assert!(matches!(
            result.expect("bad trace scalars must refuse without panicking"),
            Err(FreeSurfaceRheologyError::Invalid(_)),
        ));
        same_native(&invalid, &saved);
        assert_eq!(*rheology.material_at([0, 0]).unwrap(), receipt);
        assert_eq!(*rheology.rheology_at([0, 0]).unwrap(), sample);
    }
}
