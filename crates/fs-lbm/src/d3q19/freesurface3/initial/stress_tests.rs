//! Independent physical population moments of the shared initial strain owner.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

#[test]
fn g1_initial_density_closes_authored_mass_without_a_force_impulse_or_conversion() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 28,
                kernel_id: 69,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        for scale in [2.0_f64.powi(-40), 1.0, 2.0_f64.powi(40)] {
            for velocity in [[0.0; 3], [0.025, -0.012, 0.018]] {
                for force in [[0.0; 3], [0.0002, -0.0001, 0.0003].map(|f| f * scale)] {
                    let initial = [
                        FreeSurfaceInitialCell3 {
                            coordinates: [1, 1, 1],
                            fill_fraction: 1.0,
                            // Hydrostatic density whose raw equilibrium sum
                            // falls below its authored value by two ulps.
                            density: 1.0015000000000001 * scale,
                            velocity,
                        },
                        FreeSurfaceInitialCell3 {
                            coordinates: [2, 1, 1],
                            fill_fraction: 0.35,
                            density: 0.99775 * scale,
                            velocity,
                        },
                    ];
                    let mut flow = FreeSurface3::from_cells(
                        [4; 3],
                        0.8,
                        force,
                        0.0,
                        ContactModel3::Neutral,
                        &[],
                        &initial,
                        &cx,
                    )
                    .unwrap();
                    for state in initial {
                        let [x, y, z] = state.coordinates;
                        let (slot, lane) = flow
                            .grid
                            .resolve_source(x as i64, y as i64, z as i64)
                            .unwrap();
                        let key = flow.grid.active_keys()[slot];
                        let populations = flow.grid.populations(slot, lane);
                        let rho: f64 = populations.iter().sum();
                        let mass = flow.tiles[&key].mass[lane];
                        assert_eq!(rho.to_bits(), state.density.to_bits());
                        assert_eq!(
                            mass.to_bits(),
                            (state.fill_fraction * state.density).to_bits()
                        );
                        assert!(mass <= rho);
                        for axis in 0..3 {
                            let raw: f64 = populations
                                .iter()
                                .zip(E3)
                                .map(|(f, e)| f * f64::from([e.0, e.1, e.2][axis]))
                                .sum();
                            let physical = raw + 0.5 * force[axis];
                            assert!(
                                (physical - state.density * velocity[axis]).abs()
                                    < 64.0 * f64::EPSILON * state.density
                            );
                        }
                    }
                    let bits = flow.grid.state_bits();
                    let mass = flow.mass_snapshot();
                    flow.apply_conversions(None).unwrap();
                    assert_eq!(flow.grid.state_bits(), bits);
                    assert_eq!(flow.mass_snapshot(), mass);
                    assert_eq!(flow.conversions, ConversionStats3::default());
                    assert_eq!(flow.conversion_momentum_report().unwrap().routes, 0);
                }
            }
        }
    });
}

#[test]
fn g1_initial_strain_preserves_forced_momentum_and_recovers_molecular_stress() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 28,
                kernel_id: 69,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let coordinate = [1, 1, 1];
        let force = [0.0002, -0.0001, 0.0003];
        let tau = 0.83;
        let mut flow = FreeSurface3::from_cells(
            [4; 3],
            tau,
            force,
            0.0,
            ContactModel3::Neutral,
            &[],
            &[FreeSurfaceInitialCell3 {
                coordinates: coordinate,
                fill_fraction: 0.6,
                density: 1.2,
                velocity: [0.03, -0.02, 0.01],
            }],
            &cx,
        )
        .unwrap();
        let key = morton3(0, 0, 0);
        let slot = flow.grid.slot_of(key).unwrap();
        let lane = 1 + TILE * (1 + TILE);
        let original = flow.grid.populations(slot, lane);
        let before = flow.sample_at(1, 1, 1).unwrap();
        let mass = flow.ledger_mass();
        flow.initialize_strain_at(coordinate, [0.0; 6], tau)
            .unwrap();
        assert_eq!(
            flow.grid.populations(slot, lane).map(f64::to_bits),
            original.map(f64::to_bits)
        );

        let strain = [0.002, -0.001, 0.0005, 0.0007, -0.0004, 0.0009];
        flow.initialize_strain_at(coordinate, strain, tau).unwrap();
        let populations = flow.grid.populations(slot, lane);
        let after = flow.sample_at(1, 1, 1).unwrap();
        assert_eq!(before.density.to_bits(), after.density.to_bits());
        assert_eq!(mass.to_bits(), flow.ledger_mass().to_bits());
        for axis in 0..3 {
            assert!((before.velocity[axis] - after.velocity[axis]).abs() < 4e-16);
        }
        let strain_tensor = [
            [strain[0], strain[3], strain[4]],
            [strain[3], strain[1], strain[5]],
            [strain[4], strain[5], strain[2]],
        ];
        for a in 0..3 {
            for b in 0..3 {
                let second = populations
                    .iter()
                    .enumerate()
                    .map(|(q, &f)| {
                        let e = [E3[q].0, E3[q].1, E3[q].2];
                        f * f64::from(e[a]) * f64::from(e[b])
                    })
                    .sum::<f64>();
                // Independent physical second moment: restore the half force,
                // remove equilibrium pressure and convective momentum flux.
                let pi = second
                    + 0.5 * (after.velocity[a] * force[b] + force[a] * after.velocity[b])
                    - after.density * after.velocity[a] * after.velocity[b]
                    - if a == b { after.density / 3.0 } else { 0.0 };
                let stress = -(1.0 - 0.5 / tau) * pi;
                let expected = 2.0 * after.density * (tau - 0.5) / 3.0 * strain_tensor[a][b];
                assert!(
                    (stress - expected).abs() < 2e-16,
                    "({a},{b}): {stress} != {expected}"
                );
            }
        }
        // A nonrepresentable positive-population state cannot partially modify
        // even this private cell initializer, and the outer constructor returns
        // no candidate on the same failure.
        let accepted = populations.map(f64::to_bits);
        assert!(
            flow.initialize_strain_at(coordinate, [0.0, 0.0, 0.0, 100.0, 0.0, 0.0], tau)
                .is_err()
        );
        assert_eq!(
            flow.grid.populations(slot, lane).map(f64::to_bits),
            accepted
        );
        assert_eq!(mass.to_bits(), flow.ledger_mass().to_bits());
    });
}

#[test]
fn g1_initial_strain_retains_independent_shear_and_bulk_stress() {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey { seed: 28, kernel_id: 69, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let coordinate = [1, 1, 1];
        let force = [0.0002, -0.0001, 0.0003];
        // A spatial material solve can supply a local shear relaxation that
        // differs from the reference collision law. Neither changes tau_bulk.
        let local_tau = 0.83;
        for bulk_tau in [0.58, 0.83, 1.37] {
            for strain in [
                [0.0015, 0.0015, 0.0015, 0.0, 0.0, 0.0],
                [0.002, -0.001, 0.0005, 0.0007, -0.0004, 0.0009],
            ] {
                let mut flow = FreeSurface3::from_cells(
                    [4; 3], 0.91, force, 0.0, ContactModel3::Neutral, &[],
                    &[FreeSurfaceInitialCell3 {
                        coordinates: coordinate,
                        fill_fraction: 0.6,
                        density: 1.2,
                        velocity: [0.03, -0.02, 0.01],
                    }], &cx,
                ).unwrap();
                flow.grid.set_collision_model(super::super::super::CollisionModel3::CentralMomentBulk {
                    shear_rate: 1.0 / 0.91,
                    bulk_rate: 1.0 / bulk_tau,
                    higher_order_rate: 1.2,
                }).unwrap();
                let before = flow.sample_at(1, 1, 1).unwrap();
                let mass = flow.ledger_mass();
                flow.initialize_strain_at(coordinate, strain, local_tau).unwrap();
                let after = flow.sample_at(1, 1, 1).unwrap();
                assert_eq!(before.density.to_bits(), after.density.to_bits());
                assert_eq!(mass.to_bits(), flow.ledger_mass().to_bits());
                for axis in 0..3 {
                    assert!((before.velocity[axis] - after.velocity[axis]).abs() < 4e-16);
                }
                let slot = flow.grid.slot_of(morton3(0, 0, 0)).unwrap();
                let lane = 1 + TILE * (1 + TILE);
                let populations = flow.grid.populations(slot, lane);
                let mut pi = [[0.0; 3]; 3];
                for a in 0..3 {
                    for b in 0..3 {
                        pi[a][b] = populations.iter().enumerate().map(|(q, &f)| {
                            let e = [E3[q].0, E3[q].1, E3[q].2];
                            f * f64::from(e[a]) * f64::from(e[b])
                        }).sum::<f64>()
                            + 0.5 * (after.velocity[a] * force[b] + force[a] * after.velocity[b])
                            - after.density * after.velocity[a] * after.velocity[b]
                            - if a == b { after.density / 3.0 } else { 0.0 };
                    }
                }
                let mean_pi = (pi[0][0] + pi[1][1] + pi[2][2]) / 3.0;
                let theta = strain[0] + strain[1] + strain[2];
                let tensor = [
                    [strain[0], strain[3], strain[4]],
                    [strain[3], strain[1], strain[5]],
                    [strain[4], strain[5], strain[2]],
                ];
                let mu = after.density * (local_tau - 0.5) / 3.0;
                let zeta = 2.0 * after.density * (bulk_tau - 0.5) / 9.0;
                for a in 0..3 {
                    for b in 0..3 {
                        let spherical = if a == b { mean_pi } else { 0.0 };
                        let stress = -(1.0 - 0.5 / local_tau) * (pi[a][b] - spherical)
                            - (1.0 - 0.5 / bulk_tau) * spherical;
                        let expected = 2.0 * mu
                            * (tensor[a][b] - if a == b { theta / 3.0 } else { 0.0 })
                            + if a == b { zeta * theta } else { 0.0 };
                        assert!((stress - expected).abs() < 8e-16,
                            "tau_bulk={bulk_tau} ({a},{b}): {stress} != {expected}");
                    }
                }
                let accepted = populations.map(f64::to_bits);
                assert!(flow.initialize_strain_at(coordinate,
                    [100.0, 100.0, 100.0, 0.0, 0.0, 0.0], local_tau).is_err());
                assert_eq!(flow.grid.populations(slot, lane).map(f64::to_bits), accepted);
                assert_eq!(mass.to_bits(), flow.ledger_mass().to_bits());
            }
        }
    });
}
