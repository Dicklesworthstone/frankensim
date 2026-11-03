//! Independent physical population moments of the shared initial strain owner.

use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

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
