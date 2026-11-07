//! Wall-aware pathline history on the original spatial force and wall owners.

use super::*;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::lattice3::wall_loads::WallMomentumProbe3;

const WALLS: [[f64; 3]; 2] = [[0.01, -0.004, 0.0], [0.018, 0.003, 0.0]];

fn noise(cell: usize, salt: u64) -> f64 {
    let mut h = (cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

fn channel_flow<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) -> ForcedGrid3<L, Q> {
    let incoming = force_field(0.0);
    let grid = Grid3::<L, Q>::try_new(
        [5; 3],
        collision,
        [0.0; 3],
        [true, true, false],
        GridBudget3 {
            maximum_cells: 125,
            maximum_mach: 0.2,
        },
        |_| false,
        |[x, y, z]| {
            let cell = (z * 5 + y) * 5 + x;
            let h = (z as f64 + 0.5) / 5.0;
            let envelope = 4.0 * h * (1.0 - h);
            let velocity = core::array::from_fn(|axis| {
                WALLS[0][axis]
                    + h * (WALLS[1][axis] - WALLS[0][axis])
                    + 0.02 * envelope * noise(cell, 31 * axis as u64 + 7)
                    - 0.5 * incoming[cell][axis]
            });
            equilibrium(1.0, velocity)
        },
    )
    .unwrap();
    ForcedGrid3::try_new(grid, cx, incoming).unwrap()
}

fn drives<L: Lattice3<Q>, const Q: usize>(probe: &WallMomentumProbe3<L, Q>) -> Vec<[f64; 3]> {
    probe
        .wall_links()
        .map(|link| WALLS[usize::from(link.exterior_faces[5])])
        .collect()
}

fn bits<const Q: usize>(populations: &[[f64; Q]]) -> Vec<u64> {
    populations
        .iter()
        .flatten()
        .map(|value| value.to_bits())
        .collect()
}

fn channel_reference<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) where
    L::Collision: core::fmt::Debug + PartialEq,
{
    with_cx(|cx| {
        let mut native = channel_flow::<L, Q>(cx, collision, equilibrium);
        let mut independent = channel_flow::<L, Q>(cx, collision, equilibrium);
        let mut pooled = channel_flow::<L, Q>(cx, collision, equilibrium);
        let probe = native
            .wall_momentum_probe(cx, [0.5, -0.3, 0.2], |_| true)
            .unwrap();
        let velocities = drives(&probe);
        let mut work = native.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut pooled_work = pooled.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
        assert!(
            native
                .channel_lagrangian_sgs_workspace(work.workspace_bytes() - 1)
                .is_err()
        );
        let mut model = LagrangianDynamicSgs::new_channel([5; 3], [1.0; 3], 2, usize::MAX).unwrap();
        let molecular: Vec<_> = (0..125)
            .map(|cell| L::with_shear_relaxation(collision, 0.7 + 0.0005 * cell as f64).unwrap())
            .collect();
        let original = molecular.clone();
        let pool = TilePool::for_host(2, 501);
        let gate = CancelGate::new_clock_free();
        let mut saw_modeled_transport = false;
        let mut saw_wall_work = false;
        for step in 0..4 {
            let macros: Vec<_> = (0..125)
                .map(|cell| independent.macroscopic(cell).unwrap())
                .collect();
            let density: Vec<_> = macros.iter().map(|m| m.0).collect();
            let velocity: Vec<_> = macros.iter().map(|m| m.1).collect();
            let advance = step % 2 == 0;
            let next = if advance {
                force_field(step as f64 + 1.0)
            } else {
                independent.forces().to_vec()
            };
            let trial = model
                .prepare_with_walls(config(), step as f64, &density, &velocity, WALLS, || false)
                .unwrap();
            if step == 0 {
                // A distinct endpoint force must not reinterpret the incoming
                // velocity used for either the fit or pathline departure.
                let wrong_velocity: Vec<_> = (0..125)
                    .map(|cell| {
                        core::array::from_fn(|axis| {
                            velocity[cell][axis]
                                + 0.5 * (next[cell][axis] - independent.forces()[cell][axis])
                                    / density[cell]
                        })
                    })
                    .collect();
                let mut wrong =
                    LagrangianDynamicSgs::new_channel([5; 3], [1.0; 3], 2, usize::MAX).unwrap();
                let wrong_trial = wrong
                    .prepare_with_walls(config(), 0.0, &density, &wrong_velocity, WALLS, || false)
                    .unwrap();
                assert!(
                    trial
                        .viscosities()
                        .iter()
                        .zip(wrong_trial.viscosities())
                        .any(|(a, b)| (a - b).abs() > 1e-12)
                );
                let mut stationary =
                    LagrangianDynamicSgs::new_channel([5; 3], [1.0; 3], 2, usize::MAX).unwrap();
                let stationary_trial = stationary
                    .prepare_with_walls(config(), 0.0, &density, &velocity, [[0.0; 3]; 2], || false)
                    .unwrap();
                assert!(
                    trial
                        .viscosities()
                        .iter()
                        .zip(stationary_trial.viscosities())
                        .any(|(a, b)| (a - b).abs() > 1e-12)
                );
            }
            let rates: Vec<_> = trial
                .viscosities()
                .iter()
                .enumerate()
                .map(|(cell, &eddy)| {
                    let source = molecular[cell];
                    if eddy == 0.0 {
                        source
                    } else {
                        let nu = L::molecular_viscosity(source).unwrap();
                        L::with_shear_relaxation(source, 0.5 + (nu + eddy) / L::CS2).unwrap()
                    }
                })
                .collect();
            let expected_fit = trial.report();
            let expected_wall = if advance {
                independent
                    .advance_serial_with_moving_wall_momentum(
                        cx,
                        &probe,
                        &velocities,
                        &rates,
                        next.clone(),
                    )
                    .unwrap()
            } else {
                independent
                    .step_serial_with_moving_wall_momentum(cx, &probe, &velocities, &rates)
                    .unwrap()
            };
            let (serial, parallel) = if advance {
                (
                    native
                        .advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                            cx,
                            config(),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut work,
                            next.clone(),
                        )
                        .unwrap(),
                    pooled
                        .advance_pooled_with_channel_lagrangian_sgs_and_moving_walls(
                            &pool,
                            &gate,
                            config(),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut pooled_work,
                            next.clone(),
                        )
                        .unwrap(),
                )
            } else {
                (
                    native
                        .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                            cx,
                            config(),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut work,
                        )
                        .unwrap(),
                    pooled
                        .step_pooled_with_channel_lagrangian_sgs_and_moving_walls(
                            &pool,
                            &gate,
                            config(),
                            Some(&molecular),
                            &probe,
                            &velocities,
                            &mut pooled_work,
                        )
                        .unwrap(),
                )
            };
            trial.commit();
            assert_eq!(serial.0.fit, expected_fit);
            assert_eq!(serial.1, expected_wall);
            assert_eq!(serial, parallel);
            saw_modeled_transport |= serial.0.fit.maximum_eddy_viscosity > 0.0;
            saw_wall_work |= serial.1.wall_work.abs() > 1e-10;
            assert_eq!(work.collisions(), rates);
            assert_eq!(work.collisions(), pooled_work.collisions());
            assert_eq!(work.accepted_fit_time(), Some(step as f64));
            assert_eq!(pooled_work.accepted_fit_time(), Some(step as f64));
            assert_eq!(bits(native.populations()), bits(independent.populations()));
            assert_eq!(bits(native.populations()), bits(pooled.populations()));
            assert_eq!(native.forces(), next);
            assert_eq!(native.forces(), pooled.forces());
            assert_eq!(native.completed_steps(), step + 1);
            assert_eq!(native.completed_steps(), pooled.completed_steps());
        }
        assert!(saw_modeled_transport);
        assert!(saw_wall_work);
        assert_eq!(molecular, original);
    });
}

#[test]
fn d3q19_moving_channel_force_epochs_match_independent_pathline_history() {
    channel_reference::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
}

#[test]
fn d3q27_moving_channel_preserves_cumulants_wall_receipts_and_worker_identity() {
    channel_reference::<D3q27, 27>(cumulant(), eq27);
}

fn channel_rollback<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) {
    with_cx(|cx| {
        let mut flow = channel_flow::<L, Q>(cx, collision, equilibrium);
        let mut control = channel_flow::<L, Q>(cx, collision, equilibrium);
        let probe = flow.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let velocities = drives(&probe);
        let mut work = flow.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
        let mut control_work = control
            .channel_lagrangian_sgs_workspace(usize::MAX)
            .unwrap();
        flow.advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
            cx,
            config(),
            None,
            &probe,
            &velocities,
            &mut work,
            force_field(1.0),
        )
        .unwrap();
        control
            .advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocities,
                &mut control_work,
                force_field(1.0),
            )
            .unwrap();
        let populations = bits(flow.populations());
        let forces = flow.forces().to_vec();
        assert!(
            flow.advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocities,
                &mut work,
                vec![[1.0; 3]; 125],
            )
            .is_err()
        );
        assert_eq!(bits(flow.populations()), populations);
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        let gate = CancelGate::new_clock_free();
        let refusing = RefuseAfterStreaming {
            pool: TilePool::for_host(2, 501),
            calls: AtomicUsize::new(0),
        };
        assert!(
            flow.advance_pooled_with_channel_lagrangian_sgs_and_moving_walls(
                &refusing,
                &gate,
                config(),
                None,
                &probe,
                &velocities,
                &mut work,
                force_field(2.0),
            )
            .is_err()
        );
        assert_eq!(refusing.calls.load(Ordering::Acquire), 2);
        assert_eq!(bits(flow.populations()), populations);
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(work.accepted_fit_time(), Some(0.0));
        let retry = flow
            .advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocities,
                &mut work,
                force_field(2.0),
            )
            .unwrap();
        let expected = control
            .advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                cx,
                config(),
                None,
                &probe,
                &velocities,
                &mut control_work,
                force_field(2.0),
            )
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(bits(flow.populations()), bits(control.populations()));
        assert_eq!(flow.forces(), control.forces());
        assert_eq!(flow.completed_steps(), 2);
        assert_eq!(work.accepted_fit_time(), Some(1.0));
    });
}

#[test]
fn refused_force_or_post_stream_cancellation_preserves_channel_history_for_retry() {
    channel_rollback::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    channel_rollback::<D3q27, 27>(cumulant(), eq27);
}

fn zero_drive<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    equilibrium: fn(f64, [f64; 3]) -> [f64; Q],
) {
    with_cx(|cx| {
        let mut stationary = channel_flow::<L, Q>(cx, collision, equilibrium);
        let mut moving = channel_flow::<L, Q>(cx, collision, equilibrium);
        let probe = moving.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap();
        let zeros = vec![[0.0; 3]; probe.measured_links()];
        let mut stationary_work = stationary
            .channel_lagrangian_sgs_workspace(usize::MAX)
            .unwrap();
        let mut moving_work = moving.channel_lagrangian_sgs_workspace(usize::MAX).unwrap();
        for step in 0..3 {
            let (expected, actual) = if step % 2 == 0 {
                let next = force_field(step as f64 + 1.0);
                (
                    stationary
                        .advance_serial_with_lagrangian_sgs(
                            cx,
                            config(),
                            None,
                            &mut stationary_work,
                            next.clone(),
                        )
                        .unwrap(),
                    moving
                        .advance_serial_with_channel_lagrangian_sgs_and_moving_walls(
                            cx,
                            config(),
                            None,
                            &probe,
                            &zeros,
                            &mut moving_work,
                            next,
                        )
                        .unwrap(),
                )
            } else {
                (
                    stationary
                        .step_serial_with_lagrangian_sgs(cx, config(), None, &mut stationary_work)
                        .unwrap(),
                    moving
                        .step_serial_with_channel_lagrangian_sgs_and_moving_walls(
                            cx,
                            config(),
                            None,
                            &probe,
                            &zeros,
                            &mut moving_work,
                        )
                        .unwrap(),
                )
            };
            assert_eq!(expected, actual.0);
            assert_eq!(actual.1.wall_work, 0.0);
            assert_eq!(bits(stationary.populations()), bits(moving.populations()));
            assert_eq!(stationary.forces(), moving.forces());
            assert_eq!(stationary.completed_steps(), moving.completed_steps());
            assert_eq!(
                stationary_work.accepted_fit_time(),
                moving_work.accepted_fit_time()
            );
        }
    });
}

#[test]
fn zero_drive_keeps_stationary_channel_force_and_history_bits() {
    zero_drive::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, equilibrium3);
    zero_drive::<D3q27, 27>(cumulant(), eq27);
}
