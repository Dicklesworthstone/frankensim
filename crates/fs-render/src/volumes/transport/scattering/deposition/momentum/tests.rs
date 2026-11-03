//! G1/G3 photon momentum, angular-force and native energy-account regressions.

use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    CellPhaseFunction, MassScatteringSpectrum, PhaseFunction,
};
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 17,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn medium<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    a: f64,
    s: f64,
    phase: PhaseFunction,
) -> ScatteringGrid<'a> {
    let absorption = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], [1.0; 3]),
        MassAbsorptionSpectrum::constant(a).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 64,
        },
    )
    .unwrap();
    ScatteringGrid::try_new(
        cx,
        absorption,
        MassScatteringSpectrum::constant(s).unwrap(),
        phase,
    )
    .unwrap()
}

fn config(histories: u32, collisions: u32) -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: histories,
        max_collisions_per_history: collisions,
        max_total_segments: u64::from(histories) * u64::from(collisions) * 4,
        max_deposition_cells: 8,
        seed: 1927,
        source_id: 11,
    }
}

fn trace<const N: usize>(
    cx: &Cx<'_>,
    m: &ScatteringGrid<'_>,
    power: &[f64; N],
    options: ScatteringDepositionConfig,
) -> ScatteringMomentumDeposition<N> {
    m.deposit_multiple_scattering_momentum(
        cx,
        Point3::new(-1.0, 0.5, 0.5),
        Vec3::new(1.0, 0.0, 0.0),
        &[550.0; N],
        power,
        options,
    )
    .unwrap()
}

fn balance<const N: usize>(report: &ScatteringMomentumDeposition<N>) {
    for lane in 0..N {
        let scale = report.power.incident_power_w[lane] / LIGHT_M_S;
        for axis in 0..3 {
            assert!(
                report.momentum.balance_residual_n[lane][axis].abs() <= 2e-12 * scale,
                "momentum residual {:?}",
                report.momentum.balance_residual_n
            );
            let cells: f64 = report
                .momentum
                .cells
                .iter()
                .map(|cell| cell.force.mean_n[lane][axis])
                .sum();
            assert!((cells - report.momentum.force.mean_n[lane][axis]).abs() <= 2e-12 * scale);
        }
    }
}

#[test]
fn g1_absorption_transfers_exact_beer_momentum_without_sample_noise() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0], 0.4, 0.0, PhaseFunction::Isotropic);
        let options = config(2, 1);
        let report = trace(cx, &m, &[3.0, 7.0], options);
        let old = m
            .deposit_multiple_scattering_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0; 2],
                &[3.0, 7.0],
                options,
            )
            .unwrap();
        assert_eq!(report.power, old);
        for (lane, power) in [3.0, 7.0].into_iter().enumerate() {
            let expected = power * -(-1.2_f64).exp_m1() / LIGHT_M_S;
            assert!((report.momentum.force.mean_n[lane][0] / expected - 1.0).abs() < 2e-13);
            assert_eq!(report.momentum.force.mean_n[lane][1..], [0.0; 2]);
            assert_eq!(report.momentum.force.standard_error_n[lane], [0.0; 3]);
            for cell in &report.momentum.cells {
                let heat = report
                    .power
                    .cells
                    .iter()
                    .find(|heat| heat.cell_index == cell.cell_index)
                    .unwrap();
                assert_eq!(
                    cell.force.mean_n[lane][0],
                    heat.absorbed.mean_w[lane] / LIGHT_M_S
                );
            }
        }
        balance(&report);
    });
}

#[test]
fn g1_elastic_isotropic_scattering_exerts_force_and_retains_unresolved_momentum() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1e6], 0.0, 1.0, PhaseFunction::Isotropic);
        let options = config(4096, 1);
        let report = trace(cx, &m, &[LIGHT_M_S], options);
        assert_eq!(report.power.absorbed.mean_w, [0.0]);
        assert!(report.power.cells.is_empty());
        assert_eq!(report.power.unresolved.mean_w, [LIGHT_M_S]);
        assert_eq!(report.momentum.cells.len(), 1);
        for axis in 0..3 {
            let expected = if axis == 0 { 1.0 } else { 0.0 };
            let observed = report.momentum.force.mean_n[0][axis];
            assert!(
                (observed - expected).abs() < 6.0 * report.momentum.force.standard_error_n[0][axis]
            );
            assert!(
                report.momentum.unresolved.mean_n[0][axis].abs()
                    < 6.0 * report.momentum.unresolved.standard_error_n[0][axis]
            );
        }
        balance(&report);
        assert_eq!(report, trace(cx, &m, &[LIGHT_M_S], options));
    });
}

#[test]
fn g1_forward_and_backward_phase_transfer_their_actual_photon_recoil() {
    with_cx(|cx, _| {
        for g in [-0.9, 0.9] {
            let m = medium(
                cx,
                &[1e6],
                0.25,
                0.75,
                PhaseFunction::henyey_greenstein(g).unwrap(),
            );
            let options = config(4096, 1);
            let report = trace(cx, &m, &[LIGHT_M_S], options);
            let expected = 0.25 + 0.75 * (1.0 - g);
            assert!(
                (report.momentum.force.mean_n[0][0] - expected).abs()
                    < 6.0 * report.momentum.force.standard_error_n[0][0]
            );
            assert_eq!(report.power.absorbed.mean_w, [0.25 * LIGHT_M_S]);
            let old = m
                .deposit_multiple_scattering_power(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    Vec3::new(1.0, 0.0, 0.0),
                    &[550.0],
                    &[LIGHT_M_S],
                    options,
                )
                .unwrap();
            assert_eq!(
                report.power, old,
                "last-collision phase draws must not change energy histories"
            );
            balance(&report);
        }
    });
}

#[test]
fn g1_each_multiple_collision_history_telescopes_before_statistical_averaging() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1e6], 0.2, 0.8, PhaseFunction::Isotropic);
        let wavelengths = [550.0];
        let (absorption, scattering, extinction) = m.coefficients(&wavelengths).unwrap();
        let prepared = UniformOptics {
            medium: &m,
            wavelengths,
            absorption,
            scattering,
            extinction,
        };
        let mut tallies = MomentumTallies::<1>::new(1);
        let mut history = HistoryMomentum::new(&mut tallies, 0, 0);
        let mut segments = 0;
        let mut draws = [0.5, 0.2, 0.3, 0.5, 0.8, 0.7].into_iter();
        let result = trace_history_observed(
            &prepared,
            &mut || Ok(()),
            (
                Point3::new(0.5, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                LIGHT_M_S,
            ),
            (2, 2),
            &mut segments,
            &mut || draws.next().unwrap(),
            &mut |_, _| Ok(()),
            Some(&mut history),
        )
        .unwrap();
        assert_eq!(segments, 2);
        assert!(result.values[3] > 0.0);
        for axis in 0..3 {
            assert!(
                (history.force[axis].value() + history.escape[2][axis]
                    - if axis == 0 { 1.0 } else { 0.0 })
                .abs()
                    < 2e-14
            );
        }
        history.finish().unwrap();
    });
}

#[test]
fn g3_local_spectral_phase_and_transparent_cells_reach_the_same_force_owner() {
    with_cx(|cx, _| {
        let base = medium(cx, &[1e6; 2], 0.0, 0.0, PhaseFunction::Isotropic);
        let phases = [
            PhaseFunction::Isotropic,
            PhaseFunction::spectral_table(&[
                (380.0, PhaseFunction::henyey_greenstein(-0.8).unwrap()),
                (780.0, PhaseFunction::henyey_greenstein(0.8).unwrap()),
            ])
            .unwrap(),
        ];
        let field = CellScatteringGrid::try_new(
            cx,
            base.absorption_grid(),
            [380.0, 780.0],
            &[[0.0; 2]; 2],
            &[[0.0; 2], [1.0; 2]],
            CellPhaseFunction::Indexed {
                palette: &phases,
                cell_indices: &[0, 1],
            },
            (380.0, 780.0),
            ContentHash([41; 32]),
        )
        .unwrap();
        let options = config(4096, 1);
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let report = field
            .deposit_multiple_scattering_momentum(cx, origin, direction, &[LIGHT_M_S; 2], options)
            .unwrap();
        assert_eq!(report.momentum.cells.len(), 1);
        assert_eq!(report.momentum.cells[0].cell_index, 1);
        assert_eq!(
            report.power,
            field
                .deposit_multiple_scattering_power(cx, origin, direction, &[LIGHT_M_S; 2], options)
                .unwrap()
        );
        for (lane, expected) in [1.8, 0.2].into_iter().enumerate() {
            assert!(
                (report.momentum.force.mean_n[lane][0] - expected).abs()
                    < 6.0 * report.momentum.force.standard_error_n[lane][0]
            );
        }
        balance(&report);
    });
}

#[test]
fn g3_diffuse_escape_and_repeated_cell_visits_close_the_joint_account() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0, 0.5], 0.15, 0.85, PhaseFunction::Rayleigh);
        let options = config(512, 32);
        let report = trace(cx, &m, &[LIGHT_M_S], options);
        assert!(report.power.diffuse_escape.mean_w[0] > 0.0);
        assert!(report.power.absorbed.mean_w[0] > 0.0);
        let old = m
            .deposit_multiple_scattering_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[LIGHT_M_S],
                options,
            )
            .unwrap();
        assert_eq!(report.power, old);
        balance(&report);
    });
}

#[test]
fn g4_force_budget_and_cancellation_publish_no_partial_account() {
    with_cx(|cx, gate| {
        let m = medium(cx, &[1e6], 0.0, 1.0, PhaseFunction::Isotropic);
        let identity = m.content_identity();
        let options = ScatteringDepositionConfig {
            max_deposition_cells: 0,
            ..config(2, 1)
        };
        let query = || {
            m.deposit_multiple_scattering_momentum(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[1.0],
                options,
            )
        };
        assert_eq!(query(), Err(VolumeTransportError::Budget));
        assert_eq!(m.content_identity(), identity);
        gate.request();
        assert!(query().is_err());
    });
}

fn angular_balance<const N: usize>(report: &ScatteringWrenchDeposition<N>) {
    for lane in 0..N {
        let scale = report.power.incident_power_w[lane] / LIGHT_M_S;
        for axis in 0..3 {
            assert!(
                report.orbital.balance_residual_n_m[lane][axis].abs() <= 3e-12 * scale,
                "orbital residual {:?}",
                report.orbital.balance_residual_n_m
            );
            let cells: f64 = report
                .orbital
                .cells
                .iter()
                .map(|cell| cell.torque.mean_n_m[lane][axis])
                .sum();
            assert!((cells - report.orbital.torque.mean_n_m[lane][axis]).abs() <= 3e-12 * scale);
        }
    }
}

#[test]
fn g1_off_axis_absorption_preserves_beam_line_of_action_in_every_cell() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0], 0.4, 0.0, PhaseFunction::Isotropic);
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let options = config(2, 1);
        let report = m
            .deposit_multiple_scattering_wrench(
                cx,
                origin,
                direction,
                &[550.0; 2],
                &[3.0, 7.0],
                Point3::new(0.0, 0.0, 0.0),
                options,
            )
            .unwrap();
        let linear = m
            .deposit_multiple_scattering_momentum(
                cx,
                origin,
                direction,
                &[550.0; 2],
                &[3.0, 7.0],
                options,
            )
            .unwrap();
        assert_eq!(report.power, linear.power);
        assert_eq!(report.momentum, linear.momentum);
        assert_eq!(report.orbital.cells.len(), 2);
        for lane in 0..2 {
            let force = report.momentum.force.mean_n[lane][0];
            let expected = [0.0, 0.8 * force, -0.2 * force];
            for axis in 0..3 {
                assert!(
                    (report.orbital.torque.mean_n_m[lane][axis] - expected[axis]).abs()
                        < 2e-14 * force
                );
            }
            assert_eq!(report.orbital.torque.standard_error_n_m[lane], [0.0; 3]);
            for cell in &report.orbital.cells {
                let local_force = report
                    .momentum
                    .cells
                    .iter()
                    .find(|force| force.cell_index == cell.cell_index)
                    .unwrap()
                    .force
                    .mean_n[lane][0];
                assert_eq!(
                    cell.torque.mean_n_m[lane],
                    [0.0, 0.8 * local_force, -0.2 * local_force]
                );
                assert_eq!(cell.torque.standard_error_n_m[lane], [0.0; 3]);
            }
        }
        angular_balance(&report);
        let on_line = m
            .deposit_multiple_scattering_wrench(
                cx,
                origin,
                direction,
                &[550.0],
                &[3.0],
                Point3::new(-2.0, 0.2, 0.8),
                options,
            )
            .unwrap();
        assert_eq!(on_line.orbital.torque.mean_n_m, [[0.0; 3]]);
        assert!(on_line.orbital.cells.is_empty());
    });
}

#[test]
fn g1_native_collision_torque_uses_event_position_and_each_history_telescopes() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1e6], 0.2, 0.8, PhaseFunction::Isotropic);
        let wavelengths = [550.0];
        let (absorption, scattering, extinction) = m.coefficients(&wavelengths).unwrap();
        let prepared = UniformOptics {
            medium: &m,
            wavelengths,
            absorption,
            scattering,
            extinction,
        };
        let origin = Point3::new(0.2, 0.3, 0.4);
        let reference = Point3::new(-0.4, 0.1, 0.2);
        for collisions in [1, 2] {
            let mut tallies = MomentumTallies::<1>::with_orbital(1, reference).unwrap();
            let mut history = HistoryMomentum::new(&mut tallies, 0, 0);
            let mut segments = 0;
            let mut draws = [0.5, 0.2, 0.3, 0.4, 0.8, 0.7].into_iter();
            let result = trace_history_observed(
                &prepared,
                &mut || Ok(()),
                (origin, Vec3::new(1.0, 0.0, 0.0), LIGHT_M_S),
                (collisions, u64::from(collisions)),
                &mut segments,
                &mut || draws.next().unwrap(),
                &mut |_, _| Ok(()),
                Some(&mut history),
            )
            .unwrap();
            assert_eq!(segments, u64::from(collisions));
            assert!(result.values[3] > 0.0);
            let orbital = history.orbital.as_ref().unwrap();
            let incident = [0.0, origin.z - reference.z, -(origin.y - reference.y)];
            for axis in 0..3 {
                assert!(
                    (orbital.torque[axis].value() + orbital.escape[2][axis] - incident[axis]).abs()
                        < 3e-14
                );
            }
            if collisions == 1 {
                // Homogeneous exponential free flight: tau = -ln(1-u).
                let event_x = origin.x + std::f64::consts::LN_2 / 1e6;
                let force = history.force.map(PowerSum::value);
                let r = [
                    event_x - reference.x,
                    origin.y - reference.y,
                    origin.z - reference.z,
                ];
                let expected = [
                    r[1] * force[2] - r[2] * force[1],
                    r[2] * force[0] - r[0] * force[2],
                    r[0] * force[1] - r[1] * force[0],
                ];
                for axis in 0..3 {
                    assert!((orbital.torque[axis].value() - expected[axis]).abs() < 3e-14);
                }
                let cell_center =
                    orbital::moment_about(reference, Point3::new(0.5, 0.5, 0.5), force).unwrap();
                assert!((orbital.torque[1].value() - cell_center[1]).abs() > 1e-2);
            }
            history.finish().unwrap();
        }
    });
}

#[test]
fn g3_native_escape_and_unresolved_orbital_flux_transform_with_reference_point() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0, 0.5], 0.15, 0.85, PhaseFunction::Rayleigh);
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let options = config(512, 3);
        let query = |reference| {
            m.deposit_multiple_scattering_wrench(
                cx,
                origin,
                direction,
                &[550.0],
                &[LIGHT_M_S],
                reference,
                options,
            )
            .unwrap()
        };
        let a = query(Point3::new(0.0, 0.0, 0.0));
        let shift = Point3::new(0.3, -0.4, 0.2);
        let b = query(shift);
        assert_eq!(a.power, b.power);
        assert_eq!(a.momentum, b.momentum);
        assert!(a.power.diffuse_escape.mean_w[0] > 0.0);
        assert!(a.power.unresolved.mean_w[0] > 0.0);
        assert_eq!(a, query(Point3::new(0.0, 0.0, 0.0)));
        for (before, after, rate) in [
            (
                a.orbital.torque.mean_n_m[0],
                b.orbital.torque.mean_n_m[0],
                a.momentum.force.mean_n[0],
            ),
            (
                a.orbital.ballistic_escape.mean_n_m[0],
                b.orbital.ballistic_escape.mean_n_m[0],
                a.momentum.ballistic_escape.mean_n[0],
            ),
            (
                a.orbital.diffuse_escape.mean_n_m[0],
                b.orbital.diffuse_escape.mean_n_m[0],
                a.momentum.diffuse_escape.mean_n[0],
            ),
            (
                a.orbital.unresolved.mean_n_m[0],
                b.orbital.unresolved.mean_n_m[0],
                a.momentum.unresolved.mean_n[0],
            ),
            (
                a.orbital.incident_angular_momentum_rate_n_m[0],
                b.orbital.incident_angular_momentum_rate_n_m[0],
                a.momentum.incident_momentum_rate_n[0],
            ),
        ] {
            let reference_moment =
                orbital::moment_about(Point3::new(0.0, 0.0, 0.0), shift, rate).unwrap();
            for axis in 0..3 {
                assert!((after[axis] - (before[axis] - reference_moment[axis])).abs() < 3e-12);
            }
        }
        angular_balance(&a);
        angular_balance(&b);
    });
}

#[test]
fn g3_local_spectral_phase_reaches_original_cell_torque_owner() {
    with_cx(|cx, _| {
        let base = medium(cx, &[1e6; 2], 0.0, 0.0, PhaseFunction::Isotropic);
        let phases = [
            PhaseFunction::Isotropic,
            PhaseFunction::spectral_table(&[
                (380.0, PhaseFunction::henyey_greenstein(-0.8).unwrap()),
                (780.0, PhaseFunction::henyey_greenstein(0.8).unwrap()),
            ])
            .unwrap(),
        ];
        let field = CellScatteringGrid::try_new(
            cx,
            base.absorption_grid(),
            [380.0, 780.0],
            &[[0.0; 2], [0.3; 2]],
            &[[0.0; 2], [0.7; 2]],
            CellPhaseFunction::Indexed {
                palette: &phases,
                cell_indices: &[0, 1],
            },
            (380.0, 780.0),
            ContentHash([41; 32]),
        )
        .unwrap();
        let origin = Point3::new(-1.0, 0.2, 0.8);
        let direction = Vec3::new(1.0, 0.0, 0.0);
        let options = config(512, 1);
        let report = field
            .deposit_multiple_scattering_wrench(
                cx,
                origin,
                direction,
                &[LIGHT_M_S; 2],
                Point3::new(0.0, 0.0, 0.0),
                options,
            )
            .unwrap();
        assert_eq!(report.orbital.cells.len(), 1);
        assert_eq!(report.orbital.cells[0].cell_index, 1);
        let linear = field
            .deposit_multiple_scattering_momentum(cx, origin, direction, &[LIGHT_M_S; 2], options)
            .unwrap();
        assert_eq!(report.power, linear.power);
        assert_eq!(report.momentum, linear.momentum);
        for lane in 0..2 {
            let force = report.momentum.force.mean_n[lane];
            assert!(
                (report.orbital.torque.mean_n_m[lane][0] - (0.2 * force[2] - 0.8 * force[1])).abs()
                    < 3e-13
            );
        }
        angular_balance(&report);
    });
}

#[test]
fn g1_orbital_cross_product_retains_small_cancellation_residual() {
    let offset = 2.0_f64.powi(-27);
    let moment = orbital::moment_about(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 1.0 + offset, 1.0),
        [0.0, 1.0, 1.0 - offset],
    )
    .unwrap();
    assert_eq!(moment, [-2.0_f64.powi(-54), 0.0, 0.0]);
}

#[test]
fn g4_wrench_refuses_invalid_reference_overflow_budget_and_cancellation_atomically() {
    with_cx(|cx, gate| {
        let m = medium(cx, &[1.0], 0.4, 0.0, PhaseFunction::Isotropic);
        let identity = m.content_identity();
        let options = config(2, 1);
        let query = |reference, power, config| {
            m.deposit_multiple_scattering_wrench(
                cx,
                Point3::new(-1.0, 0.2, 0.8),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[power],
                reference,
                config,
            )
        };
        assert_eq!(
            query(Point3::new(f64::NAN, 0.0, 0.0), 1.0, options),
            Err(VolumeTransportError::Ray)
        );
        assert_eq!(
            query(Point3::new(0.0, f64::MAX, 0.0), 1e308, options),
            Err(VolumeTransportError::Unrepresentable)
        );
        assert_eq!(
            query(
                Point3::new(0.0, 0.0, 0.0),
                1.0,
                ScatteringDepositionConfig {
                    max_deposition_cells: 0,
                    ..options
                }
            ),
            Err(VolumeTransportError::Budget)
        );
        assert_eq!(m.content_identity(), identity);
        gate.request();
        assert!(query(Point3::new(0.0, 0.0, 0.0), 1.0, options).is_err());
        assert_eq!(m.content_identity(), identity);
    });
}
