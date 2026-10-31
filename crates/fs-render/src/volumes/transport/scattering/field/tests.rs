//! Layered-media references independent of the prepared coefficient implementation.

use super::*;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringGrid};
use crate::volumes::transport::{AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
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
    })
}

fn geometry<'a>(
    cx: &Cx<'_>,
    density: &'a [f64],
    temperature: Option<&'a [f64]>,
) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new([density.len(), 1, 1], density, [0.0; 3], [1.0, 4.0, 4.0]),
        MassAbsorptionSpectrum::constant(0.0).unwrap(),
        temperature,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 64,
        },
    )
    .unwrap()
}

fn config(histories: u32, collisions: u32) -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: histories,
        max_collisions_per_history: collisions,
        max_total_segments: u64::from(histories) * u64::from(collisions) * 4,
        max_deposition_cells: 16,
        seed: 1927,
        source_id: 11,
    }
}
fn origin() -> Point3 {
    Point3::new(-1.0, 2.0, 2.0)
}
fn direction() -> Vec3 {
    Vec3::new(1.0, 0.0, 0.0)
}
fn balance<const N: usize>(r: &ScatteringPowerDeposition<N>) {
    for lane in 0..N {
        let sum = r.ballistic_escape.mean_w[lane]
            + r.diffuse_escape.mean_w[lane]
            + r.absorbed.mean_w[lane]
            + r.unresolved.mean_w[lane];
        assert!((sum - r.incident_power_w[lane]).abs() < 3e-12 * r.incident_power_w[lane]);
        let cells: f64 = r.cells.iter().map(|c| c.absorbed.mean_w[lane]).sum();
        assert!((cells - r.absorbed.mean_w[lane]).abs() < 3e-12 * r.incident_power_w[lane]);
    }
}

#[test]
fn g1_uniform_cells_preserve_free_flight_beam_and_complete_history_results() {
    with_cx(|cx, _| {
        let density = [0.5, 1.5];
        let base = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0, 4.0, 4.0]),
            MassAbsorptionSpectrum::constant(0.3).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 2,
                max_steps_per_segment: 64,
            },
        )
        .unwrap();
        let uniform = ScatteringGrid::try_new(
            cx,
            base,
            MassScatteringSpectrum::constant(1.7).unwrap(),
            PhaseFunction::HenyeyGreenstein { g: -0.4 },
        )
        .unwrap();
        let a = [[0.3; 2]; 2];
        let s = [[1.7; 2]; 2];
        let field = CellScatteringGrid::try_new(
            cx,
            uniform.absorption_grid(),
            [450.0, 650.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(uniform.phase_ref().unwrap()),
            (380.0, 780.0),
            ContentHash([31; 32]),
        )
        .unwrap();
        for u in [0.001, 0.2, 0.7, 0.999] {
            assert_eq!(
                field
                    .sample_segment(cx, origin(), direction(), f64::INFINITY, 0.6, u)
                    .unwrap(),
                uniform
                    .sample_segment(
                        cx,
                        origin(),
                        direction(),
                        f64::INFINITY,
                        &[450.0, 650.0],
                        0.6,
                        u
                    )
                    .unwrap()
            );
        }
        let mut beam = field
            .deposit_beam_power(cx, origin(), direction(), f64::INFINITY, &[3.0, 7.0], 2)
            .unwrap();
        beam.medium_identity = uniform.content_identity();
        assert_eq!(
            beam,
            uniform
                .deposit_beam_power(
                    cx,
                    origin(),
                    direction(),
                    f64::INFINITY,
                    &[450.0, 650.0],
                    &[3.0, 7.0],
                    2
                )
                .unwrap()
        );
        let mut sampled = field
            .deposit_multiple_scattering_power(
                cx,
                origin(),
                direction(),
                &[3.0, 7.0],
                config(4000, 64),
            )
            .unwrap();
        sampled.medium_identity = uniform.content_identity();
        assert_eq!(
            sampled,
            uniform
                .deposit_multiple_scattering_power(
                    cx,
                    origin(),
                    direction(),
                    &[450.0, 650.0],
                    &[3.0, 7.0],
                    config(4000, 64)
                )
                .unwrap()
        );
        balance(&sampled);
    });
}

#[test]
fn g1_layered_ir_first_collision_uses_local_spectral_albedo_and_mixture_density() {
    with_cx(|cx, _| {
        let density = [0.5, 1.5];
        let base = geometry(cx, &density, None);
        let a = [[0.1, 0.6], [1.3, 0.2]];
        let s = [[0.9, 0.4], [0.2, 1.8]];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [1000.0, 2000.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (900.0, 2200.0),
            ContentHash([4; 32]),
        )
        .unwrap();
        // Hero lane 0 collides a quarter metre into the second cell.
        let sampled_tau: f64 = 0.5 * 1.0 + 1.5 * 1.5 * 0.25;
        let flight = field
            .sample_segment(
                cx,
                origin(),
                direction(),
                f64::INFINITY,
                0.1,
                (-sampled_tau).exp(),
            )
            .unwrap();
        let event = flight.event.unwrap();
        assert_eq!(event.cell_index, 1);
        assert!((event.position.x - 1.25).abs() < 2e-14);
        let t = [
            (-sampled_tau).exp(),
            (-(0.5 * 1.0 + 1.5 * 2.0 * 0.25_f64)).exp(),
        ];
        let mixture = 0.5 * (t[0] * 1.5 + t[1] * 2.0);
        for lane in 0..2 {
            assert!((flight.weight[lane] - t[lane] * s[1][lane] / mixture).abs() < 2e-14);
        }
        let report = field
            .deposit_multiple_scattering_power(
                cx,
                origin(),
                direction(),
                &[3.0, 7.0],
                config(40_000, 1),
            )
            .unwrap();
        assert_eq!(report.wavelength_bounds_nm, (900.0, 2200.0));
        for lane in 0..2 {
            let mut transmission = 1.0;
            for cell in 0..2 {
                let total = a[cell][lane] + s[cell][lane];
                let opacity = -(-density[cell] * total).exp_m1();
                let expected =
                    report.incident_power_w[lane] * transmission * opacity * a[cell][lane] / total;
                let actual = &report.cells[cell].absorbed;
                assert!(
                    (actual.mean_w[lane] - expected).abs()
                        < 6.0 * actual.standard_error_w[lane] + 1e-12
                );
                transmission *= (-density[cell] * total).exp();
            }
            assert!(
                (report.ballistic_escape.mean_w[lane]
                    - report.incident_power_w[lane] * transmission)
                    .abs()
                    < 6.0 * report.ballistic_escape.standard_error_w[lane] + 1e-12
            );
            assert!(report.unresolved.mean_w[lane] > 0.0);
        }
        balance(&report);
    });
}

#[test]
fn g1_local_absorption_only_lane_is_exact_even_when_other_lane_scatters() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0];
        let base = geometry(cx, &density, None);
        let a = [[0.2, 0.1], [0.7, 1.3]];
        let s = [[0.8, 0.0], [0.3, 0.0]];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [1100.0, 1800.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (1000.0, 2000.0),
            ContentHash([5; 32]),
        )
        .unwrap();
        let r = field
            .deposit_multiple_scattering_power(
                cx,
                origin(),
                direction(),
                &[1.0, 4.0],
                config(1000, 32),
            )
            .unwrap();
        let front = 4.0 * -(-0.1_f64).exp_m1();
        let rear = 4.0 * (-0.1_f64).exp() * -(-2.6_f64).exp_m1();
        assert_eq!(r.sampled_histories[1], 0);
        assert_eq!(r.absorbed.standard_error_w[1], 0.0);
        assert!((r.cells[0].absorbed.mean_w[1] - front).abs() < 2e-14);
        assert!((r.cells[1].absorbed.mean_w[1] - rear).abs() < 2e-14);
        assert!((r.ballistic_escape.mean_w[1] - 4.0 * (-2.7_f64).exp()).abs() < 2e-14);
        balance(&r);
    });
}

#[test]
fn g1_cell_thermal_source_uses_local_absorption_and_full_downstream_extinction() {
    with_cx(|cx, _| {
        let density = [0.5, 1.5];
        let temperatures = [700.0, 1500.0];
        let base = geometry(cx, &density, Some(&temperatures));
        assert!(!base.has_thermal_emission()); // The replaced uniform law is zero.
        let a = [[0.2], [0.8]];
        let s = [[0.8], [0.2]];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [1600.0],
            &a,
            &s,
            CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
            (1000.0, 2000.0),
            ContentHash([6; 32]),
        )
        .unwrap();
        let r = field
            .transport(cx, origin(), direction(), f64::INFINITY)
            .unwrap();
        let expected =
            -(-0.5_f64).exp_m1() * 0.2 * planck_radiance_w_m2_sr_nm(700.0, 1600.0).unwrap()
                + (-0.5_f64).exp()
                    * -(-1.5_f64).exp_m1()
                    * 0.8
                    * planck_radiance_w_m2_sr_nm(1500.0, 1600.0).unwrap();
        assert!((r.emission_w_m2_sr_nm[0] - expected).abs() < 2e-13 * expected);
        let sampled = field
            .sample_segment(cx, origin(), direction(), f64::INFINITY, 0.0, 0.99)
            .unwrap();
        assert_eq!(sampled.event.unwrap().cell_index, 0);
        assert_eq!(sampled.emission_w_m2_sr_nm, r.emission_w_m2_sr_nm);
    });
}

#[test]
fn g0_admission_refuses_unbound_fields_and_callback_errors_are_preserved() {
    with_cx(|cx, _| {
        let density = [1.0, 1.0];
        let base = geometry(cx, &density, None);
        let a = [[0.2], [0.7]];
        let s = [[0.8], [0.3]];
        let new = |a, s, wavelength, bounds, identity| {
            CellScatteringGrid::try_new(
                cx,
                &base,
                [wavelength],
                a,
                s,
                CellPhaseFunction::Uniform(&PhaseFunction::Isotropic),
                bounds,
                identity,
            )
        };
        assert!(matches!(
            new(&a[..1], &s, 1500.0, (1000.0, 2000.0), ContentHash([1; 32])),
            Err(VolumeTransportError::Grid)
        ));
        assert!(matches!(
            new(&a, &s, 999.0, (1000.0, 2000.0), ContentHash([1; 32])),
            Err(VolumeTransportError::Wavelength)
        ));
        assert!(matches!(
            new(&a, &s, 1500.0, (1000.0, 2000.0), ContentHash([0; 32])),
            Err(VolumeTransportError::SourceIdentity)
        ));
        assert!(matches!(
            new(
                &[[f64::NAN], [0.2]],
                &s,
                1500.0,
                (1000.0, 2000.0),
                ContentHash([1; 32])
            ),
            Err(VolumeTransportError::AbsorptionSpectrum)
        ));
        assert!(matches!(
            new(
                &a,
                &[[-0.1], [0.2]],
                1500.0,
                (1000.0, 2000.0),
                ContentHash([1; 32])
            ),
            Err(VolumeTransportError::ScatteringSpectrum)
        ));
        let phases = [PhaseFunction::Isotropic];
        assert!(matches!(
            CellScatteringGrid::try_new(
                cx,
                &base,
                [1500.0],
                &a,
                &s,
                CellPhaseFunction::Indexed {
                    palette: &phases,
                    cell_indices: &[0, 1]
                },
                (1000.0, 2000.0),
                ContentHash([1; 32])
            ),
            Err(VolumeTransportError::Phase)
        ));
        let field = new(&a, &s, 1500.0, (1000.0, 2000.0), ContentHash([1; 32])).unwrap();
        let mut checkpoints = 0;
        let mut check = || {
            checkpoints += 1;
            if checkpoints == 9 {
                Err(VolumeTransportError::Budget)
            } else {
                Ok(())
            }
        };
        let dynamic: &mut dyn FnMut() -> Result<(), VolumeTransportError> = &mut check;
        assert!(matches!(
            field.deposit_multiple_scattering_power_with_checkpoint(
                dynamic,
                origin(),
                direction(),
                &[1.0],
                config(100, 16)
            ),
            Err(VolumeTransportError::Budget)
        ));
        assert!(matches!(
            field.deposit_multiple_scattering_power_with_checkpoint(
                &mut || Err(VolumeTransportError::Cancelled),
                origin(),
                direction(),
                &[1.0],
                config(100, 16)
            ),
            Err(VolumeTransportError::Cancelled)
        ));
        let mut bounded = config(100, 16);
        bounded.max_total_segments = 1;
        assert!(matches!(
            field.deposit_multiple_scattering_power(cx, origin(), direction(), &[1.0], bounded),
            Err(VolumeTransportError::Budget)
        ));
        let repeated = field
            .deposit_multiple_scattering_power(cx, origin(), direction(), &[1.0], config(100, 16))
            .unwrap();
        assert_eq!(
            repeated,
            field
                .deposit_multiple_scattering_power(
                    cx,
                    origin(),
                    direction(),
                    &[1.0],
                    config(100, 16)
                )
                .unwrap()
        );
    });
}

// Independent analog random walk in a two-layer rectangular box. It uses
// elementary plane intersections, one scalar exponential optical depth and
// binary absorption (not implicit capture), standard-library elementary
// functions, a separate LCG, and a directly inverted HG angular CDF.
fn analog_two_layer(histories: usize, a: [f64; 2], s: [f64; 2], g: [f64; 2]) -> [f64; 2] {
    let mut state = 7919_u64;
    let mut uniform = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64 + 0.5) / 9_007_199_254_740_992.0
    };
    let mut captured = [0_usize; 2];
    for _ in 0..histories {
        let mut p = [0.0, 2.0, 2.0];
        let mut d = [1.0, 0.0, 0.0];
        let mut cell = 0;
        'history: for _ in 0..512 {
            let mut tau = -uniform().ln();
            loop {
                let face = [
                    if d[0] >= 0.0 {
                        cell as f64 + 1.0
                    } else {
                        cell as f64
                    },
                    if d[1] >= 0.0 { 4.0 } else { 0.0 },
                    if d[2] >= 0.0 { 4.0 } else { 0.0 },
                ];
                let distance: [f64; 3] = std::array::from_fn(|axis| {
                    if d[axis] == 0.0 {
                        f64::INFINITY
                    } else {
                        (face[axis] - p[axis]) / d[axis]
                    }
                });
                let axis = (0..3)
                    .min_by(|&i, &j| distance[i].total_cmp(&distance[j]))
                    .unwrap();
                let extinction = a[cell] + s[cell];
                if tau < extinction * distance[axis] {
                    for k in 0..3 {
                        p[k] += d[k] * (tau / extinction);
                    }
                    if uniform() < a[cell] / extinction {
                        captured[cell] += 1;
                        break 'history;
                    }
                    let u = uniform();
                    let cosine = if g[cell].abs() < 1e-12 {
                        2.0 * u - 1.0
                    } else {
                        let ratio = (1.0 - g[cell] * g[cell]) / (1.0 - g[cell] + 2.0 * g[cell] * u);
                        ((1.0 + g[cell] * g[cell] - ratio * ratio) / (2.0 * g[cell]))
                            .clamp(-1.0, 1.0)
                    };
                    let sine = (1.0 - cosine * cosine).max(0.0).sqrt();
                    let phi = 2.0 * std::f64::consts::PI * uniform();
                    let mut tangent = if d[2].abs() < 0.9 {
                        [d[1], -d[0], 0.0]
                    } else {
                        [-d[2], 0.0, d[0]]
                    };
                    let length = tangent.iter().map(|v| v * v).sum::<f64>().sqrt();
                    for value in &mut tangent {
                        *value /= length;
                    }
                    let binormal = [
                        d[1] * tangent[2] - d[2] * tangent[1],
                        d[2] * tangent[0] - d[0] * tangent[2],
                        d[0] * tangent[1] - d[1] * tangent[0],
                    ];
                    d = std::array::from_fn(|k| {
                        cosine * d[k] + sine * (phi.cos() * tangent[k] + phi.sin() * binormal[k])
                    });
                    break;
                }
                tau -= extinction * distance[axis];
                for k in 0..3 {
                    p[k] += d[k] * distance[axis];
                }
                p[axis] = face[axis];
                if axis != 0 || (cell == 0 && d[0] < 0.0) || (cell == 1 && d[0] > 0.0) {
                    break 'history;
                }
                cell = 1 - cell;
            }
        }
    }
    captured.map(|count| count as f64 / histories as f64)
}

#[test]
fn g1_cell_phase_palette_and_repeat_visits_match_independent_analog_transport() {
    with_cx(|cx, _| {
        let density = [1.0, 1.0];
        let base = geometry(cx, &density, None);
        let a = [[0.2], [0.9]];
        let s = [[1.8], [1.1]];
        let phases = [
            PhaseFunction::henyey_greenstein(-0.6).unwrap(),
            PhaseFunction::henyey_greenstein(0.6).unwrap(),
        ];
        let field = CellScatteringGrid::try_new(
            cx,
            &base,
            [550.0],
            &a,
            &s,
            CellPhaseFunction::Indexed {
                palette: &phases,
                cell_indices: &[0, 1],
            },
            (380.0, 780.0),
            ContentHash([27; 32]),
        )
        .unwrap();
        let report = field
            .deposit_multiple_scattering_power(
                cx,
                origin(),
                direction(),
                &[1.0],
                config(30_000, 128),
            )
            .unwrap();
        let reference_count = 120_000;
        let reference = analog_two_layer(reference_count, [0.2, 0.9], [1.8, 1.1], [-0.6, 0.6]);
        for (cell, expected) in report.cells.iter().zip(reference) {
            let reference_error = (expected * (1.0 - expected) / reference_count as f64).sqrt();
            let combined_error = (reference_error * reference_error
                + cell.absorbed.standard_error_w[0].powi(2))
            .sqrt();
            assert!(
                (cell.absorbed.mean_w[0] - expected).abs() < 6.0 * combined_error + 1e-10,
                "cell {}: native {}, analog {}, combined SE {}",
                cell.cell_index,
                cell.absorbed.mean_w[0],
                expected,
                combined_error
            );
        }
        assert!(report.traced_segments > 2 * 30_000);
        assert!(report.unresolved.mean_w[0] < 1e-10);
        assert!(report.absorbed.standard_error_w[0] > 0.0);
        balance(&report);
    });
}
