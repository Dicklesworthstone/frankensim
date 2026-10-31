//! Analytic first-collision and independent analog-transport references.

use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
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

fn config(samples: u32, depth: u32) -> ScatteringDepositionConfig {
    ScatteringDepositionConfig {
        histories_per_lane: samples,
        max_collisions_per_history: depth,
        max_total_segments: u64::from(samples) * u64::from(depth) * 4,
        max_deposition_cells: 16,
        seed: 1927,
        source_id: 11,
    }
}

fn trace<const N: usize>(
    cx: &Cx<'_>,
    medium: &ScatteringGrid<'_>,
    wavelengths: &[f64; N],
    power: &[f64; N],
    options: ScatteringDepositionConfig,
) -> ScatteringPowerDeposition<N> {
    medium
        .deposit_multiple_scattering_power(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            wavelengths,
            power,
            options,
        )
        .unwrap()
}

fn balance<const N: usize>(report: &ScatteringPowerDeposition<N>) {
    for lane in 0..N {
        let total = report.ballistic_escape.mean_w[lane]
            + report.diffuse_escape.mean_w[lane]
            + report.absorbed.mean_w[lane]
            + report.unresolved.mean_w[lane];
        assert!(
            (total - report.incident_power_w[lane]).abs()
                < 2e-12 * report.incident_power_w[lane].max(1e-300)
        );
        let cell_sum: f64 = report.cells.iter().map(|c| c.absorbed.mean_w[lane]).sum();
        assert!(
            (cell_sum - report.absorbed.mean_w[lane]).abs()
                < 2e-12 * report.incident_power_w[lane].max(1e-300)
        );
    }
}

#[test]
fn g1_absorption_only_is_exact_spatial_beer_transport_without_sample_noise() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0], 0.4, 0.0, PhaseFunction::Isotropic);
        let report = trace(cx, &m, &[530.0, 680.0], &[3.0, 7.0], config(2, 1));
        let exact = m
            .deposit_beam_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
                &[530.0, 680.0],
                &[3.0, 7.0],
                2,
            )
            .unwrap();
        assert_eq!(report.absorbed.mean_w, exact.absorbed_power_w);
        assert_eq!(report.ballistic_escape.mean_w, exact.transmitted_power_w);
        assert_eq!(report.sampled_histories, [0; 2]);
        assert_eq!(report.absorbed.standard_error_w, [0.0; 2]);
        assert_eq!(report.traced_segments, 2);
        for (sampled, deterministic) in report.cells.iter().zip(exact.cells) {
            assert_eq!(sampled.absorbed.mean_w, deterministic.absorbed_power_w);
        }
        balance(&report);
    });
}

#[test]
fn g1_first_collision_matches_heterogeneous_spectral_absorption_and_retains_tail() {
    with_cx(|cx, _| {
        let density = [0.5, 1.5];
        let a = MassAbsorptionSpectrum::table(&[(380.0, 0.2), (780.0, 1.0)]).unwrap();
        let base = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0; 3]),
            a,
            None,
            AbsorbingGridBudget {
                max_grid_cells: 2,
                max_steps_per_segment: 16,
            },
        )
        .unwrap();
        let m = ScatteringGrid::try_new(
            cx,
            base,
            MassScatteringSpectrum::constant(1.0).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let power = [3.0, 7.0];
        let report = trace(cx, &m, &[380.0, 780.0], &power, config(40_000, 1));
        for lane in 0..2 {
            let absorption = [0.2_f64, 1.0][lane];
            let extinction = absorption + 1.0;
            let ballistic = power[lane] * (-2.0 * extinction).exp();
            assert!(
                (report.ballistic_escape.mean_w[lane] - ballistic).abs()
                    < 5.0 * report.ballistic_escape.standard_error_w[lane]
            );
            assert_eq!(report.diffuse_escape.mean_w[lane], 0.0);
            let expected_tail = (power[lane] - ballistic) / extinction;
            assert!(
                (report.unresolved.mean_w[lane] - expected_tail).abs()
                    < 5.0 * report.unresolved.standard_error_w[lane]
            );
            for cell in &report.cells {
                let before = if cell.cell_index == 0 { 0.0 } else { 0.5 };
                let expected = power[lane]
                    * (-extinction * before).exp()
                    * -(-extinction * density[cell.cell_index]).exp_m1()
                    * absorption
                    / extinction;
                assert!(
                    (cell.absorbed.mean_w[lane] - expected).abs()
                        < 5.0 * cell.absorbed.standard_error_w[lane]
                );
            }
        }
        balance(&report);
    });
}

#[test]
fn source_backed_infrared_transport_uses_the_same_scattering_and_heat_owner() {
    with_cx(|cx, _| {
        let mut provenance = DomainHasher::new("absorption-deposition-infrared-fixture");
        provenance.update(b"declared test coefficients, not a specimen measurement");
        let identity = provenance.finalize();
        let a =
            MassAbsorptionSpectrum::source_table(&[(200.0, 0.4), (2500.0, 0.4)], identity).unwrap();
        let s =
            MassScatteringSpectrum::source_table(&[(200.0, 1.6), (2500.0, 1.6)], identity).unwrap();
        let base = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &[1.0], [0.0; 3], [1.0; 3]),
            a,
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 8,
            },
        )
        .unwrap();
        let medium = ScatteringGrid::try_new(cx, base, s, PhaseFunction::Rayleigh).unwrap();
        let report = trace(cx, &medium, &[1064.0], &[10.0], config(10_000, 1));
        let expected = 2.0 * -(-2.0_f64).exp_m1();
        assert!(
            (report.absorbed.mean_w[0] - expected).abs()
                < 5.0 * report.absorbed.standard_error_w[0]
        );
        assert_eq!(report.wavelength_bounds_nm, (200.0, 2500.0));
        balance(&report);
    });
}

/// Independent analog simulation of a homogeneous unit cube. This uses direct
/// ray/box distances, standard-library transcendentals, a separate RNG stream,
/// an absorption Bernoulli draw, and a separately derived HG frame/inverse CDF.
/// It shares no production DDA, implicit-capture, phase or tally implementation.
fn analog_absorption(a: f64, s: f64, g: f64, histories: u32) -> (f64, f64) {
    let mut rng = StreamKey {
        seed: 82719,
        kernel: 71,
        tile: 91,
    }
    .stream();
    let mut absorbed = 0_u32;
    for _ in 0..histories {
        let mut p = [0.0, 0.5, 0.5];
        let mut d = [1.0, 0.0, 0.0];
        loop {
            let exit = (0..3)
                .filter(|&axis| d[axis] != 0.0)
                .map(|axis| {
                    (if d[axis] > 0.0 {
                        1.0 - p[axis]
                    } else {
                        -p[axis]
                    }) / d[axis]
                })
                .fold(f64::INFINITY, f64::min);
            let distance = -(1.0 - rng.next_f64()).ln() / (a + s);
            if distance >= exit {
                break;
            }
            for axis in 0..3 {
                p[axis] += distance * d[axis];
            }
            if rng.next_f64() < a / (a + s) {
                absorbed += 1;
                break;
            }
            let u = rng.next_f64();
            let cosine = if g == 0.0 {
                2.0 * u - 1.0
            } else {
                let ratio = (1.0 - g * g) / (1.0 - g + 2.0 * g * u);
                (1.0 + g * g - ratio * ratio) / (2.0 * g)
            };
            let sine = (1.0 - cosine * cosine).max(0.0).sqrt();
            let angle = core::f64::consts::TAU * rng.next_f64();
            let helper = if d[0].abs() < 0.8 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            let mut tangent = cross(d, helper);
            let norm = tangent.iter().map(|v| v * v).sum::<f64>().sqrt();
            for v in &mut tangent {
                *v /= norm;
            }
            let bitangent = cross(d, tangent);
            d = std::array::from_fn(|axis| {
                cosine * d[axis]
                    + sine * (angle.cos() * tangent[axis] + angle.sin() * bitangent[axis])
            });
        }
    }
    let probability = f64::from(absorbed) / f64::from(histories);
    (
        probability,
        (probability * (1.0 - probability) / f64::from(histories)).sqrt(),
    )
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn g2_multiple_scattering_matches_independent_analog_cube_for_anisotropic_phases() {
    with_cx(|cx, _| {
        for g in [-0.7, 0.0, 0.7] {
            let m = medium(
                cx,
                &[1.0],
                0.4,
                3.0,
                PhaseFunction::henyey_greenstein(g).unwrap(),
            );
            let report = trace(cx, &m, &[550.0], &[1.0], config(30_000, 128));
            let (reference, reference_error) = analog_absorption(0.4, 3.0, g, 40_000);
            let combined =
                (reference_error.powi(2) + report.absorbed.standard_error_w[0].powi(2)).sqrt();
            assert!(
                (report.absorbed.mean_w[0] - reference).abs() < 6.0 * combined,
                "g={g}: {:?} vs {reference}, combined error {combined}",
                report.absorbed
            );
            assert_eq!(report.unresolved.mean_w[0], 0.0);
            assert!((report.cells[0].absorbed.mean_w[0] - report.absorbed.mean_w[0]).abs() < 2e-14);
            assert!(
                (report.cells[0].absorbed.standard_error_w[0]
                    - report.absorbed.standard_error_w[0])
                    .abs()
                    < 2e-14,
                "one-cell repeated collisions must preserve total variance"
            );
            let first_capture = (0.4 / 3.4) * -(-3.4_f64).exp_m1();
            assert!(report.absorbed.mean_w[0] > first_capture + 0.02);
            balance(&report);
        }
    });
}

#[test]
fn g3_conservative_scattering_and_high_albedo_caps_never_become_heat_or_disappear() {
    with_cx(|cx, _| {
        for a in [0.0, 0.001] {
            let m = medium(cx, &[1.0], a, 20.0, PhaseFunction::Isotropic);
            let report = trace(cx, &m, &[550.0], &[10.0], config(2_000, 5));
            assert!(report.unresolved.mean_w[0] > 1.0);
            if a == 0.0 {
                assert!(report.cells.is_empty());
                assert_eq!(report.absorbed.mean_w, [0.0]);
            } else {
                assert!(report.absorbed.mean_w[0] > 0.0);
            }
            balance(&report);
        }
    });
}

#[test]
fn g0_sparse_history_moments_include_missing_histories_and_repeat_visit_covariance() {
    let mut cell = CellStatistics::default();
    cell.capture(0, 1.0e280).unwrap();
    cell.capture(0, 2.0e280).unwrap();
    cell.capture(2, 6.0e280).unwrap();
    cell.flush().unwrap();
    cell.moments.zeroes_until(4);
    let (mean, se) = cell.moments.estimate().unwrap();
    assert!((mean / 1e280 - 2.25).abs() < 1e-14);
    assert!((se / 1e280 - 2.0625_f64.sqrt()).abs() < 1e-14);
}

#[test]
fn g0_near_unit_albedo_retains_cumulative_sub_ulp_absorption() {
    let mut surviving = PowerSum {
        sum: 1.0,
        correction: 0.0,
    };
    let mut heat = PowerSum::default();
    let a = 1.0e-16;
    for _ in 0..100_000 {
        let power = surviving.value();
        let captured = split_power(power, a, 1.0).unwrap().unwrap();
        let continuation = split_power(power, 1.0, 1.0).unwrap().unwrap();
        heat.add(captured).unwrap();
        surviving
            .continue_after_capture(captured, continuation)
            .unwrap();
    }
    let expected = -(100_000.0 * (-a).ln_1p()).exp_m1();
    assert!((heat.value() - expected).abs() < 2e-25);
    assert!((1.0 - surviving.value() - heat.value()).abs() < f64::EPSILON);
    assert!(surviving.value() < 1.0 - 9e-12);
}

#[test]
fn g3_underflowing_continuation_retains_the_packet_as_unresolved_radiation() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0], 1.0e6, 1.0, PhaseFunction::Isotropic);
        let report = m
            .deposit_multiple_scattering_power(
                cx,
                Point3::new(0.5, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[1.0e200],
                config(2, 512),
            )
            .unwrap();
        assert_eq!(report.precision_limited_histories, [2]);
        assert!(report.unresolved.mean_w[0] > 0.0);
        assert!(report.absorbed.mean_w[0] > 0.99e200);
        balance(&report);
    });
}

#[test]
fn g4_work_storage_and_cancellation_refuse_and_g5_source_streams_replay() {
    with_cx(|cx, gate| {
        let m = medium(cx, &[1.0], 1.0, 2.0, PhaseFunction::Isotropic);
        let options = config(100, 16);
        let first = trace(cx, &m, &[550.0], &[1.0], options);
        assert_eq!(first, trace(cx, &m, &[550.0], &[1.0], options));
        let different = trace(
            cx,
            &m,
            &[550.0],
            &[1.0],
            ScatteringDepositionConfig {
                source_id: 12,
                ..options
            },
        );
        assert_ne!(first.absorbed.mean_w, different.absorbed.mean_w);
        for refusal in [
            ScatteringDepositionConfig {
                max_total_segments: 1,
                ..options
            },
            ScatteringDepositionConfig {
                max_deposition_cells: 0,
                ..options
            },
        ] {
            assert_eq!(
                m.deposit_multiple_scattering_power(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    Vec3::new(1.0, 0.0, 0.0),
                    &[550.0],
                    &[1.0],
                    refusal
                ),
                Err(VolumeTransportError::Budget)
            );
        }
        gate.request();
        assert_eq!(
            m.deposit_multiple_scattering_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[1.0],
                options
            ),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
