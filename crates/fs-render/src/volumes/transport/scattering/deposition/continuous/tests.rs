//! Analytic absorption and an independent analog multiple-scattering reference.

use super::*;
use super::super::tests::{
    analog_absorption, balance, config, indexed_medium, medium, palette, with_cx,
};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use crate::volumes::transport::scattering::{
    CellPhaseFunction, MassScatteringSpectrum, ScatteringOptics,
};

fn trace<const N: usize>(
    cx: &Cx<'_>,
    medium: &ScatteringGrid<'_>,
    wavelengths: &[f64; N],
    power: &[f64; N],
    options: ScatteringDepositionConfig,
) -> ScatteringPowerDeposition<N> {
    medium
        .deposit_continuous_absorption_power(
            cx,
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            wavelengths,
            power,
            options,
        )
        .unwrap()
}

#[test]
fn g1_pure_absorption_is_the_original_exact_zero_variance_solve() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0], 0.4, 0.0, PhaseFunction::Isotropic);
        let mut report = trace(cx, &m, &[530.0, 680.0], &[3.0, 7.0], config(2, 1));
        assert_eq!(report.absorption_estimator, AbsorptionEstimator::Continuous);
        assert_eq!(report.sampled_histories, [0; 2]);
        assert_eq!(report.absorbed.standard_error_w, [0.0; 2]);
        balance(&report);
        report.absorption_estimator = AbsorptionEstimator::Collision;
        let original = m
            .deposit_multiple_scattering_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[530.0, 680.0],
                &[3.0, 7.0],
                config(2, 1),
            )
            .unwrap();
        assert_eq!(report, original);
    });
}

#[test]
fn g1_thin_absorbers_heat_even_when_no_collision_is_sampled() {
    with_cx(|cx, _| {
        let m = medium(cx, &[1.0, 2.0], 1e-7, 1e-12, PhaseFunction::Isotropic);
        let options = config(32, 8);
        let collision = m
            .deposit_multiple_scattering_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[10.0],
                options,
            )
            .unwrap();
        let report = trace(cx, &m, &[550.0], &[10.0], options);
        let expected = 10.0 * -(-3e-7_f64).exp_m1();
        assert_eq!(collision.absorbed.mean_w, [0.0]);
        assert!((report.absorbed.mean_w[0] - expected).abs() < 2e-13 * expected);
        assert_eq!(report.diffuse_escape.mean_w, [0.0]);
        // This fixed run has identical ballistic histories. Its empirical zero
        // error does not bound the effect of an unobserved rare scattering event.
        assert_eq!(report.absorbed.standard_error_w, [0.0]);
        assert_eq!(report.cells.len(), 2);
        assert!(report.cells.iter().all(|c| c.absorbed.mean_w[0] > 0.0));
        balance(&report);
    });
}

#[test]
fn g1_first_scatter_cap_matches_independent_layered_heat_and_tail_integrals() {
    with_cx(|cx, _| {
        let m = indexed_medium(
            cx,
            &[1.0; 4],
            palette(
                cx,
                [0.2, 0.8]
                    .map(|a| MassAbsorptionSpectrum::constant(a).unwrap())
                    .to_vec(),
                &[0, 0, 1, 1],
            ),
            palette(
                cx,
                [0.0, 0.3, 1.2]
                    .map(|s| ScatteringOptics {
                        scattering: MassScatteringSpectrum::constant(s).unwrap(),
                        phase: PhaseFunction::Isotropic,
                    })
                    .to_vec(),
                &[0, 2, 1, 0],
            ),
        );
        let report = trace(cx, &m, &[550.0], &[7.0], config(20_000, 1));
        let mut before = 7.0;
        let mut total_heat = 0.0;
        let mut tail = 0.0;
        for (index, (a, s)) in [(0.2_f64, 0.0), (0.2, 1.2), (0.8, 0.3), (0.8, 0.0)]
            .into_iter()
            .enumerate()
        {
            // Zero-scattering first/last cells absorb without drawing or
            // consuming scattering optical depth. E[segment absorption] is
            // integral P * mu_a * exp(-(mu_a + mu_s) * x) dx.
            let removed = before * -(-(a + s)).exp_m1();
            let expected = removed * a / (a + s);
            total_heat += expected;
            tail += removed * s / (a + s);
            before *= (-(a + s)).exp();
            let cell = report.cells.iter().find(|c| c.cell_index == index).unwrap();
            assert!(
                (cell.absorbed.mean_w[0] - expected).abs()
                    <= 6.0 * cell.absorbed.standard_error_w[0] + 2e-12
            );
        }
        for (estimate, expected) in [
            (report.ballistic_escape, before),
            (report.absorbed, total_heat),
            (report.unresolved, tail),
        ] {
            assert!(
                (estimate.mean_w[0] - expected).abs()
                    <= 6.0 * estimate.standard_error_w[0] + 2e-12
            );
        }
        assert_eq!(report.diffuse_escape.mean_w, [0.0]);
        balance(&report);
    });
}

#[test]
fn g2_multiple_scattering_agrees_with_independent_analog_photon_transport() {
    with_cx(|cx, _| {
        for g in [-0.7, 0.0, 0.7] {
            let m = medium(
                cx,
                &[1.0],
                0.4,
                3.0,
                PhaseFunction::henyey_greenstein(g).unwrap(),
            );
            let report = trace(cx, &m, &[550.0], &[1.0], config(15_000, 128));
            let (reference, reference_error) = analog_absorption(0.4, 3.0, g, 40_000);
            let combined =
                (reference_error.powi(2) + report.absorbed.standard_error_w[0].powi(2)).sqrt();
            assert!(
                (report.absorbed.mean_w[0] - reference).abs() < 6.0 * combined,
                "g={g}: {:?} vs {reference}, combined error {combined}",
                report.absorbed
            );
            assert_eq!(report.unresolved.mean_w, [0.0]);
            assert_eq!(report.cells.len(), 1);
            assert!(
                (report.cells[0].absorbed.standard_error_w[0]
                    - report.absorbed.standard_error_w[0])
                    .abs()
                    < 2e-14,
                "repeated visits must be summed within a history"
            );
            balance(&report);
        }
    });
}

fn subdivided<'a>(cx: &Cx<'_>, density: &'a [f64]) -> ScatteringGrid<'a> {
    let absorption = AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(
            [density.len(), 1, 1],
            density,
            [0.0; 3],
            [1.0 / density.len() as f64, 1.0, 1.0],
        ),
        MassAbsorptionSpectrum::constant(1e-4).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: density.len(),
            max_steps_per_segment: 128,
        },
    )
    .unwrap();
    ScatteringGrid::try_new(
        cx,
        absorption,
        MassScatteringSpectrum::constant(0.7).unwrap(),
        PhaseFunction::henyey_greenstein(0.6).unwrap(),
    )
    .unwrap()
}

#[test]
fn g3_native_subdivision_preserves_paths_and_integrated_heat_with_replay() {
    with_cx(|cx, _| {
        let coarse_density: &[f64] = &[1.0];
        let fine_density: &[f64] = &[1.0; 16];
        let coarse = subdivided(cx, coarse_density);
        let fine = subdivided(cx, fine_density);
        let options = config(1_024, 32);
        let a = trace(cx, &coarse, &[550.0], &[10.0], options);
        let b = trace(cx, &fine, &[550.0], &[10.0], options);
        assert_eq!(b, trace(cx, &fine, &[550.0], &[10.0], options));
        assert_eq!(a.traced_segments, b.traced_segments);
        for (left, right) in [
            (a.ballistic_escape, b.ballistic_escape),
            (a.diffuse_escape, b.diffuse_escape),
            (a.absorbed, b.absorbed),
            (a.unresolved, b.unresolved),
        ] {
            for (x, y) in [
                (left.mean_w[0], right.mean_w[0]),
                (left.standard_error_w[0], right.standard_error_w[0]),
            ] {
                assert!((x - y).abs() <= 3e-12 * x.abs().max(y.abs()).max(1e-300));
            }
        }
        balance(&a);
        balance(&b);
    });
}

#[test]
fn g3_prepared_cell_sources_use_original_coefficients_and_phase() {
    with_cx(|cx, _| {
        let phase = PhaseFunction::henyey_greenstein(-0.4).unwrap();
        let m = medium(cx, &[1.0, 2.0], 0.4, 2.0, phase.clone());
        let absorption = [[0.4]; 2];
        let scattering = [[2.0]; 2];
        let field = CellScatteringGrid::try_new(
            cx,
            m.absorption_grid(),
            [550.0],
            &absorption,
            &scattering,
            CellPhaseFunction::Uniform(&phase),
            (380.0, 780.0),
            m.content_identity(),
        )
        .unwrap();
        let options = config(512, 64);
        let actual = field
            .deposit_continuous_absorption_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[10.0],
                options,
            )
            .unwrap();
        let mut expected = trace(cx, &m, &[550.0], &[10.0], options);
        expected.medium_identity = actual.medium_identity;
        assert_eq!(actual, expected);
        let mut polls = 0;
        let error = field.deposit_continuous_absorption_power_with_checkpoint(
            &mut || {
                polls += 1;
                if polls == 5 {
                    Err(VolumeTransportError::Cancelled)
                } else {
                    Ok(())
                }
            },
            Point3::new(-1.0, 0.5, 0.5),
            Vec3::new(1.0, 0.0, 0.0),
            &[10.0],
            options,
        );
        assert_eq!(error, Err(VolumeTransportError::Cancelled));
    });
}

#[test]
fn g4_scattering_caps_and_budget_refusals_never_create_heat() {
    with_cx(|cx, gate| {
        for a in [0.0, 0.001] {
            let m = medium(cx, &[1.0], a, 20.0, PhaseFunction::Isotropic);
            let report = trace(cx, &m, &[550.0], &[10.0], config(256, 5));
            assert!(report.unresolved.mean_w[0] > 1.0);
            if a == 0.0 {
                assert_eq!(report.absorbed.mean_w, [0.0]);
                assert!(report.cells.is_empty());
            } else {
                assert!(report.absorbed.mean_w[0] > 0.0);
            }
            balance(&report);
        }
        let tiny = medium(cx, &[1.0], 0.01, 1e-12, PhaseFunction::Isotropic);
        let tiny_options = config(32, 8);
        let tiny_source = f64::from_bits(1);
        let tiny_report = trace(cx, &tiny, &[550.0], &[tiny_source], tiny_options);
        assert_eq!(tiny_report.absorbed.mean_w, [0.0]);
        assert!(tiny_report.cells.is_empty());
        assert_eq!(tiny_report.ballistic_escape.mean_w, [0.0]);
        assert_eq!(tiny_report.diffuse_escape.mean_w, [0.0]);
        assert_eq!(tiny_report.unresolved.mean_w, [tiny_source]);
        assert_eq!(
            tiny_report.precision_limited_histories,
            [tiny_options.histories_per_lane]
        );
        balance(&tiny_report);

        let m = medium(cx, &[1.0], 0.4, 2.0, PhaseFunction::Isotropic);
        let options = config(2, 8);
        for refused in [
            ScatteringDepositionConfig { max_total_segments: 1, ..options },
            ScatteringDepositionConfig { max_deposition_cells: 0, ..options },
        ] {
            assert_eq!(
                m.deposit_continuous_absorption_power(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    Vec3::new(1.0, 0.0, 0.0),
                    &[550.0],
                    &[10.0],
                    refused,
                ),
                Err(VolumeTransportError::Budget)
            );
        }
        gate.request();
        assert_eq!(
            m.deposit_continuous_absorption_power(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                &[550.0],
                &[10.0],
                options,
            ),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
