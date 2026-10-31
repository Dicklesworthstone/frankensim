//! Independent Beer–Lambert energy and spatial-deposition checks.

use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::AbsorbingGridBudget;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction, ScatteringGrid};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 19,
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

fn budget(steps: u32) -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 4096,
        max_steps_per_segment: steps,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

fn balance<const N: usize>(report: &BeamPowerDeposition<N>) {
    for lane in 0..N {
        let absorbed: f64 = report.cells.iter().map(|c| c.absorbed_power_w[lane]).sum();
        let scattered: f64 = report.cells.iter().map(|c| c.scattered_power_w[lane]).sum();
        close(absorbed, report.absorbed_power_w[lane]);
        close(scattered, report.scattered_power_w[lane]);
        close(
            report.transmitted_power_w[lane] + absorbed + scattered,
            report.incident_power_w[lane],
        );
        assert!(report.cells.iter().all(|c| {
            c.absorbed_power_w[lane].is_finite()
                && c.absorbed_power_w[lane] >= 0.0
                && c.scattered_power_w[lane].is_finite()
                && c.scattered_power_w[lane] >= 0.0
        }));
    }
}

#[test]
fn g1_heterogeneous_cells_deposit_the_local_spectral_beer_lambert_loss() {
    with_cx(|cx, _| {
        let densities = [1.0, 3.0, 0.0, 2.0];
        let wavelengths = [380.0, 580.0, 780.0];
        let coefficients = [0.0, 0.4, 1.2];
        let incident = [3.0, 5.0, 7.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 4], &densities, [0.0; 3], [1.0, 1.0, 0.25]),
            super::super::MassAbsorptionSpectrum::table(&[
                (380.0, coefficients[0]),
                (580.0, coefficients[1]),
                (780.0, coefficients[2]),
            ])
            .unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        let report = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, -1.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
                &wavelengths,
                &incident,
                3,
            )
            .unwrap();
        assert_eq!(report.medium_identity, medium.content_identity());
        assert_eq!(report.wavelengths_nm, wavelengths);
        assert_eq!(report.incident_power_w, incident);
        assert_eq!(
            report
                .cells
                .iter()
                .map(|c| c.cell_index)
                .collect::<Vec<_>>(),
            vec![0, 1, 3]
        );
        let columns = [0.0, 0.25, 1.0, 1.0];
        for cell in &report.cells {
            let index = cell.cell_index;
            for lane in 0..3 {
                let before = (-coefficients[lane] * columns[index]).exp();
                let loss = -(-coefficients[lane] * densities[index] * 0.25).exp_m1();
                close(cell.absorbed_power_w[lane], incident[lane] * before * loss);
                assert_eq!(cell.scattered_power_w[lane], 0.0);
            }
        }
        for lane in 0..3 {
            let depth = coefficients[lane] * 1.5;
            close(
                report.transmitted_power_w[lane],
                incident[lane] * (-depth).exp(),
            );
            close(
                report.absorbed_power_w[lane],
                incident[lane] * -(-depth).exp_m1(),
            );
        }
        balance(&report);
    });
}

#[test]
fn g1_scattering_removes_beam_power_without_creating_heat() {
    with_cx(|cx, _| {
        let densities = [2.0];
        let wavelengths = [380.0, 580.0, 780.0];
        let incident = [10.0, 20.0, 30.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &densities, [0.0; 3], [1.0; 3]),
            super::super::MassAbsorptionSpectrum::table(&[
                (380.0, 0.0),
                (580.0, 0.25),
                (780.0, 1.0),
            ])
            .unwrap(),
            None,
            budget(8),
        )
        .unwrap();
        let absorption_identity = absorption.content_identity();
        let medium = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::table(&[(380.0, 1.0), (580.0, 0.75), (780.0, 0.0)]).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let report = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.0,
                &wavelengths,
                &incident,
                1,
            )
            .unwrap();
        assert_eq!(report.medium_identity, medium.content_identity());
        assert_ne!(report.medium_identity, absorption_identity);
        assert_eq!(report.cells.len(), 1);
        let loss = -(-2.0_f64).exp_m1();
        for lane in 0..3 {
            close(
                report.transmitted_power_w[lane],
                incident[lane] * (-2.0_f64).exp(),
            );
            close(report.absorbed_power_w[lane], [0.0, 5.0, 30.0][lane] * loss);
            close(
                report.scattered_power_w[lane],
                [10.0, 15.0, 0.0][lane] * loss,
            );
        }
        assert_eq!(
            report.absorbed_power_w[0], 0.0,
            "elastic scattering is not heat"
        );
        balance(&report);
    });
}

#[test]
fn g1_oblique_beam_heat_is_invariant_under_grid_and_segment_subdivision() {
    with_cx(|cx, _| {
        let coarse_density = [2.0];
        let fine_density = [2.0; 4];
        let coefficient = super::super::MassAbsorptionSpectrum::constant(0.8).unwrap();
        let coarse = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &coarse_density, [0.0; 3], [1.0; 3]),
            coefficient,
            None,
            budget(16),
        )
        .unwrap();
        let fine = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 4], &fine_density, [0.0; 3], [1.0, 1.0, 0.25]),
            coefficient,
            None,
            budget(16),
        )
        .unwrap();
        let origin = Point3::new(0.0, 0.5, 0.0);
        let direction = Vec3::new(0.6, 0.0, 0.8);
        let a = coarse
            .deposit_beam_power(cx, origin, direction, 2.0, &[550.0], &[4.2], 1)
            .unwrap();
        let b = fine
            .deposit_beam_power(cx, origin, direction, 2.0, &[550.0], &[4.2], 4)
            .unwrap();
        // Physical exit length is 1/0.8 = 1.25 m, not the one-metre z extent.
        close(a.transmitted_power_w[0], 4.2 * (-2.0_f64).exp());
        close(a.absorbed_power_w[0], b.absorbed_power_w[0]);
        close(a.transmitted_power_w[0], b.transmitted_power_w[0]);
        let first = fine
            .deposit_beam_power(cx, origin, direction, 0.625, &[550.0], &[4.2], 2)
            .unwrap();
        let second = fine
            .deposit_beam_power(
                cx,
                Point3::new(0.375, 0.5, 0.5),
                direction,
                0.625,
                &[550.0],
                &first.transmitted_power_w,
                2,
            )
            .unwrap();
        close(
            first.absorbed_power_w[0] + second.absorbed_power_w[0],
            b.absorbed_power_w[0],
        );
        close(second.transmitted_power_w[0], b.transmitted_power_w[0]);
        for (expected, actual) in b.cells.iter().zip(first.cells.iter().chain(&second.cells)) {
            assert_eq!(expected.cell_index, actual.cell_index);
            close(expected.absorbed_power_w[0], actual.absorbed_power_w[0]);
        }
        let reverse = fine
            .deposit_beam_power(
                cx,
                Point3::new(0.75, 0.5, 1.0),
                Vec3::new(-0.6, 0.0, -0.8),
                2.0,
                &[550.0],
                &[4.2],
                4,
            )
            .unwrap();
        assert_eq!(
            reverse
                .cells
                .iter()
                .map(|c| c.cell_index)
                .collect::<Vec<_>>(),
            vec![3, 2, 1, 0]
        );
        close(reverse.absorbed_power_w[0], b.absorbed_power_w[0]);
        for report in [&a, &b, &first, &second, &reverse] {
            balance(report);
        }
    });
}

#[test]
fn g0_transparent_missed_and_zero_power_beams_need_no_deposition_storage() {
    with_cx(|cx, _| {
        let densities = [1.0];
        for coefficient in [0.0, 1.0] {
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &densities, [0.0; 3], [1.0; 3]),
                super::super::MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                None,
                budget(8),
            )
            .unwrap();
            let direction = Vec3::new(0.0, 0.0, 1.0);
            for (origin, distance, power) in [
                (Point3::new(0.5, 0.5, 0.0), 1.0, 0.0),
                (Point3::new(2.0, 0.5, 0.0), 1.0, 3.0),
                (Point3::new(0.5, 0.5, 0.0), 0.0, 3.0),
            ] {
                let report = medium
                    .deposit_beam_power(cx, origin, direction, distance, &[550.0], &[power], 0)
                    .unwrap();
                assert!(report.cells.is_empty());
                assert_eq!(report.transmitted_power_w, [power]);
                balance(&report);
            }
            if coefficient == 0.0 {
                let report = medium
                    .deposit_beam_power(
                        cx,
                        Point3::new(0.5, 0.5, 0.0),
                        direction,
                        1.0,
                        &[550.0],
                        &[3.0],
                        0,
                    )
                    .unwrap();
                assert!(report.cells.is_empty());
                balance(&report);
            }
            let origin = Point3::new(0.5, 0.5, 0.0);
            for invalid in [-1.0, f64::INFINITY, f64::NAN] {
                assert_eq!(
                    medium.deposit_beam_power(cx, origin, direction, 1.0, &[550.0], &[invalid], 0),
                    Err(VolumeTransportError::BeamPower { lane: 0 })
                );
            }
            assert_eq!(
                medium.deposit_beam_power(cx, origin, direction, 1.0, &[f64::NAN], &[0.0], 0),
                Err(VolumeTransportError::Wavelength)
            );
            assert_eq!(
                medium.deposit_beam_power(
                    cx,
                    origin,
                    Vec3::new(0.0, 0.0, 2.0),
                    1.0,
                    &[550.0],
                    &[0.0],
                    0
                ),
                Err(VolumeTransportError::Ray)
            );
            assert_eq!(
                medium.deposit_beam_power(cx, origin, direction, 1.0, &[], &[], 0),
                Err(VolumeTransportError::Sample)
            );
        }
    });
}

#[test]
fn g4_deposition_refuses_output_traversal_and_cancellation_without_partial_results() {
    with_cx(|cx, gate| {
        let densities = [1.0; 2];
        for (steps, limit) in [(1, 2), (8, 1)] {
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], &densities, [0.0; 3], [1.0; 3]),
                super::super::MassAbsorptionSpectrum::constant(1.0).unwrap(),
                None,
                budget(steps),
            )
            .unwrap();
            assert_eq!(
                medium.deposit_beam_power(
                    cx,
                    Point3::new(0.5, 0.5, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    2.0,
                    &[550.0],
                    &[1.0],
                    limit,
                ),
                Err(VolumeTransportError::Budget)
            );
        }
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &densities, [0.0; 3], [1.0; 3]),
            super::super::MassAbsorptionSpectrum::constant(1.0).unwrap(),
            None,
            budget(8),
        )
        .unwrap();
        gate.request();
        assert_eq!(
            medium.deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                2.0,
                &[550.0],
                &[1.0],
                2,
            ),
            Err(VolumeTransportError::Cancelled)
        );
        assert_eq!(densities, [1.0; 2]);
    });
}

#[test]
fn g1_representable_watts_survive_underflow_in_fractional_opacity() {
    with_cx(|cx, _| {
        let densities = [1.0];
        let tiny = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &densities, [0.0; 3], [1.0, 1.0, 1e-30]),
            super::super::MassAbsorptionSpectrum::constant(1e-300).unwrap(),
            None,
            budget(8),
        )
        .unwrap();
        let tiny_report = tiny
            .deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                f64::INFINITY,
                &[550.0],
                &[1e300],
                1,
            )
            .unwrap();
        close(tiny_report.absorbed_power_w[0], 1e-30);
        let opaque = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &densities, [0.0; 3], [1.0; 3]),
            super::super::MassAbsorptionSpectrum::constant(1000.0).unwrap(),
            None,
            budget(8),
        )
        .unwrap();
        let opaque_report = opaque
            .deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.0,
                &[550.0],
                &[1e300],
                1,
            )
            .unwrap();
        // exp(-1000) underflows, while 1e300 * exp(-1000) is representable.
        close(
            opaque_report.transmitted_power_w[0],
            (1e300_f64.ln() - 1000.0).exp(),
        );
        assert!(opaque_report.transmitted_power_w[0] > 0.0);
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &densities, [0.0; 3], [1.0; 3]),
            super::super::MassAbsorptionSpectrum::constant(1e-320).unwrap(),
            None,
            budget(8),
        )
        .unwrap();
        let scattering = ScatteringGrid::try_new(
            cx,
            absorption,
            MassScatteringSpectrum::constant(1.0).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let report = scattering
            .deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.0,
                &[550.0],
                &[1e300],
                1,
            )
            .unwrap();
        close(
            report.absorbed_power_w[0],
            (1e300 * 1e-320) * -(-1.0_f64).exp_m1(),
        );
        balance(&report);
    });
}

#[test]
fn g1_parallel_decimal_cell_faces_keep_the_correct_heat_receiver() {
    with_cx(|cx, _| {
        let density: Vec<_> = (0..24).map(|i| 0.2 * (i + 1) as f64).collect();
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([3, 2, 4], &density, [-0.4, 0.7, -1.2], [0.3, 0.45, 0.22]),
            super::super::MassAbsorptionSpectrum::constant(0.3).unwrap(),
            None,
            budget(16),
        )
        .unwrap();
        // These are exact stored cell faces, but subtract/divide can put
        // their computed lattice coordinate just below the integer.
        for (origin, direction, length, expected_cell) in [
            (
                Point3::new(0.1, 1.15, -0.76),
                Vec3::new(-1.0, 0.0, 0.0),
                0.05,
                16,
            ),
            (
                Point3::new(0.05, 1.15, -0.76),
                Vec3::new(1.0, 0.0, 0.0),
                0.05,
                16,
            ),
            (
                Point3::new(0.05, 1.15, -0.76),
                Vec3::new(0.0, 0.0, -1.0),
                0.1,
                10,
            ),
            (
                Point3::new(0.05, 1.15, -0.76),
                Vec3::new(0.0, 0.0, 1.0),
                0.1,
                16,
            ),
        ] {
            let report = medium
                .deposit_beam_power(cx, origin, direction, length, &[550.0], &[2.0], 4)
                .unwrap();
            assert_eq!(report.cells.len(), 1);
            assert_eq!(report.cells[0].cell_index, expected_cell);
            let depth = 0.3 * density[expected_cell] * length;
            close(report.absorbed_power_w[0], 2.0 * -(-depth).exp_m1());
            close(report.transmitted_power_w[0], 2.0 * (-depth).exp());
            close(
                medium
                    .transmittance(cx, origin, direction, length, &[550.0])
                    .unwrap()[0],
                (-depth).exp(),
            );
            balance(&report);
        }
    });
}
