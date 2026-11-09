//! Prepared-cell Stokes histories retain the original native angular owner.

use super::*;
use crate::volumes::transport::scattering::field::{CellPhaseFunction, CellScatteringGrid};
use fs_blake3::ContentHash;

#[test]
fn g1_prepared_polarized_cells_replay_the_original_uniform_histories_for_both_estimators() {
    with_cx(|cx, _| {
        let density = [1.0, 2.0, 0.5];
        let phase = PhaseFunction::spectral_table(&[
            (380.0, PhaseFunction::Rayleigh),
            (780.0, PhaseFunction::rayleigh_depolarized(0.08).unwrap()),
        ])
        .unwrap();
        let medium = medium(cx, &density, 0.15, 0.85, phase);
        let absorption = [[0.15; 2]; 3];
        let scattering = [[0.85; 2]; 3];
        let wavelengths = [450.0, 680.0];
        let field = CellScatteringGrid::try_new(
            cx,
            medium.absorption_grid(),
            wavelengths,
            &absorption,
            &scattering,
            CellPhaseFunction::Uniform(medium.phase_ref().unwrap()),
            (380.0, 780.0),
            ContentHash([91; 32]),
        )
        .unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let frame =
            PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let source = [
            Stokes {
                i: 3.0,
                q: 3.0,
                u: 0.0,
                v: 0.0,
            },
            Stokes {
                i: 7.0,
                q: 0.0,
                u: 0.0,
                v: 7.0,
            },
        ];
        let config = options(256, 16);
        for estimator in [
            AbsorptionEstimator::Collision,
            AbsorptionEstimator::Continuous,
        ] {
            let mut original = match estimator {
                AbsorptionEstimator::Collision => medium
                    .deposit_polarized_multiple_scattering_power(
                        cx,
                        origin,
                        frame,
                        &wavelengths,
                        &source,
                        config,
                    ),
                AbsorptionEstimator::Continuous => medium
                    .deposit_polarized_continuous_absorption_power(
                        cx,
                        origin,
                        frame,
                        &wavelengths,
                        &source,
                        config,
                    ),
            }
            .unwrap();
            let prepared = field
                .deposit_polarized_power(cx, origin, frame, &source, config, estimator)
                .unwrap();
            assert_eq!(prepared.sampled_histories, [256; 2]);
            assert_ne!(original.medium_identity, prepared.medium_identity);
            original.medium_identity = prepared.medium_identity;
            assert_eq!(prepared, original);
            balance(&prepared);
        }
    });
}

#[test]
fn g0_prepared_polarization_admits_only_active_phase_laws_and_preserves_refusals() {
    with_cx(|cx, _| {
        let density = [1.0, 0.0, 1.0];
        let medium = medium(cx, &density, 0.4, 0.0, PhaseFunction::Rayleigh);
        let absorption = [[0.4]; 3];
        let scattering = [[0.0], [5.0], [1.0]];
        let phases = [
            PhaseFunction::Isotropic,
            PhaseFunction::HenyeyGreenstein { g: 0.5 },
            PhaseFunction::Rayleigh,
        ];
        let indices = [0, 1, 2];
        let make = |coefficients| {
            CellScatteringGrid::try_new(
                cx,
                medium.absorption_grid(),
                [550.0],
                &absorption,
                coefficients,
                CellPhaseFunction::Indexed {
                    palette: &phases,
                    cell_indices: &indices,
                },
                (380.0, 780.0),
                ContentHash([92; 32]),
            )
        };
        let field = make(&scattering).unwrap();
        let origin = Point3::new(-1.0, 0.5, 0.5);
        let frame =
            PolarizationFrame::try_new(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap();
        let source = [Stokes::unpolarized(1.0)];
        let config = options(32, 16);
        let report = field
            .deposit_polarized_power(
                cx,
                origin,
                frame,
                &source,
                config,
                AbsorptionEstimator::Continuous,
            )
            .unwrap();
        balance(&report);
        assert!(report.cells.iter().any(|cell| cell.cell_index == 0));
        assert!(report.cells.iter().all(|cell| cell.cell_index != 1));
        let active_scalar = [[1.0], [0.0], [0.0]];
        assert_eq!(
            make(&active_scalar).unwrap().deposit_polarized_power(
                cx,
                origin,
                frame,
                &source,
                config,
                AbsorptionEstimator::Continuous,
            ),
            Err(PolarizedScatteringError::UnsupportedPhase),
        );
        assert_eq!(
            field.deposit_polarized_power_with_checkpoint(
                &mut || Err(VolumeTransportError::Cancelled),
                origin,
                frame,
                &source,
                config,
                AbsorptionEstimator::Continuous,
            ),
            Err(VolumeTransportError::Cancelled.into()),
        );
        assert_eq!(
            field.deposit_polarized_power(
                cx,
                origin,
                frame,
                &[Stokes {
                    i: 1.0,
                    q: 1.1,
                    u: 0.0,
                    v: 0.0
                }],
                config,
                AbsorptionEstimator::Continuous,
            ),
            Err(BeamError::InvalidStokes.into()),
        );
    });
}
