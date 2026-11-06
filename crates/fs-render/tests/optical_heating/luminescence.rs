//! G1/G3 physical line radiance and G0/G4 admission/cancellation boundaries.

use super::*;
use fs_blake3::ContentHash;
use fs_render::volumes::heating::fluorescence::{FluorescenceSource, FluorescentOpticalHeating};
use fs_render::volumes::transport::VolumeTransportError;
use fs_render::volumes::transport::luminescence::LuminescenceGrid;

const FOUR_PI: f64 = 4.0 * core::f64::consts::PI;
const SOURCE: ContentHash = ContentHash([109; 32]);

fn relative(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn transparent_emission_uses_original_volume_solid_angle_and_integrated_bin_units() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 1], &density, [0.0; 3], [2.0, 3.0, 4.0]),
            MassAbsorptionSpectrum::constant(0.0).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 1,
            },
        )
        .unwrap();
        assert!(medium.is_vacuum());
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            SOURCE,
            [510.0, 650.0],
            &[(0, [24.0 * FOUR_PI, 48.0 * FOUR_PI])],
            1,
        )
        .unwrap();
        let ray = view
            .transport(
                cx,
                Point3::new(-1.0, 1.5, 2.0),
                Vec3::new(1.0, 0.0, 0.0),
                f64::INFINITY,
            )
            .unwrap();
        assert_eq!(ray.transmittance, [1.0; 2]);
        relative(ray.emitted_radiance_w_m2_sr[0], 2.0);
        relative(ray.emitted_radiance_w_m2_sr[1], 4.0);
        let partial = view
            .transport(
                cx,
                Point3::new(1.0, 1.5, 2.0),
                Vec3::new(1.0, 0.0, 0.0),
                0.25,
            )
            .unwrap();
        relative(partial.emitted_radiance_w_m2_sr[0], 0.25);
        assert_eq!(view.source_identity(), SOURCE);
        assert_eq!(view.wavelengths_nm(), &[510.0, 650.0]);
    });
}

#[test]
fn homogeneous_slab_agrees_with_independent_source_integral_at_all_optical_depths() {
    with_cx(|cx| {
        let density = [2.0];
        for alpha in [0.0_f64, 1e-20, core::f64::consts::LN_2, 1e4] {
            let medium = grid(
                cx,
                &density,
                MassAbsorptionSpectrum::constant(alpha / 2.0).unwrap(),
            );
            let view =
                LuminescenceGrid::try_new(cx, &medium, SOURCE, [600.0], &[(0, [FOUR_PI])], 1)
                    .unwrap();
            let ray = view
                .transport(
                    cx,
                    Point3::new(-1.0, 0.5, 0.5),
                    Vec3::new(1.0, 0.0, 0.0),
                    3.0,
                )
                .unwrap();
            let integral = if alpha == 0.0 {
                1.0
            } else {
                -(-alpha).exp_m1() / alpha
            };
            relative(ray.emitted_radiance_w_m2_sr[0], integral);
            relative(ray.transmittance[0], (-alpha).exp());
        }
    });
}

#[test]
fn observer_direction_changes_foreground_shielding_without_changing_source_power() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let alpha = core::f64::consts::LN_2;
        let medium = grid(
            cx,
            &density,
            MassAbsorptionSpectrum::constant(alpha).unwrap(),
        );
        let powers = [(1, [FOUR_PI]), (0, [0.0])];
        let view = LuminescenceGrid::try_new(cx, &medium, SOURCE, [600.0], &powers, 2).unwrap();
        let canonical = LuminescenceGrid::try_new(
            cx,
            &medium,
            SOURCE,
            [600.0],
            &[(0, [0.0]), (1, [FOUR_PI])],
            2,
        )
        .unwrap();
        assert_eq!(view.content_identity(), canonical.content_identity());
        let left = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                4.0,
            )
            .unwrap();
        let right = view
            .transport(
                cx,
                Point3::new(3.0, 0.5, 0.5),
                Vec3::new(-1.0, 0.0, 0.0),
                4.0,
            )
            .unwrap();
        relative(left.transmittance[0], 0.25);
        relative(right.transmittance[0], 0.25);
        relative(left.emitted_radiance_w_m2_sr[0], 0.25 / alpha);
        relative(right.emitted_radiance_w_m2_sr[0], 0.5 / alpha);
        let changed = LuminescenceGrid::try_new(
            cx,
            &medium,
            SOURCE,
            [600.0],
            &[(0, [0.0]), (1, [2.0 * FOUR_PI])],
            2,
        )
        .unwrap();
        assert_ne!(view.content_identity(), changed.content_identity());
    });
}

#[test]
fn cell_refinement_preserves_uniform_emissivity_and_observed_radiance() {
    with_cx(|cx| {
        let coarse_density = [1.0];
        let fine_density = [1.0, 1.0];
        let absorption = MassAbsorptionSpectrum::constant(0.7).unwrap();
        let coarse = grid(cx, &coarse_density, absorption);
        let fine = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2, 1, 1], &fine_density, [0.0; 3], [0.5, 1.0, 1.0]),
            absorption,
            None,
            AbsorbingGridBudget {
                max_grid_cells: 2,
                max_steps_per_segment: 2,
            },
        )
        .unwrap();
        let coarse_view =
            LuminescenceGrid::try_new(cx, &coarse, SOURCE, [550.0], &[(0, [FOUR_PI])], 1).unwrap();
        let fine_view = LuminescenceGrid::try_new(
            cx,
            &fine,
            SOURCE,
            [550.0],
            &[(0, [0.5 * FOUR_PI]), (1, [0.5 * FOUR_PI])],
            2,
        )
        .unwrap();
        let observe = |view: &LuminescenceGrid<'_, '_, 1>| {
            view.transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
            )
            .unwrap()
        };
        let a = observe(&coarse_view);
        let b = observe(&fine_view);
        relative(a.emitted_radiance_w_m2_sr[0], b.emitted_radiance_w_m2_sr[0]);
        relative(a.transmittance[0], b.transmittance[0]);
    });
}

#[test]
fn native_fluorescence_powers_reach_physical_radiance_without_double_counting_heat() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let spectrum = MassAbsorptionSpectrum::source_table(
            &[
                (380.0, core::f64::consts::LN_2),
                (500.0, core::f64::consts::LN_2),
                (600.0, 0.0),
                (800.0, 0.0),
            ],
            ContentHash([110; 32]),
        )
        .unwrap();
        let medium = grid(cx, &density, spectrum);
        let deposited = medium
            .deposit_beam_power(
                cx,
                Point3::new(3.0, 0.5, 0.5),
                Vec3::new(-1.0, 0.0, 0.0),
                4.0,
                &[400.0],
                &[16.0],
                2,
            )
            .unwrap();
        let source = FluorescenceSource::try_new(
            SOURCE,
            (380.0, 500.0),
            (290.0, 800.0),
            0.75,
            [600.0, 800.0],
            [1.0, 1.0],
        )
        .unwrap();
        let thermal = thermal(
            vec![cell(1, 1.0, 1.0, 0.0), cell(2, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let heating = FluorescentOpticalHeating::try_new(
            &medium,
            source,
            deposited,
            &[(0, CellId(1)), (1, CellId(2))],
            &thermal,
            2,
            1e-12,
            || Ok(()),
        )
        .unwrap();
        let view = LuminescenceGrid::try_new(
            cx,
            &medium,
            source.source_identity(),
            *source.wavelengths_nm(),
            heating.cell_emission_w(),
            2,
        )
        .unwrap();
        let ray = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                4.0,
            )
            .unwrap();
        // Pump absorption is 12 W. Half of emitted photons are at each node;
        // their radiant power includes the 400/lambda Stokes energy ratio.
        relative(
            ray.emitted_radiance_w_m2_sr[0],
            12.0 * 0.75 * 0.5 * 400.0 / 600.0 / FOUR_PI,
        );
        relative(
            ray.emitted_radiance_w_m2_sr[1],
            12.0 * 0.75 * 0.5 * 400.0 / 800.0 / FOUR_PI,
        );
        let observed_power = FOUR_PI * ray.emitted_radiance_w_m2_sr.iter().sum::<f64>();
        relative(observed_power, heating.power().escaped_w.iter().sum());
        relative(observed_power + heating.pump_power().thermal_input_w, 12.0);
    });
}

#[test]
fn bright_emission_survives_underflowing_transmittance_and_overflowing_optical_depth() {
    with_cx(|cx| {
        let density = [1.0, 1.0];
        let medium = grid(
            cx,
            &density,
            MassAbsorptionSpectrum::constant(800.0).unwrap(),
        );
        let view =
            LuminescenceGrid::try_new(cx, &medium, SOURCE, [600.0], &[(1, [FOUR_PI * 1e300])], 1)
                .unwrap();
        let observed = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                4.0,
            )
            .unwrap();
        relative(
            observed.emitted_radiance_w_m2_sr[0],
            (1e300_f64.ln() - 800.0_f64.ln() - 800.0).exp(),
        );
        assert_eq!(observed.transmittance[0], 0.0);

        let dense = [1e200];
        let opaque = grid(cx, &dense, MassAbsorptionSpectrum::constant(1e200).unwrap());
        let view =
            LuminescenceGrid::try_new(cx, &opaque, SOURCE, [600.0], &[(0, [FOUR_PI * 1e300])], 1)
                .unwrap();
        let observed = view
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
            )
            .unwrap();
        relative(observed.emitted_radiance_w_m2_sr[0], 1e-100);
        assert_eq!(observed.transmittance[0], 0.0);
    });
}

#[test]
fn malformed_mapping_power_spectrum_and_preparation_budget_refuse() {
    with_cx(|cx| {
        let density = [1.0];
        let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(0.0).unwrap());
        let make = |powers: &[(usize, [f64; 1])], max_cells| {
            LuminescenceGrid::try_new(cx, &medium, SOURCE, [600.0], powers, max_cells)
        };
        for powers in [&[(1, [1.0])][..], &[(0, [1.0]), (0, [2.0])][..]] {
            assert!(matches!(
                make(powers, 2),
                Err(VolumeTransportError::Palette)
            ));
        }
        for power in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                make(&[(0, [power])], 1),
                Err(VolumeTransportError::BeamPower { lane: 0 })
            ));
        }
        assert!(matches!(
            make(&[(0, [1.0])], 0),
            Err(VolumeTransportError::Budget)
        ));
        assert!(matches!(
            LuminescenceGrid::try_new(cx, &medium, ContentHash([0; 32]), [600.0], &[(0, [1.0])], 1),
            Err(VolumeTransportError::SourceIdentity)
        ));
        for wavelengths in [
            [600.0, 500.0],
            [500.0, 500.0],
            [500.0, 900.0],
            [f64::NAN, 600.0],
        ] {
            assert!(matches!(
                LuminescenceGrid::try_new(cx, &medium, SOURCE, wavelengths, &[(0, [1.0; 2])], 1),
                Err(VolumeTransportError::Wavelength)
            ));
        }
        let empty = make(&[], 0).unwrap();
        let ray = empty
            .transport(
                cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                3.0,
            )
            .unwrap();
        assert_eq!(ray.emitted_radiance_w_m2_sr, [0.0]);
        assert_eq!(ray.transmittance, [1.0]);
    });
}

#[test]
fn native_ray_step_budget_and_cancellation_publish_no_partial_observation() {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 109,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let density = [1.0, 1.0];
        let medium = AbsorbingGrid::try_new(
            &cx,
            VolumeGrid::new([2, 1, 1], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(0.0).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 2,
                max_steps_per_segment: 1,
            },
        )
        .unwrap();
        let view =
            LuminescenceGrid::try_new(&cx, &medium, SOURCE, [600.0], &[(0, [1.0]), (1, [2.0])], 2)
                .unwrap();
        let trace = || {
            view.transport(
                &cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(1.0, 0.0, 0.0),
                4.0,
            )
        };
        assert_eq!(trace(), Err(VolumeTransportError::Budget));
        assert_eq!(
            view.transport(
                &cx,
                Point3::new(-1.0, 0.5, 0.5),
                Vec3::new(2.0, 0.0, 0.0),
                4.0
            ),
            Err(VolumeTransportError::Ray)
        );
        gate.request();
        assert_eq!(trace(), Err(VolumeTransportError::Cancelled));
        assert!(matches!(
            LuminescenceGrid::try_new(&cx, &medium, SOURCE, [600.0], &[], 0),
            Err(VolumeTransportError::Cancelled)
        ));
    });
}
