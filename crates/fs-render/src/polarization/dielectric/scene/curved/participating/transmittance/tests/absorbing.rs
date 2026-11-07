//! Deterministic native true absorption on the original source replay.

use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

#[test]
fn finite_absorbing_source_replay_combines_native_cells_and_scattering_once() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, Point3::new(0.0, 0.0, 4.0))]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let density = [2.0, 1.0, 3.0];
        let temperatures = [1000.0, 1300.0, 1600.0];
        let absorber = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 3], &density, [-1.0, -1.0, 0.2], [2.0, 2.0, 0.25]),
            MassAbsorptionSpectrum::constant(0.5).unwrap(),
            Some(&temperatures),
            AbsorbingGridBudget {
                max_grid_cells: 3,
                max_steps_per_segment: 3,
            },
        )
        .unwrap();
        let grids = [Some(&absorber)];
        let mass = MassScatteringSpectrum::constant(0.4).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let receiver = BundleConnectionTarget::try_volume_receiver(
            0,
            Point3::new(0.0, 0.0, 1.1),
            Vec3::new(0.0, 0.0, -1.0),
            2e-7,
        )
        .unwrap();
        for law in [
            None,
            Some(CurvedRayleighMedium::try_new(&mass, 1.0, &phase).unwrap()),
        ] {
            let laws = [law];
            let transport = CurvedRayleighScene::try_new(cx, &scene, &laws)
                .unwrap()
                .with_absorbing_grids(cx, &grids, 3)
                .unwrap();
            let source = BundleConnectionLaunch {
                normalized_source_stokes: Stokes {
                    i: 1.0,
                    q: 0.2,
                    u: -0.3,
                    v: 0.4,
                },
                ..launch(Vec3::new(0.0, 0.0, 1.0))
            };
            let measured = transport
                .connect_bundle_with_transmittance_samples(
                    cx,
                    source,
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                    CurvedTransmittanceOptions {
                        max_free_flight_draws: 0,
                        max_tracking_candidates: 0,
                    },
                    |_| panic!("uniform scattering and deterministic absorption need no variates"),
                )
                .unwrap();
            let scatter_depth = if law.is_some() { 0.44 } else { 0.0 };
            near(measured.log_absorption_survival(), -0.75, 1e-12);
            near(measured.log_scattering_survival(), -scatter_depth, 2e-8);
            near(measured.log_source_survival(), -0.75 - scatter_depth, 2e-8);
            assert_eq!(measured.absorption_cell_visits(), 3);
            assert_eq!(measured.free_flight_draws(), 0);
            assert_eq!(measured.tracking_candidates(), 0);
            assert_eq!(
                measured.trace_evaluations(),
                measured.connection().trace_evaluations() + 1
            );
            assert_eq!(
                measured
                    .connection()
                    .path()
                    .bundle()
                    .propagation()
                    .propagation()
                    .optical_depth,
                0.0,
                "finite absorption is retained separately, not applied again by the geometry solve"
            );
            assert!(measured.connection().path().is_volume_receiver());
        }
    });
}

#[test]
fn finite_absorbing_source_keeps_large_finite_depth_and_enforces_actual_cell_work() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, Point3::new(0.0, 0.0, 4.0))]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let density = [1.0, 1.0];
        let absorber = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.5]),
            MassAbsorptionSpectrum::constant(800.0).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 2,
                max_steps_per_segment: 2,
            },
        )
        .unwrap();
        let grids = [Some(&absorber)];
        let laws = [None];
        let receiver = BundleConnectionTarget::try_volume_receiver(
            0,
            Point3::new(0.0, 0.0, 1.5),
            Vec3::new(0.0, 0.0, -1.0),
            2e-7,
        )
        .unwrap();
        let connect = |cap| {
            CurvedRayleighScene::try_new(cx, &scene, &laws)
                .unwrap()
                .with_absorbing_grids(cx, &grids, cap)
                .unwrap()
                .connect_bundle_with_transmittance_samples(
                    cx,
                    launch(Vec3::new(0.0, 0.0, 1.0)),
                    receiver,
                    options(),
                    BundleBranchPolicy::Sequence(&[]),
                    CurvedTransmittanceOptions::default(),
                    |_| panic!("no sampled source extinction"),
                )
        };
        assert!(connect(1).is_err());
        let measured = connect(2).unwrap();
        near(measured.log_source_survival(), -800.0, 2e-10);
        assert!(!measured.ratio_tracking_zero());
        assert_eq!(measured.log_source_survival().exp(), 0.0);
        let bright_source = (1e300_f64.ln() + measured.log_source_survival()).exp();
        assert!(bright_source > 0.0 && bright_source.is_finite());
        assert_eq!(measured.absorption_cell_visits(), 2);
    });
}
