use super::*;
use crate::spectral::lift_rgb;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget, MassAbsorptionSpectrum};
use fs_exec::{CancelGate, StreamKey};

#[test]
fn reflected_thermal_nee_obeys_calibrated_physical_range() {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let density = [1.0];
        let temperature = [3000.0];
        let grid = AbsorbingGrid::try_new(
            &cx,
            VolumeGrid::new([1; 3], &density, [-0.5, -0.5, 2.0], [1.0; 3]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            Some(&temperature),
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 8,
            },
        )
        .unwrap();
        let point = Point3::new(0.0, 0.0, 0.0);
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let scene = Scene {
            primitives: Vec::new(),
            lights: Vec::new(),
            environment: None,
            camera: Camera {
                eye: point,
                forward: normal,
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.0,
            },
        };
        let material = Material::Lambertian {
            reflectance: lift_rgb([0.5; 3]),
        };
        let settings = Settings {
            width: 1,
            height: 1,
            spp: 1,
            max_depth: 1,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            seed: 31,
        };
        let rng = PathRng {
            pixel: 0,
            sample: 0,
            dim: 1,
            key: [31, 0],
            qmc_depth: 0,
            qmc_values: [0.0; PATH_QMC_DIMENSIONS],
            qmc_valid: false,
        };
        let evaluate = |fixtures: Option<&FixtureLighting<'_>>| {
            volume_emission::direct_contribution(
                &scene,
                &cx,
                None,
                &settings,
                0,
                PathMedium::Absorbing(&grid),
                None,
                point,
                volume_emission::Receiver::Surface {
                    material: &material,
                    frame: SurfaceFrame {
                        oriented: normal,
                        geometric: normal,
                        entering: true,
                    },
                    outgoing: normal,
                    brushed: None,
                    incident_glass: None,
                },
                &[450.0, 550.0, 650.0, 750.0],
                &rng,
                fixtures,
            )
        };
        let legacy = evaluate(None).unwrap();
        assert!(legacy.iter().all(|value| value.is_finite() && *value > 0.0));

        let mut calibrated = FixtureLighting {
            sources: &[],
            areas: area::AreaLighting::empty(),
            options: TransportOptions {
                max_distance_m: 10.0,
                ..TransportOptions::default()
            },
            medium_bounds_m: Some(grid.bounds_m()),
            scattering: false,
        };
        let admitted = evaluate(Some(&calibrated)).unwrap();
        for (actual, expected) in admitted.into_iter().zip(legacy) {
            assert!((actual - expected).abs() <= 1.0e-13 * expected);
        }

        calibrated.options.max_distance_m = 1.0;
        assert_eq!(
            evaluate(Some(&calibrated)),
            Err(TracerError::FixturePathRange),
            "reflected thermal source must not bypass the declared metre limit"
        );
        assert_eq!(evaluate(None).unwrap(), legacy);
    });
}
