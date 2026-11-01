use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, DielectricSurface,
    GlassProvenance,
};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::bundle::PlanarIrradianceReceiver;
use crate::photometry::source::bundle::connection::{
    ConnectedIrradianceError, ConnectedIrradianceRequest,
};
use crate::photometry::source::{
    AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, PhotometricSource, SPECTRAL_BINS,
    SourceSpectrum, photopic_weight_nm, wavelength_nm,
};
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::bundle::connection::{
    BundleBranchPolicy, BundleConnectionError, BundleConnectionOptions,
};
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::photometry::SpectralOptics;
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

const DIRECT: [DielectricEvent; 2] = [DielectricEvent::Transmission; 2];
const FIRST_GHOST: [DielectricEvent; 4] = [
    DielectricEvent::Transmission,
    DielectricEvent::Reflection,
    DielectricEvent::Reflection,
    DielectricEvent::Transmission,
];
const SECOND_GHOST: [DielectricEvent; 6] = [
    DielectricEvent::Transmission,
    DielectricEvent::Reflection,
    DielectricEvent::Reflection,
    DielectricEvent::Reflection,
    DielectricEvent::Reflection,
    DielectricEvent::Transmission,
];

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn source() -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0]).unwrap(),
        1200.0,
        0.75,
        AngularDistribution::cosine_power(2.0, core::f64::consts::FRAC_PI_2).unwrap(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn glass(ior: CauchyIor, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        ior,
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn plane(id: u64, z: f64, positive: bool) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("multipath-slab-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-6.0, -6.0, z], [6.0, -6.0, z], [0.0, 6.0, z]],
            vec![if positive { [0, 1, 2] } else { [0, 2, 1] }],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn geometry() -> InstanceScene {
    InstanceScene::try_new(vec![
        plane(1, 0.0, false),
        plane(2, 1.0, true),
        plane(9, 2.0, true),
    ])
    .unwrap()
}

fn bindings() -> Vec<SurfaceBinding> {
    [1, 2]
        .map(|object_id| SurfaceBinding {
            object_id,
            interaction: Interaction::Dielectric {
                outside: 1,
                inside: 0,
                surface: DielectricSurface::SMOOTH,
            },
        })
        .into_iter()
        .chain([SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Detector {
                reference_axis: Vec3::new(0.0, 1.0, 0.0),
            },
        }])
        .collect()
}

fn request() -> ConnectedIrradianceRequest {
    ConnectedIrradianceRequest {
        detector_object_id: 9,
        receiver: PlanarIrradianceReceiver::try_new(
            Point3::new(0.2, 0.0, 2.0),
            Vec3::new(0.0, 0.0, -1.0),
            2e-7,
        )
        .unwrap(),
        initial_direction: Vec3::new(0.06, 0.0, 1.0),
        source_transverse_axis: Vec3::new(1.0, 0.0, 0.0),
        initial_medium: 1,
        connection: BundleConnectionOptions {
            scene: CurvedBundleSettings {
                trace: CurvedTraceSettings {
                    max_surface_events: 12,
                    max_distance_m: 10.0,
                    max_total_steps: 16_384,
                    max_query_work: 1_000_000,
                    event: GradientEventOptions {
                        integration: GradientIndexOptions {
                            position_tolerance_m: 1e-11,
                            direction_tolerance: 1e-12,
                            optical_path_tolerance_m: 1e-11,
                            optical_depth_tolerance: 1e-12,
                            relative_tolerance: 1e-11,
                            ..GradientIndexOptions::default()
                        },
                        surface_tolerance_m: 1e-8,
                        ..GradientEventOptions::default()
                    },
                },
                ..CurvedBundleSettings::default()
            },
            ..BundleConnectionOptions::default()
        },
    }
}

fn family() -> [ConnectedIrradianceBranch<'static>; 3] {
    [
        ConnectedIrradianceBranch {
            initial_direction: Vec3::new(0.06, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&DIRECT),
        },
        ConnectedIrradianceBranch {
            initial_direction: Vec3::new(0.04, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&FIRST_GHOST),
        },
        ConnectedIrradianceBranch {
            initial_direction: Vec3::new(0.025, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&SECOND_GHOST),
        },
    ]
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-15),
        "{actual:.16e} != {expected:.16e}"
    );
}

/// Independent parallel-plate optics. Bisection solves Snell's displacement
/// equation; its closed-form derivative gives dA/dOmega. No production event,
/// Newton, derivative or Fresnel code supplies this reference.
fn plate_reference(index: f64, passes: u32) -> (Vec3, [f64; 2]) {
    let mut lo = 0.0_f64;
    let mut hi = 0.4_f64;
    let m = f64::from(passes);
    for _ in 0..80 {
        let theta = 0.5 * (lo + hi);
        let st = theta.sin() / index;
        let ct = (1.0 - st * st).sqrt();
        let radius = 2.0 * theta.tan() + m * st / ct;
        if radius < 0.2 {
            lo = theta;
        } else {
            hi = theta;
        }
    }
    let theta = 0.5 * (lo + hi);
    let si = theta.sin();
    let ci = theta.cos();
    let ct = (1.0 - (si / index).powi(2)).sqrt();
    let radial_derivative = 2.0 / (ci * ci) + m * ci / (index * ct.powi(3));
    let area = 0.2 * radial_derivative / si;
    let rs = ((ci - index * ct) / (ci + index * ct)).powi(2);
    let rp = ((index * ci - ct) / (index * ci + ct)).powi(2);
    let beer_over_area = (-0.15 * m / ct).exp() / area;
    let reflections = (passes - 1) as i32;
    (
        Vec3::new(si, 0.0, ci),
        [
            (1.0 - rp).powi(2) * rp.powi(reflections) * beer_over_area,
            (1.0 - rs).powi(2) * rs.powi(reflections) * beer_over_area,
        ],
    )
}

#[test]
fn g1_eighty_bin_direct_and_internal_reflection_images_match_independent_plate_optics() {
    with_cx(|cx, _| {
        let source = source();
        let geometry = geometry();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let slab = glass(CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&slab),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let controls = request();
        let result = source
            .connect_irradiance_branches(cx, &scene, controls, &family())
            .unwrap();
        let mut expected_radiant = [0.0; 3];
        let mut expected_lux = [0.0; 3];
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let mut expected = 0.0;
            for (branch, passes) in [1, 3, 5].into_iter().enumerate() {
                let (direction, channels) = plate_reference(index, passes);
                let value = source
                    .spectral_intensity_w_sr_nm(direction, wavelength)
                    .unwrap()
                    * 0.5
                    * (channels[0] + channels[1]);
                expected += value;
                expected_radiant[branch] += value * BIN_WIDTH_NM;
                expected_lux[branch] += value * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
            }
            close(result.spectral_irradiance_w_m2_nm()[bin], expected, 5e-6);
        }
        assert_eq!(result.branches().len(), 3);
        for (index, branch) in result.branches().iter().enumerate() {
            close(
                branch.reading.radiant_irradiance_w_m2,
                expected_radiant[index],
                5e-6,
            );
            close(branch.reading.illuminance_lux, expected_lux[index], 5e-6);
            assert!(branch.max_endpoint_error_m <= controls.receiver.position_tolerance_m());
            assert!(branch.trace_evaluations >= SPECTRAL_BINS as u64);
            assert!(branch.max_contact_resolution_m > 0.0);
        }
        assert_ne!(
            result.branches()[0].signature,
            result.branches()[1].signature
        );
        assert_ne!(
            result.branches()[1].signature,
            result.branches()[2].signature
        );
        close(
            result.reading().radiant_irradiance_w_m2,
            expected_radiant.iter().sum(),
            5e-6,
        );
        close(
            result.reading().illuminance_lux,
            expected_lux.iter().sum(),
            5e-6,
        );
        assert!(result.reading().illuminance_lux > result.branches()[0].reading.illuminance_lux);
        assert!(expected_lux[1] > 100.0 * expected_lux[2]);
        assert_eq!(
            result.trace_evaluations(),
            result.branches().iter().map(|b| b.trace_evaluations).sum()
        );
        assert_eq!(
            result.field_evaluations(),
            result.branches().iter().map(|b| b.field_evaluations).sum()
        );
        assert_eq!(
            result.query_work(),
            result.branches().iter().map(|b| b.query_work).sum()
        );
    });
}

#[test]
fn g1_polarized_multipath_sum_preserves_upstream_filter_losses_and_distinct_sp_reflections() {
    with_cx(|cx, _| {
        let source = source();
        let geometry = geometry();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let slab = glass(CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&slab),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let parallel =
            core::array::from_fn(|bin| 0.2 + 0.5 * bin as f64 / (SPECTRAL_BINS - 1) as f64);
        let perpendicular = [0.1; SPECTRAL_BINS];
        let neutral =
            core::array::from_fn(|bin| 0.6 - 0.2 * bin as f64 / (SPECTRAL_BINS - 1) as f64);
        let filtered = source
            .clone()
            .with_unpolarized_optics()
            .through(&SpectralOptics::linear_diattenuator(0.0, parallel, perpendicular).unwrap())
            .unwrap()
            .through(&SpectralOptics::neutral(neutral).unwrap())
            .unwrap();
        let result = filtered
            .connect_irradiance_branches(cx, &scene, request(), &family())
            .unwrap();
        let mut expected_lux = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let wavelength = wavelength_nm(bin);
            let index = 1.45 + 0.012 / (wavelength * 1e-3).powi(2);
            let expected: f64 = [1, 3, 5]
                .into_iter()
                .map(|passes| {
                    let (direction, channels) = plate_reference(index, passes);
                    // Projected world-x is p for every xz-plane launch. Each ghost
                    // retains separate reflected channel powers and source losses.
                    source
                        .spectral_intensity_w_sr_nm(direction, wavelength)
                        .unwrap()
                        * 0.5
                        * neutral[bin]
                        * (parallel[bin] * channels[0] + perpendicular[bin] * channels[1])
                })
                .sum();
            close(result.spectral_irradiance_w_m2_nm()[bin], expected, 5e-6);
            expected_lux += expected * photopic_weight_nm(bin) * LUMINOUS_EFFICACY_LM_W;
        }
        close(result.reading().illuminance_lux, expected_lux, 5e-6);
        assert!(filtered.intensity_source().flux_lm() < 0.3 * source.flux_lm());
    });
}

#[test]
fn g4_multipath_admission_and_failed_members_never_publish_partial_or_duplicate_light() {
    with_cx(|cx, gate| {
        let source = source();
        let geometry = geometry();
        let air = glass(CauchyIor::try_constant(1.0).unwrap(), 0.0);
        let slab = glass(CauchyIor::try_new(1.45, 0.012, 0.0).unwrap(), 0.15);
        let media = [
            SceneMedium::Homogeneous(&slab),
            SceneMedium::Homogeneous(&air),
        ];
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        assert!(matches!(
            source.connect_irradiance_branches(cx, &scene, request(), &[]),
            Err(MultipathIrradianceError::EmptyBranches)
        ));
        let too_many = [family()[0]; MAX_IRRADIANCE_BRANCHES + 1];
        assert!(matches!(
            source.connect_irradiance_branches(cx, &scene, request(), &too_many),
            Err(MultipathIrradianceError::TooManyBranches)
        ));
        let duplicate = [
            ConnectedIrradianceBranch {
                initial_direction: Vec3::new(0.06, 0.0, 1.0),
                policy: BundleBranchPolicy::TransmitOrReflectTir,
            },
            ConnectedIrradianceBranch {
                initial_direction: Vec3::new(0.08, 0.0, 1.0),
                policy: BundleBranchPolicy::Sequence(&DIRECT),
            },
        ];
        assert!(matches!(
            source.connect_irradiance_branches(cx, &scene, request(), &duplicate),
            Err(MultipathIrradianceError::DuplicateBranch {
                first: 0,
                second: 1
            })
        ));
        let failed_member = [
            family()[0],
            ConnectedIrradianceBranch {
                initial_direction: Vec3::new(0.04, 0.0, 1.0),
                policy: BundleBranchPolicy::Sequence(&[]),
            },
        ];
        assert!(matches!(
            source.connect_irradiance_branches(cx, &scene, request(), &failed_member),
            Err(MultipathIrradianceError::Branch {
                branch: 1,
                error: ConnectedIrradianceError::Connection { bin: 0, .. }
            })
        ));
        gate.request();
        assert!(matches!(
            source.connect_irradiance_branches(cx, &scene, request(), &family()),
            Err(MultipathIrradianceError::Branch {
                branch: 0,
                error: ConnectedIrradianceError::Connection {
                    bin: 0,
                    error: BundleConnectionError::Cancelled,
                }
            })
        ));
    });
}
