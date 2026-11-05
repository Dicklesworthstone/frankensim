use super::*;
use crate::camera::PhysicalCamera;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexOptions};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::photometry::scene::{
    NativeGeometryBinding, NativeOpticalBinding, NativeOpticalRole,
    camera::NativeRadianceCalibration,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::scene::curved::bundle::{
    CurvedBundleSceneError, CurvedBundleSettings,
};
use crate::polarization::dielectric::scene::curved::camera::{CurvedCameraScene, RadianceEndpoint};
use crate::spectral::{LiftedSpectrum, SpectralReflectance};
use crate::tracer::{Camera, Primitive, Scene, Shape};
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3, fixtures::SphereChart};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};

const MEDIUM: NativeOpticalMedium = NativeOpticalMedium::InteriorOf(0);
const SPECTRUM: LiftedSpectrum = LiftedSpectrum {
    c: [3.0, -2.0, -0.25],
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn glass(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn owner(cx: &Cx<'_>, center_x: f64, width: f64, detector: bool) -> NativeOpticalScene {
    let scene = Scene {
        primitives: vec![
            Primitive {
                shape: Shape::Chart(Box::new(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                })),
                material: Material::Dielectric {
                    glass: glass(1.5, 0.2),
                    surface: DielectricSurface::SMOOTH,
                },
                emission: None,
            },
            Primitive {
                shape: Shape::Mesh(TriMesh::new(
                    vec![
                        [center_x - width, -width, 2.0],
                        [center_x + width, -width, 2.0],
                        [center_x, 2.0 * width, 2.0],
                    ],
                    vec![[0, 1, 2]],
                )),
                material: Material::SpectralDiffuse {
                    reflectance: SpectralReflectance::constant(0.0).unwrap(),
                },
                emission: (!detector).then_some((SPECTRUM, 2.0)),
            },
        ],
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 0.0),
            forward: Vec3::new(0.0, 0.0, 1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.05,
        },
    };
    let bindings = [1_u64, 9].map(|id| NativeOpticalBinding {
        geometry: NativeGeometryBinding::Unplaced {
            object_id: id,
            geometry_identity: hash_domain("native-thermo-optic-test", &id.to_le_bytes()),
        },
        role: if id == 9 && detector {
            NativeOpticalRole::Detector {
                reference_axis: Vec3::new(1.0, 0.0, 0.0),
            }
        } else {
            NativeOpticalRole::Material {
                outside: NativeOpticalMedium::Ambient,
            }
        },
    });
    NativeOpticalScene::try_from_scene(cx, scene, glass(1.0, 0.1), &bindings, 2).unwrap()
}

fn thermal_law(d0: f64, e0: f64) -> SchottThermoOpticLaw {
    // Synthetic coefficients test exact transport, not a measured glass fit.
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [d0, 0.0, 0.0],
            e: [e0, 0.0],
            resonance_um: 0.0,
        },
        300.0,
        [200.0, 450.0],
        [380.0, 780.0],
        hash_domain("native-thermo-optic-test", b"synthetic-law"),
    )
    .unwrap()
}

struct Temperature {
    center_k: f64,
    gradient: [f64; 3],
    diagonal: [f64; 3],
    bound_m: f64,
}

impl Temperature {
    fn uniform() -> Self {
        Self {
            center_k: 300.0,
            gradient: [0.0; 3],
            diagonal: [0.0; 3],
            bound_m: 2.0,
        }
    }
}

impl TemperatureField for Temperature {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.temperature_gradient(p)?.temperature_k)
    }
}
impl TemperatureGradientField for Temperature {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        let coordinates = [p.x, p.y, p.z];
        if coordinates
            .iter()
            .any(|value| !value.is_finite() || value.abs() > self.bound_m)
        {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let mut temperature_k = self.center_k;
        let gradient = core::array::from_fn::<_, 3, _>(|i| {
            temperature_k += self.gradient[i] * coordinates[i]
                + 0.5 * self.diagonal[i] * coordinates[i] * coordinates[i];
            self.gradient[i] + self.diagonal[i] * coordinates[i]
        });
        Ok(TemperatureGradientSample {
            temperature_k,
            gradient_k_per_m: Vec3::new(gradient[0], gradient[1], gradient[2]),
        })
    }
}
impl TemperatureSecondOrderField for Temperature {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: self.diagonal,
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 12,
        max_distance_m: 12.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                max_step_m: 0.04,
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
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (tolerance {tolerance:.3e})"
    );
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

/// Independent conserved-tangential-momentum solution for n = 1.5 + a*x,
/// then scalar bisection of the ORIGINAL unit sphere and vector Snell into air.
/// No production propagation, intersection or prescribed optical leg is used.
fn linear_reference(wavelength_nm: f64) -> (f64, Vec3, f64, f64, f64, f64) {
    let spectral = 0.0045 / (wavelength_nm * 1e-3).powi(2);
    let a = (5.0 / 12.0) * (0.012 + spectral) * 10.0;
    let ag = (5.0 / 12.0) * (0.012 + 3.0 * spectral) * 10.0;
    let position = |s: f64| {
        let n = (2.25 + (a * s).powi(2)).sqrt();
        ((n - 1.5) / a, 1.5 / a * (a * s / 1.5).asinh(), n)
    };
    let (mut lo, mut hi) = (0.0, 2.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let (x, z, _) = position(mid);
        if x * x + z * z > 1.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let s = 0.5 * (lo + hi);
    let (x, z, n) = position(s);
    let u = Vec3::new(a * s / n, 0.0, 1.5 / n);
    let ci = u.x * x + u.z * z;
    let ct = (1.0 - n * n * (1.0 - ci * ci)).sqrt();
    let outgoing = Vec3::new(
        n * u.x + (ct - n * ci) * x,
        0.0,
        n * u.z + (ct - n * ci) * z,
    );
    let outside = (2.0 - z) / outgoing.z;
    let phase_inside = 0.5 * (s * n + 1.5 * z);
    let group_inside = 1.5 * s + (ag / a) * (phase_inside - 1.5 * s);
    (
        x + outside * outgoing.x,
        outgoing,
        s + outside,
        phase_inside + outside,
        group_inside + outside,
        0.2 * s + 0.1 * outside,
    )
}

#[test]
fn g1_original_thermal_field_bends_spectral_camera_rays_to_original_emitter() {
    with_cx(|cx, _| {
        let owner = owner(cx, 0.15, 0.15, false);
        let medium = owner.medium_index(MEDIUM).unwrap();
        let law = thermal_law(0.012, 0.0045);
        let temperature = Temperature {
            gradient: [10.0, 0.0, 0.0],
            ..Temperature::uniform()
        };
        let field = MaterialThermoOpticField::try_new(
            owner.medium(medium).unwrap(),
            &law,
            &temperature,
            hash_domain("native-thermo-optic-test", b"linear-snapshot"),
        )
        .unwrap();
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        owner
            .with_curved_camera(cx, 100_000, calibration, |native, _| {
                let cold = native
                    .sample_ray_with_samples(cx, ray(), 550.0, medium, settings(), |_| {
                        [0.0, 0.0, 0.99]
                    })
                    .unwrap();
                assert_eq!(cold.radiance_w_m2_sr_nm(), 0.0);
            })
            .unwrap();
        for binding in [
            NativeThermoOpticBinding::from_gradient_field(MEDIUM, &field),
            NativeThermoOpticBinding::from_differentiable_field(MEDIUM, &field),
        ] {
            owner
                .with_thermo_optic_camera(cx, 100_000, calibration, &[binding], |native, camera| {
                    assert_eq!(camera.eye(), owner.render_scene().camera.eye);
                    for wavelength in [401.25, 550.0, 699.125] {
                        let sample = native
                            .sample_ray_with_samples(
                                cx,
                                ray(),
                                wavelength,
                                medium,
                                settings(),
                                |_| [0.0, 0.0, 0.99],
                            )
                            .unwrap();
                        let (x, direction, length, phase, group, depth) =
                            linear_reference(wavelength);
                        let path = sample.path();
                        assert_eq!(
                            sample.endpoint(),
                            RadianceEndpoint::Surface { object_id: 9 }
                        );
                        assert_eq!(path.last_contact().unwrap().object_id, 9);
                        assert_eq!(path.surface_events(), 2);
                        close(path.beam().ray().origin.x, x, 3e-7);
                        close(path.beam().ray().origin.z, 2.0, 3e-7);
                        close(path.beam().ray().dir.x, direction.x, 3e-7);
                        close(path.beam().ray().dir.z, direction.z, 3e-7);
                        close(path.traveled_m(), length, 3e-7);
                        close(path.optical_path_length_m(), phase, 4e-7);
                        close(path.group_optical_path_length_m().unwrap(), group, 4e-7);
                        close(path.optical_depth(), depth, 1e-7);
                        let expected = 6.0 * SPECTRUM.eval(wavelength) * (-depth).exp() * 2.25;
                        close(sample.radiance_w_m2_sr_nm(), expected, expected * 3e-7);
                    }
                })
                .unwrap();
        }
    });
}

#[test]
fn g1_original_thermal_hessian_focuses_bundle_through_actual_curved_interface() {
    with_cx(|cx, _| {
        let owner = owner(cx, 0.0, 4.0, true);
        let medium = owner.medium_index(MEDIUM).unwrap();
        let law = thermal_law(0.02, 0.0);
        let temperature = Temperature {
            diagonal: [-40.0, -40.0, 0.0],
            ..Temperature::uniform()
        };
        let field = MaterialThermoOpticField::try_new(
            owner.medium(medium).unwrap(),
            &law,
            &temperature,
            hash_domain("native-thermo-optic-test", b"quadratic-snapshot"),
        )
        .unwrap();
        let first = NativeThermoOpticBinding::from_gradient_field(MEDIUM, &field);
        let second = NativeThermoOpticBinding::from_differentiable_field(MEDIUM, &field);
        assert_eq!(first.field_identity(), second.field_identity());
        owner.with_thermo_optic_scene(cx, 100_000, &[first], |native| {
            assert!(matches!(
                native.launch_bundle(ray(), Vec3::new(1.0, 0.0, 0.0), 550.0, Stokes::unpolarized(1.0), medium),
                Err(CurvedBundleSceneError::MissingSecondOrderField { medium: index }) if index == medium
            ));
        }).unwrap();
        owner
            .with_thermo_optic_scene(cx, 100_000, &[second], |native| {
                let launch = native
                    .launch_bundle(
                        ray(),
                        Vec3::new(1.0, 0.0, 0.0),
                        550.0,
                        Stokes::unpolarized(1.0),
                        medium,
                    )
                    .unwrap();
                let path = native
                    .trace_bundle_transmitted(
                        cx,
                        launch,
                        medium,
                        CurvedBundleSettings {
                            trace: settings(),
                            position_derivative_tolerance_m: 1e-10,
                            direction_derivative_tolerance: 1e-11,
                            ..CurvedBundleSettings::default()
                        },
                    )
                    .unwrap();
                assert_eq!(path.termination(), PathTermination::Detected);
                assert_eq!(path.last_contact().unwrap().object_id, 9);
                let propagation = path.bundle().propagation();
                let footprint = propagation
                    .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
                    .unwrap();
                // On the unchanged central ray: x'' + gamma²*x = 0. Native sphere
                // curvature then contributes (1-n)*x/R at its Snell interface.
                let gamma = (2.0_f64 / 9.0).sqrt();
                let area = (0.5 * gamma.sin() / gamma + 1.5 * gamma.cos()).powi(2);
                assert!(area < 3.5); // The homogeneous original has area 4 m²/sr.
                close(footprint.area_per_solid_angle_m2(), area, 5e-7);
                close(propagation.propagation().optical_path_length_m, 2.5, 4e-7);
                close(
                    propagation
                        .propagation()
                        .group_optical_path_length_m
                        .unwrap(),
                    2.5,
                    4e-7,
                );
                close(propagation.propagation().optical_depth, 0.3, 1e-7);
                close(path.bundle().interface_log_power().exp(), 0.96, 1e-10);
            })
            .unwrap();
    });
}

#[test]
fn g0_original_reference_limit_and_frozen_physical_identity() {
    with_cx(|cx, _| {
        let owner = owner(cx, 0.0, 4.0, false);
        let medium = owner.medium_index(MEDIUM).unwrap();
        let reference = owner.medium(medium).unwrap();
        let law = thermal_law(0.012, 0.0045);
        let temperature = Temperature::uniform();
        let snapshot = hash_domain("native-thermo-optic-test", b"reference-snapshot");
        let field =
            MaterialThermoOpticField::try_new(reference, &law, &temperature, snapshot).unwrap();
        let binding = NativeThermoOpticBinding::from_differentiable_field(MEDIUM, &field);
        assert_eq!(binding.medium(), MEDIUM);
        assert_eq!(binding.source_snapshot_identity(), snapshot);
        let changed_snapshot = MaterialThermoOpticField::try_new(
            reference,
            &law,
            &temperature,
            hash_domain("native-thermo-optic-test", b"other-snapshot"),
        )
        .unwrap();
        assert_ne!(
            binding.field_identity(),
            NativeThermoOpticBinding::from_gradient_field(MEDIUM, &changed_snapshot)
                .field_identity()
        );
        let other_law = thermal_law(0.013, 0.0045);
        let changed_law =
            MaterialThermoOpticField::try_new(reference, &other_law, &temperature, snapshot)
                .unwrap();
        assert_ne!(
            binding.field_identity(),
            NativeThermoOpticBinding::from_gradient_field(MEDIUM, &changed_law).field_identity()
        );
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        let observe = |native: &CurvedCameraScene<'_, '_>, _: &PhysicalCamera| {
            native
                .sample_ray_with_samples(cx, ray(), 550.0, medium, settings(), |_| [0.0, 0.0, 0.99])
                .unwrap()
        };
        let original = owner
            .with_curved_camera(cx, 100_000, calibration, observe)
            .unwrap();
        let empty = owner
            .with_thermo_optic_camera(cx, 100_000, calibration, &[], observe)
            .unwrap();
        assert_eq!(empty, original);
        let thermal = owner
            .with_thermo_optic_camera(cx, 100_000, calibration, &[binding], observe)
            .unwrap();
        close(
            thermal.radiance_w_m2_sr_nm(),
            original.radiance_w_m2_sr_nm(),
            original.radiance_w_m2_sr_nm() * 3e-7,
        );
        close(
            thermal.path().optical_path_length_m(),
            original.path().optical_path_length_m(),
            4e-7,
        );
    });
}

#[test]
fn g3_original_field_binding_refuses_wrong_owner_domain_and_work() {
    with_cx(|cx, gate| {
        let owner = owner(cx, 0.0, 4.0, false);
        let medium = owner.medium_index(MEDIUM).unwrap();
        let law = thermal_law(0.01, 0.0);
        let temperature = Temperature::uniform();
        let snapshot = hash_domain("native-thermo-optic-test", b"refusal-snapshot");
        let wrong_absorption = glass(1.5, 0.3);
        let wrong =
            MaterialThermoOpticField::try_new(&wrong_absorption, &law, &temperature, snapshot)
                .unwrap();
        assert_eq!(
            owner.with_thermo_optic_scene(
                cx,
                100_000,
                &[NativeThermoOpticBinding::from_gradient_field(
                    MEDIUM, &wrong
                )],
                |_| ()
            ),
            Err(NativeOpticalSceneError::ReferenceGlassMismatch { medium: MEDIUM })
        );
        let field = MaterialThermoOpticField::try_new(
            owner.medium(medium).unwrap(),
            &law,
            &temperature,
            snapshot,
        )
        .unwrap();
        let binding = NativeThermoOpticBinding::from_gradient_field(MEDIUM, &field);
        assert_eq!(
            owner.with_thermo_optic_scene(cx, 100_000, &[binding, binding], |_| ()),
            Err(NativeOpticalSceneError::DuplicateFieldBinding { medium: MEDIUM })
        );
        assert_eq!(
            owner.with_thermo_optic_scene(cx, 100_000, &[binding; 3], |_| ()),
            Err(NativeOpticalSceneError::TooManyFieldBindings {
                bindings: 3,
                media: 2
            })
        );
        let invalid = NativeOpticalMedium::InteriorOf(1);
        assert_eq!(
            owner.with_thermo_optic_scene(
                cx,
                100_000,
                &[NativeThermoOpticBinding::from_gradient_field(
                    invalid, &field
                )],
                |_| ()
            ),
            Err(NativeOpticalSceneError::InvalidMediumReference { medium: invalid })
        );
        assert!(
            owner
                .with_thermo_optic_scene(cx, 0, &[binding], |_| ())
                .is_err()
        );
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        for temperature in [
            Temperature {
                bound_m: 0.2,
                ..Temperature::uniform()
            },
            Temperature {
                center_k: 451.0,
                ..Temperature::uniform()
            },
        ] {
            let unsupported = MaterialThermoOpticField::try_new(
                owner.medium(medium).unwrap(),
                &law,
                &temperature,
                snapshot,
            )
            .unwrap();
            owner
                .with_thermo_optic_camera(
                    cx,
                    100_000,
                    calibration,
                    &[NativeThermoOpticBinding::from_gradient_field(
                        MEDIUM,
                        &unsupported,
                    )],
                    |native, _| {
                        assert!(
                            native
                                .sample_ray_with_samples(
                                    cx,
                                    ray(),
                                    550.0,
                                    medium,
                                    settings(),
                                    |_| [0.0, 0.0, 0.99]
                                )
                                .is_err()
                        );
                    },
                )
                .unwrap();
        }
        gate.request();
        assert_eq!(
            owner.with_thermo_optic_scene(cx, 100_000, &[binding], |_| ()),
            Err(NativeOpticalSceneError::Cancelled)
        );
    });
}
