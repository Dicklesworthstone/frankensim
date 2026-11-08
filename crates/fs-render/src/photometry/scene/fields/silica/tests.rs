use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexOptions};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::photometry::scene::{
    NativeGeometryBinding, NativeOpticalRole, NativeRadianceCalibration,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::scene::curved::bundle::CurvedBundleSettings;
use crate::polarization::dielectric::scene::curved::camera::RadianceEndpoint;
use crate::spectral::{LiftedSpectrum, SpectralReflectance};
use crate::tracer::{Camera, Primitive, Shape};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_material::silica_refraction::Corning7980Silica;

const AMBIENT: NativeOpticalMedium = NativeOpticalMedium::Ambient;
const SPECTRUM: LiftedSpectrum = LiftedSpectrum {
    c: [3.0, -2.0, -0.25],
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 161,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

struct Temperature {
    base: f64,
    slope: f64,
    curvature: f64,
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
        if ![p.x, p.y, p.z].into_iter().all(f64::is_finite) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        Ok(TemperatureGradientSample {
            temperature_k: self.base
                + self.slope * p.x
                + 0.5 * self.curvature * (p.x * p.x + p.y * p.y),
            gradient_k_per_m: Vec3::new(
                self.slope + self.curvature * p.x,
                self.curvature * p.y,
                0.0,
            ),
        })
    }
    fn is_spatially_uniform(&self) -> bool {
        self.slope == 0.0 && self.curvature == 0.0
    }
}
impl TemperatureSecondOrderField for Temperature {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [self.curvature, self.curvature, 0.0],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 4,
        max_distance_m: 2.0,
        max_total_steps: 16_384,
        max_query_work: 100_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11,
                optical_depth_tolerance: 1e-12,
                relative_tolerance: 1e-11,
                max_step_m: 0.025,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-9,
            ..GradientEventOptions::default()
        },
    }
}
fn owner(
    cx: &Cx<'_>,
    reference: DielectricGlass,
    center: f64,
    width: f64,
    detector: bool,
) -> NativeOpticalScene {
    let scene = Scene {
        primitives: vec![Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [center - width, -width, 1.0],
                    [center + width, -width, 1.0],
                    [center, 2.0 * width, 1.0],
                ],
                vec![[0, 1, 2]],
            )),
            material: Material::SpectralDiffuse {
                reflectance: SpectralReflectance::constant(0.0).unwrap(),
            },
            emission: (!detector).then_some((SPECTRUM, 2.0)),
        }],
        lights: vec![],
        environment: None,
        camera: Camera {
            eye: ray().origin,
            forward: ray().dir,
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.05,
        },
    };
    let bindings = [NativeOpticalBinding {
        geometry: NativeGeometryBinding::Unplaced {
            object_id: 9,
            geometry_identity: ContentHash([162; 32]),
        },
        role: if detector {
            NativeOpticalRole::Detector {
                reference_axis: Vec3::new(1.0, 0.0, 0.0),
            }
        } else {
            NativeOpticalRole::Material { outside: AMBIENT }
        },
    }];
    NativeOpticalScene::try_from_scene(cx, scene, reference, &bindings, 1).unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_silica_scene_uses_actual_temperature_and_analytic_bundle_curvature() {
    with_cx(|cx, _| {
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let reference = corning7980_glass(180.0, absorption).unwrap();
        let owner = owner(cx, reference, 0.0, 1.0, true);
        for curvature in [0.0, -40.0] {
            let temperature = Temperature {
                base: 220.0,
                slope: 0.0,
                curvature,
            };
            let field = MaterialSilicaRefractiveField::try_new(
                &temperature,
                &absorption,
                ContentHash([if curvature == 0.0 { 163 } else { 167 }; 32]),
            )
            .unwrap();
            let binding = NativeMaterialFieldBinding::from_differentiable_silica_field(
                AMBIENT, &reference, 180.0, &field,
            )
            .unwrap();
            for wavelength in [550.0, 1550.0] {
                let source = Corning7980Silica::at_temperature(220.0)
                    .unwrap()
                    .evaluate(wavelength)
                    .unwrap();
                let trace = |native: &CurvedPolarizedScene<'_>| {
                    let launch = native
                        .launch_bundle(
                            ray(),
                            Vec3::new(1.0, 0.0, 0.0),
                            wavelength,
                            Stokes::unpolarized(1.0),
                            0,
                        )
                        .unwrap();
                    native
                        .trace_bundle_transmitted(
                            cx,
                            launch,
                            0,
                            CurvedBundleSettings {
                                trace: settings(),
                                position_derivative_tolerance_m: 1e-10,
                                direction_derivative_tolerance: 1e-11,
                                ..CurvedBundleSettings::default()
                            },
                        )
                        .unwrap()
                };
                let path = owner
                    .with_thermo_optic_scene(cx, 100_000, &[binding], trace)
                    .unwrap();
                let prepared = field.prepare(wavelength).unwrap();
                assert!(core::ptr::eq(prepared.source_field(), &field));
                let selected =
                    NativeMaterialFieldBinding::from_prepared_differentiable_silica_field(
                        AMBIENT, &reference, 180.0, &prepared,
                    )
                    .unwrap();
                assert_eq!(binding.field_identity(), selected.field_identity());
                let repeat = owner
                    .with_thermo_optic_scene(cx, 100_000, &[selected], trace)
                    .unwrap();
                assert_eq!(path, repeat);
                assert_eq!(path.termination(), PathTermination::Detected);
                let bundle = path.bundle().propagation();
                let propagation = bundle.propagation();
                close(propagation.optical_path_length_m, source.phase_index, 5e-8);
                close(
                    propagation.group_optical_path_length_m.unwrap(),
                    source.group_index,
                    5e-8,
                );
                close(propagation.optical_depth, 0.2, 1e-8);
                // Along the axis, the original T is constant and grad(T)=0.
                // The exact transverse variational equation is x''=-gamma²*x.
                let derivative = if curvature == 0.0 {
                    1.0
                } else {
                    let gamma = (-source.derivative_temperature_per_k * curvature
                        / source.phase_index)
                        .sqrt();
                    assert!(gamma > 0.0);
                    gamma.sin() / gamma
                };
                let projection = bundle
                    .project_at_endpoint(ray().dir, ProjectionOptions::default())
                    .unwrap();
                close(
                    projection.area_per_solid_angle_m2(),
                    derivative * derivative,
                    1e-7,
                );
                if curvature != 0.0 {
                    assert!(projection.area_per_solid_angle_m2() < 1.0 - 1e-5);
                }
            }
        }
    });
}

#[test]
fn g1_original_silica_gradient_changes_the_native_camera_observation() {
    with_cx(|cx, _| {
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let reference = corning7980_glass(180.0, absorption).unwrap();
        let temperature = Temperature {
            base: 220.0,
            slope: 40.0,
            curvature: 0.0,
        };
        let field = MaterialSilicaRefractiveField::try_new(
            &temperature,
            &absorption,
            ContentHash([164; 32]),
        )
        .unwrap();
        let initial = Corning7980Silica::at_temperature(220.0)
            .unwrap()
            .evaluate(550.0)
            .unwrap();
        // A weak-gradient paraxial estimate only places a finite target. Actual
        // camera transport still uses the original nonlinear source and event walk.
        let center =
            0.5 * initial.derivative_temperature_per_k * temperature.slope / initial.phase_index;
        assert!(center > 1e-5);
        let owner = owner(cx, reference, center, 0.1 * center, false);
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        owner
            .with_curved_camera(cx, 100_000, calibration, |camera, _| {
                let baseline = camera
                    .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |_| [0.0, 0.0, 0.99])
                    .unwrap();
                assert_eq!(baseline.radiance_w_m2_sr_nm(), 0.0);
            })
            .unwrap();
        let prepared = field.prepare(550.0).unwrap();
        let binding = NativeMaterialFieldBinding::from_prepared_silica_field(
            AMBIENT, &reference, 180.0, &prepared,
        )
        .unwrap();
        owner
            .with_thermo_optic_camera(cx, 100_000, calibration, &[binding], |camera, physical| {
                assert_eq!(physical.eye(), owner.render_scene().camera.eye);
                let sample = camera
                    .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |_| [0.0, 0.0, 0.99])
                    .unwrap();
                assert_eq!(
                    sample.endpoint(),
                    RadianceEndpoint::Surface { object_id: 9 }
                );
                let path = sample.path();
                assert!((path.beam().ray().origin.x - center).abs() < center * 0.05);
                let endpoint = prepared.optical_sample(path.beam().ray().origin).unwrap();
                let expected = 6.0
                    * SPECTRUM.eval(550.0)
                    * (-0.2 * path.traveled_m()).exp()
                    * (initial.phase_index / endpoint.phase_index()).powi(2);
                close(sample.radiance_w_m2_sr_nm(), expected, 1e-7 * expected);
                assert!(
                    camera
                        .sample_ray_with_samples(cx, ray(), 551.0, 0, settings(), |_| [
                            0.0, 0.0, 0.99
                        ])
                        .is_err()
                );
            })
            .unwrap();
    });
}

#[test]
fn g3_silica_binding_retains_original_reference_absorption_snapshot_and_cancellation() {
    with_cx(|cx, gate| {
        let temperature = Temperature {
            base: 220.0,
            slope: 0.0,
            curvature: 0.0,
        };
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let reference = corning7980_glass(180.0, absorption).unwrap();
        let field = MaterialSilicaRefractiveField::try_new(
            &temperature,
            &absorption,
            ContentHash([165; 32]),
        )
        .unwrap();
        let binding =
            NativeMaterialFieldBinding::from_silica_field(AMBIENT, &reference, 180.0, &field)
                .unwrap();
        let snapshot = MaterialSilicaRefractiveField::try_new(
            &temperature,
            &absorption,
            ContentHash([166; 32]),
        )
        .unwrap();
        let changed =
            NativeMaterialFieldBinding::from_silica_field(AMBIENT, &reference, 180.0, &snapshot)
                .unwrap();
        assert_ne!(binding.field_identity(), changed.field_identity());
        for invalid_temperature in [181.0, 300.1, f64::NAN] {
            assert!(matches!(
                NativeMaterialFieldBinding::from_silica_field(
                    AMBIENT,
                    &reference,
                    invalid_temperature,
                    &field
                ),
                Err(NativeOpticalSceneError::SilicaReferenceRequired { .. })
            ));
        }
        let lookalike = DielectricGlass::new(
            CauchyIor::try_constant(reference.ior().eval(550.0).unwrap()).unwrap(),
            absorption,
            GlassProvenance::CharmsCorning7980V1,
        );
        assert!(matches!(
            NativeMaterialFieldBinding::from_silica_field(AMBIENT, &lookalike, 180.0, &field),
            Err(NativeOpticalSceneError::SilicaReferenceRequired { .. })
        ));
        let different_absorption = corning7980_glass(180.0, BeerLambertAbsorption::CLEAR).unwrap();
        assert!(matches!(
            NativeMaterialFieldBinding::from_silica_field(
                AMBIENT,
                &different_absorption,
                180.0,
                &field
            ),
            Err(NativeOpticalSceneError::SilicaAbsorptionMismatch { .. })
        ));
        let owner = owner(cx, reference, 0.0, 1.0, true);
        gate.request();
        let called = std::cell::Cell::new(false);
        let result = owner.with_thermo_optic_scene(cx, 100_000, &[binding], |_| called.set(true));
        assert!(matches!(result, Err(NativeOpticalSceneError::Cancelled)));
        assert!(!called.get());
    });
}
