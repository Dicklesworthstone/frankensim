use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, GlassProvenance, IapwsWaterIor, SellmeierIor,
};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexOptions};
use crate::material_state::gas_scattering::isobaric::IsobaricGasField;
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::material_state::water_refraction::isobaric::IsobaricWaterField;
use crate::material_state::water_refraction::pressure_temperature::{
    PressureTemperatureGradientField, PressureTemperatureGradientSample, PressureTemperatureJet,
    PressureTemperatureSecondOrderField, PressureTemperatureWaterField,
};
use crate::photometry::scene::{
    NativeGeometryBinding, NativeOpticalRole, NativeRadianceCalibration,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::scene::curved::bundle::{
    CurvedBundlePathResult, CurvedBundleSceneError, CurvedBundleSettings,
};
use crate::polarization::dielectric::scene::curved::camera::RadianceEndpoint;
use crate::spectral::{LiftedSpectrum, SpectralReflectance};
use crate::tracer::{Camera, Primitive, Shape};
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3, fixtures::SphereChart};
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::{MolecularRayleigh, MolecularRayleighMixture};
use fs_material::water_refraction::IapwsWaterRefraction;
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::{If97Region, If97State};

const WATER: NativeOpticalMedium = NativeOpticalMedium::InteriorOf(0);
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
                seed: 97,
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

fn state(pressure: f64, temperature: f64) -> If97State {
    If97State::at_pressure_temperature(
        Pressure::new(pressure),
        Temperature::new(temperature),
        If97Region::Region1,
    )
    .unwrap()
}

fn water_reference(state: If97State, absorption: BeerLambertAbsorption) -> DielectricGlass {
    DielectricGlass::new(
        IapwsWaterIor::try_new(state.temperature().value(), state.density().value()).unwrap(),
        absorption,
        GlassProvenance::Custom,
    )
}

fn ambient() -> DielectricGlass {
    DielectricGlass::new(
        SellmeierIor::try_new_in_band([0.0; 3], [0.0; 3], (200.0, 1100.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.1).unwrap(),
        GlassProvenance::Custom,
    )
}

fn scene(
    water: DielectricGlass,
    center_x: f64,
    width: f64,
    detector: bool,
) -> (Scene, [NativeOpticalBinding; 2]) {
    let original = Scene {
        primitives: vec![
            Primitive {
                shape: Shape::Chart(Box::new(SphereChart {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                })),
                material: Material::Dielectric {
                    glass: water,
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
            geometry_identity: hash_domain("native-water-scene-test", &id.to_le_bytes()),
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
    (original, bindings)
}

fn admit(cx: &Cx<'_>, original: (Scene, [NativeOpticalBinding; 2])) -> NativeOpticalScene {
    NativeOpticalScene::try_from_scene(cx, original.0, ambient(), &original.1, 2).unwrap()
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 8.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11,
                optical_depth_tolerance: 1e-12,
                relative_tolerance: 1e-11,
                max_step_m: 0.04,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}
fn bundle(
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    wavelength: f64,
    medium: usize,
) -> CurvedBundlePathResult {
    let launch = scene
        .launch_bundle(
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            wavelength,
            Stokes::unpolarized(1.0),
            medium,
        )
        .unwrap();
    scene
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
        .unwrap()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (tolerance {tolerance:.3e})"
    );
}

// Smooth manufactured snapshots. No hydrostatic/momentum equilibrium or moving
// fluid geometry is asserted; the original single-phase EOS admits each state.
struct Joint {
    temperature_slope: f64,
    temperature_curvature: f64,
    pressure: f64,
    pressure_curvature: f64,
}
impl Joint {
    fn uniform(pressure: f64) -> Self {
        Self {
            temperature_slope: 0.0,
            temperature_curvature: 0.0,
            pressure,
            pressure_curvature: 0.0,
        }
    }
}
impl TemperatureField for Joint {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.temperature_gradient(p)?.temperature_k)
    }
}
impl TemperatureGradientField for Joint {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        if ![p.x, p.y, p.z].into_iter().all(f64::is_finite)
            || p.x.abs() > 0.25
            || p.y.abs() > 0.25
            || p.z.abs() > 3.0
        {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        Ok(TemperatureGradientSample {
            temperature_k: 300.0
                + self.temperature_slope * p.x
                + 0.5 * self.temperature_curvature * (p.x * p.x + p.y * p.y),
            gradient_k_per_m: Vec3::new(
                self.temperature_slope + self.temperature_curvature * p.x,
                self.temperature_curvature * p.y,
                0.0,
            ),
        })
    }
}
impl TemperatureSecondOrderField for Joint {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [self.temperature_curvature, self.temperature_curvature, 0.0],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}
impl PressureTemperatureGradientField for Joint {
    fn pressure_temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        Ok(PressureTemperatureGradientSample {
            temperature: self.temperature_gradient(p)?,
            pressure_pa: self.pressure + 0.5 * self.pressure_curvature * (p.x * p.x + p.y * p.y),
            pressure_gradient_pa_per_m: Vec3::new(
                self.pressure_curvature * p.x,
                self.pressure_curvature * p.y,
                0.0,
            ),
        })
    }
}
impl PressureTemperatureSecondOrderField for Joint {
    fn pressure_temperature_jet(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureJet, GradientIndexError> {
        let sample = self.pressure_temperature_gradient(p)?;
        Ok(PressureTemperatureJet {
            temperature: self.temperature_jet(p)?,
            pressure_pa: sample.pressure_pa,
            pressure_gradient_pa_per_m: sample.pressure_gradient_pa_per_m,
            pressure_diagonal_pa_per_m2: [self.pressure_curvature, self.pressure_curvature, 0.0],
            pressure_off_diagonal_pa_per_m2: [0.0; 3],
        })
    }
}

#[test]
fn g1_original_joint_water_state_focuses_through_the_original_sphere_and_detector() {
    with_cx(|cx, _| {
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let original_state = state(40e6, 300.0);
        let owner = admit(
            cx,
            scene(water_reference(original_state, absorption), 0.0, 4.0, true),
        );
        let medium = owner.medium_index(WATER).unwrap();
        let reference = owner.medium(medium).unwrap();
        let source = Joint {
            temperature_curvature: 60.0,
            pressure_curvature: -300e6,
            ..Joint::uniform(40e6)
        };
        let water = PressureTemperatureWaterField::try_new(
            &source,
            If97Region::Region1,
            ContentHash([131; 32]),
        )
        .unwrap();
        let optics = water.bind_refraction(&absorption).unwrap();
        let binding =
            NativeMaterialFieldBinding::from_differentiable_water_field(WATER, reference, &optics)
                .unwrap();
        let old_name: NativeThermoOpticBinding<'_> = binding;
        assert_eq!(old_name.field_identity(), binding.field_identity());
        for wavelength in [401.25, 589.0, 1013.98] {
            let material = IapwsWaterRefraction::prepare(wavelength)
                .unwrap()
                .sample_if97_state(&original_state)
                .unwrap();
            let rho = original_state.density().value();
            let n_t = material.temperature_derivative_per_k
                - material.density_derivative_m3_per_kg
                    * rho
                    * original_state.thermal_expansion().value();
            let n_p = material.density_derivative_m3_per_kg
                * rho
                * original_state.isothermal_compressibility().value();
            let curvature = n_t * 60.0 - n_p * 300e6;
            let gamma = (-curvature / material.phase_index).sqrt();
            // Exact axial variational solution in water, then the ORIGINAL
            // unit sphere's Snell normal derivative and one metre outside.
            let derivative = (2.0 - material.phase_index) * gamma.sin() / gamma
                + material.phase_index * gamma.cos();
            assert!(derivative * derivative < 3.95);
            let path = owner
                .with_thermo_optic_scene(cx, 100_000, &[binding], |native| {
                    bundle(cx, native, wavelength, medium)
                })
                .unwrap();
            let prepared = optics.prepare(wavelength).unwrap();
            assert!(core::ptr::eq(prepared.source_field(), &optics));
            let selected = NativeMaterialFieldBinding::from_prepared_differentiable_water_field(
                WATER, reference, &prepared,
            )
            .unwrap();
            assert_eq!(binding.field_identity(), selected.field_identity());
            let repeat = owner
                .with_thermo_optic_scene(cx, 100_000, &[selected], |native| {
                    bundle(cx, native, wavelength, medium)
                })
                .unwrap();
            assert_eq!(path, repeat);
            assert_eq!(path.termination(), PathTermination::Detected);
            assert_eq!(path.last_contact().unwrap().object_id, 9);
            let propagation = path.bundle().propagation();
            let projection = propagation
                .project_at_endpoint(ray().dir, ProjectionOptions::default())
                .unwrap();
            close(
                projection.area_per_solid_angle_m2(),
                derivative * derivative,
                1e-6,
            );
            close(
                propagation.propagation().optical_path_length_m,
                material.phase_index + 1.0,
                4e-7,
            );
            close(
                propagation
                    .propagation()
                    .group_optical_path_length_m
                    .unwrap(),
                material.group_index + 1.0,
                4e-7,
            );
            close(propagation.propagation().optical_depth, 0.3, 1e-7);
            let transmission = 4.0 * material.phase_index / (material.phase_index + 1.0).powi(2);
            close(
                path.bundle().interface_log_power().exp(),
                transmission,
                1e-10,
            );
        }
        let baseline = owner
            .with_curved_scene(cx, 100_000, |native| bundle(cx, native, 589.0, medium))
            .unwrap();
        close(
            baseline
                .bundle()
                .propagation()
                .project_at_endpoint(ray().dir, ProjectionOptions::default())
                .unwrap()
                .area_per_solid_angle_m2(),
            4.0,
            1e-6,
        );
        let first =
            NativeMaterialFieldBinding::from_water_field(WATER, reference, &optics).unwrap();
        owner.with_thermo_optic_scene(cx, 100_000, &[first], |native| {
            assert!(matches!(native.launch_bundle(ray(), Vec3::new(1.0, 0.0, 0.0), 589.0, Stokes::unpolarized(1.0), medium),
                Err(CurvedBundleSceneError::MissingSecondOrderField { medium: index }) if index == medium));
        }).unwrap();
    });
}

struct CameraReference {
    x: f64,
    direction: Vec3,
    length: f64,
    phase: f64,
    group: f64,
    depth: f64,
    launch_index: f64,
}

// Independent conserved n*u_z for the original IF97/IAPWS isobaric state.
// x=-q² removes the axial turning-point singularity. Fixed one-dimensional
// Simpson quadrature supplies z, physical/phase/group length; bisection finds
// the actual original sphere. No native ODE, event walker or Hessian is used.
fn camera_reference(wavelength: f64) -> CameraReference {
    let source = IapwsWaterRefraction::prepare(wavelength).unwrap();
    let initial = state(3e6, 300.0);
    let first = source.sample_if97_state(&initial).unwrap();
    let n0 = first.phase_index;
    let nx = 100.0
        * (first.temperature_derivative_per_k
            - first.density_derivative_m3_per_kg
                * initial.density().value()
                * initial.thermal_expansion().value());
    assert!(nx < 0.0);
    let sample = |q: f64| {
        source
            .sample_if97_state(&state(3e6, 300.0 - 100.0 * q * q))
            .unwrap()
    };
    let integrate = |q: f64| {
        let mut integral = [0.0; 4];
        let step = q / 96.0;
        for i in 0..=96 {
            let local = sample(step * i as f64);
            let n = local.phase_index;
            let factor = if i == 0 {
                2.0 / (-2.0 * n0 * nx).sqrt()
            } else {
                2.0 * step * i as f64 / ((n - n0) * (n + n0)).sqrt()
            };
            let weight = if i == 0 || i == 96 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for (sum, value) in integral
                .iter_mut()
                .zip([n0, n, n * n, n * local.group_index])
            {
                *sum += weight * factor * value;
            }
        }
        integral.map(|value| value * step / 3.0)
    };
    let (mut lo, mut hi) = (0.0, 0.12);
    assert!(integrate(hi)[0] > 1.0);
    for _ in 0..48 {
        let q = 0.5 * (lo + hi);
        if q.powi(4) + integrate(q)[0].powi(2) > 1.0 {
            hi = q;
        } else {
            lo = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let x = -q * q;
    let [z, length, phase, group] = integrate(q);
    let n = sample(q).phase_index;
    let incoming = Vec3::new(-(1.0 - (n0 / n).powi(2)).sqrt(), 0.0, n0 / n);
    let ci = incoming.x * x + incoming.z * z;
    let ct = (1.0 - n * n * (1.0 - ci * ci)).sqrt();
    let direction = Vec3::new(
        n * incoming.x + (ct - n * ci) * x,
        0.0,
        n * incoming.z + (ct - n * ci) * z,
    );
    let outside = (2.0 - z) / direction.z;
    CameraReference {
        x: x + outside * direction.x,
        direction,
        length: length + outside,
        phase: phase + outside,
        group: group + outside,
        depth: 0.2 * length + 0.1 * outside,
        launch_index: n0,
    }
}

#[test]
fn g1_original_isobaric_water_bends_camera_rays_to_the_original_spectral_source() {
    with_cx(|cx, _| {
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let center = camera_reference(589.0).x;
        assert!(center < -0.003);
        let owner = admit(
            cx,
            scene(
                water_reference(state(3e6, 300.0), absorption),
                center,
                0.7 * center.abs(),
                false,
            ),
        );
        let medium = owner.medium_index(WATER).unwrap();
        let source = Joint {
            temperature_slope: 100.0,
            ..Joint::uniform(3e6)
        };
        let water = IsobaricWaterField::try_new(
            &source,
            Pressure::new(3e6),
            If97Region::Region1,
            ContentHash([132; 32]),
        )
        .unwrap();
        let optics = water.bind_refraction(&absorption).unwrap();
        let binding = NativeMaterialFieldBinding::from_water_field(
            WATER,
            owner.medium(medium).unwrap(),
            &optics,
        )
        .unwrap();
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        owner
            .with_curved_camera(cx, 100_000, calibration, |camera, _| {
                let cold = camera
                    .sample_ray_with_samples(cx, ray(), 589.0, medium, settings(), |_| {
                        [0.0, 0.0, 0.99]
                    })
                    .unwrap();
                assert_eq!(cold.radiance_w_m2_sr_nm(), 0.0);
            })
            .unwrap();
        owner
            .with_thermo_optic_camera(cx, 100_000, calibration, &[binding], |camera, physical| {
                assert_eq!(physical.eye(), owner.render_scene().camera.eye);
                for wavelength in [401.25, 589.0, 699.125] {
                    let expected = camera_reference(wavelength);
                    let sample = camera
                        .sample_ray_with_samples(cx, ray(), wavelength, medium, settings(), |_| {
                            [0.0, 0.0, 0.99]
                        })
                        .unwrap();
                    let path = sample.path();
                    assert_eq!(
                        sample.endpoint(),
                        RadianceEndpoint::Surface { object_id: 9 }
                    );
                    assert_eq!(path.surface_events(), 2);
                    close(path.beam().ray().origin.x, expected.x, 2e-6);
                    close(path.beam().ray().dir.x, expected.direction.x, 2e-6);
                    close(path.traveled_m(), expected.length, 2e-6);
                    close(path.optical_path_length_m(), expected.phase, 3e-6);
                    close(
                        path.group_optical_path_length_m().unwrap(),
                        expected.group,
                        3e-6,
                    );
                    close(path.optical_depth(), expected.depth, 5e-7);
                    let radiance = 6.0
                        * SPECTRUM.eval(wavelength)
                        * (-expected.depth).exp()
                        * expected.launch_index.powi(2);
                    close(sample.radiance_w_m2_sr_nm(), radiance, radiance * 8e-7);
                }
            })
            .unwrap();
    });
}

#[test]
fn g1_original_water_and_molecular_ambient_compose_in_the_same_camera() {
    with_cx(|cx, _| {
        let gas_source = MolecularRayleighMixture::try_new(
            &[GasSpec::dry_air_ussa1976()],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
            &[
                MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([133; 32]))
                    .unwrap(),
            ],
            None,
        )
        .unwrap();
        let source = Joint::uniform(3e6);
        let gas = IsobaricGasField::try_new(
            &source,
            &gas_source,
            Pressure::new(1e5),
            Temperature::new(300.0),
            ContentHash([134; 32]),
        )
        .unwrap();
        let gas_optics = gas.bind_refraction().unwrap();
        let gas_binding = NativeGasAmbient::from_density_field(&gas_optics);
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let original_state = state(3e6, 300.0);
        let (original, bindings) =
            scene(water_reference(original_state, absorption), 0.0, 4.0, false);
        let owner =
            NativeGasOpticalScene::try_from_scene(cx, original, gas_binding, &bindings, 2).unwrap();
        let medium = owner.medium_index(WATER).unwrap();
        let water = PressureTemperatureWaterField::try_new(
            &source,
            If97Region::Region1,
            ContentHash([135; 32]),
        )
        .unwrap();
        let optics = water.bind_refraction(&absorption).unwrap();
        let binding = NativeMaterialFieldBinding::from_water_field(
            WATER,
            owner.medium(medium).unwrap(),
            &optics,
        )
        .unwrap();
        owner
            .with_curved_camera(
                cx,
                100_000,
                NativeRadianceCalibration::try_new(3.0, 0.0).unwrap(),
                gas_binding,
                &[binding],
                |camera, _| {
                    let sample = camera
                        .sample_ray_with_samples(cx, ray(), 589.0, medium, settings(), |_| {
                            [0.0, 0.0, 0.99]
                        })
                        .unwrap();
                    let water_index = IapwsWaterRefraction::prepare(589.0)
                        .unwrap()
                        .sample_if97_state(&original_state)
                        .unwrap();
                    let gas_index = gas_optics
                        .sample(ray().origin, 589.0)
                        .unwrap()
                        .phase_index();
                    assert_eq!(
                        sample.endpoint(),
                        RadianceEndpoint::Surface { object_id: 9 }
                    );
                    assert_eq!(
                        sample.path().medium_snapshot_identity(),
                        Some(gas_optics.content_identity())
                    );
                    close(sample.path().optical_depth(), 0.2, 1e-7);
                    close(
                        sample.path().optical_path_length_m(),
                        water_index.phase_index + gas_index,
                        4e-7,
                    );
                    close(
                        sample.path().group_optical_path_length_m().unwrap(),
                        water_index.group_index + gas_index,
                        4e-7,
                    );
                    let expected = 6.0
                        * SPECTRUM.eval(589.0)
                        * (-0.2_f64).exp()
                        * (water_index.phase_index / gas_index).powi(2);
                    close(sample.radiance_w_m2_sr_nm(), expected, expected * 3e-7);
                },
            )
            .unwrap();
    });
}

#[test]
fn g3_water_binding_preserves_reference_absorption_preparation_and_physical_domains() {
    with_cx(|cx, gate| {
        let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
        let reference = water_reference(state(3e6, 300.0), absorption);
        let owner = admit(cx, scene(reference, 0.0, 4.0, true));
        let medium = owner.medium_index(WATER).unwrap();
        let source = Joint::uniform(3e6);
        let water = PressureTemperatureWaterField::try_new(
            &source,
            If97Region::Region1,
            ContentHash([136; 32]),
        )
        .unwrap();
        let optics = water.bind_refraction(&absorption).unwrap();
        let binding =
            NativeMaterialFieldBinding::from_differentiable_water_field(WATER, &reference, &optics)
                .unwrap();
        let fit = DielectricGlass::new(
            CauchyIor::try_constant(reference.ior().eval(589.0).unwrap()).unwrap(),
            absorption,
            GlassProvenance::Custom,
        );
        assert!(matches!(
            NativeMaterialFieldBinding::from_water_field(WATER, &fit, &optics),
            Err(NativeOpticalSceneError::WaterReferenceRequired { .. })
        ));
        let changed_absorption = BeerLambertAbsorption::try_constant(0.3).unwrap();
        let other_optics = water.bind_refraction(&changed_absorption).unwrap();
        assert_eq!(optics.snapshot_identity(), other_optics.snapshot_identity());
        let other_prepared = other_optics.prepare(589.0).unwrap();
        assert!(matches!(
            NativeMaterialFieldBinding::from_prepared_water_field(
                WATER,
                &reference,
                &other_prepared
            ),
            Err(NativeOpticalSceneError::WaterAbsorptionMismatch { .. })
        ));
        let changed_reference = water_reference(state(3e6, 310.0), absorption);
        let mismatch =
            NativeMaterialFieldBinding::from_water_field(WATER, &changed_reference, &optics)
                .unwrap();
        assert_eq!(
            owner.with_thermo_optic_scene(cx, 100_000, &[mismatch], |_| ()),
            Err(NativeOpticalSceneError::ReferenceGlassMismatch { medium: WATER })
        );
        assert_eq!(
            owner.with_thermo_optic_scene(cx, 100_000, &[binding, binding], |_| ()),
            Err(NativeOpticalSceneError::DuplicateFieldBinding { medium: WATER })
        );
        let changed = PressureTemperatureWaterField::try_new(
            &source,
            If97Region::Region1,
            ContentHash([137; 32]),
        )
        .unwrap();
        let changed_optics = changed.bind_refraction(&absorption).unwrap();
        let fresh =
            NativeMaterialFieldBinding::from_water_field(WATER, &reference, &changed_optics)
                .unwrap();
        assert_ne!(fresh.field_identity(), binding.field_identity());
        let prepared = optics.prepare(589.0).unwrap();
        let selected = NativeMaterialFieldBinding::from_prepared_differentiable_water_field(
            WATER, &reference, &prepared,
        )
        .unwrap();
        for (field, wavelength, origin) in [
            (selected, 600.0, ray().origin),
            (binding, 1100.0_f64.next_up(), ray().origin),
            (binding, 589.0, Point3::new(0.3, 0.0, 0.0)),
        ] {
            owner
                .with_thermo_optic_scene(cx, 100_000, &[field], |native| {
                    assert!(
                        native
                            .launch_bundle(
                                Ray { origin, ..ray() },
                                Vec3::new(1.0, 0.0, 0.0),
                                wavelength,
                                Stokes::unpolarized(1.0),
                                medium
                            )
                            .is_err()
                    );
                })
                .unwrap();
        }
        let invalid = Joint::uniform(-1.0);
        let invalid_water = PressureTemperatureWaterField::try_new(
            &invalid,
            If97Region::Region1,
            ContentHash([138; 32]),
        )
        .unwrap();
        let invalid_optics = invalid_water.bind_refraction(&absorption).unwrap();
        let invalid_binding = NativeMaterialFieldBinding::from_differentiable_water_field(
            WATER,
            &reference,
            &invalid_optics,
        )
        .unwrap();
        owner
            .with_thermo_optic_scene(cx, 100_000, &[invalid_binding], |native| {
                assert!(
                    native
                        .launch_bundle(
                            ray(),
                            Vec3::new(1.0, 0.0, 0.0),
                            589.0,
                            Stokes::unpolarized(1.0),
                            medium
                        )
                        .is_err()
                );
            })
            .unwrap();
        assert!(
            owner
                .with_thermo_optic_scene(cx, 0, &[binding], |_| ())
                .is_err()
        );
        gate.request();
        assert_eq!(
            owner.with_thermo_optic_scene(cx, 100_000, &[binding], |_| ()),
            Err(NativeOpticalSceneError::Cancelled)
        );
    });
}
