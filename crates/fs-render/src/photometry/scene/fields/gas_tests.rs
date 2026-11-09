use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{BeerLambertAbsorption, GlassProvenance, SellmeierIor};
use crate::gradient_index::bundle::ProjectionOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::gradient_index::thermal::TemperatureField;
use crate::gradient_index::{GradientIndexError, GradientIndexOptions};
use crate::material_state::gas_scattering::isobaric::IsobaricGasField;
use crate::material_state::gas_scattering::{GasPartialDensityJet, GasPartialDensitySample};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::material_state::water_refraction::pressure_temperature::{
    PressureTemperatureGradientSample, PressureTemperatureJet,
};
use crate::photometry::scene::{
    NativeGeometryBinding, NativeOpticalRole, NativeRadianceCalibration,
};
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::CurvedTraceSettings;
use crate::polarization::dielectric::scene::curved::bundle::{
    CurvedBundleSceneError, CurvedBundleSettings,
};
use crate::polarization::dielectric::scene::curved::camera::RadianceEndpoint;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighBinding, CurvedRayleighField, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::spectral::{LiftedSpectrum, SpectralReflectance};
use crate::tracer::{Camera, Primitive, Shape};
use fs_blake3::hash_domain;
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3, fixtures::SphereChart};
use fs_material::air_refraction::CiddorAir;
use fs_material::gas::{GasSpec, R_USSA_1976};
use fs_material::gas_scattering::composition::{
    MolecularCauchyDispersion, MolecularRefractionModel, PartialDensityBasis,
};
use fs_material::gas_scattering::{AVOGADRO_PER_MOL, MolecularRayleigh, MolecularRayleighMixture};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticCoefficients};
use fs_qty::semantic::{Composition, CompositionBasis};
use fs_qty::{Pressure, Temperature};

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
                seed: 89,
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

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        SellmeierIor::try_new_in_band([1.25, 0.0, 0.0], [0.0; 3], (300.0, 2500.0)).unwrap(),
        BeerLambertAbsorption::try_constant(0.2).unwrap(),
        GlassProvenance::Custom,
    )
}

fn scene(
    z: f64,
    center_x: f64,
    width: f64,
    lens: bool,
    detector: bool,
) -> (Scene, Vec<NativeOpticalBinding>) {
    let mut primitives = Vec::new();
    if lens {
        primitives.push(Primitive {
            shape: Shape::Chart(Box::new(SphereChart {
                center: Point3::new(0.0, 0.0, 4.0),
                radius: 1.0,
            })),
            material: Material::Dielectric {
                glass: glass(),
                surface: DielectricSurface::SMOOTH,
            },
            emission: None,
        });
    }
    primitives.push(Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [center_x - width, -width, z],
                [center_x + width, -width, z],
                [center_x, 2.0 * width, z],
            ],
            vec![[0, 1, 2]],
        )),
        material: Material::SpectralDiffuse {
            reflectance: SpectralReflectance::constant(0.0).unwrap(),
        },
        emission: (!detector).then_some((SPECTRUM, 2.0)),
    });
    let bindings = (0..primitives.len())
        .map(|i| {
            let id = if lens && i == 0 { 1_u64 } else { 9 };
            NativeOpticalBinding {
                geometry: NativeGeometryBinding::Unplaced {
                    object_id: id,
                    geometry_identity: hash_domain("native-gas-scene-test", &id.to_le_bytes()),
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
            }
        })
        .collect();
    (
        Scene {
            primitives,
            lights: Vec::new(),
            environment: None,
            camera: Camera {
                eye: Point3::new(0.0, 0.0, 0.0),
                forward: Vec3::new(0.0, 0.0, 1.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                half_tan: 0.05,
            },
        },
        bindings,
    )
}

fn admit(
    cx: &Cx<'_>,
    ambient: NativeGasAmbient<'_>,
    original: (Scene, Vec<NativeOpticalBinding>),
) -> NativeGasOpticalScene {
    NativeGasOpticalScene::try_from_scene(cx, original.0, ambient, &original.1, 2).unwrap()
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 64.0,
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
    }
}

fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (tolerance {tolerance:.3e})"
    );
}

// Declared synthetic molecular optical sources, not measured atmosphere claims.
fn density_model(source: u8) -> MolecularRayleighMixture {
    MolecularRayleighMixture::try_new(
        &[GasSpec::dry_air_ussa1976()],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        &[
            MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([source; 32]))
                .unwrap(),
        ],
        None,
    )
    .unwrap()
}

struct ThermalProfile {
    slope: f64,
    curvature: f64,
}
impl TemperatureField for ThermalProfile {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.temperature_gradient(p)?.temperature_k)
    }
}
impl TemperatureGradientField for ThermalProfile {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        if ![p.x, p.y, p.z].into_iter().all(f64::is_finite) || p.x.abs() > 0.5 || p.y.abs() > 0.5 {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        Ok(TemperatureGradientSample {
            temperature_k: 300.0
                + self.slope * p.x
                + 0.5 * self.curvature * (p.x * p.x + p.y * p.y),
            gradient_k_per_m: Vec3::new(
                self.slope + self.curvature * p.x,
                self.curvature * p.y,
                0.0,
            ),
        })
    }
}
impl TemperatureSecondOrderField for ThermalProfile {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: [self.curvature, self.curvature, 0.0],
            off_diagonal_k_per_m2: [0.0; 3],
        })
    }
}

// Independent exact isobaric solution. With T=300+g*x and L=C/T,
// n²=(T+2C)/(T-C), and conserved n*u_z=n(0). Substitute x=-q² to
// integrate the endpoint, physical length and phase without an ODE solver.
fn isobaric_reference() -> (f64, f64, f64, f64, f64) {
    let r = 0.0003;
    let volume = (r * (2.0 + r) / (3.0 + r * (2.0 + r))) / 2.5e25;
    let c = volume * 100_000.0 * AVOGADRO_PER_MOL / R_USSA_1976;
    let (a, b, g) = (300.0 - c, 300.0 + 2.0 * c, 300.0_f64);
    let n0 = (b / a).sqrt();
    let scale = (a / (3.0 * c * g)).sqrt();
    let integral = |q: f64, t: f64| {
        0.5 * (q * (t - g * q * q).sqrt() + t / g.sqrt() * (q * (g / t).sqrt()).asin())
    };
    let (mut lo, mut hi) = (0.0, 0.65);
    for _ in 0..90 {
        let q = 0.5 * (lo + hi);
        if 2.0 * n0 * scale * integral(q, a) < 40.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let n1 = ((b - g * q * q) / (a - g * q * q)).sqrt();
    let distance = 2.0 * scale * integral(q, b);
    let phase = 2.0 * scale * (integral(q, a) + 3.0 * c / g.sqrt() * (q * (g / a).sqrt()).asin());
    (-q * q, distance, phase, n0, n1)
}

#[test]
fn g1_original_isobaric_gas_bends_camera_to_original_off_axis_source() {
    with_cx(|cx, _| {
        let model = density_model(91);
        let temperature = ThermalProfile {
            slope: 300.0,
            curvature: 0.0,
        };
        let gas = IsobaricGasField::try_new(
            &temperature,
            &model,
            Pressure::new(100_000.0),
            Temperature::new(300.0),
            ContentHash([92; 32]),
        )
        .unwrap();
        let field = gas.bind_refraction().unwrap();
        let first = NativeGasAmbient::from_density_field(&field);
        let second = NativeGasAmbient::from_differentiable_density_field(&field);
        let owner = admit(cx, first, scene(40.0, -0.25, 0.12, false, false));
        assert!(owner.medium(0).is_none());
        assert_eq!(
            first.material_identity(),
            model.evaluate(300.0, 100_000.0).unwrap().model_identity()
        );
        assert_eq!(first.field_identity(), second.field_identity());
        let (x, distance, phase, n0, n1) = isobaric_reference();
        assert!((-0.33..-0.17).contains(&x)); // The straight axial ray misses this original triangle.
        for ambient in [first, second] {
            owner
                .with_curved_camera(
                    cx,
                    100_000,
                    NativeRadianceCalibration::try_new(3.0, 0.0).unwrap(),
                    ambient,
                    &[],
                    |native, camera| {
                        assert_eq!(camera.eye(), ray().origin);
                        for wavelength in [401.25, 550.0, 699.125] {
                            let result = native
                                .sample_ray_with_samples(
                                    cx,
                                    ray(),
                                    wavelength,
                                    0,
                                    settings(),
                                    |_| panic!("black source has no continuation"),
                                )
                                .unwrap();
                            let path = result.path();
                            assert_eq!(
                                result.endpoint(),
                                RadianceEndpoint::Surface { object_id: 9 }
                            );
                            assert_eq!(
                                path.medium_snapshot_identity(),
                                Some(field.content_identity())
                            );
                            close(path.beam().ray().origin.x, x, 5e-7);
                            close(path.beam().ray().origin.z, 40.0, 5e-7);
                            close(path.beam().ray().dir.z, n0 / n1, 2e-9);
                            close(path.traveled_m(), distance, 1e-6);
                            close(path.optical_path_length_m(), phase, 1e-6);
                            close(path.group_optical_path_length_m().unwrap(), phase, 1e-6);
                            assert_eq!(path.optical_depth(), 0.0);
                            let expected = 6.0 * SPECTRUM.eval(wavelength) * (n0 / n1).powi(2);
                            close(result.radiance_w_m2_sr_nm(), expected, expected * 1e-8);
                        }
                    },
                )
                .unwrap();
        }
    });
}

struct AtmosphericProfile(ThermalProfile);
impl PressureTemperatureGradientField for AtmosphericProfile {
    fn pressure_temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        Ok(PressureTemperatureGradientSample {
            temperature: self.0.temperature_gradient(p)?,
            pressure_pa: 101_325.0 + 2000.0 * p.x,
            pressure_gradient_pa_per_m: Vec3::new(2000.0, 0.0, 0.0),
        })
    }
}
impl PressureTemperatureSecondOrderField for AtmosphericProfile {
    fn pressure_temperature_jet(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureJet, GradientIndexError> {
        let state = self.pressure_temperature_gradient(p)?;
        Ok(PressureTemperatureJet {
            temperature: self.0.temperature_jet(p)?,
            pressure_pa: state.pressure_pa,
            pressure_gradient_pa_per_m: state.pressure_gradient_pa_per_m,
            pressure_diagonal_pa_per_m2: [0.0; 3],
            pressure_off_diagonal_pa_per_m2: [0.0; 3],
        })
    }
}

#[test]
fn g1_ciddor_atmosphere_bends_native_camera_to_original_source_with_group_delay() {
    with_cx(|cx, _| {
        let state = AtmosphericProfile(ThermalProfile {
            slope: 100.0,
            curvature: 0.0,
        });
        let composition = CiddorAir::try_new(0.01, 450.0).unwrap();
        let field = MaterialAirRefractiveField::try_nonabsorbing(
            &state,
            composition,
            ContentHash([141; 32]),
        )
        .unwrap();
        let first = NativeGasAmbient::from_ciddor_field(&field);
        let second = NativeGasAmbient::from_differentiable_ciddor_field(&field);
        // A straight axial ray misses the original off-axis emitting triangle.
        let owner = admit(cx, first, scene(40.0, -0.07, 0.025, false, false));
        assert_eq!(first.material_identity(), composition.identity());
        assert_eq!(first.wavelength_range_nm(), [300.0, 1700.0]);
        for wavelength in [450.0, 650.0] {
            let prepared = field.prepare(wavelength).unwrap();
            let third = NativeGasAmbient::from_prepared_ciddor_field(&prepared);
            let fourth = NativeGasAmbient::from_prepared_differentiable_ciddor_field(&prepared);
            for ambient in [first, second, third, fourth] {
                assert_eq!(ambient.field_identity(), field.content_identity());
                owner
                    .with_curved_camera(
                        cx,
                        100_000,
                        NativeRadianceCalibration::try_new(3.0, 0.0).unwrap(),
                        ambient,
                        &[],
                        |native, _| {
                            let result = native
                                .sample_ray_with_samples(
                                    cx,
                                    ray(),
                                    wavelength,
                                    0,
                                    settings(),
                                    |_| panic!("black source has no continuation"),
                                )
                                .unwrap();
                            assert_eq!(
                                result.endpoint(),
                                RadianceEndpoint::Surface { object_id: 9 }
                            );
                            let path = result.path();
                            let end = path.beam().ray();
                            assert!((-0.095..-0.045).contains(&end.origin.x));
                            close(end.origin.z, 40.0, 5e-7);
                            assert_eq!(
                                path.medium_snapshot_identity(),
                                Some(field.content_identity())
                            );
                            let n0 = prepared.optical_sample(ray().origin).unwrap().phase_index();
                            let n1 = prepared.optical_sample(end.origin).unwrap().phase_index();
                            close(n1 * end.dir.z, n0, 2e-9);
                            assert!(
                                path.group_optical_path_length_m().unwrap()
                                    > path.optical_path_length_m() + 1e-4
                            );
                            assert_eq!(path.optical_depth(), 0.0);
                            let expected = 6.0 * SPECTRUM.eval(wavelength) * (n0 / n1).powi(2);
                            close(result.radiance_w_m2_sr_nm(), expected, expected * 1e-8);
                        },
                    )
                    .unwrap();
            }
        }
    });
}

fn composition_model(source: u8) -> MolecularRefractionModel {
    let mut specs = [GasSpec::dry_air_ussa1976(); 2];
    specs[1].molar_mass *= 2.0;
    let optics = [0.00015, 0.0006].map(|r| {
        MolecularRayleigh::try_new(r, 2.5e25, 550.0, 0.03, ContentHash([source; 32])).unwrap()
    });
    let original = MolecularRayleighMixture::try_new(
        &specs,
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
        &optics,
        None,
    )
    .unwrap();
    let reference = MolecularRefractionModel::from_mixture(
        &original,
        original.evaluate(300.0, 100_000.0).unwrap().gas(),
    )
    .unwrap();
    let sources: Vec<_> = reference
        .components()
        .iter()
        .enumerate()
        .map(|(i, row)| {
            MolecularCauchyDispersion::try_new_in_band(
                row.optical(),
                [2e-6, 8e-6][i],
                [1e-7, 3e-7][i],
                [300.0, 2500.0],
                ContentHash([source + i as u8 + 1; 32]),
            )
            .unwrap()
        })
        .collect();
    reference.with_dispersion(&sources).unwrap()
}

struct PartialDensities {
    model: ContentHash,
    curvature: f64,
}
impl PartialDensities {
    fn jet(&self, p: Point3, i: usize) -> Result<GasPartialDensityJet, GradientIndexError> {
        if ![p.x, p.y, p.z].into_iter().all(f64::is_finite) || p.x.abs() > 1.0 || p.y.abs() > 1.0 {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let k = if i == 0 {
            self.curvature
        } else {
            -self.curvature
        };
        Ok(GasPartialDensityJet {
            sample: GasPartialDensitySample {
                density: 1.0 + 0.5 * k * (p.x * p.x + p.y * p.y),
                gradient_per_m: Vec3::new(k * p.x, k * p.y, 0.0),
            },
            hessian_diagonal_per_m2: [k, k, 0.0],
            hessian_off_diagonal_per_m2: [0.0; 3],
        })
    }
}
impl GasCompositionGradientField for PartialDensities {
    fn model_identity(&self) -> ContentHash {
        self.model
    }
    fn component_count(&self) -> usize {
        2
    }
    fn density_basis(&self) -> PartialDensityBasis {
        PartialDensityBasis::MassKgPerM3
    }
    fn sample_partial_densities(
        &self,
        p: Point3,
        output: &mut [GasPartialDensitySample],
    ) -> Result<(), GradientIndexError> {
        for (i, row) in output.iter_mut().enumerate() {
            *row = self.jet(p, i)?.sample;
        }
        Ok(())
    }
}
impl GasCompositionSecondOrderField for PartialDensities {
    fn sample_partial_density_jets(
        &self,
        p: Point3,
        output: &mut [GasPartialDensityJet],
    ) -> Result<(), GradientIndexError> {
        for (i, row) in output.iter_mut().enumerate() {
            *row = self.jet(p, i)?;
        }
        Ok(())
    }
}

// Original source coefficients, independently differentiated before the mixture
// Lorentz--Lorenz transform. No production refractive sample/jet is an oracle.
fn composition_reference(model: &MolecularRefractionModel, wavelength: f64) -> (f64, f64, f64) {
    let lambda = wavelength * 1e-3;
    let mut coefficient = [0.0; 2];
    let mut derivative = 0.0;
    for i in 0..2 {
        let (b, c) = ([2e-6, 8e-6][i], [1e-7, 3e-7][i]);
        let r = [0.00015, 0.0006][i]
            + b * (lambda.powi(-2) - 0.55_f64.powi(-2))
            + c * (lambda.powi(-4) - 0.55_f64.powi(-4));
        let dr = -2.0 * b / lambda.powi(3) - 4.0 * c / lambda.powi(5);
        let delta = r * (2.0 + r);
        let scale = AVOGADRO_PER_MOL / (2.5e25 * model.components()[i].molar_mass_kg_mol());
        coefficient[i] = scale * delta / (3.0 + delta);
        derivative += scale * 6.0 * (1.0 + r) * dr / (3.0 + delta).powi(2);
    }
    let l = coefficient[0] + coefficient[1];
    let n = ((1.0 + 2.0 * l) / (1.0 - l)).sqrt();
    let dn = 3.0 / (2.0 * n * (1.0 - l).powi(2));
    (
        n,
        n - lambda * dn * derivative,
        dn * (coefficient[0] - coefficient[1]) * 0.5,
    )
}

fn focus([x, u]: [f64; 2], n: f64, curvature: f64, distance: f64) -> [f64; 2] {
    let gamma = (-curvature / n).sqrt();
    let angle = gamma * distance;
    [
        x * angle.cos() + u * angle.sin() / gamma,
        u * angle.cos() - x * gamma * angle.sin(),
    ]
}

#[test]
fn g1_original_chromatic_gas_and_thermal_glass_share_actual_bundle_geometry() {
    with_cx(|cx, _| {
        let model = composition_model(93);
        let density = PartialDensities {
            model: model.content_identity(),
            curvature: 0.5,
        };
        let field =
            MaterialGasCompositionRefractiveField::try_new(&model, &density, ContentHash([96; 32]))
                .unwrap();
        let ambient = NativeGasAmbient::from_differentiable_composition_field(&field);
        let owner = admit(cx, ambient, scene(8.0, 0.0, 4.0, true, true));
        let interior = owner
            .medium_index(NativeOpticalMedium::InteriorOf(0))
            .unwrap();
        assert_eq!(owner.medium(interior), Some(&glass()));
        let temperature = ThermalProfile {
            slope: 0.0,
            curvature: -12.0,
        };
        let law = SchottThermoOpticLaw::try_new(
            ThermoOpticCoefficients {
                d: [0.0003, 0.0, 0.0],
                e: [0.0; 2],
                resonance_um: 0.0,
            },
            300.0,
            [200.0, 450.0],
            [300.0, 2500.0],
            ContentHash([97; 32]),
        )
        .unwrap();
        let glass_field = MaterialThermoOpticField::try_new(
            owner.medium(interior).unwrap(),
            &law,
            &temperature,
            ContentHash([98; 32]),
        )
        .unwrap();
        let thermal = [NativeThermoOpticBinding::from_differentiable_field(
            NativeOpticalMedium::InteriorOf(0),
            &glass_field,
        )];
        let mut areas = Vec::new();
        for wavelength in [401.25, 699.125, 1550.0, 2500.0] {
            let observe = |native: &CurvedPolarizedScene<'_>| {
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
                .with_curved_scene(cx, 100_000, ambient, &thermal, observe)
                .unwrap();
            let prepared = field.at_wavelength(wavelength).unwrap();
            let prepared_ambient =
                NativeGasAmbient::from_prepared_differentiable_composition_field(&field, &prepared)
                    .unwrap();
            let repeated = owner
                .with_curved_scene(cx, 100_000, prepared_ambient, &thermal, observe)
                .unwrap();
            assert_eq!(path, repeated);
            assert_eq!(path.termination(), PathTermination::Detected);
            assert_eq!(path.last_contact().unwrap().object_id, 9);
            assert_eq!(
                path.medium_snapshot_identity(),
                Some(field.content_identity())
            );
            let (n, ng, curvature) = composition_reference(&model, wavelength);
            let mut column = focus([0.0, 1.0], n, curvature, 3.0);
            let ratio = n / 1.5;
            column[1] = ratio * column[1] - (1.0 - ratio) * column[0];
            column = focus(column, 1.5, -0.0015, 2.0);
            let ratio = 1.5 / n;
            column[1] = ratio * column[1] + (1.0 - ratio) * column[0];
            column = focus(column, n, curvature, 3.0);
            let bundle = path.bundle().propagation();
            let area = bundle
                .project_at_endpoint(ray().dir, ProjectionOptions::default())
                .unwrap()
                .area_per_solid_angle_m2();
            close(area, column[0] * column[0], 5e-6);
            close(
                bundle.propagation().optical_path_length_m,
                6.0 * n + 3.0,
                1e-6,
            );
            close(
                bundle.propagation().group_optical_path_length_m.unwrap(),
                6.0 * ng + 3.0,
                1e-6,
            );
            close(bundle.propagation().optical_depth, 0.4, 1e-7);
            let transmission = 4.0 * n * 1.5 / (n + 1.5).powi(2);
            close(
                path.bundle().interface_log_power().exp(),
                transmission * transmission,
                1e-9,
            );
            areas.push(area);
        }
        assert!((areas[0] - areas[3]).abs() > 1e-4);
        owner
            .with_curved_scene(
                cx,
                100_000,
                NativeGasAmbient::from_composition_field(&field),
                &thermal,
                |native| {
                    assert!(matches!(
                        native.launch_bundle(
                            ray(),
                            Vec3::new(1.0, 0.0, 0.0),
                            550.0,
                            Stokes::unpolarized(1.0),
                            0
                        ),
                        Err(CurvedBundleSceneError::MissingSecondOrderField { medium: 0 })
                    ));
                },
            )
            .unwrap();
    });
}

#[test]
fn g1_original_gas_camera_keeps_source_coupled_rayleigh_collisions_elastic() {
    with_cx(|cx, _| {
        let model = composition_model(101);
        let density = PartialDensities {
            model: model.content_identity(),
            curvature: 0.0,
        };
        let field = MaterialGasCompositionRefractiveField::try_new(
            &model,
            &density,
            ContentHash([104; 32]),
        )
        .unwrap();
        let scattering = field.bind_scattering_field(&[1.0, 1.0]).unwrap();
        let ambient = NativeGasAmbient::from_composition_field(&field);
        let owner = admit(cx, ambient, scene(1.0, 0.0, 4.0, false, false));
        let calibration = NativeRadianceCalibration::try_new(3.0, 0.0).unwrap();
        owner
            .with_curved_camera(cx, 100_000, calibration, ambient, &[], |camera, _| {
                let fields = [Some(CurvedRayleighBinding::Field(&scattering))];
                let transport =
                    CurvedRayleighScene::try_new_with_fields(cx, camera.source_scene(), &fields)
                        .unwrap();
                let molecular = CurvedRayleighCameraScene::try_new(cx, camera, &transport).unwrap();
                let rate = scattering.scattering_majorant_per_m(550.0).unwrap();
                let first_cdf = -(-0.25 * rate).exp_m1();
                let sample = molecular
                    .sample_ray_with_tracking_samples(
                        cx,
                        ray(),
                        550.0,
                        0,
                        CurvedRayleighTraceSettings {
                            surfaces: settings(),
                            max_volume_events: 4,
                        },
                        |_| panic!("black original emitter has no continuation"),
                        |flight| {
                            [
                                if flight == 0 { first_cdf } else { 0.999999 },
                                0.999999,
                                0.3,
                                0.0,
                            ]
                        },
                    )
                    .unwrap();
                assert_eq!(sample.volume_events(), 1);
                assert_eq!(sample.free_flight_draws(), 2);
                let radiance = sample.radiance_sample();
                assert_eq!(
                    radiance.endpoint(),
                    RadianceEndpoint::Surface { object_id: 9 }
                );
                assert_eq!(
                    radiance.path().medium_snapshot_identity(),
                    Some(field.content_identity())
                );
                assert_eq!(radiance.path().optical_depth(), 0.0);
                // Collision/survival probabilities are in the sampler, not an
                // additional Beer absorption of these elastic molecular events.
                close(
                    radiance.radiance_w_m2_sr_nm(),
                    6.0 * SPECTRUM.eval(550.0),
                    1e-10,
                );
            })
            .unwrap();
        let changed = MaterialGasCompositionRefractiveField::try_new(
            &model,
            &density,
            ContentHash([105; 32]),
        )
        .unwrap();
        owner
            .with_curved_camera(
                cx,
                100_000,
                calibration,
                NativeGasAmbient::from_composition_field(&changed),
                &[],
                |camera, _| {
                    let wrong_fields = [Some(CurvedRayleighBinding::Field(&scattering))];
                    assert!(
                        CurvedRayleighScene::try_new_with_fields(
                            cx,
                            camera.source_scene(),
                            &wrong_fields
                        )
                        .is_err()
                    );
                },
            )
            .unwrap();
    });
}

#[test]
fn g3_gas_ambient_requires_original_material_snapshot_domain_and_budget() {
    with_cx(|cx, gate| {
        let model = composition_model(111);
        let density = PartialDensities {
            model: model.content_identity(),
            curvature: 0.0,
        };
        let field = MaterialGasCompositionRefractiveField::try_new(
            &model,
            &density,
            ContentHash([114; 32]),
        )
        .unwrap();
        let ambient = NativeGasAmbient::from_differentiable_composition_field(&field);
        let owner = admit(cx, ambient, scene(1.0, 0.0, 4.0, false, true));
        assert_eq!(ambient.wavelength_range_nm(), [300.0, 2500.0]);
        assert_eq!(
            owner.scene.with_curved_scene(cx, 100_000, |_| ()),
            Err(NativeOpticalSceneError::MissingGasAmbient)
        );
        let foreign_model = composition_model(121); // Same numerical optics, different original sources.
        let foreign_density = PartialDensities {
            model: foreign_model.content_identity(),
            curvature: 0.0,
        };
        let foreign = MaterialGasCompositionRefractiveField::try_new(
            &foreign_model,
            &foreign_density,
            ContentHash([114; 32]),
        )
        .unwrap();
        assert_eq!(
            owner.with_curved_scene(
                cx,
                100_000,
                NativeGasAmbient::from_composition_field(&foreign),
                &[],
                |_| ()
            ),
            Err(NativeOpticalSceneError::GasMaterialMismatch)
        );
        let changed = MaterialGasCompositionRefractiveField::try_new(
            &model,
            &density,
            ContentHash([115; 32]),
        )
        .unwrap();
        let prepared = changed.at_wavelength(550.0).unwrap();
        assert!(matches!(
            NativeGasAmbient::from_prepared_composition_field(&field, &prepared),
            Err(NativeOpticalSceneError::GasSnapshotMismatch)
        ));
        let prepared = field.at_wavelength(550.0).unwrap();
        let selected =
            NativeGasAmbient::from_prepared_differentiable_composition_field(&field, &prepared)
                .unwrap();
        for (gas, wavelength, origin) in [
            (selected, 551.0, ray().origin),
            (ambient, 2500.0_f64.next_up(), ray().origin),
            (ambient, 550.0, Point3::new(2.0, 0.0, 0.0)),
        ] {
            owner
                .with_curved_scene(cx, 100_000, gas, &[], |native| {
                    assert!(
                        native
                            .launch_bundle(
                                Ray { origin, ..ray() },
                                Vec3::new(1.0, 0.0, 0.0),
                                wavelength,
                                Stokes::unpolarized(1.0),
                                0
                            )
                            .is_err()
                    );
                })
                .unwrap();
        }
        let reference = glass();
        let temperature = ThermalProfile {
            slope: 0.0,
            curvature: 0.0,
        };
        let law = SchottThermoOpticLaw::try_new(
            ThermoOpticCoefficients {
                d: [0.0; 3],
                e: [0.0; 2],
                resonance_um: 0.0,
            },
            300.0,
            [200.0, 450.0],
            [300.0, 2500.0],
            ContentHash([116; 32]),
        )
        .unwrap();
        let thermal = MaterialThermoOpticField::try_new(
            &reference,
            &law,
            &temperature,
            ContentHash([117; 32]),
        )
        .unwrap();
        assert_eq!(
            owner.with_curved_scene(
                cx,
                100_000,
                ambient,
                &[NativeThermoOpticBinding::from_gradient_field(
                    NativeOpticalMedium::Ambient,
                    &thermal
                )],
                |_| ()
            ),
            Err(NativeOpticalSceneError::DuplicateFieldBinding {
                medium: NativeOpticalMedium::Ambient
            })
        );
        assert!(
            owner
                .with_curved_scene(cx, 0, ambient, &[], |_| ())
                .is_err()
        );
        gate.request();
        assert_eq!(
            owner.with_curved_scene(cx, 100_000, ambient, &[], |_| ()),
            Err(NativeOpticalSceneError::Cancelled)
        );
    });
}
