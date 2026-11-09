use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::material_state::water_refraction::pressure_temperature::PressureTemperatureJet;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:e} != {expected:e}; tolerance {tolerance:e}");
}
fn point(p: [f64; 3]) -> Point3 { Point3::new(p[0], p[1], p[2]) }
fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.5, max_step_m: 2.0,
        position_tolerance_m: 1e-11, direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-11,
        ..GradientIndexOptions::default()
    }
}

struct JointState;
impl PressureTemperatureGradientField for JointState {
    fn pressure_temperature_gradient(&self, p: Point3) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        Ok(self.pressure_temperature_jet(p)?.sample())
    }
}
impl PressureTemperatureSecondOrderField for JointState {
    fn pressure_temperature_jet(&self, p: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        let (x, y, z) = (p.x, p.y, p.z);
        Ok(PressureTemperatureJet {
            temperature: TemperatureJet {
                sample: TemperatureGradientSample {
                    temperature_k: 293.15 + 3.0*x + 2.0*y + 0.5*x*x + 0.4*x*y + 0.3*z*z,
                    gradient_k_per_m: Vec3::new(3.0 + x + 0.4*y, 2.0 + 0.4*x, 0.6*z),
                },
                diagonal_k_per_m2: [1.0, 0.0, 0.6],
                off_diagonal_k_per_m2: [0.4, 0.0, 0.0],
            },
            pressure_pa: 101_325.0 + 1000.0*x - 700.0*y + 500.0*z + 100.0*x*y + 50.0*y*z + 80.0*x*x,
            pressure_gradient_pa_per_m: Vec3::new(1000.0 + 100.0*y + 160.0*x, -700.0 + 100.0*x + 50.0*z, 500.0 + 50.0*y),
            pressure_diagonal_pa_per_m2: [160.0, 0.0, 0.0],
            pressure_off_diagonal_pa_per_m2: [100.0, 0.0, 50.0],
        })
    }
}

#[test]
fn g1_original_pressure_temperature_gradients_and_hessians_reach_air_optics() {
    let composition = CiddorAir::try_new(0.01, 700.0).unwrap();
    let field = MaterialAirRefractiveField::try_nonabsorbing(&JointState, composition, ContentHash([31; 32])).unwrap();
    let p = [0.1, -0.2, 0.3];
    assert!(!field.is_spatially_uniform());
    for wavelength in [300.0, 633.0, 1700.0] {
        let prepared = field.prepare(wavelength).unwrap();
        let jet = prepared.optical_jet(point(p)).unwrap();
        assert_eq!(jet.optical(), prepared.optical_sample(point(p)).unwrap());
        assert_eq!(jet.optical(), field.sample(point(p), wavelength).unwrap());
        let value = |position| {
            let state = JointState.pressure_temperature_gradient(point(position)).unwrap();
            composition.prepare(wavelength).unwrap().sample(state.temperature.temperature_k, state.pressure_pa).unwrap().phase_index_minus_one
        };
        let g = jet.optical().phase_gradient_per_m();
        let gradient = [g.x, g.y, g.z];
        let (h, center) = (0.002, value(p));
        for i in 0..3 {
            let (mut plus, mut minus) = (p, p);
            plus[i] += h; minus[i] -= h;
            close(gradient[i], (value(plus) - value(minus)) / (2.0*h), 2e-13);
            close(jet.hessian().diagonal()[i], (value(plus) - 2.0*center + value(minus)) / (h*h), 2e-13);
        }
        for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            let (mut pp, mut pm, mut mp, mut mm) = (p, p, p, p);
            pp[i] += h; pp[j] += h; pm[i] += h; pm[j] -= h;
            mp[i] -= h; mp[j] += h; mm[i] -= h; mm[j] -= h;
            close(jet.hessian().off_diagonal()[k], (value(pp)-value(pm)-value(mp)+value(mm))/(4.0*h*h), 2e-13);
        }
    }
}

struct PressureLens { curvature: f64 }
impl PressureTemperatureGradientField for PressureLens {
    fn pressure_temperature_gradient(&self, p: Point3) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        if p.x.abs() > 1.0 || p.y.abs() > 1.0 { return Err(GradientIndexError::OutsideFieldDomain); }
        Ok(PressureTemperatureGradientSample {
            temperature: TemperatureGradientSample { temperature_k: 293.15, gradient_k_per_m: Vec3::new(0.0, 0.0, 0.0) },
            pressure_pa: 100_000.0 + 0.5*self.curvature*(p.x*p.x + p.y*p.y),
            pressure_gradient_pa_per_m: Vec3::new(self.curvature*p.x, self.curvature*p.y, 0.0),
        })
    }
}
impl PressureTemperatureSecondOrderField for PressureLens {
    fn pressure_temperature_jet(&self, p: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        let state = self.pressure_temperature_gradient(p)?;
        Ok(PressureTemperatureJet {
            temperature: TemperatureJet { sample: state.temperature, diagonal_k_per_m2: [0.0; 3], off_diagonal_k_per_m2: [0.0; 3] },
            pressure_pa: state.pressure_pa, pressure_gradient_pa_per_m: state.pressure_gradient_pa_per_m,
            pressure_diagonal_pa_per_m2: [self.curvature, self.curvature, 0.0], pressure_off_diagonal_pa_per_m2: [0.0; 3],
        })
    }
}

#[test]
fn g1_pressure_alone_bends_air_rays_and_focuses_native_bundles() {
    let source = PressureLens { curvature: -10_000.0 };
    let air = CiddorAir::try_new(0.0, 450.0).unwrap();
    let field = MaterialAirRefractiveField::try_nonabsorbing(&source, air, ContentHash([32; 32])).unwrap();
    for wavelength in [400.0, 633.0, 1550.0] {
        let prepared = field.prepare(wavelength).unwrap();
        let original = air.prepare(wavelength).unwrap().sample(293.15, 100_000.0).unwrap();
        let omega = (-source.curvature * original.pressure_derivative_per_pa / original.phase_index).sqrt();
        let length = 40.0;
        let ray = Ray { origin: point([0.0; 3]), dir: Vec3::new(0.0, 0.0, 1.0) };
        let bundle = propagate_bundle(&prepared, ray, Vec3::new(1.0, 0.0, 0.0), wavelength, length,
            BundleOptions { integration: options(), position_derivative_tolerance_m: 1e-10, direction_derivative_tolerance: 1e-12 }).unwrap();
        let footprint = bundle.project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default()).unwrap();
        let expected = (omega * length).sin() / omega;
        close(footprint.area_per_solid_angle_m2(), expected*expected, 1e-7);
        assert!(footprint.area_per_solid_angle_m2() < length*length - 10.0);
        close(bundle.propagation().optical_path_length_m, length*original.phase_index, 2e-11);
        close(bundle.propagation().group_optical_path_length_m.unwrap(), length*original.group_index, 2e-11);
        close(bundle.propagation().power_transmittance(), 1.0, 0.0);
        let off_axis = Ray { origin: point([0.02, 0.0, 0.0]), ..ray };
        let initial = prepared.sample(off_axis.origin, wavelength).unwrap();
        let bent = propagate(&prepared, off_axis, Vec3::new(1.0, 0.0, 0.0), wavelength, 20.0, options()).unwrap();
        assert!(bent.ray.dir.x < -1e-5);
        // Translation symmetry in z conserves the original optical momentum.
        close(prepared.sample(bent.ray.origin, wavelength).unwrap().phase_index()*bent.ray.dir.z, initial.phase_index(), 2e-12);
    }
}

#[test]
fn g0_air_snapshot_and_source_domain_refusals_are_preserved() {
    let source = PressureLens { curvature: -10_000.0 };
    let air = CiddorAir::try_new(0.0, 450.0).unwrap();
    let bind = |id, composition| MaterialAirRefractiveField::try_nonabsorbing(&source, composition, ContentHash([id; 32]));
    assert!(matches!(bind(0, air), Err(AirFieldError::MissingSnapshotIdentity)));
    let field = bind(33, air).unwrap();
    assert_ne!(field.content_identity(), bind(34, air).unwrap().content_identity());
    assert_ne!(field.content_identity(), bind(33, CiddorAir::try_new(0.01, 450.0).unwrap()).unwrap().content_identity());
    assert_eq!(field.snapshot_identity(), ContentHash([33; 32]));
    assert!(core::ptr::eq(field.state_source(), &source));
    for p in [[2.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(field.sample(point(p), 550.0).is_err());
        assert!(field.sample_second_order(point(p), 550.0).is_err());
    }
    assert!(field.prepare(1700.1).is_err());
    let prepared = field.prepare(550.0).unwrap();
    assert!(prepared.sample(point([0.0; 3]), 550.1).is_err());
    assert!(prepared.sample_second_order(point([0.0; 3]), 550.1).is_err());
}
