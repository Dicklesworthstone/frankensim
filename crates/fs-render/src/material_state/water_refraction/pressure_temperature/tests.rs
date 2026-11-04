use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexField, GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::{TemperatureGradientField, TemperatureSecondOrderField};
use crate::material_state::water_refraction::isobaric::IsobaricWaterField;
use fs_material::water_refraction::IapwsWaterRefraction;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}; tolerance {tolerance:.3e}");
}
fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}
#[derive(Clone, Copy)]
struct Polynomial {
    value: f64,
    gradient: [f64; 3],
    hessian: [f64; 6],
}
impl Polynomial {
    fn uniform(value: f64) -> Self {
        Self { value, gradient: [0.0; 3], hessian: [0.0; 6] }
    }
    fn sample(self, p: Point3) -> Result<(f64, Vec3), GradientIndexError> {
        let x = [p.x, p.y, p.z];
        if !x.into_iter().all(|v| v.is_finite() && v.abs() <= 4.0) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let (mut value, mut gradient) = (self.value, self.gradient);
        for i in 0..3 {
            value += self.gradient[i] * x[i] + 0.5 * self.hessian[i] * x[i] * x[i];
            gradient[i] += self.hessian[i] * x[i];
        }
        for ((i, j), h) in [(0, 1), (0, 2), (1, 2)].into_iter().zip(&self.hessian[3..]) {
            value += h * x[i] * x[j];
            gradient[i] += h * x[j];
            gradient[j] += h * x[i];
        }
        Ok((value, Vec3::new(gradient[0], gradient[1], gradient[2])))
    }
}
impl TemperatureField for Polynomial {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.sample(p)?.0)
    }
}
impl TemperatureGradientField for Polynomial {
    fn temperature_gradient(&self, p: Point3) -> Result<TemperatureGradientSample, GradientIndexError> {
        let (temperature_k, gradient_k_per_m) = self.sample(p)?;
        Ok(TemperatureGradientSample { temperature_k, gradient_k_per_m })
    }
}
impl TemperatureSecondOrderField for Polynomial {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        Ok(TemperatureJet {
            sample: self.temperature_gradient(p)?,
            diagonal_k_per_m2: self.hessian[..3].try_into().unwrap(),
            off_diagonal_k_per_m2: self.hessian[3..].try_into().unwrap(),
        })
    }
}
struct Joint {
    temperature: Polynomial,
    pressure: Polynomial,
}
impl PressureTemperatureGradientField for Joint {
    fn pressure_temperature_gradient(&self, p: Point3) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        Ok(self.pressure_temperature_jet(p)?.sample())
    }
}
impl PressureTemperatureSecondOrderField for Joint {
    fn pressure_temperature_jet(&self, p: Point3) -> Result<PressureTemperatureJet, GradientIndexError> {
        let (pressure_pa, pressure_gradient_pa_per_m) = self.pressure.sample(p)?;
        Ok(PressureTemperatureJet {
            temperature: self.temperature.temperature_jet(p)?,
            pressure_pa,
            pressure_gradient_pa_per_m,
            pressure_diagonal_pa_per_m2: self.pressure.hessian[..3].try_into().unwrap(),
            pressure_off_diagonal_pa_per_m2: self.pressure.hessian[3..].try_into().unwrap(),
        })
    }
}
fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.005,
        max_step_m: 0.02,
        position_tolerance_m: 1e-12,
        direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-12,
        ..GradientIndexOptions::default()
    }
}

#[test]
fn g1_joint_pressure_temperature_optics_match_original_scalar_sources() {
    for (pressure, temperature, region, pressure_scale) in [
        (3e6, 300.0, If97Region::Region1, 1.0),
        (8e4, 400.0, If97Region::Region2, 0.02),
    ] {
        let source = Joint {
            temperature: Polynomial {
                value: temperature, gradient: [5.0, -3.0, 2.0],
                hessian: [8.0, -7.0, 2.0, 3.0, -4.0, 1.0],
            },
            pressure: Polynomial {
                value: pressure,
                gradient: [2e6, -1.3e6, 0.7e6].map(|v| v * pressure_scale),
                hessian: [7e5, -9e5, 2e5, 8e5, -6e5, 4e5].map(|v| v * pressure_scale),
            },
        };
        let field = PressureTemperatureWaterField::try_new(&source, region, ContentHash([91; 32])).unwrap();
        let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
        let optics = field.bind_refraction(&absorption).unwrap();
        assert!(core::ptr::eq(field.state_source(), &source));
        assert!(core::ptr::eq(optics.state_source(), &field));
        let p = [0.02, -0.04, 0.03];
        assert_eq!(field.sample_gradient(point(p)).unwrap(), field.sample_jet(point(p)).unwrap().sample());
        for lambda in [226.5, 589.0, 1013.98] {
            let original = IapwsWaterRefraction::prepare(lambda).unwrap();
            // Only original scalar pressure, temperature, IF97 and IAPWS
            // values enter this independent finite-difference reference.
            let phase = |x| {
                let state = If97State::at_pressure_temperature(
                    Pressure::new(source.pressure.sample(point(x)).unwrap().0),
                    Temperature::new(source.temperature.sample(point(x)).unwrap().0), region,
                ).unwrap();
                original.sample_if97_state(&state).unwrap().phase_index
            };
            let jet = optics.prepare(lambda).unwrap().optical_jet(point(p)).unwrap();
            let g = jet.optical().phase_gradient_per_m();
            let gradient = [g.x, g.y, g.z];
            let (n, h) = (phase(p), 0.002);
            for i in 0..3 {
                let (mut a, mut b) = (p, p);
                a[i] += h; b[i] -= h;
                near(gradient[i], (phase(a) - phase(b)) / (2.0 * h), 2e-9);
                near(jet.hessian().diagonal()[i], (phase(a) - 2.0 * n + phase(b)) / (h * h), 3e-9);
            }
            for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                let (mut pp, mut pm, mut mp, mut mm) = (p, p, p, p);
                pp[i] += h; pp[j] += h; pm[i] += h; pm[j] -= h;
                mp[i] -= h; mp[j] += h; mm[i] -= h; mm[j] -= h;
                near(jet.hessian().off_diagonal()[k],
                    (phase(pp) - phase(pm) - phase(mp) + phase(mm)) / (4.0 * h * h), 3e-9);
            }
        }
    }
}

#[test]
fn g1_constant_pressure_retains_the_existing_isobaric_water_jet() {
    let source = Joint {
        temperature: Polynomial {
            value: 300.0, gradient: [5.0, -3.0, 2.0],
            hessian: [8.0, -7.0, 2.0, 3.0, -4.0, 1.0],
        },
        pressure: Polynomial::uniform(3e6),
    };
    let field = PressureTemperatureWaterField::try_new(&source, If97Region::Region1, ContentHash([92; 32])).unwrap();
    let isobaric = IsobaricWaterField::try_new(&source.temperature, Pressure::new(3e6), If97Region::Region1, ContentHash([92; 32])).unwrap();
    for p in [[0.0; 3], [0.2, -0.3, 0.1]] {
        assert_eq!(field.sample_state(point(p)).unwrap(), isobaric.sample_state(point(p)).unwrap());
        assert_eq!(field.sample_gradient(point(p)).unwrap(), isobaric.sample_gradient(point(p)).unwrap());
        assert_eq!(field.sample_jet(point(p)).unwrap(), isobaric.sample_jet(point(p)).unwrap());
    }
}

#[test]
fn g1_pressure_lens_drives_native_bending_and_analytic_bundle_focusing() {
    let curvature = -300e6;
    let source = Joint {
        temperature: Polynomial::uniform(300.0),
        pressure: Polynomial {
            hessian: [curvature, curvature, 0.0, 0.0, 0.0, 0.0],
            ..Polynomial::uniform(40e6)
        },
    };
    let field = PressureTemperatureWaterField::try_new(&source, If97Region::Region1, ContentHash([93; 32])).unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let optics = field.bind_refraction(&absorption).unwrap();
    let state = field.sample_state(point([0.0; 3])).unwrap();
    for lambda in [400.0, 589.0, 900.0] {
        let original = IapwsWaterRefraction::prepare(lambda).unwrap().sample_if97_state(&state).unwrap();
        let n_p = original.density_derivative_m3_per_kg * state.density().value()
            * state.isothermal_compressibility().value();
        assert!(n_p > 0.0);
        let omega = (-curvature * n_p / original.phase_index).sqrt();
        let length = 2.0;
        let prepared = optics.prepare(lambda).unwrap();
        let bundle = propagate_bundle(
            &prepared, Ray { origin: point([0.0; 3]), dir: Vec3::new(0.0, 0.0, 1.0) },
            Vec3::new(1.0, 0.0, 0.0), lambda, length,
            BundleOptions { integration: options(), position_derivative_tolerance_m: 1e-12, direction_derivative_tolerance: 1e-12 },
        ).unwrap();
        let footprint = bundle.project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default()).unwrap();
        let expected = (omega * length).sin() / omega;
        near(footprint.area_per_solid_angle_m2(), expected * expected, 2e-9);
        assert!(footprint.area_per_solid_angle_m2() < length * length - 0.1);
        near(bundle.propagation().optical_path_length_m, length * original.phase_index, 1e-11);
        near(bundle.propagation().group_optical_path_length_m.unwrap(), length * original.group_index, 1e-11);
        near(bundle.propagation().power_transmittance(), (-0.02 * length).exp(), 1e-12);

        let start = point([0.02, 0.0, 0.0]);
        let initial = prepared.sample(start, lambda).unwrap();
        let bent = propagate(&prepared,
            Ray { origin: start, dir: Vec3::new(0.0, 0.0, 1.0) },
            Vec3::new(1.0, 0.0, 0.0), lambda, 1.0, options(),
        ).unwrap();
        assert!(bent.ray.dir.x < -1e-5);
        // The original pressure lens has no z dependence, so n*d_z is conserved.
        near(prepared.sample(bent.ray.origin, lambda).unwrap().phase_index() * bent.ray.dir.z,
            initial.phase_index(), 2e-11);
    }
}

#[test]
fn g3_joint_snapshot_phase_and_physical_domains_remain_explicit() {
    let source = Joint { temperature: Polynomial::uniform(300.0), pressure: Polynomial::uniform(1e5) };
    let make = |region, id| PressureTemperatureWaterField::try_new(&source, region, ContentHash([id; 32]));
    let field = make(If97Region::Region1, 94).unwrap();
    assert!(!field.is_spatially_uniform());
    assert_ne!(field.content_identity(), make(If97Region::Region1, 95).unwrap().content_identity());
    assert_ne!(field.content_identity(), make(If97Region::Region2, 94).unwrap().content_identity());
    assert!(matches!(make(If97Region::Region1, 0), Err(PressureTemperatureWaterError::MissingSnapshotIdentity)));
    for p in [[5.0, 0.0, 0.0], [f64::NAN, 0.0, 0.0]] {
        assert!(matches!(field.water_state_gradient(point(p)), Err(GradientIndexError::OutsideFieldDomain)));
    }
    assert!(matches!(make(If97Region::Region2, 94).unwrap().water_state_gradient(point([0.0; 3])), Err(GradientIndexError::OutsideFieldDomain)));
    for p in [0.0, -1.0, f64::NAN, 1.1e8] {
        let invalid = Joint { pressure: Polynomial::uniform(p), temperature: source.temperature };
        let water = PressureTemperatureWaterField::try_new(&invalid, If97Region::Region1, ContentHash([94; 32])).unwrap();
        assert!(water.water_state_gradient(point([0.0; 3])).is_err());
    }
    let crossing = Joint {
        temperature: Polynomial::uniform(373.15),
        pressure: Polynomial { gradient: [-1.5e5, 0.0, 0.0], ..Polynomial::uniform(2e5) },
    };
    let liquid = PressureTemperatureWaterField::try_new(&crossing, If97Region::Region1, ContentHash([96; 32])).unwrap();
    assert!(liquid.sample_state(point([0.0; 3])).is_ok());
    assert!(matches!(liquid.water_state_gradient(point([1.0, 0.0, 0.0])), Err(GradientIndexError::OutsideFieldDomain)));
    let hot = Joint { temperature: Polynomial::uniform(800.0), pressure: source.pressure };
    let vapor = PressureTemperatureWaterField::try_new(&hot, If97Region::Region2, ContentHash([97; 32])).unwrap();
    assert!(vapor.sample_state(point([0.0; 3])).is_ok());
    let absorption = BeerLambertAbsorption::try_constant(0.0).unwrap();
    let optics = vapor.bind_refraction(&absorption).unwrap();
    assert!(optics.sample(point([0.0; 3]), 589.0).is_err());
    assert!(optics.prepare(1101.0).is_err());
}

#[test]
fn g3_unrepresentable_pressure_hessian_refuses_without_poisoning_the_ray_sample() {
    let source = Joint {
        temperature: Polynomial::uniform(300.0),
        pressure: Polynomial { gradient: [1e308, 0.0, 0.0], ..Polynomial::uniform(1e5) },
    };
    let field = PressureTemperatureWaterField::try_new(&source, If97Region::Region1, ContentHash([98; 32])).unwrap();
    assert!(field.water_state_gradient(point([0.0; 3])).is_ok());
    assert!(matches!(field.water_state_jet(point([0.0; 3])), Err(GradientIndexError::NumericalRange)));
}
