use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:.3e}"
    );
}
fn xyz(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn point(x: [f64; 3]) -> Point3 {
    Point3::new(x[0], x[1], x[2])
}

#[derive(Clone, Copy)]
struct Polynomial {
    t: f64,
    r: f64,
    gt: [f64; 3],
    gr: [f64; 3],
    ht: [f64; 6],
    hr: [f64; 6],
}
impl Polynomial {
    fn uniform() -> Self {
        Self {
            t: 300.0,
            r: 995.0,
            gt: [0.0; 3],
            gr: [0.0; 3],
            ht: [0.0; 6],
            hr: [0.0; 6],
        }
    }
    fn scalar(base: f64, g: [f64; 3], h: [f64; 6], p: [f64; 3]) -> (f64, Vec3) {
        let (mut value, mut gradient) = (base, g);
        for i in 0..3 {
            value += g[i] * p[i] + 0.5 * h[i] * p[i] * p[i];
            gradient[i] += h[i] * p[i];
        }
        for ((i, j), v) in [(0, 1), (0, 2), (1, 2)].into_iter().zip(&h[3..]) {
            value += v * p[i] * p[j];
            gradient[i] += v * p[j];
            gradient[j] += v * p[i];
        }
        (value, Vec3::new(gradient[0], gradient[1], gradient[2]))
    }
}
impl WaterStateGradientField for Polynomial {
    fn water_state_gradient(
        &self,
        p: Point3,
    ) -> Result<WaterStateGradientSample, GradientIndexError> {
        Ok(self.water_state_jet(p)?.sample())
    }
}
impl WaterStateSecondOrderField for Polynomial {
    fn water_state_jet(&self, p: Point3) -> Result<WaterStateJet, GradientIndexError> {
        if ![p.x, p.y, p.z]
            .into_iter()
            .all(|x| x.is_finite() && x.abs() <= 4.0)
        {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let (temperature_k, gradient_k_per_m) =
            Self::scalar(self.t, self.gt, self.ht, [p.x, p.y, p.z]);
        let (density_kg_m3, density_gradient_kg_m4) =
            Self::scalar(self.r, self.gr, self.hr, [p.x, p.y, p.z]);
        Ok(WaterStateJet {
            temperature: TemperatureJet {
                sample: TemperatureGradientSample {
                    temperature_k,
                    gradient_k_per_m,
                },
                diagonal_k_per_m2: self.ht[..3].try_into().unwrap(),
                off_diagonal_k_per_m2: self.ht[3..].try_into().unwrap(),
            },
            density_kg_m3,
            density_gradient_kg_m4,
            density_diagonal_kg_m5: self.hr[..3].try_into().unwrap(),
            density_off_diagonal_kg_m5: self.hr[3..].try_into().unwrap(),
        })
    }
}

fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.005,
        max_step_m: 0.02,
        position_tolerance_m: 1.0e-12,
        direction_tolerance: 1.0e-12,
        optical_path_tolerance_m: 1.0e-12,
        ..GradientIndexOptions::default()
    }
}
fn ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}

#[test]
fn full_temperature_density_chain_rule_matches_independent_spatial_differences() {
    let state = Polynomial {
        gt: [3.0, -2.0, 1.0],
        gr: [-5.0, 4.0, 2.0],
        ht: [0.7, -0.4, 0.2, 0.3, -0.2, 0.1],
        hr: [1.0, -0.6, 0.3, -0.4, 0.2, 0.5],
        ..Polynomial::uniform()
    };
    let absorption = BeerLambertAbsorption::try_constant(0.1).unwrap();
    let field =
        MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([74; 32])).unwrap();
    let p = [0.2, -0.1, 0.3];
    for lambda in [226.5, 589.0, 1013.98] {
        let prepared = field.prepare(lambda).unwrap();
        let jet = prepared.optical_jet(point(p)).unwrap();
        assert_eq!(jet.optical(), field.sample(point(p), lambda).unwrap());
        let phase = |p| prepared.optical_sample(point(p)).unwrap().phase_index();
        let gradient = xyz(jet.optical().phase_gradient_per_m());
        let n = phase(p);
        let h = 0.002;
        for i in 0..3 {
            let mut a = p;
            let mut b = p;
            a[i] += h;
            b[i] -= h;
            near(gradient[i], (phase(a) - phase(b)) / (2.0 * h), 2.0e-10);
            near(
                jet.hessian().diagonal()[i],
                (phase(a) - 2.0 * n + phase(b)) / (h * h),
                3.0e-9,
            );
        }
        for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            let mut pp = p;
            let mut pm = p;
            let mut mp = p;
            let mut mm = p;
            pp[i] += h;
            pp[j] += h;
            pm[i] += h;
            pm[j] -= h;
            mp[i] -= h;
            mp[j] += h;
            mm[i] -= h;
            mm[j] -= h;
            near(
                jet.hessian().off_diagonal()[k],
                (phase(pp) - phase(pm) - phase(mp) + phase(mm)) / (4.0 * h * h),
                2.0e-9,
            );
        }
    }
}

#[test]
fn density_curvature_focuses_native_water_bundles_even_at_zero_local_gradient() {
    let state = Polynomial {
        hr: [-40.0, -40.0, 0.0, 0.0, 0.0, 0.0],
        ..Polynomial::uniform()
    };
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let field =
        MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([75; 32])).unwrap();
    assert!(!field.is_spatially_uniform());
    for lambda in [226.5, 589.0, 1013.98] {
        let prepared = field.prepare(lambda).unwrap();
        let sample = IapwsWaterRefraction::prepare(lambda)
            .unwrap()
            .sample(state.t, state.r)
            .unwrap();
        let omega = (40.0 * sample.density_derivative_m3_per_kg / sample.phase_index).sqrt();
        let length = 2.0;
        let bundle = propagate_bundle(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            lambda,
            length,
            BundleOptions {
                integration: options(),
                position_derivative_tolerance_m: 1.0e-12,
                direction_derivative_tolerance: 1.0e-12,
            },
        )
        .unwrap();
        let footprint = bundle
            .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
            .unwrap();
        let expected = (omega * length).sin() / omega;
        near(
            footprint.area_per_solid_angle_m2(),
            expected * expected,
            1.0e-9,
        );
        assert!(footprint.area_per_solid_angle_m2() < length * length - 0.01);
        near(
            bundle.propagation().optical_path_length_m,
            length * sample.phase_index,
            1.0e-11,
        );
        near(
            bundle.propagation().group_optical_path_length_m.unwrap(),
            length * sample.group_index,
            1.0e-11,
        );
        near(
            bundle.propagation().power_transmittance(),
            (-0.02 * length).exp(),
            1.0e-12,
        );
    }
}

#[test]
fn native_curved_water_ray_preserves_tangential_optical_momentum_and_beer_length() {
    let state = Polynomial {
        gt: [1.0, 0.0, 0.0],
        gr: [20.0, 0.0, 0.0],
        ..Polynomial::uniform()
    };
    let absorption = BeerLambertAbsorption::try_constant(0.2).unwrap();
    let field =
        MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([76; 32])).unwrap();
    for lambda in [300.0, 589.0, 1000.0] {
        let prepared = field.prepare(lambda).unwrap();
        let initial = prepared.optical_sample(ray().origin).unwrap();
        let result = propagate(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            lambda,
            1.5,
            options(),
        )
        .unwrap();
        let final_index = prepared.optical_sample(result.ray.origin).unwrap();
        // Translational symmetry in y,z conserves n*u_y and n*u_z. This
        // independently checks bending from both original state gradients.
        assert!(result.ray.origin.x > 0.003 && result.ray.dir.x > 0.004);
        near(
            final_index.phase_index() * result.ray.dir.z,
            initial.phase_index(),
            1.0e-10,
        );
        near(result.ray.dir.y, 0.0, 1.0e-15);
        near(result.power_transmittance(), (-0.3_f64).exp(), 1.0e-12);
        assert!(result.group_optical_path_length_m.unwrap() > result.optical_path_length_m);
        assert_eq!(
            result,
            propagate(
                &field,
                ray(),
                Vec3::new(1.0, 0.0, 0.0),
                lambda,
                1.5,
                options()
            )
            .unwrap()
        );
    }
}

#[test]
fn original_snapshot_domains_absorption_support_and_prepared_wavelength_are_enforced() {
    let state = Polynomial::uniform();
    let absorption =
        BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[(400.0, 0.1), (900.0, 0.3)])
            .unwrap();
    let field =
        MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([77; 32])).unwrap();
    assert_eq!(field.wavelength_range_nm(), [400.0, 900.0]);
    assert!(core::ptr::eq(field.state_source(), &state));
    assert!(core::ptr::eq(field.absorption_source(), &absorption));
    assert_eq!(field.snapshot_identity(), ContentHash([77; 32]));
    assert!(matches!(
        MaterialWaterRefractiveField::try_new(&state, &absorption, ContentHash([0; 32])),
        Err(WaterRefractiveFieldError::MissingSnapshotIdentity)
    ));
    assert!(field.prepare(399.0).is_err());
    let prepared = field.prepare(600.0).unwrap();
    assert!(prepared.sample(ray().origin, 601.0).is_err());
    assert!(prepared.optical_sample(Point3::new(5.0, 0.0, 0.0)).is_err());
    assert!(
        prepared
            .optical_sample(Point3::new(f64::NAN, 0.0, 0.0))
            .is_err()
    );
    near(
        prepared
            .optical_sample(ray().origin)
            .unwrap()
            .absorption_per_m(),
        0.18,
        1.0e-15,
    );
    let invalid = Polynomial { r: 1070.0, ..state };
    let invalid =
        MaterialWaterRefractiveField::try_new(&invalid, &absorption, ContentHash([78; 32]))
            .unwrap();
    assert_eq!(
        invalid.sample(ray().origin, 600.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    let other = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (1500.0, 0.1),
        (2000.0, 0.1),
    ])
    .unwrap();
    assert!(matches!(
        MaterialWaterRefractiveField::try_new(&state, &other, ContentHash([79; 32])),
        Err(WaterRefractiveFieldError::DisjointSpectralSupport)
    ));
}
