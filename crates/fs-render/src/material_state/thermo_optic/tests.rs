use super::*;
use crate::charts::Ray;
use crate::dielectric::{
    BeerLambertAbsorption, DielectricIorSample, GlassProvenance, TabulatedIor,
};
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate, propagate_with_control};
use fs_material::thermo_optic::ThermoOpticCoefficients;

struct PolynomialTemperature {
    centre_k: f64,
    gradient: Vec3,
    diagonal: [f64; 3],
    off_diagonal: [f64; 3],
}

impl PolynomialTemperature {
    fn jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        if ![p.x, p.y, p.z]
            .into_iter()
            .all(|coordinate| coordinate.is_finite() && coordinate.abs() <= 2.0)
        {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let [xx, yy, zz] = self.diagonal;
        let [xy, xz, yz] = self.off_diagonal;
        let h = Vec3::new(
            xx * p.x + xy * p.y + xz * p.z,
            xy * p.x + yy * p.y + yz * p.z,
            xz * p.x + yz * p.y + zz * p.z,
        );
        let displacement = Vec3::new(p.x, p.y, p.z);
        Ok(TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: self.centre_k
                    + self.gradient.dot(displacement)
                    + 0.5 * h.dot(displacement),
                gradient_k_per_m: Vec3::new(
                    self.gradient.x + h.x,
                    self.gradient.y + h.y,
                    self.gradient.z + h.z,
                ),
            },
            diagonal_k_per_m2: self.diagonal,
            off_diagonal_k_per_m2: self.off_diagonal,
        })
    }
}
impl TemperatureField for PolynomialTemperature {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        self.jet(p).map(|jet| jet.sample.temperature_k)
    }
}
impl TemperatureGradientField for PolynomialTemperature {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        self.jet(p).map(|jet| jet.sample)
    }
    fn is_spatially_uniform(&self) -> bool {
        // These bounded fixtures deliberately cannot use a uniform shortcut.
        false
    }
}
impl TemperatureSecondOrderField for PolynomialTemperature {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        self.jet(p)
    }
}

fn temperature(gradient: Vec3, diagonal: [f64; 3]) -> PolynomialTemperature {
    PolynomialTemperature {
        centre_k: 300.0,
        gradient,
        diagonal,
        off_diagonal: [0.0; 3],
    }
}
fn law(d: [f64; 3], e: [f64; 2]) -> SchottThermoOpticLaw {
    // Synthetic coefficient magnitudes make numerical effects observable;
    // this fixture does not assert a measured material or fit accuracy.
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d,
            e,
            resonance_um: 0.0,
        },
        300.0,
        [200.0, 450.0],
        [1000.0, 2000.0],
        ContentHash([1; 32]),
    )
    .unwrap()
}
fn glass(absorption: BeerLambertAbsorption) -> DielectricGlass {
    DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(1000.0 + i as f64 * 125.0, 1.5).unwrap()
            }),
            125.0,
        )
        .unwrap(),
        absorption,
        GlassProvenance::Custom,
    )
}
fn origin() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}
fn ray() -> Ray {
    Ray {
        origin: origin(),
        dir: Vec3::new(0.0, 0.0, 1.0),
    }
}
fn tight() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.01,
        max_step_m: 0.04,
        relative_tolerance: 2e-12,
        position_tolerance_m: 2e-12,
        direction_tolerance: 2e-12,
        optical_path_tolerance_m: 2e-12,
        optical_depth_tolerance: 2e-12,
        ..GradientIndexOptions::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}; expected={expected:.17e}; tolerance={tolerance:.3e}"
    );
}

#[test]
fn g1_original_temperature_drives_native_infrared_bending_group_delay_and_absorption() {
    let glass = glass(BeerLambertAbsorption::try_constant(0.4).unwrap());
    let law = law([0.012, 0.0, 0.0], [0.0045, 0.0]);
    let temperature = temperature(Vec3::new(40.0, 0.0, 0.0), [0.0; 3]);
    let field = MaterialThermoOpticField::try_new(&glass, &law, &temperature, ContentHash([2; 32]))
        .unwrap();
    let prepared = field.prepare(1500.0).unwrap();
    let length = 0.8;
    let path = propagate(
        &prepared,
        ray(),
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        tight(),
    )
    .unwrap();
    // At 1.5um, A=5/12, F=.014*(T-300), and group F=.018*(T-300).
    // Therefore the original linear thermal profile gives n=1.5+(7/30)*x.
    let gradient = 7.0 / 30.0;
    let n_end = (2.25 + (gradient * length).powi(2)).sqrt();
    let z = (1.5 / gradient) * (gradient * length / 1.5).asinh();
    let phase = 0.5 * (length * n_end + 1.5 * z);
    near(path.ray.origin.x, (n_end - 1.5) / gradient, 2e-11);
    near(path.ray.origin.z, z, 2e-11);
    near(path.ray.dir.z, 1.5 / n_end, 2e-11);
    near(path.optical_path_length_m, phase, 2e-11);
    near(
        path.group_optical_path_length_m.unwrap(),
        1.5 * length + (9.0 / 7.0) * (phase - 1.5 * length),
        2e-11,
    );
    near(path.power_transmittance(), (-0.4 * length).exp(), 2e-13);
    assert_eq!(
        prepared.temperature_k(path.ray.origin),
        temperature.temperature_k(path.ray.origin)
    );
    assert_eq!(prepared.snapshot_identity(), ContentHash([2; 32]));
    assert_eq!(
        field.sample(path.ray.origin, 1500.0).unwrap(),
        prepared.sample(path.ray.origin, 1500.0).unwrap(),
    );
}

#[test]
fn g1_thermal_hessian_focuses_native_ray_bundles_at_zero_local_gradient() {
    let glass = glass(BeerLambertAbsorption::try_constant(0.25).unwrap());
    let law = law([0.12, 0.0, 0.0], [0.0; 2]);
    let mut temperature = temperature(Vec3::new(0.0, 0.0, 0.0), [-2.0, -3.0, 0.0]);
    temperature.centre_k = 302.0;
    let field = MaterialThermoOpticField::try_new(&glass, &law, &temperature, ContentHash([3; 32]))
        .unwrap();
    let prepared = field.prepare(1500.0).unwrap();
    let length = 0.8;
    let bundle = propagate_bundle(
        &prepared,
        ray(),
        Vec3::new(1.0, 0.0, 0.0),
        1500.0,
        length,
        BundleOptions {
            integration: tight(),
            ..BundleOptions::default()
        },
    )
    .unwrap();
    let projection = bundle
        .project_at_endpoint(ray().dir, ProjectionOptions::default())
        .unwrap();
    let omega_x = (0.1_f64 / 1.6).sqrt();
    let omega_y = (0.15_f64 / 1.6).sqrt();
    let area = ((omega_x * length).sin() / omega_x) * ((omega_y * length).sin() / omega_y);
    near(projection.area_per_solid_angle_m2(), area, 2e-11);
    near(
        bundle.propagation().optical_path_length_m,
        1.6 * length,
        2e-12,
    );
    near(
        bundle.propagation().power_transmittance(),
        (-0.25 * length).exp(),
        2e-13,
    );
    assert_eq!(
        prepared
            .optical_sample(origin())
            .unwrap()
            .phase_gradient_per_m(),
        Vec3::new(0.0, 0.0, 0.0)
    );
    assert_ne!(
        prepared.optical_jet(origin()).unwrap().hessian(),
        IndexHessian::ZERO
    );
}

#[test]
fn g1_nonlinear_thermal_law_retains_all_hessian_chain_rule_terms() {
    let glass = glass(BeerLambertAbsorption::CLEAR);
    let law = law([0.006, 0.0002, 0.000003], [0.0; 2]);
    let mut temperature = temperature(Vec3::new(2.0, -3.0, 4.0), [1.0, -2.0, 3.0]);
    temperature.centre_k = 310.0;
    temperature.off_diagonal = [0.5, -0.75, 1.25];
    let field = MaterialThermoOpticField::try_new(&glass, &law, &temperature, ContentHash([4; 32]))
        .unwrap();
    let prepared = field.prepare(1500.0).unwrap();
    let jet = prepared.optical_jet(origin()).unwrap();
    let first = (5.0 / 12.0) * (0.006 + 2.0 * 0.0002 * 10.0 + 3.0 * 0.000003 * 100.0);
    let second = (5.0 / 12.0) * (2.0 * 0.0002 + 6.0 * 0.000003 * 10.0);
    near(
        jet.optical().phase_index(),
        1.5 + (5.0 / 12.0) * (0.06 + 0.02 + 0.003),
        2e-15,
    );
    for (actual, expected) in jet.hessian().diagonal().into_iter().zip([
        first + second * 4.0,
        -2.0 * first + second * 9.0,
        3.0 * first + second * 16.0,
    ]) {
        near(actual, expected, 2e-16);
    }
    for (actual, expected) in jet.hessian().off_diagonal().into_iter().zip([
        0.5 * first - 6.0 * second,
        -0.75 * first + 8.0 * second,
        1.25 * first - 12.0 * second,
    ]) {
        near(actual, expected, 2e-16);
    }
    assert_eq!(jet.optical(), prepared.optical_sample(origin()).unwrap());
}

#[test]
fn g0_common_source_support_and_unknown_group_remain_explicit() {
    let absorption = BeerLambertAbsorption::try_from_spectral_extinction_in_band(&[
        (1400.0, 0.2),
        (1800.0, 0.6),
    ])
    .unwrap();
    let glass = DielectricGlass::new(
        TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(
                    1000.0 + 125.0 * i as f64,
                    if i == 4 { 1.55 } else { 1.5 },
                )
                .unwrap()
            }),
            125.0,
        )
        .unwrap(),
        absorption,
        GlassProvenance::Custom,
    );
    let law = law([0.0001, 0.0, 0.0], [0.0; 2]);
    let temperature = temperature(Vec3::new(0.0, 0.0, 0.0), [0.0; 3]);
    let field = MaterialThermoOpticField::try_new(&glass, &law, &temperature, ContentHash([5; 32]))
        .unwrap();
    assert_eq!(field.wavelength_range_nm(), [1400.0, 1800.0]);
    let prepared = field.prepare(1500.0).unwrap();
    let sample = prepared.optical_sample(origin()).unwrap();
    assert_eq!(sample.group_index(), None);
    assert_eq!(sample.phase_index(), 1.55);
    near(sample.absorption_per_m(), 0.3, 2e-16);
    for wavelength in [1399.0, 1801.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            field.prepare(wavelength),
            Err(ThermoOpticFieldError::Optics(
                DielectricError::InvalidWavelength
            ))
        ));
    }
    assert_eq!(
        prepared.sample(origin(), 1501.0),
        Err(DielectricError::InvalidWavelength.into())
    );
    assert!(matches!(
        MaterialThermoOpticField::try_new(&glass, &law, &temperature, ContentHash([0; 32])),
        Err(ThermoOpticFieldError::MissingSnapshotIdentity)
    ));
    let disjoint = SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [0.0; 3],
            e: [0.0; 2],
            resonance_um: 0.0,
        },
        300.0,
        [200.0, 400.0],
        [2100.0, 2500.0],
        ContentHash([6; 32]),
    )
    .unwrap();
    assert!(matches!(
        MaterialThermoOpticField::try_new(&glass, &disjoint, &temperature, ContentHash([5; 32])),
        Err(ThermoOpticFieldError::DisjointSpectralSupport)
    ));
}

#[test]
fn g0_thermal_domains_and_native_cancellation_do_not_extrapolate_or_mutate() {
    let glass = glass(BeerLambertAbsorption::CLEAR);
    let law = law([0.0001, 0.0, 0.0], [0.0; 2]);
    let temperature = temperature(Vec3::new(100.0, 0.0, 0.0), [0.0; 3]);
    let field = MaterialThermoOpticField::try_new(&glass, &law, &temperature, ContentHash([7; 32]))
        .unwrap();
    let prepared = field.prepare(1500.0).unwrap();
    for p in [
        Point3::new(1.6, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
        Point3::new(f64::NAN, 0.0, 0.0),
    ] {
        assert_eq!(
            prepared.sample(p, 1500.0),
            Err(GradientIndexError::OutsideFieldDomain)
        );
    }
    let before = prepared.sample(origin(), 1500.0).unwrap();
    assert_eq!(
        propagate_with_control(
            &prepared,
            ray(),
            Vec3::new(1.0, 0.0, 0.0),
            1500.0,
            0.5,
            tight(),
            || false
        ),
        Err(GradientIndexError::Cancelled)
    );
    assert_eq!(prepared.sample(origin(), 1500.0).unwrap(), before);
    assert_eq!(field.temperature_source().centre_k, 300.0);
}
