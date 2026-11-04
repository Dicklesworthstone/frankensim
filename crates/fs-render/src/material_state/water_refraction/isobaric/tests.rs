use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexField, GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::conduction::PlanarConductionTemperatureField;
use fs_geom::Vec3;
use fs_material::water_refraction::IapwsWaterRefraction;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:.3e}"
    );
}
fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}
#[derive(Clone, Copy)]
struct Polynomial {
    t: f64,
    g: [f64; 3],
    h: [f64; 6],
}
impl Polynomial {
    fn sample(self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        let x = [p.x, p.y, p.z];
        if !x.into_iter().all(|v| v.is_finite() && v.abs() <= 4.0) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let (mut t, mut g) = (self.t, self.g);
        for i in 0..3 {
            t += self.g[i] * x[i] + 0.5 * self.h[i] * x[i] * x[i];
            g[i] += self.h[i] * x[i];
        }
        for ((i, j), h) in [(0, 1), (0, 2), (1, 2)].into_iter().zip(&self.h[3..]) {
            t += h * x[i] * x[j];
            g[i] += h * x[j];
            g[j] += h * x[i];
        }
        Ok(TemperatureJet {
            sample: TemperatureGradientSample {
                temperature_k: t,
                gradient_k_per_m: Vec3::new(g[0], g[1], g[2]),
            },
            diagonal_k_per_m2: self.h[..3].try_into().unwrap(),
            off_diagonal_k_per_m2: self.h[3..].try_into().unwrap(),
        })
    }
}
impl TemperatureField for Polynomial {
    fn temperature_k(&self, p: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.sample(p)?.sample.temperature_k)
    }
}
impl TemperatureGradientField for Polynomial {
    fn temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        Ok(self.sample(p)?.sample)
    }
}
impl TemperatureSecondOrderField for Polynomial {
    fn temperature_jet(&self, p: Point3) -> Result<TemperatureJet, GradientIndexError> {
        self.sample(p)
    }
}
fn uniform(t: f64) -> Polynomial {
    Polynomial {
        t,
        g: [0.0; 3],
        h: [0.0; 6],
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
fn g1_spatial_water_jet_differentiates_both_original_eos_and_optical_sources() {
    let thermal = Polynomial {
        g: [5.0, -3.0, 2.0],
        h: [8.0, -7.0, 2.0, 3.0, -4.0, 1.0],
        ..uniform(300.0)
    };
    let field = IsobaricWaterField::try_new(
        &thermal,
        Pressure::new(3e6),
        If97Region::Region1,
        ContentHash([81; 32]),
    )
    .unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let optics = field.bind_refraction(&absorption).unwrap();
    assert!(core::ptr::eq(field.temperature_source(), &thermal));
    assert!(core::ptr::eq(optics.state_source(), &field));
    let p = [0.02, -0.04, 0.03];
    assert_eq!(
        field.sample_gradient(point(p)).unwrap(),
        field.sample_jet(point(p)).unwrap().sample()
    );
    for lambda in [226.5, 589.0, 1013.98] {
        let original = IapwsWaterRefraction::prepare(lambda).unwrap();
        // Reference uses only the original scalar temperature, native IF97
        // state value, and original optical source; no adapter derivatives.
        let phase = |x| {
            let t = thermal.temperature_k(point(x)).unwrap();
            let s = If97State::at_pressure_temperature(
                Pressure::new(3e6),
                Temperature::new(t),
                If97Region::Region1,
            )
            .unwrap();
            original.sample_if97_state(&s).unwrap().phase_index
        };
        let jet = optics
            .prepare(lambda)
            .unwrap()
            .optical_jet(point(p))
            .unwrap();
        let g = jet.optical().phase_gradient_per_m();
        let gradient = [g.x, g.y, g.z];
        let n = phase(p);
        let h = 0.002;
        for i in 0..3 {
            let (mut a, mut b) = (p, p);
            a[i] += h;
            b[i] -= h;
            near(gradient[i], (phase(a) - phase(b)) / (2.0 * h), 2e-9);
            near(
                jet.hessian().diagonal()[i],
                (phase(a) - 2.0 * n + phase(b)) / (h * h),
                3e-9,
            );
        }
        for (k, (i, j)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
            let (mut pp, mut pm, mut mp, mut mm) = (p, p, p, p);
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
                3e-9,
            );
        }
    }
}

#[test]
fn g1_isobaric_thermal_expansion_drives_analytic_native_bundle_focusing() {
    let curvature = 800.0;
    let thermal = Polynomial {
        h: [curvature, curvature, 0.0, 0.0, 0.0, 0.0],
        ..uniform(333.15)
    };
    let field = IsobaricWaterField::try_new(
        &thermal,
        Pressure::new(1e5),
        If97Region::Region1,
        ContentHash([82; 32]),
    )
    .unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.02).unwrap();
    let optics = field.bind_refraction(&absorption).unwrap();
    let state = field.sample_state(point([0.0; 3])).unwrap();
    for lambda in [400.0, 589.0, 900.0] {
        let original = IapwsWaterRefraction::prepare(lambda)
            .unwrap()
            .sample_if97_state(&state)
            .unwrap();
        let rho_t = -state.density().value() * state.thermal_expansion().value();
        let total_n_t =
            original.temperature_derivative_per_k + original.density_derivative_m3_per_kg * rho_t;
        assert!(total_n_t < 0.0);
        assert!((total_n_t - original.temperature_derivative_per_k).abs() > total_n_t.abs() * 0.5);
        let omega = (-curvature * total_n_t / original.phase_index).sqrt();
        let length = 2.0;
        let prepared = optics.prepare(lambda).unwrap();
        let bundle = propagate_bundle(
            &prepared,
            Ray {
                origin: point([0.0; 3]),
                dir: Vec3::new(0.0, 0.0, 1.0),
            },
            Vec3::new(1.0, 0.0, 0.0),
            lambda,
            length,
            BundleOptions {
                integration: options(),
                position_derivative_tolerance_m: 1e-12,
                direction_derivative_tolerance: 1e-12,
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
            2e-9,
        );
        assert!(footprint.area_per_solid_angle_m2() < length * length - 0.1);
        near(
            bundle.propagation().optical_path_length_m,
            length * original.phase_index,
            1e-11,
        );
        near(
            bundle.propagation().group_optical_path_length_m.unwrap(),
            length * original.group_index,
            1e-11,
        );
        near(
            bundle.propagation().power_transmittance(),
            (-0.02 * length).exp(),
            1e-12,
        );
    }
}

#[test]
fn g1_native_conduction_solution_drives_water_state_and_curved_transport() {
    use fs_conduction::ConductivityTable;
    use fs_conduction::layered::{LayeredWall, PlanarThermalLayer};
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    let wall = LayeredWall::new(
        vec![
            PlanarThermalLayer::new(
                "declared thermal specimen",
                2.0,
                ConductivityTable::declared(4.0).unwrap(),
            )
            .unwrap(),
        ],
        vec![],
        [0.0; 2],
    )
    .unwrap();
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let solution = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 83,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        wall.solve_with_generation(&cx, [300.0; 2], 1.0, &[8.0])
            .unwrap()
    });
    let thermal = PlanarConductionTemperatureField::try_new(
        &solution,
        0,
        point([0.0; 3]),
        Vec3::new(1.0, 0.0, 0.0),
    )
    .unwrap();
    let field = IsobaricWaterField::try_new(
        &thermal,
        Pressure::new(3e6),
        If97Region::Region1,
        ContentHash([83; 32]),
    )
    .unwrap();
    let absorption = BeerLambertAbsorption::try_constant(0.03).unwrap();
    let optics = field.bind_refraction(&absorption).unwrap();
    let prepared = optics.prepare(589.0).unwrap();
    let start = point([0.4, 0.0, 0.0]);
    let state = field.sample_state(start).unwrap();
    near(
        state.temperature().value(),
        300.0 + 2.0 * 0.4 - 0.4 * 0.4,
        1e-13,
    );
    assert!(core::ptr::eq(thermal.solution(), &solution));
    let initial = prepared.sample(start, 589.0).unwrap();
    let path = propagate(
        &prepared,
        Ray {
            origin: start,
            dir: Vec3::new(0.0, 0.0, 1.0),
        },
        Vec3::new(1.0, 0.0, 0.0),
        589.0,
        0.5,
        options(),
    )
    .unwrap();
    let final_sample = prepared.sample(path.ray.origin, 589.0).unwrap();
    assert!(path.ray.dir.x < -1e-5);
    near(
        final_sample.phase_index() * path.ray.dir.z,
        initial.phase_index(),
        2e-11,
    );
    near(path.power_transmittance(), (-0.03_f64 * 0.5).exp(), 1e-12);
    assert!(matches!(
        field.water_state_gradient(point([-0.01, 0.0, 0.0])),
        Err(GradientIndexError::OutsideFieldDomain)
    ));
}

#[test]
fn g3_source_binding_and_phase_domains_remain_explicit() {
    let thermal = uniform(300.0);
    let make = |p, region, id| {
        IsobaricWaterField::try_new(&thermal, Pressure::new(p), region, ContentHash([id; 32]))
    };
    let field = make(1e5, If97Region::Region1, 84).unwrap();
    assert!(!field.is_spatially_uniform());
    assert_ne!(
        field.content_identity(),
        make(2e5, If97Region::Region1, 84)
            .unwrap()
            .content_identity()
    );
    assert_ne!(
        field.content_identity(),
        make(1e5, If97Region::Region1, 85)
            .unwrap()
            .content_identity()
    );
    assert_ne!(
        field.content_identity(),
        make(1e5, If97Region::Region2, 84)
            .unwrap()
            .content_identity()
    );
    assert!(
        make(1e5, If97Region::Region2, 84)
            .unwrap()
            .sample_state(point([0.0; 3]))
            .is_err()
    );
    for p in [0.0, -1.0, f64::NAN, 1.1e8] {
        assert!(matches!(
            make(p, If97Region::Region1, 84),
            Err(IsobaricWaterError::InvalidPressure)
        ));
    }
    assert!(matches!(
        make(1e5, If97Region::Region1, 0),
        Err(IsobaricWaterError::MissingSnapshotIdentity)
    ));
    assert!(matches!(
        field.water_state_gradient(point([5.0, 0.0, 0.0])),
        Err(GradientIndexError::OutsideFieldDomain)
    ));
    let boiling = Polynomial {
        g: [100.0, 0.0, 0.0],
        ..uniform(350.0)
    };
    let liquid = IsobaricWaterField::try_new(
        &boiling,
        Pressure::new(1e5),
        If97Region::Region1,
        ContentHash([86; 32]),
    )
    .unwrap();
    assert!(liquid.sample_state(point([0.0; 3])).is_ok());
    assert!(matches!(
        liquid.water_state_gradient(point([0.3, 0.0, 0.0])),
        Err(GradientIndexError::OutsideFieldDomain)
    ));
    let hot = uniform(800.0);
    let vapor = IsobaricWaterField::try_new(
        &hot,
        Pressure::new(1e5),
        If97Region::Region2,
        ContentHash([87; 32]),
    )
    .unwrap();
    assert!(vapor.sample_state(point([0.0; 3])).is_ok());
    let absorption = BeerLambertAbsorption::try_constant(0.0).unwrap();
    let optics = vapor.bind_refraction(&absorption).unwrap();
    assert!(optics.sample(point([0.0; 3]), 589.0).is_err());
    assert!(optics.prepare(1101.0).is_err());
}
