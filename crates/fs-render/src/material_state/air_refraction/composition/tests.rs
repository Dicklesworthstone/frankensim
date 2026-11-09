use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};
use crate::material_state::thermo_optic::{TemperatureGradientSample, TemperatureJet};
use crate::material_state::water_refraction::pressure_temperature::PressureTemperatureJet;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn point(p: [f64; 3]) -> Point3 {
    Point3::new(p[0], p[1], p[2])
}
fn options() -> GradientIndexOptions {
    GradientIndexOptions {
        initial_step_m: 0.5,
        max_step_m: 2.0,
        position_tolerance_m: 1e-11,
        direction_tolerance: 1e-12,
        optical_path_tolerance_m: 1e-11,
        ..GradientIndexOptions::default()
    }
}

struct State(bool);
impl PressureTemperatureGradientField for State {
    fn pressure_temperature_gradient(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureGradientSample, GradientIndexError> {
        Ok(self.pressure_temperature_jet(p)?.sample())
    }
    fn is_spatially_uniform(&self) -> bool {
        !self.0
    }
}
impl PressureTemperatureSecondOrderField for State {
    fn pressure_temperature_jet(
        &self,
        p: Point3,
    ) -> Result<PressureTemperatureJet, GradientIndexError> {
        if ![p.x, p.y, p.z].into_iter().all(f64::is_finite) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let s = if self.0 { 1.0 } else { 0.0 };
        Ok(PressureTemperatureJet {
            temperature: TemperatureJet {
                sample: TemperatureGradientSample {
                    temperature_k: 293.15
                        + s * (3.0 * p.x
                            + 2.0 * p.y
                            + 0.5 * p.x * p.x
                            + 0.4 * p.x * p.y
                            + 0.3 * p.z * p.z),
                    gradient_k_per_m: Vec3::new(
                        s * (3.0 + p.x + 0.4 * p.y),
                        s * (2.0 + 0.4 * p.x),
                        s * 0.6 * p.z,
                    ),
                },
                diagonal_k_per_m2: [s, 0.0, 0.6 * s],
                off_diagonal_k_per_m2: [0.4 * s, 0.0, 0.0],
            },
            pressure_pa: 101_325.0
                + s * (1000.0 * p.x - 700.0 * p.y
                    + 500.0 * p.z
                    + 100.0 * p.x * p.y
                    + 50.0 * p.y * p.z),
            pressure_gradient_pa_per_m: Vec3::new(
                s * (1000.0 + 100.0 * p.y),
                s * (-700.0 + 100.0 * p.x + 50.0 * p.z),
                s * (500.0 + 50.0 * p.y),
            ),
            pressure_diagonal_pa_per_m2: [0.0; 3],
            pressure_off_diagonal_pa_per_m2: [100.0 * s, 0.0, 50.0 * s],
        })
    }
}

struct MixedComposition;
impl AirCompositionGradientField for MixedComposition {
    fn air_composition_gradient(
        &self,
        p: Point3,
    ) -> Result<AirCompositionGradientSample, GradientIndexError> {
        Ok(self.air_composition_jet(p)?.sample)
    }
}
impl AirCompositionSecondOrderField for MixedComposition {
    fn air_composition_jet(&self, p: Point3) -> Result<AirCompositionJet, GradientIndexError> {
        let composition = CiddorAir::try_new(
            0.02 + 0.005 * p.x - 0.003 * p.y + 0.001 * p.z + 0.001 * p.x * p.x + 0.002 * p.x * p.y,
            600.0 + 100.0 * p.x + 200.0 * p.y - 50.0 * p.z + 20.0 * p.y * p.z + 30.0 * p.z * p.z,
        )
        .map_err(|_| GradientIndexError::OutsideFieldDomain)?;
        Ok(AirCompositionJet {
            sample: AirCompositionGradientSample {
                composition,
                water_gradient_per_m: Vec3::new(
                    0.005 + 0.002 * p.x + 0.002 * p.y,
                    -0.003 + 0.002 * p.x,
                    0.001,
                ),
                co2_gradient_micromol_per_mol_m: Vec3::new(
                    100.0,
                    200.0 + 20.0 * p.z,
                    -50.0 + 20.0 * p.y + 60.0 * p.z,
                ),
            },
            water_diagonal_per_m2: [0.002, 0.0, 0.0],
            water_off_diagonal_per_m2: [0.002, 0.0, 0.0],
            co2_diagonal_micromol_per_mol_m2: [0.0, 0.0, 60.0],
            co2_off_diagonal_micromol_per_mol_m2: [0.0, 0.0, 20.0],
        })
    }
}

#[test]
fn g1_full_original_state_chain_rule_retains_all_spatial_mixed_responses() {
    let state = State(true);
    let field = MaterialAirCompositionRefractiveField::try_nonabsorbing(
        &state,
        &MixedComposition,
        ContentHash([151; 32]),
    )
    .unwrap();
    let p = [0.1, -0.2, 0.3];
    for wavelength in [300.0, 633.0, 1700.0] {
        let prepared = field.prepare(wavelength).unwrap();
        let jet = prepared.optical_jet(point(p)).unwrap();
        assert_eq!(jet.optical(), prepared.optical_sample(point(p)).unwrap());
        assert_eq!(jet.optical(), field.sample(point(p), wavelength).unwrap());
        // Independent spatial differences query the original frozen-composition
        // phase law, never the full-state derivative implementation under test.
        let value = |p| {
            let state = state.pressure_temperature_gradient(point(p)).unwrap();
            MixedComposition
                .air_composition_gradient(point(p))
                .unwrap()
                .composition
                .prepare(wavelength)
                .unwrap()
                .sample(state.temperature.temperature_k, state.pressure_pa)
                .unwrap()
                .phase_index_minus_one
        };
        let (h, center) = (0.001, value(p));
        let g = jet.optical().phase_gradient_per_m();
        for (i, gradient) in [g.x, g.y, g.z].into_iter().enumerate() {
            let (mut plus, mut minus) = (p, p);
            plus[i] += h;
            minus[i] -= h;
            close(gradient, (value(plus) - value(minus)) / (2.0 * h), 3e-13);
            close(
                jet.hessian().diagonal()[i],
                (value(plus) - 2.0 * center + value(minus)) / (h * h),
                3e-13,
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
            close(
                jet.hessian().off_diagonal()[k],
                (value(pp) - value(pm) - value(mp) + value(mm)) / (4.0 * h * h),
                3e-13,
            );
        }
    }
}

struct CompositionLens {
    water_curvature: f64,
    co2_curvature: f64,
}
impl AirCompositionGradientField for CompositionLens {
    fn air_composition_gradient(
        &self,
        p: Point3,
    ) -> Result<AirCompositionGradientSample, GradientIndexError> {
        if p.x.abs() > 0.1 || p.y.abs() > 0.1 {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let r2 = p.x * p.x + p.y * p.y;
        let composition = CiddorAir::try_new(
            0.02 + 0.5 * self.water_curvature * r2,
            600.0 + 0.5 * self.co2_curvature * r2,
        )
        .map_err(|_| GradientIndexError::OutsideFieldDomain)?;
        Ok(AirCompositionGradientSample {
            composition,
            water_gradient_per_m: Vec3::new(
                self.water_curvature * p.x,
                self.water_curvature * p.y,
                0.0,
            ),
            co2_gradient_micromol_per_mol_m: Vec3::new(
                self.co2_curvature * p.x,
                self.co2_curvature * p.y,
                0.0,
            ),
        })
    }
}
impl AirCompositionSecondOrderField for CompositionLens {
    fn air_composition_jet(&self, p: Point3) -> Result<AirCompositionJet, GradientIndexError> {
        Ok(AirCompositionJet {
            sample: self.air_composition_gradient(p)?,
            water_diagonal_per_m2: [self.water_curvature, self.water_curvature, 0.0],
            water_off_diagonal_per_m2: [0.0; 3],
            co2_diagonal_micromol_per_mol_m2: [self.co2_curvature, self.co2_curvature, 0.0],
            co2_off_diagonal_micromol_per_mol_m2: [0.0; 3],
        })
    }
}

#[test]
fn g1_water_and_co2_each_bend_and_focus_with_constant_temperature_and_pressure() {
    let state = State(false);
    for composition in [
        CompositionLens {
            water_curvature: 1.0,
            co2_curvature: 0.0,
        },
        CompositionLens {
            water_curvature: 0.0,
            co2_curvature: -20_000.0,
        },
    ] {
        let field = MaterialAirCompositionRefractiveField::try_nonabsorbing(
            &state,
            &composition,
            ContentHash([152; 32]),
        )
        .unwrap();
        assert!(!field.is_spatially_uniform());
        for wavelength in [400.0, 633.0, 1550.0] {
            let prepared = field.prepare(wavelength).unwrap();
            let origin = point([0.0; 3]);
            let local = prepared.optical_jet(origin).unwrap();
            assert_eq!(
                local.optical().phase_gradient_per_m(),
                Vec3::new(0.0, 0.0, 0.0)
            );
            let omega = (-local.hessian().diagonal()[0] / local.optical().phase_index()).sqrt();
            assert!(omega > 0.0);
            let length = 80.0;
            let ray = Ray {
                origin,
                dir: Vec3::new(0.0, 0.0, 1.0),
            };
            let bundle = propagate_bundle(
                &prepared,
                ray,
                Vec3::new(1.0, 0.0, 0.0),
                wavelength,
                length,
                BundleOptions {
                    integration: options(),
                    position_derivative_tolerance_m: 1e-10,
                    direction_derivative_tolerance: 1e-12,
                },
            )
            .unwrap();
            let footprint = bundle
                .project_at_endpoint(Vec3::new(0.0, 0.0, 1.0), ProjectionOptions::default())
                .unwrap();
            let expected = (omega * length).sin() / omega;
            close(
                footprint.area_per_solid_angle_m2(),
                expected * expected,
                1e-6,
            );
            assert!(footprint.area_per_solid_angle_m2() < length * length - 0.1);
            close(
                bundle.propagation().optical_path_length_m,
                length * local.optical().phase_index(),
                5e-11,
            );
            close(
                bundle.propagation().group_optical_path_length_m.unwrap(),
                length * local.optical().group_index().unwrap(),
                5e-11,
            );
            let off_axis = Ray {
                origin: point([0.02, 0.0, 0.0]),
                ..ray
            };
            let initial = prepared.optical_sample(off_axis.origin).unwrap();
            let bent = propagate(
                &prepared,
                off_axis,
                Vec3::new(1.0, 0.0, 0.0),
                wavelength,
                20.0,
                options(),
            )
            .unwrap();
            assert!(bent.ray.dir.x < -1e-7);
            close(
                prepared
                    .optical_sample(bent.ray.origin)
                    .unwrap()
                    .phase_index()
                    * bent.ray.dir.z,
                initial.phase_index(),
                2e-12,
            );
            assert_eq!(bent.power_transmittance(), 1.0);
        }
    }
}

#[test]
fn g0_composition_support_snapshot_and_prepared_wavelength_are_retained() {
    let state = State(false);
    let composition = CompositionLens {
        water_curvature: 1.0,
        co2_curvature: 0.0,
    };
    let bind = |id| {
        MaterialAirCompositionRefractiveField::try_nonabsorbing(
            &state,
            &composition,
            ContentHash([id; 32]),
        )
    };
    assert!(matches!(
        bind(0),
        Err(AirFieldError::MissingSnapshotIdentity)
    ));
    let field = bind(153).unwrap();
    assert_ne!(
        field.content_identity(),
        bind(154).unwrap().content_identity()
    );
    assert_eq!(field.model_identity(), ciddor_model_identity());
    assert!(core::ptr::eq(field.state_source(), &state));
    assert!(core::ptr::eq(field.composition_source(), &composition));
    for p in [point([0.2, 0.0, 0.0]), point([f64::NAN, 0.0, 0.0])] {
        assert!(field.sample(p, 550.0).is_err());
        assert!(field.sample_second_order(p, 550.0).is_err());
    }
    assert!(field.prepare(299.9).is_err());
    let prepared = field.prepare(550.0).unwrap();
    assert!(prepared.sample(point([0.0; 3]), 550.1).is_err());
    assert!(
        prepared
            .sample_second_order(point([0.0; 3]), 550.1)
            .is_err()
    );
}
