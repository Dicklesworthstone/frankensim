use super::*;
use crate::charts::Ray;
use crate::gradient_index::bundle::{BundleOptions, ProjectionOptions, propagate_bundle};
use crate::gradient_index::{GradientIndexOptions, propagate};
use fs_material::gas::GasSpec;
use fs_material::gas_scattering::{MolecularRayleigh, MolecularRayleighMixture};
use fs_qty::semantic::{Composition, CompositionBasis};

fn binding() -> MaterialGasScatteringBinding {
    let model = MolecularRayleighMixture::try_new(
        &[GasSpec::dry_air_ussa1976()],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        &[MolecularRayleigh::try_new(0.0003, 2.5e25, 550.0, 0.03, ContentHash([1; 32])).unwrap()],
        None,
    )
    .unwrap();
    MaterialGasScatteringBinding::try_new(model.evaluate(300.0, 100_000.0).unwrap()).unwrap()
}

struct Density {
    rho: f64,
    gradient: [f64; 3],
    diagonal: [f64; 3],
    off_diagonal: [f64; 3],
}

impl GasDensityGradientField for Density {
    fn sample_density(&self, p: Point3) -> Result<GasDensityGradientSample, GradientIndexError> {
        let [xx, yy, zz] = self.diagonal;
        let [xy, xz, yz] = self.off_diagonal;
        let [gx, gy, gz] = self.gradient;
        Ok(GasDensityGradientSample {
            density_kg_m3: self.rho
                + gx * p.x
                + gy * p.y
                + gz * p.z
                + 0.5 * (xx * p.x * p.x + yy * p.y * p.y + zz * p.z * p.z)
                + xy * p.x * p.y
                + xz * p.x * p.z
                + yz * p.y * p.z,
            gradient_kg_m4: Vec3::new(
                gx + xx * p.x + xy * p.y + xz * p.z,
                gy + xy * p.x + yy * p.y + yz * p.z,
                gz + xz * p.x + yz * p.y + zz * p.z,
            ),
        })
    }
}

impl GasDensitySecondOrderField for Density {
    fn sample_density_second_order(&self, p: Point3) -> Result<GasDensityJet, GradientIndexError> {
        Ok(GasDensityJet {
            sample: self.sample_density(p)?,
            hessian_diagonal_kg_m5: self.diagonal,
            hessian_off_diagonal_kg_m5: self.off_diagonal,
        })
    }
}

fn close(got: f64, expected: f64, tolerance: f64) {
    assert!(
        (got - expected).abs() <= tolerance,
        "{got:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_material_index_hessian_matches_spatial_gradient_differences() {
    let binding = binding();
    for (diagonal, off_diagonal) in [
        ([0.0; 3], [0.0; 3]),
        ([0.04, 0.03, -0.01], [0.02, -0.03, 0.015]),
    ] {
        let density = Density {
            rho: binding.state().gas().density,
            gradient: [0.2, -0.1, 0.05],
            diagonal,
            off_diagonal,
        };
        let field = binding
            .bind_refractive_field(&density, ContentHash([2; 32]))
            .unwrap();
        let p = Point3::new(0.1, 0.3, -0.2);
        for lambda in [380.0, 550.0, 780.0] {
            let jet = field.sample_second_order(p, lambda).unwrap();
            assert_eq!(jet.optical(), field.sample(p, lambda).unwrap());
            let [xx, yy, zz] = jet.hessian().diagonal();
            let [xy, xz, yz] = jet.hessian().off_diagonal();
            let columns = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
            let h = 0.01;
            for (axis, expected) in columns.into_iter().enumerate() {
                let mut a = [p.x, p.y, p.z];
                let mut b = a;
                a[axis] -= h;
                b[axis] += h;
                let ga = field
                    .sample(Point3::new(a[0], a[1], a[2]), lambda)
                    .unwrap()
                    .phase_gradient_per_m();
                let gb = field
                    .sample(Point3::new(b[0], b[1], b[2]), lambda)
                    .unwrap()
                    .phase_gradient_per_m();
                for ((a, b), expected) in [ga.x, ga.y, ga.z]
                    .into_iter()
                    .zip([gb.x, gb.y, gb.z])
                    .zip(expected)
                {
                    close(
                        (b - a) / (2.0 * h),
                        expected,
                        2.0e-8 * expected.abs() + 2.0e-15,
                    );
                }
            }
            if diagonal == [0.0; 3] {
                // An affine DENSITY is not an affine phase index.
                assert!(xx > 0.0 && yy > 0.0 && zz > 0.0);
                assert!(xy < 0.0 && xz > 0.0 && yz < 0.0);
            }
        }
    }
}

#[test]
fn g1_gas_density_curvature_focuses_native_bundle_and_detector_power() {
    let binding = binding();
    let rho = binding.state().gas().density;
    let ray = Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Vec3::new(1.0, 0.0, 0.0),
    };
    let length = 40.0;
    let options = BundleOptions {
        integration: GradientIndexOptions {
            initial_step_m: 0.1,
            max_step_m: 0.5,
            relative_tolerance: 1.0e-11,
            ..Default::default()
        },
        ..Default::default()
    };
    for sign in [-1.0, 0.0, 1.0] {
        let h = sign * 0.5 * rho;
        let density = Density {
            rho,
            gradient: [0.0; 3],
            diagonal: [0.0, h, h],
            off_diagonal: [0.0; 3],
        };
        let field = binding
            .bind_refractive_field(&density, ContentHash([2; 32]))
            .unwrap();
        let result = propagate_bundle(
            &field,
            ray,
            Vec3::new(0.0, 1.0, 0.0),
            550.0,
            length,
            options,
        )
        .unwrap();
        let index = binding.state().refractive_index(550.0).unwrap();
        // Along the symmetry axis n is constant and the exact angular
        // variations obey X'' = (n'(rho) h / n) X, X(0)=0, X'(0)=1.
        let rate = (index.derivative_density_m3_kg * h.abs() / index.phase_index).sqrt();
        let expected = if sign < 0.0 {
            (rate * length).sin() / rate
        } else if sign > 0.0 {
            (rate * length).sinh() / rate
        } else {
            length
        };
        let projection = result
            .project_at_endpoint(Vec3::new(-1.0, 0.0, 0.0), ProjectionOptions::default())
            .unwrap();
        close(
            projection.area_per_solid_angle_m2(),
            expected * expected,
            2.0e-8,
        );
        close(
            projection.point_source_irradiance(7.0).unwrap(),
            7.0 / (expected * expected),
            1.0e-13,
        );
        let central = result.propagation();
        assert_eq!(central.ray.dir, ray.dir);
        assert_eq!(central.ray.origin.y, 0.0);
        assert_eq!(central.ray.origin.z, 0.0);
        close(
            central.optical_path_length_m,
            index.phase_index * length,
            2.0e-11,
        );
        assert_eq!(
            central.group_optical_path_length_m,
            Some(central.optical_path_length_m)
        );
        assert_eq!(central.power_transmittance(), 1.0);
        // Independently launched neighboring rays recover the same focusing;
        // the reference never calls the second-order bundle implementation.
        let epsilon = 1.0e-5;
        let displaced = propagate(
            &field,
            Ray {
                origin: ray.origin,
                dir: Vec3::new(1.0, epsilon, 0.0),
            },
            Vec3::new(0.0, 0.0, 1.0),
            550.0,
            length,
            options.integration,
        )
        .unwrap();
        close(displaced.ray.origin.y / epsilon, expected, 2.0e-6);
    }
}

#[test]
fn g0_second_order_gas_field_preserves_provider_and_numeric_refusals() {
    let binding = binding();
    let mut density = Density {
        rho: 1.0,
        gradient: [0.0; 3],
        diagonal: [0.0; 3],
        off_diagonal: [0.0; 3],
    };
    let p = Point3::new(0.0, 0.0, 0.0);
    let field = binding
        .bind_refractive_field(&density, ContentHash([2; 32]))
        .unwrap();
    assert_eq!(
        field.sample_second_order(Point3::new(f64::NAN, 0.0, 0.0), 550.0),
        Err(GradientIndexError::InvalidSample)
    );
    density.rho = -1.0;
    let field = binding
        .bind_refractive_field(&density, ContentHash([2; 32]))
        .unwrap();
    assert_eq!(
        field.sample_second_order(p, 550.0),
        Err(GradientIndexError::OutsideFieldDomain)
    );
    density.rho = 1.0;
    density.off_diagonal[0] = f64::NAN;
    let field = binding
        .bind_refractive_field(&density, ContentHash([2; 32]))
        .unwrap();
    assert_eq!(
        field.sample_second_order(p, 550.0),
        Err(GradientIndexError::InvalidSample)
    );

    struct Refusing;
    impl GasDensityGradientField for Refusing {
        fn sample_density(
            &self,
            _: Point3,
        ) -> Result<GasDensityGradientSample, GradientIndexError> {
            Err(GradientIndexError::OutsideFieldDomain)
        }
    }
    impl GasDensitySecondOrderField for Refusing {
        fn sample_density_second_order(
            &self,
            _: Point3,
        ) -> Result<GasDensityJet, GradientIndexError> {
            Err(GradientIndexError::Cancelled)
        }
    }
    let field = binding
        .bind_refractive_field(&Refusing, ContentHash([2; 32]))
        .unwrap();
    assert_eq!(
        field.sample_second_order(p, 550.0),
        Err(GradientIndexError::Cancelled)
    );
}
