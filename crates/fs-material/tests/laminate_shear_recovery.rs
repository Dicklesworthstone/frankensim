use fs_material::laminate::{
    LaminaLayer, LaminateSection, LaminateShearRecovery, LayerLaw, ReferenceSurface,
};
use std::sync::Arc;

fn layer(h: f64, angle: f64, scale: f64) -> LaminaLayer {
    LaminaLayer {
        thickness_m: h,
        orientation_rad: angle,
        density_kg_m3: 1400.0,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [
                12e9 * scale,
                0.8e9 * scale,
                0.0,
                0.8e9 * scale,
                4e9 * scale,
                0.0,
                0.0,
                0.0,
                1.4e9 * scale,
            ],
            transverse_shear_pa: [1.2e9 * scale, 0.25e9 * scale],
        },
    }
}
fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn recovery(
    layers: Vec<LaminaLayer>,
    reference: ReferenceSurface,
    axis: f64,
) -> LaminateShearRecovery {
    LaminateShearRecovery::try_new(Arc::new(
        LaminateSection::try_new_equilibrated_shear(layers, reference, axis).unwrap(),
    ))
    .unwrap()
}

#[test]
fn homogeneous_profile_is_physical_parabola_with_uncorrected_local_compliance() {
    let h = 0.008;
    let angle = 0.41;
    let source = Arc::new(
        LaminateSection::try_new_equilibrated_shear(
            vec![layer(h, angle, 1.0)],
            ReferenceSurface::Midsurface,
            -0.17,
        )
        .unwrap(),
    );
    let response = LaminateShearRecovery::try_new(source.clone()).unwrap();
    assert!(std::ptr::eq(response.section(), source.as_ref()));
    let q = [2100.0, -800.0];
    for s in [0.0, 0.07, 0.25, 0.5, 0.91, 1.0] {
        let point = response.evaluate(0, s, q).unwrap();
        for k in 0..2 {
            near(point.stress_pa[k], 6.0 * s * (1.0 - s) * q[k] / h, 1e6);
        }
        near(
            point.material_engineering_strain[0],
            point.material_stress_pa[0] / 1.2e9,
            1e-3,
        );
        near(
            point.material_engineering_strain[1],
            point.material_stress_pa[1] / 0.25e9,
            1e-3,
        );
        near(point.z_m, h * (s - 0.5), h);
        near(
            point.energy_density_j_m3,
            0.5 * (point.stress_pa[0] * point.engineering_strain[0]
                + point.stress_pa[1] * point.engineering_strain[1]),
            1.0,
        );
    }
    // Peak stress is 3Q/(2h), not the first-order constant G*gamma.
    near(
        response.evaluate(0, 0.5, q).unwrap().stress_pa[0],
        1.5 * q[0] / h,
        1e6,
    );
}

#[test]
fn unequal_off_axis_plies_close_resultant_energy_and_interface_traction() {
    let source = recovery(
        vec![
            layer(0.002, 0.37, 1.0),
            layer(0.005, -0.63, 0.13),
            layer(0.001, 0.82, 2.4),
        ],
        ReferenceSurface::Offset(0.0013),
        0.21,
    );
    let h = source.section().shear_pa_m();
    assert!(h[1].abs() > 1e-3 * h[0]);
    for gamma in [[0.002, 0.003], [0.002, -0.003]] {
        let q = [
            h[0] * gamma[0] + h[1] * gamma[1],
            h[2] * gamma[0] + h[3] * gamma[1],
        ];
        let mut resultant = [0.0; 2];
        let mut energy = 0.0;
        let g = (3.0_f64 / 5.0).sqrt();
        for (i, layer) in source.section().layers().iter().enumerate() {
            for (s, w) in [
                (0.5 * (1.0 - g), 5.0 / 18.0),
                (0.5, 4.0 / 9.0),
                (0.5 * (1.0 + g), 5.0 / 18.0),
            ] {
                let point = source.from_section_strain(i, s, gamma).unwrap();
                for k in 0..2 {
                    resultant[k] += w * layer.thickness_m * point.stress_pa[k];
                }
                energy += w * layer.thickness_m * point.energy_density_j_m3;
                assert!(point.energy_density_j_m3 >= 0.0);
            }
        }
        for k in 0..2 {
            near(resultant[k], q[k], q[0].abs().max(q[1].abs()));
        }
        near(energy, 0.5 * (q[0] * gamma[0] + q[1] * gamma[1]), energy);
        for i in 0..2 {
            let lower = source.evaluate(i, 1.0, q).unwrap();
            let upper = source.evaluate(i + 1, 0.0, q).unwrap();
            for k in 0..2 {
                near(lower.stress_pa[k], upper.stress_pa[k], 1e7);
            }
            assert!((lower.engineering_strain[0] - upper.engineering_strain[0]).abs() > 1e-5);
        }
        for (i, s) in [(0, 0.0), (2, 1.0)] {
            for tau in source.evaluate(i, s, q).unwrap().stress_pa {
                near(tau, 0.0, 1e7);
            }
        }
    }
}

#[test]
fn physical_axis_rotation_and_reference_shift_preserve_recovered_ply_work() {
    let layers = vec![layer(0.003, 0.2, 1.0), layer(0.001, -0.5, 0.2)];
    let a = recovery(layers.clone(), ReferenceSurface::Midsurface, 0.13);
    let phi: f64 = -0.74;
    let (s, c) = phi.sin_cos();
    let rotate = |v: [f64; 2]| [c * v[0] - s * v[1], s * v[0] + c * v[1]];
    let b = recovery(
        layers
            .into_iter()
            .map(|mut l| {
                l.orientation_rad += phi;
                l
            })
            .collect(),
        ReferenceSurface::Offset(0.009),
        0.13 + phi,
    );
    let q = [-2000.0, 1700.0];
    for i in 0..2 {
        for fraction in [0.1, 0.5, 0.9] {
            let x = a.evaluate(i, fraction, q).unwrap();
            let y = b.evaluate(i, fraction, rotate(q)).unwrap();
            for k in 0..2 {
                near(y.stress_pa[k], rotate(x.stress_pa)[k], 1e6);
                near(
                    y.engineering_strain[k],
                    rotate(x.engineering_strain)[k],
                    1e-2,
                );
                near(y.material_stress_pa[k], x.material_stress_pa[k], 1e6);
                near(
                    y.material_engineering_strain[k],
                    x.material_engineering_strain[k],
                    1e-2,
                );
            }
            near(y.energy_density_j_m3, x.energy_density_j_m3, 1.0);
            near(y.z_m, x.z_m - 0.007, 0.01);
        }
    }
}

#[test]
fn wrong_closure_and_invalid_points_refuse_without_changing_sources() {
    let whitney = Arc::new(
        LaminateSection::try_new(vec![layer(0.01, 0.0, 1.0)], ReferenceSurface::Midsurface)
            .unwrap(),
    );
    assert!(LaminateShearRecovery::try_new(whitney).is_err());
    let source = recovery(
        vec![layer(0.01, 0.0, 1.0)],
        ReferenceSurface::Midsurface,
        0.0,
    );
    for s in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(source.evaluate(0, s, [1.0; 2]).is_err());
    }
    assert!(source.evaluate(1, 0.5, [1.0; 2]).is_err());
    assert!(source.evaluate(0, 0.5, [f64::NAN, 0.0]).is_err());
    assert!(
        source
            .from_section_strain(0, 0.5, [f64::INFINITY, 0.0])
            .is_err()
    );
    let zero = source.evaluate(0, 0.5, [0.0; 2]).unwrap();
    assert_eq!(zero.stress_pa, [0.0; 2]);
    assert_eq!(zero.energy_density_j_m3, 0.0);
}
