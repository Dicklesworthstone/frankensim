//! G1/G3: resolved complex-material pressure waves, interfaces and local heat.

use fs_feec::acoustics::helmholtz::heterogeneous::{
    BoundaryAdmittance, EquivalentFluidCoefficients, EquivalentFluidControls,
    HarmonicLimits, HeterogeneousError, HeterogeneousHelmholtzFe, PreparedEquivalentFluid,
    VelocityTrace,
};
use fs_feec::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_feec::fixtures::{masked_cube_grid, single_tet};
use fs_math::c64::C64;
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.3,
    sound_speed_m_s: 2.0,
};
const MATERIAL: EquivalentFluidCoefficients = EquivalentFluidCoefficients {
    dynamic_density_kg_m3: C64::new(1.3, -0.4),
    dynamic_bulk_modulus_pa: C64::new(5.2, 0.8),
};

fn brick(nx: usize) -> (TetComplex, Vec<[f64; 3]>) {
    let (complex, mut positions) = masked_cube_grid(nx, 1, 1, &|_, _, _| true);
    for point in &mut positions {
        point[0] /= nx as f64;
        point[1] *= 0.2;
        point[2] *= 0.3;
    }
    (complex, positions)
}

fn plane(fe: &HelmholtzFe<'_>, coordinate: f64) -> Vec<[u32; 3]> {
    fe.boundary_faces()
        .into_iter()
        .filter(|face| face.iter().all(|&v| fe.positions()[v as usize][0] == coordinate))
        .collect()
}

fn velocity(faces: &[[u32; 3]], value: C64) -> Vec<VelocityTrace> {
    faces.iter().map(|&face| (face, [value; 3])).collect()
}

fn wave(k: C64, x: f64) -> C64 {
    let phase = -k.re * x;
    C64::new(phase.cos(), phase.sin()).scale((k.im * x).exp())
}

fn assert_relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}",
    );
}

#[test]
fn lossy_plane_wave_converges_with_separate_original_volume_losses() {
    let omega = 2.2;
    let rho = MATERIAL.dynamic_density_kg_m3;
    let bulk = MATERIAL.dynamic_bulk_modulus_pa;
    let k = (rho / bulk).sqrt().scale(omega);
    let impedance = (rho * bulk).sqrt();
    let area = 0.06;
    let integral = (2.0 * k.im).exp_m1() / (2.0 * k.im);
    let expected_density = 0.5 / omega * rho.recip().im
        * k.norm_sq() * impedance.norm_sq() * area * integral;
    let expected_compression = -0.5 * omega * bulk.recip().im
        * impedance.norm_sq() * area * integral;
    let expected_input = 0.5 * impedance.re * area;
    let expected_outflow = expected_input * (2.0 * k.im).exp();
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (complex, positions) = brick(nx);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let boundaries = [BoundaryAdmittance {
            faces: plane(&geometry, 1.0),
            admittance_m_pa_s: impedance.recip(),
        }];
        let prepared = PreparedEquivalentFluid::new(
            &geometry,
            omega / core::f64::consts::TAU,
            &boundaries,
            &[],
            EquivalentFluidControls::default(),
            &mut |_, _| Ok::<_, HeterogeneousError>(MATERIAL),
            &mut || false,
        ).unwrap();
        let response = prepared
            .solve(&velocity(&plane(&geometry, 0.0), C64::ONE), &mut || false)
            .unwrap();
        let error = positions.iter().zip(&response.pressure_pa)
            .map(|(point, &p)| {
                (p - impedance * wave(k, point[0])).abs() / impedance.abs()
            }).fold(0.0, f64::max);
        errors.push(error);
        assert!(response.density_dissipation_w > 0.0
            && response.compressional_dissipation_w > 0.0);
        assert!(response.power_balance_residual_w().abs() < 2e-11 * expected_input);
        assert!(response.backward_error < 5e-13);
        assert_eq!(response.element_power.len(), complex.tets.len());
        assert_relative(
            response.element_power.iter().map(|p| p.volume_m3).sum(), 0.06, 2e-13,
        );
        for element in &response.element_power {
            assert!(element.density_dissipation_w >= 0.0
                && element.compressional_dissipation_w >= 0.0);
        }
        if nx == 8 {
            assert_relative(response.density_dissipation_w, expected_density, 2e-4);
            assert_relative(response.compressional_dissipation_w, expected_compression, 2e-4);
            assert_relative(response.input_power_w, expected_input, 2e-4);
            assert_relative(response.boundary_dissipation_w, expected_outflow, 2e-4);
        }
    }
    assert!(errors[1] < errors[0] / 4.0 && errors[2] < errors[1] / 4.0, "{errors:?}");
    assert!(errors[2] < 1e-4, "{errors:?}");
}

#[test]
fn conforming_free_fluid_and_lossy_volume_preserve_interface_reflection() {
    let (complex, positions) = brick(8);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let left = EquivalentFluidCoefficients {
        dynamic_density_kg_m3: C64::from_re(MEDIUM.density_kg_m3),
        dynamic_bulk_modulus_pa: C64::from_re(
            MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2),
        ),
    };
    let right = EquivalentFluidCoefficients {
        dynamic_density_kg_m3: C64::new(4.6, -1.1),
        dynamic_bulk_modulus_pa: C64::new(8.0, 1.3),
    };
    let omega = 1.4;
    let z1 = (left.dynamic_density_kg_m3 * left.dynamic_bulk_modulus_pa).sqrt();
    let z2 = (right.dynamic_density_kg_m3 * right.dynamic_bulk_modulus_pa).sqrt();
    let k1 = (left.dynamic_density_kg_m3 / left.dynamic_bulk_modulus_pa).sqrt().scale(omega);
    let k2 = (right.dynamic_density_kg_m3 / right.dynamic_bulk_modulus_pa).sqrt().scale(omega);
    let reflection = (z2 - z1) / (z2 + z1);
    let inlet = plane(&geometry, 0.0);
    let boundaries = [
        BoundaryAdmittance { faces: inlet.clone(), admittance_m_pa_s: z1.recip() },
        BoundaryAdmittance { faces: plane(&geometry, 1.0), admittance_m_pa_s: z2.recip() },
    ];
    let model = PreparedEquivalentFluid::new(
        &geometry,
        omega / core::f64::consts::TAU,
        &boundaries,
        &[],
        EquivalentFluidControls::default(),
        &mut |_, p| Ok::<_, HeterogeneousError>(if p[0] < 0.5 { left } else { right }),
        &mut || false,
    ).unwrap();
    let response = model
        .solve(&velocity(&inlet, C64::from_re(2.0) / z1), &mut || false)
        .unwrap();
    for (&p, point) in response.pressure_pa.iter().zip(&positions) {
        let x = point[0];
        let exact = if x <= 0.5 {
            wave(k1, x) + reflection * wave(k1, 1.0 - x)
        } else {
            (C64::ONE + reflection) * wave(k1, 0.5) * wave(k2, x - 0.5)
        };
        assert!((p - exact).abs() < 7e-5, "x={x}, p={p:?}, exact={exact:?}");
    }
    for (tet, power) in complex.tets.iter().zip(&response.element_power) {
        if tet.iter().all(|&v| positions[v as usize][0] <= 0.5) {
            assert_eq!(power.density_dissipation_w, 0.0);
            assert_eq!(power.compressional_dissipation_w, 0.0);
        }
    }
    assert!(response.power_balance_residual_w().abs() < 2e-12);
}

#[test]
fn local_heat_uses_original_complex_gradient_and_pressure_in_si_units() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let omega = 3.0;
    let prepared = PreparedEquivalentFluid::new(
        &geometry,
        omega / core::f64::consts::TAU,
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut |_, _| Ok::<_, HeterogeneousError>(MATERIAL),
        &mut || false,
    ).unwrap();
    let mut pressure = vec![C64::ZERO; geometry.ndof()];
    for (i, p) in positions.iter().enumerate() {
        pressure[i] = C64::new(2.0 + p[0], 2.0 * p[1]);
    }
    let losses = prepared.element_dissipation(&pressure, &mut || false).unwrap();
    let volume = 1.0 / 6.0;
    // On the unit tetrahedron E[x]=1/4 and E[x²]=E[y²]=1/10.
    // |grad(2+x+2iy)|²=5, E[|2+x+2iy|²]=5.5.
    assert_relative(losses[0].volume_m3, volume, 2e-14);
    assert_relative(
        losses[0].density_dissipation_w,
        volume * MATERIAL.dynamic_density_kg_m3.recip().im * 5.0 / (2.0 * omega),
        2e-13,
    );
    assert_relative(
        losses[0].compressional_dissipation_w,
        -volume * MATERIAL.dynamic_bulk_modulus_pa.recip().im * 5.5 * omega / 2.0,
        2e-13,
    );
    let rotated: Vec<_> = pressure.iter()
        .map(|&p| C64::new(-p.im, p.re).scale(3.0)).collect();
    let amplified = prepared.element_dissipation(&rotated, &mut || false).unwrap();
    assert_relative(
        amplified[0].total_dissipation_w(),
        9.0 * losses[0].total_dissipation_w(),
        2e-14,
    );
    assert!(prepared.element_dissipation(&pressure[..pressure.len() - 1], &mut || false).is_err());
    assert!(matches!(
        prepared.element_dissipation(&pressure, &mut || true),
        Err(HeterogeneousError::Cancelled),
    ));
}

#[test]
fn lossless_limit_reuses_the_original_full_pressure_solution() {
    let (complex, positions) = brick(3);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let original = HeterogeneousHelmholtzFe::new(
        &geometry, 4, &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM), &mut || false,
    ).unwrap();
    let coefficients = EquivalentFluidCoefficients {
        dynamic_density_kg_m3: C64::from_re(MEDIUM.density_kg_m3),
        dynamic_bulk_modulus_pa: C64::from_re(
            MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2),
        ),
    };
    let end = plane(&geometry, 1.0);
    let old = original.prepare_with_pressure_release(
        0.2, &[], &end, HarmonicLimits::default(), &mut || false,
    ).unwrap();
    let prepared = PreparedEquivalentFluid::new(
        &geometry,
        0.2,
        &[],
        &end,
        EquivalentFluidControls::default(),
        &mut |_, _| Ok::<_, HeterogeneousError>(coefficients),
        &mut || false,
    ).unwrap();
    let drive = velocity(&plane(&geometry, 0.0), C64::new(0.3, -0.4));
    let expected = old.solve(&drive, &mut || false).unwrap();
    let actual = prepared.solve(&drive, &mut || false).unwrap();
    for (&a, &b) in actual.pressure_pa.iter().zip(&expected.pressure_pa) {
        assert!((a - b).abs() < 2e-12, "{a:?} versus {b:?}");
    }
    assert_eq!(actual.density_dissipation_w, 0.0);
    assert_eq!(actual.compressional_dissipation_w, 0.0);
    assert_eq!(actual.boundary_dissipation_w, 0.0);
    assert!(prepared.solve(&velocity(&end, C64::ONE), &mut || false).is_err());
    assert!(matches!(
        prepared.solve(&drive, &mut || true),
        Err(HeterogeneousError::Cancelled),
    ));
    assert_eq!(
        prepared.solve(&drive, &mut || false).unwrap().pressure_pa,
        actual.pressure_pa,
    );
}

#[test]
fn preparation_refuses_active_material_and_limits_before_sampling() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    for controls in [
        EquivalentFluidControls { max_volume_points: 0, ..EquivalentFluidControls::default() },
        EquivalentFluidControls { max_element_loss_bytes: 0, ..EquivalentFluidControls::default() },
        EquivalentFluidControls {
            harmonic: HarmonicLimits { max_dofs: 1, ..HarmonicLimits::default() },
            ..EquivalentFluidControls::default()
        },
    ] {
        let calls = Cell::new(0);
        let result = PreparedEquivalentFluid::new(
            &geometry,
            0.3,
            &[],
            &[],
            controls,
            &mut |_, _| {
                calls.set(calls.get() + 1);
                Ok::<_, HeterogeneousError>(MATERIAL)
            },
            &mut || false,
        );
        assert!(matches!(result, Err(HeterogeneousError::Limit { .. })));
        assert_eq!(calls.get(), 0);
    }
    for bad in [
        EquivalentFluidCoefficients {
            dynamic_density_kg_m3: MATERIAL.dynamic_density_kg_m3.conj(), ..MATERIAL
        },
        EquivalentFluidCoefficients {
            dynamic_bulk_modulus_pa: MATERIAL.dynamic_bulk_modulus_pa.conj(), ..MATERIAL
        },
    ] {
        assert!(PreparedEquivalentFluid::new(
            &geometry,
            0.3,
            &[],
            &[],
            EquivalentFluidControls::default(),
            &mut |_, _| Ok::<_, HeterogeneousError>(bad),
            &mut || false,
        ).is_err());
    }
    let calls = Cell::new(0);
    let result = PreparedEquivalentFluid::new(
        &geometry,
        0.3,
        &[],
        &[],
        EquivalentFluidControls::default(),
        &mut |_, _| {
            calls.set(calls.get() + 1);
            Ok::<_, HeterogeneousError>(MATERIAL)
        },
        &mut || calls.get() >= 3,
    );
    assert!(matches!(result, Err(HeterogeneousError::Cancelled)));
    assert_eq!(calls.get(), 3);
}
