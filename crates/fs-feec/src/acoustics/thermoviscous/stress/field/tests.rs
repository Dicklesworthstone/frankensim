use super::*;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

fn complex_close(a: C64, b: C64) {
    close(a.re, b.re);
    close(a.im, b.im);
}

// Unit tetrahedron, original hierarchical vertex functions and edge bubbles.
// Coefficients are deliberately not midpoint nodal values.
fn point(bary: [f64; 4]) -> PreparedThermoviscousPoint {
    let gradients = [[-1.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut basis_values = [0.0; 10];
    let mut physical_gradients = [[0.0; 3]; 10];
    basis_values[..4].copy_from_slice(&bary);
    physical_gradients[..4].copy_from_slice(&gradients);
    for (i, &(a, b)) in edges.iter().enumerate() {
        basis_values[i + 4] = bary[a] * bary[b];
        physical_gradients[i + 4] =
            std::array::from_fn(|j| bary[a] * gradients[b][j] + bary[b] * gradients[a][j]);
    }
    PreparedThermoviscousPoint {
        stress: PreparedThermoviscousStress {
            pressure_vertices: [0, 1, 2, 3],
            barycentric: bary,
            velocity_dofs: std::array::from_fn(|i| i),
            physical_gradients,
            pressure_count: 4,
            velocity_count: 10,
            mu: 0.03,
            bulk: 0.02,
        },
        basis_values,
        equilibrium_temperature_k: 300.0,
        conductivity_w_m_k: 0.9,
    }
}

fn coefficients() -> ([C64; 4], [[C64; 3]; 10], [C64; 10]) {
    let p = [1.0, 3.0, 4.0, 5.0].map(|x| C64::new(x, 0.3 * x));
    let real_v = [
        [2.0, -3.0, 0.0],
        [3.0, -3.0, 0.0],
        [2.0, -3.0, -1.0],
        [2.0, -1.0, 0.0],
        [0.0; 3],
        [0.0; 3],
        [0.0; 3],
        [6.0, 0.0, 0.0],
        [0.0, 0.0, 4.0],
        [0.0; 3],
    ];
    let v = real_v.map(|row| row.map(|x| C64::new(x, -0.4 * x)));
    let t = [5.0, 8.0, 3.0, 6.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0].map(|x| C64::new(0.7 * x, 1.2 * x));
    (p, v, t)
}

fn state() -> ThermoviscousTimeState {
    ThermoviscousTimeState {
        pressure_pa: vec![0.0; 4],
        velocity_m_s: vec![[0.0; 3]; 10],
        temperature_k: vec![0.0; 10],
        molecular_temperature_k: vec![],
        conductive_wall_temperature_k: vec![],
        wall_displacement_m_sqrt_kg: vec![],
        wall_velocity_m_sqrt_kg_per_s: vec![],
        wall_relaxation_coordinates_sqrt_j: vec![],
    }
}

#[test]
fn original_p1_p2_fields_gradients_and_fourier_heat_current_are_exact() {
    let observer = point([0.1, 0.2, 0.3, 0.4]);
    let (p, v, t) = coefficients();
    let sample = observer
        .evaluate([4, 10, 10], |i| p[i], |i| v[i], |i| t[i], &mut || false)
        .unwrap();
    complex_close(sample.pressure_pa, C64::new(3.9, 1.17));
    for (actual, expected) in sample.velocity_m_s.into_iter().zip([2.56, -2.2, 0.02]) {
        complex_close(actual, C64::new(expected, -0.4 * expected));
    }
    complex_close(sample.temperature_k, C64::new(0.7 * 5.82, 1.2 * 5.82));
    for (j, expected) in [5.1, -0.6, 1.0].into_iter().enumerate() {
        let gradient = C64::new(0.7 * expected, 1.2 * expected);
        complex_close(sample.temperature_gradient_k_m[j], gradient);
        complex_close(sample.heat_flux_w_m2[j], gradient.scale(-0.9));
    }
    // Analytic off-diagonal stress: mu*(dv_x/dy + dv_y/dx).
    complex_close(
        sample.stress.deviatoric_pa[3],
        C64::new(0.03 * 1.2, -0.012 * 1.2),
    );
}

#[test]
fn harmonic_flux_equals_independent_period_average_of_original_time_states() {
    let observer = point([0.1, 0.2, 0.3, 0.4]);
    let (p, v, t) = coefficients();
    let harmonic = observer
        .evaluate([4, 10, 10], |i| p[i], |i| v[i], |i| t[i], &mut || false)
        .unwrap();
    let expected = harmonic.mean_energy_flux().unwrap();
    let mut sums = [[0.0; 3]; 3];
    let mut shear = 0.0;
    let mut bulk = 0.0;
    let mut instant = state();
    for step in 0..1024 {
        let phase = std::f64::consts::TAU * step as f64 / 1024.0;
        let at = |z: C64| z.re * phase.cos() - z.im * phase.sin();
        for i in 0..4 {
            instant.pressure_pa[i] = at(p[i]);
        }
        for i in 0..10 {
            instant.velocity_m_s[i] = v[i].map(at);
            instant.temperature_k[i] = at(t[i]);
        }
        let sample = observer.instantaneous(&instant, &mut || false).unwrap();
        let flux = sample.energy_flux().unwrap();
        for (sum, channel) in sums.iter_mut().zip([
            flux.pressure_work_w_m2,
            flux.viscous_work_w_m2,
            flux.conductive_work_w_m2,
        ]) {
            for i in 0..3 {
                sum[i] += channel[i] / 1024.0;
            }
        }
        shear += sample.stress.shear_heating_w_m3() / 1024.0;
        bulk += sample.stress.bulk_heating_w_m3() / 1024.0;
    }
    for (actual, expected) in sums.into_iter().zip([
        expected.pressure_work_w_m2,
        expected.viscous_work_w_m2,
        expected.conductive_work_w_m2,
    ]) {
        for i in 0..3 {
            close(actual[i], expected[i]);
        }
    }
    close(shear, harmonic.stress.mean_shear_heating_w_m3());
    close(bulk, harmonic.stress.mean_bulk_heating_w_m3());
}

#[test]
fn viscous_and_thermal_flux_divergence_balance_their_local_dissipation() {
    // Steady Couette v=(2y,0,0) and affine T'=4+3x: exact local energy
    // identities, independent of pressure or a plane-wave approximation.
    let mut instant = state();
    instant.velocity_m_s[2][0] = 2.0;
    instant.temperature_k[..4].copy_from_slice(&[4.0, 7.0, 4.0, 4.0]);
    let evaluate = |x, y| {
        point([1.0 - x - y - 0.1, x, y, 0.1])
            .instantaneous(&instant, &mut || false)
            .unwrap()
    };
    let center = evaluate(0.2, 0.3);
    let a = evaluate(0.15, 0.25).energy_flux().unwrap();
    let b = evaluate(0.25, 0.35).energy_flux().unwrap();
    close(
        -(b.viscous_work_w_m2[1] - a.viscous_work_w_m2[1]) / 0.1,
        center.stress.shear_heating_w_m3(),
    );
    close(
        -(b.conductive_work_w_m2[0] - a.conductive_work_w_m2[0]) / 0.1,
        0.9 / 300.0 * 9.0,
    );
    assert_eq!(center.energy_flux().unwrap().pressure_work_w_m2, [0.0; 3]);
}

#[test]
fn standing_wave_quadrature_retains_velocity_without_inventing_active_intensity() {
    let observer = point([0.25; 4]);
    let sample = observer
        .evaluate(
            [4, 10, 10],
            |_| C64::from_re(5.0),
            |i| {
                if i < 4 {
                    [C64::new(0.0, 2.0), C64::ZERO, C64::ZERO]
                } else {
                    [C64::ZERO; 3]
                }
            },
            |_| C64::ZERO,
            &mut || false,
        )
        .unwrap();
    assert_eq!(sample.velocity_m_s[0], C64::new(0.0, 2.0));
    assert_eq!(
        sample.mean_energy_flux().unwrap().total_w_m2().unwrap(),
        [0.0; 3]
    );
}

#[test]
fn signed_flux_retains_quiet_mechanisms_and_oriented_normal() {
    let flux = ThermoviscousEnergyFlux {
        pressure_work_w_m2: [1e16, 2.0, 0.0],
        viscous_work_w_m2: [3.0, 0.0, 0.0],
        conductive_work_w_m2: [-1e16, 0.0, 0.0],
    };
    assert_eq!(flux.total_w_m2().unwrap(), [3.0, 2.0, 0.0]);
    close(flux.normal_w_m2([3e200, 4e200, 0.0]).unwrap(), 3.4);
    close(flux.normal_w_m2([-3e-200, -4e-200, 0.0]).unwrap(), -3.4);
    assert!(flux.normal_w_m2([0.0; 3]).is_err());
}

#[test]
fn shape_cancellation_and_nonfinite_temperature_refuse_before_returning_fields() {
    let observer = point([0.25; 4]);
    for counts in [[3, 10, 10], [4, 9, 10], [4, 10, 9]] {
        assert!(
            observer
                .evaluate(
                    counts,
                    |_| panic!("shape first"),
                    |_| panic!("shape first"),
                    |_| panic!("shape first"),
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(matches!(
        observer.evaluate(
            [4, 10, 10],
            |_| panic!("cancel first"),
            |_| panic!("cancel first"),
            |_| panic!("cancel first"),
            &mut || true
        ),
        Err(ThermoviscousError::Cancelled)
    ));
    assert!(matches!(
        observer.evaluate(
            [4, 10, 10],
            |_| C64::ZERO,
            |_| [C64::ZERO; 3],
            |i| C64::from_re(if i == 8 { f64::NAN } else { 0.0 }),
            &mut || false
        ),
        Err(ThermoviscousError::NonFinite)
    ));
    let mut calls = 0;
    let (p, v, t) = coefficients();
    assert!(matches!(
        observer.evaluate([4, 10, 10], |i| p[i], |i| v[i], |i| t[i], &mut || {
            calls += 1;
            calls >= 3
        }),
        Err(ThermoviscousError::Cancelled)
    ));
    assert!(
        observer
            .evaluate([4, 10, 10], |i| p[i], |i| v[i], |i| t[i], &mut || false)
            .is_ok()
    );
}
