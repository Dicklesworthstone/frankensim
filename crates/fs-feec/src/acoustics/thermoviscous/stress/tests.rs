use super::*;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(1e-28),
        "{a:e} != {b:e}"
    );
}

fn complex_close(a: C64, b: C64) {
    close(a.re, b.re);
    close(a.im, b.im);
}

#[test]
fn g1_complex_couette_stress_has_tensor_shear_and_cycle_mean_heating() {
    let mu = 0.007;
    let shear = C64::new(3.0, -4.0);
    let mut gradient = [[C64::ZERO; 3]; 3];
    gradient[0][1] = shear;
    let value = ThermoviscousStress::from_gradient(C64::new(7.0, 2.0), gradient, mu, 0.02).unwrap();
    complex_close(value.deviatoric_pa[3], shear.scale(mu));
    close(value.mean_shear_heating_w_m3(), 0.5 * mu * 25.0);
    assert_eq!(value.mean_bulk_heating_w_m3(), 0.0);
    let traction = value.traction([0.0, 19.0, 0.0]).unwrap();
    assert_eq!(traction.shear_pa, [shear.scale(mu), C64::ZERO, C64::ZERO]);
    assert_eq!(traction.normal_pa, C64::new(-7.0, -2.0));
    let opposite = value.traction([0.0, -1.0, 0.0]).unwrap();
    for i in 0..3 {
        complex_close(opposite.total_pa[i], traction.total_pa[i].scale(-1.0));
    }
    assert_eq!(opposite.normal_pa, traction.normal_pa);
}

#[test]
fn g3_rigid_spin_and_isotropic_expansion_have_no_shear() {
    let expansion = C64::new(2.0, -3.0);
    let spin = C64::new(4.0, 5.0);
    let gradient = [
        [expansion, spin, C64::ZERO],
        [spin.scale(-1.0), expansion, C64::ZERO],
        [C64::ZERO, C64::ZERO, expansion],
    ];
    let value = ThermoviscousStress::from_gradient(C64::ZERO, gradient, 0.17, 0.23).unwrap();
    assert_eq!(value.deviatoric_pa, [C64::ZERO; 6]);
    assert_eq!(value.mean_shear_heating_w_m3(), 0.0);
    complex_close(value.bulk_pa, expansion.scale(3.0 * 0.23));
    close(value.mean_bulk_heating_w_m3(), 0.5 * 0.23 * 9.0 * 13.0);
    assert_eq!(
        value.traction([1.0, 2.0, 3.0]).unwrap().shear_pa,
        [C64::ZERO; 3]
    );
}

#[test]
fn g3_pressure_cannot_hide_wall_shear_and_normal_scaling_is_invariant() {
    let mut gradient = [[C64::ZERO; 3]; 3];
    gradient[0][1] = C64::new(0.2, -0.1);
    let low = ThermoviscousStress::from_gradient(C64::ZERO, gradient, 1e-5, 0.0).unwrap();
    let high =
        ThermoviscousStress::from_gradient(C64::new(1e30, -1e25), gradient, 1e-5, 0.0).unwrap();
    let expected = low.traction([1.0, 2.0, 3.0]).unwrap().shear_pa;
    assert_eq!(high.traction([1.0, 2.0, 3.0]).unwrap().shear_pa, expected);
    assert_eq!(
        high.traction([1e300, 2e300, 3e300]).unwrap().shear_pa,
        expected
    );
    assert_eq!(
        high.traction([1e-300, 2e-300, 3e-300]).unwrap().shear_pa,
        expected
    );
}

#[test]
fn g3_stress_and_traction_rotate_with_the_physical_axes() {
    let g = [[0.3, -1.2, 0.7], [0.8, -0.4, 0.2], [-0.5, 1.1, 0.6]];
    let c = fs_math::det::cos(0.61);
    let s = fs_math::det::sin(0.61);
    let rotation = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
    let rotate_vector = |v: [f64; 3]| -> [f64; 3] {
        std::array::from_fn(|i| (0..3).map(|j| rotation[i][j] * v[j]).sum())
    };
    let rotated: [[C64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            C64::from_re(
                (0..3)
                    .flat_map(|a| (0..3).map(move |b| rotation[i][a] * g[a][b] * rotation[j][b]))
                    .sum(),
            )
        })
    });
    let original = ThermoviscousStress::from_gradient(
        C64::from_re(3.4),
        g.map(|r| r.map(C64::from_re)),
        0.6,
        0.8,
    )
    .unwrap();
    let transformed =
        ThermoviscousStress::from_gradient(C64::from_re(3.4), rotated, 0.6, 0.8).unwrap();
    close(
        original.mean_shear_heating_w_m3(),
        transformed.mean_shear_heating_w_m3(),
    );
    close(
        original.mean_bulk_heating_w_m3(),
        transformed.mean_bulk_heating_w_m3(),
    );
    let normal = [1.0, -2.0, 3.0];
    let a = original.traction(normal).unwrap();
    let b = transformed.traction(rotate_vector(normal)).unwrap();
    close(a.normal_pa.re, b.normal_pa.re);
    let expected = rotate_vector(a.shear_pa.map(|v| v.re));
    for i in 0..3 {
        close(expected[i], b.shear_pa[i].re);
    }
}

// One exact P2 element observation: v_x = 6 lambda_1 lambda_2,
// with lambda_1=x, lambda_2=y on the reference tetrahedron.
fn quadratic_observer() -> PreparedThermoviscousStress {
    let mut gradients = [[0.0; 3]; 10];
    gradients[0] = [-1.0; 3];
    gradients[1] = [1.0, 0.0, 0.0];
    gradients[2] = [0.0, 1.0, 0.0];
    gradients[3] = [0.0, 0.0, 1.0];
    gradients[7] = [0.3, 0.2, 0.0];
    PreparedThermoviscousStress {
        pressure_vertices: [0, 1, 2, 3],
        barycentric: [0.1, 0.2, 0.3, 0.4],
        velocity_dofs: std::array::from_fn(|i| i),
        physical_gradients: gradients,
        pressure_count: 4,
        velocity_count: 10,
        mu: 0.5,
        bulk: 0.7,
    }
}

#[test]
fn g1_original_hierarchical_coefficients_and_instantaneous_heating_are_retained() {
    let observer = quadratic_observer();
    let mut state = ThermoviscousTimeState {
        pressure_pa: vec![1.0, 2.0, 3.0, 4.0],
        velocity_m_s: vec![[0.0; 3]; 10],
        temperature_k: vec![0.0; 10],
        molecular_temperature_k: vec![],
        conductive_wall_temperature_k: vec![],
        wall_displacement_m_sqrt_kg: vec![],
        wall_velocity_m_sqrt_kg_per_s: vec![],
    };
    state.velocity_m_s[7][0] = 6.0;
    let real = observer.instantaneous(&state, &mut || false).unwrap();
    close(real.pressure_pa, 3.0);
    close(real.deviatoric_pa[3], 0.6);
    close(real.bulk_pa, 0.7 * 1.8);
    // Dxx=1.8, Dxy=0.6, dev(D)=(1.2,-0.6,-0.6,0.6,0,0).
    close(
        real.shear_heating_w_m3(),
        2.0 * 0.5 * (1.44 + 0.36 + 0.36 + 2.0 * 0.36),
    );
    close(real.bulk_heating_w_m3(), 0.7 * 1.8 * 1.8);
    let complex = observer
        .evaluate(
            4,
            10,
            |i| C64::from_re(state.pressure_pa[i]),
            |i| state.velocity_m_s[i].map(C64::from_re),
            &mut || false,
        )
        .unwrap();
    close(
        real.shear_heating_w_m3(),
        2.0 * complex.mean_shear_heating_w_m3(),
    );
    close(
        real.bulk_heating_w_m3(),
        2.0 * complex.mean_bulk_heating_w_m3(),
    );
    let a = real.traction([0.0, 1.0, 0.0]).unwrap();
    let b = complex.traction([0.0, 1.0, 0.0]).unwrap();
    assert_eq!(a.shear_pa, b.shear_pa.map(|v| v.re));
}

#[test]
fn invalid_normals_fields_and_cancelled_samples_refuse() {
    let value =
        ThermoviscousStress::from_gradient(C64::ZERO, [[C64::ZERO; 3]; 3], 1.0, 0.0).unwrap();
    for n in [[0.0; 3], [f64::NAN, 0.0, 1.0], [0.0, f64::INFINITY, 0.0]] {
        assert!(value.traction(n).is_err());
    }
    assert!(
        ThermoviscousStress::from_gradient(C64::new(f64::NAN, 0.0), [[C64::ZERO; 3]; 3], 1.0, 0.0)
            .is_err()
    );
    let observer = quadratic_observer();
    assert!(
        observer
            .evaluate(
                3,
                10,
                |_| panic!("must check dimensions first"),
                |_| panic!("must check dimensions first"),
                &mut || false
            )
            .is_err()
    );
    assert!(matches!(
        observer.evaluate(
            4,
            10,
            |_| panic!("cancelled"),
            |_| panic!("cancelled"),
            &mut || true
        ),
        Err(ThermoviscousError::Cancelled)
    ));
}
