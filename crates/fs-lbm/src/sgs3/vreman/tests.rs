use super::*;

fn close(a: f64, b: f64) {
    assert!((a / b - 1.0).abs() < 3e-14, "{a:e} != {b:e}");
}

#[test]
fn all_thirteen_laminar_flow_types_and_exact_rank_one_have_zero_viscosity() {
    let widths = [0.125, 2.0, 8.0];
    assert_eq!(vreman_eddy_viscosity(0.07, widths, [[0.0; 3]; 3]), Ok(0.0));
    let mut cases = 1;
    for i in 0..3 {
        for j in 0..3 {
            if i != j {
                let mut g = [[0.0; 3]; 3];
                g[i][j] = -3.0;
                assert_eq!(vreman_eddy_viscosity(0.07, widths, g), Ok(0.0));
                cases += 1;
            }
        }
        let mut row = [[0.0; 3]; 3];
        row[i][(i + 1) % 3] = 2.0;
        row[i][(i + 2) % 3] = -4.0;
        let column = core::array::from_fn(|a| core::array::from_fn(|b| row[b][a]));
        for g in [row, column] {
            assert_eq!(vreman_eddy_viscosity(0.07, widths, g), Ok(0.0));
            cases += 1;
        }
    }
    assert_eq!(cases, 13);
    let a = [1.0, 2.0, -4.0];
    let b = [0.5, -0.25, 1.0];
    let g = core::array::from_fn(|i| core::array::from_fn(|j| a[i] * b[j]));
    assert_eq!(vreman_eddy_viscosity(0.07, widths, g), Ok(0.0));
}

#[test]
fn directional_widths_use_derivative_axes_and_match_an_independent_beta_oracle() {
    let g = [[0.0, 2.0, 0.0], [0.0, 0.0, 3.0], [5.0, 0.0, 0.0]];
    let d = [0.25, 2.0, 4.0];
    let beta = [4.0 * d[1] * d[1], 9.0 * d[2] * d[2], 25.0 * d[0] * d[0]];
    let expected = 0.07 * ((beta[0] * beta[1] + beta[0] * beta[2] + beta[1] * beta[2]) / 38.0_f64).sqrt();
    let value = vreman_eddy_viscosity(0.07, d, g).unwrap();
    close(value, expected);
    assert!((vreman_eddy_viscosity(0.07, [d[2], d[1], d[0]], g).unwrap() / value - 1.0).abs() > 0.1);
    for p in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
        let permuted = core::array::from_fn(|i| core::array::from_fn(|j| g[p[i]][p[j]]));
        close(vreman_eddy_viscosity(0.07, p.map(|i| d[i]), permuted).unwrap(), value);
    }
    close(vreman_eddy_viscosity(0.14, d, g).unwrap(), 2.0 * value);
}

#[test]
fn isotropic_frame_invariance_and_quadratic_filter_refinement() {
    let g = [[0.3, 0.7, -0.2], [-0.1, -0.4, 0.6], [0.8, -0.5, 0.1]];
    let r = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    let rotated = core::array::from_fn(|i| core::array::from_fn(|j| {
        (0..3).flat_map(|a| (0..3).map(move |b| r[i][a] * g[a][b] * r[j][b])).sum()
    }));
    let base = vreman_eddy_viscosity(0.07, [1.0; 3], g).unwrap();
    close(vreman_eddy_viscosity(0.07, [1.0; 3], rotated).unwrap(), base);
    close(vreman_eddy_viscosity(0.07, [0.5; 3], g).unwrap(), base / 4.0);
    close(vreman_eddy_viscosity(0.07, [0.25; 3], g).unwrap(), base / 16.0);
}

#[test]
fn nearly_rank_one_flow_keeps_small_positive_transport_without_clipping() {
    let epsilon = 2.0_f64.powi(-40);
    let g = [[1.0, 1.0, 0.0], [1.0, 1.0 + epsilon, 0.0], [0.0; 3]];
    let expected = 0.07 * epsilon / (3.0 + (1.0 + epsilon).powi(2)).sqrt();
    let result = vreman_eddy_viscosity(0.07, [1.0; 3], g).unwrap();
    assert!(result > 0.0);
    close(result, expected);
}

#[test]
fn dimensional_reconstruction_avoids_intermediate_overflow_and_underflow() {
    let g = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]];
    for (gradient, width, coefficient, expected) in [
        (1e300, 1e-100, 0.07, 7e98),
        (1e-300, 1e200, 0.07, 7e98),
        (1.0, 1e200, 1e-300, 1e100),
        (1e-200, 1.0, 0.07, 7e-202),
        (1e200, 1.0, 0.07, 7e198),
    ] {
        close(vreman_eddy_viscosity(coefficient, [width; 3],
            g.map(|row| row.map(|v| v * gradient))).unwrap(), expected);
    }
    let extreme_aspect = vreman_eddy_viscosity(1.0, [1e-150, 1.0, 1e150], g).unwrap();
    close(extreme_aspect, 1e150 / 3.0_f64.sqrt());
    assert_eq!(vreman_eddy_viscosity(1.0, [1e200; 3], g), Err(SgsError::NonFiniteViscosity));
    assert_eq!(vreman_eddy_viscosity(1.0, [1e-200; 3], g), Err(SgsError::NonFiniteViscosity));
}

#[test]
fn malformed_inputs_refuse_even_for_disabled_control() {
    let g = [[0.0; 3]; 3];
    for coefficient in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(vreman_eddy_viscosity(coefficient, [1.0; 3], g), Err(SgsError::InvalidCoefficient));
    }
    for width in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        for axis in 0..3 {
            let mut d = [1.0; 3];
            d[axis] = width;
            assert_eq!(vreman_eddy_viscosity(0.0, d, g), Err(SgsError::InvalidFilterWidth));
        }
    }
    for value in [f64::NAN, f64::INFINITY] {
        assert_eq!(vreman_eddy_viscosity(0.0, [1.0; 3], [[value; 3]; 3]), Err(SgsError::NonFiniteGradient));
    }
}
