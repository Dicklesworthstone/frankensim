use super::*;

const G: [[f64; 3]; 3] = [[0.3, 0.7, -0.2], [-0.1, -0.4, 0.6], [0.8, -0.5, 0.1]];
const MODELS: [SgsModel; 3] = [
    SgsModel::Smagorinsky { cs: 0.17 },
    SgsModel::Wale { cw: 0.5 },
    SgsModel::Vreman { cv: 0.07 },
];

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:e}, expected={expected:e}"
    );
}

fn contract(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> f64 {
    (0..3)
        .flat_map(|i| (0..3).map(move |j| a[i][j] * b[i][j]))
        .sum()
}

fn finite_differences(
    gradient: [[f64; 3]; 3],
    response: SgsViscosityTangent,
    primal: impl Fn([[f64; 3]; 3]) -> Result<f64, SgsError>,
) {
    assert_eq!(
        response.viscosity.to_bits(),
        primal(gradient).unwrap().to_bits()
    );
    let h = 2.0_f64.powi(-18);
    for i in 0..3 {
        for j in 0..3 {
            let mut plus = gradient;
            let mut minus = gradient;
            plus[i][j] += h;
            minus[i][j] -= h;
            let difference = (primal(plus).unwrap() - primal(minus).unwrap()) / (2.0 * h);
            close(response.derivative[i][j], difference, 2e-8, 2e-11);
        }
    }
    // A degree-one viscosity satisfies Euler's identity at smooth gradients.
    close(
        contract(response.derivative, gradient),
        response.viscosity,
        4e-14,
        2e-16,
    );
}

#[test]
fn all_nine_entries_match_independent_centered_perturbations() {
    for g in [
        G,
        [[0.7, 1.2, 0.0], [-0.4, -0.7, 0.0], [0.0; 3]],
        [[0.8, -0.2, 0.3], [0.4, 0.1, -0.5], [-0.2, 0.6, -0.9]],
    ] {
        for model in MODELS {
            finite_differences(g, eddy_viscosity_tangent(model, 1.3, g).unwrap(), |trial| {
                eddy_viscosity(model, 1.3, trial)
            });
        }
    }
}

#[test]
fn anisotropic_vreman_differentiates_column_widths_and_preserves_primal() {
    let widths = [0.125, 2.0, 8.0];
    finite_differences(
        G,
        vreman_eddy_viscosity_tangent(0.07, widths, G).unwrap(),
        |trial| vreman_eddy_viscosity(0.07, widths, trial),
    );
    let response = vreman_eddy_viscosity_tangent(0.07, widths, G).unwrap();
    for p in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let permuted = core::array::from_fn(|i| core::array::from_fn(|j| G[p[i]][p[j]]));
        let actual = vreman_eddy_viscosity_tangent(0.07, p.map(|i| widths[i]), permuted).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(
                    actual.derivative[i][j],
                    response.derivative[p[i]][p[j]],
                    4e-14,
                    1e-15,
                );
            }
        }
    }
    let wrong = vreman_eddy_viscosity_tangent(0.07, [8.0, 2.0, 0.125], G).unwrap();
    let difference = core::array::from_fn(|i| {
        core::array::from_fn(|j| wrong.derivative[i][j] - response.derivative[i][j])
    });
    assert!(maximum(difference) > 0.1 * maximum(response.derivative));
}

#[test]
fn isotropic_derivatives_rotate_as_tensors_and_are_odd_in_gradient() {
    let r = [[0.6, -0.8, 0.0], [0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
    let rotate = |g: [[f64; 3]; 3]| {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..3)
                    .flat_map(|a| (0..3).map(move |b| r[i][a] * g[a][b] * r[j][b]))
                    .sum()
            })
        })
    };
    for model in MODELS {
        let original = eddy_viscosity_tangent(model, 1.3, G).unwrap();
        let rotated = eddy_viscosity_tangent(model, 1.3, rotate(G)).unwrap();
        let expected = rotate(original.derivative);
        let negative = eddy_viscosity_tangent(model, 1.3, G.map(|row| row.map(|v| -v))).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(rotated.derivative[i][j], expected[i][j], 5e-14, 2e-16);
                assert_eq!(negative.derivative[i][j], -original.derivative[i][j]);
            }
        }
    }
}

#[test]
fn physical_gradient_and_filter_scaling_do_not_change_the_tangent_law() {
    for model in MODELS {
        let base = eddy_viscosity_tangent(model, 1.0, G).unwrap();
        for scale in [1e-150, 1.0, 1e150] {
            let g = G.map(|row| row.map(|v| scale * v));
            for width in [1e-50, 0.5, 1.0, 1e50] {
                let actual = eddy_viscosity_tangent(model, width, g).unwrap();
                assert_eq!(
                    actual.viscosity.to_bits(),
                    eddy_viscosity(model, width, g).unwrap().to_bits()
                );
                for i in 0..3 {
                    for j in 0..3 {
                        close(
                            actual.derivative[i][j] / width / width,
                            base.derivative[i][j],
                            8e-14,
                            2e-16,
                        );
                    }
                }
            }
        }
    }
    for (model, width, reference) in [
        (
            SgsModel::Smagorinsky { cs: 1e-200 },
            1e200,
            SgsModel::Smagorinsky { cs: 1.0 },
        ),
        (
            SgsModel::Wale { cw: 1e-200 },
            1e200,
            SgsModel::Wale { cw: 1.0 },
        ),
        (
            SgsModel::Vreman { cv: 1e-300 },
            1e150,
            SgsModel::Vreman { cv: 1.0 },
        ),
    ] {
        let actual = eddy_viscosity_tangent(model, width, G).unwrap();
        let expected = eddy_viscosity_tangent(reference, 1.0, G).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                close(
                    actual.derivative[i][j],
                    expected.derivative[i][j],
                    8e-14,
                    2e-16,
                );
            }
        }
    }
}

#[test]
fn near_wall_wale_derivative_retains_the_cubic_asymptote() {
    let dy = [[0.2, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, 0.0, 0.3]];
    let mut previous = None;
    for y in [0.01, 0.005, 0.0025] {
        let g = [
            [0.2 * y, 1.0, 0.0],
            [0.0, -0.5 * y, 0.0],
            [0.0, 0.4, 0.3 * y],
        ];
        let actual = eddy_viscosity_tangent(SgsModel::Wale { cw: 0.5 }, 1.0, g).unwrap();
        let derivative = contract(actual.derivative, dy);
        close(y * derivative / actual.viscosity, 3.0, 0.0, 0.002);
        if let Some(prior) = previous {
            close(prior / derivative, 4.0, 0.0, 0.002);
        }
        previous = Some(derivative);
    }
}

#[test]
fn zero_selections_are_explicit_and_rank_one_vreman_remains_nonsmooth() {
    let zero = [[0.0; 3]; 3];
    for model in MODELS {
        assert_eq!(
            eddy_viscosity_tangent(model, 1.0, zero).unwrap(),
            SgsViscosityTangent {
                viscosity: 0.0,
                derivative: zero
            }
        );
    }
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        assert_eq!(
            eddy_viscosity_tangent(model, 1.0, G).unwrap(),
            SgsViscosityTangent {
                viscosity: 0.0,
                derivative: zero
            }
        );
    }
    let rotation = [[0.0, -0.2, 0.0], [0.2, 0.0, 0.0], [0.0; 3]];
    assert_eq!(
        eddy_viscosity_tangent(MODELS[0], 1.0, rotation)
            .unwrap()
            .derivative,
        zero
    );
    let shear = [[0.0, 1.0, 0.0], [0.0; 3], [0.0; 3]];
    assert_eq!(
        eddy_viscosity_tangent(MODELS[1], 1.0, shear)
            .unwrap()
            .derivative,
        zero
    );
    let rank_one = [[1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    assert_eq!(
        eddy_viscosity_tangent(MODELS[2], 1.0, rank_one)
            .unwrap()
            .derivative,
        zero
    );
    for h in [-1e-6, 1e-6] {
        let mut perturbed = rank_one;
        perturbed[1][1] = h;
        let slope = eddy_viscosity(MODELS[2], 1.0, perturbed).unwrap() / h;
        close(slope, 0.07_f64.copysign(h), 1e-12, 0.0);
    }
}

#[test]
fn rotation_perturbations_and_complete_stress_response_need_the_full_gradient() {
    let spin = [[0.0, 0.8, -0.2], [-0.8, 0.0, 0.3], [0.2, -0.3, 0.0]];
    let h = 2.0_f64.powi(-18);
    let plus = core::array::from_fn(|i| core::array::from_fn(|j| G[i][j] + h * spin[i][j]));
    let minus = core::array::from_fn(|i| core::array::from_fn(|j| G[i][j] - h * spin[i][j]));
    for model in MODELS {
        let response = eddy_viscosity_tangent(model, 1.3, G).unwrap();
        let derivative = contract(response.derivative, spin);
        if matches!(model, SgsModel::Smagorinsky { .. }) {
            close(derivative, 0.0, 0.0, 1e-17);
        } else {
            assert!(derivative.abs() > 1e-4);
        }
        let nu_plus = eddy_viscosity(model, 1.3, plus).unwrap();
        let nu_minus = eddy_viscosity(model, 1.3, minus).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let analytic = derivative * (G[i][j] + G[j][i])
                    + response.viscosity * (spin[i][j] + spin[j][i]);
                let difference = (nu_plus * (plus[i][j] + plus[j][i])
                    - nu_minus * (minus[i][j] + minus[j][i]))
                    / (2.0 * h);
                close(analytic, difference, 2e-8, 2e-11);
            }
        }
    }
}

#[test]
fn anisotropic_extreme_widths_and_nearly_rank_one_gradients_are_not_clipped() {
    let widths = [1e-150, 1.0, 1e150];
    let response = vreman_eddy_viscosity_tangent(0.07, widths, G).unwrap();
    assert_eq!(
        response.viscosity.to_bits(),
        vreman_eddy_viscosity(0.07, widths, G).unwrap().to_bits()
    );
    close(
        contract(response.derivative, G),
        response.viscosity,
        8e-14,
        0.0,
    );
    let epsilon = 2.0_f64.powi(-40);
    let g = [[1.0, 1.0, 0.0], [1.0, 1.0 + epsilon, 0.0], [0.0; 3]];
    let actual = vreman_eddy_viscosity_tangent(0.07, [1.0; 3], g).unwrap();
    let norm = (3.0 + (1.0 + epsilon).powi(2)).sqrt();
    let expected = 0.07 * (1.0 / norm - epsilon * (1.0 + epsilon) / norm.powi(3));
    close(actual.derivative[1][1], expected, 4e-14, 0.0);
    assert!(actual.derivative[1][1] > 0.03);
}

#[test]
fn unrepresentable_derivatives_refuse_even_when_primal_viscosity_fits() {
    for model in [
        SgsModel::Smagorinsky { cs: 1.0 },
        SgsModel::Wale { cw: 1.0 },
        SgsModel::Vreman { cv: 1.0 },
    ] {
        for (width, rate) in [(1e200, 1e-200), (1e-200, 1e200)] {
            let g = G.map(|row| row.map(|v| v * rate));
            assert!(eddy_viscosity(model, width, g).unwrap().is_finite());
            assert_eq!(
                eddy_viscosity_tangent(model, width, g),
                Err(SgsError::NonFiniteViscosity)
            );
        }
    }
}

#[test]
fn malformed_inputs_keep_original_admission_even_when_disabled() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        assert_eq!(
            eddy_viscosity_tangent(model, 0.0, G),
            Err(SgsError::InvalidFilterWidth)
        );
        assert_eq!(
            eddy_viscosity_tangent(model, 1.0, [[f64::NAN; 3]; 3]),
            Err(SgsError::NonFiniteGradient)
        );
    }
    for model in [
        SgsModel::Smagorinsky { cs: -1.0 },
        SgsModel::Wale { cw: f64::INFINITY },
        SgsModel::Vreman { cv: f64::NAN },
    ] {
        assert_eq!(
            eddy_viscosity_tangent(model, 1.0, G),
            Err(SgsError::InvalidCoefficient)
        );
    }
    assert_eq!(
        vreman_eddy_viscosity_tangent(0.0, [1.0, -1.0, 1.0], G),
        Err(SgsError::InvalidFilterWidth)
    );
}
