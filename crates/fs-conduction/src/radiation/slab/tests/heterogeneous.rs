use super::*;

fn prepared(depths: &[f64], indices: &[f64]) -> AbsorbingSlab {
    AbsorbingSlab::new(
        depths
            .iter()
            .zip(indices)
            .enumerate()
            .map(|(i, (&depth, &index))| {
                AbsorbingLayer::new(1.0, depth, [100.0, 3000.0], format!("absorption {i}"))
                    .unwrap()
                    .with_refractive_index(index, format!("index {i}"))
                    .unwrap()
            })
            .collect(),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

// Independent physical two-stream boundary equations at each ordinate. Dense
// pivoted elimination retains both incoming directions of every layer; it does
// not use the production prefix/suffix elimination or exchange coefficients.
fn ray(
    depths: &[f64],
    indices: &[f64],
    p: f64,
    pol: usize,
    b: &[f64],
    h: [f64; 2],
) -> (Vec<f64>, [f64; 2]) {
    let n = depths.len();
    let mut all = vec![1.0];
    all.extend_from_slice(indices);
    all.push(1.0);
    let mu: Vec<_> = all
        .iter()
        .map(|&index| {
            if p < index {
                (1.0 - (p / index).powi(2)).sqrt()
            } else {
                0.0
            }
        })
        .collect();
    let mut r = vec![1.0; n + 1];
    for i in 0..=n {
        if p < all[i] && p < all[i + 1] {
            let (a, z) = if pol == 0 {
                (all[i] * mu[i], all[i + 1] * mu[i + 1])
            } else {
                (all[i] * mu[i + 1], all[i + 1] * mu[i])
            };
            r[i] = ((a - z) / (a + z)).powi(2);
        }
    }
    let t: Vec<_> = (0..n)
        .map(|i| {
            if mu[i + 1] > 0.0 {
                (-depths[i] / mu[i + 1]).exp()
            } else {
                0.0
            }
        })
        .collect();
    let a: Vec<_> = (0..n)
        .map(|i| {
            if mu[i + 1] > 0.0 {
                -(-depths[i] / mu[i + 1]).exp_m1()
            } else {
                0.0
            }
        })
        .collect();
    let size = 2 * n;
    let mut matrix = vec![vec![0.0; size + 1]; size];
    for i in 0..n {
        let x = 2 * i;
        let y = x + 1;
        matrix[x][x] = 1.0;
        matrix[x][y] = -r[i] * t[i];
        matrix[x][size] = r[i] * a[i] * b[i];
        if i == 0 {
            matrix[x][size] += (1.0 - r[i]) * h[0];
        } else {
            matrix[x][2 * (i - 1)] = -(1.0 - r[i]) * t[i - 1];
            matrix[x][size] += (1.0 - r[i]) * a[i - 1] * b[i - 1];
        }
        matrix[y][y] = 1.0;
        matrix[y][x] = -r[i + 1] * t[i];
        matrix[y][size] = r[i + 1] * a[i] * b[i];
        if i + 1 == n {
            matrix[y][size] += (1.0 - r[i + 1]) * h[1];
        } else {
            matrix[y][2 * (i + 1) + 1] = -(1.0 - r[i + 1]) * t[i + 1];
            matrix[y][size] += (1.0 - r[i + 1]) * a[i + 1] * b[i + 1];
        }
    }
    for k in 0..size {
        let pivot = (k..size)
            .max_by(|&i, &j| matrix[i][k].abs().total_cmp(&matrix[j][k].abs()))
            .unwrap();
        matrix.swap(k, pivot);
        let divisor = matrix[k][k];
        assert!(divisor.abs() > 1e-13);
        for j in k..=size {
            matrix[k][j] /= divisor;
        }
        for i in 0..size {
            if i != k {
                let factor = matrix[i][k];
                for j in k..=size {
                    matrix[i][j] -= factor * matrix[k][j];
                }
            }
        }
    }
    let heat = (0..n)
        .map(|i| a[i] * (matrix[2 * i][size] + matrix[2 * i + 1][size] - 2.0 * b[i]))
        .collect();
    let outgoing = if p < 1.0 {
        [
            r[0] * h[0] + (1.0 - r[0]) * (t[0] * matrix[1][size] + a[0] * b[0]),
            r[n] * h[1]
                + (1.0 - r[n]) * (t[n - 1] * matrix[2 * (n - 1)][size] + a[n - 1] * b[n - 1]),
        ]
    } else {
        [0.0; 2]
    };
    (heat, outgoing)
}

fn reference(
    depths: &[f64],
    indices: &[f64],
    b: &[f64],
    h: [f64; 2],
    exterior_only: bool,
) -> (Vec<f64>, [f64; 2]) {
    let mut edges = vec![0.0, 1.0];
    if !exterior_only {
        edges.extend_from_slice(indices);
    }
    edges.sort_by(f64::total_cmp);
    edges.dedup();
    let panels = 4096;
    let step = core::f64::consts::FRAC_PI_2 / panels as f64;
    let mut heat = vec![0.0; depths.len()];
    let mut out = [0.0; 2];
    for band in edges.windows(2) {
        let width = band[1] * band[1] - band[0] * band[0];
        for j in 0..panels {
            let theta = (j as f64 + 0.5) * step;
            let (sin, cos) = theta.sin_cos();
            let p = (band[0] * band[0] + width * sin * sin).sqrt();
            let weight = width * sin * cos * step;
            for pol in 0..2 {
                let (q, j) = ray(depths, indices, p, pol, b, h);
                for i in 0..heat.len() {
                    heat[i] += weight * q[i];
                }
                for i in 0..2 {
                    out[i] += weight * j[i];
                }
            }
        }
    }
    (heat, out)
}

#[test]
fn g1_heterogeneous_heat_matches_independent_boundary_equations() {
    let depths = [0.2, 0.9, 0.4];
    let indices = [1.4, 2.1, 1.2];
    let slab = prepared(&depths, &indices);
    let b = [3.0, 7.0, 2.0];
    let h = [11.0, 5.0];
    let mut heat = [0.0; 3];
    slab.layer_net_fluxes(h, &b, &mut heat).unwrap();
    let (expected, out) = reference(&depths, &indices, &b, h, false);
    let emitted = slab.emission_to_faces(&b).unwrap();
    for i in 0..3 {
        close(heat[i], expected[i], 2e-6);
    }
    let actual: [f64; 2] = std::array::from_fn(|i| {
        emitted[i] + slab.reflectances()[i] * h[i] + slab.transmittance() * h[1 - i]
    });
    for i in 0..2 {
        close(actual[i], out[i], 2e-6);
    }
    close(
        heat.iter().sum::<f64>() + actual.iter().sum::<f64>(),
        h.iter().sum(),
        2e-13,
    );
    assert!((slab.reflectances()[0] - slab.reflectances()[1]).abs() > 0.01);
    assert_eq!(slab.refractive_index(), None);
    assert_eq!(slab.layers()[1].refractive_index(), 2.1);
    assert_eq!(slab.layers()[1].refractive_index_source(), Some("index 1"));
    assert_eq!(slab.refractive_index_source(), None);
}

#[test]
fn g3_heterogeneous_reversal_subdivision_and_equilibrium() {
    let original = prepared(&[0.2, 0.9, 0.4], &[1.4, 2.1, 1.2]);
    let reverse = prepared(&[0.4, 0.9, 0.2], &[1.2, 2.1, 1.4]);
    let split = prepared(&[0.2, 0.3, 0.6, 0.4], &[1.4, 2.1, 2.1, 1.2]);
    close(original.transmittance(), reverse.transmittance(), 2e-12);
    close(original.transmittance(), split.transmittance(), 2e-12);
    for side in 0..2 {
        close(
            original.reflectances()[side],
            reverse.reflectances()[1 - side],
            2e-12,
        );
        close(
            original.reflectances()[side],
            split.reflectances()[side],
            2e-12,
        );
    }
    let mut heat = [0.0; 3];
    let mut reversed = [0.0; 3];
    let mut divided = [0.0; 4];
    original
        .layer_net_fluxes([11.0, 5.0], &[3.0, 7.0, 2.0], &mut heat)
        .unwrap();
    reverse
        .layer_net_fluxes([5.0, 11.0], &[2.0, 7.0, 3.0], &mut reversed)
        .unwrap();
    split
        .layer_net_fluxes([11.0, 5.0], &[3.0, 7.0, 7.0, 2.0], &mut divided)
        .unwrap();
    for i in 0..3 {
        close(heat[i], reversed[2 - i], 2e-11);
    }
    close(heat[0], divided[0], 2e-11);
    close(heat[1], divided[1] + divided[2], 2e-11);
    close(heat[2], divided[3], 2e-11);
    original
        .layer_net_fluxes([7.0; 2], &[7.0; 3], &mut heat)
        .unwrap();
    assert_eq!(heat, [0.0; 3]);
    for i in 0..3 {
        for j in 0..3 {
            assert_eq!(
                original.layer_coupling()[3 * i + j],
                original.layer_coupling()[3 * j + i]
            );
        }
    }
}

#[test]
fn g1_lossless_layers_and_low_index_gap_retain_internal_reflection() {
    let indices = [1.4, 2.1, 1.2];
    let clear = prepared(&[0.0; 3], &indices);
    let (_, out) = reference(&[0.0; 3], &indices, &[0.0; 3], [1.0, 0.0], true);
    close(clear.reflectances()[0], out[0], 2e-7);
    close(clear.transmittance(), out[1], 2e-7);
    close(clear.reflectances()[0], clear.reflectances()[1], 2e-13);
    assert!(clear.layer_coupling().iter().all(|&x| x == 0.0));
    assert!(clear.face_coupling().iter().flatten().all(|&x| x == 0.0));
    let gap = prepared(&[f64::MAX, 0.0, f64::MAX], &[2.0, 1.0, 2.0]);
    // Only p<1 propagates across the real low-index gap. The transparent
    // two-interface transmission has the same reciprocal Fresnel series.
    let expected = super::refractive(&[0.0], 2.0).transmittance();
    close(gap.layer_coupling()[2], expected, 2e-12);
    assert!(gap.layer_coupling()[2] < 1.0);
    let common = super::refractive(&[f64::MAX, 0.0, f64::MAX], 2.0);
    close(common.layer_coupling()[2], 4.0, 2e-12);
}

#[test]
fn g0_trapped_weak_absorption_keeps_positive_thermal_exchange() {
    let small = prepared(&[1e-20, 2e-20, 3e-20], &[2.0, 3.0, 2.0]);
    let doubled = prepared(&[2e-20, 4e-20, 6e-20], &[2.0, 3.0, 2.0]);
    for (i, j) in [(0, 1), (0, 2), (1, 2)] {
        let value = small.layer_coupling()[3 * i + j];
        assert!(value.is_finite() && value > 1e-21);
        assert!((doubled.layer_coupling()[3 * i + j] / value - 2.0).abs() < 2e-12);
    }
    let opaque = prepared(&[f64::MAX, f64::MAX], &[1.5, 2.0]);
    assert_eq!(opaque.transmittance(), 0.0);
    assert!(opaque.layer_coupling()[1] > 1.0 && opaque.layer_coupling()[1] < 2.25);
}

#[test]
fn g0_index_admission_and_uniform_constructor_reuse() {
    let layer = AbsorbingLayer::new(1.0, 0.2, [100.0, 3000.0], "absorption").unwrap();
    for index in [0.0, 0.99, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            layer
                .clone()
                .with_refractive_index(index, "source")
                .is_err()
        );
    }
    assert!(layer.clone().with_refractive_index(1.5, " ").is_err());
    let selected = layer.clone().with_refractive_index(1.5, "source").unwrap();
    let direct = AbsorbingSlab::new(vec![selected.clone(); 2], SlabQuadrature::default(), || {
        false
    })
    .unwrap();
    let common = AbsorbingSlab::new_refractive(
        vec![layer.clone(); 2],
        1.5,
        "source",
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    assert_eq!(direct, common);
    assert_eq!(direct.refractive_index(), Some(1.5));
    assert_eq!(direct.refractive_index_source(), Some("source"));
    let layers = vec![
        selected,
        layer.with_refractive_index(2.0, "second source").unwrap(),
    ];
    assert!(
        AbsorbingSlab::new(
            layers.clone(),
            SlabQuadrature {
                maximum_evaluations: 24,
                ..SlabQuadrature::default()
            },
            || false
        )
        .is_err()
    );
    assert!(matches!(
        AbsorbingSlab::new(layers, SlabQuadrature::default(), || true),
        Err(ConductionError::Cancelled { .. })
    ));
}
