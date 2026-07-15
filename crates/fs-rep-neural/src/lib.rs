//! fs-rep-neural — neural implicit charts. Layer: L2.
//!
//! Small coordinate MLPs as shapes (DeepSDF-style). The load-bearing property
//! is NOT the fit — it is the CERTIFICATE: with SPECTRAL-NORM-constrained
//! layers and 1-Lipschitz activations, the network carries a certified global
//! LIPSCHITZ CONSTANT `L = Π Uᵢ`, where every `Uᵢ` is an outward-rounded upper
//! bound on a layer's largest singular value. That bound holds for every input
//! regardless of training. That is what keeps a neural shape inside the
//! certificate regime:
//!
//! - sphere tracing with PROVABLE step safety — a step of `|f(x)|/L` can never
//!   tunnel through the surface ([`safe_step_radius`]);
//! - INTERVAL evaluation via layer-wise bound propagation (IBP) — a guaranteed
//!   output enclosure over an input box ([`MlpSdf::eval_interval`]);
//! - the gradient is bounded by `L` (`‖∇f‖ ≤ L`).
//!
//! Neural charts PROPOSE geometry; certification (watertightness, Hausdorff
//! agreement, topology) comes from the certificate machinery, NEVER the loss
//! curve — so [`TopologyHint`] is honestly `Unknown` here. Deterministic; no
//! dependencies (in-house power iteration for diagnostics and scaled
//! Frobenius/induced-norm upper bounds for certificates).

/// A dense affine layer (`out × in` weights + bias).
#[derive(Debug, Clone, PartialEq)]
pub struct Layer {
    /// The weight matrix (`out` rows, each of length `in`).
    pub weights: Vec<Vec<f64>>,
    /// The bias (`out`).
    pub bias: Vec<f64>,
}

impl Layer {
    /// A layer.
    ///
    /// # Panics
    /// If the weight rows disagree with the bias length, are ragged, or contain
    /// non-finite values. Biases must also be finite.
    #[must_use]
    pub fn new(weights: Vec<Vec<f64>>, bias: Vec<f64>) -> Layer {
        assert!(weights.len() == bias.len(), "weights/bias length mismatch");
        if let Some(first) = weights.first() {
            assert!(
                weights.iter().all(|r| r.len() == first.len()),
                "ragged weights"
            );
        }
        assert!(
            weights.iter().flatten().all(|w| w.is_finite()),
            "weights must be finite"
        );
        assert!(bias.iter().all(|b| b.is_finite()), "bias must be finite");
        Layer { weights, bias }
    }

    fn in_dim(&self) -> usize {
        self.weights.first().map_or(0, Vec::len)
    }
}

/// A deterministic power-iteration estimate of the spectral norm (largest
/// singular value) of a matrix.
///
/// This estimate is retained for diagnostics and backwards compatibility. It
/// can converge from below or miss a singular direction orthogonal to its
/// fixed starting vector, so it MUST NOT be used as a safety certificate. Use
/// [`spectral_norm_upper_bound`] for admission, normalization, and Lipschitz
/// claims.
#[must_use]
pub fn spectral_norm(weights: &[Vec<f64>]) -> f64 {
    if weights.is_empty() || weights[0].is_empty() {
        return 0.0;
    }
    let n_in = weights[0].len();
    // deterministic, non-degenerate initial vector.
    let mut v: Vec<f64> = (0..n_in).map(|i| 1.0 + 0.1 * i as f64).collect();
    normalize(&mut v);
    for _ in 0..200 {
        let u = matvec(weights, &v); // W v
        let mut z = matvec_t(weights, &u); // Wᵀ (W v)
        let nz = norm(&z);
        if nz < 1e-30 {
            return 0.0;
        }
        for zi in &mut z {
            *zi /= nz;
        }
        v = z;
    }
    norm(&matvec(weights, &v)) // ‖W v‖ estimates σ_max from below.
}

/// A guaranteed finite upper bound on a matrix's spectral norm.
///
/// The result is the tighter of the scaled Frobenius bound and
/// `sqrt(||W||_1 ||W||_∞)`, each evaluated with outward rounding at every
/// positive arithmetic operation. The result can conservatively exceed the
/// true spectral norm by as much as `sqrt(rank)`, but unlike power iteration it
/// cannot miss a singular direction. Scaling by the largest entry avoids
/// intermediate square/sum overflow.
///
/// # Panics
///
/// Panics for ragged or non-finite matrices, or when no finite `f64` can
/// represent the certified upper bound. The latter is a fail-closed admission
/// outcome: returning `f64::MAX` would be unsound when the mathematical norm is
/// larger than `f64::MAX`.
#[must_use]
pub fn spectral_norm_upper_bound(weights: &[Vec<f64>]) -> f64 {
    let (max_abs, unit_bound) = scaled_spectral_upper_parts(weights);
    finite_scaled_upper(max_abs, unit_bound)
        .unwrap_or_else(|| panic!("spectral-norm upper bound is not representable as a finite f64"))
}

/// Normalize a layer so its weight matrix has spectral norm at most `bound`.
///
/// This uses [`spectral_norm_upper_bound`], not the power-iteration estimate.
/// The stored matrix is re-certified after floating-point scaling; if an
/// extreme rounding case cannot be corrected in a few deterministic passes,
/// the weights are collapsed to zero as a fail-closed fallback.
///
/// # Panics
///
/// Panics if `bound` is negative or non-finite, or if the layer was mutated
/// after construction into a ragged or non-finite state.
#[must_use]
pub fn spectral_normalize(mut layer: Layer, bound: f64) -> Layer {
    assert!(
        bound.is_finite() && bound >= 0.0,
        "spectral bound must be finite and non-negative"
    );
    assert!(
        layer.bias.iter().all(|b| b.is_finite()),
        "bias must be finite"
    );
    assert!(
        layer.weights.len() == layer.bias.len(),
        "weights/bias length mismatch"
    );

    let (max_abs, unit_bound) = scaled_spectral_upper_parts(&layer.weights);
    if max_abs == 0.0 {
        return layer;
    }
    if bound == 0.0 {
        zero_weights(&mut layer.weights);
        return layer;
    }

    // Form each normalized value relative to max_abs instead of first forming
    // bound / (max_abs * unit_bound). This remains useful when the original
    // matrix's norm exceeds f64::MAX or max_abs is subnormal.
    let target = next_down_nonnegative(bound / unit_bound);
    for row in &mut layer.weights {
        for w in row {
            *w = (*w / max_abs) * target;
        }
    }

    // Re-certify the values actually stored after rounded division and
    // multiplication. Ordinarily one pass suffices; the fixed limit keeps the
    // operation deterministic even at subnormal boundaries.
    for _ in 0..4 {
        match try_spectral_norm_upper_bound(&layer.weights) {
            Some(certified) if certified <= bound => return layer,
            Some(certified) => {
                let correction = next_down_nonnegative(bound / certified);
                for row in &mut layer.weights {
                    for w in row {
                        *w *= correction;
                    }
                }
            }
            None => {
                // The current certificate is not finitely representable yet.
                for row in &mut layer.weights {
                    for w in row {
                        *w *= 0.5;
                    }
                }
            }
        }
    }

    // A zero matrix is always within the requested bound. Reaching this path
    // requires an extreme rounding plateau; safety takes precedence over
    // preserving a non-certifiable parameterization.
    zero_weights(&mut layer.weights);
    layer
}

/// An honest topology hint: never claimed without a verifying check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TopologyHint {
    /// Topology is unknown until certified (persistent-homology check).
    Unknown,
}

/// A neural implicit chart: a coordinate MLP with a certified Lipschitz constant.
#[derive(Debug, Clone, PartialEq)]
pub struct MlpSdf {
    layers: Vec<Layer>,
    lipschitz: f64,
}

impl MlpSdf {
    /// Build a chart from raw layers, spectrally normalizing each to at most
    /// `bound` (so `L ≤ boundᵏ` with 1-Lipschitz `tanh` activations). The last
    /// layer must map to a scalar.
    ///
    /// # Panics
    /// If `layers` is empty, dimensions do not chain, the output is not scalar,
    /// `bound` is invalid, or the global certificate is not finitely
    /// representable.
    #[must_use]
    pub fn new(layers: Vec<Layer>, bound: f64) -> MlpSdf {
        assert!(!layers.is_empty(), "need at least one layer");
        assert!(
            bound.is_finite() && bound >= 0.0,
            "spectral bound must be finite and non-negative"
        );
        assert!(
            layers.last().unwrap().weights.len() == 1,
            "output must be scalar"
        );
        for pair in layers.windows(2) {
            assert!(
                pair[0].weights.len() == pair[1].in_dim(),
                "layer dims do not chain"
            );
        }
        let normalized: Vec<Layer> = layers
            .into_iter()
            .map(|l| spectral_normalize(l, bound))
            .collect();
        // Certified global Lipschitz constant = Π Uᵢ, where each Uᵢ is a
        // guaranteed spectral-norm upper bound (tanh is 1-Lipschitz).
        let lipschitz = normalized.iter().fold(1.0, |product, layer| {
            product_up(product, spectral_norm_upper_bound(&layer.weights))
        });
        assert!(
            lipschitz.is_finite(),
            "global Lipschitz certificate is not representable as a finite f64"
        );
        MlpSdf {
            layers: normalized,
            lipschitz,
        }
    }

    /// The certified global Lipschitz constant.
    #[must_use]
    pub fn lipschitz(&self) -> f64 {
        self.lipschitz
    }

    /// The topology hint — honestly `Unknown` (never inferred from the fit).
    #[must_use]
    pub fn topology_hint(&self) -> TopologyHint {
        TopologyHint::Unknown
    }

    /// Evaluate the implicit value `f(x)` (tanh between layers, linear output).
    #[must_use]
    pub fn eval(&self, x: &[f64]) -> f64 {
        let mut a = x.to_vec();
        let last = self.layers.len() - 1;
        for (li, layer) in self.layers.iter().enumerate() {
            let mut z = layer.bias.clone();
            for (j, row) in layer.weights.iter().enumerate() {
                z[j] += dot(row, &a);
            }
            if li < last {
                for zi in &mut z {
                    *zi = zi.tanh();
                }
            }
            a = z;
        }
        a[0]
    }

    /// The gradient `∇f(x)` by central finite differences (its norm is `≤ L`).
    #[must_use]
    pub fn eval_grad(&self, x: &[f64]) -> Vec<f64> {
        let h = 1e-6;
        (0..x.len())
            .map(|i| {
                let mut xp = x.to_vec();
                let mut xm = x.to_vec();
                xp[i] += h;
                xm[i] -= h;
                (self.eval(&xp) - self.eval(&xm)) / (2.0 * h)
            })
            .collect()
    }

    /// A GUARANTEED output enclosure of `f` over the input box `[lo, hi]`, by
    /// interval bound propagation (IBP) — sound for sphere-tracing sub-boxes.
    #[must_use]
    pub fn eval_interval(&self, lo: &[f64], hi: &[f64]) -> (f64, f64) {
        let (mut lo, mut hi) = (lo.to_vec(), hi.to_vec());
        let last = self.layers.len() - 1;
        for (li, layer) in self.layers.iter().enumerate() {
            let mut nlo = layer.bias.clone();
            let mut nhi = layer.bias.clone();
            for (j, row) in layer.weights.iter().enumerate() {
                for (i, &w) in row.iter().enumerate() {
                    if w >= 0.0 {
                        nlo[j] += w * lo[i];
                        nhi[j] += w * hi[i];
                    } else {
                        nlo[j] += w * hi[i];
                        nhi[j] += w * lo[i];
                    }
                }
            }
            if li < last {
                for k in 0..nlo.len() {
                    // tanh is monotone increasing.
                    nlo[k] = nlo[k].tanh();
                    nhi[k] = nhi[k].tanh();
                }
            }
            lo = nlo;
            hi = nhi;
        }
        (lo[0], hi[0])
    }
}

/// The provably-safe sphere-tracing step radius: with SDF value `value` and
/// Lipschitz constant `lipschitz`, `f` cannot change sign within `|value|/L`, so
/// a step of that size never tunnels through the surface. The returned finite
/// quotient is rounded DOWN; a nearest-rounded quotient can exceed the exact
/// safe radius by one ulp. Invalid inputs fail closed to a zero step.
#[must_use]
pub fn safe_step_radius(value: f64, lipschitz: f64) -> f64 {
    if !value.is_finite() || !lipschitz.is_finite() || lipschitz < 0.0 {
        return 0.0;
    }
    if lipschitz == 0.0 {
        return f64::INFINITY;
    }
    next_down_nonnegative(value.abs() / lipschitz)
}

// -- linear-algebra helpers -------------------------------------------------

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn matvec(w: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    w.iter().map(|row| dot(row, v)).collect()
}

fn matvec_t(w: &[Vec<f64>], u: &[f64]) -> Vec<f64> {
    let n_in = w[0].len();
    let mut z = vec![0.0; n_in];
    for (row, &ui) in w.iter().zip(u) {
        for (zi, &wij) in z.iter_mut().zip(row) {
            *zi += wij * ui;
        }
    }
    z
}

fn norm(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

fn normalize(v: &mut [f64]) {
    let n = norm(v);
    if n > 1e-30 {
        for vi in v {
            *vi /= n;
        }
    }
}

/// Return `(max_abs, unit_bound)` such that the spectral norm is at most
/// `max_abs * unit_bound`, without forming that potentially overflowing
/// product. All inputs are validated here so every certificate-bearing caller
/// shares the same fail-closed boundary.
fn scaled_spectral_upper_parts(weights: &[Vec<f64>]) -> (f64, f64) {
    let Some(first) = weights.first() else {
        return (0.0, 0.0);
    };
    let n_in = first.len();
    assert!(
        weights.iter().all(|row| row.len() == n_in),
        "ragged weights"
    );
    assert!(
        weights.iter().flatten().all(|w| w.is_finite()),
        "weights must be finite"
    );

    let max_abs = weights
        .iter()
        .flatten()
        .map(|w| w.abs())
        .fold(0.0_f64, f64::max);
    if max_abs == 0.0 {
        return (0.0, 0.0);
    }

    let mut frobenius_sum = 0.0f64;
    let mut max_row_sum = 0.0f64;
    let mut column_sums = vec![0.0f64; n_in];
    for row in weights {
        let mut row_sum = 0.0;
        for (column_sum, magnitude) in column_sums.iter_mut().zip(row.iter().map(|w| w.abs())) {
            if magnitude == 0.0 {
                continue;
            }
            let ratio = if magnitude == max_abs {
                1.0
            } else {
                next_up_nonnegative(magnitude / max_abs)
            };
            frobenius_sum = add_up_nonnegative(frobenius_sum, mul_up_nonnegative(ratio, ratio));
            row_sum = add_up_nonnegative(row_sum, ratio);
            *column_sum = add_up_nonnegative(*column_sum, ratio);
        }
        max_row_sum = max_row_sum.max(row_sum);
    }

    let frobenius_bound = if frobenius_sum == 1.0 {
        1.0
    } else {
        next_up_nonnegative(frobenius_sum.sqrt())
    };
    let max_column_sum = column_sums.into_iter().fold(0.0_f64, f64::max);
    let induced_product = mul_up_nonnegative(max_column_sum, max_row_sum);
    let induced_bound = if induced_product == 1.0 {
        1.0
    } else {
        next_up_nonnegative(induced_product.sqrt())
    };
    let unit_bound = frobenius_bound.min(induced_bound);
    (max_abs, unit_bound)
}

fn try_spectral_norm_upper_bound(weights: &[Vec<f64>]) -> Option<f64> {
    let (max_abs, unit_bound) = scaled_spectral_upper_parts(weights);
    finite_scaled_upper(max_abs, unit_bound)
}

fn finite_scaled_upper(max_abs: f64, unit_bound: f64) -> Option<f64> {
    if max_abs == 0.0 || unit_bound == 0.0 {
        return Some(0.0);
    }
    if unit_bound == 1.0 {
        // The induced certificate is exactly max_abs (for example, a single
        // nonzero or a scaled permutation matrix). Keeping this exact also
        // admits a 1×1 matrix containing f64::MAX.
        return Some(max_abs);
    }
    let product = mul_up_nonnegative(max_abs, unit_bound);
    product.is_finite().then_some(product)
}

fn add_up_nonnegative(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        return b;
    }
    if b == 0.0 {
        return a;
    }
    next_up_nonnegative(a + b)
}

fn mul_up_nonnegative(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 {
        return 0.0;
    }
    if a == 1.0 {
        return b;
    }
    if b == 1.0 {
        return a;
    }
    next_up_nonnegative(a * b)
}

fn product_up(a: f64, b: f64) -> f64 {
    mul_up_nonnegative(a, b)
}

fn next_up_nonnegative(value: f64) -> f64 {
    if value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    f64::from_bits(value.to_bits() + 1)
}

fn next_down_nonnegative(value: f64) -> f64 {
    if value <= 0.0 {
        return 0.0;
    }
    if value == f64::INFINITY {
        return f64::MAX;
    }
    f64::from_bits(value.to_bits() - 1)
}

fn zero_weights(weights: &mut [Vec<f64>]) {
    for row in weights {
        row.fill(0.0);
    }
}
