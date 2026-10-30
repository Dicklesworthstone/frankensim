//! Preserve the positive factor of the physical resistance integral.
//!
//! Expanding W into z^T W z can erase the heat of a tiny vibration about a
//! force-supported equilibrium. We instead retain W = L^T L and evaluate
//! ||L z||^2, using initial loss effort and the actual force imbalance as
//! coordinates. Positive Gauss weights on a norm-bounded short interval and
//! QR doubling preserve this factor without clipping heat or resistance modes.
use super::{LinearSplit, checked_eigh, exponential, finite, product, refuse};
use crate::{PhsError, PortHamiltonian};
use fs_la::factor::qr;
use fs_math::{det, eft::two_sum};

pub(super) struct PositiveHeat {
    n: usize,
    m: usize,
    d: usize,
    forward: Vec<f64>,
    resistance_root: Vec<f64>,
    rate: f64,
    root: Vec<f64>,
}

impl PositiveHeat {
    /// Same physical generator and energy transform as PreparedQuadraticFlow.
    /// The caller has already admitted the complete state/port capacity.
    pub(super) fn new(
        sys: &PortHamiltonian,
        split: &LinearSplit,
        dt: f64,
        resistance: &[f64],
    ) -> Result<Self, PhsError> {
        let (n, m) = (split.n, split.m);
        let d = 2 * n;
        if resistance.len() != n * n || !finite(resistance) || !dt.is_finite() || dt <= 0.0 {
            return Err(refuse(
                "positive heat requires finite resistance and a positive clock",
            ));
        }
        let mut forward = vec![0.0; n * n];
        let mut inverse = vec![0.0; n * n];
        for i in 0..n {
            let weight = if split.values[i] > 0.0 {
                det::sqrt(split.values[i])
            } else {
                1.0
            };
            for j in 0..n {
                forward[i * n + j] = weight * split.basis[j * n + i];
                inverse[j * n + i] = split.basis[j * n + i] / weight;
            }
        }
        let scale = resistance.iter().fold(0.0_f64, |s, r| s.max(r.abs()));
        if scale == 0.0 {
            return Ok(Self {
                n,
                m,
                d,
                forward,
                resistance_root: Vec::new(),
                rate: 1.0,
                root: Vec::new(),
            });
        }
        // Passive resistance graphs have an explicit positive decomposition.
        // Their exact effort nullspaces must not depend on the rounded sign of
        // a numerical eigenvalue, especially beside much stronger loss blocks.
        let resistance_root = factor_resistance(resistance, n)?;
        let mut original = vec![0.0; n * n];
        for i in 0..n {
            for k in 0..n {
                for j in 0..n {
                    original[i * n + j] +=
                        (sys.j[i * n + k] - sys.r[i * n + k]) * split.quadratic.matrix[k * n + j];
                }
            }
        }
        let balanced = product(&product(&forward, &original, n), &inverse, n);
        let effort = product(&split.quadratic.matrix, &inverse, n);
        let loss_map = product(&resistance_root, &effort, n);
        let mut rate = 1.0_f64;
        for i in 0..n {
            rate = rate
                .max((0..n).map(|j| balanced[i * n + j].abs()).sum())
                .max((0..n).map(|j| balanced[j * n + i].abs()).sum());
        }
        if !rate.is_finite() {
            return Err(refuse("positive heat coordinate rate is nonfinite"));
        }
        let root_rate = det::sqrt(rate);
        let mut generator = vec![0.0; d * d];
        let mut factor = vec![0.0; d * d];
        for i in 0..n {
            factor[i * d + i] = root_rate;
            for j in 0..n {
                generator[i * d + n + j] = root_rate * loss_map[i * n + j];
                generator[(n + i) * d + n + j] = balanced[i * n + j];
            }
        }
        // For held u, f=x' obeys f'=A f, even when A is singular. With
        // a=sqrt(R) grad(H)/sqrt(rate), b=T f/rate, a'=sqrt(rate) L b
        // and b'=T A T^-1 b. No static solve or equilibrium shift is needed.
        // This retains the complete trajectory but forms source/internal-force
        // cancellation BEFORE the heat square, keeping static modes silent.
        let root = integral_root(&generator, &factor, d, dt)?;
        Ok(Self {
            n,
            m,
            d,
            forward,
            resistance_root,
            rate,
            root,
        })
    }

    pub(super) fn evaluate(
        &self,
        sys: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
    ) -> Result<f64, PhsError> {
        self.evaluate_into(
            sys,
            state,
            input,
            &mut vec![0.0; self.n],
            &mut vec![0.0; self.n],
            &mut vec![0.0; self.d],
        )
    }

    /// No allocation, exponential, equilibrium solve or expanded heat form.
    pub(super) fn evaluate_into(
        &self,
        sys: &PortHamiltonian,
        state: &[f64],
        input: &[f64],
        effort: &mut [f64],
        derivative: &mut [f64],
        coordinates: &mut [f64],
    ) -> Result<f64, PhsError> {
        let n = self.n;
        if state.len() != n
            || input.len() != self.m
            || effort.len() != n
            || derivative.len() != n
            || coordinates.len() != self.d
            || !finite(state)
            || !finite(input)
        {
            return Err(refuse(
                "positive heat requires complete finite state/input/scratch",
            ));
        }
        if self.root.is_empty() {
            return Ok(0.0);
        }
        sys.storage.gradient(state, effort);
        let root_rate = det::sqrt(self.rate);
        for i in 0..n {
            let force: f64 = (0..self.m).map(|p| sys.g[i * self.m + p] * input[p]).sum();
            let conservative: f64 = (0..n).map(|j| sys.j[i * n + j] * effort[j]).sum();
            let resistive: f64 = (0..n).map(|j| sys.r[i * n + j] * effort[j]).sum();
            derivative[i] = (force + conservative) - resistive;
            coordinates[i] = (0..n)
                .map(|j| self.resistance_root[i * n + j] * effort[j])
                .sum::<f64>()
                / root_rate;
        }
        for i in 0..n {
            coordinates[n + i] = (0..n)
                .map(|j| self.forward[i * n + j] * derivative[j])
                .sum::<f64>()
                / self.rate;
        }
        if !finite(coordinates) {
            return Err(refuse("positive heat initial effort/rate is nonfinite"));
        }
        let mut value = 0.0;
        for row in self.root.chunks_exact(self.d) {
            let loss: f64 = row.iter().zip(coordinates.iter()).map(|(a, b)| a * b).sum();
            value += loss * loss;
        }
        if !value.is_finite() {
            return Err(refuse("physical resistance heat is unrepresentable"));
        }
        Ok(value)
    }
}

// Exact zero entries separate independent resistance components. A general
// full-C material block must not force a Foster graph's nullspace through an
// eigensolve, or set the absolute scale of a physically weak thermal block.
fn factor_resistance(resistance: &[f64], n: usize) -> Result<Vec<f64>, PhsError> {
    let mut root = vec![0.0; n * n];
    let mut seen = vec![false; n];
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut indices = vec![start];
        seen[start] = true;
        let mut next = 0;
        while next < indices.len() {
            let i = indices[next];
            next += 1;
            for j in 0..n {
                if !seen[j] && (resistance[i * n + j] != 0.0 || resistance[j * n + i] != 0.0) {
                    indices.push(j);
                    seen[j] = true;
                }
            }
        }
        indices.sort_unstable();
        let count = indices.len();
        let matrix: Vec<_> = indices
            .iter()
            .flat_map(|&i| indices.iter().map(move |&j| resistance[i * n + j]))
            .collect();
        let factor = if let Some(factor) = graph_root(&matrix, count)? {
            factor
        } else {
            let scale = matrix.iter().fold(0.0_f64, |s, r| s.max(r.abs()));
            let normalized: Vec<_> = matrix.iter().map(|r| r / scale).collect();
            let (values, vectors) = checked_eigh(&normalized, count)?;
            if values.iter().any(|&value| value < 0.0) {
                return Err(refuse(
                    "physical heat resistance has an unresolved negative eigenvalue",
                ));
            }
            let mut factor = vec![0.0; count * count];
            for i in 0..count {
                let weight = det::sqrt(values[i]) * det::sqrt(scale);
                for j in 0..count {
                    factor[i * count + j] = weight * vectors[j * count + i];
                }
            }
            factor
        };
        check_resistance_root(&matrix, &factor, count)?;
        for (i, &row) in indices.iter().enumerate() {
            for (j, &column) in indices.iter().enumerate() {
                root[row * n + column] = factor[i * count + j];
            }
        }
    }
    Ok(root)
}

// Exact signed diagonal dominance supplies physical edge/ground squares:
// |R_ij| (e_i + sign(R_ij) e_j)^2 and (R_ii-sum_j |R_ij|) e_i^2.
// Every positive coefficient survives at its own absolute scale. No tolerance
// selects rank, rescues a negative remainder or changes the original R.
fn graph_root(resistance: &[f64], n: usize) -> Result<Option<Vec<f64>>, PhsError> {
    let mut grounds = vec![0.0; n];
    let mut count = 0;
    for i in 0..n {
        for j in 0..n {
            if resistance[i * n + j] != resistance[j * n + i] {
                return Ok(None);
            }
            if i != j {
                if i < j && resistance[i * n + j] != 0.0 {
                    count += 1;
                }
            }
        }
        let Some(ground) = nonnegative_ground(resistance, n, i) else {
            return Ok(None);
        };
        grounds[i] = ground;
        if grounds[i] > 0.0 {
            count += 1;
        }
    }
    // Forests and diagonal damping already fit in n rows. Denser passive
    // graphs use the existing QR owner to retain every edge in n coordinates.
    let mut rows = vec![0.0; count.max(n) * n];
    let mut row = 0;
    for i in 0..n {
        if grounds[i] > 0.0 {
            rows[row * n + i] = det::sqrt(grounds[i]);
            row += 1;
        }
        for j in i + 1..n {
            let value = resistance[i * n + j];
            if value != 0.0 {
                let weight = det::sqrt(value.abs());
                rows[row * n + i] = weight;
                rows[row * n + j] = if value < 0.0 { -weight } else { weight };
                row += 1;
            }
        }
    }
    let root = if count > n {
        compress(&rows, count, n)?
    } else {
        rows
    };
    Ok(Some(root))
}

// Grow an exact floating-point expansion of R_ii-sum_j |R_ij|. A rounded
// off-diagonal sum can otherwise conceal a negative ground below one ulp of
// the diagonal. Nonoverlapping expansion components are ordered by magnitude,
// so the last nonzero component supplies the exact represented residual sign.
fn nonnegative_ground(resistance: &[f64], n: usize, i: usize) -> Option<f64> {
    let mut expansion = vec![resistance[i * n + i]];
    for j in 0..n {
        if i == j || resistance[i * n + j] == 0.0 {
            continue;
        }
        let mut value = -resistance[i * n + j].abs();
        let mut next = Vec::with_capacity(expansion.len() + 1);
        for &part in &expansion {
            let (sum, error) = two_sum(value, part);
            if !sum.is_finite() || !error.is_finite() {
                return None;
            }
            if error != 0.0 {
                next.push(error);
            }
            value = sum;
        }
        if value != 0.0 {
            next.push(value);
        }
        expansion = next;
    }
    if expansion.last().is_some_and(|&value| value < 0.0) {
        return None;
    }
    let ground: f64 = expansion.iter().sum();
    (ground.is_finite() && ground >= 0.0).then_some(ground)
}

// Both graph and eigen factors must resolve each represented nonzero entry;
// a large unrelated diagonal or cancellation inside QR cannot hide its loss.
// At structural zeros only, use the local product sum as the roundoff scale.
fn check_resistance_root(resistance: &[f64], root: &[f64], n: usize) -> Result<(), PhsError> {
    for i in 0..n {
        for j in i..n {
            let mut gram = 0.0;
            let mut products = 0.0;
            for k in 0..n {
                let product = root[k * n + i] * root[k * n + j];
                gram += product;
                products += product.abs();
            }
            let magnitude = if resistance[i * n + j] != 0.0 {
                resistance[i * n + j].abs()
            } else {
                products
            };
            if !gram.is_finite()
                || !magnitude.is_finite()
                || (gram - resistance[i * n + j]).abs()
                    > 128.0 * n as f64 * f64::EPSILON * magnitude
            {
                return Err(refuse(
                    "physical resistance factor exceeds componentwise roundoff",
                ));
            }
        }
    }
    Ok(())
}

// QR's factors retain exact zero rows/columns without asking for numerical
// rank. A common scale protects its squared column norms from tiny absolute
// loss magnitudes; it does not discard small directions relative to others.
fn compress(rows: &[f64], height: usize, n: usize) -> Result<Vec<f64>, PhsError> {
    let scale = rows.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if !scale.is_finite() {
        return Err(refuse("positive heat factor is nonfinite"));
    }
    if scale == 0.0 {
        return Ok(vec![0.0; n * n]);
    }
    let normalized: Vec<_> = rows.iter().map(|x| x / scale).collect();
    let decomposition = qr(&normalized, height, n);
    let mut root = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            root[i * n + j] = scale * decomposition.r(i, j);
        }
    }
    if !finite(&root) {
        return Err(refuse("positive heat QR factor is nonfinite"));
    }
    Ok(root)
}

fn integral_root(
    generator: &[f64],
    factor: &[f64],
    n: usize,
    dt: f64,
) -> Result<Vec<f64>, PhsError> {
    let mut norm = 0.0_f64;
    for i in 0..n {
        let row: f64 = (0..n).map(|j| (dt * generator[i * n + j]).abs()).sum();
        let column: f64 = (0..n).map(|j| (dt * generator[j * n + i]).abs()).sum();
        norm = norm.max(row).max(column);
    }
    if !norm.is_finite() || !finite(factor) {
        return Err(refuse("positive heat operator is nonfinite"));
    }
    let mut scale = 1.0;
    let mut doublings = 0;
    while norm * scale > 0.25 {
        if doublings == 48 {
            return Err(refuse("positive heat scaling exceeds 48 doublings"));
        }
        scale *= 0.5;
        doublings += 1;
    }
    let h = dt * scale;
    let a: Vec<_> = generator.iter().map(|x| h * x).collect();
    let mut transition = exponential(&a, n)?;
    // Eight-point Gauss on ||A h|| <= 1/4: the first omitted degree is 16.
    // The exponential product's truncation is below f64 roundoff at this norm;
    // quadrature weights remain positive, unlike an expanded Lyapunov series.
    const NODES: [f64; 8] = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.525532409916329,
        -0.1834346424956498,
        0.1834346424956498,
        0.525532409916329,
        0.7966664774136267,
        0.9602898564975363,
    ];
    const WEIGHTS: [f64; 8] = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.362683783378362,
        0.362683783378362,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ];
    let mut samples = vec![0.0; 8 * n * n];
    for k in 0..8 {
        let node = 0.5 * (1.0 + NODES[k]);
        let at: Vec<_> = a.iter().map(|x| node * x).collect();
        let map = exponential(&at, n)?;
        let output = product(factor, &map, n);
        let weight = det::sqrt(h) * det::sqrt(0.5 * WEIGHTS[k]);
        for (target, value) in samples[k * n * n..(k + 1) * n * n].iter_mut().zip(output) {
            *target = weight * value;
        }
    }
    let mut root = compress(&samples, 8 * n, n)?;
    for _ in 0..doublings {
        // W(2h)=W(h)+E(h)^T W(h) E(h): QR of [L; L E], never a
        // Gram matrix or exp(-A h). Stiff passive relaxation stays bounded.
        let advanced = product(&root, &transition, n);
        let mut stacked = Vec::with_capacity(2 * n * n);
        stacked.extend_from_slice(&root);
        stacked.extend_from_slice(&advanced);
        root = compress(&stacked, 2 * n, n)?;
        transition = product(&transition, &transition, n);
        if !finite(&transition) {
            return Err(refuse("positive heat transition overflows"));
        }
    }
    Ok(root)
}
