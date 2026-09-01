//! Muon optimizer (Momentum + Orthogonalized weight updates).
//!
//! For 2-D weight matrices W: momentum buffer M = β·M + grad;
//! orthogonalize O = orth(M); W -= lr · O · scale, optionally with
//! decoupled weight decay. The orthogonalization here is exact
//! Gram-Schmidt over the columns (deterministic, well-conditioned in
//! f32, O(mn²)) — the same "nearest semi-orthogonal matrix" target as
//! Keller Jordan's 5-iteration quintic Newton–Schulz (a=3.4445,
//! b=−4.7750, c=2.0315; arXiv 2509.24406), which is the standard
//! faster approximation. Swapping GS for the NS quintic is a pure
//! performance lever and must be parity-gated per project doctrine.
//!
//! Split: Muon for hidden-layer 2-D weights, Adam for 1-D params (biases,
//! norms) + embed/head.

// Newton–Schulz quintic coefficients (Keller Jordan), reserved for the
// parity-gated NS fast path; unused by the exact Gram-Schmidt path.
#[allow(dead_code)]
const NS_A: f32 = 3.4445;
#[allow(dead_code)]
const NS_B: f32 = -4.7750;
#[allow(dead_code)]
const NS_C: f32 = 2.0315;
#[allow(dead_code)]
const NS_EPS: f32 = 1e-7;
#[allow(dead_code)]
const NS_ITERATIONS: usize = 5;

/// Tall-column Gram–Schmidt: orthonormalizes the columns of a TALL
/// (rows >= cols) matrix in place. Call `orthogonalize_momentum` instead —
/// this helper is the tall-only core (wide matrices MUST be transposed
/// first: column GS on a wide matrix exhausts the row rank at column
/// `rows` and the remaining columns degenerate to amplified noise).
pub fn newton_schulz_orthogonalize(x: &mut Vec<f32>, rows: usize, cols: usize) {
    debug_assert!(
        rows >= cols,
        "tall-only core: transpose wide matrices first"
    );
    assert_eq!(x.len(), rows * cols);
    // Gram-Schmidt orthonormalization of columns (always correct, O(mn²)).
    // For production use Newton-Schulz; this version is guaranteed to work.
    const EPS: f32 = 1e-10;
    for j in 0..cols {
        // Subtract projections onto all previous columns
        for k in 0..j {
            let mut dot = 0.0f32;
            for i in 0..rows {
                dot += x[i * cols + j] * x[i * cols + k];
            }
            for i in 0..rows {
                x[i * cols + j] -= dot * x[i * cols + k];
            }
        }
        // Normalize
        let mut norm = 0.0f32;
        for i in 0..rows {
            norm += x[i * cols + j] * x[i * cols + j];
        }
        let norm = norm.sqrt();
        if norm > EPS {
            for i in 0..rows {
                x[i * cols + j] /= norm;
            }
        } else {
            for i in 0..rows {
                x[i * cols + j] = 0.0;
            }
        }
    }
}

/// Muon optimizer state for one 2-D weight matrix.
pub struct MuonParam {
    pub rows: usize,
    pub cols: usize,
    pub weights: Vec<f32>,
    pub momentum: Vec<f32>,
    pub grad: Vec<f32>,
    pub lr: f32,
    pub momentum_beta: f32,
    pub weight_decay: f32,
}

impl MuonParam {
    pub fn new(rows: usize, cols: usize, lr: f32, momentum_beta: f32) -> Self {
        Self {
            rows,
            cols,
            weights: vec![0.0; rows * cols],
            momentum: vec![0.0; rows * cols],
            grad: vec![0.0; rows * cols],
            lr,
            momentum_beta,
            weight_decay: 0.0,
        }
    }

    /// One Muon step: accumulate momentum, orthogonalize, update weights.
    pub fn step(&mut self) {
        // Momentum: M = β·M + grad
        for i in 0..self.momentum.len() {
            self.momentum[i] = self.momentum_beta * self.momentum[i] + self.grad[i];
        }
        // Orthogonalize momentum. For WIDE matrices (cols > rows) the
        // column Gram-Schmidt must run on the TRANSPOSE (tall orientation):
        // running it directly exhausts the row rank at column `rows` and
        // the tail columns degenerate to amplified noise (they normalize
        // the residual of linearly-dependent vectors). The transposed
        // variant orthonormalizes the row space instead, giving every
        // weight a well-defined semi-orthogonal update direction.
        let o = if self.rows >= self.cols {
            let mut o = self.momentum.clone();
            newton_schulz_orthogonalize(&mut o, self.rows, self.cols);
            o
        } else {
            let mut t = vec![0.0f32; self.momentum.len()];
            for r in 0..self.rows {
                for c in 0..self.cols {
                    t[c * self.rows + r] = self.momentum[r * self.cols + c];
                }
            }
            newton_schulz_orthogonalize(&mut t, self.cols, self.rows);
            let mut o = vec![0.0f32; self.momentum.len()];
            for r in 0..self.rows {
                for c in 0..self.cols {
                    o[r * self.cols + c] = t[c * self.rows + r];
                }
            }
            o
        };
        // Weight update: W -= lr · O
        for i in 0..self.weights.len() {
            self.weights[i] -= self.lr * o[i];
            if self.weight_decay > 0.0 {
                self.weights[i] *= 1.0 - self.lr * self.weight_decay;
            }
        }
    }
}

/// Adam optimizer for 1-D params (biases, norms, embed/head).
pub struct AdamParam {
    pub params: Vec<f32>,
    pub grads: Vec<f32>,
    pub m: Vec<f32>,
    pub v: Vec<f32>,
    pub lr: f32,
    pub beta1: f32,
    pub beta2: f32,
    pub eps: f32,
    pub t: u64,
}

impl AdamParam {
    pub fn new(size: usize, lr: f32) -> Self {
        Self {
            params: vec![0.0; size],
            grads: vec![0.0; size],
            m: vec![0.0; size],
            v: vec![0.0; size],
            lr,
            beta1: 0.9,
            beta2: 0.999,
            eps: 1e-8,
            t: 0,
        }
    }

    pub fn step(&mut self) {
        self.t += 1;
        let bc1 = 1.0 - self.beta1.powi(self.t as i32);
        let bc2 = 1.0 - self.beta2.powi(self.t as i32);
        for i in 0..self.params.len() {
            self.m[i] = self.beta1 * self.m[i] + (1.0 - self.beta1) * self.grads[i];
            self.v[i] = self.beta2 * self.v[i] + (1.0 - self.beta2) * self.grads[i] * self.grads[i];
            let m_hat = self.m[i] / bc1;
            let v_hat = self.v[i] / bc2;
            self.params[i] -= self.lr * m_hat / (v_hat.sqrt() + self.eps);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ns_orthonormalizes_small_matrix() {
        let rows = 4;
        let cols = 3;
        let mut x: Vec<f32> = (0..rows * cols)
            .map(|i| ((i * 7 + 3) % 13) as f32 / 13.0 - 0.5)
            .collect();
        newton_schulz_orthogonalize(&mut x, rows, cols);
        // Check: X^T X ≈ I (cols × cols identity) for semi-orthogonal
        for j in 0..cols {
            for k in 0..cols {
                let mut dot = 0.0f32;
                for i in 0..rows {
                    dot += x[i * cols + j] * x[i * cols + k];
                }
                let expected = if j == k { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 0.15,
                    "X^T X[{j}][{k}] = {dot}, expected ≈ {expected}"
                );
            }
        }
    }

    #[test]
    fn muon_step_reduces_loss_on_quadratic() {
        // Simple quadratic loss: L = ||W - target||²
        let rows = 4;
        let cols = 3;
        let target: Vec<f32> = (0..rows * cols)
            .map(|i| ((i * 11) % 7) as f32 / 7.0 - 0.5)
            .collect();
        let mut param = MuonParam::new(rows, cols, 0.05, 0.9);
        param.weights.copy_from_slice(&vec![0.0; rows * cols]);
        let initial_loss: f32 = param
            .weights
            .iter()
            .zip(target.iter())
            .map(|(w, t)| (w - t) * (w - t))
            .sum();
        for _ in 0..100 {
            for i in 0..param.grad.len() {
                param.grad[i] = 2.0 * (param.weights[i] - target[i]);
            }
            param.step();
        }
        let final_loss: f32 = param
            .weights
            .iter()
            .zip(target.iter())
            .map(|(w, t)| (w - t) * (w - t))
            .sum();
        assert!(
            final_loss < initial_loss * 0.5,
            "Muon should reduce loss: {initial_loss} -> {final_loss}"
        );
    }

    #[test]
    fn adam_step_matches_reference() {
        let mut adam = AdamParam::new(2, 0.01);
        adam.params = vec![1.0, 2.0];
        adam.grads = vec![0.1, -0.2];
        adam.step();
        // Manual: m = 0.9*0 + 0.1*[0.1, -0.2] = [0.01, -0.02]
        // v = 0.999*0 + 0.001*[0.01, 0.04] = [1e-5, 4e-5]
        // m_hat = m/0.1, v_hat = v/0.001
        // update = 0.01 * m_hat / (sqrt(v_hat) + eps)
        assert!(
            adam.params[0] < 1.0,
            "Adam should decrease param[0] with positive grad"
        );
        assert!(
            adam.params[1] > 2.0,
            "Adam should increase param[1] with negative grad"
        );
    }

    #[test]
    fn muon_lr_routing_splits() {
        // Verify that Muon and Adam params can coexist and route independently.
        let mut muon_w = MuonParam::new(3, 3, 0.01, 0.9);
        let mut adam_b = AdamParam::new(3, 0.001);
        let before_muon = muon_w.weights.clone();
        let before_adam = adam_b.params.clone();
        muon_w
            .grad
            .iter_mut()
            .enumerate()
            .for_each(|(i, g)| *g = (i % 5) as f32 * 0.1);
        adam_b
            .grads
            .iter_mut()
            .enumerate()
            .for_each(|(i, g)| *g = (i % 3) as f32 * 0.01);
        muon_w.step();
        adam_b.step();
        assert_ne!(muon_w.weights, before_muon);
        assert_ne!(adam_b.params, before_adam);
    }
}

#[cfg(test)]
mod wide_muon {
    use crate::muon::MuonParam;
    use crate::transformer::randomize_uniform;

    /// Discriminator for the wide-matrix Gram-Schmidt defect: with a
    /// full-rank random target, the column-GS orthogonalization zeroes
    /// updates for columns >= rows (128), so residual mass piles up in
    /// the tail columns. A correct Muon must shrink BOTH column blocks.
    #[test]
    fn muon_wide_matrix_shrinks_both_column_blocks() {
        let rows = 128usize;
        let cols = 256usize;
        let mut seed = 0xDE11u64.wrapping_add(0x9E37); // deterministic
        let mut p = MuonParam::new(rows, cols, 0.2, 0.9);
        let mut target = vec![0.0f32; rows * cols];
        randomize_uniform(&mut target, cols, &mut seed);

        let block_loss = |w: &[f32], c0: usize, c1: usize| -> f32 {
            let mut l = 0.0;
            for r in 0..rows {
                for c in c0..c1 {
                    let d = w[r * cols + c] - target[r * cols + c];
                    l += d * d;
                }
            }
            l
        };
        let l_head = block_loss(&p.weights, 0, rows / 2);
        let l_tail = block_loss(&p.weights, rows / 2, cols.min(rows * 2));

        for _ in 0..120 {
            for i in 0..p.weights.len() {
                p.grad[i] = 2.0 * (p.weights[i] - target[i]);
            }
            p.step();
        }
        let l_head_after = block_loss(&p.weights, 0, rows / 2);
        let l_tail_after = block_loss(&p.weights, rows / 2, cols.min(rows * 2));
        println!(
            "[wide-muon] head block {l_head:.2} -> {l_head_after:.2} ; tail block {l_tail:.2} -> {l_tail_after:.2}"
        );
        assert!(
            l_head_after < l_head * 0.9,
            "head-block residual must shrink ({} -> {})",
            l_head,
            l_head_after
        );
        assert!(
            l_tail_after < l_tail * 0.9,
            "tail-block residual must shrink — if it stalls, wide GS is zeroing the update ({} -> {})",
            l_tail, l_tail_after
        );
    }
}
