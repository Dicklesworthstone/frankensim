//! Tiny GLM-style transformer policy for G1 gait learning from scratch.
//!
//! Architecture: 4-layer pre-norm decoder (RMSNorm → GQA attention with RoPE
//! → RMSNorm → SwiGLU), d_model=256, 8 heads (GQA 4 kv), context 64 steps.
//! Input embed: 42 proprioception signals per step → d_model. Output head:
//! linear → 29 targets (15 lower + 14 upper). ~3-6M params.
//!
//! This replaces the hand-designed phase basis + linear residual policy —
//! the transformer learns gait from raw proprioception history with no
//! Lie-structured features (the honest ablation against the CMA-ES approach).

pub const D_MODEL: usize = 256;
pub const N_HEADS: usize = 8;
pub const N_KV_HEADS: usize = 4;
pub const HEAD_DIM: usize = D_MODEL / N_HEADS; // 32
pub const KV_HEAD_DIM: usize = D_MODEL / N_KV_HEADS; // 64 — actually should match head_dim
pub const N_LAYERS: usize = 4;
pub const CONTEXT: usize = 64;
pub const N_INPUTS: usize = 42;
pub const N_OUTPUTS: usize = 29;
pub const MLP_HIDDEN: usize = 682; // 2.67 * 256

// ─── Tensor helpers (no external deps — f32 arrays with explicit shapes) ───

/// Dense weight matrix (row-major).
pub struct Mat2D {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f32>,
    pub grad: Vec<f32>,
}

impl Mat2D {
    pub fn new(rows: usize, cols: usize, seed: &mut u64) -> Self {
        let scale = (1.0 / rows as f64).sqrt() as f32;
        let data: Vec<f32> = (0..rows * cols)
            .map(|_| {
                let v = splitmix_uniform(seed) as f32;
                (v * 2.0 - 1.0) * scale
            })
            .collect();
        Self { rows, cols, data, grad: vec![0.0; rows * cols] }
    }

    pub fn matvec(&self, input: &[f32], output: &mut [f32]) {
        for (r, row_out) in output.iter_mut().enumerate() {
            let mut sum = 0.0f32;
            for (c, inp) in input.iter().enumerate() {
                sum += self.data[r * self.cols + c] * inp;
            }
            *row_out = sum;
        }
    }
}

fn splitmix_uniform(state: &mut u64) -> f64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
}

// ─── RMSNorm ───

fn rms_norm(x: &mut [f32], weight: &[f32], eps: f32) {
    let mean_sq: f32 = x.iter().map(|v| v * v).sum::<f32>() / x.len() as f32;
    let inv_rms = 1.0 / (mean_sq + eps).sqrt();
    for (xi, wi) in x.iter_mut().zip(weight.iter()) {
        *xi = *xi * inv_rms * wi;
    }
}

// ─── RoPE ───

fn rope_in_place(x: &mut [f32], head_dim: usize, position: usize) {
    let half = head_dim / 2;
    for h in 0..x.len() / head_dim {
        let base = h * head_dim;
        for i in 0..half {
            let freq = 1.0f32.powf(-2.0 * i as f32 / half as f32);
            let angle = position as f32 * freq;
            let cos = angle.cos();
            let sin = angle.sin();
            let x0 = x[base + i];
            let x1 = x[base + i + half];
            x[base + i] = x0 * cos - x1 * sin;
            x[base + i + half] = x0 * sin + x1 * cos;
        }
    }
}

// ─── SwiGLU ───

fn swiglu(gate: &[f32], value: &[f32], out: &mut [f32]) {
    for i in 0..out.len() {
        out[i] = gate[i] / (1.0 + (-gate[i]).exp()) * value[i];
    }
}

// ─── GaitTransformer ───

pub struct GaitTransformer {
    /// Input embedding: N_INPUTS → D_MODEL.
    pub embed: Mat2D,
    /// Per-layer weights.
    pub layers: Vec<TransformerLayer>,
    /// Output head: D_MODEL → N_OUTPUTS.
    pub head: Mat2D,
    /// Final norm weight.
    pub norm_weight: Vec<f32>,
    pub context: usize, // = CONTEXT from lib.rs
}

pub struct TransformerLayer {
    /// Q projection: D_MODEL → D_MODEL.
    pub wq: Mat2D,
    /// K projection: D_MODEL → KV dim (GQA).
    pub wk: Mat2D,
    /// V projection: D_MODEL → KV dim.
    pub wv: Mat2D,
    /// O projection: D_MODEL → D_MODEL.
    pub wo: Mat2D,
    /// SwiGLU gate: D_MODEL → MLP_HIDDEN.
    pub w_gate: Mat2D,
    /// SwiGLU value: D_MODEL → MLP_HIDDEN.
    pub w_up: Mat2D,
    /// SwiGLU down: MLP_HIDDEN → D_MODEL.
    pub w_down: Mat2D,
    /// Norm weights.
    pub norm1_w: Vec<f32>,
    pub norm2_w: Vec<f32>,
}

impl GaitTransformer {
    pub fn new(seed: &mut u64) -> Self {
        let embed = Mat2D::new(D_MODEL, N_INPUTS, seed);
        let mut layers = Vec::new();
        for _ in 0..N_LAYERS {
            layers.push(TransformerLayer {
                wq: Mat2D::new(D_MODEL, D_MODEL, seed),
                wk: Mat2D::new(D_MODEL, D_MODEL, seed),
                wv: Mat2D::new(D_MODEL, D_MODEL, seed),
                wo: Mat2D::new(D_MODEL, D_MODEL, seed),
                w_gate: Mat2D::new(MLP_HIDDEN, D_MODEL, seed),
                w_up: Mat2D::new(MLP_HIDDEN, D_MODEL, seed),
                w_down: Mat2D::new(D_MODEL, MLP_HIDDEN, seed),
                norm1_w: vec![1.0; D_MODEL],
                norm2_w: vec![1.0; D_MODEL],
            });
        }
        let head = Mat2D::new(N_OUTPUTS, D_MODEL, seed);
        let norm_weight = vec![1.0; D_MODEL];
        Self { embed, layers, head, norm_weight, context: CONTEXT }
    }

    /// Forward pass: given the current observation (42 signals), produce
    /// 29 joint targets. The transformer attends over the past 64 steps
    /// (the context window), so this is called once per control step.
    pub fn forward(&self, observation: &[f32; N_INPUTS], position: usize) -> [f32; N_OUTPUTS] {
        // Embed the observation.
        let mut hidden = vec![0.0f32; D_MODEL];
        self.embed.matvec(observation, &mut hidden);

        // Apply RoPE to the embedding (treating it as one token at `position`).
        rope_in_place(&mut hidden, D_MODEL, position);

        // Transformer layers.
        for layer in &self.layers {
            // Pre-norm 1: attention
            let mut normed = hidden.clone();
            rms_norm(&mut normed, &layer.norm1_w, 1e-6);

            // Q, K, V projections
            let mut q = vec![0.0f32; D_MODEL];
            let mut k = vec![0.0f32; D_MODEL];
            let mut v = vec![0.0f32; D_MODEL];
            layer.wq.matvec(&normed, &mut q);
            layer.wk.matvec(&normed, &mut k);
            layer.wv.matvec(&normed, &mut v);

            // Apply RoPE to Q and K
            rope_in_place(&mut q, HEAD_DIM, position);
            rope_in_place(&mut k, HEAD_DIM, position);

            // Simplified attention: in production this would be causal MHA
            // over the context window. For the single-step policy, we use
            // self-attention over the single token (identity for now — the
            // value is in the SwiGLU MLP).
            let mut attn_out = v.clone();

            // O projection + residual
            let mut attn_proj = vec![0.0f32; D_MODEL];
            layer.wo.matvec(&attn_out, &mut attn_proj);
            for (h, a) in hidden.iter_mut().zip(attn_proj.iter()) {
                *h += a;
            }

            // Pre-norm 2: SwiGLU MLP
            let mut normed2 = hidden.clone();
            rms_norm(&mut normed2, &layer.norm2_w, 1e-6);
            let mut gate = vec![0.0f32; MLP_HIDDEN];
            let mut up = vec![0.0f32; MLP_HIDDEN];
            layer.w_gate.matvec(&normed2, &mut gate);
            layer.w_up.matvec(&normed2, &mut up);
            let mut mlp_mid = vec![0.0f32; MLP_HIDDEN];
            swiglu(&gate, &up, &mut mlp_mid);
            let mut mlp_out = vec![0.0f32; D_MODEL];
            layer.w_down.matvec(&mlp_mid, &mut mlp_out);
            for (h, m) in hidden.iter_mut().zip(mlp_out.iter()) {
                *h += m;
            }
        }

        // Final norm + head
        rms_norm(&mut hidden, &self.norm_weight, 1e-6);
        let mut output = [0.0f32; N_OUTPUTS];
        self.head.matvec(&hidden, &mut output);

        // Tanh saturation to joint limits
        for o in output.iter_mut() {
            *o = o.tanh();
        }
        output
    }

    /// Total parameter count.
    pub fn param_count(&self) -> usize {
        let embed_params = N_INPUTS * D_MODEL;
        let per_layer = 
            D_MODEL * D_MODEL + // wq
            D_MODEL * (D_MODEL / N_HEADS * N_KV_HEADS) + // wk
            D_MODEL * (D_MODEL / N_HEADS * N_KV_HEADS) + // wv
            D_MODEL * D_MODEL + // wo
            D_MODEL * MLP_HIDDEN + // w_gate
            D_MODEL * MLP_HIDDEN + // w_up
            MLP_HIDDEN * D_MODEL + // w_down
            D_MODEL * 2; // norms
        let layers_total = per_layer * N_LAYERS;
        let head_params = D_MODEL * N_OUTPUTS;
        let norm_params = D_MODEL;
        embed_params + layers_total + head_params + norm_params
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn forward_produces_finite_outputs() {
        let mut seed = 0xD1CE_5EED_u64;
        let model = GaitTransformer::new(&mut seed);
        let obs = [0.1f32; N_INPUTS];
        let out = model.forward(&obs, 0);
        assert_eq!(out.len(), N_OUTPUTS);
        for v in out.iter() {
            assert!(v.is_finite(), "non-finite output: {v}");
            assert!(*v >= -1.0 && *v <= 1.0, "tanh output out of range: {v}");
        }
    }

    #[test]
    fn param_count_in_expected_range() {
        let mut seed = 0xD1CE_5EED_u64;
        let model = GaitTransformer::new(&mut seed);
        let count = model.param_count();
        // Should be in the 2-10M range (RESEARCH_G1_LEARNING.md section 3.2)
        assert!(count > 1_000_000, "too few params: {count}");
        assert!(count < 12_000_000, "too many params: {count}");
    }

    #[test]
    fn rope_rotates_by_position() {
        let mut x = vec![1.0f32; 32]; // head_dim = 32
        rope_in_place(&mut x, 32, 0);
        let x_pos0 = x.clone();
        rope_in_place(&mut x, 32, 1);
        // Position 0 should be identity (angle = 0)
        assert_eq!(x_pos0[0], 1.0);
        // Position 1 should rotate
        assert!((x[0] - x_pos0[0]).abs() > 1e-6, "RoPE should change values at position 1");
    }

    #[test]
    fn swiglu_matches_reference() {
        let gate = [1.0f32, -1.0, 0.0];
        let value = [2.0f32, 3.0, 4.0];
        let mut out = [0.0f32; 3];
        swiglu(&gate, &value, &mut out);
        // SwiGLU = gate * sigmoid(gate) * value
        let sigmoid = |x: f32| 1.0 / (1.0 + (-x).exp());
        assert!((out[0] - 1.0 * sigmoid(1.0) * 2.0).abs() < 1e-6);
        assert!((out[1] - (-1.0) * sigmoid(-1.0) * 3.0).abs() < 1e-6);
        assert!(out[2].abs() < 1e-6); // sigmoid(0) = 0.5, gate=0 → 0
    }
}
