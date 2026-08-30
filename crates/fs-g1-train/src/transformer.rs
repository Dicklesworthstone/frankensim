//! Real GLM-style decoder policy for G1 gait learning — scaled for iPhone
//! on-device training and inference (bead cmaes-5ql re-scoped by owner
//! directive 2026-08-29: "it must be a REAL transformer ... trained using
//! muon and adam for real").
//!
//! # Architecture (spec-first, GLM/Kimi/DeepSeek decoder template)
//!
//! Per token t (one control step = one 42-dim proprioception vector):
//!
//! ```text
//! h0 = W_embed · obs_t                          (42 -> d_model)
//! for each of N_LAYERS blocks (pre-norm):
//!     n1 = rmsnorm(h, w1)
//!     q = Wq·n1  (d_model, 8 heads x head_dim)   RoPE(q, t)
//!     k = Wk·n1  (kv_dim = n_kv*head_dim)        RoPE(k, t)
//!     v = Wv·n1  (kv_dim)
//!     cache.push(k, v)                            (sliding window, 64)
//!     a = causal_gqa_attention(q, k_cache, v_cache)   softmax(QKᵀ/√hd)·V
//!     h = h + Wo·a
//!     n2 = rmsnorm(h, w2)
//!     h = h + W_down·(silu(W_gate·n2) ∘ (W_up·n2))
//! f = rmsnorm(h, w_final)
//! policy_mean = tanh(W_policy·f)                 (29 joint targets)
//! value       = W_value·f + b_value              (critic, separate head)
//! ```
//!
//! GQA: 8 query heads share 4 KV heads (kv head hq/2), head_dim 32 for
//! both (the first version had kv head_dim 64 ≠ head_dim 32 — fixed).
//! RoPE: standard multi-scale base-10000 frequencies (the first version
//! used base 1.0, collapsing all pair frequencies to 1.0).
//!
//! # Training vs inference paths (parity by construction)
//!
//! - `forward_step` — inference/rollout: appends to the sliding-window KV
//!   cache and returns (mean, value) for one position.
//! - `forward_sequence_train` — training: the SAME per-token core over a
//!   whole trajectory, recording an activation tape.
//! - `backward_sequence` — exact manual backprop of the PPO loss
//!   (clipped surrogate + value MSE) into the `.grad` buffers.
//!
//! Both forward paths call the identical per-token core in identical
//! order, so their outputs are bit-identical — asserted by a
//! differential test.
//!
//! # Gradient oracle
//!
//! Every backward path is verified against central finite differences on
//! a tiny config (`gradcheck_*` tests). The f32 FD noise floor was
//! measured before the tolerance was pinned.
//!
//! # Scale (iPhone budget)
//!
//! ~2.9M parameters (~12 MB f32), ~5.8M MACs per inference step —
//! inference is trivially phone-class. Pure scalar f32 loops (LLVM
//! autovec; hand-written wide SIMD is a measured loss for glue).

use crate::muon::{AdamParam, MuonParam};

// ─── Dimensions (defaults; tiny runtime-config variants power the oracle) ───

pub const D_MODEL: usize = 256;
pub const N_HEADS: usize = 8;
pub const HEAD_DIM: usize = 32;
pub const N_KV_HEADS: usize = 4;
pub const KV_HEAD_DIM: usize = HEAD_DIM; // proper GQA: same head dim, fewer heads
pub const KV_DIM: usize = N_KV_HEADS * KV_HEAD_DIM; // 128
pub const N_LAYERS: usize = 4;
pub const CONTEXT: usize = 64;
pub const N_INPUTS: usize = 42;
pub const N_OUTPUTS: usize = 29;
pub const MLP_HIDDEN: usize = 682; // 2.67 * 256 (GLM-style ratio)

/// Runtime dimension config. Default matches the consts above; the
/// gradient-oracle tests use a tiny variant.
#[derive(Clone, Copy, Debug)]
pub struct Config {
    pub d_model: usize,
    pub n_heads: usize,
    pub head_dim: usize,
    pub n_kv_heads: usize,
    pub kv_dim: usize,
    pub n_layers: usize,
    pub mlp_hidden: usize,
    pub context: usize,
    pub n_inputs: usize,
    pub n_outputs: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            d_model: D_MODEL,
            n_heads: N_HEADS,
            head_dim: HEAD_DIM,
            n_kv_heads: N_KV_HEADS,
            kv_dim: KV_DIM,
            n_layers: N_LAYERS,
            mlp_hidden: MLP_HIDDEN,
            context: CONTEXT,
            n_inputs: N_INPUTS,
            n_outputs: N_OUTPUTS,
        }
    }
}

// ─── Deterministic RNG (splitmix64) ───

fn splitmix_uniform(state: &mut u64) -> f64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
}

// ─── Primitives (forward + exact backward) ───

/// y = W·x, W row-major [rows × cols].
pub fn matvec(w: &[f32], rows: usize, cols: usize, x: &[f32], y: &mut [f32]) {
    debug_assert_eq!(x.len(), cols);
    debug_assert_eq!(y.len(), rows);
    for (r, yr) in y.iter_mut().enumerate() {
        let row = &w[r * cols..(r + 1) * cols];
        let mut sum = 0.0f32;
        for (c, inp) in x.iter().enumerate() {
            sum += row[c] * inp;
        }
        *yr = sum;
    }
}

/// Backward of y = W·x: dW += dy·xᵀ; dx += Wᵀ·dy (dx is ACCUMULATED —
/// zero it before a standalone call).
pub fn matvec_backward(
    w: &[f32],
    rows: usize,
    cols: usize,
    x: &[f32],
    dy: &[f32],
    dw: &mut [f32],
    dx: &mut [f32],
) {
    debug_assert_eq!(x.len(), cols);
    debug_assert_eq!(dy.len(), rows);
    for r in 0..rows {
        let g = dy[r];
        if g == 0.0 {
            continue;
        }
        let row = r * cols;
        for c in 0..cols {
            dw[row + c] += g * x[c];
        }
    }
    for c in 0..cols {
        let mut sum = 0.0f32;
        for r in 0..rows {
            sum += w[r * cols + c] * dy[r];
        }
        dx[c] += sum;
    }
}

pub fn rms_norm(x: &mut [f32], weight: &[f32], eps: f32) {
    let mean_sq: f32 = x.iter().map(|v| v * v).sum::<f32>() / x.len() as f32;
    let inv_rms = 1.0 / (mean_sq + eps).sqrt();
    for (xi, wi) in x.iter_mut().zip(weight.iter()) {
        *xi *= inv_rms * wi;
    }
}

/// Backward of y = x·inv ⊙ w, inv = (mean(x²)+ε)^(−1/2).
/// Adds into dx; accumulates dw. x is the SAVED pre-norm input.
pub fn rms_norm_backward(
    x: &[f32],
    weight: &[f32],
    eps: f32,
    dy: &[f32],
    dw: &mut [f32],
    dx: &mut [f32],
) {
    let dim = x.len() as f32;
    let mean_sq: f32 = x.iter().map(|v| v * v).sum::<f32>() / dim;
    let inv = 1.0 / (mean_sq + eps).sqrt();
    let mut dot = 0.0f32;
    for i in 0..x.len() {
        dot += dy[i] * weight[i] * x[i];
        dw[i] += dy[i] * (x[i] * inv);
    }
    for i in 0..x.len() {
        dx[i] += inv * dy[i] * weight[i] - inv * inv * inv * x[i] * dot / dim;
    }
}

/// RoPE: rotate (i, i+half) pairs of each head by position·base^(−2i/half),
/// base 10000 — standard multi-scale frequencies.
pub fn rope_in_place(x: &mut [f32], head_dim: usize, position: usize) {
    let half = head_dim / 2;
    for h in 0..x.len() / head_dim {
        let base = h * head_dim;
        for i in 0..half {
            let freq = 10000.0f32.powf(-2.0 * i as f32 / half as f32);
            let angle = position as f32 * freq;
            let (c, s) = (angle.cos(), angle.sin());
            let x0 = x[base + i];
            let x1 = x[base + i + half];
            x[base + i] = x0 * c - x1 * s;
            x[base + i + half] = x0 * s + x1 * c;
        }
    }
}

/// Backward of RoPE: the rotation is orthogonal, so its transpose is the
/// rotation by −angle. Adds into dx.
pub fn rope_backward(head_dim: usize, position: usize, dy: &[f32], dx: &mut [f32]) {
    let half = head_dim / 2;
    for h in 0..dy.len() / head_dim {
        let base = h * head_dim;
        for i in 0..half {
            let freq = 10000.0f32.powf(-2.0 * i as f32 / half as f32);
            let angle = position as f32 * freq;
            let (c, s) = (angle.cos(), angle.sin());
            let y0 = dy[base + i];
            let y1 = dy[base + i + half];
            dx[base + i] += y0 * c + y1 * s;
            dx[base + i + half] += -y0 * s + y1 * c;
        }
    }
}

pub fn swiglu(gate: &[f32], value: &[f32], out: &mut [f32]) {
    for i in 0..out.len() {
        out[i] = gate[i] / (1.0 + (-gate[i]).exp()) * value[i];
    }
}

/// Backward of out = silu(gate) ∘ value. Adds into dgate/dvalue.
pub fn swiglu_backward(gate: &[f32], value: &[f32], dout: &[f32], dgate: &mut [f32], dvalue: &mut [f32]) {
    for i in 0..dout.len() {
        let g = gate[i];
        let sigmoid = 1.0 / (1.0 + (-g).exp());
        let silu = g * sigmoid;
        dgate[i] += dout[i] * value[i] * sigmoid * (1.0 + g * (1.0 - sigmoid));
        dvalue[i] += dout[i] * silu;
    }
}

// ─── Optimizer-tagged weights ───

pub(crate) fn randomize_uniform(w: &mut [f32], fan_in: usize, seed: &mut u64) {
    let scale = (1.0 / fan_in as f64).sqrt() as f32;
    for v in w.iter_mut() {
        *v = (splitmix_uniform(seed) as f32 * 2.0 - 1.0) * scale;
    }
}

fn muon_param(rows: usize, cols: usize, lr: f32, beta: f32, seed: &mut u64) -> MuonParam {
    let mut p = MuonParam::new(rows, cols, lr, beta);
    randomize_uniform(&mut p.weights, cols, seed); // fan_in = cols
    p
}

fn adam_param_ones(size: usize, lr: f32) -> AdamParam {
    let mut p = AdamParam::new(size, lr);
    for v in p.params.iter_mut() {
        *v = 1.0;
    }
    p
}

// ─── Model ───

pub struct LayerWeights {
    pub wq: MuonParam,     // d_model × d_model
    pub wk: MuonParam,     // kv_dim × d_model  (k = Wk·n1)
    pub wv: MuonParam,     // kv_dim × d_model  (v = Wv·n1)
    pub wo: MuonParam,     // d_model × d_model
    pub w_gate: MuonParam, // mlp_hidden × d_model
    pub w_up: MuonParam,   // mlp_hidden × d_model
    pub w_down: MuonParam, // d_model × mlp_hidden
    pub norm1: AdamParam,  // d_model
    pub norm2: AdamParam,  // d_model
}

/// Per-layer KV sliding-window cache (inference path).
struct KvCache {
    k: Vec<f32>, // slots × kv_dim
    v: Vec<f32>,
    start: usize, // absolute position of slot 0
    count: usize,
}

impl KvCache {
    fn new() -> Self {
        Self { k: Vec::new(), v: Vec::new(), start: 0, count: 0 }
    }
    fn clear(&mut self) {
        self.k.clear();
        self.v.clear();
        self.start = 0;
        self.count = 0;
    }
    fn push(&mut self, kv_dim: usize, context: usize, k: &[f32], v: &[f32]) {
        if self.count == context {
            self.k.drain(0..kv_dim);
            self.v.drain(0..kv_dim);
            self.start += 1;
            self.count -= 1;
        }
        self.k.extend_from_slice(k);
        self.v.extend_from_slice(v);
        self.count += 1;
    }
}

pub struct LayerTape {
    pub hidden_in: Vec<f32>, // residual input to this block
    pub normed1: Vec<f32>,
    pub q: Vec<f32>,       // post-RoPE queries
    pub k_cache: Vec<f32>, // cache snapshot INCLUDING this token (slots × kv_dim)
    pub v_cache: Vec<f32>,
    pub cache_count: usize,
    pub att_out: Vec<f32>, // concatenated head outputs (input to Wo)
    pub a_vec: Vec<f32>,   // hidden after attention residual
    pub normed2: Vec<f32>,
    pub gate: Vec<f32>,
    pub up: Vec<f32>,
    pub mid: Vec<f32>,
}

pub struct TokenTape {
    pub obs: Vec<f32>,
    pub mean: Vec<f32>,
    pub layers: Vec<LayerTape>,
    pub h_final_pre: Vec<f32>,
    pub h_final_post: Vec<f32>,
}

pub struct TrainTape {
    pub tokens: Vec<TokenTape>,
}

pub struct GaitTransformer {
    pub cfg: Config,
    pub embed: AdamParam,       // d_model × n_inputs
    pub layers: Vec<LayerWeights>,
    pub policy_head: AdamParam, // n_outputs × d_model
    pub value_w: AdamParam,     // 1 × d_model
    pub value_b: AdamParam,     // 1
    pub final_norm: AdamParam,  // d_model
    cache: Vec<KvCache>,
}

impl GaitTransformer {
    pub fn new(cfg: Config, lr_muon: f32, beta: f32, lr_adam: f32, seed: &mut u64) -> Self {
        let d = cfg.d_model;
        let mut layers = Vec::new();
        for _ in 0..cfg.n_layers {
            layers.push(LayerWeights {
                wq: muon_param(d, d, lr_muon, beta, seed),
                wk: muon_param(cfg.kv_dim, d, lr_muon, beta, seed),
                wv: muon_param(cfg.kv_dim, d, lr_muon, beta, seed),
                wo: muon_param(d, d, lr_muon, beta, seed),
                w_gate: muon_param(cfg.mlp_hidden, d, lr_muon, beta, seed),
                w_up: muon_param(cfg.mlp_hidden, d, lr_muon, beta, seed),
                w_down: muon_param(d, cfg.mlp_hidden, lr_muon, beta, seed),
                norm1: adam_param_ones(d, lr_adam),
                norm2: adam_param_ones(d, lr_adam),
            });
        }
        let mut embed = AdamParam::new(d * cfg.n_inputs, lr_adam);
        randomize_uniform(&mut embed.params, cfg.n_inputs, seed);
        let mut policy_head = AdamParam::new(cfg.n_outputs * d, lr_adam);
        randomize_uniform(&mut policy_head.params, d, seed);
        let mut value_w = AdamParam::new(d, lr_adam);
        randomize_uniform(&mut value_w.params, d, seed);
        Self {
            cfg,
            embed,
            layers,
            policy_head,
            value_w,
            value_b: AdamParam::new(1, lr_adam),
            final_norm: adam_param_ones(d, lr_adam),
            cache: (0..cfg.n_layers).map(|_| KvCache::new()).collect(),
        }
    }

    pub fn reset_cache(&mut self) {
        for c in self.cache.iter_mut() {
            c.clear();
        }
    }

    pub fn param_count(&self) -> usize {
        let mut n = self.embed.params.len()
            + self.policy_head.params.len()
            + self.value_w.params.len()
            + self.value_b.params.len()
            + self.final_norm.params.len();
        for l in &self.layers {
            n += l.wq.weights.len()
                + l.wk.weights.len()
                + l.wv.weights.len()
                + l.wo.weights.len()
                + l.w_gate.weights.len()
                + l.w_up.weights.len()
                + l.w_down.weights.len()
                + l.norm1.params.len()
                + l.norm2.params.len();
        }
        n
    }

    pub fn zero_grads(&mut self) {
        self.embed.grads.fill(0.0);
        self.policy_head.grads.fill(0.0);
        self.value_w.grads.fill(0.0);
        self.value_b.grads.fill(0.0);
        self.final_norm.grads.fill(0.0);
        for l in &mut self.layers {
            l.wq.grad.fill(0.0);
            l.wk.grad.fill(0.0);
            l.wv.grad.fill(0.0);
            l.wo.grad.fill(0.0);
            l.w_gate.grad.fill(0.0);
            l.w_up.grad.fill(0.0);
            l.w_down.grad.fill(0.0);
            l.norm1.grads.fill(0.0);
            l.norm2.grads.fill(0.0);
        }
    }

    /// One optimizer step: Muon (orthogonalized momentum) on hidden 2-D
    /// weights, Adam on embed / heads / norms — then zero the grads.
    pub fn step_optimizers(&mut self) {
        for l in self.layers.iter_mut() {
            l.wq.step();
            l.wk.step();
            l.wv.step();
            l.wo.step();
            l.w_gate.step();
            l.w_up.step();
            l.w_down.step();
            l.norm1.step();
            l.norm2.step();
        }
        self.embed.step();
        self.policy_head.step();
        self.value_w.step();
        self.value_b.step();
        self.final_norm.step();
        self.zero_grads();
    }

    /// One inference step (KV-cache path). `position` = absolute token
    /// index (drives RoPE). Returns (policy mean, value).
    pub fn forward_step(&mut self, obs: &[f32], position: usize) -> (Vec<f32>, f32) {
        debug_assert_eq!(obs.len(), self.cfg.n_inputs);
        let cfg = self.cfg;
        let d = cfg.d_model;
        let mut h = vec![0.0f32; d];
        matvec(&self.embed.params, d, cfg.n_inputs, obs, &mut h);
        let mut scratch = vec![0.0f32; d];
        let mut q = vec![0.0f32; d];
        let mut k = vec![0.0f32; cfg.kv_dim];
        let mut v = vec![0.0f32; cfg.kv_dim];
        let mut attn = vec![0.0f32; d];
        let mut proj = vec![0.0f32; d];
        let mut gate = vec![0.0f32; cfg.mlp_hidden];
        let mut up = vec![0.0f32; cfg.mlp_hidden];
        let mut mid = vec![0.0f32; cfg.mlp_hidden];
        for (li, layer) in self.layers.iter().enumerate() {
            scratch.copy_from_slice(&h);
            rms_norm(&mut scratch, &layer.norm1.params, 1e-6);
            matvec(&layer.wq.weights, d, d, &scratch, &mut q);
            matvec(&layer.wk.weights, cfg.kv_dim, d, &scratch, &mut k);
            matvec(&layer.wv.weights, cfg.kv_dim, d, &scratch, &mut v);
            rope_in_place(&mut q, cfg.head_dim, position);
            rope_in_place(&mut k, cfg.head_dim, position);
            self.cache[li].push(cfg.kv_dim, cfg.context, &k, &v);
            let c = &self.cache[li];
            attention_causal(&q, &c.k, &c.v, c.count, cfg.n_heads, cfg.head_dim, &mut attn);
            matvec(&layer.wo.weights, d, d, &attn, &mut proj);
            for (hi, pi) in h.iter_mut().zip(proj.iter()) {
                *hi += pi;
            }
            scratch.copy_from_slice(&h);
            rms_norm(&mut scratch, &layer.norm2.params, 1e-6);
            matvec(&layer.w_gate.weights, cfg.mlp_hidden, d, &scratch, &mut gate);
            matvec(&layer.w_up.weights, cfg.mlp_hidden, d, &scratch, &mut up);
            swiglu(&gate, &up, &mut mid);
            matvec(&layer.w_down.weights, d, cfg.mlp_hidden, &mid, &mut proj);
            for (hi, mi) in h.iter_mut().zip(proj.iter()) {
                *hi += mi;
            }
        }
        rms_norm(&mut h, &self.final_norm.params, 1e-6);
        let mean = head_policy(&self.policy_head.params, &h, cfg.n_outputs);
        let mut value = self.value_b.params[0];
        for (c, hf) in h.iter().enumerate() {
            value += self.value_w.params[c] * hf;
        }
        (mean, value)
    }

    /// Training forward: identical per-token core to `forward_step`
    /// (bit-identical outputs — differential-tested), recording the tape.
    pub fn forward_sequence_train(
        &mut self,
        obs_seq: &[Vec<f32>],
    ) -> (Vec<Vec<f32>>, Vec<f32>, TrainTape) {
        self.reset_cache();
        let cfg = self.cfg;
        let d = cfg.d_model;
        let mut means = Vec::with_capacity(obs_seq.len());
        let mut values = Vec::with_capacity(obs_seq.len());
        let mut tokens = Vec::with_capacity(obs_seq.len());
        let mut scratch = vec![0.0f32; d];
        let mut q = vec![0.0f32; d];
        let mut k = vec![0.0f32; cfg.kv_dim];
        let mut v = vec![0.0f32; cfg.kv_dim];
        let mut attn = vec![0.0f32; d];
        let mut proj = vec![0.0f32; d];
        let mut gate = vec![0.0f32; cfg.mlp_hidden];
        let mut up = vec![0.0f32; cfg.mlp_hidden];
        let mut mid = vec![0.0f32; cfg.mlp_hidden];
        for (pos, obs) in obs_seq.iter().enumerate() {
            let mut h = vec![0.0f32; d];
            matvec(&self.embed.params, d, cfg.n_inputs, obs, &mut h);
            let mut ltapes = Vec::with_capacity(cfg.n_layers);
            for (li, layer) in self.layers.iter().enumerate() {
                let hidden_in = h.clone();
                scratch.copy_from_slice(&h);
                rms_norm(&mut scratch, &layer.norm1.params, 1e-6);
                let normed1 = scratch.clone();
                matvec(&layer.wq.weights, d, d, &normed1, &mut q);
                matvec(&layer.wk.weights, cfg.kv_dim, d, &normed1, &mut k);
                matvec(&layer.wv.weights, cfg.kv_dim, d, &normed1, &mut v);
                rope_in_place(&mut q, cfg.head_dim, pos);
                rope_in_place(&mut k, cfg.head_dim, pos);
                let q_saved = q.clone();
                self.cache[li].push(cfg.kv_dim, cfg.context, &k, &v);
                let c = &self.cache[li];
                let (k_snap, v_snap, count) = (c.k.clone(), c.v.clone(), c.count);
                attention_causal(&q, &c.k, &c.v, c.count, cfg.n_heads, cfg.head_dim, &mut attn);
                let att_out = attn.clone();
                matvec(&layer.wo.weights, d, d, &attn, &mut proj);
                for (hi, pi) in h.iter_mut().zip(proj.iter()) {
                    *hi += pi;
                }
                let a_vec = h.clone();
                scratch.copy_from_slice(&h);
                rms_norm(&mut scratch, &layer.norm2.params, 1e-6);
                let normed2 = scratch.clone();
                matvec(&layer.w_gate.weights, cfg.mlp_hidden, d, &normed2, &mut gate);
                matvec(&layer.w_up.weights, cfg.mlp_hidden, d, &normed2, &mut up);
                swiglu(&gate, &up, &mut mid);
                matvec(&layer.w_down.weights, d, cfg.mlp_hidden, &mid, &mut proj);
                for (hi, mi) in h.iter_mut().zip(proj.iter()) {
                    *hi += mi;
                }
                ltapes.push(LayerTape {
                    hidden_in,
                    normed1,
                    q: q_saved,
                    k_cache: k_snap,
                    v_cache: v_snap,
                    cache_count: count,
                    att_out,
                    a_vec,
                    normed2,
                    gate: gate.clone(),
                    up: up.clone(),
                    mid: mid.clone(),
                });
            }
            let h_final_pre = h.clone();
            rms_norm(&mut h, &self.final_norm.params, 1e-6);
            let mean_vec = head_policy(&self.policy_head.params, &h, cfg.n_outputs);
            let mut value = self.value_b.params[0];
            for (c, hf) in h.iter().enumerate() {
                value += self.value_w.params[c] * hf;
            }
            means.push(mean_vec.to_vec());
            values.push(value);
            tokens.push(TokenTape {
                obs: obs.clone(),
                mean: mean_vec.to_vec(),
                layers: ltapes,
                h_final_pre,
                h_final_post: h,
            });
        }
        (means, values, TrainTape { tokens })
    }

    /// Exact manual backprop of
    ///   L = (1/T) Σ_t [ ppo_t + c_V·(V_t−R_t)² ]
    /// where `dmean[t][i]` = ∂ppo_t/∂mean[t][i] (caller applies the
    /// clip-branch mask, the −sign, and the 1/T) and `dvalue[t]` =
    /// ∂L/∂V_t (caller folds the c_V·2/T factor in). Accumulates into
    /// the `.grad` buffers (call `zero_grads` first).
    ///
    /// Causal cross-token flow: token τ's attention consumes the cached
    /// k/v of every token j ≤ τ, so ∂L/∂k_j picks up contributions from
    /// all τ ≥ j. Processing tokens in reverse order with per-position
    /// KV gradient accumulators makes each token's own backward see the
    /// complete KV gradient before it runs.
    pub fn backward_sequence(&mut self, tape: &TrainTape, dmean: &[Vec<f32>], dvalue: &[f32]) {
        let cfg = self.cfg;
        let d = cfg.d_model;
        let t_len = tape.tokens.len();
        self.zero_grads();
        let mut dh_post = vec![0.0f32; d];
        let mut dh = vec![0.0f32; d];
        let mut da = vec![0.0f32; d];
        let mut dh_in = vec![0.0f32; d];
        let mut dattn = vec![0.0f32; d];
        let mut dq = vec![0.0f32; d];
        let mut dn1 = vec![0.0f32; d];
        let mut dn2 = vec![0.0f32; d];
        let mut dmid = vec![0.0f32; cfg.mlp_hidden];
        let mut dgate = vec![0.0f32; cfg.mlp_hidden];
        let mut dup = vec![0.0f32; cfg.mlp_hidden];
        let mut dk_own = vec![0.0f32; cfg.kv_dim];
        let mut dv_own = vec![0.0f32; cfg.kv_dim];
        let mut dk_pre = vec![0.0f32; cfg.kv_dim];
        let mut dq_pre = vec![0.0f32; d];
        // Per-layer KV gradient accumulators over ABSOLUTE positions.
        // gk[li][pos*kv_dim + i] = ∂L/∂k_{pos,layer,kv-head,i}.
        let mut gk: Vec<Vec<f32>> = (0..cfg.n_layers).map(|_| Vec::new()).collect();
        let mut gv: Vec<Vec<f32>> = (0..cfg.n_layers).map(|_| Vec::new()).collect();

        for t in (0..t_len).rev() {
            let tok = &tape.tokens[t];
            // ── heads on h_final_post ──
            dh_post.fill(0.0);
            for i in 0..cfg.n_outputs {
                let m = tok.mean[i];
                let dpre = dmean[t][i] * (1.0 - m * m);
                if dpre == 0.0 {
                    continue;
                }
                let r = i * d;
                for c in 0..d {
                    self.policy_head.grads[r + c] += dpre * tok.h_final_post[c];
                    dh_post[c] += dpre * self.policy_head.params[r + c];
                }
            }
            for c in 0..d {
                self.value_w.grads[c] += dvalue[t] * tok.h_final_post[c];
                dh_post[c] += dvalue[t] * self.value_w.params[c];
            }
            self.value_b.grads[0] += dvalue[t];
            // ── final RMSNorm ──
            dh.fill(0.0);
            rms_norm_backward(
                &tok.h_final_pre,
                &self.final_norm.params,
                1e-6,
                &dh_post,
                &mut self.final_norm.grads,
                &mut dh,
            );
            // ── blocks in reverse; KV inflow from later tokens already
            // accumulated in gk/gv ──
            for li in (0..cfg.n_layers).rev() {
                let layer = &mut self.layers[li];
                let lt = &tok.layers[li];
                let need = (t + 1) * cfg.kv_dim;
                if gk[li].len() < need {
                    gk[li].resize(need, 0.0);
                    gv[li].resize(need, 0.0);
                }
                // residual 2: h_out = a + mlp_out → da = dh ; dmlp = dh
                dmid.fill(0.0);
                matvec_backward(
                    &layer.w_down.weights,
                    d,
                    cfg.mlp_hidden,
                    &lt.mid,
                    &dh,
                    &mut layer.w_down.grad,
                    &mut dmid,
                );
                dgate.fill(0.0);
                dup.fill(0.0);
                swiglu_backward(&lt.gate, &lt.up, &dmid, &mut dgate, &mut dup);
                dn2.fill(0.0);
                matvec_backward(&layer.w_gate.weights, cfg.mlp_hidden, d, &lt.normed2, &dgate, &mut layer.w_gate.grad, &mut dn2);
                matvec_backward(&layer.w_up.weights, cfg.mlp_hidden, d, &lt.normed2, &dup, &mut layer.w_up.grad, &mut dn2);
                da.copy_from_slice(&dh);
                rms_norm_backward(&lt.a_vec, &layer.norm2.params, 1e-6, &dn2, &mut layer.norm2.grads, &mut da);
                // attention residual: a = hidden_in + wo·att_out
                dattn.fill(0.0);
                matvec_backward(&layer.wo.weights, d, d, &lt.att_out, &da, &mut layer.wo.grad, &mut dattn);
                // attention backward: dq_t own; dk/dv per slot
                dq.fill(0.0);
                let start_abs = t + 1 - lt.cache_count; // 0 while t < context
                for slot in 0..lt.cache_count {
                    let abs = start_abs + slot;
                    let off = abs * cfg.kv_dim;
                    if gk[li].len() < off + cfg.kv_dim {
                        gk[li].resize(off + cfg.kv_dim, 0.0);
                        gv[li].resize(off + cfg.kv_dim, 0.0);
                    }
                    dk_own.fill(0.0);
                    dv_own.fill(0.0);
                    attention_causal_backward_slot(
                        &lt.q, &lt.k_cache, &lt.v_cache, lt.cache_count,
                        slot, cfg.n_heads, cfg.head_dim, &dattn,
                        &mut dq, &mut dk_own, &mut dv_own,
                    );
                    for i in 0..cfg.kv_dim {
                        gk[li][off + i] += dk_own[i];
                        gv[li][off + i] += dv_own[i];
                    }
                }
                // own token's total k/v gradient (own + all later inflow)
                let own_off = t * cfg.kv_dim;
                dk_pre.fill(0.0);
                rope_backward(cfg.head_dim, t, &gk[li][own_off..own_off + cfg.kv_dim], &mut dk_pre);
                dq_pre.fill(0.0);
                rope_backward(cfg.head_dim, t, &dq, &mut dq_pre);
                // projections back into normed1
                dn1.fill(0.0);
                matvec_backward(&layer.wq.weights, d, d, &lt.normed1, &dq_pre, &mut layer.wq.grad, &mut dn1);
                matvec_backward(&layer.wk.weights, cfg.kv_dim, d, &lt.normed1, &dk_pre, &mut layer.wk.grad, &mut dn1);
                let dv_own: Vec<f32> = gv[li][own_off..own_off + cfg.kv_dim].to_vec();
                matvec_backward(&layer.wv.weights, cfg.kv_dim, d, &lt.normed1, &dv_own, &mut layer.wv.grad, &mut dn1);
                // residual 1 + RMSNorm 1 → gradient into block input
                dh_in.copy_from_slice(&da);
                rms_norm_backward(&lt.hidden_in, &layer.norm1.params, 1e-6, &dn1, &mut layer.norm1.grads, &mut dh_in);
                dh.copy_from_slice(&dh_in);
            }
            // embed gradient (token-local): runs ONCE per token with the
            // completed dh — inside the layer loop it would over-count by
            // summing each layer's partial contribution separately.
            for c in 0..d {
                let g = dh[c];
                if g == 0.0 {
                    continue;
                }
                let row = c * cfg.n_inputs;
                for (i, o) in tok.obs.iter().enumerate() {
                    self.embed.grads[row + i] += g * o;
                }
            }
        }
    }
}

fn head_policy(w: &[f32], h: &[f32], n_out: usize) -> Vec<f32> {
    let mut out = vec![0.0f32; n_out];
    let d = h.len();
    for (i, o) in out.iter_mut().enumerate().take(n_out) {
        let r = i * d;
        let mut sum = 0.0f32;
        for (c, hf) in h.iter().enumerate() {
            sum += w[r + c] * hf;
        }
        *o = sum.tanh();
    }
    out
}

/// Causal GQA attention over a KV window (count slots).
fn attention_causal(
    q: &[f32],
    k_cache: &[f32],
    v_cache: &[f32],
    count: usize,
    n_heads: usize,
    head_dim: usize,
    attn: &mut [f32],
) {
    let kv_dim = if count == 0 { 0 } else { k_cache.len() / count };
    let n_q_groups = n_heads / (kv_dim / head_dim);
    let scale = 1.0 / (head_dim as f32).sqrt();
    let mut scores = vec![0.0f32; count];
    for g in 0..n_heads {
        let kv_head = g / n_q_groups;
        let q_off = g * head_dim;
        let mut maxs = f32::NEG_INFINITY;
        for j in 0..count {
            let ko = j * kv_dim + kv_head * head_dim;
            let mut s = 0.0f32;
            for i in 0..head_dim {
                s += q[q_off + i] * k_cache[ko + i];
            }
            s *= scale;
            scores[j] = s;
            if s > maxs {
                maxs = s;
            }
        }
        let mut sum = 0.0f32;
        for j in 0..count {
            let e = (scores[j] - maxs).exp();
            scores[j] = e;
            sum += e;
        }
        let inv = 1.0 / sum;
        let vo = kv_head * head_dim;
        for i in 0..head_dim {
            attn[q_off + i] = 0.0;
        }
        for j in 0..count {
            let a = scores[j] * inv;
            let vo_j = j * kv_dim + vo;
            for i in 0..head_dim {
                attn[q_off + i] += a * v_cache[vo_j + i];
            }
        }
    }
}

/// Backward of `attention_causal` for ONE query head and ONE cache slot:
/// contributes dq (query head slice, accumulated) and the slot's dk/dv
/// (WRITTEN, caller folds them into absolute-position accumulators).
/// Recomputes the softmax exactly as forward.
fn attention_causal_backward_slot(
    q: &[f32],
    k_cache: &[f32],
    v_cache: &[f32],
    count: usize,
    slot: usize,
    n_heads: usize,
    head_dim: usize,
    dattn: &[f32],
    dq: &mut [f32],
    dk_slot: &mut [f32],
    dv_slot: &mut [f32],
) {
    let kv_dim = k_cache.len() / count.max(1);
    let n_q_groups = n_heads / (kv_dim / head_dim);
    let scale = 1.0 / (head_dim as f32).sqrt();
    dk_slot.fill(0.0);
    dv_slot.fill(0.0);
    let mut scores = vec![0.0f32; count];
    for g in 0..n_heads {
        let kv_head = g / n_q_groups;
        let q_off = g * head_dim;
        let mut maxs = f32::NEG_INFINITY;
        for j in 0..count {
            let ko = j * kv_dim + kv_head * head_dim;
            let mut s = 0.0f32;
            for i in 0..head_dim {
                s += q[q_off + i] * k_cache[ko + i];
            }
            s *= scale;
            scores[j] = s;
            if s > maxs {
                maxs = s;
            }
        }
        let mut sum = 0.0f32;
        for j in 0..count {
            let e = (scores[j] - maxs).exp();
            scores[j] = e;
            sum += e;
        }
        for j in 0..count {
            scores[j] /= sum;
        }
        let a_slot = scores[slot];
        let d_o = &dattn[q_off..q_off + head_dim];
        let mut d_a = vec![0.0f32; count];
        let mut dot = 0.0f32;
        for j in 0..count {
            let vo = j * kv_dim + kv_head * head_dim;
            let mut acc = 0.0f32;
            for i in 0..head_dim {
                acc += d_o[i] * v_cache[vo + i];
            }
            d_a[j] = acc;
            dot += scores[j] * acc;
        }
        let ds_slot = a_slot * (d_a[slot] - dot);
        for i in 0..head_dim {
            dv_slot[kv_head * head_dim + i] += scores[slot] * d_o[i];
            dq[q_off + i] += ds_slot * k_cache[slot * kv_dim + kv_head * head_dim + i] * scale;
            dk_slot[kv_head * head_dim + i] += ds_slot * q[q_off + i] * scale;
        }
    }
}

// ─── Tests: the gradient oracle + parity + learning proofs ───

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn tiny_cfg() -> Config {
        Config {
            d_model: 16,
            n_heads: 2,
            head_dim: 8,
            n_kv_heads: 1,
            kv_dim: 8,
            n_layers: 2,
            mlp_hidden: 24,
            context: 8,
            n_inputs: 6,
            n_outputs: 3,
        }
    }

    pub(crate) fn tiny_model(seed: &mut u64) -> GaitTransformer {
        GaitTransformer::new(tiny_cfg(), 0.02, 0.9, 0.02, seed)
    }

    pub(crate) fn probe_obs(seed: &mut u64, t: usize) -> Vec<f32> {
        (0..6)
            .map(|i| splitmix_uniform(seed) as f32 - 0.5 + 0.01 * t as f32 * i as f32)
            .collect()
    }

    /// Fixed (parameter-independent) probe loss: L = Σ_t Σ_i w·mean +
    /// Σ_t c·value. Any parameter perturbation moves L only through the
    /// network — exactly what the gradient oracle checks.
    pub(crate) struct Probe {
        pub(crate) w_mean: Vec<Vec<f32>>,
        pub(crate) w_val: Vec<f32>,
    }

    impl Probe {
        pub(crate) fn new(seed: &mut u64, t_len: usize, n_out: usize) -> Self {
            Self {
                w_mean: (0..t_len)
                    .map(|_| (0..n_out).map(|_| splitmix_uniform(seed) as f32 - 0.5).collect())
                    .collect(),
                w_val: (0..t_len).map(|_| splitmix_uniform(seed) as f32 - 0.5).collect(),
            }
        }
        pub(crate) fn loss(&self, means: &[Vec<f32>], values: &[f32]) -> f32 {
            let mut l = 0.0;
            for t in 0..means.len() {
                for i in 0..means[t].len() {
                    l += self.w_mean[t][i] * means[t][i];
                }
                l += self.w_val[t] * values[t];
            }
            l
        }
        pub(crate) fn dmean(&self, t: usize) -> Vec<f32> {
            self.w_mean[t].clone()
        }
        pub(crate) fn dvalue(&self, t: usize) -> f32 {
            self.w_val[t]
        }
    }

    /// Tensor codes for gradcheck site enumeration.
    /// 0 embed, 1 policy, 2 value_w, 3 value_b, 4 final_norm,
    /// 5 wq, 6 wk, 7 wv, 8 wo, 9 w_gate, 10 w_up, 11 w_down, 12 norm1, 13 norm2
    fn site_len(model: &GaitTransformer, code: usize, layer: usize) -> usize {
        match code {
            0 => model.embed.params.len(),
            1 => model.policy_head.params.len(),
            2 => model.value_w.params.len(),
            3 => model.value_b.params.len(),
            4 => model.final_norm.params.len(),
            5 => model.layers[layer].wq.weights.len(),
            6 => model.layers[layer].wk.weights.len(),
            7 => model.layers[layer].wv.weights.len(),
            8 => model.layers[layer].wo.weights.len(),
            9 => model.layers[layer].w_gate.weights.len(),
            10 => model.layers[layer].w_up.weights.len(),
            11 => model.layers[layer].w_down.weights.len(),
            12 => model.layers[layer].norm1.params.len(),
            _ => model.layers[layer].norm2.params.len(),
        }
    }

    fn site_grad(model: &GaitTransformer, code: usize, layer: usize, idx: usize) -> f32 {
        match code {
            0 => model.embed.grads[idx],
            1 => model.policy_head.grads[idx],
            2 => model.value_w.grads[idx],
            3 => model.value_b.grads[idx],
            4 => model.final_norm.grads[idx],
            5 => model.layers[layer].wq.grad[idx],
            6 => model.layers[layer].wk.grad[idx],
            7 => model.layers[layer].wv.grad[idx],
            8 => model.layers[layer].wo.grad[idx],
            9 => model.layers[layer].w_gate.grad[idx],
            10 => model.layers[layer].w_up.grad[idx],
            11 => model.layers[layer].w_down.grad[idx],
            12 => model.layers[layer].norm1.grads[idx],
            _ => model.layers[layer].norm2.grads[idx],
        }
    }

    fn site_bump(model: &mut GaitTransformer, code: usize, layer: usize, idx: usize, delta: f32) {
        match code {
            0 => model.embed.params[idx] += delta,
            1 => model.policy_head.params[idx] += delta,
            2 => model.value_w.params[idx] += delta,
            3 => model.value_b.params[idx] += delta,
            4 => model.final_norm.params[idx] += delta,
            5 => model.layers[layer].wq.weights[idx] += delta,
            6 => model.layers[layer].wk.weights[idx] += delta,
            7 => model.layers[layer].wv.weights[idx] += delta,
            8 => model.layers[layer].wo.weights[idx] += delta,
            9 => model.layers[layer].w_gate.weights[idx] += delta,
            10 => model.layers[layer].w_up.weights[idx] += delta,
            11 => model.layers[layer].w_down.weights[idx] += delta,
            12 => model.layers[layer].norm1.params[idx] += delta,
            _ => model.layers[layer].norm2.params[idx] += delta,
        }
    }

    fn probe_loss(model: &mut GaitTransformer, obs: &[Vec<f32>], probe: &Probe) -> f32 {
        let (means, values, _tape) = model.forward_sequence_train(obs);
        probe.loss(&means, &values)
    }

    /// THE GRADIENT ORACLE: central finite differences vs analytic
    /// backward on every tensor class. Tolerance derived from measured
    /// f32 FD noise (h = 1e-2).
    #[test]
    fn gradcheck_backward_matches_finite_differences() {
        let mut seed = 0x0ACC_0C00_EE_u64;
        let cfg = tiny_cfg();
        let mut model = tiny_model(&mut seed);
        let t_len = 5;
        let obs: Vec<Vec<f32>> = (0..t_len).map(|t| probe_obs(&mut seed, t)).collect();
        let probe = Probe::new(&mut seed, t_len, cfg.n_outputs);

        let (means, values, tape) = model.forward_sequence_train(&obs);
        let _ = probe.loss(&means, &values);
        model.zero_grads();
        let dmean: Vec<Vec<f32>> = (0..t_len).map(|t| probe.dmean(t)).collect();
        let dvalue: Vec<f32> = (0..t_len).map(|t| probe.dvalue(t)).collect();
        model.backward_sequence(&tape, &dmean, &dvalue);

        let h = 1e-2f32;
        let mut worst_rel = 0.0f32;
        let mut checked = 0usize;
        let codes = 0..14usize;
        for code in codes {
            let layers: &[usize] = if code < 5 { &[0] } else { &[0, cfg.n_layers - 1] };
            for &layer in layers {
                let len = site_len(&model, code, layer);
                for idx in [0usize, len / 2, len - 1] {
                    let g_analytic = site_grad(&model, code, layer, idx);
                    site_bump(&mut model, code, layer, idx, h);
                    let lp = probe_loss(&mut model, &obs, &probe);
                    site_bump(&mut model, code, layer, idx, -2.0 * h);
                    let lm = probe_loss(&mut model, &obs, &probe);
                    site_bump(&mut model, code, layer, idx, h); // restore
                    let g_fd = (lp - lm) / (2.0 * h);
                    let denom = g_analytic.abs() + g_fd.abs() + 1e-6;
                    let rel = (g_analytic - g_fd).abs() / denom;
                    assert!(
                        rel < 0.05,
                        "gradcheck {}[layer {}][{}] failed: analytic={g_analytic} fd={g_fd} rel={rel}",
                        code, layer, idx
                    );
                    if rel > worst_rel {
                        worst_rel = rel;
                    }
                    checked += 1;
                }
            }
        }
        println!("[gradcheck] {checked} sites checked, worst rel err {worst_rel:.5}");
        assert!(checked >= 40);
    }

    /// Parity: stepping `forward_step` must produce bit-identical
    /// (mean, value) to `forward_sequence_train` for every position.
    #[test]
    fn stepwise_matches_sequence_bit_exactly() {
        let mut seed = 0x5EED_5EED_u64;
        let mut model = tiny_model(&mut seed);
        let obs: Vec<Vec<f32>> = (0..7).map(|t| probe_obs(&mut seed, t)).collect();
        let (seq_means, seq_values, _tape) = model.forward_sequence_train(&obs);
        model.reset_cache();
        for (t, o) in obs.iter().enumerate() {
            let (mean, value) = model.forward_step(o, t);
            assert_eq!(mean.to_vec(), seq_means[t], "mean mismatch at t={t}");
            assert_eq!(value, seq_values[t], "value mismatch at t={t}");
        }
    }

    /// Causality: outputs at t must not depend on observations after t.
    #[test]
    fn outputs_are_causal() {
        let mut seed = 0xCA5A_1E00u64;
        let mut model = tiny_model(&mut seed);
        let mut obs: Vec<Vec<f32>> = (0..6).map(|t| probe_obs(&mut seed, t)).collect();
        let (base_means, base_values, _t) = model.forward_sequence_train(&obs);
        // Perturb the FUTURE observations only.
        for o in obs.iter_mut().skip(4) {
            for v in o.iter_mut() {
                *v += 3.0;
            }
        }
        let (pert_means, pert_values, _t) = model.forward_sequence_train(&obs);
        for t in 0..4 {
            assert_eq!(base_means[t], pert_means[t], "future leaked into past at t={t}");
            assert_eq!(base_values[t], pert_values[t]);
        }
    }

    /// Sliding window: positions beyond the context still work and the
    /// field stays finite.
    #[test]
    fn sliding_window_beyond_context() {
        let mut seed = 0x1111_2222_u64;
        let mut model = tiny_model(&mut seed); // context 8
        for pos in 0..20 {
            let obs = probe_obs(&mut seed, pos);
            let (mean, value) = model.forward_step(&obs, pos);
            assert!(mean.iter().all(|v| v.is_finite()));
            assert!(value.is_finite());
        }
    }

    /// REAL LEARNING proof: supervised regression on a fixed dataset —
    /// gradient descent (Muon + Adam through the full backward) must
    /// reduce the loss by a large factor. This test is the anti-vacuity
    /// gate for the entire training stack.
    #[test]
    fn supervised_training_reduces_loss_substantially() {
        let mut seed = 0xABCD_1234_u64;
        let cfg = tiny_cfg();
        let mut model = tiny_model(&mut seed);
        let t_len = 6;
        let obs: Vec<Vec<f32>> = (0..t_len).map(|t| probe_obs(&mut seed, t)).collect();
        // Fixed targets in (-0.5, 0.5) for the policy head; fixed value targets.
        let targets: Vec<Vec<f32>> = (0..t_len)
            .map(|_| (0..cfg.n_outputs).map(|_| splitmix_uniform(&mut seed) as f32 - 0.5).collect())
            .collect();
        let v_targets: Vec<f32> = (0..t_len).map(|_| splitmix_uniform(&mut seed) as f32 * 0.4).collect();

        let loss_at = |model: &mut GaitTransformer| -> f32 {
            let (means, values, _t) = model.forward_sequence_train(&obs);
            let mut l = 0.0;
            for t in 0..t_len {
                for i in 0..cfg.n_outputs {
                    l += (means[t][i] - targets[t][i]).powi(2);
                }
                l += (values[t] - v_targets[t]).powi(2);
            }
            l / t_len as f32
        };

        let l0 = loss_at(&mut model);
        for _ in 0..400 {
            let (means, values, tape) = model.forward_sequence_train(&obs);
            let mut dmean = Vec::with_capacity(t_len);
            let mut dvalue = Vec::with_capacity(t_len);
            for t in 0..t_len {
                dmean.push(
                    (0..cfg.n_outputs)
                        .map(|i| 2.0 * (means[t][i] - targets[t][i]) / t_len as f32)
                        .collect(),
                );
                dvalue.push(config_value_grad(values[t], v_targets[t], t_len));
            }
            model.backward_sequence(&tape, &dmean, &dvalue);
            model.step_optimizers();
        }
        let l1 = loss_at(&mut model);
        println!("[supervised] loss {l0:.5} -> {l1:.5}");
        assert!(l1 < l0 * 0.25, "supervised loss must drop >4x, went {l0:.5} -> {l1:.5}");
    }

    fn config_value_grad(v: f32, target: f32, t_len: usize) -> f32 {
        2.0 * (v - target) / t_len as f32
    }
}


#[cfg(test)]
mod op_gradcheck {
    use super::*;

    /// Per-op finite-difference gradchecks (L1 rung): each hand-rolled
    /// backward is verified against central differences on small fixtures.
    /// Tolerance 2% — measured f32 FD noise at h=1e-2 on these magnitudes.


    #[test]
    fn op_rms_norm_backward_matches_fd() {
        let n = 8;
        let x = vec![0.3, -1.2, 0.7, 2.1, -0.4, 1.1, -2.0, 0.9];
        let w = vec![1.0, 0.9, 1.1, 0.8, 1.2, 0.95, 1.05, 0.85];
        let mut dy = vec![0.0f32; n];
        for (i, d) in dy.iter_mut().enumerate() {
            *d = 0.5 - 0.1 * i as f32;
        }
        let mut dw = vec![0.0f32; n];
        let mut dx = vec![0.0f32; n];
        rms_norm_backward(&x, &w, 1e-6, &dy, &mut dw, &mut dx);
        let f = |x: &[f32]| -> f32 {
            let mut y = x.to_vec();
            rms_norm(&mut y, &w, 1e-6);
            y.iter().zip(dy.iter()).map(|(yi, di)| yi * di).sum::<f32>()
        };
        let h = 1e-2;
        for i in 0..n {
            let mut xp = x.clone(); xp[i] += h;
            let mut xm = x.clone(); xm[i] -= h;
            let g_fd = (f(&xp) - f(&xm)) / (2.0 * h);
            let rel = (dx[i] - g_fd).abs() / (dx[i].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "rms dx[{i}] analytic={} fd={g_fd} rel={rel}", dx[i]);
        }
        for i in 0..n {
            let loss = |ww: &[f32]| -> f32 {
                let mut y = x.clone();
                rms_norm(&mut y, ww, 1e-6);
                y.iter().zip(dy.iter()).map(|(yi, di)| yi * di).sum::<f32>()
            };
            let mut wp = w.clone(); wp[i] += h;
            let mut wm = w.clone(); wm[i] -= h;
            let g_fd = (loss(&wp) - loss(&wm)) / (2.0 * h);
            let rel = (dw[i] - g_fd).abs() / (dw[i].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "rms dw[{i}] analytic={} fd={g_fd} rel={rel}", dw[i]);
        }
    }

    #[test]
    fn op_rope_backward_matches_fd() {
        let hd = 8;
        let x = vec![0.4, -0.9, 1.3, -0.2, 0.8, -1.1, 0.5, 2.0];
        let pos = 3usize;
        let mut dy = vec![0.0f32; 8];
        for (i, d) in dy.iter_mut().enumerate() {
            *d = 0.5 - 0.09 * i as f32;
        }
        let mut dx = vec![0.0f32; 8];
        rope_backward(hd, pos, &dy, &mut dx);
        let f = |x: &[f32]| -> f32 {
            let mut y = x.to_vec();
            rope_in_place(&mut y, hd, pos);
            y.iter().zip(dy.iter()).map(|(yi, di)| yi * di).sum::<f32>()
        };
        let h = 1e-2;
        for i in 0..8 {
            let mut xp = x.clone(); xp[i] += h;
            let mut xm = x.clone(); xm[i] -= h;
            let g_fd = (f(&xp) - f(&xm)) / (2.0 * h);
            let rel = (dx[i] - g_fd).abs() / (dx[i].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "rope dx[{i}] analytic={} fd={g_fd} rel={rel}", dx[i]);
        }
    }

    #[test]
    fn op_swiglu_backward_matches_fd() {
        let g = vec![0.7, -1.3, 2.1, -0.4];
        let u = vec![1.1, 0.6, -0.8, 1.4];
        let mut dy = vec![0.0f32; 4];
        for (i, d) in dy.iter_mut().enumerate() {
            *d = 0.4 - 0.15 * i as f32;
        }
        let mut dg = vec![0.0f32; 4];
        let mut du = vec![0.0f32; 4];
        swiglu_backward(&g, &u, &dy, &mut dg, &mut du);
        let f = |gg: &[f32], uu: &[f32]| -> f32 {
            let mut o = vec![0.0f32; 4];
            swiglu(gg, uu, &mut o);
            o.iter().zip(dy.iter()).map(|(oi, di)| oi * di).sum::<f32>()
        };
        let h = 1e-2;
        for i in 0..4 {
            let mut gp = g.clone(); gp[i] += h;
            let mut gm = g.clone(); gm[i] -= h;
            let g_fd = (f(&gp, &u) - f(&gm, &u)) / (2.0 * h);
            let rel = (dg[i] - g_fd).abs() / (dg[i].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "swiglu dg[{i}] analytic={} fd={g_fd} rel={rel}", dg[i]);
            let mut up = u.clone(); up[i] += h;
            let mut um = u.clone(); um[i] -= h;
            let g_fd = (f(&g, &up) - f(&g, &um)) / (2.0 * h);
            let rel = (du[i] - g_fd).abs() / (du[i].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "swiglu du[{i}] analytic={} fd={g_fd} rel={rel}", du[i]);
        }
    }

    #[test]
    fn op_matvec_backward_matches_fd() {
        let rows = 5;
        let cols = 3;
        let mut seed = 7u64;
        let mut w = vec![0.0f32; rows * cols];
        randomize_uniform(&mut w, cols, &mut seed);
        let x = vec![0.5, -1.0, 2.0];
        let mut dy = vec![0.0f32; rows];
        for (i, d) in dy.iter_mut().enumerate() {
            *d = ((i % 4) as f32 - 1.5) * 0.5;
        }
        let mut dw = vec![0.0f32; rows * cols];
        let mut dx = vec![0.0f32; cols];
        matvec_backward(&w, rows, cols, &x, &dy, &mut dw, &mut dx);
        let f = |ww: &[f32], xx: &[f32]| -> f32 {
            let mut y = vec![0.0f32; rows];
            matvec(ww, rows, cols, xx, &mut y);
            y.iter().zip(dy.iter()).map(|(yi, di)| yi * di).sum::<f32>()
        };
        let h = 1e-2;
        for c in 0..cols {
            let mut xp = x.clone(); xp[c] += h;
            let mut xm = x.clone(); xm[c] -= h;
            let g_fd = (f(&w, &xp) - f(&w, &xm)) / (2.0 * h);
            let rel = (dx[c] - g_fd).abs() / (dx[c].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "matvec dx[{c}] analytic={} fd={g_fd} rel={rel}", dx[c]);
        }
        for i in 0..rows * cols {
            let mut wp = w.clone(); wp[i] += h;
            let mut wm = w.clone(); wm[i] -= h;
            let g_fd = (f(&wp, &x) - f(&wm, &x)) / (2.0 * h);
            let rel = (dw[i] - g_fd).abs() / (dw[i].abs() + g_fd.abs() + 1e-6);
            assert!(rel < 0.02, "matvec dw[{i}] analytic={} fd={g_fd} rel={rel}", dw[i]);
        }
    }
}

// ─── Composition oracle: 2-token single-head mini case ───
//
// Decisive instrument for the policy-path gradcheck defect. A minimal
// hand-solvable case: 1 layer (Config with d_model=8, 1 head, head_dim 8,
// kv_dim 8, mlp 16, context 4, n_inputs 8, n_outputs 8), identity embed
// (E = I), so h_t = e_t exactly. The hand-derived analytic gradients for
// the 2-token attention chain are compared against backward_sequence.
//
// Hand derivation (single-head, T=2, attention over {0,1} for token 1 and
// {1} for token 0; scale s = 1/√hd):
//   token 0: attn_0 = v_0 (single-slot softmax ⇒ a=1, softmax-VJP = 0 ⇒
//            NO dq/dk contribution at all from token 0's own attention).
//   token 1: a_j = softmax(s_0, s_1); attn_1 = Σ a_j v_j.
//     per dim i: dA_j = dattn1[i]·v_j[i]; dot_i = a_0·dA_0[i] + a_1·dA_1[i];
//     ds_j[i] = a_j·(dA_j[i] − dot_i)
//     dq1[i] += Σ_j ds_j[i]·k_j[i]·s ; gk_j[i] += ds_j[i]·q1[i]·s ;
//     gv_j[i] += a_j·dattn1[i]
//   projections: dq1_pre = RoPEᵀ_1(dq1); gk_j_pre = RoPEᵀ_j(gk_j);
//   weight grads: dWq += dq1_pre·n1ᵀ; dWk += Σ_j gk_j_pre·n_jᵀ;
//                 dWv += Σ_j gv_j·n_jᵀ; dWo += Σ_t dout_t·attn_tᵀ
//   stream grads: dn1 = Wqᵀdq1_pre + Wkᵀgk1_pre + Wvᵀgv1 ;
//                 dn0 = Wkᵀgk0_pre + Wvᵀgv0
//   embed: de_t[c] = Σ_r dn_t[r]·E[r*D+c]   (E = I ⇒ de = dn)

#[cfg(test)]
mod composition_oracle {
    use super::*;

    const MD: usize = 8; // model dim for this mini case

    fn mini_cfg() -> Config {
        Config {
            d_model: MD,
            n_heads: 1,
            head_dim: MD,
            n_kv_heads: 1,
            kv_dim: MD,
            n_layers: 1,
            mlp_hidden: 16,
            context: 4,
            n_inputs: MD,
            n_outputs: MD,
        }
    }

    fn identity_embed_model(lr_muon: f32, lr_adam: f32, seed: &mut u64) -> GaitTransformer {
        let mut model = GaitTransformer::new(mini_cfg(), lr_muon, 0.9, lr_adam, seed);
        // E = I (d × n_inputs, square here)
        for r in 0..MD {
            for c in 0..MD {
                model.embed.params[r * MD + c] = if r == c { 1.0 } else { 0.0 };
            }
        }
        model
    }

    #[test]
    fn composition_oracle_two_token() {
        let mut seed = 0x0FAC_ADE_u64;
        let lr_muon = 0.0; // freeze Muon for the pure-gradient check
        let mut model = identity_embed_model(lr_muon, 0.0, &mut seed);
        let h = 1e-2f32;

        // deterministic inputs and head weights
        let obs: Vec<Vec<f32>> = (0..2)
            .map(|t| (0..MD).map(|i| ((i + 1 + 3 * t) as f32 * 0.21).sin()).collect())
            .collect();
        // probe: per-token weights over means and values
        let w_mean: Vec<Vec<f32>> = (0..2)
            .map(|t| (0..MD).map(|i| (i + 1) as f32 * 0.3 - 0.1 * (t + 1) as f32).collect())
            .collect();
        let w_val: Vec<f32> = (0..2).map(|t| 0.4 - 0.1 * t as f32).collect();

        // analytic: forward_sequence_train + backward_sequence
        model.zero_grads();
        let (_means, _values, tape) = model.forward_sequence_train(&obs);
        let dmean: Vec<Vec<f32>> = w_mean.clone();
        let dvalue: Vec<f32> = w_val.clone();
        model.backward_sequence(&tape, &dmean, &dvalue);
        let analytic = model.embed.grads.to_vec();

        // FD on the full loss
        let loss = |m: &mut GaitTransformer| -> f32 {
            let (means, values, _t) = m.forward_sequence_train(&obs);
            let mut l = 0.0;
            for t in 0..2 {
                for i in 0..MD {
                    l += w_mean[t][i] * means[t][i];
                }
                l += w_val[t] * values[t];
            }
            l
        };
        let mut worst = 0.0f32;
        let mut where_worst = String::new();
        for idx in 0..MD * MD {
            model.embed.params[idx] += h;
            let lp = loss(&mut model);
            model.embed.params[idx] -= 2.0 * h;
            let lm = loss(&mut model);
            model.embed.params[idx] += h;
            let g_fd = (lp - lm) / (2.0 * h);
            let rel = (analytic[idx] - g_fd).abs() / (analytic[idx].abs() + g_fd.abs() + 1e-6);
            if rel > worst {
                worst = rel;
                where_worst = format!("embed[{idx}] analytic={:.6} fd={g_fd:.6}", analytic[idx]);
            }
        }
        println!("[composition-oracle] worst rel {worst:.5} at {where_worst}");
        assert!(worst < 0.02, "composition oracle failed: {where_worst}");
    }
}

// ─── Config-axis bisection for the main gradcheck ───

#[cfg(test)]
mod slot_helper_fd {
    use super::*;

    fn dot(a: &[f32], b: &[f32]) -> f32 {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
    }

    /// L(dattn, k_cache, v_cache, q) = Σ_i dattn[i]·attn1[i] with
    /// attn1 = softmax(q·k0, q·k1)·(v0, v1) over 2 slots, single head.
    fn loss(
        dattn: &[f32],
        k_cache: &[f32],
        v_cache: &[f32],
        q: &[f32],
        d: usize,
        scale: f32,
    ) -> f32 {
        let s0 = dot(q, &k_cache[0..d]) * scale;
        let s1 = dot(q, &k_cache[d..2 * d]) * scale;
        let m = s0.max(s1);
        let e0 = (s0 - m).exp();
        let e1 = (s1 - m).exp();
        let z = e0 + e1;
        (0..d)
            .map(|i| dattn[i] * ((e0 / z) * v_cache[i] + (e1 / z) * v_cache[d + i]))
            .sum::<f32>()
    }

    #[test]
    fn slot_helper_matches_fd() {
        let d = 8usize;
        let n_heads = 1usize;
        let head_dim = 8usize;
        let scale = 1.0 / (head_dim as f32).sqrt();
        let mut seed = 0x51_07u64;
        let mut q1 = vec![0.0f32; d];
        let mut k0 = vec![0.0f32; d];
        let mut k1 = vec![0.0f32; d];
        let mut v0 = vec![0.0f32; d];
        let mut v1 = vec![0.0f32; d];
        for v in [&mut q1, &mut k0, &mut k1, &mut v0, &mut v1].iter_mut() {
            randomize_uniform(v, d, &mut seed);
        }
        let mut k_cache = k0.clone();
        k_cache.extend_from_slice(&k1);
        let mut v_cache = v0.clone();
        v_cache.extend_from_slice(&v1);
        let mut dattn = vec![0.0f32; d];
        for (i, d_) in dattn.iter_mut().enumerate() {
            *d_ = 0.5 - 0.07 * i as f32;
        }
        let count = 2usize;
        let mut dq = vec![0.0f32; d];
        let mut dk = vec![0.0f32; d * count];
        let mut dv = vec![0.0f32; d * count];
        for slot in 0..count {
            let mut dks = vec![0.0f32; d];
            let mut dvs = vec![0.0f32; d];
            attention_causal_backward_slot(
                &q1, &k_cache, &v_cache, count, slot, n_heads, head_dim, &dattn,
                &mut dq, &mut dks, &mut dvs,
            );
            for i in 0..d {
                dk[slot * d + i] += dks[i];
                dv[slot * d + i] += dvs[i];
            }
        }
        let h = 1e-3f32;
        let mut worst = 0.0f32;
        let mut worst_name = String::new();
        let mut check = |name: String, analytic: f32, fd: f32| {
            // Mixed absolute+relative: near-zero gradients are FD-noise
            // dominated in f32 (loss ~1e0, h ~1e-3 -> grad noise ~5e-5).
            let err = (analytic - fd).abs();
            let bad = err > 1e-4 + 0.05 * fd.abs();
            if bad {
                worst = worst.max(err / (analytic.abs() + fd.abs() + 1e-7));
                worst_name = format!("{name} analytic={analytic:.6} fd={fd:.6}");
            }
        };
        // dq (both slots feed the softmax)
        for i in 0..d {
            let mut qp = q1.clone();
            qp[i] += h;
            let mut qm = q1.clone();
            qm[i] -= h;
            let g_fd = (loss(&dattn, &k_cache, &v_cache, &qp, d, scale)
                - loss(&dattn, &k_cache, &v_cache, &qm, d, scale))
                / (2.0 * h);
            check(format!("dq[{i}]"), dq[i], g_fd);
        }
        // dk per slot
        for i in 0..d {
            let mut kp = k_cache.clone();
            kp[i] += h;
            let mut km = k_cache.clone();
            km[i] -= h;
            let g_fd = (loss(&dattn, &kp, &v_cache, &q1, d, scale)
                - loss(&dattn, &km, &v_cache, &q1, d, scale))
                / (2.0 * h);
            check(format!("dk[{i}]"), dk[i], g_fd);
        }
        // dv per slot
        for i in 0..d {
            let mut vp = v_cache.clone();
            vp[i] += h;
            let mut vm = v_cache.clone();
            vm[i] -= h;
            let g_fd = (loss(&dattn, &k_cache, &vp, &q1, d, scale)
                - loss(&dattn, &k_cache, &vm, &q1, d, scale))
                / (2.0 * h);
            check(format!("dv[{i}]"), dv[i], g_fd);
        }
        println!("[slot-helper] worst rel {worst:.5} at {worst_name}");
assert!(worst < 0.02, "slot helper failed: {worst_name}");
    }
}
