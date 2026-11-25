//! Transactional fixed-point interface accelerators (bounded Aitken and vector
//! IQN-ILS) and controller memory stored in the System's two-generation store.
//!
//! Every workspace is pre-allocated at seal time so warm coupled steps perform
//! zero heap allocation. Accelerator and step-size controller state live in the
//! accepted generation, swap atomically at publish, roll back bitwise on trial
//! rejection or cancellation, and round-trip verbatim through `Checkpoint` and
//! `SystemImage`.

#![deny(unsafe_code)]

/// Maximum retained secant columns in an IQN-ILS block accelerator.
pub const MAX_IQN_HISTORY: usize = 32;

/// Configuration for bounded Aitken dynamic under-relaxation on torn port
/// variables of a looped block.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AitkenConfig {
    /// Initial relaxation factor `omega_0` in `(0, max_relaxation]`.
    pub initial_relaxation: f64,
    /// Lower bound on dynamic relaxation `omega_k` in `(0, initial_relaxation]`.
    pub min_relaxation: f64,
    /// Upper bound on dynamic relaxation `omega_k` in `[initial_relaxation, 2.0]`.
    pub max_relaxation: f64,
}

impl Default for AitkenConfig {
    fn default() -> Self {
        Self {
            initial_relaxation: 0.5,
            min_relaxation: 0.05,
            max_relaxation: 1.25,
        }
    }
}

impl AitkenConfig {
    /// Validate bounds.
    #[must_use]
    pub fn is_valid(self) -> bool {
        self.initial_relaxation.is_finite()
            && self.min_relaxation.is_finite()
            && self.max_relaxation.is_finite()
            && self.min_relaxation > 0.0
            && self.initial_relaxation >= self.min_relaxation
            && self.max_relaxation >= self.initial_relaxation
            && self.max_relaxation <= 2.0
    }
}

/// Configuration for bounded vector IQN-ILS (Anderson acceleration with
/// twice-orthogonalized modified Gram--Schmidt rank filtering) on torn port
/// variables of a looped block.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IqnIlsConfig {
    /// Maximum retained secant pairs in `1..=MAX_IQN_HISTORY`.
    pub max_history: usize,
    /// Drop a secant column when its orthogonalized norm is at most this
    /// fraction of its original norm, in `(0, 1)`.
    pub relative_rank_tolerance: f64,
    /// Fallback Aitken relaxation parameters used before an independent secant
    /// column is available.
    pub fallback: AitkenConfig,
    /// Whether to retain accepted secant columns across consecutive time steps
    /// (subject to generation rollback and snapshots).
    pub reuse_history_across_steps: bool,
}

impl Default for IqnIlsConfig {
    fn default() -> Self {
        Self {
            max_history: 8,
            relative_rank_tolerance: 1.0e-10,
            fallback: AitkenConfig::default(),
            reuse_history_across_steps: true,
        }
    }
}

impl IqnIlsConfig {
    /// Validate configuration parameters.
    #[must_use]
    pub fn is_valid(self) -> bool {
        (1..=MAX_IQN_HISTORY).contains(&self.max_history)
            && self.relative_rank_tolerance.is_finite()
            && self.relative_rank_tolerance > 0.0
            && self.relative_rank_tolerance < 1.0
            && self.fallback.is_valid()
    }
}

/// Pre-allocated scratch buffers for zero-allocation IQN-ILS QR least-squares
/// solves on one looped block.
#[derive(Debug, Clone, Default)]
pub(crate) struct IqnIlsScratch {
    residual: Vec<f64>,
    cand_dr: Vec<f64>,
    cand_dg: Vec<f64>,
    /// Q orthonormal basis columns, row-major `max_cols * dim`.
    q: Vec<f64>,
    /// Upper triangular R, column-major `max_cols * max_cols`.
    r: Vec<f64>,
    /// Scaled W (image difference) columns, `max_cols * dim`.
    w: Vec<f64>,
    /// Temporary column vector during twice-orthogonalized MGS.
    v_col: Vec<f64>,
    /// Temporary right-hand side and solution coefficients.
    rhs: Vec<f64>,
    coeffs: Vec<f64>,
    /// Per-column scale factors and active-column ordering.
    tri_col: Vec<f64>,
}

impl IqnIlsScratch {
    pub(crate) fn new(dim: usize, max_history: usize) -> Self {
        let cap = max_history.clamp(1, MAX_IQN_HISTORY);
        let max_cols = cap.min(dim.max(1));
        Self {
            residual: vec![0.0; dim],
            cand_dr: vec![0.0; dim],
            cand_dg: vec![0.0; dim],
            q: vec![0.0; max_cols * dim],
            r: vec![0.0; max_cols * max_cols],
            w: vec![0.0; max_cols * dim],
            v_col: vec![0.0; dim],
            rhs: vec![0.0; dim],
            coeffs: vec![0.0; max_cols],
            tri_col: vec![0.0; max_cols],
        }
    }

    fn ensure_capacity(&mut self, dim: usize, max_history: usize) {
        let cap = max_history.clamp(1, MAX_IQN_HISTORY);
        let max_cols = cap.min(dim.max(1));
        if self.residual.len() != dim || self.coeffs.len() < max_cols {
            *self = Self::new(dim, cap);
        }
    }
}

/// Transactional fixed-point accelerator state for one looped block.
///
/// Stored in both accepted and candidate generations of [`crate::System`].
/// On commit, accepted and candidate swap; on rejection or cancellation, the
/// accepted instance is bitwise untouched.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockAcceleratorState {
    dim: usize,
    max_history: usize,
    aitken_omega: f64,
    has_prev: bool,
    prev_residual: Vec<f64>,
    prev_image: Vec<f64>,
    /// Number of stored secant columns in `0..=max_history`, newest at index 0.
    history_len: usize,
    /// Residual difference columns `dR`, flat `max_history * dim`, newest first.
    history_v: Vec<f64>,
    /// Image difference columns `dG`, flat `max_history * dim`, newest first.
    history_w: Vec<f64>,
}

impl BlockAcceleratorState {
    /// Allocate fixed-capacity accelerator state for a block of `dim` unknowns.
    #[must_use]
    pub fn new(dim: usize, max_history: usize, initial_omega: f64) -> Self {
        let cap = max_history.min(MAX_IQN_HISTORY);
        Self {
            dim,
            max_history: cap,
            aitken_omega: initial_omega,
            has_prev: false,
            prev_residual: vec![0.0; dim],
            prev_image: vec![0.0; dim],
            history_len: 0,
            history_v: vec![0.0; cap * dim],
            history_w: vec![0.0; cap * dim],
        }
    }

    /// Unknown dimension of this block accelerator (`0` for non-looped blocks).
    #[must_use]
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// Current dynamic Aitken relaxation weight.
    #[must_use]
    pub fn relaxation_omega(&self) -> f64 {
        self.aitken_omega
    }

    /// Number of currently retained secant columns.
    #[must_use]
    pub fn history_len(&self) -> usize {
        self.history_len
    }

    /// Zero-allocation copy from another block accelerator of matching shape.
    pub fn copy_from(&mut self, other: &Self) {
        if self.dim != other.dim || self.max_history != other.max_history {
            *self = other.clone();
            return;
        }
        self.aitken_omega = other.aitken_omega;
        self.has_prev = other.has_prev;
        self.prev_residual.copy_from_slice(&other.prev_residual);
        self.prev_image.copy_from_slice(&other.prev_image);
        self.history_len = other.history_len;
        let used = other.history_len * other.dim;
        self.history_v[..used].copy_from_slice(&other.history_v[..used]);
        self.history_w[..used].copy_from_slice(&other.history_w[..used]);
    }

    /// Clear previous-iterate sample at the start of a new time step while
    /// optionally keeping cross-step secant history and dynamic relaxation.
    pub fn begin_step(&mut self, keep_history: bool, initial_omega: f64) {
        self.has_prev = false;
        if !keep_history {
            self.history_len = 0;
            self.aitken_omega = initial_omega;
        }
    }

    /// Perform one bounded Aitken dynamic under-relaxation update:
    /// `out = current + omega * (image - current)`.
    ///
    /// Leaves `self` untouched on non-finite arithmetic.
    pub(crate) fn aitken_step(
        &mut self,
        current: &[f64],
        image: &[f64],
        cfg: AitkenConfig,
        out: &mut [f64],
    ) -> Result<f64, &'static str> {
        let n = self.dim;
        if current.len() != n || image.len() != n || out.len() != n {
            return Err("aitken dimension mismatch");
        }
        if !cfg.is_valid() {
            return Err("invalid aitken configuration");
        }
        let mut omega = if self.aitken_omega.is_finite() && self.aitken_omega > 0.0 {
            self.aitken_omega
                .clamp(cfg.min_relaxation, cfg.max_relaxation)
        } else {
            cfg.initial_relaxation
        };

        // Compute residual norm and Aitken numerator/denominator without mutating self yet.
        // Scale prev_residual and dr by their maximum magnitude so `dr * dr` cannot overflow
        // or underflow for extreme physical units.
        let mut max_r = 0.0_f64;
        let mut pair_scale = 0.0_f64;
        for i in 0..n {
            let r = image[i] - current[i];
            if !r.is_finite() {
                return Err("non-finite interface residual");
            }
            max_r = max_r.max(r.abs());
            if self.has_prev {
                let dr = r - self.prev_residual[i];
                if !dr.is_finite() {
                    return Err("non-finite residual difference");
                }
                pair_scale = pair_scale.max(self.prev_residual[i].abs()).max(dr.abs());
            }
        }
        if max_r == 0.0 {
            out.copy_from_slice(image);
            return Ok(omega);
        }
        if self.has_prev && pair_scale > 0.0 {
            let mut num = 0.0_f64;
            let mut den = 0.0_f64;
            for i in 0..n {
                let r = image[i] - current[i];
                let r_prev_s = self.prev_residual[i] / pair_scale;
                let dr_s = (r - self.prev_residual[i]) / pair_scale;
                num += r_prev_s * dr_s;
                den += dr_s * dr_s;
            }
            if den.is_finite() && den > 0.0 && num.is_finite() {
                let candidate = -omega * (num / den);
                if candidate.is_finite() {
                    omega = candidate.clamp(cfg.min_relaxation, cfg.max_relaxation);
                }
            }
        }
        for i in 0..n {
            let r = image[i] - current[i];
            let next = current[i] + omega * r;
            if !next.is_finite() {
                return Err("non-finite aitken update");
            }
            out[i] = next;
        }
        // Commit sample to this candidate accelerator state only after all checks succeed.
        for i in 0..n {
            self.prev_residual[i] = image[i] - current[i];
            self.prev_image[i] = image[i];
        }
        self.has_prev = true;
        self.aitken_omega = omega;
        Ok(omega)
    }

    /// Perform one bounded vector IQN-ILS update with twice-orthogonalized MGS
    /// rank filtering, falling back to bounded Aitken relaxation when no
    /// independent secant column is available.
    ///
    /// Leaves `self` untouched on non-finite arithmetic.
    #[allow(clippy::too_many_lines)]
    pub(crate) fn iqn_ils_step(
        &mut self,
        current: &[f64],
        image: &[f64],
        cfg: IqnIlsConfig,
        scratch: &mut IqnIlsScratch,
        out: &mut [f64],
    ) -> Result<(usize, f64), &'static str> {
        let n = self.dim;
        if current.len() != n || image.len() != n || out.len() != n {
            return Err("iqn-ils dimension mismatch");
        }
        if !cfg.is_valid() {
            return Err("invalid iqn-ils configuration");
        }
        let max_hist = cfg.max_history.min(self.max_history).max(1);
        scratch.ensure_capacity(n, max_hist);

        let mut residual_scale = 0.0_f64;
        for i in 0..n {
            let r = image[i] - current[i];
            if !r.is_finite() {
                return Err("non-finite interface residual");
            }
            scratch.residual[i] = r;
            residual_scale = residual_scale.max(r.abs());
        }
        if residual_scale == 0.0 {
            out.copy_from_slice(image);
            return Ok((0, self.aitken_omega));
        }

        // Build candidate newest secant column if a previous sample exists.
        let mut has_candidate = false;
        let mut cand_scale = 0.0_f64;
        if self.has_prev {
            for i in 0..n {
                let dr = scratch.residual[i] - self.prev_residual[i];
                let dg = image[i] - self.prev_image[i];
                if !(dr.is_finite() && dg.is_finite()) {
                    return Err("non-finite secant difference");
                }
                scratch.cand_dr[i] = dr;
                scratch.cand_dg[i] = dg;
                cand_scale = cand_scale.max(dr.abs());
            }
            has_candidate = cand_scale > 0.0;
        }

        // Factor [candidate, history[0..]] with twice-orthogonalized MGS.
        let max_cols = max_hist.min(n);
        let total_available = usize::from(has_candidate)
            + self
                .history_len
                .min(max_hist.saturating_sub(usize::from(has_candidate)));
        let mut q_cols = 0usize;

        for idx in 0..total_available {
            if q_cols == max_cols {
                break;
            }
            let (dr_slice, dg_slice): (&[f64], &[f64]) = if has_candidate && idx == 0 {
                (&scratch.cand_dr[..n], &scratch.cand_dg[..n])
            } else {
                let h_idx = idx - usize::from(has_candidate);
                let start = h_idx * n;
                (
                    &self.history_v[start..start + n],
                    &self.history_w[start..start + n],
                )
            };
            let col_scale = dr_slice.iter().fold(0.0_f64, |acc, &v| acc.max(v.abs()));
            if !(col_scale.is_finite() && col_scale > 0.0) {
                continue;
            }
            let mut orig_norm2 = 0.0_f64;
            for i in 0..n {
                let v = dr_slice[i] / col_scale;
                scratch.v_col[i] = v;
                orig_norm2 += v * v;
            }
            let orig_norm = fs_math::det::sqrt(orig_norm2);
            if !(orig_norm.is_finite() && orig_norm > 0.0) {
                continue;
            }
            scratch.tri_col[..q_cols].fill(0.0);
            for _reorth in 0..2 {
                for j in 0..q_cols {
                    let qj = &scratch.q[j * n..(j + 1) * n];
                    let mut proj = 0.0_f64;
                    for i in 0..n {
                        proj += qj[i] * scratch.v_col[i];
                    }
                    if !proj.is_finite() {
                        return Err("non-finite MGS projection");
                    }
                    scratch.tri_col[j] += proj;
                    for i in 0..n {
                        scratch.v_col[i] -= proj * qj[i];
                    }
                }
            }
            let mut orth_norm2 = 0.0_f64;
            for &v in &scratch.v_col[..n] {
                orth_norm2 += v * v;
            }
            let orth_norm = fs_math::det::sqrt(orth_norm2);
            if !orth_norm.is_finite() || orth_norm <= cfg.relative_rank_tolerance * orig_norm {
                continue;
            }
            let q_dst = &mut scratch.q[q_cols * n..(q_cols + 1) * n];
            let w_dst = &mut scratch.w[q_cols * n..(q_cols + 1) * n];
            for i in 0..n {
                q_dst[i] = scratch.v_col[i] / orth_norm;
                let wv = dg_slice[i] / col_scale;
                if !wv.is_finite() {
                    return Err("non-finite scaled map difference");
                }
                w_dst[i] = wv;
            }
            for j in 0..q_cols {
                scratch.r[q_cols * max_cols + j] = scratch.tri_col[j];
            }
            scratch.r[q_cols * max_cols + q_cols] = orth_norm;
            q_cols += 1;
        }

        if q_cols == 0 {
            // No independent secant available yet -> bounded Aitken step.
            let omega = self.aitken_step(current, image, cfg.fallback, out)?;
            if has_candidate {
                self.push_secant(&scratch.cand_dr[..n], &scratch.cand_dg[..n], max_hist);
            }
            return Ok((0, omega));
        }

        // Solve R c = Q^T (-r / residual_scale).
        for i in 0..n {
            scratch.rhs[i] = -scratch.residual[i] / residual_scale;
        }
        for j in 0..q_cols {
            let qj = &scratch.q[j * n..(j + 1) * n];
            let mut sum = 0.0_f64;
            for i in 0..n {
                sum += qj[i] * scratch.rhs[i];
            }
            if !sum.is_finite() {
                return Err("non-finite IQN-ILS RHS dot");
            }
            scratch.coeffs[j] = sum;
        }
        for row in (0..q_cols).rev() {
            let mut val = scratch.coeffs[row];
            for col in (row + 1)..q_cols {
                val -= scratch.r[col * max_cols + row] * scratch.coeffs[col];
            }
            let diag = scratch.r[row * max_cols + row];
            let solved = val / diag;
            if !solved.is_finite() {
                return Err("non-finite IQN-ILS back substitution");
            }
            scratch.coeffs[row] = solved;
        }

        // Form x_next = G(x) + (W c) * residual_scale.
        for i in 0..n {
            let mut corr = 0.0_f64;
            for j in 0..q_cols {
                corr += scratch.w[j * n + i] * scratch.coeffs[j];
            }
            let next = image[i] + corr * residual_scale;
            if !next.is_finite() {
                return Err("non-finite IQN-ILS update");
            }
            out[i] = next;
        }

        // Commit candidate secant and latest sample only after full update succeeded.
        if has_candidate {
            self.push_secant(&scratch.cand_dr[..n], &scratch.cand_dg[..n], max_hist);
        }
        self.prev_residual.copy_from_slice(&scratch.residual[..n]);
        self.prev_image.copy_from_slice(image);
        self.has_prev = true;
        Ok((q_cols, 1.0))
    }

    fn push_secant(&mut self, dr: &[f64], dg: &[f64], max_hist: usize) {
        let n = self.dim;
        let cap = max_hist.min(self.max_history);
        if cap == 0 || n == 0 {
            return;
        }
        let keep = self.history_len.min(cap.saturating_sub(1));
        for col in (0..keep).rev() {
            let src = col * n;
            let dst = (col + 1) * n;
            self.history_v.copy_within(src..src + n, dst);
            self.history_w.copy_within(src..src + n, dst);
        }
        self.history_v[..n].copy_from_slice(dr);
        self.history_w[..n].copy_from_slice(dg);
        self.history_len = (keep + 1).min(cap);
    }

    /// Encode this block accelerator state into canonical `u64` words.
    pub fn push_words(&self, words: &mut Vec<u64>) {
        words.push(self.dim as u64);
        words.push(self.max_history as u64);
        words.push(self.aitken_omega.to_bits());
        words.push(u64::from(self.has_prev));
        if self.has_prev {
            for &v in &self.prev_residual {
                words.push(v.to_bits());
            }
            for &v in &self.prev_image {
                words.push(v.to_bits());
            }
        }
        words.push(self.history_len as u64);
        let used = self.history_len * self.dim;
        for &v in &self.history_v[..used] {
            words.push(v.to_bits());
        }
        for &v in &self.history_w[..used] {
            words.push(v.to_bits());
        }
    }

    /// Decode a block accelerator state from canonical `u64` words, verifying
    /// expected shape and finiteness.
    pub fn read_words(
        words: &[u64],
        cursor: &mut usize,
        expected_dim: usize,
        expected_max_history: usize,
    ) -> Option<Self> {
        let next = |cursor: &mut usize| -> Option<u64> {
            let word = *words.get(*cursor)?;
            *cursor += 1;
            Some(word)
        };
        let dim = usize::try_from(next(cursor)?).ok()?;
        let max_history = usize::try_from(next(cursor)?).ok()?;
        if dim != expected_dim || max_history != expected_max_history {
            return None;
        }
        let aitken_omega = f64::from_bits(next(cursor)?);
        if !(aitken_omega.is_finite() && aitken_omega > 0.0) {
            return None;
        }
        let has_prev_flag = next(cursor)?;
        if has_prev_flag > 1 {
            return None;
        }
        let has_prev = has_prev_flag == 1;
        let mut prev_residual = vec![0.0; dim];
        let mut prev_image = vec![0.0; dim];
        if has_prev {
            for slot in &mut prev_residual {
                let v = f64::from_bits(next(cursor)?);
                if !v.is_finite() {
                    return None;
                }
                *slot = v;
            }
            for slot in &mut prev_image {
                let v = f64::from_bits(next(cursor)?);
                if !v.is_finite() {
                    return None;
                }
                *slot = v;
            }
        }
        let history_len = usize::try_from(next(cursor)?).ok()?;
        if history_len > max_history {
            return None;
        }
        let used = history_len.checked_mul(dim)?;
        let cap_len = max_history.checked_mul(dim)?;
        let mut history_v = vec![0.0; cap_len];
        let mut history_w = vec![0.0; cap_len];
        for slot in &mut history_v[..used] {
            let v = f64::from_bits(next(cursor)?);
            if !v.is_finite() {
                return None;
            }
            *slot = v;
        }
        for slot in &mut history_w[..used] {
            let v = f64::from_bits(next(cursor)?);
            if !v.is_finite() {
                return None;
            }
            *slot = v;
        }
        Some(Self {
            dim,
            max_history,
            aitken_omega,
            has_prev,
            prev_residual,
            prev_image,
            history_len,
            history_v,
            history_w,
        })
    }
}

/// Complete numerical controller and accelerator state owned by the system's
/// two-generation store.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ControllerState {
    /// Step-size memory carried across accepted adaptive intervals [s].
    pub suggested_dt_s: Option<f64>,
    /// Per-block fixed-point accelerator state (length == number of blocks).
    pub blocks: Vec<BlockAcceleratorState>,
}

impl ControllerState {
    /// Zero-allocation copy from another controller state of matching layout.
    pub fn copy_from(&mut self, other: &Self) {
        self.suggested_dt_s = other.suggested_dt_s;
        if self.blocks.len() != other.blocks.len() {
            self.blocks.clone_from(&other.blocks);
            return;
        }
        for (dst, src) in self.blocks.iter_mut().zip(&other.blocks) {
            dst.copy_from(src);
        }
    }

    /// Encode into canonical `u64` words for durable snapshots.
    pub fn push_words(&self, words: &mut Vec<u64>) {
        match self.suggested_dt_s {
            Some(dt) => {
                words.push(1);
                words.push(dt.to_bits());
            }
            None => words.push(0),
        }
        words.push(self.blocks.len() as u64);
        for block in &self.blocks {
            block.push_words(words);
        }
    }

    /// Decode from canonical `u64` words, validating shape against `template`.
    #[must_use]
    pub fn read_words(&self, words: &[u64]) -> Option<Self> {
        let mut cursor = 0usize;
        let has_dt = *words.get(cursor)?;
        cursor += 1;
        let suggested_dt_s = match has_dt {
            0 => None,
            1 => {
                let dt = f64::from_bits(*words.get(cursor)?);
                cursor += 1;
                if !(dt.is_finite() && dt > 0.0) {
                    return None;
                }
                Some(dt)
            }
            _ => return None,
        };
        let count = usize::try_from(*words.get(cursor)?).ok()?;
        cursor += 1;
        if count != self.blocks.len() {
            return None;
        }
        let mut blocks = Vec::with_capacity(count);
        for template in &self.blocks {
            blocks.push(BlockAcceleratorState::read_words(
                words,
                &mut cursor,
                template.dim,
                template.max_history,
            )?);
        }
        if cursor != words.len() {
            return None;
        }
        Some(Self {
            suggested_dt_s,
            blocks,
        })
    }
}

/// Identity key for caches bound to one accepted system generation and step size.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StateCacheKey {
    /// Monotonic accepted-state generation counter of the owning [`crate::System`].
    pub accepted_generation: u64,
    /// Exact IEEE-754 bit pattern of the step duration `dt_s`.
    pub dt_bits: u64,
}

impl StateCacheKey {
    /// Construct a cache key from an accepted generation and step duration.
    #[must_use]
    pub fn new(accepted_generation: u64, dt_s: f64) -> Self {
        Self {
            accepted_generation,
            dt_bits: dt_s.to_bits(),
        }
    }
}

/// Single-slot cache keyed strictly by [`StateCacheKey`].
///
/// Any accepted-state change (step publish, `restore`, `resume`) advances
/// `accepted_generation`, invalidating stale entries automatically so cached
/// values are never observable across differing accepted states.
#[derive(Debug, Clone)]
pub struct AcceptedStateCache<V> {
    entry: Option<(StateCacheKey, V)>,
}

impl<V> Default for AcceptedStateCache<V> {
    fn default() -> Self {
        Self { entry: None }
    }
}

impl<V> AcceptedStateCache<V> {
    /// Empty cache.
    #[must_use]
    pub fn new() -> Self {
        Self { entry: None }
    }

    /// Borrow the cached value only if keyed by the exact `key`.
    #[must_use]
    pub fn get(&self, key: StateCacheKey) -> Option<&V> {
        self.entry
            .as_ref()
            .and_then(|(k, v)| (*k == key).then_some(v))
    }

    /// Store a value under `key`, replacing any prior entry.
    pub fn insert(&mut self, key: StateCacheKey, value: V) {
        self.entry = Some((key, value));
    }

    /// Drop any cached entry.
    pub fn invalidate(&mut self) {
        self.entry = None;
    }
}
