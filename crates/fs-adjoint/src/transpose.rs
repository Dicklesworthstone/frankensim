//! TRANSPOSE THE LEDGER (addendum Proposal 1, bead bk0o.1; [F] — behind
//! the `ledger-transpose` feature until its Gauntlet tier + kill metric
//! are green): the ledger composes error FORWARD; the same DAG,
//! transposed, composes sensitivity BACKWARD — ∂(lift)/∂(control point)
//! THROUGH the conversion, THROUGH the mesh, THROUGH the solve. Per-op
//! adjoints exist today and die at every seam; a ledger-shaped system
//! gets the chain almost by transposition, because restriction and
//! conversion maps are linear operators whose adjoints are free.
//!
//! The AMENDMENT this module enforces: every op REGISTERS a VJP or an
//! explicit non-differentiable declaration with color consequences. A
//! missing VJP inside a differentiation path is a STRUCTURED, LOUD
//! error that blocks the gradient — never a silent zero (a silently
//! zero seam gradient is a Goodhart trap).
//!
//! Boundary vs the base crate: fs-adjoint's other modules own per-op
//! discrete adjoints (IFT, revolve, Hadamard); this module owns ONLY
//! the DAG transposition that chains those VJPs across seams, plus the
//! content-addressed checkpoint SPILL contract shared with Proposal 2's
//! store discipline.

use std::collections::BTreeMap;
use std::sync::Arc;

/// One op's vector-Jacobian product: given the primal inputs it saw and
/// the cotangent arriving at its output, produce the cotangents for
/// each input.
///
/// SHAPE CONTRACT, enforced by [`Tape::transpose`]: the returned vector
/// has the op's arity, and cotangent `i` has exactly the length of the
/// recorded primal value of input `i`. Both are checked on every sweep
/// and a violation PANICS — the accumulation used to `zip`, which
/// truncates to the shorter vector and turns a masking/slicing bug into
/// a partly-dropped gradient indistinguishable from a genuine zero
/// cotangent (the silent-zero Goodhart trap this module exists to
/// prevent). The seed cotangent handed to `transpose` is checked the
/// same way against the output node's recorded value.
pub trait Vjp: Send + Sync {
    /// Pull the output cotangent back through the op.
    fn vjp(&self, primal_inputs: &[&[f64]], out_cotangent: &[f64]) -> Vec<Vec<f64>>;
}

/// Registry entry: differentiable, or DECLARED non-differentiable with
/// the color consequence spelled out.
#[derive(Clone)]
pub enum VjpEntry {
    /// A registered VJP.
    Differentiable(Arc<dyn Vjp>),
    /// An explicit refusal: gradients through this op are blocked, and
    /// downstream claims degrade to the named color at best.
    NonDifferentiable {
        /// Why (teaches the caller).
        reason: String,
        /// The color consequence (e.g. "estimated at best").
        color_consequence: String,
    },
}

impl std::fmt::Debug for VjpEntry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VjpEntry::Differentiable(_) => f.write_str("Differentiable"),
            VjpEntry::NonDifferentiable { reason, .. } => {
                write!(f, "NonDifferentiable({reason})")
            }
        }
    }
}

/// The per-op VJP registry (the op-spec amendment made executable).
#[derive(Debug, Default)]
pub struct VjpRegistry {
    entries: BTreeMap<String, VjpEntry>,
}

/// Structured transposition failures — loud, teaching, never silent.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TransposeError {
    /// An op in the differentiation path has NO registration at all.
    MissingVjp {
        /// The offending op kind.
        op: String,
    },
    /// An op in the path is declared non-differentiable.
    NonDifferentiableInPath {
        /// The offending op kind.
        op: String,
        /// The declared reason.
        reason: String,
        /// The declared color consequence.
        color_consequence: String,
    },
}

impl std::fmt::Display for TransposeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TransposeError::MissingVjp { op } => write!(
                f,
                "op '{op}' has no registered VJP and no non-differentiable declaration: \
                 the gradient is BLOCKED (register a VJP or declare the op, with color \
                 consequences — a silent zero here would be a Goodhart trap)"
            ),
            TransposeError::NonDifferentiableInPath {
                op,
                reason,
                color_consequence,
            } => write!(
                f,
                "op '{op}' is declared non-differentiable ({reason}); the gradient is \
                 blocked and downstream claims are {color_consequence}"
            ),
        }
    }
}

impl std::error::Error for TransposeError {}

impl VjpRegistry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        VjpRegistry::default()
    }

    /// Register an op's VJP.
    pub fn register(&mut self, op: &str, vjp: Arc<dyn Vjp>) {
        self.entries
            .insert(op.to_string(), VjpEntry::Differentiable(vjp));
    }

    /// Declare an op non-differentiable (the honest alternative).
    pub fn declare_non_differentiable(&mut self, op: &str, reason: &str, consequence: &str) {
        self.entries.insert(
            op.to_string(),
            VjpEntry::NonDifferentiable {
                reason: reason.to_string(),
                color_consequence: consequence.to_string(),
            },
        );
    }

    /// Coverage report: (registered, declared-non-differentiable) names.
    #[must_use]
    pub fn coverage(&self) -> (Vec<&str>, Vec<&str>) {
        let mut diff = Vec::new();
        let mut nondiff = Vec::new();
        for (k, v) in &self.entries {
            match v {
                VjpEntry::Differentiable(_) => diff.push(k.as_str()),
                VjpEntry::NonDifferentiable { .. } => nondiff.push(k.as_str()),
            }
        }
        (diff, nondiff)
    }

    fn lookup(&self, op: &str) -> Option<&VjpEntry> {
        self.entries.get(op)
    }
}

/// One recorded op application on the tape.
#[derive(Debug, Clone)]
pub struct TapeNode {
    /// The op kind (registry key).
    pub op: String,
    /// Input node ids (leaves are inputs pushed via [`Tape::leaf`]).
    pub inputs: Vec<usize>,
    /// The value this node produced.
    pub value: Vec<f64>,
}

/// The reserved op name marking a tape LEAF. [`Tape::transpose`] stops
/// the sweep at these nodes without consulting the registry, so it may
/// never name a real op (see [`Tape::apply`]).
const LEAF_OP: &str = "leaf";

/// The forward recording of a DAG execution.
#[derive(Debug, Default)]
pub struct Tape {
    nodes: Vec<TapeNode>,
}

impl Tape {
    /// Empty tape.
    #[must_use]
    pub fn new() -> Self {
        Tape::default()
    }

    /// Record a LEAF (an input the caller wants gradients for).
    pub fn leaf(&mut self, value: Vec<f64>) -> usize {
        self.nodes.push(TapeNode {
            op: LEAF_OP.to_string(),
            inputs: Vec::new(),
            value,
        });
        self.nodes.len() - 1
    }

    /// Record an op application (the value was computed by the caller's
    /// forward code — the tape only remembers structure + primals).
    ///
    /// # Panics
    /// If `op` is the reserved leaf sentinel `"leaf"`: [`Tape::transpose`]
    /// terminates the chain at any node whose op is that name WITHOUT a
    /// registry lookup, so recording a real op under it would silently
    /// drop the rest of the differentiation path instead of raising
    /// [`TransposeError::MissingVjp`].
    pub fn apply(&mut self, op: &str, inputs: &[usize], value: Vec<f64>) -> usize {
        assert_ne!(
            op, LEAF_OP,
            "'{LEAF_OP}' is the reserved leaf sentinel: an op recorded under it would terminate \
             the transposed sweep silently instead of raising MissingVjp"
        );
        self.nodes.push(TapeNode {
            op: op.to_string(),
            inputs: inputs.to_vec(),
            value,
        });
        self.nodes.len() - 1
    }

    /// A node's recorded value.
    #[must_use]
    pub fn value(&self, id: usize) -> &[f64] {
        &self.nodes[id].value
    }

    /// TRANSPOSE the DAG: pull `seed` (the cotangent at `output`) back
    /// to every leaf. Deterministic: reverse node order, accumulation
    /// in fixed index order — re-runs are bit-equal.
    ///
    /// # Errors
    /// [`TransposeError`] when any op on the path lacks a VJP or is
    /// declared non-differentiable — the gradient is blocked, loudly.
    ///
    /// # Panics
    /// If `seed` does not have the length of the output node's recorded
    /// value, or a registered VJP returns the wrong arity or a cotangent
    /// whose length differs from its input's recorded primal. These are
    /// shape-contract violations (see [`Vjp`]); accumulating them would
    /// truncate to the shorter vector and publish a partly-dropped
    /// gradient as if it were complete.
    pub fn transpose(
        &self,
        registry: &VjpRegistry,
        output: usize,
        seed: &[f64],
    ) -> Result<BTreeMap<usize, Vec<f64>>, TransposeError> {
        assert_eq!(
            seed.len(),
            self.nodes[output].value.len(),
            "seed cotangent length must match the output node's value"
        );
        let mut cotangents: Vec<Option<Vec<f64>>> = vec![None; self.nodes.len()];
        cotangents[output] = Some(seed.to_vec());
        for id in (0..self.nodes.len()).rev() {
            let Some(bar) = cotangents[id].clone() else {
                continue;
            };
            let node = &self.nodes[id];
            if node.op == LEAF_OP {
                continue;
            }
            let entry = registry
                .lookup(&node.op)
                .ok_or_else(|| TransposeError::MissingVjp {
                    op: node.op.clone(),
                })?;
            let vjp = match entry {
                VjpEntry::Differentiable(v) => v,
                VjpEntry::NonDifferentiable {
                    reason,
                    color_consequence,
                } => {
                    return Err(TransposeError::NonDifferentiableInPath {
                        op: node.op.clone(),
                        reason: reason.clone(),
                        color_consequence: color_consequence.clone(),
                    });
                }
            };
            let primal_inputs: Vec<&[f64]> = node
                .inputs
                .iter()
                .map(|&i| self.nodes[i].value.as_slice())
                .collect();
            let input_bars = vjp.vjp(&primal_inputs, &bar);
            assert_eq!(
                input_bars.len(),
                node.inputs.len(),
                "op '{}' VJP arity",
                node.op
            );
            for (input, (&src, ib)) in node.inputs.iter().zip(input_bars).enumerate() {
                // The LENGTH check the trait doc advertises. Without it the
                // accumulation below zips against the shorter vector and a
                // masked/sliced cotangent silently drops its tail — the
                // dropped components are indistinguishable from genuine
                // zeros in the returned map.
                assert_eq!(
                    ib.len(),
                    self.nodes[src].value.len(),
                    "op '{}' VJP cotangent length for input {input}",
                    node.op
                );
                match &mut cotangents[src] {
                    Some(acc) => {
                        for (a, b) in acc.iter_mut().zip(&ib) {
                            *a += b;
                        }
                    }
                    slot @ None => *slot = Some(ib),
                }
            }
        }
        let mut grads = BTreeMap::new();
        for (id, node) in self.nodes.iter().enumerate() {
            if node.op == LEAF_OP
                && let Some(g) = &cotangents[id]
            {
                grads.insert(id, g.clone());
            }
        }
        Ok(grads)
    }
}

/// Why a transpose-consistency probe request carries no evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TransposeProbeError {
    /// Zero probes were requested. The fold is seeded at `0.0`, so the
    /// call would return the STRONGEST possible score having evaluated
    /// neither operator.
    NoProbes,
    /// An operator dimension is empty: `⟨Av, w⟩` and `⟨v, Aᵀw⟩` are then
    /// both the empty sum for every probe and agree vacuously.
    EmptyDimension {
        /// Input dimension.
        n_in: usize,
        /// Output dimension.
        n_out: usize,
    },
    /// An operator returned a wrong-length image; the inner products
    /// would `zip`-truncate and compare a prefix.
    ImageLength {
        /// Probe index.
        probe: usize,
        /// `"apply"` or `"apply_transpose"`.
        operator: &'static str,
        /// Length the declared dimensions require.
        expected: usize,
        /// Length the operator returned.
        got: usize,
    },
    /// An operator produced a non-finite image or inner product, so the
    /// residual is not a number the caller can threshold.
    NonFinite {
        /// Probe index.
        probe: usize,
        /// `"apply"` or `"apply_transpose"`.
        operator: &'static str,
    },
}

impl std::fmt::Display for TransposeProbeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TransposeProbeError::NoProbes => f.write_str(
                "transpose-consistency check requires at least one probe: zero probes would \
                 return the perfect residual 0.0 without evaluating either operator",
            ),
            TransposeProbeError::EmptyDimension { n_in, n_out } => write!(
                f,
                "transpose-consistency check requires non-empty dimensions (n_in = {n_in}, \
                 n_out = {n_out}): empty sums agree vacuously"
            ),
            TransposeProbeError::ImageLength {
                probe,
                operator,
                expected,
                got,
            } => write!(
                f,
                "probe {probe}: '{operator}' returned {got} components where {expected} were \
                 declared; the inner product would compare a truncated prefix"
            ),
            TransposeProbeError::NonFinite { probe, operator } => write!(
                f,
                "probe {probe}: '{operator}' produced a non-finite value, so the transpose \
                 residual is not a threshold-able number"
            ),
        }
    }
}

impl std::error::Error for TransposeProbeError {}

/// The evidence a transpose-consistency check actually produced.
///
/// `max_abs_residual` alone is not thresholdable: it is an ABSOLUTE
/// quantity, so an operator whose entries are ~1e-20 scores ~1e-40
/// however wrong its transpose is. The probe count and the pairing
/// scale are reported so a caller's threshold can be relative.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransposeProbeReport {
    /// Probes actually evaluated (always ≥ 1).
    pub probes: usize,
    /// `max |⟨Av, w⟩ − ⟨v, Aᵀw⟩|` over those probes.
    pub max_abs_residual: f64,
    /// `max(|⟨Av, w⟩|, |⟨v, Aᵀw⟩|)` over the same probes — what the
    /// absolute residual must be read against.
    pub scale: f64,
    /// `max_abs_residual / scale`, or `+∞` when `scale` is `0.0` (both
    /// pairings vanished on every probe, so the check discriminated
    /// nothing).
    pub rel_residual: f64,
}

/// Transpose-consistency check `max |⟨Av, w⟩ − ⟨v, Aᵀw⟩|` over seeded
/// deterministic probes — the G0 suite every registered linear op runs.
///
/// # Errors
/// [`TransposeProbeError`] when the request cannot produce evidence
/// (no probes, an empty dimension) or an operator answers with a
/// wrong-length or non-finite image. Refusing beats returning the
/// perfect score `0.0` from an experiment that never ran.
pub fn check_transpose(
    apply: &dyn Fn(&[f64]) -> Vec<f64>,
    apply_t: &dyn Fn(&[f64]) -> Vec<f64>,
    n_in: usize,
    n_out: usize,
    probes: usize,
) -> Result<TransposeProbeReport, TransposeProbeError> {
    if probes == 0 {
        return Err(TransposeProbeError::NoProbes);
    }
    if n_in == 0 || n_out == 0 {
        return Err(TransposeProbeError::EmptyDimension { n_in, n_out });
    }
    let mut state = 0x7ea5_e11e_u64;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    let mut worst = 0.0f64;
    let mut scale = 0.0f64;
    for probe in 0..probes {
        let v: Vec<f64> = (0..n_in).map(|_| lcg()).collect();
        let w: Vec<f64> = (0..n_out).map(|_| lcg()).collect();
        let av = apply(&v);
        if av.len() != n_out {
            return Err(TransposeProbeError::ImageLength {
                probe,
                operator: "apply",
                expected: n_out,
                got: av.len(),
            });
        }
        let atw = apply_t(&w);
        if atw.len() != n_in {
            return Err(TransposeProbeError::ImageLength {
                probe,
                operator: "apply_transpose",
                expected: n_in,
                got: atw.len(),
            });
        }
        let lhs: f64 = av.iter().zip(&w).map(|(a, b)| a * b).sum();
        let rhs: f64 = v.iter().zip(&atw).map(|(a, b)| a * b).sum();
        if !lhs.is_finite() {
            return Err(TransposeProbeError::NonFinite {
                probe,
                operator: "apply",
            });
        }
        if !rhs.is_finite() {
            return Err(TransposeProbeError::NonFinite {
                probe,
                operator: "apply_transpose",
            });
        }
        worst = worst.max((lhs - rhs).abs());
        scale = scale.max(lhs.abs()).max(rhs.abs());
    }
    let rel_residual = if scale > 0.0 {
        worst / scale
    } else {
        f64::INFINITY
    };
    Ok(TransposeProbeReport {
        probes,
        max_abs_residual: worst,
        scale,
        rel_residual,
    })
}

/// The conditioning-aware FD falsifier verdict (review-round-3
/// hardening: an ill-conditioned seam where adjoint and FD legitimately
/// diverge must NOT fire a false falsifier hit).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FdVerdict {
    /// The adjoint directional derivative under test.
    pub adjoint_dd: f64,
    /// Central FD at step h.
    pub fd_coarse: f64,
    /// Central FD at step h/2 (the Richardson probe).
    pub fd_fine: f64,
    /// The conditioning-scaled tolerance actually used.
    pub tolerance: f64,
    /// The central-difference ROUNDING floor at the fine step,
    /// `ε·(|f(x+½hd)| + |f(x−½hd)|)/|h|`. A difference quotient cannot
    /// resolve derivatives smaller than the cancellation noise of its
    /// own numerator, so this is the smallest directional derivative
    /// the experiment can see at this step.
    pub noise_floor: f64,
    /// Whether the experiment DISCRIMINATES: the larger of
    /// `|adjoint_dd|` and `|fd_fine|` clears `noise_floor`. When false
    /// the honest verdict is "no signal at this step" — every value in
    /// the band, including `0.0` and the sign-flipped truth, satisfies
    /// the comparison, so agreement would be vacuous.
    pub falsifiable: bool,
    /// True when the probe was falsifiable AND the adjoint agrees
    /// within the scaled tolerance.
    pub consistent: bool,
}

/// Compare an adjoint directional derivative against central finite
/// differences with a conditioning-aware tolerance: the FD self-error
/// |FD(h) − FD(h/2)| estimates how much the seam itself wobbles, and
/// the acceptance band is `base_tol · scale + 3·self_error`.
///
/// `scale` is `max(|adjoint_dd|, |fd_fine|, noise_floor)` — homogeneous
/// in the derivative being tested. It deliberately does NOT carry an
/// absolute floor: a `max(…, 1.0)` term makes the band never tighter
/// than `base_tol` in absolute units, so on any problem whose true
/// directional derivative is smaller than `base_tol` every value in
/// `[−base_tol, base_tol]` passes — including a silently zero adjoint
/// and the exact sign flip the falsifier exists to catch.
pub fn fd_falsifier(
    f: &dyn Fn(&[f64]) -> f64,
    x: &[f64],
    dir: &[f64],
    adjoint_dd: f64,
    h: f64,
    base_tol: f64,
) -> FdVerdict {
    let eval = |step: f64| {
        let xp: Vec<f64> = x.iter().zip(dir).map(|(a, d)| a + step * d).collect();
        let xm: Vec<f64> = x.iter().zip(dir).map(|(a, d)| a - step * d).collect();
        let f_plus = f(&xp);
        let f_minus = f(&xm);
        (
            (f_plus - f_minus) / (2.0 * step),
            f_plus.abs() + f_minus.abs(),
        )
    };
    let (fd_coarse, _) = eval(h);
    let (fd_fine, magnitude) = eval(h / 2.0);
    let self_error = (fd_coarse - fd_fine).abs();
    // The quotient's own rounding floor: the numerator carries roughly
    // ε·(|f₊| + |f₋|) of absolute error and is divided by 2·(h/2) = h.
    // The floor is homogeneous under the paired rescaling (d → s·d,
    // h → h/s) that leaves the evaluation points fixed, exactly like the
    // sibling gate's `1e-12·‖d‖∞`.
    let noise_floor = f64::EPSILON * magnitude / h.abs();
    let signal = adjoint_dd.abs().max(fd_fine.abs());
    let falsifiable = signal > noise_floor;
    let scale = signal.max(noise_floor);
    let tolerance = base_tol * scale + 3.0 * self_error;
    FdVerdict {
        adjoint_dd,
        fd_coarse,
        fd_fine,
        tolerance,
        noise_floor,
        falsifiable,
        consistent: falsifiable && (adjoint_dd - fd_fine).abs() <= tolerance,
    }
}

/// The content-addressed checkpoint contract (shared storage discipline
/// with Proposal 2's incremental cache): `put` returns a stable key for
/// the bytes; `get` returns exactly those bytes.
pub trait CheckpointStore {
    /// Store bytes, returning the content key.
    fn put(&mut self, bytes: &[u8]) -> Vec<u8>;
    /// Fetch by key (panics on unknown keys — a checkpointing logic
    /// bug, not a runtime condition).
    fn get(&self, key: &[u8]) -> Vec<u8>;
}

/// The trivial in-memory store (the no-spill baseline).
#[derive(Debug, Default)]
pub struct MemStore {
    items: BTreeMap<Vec<u8>, Vec<u8>>,
    counter: u64,
}

impl CheckpointStore for MemStore {
    fn put(&mut self, bytes: &[u8]) -> Vec<u8> {
        self.counter += 1;
        let key = self.counter.to_le_bytes().to_vec();
        self.items.insert(key.clone(), bytes.to_vec());
        key
    }

    fn get(&self, key: &[u8]) -> Vec<u8> {
        self.items.get(key).expect("checkpoint present").clone()
    }
}

fn state_to_bytes(u: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(u.len() * 8);
    for v in u {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}

fn state_from_bytes(b: &[u8]) -> Vec<f64> {
    let (chunks, rest) = b.as_chunks::<8>();
    assert!(rest.is_empty(), "checkpoint bytes are whole f64s");
    chunks.iter().map(|c| f64::from_le_bytes(*c)).collect()
}

/// A uniform-checkpoint adjoint sweep with checkpoints SPILLED through
/// a [`CheckpointStore`]: states at every `every`-th step round-trip
/// through the store (bytes → key → bytes), segments are recomputed
/// from the fetched checkpoints, and the reverse sweep runs the same
/// deterministic step sequence — so gradients are BIT-EQUAL with or
/// without spill (the f64 ↔ bytes round-trip is exact).
///
/// Returns (gradient, checkpoints stored, forward step evaluations).
pub fn spilled_adjoint(
    u0: &[f64],
    steps: usize,
    every: usize,
    store: &mut dyn CheckpointStore,
    step_forward: &dyn Fn(&[f64]) -> Vec<f64>,
    step_reverse: &dyn Fn(&[f64]) -> Vec<f64>,
    terminal_seed: &dyn Fn(&[f64]) -> Vec<f64>,
) -> (Vec<f64>, usize, u64) {
    assert!(every >= 1, "checkpoint stride");
    // Forward: spill checkpoints at stride boundaries.
    let mut keys = Vec::new();
    let mut offsets = Vec::new();
    let mut u = u0.to_vec();
    let mut fwd_evals = 0u64;
    for k in 0..steps {
        if k % every == 0 {
            keys.push(store.put(&state_to_bytes(&u)));
            offsets.push(k);
        }
        u = step_forward(&u);
        fwd_evals += 1;
    }
    let mut bar = terminal_seed(&u);
    // Reverse by segments, newest checkpoint first.
    for (seg, key) in keys.iter().enumerate().rev() {
        let seg_start = offsets[seg];
        let seg_end = if seg + 1 < offsets.len() {
            offsets[seg + 1]
        } else {
            steps
        };
        // Recompute the segment's states from the SPILLED checkpoint.
        let mut states = Vec::with_capacity(seg_end - seg_start);
        let mut s = state_from_bytes(&store.get(key));
        for _ in seg_start..seg_end {
            states.push(s.clone());
            s = step_forward(&s);
            fwd_evals += 1;
        }
        // Reverse through the segment (states are not needed by the
        // linear reverse step here, but the recompute pattern is the
        // contract nonlinear steps rely on).
        for _ in (seg_start..seg_end).rev() {
            bar = step_reverse(&bar);
        }
        let _ = states;
    }
    (bar, keys.len(), fwd_evals)
}
