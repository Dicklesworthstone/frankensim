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
/// each input (same arity, same lengths — checked by the sweep).
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
            op: "leaf".to_string(),
            inputs: Vec::new(),
            value,
        });
        self.nodes.len() - 1
    }

    /// Record an op application (the value was computed by the caller's
    /// forward code — the tape only remembers structure + primals).
    pub fn apply(&mut self, op: &str, inputs: &[usize], value: Vec<f64>) -> usize {
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
    pub fn transpose(
        &self,
        registry: &VjpRegistry,
        output: usize,
        seed: &[f64],
    ) -> Result<BTreeMap<usize, Vec<f64>>, TransposeError> {
        let mut cotangents: Vec<Option<Vec<f64>>> = vec![None; self.nodes.len()];
        cotangents[output] = Some(seed.to_vec());
        for id in (0..self.nodes.len()).rev() {
            let Some(bar) = cotangents[id].clone() else {
                continue;
            };
            let node = &self.nodes[id];
            if node.op == "leaf" {
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
            for (&src, ib) in node.inputs.iter().zip(input_bars) {
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
            if node.op == "leaf"
                && let Some(g) = &cotangents[id]
            {
                grads.insert(id, g.clone());
            }
        }
        Ok(grads)
    }
}

/// Transpose-consistency check `max |⟨Av, w⟩ − ⟨v, Aᵀw⟩|` over seeded
/// deterministic probes — the G0 suite every registered linear op runs.
#[must_use]
pub fn check_transpose(
    apply: &dyn Fn(&[f64]) -> Vec<f64>,
    apply_t: &dyn Fn(&[f64]) -> Vec<f64>,
    n_in: usize,
    n_out: usize,
    probes: usize,
) -> f64 {
    let mut state = 0x7ea5_e11e_u64;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    let mut worst = 0.0f64;
    for _ in 0..probes {
        let v: Vec<f64> = (0..n_in).map(|_| lcg()).collect();
        let w: Vec<f64> = (0..n_out).map(|_| lcg()).collect();
        let av = apply(&v);
        let atw = apply_t(&w);
        let lhs: f64 = av.iter().zip(&w).map(|(a, b)| a * b).sum();
        let rhs: f64 = v.iter().zip(&atw).map(|(a, b)| a * b).sum();
        worst = worst.max((lhs - rhs).abs());
    }
    worst
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
    /// True when the adjoint agrees within the scaled tolerance.
    pub consistent: bool,
}

/// Compare an adjoint directional derivative against central finite
/// differences with a conditioning-aware tolerance: the FD self-error
/// |FD(h) − FD(h/2)| estimates how much the seam itself wobbles, and
/// the acceptance band is `base_tol · scale + 3·self_error`.
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
        (f(&xp) - f(&xm)) / (2.0 * step)
    };
    let fd_coarse = eval(h);
    let fd_fine = eval(h / 2.0);
    let self_error = (fd_coarse - fd_fine).abs();
    let scale = adjoint_dd.abs().max(fd_fine.abs()).max(1.0);
    let tolerance = base_tol * scale + 3.0 * self_error;
    FdVerdict {
        adjoint_dd,
        fd_coarse,
        fd_fine,
        tolerance,
        consistent: (adjoint_dd - fd_fine).abs() <= tolerance,
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
    b.chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().expect("8 bytes")))
        .collect()
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
