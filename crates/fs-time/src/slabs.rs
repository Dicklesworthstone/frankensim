//! TIME SLABS AS CELLS (addendum Proposal 4, bead bk0o.7; [F] — behind
//! the `time-slabs` feature): coupled multiphysics is glued in
//! SPACETIME — operator splitting hands data across temporal interfaces
//! exactly as patches hand traces across spatial ones. Making slabs
//! cells of the complex puts SPLITTING ERROR IN THE LEDGER: measurable,
//! localized to specific coupling steps, adaptively controllable.
//!
//! HONEST SCOPE (verbatim from the proposal): this measures and
//! controls splitting error; it does NOT claim coupling STABILITY —
//! added-mass FSI instabilities are analysis problems solved
//! per-coupling. The machinery guarantees consistent data TRANSFER
//! and visible defects, not convergent iteration.
//!
//! ACTIVATION DISCIPLINE: the controller is built, but
//! [`activation_report`] encodes the gate — when splitting error is
//! under 20% of the budget the recommendation is INSTRUMENT-ONLY
//! (measure it, don't control it).

/// The coupled two-field linear fixture: `x′ = −x + c(t)·y`,
/// `y′ = −2y + c(t)·x` — the canonical splitting testbed (the
/// commutator of the split operators is what the defect measures).
#[derive(Clone)]
pub struct CoupledFixture {
    /// Time-dependent coupling strength.
    pub coupling: fn(f64) -> f64,
}

impl CoupledFixture {
    /// One MONOLITHIC reference step over `[t0, t1]`: RK4 on the full
    /// system at `fine` substeps — the slab's monolithic residual
    /// reference.
    #[must_use]
    pub fn monolithic(&self, state: [f64; 2], t0: f64, t1: f64, fine: usize) -> [f64; 2] {
        let mut u = state;
        #[allow(clippy::cast_precision_loss)]
        let dt = (t1 - t0) / fine as f64;
        let f = |t: f64, u: [f64; 2]| -> [f64; 2] {
            let c = (self.coupling)(t);
            [-u[0] + c * u[1], -2.0 * u[1] + c * u[0]]
        };
        for k in 0..fine {
            #[allow(clippy::cast_precision_loss)]
            let t = t0 + k as f64 * dt;
            let k1 = f(t, u);
            let k2 = f(
                t + dt / 2.0,
                [u[0] + dt / 2.0 * k1[0], u[1] + dt / 2.0 * k1[1]],
            );
            let k3 = f(
                t + dt / 2.0,
                [u[0] + dt / 2.0 * k2[0], u[1] + dt / 2.0 * k2[1]],
            );
            let k4 = f(t + dt, [u[0] + dt * k3[0], u[1] + dt * k3[1]]);
            for i in 0..2 {
                u[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
        }
        u
    }

    /// One LIE-SPLIT coupling step over `[t0, t1]` with `subcycles`
    /// sub-slabs: within each sub-slab, field x advances with y FROZEN
    /// at the interface value (the temporal trace handoff), then y
    /// advances with the new x frozen — exactly the splitting used when
    /// two physics codes exchange boundary data once per step.
    #[must_use]
    pub fn split_step(&self, state: [f64; 2], t0: f64, t1: f64, subcycles: usize) -> [f64; 2] {
        let mut u = state;
        #[allow(clippy::cast_precision_loss)]
        let dt = (t1 - t0) / subcycles as f64;
        for k in 0..subcycles {
            #[allow(clippy::cast_precision_loss)]
            let ta = t0 + k as f64 * dt;
            let c = (self.coupling)(ta + dt / 2.0);
            // Field x: x′ = −x + c·y_frozen (exact 1-D solve).
            let y_frozen = u[1];
            let x_new = (u[0] - c * y_frozen) * (-dt).exp() + c * y_frozen;
            // Field y: y′ = −2y + c·x_frozen (x handed across).
            let x_frozen = x_new;
            let y_new = (u[1] - c * x_frozen / 2.0) * (-2.0 * dt).exp() + c * x_frozen / 2.0;
            u = [x_new, y_new];
        }
        u
    }
}

/// One time slab — a CELL of the temporal complex — with its ledgered
/// coupling defect.
#[derive(Debug, Clone, PartialEq)]
pub struct SlabEntry {
    /// Slab start.
    pub t0: f64,
    /// Slab end.
    pub t1: f64,
    /// Subcycles the controller used.
    pub subcycles: usize,
    /// THE TEMPORAL COCYCLE: ‖split − monolithic‖ over the slab — the
    /// defect of the coupling handoff against the monolithic residual.
    pub defect: f64,
}

/// The per-coupling-step ledger (Proposal 3's budget pie, pointed at
/// time).
#[derive(Debug, Clone, Default)]
pub struct SlabLedger {
    /// Entries in slab order.
    pub entries: Vec<SlabEntry>,
}

impl SlabLedger {
    /// Total splitting-error mass.
    #[must_use]
    pub fn total_defect(&self) -> f64 {
        self.entries.iter().map(|e| e.defect).sum()
    }

    /// THE BUDGET PIE: intervals ranked by defect share — "your error
    /// is in the coupling handoff at t ∈ [2.1, 2.3]".
    #[must_use]
    pub fn attribute(&self) -> Vec<(f64, f64, f64)> {
        let total = self.total_defect().max(1e-300);
        let mut shares: Vec<(f64, f64, f64)> = self
            .entries
            .iter()
            .map(|e| (e.t0, e.t1, e.defect / total))
            .collect();
        shares.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.0.total_cmp(&b.0)));
        shares
    }

    /// Structured line for the study log.
    #[must_use]
    pub fn to_json(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::from("{\"slabs\":[");
        for (k, e) in self.entries.iter().enumerate() {
            if k > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"t0\":{:.3},\"t1\":{:.3},\"sub\":{},\"defect\":{:.3e}}}",
                e.t0, e.t1, e.subcycles, e.defect
            );
        }
        let _ = write!(out, "],\"total\":{:.3e}}}", self.total_defect());
        out
    }
}

/// March `[0, t_end]` in `slabs` uniform slabs with FIXED subcycling,
/// measuring every slab's temporal cocycle against the monolithic
/// reference. Returns (final state, ledger).
#[must_use]
pub fn march_instrumented(
    fixture: &CoupledFixture,
    state0: [f64; 2],
    t_end: f64,
    slabs: usize,
    subcycles: usize,
) -> ([f64; 2], SlabLedger) {
    let mut u = state0;
    let mut ledger = SlabLedger::default();
    #[allow(clippy::cast_precision_loss)]
    let dt = t_end / slabs as f64;
    for k in 0..slabs {
        #[allow(clippy::cast_precision_loss)]
        let (t0, t1) = (k as f64 * dt, (k as f64 + 1.0) * dt);
        let split = fixture.split_step(u, t0, t1, subcycles);
        let mono = fixture.monolithic(u, t0, t1, 64);
        let defect = ((split[0] - mono[0]).powi(2) + (split[1] - mono[1]).powi(2)).sqrt();
        ledger.entries.push(SlabEntry {
            t0,
            t1,
            subcycles,
            defect,
        });
        u = split;
    }
    (u, ledger)
}

/// THE ADAPTIVE CONTROLLER: per slab, double the subcycles until the
/// temporal cocycle is under `tol` (cap 64) — tightening the coupling
/// exactly where the cocycle is large. Returns (final state, ledger,
/// total substeps spent — the cost the payoff test compares).
#[must_use]
pub fn march_adaptive(
    fixture: &CoupledFixture,
    state0: [f64; 2],
    t_end: f64,
    slabs: usize,
    tol: f64,
) -> ([f64; 2], SlabLedger, usize) {
    let mut u = state0;
    let mut ledger = SlabLedger::default();
    let mut spent = 0usize;
    #[allow(clippy::cast_precision_loss)]
    let dt = t_end / slabs as f64;
    for k in 0..slabs {
        #[allow(clippy::cast_precision_loss)]
        let (t0, t1) = (k as f64 * dt, (k as f64 + 1.0) * dt);
        let mono = fixture.monolithic(u, t0, t1, 64);
        let mut subcycles = 1usize;
        let (split, defect) = loop {
            let s = fixture.split_step(u, t0, t1, subcycles);
            let d = ((s[0] - mono[0]).powi(2) + (s[1] - mono[1]).powi(2)).sqrt();
            if d <= tol || subcycles >= 64 {
                break (s, d);
            }
            subcycles *= 2;
        };
        spent += subcycles;
        ledger.entries.push(SlabEntry {
            t0,
            t1,
            subcycles,
            defect,
        });
        u = split;
    }
    (u, ledger, spent)
}

/// The ACTIVATION verdict (the Proposal-4 gate as code).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Activation {
    /// Splitting error dominates: build/enable the controller.
    ControlJustified,
    /// Splitting error is a minor budget line: keep it INSTRUMENTED
    /// but uncontrolled — measure, don't build.
    InstrumentOnly,
}

/// Gate activation on demand: control is justified only when the
/// ledgered splitting error exceeds 20% of the total error budget.
#[must_use]
pub fn activation_report(ledger: &SlabLedger, total_error_budget: f64) -> (f64, Activation) {
    let fraction = ledger.total_defect() / total_error_budget.max(1e-300);
    let verdict = if fraction >= 0.2 {
        Activation::ControlJustified
    } else {
        Activation::InstrumentOnly
    };
    (fraction, verdict)
}
