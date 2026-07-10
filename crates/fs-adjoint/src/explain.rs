//! EXPLANATION OBJECTS (addendum Proposal B, bead knh1.5; [F] — behind
//! the `explanation-objects` feature): when a QoI moves, return not a
//! NUMBER but a CAUSAL DECOMPOSITION along named physical channels,
//! each term with a bound, as a first-class, re-derivable, certified
//! artifact. The difference between confabulation and understanding is
//! whether the system can CHECK the story — so the explanation is an
//! OBJECT the system checks, and the natural-language rendering on top
//! is explicitly NON-AUTHORITATIVE.
//!
//! Three attribution engines feed one tree:
//! 1. ADJOINT attribution — for elliptic compliance the bilinear trick
//!    gives the EXACT identity `J₁ − J₀ = −∫ Δa·u₀′·u₁′`, so channel
//!    masks decompose ΔJ exactly (quadrature-level bounds).
//! 2. PROVENANCE attribution — which EDIT moved the number, by exact
//!    telescoping over replayed ledger states.
//! 3. PHYSICAL decomposition — the far-field drag FLAGSHIP: induced
//!    drag via the Trefftz-plane wake integral on a lifting-line
//!    fixture (reconciling with the analytic `C_L²/(π·AR)`), a
//!    viscous strip channel, and the wave channel DECLARED zero in the
//!    subsonic regime rather than silently omitted.
//!
//! THE HONESTY GATE (a permanent runtime invariant): if the
//! unattributed residual exceeds its threshold the system REFUSES to
//! explain rather than smearing the residual across plausible
//! channels. A partial explanation with a declared gap beats a
//! complete story with a hidden one.

use fs_evidence::Color;

/// Local deterministic FNV-1a fingerprint (keeps this L3 module free
/// of ledger-layer (L6) dependencies; the ledger stores the same hex).
fn fingerprint_hex(bytes: &[u8]) -> String {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{h:016x}")
}

fn assert_node_payload(
    channel: &str,
    contribution: f64,
    bound: f64,
    color: &Color,
    evidence: &[String],
) {
    assert!(!channel.is_empty(), "explanation channel must not be empty");
    assert!(
        contribution.is_finite(),
        "explanation contribution must be finite"
    );
    assert!(
        bound.is_finite() && bound >= 0.0,
        "explanation bound must be finite and non-negative"
    );
    assert!(
        !evidence.is_empty(),
        "explanation nodes must carry evidence links"
    );
    match color {
        Color::Verified { lo, hi } => {
            assert!(
                lo.is_finite() && hi.is_finite() && lo <= hi,
                "verified explanation color must carry finite ordered bounds"
            );
        }
        Color::Estimated { dispersion, .. } => {
            assert!(
                dispersion.is_finite() && *dispersion >= 0.0,
                "estimated explanation color must carry finite non-negative dispersion"
            );
        }
        Color::Validated { .. } => {}
    }
}

/// One attribution node: a named channel's contribution with its
/// bound, evidence links, and a re-derivation fingerprint.
#[derive(Debug, Clone, PartialEq)]
pub struct ExplanationNode {
    /// The named physical/provenance channel.
    pub channel: String,
    /// Signed contribution to the observed ΔQoI.
    pub contribution: f64,
    /// Certified half-width on the contribution.
    pub bound: f64,
    /// The evidence color of this term.
    pub color: Color,
    /// Ledger keys backing the term (solves, diffs, integrals).
    pub evidence: Vec<String>,
    /// Deterministic fingerprint of (channel, inputs) — the
    /// re-derivability witness: recomputing the node from the ledger
    /// must reproduce this exactly.
    pub fingerprint: String,
}

impl ExplanationNode {
    /// Build a node with its fingerprint derived from the payload.
    #[must_use]
    pub fn new(
        channel: &str,
        contribution: f64,
        bound: f64,
        color: Color,
        evidence: Vec<String>,
    ) -> ExplanationNode {
        use std::fmt::Write as _;

        assert_node_payload(channel, contribution, bound, &color, &evidence);
        let mut canon = String::new();
        let _ = write!(
            canon,
            "channel:{}:{channel};contrib:{:x};bound:{:x};color:{}:{};evidence:{}",
            channel.len(),
            contribution.to_bits(),
            bound.to_bits(),
            color.name(),
            color.payload_json(),
            evidence.len()
        );
        for item in &evidence {
            let _ = write!(canon, ";{}:{item}", item.len());
        }
        ExplanationNode {
            channel: channel.to_string(),
            contribution,
            bound,
            color,
            evidence,
            fingerprint: fingerprint_hex(canon.as_bytes()),
        }
    }
}

/// A finalized explanation, or the refusal that keeps it honest.
#[derive(Debug, Clone, PartialEq)]
pub enum Explanation {
    /// The tree reconciles: channels + declared residual = observed.
    Explained {
        /// Channel nodes.
        nodes: Vec<ExplanationNode>,
        /// The observed ΔQoI being explained.
        observed: f64,
        /// The declared unattributed residual (within threshold).
        residual: f64,
    },
    /// The residual exceeded its threshold: NO explanation is issued.
    /// The partial tree is returned as forensics, clearly not a claim.
    Refused {
        /// The partial (non-authoritative) attribution.
        partial: Vec<ExplanationNode>,
        /// The unattributed residual that triggered refusal.
        residual: f64,
        /// The threshold it exceeded.
        threshold: f64,
    },
}

impl Explanation {
    /// THE PERMANENT INVARIANT (the Proposal-B kill criterion): the
    /// certified channel bounds must COVER the unattributed residual.
    /// An engine failing this on any case is lying and ships nothing.
    ///
    /// The earlier form checked `attributed + residual == observed`,
    /// which `finalize` makes true BY CONSTRUCTION (residual is
    /// defined as observed − attributed) — a kill criterion that can
    /// never fire is not a kill criterion (bead 9sf6 F3). Coverage is
    /// the non-vacuous claim: bad channel math leaves a residual the
    /// bounds cannot absorb, and THIS fires.
    #[must_use]
    pub fn reconciles(&self) -> bool {
        match self {
            Explanation::Explained {
                nodes,
                observed,
                residual,
            } => {
                let attributed: f64 = nodes.iter().map(|n| n.contribution).sum();
                let bounds: f64 = nodes.iter().map(|n| n.bound).sum();
                // Construction sanity (cheap, and guards hand-built trees).
                let identity_ok = (attributed + residual - observed).abs() <= bounds.max(1e-14);
                let covered = residual.abs() <= bounds.max(1e-14);
                identity_ok && covered
            }
            Explanation::Refused { .. } => true, // a refusal claims nothing
        }
    }

    /// NON-AUTHORITATIVE natural-language rendering. The TREE is the
    /// artifact; this string is for humans skimming and says so.
    #[must_use]
    pub fn render_narrative(&self) -> String {
        use std::fmt::Write as _;
        let mut out =
            String::from("NON-AUTHORITATIVE RENDERING (the explanation tree is the artifact):\n");
        match self {
            Explanation::Explained {
                nodes,
                observed,
                residual,
            } => {
                let _ = writeln!(out, "observed change {observed:+.6e}");
                for n in nodes {
                    let _ = writeln!(
                        out,
                        "  {} contributed {:+.6e} (± {:.1e})",
                        n.channel, n.contribution, n.bound
                    );
                }
                let _ = writeln!(out, "  unattributed residual {residual:+.6e}");
            }
            Explanation::Refused {
                residual,
                threshold,
                ..
            } => {
                let _ = writeln!(
                    out,
                    "REFUSED: unattributed residual {residual:.3e} exceeds the honesty \
                     threshold {threshold:.3e}; no causal story is issued."
                );
            }
        }
        out
    }
}

/// Assemble + gate: compute the residual against the observed change
/// and REFUSE when it exceeds `threshold`.
#[must_use]
pub fn finalize(nodes: Vec<ExplanationNode>, observed: f64, threshold: f64) -> Explanation {
    assert!(observed.is_finite(), "observed change must be finite");
    assert!(
        threshold.is_finite() && threshold >= 0.0,
        "explanation threshold must be finite and non-negative"
    );
    for node in &nodes {
        assert_node_payload(
            &node.channel,
            node.contribution,
            node.bound,
            &node.color,
            &node.evidence,
        );
    }
    let attributed: f64 = nodes.iter().map(|n| n.contribution).sum();
    let residual = observed - attributed;
    if residual.abs() > threshold {
        Explanation::Refused {
            partial: nodes,
            residual,
            threshold,
        }
    } else {
        Explanation::Explained {
            nodes,
            observed,
            residual,
        }
    }
}

// ---------------------------------------------------------------------------
// Engine 1: ADJOINT attribution on the elliptic compliance fixture.
// ---------------------------------------------------------------------------

/// The 1-D elliptic fixture: `−(a u′)′ = 1`, u(0)=u(1)=0, P1 elements;
/// compliance `J = ∫ u`. Channel masks partition the elements.
#[derive(Debug, Clone)]
pub struct Elliptic1d {
    /// Interior nodes.
    pub n: usize,
}

impl Elliptic1d {
    /// Solve with per-element conductivity `a` (length n+1).
    #[must_use]
    pub fn solve(&self, a: &[f64]) -> Vec<f64> {
        let n = self.n;
        assert!(n > 0, "Elliptic1d requires at least one interior node");
        assert_eq!(
            a.len(),
            n + 1,
            "conductivity length must equal n + 1 elements"
        );
        assert!(
            a.iter().all(|v| v.is_finite() && *v > 0.0),
            "conductivity values must be finite and positive"
        );
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let mut diag = vec![0.0f64; n];
        let mut off = vec![0.0f64; n.saturating_sub(1)];
        for (e, &ae) in a.iter().enumerate() {
            let w = ae / h;
            if e < n {
                diag[e] += w;
            }
            if e > 0 {
                diag[e - 1] += w;
            }
            if e > 0 && e < n {
                off[e - 1] -= w;
            }
        }
        let mut c = off.clone();
        let mut d = vec![h; n];
        if n > 1 {
            c[0] /= diag[0];
        }
        d[0] /= diag[0];
        for i in 1..n {
            let m = diag[i] - off[i - 1] * c[i - 1];
            if i < n - 1 {
                c[i] = off[i] / m;
            }
            d[i] = (d[i] - off[i - 1] * d[i - 1]) / m;
        }
        for i in (0..n - 1).rev() {
            let t = c[i] * d[i + 1];
            d[i] -= t;
        }
        d
    }

    /// Compliance `J = h Σ u`.
    #[must_use]
    pub fn compliance(&self, u: &[f64]) -> f64 {
        assert_eq!(
            u.len(),
            self.n,
            "state length must equal the number of interior nodes"
        );
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (self.n as f64 + 1.0);
        h * u.iter().sum::<f64>()
    }

    /// Element slope of the P1 solution.
    fn slope(&self, u: &[f64], e: usize) -> f64 {
        let n = self.n;
        #[allow(clippy::cast_precision_loss)]
        let h = 1.0 / (n as f64 + 1.0);
        let lo = if e == 0 { 0.0 } else { u[e - 1] };
        let hi = if e == n { 0.0 } else { u[e] };
        (hi - lo) / h
    }
}

/// ADJOINT attribution of a conductivity edit `a0 → a1` over named
/// channel masks (element index sets). Uses the EXACT bilinear
/// identity `J(a1) − J(a0) = −∫ Δa · u0′ · u1′` (compliance is
/// self-adjoint; both states enter, no linearization error), so the
/// channel terms sum to the observed change EXACTLY up to rounding —
/// the certified bound per node is a rounding allowance.
#[must_use]
pub fn adjoint_attribution(
    fixture: &Elliptic1d,
    a0: &[f64],
    a1: &[f64],
    channels: &[(&str, Vec<usize>)],
) -> Vec<ExplanationNode> {
    assert_eq!(a0.len(), a1.len(), "conductivity edits must align");
    let u0 = fixture.solve(a0);
    let u1 = fixture.solve(a1);
    #[allow(clippy::cast_precision_loss)]
    let h = 1.0 / (fixture.n as f64 + 1.0);
    channels
        .iter()
        .map(|(name, elems)| {
            let mut acc = 0.0f64;
            for &e in elems {
                assert!(
                    e <= fixture.n,
                    "channel element index {e} exceeds the fixture element count"
                );
                let da = a1[e] - a0[e];
                acc -= da * fixture.slope(&u0, e) * fixture.slope(&u1, e) * h;
            }
            // The Verified interval must use the SAME (relative) bound as the
            // node's `bound` field — a fixed absolute ±1e-12 half-width is NOT a
            // rigorous enclosure once |acc| is large (the fp summation error of
            // the discrete bilinear form is ~n·eps·|acc| ≫ 1e-12), so it would
            // EXCLUDE the exact discrete value it claims to enclose (bead 9sf6
            // F2). The relative bound covers the fp error for reasonable n.
            let bound = 1e-12 * acc.abs().max(1.0);
            ExplanationNode::new(
                name,
                acc,
                bound,
                Color::Verified {
                    lo: acc - bound,
                    hi: acc + bound,
                },
                vec![
                    "solve(a0)".to_string(),
                    "solve(a1)".to_string(),
                    format!("mask:{name}"),
                ],
            )
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Engine 2: PROVENANCE attribution (which edit moved the number).
// ---------------------------------------------------------------------------

/// Telescoping edit attribution: for a replayed ledger sequence of
/// states with QoI values, each edit's contribution is the exact
/// difference it produced. Exact by construction; re-derivable by
/// replaying the same sequence.
#[must_use]
pub fn provenance_attribution(edits: &[(String, f64, f64)]) -> Vec<ExplanationNode> {
    edits
        .iter()
        .map(|(name, before, after)| {
            ExplanationNode::new(
                &format!("edit:{name}"),
                after - before,
                0.0,
                Color::Verified {
                    lo: after - before,
                    hi: after - before,
                },
                vec![format!("replay:{name}")],
            )
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Engine 3: PHYSICAL decomposition — the far-field drag flagship.
// ---------------------------------------------------------------------------

/// Lifting-line wing fixture: span-stations with circulation Γ(y) on a
/// span `b`, freestream `v_inf`, reference area `s_ref`.
#[derive(Debug, Clone)]
pub struct LiftingLine {
    /// Station midpoint circulations.
    pub gamma: Vec<f64>,
    /// Span.
    pub b: f64,
    /// Freestream speed.
    pub v_inf: f64,
    /// Reference area.
    pub s_ref: f64,
}

impl LiftingLine {
    /// Elliptic distribution `Γ = Γ0 √(1 − (2y/b)²)` at `n` stations.
    #[must_use]
    pub fn elliptic(gamma0: f64, b: f64, v_inf: f64, s_ref: f64, n: usize) -> LiftingLine {
        assert!(n > 0, "lifting-line station count must be positive");
        assert!(
            gamma0.is_finite()
                && b.is_finite()
                && b > 0.0
                && v_inf.is_finite()
                && v_inf > 0.0
                && s_ref.is_finite()
                && s_ref > 0.0,
            "lifting-line parameters must be finite with positive span, speed, and reference area"
        );
        let gamma = (0..n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let y = -0.5 + (i as f64 + 0.5) / n as f64; // 2y/b in (−1,1)
                gamma0 * (1.0 - (2.0 * y) * (2.0 * y)).max(0.0).sqrt()
            })
            .collect();
        LiftingLine {
            gamma,
            b,
            v_inf,
            s_ref,
        }
    }

    fn assert_valid(&self) {
        assert!(
            !self.gamma.is_empty(),
            "lifting-line circulation stations must not be empty"
        );
        assert!(
            self.gamma.iter().all(|g| g.is_finite())
                && self.b.is_finite()
                && self.b > 0.0
                && self.v_inf.is_finite()
                && self.v_inf > 0.0
                && self.s_ref.is_finite()
                && self.s_ref > 0.0,
            "lifting-line state must be finite with positive span, speed, and reference area"
        );
    }

    /// Lift coefficient from the circulation integral (KJ theorem).
    #[must_use]
    pub fn cl(&self) -> f64 {
        self.assert_valid();
        #[allow(clippy::cast_precision_loss)]
        let dy = self.b / self.gamma.len() as f64;
        let lift_per_rho = self.v_inf * self.gamma.iter().sum::<f64>() * dy;
        lift_per_rho / (0.5 * self.v_inf * self.v_inf * self.s_ref)
    }

    /// INDUCED drag by the TREFFTZ-PLANE wake integral: the shed
    /// vorticity sheet's kinetic energy,
    /// `D_i/ρ = (1/4π) ΣΣ γ_i γ_j ln|y_i − y_j|`-free discrete form via
    /// downwash: `w(y_i) = Σ_j γ'_j / (4π (y_i − y_j))`,
    /// `D_i/ρ = Σ_i Γ_i w_i dy`. Deterministic midpoint discretization.
    #[must_use]
    pub fn induced_drag_coefficient(&self) -> f64 {
        self.assert_valid();
        let n = self.gamma.len();
        #[allow(clippy::cast_precision_loss)]
        let dy = self.b / n as f64;
        // Shed vorticity between stations: γ_shed = −dΓ/dy at panel
        // edges (n+1 trailing vortices including tips).
        // Downwash convention (Katz & Plotkin): w(y) =
        // −(1/4π)∫(dΓ/dy′)/(y−y′) dy′; the discrete jump ΔΓ_j enters
        // NEGATED, i.e. as (Γ_j − Γ_{j−1}) with the sign folded here —
        // the original left−right form double-negated and produced a
        // wake integral of the right magnitude and wrong sign (caught
        // by the analytic envelope in conformance).
        let mut shed = Vec::with_capacity(n + 1);
        for j in 0..=n {
            let left = if j == 0 { 0.0 } else { self.gamma[j - 1] };
            let right = if j == n { 0.0 } else { self.gamma[j] };
            shed.push(right - left);
        }
        let mut drag_per_rho = 0.0f64;
        for (i, &g) in self.gamma.iter().enumerate() {
            #[allow(clippy::cast_precision_loss)]
            let yi = (i as f64 + 0.5) * dy;
            // Downwash at station i from every trailing vortex edge.
            let mut w = 0.0f64;
            for (j, &sv) in shed.iter().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let yj = j as f64 * dy;
                let r = yi - yj;
                w += sv / (4.0 * std::f64::consts::PI * r);
            }
            drag_per_rho += g * w * dy;
        }
        drag_per_rho / (0.5 * self.v_inf * self.v_inf * self.s_ref)
    }

    /// Aspect ratio.
    #[must_use]
    pub fn aspect_ratio(&self) -> f64 {
        self.assert_valid();
        self.b * self.b / self.s_ref
    }
}

/// The FLAGSHIP: decompose total drag into (induced, viscous, wave)
/// with bounds, via the wake integral + a strip-friction model + the
/// declared-zero subsonic wave channel. `cd_total_observed` is the
/// near-field measurement being explained.
#[must_use]
pub fn drag_decomposition(
    wing: &LiftingLine,
    cf_strip: f64,
    wetted_over_sref: f64,
    cd_total_observed: f64,
    threshold: f64,
) -> Explanation {
    assert!(
        cf_strip.is_finite()
            && cf_strip >= 0.0
            && wetted_over_sref.is_finite()
            && wetted_over_sref >= 0.0,
        "drag decomposition inputs must be finite and non-negative"
    );
    let n_stations = wing.gamma.len();
    let cdi = wing.induced_drag_coefficient();
    // Discretization bound for the wake integral: the midpoint panel
    // scheme converges O(1/n); the certified allowance is the analytic
    // envelope C/n with C = CDi itself (conservative on the elliptic
    // family, verified in conformance).
    #[allow(clippy::cast_precision_loss)]
    let cdi_bound = cdi.abs() / n_stations as f64 * 4.0;
    let cdv = cf_strip * wetted_over_sref;
    let nodes = vec![
        ExplanationNode::new(
            "induced (Trefftz wake integral)",
            cdi,
            cdi_bound,
            Color::Verified {
                lo: cdi - cdi_bound,
                hi: cdi + cdi_bound,
            },
            vec![
                "wake-integral".to_string(),
                format!("stations:{n_stations}"),
            ],
        ),
        ExplanationNode::new(
            "viscous (strip friction)",
            cdv,
            0.15 * cdv,
            Color::Estimated {
                estimator: "strip-friction".to_string(),
                dispersion: 0.15 * cdv,
            },
            vec![
                format!("cf:{cf_strip}"),
                format!("swet/sref:{wetted_over_sref}"),
            ],
        ),
        ExplanationNode::new(
            "wave (declared zero: subsonic regime)",
            0.0,
            0.0,
            Color::Verified { lo: 0.0, hi: 0.0 },
            vec!["regime:subsonic".to_string()],
        ),
    ];
    finalize(nodes, cd_total_observed, threshold)
}
