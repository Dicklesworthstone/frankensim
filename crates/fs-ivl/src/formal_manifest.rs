//! Formal statement, model, and proof-TCB freeze for selected interval primitives
//! (bead `frankensim-extreal-program-f85xj.3.8.1`).
//!
//! Freezes the exact formal specification for:
//! 1. `next_up` (directed successor)
//! 2. `next_down` (directed predecessor)
//! 3. Outward-rounded interval addition (`Interval::add`)
//! 4. Outward-rounded interval multiplication (`Interval::mul`)
//!
//! Includes toolchain, axioms, Trusted Computing Base (TCB), extraction boundary,
//! exceptional value policies, stretch lemmas, and residual no-claim surfaces.

/// Schema version for the formal proof manifest.
pub const FORMAL_MANIFEST_SCHEMA_VERSION: u32 = 1;

/// Domain separator for formal proof manifest identity.
pub const FORMAL_MANIFEST_DOMAIN: &str = "org.frankensim.fs-ivl.formal-manifest.v1";

/// A 64-bit deterministic content hash.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ManifestFingerprint(pub u64);

impl ManifestFingerprint {
    /// Convert to hex string.
    #[must_use]
    pub fn to_hex(self) -> String {
        format!("{:016x}", self.0)
    }
}

fn fnv1a_hash(bytes: &[u8]) -> ManifestFingerprint {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    ManifestFingerprint(h)
}

/// Classification of a theorem in the formal proof program.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TheoremClass {
    /// Minimum required core theorem (blocking 3.8.2/3.8.3).
    MinimumCore,
    /// Stretch lemma (non-blocking).
    Stretch,
}

/// Formal statement specification for one certified arithmetic primitive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormalTheoremSpec {
    /// Unique theorem identifier in the formal proof database.
    pub theorem_id: &'static str,
    /// Classification (MinimumCore vs Stretch).
    pub class: TheoremClass,
    /// Rust source symbol and crate path.
    pub rust_symbol: &'static str,
    /// Mathematical quantifier statement in formal logic.
    pub statement_logic: &'static str,
    /// IEEE-754 / extended real domain and exceptional cases handling.
    pub domain_and_exceptions: &'static str,
    /// Target proof vehicle module and lemma name.
    pub formal_lemma_target: &'static str,
}

/// Specification of the formal proof toolchain and Trusted Computing Base (TCB).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProofTcbSpec {
    /// Proof assistant vehicle (e.g. "Coq 8.18 / Flocq 4.1.0" or "Lean 4 / Mathlib").
    pub proof_vehicle: &'static str,
    /// Formal IEEE-754 floating-point model reference.
    pub ieee_model_ref: &'static str,
    /// Trusted axioms and base theories.
    pub trusted_axioms: &'static [&'static str],
    /// Extraction boundary and translation model.
    pub extraction_boundary: &'static str,
}

/// Residual no-claim boundaries explicitly disclaimed by this formal freeze.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResidualNoClaim {
    /// Disclaimed boundary name.
    pub boundary_name: &'static str,
    /// Plain language explanation of what is NOT certified.
    pub explanation: &'static str,
}

/// The immutable formal manifest for fs-ivl certified primitives.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FormalProofManifest<'a> {
    /// Version of the manifest schema.
    pub schema_version: u32,
    /// Trusted Computing Base specification.
    pub tcb: ProofTcbSpec,
    /// Minimum required core theorems for the proof program.
    pub minimum_theorems: &'a [FormalTheoremSpec],
    /// Optional stretch lemmas.
    pub stretch_theorems: &'a [FormalTheoremSpec],
    /// Explicit residual no-claim statements.
    pub no_claims: &'a [ResidualNoClaim],
}

/// The frozen proof TCB specification.
pub const FROZEN_PROOF_TCB: ProofTcbSpec = ProofTcbSpec {
    proof_vehicle: "Coq 8.18 / Flocq 4.1.0 (IEEE-754 binary64 formalization)",
    ieee_model_ref: "IEEE 754-2008 Standard for Floating-Point Arithmetic (binary64, roundNearestTiesToEven)",
    trusted_axioms: &[
        "Flocq Core IEEE-754 binary64 definition and rounding operator round(RNE, x)",
        "Correctly rounded IEEE-754 basic operations err <= 0.5 ULP in absence of overflow",
        "Bit-level sign-magnitude integer bijection on finite binary64 floats",
    ],
    extraction_boundary: "Pure computational model verified against bitwise MODEL v1; no runtime FPU mode alteration assumed",
};

/// The frozen minimum core theorems.
pub const FROZEN_MINIMUM_THEOREMS: [FormalTheoremSpec; 4] = [
    FormalTheoremSpec {
        theorem_id: "thm_next_up_enclosure",
        class: TheoremClass::MinimumCore,
        rust_symbol: "fs_math::next_up",
        statement_logic: "forall (x : binary64), x <> +inf -> ~isnan(x) -> Real(x) < Real(next_up(x))",
        domain_and_exceptions: "x in [-inf, +inf); next_up(+inf) = +inf; next_up(-0.0) = +min_subnormal; next_up(+0.0) = +min_subnormal; NaN preserved",
        formal_lemma_target: "fs_ivl_formal.primitives.thm_next_up_sound",
    },
    FormalTheoremSpec {
        theorem_id: "thm_next_down_enclosure",
        class: TheoremClass::MinimumCore,
        rust_symbol: "fs_math::next_down",
        statement_logic: "forall (x : binary64), x <> -inf -> ~isnan(x) -> Real(next_down(x)) < Real(x)",
        domain_and_exceptions: "x in (-inf, +inf]; next_down(-inf) = -inf; next_down(+0.0) = -min_subnormal; next_down(-0.0) = -min_subnormal; NaN preserved",
        formal_lemma_target: "fs_ivl_formal.primitives.thm_next_down_sound",
    },
    FormalTheoremSpec {
        theorem_id: "thm_interval_add_enclosure",
        class: TheoremClass::MinimumCore,
        rust_symbol: "fs_ivl::Interval::add",
        statement_logic: "forall (I1 I2 : Interval) (x y : Real), in_interval(x, I1) -> in_interval(y, I2) -> in_interval(x + y, I1 + I2)",
        domain_and_exceptions: "Finite and extended real endpoints; overflow bounded by [-inf, -f64::MAX] or [f64::MAX, +inf]; NaN rejected at construction",
        formal_lemma_target: "fs_ivl_formal.interval.thm_add_enclosure",
    },
    FormalTheoremSpec {
        theorem_id: "thm_interval_mul_enclosure",
        class: TheoremClass::MinimumCore,
        rust_symbol: "fs_ivl::Interval::mul",
        statement_logic: "forall (I1 I2 : Interval) (x y : Real), in_interval(x, I1) -> in_interval(y, I2) -> in_interval(x * y, I1 * I2)",
        domain_and_exceptions: "Finite and extended real endpoints; 0 * inf indeterminacies yield Interval::WHOLE [-inf, +inf]; finite overflow yields one-sided enclosure",
        formal_lemma_target: "fs_ivl_formal.interval.thm_mul_enclosure",
    },
];

/// The frozen stretch theorems.
pub const FROZEN_STRETCH_THEOREMS: [FormalTheoremSpec; 3] = [
    FormalTheoremSpec {
        theorem_id: "thm_interval_sub_enclosure",
        class: TheoremClass::Stretch,
        rust_symbol: "fs_ivl::Interval::sub",
        statement_logic: "forall (I1 I2 : Interval) (x y : Real), in_interval(x, I1) -> in_interval(y, I2) -> in_interval(x - y, I1 - I2)",
        domain_and_exceptions: "Defined as I1 + (-I2); negation is exact on binary64 floats (sign-bit flip)",
        formal_lemma_target: "fs_ivl_formal.interval.thm_sub_enclosure",
    },
    FormalTheoremSpec {
        theorem_id: "thm_interval_div_enclosure",
        class: TheoremClass::Stretch,
        rust_symbol: "fs_ivl::Interval::div",
        statement_logic: "forall (I1 I2 : Interval) (x y : Real), in_interval(x, I1) -> in_interval(y, I2) -> y <> 0 -> in_interval(x / y, I1 / I2)",
        domain_and_exceptions: "If 0 in I2, returns Interval::WHOLE [-inf, +inf]; 0/0 and inf/inf indeterminacies return Interval::WHOLE",
        formal_lemma_target: "fs_ivl_formal.interval.thm_div_enclosure",
    },
    FormalTheoremSpec {
        theorem_id: "thm_interval_sqrt_enclosure",
        class: TheoremClass::Stretch,
        rust_symbol: "fs_ivl::Interval::sqrt",
        statement_logic: "forall (I : Interval) (x : Real), in_interval(x, I) -> x >= 0 -> in_interval(sqrt(x), I.sqrt())",
        domain_and_exceptions: "Requires I.hi >= 0; negative lower bound clamped to 0; correctly rounded IEEE-754 sqrt with next_down/next_up nudge",
        formal_lemma_target: "fs_ivl_formal.interval.thm_sqrt_sound",
    },
];

/// The frozen residual no-claim boundary statements.
pub const FROZEN_RESIDUAL_NO_CLAIMS: [ResidualNoClaim; 4] = [
    ResidualNoClaim {
        boundary_name: "non_compliant_fpu_modes",
        explanation: "Proofs assume standard IEEE-754 binary64 operation without flush-to-zero (FTZ) or denormals-are-zero (DAZ) hardware flags active.",
    },
    ResidualNoClaim {
        boundary_name: "compiler_fast_math",
        explanation: "Proofs require strict IEEE-754 associativity; compiling with unsafe aggressive fast-math or associative-math flags breaks formal validity.",
    },
    ResidualNoClaim {
        boundary_name: "transcendental_functions",
        explanation: "Transcendental functions (exp, ln, sin, cos, tanh) rely on declared ULP budgets from fs-math rather than formal proof under this minimal TCB.",
    },
    ResidualNoClaim {
        boundary_name: "multivariate_taylor_models",
        explanation: "This formal manifest covers basic scalar and interval arithmetic primitives; multivariate Taylor models remain separate future scope.",
    },
];

/// The canonical frozen formal proof manifest.
pub const FROZEN_FORMAL_MANIFEST: FormalProofManifest<'static> = FormalProofManifest {
    schema_version: FORMAL_MANIFEST_SCHEMA_VERSION,
    tcb: FROZEN_PROOF_TCB,
    minimum_theorems: &FROZEN_MINIMUM_THEOREMS,
    stretch_theorems: &FROZEN_STRETCH_THEOREMS,
    no_claims: &FROZEN_RESIDUAL_NO_CLAIMS,
};

impl<'a> FormalProofManifest<'a> {
    /// Compute the immutable content hash of this formal manifest.
    #[must_use]
    pub fn content_hash(&self) -> ManifestFingerprint {
        let mut buf = Vec::new();
        buf.extend_from_slice(FORMAL_MANIFEST_DOMAIN.as_bytes());
        buf.extend_from_slice(&self.schema_version.to_le_bytes());
        buf.extend_from_slice(self.tcb.proof_vehicle.as_bytes());
        buf.extend_from_slice(self.tcb.ieee_model_ref.as_bytes());
        for ax in self.tcb.trusted_axioms {
            buf.extend_from_slice(ax.as_bytes());
        }
        buf.extend_from_slice(self.tcb.extraction_boundary.as_bytes());
        for thm in self.minimum_theorems {
            buf.extend_from_slice(thm.theorem_id.as_bytes());
            buf.extend_from_slice(thm.rust_symbol.as_bytes());
            buf.extend_from_slice(thm.statement_logic.as_bytes());
            buf.extend_from_slice(thm.domain_and_exceptions.as_bytes());
            buf.extend_from_slice(thm.formal_lemma_target.as_bytes());
        }
        for thm in self.stretch_theorems {
            buf.extend_from_slice(thm.theorem_id.as_bytes());
            buf.extend_from_slice(thm.rust_symbol.as_bytes());
            buf.extend_from_slice(thm.statement_logic.as_bytes());
            buf.extend_from_slice(thm.domain_and_exceptions.as_bytes());
            buf.extend_from_slice(thm.formal_lemma_target.as_bytes());
        }
        for nc in self.no_claims {
            buf.extend_from_slice(nc.boundary_name.as_bytes());
            buf.extend_from_slice(nc.explanation.as_bytes());
        }
        fnv1a_hash(&buf)
    }

    /// Validate manifest completeness against required minimum theorems.
    ///
    /// # Errors
    /// Returns an error string if any minimum theorem or required field is missing.
    pub fn validate(&self) -> Result<(), &'static str> {
        if self.schema_version != FORMAL_MANIFEST_SCHEMA_VERSION {
            return Err("schema version mismatch");
        }
        if self.tcb.proof_vehicle.is_empty() || self.tcb.ieee_model_ref.is_empty() {
            return Err("proof TCB missing required vehicle or model reference");
        }
        if self.tcb.trusted_axioms.is_empty() {
            return Err("proof TCB has no trusted axioms declared");
        }
        let required_ids = [
            "thm_next_up_enclosure",
            "thm_next_down_enclosure",
            "thm_interval_add_enclosure",
            "thm_interval_mul_enclosure",
        ];
        for req in required_ids {
            if !self.minimum_theorems.iter().any(|t| t.theorem_id == req) {
                return Err("missing required minimum core theorem");
            }
        }
        for thm in self.minimum_theorems {
            if thm.rust_symbol.is_empty()
                || thm.statement_logic.is_empty()
                || thm.domain_and_exceptions.is_empty()
                || thm.formal_lemma_target.is_empty()
            {
                return Err("theorem contains empty required field");
            }
        }
        if self.no_claims.is_empty() {
            return Err("residual no-claims list must not be empty");
        }
        Ok(())
    }
}
