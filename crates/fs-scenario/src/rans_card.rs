//! Low-Re forced-convection RANS model card (bead
//! frankensim-extreal-program-f85xj.5.8.1).
//!
//! This module freezes, BEFORE any implementation claim (.5.8.2), every
//! governing term, closure coefficient, wall treatment, scalar model,
//! porous/buoyancy option, boundary condition, unit, regime, determinism
//! and cancellation semantic, budget, validation-plan family, falsifier,
//! and no-claim surface of the low-Re RANS fidelity-graph node (e10).
//!
//! Design contract:
//! - Cards are built through [`RansCardDraft`] and admitted once by
//!   [`RansCardDraft::freeze`]; the resulting [`RansModelCard`] is opaque,
//!   so a frozen card cannot be mutated in place.
//! - [`RansCardDraft::freeze`] refuses malformed drafts, out-of-bound
//!   coefficients, missing governing terms, an empty or transition-claiming
//!   no-claim surface, oversized manifests, and unavailable capabilities.
//! - The canonical statement manifest is content-hashed so the independent
//!   adjudicator (.5.8.4) binds against exact bytes.
//!
//! Citations recorded on the card are bibliographic strings; they justify
//! the closure CHOICE, not any measured agreement (that is .5.8.3's lane).

use std::collections::BTreeMap;

/// Frozen schema identity of the card family.
pub const RANS_MODEL_CARD_SCHEMA: &str = "frankensim.rans-model-card.v1";

/// Canonical maximum serialized manifest size at freeze time.
pub const MAX_MANIFEST_BYTES: usize = 32 * 1024;

/// Governing mean-flow/scalar terms that MUST be present on a card.
pub const REQUIRED_TERMS: [&str; 8] = [
    "continuity-incompressible",
    "momentum-mean-x",
    "momentum-mean-y",
    "momentum-mean-z",
    "turbulent-kinetic-energy",
    "dissipation-rate",
    "scalar-temperature",
    "boussinesq-source-optional",
];

/// Closure coefficient set of Launder–Sharma (1974) low-Re k–epsilon,
/// with admissible bounds enforced at freeze time.
#[derive(Clone, Copy, Debug)]
pub struct LaunderSharmaCoefficients {
    /// C_mu.
    pub c_mu: f64,
    /// C_eps_1.
    pub c_eps_1: f64,
    /// C_eps_2.
    pub c_eps_2: f64,
    /// sigma_k.
    pub sigma_k: f64,
    /// sigma_eps.
    pub sigma_eps: f64,
}

impl LaunderSharmaCoefficients {
    /// Canonical LS74 values.
    #[must_use]
    pub const fn launder_sharma_1974() -> Self {
        Self {
            c_mu: 0.09,
            c_eps_1: 1.44,
            c_eps_2: 1.92,
            sigma_k: 1.0,
            sigma_eps: 1.3,
        }
    }

    fn bounds_ok(&self) -> bool {
        (0.08..=0.10).contains(&self.c_mu)
            && (1.30..=1.55).contains(&self.c_eps_1)
            && (1.80..=2.00).contains(&self.c_eps_2)
            && (0.9..=1.3).contains(&self.sigma_k)
            && (1.1..=1.5).contains(&self.sigma_eps)
    }
}

/// Wall treatment declaration. The frozen choice resolves the mean flow to
/// the viscous sublayer (low-Re integration); wall functions are refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WallTreatment {
    /// Integrate to the wall with damping functions; first cell y+ < 1.
    ResolveToViscousSublayer,
}

/// Buoyancy option: Boussinesq source, feature-gated inside the solver and
/// bounded by a physically plausible thermal expansion coefficient.
#[derive(Clone, Copy, Debug)]
pub struct BoussinesqOption {
    /// Enabled flag (solver-side gate must match).
    pub enabled: bool,
    /// Thermal expansion coefficient beta [1/K] bound-checked when enabled.
    pub beta_per_k: Option<f64>,
    /// Reference temperature [K].
    pub reference_temperature_k: f64,
}

/// Porous-media momentum sink for fin arrays (Darcy–Forchheimer form),
/// carried by the card and gated off until stage (d).
#[derive(Clone, Copy, Debug)]
pub struct PorousFinSink {
    /// Enabled flag.
    pub enabled: bool,
    /// Permeability K [m^2].
    pub permeability_m2: Option<f64>,
    /// Forchheimer coefficient c_F [-].
    pub forchheimer_c_f: Option<f64>,
}

/// Typed admission failure with the violated clause.
#[derive(Clone, Debug, PartialEq)]
pub enum AdmissionError {
    /// A required governing term is missing from the term list.
    MissingTerm {
        /// Missing term id.
        term: &'static str,
    },
    /// A closure coefficient is outside its admissible bound.
    CoefficientOutOfBounds {
        /// Coefficient name.
        name: &'static str,
    },
    /// The no-claim surface is empty or claims something forbidden.
    NoClaimViolation {
        /// Human-readable cause.
        what: String,
    },
    /// A numeric field was non-finite.
    NonFinite {
        /// Field name.
        field: &'static str,
    },
    /// Regime bounds are inverted or non-physical.
    InvalidRegime {
        /// Human-readable cause.
        what: String,
    },
    /// Canonical manifest exceeds the size cap.
    ManifestTooLarge {
        /// Actual byte count.
        actual: usize,
    },
    /// A required capability is unavailable.
    CapabilityUnavailable {
        /// Capability name.
        capability: &'static str,
    },
}

/// Builder draft; every field public, admitted only at [`Self::freeze`].
#[derive(Clone, Debug)]
pub struct RansCardDraft {
    /// Instrumented system family this card applies to.
    pub system_family: String,
    /// Required governing terms (subset superset of [`REQUIRED_TERMS`]).
    pub governing_terms: Vec<String>,
    /// Closure coefficients.
    pub coefficients: LaunderSharmaCoefficients,
    /// Damping-function formulas, verbatim (frozen text authority).
    pub damping_formulas: BTreeMap<String, String>,
    /// Wall treatment.
    pub wall_treatment: WallTreatment,
    /// Turbulent Prandtl number (constant model).
    pub turbulent_prandtl: f64,
    /// Boussinesq buoyancy option.
    pub boussinesq: BoussinesqOption,
    /// Porous fin sink option.
    pub porous_fin: PorousFinSink,
    /// Reynolds-number applicability band (D_h based).
    pub reynolds_band: (f64, f64),
    /// Boundary condition inventory (named BCs with units strings).
    pub boundary_conditions: BTreeMap<String, String>,
    /// Discretization targets (named quantity -> target string).
    pub discretization_targets: BTreeMap<String, String>,
    /// Solver iteration cap.
    pub max_iterations: u32,
    /// Relative residual tolerance at which the solve may stop.
    pub residual_tolerance_rel: f64,
    /// Validation-case families (exact envelope numbers defer to .5.8.3).
    pub validation_case_families: Vec<String>,
    /// Executable falsifier descriptions.
    pub falsifiers: Vec<String>,
    /// Explicit exclusion statements; MUST include the transition refusal.
    pub exclusions: Vec<String>,
    /// Feature gate name in the owning crate.
    pub feature_gate: String,
}

impl RansCardDraft {
    /// A draft pre-filled with the frozen LS74 choice and every mandatory
    /// clause, ready for caller review and [`Self::freeze`].
    #[must_use]
    pub fn launder_sharma_channel(system_family: impl Into<String>) -> Self {
        let mut damping = BTreeMap::new();
        damping.insert(
            "f_mu".to_string(),
            "exp(-3.4 / (1 + Re_t/50)^2)".to_string(),
        );
        damping.insert("f_2".to_string(), "1 - 0.3 exp(-(Re_t/6.5)^2)".to_string());
        let mut bcs = BTreeMap::new();
        bcs.insert("inlet".to_string(), "velocity [m/s], uniform".to_string());
        bcs.insert("walls".to_string(), "no-slip, resolved y+ < 1".to_string());
        bcs.insert(
            "thermal-wall".to_string(),
            "fixed heat flux [W/m^2] or fixed temperature [K]".to_string(),
        );
        bcs.insert(
            "outlet".to_string(),
            "zero-gradient pressure [Pa]".to_string(),
        );
        let mut targets = BTreeMap::new();
        targets.insert(
            "wall-resolution".to_string(),
            "first-cell y+ <= 1 on all resolved walls".to_string(),
        );
        targets.insert(
            "channel-modes".to_string(),
            ">= 40 cells per channel half-width at Re mid-band".to_string(),
        );
        Self {
            system_family: system_family.into(),
            governing_terms: REQUIRED_TERMS.iter().map(|t| (*t).to_string()).collect(),
            coefficients: LaunderSharmaCoefficients::launder_sharma_1974(),
            damping_formulas: damping,
            wall_treatment: WallTreatment::ResolveToViscousSublayer,
            turbulent_prandtl: 0.85,
            boussinesq: BoussinesqOption {
                enabled: false,
                beta_per_k: None,
                reference_temperature_k: 300.0,
            },
            porous_fin: PorousFinSink {
                enabled: false,
                permeability_m2: None,
                forchheimer_c_f: None,
            },
            reynolds_band: (500.0, 20_000.0),
            boundary_conditions: bcs,
            discretization_targets: targets,
            max_iterations: 20_000,
            residual_tolerance_rel: 1.0e-8,
            validation_case_families: vec![
                "vvreg:thermal-level-a:laminar-channel-friction".to_string(),
                "vvreg:thermal-level-b:cross-code-cavity".to_string(),
                "vvreg:thermal-level-b:cross-code-heated-channel".to_string(),
            ],
            falsifiers: vec![
                "laminar-channel friction factor within 5% of 64/Re_Dh at Re<=2000 equivalent"
                    .to_string(),
                "log-region mean-profile slope kappa within 0.38..0.44 where y+ in 30..300"
                    .to_string(),
                "decaying-HIT k(t) monotone decay under periodic reset seeds".to_string(),
            ],
            exclusions: vec![
                "NO transitional-flow validity between laminar onset and fully turbulent;"
                    .to_string(),
                "NO unvalidated turbulence-model authority; discrepancy vs correlations is "
                    .to_string()
                    + "fidelity-graph edge data, never an upgrade.",
                "NO compressibility, shocks, or high-Mach effects.".to_string(),
            ],
            feature_gate: "rans-rung".to_string(),
        }
    }

    /// Admit the draft into an immutable card.
    ///
    /// # Errors
    /// Every clause documented on [`AdmissionError`], plus the transition
    /// no-claim requirement and the manifest size cap.
    pub fn freeze(self) -> Result<RansModelCard, AdmissionError> {
        // Terms: every required term present.
        for term in REQUIRED_TERMS {
            if !self.governing_terms.iter().any(|g| g == term) {
                return Err(AdmissionError::MissingTerm { term });
            }
        }
        // Coefficients finite and bounded.
        let c = self.coefficients;
        for (name, value) in [
            ("c_mu", c.c_mu),
            ("c_eps_1", c.c_eps_1),
            ("c_eps_2", c.c_eps_2),
            ("sigma_k", c.sigma_k),
            ("sigma_eps", c.sigma_eps),
        ] {
            if !value.is_finite() {
                return Err(AdmissionError::NonFinite { field: name });
            }
        }
        if !c.bounds_ok() {
            return Err(AdmissionError::CoefficientOutOfBounds { name: "closure" });
        }
        if !self.turbulent_prandtl.is_finite() || !(0.5..=1.2).contains(&self.turbulent_prandtl) {
            return Err(AdmissionError::CoefficientOutOfBounds {
                name: "turbulent_prandtl",
            });
        }
        // Boussinesq bounds when enabled.
        if self.boussinesq.enabled {
            match self.boussinesq.beta_per_k {
                Some(beta) if (1.0e-4..=1.0e-2).contains(&beta) && beta.is_finite() => {}
                _ => return Err(AdmissionError::CoefficientOutOfBounds { name: "beta_per_k" }),
            }
            if !(1.0..=2_000.0).contains(&self.boussinesq.reference_temperature_k)
                || !self.boussinesq.reference_temperature_k.is_finite()
            {
                return Err(AdmissionError::InvalidRegime {
                    what: "reference temperature out of physical range".to_string(),
                });
            }
        }
        // Porous bounds when enabled.
        if self.porous_fin.enabled {
            match (
                self.porous_fin.permeability_m2,
                self.porous_fin.forchheimer_c_f,
            ) {
                (Some(k), Some(cf))
                    if (1.0e-12..=1.0e-6).contains(&k) && (0.0..=10.0).contains(&cf) => {}
                _ => return Err(AdmissionError::CoefficientOutOfBounds { name: "porous_fin" }),
            }
        }
        // Regime sanity.
        let (re_lo, re_hi) = self.reynolds_band;
        if !(re_lo.is_finite() && re_hi.is_finite() && re_hi > re_lo && re_lo > 0.0) {
            return Err(AdmissionError::InvalidRegime {
                what: "reynolds band must be positive and increasing".to_string(),
            });
        }
        if !(self.residual_tolerance_rel > 0.0 && self.residual_tolerance_rel < 1.0) {
            return Err(AdmissionError::InvalidRegime {
                what: "residual tolerance must lie in (0,1)".to_string(),
            });
        }
        // No-claim law: exclusions non-empty AND transition refusal present
        // verbatim enough to be mechanically checked.
        if self.exclusions.is_empty() {
            return Err(AdmissionError::NoClaimViolation {
                what: "exclusions empty".to_string(),
            });
        }
        // Mechanical no-claim law: at least one exclusion must be phrased
        // as a refusal ("NO ...") naming transition. A card that merely
        // mentions transition while claiming validity refuses here.
        let ok = self.exclusions.iter().any(|e| {
            let l = e.to_ascii_lowercase();
            l.starts_with("no ") && l.contains("transition")
        });
        if !ok {
            return Err(AdmissionError::NoClaimViolation {
                what: "the transition NO-claim is mandatory and cannot be erased or inverted"
                    .to_string(),
            });
        }
        if self.falsifiers.is_empty() {
            return Err(AdmissionError::NoClaimViolation {
                what: "a card without falsifiers cannot be frozen".to_string(),
            });
        }
        // Capabilities: freezing declares the solver-side gates it needs.
        // Solver-side capabilities that do not exist yet (.5.8.2 lands
        // them): freezing may DECLARE such options only while disabled.
        const AVAILABLE_CAPABILITIES: [&str; 3] =
            ["mean-flow", "scalar-temperature", "conjugate-wall"];
        if self.boussinesq.enabled && !AVAILABLE_CAPABILITIES.contains(&"buoyancy-source") {
            return Err(AdmissionError::CapabilityUnavailable {
                capability: "buoyancy-source",
            });
        }
        if self.porous_fin.enabled && !AVAILABLE_CAPABILITIES.contains(&"porous-sink") {
            return Err(AdmissionError::CapabilityUnavailable {
                capability: "porous-sink",
            });
        }

        let card = RansModelCard {
            schema: RANS_MODEL_CARD_SCHEMA.to_string(),
            system_family: self.system_family,
            governing_terms: self.governing_terms,
            coefficients: self.coefficients,
            damping_formulas: self.damping_formulas,
            wall_treatment: self.wall_treatment,
            turbulent_prandtl: self.turbulent_prandtl,
            boussinesq: self.boussinesq,
            porous_fin: self.porous_fin,
            reynolds_band: self.reynolds_band,
            boundary_conditions: self.boundary_conditions,
            discretization_targets: self.discretization_targets,
            max_iterations: self.max_iterations,
            residual_tolerance_rel: self.residual_tolerance_rel,
            validation_case_families: self.validation_case_families,
            falsifiers: self.falsifiers,
            exclusions: self.exclusions,
            feature_gate: self.feature_gate,
            citations: vec![
                "Launder & Spalding (1974), The numerical computation of turbulent flows, "
                    .to_string()
                    + "Comput. Methods Appl. Mech. Eng. 3:269-289",
                "Kays (1994), Turbulent Prandtl number - where are we?, J. Heat Transfer "
                    .to_string()
                    + "116:284-295",
                "Boussinesq (1903), Theorie de l'ecoulement tourbillonnant, Paris".to_string(),
            ],
        };
        // Size cap on the canonical manifest bytes.
        let manifest_len = card
            .statement_manifest()
            .iter()
            .map(|(k, v)| k.len() + v.len())
            .sum::<usize>();
        if manifest_len > MAX_MANIFEST_BYTES {
            return Err(AdmissionError::ManifestTooLarge {
                actual: manifest_len,
            });
        }
        Ok(card)
    }
}

/// An admitted, immutable model card.
#[derive(Clone, Debug)]
pub struct RansModelCard {
    schema: String,
    system_family: String,
    governing_terms: Vec<String>,
    coefficients: LaunderSharmaCoefficients,
    damping_formulas: BTreeMap<String, String>,
    wall_treatment: WallTreatment,
    turbulent_prandtl: f64,
    boussinesq: BoussinesqOption,
    porous_fin: PorousFinSink,
    reynolds_band: (f64, f64),
    boundary_conditions: BTreeMap<String, String>,
    discretization_targets: BTreeMap<String, String>,
    max_iterations: u32,
    residual_tolerance_rel: f64,
    validation_case_families: Vec<String>,
    falsifiers: Vec<String>,
    exclusions: Vec<String>,
    feature_gate: String,
    citations: Vec<String>,
}

impl RansModelCard {
    /// Schema identity.
    #[must_use]
    pub fn schema(&self) -> &str {
        &self.schema
    }

    /// Closure coefficients.
    #[must_use]
    pub const fn coefficients(&self) -> &LaunderSharmaCoefficients {
        &self.coefficients
    }

    /// Turbulent Prandtl number.
    #[must_use]
    pub const fn turbulent_prandtl(&self) -> f64 {
        self.turbulent_prandtl
    }

    /// Reynolds applicability band.
    #[must_use]
    pub const fn reynolds_band(&self) -> (f64, f64) {
        self.reynolds_band
    }

    /// Solver iteration cap.
    #[must_use]
    pub const fn max_iterations(&self) -> u32 {
        self.max_iterations
    }

    /// Relative residual stop tolerance.
    #[must_use]
    pub const fn residual_tolerance_rel(&self) -> f64 {
        self.residual_tolerance_rel
    }

    /// Exclusion statements (the no-claim surface).
    #[must_use]
    pub fn exclusions(&self) -> &[String] {
        &self.exclusions
    }

    /// Bibliographic citations justifying the closure choice.
    #[must_use]
    pub fn citations(&self) -> &[String] {
        &self.citations
    }

    /// Whether the Boussinesq source is admitted on this card.
    #[must_use]
    pub const fn boussinesq_enabled(&self) -> bool {
        self.boussinesq.enabled
    }

    /// Whether the porous fin sink is admitted on this card.
    #[must_use]
    pub const fn porous_fin_enabled(&self) -> bool {
        self.porous_fin.enabled
    }

    /// Canonical `(key, value)` statement pairs, sorted; the exact bytes of
    /// this listing are what [`Self::manifest_hash`] binds.
    #[must_use]
    pub fn statement_manifest(&self) -> Vec<(String, String)> {
        let mut m: Vec<(String, String)> = Vec::new();
        m.push(("schema".into(), self.schema.clone()));
        m.push(("system_family".into(), self.system_family.clone()));
        m.push(("governing_terms".into(), self.governing_terms.join(",")));
        m.push(("c_mu".into(), format!("{}", self.coefficients.c_mu)));
        m.push(("c_eps_1".into(), format!("{}", self.coefficients.c_eps_1)));
        m.push(("c_eps_2".into(), format!("{}", self.coefficients.c_eps_2)));
        m.push(("sigma_k".into(), format!("{}", self.coefficients.sigma_k)));
        m.push((
            "sigma_eps".into(),
            format!("{}", self.coefficients.sigma_eps),
        ));
        for (k, v) in &self.damping_formulas {
            m.push((format!("damping/{k}"), v.clone()));
        }
        m.push((
            "wall_treatment".into(),
            format!("{:?}", self.wall_treatment),
        ));
        m.push((
            "turbulent_prandtl".into(),
            format!("{}", self.turbulent_prandtl),
        ));
        m.push((
            "boussinesq".into(),
            format!(
                "enabled={} beta={:?} T_ref={}",
                self.boussinesq.enabled,
                self.boussinesq.beta_per_k,
                self.boussinesq.reference_temperature_k
            ),
        ));
        m.push((
            "porous_fin".into(),
            format!(
                "enabled={} K={:?} cF={:?}",
                self.porous_fin.enabled,
                self.porous_fin.permeability_m2,
                self.porous_fin.forchheimer_c_f
            ),
        ));
        m.push((
            "reynolds_band".into(),
            format!("{}..{}", self.reynolds_band.0, self.reynolds_band.1),
        ));
        for (k, v) in &self.boundary_conditions {
            m.push((format!("bc/{k}"), v.clone()));
        }
        for (k, v) in &self.discretization_targets {
            m.push((format!("target/{k}"), v.clone()));
        }
        m.push(("max_iterations".into(), format!("{}", self.max_iterations)));
        m.push((
            "residual_tolerance_rel".into(),
            format!("{}", self.residual_tolerance_rel),
        ));
        m.push((
            "validation_families".into(),
            self.validation_case_families.join(","),
        ));
        m.push(("falsifiers".into(), self.falsifiers.join(";")));
        m.push(("exclusions".into(), self.exclusions.join(";")));
        m.push(("feature_gate".into(), self.feature_gate.clone()));
        m.push(("citations".into(), self.citations.join(";")));
        m.sort();
        m
    }

    /// BLAKE3 digest over the canonical manifest; the adjudicator's binding.
    #[must_use]
    pub fn manifest_hash(&self) -> fs_blake3::ContentHash {
        let mut bytes = Vec::new();
        for (k, v) in self.statement_manifest() {
            bytes.extend_from_slice(k.as_bytes());
            bytes.push(0);
            bytes.extend_from_slice(v.as_bytes());
            bytes.push(0);
        }
        fs_blake3::hash_bytes(&bytes)
    }
}
