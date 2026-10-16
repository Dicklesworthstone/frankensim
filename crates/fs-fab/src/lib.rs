//! fs-fab — manufacturing, fabrication & code compliance as its own layer.
//! Layer: L4.
//!
//! Optimization without fabrication semantics produces FANTASY ARTIFACTS: a
//! beautiful topology no one can build. This crate makes manufacturability a
//! first-class, CHECKABLE constraint layer — differentiable or certifiable
//! wherever possible — across additive, subtractive, casting, structural-steel,
//! and reinforced-concrete processes.
//!
//! Every [`FabConstraint`] declares its DIFFERENTIABILITY class (so the
//! optimizer knows whether it can take a gradient), its CERTIFICATE
//! availability, and its [`ConstraintKind`] (`Fabrication{process}` or
//! `Code{standard}`). Checking a design yields per-constraint detection +
//! localization + a REPAIR SUGGESTION for each violation ([`check_all`]). Cost
//! and embodied carbon are modeled quantities with uncertainty envelopes
//! ([`Estimate`], Evidence-typed) so "minimize embodied carbon subject to code
//! compliance" is a first-class study. Deterministic.

pub use fs_evidence::Color;

/// A fabrication process family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Process {
    /// Additive (powder-bed / deposition).
    Additive,
    /// Subtractive (CNC milling / turning).
    Subtractive,
    /// Casting / molding.
    Casting,
    /// Structural steel fabrication.
    StructuralSteel,
    /// Reinforced concrete.
    ReinforcedConcrete,
}

/// Whether a constraint is a process capability or a code/standard rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstraintKind {
    /// A fabrication capability of a process.
    Fabrication(Process),
    /// A code/standard rule (e.g. `"AISC-360"`, `"ACI-318"`).
    Code(&'static str),
}

/// The differentiability class of a constraint's margin (routes the optimizer).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Differentiability {
    /// Smooth in the design parameter (a gradient exists).
    Differentiable,
    /// Piecewise / non-smooth but has subgradients.
    Subdifferentiable,
    /// Discrete (catalog snapping, no gradient).
    Discrete,
}

/// Whether the constraint carries a certificate or is heuristic.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CertAvailability {
    /// A certified check.
    Certified,
    /// A heuristic check.
    Heuristic,
}

/// The direction of a one-sided constraint on a scalar feature value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sense {
    /// `value >= limit`.
    AtLeast,
    /// `value <= limit`.
    AtMost,
}

/// A manufacturing / code constraint on a scalar geometric feature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FabConstraint {
    /// A stable name.
    pub name: &'static str,
    /// Fabrication or code.
    pub kind: ConstraintKind,
    /// The differentiability class.
    pub differentiability: Differentiability,
    /// Certificate availability.
    pub cert: CertAvailability,
    /// One-sided direction.
    pub sense: Sense,
    /// The limit.
    pub limit: f64,
    /// The feature units (for reports).
    pub units: &'static str,
}

/// A repair suggestion for a violated constraint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Repair {
    /// The feature value that would just satisfy the constraint.
    pub target_value: f64,
    /// The change from the current value (`target − value`).
    pub delta: f64,
}

impl FabConstraint {
    /// The signed margin (`>= 0` ⇒ satisfied): distance into the feasible side.
    #[must_use]
    pub fn margin(&self, value: f64) -> f64 {
        match self.sense {
            Sense::AtLeast => value - self.limit,
            Sense::AtMost => self.limit - value,
        }
    }

    /// Is the feature value feasible?
    #[must_use]
    pub fn satisfied(&self, value: f64) -> bool {
        self.margin(value) >= 0.0
    }

    /// `d(margin)/d(value)` for a differentiable constraint (`±1`); `None` for a
    /// discrete one.
    #[must_use]
    pub fn margin_gradient(&self) -> Option<f64> {
        match self.differentiability {
            Differentiability::Discrete => None,
            _ => Some(match self.sense {
                Sense::AtLeast => 1.0,
                Sense::AtMost => -1.0,
            }),
        }
    }

    /// A repair suggestion if the constraint is violated: move the feature to
    /// the limit.
    #[must_use]
    pub fn repair(&self, value: f64) -> Option<Repair> {
        if self.satisfied(value) {
            None
        } else {
            Some(Repair {
                target_value: self.limit,
                delta: self.limit - value,
            })
        }
    }
}

// -- Constraint catalog -----------------------------------------------------

/// Additive: the self-supporting overhang angle from vertical must not exceed
/// `max_overhang_deg` (differentiable via surface-normal fields).
#[must_use]
pub fn overhang_angle(max_overhang_deg: f64) -> FabConstraint {
    FabConstraint {
        name: "additive-overhang-angle",
        kind: ConstraintKind::Fabrication(Process::Additive),
        differentiability: Differentiability::Differentiable,
        cert: CertAvailability::Heuristic,
        sense: Sense::AtMost,
        limit: max_overhang_deg,
        units: "deg",
    }
}

/// Additive: minimum printable feature size (medial-oracle bound).
#[must_use]
pub fn min_feature_size(min_mm: f64) -> FabConstraint {
    FabConstraint {
        name: "additive-min-feature",
        kind: ConstraintKind::Fabrication(Process::Additive),
        differentiability: Differentiability::Subdifferentiable,
        cert: CertAvailability::Certified,
        sense: Sense::AtLeast,
        limit: min_mm,
        units: "mm",
    }
}

/// Subtractive: a concave feature's radius must admit the tool (radius ≥ tool
/// radius) — CNC accessibility as a curvature bound.
#[must_use]
pub fn cnc_tool_radius(tool_radius_mm: f64) -> FabConstraint {
    FabConstraint {
        name: "cnc-tool-reachability",
        kind: ConstraintKind::Fabrication(Process::Subtractive),
        differentiability: Differentiability::Differentiable,
        cert: CertAvailability::Certified,
        sense: Sense::AtLeast,
        limit: tool_radius_mm,
        units: "mm",
    }
}

/// Casting: draft angle must be at least `min_draft_deg` for mold release.
#[must_use]
pub fn draft_angle(min_draft_deg: f64) -> FabConstraint {
    FabConstraint {
        name: "casting-draft-angle",
        kind: ConstraintKind::Fabrication(Process::Casting),
        differentiability: Differentiability::Differentiable,
        cert: CertAvailability::Certified,
        sense: Sense::AtLeast,
        limit: min_draft_deg,
        units: "deg",
    }
}

/// Structural steel (AISC-360 §J3.3): minimum bolt center-to-center spacing is
/// `2⅔ · bolt_diameter` (preferred `3d`).
#[must_use]
pub fn bolt_spacing_aisc(bolt_diameter_mm: f64) -> FabConstraint {
    FabConstraint {
        name: "bolt-spacing-min",
        kind: ConstraintKind::Code("AISC-360"),
        differentiability: Differentiability::Differentiable,
        cert: CertAvailability::Certified,
        sense: Sense::AtLeast,
        limit: (8.0 / 3.0) * bolt_diameter_mm,
        units: "mm",
    }
}

/// Assembly: a member's length must not exceed the transport limit — a discrete
/// catalog/standardization constraint.
#[must_use]
pub fn member_length_transport(max_length_m: f64) -> FabConstraint {
    FabConstraint {
        name: "member-transport-length",
        kind: ConstraintKind::Fabrication(Process::StructuralSteel),
        differentiability: Differentiability::Discrete,
        cert: CertAvailability::Certified,
        sense: Sense::AtMost,
        limit: max_length_m,
        units: "m",
    }
}

/// Reinforced concrete (ACI-318 §25.2.1): minimum clear bar spacing is
/// `max(25 mm, bar_diameter)`.
#[must_use]
pub fn rebar_spacing_aci(bar_diameter_mm: f64) -> FabConstraint {
    FabConstraint {
        name: "rebar-clear-spacing-min",
        kind: ConstraintKind::Code("ACI-318"),
        differentiability: Differentiability::Differentiable,
        cert: CertAvailability::Certified,
        sense: Sense::AtLeast,
        limit: 25.0_f64.max(bar_diameter_mm),
        units: "mm",
    }
}

// -- Evaluation + report ----------------------------------------------------

/// The result of checking one constraint against a feature value.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintResult {
    /// The constraint name.
    pub name: &'static str,
    /// The checked feature value.
    pub value: f64,
    /// The signed margin.
    pub margin: f64,
    /// Feasible?
    pub satisfied: bool,
    /// A repair suggestion if violated.
    pub repair: Option<Repair>,
}

/// Check one constraint against a value.
#[must_use]
pub fn evaluate(constraint: &FabConstraint, value: f64) -> ConstraintResult {
    ConstraintResult {
        name: constraint.name,
        value,
        margin: constraint.margin(value),
        satisfied: constraint.satisfied(value),
        repair: constraint.repair(value),
    }
}

/// A manufacturability report over many constraints.
#[derive(Debug, Clone, PartialEq)]
pub struct FabReport {
    /// Per-constraint results.
    pub results: Vec<ConstraintResult>,
    /// Is the whole design manufacturable (no fantasy artifact)?
    pub feasible: bool,
    /// The names of violated constraints (localization).
    pub violations: Vec<&'static str>,
}

/// Check a design: each `(constraint, feature value)` pair is evaluated;
/// violations are detected + localized + given a repair suggestion.
#[must_use]
pub fn check_all(pairs: &[(FabConstraint, f64)]) -> FabReport {
    let mut results = Vec::with_capacity(pairs.len());
    let mut violations = Vec::new();
    for (c, v) in pairs {
        let r = evaluate(c, *v);
        if !r.satisfied {
            violations.push(r.name);
        }
        results.push(r);
    }
    FabReport {
        feasible: violations.is_empty(),
        violations,
        results,
    }
}

// -- Cost & carbon (Evidence-typed) -----------------------------------------

/// A modeled quantity with an uncertainty envelope (Evidence-typed).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Estimate {
    /// The mean estimate.
    pub mean: f64,
    /// The relative standard deviation (`std / mean`).
    pub rel_std: f64,
}

impl Estimate {
    /// The absolute standard deviation.
    #[must_use]
    pub fn std(&self) -> f64 {
        self.mean.abs() * self.rel_std
    }

    /// A cost/carbon model is a MODELED quantity — estimated color.
    #[must_use]
    pub fn color(&self) -> Color {
        Color::Estimated {
            estimator: "fab-cost-model".to_string(),
            dispersion: self.std(),
        }
    }
}

/// Parametric process cost: `quantity · unit_rate`, with a relative uncertainty.
#[must_use]
pub fn process_cost(quantity: f64, unit_rate: f64, rel_uncertainty: f64) -> Estimate {
    Estimate {
        mean: quantity * unit_rate,
        rel_std: rel_uncertainty.max(0.0),
    }
}

/// Embodied carbon: `mass_kg · carbon_per_kg`, with a relative uncertainty.
#[must_use]
pub fn embodied_carbon(mass_kg: f64, carbon_per_kg: f64, rel_uncertainty: f64) -> Estimate {
    Estimate {
        mean: mass_kg * carbon_per_kg,
        rel_std: rel_uncertainty.max(0.0),
    }
}
