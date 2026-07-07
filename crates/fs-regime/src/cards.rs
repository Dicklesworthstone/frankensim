//! The solver/model registry and the ADMISSION gate: every FLUX/solid
//! model registers its validity domain in dimensionless-group space
//! (fs-evidence `ModelCard`s), and `admit` turns "is this solver even
//! valid here?" into a checkable, ledgerable verdict with RANKED
//! alternatives instead of a silent wrong answer.

use crate::RegimeError;
use fs_evidence::{Ambition, ModelCard, ProvenanceHash, ValidityDomain};
use fs_math::det;
use std::collections::BTreeMap;

fn log10(x: f64) -> f64 {
    det::ln(x) / core::f64::consts::LN_10
}

/// The built-in registry of FLUX/solid model cards, validity boxes in
/// dimensionless-group space.
#[must_use]
pub fn flux_model_cards() -> Vec<ModelCard> {
    let card = |name: &str,
                ambition: Ambition,
                assumptions: Vec<&str>,
                validity: ValidityDomain,
                failures: Vec<&str>,
                disc: f64| {
        ModelCard::new(
            name,
            "0.1.0",
            ambition,
            assumptions.into_iter().map(str::to_string).collect(),
            validity,
            failures.into_iter().map(str::to_string).collect(),
            disc,
        )
    };
    vec![
        card(
            "flux.stokes-creeping",
            Ambition::Solid,
            vec!["inertialess (Re << 1)", "incompressible", "Newtonian"],
            ValidityDomain::unconstrained().with("Re", 0.0, 1.0),
            vec!["silently wrong wakes/separation once inertia matters"],
            0.02,
        ),
        card(
            "flux.laminar-ns",
            Ambition::Solid,
            vec!["laminar (no closure)", "incompressible", "Newtonian"],
            ValidityDomain::unconstrained()
                .with("Re", 0.0, 2300.0)
                .with("Ma", 0.0, 0.3),
            vec!["under-resolved transition; no turbulence closure"],
            0.03,
        ),
        card(
            "flux.les-ns",
            Ambition::Frontier,
            vec![
                "LES closure (labeled)",
                "incompressible",
                "resolved large scales",
            ],
            ValidityDomain::unconstrained()
                .with("Re", 1.0e3, 1.0e9)
                .with("Ma", 0.0, 0.3),
            vec!["near-wall modeling error", "closure discrepancy"],
            0.10,
        ),
        card(
            "flux.free-surface-lbm",
            Ambition::Solid,
            vec![
                "weakly compressible lattice (Ma_lattice <= 0.3)",
                "free-surface capture",
            ],
            ValidityDomain::unconstrained()
                .with("Re", 1.0, 1.0e5)
                .with("Ma", 0.0, 0.3),
            vec![
                "lattice Mach artifacts at high speed",
                "thin-film breakup under-resolution",
            ],
            0.05,
        ),
        card(
            "flux.potential-flow",
            Ambition::Solid,
            vec!["inviscid", "irrotational", "attached flow"],
            ValidityDomain::unconstrained().with("Re", 1.0e4, f64::INFINITY),
            vec!["no separation, no drag from viscosity (d'Alembert)"],
            0.15,
        ),
        card(
            "solid.euler-bernoulli",
            Ambition::Solid,
            vec![
                "plane sections remain plane",
                "shear deformation negligible",
            ],
            ValidityDomain::unconstrained().with("slenderness", 20.0, f64::INFINITY),
            vec!["stiff-in-error for deep beams (shear ignored)"],
            0.02,
        ),
        card(
            "solid.timoshenko",
            Ambition::Solid,
            vec!["first-order shear deformation"],
            ValidityDomain::unconstrained().with("slenderness", 5.0, f64::INFINITY),
            vec!["cross-section warping ignored"],
            0.02,
        ),
    ]
}

/// The admission verdict for one model at one group point.
#[derive(Debug, Clone, PartialEq)]
pub struct Admission {
    /// The queried model.
    pub model: String,
    /// Whether the model's validity domain contains the group point.
    pub allowed: bool,
    /// Human-readable reasons when refused (one per violated bound).
    pub reasons: Vec<String>,
    /// Alternatives ranked by log-space distance to their validity boxes
    /// (admissible ones first, distance 0).
    pub alternatives: Vec<(String, f64)>,
    /// Provenance of the verdict (hash of model + point).
    pub provenance: ProvenanceHash,
}

/// Log-space distance from a group point to a card's validity box: 0
/// inside; the sum over violated bounds of `|log10(value/bound)|`;
/// a unit penalty per constrained-but-unavailable group.
#[must_use]
pub fn distance_to_validity(card: &ModelCard, groups: &BTreeMap<String, f64>) -> f64 {
    let mut d = 0.0f64;
    for param in card.validity.param_names() {
        let Some((lo, hi)) = card.validity.bound(&param) else {
            continue;
        };
        match groups.get(&param) {
            None => d += 1.0,
            Some(&v) => {
                if v < lo && lo > 0.0 && v > 0.0 {
                    d += log10(lo / v).abs();
                } else if v < lo {
                    d += 1.0;
                }
                if v > hi && hi.is_finite() && hi > 0.0 && v > 0.0 {
                    d += log10(v / hi).abs();
                }
            }
        }
    }
    d
}

/// Gate a model choice against the computed groups.
///
/// # Errors
/// [`RegimeError::UnknownModel`] when the name is not in the registry.
pub fn admit(
    registry: &[ModelCard],
    groups: &BTreeMap<String, f64>,
    model: &str,
) -> Result<Admission, RegimeError> {
    let card =
        registry
            .iter()
            .find(|c| c.name == model)
            .ok_or_else(|| RegimeError::UnknownModel {
                name: model.to_string(),
            })?;
    let allowed = card.validity.contains(groups);
    let mut reasons = Vec::new();
    if !allowed {
        for param in card.validity.param_names() {
            let (lo, hi) = card.validity.bound(&param).unwrap_or((0.0, f64::INFINITY));
            match groups.get(&param) {
                None => reasons.push(format!(
                    "{model} constrains {param} in [{lo}, {hi}] but the regime report has \
                     no {param} (missing role inputs)"
                )),
                Some(&v) if v < lo || v > hi => reasons.push(format!(
                    "{param} = {v:.4e} outside {model}'s validity [{lo}, {hi}]"
                )),
                Some(_) => {}
            }
        }
    }
    let mut alternatives: Vec<(String, f64)> = registry
        .iter()
        .filter(|c| c.name != model)
        .map(|c| (c.name.clone(), distance_to_validity(c, groups)))
        .collect();
    alternatives.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
    let mut provenance_text = format!("admit:{model}");
    for (k, v) in groups {
        let _ = core::fmt::Write::write_fmt(&mut provenance_text, format_args!(";{k}={v}"));
    }
    Ok(Admission {
        model: model.to_string(),
        allowed,
        reasons,
        alternatives,
        provenance: ProvenanceHash::of_bytes(provenance_text.as_bytes()),
    })
}
