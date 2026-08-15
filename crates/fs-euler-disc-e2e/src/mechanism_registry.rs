//! Version-pinned source registry and rival mechanism cards for the
//! Euler-disc campaign (bead frankensim-euler-disc-emergent-flagship-t6314.1.2).
//!
//! Euler-disc dissipation is regime- and apparatus-dependent: published
//! models variously emphasize whole-gap viscosity, oscillatory air boundary
//! layers, rolling/contour resistance, slip, contact deformation,
//! vibration, and impacts — and several predict OVERLAPPING decay-exponent
//! families, so an exponent alone cannot identify a mechanism. This module
//! makes that epistemic state machine-checkable:
//!
//! - [`RegisteredSource`] rows are version-pinned and fail closed when the
//!   source type, version/date, locator, license, equations (for analytic
//!   models), assumptions, or transfer limitations are missing. A row's
//!   identity moves with every semantic field, so a version or assumption
//!   edit is a NEW source declaration.
//! - [`MechanismCard`]s state observables and the orthogonal interventions
//!   that could falsify or distinguish each rival; a card resolves only
//!   against registered rows, so a citation NAME alone cannot authorize a
//!   model.
//! - [`MechanismRegistry::mechanisms_matching_exponent`] returns EVERY card
//!   whose predicted exponent family covers an observation: the ambiguity
//!   is the result, not an error to hide.
//!
//! No-claims: registration is declaration bookkeeping, not validation. It
//! does not encode video rankings or an ~1 mm optimum as expected truth,
//! does not promote a preprint to consensus, and does not make any paper's
//! material/surface result transferable to a different finish, support,
//! gas, load, or scale — each row must SAY its transfer limitations.

use std::collections::BTreeMap;

use fs_blake3::ContentHash;

/// Versioned identity domain for source declarations.
pub const SOURCE_DECLARATION_DOMAIN: &str = "org.frankensim.euler-disc.registered-source.v1";
/// Versioned identity domain for mechanism cards.
pub const MECHANISM_CARD_DOMAIN: &str = "org.frankensim.euler-disc.mechanism-card.v1";
/// Bounded text field size.
pub const MAX_REGISTRY_TEXT_BYTES: usize = 4_096;
/// Bounded list size for any registry collection.
pub const MAX_REGISTRY_ITEMS: usize = 64;

/// Typed refusal of the registry boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegistryError {
    /// Stable machine rule.
    pub rule: &'static str,
    /// Field that refused.
    pub field: String,
    /// Human diagnosis.
    pub detail: String,
}

impl core::fmt::Display for RegistryError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}: {}: {}", self.rule, self.field, self.detail)
    }
}

impl std::error::Error for RegistryError {}

fn refuse(
    rule: &'static str,
    field: impl Into<String>,
    detail: impl Into<String>,
) -> RegistryError {
    RegistryError {
        rule,
        field: field.into(),
        detail: detail.into(),
    }
}

fn checked_text(field: &str, value: &str) -> Result<(), RegistryError> {
    if value.is_empty() || value.len() > MAX_REGISTRY_TEXT_BYTES {
        return Err(refuse(
            "euler-registry-text-bounds",
            field,
            format!("text must be 1..={MAX_REGISTRY_TEXT_BYTES} bytes"),
        ));
    }
    Ok(())
}

fn checked_list(field: &str, values: &[String], minimum: usize) -> Result<(), RegistryError> {
    if values.len() < minimum || values.len() > MAX_REGISTRY_ITEMS {
        return Err(refuse(
            "euler-registry-list-bounds",
            field,
            format!("must declare {minimum}..={MAX_REGISTRY_ITEMS} entries"),
        ));
    }
    for value in values {
        checked_text(field, value)?;
    }
    Ok(())
}

/// What kind of source a registry row pins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceType {
    /// Closed-form or asymptotic analytic model.
    AnalyticModel,
    /// Controlled measurement campaign.
    ExperimentalStudy,
    /// Preprint: explicitly NOT consensus; rows carry that status.
    Preprint,
    /// Video/observational material: motivating, never validating.
    Observational,
}

impl SourceType {
    const fn tag(self) -> u8 {
        match self {
            Self::AnalyticModel => 0,
            Self::ExperimentalStudy => 1,
            Self::Preprint => 2,
            Self::Observational => 3,
        }
    }
}

/// One version-pinned source declaration. Constructed only through
/// [`RegisteredSource::try_new`], which fails closed on every missing
/// semantic field.
#[derive(Debug, Clone, PartialEq)]
pub struct RegisteredSource {
    /// Stable registry id (e.g. `"moffatt-2000-nature"`).
    pub id: String,
    /// Source class.
    pub source_type: SourceType,
    /// Exact version or publication date being pinned.
    pub version_or_date: String,
    /// Resolvable locator (doi, arxiv id, or full bibliographic citation).
    pub locator: String,
    /// License or access status.
    pub license: String,
    /// Model equations (nonempty for analytic models; may be empty for
    /// purely observational rows).
    pub model_equations: Vec<String>,
    /// Stated assumptions (always nonempty: even observation assumes).
    pub assumptions: Vec<String>,
    /// Measured configurations, where the source reports them.
    pub measured_configurations: Vec<String>,
    /// QoIs the source reports or predicts.
    pub qois: Vec<String>,
    /// Reported uncertainty/repeats; `"unreported"` must be said, not
    /// implied by absence.
    pub reported_uncertainty: String,
    /// Explicit transfer limitations (always nonempty: no result
    /// transfers unconditionally).
    pub transfer_limitations: Vec<String>,
}

impl RegisteredSource {
    /// Validate a complete row.
    ///
    /// # Errors
    /// Fails closed when any semantic field is missing or out of bounds;
    /// analytic models must declare at least one equation.
    pub fn try_new(row: RegisteredSource) -> Result<RegisteredSource, RegistryError> {
        checked_text("source.id", &row.id)?;
        checked_text("source.version_or_date", &row.version_or_date)?;
        checked_text("source.locator", &row.locator)?;
        checked_text("source.license", &row.license)?;
        checked_text("source.reported_uncertainty", &row.reported_uncertainty)?;
        let minimum_equations = usize::from(row.source_type == SourceType::AnalyticModel);
        checked_list(
            "source.model_equations",
            &row.model_equations,
            minimum_equations,
        )?;
        checked_list("source.assumptions", &row.assumptions, 1)?;
        checked_list(
            "source.measured_configurations",
            &row.measured_configurations,
            0,
        )?;
        checked_list("source.qois", &row.qois, 1)?;
        checked_list("source.transfer_limitations", &row.transfer_limitations, 1)?;
        Ok(row)
    }

    /// Declaration identity: moves with EVERY semantic field, so a version
    /// or assumption edit is a new declaration.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        fn push_text(payload: &mut Vec<u8>, text: &str) {
            payload.extend_from_slice(text.as_bytes());
            payload.push(0);
        }
        let mut payload = Vec::new();
        push_text(&mut payload, &self.id);
        payload.push(self.source_type.tag());
        push_text(&mut payload, &self.version_or_date);
        push_text(&mut payload, &self.locator);
        push_text(&mut payload, &self.license);
        for list in [
            &self.model_equations,
            &self.assumptions,
            &self.measured_configurations,
            &self.qois,
        ] {
            payload.extend_from_slice(&(list.len() as u32).to_le_bytes());
            for entry in list.iter() {
                payload.extend_from_slice(entry.as_bytes());
                payload.push(0);
            }
        }
        push_text(&mut payload, &self.reported_uncertainty);
        payload.extend_from_slice(&(self.transfer_limitations.len() as u32).to_le_bytes());
        for entry in &self.transfer_limitations {
            payload.extend_from_slice(entry.as_bytes());
            payload.push(0);
        }
        fs_blake3::hash_domain(SOURCE_DECLARATION_DOMAIN, &payload)
    }
}

/// The rival mechanism classes named by the campaign.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum MechanismClass {
    /// Ideal conservative rolling (the null: no dissipation).
    IdealConservativeRolling,
    /// Phenomenological dry contour moment.
    DryContourMoment,
    /// Phenomenological viscous contour moment.
    ViscousContourMoment,
    /// Finite-patch elastic/viscoelastic/plastic/adhesive contact loss.
    FinitePatchContactLoss,
    /// Tangential microslip/creep in the contact patch.
    TangentialMicroslip,
    /// Gross slip at the contact.
    GrossSlip,
    /// Support/base modal loss.
    BaseModalLoss,
    /// Runout / contact-loss impact cascade.
    RunoutContactImpacts,
    /// Bulk aerodynamic drag on the disc body.
    BulkDrag,
    /// Wedge / thin-gap viscous air flow under the disc.
    WedgeThinGapFlow,
    /// Oscillatory air boundary layer.
    OscillatoryBoundaryLayer,
}

/// Closed interval of decay exponents a mechanism family predicts for the
/// terminal power law (in `(t_f - t)^n` conventions declared by the card).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExponentFamily {
    /// Inclusive lower bound.
    pub lo: f64,
    /// Inclusive upper bound.
    pub hi: f64,
}

impl ExponentFamily {
    fn admitted(&self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }

    /// Whether an observed exponent lies in this family.
    #[must_use]
    pub fn covers(&self, exponent: f64) -> bool {
        exponent.is_finite() && exponent >= self.lo && exponent <= self.hi
    }
}

/// One rival mechanism card.
#[derive(Debug, Clone, PartialEq)]
pub struct MechanismCard {
    /// Stable card id.
    pub id: String,
    /// Mechanism class.
    pub class: MechanismClass,
    /// Observables this mechanism predicts or shapes.
    pub observables: Vec<String>,
    /// Orthogonal interventions that could falsify or distinguish it
    /// (vacuum runs, surface swaps, base isolation, load scaling, ...).
    pub distinguishing_interventions: Vec<String>,
    /// Predicted terminal decay-exponent family (with its convention
    /// stated in the card's observables), when the mechanism predicts one.
    pub exponent_family: Option<ExponentFamily>,
    /// Registered source ids informing this card (informing, never
    /// validating: resolution checks registration, nothing more).
    pub source_ids: Vec<String>,
}

impl MechanismCard {
    /// Validate a card.
    ///
    /// # Errors
    /// Fails closed on missing observables or interventions, malformed
    /// exponent families, or unbounded lists.
    pub fn try_new(card: MechanismCard) -> Result<MechanismCard, RegistryError> {
        checked_text("card.id", &card.id)?;
        checked_list("card.observables", &card.observables, 1)?;
        checked_list(
            "card.distinguishing_interventions",
            &card.distinguishing_interventions,
            1,
        )?;
        if let Some(family) = &card.exponent_family
            && !family.admitted()
        {
            return Err(refuse(
                "euler-registry-exponent-family",
                "card.exponent_family",
                "exponent bounds must be finite with lo <= hi",
            ));
        }
        checked_list("card.source_ids", &card.source_ids, 0)?;
        Ok(card)
    }

    /// Card identity: moves with every semantic field.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut payload = Vec::new();
        payload.extend_from_slice(self.id.as_bytes());
        payload.push(0);
        payload.extend_from_slice(format!("{:?}", self.class).as_bytes());
        payload.push(0);
        for list in [
            &self.observables,
            &self.distinguishing_interventions,
            &self.source_ids,
        ] {
            payload.extend_from_slice(&(list.len() as u32).to_le_bytes());
            for entry in list.iter() {
                payload.extend_from_slice(entry.as_bytes());
                payload.push(0);
            }
        }
        match &self.exponent_family {
            None => payload.push(0),
            Some(family) => {
                payload.push(1);
                payload.extend_from_slice(&family.lo.to_le_bytes());
                payload.extend_from_slice(&family.hi.to_le_bytes());
            }
        }
        fs_blake3::hash_domain(MECHANISM_CARD_DOMAIN, &payload)
    }
}

/// The campaign registry: sources by id, cards resolved against them.
#[derive(Debug, Default)]
pub struct MechanismRegistry {
    sources: BTreeMap<String, RegisteredSource>,
    cards: BTreeMap<String, MechanismCard>,
}

impl MechanismRegistry {
    /// Empty registry.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a validated source row; duplicate ids refuse.
    ///
    /// # Errors
    /// Row validation refusals plus duplicate-id refusal.
    pub fn register_source(&mut self, row: RegisteredSource) -> Result<(), RegistryError> {
        let row = RegisteredSource::try_new(row)?;
        if self.sources.contains_key(&row.id) {
            return Err(refuse(
                "euler-registry-duplicate-source",
                row.id.clone(),
                "source id is already registered; a revision is a NEW id or an explicit replacement",
            ));
        }
        self.sources.insert(row.id.clone(), row);
        Ok(())
    }

    /// Register a validated card, resolving every named source id: a
    /// citation NAME that is not a registered row refuses — names alone
    /// never authorize a model.
    ///
    /// # Errors
    /// Card validation refusals, duplicate card id, unresolved source id.
    pub fn register_card(&mut self, card: MechanismCard) -> Result<(), RegistryError> {
        let card = MechanismCard::try_new(card)?;
        if self.cards.contains_key(&card.id) {
            return Err(refuse(
                "euler-registry-duplicate-card",
                card.id.clone(),
                "card id is already registered",
            ));
        }
        for source_id in &card.source_ids {
            if !self.sources.contains_key(source_id) {
                return Err(refuse(
                    "euler-registry-unresolved-source",
                    source_id.clone(),
                    "a citation name alone cannot authorize a mechanism card; \
                     register the version-pinned source row first",
                ));
            }
        }
        self.cards.insert(card.id.clone(), card);
        Ok(())
    }

    /// Registered source count.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.sources.len()
    }

    /// Registered card count.
    #[must_use]
    pub fn card_count(&self) -> usize {
        self.cards.len()
    }

    /// Every card whose predicted exponent family covers the observation,
    /// in canonical id order. The AMBIGUITY is the scientific result: an
    /// exponent-only observation selecting more than one card cannot
    /// identify a mechanism, and this query makes that undeniable.
    #[must_use]
    pub fn mechanisms_matching_exponent(&self, observed: f64) -> Vec<&MechanismCard> {
        self.cards
            .values()
            .filter(|card| {
                card.exponent_family
                    .as_ref()
                    .is_some_and(|family| family.covers(observed))
            })
            .collect()
    }
}

/// The campaign's principal pinned sources and rival cards, exactly as the
/// bead scopes them. The March 2026 air-drag preprint is deliberately NOT
/// registered here: its exact locator is custodial knowledge this module
/// must not guess, and the fail-closed test proves an incomplete row
/// refuses rather than being quietly admitted.
///
/// # Errors
/// Propagates registration refusals (none for the shipped rows; the
/// batteries prove the set registers clean).
pub fn campaign_registry() -> Result<MechanismRegistry, RegistryError> {
    let mut registry = MechanismRegistry::new();
    let text = |values: &[&str]| values.iter().map(|&v| v.to_string()).collect::<Vec<_>>();

    registry.register_source(RegisteredSource {
        id: "moffatt-2000-nature".into(),
        source_type: SourceType::AnalyticModel,
        version_or_date: "2000-04-20".into(),
        locator: "Moffatt, H.K., Euler's disk and its finite-time singularity, Nature 404, 833-834 (2000), doi:10.1038/35009017".into(),
        license: "publisher-copyright; cited under quotation".into(),
        model_equations: text(&[
            "whole-gap laminar viscous air dissipation Phi ~ pi*mu*g^2*a^4/(2*h*Omega^2*alpha^2) with alpha(t) ~ (t_f - t)^(1/3)",
        ]),
        assumptions: text(&[
            "thin whole-gap laminar air film under the disc dominates dissipation",
            "small inclination asymptotics; rolling without slip; rigid disc and support",
        ]),
        measured_configurations: text(&["order-of-magnitude check against a commercial Euler disc toy"]),
        qois: text(&["inclination alpha(t)", "precession rate Omega(t)", "finite settling time t_f"]),
        reported_uncertainty: "order-of-magnitude agreement only; no formal uncertainty".into(),
        transfer_limitations: text(&[
            "gap-flow dominance is contested by vacuum experiments; not transferable to rough or grooved supports",
            "asymptotic small-alpha regime only",
        ]),
    })?;

    registry.register_source(RegisteredSource {
        id: "van-den-engh-2000-nature".into(),
        source_type: SourceType::ExperimentalStudy,
        version_or_date: "2000-11-30".into(),
        locator: "van den Engh, G., Nelson, P., Roach, J., Numismatic gyrations, Nature 408, 540 (2000), doi:10.1038/35046239".into(),
        license: "publisher-copyright; cited under quotation".into(),
        model_equations: Vec::new(),
        assumptions: text(&["coin spin-down in reduced pressure comparable to air"]),
        measured_configurations: text(&["coins spun under vacuum-jar reduced pressure"]),
        qois: text(&["settling time vs ambient pressure"]),
        reported_uncertainty: "qualitative comparison; no formal uncertainty".into(),
        transfer_limitations: text(&[
            "coins are not machined Euler discs; support and edge finish uncontrolled",
        ]),
    })?;

    registry.register_source(RegisteredSource {
        id: "bildsten-2002-pre".into(),
        source_type: SourceType::AnalyticModel,
        version_or_date: "2002".into(),
        locator: "Bildsten, L., Viscous dissipation for Euler's disk, Phys. Rev. E 66, 056309 (2002), doi:10.1103/PhysRevE.66.056309".into(),
        license: "publisher-copyright; cited under quotation".into(),
        model_equations: text(&[
            "oscillatory (Stokes) boundary-layer dissipation scaling with delta ~ sqrt(nu/Omega); modified finite-time exponent",
        ]),
        assumptions: text(&[
            "boundary-layer thickness smaller than gap height near the end; laminar oscillatory flow",
        ]),
        measured_configurations: Vec::new(),
        qois: text(&["dissipation-rate scaling", "terminal exponent"]),
        reported_uncertainty: "unreported (theory)".into(),
        transfer_limitations: text(&[
            "scaling regime depends on gap-to-boundary-layer ratio; not valid where contact mechanics dominates",
        ]),
    })?;

    registry.register_source(RegisteredSource {
        id: "mcdonald-2001-ajp".into(),
        source_type: SourceType::ExperimentalStudy,
        version_or_date: "2001 (arXiv physics/0008227)".into(),
        locator: "McDonald, A.J., McDonald, K.T., The rolling motion of a disk on a horizontal plane, arXiv:physics/0008227".into(),
        license: "arXiv non-exclusive license".into(),
        model_equations: text(&["rolling-friction dissipation with empirical coefficient"]),
        assumptions: text(&["rolling friction dominates for their disc/surface pair"]),
        measured_configurations: text(&["steel disc on glass and other supports; video timing"]),
        qois: text(&["settling time", "precession-rate growth law"]),
        reported_uncertainty: "run-to-run scatter reported qualitatively".into(),
        transfer_limitations: text(&[
            "friction coefficient is surface-pair-specific; not transferable across finish or cleanliness",
        ]),
    })?;

    registry.register_source(RegisteredSource {
        id: "kessler-oreilly-2002".into(),
        source_type: SourceType::AnalyticModel,
        version_or_date: "2002".into(),
        locator: "Kessler, P., O'Reilly, O.M., The ringing of Euler's disk, Regular and Chaotic Dynamics 7(1), 49-60 (2002)".into(),
        license: "publisher-copyright; cited under quotation".into(),
        model_equations: text(&[
            "rigid-body rolling with dissipation from slip and vibration; contact-loss (ringing) episodes",
        ]),
        assumptions: text(&["dissipation split among slip, rolling resistance, and impacts"]),
        measured_configurations: Vec::new(),
        qois: text(&["energy partition among mechanisms", "contact-loss episodes"]),
        reported_uncertainty: "unreported (theory)".into(),
        transfer_limitations: text(&[
            "mechanism weights are apparatus-dependent by the paper's own construction",
        ]),
    })?;

    registry.register_source(RegisteredSource {
        id: "leine-2009-aam".into(),
        source_type: SourceType::ExperimentalStudy,
        version_or_date: "2009".into(),
        locator: "Leine, R.I., Experimental and theoretical investigation of the energy dissipation of a rolling disk during its final stage of motion, Archive of Applied Mechanics 79, 1063-1082 (2009), doi:10.1007/s00419-008-0278-6".into(),
        license: "publisher-copyright; cited under quotation".into(),
        model_equations: text(&["contour-friction moment models (dry and viscous) fitted to decay data"]),
        assumptions: text(&["contour friction dominates the final stage for the measured disc"]),
        measured_configurations: text(&["machined disc on glass plate; high-speed measurement of the final stage"]),
        qois: text(&["inclination decay law", "fitted contour-moment parameters"]),
        reported_uncertainty: "fit residuals reported per model".into(),
        transfer_limitations: text(&[
            "fitted moments are specific to the measured disc/plate pair and load",
        ]),
    })?;

    registry.register_source(RegisteredSource {
        id: "caps-2004-pre".into(),
        source_type: SourceType::ExperimentalStudy,
        version_or_date: "2004".into(),
        locator: "Caps, H., Dorbolo, S., Ponte, S., Croisier, H., Vandewalle, N., Rolling and slipping motion of Euler's disk, Phys. Rev. E 69, 056610 (2004), doi:10.1103/PhysRevE.69.056610".into(),
        license: "publisher-copyright; cited under quotation".into(),
        model_equations: Vec::new(),
        assumptions: text(&["slip observable via marker tracking through the final stage"]),
        measured_configurations: text(&["disc with optical markers; slip velocity measured directly"]),
        qois: text(&["slip velocity", "rolling/slipping regime boundaries"]),
        reported_uncertainty: "measurement scatter shown in figures".into(),
        transfer_limitations: text(&[
            "slip onset depends on surface pair and spin history; regime map is apparatus-specific",
        ]),
    })?;

    // Rival mechanism cards. Exponent families use the precession-rate
    // convention Omega ~ (t_f - t)^(-n) stated per-card in observables;
    // the OVERLAP among families is deliberate and load-bearing.
    registry.register_card(MechanismCard {
        id: "card-ideal-conservative".into(),
        class: MechanismClass::IdealConservativeRolling,
        observables: text(&["no settling in finite time; constant energy"]),
        distinguishing_interventions: text(&["any observed finite settling falsifies it outright"]),
        exponent_family: None,
        source_ids: Vec::new(),
    })?;
    registry.register_card(MechanismCard {
        id: "card-whole-gap-viscous".into(),
        class: MechanismClass::WedgeThinGapFlow,
        observables: text(&[
            "Omega ~ (t_f - t)^(-n) with n near 1/3 under gap-flow scaling; gap-height sensitivity",
        ]),
        distinguishing_interventions: text(&[
            "vacuum or reduced-pressure runs (kills air terms)",
            "grooved support breaking the thin-gap film",
        ]),
        exponent_family: Some(ExponentFamily { lo: 0.30, hi: 0.40 }),
        source_ids: text(&["moffatt-2000-nature"]),
    })?;
    registry.register_card(MechanismCard {
        id: "card-oscillatory-boundary-layer".into(),
        class: MechanismClass::OscillatoryBoundaryLayer,
        observables: text(&[
            "modified terminal exponent from delta ~ sqrt(nu/Omega) scaling; pressure dependence",
        ]),
        distinguishing_interventions: text(&[
            "gas-species swap changing nu at fixed pressure",
            "reduced-pressure sweep",
        ]),
        exponent_family: Some(ExponentFamily { lo: 0.35, hi: 0.50 }),
        source_ids: text(&["bildsten-2002-pre"]),
    })?;
    registry.register_card(MechanismCard {
        id: "card-dry-contour-moment".into(),
        class: MechanismClass::DryContourMoment,
        observables: text(&[
            "decay law consistent with constant-moment contour friction; load dependence",
        ]),
        distinguishing_interventions: text(&[
            "normal-load scaling via disc mass at fixed geometry",
            "edge-finish swap at fixed air environment",
        ]),
        exponent_family: Some(ExponentFamily { lo: 0.30, hi: 0.45 }),
        source_ids: text(&["leine-2009-aam", "mcdonald-2001-ajp"]),
    })?;
    registry.register_card(MechanismCard {
        id: "card-viscous-contour-moment".into(),
        class: MechanismClass::ViscousContourMoment,
        observables: text(&["rate-dependent contour moment; distinct curvature of the decay tail"]),
        distinguishing_interventions: text(&["spin-history sweep at fixed load and finish"]),
        exponent_family: Some(ExponentFamily { lo: 0.40, hi: 0.60 }),
        source_ids: text(&["leine-2009-aam"]),
    })?;
    registry.register_card(MechanismCard {
        id: "card-microslip".into(),
        class: MechanismClass::TangentialMicroslip,
        observables: text(&["marker-measurable slip velocity in the final stage"]),
        distinguishing_interventions: text(&["optical marker tracking (direct observation)"]),
        exponent_family: None,
        source_ids: text(&["caps-2004-pre"]),
    })?;
    registry.register_card(MechanismCard {
        id: "card-base-modal-loss".into(),
        class: MechanismClass::BaseModalLoss,
        observables: text(&[
            "support-plate ringing correlated with decay; base-material sensitivity",
        ]),
        distinguishing_interventions: text(&[
            "base isolation / massive-rigid support swap",
            "support-material swap at fixed disc",
        ]),
        exponent_family: None,
        source_ids: text(&["kessler-oreilly-2002"]),
    })?;
    registry.register_card(MechanismCard {
        id: "card-contact-loss-impacts".into(),
        class: MechanismClass::RunoutContactImpacts,
        observables: text(&["audible/high-speed-visible contact-loss episodes near termination"]),
        distinguishing_interventions: text(&[
            "high-speed imaging of the final tens of milliseconds",
        ]),
        exponent_family: None,
        source_ids: text(&["kessler-oreilly-2002"]),
    })?;

    Ok(registry)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_source(id: &str) -> RegisteredSource {
        RegisteredSource {
            id: id.into(),
            source_type: SourceType::ExperimentalStudy,
            version_or_date: "2020".into(),
            locator: "Journal 1, 1-2 (2020)".into(),
            license: "cited".into(),
            model_equations: Vec::new(),
            assumptions: vec!["assumption".into()],
            measured_configurations: Vec::new(),
            qois: vec!["qoi".into()],
            reported_uncertainty: "unreported".into(),
            transfer_limitations: vec!["apparatus-specific".into()],
        }
    }

    #[test]
    fn the_campaign_registry_registers_clean_and_is_deterministic() {
        let registry = campaign_registry().expect("shipped rows register");
        assert_eq!(registry.source_count(), 7);
        assert_eq!(registry.card_count(), 8);
        // Identity determinism across rebuilds.
        let again = campaign_registry().expect("registers again");
        let ids: Vec<_> = registry
            .mechanisms_matching_exponent(0.38)
            .iter()
            .map(|card| card.identity())
            .collect();
        let ids_again: Vec<_> = again
            .mechanisms_matching_exponent(0.38)
            .iter()
            .map(|card| card.identity())
            .collect();
        assert_eq!(ids, ids_again);
    }

    #[test]
    fn missing_fields_fail_closed() {
        // Analytic model without equations.
        let mut row = minimal_source("no-eq");
        row.source_type = SourceType::AnalyticModel;
        assert_eq!(
            RegisteredSource::try_new(row)
                .expect_err("analytic needs equations")
                .rule,
            "euler-registry-list-bounds"
        );
        // Missing assumptions.
        let mut row = minimal_source("no-assume");
        row.assumptions.clear();
        assert!(RegisteredSource::try_new(row).is_err());
        // Missing transfer limitations (the no-claim carrier).
        let mut row = minimal_source("no-transfer");
        row.transfer_limitations.clear();
        assert!(RegisteredSource::try_new(row).is_err());
        // Missing version pin.
        let mut row = minimal_source("no-version");
        row.version_or_date = String::new();
        assert!(RegisteredSource::try_new(row).is_err());
        // Card without interventions.
        let card = MechanismCard {
            id: "card-x".into(),
            class: MechanismClass::BulkDrag,
            observables: vec!["obs".into()],
            distinguishing_interventions: Vec::new(),
            exponent_family: None,
            source_ids: Vec::new(),
        };
        assert!(MechanismCard::try_new(card).is_err());
        // Malformed exponent family.
        let card = MechanismCard {
            id: "card-y".into(),
            class: MechanismClass::BulkDrag,
            observables: vec!["obs".into()],
            distinguishing_interventions: vec!["vacuum".into()],
            exponent_family: Some(ExponentFamily { lo: 0.5, hi: 0.2 }),
            source_ids: Vec::new(),
        };
        assert_eq!(
            MechanismCard::try_new(card)
                .expect_err("lo > hi refuses")
                .rule,
            "euler-registry-exponent-family"
        );
    }

    #[test]
    fn a_citation_name_alone_cannot_authorize_a_card() {
        let mut registry = MechanismRegistry::new();
        let card = MechanismCard {
            id: "card-unbacked".into(),
            class: MechanismClass::BulkDrag,
            observables: vec!["obs".into()],
            distinguishing_interventions: vec!["vacuum".into()],
            exponent_family: None,
            source_ids: vec!["moffatt-2000-nature".into()], // famous NAME, unregistered
        };
        let error = registry
            .register_card(card)
            .expect_err("name alone refuses");
        assert_eq!(error.rule, "euler-registry-unresolved-source");
        // Registering the pinned row first admits the card.
        let mut registry = campaign_registry().expect("registry");
        let card = MechanismCard {
            id: "card-backed".into(),
            class: MechanismClass::BulkDrag,
            observables: vec!["obs".into()],
            distinguishing_interventions: vec!["vacuum".into()],
            exponent_family: None,
            source_ids: vec!["moffatt-2000-nature".into()],
        };
        registry
            .register_card(card)
            .expect("registered row authorizes resolution");
    }

    #[test]
    fn version_or_assumption_changes_move_identity() {
        let base = RegisteredSource::try_new(minimal_source("id-move")).expect("admits");
        let base_id = base.identity();
        let mut versioned = base.clone();
        versioned.version_or_date = "2021".into();
        assert_ne!(versioned.identity(), base_id, "version edit moves identity");
        let mut assumed = base.clone();
        assumed.assumptions.push("extra assumption".into());
        assert_ne!(
            assumed.identity(),
            base_id,
            "assumption edit moves identity"
        );
        let mut relabeled = base;
        relabeled.reported_uncertainty = "now reported".into();
        assert_ne!(
            relabeled.identity(),
            base_id,
            "uncertainty edit moves identity"
        );
    }

    #[test]
    fn an_exponent_only_observation_cannot_select_a_unique_mechanism() {
        let registry = campaign_registry().expect("registry");
        // 0.38 sits inside gap-flow, boundary-layer, AND dry-contour
        // families: three rivals, one observation. The ambiguity is the
        // scientific result this registry exists to make undeniable.
        let matches = registry.mechanisms_matching_exponent(0.38);
        assert!(
            matches.len() >= 3,
            "overlapping exponent families must produce ambiguity, got {}",
            matches.len()
        );
        // And an out-of-family observation matches nothing power-law.
        assert!(registry.mechanisms_matching_exponent(5.0).is_empty());
        // Non-finite observations match nothing.
        assert!(registry.mechanisms_matching_exponent(f64::NAN).is_empty());
    }

    #[test]
    fn duplicates_refuse() {
        let mut registry = campaign_registry().expect("registry");
        assert_eq!(
            registry
                .register_source(minimal_source("moffatt-2000-nature"))
                .expect_err("duplicate id refuses")
                .rule,
            "euler-registry-duplicate-source"
        );
    }
}
