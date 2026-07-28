//! Evidence-bearing quantities of interest for the thermal vertical.
//!
//! This module is the narrow consumer seam between the steady conduction
//! field and the enclosure-airflow operating point. It deliberately does not
//! infer temperature error from a residual, treat an estimated airflow band
//! as a physical certificate, or turn a missing receipt into zero. Every QoI
//! carries an exactly-eight-source [`EngineeringUncertaintyBudget`]; sources
//! without a valid propagation path are explicit [`TermValue::Unknown`].

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};
use fs_conduction::{ConductionMesh, ConductionSolution, ProvenanceClass};
use fs_evidence::uncertainty::{
    BudgetTotal, EngineeringUncertaintyBudget, EngineeringUncertaintyKind,
    EngineeringUncertaintyTerm, TermValue, UncertaintyArtifactRef, UncertaintyError,
};
use fs_evidence::{
    Evidence, ModelCard, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate, ValidityDomain, color_of,
};
use fs_qty::{Power, Pressure, Temperature};
use fs_regime::{
    OperatingPoint as RegimeOperatingPoint, OutputAuditBudgetError, OutputAuditError,
    OverrideAcknowledgement, ProductOutputAudit, QoiClaim, RegimeAuditCard,
    apply_output_audit_to_budget, audit_product_output_with_cards,
};

use crate::{OperatingPoint, SourceProvenance};

const QOI_IDENTITY_DOMAIN: &str = "org.frankensim.fs-airflow.thermal-qoi.v1";
const QOI_TERM_DOMAIN: &str = "org.frankensim.fs-airflow.thermal-qoi.term.v1";
const MAX_REGION_ENTRIES: usize = 1_000_000;

/// A declared set of conduction vertices representing one component region.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct JunctionRegion {
    name: String,
    vertices: Vec<usize>,
}

impl JunctionRegion {
    /// Admit a non-empty component region in canonical vertex order.
    pub fn try_new(name: impl Into<String>, mut vertices: Vec<usize>) -> Result<Self, QoiError> {
        let name = admit_name("junction region", name.into())?;
        canonicalize_indices("junction vertices", &mut vertices)?;
        Ok(Self { name, vertices })
    }

    /// Stable region name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Canonically ordered vertex indices.
    #[must_use]
    pub fn vertices(&self) -> &[usize] {
        &self.vertices
    }
}

/// A declared set of boundary faces representing one reporting surface.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SurfaceRegion {
    name: String,
    boundary_faces: Vec<usize>,
}

impl SurfaceRegion {
    /// Admit a non-empty surface region in canonical boundary-face order.
    pub fn try_new(
        name: impl Into<String>,
        mut boundary_faces: Vec<usize>,
    ) -> Result<Self, QoiError> {
        let name = admit_name("surface region", name.into())?;
        canonicalize_indices("surface boundary faces", &mut boundary_faces)?;
        Ok(Self {
            name,
            boundary_faces,
        })
    }

    /// Stable region name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Canonically ordered boundary-face slots.
    #[must_use]
    pub fn boundary_faces(&self) -> &[usize] {
        &self.boundary_faces
    }
}

/// Total fan-efficiency input needed to turn air power into fan input power.
#[derive(Debug, Clone, PartialEq)]
pub struct FanPowerSpec {
    /// Nominal total efficiency in `(0, 1]`.
    total_efficiency: f64,
    /// Declared absolute efficiency half-width. The full interval must remain
    /// inside `(0, 1]`.
    efficiency_half_width: f64,
    /// Authority for the efficiency value and its allowance.
    source: SourceProvenance,
}

impl FanPowerSpec {
    /// Validate a fan-power conversion input.
    pub fn try_new(
        total_efficiency: f64,
        efficiency_half_width: f64,
        source: SourceProvenance,
    ) -> Result<Self, QoiError> {
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(QoiError::invalid(
                "fan efficiency source",
                "citation and stable identifier must both be non-empty",
            ));
        }
        let low = total_efficiency - efficiency_half_width;
        let high = total_efficiency + efficiency_half_width;
        if !(total_efficiency.is_finite()
            && efficiency_half_width.is_finite()
            && efficiency_half_width >= 0.0
            && low > 0.0
            && high <= 1.0)
        {
            return Err(QoiError::invalid(
                "fan efficiency",
                "nominal +/- half-width must be a finite interval inside (0, 1]",
            ));
        }
        Ok(Self {
            total_efficiency,
            efficiency_half_width,
            source,
        })
    }

    /// Nominal total efficiency.
    #[must_use]
    pub const fn total_efficiency(&self) -> f64 {
        self.total_efficiency
    }

    /// Declared absolute efficiency half-width.
    #[must_use]
    pub const fn efficiency_half_width(&self) -> f64 {
        self.efficiency_half_width
    }

    /// Authority for the efficiency value and allowance.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
}

/// Safety-factor authority already reflected in an effective requirement.
///
/// No universal multiply/divide rule is inferred here. The declaring source
/// owns how the finite factor was applied, and
/// [`ThermalRequirement::maximum_temperature`] stores the effective value
/// consumed by margin evaluation. This record exists so the factor and its
/// authority bind into requirement identity; nothing in this crate re-applies
/// it. Applying it a second time here would be a double derating.
#[derive(Debug, Clone, PartialEq)]
pub struct SafetyFactorAuthority {
    /// Applied factor, finite and at least one.
    factor: f64,
    /// Versioned policy defining the application rule.
    source: SourceProvenance,
}

impl SafetyFactorAuthority {
    /// Admit a finite factor of at least one with a cited application policy.
    ///
    /// # Errors
    /// Refuses a non-finite factor, a factor below one, or a blank citation
    /// or identifier.
    pub fn try_new(factor: f64, source: SourceProvenance) -> Result<Self, QoiError> {
        if !(factor.is_finite() && factor >= 1.0) {
            return Err(QoiError::invalid(
                "safety factor",
                "factor must be finite and at least one",
            ));
        }
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(QoiError::invalid(
                "safety factor source",
                "citation and stable identifier must both be non-empty",
            ));
        }
        Ok(Self { factor, source })
    }

    /// Applied factor, retained for identity only.
    #[must_use]
    pub const fn factor(&self) -> f64 {
        self.factor
    }

    /// Authority defining how the factor was applied upstream.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
}

/// Maximum admissible component temperature from a cited requirement.
///
/// The temperature is the *effective* limit: any safety factor was already
/// applied by the declaring authority. The retained
/// [`SafetyFactorAuthority`] records which factor produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalRequirement {
    /// Effective maximum permitted temperature in kelvin, post-factor.
    maximum_temperature: Temperature,
    /// Safety-factor authority already reflected in the value above.
    safety_factor: SafetyFactorAuthority,
    /// Requirement/standard authority.
    source: SourceProvenance,
}

impl ThermalRequirement {
    /// Admit a finite positive effective absolute-temperature requirement.
    ///
    /// # Errors
    /// Refuses a non-finite or non-positive kelvin limit, or a blank
    /// requirement citation or identifier.
    pub fn try_new(
        maximum_temperature: Temperature,
        safety_factor: SafetyFactorAuthority,
        source: SourceProvenance,
    ) -> Result<Self, QoiError> {
        if !(maximum_temperature.value().is_finite() && maximum_temperature.value() > 0.0) {
            return Err(QoiError::invalid(
                "thermal requirement",
                "maximum temperature must be finite and positive in kelvin",
            ));
        }
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(QoiError::invalid(
                "thermal requirement source",
                "citation and stable identifier must both be non-empty",
            ));
        }
        Ok(Self {
            maximum_temperature,
            safety_factor,
            source,
        })
    }

    /// Effective maximum permitted temperature, with the factor already applied.
    #[must_use]
    pub const fn maximum_temperature(&self) -> Temperature {
        self.maximum_temperature
    }

    /// Safety-factor authority behind the effective limit.
    #[must_use]
    pub const fn safety_factor(&self) -> &SafetyFactorAuthority {
        &self.safety_factor
    }

    /// Requirement/standard authority.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
}

/// A caller-declared refinement-ladder half-width for the junction maximum.
///
/// This is an *observation the caller retains*, in kelvin, from an actual
/// mesh-refinement ladder or goal-oriented estimate produced elsewhere. This
/// crate neither computes it nor derives it from a solver residual: a
/// conduction residual measured in watts is not converted into kelvin here.
///
/// A declared ladder observation is not a certificate, so admitting one leaves
/// [`Evidence::numerical`](fs_evidence::Evidence) at
/// [`NumericalKind::NoClaim`](fs_evidence::NumericalKind); it populates only
/// the `Discretization` term of the engineering budget, exactly as a declared
/// fan-efficiency band populates `Parameters`.
#[derive(Debug, Clone, PartialEq)]
pub struct DiscretizationReceipt {
    /// Conservative half-width in kelvin, finite and non-negative.
    half_width: f64,
    /// Authority for the retained ladder or estimator run.
    source: SourceProvenance,
}

impl DiscretizationReceipt {
    /// Admit a finite non-negative kelvin half-width with a cited source.
    ///
    /// # Errors
    /// Refuses a non-finite or negative half-width, negative zero, or a blank
    /// citation or identifier.
    pub fn try_new(half_width: f64, source: SourceProvenance) -> Result<Self, QoiError> {
        if !(half_width.is_finite() && half_width >= 0.0)
            || (half_width == 0.0 && half_width.is_sign_negative())
        {
            return Err(QoiError::invalid(
                "discretization receipt",
                "half-width must be finite, non-negative, and not negative zero",
            ));
        }
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(QoiError::invalid(
                "discretization receipt source",
                "citation and stable identifier must both be non-empty",
            ));
        }
        Ok(Self { half_width, source })
    }

    /// Conservative kelvin half-width.
    #[must_use]
    pub const fn half_width(&self) -> f64 {
        self.half_width
    }

    /// Authority for the retained ladder or estimator run.
    #[must_use]
    pub const fn source(&self) -> &SourceProvenance {
        &self.source
    }
}

/// Everything the caller declares for one thermal QoI extraction.
///
/// Grouping these keeps the extraction entry point narrow as QoI families grow,
/// and keeps "what the caller declared" separable from "what the solvers
/// produced" (the mesh, field, and operating point passed alongside it).
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalQoiDeclarations<'a> {
    /// Component region the junction maximum is taken over.
    pub junction_region: &'a JunctionRegion,
    /// Reporting surface the uniformity family is taken over.
    pub surface_region: &'a SurfaceRegion,
    /// Cited fan total-efficiency interval.
    pub fan_power: &'a FanPowerSpec,
    /// Cited effective temperature requirement. Absence refuses; there is no
    /// default limit.
    pub requirement: Option<&'a ThermalRequirement>,
    /// Optional declared refinement-ladder half-width for the junction
    /// maximum. Absence leaves the discretization term a named `Unknown`.
    pub discretization: Option<&'a DiscretizationReceipt>,
}

/// One typed QoI with its legacy evidence carrier and rich engineering budget.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalQoi<T> {
    /// Typed value plus the existing numerical/statistical/model evidence.
    pub evidence: Evidence<T>,
    /// Exactly-eight-source engineering uncertainty budget.
    pub uncertainty: EngineeringUncertaintyBudget,
}

/// Deterministically selected maximum and its canonical vertex witness.
#[derive(Debug, Clone, PartialEq)]
pub struct JunctionMaximum {
    /// Evidence-bearing maximum temperature.
    pub qoi: ThermalQoi<Temperature>,
    /// Lowest canonical vertex index among equal maxima.
    pub vertex: usize,
}

/// Area-weighted surface temperature and two uniformity diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceUniformity {
    /// Exact P1 face-integral mean divided by selected area.
    pub mean_temperature: ThermalQoi<Temperature>,
    /// Maximum minus minimum selected surface-vertex temperature.
    pub spread: ThermalQoi<Temperature>,
    /// Area-weighted standard deviation of face-mean temperatures.
    pub face_mean_standard_deviation: ThermalQoi<Temperature>,
}

/// The five E05.10 QoI families emitted together from one source snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalQoiSet {
    /// Maximum temperature over the declared component region.
    pub junction_maximum: JunctionMaximum,
    /// Enclosure pressure drop at the solved operating point.
    pub pressure_drop: ThermalQoi<Pressure>,
    /// Fan input power `Delta p * Q / eta_total`.
    pub fan_power: ThermalQoi<Power>,
    /// Surface mean/spread/std family.
    pub uniformity: SurfaceUniformity,
    /// `requirement maximum - junction maximum`.
    pub thermal_margin: ThermalQoi<Temperature>,
}

/// Explicit model-card use declaration for one thermal QoI.
///
/// The declaration supplies card identities and an optional non-restoring
/// override acknowledgement. The QoI color is always derived from the actual
/// [`Evidence`] receipt and cannot be supplied by the caller.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThermalQoiCardUse {
    /// Exact QoI identity, matching its engineering uncertainty budget.
    pub qoi: String,
    /// Every model card consumed by this QoI.
    pub model_cards: Vec<String>,
    /// Explicit authorization to proceed despite demotion, if any.
    pub override_acknowledgement: Option<OverrideAcknowledgement>,
}

/// The final E05.10 output: values/budgets plus the complete regime audit.
#[derive(Debug, Clone, PartialEq)]
pub struct AuditedThermalQoiSet {
    /// Original QoI values with every demotion applied to its rich budget.
    pub qois: ThermalQoiSet,
    /// One final-envelope receipt per QoI.
    pub audit: ProductOutputAudit,
}

/// Refusal from final operating-envelope admission of the thermal QoI set.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermalOutputAuditError {
    /// Two card-use declarations named the same QoI.
    DuplicateCardUse {
        /// Repeated QoI identity.
        qoi: String,
    },
    /// One emitted QoI had no consumed-card declaration.
    MissingCardUse {
        /// Emitted QoI identity.
        qoi: String,
    },
    /// A declaration named no emitted thermal QoI.
    UnknownQoi {
        /// Foreign QoI identity.
        qoi: String,
    },
    /// The shared audit returned no receipt for an emitted QoI.
    MissingReceipt {
        /// Emitted QoI identity.
        qoi: String,
    },
    /// The shared `fs-regime` product audit refused the envelope.
    OutputAudit(OutputAuditError),
    /// A receipt could not be applied to its matching rich budget.
    Budget(OutputAuditBudgetError),
}

impl fmt::Display for ThermalOutputAuditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateCardUse { qoi } => {
                write!(f, "duplicate thermal QoI model-card declaration {qoi:?}")
            }
            Self::MissingCardUse { qoi } => {
                write!(f, "thermal QoI {qoi:?} has no model-card declaration")
            }
            Self::UnknownQoi { qoi } => {
                write!(
                    f,
                    "model-card declaration names unknown thermal QoI {qoi:?}"
                )
            }
            Self::MissingReceipt { qoi } => {
                write!(
                    f,
                    "thermal product audit returned no receipt for QoI {qoi:?}"
                )
            }
            Self::OutputAudit(error) => write!(f, "thermal product audit refused: {error}"),
            Self::Budget(error) => write!(f, "thermal budget demotion refused: {error}"),
        }
    }
}

impl std::error::Error for ThermalOutputAuditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::DuplicateCardUse { .. }
            | Self::MissingCardUse { .. }
            | Self::UnknownQoi { .. }
            | Self::MissingReceipt { .. } => None,
            Self::OutputAudit(error) => Some(error),
            Self::Budget(error) => Some(error),
        }
    }
}

impl From<OutputAuditError> for ThermalOutputAuditError {
    fn from(error: OutputAuditError) -> Self {
        Self::OutputAudit(error)
    }
}

impl From<OutputAuditBudgetError> for ThermalOutputAuditError {
    fn from(error: OutputAuditBudgetError) -> Self {
        Self::Budget(error)
    }
}

impl ThermalQoiSet {
    /// Every emitted budget, including all three uniformity records.
    #[must_use]
    pub fn budgets(&self) -> [&EngineeringUncertaintyBudget; 7] {
        [
            &self.junction_maximum.qoi.uncertainty,
            &self.pressure_drop.uncertainty,
            &self.fan_power.uncertainty,
            &self.uniformity.mean_temperature.uncertainty,
            &self.uniformity.spread.uncertainty,
            &self.uniformity.face_mean_standard_deviation.uncertainty,
            &self.thermal_margin.uncertainty,
        ]
    }

    /// True only when every QoI retains an explicit unknown source. This is a
    /// useful integration assertion while the vertical lacks complete DWR,
    /// model-validation, and measurement propagation.
    #[must_use]
    pub fn all_totals_are_honestly_unknown(&self) -> bool {
        self.budgets()
            .iter()
            .all(|budget| matches!(budget.total(), BudgetTotal::Unknown { .. }))
    }

    /// Run the mandatory final operating-envelope audit for all seven records.
    ///
    /// Every emitted QoI must have exactly one [`ThermalQoiCardUse`]. The
    /// incoming color is derived from the actual evidence receipt, all claims
    /// are audited together against every supplied point/card, and each exact
    /// receipt is applied to the matching eight-term budget before return.
    ///
    /// # Errors
    ///
    /// Refuses missing, duplicate, or foreign card-use declarations, any
    /// shared `fs-regime` audit refusal, or a receipt/budget mismatch.
    pub fn audit_operating_envelope(
        self,
        registry: &[ModelCard],
        operating_points: &[RegimeOperatingPoint],
        card_uses: &[ThermalQoiCardUse],
    ) -> Result<AuditedThermalQoiSet, ThermalOutputAuditError> {
        let audit_cards = registry
            .iter()
            .map(RegimeAuditCard::from)
            .collect::<Vec<_>>();
        self.audit_operating_envelope_with_cards(&audit_cards, operating_points, card_uses)
    }

    /// Run the final audit with owner-neutral identity/validity projections.
    ///
    /// This entry point admits material and other domain-native card schemas
    /// without inventing the extra authority fields of an evidence
    /// [`ModelCard`]. Each projection must preserve its owner's exact identity,
    /// version, and shared validity box. The ordinary
    /// [`Self::audit_operating_envelope`] method is an exact wrapper over this
    /// path for existing evidence cards.
    ///
    /// # Errors
    ///
    /// Refuses missing, duplicate, or foreign card-use declarations, any
    /// shared `fs-regime` audit refusal, or a receipt/budget mismatch.
    pub fn audit_operating_envelope_with_cards(
        mut self,
        registry: &[RegimeAuditCard],
        operating_points: &[RegimeOperatingPoint],
        card_uses: &[ThermalQoiCardUse],
    ) -> Result<AuditedThermalQoiSet, ThermalOutputAuditError> {
        let mut uses = BTreeMap::new();
        for declaration in card_uses {
            if uses.insert(declaration.qoi.clone(), declaration).is_some() {
                return Err(ThermalOutputAuditError::DuplicateCardUse {
                    qoi: declaration.qoi.clone(),
                });
            }
        }

        let mut claims = Vec::with_capacity(7);
        push_regime_claim(&mut claims, &mut uses, &self.junction_maximum.qoi)?;
        push_regime_claim(&mut claims, &mut uses, &self.pressure_drop)?;
        push_regime_claim(&mut claims, &mut uses, &self.fan_power)?;
        push_regime_claim(&mut claims, &mut uses, &self.uniformity.mean_temperature)?;
        push_regime_claim(&mut claims, &mut uses, &self.uniformity.spread)?;
        push_regime_claim(
            &mut claims,
            &mut uses,
            &self.uniformity.face_mean_standard_deviation,
        )?;
        push_regime_claim(&mut claims, &mut uses, &self.thermal_margin)?;
        if let Some(qoi) = uses.keys().next() {
            return Err(ThermalOutputAuditError::UnknownQoi { qoi: qoi.clone() });
        }

        let audit = audit_product_output_with_cards(registry, operating_points, &claims)?;
        apply_receipt_to_budget(&audit, &mut self.junction_maximum.qoi.uncertainty)?;
        apply_receipt_to_budget(&audit, &mut self.pressure_drop.uncertainty)?;
        apply_receipt_to_budget(&audit, &mut self.fan_power.uncertainty)?;
        apply_receipt_to_budget(&audit, &mut self.uniformity.mean_temperature.uncertainty)?;
        apply_receipt_to_budget(&audit, &mut self.uniformity.spread.uncertainty)?;
        apply_receipt_to_budget(
            &audit,
            &mut self.uniformity.face_mean_standard_deviation.uncertainty,
        )?;
        apply_receipt_to_budget(&audit, &mut self.thermal_margin.uncertainty)?;

        Ok(AuditedThermalQoiSet { qois: self, audit })
    }
}

fn push_regime_claim<'a, T>(
    claims: &mut Vec<QoiClaim>,
    uses: &mut BTreeMap<String, &'a ThermalQoiCardUse>,
    qoi: &ThermalQoi<T>,
) -> Result<(), ThermalOutputAuditError> {
    let qoi_id = qoi.uncertainty.qoi();
    let declaration =
        uses.remove(qoi_id)
            .ok_or_else(|| ThermalOutputAuditError::MissingCardUse {
                qoi: qoi_id.to_string(),
            })?;
    claims.push(QoiClaim {
        qoi: qoi_id.to_string(),
        color: color_of(&qoi.evidence.numerical, &qoi.evidence.model),
        model_cards: declaration.model_cards.clone(),
        override_acknowledgement: declaration.override_acknowledgement.clone(),
    });
    Ok(())
}

fn apply_receipt_to_budget(
    audit: &ProductOutputAudit,
    budget: &mut EngineeringUncertaintyBudget,
) -> Result<(), ThermalOutputAuditError> {
    let receipt = audit
        .receipts
        .iter()
        .find(|receipt| receipt.qoi == budget.qoi())
        .ok_or_else(|| ThermalOutputAuditError::MissingReceipt {
            qoi: budget.qoi().to_string(),
        })?;
    *budget = apply_output_audit_to_budget(receipt, budget)?;
    Ok(())
}

/// Deterministic refusal from thermal QoI extraction.
#[derive(Debug, Clone, PartialEq)]
pub enum QoiError {
    /// A caller declaration or upstream record is malformed.
    InvalidInput {
        /// Stable input category.
        field: &'static str,
        /// Actionable diagnosis.
        detail: String,
    },
    /// Margin extraction has no cited requirement.
    MissingRequirement,
    /// The eight-term evidence layer refused a malformed budget.
    Uncertainty(UncertaintyError),
}

impl QoiError {
    fn invalid(field: &'static str, detail: impl Into<String>) -> Self {
        Self::InvalidInput {
            field,
            detail: detail.into(),
        }
    }
}

impl fmt::Display for QoiError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput { field, detail } => write!(f, "invalid {field}: {detail}"),
            Self::MissingRequirement => write!(
                f,
                "thermal margin requires a cited maximum-temperature requirement; no default is admissible"
            ),
            Self::Uncertainty(error) => write!(f, "uncertainty budget refused: {error}"),
        }
    }
}

impl std::error::Error for QoiError {}

impl From<UncertaintyError> for QoiError {
    fn from(error: UncertaintyError) -> Self {
        Self::Uncertainty(error)
    }
}

/// Extract all E05.10 QoI families from one conduction/airflow snapshot.
///
/// The field values are useful estimates, but temperature extrema and surface
/// metrics receive `NumericalKind::NoClaim` until a DWR/refinement path binds
/// an actual error enclosure. Pressure and fan power retain the upstream
/// estimate bands. Every budget separately retains an unknown model-form term,
/// so conditional parameter/BC envelopes cannot become product authority.
///
/// An optional [`DiscretizationReceipt`] is the only route by which the
/// junction maximum acquires a known `Discretization` term; the thermal margin
/// then inherits it, because `margin = limit - maximum` has unit sensitivity to
/// the maximum in the same kelvin. Admitting a receipt does not upgrade
/// `NumericalKind::NoClaim`: a declared ladder observation is not a certificate.
///
/// The requirement's [`SafetyFactorAuthority`] is retained for identity only.
/// The supplied limit is already the effective post-factor value, so this
/// function never re-applies the factor.
///
/// # Errors
/// Refuses malformed regions, mismatched/non-finite fields, invalid operating
/// evidence, invalid fan efficiency, a non-finite margin, or an absent
/// requirement.
pub fn extract_thermal_qois(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    operating_point: &OperatingPoint,
    declarations: &ThermalQoiDeclarations<'_>,
) -> Result<ThermalQoiSet, QoiError> {
    let &ThermalQoiDeclarations {
        junction_region,
        surface_region,
        fan_power,
        requirement,
        discretization,
    } = declarations;
    let requirement = requirement.ok_or(QoiError::MissingRequirement)?;
    validate_solution(mesh, solution)?;
    validate_operating_point(operating_point)?;
    validate_region_indices(mesh, junction_region, surface_region)?;

    let solution_id = solution_identity(mesh, solution);
    let operating_id = operating_identity(operating_point);
    let fan_id = fan_power_identity(fan_power);
    let requirement_id = requirement_identity(requirement);
    let discretization_id = discretization.map(discretization_identity);

    let (maximum, maximum_vertex) = junction_maximum(solution, junction_region);
    let temperature_model = conduction_model(solution);
    let mut maximum_parents = vec![solution_id];
    if let Some(identity) = discretization_id {
        maximum_parents.push(identity);
    }
    let maximum_identity = qoi_identity(
        "junction-maximum",
        &maximum_parents,
        &region_identity(junction_region.name(), junction_region.vertices()),
    );
    // A declared refinement-ladder half-width is the ONLY route by which the
    // junction maximum acquires a known Discretization term. Absent a receipt
    // the term stays a named Unknown; nothing here manufactures one.
    // Every temperature QoI carries the same provenance-specific Parameters
    // gap, so a caller-declared conductivity can never read like a
    // receipt-backed one.
    let parameter_term = material_parameter_term(solution, solution_id)?;
    let mut maximum_known = vec![parameter_term.clone()];
    if let (Some(receipt), Some(identity)) = (discretization, discretization_id) {
        maximum_known.push((
            EngineeringUncertaintyKind::Discretization,
            TermValue::interval(0.0, receipt.half_width())?,
            "thermal-qoi-discretization-receipt",
            identity,
        ));
    }
    let maximum_budget = unknown_budget(
        "thermal-junction-maximum",
        "kelvin",
        maximum_identity,
        &maximum_known,
    )?;
    let maximum_evidence =
        no_claim_temperature(maximum, temperature_model.clone(), maximum_identity);
    let junction_maximum = JunctionMaximum {
        qoi: ThermalQoi {
            evidence: maximum_evidence,
            uncertainty: maximum_budget,
        },
        vertex: maximum_vertex,
    };

    let surface = surface_statistics(mesh, solution, surface_region)?;
    let surface_region_id = region_identity(surface_region.name(), surface_region.boundary_faces());
    let mean_identity = qoi_identity("surface-mean", &[solution_id], &surface_region_id);
    let spread_identity = qoi_identity("surface-spread", &[solution_id], &surface_region_id);
    let std_identity = qoi_identity("surface-face-mean-std", &[solution_id], &surface_region_id);
    let uniformity = SurfaceUniformity {
        mean_temperature: ThermalQoi {
            evidence: no_claim_temperature(
                Temperature::new(surface.mean),
                temperature_model.clone(),
                mean_identity,
            ),
            uncertainty: unknown_budget(
                "thermal-surface-mean",
                "kelvin",
                mean_identity,
                std::slice::from_ref(&parameter_term),
            )?,
        },
        spread: ThermalQoi {
            evidence: no_claim_temperature(
                Temperature::new(surface.spread),
                temperature_model.clone(),
                spread_identity,
            ),
            uncertainty: unknown_budget(
                "thermal-surface-spread",
                "kelvin",
                spread_identity,
                std::slice::from_ref(&parameter_term),
            )?,
        },
        face_mean_standard_deviation: ThermalQoi {
            evidence: no_claim_temperature(
                Temperature::new(surface.standard_deviation),
                temperature_model.clone(),
                std_identity,
            ),
            uncertainty: unknown_budget(
                "thermal-surface-face-mean-std",
                "kelvin",
                std_identity,
                std::slice::from_ref(&parameter_term),
            )?,
        },
    };

    let pressure_value = operating_point.pressure.value.value();
    let pressure_half_width = interval_half_width(
        pressure_value,
        operating_point.pressure.numerical.lo,
        operating_point.pressure.numerical.hi,
    )?;
    let pressure_identity = qoi_identity("pressure-drop", &[operating_id], b"operating-pressure");
    let pressure_term = TermValue::interval(0.0, pressure_half_width)?;
    let pressure_drop = ThermalQoi {
        evidence: operating_point.pressure.clone(),
        uncertainty: unknown_budget(
            "thermal-pressure-drop",
            "pascal",
            pressure_identity,
            &[(
                EngineeringUncertaintyKind::BoundaryConditions,
                pressure_term,
                "thermal-qoi-operating-envelope",
                operating_id,
            )],
        )?,
    };

    let fan_power_qoi = fan_power_qoi(operating_point, fan_power, operating_id, fan_id)?;

    // The effective limit already reflects the declared safety factor; this
    // subtraction is the ONLY place the factor could be applied a second time,
    // and it deliberately is not. `safety_factor` reaches `requirement_id`
    // only, so changing the factor rebinds identity without moving the value.
    let margin = requirement.maximum_temperature.value() - maximum.value();
    if !margin.is_finite() {
        return Err(QoiError::invalid(
            "thermal margin",
            "requirement minus junction maximum was non-finite",
        ));
    }
    let margin_identity = qoi_identity(
        "thermal-margin",
        &[maximum_identity, requirement_id],
        requirement.source.identifier.as_bytes(),
    );
    // margin = L_eff - T_max has unit sensitivity to the junction maximum and
    // is expressed in the same kelvin, so every junction-maximum term
    // transfers 1:1. Propagating the term VALUES (rather than re-deriving
    // independent Unknowns) is what makes the margin a functional of its input
    // instead of a scalar wearing a budget. The declared limit contributes no
    // term of its own: it is exact by declaration, not a measured quantity.
    let thermal_margin = ThermalQoi {
        evidence: no_claim_temperature(
            Temperature::new(margin),
            temperature_model,
            margin_identity,
        ),
        uncertainty: propagate_budget(
            &junction_maximum.qoi.uncertainty,
            "thermal-margin",
            "kelvin",
            margin_identity,
            "thermal-qoi-margin-from-junction-maximum",
        )?,
    };

    Ok(ThermalQoiSet {
        junction_maximum,
        pressure_drop,
        fan_power: fan_power_qoi,
        uniformity,
        thermal_margin,
    })
}

fn fan_power_qoi(
    operating_point: &OperatingPoint,
    spec: &FanPowerSpec,
    operating_id: ContentHash,
    fan_id: ContentHash,
) -> Result<ThermalQoi<Power>, QoiError> {
    let flow = &operating_point.flow;
    let pressure = &operating_point.pressure;
    let eta = spec.total_efficiency;
    let eta_low = eta - spec.efficiency_half_width;
    let eta_high = eta + spec.efficiency_half_width;
    let nominal_air_power = pressure.value.value() * flow.value.value();
    let nominal = nominal_air_power / eta;
    let low = pressure.numerical.lo.max(0.0) * flow.numerical.lo.max(0.0) / eta_high;
    let high = pressure.numerical.hi * flow.numerical.hi / eta_low;
    if !(nominal.is_finite()
        && low.is_finite()
        && high.is_finite()
        && 0.0 <= low
        && low <= nominal
        && nominal <= high)
    {
        return Err(QoiError::invalid(
            "fan power",
            "pressure/flow/efficiency envelope did not produce an ordered finite non-negative power interval",
        ));
    }

    let boundary_low = pressure.numerical.lo.max(0.0) * flow.numerical.lo.max(0.0) / eta;
    let boundary_high = pressure.numerical.hi * flow.numerical.hi / eta;
    let boundary_half_width = interval_half_width(nominal, boundary_low, boundary_high)?;
    let efficiency_low_power = nominal_air_power / eta_high;
    let efficiency_high_power = nominal_air_power / eta_low;
    let efficiency_half_width =
        interval_half_width(nominal, efficiency_low_power, efficiency_high_power)?;
    let identity = qoi_identity(
        "fan-power",
        &[operating_id, fan_id],
        b"pressure-flow-efficiency",
    );

    let mut model = pressure.model.clone();
    model
        .cards
        .push(format!("fan-efficiency:{}", spec.source.identifier));
    model.cards.sort_unstable();
    model.cards.dedup();
    model.assumptions.push(
        "fan input power uses Delta-p times volume-flow divided by declared total efficiency"
            .to_string(),
    );
    model.assumptions.sort_unstable();
    model.assumptions.dedup();
    let leaf = ProvenanceHash::of_bytes(identity.as_bytes());
    let provenance = ProvenanceHash::chain(
        "thermal-qoi-fan-power-v1",
        &[pressure.provenance, flow.provenance, leaf],
    );
    let evidence = Evidence {
        value: Power::new(nominal),
        qoi: nominal,
        numerical: NumericalCertificate::estimate(low, high),
        statistical: StatisticalCertificate::None,
        model,
        sensitivity: SensitivitySummary::default(),
        provenance,
        adjoint_ref: None,
    };
    let uncertainty = unknown_budget(
        "thermal-fan-input-power",
        "watt",
        identity,
        &[
            (
                EngineeringUncertaintyKind::Parameters,
                TermValue::interval(0.0, efficiency_half_width)?,
                "thermal-qoi-efficiency-envelope",
                fan_id,
            ),
            (
                EngineeringUncertaintyKind::BoundaryConditions,
                TermValue::interval(0.0, boundary_half_width)?,
                "thermal-qoi-operating-envelope",
                operating_id,
            ),
        ],
    )?;
    Ok(ThermalQoi {
        evidence,
        uncertainty,
    })
}

fn no_claim_temperature(
    value: Temperature,
    model: ModelEvidence,
    identity: ContentHash,
) -> Evidence<Temperature> {
    Evidence {
        value,
        qoi: value.value(),
        numerical: NumericalCertificate::no_claim(),
        statistical: StatisticalCertificate::None,
        model,
        sensitivity: SensitivitySummary::default(),
        provenance: ProvenanceHash::of_bytes(identity.as_bytes()),
        adjoint_ref: None,
    }
}

/// The card naming where this solve's conductivity numbers came from.
///
/// Declared and receipt-backed material provenance must not present the same
/// machine-readable evidence: `ProvenanceClass::Declared` says there is no
/// material provenance at all, which is strictly weaker than a retained
/// `fs-matdb` receipt. Before this existed the class reached only the identity
/// hash, so two solves with different material authority were distinguishable
/// by digest but read identically to a reviewer.
fn material_provenance_card(provenance: ProvenanceClass) -> &'static str {
    match provenance {
        ProvenanceClass::MatdbReceipts => "fs-conduction:material-matdb-receipts",
        ProvenanceClass::Declared => "fs-conduction:material-declared",
    }
}

/// The `Parameters` gap reason, specialized by material authority.
fn material_parameter_gap(provenance: ProvenanceClass) -> &'static str {
    match provenance {
        ProvenanceClass::MatdbReceipts => {
            "material receipts are retained but no complete parameter propagation reaches this QoI"
        }
        ProvenanceClass::Declared => {
            "the conductivity values are caller-declared with no material receipt, so there is no parameter authority to propagate"
        }
    }
}

fn conduction_model(solution: &ConductionSolution) -> ModelEvidence {
    let provenance = solution.report.material_provenance;
    ModelEvidence {
        cards: vec![
            "fs-conduction:steady-p1".to_string(),
            material_provenance_card(provenance).to_string(),
        ],
        assumptions: vec![format!(
            "temperature QoI is a direct functional of a steady P1 field whose conductivity provenance is {} with {} material receipt(s); no DWR-to-QoI bound is attached",
            provenance.tag(),
            solution.report.material_receipts
        )],
        validity: ValidityDomain::unconstrained(),
        discrepancy_rel: f64::INFINITY,
        in_domain: true,
    }
}

/// The provenance-specific `Parameters` term every temperature QoI carries.
///
/// It stays `Unknown` in both cases — neither authority supplies a propagation
/// theorem — but the named reason distinguishes "receipts exist, propagation
/// does not" from "no material authority exists at all".
fn material_parameter_term(
    solution: &ConductionSolution,
    solution_id: ContentHash,
) -> Result<
    (
        EngineeringUncertaintyKind,
        TermValue,
        &'static str,
        ContentHash,
    ),
    QoiError,
> {
    Ok((
        EngineeringUncertaintyKind::Parameters,
        TermValue::unknown(material_parameter_gap(solution.report.material_provenance))?,
        "thermal-qoi-material-provenance",
        solution_id,
    ))
}

fn unknown_budget(
    qoi: &str,
    unit: &str,
    qoi_identity: ContentHash,
    known: &[(
        EngineeringUncertaintyKind,
        TermValue,
        &'static str,
        ContentHash,
    )],
) -> Result<EngineeringUncertaintyBudget, QoiError> {
    let mut terms = Vec::with_capacity(EngineeringUncertaintyKind::ALL.len());
    for kind in EngineeringUncertaintyKind::ALL {
        let (value, role, digest) = known
            .iter()
            .find(|(candidate, _, _, _)| *candidate == kind)
            .map(|(_, value, role, digest)| (value.clone(), *role, *digest))
            .unwrap_or((
                TermValue::unknown(unknown_reason(kind))?,
                "thermal-qoi-evidence-gap",
                term_identity(qoi_identity, kind),
            ));
        let provenance = UncertaintyArtifactRef::new(role, digest)?;
        terms.push(EngineeringUncertaintyTerm::try_new(
            kind, value, provenance,
        )?);
    }
    EngineeringUncertaintyBudget::try_new(qoi, unit, terms).map_err(Into::into)
}

/// Re-express a source budget as a downstream QoI with unit sensitivity.
///
/// Every term VALUE transfers unchanged, which is exact when the downstream
/// QoI is `constant ± source` in the same unit. Provenance is rebound to the
/// downstream identity under `role` so that a change anywhere in the
/// downstream identity chain (for a margin: the requirement and its safety
/// factor) still rebinds the budget, while lineage back to the source stays
/// carried by `identity`'s parent hashes.
///
/// This deliberately cannot strengthen a term: it copies representations
/// rather than composing them, so an `Unknown` source term stays `Unknown`.
fn propagate_budget(
    source: &EngineeringUncertaintyBudget,
    qoi: &str,
    unit: &str,
    identity: ContentHash,
    role: &'static str,
) -> Result<EngineeringUncertaintyBudget, QoiError> {
    let mut terms = Vec::with_capacity(EngineeringUncertaintyKind::ALL.len());
    for kind in EngineeringUncertaintyKind::ALL {
        let provenance = UncertaintyArtifactRef::new(role, term_identity(identity, kind))?;
        terms.push(EngineeringUncertaintyTerm::try_new(
            kind,
            source.term(kind).value().clone(),
            provenance,
        )?);
    }
    EngineeringUncertaintyBudget::try_new(qoi, unit, terms).map_err(Into::into)
}

fn unknown_reason(kind: EngineeringUncertaintyKind) -> &'static str {
    match kind {
        EngineeringUncertaintyKind::Roundoff => {
            "no retained outward-roundoff propagation receipt reaches this QoI"
        }
        EngineeringUncertaintyKind::SolverAlgebraic => {
            "the recomputed residual has no admitted inverse-stability map into this QoI unit"
        }
        EngineeringUncertaintyKind::Discretization => {
            "no retained DWR or refinement-ladder bound is attached to this QoI"
        }
        EngineeringUncertaintyKind::Geometry => {
            "no import, meshing, registration, or as-built geometry uncertainty is propagated"
        }
        EngineeringUncertaintyKind::Parameters => {
            "no complete material and manufacturing parameter propagation is attached"
        }
        EngineeringUncertaintyKind::BoundaryConditions => {
            "no complete boundary and operating-condition propagation is attached"
        }
        EngineeringUncertaintyKind::ModelForm => {
            "no externally validated model-form discrepancy envelope is attached"
        }
        EngineeringUncertaintyKind::Measurement => {
            "no calibrated observation or comparison-data uncertainty is attached"
        }
        _ => "no authoritative propagation is attached for this uncertainty family",
    }
}

fn validate_solution(mesh: &ConductionMesh, solution: &ConductionSolution) -> Result<(), QoiError> {
    if solution.temperature.len() != mesh.vertex_count() {
        return Err(QoiError::invalid(
            "temperature field",
            format!(
                "{} nodal values for {} mesh vertices",
                solution.temperature.len(),
                mesh.vertex_count()
            ),
        ));
    }
    if solution.temperature.iter().any(|value| !value.is_finite()) {
        return Err(QoiError::invalid(
            "temperature field",
            "every nodal temperature must be finite",
        ));
    }
    // `ConductionSolution` has public fields, so this seam validates the
    // report the same way it already validates the field. The production
    // solver fails closed on a non-converged linear solve
    // (`ConductionError::LinearSolveFailed`), so a caller reaching here with
    // one has assembled the report by hand; a QoI taken from an algebraically
    // unconverged field is unsupported, not merely weaker.
    if let Some(record) = solution
        .report
        .linear
        .iter()
        .find(|record| !record.converged_true)
    {
        return Err(QoiError::invalid(
            "conduction report",
            format!(
                "linear solve at nonlinear iteration {} ({}) reported true relative residual {} without converging{}",
                record.nonlinear_iteration,
                record.method,
                record.true_relative_residual,
                match &record.stall {
                    Some(stall) => format!("; retained stall diagnosis {stall:?}"),
                    None => String::new(),
                }
            ),
        ));
    }
    // A claim of retained receipts with no receipts is self-contradictory.
    if solution.report.material_provenance == ProvenanceClass::MatdbReceipts
        && solution.report.material_receipts == 0
    {
        return Err(QoiError::invalid(
            "conduction report",
            "material provenance claims retained matdb receipts while reporting zero of them",
        ));
    }
    Ok(())
}

fn validate_operating_point(point: &OperatingPoint) -> Result<(), QoiError> {
    validate_estimate(
        "operating pressure",
        point.pressure.value.value(),
        point.pressure.numerical.lo,
        point.pressure.numerical.hi,
    )?;
    validate_estimate(
        "operating flow",
        point.flow.value.value(),
        point.flow.numerical.lo,
        point.flow.numerical.hi,
    )?;
    if point.pressure.value.value() < 0.0
        || point.pressure.numerical.lo < 0.0
        || point.flow.value.value() < 0.0
        || point.flow.numerical.lo < 0.0
    {
        return Err(QoiError::invalid(
            "operating point",
            "pressure and volume-flow estimates must be non-negative",
        ));
    }
    Ok(())
}

fn validate_estimate(field: &'static str, value: f64, lo: f64, hi: f64) -> Result<(), QoiError> {
    if !(value.is_finite() && lo.is_finite() && hi.is_finite() && lo <= value && value <= hi) {
        return Err(QoiError::invalid(
            field,
            "value and bounds must be finite, ordered, and self-containing",
        ));
    }
    Ok(())
}

fn validate_region_indices(
    mesh: &ConductionMesh,
    junction: &JunctionRegion,
    surface: &SurfaceRegion,
) -> Result<(), QoiError> {
    if let Some(&vertex) = junction
        .vertices()
        .iter()
        .find(|&&vertex| vertex >= mesh.vertex_count())
    {
        return Err(QoiError::invalid(
            "junction region",
            format!(
                "vertex {vertex} is outside {} vertices",
                mesh.vertex_count()
            ),
        ));
    }
    if let Some(&face) = surface
        .boundary_faces()
        .iter()
        .find(|&&face| face >= mesh.boundary().len())
    {
        return Err(QoiError::invalid(
            "surface region",
            format!(
                "boundary-face slot {face} is outside {} boundary faces",
                mesh.boundary().len()
            ),
        ));
    }
    Ok(())
}

fn junction_maximum(
    solution: &ConductionSolution,
    region: &JunctionRegion,
) -> (Temperature, usize) {
    // `JunctionRegion` admission sorts vertices ascending and rejects
    // duplicates, so every candidate is strictly greater than the incumbent.
    // The strict `>` is therefore the whole tie-break: an equal value never
    // displaces the incumbent, so the lowest canonical index wins. A
    // `candidate < vertex` guard would be unreachable here, not a safeguard.
    // `validate_solution` has already rejected non-finite nodal values, so no
    // NaN can reach this comparison.
    let mut vertex = region.vertices()[0];
    let mut maximum = solution.temperature[vertex];
    for &candidate in &region.vertices()[1..] {
        let value = solution.temperature[candidate];
        if value > maximum {
            maximum = value;
            vertex = candidate;
        }
    }
    (Temperature::new(maximum), vertex)
}

struct SurfaceStatistics {
    mean: f64,
    spread: f64,
    standard_deviation: f64,
}

fn surface_statistics(
    mesh: &ConductionMesh,
    solution: &ConductionSolution,
    region: &SurfaceRegion,
) -> Result<SurfaceStatistics, QoiError> {
    let mut area = 0.0;
    let mut integral = 0.0;
    let mut face_means = Vec::with_capacity(region.boundary_faces().len());
    let mut vertices = BTreeSet::new();
    for &slot in region.boundary_faces() {
        let face = &mesh.boundary()[slot];
        let mean = face
            .vertices
            .iter()
            .map(|&vertex| solution.temperature[vertex as usize])
            .sum::<f64>()
            / 3.0;
        area += face.area;
        integral = face.area.mul_add(mean, integral);
        face_means.push((face.area, mean));
        vertices.extend(face.vertices.iter().map(|&vertex| vertex as usize));
    }
    if !(area.is_finite() && area > 0.0 && integral.is_finite()) {
        return Err(QoiError::invalid(
            "surface region",
            "selected area and temperature integral must be finite with positive area",
        ));
    }
    let mean = integral / area;
    let mut variance_integral = 0.0;
    for (face_area, face_mean) in face_means {
        let delta = face_mean - mean;
        variance_integral = face_area.mul_add(delta * delta, variance_integral);
    }
    let variance = (variance_integral / area).max(0.0);
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for vertex in vertices {
        let temperature = solution.temperature[vertex];
        minimum = minimum.min(temperature);
        maximum = maximum.max(temperature);
    }
    let spread = maximum - minimum;
    let standard_deviation = fs_math::det::sqrt(variance);
    if !(mean.is_finite() && spread.is_finite() && standard_deviation.is_finite()) {
        return Err(QoiError::invalid(
            "surface statistics",
            "mean, spread, or standard deviation was non-finite",
        ));
    }
    Ok(SurfaceStatistics {
        mean,
        spread,
        standard_deviation,
    })
}

fn interval_half_width(value: f64, lo: f64, hi: f64) -> Result<f64, QoiError> {
    validate_estimate("QoI envelope", value, lo, hi)?;
    let half_width = (value - lo).abs().max((hi - value).abs());
    if half_width.is_finite() {
        Ok(half_width)
    } else {
        Err(QoiError::invalid(
            "QoI envelope",
            "half-width arithmetic overflowed",
        ))
    }
}

fn admit_name(field: &'static str, value: String) -> Result<String, QoiError> {
    if value.trim().is_empty() || value.len() > 1024 {
        Err(QoiError::invalid(field, "name must contain 1..=1024 bytes"))
    } else {
        Ok(value)
    }
}

fn canonicalize_indices(field: &'static str, values: &mut Vec<usize>) -> Result<(), QoiError> {
    if values.is_empty() || values.len() > MAX_REGION_ENTRIES {
        return Err(QoiError::invalid(
            field,
            format!("entry count must lie in 1..={MAX_REGION_ENTRIES}"),
        ));
    }
    values.sort_unstable();
    if values.windows(2).any(|pair| pair[0] == pair[1]) {
        return Err(QoiError::invalid(
            field,
            "duplicate indices are not admissible",
        ));
    }
    Ok(())
}

fn solution_identity(mesh: &ConductionMesh, solution: &ConductionSolution) -> ContentHash {
    let mut bytes = Vec::new();
    push_usize(&mut bytes, mesh.vertex_count());
    push_usize(&mut bytes, mesh.element_count());
    for position in mesh.positions() {
        for coordinate in position {
            bytes.extend_from_slice(&coordinate.to_bits().to_le_bytes());
        }
    }
    for tet in &mesh.complex().tets {
        for vertex in tet {
            bytes.extend_from_slice(&vertex.to_le_bytes());
        }
    }
    for value in &solution.temperature {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(&solution.report.final_residual.to_bits().to_le_bytes());
    bytes.extend_from_slice(&solution.report.residual_threshold.to_bits().to_le_bytes());
    bytes.extend_from_slice(&solution.report.energy.closure_w.to_bits().to_le_bytes());
    bytes.extend_from_slice(solution.report.material_provenance.tag().as_bytes());
    push_usize(&mut bytes, solution.report.material_receipts);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn operating_identity(point: &OperatingPoint) -> ContentHash {
    let mut bytes = Vec::new();
    for value in [
        point.flow.value.value(),
        point.flow.numerical.lo,
        point.flow.numerical.hi,
        point.pressure.value.value(),
        point.pressure.numerical.lo,
        point.pressure.numerical.hi,
        point.nominal_root.flow.lo(),
        point.nominal_root.flow.hi(),
        point.leakage_fraction,
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    bytes.extend_from_slice(&point.flow.provenance.0.to_le_bytes());
    bytes.extend_from_slice(&point.pressure.provenance.0.to_le_bytes());
    for card in &point.pressure.model.cards {
        push_string(&mut bytes, card);
    }
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn fan_power_identity(spec: &FanPowerSpec) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&spec.total_efficiency.to_bits().to_le_bytes());
    bytes.extend_from_slice(&spec.efficiency_half_width.to_bits().to_le_bytes());
    push_string(&mut bytes, &spec.source.citation);
    push_string(&mut bytes, &spec.source.identifier);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn requirement_identity(requirement: &ThermalRequirement) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(
        &requirement
            .maximum_temperature
            .value()
            .to_bits()
            .to_le_bytes(),
    );
    push_string(&mut bytes, &requirement.source.citation);
    push_string(&mut bytes, &requirement.source.identifier);
    // The already-applied factor and its policy bind into identity so two runs
    // that share an effective limit but derated it under different authority
    // are distinguishable. The factor never reaches the margin arithmetic.
    bytes.extend_from_slice(&requirement.safety_factor.factor.to_bits().to_le_bytes());
    push_string(&mut bytes, &requirement.safety_factor.source.citation);
    push_string(&mut bytes, &requirement.safety_factor.source.identifier);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn discretization_identity(receipt: &DiscretizationReceipt) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&receipt.half_width.to_bits().to_le_bytes());
    push_string(&mut bytes, &receipt.source.citation);
    push_string(&mut bytes, &receipt.source.identifier);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn region_identity(name: &str, indices: &[usize]) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_string(&mut bytes, name);
    for &index in indices {
        push_usize(&mut bytes, index);
    }
    bytes
}

fn qoi_identity(label: &str, parents: &[ContentHash], context: &[u8]) -> ContentHash {
    let mut bytes = Vec::new();
    push_string(&mut bytes, label);
    for parent in parents {
        bytes.extend_from_slice(parent.as_bytes());
    }
    bytes.extend_from_slice(context);
    hash_domain(QOI_IDENTITY_DOMAIN, &bytes)
}

fn term_identity(qoi: ContentHash, kind: EngineeringUncertaintyKind) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(qoi.as_bytes());
    bytes.extend_from_slice(kind.name().as_bytes());
    hash_domain(QOI_TERM_DOMAIN, &bytes)
}

fn push_string(bytes: &mut Vec<u8>, value: &str) {
    push_usize(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_usize(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(&(value as u64).to_le_bytes());
}
