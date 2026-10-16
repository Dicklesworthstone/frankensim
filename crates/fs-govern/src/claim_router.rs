//! Executable claim-to-evidence routing for the E09 certificate doctrine.
//!
//! Routing is a governance decision, not evidence. A successful route names the
//! machinery and assumptions a claim would need; it does not execute that
//! machinery, authenticate an artifact, or mint scientific authority.

use core::fmt;

use crate::{
    CERTIFICATE_REGIME_NO_CLAIM, CertificateRegimeRow, ClaimClass, EvidenceRegime,
    certificate_regime,
};

/// Schema version for typed claim requests and routing decisions.
pub const CLAIM_ROUTER_SCHEMA_VERSION: u16 = 1;

/// Maximum UTF-8 bytes in one caller-supplied routing field.
pub const MAX_CLAIM_ROUTER_FIELD_BYTES: usize = 512;

/// Maximum explicit assumptions retained by one claim request.
pub const MAX_CLAIM_ROUTER_ASSUMPTIONS: usize = 64;

/// Stable no-authority boundary for every route and refusal.
pub const CLAIM_ROUTER_NO_CLAIM: &str = "routing selects required evidence machinery and records assumptions; it does not execute a solver, mint evidence, authenticate artifacts, establish scientific truth, or admit a runtime claim";

fn valid_field(value: &str) -> bool {
    !value.trim().is_empty()
        && value.len() <= MAX_CLAIM_ROUTER_FIELD_BYTES
        && !value.chars().any(char::is_control)
}

fn validate_field(field: &'static str, value: &str) -> Result<(), ClaimRouterError> {
    if valid_field(value) {
        Ok(())
    } else {
        Err(ClaimRouterError::InvalidField {
            field,
            bytes: value.len(),
        })
    }
}

/// Decision-relevant extent of a requested claim.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimExtent {
    /// A local statement at one declared state, equilibrium, or orbit.
    Local,
    /// A finite time or parameter horizon.
    FiniteHorizon {
        /// Positive finite horizon length.
        duration: f64,
        /// Horizon unit.
        unit: String,
    },
    /// A finite non-temporal parameter domain.
    FiniteParameterDomain {
        /// Positive finite parameter-domain span.
        span: f64,
        /// Parameter unit.
        unit: String,
    },
    /// A declared long-run horizon.
    LongHorizon {
        /// Positive finite requested horizon length.
        duration: f64,
        /// Horizon unit.
        unit: String,
    },
    /// A population or ensemble rather than one pointwise trajectory.
    Population {
        /// Declared population.
        population: String,
    },
}

impl ClaimExtent {
    /// Construct a finite-horizon extent.
    pub fn try_finite_horizon(
        duration: f64,
        unit: impl Into<String>,
    ) -> Result<Self, ClaimRouterError> {
        let unit = unit.into();
        validate_duration("finite_horizon.duration", duration)?;
        validate_field("finite_horizon.unit", &unit)?;
        Ok(Self::FiniteHorizon { duration, unit })
    }

    /// Construct a long-horizon extent.
    pub fn try_long_horizon(
        duration: f64,
        unit: impl Into<String>,
    ) -> Result<Self, ClaimRouterError> {
        let unit = unit.into();
        validate_duration("long_horizon.duration", duration)?;
        validate_field("long_horizon.unit", &unit)?;
        Ok(Self::LongHorizon { duration, unit })
    }

    /// Construct a finite parameter-domain extent.
    pub fn try_finite_parameter_domain(
        span: f64,
        unit: impl Into<String>,
    ) -> Result<Self, ClaimRouterError> {
        let unit = unit.into();
        validate_duration("finite_parameter_domain.span", span)?;
        validate_field("finite_parameter_domain.unit", &unit)?;
        Ok(Self::FiniteParameterDomain { span, unit })
    }

    /// Construct a population extent.
    pub fn try_population(population: impl Into<String>) -> Result<Self, ClaimRouterError> {
        let population = population.into();
        validate_field("population", &population)?;
        Ok(Self::Population { population })
    }

    /// Stable extent code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::Local => "local",
            Self::FiniteHorizon { .. } => "finite-horizon",
            Self::FiniteParameterDomain { .. } => "finite-parameter-domain",
            Self::LongHorizon { .. } => "long-horizon",
            Self::Population { .. } => "population",
        }
    }

    fn temporal_duration(&self) -> Option<(f64, &str)> {
        match self {
            Self::FiniteHorizon { duration, unit } | Self::LongHorizon { duration, unit } => {
                Some((*duration, unit))
            }
            Self::Local | Self::FiniteParameterDomain { .. } | Self::Population { .. } => None,
        }
    }

    fn validate(&self) -> Result<(), ClaimRouterError> {
        match self {
            Self::Local => Ok(()),
            Self::FiniteHorizon { duration, unit } => {
                validate_duration("finite_horizon.duration", *duration)?;
                validate_field("finite_horizon.unit", unit)
            }
            Self::FiniteParameterDomain { span, unit } => {
                validate_duration("finite_parameter_domain.span", *span)?;
                validate_field("finite_parameter_domain.unit", unit)
            }
            Self::LongHorizon { duration, unit } => {
                validate_duration("long_horizon.duration", *duration)?;
                validate_field("long_horizon.unit", unit)
            }
            Self::Population { population } => validate_field("population", population),
        }
    }
}

fn validate_duration(field: &'static str, duration: f64) -> Result<(), ClaimRouterError> {
    if duration.is_finite() && duration > 0.0 {
        Ok(())
    } else {
        Err(ClaimRouterError::InvalidPositiveFinite {
            field,
            found: duration,
        })
    }
}

/// Basis for treating a system as chaotic in a routing request.
#[derive(Debug, Clone, PartialEq)]
pub enum ChaosBasis {
    /// No chaotic-system assumption or probe result is asserted.
    NotIndicated,
    /// Caller-declared chaotic behavior with a declared predictability horizon.
    Declared {
        /// Positive finite predictability horizon.
        predictability_horizon: f64,
        /// Horizon unit.
        unit: String,
    },
    /// A positive local-Lyapunov probe plus a declared predictability horizon.
    Probed {
        /// Positive finite lower estimate for the local Lyapunov exponent.
        local_lyapunov_lower: f64,
        /// Positive finite predictability horizon.
        predictability_horizon: f64,
        /// Horizon unit.
        unit: String,
    },
}

impl ChaosBasis {
    /// No chaotic-system classification.
    #[must_use]
    pub const fn not_indicated() -> Self {
        Self::NotIndicated
    }

    /// Construct a declared chaotic-system assumption.
    pub fn try_declared(
        predictability_horizon: f64,
        unit: impl Into<String>,
    ) -> Result<Self, ClaimRouterError> {
        let unit = unit.into();
        validate_duration("chaos.predictability_horizon", predictability_horizon)?;
        validate_field("chaos.unit", &unit)?;
        Ok(Self::Declared {
            predictability_horizon,
            unit,
        })
    }

    /// Construct a probe-backed chaotic-system assumption.
    pub fn try_probed(
        local_lyapunov_lower: f64,
        predictability_horizon: f64,
        unit: impl Into<String>,
    ) -> Result<Self, ClaimRouterError> {
        let unit = unit.into();
        validate_duration("chaos.local_lyapunov_lower", local_lyapunov_lower)?;
        validate_duration("chaos.predictability_horizon", predictability_horizon)?;
        validate_field("chaos.unit", &unit)?;
        Ok(Self::Probed {
            local_lyapunov_lower,
            predictability_horizon,
            unit,
        })
    }

    /// Whether the request explicitly treats the system as chaotic.
    #[must_use]
    pub const fn is_chaotic(&self) -> bool {
        !matches!(self, Self::NotIndicated)
    }

    fn predictability_horizon(&self) -> Option<(f64, &str)> {
        match self {
            Self::NotIndicated => None,
            Self::Declared {
                predictability_horizon,
                unit,
            }
            | Self::Probed {
                predictability_horizon,
                unit,
                ..
            } => Some((*predictability_horizon, unit)),
        }
    }

    fn code(&self) -> &'static str {
        match self {
            Self::NotIndicated => "not-indicated",
            Self::Declared { .. } => "declared-chaotic",
            Self::Probed { .. } => "local-lyapunov-probed",
        }
    }

    fn validate(&self) -> Result<(), ClaimRouterError> {
        match self {
            Self::NotIndicated => Ok(()),
            Self::Declared {
                predictability_horizon,
                unit,
            } => {
                validate_duration("chaos.predictability_horizon", *predictability_horizon)?;
                validate_field("chaos.unit", unit)
            }
            Self::Probed {
                local_lyapunov_lower,
                predictability_horizon,
                unit,
            } => {
                validate_duration("chaos.local_lyapunov_lower", *local_lyapunov_lower)?;
                validate_duration("chaos.predictability_horizon", *predictability_horizon)?;
                validate_field("chaos.unit", unit)
            }
        }
    }
}

/// Declared dynamics properties used by the router.
///
/// The Boolean characteristics are deliberately non-exclusive: a coupled model
/// can contain conservative, dissipative, and stiff subsystems at once.
#[derive(Debug, Clone, PartialEq)]
pub struct DynamicsProfile {
    conservative: bool,
    dissipative: bool,
    stiff: bool,
    chaos: ChaosBasis,
}

impl DynamicsProfile {
    /// Construct a profile from explicit characteristics.
    #[must_use]
    pub const fn new(
        conservative: bool,
        dissipative: bool,
        stiff: bool,
        chaos: ChaosBasis,
    ) -> Self {
        Self {
            conservative,
            dissipative,
            stiff,
            chaos,
        }
    }

    /// A system with no stronger classification asserted.
    #[must_use]
    pub const fn general() -> Self {
        Self::new(false, false, false, ChaosBasis::NotIndicated)
    }

    /// Whether a conservative characteristic is declared.
    #[must_use]
    pub const fn conservative(&self) -> bool {
        self.conservative
    }

    /// Whether a dissipative characteristic is declared.
    #[must_use]
    pub const fn dissipative(&self) -> bool {
        self.dissipative
    }

    /// Whether a stiff characteristic is declared.
    #[must_use]
    pub const fn stiff(&self) -> bool {
        self.stiff
    }

    /// Chaotic-system assumption or probe basis.
    #[must_use]
    pub const fn chaos(&self) -> &ChaosBasis {
        &self.chaos
    }
}

/// Caller-declared decision need.
#[derive(Debug, Clone, PartialEq)]
pub struct DecisionNeed {
    context: String,
    unit: String,
    tolerance: f64,
}

impl DecisionNeed {
    /// Construct a positive finite decision tolerance.
    pub fn try_new(
        context: impl Into<String>,
        unit: impl Into<String>,
        tolerance: f64,
    ) -> Result<Self, ClaimRouterError> {
        let context = context.into();
        let unit = unit.into();
        validate_field("decision.context", &context)?;
        validate_field("decision.unit", &unit)?;
        validate_duration("decision.tolerance", tolerance)?;
        Ok(Self {
            context,
            unit,
            tolerance,
        })
    }

    /// Decision context.
    #[must_use]
    pub fn context(&self) -> &str {
        &self.context
    }

    /// Tolerance unit.
    #[must_use]
    pub fn unit(&self) -> &str {
        &self.unit
    }

    /// Positive decision tolerance.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

/// Typed request presented to the claim router.
#[derive(Debug, Clone, PartialEq)]
pub struct ClaimRequest {
    request_id: String,
    claim: ClaimClass,
    quantity: String,
    extent: ClaimExtent,
    decision: DecisionNeed,
    system: DynamicsProfile,
    assumptions: Vec<String>,
}

impl ClaimRequest {
    /// Construct a claim request with canonical, explicit assumptions.
    pub fn try_new(
        request_id: impl Into<String>,
        claim: ClaimClass,
        quantity: impl Into<String>,
        extent: ClaimExtent,
        decision: DecisionNeed,
        system: DynamicsProfile,
        mut assumptions: Vec<String>,
    ) -> Result<Self, ClaimRouterError> {
        let request_id = request_id.into();
        let quantity = quantity.into();
        validate_field("request_id", &request_id)?;
        validate_field("quantity", &quantity)?;
        extent.validate()?;
        system.chaos().validate()?;
        if assumptions.len() > MAX_CLAIM_ROUTER_ASSUMPTIONS {
            return Err(ClaimRouterError::TooManyAssumptions {
                found: assumptions.len(),
            });
        }
        for assumption in &assumptions {
            validate_field("assumption", assumption)?;
        }
        assumptions.sort_unstable();
        assumptions.dedup();
        if assumptions.is_empty() {
            return Err(ClaimRouterError::MissingAssumptions);
        }
        Ok(Self {
            request_id,
            claim,
            quantity,
            extent,
            decision,
            system,
            assumptions,
        })
    }

    /// Stable request id.
    #[must_use]
    pub fn request_id(&self) -> &str {
        &self.request_id
    }

    /// Requested claim class.
    #[must_use]
    pub const fn claim(&self) -> ClaimClass {
        self.claim
    }

    /// Named quantity.
    #[must_use]
    pub fn quantity(&self) -> &str {
        &self.quantity
    }

    /// Decision-relevant extent.
    #[must_use]
    pub const fn extent(&self) -> &ClaimExtent {
        &self.extent
    }

    /// Caller-declared decision need.
    #[must_use]
    pub const fn decision(&self) -> &DecisionNeed {
        &self.decision
    }

    /// Dynamics classification and its chaos basis.
    #[must_use]
    pub const fn system(&self) -> &DynamicsProfile {
        &self.system
    }

    /// Canonical explicit assumptions.
    #[must_use]
    pub fn assumptions(&self) -> &[String] {
        &self.assumptions
    }
}

/// Why a valid request was refused before computation.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimRouteRefusalCause {
    /// The claim class and extent do not match the doctrine.
    ExtentMismatch {
        /// Extent required by this claim class.
        required: &'static str,
        /// Extent presented by the caller.
        found: &'static str,
    },
    /// A chaotic finite-horizon request exceeds its own predictability horizon.
    PredictabilityHorizonExceeded {
        /// Requested horizon.
        requested: f64,
        /// Declared predictability horizon.
        predictability_horizon: f64,
        /// Shared horizon unit.
        unit: String,
    },
    /// An exact-chaotic-trajectory request supplied no chaotic classification.
    MissingChaoticClassification,
    /// An exact-chaotic-trajectory request is still inside its predictability horizon.
    InsidePredictabilityHorizon {
        /// Requested horizon.
        requested: f64,
        /// Declared predictability horizon.
        predictability_horizon: f64,
        /// Shared horizon unit.
        unit: String,
    },
    /// The requested exact trajectory is beyond the admitted predictability horizon.
    ExactLongChaoticTrajectoryHasNoUsefulRoute {
        /// Requested horizon.
        requested: f64,
        /// Declared predictability horizon.
        predictability_horizon: f64,
        /// Shared horizon unit.
        unit: String,
    },
    /// Horizon units cannot be compared.
    IncompatibleHorizonUnits {
        /// Request extent unit.
        requested_unit: String,
        /// Predictability-horizon unit.
        predictability_unit: String,
    },
}

impl ClaimRouteRefusalCause {
    /// Stable machine code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::ExtentMismatch { .. } => "extent-mismatch",
            Self::PredictabilityHorizonExceeded { .. } => "predictability-horizon-exceeded",
            Self::MissingChaoticClassification => "missing-chaotic-classification",
            Self::InsidePredictabilityHorizon { .. } => "inside-predictability-horizon",
            Self::ExactLongChaoticTrajectoryHasNoUsefulRoute { .. } => {
                "exact-long-chaotic-trajectory-no-useful-route"
            }
            Self::IncompatibleHorizonUnits { .. } => "incompatible-horizon-units",
        }
    }
}

/// Successful doctrine route.
#[derive(Debug, Clone, PartialEq)]
pub struct RoutedClaim {
    request: ClaimRequest,
    row_id: &'static str,
    evidence: EvidenceRegime,
}

impl RoutedClaim {
    /// Exact request retained as routing provenance.
    #[must_use]
    pub const fn request(&self) -> &ClaimRequest {
        &self.request
    }

    /// Canonical doctrine row id.
    #[must_use]
    pub const fn row_id(&self) -> &'static str {
        self.row_id
    }

    /// Required evidence regime.
    #[must_use]
    pub const fn evidence(&self) -> EvidenceRegime {
        self.evidence
    }
}

/// Request-time refusal with doctrine reasoning and a reformulation.
#[derive(Debug, Clone, PartialEq)]
pub struct ClaimRouteRefusal {
    request: ClaimRequest,
    row_id: &'static str,
    required_evidence: EvidenceRegime,
    cause: ClaimRouteRefusalCause,
    suggested_reformulation: ClaimClass,
}

impl ClaimRouteRefusal {
    /// Exact request retained as routing provenance.
    #[must_use]
    pub const fn request(&self) -> &ClaimRequest {
        &self.request
    }

    /// Canonical doctrine row id.
    #[must_use]
    pub const fn row_id(&self) -> &'static str {
        self.row_id
    }

    /// Evidence regime the requested claim would have required.
    #[must_use]
    pub const fn required_evidence(&self) -> EvidenceRegime {
        self.required_evidence
    }

    /// Typed refusal cause.
    #[must_use]
    pub const fn cause(&self) -> &ClaimRouteRefusalCause {
        &self.cause
    }

    /// Doctrine-supported claim reformulation.
    #[must_use]
    pub const fn suggested_reformulation(&self) -> ClaimClass {
        self.suggested_reformulation
    }
}

/// Executable result of routing one valid claim request.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimRouteDecision {
    /// The doctrine selected required evidence machinery.
    Routed(RoutedClaim),
    /// The request was refused before computation.
    Refused(ClaimRouteRefusal),
}

impl ClaimRouteDecision {
    /// Exact request retained in either outcome.
    #[must_use]
    pub const fn request(&self) -> &ClaimRequest {
        match self {
            Self::Routed(route) => route.request(),
            Self::Refused(refusal) => refusal.request(),
        }
    }

    /// Canonical doctrine row used by the decision.
    #[must_use]
    pub fn doctrine_row(&self) -> &'static CertificateRegimeRow {
        certificate_regime(self.request().claim())
    }

    /// Successful route, if any.
    #[must_use]
    pub const fn routed(&self) -> Option<&RoutedClaim> {
        match self {
            Self::Routed(route) => Some(route),
            Self::Refused(_) => None,
        }
    }

    /// Request-time refusal, if any.
    #[must_use]
    pub const fn refusal(&self) -> Option<&ClaimRouteRefusal> {
        match self {
            Self::Routed(_) => None,
            Self::Refused(refusal) => Some(refusal),
        }
    }

    /// Deterministic line-oriented routing provenance.
    ///
    /// This rendering is inspectable governance data, not a durable identity
    /// transport. A future schema owner must register any content-addressed
    /// encoding separately.
    #[must_use]
    pub fn render_record(&self) -> String {
        let request = self.request();
        let row = self.doctrine_row();
        let mut output = format!(
            "claim-router-schema={CLAIM_ROUTER_SCHEMA_VERSION}\nrequest-id={}\noutcome={}\nclaim={}\nquantity={}\nextent={}\ndecision-context={}\ndecision-tolerance={}\ndecision-unit={}\nsystem-conservative={}\nsystem-dissipative={}\nsystem-stiff={}\nchaos-basis={}\ndoctrine-row={}\nrequired-evidence={}\n",
            request.request_id(),
            if self.routed().is_some() {
                "routed"
            } else {
                "refused"
            },
            request.claim().code(),
            request.quantity(),
            request.extent().code(),
            request.decision().context(),
            request.decision().tolerance(),
            request.decision().unit(),
            request.system().conservative(),
            request.system().dissipative(),
            request.system().stiff(),
            request.system().chaos().code(),
            row.id,
            row.evidence.code(),
        );
        match request.extent() {
            ClaimExtent::Local => {}
            ClaimExtent::FiniteHorizon { duration, unit }
            | ClaimExtent::LongHorizon { duration, unit } => {
                output.push_str(&format!("extent-duration={duration}\nextent-unit={unit}\n"));
            }
            ClaimExtent::FiniteParameterDomain { span, unit } => {
                output.push_str(&format!(
                    "extent-span={span}\nextent-unit={unit}\nextent-axis=parameter\n"
                ));
            }
            ClaimExtent::Population { population } => {
                output.push_str(&format!("population={population}\n"));
            }
        }
        match request.system().chaos() {
            ChaosBasis::NotIndicated => {}
            ChaosBasis::Declared {
                predictability_horizon,
                unit,
            } => {
                output.push_str(&format!(
                    "predictability-horizon={predictability_horizon}\npredictability-unit={unit}\n"
                ));
            }
            ChaosBasis::Probed {
                local_lyapunov_lower,
                predictability_horizon,
                unit,
            } => {
                output.push_str(&format!(
                    "local-lyapunov-lower={local_lyapunov_lower}\npredictability-horizon={predictability_horizon}\npredictability-unit={unit}\n"
                ));
            }
        }
        if let Some(refusal) = self.refusal() {
            output.push_str(&format!(
                "refusal-cause={}\nsuggested-reformulation={}\n",
                refusal.cause().code(),
                refusal.suggested_reformulation().code(),
            ));
        }
        for capability in row.capabilities {
            output.push_str(&format!(
                "capability={}/{}:{}\n",
                capability.crate_name,
                capability.capability,
                capability.status.code(),
            ));
        }
        for assumption in request.assumptions() {
            output.push_str(&format!("assumption={assumption}\n"));
        }
        output.push_str(&format!(
            "row-no-claim={}\nrouter-no-claim={CLAIM_ROUTER_NO_CLAIM}\ndoctrine-no-claim={CERTIFICATE_REGIME_NO_CLAIM}\n",
            row.no_claim,
        ));
        output
    }
}

fn routed(request: ClaimRequest, row: &'static CertificateRegimeRow) -> ClaimRouteDecision {
    ClaimRouteDecision::Routed(RoutedClaim {
        request,
        row_id: row.id,
        evidence: row.evidence,
    })
}

fn refused(
    request: ClaimRequest,
    row: &'static CertificateRegimeRow,
    cause: ClaimRouteRefusalCause,
    suggested_reformulation: ClaimClass,
) -> ClaimRouteDecision {
    ClaimRouteDecision::Refused(ClaimRouteRefusal {
        request,
        row_id: row.id,
        required_evidence: row.evidence,
        cause,
        suggested_reformulation,
    })
}

fn extent_mismatch(
    request: ClaimRequest,
    row: &'static CertificateRegimeRow,
    required: &'static str,
    suggested_reformulation: ClaimClass,
) -> ClaimRouteDecision {
    let found = request.extent().code();
    refused(
        request,
        row,
        ClaimRouteRefusalCause::ExtentMismatch { required, found },
        suggested_reformulation,
    )
}

fn compare_to_predictability_horizon(
    request: &ClaimRequest,
) -> Result<Option<(f64, f64, String)>, ClaimRouteRefusalCause> {
    let Some((requested, requested_unit)) = request.extent().temporal_duration() else {
        return Ok(None);
    };
    let Some((predictability_horizon, predictability_unit)) =
        request.system().chaos().predictability_horizon()
    else {
        return Ok(None);
    };
    if requested_unit != predictability_unit {
        return Err(ClaimRouteRefusalCause::IncompatibleHorizonUnits {
            requested_unit: requested_unit.to_string(),
            predictability_unit: predictability_unit.to_string(),
        });
    }
    Ok(Some((
        requested,
        predictability_horizon,
        requested_unit.to_string(),
    )))
}

/// Route one valid request through the closed E09 doctrine.
///
/// Decision tolerance is retained as context but never used to widen the set of
/// admissible claim/evidence pairs. Relaxing it therefore cannot turn a route
/// into a refusal.
#[must_use]
pub fn route_claim(request: ClaimRequest) -> ClaimRouteDecision {
    let row = certificate_regime(request.claim());
    match request.claim() {
        ClaimClass::RootOrEventTime => {
            if !matches!(
                request.extent(),
                ClaimExtent::FiniteHorizon { .. } | ClaimExtent::FiniteParameterDomain { .. }
            ) {
                return extent_mismatch(
                    request,
                    row,
                    "finite-time-or-parameter-domain",
                    ClaimClass::LongHorizonMeanLoad,
                );
            }
            match compare_to_predictability_horizon(&request) {
                Err(cause) => refused(request, row, cause, ClaimClass::LongHorizonMeanLoad),
                Ok(Some((requested, predictability_horizon, unit)))
                    if requested > predictability_horizon =>
                {
                    refused(
                        request,
                        row,
                        ClaimRouteRefusalCause::PredictabilityHorizonExceeded {
                            requested,
                            predictability_horizon,
                            unit,
                        },
                        ClaimClass::LongHorizonMeanLoad,
                    )
                }
                Ok(_) => routed(request, row),
            }
        }
        ClaimClass::ShortHorizonReachability => {
            if !matches!(request.extent(), ClaimExtent::FiniteHorizon { .. }) {
                return extent_mismatch(
                    request,
                    row,
                    "finite-horizon",
                    ClaimClass::LongHorizonMeanLoad,
                );
            }
            match compare_to_predictability_horizon(&request) {
                Err(cause) => refused(request, row, cause, ClaimClass::LongHorizonMeanLoad),
                Ok(Some((requested, predictability_horizon, unit)))
                    if requested > predictability_horizon =>
                {
                    refused(
                        request,
                        row,
                        ClaimRouteRefusalCause::PredictabilityHorizonExceeded {
                            requested,
                            predictability_horizon,
                            unit,
                        },
                        ClaimClass::LongHorizonMeanLoad,
                    )
                }
                Ok(_) => routed(request, row),
            }
        }
        ClaimClass::ConservedQuantity => {
            if matches!(request.extent(), ClaimExtent::FiniteHorizon { .. }) {
                routed(request, row)
            } else {
                extent_mismatch(
                    request,
                    row,
                    "finite-horizon",
                    ClaimClass::ConservedQuantity,
                )
            }
        }
        ClaimClass::LocalStability => {
            if matches!(request.extent(), ClaimExtent::Local) {
                routed(request, row)
            } else {
                extent_mismatch(request, row, "local", ClaimClass::LocalStability)
            }
        }
        ClaimClass::LongHorizonMeanLoad | ClaimClass::BroadbandSpectrum => {
            if matches!(request.extent(), ClaimExtent::LongHorizon { .. }) {
                routed(request, row)
            } else {
                extent_mismatch(
                    request,
                    row,
                    "long-horizon",
                    ClaimClass::LongHorizonMeanLoad,
                )
            }
        }
        ClaimClass::DutyCycleReliability => {
            if matches!(request.extent(), ClaimExtent::Population { .. }) {
                routed(request, row)
            } else {
                extent_mismatch(request, row, "population", ClaimClass::DutyCycleReliability)
            }
        }
        ClaimClass::ExactLongChaoticTrajectory => {
            if !matches!(request.extent(), ClaimExtent::LongHorizon { .. }) {
                return extent_mismatch(
                    request,
                    row,
                    "long-horizon-beyond-predictability",
                    ClaimClass::RootOrEventTime,
                );
            }
            if !request.system().chaos().is_chaotic() {
                return refused(
                    request,
                    row,
                    ClaimRouteRefusalCause::MissingChaoticClassification,
                    ClaimClass::LocalStability,
                );
            }
            match compare_to_predictability_horizon(&request) {
                Err(cause) => refused(request, row, cause, ClaimClass::LongHorizonMeanLoad),
                Ok(Some((requested, predictability_horizon, unit)))
                    if requested > predictability_horizon =>
                {
                    refused(
                        request,
                        row,
                        ClaimRouteRefusalCause::ExactLongChaoticTrajectoryHasNoUsefulRoute {
                            requested,
                            predictability_horizon,
                            unit,
                        },
                        ClaimClass::LongHorizonMeanLoad,
                    )
                }
                Ok(Some((requested, predictability_horizon, unit))) => refused(
                    request,
                    row,
                    ClaimRouteRefusalCause::InsidePredictabilityHorizon {
                        requested,
                        predictability_horizon,
                        unit,
                    },
                    ClaimClass::RootOrEventTime,
                ),
                Ok(None) => refused(
                    request,
                    row,
                    ClaimRouteRefusalCause::MissingChaoticClassification,
                    ClaimClass::LocalStability,
                ),
            }
        }
    }
}

/// Malformed claim-request input.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimRouterError {
    /// A bounded text field was empty, oversized, or contained control bytes.
    InvalidField {
        /// Rejected field.
        field: &'static str,
        /// Presented UTF-8 byte count.
        bytes: usize,
    },
    /// A numeric policy field was not positive and finite.
    InvalidPositiveFinite {
        /// Rejected field.
        field: &'static str,
        /// Presented value.
        found: f64,
    },
    /// More assumptions were presented than the schema permits.
    TooManyAssumptions {
        /// Presented count.
        found: usize,
    },
    /// No explicit model/system assumption was retained.
    MissingAssumptions,
}

impl fmt::Display for ClaimRouterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for ClaimRouterError {}
