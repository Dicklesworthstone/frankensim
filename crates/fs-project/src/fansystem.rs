//! Versioned Cooling fan-system declaration and fail-closed lowering to
//! production `fs-airflow` evidence (bead frn2i.1).
//!
//! A [`FanSystemDecl`] carries every parameter the solver consumes: stable
//! bank identities, explicit counts, series/parallel arrangements, declared
//! speed ratios inside declared domains, monotone curves with stall
//! boundaries, pressure tolerance and basis, source citation/ID, and an
//! explicit system topology whenever more than one bank is declared. The
//! current rated point is either checked against the declared curve within
//! the declared tolerance or typed as correlation-only and excluded from
//! the network solve.
//!
//! Nothing is inferred: no first-fan choice, no implicit count or speed,
//! no midpoint, no identical-fan collapse. [`lower_fan_system`] consumes
//! every declared bank exactly once and constructs the production fan
//! system; zero/duplicate/orphan banks, ambiguous topology, and unsupported
//! composition refuse with structured [`ProjectError`]s.

use fs_airflow::{
    FanArrangement, FanBank, FanCurve, FanPoint, SourceProvenance, ToleranceBasis,
    composite::{compose_parallel, compose_series},
};
use fs_qty::{Dims, QtyAny};

use crate::spec::{FanCurveDecl, FanToleranceBasis};
use crate::wire::ProjectError;

/// The fan-system declaration schema version.
pub const FAN_SYSTEM_DECL_VERSION: u32 = 1;

/// Wire prefix of the fan-system declaration identity.
pub const FAN_SYSTEM_IDENTITY_PREFIX: &str = "fan-system-decl:v1:";

const IDENTITY_DOMAIN: &str = "org.frankensim.fs-project.fan-system-decl.v1";

const FLOW_DIMS: Dims = Dims([3, 0, -1, 0, 0, 0]);
const PRESSURE_DIMS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn fan_error(code: &'static str, detail: String, hint: impl Into<String>) -> ProjectError {
    ProjectError {
        code,
        detail,
        hint: hint.into(),
    }
}

/// How the rated (nominal) operating point participates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RatedPointAdmission {
    /// The rated point was checked against the declared curve within the
    /// declared relative pressure tolerance.
    CheckedWithinDeclaredTolerance,
    /// Correlation-rung evidence only: excluded from the network solve.
    CorrelationOnly,
}

/// A declared rated operating point with its admission classification.
#[derive(Debug, Clone, PartialEq)]
pub struct RatedPointDecl {
    /// Rated volumetric flow (m^3/s).
    pub flow: QtyAny,
    /// Rated static pressure (Pa).
    pub static_pressure: QtyAny,
    /// How this point participates in the solve.
    pub admission: RatedPointAdmission,
}

/// One fan bank declaration: every solver-consumed parameter, explicit.
#[derive(Debug, Clone, PartialEq)]
pub struct FanBankDecl {
    /// Stable bank identity.
    pub bank_id: String,
    /// Declared pressure-flow curve with provenance and stall boundary.
    pub curve: FanCurveDecl,
    /// Identical fans in the bank; explicit, never defaulted.
    pub count: usize,
    /// Arrangement of the identical fans inside the bank.
    pub arrangement: FanArrangement,
    /// Declared fan-law speed ratio, inside the declared domain.
    pub speed_ratio: f64,
    /// Admitted speed-ratio domain: `0 < low <= 1 <= high`.
    pub speed_ratio_domain: (f64, f64),
    /// Optional rated point with its admission classification.
    pub rated_point: Option<RatedPointDecl>,
}

/// The explicit system topology over declared banks.
#[derive(Debug, Clone, PartialEq)]
pub enum FanSystemTopology {
    /// Exactly one bank; no composition.
    Single,
    /// Named banks in series: equal total flow, pressures add.
    Series(Vec<String>),
    /// Named banks in parallel: shared pressure, flows add.
    Parallel(Vec<String>),
}

/// The versioned fan-system declaration.
#[derive(Debug, Clone, PartialEq)]
pub struct FanSystemDecl {
    /// Declaration schema version ([`FAN_SYSTEM_DECL_VERSION`]).
    pub version: u32,
    /// Declared banks.
    pub banks: Vec<FanBankDecl>,
    /// Explicit system topology. Required when more than one bank is
    /// declared; `Single` is the only legal form for one bank.
    pub topology: FanSystemTopology,
}

fn validate_bank_identity(bank_id: &str) -> Result<(), ProjectError> {
    let valid = !bank_id.is_empty()
        && bank_id.len() <= 128
        && bank_id
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b':'));
    if !valid {
        return Err(fan_error(
            "fan-bank-identity",
            format!("bank identity {bank_id:?} is not a bounded token"),
            "use 1..=128 ASCII alphanumeric, '-', '_', '.', or ':' bytes",
        ));
    }
    Ok(())
}

fn validate_curve(curve: &FanCurveDecl, bank_id: &str) -> Result<(), ProjectError> {
    if curve.points.len() < 2 {
        return Err(fan_error(
            "fan-curve-points",
            format!(
                "bank {bank_id:?} declares {} curve points; at least two are required",
                curve.points.len()
            ),
            "declare at least a two-point curve",
        ));
    }
    for point in &curve.points {
        if point.flow.dims != FLOW_DIMS {
            return Err(fan_error(
                "fan-curve-units",
                format!(
                    "bank {bank_id:?} curve point flow carries dims {}",
                    point.flow.dims.unit_string()
                ),
                "curve flow must carry m^3/s dimensions",
            ));
        }
        if point.static_pressure.dims != PRESSURE_DIMS {
            return Err(fan_error(
                "fan-curve-units",
                format!(
                    "bank {bank_id:?} curve point pressure carries dims {}",
                    point.static_pressure.dims.unit_string()
                ),
                "curve pressure must carry Pa dimensions",
            ));
        }
        let q = point.flow.value;
        let p = point.static_pressure.value;
        if !(q.is_finite() && q >= 0.0 && p.is_finite() && p >= 0.0) {
            return Err(fan_error(
                "fan-curve-value",
                format!("bank {bank_id:?} declares a non-finite or negative curve point"),
                "curve points must be finite and non-negative",
            ));
        }
    }
    for pair in curve.points.windows(2) {
        if pair[1].flow.value <= pair[0].flow.value {
            return Err(fan_error(
                "fan-curve-monotonicity",
                format!("bank {bank_id:?} curve flow does not strictly increase"),
                "declare knots in strictly increasing flow order",
            ));
        }
        if pair[1].static_pressure.value > pair[0].static_pressure.value {
            return Err(fan_error(
                "fan-curve-pressure-rise",
                format!("bank {bank_id:?} curve pressure rises with flow"),
                "a fan curve's pressure must not increase with flow",
            ));
        }
    }
    if !(curve.pressure_tolerance_rel.is_finite()
        && (0.0..1.0).contains(&curve.pressure_tolerance_rel))
    {
        return Err(fan_error(
            "fan-curve-tolerance",
            format!(
                "bank {bank_id:?} tolerance {} is outside [0, 1)",
                curve.pressure_tolerance_rel
            ),
            "declare a finite relative pressure tolerance below 1",
        ));
    }
    if curve.source.trim().is_empty() || curve.source_id.trim().is_empty() {
        return Err(fan_error(
            "fan-curve-provenance",
            format!("bank {bank_id:?} omits the source citation or identifier"),
            "both the citation and the stable source identifier are required",
        ));
    }
    let q_min = curve.points[0].flow.value;
    let q_max = curve.points.last().expect("len checked").flow.value;
    let admissible = curve.min_flow.value;
    if curve.min_flow.dims != FLOW_DIMS {
        return Err(fan_error(
            "fan-curve-units",
            format!(
                "bank {bank_id:?} stall boundary carries dims {}",
                curve.min_flow.dims.unit_string()
            ),
            "the stall boundary must carry m^3/s dimensions",
        ));
    }
    if !(admissible.is_finite() && admissible >= q_min && admissible < q_max) {
        return Err(fan_error(
            "fan-stall-boundary",
            format!("bank {bank_id:?} stall boundary {admissible} is outside [{q_min}, {q_max})"),
            "the stall boundary must cover the first knot and stay below the last",
        ));
    }
    Ok(())
}

impl FanSystemDecl {
    /// Validate a fan-system declaration: every solver-consumed parameter
    /// present, explicit, and consistent. Nothing is inferred or defaulted.
    ///
    /// # Errors
    /// Returns a structured [`ProjectError`] for empty/duplicate/orphan
    /// banks, ambiguous or incomplete topology, out-of-domain speed ratios,
    /// malformed curves, or a rated point that fails its declared check.
    #[allow(clippy::too_many_lines)] // one admission protocol; splitting scatters the refusal taxonomy
    pub fn validate(&self) -> Result<(), ProjectError> {
        if self.version != FAN_SYSTEM_DECL_VERSION {
            return Err(fan_error(
                "fan-system-version",
                format!(
                    "declaration version {} is not {FAN_SYSTEM_DECL_VERSION}",
                    self.version
                ),
                "migrate the declaration explicitly; there is no compatibility default",
            ));
        }
        if self.banks.is_empty() {
            return Err(fan_error(
                "fan-system-empty",
                "the fan system declares no banks".to_string(),
                "declare at least one fan bank",
            ));
        }
        for (index, bank) in self.banks.iter().enumerate() {
            validate_bank_identity(&bank.bank_id)?;
            if self.banks[..index]
                .iter()
                .any(|prior| prior.bank_id == bank.bank_id)
            {
                return Err(fan_error(
                    "fan-bank-duplicate",
                    format!("bank identity {:?} is declared twice", bank.bank_id),
                    "bank identities must be unique",
                ));
            }
            if bank.count == 0 {
                return Err(fan_error(
                    "fan-bank-count",
                    format!("bank {:?} declares zero fans", bank.bank_id),
                    "the fan count is explicit and must be at least one",
                ));
            }
            let (low, high) = bank.speed_ratio_domain;
            if !(low.is_finite() && high.is_finite() && low > 0.0 && low <= 1.0 && high >= 1.0) {
                return Err(fan_error(
                    "fan-speed-domain",
                    format!(
                        "bank {:?} declares speed-ratio domain [{low}, {high}]",
                        bank.bank_id
                    ),
                    "the admitted domain must satisfy 0 < low <= 1 <= high",
                ));
            }
            if !(bank.speed_ratio.is_finite()
                && bank.speed_ratio >= low
                && bank.speed_ratio <= high)
            {
                return Err(fan_error(
                    "fan-speed-out-of-domain",
                    format!(
                        "bank {:?} declares speed ratio {} outside [{low}, {high}]",
                        bank.bank_id, bank.speed_ratio
                    ),
                    "the declared speed ratio must sit inside the declared domain",
                ));
            }
            validate_curve(&bank.curve, &bank.bank_id)?;
            if let Some(rated) = &bank.rated_point {
                Self::validate_rated_point(bank, rated)?;
            }
        }
        match &self.topology {
            FanSystemTopology::Single => {
                if self.banks.len() != 1 {
                    return Err(fan_error(
                        "fan-system-topology",
                        format!("Single topology with {} banks declared", self.banks.len()),
                        "a single-bank system declares exactly one bank; more than one bank requires an explicit Series or Parallel member list",
                    ));
                }
            }
            FanSystemTopology::Series(members) | FanSystemTopology::Parallel(members) => {
                if members.len() < 2 {
                    return Err(fan_error(
                        "fan-system-topology",
                        "a composite topology names fewer than two members".to_string(),
                        "name at least two member bank identities",
                    ));
                }
                for member in members {
                    if !self.banks.iter().any(|bank| &bank.bank_id == member) {
                        return Err(fan_error(
                            "fan-system-orphan",
                            format!("topology member {member:?} is not a declared bank"),
                            "every topology member must be a declared bank identity",
                        ));
                    }
                }
                for bank in &self.banks {
                    if !members.contains(&bank.bank_id) {
                        return Err(fan_error(
                            "fan-system-unreferenced",
                            format!(
                                "declared bank {:?} never appears in the topology",
                                bank.bank_id
                            ),
                            "every declared bank must be consumed by the topology exactly once",
                        ));
                    }
                }
                let mut sorted = members.clone();
                sorted.sort();
                for window in sorted.windows(2) {
                    if window[0] == window[1] {
                        return Err(fan_error(
                            "fan-system-duplicate-member",
                            format!("topology member {:?} appears twice", window[0]),
                            "each member may appear exactly once",
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    fn validate_rated_point(
        bank: &FanBankDecl,
        rated: &RatedPointDecl,
    ) -> Result<(), ProjectError> {
        if rated.flow.dims != FLOW_DIMS || rated.static_pressure.dims != PRESSURE_DIMS {
            return Err(fan_error(
                "fan-rated-units",
                format!(
                    "bank {:?} rated point carries non-flow or non-pressure units",
                    bank.bank_id
                ),
                "rated flow must be m^3/s and rated pressure Pa",
            ));
        }
        match rated.admission {
            RatedPointAdmission::CorrelationOnly => Ok(()),
            RatedPointAdmission::CheckedWithinDeclaredTolerance => {
                let curve = &bank.curve;
                let q = rated.flow.value;
                let q_min = curve.points[0].flow.value;
                let q_max = curve.points.last().expect("len checked").flow.value;
                if !(q.is_finite() && q >= q_min && q <= q_max) {
                    return Err(fan_error(
                        "fan-rated-off-curve",
                        format!(
                            "bank {:?} rated flow {q} is off the declared curve",
                            bank.bank_id
                        ),
                        "the rated point must sit on the declared curve's flow range",
                    ));
                }
                let expected = interpolate_curve(curve, q);
                let tolerance = curve.pressure_tolerance_rel * expected;
                let deviation = (rated.static_pressure.value - expected).abs();
                if deviation > tolerance {
                    return Err(fan_error(
                        "fan-rated-mismatch",
                        format!(
                            "bank {:?} rated pressure {} deviates from the declared curve value {expected} by {deviation}, beyond the declared relative tolerance {}",
                            bank.bank_id, rated.static_pressure.value, curve.pressure_tolerance_rel
                        ),
                        "check the curve against the rating or reclassify the point as correlation-only",
                    ));
                }
                Ok(())
            }
        }
    }

    /// Domain-separated declaration identity binding the version, every
    /// bank field, and the topology.
    #[must_use]
    pub fn identity(&self) -> String {
        let mut hasher = fs_blake3::DomainHasher::new(IDENTITY_DOMAIN);
        hasher.update(&self.version.to_le_bytes());
        hasher.update(&(self.banks.len() as u64).to_le_bytes());
        for bank in &self.banks {
            push_identity_string(&mut hasher, &bank.bank_id);
            hasher.update(&(bank.count as u64).to_le_bytes());
            hasher.update(&[match bank.arrangement {
                FanArrangement::Series => 0,
                FanArrangement::Parallel => 1,
            }]);
            hasher.update(&bank.speed_ratio.to_le_bytes());
            hasher.update(&bank.speed_ratio_domain.0.to_le_bytes());
            hasher.update(&bank.speed_ratio_domain.1.to_le_bytes());
            push_identity_string(&mut hasher, &bank.curve.source);
            push_identity_string(&mut hasher, &bank.curve.source_id);
            hasher.update(&(bank.curve.points.len() as u64).to_le_bytes());
            for point in &bank.curve.points {
                hasher.update(&point.flow.value.to_le_bytes());
                hasher.update(&point.static_pressure.value.to_le_bytes());
            }
            hasher.update(&bank.curve.pressure_tolerance_rel.to_le_bytes());
            hasher.update(&[match bank.curve.tolerance_basis {
                FanToleranceBasis::Manufacturer => 0,
                FanToleranceBasis::Engineering => 1,
                FanToleranceBasis::Analytic => 2,
            }]);
            hasher.update(&bank.curve.min_flow.value.to_le_bytes());
            match &bank.rated_point {
                Some(rated) => {
                    hasher.update(&[1]);
                    hasher.update(&rated.flow.value.to_le_bytes());
                    hasher.update(&rated.static_pressure.value.to_le_bytes());
                    hasher.update(&[match rated.admission {
                        RatedPointAdmission::CheckedWithinDeclaredTolerance => 0,
                        RatedPointAdmission::CorrelationOnly => 1,
                    }]);
                }
                None => hasher.update(&[0]),
            }
        }
        match &self.topology {
            FanSystemTopology::Single => hasher.update(&[0]),
            FanSystemTopology::Series(members) => {
                hasher.update(&[1]);
                for member in members {
                    push_identity_string(&mut hasher, member);
                }
            }
            FanSystemTopology::Parallel(members) => {
                hasher.update(&[2]);
                for member in members {
                    push_identity_string(&mut hasher, member);
                }
            }
        }
        format!("{FAN_SYSTEM_IDENTITY_PREFIX}{}", hasher.finalize())
    }
}

fn push_identity_string(hasher: &mut fs_blake3::DomainHasher, value: &str) {
    hasher.update(&(value.len() as u64).to_le_bytes());
    hasher.update(value.as_bytes());
}

/// Exact piecewise-linear curve pressure at flow `q`; `q` is inside the
/// declared range by the caller's check.
fn interpolate_curve(curve: &FanCurveDecl, q: f64) -> f64 {
    for pair in curve.points.windows(2) {
        let q0 = pair[0].flow.value;
        let q1 = pair[1].flow.value;
        if q >= q0 && q <= q1 {
            let p0 = pair[0].static_pressure.value;
            let p1 = pair[1].static_pressure.value;
            return p0 + (p1 - p0) * (q - q0) / (q1 - q0);
        }
    }
    curve.points.last().expect("nonempty").static_pressure.value
}

/// The production lowering product: the solver-ready composite fan bank
/// plus the member banks and the topology record that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct LoweredFanSystem {
    /// The production fan bank: the single declared bank, or the exact
    /// composite for multi-bank systems.
    pub system_bank: FanBank,
    /// Member banks in declaration order, one per declared bank identity.
    pub members: Vec<(String, FanBank)>,
    /// The topology the lowering consumed.
    pub topology: FanSystemTopology,
    /// The declaration identity this lowering binds.
    pub declaration_identity: String,
}

/// Lower a validated declaration into production `fs-airflow` evidence:
/// every declared bank becomes exactly one [`FanBank`], and a multi-bank
/// topology composes them exactly. Rated points typed correlation-only are
/// excluded from the solve input; checked points were verified against
/// their declared curves at declaration time.
///
/// # Errors
/// Returns a structured [`ProjectError`] from declaration validation or
/// the production curve/bank/composite constructors.
pub fn lower_fan_system(decl: &FanSystemDecl) -> Result<LoweredFanSystem, ProjectError> {
    decl.validate()?;
    let mut members = Vec::with_capacity(decl.banks.len());
    for bank in &decl.banks {
        let curve = FanCurve::new(
            bank.bank_id.clone(),
            bank.curve
                .points
                .iter()
                .map(|point| {
                    FanPoint::new(
                        fs_qty::VolumetricFlowRate::new(point.flow.value),
                        fs_qty::Pressure::new(point.static_pressure.value),
                    )
                })
                .collect(),
            SourceProvenance::new(bank.curve.source.clone(), bank.curve.source_id.clone()),
            bank.curve.pressure_tolerance_rel,
            match bank.curve.tolerance_basis {
                FanToleranceBasis::Manufacturer => ToleranceBasis::ManufacturerDeclared,
                FanToleranceBasis::Engineering => ToleranceBasis::EngineeringAllowance,
                FanToleranceBasis::Analytic => ToleranceBasis::Analytic,
            },
            fs_qty::VolumetricFlowRate::new(bank.curve.min_flow.value),
            bank.speed_ratio_domain,
        )
        .map_err(|error| {
            fan_error(
                "fan-curve-lowering",
                format!(
                    "bank {:?} curve refused by production admission: {error}",
                    bank.bank_id
                ),
                "the declaration validator and the production constructor disagree; report this",
            )
        })?;
        let member = FanBank::new(curve, bank.count, bank.arrangement, bank.speed_ratio).map_err(
            |error| {
                fan_error(
                    "fan-bank-lowering",
                    format!("bank {:?} refused by production admission: {error}", bank.bank_id),
                    "the declaration validator and the production constructor disagree; report this",
                )
            },
        )?;
        members.push((bank.bank_id.clone(), member));
    }
    let system_bank = match &decl.topology {
        FanSystemTopology::Single => members[0].1.clone(),
        FanSystemTopology::Series(member_ids) => {
            let ordered: Vec<FanBank> = member_ids
                .iter()
                .map(|id| {
                    members
                        .iter()
                        .find(|(bank_id, _)| bank_id == id)
                        .expect("validated member")
                        .1
                        .clone()
                })
                .collect();
            compose_series(&ordered).map_err(|error| {
                fan_error(
                    "fan-system-composition",
                    format!("series composition refused: {error}"),
                    "the member banks share no common flow domain",
                )
            })?
        }
        FanSystemTopology::Parallel(member_ids) => {
            let ordered: Vec<FanBank> = member_ids
                .iter()
                .map(|id| {
                    members
                        .iter()
                        .find(|(bank_id, _)| bank_id == id)
                        .expect("validated member")
                        .1
                        .clone()
                })
                .collect();
            compose_parallel(&ordered).map_err(|error| {
                fan_error(
                    "fan-system-composition",
                    format!("parallel composition refused: {error}"),
                    "the member banks share no common pressure domain",
                )
            })?
        }
    };
    Ok(LoweredFanSystem {
        system_bank,
        members,
        topology: decl.topology.clone(),
        declaration_identity: decl.identity(),
    })
}
