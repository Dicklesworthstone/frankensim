//! Measured, multi-axis convergence studies (IRP M5.5).
//!
//! Coordinates are module-reported resolutions, never inferred from requests.
//! Power-law fits use all admitted rungs and the shared convergence arithmetic.
//! A fit and a GCI are numerical estimates, not continuum certificates.

use fs_blake3::{ContentHash, hash_domain};

use crate::convergence::{PowerFit, fit_power_law};

mod equivalence;
pub use equivalence::{
    BoundSource, BoundedSolution, Coverage, EquivalenceDomain, EquivalenceReport,
    EquivalenceVerdict, Metric, compare_equivalence, qoi_bound,
};

/// Applicable numerical error model; no model is selected implicitly.
#[derive(Debug, Clone, PartialEq)]
pub enum ErrorModel {
    /// A leading power of the resolution, with a declared scheme order.
    Algebraic { formal_order: f64 },
    /// Declared truncation/tail estimate, in each QoI's units.
    SpectralTail,
    /// Externally supplied operator-error bound.
    OperatorError,
    /// Externally supplied residual-to-QoI bound.
    ResidualBound,
    /// Externally supplied statistical uncertainty.
    Statistical,
    /// Missing error model. Never qualifies a convergence claim.
    Unknown,
}

/// Map a realized resolution into a decreasing error coordinate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Coordinate {
    /// Mesh spacing, time step, or tolerance.
    Direct,
    /// Mode/sample/panel count: the error coordinate is its reciprocal.
    Reciprocal,
}

/// One explicitly declared refinement axis.
#[derive(Debug, Clone, PartialEq)]
pub struct RefinableAxis {
    /// Unique name also used in measured dimension records.
    pub name: String,
    /// Unit of the realized coordinate.
    pub unit: String,
    /// Declared numerical error model.
    pub model: ErrorModel,
    /// Coupled group name, if this axis must be refined with others.
    pub group: Option<String>,
    /// Meaning of increasing/decreasing realized values.
    pub coordinate: Coordinate,
    /// Exact observables bounded by a non-algebraic model. Other observables
    /// remain measurements only; a residual cannot silently bound pitch/dB.
    pub bound_targets: Vec<String>,
}

/// Named physical observable. Units and order must agree across rungs.
#[derive(Debug, Clone, PartialEq)]
pub struct Observable {
    /// Stable observable name.
    pub name: String,
    /// Measured value, before any audio normalization.
    pub value: f64,
    /// Physical/reporting unit.
    pub unit: String,
}

/// A dimension or tolerance actually used by the module.
#[derive(Debug, Clone, PartialEq)]
pub struct Dimension {
    /// Axis name.
    pub name: String,
    /// Realized positive finite value.
    pub value: f64,
    /// Axis unit.
    pub unit: String,
}

/// One completed measurement. No authored residual or time defaults exist.
#[derive(Debug, Clone, PartialEq)]
pub struct RungRecord {
    /// Requested axis value, retained to expose clamping.
    pub requested: f64,
    /// Realized dimensions/tolerances supplied by the producer.
    pub dimensions: Vec<Dimension>,
    /// Raw observables in declared order.
    pub qoi: Vec<Observable>,
    /// Per-observable declared bounds, in the matching observable's units.
    pub bounds: Vec<Observable>,
    /// Register/slot/regime identity; changes refuse fitting.
    pub regime: String,
    /// Measured solver residual (its norm and unit belong in producer provenance).
    pub solver_residual: f64,
    /// Measured elapsed seconds.
    pub wall_s: f64,
}

/// A separate physical fixture demonstrating that a resolution can affect a QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct SensitivityFixture {
    /// Reproducible fixture identifier, including its source description.
    pub fixture: String,
    /// Two actual, distinct module resolutions in axis units.
    pub realized: [f64; 2],
    /// Same named observable evaluated at those resolutions.
    pub qoi: [Observable; 2],
    /// Declared absolute difference needed to resolve sensitivity.
    pub minimum_change: f64,
}

/// Explicit analysis thresholds. Callers must declare these before a run.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisPolicy {
    /// Largest admissible power in the bounded nonlinear fit.
    pub maximum_order: f64,
    /// Allowed absolute deviation from the declared formal order.
    pub order_tolerance: f64,
    /// Maximum least-squares RMS residual divided by the data range.
    pub relative_fit_tolerance: f64,
    /// Maximum relative spread of adjacent power-law coefficients.
    pub asymptotic_tolerance: f64,
    /// Minimum adjacent ratio in decreasing error coordinates.
    pub minimum_refinement_ratio: f64,
}

impl AnalysisPolicy {
    pub(crate) fn valid(&self) -> bool {
        self.maximum_order.is_finite()
            && self.maximum_order > 0.05
            && self.maximum_order <= 64.0
            && self.minimum_refinement_ratio.is_finite()
            && self.minimum_refinement_ratio > 1.0
            && [
                self.order_tolerance,
                self.relative_fit_tolerance,
                self.asymptotic_tolerance,
            ]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
    }
}

/// Per-observable interpretation, with no fallback GCI.
#[derive(Debug, Clone, PartialEq)]
pub enum AxisVerdict {
    /// All-rung fit, order and coefficient consistency passed.
    Asymptotic {
        observed_order: f64,
        gci: Option<f64>,
    },
    /// A monotone power-law fit exists, but asymptotic assumptions did not pass.
    PreAsymptotic { observed_order: f64 },
    /// All-rung least-squares residual is too large for the declared fit model.
    Scattered,
    /// Successive observable differences change sign.
    Oscillatory,
    /// A declared bound for the finest admitted rung, not an inferred GCI.
    Bounded { bound: f64 },
    /// Exact/roundoff-flat values: no observable order and no GCI.
    Unresolved,
    /// Measured QoI outside the declared bound's target; no accuracy claim.
    ObservedOnly {
        /// Missing map from the model's diagnostic norm to this observable.
        reason: String,
    },
    /// Fewer than three distinct, valid realized resolutions.
    InsufficientRungs,
    /// Missing model, malformed measurements, sensitivity or policy.
    Refused { reason: String },
}

/// All-rung result for a named observable.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiReport {
    /// Name and unit of the observable.
    pub name: String,
    /// Observable unit.
    pub unit: String,
    /// Measured max-minus-min spread (data, not an acceptance gate).
    pub spread: f64,
    /// Applicable error-model result.
    pub verdict: AxisVerdict,
    /// Full fit diagnostics for algebraic data.
    pub fit: Option<PowerFit>,
}

/// Axis report retaining every input, including dead and failed rungs.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisReport {
    /// Declared axis.
    pub axis: RefinableAxis,
    /// Unmodified input measurements.
    pub rungs: Vec<RungRecord>,
    /// Indices of repeated realized resolutions (not counted as live).
    pub dead_rungs: Vec<usize>,
    /// Indices of invalid or changed-regime measurements.
    pub failed_rungs: Vec<usize>,
    /// Number of distinct, valid resolutions.
    pub live_rungs: usize,
    /// Independent liveness evidence, retained even on refusal.
    pub sensitivity: SensitivityFixture,
    /// One result for every declared observable.
    pub qois: Vec<QoiReport>,
}

impl AxisReport {
    /// Complete measurements are distinct from a Q2 accuracy acceptance.
    #[must_use]
    pub fn complete(&self) -> bool {
        self.live_rungs >= 3
            && self.failed_rungs.is_empty()
            && !self.qois.is_empty()
            && self.qois.iter().all(|q| {
                !matches!(
                    q.verdict,
                    AxisVerdict::InsufficientRungs | AxisVerdict::Refused { .. }
                )
            })
    }
}

fn dimension<'a>(rung: &'a RungRecord, axis: &RefinableAxis) -> Option<&'a Dimension> {
    let mut matches = rung.dimensions.iter().filter(|d| d.name == axis.name);
    let d = matches.next()?;
    (matches.next().is_none() && d.unit == axis.unit && d.value.is_finite() && d.value > 0.0)
        .then_some(d)
}

fn compatible(a: &[Observable], b: &[Observable]) -> bool {
    !a.is_empty()
        && a.len() == b.len()
        && a.iter()
            .enumerate()
            .all(|(i, x)| !a[..i].iter().any(|q| q.name == x.name))
        && a.iter().zip(b).all(|(x, y)| {
            !x.name.is_empty()
                && !x.unit.is_empty()
                && x.name == y.name
                && x.unit == y.unit
                && x.value.is_finite()
                && y.value.is_finite()
        })
}

fn sensitivity_valid(s: &SensitivityFixture) -> bool {
    !s.fixture.is_empty()
        && s.realized.iter().all(|v| v.is_finite() && *v > 0.0)
        && s.realized[0] != s.realized[1]
        && compatible(&s.qoi[..1], &s.qoi[1..])
        && s.minimum_change.is_finite()
        && s.minimum_change >= 0.0
        && (s.qoi[0].value - s.qoi[1].value).abs() > s.minimum_change
}

/// Analyze every valid rung, counting a repeated realized dimension only once.
/// A regime change or malformed rung refuses the study rather than being dropped.
#[must_use]
pub fn evaluate_axis(
    axis: &RefinableAxis,
    rungs: &[RungRecord],
    sensitivity: &SensitivityFixture,
    policy: &AnalysisPolicy,
) -> AxisReport {
    let mut report = AxisReport {
        axis: axis.clone(),
        rungs: rungs.to_vec(),
        dead_rungs: vec![],
        failed_rungs: vec![],
        live_rungs: 0,
        sensitivity: sensitivity.clone(),
        qois: vec![],
    };
    let Some(first) = rungs.first() else {
        return report;
    };
    let mut admitted = Vec::new();
    let mut seen = Vec::new();
    for (i, r) in rungs.iter().enumerate() {
        let Some(d) = dimension(r, axis) else {
            report.failed_rungs.push(i);
            continue;
        };
        if !r.requested.is_finite()
            || r.requested <= 0.0
            || r.regime.is_empty()
            || r.regime != first.regime
            || !compatible(&first.qoi, &r.qoi)
            || !r.solver_residual.is_finite()
            || r.solver_residual < 0.0
            || !r.wall_s.is_finite()
            || r.wall_s < 0.0
        {
            report.failed_rungs.push(i);
            continue;
        }
        if seen.contains(&d.value) {
            report.dead_rungs.push(i);
            continue;
        }
        seen.push(d.value);
        let h = match axis.coordinate {
            Coordinate::Direct => d.value,
            Coordinate::Reciprocal => 1.0 / d.value,
        };
        if !h.is_finite() || h <= 0.0 {
            report.failed_rungs.push(i);
            continue;
        }
        admitted.push((h, r));
    }
    admitted.sort_by(|a, b| b.0.total_cmp(&a.0));
    report.live_rungs = admitted.len();
    let refusal = if axis.name.is_empty() || axis.unit.is_empty() || !policy.valid() {
        Some("axis declaration or analysis policy is invalid")
    } else if !report.failed_rungs.is_empty() {
        Some("invalid measurement or regime change")
    } else if !sensitivity_valid(sensitivity) {
        Some("independent sensitivity fixture did not resolve a QoI change")
    } else if matches!(axis.model, ErrorModel::Unknown) {
        Some("unknown error model")
    } else if !matches!(axis.model, ErrorModel::Algebraic { .. })
        && (axis.bound_targets.is_empty()
            || axis.bound_targets.iter().enumerate().any(|(i, n)| {
                axis.bound_targets[..i].contains(n) || !first.qoi.iter().any(|q| q.name == *n)
            }))
    {
        Some("non-algebraic bound targets are missing, duplicated or unknown")
    } else {
        None
    };
    for (j, q) in first.qoi.iter().enumerate() {
        let values: Vec<_> = admitted.iter().map(|(h, r)| (*h, r.qoi[j].value)).collect();
        let spread = values
            .iter()
            .map(|(_, y)| *y)
            .fold(f64::NEG_INFINITY, f64::max)
            - values.iter().map(|(_, y)| *y).fold(f64::INFINITY, f64::min);
        let mut fit = None;
        let verdict = if let Some(reason) = refusal {
            AxisVerdict::Refused {
                reason: reason.into(),
            }
        } else if admitted.len() < 3 {
            AxisVerdict::InsufficientRungs
        } else if !spread.is_finite() {
            AxisVerdict::Refused {
                reason: "observable spread is not representable".into(),
            }
        } else if admitted
            .windows(2)
            .any(|w| w[0].0 / w[1].0 < policy.minimum_refinement_ratio)
        {
            AxisVerdict::Refused {
                reason: "realized refinement ratios are too small".into(),
            }
        } else if let ErrorModel::Algebraic { formal_order } = axis.model {
            let analysis = fit_power_law(&values, formal_order, policy);
            let v = analysis.verdict.clone();
            fit = Some(analysis);
            v
        } else if !axis.bound_targets.contains(&q.name) {
            AxisVerdict::ObservedOnly {
                reason: "accuracy bound unavailable: no declared conditioning/stability map from this axis diagnostic into this QoI".into(),
            }
        } else {
            let bounds: Option<Vec<f64>> = admitted
                .iter()
                .map(|(_, r)| {
                    let mut bs = r
                        .bounds
                        .iter()
                        .filter(|b| b.name == q.name && b.unit == q.unit);
                    let b = bs.next()?;
                    (bs.next().is_none() && b.value.is_finite() && b.value >= 0.0)
                        .then_some(b.value)
                })
                .collect();
            match bounds {
                Some(bs) => AxisVerdict::Bounded {
                    bound: bs[bs.len() - 1],
                },
                None => AxisVerdict::Refused {
                    reason: "every non-algebraic rung needs a bound in this QoI's units".into(),
                },
            }
        };
        report.qois.push(QoiReport {
            name: q.name.clone(),
            unit: q.unit.clone(),
            spread: if spread.is_finite() { spread } else { 0.0 },
            verdict,
            fit,
        });
    }
    report
}

/// Coupled group with a predeclared product invariant, e.g. Δt * N_FFT = seconds.
#[derive(Debug, Clone, PartialEq)]
pub struct CoupledGroup {
    /// Group identity referenced by every member axis.
    pub name: String,
    /// Axis names and powers in the product invariant.
    pub members: Vec<(String, f64)>,
    /// Fixed physical product, declared before measurement.
    pub invariant: f64,
    /// Relative tolerance for the product invariant.
    pub relative_tolerance: f64,
}

/// Refuse missing/duplicate group members, dead members or a moving invariant.
pub fn validate_group(
    group: &CoupledGroup,
    axes: &[RefinableAxis],
    rungs: &[RungRecord],
) -> Result<(), String> {
    if group.name.is_empty()
        || group.members.len() < 2
        || rungs.len() < 3
        || !group.invariant.is_finite()
        || group.invariant <= 0.0
        || !group.relative_tolerance.is_finite()
        || group.relative_tolerance < 0.0
    {
        return Err("invalid coupled group declaration or fewer than three joint rungs".into());
    }
    let mut names = Vec::new();
    let mut previous: Vec<Option<f64>> = vec![None; group.members.len()];
    for (name, exponent) in &group.members {
        if !exponent.is_finite() || *exponent == 0.0 || names.contains(name) {
            return Err("duplicate coupled member or invalid exponent".into());
        }
        names.push(name.clone());
        if axes
            .iter()
            .filter(|a| a.name == *name && a.group.as_deref() == Some(&group.name))
            .count()
            != 1
        {
            return Err(format!(
                "coupled axis {name} is missing or not in {}",
                group.name
            ));
        }
    }
    for r in rungs {
        let mut product = 1.0;
        for (j, (name, exponent)) in group.members.iter().enumerate() {
            let axis = axes
                .iter()
                .find(|a| a.name == *name)
                .ok_or("missing group axis")?;
            let d = dimension(r, axis).ok_or("missing realized coupled dimension")?;
            if previous[j] == Some(d.value) {
                return Err(format!("coupled axis {name} is dead"));
            }
            previous[j] = Some(d.value);
            product *= fs_math::det::pow(d.value, *exponent);
        }
        if !product.is_finite()
            || (product / group.invariant - 1.0).abs() > group.relative_tolerance
        {
            return Err(format!(
                "coupled group {} violated its physical invariant",
                group.name
            ));
        }
    }
    Ok(())
}

/// Seed ensemble is mandatory even for a deterministic source.
#[derive(Debug, Clone, PartialEq)]
pub enum EnsembleDecl {
    /// No stochastic source, with an explicit reason.
    Deterministic { reason: String },
    /// Actual seeds and the independently declared statistical interpretation.
    Seeds {
        seeds: Vec<u64>,
        uncertainty: String,
    },
}

/// Signed finest-rung-reference differences, in observable units.
#[derive(Debug, Clone, PartialEq)]
pub struct Interaction {
    /// Coarse-to-fine joint-rung index.
    pub rung: usize,
    /// Observable name/unit and joint difference from its finest rung.
    pub joint_error: Observable,
    /// Sum of the corresponding single-axis differences from their finest rungs.
    pub single_axis_sum: f64,
    /// Joint difference minus the single-axis sum. Not a rigorous error bound.
    pub interaction: f64,
}

/// Joint run with retained realized dimensions and interaction measurements.
#[derive(Debug, Clone, PartialEq)]
pub struct JointReport {
    /// All joint rung records.
    pub rungs: Vec<RungRecord>,
    /// Every rung/observable interaction (finest rung is exactly zero).
    pub interactions: Vec<Interaction>,
}

/// Complete reproducible input and analysis, in declaration order.
#[derive(Debug, Clone, PartialEq)]
pub struct LadderReport {
    /// Predeclared thresholds.
    pub policy: AnalysisPolicy,
    /// Every axis and its raw measurements.
    pub axes: Vec<AxisReport>,
    /// Coupled groups checked against the joint measurements.
    pub groups: Vec<CoupledGroup>,
    /// Jointly refined end-to-end run.
    pub joint: JointReport,
    /// Explicit ensemble/absence declaration.
    pub ensemble: EnsembleDecl,
    /// Full-width content identity of the retained report bytes.
    pub hash: ContentHash,
}

impl LadderReport {
    /// Validate complete axis coverage, matched observables and simultaneous
    /// refinement, then measure joint interactions using every rung.
    pub fn new(
        policy: AnalysisPolicy,
        axes: Vec<AxisReport>,
        groups: Vec<CoupledGroup>,
        joint: Vec<RungRecord>,
        ensemble: EnsembleDecl,
    ) -> Result<Self, String> {
        if !policy.valid()
            || axes.is_empty()
            || axes.iter().any(|a| !a.complete())
            || joint.len() < 3
        {
            return Err("complete axes and at least three joint rungs are required".into());
        }
        let ensemble_valid = match &ensemble {
            EnsembleDecl::Deterministic { reason } => !reason.is_empty(),
            EnsembleDecl::Seeds { seeds, uncertainty } => {
                !seeds.is_empty()
                    && !uncertainty.is_empty()
                    && seeds
                        .iter()
                        .enumerate()
                        .all(|(i, s)| !seeds[..i].contains(s))
            }
        };
        if !ensemble_valid {
            return Err("ensemble declaration is empty or has duplicate seeds".into());
        }
        let declarations: Vec<_> = axes.iter().map(|a| a.axis.clone()).collect();
        for (i, a) in declarations.iter().enumerate() {
            if declarations[..i].iter().any(|b| b.name == a.name) {
                return Err("duplicate axis".into());
            }
            if let Some(g) = &a.group {
                if groups
                    .iter()
                    .filter(|x| x.name == *g && x.members.iter().any(|(n, _)| *n == a.name))
                    .count()
                    != 1
                {
                    return Err(format!("axis {} lacks its coupled group", a.name));
                }
            }
        }
        for g in &groups {
            validate_group(g, &declarations, &joint)?;
        }
        let finest = joint.last().ok_or("joint measurements absent")?;
        let mut interactions = Vec::new();
        for (i, r) in joint.iter().enumerate() {
            if !compatible(&finest.qoi, &r.qoi)
                || r.regime != finest.regime
                || r.regime.is_empty()
                || !r.solver_residual.is_finite()
                || r.solver_residual < 0.0
                || !r.wall_s.is_finite()
                || r.wall_s < 0.0
            {
                return Err("joint rung is invalid or changes regime".into());
            }
            for a in &axes {
                let d = dimension(r, &a.axis).ok_or("joint rung omits an axis")?;
                if i > 0 {
                    let prev = dimension(&joint[i - 1], &a.axis)
                        .ok_or("previous joint dimension absent")?;
                    let ratio = match a.axis.coordinate {
                        Coordinate::Direct => prev.value / d.value,
                        Coordinate::Reciprocal => d.value / prev.value,
                    };
                    if ratio < policy.minimum_refinement_ratio {
                        return Err(format!("joint axis {} did not refine", a.axis.name));
                    }
                }
                if a.rungs.len() != joint.len() || !compatible(&r.qoi, &a.rungs[i].qoi) {
                    return Err("joint and single-axis rung/QoI sets differ".into());
                }
                if dimension(&a.rungs[i], &a.axis) != Some(d) {
                    return Err("joint resolution does not match its single-axis rung".into());
                }
            }
            for (j, q) in r.qoi.iter().enumerate() {
                let error = q.value - finest.qoi[j].value;
                let sum: f64 = axes
                    .iter()
                    .map(|a| a.rungs[i].qoi[j].value - a.rungs[joint.len() - 1].qoi[j].value)
                    .sum();
                if !error.is_finite() || !sum.is_finite() || !(error - sum).is_finite() {
                    return Err("nonfinite joint interaction".into());
                }
                interactions.push(Interaction {
                    rung: i,
                    joint_error: Observable {
                        name: q.name.clone(),
                        value: error,
                        unit: q.unit.clone(),
                    },
                    single_axis_sum: sum,
                    interaction: error - sum,
                });
            }
        }
        let mut report = Self {
            policy,
            axes,
            groups,
            joint: JointReport {
                rungs: joint,
                interactions,
            },
            ensemble,
            hash: ContentHash([0; 32]),
        };
        report.hash = hash_domain("org.frankensim.axis-ladder.v1", &report.canonical_bytes());
        Ok(report)
    }

    /// Canonical bytes of the retained report, excluding the hash itself.
    /// Rust Debug records are length-framed here; this is an internal v1
    /// transport, not a promised public interchange schema.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = b"frankensim.axis-ladder.internal.v1\0".to_vec();
        for s in [
            format!("{:?}", self.policy),
            format!("{:?}", self.axes),
            format!("{:?}", self.groups),
            format!("{:?}", self.joint),
            format!("{:?}", self.ensemble),
        ] {
            bytes.extend_from_slice(&(s.len() as u64).to_le_bytes());
            bytes.extend_from_slice(s.as_bytes());
        }
        bytes
    }
}
