//! Steady countercurrent extraction or stripping through ideal equilibrium stages.
//!
//! Aqueous flow visits stages `1..=N`; organic flow visits them in reverse.
//! With fixed positive volumetric flows `A` and `O`, every stage satisfies
//! `A*x[i-1] + O*y[i+1] = A*x[i] + O*y[i]` and `y[i] = F_i(x[i])`.
//! These are the single-solute, constant-phase-flow balances underlying the
//! McCabe–Thiele construction; see Srivastava, Werner and Honaker, *Mining*
//! 3 (2023), equation 4: <https://doi.org/10.3390/mining3030031>.
//!
//! Summing balances from stage i to N gives a backward shooting recurrence
//! `x[i-1] = r + (O/A)*(F_i(x[i]) - y_feed)`, starting with `x[N] = r`.
//! Each supplied isotherm is nondecreasing, so the predicted aqueous feed
//! increases strictly with r, even when a loading curve has flat segments.
//! Bounded bisection finds r without extrapolating any stage's original table.
//! A trial outside a table narrows the feasible bracket; it is never accepted
//! as a physical stream. Every reported stage and the complete train must meet
//! the caller's absolute/relative mol/s balance tolerance independently.
//!
//! Each stage retains its original isotherm, temperature and chemical context.
//! Its original O/A validity coordinate is evaluated at the actual phase-flow
//! ratio (equal to the contacted-volume ratio over any common time interval).
//! The reported flows are not stage hold-ups. Sources must justify applying
//! their equilibrium relation at each stage's conditions; this owner does not
//! infer changing acid, pH, extractant concentration or competing species.
//!
//! This is steady ideal-stage performance, not a residence-time, kinetic,
//! stage-efficiency, entrainment, heat, variable-volume, recycle, or complete
//! SX–EW plant model. One analytical element is conserved, without chemical
//! species or charge closure. Numerical residuals are not interval certificates
//! or experimental calibration. Fixed traversal/order is deterministic on one
//! ISA. Cancellation and exhausted budgets publish no partial solution.

use super::{
    AQUEOUS_CONCENTRATION_AXIS, ORGANIC_CONCENTRATION_AXIS, OWNED_AXES, SolventExtractionError,
    SolventExtractionIsotherm, TEMPERATURE_AXIS, VOLUME_RATIO_AXIS, product, represented, sum,
};
use fs_matdb::QueryPoint;
use fs_math::dd::Dd;
use fs_qty::{
    AmountConcentration, Dimensionless, Qty, QuantitySpec, Temperature, VolumetricFlowRate,
    chemistry::ElementId,
};
use std::sync::Arc;

/// Analytical-element flow in mol/s, not a mass flow or electrical current.
pub type AmountFlowRate = Qty<0, 0, -1, 0, 0, 1>;
/// Fixed memory and per-trial stage-count limit.
pub const MAX_STAGES: usize = 64;

/// Invalid train, original-source refusal, or an explicitly bounded solve failure.
#[derive(Clone, Debug, PartialEq)]
pub enum CountercurrentError {
    /// Original isotherm, unit, source-domain, or arithmetic refusal.
    Source(SolventExtractionError),
    /// Invalid flow, stage arrangement, or requested numerical allowance.
    Invalid(&'static str),
    /// No supported countercurrent solution was found on the original tables.
    OutsideIsotherms,
    /// Explicit root-evaluation or stage-evaluation budget was exhausted.
    Budget,
    /// Finite arithmetic cannot resolve the requested balance allowance.
    UnresolvedBalance,
    /// Cancellation was observed before publishing a complete solution.
    Cancelled,
}
impl core::fmt::Display for CountercurrentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "countercurrent solvent extraction: {self:?}")
    }
}
impl std::error::Error for CountercurrentError {}
impl From<SolventExtractionError> for CountercurrentError {
    fn from(error: SolventExtractionError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, CountercurrentError>;

/// Constant phase flow and its named analytical-element concentration.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseFlow {
    solute: ElementId,
    volume_flow: VolumetricFlowRate,
    concentration: AmountConcentration,
    amount_flow: AmountFlowRate,
}
impl PhaseFlow {
    /// Positive phase flow, nonnegative concentration, and representable mol/s.
    /// Feed concentrations need not themselves be equilibrium concentrations.
    pub fn new(
        solute: ElementId,
        volume_flow: VolumetricFlowRate,
        concentration: AmountConcentration,
    ) -> Result<Self> {
        require(
            volume_flow.value().is_finite() && volume_flow.value() > 0.0,
            "phase flow must be finite positive m3/s",
        )?;
        require(
            concentration.value().is_finite() && concentration.value() >= 0.0,
            "phase concentration must be finite nonnegative mol/m3",
        )?;
        let amount_flow =
            AmountFlowRate::new(product(volume_flow.value(), concentration.value())?.to_f64());
        Ok(Self {
            solute,
            volume_flow,
            concentration,
            amount_flow,
        })
    }
    pub fn solute(&self) -> &ElementId {
        &self.solute
    }
    pub const fn volume_flow(&self) -> VolumetricFlowRate {
        self.volume_flow
    }
    pub const fn concentration(&self) -> AmountConcentration {
        self.concentration
    }
    pub const fn amount_flow(&self) -> AmountFlowRate {
        self.amount_flow
    }
}

/// One ideal stage's immutable source and explicitly supplied chemical context.
#[derive(Clone, Debug, PartialEq)]
pub struct StageSpecification {
    isotherm: SolventExtractionIsotherm,
    temperature: Temperature,
    conditions: QueryPoint,
}
impl StageSpecification {
    pub fn new(
        isotherm: SolventExtractionIsotherm,
        temperature: Temperature,
        conditions: QueryPoint,
    ) -> Result<Self> {
        require(
            temperature.value().is_finite() && temperature.value() > 0.0,
            "stage temperature must be finite positive kelvin",
        )?;
        require(
            OWNED_AXES
                .iter()
                .all(|axis| !conditions.axes().contains_key(*axis)),
            "stage conditions cannot override owned equilibrium coordinates",
        )?;
        Ok(Self {
            isotherm,
            temperature,
            conditions,
        })
    }
    pub fn isotherm(&self) -> &SolventExtractionIsotherm {
        &self.isotherm
    }
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub fn conditions(&self) -> &QueryPoint {
        &self.conditions
    }
}

/// Explicit balance accuracy and bounded scalar work for one complete solve.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CascadeNumerics {
    /// Nonnegative mol/s allowance applied to each stage and the whole train.
    pub absolute_flow_tolerance: AmountFlowRate,
    /// Nonnegative fraction of the larger incoming/outgoing mol/s account.
    /// Must be below one; at least one of the two tolerances must be positive.
    pub relative_flow_tolerance: Dimensionless,
    /// Includes the two bracket endpoint trials and all interior trials.
    pub max_root_evaluations: usize,
    /// Includes admission, interpolation trials and final stage accounting.
    pub max_stage_evaluations: usize,
}
impl CascadeNumerics {
    fn validate(self) -> Result<()> {
        let abs = self.absolute_flow_tolerance.value();
        let rel = self.relative_flow_tolerance.value();
        require(
            abs.is_finite()
                && abs >= 0.0
                && rel.is_finite()
                && (0.0..1.0).contains(&rel)
                && (abs > 0.0 || rel > 0.0)
                && self.max_root_evaluations > 0
                && self.max_stage_evaluations > 0,
            "require finite nonnegative balance tolerances, positive accuracy and work budgets",
        )
    }
    fn admits(self, residual: f64, incoming: f64, outgoing: f64) -> bool {
        // Divide after the absolute allowance, avoiding overflow of abs + rel*scale.
        residual.is_finite()
            && (residual.abs() <= self.absolute_flow_tolerance.value()
                || (incoming.max(outgoing) > 0.0
                    && (residual.abs() - self.absolute_flow_tolerance.value())
                        / incoming.max(outgoing)
                        <= self.relative_flow_tolerance.value()))
    }
}

/// Actual inlet/outlet streams for one stage, ordered along aqueous flow.
#[derive(Clone, Debug, PartialEq)]
pub struct StageSolution {
    pub aqueous_inlet: PhaseFlow,
    pub organic_inlet: PhaseFlow,
    pub aqueous_outlet: PhaseFlow,
    pub organic_outlet: PhaseFlow,
    /// Positive for extraction, negative for stripping, in mol/s.
    pub transfer_to_organic: AmountFlowRate,
    /// Outgoing minus incoming analytical-element mol/s, independently summed.
    pub conservation_residual: AmountFlowRate,
    /// Lower-knot index, with the left segment chosen at an exact internal knot.
    pub segment_index: usize,
    /// Actual equilibrium coordinates and the unchanged additional source context.
    pub conditions: QueryPoint,
}

/// Complete retained-source solution; no stream is published before all checks pass.
#[derive(Clone, Debug, PartialEq)]
pub struct CountercurrentSolution {
    specifications: Arc<[StageSpecification]>,
    pub aqueous_feed: PhaseFlow,
    pub organic_feed: PhaseFlow,
    pub aqueous_outlet: PhaseFlow,
    pub organic_outlet: PhaseFlow,
    pub stages: Vec<StageSolution>,
    pub transfer_to_organic: AmountFlowRate,
    /// External outlets minus feeds; internal transfers are excluded.
    pub conservation_residual: AmountFlowRate,
    pub numerics: CascadeNumerics,
    pub root_evaluations: usize,
    pub stage_evaluations: usize,
}
impl CountercurrentSolution {
    pub fn specifications(&self) -> &[StageSpecification] {
        &self.specifications
    }
}

/// An ordered train of 1..=64 ideal contacts, all tracking the same element.
#[derive(Clone, Debug, PartialEq)]
pub struct CountercurrentCascade {
    specifications: Arc<[StageSpecification]>,
}
impl CountercurrentCascade {
    pub fn new(specifications: Vec<StageSpecification>) -> Result<Self> {
        require(
            !specifications.is_empty() && specifications.len() <= MAX_STAGES,
            "require 1..=64 countercurrent stages",
        )?;
        require(
            specifications
                .iter()
                .all(|s| s.isotherm.solute() == specifications[0].isotherm.solute()),
            "all stages must track the same analytical element",
        )?;
        Ok(Self {
            specifications: specifications.into(),
        })
    }
    pub fn specifications(&self) -> &[StageSpecification] {
        &self.specifications
    }

    /// Solve with aqueous feed entering stage 1 and organic feed entering stage N.
    /// Neither input, source cards nor this immutable train is mutated on failure.
    /// Cancellation is checked before every bounded stage evaluation.
    pub fn solve(
        &self,
        aqueous_feed: &PhaseFlow,
        organic_feed: &PhaseFlow,
        numerics: CascadeNumerics,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CountercurrentSolution> {
        numerics.validate()?;
        require(
            aqueous_feed.solute() == self.specifications[0].isotherm.solute()
                && organic_feed.solute() == aqueous_feed.solute(),
            "both feeds must track the isotherms' analytical element",
        )?;
        let ratio = represented(
            organic_feed.volume_flow.value() / aqueous_feed.volume_flow.value(),
            true,
        )?;
        let mut work = Work {
            numerics,
            roots: 0,
            stages: 0,
            cancelled: &mut cancelled,
        };
        for stage in self.specifications.iter() {
            work.stage()?;
            // Every table knot already lies inside its concentration bounds.
            // This first checks the fixed flow ratio, temperature and chemistry.
            let first = stage.isotherm.points()[0];
            actual_point(stage, first.aqueous.value(), first.organic.value(), ratio)?;
        }
        let last = self.specifications.last().expect("admitted nonempty train");
        let mut lower = last.isotherm.points()[0].aqueous.value();
        let mut upper = last.isotherm.points().last().unwrap().aqueous.value();
        let incoming = flow_sum(aqueous_feed, organic_feed)?;
        let mut lower_supported = false;
        let mut upper_supported = false;

        for endpoint in [lower, upper] {
            let trial = self.shoot(
                endpoint,
                organic_feed.concentration.value(),
                ratio,
                &mut work,
            )?;
            if let Trial::Supported(ref path) = trial {
                if self.root_admitted(path, aqueous_feed, organic_feed, incoming, numerics)? {
                    match self.finish(path, aqueous_feed, organic_feed, ratio, &mut work) {
                        // A smaller stage account can impose a tighter relative
                        // allowance than the external account. Keep refining.
                        Err(CountercurrentError::UnresolvedBalance) => {}
                        result => return result,
                    }
                }
            }
            let ordering = trial.ordering(aqueous_feed.concentration.value());
            if endpoint == lower {
                if ordering.is_gt() {
                    return Err(CountercurrentError::OutsideIsotherms);
                }
                lower_supported = matches!(trial, Trial::Supported(_));
            } else {
                if ordering.is_lt() {
                    return Err(CountercurrentError::OutsideIsotherms);
                }
                upper_supported = matches!(trial, Trial::Supported(_));
            }
        }
        loop {
            let middle = lower + (upper - lower) * 0.5;
            if middle == lower || middle == upper {
                return Err(if lower_supported && upper_supported {
                    CountercurrentError::UnresolvedBalance
                } else {
                    CountercurrentError::OutsideIsotherms
                });
            }
            let trial = self.shoot(middle, organic_feed.concentration.value(), ratio, &mut work)?;
            if let Trial::Supported(ref path) = trial {
                if self.root_admitted(path, aqueous_feed, organic_feed, incoming, numerics)? {
                    match self.finish(path, aqueous_feed, organic_feed, ratio, &mut work) {
                        Err(CountercurrentError::UnresolvedBalance) => {}
                        result => return result,
                    }
                }
            }
            if trial.ordering(aqueous_feed.concentration.value()).is_lt() {
                lower = middle;
                lower_supported = matches!(trial, Trial::Supported(_));
            } else {
                upper = middle;
                upper_supported = matches!(trial, Trial::Supported(_));
            }
        }
    }

    fn shoot(
        &self,
        raffinate: f64,
        organic_feed: f64,
        ratio: f64,
        work: &mut Work<'_>,
    ) -> Result<Trial> {
        work.root()?;
        let mut points = vec![
            RawPoint {
                aqueous: 0.0,
                organic: 0.0,
                segment: 0
            };
            self.specifications.len()
        ];
        let mut aqueous = Dd::from_f64(raffinate);
        for (index, stage) in self.specifications.iter().enumerate().rev() {
            work.stage()?;
            let knots = stage.isotherm.points();
            if aqueous.lt(Dd::from_f64(knots[0].aqueous.value())) {
                return Ok(Trial::Below);
            }
            if Dd::from_f64(knots.last().unwrap().aqueous.value()).lt(aqueous) {
                return Ok(Trial::Above);
            }
            let x = represented(aqueous.to_f64(), Dd::ZERO.lt(aqueous))?;
            let segment = knots
                .windows(2)
                .position(|p| x <= p[1].aqueous.value())
                .unwrap();
            let lo = knots[segment];
            let hi = knots[segment + 1];
            let fraction = (Dd::from_f64(x) - Dd::from_f64(lo.aqueous.value()))
                / (Dd::from_f64(hi.aqueous.value()) - Dd::from_f64(lo.aqueous.value()));
            let y_dd = Dd::from_f64(lo.organic.value())
                + fraction * (Dd::from_f64(hi.organic.value()) - Dd::from_f64(lo.organic.value()));
            let y = represented(y_dd.to_f64(), Dd::ZERO.lt(y_dd))?;
            points[index] = RawPoint {
                aqueous: x,
                organic: y,
                segment,
            };
            let difference = Dd::from_f64(y) - Dd::from_f64(organic_feed);
            // An unrepresentably distant trial can only be outside the next
            // nonnegative finite source interval; do not extrapolate through it.
            let shift = ratio * difference.to_f64();
            if shift == f64::INFINITY {
                return Ok(Trial::Above);
            }
            if shift == f64::NEG_INFINITY {
                return Ok(Trial::Below);
            }
            aqueous = Dd::from_f64(raffinate) + Dd::from_f64(ratio) * difference;
            if !aqueous.to_f64().is_finite() {
                return Err(CountercurrentError::UnresolvedBalance);
            }
        }
        Ok(Trial::Supported(Path {
            points,
            predicted_feed: aqueous,
        }))
    }

    fn root_admitted(
        &self,
        path: &Path,
        aqueous: &PhaseFlow,
        organic: &PhaseFlow,
        incoming: Dd,
        numerics: CascadeNumerics,
    ) -> Result<bool> {
        // Check the actual external account, not a concentration-only stopping rule.
        let outgoing = sum(
            Dd::from_f64(
                product(
                    aqueous.volume_flow.value(),
                    path.points.last().unwrap().aqueous,
                )?
                .to_f64(),
            ),
            Dd::from_f64(product(organic.volume_flow.value(), path.points[0].organic)?.to_f64()),
        )?;
        let residual = (outgoing - incoming).to_f64();
        Ok(numerics.admits(residual, incoming.to_f64(), outgoing.to_f64()))
    }

    fn finish(
        &self,
        path: &Path,
        aqueous_feed: &PhaseFlow,
        organic_feed: &PhaseFlow,
        ratio: f64,
        work: &mut Work<'_>,
    ) -> Result<CountercurrentSolution> {
        let phase = |volume_flow, concentration| {
            PhaseFlow::new(
                aqueous_feed.solute.clone(),
                volume_flow,
                AmountConcentration::new(concentration),
            )
        };
        let mut stages = Vec::with_capacity(path.points.len());
        for (i, point) in path.points.iter().enumerate() {
            work.stage()?;
            let aqueous_inlet = if i == 0 {
                aqueous_feed.clone()
            } else {
                phase(aqueous_feed.volume_flow, path.points[i - 1].aqueous)?
            };
            let organic_inlet = if i + 1 == path.points.len() {
                organic_feed.clone()
            } else {
                phase(organic_feed.volume_flow, path.points[i + 1].organic)?
            };
            let aqueous_outlet = phase(aqueous_feed.volume_flow, point.aqueous)?;
            let organic_outlet = phase(organic_feed.volume_flow, point.organic)?;
            let residual = account(
                &aqueous_inlet,
                &organic_inlet,
                &aqueous_outlet,
                &organic_outlet,
                work.numerics,
            )?;
            let conditions =
                actual_point(&self.specifications[i], point.aqueous, point.organic, ratio)?;
            let transfer_to_organic = AmountFlowRate::new(
                organic_outlet.amount_flow.value() - organic_inlet.amount_flow.value(),
            );
            stages.push(StageSolution {
                aqueous_inlet,
                organic_inlet,
                aqueous_outlet,
                organic_outlet,
                transfer_to_organic,
                conservation_residual: AmountFlowRate::new(residual),
                segment_index: point.segment,
                conditions,
            });
        }
        let aqueous_outlet = stages.last().unwrap().aqueous_outlet.clone();
        let organic_outlet = stages[0].organic_outlet.clone();
        let residual = account(
            aqueous_feed,
            organic_feed,
            &aqueous_outlet,
            &organic_outlet,
            work.numerics,
        )?;
        work.check_cancelled()?;
        Ok(CountercurrentSolution {
            specifications: self.specifications.clone(),
            aqueous_feed: aqueous_feed.clone(),
            organic_feed: organic_feed.clone(),
            transfer_to_organic: AmountFlowRate::new(
                organic_outlet.amount_flow.value() - organic_feed.amount_flow.value(),
            ),
            aqueous_outlet,
            organic_outlet,
            stages,
            conservation_residual: AmountFlowRate::new(residual),
            numerics: work.numerics,
            root_evaluations: work.roots,
            stage_evaluations: work.stages,
        })
    }
}

#[derive(Clone, Copy)]
struct RawPoint {
    aqueous: f64,
    organic: f64,
    segment: usize,
}
struct Path {
    points: Vec<RawPoint>,
    predicted_feed: Dd,
}
enum Trial {
    Below,
    Above,
    Supported(Path),
}
impl Trial {
    fn ordering(&self, feed: f64) -> core::cmp::Ordering {
        match self {
            Self::Below => core::cmp::Ordering::Less,
            Self::Above => core::cmp::Ordering::Greater,
            Self::Supported(path) => {
                if path.predicted_feed.lt(Dd::from_f64(feed)) {
                    core::cmp::Ordering::Less
                } else if Dd::from_f64(feed).lt(path.predicted_feed) {
                    core::cmp::Ordering::Greater
                } else {
                    core::cmp::Ordering::Equal
                }
            }
        }
    }
}
struct Work<'a> {
    numerics: CascadeNumerics,
    roots: usize,
    stages: usize,
    cancelled: &'a mut dyn FnMut() -> bool,
}
impl Work<'_> {
    fn check_cancelled(&mut self) -> Result<()> {
        if (self.cancelled)() {
            Err(CountercurrentError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn root(&mut self) -> Result<()> {
        self.check_cancelled()?;
        if self.roots >= self.numerics.max_root_evaluations {
            return Err(CountercurrentError::Budget);
        }
        self.roots += 1;
        Ok(())
    }
    fn stage(&mut self) -> Result<()> {
        self.check_cancelled()?;
        if self.stages >= self.numerics.max_stage_evaluations {
            return Err(CountercurrentError::Budget);
        }
        self.stages += 1;
        Ok(())
    }
}
fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(CountercurrentError::Invalid(message))
    }
}
fn actual_point(stage: &StageSpecification, x: f64, y: f64, ratio: f64) -> Result<QueryPoint> {
    let mut point = stage.conditions.clone();
    for (axis, dims, value) in [
        (AQUEOUS_CONCENTRATION_AXIS, AmountConcentration::DIMS, x),
        (ORGANIC_CONCENTRATION_AXIS, AmountConcentration::DIMS, y),
        (VOLUME_RATIO_AXIS, Dimensionless::DIMS, ratio),
        (
            TEMPERATURE_AXIS,
            Temperature::DIMS,
            stage.temperature.value(),
        ),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(dims), value)
            .map_err(SolventExtractionError::from)?;
    }
    if !stage
        .isotherm
        .source()
        .validity
        .contains_typed(point.axes(), point.axis_quantities())
    {
        return Err(SolventExtractionError::OutsideValidity.into());
    }
    Ok(point)
}
fn flow_sum(aqueous: &PhaseFlow, organic: &PhaseFlow) -> Result<Dd> {
    Ok(sum(
        Dd::from_f64(aqueous.amount_flow.value()),
        Dd::from_f64(organic.amount_flow.value()),
    )?)
}
fn account(
    ai: &PhaseFlow,
    oi: &PhaseFlow,
    ao: &PhaseFlow,
    oo: &PhaseFlow,
    numerics: CascadeNumerics,
) -> Result<f64> {
    let incoming = sum(
        Dd::from_f64(ai.amount_flow.value()),
        Dd::from_f64(oi.amount_flow.value()),
    )?;
    let outgoing = sum(
        Dd::from_f64(ao.amount_flow.value()),
        Dd::from_f64(oo.amount_flow.value()),
    )?;
    let residual = (outgoing - incoming).to_f64();
    if !numerics.admits(residual, incoming.to_f64(), outgoing.to_f64()) {
        return Err(CountercurrentError::UnresolvedBalance);
    }
    Ok(residual)
}

#[cfg(test)]
mod tests;
