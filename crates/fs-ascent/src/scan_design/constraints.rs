//! Constrained scan design using the original physical objective and native SQP.
//!
//! Limits are explicit inequalities, not clipped sources or fictitious objective
//! penalties. SQP may evaluate infeasible designs; original source/material/solve
//! failures still propagate. A successful simulation or exhausted optimizer is
//! not a feasible optimum: inspect the returned KKT residuals and physical rows.

use super::{ScanGeometryControl,ScanPowerStudy,ScanPowerEvaluation,ScanSensitivityError,
    MovingHeatSource,ScanningEnthalpyRun,ScanObservation,ScanObservationSpec,finite,invalid,poll};
use crate::{SqpState,SqpSample,SqpError,SqpRunReport};
use fs_exec::Cx;
use std::collections::{BTreeMap,BTreeSet};
pub use fs_conduction::bc::history::march::enthalpy::scanning::sensitivity::objective::responses::{
    ScanResponse, ScanResponseBudget, peaks::ScanPeakResponse};

/// A bound on one selected control in its original physical units, not its
/// transformed optimizer coordinate. Powers come first [W], then geometry [m].
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ScanControlBound {
    pub coordinate: usize,
    pub lower: Option<f64>,
    pub upper: Option<f64>,
    /// Positive W or m scale. A feasibility tolerance tau permits tau*scale
    /// physical violation. It is not a numerical solver or measurement error.
    pub scale: f64,
}

/// Total authored domain-absorbed energy of ALL attached moving sources over
/// this absolute-time window. Unselected power values are included. Background
/// volumetric heating, radiation and electricity consumed by a heater are not.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ScanEnergyCap {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub maximum_j: f64,
    pub scale_j: f64,
}

/// A lower/upper limit on actual temperature [K] or specific enthalpy [J/kg]
/// at one fixed-grid endpoint and mesh vertex. These are separate constraints,
/// not fitting targets. One response sweep supplies both sides of an interval.
/// No claim is made about unsampled vertices, times or a continuous maximum.
#[derive(Debug,Clone,Copy,PartialEq)]
pub struct ScanStateBound {
    pub response: ScanResponse,
    pub lower: Option<f64>,
    pub upper: Option<f64>,
    /// Same physical units as the response. SQP feasibility tolerance tau
    /// permits tau*scale violation, independently of the primal solve error.
    pub scale: f64,
}

/// Upper limit on a smooth envelope of selected nodal temperature/enthalpy
/// samples. Unlike bounding their mean, this constrains every selected sample
/// in exact arithmetic. The explicit smoothing bias makes it conservative,
/// not an interval/continuum certificate. A lower bound on this envelope would
/// NOT establish a lower bound on every sample and is intentionally absent.
#[derive(Debug,Clone,PartialEq)]
pub struct ScanPeakBound {
    pub response: ScanPeakResponse,
    pub maximum: f64,
    /// Same K or J/kg units as the response; SQP tolerance tau permits a
    /// tau*scale violation independently of smoothing and physical solve error.
    pub scale: f64,
}

/// Constraint observations and derivatives at exactly one design point.
#[derive(Debug,Clone,PartialEq)]
pub struct ScanConstraintEvaluation {
    /// Selected physical watts for automatically inserted nonnegative-domain
    /// rows. Empty under the original logarithmic map. These rows are first.
    pub nonnegative_power_values_w: Vec<f64>,
    /// One value in watts or metres per declared bound, in declaration order.
    pub bound_values: Vec<f64>,
    /// Joules per declared energy window, including fixed source contributions.
    pub absorbed_energy_j: Vec<f64>,
    /// Actual kelvin/J-per-kg responses in state-bound declaration order.
    pub state_values: Vec<f64>,
    /// Smooth upper envelopes, in peak-bound declaration order.
    pub peak_values: Vec<f64>,
    /// Actual selected maxima from the same trajectories, before smoothing.
    pub sampled_peak_values: Vec<f64>,
    /// Added smoothing margins in physical units, not solver-error estimates.
    pub peak_smoothing_gaps: Vec<f64>,
    /// <= 0 is feasible. Nonnegative-domain rows, control lower/upper rows, energy rows, then state
    /// lower/upper rows, then peak upper rows, in original declaration order.
    pub residuals: Vec<f64>,
    /// Row-major derivatives with respect to original optimizer coordinates.
    pub jacobian: Vec<f64>,
    /// Additional response reverse work, excluding the fitting-objective sweep.
    pub state_transpose_iterations: usize,
    pub max_state_transpose_relative_residual: f64,
}
impl ScanConstraintEvaluation {
    pub fn max_scaled_violation(&self) -> f64 {
        self.residuals.iter().map(|v|v.max(0.0)).fold(0.0,f64::max)
    }
}

struct EnergyRow { cap: ScanEnergyCap, fixed_j: f64, power_weights_s: Vec<f64> }
struct Limits { bounds: Vec<ScanControlBound>, energy: Vec<EnergyRow>,
    states: Vec<ScanStateBound>, peaks: Vec<ScanPeakBound>, max_peak_samples: usize,
    response_budget: ScanResponseBudget, rows: usize }

/// One study, immutable constraints and its own accepted SQP checkpoint.
/// This owner prevents resuming a checkpoint with a different constraint set.
/// The native dense SQP is for bounded small design spaces, not large sparse
/// constrained optimization. All physical evaluations retain the study budget.
pub struct ConstrainedScanStudy<F, O = ScanObservation> {
    study: ScanPowerStudy<F, O>,
    limits: Limits,
    state: SqpState,
}

impl<F, O: ScanObservationSpec> ConstrainedScanStudy<F, O> {
    /// Admit all constraint shapes before allocating a dense Jacobian or doing
    /// physics, then evaluate the original initial design once. Starts may be
    /// constraint-infeasible but must be physically evaluable. Equal lower/upper
    /// bounds are refused: keep a fixed coordinate unselected instead.
    /// All energy windows require constant/ramp/held-or-linear power schedules.
    /// Nonnegative-linear powers automatically add one `-x <= 0` domain row
    /// per selected power, counted in the KKT cap. No user bounds are needed
    /// just to admit the off boundary. Other design limits may remain infeasible
    /// during the search; negative heater inputs never enter a physics solve.
    pub fn new<'h,'m,'c>(cx: &Cx<'_>, study: ScanPowerStudy<F, O>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, max_kkt_dimension: usize)
        -> Result<Self,SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        Self::new_with_state_bounds(cx,study,bounds,energy,Vec::new(),
            ScanResponseBudget {max_responses:0,max_gradient_entries:0,max_transpose_solves:0},max_kkt_dimension)
    }

    /// Add limits on the solved trajectory, sharing each original forward run
    /// with the fitting objective. All response gradients pass through full
    /// previous enthalpy, nonlinear conductivity and optional radiation.
    /// State-only constrained studies are valid. Bounds stay immutable across
    /// SQP continuation, and every extra reverse sweep has an explicit budget.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_state_bounds<'h,'m,'c>(cx: &Cx<'_>, study: ScanPowerStudy<F, O>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>,
        response_budget: ScanResponseBudget, max_kkt_dimension: usize)
        -> Result<Self,SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        Self::new_with_peak_bounds(cx,study,bounds,energy,states,Vec::new(),
            response_budget,0,max_kkt_dimension)
    }

    /// Constrain hot spots across whole selected regions and time windows with
    /// one response sweep per peak. Ordinary point/input constraints can coexist.
    /// All rows share one original forward trajectory. Each peak adds just one
    /// dense SQP row irrespective of its number of nodal/time samples.
    ///
    /// Combined point/peak reverse work uses response_budget; max_peak_samples
    /// bounds total selected vertex/time pairs. Peak selection, smoothing and
    /// limits stay immutable across optimizer continuation and physical replay.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_peak_bounds<'h,'m,'c>(cx: &Cx<'_>, mut study: ScanPowerStudy<F, O>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>,
        peaks: Vec<ScanPeakBound>, response_budget: ScanResponseBudget,
        max_peak_samples: usize, max_kkt_dimension: usize)
        -> Result<Self,SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        let limits=Limits::new(cx,&study,bounds,energy,states,peaks,response_budget,max_peak_samples,max_kkt_dimension)
            .map_err(SqpError::Evaluation)?;
        let initial=study.initial_coordinates().to_vec();
        let state=SqpState::try_new(&initial,max_kkt_dimension,
            &mut |point| limits.sample(cx,&mut study,point),Some(cx))?;
        Ok(Self {study,limits,state})
    }

    pub fn point(&self) -> &[f64] { self.state.point() }
    pub fn sample(&self) -> &SqpSample { self.state.sample() }
    pub fn evaluations(&self) -> usize { self.study.evaluations() }
    pub fn evaluation_limit(&self) -> usize { self.study.evaluation_limit() }
    pub fn iterations(&self) -> usize { self.state.iterations() }

    /// Check input-only constraints without physics. With state bounds this
    /// refuses rather than silently omitting them; use evaluate_at or
    /// replay_with_constraints for all physical response values and gradients.
    pub fn constraints_at(&self,cx: &Cx<'_>,point: &[f64])
        -> Result<ScanConstraintEvaluation,ScanSensitivityError>
    {
        poll(cx)?;
        if !self.limits.states.is_empty() || !self.limits.peaks.is_empty() {
            return Err(invalid("state constraints need a physical evaluation; use evaluate_at or replay_with_constraints"));
        }
        self.limits.at(cx,&self.study,point)
    }

    /// Original immutable sampled-state requirements, in report order.
    pub fn state_bounds(&self) -> &[ScanStateBound] { &self.limits.states }

    /// Original peak selections and physical smoothing choices, in report order.
    pub fn peak_bounds(&self) -> &[ScanPeakBound] { &self.limits.peaks }

    /// Evaluate the objective AND all constraints at a supplied candidate.
    /// Costs one physical evaluation plus the admitted response sweeps; never
    /// changes the accepted SQP point, multipliers or curvature history.
    pub fn evaluate_at<'h,'m,'c>(&mut self,cx: &Cx<'_>,point: &[f64])
        -> Result<(ScanPowerEvaluation,ScanConstraintEvaluation),ScanSensitivityError>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    { self.limits.evaluate(cx,&mut self.study,point) }

    /// Re-run the accepted design with every attached state requirement. The
    /// values come from that same trajectory, not reconstructed penalty terms.
    pub fn replay_with_constraints<'h,'m,'c>(&mut self,cx: &Cx<'_>)
        -> Result<(ScanPowerEvaluation,ScanConstraintEvaluation),ScanSensitivityError>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    { self.limits.evaluate(cx,&mut self.study,self.state.point()) }

    /// Export the current candidate. It may be infeasible or unconverged after
    /// an interrupted run; export alone makes no constraint-satisfaction claim.
    pub fn sources(&self,cx: &Cx<'_>) -> Result<Vec<MovingHeatSource>,ScanSensitivityError> {
        self.study.sources_at(cx,self.state.point())
    }

    pub fn extend_evaluation_budget(&mut self,total: usize) -> Result<(),ScanSensitivityError> {
        self.study.extend_evaluation_budget(total)
    }

    /// Continue native SQP. Reports retain the original objective, inequality
    /// multipliers and all four KKT residuals. A stalled search is not an
    /// infeasibility certificate; a feasible optimum need not match the targets.
    /// Cancellation/failure retains the accepted state while counting spent work.
    pub fn advance_optimizer<'h,'m,'c>(&mut self,cx: &Cx<'_>,tolerance: f64,additional_iterations: usize)
        -> Result<SqpRunReport,SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        let remaining=self.study.evaluation_limit().saturating_sub(self.study.evaluations());
        let ceiling=self.state.evaluations().checked_add(remaining)
            .ok_or(SqpError::Invalid("constrained scan evaluation ceiling overflow"))?;
        let Self {study,limits,state}=self;
        state.try_run(&mut |point|limits.sample(cx,study,point),tolerance,
            additional_iterations,ceiling,Some(cx))
    }

    /// Re-run the accepted candidate through the original objective/adjoint,
    /// consuming one physical evaluation. This does not reset the SQP checkpoint.
    pub fn replay<'h,'m,'c>(&mut self,cx: &Cx<'_>) -> Result<ScanPowerEvaluation,ScanSensitivityError>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    { self.replay_with_constraints(cx).map(|(physical,_)|physical) }
}

impl Limits {
    #[allow(clippy::too_many_arguments)]
    fn new<F, O: ScanObservationSpec>(cx: &Cx<'_>, study: &ScanPowerStudy<F, O>, bounds: Vec<ScanControlBound>,
        caps: Vec<ScanEnergyCap>, states: Vec<ScanStateBound>, peaks: Vec<ScanPeakBound>,
        response_budget: ScanResponseBudget, max_peak_samples: usize, maximum: usize)
        -> Result<Self,ScanSensitivityError>
    {
        poll(cx)?;
        let n=study.initial.len();
        let mut rows=caps.len().checked_add(if study.has_nonnegative_powers() {study.controls.len()} else {0})
            .ok_or_else(||invalid("scan constraint count overflow"))?;
        let mut seen=BTreeSet::new();
        if states.len().checked_add(peaks.len()).is_none_or(|n|n>response_budget.max_responses) {
            return Err(invalid("state bounds exceed their response allowance"));
        }
        let mut response_solves=study.objective.endpoints_s.len();
        for state in &states {
            poll(cx)?;
            if state.response.endpoint>=study.objective.endpoints_s.len()
                || !state.scale.is_finite() || state.scale<=0.0
                || (state.lower.is_none() && state.upper.is_none())
                || state.lower.into_iter().chain(state.upper).any(|v|!v.is_finite())
                || state.lower.zip(state.upper).is_some_and(|(a,b)|a>=b)
            { return Err(invalid("state bounds require an existing endpoint, finite ordered limits and a positive physical scale")); }
            rows=rows.checked_add(usize::from(state.lower.is_some())+usize::from(state.upper.is_some()))
                .ok_or_else(||invalid("state constraint count overflow"))?;
            response_solves=response_solves.checked_add(state.response.endpoint+1)
                .ok_or_else(||invalid("state response work count overflow"))?;
        }
        let mut peak_samples=0usize;
        for peak in &peaks {
            poll(cx)?;
            if !peak.maximum.is_finite() || !peak.scale.is_finite() || peak.scale<=0.0 {
                return Err(invalid("scan peak bounds need a finite physical maximum and positive scale"));
            }
            // The factory supplies the actual mesh. The physical evaluator
            // additionally validates every vertex ID before the forward solve.
            peak_samples+=peak.response.validate(cx,study.objective.endpoints_s.len(),
                usize::MAX,max_peak_samples-peak_samples)?;
            rows=rows.checked_add(1).ok_or_else(||invalid("scan peak row count overflow"))?;
            response_solves=response_solves.checked_add(peak.response.last_endpoint+1)
                .ok_or_else(||invalid("scan peak response work count overflow"))?;
        }
        if (!states.is_empty() || !peaks.is_empty()) && response_solves>response_budget.max_transpose_solves {
            return Err(invalid("state bounds and objective exceed their transpose allowance"));
        }
        for bound in &bounds {
            poll(cx)?;
            if bound.coordinate>=n || !seen.insert(bound.coordinate) || !bound.scale.is_finite()
                || bound.scale<=0.0 || (bound.lower.is_none() && bound.upper.is_none())
                || bound.lower.into_iter().chain(bound.upper).any(|v|!v.is_finite())
                || bound.lower.zip(bound.upper).is_some_and(|(a,b)|a>=b)
            { return Err(invalid("scan bounds need distinct existing coordinates, finite ordered limits and positive physical scales")); }
            let positive=match bound.coordinate.checked_sub(study.controls.len()) {
                None=>!study.has_nonnegative_powers(),
                Some(index)=>matches!(study.geometry_controls[index],ScanGeometryControl::HalfWidth {..}),
            };
            if positive && (bound.lower.is_some_and(|v|v<0.0) || bound.upper.is_some_and(|v|v<=0.0)) {
                return Err(invalid("power/half-width bounds must admit strictly positive physical values"));
            }
            if bound.coordinate < study.controls.len() && study.has_nonnegative_powers()
                && bound.lower.into_iter().chain(bound.upper).any(|v| v < 0.0) {
                return Err(invalid("nonnegative heater bounds cannot request negative watts"));
            }
            rows=rows.checked_add(usize::from(bound.lower.is_some())+usize::from(bound.upper.is_some()))
                .ok_or_else(||invalid("scan constraint count overflow"))?;
        }
        let dim=n.checked_add(rows).filter(|&d|d<=maximum)
            .ok_or_else(||invalid("scan constraints exceed their explicit dense KKT dimension allowance"))?;
        if rows==0 || dim.checked_mul(dim).is_none() {return Err(invalid("constrained scan needs bounded nonempty constraints"));}
        let selected: BTreeMap<_,_>=study.controls.iter().enumerate().map(|(i,c)|((c.source,c.value),i)).collect();
        let mut energy=Vec::with_capacity(caps.len());
        for cap in caps {
            poll(cx)?;
            if !cap.start_time_s.is_finite() || !cap.end_time_s.is_finite() || cap.end_time_s<=cap.start_time_s
                || !(cap.end_time_s-cap.start_time_s).is_finite() || !cap.maximum_j.is_finite()
                || cap.maximum_j<0.0 || !cap.scale_j.is_finite() || cap.scale_j<=0.0 {
                return Err(invalid("scan energy caps need finite positive-duration windows, nonnegative joules and positive joule scales"));
            }
            let mut row=EnergyRow {cap,fixed_j:0.0,power_weights_s:vec![0.0;study.controls.len()]};
            for (source_index,source) in study.sources.iter().enumerate() {
                poll(cx)?;
                let integral=source.integrate_absorbed_power(cx,cap.start_time_s,cap.end_time_s,study.objective.budget.max_schedule_controls)?;
                let values=source.absorbed_power_values(study.objective.budget.max_schedule_controls)?;
                // Sum fixed terms directly: subtracting selected energy from a
                // large total could erase a small but real fixed heat supply.
                for (value_index,(&power,&weight)) in values.iter().zip(&integral.power_value_weights_s).enumerate() {
                    if value_index%256==0 {poll(cx)?;}
                    if let Some(&index)=selected.get(&(source_index,value_index)) {row.power_weights_s[index]=weight;}
                    else {row.fixed_j=finite(row.fixed_j+finite(power*weight)?)?;}
                }
            }
            energy.push(row);
        }
        poll(cx)?;
        Ok(Self {bounds,energy,states,peaks,max_peak_samples,response_budget,rows})
    }

    fn at<F, O: ScanObservationSpec>(&self,cx: &Cx<'_>,study: &ScanPowerStudy<F, O>,point: &[f64])
        -> Result<ScanConstraintEvaluation,ScanSensitivityError>
    {
        poll(cx)?;
        let powers=study.powers_at(point)?;let geometry=study.geometry_at(point)?;let n=point.len();
        let mut out=ScanConstraintEvaluation {nonnegative_power_values_w:Vec::new(),
            bound_values:Vec::with_capacity(self.bounds.len()),
            absorbed_energy_j:Vec::with_capacity(self.energy.len()),residuals:Vec::with_capacity(self.rows),
            state_values:Vec::with_capacity(self.states.len()),jacobian:vec![0.0;self.rows*n],
            peak_values:Vec::with_capacity(self.peaks.len()),sampled_peak_values:Vec::with_capacity(self.peaks.len()),
            peak_smoothing_gaps:Vec::with_capacity(self.peaks.len()),
            state_transpose_iterations:0,max_state_transpose_relative_residual:0.0};
        if study.has_nonnegative_powers() {
            out.nonnegative_power_values_w = powers.clone();
            for (i, &x) in point[..powers.len()].iter().enumerate() {
                if i % 256 == 0 { poll(cx)?; }
                out.jacobian[out.residuals.len()*n+i] = -1.0;
                out.residuals.push(-x);
            }
        }
        for b in &self.bounds {
            poll(cx)?;
            let (value,chain)=if b.coordinate<powers.len() {
                (powers[b.coordinate],study.power_derivative(b.coordinate,powers[b.coordinate]))}
                else {let index=b.coordinate-powers.len();let value=geometry[index];
                    (value,match study.geometry_controls[index] {
                        ScanGeometryControl::CenterValue {scale_m,..}=>scale_m,
                        ScanGeometryControl::HalfWidth {..}=>value})};
            out.bound_values.push(value);
            for (limit,sign) in [(b.lower,-1.0),(b.upper,1.0)] {
                if let Some(limit)=limit {
                    out.jacobian[out.residuals.len()*n+b.coordinate]=finite(sign*(chain/b.scale))?;
                    out.residuals.push(finite(sign*finite(value-limit)?/b.scale)?);
                }
            }
        }
        for row in &self.energy {
            poll(cx)?;
            let mut joules=row.fixed_j;let offset=out.residuals.len()*n;
            for (i,(&weight,&power)) in row.power_weights_s.iter().zip(&powers).enumerate() {
                if i%256==0 {poll(cx)?;}
                let term=finite(weight*power)?;
                joules=finite(joules+term)?;
                out.jacobian[offset+i]=finite(finite(weight*study.power_derivative(i,power))?/row.cap.scale_j)?;
            }
            out.absorbed_energy_j.push(joules);
            out.residuals.push(finite(finite(joules-row.cap.maximum_j)?/row.cap.scale_j)?);
        }
        poll(cx)?;
        Ok(out)
    }

    fn sample<'h,'m,'c,F,O: ScanObservationSpec>(&self,cx: &Cx<'_>,study: &mut ScanPowerStudy<F, O>,point: &[f64])
        -> Result<Option<SqpSample>,ScanSensitivityError>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        poll(cx)?;
        if !study.power_domain_contains(point)? {
            // Only this analytic domain refusal becomes an unavailable SQP
            // trial. Phase, source, arithmetic and solver errors still propagate.
            // The SQP callback attempt is counted, but no physical reservation
            // or source reconstruction is spent on negative heater power.
            return Ok(None);
        }
        let (physical,constraints)=self.evaluate(cx,study,point)?;
        Ok(Some(SqpSample {f:physical.physical.value,gradient:physical.gradient,
            ce:Vec::new(),je:Vec::new(),ci:constraints.residuals,ji:constraints.jacobian}))
    }

    fn evaluate<'h,'m,'c,F,O: ScanObservationSpec>(&self,cx: &Cx<'_>,study: &mut ScanPowerStudy<F, O>,point: &[f64])
        -> Result<(ScanPowerEvaluation,ScanConstraintEvaluation),ScanSensitivityError>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        let mut out=self.at(cx,study,point)?;
        let requests:Vec<_>=self.states.iter().map(|state|state.response).collect();
        let peaks:Vec<_>=self.peaks.iter().map(|peak|peak.response.clone()).collect();
        let (physical,responses,peak_responses)=study.evaluate_with_responses(cx,point,
            if requests.is_empty() && peaks.is_empty() {None}
            else {Some((&requests,&peaks,self.response_budget,self.max_peak_samples))})?;
        if responses.len()!=self.states.len() {return Err(invalid("state response layout changed"));}
        for (state,response) in self.states.iter().zip(responses) {
            poll(cx)?;
            if response.response!=state.response {return Err(invalid("state response identity changed"));}
            let gradient=study.chain_gradient(cx,&physical.powers_w,&physical.geometry_m,&response.sources)?;
            out.state_values.push(response.value);
            out.state_transpose_iterations=out.state_transpose_iterations.checked_add(response.transpose_iterations)
                .ok_or_else(||invalid("state response iteration count overflow"))?;
            out.max_state_transpose_relative_residual=out.max_state_transpose_relative_residual
                .max(response.max_transpose_relative_residual);
            for (limit,sign) in [(state.lower,-1.0),(state.upper,1.0)] {
                if let Some(limit)=limit {
                    let offset=out.residuals.len()*point.len();
                    for (i,&g) in gradient.iter().enumerate() {
                        if i%256==0 {poll(cx)?;}
                        out.jacobian[offset+i]=finite(sign*(g/state.scale))?;
                    }
                    out.residuals.push(finite(sign*finite(response.value-limit)?/state.scale)?);
                }
            }
        }
        if peak_responses.len()!=self.peaks.len() {return Err(invalid("scan peak response layout changed"));}
        for (bound,response) in self.peaks.iter().zip(peak_responses) {
            poll(cx)?;
            if response.peak!=bound.response {return Err(invalid("scan peak response identity changed"));}
            let gradient=study.chain_gradient(cx,&physical.powers_w,&physical.geometry_m,&response.sources)?;
            out.peak_values.push(response.value);
            out.sampled_peak_values.push(response.sampled_maximum);
            out.peak_smoothing_gaps.push(response.smoothing_gap);
            out.state_transpose_iterations=out.state_transpose_iterations.checked_add(response.transpose_iterations)
                .ok_or_else(||invalid("scan peak iteration count overflow"))?;
            out.max_state_transpose_relative_residual=out.max_state_transpose_relative_residual
                .max(response.max_transpose_relative_residual);
            let offset=out.residuals.len()*point.len();
            for (i,&g) in gradient.iter().enumerate() {
                if i%256==0 {poll(cx)?;}
                out.jacobian[offset+i]=finite(g/bound.scale)?;
            }
            out.residuals.push(finite(finite(response.value-bound.maximum)?/bound.scale)?);
        }
        poll(cx)?;
        Ok((physical,out))
    }
}
