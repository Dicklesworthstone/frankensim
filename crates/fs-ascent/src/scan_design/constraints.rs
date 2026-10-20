//! Constrained scan design using the original physical objective and native SQP.
//!
//! Limits are explicit inequalities, not clipped sources or fictitious objective
//! penalties. SQP may evaluate infeasible designs; original source/material/solve
//! failures still propagate. A successful simulation or exhausted optimizer is
//! not a feasible optimum: inspect the returned KKT residuals and physical rows.

use super::{ScanGeometryControl,ScanPowerStudy,ScanPowerEvaluation,ScanSensitivityError,
    MovingHeatSource,ScanningEnthalpyRun,finite,invalid,poll};
use crate::{SqpState,SqpSample,SqpError,SqpRunReport};
use fs_exec::Cx;
use std::collections::{BTreeMap,BTreeSet};

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

/// Constraint observations and derivatives at exactly one design point.
#[derive(Debug,Clone,PartialEq)]
pub struct ScanConstraintEvaluation {
    /// One value in watts or metres per declared bound, in declaration order.
    pub bound_values: Vec<f64>,
    /// Joules per declared energy window, including fixed source contributions.
    pub absorbed_energy_j: Vec<f64>,
    /// <= 0 is feasible. Bound lower then upper rows, followed by energy rows.
    pub residuals: Vec<f64>,
    /// Row-major derivatives with respect to original optimizer coordinates.
    pub jacobian: Vec<f64>,
}
impl ScanConstraintEvaluation {
    pub fn max_scaled_violation(&self) -> f64 {
        self.residuals.iter().map(|v|v.max(0.0)).fold(0.0,f64::max)
    }
}

struct EnergyRow { cap: ScanEnergyCap, fixed_j: f64, power_weights_s: Vec<f64> }
struct Limits { bounds: Vec<ScanControlBound>, energy: Vec<EnergyRow>, rows: usize }

/// One study, immutable constraints and its own accepted SQP checkpoint.
/// This owner prevents resuming a checkpoint with a different constraint set.
/// The native dense SQP is for bounded small design spaces, not large sparse
/// constrained optimization. All physical evaluations retain the study budget.
pub struct ConstrainedScanStudy<F> {
    study: ScanPowerStudy<F>,
    limits: Limits,
    state: SqpState,
}

impl<F> ConstrainedScanStudy<F> {
    /// Admit all constraint shapes before allocating a dense Jacobian or doing
    /// physics, then evaluate the original initial design once. Starts may be
    /// constraint-infeasible but must be physically evaluable. Equal lower/upper
    /// bounds are refused: keep a fixed coordinate unselected instead.
    /// All energy windows require constant/ramp/held-or-linear power schedules.
    pub fn new<'h,'m,'c>(cx: &Cx<'_>, mut study: ScanPowerStudy<F>,
        bounds: Vec<ScanControlBound>, energy: Vec<ScanEnergyCap>, max_kkt_dimension: usize)
        -> Result<Self,SqpError<ScanSensitivityError>>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        let limits=Limits::new(cx,&study,bounds,energy,max_kkt_dimension).map_err(SqpError::Evaluation)?;
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

    /// Check constraints without a physics evaluation or spent-budget change.
    pub fn constraints_at(&self,cx: &Cx<'_>,point: &[f64])
        -> Result<ScanConstraintEvaluation,ScanSensitivityError>
    { self.limits.at(cx,&self.study,point) }

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
    { self.study.evaluate(cx,self.state.point()) }
}

impl Limits {
    fn new<F>(cx: &Cx<'_>, study: &ScanPowerStudy<F>, bounds: Vec<ScanControlBound>,
        caps: Vec<ScanEnergyCap>, maximum: usize) -> Result<Self,ScanSensitivityError>
    {
        poll(cx)?;
        let n=study.initial.len();let mut rows=caps.len();let mut seen=BTreeSet::new();
        for bound in &bounds {
            poll(cx)?;
            if bound.coordinate>=n || !seen.insert(bound.coordinate) || !bound.scale.is_finite()
                || bound.scale<=0.0 || (bound.lower.is_none() && bound.upper.is_none())
                || bound.lower.into_iter().chain(bound.upper).any(|v|!v.is_finite())
                || bound.lower.zip(bound.upper).is_some_and(|(a,b)|a>=b)
            { return Err(invalid("scan bounds need distinct existing coordinates, finite ordered limits and positive physical scales")); }
            let positive=match bound.coordinate.checked_sub(study.controls.len()) {
                None=>true,
                Some(index)=>matches!(study.geometry_controls[index],ScanGeometryControl::HalfWidth {..}),
            };
            if positive && (bound.lower.is_some_and(|v|v<0.0) || bound.upper.is_some_and(|v|v<=0.0)) {
                return Err(invalid("power/half-width bounds must admit strictly positive physical values"));
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
        Ok(Self {bounds,energy,rows})
    }

    fn at<F>(&self,cx: &Cx<'_>,study: &ScanPowerStudy<F>,point: &[f64])
        -> Result<ScanConstraintEvaluation,ScanSensitivityError>
    {
        poll(cx)?;
        let powers=study.powers_at(point)?;let geometry=study.geometry_at(point)?;let n=point.len();
        let mut out=ScanConstraintEvaluation {bound_values:Vec::with_capacity(self.bounds.len()),
            absorbed_energy_j:Vec::with_capacity(self.energy.len()),residuals:Vec::with_capacity(self.rows),
            jacobian:vec![0.0;self.rows*n]};
        for b in &self.bounds {
            poll(cx)?;
            let (value,chain)=if b.coordinate<powers.len() {(powers[b.coordinate],powers[b.coordinate])}
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
                joules=finite(joules+term)?;out.jacobian[offset+i]=finite(term/row.cap.scale_j)?;
            }
            out.absorbed_energy_j.push(joules);
            out.residuals.push(finite(finite(joules-row.cap.maximum_j)?/row.cap.scale_j)?);
        }
        poll(cx)?;
        Ok(out)
    }

    fn sample<'h,'m,'c,F>(&self,cx: &Cx<'_>,study: &mut ScanPowerStudy<F>,point: &[f64])
        -> Result<Option<SqpSample>,ScanSensitivityError>
    where F: FnMut(&Cx<'_>,Vec<MovingHeatSource>) -> Result<ScanningEnthalpyRun<'h,'m,'c>,ScanSensitivityError>
    {
        let constraints=self.at(cx,study,point)?;
        let physical=study.evaluate(cx,point)?;
        Ok(Some(SqpSample {f:physical.physical.value,gradient:physical.gradient,
            ce:Vec::new(),je:Vec::new(),ci:constraints.residuals,ji:constraints.jacobian}))
    }
}
