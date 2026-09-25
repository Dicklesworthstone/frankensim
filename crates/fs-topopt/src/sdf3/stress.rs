//! Exact discrete qp-relaxed stress sensitivities on retained 3-D cut quadrature.
//!
//! `s_lq = rho_bar[cell]^q * vm(C_ref eps(u_l))`, with `1 <= q < penal`.
//! The aggregate is the p-th root of the volume-normalized integral of `s^p`,
//! averaged with normalized positive load weights. Loads are solved separately;
//! opposite loads cannot cancel stresses. This is a smooth quadrature aggregate,
//! not a continuum maximum or an upper bound on any sampled maximum.
//!
//! The reference material in the stress measure is deliberate: multiplying the
//! already density-softened stress by rho^q can make a force-loaded void appear
//! stronger. The actual SIMP stress is reported separately. The derivative adds
//! the direct relaxation term to `-z^T (dK/dscale) u * dscale/drho_bar`, then
//! follows the same projection and graph-filter transpose as compliance.
//! With a nonzero ersatz modulus, the algebraic qp measure eventually folds
//! back as density approaches zero. This evaluator exposes that declared model;
//! the minimum-volume driver separately admits only densities above its turnover.

use super::*;
use crate::SolveWork;
use fs_cutfem::elastic3::{ElasticityError3, stress::BulkStressPoint3};
use std::ops::ControlFlow;

#[derive(Debug, Clone, Copy)]
pub struct StressOptions3 {
    /// qp-relaxation exponent. Values below one are not admitted because their
    /// derivative is singular at an exactly zero projected density.
    pub relaxation_power: f64,
    /// Fixed p in [2, 64]. No iteration-dependent maximum rescaling is applied.
    pub aggregation_power: f64,
    pub max_cases: usize,
    /// Total retained bulk points across ALL independent load cases, including
    /// zero-weight cases. Checked before allocating each case's point vector.
    pub max_points: usize,
}
impl Default for StressOptions3 {
    fn default() -> Self {
        Self {
            relaxation_power: 1.0,
            aggregation_power: 8.0,
            max_cases: 32,
            max_points: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum StressError3 {
    Invalid(&'static str),
    Evaluation(EvaluationStop),
    Physics(ElasticityError3),
    /// No partial point family or aggregate is returned on an allowance stop.
    PointBudget,
}
impl From<EvaluationStop> for StressError3 {
    fn from(e: EvaluationStop) -> Self {
        Self::Evaluation(e)
    }
}
impl From<ElasticityError3> for StressError3 {
    fn from(e: ElasticityError3) -> Self {
        match e {
            ElasticityError3::Cancelled => Self::Evaluation(EvaluationStop::Cancelled),
            ElasticityError3::Invalid("bulk stress point allowance exhausted") => Self::PointBudget,
            e => Self::Physics(e),
        }
    }
}
impl std::fmt::Display for StressError3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "3-D stress evaluation refused: {self:?}")
    }
}
impl std::error::Error for StressError3 {}

/// Complete independently loaded equilibrium and stress derivative at one raw
/// design. Sampled maxima are diagnostics, not certified continuum maxima.
#[derive(Debug, Clone)]
pub struct StressEvaluation3 {
    pub rho: Vec<f64>,
    pub projected_rho: Vec<f64>,
    pub scales: Vec<f64>,
    pub aggregate: f64,
    pub gradient: Vec<f64>,
    /// Maximum over all quadrature points of all positive-weight load cases.
    pub sampled_relaxed_max: f64,
    /// Maximum of vm(scale * C_ref eps(u)), separately from the relaxed measure.
    pub sampled_physical_max: f64,
    pub case_relaxed_max: Vec<f64>,
    pub case_physical_max: Vec<f64>,
    /// Per-case active-cell maxima, in the operator's exact density order.
    pub cell_relaxed_max: Vec<Vec<f64>>,
    pub normalized_load_weights: Vec<f64>,
    pub volume_fraction: f64,
    pub volume_gradient: Vec<f64>,
    pub displacements: Vec<Vec<f64>>,
    /// One adjoint of the whole aggregate per independent load case.
    pub adjoints: Vec<Vec<f64>>,
    pub case_compliances: Vec<f64>,
    pub point_count: usize,
    pub work: SolveWork,
}

pub(super) fn admit<O: Sdf3Elasticity>(
    study: &CutDensityStudy3<O>,
    rho: &[f64],
    loads: &[LoadCase<'_>],
    options: StressOptions3,
) -> Result<(), StressError3> {
    if !options.relaxation_power.is_finite()
        || options.relaxation_power < 1.0
        || options.relaxation_power >= study.params.penal
        || !options.aggregation_power.is_finite()
        || !(2.0..=64.0).contains(&options.aggregation_power)
        || options.max_cases == 0
        || options.max_points == 0
        || loads.is_empty()
        || loads.len() > options.max_cases
        || !loads.iter().any(|l| l.weight > 0.0)
    {
        return Err(StressError3::Invalid(
            "invalid qp stress policy or load count",
        ));
    }
    if rho.len() != study.cells()
        || rho
            .iter()
            .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
        || loads.iter().any(|l| {
            !l.weight.is_finite()
                || l.weight < 0.0
                || l.force.len() != study.operator.n()
                || l.force.iter().any(|f| !f.is_finite())
        })
    {
        return Err(StressError3::Invalid(
            "invalid stress design or independent loads",
        ));
    }
    Ok(())
}
fn poll(control: &mut SolveControl<'_>) -> ControlFlow<()> {
    if control.checkpoint("sdf3-stress-physics").is_ok() {
        ControlFlow::Continue(())
    } else {
        ControlFlow::Break(())
    }
}
fn finite(value: f64) -> Result<f64, StressError3> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(failure("sdf3-stress-arithmetic").into())
    }
}
fn power(r: f64, p: f64) -> f64 {
    if p == 0.0 {
        1.0
    } else if r == 0.0 {
        0.0
    } else {
        fs_math::det::pow(r, p)
    }
}
// Scale first so a finite large stress does not overflow merely on squaring.
// Hydrostatic and zero stress have the zero subgradient of the Euclidean norm;
// their aggregate contribution is differentiable for every admitted p >= 2.
fn von_mises(s: &[f64; 6]) -> (f64, [f64; 6]) {
    let scale = s.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return (0.0, [0.0; 6]);
    }
    let t = s.map(|v| v / scale);
    let v = fs_math::det::sqrt(
        0.5 * ((t[0] - t[1]).powi(2) + (t[1] - t[2]).powi(2) + (t[2] - t[0]).powi(2))
            + 3.0 * (t[3] * t[3] + t[4] * t[4] + t[5] * t[5]),
    );
    if v == 0.0 {
        return (0.0, [0.0; 6]);
    }
    (
        scale * v,
        [
            (2.0 * t[0] - t[1] - t[2]) / (2.0 * v),
            (2.0 * t[1] - t[0] - t[2]) / (2.0 * v),
            (2.0 * t[2] - t[0] - t[1]) / (2.0 * v),
            3.0 * t[3] / v,
            3.0 * t[4] / v,
            3.0 * t[5] / v,
        ],
    )
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Evaluate real equilibria, the reference-material stress quadrature and
    /// exact non-self-adjoint density sensitivities. All primal, adjoint and
    /// filter solves share one cumulative control. The preconditioner is built
    /// once at the evaluated scales and shared by all solves at this design.
    ///
    /// Every result, including success, restores incoming operator scales. Only
    /// the owning constrained optimizer installs accepted scales. Loads must be
    /// fixed external forces on homogeneous supports; nonzero prescribed-motion
    /// lifting is not a density-independent load and is outside this API.
    /// A zero aggregate under nonzero applied load is refused because no
    /// differentiable stress branch has been established there.
    pub fn evaluate_stress(
        &mut self,
        rho: &[f64],
        loads: &[LoadCase<'_>],
        options: StressOptions3,
        control: &mut SolveControl<'_>,
    ) -> Result<StressEvaluation3, StressError3> {
        admit(self, rho, loads, options)?;
        control.checkpoint("sdf3-stress-evaluation")?;
        let design = self.design(rho, control)?;
        let previous = self.operator.scales().to_vec();
        self.operator.set_scales(&design.scales)?;
        let outcome = (|| {
            let prepared = self.operator.prepare_elasticity(control)?;
            let weight_scale = loads.iter().map(|l| l.weight).fold(0.0_f64, f64::max);
            let total_weight: f64 = loads.iter().map(|l| l.weight / weight_scale).sum();
            let weights: Vec<f64> = loads
                .iter()
                .map(|l| (l.weight / weight_scale) / total_weight)
                .collect();
            let volume = finite(self.operator.volumes().iter().sum())?;
            let relaxation: Vec<f64> = design
                .projected
                .iter()
                .map(|r| power(*r, options.relaxation_power))
                .collect();
            let mut points: Vec<Vec<BulkStressPoint3>> = Vec::with_capacity(loads.len());
            let mut displacements = Vec::with_capacity(loads.len());
            let mut compliances = Vec::with_capacity(loads.len());
            let mut case_relaxed_max = Vec::with_capacity(loads.len());
            let mut case_physical_max = Vec::with_capacity(loads.len());
            let mut cell_relaxed_max = Vec::with_capacity(loads.len());
            let mut point_count = 0usize;
            let mut maximum = 0.0_f64;
            let mut physical_maximum = 0.0_f64;
            let mut nonzero_weighted_force = false;
            for (case, load) in loads.iter().enumerate() {
                control.checkpoint("sdf3-stress-primal")?;
                let rhs: Vec<f64> = load
                    .force
                    .iter()
                    .enumerate()
                    .map(|(i, &f)| if self.operator.fixed()[i / 3] { 0.0 } else { f })
                    .collect();
                nonzero_weighted_force |= weights[case] > 0.0 && rhs.iter().any(|f| *f != 0.0);
                let u = checked_solve_preconditioned(
                    &self.operator,
                    &prepared,
                    &rhs,
                    1e-11,
                    "sdf3-stress-primal",
                    control,
                )?;
                let compliance = finite(rhs.iter().zip(&u).map(|(f, u)| f * u).sum())?;
                let data =
                    self.operator
                        .bulk_stress(&u, options.max_points - point_count, || poll(control))?;
                point_count = point_count
                    .checked_add(data.len())
                    .ok_or(StressError3::PointBudget)?;
                if point_count > options.max_points {
                    return Err(StressError3::PointBudget);
                }
                let mut cm = 0.0_f64;
                let mut pm = 0.0_f64;
                let mut cells = vec![0.0_f64; self.cells()];
                for (index, point) in data.iter().enumerate() {
                    if index % 256 == 0 {
                        control.checkpoint("sdf3-stress-measures")?;
                    }
                    let stress =
                        finite(relaxation[point.cell] * von_mises(&point.reference_stress).0)?;
                    cm = cm.max(stress);
                    pm = pm.max(finite(von_mises(&point.stress).0)?);
                    cells[point.cell] = cells[point.cell].max(stress);
                }
                if weights[case] > 0.0 {
                    maximum = maximum.max(cm);
                    physical_maximum = physical_maximum.max(pm);
                }
                points.push(data);
                displacements.push(u);
                compliances.push(compliance);
                case_relaxed_max.push(cm);
                case_physical_max.push(pm);
                cell_relaxed_max.push(cells);
            }
            // Normalize by the observed maximum only to avoid overflow. This
            // is an algebraically equivalent p-norm, not a stop-gradient scale.
            let mut sum = 0.0;
            if maximum > 0.0 {
                for (case, data) in points.iter().enumerate() {
                    if weights[case] == 0.0 {
                        continue;
                    }
                    for (index, point) in data.iter().enumerate() {
                        if index % 256 == 0 {
                            control.checkpoint("sdf3-stress-aggregate")?;
                        }
                        let ratio =
                            relaxation[point.cell] * von_mises(&point.reference_stress).0 / maximum;
                        sum += weights[case]
                            * (point.weight / volume)
                            * power(ratio, options.aggregation_power);
                    }
                }
            }
            let root = power(finite(sum)?, 1.0 / options.aggregation_power);
            let aggregate = finite(maximum * root)?;
            if aggregate == 0.0 && nonzero_weighted_force {
                return Err(StressError3::Invalid(
                    "zero stress aggregate has no admitted differentiable branch",
                ));
            }
            let mut local = vec![0.0; self.cells()];
            let mut adjoints = Vec::with_capacity(loads.len());
            for (case, data) in points.iter().enumerate() {
                control.checkpoint("sdf3-stress-adjoint-rhs")?;
                let mut derivatives = Vec::with_capacity(data.len());
                for (index, point) in data.iter().enumerate() {
                    if index % 256 == 0 {
                        control.checkpoint("sdf3-stress-derivative")?;
                    }
                    let (vm, dvm) = von_mises(&point.reference_stress);
                    let derivative = if aggregate == 0.0 || weights[case] == 0.0 {
                        0.0
                    } else {
                        let ratio = relaxation[point.cell] * vm / maximum;
                        finite(
                            weights[case]
                                * (point.weight / volume)
                                * power(ratio / root, options.aggregation_power - 1.0),
                        )?
                    };
                    let r = design.projected[point.cell];
                    local[point.cell] += derivative
                        * vm
                        * options.relaxation_power
                        * power(r, options.relaxation_power - 1.0);
                    let scale = derivative * relaxation[point.cell];
                    derivatives.push(dvm.map(|v| scale * v));
                }
                let rhs = self
                    .operator
                    .reference_bulk_stress_pullback(&derivatives, || poll(control))?;
                let z = checked_solve_preconditioned(
                    &self.operator,
                    &prepared,
                    &rhs,
                    1e-11,
                    "sdf3-stress-adjoint",
                    control,
                )?;
                let cross = self
                    .operator
                    .scale_bilinear_forms(&z, &displacements[case], || poll(control))?;
                for ((value, &r), contraction) in local.iter_mut().zip(&design.projected).zip(cross)
                {
                    let dsimp = (1.0 - self.params.e_min)
                        * self.params.penal
                        * power(r, self.params.penal - 1.0);
                    *value -= dsimp * contraction;
                }
                adjoints.push(z);
            }
            for (value, slope) in local.iter_mut().zip(&design.slope) {
                *value = finite(*value * slope)?;
            }
            let gradient = self.pullback(&local, control)?;
            let local_volume: Vec<f64> = self
                .mass
                .iter()
                .zip(&design.slope)
                .map(|(m, s)| m * s)
                .collect();
            let volume_gradient = self.pullback(&local_volume, control)?;
            control.checkpoint("sdf3-stress-publish")?;
            Ok(StressEvaluation3 {
                rho: rho.to_vec(),
                projected_rho: design.projected,
                scales: design.scales,
                aggregate,
                gradient,
                sampled_relaxed_max: maximum,
                sampled_physical_max: physical_maximum,
                case_relaxed_max,
                case_physical_max,
                cell_relaxed_max,
                normalized_load_weights: weights,
                volume_fraction: design.volume,
                volume_gradient,
                displacements,
                adjoints,
                case_compliances: compliances,
                point_count,
                work: control.work(),
            })
        })();
        self.operator
            .set_scales(&previous)
            .expect("incoming stress scales remain valid");
        outcome
    }
}
