//! Accuracy-controlled physical intervals of the original mixed viscoelastic motion.
//!
//! One call accepts one pair of native backward-Euler half intervals toward a
//! target. The full interval is an error comparison only. All original geometry,
//! pressure equations, mass and Maxwell histories remain with their native owners.
//! No Richardson extrapolation, pressure projection or history interpolation occurs.
//! This local consistency estimate is not a global error bound, a frequency gate,
//! or a zero-algorithmic-damping integrator.

use fs_material::hyper::mixed::MixedFiniteStrainPotential;

use super::{
    Cx, FiniteViscoMotionControls, MixedUpSettings, MixedViscoMotionStep, MixedViscoTetError,
    MixedViscoTetMotion, checkpoint, gate,
};
use crate::hyper3::mixed::MixedUpError;
use crate::hyper3::response::{green_lagrange_strain, stress_strain_fields};

/// One explicit physical source snapshot at the requested absolute time.
///
/// Supports are absolute displacements. Forces are original nodal dead loads
/// and force per reference volume. The callback owns the time dependence; no
/// midpoint load or support position is inferred from an end-only snapshot.
#[derive(Clone, Debug)]
pub struct MixedViscoPairLoads {
    pub prescribed_m: Vec<(usize, f64)>,
    pub nodal_forces_n: Vec<(usize, f64)>,
    pub body_force_n_m3: [f64; 3],
    pub pressure_gauge: Option<(usize, f64)>,
}

/// Explicit local error tolerances and bounded pair retries; no defaults.
///
/// Vector/tensor errors use Euclidean/Frobenius norms, taking the maximum over
/// original nodes, quadrature points and branches. Relative motion, Green strain
/// and metric scales use changes from the accepted state, so a large initial
/// displacement or prestrain cannot hide a small evolving error. Cauchy stress
/// uses the larger endpoint norm, retaining physical hydrostatic pressure too.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedViscoAdaptiveControls {
    pub absolute_displacement_m: f64,
    pub absolute_velocity_m_s: f64,
    /// Frobenius norm of the total Green-Lagrange strain discrepancy.
    pub absolute_green_strain: f64,
    /// Frobenius norm of each original symmetric C_v^-1 discrepancy.
    pub absolute_inverse_viscous_metric: f64,
    pub absolute_cauchy_stress_pa: f64,
    pub relative: f64,
    /// Retry floor for a PAIR span, not its individual half intervals.
    /// An explicit proposal or the target may initially shorten the pair below this floor.
    pub minimum_pair_span_s: f64,
    pub maximum_pair_span_s: f64,
    /// Includes the first attempt; 1..=32, hence at most 96 native solves.
    /// A successful call accepts exactly two physical history intervals.
    pub maximum_attempts: usize,
}

impl MixedViscoAdaptiveControls {
    fn validate(self) -> Result<(), MixedViscoAdaptiveError> {
        if [
            self.absolute_displacement_m,
            self.absolute_velocity_m_s,
            self.absolute_green_strain,
            self.absolute_inverse_viscous_metric,
            self.absolute_cauchy_stress_pa,
            self.minimum_pair_span_s,
            self.maximum_pair_span_s,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
            || !self.relative.is_finite()
            || !(0.0..1.0).contains(&self.relative)
            || self.minimum_pair_span_s > self.maximum_pair_span_s
            || !(1..=32).contains(&self.maximum_attempts)
        {
            return Err(MixedViscoAdaptiveError::Invalid(
                "positive finite absolute tolerances and pair spans, relative tolerance in [0,1), and 1..=32 attempts required",
            ));
        }
        Ok(())
    }
}

/// Raw full/two-half disagreement divided by the declared local allowances.
/// No order is assumed for prescribed-support derivative recurrence. For smooth
/// first-order native evolution the asymptotic step-doubling divisor is one.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MixedViscoPairEstimate {
    pub displacement_ratio: f64,
    pub velocity_ratio: f64,
    pub green_strain_ratio: f64,
    pub inverse_viscous_metric_ratio: f64,
    pub cauchy_stress_ratio: f64,
}

impl MixedViscoPairEstimate {
    pub fn maximum_ratio(self) -> f64 {
        [
            self.displacement_ratio,
            self.velocity_ratio,
            self.green_strain_ratio,
            self.inverse_viscous_metric_ratio,
            self.cauchy_stress_ratio,
        ]
        .into_iter()
        .fold(0.0, f64::max)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum MixedViscoAdaptiveError {
    Native(MixedViscoTetError),
    Invalid(&'static str),
    Rejected {
        attempts: usize,
        attempted_pair_span_s: f64,
        estimate: MixedViscoPairEstimate,
    },
}

impl From<MixedViscoTetError> for MixedViscoAdaptiveError {
    fn from(value: MixedViscoTetError) -> Self {
        Self::Native(value)
    }
}
impl From<MixedUpError> for MixedViscoAdaptiveError {
    fn from(value: MixedUpError) -> Self {
        Self::Native(value.into())
    }
}
impl core::fmt::Display for MixedViscoAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive mixed viscoelastic pair refused: {self:?}")
    }
}
impl core::error::Error for MixedViscoAdaptiveError {}

/// Both accepted physical intervals, including their actual independent work
/// and loss ledgers. Continue from second.motion; first is the accepted midpoint.
/// The discarded coarse trajectory contributes no accepted history or work.
#[derive(Debug)]
pub struct MixedViscoAdaptivePair<'a> {
    pub first: MixedViscoMotionStep<'a>,
    pub second: MixedViscoMotionStep<'a>,
    pub estimate: MixedViscoPairEstimate,
    pub accepted_pair_span_s: f64,
    pub rejected_attempts: usize,
    pub native_solves: usize,
    pub suggested_pair_span_s: f64,
    /// Sum of ABSOLUTE signed material endpoint-work defects [J], not heat.
    pub endpoint_work_defect_magnitude_j: f64,
    /// Sum of actual backward-Euler kinetic losses [J], not material heat.
    pub numerical_kinetic_loss_j: f64,
    pub balance_error_magnitude_j: f64,
}

impl<'a> MixedViscoTetMotion<'a> {
    /// Advance ONE accepted pair toward target_time_s, possibly ending earlier.
    ///
    /// A pure source callback is evaluated at each actual candidate midpoint
    /// and endpoint. The captured endpoint is shared by the coarse and fine
    /// solves. The caller must bound source bandwidth and end the target at
    /// discontinuities; two samples cannot detect a skipped source cycle.
    ///
    /// All nodes, including prescribed supports, participate in motion checks.
    /// Green strain and natural-frame viscous-metric errors are objective.
    /// Cauchy stress is evaluated from each accepted history at ZERO duration,
    /// with its original pressure and material owner; positive-duration
    /// incremental stress is never used as the physical comparison field.
    ///
    /// Native equilibrium, source, material and work refusals propagate without
    /// retry. Only an excessive local discrepancy shrinks the pair. Both fine
    /// intervals retain all native gates, and their summed absolute endpoint
    /// defects, numerical kinetic loss and balance defects must fit the SAME
    /// supplied work allowances for the pair. Cancellation or any refusal
    /// publishes neither half and leaves this accepted owner unchanged.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_adaptive_pair_toward(
        &self,
        target_time_s: f64,
        proposed_pair_span_s: f64,
        mut endpoint_source: impl FnMut(f64) -> Result<MixedViscoPairLoads, MixedViscoTetError>,
        settings: MixedUpSettings,
        work: FiniteViscoMotionControls,
        control: MixedViscoAdaptiveControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoAdaptivePair<'a>, MixedViscoAdaptiveError> {
        checkpoint(cx)?;
        control.validate()?;
        let start = self.solid.time_s;
        if !target_time_s.is_finite()
            || target_time_s <= start
            || !proposed_pair_span_s.is_finite()
            || proposed_pair_span_s <= 0.0
            || settings.load_steps != 1
            || [
                work.material.equilibrium_work_tolerance_j,
                work.material.maximum_endpoint_work_defect_j,
                work.maximum_numerical_kinetic_loss_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(MixedViscoAdaptiveError::Invalid(
                "a later finite target, positive pair proposal, one native endpoint and finite nonnegative work allowances required",
            ));
        }
        let mut span_cap = proposed_pair_span_s.min(control.maximum_pair_span_s);
        for attempt in 1..=control.maximum_attempts {
            checkpoint(cx)?;
            let end = (start + span_cap).min(target_time_s);
            let span = end - start;
            let middle = start + 0.5 * span;
            if !span.is_finite() || middle <= start || middle >= end {
                return Err(MixedViscoAdaptiveError::Invalid(
                    "adaptive pair has no distinct finite representable midpoint",
                ));
            }
            let middle_load = endpoint_source(middle)?;
            checkpoint(cx)?;
            let end_load = endpoint_source(end)?;
            checkpoint(cx)?;
            let solve = |motion: &Self, dt: f64, load: &MixedViscoPairLoads| {
                motion.advance(
                    dt,
                    &load.prescribed_m,
                    &load.nodal_forces_n,
                    load.body_force_n_m3,
                    load.pressure_gauge,
                    settings,
                    work,
                    cx,
                )
            };
            let coarse = solve(self, span, &end_load)?;
            let first = solve(self, middle - start, &middle_load)?;
            let second = solve(&first.motion, end - middle, &end_load)?;
            if coarse.motion.solid.time_s != end
                || first.motion.solid.time_s != middle
                || second.motion.solid.time_s != end
            {
                return Err(MixedViscoAdaptiveError::Invalid(
                    "native comparison did not reach the same physical endpoint",
                ));
            }
            let estimate = discrepancy(self, &coarse.motion, &second.motion, control, cx)?;
            let error = estimate.maximum_ratio();
            if error <= 1.0 {
                let endpoint_work_defect_magnitude_j =
                    first.material_energy.endpoint_work_defect_j.abs()
                        + second.material_energy.endpoint_work_defect_j.abs();
                let numerical_kinetic_loss_j =
                    first.inertia.numerical_loss_j + second.inertia.numerical_loss_j;
                let balance_error_magnitude_j =
                    first.balance_error_j.abs() + second.balance_error_j.abs();
                gate(
                    "adaptive pair absolute endpoint-work defect",
                    endpoint_work_defect_magnitude_j,
                    work.material.maximum_endpoint_work_defect_j,
                )?;
                gate(
                    "adaptive pair numerical kinetic loss",
                    numerical_kinetic_loss_j,
                    work.maximum_numerical_kinetic_loss_j,
                )?;
                gate(
                    "adaptive pair absolute balance defect",
                    balance_error_magnitude_j,
                    work.material.equilibrium_work_tolerance_j,
                )?;
                let factor = if error == 0.0 {
                    2.0
                } else {
                    (0.8 / fs_math::det::sqrt(error)).clamp(0.5, 2.0)
                };
                let suggested_pair_span_s = (span * factor)
                    .min(control.maximum_pair_span_s)
                    .max(control.minimum_pair_span_s);
                checkpoint(cx)?;
                return Ok(MixedViscoAdaptivePair {
                    first,
                    second,
                    estimate,
                    accepted_pair_span_s: span,
                    rejected_attempts: attempt - 1,
                    native_solves: 3 * attempt,
                    suggested_pair_span_s,
                    endpoint_work_defect_magnitude_j,
                    numerical_kinetic_loss_j,
                    balance_error_magnitude_j,
                });
            }
            let next = (span * (0.8 / fs_math::det::sqrt(error)).clamp(0.2, 0.8))
                .max(control.minimum_pair_span_s);
            if attempt == control.maximum_attempts || next >= span || start + next <= start {
                return Err(MixedViscoAdaptiveError::Rejected {
                    attempts: attempt,
                    attempted_pair_span_s: span,
                    estimate,
                });
            }
            span_cap = next;
        }
        Err(MixedViscoAdaptiveError::Invalid(
            "empty adaptive attempt budget",
        ))
    }
}

fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0, fs_math::det::hypot)
}

fn difference<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    norm(a.iter().zip(b).map(|(a, b)| a - b))
}

fn ratio(
    discrepancy: f64,
    scale: f64,
    absolute: f64,
    relative: f64,
) -> Result<f64, MixedViscoAdaptiveError> {
    let tolerance = relative.mul_add(scale, absolute);
    let value = discrepancy / tolerance;
    if !scale.is_finite() || !tolerance.is_finite() || !value.is_finite() {
        return Err(MixedViscoAdaptiveError::Invalid(
            "unrepresentable physical error comparison",
        ));
    }
    Ok(value)
}

fn change_ratio<const N: usize>(
    old: &[f64; N],
    coarse: &[f64; N],
    fine: &[f64; N],
    absolute: f64,
    relative: f64,
) -> Result<f64, MixedViscoAdaptiveError> {
    ratio(
        difference(coarse, fine),
        difference(old, coarse).max(difference(old, fine)),
        absolute,
        relative,
    )
}

// Material storage packs [xx,yy,zz,xy,xz,yz]. Expansion gives the true
// Frobenius metric, with both occurrences of every off-diagonal component.
fn metric(a: &[f64; 6]) -> [f64; 9] {
    [a[0], a[3], a[4], a[3], a[1], a[5], a[4], a[5], a[2]]
}

fn physical_fields(
    motion: &MixedViscoTetMotion<'_>,
    element: usize,
    local_point: usize,
) -> Result<([f64; 9], [f64; 9]), MixedViscoAdaptiveError> {
    let solid = &motion.solid;
    let points = solid.mesh.points_per_element;
    let geometry = &solid.mesh.elements[element];
    let point = &geometry.points[local_point];
    let accepted = &solid.states[element * points + local_point];
    let pressure = (0..4)
        .map(|a| solid.pressure_pa[geometry.pressure_indices[a]] * point.pressure_shapes[a])
        .sum();
    let material_error = |error| MixedUpError::Material {
        element,
        point: local_point,
        error,
    };
    let frozen = solid
        .material
        .at(element, local_point, points)
        .freeze_remainder(accepted, 0.0)
        .map_err(material_error)?;
    let physical = frozen
        .mixed_state(&accepted.deformation_gradient, pressure, solid.volumetric)
        .map_err(material_error)?;
    let fields = stress_strain_fields(
        &accepted.deformation_gradient,
        &physical.first_piola_pa,
        physical.volume_ratio,
    );
    if !physical.volume_ratio.is_finite()
        || physical.volume_ratio <= 0.0
        || fields.0.iter().chain(&fields.1).any(|x| !x.is_finite())
    {
        return Err(MixedViscoAdaptiveError::Invalid(
            "unrepresentable accepted strain or physical stress",
        ));
    }
    Ok(fields)
}

fn discrepancy(
    old: &MixedViscoTetMotion<'_>,
    coarse: &MixedViscoTetMotion<'_>,
    fine: &MixedViscoTetMotion<'_>,
    control: MixedViscoAdaptiveControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoPairEstimate, MixedViscoAdaptiveError> {
    let mut result = MixedViscoPairEstimate::default();
    for i in 0..old.solid.nodes_m.len() {
        if i % 256 == 0 {
            checkpoint(cx)?;
        }
        result.displacement_ratio = result.displacement_ratio.max(change_ratio(
            &old.solid.displacement_m[i],
            &coarse.solid.displacement_m[i],
            &fine.solid.displacement_m[i],
            control.absolute_displacement_m,
            control.relative,
        )?);
        result.velocity_ratio = result.velocity_ratio.max(change_ratio(
            &old.velocity_mps[i],
            &coarse.velocity_mps[i],
            &fine.velocity_mps[i],
            control.absolute_velocity_m_s,
            control.relative,
        )?);
    }
    let points = old.solid.mesh.points_per_element;
    for (index, accepted) in old.solid.states.iter().enumerate() {
        checkpoint(cx)?;
        let (e, q) = (index / points, index % points);
        let full = &coarse.solid.states[index];
        let half = &fine.solid.states[index];
        let (cauchy_full, green_full) = physical_fields(coarse, e, q)?;
        let (cauchy_half, green_half) = physical_fields(fine, e, q)?;
        result.green_strain_ratio = result.green_strain_ratio.max(change_ratio(
            &green_lagrange_strain(&accepted.deformation_gradient),
            &green_full,
            &green_half,
            control.absolute_green_strain,
            control.relative,
        )?);
        result.cauchy_stress_ratio = result.cauchy_stress_ratio.max(ratio(
            difference(&cauchy_full, &cauchy_half),
            norm(cauchy_full.into_iter()).max(norm(cauchy_half.into_iter())),
            control.absolute_cauchy_stress_pa,
            control.relative,
        )?);
        for ((base, full), half) in accepted
            .inverse_viscous_cauchy_green
            .iter()
            .zip(&full.inverse_viscous_cauchy_green)
            .zip(&half.inverse_viscous_cauchy_green)
        {
            result.inverse_viscous_metric_ratio =
                result.inverse_viscous_metric_ratio.max(change_ratio(
                    &metric(base),
                    &metric(full),
                    &metric(half),
                    control.absolute_inverse_viscous_metric,
                    control.relative,
                )?);
        }
    }
    checkpoint(cx)?;
    Ok(result)
}
