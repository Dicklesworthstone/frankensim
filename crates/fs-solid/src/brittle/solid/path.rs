//! Adaptive physical load intervals for the original retained brittle solid.
//!
//! Prescribed coordinates interpolate from the accepted displacement to the
//! requested endpoint; dead forces interpolate between their explicit initial
//! and final values. The original native force assembly checks that the initial
//! free equations are in equilibrium before any continuation trial.
//!
//! Rejected trials never change retained history. Subdivision controls the
//! largest element damage increment and the sum of absolute endpoint-work
//! defects. Signed cancellation cannot hide numerical work error. These gates
//! do not estimate the full displacement error or establish crack-path
//! convergence. Loading reversals require separate intervals.
//!
//! This is quasistatic affine load continuation at fixed material conditions,
//! not pseudo-arclength continuation through a force-controlled limit point,
//! dynamics, a physical clock, moving crack planes or fragment contact.

use fs_exec::Cx;

use super::{BrittleTetControls, BrittleTetEnergy, BrittleTetError, BrittleTetSolid, checkpoint};
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
};

/// Explicit physical-subdivision, memory, work and initial-equilibrium bounds.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetPathControls {
    pub initial_fraction: f64,
    /// Smallest attempted interval, except a final shorter remainder.
    pub minimum_fraction: f64,
    pub maximum_fraction: f64,
    /// Maximum number of compact accepted-step records.
    pub maximum_accepted_steps: usize,
    /// Maximum native equilibrium solves, including rejected trials.
    pub maximum_trials: usize,
    /// Largest admitted damage increase of any element in an accepted substep.
    pub maximum_damage_increment: f64,
    /// Absolute infinity norm of the initial free-force residual, N.
    pub initial_equilibrium_tolerance_n: f64,
    /// Material/equilibrium tolerances apply to each endpoint and the complete
    /// interval. The endpoint-defect allowance is a WHOLE-INTERVAL budget:
    /// each trial receives this value times its actual interval fraction.
    pub work: BrittleTetControls,
}

/// A compact observation of one accepted endpoint, without copying its fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetPathRecord {
    /// Cumulative fraction reached within the requested affine interval.
    pub fraction: f64,
    pub maximum_damage: f64,
    pub maximum_damage_increment: f64,
    pub energy: BrittleTetEnergy,
}

/// A complete accepted interval. The native solution retained by `solid`
/// describes its final substep; `records` covers every accepted substep.
#[derive(Debug)]
pub struct BrittleTetPathStep<'a> {
    pub solid: BrittleTetSolid<'a>,
    pub records: Vec<BrittleTetPathRecord>,
    pub trials: usize,
    pub rejected_trials: usize,
    pub energy: BrittleTetEnergy,
    /// Sum of absolute local endpoint-work defects, never a physical loss.
    pub absolute_endpoint_work_defect_j: f64,
}

/// A refinable numerical endpoint or a resolved damage increment that is too large.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetRefinementCause {
    Endpoint(BrittleTetError),
    DamageIncrement { observed: f64, maximum: f64 },
}

/// Refusal publishes no partially advanced solid.
#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetPathError {
    Invalid(&'static str),
    Allocation,
    Endpoint(BrittleTetError),
    InitialEquilibrium {
        maximum_residual_n: f64,
        tolerance_n: f64,
    },
    WorkLimit {
        accepted_steps: usize,
        trials: usize,
    },
    MinimumInterval {
        accepted_fraction: f64,
        attempted_fraction: f64,
        cause: BrittleTetRefinementCause,
    },
}

impl core::fmt::Display for BrittleTetPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "brittle physical load path refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetPathError {}
impl From<BrittleTetError> for BrittleTetPathError {
    fn from(error: BrittleTetError) -> Self {
        Self::Endpoint(error)
    }
}

impl<'a> BrittleTetSolid<'a> {
    /// Resolve one affine physical load interval using the original `advance`.
    ///
    /// `prescribed_end_m` defines a fixed support set; its initial coordinates
    /// are the retained displacement at those DOFs. `force_intervals_n` supplies
    /// each dead force as `(dof, [initial_n, final_n])`; repeated force entries
    /// add, as in the native solver. The initial free equations must satisfy the
    /// explicit force tolerance, so releasing a support requires a compatible
    /// retained equilibrium. New source conditions require the temperature API.
    ///
    /// Only typed Newton/linear-solve failures, finite endpoint-work defects and
    /// excessive damage increments trigger subdivision. Material domains,
    /// identity, geometry, other energy gates and cancellation are terminal.
    /// Every failed complete request leaves `self` unchanged.
    pub fn advance_path(
        &self,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleTetPathControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetPathStep<'a>, BrittleTetPathError> {
        checkpoint(cx)?;
        validate_controls(controls)?;
        let old_u = self.solution.displacement_m.as_flattened();
        if prescribed_end_m.len() > old_u.len() || force_intervals_n.len() > old_u.len() {
            return Err(BrittleTetPathError::Invalid(
                "load entries exceed the original DOF count",
            ));
        }
        let mut prescribed = buffer(prescribed_end_m.len())?;
        for &(dof, end) in prescribed_end_m {
            checkpoint(cx)?;
            let Some(&start) = old_u.get(dof) else {
                return Err(BrittleTetPathError::Invalid(
                    "prescribed DOF is outside the mesh",
                ));
            };
            if !end.is_finite() {
                return Err(BrittleTetPathError::Invalid(
                    "finite prescribed endpoint required",
                ));
            }
            prescribed.push((dof, start));
        }
        let mut forces = buffer(force_intervals_n.len())?;
        for &(dof, values) in force_intervals_n {
            checkpoint(cx)?;
            if dof >= old_u.len() || values.iter().any(|v| !v.is_finite()) {
                return Err(BrittleTetPathError::Invalid(
                    "finite forces on existing DOFs required",
                ));
            }
            forces.push((dof, values[0]));
        }
        // Reuse native shape gradients, constitutive stress and DOF admission;
        // no second FE residual or material law is reconstructed here.
        {
            let mut laws = buffer(self.materials.len())?;
            for source in self.materials {
                checkpoint(cx)?;
                laws.push(source.small_strain_law());
            }
            let problem = SmallStrainTetProblem {
                nodes_m: self.geometry.nodes_m(),
                tetrahedra: self.geometry.tetrahedra(),
                laws: SmallStrainLawField::PerElement(&laws),
                prescribed_m: &prescribed,
                nodal_forces_n: &forces,
                load_path: &[1.0],
                budget: self.budget,
                settings,
            };
            let mut residual = problem
                .internal_force(old_u, &self.solution.element_states, cx)
                .map_err(BrittleTetError::from)?;
            for (index, &(dof, force)) in forces.iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                residual[dof] -= force;
                if !residual[dof].is_finite() {
                    return Err(BrittleTetPathError::Invalid(
                        "initial force sum is not finite",
                    ));
                }
            }
            for (index, &(dof, _)) in prescribed.iter().enumerate() {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                residual[dof] = 0.0;
            }
            let mut maximum = 0.0_f64;
            for (dof, value) in residual.iter().enumerate() {
                if dof % 256 == 0 {
                    checkpoint(cx)?;
                }
                if !value.is_finite() {
                    return Err(BrittleTetPathError::Invalid(
                        "initial force residual is not finite",
                    ));
                }
                maximum = maximum.max(value.abs());
            }
            if maximum > controls.initial_equilibrium_tolerance_n {
                return Err(BrittleTetPathError::InitialEquilibrium {
                    maximum_residual_n: maximum,
                    tolerance_n: controls.initial_equilibrium_tolerance_n,
                });
            }
        }

        let mut current: Option<BrittleTetSolid<'a>> = None;
        let mut records = buffer(controls.maximum_accepted_steps.min(64))?;
        let mut sums = [Sum::default(); 6];
        let mut progress = 0.0_f64;
        let mut fraction = controls.initial_fraction;
        let mut trials = 0;
        let mut rejected = 0;
        while progress < 1.0 {
            checkpoint(cx)?;
            if trials >= controls.maximum_trials || records.len() >= controls.maximum_accepted_steps
            {
                return Err(BrittleTetPathError::WorkLimit {
                    accepted_steps: records.len(),
                    trials,
                });
            }
            let end = (progress + fraction).min(1.0);
            let span = end - progress;
            if end <= progress || !span.is_finite() {
                return Err(BrittleTetPathError::Invalid(
                    "physical interval cannot be represented",
                ));
            }
            for (index, (entry, &(dof, target))) in
                prescribed.iter_mut().zip(prescribed_end_m).enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = interpolate(old_u[dof], target, end);
            }
            for (index, (entry, &(_, [start, target]))) in
                forces.iter_mut().zip(force_intervals_n).enumerate()
            {
                if index % 256 == 0 {
                    checkpoint(cx)?;
                }
                entry.1 = interpolate(start, target, end);
            }
            let allowance = controls.work.maximum_endpoint_work_defect_j * span;
            let endpoint_controls = BrittleTetControls {
                maximum_endpoint_work_defect_j: allowance,
                ..controls.work
            };
            let accepted = current.as_ref().unwrap_or(self);
            trials += 1;
            let candidate =
                match accepted.advance(&prescribed, &forces, settings, endpoint_controls, cx) {
                    Ok(candidate) => candidate,
                    Err(error) if retryable(&error) => {
                        fraction = refine(
                            progress,
                            span,
                            controls.minimum_fraction,
                            BrittleTetRefinementCause::Endpoint(error),
                        )?;
                        rejected += 1;
                        continue;
                    }
                    Err(error) => return Err(error.into()),
                };
            let mut maximum_damage = 0.0_f64;
            let mut maximum_increment = 0.0_f64;
            for (element, (old, new)) in accepted
                .solution
                .element_states
                .iter()
                .zip(&candidate.solid.solution.element_states)
                .enumerate()
            {
                if element % 256 == 0 {
                    checkpoint(cx)?;
                }
                maximum_damage = maximum_damage.max(new.damage());
                maximum_increment = maximum_increment.max(new.damage() - old.damage());
            }
            if maximum_increment > controls.maximum_damage_increment {
                fraction = refine(
                    progress,
                    span,
                    controls.minimum_fraction,
                    BrittleTetRefinementCause::DamageIncrement {
                        observed: maximum_increment,
                        maximum: controls.maximum_damage_increment,
                    },
                )?;
                rejected += 1;
                continue;
            }
            let energy = candidate.energy;
            let defect = energy.endpoint_work_defect_j.abs();
            for (sum, value) in sums.iter_mut().zip([
                energy.material_work_j,
                energy.fracture_dissipation_j,
                energy.internal_endpoint_work_j,
                energy.dead_load_endpoint_work_j,
                energy.support_endpoint_work_j,
                defect,
            ]) {
                sum.add(value);
                if !sum.value().is_finite() {
                    return Err(BrittleTetPathError::Invalid(
                        "interval work cannot be represented",
                    ));
                }
            }
            if records.len() == records.capacity() {
                records
                    .try_reserve(1)
                    .map_err(|_| BrittleTetPathError::Allocation)?;
            }
            records.push(BrittleTetPathRecord {
                fraction: end,
                maximum_damage,
                maximum_damage_increment: maximum_increment,
                energy,
            });
            current = Some(candidate.solid);
            progress = end;
            let work_growth = if defect > 0.0 {
                (0.8 * allowance / defect).clamp(0.5, 2.0)
            } else {
                2.0
            };
            let damage_growth = if maximum_increment > 0.0 {
                (0.8 * controls.maximum_damage_increment / maximum_increment).clamp(0.5, 2.0)
            } else {
                2.0
            };
            fraction = (span * work_growth.min(damage_growth))
                .clamp(controls.minimum_fraction, controls.maximum_fraction);
        }
        let solid = current.ok_or(BrittleTetPathError::Invalid(
            "no physical endpoint accepted",
        ))?;
        let [
            material_work,
            fracture,
            internal,
            dead,
            support,
            absolute_defect,
        ] = sums.map(Sum::value);
        let energy = BrittleTetEnergy {
            stored_before_j: self.stored_energy_j,
            stored_after_j: solid.stored_energy_j,
            material_work_j: material_work,
            fracture_dissipation_j: fracture,
            material_balance_error_j: (material_work
                - (solid.stored_energy_j - self.stored_energy_j))
                - fracture,
            internal_endpoint_work_j: internal,
            dead_load_endpoint_work_j: dead,
            support_endpoint_work_j: support,
            equilibrium_work_residual_j: (dead + support) - internal,
            endpoint_work_defect_j: internal - material_work,
        };
        super::energy_gate(
            "path material work/storage/fracture",
            energy.material_balance_error_j,
            controls.work.material_energy_tolerance_j,
        )?;
        super::energy_gate(
            "path cumulative fracture loss",
            (solid.fracture_dissipation_j - self.fracture_dissipation_j) - fracture,
            controls.work.material_energy_tolerance_j,
        )?;
        super::energy_gate(
            "path support/dead-load equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.work.equilibrium_work_tolerance_j,
        )?;
        super::energy_gate(
            "path absolute endpoint integration defect",
            absolute_defect,
            controls.work.maximum_endpoint_work_defect_j,
        )?;
        checkpoint(cx)?;
        Ok(BrittleTetPathStep {
            solid,
            records,
            trials,
            rejected_trials: rejected,
            energy,
            absolute_endpoint_work_defect_j: absolute_defect,
        })
    }
}

fn validate_controls(c: BrittleTetPathControls) -> Result<(), BrittleTetPathError> {
    if ![
        c.minimum_fraction,
        c.initial_fraction,
        c.maximum_fraction,
        c.maximum_damage_increment,
    ]
    .iter()
    .all(|x| x.is_finite() && *x > 0.0)
        || c.minimum_fraction > c.initial_fraction
        || c.initial_fraction > c.maximum_fraction
        || c.maximum_fraction > 1.0
        || c.maximum_damage_increment > 1.0
        || c.maximum_accepted_steps == 0
        || c.maximum_trials == 0
        || ![
            c.initial_equilibrium_tolerance_n,
            c.work.material_energy_tolerance_j,
            c.work.equilibrium_work_tolerance_j,
            c.work.maximum_endpoint_work_defect_j,
        ]
        .iter()
        .all(|x| x.is_finite() && *x >= 0.0)
    {
        return Err(BrittleTetPathError::Invalid(
            "ordered positive fractions, positive work counts and finite nonnegative tolerances required",
        ));
    }
    Ok(())
}
fn buffer<T>(capacity: usize) -> Result<Vec<T>, BrittleTetPathError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| BrittleTetPathError::Allocation)?;
    Ok(result)
}
pub(super) fn interpolate(start: f64, end: f64, fraction: f64) -> f64 {
    // Preserve the caller's exact endpoint and avoid forming end-start,
    // which can overflow for finite endpoints of opposite sign.
    if fraction == 1.0 {
        end
    } else {
        fraction.mul_add(end, (1.0 - fraction) * start)
    }
}
fn retryable(error: &BrittleTetError) -> bool {
    matches!(
        error,
        BrittleTetError::Solve(
            SmallStrainTetError::NewtonStalled { .. }
                | SmallStrainTetError::LinearSolveFailed { .. }
        )
    ) || matches!(error, BrittleTetError::EnergyGate {
            account: "signed endpoint integration defect", defect_j, ..
        } if defect_j.is_finite())
}
fn refine(
    progress: f64,
    span: f64,
    minimum: f64,
    cause: BrittleTetRefinementCause,
) -> Result<f64, BrittleTetPathError> {
    let smaller = 0.5 * span;
    if smaller < minimum || progress + smaller <= progress {
        Err(BrittleTetPathError::MinimumInterval {
            accepted_fraction: progress,
            attempted_fraction: span,
            cause,
        })
    } else {
        Ok(smaller)
    }
}
#[derive(Clone, Copy, Default)]
pub(super) struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    pub(super) fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    pub(super) fn value(self) -> f64 {
        self.sum + self.correction
    }
}
