//! Instantaneous impulses in the complete solid's original free coordinates.
//!
//! Both propagation owners use their prepared mass response, then audit the
//! actual represented velocity against the original consistent mass. Positions,
//! clocks and material/source preparation stay fixed. The reduced pencil has
//! no discarded support rows: this operation does not infer support reactions
//! or claim a full-body angular-momentum inventory from an incomplete map.

use super::{ExactLinearSolid, ExactLinearSolidEnergy, quadratic};
use fs_exec::Cx;
use fs_math::dd::Dd;

/// Failure leaves the accepted motion and all cumulative accounts unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum LinearSolidImpulseError {
    Invalid(&'static str),
    Momentum {
        defect_kg_m_s: f64,
        allowance_kg_m_s: f64,
    },
    Energy {
        defect_j: f64,
        allowance_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for LinearSolidImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "complete solid impulse refused: {self:?}")
    }
}
impl core::error::Error for LinearSolidImpulseError {}

/// Every vector follows the owner's original `free_dofs()` order.
#[derive(Clone, Debug)]
pub struct LinearSolidImpulseReport {
    pub time_s: f64,
    /// M times the actual representable velocity change, not subtraction of
    /// potentially much larger before/after momentum totals [kg m/s].
    pub momentum_change_kg_m_s: Vec<f64>,
    /// Applied impulse minus the actual momentum change [kg m/s].
    pub impulse_residual_kg_m_s: Vec<f64>,
    pub maximum_impulse_residual_kg_m_s: f64,
    /// Average pre/post velocity dotted with the actual momentum change [J].
    pub kinetic_change_j: f64,
    /// Supplied work is J dot average(v). No elapsed time or physical loss is
    /// introduced; earlier viscosity losses and cumulative work are retained.
    pub energy: ExactLinearSolidEnergy,
}

type Result<T> = core::result::Result<T, LinearSolidImpulseError>;

pub(super) fn poll(cx: &Cx<'_>) -> Result<()> {
    cx.checkpoint()
        .map_err(|_| LinearSolidImpulseError::Cancelled)
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LinearSolidImpulseError::Invalid(
            "unrepresentable impulse, velocity or work",
        ))
    }
}

pub(super) fn admit(impulse: &[f64], n: usize, allowance: f64, cx: &Cx<'_>) -> Result<()> {
    poll(cx)?;
    if impulse.len() != n
        || impulse.iter().any(|j| !j.is_finite())
        || !allowance.is_finite()
        || allowance < 0.0
    {
        return Err(LinearSolidImpulseError::Invalid(
            "complete finite free-coordinate impulses and a finite nonnegative impulse allowance required",
        ));
    }
    Ok(())
}

pub(super) fn action(matrix: &[f64], x: &[f64], cx: &Cx<'_>) -> Result<Vec<f64>> {
    let n = x.len();
    let mut result = vec![0.0; n];
    for i in 0..n {
        poll(cx)?;
        result[i] = finite((0..n).fold(0.0, |s, j| matrix[i * n + j].mul_add(x[j], s)))?;
    }
    Ok(result)
}

fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    // Opposing impulses may do nearly cancelling signed work. This O(n)
    // account keeps the product roundoff, without changing the O(n^2) mass map.
    let value = a.iter().zip(b).fold(Dd::ZERO, |s, (&a, &b)| {
        s + Dd::from_f64(a) * Dd::from_f64(b)
    });
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())
}

pub(super) struct History {
    pub time_s: f64,
    pub stored_j: f64,
    pub kinetic_before_j: f64,
    pub initial_energy_j: f64,
    pub cumulative_work_j: f64,
    pub work_correction_j: f64,
    pub cumulative_dissipated_j: f64,
    pub maximum_work_error_j: f64,
}

pub(super) struct Candidate {
    pub report: LinearSolidImpulseReport,
    pub cumulative_work_j: f64,
    pub work_correction_j: f64,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn audit(
    mass: &[f64],
    before: &[f64],
    after: &[f64],
    impulse: &[f64],
    allowance: f64,
    kinetic_after_j: f64,
    history: History,
    cx: &Cx<'_>,
) -> Result<Candidate> {
    let n = before.len();
    let mut delta = Vec::with_capacity(n);
    let mut mean = Vec::with_capacity(n);
    for i in 0..n {
        delta.push(finite(after[i] - before[i])?);
        mean.push(finite(0.5 * before[i] + 0.5 * after[i])?);
    }
    let momentum_change_kg_m_s = action(mass, &delta, cx)?;
    let mut impulse_residual_kg_m_s = Vec::with_capacity(n);
    let mut maximum_impulse_residual_kg_m_s = 0.0_f64;
    for i in 0..n {
        let residual = finite(impulse[i] - momentum_change_kg_m_s[i])?;
        maximum_impulse_residual_kg_m_s = maximum_impulse_residual_kg_m_s.max(residual.abs());
        impulse_residual_kg_m_s.push(residual);
    }
    if maximum_impulse_residual_kg_m_s > allowance {
        return Err(LinearSolidImpulseError::Momentum {
            defect_kg_m_s: maximum_impulse_residual_kg_m_s,
            allowance_kg_m_s: allowance,
        });
    }
    let supplied_work_j = dot(&mean, impulse)?;
    let kinetic_change_j = dot(&mean, &momentum_change_kg_m_s)?;
    // Evaluate this from the original residual, not two large work totals.
    let balance_error_j = dot(&mean, &impulse_residual_kg_m_s)?;
    let (cumulative_work_j, work_correction_j) = if supplied_work_j == 0.0 {
        (history.cumulative_work_j, history.work_correction_j)
    } else {
        let increment = finite(supplied_work_j - history.work_correction_j)?;
        let work = finite(history.cumulative_work_j + increment)?;
        (
            work,
            finite((work - history.cumulative_work_j) - increment)?,
        )
    };
    let energy = ExactLinearSolidEnergy {
        stored_before_j: history.stored_j,
        stored_after_j: history.stored_j,
        kinetic_before_j: history.kinetic_before_j,
        kinetic_after_j: finite(kinetic_after_j)?,
        supplied_work_j,
        dissipated_j: 0.0,
        cumulative_dissipated_j: history.cumulative_dissipated_j,
        balance_error_j,
        cumulative_balance_error_j: cumulative_work_j
            - ((history.stored_j + kinetic_after_j) - history.initial_energy_j)
            - history.cumulative_dissipated_j,
    };
    for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
        if !defect_j.is_finite() || defect_j.abs() > history.maximum_work_error_j {
            return Err(LinearSolidImpulseError::Energy {
                defect_j,
                allowance_j: history.maximum_work_error_j,
            });
        }
    }
    poll(cx)?;
    Ok(Candidate {
        report: LinearSolidImpulseReport {
            time_s: history.time_s,
            momentum_change_kg_m_s,
            impulse_residual_kg_m_s,
            maximum_impulse_residual_kg_m_s,
            kinetic_change_j,
            energy,
        },
        cumulative_work_j,
        work_correction_j,
    })
}

impl ExactLinearSolid {
    /// Apply a resolved instantaneous impulse in `free_dofs()` order [kg m/s].
    /// The cached, admitted inverse supplies the velocity response, including
    /// zero-stiffness rigid modes. The original consistent M then checks the
    /// actual rounded change against the explicit component allowance.
    ///
    /// Displacements, accepted steps, clock, viscosity and force preparation
    /// stay fixed. Signed J dot average(v) enters the existing cumulative work;
    /// no constitutive heat or elapsed-time dissipation occurs. Homogeneous
    /// supports remain held; their discarded reaction rows are unavailable.
    /// This is an impulse actuator, not a contact or impact constitutive law.
    pub fn apply_momentum_impulse(
        &mut self,
        impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> Result<LinearSolidImpulseReport> {
        let n = self.free_dofs.len();
        admit(impulse_kg_m_s, n, maximum_impulse_error_kg_m_s, cx)?;
        let mut velocity = action(&self.inverse_mass, impulse_kg_m_s, cx)?;
        for (next, &old) in velocity.iter_mut().zip(self.velocity_mps()) {
            *next = if *next == 0.0 {
                old
            } else {
                finite(old + *next)?
            };
        }
        let candidate = audit(
            &self.mass,
            self.velocity_mps(),
            &velocity,
            impulse_kg_m_s,
            maximum_impulse_error_kg_m_s,
            quadratic(&self.mass, &velocity),
            History {
                time_s: self.time_s(),
                stored_j: self.stored_energy_j(),
                kinetic_before_j: self.kinetic_energy_j(),
                initial_energy_j: self.initial_energy_j,
                cumulative_work_j: self.cumulative_work_j,
                work_correction_j: self.work_correction_j,
                cumulative_dissipated_j: self.cumulative_dissipated_j,
                maximum_work_error_j: self.controls.maximum_work_error_j,
            },
            cx,
        )?;
        self.state[n..].copy_from_slice(&velocity);
        self.cumulative_work_j = candidate.cumulative_work_j;
        self.work_correction_j = candidate.work_correction_j;
        Ok(candidate.report)
    }
}
