//! Bounded sparse Newton solution of the original material-backed heat step.

use super::{LayeredThermalLinearization, LayeredThermalStep};
use crate::thermal::transport::stack::transient::MAX_LAYERED_THERMAL_CELLS;
use fs_la::tridiagonal::{TridiagonalError, TridiagonalWorkspace};
use fs_thermal::enthalpy::{
    DepositionError, DepositionSimulation, EnergyAudit, entropy::EntropyAudit,
};

/// Nonlinear work and independent local/whole-body energy acceptance limits.
/// These control algebraic error; temporal and spatial errors need refinement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredThermalSolveControls {
    pub maximum_newton_iterations: usize,
    /// Number of halvings after an unsuccessful full Newton proposal.
    pub maximum_backtracks: usize,
    /// Includes the initial residual, rejected trials and final entropy query.
    pub maximum_residual_evaluations: usize,
    pub cell_tolerance_j: f64,
    pub total_tolerance_j: f64,
}

impl Default for LayeredThermalSolveControls {
    fn default() -> Self {
        Self {
            maximum_newton_iterations: 32,
            maximum_backtracks: 24,
            maximum_residual_evaluations: 256,
            cell_tolerance_j: 1e-8,
            total_tolerance_j: 1e-8,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum LayeredThermalSolveError {
    Thermal(DepositionError),
    Linear(TridiagonalError),
    NotConverged {
        maximum_cell_residual_j: f64,
        total_residual_j: f64,
    },
}
impl From<DepositionError> for LayeredThermalSolveError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<TridiagonalError> for LayeredThermalSolveError {
    fn from(error: TridiagonalError) -> Self {
        Self::Linear(error)
    }
}
impl core::fmt::Display for LayeredThermalSolveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit layered thermal solve: {self:?}")
    }
}
impl std::error::Error for LayeredThermalSolveError {}

/// Accounts describe only the accepted endpoint; work counts include rejections.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayeredThermalSolveReport {
    pub energy: EnergyAudit,
    pub entropy: Option<EntropyAudit>,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub backtracks: usize,
    pub maximum_cell_residual_j: f64,
}

/// Reusable O(cells) Newton vectors and adjacent-pivot tridiagonal factors.
/// The existing material/residual evaluation still allocates O(cells) trials.
/// No dense matrix, independent thermal state or numerical material fit is kept.
#[derive(Debug)]
pub struct LayeredThermalWorkspace {
    linear: TridiagonalWorkspace,
    energy: Vec<f64>,
    proposal: Vec<f64>,
    direction: Vec<f64>,
    rhs: Vec<f64>,
}
impl LayeredThermalWorkspace {
    pub fn new(cells: usize) -> Result<Self, LayeredThermalSolveError> {
        if cells == 0 || cells > MAX_LAYERED_THERMAL_CELLS {
            return Err(DepositionError::Budget("layered thermal workspace cells").into());
        }
        let allocate = || {
            let mut values = Vec::new();
            values
                .try_reserve_exact(cells)
                .map_err(|_| TridiagonalError::Capacity)?;
            values.resize(cells, 0.0);
            Ok::<_, TridiagonalError>(values)
        };
        Ok(Self {
            linear: TridiagonalWorkspace::new(cells)?,
            energy: allocate()?,
            proposal: allocate()?,
            direction: allocate()?,
            rhs: allocate()?,
        })
    }
    pub fn cells(&self) -> usize {
        self.energy.len()
    }
}

impl<'m> LayeredThermalStep<'_, 'm, '_> {
    /// Solve and accept this complete backward-Euler interval, without an
    /// explicit diffusion CFL restriction or silently subdividing physical time.
    /// The exact source-integrated flux Jacobian is reevaluated at each Newton
    /// iterate. Backtracking enforces the original material domain and residual
    /// decrease; no out-of-domain conductivity or capacity is extrapolated.
    /// All failures leave the original timeline unchanged; scratch may change.
    /// Events exactly at the target remain pending for the next owner action.
    pub fn solve(
        &self,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            LayeredThermalLinearization,
            LayeredThermalSolveReport,
        ),
        LayeredThermalSolveError,
    > {
        self.solve_inner(controls, workspace, false, checkpoint)
    }

    /// Same solve, with the original caloric/heat owner's endpoint entropy audit.
    /// Entropy is queried only at the converged candidate and requires the
    /// source laws' explicit sensible-entropy opt-in. It is never guessed.
    pub fn solve_with_entropy(
        &self,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            LayeredThermalLinearization,
            LayeredThermalSolveReport,
        ),
        LayeredThermalSolveError,
    > {
        self.solve_inner(controls, workspace, true, checkpoint)
    }

    fn solve_inner(
        &self,
        controls: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        with_entropy: bool,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            DepositionSimulation<'m>,
            LayeredThermalLinearization,
            LayeredThermalSolveReport,
        ),
        LayeredThermalSolveError,
    > {
        checkpoint()?;
        if workspace.cells() != self.cell_ids().len()
            || controls.maximum_newton_iterations == 0
            || controls.maximum_residual_evaluations == 0
            || ![controls.cell_tolerance_j, controls.total_tolerance_j]
                .into_iter()
                .all(|x| x.is_finite() && x >= 0.0)
        {
            return Err(DepositionError::InvalidInput(
                "invalid layered thermal solve controls or workspace",
            )
            .into());
        }
        let mut evaluations = 0;
        let mut evaluate =
            |energy: &[f64],
             entropy: bool,
             checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>| {
                if evaluations >= controls.maximum_residual_evaluations {
                    return Err(DepositionError::Budget(
                        "layered thermal residual evaluations",
                    ));
                }
                evaluations += 1;
                if entropy {
                    self.evaluate_with_entropy(energy, checkpoint)
                } else {
                    self.evaluate(energy, checkpoint)
                }
            };
        workspace.energy.fill(0.0);
        let mut trial = evaluate(&workspace.energy, false, &mut checkpoint)?;
        let mut backtracks = 0;
        for iteration in 0..=controls.maximum_newton_iterations {
            checkpoint()?;
            let maximum = trial.maximum_cell_residual_j();
            let total = trial.audit().residual_j;
            let converged = |local: f64, total: f64| {
                local <= controls.cell_tolerance_j && total.abs() <= controls.total_tolerance_j
            };
            if converged(maximum, total) {
                if with_entropy {
                    trial = evaluate(&workspace.energy, true, &mut checkpoint)?;
                }
                let maximum_cell_residual_j = trial.maximum_cell_residual_j();
                let entropy = trial.entropy_audit().copied();
                let (next, linearization, energy) = trial.accept(
                    controls.cell_tolerance_j,
                    controls.total_tolerance_j,
                    &mut checkpoint,
                )?;
                return Ok((
                    next,
                    linearization,
                    LayeredThermalSolveReport {
                        energy,
                        entropy,
                        newton_iterations: iteration,
                        residual_evaluations: evaluations,
                        backtracks,
                        maximum_cell_residual_j,
                    },
                ));
            }
            let not_converged = || LayeredThermalSolveError::NotConverged {
                maximum_cell_residual_j: maximum,
                total_residual_j: total,
            };
            if iteration == controls.maximum_newton_iterations {
                return Err(not_converged());
            }
            for (rhs, &residual) in workspace.rhs.iter_mut().zip(trial.residual_j()) {
                *rhs = -residual;
            }
            let jacobian = trial.payload();
            let mut cancellation_error = None;
            let solved = workspace.linear.solve_into(
                jacobian.lower_diagonal(),
                jacobian.diagonal(),
                jacobian.upper_diagonal(),
                &workspace.rhs,
                &mut workspace.direction,
                || match checkpoint() {
                    Ok(()) => false,
                    Err(error) => {
                        cancellation_error = Some(error);
                        true
                    }
                },
            );
            if let Some(error) = cancellation_error {
                return Err(error.into());
            }
            solved?;
            let mut fraction = 1.0;
            let mut accepted = None;
            for attempt in 0..=controls.maximum_backtracks {
                checkpoint()?;
                let mut representable = false;
                for (i, &id) in self.cell_ids().iter().enumerate() {
                    let old = &self.initial.cells()[&id];
                    let mut value = workspace.energy[i] + fraction * workspace.direction[i];
                    // Tiny far-field increments may not change the actual
                    // storage coordinate. Submit their representable zero;
                    // the original residual still charges every heat flux.
                    if old.specific_h_j_kg + value / old.mass_kg == old.specific_h_j_kg {
                        value = 0.0;
                    }
                    representable |= value != workspace.energy[i];
                    workspace.proposal[i] = value;
                }
                if representable {
                    match evaluate(&workspace.proposal, false, &mut checkpoint) {
                        Ok(candidate) => {
                            let residual = candidate.maximum_cell_residual_j();
                            if residual <= (1.0 - 1e-4 * fraction) * maximum
                                || converged(residual, candidate.audit().residual_j)
                            {
                                accepted = Some(candidate);
                                break;
                            }
                        }
                        Err(DepositionError::MaterialDomain | DepositionError::Arithmetic) => {}
                        Err(error) => return Err(error.into()),
                    }
                }
                if attempt == controls.maximum_backtracks {
                    break;
                }
                backtracks += 1;
                fraction *= 0.5;
                if fraction == 0.0 {
                    break;
                }
            }
            let Some(candidate) = accepted else {
                return Err(not_converged());
            };
            workspace.energy.copy_from_slice(&workspace.proposal);
            trial = candidate;
        }
        unreachable!("bounded Newton loop returns at its limit")
    }
}
