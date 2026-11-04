//! Backward-Euler heat residuals on original kinetic crystallization cells.
//!
//! Unknowns remain changes in each original cell's energy, in joules. Every
//! candidate replays the same physical interval from its FROZEN accepted phase
//! histories, then the existing thermal owner evaluates all links, covered
//! boundaries and external powers at those endpoints. Nonlinear iterations
//! never accumulate extra kinetic age, and latent heat is never a second source.
//!
//! All active kinetic materials require explicitly selected phase adaptivity.
//! Ordinary mold/substrate cells retain their original caloric/storage laws.
//! The local phase solver follows the candidate's linear enthalpy path; its
//! tolerances do not control backward-Euler transport time error. This remains
//! fixed reference geometry and fixed cell mass, without a pressure-work,
//! moving-interface or mixture-entropy model inferred from the parent law.

use super::super::backward_euler::{
    BackwardEulerStep, BackwardEulerTrial, BoundaryStepCoefficients, StepCoefficients,
};
use super::*;
use fs_math::dd::Dd;
use std::cell::RefCell;

/// Original Newton-Krylov solve and event-aware advancing of these residuals.
#[path = "implicit/solve.rs"]
pub mod solve;
pub use solve::{
    CrystallizationImplicitControls, CrystallizationImplicitError, CrystallizationImplicitReport,
};

/// One physical interval and immutable accepted phase inventory.
pub struct CrystallizationBackwardEulerStep<'s, 'a> {
    initial: &'s CrystallizationSimulation<'a>,
    thermal: BackwardEulerStep<'s, 'a>,
    controls: CrystallizationControls,
    initial_caloric_defects_j: Vec<f64>,
}

impl<'s, 'a> CrystallizationBackwardEulerStep<'s, 'a> {
    /// Freeze one interval without processing events or changing any history.
    /// Start-time events must already be processed. Interior events refuse;
    /// target-time events remain pending in an accepted left-limit endpoint,
    /// exactly as on the original backward-Euler thermal owner.
    ///
    /// The supplied endpoint determines this single interval. `maximum_dt_s`
    /// is validated here, while an advancing/marching owner applies that cap.
    pub fn new(
        initial: &'s CrystallizationSimulation<'a>,
        target_s: f64,
        controls: CrystallizationControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, CrystallizationTimelineError> {
        checkpoint()?;
        controls.caloric.validate()?;
        if !positive(controls.maximum_dt_s)
            || !nonnegative(controls.cell_energy_tolerance_j)
            || !nonnegative(controls.total_energy_tolerance_j)
        {
            return Err(DepositionError::InvalidInput(
                "invalid implicit crystallization phase or energy controls",
            )
            .into());
        }
        let thermal = BackwardEulerStep::new(&initial.thermal, target_s, &mut checkpoint)?;
        let mut initial_caloric_defects_j = Vec::with_capacity(thermal.cell_ids().len());
        for &id in thermal.cell_ids() {
            checkpoint()?;
            let cell = &initial.thermal.cells[&id];
            let defect = if let Some(material) = initial.materials.get(&cell.material_id) {
                let adaptivity = initial.adaptivity.get(&cell.material_id).ok_or(
                    DepositionError::InvalidInput(
                        "implicit crystallization requires selected adaptivity for every active kinetic material",
                    ),
                )?;
                adaptivity.validate()?;
                let state = initial
                    .states
                    .get(&id)
                    .ok_or(CrystallizationError::StateMismatch)?;
                if state.specific_enthalpy_j_kg() != cell.specific_h_j_kg {
                    return Err(CrystallizationError::StateMismatch.into());
                }
                material.law.validate_state(*state, controls.caloric)?;
                let reconstructed = material.law.specific_enthalpy_j_kg(
                    state.temperature_k(),
                    state.crystalline_mass_fraction(),
                )?;
                caloric_joules(
                    cell.mass_kg,
                    finite(reconstructed - cell.specific_h_j_kg)?.abs(),
                )?
            } else {
                0.0
            };
            initial_caloric_defects_j.push(defect);
        }
        checkpoint()?;
        Ok(Self {
            initial,
            thermal,
            controls,
            initial_caloric_defects_j,
        })
    }

    pub fn cell_ids(&self) -> &[CellId] {
        self.thermal.cell_ids()
    }
    pub fn duration_s(&self) -> f64 {
        self.thermal.duration_s()
    }
    pub fn target_s(&self) -> f64 {
        self.thermal.target_s()
    }

    /// The original energy-coordinate resolution projection, for solver iterates
    /// only. Every source/link joule still enters the actual endpoint residual.
    pub fn resolvable_increments_j(&self, energy_increments_j: &[f64]) -> Vec<f64> {
        self.thermal.resolvable_increments_j(energy_increments_j)
    }

    /// Evaluate one complete original-state proposal without publishing it.
    /// The source callback sees ALL actual kinetic endpoints simultaneously and
    /// may supply original source/link/Robin conductances. It must remain pure
    /// with respect to physical history across speculative evaluations.
    /// No derivative, phase mixture conductivity or entropy chart is invented.
    pub fn evaluate<T>(
        &self,
        energy_increments_j: &[f64],
        resolve: impl FnOnce(
            &CrystallizationSimulation<'a>,
        )
            -> Result<BoundaryStepCoefficients<T>, CrystallizationTimelineError>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<CrystallizationBackwardEulerTrial<'a, T>, CrystallizationTimelineError> {
        let checkpoint = RefCell::new(checkpoint);
        let trial = self.thermal.evaluate_with_state_calorics(
            energy_increments_j,
            |thermal| self.endpoint(thermal, resolve, &mut || (checkpoint.borrow_mut())()),
            || (checkpoint.borrow_mut())(),
        )?;
        let endpoint = trial.payload();
        let mut local = Dd::ZERO;
        let mut signed = Dd::ZERO;
        let mut maximum: f64 = 0.0;
        for (&residual, &caloric) in trial.residual_j().iter().zip(&endpoint.caloric_defects_j) {
            (checkpoint.borrow_mut())()?;
            let defect = finite(residual.abs() + caloric)?;
            local = local + Dd::from_f64(defect);
            signed = signed + Dd::from_f64(residual);
            maximum = maximum.max(defect);
        }
        // Independent cell-relative storage reduction detects heat aggregation
        // and stored-total roundoff, without subtracting large formation totals.
        let thermal = *trial.audit();
        let physical_residual = endpoint.stored_change_j - Dd::from_f64(thermal.external_heat_j);
        let storage_roundoff =
            (Dd::from_f64(thermal.stored_energy_change_j) - endpoint.stored_change_j).abs();
        let exchange_roundoff = (physical_residual - signed).abs();
        let absolute = finite((local + storage_roundoff + exchange_roundoff).to_f64())?
            .max(finite(physical_residual.to_f64())?.abs())
            .max(thermal.residual_j.abs());
        let audit = CrystallizationAudit {
            thermal,
            crystallized_mass_kg: endpoint.crystallized_mass_kg,
            maximum_cell_energy_residual_j: maximum,
            absolute_energy_residual_j: absolute,
            adaptive_phase: endpoint.adaptive_phase,
        };
        (checkpoint.borrow_mut())()?;
        Ok(CrystallizationBackwardEulerTrial {
            thermal: trial,
            audit,
            controls: self.controls,
        })
    }

    #[allow(clippy::type_complexity)]
    fn endpoint<T>(
        &self,
        thermal: &DepositionSimulation<'a>,
        resolve: impl FnOnce(
            &CrystallizationSimulation<'a>,
        )
            -> Result<BoundaryStepCoefficients<T>, CrystallizationTimelineError>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (
            BTreeMap<CellId, (f64, f64)>,
            BoundaryStepCoefficients<Endpoint<'a, T>>,
        ),
        CrystallizationTimelineError,
    > {
        checkpoint()?;
        let mut next = CrystallizationSimulation {
            thermal: thermal.clone(),
            materials: self.initial.materials.clone(),
            states: BTreeMap::new(),
            remelting: self.initial.remelting.clone(),
            adaptivity: self.initial.adaptivity.clone(),
        };
        let mut calorics = BTreeMap::new();
        let mut caloric_defects_j = Vec::with_capacity(self.cell_ids().len());
        let mut stored_change_j = Dd::ZERO;
        let mut crystals = Dd::ZERO;
        let mut adaptive_phase = CrystallizationAdaptiveAudit::default();
        for (i, &id) in self.cell_ids().iter().enumerate() {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let old = &self.initial.thermal.cells[&id];
            let represented =
                finite(old.mass_kg * finite(cell.specific_h_j_kg - old.specific_h_j_kg)?)?;
            stored_change_j = stored_change_j + Dd::from_f64(represented);
            let (temperature, slope, final_caloric) =
                if let Some(material) = self.initial.materials.get(&cell.material_id) {
                    let previous = self.initial.states[&id];
                    let mut refused = None;
                    let result = material.law.advance_specific_enthalpy_adaptive(
                        previous,
                        cell.specific_h_j_kg,
                        self.duration_s(),
                        self.controls.caloric,
                        self.initial.adaptivity[&cell.material_id],
                        self.initial.remelting.get(&cell.material_id).copied(),
                        || match checkpoint() {
                            Ok(()) => false,
                            Err(error) => {
                                refused = Some(error);
                                true
                            }
                        },
                    );
                    if let Some(error) = refused {
                        return Err(error.into());
                    }
                    let step = result?;
                    adaptive_phase.include(&step)?;
                    crystals = crystals
                        + Dd::from_f64(finite(
                            cell.mass_kg
                                * (step.state.crystalline_mass_fraction()
                                    - previous.crystalline_mass_fraction()),
                        )?);
                    next.states.insert(id, step.state);
                    (
                        step.state.temperature_k(),
                        material.maximum_temperature_slope,
                        caloric_joules(cell.mass_kg, step.caloric_residual_j_kg.abs())?,
                    )
                } else {
                    (
                        thermal.temperature_k(id)?,
                        thermal.material_temperature_slopes[&cell.material_id],
                        0.0,
                    )
                };
            calorics.insert(id, (temperature, slope));
            caloric_defects_j.push(finite(self.initial_caloric_defects_j[i] + final_caloric)?);
        }
        finite(stored_change_j.to_f64())?;
        let crystallized_mass_kg = finite(crystals.to_f64())?;
        checkpoint()?;
        let BoundaryStepCoefficients {
            bulk:
                StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload,
                },
            boundaries,
        } = resolve(&next)?;
        checkpoint()?;
        Ok((
            calorics,
            BoundaryStepCoefficients {
                bulk: StepCoefficients {
                    sources_w,
                    conductances_w_k,
                    payload: Endpoint {
                        next,
                        payload,
                        caloric_defects_j,
                        stored_change_j,
                        crystallized_mass_kg,
                        adaptive_phase,
                    },
                },
                boundaries,
            },
        ))
    }
}

struct Endpoint<'a, T> {
    next: CrystallizationSimulation<'a>,
    payload: T,
    caloric_defects_j: Vec<f64>,
    stored_change_j: Dd,
    crystallized_mass_kg: f64,
    adaptive_phase: CrystallizationAdaptiveAudit,
}

/// One tentative residual and the original phase inventory that generated it.
pub struct CrystallizationBackwardEulerTrial<'a, T> {
    thermal: BackwardEulerTrial<'a, Endpoint<'a, T>>,
    audit: CrystallizationAudit,
    controls: CrystallizationControls,
}

impl<'a, T> CrystallizationBackwardEulerTrial<'a, T> {
    pub fn residual_j(&self) -> &[f64] {
        self.thermal.residual_j()
    }
    /// Mathematical heat-residual maximum for the nonlinear solver. Admission
    /// also checks the separate caloric and noncancelling defects in `audit()`.
    pub fn maximum_cell_residual_j(&self) -> f64 {
        self.thermal.maximum_cell_residual_j()
    }
    pub fn payload(&self) -> &T {
        &self.thermal.payload().payload
    }
    pub const fn audit(&self) -> &CrystallizationAudit {
        &self.audit
    }

    /// Accept original enthalpy, all phase histories and physical time together.
    /// The controls frozen at step construction gate every local caloric/heat
    /// defect and their noncancelling total, in addition to the original thermal
    /// owner's local and global residual gates. Endpoint events remain pending.
    pub fn accept(
        self,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (CrystallizationSimulation<'a>, T, CrystallizationAudit),
        CrystallizationTimelineError,
    > {
        checkpoint()?;
        if self.audit.maximum_cell_energy_residual_j > self.controls.cell_energy_tolerance_j {
            return Err(
                DepositionError::EnergyBalance(self.audit.maximum_cell_energy_residual_j).into(),
            );
        }
        if self.audit.absolute_energy_residual_j > self.controls.total_energy_tolerance_j {
            return Err(
                DepositionError::EnergyBalance(self.audit.absolute_energy_residual_j).into(),
            );
        }
        let (thermal, endpoint, _) = self.thermal.accept(
            self.controls.cell_energy_tolerance_j,
            self.controls.total_energy_tolerance_j,
            &mut checkpoint,
        )?;
        let mut next = endpoint.next;
        next.thermal = thermal;
        checkpoint()?;
        Ok((next, endpoint.payload, self.audit))
    }
}

fn caloric_joules(mass_kg: f64, residual_j_kg: f64) -> Result<f64, DepositionError> {
    let result = finite(mass_kg * residual_j_kg)?;
    if mass_kg != 0.0 && residual_j_kg != 0.0 && result == 0.0 {
        return Err(DepositionError::Arithmetic);
    }
    Ok(result)
}
