//! Finite-rate dissolution and precipitation in one shared finite thermal bath.
//!
//! In the independent ideal limit, each stage uses the original exact fixed-T
//! kinetic trajectory at its unknown
//! endpoint temperature. A scalar first-law solve closes that temperature with
//! the SAME solubility law, dissolution enthalpy and phase Cp sources. On an
//! active branch its derivative is positive phase-weighted sensible capacity
//! plus `w * n_sat * delta_H(T)^2 / (R*T^2)`. Here `w` is the nonnegative
//! sensitivity of final dissolved amount to saturation capacity. Fixed area
//! gives `w = 1-exp(-G*h/V)`. For particle area, scalar comparison of the
//! original kinetic ODE preserves `w >= 0`: its sensitivity has positive
//! forcing `G(s)/V` and starts at zero. No frozen-area derivative is assumed.
//! Exhaustion and the declared growth policy retain their finite-stock limits.
//!
//! The thermal integration is FIRST ORDER. Bounded full-step/two-half-step
//! comparison controls temporal disagreement; the accepted state uses the two
//! half steps without extrapolating phase stocks. This estimate is not a rigorous
//! error bound. A frozen-stage exhaustion time is not a continuous thermal event
//! time and is deliberately not exposed as one.
//!
//! Solvent volume, reference conductance and growth policy stay fixed over a window.
//! Selected particle area evolves through each original isothermal phase solve;
//! its original reference stock is never reset at a temperature trial.
//! Specified heat is delivered uniformly during that window. Alternatively a
//! prescribed-temperature reservoir exchanges heat through a finite thermal
//! conductance, using the same implicit stages and accepted heat account. There is no
//! simultaneous transport, evolving film coefficient, nucleation, fragmentation
//! or population balance. An explicitly bound [`SolubilityThermalInteractionV1`]
//! instead supplies conditional activity equilibrium and its consistent excess
//! calorics. Its phase stage uses the CURRENT temperature before a frozen-stock
//! heat projection. This first-order split retains the same native kinetics,
//! state, adaptive controller and whole-window energy account, without assuming
//! the ideal endpoint-phase monotonicity proof applies to interacting solutes.

use super::{
    MAX_SOLUBILITY_CALORIC_COMPONENTS_V1, PhaseIncrement, SolubilityCaloricErrorV1,
    SolubilityCaloricModelV1, SolubilityCaloricSolveV1, checked, poll, positive,
};
use crate::solution::solubility::kinetics::{SolubilityKineticStateV1, SolubilityKineticsV1};
use fs_qty::{Amount, Temperature, Time, Volume};
use std::sync::Arc;

mod heat_transfer;
use heat_transfer::{HeatReceipt, StageHeat, StageResult, WindowHeat};

mod finite_reservoir;
pub use finite_reservoir::{SensibleHeatReservoirV1, SolubilityThermalReservoirOutcomeV1};

pub mod interaction;
pub use interaction::SolubilityThermalInteractionV1;

#[derive(Clone, Debug)]
struct ThermalKineticOwner {
    caloric: SolubilityCaloricModelV1,
    kinetics: Vec<SolubilityKineticsV1>,
    interactions: Option<Arc<dyn SolubilityThermalInteractionV1>>,
}
impl PartialEq for ThermalKineticOwner {
    fn eq(&self, other: &Self) -> bool {
        self.caloric == other.caloric
            && self.kinetics == other.kinetics
            && match (&self.interactions, &other.interactions) {
                (None, None) => true,
                (Some(left), Some(right)) => Arc::ptr_eq(left, right),
                _ => false,
            }
    }
}

/// One immutable original caloric/kinetic source binding for a common bath.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityThermalKineticsV1 {
    owner: Arc<ThermalKineticOwner>,
}
impl SolubilityThermalKineticsV1 {
    /// Every ordered kinetic law must equal its original caloric law, including
    /// attached phase curves. Matching species names alone is insufficient.
    pub fn new(
        caloric: SolubilityCaloricModelV1,
        kinetics: Vec<SolubilityKineticsV1>,
    ) -> Result<Self, SolubilityCaloricErrorV1> {
        if kinetics.len() != caloric.components().len()
            || kinetics
                .iter()
                .zip(caloric.components())
                .any(|(kinetic, thermal)| kinetic.law() != thermal.law())
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal kinetics require the same ordered original solubility laws",
            ));
        }
        Ok(Self {
            owner: Arc::new(ThermalKineticOwner {
                caloric,
                kinetics,
                interactions: None,
            }),
        })
    }

    /// Bind explicit interaction thermodynamics and its excess calorics to the
    /// SAME original laws and thermal integrator. Preparation must qualify a
    /// positive frozen-stock heat capacity over the source's temperature support.
    /// The original fixed-temperature interaction fit alone is not enough.
    ///
    /// The source Arc remains part of this owner and every state/outcome. A
    /// different interaction source cannot replace an existing binding.
    pub fn with_interactions(
        &self,
        interactions: Arc<dyn SolubilityThermalInteractionV1>,
    ) -> Result<Self, SolubilityCaloricErrorV1> {
        if let Some(original) = &self.owner.interactions {
            return if Arc::ptr_eq(original, &interactions) {
                Ok(self.clone())
            } else {
                Err(SolubilityCaloricErrorV1::Invalid(
                    "thermal interaction source cannot replace an existing binding",
                ))
            };
        }
        let domain = interactions.temperature_domain();
        let caloric = self.caloric().temperature_domain();
        if !positive(domain[0].value())
            || !positive(domain[1].value())
            || domain[0].value() >= domain[1].value()
            || domain[0].value().max(caloric[0].value())
                >= domain[1].value().min(caloric[1].value())
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal interaction source needs overlapping positive temperature support",
            ));
        }
        interactions.admit_binding(self.caloric(), self.kinetics())?;
        Ok(Self {
            owner: Arc::new(ThermalKineticOwner {
                caloric: self.caloric().clone(),
                kinetics: self.kinetics().to_vec(),
                interactions: Some(interactions),
            }),
        })
    }

    /// Exact immutable interaction source, if this owner explicitly selected one.
    pub fn interaction_source(&self) -> Option<&Arc<dyn SolubilityThermalInteractionV1>> {
        self.owner.interactions.as_ref()
    }

    /// Common caloric and interaction-temperature support; no extrapolation.
    pub fn temperature_domain(&self) -> [Temperature; 2] {
        let mut domain = self.caloric().temperature_domain();
        if let Some(interactions) = &self.owner.interactions {
            let support = interactions.temperature_domain();
            domain[0] = Temperature::new(domain[0].value().max(support[0].value()));
            domain[1] = Temperature::new(domain[1].value().min(support[1].value()));
        }
        domain
    }

    pub fn caloric(&self) -> &SolubilityCaloricModelV1 {
        &self.owner.caloric
    }
    pub fn kinetics(&self) -> &[SolubilityKineticsV1] {
        &self.owner.kinetics
    }

    /// Preserve the original material volume and BOTH actual phase stocks.
    /// Preparation neither equilibrates the phases nor introduces heat.
    pub fn prepare_at_solvent_volume(
        &self,
        temperature: Temperature,
        solvent_volume: Volume,
        phases: &[SolubilityKineticStateV1],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityThermalKineticStateV1, SolubilityCaloricErrorV1> {
        poll(&mut cancelled)?;
        let domain = self.temperature_domain();
        let moles = solvent_volume.value() / self.caloric().solvent_molar_volume_m3_mol();
        if !temperature.value().is_finite()
            || !(domain[0].value()..=domain[1].value()).contains(&temperature.value())
            || !positive(solvent_volume.value())
            || !positive(moles)
            || phases.len() != self.kinetics().len()
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal kinetic temperature, solvent volume or phase layout",
            ));
        }
        for (index, (kinetic, phase)) in self.kinetics().iter().zip(phases).enumerate() {
            poll(&mut cancelled)?;
            if let Some(interactions) = &self.owner.interactions {
                kinetic.advance_at_saturation(
                    temperature,
                    solvent_volume,
                    phase,
                    Time::new(0.0),
                    interactions.conditional_saturation(
                        temperature,
                        solvent_volume,
                        phases,
                        index,
                    )?,
                )?;
            } else {
                kinetic.advance(temperature, solvent_volume, phase, Time::new(0.0))?;
            }
        }
        let result = SolubilityThermalKineticStateV1 {
            owner: Arc::clone(&self.owner),
            temperature,
            solvent_moles: Amount::new(moles),
            solvent_volume,
            phases: Arc::from(phases),
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Direct original-to-final energy increment, with no reference-enthalpy
    /// subtraction and no accumulation of individual stage receipts.
    pub fn enthalpy_change_j(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        final_state: &SolubilityThermalKineticStateV1,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        if !self.owns(initial)
            || !self.owns(final_state)
            || initial.solvent_moles != final_state.solvent_moles
            || initial.solvent_volume != final_state.solvent_volume
            || initial
                .phases
                .iter()
                .zip(final_state.phases.iter())
                .any(|(a, b)| a.total_solute() != b.total_solute())
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal kinetic energy requires the same source owner and inventories",
            ));
        }
        self.change(
            initial,
            final_state.temperature,
            &final_state.phases,
            &mut || false,
        )
    }

    /// Integrate a closed finite bath with uniform external heating/cooling.
    ///
    /// Every stage is anchored to the prescribed cumulative window heat using
    /// direct original-to-current energy. This prevents stage root tolerances
    /// from accumulating as heat drift. Final accounting independently uses
    /// actual initial/final stocks, not the accepted stages' transfer receipts.
    /// Trials allocate nothing; cancellation or any refusal leaves input intact.
    pub fn advance(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        duration: Time,
        heat_input_j: f64,
        controls: SolubilityThermalKineticControlsV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityThermalKineticOutcomeV1, SolubilityCaloricErrorV1> {
        self.advance_with_source(
            initial,
            duration,
            WindowHeat::Uniform(heat_input_j),
            controls,
            cancelled,
        )
        .map(|(balance, _)| balance)
    }

    fn advance_with_source(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        duration: Time,
        heat_source: WindowHeat,
        controls: SolubilityThermalKineticControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(SolubilityThermalKineticOutcomeV1, HeatReceipt), SolubilityCaloricErrorV1> {
        poll(&mut cancelled)?;
        self.admit(
            initial,
            duration,
            match heat_source {
                WindowHeat::Uniform(heat) => heat,
                WindowHeat::Reservoir(_) | WindowHeat::Finite(_) => 0.0,
            },
            controls,
        )?;
        let count = self.kinetics().len();
        let mut current = WorkingState::from(initial);
        let mut accepted_heat = HeatReceipt::default();
        let mut elapsed = 0.0;
        let mut interval = controls.initial_step.value();
        let mut accepted_steps = 0;
        let mut rejected_steps = 0;
        let mut root_iterations = 0;
        let mut maximum_accepted_error_norm: f64 = 0.0;
        let mut root = controls.root;
        // Reserve room for independent final-state rounding in the window's
        // requested energy tolerance. No tolerance is added to the heat input.
        root.energy_tolerance_j *= 0.25;
        while elapsed < duration.value() {
            poll(&mut cancelled)?;
            if accepted_steps + rejected_steps == controls.maximum_attempts {
                return Err(SolubilityCaloricErrorV1::IntegrationLimit);
            }
            let remaining = duration.value() - elapsed;
            let h = interval.min(remaining);
            let end = if h == remaining {
                duration.value()
            } else {
                elapsed + h
            };
            let h = end - elapsed;
            if h <= 0.0 || 0.5 * h == 0.0 {
                return Err(SolubilityCaloricErrorV1::Arithmetic(
                    "thermal kinetic time increment is unrepresentable",
                ));
            }
            let (full_heat, middle_heat) = match heat_source {
                WindowHeat::Uniform(heat_input_j) => {
                    // Preserve the original prescribed-heat arithmetic exactly.
                    let target_heat = heat_input_j * (end / duration.value());
                    let start_heat = heat_input_j * (elapsed / duration.value());
                    let middle_heat = 0.5 * start_heat + 0.5 * target_heat;
                    (
                        StageHeat::Prescribed(target_heat),
                        StageHeat::Prescribed(middle_heat),
                    )
                }
                WindowHeat::Reservoir(source) => {
                    let heat = StageHeat::Reservoir {
                        source,
                        prefix: accepted_heat,
                    };
                    (heat, heat)
                }
                WindowHeat::Finite(source) => {
                    let heat = StageHeat::Finite {
                        source,
                        prefix: accepted_heat,
                    };
                    (heat, heat)
                }
            };
            let trial = (|| -> Result<_, SolubilityCaloricErrorV1> {
                let full = self.stage(
                    initial,
                    &current,
                    Time::new(h),
                    full_heat,
                    root,
                    &mut root_iterations,
                    &mut cancelled,
                )?;
                let middle = self.stage(
                    initial,
                    &current,
                    Time::new(0.5 * h),
                    middle_heat,
                    root,
                    &mut root_iterations,
                    &mut cancelled,
                )?;
                let fine_heat = match heat_source {
                    WindowHeat::Uniform(_) => full_heat,
                    WindowHeat::Reservoir(source) => StageHeat::Reservoir {
                        source,
                        prefix: middle.heat,
                    },
                    WindowHeat::Finite(source) => StageHeat::Finite {
                        source,
                        prefix: middle.heat,
                    },
                };
                let fine = self.stage(
                    initial,
                    &middle.state,
                    Time::new(h - 0.5 * h),
                    fine_heat,
                    root,
                    &mut root_iterations,
                    &mut cancelled,
                )?;
                Ok((full, fine))
            })();
            let (full, fine) = match trial {
                Ok(pair) => pair,
                Err(error @ SolubilityCaloricErrorV1::HeatOutsideBracket { .. }) => {
                    rejected_steps += 1;
                    interval = match half_interval(h, controls.minimum_step.value()) {
                        Some(smaller) => smaller,
                        None => return Err(error),
                    };
                    continue;
                }
                Err(error) => return Err(error),
            };
            let mut error = temporal_error(&full.state, &fine.state, count, controls);
            if let WindowHeat::Finite(source) = heat_source {
                error = error.max(source.temporal_error(full.heat, fine.heat, controls)?);
            }
            if error > 1.0 {
                rejected_steps += 1;
                interval = half_interval(h, controls.minimum_step.value()).ok_or(
                    SolubilityCaloricErrorV1::TemporalTolerance {
                        step_seconds: h,
                        estimated_error_norm: error,
                    },
                )?;
                continue;
            }
            current = fine.state;
            accepted_heat = fine.heat;
            elapsed = end;
            accepted_steps += 1;
            maximum_accepted_error_norm = maximum_accepted_error_norm.max(error);
            interval = if error <= 0.25 {
                if h >= 0.5 * controls.maximum_step.value() {
                    controls.maximum_step.value()
                } else {
                    h + h
                }
            } else {
                h
            };
        }
        let enthalpy_change_j = self.change(
            initial,
            current.temperature,
            &current.phases[..count],
            &mut cancelled,
        )?;
        let (heat_input_j, energy_residual_j) = match heat_source {
            WindowHeat::Uniform(heat_input_j) => (
                heat_input_j,
                checked(
                    enthalpy_change_j - heat_input_j,
                    "thermal kinetic window residual",
                )?,
            ),
            WindowHeat::Reservoir(_) | WindowHeat::Finite(_) => (
                accepted_heat.value()?,
                accepted_heat.residual(enthalpy_change_j)?,
            ),
        };
        if energy_residual_j.abs() > controls.root.energy_tolerance_j {
            return Err(SolubilityCaloricErrorV1::EnergyResidual {
                residual_j: energy_residual_j,
            });
        }
        let final_state = if accepted_steps == 0 {
            initial.clone()
        } else {
            SolubilityThermalKineticStateV1 {
                owner: Arc::clone(&self.owner),
                temperature: current.temperature,
                solvent_moles: initial.solvent_moles,
                solvent_volume: initial.solvent_volume,
                phases: Arc::from(&current.phases[..count]),
            }
        };
        poll(&mut cancelled)?;
        Ok((
            SolubilityThermalKineticOutcomeV1 {
                final_state,
                duration,
                heat_input_j,
                enthalpy_change_j,
                energy_residual_j,
                accepted_steps,
                rejected_steps,
                root_iterations,
                maximum_accepted_error_norm,
            },
            accepted_heat,
        ))
    }

    fn owns(&self, state: &SolubilityThermalKineticStateV1) -> bool {
        Arc::ptr_eq(&self.owner, &state.owner) || self.owner == state.owner
    }

    fn admit(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        duration: Time,
        heat: f64,
        controls: SolubilityThermalKineticControlsV1,
    ) -> Result<(), SolubilityCaloricErrorV1> {
        let [lo, hi] = controls.root.temperature_bracket.map(|t| t.value());
        let domain = self.temperature_domain();
        let relative = controls.relative_tolerance;
        let amount = controls.amount_absolute_tolerance_mol;
        let temperature = controls.temperature_absolute_tolerance_k;
        if !self.owns(initial)
            || !duration.value().is_finite()
            || duration.value() < 0.0
            || !heat.is_finite()
            || (duration.value() == 0.0 && heat != 0.0)
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo < domain[0].value()
            || hi > domain[1].value()
            || !(lo..=hi).contains(&initial.temperature.value())
            || !positive(controls.root.temperature_tolerance_k)
            || !controls.root.energy_tolerance_j.is_finite()
            || controls.root.energy_tolerance_j < 0.0
            || !(1..=256).contains(&controls.root.maximum_iterations)
            || !positive(controls.minimum_step.value())
            || !positive(controls.maximum_step.value())
            || !positive(controls.initial_step.value())
            || controls.minimum_step.value() > controls.initial_step.value()
            || controls.initial_step.value() > controls.maximum_step.value()
            || !(0.0..=1.0).contains(&relative)
            || !amount.is_finite()
            || amount < 0.0
            || !temperature.is_finite()
            || temperature < 0.0
            || (relative == 0.0 && (amount == 0.0 || temperature == 0.0))
            || !(1..=65_536).contains(&controls.maximum_attempts)
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "thermal kinetic owner, duration, heat or integration controls",
            ));
        }
        Ok(())
    }

    fn change(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        temperature: Temperature,
        phases: &[SolubilityKineticStateV1],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        let base = self.caloric().phase_enthalpy_increment(
            initial.temperature,
            initial.solvent_moles,
            temperature,
            |index| {
                let old = &initial.phases[index];
                Ok(PhaseIncrement {
                    total: old.total_solute().value(),
                    solid: old.solid().value(),
                    dissolved: old.dissolved().value(),
                    transferred: resolved_transfer(old, &phases[index]),
                })
            },
            cancelled,
        )?;
        match &self.owner.interactions {
            Some(interactions) => checked(
                base + interactions.enthalpy_increment_j(
                    initial.solvent_volume,
                    initial.temperature,
                    &initial.phases,
                    temperature,
                    phases,
                )?,
                "total interacting phase enthalpy increment",
            ),
            None => Ok(base),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn stage(
        &self,
        original: &SolubilityThermalKineticStateV1,
        current: &WorkingState,
        duration: Time,
        heat: StageHeat,
        controls: SolubilityCaloricSolveV1,
        root_iterations: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<StageResult, SolubilityCaloricErrorV1> {
        let baseline = self.change(
            original,
            current.temperature,
            &current.phases[..self.kinetics().len()],
            cancelled,
        )?;
        // The prescribed path keeps its original subtraction and evaluation
        // order. Reservoir heat instead depends on the unknown endpoint.
        let prescribed_required = match heat {
            StageHeat::Prescribed(target) => {
                Some(checked(target - baseline, "thermal kinetic stage heat")?)
            }
            StageHeat::Reservoir { .. } | StageHeat::Finite { .. } => None,
        };
        // Interacting phase evolution freezes the stage's initial temperature.
        // The subsequent heat solve changes only T, so the source-qualified
        // positive frozen-stock heat capacity establishes a monotone root.
        // Native phase trajectories and ordered cross-component updates remain
        // solely owned here; no additional integrator or phase state is created.
        let count = self.kinetics().len();
        let mut frozen = *current;
        let mut transfers = [0.0; MAX_SOLUBILITY_CALORIC_COMPONENTS_V1];
        if let Some(interactions) = &self.owner.interactions {
            for (index, kinetic) in self.kinetics().iter().enumerate() {
                poll(cancelled)?;
                let saturation = interactions.conditional_saturation(
                    current.temperature,
                    original.solvent_volume,
                    &frozen.phases[..count],
                    index,
                )?;
                let step = kinetic.advance_at_saturation(
                    current.temperature,
                    original.solvent_volume,
                    &current.phases[index],
                    duration,
                    saturation,
                )?;
                frozen.phases[index] = *step.state();
                transfers[index] = step.transferred_to_solution().value();
            }
        }
        let mut sample = |t| -> Result<Sample, SolubilityCaloricErrorV1> {
            let mut state = frozen;
            state.temperature = Temperature::new(t);
            let mut change = self.caloric().phase_enthalpy_increment(
                current.temperature,
                original.solvent_moles,
                state.temperature,
                |index| {
                    let old = &current.phases[index];
                    let transferred = if self.owner.interactions.is_some() {
                        transfers[index]
                    } else {
                        let step = self.kinetics()[index].advance(
                            state.temperature,
                            original.solvent_volume,
                            old,
                            duration,
                        )?;
                        state.phases[index] = *step.state();
                        step.transferred_to_solution().value()
                    };
                    Ok(PhaseIncrement {
                        total: old.total_solute().value(),
                        solid: old.solid().value(),
                        dissolved: old.dissolved().value(),
                        transferred,
                    })
                },
                cancelled,
            )?;
            if let Some(interactions) = &self.owner.interactions {
                change = checked(
                    change
                        + interactions.enthalpy_increment_j(
                            original.solvent_volume,
                            current.temperature,
                            &current.phases[..count],
                            state.temperature,
                            &state.phases[..count],
                        )?,
                    "interacting thermal stage enthalpy increment",
                )?;
                poll(cancelled)?;
            }
            Ok(Sample {
                state,
                residual: checked(
                    change
                        - match prescribed_required {
                            Some(required) => required,
                            None => heat.required(baseline, duration, state.temperature)?,
                        },
                    "thermal kinetic stage residual",
                )?,
            })
        };
        let center = sample(current.temperature.value())?;
        if center.residual == 0.0 {
            return heat.finish(center.state, duration);
        }
        let [lo, hi] = controls.temperature_bracket.map(|t| t.value());
        let (mut lower, mut upper) = if center.residual < 0.0 {
            (center, sample(hi)?)
        } else {
            (sample(lo)?, center)
        };
        if lower.residual > 0.0 || upper.residual < 0.0 {
            return Err(SolubilityCaloricErrorV1::HeatOutsideBracket {
                residuals_j: [lower.residual, upper.residual],
            });
        }
        for iteration in 0..=controls.maximum_iterations {
            let lo = lower.state.temperature.value();
            let hi = upper.state.temperature.value();
            let best = if lower.residual.abs() <= upper.residual.abs() {
                &lower
            } else {
                &upper
            };
            if best.residual == 0.0
                || (hi - lo <= controls.temperature_tolerance_k
                    && best.residual.abs() <= controls.energy_tolerance_j)
            {
                return heat.finish(best.state, duration);
            }
            if iteration == controls.maximum_iterations {
                break;
            }
            let midpoint = 0.5 * lo + 0.5 * hi;
            if midpoint == lo || midpoint == hi {
                break;
            }
            *root_iterations += 1;
            let middle = sample(midpoint)?;
            if middle.residual <= 0.0 {
                lower = middle;
            } else {
                upper = middle;
            }
        }
        Err(SolubilityCaloricErrorV1::IterationLimit)
    }
}

/// Explicit scalar-root, temporal-accuracy and work bounds for one window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolubilityThermalKineticControlsV1 {
    /// The bracket applies to every stage. Energy tolerance applies to the
    /// independently checked WHOLE window, not a sum of per-stage tolerances.
    pub root: SolubilityCaloricSolveV1,
    pub initial_step: Time,
    /// Smallest attempted comparison interval; its accepted half stages may
    /// be shorter. A final remainder below this value is attempted once.
    pub minimum_step: Time,
    pub maximum_step: Time,
    /// Relative temporal tolerance in [0, 1], applied to each phase separately.
    pub relative_tolerance: f64,
    pub amount_absolute_tolerance_mol: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Total accepted plus rejected comparison intervals, 1..=65_536.
    pub maximum_attempts: usize,
}

/// Immutable actual phase stocks and common temperature; no equilibrium claim.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityThermalKineticStateV1 {
    owner: Arc<ThermalKineticOwner>,
    temperature: Temperature,
    solvent_moles: Amount,
    solvent_volume: Volume,
    phases: Arc<[SolubilityKineticStateV1]>,
}
impl SolubilityThermalKineticStateV1 {
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn solvent_moles(&self) -> Amount {
        self.solvent_moles
    }
    pub const fn solvent_volume(&self) -> Volume {
        self.solvent_volume
    }
    pub fn phases(&self) -> &[SolubilityKineticStateV1] {
        &self.phases
    }
}

/// Accepted complete window with independent conservation and temporal receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityThermalKineticOutcomeV1 {
    final_state: SolubilityThermalKineticStateV1,
    duration: Time,
    heat_input_j: f64,
    enthalpy_change_j: f64,
    energy_residual_j: f64,
    accepted_steps: usize,
    rejected_steps: usize,
    root_iterations: usize,
    maximum_accepted_error_norm: f64,
}
impl SolubilityThermalKineticOutcomeV1 {
    pub const fn final_state(&self) -> &SolubilityThermalKineticStateV1 {
        &self.final_state
    }
    pub const fn duration(&self) -> Time {
        self.duration
    }
    /// Signed heat delivered to the bath [J]. For uniform heating this is the
    /// prescribed window input; for reservoir exchange it is the independently
    /// accumulated heat of the accepted implicit half stages.
    pub const fn heat_input_j(&self) -> f64 {
        self.heat_input_j
    }
    pub const fn enthalpy_change_j(&self) -> f64 {
        self.enthalpy_change_j
    }
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
    /// Accepted comparison intervals, each integrated with two half stages.
    pub const fn accepted_steps(&self) -> usize {
        self.accepted_steps
    }
    pub const fn rejected_steps(&self) -> usize {
        self.rejected_steps
    }
    /// Total bisections across all accepted and rejected trials.
    pub const fn root_iterations(&self) -> usize {
        self.root_iterations
    }
    pub const fn maximum_accepted_error_norm(&self) -> f64 {
        self.maximum_accepted_error_norm
    }
}

// All trial storage is bounded independently of the requested duration or work.
#[derive(Clone, Copy)]
struct WorkingState {
    temperature: Temperature,
    phases: [SolubilityKineticStateV1; MAX_SOLUBILITY_CALORIC_COMPONENTS_V1],
}
impl WorkingState {
    fn from(state: &SolubilityThermalKineticStateV1) -> Self {
        let mut phases = [state.phases[0]; MAX_SOLUBILITY_CALORIC_COMPONENTS_V1];
        phases[..state.phases.len()].copy_from_slice(&state.phases);
        Self {
            temperature: state.temperature,
            phases,
        }
    }
}
struct Sample {
    state: WorkingState,
    residual: f64,
}

fn half_interval(interval: f64, minimum: f64) -> Option<f64> {
    let half = 0.5 * interval;
    (half >= minimum && half > 0.0 && half < interval).then_some(half)
}

// Choose the phase whose subtraction has the smaller absolute rounding scale.
// In particular a decaying solid tail remains visible after dissolved==total.
fn resolved_transfer(
    initial: &SolubilityKineticStateV1,
    final_state: &SolubilityKineticStateV1,
) -> f64 {
    if initial
        .dissolved()
        .value()
        .max(final_state.dissolved().value())
        <= initial.solid().value().max(final_state.solid().value())
    {
        final_state.dissolved().value() - initial.dissolved().value()
    } else {
        initial.solid().value() - final_state.solid().value()
    }
}

fn temporal_error(
    full: &WorkingState,
    fine: &WorkingState,
    count: usize,
    controls: SolubilityThermalKineticControlsV1,
) -> f64 {
    let mut error = normalized_difference(
        full.temperature.value(),
        fine.temperature.value(),
        controls.temperature_absolute_tolerance_k,
        controls.relative_tolerance,
    );
    for (a, b) in full.phases[..count].iter().zip(&fine.phases[..count]) {
        for (first, second) in [(a.solid(), b.solid()), (a.dissolved(), b.dissolved())] {
            error = error.max(normalized_difference(
                first.value(),
                second.value(),
                controls.amount_absolute_tolerance_mol,
                controls.relative_tolerance,
            ));
        }
    }
    error
}

fn normalized_difference(first: f64, second: f64, absolute: f64, relative: f64) -> f64 {
    if first == second {
        return 0.0;
    }
    let magnitude = first.max(second);
    let scale = magnitude.max(absolute);
    // Scale before adding tolerances so a large valid inventory cannot overflow
    // the denominator and silently report zero temporal disagreement.
    ((first - second).abs() / scale) / (absolute / scale + relative * (magnitude / scale))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod interaction_tests;
