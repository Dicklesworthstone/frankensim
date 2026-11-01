//! Adiabatic equilibrium compressor/expander outlets over the original mixture.
//!
//! A reversible pressure/entropy solve establishes the reference work. An
//! authored isentropic efficiency then selects an actual pressure/enthalpy
//! endpoint. Both use the same source and phase policy; phase amounts, latent
//! energy and species reference datums never become a fitted heat capacity.

use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError,
    equilibrium::{
        PhaseEnthalpyControls, ResolvedThermochemicalPhaseState, entropy::PhaseEntropyControls,
    },
};
use fs_qty::{MassFlowRate, Power, Pressure, Temperature};

/// Whole-machine ceiling: at most two bounded original caloric inversions.
pub const MAX_MACHINE_PT_EVALUATIONS: usize = 256;

/// Pressure direction and convention for the supplied isentropic efficiency.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixtureMachineKind {
    /// Pressure rises; actual shaft work into the fluid is ideal work / efficiency.
    Compressor,
    /// Pressure falls; actual work out of the fluid is efficiency * ideal work.
    Expander,
}

/// Source domains, physical residuals and a shared bound on all attempted work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureMachineControl {
    /// Reversible outlet at the original inlet entropy.
    pub isentropic: PhaseEntropyControls,
    /// Actual outlet at the efficiency-selected enthalpy increment.
    /// Its phase policy must equal `isentropic.phase_options()` exactly.
    pub outlet: PhaseEnthalpyControls,
    /// Allowed numerical shortfall from nonnegative actual Delta s [J/(kg K)].
    /// The actual value is reported; a negative residual is never clamped away.
    pub entropy_tolerance_j_kg_k: f64,
    /// Total original PT calls in both stages, including both reference queries.
    pub maximum_pt_evaluations: usize,
}
impl MixtureMachineControl {
    /// Start with the same explicit temperature domain for both original solves.
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2]) -> Self {
        Self {
            isentropic: PhaseEntropyControls::new(temperature_bracket),
            outlet: PhaseEnthalpyControls::new(temperature_bracket),
            entropy_tolerance_j_kg_k: 1e-7,
            maximum_pt_evaluations: MAX_MACHINE_PT_EVALUATIONS,
        }
    }

    /// Admit the complete numerical and phase policy before source work begins.
    pub fn validate(self) -> Result<(), MixtureMachineError> {
        self.isentropic
            .validate()
            .map_err(|error| phase_error(MixtureMachineStage::Isentropic, 0, error))?;
        self.outlet
            .validate()
            .map_err(|error| phase_error(MixtureMachineStage::ActualOutlet, 0, error))?;
        if self.isentropic.phase_options() != self.outlet.phase_options() {
            return Err(MixtureMachineError::Invalid(
                "machine stages require identical original phase policies",
            ));
        }
        if !positive(self.entropy_tolerance_j_kg_k)
            || self.maximum_pt_evaluations == 0
            || self.maximum_pt_evaluations > MAX_MACHINE_PT_EVALUATIONS
        {
            return Err(MixtureMachineError::Invalid(
                "positive entropy tolerance and bounded total PT work required",
            ));
        }
        Ok(())
    }
}

/// Original caloric stage responsible for a source or numerical refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixtureMachineStage {
    /// Constant-inlet-entropy pressure change.
    Isentropic,
    /// Efficiency-selected enthalpy change at the final pressure.
    ActualOutlet,
}

/// No partially resolved machine outlet is returned on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum MixtureMachineError {
    /// Invalid physical input, conflicting phase policy or unrepresentable work.
    Invalid(&'static str),
    /// The original phase/caloric owner refused a state or bounded solve.
    Equilibrium {
        /// Stage containing the refusal.
        stage: MixtureMachineStage,
        /// Successful PT evaluations completed before this stage began.
        completed_pt_evaluations: usize,
        /// Original refusal, including this stage's attempted-work count.
        source: ThermochemicalFlashError,
    },
    /// The shared whole-machine PT allowance was exhausted between stages.
    WorkLimit {
        /// Original PT calls completed before refusing another stage.
        pt_evaluations: usize,
    },
    /// Actual source entropy change violates the declared adiabatic allowance.
    Entropy {
        /// Actual integrated source Delta s [J/(kg K)].
        change_j_kg_k: f64,
        /// Admitted numerical tolerance [J/(kg K)].
        tolerance_j_kg_k: f64,
    },
    /// Cancellation before complete publication.
    Cancelled,
}
impl core::fmt::Display for MixtureMachineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "equilibrium mixture machine refusal: {self:?}")
    }
}
impl std::error::Error for MixtureMachineError {}

/// Accepted steady, adiabatic outlet with source-integrated work and entropy.
#[derive(Clone, Debug, PartialEq)]
pub struct MixtureMachineStep {
    isentropic_outlet: ResolvedThermochemicalPhaseState,
    outlet: ResolvedThermochemicalPhaseState,
    /// Declared machine and signed pressure direction.
    pub kind: MixtureMachineKind,
    /// Authored efficiency, strictly positive and at most one.
    pub efficiency: f64,
    /// Original-source reversible Delta h [J/kg].
    pub isentropic_enthalpy_change_j_kg: f64,
    /// Actual original-source Delta h [J/kg], positive for work into the fluid.
    /// Stream enthalpy already includes pressure-flow work; do not add it again.
    pub shaft_work_to_fluid_j_kg: f64,
    /// Actual original-source Delta s [J/(kg K)], without clipping roundoff.
    pub entropy_production_j_kg_k: f64,
    /// Actual Delta h minus the efficiency-selected target [J/kg].
    pub enthalpy_residual_j_kg: f64,
    /// Actual reversible-reference Delta s minus zero [J/(kg K)].
    pub isentropic_entropy_residual_j_kg_k: f64,
    /// Total original PT calls in all stages, including repeated references.
    pub pt_evaluations: usize,
    /// Original observed screen work plus optional split-work ceilings.
    /// This is an upper bound, not a measured EOS evaluation count.
    pub eos_evaluation_upper_bound: usize,
}
impl MixtureMachineStep {
    /// Reference outlet at the original inlet entropy and declared final pressure.
    #[must_use]
    pub const fn isentropic_outlet(&self) -> &ResolvedThermochemicalPhaseState {
        &self.isentropic_outlet
    }
    /// Actual outlet with its original source and finite phase-screen evidence.
    #[must_use]
    pub const fn outlet(&self) -> &ResolvedThermochemicalPhaseState {
        &self.outlet
    }
    /// Move the accepted state without reevaluating its source or phase split.
    #[must_use]
    pub fn into_outlet(self) -> ResolvedThermochemicalPhaseState {
        self.outlet
    }
    /// Signed shaft power for a separately supplied nonnegative mass flow.
    /// This conversion does not predict flow from geometry or pressure ratio.
    pub fn shaft_power_to_fluid(&self, flow: MassFlowRate) -> Result<Power, MixtureMachineError> {
        let mass_flow = flow.value();
        let power = mass_flow * self.shaft_work_to_fluid_j_kg;
        if !mass_flow.is_finite() || mass_flow < 0.0 || !power.is_finite() {
            return Err(MixtureMachineError::Invalid(
                "nonnegative finite mass flow and representable shaft power required",
            ));
        }
        Ok(Power::new(power))
    }
}

/// Resolve a compressor or expander using one immutable feed and original EOS.
///
/// The steady-flow energy balance assumes no heat transfer and negligible
/// inlet/outlet kinetic and potential energy differences: shaft work into the
/// fluid equals its integrated enthalpy rise. The efficiency is supplied by the
/// caller; this API provides no empirical machine map, shaft dynamics, phase
/// kinetics, species reactions, flow-rate law or global phase-stability proof.
///
/// Efficiency one reuses the accepted isentropic endpoint without a PH solve.
/// Equal pressure is exact zero work after one original inlet PT selection,
/// provided the reference lies inside the declared isentropic bracket.
/// All source states remain unpublished until the final cancellation poll.
#[allow(clippy::too_many_arguments)]
pub fn adiabatic_mixture_machine(
    model: &ThermochemicalPengRobinsonModel,
    kind: MixtureMachineKind,
    inlet_temperature: Temperature,
    inlet_pressure: Pressure,
    outlet_pressure: Pressure,
    efficiency: f64,
    controls: MixtureMachineControl,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixtureMachineStep, MixtureMachineError> {
    poll(&mut cancelled)?;
    controls.validate()?;
    let inlet_p = inlet_pressure.value();
    let outlet_p = outlet_pressure.value();
    if !positive(inlet_temperature.value())
        || !positive(inlet_p)
        || !positive(outlet_p)
        || !positive(efficiency)
        || efficiency > 1.0
        || matches!(kind, MixtureMachineKind::Compressor) && outlet_p < inlet_p
        || matches!(kind, MixtureMachineKind::Expander) && outlet_p > inlet_p
    {
        return Err(MixtureMachineError::Invalid(
            "positive inlet state, consistent pressure direction and efficiency in (0,1] required",
        ));
    }
    let mut isentropic_control = controls.isentropic;
    isentropic_control.maximum_pt_evaluations = isentropic_control
        .maximum_pt_evaluations
        .min(controls.maximum_pt_evaluations);
    let reversible = model
        .phase_equilibrium_pressure_entropy_increment(
            inlet_temperature,
            inlet_pressure,
            outlet_pressure,
            0.0,
            isentropic_control,
            &mut cancelled,
        )
        .map_err(|error| phase_error(MixtureMachineStage::Isentropic, 0, error))?;
    let ideal_h = reversible
        .enthalpy_change_j_kg()
        .ok_or(MixtureMachineError::Invalid(
            "isentropic source did not retain its enthalpy increment",
        ))?;
    admit_work_direction(kind, inlet_p == outlet_p, ideal_h)?;
    let requested_h = match kind {
        MixtureMachineKind::Compressor => ideal_h / efficiency,
        MixtureMachineKind::Expander => ideal_h * efficiency,
    };
    if !requested_h.is_finite() || ideal_h != 0.0 && requested_h == 0.0 {
        return Err(MixtureMachineError::Invalid(
            "efficiency-selected enthalpy increment is unrepresentable",
        ));
    }
    let mut pt_evaluations = reversible.pt_evaluations();
    let mut eos_evaluation_upper_bound = reversible.eos_evaluation_upper_bound();
    let entropy_residual = reversible.residual_j_kg_k();
    let (outlet, actual_h, actual_s, enthalpy_residual) =
        if efficiency == 1.0 || inlet_p == outlet_p {
            (
                reversible.state().clone(),
                ideal_h,
                reversible
                    .accepted_increment_j_kg_k()
                    .ok_or(MixtureMachineError::Invalid(
                        "isentropic source did not retain its entropy increment",
                    ))?,
                0.0,
            )
        } else {
            let remaining = controls.maximum_pt_evaluations - pt_evaluations;
            if remaining == 0 {
                return Err(MixtureMachineError::WorkLimit { pt_evaluations });
            }
            let mut actual_control = controls.outlet;
            actual_control.enthalpy.maximum_flash_evaluations = actual_control
                .enthalpy
                .maximum_flash_evaluations
                .min(remaining);
            // Resolve small nonzero work through original source increments,
            // rather than accepting it within a larger absolute heat tolerance.
            actual_control.enthalpy.enthalpy_tolerance_j_kg = actual_control
                .enthalpy
                .enthalpy_tolerance_j_kg
                .min(0.25 * requested_h.abs());
            let actual = model
                .phase_equilibrium_pressure_enthalpy_increment(
                    inlet_temperature,
                    inlet_pressure,
                    outlet_pressure,
                    requested_h,
                    actual_control,
                    &mut cancelled,
                )
                .map_err(|error| {
                    phase_error(MixtureMachineStage::ActualOutlet, pt_evaluations, error)
                })?;
            pt_evaluations += actual.pt_evaluations();
            eos_evaluation_upper_bound += actual.eos_evaluation_upper_bound();
            let actual_h = actual
                .accepted_increment_j_kg()
                .ok_or(MixtureMachineError::Invalid(
                    "actual source did not retain its enthalpy increment",
                ))?;
            let actual_s = actual
                .entropy_change_j_kg_k()
                .ok_or(MixtureMachineError::Invalid(
                    "actual source did not retain its entropy increment",
                ))?;
            let residual = actual.residual_j_kg();
            (actual.into_state(), actual_h, actual_s, residual)
        };
    admit_work_direction(kind, inlet_p == outlet_p, actual_h)?;
    if !actual_s.is_finite() || actual_s < -controls.entropy_tolerance_j_kg_k {
        return Err(MixtureMachineError::Entropy {
            change_j_kg_k: actual_s,
            tolerance_j_kg_k: controls.entropy_tolerance_j_kg_k,
        });
    }
    poll(&mut cancelled)?;
    Ok(MixtureMachineStep {
        isentropic_outlet: reversible.into_state(),
        outlet,
        kind,
        efficiency,
        isentropic_enthalpy_change_j_kg: ideal_h,
        shaft_work_to_fluid_j_kg: actual_h,
        entropy_production_j_kg_k: actual_s,
        enthalpy_residual_j_kg: enthalpy_residual,
        isentropic_entropy_residual_j_kg_k: entropy_residual,
        pt_evaluations,
        eos_evaluation_upper_bound,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MixtureMachineError> {
    if cancelled() {
        Err(MixtureMachineError::Cancelled)
    } else {
        Ok(())
    }
}
fn phase_error(
    stage: MixtureMachineStage,
    completed_pt_evaluations: usize,
    source: ThermochemicalFlashError,
) -> MixtureMachineError {
    if source == ThermochemicalFlashError::Cancelled {
        MixtureMachineError::Cancelled
    } else {
        MixtureMachineError::Equilibrium {
            stage,
            completed_pt_evaluations,
            source,
        }
    }
}
fn admit_work_direction(
    kind: MixtureMachineKind,
    equal_pressure: bool,
    change_j_kg: f64,
) -> Result<(), MixtureMachineError> {
    let valid = change_j_kg.is_finite()
        && if equal_pressure {
            change_j_kg == 0.0
        } else {
            match kind {
                MixtureMachineKind::Compressor => change_j_kg > 0.0,
                MixtureMachineKind::Expander => change_j_kg < 0.0,
            }
        };
    if valid {
        Ok(())
    } else {
        Err(MixtureMachineError::Invalid(
            "original enthalpy change does not resolve the declared work direction",
        ))
    }
}
