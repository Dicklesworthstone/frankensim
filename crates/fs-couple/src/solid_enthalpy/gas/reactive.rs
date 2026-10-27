//! Declared gas chemistry heats the original spatial solid and material fields.
//!
//! The existing reaction owner first advances composition at fixed volume and
//! TOTAL internal energy, including formation energy. The existing joint
//! gas/cell enthalpy owner then exchanges heat over the SAME physical interval.
//! This is a first-order reaction -> thermal split, not a new reaction or heat
//! solver. Continuing the accepted gas makes later rates see its actual cooled
//! temperature, composition and pressure. Time refinement remains necessary.

/// Shared implicit endpoint chemistry and spatial gas/solid heat exchange.
pub mod coupled;

use super::{
    Cx, DepositionSimulation, IsochoricGasParcel, PreparedMaterialNetwork, SolidEnthalpy,
    SolidEnthalpyState, SolidGasContact, SolidGasError, SolidGasReport, SolidGasSolve,
};
use crate::cavity_phs::piston::reactive::{
    GasReaction, GasReactionError, GasReactionSolve, GasReactionStep,
};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// Original chemical and thermal refusals remain distinguishable.
#[derive(Debug)]
pub enum ReactiveSolidGasError {
    Reaction(GasReactionError),
    Thermal(SolidGasError),
    Invalid(&'static str),
    Cancelled,
}
impl std::fmt::Display for ReactiveSolidGasError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "reacting gas/spatial solid coupling refused: {self:?}")
    }
}
impl std::error::Error for ReactiveSolidGasError {}
impl From<GasReactionError> for ReactiveSolidGasError {
    fn from(error: GasReactionError) -> Self {
        match error {
            GasReactionError::Cancelled => Self::Cancelled,
            other => Self::Reaction(other),
        }
    }
}
impl From<SolidGasError> for ReactiveSolidGasError {
    fn from(error: SolidGasError) -> Self {
        match error {
            SolidGasError::Cancelled => Self::Cancelled,
            other => Self::Thermal(other),
        }
    }
}

/// Stage receipts and an independently measured account across BOTH stages.
#[derive(Clone, Debug)]
pub struct ReactiveSolidGasReport {
    /// One interval shared by chemistry and thermal exchange, not their sum.
    pub duration_s: f64,
    /// The original isolated chemical stage, including its own entropy gate.
    /// Formation energy stays inside gas storage and is never credited as heat.
    pub reaction: GasReactionStep,
    /// Original spatial heat solve; its gas account starts AFTER chemistry.
    pub thermal: SolidGasReport,
    /// Actual m*(u_final-u_before_reaction) [J], on the same chemical datum.
    pub gas_energy_change_j: f64,
    /// Gas change + cell enthalpy change - authored external heat [J].
    /// Internal patch transfers cancel; reaction energy is not an extra source.
    pub combined_energy_residual_j: f64,
}

/// No stage mutates its inputs. All original endpoints publish together.
pub struct ReactiveSolidGasCandidate<'m, 't> {
    thermal: DepositionSimulation<'t>,
    gas: IsochoricGasParcel,
    state: SolidEnthalpyState<'m>,
    report: ReactiveSolidGasReport,
}
impl<'m, 't> ReactiveSolidGasCandidate<'m, 't> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    pub const fn state(&self) -> &SolidEnthalpyState<'m> {
        &self.state
    }
    pub const fn report(&self) -> &ReactiveSolidGasReport {
        &self.report
    }
    /// The returned solid state feeds the existing load, deformation, stress,
    /// prestress and transient-solid owners; no independently animated shape.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        IsochoricGasParcel,
        SolidEnthalpyState<'m>,
        ReactiveSolidGasReport,
    ) {
        (self.thermal, self.gas, self.state, self.report)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ReactiveSolidGasError> {
    cx.checkpoint()
        .map_err(|_| ReactiveSolidGasError::Cancelled)
}

impl<'m> SolidEnthalpy<'m> {
    /// React one finite gas and exchange its heat with the SAME spatial cells.
    ///
    /// `reaction` retains the existing declared mechanism, source rate scales,
    /// exact species/element conventions and optional implicit/subcycling policy.
    /// The accepted gas is cloned into that owner without a composition reset.
    /// Products, NASA-9 caloric data and material transport then enter the original
    /// `step_gas_with_transport` solve at every actual thermal endpoint.
    ///
    /// `transport` selects the original prepared material network, when needed,
    /// so internal k(T) and gas-facing material layers share pinned sources.
    /// None retains authored internal conductances. Patch geometry, local cell
    /// temperatures and complete thermal/source-domain checks remain unchanged.
    /// `external_gas_heat_j` is external heat ONLY, never chemical heat release.
    /// The total tolerance encloses both stages; the temperature-increment limit
    /// also bounds the chemical intermediate and final gas from the original T.
    ///
    /// The split is first order. Gas stays well mixed at fixed mass and volume;
    /// there is no flame/shock, gas traction, moving boundary or melting model.
    /// The returned solid field retains the existing sequential thermoelastic
    /// interpretation: its caloric audit does not include elastic energy or
    /// reciprocal thermoelastic heat. Reaction entropy belongs to the isolated
    /// chemical stage, not a claimed whole-system entropy-production theorem.
    ///
    /// Any reaction, thermal, source, work, interval or cancellation refusal
    /// publishes neither stage. Process/split thermal events as required by the
    /// original backward-Euler owner before requesting the interval.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_reacting_isochoric_gas<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        gas: &IsochoricGasParcel,
        reaction: &GasReaction,
        target_s: f64,
        external_gas_heat_j: f64,
        contacts: &[SolidGasContact],
        transport: Option<&PreparedMaterialNetwork<'_, '_>>,
        reaction_controls: GasReactionSolve,
        controls: SolidGasSolve,
        initial_energy_increments_j: &[f64],
        cx: &Cx<'_>,
    ) -> Result<ReactiveSolidGasCandidate<'m, 't>, ReactiveSolidGasError> {
        checkpoint(cx)?;
        let duration_s = target_s - thermal.time_s();
        if !target_s.is_finite()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !external_gas_heat_j.is_finite()
            || thermal.cells().len() > controls.maximum_cells
            || contacts.is_empty()
            || contacts.len() > controls.maximum_contacts
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || initial_energy_increments_j.len() != thermal.cells().len().saturating_add(1)
            || initial_energy_increments_j
                .iter()
                .any(|value| !value.is_finite())
            || !controls.maximum_gas_temperature_increment_k.is_finite()
            || controls.maximum_gas_temperature_increment_k <= 0.0
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j < 0.0
        {
            return Err(ReactiveSolidGasError::Invalid(
                "positive represented interval and valid coupled temperature/energy bounds required",
            ));
        }
        // Let the original timeline reject pending events and exhausted thermal
        // work before spending a chemical stage. The heat owner repeats this
        // admission when it constructs its own unchanged residual.
        fs_thermal::enthalpy::backward_euler::BackwardEulerStep::new(thermal, target_s, || {
            cx.checkpoint()
                .map_err(|_| fs_thermal::enthalpy::DepositionError::Cancelled)
        })
        .map_err(SolidGasError::from)?;
        let before = ClosedGasParcel::from_isochoric(gas.clone());
        let chemistry =
            reaction.parcel_candidate(&before, duration_s, reaction_controls, || {
                cx.checkpoint().is_err()
            })?;
        let initial_temperature = gas.state().state().temperature;
        if (chemistry.gas().state().state().temperature - initial_temperature).abs()
            > controls.maximum_gas_temperature_increment_k
        {
            return Err(ReactiveSolidGasError::Invalid(
                "chemical intermediate exceeds the coupled gas temperature-increment allowance",
            ));
        }
        let (reacted, reaction) = chemistry.into_parts();
        checkpoint(cx)?;
        let candidate = self.step_gas_with_transport(
            thermal,
            &reacted.into_isochoric(),
            target_s,
            external_gas_heat_j,
            contacts,
            transport,
            controls,
            initial_energy_increments_j,
            cx,
        )?;
        let gas_energy_change_j = gas.mass_kg()
            * (candidate.gas().state().specific_internal_energy_j_kg()
                - gas.state().specific_internal_energy_j_kg());
        let heat = candidate.report();
        let total = gas_energy_change_j + heat.thermal.stored_energy_change_j
            - heat.external_thermal_heat_j
            - heat.external_gas_heat_j;
        if !gas_energy_change_j.is_finite()
            || !total.is_finite()
            || total.abs() > controls.total_energy_tolerance_j
            || candidate.gas().mass_kg().to_bits() != gas.mass_kg().to_bits()
            || candidate.gas().volume_m3().to_bits() != gas.volume_m3().to_bits()
            || (candidate.gas().state().state().temperature - initial_temperature).abs()
                > controls.maximum_gas_temperature_increment_k
        {
            return Err(ReactiveSolidGasError::Invalid(
                "reaction and spatial heat do not close the original gas/solid inventories",
            ));
        }
        checkpoint(cx)?;
        let (thermal, gas, state, thermal_report) = candidate.into_parts();
        Ok(ReactiveSolidGasCandidate {
            thermal,
            gas,
            state,
            report: ReactiveSolidGasReport {
                duration_s,
                reaction,
                thermal: thermal_report,
                gas_energy_change_j,
                combined_energy_residual_j: total,
            },
        })
    }
}
