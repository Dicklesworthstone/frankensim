//! Complete heat transport on the original squirrel-cage induction machine.
//!
//! Thermal contacts, layers, radiation and finite coolant laws resolve on the
//! same endpoint as both winding resistances. The native flux/shaft evolution,
//! prepared workspace and three physical heat receivers remain unchanged.

use super::*;
use crate::circuit_enthalpy::contacts::heat::{
    HeatAssembly, HeatNetworkError, closed_balance,
};
use fs_exec::Cx;

/// Preserve native induction failures and the complete assembly's typed refusal.
#[derive(Debug)]
pub enum InductionHeatError<E> {
    Induction(InductionThermalError),
    Heat(E),
}
impl<E> From<InductionThermalError> for InductionHeatError<E> {
    fn from(error: InductionThermalError) -> Self {
        Self::Induction(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for InductionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "induction heat assembly: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for InductionHeatError<E> {}

type Error<E> = CircuitEnthalpyError<InductionHeatError<E>>;

fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(InductionHeatError::Heat(other.into())),
    }
}

/// Two source-bound windings and explicit bearing heat on one complete assembly.
///
/// Cells and original caloric/material bindings come from this assembly alone.
/// Stator and rotor cells remain distinct. The bearing receiver may coincide
/// with either winding or be a separate bound cell. Graph resistors, electrical
/// contacts and circuit-basis bindings refuse because this reduced machine
/// does not own their work. Keep the same assembly across accepted steps to
/// retain its thermal source selections.
pub struct InductionHeatNetwork<'a, A: HeatAssembly + ?Sized> {
    pub assembly: &'a A,
    pub stator: InductionWinding,
    pub rotor: InductionWinding,
    pub bearing_cell: CellId,
}

impl<A: HeatAssembly + ?Sized> InductionHeatNetwork<'_, A> {
    fn network(&self, cx: &Cx<'_>) -> Result<InductionThermalNetwork<'_>, Error<A::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = self.assembly.prepared_heat_network();
        let bulk = prepared.thermal_material_network().map_err(heat_error)?;
        prepared.admit_field_independent_boundaries(cx, &[]).map_err(heat_error)?;
        Ok(InductionThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            stator: self.stator.clone(),
            rotor: self.rotor.clone(),
            bearing_cell: self.bearing_cell,
        })
    }

    fn resolve(
        &self,
        cx: &Cx<'_>,
        network: &InductionThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        original: Option<&[Arc<ResolvedConductor>; 2]>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(CircuitThermalCoefficients<A::Sample>, [ResolvedConductor; 2]), Error<A::Error>> {
        let coefficients = self.assembly.resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(InductionHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(heat_error(HeatNetworkError::Invalid(
                "induction heat assembly cannot discard graph resistor or electrical-contact work",
            )));
        }
        let windings = network.resolve_windings(
            thermal, &coefficients.bulk.source.as_ref().material, original, cancel,
        ).map_err(induction_error)?;
        Ok((coefficients.bulk, windings))
    }

    /// Bind both original conductor sources without changing fluxes, angle,
    /// momenta, finite enthalpy or time. Numeric winding resistance placeholders
    /// are replaced in the explicitly authored native winding basis.
    pub fn bind(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        state: [f64; 7],
        parameters: InductionMotorParameters,
    ) -> Result<MaterialInductionMotor, Error<A::Error>> {
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.bind_resolved(
            thermal,
            state,
            parameters,
            || self.resolve(cx, &network, thermal, None, &cancel).map(|(_, windings)| windings),
            &cancel,
        )
    }

    /// Jointly solve native induction, both winding resistances and the entire
    /// heat assembly at one endpoint. Inputs are held
    /// `[v_alpha, v_beta, tau_rotor, tau_stator]`; the shorted cage remains an
    /// original flux state. Native stator/rotor copper and bearing heat enter
    /// their receivers once, while additional transport keeps its original
    /// source payload in `report.source.network`.
    ///
    /// Closed-radiation roundoff is internal: external heat is
    /// `report.non_electrical_heat_j` minus duration times
    /// `network.as_ref().radiation_net_power_w()`. The combined balance is
    /// independently admitted after that correction without changing a state
    /// or depositing another source. Cancellation or any source, work-budget
    /// or solve refusal publishes neither machine nor thermal history.
    ///
    /// Magnetic inductances, geometry, inertia and bearing drag remain fixed.
    /// This preserves the native reduced machine's generally first-order
    /// thermal feedback; it adds no saturation, skin effect or moving geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialInductionMotor,
        inputs: &[f64; 4],
        controls: InductionEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<InductionEnthalpyCandidate<'t, A::Sample>, Error<A::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        admit_driver_unknowns(thermal, controls.thermal, initial_energy_increments_j)?;
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let candidate = network.step_resolved_candidate(
            thermal,
            target_s,
            motor,
            inputs,
            controls,
            initial_energy_increments_j,
            |state, original| self.resolve(cx, &network, state, original, &cancel),
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(), candidate.report().total_energy_residual_j,
            candidate.report().source.network.as_ref(), controls.thermal.total_energy_tolerance_j,
        ).map_err(heat_error)?;
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }
}
