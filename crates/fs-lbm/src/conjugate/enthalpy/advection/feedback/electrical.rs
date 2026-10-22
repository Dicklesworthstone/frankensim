//! The exact live-flow material epoch supplies current and irreversible heat.
//!
//! Momentum, advection and fluid-face conduction remain in MaterialFlowFeedback.
//! The original Ohmic Account owns current, contact deposition and electrical
//! work. No second material evaluation or thermal stepping loop is introduced.

use super::*;
use crate::conjugate::enthalpy::{
    Account, EnthalpyCurrentState, JouleEnthalpyError, JouleEnthalpyInput,
    JouleEnthalpyReport, StepConfig, admit_map as admit_electrical_map, checkpoint,
};

/// Original flow/history and electrical failures remain distinguishable.
#[derive(Debug)]
pub enum MaterialFlowElectricalError<E> {
    Flow(MaterialFlowError<E>),
    Electrical(JouleEnthalpyError),
}
impl<E> From<DepositionError> for MaterialFlowElectricalError<E> {
    fn from(error: DepositionError) -> Self { Self::Electrical(error.into()) }
}
impl<E> From<JouleEnthalpyError> for MaterialFlowElectricalError<E> {
    fn from(error: JouleEnthalpyError) -> Self { Self::Electrical(error) }
}
impl<E> From<MaterialFlowError<E>> for MaterialFlowElectricalError<E> {
    fn from(error: MaterialFlowError<E>) -> Self { Self::Flow(error) }
}
impl<E: core::fmt::Display> core::fmt::Display for MaterialFlowElectricalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(error) => write!(f, "electrical material flow: {error}"),
            Self::Electrical(error) => write!(f, "electrical material flow: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MaterialFlowElectricalError<E> {}

type ElectricalResponse<'a, F> =
    dyn Fn(&FlowMaterialSample<F>) -> Result<Option<f64>, ChtError> + 'a;

#[derive(Clone)]
struct CurrentFlow<'a, 'f, F, R, O> {
    flow: MaterialFlowFeedback<'f, F, R, O>,
    account: Account<'a>,
    response: &'a ElectricalResponse<'a, F>,
}

// Require the complete actual epoch, including ideal electrode voxels. The
// getter below reads a retained material object; it never re-resolves its T/h.
fn evaluate<F>(
    account: &mut Account<'_>,
    thermal: &DepositionSimulation<'_>,
    samples: &[FlowMaterialSample<F>],
    response: &ElectricalResponse<'_, F>,
) -> Result<EnthalpyCurrentState, JouleEnthalpyError> {
    admit_electrical_map(thermal, account.input, account.gate)?;
    if samples.len() != account.input.cell_ids.len() {
        return Err(JouleEnthalpyError::Invalid("current needs every original flow material sample"));
    }
    for (voxel, (sample, &id)) in samples.iter().zip(account.input.cell_ids).enumerate() {
        checkpoint(account.gate)?;
        if sample.voxel != voxel || sample.cell != id
            || sample.specific_enthalpy_j_kg.to_bits() != thermal.cells()[&id].specific_h_j_kg.to_bits()
            || sample.temperature_k.to_bits() != thermal.temperature_k(id)?.to_bits()
        {
            return Err(JouleEnthalpyError::Invalid(
                "electrical response must use the exact receiving flow material epoch",
            ));
        }
    }
    account.evaluate(thermal, &|voxel, _, _| response(&samples[voxel]))
}

impl<F: FluidTransport + Clone, R, O: DepositionObserver> DepositionObserver
    for CurrentFlow<'_, '_, F, R, O>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    type Error = MaterialFlowElectricalError<O::Error>;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.flow.begin()?;
        self.account = Account::new(self.account.input, self.account.gate);
        Ok(())
    }

    fn sources(
        &mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        check()?;
        // Refuse before another expensive momentum solve as well as before
        // current assembly. The same ceiling also reserves no free endpoint.
        if self.account.solves >= self.account.input.max_current_solves {
            return Err(JouleEnthalpyError::CurrentBudget.into());
        }
        let mut sources = self.flow.sources(thermal, check)?;
        let epoch = self.flow.applied().ok_or(JouleEnthalpyError::Invalid(
            "current has no applied flow material epoch",
        ))?;
        if epoch.time_s.to_bits() != thermal.time_s().to_bits() {
            return Err(JouleEnthalpyError::Invalid("stale flow source time").into());
        }
        let state = evaluate(&mut self.account, thermal, &epoch.samples, self.response)?;
        for (cell, watts) in self.account.apply(state)? {
            check()?;
            let value = sources.entry(cell).or_insert(0.0);
            *value += watts;
            if !value.is_finite() { return Err(DepositionError::Arithmetic.into()); }
        }
        Ok(sources)
    }

    fn conductances(
        &mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.flow.conductances(thermal, check).map_err(Into::into)
    }

    fn heat_transfers(
        &mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.flow.heat_transfers(thermal, check).map_err(Into::into)
    }

    fn limit_timestep_s(
        &mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.flow.limit_timestep_s(thermal, proposed_s, check).map_err(Into::into)
    }

    fn event(
        &mut self, event: &Event, thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.flow.event(event, thermal)?;
        admit_electrical_map(thermal, self.account.input, self.account.gate)?;
        Ok(())
    }

    fn step(
        &mut self, start_s: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.flow.step(start_s, dt, old, thermal, check)?;
        self.account.record(start_s, dt)?;
        Ok(())
    }
}

/// Advance current, local material flow and total enthalpy from ONE material epoch.
///
/// Each native thermal substep first resolves the flow's actual h/T/material
/// objects. Their viscosity and any declared buoyancy/mushy mobility drive the
/// original momentum solve; its face mass fluxes advect the original total h.
/// `resistivity` reads those SAME immutable samples for the original current
/// solve. It must return a retained response, not query another material or
/// mutate external physical history. `None` declares insulation; missing source
/// support must return an error. Ideal electrodes retain the existing current
/// owner's convention: their bulk electrical getter is unused.
///
/// Attach `MaterialFlowFeedback::with_fluid_conduction` before this call to use
/// those samples' k for its owned faces as well. Existing inner source/contact
/// history and admitted molecular/Darcy heat continue on the same actual CFL
/// steps. Do not wrap an independent owner that deposits this same Joule loss.
/// Fixed electrical contact laws retain their original area and heat partition.
///
/// The current and flow domains must be the SAME borrowed object and every
/// voxel must be a mapped fluid cell, including ideal electrical electrodes.
/// A potential constraint does not turn that cell into a solid obstacle. The
/// original caloric object, fixed masses and control volumes remain owners.
/// The full flow epoch is checked against the actual receiving cells before
/// each current solve. Endpoint current uses the flow's already resolved
/// endpoint, without another material query or an unreported momentum solve.
///
/// Returns thermal state, flow/inner history and the original electrical report
/// together. Electrode work, contact/bulk heat and mapped deposition use their
/// actual accepted durations; Joule heat is ALREADY in the thermal audit.
/// Current budgets include one solve per substep and one final support solve.
/// Any map, material, current, flow, heat, budget, endpoint or cancellation
/// refusal leaves both inputs untouched. Mutable data captured by callbacks
/// remains caller-owned and must obey the original observer cloning contract.
///
/// This is explicit quasisteady fixed-volume Ohmic/enthalpy coupling. Thermal
/// refinement and rapid momentum equilibration still require justification.
/// No magnetic field/motional EMF, induction, charge convection, density-change
/// volume source, free surface or total kinetic/thermal energy claim is added.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_joule_enthalpy_with_material_flow<'t, 'f, F, R, O>(
    thermal: &DepositionSimulation<'t>, target_s: f64, config: StepConfig,
    input: JouleEnthalpyInput<'_>, flow: &MaterialFlowFeedback<'f, F, R, O>,
    resistivity: impl Fn(&FlowMaterialSample<F>) -> Result<Option<f64>, ChtError>,
) -> Result<
    (DepositionSimulation<'t>, MaterialFlowFeedback<'f, F, R, O>, JouleEnthalpyReport),
    MaterialFlowElectricalError<O::Error>,
>
where
    F: FluidTransport + Clone,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    let gate = flow.gate;
    checkpoint(gate)?;
    admit_electrical_map(thermal, input, gate)?;
    flow.check_map(thermal, &mut || Ok(()))?;
    if !core::ptr::eq(input.domain, flow.domain)
        || flow.domain.fluid_count() != flow.domain.cell_count()
        || input.cell_ids.len() != flow.input.cell_ids.len()
        || input.cell_ids.iter().zip(flow.input.cell_ids)
            .any(|(&electric, &fluid)| fluid != Some(electric))
    {
        return Err(JouleEnthalpyError::Invalid(
            "current and live flow need the same all-fluid domain and exact cell map",
        ).into());
    }
    let source = CurrentFlow {
        flow: flow.clone(), account: Account::new(input, gate), response: &resistivity,
    };
    let (next, mut accepted, audit) =
        thermal.advance_observed(&source, target_s, config, || checkpoint(gate))?;
    let endpoint = evaluate(
        &mut accepted.account, &next, accepted.flow.endpoint(), accepted.response,
    )?;
    let report = accepted.account.finish(audit, endpoint)?;
    checkpoint(gate)?;
    Ok((next, accepted.flow, report))
}

