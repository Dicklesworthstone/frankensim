//! Finite solid-layer heating on the complete structural/fluid trajectory.
//!
//! Reuses the original frame-layer/card/caloric binding and the native cell
//! enthalpy/Newton owner. Every thermal trial resolves the original solid
//! viscosities and reruns the SAME prepared mechanical state with ALL fluid
//! momenta intact. Layer heat uses that trial's actual discrete-effort losses
//! and fresh eta*A / eta*I partitions. No coefficient-lag fallback or second
//! thermal/mechanical integrator is introduced.
//!
//! Solid density, elastic energy and reference geometry stay fixed. The fluid
//! remains at its prepared externally held conditions; its viscous heat leaves
//! through an explicit exterior sink. An unconfined fluid is not assigned an
//! invented finite thermal volume. Endpoint viscosity plus thermal backward
//! Euler is first order in general; refine time, dry modes and fluid spectrum.

use super::*;
use crate::air_path::cylinder::frame::{FrameFluidError, MaterialFluidFrame};
use crate::air_path::cylinder::frame::viscosity::{PreparedFluidFrame, PreparedFrameFluidStep};
use std::cell::RefCell;

fn fluid(error: FrameFluidError) -> FrameThermalError {
    match error {
        FrameFluidError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(FrameThermalSourceError::Fluid(other)),
    }
}

/// One accepted frame/fluid history, source-pinned finite solid cells, and
/// reusable native solver scratch. The original thermal owner keeps its cells,
/// links, reservoirs and event cursor; `step` publishes both histories together.
pub struct FrameFluidThermalSystem<'f, 's, 'n, 'a, F> {
    mechanical: PreparedFluidFrame<'f, 's, F>,
    state: Vec<f64>,
    time_s: f64,
    prepared: PreparedMaterialNetwork<'n, 'a>,
    regions: Vec<Region<'a>>,
    viscosity: Vec<(usize, FrameViscousSection)>,
}

/// Actual source viscosities and heat/work accounts for ONE accepted interval.
pub struct FrameFluidThermalReport {
    pub mechanical: PreparedFrameFluidStep,
    pub thermal: EnergyAudit,
    pub source: MaterialNetworkSample,
    pub layer_heat: Vec<FrameLayerHeat>,
    /// Original sections at the accepted thermal endpoint, used during this step.
    pub viscosity: Vec<(usize, FrameViscousSection)>,
    pub deposited_solid_heat_j: f64,
    /// Unconfined fluid loss goes to its declared external thermal reservoir.
    /// It is not deposited a second time into the solid caloric cells.
    pub exported_fluid_heat_j: f64,
    pub external_thermal_heat_j: f64,
    pub combined_energy_residual_j: f64,
    /// Entropy of the finite solid thermal network only; no exterior-fluid
    /// entropy or global second-law claim is inferred from the exported heat.
    pub entropy: FrameHeatEntropy,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    pub material_trials: usize,
}

impl<'f, 's, 'n, 'a, F> FrameFluidThermalSystem<'f, 's, 'n, 'a, F> {
    /// Bind EVERY original viscous solid layer to its actual whole caloric cell.
    /// Reuses `bind_regions`: same immutable material card, exact source context,
    /// mass, volume, phase chart/law instance and constant density/E admission.
    /// A pre-existing nonzero mechanical/fluid state is retained, not reinitialized.
    pub fn new(
        model: &'f MaterialFluidFrame<'s, F>,
        state: Vec<f64>,
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        destinations: &[FrameLayerCell],
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        model.modal_states(&state).map_err(fluid)?;
        let (prepared, regions, viscosity) =
            bind_regions(model.basis(), thermal, source, destinations, true, gate)?;
        let mechanical = model.prepare(gate).map_err(fluid)?;
        check(gate)?;
        Ok(Self {mechanical, state, time_s: thermal.time_s(), prepared, regions, viscosity})
    }

    pub fn model(&self) -> &'f MaterialFluidFrame<'s, F> { self.mechanical.model() }
    pub fn state(&self) -> &[f64] { &self.state }
    pub const fn time_s(&self) -> f64 { self.time_s }
    pub fn viscous_sections(&self) -> &[(usize, FrameViscousSection)] { &self.viscosity }
    pub fn condensed_dimension(&self) -> Option<usize> { self.mechanical.condensed_dimension() }

    /// Solve actual solid heat, source viscosity and loaded motion on one clock.
    ///
    /// Every Newton/stencil/final trial starts from identical accepted mechanics,
    /// fluid memory and cell enthalpy. All source, phase, inventory, local/global
    /// energy and cancellation checks precede joint publication. Source/event
    /// discontinuities require interval cuts; target events stay pending under
    /// the existing thermal residual's rules. Failed attempts change scratch only.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        target_s: f64,
        nodal_forces_and_moments: &[f64],
        controls: FrameThermalSolve,
        initial_energy_increments_j: &[f64],
        gate: &CancelGate,
    ) -> Result<FrameFluidThermalReport> {
        check(gate)?;
        if thermal.time_s().to_bits() != self.time_s.to_bits()
            || !target_s.is_finite() || target_s <= self.time_s
            || !controls.total_energy_tolerance_j.is_finite() || controls.total_energy_tolerance_j <= 0.0 {
            return Err(invalid("matching frame/fluid/thermal clocks, positive interval and work tolerance required"));
        }
        let dt = finite(target_s-self.time_s)?;
        // Reject an independently altered caloric history, even when its clock
        // matches. The accepted source snapshots must describe its actual T.
        for (region, (_, section)) in self.regions.iter().zip(&self.viscosity) {
            for (receiver, layer) in region.receivers.iter().zip(section.elastic().layers()) {
                check(gate)?;
                let t = thermal.temperature_k(receiver.binding.cell)?;
                let source_t = layer.material.query_point().iter()
                    .find(|(axis,_)| axis == &receiver.temperature_axis).map(|(_,v)| *v)
                    .ok_or_else(|| invalid("accepted layer viscosity lost its original temperature"))?;
                if t.to_bits() != source_t.to_bits() {
                    return Err(invalid("accepted frame layer and thermal history do not share one source state"));
                }
            }
        }
        let regions = &self.regions;
        let previous = &self.viscosity;
        let prepared = &self.prepared;
        let state = &self.state;
        let scratch = RefCell::new(&mut self.mechanical);
        let (next_thermal, (), audit, entropy) = cell_driver_candidate_entropy_option(
            thermal, target_s, controls.driver(), initial_energy_increments_j,
            |trial, _remaining| {
                check(gate)?;
                let source = prepared.resolve(trial, || gate.is_requested()).map_err(network)?;
                let viscosity = resolve_sections(regions, previous, trial, dt, gate)?;
                let resistances = original_resistances(regions, &viscosity)?;
                let mechanical = scratch.borrow_mut().step_with_material_resistances(
                    &resistances, state, nodal_forces_and_moments, dt,
                    controls.total_energy_tolerance_j, gate).map_err(fluid)?;
                let layer_heat = partition_heat(regions, &viscosity,
                    &mechanical.material_port_dissipation_j)?;
                let material_heat_j = finite(mechanical.material_port_dissipation_j.iter().sum())?;
                let outside_loss_j = finite(mechanical.fluid_member_dissipation_j.iter().sum())?;
                Ok(DrivenCoefficients {
                    circuit: FluidMechanicalCandidate {mechanical, viscosity, layer_heat,
                        material_heat_j, outside_loss_j},
                    sources_w: source.sources_w,
                    conductances_w_k: source.conductances_w_k,
                    source: source.source,
                    circuit_trials: 1,
                })
            },
            |candidate, duration, checkpoint| {
                let mut powers = BTreeMap::new();
                for heat in &candidate.layer_heat {
                    checkpoint()?;
                    let power = finite(heat.total_j()/duration)?;
                    if heat.total_j() != 0.0 && power == 0.0 {
                        return Err(invalid("solid heat disappeared during conversion to cell power"));
                    }
                    // The original binding owns distinct whole-layer cells.
                    if powers.insert(heat.destination.cell, power).is_some() {
                        return Err(invalid("duplicate original solid heat destination"));
                    }
                }
                Ok(powers)
            },
            true,
            || gate.is_requested(),
        )?;
        let entropy = entropy.ok_or_else(|| invalid("joint frame/fluid heat entropy account is absent"))?;
        let driver = audit.electrical;
        // Prepare all fallible/allocation work before the final publication gate.
        let next_state = driver.mechanical.state.clone();
        let next_viscosity = driver.viscosity.clone();
        let report = FrameFluidThermalReport {
            mechanical: driver.mechanical,
            thermal: audit.thermal,
            source: audit.source,
            layer_heat: driver.layer_heat,
            viscosity: driver.viscosity,
            deposited_solid_heat_j: audit.deposited_joule_j,
            exported_fluid_heat_j: driver.outside_loss_j,
            external_thermal_heat_j: audit.non_electrical_heat_j,
            combined_energy_residual_j: audit.total_energy_residual_j,
            entropy: entropy.into(),
            nonlinear: audit.nonlinear,
            maximum_cell_residual_j: audit.maximum_cell_residual_j,
            material_trials: audit.circuit_trials,
        };
        drop(scratch);
        check(gate)?;
        self.state = next_state;
        self.viscosity = next_viscosity;
        self.time_s = target_s;
        *thermal = next_thermal;
        Ok(report)
    }
}

fn resolve_sections(regions: &[Region<'_>], previous: &[(usize,FrameViscousSection)],
    state: &DepositionSimulation<'_>, dt: f64, gate: &CancelGate)
    -> Result<Vec<(usize,FrameViscousSection)>> {
    if previous.len() != regions.len() { return Err(invalid("original solid region layout changed")); }
    let mut sections = Vec::with_capacity(regions.len());
    for (region, (element, old)) in regions.iter().zip(previous) {
        check(gate)?;
        if *element != region.element { return Err(invalid("original solid member order changed")); }
        let mut temperatures = Vec::with_capacity(region.receivers.len());
        for receiver in &region.receivers {
            let original = &receiver.original;
            let cell = state.cells().get(&receiver.binding.cell)
                .ok_or_else(|| invalid("frame heat receiver disappeared"))?;
            if cell.id != original.id || cell.material_id != original.material_id
                || cell.segment_id != original.segment_id || cell.body_id != original.body_id
                || cell.mass_kg.to_bits() != original.mass_kg.to_bits() || cell.geometry != original.geometry {
                return Err(invalid("frame heat receiver changed original material inventory or geometry"));
            }
            temperatures.push(state.temperature_k(cell.id)?);
        }
        let section = region.law.at_temperatures_k(&temperatures).map_err(material)?;
        let section = FrameViscousSection::try_new(section,
            old.axial().viscosity_property(), old.bending().viscosity_property())
            .map_err(|e| CircuitEnthalpyError::Source(FrameThermalSourceError::Frame(e.into())))?;
        let band = section.omega_band_rad_s();
        if band.0 != 0.0 || band.1 < core::f64::consts::PI/dt {
            return Err(invalid("refreshed solid viscosity does not admit the actual sampled band"));
        }
        sections.push((region.element,section));
    }
    Ok(sections)
}

// SAME axial / symmetric-bending / antisymmetric-bending factors and original
// member length as FrameMaterialDamping and the existing thermal frame owner.
fn original_resistances(regions: &[Region<'_>], sections: &[(usize,FrameViscousSection)])
    -> Result<Vec<f64>> {
    let mut values = Vec::with_capacity(3*regions.len());
    for (region, (_, section)) in regions.iter().zip(sections) {
        let axial = finite(section.axial().viscous_rigidity_n_s()/region.length_m)?;
        let bending = finite(section.bending().viscous_stiffness_n_m2_s()
            /region.length_m/region.length_m/region.length_m)?;
        if section.axial().viscous_rigidity_n_s()>0.0 && axial==0.0
            || section.bending().viscous_stiffness_n_m2_s()>0.0 && bending==0.0 {
            return Err(invalid("original source viscosity underflow on the physical member geometry"));
        }
        values.extend([axial,finite(3.0*bending)?,bending]);
    }
    Ok(values)
}

fn partition_heat(regions: &[Region<'_>], sections: &[(usize,FrameViscousSection)],
    port_heat: &[f64]) -> Result<Vec<FrameLayerHeat>> {
    if port_heat.len() != 3*regions.len() || sections.len() != regions.len()
        || port_heat.iter().any(|h| !h.is_finite() || *h<0.0) {
        return Err(invalid("complete nonnegative original solid-port losses required"));
    }
    let mut result = Vec::new();
    for (index,(region,(_,section))) in regions.iter().zip(sections).enumerate() {
        for (layer,receiver) in region.receivers.iter().enumerate() {
            let a = share(section.axial().layer_viscous_rigidity_n_s()[layer],
                section.axial().viscous_rigidity_n_s())?;
            let b = share(section.bending().layer_viscous_stiffness_n_m2_s()[layer],
                section.bending().viscous_stiffness_n_m2_s())?;
            let partition = |fraction:f64,energy:f64| -> Result<f64> {
                let value = finite(fraction*energy)?;
                if fraction>0.0 && energy>0.0 && value==0.0 {
                    return Err(invalid("physical solid-layer heat underflow"));
                }
                Ok(value)
            };
            result.push(FrameLayerHeat {destination:receiver.binding,
                axial_heat_j:partition(a,port_heat[3*index])?,
                bending_heat_j:partition(b,finite(port_heat[3*index+1]+port_heat[3*index+2])?)?});
        }
    }
    Ok(result)
}

struct FluidMechanicalCandidate {
    mechanical: PreparedFrameFluidStep,
    viscosity: Vec<(usize,FrameViscousSection)>,
    layer_heat: Vec<FrameLayerHeat>,
    material_heat_j: f64,
    outside_loss_j: f64,
}
impl ThermalDriver for FluidMechanicalCandidate {
    type Model = ();
    type Step = Self;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {resistors:&[],contacts:&[],heat_j:self.material_heat_j,
            numerical_dissipation_j:0.0,
            stored_change_j:self.mechanical.stored_after_j-self.mechanical.stored_before_j,
            supplied_j:self.mechanical.nodal_work_j-self.outside_loss_j}
    }
    fn into_parts(self) -> ((),Self) { ((),self) }
}

#[cfg(test)]
mod tests;
