//! Material-driven momentum supplies the actual mass flux of enthalpy transport.
//!
//! The original SIMPLEC owner solves flow; the original `Advected` observer
//! transports h and sets the combined conduction/advection CFL. This boundary
//! only maps their shared control volumes and converts volume flux to mass.
//! No cell-centred velocity is interpolated into a second, different flux.

/// Source-pinned liquid admission on the existing native advection owner.
pub mod liquid;

/// Local material feedback through the existing momentum and enthalpy owners.
pub mod feedback;

use std::collections::{BTreeMap, BTreeSet};

use fs_conduction::transient::deposition::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, MaterialId, ThermalLaw,
};
use fs_conduction::transient::deposition::advection::{Advected, MassBoundary, MassFlow};
use fs_exec::CancelGate;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

use crate::conjugate::{FlowField, VoxelDomain};
use crate::material::conjugate::MaterialChtFlow;

/// Fixed control-volume and caloric assignment for one frozen momentum field.
#[derive(Clone, Copy)]
pub struct FlowEnthalpyInput<'a> {
    /// One unique active thermal cell for every fluid voxel, in x-fastest order.
    /// Solid voxels must be None. Other thermal cells may exchange heat through
    /// their original contacts but do not receive advective mass or enthalpy.
    pub cell_ids: &'a [Option<CellId>],
    /// Common immutable caloric law and enthalpy datum of the moving material.
    pub material_id: MaterialId,
    /// The actual borrowed law registered for `material_id`, not a second Cp
    /// model. Storage and T(h), including latent heat, remain with that owner.
    pub caloric: &'a dyn ThermalLaw,
    /// Incoming specific enthalpy [J/kg] on each Face3::ALL domain face.
    /// Every face with local inward flux requires Some, including outlet
    /// backflow. Outgoing flux always carries CURRENT local h. Unused Some
    /// values are checked for support; None never fabricates ambient enthalpy.
    pub inlet_specific_h_j_kg: [Option<f64>; 6],
    /// Bound on all voxels, including solids, checked before map allocation.
    pub max_cells: usize,
    /// Bound on nonzero interior AND boundary connections. The thermal owner's
    /// max_links must separately cover interior heat transfers plus other owners.
    pub max_connections: usize,
}

/// Original material source, solved face fluxes and their enthalpy observer.
///
/// This implements `DepositionObserver`, so it composes directly with existing
/// Joule/conduction/material-interface owners. State and audit publication use
/// the original thermal timeline's all-or-nothing candidate semantics.
pub struct MaterialFlowAdvection<'a, F: FluidTransport + ?Sized, O> {
    material: &'a F,
    transport: FluidTransportState,
    domain: &'a VoxelDomain,
    field: &'a FlowField,
    input: FlowEnthalpyInput<'a>,
    advected: Advected<O>,
}

impl<F: FluidTransport + ?Sized, O: Clone> Clone for MaterialFlowAdvection<'_, F, O> {
    fn clone(&self) -> Self {
        Self {
            material: self.material, transport: self.transport, domain: self.domain,
            field: self.field, input: self.input, advected: self.advected.clone(),
        }
    }
}

impl<'a, F: FluidTransport + ?Sized, O> MaterialFlowAdvection<'a, F, O> {
    /// Original concrete source owner, without resampling or type erasure.
    #[must_use]
    pub const fn material(&self) -> &'a F { self.material }
    /// SAME density/viscosity snapshot used to solve momentum.
    #[must_use]
    pub const fn transport(&self) -> &FluidTransportState { &self.transport }
    /// Actual face volume fluxes used to construct the mass-flow graph.
    #[must_use]
    pub const fn field(&self) -> &'a FlowField { self.field }
    /// Existing mass/enthalpy graph, actual transport audit and inner history.
    #[must_use]
    pub const fn advection(&self) -> &Advected<O> { &self.advected }
    /// Return the accepted inner physical history without cloning it.
    pub fn into_inner(self) -> O { self.advected.into_inner() }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        admit_map(thermal, self.domain, self.transport.density_kg_m3(), self.input, check)
    }
}

impl<F: FluidTransport + ?Sized> MaterialChtFlow<'_, F> {
    /// Connect this material-driven momentum solution to the shared enthalpy
    /// timeline. Each UNIQUE fluid face contributes `rho * volume_flux`; the
    /// original upwind observer supplies internal `mass_rate * donor_h` and
    /// boundary enthalpy exactly once. No Cp*T approximation or latent-heat
    /// source is introduced. Volume, mass and caloric-law identity are checked
    /// at construction, every source epoch and every accepted step/event.
    ///
    /// Source callbacks, conduction links and additional histories stay on the
    /// same actual thermal substeps. Use the returned observer directly with
    /// `advance_observed`, or compose it with an existing electrical driver.
    ///
    /// This is one-way frozen, constant-density momentum coupling on stationary
    /// cells. A caloric plateau does not supply phase-dependent mobility, a
    /// solid skeleton, density-change volume flux, or melt/free-surface motion.
    /// Viscous/porous dissipation, pressure work and turbulence heat transport
    /// are not inferred; those require their own admitted thermal owners.
    ///
    /// # Errors
    /// Refuses incompatible geometry/mass/law, missing inlet or backflow h,
    /// unsupported caloric states, nonfinite flux, budgets, cancellation, or a
    /// momentum residual larger than the shared fixed-mass roundoff admission.
    /// Tighten/refine the momentum solve when necessary; no flux is silently
    /// repaired and no mass source is invented to absorb a divergence residual.
    pub fn advect_enthalpy<'s, O>(
        &'s self,
        thermal: &DepositionSimulation<'_>,
        input: FlowEnthalpyInput<'s>,
        inner: O,
        gate: &CancelGate,
    ) -> Result<MaterialFlowAdvection<'s, F, O>, DepositionError> {
        let mut check = || {
            if gate.is_requested() { Err(DepositionError::Cancelled) } else { Ok(()) }
        };
        admit_map(thermal, self.domain(), self.transport().density_kg_m3(), input, &mut check)?;
        let (flows, boundaries) = mass_graph(
            self.domain(), &self.flow().field, self.transport().density_kg_m3(), input,
            &mut check,
        )?;
        let advected = Advected::new(inner, flows, input.max_connections)?
            .with_boundaries(boundaries, input.max_connections)?;
        check()?;
        Ok(MaterialFlowAdvection {
            material: self.material(), transport: *self.transport(), domain: self.domain(),
            field: &self.flow().field, input, advected,
        })
    }
}

impl<F: FluidTransport + ?Sized, O: DepositionObserver> DepositionObserver
    for MaterialFlowAdvection<'_, F, O>
{
    type Error = O::Error;

    fn requires_step_interval(&self) -> bool { self.advected.requires_step_interval() }

    fn step_sources(&mut self, thermal_start: &DepositionSimulation<'_>,
        start_s: f64, dt_s: f64, end_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.advected.step_sources(thermal_start, start_s, dt_s, end_s, check)
    }

    fn begin(&mut self) -> Result<(), Self::Error> { self.advected.begin() }

    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error>
    {
        self.admit(thermal, check)?;
        self.advected.sources(thermal, check)
    }

    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<u64, f64>, Self::Error>
    {
        self.advected.conductances(thermal, check)
    }

    fn heat_transfers(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<Vec<HeatTransfer>, Self::Error>
    {
        self.advected.heat_transfers(thermal, check)
    }

    fn reservoir_conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error>
    {
        self.advected.reservoir_conductances(thermal, check)
    }

    fn limit_timestep_s(&mut self, thermal: &DepositionSimulation<'_>, proposed_s: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<f64, Self::Error>
    {
        self.advected.limit_timestep_s(thermal, proposed_s, check)
    }

    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>)
        -> Result<(), Self::Error>
    {
        self.admit(thermal, &mut || Ok(()))?;
        self.advected.event(event, thermal)
    }

    fn step(&mut self, start_s: f64, duration_s: f64,
        old_temperatures_k: &BTreeMap<CellId, f64>, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error>
    {
        self.admit(thermal, check)?;
        self.advected.step(start_s, duration_s, old_temperatures_k, thermal, check)
    }
}

fn invalid(reason: &'static str) -> DepositionError { DepositionError::InvalidInput(reason) }

fn admit_map(
    thermal: &DepositionSimulation<'_>, domain: &VoxelDomain, density: f64,
    input: FlowEnthalpyInput<'_>, check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), DepositionError> {
    check()?;
    if domain.cell_count() > input.max_cells {
        return Err(DepositionError::Budget("flow/enthalpy voxels"));
    }
    if input.cell_ids.len() != domain.cell_count()
        || !density.is_finite() || density <= 0.0
        || !thermal.uses_material_law(input.material_id, input.caloric)
    {
        return Err(invalid("flow/enthalpy map or original caloric law does not match"));
    }
    for h in input.inlet_specific_h_j_kg.into_iter().flatten() {
        if !h.is_finite() { return Err(DepositionError::Arithmetic); }
        let t = input.caloric.temperature_k(h)?;
        if !t.is_finite() || t <= 0.0 { return Err(DepositionError::MaterialDomain); }
    }
    let mut seen = BTreeSet::new();
    for (c, &id) in input.cell_ids.iter().enumerate() {
        if c % 256 == 0 { check()?; }
        if !domain.is_fluid(c) {
            if id.is_some() { return Err(invalid("solid voxel cannot receive advective mass")); }
            continue;
        }
        let id = id.ok_or(invalid("every fluid voxel needs a thermal cell"))?;
        if !seen.insert(id) { return Err(invalid("fluid voxels cannot share a thermal cell")); }
        let cell = thermal.cells().get(&id).ok_or(invalid("inactive flow/enthalpy cell"))?;
        if cell.material_id != input.material_id {
            return Err(invalid("flowing cells must share one caloric law and enthalpy datum"));
        }
        let volume = domain.volume(c);
        let mass = density * volume;
        if !same_positive(cell.geometry.volume_m3, volume) || !same_positive(cell.mass_kg, mass) {
            return Err(invalid("flow/enthalpy geometry or mass differs from sourced density times volume"));
        }
        let [x, y, z] = domain.coords(c);
        let centre = domain.center(x, y, z);
        let widths = domain.widths(c);
        for a in 0..3 {
            let mid = 0.5 * cell.geometry.start_m[a] + 0.5 * cell.geometry.end_m[a];
            if !mid.is_finite() || !centre[a].is_finite()
                || (mid - centre[a]).abs() > 64.0 * f64::EPSILON * centre[a].abs().max(widths[a])
            {
                return Err(invalid("flow and thermal cell centres do not match"));
            }
        }
    }
    check()
}

fn same_positive(a: f64, b: f64) -> bool {
    a.is_finite() && a > 0.0 && b.is_finite() && b > 0.0
        && (a - b).abs() <= 64.0 * f64::EPSILON * a.max(b)
}

// Unique positive-side interior faces, then actual boundary signs. Stable IDs
// use the owner's existing face-array layout, so zero/reversed flux never
// changes another physical connection's identity.
fn mass_graph(
    domain: &VoxelDomain, field: &FlowField, density: f64, input: FlowEnthalpyInput<'_>,
    check: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(Vec<MassFlow>, Vec<MassBoundary>), DepositionError> {
    if field.dims() != domain.dims() { return Err(invalid("flow field domain differs from thermal map")); }
    let mut flows = Vec::new();
    let mut boundaries = Vec::new();
    let mut rates = vec![[0.0_f64; 2]; domain.cell_count()];
    let arrays = field.face_arrays();
    let offsets = [0, arrays[0].len(), arrays[0].len().checked_add(arrays[1].len())
        .ok_or(DepositionError::Arithmetic)?];
    for c in 0..domain.cell_count() {
        if c % 256 == 0 { check()?; }
        for f in 0..6 {
            let flux = field.outward(domain, c, f);
            if !flux.is_finite() { return Err(DepositionError::Arithmetic); }
            if !domain.is_fluid(c) {
                if flux != 0.0 { return Err(invalid("solved flow crosses a solid voxel")); }
                continue;
            }
            let neighbor = domain.neighbor(c, f);
            if let Some(n) = neighbor {
                if !domain.is_fluid(n) {
                    if flux != 0.0 { return Err(invalid("solved flow crosses a solid boundary")); }
                    continue;
                }
                if f % 2 == 0 { continue; }
            }
            if flux == 0.0 { continue; }
            if flows.len().checked_add(boundaries.len()).is_none_or(|n| n >= input.max_connections) {
                return Err(DepositionError::Budget("flow/enthalpy connections"));
            }
            let rate = density * flux;
            if !rate.is_finite() || rate == 0.0 { return Err(DepositionError::Arithmetic); }
            let (slot, _) = field.slot(domain, c, f);
            let id = u64::try_from(offsets[f / 2].checked_add(slot)
                .ok_or(DepositionError::Arithmetic)?).map_err(|_| DepositionError::Arithmetic)?;
            let cell = input.cell_ids[c].ok_or(invalid("unmapped flowing cell"))?;
            if let Some(n) = neighbor {
                let b = input.cell_ids[n].ok_or(invalid("unmapped flowing neighbor"))?;
                flows.push(MassFlow { id, a: cell, b, kg_s: rate });
                let (up, down) = if rate > 0.0 { (c, n) } else { (n, c) };
                rates[up][1] += rate.abs();
                rates[down][0] += rate.abs();
            } else if rate > 0.0 {
                boundaries.push(MassBoundary::Outlet { id, cell, kg_s: rate });
                rates[c][1] += rate;
            } else {
                let h = input.inlet_specific_h_j_kg[f]
                    .ok_or(invalid("inward flow or outlet backflow needs supported inlet enthalpy"))?;
                boundaries.push(MassBoundary::Inlet {
                    id, cell, kg_s: -rate, material_id: input.material_id, specific_h_j_kg: h,
                });
                rates[c][0] -= rate;
            }
        }
    }
    for [incoming, outgoing] in rates {
        if !incoming.is_finite() || !outgoing.is_finite() { return Err(DepositionError::Arithmetic); }
        if (incoming - outgoing).abs() > 64.0 * f64::EPSILON * incoming.max(outgoing) {
            return Err(invalid("solved flow must conserve fixed cell mass to roundoff; refine momentum solve"));
        }
    }
    check()?;
    Ok((flows, boundaries))
}

#[cfg(test)]
mod tests;
