//! Complete source-bound wet frames in the ORIGINAL shared thermal transaction.
//!
//! Material bindings protect their solid source ports. This extension additionally
//! protects every cold fluid component and reciprocal exterior connection before
//! handing the whole assembly to that same owner. No thermal equation, integrator
//! or loss-to-heat rule is added. Placements retain the physical coordinate maps.

use super::*;
use crate::air_path::cylinder::frame::modal::{FrameFluidConnectionLayout, FrameFluidPlacement};
use core::ops::Range;

/// Cold parts for one source frame at explicit native assembly addresses.
/// Append ALL components/connections unchanged. Retain `placement` for physical
/// ports/observations; raw native states are not dry-basis coordinates. Bind all
/// returned bindings through `FrameFluidThermoelasticBinding::into_shared_system`.
pub struct FrameFluidThermoelasticParts<'m, 'f, 'n, 'a, F> {
    pub placement: FrameFluidPlacement<'m, 'f, F>,
    pub binding: FrameFluidThermoelasticBinding<'n, 'a>,
    pub components: Vec<ModalAcousticTimeModel>,
    pub connections: Vec<ModalConnection>,
}

/// Original solid thermal binding plus its unchanged cold exterior-fluid state.
/// The source-typed placement stays with the caller. The binding has no fluid
/// type parameter, so different original source types can share one assembly.
/// Exterior dissipation remains outside finite solid heat unless a separate
/// physical owner explicitly supplies a justified receiving thermal inventory.
pub struct FrameFluidThermoelasticBinding<'n, 'a> {
    material: FrameThermoelasticBinding<'n, 'a>,
    components: Range<usize>,
    connections: Range<usize>,
    initial_memory: Option<ModalAcousticTimeModel>,
    fluid_start: usize,
    fluid_definitions: Vec<ModalConnection>,
}

impl<'f, F> MaterialFluidFrame<'f, F> {
    /// Prepare the SAME thermoelastic fluid-loaded frame for a larger shared
    /// assembly. Source stiffness, expansion and damping transform together;
    /// original inertia, all fluid momenta, regions and source references survive.
    ///
    /// The component and connection addresses are offsets, not new physics.
    /// Joint springs/contact ports may precede or follow the exported ranges.
    /// Whole-layer thermal destinations remain disjoint between source frames.
    /// `config.prepared` selects only the standalone convenience constructor;
    /// the caller chooses preparation for the final COMPLETE native assembly.
    ///
    /// This cold operation does not advance accepted mechanics or thermal time.
    /// Constant reference mass/geometry, held exterior transport, complete modes
    /// and the existing thermoelastic source/strain/capacity restrictions apply.
    #[allow(clippy::too_many_arguments)]
    pub fn thermoelastic_parts<'m, 'n, 'a>(
        &'m self,
        initial_state: &[f64],
        thermal: &DepositionSimulation<'_>,
        source: &'n MaterialNetwork<'a>,
        destinations: &[FrameLayerCell],
        config: FrameFluidThermoelasticConfig,
        component_start: usize,
        connection_start: usize,
        gate: &CancelGate,
    ) -> FluidResult<FrameFluidThermoelasticParts<'m, 'f, 'n, 'a, F>> {
        check(gate)?;
        positive(config.initial_energy_tolerance_j)?;
        let count = self.basis().modes().len();
        let mut component = self.basis().coupled_component(config.sample_rate_hz,
            &vec![0.0; count], config.modal, config.maximum_projection_terms,
            || gate.is_requested()).map_err(frame_error)?;
        component.restore_states(&self.modal_states(initial_state)?).map_err(frame_error)?;
        let (mut material, _, mut definitions) = component.into_thermoelastic_parts(
            thermal, source, destinations, config.maximum_absolute_strain,
            component_start, connection_start, config.coupling, gate)?;
        let FrameFluidModalParts { coordinates, structural, memory, mut connections } =
            self.modal_parts(initial_state, config.sample_rate_hz, config.modal,
                definitions.len(), config.coupling, gate)?;
        let component_count = 1 + usize::from(memory.is_some());
        let component_end = component_start.checked_add(component_count)
            .ok_or_else(|| invalid("complete wet-frame component address overflow"))?;
        let fluid_start = connection_start.checked_add(definitions.len())
            .ok_or_else(|| invalid("complete wet-frame connection address overflow"))?;
        let connection_end = fluid_start.checked_add(connections.len())
            .ok_or_else(|| invalid("complete wet-frame connection extent overflow"))?;
        if component_end > config.coupling.max_modes || connection_end > config.coupling.max_connections {
            return Err(invalid("complete wet-frame placement exceeds its declared assembly budget").into());
        }
        for definition in &mut definitions {
            check(gate)?;
            definition.left.shapes = coordinates.project_column(&definition.left.shapes)?;
            definition.right.shapes = coordinates.project_column(&definition.right.shapes)?;
        }
        for member in &mut material.members {
            check(gate)?;
            for column in &mut member.columns { *column = coordinates.project_column(column)?; }
        }
        let mut states = Vec::with_capacity(material.members.len());
        for member in &material.members {
            check(gate)?;
            states.push(resolve_member(member, &temperatures(&member.region, thermal)?,
                structural.states(), structural.sample_period_s())?);
        }
        material.states = states;
        material.initial_model = structural.clone();
        material.definitions = definitions.clone();
        for connection in &mut connections {
            check(gate)?;
            connection.left.component += component_start;
            connection.right.component += component_start;
        }
        // The original constructor orders all member viscosities first, then
        // all layer springs. Source spring addresses retain that exact layout.
        let elastic_start = material.members.first()
            .ok_or_else(|| invalid("complete thermoelastic frame has no source members"))?.spring_start;
        let mut region_start = fluid_start;
        let mut fluid_regions = Vec::with_capacity(self.members().len());
        for region in self.members() {
            check(gate)?;
            // The existing cylinder assembly has exactly four original spatial
            // stations per region; each exports its tail plus every memory arm.
            let links = region.flow().branches().len().checked_add(1)
                .and_then(|n| n.checked_mul(4))
                .ok_or_else(|| invalid("regional fluid port extent overflow"))?;
            let end = region_start.checked_add(links)
                .ok_or_else(|| invalid("regional fluid port address overflow"))?;
            fluid_regions.push(region_start..end);
            region_start = end;
        }
        if region_start != connection_end {
            return Err(invalid("source region and complete fluid connection layout differ").into());
        }
        let period = structural.sample_period_s();
        let initial_memory = memory.clone();
        let mut components = vec![structural];
        if let Some(memory) = memory { components.push(memory); }
        // Independently check the source's original canonical H before any
        // external joints are introduced. This includes ALL fluid kinetics.
        let mut composed = sum(components.iter().flat_map(|m| m.states())
            .map(|s| 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2)))?;
        for definition in &definitions {
            check(gate)?;
            let extension = finite(sum(definition.left.shapes.iter()
                .zip(&definition.right.shapes).zip(components[0].states())
                .map(|((&a, &b), s)| (a-b)*s.displacement_m_sqrt_kg))?
                - definition.rest_extension_m)?;
            composed = finite(composed + 0.5*definition.stiffness_n_m*extension*extension)?;
        }
        if finite(composed-self.system().hamiltonian(initial_state))?.abs() > config.initial_energy_tolerance_j {
            return Err(invalid("complete wet thermoelastic parts changed original initial storage").into());
        }
        let placement = coordinates.into_placement(component_start,
            FrameFluidConnectionLayout { material: connection_start..elastic_start,
                elastic: elastic_start..fluid_start, fluid_regions }, period)?;
        let binding = FrameFluidThermoelasticBinding { material,
            components: component_start..component_end,
            connections: connection_start..connection_end,
            initial_memory, fluid_start, fluid_definitions: connections.clone() };
        definitions.extend(connections);
        check(gate)?;
        Ok(FrameFluidThermoelasticParts { placement, binding, components, connections: definitions })
    }
}

impl<'n, 'a> FrameFluidThermoelasticBinding<'n, 'a> {
    /// One-frame specialization of the SAME shared binding/admission path.
    pub fn into_system(self, mechanical: CoupledModalSystem,
        thermal: &DepositionSimulation<'_>, gate: &CancelGate)
        -> FluidResult<FrameThermoelasticSystem<'n, 'a>> {
        Self::into_shared_system(vec![self], mechanical, thermal, gate)
    }

    /// Admit all wet/dry frames to one existing mechanics/finite-heat owner.
    /// Different fluid source types are allowed: each caller-retained placement
    /// still borrows its exact original typed material model. No source object is
    /// replaced by an averaged fluid or a new material database.
    ///
    /// Every protected fluid momentum, gauge, loss column and coefficient must
    /// match its original cold export. Solid source pins, complete layer cells,
    /// common material-network identity and clocks are admitted by the original
    /// shared thermoelastic owner. Overlap or an omitted/changed fluid refuses.
    /// Subsequent native steps, electrical/field/heat-network bindings, entropy,
    /// reversible heat and atomic publication all remain on that original owner.
    pub fn into_shared_system(mut bindings: Vec<Self>, mechanical: CoupledModalSystem,
        thermal: &DepositionSimulation<'_>, gate: &CancelGate)
        -> FluidResult<FrameThermoelasticSystem<'n, 'a>> {
        check(gate)?;
        if bindings.is_empty() || bindings.len() > mechanical.components().len() {
            return Err(invalid("nonempty bounded fluid-frame thermal bindings required").into());
        }
        bindings.sort_unstable_by_key(|b| b.components.start);
        if bindings.windows(2).any(|p| p[0].components.end > p[1].components.start) {
            return Err(invalid("complete fluid-frame component ranges overlap").into());
        }
        let mut ranges: Vec<_> = bindings.iter().map(|b| b.connections.clone()).collect();
        ranges.sort_unstable_by_key(|r| r.start);
        if ranges.windows(2).any(|p| p[0].end > p[1].start) {
            return Err(invalid("complete fluid-frame connection ranges overlap").into());
        }
        for binding in &bindings { binding.admit_fluid(&mechanical, gate)?; }
        let material = bindings.into_iter().map(|b| b.material).collect();
        Ok(FrameThermoelasticBinding::into_shared_system(material, mechanical, thermal, gate)?)
    }

    fn admit_fluid(&self, mechanical: &CoupledModalSystem, gate: &CancelGate) -> FluidResult<()> {
        check(gate)?;
        if mechanical.samples_rendered() != 0 {
            return Err(invalid("fluid thermal binding cannot discard unobserved warm history").into());
        }
        let local = mechanical.components().get(self.components.clone())
            .ok_or_else(|| invalid("complete source fluid component range is absent"))?;
        if let Some(expected) = &self.initial_memory {
            let actual = &local[1];
            if actual.sample_period_s().to_bits() != expected.sample_period_s().to_bits()
                || actual.modes() != expected.modes() || actual.states() != expected.states()
                || actual.budget() != expected.budget() {
                return Err(invalid("original fluid momentum, free gauge, clock or limits changed before binding").into());
            }
        }
        let mut offsets = Vec::with_capacity(mechanical.components().len()+1);
        offsets.push(0usize);
        for component in mechanical.components() {
            check(gate)?;
            let end = offsets.last().unwrap().checked_add(component.modes().len())
                .ok_or_else(|| invalid("complete native fluid column extent overflow"))?;
            offsets.push(end);
        }
        if mechanical.constant_loss_connections().iter()
            .any(|c| (self.fluid_start..self.connections.end).contains(&c.connection)) {
            return Err(invalid("original fluid resistance cannot be replaced by a fitted loss factor").into());
        }
        for (i, expected) in self.fluid_definitions.iter().enumerate() {
            check(gate)?;
            let index = self.fluid_start+i;
            let actual = mechanical.connection_coefficients(index)
                .ok_or_else(|| invalid("source fluid connection is absent"))?;
            let column = mechanical.connection_column(index)
                .ok_or_else(|| invalid("source fluid force column is absent"))?;
            if actual.stiffness_n_m.to_bits() != expected.stiffness_n_m.to_bits()
                || actual.damping_n_s_m.to_bits() != expected.damping_n_s_m.to_bits()
                || actual.rest_extension_m.to_bits() != expected.rest_extension_m.to_bits()
                || column.len() != mechanical.mode_count() {
                return Err(invalid("original fluid connection coefficients or dimensions changed").into());
            }
            let left = offsets[expected.left.component];
            let right = offsets[expected.right.component];
            for (j, &value) in column.iter().enumerate() {
                if j % 128 == 0 { check(gate)?; }
                let a = j.checked_sub(left).and_then(|k| expected.left.shapes.get(k)).copied().unwrap_or(0.0);
                let b = j.checked_sub(right).and_then(|k| expected.right.shapes.get(k)).copied().unwrap_or(0.0);
                if value != a-b {
                    return Err(invalid("original fluid connection was rewired or lost reciprocal force").into());
                }
            }
        }
        check(gate)?;
        Ok(())
    }
}
