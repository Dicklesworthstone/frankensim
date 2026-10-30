//! Complete frozen material frames in the original shared mechanical owner.
//!
//! The existing loaded-inertia chart supplies q=L^-T r. Dry stiffness is
//! sum_i k_i (e_i^T q)^2/2; transform each e_i with L^-1, exactly as every
//! material-viscosity and fluid-traction column. The native components are
//! FREE coordinates, so these springs are the only copy of solid elasticity.
//! Modal spring factors use an explicit one-kilogram coordinate reference:
//! x_i=q_i/sqrt(m_ref), k_i=m_ref*omega_i^2. It is only a units convention,
//! not an additional mass or material property; k_i*x_i^2 is independent of it.
//! No new integrator, mode retuning or diagonal damping approximation occurs.

use super::*;
use crate::render::schedule::force::coupled::{CoupledModalSystem, PreparedModalBudget};
use core::ops::Range;

const MODAL_SPRING_REFERENCE_MASS_KG: f64 = 1.0;

mod ports;
pub use ports::FrameFluidPortMap;

/// Complete original source, normalized motion and every physical connection.
/// Unlike `modal_parts`, this includes solid elasticity AND source viscosity.
/// `into_parts` supports arbitrary native assemblies; `into_system` prepares a
/// standalone instance through that same path. No thermal cells are required.
pub struct FrameFluidComponent<'m, 'f, F> {
    coordinates: FrameFluidCoordinates<'m, 'f, F>,
    components: Vec<ModalAcousticTimeModel>,
    connections: Vec<ModalConnection>,
    layout: FrameFluidConnectionLayout,
    config: ModalCouplingConfig,
}

/// Disjoint ORIGINAL connection accounts. Solid springs store energy; source
/// material ports dissipate it; each exterior region has its own loss range.
/// External assembly connectors are not silently included in these accounts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameFluidConnectionLayout {
    pub elastic: Range<usize>,
    pub material: Range<usize>,
    /// Same order as source().members(); repeated element IDs remain separate.
    pub fluid_regions: Vec<Range<usize>>,
}

/// Placement and coordinate map retained with a complete native assembly.
/// The assembly owns motion. This map owns neither an independent clock nor
/// another accepted history. Preserve EVERY exported component and connection.
/// Changing their coefficients changes the physical model; this map is not
/// authority for arbitrary rewiring or a live material/geometry replacement.
pub struct FrameFluidPlacement<'m, 'f, F> {
    coordinates: FrameFluidCoordinates<'m, 'f, F>,
    components: Range<usize>,
    layout: FrameFluidConnectionLayout,
    period_s: f64,
}

impl<'m, 'f, F> FrameFluidCoordinates<'m, 'f, F> {
    /// Shared cold placement seam for source constitutive owners. The caller
    /// supplies its actual connection layout; the physical coordinate chart
    /// and complete fluid state extent still come from the original source.
    pub(crate) fn into_placement(self, component_start: usize,
        layout: FrameFluidConnectionLayout, period_s: f64)
        -> Result<FrameFluidPlacement<'m, 'f, F>> {
        positive(period_s)?;
        let count = 1 + usize::from(!self.source().memory_masses.is_empty());
        let end = component_start.checked_add(count)
            .ok_or_else(|| invalid("complete fluid placement extent overflow"))?;
        Ok(FrameFluidPlacement { coordinates: self, components: component_start..end,
            layout, period_s })
    }
}

impl<'f, F> MaterialFluidFrame<'f, F> {
    /// Prepare the complete frozen frame/fluid component for the existing
    /// simultaneous time-domain/contact/radiation assembly. Nonzero structural
    /// and fluid velocities survive. The four regional quadrature stations,
    /// loaded inertia, eigenvectors and constitutive sources are reused.
    ///
    /// Retained-mode truncation is the original Galerkin approximation. Native
    /// free fluid position gauges are kept and budgeted, never reset per step.
    /// Original source viscosity must cover the ACTUAL sampled frequency band.
    /// Complete native assembly, temporal and physical-model limits still apply.
    pub fn coupled_component<'m>(
        &'m self,
        initial_state: &[f64],
        sample_rate_hz: u32,
        modal_budget: ModalAcousticTimeBudget,
        config: ModalCouplingConfig,
        gate: &CancelGate,
    ) -> Result<FrameFluidComponent<'m, 'f, F>> {
        poll(gate)?;
        if sample_rate_hz == 0
            || !config.energy_absolute_tolerance_j.is_finite()
            || config.energy_absolute_tolerance_j <= 0.0
            || !config.energy_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&config.energy_relative_tolerance)
        {
            return Err(invalid("complete frame component requires a finite clock and energy allowance"));
        }
        if let Some(damping) = self.basis.material_damping() {
            damping.admit_timestep(1.0 / f64::from(sample_rate_hz))?;
        }
        let m = self.basis.modes().len();
        let solid_count = m.checked_add(self.material_ports.len())
            .ok_or_else(|| invalid("complete solid connection count overflow"))?;
        let FrameFluidModalParts { coordinates, structural, memory, connections: fluid } =
            self.modal_parts(initial_state, sample_rate_hz, modal_budget, solid_count, config, gate)?;
        let capacity = solid_count.checked_add(fluid.len())
            .ok_or_else(|| invalid("complete frame connection count overflow"))?;
        let mut connections = Vec::with_capacity(capacity);
        let mut unit = vec![0.0; m];
        for (i, mode) in self.basis.modes().iter().enumerate() {
            poll(gate)?;
            unit[i] = 1.0;
            let mut shapes = coordinates.project_column(&unit)?;
            for value in &mut shapes { *value /= MODAL_SPRING_REFERENCE_MASS_KG.sqrt(); }
            connections.push(ground(shapes,
                positive(MODAL_SPRING_REFERENCE_MASS_KG * mode.eigenvalue_per_s2)?, 0.0));
            unit[i] = 0.0;
        }
        for port in &self.material_ports {
            poll(gate)?;
            connections.push(ground(coordinates.project_column(&port.projection)?,
                0.0, port.resistance));
        }
        // The original modal export emits one low-tail connection then every
        // spectral branch per station. Retain that order and source partition.
        let mut region_counts = vec![0usize; self.members.len()];
        for station in &self.stations {
            poll(gate)?;
            region_counts[station.member] = region_counts[station.member]
                .checked_add(1 + self.members[station.member].flow.branches().len())
                .ok_or_else(|| invalid("regional connection count overflow"))?;
        }
        let mut end = solid_count;
        let mut fluid_regions = Vec::with_capacity(region_counts.len());
        for count in region_counts {
            let next = end.checked_add(count)
                .ok_or_else(|| invalid("regional connection range overflow"))?;
            fluid_regions.push(end..next);
            end = next;
        }
        if end != capacity {
            return Err(invalid("original fluid ports and regional history layout disagree"));
        }
        connections.extend(fluid);
        let mut components = vec![structural];
        if let Some(memory) = memory { components.push(memory); }
        // Independent initial energy in exported physical springs and free
        // kinetics, compared with the original canonical frame/fluid H.
        let kinetic = components.iter().flat_map(|c| c.states()).try_fold(0.0,
            |sum, s| finite(sum + 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2)))?;
        let mut composed = kinetic;
        for connection in &connections[..m] {
            poll(gate)?;
            let extension = connection.left.shapes.iter().zip(components[0].states())
                .try_fold(0.0, |sum, (&b, s)| finite(sum + b * s.displacement_m_sqrt_kg))?;
            composed = finite(composed + 0.5 * connection.stiffness_n_m * extension * extension)?;
        }
        let original = finite(self.system.hamiltonian(initial_state))?;
        let allowance = finite(config.energy_absolute_tolerance_j
            + config.energy_relative_tolerance * original.abs().max(composed.abs()))?;
        if finite(composed - original)?.abs() > allowance {
            return Err(invalid("complete frame component changed original initial stored energy"));
        }
        poll(gate)?;
        Ok(FrameFluidComponent { coordinates, components, connections,
            layout: FrameFluidConnectionLayout { elastic: 0..m, material: m..solid_count, fluid_regions },
            config })
    }
}

fn ground(shapes: Vec<f64>, stiffness: f64, damping: f64) -> ModalConnection {
    let n = shapes.len();
    ModalConnection {
        left: ModalAttachment { component: 0, shapes },
        right: ModalAttachment { component: 0, shapes: vec![0.0; n] },
        stiffness_n_m: stiffness, damping_n_s_m: damping, rest_extension_m: 0.0,
    }
}

impl<'m, 'f, F> FrameFluidComponent<'m, 'f, F> {
    pub const fn source(&self) -> &'m MaterialFluidFrame<'f, F> { self.coordinates.source() }
    pub fn component_count(&self) -> usize { self.components.len() }
    pub fn connection_count(&self) -> usize { self.connections.len() }

    /// Append the returned components/links at the declared offsets, preserving
    /// their full order. Other components and connectors may precede/follow.
    /// No warm rebind, gauge reset or source transition is performed. Final
    /// native construction admits the WHOLE assembly and its separate budget.
    pub fn into_parts(mut self, component_start: usize, connection_start: usize,
        gate: &CancelGate) -> Result<(FrameFluidPlacement<'m, 'f, F>,
            Vec<ModalAcousticTimeModel>, Vec<ModalConnection>)> {
        poll(gate)?;
        let component_end = component_start.checked_add(self.components.len())
            .ok_or_else(|| invalid("frame component address overflow"))?;
        let connection_end = connection_start.checked_add(self.connections.len())
            .ok_or_else(|| invalid("frame connection address overflow"))?;
        // A component has at least one mode, so max_modes also bounds its
        // possible component address. This does not replace whole-assembly admission.
        if component_end > self.config.max_modes || connection_end > self.config.max_connections {
            return Err(invalid("frame placement exceeds its declared assembly limits"));
        }
        for connection in &mut self.connections {
            poll(gate)?;
            connection.left.component += component_start;
            connection.right.component += component_start;
        }
        let shift = |range: &mut Range<usize>| {
            range.start += connection_start;
            range.end += connection_start;
        };
        shift(&mut self.layout.elastic);
        shift(&mut self.layout.material);
        for region in &mut self.layout.fluid_regions { shift(region); }
        let period_s = self.components[0].sample_period_s();
        let placement = self.coordinates.into_placement(component_start, self.layout, period_s)?;
        poll(gate)?;
        Ok((placement, self.components, self.connections))
    }

    /// Standalone use of the SAME complete native owner. For cross-body,
    /// contact or radiation composition use `into_parts`, not separate time
    /// steps followed by exchanging endpoint forces. Prepared propagation is
    /// optional and retains the original native work/capacity limits.
    pub fn into_system(self, prepared: Option<PreparedModalBudget>, gate: &CancelGate)
        -> Result<(FrameFluidPlacement<'m, 'f, F>, CoupledModalSystem)> {
        let config = self.config;
        let (placement, components, connections) = self.into_parts(0, 0, gate)?;
        let system = match prepared {
            Some(budget) => CoupledModalSystem::new_prepared(components, connections, config, budget, gate),
            None => CoupledModalSystem::new(components, connections, config, gate),
        }.map_err(FrameVibrationError::from)?;
        poll(gate)?;
        Ok((placement, system))
    }
}

impl<'m, 'f, F> FrameFluidPlacement<'m, 'f, F> {
    pub(crate) fn into_coordinates(self) -> FrameFluidCoordinates<'m, 'f, F> { self.coordinates }
    pub const fn source(&self) -> &'m MaterialFluidFrame<'f, F> { self.coordinates.source() }
    pub fn component_range(&self) -> Range<usize> { self.components.clone() }
    pub fn structural_component(&self) -> usize { self.components.start }
    pub fn fluid_component(&self) -> Option<usize> {
        (self.components.len() == 2).then_some(self.components.start + 1)
    }
    pub const fn connections(&self) -> &FrameFluidConnectionLayout { &self.layout }

    /// Physical nodal force/moment -> this placement's structural force slots.
    /// Append zeros for its fluid memory in a complete native forcing vector.
    /// No force is applied to the private fluid position gauges.
    pub fn project_nodal_force(&self, nodal: &[f64], gate: &CancelGate) -> Result<Vec<f64>> {
        self.coordinates.project_nodal_force(nodal, gate)
    }

    /// Convert this placement's complete state from a native component slice.
    /// Pass the ORIGINAL exported component order, as retained by the assembly
    /// or its contact wrapper. Dimensions/clock/free charts are checked; this
    /// coordinate map does not authorize changing the exported connections.
    pub fn physical_state(&self, components: &[ModalAcousticTimeModel], gate: &CancelGate) -> Result<Vec<f64>> {
        poll(gate)?;
        let local = self.local_components(components)?;
        let motion = self.coordinates.original_states(local[0].states())?;
        let mut velocities = Vec::with_capacity(self.source().memory_masses.len());
        if let Some(memory) = local.get(1) {
            for (state, &mass) in memory.states().iter().zip(&self.source().memory_masses) {
                poll(gate)?;
                velocities.push(finite(state.velocity_m_sqrt_kg_per_s / mass.sqrt())?);
            }
        }
        let state = self.source().state_from_motion(&motion, &velocities)?;
        poll(gate)?;
        Ok(state)
    }

    /// Original node-major physical observations, with the source basis's
    /// reference displacement convention. Output arrays are scratch: a failed
    /// reconstruction may overwrite them, but never accepted physical history.
    pub fn reconstruct_into(&self, components: &[ModalAcousticTimeModel], displacement: &mut [f64],
        velocity: &mut [f64], gate: &CancelGate) -> Result<()> {
        poll(gate)?;
        let local = self.local_components(components)?;
        let motion = self.coordinates.original_states(local[0].states())?;
        self.source().basis.reconstruct_into(&motion, displacement, velocity,
            || gate.is_requested())?;
        poll(gate)
    }

    fn local_components<'c>(&self, components: &'c [ModalAcousticTimeModel])
        -> Result<&'c [ModalAcousticTimeModel]> {
        let local = components.get(self.components.clone())
            .ok_or_else(|| invalid("original complete frame placement is absent"))?;
        for (index, component) in local.iter().enumerate() {
            let count = if index == 0 { self.source().basis.modes().len() }
                else { self.source().memory_masses.len() };
            if component.modes().len() != count || component.states().len() != count
                || component.sample_period_s().to_bits() != self.period_s.to_bits()
                || component.modes().iter().any(|m| m.angular_frequency_rad_s != 0.0
                    || m.damping_ratio != 0.0 || m.pressure_per_modal_velocity != C64::ZERO)
            {
                return Err(invalid("frame placement lost its original free coordinates, clock or state dimension"));
            }
        }
        Ok(local)
    }
}

#[cfg(test)]
mod tests;
