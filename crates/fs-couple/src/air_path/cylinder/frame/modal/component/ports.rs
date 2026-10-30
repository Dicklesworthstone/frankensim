//! The SAME original physical ports in harmonic and transient assemblies.
use super::*;
use crate::air_path::cylinder::frame::harmonic_ports::FrameHarmonicPorts;

/// Cached full-inertia maps for one original physical port set and placement.
/// Repeated, signed and dependent columns are kept. No new port geometry,
/// scalar material coefficient, observer gain or pressure transfer is inferred.
pub struct FrameFluidPortMap<'p, 'm, 'f, F> {
    placement: &'p FrameFluidPlacement<'m, 'f, F>,
    original: &'p FrameHarmonicPorts<'m, 'f, F>,
    shapes: Vec<Vec<f64>>,
    reference: Vec<f64>,
}

impl<'m, 'f, F> FrameFluidPlacement<'m, 'f, F> {
    /// Reuse a harmonic port definition for physical time-domain forcing,
    /// observation and native reciprocal attachments. Its source must be the
    /// EXACT same frame/fluid model, not merely a same-sized basis.
    ///
    /// `maximum_projection_terms` bounds ports*(modes^2+nodal_dofs) before
    /// allocating the transformed maps. Cold triangular transforms occur once;
    /// subsequent force and observation calls allocate no owner workspace.
    pub fn bind_ports<'p>(&'p self, ports: &'p FrameHarmonicPorts<'m, 'f, F>,
        maximum_projection_terms: usize, gate: &CancelGate)
        -> Result<FrameFluidPortMap<'p, 'm, 'f, F>> {
        poll(gate)?;
        if !core::ptr::eq(self.source(), ports.model()) {
            return Err(invalid("time and harmonic ports must retain the SAME original source model"));
        }
        let m = self.source().basis.modes().len();
        let reference = self.source().basis.reference_displacement();
        if m.checked_mul(m).and_then(|n| n.checked_add(reference.len()))
            .and_then(|n| n.checked_mul(ports.port_count()))
            .is_none_or(|n| n > maximum_projection_terms) {
            return Err(invalid("complete source port coordinate projection exceeds budget"));
        }
        let mut shapes = Vec::with_capacity(ports.port_count());
        let mut datum = Vec::with_capacity(ports.port_count());
        for (modal, nodal) in ports.modal_columns().iter().zip(ports.nodal_columns()) {
            poll(gate)?;
            shapes.push(self.coordinates.project_column(modal)?);
            datum.push(dot(nodal, reference)?);
        }
        poll(gate)?;
        Ok(FrameFluidPortMap { placement: self, original: ports, shapes, reference: datum })
    }
}

impl<'p, 'm, 'f, F> FrameFluidPortMap<'p, 'm, 'f, F> {
    pub fn port_count(&self) -> usize { self.shapes.len() }
    pub const fn original(&self) -> &'p FrameHarmonicPorts<'m, 'f, F> { self.original }
    /// Fixed reference coordinates in each original port's declared units.
    /// Native attachments act on incremental motion; include these datums when
    /// specifying a physical spring rest extension, support position or gap.
    pub fn reference_coordinates(&self) -> &[f64] { &self.reference }

    /// Native signed attachment on the structural component, never on a fluid
    /// gauge. Source transverse ports have force/displacement SI units. A
    /// general user column retains its declared conjugate effort/motion units:
    /// only connect compatible coordinates or an explicit physical transformer.
    /// In particular, an unscaled angle is not a displacement in metres.
    pub fn attachment(&self, port: usize) -> Result<ModalAttachment> {
        let shapes = self.shapes.get(port)
            .ok_or_else(|| invalid("original physical port index is absent"))?;
        Ok(ModalAttachment { component: self.placement.structural_component(), shapes: shapes.clone() })
    }

    /// Add coherent physical efforts into the COMPLETE native force vector.
    /// Existing other-body inputs are retained; fluid-memory slots are untouched.
    /// Call once per body after initializing the shared vector. No accepted
    /// state is changed; this force buffer is scratch and may be partial on Err.
    pub fn accumulate_forces(&self, components: &[ModalAcousticTimeModel], efforts: &[f64],
        force: &mut [f64], gate: &CancelGate) -> Result<()> {
        poll(gate)?;
        self.placement.local_components(components)?;
        if efforts.len() != self.port_count() || efforts.iter().any(|x| !x.is_finite())
            || force.iter().any(|x| !x.is_finite()) {
            return Err(invalid("complete finite physical efforts and native force scratch required"));
        }
        let mut offset = 0usize;
        let mut total = 0usize;
        for (i, component) in components.iter().enumerate() {
            poll(gate)?;
            if i == self.placement.structural_component() { offset = total; }
            total = total.checked_add(component.modes().len())
                .ok_or_else(|| invalid("complete native force extent overflow"))?;
        }
        if force.len() != total { return Err(invalid("native forcing must retain ALL component slots")); }
        let target = &mut force[offset..offset+self.source_modes()];
        for (shape, &effort) in self.shapes.iter().zip(efforts) {
            poll(gate)?;
            for (value, &b) in target.iter_mut().zip(shape) { *value = finite(*value + b*effort)?; }
        }
        poll(gate)
    }

    /// Observe the SAME conjugate coordinates without reconstructing every
    /// structural node. Displacements include the retained reference datum;
    /// velocities do not. Output arrays are scratch on refusal/cancellation.
    pub fn observe_into(&self, components: &[ModalAcousticTimeModel], displacement: &mut [f64],
        velocity: &mut [f64], gate: &CancelGate) -> Result<()> {
        poll(gate)?;
        let local = self.placement.local_components(components)?;
        if displacement.len() != self.port_count() || velocity.len() != self.port_count() {
            return Err(invalid("one displacement and velocity output per original physical port required"));
        }
        for (i, shape) in self.shapes.iter().enumerate() {
            poll(gate)?;
            let mut q = self.reference[i];
            let mut v = 0.0;
            for (&b, state) in shape.iter().zip(local[0].states()) {
                q = finite(q + b*state.displacement_m_sqrt_kg)?;
                v = finite(v + b*state.velocity_m_sqrt_kg_per_s)?;
            }
            displacement[i] = q;
            velocity[i] = v;
        }
        poll(gate)
    }

    fn source_modes(&self) -> usize { self.placement.source().basis.modes().len() }
}
