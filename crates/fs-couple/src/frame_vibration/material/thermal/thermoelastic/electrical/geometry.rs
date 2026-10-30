//! Accepted frame kinematics lower original resistor material through actual geometry.
//!
//! The frame's one-dimensional thermoelastic model owns axial extension. Each
//! concentric layer is traction-free transversely in that model, so its area
//! follows only the same sourced scalar free-expansion strain. Poisson
//! contraction, bending-current paths, changing contacts and motional EMF are
//! deliberately not inferred. Geometry is frozen over one existing first-order
//! mechanical/electrical interval and refreshed from the next accepted state.
use super::*;
use crate::material_circuit::{MaterialCircuitError, MaterialResistor};
use fs_material::conductor::ResolvedConductorState;
use std::collections::BTreeSet;

/// One whole-layer current path through an original thermoelastic frame member.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameConductorPath {
    /// Original circuit branch carrying this whole-layer uniform resistor.
    pub branch: usize,
    /// Index in FrameThermoelasticSystem::members().
    pub member: usize,
    /// Ordered concentric layer in that member.
    pub layer: usize,
}

/// Geometry actually used by one electrical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrameConductorGeometry {
    pub branch: usize,
    pub member: usize,
    pub layer: usize,
    pub length_m: f64,
    pub area_m2: f64,
    pub axial_extension_m: f64,
    pub transverse_free_strain: f64,
}

#[derive(Clone, Debug)]
pub(super) struct BoundConductorPath {
    branch: usize,
    member: usize,
    layer: usize,
    component: usize,
    axial_column: Vec<f64>,
    reference_length_m: f64,
    reference_area_m2: f64,
    linear_strain_limit: f64,
}

fn close(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a - b).abs() <= 128.0 * f64::EPSILON * a.abs().max(b.abs())
}

fn sample(
    path: &BoundConductorPath,
    states: &[FrameThermoelasticMemberState],
    models: &[ModalAcousticTimeModel],
) -> Result<FrameConductorGeometry> {
    let material = states.get(path.member)
        .ok_or_else(|| invalid("deforming conductor lost its original frame member"))?;
    if material.component != path.component {
        return Err(invalid("deforming conductor member changed native component"));
    }
    let model = models.get(path.component)
        .ok_or_else(|| invalid("deforming conductor native component is absent"))?;
    if path.axial_column.len() != model.states().len() {
        return Err(invalid("deforming conductor axial coordinate changed dimension"));
    }
    let extension = sum(path.axial_column.iter().zip(model.states())
        .map(|(&shape, state)| shape * state.displacement_m_sqrt_kg))?;
    let strain = finite(extension / path.reference_length_m)?;
    if strain.abs() > path.linear_strain_limit {
        return Err(invalid("deforming conductor exceeds the frame axial strain envelope"));
    }
    let length = finite(path.reference_length_m + extension)?;
    let free = material.source.section.expansion_paths()
        .get(path.layer)
        .ok_or_else(|| invalid("deforming conductor lost its original expansion layer"))?
        .free_linear_strain();
    let transverse = finite(1.0 + free)?;
    let area = finite(path.reference_area_m2 * transverse * transverse)?;
    if length <= 0.0 || transverse <= 0.0 || area <= 0.0 {
        return Err(invalid("deforming conductor geometry must remain positive"));
    }
    Ok(FrameConductorGeometry {
        branch: path.branch, member: path.member, layer: path.layer,
        length_m: length, area_m2: area, axial_extension_m: extension,
        transverse_free_strain: free,
    })
}

pub(super) fn apply(
    geometry: &[FrameConductorGeometry],
    resistors: &mut [MaterialResistor],
    gate: &CancelGate,
) -> Result<()> {
    for geometry in geometry {
        check(gate)?;
        let resistor = resistors.iter_mut().find(|r| r.branch == geometry.branch)
            .ok_or_else(|| invalid("deforming conductor disappeared from endpoint material state"))?;
        resistor.conductor = ResolvedConductorState::try_from_material_with_property(
            resistor.conductor.material().clone(),
            resistor.conductor.source_basis(),
            resistor.conductor.source_key(),
        )
        .and_then(|state| state.into_uniform_conductor(geometry.length_m, geometry.area_m2))
        .map_err(|error| CircuitEnthalpyError::Circuit(MaterialCircuitError::from(error)))?;
    }
    check(gate)
}

pub(super) fn admit_heat_resistors(
    binding: &Binding<'_, '_>,
    resistors: &mut [MaterialResistor],
    gate: &CancelGate,
) -> Result<()> {
    apply(&binding.conductor_geometry, resistors, gate)?;
    resistors.sort_unstable_by_key(|r| r.branch);
    if resistors.len() != binding.circuit.resistors().len() {
        return Err(invalid("complete heat network changed resistor coverage"));
    }
    let deformed: BTreeSet<_> = binding.conductor_geometry.iter().map(|g| g.branch).collect();
    for (fresh, accepted) in resistors.iter().zip(binding.circuit.resistors()) {
        check(gate)?;
        let a = &fresh.conductor;
        let b = &accepted.conductor;
        let same_source = fresh.branch == accepted.branch
            && a.source_basis() == b.source_basis()
            && a.source_key() == b.source_key()
            && a.material().card_identity() == b.material().card_identity()
            && a.material().query() == b.material().query()
            && a.source_property().answer().receipt.selected
                == b.source_property().answer().receipt.selected;
        if !same_source || (!deformed.contains(&fresh.branch) && fresh != accepted) {
            return Err(invalid("complete heat network changed original resistor source or unowned geometry"));
        }
    }
    check(gate)
}

impl FrameThermoelasticSystem<'_, '_> {
    /// Couple selected whole-layer resistors to accepted frame kinematics.
    ///
    /// Axial length follows the member's physical extension. Cross-sectional
    /// area follows the same layer's sourced scalar free thermal expansion in
    /// the two transverse directions, consistent with this 1D frame law. The
    /// electrical material bundle still supplies resistivity/conductivity at
    /// the thermal endpoint. No Poisson strain or electromagnetic force is
    /// invented. Geometry is frozen for one interval and refreshed from the
    /// accepted start state on the next call to step_electrical/step_heat_network.
    pub fn with_deforming_conductors(
        mut self,
        paths: Vec<FrameConductorPath>,
        gate: &CancelGate,
    ) -> Result<Self> {
        check(gate)?;
        let Some(binding) = self.electrical.as_ref() else {
            return Err(invalid("deforming conductors require an existing material circuit"));
        };
        if paths.is_empty() || paths.len() > binding.circuit.resistors().len()
            || !binding.conductor_paths.is_empty() {
            return Err(invalid("deforming conductors require one bounded, previously unbound path set"));
        }
        let mut seen = BTreeSet::new();
        let mut bound = Vec::with_capacity(paths.len());
        for path in paths {
            check(gate)?;
            if !seen.insert(path.branch) {
                return Err(invalid("deforming conductor branches must be unique"));
            }
            let member = self.members.get(path.member)
                .ok_or_else(|| invalid("deforming conductor names no frame member"))?;
            let state = self.states.get(path.member)
                .ok_or_else(|| invalid("deforming conductor member has no material state"))?;
            let layer = state.source.section.section().layers().get(path.layer)
                .ok_or_else(|| invalid("deforming conductor names no concentric layer"))?;
            let integral = state.source.section.section().layer_integrals().get(path.layer)
                .ok_or_else(|| invalid("deforming conductor layer has no reference geometry"))?;
            let resistor = binding.circuit.resistors().iter().find(|r| r.branch == path.branch)
                .ok_or_else(|| invalid("deforming conductor names no material resistor"))?;
            let destination = binding.destinations.get(&path.branch)
                .ok_or_else(|| invalid("deforming conductor lacks a caloric destination"))?;
            if *destination != state.layers[path.layer].destination.cell
                || resistor.conductor.material().card_identity() != layer.material.card_identity()
                || resistor.conductor.material().query() != layer.material.query() {
                return Err(invalid("deforming conductor must share its frame layer material and caloric cell"));
            }
            bound.push(BoundConductorPath {
                branch: path.branch, member: path.member, layer: path.layer,
                component: member.component, axial_column: member.columns[0].clone(),
                reference_length_m: member.region.length_m,
                reference_area_m2: integral.area_m2(),
                linear_strain_limit: member.linear_strain_limit,
            });
        }
        let samples: Vec<_> = bound.iter().map(|path|
            sample(path, &self.states, self.mechanical.components())).collect::<Result<_>>()?;
        for geometry in &samples {
            let resistor = binding.circuit.resistors().iter().find(|r| r.branch == geometry.branch)
                .expect("validated branch");
            if !close(resistor.conductor.length_m(), geometry.length_m)
                || !close(resistor.conductor.area_m2(), geometry.area_m2) {
                return Err(invalid("initial material circuit geometry must match the accepted frame layer"));
            }
        }
        let binding = self.electrical.as_mut().expect("validated material circuit");
        binding.conductor_paths = bound;
        binding.conductor_geometry = samples;
        check(gate)?;
        Ok(self)
    }

    /// Geometry used by the last attempted accepted electrical interval.
    pub fn electrical_conductor_geometry(&self) -> &[FrameConductorGeometry] {
        self.electrical.as_ref().map_or(&[], |binding| binding.conductor_geometry.as_slice())
    }

    pub(in super::super) fn refresh_conductor_geometry(
        &mut self, gate: &CancelGate,
    ) -> Result<Option<Vec<FrameConductorGeometry>>> {
        check(gate)?;
        let Some(binding) = self.electrical.as_ref() else { return Ok(None); };
        let previous = binding.conductor_geometry.clone();
        if binding.conductor_paths.is_empty() { return Ok(Some(previous)); }
        let next: Vec<_> = binding.conductor_paths.iter().map(|path|
            sample(path, &self.states, self.mechanical.components())).collect::<Result<_>>()?;
        self.electrical.as_mut().expect("binding exists").conductor_geometry = next;
        check(gate)?;
        Ok(Some(previous))
    }

    pub(in super::super) fn restore_conductor_geometry(
        &mut self, previous: Option<Vec<FrameConductorGeometry>>,
    ) {
        if let (Some(binding), Some(previous)) = (&mut self.electrical, previous) {
            binding.conductor_geometry = previous;
        }
    }
}
