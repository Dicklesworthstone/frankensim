//! Regional conductivity and named contact controls of the accepted trajectory.
//!
//! Enthalpy conductivity bars already include dt; its contact dual is the
//! shared masked nodal-load pullback, with dt included exactly once. Capacity
//! storage also supplies watt-coordinate bars. Both use the SAME complete dual
//! and original material/contact owners. Neither path resolves a primal or
//! another transpose.

use std::collections::BTreeMap;
use std::time::Instant;

use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ElementMaterials, ThermalInterfaces,
};
use fs_exec::Cx;
use fs_project::ProjectSpec;

use super::{SolveRefusal, bad, finite, json_string, lower, number, poll, zeros};

struct ConductivityRegion {
    region: String,
    material_card: String,
    material_state: String,
    elements: usize,
    derivative: f64,
}

struct ContactInterface {
    interface: String,
    interface_card: String,
    mapped: bool,
    face_pairs: usize,
    derivative: f64,
}

/// Extra accepted-step conductive terms (32 bytes/element), returned element
/// bars (8), regional owner indices, and one nodal contact-load dual. Metadata,
/// checked aggregation and output rows have a bounded allowance per declaration.
/// The estimate describes logical storage, not allocator overhead or peak RSS.
pub(in crate::solve::transient) fn storage_bytes(
    n: usize,
    elements: usize,
    regions: usize,
    contacts: usize,
) -> Result<usize, SolveRefusal> {
    if n == 0 || elements == 0 || regions > 64 || contacts > 64 {
        return Err(bad(
            "enthalpy material controls need a nonempty mesh and at most 64 regions and contacts",
        ));
    }
    let per_element = 5 * std::mem::size_of::<f64>() + std::mem::size_of::<usize>();
    elements.checked_mul(per_element)
        .and_then(|v| n.checked_mul(std::mem::size_of::<f64>())
            .and_then(|b| v.checked_add(b)))
        .and_then(|v| regions.checked_add(contacts)
            .and_then(|rows| rows.checked_mul(8192))
            .and_then(|b| v.checked_add(b)))
        .filter(|&v| isize::try_from(v).is_ok())
        .ok_or_else(|| bad("enthalpy material-control storage allowance overflows"))
}

/// Immutable material/contact bindings and controls held throughout the
/// selected prefix. Sharing a card never merges distinct declared coordinates.
pub(in crate::solve::transient::adjoint) struct Controls<'m> {
    mesh: &'m ConductionMesh,
    fallback: &'m ConductivityModel,
    materials: Option<&'m ElementMaterials>,
    interfaces: Option<&'m ThermalInterfaces>,
    element_regions: Vec<usize>,
    conductivity: Vec<ConductivityRegion>,
    contacts: Vec<ContactInterface>,
    steps: usize,
}

impl<'m> Controls<'m> {
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::solve::transient::adjoint) fn new(
        cx: &Cx<'_>,
        spec: &ProjectSpec,
        problem: ConductionProblem<'m>,
        interfaces: Option<&'m ThermalInterfaces>,
        labels: &[u32],
        ids: &BTreeMap<String, u32>,
        initial_temperature: &[f64],
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        poll(cx, deadline)?;
        let bindings = spec.materials.as_deref().unwrap_or(&[]);
        let declared_contacts = spec.interface_cards.as_deref().unwrap_or(&[]);
        let n = problem.mesh.vertex_count();
        let elements = problem.mesh.element_count();
        storage_bytes(n, elements, bindings.len(), declared_contacts.len())?;
        if labels.len() != elements || initial_temperature.len() != n
            || problem.element_materials.is_some_and(|value| value.of_element().len() != elements)
            || interfaces.map_or(0, ThermalInterfaces::surface_count) != declared_contacts.len()
        {
            return Err(bad(
                "enthalpy material controls require the original element and contact layout",
            ));
        }
        for (vertex, &temperature) in initial_temperature.iter().enumerate() {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(bad("enthalpy material-control initial temperatures must be positive"));
            }
        }
        let mut region_owners = BTreeMap::new();
        let mut conductivity = Vec::new();
        conductivity.try_reserve_exact(bindings.len())
            .map_err(|_| bad("enthalpy conductivity-control metadata allocation refused"))?;
        for (index, binding) in bindings.iter().enumerate() {
            poll(cx, deadline)?;
            let id = *ids.get(&binding.region)
                .ok_or_else(|| bad("enthalpy conductivity control has no native region label"))?;
            if region_owners.insert(id, index).is_some() {
                return Err(bad("enthalpy conductivity controls need one binding per native region"));
            }
            conductivity.push(ConductivityRegion {
                region: binding.region.clone(),
                material_card: binding.card.clone(),
                material_state: binding.state.clone(),
                elements: 0,
                derivative: 0.0,
            });
        }
        let mut element_regions = Vec::new();
        if !bindings.is_empty() {
            element_regions.try_reserve_exact(elements)
                .map_err(|_| bad("enthalpy conductivity owner-map allocation refused"))?;
            for (element, id) in labels.iter().enumerate() {
                if element % 256 == 0 { poll(cx, deadline)?; }
                let owner = *region_owners.get(id)
                    .ok_or_else(|| bad("a retained enthalpy element has no declared material control"))?;
                element_regions.push(owner);
                conductivity[owner].elements = conductivity[owner].elements.checked_add(1)
                    .ok_or_else(|| bad("enthalpy conductivity-control element count overflows"))?;
            }
            if conductivity.iter().any(|row| row.elements == 0) {
                return Err(bad("declared conductivity control has no retained volume elements"));
            }
        }
        let mut contacts = Vec::new();
        contacts.try_reserve_exact(declared_contacts.len())
            .map_err(|_| bad("enthalpy contact-control metadata allocation refused"))?;
        if let Some(interfaces) = interfaces {
            // A constant objective has an exactly zero dual. The physical
            // initial temperature retains both contact traces while the shared
            // owner supplies original names, cards, mapping and face counts.
            let zero_dual = zeros(n)?;
            let gradients = interfaces.matching_resistance_scale_pullback(
                cx, initial_temperature, &zero_dual, problem.mesh.boundary().len(),
            ).map_err(lower)?;
            let mut by_name = BTreeMap::new();
            for gradient in gradients {
                poll(cx, deadline)?;
                if by_name.insert(gradient.interface.clone(), gradient).is_some() {
                    return Err(bad("enthalpy contact-control names are duplicated"));
                }
            }
            for binding in declared_contacts {
                poll(cx, deadline)?;
                let gradient = by_name.remove(&binding.interface)
                    .ok_or_else(|| bad("declared enthalpy contact control has no unique bound identity"))?;
                let card = gradient.card_identity.to_hex();
                if card != binding.card || gradient.face_pairs == 0 || gradient.derivative != 0.0 {
                    return Err(bad("enthalpy contact-control card, paired faces or zero-dual binding differs"));
                }
                contacts.push(ContactInterface {
                    interface: gradient.interface,
                    interface_card: card,
                    mapped: gradient.mapped,
                    face_pairs: gradient.face_pairs,
                    derivative: 0.0,
                });
            }
            if !by_name.is_empty() {
                return Err(bad("a bound enthalpy contact control has no declared owner"));
            }
        }
        poll(cx, deadline)?;
        Ok(Self {
            mesh: problem.mesh,
            fallback: problem.material,
            materials: problem.element_materials,
            interfaces,
            element_regions,
            conductivity,
            contacts,
            steps: 0,
        })
    }

    /// Aggregate one accepted endpoint after its complete multiplier passes the
    /// original true-transpose gate. Conductivity already includes dt; contact
    /// receives the masked watt-coordinate dual. Source and boundary values
    /// may vary by step.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn add(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        temperature: &[f64],
        element_conductivity: &[f64],
        nodal_load_adjoint: &[f64],
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        self.add_scaled(
            cx, problem, interfaces, temperature, element_conductivity,
            nodal_load_adjoint, 1.0, deadline,
        )
    }

    /// Capacity storage already exposes dJ/d(nodal load in W). Contract the
    /// production K(T) at the unchanged endpoint, including anisotropy and
    /// prescribed-temperature lifts. K'(T) is in the complete state tangent,
    /// not another factor in this fixed-state material partial.
    ///
    /// The caller must first gate the complete transpose, including all air,
    /// natural and radiative feedback. No new primal or transpose is solved.
    /// A scale of one below is the identity change of residual coordinates,
    /// NOT a one-second physical timestep. Neither control receives another dt.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::solve::transient::adjoint) fn add_capacity(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        temperature: &[f64],
        nodal_load_adjoint: &[f64],
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        let conductivity = fs_conduction::adjoint::RobinResponse::conductivity_scale_pullback(
            cx, problem, temperature, nodal_load_adjoint, self.mesh.element_count(),
        ).map_err(lower)?;
        self.add_scaled(
            cx, problem, interfaces, temperature, &conductivity,
            nodal_load_adjoint, 1.0, deadline,
        )
    }

    /// Element bars are already objective derivatives. Contact uses the
    /// supplied masked nodal-load dual with its explicit coordinate scale.
    /// Both storage formulations retain the same owner checks and aggregation.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn add_scaled(
        &mut self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        temperature: &[f64],
        element_conductivity: &[f64],
        full_residual_multiplier: &[f64],
        load_scale: f64,
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        let same_materials = match (self.materials, problem.element_materials) {
            (Some(original), Some(current)) => std::ptr::eq(original, current),
            (None, None) => true,
            _ => false,
        };
        let same_interfaces = match (self.interfaces, interfaces) {
            (Some(original), Some(current)) => std::ptr::eq(original, current),
            (None, None) => true,
            _ => false,
        };
        let n = self.mesh.vertex_count();
        if !std::ptr::eq(self.mesh, problem.mesh)
            || !std::ptr::eq(self.fallback, problem.material)
            || !same_materials || !same_interfaces
            || temperature.len() != n || full_residual_multiplier.len() != n
            || element_conductivity.len() != self.mesh.element_count()
            || (!self.conductivity.is_empty()
                && self.element_regions.len() != element_conductivity.len())
            || !(load_scale.is_finite() && load_scale > 0.0)
        {
            return Err(bad("enthalpy material controls lost their accepted physical owner or layout"));
        }
        for (vertex, (&temperature, &lambda)) in temperature.iter()
            .zip(full_residual_multiplier).enumerate()
        {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            if !temperature.is_finite() || temperature <= 0.0 || !lambda.is_finite() {
                return Err(bad("enthalpy material controls require positive temperatures and finite full duals"));
            }
        }
        let mut regional = zeros(self.conductivity.len())?;
        for (element, &derivative) in element_conductivity.iter().enumerate() {
            if element % 256 == 0 { poll(cx, deadline)?; }
            let derivative = finite(derivative)?;
            if !self.conductivity.is_empty() {
                let owner = self.element_regions[element];
                regional[owner] = finite(regional[owner] + derivative)?;
            }
        }
        let mut contact_values = zeros(self.contacts.len())?;
        if let Some(interfaces) = interfaces {
            let mut load_dual = zeros(n)?;
            for (vertex, (weight, &lambda)) in load_dual.iter_mut()
                .zip(full_residual_multiplier).enumerate()
            {
                if vertex % 256 == 0 { poll(cx, deadline)?; }
                *weight = finite(load_scale * lambda)?;
            }
            let gradients = interfaces.matching_resistance_scale_pullback(
                cx, temperature, &load_dual, self.mesh.boundary().len(),
            ).map_err(lower)?;
            if gradients.len() != self.contacts.len() {
                return Err(bad("enthalpy contact controls changed their bound surface count"));
            }
            let mut by_name = BTreeMap::new();
            for gradient in gradients {
                poll(cx, deadline)?;
                if by_name.insert(gradient.interface.clone(), gradient).is_some() {
                    return Err(bad("enthalpy contact-gradient identities are duplicated"));
                }
            }
            for (index, original) in self.contacts.iter().enumerate() {
                poll(cx, deadline)?;
                let gradient = by_name.remove(&original.interface)
                    .ok_or_else(|| bad("enthalpy contact gradient lost a bound interface"))?;
                if gradient.card_identity.to_hex() != original.interface_card
                    || gradient.mapped != original.mapped
                    || gradient.face_pairs != original.face_pairs
                {
                    return Err(bad("enthalpy contact-gradient card, map or paired-face identity changed"));
                }
                contact_values[index] = finite(gradient.derivative)?;
            }
            if !by_name.is_empty() {
                return Err(bad("enthalpy contact gradient has an undeclared interface"));
            }
        }
        for (control, value) in self.conductivity.iter_mut().zip(regional) {
            control.derivative = finite(control.derivative + value)?;
        }
        for (control, value) in self.contacts.iter_mut().zip(contact_values) {
            control.derivative = finite(control.derivative + value)?;
        }
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("enthalpy material-control step count overflows"))?;
        poll(cx, deadline)
    }

    pub(in crate::solve::transient::adjoint) fn fragment(&self, expected_steps: usize) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("enthalpy material controls do not cover the selected reverse prefix"));
        }
        let conductivity = self.conductivity.iter().map(|control| Ok(format!(
            "{{\"region\":{},\"material_card\":{},\"material_state\":{},\"elements\":{},\"reference_value\":1,\"log_reference_value\":0,\"dtemperature_dconductivity_multiplier_k\":{}}}",
            json_string(&control.region), json_string(&control.material_card),
            json_string(&control.material_state), control.elements, number(control.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let contacts = self.contacts.iter().map(|control| Ok(format!(
            "{{\"interface\":{},\"interface_card\":{},\"mapped\":{},\"face_pairs\":{},\"reference_value\":1,\"log_reference_value\":0,\"dtemperature_dresistance_multiplier_k\":{}}}",
            json_string(&control.interface), json_string(&control.interface_card),
            control.mapped, control.face_pairs, number(control.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"material_controls\":{{\"control_steps\":{},\"parameter_scope\":\"constant-through-selected-trajectory-independent-regional-effective-conductivity-and-named-matching-resistance-multipliers\",\"units\":\"K\",\"conductivity_regions\":[{}],\"contact_interfaces\":[{}]}}",
            self.steps, conductivity.join(","), contacts.join(","),
        ))
    }
}
