//! Source-resolved electrical contacts on matching original exterior traces.
//!
//! The surface term integrates (V_a-V_b)(test_a-test_b)/R'' with the same P1
//! trace basis as the bulk FEEC operator. Duplicate, coincident vertices retain
//! the physical voltage jump; a small fictitious conducting layer is unnecessary.

use super::*;
use crate::planar::{MagneticError, PlanarMagnetostatics};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ResolvedElectricalInterface,
};
use fs_qty::{Dims, QuantitySpec};
use fs_sparse::Coo;

/// One ordered, zero-thickness material interface. The two original exterior
/// edges must have distinct vertices with matching endpoint coordinates, in
/// this exact pairing order, and their adjacent cells must lie on opposite sides.
/// Source surface A belongs to `first_vertices`; B belongs to `second_vertices`.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionContact {
    pub first_vertices: [usize; 2],
    pub second_vertices: [usize; 2],
    pub source: ResolvedElectricalInterface,
}

/// Physical contact response, oriented from the first surface to the second.
/// Current density and voltage jump vary linearly along the paired edge;
/// current and irreversible heat are their exact surface integrals.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionContactResponse {
    pub voltage_jump_v: [f64; 2],
    pub current_density_a_m2: [f64; 2],
    pub current_a: f64,
    pub power_w: f64,
}

fn edge_key(mut vertices: [usize; 2]) -> [usize; 2] {
    vertices.sort_unstable();
    vertices
}

impl<C: PlanarConductorMaterial> PlanarConduction<C> {
    /// Assemble bulk and source-backed surface conduction in one native field.
    /// The specific contact resistance R'' [ohm m2] is divided by actual area,
    /// with the exact P1 surface mass matrix, not a lumped endpoint resistor.
    /// Material, pressure or temperature changes require resolving the original
    /// interface source again and preparing a new operator.
    ///
    /// Finite positive, static R'' is required. Explicit zero resistance belongs
    /// on a conforming shared-node mesh; this owner does not invent a penalty.
    /// No surface capacitance, moving contact, or thermal heat partition is inferred.
    pub fn new_with_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if contacts.len() > limits.max_triplets / 16 {
            return Err(ConductionError::Invalid(
                "electrical contacts exceed the declared triplet capacity",
            ));
        }
        let axes = vec![[None, None]; contacts.len()];
        Self::new_with_phase_contacts(mesh, conductors, contacts, axes, limits, cancelled)
    }

    /// The SAME native surface operator with explicit externally owned liquid
    /// mass fractions on either ordered side. Each declared coordinate must be
    /// present as a dimensional scalar in [0,1]. Every other source response
    /// and applicability coordinate retains the fixed-external admission rule.
    ///
    /// A thermal/phase owner must resolve these exact values from its accepted
    /// caloric state before assembly, then refresh the same sources on every
    /// proposed change. This method neither determines phase from temperature
    /// nor adds phase mixing, motion, storage or an electrical source policy.
    /// Axis rows correspond to INPUT contacts and are sorted with their traces;
    /// returned axes use the same canonical order as `contacts()`.
    pub fn new_with_phase_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        Self::assemble_with_phase_contacts(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            |source| source.has_admitted_static_sources(),
        )
    }

    // Constitutive response owners admit their own bulk coordinates; geometry,
    // interface support, surface quadrature and current/heat stay native here.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn assemble_with_phase_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
        admit: impl Fn(&C) -> bool,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if liquid_mass_fraction_axes.len() != contacts.len() {
            return Err(ConductionError::Invalid(
                "one phase-coordinate row per electrical contact required",
            ));
        }
        let triplets = mesh.faces().len().checked_mul(9).and_then(|bulk| {
            contacts
                .len()
                .checked_mul(16)
                .and_then(|surface| bulk.checked_add(surface))
        });
        if triplets.is_none_or(|count| count > limits.max_triplets) {
            return Err(ConductionError::Invalid(
                "bulk and contact assembly exceed the declared triplet capacity",
            ));
        }
        let mut problem = Self::assemble(mesh, conductors, limits, cancelled, admit)?;
        if contacts.is_empty() {
            return Ok(problem);
        }
        let trace = problem.mesh.boundary_trace().map_err(|_| {
            ConductionError::Invalid("electrical contacts require original exterior traces")
        })?;
        let exterior: BTreeMap<_, _> = trace
            .edges
            .iter()
            .map(|edge| (edge_key(edge.oriented_vertices.map(|v| v as usize)), edge))
            .collect();
        let mut bound: Vec<_> = contacts
            .into_iter()
            .zip(liquid_mass_fraction_axes)
            .collect();
        bound.sort_by_key(|(contact, _)| {
            (
                edge_key(contact.first_vertices),
                edge_key(contact.second_vertices),
            )
        });
        let (contacts, liquid_mass_fraction_axes): (Vec<_>, Vec<_>) = bound.into_iter().unzip();
        let n = problem.mesh.vertices().len();
        let mut matrix = Coo::new(n, n);
        for row in 0..n {
            poll(cancelled)?;
            let (columns, values) = problem.stiffness.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                matrix.push(row, column, value);
            }
        }
        for (contact, axes) in contacts.iter().zip(&liquid_mass_fraction_axes) {
            poll(cancelled)?;
            let first = edge_key(contact.first_vertices);
            let second = edge_key(contact.second_vertices);
            let Some(first_edge) = exterior.get(&first) else {
                return Err(ConductionError::Invalid(
                    "first contact trace is not an exterior edge",
                ));
            };
            let Some(second_edge) = exterior.get(&second) else {
                return Err(ConductionError::Invalid(
                    "second contact trace is not an exterior edge",
                ));
            };
            if first.iter().any(|vertex| second.contains(vertex))
                || !problem.contact_edges.insert(first)
                || !problem.contact_edges.insert(second)
            {
                return Err(ConductionError::Invalid(
                    "contact traces require distinct original vertices and one interface per edge",
                ));
            }
            for endpoint in 0..2 {
                if problem.mesh.vertices()[contact.first_vertices[endpoint]]
                    != problem.mesh.vertices()[contact.second_vertices[endpoint]]
                {
                    return Err(ConductionError::Invalid(
                        "contact endpoint pairs must have identical physical coordinates",
                    ));
                }
            }
            let [a, b] = contact.first_vertices.map(|v| problem.mesh.vertices()[v]);
            let side = |face: usize| -> Result<f64, ConductionError> {
                let vertex = problem.mesh.faces()[face]
                    .iter()
                    .copied()
                    .find(|&v| {
                        let p = problem.mesh.vertices()[v as usize];
                        p != a && p != b
                    })
                    .ok_or(ConductionError::Invalid(
                        "contact has no adjacent cell interior",
                    ))?;
                let p = problem.mesh.vertices()[vertex as usize];
                finite((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]))
            };
            let first_side = side(first_edge.source_face)?;
            let second_side = side(second_edge.source_face)?;
            if first_side == 0.0
                || second_side == 0.0
                || first_side.is_sign_positive() == second_side.is_sign_positive()
            {
                return Err(ConductionError::Invalid(
                    "contact cells must lie on opposite sides of the matched trace",
                ));
            }
            contact
                .source
                .require_surface_materials(
                    problem.conductors[first_edge.source_face]
                        .material()
                        .material(),
                    problem.conductors[second_edge.source_face]
                        .material()
                        .material(),
                )
                .map_err(|_| {
                    ConductionError::Invalid(
                        "contact source must match both ordered adjacent material states",
                    )
                })?;
            let mut owned: Vec<(&str, QuantitySpec)> = Vec::new();
            for axis in axes.iter().flatten() {
                if axis.trim().is_empty()
                    || owned.iter().any(|(previous, _)| *previous == axis.as_str())
                    || contact.source.resolved().axis_quantities().get(axis)
                        != Some(&QuantitySpec::dimensional(Dims::NONE))
                    || !contact
                        .source
                        .resolved()
                        .state_coordinate(axis)
                        .is_some_and(|value| value.is_finite() && (0.0..=1.0).contains(&value))
                {
                    return Err(ConductionError::Invalid(
                        "owned contact phase needs a unique existing dimensional mass fraction in [0,1]",
                    ));
                }
                owned.push((axis.as_str(), QuantitySpec::dimensional(Dims::NONE)));
            }
            if !contact
                .source
                .resolved()
                .property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
                .is_some_and(|property| {
                    if owned.is_empty() {
                        property.is_static_at_fixed_external_conditions()
                    } else {
                        property.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
                    }
                })
            {
                return Err(ConductionError::Invalid(
                    "contact resistance requires static source support at fixed external conditions",
                ));
            }
            let resistance = contact.source.resistance_ohm_m2();
            if !positive(resistance) {
                return Err(ConductionError::Invalid(
                    "positive specific contact resistance required; perfect contact needs shared mesh vertices",
                ));
            }
            let area = problem.mesh.edge_measure(first_edge.global_edge).ok_or(
                ConductionError::Invalid("contact trace lacks physical area"),
            )?;
            let conductance = finite(area / resistance)?;
            let sixth = conductance / 6.0;
            if !positive(sixth) {
                return Err(ConductionError::Invalid(
                    "contact conductance is not representable",
                ));
            }
            for i in 0..2 {
                for j in 0..2 {
                    let value = sixth * if i == j { 2.0 } else { 1.0 };
                    let ai = contact.first_vertices[i];
                    let aj = contact.first_vertices[j];
                    let bi = contact.second_vertices[i];
                    let bj = contact.second_vertices[j];
                    matrix.push(ai, aj, value);
                    matrix.push(bi, bj, value);
                    matrix.push(ai, bj, -value);
                    matrix.push(bi, aj, -value);
                }
            }
            problem.contact_conductance_s.push(conductance);
        }
        poll(cancelled)?;
        problem.stiffness = matrix.assemble();
        for row in 0..n {
            poll(cancelled)?;
            if problem
                .stiffness
                .row(row)
                .1
                .iter()
                .any(|value| !value.is_finite())
            {
                return Err(ConductionError::Invalid(
                    "unrepresentable bulk/contact stiffness",
                ));
            }
        }
        problem.contacts = contacts;
        problem.contact_liquid_mass_fraction_axes = liquid_mass_fraction_axes;
        poll(cancelled)?;
        Ok(problem)
    }
}

impl<C: PlanarOhmicResponse> PlanarConduction<C> {
    /// Original ordered source bindings in canonical first/second edge order.
    pub fn contacts(&self) -> &[ConductionContact] {
        &self.contacts
    }

    /// Explicit phase-owner coordinates in the same canonical order as contacts.
    /// The fixed-external constructor records two `None` entries per contact.
    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        &self.contact_liquid_mass_fraction_axes
    }

    pub(super) fn is_contact_edge(&self, vertices: [usize; 2]) -> bool {
        self.contact_edges.contains(&edge_key(vertices))
    }

    pub(super) fn prescribed_boundary(
        &self,
        boundary: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>, ConductionError> {
        if self.contacts.is_empty() {
            return Ok(scalar_potential::prescribed_boundary(
                &self.mesh, boundary, cancelled,
            )?);
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_traces(&self.mesh, boundary, cancelled)
            .map_err(|error| match error {
            MagneticError::Cancelled => ConductionError::Cancelled,
            MagneticError::Invalid(reason) => ConductionError::Invalid(reason),
            _ => ConductionError::Invalid("invalid conductor boundary trace"),
        })?;
        for &[a, b] in &self.contact_edges {
            poll(cancelled)?;
            if prescribed[a].is_some() && prescribed[b].is_some() {
                return Err(ConductionError::Invalid(
                    "internal contact traces cannot also be voltage electrodes",
                ));
            }
        }
        let n = prescribed.len();
        let mut neighbors = vec![Vec::new(); n];
        for &[a, b] in self.mesh.edges() {
            poll(cancelled)?;
            neighbors[a as usize].push(b as usize);
            neighbors[b as usize].push(a as usize);
        }
        for contact in &self.contacts {
            poll(cancelled)?;
            for endpoint in 0..2 {
                let a = contact.first_vertices[endpoint];
                let b = contact.second_vertices[endpoint];
                neighbors[a].push(b);
                neighbors[b].push(a);
            }
        }
        let mut visited = vec![false; n];
        let mut stack = Vec::new();
        for start in 0..n {
            if visited[start] {
                continue;
            }
            visited[start] = true;
            stack.push(start);
            let mut anchored = false;
            while let Some(vertex) = stack.pop() {
                poll(cancelled)?;
                anchored |= prescribed[vertex].is_some();
                for &next in &neighbors[vertex] {
                    if !visited[next] {
                        visited[next] = true;
                        stack.push(next);
                    }
                }
            }
            if !anchored {
                return Err(ConductionError::Invalid(
                    "every contact-connected component needs a prescribed exterior edge",
                ));
            }
        }
        Ok(prescribed)
    }

    pub(super) fn contact_response(
        &self,
        voltage: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<ConductionContactResponse>, ConductionError> {
        let mut response = Vec::with_capacity(self.contacts.len());
        for (contact, &conductance) in self.contacts.iter().zip(&self.contact_conductance_s) {
            poll(cancelled)?;
            let [a0, a1] = contact.first_vertices;
            let [b0, b1] = contact.second_vertices;
            let d0 = finite(voltage[a0] - voltage[b0])?;
            let d1 = finite(voltage[a1] - voltage[b1])?;
            let resistance = contact.source.resistance_ohm_m2();
            let sum = finite(d0 + d1)?;
            let power = finite((conductance / 6.0) * (d0 * d0 + d1 * d1 + sum * sum))?;
            if power < 0.0 || ((d0 != 0.0 || d1 != 0.0) && power == 0.0) {
                return Err(ConductionError::Invalid(
                    "contact Joule heat is not representably positive",
                ));
            }
            response.push(ConductionContactResponse {
                voltage_jump_v: [d0, d1],
                current_density_a_m2: [finite(d0 / resistance)?, finite(d1 / resistance)?],
                current_a: finite(conductance * (0.5 * d0 + 0.5 * d1))?,
                power_w: power,
            });
        }
        Ok(response)
    }
}
