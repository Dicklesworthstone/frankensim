//! Source-bound voltage jumps across matching original triangular traces.
//! The P1 contact term is integral([test][V]/R'' dA); its exact surface mass
//! matrix is A/12 times the matrix with diagonal 2 and off-diagonal 1.
//! Bulk Hodge, boundary/gauge preparation and sparse solve retain their owners.

use super::*;
use crate::conduction::{ConductionError, admit_contact_source};
use fs_material::electrical_interface::ResolvedElectricalInterface;
use fs_sparse::Coo;
use std::collections::BTreeMap;

/// Ordered zero-thickness interface between two matching exterior triangles.
/// Vertices are distinct across sides and pair by their EXACT coordinates in
/// the supplied order. Opposite adjacent volumes keep the physical voltage jump.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialConductionContact {
    pub first_vertices: [usize; 3],
    pub second_vertices: [usize; 3],
    pub source: ResolvedElectricalInterface,
}

/// P1 trace response oriented from source surface A to source surface B.
/// Current and heat are exact triangle integrals, not lumped nodal estimates.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialConductionContactResponse {
    pub voltage_jump_v: [f64; 3],
    pub current_density_a_m2: [f64; 3],
    pub current_a: f64,
    pub power_w: f64,
}

fn key(mut vertices: [usize; 3]) -> [usize; 3] {
    vertices.sort_unstable();
    vertices
}

impl<C: SpatialConductorMaterial> SpatialConduction<C> {
    /// Prepare finite source-backed contact resistance on the original 3D
    /// field. R'' has units ohm m2; actual triangle area determines conductance.
    /// Explicit zero resistance requires shared mesh vertices, not a penalty.
    /// No opening, spreading, arcing, interface storage or heat split is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_contacts(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        conductors: Vec<C>,
        contacts: Vec<SpatialConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if contacts.len() > limits.max_triplets / 36 {
            return Err(SpatialConductionError::Invalid(
                "3D contacts exceed the assembly budget",
            ));
        }
        let axes = vec![[None, None]; contacts.len()];
        Self::new_with_phase_contacts(
            positions_m,
            tetrahedra,
            conductors,
            contacts,
            axes,
            limits,
            cancelled,
        )
    }

    /// The same native surface operator with explicitly owned liquid MASS
    /// fractions on either ordered side. Axis rows follow input contacts and
    /// are canonicalized with them. A thermal owner must refresh these source
    /// coordinates from its caloric states; no mixing law follows from a name.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_phase_contacts(
        positions_m: Vec<[f64; 3]>,
        tetrahedra: Vec<[u32; 4]>,
        conductors: Vec<C>,
        contacts: Vec<SpatialConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if contacts.len() > limits.max_triplets / 36 {
            return Err(SpatialConductionError::Invalid(
                "3D contacts exceed the assembly budget",
            ));
        }
        Self::new(positions_m, tetrahedra, conductors, limits, cancelled)?.bind_contacts(
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
        )
    }

    /// Refresh original bulk and interface sources on unchanged matching traces.
    /// Contacts must remain in `contacts()` order with the same paired vertices;
    /// phase ownership is retained. Thermal history and source pinning belong
    /// to the coupling owner, not to this immutable field preparation.
    pub fn with_conductors_and_contacts(
        &self,
        conductors: Vec<C>,
        contacts: Vec<SpatialConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        if contacts.len() != self.contacts.len()
            || contacts.iter().zip(&self.contacts).any(|(a, b)| {
                a.first_vertices != b.first_vertices || a.second_vertices != b.second_vertices
            })
        {
            return Err(SpatialConductionError::Invalid(
                "3D contact refresh must preserve original paired traces",
            ));
        }
        Self::assemble(Arc::clone(&self.geometry), conductors, limits, cancelled)?.bind_contacts(
            contacts,
            self.contact_liquid_mass_fraction_axes.clone(),
            limits,
            cancelled,
        )
    }

    pub(super) fn bind_contacts(
        mut self,
        contacts: Vec<SpatialConductionContact>,
        axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        check(cancelled)?;
        let triplets = self
            .geometry
            .mesh
            .tets
            .len()
            .checked_mul(36)
            .and_then(|bulk| {
                contacts
                    .len()
                    .checked_mul(36)
                    .and_then(|surface| bulk.checked_add(surface))
            });
        if axes.len() != contacts.len() || triplets.is_none_or(|n| n > limits.max_triplets) {
            return Err(SpatialConductionError::Invalid(
                "complete contact phase roles and bounded bulk/surface assembly required",
            ));
        }
        if contacts.is_empty() {
            return Ok(self);
        }
        let g = &self.geometry;
        let faces: BTreeMap<_, _> = g
            .mesh
            .faces
            .iter()
            .enumerate()
            .filter(|(face, _)| g.faces[*face].cells.len() == 1)
            .map(|(face, vertices)| (vertices.map(|v| v as usize), face))
            .collect();
        let mut bound: Vec<_> = contacts.into_iter().zip(axes).collect();
        bound.sort_by_key(|(contact, _)| {
            (key(contact.first_vertices), key(contact.second_vertices))
        });
        let (contacts, axes): (Vec<_>, Vec<_>) = bound.into_iter().unzip();
        let n = g.mesh.vertex_count;
        let mut matrix = Coo::new(n, n);
        for row in 0..n {
            check(cancelled)?;
            let (columns, values) = self.stiffness.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                matrix.push(row, column, value);
            }
        }
        let mut parent = self.components.clone();
        for (contact, axes) in contacts.iter().zip(&axes) {
            check(cancelled)?;
            let first = key(contact.first_vertices);
            let second = key(contact.second_vertices);
            let a = *faces.get(&first).ok_or(SpatialConductionError::Invalid(
                "first contact trace is not an original exterior triangle",
            ))?;
            let b = *faces.get(&second).ok_or(SpatialConductionError::Invalid(
                "second contact trace is not an original exterior triangle",
            ))?;
            if first.iter().any(|v| second.contains(v))
                || !self.contact_faces.insert(a)
                || !self.contact_faces.insert(b)
            {
                return Err(SpatialConductionError::Invalid(
                    "3D contacts need distinct vertices and exactly one interface per triangle",
                ));
            }
            for i in 0..3 {
                if g.positions_m[contact.first_vertices[i]]
                    != g.positions_m[contact.second_vertices[i]]
                {
                    return Err(SpatialConductionError::Invalid(
                        "paired contact vertices must have identical physical coordinates",
                    ));
                }
            }
            if dot(g.faces[a].outward_normal, g.faces[b].outward_normal) > -1.0 + 1e-12 {
                return Err(SpatialConductionError::Invalid(
                    "contact tetrahedra must lie on opposite sides of the matched face",
                ));
            }
            let cells = [g.faces[a].cells[0], g.faces[b].cells[0]];
            contact
                .source
                .require_surface_materials(
                    self.conductors[cells[0]].material().material(),
                    self.conductors[cells[1]].material().material(),
                )
                .map_err(|_| {
                    SpatialConductionError::Invalid(
                        "3D contact source must match the ordered adjacent materials",
                    )
                })?;
            let resistance =
                admit_contact_source(&contact.source, axes).map_err(|error| match error {
                    ConductionError::Invalid(message) => SpatialConductionError::Invalid(message),
                    _ => SpatialConductionError::Invalid("3D contact source admission refused"),
                })?;
            let conductance = finite(g.faces[a].area_m2 / resistance)?;
            let twelfth = conductance / 12.0;
            if !twelfth.is_finite() || twelfth <= 0.0 {
                return Err(SpatialConductionError::Invalid(
                    "3D contact conductance is not representably positive",
                ));
            }
            for i in 0..3 {
                for j in 0..3 {
                    let value = twelfth * if i == j { 2.0 } else { 1.0 };
                    let (ai, aj) = (contact.first_vertices[i], contact.first_vertices[j]);
                    let (bi, bj) = (contact.second_vertices[i], contact.second_vertices[j]);
                    matrix.push(ai, aj, value);
                    matrix.push(bi, bj, value);
                    matrix.push(ai, bj, -value);
                    matrix.push(bi, aj, -value);
                }
                let a = super::super::root(&mut parent, contact.first_vertices[i]);
                let b = super::super::root(&mut parent, contact.second_vertices[i]);
                if a != b {
                    parent[a.max(b)] = a.min(b);
                }
            }
            self.contact_cells.push(cells);
            self.contact_conductance_s.push(conductance);
        }
        check(cancelled)?;
        self.stiffness = matrix.assemble();
        for row in 0..n {
            check(cancelled)?;
            for &value in self.stiffness.row(row).1 {
                finite(value)?;
            }
        }
        self.components = (0..n).map(|v| super::super::root(&mut parent, v)).collect();
        self.contacts = contacts;
        self.contact_liquid_mass_fraction_axes = axes;
        check(cancelled)?;
        Ok(self)
    }

    /// Original source bindings in canonical first/second triangular order.
    pub fn contacts(&self) -> &[SpatialConductionContact] {
        &self.contacts
    }
    /// Original adjacent tetrahedra in contact order and source A/B order.
    pub fn contact_cells(&self) -> &[[usize; 2]] {
        &self.contact_cells
    }
    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        &self.contact_liquid_mass_fraction_axes
    }

    pub(super) fn contact_response(
        &self,
        voltage: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<SpatialConductionContactResponse>> {
        let mut result = Vec::with_capacity(self.contacts.len());
        for (contact, &conductance) in self.contacts.iter().zip(&self.contact_conductance_s) {
            check(cancelled)?;
            let jump = [0, 1, 2].map(|i| {
                finite(voltage[contact.first_vertices[i]] - voltage[contact.second_vertices[i]])
            });
            let [a, b, c] = jump;
            let jump = [a?, b?, c?];
            let sum = finite(jump.iter().sum())?;
            let power = finite(
                (conductance / 12.0) * (jump.iter().map(|v| v * v).sum::<f64>() + sum * sum),
            )?;
            if power < 0.0 || (jump.iter().any(|&v| v != 0.0) && power == 0.0) {
                return Err(SpatialConductionError::Invalid(
                    "3D contact Joule heat is not representably positive",
                ));
            }
            let density = jump.map(|v| finite(v / contact.source.resistance_ohm_m2()));
            let [a, b, c] = density;
            result.push(SpatialConductionContactResponse {
                voltage_jump_v: jump,
                current_density_a_m2: [a?, b?, c?],
                current_a: finite((conductance / 3.0) * sum)?,
                power_w: power,
            });
        }
        Ok(result)
    }
}
