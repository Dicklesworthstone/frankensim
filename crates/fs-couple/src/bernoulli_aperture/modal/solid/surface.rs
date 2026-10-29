//! Original three-dimensional exterior traces in the shared spatial closure.
use super::super::{ModalAperture, invalid};
use super::{poll, scaled};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::modal::closure::{ModalPlateClosure, displacements};
use crate::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_exec::CancelGate;
use std::collections::{BTreeMap, BTreeSet};

/// Fixed planar exterior patches on one original Tet4 material solid.
/// `contact.lay_triangles` names original `surface_geometry().faces()` indices;
/// clearances cover every original volume node. Only the selected face nodes
/// carry contact. Slit edges name original nodes on the patch boundary.
///
/// Positive gap motion is `opening_sign` times the original outward normal.
/// All selected faces must have the same normal; curved/moving normal contact
/// requires another geometric owner and is refused here. Clearances can vary
/// spatially, including explicit installed interference. This does not change
/// the solid's independently integrated pressure/swept-volume port.
#[derive(Clone, Debug)]
pub struct SolidSurfaceClosureSpec {
    pub contact: PlateClosureSpec,
    pub slit_edges: Vec<[usize; 2]>,
    /// Exactly +1 (outward opening) or -1 (inward opening).
    pub opening_sign: f64,
    /// Upper bound on 6*(nodes+free_dofs)*modes + 16*(faces+slit_edges).
    pub max_projection_terms: usize,
}

impl ModalAperture {
    /// Replace mean-slit contact by the original solid's local surface gaps.
    /// Source face areas, full signed nodal modes and their transposes feed the
    /// existing point-contact potential and coupled spatial pressure solve.
    /// Zero-mean modes, full viscosity and hereditary states remain unchanged.
    /// Geometry is prepared once; no surface construction occurs while playing.
    pub fn with_solid_surface_closure(
        mut self,
        spec: SolidSurfaceClosureSpec,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        if self.samples != 0 || self.closure.is_some() || self.solid_closure.is_some() {
            return Err(invalid(
                "solid surface contact requires one cold attachment",
            ));
        }
        let (surface, modes, free, root_mass) = if let Some(solid) = &self.hereditary_solid {
            let source = solid.source();
            (
                source.surface_geometry(),
                source.modes(),
                source.source().equilibrium().free_dofs.as_slice(),
                solid.coordinate_mass_kg().sqrt(),
            )
        } else if let Some(solid) = &self.material_solid {
            let source = solid.source();
            (
                source.surface_geometry(),
                source.modes(),
                source.source().equilibrium().free_dofs.as_slice(),
                solid.coordinate_mass_kg().sqrt(),
            )
        } else {
            return Err(invalid(
                "surface contact requires its original material solid",
            ));
        };
        let surface = surface.ok_or_else(|| {
            invalid("solid surface contact requires original P1 exterior geometry")
        })?;
        let n = modes.len();
        let nodes = surface.nodes_m();
        let faces = surface.faces();
        let terms = nodes
            .len()
            .checked_add(free.len())
            .and_then(|v| v.checked_mul(n))
            .and_then(|v| v.checked_mul(6))
            .and_then(|v| {
                spec.contact
                    .lay_triangles
                    .len()
                    .checked_add(spec.slit_edges.len())
                    .and_then(|s| s.checked_mul(16))
                    .and_then(|s| v.checked_add(s))
            })
            .ok_or_else(|| invalid("solid surface projection budget overflow"))?;
        if terms > spec.max_projection_terms
            || spec.max_projection_terms == 0
            || spec.opening_sign != 1.0 && spec.opening_sign != -1.0
            || spec.contact.nodal_rest_gap_m.len() != nodes.len()
            || spec.contact.lay_triangles.is_empty()
            || spec.slit_edges.is_empty()
        {
            return Err(invalid(
                "solid surface contact requires complete bounded original geometry",
            ));
        }
        let first = faces
            .get(spec.contact.lay_triangles[0])
            .ok_or_else(|| invalid("solid contact names an absent exterior face"))?;
        let normal = first.outward_normal.map(|v| spec.opening_sign * v);
        let mut weights = vec![0.0; nodes.len()];
        let mut seen = BTreeSet::new();
        let mut incidence = BTreeMap::<[usize; 2], usize>::new();
        for &id in &spec.contact.lay_triangles {
            poll(gate)?;
            let face = faces
                .get(id)
                .ok_or_else(|| invalid("solid contact names an absent exterior face"))?;
            if !seen.insert(id) {
                return Err(invalid("solid contact repeats a physical face"));
            }
            if face
                .outward_normal
                .iter()
                .zip(first.outward_normal)
                .any(|(a, b)| (a - b).abs() > 128.0 * f64::EPSILON)
            {
                return Err(invalid(
                    "solid surface contact requires parallel original outward normals",
                ));
            }
            for &node in &face.nodes {
                weights[node] += face.area_m2 / 3.0;
            }
            for mut edge in [
                [face.nodes[0], face.nodes[1]],
                [face.nodes[1], face.nodes[2]],
                [face.nodes[2], face.nodes[0]],
            ] {
                edge.sort_unstable();
                *incidence.entry(edge).or_default() += 1;
            }
        }
        let mut edges = Vec::with_capacity(spec.slit_edges.len());
        let mut seen_edges = BTreeSet::new();
        for &original in &spec.slit_edges {
            poll(gate)?;
            let mut edge = original;
            edge.sort_unstable();
            if incidence.get(&edge) != Some(&1) || !seen_edges.insert(edge) {
                return Err(invalid(
                    "solid slit edges must be unique original patch boundary edges",
                ));
            }
            let a = nodes[edge[0]];
            let b = nodes[edge[1]];
            let length = (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
            edges.push((original, length));
        }
        let mut shapes = vec![vec![0.0; n]; nodes.len()];
        for (row, &dof) in free.iter().enumerate() {
            poll(gate)?;
            for (mode, shape) in modes.iter().zip(&mut shapes[dof / 3]) {
                *shape += scaled(mode.phi[row], root_mass)? * normal[dof % 3];
            }
        }
        let contact_nodes = weights
            .iter()
            .enumerate()
            .filter(|(_, w)| **w > 0.0)
            .map(|(node, _)| node)
            .collect();
        let closure = ModalPlateClosure::from_traces(
            shapes,
            spec.contact.clone(),
            weights,
            contact_nodes,
            edges,
        )?;
        closure.validate_state(&self.state)?;
        if let Some(law) = self.viscous_slit {
            let q = displacements(&self.state);
            closure.viscous_flow(law, &q[..n], self.spec.density_kg_m3, 0.0)?;
        }
        // The physical total slit length now comes from the actual boundary.
        // Pressure force and swept-volume maps remain the original face integral.
        self.spec.width_m = closure.width_m();
        self.closure = Some(closure);
        self.solid_closure = Some(spec);
        self.validate_state(&self.state)?;
        let energy = self.stored_energy_j();
        if !energy.is_finite() || energy > self.spec.limits.maximum_energy_j {
            return Err(invalid(
                "solid surface contact exceeds the complete energy allowance",
            ));
        }
        poll(gate)?;
        Ok(self)
    }

    /// Original source-bound face/edge selection and contact data.
    pub fn solid_surface_closure(&self) -> Option<&SolidSurfaceClosureSpec> {
        self.solid_closure.as_ref()
    }
}
