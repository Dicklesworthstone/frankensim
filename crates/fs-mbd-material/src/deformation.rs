//! Transport a reference material inventory through a compatible P1 tet map.
//!
//! This is geometry/mass transport, not another constitutive or dynamics solver.
//! One reference region is one tetrahedron; shared vertex indices enforce the
//! same displacement on adjoining cells. Every cell retains its original mass
//! while its volume, density, centroid and complete inertia follow actual nodes.
//! Overlap detection, stress validity and momentum/work during a live change
//! remain with the geometry, constitutive and time-integration owners.

/// Material-driven constrained thermal deformation through the existing FEM owner.
#[cfg(feature = "thermoelastic-coupling")]
pub mod thermoelastic;

use std::collections::BTreeSet;
use fs_ga::Vec3;
use fs_mbd::articulated::SpatialInertia;
use crate::{MaterialInertiaError, ResolvedBody, VolumeMoments, merge_inertias, poll, positive};

/// A fixed-topology P1 mesh bound to the exact reference mass inventory.
/// Each region must have been resolved using `VolumeMoments::tetrahedron` on
/// the supplied vertices, in the same order. No tolerance silently changes the
/// original geometry. One region per cell also permits nonuniform temperatures.
#[derive(Clone, Debug)]
pub struct TetrahedralMass<'a> {
    reference: &'a ResolvedBody,
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    // Canonical region -> element, independent of assignment input order.
    elements: Vec<usize>,
    determinants: Vec<f64>,
}

/// Actual deformed geometry and the SAME material mass for a reference region.
#[derive(Clone, Debug, PartialEq)]
pub struct DeformedRegion {
    geometry: VolumeMoments,
    inertia: SpatialInertia,
}
impl DeformedRegion {
    /// Current geometric moments, in the original body's coordinate frame.
    pub const fn geometry(&self) -> VolumeMoments { self.geometry }
    /// Canonical spatial inertia with reference mass copied, not re-evaluated.
    pub const fn inertia(&self) -> SpatialInertia { self.inertia }
    /// Current bulk density from reference mass divided by current volume.
    pub fn density_kg_m3(&self) -> f64 { self.inertia.mass() / self.geometry.volume_m3() }
}

/// A complete immutable geometry/mass configuration, not a live state mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct DeformedBody<'a> {
    reference: &'a ResolvedBody,
    nodes_m: Vec<[f64; 3]>,
    tetrahedra: Vec<[usize; 4]>,
    regions: Vec<DeformedRegion>,
    inertia: SpatialInertia,
}
impl<'a> DeformedBody<'a> {
    /// Original material cards, density receipts and geometry are retained here.
    pub const fn reference(&self) -> &'a ResolvedBody { self.reference }
    /// Actual updated shared mesh vertices [m].
    pub fn nodes_m(&self) -> &[[f64; 3]] { &self.nodes_m }
    /// Original fixed-topology element order.
    pub fn tetrahedra(&self) -> &[[usize; 4]] { &self.tetrahedra }
    /// Current regions in the reference body's canonical region-key order.
    pub fn regions(&self) -> &[DeformedRegion] { &self.regions }
    /// Feed directly to the existing body or articulated dynamics owner.
    pub const fn inertia(&self) -> SpatialInertia { self.inertia }
}

impl<'a> TetrahedralMass<'a> {
    /// Bind a complete fixed-topology mesh before any deformation.
    ///
    /// `region_keys[e]` identifies the reference region occupied by element e.
    /// All regions are required once, and duplicate cells refuse. `max_nodes`
    /// bounds the vertex payload before allocation; cells inherit the reference
    /// body's existing region cap. Coordinates are metres in ONE body frame.
    /// Exact reference moments prevent attaching a displaced/foreign mesh to a
    /// valid mass receipt. Global geometric overlap is not certified here.
    pub fn bind(
        reference: &'a ResolvedBody,
        nodes_m: &[[f64; 3]],
        tetrahedra: &[[usize; 4]],
        region_keys: &[&str],
        max_nodes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialInertiaError> {
        poll(&mut cancelled)?;
        if nodes_m.is_empty() || nodes_m.len() > max_nodes
            || tetrahedra.len() != reference.regions().len()
            || region_keys.len() != tetrahedra.len() {
            return Err(invalid("complete bounded reference mesh and one region key per cell required"));
        }
        for node in nodes_m {
            poll(&mut cancelled)?;
            if !node.iter().all(|v| v.is_finite()) { return Err(invalid("reference nodes must be finite")); }
        }
        let mut elements = vec![usize::MAX; reference.regions().len()];
        let mut determinants = Vec::with_capacity(tetrahedra.len());
        let mut cells = BTreeSet::new();
        for (element, (&tet, &key)) in tetrahedra.iter().zip(region_keys).enumerate() {
            poll(&mut cancelled)?;
            if tet.iter().any(|&v| v >= nodes_m.len()) { return Err(invalid("tetrahedron node is out of range")); }
            let mut canonical = tet;
            canonical.sort_unstable();
            if canonical.windows(2).any(|w| w[0] == w[1]) || !cells.insert(canonical) {
                return Err(invalid("tetrahedra must have four distinct nodes and cannot be duplicated"));
            }
            let index = reference.regions().binary_search_by(|r| r.key().cmp(key))
                .map_err(|_| invalid("unknown reference region"))?;
            if elements[index] != usize::MAX { return Err(invalid("reference region is assigned more than once")); }
            let vertices = tet.map(|v| vector(nodes_m[v]));
            let determinant = signed_determinant(vertices);
            if !determinant.is_finite() || determinant == 0.0 {
                return Err(invalid("reference tetrahedron has unrepresentable signed volume"));
            }
            let moments = VolumeMoments::tetrahedron(vertices)?;
            if moments != reference.regions()[index].geometry() {
                return Err(invalid("tetrahedron geometry differs from its exact reference mass inventory"));
            }
            elements[index] = element;
            determinants.push(determinant);
        }
        poll(&mut cancelled)?;
        Ok(Self { reference, nodes_m: nodes_m.to_vec(), tetrahedra: tetrahedra.to_vec(), elements, determinants })
    }

    /// Unchanged reference inventory, including original density/source records.
    pub const fn reference(&self) -> &'a ResolvedBody { self.reference }
    /// Reference coordinates used by the physical displacement solver.
    pub fn nodes_m(&self) -> &[[f64; 3]] { &self.nodes_m }
    /// Fixed topology and element order expected by the physical solver.
    pub fn tetrahedra(&self) -> &[[usize; 4]] { &self.tetrahedra }
    /// Canonical reference region -> solver element map.
    pub fn region_elements(&self) -> &[usize] { &self.elements }

    /// Apply actual current node coordinates, conserving each reference mass.
    ///
    /// A P1 tetrahedron has an affine map and constant Jacobian, so its uniform
    /// reference density remains uniform in that cell. Analytic current volume
    /// moments come from the existing geometry owner. No isotropic scale is
    /// substituted for shear, bending or incompatible cellwise free strains.
    /// Inverted/collapsed cells, nonfinite arithmetic and cancellation refuse
    /// the entire candidate. Repeated calls always start from reference matter.
    pub fn deform(
        &self, nodes_m: &[[f64; 3]], mut cancelled: impl FnMut() -> bool,
    ) -> Result<DeformedBody<'a>, MaterialInertiaError> {
        poll(&mut cancelled)?;
        if nodes_m.len() != self.nodes_m.len() { return Err(invalid("deformation must retain every reference vertex")); }
        for node in nodes_m {
            poll(&mut cancelled)?;
            if !node.iter().all(|v| v.is_finite()) { return Err(invalid("deformed nodes must be finite")); }
        }
        let mut regions = Vec::with_capacity(self.elements.len());
        let mut combined = None;
        for (index, &element) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let tet = self.tetrahedra[element];
            let vertices = tet.map(|v| vector(nodes_m[v]));
            let determinant = signed_determinant(vertices);
            if !determinant.is_finite() || determinant == 0.0
                || determinant.is_sign_positive() != self.determinants[element].is_sign_positive() {
                return Err(invalid("deformation inverts, collapses or overflows a material cell"));
            }
            let original = &self.reference.regions()[index];
            let geometry = VolumeMoments::tetrahedron(vertices)?;
            let density = original.inertia().mass() / geometry.volume_m3();
            if !positive(density) { return Err(invalid("transported density is unrepresentable")); }
            let inertia = if tet.iter().all(|&v| nodes_m[v] == self.nodes_m[v]) {
                original.inertia()
            } else {
                let mut tensor = geometry.inertia_m5();
                for entry in &mut tensor.m { *entry *= density; }
                // Do not multiply density by volume again: preserve mass bits.
                SpatialInertia::new(original.inertia().mass(), geometry.centroid_m(), tensor)?
            };
            combined = Some(match combined { None => inertia, Some(previous) => merge_inertias(previous, inertia)? });
            regions.push(DeformedRegion { geometry, inertia });
        }
        let inertia = combined.expect("the reference inventory is nonempty");
        inertia.matrix()?;
        poll(&mut cancelled)?;
        Ok(DeformedBody { reference: self.reference, nodes_m: nodes_m.to_vec(),
            tetrahedra: self.tetrahedra.clone(), regions, inertia })
    }
}

fn invalid(what: &'static str) -> MaterialInertiaError { MaterialInertiaError::Invalid(what) }
fn vector([x,y,z]: [f64;3]) -> Vec3 { Vec3::new(x,y,z) }
fn signed_determinant(p: [Vec3;4]) -> f64 { (p[1]-p[0]).dot((p[2]-p[0]).cross(p[3]-p[0])) }
