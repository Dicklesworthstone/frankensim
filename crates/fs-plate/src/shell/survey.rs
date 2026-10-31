//! Explicit, non-axisymmetric midsurface samples and physical thickness.
//!
//! This is a cold input adapter into the existing shell and section owners, not
//! a second FEM or a reconstruction algorithm. It retains the supplied vertices,
//! triangles, nodal thickness and per-facet isotropic material. No smoothing,
//! welding, hole filling, resampling or inferred manufacturing data is applied.
use super::{ShellMesh, ShellModel, ShellSupport, FacetGeometry, CstDktMass, assemble_shell_sections_with_mass};
use super::profile::{ProfileBudget, ProfileShell};
use crate::{PlateError, PlateSection};
use std::collections::{BTreeMap, BTreeSet};

/// Isotropic physical material for one facet; thickness is a separate nodal field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IsotropicMaterial {
    /// Young's modulus [Pa].
    pub young_pa: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Volume density [kg/m^3].
    pub density_kg_m3: f64,
}

/// Explicit facet material. Orthotropic axis 1 is a supplied GLOBAL direction:
/// it is projected into each actual facet tangent before stiffness rotation.
/// One material ID can therefore cross mesh diagonals without changing its
/// physical grain/rolling direction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SurveyMaterial {
    Isotropic(IsotropicMaterial),
    Orthotropic {
        e1_pa: f64,
        e2_pa: f64,
        nu12: f64,
        g12_pa: f64,
        density_kg_m3: f64,
        axis1_global: [f64; 3],
    },
}
impl From<IsotropicMaterial> for SurveyMaterial {
    fn from(value:IsotropicMaterial)->Self { Self::Isotropic(value) }
}
impl SurveyMaterial {
    fn section(self,h:f64,facet:&FacetGeometry)->Result<PlateSection,PlateError> {
        match self {
            Self::Isotropic(m)=>PlateSection::isotropic(m.young_pa,m.poisson,h,m.density_kg_m3),
            Self::Orthotropic{e1_pa,e2_pa,nu12,g12_pa,density_kg_m3,axis1_global}=>{
                let angle=material_axis_angle(facet,axis1_global)?;
                PlateSection::orthotropic_plane_stress_at_angle(
                    e1_pa,e2_pa,nu12,g12_pa,h,density_kg_m3,angle)
            }
        }
    }
}

/// Project a physical GLOBAL material direction onto an actual facet tangent.
/// This preserves grain across differently ordered triangle edges. A normal or
/// unresolved direction refuses; no fallback direction or smoothing is inferred.
pub fn material_axis_angle(facet: &FacetGeometry, axis1_global: [f64;3]) -> Result<f64,PlateError> {
    if axis1_global.iter().any(|v|!v.is_finite()) {
        return Err(bad("orthotropic shell material axis must be finite"));
    }
    let axis_norm=axis1_global.iter().map(|v|v*v).sum::<f64>().sqrt();
    if !axis_norm.is_finite() || axis_norm<=0.0 {
        return Err(bad("orthotropic shell material needs a nonzero global axis"));
    }
    let dot=|a:[f64;3],b:[f64;3]|a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    let x=dot(axis1_global,facet.frame[0]);
    let y=dot(axis1_global,facet.frame[1]);
    let tangent=x.hypot(y);
    if !tangent.is_finite() || tangent<=1e-10*axis_norm {
        return Err(bad("orthotropic shell material axis is unresolved on a facet tangent"));
    }
    Ok(y.atan2(x))
}

/// One connected, consistently oriented, manifold triangular shell.
/// Admitted topology is not a global geometric self-intersection certificate.
#[derive(Debug, Clone)]
pub struct MeshShell {
    /// Original reference midsurface positions and connectivity, in input order.
    pub mesh: ShellMesh,
    /// Existing section law, with nodal thickness averaged within each facet.
    pub sections: Vec<PlateSection>,
    /// Original physical thickness [m], also consumed by exterior-surface lifting.
    pub nodal_thickness_m: Vec<f64>,
    /// Sum rho*h*A over actual three-dimensional facets [kg].
    pub mass_kg: f64,
    /// Longest actual three-dimensional facet edge [m].
    pub max_edge_m: f64,
}
fn bad(what: &'static str) -> PlateError { PlateError::BadSection { what } }

impl MeshShell {
    /// Admit explicit samples without changing their geometry or indexing.
    /// Each node must be used, every vertex link must be one fan, and triangles
    /// must connect through consistently oriented edges. Boundaries and multiple
    /// holes are allowed; disconnected shells must be separate physical bodies.
    ///
    /// # Errors
    /// Refuses malformed geometry, thickness/material counts, nonmanifold or
    /// disconnected topology, exact duplicate positions/facets and mesh budgets.
    /// Near duplicates and global self-intersections are not certified here.
    pub fn new(mesh: ShellMesh, nodal_thickness_m: Vec<f64>,
        materials: &[IsotropicMaterial], budget: ProfileBudget) -> Result<Self, PlateError>
    {
        let materials=materials.iter().copied().map(SurveyMaterial::from).collect::<Vec<_>>();
        Self::new_with_materials(mesh,nodal_thickness_m,&materials,budget)
    }

    /// Same explicit geometry with isotropic or physically oriented orthotropic
    /// material per facet. The global material direction is never inferred from
    /// triangle edge order.
    pub fn new_with_materials(mesh: ShellMesh, nodal_thickness_m: Vec<f64>,
        materials: &[SurveyMaterial], budget: ProfileBudget) -> Result<Self, PlateError>
    {
        let n = mesh.nodes.len(); let count = mesh.tris.len();
        if n < 3 || count == 0 || n > budget.max_nodes || count > budget.max_triangles
            || nodal_thickness_m.len() != n || materials.len() != count
        { return Err(bad("shell survey exceeds mesh budgets or has incomplete physical fields")); }
        if nodal_thickness_m.iter().any(|h| !h.is_finite() || *h <= 0.0)
            || mesh.nodes.iter().flatten().any(|x| !x.is_finite())
        { return Err(bad("shell survey requires finite positions and positive finite thickness")); }
        validate_midsurface(&mesh, budget)?;
        let mut sections = Vec::with_capacity(count);
        let (mut mass, mut max_edge) = (0.0, 0.0_f64);
        for (e, tri) in mesh.tris.iter().enumerate() {
            let facet = mesh.facet(e)?;
            let h = (nodal_thickness_m[tri[0]]+nodal_thickness_m[tri[1]]+nodal_thickness_m[tri[2]])/3.0;
            let section = materials[e].section(h,&facet)?;
            mass += section.density*section.thickness*facet.area_m2;
            sections.push(section);
            for k in 0..3 {
                let (a,b) = (tri[k],tri[(k+1)%3]);
                let p = mesh.nodes[a]; let q = mesh.nodes[b];
                max_edge = max_edge.max((p[0]-q[0]).hypot(p[1]-q[1]).hypot(p[2]-q[2]));
            }
        }
        if !mass.is_finite() || mass <= 0.0 || !max_edge.is_finite() {
            return Err(bad("shell survey mass or geometric scale is not representable"));
        }
        Ok(Self { mesh, sections, nodal_thickness_m, mass_kg: mass, max_edge_m: max_edge })
    }

    /// Assemble the existing CST/DKT stiffness with consistent P1/director
    /// inertia, including the full orientation-dependent through-thickness tensor.
    /// # Errors
    /// Propagates the existing physical-section, support and shell assembly errors.
    pub fn assemble(&self, nodes: &[usize], support: ShellSupport) -> Result<ShellModel, PlateError> {
        self.assemble_with_mass(nodes, support, CstDktMass::Consistent)
    }
    /// Explicit inertia choice for comparison with an original lumped pencil.
    /// # Errors
    /// The same physical-section and support refusals as [`Self::assemble`].
    pub fn assemble_with_mass(&self, nodes: &[usize], support: ShellSupport, mass: CstDktMass)
        -> Result<ShellModel, PlateError> {
        assemble_shell_sections_with_mass(&self.mesh, &self.sections, nodes, support, mass)
    }
}

/// Admit the original survey's connected, oriented manifold midsurface without
/// requiring fictitious material or thickness fields. Explicit mesh consumers
/// share these SAME edge/fan checks before constructing their physical sections.
/// Global self-intersections and near-coincident vertices remain outside scope.
/// ProfileBudget::max_feature_evaluations is irrelevant to an explicit mesh.
pub fn validate_midsurface(mesh: &ShellMesh, budget: ProfileBudget) -> Result<(), PlateError> {
    let n = mesh.nodes.len(); let count = mesh.tris.len();
    if n < 3 || count == 0 || n > budget.max_nodes || count > budget.max_triangles {
        return Err(bad("shell survey exceeds explicit mesh budgets"));
    }
    if mesh.nodes.iter().flatten().any(|x| !x.is_finite()) {
        return Err(bad("shell survey requires finite positions"));
    }
    let mut positions = BTreeSet::new();
    for node in &mesh.nodes {
        // +0/-0 are the same position, not two physical vertices.
        let key = node.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
        if !positions.insert(key) { return Err(bad("shell survey contains duplicate positions; explicit welding is required")); }
    }
    let mut edges = BTreeMap::<(usize, usize), (usize, usize, usize)>::new();
    let mut facets = BTreeSet::new();
    let mut adjacency = vec![Vec::new(); count];
    let mut links = vec![Vec::<(usize, usize)>::new(); n];
    for (e, tri) in mesh.tris.iter().enumerate() {
        mesh.facet(e)?;
        let mut key = *tri; key.sort_unstable();
        if !facets.insert(key) { return Err(bad("shell survey contains duplicate facets")); }
        for k in 0..3 {
            let (a, b, c) = (tri[k], tri[(k+1)%3], tri[(k+2)%3]);
            links[a].push((b, c));
            let key = (a.min(b), a.max(b));
            if let Some(&(first, from, to)) = edges.get(&key) {
                if from != b || to != a || first == usize::MAX {
                    return Err(bad("shell survey has inconsistent winding or a nonmanifold edge"));
                }
                adjacency[e].push(first); adjacency[first].push(e);
                edges.insert(key, (usize::MAX, a, b));
            } else { edges.insert(key, (e, a, b)); }
        }
    }
    let mut seen = vec![false; count]; let mut stack = vec![0]; seen[0] = true;
    while let Some(e) = stack.pop() {
        for &next in &adjacency[e] { if !seen[next] { seen[next] = true; stack.push(next); } }
    }
    if seen.iter().any(|v| !v) { return Err(bad("shell survey must contain one edge-connected physical body")); }
    for link in links {
        if link.is_empty() { return Err(bad("shell survey contains an unused vertex")); }
        let mut graph = BTreeMap::<usize, Vec<usize>>::new();
        for (a, b) in link { graph.entry(a).or_default().push(b); graph.entry(b).or_default().push(a); }
        let ends = graph.values().filter(|v| v.len() == 1).count();
        if graph.values().any(|v| v.len() > 2) || (ends != 0 && ends != 2) {
            return Err(bad("shell survey vertex is not a manifold fan"));
        }
        let mut reached = BTreeSet::new(); let mut pending = vec![*graph.keys().next().unwrap()];
        while let Some(a) = pending.pop() {
            if reached.insert(a) { pending.extend(graph[&a].iter().copied()); }
        }
        if reached.len() != graph.len() { return Err(bad("shell survey vertex joins disconnected fans")); }
    }
    Ok(())
}

impl From<ProfileShell> for MeshShell {
    /// Retain a previously generated profile without recomputation or rounding.
    fn from(s: ProfileShell) -> Self {
        Self { mesh: s.mesh, sections: s.sections, nodal_thickness_m: s.nodal_thickness_m,
            mass_kg: s.mass_kg, max_edge_m: s.max_edge_m }
    }
}
