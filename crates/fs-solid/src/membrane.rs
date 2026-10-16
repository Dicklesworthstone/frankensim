//! Total-Lagrangian P1 membranes embedded in three dimensions.
//!
//! The original `Mesh2` owns the flat reference mesh, `fs-material::sheet`
//! owns the constitutive response, and `fs-solver` owns nonlinear iteration.
//! This module supplies only kinematics and reference-area virtual work.
//! It is not an edge-spring surrogate or a t-cubed solid-plate conversion.
//!
//! With F a 3x2 surface deformation gradient in material axes,
//! E = [(F1.F1-1)/2, (F2.F2-1)/2, F1.F2]. For work-conjugate resultants S,
//! the residual is A B^t S. The tangent action differentiates BOTH B and S;
//! the geometric-stiffness term is essential under transverse loading.
//! No bending, self-contact, wrinkling relaxation, damage, moisture evolution
//! or stability certificate is implied. A flat unstressed membrane has zero
//! transverse tangent: impose physical prestretch or an admissible initial
//! deformation, not an artificial stiffness floor.

use crate::mesh2::{Mesh2, shapes_at};
use fs_material::sheet::{MembraneLaw, ResolvedSheet, SheetError, SheetLaw, SheetResponse};
use std::collections::BTreeSet;

/// Allocation/topology limits checked before allocating element data.
#[derive(Clone, Copy, Debug)]
pub struct MembraneBudget {
    /// Maximum reference nodes (three displacement DOFs per node).
    pub max_nodes: usize,
    /// Maximum P1 triangles and source bindings.
    pub max_triangles: usize,
}

/// Caller-declared operating envelope, not an empirical source qualification.
#[derive(Clone, Copy, Debug)]
pub struct MembraneDomain {
    /// Lower bound on deformed/reference surface area, strictly positive.
    pub minimum_area_ratio: f64,
    /// Upper bound on each absolute material engineering Green-strain component.
    pub maximum_abs_strain: f64,
}

/// A refused source, mesh, state, work budget, or equilibrium solve.
#[derive(Debug)]
pub enum MembraneError {
    /// Invalid finite data, topology, dimensions, or boundary condition.
    Invalid(&'static str),
    /// Original material refusal with reference-element attribution.
    Material { element: usize, source: SheetError },
    /// The trial state left the declared geometric/strain domain.
    State { element: usize, reason: &'static str },
    /// No candidate is published after cancellation.
    Cancelled,
    /// The explicitly bounded work or topology was exhausted.
    Budget(&'static str),
    /// Original Newton configuration/initial-state refusal.
    Newton(fs_solver::nonlinear::NewtonError),
    /// The original nonlinear solver did not establish equilibrium.
    Unconverged(fs_solver::nonlinear::NewtonReport),
}
impl core::fmt::Display for MembraneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "membrane: {self:?}")
    }
}
impl std::error::Error for MembraneError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material { source, .. } => Some(source),
            Self::Newton(source) => Some(source),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
struct Element {
    nodes: [usize; 3],
    // Reference shape gradients rotated ONCE into the material axes.
    gradients: [[f64; 2]; 3],
    area: f64,
    law: MembraneLaw,
}

/// Immutable reference geometry, per-element grain and complete source bundles.
#[derive(Clone, Debug)]
pub struct MembranePatch {
    mesh: Mesh2,
    sheets: Vec<ResolvedSheet>,
    grain_rad: Vec<f64>,
    elements: Vec<Element>,
    nodal_mass_kg: Vec<f64>,
    domain: MembraneDomain,
}

/// Full-space reference-area energy and internal virtual-work residual.
#[derive(Clone, Debug)]
pub struct MembraneEvaluation {
    /// Stored strain energy, in joules.
    pub strain_energy_j: f64,
    /// Gradient of strain energy with respect to displacement, in newtons.
    /// Physical elastic nodal force is its negative.
    pub internal_gradient_n: Vec<[f64; 3]>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MembraneError> {
    if cancelled() { Err(MembraneError::Cancelled) } else { Ok(()) }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn finite_vectors(v: &[[f64; 3]], n: usize) -> bool {
    v.len() == n && v.iter().flatten().all(|x| x.is_finite())
}

impl MembranePatch {
    /// Bind complete membrane sources to the EXISTING P1 reference mesh.
    /// Areal mass remains independent of thickness, stiffness and current area.
    /// Only triangles are admitted; quadrilaterals are not silently remeshed.
    pub fn try_new(
        mesh: Mesh2, sheets: Vec<ResolvedSheet>, grain_rad: Vec<f64>,
        domain: MembraneDomain, budget: MembraneBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MembraneError> {
        poll(&mut cancelled)?;
        if mesh.nodes.len() > budget.max_nodes || mesh.elems.len() > budget.max_triangles {
            return Err(MembraneError::Budget("reference topology"));
        }
        if mesh.nodes.len() < 3 || mesh.elems.is_empty()
            || mesh.nodes.len().checked_mul(3).is_none()
            || sheets.len() != mesh.elems.len() || grain_rad.len() != mesh.elems.len()
            || mesh.nodes.iter().flatten().any(|x| !x.is_finite())
            || !domain.minimum_area_ratio.is_finite() || domain.minimum_area_ratio <= 0.0
            || domain.minimum_area_ratio > 1.0
            || !domain.maximum_abs_strain.is_finite() || domain.maximum_abs_strain <= 0.0
        {
            return Err(MembraneError::Invalid("complete finite P1 mesh, sources and positive operating envelope required"));
        }
        let mut elements = Vec::with_capacity(mesh.elems.len());
        let mut nodal_mass_kg = vec![0.0; mesh.nodes.len()];
        let mut seen = BTreeSet::new();
        for (index, conn) in mesh.elems.iter().enumerate() {
            poll(&mut cancelled)?;
            if conn.len() != 3 || conn.iter().any(|&i| i >= mesh.nodes.len()) {
                return Err(MembraneError::Invalid("in-range P1 triangles required"));
            }
            let nodes = [conn[0], conn[1], conn[2]];
            let mut key = nodes;
            key.sort_unstable();
            if key[0] == key[1] || key[1] == key[2] || !seen.insert(key) {
                return Err(MembraneError::Invalid("duplicate node or duplicate triangle"));
            }
            // Validate before entering the existing shape owner, which assumes
            // a nondegenerate element. Its Jacobian is twice the reference area.
            let [a, b, c] = nodes.map(|i| mesh.nodes[i]);
            let determinant = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
            if !determinant.is_finite() || determinant <= 0.0 || !grain_rad[index].is_finite() {
                return Err(MembraneError::Invalid("finite positive reference orientation and grain required"));
            }
            // Shift to the local origin to avoid translation-dependent loss of
            // precision in the shared isoparametric reference transformation.
            let local = [[0.0, 0.0], [b[0]-a[0], b[1]-a[1]], [c[0]-a[0], c[1]-a[1]]];
            let (_, gradients, jacobian) = shapes_at(&local, &[0, 1, 2], 0.0, 0.0);
            let (sin, cos) = grain_rad[index].sin_cos();
            let gradients: [[f64; 2]; 3] = std::array::from_fn(|i| [
                cos*gradients[i][0] + sin*gradients[i][1],
                -sin*gradients[i][0] + cos*gradients[i][1],
            ]);
            let area = 0.5*jacobian;
            let mass = sheets[index].areal_density_kg_m2()*area;
            if !area.is_finite() || area <= 0.0 || !mass.is_finite() || mass/3.0 <= 0.0
                || gradients.iter().flatten().any(|x| !x.is_finite()) {
                return Err(MembraneError::Invalid("unrepresentable reference shape gradients or mass"));
            }
            let law = sheets[index].membrane_law().map_err(|source| MembraneError::Material { element: index, source })?;
            for &i in &nodes { nodal_mass_kg[i] += mass/3.0; }
            elements.push(Element { nodes, gradients, area, law });
        }
        if nodal_mass_kg.iter().any(|m| !m.is_finite() || *m <= 0.0)
            || !nodal_mass_kg.iter().sum::<f64>().is_finite() {
            return Err(MembraneError::Invalid("unused node or unrepresentable reference mass"));
        }
        poll(&mut cancelled)?;
        Ok(Self { mesh, sheets, grain_rad, elements, nodal_mass_kg, domain })
    }
    /// Original reference mesh, including its named boundary patches.
    pub fn mesh(&self) -> &Mesh2 { &self.mesh }
    /// Complete original source objects in original triangle order.
    pub fn sheets(&self) -> &[ResolvedSheet] { &self.sheets }
    /// Material axis angles in the original reference plane, in radians.
    pub fn grain_rad(&self) -> &[f64] { &self.grain_rad }
    /// Lumped reference masses; current deformation cannot create mass.
    pub fn nodal_mass_kg(&self) -> &[f64] { &self.nodal_mass_kg }
    /// Conserved total reference mass, in kilograms.
    pub fn reference_mass_kg(&self) -> f64 { self.nodal_mass_kg.iter().sum() }

    // Columns of the surface gradient. Displacement differences (rather than
    // global current positions) make rigid translation drop out algebraically.
    fn deformation(&self, e: &Element, u: &[[f64; 3]]) -> [[f64; 3]; 2] {
        let mut f = [[0.0; 3]; 2];
        let origin = e.nodes[0];
        for a in 1..3 {
            let node = e.nodes[a];
            for i in 0..3 {
                let reference = if i < 2 { self.mesh.nodes[node][i]-self.mesh.nodes[origin][i] } else { 0.0 };
                let edge = reference + (u[node][i]-u[origin][i]);
                for (alpha, col) in f.iter_mut().enumerate() { col[i] += edge*e.gradients[a][alpha]; }
            }
        }
        f
    }
    fn response(&self, index: usize, f: [[f64; 3]; 2]) -> Result<SheetResponse, MembraneError> {
        let [a, b] = f;
        let cross = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
        let area_ratio = cross[0].hypot(cross[1]).hypot(cross[2]);
        let strain = [0.5*(dot(a,a)-1.0), 0.5*(dot(b,b)-1.0), dot(a,b)];
        if !area_ratio.is_finite() || area_ratio < self.domain.minimum_area_ratio
            || strain.iter().any(|e| !e.is_finite() || e.abs() > self.domain.maximum_abs_strain) {
            return Err(MembraneError::State { element: index, reason: "surface collapse or material strain outside declared operating envelope" });
        }
        self.elements[index].law.response(strain)
            .map_err(|source| MembraneError::Material { element: index, source })
    }
    // A semismooth equilibrium tangent is not automatically a unique vibration
    // tangent. Preserve material attribution and the same surface kinematics.
    fn admit_linearization(&self, displacement_m:&[[f64;3]],
        mut cancelled:impl FnMut()->bool) -> Result<(),MembraneError> {
        if !finite_vectors(displacement_m,self.mesh.nodes.len()) {
            return Err(MembraneError::Invalid("finite nodal displacement required"));
        }
        for (index,e) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let [a,b]=self.deformation(e,displacement_m);
            e.law.admit_linearization([0.5*(dot(a,a)-1.0),0.5*(dot(b,b)-1.0),dot(a,b)])
                .map_err(|source|MembraneError::Material {element:index,source})?;
        }
        poll(&mut cancelled)
    }
    /// Evaluate strain energy and its full displacement gradient. Outputs are
    /// staged, so cancellation/refusal cannot publish a partial force array.
    pub fn evaluate(&self, displacement_m: &[[f64; 3]], mut cancelled: impl FnMut() -> bool)
        -> Result<MembraneEvaluation, MembraneError> {
        poll(&mut cancelled)?;
        if !finite_vectors(displacement_m, self.mesh.nodes.len()) {
            return Err(MembraneError::Invalid("one finite displacement vector per reference node required"));
        }
        let mut out = MembraneEvaluation { strain_energy_j: 0.0, internal_gradient_n: vec![[0.0;3]; self.mesh.nodes.len()] };
        for (index, e) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let f = self.deformation(e, displacement_m);
            let r = self.response(index, f)?;
            out.strain_energy_j += e.area*r.energy_j_m2;
            let [s11,s22,s12] = r.resultants_n_m;
            for (a, &node) in e.nodes.iter().enumerate() {
                let [gx,gy] = e.gradients[a];
                for i in 0..3 {
                    out.internal_gradient_n[node][i] += e.area*(
                        (s11*gx+s12*gy)*f[0][i] + (s12*gx+s22*gy)*f[1][i]);
                }
            }
        }
        if !out.strain_energy_j.is_finite() || !finite_vectors(&out.internal_gradient_n, self.mesh.nodes.len()) {
            return Err(MembraneError::Invalid("unrepresentable assembled energy or force"));
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
    /// Exact matrix-free tangent action, including initial-stress/geometric
    /// stiffness. No finite differencing and no dense global Hessian allocation.
    pub fn tangent_apply(&self, displacement_m: &[[f64; 3]], direction_m: &[[f64; 3]],
        mut cancelled: impl FnMut() -> bool) -> Result<Vec<[f64; 3]>, MembraneError> {
        poll(&mut cancelled)?;
        let n = self.mesh.nodes.len();
        if !finite_vectors(displacement_m,n) || !finite_vectors(direction_m,n) {
            return Err(MembraneError::Invalid("finite displacement and direction vectors required"));
        }
        let mut out = vec![[0.0;3];n];
        for (index,e) in self.elements.iter().enumerate() {
            poll(&mut cancelled)?;
            let f = self.deformation(e,displacement_m);
            let r = self.response(index,f)?;
            let mut df = [[0.0;3];2];
            for a in 1..3 { for i in 0..3 { for alpha in 0..2 {
                df[alpha][i] += (direction_m[e.nodes[a]][i]-direction_m[e.nodes[0]][i])*e.gradients[a][alpha];
            } } }
            let de = [dot(f[0],df[0]),dot(f[1],df[1]),dot(f[0],df[1])+dot(f[1],df[0])];
            let ds: [f64;3] = std::array::from_fn(|i| (0..3).map(|j|r.tangent_n_m[3*i+j]*de[j]).sum());
            let [s11,s22,s12] = r.resultants_n_m;
            let [ds11,ds22,ds12] = ds;
            for (a,&node) in e.nodes.iter().enumerate() {
                let [gx,gy] = e.gradients[a];
                for i in 0..3 {
                    out[node][i] += e.area*((ds11*gx+ds12*gy)*f[0][i] + (ds12*gx+ds22*gy)*f[1][i]
                        + (s11*gx+s12*gy)*df[0][i] + (s12*gx+s22*gy)*df[1][i]);
                }
            }
        }
        if !finite_vectors(&out,n) { return Err(MembraneError::Invalid("unrepresentable tangent action")); }
        poll(&mut cancelled)?;
        Ok(out)
    }
    /// Uniform body acceleration to consistent P1 nodal dead forces, in newtons.
    /// No mass is inferred from bulk-fiber density or current sheet volume.
    pub fn body_load_n(&self, acceleration_m_s2: [f64;3]) -> Result<Vec<[f64;3]>,MembraneError> {
        if acceleration_m_s2.iter().any(|x|!x.is_finite()) { return Err(MembraneError::Invalid("finite acceleration required")); }
        let out: Vec<_> = self.nodal_mass_kg.iter().map(|m|acceleration_m_s2.map(|a|m*a)).collect();
        if !finite_vectors(&out,self.mesh.nodes.len()) { return Err(MembraneError::Invalid("unrepresentable body load")); }
        Ok(out)
    }
}

mod equilibrium;
pub use equilibrium::{MembraneSolveSettings, MembraneSolution};

/// Prestressed vibration and work-conjugate force/volume ports.
pub mod vibration;
