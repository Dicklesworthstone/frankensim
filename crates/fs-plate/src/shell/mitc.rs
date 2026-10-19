//! Linear shear-deformable triangles with two condensed director rotations.
//!
//! Rotation interpolation and the six interior shear tying points follow
//! Lee, Lee & Bathe (2014), doi:10.1016/j.compstruc.2014.02.005, Eqs. 6–8,
//! 15–17 and Table 1. Their nonzero twisting distance is 1/10000. Geometry
//! remains a convergent faceted midsurface with independently supplied nodal
//! directors. CST/DKT operators are not called or modified by this path.
//!
//! Local support coordinates travel with the pencil: `ShellModel::dof_map`
//! alone cannot reconstruct an oblique support. Use `MitcShellModel`'s load
//! projection and displacement expansion for physical work and observation.

use super::{ShellMesh, ShellModel, cross, dot};
use crate::PlateError;
use fs_material::sheet::SheetMatrix;
use fs_math::det;
use fs_sparse::{Coo, DirectOrdering, LdltOptions, SymbolicLdlt};
use std::collections::BTreeMap;

mod element;
/// Full membrane, bending and transverse-shear material dissipation.
pub mod viscosity;

/// Explicit displacement/director inertia discretization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellMass {
    /// Exact polynomial mass, congruently condensed with the bubble rotations.
    Consistent,
    /// Nodal translational and rotary inertia, without a bubble inertia field.
    Lumped,
}

/// Caller-owned cold-work ceilings, checked before assembly allocation.
#[derive(Debug, Clone, Copy)]
pub struct MitcBudget {
    /// Maximum physical midsurface nodes.
    pub max_nodes: usize,
    /// Maximum triangular elements.
    pub max_triangles: usize,
}

/// Rows are local support axes expressed in global coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConstraintFrame {
    /// Cartesian world axes.
    Global,
    /// Explicit proper orthonormal frame, shared by translations and rotations.
    Local([[f64; 3]; 3]),
}

/// One homogeneous support declaration per node; combine masks before calling.
#[derive(Debug, Clone, Copy)]
pub struct ShellConstraint {
    /// Physical node.
    pub node: usize,
    /// Fixed `(u,v,w,theta_x,theta_y,theta_z)` components in `frame`.
    pub mask: [bool; 6],
    /// Work-conjugate support axes.
    pub frame: ConstraintFrame,
}

/// Unit directors and the explicit component-surface angle envelope.
#[derive(Debug, Clone)]
pub struct NodalDirectors {
    vectors: Vec<[f64; 3]>,
    minimum_facet_dot: f64,
}

impl NodalDirectors {
    /// Admit unit vectors, their alignment with every incident facet, and
    /// pairwise incident-facet normals. Junctions must use separate components.
    pub fn new(
        mesh: &ShellMesh,
        vectors: Vec<[f64; 3]>,
        minimum_facet_dot: f64,
    ) -> Result<Self, PlateError> {
        let out = Self {
            vectors,
            minimum_facet_dot,
        };
        out.validate(mesh)?;
        Ok(out)
    }

    /// Source-supplied or reconstructed unit directors, without hidden smoothing.
    #[must_use]
    pub fn vectors(&self) -> &[[f64; 3]] {
        &self.vectors
    }

    fn validate(&self, mesh: &ShellMesh) -> Result<(), PlateError> {
        if self.vectors.len() != mesh.nodes.len()
            || !self.minimum_facet_dot.is_finite()
            || !(0.0..=1.0).contains(&self.minimum_facet_dot)
            || self.minimum_facet_dot == 0.0
        {
            return Err(PlateError::BadDirector {
                node: self.vectors.len().min(mesh.nodes.len()),
            });
        }
        mesh.validate()?;
        let mut incident = vec![Vec::new(); mesh.nodes.len()];
        for (node, n) in self.vectors.iter().enumerate() {
            if n.iter().any(|v| !v.is_finite()) || (dot(*n, *n) - 1.0).abs() > 1e-12 {
                return Err(PlateError::BadDirector { node });
            }
        }
        for (e, t) in mesh.tris.iter().enumerate() {
            let normal = mesh.facet(e)?.frame[2];
            for &node in t {
                if dot(normal, self.vectors[node]) < self.minimum_facet_dot
                    || incident[node]
                        .iter()
                        .any(|&old| dot(normal, old) < self.minimum_facet_dot)
                {
                    return Err(PlateError::BadDirector { node });
                }
                incident[node].push(normal);
            }
        }
        if let Some(node) = incident.iter().position(Vec::is_empty) {
            return Err(PlateError::BadDirector { node });
        }
        Ok(())
    }
}

/// Independent membrane, bending and transverse-shear resultants. There is no
/// membrane/bending coupling B and no loss or prestress in this linear section.
#[derive(Debug, Clone, Copy)]
pub struct ShearShellSection {
    membrane: SheetMatrix,
    bending: SheetMatrix,
    shear_n_m: [f64; 3],
    shear_correction: f64,
    areal_mass_kg_m2: f64,
    rotary_inertia_kg: f64,
    material_angle_rad: f64,
}

impl ShearShellSection {
    /// Independent principal-frame A and D, symmetric shear `[K11,K12,K22]`,
    /// explicit correction, areal mass, rotary inertia and in-plane angle.
    /// `shear_n_m` is BEFORE multiplication by the correction factor.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        membrane: SheetMatrix,
        bending: SheetMatrix,
        shear_n_m: [f64; 3],
        shear_correction: f64,
        areal_mass_kg_m2: f64,
        rotary_inertia_kg: f64,
        material_angle_rad: f64,
    ) -> Result<Self, PlateError> {
        let out = Self {
            membrane,
            bending,
            shear_n_m,
            shear_correction,
            areal_mass_kg_m2,
            rotary_inertia_kg,
            material_angle_rad,
        };
        out.validate()?;
        Ok(out)
    }

    /// Homogeneous orthotropic solid, SI units. Both transverse shear moduli
    /// and the correction factor are required; isotropic values are not inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn from_orthotropic_solid(
        e1: f64,
        e2: f64,
        nu12: f64,
        g12: f64,
        g13: f64,
        g23: f64,
        h: f64,
        rho: f64,
        angle: f64,
        kappa: f64,
    ) -> Result<Self, PlateError> {
        if [e1, e2, g12, g13, g23, h, rho]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
            || !nu12.is_finite()
        {
            return Err(PlateError::BadSection {
                what: "positive finite moduli, thickness and density required",
            });
        }
        let denominator = 1.0 - nu12 * nu12 * e2 / e1;
        if !(denominator > 0.0) {
            return Err(PlateError::BadSection {
                what: "orthotropic plane-stress compliance is not positive definite",
            });
        }
        let sheet = |scale| {
            SheetMatrix::orthotropic(
                e1 * scale / denominator,
                e2 * scale / denominator,
                nu12 * e2 * scale / denominator,
                g12 * scale,
            )
            .map_err(|_| PlateError::BadSection {
                what: "unrepresentable orthotropic shell section",
            })
        };
        Self::new(
            sheet(h)?,
            sheet(h * h * h / 12.0)?,
            [g13 * h, 0.0, g23 * h],
            kappa,
            rho * h,
            rho * h * h * h / 12.0,
            angle,
        )
    }

    fn validate(&self) -> Result<(), PlateError> {
        let [a, b, c] = self.shear_n_m;
        if [
            a,
            b,
            c,
            self.shear_correction,
            self.areal_mass_kg_m2,
            self.rotary_inertia_kg,
            self.material_angle_rad,
        ]
        .iter()
        .any(|v| !v.is_finite())
            || a <= 0.0
            || c <= 0.0
            || (b / det::sqrt(a) / det::sqrt(c)).abs() >= 1.0
            || self.shear_correction <= 0.0
            || self.areal_mass_kg_m2 <= 0.0
            || self.rotary_inertia_kg <= 0.0
        {
            return Err(PlateError::BadSection {
                what: "positive-definite shear, explicit positive correction and inertia required",
            });
        }
        Ok(())
    }

    /// Preserve the physical section while expressing its material axes in
    /// another facet frame. Geometry-owned global-axis projection supplies angle.
    pub fn with_material_angle_rad(mut self, angle: f64) -> Result<Self, PlateError> {
        self.material_angle_rad=angle;
        self.validate()?;
        Ok(self)
    }

    /// Physical areal mass [kg/m²].
    #[must_use]
    pub const fn areal_mass_kg_m2(&self) -> f64 {
        self.areal_mass_kg_m2
    }

    /// Principal A and D, prior to the explicit facet material-angle transform.
    #[must_use]
    pub const fn principal_matrices(&self) -> (&SheetMatrix, &SheetMatrix) {
        (&self.membrane, &self.bending)
    }

    fn matrices(&self) -> ([f64; 9], [f64; 9], [[f64; 2]; 2]) {
        let s = det::sin(self.material_angle_rad);
        let c = det::cos(self.material_angle_rad);
        let t = [
            [c * c, s * s, c * s],
            [s * s, c * c, -c * s],
            [-2.0 * c * s, 2.0 * c * s, c * c - s * s],
        ];
        let rotate = |matrix: &[f64; 9]| {
            let mut out = [0.0; 9];
            for i in 0..3 {
                for j in i..3 {
                    for a in 0..3 {
                        for b in 0..3 {
                            out[3 * i + j] += t[a][i] * matrix[3 * a + b] * t[b][j];
                        }
                    }
                    out[3 * j + i] = out[3 * i + j];
                }
            }
            out
        };
        let [a, b, d] = self.shear_n_m;
        let r = [[c, s], [-s, c]];
        let ks = [[a, b], [b, d]];
        let mut shear = [[0.0; 2]; 2];
        for i in 0..2 {
            for j in 0..2 {
                for a in 0..2 {
                    for b in 0..2 {
                        shear[i][j] += self.shear_correction * r[a][i] * ks[a][b] * r[b][j];
                    }
                }
            }
        }
        (
            rotate(self.membrane.coefficients()),
            rotate(self.bending.coefficients()),
            shear,
        )
    }
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Ordinary sparse pencil plus the indispensable physical-coordinate map.
#[derive(Debug, Clone)]
pub struct MitcShellModel {
    pencil: ShellModel,
    frames: Vec<[[f64; 3]; 3]>,
    mass_model: ShellMass,
}
impl MitcShellModel {
    /// The existing eigenproblem carrier, in the declared support coordinates.
    #[must_use]
    pub const fn pencil(&self) -> &ShellModel {
        &self.pencil
    }
    /// Retained caller mass choice.
    #[must_use]
    pub const fn mass_model(&self) -> ShellMass {
        self.mass_model
    }
    /// Pull back a physical global force/moment vector by the same frame map
    /// used in assembly, preserving `f_global · u_global = f_free · u_free`.
    pub fn project_load(&self, full: &[f64]) -> Result<Vec<f64>, PlateError> {
        if full.len() != 6 * self.frames.len() || full.iter().any(|x| !x.is_finite()) {
            return Err(PlateError::BadConstraint {
                what: "one finite physical force/moment per shell DOF required",
            });
        }
        let mut out = vec![0.0; self.pencil.free];
        for (node, r) in self.frames.iter().enumerate() {
            for block in [0, 3] {
                for a in 0..3 {
                    if let Some(i) = self.pencil.dof_map[6 * node + block + a] {
                        out[i] = (0..3).map(|k| r[a][k] * full[6 * node + block + k]).sum();
                    }
                }
            }
        }
        Ok(out)
    }
    /// Expand a reduced vector into physical global displacements/rotations.
    pub fn expand(&self, reduced: &[f64]) -> Result<Vec<f64>, PlateError> {
        if reduced.len() != self.pencil.free || reduced.iter().any(|x| !x.is_finite()) {
            return Err(PlateError::BadConstraint {
                what: "one finite value per free shell coordinate required",
            });
        }
        let mut full = vec![0.0; 6 * self.frames.len()];
        for (node, r) in self.frames.iter().enumerate() {
            for block in [0, 3] {
                for a in 0..3 {
                    if let Some(i) = self.pencil.dof_map[6 * node + block + a] {
                        for k in 0..3 {
                            full[6 * node + block + k] += r[a][k] * reduced[i];
                        }
                    }
                }
            }
        }
        Ok(full)
    }
}

fn valid_frame(r: [[f64; 3]; 3]) -> bool {
    r.iter().flatten().all(|v| v.is_finite())
        && (0..3).all(|i| {
            (0..3).all(|j| (dot(r[i], r[j]) - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12)
        })
        && dot(cross(r[0], r[1]), r[2]) > 0.0
}

/// Assemble only the shear-deformable element, with explicit drilling strength.
/// All supports are homogeneous, and all sections are supplied per triangle.
/// Cancellation is polled before admission, per element and per quadrature point.
#[allow(clippy::too_many_arguments)]
pub fn assemble_mitc3p(
    mesh: &ShellMesh,
    directors: &NodalDirectors,
    sections: &[ShearShellSection],
    constraints: &[ShellConstraint],
    mass: ShellMass,
    budget: MitcBudget,
    drilling_alpha: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MitcShellModel, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    for (what, requested, limit) in [
        ("nodes", mesh.nodes.len(), budget.max_nodes),
        ("triangles", mesh.tris.len(), budget.max_triangles),
    ] {
        if requested > limit {
            return Err(PlateError::Budget {
                what,
                requested,
                limit,
            });
        }
    }
    if !drilling_alpha.is_finite() || drilling_alpha <= 0.0 {
        return Err(PlateError::BadSection {
            what: "positive finite explicit drilling strength required",
        });
    }
    directors.validate(mesh)?;
    if sections.len() != mesh.tris.len() {
        return Err(PlateError::SectionCount {
            expected: mesh.tris.len(),
            actual: sections.len(),
        });
    }
    for section in sections {
        section.validate()?;
    }
    let full = mesh
        .nodes
        .len()
        .checked_mul(6)
        .ok_or(PlateError::BadSection {
            what: "shell DOF extent overflow",
        })?;
    let mut frames = vec![IDENTITY; mesh.nodes.len()];
    let mut fixed = vec![false; full];
    let mut seen = vec![false; mesh.nodes.len()];
    for constraint in constraints {
        if constraint.node >= mesh.nodes.len() || seen[constraint.node] {
            return Err(PlateError::BadConstraint {
                what: "support node out of range or declared more than once",
            });
        }
        let r = match constraint.frame {
            ConstraintFrame::Global => IDENTITY,
            ConstraintFrame::Local(r) => r,
        };
        if !valid_frame(r) {
            return Err(PlateError::BadConstraint {
                what: "support frame must be proper orthonormal",
            });
        }
        seen[constraint.node] = true;
        frames[constraint.node] = r;
        fixed[6 * constraint.node..6 * constraint.node + 6].copy_from_slice(&constraint.mask);
    }
    let mut free = 0;
    let dof_map: Vec<_> = fixed
        .iter()
        .map(|fixed| {
            if *fixed {
                None
            } else {
                let i = free;
                free += 1;
                Some(i)
            }
        })
        .collect();
    if free == 0 {
        return Err(PlateError::BadConstraint {
            what: "no free shell coordinates",
        });
    }
    let mut k = Coo::new(free, free);
    let mut m = Coo::new(free, free);
    let mut drill = vec![BTreeMap::<usize, f64>::new(); mesh.nodes.len()];
    let mut node_area = vec![0.0; mesh.nodes.len()];
    let mut drill_scale = vec![0.0; mesh.nodes.len()];
    for (e, tri) in mesh.tris.iter().enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let g = mesh.facet(e)?;
        let (ke, me) = element::matrices(
            mesh,
            tri,
            &g,
            tri.map(|i| directors.vectors[i]),
            &sections[e],
            mass,
            &mut cancelled,
        )?;
        for i in 0..18 {
            let ni = tri[i / 6];
            let bi = i % 6;
            let Some(row) = dof_map[6 * ni + bi] else {
                continue;
            };
            for j in 0..18 {
                let nj = tri[j / 6];
                let bj = j % 6;
                let Some(col) = dof_map[6 * nj + bj] else {
                    continue;
                };
                let mut kv = 0.0;
                let mut mv = 0.0;
                for p in 0..3 {
                    for q in 0..3 {
                        let li = 6 * (i / 6) + 3 * (bi / 3) + p;
                        let lj = 6 * (j / 6) + 3 * (bj / 3) + q;
                        let rotation = frames[ni][bi % 3][p] * frames[nj][bj % 3][q];
                        kv += rotation * ke[li][lj];
                        mv += rotation * me[li][lj];
                    }
                }
                if !kv.is_finite() || !mv.is_finite() {
                    return Err(PlateError::BadSection {
                        what: "nonfinite MITC assembly entry",
                    });
                }
                if kv != 0.0 {
                    k.push(row, col, kv);
                }
                if mv != 0.0 {
                    m.push(row, col, mv);
                }
            }
        }
        // One relative-spin equation per physical node, area averaged over
        // its fan. Physical stiffness is exactly blind to each node's director
        // drilling rotation, so this fixes only that gauge; averaging avoids
        // multiple inconsistent drilling penalties on a curved fan.
        let a66 = sections[e].matrices().0[8];
        for &node in tri {
            let weight = g.area_m2 / 3.0;
            node_area[node] += weight;
            drill_scale[node] += weight * a66;
            for c in 0..3 {
                *drill[node].entry(6 * node + 3 + c).or_default() += weight * g.frame[2][c];
            }
            for (a, &source) in tri.iter().enumerate() {
                for c in 0..3 {
                    let spin =
                        0.5 * (g.gradient[a][0] * g.frame[1][c] - g.gradient[a][1] * g.frame[0][c]);
                    *drill[node].entry(6 * source + c).or_default() -= weight * spin;
                }
            }
        }
    }
    for node in 0..mesh.nodes.len() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let mut row = BTreeMap::<usize, f64>::new();
        for (&full, &value) in &drill[node] {
            let n = full / 6;
            let block = 3 * ((full % 6) / 3);
            let global = full % 3;
            for a in 0..3 {
                if let Some(i) = dof_map[6 * n + block + a] {
                    *row.entry(i).or_default() += frames[n][a][global] * value / node_area[node];
                }
            }
        }
        let scale = drilling_alpha * drill_scale[node];
        for (&i, &a) in &row {
            for (&j, &b) in &row {
                k.push(i, j, scale * a * b);
            }
        }
    }
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    let pencil = ShellModel {
        k: k.assemble(),
        m: m.assemble(),
        dof_map,
        free,
    };
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    Ok(MitcShellModel {
        pencil,
        frames,
        mass_model: mass,
    })
}

/// Positive-definite static solve with physical world loads and output. A
/// free body, an unstable stiffness or a nonfinite/large true residual refuses.
pub fn solve_static(model: &MitcShellModel, full_load_6n: &[f64]) -> Result<Vec<f64>, PlateError> {
    let rhs = model.project_load(full_load_6n)?;
    let symbolic = SymbolicLdlt::analyze(&model.pencil.k, DirectOrdering::Amd)
        .map_err(PlateError::StaticFactor)?;
    let factor = symbolic
        .factor(&model.pencil.k, &LdltOptions::default())
        .map_err(PlateError::StaticFactor)?;
    if factor.inertia().negative != 0 || factor.inertia().positive != model.pencil.free {
        return Err(PlateError::StaticSolve {
            what: "static shell stiffness is not positive definite",
        });
    }
    let values = factor.solve(&rhs);
    let mut ku = vec![0.0; rhs.len()];
    model.pencil.k.spmv(&values, &mut ku);
    let load = rhs.iter().fold(0.0_f64, |r, v| r.hypot(*v));
    let residual = ku
        .iter()
        .zip(&rhs)
        .fold(0.0_f64, |r, (a, b)| r.hypot(a - b));
    if !residual.is_finite() || residual > 1e-8 * load.max(f64::MIN_POSITIVE) {
        return Err(PlateError::StaticSolve {
            what: "static shell true residual exceeds 1e-8 relative",
        });
    }
    model.expand(&values)
}

#[cfg(test)]
mod tests;
