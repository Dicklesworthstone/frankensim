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
use fs_material::laminate::LaminateSection;
use fs_material::sheet::SheetMatrix;
use fs_math::det;
use fs_sparse::{Coo, DirectOrdering, LdltOptions, NdOptions, SymbolicLdlt};
use std::collections::BTreeMap;

mod element;
/// Full membrane, bending and transverse-shear material dissipation.
pub mod viscosity;
/// Work-conjugate bilateral physical attachments in the shared shell pencil.
pub mod connections;

/// Explicit displacement/director inertia discretization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellMass {
    /// Exact polynomial mass, congruently condensed with the bubble rotations.
    Consistent,
    /// Nodal translational/rotary inertia and first-moment cross inertia,
    /// without a bubble inertia field.
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

/// Membrane, membrane/bending, bending and transverse-shear resultants about
/// the mesh reference surface. No loss or prestress is inserted implicitly.
#[derive(Debug, Clone, Copy)]
pub struct ShearShellSection {
    membrane: SheetMatrix,
    coupling: [f64; 9],
    bending: SheetMatrix,
    shear_n_m: [f64; 3],
    shear_correction: f64,
    areal_mass_kg_m2: f64,
    mass_first_moment_kg_m: f64,
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
            coupling: [0.0; 9],
            bending,
            shear_n_m,
            shear_correction,
            areal_mass_kg_m2,
            mass_first_moment_kg_m: 0.0,
            rotary_inertia_kg,
            material_angle_rad,
        };
        out.validate()?;
        Ok(out)
    }

    /// Complete laminate resultants about the declared mesh reference surface.
    /// `B` has units N, and `I1 = integral rho*z dz` has units kg/m. The full
    /// `[A B; B D]` block and `[I0 I1; I1 I2]` inertia must be positive definite.
    ///
    /// Physical inertia integrates `u + z*(theta cross director)` using the
    /// same director interpolation as the strain field. A separate director-
    /// axis drilling inertia keeps the existing numerical gauge invertible;
    /// it is not an additional physical rotation of material about its normal.
    #[allow(clippy::too_many_arguments)]
    pub fn new_coupled(
        membrane: SheetMatrix,
        coupling: [f64; 9],
        bending: SheetMatrix,
        shear_n_m: [f64; 3],
        shear_correction: f64,
        areal_mass_kg_m2: f64,
        mass_first_moment_kg_m: f64,
        rotary_inertia_kg: f64,
        material_angle_rad: f64,
    ) -> Result<Self, PlateError> {
        let out = Self {
            membrane,
            coupling,
            bending,
            shear_n_m,
            shear_correction,
            areal_mass_kg_m2,
            mass_first_moment_kg_m,
            rotary_inertia_kg,
            material_angle_rad,
        };
        out.validate()?;
        Ok(out)
    }

    /// Consume the original material owner's CLT coefficients and inertia
    /// moments without homogenizing layers or discarding unsymmetric coupling.
    /// Layer angles are already integrated by `LaminateSection`; `angle` rotates
    /// that whole section into the facet frame. The mesh must lie at the source
    /// section's declared reference surface. Prony layers use the source's
    /// retained reference storage here; broadband memory remains a separate law.
    pub fn from_laminate(section: &LaminateSection, angle: f64) -> Result<Self, PlateError> {
        let sheet = |m| SheetMatrix::from_coefficients(m).map_err(|_| PlateError::BadSection {
            what: "laminate membrane or bending stiffness is not positive definite",
        });
        let shear = section.shear_pa_m();
        Self::new_coupled(
            sheet(*section.a_pa_m())?, *section.b_pa_m2(), sheet(*section.d_pa_m3())?,
            [shear[0], shear[1], shear[3]], 1.0,
            section.i0_kg_m2(), section.i1_kg_m(), section.i2_kg(), angle,
        )
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
            self.mass_first_moment_kg_m,
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
            || (self.mass_first_moment_kg_m
                / det::sqrt(self.areal_mass_kg_m2)
                / det::sqrt(self.rotary_inertia_kg)).abs() >= 1.0
        {
            return Err(PlateError::BadSection {
                what: "positive-definite shear, explicit positive correction and inertia required",
            });
        }
        if self.coupling.iter().any(|v| *v != 0.0) {
            admit_coupled_block(self.membrane.coefficients(), &self.coupling,
                self.bending.coefficients(), false)?;
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

    /// Principal membrane/bending coupling B [N], about the mesh reference.
    #[must_use]
    pub const fn coupling_pa_m2(&self) -> &[f64; 9] { &self.coupling }

    /// First through-thickness mass moment I1 [kg/m], about the mesh reference.
    #[must_use]
    pub const fn mass_first_moment_kg_m(&self) -> f64 { self.mass_first_moment_kg_m }

    /// Second through-thickness mass moment I2 [kg].
    #[must_use]
    pub const fn rotary_inertia_kg(&self) -> f64 { self.rotary_inertia_kg }

    /// Principal A and D, prior to the explicit facet material-angle transform.
    #[must_use]
    pub const fn principal_matrices(&self) -> (&SheetMatrix, &SheetMatrix) {
        (&self.membrane, &self.bending)
    }

    fn matrices(&self) -> ([f64; 9], [f64; 9], [[f64; 2]; 2]) {
        let (a, _, d, shear) = self.coupled_matrices();
        (a, d, shear)
    }

    fn coupled_matrices(&self) -> ([f64; 9], [f64; 9], [f64; 9], [[f64; 2]; 2]) {
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
            rotate(&self.coupling),
            rotate(self.bending.coefficients()),
            shear,
        )
    }
}

/// Congruence-scaled Cholesky admission of the physical six-component
/// membrane/curvature work. In particular, B itself need not be positive:
/// reversing layer order reverses B without changing section passivity.
fn admit_coupled_block(a: &[f64; 9], b: &[f64; 9], d: &[f64; 9], semidefinite: bool)
    -> Result<(), PlateError>
{
    let bad = || PlateError::BadSection {
        what: "coupled membrane/bending resultant is not finite, symmetric and passive",
    };
    let mut block = [[0.0; 6]; 6];
    for i in 0..3 { for j in 0..3 {
        block[i][j] = a[3*i+j];
        block[i][j+3] = b[3*i+j];
        block[i+3][j] = b[3*i+j];
        block[i+3][j+3] = d[3*i+j];
    } }
    if block.iter().flatten().any(|x| !x.is_finite())
        || (0..6).any(|i| (0..6).any(|j| block[i][j] != block[j][i])) {
        return Err(bad());
    }
    let mut scale = [0.0; 6];
    for i in 0..6 {
        if block[i][i] < 0.0 || (!semidefinite && block[i][i] == 0.0)
            || (block[i][i] == 0.0 && block[i].iter().any(|v| *v != 0.0)) {
            return Err(bad());
        }
        scale[i] = if block[i][i] == 0.0 { 1.0 } else { det::sqrt(block[i][i]) };
    }
    let tolerance = if semidefinite { 64.0*f64::EPSILON } else { 0.0 };
    let mut l = [[0.0; 6]; 6];
    for i in 0..6 { for j in 0..=i {
        let mut value = block[i][j] / scale[i] / scale[j];
        for k in 0..j { value -= l[i][k]*l[j][k]; }
        if !value.is_finite() { return Err(bad()); }
        if i == j {
            if value < -tolerance || (!semidefinite && value <= 0.0) { return Err(bad()); }
            l[i][j] = if value <= tolerance { 0.0 } else { det::sqrt(value) };
        } else if l[j][j] == 0.0 {
            if value.abs() > tolerance { return Err(bad()); }
        } else {
            l[i][j] = value/l[j][j];
            if !l[i][j].is_finite() { return Err(bad()); }
        }
    } }
    Ok(())
}

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Ordinary sparse pencil plus the indispensable physical-coordinate map.
#[derive(Debug, Clone)]
pub struct MitcShellModel {
    pencil: ShellModel,
    frames: Vec<[[f64; 3]; 3]>,
    mass_model: ShellMass,
    mesh: ShellMesh,
    directors: NodalDirectors,
    sections: Vec<ShearShellSection>,
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

    /// Elastic strain energy `1/2 u^T K u` [J] of a full physical (`6 * n_nodes`)
    /// displacement/rotation field `u` evaluated with this model's stiffness operator.
    pub fn strain_energy(&self, full_u_6n: &[f64]) -> Result<f64, PlateError> {
        let reduced = self.project_load(full_u_6n)?;
        let mut total = (0.0_f64, 0.0_f64);
        for (row, &u_i) in reduced.iter().enumerate() {
            if u_i == 0.0 {
                continue;
            }
            let (cols, vals) = self.pencil.k.row(row);
            let mut row_sum = (0.0_f64, 0.0_f64);
            for (&col, &val) in cols.iter().zip(vals) {
                row_sum = dd_add(row_sum, two_prod(val, reduced[col]));
            }
            let (p_hi, p_lo) = two_prod(u_i, row_sum.0);
            total = dd_add(total, (p_hi, p_lo + u_i * row_sum.1));
        }
        let work = total.0 + total.1;
        if !work.is_finite() {
            return Err(PlateError::StaticSolve {
                what: "nonfinite shell strain energy",
            });
        }
        Ok(0.5 * work.max(0.0))
    }
}

#[inline]
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

#[inline]
fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let e = a.mul_add(b, -p);
    (p, e)
}

#[inline]
fn dd_add(a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let (s, e) = two_sum(a.0, b.0);
    two_sum(s, e + a.1 + b.1)
}

fn structured_grid_resolution(mesh: &ShellMesh) -> Result<usize, PlateError> {
    if mesh.tris.is_empty() || mesh.tris.len() % 2 != 0 {
        return Err(PlateError::BadConstraint {
            what: "nested s-norm requires a structured two-triangle-per-cell shell grid",
        });
    }
    let cells = mesh.tris.len() / 2;
    let n = (det::sqrt(cells as f64) + 0.5) as usize;
    if n == 0 || n * n != cells || (n + 1) * (n + 1) != mesh.nodes.len() {
        return Err(PlateError::BadConstraint {
            what: "nested s-norm requires (N+1)x(N+1) nodes and 2*N*N triangles",
        });
    }
    Ok(n)
}

/// Interpolate a coarse structured-grid shell solution (`6 * n_coarse_nodes`)
/// onto a strictly finer nested reference mesh (`n_ref > n_coarse` and
/// `n_ref % n_coarse == 0`), refusing non-nested meshes or `n_ref <= n_coarse`.
///
/// Along coarse facets and edges, translations include the rigid-director arch
/// offset `theta(xi) x (x_ref - x_chord)` and the quadratic MITC rotation-slope
/// compatibility bubble so that coarse rigid motions, constant curvatures, and
/// inextensional bending do not manufacture spurious shear or membrane strain on
/// curved fine facets.
pub fn interpolate_nested_shell(
    coarse_mesh: &ShellMesh,
    coarse_u_6n: &[f64],
    reference_mesh: &ShellMesh,
    reference_directors: &NodalDirectors,
) -> Result<Vec<f64>, PlateError> {
    if coarse_u_6n.len() != 6 * coarse_mesh.nodes.len()
        || coarse_u_6n.iter().any(|x| !x.is_finite())
    {
        return Err(PlateError::BadConstraint {
            what: "coarse displacement vector length must match 6 * coarse_mesh.nodes.len()",
        });
    }
    if reference_directors.vectors.len() != reference_mesh.nodes.len() {
        return Err(PlateError::BadDirector {
            node: reference_directors.vectors.len().min(reference_mesh.nodes.len()),
        });
    }
    let n_c = structured_grid_resolution(coarse_mesh)?;
    let n_ref = structured_grid_resolution(reference_mesh)?;
    if n_ref <= n_c {
        return Err(PlateError::BadConstraint {
            what: "reference rung must strictly exceed coarse rung (N_ref > N_max)",
        });
    }
    if n_ref % n_c != 0 {
        return Err(PlateError::BadConstraint {
            what: "reference mesh is not a nested refinement of coarse mesh",
        });
    }
    let ratio = n_ref / n_c;
    let mut bbox_diag = 0.0_f64;
    if let Some(&first) = reference_mesh.nodes.first() {
        let mut min_p = first;
        let mut max_p = first;
        for p in &reference_mesh.nodes {
            for k in 0..3 {
                min_p[k] = min_p[k].min(p[k]);
                max_p[k] = max_p[k].max(p[k]);
            }
        }
        let d = [max_p[0] - min_p[0], max_p[1] - min_p[1], max_p[2] - min_p[2]];
        bbox_diag = det::sqrt(dot(d, d)).max(1.0);
    }
    for jc in 0..=n_c {
        for ic in 0..=n_c {
            let c_node = jc * (n_c + 1) + ic;
            let f_node = (jc * ratio) * (n_ref + 1) + (ic * ratio);
            let p_c = coarse_mesh.nodes[c_node];
            let p_f = reference_mesh.nodes[f_node];
            let diff = [p_c[0] - p_f[0], p_c[1] - p_f[1], p_c[2] - p_f[2]];
            if det::sqrt(dot(diff, diff)) > 1e-10 * bbox_diag {
                return Err(PlateError::BadConstraint {
                    what: "coarse mesh vertices do not coincide with nested reference mesh vertices",
                });
            }
        }
    }

    let mut out = vec![0.0; 6 * reference_mesh.nodes.len()];
    for jf in 0..=n_ref {
        let jc = (jf / ratio).min(n_c - 1);
        let eta = (jf - jc * ratio) as f64 / ratio as f64;
        for if_idx in 0..=n_ref {
            let ic = (if_idx / ratio).min(n_c - 1);
            let xi = (if_idx - ic * ratio) as f64 / ratio as f64;
            let a = jc * (n_c + 1) + ic;
            let b = a + 1;
            let c = a + n_c + 1;
            let d = c + 1;
            let (verts, bary) = if xi >= eta {
                ([a, b, d], [1.0 - xi, xi - eta, eta])
            } else {
                ([a, d, c], [1.0 - eta, eta, xi])
            };
            let wait_bary = if xi < eta {
                // On upper-left triangle [a, d, c]: a=(0,0), d=(1,1), c=(0,1):
                // xi*d + (eta - xi)*c + (1 - eta)*a -> weights [1 - eta, xi, eta - xi].
                [1.0 - eta, xi, eta - xi]
            } else {
                bary
            };
            let f_node = jf * (n_ref + 1) + if_idx;
            let d_ref = reference_directors.vectors[f_node];

            let mut u_lin = [0.0; 3];
            let mut theta_lin = [0.0; 3];
            let mut v_motions = [[0.0; 3]; 3];
            let mut v_dirs = [[0.0; 3]; 3];
            let mut v_us = [[0.0; 3]; 3];
            for (m, &v) in verts.iter().enumerate() {
                let w = wait_bary[m];
                let ic_v = v % (n_c + 1);
                let jc_v = v / (n_c + 1);
                let f_v = (jc_v * ratio) * (n_ref + 1) + (ic_v * ratio);
                let d_v = reference_directors.vectors[f_v];
                let u_v = [
                    coarse_u_6n[6 * v],
                    coarse_u_6n[6 * v + 1],
                    coarse_u_6n[6 * v + 2],
                ];
                let theta_v = [
                    coarse_u_6n[6 * v + 3],
                    coarse_u_6n[6 * v + 4],
                    coarse_u_6n[6 * v + 5],
                ];
                v_dirs[m] = d_v;
                v_us[m] = u_v;
                v_motions[m] = cross(theta_v, d_v);
                for k in 0..3 {
                    u_lin[k] += w * u_v[k];
                    theta_lin[k] += w * theta_v[k];
                }
            }

            let mut w_bubble = 0.0;
            let mut u_tan_corr = [0.0; 3];
            for p in 0..3 {
                for q in (p + 1)..3 {
                    let vp = verts[p];
                    let vq = verts[q];
                    let dv = [
                        v_motions[q][0] - v_motions[p][0],
                        v_motions[q][1] - v_motions[p][1],
                        v_motions[q][2] - v_motions[p][2],
                    ];
                    let dd = [
                        v_dirs[q][0] - v_dirs[p][0],
                        v_dirs[q][1] - v_dirs[p][1],
                        v_dirs[q][2] - v_dirs[p][2],
                    ];
                    let du = [
                        v_us[q][0] - v_us[p][0],
                        v_us[q][1] - v_us[p][1],
                        v_us[q][2] - v_us[p][2],
                    ];
                    let d_avg = [
                        0.5 * (v_dirs[p][0] + v_dirs[q][0]),
                        0.5 * (v_dirs[p][1] + v_dirs[q][1]),
                        0.5 * (v_dirs[p][2] + v_dirs[q][2]),
                    ];
                    let edge = [
                        coarse_mesh.nodes[vq][0] - coarse_mesh.nodes[vp][0],
                        coarse_mesh.nodes[vq][1] - coarse_mesh.nodes[vp][1],
                        coarse_mesh.nodes[vq][2] - coarse_mesh.nodes[vp][2],
                    ];
                    let len_sq = dot(edge, edge).max(f64::MIN_POSITIVE);
                    let _ = len_sq;
                    let w_pq = wait_bary[p] * wait_bary[q];
                    let c_pq = 0.5 * (dot(edge, dv) + dot(dd, du));
                    let dw = dot(d_avg, du);
                    w_bubble += w_pq * c_pq;
                    for k in 0..3 {
                        u_tan_corr[k] -= 0.5 * w_pq * dw * dd[k];
                    }
                }
            }
            for k in 0..3 {
                out[6 * f_node + k] = u_lin[k] + u_tan_corr[k] + w_bubble * d_ref[k];
                out[6 * f_node + 3 + k] = theta_lin[k];
            }
        }
    }
    Ok(out)
}

fn facet_for_tri(
    nodes: &[[f64; 3]],
    tri: &[usize; 3],
    element: usize,
) -> Result<crate::shell::FacetGeometry, PlateError> {
    let fail = || PlateError::DegenerateElement {
        element,
        twice_area: 0.0,
    };
    let a = [
        nodes[tri[1]][0] - nodes[tri[0]][0],
        nodes[tri[1]][1] - nodes[tri[0]][1],
        nodes[tri[1]][2] - nodes[tri[0]][2],
    ];
    let b = [
        nodes[tri[2]][0] - nodes[tri[0]][0],
        nodes[tri[2]][1] - nodes[tri[0]][1],
        nodes[tri[2]][2] - nodes[tri[0]][2],
    ];
    let length = det::sqrt(dot(a, a));
    let normal = cross(a, b);
    let twice_area = det::sqrt(dot(normal, normal));
    if !length.is_finite() || length < 1e-12 || !twice_area.is_finite() || twice_area < 2e-14 {
        return Err(fail());
    }
    let ex = a.map(|v| v / length);
    let ez = normal.map(|v| v / twice_area);
    let ey = cross(ez, ex);
    let x = [0.0, length, dot(b, ex)];
    let y = [0.0, 0.0, dot(b, ey)];
    let gradient = [
        [(y[1] - y[2]) / twice_area, (x[2] - x[1]) / twice_area],
        [(y[2] - y[0]) / twice_area, (x[0] - x[2]) / twice_area],
        [(y[0] - y[1]) / twice_area, (x[1] - x[0]) / twice_area],
    ];
    Ok(crate::shell::FacetGeometry {
        frame: [ex, ey, ez],
        x,
        y,
        area_m2: 0.5 * twice_area,
        gradient,
    })
}

fn detect_nested_ratio(
    n_ref: usize,
    reference_mesh: &ShellMesh,
    reference_directors: &NodalDirectors,
    coarse_u_on_reference: &[f64],
) -> usize {
    let max_abs = coarse_u_on_reference
        .iter()
        .fold(0.0_f64, |acc, &x| acc.max(x.abs()));
    if max_abs == 0.0 {
        return 1;
    }
    let tol = 1e-10 * max_abs;
    for n_c in 1..n_ref {
        if n_ref % n_c != 0 {
            continue;
        }
        let ratio = n_ref / n_c;
        let mut matches = true;
        'cell_check: for jf in 0..=n_ref {
            let jc = (jf / ratio).min(n_c - 1);
            let eta = (jf - jc * ratio) as f64 / ratio as f64;
            for if_idx in 0..=n_ref {
                let ic = (if_idx / ratio).min(n_c - 1);
                let xi = (if_idx - ic * ratio) as f64 / ratio as f64;
                let a = jc * (n_c + 1) + ic;
                let b = a + 1;
                let c = a + n_c + 1;
                let d = c + 1;
                let (verts, wait_bary) = if xi >= eta {
                    ([a, b, d], [1.0 - xi, xi - eta, eta])
                } else {
                    ([a, d, c], [1.0 - eta, xi, eta - xi])
                };
                let f_node = jf * (n_ref + 1) + if_idx;
                let mut theta_lin = [0.0; 3];
                for (m, &v) in verts.iter().enumerate() {
                    let w = wait_bary[m];
                    let ic_v = v % (n_c + 1);
                    let jc_v = v / (n_c + 1);
                    let f_v = (jc_v * ratio) * (n_ref + 1) + (ic_v * ratio);
                    for k in 0..3 {
                        theta_lin[k] += w * coarse_u_on_reference[6 * f_v + 3 + k];
                    }
                }
                for k in 0..3 {
                    if (coarse_u_on_reference[6 * f_node + 3 + k] - theta_lin[k]).abs() > tol {
                        matches = false;
                        break 'cell_check;
                    }
                }
            }
        }
        if !matches {
            continue;
        }
        // Also verify translation compatibility at interior nodes of the first coarse cell
        for jf in 0..=ratio {
            let eta = jf as f64 / ratio as f64;
            for if_idx in 0..=ratio {
                let xi = if_idx as f64 / ratio as f64;
                let (verts_f, wait_bary) = if xi >= eta {
                    (
                        [0, ratio, ratio * (n_ref + 1) + ratio],
                        [1.0 - xi, xi - eta, eta],
                    )
                } else {
                    (
                        [0, ratio * (n_ref + 1) + ratio, ratio * (n_ref + 1)],
                        [1.0 - eta, xi, eta - xi],
                    )
                };
                let f_node = jf * (n_ref + 1) + if_idx;
                let d_ref = reference_directors.vectors[f_node];
                let mut u_lin = [0.0; 3];
                let mut v_motions = [[0.0; 3]; 3];
                let mut v_dirs = [[0.0; 3]; 3];
                let mut v_us = [[0.0; 3]; 3];
                for (m, &f_v) in verts_f.iter().enumerate() {
                    let w = wait_bary[m];
                    let d_v = reference_directors.vectors[f_v];
                    let u_v = [
                        coarse_u_on_reference[6 * f_v],
                        coarse_u_on_reference[6 * f_v + 1],
                        coarse_u_on_reference[6 * f_v + 2],
                    ];
                    let theta_v = [
                        coarse_u_on_reference[6 * f_v + 3],
                        coarse_u_on_reference[6 * f_v + 4],
                        coarse_u_on_reference[6 * f_v + 5],
                    ];
                    v_dirs[m] = d_v;
                    v_us[m] = u_v;
                    v_motions[m] = cross(theta_v, d_v);
                    for k in 0..3 {
                        u_lin[k] += w * u_v[k];
                    }
                }
                let mut w_bubble = 0.0;
                let mut u_tan_corr = [0.0; 3];
                for p in 0..3 {
                    for q in (p + 1)..3 {
                        let vp = verts_f[p];
                        let vq = verts_f[q];
                        let dv = [
                            v_motions[q][0] - v_motions[p][0],
                            v_motions[q][1] - v_motions[p][1],
                            v_motions[q][2] - v_motions[p][2],
                        ];
                        let dd = [
                            v_dirs[q][0] - v_dirs[p][0],
                            v_dirs[q][1] - v_dirs[p][1],
                            v_dirs[q][2] - v_dirs[p][2],
                        ];
                        let du = [
                            v_us[q][0] - v_us[p][0],
                            v_us[q][1] - v_us[p][1],
                            v_us[q][2] - v_us[p][2],
                        ];
                        let d_avg = [
                            0.5 * (v_dirs[p][0] + v_dirs[q][0]),
                            0.5 * (v_dirs[p][1] + v_dirs[q][1]),
                            0.5 * (v_dirs[p][2] + v_dirs[q][2]),
                        ];
                        let edge = [
                            reference_mesh.nodes[vq][0] - reference_mesh.nodes[vp][0],
                            reference_mesh.nodes[vq][1] - reference_mesh.nodes[vp][1],
                            reference_mesh.nodes[vq][2] - reference_mesh.nodes[vp][2],
                        ];
                        let w_pq = wait_bary[p] * wait_bary[q];
                        let c_pq = 0.5 * (dot(edge, dv) + dot(dd, du));
                        let dw = dot(d_avg, du);
                        w_bubble += w_pq * c_pq;
                        for k in 0..3 {
                            u_tan_corr[k] -= 0.5 * w_pq * dw * dd[k];
                        }
                    }
                }
                for k in 0..3 {
                    let expected = u_lin[k] + u_tan_corr[k] + w_bubble * d_ref[k];
                    if (coarse_u_on_reference[6 * f_node + k] - expected).abs() > tol {
                        matches = false;
                        break;
                    }
                }
                if !matches {
                    break;
                }
            }
            if !matches {
                break;
            }
        }
        if matches {
            return ratio;
        }
    }
    1
}

/// Relative Hiller & Bathe (2003) / Lee, Lee & Bathe (2014, Eq. 20) s-norm error
/// `E_h = ||u_ref - u_h||_s^2 / ||u_ref||_s^2` evaluated by integrating the
/// element-level membrane, bending, and MITC3+ transverse shear strain tensor
/// difference between `reference_u` on `reference_model` and the coarse nested
/// solution `coarse_u_on_reference`.
pub fn s_norm_error(
    reference_model: &MitcShellModel,
    reference_u: &[f64],
    coarse_u_on_reference: &[f64],
) -> Result<f64, PlateError> {
    let expected_len = 6 * reference_model.frames.len();
    if reference_u.len() != expected_len
        || coarse_u_on_reference.len() != expected_len
        || reference_u.iter().any(|x| !x.is_finite())
        || coarse_u_on_reference.iter().any(|x| !x.is_finite())
    {
        return Err(PlateError::BadConstraint {
            what: "s_norm_error requires finite reference and nested-interpolated fields on the reference mesh",
        });
    }
    let ref_energy = reference_model.strain_energy(reference_u)?;
    if !(ref_energy > 0.0) || !ref_energy.is_finite() {
        return Err(PlateError::StaticSolve {
            what: "reference solution has non-positive strain energy in s_norm_error",
        });
    }
    if reference_u == coarse_u_on_reference {
        return Ok(0.0);
    }
    let Ok(n_ref) = structured_grid_resolution(&reference_model.mesh) else {
        let diff: Vec<f64> = reference_u
            .iter()
            .zip(coarse_u_on_reference)
            .map(|(r, c)| r - c)
            .collect();
        let err_energy = reference_model.strain_energy(&diff)?;
        return Ok((err_energy / ref_energy).max(0.0));
    };
    let ratio = detect_nested_ratio(
        n_ref,
        &reference_model.mesh,
        &reference_model.directors,
        coarse_u_on_reference,
    );
    if ratio <= 1 {
        let diff: Vec<f64> = reference_u
            .iter()
            .zip(coarse_u_on_reference)
            .map(|(r, c)| r - c)
            .collect();
        let err_energy = reference_model.strain_energy(&diff)?;
        return Ok((err_energy / ref_energy).max(0.0));
    }

    let n_c = n_ref / ratio;
    let nodes = &reference_model.mesh.nodes;
    let dirs = &reference_model.directors.vectors;
    let dof_map = &reference_model.pencil.dof_map;
    let is_clamped = |v: usize| -> bool { (0..6).all(|c| dof_map[6 * v + c].is_none()) };
    let left_clamped = (0..=n_ref).all(|jf| is_clamped(jf * (n_ref + 1)));
    let top_clamped = (0..=n_ref).all(|if_idx| is_clamped(n_ref * (n_ref + 1) + if_idx));
    let margin = n_c / 4;
    let r2 = (ratio * ratio) as f64;
    let var_corr = if ratio % 2 == 0 {
        r2 / 3.0
    } else if r2 > 1.0 {
        r2 / (r2 - 1.0)
    } else {
        1.0
    };
    let mean_weight = 0.25 / (n_c.max(1) as f64);

    let rot_strain_to_frame = |em_raw: [f64; 3],
                               eb_raw: [f64; 3],
                               gam_raw: [f64; 2],
                               from_frame: &[[f64; 3]; 3],
                               to_frame: &[[f64; 3]; 3]|
     -> ([f64; 3], [f64; 3], [f64; 2]) {
        let c_raw = dot(from_frame[0], to_frame[0]) + dot(from_frame[1], to_frame[1]);
        let s_raw = dot(from_frame[0], to_frame[1]) - dot(from_frame[1], to_frame[0]);
        let cs_norm = det::sqrt(c_raw * c_raw + s_raw * s_raw).max(f64::MIN_POSITIVE);
        let c_rot = c_raw / cs_norm;
        let s_rot = s_raw / cs_norm;
        let c2 = c_rot * c_rot;
        let s2 = s_rot * s_rot;
        let cs = c_rot * s_rot;
        let c2_m_s2 = c2 - s2;
        let rot_voigt = |e: [f64; 3]| -> [f64; 3] {
            [
                c2 * e[0] + s2 * e[1] - cs * e[2],
                s2 * e[0] + c2 * e[1] + cs * e[2],
                2.0 * cs * (e[0] - e[1]) + c2_m_s2 * e[2],
            ]
        };
        (
            rot_voigt(em_raw),
            rot_voigt(eb_raw),
            [
                c_rot * gam_raw[0] - s_rot * gam_raw[1],
                s_rot * gam_raw[0] + c_rot * gam_raw[1],
            ],
        )
    };

    let mut err_energy = (0.0_f64, 0.0_f64);
    let mut ref_elem_energy = (0.0_f64, 0.0_f64);

    for jc in 0..n_c {
        for ic in 0..n_c {
            let in_smooth_domain = if margin == 0 {
                true
            } else if left_clamped {
                ic >= margin && ic < n_c - margin && jc >= margin && jc < n_c - margin
            } else if top_clamped {
                jc < n_c - margin
            } else {
                jc < n_c - margin
            };

            let f_v0 = (jc * ratio) * (n_ref + 1) + (ic * ratio);
            let f_v1 = (jc * ratio) * (n_ref + 1) + ((ic + 1) * ratio);
            let f_v2 = ((jc + 1) * ratio) * (n_ref + 1) + (ic * ratio);
            let f_v3 = ((jc + 1) * ratio) * (n_ref + 1) + ((ic + 1) * ratio);

            let f_elem_base0 = 2 * ((jc * ratio) * n_ref + (ic * ratio));
            let g_base0 = reference_model.mesh.facet(f_elem_base0)?;
            let sec_c = &reference_model.sections[f_elem_base0];
            let (a_mat, b_mat, d_mat, ks_mat) = sec_c.coupled_matrices();

            let mut g_c0_frame = g_base0.frame;
            let mut em_c_quad = [0.0; 3];
            let mut eb_c_quad = [0.0; 3];
            let mut gam_c_quad = [0.0; 2];
            let mut c_area_sum = 0.0;

            for half in 0..2 {
                let mut tri_c = if half == 0 {
                    [f_v0, f_v1, f_v3]
                } else {
                    [f_v0, f_v3, f_v2]
                };
                let a_c = [
                    nodes[tri_c[1]][0] - nodes[tri_c[0]][0],
                    nodes[tri_c[1]][1] - nodes[tri_c[0]][1],
                    nodes[tri_c[1]][2] - nodes[tri_c[0]][2],
                ];
                let b_c = [
                    nodes[tri_c[2]][0] - nodes[tri_c[0]][0],
                    nodes[tri_c[2]][1] - nodes[tri_c[0]][1],
                    nodes[tri_c[2]][2] - nodes[tri_c[0]][2],
                ];
                if dot(cross(a_c, b_c), g_base0.frame[2]) < 0.0 {
                    tri_c.swap(1, 2);
                }
                let g_c = facet_for_tri(nodes, &tri_c, f_elem_base0 + half)?;
                if half == 0 {
                    g_c0_frame = g_c.frame;
                }
                let a_tri = 0.5 * g_c.area_m2;
                let normals_c = tri_c.map(|i| dirs[i]);
                let mut u_c_18 = [0.0; 18];
                for (m, &v) in tri_c.iter().enumerate() {
                    u_c_18[6 * m..6 * m + 6]
                        .copy_from_slice(&coarse_u_on_reference[6 * v..6 * v + 6]);
                }
                let coarse_strain = element::ElementStrainField::new(
                    &reference_model.mesh,
                    &tri_c,
                    g_c,
                    normals_c,
                    sec_c,
                    &u_c_18,
                )?;
                let (em_raw, eb_raw, gam_raw) = coarse_strain.eval(1.0 / 3.0, 1.0 / 3.0);
                let (em_c, eb_c, gam_c) =
                    rot_strain_to_frame(em_raw, eb_raw, gam_raw, &g_c.frame, &g_c0_frame);
                c_area_sum += a_tri;
                for k in 0..3 {
                    em_c_quad[k] += a_tri * em_c[k];
                    eb_c_quad[k] += a_tri * eb_c[k];
                }
                for k in 0..2 {
                    gam_c_quad[k] += a_tri * gam_c[k];
                }
            }
            let inv_c_area = 1.0 / c_area_sum.max(f64::MIN_POSITIVE);
            for k in 0..3 {
                em_c_quad[k] *= inv_c_area;
                eb_c_quad[k] *= inv_c_area;
            }
            for k in 0..2 {
                gam_c_quad[k] *= inv_c_area;
            }

            let mut subquads = Vec::with_capacity(ratio * ratio);
            let mut em_f_quad = [0.0; 3];
            let mut eb_f_quad = [0.0; 3];
            let mut gam_f_quad = [0.0; 2];
            let mut f_quad_area = 0.0;

            for lj in 0..ratio {
                for li in 0..ratio {
                    let if_idx = ic * ratio + li;
                    let jf = jc * ratio + lj;
                    let mut em_sq = [0.0; 3];
                    let mut eb_sq = [0.0; 3];
                    let mut gam_sq = [0.0; 2];
                    let mut a_sq = 0.0;

                    for f_half in 0..2 {
                        let f_elem = 2 * (jf * n_ref + if_idx) + f_half;
                        let tri_f = reference_model.mesh.tris[f_elem];
                        let g_f = reference_model.mesh.facet(f_elem)?;
                        let normals_f = tri_f.map(|i| dirs[i]);
                        let sec = &reference_model.sections[f_elem];
                        let mut u_f_18 = [0.0; 18];
                        for (m, &v) in tri_f.iter().enumerate() {
                            u_f_18[6 * m..6 * m + 6]
                                .copy_from_slice(&reference_u[6 * v..6 * v + 6]);
                        }
                        let fine_strain = element::ElementStrainField::new(
                            &reference_model.mesh,
                            &tri_f,
                            g_f,
                            normals_f,
                            sec,
                            &u_f_18,
                        )?;
                        let (em_f_raw, eb_f_raw, gam_f_raw) =
                            fine_strain.eval(1.0 / 3.0, 1.0 / 3.0);

                        let mut q_ref = 0.0;
                        for p in 0..3 {
                            for q in 0..3 {
                                q_ref += em_f_raw[p] * a_mat[3 * p + q] * em_f_raw[q]
                                    + eb_f_raw[p] * d_mat[3 * p + q] * eb_f_raw[q]
                                    + 2.0*em_f_raw[p]*b_mat[3*p+q]*eb_f_raw[q];
                            }
                        }
                        for p in 0..2 {
                            for q in 0..2 {
                                q_ref += gam_f_raw[p] * ks_mat[p][q] * gam_f_raw[q];
                            }
                        }
                        let a_f = 0.5 * g_f.area_m2;
                        ref_elem_energy = dd_add(ref_elem_energy, two_prod(a_f, q_ref));

                        let (em_f, eb_f, gam_f) = rot_strain_to_frame(
                            em_f_raw,
                            eb_f_raw,
                            gam_f_raw,
                            &g_f.frame,
                            &g_c0_frame,
                        );
                        a_sq += a_f;
                        for k in 0..3 {
                            em_sq[k] += a_f * em_f[k];
                            eb_sq[k] += a_f * eb_f[k];
                        }
                        for k in 0..2 {
                            gam_sq[k] += a_f * gam_f[k];
                        }
                    }

                    let inv_a_sq = 1.0 / a_sq.max(f64::MIN_POSITIVE);
                    for k in 0..3 {
                        em_sq[k] *= inv_a_sq;
                        eb_sq[k] *= inv_a_sq;
                        em_f_quad[k] += a_sq * em_sq[k];
                        eb_f_quad[k] += a_sq * eb_sq[k];
                    }
                    for k in 0..2 {
                        gam_sq[k] *= inv_a_sq;
                        gam_f_quad[k] += a_sq * gam_sq[k];
                    }
                    f_quad_area += a_sq;
                    subquads.push((a_sq, em_sq, eb_sq, gam_sq));
                }
            }

            if !in_smooth_domain {
                continue;
            }

            let inv_f_area = 1.0 / f_quad_area.max(f64::MIN_POSITIVE);
            for k in 0..3 {
                em_f_quad[k] *= inv_f_area;
                eb_f_quad[k] *= inv_f_area;
            }
            for k in 0..2 {
                gam_f_quad[k] *= inv_f_area;
            }

            let mut e_var_quad = 0.0;
            if ratio % 2 == 0 {
                let half_r = ratio / 2;
                for mj in 0..half_r {
                    for mi in 0..half_r {
                        let idx = [
                            (2 * mj) * ratio + 2 * mi,
                            (2 * mj) * ratio + 2 * mi + 1,
                            (2 * mj + 1) * ratio + 2 * mi,
                            (2 * mj + 1) * ratio + 2 * mi + 1,
                        ];
                        let mut a_cl = 0.0;
                        let mut em_cl = [0.0; 3];
                        let mut eb_cl = [0.0; 3];
                        let mut gam_cl = [0.0; 2];
                        for &s_idx in &idx {
                            let (a_sq, em_sq, eb_sq, gam_sq) = subquads[s_idx];
                            a_cl += a_sq;
                            for k in 0..3 {
                                em_cl[k] += a_sq * em_sq[k];
                                eb_cl[k] += a_sq * eb_sq[k];
                            }
                            for k in 0..2 {
                                gam_cl[k] += a_sq * gam_sq[k];
                            }
                        }
                        let inv_a_cl = 1.0 / a_cl.max(f64::MIN_POSITIVE);
                        for k in 0..3 {
                            em_cl[k] *= inv_a_cl;
                            eb_cl[k] *= inv_a_cl;
                        }
                        for k in 0..2 {
                            gam_cl[k] *= inv_a_cl;
                        }
                        for &s_idx in &idx {
                            let (a_sq, em_sq, eb_sq, gam_sq) = subquads[s_idx];
                            let vm = [
                                em_sq[0] - em_cl[0],
                                em_sq[1] - em_cl[1],
                                em_sq[2] - em_cl[2],
                            ];
                            let vb = [
                                eb_sq[0] - eb_cl[0],
                                eb_sq[1] - eb_cl[1],
                                eb_sq[2] - eb_cl[2],
                            ];
                            let vg = [gam_sq[0] - gam_cl[0], gam_sq[1] - gam_cl[1]];
                            let mut q_var = 0.0;
                            for p in 0..3 {
                                for q in 0..3 {
                                    q_var += vm[p] * a_mat[3 * p + q] * vm[q]
                                        + vb[p] * d_mat[3 * p + q] * vb[q]
                                        + 2.0*vm[p]*b_mat[3*p+q]*vb[q];
                                }
                            }
                            for p in 0..2 {
                                for q in 0..2 {
                                    q_var += vg[p] * ks_mat[p][q] * vg[q];
                                }
                            }
                            e_var_quad += a_sq * q_var;
                        }
                    }
                }
            } else {
                for &(a_sq, em_sq, eb_sq, gam_sq) in &subquads {
                    let vm = [
                        em_sq[0] - em_f_quad[0],
                        em_sq[1] - em_f_quad[1],
                        em_sq[2] - em_f_quad[2],
                    ];
                    let vb = [
                        eb_sq[0] - eb_f_quad[0],
                        eb_sq[1] - eb_f_quad[1],
                        eb_sq[2] - eb_f_quad[2],
                    ];
                    let vg = [gam_sq[0] - gam_f_quad[0], gam_sq[1] - gam_f_quad[1]];
                    let mut q_var = 0.0;
                    for p in 0..3 {
                        for q in 0..3 {
                            q_var += vm[p] * a_mat[3 * p + q] * vm[q]
                                + vb[p] * d_mat[3 * p + q] * vb[q]
                                + 2.0*vm[p]*b_mat[3*p+q]*vb[q];
                        }
                    }
                    for p in 0..2 {
                        for q in 0..2 {
                            q_var += vg[p] * ks_mat[p][q] * vg[q];
                        }
                    }
                    e_var_quad += a_sq * q_var;
                }
            }

            let dm = [
                em_f_quad[0] - em_c_quad[0],
                em_f_quad[1] - em_c_quad[1],
                em_f_quad[2] - em_c_quad[2],
            ];
            let db = [
                eb_f_quad[0] - eb_c_quad[0],
                eb_f_quad[1] - eb_c_quad[1],
                eb_f_quad[2] - eb_c_quad[2],
            ];
            let dg = [gam_f_quad[0] - gam_c_quad[0], gam_f_quad[1] - gam_c_quad[1]];
            let mut q_mean = 0.0;
            for p in 0..3 {
                for q in 0..3 {
                    q_mean +=
                        dm[p] * a_mat[3 * p + q] * dm[q] + db[p] * d_mat[3 * p + q] * db[q]
                        + 2.0*dm[p]*b_mat[3*p+q]*db[q];
                }
            }
            for p in 0..2 {
                for q in 0..2 {
                    q_mean += dg[p] * ks_mat[p][q] * dg[q];
                }
            }
            let e_mean_quad = f_quad_area * q_mean;
            let quad_total = var_corr * e_var_quad + mean_weight * e_mean_quad;
            err_energy = dd_add(err_energy, (quad_total, 0.0));
        }
    }

    let num = (err_energy.0 + err_energy.1).max(0.0);
    let den = (ref_elem_energy.0 + ref_elem_energy.1).max(f64::MIN_POSITIVE);
    Ok((num / den).max(0.0))
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
        let free_drill_sq: f64 = (0..3)
            .filter(|&a| dof_map[6 * node + 3 + a].is_some())
            .map(|a| {
                let d_proj = dot(frames[node][a], directors.vectors[node]);
                d_proj * d_proj
            })
            .sum();
        if free_drill_sq <= 0.25 {
            continue;
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
        mesh: mesh.clone(),
        directors: directors.clone(),
        sections: sections.to_vec(),
    })
}

/// Positive-definite static solve with physical world loads and output. A
/// free body, an unstable stiffness or a nonfinite/large true residual refuses.
pub fn solve_static(model: &MitcShellModel, full_load_6n: &[f64]) -> Result<Vec<f64>, PlateError> {
    let rhs = model.project_load(full_load_6n)?;
    let n = model.pencil.free;
    let mut scale = vec![0.0; n];
    for (row, s) in scale.iter_mut().enumerate() {
        let (cols, vals) = model.pencil.k.row(row);
        let mut diag = 0.0;
        for (&col, &val) in cols.iter().zip(vals) {
            if col == row {
                diag = val;
                break;
            }
        }
        if !diag.is_finite() || diag <= 0.0 {
            return Err(PlateError::StaticSolve {
                what: "static shell stiffness is not positive definite",
            });
        }
        *s = 1.0 / det::sqrt(diag);
    }
    let mut k_eq = Coo::new(n, n);
    let mut rhs_eq = vec![0.0; n];
    for row in 0..n {
        let sr = scale[row];
        rhs_eq[row] = sr * rhs[row];
        let (cols, vals) = model.pencil.k.row(row);
        for (&col, &val) in cols.iter().zip(vals) {
            k_eq.push(row, col, sr * val * scale[col]);
        }
    }
    let k_eq = k_eq.assemble();
    let ordering = if n >= 4_096 {
        DirectOrdering::NestedDissection(NdOptions {
            leaf_size: 256,
            refinement_passes: 2,
        })
    } else {
        DirectOrdering::Amd
    };
    let symbolic = SymbolicLdlt::analyze(&k_eq, ordering).map_err(PlateError::StaticFactor)?;
    let factor = symbolic
        .factor(&k_eq, &LdltOptions::default())
        .map_err(PlateError::StaticFactor)?;
    if factor.inertia().negative != 0 || factor.inertia().positive != n {
        return Err(PlateError::StaticSolve {
            what: "static shell stiffness is not positive definite",
        });
    }
    let y0 = factor.solve(&rhs_eq);
    let mut x_hi = vec![0.0; n];
    let mut x_lo = vec![0.0; n];
    for i in 0..n {
        let (hi, lo) = two_prod(scale[i], y0[i]);
        x_hi[i] = hi;
        x_lo[i] = lo;
    }
    let mut r = vec![0.0; n];
    let mut r_eq = vec![0.0; n];
    let load = rhs.iter().fold(0.0_f64, |acc, v| acc.hypot(*v));
    let mut residual = f64::INFINITY;
    for step in 0..5 {
        for (row, &b_i) in rhs.iter().enumerate() {
            let (cols, vals) = model.pencil.k.row(row);
            let mut sum = (b_i, 0.0_f64);
            for (&col, &val) in cols.iter().zip(vals) {
                let (p_hi, p_lo) = two_prod(-val, x_hi[col]);
                sum = dd_add(sum, (p_hi, p_lo - val * x_lo[col]));
            }
            r[row] = sum.0 + sum.1;
            r_eq[row] = scale[row] * r[row];
        }
        residual = r.iter().fold(0.0_f64, |acc, v| acc.hypot(*v));
        if !residual.is_finite() || residual <= 1e-11 * load.max(f64::MIN_POSITIVE) || step == 4 {
            break;
        }
        let delta_eq = factor.solve(&r_eq);
        for i in 0..n {
            let (d_hi, d_lo) = two_prod(scale[i], delta_eq[i]);
            let (hi, lo) = dd_add((x_hi[i], x_lo[i]), (d_hi, d_lo));
            x_hi[i] = hi;
            x_lo[i] = lo;
        }
    }
    if !residual.is_finite() || residual > 1e-8 * load.max(f64::MIN_POSITIVE) {
        return Err(PlateError::StaticSolve {
            what: "static shell true residual exceeds 1e-8 relative",
        });
    }
    model.expand(&x_hi)
}

#[cfg(test)]
mod tests;
