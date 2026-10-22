//! Work-conjugate dual-mortar ties and eccentric junctions between shell surfaces.
//!
//! Each curved or flat shell piece is assembled on its own smooth component chart
//! (`ShellComponent`) with well-defined facet frames and nodal directors. Multi-surface
//! instrument assemblies (guitar top/back/rim/kerfing, braced soundboards, violin corpus,
//! piano rim/belly) couple components along shared 1-D edge polylines (`EdgeTie`) or
//! 2-D overlapping skin facets (`FaceTie`) without averaging directors across fold seams.
//!
//! Along a 1-D seam $\Gamma$, the Wohlmuth (2000) linear biorthogonal dual basis
//! $\psi_0(\xi) = 2 - 3\xi$, $\psi_1(\xi) = 3\xi - 1$ on $\xi \in [0, 1]$ satisfies
//! $\int_0^1 \psi_a(\xi)\,\phi_b(\xi)\,d\xi = \frac{1}{2}\delta_{ab}$ and
//! $\psi_0(\xi) + \psi_1(\xi) = 1$, making the slave mortar matrix $D$ strictly diagonal
//! ($D_{ii} = \sum_{e \ni i} L_e / 2 > 0$) and the master-to-slave transfer operator
//! $C = D^{-1} M$ locally supported without dense matrix inversion. On curved seams,
//! a minimum-norm local affine consistency projection onto the master polyline's affine
//! hull restores exact 3-D linear coordinate reproduction $\sum_j C_{ij} \mathbf{x}_{m,j} = \mathbf{x}_{s,i}$
//! (Puso 2004), ensuring zero-energy rigid-body rotations across curved non-matching seams.
//!
//! Slave degrees of freedom are eliminated by congruence transformation
//! $K_r = T^T K T$, $M_r = T^T M T$, preserving exact symmetry, positive-semidefinite
//! energy, and zero net virtual interface power $P_m + P_s = 0$.

use super::{ShellMesh, ShellModel, dot, sub};
use crate::PlateError;
use fs_math::det;
use fs_sparse::{Coo, Csr};
use std::collections::{BTreeMap, BTreeSet};

/// One smooth shell component participating in a multi-surface junction assembly.
#[derive(Debug, Clone, Copy)]
pub struct ShellComponent<'a> {
    /// Canonical component index (`0..components.len()`).
    pub id: usize,
    /// Component midsurface mesh.
    pub mesh: &'a ShellMesh,
    /// Assembled component stiffness/mass pencil (either CST/DKT or MITC3+).
    pub model: &'a ShellModel,
    /// Positive finite physical thickness [m] at each node of `mesh`.
    pub nodal_thickness_m: &'a [f64],
}

/// Ordered node sequence on one side of an [`EdgeTie`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TieSide {
    /// Index of the [`ShellComponent`] in the slice passed to [`tie_shells`].
    pub component: usize,
    /// Ordered polyline of node indices on `component`.
    pub nodes: Vec<usize>,
}

/// Through-thickness kinematic attachment offset relative to the master midsurface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TieOffset {
    /// Midsurface-to-midsurface tie ($\mathbf{d} = \mathbf{0}$).
    Midsurface,
    /// Attachment on the master upper (`side = +1`) or lower (`side = -1`) skin,
    /// shifting the slave attachment arm by $\mathbf{d} = \text{side} \cdot \frac{1}{2} h_m \mathbf{n}_m$
    /// and generating $\mathbf{u}_s = \mathbf{u}_m + \boldsymbol{\theta}_m \times \mathbf{d}$.
    Skin {
        /// Skin sign in `{-1, +1}`.
        side: i8,
    },
}

/// Rotational kinematic coupling mode across a tie interface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RotationCoupling {
    /// Couple both translations $(u, v, w)$ and physical axial rotations $(\theta_x, \theta_y, \theta_z)$.
    Full,
    /// Couple only translations $(u, v, w)$ (hinged seam; slave rotations remain independent free DOFs).
    TranslationsOnly,
    /// Couple translations and only the rotation about the seam tangent,
    /// $\boldsymbol\theta_s\cdot\mathbf t = \mathbf t\cdot\sum_j C_{ij}\boldsymbol\theta_{m,j}$.
    /// At a fold the other two rotations are each side's slope along the seam or
    /// its drilling spin. The translational ties already make the slopes
    /// continuous, and tying one side's slope to the other's drilling spin would
    /// put the drilling penalty into the junction. Edge ties only.
    SeamTangent,
}

/// Area-specific Kelvin–Voigt adhesive compliance along a 1-D glue line or 2-D glue patch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GlueLine {
    /// Area-specific normal stiffness $k_n$ [$\text{N/m}^3$].
    pub normal_n_m3: f64,
    /// Area-specific tangential shear stiffnesses $[k_{s1}, k_{s2}]$ [$\text{N/m}^3$]
    /// along the local seam tangent $\mathbf{t}_1$ and transverse in-plane direction $\mathbf{t}_2$.
    pub shear_n_m3: [f64; 2],
    /// Optional area-specific Kelvin–Voigt viscosities $[\eta_n, \eta_{s1}, \eta_{s2}]$ [$\text{N}\cdot\text{s/m}^3$].
    pub viscosity_n_s_m3: Option<[f64; 3]>,
    /// Physical glue-line width $w_g$ [m] for 1-D [`EdgeTie`]s ($dA = w_g\,ds$).
    pub width_m: f64,
}

impl GlueLine {
    /// Validate that all stiffnesses and width are finite and strictly positive, and optional
    /// Kelvin–Voigt viscosities are finite and non-negative.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if any parameter is non-finite or out of range.
    pub fn validate(&self) -> Result<(), PlateError> {
        if !self.normal_n_m3.is_finite()
            || self.normal_n_m3 <= 0.0
            || !self.shear_n_m3[0].is_finite()
            || self.shear_n_m3[0] <= 0.0
            || !self.shear_n_m3[1].is_finite()
            || self.shear_n_m3[1] <= 0.0
            || !self.width_m.is_finite()
            || self.width_m <= 0.0
        {
            return Err(PlateError::BadConstraint {
                what: "GlueLine requires finite positive normal_n_m3, shear_n_m3, and width_m",
            });
        }
        if let Some(visc) = self.viscosity_n_s_m3
            && visc.iter().any(|&v| !v.is_finite() || v < 0.0)
        {
            return Err(PlateError::BadConstraint {
                what: "GlueLine viscosity_n_s_m3 entries must be finite and non-negative",
            });
        }
        Ok(())
    }
}

/// Constitutive nature of the interface tie.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TieKind {
    /// Rigid work-conjugate master-slave elimination ($K_r = T^T K T$, $M_r = T^T M T$).
    Rigid,
    /// Compliant Kelvin–Voigt adhesive tie adding positive interface stiffness
    /// $\int_\Gamma (\mathbf{u}_a - \mathbf{u}_b)^T K_g (\mathbf{u}_a - \mathbf{u}_b)\,d\Gamma$
    /// and optional viscous damping without slave-DOF elimination.
    Compliant(GlueLine),
}

/// 1-D polyline mortar tie between master side `a` and slave side `b`.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeTie {
    /// Master component and ordered edge polyline.
    pub a: TieSide,
    /// Slave component and ordered edge polyline.
    pub b: TieSide,
    /// Midsurface or $\pm h_m / 2$ skin offset on the master surface.
    pub offset: TieOffset,
    /// Full 6-DOF vs 3-DOF translational coupling.
    pub rotations: RotationCoupling,
    /// Rigid elimination or compliant Kelvin–Voigt glue line.
    pub kind: TieKind,
}

/// 2-D overlapping facet tie (e.g., bridge plate or lap joint) between master and slave components.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceTie {
    /// Master component index.
    pub master_component: usize,
    /// Triangle indices on the master component forming the master patch.
    pub master_tris: Vec<usize>,
    /// Slave component index.
    pub slave_component: usize,
    /// Triangle indices on the slave component forming the slave patch.
    pub slave_tris: Vec<usize>,
    /// Midsurface or $\pm h_m / 2$ skin offset on the master surface.
    pub offset: TieOffset,
    /// Full 6-DOF vs 3-DOF translational coupling.
    pub rotations: RotationCoupling,
    /// Rigid elimination or compliant Kelvin–Voigt glue patch.
    pub kind: TieKind,
}

/// Explicit geometric tolerances and work ceilings for [`tie_shells`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TieBudget {
    /// Maximum mortar integration segments across all ties.
    pub max_segments: usize,
    /// Maximum eliminated slave (or compliant interface) scalar DOFs across all ties.
    pub max_slave_dofs: usize,
    /// Maximum admissible normal/interface gap [m] between slave and master geometry.
    pub gap_tolerance_m: f64,
    /// Maximum admissible tangent-frame (or facet-normal) angle mismatch [rad].
    pub angle_tolerance_rad: f64,
}

/// Auditable geometric and algebraic receipt emitted by [`tie_shells`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TieReceipt {
    /// Maximum interface gap [m] observed across all quadrature points and nodes.
    pub max_gap_m: f64,
    /// Maximum tangent/normal angle mismatch [rad] observed across all segments.
    pub max_angle_rad: f64,
    /// Total number of mortar integration segments processed.
    pub segments: usize,
    /// Total number of eliminated slave (or compliant interface) scalar DOFs.
    pub slave_dofs: usize,
}

/// One eliminated slave scalar DOF and its linear combination of retained reduced DOFs.
#[derive(Debug, Clone)]
pub struct SlaveDofCoupling {
    /// Slave component index.
    pub slave_component: usize,
    /// Slave node index on `slave_component`.
    pub slave_node: usize,
    /// Local DOF index in `0..6`.
    pub slave_dof: usize,
    /// Component-level free DOF index in `components[slave_component].model`.
    pub slave_component_free: usize,
    /// Pairs `(reduced_dof, weight)` such that $q_{\text{slave}} = \sum (w_k \, q_{r,k})$.
    pub master_terms: Vec<(usize, f64)>,
}

/// Coupled multi-component shell pencil produced by [`tie_shells`].
#[derive(Debug, Clone)]
pub struct TiedShellModel {
    /// Congruence-reduced stiffness and mass pencil $(K_r = T^T K T + K_{\text{glue}},\ M_r = T^T M T)$.
    pub model: ShellModel,
    /// Per-component `[component][node][dof_0..6]` map to retained reduced coordinate in `0..model.free`.
    /// Eliminated slave DOFs and essential-BC fixed DOFs have `None`.
    pub map: Vec<Vec<[Option<usize>; 6]>>,
    /// Auditable tie geometry and budget receipt.
    pub receipt: TieReceipt,
    /// Pure compliant interface stiffness matrix $K_{\text{glue}}$ (`model.free × model.free`).
    pub glue_stiffness: Csr,
    /// Optional compliant Kelvin–Voigt interface damping matrix $C_{\text{glue}}$ (`model.free × model.free`).
    pub glue_damping: Option<Csr>,
    /// Per-component sparse rows of $T_c$: `comp_t[c][local_free]` is a list of `(reduced_dof, weight)`.
    comp_t: Vec<Vec<Vec<(usize, f64)>>>,
    /// Per-component full-DOF (`6 * node + d`) to component-free index (`components[c].model.dof_map`).
    comp_dof_maps: Vec<Vec<Option<usize>>>,
    /// Explicit list of eliminated slave DOFs in deterministic order.
    slave_couplings: Vec<SlaveDofCoupling>,
}

impl TiedShellModel {
    /// Eliminated slave DOF coupling rows in deterministic order.
    #[must_use]
    pub fn slave_couplings(&self) -> &[SlaveDofCoupling] {
        &self.slave_couplings
    }

    /// Evaluate the compliant glue strain energy $\frac{1}{2} \mathbf{q}_r^T K_{\text{glue}} \mathbf{q}_r$ [J].
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if `reduced.len() != self.model.free` or contains non-finite entries.
    pub fn glue_strain_energy(&self, reduced: &[f64]) -> Result<f64, PlateError> {
        if reduced.len() != self.model.free || reduced.iter().any(|x| !x.is_finite()) {
            return Err(PlateError::BadConstraint {
                what: "reduced vector length must equal model.free and be finite",
            });
        }
        let mut kq = vec![0.0; self.model.free];
        self.glue_stiffness.spmv(reduced, &mut kq);
        Ok(0.5 * dot_slice(reduced, &kq))
    }

    /// Evaluate the compliant glue viscous dissipation rate $\dot{\mathbf{q}}_r^T C_{\text{glue}} \dot{\mathbf{q}}_r$ [W].
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if `reduced_velocity.len() != self.model.free` or contains non-finite entries.
    pub fn glue_dissipation_rate(&self, reduced_velocity: &[f64]) -> Result<f64, PlateError> {
        if reduced_velocity.len() != self.model.free
            || reduced_velocity.iter().any(|x| !x.is_finite())
        {
            return Err(PlateError::BadConstraint {
                what: "reduced_velocity length must equal model.free and be finite",
            });
        }
        let Some(c_glue) = &self.glue_damping else {
            return Ok(0.0);
        };
        let mut cq = vec![0.0; self.model.free];
        c_glue.spmv(reduced_velocity, &mut cq);
        Ok(dot_slice(reduced_velocity, &cq))
    }

    /// Project the compliant glue damping matrix $C_{\text{glue}}$ onto a set of reduced modal
    /// vectors $\Phi = [\boldsymbol{\phi}_0, \dots, \boldsymbol{\phi}_{m-1}]$, returning the
    /// row-major $m \times m$ symmetric positive-semidefinite matrix $C_{\text{modal}} = \Phi^T C_{\text{glue}} \Phi$.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if any mode vector length mismatches `model.free` or is non-finite.
    pub fn project_glue_damping(&self, modes: &[Vec<f64>]) -> Result<Vec<f64>, PlateError> {
        let m = modes.len();
        if m == 0
            || modes
                .iter()
                .any(|phi| phi.len() != self.model.free || phi.iter().any(|x| !x.is_finite()))
        {
            return Err(PlateError::BadConstraint {
                what: "modes must be non-empty, finite, and have length equal to model.free",
            });
        }
        let Some(c_glue) = &self.glue_damping else {
            return Ok(vec![0.0; m * m]);
        };
        let mut out = vec![0.0; m * m];
        let mut c_phi_j = vec![0.0; self.model.free];
        for j in 0..m {
            c_glue.spmv(&modes[j], &mut c_phi_j);
            for i in 0..=j {
                let val = dot_slice(&modes[i], &c_phi_j);
                out[i * m + j] = val;
                out[j * m + i] = val;
            }
        }
        Ok(out)
    }

    /// Expand a reduced coordinate vector `reduced` (`len == self.model.free`) into each
    /// component's free coordinate vector $\mathbf{q}_c = T_c \mathbf{q}_r$.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if `reduced.len() != self.model.free` or contains non-finite values.
    pub fn expand_component_free(&self, reduced: &[f64]) -> Result<Vec<Vec<f64>>, PlateError> {
        if reduced.len() != self.model.free || reduced.iter().any(|x| !x.is_finite()) {
            return Err(PlateError::BadConstraint {
                what: "reduced vector length must equal model.free and be finite",
            });
        }
        let mut out = Vec::with_capacity(self.comp_t.len());
        for comp_rows in &self.comp_t {
            let mut q_c = vec![0.0; comp_rows.len()];
            for (i, terms) in comp_rows.iter().enumerate() {
                let mut sum = 0.0;
                for &(r_idx, w) in terms {
                    sum += w * reduced[r_idx];
                }
                q_c[i] = sum;
            }
            out.push(q_c);
        }
        Ok(out)
    }

    /// Expand a reduced coordinate vector `reduced` (`len == self.model.free`) into each
    /// component's full `6 * node_count` physical displacement/rotation vector.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if `reduced.len() != self.model.free` or contains non-finite values.
    pub fn expand_components(&self, reduced: &[f64]) -> Result<Vec<Vec<f64>>, PlateError> {
        let free_vecs = self.expand_component_free(reduced)?;
        let mut out = Vec::with_capacity(free_vecs.len());
        for (c, q_free) in free_vecs.iter().enumerate() {
            let dof_map = &self.comp_dof_maps[c];
            let mut full = vec![0.0; dof_map.len()];
            for (full_idx, maybe_free) in dof_map.iter().enumerate() {
                if let Some(f_idx) = *maybe_free {
                    full[full_idx] = q_free[f_idx];
                }
            }
            out.push(full);
        }
        Ok(out)
    }

    /// Project full `6 * node_count` physical load vectors for each component onto the
    /// reduced space via $\mathbf{f}_r = \sum_c T_c^T \mathbf{f}_c$.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] if component counts or vector lengths mismatch or contain non-finite values.
    pub fn project_loads(&self, component_full_loads: &[&[f64]]) -> Result<Vec<f64>, PlateError> {
        if component_full_loads.len() != self.comp_t.len() {
            return Err(PlateError::BadConstraint {
                what: "component load slice count must match component count",
            });
        }
        let mut f_r = vec![0.0; self.model.free];
        for (c, &full_f) in component_full_loads.iter().enumerate() {
            let dof_map = &self.comp_dof_maps[c];
            if full_f.len() != dof_map.len() || full_f.iter().any(|x| !x.is_finite()) {
                return Err(PlateError::BadConstraint {
                    what: "component full load vector length mismatch or non-finite entry",
                });
            }
            let comp_rows = &self.comp_t[c];
            for (full_idx, maybe_free) in dof_map.iter().enumerate() {
                if let Some(f_idx) = *maybe_free {
                    let f_val = full_f[full_idx];
                    if f_val != 0.0 {
                        for &(r_idx, w) in &comp_rows[f_idx] {
                            f_r[r_idx] += w * f_val;
                        }
                    }
                }
            }
        }
        Ok(f_r)
    }

    /// Evaluate the master-side and slave-side virtual interface powers `(P_m, P_s)` for a
    /// reduced velocity field `reduced_velocity` (`len == self.model.free`) and an arbitrary
    /// discrete Lagrange multiplier / interface traction vector `lambda` (`len == self.receipt.slave_dofs`).
    ///
    /// Under the kinematic tie constraint $\mathbf{v}_s = T_{sm} \mathbf{v}_m$, work conjugacy
    /// $\mathbf{f}_m = T_{sm}^T \boldsymbol{\lambda}$, $\mathbf{f}_s = -\boldsymbol{\lambda}$
    /// guarantees $P_m + P_s = \mathbf{f}_m \cdot \mathbf{v}_m + \mathbf{f}_s \cdot \mathbf{v}_s = 0$.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] on length mismatch or non-finite inputs.
    pub fn interface_virtual_power(
        &self,
        reduced_velocity: &[f64],
        lambda: &[f64],
    ) -> Result<(f64, f64), PlateError> {
        if reduced_velocity.len() != self.model.free
            || reduced_velocity.iter().any(|x| !x.is_finite())
        {
            return Err(PlateError::BadConstraint {
                what: "reduced_velocity must match model.free and be finite",
            });
        }
        if lambda.len() != self.slave_couplings.len() || lambda.iter().any(|x| !x.is_finite()) {
            return Err(PlateError::BadConstraint {
                what: "lambda length must equal receipt.slave_dofs and be finite",
            });
        }
        let comp_v = self.expand_component_free(reduced_velocity)?;
        // Accumulate master reaction forces f_m on the reduced master DOFs: f_m = T_{sm}^T lambda
        let mut f_master_reduced = vec![0.0; self.model.free];
        let mut p_slave = 0.0;
        for (s_idx, coupling) in self.slave_couplings.iter().enumerate() {
            let lam = lambda[s_idx];
            let v_s = comp_v[coupling.slave_component][coupling.slave_component_free];
            p_slave += -lam * v_s;
            for &(r_idx, w) in &coupling.master_terms {
                f_master_reduced[r_idx] += w * lam;
            }
        }
        let mut p_master = 0.0;
        for (f_m, &v_m) in f_master_reduced.iter().zip(reduced_velocity.iter()) {
            p_master += f_m * v_m;
        }
        Ok((p_master, p_slave))
    }
}

/// Wohlmuth (2000) biorthogonal linear dual basis $(\psi_0(\xi), \psi_1(\xi)) = (2 - 3\xi,\ 3\xi - 1)$
/// on the reference segment $\xi \in [0, 1]$.
#[must_use]
pub const fn dual_basis(xi: f64) -> [f64; 2] {
    [2.0 - 3.0 * xi, 3.0 * xi - 1.0]
}

/// Standard P1 primal basis $(\phi_0(\xi), \phi_1(\xi)) = (1 - \xi,\ \xi)$ on $\xi \in [0, 1]$.
#[must_use]
pub const fn primal_basis(xi: f64) -> [f64; 2] {
    [1.0 - xi, xi]
}

/// Raw 1-D dual-mortar matrices along an interface between `master_points` and `slave_points`.
#[derive(Debug, Clone)]
pub struct PolylineMortarOperators {
    /// Diagonal slave mass entries $D_{ii} = \int_\Gamma \psi_i^s \phi_i^s\,ds$.
    pub d_diag: Vec<f64>,
    /// Full slave dual-primal mass matrix $D_{ij} = \int_\Gamma \psi_i^s \phi_j^s\,ds$ (diagonal to machine precision).
    pub d_full: Vec<Vec<f64>>,
    /// Raw mortar coupling matrix $M_{ij} = \int_\Gamma \psi_i^s \phi_j^m\,ds$.
    pub m: Vec<Vec<f64>>,
    /// Coordinate-consistent mortar transfer matrix $C_{ij}$ (`n_slave × n_master`),
    /// satisfying $\sum_j C_{ij} = 1$ and $\sum_j C_{ij} \mathbf{p}_{m,j} = \mathbf{p}_{s,i}$.
    pub c: Vec<Vec<f64>>,
    /// Number of mortar sub-segments integrated.
    pub segments: usize,
    /// Maximum interface gap [m].
    pub max_gap_m: f64,
    /// Maximum segment tangent angle mismatch [rad].
    pub max_angle_rad: f64,
}

fn norm3(v: [f64; 3]) -> f64 {
    det::sqrt(dot(v, v))
}

fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        (1.0 - t) * a[0] + t * b[0],
        (1.0 - t) * a[1] + t * b[1],
        (1.0 - t) * a[2] + t * b[2],
    ]
}

/// Cross-product matrix $S(\mathbf{d})$ such that $S(\mathbf{d})\boldsymbol{\theta} = \boldsymbol{\theta} \times \mathbf{d}$.
const fn skew_theta_cross_d(d: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, d[2], -d[1]], [-d[2], 0.0, d[0]], [d[1], -d[0], 0.0]]
}

/// Assemble the 1-D dual-mortar operators between a master polyline and a slave polyline.
///
/// Both `master_midsurface` and `master_offsets` have length `n_master >= 2`; `slave_points`
/// has length `n_slave >= 2`.
///
/// # Errors
/// Returns [`PlateError::BadConstraint`], [`PlateError::Budget`], or [`PlateError::Cancelled`]
/// on degenerate segments, excessive gap/angle mismatch, budget overflow, or cancellation.
pub fn assemble_polyline_dual_mortar(
    master_midsurface: &[[f64; 3]],
    master_offsets: &[[f64; 3]],
    slave_points: &[[f64; 3]],
    budget: &TieBudget,
    segments_used: &mut usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PolylineMortarOperators, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    let n_m = master_midsurface.len();
    let n_s = slave_points.len();
    if n_m < 2 || n_s < 2 || master_offsets.len() != n_m {
        return Err(PlateError::BadConstraint {
            what: "polylines require at least two nodes on master and slave sides",
        });
    }

    // Master segment lengths, cumulative arc-lengths, and unit tangents.
    let mut m_seg_len = Vec::with_capacity(n_m - 1);
    let mut m_seg_tan = Vec::with_capacity(n_m - 1);
    let mut m_arc = Vec::with_capacity(n_m);
    m_arc.push(0.0);
    for j in 0..(n_m - 1) {
        let chord = sub(master_midsurface[j + 1], master_midsurface[j]);
        let len = norm3(chord);
        if !len.is_finite() || len <= 1e-12 {
            return Err(PlateError::BadConstraint {
                what: "degenerate master tie polyline segment",
            });
        }
        m_seg_len.push(len);
        m_seg_tan.push(scale3(chord, 1.0 / len));
        m_arc.push(m_arc[j] + len);
    }
    let total_m_len = m_arc[n_m - 1];

    // Master offset skin attachment points p_m = x_m + d_m.
    let master_skin: Vec<[f64; 3]> = (0..n_m)
        .map(|j| add3(master_midsurface[j], master_offsets[j]))
        .collect();

    // Slave segment lengths and unit tangents.
    let mut s_seg_len = Vec::with_capacity(n_s - 1);
    let mut s_seg_tan = Vec::with_capacity(n_s - 1);
    let mut s_arc = Vec::with_capacity(n_s);
    s_arc.push(0.0);
    for k in 0..(n_s - 1) {
        let chord = sub(slave_points[k + 1], slave_points[k]);
        let len = norm3(chord);
        if !len.is_finite() || len <= 1e-12 {
            return Err(PlateError::BadConstraint {
                what: "degenerate slave tie polyline segment",
            });
        }
        s_seg_len.push(len);
        s_seg_tan.push(scale3(chord, 1.0 / len));
        s_arc.push(s_arc[k] + len);
    }

    // Project each slave node onto the master polyline to determine its master arc coordinate
    // and check Whether the slave nodes sit closer to master_skin or master_midsurface.
    let mut s_on_m_arc = vec![0.0; n_s];
    let mut max_gap_m = 0.0_f64;
    let mut sum_dist_skin = 0.0_f64;
    let mut sum_dist_mid = 0.0_f64;

    for (i, &sp) in slave_points.iter().enumerate() {
        let mut best_arc = 0.0;
        let mut best_dist_mid = f64::INFINITY;
        let mut best_dist_skin = f64::INFINITY;
        for j in 0..(n_m - 1) {
            // Project onto midsurface chord
            let a_mid = master_midsurface[j];
            let ab_mid = sub(master_midsurface[j + 1], a_mid);
            let xi_mid =
                (dot(sub(sp, a_mid), ab_mid) / (m_seg_len[j] * m_seg_len[j])).clamp(0.0, 1.0);
            let proj_mid = lerp3(a_mid, master_midsurface[j + 1], xi_mid);
            let d_mid = norm3(sub(sp, proj_mid));

            // Project onto skin chord
            let a_skin = master_skin[j];
            let ab_skin = sub(master_skin[j + 1], a_skin);
            let len_skin_sq = dot(ab_skin, ab_skin);
            let xi_skin = if len_skin_sq > 1e-24 {
                (dot(sub(sp, a_skin), ab_skin) / len_skin_sq).clamp(0.0, 1.0)
            } else {
                xi_mid
            };
            let proj_skin = lerp3(a_skin, master_skin[j + 1], xi_skin);
            let d_skin = norm3(sub(sp, proj_skin));

            let d_eff = d_mid.min(d_skin);
            let best_eff = best_dist_mid.min(best_dist_skin);
            if d_eff < best_eff {
                best_dist_mid = d_mid;
                best_dist_skin = d_skin;
                let xi = if d_skin < d_mid { xi_skin } else { xi_mid };
                best_arc = m_arc[j] + xi * m_seg_len[j];
            }
        }
        let gap_node = best_dist_mid.min(best_dist_skin);
        if gap_node > budget.gap_tolerance_m {
            return Err(PlateError::BadConstraint {
                what: "tie geometric gap exceeds budget.gap_tolerance_m",
            });
        }
        max_gap_m = max_gap_m.max(gap_node);
        sum_dist_skin += best_dist_skin;
        sum_dist_mid += best_dist_mid;
        s_on_m_arc[i] = best_arc;
    }

    // Ensure monotonic master arc-length parameterization along slave segments.
    // If projection of interior nodes on curved arcs slightly perturbs from proportional arc-length,
    // s_on_m_arc must remain strictly increasing.
    for k in 0..(n_s - 1) {
        if s_on_m_arc[k + 1] <= s_on_m_arc[k] + 1e-12 * total_m_len {
            // Fall back to normalized arc-length if projection is non-monotonic;
            // if the polylines are oppositely oriented, the tangent angle check below will catch and refuse it.
            let total_s_len = s_arc[n_s - 1];
            for i in 0..n_s {
                s_on_m_arc[i] = (s_arc[i] / total_s_len) * total_m_len;
            }
            break;
        }
    }

    let mut d_full = vec![vec![0.0; n_s]; n_s];
    let mut m_mat = vec![vec![0.0; n_m]; n_s];
    let mut max_angle_rad = 0.0_f64;
    let mut local_segments = 0usize;

    // 2-point Gauss-Legendre rule on [0, 1].
    let inv_two_sqrt3 = 0.5 / det::sqrt(3.0);
    let gauss_u = [0.5 - inv_two_sqrt3, 0.5 + inv_two_sqrt3];
    let gauss_w = [0.5, 0.5];

    for k in 0..(n_s - 1) {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let s_lo_arc = s_on_m_arc[k];
        let s_hi_arc = s_on_m_arc[k + 1];
        let ds_arc = s_hi_arc - s_lo_arc;
        if ds_arc <= 1e-14 * total_m_len {
            return Err(PlateError::BadConstraint {
                what: "degenerate slave-to-master arc projection on tie segment",
            });
        }

        for j in 0..(n_m - 1) {
            let m_lo_arc = m_arc[j];
            let m_hi_arc = m_arc[j + 1];
            let clip_lo = s_lo_arc.max(m_lo_arc);
            let clip_hi = s_hi_arc.min(m_hi_arc);
            if clip_hi - clip_lo <= 1e-12 * total_m_len {
                continue;
            }

            // Count mortar segment and enforce budget.
            *segments_used += 1;
            local_segments += 1;
            if *segments_used > budget.max_segments {
                return Err(PlateError::Budget {
                    what: "tie_segments",
                    requested: *segments_used,
                    limit: budget.max_segments,
                });
            }

            // Check tangent angle mismatch between slave segment k and master segment j.
            let cos_theta = dot(s_seg_tan[k], m_seg_tan[j]).clamp(-1.0, 1.0);
            let angle = det::acos(cos_theta);
            max_angle_rad = max_angle_rad.max(angle);
            if angle > budget.angle_tolerance_rad {
                return Err(PlateError::BadConstraint {
                    what: "tie tangent-frame angle mismatch exceeds budget.angle_tolerance_rad",
                });
            }

            let xi_s_lo = (clip_lo - s_lo_arc) / ds_arc;
            let xi_s_hi = (clip_hi - s_lo_arc) / ds_arc;
            let xi_m_lo = (clip_lo - m_lo_arc) / m_seg_len[j];
            let xi_m_hi = (clip_hi - m_lo_arc) / m_seg_len[j];
            let sub_len_s = (xi_s_hi - xi_s_lo) * s_seg_len[k];

            for q in 0..2 {
                let u = gauss_u[q];
                let xi_s = xi_s_lo + u * (xi_s_hi - xi_s_lo);
                let xi_m = xi_m_lo + u * (xi_m_hi - xi_m_lo);
                let dw = gauss_w[q] * sub_len_s;

                let x_s = lerp3(slave_points[k], slave_points[k + 1], xi_s);
                let x_m_mid = lerp3(master_midsurface[j], master_midsurface[j + 1], xi_m);
                let x_m_skin = lerp3(master_skin[j], master_skin[j + 1], xi_m);
                let gap = norm3(sub(x_s, x_m_mid)).min(norm3(sub(x_s, x_m_skin)));
                max_gap_m = max_gap_m.max(gap);
                if gap > budget.gap_tolerance_m {
                    return Err(PlateError::BadConstraint {
                        what: "tie quadrature gap exceeds budget.gap_tolerance_m",
                    });
                }

                let psi = dual_basis(xi_s);
                let phi_s = primal_basis(xi_s);
                let phi_m = primal_basis(xi_m);

                for a in 0..2 {
                    for b in 0..2 {
                        d_full[k + a][k + b] += dw * psi[a] * phi_s[b];
                    }
                    for c in 0..2 {
                        m_mat[k + a][j + c] += dw * psi[a] * phi_m[c];
                    }
                }
            }
        }
    }

    let mut d_diag = vec![0.0; n_s];
    let mut c_mat = vec![vec![0.0; n_m]; n_s];
    for i in 0..n_s {
        let d_ii = d_full[i][i];
        if !d_ii.is_finite() || d_ii <= 1e-15 {
            return Err(PlateError::BadConstraint {
                what: "zero or non-positive slave dual-mortar diagonal mass",
            });
        }
        d_diag[i] = d_ii;
        for j in 0..n_m {
            c_mat[i][j] = m_mat[i][j] / d_ii;
        }
    }

    // Apply Puso (2004) / Petrov-Galerkin coordinate consistency correction on curved interfaces
    // so that sum_j C_{ij} = 1 and sum_j C_{ij} p_{m,j} = p_{s,i} hold to machine precision.
    let use_skin_coords = sum_dist_skin + 1e-14 < sum_dist_mid;
    let master_ref: &[[f64; 3]] = if use_skin_coords {
        &master_skin
    } else {
        master_midsurface
    };
    enforce_linear_coordinate_reproduction(&mut c_mat, master_ref, slave_points);

    Ok(PolylineMortarOperators {
        d_diag,
        d_full,
        m: m_mat,
        c: c_mat,
        segments: local_segments,
        max_gap_m,
        max_angle_rad,
    })
}

/// Orthonormal basis of the affine hull of `points` around ` points[0]`.
fn affine_hull_basis(points: &[[f64; 3]]) -> ([f64; 3], Vec<[f64; 3]>) {
    let origin = points[0];
    let mut max_dist = 0.0_f64;
    for p in points.iter().skip(1) {
        max_dist = max_dist.max(norm3(sub(*p, origin)));
    }
    let tol = 1e-10 * max_dist.max(1.0);
    let mut basis: Vec<[f64; 3]> = Vec::with_capacity(3);
    for p in points.iter().skip(1) {
        let mut v = sub(*p, origin);
        // Two-pass Modified Gram-Schmidt for numerical stability.
        for _ in 0..2 {
            for b in &basis {
                let proj = dot(v, *b);
                v = sub(v, scale3(*b, proj));
            }
        }
        let n = norm3(v);
        if n > tol {
            basis.push(scale3(v, 1.0 / n));
            if basis.len() == 3 {
                break;
            }
        }
    }
    (origin, basis)
}

/// Solve a small dense symmetric positive-definite linear system $A x = b$ ($m \le 4$).
fn solve_small_spd(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let m = b.len();
    for k in 0..m {
        let mut piv = k;
        let mut max_v = a[k][k].abs();
        for r in (k + 1)..m {
            if a[r][k].abs() > max_v {
                max_v = a[r][k].abs();
                piv = r;
            }
        }
        if max_v <= 1e-22 {
            return None;
        }
        if piv != k {
            a.swap(k, piv);
            b.swap(k, piv);
        }
        let diag = a[k][k];
        for r in (k + 1)..m {
            let factor = a[r][k] / diag;
            for c in k..m {
                a[r][c] -= factor * a[k][c];
            }
            b[r] -= factor * b[k];
        }
    }
    let mut x = vec![0.0; m];
    for k in (0..m).rev() {
        let mut sum = b[k];
        for c in (k + 1)..m {
            sum -= a[k][c] * x[c];
        }
        x[k] = sum / a[k][k];
    }
    Some(x)
}

/// Minimum-norm local adjustment of each row $c_i$ of $C$ so that
/// $\sum_j C_{ij} = 1$ and $\sum_j C_{ij} \mathbf{p}_{m,j} = \mathbf{p}_{s,i}$ (in the affine span of the master polyline).
/// On straight seams or matching meshes, $C = D^{-1}M$ already satisfies these identities to machine precision,
/// so the adjustment is zero.
fn enforce_linear_coordinate_reproduction(
    c_mat: &mut [Vec<f64>],
    master_pts: &[[f64; 3]],
    slave_pts: &[[f64; 3]],
) {
    let n_m = master_pts.len();
    let (origin, basis) = affine_hull_basis(master_pts);
    let dim = basis.len();
    if dim == 0 {
        return;
    }
    let needed_stencil = (dim + 1).min(n_m);

    for (i, row) in c_mat.iter_mut().enumerate() {
        // Active stencil of master nodes for slave node i.
        let mut stencil: BTreeSet<usize> = (0..n_m).filter(|&j| row[j].abs() > 1e-14).collect();
        if stencil.is_empty() {
            // Find closest master node.
            let mut best_j = 0;
            let mut best_d = f64::INFINITY;
            for (j, mp) in master_pts.iter().enumerate() {
                let d = norm3(sub(slave_pts[i], *mp));
                if d < best_d {
                    best_d = d;
                    best_j = j;
                }
            }
            stencil.insert(best_j);
        }
        while stencil.len() < needed_stencil {
            let first = *stencil.iter().next().unwrap_or(&0);
            let last = *stencil.iter().next_back().unwrap_or(&0);
            let add_left = if first > 0 {
                Some((first - 1, norm3(sub(slave_pts[i], master_pts[first - 1]))))
            } else {
                None
            };
            let add_right = if last + 1 < n_m {
                Some((last + 1, norm3(sub(slave_pts[i], master_pts[last + 1]))))
            } else {
                None
            };
            match (add_left, add_right) {
                (Some((jl, dl)), Some((jr, dr))) => {
                    if dl <= dr {
                        stencil.insert(jl);
                    } else {
                        stencil.insert(jr);
                    }
                }
                (Some((jl, _)), None) => {
                    stencil.insert(jl);
                }
                (None, Some((jr, _))) => {
                    stencil.insert(jr);
                }
                (None, None) => break,
            }
        }

        let j_list: Vec<usize> = stencil.into_iter().collect();
        let m_rows = 1 + dim;
        let n_cols = j_list.len();
        // Build constraint matrix A (m_rows x n_cols) and target vector b (m_rows).
        let mut a_mat = vec![vec![0.0; n_cols]; m_rows];
        let mut b_vec = vec![0.0; m_rows];
        b_vec[0] = 1.0;
        for (col, _) in j_list.iter().enumerate() {
            a_mat[0][col] = 1.0;
        }
        let rel_s = sub(slave_pts[i], origin);
        for d in 0..dim {
            b_vec[1 + d] = dot(rel_s, basis[d]);
            for (col, &j) in j_list.iter().enumerate() {
                let rel_m = sub(master_pts[j], origin);
                a_mat[1 + d][col] = dot(rel_m, basis[d]);
            }
        }

        // Compute residual r = b - A c_0
        let mut res = b_vec;
        for r in 0..m_rows {
            for (col, &j) in j_list.iter().enumerate() {
                res[r] -= a_mat[r][col] * row[j];
            }
        }
        let max_res = res.iter().fold(0.0_f64, |acc, x| acc.max(x.abs()));
        if max_res <= 1e-15 {
            continue;
        }

        // Solve (A A^T) mu = res, then c = c_0 + A^T mu.
        let mut aat = vec![vec![0.0; m_rows]; m_rows];
        for r1 in 0..m_rows {
            for r2 in 0..m_rows {
                let mut sum = 0.0;
                for col in 0..n_cols {
                    sum += a_mat[r1][col] * a_mat[r2][col];
                }
                aat[r1][r2] = sum;
            }
        }
        if let Some(mu) = solve_small_spd(aat, res) {
            for (col, &j) in j_list.iter().enumerate() {
                let mut delta = 0.0;
                for r in 0..m_rows {
                    delta += a_mat[r][col] * mu[r];
                }
                row[j] += delta;
            }
        }
    }
}

/// Area-weighted nodal unit normals for a component mesh.
fn compute_component_normals(mesh: &ShellMesh) -> Result<Vec<[f64; 3]>, PlateError> {
    let mut normals = vec![[0.0; 3]; mesh.nodes.len()];
    for (e, tri) in mesh.tris.iter().enumerate() {
        let g = mesh.facet(e)?;
        let weighted_n = scale3(g.frame[2], g.area_m2);
        for &node in tri {
            normals[node] = add3(normals[node], weighted_n);
        }
    }
    for (node, n) in normals.iter_mut().enumerate() {
        let len = norm3(*n);
        if !len.is_finite() || len <= 1e-12 {
            return Err(PlateError::BadDirector { node });
        }
        *n = scale3(*n, 1.0 / len);
    }
    Ok(normals)
}

/// Tie multiple [`ShellComponent`] models along 1-D [`EdgeTie`] seams using work-conjugate
/// dual-mortar master-slave elimination ($K_r = T^T K T,\ M_r = T^T M T$).
///
/// # Errors
/// Returns [`PlateError::BadConstraint`] on invalid component/node references, cross-point
/// double constraints, `TieKind::Compliant`, or gap/angle tolerance violations;
/// [`PlateError::Budget`] when mortar segments or slave DOFs exceed `budget`; and
/// [`PlateError::Cancelled`] on cooperative cancellation.
pub fn tie_shells(
    components: &[ShellComponent<'_>],
    ties: &[EdgeTie],
    budget: TieBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<TiedShellModel, PlateError> {
    tie_shells_with_faces(components, ties, &[], budget, cancelled)
}

/// Tie multiple [`ShellComponent`] models along both 1-D [`EdgeTie`] polylines and 2-D [`FaceTie`]
/// overlapping patches.
///
/// # Errors
/// Same refusal conditions as [`tie_shells`].
#[allow(clippy::too_many_lines)]
pub fn tie_shells_with_faces(
    components: &[ShellComponent<'_>],
    edge_ties: &[EdgeTie],
    face_ties: &[FaceTie],
    budget: TieBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<TiedShellModel, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    if !budget.gap_tolerance_m.is_finite()
        || budget.gap_tolerance_m < 0.0
        || !budget.angle_tolerance_rad.is_finite()
        || budget.angle_tolerance_rad < 0.0
    {
        return Err(PlateError::BadConstraint {
            what: "tie budget tolerances must be finite and non-negative",
        });
    }
    if components.is_empty() {
        return Err(PlateError::BadConstraint {
            what: "at least one shell component is required",
        });
    }

    // Validate components and compute nodal normals.
    let mut comp_normals = Vec::with_capacity(components.len());
    for (idx, comp) in components.iter().enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if comp.id != idx {
            return Err(PlateError::BadConstraint {
                what: "component id must match its canonical index in components slice",
            });
        }
        comp.mesh.validate()?;
        let n_nodes = comp.mesh.nodes.len();
        if comp.nodal_thickness_m.len() != n_nodes
            || comp
                .nodal_thickness_m
                .iter()
                .any(|&h| !h.is_finite() || h <= 0.0)
        {
            return Err(PlateError::BadSection {
                what: "component nodal_thickness_m must match node count and be positive finite",
            });
        }
        if comp.model.dof_map.len() != 6 * n_nodes
            || comp.model.k.nrows() != comp.model.free
            || comp.model.k.ncols() != comp.model.free
            || comp.model.m.nrows() != comp.model.free
            || comp.model.m.ncols() != comp.model.free
        {
            return Err(PlateError::BadConstraint {
                what: "component ShellModel dimensions or dof_map length mismatch",
            });
        }
        comp_normals.push(compute_component_normals(comp.mesh)?);
    }

    let mut slave_nodes_seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut master_nodes_seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    let mut slave_edge_spans_seen: BTreeSet<(usize, usize, usize)> = BTreeSet::new();
    let mut corner_parent: BTreeMap<(usize, usize), (usize, usize)> = BTreeMap::new();

    // Count how many EdgeTies touch each (component, endpoint_node) so multi-tie 3-way corner
    // crosspoints use Wohlmuth (2000) §2.3 wirebasket corner pinning and consistent corner root resolution.
    let mut endpoint_degree: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for tie in edge_ties {
        if let (Some(&m_first), Some(&m_last)) = (tie.a.nodes.first(), tie.a.nodes.last()) {
            *endpoint_degree
                .entry((tie.a.component, m_first))
                .or_insert(0) += 1;
            if m_last != m_first {
                *endpoint_degree
                    .entry((tie.a.component, m_last))
                    .or_insert(0) += 1;
            }
        }
        if let (Some(&s_first), Some(&s_last)) = (tie.b.nodes.first(), tie.b.nodes.last()) {
            *endpoint_degree
                .entry((tie.b.component, s_first))
                .or_insert(0) += 1;
            if s_last != s_first {
                *endpoint_degree
                    .entry((tie.b.component, s_last))
                    .or_insert(0) += 1;
            }
        }
    }

    // Map from (slave_comp, slave_node, slave_dof_0..6) -> Vec<(master_comp, master_node, master_dof_0..6, weight)>
    let mut nodal_slave_constraints: BTreeMap<
        (usize, usize, usize),
        Vec<(usize, usize, usize, f64)>,
    > = BTreeMap::new();

    // Pending compliant glue nodal patches to assemble after reduced DOFs and T_c are constructed.
    let mut compliant_patches: Vec<CompliantNodalPatch> = Vec::new();
    let mut compliant_dofs_count = 0usize;

    let mut segments_used = 0usize;
    let mut max_gap_m = 0.0_f64;
    let mut max_angle_rad = 0.0_f64;

    // Process 1-D EdgeTies.
    for tie in edge_ties {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if let TieKind::Compliant(glue) = &tie.kind {
            glue.validate()?;
        }
        let side_sign = match tie.offset {
            TieOffset::Midsurface => 0.0,
            TieOffset::Skin { side } => {
                if side != -1 && side != 1 {
                    return Err(PlateError::BadConstraint {
                        what: "TieOffset::Skin side must be -1 or +1",
                    });
                }
                f64::from(side)
            }
        };

        let c_m = tie.a.component;
        let c_s = tie.b.component;
        if c_m >= components.len() || c_s >= components.len() || c_m == c_s {
            return Err(PlateError::BadConstraint {
                what: "EdgeTie component indices must be in range and distinct",
            });
        }
        let n_m = tie.a.nodes.len();
        let n_s = tie.b.nodes.len();
        if n_m < 2 || n_s < 2 {
            return Err(PlateError::BadConstraint {
                what: "EdgeTie polylines must contain at least 2 nodes",
            });
        }

        // Check in-range and duplicate nodes within each side.
        let mut local_m_seen = BTreeSet::new();
        for &mn in &tie.a.nodes {
            if mn >= components[c_m].mesh.nodes.len() || !local_m_seen.insert(mn) {
                return Err(PlateError::BadConstraint {
                    what: "EdgeTie master node out of range or duplicated within tie",
                });
            }
        }
        let mut local_s_seen = BTreeSet::new();
        for &sn in &tie.b.nodes {
            if sn >= components[c_s].mesh.nodes.len() || !local_s_seen.insert(sn) {
                return Err(PlateError::BadConstraint {
                    what: "EdgeTie slave node out of range or duplicated within tie",
                });
            }
        }

        let s_first = tie.b.nodes[0];
        let s_last = tie.b.nodes[n_s - 1];
        let s_span = (c_s, s_first.min(s_last), s_first.max(s_last));
        if !slave_edge_spans_seen.insert(s_span) {
            return Err(PlateError::BadConstraint {
                what: "cross-point double constraint: duplicate slave edge span on component",
            });
        }

        let mut skip_slave_node = vec![false; n_s];
        let corner_0_shared = endpoint_degree
            .get(&(c_m, tie.a.nodes[0]))
            .copied()
            .unwrap_or(0)
            > 1
            || endpoint_degree
                .get(&(c_s, tie.b.nodes[0]))
                .copied()
                .unwrap_or(0)
                > 1;
        let corner_last_shared = endpoint_degree
            .get(&(c_m, tie.a.nodes[n_m - 1]))
            .copied()
            .unwrap_or(0)
            > 1
            || endpoint_degree
                .get(&(c_s, tie.b.nodes[n_s - 1]))
                .copied()
                .unwrap_or(0)
                > 1;

        if matches!(tie.kind, TieKind::Rigid) {
            for (mi, &mn) in tie.a.nodes.iter().enumerate() {
                let is_end = mi == 0 || mi == n_m - 1;
                let shared_corner =
                    (mi == 0 && corner_0_shared) || (mi == n_m - 1 && corner_last_shared);
                if (!is_end || !shared_corner) && slave_nodes_seen.contains(&(c_m, mn)) {
                    return Err(PlateError::BadConstraint {
                        what: "cross-point double constraint: master node is already a slave node",
                    });
                }
            }
            for (si, &sn) in tie.b.nodes.iter().enumerate() {
                let is_end = si == 0 || si == n_s - 1;
                let shared_corner =
                    (si == 0 && corner_0_shared) || (si == n_s - 1 && corner_last_shared);
                if !is_end || !shared_corner {
                    if slave_nodes_seen.contains(&(c_s, sn))
                        || master_nodes_seen.contains(&(c_s, sn))
                    {
                        return Err(PlateError::BadConstraint {
                            what: "cross-point double constraint: slave node cannot be constrained twice or act as master",
                        });
                    }
                    if is_end {
                        let mn_end = if si == 0 {
                            tie.a.nodes[0]
                        } else {
                            tie.a.nodes[n_m - 1]
                        };
                        corner_parent.insert((c_s, sn), (c_m, mn_end));
                    }
                } else {
                    // Multi-tie 3-way corner crosspoint: resolve corner roots of (c_s, sn) and (c_m, mn_end).
                    let mn_end = if si == 0 {
                        tie.a.nodes[0]
                    } else {
                        tie.a.nodes[n_m - 1]
                    };
                    let root_s = find_corner_root(&corner_parent, (c_s, sn), components.len())?;
                    let root_m = find_corner_root(&corner_parent, (c_m, mn_end), components.len())?;
                    if root_s == root_m {
                        // This 3-way corner triangle is already closed to the same root master corner.
                        skip_slave_node[si] = true;
                    } else if root_s == (c_s, sn) && !master_nodes_seen.contains(&(c_s, sn)) {
                        corner_parent.insert((c_s, sn), (c_m, mn_end));
                    } else {
                        return Err(PlateError::BadConstraint {
                            what: "cross-point double constraint: conflicting master roots at corner",
                        });
                    }
                }
            }
            for (mi, &mn) in tie.a.nodes.iter().enumerate() {
                let shared_corner =
                    (mi == 0 && corner_0_shared) || (mi == n_m - 1 && corner_last_shared);
                if !shared_corner {
                    master_nodes_seen.insert((c_m, mn));
                }
            }
            for (si, &sn) in tie.b.nodes.iter().enumerate() {
                if !skip_slave_node[si] {
                    slave_nodes_seen.insert((c_s, sn));
                }
            }
        }

        let master_midsurface: Vec<[f64; 3]> = tie
            .a
            .nodes
            .iter()
            .map(|&mn| components[c_m].mesh.nodes[mn])
            .collect();
        let master_offsets: Vec<[f64; 3]> = tie
            .a
            .nodes
            .iter()
            .map(|&mn| {
                if side_sign == 0.0 {
                    [0.0, 0.0, 0.0]
                } else {
                    let h = components[c_m].nodal_thickness_m[mn];
                    scale3(comp_normals[c_m][mn], side_sign * 0.5 * h)
                }
            })
            .collect();
        let slave_points: Vec<[f64; 3]> = tie
            .b
            .nodes
            .iter()
            .map(|&sn| components[c_s].mesh.nodes[sn])
            .collect();

        let mut ops = assemble_polyline_dual_mortar(
            &master_midsurface,
            &master_offsets,
            &slave_points,
            &budget,
            &mut segments_used,
            &mut cancelled,
        )?;
        max_gap_m = max_gap_m.max(ops.max_gap_m);
        max_angle_rad = max_angle_rad.max(ops.max_angle_rad);

        // Wohlmuth (2000) §2.3 wirebasket crosspoint condition: at multi-tie corner crosspoints
        // where master and slave endpoints coincide, pin the corner row of C to the master corner.
        if corner_0_shared {
            for val in &mut ops.c[0] {
                *val = 0.0;
            }
            ops.c[0][0] = 1.0;
        }
        if corner_last_shared {
            for val in &mut ops.c[n_s - 1] {
                *val = 0.0;
            }
            ops.c[n_s - 1][n_m - 1] = 1.0;
        }

        // After the crosspoint pinning, so a pinned corner row takes the master
        // corner's tangent.
        let tangents = slave_seam_tangents(&master_midsurface, &ops.c)?;
        match tie.kind {
            TieKind::Rigid => {
                record_slave_nodal_couplings_filtered(
                    c_m,
                    &tie.a.nodes,
                    &master_offsets,
                    c_s,
                    &tie.b.nodes,
                    &ops.c,
                    &skip_slave_node,
                    tie.rotations,
                    Some(&tangents),
                    &mut nodal_slave_constraints,
                );
            }
            TieKind::Compliant(glue) => {
                if tie.rotations == RotationCoupling::SeamTangent {
                    return Err(PlateError::BadConstraint {
                        what: "a compliant glue tie couples Full or TranslationsOnly rotations; SeamTangent is rigid-tie kinematics",
                    });
                }
                build_compliant_polyline_patches(
                    components,
                    &comp_normals,
                    c_m,
                    &tie.a.nodes,
                    &master_offsets,
                    c_s,
                    &tie.b.nodes,
                    &slave_points,
                    &ops,
                    tie.rotations,
                    &glue,
                    &mut compliant_patches,
                    &mut compliant_dofs_count,
                )?;
            }
        }
    }

    // Process 2-D FaceTies (if any).
    for ftie in face_ties {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if let TieKind::Compliant(glue) = &ftie.kind {
            glue.validate()?;
        }
        let side_sign = match ftie.offset {
            TieOffset::Midsurface => 0.0,
            TieOffset::Skin { side } => {
                if side != -1 && side != 1 {
                    return Err(PlateError::BadConstraint {
                        what: "TieOffset::Skin side must be -1 or +1",
                    });
                }
                f64::from(side)
            }
        };
        if ftie.rotations == RotationCoupling::SeamTangent {
            return Err(PlateError::BadConstraint {
                what: "RotationCoupling::SeamTangent needs a 1-D seam; a FaceTie has none",
            });
        }
        let c_m = ftie.master_component;
        let c_s = ftie.slave_component;
        if c_m >= components.len()
            || c_s >= components.len()
            || c_m == c_s
            || ftie.master_tris.is_empty()
            || ftie.slave_tris.is_empty()
        {
            return Err(PlateError::BadConstraint {
                what: "FaceTie components and triangle lists must be valid and non-empty",
            });
        }
        let m_mesh = components[c_m].mesh;
        let s_mesh = components[c_s].mesh;
        let mut m_nodes_set = BTreeSet::new();
        for &t in &ftie.master_tris {
            let tri = *m_mesh.tris.get(t).ok_or(PlateError::BadConstraint {
                what: "FaceTie master triangle index out of range",
            })?;
            for n in tri {
                if matches!(ftie.kind, TieKind::Rigid) && slave_nodes_seen.contains(&(c_m, n)) {
                    return Err(PlateError::BadConstraint {
                        what: "cross-point double constraint in FaceTie master node",
                    });
                }
                m_nodes_set.insert(n);
            }
        }
        let mut s_nodes_set = BTreeSet::new();
        let mut s_nodal_area: BTreeMap<usize, f64> = BTreeMap::new();
        for &t in &ftie.slave_tris {
            let tri = *s_mesh.tris.get(t).ok_or(PlateError::BadConstraint {
                what: "FaceTie slave triangle index out of range",
            })?;
            segments_used += 1;
            if segments_used > budget.max_segments {
                return Err(PlateError::Budget {
                    what: "tie_segments",
                    requested: segments_used,
                    limit: budget.max_segments,
                });
            }
            let gs = s_mesh.facet(t)?;
            let mut min_angle = f64::INFINITY;
            for &mt in &ftie.master_tris {
                let gm = m_mesh.facet(mt)?;
                let ang = det::acos(dot(gs.frame[2], gm.frame[2]).clamp(-1.0, 1.0));
                min_angle = min_angle.min(ang);
            }
            max_angle_rad = max_angle_rad.max(min_angle);
            if min_angle > budget.angle_tolerance_rad {
                return Err(PlateError::BadConstraint {
                    what: "FaceTie normal angle mismatch exceeds budget.angle_tolerance_rad",
                });
            }
            for n in tri {
                s_nodes_set.insert(n);
                *s_nodal_area.entry(n).or_insert(0.0) += gs.area_m2 / 3.0;
            }
        }
        if matches!(ftie.kind, TieKind::Rigid) {
            for &n in &s_nodes_set {
                if slave_nodes_seen.contains(&(c_s, n)) || master_nodes_seen.contains(&(c_s, n)) {
                    return Err(PlateError::BadConstraint {
                        what: "cross-point double constraint in FaceTie slave node",
                    });
                }
            }
            for &n in &m_nodes_set {
                master_nodes_seen.insert((c_m, n));
            }
            for &n in &s_nodes_set {
                slave_nodes_seen.insert((c_s, n));
            }
        }

        let m_nodes: Vec<usize> = m_nodes_set.into_iter().collect();
        let m_offsets: Vec<[f64; 3]> = m_nodes
            .iter()
            .map(|&mn| {
                if side_sign == 0.0 {
                    [0.0, 0.0, 0.0]
                } else {
                    scale3(
                        comp_normals[c_m][mn],
                        side_sign * 0.5 * components[c_m].nodal_thickness_m[mn],
                    )
                }
            })
            .collect();
        let s_nodes: Vec<usize> = s_nodes_set.into_iter().collect();
        let mut c_mat = vec![vec![0.0; m_nodes.len()]; s_nodes.len()];
        for (si, &sn) in s_nodes.iter().enumerate() {
            let sp = s_mesh.nodes[sn];
            let mut best_gap = f64::INFINITY;
            let mut best_tri = ftie.master_tris[0];
            let mut best_bary = [1.0, 0.0, 0.0];
            for &mt in &ftie.master_tris {
                let tri = m_mesh.tris[mt];
                let p0 = m_mesh.nodes[tri[0]];
                let p1 = m_mesh.nodes[tri[1]];
                let p2 = m_mesh.nodes[tri[2]];
                let (bary, proj) = closest_barycentric_triangle(sp, p0, p1, p2);
                let d_mid = norm3(sub(sp, proj));
                let off_proj = add3(
                    proj,
                    scale3(
                        comp_normals[c_m][tri[0]],
                        side_sign * 0.5 * components[c_m].nodal_thickness_m[tri[0]],
                    ),
                );
                let d_eff = d_mid.min(norm3(sub(sp, off_proj)));
                if d_eff < best_gap {
                    best_gap = d_eff;
                    best_tri = mt;
                    best_bary = bary;
                }
            }
            max_gap_m = max_gap_m.max(best_gap);
            if best_gap > budget.gap_tolerance_m {
                return Err(PlateError::BadConstraint {
                    what: "FaceTie gap exceeds budget.gap_tolerance_m",
                });
            }
            let tri = m_mesh.tris[best_tri];
            for (vertex_slot, &mn) in tri.iter().enumerate() {
                if let Ok(col) = m_nodes.binary_search(&mn) {
                    c_mat[si][col] = best_bary[vertex_slot];
                }
            }
        }
        match ftie.kind {
            TieKind::Rigid => {
                let no_skip = vec![false; s_nodes.len()];
                record_slave_nodal_couplings_filtered(
                    c_m,
                    &m_nodes,
                    &m_offsets,
                    c_s,
                    &s_nodes,
                    &c_mat,
                    &no_skip,
                    ftie.rotations,
                    None,
                    &mut nodal_slave_constraints,
                );
            }
            TieKind::Compliant(glue) => {
                let g0 = m_mesh.facet(ftie.master_tris[0])?;
                for (si, &sn) in s_nodes.iter().enumerate() {
                    let area_i = s_nodal_area.get(&sn).copied().unwrap_or(0.0);
                    append_compliant_nodal_patch(
                        components,
                        c_m,
                        &m_nodes,
                        &m_offsets,
                        c_s,
                        sn,
                        &c_mat[si],
                        g0.frame[2],
                        g0.frame[0],
                        g0.frame[1],
                        area_i,
                        ftie.rotations,
                        &glue,
                        &mut compliant_patches,
                        &mut compliant_dofs_count,
                    );
                }
            }
        }
    }

    // Transitively resolve slave -> master references at multi-tie corners so every slave DOF
    // is expressed strictly in terms of unconstrained root master DOFs.
    resolve_transitive_slave_constraints(&mut nodal_slave_constraints, components.len())?;

    // Count active slave DOFs (those that were free in their component model) and check budget.
    let mut slave_dofs_count = 0usize;
    for &(c_s, sn, d_s) in nodal_slave_constraints.keys() {
        if components[c_s].model.dof_map[6 * sn + d_s].is_some() {
            slave_dofs_count += 1;
        }
    }
    let total_budget_dofs = slave_dofs_count + compliant_dofs_count;
    if total_budget_dofs > budget.max_slave_dofs {
        return Err(PlateError::Budget {
            what: "tie_slave_dofs",
            requested: total_budget_dofs,
            limit: budget.max_slave_dofs,
        });
    }

    // Assign reduced free DOF indices to all retained (non-slave, non-fixed) DOFs across components.
    let mut map: Vec<Vec<[Option<usize>; 6]>> = Vec::with_capacity(components.len());
    let mut reduced_free = 0usize;
    for (c, comp) in components.iter().enumerate() {
        let n_nodes = comp.mesh.nodes.len();
        let mut comp_map = vec![[None; 6]; n_nodes];
        for n in 0..n_nodes {
            for d in 0..6 {
                if comp.model.dof_map[6 * n + d].is_some()
                    && !nodal_slave_constraints.contains_key(&(c, n, d))
                {
                    comp_map[n][d] = Some(reduced_free);
                    reduced_free += 1;
                }
            }
        }
        map.push(comp_map);
    }
    if reduced_free == 0 {
        return Err(PlateError::BadConstraint {
            what: "tied shell assembly has zero retained free coordinates",
        });
    }

    // Build per-component transformation matrices T_c (free_c x reduced_free) and SlaveDofCoupling list.
    let mut comp_t: Vec<Vec<Vec<(usize, f64)>>> = Vec::with_capacity(components.len());
    let mut comp_dof_maps: Vec<Vec<Option<usize>>> = Vec::with_capacity(components.len());
    let mut slave_couplings: Vec<SlaveDofCoupling> = Vec::with_capacity(slave_dofs_count);

    for (c, comp) in components.iter().enumerate() {
        comp_dof_maps.push(comp.model.dof_map.clone());
        let mut rows = vec![Vec::new(); comp.model.free];
        let n_nodes = comp.mesh.nodes.len();
        for n in 0..n_nodes {
            for d in 0..6 {
                let Some(local_free) = comp.model.dof_map[6 * n + d] else {
                    continue;
                };
                if let Some(r_idx) = map[c][n][d] {
                    rows[local_free].push((r_idx, 1.0));
                } else if let Some(master_list) = nodal_slave_constraints.get(&(c, n, d)) {
                    let mut accum: BTreeMap<usize, f64> = BTreeMap::new();
                    for &(c_m, mn, d_m, w) in master_list {
                        if w.abs() <= 1e-16 {
                            continue;
                        }
                        if let Some(r_m) = map[c_m][mn][d_m] {
                            *accum.entry(r_m).or_insert(0.0) += w;
                        }
                    }
                    let terms: Vec<(usize, f64)> = accum
                        .into_iter()
                        .filter(|&(_, w)| w.abs() > 1e-16)
                        .collect();
                    rows[local_free].clone_from(&terms);
                    slave_couplings.push(SlaveDofCoupling {
                        slave_component: c,
                        slave_node: n,
                        slave_dof: d,
                        slave_component_free: local_free,
                        master_terms: terms,
                    });
                }
            }
        }
        comp_t.push(rows);
    }

    // Assemble K_r = sum_c T_c^T K_c T_c + K_glue and M_r = sum_c T_c^T M_c T_c with exact symmetry.
    let mut k_coo = Coo::new(reduced_free, reduced_free);
    let mut m_coo = Coo::new(reduced_free, reduced_free);
    for (c, comp) in components.iter().enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        accumulate_congruence(&comp.model.k, &comp_t[c], &mut k_coo, &mut cancelled)?;
        accumulate_congruence(&comp.model.m, &comp_t[c], &mut m_coo, &mut cancelled)?;
    }

    let mut k_glue_coo = Coo::new(reduced_free, reduced_free);
    let mut c_glue_coo = Coo::new(reduced_free, reduced_free);
    let mut has_glue_damping = false;
    for patch in &compliant_patches {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        accumulate_compliant_patch(
            patch,
            &comp_dof_maps,
            &comp_t,
            &mut k_coo,
            &mut k_glue_coo,
            &mut c_glue_coo,
            &mut has_glue_damping,
        );
    }

    let k_r = k_coo.assemble();
    let m_r = m_coo.assemble();
    let glue_stiffness = k_glue_coo.assemble();
    let glue_damping = if has_glue_damping {
        Some(c_glue_coo.assemble())
    } else {
        None
    };

    // Concatenate global dof_map across components.
    let total_full_dofs: usize = components.iter().map(|c| 6 * c.mesh.nodes.len()).sum();
    let mut combined_dof_map = Vec::with_capacity(total_full_dofs);
    for comp_map in &map {
        for node_dofs in comp_map {
            combined_dof_map.extend_from_slice(node_dofs);
        }
    }

    Ok(TiedShellModel {
        model: ShellModel {
            k: k_r,
            m: m_r,
            dof_map: combined_dof_map,
            free: reduced_free,
        },
        map,
        receipt: TieReceipt {
            max_gap_m,
            max_angle_rad,
            segments: segments_used,
            slave_dofs: total_budget_dofs,
        },
        glue_stiffness,
        glue_damping,
        comp_t,
        comp_dof_maps,
        slave_couplings,
    })
}

/// Unit seam tangent at each slave node: the master polyline's unit nodal
/// tangents (central differences, one-sided at the ends) carried by the mortar
/// transfer `c`, so a straight seam gets its exact direction.
fn slave_seam_tangents(master: &[[f64; 3]], c: &[Vec<f64>]) -> Result<Vec<[f64; 3]>, PlateError> {
    let degenerate = || PlateError::BadConstraint {
        what: "EdgeTie seam tangent is degenerate",
    };
    let last = master.len() - 1;
    let mut nodal = Vec::with_capacity(master.len());
    for j in 0..=last {
        let d = sub(master[(j + 1).min(last)], master[j.saturating_sub(1)]);
        let len = norm3(d);
        if !(len.is_finite() && len > 0.0) {
            return Err(degenerate());
        }
        nodal.push(scale3(d, 1.0 / len));
    }
    c.iter()
        .map(|row| {
            let t = row
                .iter()
                .zip(&nodal)
                .fold([0.0; 3], |t, (&w, &tj)| add3(t, scale3(tj, w)));
            // A seam folding back within one mortar support cancels its
            // tangents; refuse rather than pick a direction.
            let len = norm3(t);
            if len.is_finite() && len > 0.5 {
                Ok(scale3(t, 1.0 / len))
            } else {
                Err(degenerate())
            }
        })
        .collect()
}

fn find_corner_root(
    corner_parent: &BTreeMap<(usize, usize), (usize, usize)>,
    start: (usize, usize),
    max_hops: usize,
) -> Result<(usize, usize), PlateError> {
    let mut cur = start;
    for _ in 0..=max_hops {
        if let Some(&next) = corner_parent.get(&cur) {
            if next == start {
                return Err(PlateError::BadConstraint {
                    what: "cyclic corner constraint chain in EdgeTies",
                });
            }
            cur = next;
        } else {
            return Ok(cur);
        }
    }
    Err(PlateError::BadConstraint {
        what: "cyclic corner constraint chain in EdgeTies",
    })
}

fn resolve_transitive_slave_constraints(
    constraints: &mut BTreeMap<(usize, usize, usize), Vec<(usize, usize, usize, f64)>>,
    max_hops: usize,
) -> Result<(), PlateError> {
    for _ in 0..=max_hops {
        let mut any_expanded = false;
        let snapshot = constraints.clone();
        for terms in constraints.values_mut() {
            let mut needs_expand = false;
            for &(c_m, mn, d_m, _) in terms.iter() {
                if snapshot.contains_key(&(c_m, mn, d_m)) {
                    needs_expand = true;
                    break;
                }
            }
            if !needs_expand {
                continue;
            }
            any_expanded = true;
            let mut accum: BTreeMap<(usize, usize, usize), f64> = BTreeMap::new();
            for &(c_m, mn, d_m, w) in terms.iter() {
                if let Some(sub_terms) = snapshot.get(&(c_m, mn, d_m)) {
                    for &(c_r, mr, d_r, w_sub) in sub_terms {
                        *accum.entry((c_r, mr, d_r)).or_insert(0.0) += w * w_sub;
                    }
                } else {
                    *accum.entry((c_m, mn, d_m)).or_insert(0.0) += w;
                }
            }
            *terms = accum
                .into_iter()
                .filter(|&(_, w)| w.abs() > 1e-16)
                .map(|((c, n, d), w)| (c, n, d, w))
                .collect();
        }
        if !any_expanded {
            return Ok(());
        }
    }
    Err(PlateError::BadConstraint {
        what: "cyclic slave constraint chain in tied shells",
    })
}

struct CompliantNodalPatch {
    /// For each translational axis a in 0..3: list of `(comp, node, dof_0..6, weight)` for $\Delta u_a$.
    u_rows: [Vec<(usize, usize, usize, f64)>; 3],
    /// Global 3x3 translational stiffness matrix $A_i \mathbf{K}_{g,i}^u$ [N/m].
    k_u: [[f64; 3]; 3],
    /// Optional global 3x3 translational damping matrix $A_i \mathbf{C}_{g,i}^u$ [N*s/m].
    c_u: Option<[[f64; 3]; 3]>,
    /// Optional rotational jump rows and global 3x3 rotational stiffness/damping matrices when `RotationCoupling::Full`.
    rot: Option<(
        [Vec<(usize, usize, usize, f64)>; 3],
        [[f64; 3]; 3],
        Option<[[f64; 3]; 3]>,
    )>,
}

fn frame_tensor_3x3(
    n: [f64; 3],
    t1: [f64; 3],
    t2: [f64; 3],
    val_n: f64,
    val_s1: f64,
    val_s2: f64,
    area_scale: f64,
) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            out[a][b] = area_scale
                * (val_n * n[a] * n[b] + val_s1 * t1[a] * t1[b] + val_s2 * t2[a] * t2[b]);
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn build_compliant_polyline_patches(
    components: &[ShellComponent<'_>],
    comp_normals: &[Vec<[f64; 3]>],
    c_m: usize,
    m_nodes: &[usize],
    m_offsets: &[[f64; 3]],
    c_s: usize,
    s_nodes: &[usize],
    s_points: &[[f64; 3]],
    ops: &PolylineMortarOperators,
    rotations: RotationCoupling,
    glue: &GlueLine,
    out: &mut Vec<CompliantNodalPatch>,
    compliant_dofs_count: &mut usize,
) -> Result<(), PlateError> {
    let n_s = s_nodes.len();
    for (si, &sn) in s_nodes.iter().enumerate() {
        let chord = if si == 0 {
            sub(s_points[1], s_points[0])
        } else if si == n_s - 1 {
            sub(s_points[n_s - 1], s_points[n_s - 2])
        } else {
            sub(s_points[si + 1], s_points[si - 1])
        };
        let c_len = norm3(chord);
        if !c_len.is_finite() || c_len <= 1e-12 {
            return Err(PlateError::BadConstraint {
                what: "degenerate tangent on compliant glue polyline",
            });
        }
        let t1 = scale3(chord, 1.0 / c_len);
        let mut n_interp = [0.0; 3];
        for (mj_slot, &mn) in m_nodes.iter().enumerate() {
            let w = ops.c[si][mj_slot];
            let nm = comp_normals[c_m][mn];
            n_interp[0] += w * nm[0];
            n_interp[1] += w * nm[1];
            n_interp[2] += w * nm[2];
        }
        // Orthogonalize n against t1 so [n, t1, t2] is an exact right-handed orthonormal frame.
        let n_proj = sub(n_interp, scale3(t1, dot(n_interp, t1)));
        let n_len = norm3(n_proj);
        let n_unit = if n_len > 1e-12 {
            scale3(n_proj, 1.0 / n_len)
        } else {
            comp_normals[c_m][m_nodes[0]]
        };
        let t2 = super::cross(n_unit, t1);
        let area_i = glue.width_m * ops.d_diag[si];

        append_compliant_nodal_patch(
            components,
            c_m,
            m_nodes,
            m_offsets,
            c_s,
            sn,
            &ops.c[si],
            n_unit,
            t1,
            t2,
            area_i,
            rotations,
            glue,
            out,
            compliant_dofs_count,
        );
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_compliant_nodal_patch(
    components: &[ShellComponent<'_>],
    c_m: usize,
    m_nodes: &[usize],
    m_offsets: &[[f64; 3]],
    c_s: usize,
    sn: usize,
    c_row: &[f64],
    n_unit: [f64; 3],
    t1: [f64; 3],
    t2: [f64; 3],
    area_i: f64,
    rotations: RotationCoupling,
    glue: &GlueLine,
    out: &mut Vec<CompliantNodalPatch>,
    compliant_dofs_count: &mut usize,
) {
    let dof_cap = match rotations {
        RotationCoupling::Full => 6,
        // SeamTangent refuses before compliant assembly (rigid-tie kinematics).
        RotationCoupling::TranslationsOnly | RotationCoupling::SeamTangent => 3,
    };
    for d in 0..dof_cap {
        if components[c_s].model.dof_map[6 * sn + d].is_some() {
            *compliant_dofs_count += 1;
        }
    }

    let u_rows: [Vec<(usize, usize, usize, f64)>; 3] = core::array::from_fn(|a| {
        let mut row = Vec::new();
        row.push((c_s, sn, a, 1.0));
        for (mj_slot, &mn) in m_nodes.iter().enumerate() {
            let c_ij = c_row[mj_slot];
            if c_ij.abs() <= 1e-16 {
                continue;
            }
            row.push((c_m, mn, a, -c_ij));
            let skew = skew_theta_cross_d(m_offsets[mj_slot]);
            for b in 0..3 {
                let s_ab = skew[a][b];
                if s_ab != 0.0 {
                    row.push((c_m, mn, 3 + b, -c_ij * s_ab));
                }
            }
        }
        row
    });

    let k_u = frame_tensor_3x3(
        n_unit,
        t1,
        t2,
        glue.normal_n_m3,
        glue.shear_n_m3[0],
        glue.shear_n_m3[1],
        area_i,
    );
    let c_u = glue
        .viscosity_n_s_m3
        .filter(|v| v.iter().any(|&x| x > 0.0))
        .map(|v| frame_tensor_3x3(n_unit, t1, t2, v[0], v[1], v[2], area_i));

    let rot = if rotations == RotationCoupling::Full {
        let rot_scale = area_i * (glue.width_m * glue.width_m / 12.0);
        let theta_rows: [Vec<(usize, usize, usize, f64)>; 3] = core::array::from_fn(|a| {
            let mut row = Vec::new();
            row.push((c_s, sn, 3 + a, 1.0));
            for (mj_slot, &mn) in m_nodes.iter().enumerate() {
                let c_ij = c_row[mj_slot];
                if c_ij.abs() <= 1e-16 {
                    continue;
                }
                row.push((c_m, mn, 3 + a, -c_ij));
            }
            row
        });
        let k_theta = frame_tensor_3x3(
            n_unit,
            t1,
            t2,
            glue.normal_n_m3,
            glue.shear_n_m3[0],
            glue.shear_n_m3[1],
            rot_scale,
        );
        let c_theta = glue
            .viscosity_n_s_m3
            .filter(|v| v.iter().any(|&x| x > 0.0))
            .map(|v| frame_tensor_3x3(n_unit, t1, t2, v[0], v[1], v[2], rot_scale));
        Some((theta_rows, k_theta, c_theta))
    } else {
        None
    };

    out.push(CompliantNodalPatch {
        u_rows,
        k_u,
        c_u,
        rot,
    });
}

fn project_jump_row_to_reduced(
    row: &[(usize, usize, usize, f64)],
    comp_dof_maps: &[Vec<Option<usize>>],
    comp_t: &[Vec<Vec<(usize, f64)>>],
) -> Vec<(usize, f64)> {
    let mut accum: BTreeMap<usize, f64> = BTreeMap::new();
    for &(c, n, d, w) in row {
        if let Some(local_free) = comp_dof_maps[c][6 * n + d] {
            for &(r_idx, t_w) in &comp_t[c][local_free] {
                *accum.entry(r_idx).or_insert(0.0) += w * t_w;
            }
        }
    }
    accum
        .into_iter()
        .filter(|&(_, v)| v.abs() > 1e-16)
        .collect()
}

fn dot_slice(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .fold(0.0, |acc, (&x, &y)| x.mul_add(y, acc))
}

fn accumulate_quadratic_3x3(
    r_rows: &[Vec<(usize, f64)>; 3],
    mat3: &[[f64; 3]; 3],
    coo_a: &mut Coo,
    coo_b: Option<&mut Coo>,
) {
    let mut entries: BTreeMap<(usize, usize), f64> = BTreeMap::new();
    for a in 0..3 {
        for b in 0..3 {
            let k_ab = mat3[a][b];
            if k_ab.abs() <= 1e-24 {
                continue;
            }
            for &(ri, wi) in &r_rows[a] {
                for &(rj, wj) in &r_rows[b] {
                    *entries.entry((ri, rj)).or_insert(0.0) += k_ab * wi * wj;
                }
            }
        }
    }
    let mut sym_entries: Vec<(usize, usize, f64)> = Vec::new();
    for (&(ri, rj), &val) in &entries {
        if ri <= rj {
            let rev = entries.get(&(rj, ri)).copied().unwrap_or(val);
            let sym = 0.5 * (val + rev);
            if sym.abs() > 0.0 {
                sym_entries.push((ri, rj, sym));
            }
        }
    }
    for &(ri, rj, sym) in &sym_entries {
        coo_a.push(ri, rj, sym);
        if ri != rj {
            coo_a.push(rj, ri, sym);
        }
    }
    if let Some(cb) = coo_b {
        for &(ri, rj, sym) in &sym_entries {
            cb.push(ri, rj, sym);
            if ri != rj {
                cb.push(rj, ri, sym);
            }
        }
    }
}

fn accumulate_compliant_patch(
    patch: &CompliantNodalPatch,
    comp_dof_maps: &[Vec<Option<usize>>],
    comp_t: &[Vec<Vec<(usize, f64)>>],
    k_coo: &mut Coo,
    k_glue_coo: &mut Coo,
    c_glue_coo: &mut Coo,
    has_glue_damping: &mut bool,
) {
    let r_u: [Vec<(usize, f64)>; 3] = core::array::from_fn(|a| {
        project_jump_row_to_reduced(&patch.u_rows[a], comp_dof_maps, comp_t)
    });
    accumulate_quadratic_3x3(&r_u, &patch.k_u, k_coo, Some(k_glue_coo));
    if let Some(c_u) = &patch.c_u {
        *has_glue_damping = true;
        accumulate_quadratic_3x3(&r_u, c_u, c_glue_coo, None);
    }
    if let Some((theta_rows, k_theta, c_theta)) = &patch.rot {
        let r_theta: [Vec<(usize, f64)>; 3] = core::array::from_fn(|a| {
            project_jump_row_to_reduced(&theta_rows[a], comp_dof_maps, comp_t)
        });
        accumulate_quadratic_3x3(&r_theta, k_theta, k_coo, Some(k_glue_coo));
        if let Some(c_t) = c_theta {
            *has_glue_damping = true;
            accumulate_quadratic_3x3(&r_theta, c_t, c_glue_coo, None);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn record_slave_nodal_couplings_filtered(
    c_m: usize,
    m_nodes: &[usize],
    m_offsets: &[[f64; 3]],
    c_s: usize,
    s_nodes: &[usize],
    c_mat: &[Vec<f64>],
    skip_slave_node: &[bool],
    rotations: RotationCoupling,
    seam_tangents: Option<&[[f64; 3]]>,
    out: &mut BTreeMap<(usize, usize, usize), Vec<(usize, usize, usize, f64)>>,
) {
    for (si, &sn) in s_nodes.iter().enumerate() {
        if skip_slave_node.get(si).copied().unwrap_or(false) {
            continue;
        }
        // Translations a in 0..3: u_{s,a} = sum_j C_{si,j} (u_{m,a} + (theta_m x d_j)_a)
        for a in 0..3 {
            let mut terms = Vec::new();
            for (mj_slot, &mn) in m_nodes.iter().enumerate() {
                let c_ij = c_mat[si][mj_slot];
                if c_ij.abs() <= 1e-16 {
                    continue;
                }
                terms.push((c_m, mn, a, c_ij));
                let skew = skew_theta_cross_d(m_offsets[mj_slot]);
                for b in 0..3 {
                    let s_ab = skew[a][b];
                    if s_ab != 0.0 {
                        terms.push((c_m, mn, 3 + b, c_ij * s_ab));
                    }
                }
            }
            out.insert((c_s, sn, a), terms);
        }
        // Rotations a in 0..3 (if RotationCoupling::Full): theta_{s,a} = sum_j C_{si,j} theta_{m,a}
        if rotations == RotationCoupling::Full {
            for a in 0..3 {
                let mut terms = Vec::new();
                for (mj_slot, &mn) in m_nodes.iter().enumerate() {
                    let c_ij = c_mat[si][mj_slot];
                    if c_ij.abs() <= 1e-16 {
                        continue;
                    }
                    terms.push((c_m, mn, 3 + a, c_ij));
                }
                out.insert((c_s, sn, 3 + a), terms);
            }
        }
        // SeamTangent: eliminate the slave's dominant tangent component k from
        // t.theta_s = t.sum_j C_ij theta_m,j; its other two rotations stay free.
        if let (RotationCoupling::SeamTangent, Some(tangents)) = (rotations, seam_tangents) {
            let t = tangents[si];
            let k = (1..3).fold(0, |k, a| if t[a].abs() > t[k].abs() { a } else { k });
            let mut terms = Vec::new();
            for (mj_slot, &mn) in m_nodes.iter().enumerate() {
                let c_ij = c_mat[si][mj_slot];
                if c_ij.abs() <= 1e-16 {
                    continue;
                }
                for a in (0..3).filter(|&a| t[a] != 0.0) {
                    terms.push((c_m, mn, 3 + a, c_ij * t[a] / t[k]));
                }
            }
            for a in (0..3).filter(|&a| a != k && t[a] != 0.0) {
                terms.push((c_s, sn, 3 + a, -t[a] / t[k]));
            }
            out.insert((c_s, sn, 3 + k), terms);
        }
    }
}

/// Closest point and barycentric coordinates $(\lambda_0, \lambda_1, \lambda_2)$ on triangle $(p_0, p_1, p_2)$.
fn closest_barycentric_triangle(
    p: [f64; 3],
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
) -> ([f64; 3], [f64; 3]) {
    let v0 = sub(p1, p0);
    let v1 = sub(p2, p0);
    let v2 = sub(p, p0);
    let d00 = dot(v0, v0);
    let d01 = dot(v0, v1);
    let d11 = dot(v1, v1);
    let d20 = dot(v2, v0);
    let d21 = dot(v2, v1);
    let denom = d00 * d11 - d01 * d01;
    if denom.abs() <= 1e-24 {
        return ([1.0, 0.0, 0.0], p0);
    }
    let mut v = (d11 * d20 - d01 * d21) / denom;
    let mut w = (d00 * d21 - d01 * d20) / denom;
    v = v.clamp(0.0, 1.0);
    w = w.clamp(0.0, 1.0);
    if v + w > 1.0 {
        let s = v + w;
        v /= s;
        w /= s;
    }
    let u = 1.0 - v - w;
    let proj = [
        u * p0[0] + v * p1[0] + w * p2[0],
        u * p0[1] + v * p1[1] + w * p2[1],
        u * p0[2] + v * p1[2] + w * p2[2],
    ];
    ([u, v, w], proj)
}

/// Accumulate $T_c^T A T_c$ into `target` with exact upper/lower symmetry.
fn accumulate_congruence(
    source: &Csr,
    t_rows: &[Vec<(usize, f64)>],
    target: &mut Coo,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), PlateError> {
    let n = source.nrows();
    for i in 0..n {
        if (i & 63) == 0 && cancelled() {
            return Err(PlateError::Cancelled);
        }
        let ti = &t_rows[i];
        if ti.is_empty() {
            continue;
        }
        let (cols, vals) = source.row(i);
        for (&j, &val) in cols.iter().zip(vals.iter()) {
            if j < i || val == 0.0 {
                continue;
            }
            let tj = &t_rows[j];
            if tj.is_empty() {
                continue;
            }
            if !val.is_finite() {
                return Err(PlateError::BadSection {
                    what: "non-finite entry encountered during tie congruence accumulation",
                });
            }
            if i == j {
                for (a_idx, &(ra, wa)) in ti.iter().enumerate() {
                    let diag_v = wa * val * wa;
                    target.push(ra, ra, diag_v);
                    for &(rb, wb) in ti.iter().skip(a_idx + 1) {
                        let off_v = wa * val * wb;
                        target.push(ra, rb, off_v);
                        target.push(rb, ra, off_v);
                    }
                }
            } else {
                // Symmetrize pair (i, j) and (j, i) so target(p, q) and target(q, p)
                // receive identical floating-point addends in identical order.
                for &(ra, wa) in ti {
                    for &(rb, wb) in tj {
                        let prod = wa * val * wb;
                        if ra == rb {
                            target.push(ra, ra, 2.0 * prod);
                        } else {
                            target.push(ra, rb, prod);
                            target.push(rb, ra, prod);
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PlateSection;
    use crate::shell::{ShellSupport, assemble_shell};

    fn sample_rect_shell(
        x0: f64,
        x1: f64,
        y0: f64,
        y1: f64,
        nx: usize,
        ny: usize,
        z_offset: f64,
    ) -> ShellMesh {
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            let y = y0 + (y1 - y0) * (j as f64 / ny as f64);
            for i in 0..=nx {
                let x = x0 + (x1 - x0) * (i as f64 / nx as f64);
                nodes.push([x, y, z_offset]);
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut tris = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                tris.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                tris.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
            }
        }
        ShellMesh::new(nodes, tris).expect("valid rectangular shell mesh")
    }

    #[test]
    fn g0_dual_basis_biorthogonality() {
        // Verify on a single reference segment and on a multi-segment non-matching polyline
        // that \int_\Gamma \psi_i \phi_j^s ds = \delta_{ij} D_{ii} to 1e-14.
        let master: Vec<[f64; 3]> = (0..5).map(|i| [0.0, i as f64 * 0.15, 0.0]).collect();
        let offsets = vec![[0.0; 3]; master.len()];
        let slave: Vec<[f64; 3]> = (0..8)
            .map(|i| [0.0, i as f64 * (0.60 / 7.0), 0.0])
            .collect();
        let budget = TieBudget {
            max_segments: 32,
            max_slave_dofs: 128,
            gap_tolerance_m: 1e-9,
            angle_tolerance_rad: 1e-6,
        };
        let mut used = 0;
        let ops =
            assemble_polyline_dual_mortar(&master, &offsets, &slave, &budget, &mut used, || false)
                .expect("dual mortar assembly succeeds");

        for i in 0..slave.len() {
            assert!(ops.d_diag[i] > 0.0);
            for j in 0..slave.len() {
                let expected = if i == j { ops.d_diag[i] } else { 0.0 };
                let err = (ops.d_full[i][j] - expected).abs();
                assert!(
                    err < 1e-14,
                    "biorthogonality error at ({i}, {j}): got {}, expected {expected}, err {err:e}",
                    ops.d_full[i][j]
                );
            }
        }
    }

    #[test]
    fn g0_mortar_reproduces_linear_field_on_curved_interface() {
        // Circular-arc interface of radius R = 0.30 m over [0, pi/3] with 5-node master
        // and 7-node slave discretization.
        let radius = 0.30_f64;
        let arc_angle = core::f64::consts::PI / 3.0;
        let n_master = 5usize;
        let n_slave = 7usize;
        let master: Vec<[f64; 3]> = (0..n_master)
            .map(|j| {
                let theta = arc_angle * (j as f64 / (n_master - 1) as f64);
                [radius * det::cos(theta), radius * det::sin(theta), 0.0]
            })
            .collect();
        let offsets = vec![[0.0; 3]; n_master];
        let slave: Vec<[f64; 3]> = (0..n_slave)
            .map(|i| {
                let theta = arc_angle * (i as f64 / (n_slave - 1) as f64);
                [radius * det::cos(theta), radius * det::sin(theta), 0.0]
            })
            .collect();

        let budget = TieBudget {
            max_segments: 32,
            max_slave_dofs: 128,
            gap_tolerance_m: 0.01,
            angle_tolerance_rad: 0.35,
        };
        let mut used = 0;
        let ops =
            assemble_polyline_dual_mortar(&master, &offsets, &slave, &budget, &mut used, || false)
                .expect("curved dual mortar assembly succeeds");

        // Linear displacement field u(x) = A x + b.
        let a_mat = [
            [0.25, -0.40, 0.15],
            [0.70, 0.10, -0.35],
            [-0.20, 0.55, 0.30],
        ];
        let b_vec = [0.12, -0.08, 0.05];
        let eval_linear = |p: [f64; 3]| -> [f64; 3] {
            [
                a_mat[0][0] * p[0] + a_mat[0][1] * p[1] + a_mat[0][2] * p[2] + b_vec[0],
                a_mat[1][0] * p[0] + a_mat[1][1] * p[1] + a_mat[1][2] * p[2] + b_vec[1],
                a_mat[2][0] * p[0] + a_mat[2][1] * p[1] + a_mat[2][2] * p[2] + b_vec[2],
            ]
        };

        let u_master: Vec<[f64; 3]> = master.iter().map(|&p| eval_linear(p)).collect();
        let mut max_err = 0.0_f64;
        for (i, &sp) in slave.iter().enumerate() {
            let expected = eval_linear(sp);
            let mut reproduced = [0.0; 3];
            for (j, &um) in u_master.iter().enumerate() {
                let w = ops.c[i][j];
                for d in 0..3 {
                    reproduced[d] += w * um[d];
                }
            }
            for d in 0..3 {
                max_err = max_err.max((reproduced[d] - expected[d]).abs());
            }
        }
        assert!(
            max_err < 1e-12,
            "curved linear field reproduction error {max_err:e} >= 1e-12"
        );
    }

    #[test]
    fn g0_crosspoint_double_constraint_refuses() {
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 2, 0.0);
        let m1 = sample_rect_shell(0.2, 0.4, 0.0, 0.2, 2, 2, 0.0);
        let m2 = sample_rect_shell(0.2, 0.4, 0.0, 0.2, 2, 2, 0.0);
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let mod2 = assemble_shell(&m2, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let h2 = vec![0.003; m2.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
            ShellComponent {
                id: 2,
                mesh: &m2,
                model: &mod2,
                nodal_thickness_m: &h2,
            },
        ];
        // Declare component 2's seam nodes [0, 3, 6] as slave in TWO ties.
        let ties = [
            EdgeTie {
                a: TieSide {
                    component: 0,
                    nodes: vec![2, 5, 8],
                },
                b: TieSide {
                    component: 2,
                    nodes: vec![0, 3, 6],
                },
                offset: TieOffset::Midsurface,
                rotations: RotationCoupling::Full,
                kind: TieKind::Rigid,
            },
            EdgeTie {
                a: TieSide {
                    component: 1,
                    nodes: vec![0, 3, 6],
                },
                b: TieSide {
                    component: 2,
                    nodes: vec![0, 3, 6],
                },
                offset: TieOffset::Midsurface,
                rotations: RotationCoupling::Full,
                kind: TieKind::Rigid,
            },
        ];
        let budget = TieBudget {
            max_segments: 16,
            max_slave_dofs: 64,
            gap_tolerance_m: 1e-6,
            angle_tolerance_rad: 1e-4,
        };
        let err = tie_shells(&comps, &ties, budget, || false).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }

    #[test]
    fn g0_compliant_kind_refused_until_w2b() {
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 2, 0.0);
        let m1 = sample_rect_shell(0.2, 0.4, 0.0, 0.2, 2, 2, 0.0);
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        // Non-positive glue stiffness or width must refuse cleanly.
        let bad_glue = GlueLine {
            normal_n_m3: 0.0,
            shear_n_m3: [1e10, 1e10],
            viscosity_n_s_m3: None,
            width_m: 0.003,
        };
        let ties = [EdgeTie {
            a: TieSide {
                component: 0,
                nodes: vec![2, 5, 8],
            },
            b: TieSide {
                component: 1,
                nodes: vec![0, 3, 6],
            },
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Compliant(bad_glue),
        }];
        let budget = TieBudget {
            max_segments: 16,
            max_slave_dofs: 64,
            gap_tolerance_m: 1e-6,
            angle_tolerance_rad: 1e-4,
        };
        let err = tie_shells(&comps, &ties, budget, || false).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }

    #[test]
    fn g0_seam_tangent_needs_an_edge_seam() {
        // A lap joint has no seam tangent, so SeamTangent refuses rather than
        // picking an axis.
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 2, 0.0);
        let m1 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 2, 0.003);
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        let face = FaceTie {
            master_component: 0,
            master_tris: (0..m0.tris.len()).collect(),
            slave_component: 1,
            slave_tris: (0..m1.tris.len()).collect(),
            offset: TieOffset::Skin { side: 1 },
            rotations: RotationCoupling::SeamTangent,
            kind: TieKind::Rigid,
        };
        let budget = TieBudget {
            max_segments: 16,
            max_slave_dofs: 64,
            gap_tolerance_m: 1e-6,
            angle_tolerance_rad: 1e-4,
        };
        let err = tie_shells_with_faces(&comps, &[], &[face], budget, || false).unwrap_err();
        assert_eq!(
            err,
            PlateError::BadConstraint {
                what: "RotationCoupling::SeamTangent needs a 1-D seam; a FaceTie has none",
            }
        );
    }

    #[test]
    fn g0_glue_line_stiffness_spd_and_rigid_motion_free() {
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 2, 0.0);
        let m1 = sample_rect_shell(0.2, 0.4, 0.0, 0.2, 2, 3, 0.0);
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        let glue = GlueLine {
            normal_n_m3: 2.5e11,
            shear_n_m3: [1.2e11, 1.5e11],
            viscosity_n_s_m3: Some([2.5e5, 1.2e5, 1.5e5]),
            width_m: 0.003,
        };
        let ties = [EdgeTie {
            a: TieSide {
                component: 0,
                nodes: vec![2, 5, 8],
            },
            b: TieSide {
                component: 1,
                nodes: vec![0, 3, 6, 9],
            },
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Compliant(glue),
        }];
        let budget = TieBudget {
            max_segments: 32,
            max_slave_dofs: 64,
            gap_tolerance_m: 1e-6,
            angle_tolerance_rad: 1e-4,
        };
        let tied =
            tie_shells(&comps, &ties, budget, || false).expect("compliant glue tie succeeds");
        assert!(
            tied.slave_couplings().is_empty(),
            "compliant tie eliminates no DOFs"
        );

        // 1. Check symmetry and positive definiteness on any relative interface jump vector.
        let mut jump_q = vec![0.0; tied.model.free];
        for &sn in &[0usize, 3, 6, 9] {
            for d in 0..6 {
                if let Some(r_idx) = tied.map[1][sn][d] {
                    jump_q[r_idx] = 1e-3 * (1.0 + 0.25 * (d as f64));
                }
            }
        }
        let e_jump = tied.glue_strain_energy(&jump_q).unwrap();
        let d_jump = tied.glue_dissipation_rate(&jump_q).unwrap();
        assert!(
            e_jump > 0.0,
            "glue stiffness must be strictly positive on interface jumps: {e_jump}"
        );
        assert!(
            d_jump > 0.0,
            "glue viscosity must be strictly positive on interface jumps: {d_jump}"
        );
        let e_scale = 1.0e6 * e_jump;
        let d_scale = 1.0e6 * d_jump;

        // 2. Check all 6 common 3-D rigid-body motions (3 translations + 3 rotations) produce zero glue energy <= 1e-14 * e_scale.
        let motions: [([f64; 3], [f64; 3]); 7] = [
            ([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            ([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]),
            ([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]),
            ([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            ([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            ([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            ([0.37, -0.51, 0.83], [0.41, -0.29, 0.67]),
        ];
        for (u0, omega) in motions {
            let mut q_rigid = vec![0.0; tied.model.free];
            for (c, comp) in comps.iter().enumerate() {
                for (n, &pos) in comp.mesh.nodes.iter().enumerate() {
                    let u = add3(u0, crate::shell::cross(omega, pos));
                    let vals = [u[0], u[1], u[2], omega[0], omega[1], omega[2]];
                    for d in 0..6 {
                        if let Some(r_idx) = tied.map[c][n][d] {
                            q_rigid[r_idx] = vals[d];
                        }
                    }
                }
            }
            let e_rigid = tied.glue_strain_energy(&q_rigid).unwrap();
            let d_rigid = tied.glue_dissipation_rate(&q_rigid).unwrap();
            assert!(
                e_rigid.abs() <= 1e-14 * e_scale,
                "rigid motion ({u0:?}, {omega:?}) glue energy {e_rigid:.3e} > 1e-14 * {e_scale:.3e}"
            );
            assert!(
                d_rigid.abs() <= 1e-14 * d_scale,
                "rigid motion ({u0:?}, {omega:?}) glue dissipation {d_rigid:.3e} > 1e-14 * {d_scale:.3e}"
            );
        }
    }

    #[test]
    fn g0_glue_budget_cap_and_cap_plus_one() {
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 3, 0.0);
        let m1 = sample_rect_shell(0.2, 0.4, 0.0, 0.2, 2, 3, 0.0);
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        let glue = GlueLine {
            normal_n_m3: 1.0e11,
            shear_n_m3: [5.0e10, 5.0e10],
            viscosity_n_s_m3: None,
            width_m: 0.003,
        };
        // 4 seam nodes -> N = 3 segments, M = 24 compliant interface DOFs.
        let ties = [EdgeTie {
            a: TieSide {
                component: 0,
                nodes: vec![2, 5, 8, 11],
            },
            b: TieSide {
                component: 1,
                nodes: vec![0, 3, 6, 9],
            },
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Compliant(glue),
        }];
        let n_segments = 3usize;
        let n_glue_dofs = 24usize;
        let ok_budget = TieBudget {
            max_segments: n_segments,
            max_slave_dofs: n_glue_dofs,
            gap_tolerance_m: 1e-6,
            angle_tolerance_rad: 1e-4,
        };
        let tied =
            tie_shells(&comps, &ties, ok_budget, || false).expect("passes at exact glue cap");
        assert_eq!(tied.receipt.segments, n_segments);
        assert_eq!(tied.receipt.slave_dofs, n_glue_dofs);

        let over_seg_budget = TieBudget {
            max_segments: n_segments - 1,
            ..ok_budget
        };
        let err_seg = tie_shells(&comps, &ties, over_seg_budget, || false).unwrap_err();
        assert!(matches!(
            err_seg,
            PlateError::Budget {
                what: "tie_segments",
                ..
            }
        ));

        let over_dof_budget = TieBudget {
            max_slave_dofs: n_glue_dofs - 1,
            ..ok_budget
        };
        let err_dof = tie_shells(&comps, &ties, over_dof_budget, || false).unwrap_err();
        assert!(matches!(
            err_dof,
            PlateError::Budget {
                what: "tie_slave_dofs",
                ..
            }
        ));
    }

    #[test]
    fn g0_tie_budget_cap_and_cap_plus_one() {
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 3, 0.0);
        let m1 = sample_rect_shell(0.2, 0.4, 0.0, 0.2, 2, 3, 0.0);
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        // 4 seam nodes -> N = 3 segments, M = 24 slave DOFs.
        let ties = [EdgeTie {
            a: TieSide {
                component: 0,
                nodes: vec![2, 5, 8, 11],
            },
            b: TieSide {
                component: 1,
                nodes: vec![0, 3, 6, 9],
            },
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Rigid,
        }];
        let n_segments = 3usize;
        let n_slave_dofs = 24usize;
        let ok_budget = TieBudget {
            max_segments: n_segments,
            max_slave_dofs: n_slave_dofs,
            gap_tolerance_m: 1e-6,
            angle_tolerance_rad: 1e-4,
        };
        let tied = tie_shells(&comps, &ties, ok_budget, || false).expect("passes at exact cap N");
        assert_eq!(tied.receipt.segments, n_segments);
        assert_eq!(tied.receipt.slave_dofs, n_slave_dofs);

        let over_seg_budget = TieBudget {
            max_segments: n_segments - 1,
            ..ok_budget
        };
        let err_seg = tie_shells(&comps, &ties, over_seg_budget, || false).unwrap_err();
        assert!(matches!(
            err_seg,
            PlateError::Budget {
                what: "tie_segments",
                ..
            }
        ));

        let over_dof_budget = TieBudget {
            max_slave_dofs: n_slave_dofs - 1,
            ..ok_budget
        };
        let err_dof = tie_shells(&comps, &ties, over_dof_budget, || false).unwrap_err();
        assert!(matches!(
            err_dof,
            PlateError::Budget {
                what: "tie_slave_dofs",
                ..
            }
        ));
    }

    #[test]
    fn g0_gap_beyond_tolerance_refuses() {
        let sec = PlateSection::isotropic(10e9, 0.3, 0.003, 500.0).unwrap();
        let gap_tol = 1e-4;
        let m0 = sample_rect_shell(0.0, 0.2, 0.0, 0.2, 2, 2, 0.0);
        // Displace slave polyline by 2 * gap_tolerance_m in x.
        let m1 = sample_rect_shell(
            0.2 + 2.0 * gap_tol,
            0.4 + 2.0 * gap_tol,
            0.0,
            0.2,
            2,
            2,
            0.0,
        );
        let mod0 = assemble_shell(&m0, &sec, &[], ShellSupport::Free).unwrap();
        let mod1 = assemble_shell(&m1, &sec, &[], ShellSupport::Free).unwrap();
        let h0 = vec![0.003; m0.nodes.len()];
        let h1 = vec![0.003; m1.nodes.len()];
        let comps = [
            ShellComponent {
                id: 0,
                mesh: &m0,
                model: &mod0,
                nodal_thickness_m: &h0,
            },
            ShellComponent {
                id: 1,
                mesh: &m1,
                model: &mod1,
                nodal_thickness_m: &h1,
            },
        ];
        let ties = [EdgeTie {
            a: TieSide {
                component: 0,
                nodes: vec![2, 5, 8],
            },
            b: TieSide {
                component: 1,
                nodes: vec![0, 3, 6],
            },
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Rigid,
        }];
        let budget = TieBudget {
            max_segments: 16,
            max_slave_dofs: 64,
            gap_tolerance_m: gap_tol,
            angle_tolerance_rad: 1e-4,
        };
        let err = tie_shells(&comps, &ties, budget, || false).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }
}
