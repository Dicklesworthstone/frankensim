//! Weak-form assembly on the FEEC 0-form (P₁ vertex-hat) space.
//!
//! # The discrete forms
//!
//! With `λ_a` the barycentric hat functions (constant gradients `g_a` on
//! each tet), `V` the element volume, and `A` a boundary-face area:
//!
//! ```text
//!   element conduction   K^e_ab = V · gᵀ_a K(T̄_e) g_b
//!   element source       b^e_a  = Σ_b  V (1 + δ_ab)/20 · f_b
//!   face Neumann         b^F_a −= Σ_b  A (1 + δ_ab)/12 · q_b
//!   face Robin (matrix)  R^F_ab = h̄_F · A (1 + δ_ab)/12
//!   face Robin (load)    b^F_a += h̄_F · Σ_b A (1 + δ_ab)/12 · T_ref,b
//! ```
//!
//! `∫ λ_a λ_b dV = V(1+δ_ab)/20` on a tet and `∫ λ_a λ_b dA =
//! A(1+δ_ab)/12` on a triangle are exact, so the source and boundary
//! rules are EXACT for data that is linear on the element/face and
//! second-order accurate otherwise. `h̄_F` is the face mean of `h`: the
//! Robin operator block is exact for a face-constant transfer
//! coefficient and second-order otherwise. All four rules are stated in
//! `CONTRACT.md` as the quadrature claim.
//!
//! `K` is evaluated ONCE per element at `T̄_e`, the element mean of the
//! P₁ iterate — which for P₁ is the EXACT element average of `T_h`, not
//! an approximation of it. For `k` linear in `T` the resulting element
//! integral `V·k(T̄_e)` equals `∫_e k(T_h) dV` exactly; for a curved
//! `k(T)` it is the midpoint rule, second-order in the element diameter.
//!
//! # Cancellation
//!
//! Element and face loops run in tiles of [`ASSEMBLY_TILE`] and poll
//! `Cx` at every tile boundary. A cancelled assembly returns
//! [`crate::ConductionError::Cancelled`] naming the stage and the tile
//! that was about to run; the partially staged `Coo` is dropped, so no
//! half-built operator can escape.

use fs_exec::Cx;
use fs_sparse::{Coo, Csr};

use crate::ConductionError;
use crate::bc::{ThermalBc, ThermalBoundary};
use crate::field::ScalarField;
use crate::interface::ThermalInterfaces;
use crate::material::ConductivityModel;
use crate::mesh::ConductionMesh;

/// Elements (or faces) processed between cancellation polls.
pub const ASSEMBLY_TILE: usize = 512;

/// A full (pre-elimination) conduction system.
#[derive(Debug, Clone)]
pub struct AssembledSystem {
    /// `n × n` operator over ALL vertices: conduction plus Robin.
    pub operator: Csr,
    /// Load vector over all vertices: source, Neumann, and Robin terms.
    pub load: Vec<f64>,
}

/// Free/prescribed degree-of-freedom bookkeeping. Built once per
/// boundary partition and reused across nonlinear iterations, so the
/// elimination pattern cannot drift between Newton steps.
#[derive(Debug, Clone, PartialEq)]
pub struct DofMap {
    free: Vec<usize>,
    slot: Vec<usize>,
    prescribed: Vec<f64>,
    fixed: Vec<usize>,
}

impl DofMap {
    /// Build from a boundary partition.
    ///
    /// # Errors
    /// [`ConductionError::NoFreeDofs`] when every vertex is pinned.
    pub fn new(boundary: &ThermalBoundary, vertex_count: usize) -> Result<DofMap, ConductionError> {
        let mut slot = vec![usize::MAX; vertex_count];
        let mut prescribed = vec![0.0f64; vertex_count];
        let mut fixed = Vec::with_capacity(boundary.dirichlet().len());
        for &(v, value) in boundary.dirichlet() {
            slot[v] = usize::MAX;
            prescribed[v] = value;
            fixed.push(v);
        }
        let mut free = Vec::with_capacity(vertex_count - fixed.len());
        for (v, s) in slot.iter_mut().enumerate() {
            if fixed.binary_search(&v).is_err() {
                *s = free.len();
                free.push(v);
            }
        }
        if free.is_empty() {
            return Err(ConductionError::NoFreeDofs);
        }
        Ok(DofMap {
            free,
            slot,
            prescribed,
            fixed,
        })
    }

    /// Number of unknowns.
    #[must_use]
    pub fn n(&self) -> usize {
        self.free.len()
    }

    /// The free vertices, ascending.
    #[must_use]
    pub fn free(&self) -> &[usize] {
        &self.free
    }

    /// The Dirichlet-pinned vertices, ascending.
    #[must_use]
    pub fn fixed(&self) -> &[usize] {
        &self.fixed
    }

    /// The prescribed value at every vertex (`0.0` on free vertices).
    #[must_use]
    pub fn prescribed(&self) -> &[f64] {
        &self.prescribed
    }

    /// Expand a free-dof vector into a full nodal temperature field.
    #[must_use]
    pub fn scatter(&self, free_values: &[f64]) -> Vec<f64> {
        let mut full = self.prescribed.clone();
        for (i, &v) in self.free.iter().enumerate() {
            full[v] = free_values[i];
        }
        full
    }

    /// Restrict a full nodal vector to the free dofs.
    #[must_use]
    pub fn gather(&self, full: &[f64]) -> Vec<f64> {
        self.free.iter().map(|&v| full[v]).collect()
    }

    /// The free-dof index of a vertex, or `None` when it is prescribed.
    #[must_use]
    pub fn slot_of(&self, vertex: usize) -> Option<usize> {
        let s = self.slot[vertex];
        (s != usize::MAX).then_some(s)
    }
}

fn checkpoint(cx: &Cx<'_>, stage: &'static str, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint()
        .map_err(|_| ConductionError::Cancelled { stage, at })
}

/// The 4×4 element conduction matrix `V · gᵀ_a K g_b` for one tet.
#[must_use]
pub fn element_stiffness(
    mesh: &ConductionMesh,
    element: usize,
    k: &[[f64; 3]; 3],
) -> [[f64; 4]; 4] {
    let volume = mesh.element_volume(element);
    let g = &mesh.geometry().grads[element];
    let mut out = [[0.0f64; 4]; 4];
    for a in 0..4 {
        // K·g_a, then dot with g_b.
        let mut kga = [0.0f64; 3];
        for (i, kg) in kga.iter_mut().enumerate() {
            *kg = k[i][0].mul_add(g[a][0], k[i][1].mul_add(g[a][1], k[i][2] * g[a][2]));
        }
        for b in 0..4 {
            let v = g[b][0].mul_add(kga[0], g[b][1].mul_add(kga[1], g[b][2] * kga[2]));
            out[a][b] = volume * v;
        }
    }
    out
}

/// Element-mean temperature — for P₁ this IS the exact element average
/// of the discrete field, not an approximation of it.
#[must_use]
pub fn element_temperature(mesh: &ConductionMesh, element: usize, temperature: &[f64]) -> f64 {
    let tet = mesh.complex().tets[element];
    let mut sum = 0.0f64;
    for &v in &tet {
        sum += temperature[v as usize];
    }
    sum / 4.0
}

/// Assemble `A(T)` and the load `b`.
///
/// # Errors
/// [`ConductionError::Cancelled`] at a tile boundary;
/// [`ConductionError::OutsideTemperatureSpan`] when the iterate leaves
/// the sampled conductivity span; [`ConductionError::FieldLength`] for a
/// mis-sized field.
pub fn assemble_operator(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    source: &ScalarField,
    temperature: &[f64],
) -> Result<AssembledSystem, ConductionError> {
    assemble_operator_scaled_with_interfaces(
        cx,
        mesh,
        boundary,
        material,
        source,
        temperature,
        None,
        None,
    )
}

/// [`assemble_operator`] with an explicitly bound matching-face contact set.
///
/// The contact block is temperature independent and symmetric positive
/// semidefinite.  Every coincident face pair was already checked as complete by
/// [`ThermalInterfaces::new`].
pub fn assemble_operator_with_interfaces(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    source: &ScalarField,
    temperature: &[f64],
    interfaces: &ThermalInterfaces,
) -> Result<AssembledSystem, ConductionError> {
    assemble_operator_scaled_with_interfaces(
        cx,
        mesh,
        boundary,
        material,
        source,
        temperature,
        None,
        Some(interfaces),
    )
}

/// [`assemble_operator`] with an optional per-element multiplier on the
/// conduction block: `K_e ← ρ_e · K(T̄_e)`.
///
/// This is the design-parameterized form the adjoint path differentiates
/// (and the shape a density/SIMP thermal topology study wants). `ρ` does
/// NOT scale the Robin block, the Neumann load, or the source: those are
/// boundary and body data, not conductivity.
///
/// # Errors
/// As [`assemble_operator`], plus [`ConductionError::FieldLength`] when
/// the scale array does not have one entry per element.
pub fn assemble_operator_scaled(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    source: &ScalarField,
    temperature: &[f64],
    element_scale: Option<&[f64]>,
) -> Result<AssembledSystem, ConductionError> {
    assemble_operator_scaled_with_interfaces(
        cx,
        mesh,
        boundary,
        material,
        source,
        temperature,
        element_scale,
        None,
    )
}

pub(crate) fn assemble_operator_scaled_with_interfaces(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    source: &ScalarField,
    temperature: &[f64],
    element_scale: Option<&[f64]>,
    interfaces: Option<&ThermalInterfaces>,
) -> Result<AssembledSystem, ConductionError> {
    if let Some(interfaces) = interfaces {
        interfaces.validate_for(mesh, boundary)?;
    } else {
        ThermalInterfaces::require_no_undeclared(mesh)?;
    }
    let n = mesh.vertex_count();
    if let Some(scale) = element_scale
        && scale.len() != mesh.element_count()
    {
        return Err(ConductionError::FieldLength {
            field: "element conductivity scale",
            expected: mesh.element_count(),
            found: scale.len(),
        });
    }
    if temperature.len() != n {
        return Err(ConductionError::FieldLength {
            field: "temperature iterate",
            expected: n,
            found: temperature.len(),
        });
    }
    source.validate("volumetric source", n)?;

    let mut coo = Coo::new(n, n);
    let mut load = vec![0.0f64; n];
    let ne = mesh.element_count();
    let mut start = 0usize;
    while start < ne {
        checkpoint(cx, "assemble-elements", start)?;
        let end = (start + ASSEMBLY_TILE).min(ne);
        for e in start..end {
            let tet = mesh.complex().tets[e];
            let t_e = element_temperature(mesh, e, temperature);
            let k = material.tensor_at(t_e)?;
            let mut ke = element_stiffness(mesh, e, &k);
            if let Some(scale) = element_scale {
                for row in &mut ke {
                    for v in row.iter_mut() {
                        *v *= scale[e];
                    }
                }
            }
            for a in 0..4 {
                for b in 0..4 {
                    coo.push(tet[a] as usize, tet[b] as usize, ke[a][b]);
                }
            }
            let volume = mesh.element_volume(e);
            for a in 0..4 {
                let mut acc = 0.0f64;
                for (b, &vb) in tet.iter().enumerate() {
                    let w = if a == b { volume / 10.0 } else { volume / 20.0 };
                    acc = w.mul_add(source.at(vb as usize), acc);
                }
                load[tet[a] as usize] += acc;
            }
        }
        start = end;
    }

    assemble_boundary(cx, mesh, boundary, &mut coo, &mut load)?;
    if let Some(interfaces) = interfaces {
        interfaces.assemble_into(cx, &mut coo)?;
    }

    Ok(AssembledSystem {
        operator: coo.assemble(),
        load,
    })
}

fn assemble_boundary(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    coo: &mut Coo,
    load: &mut [f64],
) -> Result<(), ConductionError> {
    let nf = mesh.boundary().len();
    let mut start = 0usize;
    while start < nf {
        checkpoint(cx, "assemble-boundary", start)?;
        let end = (start + ASSEMBLY_TILE).min(nf);
        for slot in start..end {
            let face = &mesh.boundary()[slot];
            let Some(condition) = boundary.condition_for(slot) else {
                continue;
            };
            let verts = [
                face.vertices[0] as usize,
                face.vertices[1] as usize,
                face.vertices[2] as usize,
            ];
            let area = face.area;
            match condition {
                ThermalBc::Dirichlet { .. } => {}
                ThermalBc::Neumann { outward_flux } => {
                    for a in 0..3 {
                        let mut acc = 0.0f64;
                        for (b, &vb) in verts.iter().enumerate() {
                            let w = if a == b { area / 6.0 } else { area / 12.0 };
                            acc = w.mul_add(outward_flux.at(vb), acc);
                        }
                        load[verts[a]] -= acc;
                    }
                }
                ThermalBc::Robin { htc, t_ref } => {
                    let h_bar = (htc.at(verts[0]) + htc.at(verts[1]) + htc.at(verts[2])) / 3.0;
                    for a in 0..3 {
                        let mut acc = 0.0f64;
                        for (b, &vb) in verts.iter().enumerate() {
                            let w = if a == b { area / 6.0 } else { area / 12.0 };
                            coo.push(verts[a], vb, h_bar * w);
                            acc = w.mul_add(t_ref.at(vb), acc);
                        }
                        load[verts[a]] += h_bar * acc;
                    }
                }
            }
        }
        start = end;
    }
    Ok(())
}

/// Assemble the Newton Jacobian `J = ∂R/∂T` of `R(T) = A(T)T − b`.
///
/// `J^e_ac = V (gᵀ_a K g_c) + (V/4) (gᵀ_a K'(T̄_e) ∇T_h)` plus the
/// (temperature-independent) Robin block. The second term is what makes
/// `J` NONSYMMETRIC, which is why the Newton path runs FGMRES rather
/// than CG.
///
/// # Errors
/// As [`assemble_operator`].
pub fn assemble_jacobian(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    temperature: &[f64],
) -> Result<Csr, ConductionError> {
    assemble_jacobian_with_optional_interfaces(cx, mesh, boundary, material, temperature, None)
}

/// [`assemble_jacobian`] with the constant matching-face contact block.
pub fn assemble_jacobian_with_interfaces(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    temperature: &[f64],
    interfaces: &ThermalInterfaces,
) -> Result<Csr, ConductionError> {
    assemble_jacobian_with_optional_interfaces(
        cx,
        mesh,
        boundary,
        material,
        temperature,
        Some(interfaces),
    )
}

pub(crate) fn assemble_jacobian_with_optional_interfaces(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    material: &ConductivityModel,
    temperature: &[f64],
    interfaces: Option<&ThermalInterfaces>,
) -> Result<Csr, ConductionError> {
    if let Some(interfaces) = interfaces {
        interfaces.validate_for(mesh, boundary)?;
    } else {
        ThermalInterfaces::require_no_undeclared(mesh)?;
    }
    let n = mesh.vertex_count();
    if temperature.len() != n {
        return Err(ConductionError::FieldLength {
            field: "temperature iterate",
            expected: n,
            found: temperature.len(),
        });
    }
    let nonlinear = material.is_temperature_dependent();
    let mut coo = Coo::new(n, n);
    let mut discard = vec![0.0f64; n];
    let ne = mesh.element_count();
    let mut start = 0usize;
    while start < ne {
        checkpoint(cx, "assemble-jacobian", start)?;
        let end = (start + ASSEMBLY_TILE).min(ne);
        for e in start..end {
            let tet = mesh.complex().tets[e];
            let t_e = element_temperature(mesh, e, temperature);
            let k = material.tensor_at(t_e)?;
            let ke = element_stiffness(mesh, e, &k);
            for a in 0..4 {
                for b in 0..4 {
                    coo.push(tet[a] as usize, tet[b] as usize, ke[a][b]);
                }
            }
            if !nonlinear {
                continue;
            }
            let kp = material.tensor_derivative_at(t_e)?;
            let g = &mesh.geometry().grads[e];
            let volume = mesh.element_volume(e);
            // ∇T_h on this element (constant).
            let mut grad_t = [0.0f64; 3];
            for (b, gb) in g.iter().enumerate() {
                let tb = temperature[tet[b] as usize];
                for i in 0..3 {
                    grad_t[i] = tb.mul_add(gb[i], grad_t[i]);
                }
            }
            // K'·∇T_h.
            let mut kpg = [0.0f64; 3];
            for (i, v) in kpg.iter_mut().enumerate() {
                *v = kp[i][0].mul_add(grad_t[0], kp[i][1].mul_add(grad_t[1], kp[i][2] * grad_t[2]));
            }
            for a in 0..4 {
                let ga_kpg = g[a][0].mul_add(kpg[0], g[a][1].mul_add(kpg[1], g[a][2] * kpg[2]));
                let contribution = volume * ga_kpg / 4.0;
                for c in 0..4 {
                    coo.push(tet[a] as usize, tet[c] as usize, contribution);
                }
            }
        }
        start = end;
    }
    assemble_boundary(cx, mesh, boundary, &mut coo, &mut discard)?;
    if let Some(interfaces) = interfaces {
        interfaces.assemble_into(cx, &mut coo)?;
    }
    Ok(coo.assemble())
}

/// Eliminate the Dirichlet rows/columns: returns `(A_ff, b_f)` with
/// `b_f = rhs_f − A_fd · T_D`.
#[must_use]
pub fn reduce(system: &AssembledSystem, dofs: &DofMap) -> (Csr, Vec<f64>) {
    let (a, b) = reduce_matrix_and_lift(&system.operator, dofs);
    let mut rhs = Vec::with_capacity(dofs.n());
    for (i, &v) in dofs.free().iter().enumerate() {
        rhs.push(system.load[v] + b[i]);
    }
    (a, rhs)
}

/// Eliminate the Dirichlet rows/columns of a bare matrix, returning the
/// reduced matrix and the lift `−A_fd · T_D`.
#[must_use]
pub fn reduce_matrix_and_lift(matrix: &Csr, dofs: &DofMap) -> (Csr, Vec<f64>) {
    let nf = dofs.n();
    let mut coo = Coo::new(nf, nf);
    let mut lift = vec![0.0f64; nf];
    let prescribed = dofs.prescribed();
    for (i, &v) in dofs.free().iter().enumerate() {
        let (cols, vals) = matrix.row(v);
        for (&c, &value) in cols.iter().zip(vals) {
            let s = dofs.slot[c];
            if s == usize::MAX {
                lift[i] -= value * prescribed[c];
            } else {
                coo.push(i, s, value);
            }
        }
    }
    (coo.assemble(), lift)
}

/// The free-dof residual `R = (A·T_full − load)|_free`.
#[must_use]
pub fn residual(system: &AssembledSystem, dofs: &DofMap, full_temperature: &[f64]) -> Vec<f64> {
    let n = full_temperature.len();
    let mut ax = vec![0.0f64; n];
    system.operator.spmv(full_temperature, &mut ax);
    dofs.free()
        .iter()
        .map(|&v| ax[v] - system.load[v])
        .collect()
}

/// The full-vector residual `A·T_full − load`, whose Dirichlet entries
/// are the nodal reaction fluxes used by the energy balance.
#[must_use]
pub fn full_residual(system: &AssembledSystem, full_temperature: &[f64]) -> Vec<f64> {
    let n = full_temperature.len();
    let mut ax = vec![0.0f64; n];
    system.operator.spmv(full_temperature, &mut ax);
    for (i, a) in ax.iter_mut().enumerate() {
        *a -= system.load[i];
    }
    ax
}
