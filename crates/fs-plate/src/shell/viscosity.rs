//! Physical Kelvin--Voigt loss on the original CST/DKT shell coordinates.
//!
//! Integrates membrane A [N s/m], coupling B [N s], and bending D [N m s]
//! with the original CST strain and three-point DKT curvature rule. For the
//! physical thickness strain `epsilon(z) = epsilon0 + z*kappa`, the native
//! DKT field `beta=-grad(w)` already gives the signed first-moment curvature.
//! Retaining this sign preserves offset-layer B work.
//! Drilling stabilization, prestress, stiffeners, and inertia acquire no loss.
//!
//! The same constitutive Gram factors supply sparse spatial assembly and
//! direct projection. Projection forms strain before constitutive channels;
//! it never subtracts elastic matrices or erases off-diagonal modal damping.

use super::{DKT_SHELL_SIGNS, DKT_SHELL_SLOTS, FacetGeometry, ShellMesh, ShellModel};
use crate::{PlateError, dkt_curvature_quadrature};
use fs_math::{
    det,
    eft::{two_prod, two_sum},
};
use fs_sparse::{Coo, Csr};

/// Source-bound hereditary material response on the original CST/DKT shell.
pub mod harmonic;

fn bad(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// One facet's passive viscous resultants in its original tangent frame.
/// No frequency, temperature, source applicability or material fit is inferred.
#[derive(Clone, Debug)]
pub struct ShellViscosity {
    membrane: [f64; 9],
    coupling: [f64; 9],
    bending: [f64; 9],
    factor: [[f64; 6]; 6],
}

impl ShellViscosity {
    /// Independent membrane and bending loss. Exact zero and directional-only
    /// tensors are legal; no artificial positive damping floor is introduced.
    pub fn new(membrane: [f64; 9], bending: [f64; 9]) -> Result<Self, PlateError> {
        Self::new_coupled(membrane, [0.0; 9], bending)
    }

    /// Complete passive A/B/D law, with engineering shear in each third slot.
    /// B is symmetric, signed, and checked as part of the complete 6x6 block;
    /// reversing layer order reverses B without changing material passivity.
    /// Resultants follow `epsilon(z)=epsilon0+z*kappa` in the oriented facet.
    pub fn new_coupled(
        membrane: [f64; 9],
        coupling: [f64; 9],
        bending: [f64; 9],
    ) -> Result<Self, PlateError> {
        let factor = super::mitc::factor_coupled_block(&membrane, &coupling, &bending, true)?;
        Ok(Self {
            membrane,
            coupling,
            bending,
            factor,
        })
    }

    /// Original A/B/D resultants, without rotation, averaging or diagonalization.
    pub const fn matrices(&self) -> ([f64; 9], [f64; 9], [f64; 9]) {
        (self.membrane, self.coupling, self.bending)
    }

    /// Explicit source zero; numerical stabilization never changes this value.
    pub fn is_zero(&self) -> bool {
        self.membrane
            .iter()
            .chain(&self.coupling)
            .chain(&self.bending)
            .all(|x| *x == 0.0)
    }

    pub(super) fn channels(&self, strain: &[f64; 6]) -> [f64; 6] {
        core::array::from_fn(|column| {
            compensated_products((0..6).map(|row| (self.factor[row][column], strain[row])))
        })
    }
}

/// Explicit cold-work and allocation ceilings, checked before allocation.
#[derive(Clone, Copy, Debug)]
pub struct ShellViscosityBudget {
    /// Maximum triangle count.
    pub max_elements: usize,
    /// Maximum coordinate count before support elimination (six per node).
    pub max_full_dofs: usize,
    /// Maximum free spatial coordinate count.
    pub max_dofs: usize,
    /// Maximum supplied vectors in a direct modal projection.
    pub max_modes: usize,
    /// Conservative operation allowance; see [`projection_terms`]. Sparse
    /// assembly uses the same estimate with 18 local basis vectors per facet.
    pub max_work: usize,
}

/// Conservative cold-operation count including validation, original geometry,
/// strain/channel contractions and the complete symmetric output. Not a timing
/// promise. Overflow returns `None` before any output allocation.
pub fn projection_terms(full_dofs: usize, elements: usize, modes: usize) -> Option<usize> {
    let square = modes.checked_mul(modes)?;
    let local = 4096usize
        .checked_add(8192usize.checked_mul(modes)?)?
        .checked_add(256usize.checked_mul(square)?)?;
    full_dofs
        .checked_mul(modes.checked_add(1)?)?
        .checked_add(square)?
        .checked_add(elements.checked_mul(local)?)
}

pub(super) fn admit(
    mesh: &ShellMesh,
    section_count: usize,
    model: &ShellModel,
    budget: ShellViscosityBudget,
    modes: usize,
) -> Result<(), PlateError> {
    let full = mesh
        .nodes
        .len()
        .checked_mul(6)
        .ok_or_else(|| bad("shell viscosity DOF overflow"))?;
    let work = projection_terms(full, mesh.tris.len(), modes)
        .ok_or_else(|| bad("shell viscosity work overflow"))?;
    for (what, requested, limit) in [
        (
            "shell viscosity elements",
            mesh.tris.len(),
            budget.max_elements,
        ),
        ("shell viscosity full DOFs", full, budget.max_full_dofs),
        ("shell viscosity free DOFs", model.free, budget.max_dofs),
        ("shell viscosity work", work, budget.max_work),
    ] {
        if requested > limit {
            return Err(PlateError::Budget {
                what,
                requested,
                limit,
            });
        }
    }
    if section_count != mesh.tris.len() {
        return Err(PlateError::SectionCount {
            expected: mesh.tris.len(),
            actual: section_count,
        });
    }
    if mesh.nodes.len() < 3
        || mesh.tris.is_empty()
        || mesh.nodes.iter().flatten().any(|v| !v.is_finite())
        || model.free == 0
        || model.free > full
        || model.dof_map.len() != full
        || model.k.nrows() != model.free
        || model.k.ncols() != model.free
        || model.m.nrows() != model.free
        || model.m.ncols() != model.free
    {
        return Err(bad(
            "shell viscosity requires the original finite CST/DKT mesh and support map",
        ));
    }
    let mut seen = vec![false; model.free];
    for &index in model.dof_map.iter().flatten() {
        if index >= model.free || seen[index] {
            return Err(bad("invalid shell viscosity support map"));
        }
        seen[index] = true;
    }
    if seen.iter().any(|present| !present) {
        return Err(bad("incomplete shell viscosity support map"));
    }
    Ok(())
}

// The original DKT owner uses beta=-grad(w), so its returned curvature already
// has the physical first-moment sign. Reversing it would corrupt B coupling.
fn strain_rows(g: &FacetGeometry, element: usize) -> Result<[[[f64; 18]; 6]; 3], PlateError> {
    let membrane = super::membrane_strain_rows(g);
    let (curvature, _) = dkt_curvature_quadrature(&g.x, &g.y, element)?;
    let mut rows = [[[0.0; 18]; 6]; 3];
    for (point, bending) in curvature.iter().enumerate() {
        rows[point][..3].copy_from_slice(&membrane);
        for p in 0..3 {
            for j in 0..9 {
                rows[point][p + 3][DKT_SHELL_SLOTS[j]] = DKT_SHELL_SIGNS[j] * bending[p][j];
            }
        }
    }
    Ok(rows)
}

fn compensated_products(terms: impl Iterator<Item = (f64, f64)>) -> f64 {
    let (mut sum, mut correction) = (0.0, 0.0);
    for (a, b) in terms {
        let (product, product_error) = two_prod(a, b);
        let (next, sum_error) = two_sum(sum, product);
        sum = next;
        correction += product_error + sum_error;
    }
    sum + correction
}

// Shared CST/DKT kinematics for instantaneous viscosity and hereditary excess
// stiffness. The material owner only sees physical strain: drilling,
// equilibrium prestress and numerical stabilization never enter either law.
pub(super) fn projected_strains(
    g: &FacetGeometry,
    element: usize,
    tri: &[usize; 3],
    model: &ShellModel,
    modes: &[Vec<f64>],
    output: &mut [[[f64; 6]; 3]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), PlateError> {
    let rows = strain_rows(g, element)?;
    for (mode, strain) in modes.iter().zip(output) {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let nodal: [f64; 18] = core::array::from_fn(|i| {
            let base = 6 * tri[i / 6] + 3 * ((i % 6) / 3);
            compensated_products((0..3).map(|p| {
                (
                    g.frame[i % 3][p],
                    model.dof_map[base + p].map_or(0.0, |index| mode[index]),
                )
            }))
        });
        for (out, point) in strain.iter_mut().zip(&rows) {
            *out = core::array::from_fn(|r| {
                compensated_products((0..18).map(|i| (point[r][i], nodal[i])))
            });
            if out.iter().any(|value| !value.is_finite()) {
                return Err(bad("shell projected physical strain is nonfinite"));
            }
        }
    }
    Ok(())
}

/// Assemble the complete sparse physical viscous operator on a bare CST/DKT
/// shell's original supported global translation/rotation coordinates.
///
/// The caller supplies the SAME mesh used to assemble `model`. `ShellModel`
/// does not retain geometry for comparison: matching dimensions alone do not
/// establish identity. MITC oblique supports, tied assemblies and nonlinear
/// changed trial spaces must use their own physical expansion owners.
///
/// No K/M entry, support or drilling penalty is changed. Finite/budget/map
/// admission precedes allocation; cancellation is polled at every element.
pub fn assemble_constitutive_operator(
    mesh: &ShellMesh,
    sections: &[ShellViscosity],
    model: &ShellModel,
    budget: ShellViscosityBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Csr, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    admit(mesh, sections.len(), model, budget, 18)?;
    let mut coo = Coo::new(model.free, model.free);
    for (element, (tri, law)) in mesh.tris.iter().zip(sections).enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let g = mesh.facet(element)?;
        if law.is_zero() {
            continue;
        }
        let rows = strain_rows(&g, element)?;
        let scale = det::sqrt(g.area_m2) * det::sqrt(1.0 / 3.0);
        let mut local = [[0.0; 18]; 18];
        for point in rows {
            let channels: [[f64; 6]; 18] = core::array::from_fn(|i| {
                let base = 6 * (i / 6) + 3 * ((i % 6) / 3);
                let strain = core::array::from_fn(|r| {
                    compensated_products((0..3).map(|p| (g.frame[p][i % 3], point[r][base + p])))
                });
                law.channels(&strain).map(|x| scale * x)
            });
            if channels.iter().flatten().any(|v| !v.is_finite()) {
                return Err(bad("shell viscosity element channels overflow"));
            }
            for i in 0..18 {
                for j in i..18 {
                    local[i][j] +=
                        compensated_products((0..6).map(|k| (channels[i][k], channels[j][k])));
                    if !local[i][j].is_finite() {
                        return Err(bad("shell viscosity element overflow"));
                    }
                }
            }
        }
        for i in 0..18 {
            let Some(row) = model.dof_map[6 * tri[i / 6] + i % 6] else {
                continue;
            };
            for j in i..18 {
                let Some(column) = model.dof_map[6 * tri[j / 6] + j % 6] else {
                    continue;
                };
                if local[i][j] != 0.0 {
                    coo.push(row, column, local[i][j]);
                    if i != j {
                        coo.push(column, row, local[i][j]);
                    }
                }
            }
        }
    }
    let matrix = coo.assemble();
    if (0..model.free).any(|r| matrix.row(r).1.iter().any(|v| !v.is_finite())) {
        return Err(bad("shell viscosity assembly overflow"));
    }
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    Ok(matrix)
}

/// Project every material cross term on caller-supplied reduced-coordinate
/// vectors. The result is row-major `Phi^T C Phi`. It has units 1/s only for
/// mass-normalized modes; this function makes no eigenpair or mass claim.
///
/// The original global-axis support map expands each vector; mesh identity and
/// bare CST/DKT restrictions are as in [`assemble_constitutive_operator`].
/// Strain is evaluated before material factors, preserving resolved weak loss
/// and null channels without cancellation in a dense nodal matrix contraction.
/// No projection or strain rows survive a cancelled or refused call.
pub fn project_viscosity(
    mesh: &ShellMesh,
    sections: &[ShellViscosity],
    model: &ShellModel,
    modes: &[Vec<f64>],
    budget: ShellViscosityBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<f64>, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    let n = modes.len();
    if n == 0 || n > budget.max_modes || n > super::mitc::viscosity::MAX_VISCOSITY_MODES {
        return Err(bad(
            "shell viscosity requires a nonempty bounded projection basis",
        ));
    }
    admit(mesh, sections.len(), model, budget, n)?;
    if modes
        .iter()
        .any(|mode| mode.len() != model.free || mode.iter().any(|v| !v.is_finite()))
    {
        return Err(bad(
            "shell viscosity requires finite vectors on the original free coordinates",
        ));
    }
    let mut out = vec![0.0; n * n];
    let mut strains = vec![[[0.0; 6]; 3]; n];
    let mut channels = vec![[0.0; 6]; n];
    for (element, (tri, law)) in mesh.tris.iter().zip(sections).enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let g = mesh.facet(element)?;
        if law.is_zero() {
            continue;
        }
        let scale = det::sqrt(g.area_m2) * det::sqrt(1.0 / 3.0);
        projected_strains(&g, element, tri, model, modes, &mut strains, &mut cancelled)?;
        for point in 0..3 {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            for (strain, channel) in strains.iter().zip(&mut channels) {
                *channel = law.channels(&strain[point]).map(|x| scale * x);
                if channel.iter().any(|v| !v.is_finite()) {
                    return Err(bad("shell viscosity projected strain channels overflow"));
                }
            }
            for i in 0..n {
                if cancelled() {
                    return Err(PlateError::Cancelled);
                }
                for j in i..n {
                    out[n * i + j] +=
                        compensated_products((0..6).map(|k| (channels[i][k], channels[j][k])));
                    if !out[n * i + j].is_finite() {
                        return Err(bad("shell modal viscosity overflow"));
                    }
                    out[n * j + i] = out[n * i + j];
                }
            }
        }
    }
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    Ok(out)
}

// A constitutive increment or frequency derivative is signed, not a passive
// viscosity. Integrate it on the SAME physical strain rows and support map;
// never project it through a positive factor or a separate shell operator.
fn assemble_signed_constitutive_operator(
    mesh: &ShellMesh,
    sections: &[([f64; 9], [f64; 9], [f64; 9])],
    model: &ShellModel,
    budget: ShellViscosityBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Csr, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    admit(mesh, sections.len(), model, budget, 18)?;
    let mut coo = Coo::new(model.free, model.free);
    for (element, (tri, &(a, b, d))) in mesh.tris.iter().zip(sections).enumerate() {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        for matrix in [&a, &b, &d] {
            if matrix.iter().any(|v| !v.is_finite())
                || (0..3).any(|i| (0..i).any(|j| matrix[3*i+j] != matrix[3*j+i]))
            {
                return Err(bad("signed shell resultants must be finite and symmetric"));
            }
        }
        let g = mesh.facet(element)?;
        let rows = strain_rows(&g, element)?;
        let mut local = [[0.0; 18]; 18];
        for point in rows {
            let strains: [[f64; 6]; 18] = core::array::from_fn(|i| {
                let base = 6 * (i / 6) + 3 * ((i % 6) / 3);
                core::array::from_fn(|r| {
                    compensated_products((0..3).map(|p| (g.frame[p][i % 3], point[r][base + p])))
                })
            });
            for i in 0..18 {
                let stress: [f64; 6] = core::array::from_fn(|r| {
                    compensated_products((0..6).map(|s| {
                        let matrix = if r < 3 && s < 3 { &a }
                            else if r >= 3 && s >= 3 { &d } else { &b };
                        (matrix[3 * (r % 3) + s % 3], strains[i][s])
                    }))
                });
                for j in i..18 {
                    local[i][j] += (g.area_m2 / 3.0)
                        * compensated_products((0..6).map(|r| (stress[r], strains[j][r])));
                    if !local[i][j].is_finite() {
                        return Err(bad("signed shell constitutive assembly overflow"));
                    }
                }
            }
        }
        for i in 0..18 {
            let Some(row) = model.dof_map[6 * tri[i / 6] + i % 6] else { continue; };
            for j in i..18 {
                let Some(column) = model.dof_map[6 * tri[j / 6] + j % 6] else { continue; };
                if local[i][j] != 0.0 {
                    coo.push(row, column, local[i][j]);
                    if i != j { coo.push(column, row, local[i][j]); }
                }
            }
        }
    }
    let matrix = coo.assemble();
    if (0..model.free).any(|r| matrix.row(r).1.iter().any(|v| !v.is_finite())) {
        return Err(bad("signed shell constitutive assembly overflow"));
    }
    if cancelled() { return Err(PlateError::Cancelled); }
    Ok(matrix)
}
