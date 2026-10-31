//! Frozen membrane initial-stress stiffness of planar MITC panels.
//!
//! For the original P1 transverse displacement w, the added second variation
//! is integral Nxx*w_x^2 + 2*Nxy*w_x*w_y + Nyy*w_y^2 dA. The native compatible
//! thermal recovery supplies signed N, including the affine thermal director
//! bubble. P1 gradients are constant, so its native area integral is sufficient
//! even when membrane stress varies across the facet. Tension stiffens;
//! compression softens and may produce negative generalized eigenvalues.
//!
//! This chosen geometric term has no director/bubble coordinates. Its internal
//! blocks are zero, so the original elastic Schur complement and condensed
//! mass are unchanged: condense(K_material + K_geometric) = K_condensed + K_G.
//! No separately condensed geometric matrix or changed inertia is substituted.
//!
//! This is a bounded, flat-reference, small-displacement initial-stress model,
//! not a general nonlinear shell tangent. It excludes geometric bending/shear
//! stress terms, follower-load stiffness, evolving geometry and post-buckling.
//! Source thermal-strain/material validity remains with the original caller.
//! Signed matrices are exposed to the existing spectral owner; a nonpositive
//! eigenvalue is not converted into a real vibration frequency.

use super::{PreparedThermalShell, ThermalShellResponse};
use crate::shell::ShellModel;
use crate::shell::mitc::{MitcShellModel, PlateError, dot};
use fs_math::det;
use fs_sparse::{Coo, Csr};

/// Caller-selected small-kinematic envelope; both limits lie in (0, 0.1].
/// These are model limits, not measured failure strains or a stability proof.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MembranePrestressLimits {
    /// Frobenius norm of the complete P1 midsurface displacement gradient.
    pub maximum_displacement_gradient: f64,
    /// Physical nodal rotation and complete enriched director-motion bound.
    pub maximum_rotation_rad: f64,
}

/// Original supports/mass and a signed initial-stress stiffness correction.
/// The source model is retained by borrow, preventing an altered tangent from
/// being mistaken for a fresh elastic model or receiving this correction twice.
#[derive(Debug)]
pub struct MembranePrestressedShell<'a> {
    original: &'a MitcShellModel,
    pencil: ShellModel,
    geometric: Csr,
    state: ThermalShellResponse,
}
impl MembranePrestressedShell<'_> {
    /// Native sparse (K_material + K_G, M), in original support coordinates.
    pub const fn pencil(&self) -> &ShellModel {
        &self.pencil
    }
    /// Signed stress correction alone, with no positive-eigenvalue repair.
    pub const fn geometric_stiffness(&self) -> &Csr {
        &self.geometric
    }
    pub fn elastic_stiffness(&self) -> &Csr {
        &self.original.pencil().k
    }
    /// Recovered against this preparation's exact source/section/bubble.
    /// Only the solve entry point asserts equilibrium under supplied dead loads.
    pub const fn reference_state(&self) -> &ThermalShellResponse {
        &self.state
    }
    pub fn project_load(&self, full: &[f64]) -> Result<Vec<f64>, PlateError> {
        self.original.project_load(full)
    }
    pub fn expand(&self, reduced: &[f64]) -> Result<Vec<f64>, PlateError> {
        self.original.expand(reduced)
    }
}

impl<'a> PreparedThermalShell<'a> {
    /// Solve the original thermal equilibrium under explicit dead nodal loads,
    /// then construct its planar membrane initial-stress vibration/stability
    /// pencil. Supports, reference geometry and chosen mass remain unchanged.
    /// The initial-stress correction does not reapply the thermal load.
    pub fn solve_membrane_prestress(
        &self,
        external_load_6n: &[f64],
        limits: MembranePrestressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembranePrestressedShell<'a>, PlateError> {
        poll(&mut cancelled)?;
        admit_limits(limits)?;
        admit_planar(self.model, &mut cancelled)?;
        let state = self.solve_static(external_load_6n, &mut cancelled)?;
        self.assemble_membrane_prestress(state, limits, &mut cancelled)
    }

    /// Recompute actual thermal stress at a support-compatible physical state.
    /// No caller-authored stress array is accepted and no equilibrium claim is
    /// made for this read-only observer. A source consumer that already solved
    /// this exact thermal model can reuse its accepted displacement here.
    pub fn membrane_prestress_at(
        &self,
        displacement_6n: &[f64],
        limits: MembranePrestressLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembranePrestressedShell<'a>, PlateError> {
        poll(&mut cancelled)?;
        admit_limits(limits)?;
        admit_planar(self.model, &mut cancelled)?;
        admit_supports(self.model, displacement_6n)?;
        let state = self.recover(displacement_6n, &mut cancelled)?;
        self.assemble_membrane_prestress(state, limits, &mut cancelled)
    }

    fn assemble_membrane_prestress(
        &self,
        state: ThermalShellResponse,
        limits: MembranePrestressLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MembranePrestressedShell<'a>, PlateError> {
        let base = self.model.pencil();
        if base.free == 0 {
            return Err(invalid(
                "initial-stress pencil needs free incremental coordinates",
            ));
        }
        for nodal in state.displacement_6n.chunks_exact(6) {
            poll(cancelled)?;
            if norm(&nodal[3..]) > limits.maximum_rotation_rad {
                return Err(invalid(
                    "initial-stress state exceeds the physical rotation limit",
                ));
            }
        }
        let mut geometric = Coo::new(base.free, base.free);
        for (element, tri) in self.model.mesh.tris.iter().enumerate() {
            poll(cancelled)?;
            let g = self.model.mesh.facet(element)?;
            let nodal = core::array::from_fn(|i| state.displacement_6n[6 * tri[i / 6] + i % 6]);
            let gradient: [f64; 6] = core::array::from_fn(|i| {
                (0..3)
                    .map(|j| g.gradient[j][i / 3] * nodal[6 * j + i % 3])
                    .sum()
            });
            if norm(&gradient) > limits.maximum_displacement_gradient {
                return Err(invalid(
                    "initial-stress state exceeds the displacement-gradient limit",
                ));
            }
            if self.facets[element].director_motion_bound(&nodal, cancelled)?
                > limits.maximum_rotation_rad
            {
                return Err(invalid(
                    "initial-stress state exceeds the complete director limit",
                ));
            }
            let [xx, yy, xy] = state.facets[element].membrane_force_n_m;
            let mut scalar = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in i..3 {
                    let [ix, iy] = g.gradient[i];
                    let [jx, jy] = g.gradient[j];
                    let value =
                        g.area_m2 * (xx * ix * jx + xy * (ix * jy + iy * jx) + yy * iy * jy);
                    if !value.is_finite() {
                        return Err(invalid("unrepresentable membrane geometric stiffness"));
                    }
                    scalar[i][j] = value;
                    scalar[j][i] = value;
                }
            }
            // w_i = n dot u_i. Apply the ORIGINAL per-node support-frame
            // map; dof_map alone is insufficient for an oblique support.
            let mut entries = [(0, 0, 0.0); 9];
            let mut count = 0;
            for (local, &node) in tri.iter().enumerate() {
                for axis in 0..3 {
                    if let Some(row) = base.dof_map[6 * node + axis] {
                        let projection = dot(self.model.frames[node][axis], g.frame[2]);
                        if projection != 0.0 {
                            entries[count] = (local, row, projection);
                            count += 1;
                        }
                    }
                }
            }
            for (i, &(a, row, pa)) in entries[..count].iter().enumerate() {
                for &(b, column, pb) in &entries[i..count] {
                    let value = scalar[a][b] * (pa * pb);
                    if !value.is_finite() {
                        return Err(invalid(
                            "unrepresentable support-projected geometric stiffness",
                        ));
                    }
                    if value != 0.0 {
                        geometric.push(row, column, value);
                        if row != column {
                            geometric.push(column, row, value);
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        let geometric = geometric.assemble();
        let mut tangent = Coo::new(base.free, base.free);
        for row in 0..base.free {
            poll(cancelled)?;
            for matrix in [&base.k, &geometric] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    if !value.is_finite() {
                        return Err(invalid("nonfinite initial-stress stiffness"));
                    }
                    tangent.push(row, column, value);
                }
            }
        }
        let mut pencil = base.clone();
        pencil.k = tangent.assemble();
        for row in 0..pencil.free {
            poll(cancelled)?;
            if pencil.k.row(row).1.iter().any(|v| !v.is_finite()) {
                return Err(invalid("assembled initial-stress tangent overflows"));
            }
        }
        poll(cancelled)?;
        Ok(MembranePrestressedShell {
            original: self.model,
            pencil,
            geometric,
            state,
        })
    }
}

fn invalid(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PlateError> {
    if cancelled() {
        Err(PlateError::Cancelled)
    } else {
        Ok(())
    }
}
fn norm(values: &[f64]) -> f64 {
    if values.iter().any(|v| !v.is_finite()) {
        return f64::INFINITY;
    }
    let scale = values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 {
        0.0
    } else {
        scale * det::sqrt(values.iter().map(|v| (v / scale).powi(2)).sum())
    }
}
fn admit_limits(limits: MembranePrestressLimits) -> Result<(), PlateError> {
    if [
        limits.maximum_displacement_gradient,
        limits.maximum_rotation_rad,
    ]
    .iter()
    .any(|v| !v.is_finite() || *v <= 0.0 || *v > 0.1)
    {
        return Err(invalid(
            "initial-stress gradient and rotation limits must lie in (0,0.1]",
        ));
    }
    Ok(())
}
fn admit_planar(
    model: &MitcShellModel,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), PlateError> {
    let normal = model.mesh.facet(0)?.frame[2];
    let origin = model.mesh.nodes[0];
    let mut extent = 0.0_f64;
    let mut departure = 0.0_f64;
    for (point, director) in model.mesh.nodes.iter().zip(model.directors.vectors()) {
        poll(cancelled)?;
        let delta: [f64; 3] = core::array::from_fn(|i| point[i] - origin[i]);
        extent = extent.max(norm(&delta));
        departure = departure.max(dot(delta, normal).abs());
        if director
            .iter()
            .zip(normal)
            .any(|(a, b)| (a - b).abs() > 1e-12)
        {
            return Err(invalid(
                "membrane prestress requires constant facet-normal directors",
            ));
        }
    }
    if !extent.is_finite() || !departure.is_finite() || departure > 1e-12 * extent {
        return Err(invalid(
            "membrane prestress requires a planar reference panel",
        ));
    }
    for element in 0..model.mesh.tris.len() {
        poll(cancelled)?;
        if model.mesh.facet(element)?.frame[2]
            .iter()
            .zip(normal)
            .any(|(a, b)| (a - b).abs() > 1e-12)
        {
            return Err(invalid(
                "membrane prestress requires consistently oriented planar facets",
            ));
        }
    }
    Ok(())
}
fn admit_supports(model: &MitcShellModel, full: &[f64]) -> Result<(), PlateError> {
    let expanded = model.expand(&model.project_load(full)?)?;
    for (a, b) in full.chunks_exact(3).zip(expanded.chunks_exact(3)) {
        let scale = norm(a).max(norm(b));
        // Native support frames admit entrywise Gram error below 1e-12.
        // Its 3D operator norm can reach 3e-12 even for a compatible state.
        if a.iter().zip(b).any(|(a, b)| (a - b).abs() > 4e-12 * scale) {
            return Err(invalid(
                "initial-stress state violates the original homogeneous supports",
            ));
        }
    }
    Ok(())
}
