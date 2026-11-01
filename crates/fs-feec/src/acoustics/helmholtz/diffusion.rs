//! Isothermal-wall diffusion on the original acoustic FE geometry and basis.
//!
//! Solve `K theta = lambda M theta` for unit diffusivity, imposing homogeneous
//! Dirichlet data only on explicitly selected wall faces. Other boundary faces
//! have natural zero heat flux. The material owner supplies diffusivity.

use super::*;
use std::collections::BTreeSet;

/// Original-coordinate thermal modes and the unchanged eigenvalue certificate.
#[derive(Clone, Debug)]
pub struct DiffusionModes {
    /// Unit-diffusivity Laplacian eigenvalues [1/m²].
    pub eigenvalues_per_m2: Vec<f64>,
    /// M-normalized original FE vectors; constrained entries are exactly zero.
    pub shapes: Vec<Vec<f64>>,
    /// Original certified generalized-eigenpair residuals.
    pub residual: Vec<f64>,
    /// Original inertia counts for the declared eigenwindow.
    pub below_low: usize,
    pub below_high: usize,
    /// Number of unconstrained original FE degrees of freedom.
    pub free_dofs: usize,
    /// Canonical, sorted, unique isothermal wall faces.
    pub wall_faces: Vec<[u32; 3]>,
    /// True when the whole free pencil was solved once, densely, and the window only filtered the result: the
    /// eigenvalues and shapes are then bit-identical for every window that keeps them. False for a pencil above
    /// [`DENSE_FULL_PENCIL_DOFS`], which the sparse window owner solves; its shift sits at the window's midpoint,
    /// so its last bits follow the window.
    pub dense: bool,
}

/// Largest free-DOF count solved as ONE dense full pencil (a few n x n matrices; cyclic Jacobi). The unit-diffusivity
/// eigenvalues are a property of the geometry alone, and the gas only decides which of them are kept
/// (`max_decay_rate / diffusivity`, so the ceiling moves with pressure): solving the same pencil whatever the ceiling
/// makes them independent of it, bitwise. Above this size a dense solve costs more than the sparse window owner, which
/// is then used as before.
pub const DENSE_FULL_PENCIL_DOFS: usize = 512;

impl HelmholtzFe<'_> {
    /// Solve the original heat-diffusion pencil with exact zero temperature
    /// perturbation on the supplied faces. No penalty, copied mesh, lumped mass,
    /// or fitted relaxation time is introduced. Every in-window mode is kept.
    /// Cancellation is checked during assembly and around the original bounded
    /// eigensolver call, which retains its existing cancellation limitations.
    pub fn isothermal_diffusion_modes(
        &self, wall_faces: &[[u32; 3]], maximum_eigenvalue_per_m2: f64,
        maximum_modes: usize, options: &SliceOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DiffusionModes, AcousticFeError> {
        self.isothermal_diffusion_modes_with(wall_faces, maximum_eigenvalue_per_m2, maximum_modes, options,
            DENSE_FULL_PENCIL_DOFS, cancelled)
    }

    /// [`Self::isothermal_diffusion_modes`] with an explicit dense-route ceiling (0 forces the sparse window owner).
    fn isothermal_diffusion_modes_with(
        &self, wall_faces: &[[u32; 3]], maximum_eigenvalue_per_m2: f64,
        maximum_modes: usize, options: &SliceOptions, dense_limit: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DiffusionModes, AcousticFeError> {
        let invalid = |what| AcousticFeError::InvalidInput { what };
        if wall_faces.is_empty() || !maximum_eigenvalue_per_m2.is_finite()
            || maximum_eigenvalue_per_m2 <= 0.0 || !(1..=MAX_MODES).contains(&maximum_modes)
        {
            return Err(invalid("diffusion needs isothermal faces and a bounded positive eigenwindow"));
        }
        let quad = tri_quadrature(self.order + 1);
        let mut wall = BTreeSet::new();
        let mut constrained = vec![false; self.ndof()];
        for &face in wall_faces {
            if cancelled() { return Err(invalid("diffusion construction cancelled")); }
            let mut key = face; key.sort_unstable();
            if !wall.insert(key) { return Err(invalid("diffusion wall faces must be unique")); }
            // Original hierarchical basis functions with nonzero face trace.
            for dof in self.face_rule(face, &quad)?.dofs { constrained[dof] = true; }
        }
        let free: Vec<_> = constrained.iter().enumerate()
            .filter_map(|(index, &fixed)| (!fixed).then_some(index)).collect();
        if free.is_empty() {
            return Err(invalid("isothermal thermal field needs interior FE degrees of freedom"));
        }
        let mut index = vec![usize::MAX; self.ndof()];
        for (local, &global) in free.iter().enumerate() { index[global] = local; }
        let mut stiffness = Coo::new(free.len(), free.len());
        let mut mass = Coo::new(free.len(), free.len());
        for (row, &global) in free.iter().enumerate() {
            if cancelled() { return Err(invalid("diffusion construction cancelled")); }
            for (source, target) in [(&self.stiffness, &mut stiffness), (&self.mass, &mut mass)] {
                let (columns, values) = source.row(global);
                for (&column, &value) in columns.iter().zip(values) {
                    if index[column] != usize::MAX { target.push(row, index[column], value); }
                }
            }
        }
        let (stiffness, mass) = (stiffness.assemble(), mass.assemble());
        let dense = free.len() <= dense_limit;
        let (modes, below_low, below_high) = if dense {
            // The same dense pencil whatever the ceiling; the ceiling only filters, in the sparse owner's own
            // window convention (low, high].
            let n = free.len();
            let (mut k, mut m) = (vec![0.0; n * n], vec![0.0; n * n]);
            for (source, target) in [(&stiffness, &mut k), (&mass, &mut m)] {
                for row in 0..n {
                    let (columns, values) = source.row(row);
                    for (&column, &value) in columns.iter().zip(values) { target[row * n + column] += value; }
                }
            }
            if cancelled() { return Err(invalid("diffusion construction cancelled")); }
            let all = fs_modal::eigh_gen_dense(&k, &m, n)?;
            let below_low = all.iter().filter(|mode| mode.lambda <= 0.0).count();
            let kept: Vec<_> = all.into_iter()
                .filter(|mode| mode.lambda > 0.0 && mode.lambda <= maximum_eigenvalue_per_m2).collect();
            let below_high = below_low + kept.len();
            (kept.into_iter().map(|mode| (mode.lambda, mode.phi, mode.residual)).collect::<Vec<_>>(),
                below_low, below_high)
        } else {
            let report = slice_window(&stiffness, &mass, (0.0, maximum_eigenvalue_per_m2), options)?;
            (report.modes.into_iter().map(|mode| (mode.lambda, mode.phi, mode.residual)).collect::<Vec<_>>(),
                report.below_low, report.below_high)
        };
        if cancelled() { return Err(invalid("diffusion construction cancelled")); }
        if modes.len() > maximum_modes {
            return Err(AcousticFeError::TooManyModes { found: modes.len(), cap: maximum_modes });
        }
        if modes.is_empty() || below_low != 0 || below_high != modes.len()
            || modes.iter().any(|(lambda, _, _)| !lambda.is_finite() || *lambda <= 0.0)
        {
            return Err(invalid("diffusion window must retain a nonempty complete positive spectrum"));
        }
        let mut shapes = Vec::with_capacity(modes.len());
        for (_, phi, _) in &modes {
            if cancelled() { return Err(invalid("diffusion construction cancelled")); }
            let mut shape = vec![0.0; self.ndof()];
            if phi.len() != free.len() {
                return Err(invalid("diffusion eigensolver lost the constrained original basis"));
            }
            for (&global, &value) in free.iter().zip(phi) { shape[global] = value; }
            shapes.push(shape);
        }
        Ok(DiffusionModes {
            eigenvalues_per_m2: modes.iter().map(|(lambda, _, _)| *lambda).collect(),
            shapes, residual: modes.iter().map(|(_, _, residual)| *residual).collect(),
            below_low, below_high,
            free_dofs: free.len(), wall_faces: wall.into_iter().collect(), dense,
        })
    }
}

#[cfg(test)]
mod tests;
