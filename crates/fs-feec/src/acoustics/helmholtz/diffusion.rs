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
}

impl HelmholtzFe<'_> {
    /// Solve the original heat-diffusion pencil with exact zero temperature
    /// perturbation on the supplied faces. No penalty, copied mesh, lumped mass,
    /// or fitted relaxation time is introduced. Every in-window mode is kept.
    /// Cancellation is checked during assembly and around the original bounded
    /// eigensolver call, which retains its existing cancellation limitations.
    pub fn isothermal_diffusion_modes(
        &self, wall_faces: &[[u32; 3]], maximum_eigenvalue_per_m2: f64,
        maximum_modes: usize, options: &SliceOptions,
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
        let report = slice_window(&stiffness.assemble(), &mass.assemble(),
            (0.0, maximum_eigenvalue_per_m2), options)?;
        if cancelled() { return Err(invalid("diffusion construction cancelled")); }
        if report.modes.len() > maximum_modes {
            return Err(AcousticFeError::TooManyModes { found: report.modes.len(), cap: maximum_modes });
        }
        if report.modes.is_empty() || report.below_low != 0
            || report.below_high != report.modes.len()
            || report.modes.iter().any(|mode| !mode.lambda.is_finite() || mode.lambda <= 0.0)
        {
            return Err(invalid("diffusion window must retain a nonempty complete positive spectrum"));
        }
        let mut shapes = Vec::with_capacity(report.modes.len());
        for mode in &report.modes {
            if cancelled() { return Err(invalid("diffusion construction cancelled")); }
            let mut shape = vec![0.0; self.ndof()];
            if mode.phi.len() != free.len() {
                return Err(invalid("diffusion eigensolver lost the constrained original basis"));
            }
            for (&global, &value) in free.iter().zip(&mode.phi) { shape[global] = value; }
            shapes.push(shape);
        }
        Ok(DiffusionModes {
            eigenvalues_per_m2: report.modes.iter().map(|mode| mode.lambda).collect(),
            shapes, residual: report.modes.iter().map(|mode| mode.residual).collect(),
            below_low: report.below_low, below_high: report.below_high,
            free_dofs: free.len(), wall_faces: wall.into_iter().collect(),
        })
    }
}

#[cfg(test)]
mod tests;
