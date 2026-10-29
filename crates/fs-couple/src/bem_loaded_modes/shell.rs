//! Source-owned shell geometry, tangent and full material viscosity.
use super::*;
use crate::plate_viscosity::ModalPlateViscosity;
use fs_plate::shell::reduction::ShellReduction;
use fs_plate::shell::reduction::radiation::RadiationSurfaceBudget;

impl BemLoadedModes {
    /// Use the supplied shell's selected basis and its own finite two-sided
    /// surface, including every edge wall. Mass remains the source identity;
    /// K is its exact potential Hessian at zero, including cross terms. The
    /// admitted full viscosity is copied unchanged. Its producer must use this
    /// same source basis; matching dimensions cannot establish its provenance.
    ///
    /// No prestressed/deformed geometry or residual flexibility outside the
    /// selected basis is inferred. Source-mode convergence remains required.
    /// Refinement holds each original panel's constant velocity fixed; this
    /// increases acoustic resolution, not structural/source-field resolution.
    /// Both surface budgets apply to FINAL refined panel/projection counts.
    pub fn from_shell(
        shell: &ShellReduction,
        nodal_thickness_m: &[f64],
        surface_budget: RadiationSurfaceBudget,
        medium: Medium,
        viscosity: Option<&ModalPlateViscosity>,
        settings: BemModesOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, BemModesError> {
        poll(&mut cancelled)?;
        let m = shell.mode_count();
        if m == 0
            || m >= settings.max_dimension
            || m >= MAX_NONLINEAR_POLE_DIMENSION
            || settings.panel_refinements >= usize::BITS as usize
            || viscosity.is_some_and(|c| c.mode_count() != m)
        {
            return Err(bad("shell basis, damping or refinement extent is invalid"));
        }
        let skin = shell
            .radiation_surface(nodal_thickness_m, surface_budget)
            .map_err(|e| BemModesError::Source {
                what: e.to_string(),
            })?;
        let p = skin
            .triangles()
            .len()
            .checked_mul(1usize << settings.panel_refinements)
            .ok_or_else(|| bad("refined panel count overflow"))?;
        if p > surface_budget.max_panels
            || p > settings.boundary.max_panels
            || p.checked_mul(m)
                .is_none_or(|v| v > surface_budget.max_panel_modes)
        {
            return Err(bad("refined acoustic skin exceeds its declared budget"));
        }
        Self::admit_size(m, p, &settings)?;
        let mut triangles = skin.triangles().to_vec();
        let mut weights = skin.normal_velocity_weights().to_vec();
        for _ in 0..settings.panel_refinements {
            poll(&mut cancelled)?;
            let mut next = Vec::with_capacity(2 * triangles.len());
            for &t in &triangles {
                let mut edge = 0;
                let mut longest = 0.0;
                for e in 0..3 {
                    let length = (0..3)
                        .map(|c| (t[e][c] - t[(e + 1) % 3][c]).abs())
                        .fold(0.0, f64::max);
                    let scaled_sq = if length == 0.0 {
                        0.0
                    } else {
                        (0..3)
                            .map(|c| ((t[e][c] - t[(e + 1) % 3][c]) / length).powi(2))
                            .sum::<f64>()
                    };
                    let distance = length * det::sqrt(scaled_sq);
                    if distance > longest {
                        longest = distance;
                        edge = e;
                    }
                }
                let (a, b, c) = (t[edge], t[(edge + 1) % 3], t[(edge + 2) % 3]);
                let midpoint = core::array::from_fn(|j| 0.5 * a[j] + 0.5 * b[j]);
                next.push([a, midpoint, c]);
                next.push([midpoint, b, c]);
            }
            triangles = next;
            for row in &mut weights {
                *row = row.iter().flat_map(|&x| [x, x]).collect();
            }
        }
        let surface =
            SpherePanels::from_triangles(triangles).map_err(|e| BemModesError::Source {
                what: e.to_string(),
            })?;
        let mut stiffness = vec![0.0; m * m];
        shell.hessian_matrix(&vec![0.0; m], &mut stiffness, &mut vec![[0.0; 3]; m]);
        poll(&mut cancelled)?;
        let damping = viscosity.map_or_else(|| vec![0.0; m * m], |c| c.matrix_per_s().to_vec());
        Self::build(
            surface,
            weights,
            stiffness,
            damping,
            medium,
            settings,
            &mut cancelled,
        )
    }
}
