//! Small perturbations about the same finitely deformed, loaded equilibrium.
//!
//! Native hyper3 owns equilibrium, the complete dP/dF tangent, supports and
//! consistent reference inertia. fs-modal owns the generalized eigensolve.
//! No additional initial-stress matrix is added to the finite-strain tangent.
//! All free coordinates and signed eigenvalues are retained; rigid and unstable
//! modes are not clipped or converted into real vibration frequencies.
//! Density may be cellwise constant or a positive P1 reference field. Both
//! routes retain the original consistent inertia through finite deformation;
//! averaging a graded density before constructing the pencil is not equivalent.
//!
//! This bounded dense route freezes hyperelastic material conditions, dead
//! loads and optional uniform closed-surface pressure. The latter contributes
//! the original conservative follower tangent exactly once. It does not admit
//! a viscous/plastic algorithmic tangent, nonconservative pressure patches,
//! evolving gravity, damping, finite-amplitude motion or mode truncation.
//! At constitutive switching surfaces it reports the original declared tangent,
//! not a unique two-sided harmonic response. P1 volumetric locking and spatial
//! convergence remain properties of the original discretization.

use super::pressure::ClosedTetPressureSurface;
use super::{HyperTetError, HyperTetProblem, HyperTetSolution, checkpoint};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_modal::{ModalError, ModePair, eigh_gen_dense};
use fs_sparse::Csr;

/// Sparse, inertia-certified signed windows on the original loaded pencil.
pub mod windows;

/// Hard bound on the existing non-interruptible dense eigensolver leaf.
pub const MAX_DENSE_FREE_DOFS: usize = 128;

/// Explicit allocation/work bound and eigenpair accuracy gate.
#[derive(Clone, Copy, Debug)]
pub struct VibrationSettings {
    /// Caller cap, additionally bounded by [`MAX_DENSE_FREE_DOFS`].
    pub maximum_free_dofs: usize,
    /// Maximum M-inverse residual divided by the largest absolute eigenvalue
    /// in the complete spectrum; finite and strictly between zero and one.
    /// A near-zero mode retains its signed value and residual interval.
    pub relative_residual_tolerance: f64,
}

/// Original equilibrium, spectral, input or cancellation refusal.
#[derive(Debug)]
pub enum VibrationError {
    Solid(HyperTetError),
    Modal(ModalError),
    /// An eigenpair failed the declared numerical gate, independently of sign.
    UnresolvedMode {
        index: usize,
        eigenvalue_per_s2: f64,
        residual_per_s2: f64,
    },
}

impl From<HyperTetError> for VibrationError {
    fn from(error: HyperTetError) -> Self {
        Self::Solid(error)
    }
}

impl core::fmt::Display for VibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-strain solid vibration: {self:?}")
    }
}

impl std::error::Error for VibrationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Modal(error) => Some(error),
            Self::UnresolvedMode { .. } => None,
        }
    }
}

/// One immutable equilibrium and its complete original free-coordinate pencil.
///
/// Density is per original reference volume, with its original spatial basis.
/// Geometry, material/source adapters, prescribed positions and dead loads remain
/// borrowed together. The modes describe infinitesimal perturbations with zero increments
/// on prescribed DOFs; their mass-normalized shapes use the retained free map.
#[derive(Debug)]
pub struct HyperTetVibration<'a, L: HyperelasticMaterial = Hyperelastic> {
    problem: HyperTetProblem<'a, L>,
    held_pressure: Option<(&'a ClosedTetPressureSurface, f64)>,
    reference_density: RetainedReferenceDensity,
    equilibrium: HyperTetSolution,
    free_dofs: Vec<usize>,
    stiffness_n_m: Vec<f64>,
    mass_kg: Vec<f64>,
    modes: Vec<ModePair>,
}

/// Select the existing mass owner without reinterpreting one spatial basis as
/// another. The result retains cell means only for inventory reporting.
#[derive(Clone, Copy)]
enum ReferenceDensityInput<'a> {
    PerElement(&'a [f64]),
    Nodal(&'a [f64]),
}

#[derive(Debug)]
struct RetainedReferenceDensity {
    cell_means_kg_m3: Vec<f64>,
    nodal_kg_m3: Option<Vec<f64>>,
}

impl ReferenceDensityInput<'_> {
    fn mass_matrix<L: HyperelasticMaterial>(
        self,
        problem: &HyperTetProblem<'_, L>,
        cx: &Cx<'_>,
    ) -> Result<Csr, HyperTetError> {
        match self {
            Self::PerElement(density) => problem.reference_mass_matrix(density, cx),
            Self::Nodal(density) => problem.reference_mass_matrix_with_nodal_density(density, cx),
        }
    }

    /// Retain an already admitted field. The means never enter the mass solve.
    fn retain(
        self,
        tetrahedra: &[[usize; 4]],
        cx: &Cx<'_>,
    ) -> Result<RetainedReferenceDensity, HyperTetError> {
        checkpoint(cx)?;
        let retained = match self {
            Self::PerElement(density) => RetainedReferenceDensity {
                cell_means_kg_m3: density.to_vec(),
                nodal_kg_m3: None,
            },
            Self::Nodal(density) => {
                let mut means = Vec::with_capacity(tetrahedra.len());
                for nodes in tetrahedra {
                    checkpoint(cx)?;
                    let rho = nodes.map(|node| density[node]);
                    // The P1 volume mean is the four-vertex mean. Normalize
                    // before adding, just as the native inertia owner does.
                    let mean = if rho.iter().all(|value| value.to_bits() == rho[0].to_bits()) {
                        rho[0]
                    } else {
                        let scale = rho.into_iter().fold(0.0_f64, f64::max);
                        scale * (0.25 * rho.iter().map(|value| value / scale).sum::<f64>())
                    };
                    means.push(mean);
                }
                RetainedReferenceDensity {
                    cell_means_kg_m3: means,
                    nodal_kg_m3: Some(density.to_vec()),
                }
            }
        };
        checkpoint(cx)?;
        Ok(retained)
    }
}

fn invalid(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput { what: what.into() }
}

fn dense_restriction(matrix: &Csr, free: &[usize], cx: &Cx<'_>) -> Result<Vec<f64>, HyperTetError> {
    let n = free.len();
    let mut dense = vec![0.0; n * n];
    for (row, &full_row) in free.iter().enumerate() {
        checkpoint(cx)?;
        for (col, &full_col) in free.iter().enumerate() {
            dense[row * n + col] = matrix.get(full_row, full_col);
        }
    }
    admit_symmetric(&dense, n)?;
    Ok(dense)
}

fn admit_symmetric(dense: &[f64], n: usize) -> Result<(), HyperTetError> {
    if dense.iter().any(|value| !value.is_finite()) {
        return Err(invalid("nonfinite finite-strain vibration pencil"));
    }
    let scale = dense.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let tolerance = (64.0 * f64::EPSILON * n as f64) * scale;
    // The shared dense owner symmetrizes roundoff after mass whitening.
    // Admit the original physical pencil before that operation.
    for row in 0..n {
        for col in 0..row {
            if (dense[row * n + col] - dense[col * n + row]).abs() > tolerance {
                return Err(invalid(
                    "native vibration pencil is not symmetric to roundoff",
                ));
            }
        }
    }
    Ok(())
}

impl<'a, L: HyperelasticMaterial> HyperTetProblem<'a, L> {
    /// Solve the original dead-load equilibrium, then its exact signed spectrum.
    ///
    /// Coupling the operations prevents a foreign/stale accepted displacement
    /// from silently becoming a prestressed pencil. The material seam is retained
    /// throughout: reference-deformation adapters are not replaced by bare laws.
    /// Supply exactly one positive reference density per original tetrahedron.
    ///
    /// # Errors
    /// Native geometry/material/Newton/mass errors, cancellation, a free-DOF
    /// budget violation, a nonsymmetric tangent or unresolved eigenpairs.
    pub fn solve_vibration(
        &self,
        reference_density_kg_m3: &[f64],
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        self.solve_vibration_impl(
            ReferenceDensityInput::PerElement(reference_density_kg_m3),
            None,
            None,
            settings,
            cx,
        )
    }

    /// The same equilibrium and complete spectrum with a positive P1 density
    /// field on the original mesh nodes. The native mass owner integrates
    /// `rho * N_a * N_b` exactly; equal cell means do not erase different
    /// material distributions. Density is admitted before nonlinear work and
    /// neither it nor the reference mass is rescaled by the accepted deformation.
    ///
    /// # Errors
    /// Nodal-density admission and all errors of [`Self::solve_vibration`].
    pub fn solve_vibration_with_nodal_density(
        &self,
        reference_nodal_density_kg_m3: &[f64],
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        self.solve_vibration_impl(
            ReferenceDensityInput::Nodal(reference_nodal_density_kg_m3),
            None,
            None,
            settings,
            cx,
        )
    }

    /// Original pressure-loaded equilibrium and its complete signed spectrum.
    ///
    /// Positive pressure acts inward on the same complete closed surface.
    /// The linearization is K=d(f_internal-f_pressure)/du, including p*d²V/du²
    /// exactly once. Pressure is held by an external reservoir; neither a finite
    /// gas stiffness nor pressure dynamics is inferred. The exact reference
    /// surface and pressure remain borrowed in the returned spectrum.
    ///
    /// # Errors
    /// The original pressure reference/finite-value gates and all errors of
    /// [`Self::solve_vibration`]. Nonuniform/open-patch pressure is not admitted.
    pub fn solve_vibration_with_pressure(
        &self,
        surface: &'a ClosedTetPressureSurface,
        pressure_pa: f64,
        reference_density_kg_m3: &[f64],
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        self.solve_vibration_impl(
            ReferenceDensityInput::PerElement(reference_density_kg_m3),
            None,
            Some((surface, pressure_pa)),
            settings,
            cx,
        )
    }

    /// Held closed-surface pressure with the original positive P1 density.
    /// The pressure tangent, material tangent and supports are identical to
    /// [`Self::solve_vibration_with_pressure`]; only the explicit spatial
    /// density basis changes. The complete original nodal field is retained.
    ///
    /// # Errors
    /// Nodal-density and native pressure/vibration refusals.
    pub fn solve_vibration_with_nodal_density_and_pressure(
        &self,
        surface: &'a ClosedTetPressureSurface,
        pressure_pa: f64,
        reference_nodal_density_kg_m3: &[f64],
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        self.solve_vibration_impl(
            ReferenceDensityInput::Nodal(reference_nodal_density_kg_m3),
            None,
            Some((surface, pressure_pa)),
            settings,
            cx,
        )
    }

    /// Re-equilibrate one absolute dead-load endpoint from an initial iterate,
    /// then return that equilibrium's complete signed spectrum.
    ///
    /// This delegates to [`Self::solve_from`], including its one-load-step rule.
    /// The initial field is never trusted as an accepted state. It permits a
    /// caller to select a finite branch without inventing a continuation solver;
    /// convergence to an unstable equilibrium retains negative eigenvalues.
    pub fn solve_vibration_from(
        &self,
        initial_displacement_m: &[f64],
        reference_density_kg_m3: &[f64],
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        self.solve_vibration_impl(
            ReferenceDensityInput::PerElement(reference_density_kg_m3),
            Some(initial_displacement_m),
            None,
            settings,
            cx,
        )
    }

    /// Re-equilibrate from an initial iterate with the same original positive
    /// P1 reference density. The initial field is not treated as an accepted
    /// state; [`Self::solve_vibration_from`]'s one-load-step rule still applies.
    ///
    /// # Errors
    /// Nodal-density and native initial-iterate/vibration refusals.
    pub fn solve_vibration_from_with_nodal_density(
        &self,
        initial_displacement_m: &[f64],
        reference_nodal_density_kg_m3: &[f64],
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        self.solve_vibration_impl(
            ReferenceDensityInput::Nodal(reference_nodal_density_kg_m3),
            Some(initial_displacement_m),
            None,
            settings,
            cx,
        )
    }

    fn solve_vibration_impl(
        &self,
        density: ReferenceDensityInput<'_>,
        initial: Option<&[f64]>,
        held_pressure: Option<(&'a ClosedTetPressureSurface, f64)>,
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        checkpoint(cx)?;
        self.validate_inputs()?;
        if let Some((surface, pressure)) = held_pressure {
            if !pressure.is_finite() || !surface.matches_reference(self.nodes_m, self.tetrahedra) {
                return Err(invalid(
                    "finite pressure and the exact reference pressure surface are required",
                )
                .into());
            }
        }
        if !settings.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&settings.relative_residual_tolerance)
            || settings.relative_residual_tolerance == 0.0
        {
            return Err(invalid("modal residual tolerance must be finite and in (0,1)").into());
        }
        let mut fixed = vec![false; 3 * self.nodes_m.len()];
        for &(dof, _) in self.prescribed_m {
            fixed[dof] = true;
        }
        let free_dofs: Vec<_> = (0..fixed.len()).filter(|&dof| !fixed[dof]).collect();
        let n = free_dofs.len();
        if n == 0 || n > settings.maximum_free_dofs || n > MAX_DENSE_FREE_DOFS {
            return Err(invalid(
                "finite-strain dense vibration requires 1..=128 bounded free DOFs",
            )
            .into());
        }
        // Admit mass before the nonlinear work. Neither density nor reference
        // volume is updated to the deformed geometry or statically lumped.
        let mass_kg = dense_restriction(&density.mass_matrix(self, cx)?, &free_dofs, cx)?;
        let equilibrium = match held_pressure {
            Some((surface, pressure)) => self.solve_with_pressure(surface, pressure, cx)?,
            None => match initial {
                None => self.solve(cx)?,
                Some(displacement) => self.solve_from(displacement, cx)?,
            },
        };
        let displacement: Vec<_> = equilibrium
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        let mut stiffness_n_m = dense_restriction(
            &self.prepare(cx)?.tangent_matrix(&displacement, cx)?,
            &free_dofs,
            cx,
        )?;
        if let Some((surface, pressure)) = held_pressure {
            if pressure != 0.0 {
                let external = dense_restriction(
                    &surface.force_tangent(&displacement, pressure, cx)?,
                    &free_dofs,
                    cx,
                )?;
                for (stiffness, load) in stiffness_n_m.iter_mut().zip(external) {
                    *stiffness -= load;
                }
                admit_symmetric(&stiffness_n_m, n)?;
            }
        }
        checkpoint(cx)?;
        let modes = eigh_gen_dense(&stiffness_n_m, &mass_kg, n).map_err(VibrationError::Modal)?;
        checkpoint(cx)?;
        if modes.len() != n {
            return Err(invalid("modal owner returned an incomplete finite-strain basis").into());
        }
        let scale = modes
            .iter()
            .map(|mode| mode.lambda.abs())
            .fold(0.0, f64::max);
        for (index, mode) in modes.iter().enumerate() {
            if mode.phi.len() != n
                || mode.phi.iter().any(|value| !value.is_finite())
                || !mode.lambda.is_finite()
                || !mode.residual.is_finite()
                || mode.residual < 0.0
                || !mode.interval.0.is_finite()
                || !mode.interval.1.is_finite()
                || mode.residual > settings.relative_residual_tolerance * scale
            {
                return Err(VibrationError::UnresolvedMode {
                    index,
                    eigenvalue_per_s2: mode.lambda,
                    residual_per_s2: mode.residual,
                });
            }
        }
        let reference_density = density.retain(self.tetrahedra, cx)?;
        Ok(HyperTetVibration {
            problem: *self,
            held_pressure,
            reference_density,
            equilibrium,
            free_dofs,
            stiffness_n_m,
            mass_kg,
            modes,
        })
    }
}

impl<L: HyperelasticMaterial> HyperTetVibration<'_, L> {
    /// Original mesh, material/source adapters, supports and dead loads.
    pub const fn problem(&self) -> &HyperTetProblem<'_, L> {
        &self.problem
    }

    /// Exact original closed pressure surface and held signed pressure [Pa].
    /// None selects dead loads only; Some with zero retains an explicitly
    /// supplied unloaded surface. Positive pressure acts inward.
    pub const fn held_pressure(&self) -> Option<(&ClosedTetPressureSurface, f64)> {
        self.held_pressure
    }

    /// Conserved mass density per original volume, in tetrahedron order.
    /// For a graded field these are cell means for inventory reporting; the
    /// mass matrix retains the full field's first and second moments.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density.cell_means_kg_m3
    }

    /// Original P1 nodal density, or None for the cellwise-density model.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.reference_density.nodal_kg_m3.as_deref()
    }

    /// Native accepted equilibrium, including finite support reactions.
    pub const fn equilibrium(&self) -> &HyperTetSolution {
        &self.equilibrium
    }

    /// Reduced-to-full node-major xyz DOFs; prescribed increments are zero.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Row-major free-by-free tangent, N/m: native material tangent minus any
    /// held-pressure force tangent, without a second initial-stress matrix.
    pub fn stiffness_n_m(&self) -> &[f64] {
        &self.stiffness_n_m
    }

    /// Row-major free-by-free consistent reference mass, kg.
    pub fn mass_kg(&self) -> &[f64] {
        &self.mass_kg
    }

    /// All signed squared frequencies, ascending, with mass-normalized shapes
    /// and original M-inverse residual intervals. A nonpositive eigenvalue is
    /// not a real phase frequency; an interval crossing zero does not resolve
    /// stability. No modal truncation or continuum stability claim is made.
    pub fn modes(&self) -> &[ModePair] {
        &self.modes
    }
}
