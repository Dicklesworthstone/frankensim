//! Small perturbations about the same finitely deformed, loaded equilibrium.
//!
//! Native hyper3 owns equilibrium, the complete dP/dF tangent, supports and
//! consistent reference inertia. fs-modal owns the generalized eigensolve.
//! No additional initial-stress matrix is added to the finite-strain tangent.
//! All free coordinates and signed eigenvalues are retained; rigid and unstable
//! modes are not clipped or converted into real vibration frequencies.
//!
//! This bounded dense route freezes hyperelastic material conditions and dead
//! loads. It does not admit a viscous/plastic algorithmic tangent, follower
//! pressure, evolving gravity, damping, finite-amplitude motion or mode truncation.
//! At constitutive switching surfaces it reports the original declared tangent,
//! not a unique two-sided harmonic response. P1 volumetric locking and spatial
//! convergence remain properties of the original discretization.

use super::{HyperTetError, HyperTetProblem, HyperTetSolution, checkpoint};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_modal::{ModalError, ModePair, eigh_gen_dense};
use fs_sparse::Csr;

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
/// Density is per original reference volume, in element order. The geometry,
/// material/source adapters, prescribed positions and dead loads remain borrowed
/// together. The modes describe infinitesimal perturbations with zero increments
/// on prescribed DOFs; their mass-normalized shapes use the retained free map.
#[derive(Debug)]
pub struct HyperTetVibration<'a, L: HyperelasticMaterial = Hyperelastic> {
    problem: HyperTetProblem<'a, L>,
    reference_density_kg_m3: Vec<f64>,
    equilibrium: HyperTetSolution,
    free_dofs: Vec<usize>,
    stiffness_n_m: Vec<f64>,
    mass_kg: Vec<f64>,
    modes: Vec<ModePair>,
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
    Ok(dense)
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
        self.solve_vibration_impl(reference_density_kg_m3, None, settings, cx)
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
            reference_density_kg_m3,
            Some(initial_displacement_m),
            settings,
            cx,
        )
    }

    fn solve_vibration_impl(
        &self,
        reference_density_kg_m3: &[f64],
        initial: Option<&[f64]>,
        settings: VibrationSettings,
        cx: &Cx<'_>,
    ) -> Result<HyperTetVibration<'a, L>, VibrationError> {
        checkpoint(cx)?;
        self.validate_inputs()?;
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
        let mass_kg = dense_restriction(
            &self.reference_mass_matrix(reference_density_kg_m3, cx)?,
            &free_dofs,
            cx,
        )?;
        let equilibrium = match initial {
            None => self.solve(cx)?,
            Some(displacement) => self.solve_from(displacement, cx)?,
        };
        let displacement: Vec<_> = equilibrium
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        let stiffness_n_m = dense_restriction(
            &self.prepare(cx)?.tangent_matrix(&displacement, cx)?,
            &free_dofs,
            cx,
        )?;
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
        checkpoint(cx)?;
        Ok(HyperTetVibration {
            problem: *self,
            reference_density_kg_m3: reference_density_kg_m3.to_vec(),
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

    /// Conserved mass density per original volume, in tetrahedron order.
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }

    /// Native accepted equilibrium, including finite support reactions.
    pub const fn equilibrium(&self) -> &HyperTetSolution {
        &self.equilibrium
    }

    /// Reduced-to-full node-major xyz DOFs; prescribed increments are zero.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Row-major free-by-free native tangent, N/m, without an extra stress term.
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
