//! Frequency design about an original, independently equilibrated thermal case.
//!
//! The conservative initial-stress approximation uses the SAME reference mesh:
//! K_t = sum m_e K_e + G(sigma), with
//! G[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma grad(N_b).
//! Reference consistent mass is sum r_e rho_reference,e M_unit,e. The caller's
//! dimensionless mass fractions never replace the attached material density.
//!
//! For a mass-normalized mode, d lambda = phi^T(dK_t-lambda dM)phi.
//! The derivative of G includes sigma=m C(Bu-e_thermal) AND the equilibrated u.
//! A smooth modal aggregate needs just ONE original equilibrium adjoint,
//! K z = J_u; its additional modulus derivative is z^T(q_e-K_e u).
//! No global tangent or parameter Jacobian is assembled and no eigensolver is
//! introduced: fs-spectral owns the residual-qualified block inverse solve.
//!
//! This shares fs-solid::prestress's small-displacement initial-stress model.
//! It is not a finite-strain tangent, follower-load model, buckling certificate
//! or mode/state transfer. The admitted tangent must be positive definite.
//! Nonpositive diagonal/Ritz stiffness is refused, never clamped or absolutized;
//! finite search-space checks are not proof of global stability/completeness.

use super::super::{Diagonal, norm};
use super::{ComplianceError, ThermalComplianceGoal, WeightedOp, finite, gather};
use crate::SolveControl;
use crate::eigenfreq::smooth_min;
use crate::modal::{MatrixFreeEigenOptions, MatrixFreeEigenReport};
use fs_solver::LinearOp;
use fs_spectral::generalized::{
    GeneralizedError, GeneralizedOp, GeneralizedOptions, GeneralizedState,
};

/// Native physical, budget and spectral errors retain their original detail.
pub type ThermalFrequencyError = GeneralizedError<ComplianceError>;

#[derive(Debug, Clone, Copy)]
pub struct ThermalFrequencyOptions {
    pub solver: MatrixFreeEigenOptions,
    /// Positive smooth-min sharpness in inverse squared-frequency units.
    /// Keep an entire repeated cluster; truncating one is not differentiable.
    pub aggregation_beta: f64,
}
impl Default for ThermalFrequencyOptions {
    fn default() -> Self {
        Self {
            solver: MatrixFreeEigenOptions::default(),
            aggregation_beta: 1.0,
        }
    }
}

/// A frequency value and its two independent physical design derivatives.
/// The borrow binds every mode and derivative to the original thermal source.
pub struct ThermalFrequencyGoal<'g, 'a> {
    source: &'g ThermalComplianceGoal<'a>,
    case_index: usize,
    mass_fractions: Vec<f64>,
    modal: MatrixFreeEigenReport,
    value: f64,
    modulus_gradient: Vec<f64>,
    mass_fraction_gradient: Vec<f64>,
    adjoint_relative_residual: f64,
}
impl<'g, 'a> ThermalFrequencyGoal<'g, 'a> {
    /// Unnormalized smooth minimum of omega^2 [s^-2], possibly negative.
    /// This optimization aggregate is not itself a squared physical frequency.
    pub fn value(&self) -> f64 {
        self.value
    }
    pub fn modulus_gradient(&self) -> &[f64] {
        &self.modulus_gradient
    }
    pub fn mass_fraction_gradient(&self) -> &[f64] {
        &self.mass_fraction_gradient
    }
    pub fn modal(&self) -> &MatrixFreeEigenReport {
        &self.modal
    }
    pub fn source(&self) -> &'g ThermalComplianceGoal<'a> {
        self.source
    }
    pub fn case_index(&self) -> usize {
        self.case_index
    }
    pub fn mass_fractions(&self) -> &[f64] {
        &self.mass_fractions
    }
    pub fn adjoint_relative_residual(&self) -> f64 {
        self.adjoint_relative_residual
    }
}

struct PrestressedOp<'g, 'a> {
    source: &'g ThermalComplianceGoal<'a>,
    // One symmetric scalar 4x4 block per P1 cell; it acts identically on x/y/z.
    geometric: Vec<[f64; 16]>,
}
impl LinearOp for PrestressedOp<'_, '_> {
    fn n(&self) -> usize {
        self.source.base.elasticity.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let base = &self.source.base;
        let el = base.elasticity;
        WeightedOp {
            elasticity: el,
            weights: &base.moduli,
            identity: true,
        }
        .apply(x, y);
        for (tet, block) in el.tets.iter().zip(&self.geometric) {
            let local = gather(el, tet, x);
            for (a, &vertex) in tet.iter().enumerate() {
                for component in 0..3 {
                    let dof = 3 * vertex as usize + component;
                    if el.free[dof] {
                        for b in 0..4 {
                            y[dof] = block[4 * a + b].mul_add(local[3 * b + component], y[dof]);
                        }
                    }
                }
            }
        }
    }
}

struct Pencil<'g, 'a, 'c, 'callback> {
    tangent: PrestressedOp<'g, 'a>,
    mass: Vec<f64>,
    free: Vec<usize>,
    diagonal: Vec<f64>,
    tolerance: f64,
    control: &'c mut SolveControl<'callback>,
}
impl Pencil<'_, '_, '_, '_> {
    fn scatter(&self, reduced: &[f64]) -> Vec<f64> {
        let mut full = vec![0.0; self.tangent.n()];
        for (&dof, &value) in self.free.iter().zip(reduced) {
            full[dof] = value;
        }
        full
    }
    fn gather(&self, full: &[f64], reduced: &mut [f64]) {
        for (&dof, value) in self.free.iter().zip(reduced) {
            *value = full[dof];
        }
    }
}
impl GeneralizedOp for Pencil<'_, '_, '_, '_> {
    type Error = ComplianceError;
    fn dim(&self) -> usize {
        self.free.len()
    }
    fn apply_a(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.control.checkpoint("thermal-modal-tangent")?;
        let mut output = vec![0.0; self.tangent.n()];
        self.tangent.apply(&self.scatter(x), &mut output);
        self.gather(&output, y);
        self.control.checkpoint("thermal-modal-tangent")?;
        Ok(())
    }
    fn apply_b(&mut self, x: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        self.control.checkpoint("thermal-modal-mass")?;
        let mut output = vec![0.0; self.tangent.n()];
        self.tangent
            .source
            .base
            .elasticity
            .apply_mass(&self.mass, &self.scatter(x), &mut output);
        self.gather(&output, y);
        self.control.checkpoint("thermal-modal-mass")?;
        Ok(())
    }
    fn solve_a(&mut self, rhs: &[f64], y: &mut [f64]) -> Result<(), Self::Error> {
        let full = self.scatter(rhs);
        let solution = self.control.solve_preconditioned(
            &self.tangent,
            &Diagonal(&self.diagonal),
            &full,
            self.tolerance,
            self.tangent.source.base.options.max_iterations,
            "thermal-modal-inverse",
        )?;
        self.control.checkpoint("thermal-modal-inverse-residual")?;
        let mut residual = vec![0.0; full.len()];
        self.tangent.apply(&solution, &mut residual);
        for (r, &f) in residual.iter_mut().zip(&full) {
            *r -= f;
        }
        let relative = finite(
            norm(&residual) / norm(&full).max(f64::MIN_POSITIVE),
            "thermal modal inverse residual",
        )?;
        if relative >= self.tolerance {
            return Err(ComplianceError::Residual {
                stage: "thermal-modal-inverse",
                relative,
                tolerance: self.tolerance,
            });
        }
        self.gather(&solution, y);
        Ok(())
    }
    fn checkpoint(&mut self) -> Result<(), Self::Error> {
        self.control.checkpoint("thermal-modal")?;
        Ok(())
    }
}

impl<'a> ThermalComplianceGoal<'a> {
    /// Thermally prestressed frequency and exact modulus/mass-fraction gradients.
    /// Select ONE original equilibrium; case weights do not scale its prestress.
    /// Each fraction lies in [0,1] and multiplies the original reference density.
    /// No SIMP floor or void-mode mass interpolation is silently applied here.
    /// Temperature, source tensors, free strain, geometry and supports are fixed.
    ///
    /// Reuses the retained equilibrium. The mode solve uses the existing native
    /// work controller; aggregate differentiation adds one equilibrium adjoint.
    /// Every returned mode passes the original pencil-residual gate. One matrix
    /// application is not preemptible, as in the existing elasticity operator.
    pub fn prestressed_frequency_goal<'g>(
        &'g self,
        case_index: usize,
        mass_fractions: &[f64],
        options: ThermalFrequencyOptions,
        control: &mut SolveControl<'_>,
    ) -> Result<ThermalFrequencyGoal<'g, 'a>, ThermalFrequencyError> {
        let (mut pencil, block_size) = self
            .frequency_pencil(case_index, mass_fractions, options, control)
            .map_err(GeneralizedError::Operator)?;
        let settings = GeneralizedOptions {
            count: options.solver.count,
            block_size,
            relative_tolerance: options.solver.relative_tolerance,
            seed: options.solver.seed,
        };
        let mut state = GeneralizedState::new(&mut pencil, settings, None)?;
        let report = state.run(&mut pencil, options.solver.max_iterations)?;
        let modes = report.vectors.iter().map(|v| pencil.scatter(v)).collect();
        let mut modal = MatrixFreeEigenReport {
            values: report.values,
            modes,
            relative_residuals: report.relative_residuals,
            iterations: report.iterations,
            work: pencil.control.work(),
        };
        let (value, weights) = smooth_min(&modal.values, options.aggregation_beta);
        let (modulus_gradient, mass_fraction_gradient, adjoint_relative_residual) = self
            .frequency_gradient(case_index, &modal, &weights, pencil.control)
            .map_err(GeneralizedError::Operator)?;
        finite(value, "thermal frequency aggregate").map_err(GeneralizedError::Operator)?;
        pencil
            .control
            .checkpoint("thermal-frequency-publish")
            .map_err(|error| GeneralizedError::Operator(error.into()))?;
        modal.work = pencil.control.work();
        Ok(ThermalFrequencyGoal {
            source: self,
            case_index,
            mass_fractions: mass_fractions.to_vec(),
            modal,
            value,
            modulus_gradient,
            mass_fraction_gradient,
            adjoint_relative_residual,
        })
    }

    fn frequency_pencil<'g, 'c, 'callback>(
        &'g self,
        case_index: usize,
        mass_fractions: &[f64],
        options: ThermalFrequencyOptions,
        control: &'c mut SolveControl<'callback>,
    ) -> Result<(Pencil<'g, 'a, 'c, 'callback>, usize), ComplianceError> {
        control.checkpoint("thermal-frequency-admission")?;
        let el = self.base.elasticity;
        let case = self
            .cases
            .get(case_index)
            .ok_or(ComplianceError::InvalidInput(
                "thermal frequency case index out of range",
            ))?;
        if mass_fractions.len() != el.cells()
            || mass_fractions
                .iter()
                .any(|r| !r.is_finite() || !(0.0..=1.0).contains(r))
            || !options.aggregation_beta.is_finite()
            || options.aggregation_beta <= 0.0
            || !options.solver.linear_tolerance.is_finite()
            || options.solver.linear_tolerance <= 0.0
            || options.solver.linear_tolerance >= 1.0
        {
            return Err(ComplianceError::InvalidInput(
                "thermal modes require cell mass fractions in [0,1], positive beta and linear tolerance in (0,1)",
            ));
        }
        let free: Vec<usize> = el
            .free
            .iter()
            .enumerate()
            .filter_map(|(d, &f)| f.then_some(d))
            .collect();
        if options.solver.count == 0 || options.solver.count > free.len() {
            return Err(ComplianceError::InvalidInput(
                "thermal mode count must lie in 1..=free dofs",
            ));
        }
        let block_size = options
            .solver
            .count
            .checked_add(options.solver.oversampling)
            .ok_or(ComplianceError::InvalidInput(
                "thermal modal block width overflows",
            ))?
            .min(free.len());
        let mut mass = Vec::with_capacity(el.cells());
        let mut geometric = Vec::with_capacity(el.cells());
        let mut diagonal = self.base.diagonal.clone();
        for cell in 0..el.cells() {
            if cell % 256 == 0 {
                control.checkpoint("thermal-frequency-assembly")?;
            }
            mass.push(finite(
                mass_fractions[cell]
                    * el.material_density_kg_m3(cell)
                        .expect("original thermal material"),
                "thermal reference mass",
            )?);
            let s = &case.response.elements[cell].effective_stress_mandel_pa;
            let h = core::f64::consts::FRAC_1_SQRT_2;
            let sigma = [
                [s[0], h * s[3], h * s[5]],
                [h * s[3], s[1], h * s[4]],
                [h * s[5], h * s[4], s[2]],
            ];
            let mut block = [0.0; 16];
            for a in 0..4 {
                for b in a..4 {
                    let mut contraction = 0.0;
                    for i in 0..3 {
                        for j in 0..3 {
                            contraction +=
                                el.gradients[cell][a][i] * sigma[i][j] * el.gradients[cell][b][j];
                        }
                    }
                    let value = finite(
                        el.volumes[cell] * contraction,
                        "thermal geometric stiffness",
                    )?;
                    block[4 * a + b] = value;
                    block[4 * b + a] = value;
                }
                for component in 0..3 {
                    let dof = 3 * el.tets[cell][a] as usize + component;
                    if el.free[dof] {
                        diagonal[dof] += block[5 * a];
                    }
                }
            }
            geometric.push(block);
        }
        let mass_diagonal = el.mass_diagonal(&mass);
        if free.iter().any(|&d| {
            !diagonal[d].is_finite()
                || diagonal[d] <= 0.0
                || !mass_diagonal[d].is_finite()
                || mass_diagonal[d] <= 0.0
        }) {
            return Err(ComplianceError::InvalidInput(
                "every free thermal modal dof requires positive tangent stiffness and reference mass",
            ));
        }
        Ok((
            Pencil {
                tangent: PrestressedOp {
                    source: self,
                    geometric,
                },
                mass,
                free,
                diagonal,
                tolerance: options.solver.linear_tolerance,
                control,
            },
            block_size,
        ))
    }

    fn frequency_gradient(
        &self,
        case_index: usize,
        modal: &MatrixFreeEigenReport,
        weights: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(Vec<f64>, Vec<f64>, f64), ComplianceError> {
        let el = self.base.elasticity;
        let mut modulus = vec![0.0; el.cells()];
        let mut mass = vec![0.0; el.cells()];
        let mut rhs = vec![0.0; el.n];
        for ((mode, &lambda), &weight) in modal.modes.iter().zip(&modal.values).zip(weights) {
            for (cell, tet) in el.tets.iter().enumerate() {
                if cell % 256 == 0 {
                    control.checkpoint("thermal-frequency-gradient")?;
                }
                let local = gather(el, tet, mode);
                let mut elastic = 0.0;
                let mut kinetic = 0.0;
                for i in 0..12 {
                    for j in 0..12 {
                        elastic += local[i] * el.ke[cell][12 * i + j] * local[j];
                        kinetic += local[i] * el.me[cell][12 * i + j] * local[j];
                    }
                }
                // D = integral sum_component grad(phi_component) tensor grad(phi_component).
                // Mandel D is exactly the derivative of phi^T G(sigma) phi.
                let mut grad = [[0.0; 3]; 3];
                for a in 0..4 {
                    for i in 0..3 {
                        for j in 0..3 {
                            grad[i][j] += local[3 * a + i] * el.gradients[cell][a][j];
                        }
                    }
                }
                let gram = |i: usize, j: usize| -> f64 {
                    el.volumes[cell] * (0..3).map(|k| grad[k][i] * grad[k][j]).sum::<f64>()
                };
                let h = core::f64::consts::SQRT_2;
                let derivative = [
                    gram(0, 0),
                    gram(1, 1),
                    gram(2, 2),
                    h * gram(0, 1),
                    h * gram(1, 2),
                    h * gram(2, 0),
                ];
                let reference = self.reference_stress(case_index, cell);
                let explicit: f64 = derivative.iter().zip(reference).map(|(a, b)| a * b).sum();
                modulus[cell] = finite(
                    modulus[cell] + weight * (elastic + explicit),
                    "thermal frequency explicit modulus gradient",
                )?;
                mass[cell] = finite(
                    mass[cell]
                        - weight
                            * lambda
                            * kinetic
                            * el.material_density_kg_m3(cell)
                                .expect("original thermal material"),
                    "thermal frequency reference mass gradient",
                )?;
                self.scatter_stress(cell, &derivative, weight * self.base.moduli[cell], &mut rhs)?;
            }
        }
        let (adjoint, residual) =
            self.base
                .solve(&rhs, control, "thermal-frequency-equilibrium-adjoint")?;
        for (cell, tet) in el.tets.iter().enumerate() {
            if cell % 256 == 0 {
                control.checkpoint("thermal-frequency-adjoint-gradient")?;
            }
            let local = gather(el, tet, &adjoint);
            let balance = self.local_design_load(case_index, cell);
            modulus[cell] = finite(
                modulus[cell] + local.iter().zip(balance).map(|(a, b)| a * b).sum::<f64>(),
                "thermal frequency total modulus gradient",
            )?;
        }
        Ok((modulus, mass, residual))
    }
}
