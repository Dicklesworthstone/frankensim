//! Finite-bath Gibbs--Donnan equilibrium from the original convex virial potential.
//!
//! Ion exchange minimizes the sum of both compartments' original mixing free
//! energies, subject to every species total and fixed charge. The KKT force is
//! `ln(cL/cR) + B*(cL-cR) + z*psi`, with `psi=F*(phiL-phiR)/(RT)`.
//! This owner changes neither coefficients nor concentration support. The caps
//! limit applicability; they are never accepted as physical partition clamps.
//! Trapped osmotic particles are separate, ideal and noninteracting with B.
//! Fixed charge itself is not counted again as osmotic particles.
//! Each solve fixes both solvent volumes and the original model temperature.
//! The admitted virial model retains negligible solute partial volumes and its
//! single-phase concentration box; no association or thermal state is added.
//! This is fast-ion equilibrium, not a finite ionic permeability law.
//!
//! A prepared native dense LU workspace solves diagonally scaled exact Hessian
//! systems. All subsequent calls allocate nothing and restart deterministically;
//! failed or cancelled calls leave both caller output slices unchanged. Exact
//! absent totals and counterion-only limits need no artificial salt floor.

use std::sync::Arc;

use fs_la::LuWorkspace;

use super::{
    Accumulator, ChargePoint, DonnanControls, DonnanError, DonnanPotential, DonnanProblem,
    DonnanReport, FARADAY_C_PER_MOL, IdealDonnan, MAX_DONNAN_SPECIES, Result, admit, det, finite,
    poll, positive, solve_potential, split_amount,
};
use crate::osmosis::nonideal::{NonidealSolutionError, VirialSolutionModel, log_ratio};

#[cfg(test)]
mod tests;

/// Prepared fast-ion equilibrium using ONE immutable original virial model.
#[derive(Debug)]
pub struct VirialDonnan {
    ideal: IdealDonnan,
    model: Arc<VirialSolutionModel>,
    chemical_potential_tolerance_rt: f64,
    matrix: Box<[f64]>,
    factor: LuWorkspace,
}

impl Clone for VirialDonnan {
    fn clone(&self) -> Self {
        // Scratch has no physical history: a clone starts the same deterministic solve.
        Self {
            ideal: self.ideal.clone(),
            model: Arc::clone(&self.model),
            chemical_potential_tolerance_rt: self.chemical_potential_tolerance_rt,
            matrix: self.matrix.clone(),
            factor: LuWorkspace::new(self.valences().len())
                .expect("clone of an admitted bounded Donnan workspace"),
        }
    }
}

impl VirialDonnan {
    /// Every model species is an explicitly charged mobile ion, in original order.
    /// The absolute electrochemical residual tolerance is dimensionless (mu/RT),
    /// independently of the physical charge tolerance passed to each solve.
    pub fn new(
        valences: Vec<i32>,
        model: Arc<VirialSolutionModel>,
        chemical_potential_tolerance_rt: f64,
    ) -> Result<Self> {
        if model.species_count() != valences.len()
            || !positive(chemical_potential_tolerance_rt)
            || chemical_potential_tolerance_rt >= 1.0
        {
            return Err(DonnanError::Invalid(
                "virial Donnan needs the original ion layout and a chemical tolerance in (0,1)",
            ));
        }
        let ideal = IdealDonnan::new(valences)?;
        let n = ideal.valences().len();
        let factor = LuWorkspace::new(n)
            .map_err(|_| DonnanError::Invalid("virial Donnan workspace allocation failed"))?;
        let mut matrix = Vec::new();
        matrix
            .try_reserve_exact(n * n)
            .map_err(|_| DonnanError::Invalid("virial Donnan matrix allocation failed"))?;
        matrix.resize(n * n, 0.0);
        Ok(Self {
            ideal,
            model,
            chemical_potential_tolerance_rt,
            matrix: matrix.into_boxed_slice(),
            factor,
        })
    }
    pub fn model(&self) -> &VirialSolutionModel {
        &self.model
    }
    pub fn valences(&self) -> &[i32] {
        self.ideal.valences()
    }
    pub const fn chemical_potential_tolerance_rt(&self) -> f64 {
        self.chemical_potential_tolerance_rt
    }
    /// Conservative owned/shared heap footprint; counts the full shared model.
    pub fn heap_bytes(&self) -> usize {
        let n = self.valences().len();
        self.ideal.heap_bytes()
            + (2 * n * n + n) * core::mem::size_of::<f64>()
            + n * core::mem::size_of::<usize>()
            + self.model.heap_bytes()
            + core::mem::size_of::<VirialSolutionModel>()
            + 2 * core::mem::size_of::<usize>()
            + core::mem::align_of::<VirialSolutionModel>()
    }

    /// Solve with bounded native work and atomic output publication. `iterations`
    /// counts ideal initialization evaluations and all virial residual trials,
    /// including support-feasibility initialization and rejected damping trials.
    /// No state outside the original concentration box is published. A supported
    /// seed is found when the ideal partition lies outside that box.
    pub fn equilibrate_into(
        &mut self,
        problem: DonnanProblem<'_>,
        controls: DonnanControls,
        left_amounts_mol: &mut [f64],
        right_amounts_mol: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DonnanReport> {
        let n = self.valences().len();
        for actual in [
            problem.mobile_totals_mol.len(),
            left_amounts_mol.len(),
            right_amounts_mol.len(),
        ] {
            if actual != n {
                return Err(DonnanError::Shape {
                    expected: n,
                    actual,
                });
            }
        }
        let rt = admit(problem, controls)?;
        if problem.temperature_k != self.model.temperature_k() {
            return Err(DonnanError::Invalid(
                "virial Donnan temperature must match its original model",
            ));
        }
        let mut left = [0.0; MAX_DONNAN_SPECIES];
        let mut right = left;
        let mut report = self.ideal.equilibrate_into(
            problem,
            controls,
            &mut left[..n],
            &mut right[..n],
            &mut cancelled,
        )?;
        let mut work = report.iterations;
        let mut charge_scale = Accumulator::default();
        for charge in problem.fixed_charge_mol {
            charge_scale.add(charge.abs());
        }
        for (&z, &total) in self.valences().iter().zip(problem.mobile_totals_mol) {
            charge_scale.add((f64::from(z) * total).abs());
        }
        let charge_tolerance = finite(
            controls.absolute_charge_tolerance_mol
                + controls.relative_charge_tolerance * charge_scale.finish()?,
            "virial charge tolerance is not representable",
        )?;
        if let DonnanPotential::Finite {
            dimensionless_left_minus_right: mut psi,
            ..
        } = report.potential
        {
            if !self.supports(problem, &left[..n], &right[..n])? {
                self.supported_seed(
                    problem,
                    controls,
                    charge_tolerance,
                    &mut left[..n],
                    &mut right[..n],
                    &mut work,
                    &mut cancelled,
                )?;
            }
            let mut residual = [0.0; MAX_DONNAN_SPECIES];
            let mut weights = residual;
            let mut scaled_charge = residual;
            let mut chemical_solution = residual;
            let mut charge_solution = residual;
            loop {
                spend(&mut work, controls, &mut cancelled)?;
                let charge =
                    self.residual(problem, &left[..n], &right[..n], psi, &mut residual[..n])?;
                let current_merit = merit(
                    &residual[..n],
                    charge,
                    charge_tolerance,
                    self.chemical_potential_tolerance_rt,
                );
                if current_merit <= 1.0 {
                    break;
                }
                self.factor_hessian(problem, &left[..n], &right[..n], &mut weights[..n])?;
                let mut rhs = [0.0; MAX_DONNAN_SPECIES];
                for i in 0..n {
                    rhs[i] = weights[i] * residual[i];
                    scaled_charge[i] = weights[i] * f64::from(self.valences()[i]);
                }
                self.factor
                    .solve_factored_into(&rhs[..n], &mut chemical_solution[..n])
                    .map_err(|_| numerical())?;
                self.factor
                    .solve_factored_into(&scaled_charge[..n], &mut charge_solution[..n])
                    .map_err(|_| numerical())?;
                let denominator = dot(&scaled_charge[..n], &charge_solution[..n])?;
                if !positive(denominator) {
                    return Err(numerical());
                }
                let delta_psi = finite(
                    (charge - dot(&scaled_charge[..n], &chemical_solution[..n])?) / denominator,
                    "virial potential correction is not representable",
                )?;
                let mut direction = [0.0; MAX_DONNAN_SPECIES];
                let mut alpha: f64 = 1.0;
                for i in 0..n {
                    direction[i] = finite(
                        -weights[i] * (chemical_solution[i] + charge_solution[i] * delta_psi),
                        "virial amount correction is not representable",
                    )?;
                    if direction[i] > 0.0 {
                        alpha = alpha.min(0.99 * right[i] / direction[i]);
                    }
                    if direction[i] < 0.0 {
                        alpha = alpha.min(-0.99 * left[i] / direction[i]);
                    }
                }
                let mut trial_left = [0.0; MAX_DONNAN_SPECIES];
                let mut trial_right = trial_left;
                let mut trial_residual = trial_left;
                loop {
                    spend(&mut work, controls, &mut cancelled)?;
                    if !positive(alpha) {
                        return Err(DonnanError::NonConvergence);
                    }
                    for i in 0..n {
                        let change = alpha * direction[i];
                        let updated_left = left[i] + change;
                        let updated_right = right[i] - change;
                        // Preserve the minority explicitly before debiting the
                        // same amount from the majority; no trace is rounded away.
                        if updated_left <= updated_right {
                            trial_left[i] = updated_left;
                            trial_right[i] = problem.mobile_totals_mol[i] - updated_left;
                        } else {
                            trial_right[i] = updated_right;
                            trial_left[i] = problem.mobile_totals_mol[i] - updated_right;
                        }
                    }
                    let trial_psi = finite(
                        psi + alpha * delta_psi,
                        "virial trial potential is not representable",
                    )?;
                    if self.supports(problem, &trial_left[..n], &trial_right[..n])? {
                        let trial_charge = self.residual(
                            problem,
                            &trial_left[..n],
                            &trial_right[..n],
                            trial_psi,
                            &mut trial_residual[..n],
                        )?;
                        let trial_merit = merit(
                            &trial_residual[..n],
                            trial_charge,
                            charge_tolerance,
                            self.chemical_potential_tolerance_rt,
                        );
                        if trial_merit <= 1.0 || trial_merit < current_merit * (1.0 - 1e-4 * alpha)
                        {
                            left = trial_left;
                            right = trial_right;
                            psi = trial_psi;
                            break;
                        }
                    }
                    alpha *= 0.5;
                }
            }
            report.potential = DonnanPotential::Finite {
                dimensionless_left_minus_right: psi,
                voltage_left_minus_right_v: finite(
                    rt * psi / FARADAY_C_PER_MOL,
                    "virial Donnan voltage is not representable",
                )?,
            };
        }
        poll(&mut cancelled)?;
        let states = [
            self.model
                .evaluate_inventories(problem.volumes_m3[0], &left[..n])
                .map_err(material_error)?,
            self.model
                .evaluate_inventories(problem.volumes_m3[1], &right[..n])
                .map_err(material_error)?,
        ];
        for side in 0..2 {
            let amounts = if side == 0 { &left[..n] } else { &right[..n] };
            let mut osmoles = Accumulator::default();
            let mut charge = Accumulator::default();
            charge.add(problem.fixed_charge_mol[side]);
            for (&amount, &z) in amounts.iter().zip(self.valences()) {
                osmoles.add(amount);
                charge.add(f64::from(z) * amount);
            }
            report.mobile_osmoles_mol[side] = osmoles.finish()?;
            report.charge_residual_mol[side] = charge.finish()?;
            if report.charge_residual_mol[side].abs() > charge_tolerance {
                return Err(DonnanError::NonConvergence);
            }
            report.osmotic_pressure_pa[side] = finite(
                states[side].osmotic_pressure_pa
                    + rt * (problem.trapped_particles_mol[side] / problem.volumes_m3[side]),
                "virial Donnan osmotic pressure is not representable",
            )?;
        }
        report.iterations = work;
        report.osmotic_pressure_drop_slope_pa_per_m3 = self
            .pressure_slope(
                problem,
                &left[..n],
                &right[..n],
                states.map(|state| state.osmotic_stiffness_pa),
                rt,
                report.potential,
            )
            .ok();
        poll(&mut cancelled)?;
        left_amounts_mol.copy_from_slice(&left[..n]);
        right_amounts_mol.copy_from_slice(&right[..n]);
        Ok(report)
    }

    fn supports(&self, problem: DonnanProblem<'_>, left: &[f64], right: &[f64]) -> Result<bool> {
        for (side, amounts) in [left, right].into_iter().enumerate() {
            match self
                .model
                .evaluate_inventories(problem.volumes_m3[side], amounts)
            {
                Ok(_) => (),
                Err(NonidealSolutionError::OutsideDomain { .. }) => return Ok(false),
                Err(error) => return Err(material_error(error)),
            }
        }
        Ok(true)
    }

    /// Find an electrically feasible point INSIDE the original support box.
    /// This bounded logistic is only an initializer; the final KKT force uses
    /// the unclipped original chemical potentials, never this parametrization.
    #[allow(clippy::too_many_arguments)]
    fn supported_seed(
        &self,
        problem: DonnanProblem<'_>,
        controls: DonnanControls,
        tolerance: f64,
        left: &mut [f64],
        right: &mut [f64],
        work: &mut usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        let n = self.valences().len();
        let mut lower = [0.0; MAX_DONNAN_SPECIES];
        let mut upper = lower;
        let mut minimum = Accumulator::default();
        let mut maximum = Accumulator::default();
        minimum.add(problem.fixed_charge_mol[0]);
        maximum.add(problem.fixed_charge_mol[0]);
        for i in 0..n {
            let total = problem.mobile_totals_mol[i];
            let cap = self.model.maximum_concentrations_mol_m3()[i];
            lower[i] = (total - cap * problem.volumes_m3[1]).max(0.0);
            upper[i] = total.min(cap * problem.volumes_m3[0]);
            if lower[i] > upper[i] {
                return Err(domain());
            }
            let z = f64::from(self.valences()[i]);
            minimum.add(z * if z > 0.0 { lower[i] } else { upper[i] });
            maximum.add(z * if z > 0.0 { upper[i] } else { lower[i] });
        }
        let (minimum, maximum) = (minimum.finish()?, maximum.finish()?);
        if minimum > 0.0 || maximum < 0.0 {
            return Err(domain());
        }
        let log_volumes = log_ratio(problem.volumes_m3[0], problem.volumes_m3[1]);
        let mut evaluate = |parameter: f64| -> Result<ChargePoint> {
            spend(work, controls, cancelled)?;
            let mut charge = Accumulator::default();
            charge.add(problem.fixed_charge_mol[0]);
            let mut slope = Accumulator::default();
            for i in 0..n {
                let z = f64::from(self.valences()[i]);
                let span = upper[i] - lower[i];
                let (a, b) = split_amount(span, log_volumes - z * parameter);
                charge.add(z * (lower[i] + a));
                if span > 0.0 {
                    slope.add(z * z * a.min(b) * (a.max(b) / span));
                }
            }
            Ok(ChargePoint {
                psi: parameter,
                residual: charge.finish()?,
                negative_derivative: slope.finish()?,
            })
        };
        let parameter = if minimum == 0.0 {
            f64::INFINITY
        } else if maximum == 0.0 {
            f64::NEG_INFINITY
        } else {
            let scale = self
                .valences()
                .iter()
                .map(|z| f64::from(*z).abs())
                .fold(1.0, f64::max);
            solve_potential(0.25 * tolerance, 1.0 / scale, &mut evaluate)?
        };
        for i in 0..n {
            let (a, b) = split_amount(
                upper[i] - lower[i],
                log_volumes - f64::from(self.valences()[i]) * parameter,
            );
            left[i] = lower[i] + a;
            right[i] = (problem.mobile_totals_mol[i] - upper[i]) + b;
        }
        if !self.supports(problem, left, right)? {
            return Err(domain());
        }
        Ok(())
    }

    fn residual(
        &self,
        problem: DonnanProblem<'_>,
        left: &[f64],
        right: &[f64],
        psi: f64,
        residual: &mut [f64],
    ) -> Result<f64> {
        let n = self.valences().len();
        let [vl, vr] = problem.volumes_m3;
        let mut charge = Accumulator::default();
        charge.add(problem.fixed_charge_mol[0]);
        for i in 0..n {
            charge.add(f64::from(self.valences()[i]) * left[i]);
            if problem.mobile_totals_mol[i] == 0.0 {
                residual[i] = 0.0;
                continue;
            }
            if !positive(left[i]) || !positive(right[i]) {
                return Err(DonnanError::Arithmetic(
                    "finite virial equilibrium needs representable minority ions",
                ));
            }
            let mut value = Accumulator::default();
            value.add(log_ratio(left[i], right[i]) + log_ratio(vr, vl));
            for j in 0..n {
                value.add(
                    self.model.interaction_matrix_m3_mol()[i * n + j]
                        * (left[j] / vl - right[j] / vr),
                );
            }
            value.add(f64::from(self.valences()[i]) * psi);
            residual[i] = value.finish()?;
        }
        charge.finish()
    }

    /// K = D H D, H=diag(1/nL+1/nR)+B*(1/VL+1/VR).
    /// D_i=sqrt(nL*nR/N) admits traces without forming a huge reciprocal.
    fn factor_hessian(
        &mut self,
        problem: DonnanProblem<'_>,
        left: &[f64],
        right: &[f64],
        weights: &mut [f64],
    ) -> Result<()> {
        let n = self.valences().len();
        let inverse_volume = finite(
            1.0 / problem.volumes_m3[0] + 1.0 / problem.volumes_m3[1],
            "virial inverse volume is not representable",
        )?;
        for i in 0..n {
            weights[i] = if problem.mobile_totals_mol[i] == 0.0 {
                0.0
            } else {
                det::sqrt(left[i].min(right[i]))
                    * det::sqrt(left[i].max(right[i]) / problem.mobile_totals_mol[i])
            };
        }
        for i in 0..n {
            for j in 0..n {
                self.matrix[i * n + j] = finite(
                    (if i == j { 1.0 } else { 0.0 })
                        + ((self.model.interaction_matrix_m3_mol()[i * n + j] * weights[i])
                            * inverse_volume)
                            * weights[j],
                    "virial scaled Hessian is not representable",
                )?;
            }
        }
        self.factor.factor(&self.matrix).map_err(|_| numerical())
    }

    fn pressure_slope(
        &mut self,
        problem: DonnanProblem<'_>,
        left: &[f64],
        right: &[f64],
        stiffness: [f64; 2],
        rt: f64,
        potential: DonnanPotential,
    ) -> Result<f64> {
        let n = self.valences().len();
        let [vl, vr] = problem.volumes_m3;
        if !matches!(potential, DonnanPotential::Finite { .. }) {
            let mut curvature = Accumulator::default();
            for side in 0..2 {
                let volume = problem.volumes_m3[side];
                curvature.add(stiffness[side] / volume);
                curvature.add((rt * (problem.trapped_particles_mol[side] / volume)) / volume);
            }
            return curvature.finish();
        }
        let mut weights = [0.0; MAX_DONNAN_SPECIES];
        self.factor_hessian(problem, left, right, &mut weights[..n])?;
        let mut h = weights;
        let mut z = weights;
        for i in 0..n {
            let mut value = Accumulator::default();
            value.add(1.0 / vl + 1.0 / vr);
            for j in 0..n {
                value.add(
                    self.model.interaction_matrix_m3_mol()[i * n + j]
                        * ((left[j] / vl) / vl + (right[j] / vr) / vr),
                );
            }
            h[i] = finite(
                weights[i] * value.finish()?,
                "virial water sensitivity is not representable",
            )?;
            z[i] = weights[i] * f64::from(self.valences()[i]);
        }
        let mut u = [0.0; MAX_DONNAN_SPECIES];
        let mut w = u;
        self.factor
            .solve_factored_into(&h[..n], &mut u[..n])
            .map_err(|_| numerical())?;
        self.factor
            .solve_factored_into(&z[..n], &mut w[..n])
            .map_err(|_| numerical())?;
        let denominator = dot(&z[..n], &w[..n])?;
        if !positive(denominator) {
            return Err(numerical());
        }
        let projection = dot(&z[..n], &u[..n])?;
        let multiplier = finite(
            projection / denominator,
            "virial constrained water response is not representable",
        )?;
        // nL' = D*(-u + w*(z.D.u)/(z.D.w)). For x water transferred L->R,
        // the reduced curvature/RT is the sum over both sides of
        // (nL'+c)^T [diag(1/c)+B] (nL'+c)/V, plus ideal trapped n/V^2.
        // This avoids subtracting two large Schur-complement terms when the
        // relaxed osmotic stiffness is small (or exactly zero).
        let mut curvature = Accumulator::default();
        for (side, amounts) in [left, right].into_iter().enumerate() {
            let volume = problem.volumes_m3[side];
            let mut response = [0.0; MAX_DONNAN_SPECIES];
            let mut quadratic = Accumulator::default();
            for i in 0..n {
                if problem.mobile_totals_mol[i] == 0.0 {
                    continue;
                }
                response[i] = finite(
                    weights[i].mul_add(multiplier.mul_add(w[i], -u[i]), amounts[i] / volume),
                    "virial concentration response is not representable",
                )?;
                // v_i^2/(V*c_i) = (v_i/sqrt(n_i))^2 avoids forming a
                // singular reciprocal for a representable trace inventory.
                let ideal = response[i] / det::sqrt(amounts[i]);
                quadratic.add(ideal * ideal);
            }
            for i in 0..n {
                for j in 0..n {
                    quadratic.add(
                        ((self.model.interaction_matrix_m3_mol()[i * n + j] * response[i])
                            / volume)
                            * response[j],
                    );
                }
            }
            let value = quadratic.finish()?;
            // Convexity certifies the sign; a numerically unresolved quadratic
            // still refuses instead of clamping it to an invented stiffness.
            if value < 0.0 {
                return Err(numerical());
            }
            curvature.add(value);
            curvature.add((problem.trapped_particles_mol[side] / volume) / volume);
        }
        finite(
            rt * curvature.finish()?,
            "virial osmotic slope is not representable",
        )
    }
}

fn merit(residual: &[f64], charge: f64, charge_tolerance: f64, chemical_tolerance: f64) -> f64 {
    residual.iter().map(|r| r.abs() / chemical_tolerance).fold(
        if charge == 0.0 {
            0.0
        } else {
            charge.abs() / charge_tolerance
        },
        f64::max,
    )
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    let mut sum = Accumulator::default();
    for (&a, &b) in a.iter().zip(b) {
        sum.add(a * b);
    }
    sum.finish()
}
fn spend(
    work: &mut usize,
    controls: DonnanControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    if *work >= controls.maximum_iterations {
        return Err(DonnanError::NonConvergence);
    }
    *work += 1;
    Ok(())
}
fn numerical() -> DonnanError {
    DonnanError::Arithmetic("virial Donnan linear system or derivative is unresolved")
}
fn domain() -> DonnanError {
    DonnanError::Invalid("virial Donnan equilibrium is outside the original concentration support")
}
fn material_error(error: NonidealSolutionError) -> DonnanError {
    match error {
        NonidealSolutionError::OutsideDomain { .. } => domain(),
        _ => DonnanError::Arithmetic("original virial thermodynamics refused the Donnan state"),
    }
}
