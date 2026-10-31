//! Sparse real-equivalent field solves; the shared solver owns every recurrence.

use super::*;
use crate::planar::eddy::moving::MovingJacobi;
use fs_solver::{
    FgmresState, FlexiblePreconditioner, LinearOp, RealEquivalentComplexOp, ResidualClaim,
    SolverRunProgress,
};
use std::cell::{Cell, RefCell};

struct ReducedMagnetic<'a> {
    tangent: &'a Csr,
    free: &'a [usize],
    full: RefCell<Vec<f64>>,
    check: &'a dyn Fn() -> bool,
}
impl LinearOp for ReducedMagnetic<'_> {
    fn n(&self) -> usize {
        self.free.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut full = self.full.borrow_mut();
        for (&v, &value) in self.free.iter().zip(x) {
            full[v] = value;
        }
        for (row, &v) in self.free.iter().enumerate() {
            if (self.check)() {
                y.fill(f64::NAN);
                return;
            }
            let (columns, values) = self.tangent.row(v);
            y[row] = columns.iter().zip(values).map(|(&w, &a)| a * full[w]).sum();
        }
    }
    // The admitted constitutive Hessian and homogeneous restriction are symmetric.
}

struct ReducedConduction<'a> {
    conduction: &'a Conduction,
    free: &'a [usize],
    omega: f64,
    scratch: RefCell<(Vec<f64>, Vec<f64>)>,
    check: &'a dyn Fn() -> bool,
}
impl LinearOp for ReducedConduction<'_> {
    fn n(&self) -> usize {
        self.free.len()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let mut scratch = self.scratch.borrow_mut();
        let (full, applied) = &mut *scratch;
        for (&v, &value) in self.free.iter().zip(x) {
            full[v] = value;
        }
        if !self.conduction.apply(full, applied, self.check) {
            y.fill(f64::NAN);
            return;
        }
        for (&v, out) in self.free.iter().zip(y) {
            *out = -self.omega * applied[v];
        }
    }
    // -omega times a symmetric PSD projected mass. The existing real-equivalent
    // wrapper supplies the nonsymmetric block's CORRECT signed transpose.
}

struct HarmonicJacobi(MovingJacobi);
impl FlexiblePreconditioner for HarmonicJacobi {
    fn apply(&self, iteration: usize, residual: &[f64], output: &mut [f64]) {
        let n = residual.len() / 2;
        let (rr, ri) = residual.split_at(n);
        let (zr, zi) = output.split_at_mut(n);
        self.0.apply(iteration, rr, zr);
        self.0.apply(iteration, ri, zi);
    }
}

impl PlanarHarmonicWindings<'_, '_> {
    fn solve_units(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<HarmonicWindingUnitResponse>> {
        // With no conductor the DC derivative fields solve the same equations.
        // Reuse only after recomputing their residual against THIS caller's
        // allowances; a looser earlier DC solve must not bypass AC accuracy.
        if self.conduction.cells.is_empty() {
            let mut responses = Vec::with_capacity(self.basis.winding_count());
            let mut applied = vec![0.0; self.tangent.nrows()];
            for (response, load) in self.basis.unit_responses().iter().zip(&self.loads) {
                poll(cancelled)?;
                self.tangent
                    .spmv(response.potential_per_current_h_m(), &mut applied);
                let mut norm = 0.0_f64;
                let mut residual = 0.0_f64;
                for &v in &self.free {
                    poll(cancelled)?;
                    norm = norm.hypot(load[v]);
                    residual = residual.hypot(finite(load[v] - applied[v])?);
                }
                finite(norm)?;
                let tolerance = if norm == 0.0 {
                    self.controls.relative_tolerance
                } else {
                    self.controls
                        .relative_tolerance
                        .max((self.controls.absolute_residual_a_m / norm).min(f64::MAX))
                };
                let relative = if norm == 0.0 {
                    if residual == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    residual / norm
                };
                let mut history = response.report().history.clone();
                history.push(relative);
                let report = SolveReport::from_claim(
                    response.report().iters,
                    ResidualClaim::TrueEuclidean(relative),
                    tolerance,
                    history,
                );
                if !report.converged_euclidean() {
                    break;
                }
                responses.push(HarmonicWindingUnitResponse {
                    potential_per_current_h_m: response
                        .potential_per_current_h_m()
                        .iter()
                        .copied()
                        .map(C64::from_re)
                        .collect(),
                    report,
                });
            }
            if responses.len() == self.basis.winding_count() {
                return Ok(responses);
            }
        }
        let preconditioner = HarmonicJacobi(MovingJacobi::new(
            &self.tangent,
            &self.free,
            &self.conduction,
            omega.recip(),
            1.0,
            cancelled,
        )?);
        let callback = RefCell::new(&mut *cancelled);
        let interrupted = Cell::new(false);
        let check = || {
            if interrupted.get() || (*callback.borrow_mut())() {
                interrupted.set(true);
                true
            } else {
                false
            }
        };
        let n = self.tangent.nrows();
        let magnetic = ReducedMagnetic {
            tangent: &self.tangent,
            free: &self.free,
            full: RefCell::new(vec![0.0; n]),
            check: &check,
        };
        let conduction = ReducedConduction {
            conduction: &self.conduction,
            free: &self.free,
            omega,
            scratch: RefCell::new((vec![0.0; n], vec![0.0; n])),
            check: &check,
        };
        let operator = RealEquivalentComplexOp::new(&magnetic, &conduction)
            .map_err(|_| HarmonicWindingError::Invalid("harmonic block dimensions disagree"))?;
        let mut responses = Vec::with_capacity(self.loads.len());
        for load in &self.loads {
            if check() {
                return Err(HarmonicWindingError::Cancelled);
            }
            let mut rhs = vec![0.0; 2 * self.free.len()];
            for (out, &v) in rhs.iter_mut().zip(&self.free) {
                *out = load[v];
            }
            let norm = rhs.iter().fold(0.0_f64, |n, &v| n.hypot(v));
            finite(norm)?;
            let tolerance = if norm == 0.0 {
                self.controls.relative_tolerance
            } else {
                self.controls
                    .relative_tolerance
                    .max((self.controls.absolute_residual_a_m / norm).min(f64::MAX))
            };
            let mut solver = FgmresState::new(&rhs, self.restart);
            let solved = solver.run_controlled(
                &operator,
                &preconditioner,
                &rhs,
                tolerance,
                self.controls.max_cycles,
                &mut || check(),
            );
            if check() || solved.progress == SolverRunProgress::Paused {
                return Err(HarmonicWindingError::Cancelled);
            }
            if solved.progress != SolverRunProgress::Complete {
                return Err(HarmonicWindingError::Invalid(
                    "harmonic solver callback failed",
                ));
            }
            if !solved.report.converged_euclidean() {
                return Err(HarmonicWindingError::Nonconvergence(solved.report));
            }
            let mut applied = vec![0.0; rhs.len()];
            operator.apply(&solver.x, &mut applied);
            if check() {
                return Err(HarmonicWindingError::Cancelled);
            }
            let mut residual = 0.0_f64;
            for (&b, &a) in rhs.iter().zip(&applied) {
                residual = residual.hypot(finite(b - a)?);
            }
            let relative = if norm == 0.0 {
                if residual == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                residual / norm
            };
            let mut history = solved.report.history;
            history.push(relative);
            let report = SolveReport::from_claim(
                solved.report.iters,
                ResidualClaim::TrueEuclidean(relative),
                tolerance,
                history,
            );
            if !report.converged_euclidean() {
                return Err(HarmonicWindingError::Nonconvergence(report));
            }
            let mut potential = vec![C64::ZERO; n];
            for (i, &v) in self.free.iter().enumerate() {
                potential[v] = phasor(C64::new(solver.x[i], solver.x[self.free.len() + i]))?;
            }
            responses.push(HarmonicWindingUnitResponse {
                potential_per_current_h_m: potential,
                report,
            });
        }
        Ok(responses)
    }
}

impl<'b, 'm> PlanarHarmonicWindings<'b, 'm> {
    /// Solve a strictly positive angular frequency. No trial edits the admitted
    /// DC field. Exhausted solve, under-resolved skin depth or failed original
    /// field/terminal power pairing publishes no frequency matrix.
    pub fn at_frequency<'p>(
        &'p self,
        omega_rad_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HarmonicWindingFrequency<'p, 'b, 'm>> {
        poll(cancelled)?;
        if !positive(omega_rad_s) || !positive(omega_rad_s.recip()) {
            return Err(HarmonicWindingError::Invalid(
                "representably positive angular frequency required",
            ));
        }
        let mut maximum_cell_skin_ratio = 0.0_f64;
        for &factor in &self.skin_ratio_per_sqrt_rad_s {
            poll(cancelled)?;
            maximum_cell_skin_ratio =
                maximum_cell_skin_ratio.max(finite(factor * omega_rad_s.sqrt())?);
        }
        if maximum_cell_skin_ratio > self.controls.max_cell_skin_ratio {
            return Err(HarmonicWindingError::Invalid(
                "cell diameter exceeds admitted differential skin-depth ratio",
            ));
        }
        let unit_responses = self.solve_units(omega_rad_s, cancelled)?;
        let count = self.basis.winding_count();
        let mut impedance = vec![C64::ZERO; count * count];
        for (j, response) in unit_responses.iter().enumerate() {
            for (i, load) in self.loads.iter().enumerate() {
                let mut linkage = C64::ZERO;
                for (&value, &a) in load.iter().zip(&response.potential_per_current_h_m) {
                    poll(cancelled)?;
                    linkage = phasor(linkage + a.scale(value))?;
                }
                if i == j {
                    linkage.re = finite(linkage.re + self.basis.end_inductance_h()[i])?;
                }
                let mut z = phasor(C64::new(
                    omega_rad_s * linkage.im,
                    -omega_rad_s * linkage.re,
                ))?;
                if i == j {
                    z.re = finite(z.re + self.resistance_ohm[i])?;
                }
                impedance[i * count + j] = z;
            }
        }
        let mut reciprocity = 0.0_f64;
        for i in 0..count {
            for j in 0..i {
                poll(cancelled)?;
                reciprocity = reciprocity
                    .max(phasor(impedance[i * count + j] - impedance[j * count + i])?.abs());
            }
        }
        if reciprocity > self.controls.impedance_tolerance_ohm {
            return Err(HarmonicWindingError::Balance(
                "complex terminal reciprocity failed",
            ));
        }
        // Audit the COMPLETE complex power operator, including mutual terms:
        // Z_ij = R_ij + omega^2 a_i* C a_j - i omega a_i* K a_j.
        // Diagonal checks alone would miss phase-dependent cross-winding loss.
        let n = self.tangent.nrows();
        let mut k_real = vec![0.0; n];
        let mut k_imaginary = vec![0.0; n];
        let mut c_real = vec![0.0; n];
        let mut c_imaginary = vec![0.0; n];
        let mut maximum_power_pairing_defect_ohm = 0.0_f64;
        for (j, response) in unit_responses.iter().enumerate() {
            poll(cancelled)?;
            let real: Vec<_> = response
                .potential_per_current_h_m
                .iter()
                .map(|a| a.re)
                .collect();
            let imaginary: Vec<_> = response
                .potential_per_current_h_m
                .iter()
                .map(|a| a.im)
                .collect();
            self.tangent.spmv(&real, &mut k_real);
            self.tangent.spmv(&imaginary, &mut k_imaginary);
            {
                let callback = RefCell::new(&mut *cancelled);
                let check = || (*callback.borrow_mut())();
                if !self.conduction.apply(&real, &mut c_real, &check)
                    || !self.conduction.apply(&imaginary, &mut c_imaginary, &check)
                {
                    return Err(HarmonicWindingError::Cancelled);
                }
            }
            for (i, receiver) in unit_responses.iter().enumerate() {
                let mut magnetic = C64::ZERO;
                let mut loss = C64::ZERO;
                for (v, &a) in receiver.potential_per_current_h_m.iter().enumerate() {
                    poll(cancelled)?;
                    magnetic = phasor(magnetic + a.conj() * C64::new(k_real[v], k_imaginary[v]))?;
                    loss = phasor(loss + a.conj() * C64::new(c_real[v], c_imaginary[v]))?;
                }
                if i == j {
                    magnetic.re = finite(magnetic.re + self.basis.end_inductance_h()[i])?;
                }
                let mut from_field = phasor(
                    loss.scale(omega_rad_s).scale(omega_rad_s)
                        + C64::new(omega_rad_s * magnetic.im, -omega_rad_s * magnetic.re),
                )?;
                if i == j {
                    from_field.re = finite(from_field.re + self.resistance_ohm[i])?;
                }
                maximum_power_pairing_defect_ohm = maximum_power_pairing_defect_ohm
                    .max(phasor(from_field - impedance[i * count + j])?.abs());
            }
        }
        if maximum_power_pairing_defect_ohm > self.controls.impedance_tolerance_ohm {
            return Err(HarmonicWindingError::Balance(
                "terminal impedance disagrees with resolved field power",
            ));
        }
        let effective_inductance_h = impedance
            .iter()
            .map(|z| finite(-z.im / omega_rad_s))
            .collect::<Result<Vec<_>>>()?;
        poll(cancelled)?;
        Ok(HarmonicWindingFrequency {
            prepared: self,
            omega_rad_s,
            impedance_ohm: impedance,
            effective_inductance_h,
            unit_responses,
            maximum_reciprocity_defect_ohm: reciprocity,
            maximum_power_pairing_defect_ohm,
            maximum_cell_skin_ratio,
        })
    }
}
