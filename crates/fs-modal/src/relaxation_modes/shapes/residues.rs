//! Physical simple-pole force residues, retaining nonsymmetric mechanics.
//!
//! For T(s0)q=0 and T(s0)^T ell=0, the resolvent has residue
//! q ell^T / (ell^T T'(s0)q). The transpose is bilinear, not Hermitian:
//! gyroscopic and other nonsymmetric terms must not replace ell by q or q*.
//! Maxwell stress residues use their original rational branch response.

use super::{
    MaxwellMode, MaxwellModeError, PoleOptions, PreparedMaxwellModes, SourceSum, budget, finite,
    invalid, poll, quadratic_fraction, quadratic_power_scale, scale_complex,
};
use fs_math::{c64::C64, det};

/// A resolved right/transpose mechanical pair and its physical normalization.
///
/// Borrows the original prepared source; subsequent forces do not repeat an
/// eigensolve. The pair represents ONE local Laurent coefficient, not a
/// complete forced response, certified modal expansion or time realization.
#[derive(Clone, Debug)]
pub struct MaxwellModeResidue<'p, 'a> {
    source: &'p PreparedMaxwellModes<'a>,
    right: MaxwellMode,
    left: Vec<C64>,
    normalization: C64,
    alignment: f64,
    cancellation_ratio: f64,
    left_residual: (f64, f64),
}

/// Coefficients of 1/(s-s0) in the physical response to an original force f.
///
/// The force is conjugate to the original displacement coordinates. It is a
/// complex Laplace-domain input amplitude, not an assumed real time history.
/// No analytic/background response or contribution from other poles is added.
#[derive(Clone, Debug)]
pub struct MaxwellForceResidue {
    /// The local pole s0 [1/s], with convention exp(s0*t).
    pub laplace_frequency_per_s: C64,
    /// Residue of the original physical displacement transfer T(s)^-1 f.
    pub displacement: Vec<C64>,
    /// Residue of s*T(s)^-1 f, equal to s0 times displacement residue.
    pub velocity: Vec<C64>,
    /// Residue of every original branch stress, in original source order.
    pub branch_stress: Vec<C64>,
    /// ell^T f / (ell^T T'(s0) q), with the published right mode's gauge.
    pub participation: C64,
}

/// One original branch's derivative with respect to a caller's parameter.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellBranchDirection<'a> {
    /// Derivative of the source branch stiffness, in its original units.
    pub stiffness: f64,
    /// Derivative of ln(tau); equivalently d(tau)/tau.
    pub log_relaxation_time: f64,
    /// Derivative of every original strain-direction entry, in source order.
    pub strain: &'a [f64],
}

/// An infinitesimal source direction at fixed physical coordinate/time units.
///
/// All entries are derivatives, not replacement coefficients. Explicit zeros
/// represent unchanged inputs. No solver iteration, coordinate equilibration,
/// material fit or mode-shape normalization is differentiated.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellPoleDirection<'a> {
    /// Symmetric n-by-n original mass derivative, dense row-major.
    pub mass: &'a [f64],
    /// Full n-by-n original damping derivative, including gyroscopic entries.
    pub damping: &'a [f64],
    /// Full n-by-n RELAXED stiffness derivative; branch changes remain separate.
    pub relaxed_stiffness: &'a [f64],
    /// One derivative record per original Maxwell branch, in the same order.
    pub branches: &'a [MaxwellBranchDirection<'a>],
}

impl<'a> PreparedMaxwellModes<'a> {
    /// Refine a nearby simple pole, its transpose mode, and force normalization.
    ///
    /// Both seeds are ORIGINAL-coordinate displacements; the left seed targets
    /// T(s)^T ell=0 without conjugation. The two shared Newton solves retain
    /// the original C and K, including gyroscopic/active coefficients. The
    /// transpose mode is independently checked at the RIGHT pole, so a left
    /// iteration converging to another pole cannot silently define a residue.
    ///
    /// `minimum_derivative_coupling` is a caller-authored value strictly
    /// between zero and one. Both the equilibrated derivative alignment and
    /// the original-source derivative cancellation ratio must meet it.
    /// These are numerical conditioning checks, not eigenvalue enclosures.
    /// Original `PoleOptions` ceilings admit BOTH solves and retained results
    /// together; too-small paired resources refuse before either solve.
    ///
    /// # Errors
    /// All ordinary mode refusals, insufficient paired resources, unresolved
    /// left/right derivative normalization, or cancellation. Material poles,
    /// dark modes and multiple/defective modes remain outside this local API.
    pub fn refine_residue(
        &self,
        laplace_seed_per_s: C64,
        displacement_seed: &[C64],
        transpose_displacement_seed: &[C64],
        minimum_derivative_coupling: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellModeResidue<'_, 'a>, MaxwellModeError> {
        poll(&mut cancelled)?;
        if !minimum_derivative_coupling.is_finite()
            || minimum_derivative_coupling <= 0.0
            || minimum_derivative_coupling >= 1.0
        {
            return Err(invalid(
                "residue derivative-coupling floor must lie in (0,1)",
            ));
        }
        let options = self.paired_options()?;
        let right = self.refine_orientation(
            laplace_seed_per_s,
            displacement_seed,
            false,
            &options,
            &mut cancelled,
        )?;
        let left = self.refine_orientation(
            right.laplace_frequency_per_s,
            transpose_displacement_seed,
            true,
            &options,
            &mut cancelled,
        )?;
        let scaled_pole = scale_complex(right.laplace_frequency_per_s, -self.time_exponent)?;
        let (left_velocity, left_stress) = self.recover_response(
            scaled_pole,
            right.laplace_frequency_per_s,
            &left.displacement,
            &mut cancelled,
        )?;
        let left_residual = self.original_residual(
            right.laplace_frequency_per_s,
            &left.displacement,
            &left_velocity,
            &left_stress,
            true,
            &mut cancelled,
        )?;
        if left_residual.0.max(left_residual.1) > options.residual_tolerance {
            return Err(MaxwellModeError::Residual {
                mechanical: left_residual.0,
                relaxation: left_residual.1,
                limit: options.residual_tolerance,
            });
        }
        let normalization =
            self.source_derivative_pair(scaled_pole, &left.displacement, &right, &mut cancelled)?;
        let cancellation_ratio = normalization.relative()?;
        let alignment = self.derivative_alignment(
            scaled_pole,
            &left.displacement,
            &right.displacement,
            &mut cancelled,
        )?;
        if alignment < minimum_derivative_coupling
            || cancellation_ratio < minimum_derivative_coupling
        {
            return Err(MaxwellModeError::ResidueNormalization {
                alignment,
                cancellation_ratio,
                minimum: minimum_derivative_coupling,
            });
        }
        let normalization = normalization.value()?;
        if normalization == C64::ZERO {
            return Err(invalid("nonzero physical residue normalization required"));
        }
        poll(&mut cancelled)?;
        Ok(MaxwellModeResidue {
            source: self,
            right,
            left: left.displacement,
            normalization,
            alignment,
            cancellation_ratio,
            left_residual,
        })
    }

    fn paired_options(&self) -> Result<PoleOptions, MaxwellModeError> {
        let n = self.n as u128;
        let branches = self.branches.len() as u128;
        let d = n + 1;
        let result_bytes = 128 * (n + branches);
        let bytes = self.preparation_bounds.0 as u128 + result_bytes;
        budget(
            "paired numerical payload bytes",
            bytes + 128 * d * d + 256 * d,
            self.resource_limits.0 as u128,
        )?;
        // Twice the existing preparation/evaluation bound conservatively
        // covers both refinements, although preparation itself is shared.
        // The last term covers transpose traversal, source contraction and
        // the common-pole audit independently of the Newton matrices.
        let work = (2 * u128::from(self.preparation_bounds.1))
            .saturating_add(512 * (n * n + n + 1) * (branches + 1));
        budget(
            "paired arithmetic work",
            work.saturating_add(2 * d * d * (d + 16) * self.options.max_iterations as u128),
            u128::from(self.resource_limits.1),
        )?;
        let mut options = self.options;
        options.max_dense_bytes = self.resource_limits.0 - bytes as usize;
        options.max_work = (self.resource_limits.1 - work as u64) / 2;
        Ok(options)
    }

    fn source_derivative_pair(
        &self,
        scaled_pole: C64,
        left: &[C64],
        right: &MaxwellMode,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SourceSum, MaxwellModeError> {
        let mut sum = SourceSum::default();
        for (i, &ell) in left.iter().enumerate() {
            poll(cancelled)?;
            for j in 0..self.n {
                sum.add(
                    ell,
                    right.velocity[j],
                    &[2.0, self.source[0][i * self.n + j]],
                )?;
                sum.add(
                    ell,
                    right.displacement[j],
                    &[self.source[1][i * self.n + j]],
                )?;
            }
        }
        for (index, branch) in self.branches.iter().enumerate() {
            let (_, slope) = self.weights(scaled_pole, index)?;
            for (i, &ell) in left.iter().enumerate() {
                poll(cancelled)?;
                for j in 0..self.n {
                    // weights differentiates in the binary-scaled frequency.
                    // Apply 1/rho only after all source products are factored.
                    sum.add_weighted(
                        ell,
                        right.displacement[j],
                        slope,
                        &[branch.stiffness, branch.strain[i], branch.strain[j]],
                        -self.time_exponent,
                    )?;
                }
            }
        }
        Ok(sum)
    }

    fn derivative_alignment(
        &self,
        scaled_pole: C64,
        left: &[C64],
        right: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, MaxwellModeError> {
        let left = self.equilibrate_vector(left)?;
        let right = self.equilibrate_vector(right)?;
        let derivative = self.evaluate(scaled_pole, cancelled)?.derivative;
        let mut matrix_norm = 0.0;
        let mut pairing = SourceSum::default();
        for (i, &ell) in left.iter().enumerate() {
            poll(cancelled)?;
            for (j, &q) in right.iter().enumerate() {
                let value = derivative[i * self.n + j];
                matrix_norm = det::hypot(matrix_norm, value.abs());
                pairing.add_weighted(ell, q, value, &[], 0)?;
            }
        }
        let norm = |vector: &[C64]| {
            vector
                .iter()
                .fold(0.0, |n, value| det::hypot(n, value.abs()))
        };
        if matrix_norm == 0.0 {
            return Ok(0.0);
        }
        pairing.relative_to(&[norm(&left), norm(&right), matrix_norm])
    }
}

impl MaxwellModeResidue<'_, '_> {
    /// The original physical right mode and its independent source residuals.
    #[must_use]
    pub fn right_mode(&self) -> &MaxwellMode {
        &self.right
    }

    /// Unit-modal-mass transpose displacement ell, with T(s0)^T ell=0.
    /// There is no conjugation in ell^T f or in the derivative normalization.
    #[must_use]
    pub fn transpose_displacement(&self) -> &[C64] {
        &self.left
    }

    /// Original-unit bilinear derivative ell^T T'(s0) q.
    #[must_use]
    pub fn derivative_normalization(&self) -> C64 {
        self.normalization
    }

    /// |ell^T T' q|/(||ell|| ||T'||_F ||q||) in binary-equilibrated coordinates.
    #[must_use]
    pub fn derivative_alignment(&self) -> f64 {
        self.alignment
    }

    /// Original derivative contraction divided by its absolute source-term sum.
    #[must_use]
    pub fn source_cancellation_ratio(&self) -> f64 {
        self.cancellation_ratio
    }

    /// Original transposed mechanics and branch residuals at the RIGHT pole.
    #[must_use]
    pub fn transpose_source_residuals(&self) -> (f64, f64) {
        self.left_residual
    }

    /// Differentiate the local physical pole through original source inputs.
    ///
    /// ds/dp = -ell^T (dT/dp) q / (ell^T T' q). The complete mass, damping,
    /// relaxed stiffness, branch stiffness, log relaxation time and strain
    /// direction derivatives participate. Original products use factored DD
    /// sums without building a perturbed operator or repeating an eigensolve.
    /// This is an infinitesimal derivative of a simple pole; it does not bound
    /// finite parameter steps or resolve crossings and repeated modes.
    ///
    /// # Errors
    /// Invalid derivative dimensions/values or mass asymmetry, unrepresentable
    /// physical arithmetic, or cancellation. The retained pair is unchanged.
    #[allow(clippy::float_cmp, clippy::too_many_lines)] // One atomic source derivative; mass remains exactly symmetric.
    pub fn pole_derivative(
        &self,
        direction: &MaxwellPoleDirection<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<C64, MaxwellModeError> {
        poll(&mut cancelled)?;
        let n = self.source.n;
        if [
            direction.mass,
            direction.damping,
            direction.relaxed_stiffness,
        ]
        .iter()
        .any(|a| a.len() != n * n || a.iter().any(|x| !x.is_finite()))
            || direction.branches.len() != self.source.branches.len()
        {
            return Err(invalid(
                "complete finite source matrix and branch derivatives required",
            ));
        }
        let mut numerator = SourceSum::default();
        for (i, &ell) in self.left.iter().enumerate() {
            poll(&mut cancelled)?;
            for j in 0..n {
                if direction.mass[i * n + j] != direction.mass[j * n + i] {
                    return Err(invalid("mass derivative must be symmetric"));
                }
                numerator.add_weighted(
                    ell,
                    self.right.velocity[j],
                    self.right.laplace_frequency_per_s,
                    &[direction.mass[i * n + j]],
                    0,
                )?;
                numerator.add(ell, self.right.velocity[j], &[direction.damping[i * n + j]])?;
                numerator.add(
                    ell,
                    self.right.displacement[j],
                    &[direction.relaxed_stiffness[i * n + j]],
                )?;
            }
        }
        let scaled_pole = scale_complex(
            self.right.laplace_frequency_per_s,
            -self.source.time_exponent,
        )?;
        for (index, (branch, delta)) in self
            .source
            .branches
            .iter()
            .zip(direction.branches)
            .enumerate()
        {
            poll(&mut cancelled)?;
            if !delta.stiffness.is_finite()
                || !delta.log_relaxation_time.is_finite()
                || delta.strain.len() != n
                || delta.strain.iter().any(|x| !x.is_finite())
            {
                return Err(invalid(
                    "one finite stiffness, log-time and strain derivative per branch required",
                ));
            }
            let (weight, slope) = self.source.weights(scaled_pole, index)?;
            for (i, &ell) in self.left.iter().enumerate() {
                poll(&mut cancelled)?;
                for j in 0..n {
                    numerator.add_weighted(
                        ell,
                        self.right.displacement[j],
                        weight,
                        &[delta.stiffness, branch.strain[i], branch.strain[j]],
                        0,
                    )?;
                    numerator.add_weighted(
                        ell,
                        self.right.displacement[j],
                        weight,
                        &[branch.stiffness, delta.strain[i], branch.strain[j]],
                        0,
                    )?;
                    numerator.add_weighted(
                        ell,
                        self.right.displacement[j],
                        weight,
                        &[branch.stiffness, branch.strain[i], delta.strain[j]],
                        0,
                    )?;
                    // d/dln(tau) of s*tau/(1+s*tau) equals s times its
                    // physical frequency derivative. Reuse v=s*q so neither
                    // a denominator square nor an unscaled product is formed.
                    numerator.add_weighted(
                        ell,
                        self.right.velocity[j],
                        slope,
                        &[
                            branch.stiffness,
                            branch.strain[i],
                            branch.strain[j],
                            delta.log_relaxation_time,
                        ],
                        -self.source.time_exponent,
                    )?;
                }
            }
        }
        let result = -numerator.divided_by(self.normalization)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Apply one original complex force without a dense inverse or eigensolve.
    /// Work and result storage are linear in coordinates plus branch count.
    /// Fixed source order and factored DD contractions are deterministic.
    ///
    /// # Errors
    /// Invalid force dimensions/values, unrepresentable physical coefficients,
    /// or cancellation. No partial response is returned on failure.
    pub fn force_residue(
        &self,
        force: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellForceResidue, MaxwellModeError> {
        poll(&mut cancelled)?;
        if force.len() != self.source.n || force.iter().any(|&f| !finite(f)) {
            return Err(invalid(
                "one finite complex force per original coordinate required",
            ));
        }
        let mut input = SourceSum::default();
        for (&ell, &force) in self.left.iter().zip(force) {
            poll(&mut cancelled)?;
            input.add(ell, force, &[])?;
        }
        let participation = input.divided_by(self.normalization)?;
        let mut apply = |mode: &[C64]| -> Result<Vec<C64>, MaxwellModeError> {
            let mut result = Vec::with_capacity(mode.len());
            for &value in mode {
                poll(&mut cancelled)?;
                let mut product = SourceSum::default();
                product.add(value, participation, &[])?;
                result.push(product.value()?);
            }
            Ok(result)
        };
        let displacement = apply(&self.right.displacement)?;
        let velocity = apply(&self.right.velocity)?;
        let branch_stress = apply(&self.right.branch_stress)?;
        poll(&mut cancelled)?;
        Ok(MaxwellForceResidue {
            laplace_frequency_per_s: self.right.laplace_frequency_per_s,
            displacement,
            velocity,
            branch_stress,
            participation,
        })
    }
}

impl SourceSum {
    fn divided_by(&self, divisor: C64) -> Result<C64, MaxwellModeError> {
        let power = quadratic_fraction(divisor.re.abs().max(divisor.im.abs())).1;
        let divisor = scale_complex(divisor, -power)?;
        let quotient = C64::new(self.re.to_f64(), self.im.to_f64()) / divisor;
        scale_complex(quotient, self.exponent - power)
    }

    fn relative_to(&self, factors: &[f64]) -> Result<f64, MaxwellModeError> {
        let mut ratio = det::hypot(self.re.to_f64(), self.im.to_f64());
        let mut power = self.exponent;
        for &factor in factors {
            if !factor.is_finite() || factor <= 0.0 {
                return Err(invalid(
                    "positive finite residue conditioning scale required",
                ));
            }
            let (fraction, exponent) = quadratic_fraction(factor);
            ratio /= fraction;
            power -= exponent;
        }
        let ratio = quadratic_power_scale(ratio, power);
        if ratio.is_finite() {
            Ok(ratio)
        } else {
            Err(invalid("finite residue conditioning ratio required"))
        }
    }
}
