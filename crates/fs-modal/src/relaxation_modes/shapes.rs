//! Local mechanical right modes of the original rational Maxwell pencil.
//!
//! The shared bordered Newton owner receives the full analytic T(s) and T'(s).
//! No fitted damping ratio, companion eigenvector or real modal projection is
//! used. Recovery publishes original-coordinate displacement, velocity and
//! every original branch stress, then audits the original first-order source
//! equations independently of the Newton matrices.
//!
//! This is a nearby simple mechanical mode, not a complete spectrum, a decay
//! or shape-error certificate, left/right force residues or a time stepper.
//! The rational domain excludes every s=-1/tau, even for a dark strain channel.
//! Pure material modes with zero mechanical displacement cannot be normalized
//! here; use the full-state spectrum owner to retain their pole values.

use super::{MaxwellDampingBranch, coupling_parts};
use crate::{quadratic_fraction, quadratic_power_scale};
use fs_la::factor::{Cholesky, cholesky};
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, PoleEvaluation, newton_pole_with_cancel,
};
use fs_math::{c64::C64, dd::Dd, det};

/// Shared local accuracy/resource controls; spectral_scale is physical [1/s].
pub use fs_la::nonlinear_eig::PoleOptions;

/// A locally refined mechanical right mode, using exp(s*t).
#[derive(Clone, Debug)]
pub struct MaxwellMode {
    /// Original physical Laplace exponent [1/s], including active growth.
    pub laplace_frequency_per_s: C64,
    /// Original displacement with q^H M q=1. Its largest component is real
    /// and positive, with the lowest index breaking exact magnitude ties.
    pub displacement: Vec<C64>,
    /// Physical velocity s*q, in the same phase and normalization.
    pub velocity: Vec<C64>,
    /// Stress conjugate to each supplied branch strain, in source order:
    /// k_j*s*tau_j/(1+s*tau_j) * b_j^T q. These are mode amplitudes.
    pub branch_stress: Vec<C64>,
    /// Largest componentwise backward residual of M*s*v+C*v+K*q+sum(b*sigma).
    pub mechanical_relative_residual: f64,
    /// Largest componentwise backward residual of (1+s*tau)*sigma-k*tau*b^T*v.
    pub relaxation_relative_residual: f64,
    /// Shared Newton owner's residual before physical recovery.
    pub solver_relative_residual: f64,
    /// Last correction divided by max(reference rate, pole magnitude).
    pub relative_correction: f64,
    /// Number of Newton operator evaluations.
    pub iterations: usize,
}

/// A refusal publishes no partial physical mode and changes no source.
#[derive(Clone, Debug, PartialEq)]
pub enum MaxwellModeError {
    /// Invalid shape, finite domain, source data or representability.
    Invalid(&'static str),
    /// Mass is not exactly symmetric and numerically positive definite.
    MassNotPositiveDefinite,
    /// The declared memory or arithmetic ceiling cannot admit this solve.
    Budget {
        /// Admitted resource.
        quantity: &'static str,
        /// Conservative requested count.
        requested: u128,
        /// Caller ceiling.
        cap: u128,
    },
    /// A material denominator vanishes; pure material dark modes are outside
    /// this displacement-based rational realization too.
    MaterialPole {
        /// Original branch index.
        branch: usize,
        /// Physical Laplace argument [1/s].
        at: C64,
    },
    /// Shared local solve refusal, including unresolved multiplicity.
    Solve(NonlinearEigError),
    /// Published physical amplitudes failed the independent source equations.
    Residual {
        /// Mechanical source-row residual.
        mechanical: f64,
        /// Original branch-history residual.
        relaxation: f64,
        /// Caller ceiling for both residuals.
        limit: f64,
    },
    /// Cancellation observed before publication.
    Cancelled,
}
impl core::fmt::Display for MaxwellModeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-MAXWELL-MODE: {self:?}")
    }
}
impl std::error::Error for MaxwellModeError {}
impl From<NonlinearEigError> for MaxwellModeError {
    fn from(error: NonlinearEigError) -> Self {
        if error == NonlinearEigError::Cancelled {
            Self::Cancelled
        } else {
            Self::Solve(error)
        }
    }
}

#[derive(Debug)]
struct PreparedBranch {
    // sqrt(k)*D*b/rho: energy-normalized direction, no n-by-n branch cache.
    direction: Vec<f64>,
    rate: f64,
    norm_squared: f64,
}

/// Immutable source mechanics and Maxwell branches for repeated local modes.
///
/// M is SPD; C and relaxed K retain all entries, including nonsymmetric,
/// gyroscopic and active terms. Branch stiffness must not be added to K again.
/// The existing PoleOptions ceilings include branch storage and assembly work,
/// as well as one full refinement. They do not bound a caller's total number
/// of later refinements or include caller-owned sources/previous results,
/// allocator metadata or process stacks.
#[derive(Debug)]
pub struct PreparedMaxwellModes<'a> {
    n: usize,
    source: [&'a [f64]; 3],
    branches: &'a [MaxwellDampingBranch],
    scaled: [Vec<f64>; 3],
    memory: Vec<PreparedBranch>,
    mass_factor: Cholesky,
    coordinate_exponents: Vec<i32>,
    time_exponent: i32,
    coefficient_norms: [f64; 3],
    options: PoleOptions,
}

impl<'a> PreparedMaxwellModes<'a> {
    /// Admit original dense row-major M/C/K and all scalar strain channels.
    /// Binary congruence/rate scaling and one mass factor are reused. Dense
    /// Cholesky is a bounded indivisible stage; other traversals poll by row
    /// or branch. No source coefficient or relaxation time is fitted.
    ///
    /// # Errors
    /// Invalid sources/controls, insufficient budgets, failed SPD mass
    /// admission, unrepresentable scaled coefficients, or cancellation.
    #[allow(clippy::too_many_lines, clippy::float_cmp)] // One atomic admission; exact original mass symmetry.
    pub fn new(
        mass: &'a [f64],
        damping: &'a [f64],
        stiffness: &'a [f64],
        n: usize,
        branches: &'a [MaxwellDampingBranch],
        mut options: PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaxwellModeError> {
        poll(&mut cancelled)?;
        if n == 0
            || n > MAX_NONLINEAR_POLE_DIMENSION
            || n > options.max_dimension
            || options.max_iterations == 0
            || !options.spectral_scale.is_finite()
            || options.spectral_scale <= 0.0
            || [options.residual_tolerance, options.step_tolerance]
                .iter()
                .any(|t| !t.is_finite() || *t <= 0.0 || *t >= 1.0)
        {
            return Err(invalid(
                "positive bounded dimension, rate and tolerances required",
            ));
        }
        let size = n as u128;
        let count = branches.len() as u128;
        let iterations = options.max_iterations as u128;
        let d = size + 1;
        let bytes = 64 * size * size
            + 512 * size
            + (8 * size + 256) * count
            + if n > 32 { 2 * 1024 * 1024 } else { 0 };
        budget(
            "numerical payload bytes",
            bytes + 128 * d * d + 256 * d,
            options.max_dense_bytes as u128,
        )?;
        let work = (8 * size * size * size + 64 * size * size * (iterations + 1)).saturating_add(
            (256 * (size * size + size + 1))
                .saturating_mul(count)
                .saturating_mul(iterations + 1),
        );
        budget(
            "arithmetic work",
            work.saturating_add(d * d * (d + 16) * iterations),
            options.max_work as u128,
        )?;
        options.max_dense_bytes -= bytes as usize;
        options.max_work -= work as u64;
        let source = [mass, damping, stiffness];
        if source
            .iter()
            .any(|a| a.len() != n * n || a.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid(
                "complete finite n-by-n source coefficients required",
            ));
        }
        for i in 0..n {
            poll(&mut cancelled)?;
            if mass[i * n + i] <= 0.0 || (0..i).any(|j| mass[i * n + j] != mass[j * n + i]) {
                return Err(MaxwellModeError::MassNotPositiveDefinite);
            }
        }
        for branch in branches {
            poll(&mut cancelled)?;
            if branch.strain.len() != n
                || branch.strain.iter().any(|x| !x.is_finite())
                || !branch.stiffness.is_finite()
                || branch.stiffness <= 0.0
                || !branch.relaxation_time_s.is_finite()
                || branch.relaxation_time_s <= 0.0
            {
                return Err(invalid(
                    "finite strain channels and positive branch stiffness/times required",
                ));
            }
        }
        let coordinate_exponents: Vec<_> = (0..n)
            .map(|i| quadratic_fraction(mass[i * n + i]).1.div_euclid(2))
            .collect();
        let mut time_exponent = quadratic_fraction(options.spectral_scale).1;
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                let coordinate = coordinate_exponents[i] + coordinate_exponents[j];
                if damping[i * n + j] != 0.0 {
                    time_exponent =
                        time_exponent.max(quadratic_fraction(damping[i * n + j]).1 - coordinate);
                }
                if stiffness[i * n + j] != 0.0 {
                    time_exponent = time_exponent.max(
                        (quadratic_fraction(stiffness[i * n + j]).1 - coordinate + 1).div_euclid(2),
                    );
                }
            }
        }
        for branch in branches {
            poll(&mut cancelled)?;
            for (i, &strain) in branch.strain.iter().enumerate() {
                if strain != 0.0 {
                    time_exponent = time_exponent
                        .max(coupling_parts(branch.stiffness, strain).1 - coordinate_exponents[i]);
                }
            }
        }
        options.spectral_scale = scale(options.spectral_scale, -time_exponent)?;
        let mut scaled = [
            Vec::with_capacity(n * n),
            Vec::with_capacity(n * n),
            Vec::with_capacity(n * n),
        ];
        let mut coefficient_norms = [0.0; 3];
        for order in 0..3 {
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..n {
                    let value = scale(
                        source[order][i * n + j],
                        -coordinate_exponents[i]
                            - coordinate_exponents[j]
                            - order as i32 * time_exponent,
                    )?;
                    scaled[order].push(value);
                    coefficient_norms[order] = det::hypot(coefficient_norms[order], value);
                }
            }
        }
        let mass_factor =
            cholesky(&scaled[0], n).map_err(|_| MaxwellModeError::MassNotPositiveDefinite)?;
        let mut memory = Vec::with_capacity(branches.len());
        for branch in branches {
            poll(&mut cancelled)?;
            let mut direction = Vec::with_capacity(n);
            let mut norm = 0.0;
            for (i, &strain) in branch.strain.iter().enumerate() {
                let value = if strain == 0.0 {
                    0.0
                } else {
                    let (fraction, power) = coupling_parts(branch.stiffness, strain);
                    scale(fraction, power - coordinate_exponents[i] - time_exponent)?
                };
                direction.push(value);
                norm = det::hypot(norm, value);
            }
            let norm_squared = norm * norm;
            let (tau_fraction, tau_power) = quadratic_fraction(branch.relaxation_time_s);
            let rate = scale(1.0 / tau_fraction, -tau_power - time_exponent)?;
            if !norm_squared.is_finite() || (norm != 0.0 && norm_squared == 0.0) {
                return Err(invalid("scaled branch stiffness must be representable"));
            }
            memory.push(PreparedBranch {
                direction,
                rate,
                norm_squared,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            n,
            source,
            branches,
            scaled,
            memory,
            mass_factor,
            coordinate_exponents,
            time_exponent,
            coefficient_norms,
            options,
        })
    }

    /// Refine a nearby physical pole with a nonzero original displacement seed.
    /// Both residual and Newton correction must pass; the recovered physical
    /// mechanics and branch histories are independently checked before return.
    ///
    /// # Errors
    /// Invalid/zero seeds, material singularities, shared Newton refusals
    /// (including multiple/defective modes), unrepresentable recovery, failed
    /// original-equation residuals, or cancellation. No convergence from an
    /// arbitrary seed and no relative weak-decay accuracy are promised.
    pub fn refine(
        &self,
        laplace_seed_per_s: C64,
        displacement_seed: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellMode, MaxwellModeError> {
        poll(&mut cancelled)?;
        let seed = self.equilibrate_vector(displacement_seed)?;
        let value = scale_complex(laplace_seed_per_s, -self.time_exponent)?;
        let mut evaluation_error = None;
        let solved = newton_pole_with_cancel(
            self.n,
            value,
            &seed,
            &self.options,
            |s, cancel| {
                self.evaluate(s, cancel).map_err(|error| {
                    let mapped = if error == MaxwellModeError::Cancelled {
                        NonlinearEigError::Cancelled
                    } else {
                        NonlinearEigError::FunctionRefused { at: s }
                    };
                    evaluation_error = Some(error);
                    mapped
                })
            },
            &mut cancelled,
        );
        let pole = solved.map_err(|error| evaluation_error.unwrap_or_else(|| error.into()))?;
        let displacement = self.physical_displacement(&pole.vector, &mut cancelled)?;
        let laplace_frequency_per_s = scale_complex(pole.value, self.time_exponent)?;
        let velocity: Vec<_> = displacement
            .iter()
            .map(|&q| laplace_frequency_per_s * q)
            .collect();
        if velocity.iter().any(|&v| !finite(v))
            || velocity.iter().zip(&displacement).any(|(&v, &q)| {
                laplace_frequency_per_s != C64::ZERO && q != C64::ZERO && v == C64::ZERO
            })
        {
            return Err(invalid("physical velocity must be representable"));
        }
        let mut branch_stress = Vec::with_capacity(self.branches.len());
        for (index, branch) in self.branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let (weight, _) = self.weights(pole.value, index)?;
            let mut stress = SourceSum::default();
            for (&strain, &q) in branch.strain.iter().zip(&displacement) {
                stress.add(q, weight, &[branch.stiffness, strain])?;
            }
            branch_stress.push(stress.value()?);
        }
        let (mechanical, relaxation) = self.original_residual(
            laplace_frequency_per_s,
            &displacement,
            &velocity,
            &branch_stress,
            &mut cancelled,
        )?;
        if mechanical.max(relaxation) > self.options.residual_tolerance {
            return Err(MaxwellModeError::Residual {
                mechanical,
                relaxation,
                limit: self.options.residual_tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(MaxwellMode {
            laplace_frequency_per_s,
            displacement,
            velocity,
            branch_stress,
            mechanical_relative_residual: mechanical,
            relaxation_relative_residual: relaxation,
            solver_relative_residual: pole.residual,
            relative_correction: pole.relative_step,
            iterations: pole.iterations,
        })
    }

    fn weights(&self, s: C64, index: usize) -> Result<(C64, C64), MaxwellModeError> {
        let rate = self.memory[index].rate;
        let denominator = s + C64::from_re(rate);
        if denominator == C64::ZERO {
            return Err(MaxwellModeError::MaterialPole {
                branch: index,
                at: scale_complex(s, self.time_exponent)?,
            });
        }
        // r=s/(s+a), r'=a/(s+a)^2. Sequential divisions avoid forming an
        // unscaled complex denominator square or a sampled loss ratio.
        let value = s / denominator;
        let derivative = (C64::from_re(rate) / denominator) / denominator;
        if !finite(value)
            || !finite(derivative)
            || derivative == C64::ZERO
            || (s != C64::ZERO && value == C64::ZERO)
        {
            return Err(invalid(
                "rational branch value and derivative must be representable",
            ));
        }
        Ok((value, derivative))
    }

    fn evaluate(
        &self,
        s: C64,
        cancel: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, MaxwellModeError> {
        let mut matrix = vec![C64::ZERO; self.n * self.n];
        let mut derivative = vec![C64::ZERO; self.n * self.n];
        for i in 0..self.n {
            poll(cancel)?;
            for j in 0..self.n {
                let index = i * self.n + j;
                let [m, c, k] = self.scaled.each_ref().map(|a| a[index]);
                matrix[index] = s * (s.scale(m) + C64::from_re(c)) + C64::from_re(k);
                derivative[index] = s.scale(m).scale(2.0) + C64::from_re(c);
            }
        }
        let radius = s.abs().max(self.options.spectral_scale);
        let [m, c, k] = self.coefficient_norms;
        let mut residual_scale = (m * radius + c) * radius + k;
        for (index, branch) in self.memory.iter().enumerate() {
            let (weight, slope) = self.weights(s, index)?;
            for i in 0..self.n {
                poll(cancel)?;
                for j in 0..self.n {
                    let entry = branch.direction[i] * branch.direction[j];
                    if entry == 0.0 && branch.direction[i] != 0.0 && branch.direction[j] != 0.0 {
                        return Err(invalid("branch outer product must be representable"));
                    }
                    matrix[i * self.n + j] = matrix[i * self.n + j] + weight.scale(entry);
                    derivative[i * self.n + j] = derivative[i * self.n + j] + slope.scale(entry);
                }
            }
            residual_scale +=
                branch.norm_squared * (radius / (s + C64::from_re(branch.rate)).abs());
        }
        if !residual_scale.is_finite()
            || residual_scale <= 0.0
            || matrix.iter().chain(&derivative).any(|&z| !finite(z))
        {
            return Err(invalid(
                "rational operator and residual scale must be finite",
            ));
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale,
        })
    }

    fn equilibrate_vector(&self, vector: &[C64]) -> Result<Vec<C64>, MaxwellModeError> {
        if vector.len() != self.n || vector.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "one finite complex displacement per source coordinate required",
            ));
        }
        let mut largest = i32::MIN;
        for (i, z) in vector.iter().enumerate() {
            for value in [z.re, z.im] {
                if value != 0.0 {
                    largest =
                        largest.max(quadratic_fraction(value).1 + self.coordinate_exponents[i]);
                }
            }
        }
        if largest == i32::MIN {
            return Err(invalid(
                "nonzero mechanical displacement seed required; dark modes are excluded",
            ));
        }
        vector
            .iter()
            .enumerate()
            .map(|(i, &z)| scale_complex(z, self.coordinate_exponents[i] - largest))
            .collect()
    }

    fn physical_displacement(
        &self,
        vector: &[C64],
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, MaxwellModeError> {
        let mut norm = 0.0;
        for i in 0..self.n {
            poll(cancel)?;
            let mut image = C64::ZERO;
            for (j, &value) in vector.iter().enumerate().skip(i) {
                image = image + value.scale(self.mass_factor.l(j, i));
            }
            norm = det::hypot(norm, image.abs());
        }
        if !norm.is_finite() || norm <= 0.0 {
            return Err(invalid("positive finite modal mass norm required"));
        }
        let (fraction, power) = quadratic_fraction(norm);
        let mut displacement: Vec<_> = vector
            .iter()
            .enumerate()
            .map(|(i, &z)| {
                scale_complex(divide(z, fraction), -self.coordinate_exponents[i] - power)
            })
            .collect::<Result<_, _>>()?;
        let pivot = (1..self.n).fold(0, |best, i| {
            if displacement[i].abs() > displacement[best].abs() {
                i
            } else {
                best
            }
        });
        let amplitude = displacement[pivot].abs();
        if !amplitude.is_finite() || amplitude <= 0.0 {
            return Err(invalid("finite physical displacement gauge required"));
        }
        let phase = divide(displacement[pivot].conj(), amplitude);
        for value in &mut displacement {
            *value = *value * phase;
        }
        displacement[pivot] = C64::from_re(amplitude);
        Ok(displacement)
    }

    fn original_residual(
        &self,
        s: C64,
        q: &[C64],
        v: &[C64],
        stress: &[C64],
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), MaxwellModeError> {
        let mut mechanical = 0.0_f64;
        for i in 0..self.n {
            poll(cancel)?;
            let mut row = SourceSum::default();
            for j in 0..self.n {
                row.add(s, v[j], &[self.source[0][i * self.n + j]])?;
                row.add(C64::ONE, v[j], &[self.source[1][i * self.n + j]])?;
                row.add(C64::ONE, q[j], &[self.source[2][i * self.n + j]])?;
            }
            for (branch, &sigma) in self.branches.iter().zip(stress) {
                poll(cancel)?;
                row.add(C64::ONE, sigma, &[branch.strain[i]])?;
            }
            mechanical = mechanical.max(row.relative()?);
        }
        let mut relaxation = 0.0_f64;
        for (branch, &sigma) in self.branches.iter().zip(stress) {
            poll(cancel)?;
            let mut row = SourceSum::default();
            row.add(C64::ONE, sigma, &[])?;
            row.add(s, sigma, &[branch.relaxation_time_s])?;
            for (&strain, &velocity) in branch.strain.iter().zip(v) {
                row.add(
                    C64::ONE,
                    velocity,
                    &[-branch.stiffness, branch.relaxation_time_s, strain],
                )?;
            }
            relaxation = relaxation.max(row.relative()?);
        }
        Ok((mechanical, relaxation))
    }
}

// Original-unit source products share one power-of-two row scale. DD sums
// retain cancellation without ever constructing dimensional m*s*v, k*tau*b*v
// or b*sigma products that may individually overflow. This is a backward
// residual, not a forward error or a relative weak-decay certificate.
#[derive(Default)]
struct SourceSum {
    re: Dd,
    im: Dd,
    envelope: f64,
    exponent: i32,
}
impl SourceSum {
    fn add(&mut self, left: C64, right: C64, factors: &[f64]) -> Result<(), MaxwellModeError> {
        if left == C64::ZERO || right == C64::ZERO || factors.contains(&0.0) {
            return Ok(());
        }
        let parts = |z: C64| -> Result<(C64, i32), MaxwellModeError> {
            let exponent = quadratic_fraction(z.re.abs().max(z.im.abs())).1;
            Ok((scale_complex(z, -exponent)?, exponent))
        };
        let (left, le) = parts(left)?;
        let (right, re) = parts(right)?;
        let mut value = left * right;
        let mut exponent = le + re;
        for &factor in factors {
            let (fraction, power) = quadratic_fraction(factor);
            value = value.scale(fraction.copysign(factor));
            exponent += power;
        }
        let common = if self.envelope == 0.0 {
            exponent
        } else {
            self.exponent.max(exponent)
        };
        let shift_dd = |value: Dd, power| {
            Dd::from_pair(
                quadratic_power_scale(value.hi, power),
                quadratic_power_scale(value.lo, power),
            )
        };
        self.re = shift_dd(self.re, self.exponent - common);
        self.im = shift_dd(self.im, self.exponent - common);
        self.envelope = quadratic_power_scale(self.envelope, self.exponent - common);
        value = C64::new(
            quadratic_power_scale(value.re, exponent - common),
            quadratic_power_scale(value.im, exponent - common),
        );
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
        self.envelope += value.abs();
        self.exponent = common;
        Ok(())
    }
    fn value(&self) -> Result<C64, MaxwellModeError> {
        scale_complex(C64::new(self.re.to_f64(), self.im.to_f64()), self.exponent)
    }
    fn relative(&self) -> Result<f64, MaxwellModeError> {
        let value = if self.envelope == 0.0 {
            0.0
        } else {
            det::hypot(self.re.to_f64(), self.im.to_f64()) / self.envelope
        };
        if !value.is_finite() {
            Err(invalid("finite original source-row residual required"))
        } else {
            Ok(value)
        }
    }
}

fn invalid(what: &'static str) -> MaxwellModeError {
    MaxwellModeError::Invalid(what)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn divide(z: C64, value: f64) -> C64 {
    C64::new(z.re / value, z.im / value)
}
fn poll(cancel: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), MaxwellModeError> {
    if cancel() {
        Err(MaxwellModeError::Cancelled)
    } else {
        Ok(())
    }
}
fn budget(quantity: &'static str, requested: u128, cap: u128) -> Result<(), MaxwellModeError> {
    if requested > cap {
        Err(MaxwellModeError::Budget {
            quantity,
            requested,
            cap,
        })
    } else {
        Ok(())
    }
}
fn scale(value: f64, power: i32) -> Result<f64, MaxwellModeError> {
    let result = quadratic_power_scale(value, power);
    if !result.is_finite() || (value != 0.0 && result == 0.0) {
        Err(invalid(
            "binary scaling erased or overflowed a nonzero value",
        ))
    } else {
        Ok(result)
    }
}
fn scale_complex(value: C64, power: i32) -> Result<C64, MaxwellModeError> {
    Ok(C64::new(scale(value.re, power)?, scale(value.im, power)?))
}
