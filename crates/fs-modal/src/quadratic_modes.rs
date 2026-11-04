//! Physical right mode shapes of the complete mechanical quadratic pencil.
//!
//! Refine `(s^2 M + s C + K) phi = 0` with the shared bordered Newton owner.
//! Exact binary congruence and time scaling keep mixed coordinate/rate units
//! out of its numerical conditioning. No companion eigenvectors, diagonal
//! damping approximation or discarded gyroscopic entries enter this path.
//!
//! These are local right eigenpairs, not a complete spectrum or certified
//! eigenvalue/shape enclosures. General force participation and pole residues
//! require left-mode information too; right modes alone cannot diagonalize a
//! nonproportionally damped forced response. [`PreparedQuadraticModes::prepare_response`]
//! obtains those residues through the shared full bordered factor, preserving
//! that left/right coupling. Repeated/defective roots can refuse.

use crate::{quadratic_fraction, quadratic_power_scale};
use fs_la::factor::{Cholesky, cholesky};
use fs_la::nonlinear_eig::{
    MAX_NONLINEAR_POLE_DIMENSION, NonlinearEigError, PoleEvaluation, newton_pole_with_cancel,
};
use fs_math::{c64::C64, det};

/// Shared local-pole accuracy and resource options; reference rate is physical.
pub use fs_la::nonlinear_eig::PoleOptions;

mod response;
pub use response::{
    PoleResponseOptions, PreparedQuadraticResponse, PreparedQuadraticSensitivities,
    QuadraticForceResidue, QuadraticModeDerivative, QuadraticObservation,
    QuadraticParameterDirection, QuadraticResponseError, QuadraticResponseOptions,
};

/// One locally resolved right mode in the ORIGINAL generalized coordinates.
#[derive(Clone, Debug)]
pub struct QuadraticMode {
    /// Laplace exponent [1/s]; positive real parts remain physically unstable.
    pub laplace_frequency_per_s: C64,
    /// Complex displacement shape, with `phi^H M phi = 1`. The largest
    /// original-coordinate component is real and positive (lowest-index tie).
    pub displacement: Vec<C64>,
    /// Complex velocity shape `s * phi`, without frequency/sign modification.
    pub velocity: Vec<C64>,
    /// Maximum original-row residual after exact binary equilibration, divided
    /// by that row's coefficient envelope times the full vector infinity norm.
    /// The envelope uses max(pole magnitude, reference rate), resolving zero poles.
    /// Recomputed from source M/C/K and the PUBLISHED displacement and pole.
    pub relative_residual: f64,
    /// Shared Newton owner's residual before physical recovery.
    pub solver_relative_residual: f64,
    /// Last pole correction divided by max(reference rate, pole magnitude).
    pub relative_correction: f64,
    /// Number of original-pencil Newton evaluations.
    pub iterations: usize,
}

/// No partial mode or changed source coefficients accompany a refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum QuadraticModeError {
    /// Invalid shape, finite domain, symmetry or representability.
    Invalid(&'static str),
    /// The original mass failed symmetric positive-definite admission.
    MassNotPositiveDefinite,
    /// Explicit numerical payload or arithmetic budget was insufficient.
    Budget {
        /// Explicit resource being admitted.
        quantity: &'static str,
        /// Conservative requested resource count.
        requested: u128,
        /// Caller-authored resource ceiling.
        cap: u128,
    },
    /// Shared local-pole solve refusal, including unresolved multiplicity.
    Solve(NonlinearEigError),
    /// The recovered physical eigenpair failed its independent source-row gate.
    Residual {
        /// Measured original-row backward residual.
        measured: f64,
        /// Caller-authored relative residual ceiling.
        limit: f64,
    },
    /// Cancellation observed before publication.
    Cancelled,
}

impl core::fmt::Display for QuadraticModeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-COMPLEX-MODE: {self:?}")
    }
}
impl std::error::Error for QuadraticModeError {}
impl From<NonlinearEigError> for QuadraticModeError {
    fn from(error: NonlinearEigError) -> Self {
        if error == NonlinearEigError::Cancelled {
            Self::Cancelled
        } else {
            Self::Solve(error)
        }
    }
}

/// Immutable source pencil and cached mass factor for repeated local modes.
///
/// Borrowed source arrays remain the authority. M is exactly symmetric and
/// numerically SPD. C and K may be nonsymmetric or indefinite: active damping,
/// gyroscopic coupling and unstable stiffness are retained, never projected.
#[derive(Debug)]
pub struct PreparedQuadraticModes<'a> {
    n: usize,
    source: [&'a [f64]; 3],
    scaled: [Vec<f64>; 3],
    mass_factor: Cholesky,
    coordinate_exponents: Vec<i32>,
    time_exponent: i32,
    coefficient_norms: [f64; 3],
    options: PoleOptions,
}

impl<'a> PreparedQuadraticModes<'a> {
    /// Admit the complete real row-major M/C/K once, without inverting M.
    /// `options.spectral_scale` is a positive physical reference rate [1/s].
    /// The residual and correction tolerances retain the shared owner's meaning.
    ///
    /// Here `max_dense_bytes` includes this owner's cached/staged numerical
    /// arrays and the shared solver's payload; caller inputs, previously
    /// returned modes, allocator metadata and process stacks are excluded.
    /// `max_work` conservatively admits preparation and one complete refinement;
    /// every later refinement is bounded identically. Dense Cholesky is one
    /// bounded indivisible stage; other dense traversals poll once per row.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        mass: &'a [f64],
        damping: &'a [f64],
        stiffness: &'a [f64],
        n: usize,
        mut options: PoleOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, QuadraticModeError> {
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
                "positive bounded dimension, reference rate and tolerances required",
            ));
        }
        let size = n as u128;
        let d = size + 1;
        // Includes complete callback matrices, physical recovery and the
        // original Cholesky's trailing update plus serial GEMM panel scratch.
        let additional_bytes =
            64 * size * size + 512 * size + if n > 32 { 2 * 1024 * 1024 } else { 0 };
        let solver_bytes = 128 * d * d + 256 * d;
        budget(
            "numerical payload bytes",
            additional_bytes + solver_bytes,
            options.max_dense_bytes as u128,
        )?;
        let additional_work =
            8 * size * size * size + 64 * size * size * (options.max_iterations as u128 + 1);
        let solver_work = d * d * (d + 16) * options.max_iterations as u128;
        budget(
            "arithmetic work",
            additional_work + solver_work,
            options.max_work as u128,
        )?;
        options.max_dense_bytes -= additional_bytes as usize;
        options.max_work -= additional_work as u64;
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
                return Err(QuadraticModeError::MassNotPositiveDefinite);
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
        poll(&mut cancelled)?;
        let mass_factor =
            cholesky(&scaled[0], n).map_err(|_| QuadraticModeError::MassNotPositiveDefinite)?;
        poll(&mut cancelled)?;
        Ok(Self {
            n,
            source,
            scaled,
            mass_factor,
            coordinate_exponents,
            time_exponent,
            coefficient_norms,
            options,
        })
    }

    /// Refine a nearby physical Laplace pole and nonzero physical displacement
    /// seed. Complex phase and relative motion are retained. The seed's norm
    /// carries no physical amplitude; the returned mode has unit modal mass.
    ///
    /// The shared owner checks both pole correction and equation residual,
    /// including at exact roots. A singular normalized Jacobian refuses
    /// repeated or defective poles rather than publishing an arbitrary shape.
    /// This local method does not promise convergence from an arbitrary seed.
    pub fn refine(
        &self,
        laplace_seed_per_s: C64,
        displacement_seed: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<QuadraticMode, QuadraticModeError> {
        poll(&mut cancelled)?;
        let seed = self.equilibrate_vector(displacement_seed)?;
        let value = scale_complex(laplace_seed_per_s, -self.time_exponent)?;
        let pole = newton_pole_with_cancel(
            self.n,
            value,
            &seed,
            &self.options,
            |s, cancel| self.evaluate(s, cancel),
            &mut cancelled,
        )?;
        let displacement = self.physical_displacement(&pole.vector, &mut cancelled)?;
        let laplace_frequency_per_s = scale_complex(pole.value, self.time_exponent)?;
        let velocity: Vec<_> = displacement
            .iter()
            .map(|&value| laplace_frequency_per_s * value)
            .collect();
        if displacement.iter().chain(&velocity).any(|&z| !finite(z))
            || velocity.iter().zip(&displacement).any(|(&v, &q)| {
                laplace_frequency_per_s != C64::ZERO && q != C64::ZERO && v == C64::ZERO
            })
        {
            return Err(invalid("physical mode and velocity must be representable"));
        }
        let relative_residual =
            self.original_residual(laplace_frequency_per_s, &displacement, &mut cancelled)?;
        if relative_residual > self.options.residual_tolerance {
            return Err(QuadraticModeError::Residual {
                measured: relative_residual,
                limit: self.options.residual_tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(QuadraticMode {
            laplace_frequency_per_s,
            displacement,
            velocity,
            relative_residual,
            solver_relative_residual: pole.residual,
            relative_correction: pole.relative_step,
            iterations: pole.iterations,
        })
    }

    // Canonical physical gauge shared by mode discovery and response preparation.
    fn physical_displacement(
        &self,
        vector: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, QuadraticModeError> {
        let norm = self.mass_norm(vector, cancelled)?;
        let (mantissa, exponent) = quadratic_fraction(norm);
        let mut displacement = Vec::with_capacity(self.n);
        for (i, &value) in vector.iter().enumerate() {
            poll(cancelled)?;
            displacement.push(scale_complex(
                divide(value, mantissa),
                -self.coordinate_exponents[i] - exponent,
            )?);
        }
        let pivot = (1..self.n).fold(0, |best, i| {
            if displacement[i].abs() > displacement[best].abs() {
                i
            } else {
                best
            }
        });
        let amplitude = displacement[pivot].abs();
        let phase = divide(displacement[pivot].conj(), amplitude);
        for value in &mut displacement {
            *value = *value * phase;
        }
        displacement[pivot] = C64::from_re(amplitude);
        Ok(displacement)
    }

    fn evaluate(
        &self,
        s: C64,
        cancel: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, NonlinearEigError> {
        let mut matrix = vec![C64::ZERO; self.n * self.n];
        let mut derivative = vec![C64::ZERO; self.n * self.n];
        for i in 0..self.n {
            if cancel() {
                return Err(NonlinearEigError::Cancelled);
            }
            for j in 0..self.n {
                let index = i * self.n + j;
                let [m, c, k] = self.scaled.each_ref().map(|a| a[index]);
                matrix[index] = s * (s.scale(m) + C64::from_re(c)) + C64::from_re(k);
                derivative[index] = s.scale(m).scale(2.0) + C64::from_re(c);
            }
        }
        let radius = s.abs().max(self.options.spectral_scale);
        let [m, c, k] = self.coefficient_norms;
        let residual_scale = (m * radius + c) * radius + k;
        if !residual_scale.is_finite()
            || residual_scale <= 0.0
            || matrix.iter().chain(&derivative).any(|&z| !finite(z))
        {
            return Err(NonlinearEigError::FunctionRefused { at: s });
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale,
        })
    }

    // A common power of two normalizes the seed BEFORE coordinate scaling;
    // no intermediate may overflow just because its arbitrary amplitude does.
    fn equilibrate_vector(&self, x: &[C64]) -> Result<Vec<C64>, QuadraticModeError> {
        if x.len() != self.n || x.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "one finite complex displacement per original coordinate required",
            ));
        }
        let mut largest = i32::MIN;
        for (i, &z) in x.iter().enumerate() {
            for value in [z.re, z.im] {
                if value != 0.0 {
                    largest =
                        largest.max(quadratic_fraction(value).1 + self.coordinate_exponents[i]);
                }
            }
        }
        if largest == i32::MIN {
            return Err(invalid("nonzero displacement seed required"));
        }
        x.iter()
            .enumerate()
            .map(|(i, &z)| scale_complex(z, self.coordinate_exponents[i] - largest))
            .collect()
    }

    fn mass_norm(
        &self,
        y: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, QuadraticModeError> {
        // phi^H M phi = y^H (D M D) y = ||L^T y||^2. Computing
        // the norm of the Cholesky image avoids cancellation/negative energy.
        let mut norm = 0.0;
        for i in 0..self.n {
            poll(cancelled)?;
            let mut image = C64::ZERO;
            for (j, &value) in y.iter().enumerate().skip(i) {
                image = image + value.scale(self.mass_factor.l(j, i));
            }
            norm = det::hypot(norm, image.abs());
        }
        if !norm.is_finite() || norm <= 0.0 {
            return Err(invalid("positive finite physical mass norm required"));
        }
        Ok(norm)
    }

    fn original_residual(
        &self,
        physical_s: C64,
        physical_q: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, QuadraticModeError> {
        let y = self.equilibrate_vector(physical_q)?;
        let s = scale_complex(physical_s, -self.time_exponent)?;
        let magnitude = s.abs().max(self.options.spectral_scale);
        let vector_norm = y.iter().map(|z| z.abs()).fold(0.0, f64::max);
        let mut maximum = 0.0_f64;
        for i in 0..self.n {
            poll(cancelled)?;
            let mut actions = [C64::ZERO; 3];
            let mut envelopes = [0.0; 3];
            for (j, &value) in y.iter().enumerate() {
                for order in 0..3 {
                    let a = scale(
                        self.source[order][i * self.n + j],
                        -self.coordinate_exponents[i]
                            - self.coordinate_exponents[j]
                            - order as i32 * self.time_exponent,
                    )?;
                    actions[order] = actions[order] + value.scale(a);
                    envelopes[order] += a.abs();
                }
            }
            // Independent split source contractions, not the Newton matrix.
            let residual = (s * (s * actions[0] + actions[1]) + actions[2]).abs();
            let envelope = ((envelopes[0] * magnitude + envelopes[1]) * magnitude + envelopes[2])
                * vector_norm;
            let relative = if envelope == 0.0 && residual == 0.0 {
                0.0
            } else {
                residual / envelope
            };
            if !relative.is_finite() {
                return Err(invalid(
                    "physical source-row residual must be representable",
                ));
            }
            maximum = maximum.max(relative);
        }
        Ok(maximum)
    }
}

fn invalid(what: &'static str) -> QuadraticModeError {
    QuadraticModeError::Invalid(what)
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), QuadraticModeError> {
    if cancelled() {
        Err(QuadraticModeError::Cancelled)
    } else {
        Ok(())
    }
}
fn budget(quantity: &'static str, requested: u128, cap: u128) -> Result<(), QuadraticModeError> {
    if requested > cap {
        Err(QuadraticModeError::Budget {
            quantity,
            requested,
            cap,
        })
    } else {
        Ok(())
    }
}
fn scale(value: f64, exponent: i32) -> Result<f64, QuadraticModeError> {
    let result = quadratic_power_scale(value, exponent);
    if !result.is_finite() || value != 0.0 && result == 0.0 {
        Err(invalid(
            "binary coordinate/rate scaling erased or overflowed a nonzero value",
        ))
    } else {
        Ok(result)
    }
}
fn scale_complex(value: C64, exponent: i32) -> Result<C64, QuadraticModeError> {
    Ok(C64::new(
        scale(value.re, exponent)?,
        scale(value.im, exponent)?,
    ))
}
