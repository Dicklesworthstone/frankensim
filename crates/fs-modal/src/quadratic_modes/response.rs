//! Physical Laurent residues of complete mechanical quadratic pencils.
//!
//! With D = diag(2^-e_i), a = 2^t, s = a*mu and q = D*y, the prepared
//! polynomial is T_hat(mu) = D*T(a*mu)*D/a^2. The physical residue is therefore
//! D*Res_mu(T_hat^-1 * D*F)/a. Keeping that rate factor is essential: the
//! internal root variable, displacement coordinates and modal mass all differ.
//!
//! fs-la's shared bordered factor obtains the left/right participation without
//! forming an adjoint eigenvector or assuming symmetry/proportional damping.

use super::{
    C64, PreparedQuadraticModes, QuadraticMode, QuadraticModeError, divide, finite,
    quadratic_fraction, scale, scale_complex,
};
use fs_la::pole_response::{PoleResponseError, PreparedPole};

mod sensitivity;
pub use sensitivity::{
    PreparedQuadraticSensitivities, QuadraticModeDerivative, QuadraticParameterDirection,
};

/// Shared numerical factor and batch limits. The reference rate is physical [1/s].
pub use fs_la::pole_response::PoleResponseOptions;

/// Independent limits for numerical factorization and physical mapping.
/// Each limit applies to preparation and independently to EACH complete batch.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticResponseOptions {
    /// Shared bordered-factor accuracy and resource limits.
    pub numerical: PoleResponseOptions,
    /// Complete observation count per batch, in 0..=256; zero requests states only.
    pub max_observations: usize,
    /// Conservative live physical payload: evaluated pencil, mapped forces,
    /// physical shapes and largest result batch. Excludes the separately
    /// admitted fs-la factor, original source model, caller inputs and allocator.
    pub max_physical_bytes: usize,
    /// Admission units 64*(k+1)*n^2 + 64*k*(m+1)*n + 256*n, where k and m
    /// are the declared force and observation ceilings; not measured FLOPs.
    pub max_physical_work: u64,
}

/// A refused batch publishes no partial result and leaves its factor reusable.
#[derive(Clone, Debug, PartialEq)]
pub enum QuadraticResponseError {
    /// Invalid dimensions, inputs or declared limits.
    Invalid(&'static str),
    /// Original-coordinate normalization, residual, budget or scaling refusal.
    Physical(QuadraticModeError),
    /// Shared bordered-factor refusal, including a singular pole linearization.
    Numerical(PoleResponseError),
    /// Published physical parameter derivatives failed an independent source
    /// polynomial or differentiated unit-mass constraint check.
    SensitivityResidual {
        /// Maximum original-polynomial row backward residual.
        polynomial: f64,
        /// Differentiated physical unit-mass constraint residual.
        normalization: f64,
        /// Caller-authored linear residual ceiling.
        limit: f64,
    },
    /// Cancellation observed before publication.
    Cancelled,
}
impl core::fmt::Display for QuadraticResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-QUADRATIC-RESPONSE: {self:?}")
    }
}
impl std::error::Error for QuadraticResponseError {}
impl From<QuadraticModeError> for QuadraticResponseError {
    fn from(error: QuadraticModeError) -> Self {
        match error {
            QuadraticModeError::Cancelled => Self::Cancelled,
            error => Self::Physical(error),
        }
    }
}
impl From<PoleResponseError> for QuadraticResponseError {
    fn from(error: PoleResponseError) -> Self {
        match error {
            PoleResponseError::Cancelled => Self::Cancelled,
            error => Self::Numerical(error),
        }
    }
}

/// Linear physical observation y = displacement^T*q + velocity^T*v.
/// Weights are in the ORIGINAL coordinates. This is a bilinear contraction;
/// no implicit complex conjugation or energy/power interpretation is applied.
/// Supply both complete vectors, using zeros for an unused term. For analytic
/// frequency-dependent sensors, these are their weights evaluated at the pole.
#[derive(Clone, Copy, Debug)]
pub struct QuadraticObservation<'a> {
    /// Weights multiplying the original displacement coordinates.
    pub displacement: &'a [C64],
    /// Weights multiplying the original velocity coordinates.
    pub velocity: &'a [C64],
}

/// One pole's physical contribution for an analytic force evaluated at s0.
/// Q(s) = displacement/(s-s0) + regular, with zero-initial-state V(s)=s*Q(s).
/// These are Laplace residues, not driven harmonic amplitudes or powers.
/// For exp(-i*omega*t), omega0=i*s0 and the omega-residues are i times
/// these coefficients, since s-s0=-i*(omega-omega0).
#[derive(Clone, Debug)]
pub struct QuadraticForceResidue {
    /// Original physical Laplace exponent [1/s], including active growth.
    pub laplace_pole_per_s: C64,
    /// Laurent residue in each original generalized displacement coordinate.
    pub displacement: Vec<C64>,
    /// s0 times the displacement residue. A simple zero pole can cancel here.
    pub velocity: Vec<C64>,
    /// Residue of every supplied linear observation, in its input order.
    pub observations: Vec<C64>,
    /// displacement = prepared.right_displacement()*participation, in the
    /// canonical ORIGINAL-coordinate phase gauge with phi^H*M*phi=1.
    /// This includes the true left-mode force coupling; it is not phi^H*force.
    pub participation: C64,
    /// Shared owner's maximum backward residual of the original bordered rows.
    pub linear_residual: f64,
}

/// Immutable factor for one isolated simple pole of the original full M/C/K.
/// Preparation snapshots the pencil; subsequent forces and sensors need no
/// refactorization. One pole cannot reconstruct the full transfer function or
/// certify observability, stability, spectral completeness or multiplicity.
#[derive(Debug)]
pub struct PreparedQuadraticResponse {
    numerical: PreparedPole,
    coordinate_exponents: Vec<i32>,
    time_exponent: i32,
    pole: C64,
    displacement: Vec<C64>,
    velocity: Vec<C64>,
    participation_direction: C64,
    participation_exponent: i32,
    physical_residual: f64,
    options: QuadraticResponseOptions,
}

fn invalid(what: &'static str) -> QuadraticResponseError {
    QuadraticResponseError::Invalid(what)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), QuadraticResponseError> {
    if cancelled() {
        Err(QuadraticResponseError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(n: usize, o: &QuadraticResponseOptions) -> Result<(), QuadraticResponseError> {
    if n > o.numerical.max_dimension
        || !(1..=256).contains(&o.numerical.max_vectors)
        || o.max_observations > 256
        || !o.numerical.spectral_scale.is_finite()
        || o.numerical.spectral_scale <= 0.0
        || [
            o.numerical.pole_residual_tolerance,
            o.numerical.pole_step_tolerance,
            o.numerical.linear_residual_tolerance,
        ]
        .iter()
        .any(|t| !t.is_finite() || *t <= 0.0 || *t >= 1.0)
    {
        return Err(invalid(
            "complete dimension, bounded batches, positive reference rate and tolerances required",
        ));
    }
    let (n, k, m) = (
        n as u128,
        o.numerical.max_vectors as u128,
        o.max_observations as u128,
    );
    super::budget(
        "physical response bytes",
        16 * (2 * n * n + k * (5 * n + m + 2) + 16 * n + 2 * m) + 256 * (k + m),
        o.max_physical_bytes as u128,
    )?;
    super::budget(
        "physical response work",
        64 * (k + 1) * n * n + 64 * k * (m + 1) * n + 256 * n,
        o.max_physical_work as u128,
    )?;
    Ok(())
}

impl PreparedQuadraticModes<'_> {
    /// Recheck an already resolved physical pole against THIS original pencil,
    /// then factor its full bordered linearization once. Cached mode receipts,
    /// velocities and normalization are not trusted. The returned physical
    /// shape is independently mass-normalized in the canonical phase gauge.
    /// Nonsymmetric/active C/K and signed gyroscopic terms are retained. Exact
    /// repeated/defective examples refuse through the shared singular border.
    pub fn prepare_response(
        &self,
        mode: &QuadraticMode,
        options: QuadraticResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedQuadraticResponse, QuadraticResponseError> {
        poll(&mut cancelled)?;
        admit(self.n, &options)?;
        let vector = self.equilibrate_vector(&mode.displacement)?;
        let value = scale_complex(mode.laplace_frequency_per_s, -self.time_exponent)?;
        let evaluation = self
            .evaluate(value, &mut cancelled)
            .map_err(QuadraticModeError::from)?;
        let mut numerical_options = options.numerical;
        numerical_options.spectral_scale =
            scale(numerical_options.spectral_scale, -self.time_exponent)?;
        let numerical = PreparedPole::new(
            value,
            &vector,
            &evaluation,
            numerical_options,
            &mut cancelled,
        )?;
        let displacement = self.physical_displacement(numerical.right_vector(), &mut cancelled)?;
        let pole = scale_complex(numerical.value(), self.time_exponent)?;
        let velocity = displacement
            .iter()
            .map(|&q| product(pole, q, 0))
            .collect::<Result<Vec<_>, _>>()?;
        let physical_residual = self.original_residual(pole, &displacement, &mut cancelled)?;
        if physical_residual > options.numerical.pole_residual_tolerance {
            return Err(QuadraticModeError::Residual {
                measured: physical_residual,
                limit: options.numerical.pole_residual_tolerance,
            }
            .into());
        }
        let pivot = (1..self.n).fold(0, |best, i| {
            if displacement[i].abs() > displacement[best].abs() {
                i
            } else {
                best
            }
        });
        // D*v_hat = phi * (v_hat[p] * 2^-e_p / phi[p]). The physical
        // phase gauge makes phi[p] positive real. Keep its power separate so
        // neither enormous coordinate units nor tiny unit-mass shapes overflow.
        let (mantissa, exponent) = quadratic_fraction(displacement[pivot].re);
        let participation_direction = divide(numerical.right_vector()[pivot], mantissa);
        let participation_exponent =
            -self.time_exponent - self.coordinate_exponents[pivot] - exponent;
        poll(&mut cancelled)?;
        Ok(PreparedQuadraticResponse {
            numerical,
            coordinate_exponents: self.coordinate_exponents.clone(),
            time_exponent: self.time_exponent,
            pole,
            displacement,
            velocity,
            participation_direction,
            participation_exponent,
            physical_residual,
            options,
        })
    }
}

impl PreparedQuadraticResponse {
    /// Prepared original physical Laplace pole [1/s].
    pub fn laplace_pole_per_s(&self) -> C64 {
        self.pole
    }
    /// Canonical physical right shape with phi^H*M*phi=1.
    pub fn right_displacement(&self) -> &[C64] {
        &self.displacement
    }
    /// Original physical right velocity shape s0*phi.
    pub fn right_velocity(&self) -> &[C64] {
        &self.velocity
    }
    /// Independent residual recomputed from the original source M/C/K.
    pub fn physical_residual(&self) -> f64 {
        self.physical_residual
    }
    /// Rechecked pole correction relative to the declared reference rate.
    pub fn relative_correction(&self) -> f64 {
        self.numerical.relative_correction()
    }

    /// Reuse the factor for a complete batch of physical generalized forces.
    /// For arbitrary analytic F(s), supply F(s0); its derivative contributes
    /// only to regular terms. Each solve costs O(n^2), with no new eigenproblem
    /// or adjoint factor. Observation contractions cost O(n) per force/sensor.
    /// Zero force returns zero. A force may be orthogonal to the right shape
    /// and still excite a nonnormal mode; the full bordered solve decides.
    pub fn forcing_residues(
        &self,
        forces: &[&[C64]],
        observations: &[QuadraticObservation<'_>],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<QuadraticForceResidue>, QuadraticResponseError> {
        poll(&mut cancelled)?;
        let n = self.displacement.len();
        if forces.is_empty()
            || forces.len() > self.options.numerical.max_vectors
            || observations.len() > self.options.max_observations
            || forces
                .iter()
                .any(|x| x.len() != n || x.iter().any(|&z| !finite(z)))
            || observations.iter().any(|o| {
                o.displacement.len() != n
                    || o.velocity.len() != n
                    || o.displacement.iter().chain(o.velocity).any(|&z| !finite(z))
            })
        {
            return Err(invalid(
                "complete finite force and observation batches within declared ceilings required",
            ));
        }
        let mut mapped = Vec::with_capacity(forces.len());
        let mut exponents = Vec::with_capacity(forces.len());
        for &force in forces {
            poll(&mut cancelled)?;
            // Map D*F, not D*F/a^2. Normalize its common amplitude BEFORE
            // multiplying coordinate scales; recover both amplitude and 1/a
            // only after the solve. No physical source column is truncated.
            let exponent = force
                .iter()
                .enumerate()
                .filter_map(|(i, &z)| {
                    component_exponent(z).map(|e| e - self.coordinate_exponents[i])
                })
                .max()
                .unwrap_or(0);
            let column = force
                .iter()
                .enumerate()
                .map(|(i, &z)| scale_complex(z, -self.coordinate_exponents[i] - exponent))
                .collect::<Result<Vec<_>, _>>()?;
            exponents.push(exponent);
            mapped.push(column);
        }
        let refs: Vec<_> = mapped.iter().map(Vec::as_slice).collect();
        let residues = self.numerical.forcing_residues(&refs, &mut cancelled)?;
        let mut out = Vec::with_capacity(residues.len());
        for (r, exponent) in residues.into_iter().zip(exponents) {
            poll(&mut cancelled)?;
            let participation = product(
                r.coefficient,
                self.participation_direction,
                exponent + self.participation_exponent,
            )?;
            let mut displacement = Vec::with_capacity(n);
            let mut velocity = Vec::with_capacity(n);
            for i in 0..n {
                poll(&mut cancelled)?;
                // Recover from the factor's coefficient and vector separately:
                // multiplying them prematurely can erase a tiny component
                // that is finite after conversion to its physical coordinates.
                let q = product(
                    r.coefficient,
                    self.numerical.right_vector()[i],
                    exponent - self.time_exponent - self.coordinate_exponents[i],
                )?;
                displacement.push(q);
                velocity.push(product(self.pole, q, 0)?);
            }
            let mut outputs = Vec::with_capacity(observations.len());
            for observation in observations {
                poll(&mut cancelled)?;
                outputs.push(observe(
                    observation,
                    &displacement,
                    &velocity,
                    &mut cancelled,
                )?);
            }
            out.push(QuadraticForceResidue {
                laplace_pole_per_s: self.pole,
                displacement,
                velocity,
                observations: outputs,
                participation,
                linear_residual: r.linear_residual,
            });
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

fn component_exponent(z: C64) -> Option<i32> {
    [z.re, z.im]
        .into_iter()
        .filter(|x| *x != 0.0)
        .map(|x| quadratic_fraction(x).1)
        .max()
}

// Multiply with a deferred binary scale. A representable physical product
// must not overflow/underflow simply because its factors use different units.
fn product(a: C64, b: C64, exponent: i32) -> Result<C64, QuadraticResponseError> {
    let (Some(ea), Some(eb)) = (component_exponent(a), component_exponent(b)) else {
        return Ok(C64::ZERO);
    };
    let a = scale_complex(a, -ea)?;
    let b = scale_complex(b, -eb)?;
    Ok(scale_complex(a * b, exponent + ea + eb)?)
}

fn observe(
    observation: &QuadraticObservation<'_>,
    q: &[C64],
    v: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<C64, QuadraticResponseError> {
    let terms = || {
        observation
            .displacement
            .iter()
            .zip(q)
            .chain(observation.velocity.iter().zip(v))
    };
    let largest = terms()
        .filter_map(|(&a, &b)| Some(component_exponent(a)? + component_exponent(b)?))
        .max();
    let Some(largest) = largest else {
        return Ok(C64::ZERO);
    };
    let (mut sum, mut correction) = (C64::ZERO, C64::ZERO);
    for (&a, &b) in terms() {
        poll(cancelled)?;
        let term = product(a, b, -largest)?;
        let adjusted = term - correction;
        let next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    Ok(scale_complex(sum, largest)?)
}
