//! Linear rotor-bearing motion at a declared steady operating speed.
//!
//! Original generalized coordinates satisfy `M q'' + (C + Omega G) q' + K q = f`.
//! Mass is SPD and G is skew; complete bearing C/K may be nonsymmetric or active.
//! Coordinates may mix translations and angles, with conjugate forces/torques
//! such that `f dot q'` is power. All matrices must use those same coordinates.
//! The phasor convention is `q(t) = Re(q_hat exp(+i omega t))`.
//!
//! Native fs-modal owns the scaled quadratic poles; native fs-la owns complex
//! factorization and solve. These are numerical operating-point estimates,
//! without certified crossings, mode identities, continuum stability or peak
//! response claims. Changing bearings, preload or spin requires the caller's
//! actual new coefficients; a fixed-matrix speed sweep does not infer them.

use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_la::factor::cholesky;
use fs_math::c64::C64;
use fs_modal::{ModalError, quadratic_eigenvalues_with_cancel};

/// Physical rigid-disk, bearing and imbalance assembly in original coordinates.
pub mod rigid;

/// Small dense operating-point bound; the quadratic companion has twice this size.
pub const MAX_ROTOR_DOFS: usize = 32;
/// Maximum operating points in one cancellable fixed-coefficient sweep.
pub const MAX_CAMPBELL_SAMPLES: usize = 256;
/// Native complex QR sweep bound for each active eigenvalue window.
pub const ROTOR_QR_SWEEPS_PER_WINDOW: usize = 120;

/// No partial physical response or sample series is returned on refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum RotorError {
    /// Dimension, coefficient structure, load or representability admission failed.
    Invalid(&'static str),
    /// Original full mass is not numerically positive definite.
    MassNotSpd,
    /// Cancellation was observed before publication.
    Cancelled,
    /// Original quadratic-pencil solver refused.
    Modal(ModalError),
    /// Original complex factorization/solve refused; resonance is never floored.
    Linear(LuComplexError),
    /// Original dynamic equilibrium did not meet the caller's residual limit.
    Residual {
        /// Measured scale-relative original-force residual.
        relative: f64,
        /// Caller-supplied maximum.
        limit: f64,
    },
}

impl core::fmt::Display for RotorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rotor-bearing analysis refused: {self:?}")
    }
}
impl std::error::Error for RotorError {}
impl From<ModalError> for RotorError {
    fn from(error: ModalError) -> Self {
        if error == ModalError::Cancelled {
            Self::Cancelled
        } else {
            Self::Modal(error)
        }
    }
}
impl From<LuComplexError> for RotorError {
    fn from(error: LuComplexError) -> Self {
        if error == LuComplexError::Cancelled {
            Self::Cancelled
        } else {
            Self::Linear(error)
        }
    }
}

/// Immutable full matrices on one original generalized-coordinate space.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorBearingModel {
    n: usize,
    mass: Vec<f64>,
    damping: Vec<f64>,
    stiffness: Vec<f64>,
    gyro: Vec<f64>,
}

/// Fixed coefficients and assembled velocity operator at one signed spin rate.
#[derive(Debug, Clone)]
pub struct RotorOperatingPoint<'a> {
    model: &'a RotorBearingModel,
    spin_rad_s: f64,
    velocity_operator: Vec<f64>,
}

/// All complex poles of one original quadratic pencil, without mode identities.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorPoleReport {
    /// Signed steady shaft speed used in `C + Omega G` [rad/s].
    pub spin_rad_s: f64,
    /// All 2n poles [1/s], in native canonical (real, imaginary) order.
    /// Positive real parts are retained; no stability threshold is inferred.
    pub poles_per_s: Vec<C64>,
}

/// Forced complex response in the caller's original coordinates and units.
#[derive(Debug, Clone, PartialEq)]
pub struct RotorHarmonicResponse {
    /// Signed steady shaft rate [rad/s].
    pub spin_rad_s: f64,
    /// Signed forcing frequency in `exp(+i omega t)` [rad/s].
    pub frequency_rad_s: f64,
    /// Displacement/angle phasors in the original generalized coordinates.
    pub displacement: Vec<C64>,
    /// Velocity phasors `i omega q_hat`.
    pub velocity: Vec<C64>,
    /// Original equilibrium residual divided by the sum of force-component norms.
    /// This is not a condition-number or response-error bound.
    pub relative_residual: f64,
    /// Cycle-averaged supplied power, `Re(f_hat^H v_hat)/2` [W].
    pub mean_input_power_w: f64,
    /// Signed work rate against the original C term [W]. It may be negative;
    /// nonsymmetric/active bearings are not silently classified as passive heat.
    pub mean_velocity_force_power_w: f64,
    /// Signed work rate against K q [W], zero for symmetric conservative K
    /// in exact arithmetic; cross-coupled stiffness can exchange net work.
    pub mean_stiffness_force_power_w: f64,
    /// Measured gyroscopic work rate [W], zero in exact arithmetic for skew G.
    pub mean_gyroscopic_power_w: f64,
    /// Measured inertial work rate [W], zero for symmetric M in exact arithmetic.
    pub mean_inertial_power_w: f64,
    /// Supplied minus the four independently evaluated component powers [W].
    pub mean_power_residual_w: f64,
}

impl RotorBearingModel {
    /// Admit full finite row-major n-by-n coefficients. M must be exactly
    /// symmetric and numerically SPD; G must be exactly skew with zero diagonal.
    /// C/K are retained unchanged. Nothing is diagonalized, lumped or symmetrized.
    /// Mass admission uses a diagonally scaled native Cholesky, then stores the
    /// original tensor. Dense work is bounded by [`MAX_ROTOR_DOFS`].
    pub fn new(
        dofs: usize,
        mass: &[f64],
        damping: &[f64],
        stiffness: &[f64],
        gyro_per_rad_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RotorError> {
        poll(&mut cancelled)?;
        if !(1..=MAX_ROTOR_DOFS).contains(&dofs)
            || [mass, damping, stiffness, gyro_per_rad_s]
                .iter()
                .any(|a| a.len() != dofs * dofs || a.iter().any(|x| !x.is_finite()))
        {
            return Err(RotorError::Invalid(
                "matching finite bounded full matrices required",
            ));
        }
        let mut scale = Vec::with_capacity(dofs);
        for i in 0..dofs {
            poll(&mut cancelled)?;
            if mass[i * dofs + i] <= 0.0 {
                return Err(RotorError::MassNotSpd);
            }
            scale.push(fs_math::det::sqrt(mass[i * dofs + i]));
            for j in 0..dofs {
                if mass[i * dofs + j] != mass[j * dofs + i] {
                    return Err(RotorError::Invalid(
                        "mass must be symmetric in the declared coordinates",
                    ));
                }
                if gyro_per_rad_s[i * dofs + j] != -gyro_per_rad_s[j * dofs + i] {
                    return Err(RotorError::Invalid(
                        "gyroscopic coefficient must be skew-symmetric",
                    ));
                }
            }
        }
        let mut scaled = vec![0.0; dofs * dofs];
        for i in 0..dofs {
            poll(&mut cancelled)?;
            for j in 0..dofs {
                scaled[i * dofs + j] = (mass[i * dofs + j] / scale[i]) / scale[j];
            }
        }
        if scaled.iter().any(|x| !x.is_finite()) || cholesky(&scaled, dofs).is_err() {
            return Err(RotorError::MassNotSpd);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            n: dofs,
            mass: mass.to_vec(),
            damping: damping.to_vec(),
            stiffness: stiffness.to_vec(),
            gyro: gyro_per_rad_s.to_vec(),
        })
    }

    /// Original unconstrained generalized-coordinate dimension.
    pub const fn dofs(&self) -> usize {
        self.n
    }
    /// Original full row-major mass/inertia tensor.
    pub fn mass(&self) -> &[f64] {
        &self.mass
    }
    /// Original full row-major bearing velocity-force coefficients.
    pub fn damping(&self) -> &[f64] {
        &self.damping
    }
    /// Original full row-major tangent stiffness.
    pub fn stiffness(&self) -> &[f64] {
        &self.stiffness
    }
    /// Original skew gyroscopic matrix multiplying signed spin rate.
    pub fn gyro_per_rad_s(&self) -> &[f64] {
        &self.gyro
    }

    /// Assemble the velocity operator at one signed, steady spin rate.
    /// Bearings, preload and stiffness stay fixed exactly as supplied.
    pub fn at_speed(
        &self,
        spin_rad_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorOperatingPoint<'_>, RotorError> {
        poll(&mut cancelled)?;
        if !spin_rad_s.is_finite() {
            return Err(RotorError::Invalid("finite signed spin required"));
        }
        let mut velocity_operator = Vec::with_capacity(self.n * self.n);
        for (i, (&c, &g)) in self.damping.iter().zip(&self.gyro).enumerate() {
            if i % self.n == 0 {
                poll(&mut cancelled)?;
            }
            velocity_operator.push(finite(spin_rad_s.mul_add(g, c))?);
        }
        poll(&mut cancelled)?;
        Ok(RotorOperatingPoint {
            model: self,
            spin_rad_s,
            velocity_operator,
        })
    }

    /// Sample a fixed-coefficient Campbell spectrum in caller-supplied speed
    /// order. Returns complete pole sets, never persistent sorted-index modes
    /// or certified critical-speed brackets. Real bearing speed dependence
    /// requires newly supplied physical C/K matrices at each operating point.
    pub fn campbell_samples(
        &self,
        spin_rad_s: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<RotorPoleReport>, RotorError> {
        poll(&mut cancelled)?;
        if spin_rad_s.is_empty()
            || spin_rad_s.len() > MAX_CAMPBELL_SAMPLES
            || spin_rad_s.iter().any(|w| !w.is_finite())
        {
            return Err(RotorError::Invalid(
                "finite nonempty bounded speed samples required",
            ));
        }
        let mut samples = Vec::with_capacity(spin_rad_s.len());
        for &speed in spin_rad_s {
            samples.push(
                self.at_speed(speed, &mut cancelled)?
                    .poles(&mut cancelled)?,
            );
        }
        poll(&mut cancelled)?;
        Ok(samples)
    }
}

impl RotorOperatingPoint<'_> {
    /// Signed frozen operating speed [rad/s].
    pub const fn spin_rad_s(&self) -> f64 {
        self.spin_rad_s
    }
    /// Complete assembled C + Omega G, without modal decoupling.
    pub fn velocity_operator(&self) -> &[f64] {
        &self.velocity_operator
    }

    /// Original scaled quadratic poles; no eigenvectors or accuracy certificates.
    /// Cancellation reaches the native QR owner, including within an iteration.
    pub fn poles(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorPoleReport, RotorError> {
        let m = self.model;
        poll(&mut cancelled)?;
        let poles_per_s = quadratic_eigenvalues_with_cancel(
            &m.mass,
            &self.velocity_operator,
            &m.stiffness,
            m.n,
            ROTOR_QR_SWEEPS_PER_WINDOW,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(RotorPoleReport {
            spin_rad_s: self.spin_rad_s,
            poles_per_s,
        })
    }

    /// Solve the original full complex dynamic stiffness for prescribed force
    /// or torque phasors. Signed omega sets the excitation's phase convention.
    /// Row equilibration improves mixed-coordinate scaling; residual and power
    /// are recomputed from the ORIGINAL separate M/C/G/K actions. The requested
    /// residual is numerical equilibrium, not an amplitude-error certificate.
    /// This is a sinusoidal particular solution; unstable homogeneous transients
    /// need not decay, so active bearings do not imply an attained steady response.
    /// Singular resonance, nonfinite arithmetic and cancellation refuse atomically.
    /// No force floor, pole shift, damping clamp or heat allocation is introduced.
    pub fn harmonic_response(
        &self,
        frequency_rad_s: f64,
        force: &[C64],
        relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorHarmonicResponse, RotorError> {
        poll(&mut cancelled)?;
        let m = self.model;
        let w = frequency_rad_s;
        if !w.is_finite()
            || force.len() != m.n
            || !finite_complex(force)
            || !relative_tolerance.is_finite()
            || relative_tolerance <= 0.0
            || relative_tolerance >= 1.0
        {
            return Err(RotorError::Invalid(
                "finite frequency/force and relative tolerance in (0,1) required",
            ));
        }
        let mut matrix = vec![C64::ZERO; m.n * m.n];
        let mut displacement = force.to_vec();
        for i in 0..m.n {
            poll(&mut cancelled)?;
            let mut row_scale = 0.0_f64;
            for j in 0..m.n {
                let ij = i * m.n + j;
                let entry = C64::new(
                    finite((-w).mul_add(w * m.mass[ij], m.stiffness[ij]))?,
                    finite(w * self.velocity_operator[ij])?,
                );
                row_scale = row_scale.max(finite(entry.abs())?);
                matrix[ij] = entry;
            }
            if row_scale == 0.0 {
                return Err(LuComplexError::Singular { index: i }.into());
            }
            for entry in &mut matrix[i * m.n..(i + 1) * m.n] {
                *entry = C64::new(entry.re / row_scale, entry.im / row_scale);
            }
            displacement[i] = C64::new(force[i].re / row_scale, force[i].im / row_scale);
        }
        if !finite_complex(&displacement) {
            return Err(RotorError::Invalid("scaled forcing is unrepresentable"));
        }
        let factor = lu_complex_with_cancel(&matrix, m.n, &mut cancelled)?;
        factor.solve_with_cancel(&mut displacement, &mut cancelled)?;
        let velocity: Vec<_> = displacement
            .iter()
            .map(|q| C64::new(-w * q.im, w * q.re))
            .collect();
        if !finite_complex(&displacement) || !finite_complex(&velocity) {
            return Err(RotorError::Invalid("harmonic response is unrepresentable"));
        }
        let mut inertial = action(&m.mass, &displacement, &mut cancelled)?;
        for value in &mut inertial {
            *value = value.scale(-w).scale(w);
        }
        let damping = action(&m.damping, &velocity, &mut cancelled)?;
        let mut gyro = action(&m.gyro, &velocity, &mut cancelled)?;
        for value in &mut gyro {
            *value = value.scale(self.spin_rad_s);
        }
        let stiffness = action(&m.stiffness, &displacement, &mut cancelled)?;
        let residual: Vec<_> = (0..m.n)
            .map(|i| inertial[i] + damping[i] + gyro[i] + stiffness[i] - force[i])
            .collect();
        let norms = [
            norm(&inertial)?,
            norm(&damping)?,
            norm(&gyro)?,
            norm(&stiffness)?,
            norm(force)?,
        ];
        let scale = norms.into_iter().fold(0.0_f64, f64::max);
        let residual_norm = norm(&residual)?;
        let relative_residual = if scale == 0.0 {
            residual_norm
        } else {
            finite((residual_norm / scale) / norms.into_iter().map(|x| x / scale).sum::<f64>())?
        };
        if relative_residual > relative_tolerance {
            return Err(RotorError::Residual {
                relative: relative_residual,
                limit: relative_tolerance,
            });
        }
        let mean_input_power_w = power(force, &velocity)?;
        let mean_velocity_force_power_w = power(&damping, &velocity)?;
        let mean_stiffness_force_power_w = power(&stiffness, &velocity)?;
        let mean_gyroscopic_power_w = power(&gyro, &velocity)?;
        let mean_inertial_power_w = power(&inertial, &velocity)?;
        let mean_power_residual_w = finite(
            mean_input_power_w
                - mean_velocity_force_power_w
                - mean_stiffness_force_power_w
                - mean_gyroscopic_power_w
                - mean_inertial_power_w,
        )?;
        poll(&mut cancelled)?;
        Ok(RotorHarmonicResponse {
            spin_rad_s: self.spin_rad_s,
            frequency_rad_s: w,
            displacement,
            velocity,
            relative_residual,
            mean_input_power_w,
            mean_velocity_force_power_w,
            mean_stiffness_force_power_w,
            mean_gyroscopic_power_w,
            mean_inertial_power_w,
            mean_power_residual_w,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), RotorError> {
    if cancelled() {
        Err(RotorError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, RotorError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(RotorError::Invalid("unrepresentable rotor arithmetic"))
    }
}
fn finite_complex(values: &[C64]) -> bool {
    values.iter().all(|z| z.re.is_finite() && z.im.is_finite())
}
fn norm(values: &[C64]) -> Result<f64, RotorError> {
    values.iter().try_fold(0.0, |sum, z| {
        finite(fs_math::det::hypot(sum, finite(z.abs())?))
    })
}
fn action(
    matrix: &[f64],
    values: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<C64>, RotorError> {
    let n = values.len();
    let mut output = vec![C64::ZERO; n];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            output[i] = output[i] + values[j].scale(matrix[i * n + j]);
        }
    }
    if !finite_complex(&output) {
        return Err(RotorError::Invalid(
            "original force action is unrepresentable",
        ));
    }
    Ok(output)
}
fn power(force: &[C64], velocity: &[C64]) -> Result<f64, RotorError> {
    force
        .iter()
        .zip(velocity)
        .try_fold(0.0, |sum, (f, v)| finite(sum + 0.5 * (f.conj() * *v).re))
}
