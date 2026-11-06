//! Special-relativistic observation of one vacuum photon path.
//!
//! All velocities and the photon propagation direction use the SAME inertial
//! world frame. The direction points from emission toward reception, opposite
//! a camera's tracing direction. The emitter and receiver need not be at rest
//! in that frame. Their local frames use rotation-free Lorentz boosts.
//!
//! For `D = nu_receiver / nu_emitter`, invariance of `I_nu / nu^3` gives
//! `I_lambda,receiver(lambda) = D^5 I_lambda,emitter(D lambda)`; the fifth power
//! includes the frequency-to-wavelength Jacobian. The associated solid angle
//! changes by `D^-2`, already accounted for in radiance: do not multiply it
//! into radiance a second time. See Younsi, Wu & Fuerst (2012), equations 9–13,
//! <https://doi.org/10.1051/0004-6361/201219599>.
//!
//! This is a local photon-frame transform, not a spacetime intersection solver.
//! Callers supply the direction and endpoint velocities at the actual emission
//! and reception events. No instantaneous-position approximation, moving-medium
//! index, gravitational shift, scattering recoil or polarization is inferred.

use fs_math::det;

/// Exact SI speed of light in vacuum, metres per second.
pub const VACUUM_LIGHT_SPEED_M_S: f64 = 299_792_458.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DopplerError {
    InvalidDirection,
    /// Each finite velocity must have magnitude strictly below c.
    InvalidVelocity,
    InvalidWavelength,
    InvalidRadiance,
    /// A positive transform or output is not representable in binary64.
    NumericalRange,
}

impl core::fmt::Display for DopplerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "vacuum Doppler: {self:?}")
    }
}
impl core::error::Error for DopplerError {}

/// Immutable, allocation-free transform prepared once per physical ray.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VacuumDoppler {
    frequency_ratio: f64,
    source_direction: [f64; 3],
    receiver_direction: [f64; 3],
}

impl VacuumDoppler {
    /// Admit one propagation direction (normalized here) and two SI velocities.
    /// Source approach increases D; receiver motion along the photon decreases
    /// it. Identical source/receiver velocities give exactly D=1 even when their
    /// common transverse speed is relativistic.
    pub fn try_new(
        world_propagation_direction: [f64; 3],
        source_velocity_m_s: [f64; 3],
        receiver_velocity_m_s: [f64; 3],
    ) -> Result<Self, DopplerError> {
        let direction = unit(world_propagation_direction)?;
        let (source_frequency, source_direction) = in_frame(direction, source_velocity_m_s)?;
        let (receiver_frequency, receiver_direction) = in_frame(direction, receiver_velocity_m_s)?;
        let frequency_ratio = positive(receiver_frequency / source_frequency)?;
        Ok(Self {
            frequency_ratio,
            source_direction,
            receiver_direction,
        })
    }

    /// nu_receiver / nu_source. Greater than one is a blueshift.
    pub const fn frequency_ratio(&self) -> f64 {
        self.frequency_ratio
    }

    /// Photon propagation direction in the emitting material's inertial frame.
    pub const fn source_direction(&self) -> [f64; 3] {
        self.source_direction
    }

    /// Photon propagation direction in the receiver's inertial frame.
    pub const fn receiver_direction(&self) -> [f64; 3] {
        self.receiver_direction
    }

    /// Source vacuum wavelength corresponding to a receiver vacuum wavelength.
    pub fn source_wavelength_nm(&self, receiver_wavelength_nm: f64) -> Result<f64, DopplerError> {
        wavelength(receiver_wavelength_nm)?;
        positive(receiver_wavelength_nm * self.frequency_ratio)
    }

    /// Receiver vacuum wavelength corresponding to a source vacuum wavelength.
    pub fn receiver_wavelength_nm(&self, source_wavelength_nm: f64) -> Result<f64, DopplerError> {
        wavelength(source_wavelength_nm)?;
        positive(source_wavelength_nm / self.frequency_ratio)
    }

    /// Transform source radiance at the CORRESPONDING source wavelength, in
    /// W/(m² sr nm), to receiver radiance. This does not evaluate a spectrum;
    /// obtain that wavelength with `source_wavelength_nm` before source lookup.
    pub fn receiver_radiance_w_m2_sr_nm(&self, source_radiance: f64) -> Result<f64, DopplerError> {
        if !source_radiance.is_finite() || source_radiance < 0.0 {
            return Err(DopplerError::InvalidRadiance);
        }
        if source_radiance == 0.0 || self.frequency_ratio == 1.0 {
            return Ok(source_radiance);
        }
        let d2 = self.frequency_ratio * self.frequency_ratio;
        let direct = source_radiance * (d2 * d2 * self.frequency_ratio);
        if direct.is_finite() && direct > 0.0 {
            Ok(direct)
        } else {
            positive(det::exp(
                det::ln(source_radiance) + 5.0 * det::ln(self.frequency_ratio),
            ))
        }
    }

    /// dOmega_receiver / dOmega_source, for the corresponding infinitesimal
    /// ray bundle. Radiance already contains this transform.
    pub fn receiver_solid_angle_ratio(&self) -> f64 {
        let inverse = 1.0 / self.frequency_ratio;
        inverse * inverse
    }
}

fn wavelength(value: f64) -> Result<(), DopplerError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(DopplerError::InvalidWavelength)
    }
}

fn positive(value: f64) -> Result<f64, DopplerError> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(DopplerError::NumericalRange)
    }
}

fn unit(vector: [f64; 3]) -> Result<[f64; 3], DopplerError> {
    if vector.iter().any(|x| !x.is_finite()) {
        return Err(DopplerError::InvalidDirection);
    }
    let scale = vector.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    if scale == 0.0 {
        return Err(DopplerError::InvalidDirection);
    }
    let scaled = vector.map(|x| x / scale);
    let norm = det::sqrt(dot(scaled, scaled));
    Ok(scaled.map(|x| x / norm))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn in_frame(direction: [f64; 3], velocity_m_s: [f64; 3]) -> Result<(f64, [f64; 3]), DopplerError> {
    if velocity_m_s.iter().any(|x| !x.is_finite()) {
        return Err(DopplerError::InvalidVelocity);
    }
    let beta = velocity_m_s.map(|x| x / VACUUM_LIGHT_SPEED_M_S);
    let speed = det::sqrt(dot(beta, beta));
    if !speed.is_finite() || speed >= 1.0 {
        return Err(DopplerError::InvalidVelocity);
    }
    if speed == 0.0 {
        return Ok((1.0, direction));
    }
    let axis = beta.map(|x| x / speed);
    // Inputs have been normalized; this clamp removes only dot-product
    // roundoff at parallel/antiparallel endpoints, never a physical speed cap.
    let cosine = dot(direction, axis).clamp(-1.0, 1.0);
    let inverse_gamma = det::sqrt((1.0 - speed) * (1.0 + speed));
    let denominator = (1.0 - speed) + speed * (1.0 - cosine);
    let frequency = positive(denominator / inverse_gamma)?;
    let longitudinal = (cosine - speed) / denominator;
    let transverse = inverse_gamma / denominator;
    // Parallel/perpendicular decomposition avoids subtracting gamma-sized
    // vectors when a ray follows a nearly light-speed frame.
    let aberrated = core::array::from_fn(|i| {
        (direction[i] - cosine * axis[i]) * transverse + axis[i] * longitudinal
    });
    Ok((frequency, unit(aberrated)?))
}

#[cfg(test)]
mod tests;
