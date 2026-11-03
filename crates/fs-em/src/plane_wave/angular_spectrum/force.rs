//! Integrated periodic-cell Maxwell force without reconstructing a field plane.

use super::*;
use fs_math::dd::Dd;

impl PlanarAngularSpectrum<'_, '_> {
    /// Cycle-mean stack/interface electromagnetic stress resultant [N] over
    /// the original periodic cell, in the beam's x/y/z coordinates. Each
    /// exterior retains its own original positive lossless material law.
    /// Material striction and ambient mechanical pressure remain excluded.
    ///
    /// Parseval orthogonality eliminates distinct transverse-mode cross terms
    /// in this FULL-cell integral. Coherent TE/TM film phases have already
    /// determined each mode's reflection and positive physical absorption.
    /// Per-mode direction is retained; replacing it by one incidence angle or
    /// total absorbed watts would lose lateral momentum. O(N), allocation-free,
    /// no additional FFT or constitutive solve. Cancellation publishes no force.
    ///
    /// This is not a local force density, finite-edge or torque calculation,
    /// nor a moving-body power closure. The load must not be counted twice by
    /// a separate Maxwell/Lorentz mechanical consumer.
    pub fn radiation_force(
        &self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[f64; 3], AngularSpectrumError> {
        poll(&mut cancelled)?;
        self.source.admit_refracting_momentum_exterior()?;
        let omega = self.source.angular_frequency_rad_s();
        let mut sum = [Dd::ZERO; 3];
        for (index, mode) in self.modes.iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let force = super::super::momentum::one_sided(
                [mode.kx / omega, mode.ky / omega, mode.incident_kz / omega],
                self.source
                    .momentum_normal_difference(mode.incidence_cosine, mode.exit_kz),
                mode.reflected_power_w,
                mode.transmitted_power_w,
                mode.absorbed_power_w,
            );
            for (s, value) in sum.iter_mut().zip(force) {
                if !value.is_finite() {
                    return Err(bad("beam radiation force is unrepresentable"));
                }
                *s = *s + Dd::from_f64(value);
            }
        }
        let force = sum.map(Dd::to_f64);
        if force.iter().any(|f| !f.is_finite()) {
            return Err(bad("integrated radiation force is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(force)
    }
}
