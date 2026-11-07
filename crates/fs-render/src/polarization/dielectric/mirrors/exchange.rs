//! Heat and stationary radiation force from the original smooth reflection.

use super::{BeamSample, ConductorBeamError, ConductorOptics, PolarizedBeam, Vec3, phase_index, unit};

const LIGHT_M_S: f64 = 299_792_458.0;

/// Reflected ray power and material exchange at one smooth opaque boundary.
///
/// The incident Stokes intensity is a power weight in watts. Absorption covers
/// the supplied coherent coating and opaque substrate together; it is not
/// assigned to a thermal cell or penetration depth. For sampled input paths,
/// power and force are estimators with the same path weight as the input beam.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductorBeamBalance {
    reflection: BeamSample,
    absorbed_power_w: f64,
    vacuum_force_n: Option<Vec3>,
}

impl ConductorBeamBalance {
    /// The original phase-preserving reflection, with no second attenuation.
    pub const fn reflection(&self) -> &BeamSample {
        &self.reflection
    }

    /// Positive optical power entering the coating and opaque substrate [W].
    pub const fn absorbed_power_w(&self) -> f64 {
        self.absorbed_power_w
    }

    /// World-space force on stationary matter [N] in an index-one optical host.
    ///
    /// This is incident minus reflected photon momentum per unit time. No
    /// additional incidence cosine or surface-area factor applies to ray power.
    /// Returns `None` for a refractive host: its electromagnetic/material stress
    /// partition requires a separately specified momentum convention. Body
    /// motion, Doppler work, emitted recoil and polarization spin torque are not
    /// included. A mechanical owner may integrate the returned stationary force.
    pub const fn vacuum_force_n(&self) -> Option<Vec3> {
        self.vacuum_force_n
    }
}

impl PolarizedBeam {
    /// Reflect with the original optical law and retain absorbed heat and force.
    ///
    /// Absorption comes from the original positive boundary Poynting flux for
    /// each s/p component, including coherent coating interference. It is never
    /// obtained by subtracting rounded reflected power from incident power, so
    /// weak metal/film absorption survives reflectance that rounds to one.
    /// Arbitrary input reference axes are rotated into the incidence frame.
    ///
    /// The original `reflect_conductor` performs no additional field solve or
    /// allocation. This opt-in query adds the native boundary-flux solve and
    /// returns the same reflected beam. It does not advance a thermal or rigid
    /// body state. Force is available only for an index-one incident medium.
    pub fn reflect_conductor_with_balance(
        self,
        normal: Vec3,
        optics: &ConductorOptics,
    ) -> Result<ConductorBeamBalance, ConductorBeamError> {
        let reflection = self.reflect_conductor(normal, optics)?;
        let k = self.frame().direction();
        let mut normal = unit(normal)?;
        if normal.dot(k) > 0.0 {
            normal = normal.scale(-1.0);
        }
        let cosine = (-normal.dot(k)).clamp(0.0, 1.0);
        let eta = phase_index(self)?;
        let balance = optics.power_balance(self.wavelength_nm(), eta, cosine)?;
        let incident = self.with_reference_axis(reflection.beam.frame().horizontal())?.power();
        let absorbed_power_w = if incident.i == 0.0 {
            0.0
        } else {
            // Normalize Q before mixing: I+Q can overflow for a finite beam.
            // Apply the small absorption fraction to I before the polarization
            // weights, retaining weak absolute heat even for very large powers.
            let s = (0.5 + 0.5 * (incident.q / incident.i)).clamp(0.0, 1.0);
            (incident.i * balance.absorbed[0]) * s
                + (incident.i * balance.absorbed[1]) * (1.0 - s)
        };
        let vacuum_force_n = (eta == 1.0).then(|| {
            // P_in k - P_ref k_ref = P_abs k - 2 P_ref cos(theta) n.
            // This positive-loss form preserves tangential momentum transfer
            // even if P_ref rounds to P_in. Divide before doubling large power.
            let absorption = absorbed_power_w / LIGHT_M_S;
            let reflection = 2.0 * (reflection.beam.power().i / LIGHT_M_S) * cosine;
            Vec3::new(
                absorption * k.x - reflection * normal.x,
                absorption * k.y - reflection * normal.y,
                absorption * k.z - reflection * normal.z,
            )
        });
        Ok(ConductorBeamBalance { reflection, absorbed_power_w, vacuum_force_n })
    }
}

#[cfg(test)]
mod tests;
