//! Continuous sinusoidal response of the original shared gas/body/neck network.
//! All values are AC peak phasors in exp(-i omega t), about the fixed reference
//! equilibrium. This does not advance accepted histories or approximate a
//! sampled forcing sequence. The original constitutive and radiation owners
//! still determine the validity of their physical approximation.

use super::*;
use crate::render::schedule::force::coupled::harmonic::ModalHarmonicResponse;

/// Physical cavity observations from the SAME loaded harmonic solve.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedCavityHarmonicResponse {
    /// Original modal, thermal/viscous-memory and radiation response, with its
    /// complete cycle-mean power and equation residuals.
    pub network: ModalHarmonicResponse,
    /// Increase in effective cavity volume [m3], including outward neck travel.
    pub volume_change_m3: C64,
    /// Signed small-signal volume change divided by the reference volume.
    pub fractional_volume_change: C64,
    /// Uniform internal pressure increment [Pa]. Includes thermal wall memory.
    pub pressure_increment_pa: C64,
    /// Sum of original moving-surface outward volume rates [m3/s].
    pub surface_volume_rate_m3_s: C64,
    /// Physical outward neck volume flow [m3/s].
    pub neck_outflow_m3_s: C64,
    /// Exterior reaction pressure opposing neck flow [Pa], including signed
    /// surface/neck mutual radiation. Zero without an attached exterior model.
    pub neck_exterior_pressure_pa: C64,
    /// Cycle-mean gas compression availability already in network storage [J].
    pub mean_compression_energy_j: f64,
    /// Cycle-mean physical neck kinetic energy already in network storage [J].
    pub mean_neck_kinetic_energy_j: f64,
    /// Mean power supplied by original surface forces [W]. May be signed.
    pub surface_input_power_w: f64,
    /// Mean power supplied by incident exterior neck pressure [W]. May be signed.
    pub incident_pressure_input_power_w: f64,
}

impl SharedCavity {
    /// Solve the loaded continuous AC response in exp(-i omega t).
    ///
    /// Surface forces retain their ORIGINAL flattened mass-normalized basis.
    /// Positive incident pressure acts inward on the neck; the computed cavity
    /// and exterior reaction pressures are already internal to this solve.
    /// Initial displacement, accepted history and sample count do not enter a
    /// steady harmonic response, and none of them changes. DC and a singular
    /// undamped resonance refuse rather than report a finite invented peak.
    ///
    /// `maximum_work` bounds the existing network harmonic assembly/solve.
    /// Queries must stay inside each admitted viscous/thermal diffusion band.
    /// Raw exterior models do not carry a validity band: their producer's band
    /// and spatial accuracy remain binding, including for geometry-derived fits.
    pub fn harmonic_response(
        &self,
        omega_rad_s: f64,
        surface_force: &[C64],
        incident_neck_pressure_pa: C64,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<SharedCavityHarmonicResponse, ModalCouplingError> {
        let force = self.prepare_harmonic_force(
            omega_rad_s,
            surface_force,
            incident_neck_pressure_pa,
            gate,
        )?;
        let network = self
            .system
            .harmonic_response(omega_rad_s, &force, maximum_work, gate)?;
        self.observe_harmonic_network(
            omega_rad_s,
            surface_force,
            incident_neck_pressure_pa,
            network,
            gate,
        )
    }

    pub(super) fn prepare_harmonic_force(
        &self,
        omega_rad_s: f64,
        surface_force: &[C64],
        incident_neck_pressure_pa: C64,
        gate: &CancelGate,
    ) -> Result<Vec<C64>, ModalCouplingError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled);
        }
        let finite = |z: C64| z.re.is_finite() && z.im.is_finite();
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || surface_force.len() != self.surface_modes
            || surface_force.iter().any(|&z| !finite(z))
            || !finite(incident_neck_pressure_pa)
        {
            return Err(network_invalid(
                "harmonic cavity needs positive frequency and finite original physical forces",
            ));
        }
        for band in self
            .viscous_neck
            .as_ref()
            .map(|m| m.band_rad_s())
            .into_iter()
            .chain(
                self.thermal_compression
                    .as_ref()
                    .map(|m| m.band_rad_s()),
            )
        {
            if !(band[0]..=band[1]).contains(&omega_rad_s) {
                return Err(network_invalid(
                    "harmonic cavity query is outside an admitted wall-diffusion band",
                ));
            }
        }
        let projection = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let mut force = vec![C64::ZERO; self.system.mode_count()];
        force[..self.surface_modes].copy_from_slice(surface_force);
        force[self.surface_modes] = incident_neck_pressure_pa.scale(-projection);
        Ok(force)
    }

    pub(super) fn observe_harmonic_network(
        &self,
        omega_rad_s: f64,
        surface_force: &[C64],
        incident_neck_pressure_pa: C64,
        network: ModalHarmonicResponse,
        gate: &CancelGate,
    ) -> Result<SharedCavityHarmonicResponse, ModalCouplingError> {
        let finite = |z: C64| z.re.is_finite() && z.im.is_finite();
        let projection = self.source.neck_volume_shape_m2_per_sqrt_kg();
        let extension = self
            .column
            .iter()
            .zip(&network.displacement_m_sqrt_kg)
            .fold(C64::ZERO, |sum, (&c, &q)| sum + q.scale(c));
        let area = self.source.neck_area_m2();
        let volume = extension.scale(area);
        let surface_rate = self.column[..self.surface_modes]
            .iter()
            .zip(&network.velocity_m_sqrt_kg_per_s)
            .fold(C64::ZERO, |sum, (&c, &v)| sum + v.scale(c * area));
        let neck_velocity = network.velocity_m_sqrt_kg_per_s[self.surface_modes];
        let neck_flow = neck_velocity.scale(projection);
        let (compression_force, compression_energy) = match (&self.thermal_compression, &self.molecular_compression) {
            (Some(law), _) => {
                // Reuse the original constitutive law; conjugation converts its
                // documented +i omega convention to this solve's -i omega.
                let compression_force = law.stiffness_at(omega_rad_s)?.conj() * extension;
                let s = C64::new(0.0, -omega_rad_s);
                let mut energy = 0.25 * law.equilibrium_stiffness_n_m() * extension.norm_sq();
                for arm in law.arms() {
                    // Evaluate the relaxed-strain difference directly, avoiding
                    // subtracting nearly equal histories at low frequency.
                    let slip = s / (s + C64::from_re(arm.relaxation_time_s.recip())) * extension;
                    energy += 0.25 * arm.stiffness_n_m * slip.norm_sq();
                }
                (compression_force, energy)
            }
            (_, Some(law)) => law.harmonic_observation(omega_rad_s, extension)?,
            (None, None) => {
                let stiffness = self.source.compression_stiffness_n_m();
                (
                    extension.scale(stiffness),
                    0.25 * stiffness * extension.norm_sq(),
                )
            }
        };
        let mut exterior = C64::ZERO;
        if let (Some(port), Some(model)) = (self.radiation_neck_port, self.system.radiation_model())
        {
            for (pole, state) in model
                .poles
                .iter()
                .zip(network.radiation_state.chunks_exact(2))
            {
                exterior = exterior + state[1].scale(pole.coupling[port] / projection);
            }
        }
        let surface_power = surface_force
            .iter()
            .zip(&network.velocity_m_sqrt_kg_per_s)
            .map(|(&f, &v)| 0.5 * (f * v.conj()).re)
            .sum();
        let result = SharedCavityHarmonicResponse {
            volume_change_m3: volume,
            fractional_volume_change: volume.scale(self.source.geometry().volume_m3.recip()),
            pressure_increment_pa: compression_force.scale(-area.recip()),
            surface_volume_rate_m3_s: surface_rate,
            neck_outflow_m3_s: neck_flow,
            neck_exterior_pressure_pa: exterior,
            mean_compression_energy_j: compression_energy,
            mean_neck_kinetic_energy_j: 0.25 * neck_velocity.norm_sq(),
            surface_input_power_w: surface_power,
            incident_pressure_input_power_w: -0.5
                * (incident_neck_pressure_pa * neck_flow.conj()).re,
            network,
        };
        if [
            result.volume_change_m3,
            result.fractional_volume_change,
            result.pressure_increment_pa,
            result.surface_volume_rate_m3_s,
            result.neck_outflow_m3_s,
            result.neck_exterior_pressure_pa,
        ]
        .iter()
        .any(|&z| !finite(z))
            || [
                result.mean_compression_energy_j,
                result.mean_neck_kinetic_energy_j,
                result.surface_input_power_w,
                result.incident_pressure_input_power_w,
            ]
            .iter()
            .any(|x| !x.is_finite())
        {
            return Err(network_invalid(
                "harmonic cavity physical observation overflows",
            ));
        }
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled);
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
