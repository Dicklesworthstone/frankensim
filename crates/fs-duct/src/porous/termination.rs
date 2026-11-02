//! Passive rear loading of the original resolved porous volume fields.
//!
//! The load's p/U impedance uses exp(-i omega t), with U leaving the porous
//! specimen. It describes one frequency and conserved tangential mode at the
//! physical rear face. Any backing/gap/radiation model and its source support
//! remain the caller's responsibility; no material or reflection law is fitted.

use super::{DuctError, PorousLayerStack, PorousStackExcitation, bad, finite};
use crate::compose::feedback;
use fs_math::{c64::C64, dd::Dd};

/// A fixed passive terminal impedance [Pa s/m³], or an exact rigid cap.
/// The finite impedance already includes its physical face area and normal
/// modal convention. A specific p/v impedance must first be divided by area.
/// A load may absorb energy, but its heat never belongs to the porous cells.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PassivePorousTermination {
    impedance: Option<C64>,
}

impl PassivePorousTermination {
    /// Exact zero normal terminal flow; no infinite impedance is formed.
    pub const fn rigid() -> Self {
        Self { impedance: None }
    }

    /// Admit a finite passive pressure/volume-flow impedance. Its real part
    /// must be nonnegative; zero impedance is an exact pressure-release face.
    pub fn try_new_impedance(impedance_pa_s_m3: C64) -> Result<Self, DuctError> {
        if !finite(impedance_pa_s_m3) || impedance_pa_s_m3.re < 0.0 {
            return Err(bad(
                "porous termination needs a finite passive volume impedance",
            ));
        }
        Ok(Self {
            impedance: Some(impedance_pa_s_m3),
        })
    }

    /// Original finite load, or `None` for a rigid cap.
    pub const fn impedance_pa_s_m3(self) -> Option<C64> {
        self.impedance
    }

    // Normalize before adding Z+R: a finite passive load cannot have a pole,
    // and very large impedances must not overflow that sum. Keep Re(Z) itself
    // for positive dissipated work; 1-|Gamma|² erases weak load absorption.
    fn reflection(self, reference: f64) -> C64 {
        self.impedance.map_or(C64::ONE, |z| {
            let scale = z.re.abs().max(z.im.abs()).max(reference);
            let z = C64::new(z.re / scale, z.im / scale);
            let r = C64::from_re(reference / scale);
            (z - r) / (z + r)
        })
    }

    fn absorbed_power(self, outgoing: C64, reference: f64) -> Result<f64, DuctError> {
        let Some(z) = self.impedance else {
            return Ok(0.0);
        };
        if z.re == 0.0 || outgoing == C64::ZERO {
            return Ok(0.0);
        }
        let scale = z.re.abs().max(z.im.abs()).max(reference);
        let denominator = C64::new(z.re / scale + reference / scale, z.im / scale);
        // P_load = Re(Z)|U|²/2, U = 2*b_R/(Z+R). Square only after
        // power normalization. Square roots retain weak positive Re(Z) even
        // when its ratio to the reactive impedance would underflow.
        let normalized = (outgoing / denominator).scale(1.0 / scale.sqrt());
        let weighted = normalized.scale((z.re.sqrt() / scale.sqrt()) * core::f64::consts::SQRT_2);
        let power = weighted.norm_sq();
        if finite(weighted) && power.is_finite() && power > 0.0 {
            Ok(power)
        } else {
            Err(bad(
                "positive porous terminal dissipated power is unrepresentable",
            ))
        }
    }
}

/// Acoustic power supplied externally to a resolved porous specimen [W].
/// An attached load is separate from the porous thermal control volumes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PorousSourcePower {
    /// Authored source waves only; returned rear pressure is never a new source.
    pub incident_power_w: f64,
    /// Outgoing waves at prescribed-source faces (left only when terminated).
    pub outgoing_power_w: f64,
    /// Positive work exported into a passive rear load; never pore-cell heat.
    pub load_absorbed_power_w: f64,
    /// Original positive viscous plus thermal work in the porous layers.
    pub porous_absorbed_power_w: f64,
    /// Incident minus outgoing minus load work minus original pore work.
    pub residual_w: f64,
}

impl<'l, 'm> PorousLayerStack<'l, 'm> {
    /// Illuminate the left face while a passive load closes the right face.
    /// The solved counter-wave is a_R=Gamma*S21*a_L/(1-Gamma*S22).
    /// Native feedback conditioning and original decaying-wave reconstruction
    /// apply; no inverse transmission, growing wave or new material query is
    /// used. `boundary_power` remains the RAW two-face stack diagnostic;
    /// `source_power` excludes recirculation and separates terminal absorption.
    pub fn prepare_terminated_excitation(
        &self,
        incident_left_pressure_pa: C64,
        termination: PassivePorousTermination,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PorousStackExcitation<'_, 'l, 'm>, DuctError> {
        if cancelled() {
            return Err(DuctError::Cancelled { completed_steps: 0 });
        }
        if !finite(incident_left_pressure_pa) {
            return Err(bad("porous prescribed incident pressure must be finite"));
        }
        let port = self.two_port();
        let gamma = termination.reflection(port.reference_impedance);
        let (denominator, condition) = feedback(port.scattering[3], gamma);
        if !finite(gamma) || !condition.is_finite() || condition > self.maximum_condition {
            return Err(DuctError::Singular);
        }
        let returned = (gamma * port.scattering[2] / denominator) * incident_left_pressure_pa;
        if !finite(returned) {
            return Err(DuctError::Singular);
        }
        let mut excitation =
            self.prepare_excitation([incident_left_pressure_pa, returned], &mut cancelled)?;
        let boundary = excitation.boundary_power();
        let normalization = 1.0 / port.reference_impedance.sqrt();
        let incoming = 0.5 * incident_left_pressure_pa.scale(normalization).norm_sq();
        let outgoing = 0.5
            * boundary.outgoing_pressure_pa[0]
                .scale(normalization)
                .norm_sq();
        let load = termination
            .absorbed_power(boundary.outgoing_pressure_pa[1], port.reference_impedance)?;
        let porous = excitation.dissipation().absorbed_power_w;
        let residual = (Dd::from_f64(incoming)
            - Dd::from_f64(outgoing)
            - Dd::from_f64(load)
            - Dd::from_f64(porous))
        .to_f64();
        let scale = incoming.max(outgoing).max(load).max(porous);
        if ![incoming, outgoing, load, porous, residual]
            .into_iter()
            .all(f64::is_finite)
            || residual.abs() > 512.0 * f64::EPSILON * (self.layers().len() + 1) as f64 * scale
        {
            return Err(bad(
                "porous source, terminal load and original cell heat do not close",
            ));
        }
        excitation.source_power = PorousSourcePower {
            incident_power_w: incoming,
            outgoing_power_w: outgoing,
            load_absorbed_power_w: load,
            porous_absorbed_power_w: porous,
            residual_w: residual,
        };
        excitation.condition_estimate = self.condition_estimate().max(condition);
        Ok(excitation)
    }
}

#[cfg(test)]
mod tests;
