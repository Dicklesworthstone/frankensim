//! Read-only released-contact response and the original physical observation.
//! This is an AC particular solution, not a bowed operating-point linearization.

use super::*;
use crate::render::schedule::force::coupled::harmonic::ModalHarmonicResponse;
use fs_math::c64::C64;

/// Peak phasors in exp(-i omega t), with the same receiver as the time owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ReleasedHarmonicResponse {
    /// Complete original modal, Maxwell and acoustic response and power account.
    pub dynamics: ModalHarmonicResponse,
    /// Original compact acceleration or real modal observation [Pa].
    /// This does not add a finite-distance BEM receiver, delay or directivity.
    pub observer_pressure_pa: C64,
}

impl FrictionModalSystem {
    /// Forced linear response with the dry bow/contact released. Every original
    /// string, body, cavity, bilateral connection, Maxwell branch and selected
    /// radiation pole remains in the original coupled harmonic solve.
    ///
    /// No accepted state, contact branch, initial history or clock is modified.
    /// The AC particular solution excludes DC preload and initial transients.
    /// Nonlinear/adhesive storage refuses: it is never silently linearized or
    /// removed. Source frequency support remains with the coefficient producer.
    ///
    /// The physical pressure readout is the exact real observation of the time
    /// owner. For a narrow-band modal transfer T that readout is Re(T)*v plus
    /// Im(T)*omega_mode*q; multiplying T by v at an arbitrary drive frequency
    /// would be a different observer. Acceleration weights instead use -omega²q.
    pub fn released_harmonic_response(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ReleasedHarmonicResponse, ModalFrictionError> {
        poll(Some(gate))?;
        if self.nonlinear.is_some() || self.adhesive.is_some() {
            return Err(invalid(
                "released harmonic response cannot omit nonlinear or adhesive storage",
            )
            .into());
        }
        let dynamics = self
            .network
            .harmonic_response(omega_rad_s, force, maximum_work, gate)?;
        let mut pressure = C64::ZERO;
        if let Some(weights) = &self.acceleration_pressure {
            if weights.len() != dynamics.displacement_m_sqrt_kg.len() {
                return Err(
                    invalid("harmonic pressure lost its complete original modal map").into(),
                );
            }
            let minus_omega_squared = -omega_rad_s * omega_rad_s;
            for (chunk, weight_chunk) in dynamics
                .displacement_m_sqrt_kg
                .chunks(256)
                .zip(weights.chunks(256))
            {
                poll(Some(gate))?;
                for (&q, &weight) in chunk.iter().zip(weight_chunk) {
                    pressure = pressure + q.scale(minus_omega_squared * weight);
                }
            }
        } else {
            let mut start = 0;
            for component in self.network.components() {
                poll(Some(gate))?;
                let mut component_pressure = C64::ZERO;
                for (index, mode) in component.modes().iter().enumerate() {
                    let q = dynamics.displacement_m_sqrt_kg[start + index];
                    let v = dynamics.velocity_m_sqrt_kg_per_s[start + index];
                    let transfer = mode.pressure_per_modal_velocity;
                    component_pressure = component_pressure
                        + v.scale(transfer.re)
                        + q.scale(transfer.im * mode.angular_frequency_rad_s);
                }
                limit(
                    "harmonic component observer pressure",
                    component_pressure.abs(),
                    component.budget().maximum_abs_pressure_pa,
                )?;
                pressure = pressure + component_pressure;
                start += component.modes().len();
            }
        }
        if !pressure.re.is_finite() || !pressure.im.is_finite() {
            return Err(invalid("harmonic observer pressure is not representable").into());
        }
        limit(
            "harmonic observer pressure",
            pressure.abs(),
            self.network.config.maximum_abs_pressure_pa,
        )?;
        Ok(ReleasedHarmonicResponse {
            dynamics,
            observer_pressure_pa: pressure,
        })
    }
}
