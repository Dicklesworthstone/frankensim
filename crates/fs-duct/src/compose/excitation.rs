//! Coherent boundary fields and physical power for the original two-port.

use super::{DuctError, TwoPort, bad, finite};
use fs_math::c64::C64;

/// Coherent, cycle-averaged boundary work under two-sided excitation.
/// Pressure phasors are peak amplitudes, so power includes the factor 1/2.
/// Right-port volume flow points INTO the block, opposite to axial flow.
#[derive(Debug, Clone, Copy)]
pub struct TwoPortPower {
    /// Outgoing pressure-wave amplitudes at the left and right faces [Pa].
    pub outgoing_pressure_pa: [C64; 2],
    /// Total pressure at each face, including both incident fields [Pa].
    pub pressure_pa: [C64; 2],
    /// Volume flow into each face [m³/s].
    pub inward_volume_flow_m3_s: [C64; 2],
    /// Total power in the two prescribed incoming waves [W].
    pub incident_power_w: f64,
    /// Power in the two coherently combined outgoing waves [W].
    pub outgoing_power_w: f64,
    /// Signed incoming-minus-outgoing power [W], retaining roundoff.
    /// The original block owner determines its physical destination. This is
    /// not an independently integrated heat source or a passivity certificate.
    pub absorbed_power_w: f64,
}

impl TwoPort {
    /// Observe the original scattering block without allocation or a new solve.
    /// Incoming/outgoing pressure amplitudes use the same positive real volume
    /// impedance [Pa s/m³] as this block. Interfering amplitudes combine BEFORE
    /// power is evaluated; transmitted power is outgoing, not absorbed.
    pub fn power_for_incident_pressure(
        self,
        incident_pressure_pa: [C64; 2],
    ) -> Result<TwoPortPower, DuctError> {
        Self::new(self.scattering, self.reference_impedance)?;
        if !incident_pressure_pa.into_iter().all(finite) {
            return Err(bad("two-port incident pressure must be finite"));
        }
        // Preserve the original porous-layer and stack arithmetic.
        let s = self.scattering;
        let incident = incident_pressure_pa;
        let outgoing = [
            s[0] * incident[0] + s[1] * incident[1],
            s[2] * incident[0] + s[3] * incident[1],
        ];
        let reference = self.reference_impedance;
        let pressure = std::array::from_fn(|i| incident[i] + outgoing[i]);
        let flow = std::array::from_fn(|i| (incident[i] - outgoing[i]).scale(1.0 / reference));
        let amplitude_scale = 1.0 / reference.sqrt();
        let incident_power = 0.5
            * incident
                .into_iter()
                .map(|p| p.scale(amplitude_scale).norm_sq())
                .sum::<f64>();
        let outgoing_power = 0.5
            * outgoing
                .into_iter()
                .map(|p| p.scale(amplitude_scale).norm_sq())
                .sum::<f64>();
        let absorbed = incident_power - outgoing_power;
        if !outgoing.into_iter().all(finite)
            || !pressure.into_iter().all(finite)
            || !flow.into_iter().all(finite)
            || !incident_power.is_finite()
            || !outgoing_power.is_finite()
            || !absorbed.is_finite()
        {
            return Err(bad("two-port field or power is unrepresentable"));
        }
        Ok(TwoPortPower {
            outgoing_pressure_pa: outgoing,
            pressure_pa: pressure,
            inward_volume_flow_m3_s: flow,
            incident_power_w: incident_power,
            outgoing_power_w: outgoing_power,
            absorbed_power_w: absorbed,
        })
    }
}
