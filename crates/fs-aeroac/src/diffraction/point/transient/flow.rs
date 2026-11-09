//! The causal Lorentz image of the original stationary half-plane kernel.
//!
//! For an edge-parallel mean velocity U, beta²=1-(U/c)² and z'=dz/beta.
//! The original harmonic point-flow owner uses
//! G_U=exp(-i*omega*U*dz/(beta²*c²))*G_0(k/beta)/beta.
//! On the ORIGINAL real source clock this is a retarded delay
//! tau=(sqrt(beta²*r_perp²+dz²)-(U/c)*dz)/(beta²*c).
//! The same field transformation gives
//! p=p0/beta³-rho*U*u0z/beta²,
//! u_perp=u0_perp/beta, uz=u0z/beta²-U*p0/(K*beta³).
//! This is not a replacement physical sound speed or a Doppler-shifted source.
//! Reference: Griesmaier/Raumer, §2.1, Proposition 2.1 and Eq.(2.2),
//! https://www.waves.kit.edu/downloads/CRC1173_Preprint_2022-16.pdf.
//!
//! Delay evaluation rationalizes the downstream difference. The integrand
//! uses tau_edge + (L/(beta*c))*(sqrt(1+w²)-1), and resolves that last
//! difference as w²/(sqrt(1+w²)+1). No nearly cancelling global clock
//! advance is subtracted from a large stretched flight time.

use super::{AeroacError, Dd, FluidAcousticState, det, invalid};

pub(super) const SCOPE: &str = "Causal stationary volume-flow history diffracted by one infinite ideal rigid or pressure-release half-plane in original homogeneous lossless fluid with uniform strictly subsonic flow parallel to its edge. Original laboratory source clock, density and modulus; full pressure and perturbation velocity; instantaneous pressure work and total perturbation-energy transport including mean-fluid advection. C1 cubic source interpolation and complete edge memory within the declared horizon. Quadrature estimates exclude source interpolation, geometry and floating-point error. No cross-edge or sheared flow, source motion, finite barrier, independent momentum source, impedance wall, absorption, feedback or measured-SPL qualification.";

pub(super) struct FlowImage {
    pub(super) stretched_axial_separation_m: f64,
    pub(super) operator_speed_m_s: f64,
    pub(super) edge_delay_s: f64,
    pub(super) direct_delays_s: [f64; 2],
    inverse_beta: f64,
    inverse_beta_squared: f64,
    pressure_gain: f64,
    pressure_from_velocity: f64,
    velocity_from_pressure: f64,
}
impl FlowImage {
    pub(super) fn new(
        acoustic: FluidAcousticState,
        source: [f64; 3],
        point: [f64; 3],
        velocity: f64,
    ) -> Result<Self, AeroacError> {
        let c = acoustic.sound_speed_m_s();
        let mach = velocity / c;
        if !velocity.is_finite() || !mach.is_finite() || mach.abs() >= 1.0 {
            return Err(invalid(
                "causal edge flow must be finite, uniform and strictly subsonic",
            ));
        }
        let beta_squared = (1.0 - mach) * (1.0 + mach);
        let beta = det::sqrt(beta_squared);
        let inverse_beta = 1.0 / beta;
        let inverse_beta_squared = 1.0 / beta_squared;
        let pressure_gain = inverse_beta * inverse_beta_squared;
        let pressure_from_velocity = (Dd::from_f64(acoustic.density_kg_m3())
            * Dd::from_f64(velocity)
            * Dd::from_f64(inverse_beta_squared))
        .to_f64();
        let velocity_from_pressure = (Dd::from_f64(velocity) * Dd::from_f64(pressure_gain)
            / Dd::from_f64(acoustic.adiabatic_bulk_modulus_pa()))
        .to_f64();
        let dz = (Dd::from_f64(point[2]) - Dd::from_f64(source[2])).to_f64();
        let axial = dz * inverse_beta;
        let speed = c * beta;
        let delay = |transverse: f64| {
            let physical_radius = det::hypot(transverse, dz);
            if physical_radius == 0.0 {
                return 0.0; // An occluded image may coincide; visible pairs refuse elsewhere.
            }
            let mach_radius = det::hypot(beta * transverse, dz);
            let drift = Dd::from_f64(mach) * Dd::from_f64(dz);
            if mach * dz > 0.0 {
                ((Dd::from_f64(physical_radius) / Dd::from_f64(c))
                    * (Dd::from_f64(physical_radius) / (Dd::from_f64(mach_radius) + drift)))
                    .to_f64()
            } else {
                ((Dd::from_f64(mach_radius) - drift) / Dd::from_f64(beta_squared) / Dd::from_f64(c))
                    .to_f64()
            }
        };
        let edge_delay = delay(det::hypot(source[0], source[1]) + det::hypot(point[0], point[1]));
        let dx = (Dd::from_f64(point[0]) - Dd::from_f64(source[0])).to_f64();
        let direct_delays = core::array::from_fn(|j| {
            let dy = (Dd::from_f64(point[1])
                - Dd::from_f64(if j == 0 { source[1] } else { -source[1] }))
            .to_f64();
            delay(det::hypot(dx, dy))
        });
        if [
            axial,
            speed,
            inverse_beta,
            inverse_beta_squared,
            pressure_gain,
            pressure_from_velocity,
            velocity_from_pressure,
            edge_delay,
        ]
        .iter()
        .chain(&direct_delays)
        .any(|x| !x.is_finite())
            || speed <= 0.0
            || edge_delay <= 0.0
            || direct_delays.iter().any(|x| *x < 0.0)
            || (velocity != 0.0 && (pressure_from_velocity == 0.0 || velocity_from_pressure == 0.0))
        {
            return Err(invalid(
                "causal Lorentz geometry, delays or physical fields are unrepresentable",
            ));
        }
        Ok(Self {
            stretched_axial_separation_m: axial,
            operator_speed_m_s: speed,
            edge_delay_s: edge_delay,
            direct_delays_s: direct_delays,
            inverse_beta,
            inverse_beta_squared,
            pressure_gain,
            pressure_from_velocity,
            velocity_from_pressure,
        })
    }

    /// Apply the same physical p/u map as the harmonic point-flow owner to
    /// both original Q and dt*Qdot columns, BEFORE quadrature norm acceptance.
    pub(super) fn apply(&self, values: [f64; 8]) -> Result<[f64; 8], AeroacError> {
        let mut result = [0.0; 8];
        for basis in 0..2 {
            let p = values[basis];
            let velocity = &values[2 + 3 * basis..5 + 3 * basis];
            result[basis] = (Dd::from_f64(p) * Dd::from_f64(self.pressure_gain)
                - Dd::from_f64(velocity[2]) * Dd::from_f64(self.pressure_from_velocity))
            .to_f64();
            result[2 + 3 * basis] = velocity[0] * self.inverse_beta;
            result[3 + 3 * basis] = velocity[1] * self.inverse_beta;
            result[4 + 3 * basis] = (Dd::from_f64(velocity[2])
                * Dd::from_f64(self.inverse_beta_squared)
                - Dd::from_f64(p) * Dd::from_f64(self.velocity_from_pressure))
            .to_f64();
        }
        if result.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "causal convected pressure or velocity is unrepresentable",
            ));
        }
        Ok(result)
    }
}
