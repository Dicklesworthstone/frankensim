//! Finite-length, inviscid plane sound in a smooth isentropic nozzle or duct.
//!
//! The mean owner retains the original gas and geometry and resolves the
//! subsonic branch of constant mass flow, entropy and stagnation enthalpy.
//! With no entropy perturbation, linearized mass conservation and Bernoulli
//! give `m'_x = i omega A p'/c²` and `H'_x = i omega u'`. Consequently,
//! in the crate's peak-phasor exp(-i omega t) convention,
//!
//! ```text
//! d[H',m']/dx = i omega/(c²-U²) [[-U, c²/(A rho)], [A rho, -U]] [H',m'].
//! ```
//!
//! The scalar convective phase commutes with the traceless generator. The
//! existing fourth-order two-node Magnus exponential integrates the latter
//! in scaled conjugate variables [H',R0 m']. Both preserve the physical flux
//! `Re(H' conj(m'))/2` in exact arithmetic. No finite-difference mean gradient,
//! compact-junction approximation or stationary p/Q energy port is introduced.
//! See Duran & Moreau (2013), <https://doi.org/10.1017/jfm.2013.118>, and the
//! no-heat/no-entropy specialization of Yeddula et al. (2022), equations 2.3–2.11,
//! <https://doi.org/10.1017/jfm.2022.101>. Their +i omega t signs are reversed here.
//!
//! This remains a quasi-one-dimensional, small-signal, calorically-perfect
//! model. The caller must justify the slender smooth geometry and original
//! gas's frozen-Cp range. Choking, supersonic branches, shocks, separation,
//! entropy/vorticity/composition waves, wall/bulk loss and radiation are not
//! modeled. A conservative whole-profile transverse cut-on check is necessary,
//! not a certificate of three-dimensional accuracy. Refining the explicit
//! integration axis controls numerical error within this model only.

mod mean;
pub use mean::{IsentropicMeanFlow, IsentropicMeanState, MINIMUM_SUBSONIC_MARGIN};
/// Driven original conjugate loads on the variable-flow profile.
pub mod termination;

use super::{bad, finite, phase, positive};
use crate::profile::{MAX_RICCATI_STEPS, RiccatiSteps, magnus_step};
use crate::{DuctError, mul2};
use fs_math::{c64::C64, det};

#[derive(Debug)]
struct Node {
    x: f64,
    mean: IsentropicMeanState,
    // After preparation, maps incoming sqrt(W) waves to physical [H',m'].
    response: [C64; 4],
    pressure_scales: [f64; 2],
    velocity_scale: f64,
    energy_scales: [f64; 2],
}

impl Node {
    fn new(x: f64, mean: IsentropicMeanState, response: [C64; 4]) -> Result<Self, DuctError> {
        let pressure_scales =
            [1.0 + mean.mach, 1.0 - mean.mach].map(|factor| 0.5 * mean.density_kg_m3 / factor);
        let velocity_scale = (1.0 / mean.density_kg_m3) / mean.sound_speed_m_s;
        let energy_scales = [1.0 + mean.mach, 1.0 - mean.mach]
            .map(|factor| det::sqrt((0.5 / mean.density_kg_m3) * factor) / mean.sound_speed_m_s);
        if !pressure_scales
            .into_iter()
            .chain(energy_scales)
            .chain([velocity_scale])
            .all(positive)
        {
            return Err(bad(
                "isentropic local observation coefficients are unrepresentable",
            ));
        }
        Ok(Self {
            x,
            mean,
            response,
            pressure_scales,
            velocity_scale,
            energy_scales,
        })
    }
}

/// Frequency-specific scattering and interior field observations.
///
/// Preparation costs O(steps) time/storage; every coherent source evaluation
/// and retained-node observation uses fixed work, no allocation or mean solve.
/// The requested steps apply to EACH spline span and the total is capped by
/// MAX_RICCATI_STEPS. No original spline knot is skipped or decimated.
#[derive(Debug)]
pub struct PreparedIsentropicProfile<'flow, 'source> {
    mean: &'flow IsentropicMeanFlow<'source>,
    omega: f64,
    scattering: [C64; 4],
    nodes: Vec<Node>,
    reference_impedances: [f64; 2],
}

/// Physical local fields, excluding the maintained steady state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowProfileField {
    /// Local axial coordinate [m].
    pub x_m: f64,
    /// Original isentropic background at this section.
    pub mean: IsentropicMeanState,
    /// Acoustic pressure [Pa], peak phasor.
    pub pressure_pa: C64,
    /// Axial acoustic particle velocity [m/s], peak phasor.
    pub particle_velocity_m_s: C64,
    /// Axial acoustic volume flow [m³/s], peak phasor.
    pub volume_flow_m3_s: C64,
    /// Total axial mass-flow perturbation [kg/s], peak phasor.
    pub mass_flow_perturbation_kg_s: C64,
    /// Stagnation-enthalpy perturbation [J/kg], peak phasor.
    pub stagnation_enthalpy_perturbation_j_kg: C64,
    /// Cycle-mean stationary pressure work [W], distinct from physical flux.
    pub pressure_work_power_w: f64,
    /// Cycle-mean Myers energy including the mean-flow cross term [J/m³].
    pub disturbance_energy_density_j_m3: f64,
    /// Cycle-mean Re(H' conj(m'))/2 along +x [W].
    pub disturbance_power_w: f64,
}

/// Coherent two-sided incoming characteristics on the prepared profile.
#[derive(Debug)]
pub struct IsentropicProfileExcitation<'prepared, 'flow, 'source> {
    prepared: &'prepared PreparedIsentropicProfile<'flow, 'source>,
    incoming: [C64; 2],
    outgoing: [C64; 2],
    incident_power: f64,
    outgoing_power: f64,
}

impl<'source> IsentropicMeanFlow<'source> {
    /// Prepare the complete finite-length propagator with explicit resolution.
    /// Refine `steps` to establish convergence for the actual frequency/shape;
    /// neither the step count nor a small flux residual bounds phase error.
    pub fn prepare(
        &self,
        omega_rad_s: f64,
        steps: RiccatiSteps,
    ) -> Result<PreparedIsentropicProfile<'_, 'source>, DuctError> {
        self.prepare_controlled(omega_rad_s, steps, &mut || false)
    }

    /// As `prepare`, polling before every integration panel. Cancellation
    /// returns no partial propagator. Requests exceeding the cap refuse.
    pub fn prepare_controlled(
        &self,
        omega_rad_s: f64,
        steps: RiccatiSteps,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<PreparedIsentropicProfile<'_, 'source>, DuctError> {
        steps.validate()?;
        if !positive(omega_rad_s) || omega_rad_s >= self.transverse_cutoff_bound_rad_s() {
            return Err(bad(
                "isentropic frequency must be positive and below whole-profile transverse cut-on",
            ));
        }
        let spans = self.profile.stations().map_or(1, |s| s.len() - 1);
        let count = steps
            .0
            .checked_mul(spans)
            .filter(|&n| n <= MAX_RICCATI_STEPS)
            .ok_or_else(|| bad("isentropic profile exceeds total integration step cap"))?;
        let inlet = self.at_position(0.0)?;
        let outlet = self.at_position(self.profile.length())?;
        let r0 = inlet.flow_impedance_j_s_kg2;
        let r1 = outlet.flow_impedance_j_s_kg2;
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(count + 1)
            .map_err(|_| bad("isentropic profile observation storage unavailable"))?;
        let mut transfer = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
        nodes.push(Node::new(0.0, inlet, transfer)?);
        let root3 = det::sqrt(3.0);
        let c1 = 0.5 - root3 / 6.0;
        let c2 = 0.5 + root3 / 6.0;
        let mut total_phase = C64::ONE;
        for span in 0..spans {
            let (start, end) = self
                .profile
                .stations()
                .map_or((0.0, self.profile.length()), |s| {
                    (s[span].0 - s[0].0, s[span + 1].0 - s[0].0)
                });
            for panel in 0..steps.0 {
                if cancel() {
                    return Err(DuctError::Cancelled {
                        completed_steps: nodes.len() - 1,
                    });
                }
                let left = start + (end - start) * (panel as f64 / steps.0 as f64);
                let right = if panel + 1 == steps.0 {
                    end
                } else {
                    start + (end - start) * ((panel + 1) as f64 / steps.0 as f64)
                };
                let h = right - left;
                if !positive(h) {
                    return Err(bad("isentropic integration panel is unrepresentable"));
                }
                let a = self.at_position(left + c1 * h)?;
                let b = self.at_position(left + c2 * h)?;
                let (ga, da) = generator(a, omega_rad_s, r0)?;
                let (gb, db) = generator(b, omega_rad_s, r0)?;
                let angle = h * (0.5 * da + 0.5 * db);
                if !angle.is_finite() {
                    return Err(bad("isentropic convective phase is unrepresentable"));
                }
                let scalar = phase(angle);
                let step = magnus_step(ga, gb, h, root3)?.map(|z| z * scalar);
                transfer = mul2(step, transfer);
                total_phase = total_phase * scalar;
                if !transfer.into_iter().all(finite) {
                    return Err(bad(
                        "isentropic transfer is unrepresentable; refine or repartition the profile",
                    ));
                }
                nodes.push(Node::new(right, self.at_position(right)?, transfer)?);
            }
        }
        let ratio = r1 / r0;
        let scale = 0.5 * det::sqrt(r0 / r1);
        if !positive(ratio) || !positive(scale) {
            return Err(bad("isentropic end normalization is unrepresentable"));
        }
        // Map [incoming-left,outgoing-left] to [outgoing-right,incoming-right]
        // in sqrt(W) waves. Each local H characteristic is sqrt(2R)*wave.
        let t = transfer;
        let b12 = (t[0] - t[1] + (t[2] - t[3]).scale(ratio)).scale(scale);
        let b21 = (t[0] + t[1] - (t[2] + t[3]).scale(ratio)).scale(scale);
        let b22 = (t[0] - t[1] - (t[2] - t[3]).scale(ratio)).scale(scale);
        if !finite(b22) || b22.abs() == 0.0 {
            return Err(DuctError::Singular);
        }
        // det(T)=exp(-2i integral omega U/(c²-U²) dx). Using this identity
        // avoids subtracting two potentially large products for transmission.
        let scattering = [
            -b21 / b22,
            C64::ONE / b22,
            total_phase * total_phase / b22,
            b12 / b22,
        ];
        let h_scale = det::sqrt(r0) * core::f64::consts::SQRT_2;
        let inlet_response = [
            (C64::ONE + scattering[0]).scale(h_scale),
            scattering[1].scale(h_scale),
            (C64::ONE - scattering[0]).scale(h_scale),
            -scattering[1].scale(h_scale),
        ];
        if !scattering.into_iter().all(finite) || !positive(h_scale) {
            return Err(bad("isentropic scattering is unrepresentable"));
        }
        for node in &mut nodes {
            node.response = mul2(node.response, inlet_response);
            node.response[2] = node.response[2].scale(1.0 / r0);
            node.response[3] = node.response[3].scale(1.0 / r0);
            if !node.response.into_iter().all(finite) {
                return Err(bad("isentropic field response is unrepresentable"));
            }
        }
        Ok(PreparedIsentropicProfile {
            mean: self,
            omega: omega_rad_s,
            scattering,
            nodes,
            reference_impedances: [r0, r1],
        })
    }
}

fn generator(
    state: IsentropicMeanState,
    omega: f64,
    r0: f64,
) -> Result<((C64, C64), f64), DuctError> {
    let d = (1.0 - state.mach) * (1.0 + state.mach);
    let k = (omega / state.sound_speed_m_s) / d;
    let ratio = state.flow_impedance_j_s_kg2 / r0;
    let pair = (C64::new(0.0, k * ratio), C64::new(0.0, k / ratio));
    if !positive(k) || !positive(ratio) || !finite(pair.0) || !finite(pair.1) {
        return Err(bad("isentropic acoustic generator is unrepresentable"));
    }
    Ok((pair, -k * state.mach))
}

impl<'flow, 'source> PreparedIsentropicProfile<'flow, 'source> {
    /// Original background, including original geometry and material owner.
    pub const fn mean_flow(&self) -> &'flow IsentropicMeanFlow<'source> {
        self.mean
    }
    /// Authored angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    /// Actual integration panels, including every spline span.
    pub fn steps_used(&self) -> usize {
        self.nodes.len() - 1
    }
    /// Cached observation sections, including both ends and every spline knot.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Axial coordinate [m] of a retained field observation.
    pub fn node_position_m(&self, index: usize) -> Option<f64> {
        self.nodes.get(index).map(|n| n.x)
    }
    /// Local H'/m' references [left,right], in J s/kg².
    pub const fn reference_impedances_j_s_kg2(&self) -> [f64; 2] {
        self.reference_impedances
    }
    /// Row-major scattering from incoming to outgoing sqrt(W) characteristics,
    /// ordered [left,right]. Squared modulus is physical disturbance power.
    /// Mean-flow reversal exchanges directional phases; reciprocity is not assumed.
    pub const fn scattering(&self) -> [C64; 4] {
        self.scattering
    }

    /// Prescribe incoming peak pressure waves [Pa] at x=0 and x=L.
    pub fn prepare_excitation(
        &self,
        pressure_pa: [C64; 2],
    ) -> Result<IsentropicProfileExcitation<'_, 'flow, 'source>, DuctError> {
        let ends = [self.nodes[0].mean, self.nodes[self.nodes.len() - 1].mean];
        let incoming = core::array::from_fn(|i| {
            let direction = if i == 0 { 1.0 } else { -1.0 };
            pressure_pa[i].scale(
                ((1.0 + direction * ends[i].mach) / ends[i].density_kg_m3)
                    / (det::sqrt(self.reference_impedances[i]) * core::f64::consts::SQRT_2),
            )
        });
        if !pressure_pa.into_iter().all(finite) {
            return Err(bad("isentropic incoming pressures must be finite"));
        }
        self.excite_power_waves(incoming)
    }

    /// Prescribe incoming power-normalized waves [sqrt(W)] at the two ends.
    /// The original source support/small-signal limits still apply.
    pub fn excite_power_waves(
        &self,
        incoming: [C64; 2],
    ) -> Result<IsentropicProfileExcitation<'_, 'flow, 'source>, DuctError> {
        let s = self.scattering;
        let outgoing = [
            s[0] * incoming[0] + s[1] * incoming[1],
            s[2] * incoming[0] + s[3] * incoming[1],
        ];
        let incident_power = incoming[0].norm_sq() + incoming[1].norm_sq();
        let outgoing_power = outgoing[0].norm_sq() + outgoing[1].norm_sq();
        if !incoming.into_iter().chain(outgoing).all(finite)
            || !incident_power.is_finite()
            || !outgoing_power.is_finite()
        {
            return Err(bad(
                "isentropic excitation or boundary power is unrepresentable",
            ));
        }
        Ok(IsentropicProfileExcitation {
            prepared: self,
            incoming,
            outgoing,
            incident_power,
            outgoing_power,
        })
    }
}

impl<'prepared, 'flow, 'source> IsentropicProfileExcitation<'prepared, 'flow, 'source> {
    /// Original prepared propagator and its source-bound background.
    pub const fn prepared(&self) -> &'prepared PreparedIsentropicProfile<'flow, 'source> {
        self.prepared
    }
    /// Incoming sqrt(W) waves [left,right].
    pub const fn incoming_power_waves(&self) -> [C64; 2] {
        self.incoming
    }
    /// Outgoing sqrt(W) waves [left,right].
    pub const fn outgoing_power_waves(&self) -> [C64; 2] {
        self.outgoing
    }
    /// Total incident physical disturbance power [W].
    pub const fn incident_disturbance_power_w(&self) -> f64 {
        self.incident_power
    }
    /// Total outgoing physical disturbance power [W].
    pub const fn outgoing_disturbance_power_w(&self) -> f64 {
        self.outgoing_power
    }
    /// Signed numerical incoming-minus-outgoing residual [W], not absorbed heat
    /// or a bound on integration/model error.
    pub fn conservation_residual_w(&self) -> f64 {
        self.incident_power - self.outgoing_power
    }

    /// Fixed-work local physical field at a retained section. No interpolation,
    /// allocation, material query or numerical integration occurs here.
    pub fn at_node(&self, index: usize) -> Result<FlowProfileField, DuctError> {
        let n = self
            .prepared
            .nodes
            .get(index)
            .ok_or_else(|| bad("isentropic observation index outside prepared nodes"))?;
        let t = n.response;
        let enthalpy = t[0] * self.incoming[0] + t[1] * self.incoming[1];
        let mass = t[2] * self.incoming[0] + t[3] * self.incoming[1];
        let s = n.mean;
        let forward = (enthalpy + mass.scale(s.flow_impedance_j_s_kg2)).scale(n.pressure_scales[0]);
        let backward =
            (enthalpy - mass.scale(s.flow_impedance_j_s_kg2)).scale(n.pressure_scales[1]);
        let pressure = forward + backward;
        let velocity = (forward - backward).scale(n.velocity_scale);
        let volume = velocity.scale(s.area_m2);
        let energy = forward.scale(n.energy_scales[0]).norm_sq()
            + backward.scale(n.energy_scales[1]).norm_sq();
        let pressure_work = 0.5 * (pressure * volume.conj()).re;
        let power = 0.5 * (enthalpy * mass.conj()).re;
        if ![pressure, velocity, volume, enthalpy, mass]
            .into_iter()
            .all(finite)
            || ![energy, pressure_work, power]
                .into_iter()
                .all(f64::is_finite)
        {
            return Err(bad(
                "isentropic local acoustic field or work is unrepresentable",
            ));
        }
        Ok(FlowProfileField {
            x_m: n.x,
            mean: s,
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            volume_flow_m3_s: volume,
            mass_flow_perturbation_kg_s: mass,
            stagnation_enthalpy_perturbation_j_kg: enthalpy,
            pressure_work_power_w: pressure_work,
            disturbance_energy_density_j_m3: energy,
            disturbance_power_w: power,
        })
    }
}

#[cfg(test)]
mod tests;
