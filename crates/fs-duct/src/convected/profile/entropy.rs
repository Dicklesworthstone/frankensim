//! Advected entropy and noncompact entropy-to-sound conversion.
//!
//! Let sigma=s'/cp, B=p'/rho+U u', and j=A(rho u'+U p'/c²). For the
//! ORIGINAL steady homentropic mean flow, linearized Euler gives
//!
//! ```text
//! sigma_x = i omega sigma/U,
//! B_x = i omega u' + U U_x sigma,    j_x = i omega A p'/c²,
//! H' = B + cp T sigma,              m' = j - m_bar sigma.
//! ```
//!
//! This is the same acoustic generator as the parent owner, forced by the
//! original entropy characteristic. The area/Mach law gives the analytic
//! gradient `U_x=-U (A_x/A)/(1-M²)`; no differenced mean state is introduced.
//! These relations follow from Yeddula et al. (2022), equations 2.1–2.8,
//! <https://doi.org/10.1017/jfm.2022.101>, with zero heat transfer and the
//! present exp(-i omega t) convention. Uniform flow transports entropy without
//! emitting sound. Acceleration/deceleration permits signed acoustic exchange.
//!
//! The two-node fourth-order Magnus method now retains the third characteristic
//! and its coupling. Both propagation directions use the original mean flow;
//! entropy is prescribed only at the UPSTREAM end. The parent acoustic
//! scattering/field response is reused, including its normalization and bounds.
//! Zero mean flow refuses: a finite-frequency advected entropy boundary problem
//! does not have a stationary-flow continuation. The original acoustic owner
//! still handles zero flow. Entropy wavelengths can be much shorter than sound
//! wavelengths: the caller must separately refine transport phase and conversion.
//!
//! This is linear, inviscid, calorically perfect quasi-1D flow. It does not add
//! mean heat transfer, entropy generation, wall losses, shocks, separation,
//! vorticity/composition waves or a three-dimensional source qualification.

use super::{IsentropicMeanState, PreparedIsentropicProfile, generator};
use crate::DuctError;
use crate::convected::{bad, finite, positive};
use fs_math::{c64::C64, det};

mod magnus;
/// Original full stagnation-enthalpy / mass-flow boundary laws with entropy.
pub mod termination;

/// End that supplies the incoming entropy disturbance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntropyInlet {
    Left,
    Right,
}

#[derive(Debug)]
struct EntropyNode {
    // Physical B [J/kg], j [kg/s] caused by unit UPSTREAM sigma.
    acoustic: [C64; 2],
    sigma: C64,
    acceleration_m_s2: f64,
}

/// Three-characteristic observation and scattering on an original profile.
///
/// Borrows the complete acoustic cache; it does not clone geometry/material or
/// rebuild the acoustic scattering operator. Additional preparation/storage is
/// O(original integration panels), capped by the parent's MAX_RICCATI_STEPS.
/// Warm source, boundary and node observations have fixed work and no allocation.
#[derive(Debug)]
pub struct PreparedEntropyProfile<'acoustic, 'flow, 'source> {
    acoustic: &'acoustic PreparedIsentropicProfile<'flow, 'source>,
    nodes: Vec<EntropyNode>,
    entropy_to_outgoing: [C64; 2],
    inlet: EntropyInlet,
}

/// Complete physical perturbations on the maintained background.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EntropyProfileField {
    pub x_m: f64,
    pub mean: IsentropicMeanState,
    pub pressure_pa: C64,
    pub particle_velocity_m_s: C64,
    pub volume_flow_m3_s: C64,
    /// sigma=s'/cp, dimensionless peak phasor.
    pub relative_entropy: C64,
    pub entropy_perturbation_j_kg_k: C64,
    pub density_perturbation_kg_m3: C64,
    pub temperature_perturbation_k: C64,
    /// TOTAL A(rho u'+U rho'), including advected entropy.
    pub mass_flow_perturbation_kg_s: C64,
    /// TOTAL h'+Uu', including cp T sigma.
    pub stagnation_enthalpy_perturbation_j_kg: C64,
    /// Acoustic component B=p'/rho+Uu'; not total H' when sigma is nonzero.
    pub acoustic_bernoulli_j_kg: C64,
    /// Acoustic component j=A(rho u'+U p'/c²); not total mass perturbation.
    pub acoustic_mass_flow_kg_s: C64,
    /// Cycle-mean stationary pressure work Re(p' conj(Q'))/2.
    pub pressure_work_power_w: f64,
    /// Re(B conj(j))/2: acoustic component flux, NOT total entropy-wave power.
    pub acoustic_power_w: f64,
    /// d(acoustic power)/dx = U U_x Re(sigma conj(j))/2 [W/m].
    /// Signed conversion with the maintained flow, never dissipated heat.
    pub entropy_conversion_power_per_length_w_m: f64,
}

/// Coherent acoustic waves at both ends plus the original upstream entropy.
#[derive(Debug)]
pub struct EntropyProfileExcitation<'prepared, 'acoustic, 'flow, 'source> {
    prepared: &'prepared PreparedEntropyProfile<'acoustic, 'flow, 'source>,
    incoming: [C64; 2],
    outgoing: [C64; 2],
    entropy: C64,
    incident_acoustic_power: f64,
    outgoing_acoustic_power: f64,
}

impl<'flow, 'source> PreparedIsentropicProfile<'flow, 'source> {
    /// Prepare entropy transport and its acoustic forcing on every ORIGINAL
    /// integration panel. Refine the original RiccatiSteps to resolve both
    /// entropy phase and sound generation; plane-mode cut-on is not this test.
    pub fn prepare_entropy(&self) -> Result<PreparedEntropyProfile<'_, 'flow, 'source>, DuctError> {
        self.prepare_entropy_controlled(&mut || false)
    }

    /// Cancellation publishes no partial cache and leaves this owner unchanged.
    /// Each panel has two interior mean-state samples, analytic nodal geometry
    /// and bounded fixed-size exponential work. The parent's explicit panel
    /// count bounds this work.
    pub fn prepare_entropy_controlled(
        &self,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<PreparedEntropyProfile<'_, 'flow, 'source>, DuctError> {
        if cancel() {
            return Err(DuctError::Cancelled { completed_steps: 0 });
        }
        let mass = self.mean.mass_flow_kg_s();
        if mass == 0.0 {
            return Err(bad("advected entropy requires nonzero mean through-flow"));
        }
        let inlet = if mass > 0.0 {
            EntropyInlet::Left
        } else {
            EntropyInlet::Right
        };
        let r0 = self.reference_impedances[0];
        let c0 = self.nodes[0].mean.sound_speed_m_s;
        let enthalpy_scale = c0 * c0;
        if !positive(enthalpy_scale) || !positive(enthalpy_scale.recip()) {
            return Err(bad("entropy reference scale is unrepresentable"));
        }
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(self.nodes.len())
            .map_err(|_| bad("entropy observation storage unavailable"))?;
        nodes.push(EntropyNode {
            acoustic: [C64::ZERO; 2],
            sigma: C64::ONE,
            acceleration_m_s2: acceleration(self, 0.0, self.nodes[0].mean)?,
        });
        let root3 = det::sqrt(3.0);
        let gauss = [0.5 - root3 / 6.0, 0.5 + root3 / 6.0];
        let mut particular = [C64::ZERO; 2];
        let mut sigma = C64::ONE;
        for panel in 0..self.steps_used() {
            if cancel() {
                return Err(DuctError::Cancelled {
                    completed_steps: panel,
                });
            }
            let left = self.nodes[panel].x;
            let right = self.nodes[panel + 1].x;
            let h = right - left;
            let ga = entropy_generator(self, left + h * gauss[0], r0, enthalpy_scale)?;
            let gb = entropy_generator(self, left + h * gauss[1], r0, enthalpy_scale)?;
            let step = magnus::step(ga, gb, h, root3)?;
            particular = [
                step.acoustic[0] * particular[0]
                    + step.acoustic[1] * particular[1]
                    + step.forcing[0] * sigma,
                step.acoustic[2] * particular[0]
                    + step.acoustic[3] * particular[1]
                    + step.forcing[1] * sigma,
            ];
            sigma = sigma * step.entropy;
            let acoustic = [
                particular[0].scale(enthalpy_scale),
                particular[1].scale(enthalpy_scale / r0),
            ];
            if !acoustic.into_iter().chain([sigma]).all(finite) {
                return Err(bad(
                    "entropy propagation is unrepresentable; refine the profile",
                ));
            }
            nodes.push(EntropyNode {
                acoustic,
                sigma,
                acceleration_m_s2: acceleration(self, right, self.nodes[panel + 1].mean)?,
            });
        }
        // The particular solution starts with B=j=0 at x=0. Remove its
        // incoming acoustic wave at x=L with the ORIGINAL homogeneous solution.
        let last = nodes.len() - 1;
        let r1 = self.reference_impedances[1];
        let h1 = det::sqrt(r1) * core::f64::consts::SQRT_2;
        let correction =
            -(nodes[last].acoustic[0] - nodes[last].acoustic[1].scale(r1)).scale(0.5 / h1);
        let upstream_scale = match inlet {
            EntropyInlet::Left => C64::ONE,
            EntropyInlet::Right => sigma.recip(),
        };
        if !finite(upstream_scale) || !finite(correction) {
            return Err(bad(
                "entropy upstream/boundary normalization is unrepresentable",
            ));
        }
        for (index, (node, original)) in nodes.iter_mut().zip(&self.nodes).enumerate() {
            if index.is_multiple_of(128) && cancel() {
                return Err(DuctError::Cancelled {
                    completed_steps: self.steps_used(),
                });
            }
            node.acoustic[0] =
                (node.acoustic[0] + original.response[1] * correction) * upstream_scale;
            node.acoustic[1] =
                (node.acoustic[1] + original.response[3] * correction) * upstream_scale;
            node.sigma = node.sigma * upstream_scale;
            if !node.acoustic.into_iter().chain([node.sigma]).all(finite) {
                return Err(bad("entropy field response is unrepresentable"));
            }
        }
        let h0 = det::sqrt(r0) * core::f64::consts::SQRT_2;
        let entropy_to_outgoing = [
            (nodes[0].acoustic[0] - nodes[0].acoustic[1].scale(r0)).scale(0.5 / h0),
            (nodes[last].acoustic[0] + nodes[last].acoustic[1].scale(r1)).scale(0.5 / h1),
        ];
        if !entropy_to_outgoing.into_iter().all(finite) {
            return Err(bad("entropy outgoing waves are unrepresentable"));
        }
        if cancel() {
            return Err(DuctError::Cancelled {
                completed_steps: self.steps_used(),
            });
        }
        Ok(PreparedEntropyProfile {
            acoustic: self,
            nodes,
            entropy_to_outgoing,
            inlet,
        })
    }
}

fn acceleration(
    acoustic: &PreparedIsentropicProfile<'_, '_>,
    x: f64,
    mean: IsentropicMeanState,
) -> Result<f64, DuctError> {
    let geometry = acoustic.mean.radius_profile();
    let area_gradient = 2.0 * (geometry.slope(x)? / geometry.radius(x)?);
    let result = -mean.velocity_m_s * mean.velocity_m_s * area_gradient
        / ((1.0 - mean.mach) * (1.0 + mean.mach));
    if !result.is_finite() {
        return Err(bad("entropy mean acceleration is unrepresentable"));
    }
    Ok(result)
}

fn entropy_generator(
    acoustic: &PreparedIsentropicProfile<'_, '_>,
    x: f64,
    reference: f64,
    enthalpy_scale: f64,
) -> Result<magnus::Block, DuctError> {
    let mean = acoustic.mean.at_position(x)?;
    let ((a, b), d) = generator(mean, acoustic.omega, reference)?;
    let result = magnus::Block {
        acoustic: [C64::new(0.0, d), a, b, C64::new(0.0, d)],
        forcing: [
            C64::from_re(acceleration(acoustic, x, mean)? / enthalpy_scale),
            C64::ZERO,
        ],
        entropy: C64::new(0.0, acoustic.omega / mean.velocity_m_s),
    };
    if !result.is_finite() {
        return Err(bad("entropy transport generator is unrepresentable"));
    }
    Ok(result)
}

impl<'acoustic, 'flow, 'source> PreparedEntropyProfile<'acoustic, 'flow, 'source> {
    pub const fn acoustic_profile(&self) -> &'acoustic PreparedIsentropicProfile<'flow, 'source> {
        self.acoustic
    }
    pub const fn entropy_inlet(&self) -> EntropyInlet {
        self.inlet
    }
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }
    /// Acoustic sqrt(W) responses per unit dimensionless UPSTREAM s'/cp.
    pub const fn entropy_to_outgoing_power_waves(&self) -> [C64; 2] {
        self.entropy_to_outgoing
    }
    /// Upstream-to-downstream entropy transfer, with its separately integrated phase.
    pub fn entropy_transfer(&self) -> C64 {
        self.nodes[match self.inlet {
            EntropyInlet::Left => self.nodes.len() - 1,
            EntropyInlet::Right => 0,
        }]
        .sigma
    }
    /// Row-major [out-left sqrt(W), out-right sqrt(W), downstream sigma] from
    /// [in-left sqrt(W), in-right sqrt(W), upstream sigma]. Units are mixed;
    /// the three-channel matrix is NOT asserted to be unitary or passive.
    pub fn scattering(&self) -> [C64; 9] {
        let s = self.acoustic.scattering;
        [
            s[0],
            s[1],
            self.entropy_to_outgoing[0],
            s[2],
            s[3],
            self.entropy_to_outgoing[1],
            C64::ZERO,
            C64::ZERO,
            self.entropy_transfer(),
        ]
    }
    pub fn prepare_excitation(
        &self,
        incoming_pressure_pa: [C64; 2],
        upstream_relative_entropy: C64,
    ) -> Result<EntropyProfileExcitation<'_, 'acoustic, 'flow, 'source>, DuctError> {
        let incoming = self
            .acoustic
            .prepare_excitation(incoming_pressure_pa)?
            .incoming_power_waves();
        self.excite_power_waves(incoming, upstream_relative_entropy)
    }
    pub fn excite_power_waves(
        &self,
        incoming: [C64; 2],
        entropy: C64,
    ) -> Result<EntropyProfileExcitation<'_, 'acoustic, 'flow, 'source>, DuctError> {
        let s = self.acoustic.scattering;
        let outgoing = [
            s[0] * incoming[0] + s[1] * incoming[1] + self.entropy_to_outgoing[0] * entropy,
            s[2] * incoming[0] + s[3] * incoming[1] + self.entropy_to_outgoing[1] * entropy,
        ];
        let incident_acoustic_power = incoming[0].norm_sq() + incoming[1].norm_sq();
        let outgoing_acoustic_power = outgoing[0].norm_sq() + outgoing[1].norm_sq();
        if !incoming
            .into_iter()
            .chain(outgoing)
            .chain([entropy])
            .all(finite)
            || !incident_acoustic_power.is_finite()
            || !outgoing_acoustic_power.is_finite()
        {
            return Err(bad("entropy/acoustic excitation is unrepresentable"));
        }
        Ok(EntropyProfileExcitation {
            prepared: self,
            incoming,
            outgoing,
            entropy,
            incident_acoustic_power,
            outgoing_acoustic_power,
        })
    }
}

impl<'prepared, 'acoustic, 'flow, 'source>
    EntropyProfileExcitation<'prepared, 'acoustic, 'flow, 'source>
{
    pub const fn prepared(&self) -> &'prepared PreparedEntropyProfile<'acoustic, 'flow, 'source> {
        self.prepared
    }
    pub const fn incoming_power_waves(&self) -> [C64; 2] {
        self.incoming
    }
    pub const fn outgoing_power_waves(&self) -> [C64; 2] {
        self.outgoing
    }
    pub const fn upstream_relative_entropy(&self) -> C64 {
        self.entropy
    }
    pub const fn incident_acoustic_power_w(&self) -> f64 {
        self.incident_acoustic_power
    }
    pub const fn outgoing_acoustic_power_w(&self) -> f64 {
        self.outgoing_acoustic_power
    }
    /// Outgoing minus incoming ACOUSTIC component flux. Entropy conversion can
    /// supply or remove it; this is neither a conservation residual nor heat.
    pub fn net_acoustic_generation_w(&self) -> f64 {
        self.outgoing_acoustic_power - self.incident_acoustic_power
    }

    pub fn at_node(&self, index: usize) -> Result<EntropyProfileField, DuctError> {
        let original = self
            .prepared
            .acoustic
            .nodes
            .get(index)
            .ok_or_else(|| bad("entropy observation index outside original nodes"))?;
        let response = &self.prepared.nodes[index];
        let t = original.response;
        let bernoulli =
            t[0] * self.incoming[0] + t[1] * self.incoming[1] + response.acoustic[0] * self.entropy;
        let mass =
            t[2] * self.incoming[0] + t[3] * self.incoming[1] + response.acoustic[1] * self.entropy;
        let acoustic = original.field(bernoulli, mass)?;
        let sigma = response.sigma * self.entropy;
        let mean = original.mean;
        let gas = self.prepared.acoustic.mean.source();
        let entropy = sigma.scale(gas.specific_heat_cp);
        let density = acoustic
            .pressure_pa
            .scale(1.0 / mean.sound_speed_m_s / mean.sound_speed_m_s)
            - sigma.scale(mean.density_kg_m3);
        let temperature = (acoustic.pressure_pa.scale(
            (gas.gamma - 1.0) / mean.density_kg_m3 / mean.sound_speed_m_s / mean.sound_speed_m_s,
        ) + sigma)
            .scale(mean.temperature_k);
        let total_mass = mass - sigma.scale(self.prepared.acoustic.mean.mass_flow_kg_s());
        let total_enthalpy = bernoulli + sigma.scale(gas.specific_heat_cp * mean.temperature_k);
        let conversion = 0.5 * response.acceleration_m_s2 * (sigma * mass.conj()).re;
        if ![
            sigma,
            entropy,
            density,
            temperature,
            total_mass,
            total_enthalpy,
        ]
        .into_iter()
        .all(finite)
            || !conversion.is_finite()
        {
            return Err(bad(
                "entropy physical field or conversion work is unrepresentable",
            ));
        }
        Ok(EntropyProfileField {
            x_m: acoustic.x_m,
            mean,
            pressure_pa: acoustic.pressure_pa,
            particle_velocity_m_s: acoustic.particle_velocity_m_s,
            volume_flow_m3_s: acoustic.volume_flow_m3_s,
            relative_entropy: sigma,
            entropy_perturbation_j_kg_k: entropy,
            density_perturbation_kg_m3: density,
            temperature_perturbation_k: temperature,
            mass_flow_perturbation_kg_s: total_mass,
            stagnation_enthalpy_perturbation_j_kg: total_enthalpy,
            acoustic_bernoulli_j_kg: bernoulli,
            acoustic_mass_flow_kg_s: mass,
            pressure_work_power_w: acoustic.pressure_work_power_w,
            acoustic_power_w: acoustic.disturbance_power_w,
            entropy_conversion_power_per_length_w_m: conversion,
        })
    }
}

#[cfg(test)]
mod tests;
