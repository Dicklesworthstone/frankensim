//! Two acoustic boundary laws plus the independently prescribed upstream
//! entropy characteristic. The laws use TOTAL H' and TOTAL outward m'.
//! The H'/m' product is reported as an algebraic boundary pairing; with entropy
//! it is not relabelled as the full Myers flux or an acoustic passivity proof.

use super::{EntropyProfileExcitation, PreparedEntropyProfile};
use crate::DuctError;
use crate::convected::termination::FlowImpedance;
use crate::convected::{bad, finite};
use fs_math::{c64::C64, det};

#[derive(Debug)]
pub struct PreparedEntropyTerminations<'prepared, 'acoustic, 'flow, 'source> {
    profile: &'prepared PreparedEntropyProfile<'acoustic, 'flow, 'source>,
    impedances: [FlowImpedance; 2],
    // Row-major 2x3: two physical source enthalpies and upstream sigma.
    source_to_wave: [C64; 6],
    condition: f64,
}

/// Actual physical boundary variables and the authored law residual.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EntropyTerminalObservation {
    pub stagnation_enthalpy_perturbation_j_kg: C64,
    pub outward_mass_flow_perturbation_kg_s: C64,
    pub source_enthalpy_j_kg: C64,
    pub relative_entropy: C64,
    pub boundary_residual_j_kg: C64,
    /// Re(H' conj(m_out'))/2 [W]; not total entropy/acoustic energy flux.
    pub enthalpy_mass_pair_power_w: f64,
    /// Re(Z)|m_out'|²/2 [W]; resistance in the authored pairing.
    pub resistance_pair_power_w: f64,
    /// -Re(H_drive' conj(m_out'))/2 [W], in the same algebraic pairing.
    pub source_pair_work_into_duct_w: f64,
}

#[derive(Debug)]
pub struct TerminatedEntropyProfile<'prepared, 'acoustic, 'flow, 'source> {
    field: EntropyProfileExcitation<'prepared, 'acoustic, 'flow, 'source>,
    terminals: [EntropyTerminalObservation; 2],
}

impl<'acoustic, 'flow, 'source> PreparedEntropyProfile<'acoustic, 'flow, 'source> {
    /// Reuse the original two-ended acoustic feedback solve and its explicit
    /// condition ceiling. Each authored boundary obeys
    /// `H' = Z_Hm m_out' + H_drive'`, with entropy included in H' and m'.
    ///
    /// The independently prescribed entropy enters at the physical upstream
    /// end, even for negative flow. A nonzero entropy source changes the
    /// boundary forcing as well as the interior propagation. No stationary
    /// pressure/volume-flow load, supply model or thermal reservoir is inferred.
    pub fn prepare_terminations(
        &self,
        impedances: [FlowImpedance; 2],
        maximum_condition: f64,
    ) -> Result<PreparedEntropyTerminations<'_, 'acoustic, 'flow, 'source>, DuctError> {
        let original = self
            .acoustic
            .prepare_terminations(impedances, maximum_condition)?;
        let left = original
            .solve([C64::ONE, C64::ZERO])?
            .field()
            .incoming_power_waves();
        let right = original
            .solve([C64::ZERO, C64::ONE])?
            .field()
            .incoming_power_waves();
        let ends = [0, self.nodes.len() - 1];
        let mass = self.acoustic.mean.mass_flow_kg_s();
        let cp = self.acoustic.mean.source().specific_heat_cp;
        let mut effective_drive = [C64::ZERO; 2];
        for side in 0..2 {
            let index = ends[side];
            let mean = self.acoustic.nodes[index].mean;
            let z = impedances[side].value_j_s_kg2();
            let r = self.acoustic.reference_impedances[side];
            let normal = if side == 0 { -1.0 } else { 1.0 };
            let h_scale = det::sqrt(r) * core::f64::consts::SQRT_2;
            // r_reflection / g_source = sqrt(2R)*(Z/R-1).
            let returned =
                (z.scale(1.0 / r) - C64::ONE).scale(h_scale) * self.entropy_to_outgoing[side];
            let physical_entropy_effort =
                C64::from_re(cp * mean.temperature_k) + z.scale(normal * mass);
            effective_drive[side] = returned - physical_entropy_effort * self.nodes[index].sigma;
        }
        let entropy = [
            left[0] * effective_drive[0] + right[0] * effective_drive[1],
            left[1] * effective_drive[0] + right[1] * effective_drive[1],
        ];
        let source_to_wave = [left[0], right[0], entropy[0], left[1], right[1], entropy[1]];
        if !effective_drive
            .into_iter()
            .chain(source_to_wave)
            .all(finite)
        {
            return Err(bad("entropy boundary forcing is unrepresentable"));
        }
        Ok(PreparedEntropyTerminations {
            profile: self,
            impedances,
            source_to_wave,
            condition: original.condition_estimate(),
        })
    }
}

impl<'prepared, 'acoustic, 'flow, 'source>
    PreparedEntropyTerminations<'prepared, 'acoustic, 'flow, 'source>
{
    pub const fn profile(&self) -> &'prepared PreparedEntropyProfile<'acoustic, 'flow, 'source> {
        self.profile
    }
    pub const fn impedances(&self) -> [FlowImpedance; 2] {
        self.impedances
    }
    pub const fn condition_estimate(&self) -> f64 {
        self.condition
    }

    /// Fixed-work coherent solution with no allocation or integration.
    pub fn solve(
        &self,
        source_enthalpy_j_kg: [C64; 2],
        upstream_relative_entropy: C64,
    ) -> Result<TerminatedEntropyProfile<'prepared, 'acoustic, 'flow, 'source>, DuctError> {
        let h = source_enthalpy_j_kg;
        let sigma = upstream_relative_entropy;
        if !h.into_iter().chain([sigma]).all(finite) {
            return Err(bad("entropy boundary drives must be finite"));
        }
        let t = self.source_to_wave;
        let incoming = [
            t[0] * h[0] + t[1] * h[1] + t[2] * sigma,
            t[3] * h[0] + t[4] * h[1] + t[5] * sigma,
        ];
        let field = self.profile.excite_power_waves(incoming, sigma)?;
        let sections = [
            field.at_node(0)?,
            field.at_node(self.profile.node_count() - 1)?,
        ];
        let terminals = core::array::from_fn(|side| {
            let enthalpy = sections[side].stagnation_enthalpy_perturbation_j_kg;
            let mass = sections[side]
                .mass_flow_perturbation_kg_s
                .scale(if side == 0 { -1.0 } else { 1.0 });
            let z = self.impedances[side].value_j_s_kg2();
            EntropyTerminalObservation {
                stagnation_enthalpy_perturbation_j_kg: enthalpy,
                outward_mass_flow_perturbation_kg_s: mass,
                source_enthalpy_j_kg: h[side],
                relative_entropy: sections[side].relative_entropy,
                boundary_residual_j_kg: enthalpy - z * mass - h[side],
                enthalpy_mass_pair_power_w: 0.5 * (enthalpy * mass.conj()).re,
                resistance_pair_power_w: mass
                    .scale(det::sqrt(z.re) * core::f64::consts::FRAC_1_SQRT_2)
                    .norm_sq(),
                source_pair_work_into_duct_w: -0.5 * (h[side] * mass.conj()).re,
            }
        });
        if terminals.iter().any(|end| {
            !finite(end.boundary_residual_j_kg)
                || ![
                    end.enthalpy_mass_pair_power_w,
                    end.resistance_pair_power_w,
                    end.source_pair_work_into_duct_w,
                ]
                .into_iter()
                .all(f64::is_finite)
        }) {
            return Err(bad("entropy terminal law or pairing is unrepresentable"));
        }
        Ok(TerminatedEntropyProfile { field, terminals })
    }
}

impl<'prepared, 'acoustic, 'flow, 'source>
    TerminatedEntropyProfile<'prepared, 'acoustic, 'flow, 'source>
{
    pub const fn field(&self) -> &EntropyProfileExcitation<'prepared, 'acoustic, 'flow, 'source> {
        &self.field
    }
    pub const fn terminals(&self) -> &[EntropyTerminalObservation; 2] {
        &self.terminals
    }
}
