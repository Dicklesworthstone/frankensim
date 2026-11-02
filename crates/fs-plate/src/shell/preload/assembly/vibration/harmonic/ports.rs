//! Frequency-dependent impedance on the accepted assembly's physical ports.
//!
//! A real unit-effort force shape W defines BOTH the physical reaction -W e
//! and flow W^T v. The original support/tie/modal transpose projects B once.
//! The shared harmonic owner then solves M,C,K with B Z B^T before factoring.
//! No bare-body inverse, damping approximation or new radiation law is used.

use super::super::checkpoint;
use super::{AssemblyHarmonicModel, AssemblyHarmonicResponse, PreparedAssemblyHarmonicResponse};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicPortControls, HarmonicResponseControls, PreparedHarmonicResponse,
};

/// One immutable physical attachment, in original component/DOF order.
#[derive(Clone, Debug)]
pub struct AssemblyHarmonicPort {
    /// Force/moment per unit resisting effort on each physical component.
    /// The transpose acting on velocity defines the conjugate flow. For a
    /// pressure/volume-flow port these weights include the physical area;
    /// for a force/velocity port they are dimensionless. Units and physical
    /// ownership must agree with the supplied impedance, with no duplicate loss.
    pub component_force_shapes: Vec<Vec<f64>>,
}

/// Cold physical-map and one additional physical-response payload limits.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyPortProjectionBudget {
    /// Complete physical port count; dark attachments still count.
    pub maximum_ports: usize,
    /// p*(2*physical_dofs + slave_terms + modes*free_dofs + modes).
    /// Counts original support/tie work and complete modal projections.
    pub maximum_projection_terms: u64,
    /// Original owned shapes, B, projection scratch, and one reaction/flow
    /// response. The body, base harmonic response and shared factor are separate.
    pub maximum_payload_bytes: usize,
}

/// Cached physical port map on one immutable accepted harmonic assembly.
pub struct AssemblyHarmonicPorts<'m, 'v, 'body> {
    model: &'m AssemblyHarmonicModel<'v, 'body>,
    ports: Vec<AssemblyHarmonicPort>,
    force_map: Vec<f64>,
}

/// One loaded frequency factor with the original physical attachment maps.
pub struct PreparedAssemblyPortHarmonicResponse<'p, 'm, 'v, 'body> {
    ports: &'p AssemblyHarmonicPorts<'m, 'v, 'body>,
    response: PreparedAssemblyHarmonicResponse<'m, 'v, 'body>,
    relative_tolerance: f64,
}

/// Loaded perturbation and independently reconstructed physical port work.
#[derive(Debug)]
pub struct AssemblyPortHarmonicResponse {
    /// Original accepted-body motion, source work and retained-coordinate solve.
    pub assembly: AssemblyHarmonicResponse,
    /// W^T v evaluated from full physical component velocities.
    pub physical_port_flows: Vec<C64>,
    /// Force/moment ON the structure, -W e, in original component/DOF order.
    /// Supported coordinates may carry a reaction while their motion is zero.
    pub reaction_forces: Vec<Vec<C64>>,
    /// -Re(reaction^H v)/2, positive when the external load absorbs work.
    /// This is boundary work; no storage or heat certificate is inferred.
    pub average_physical_load_power_w: f64,
    /// Physical applied work minus glue heat minus physical load work [W].
    pub physical_power_balance_residual_w: f64,
}

impl<'m, 'v, 'body> AssemblyHarmonicPorts<'m, 'v, 'body> {
    /// Consume the source attachment shapes without cloning or changing them.
    /// Every port is projected through the original Phi^T T^T exactly once.
    /// A fixed or modal-dark attachment is retained, never silently discarded.
    pub fn new(
        model: &'m AssemblyHarmonicModel<'v, 'body>,
        ports: Vec<AssemblyHarmonicPort>,
        budget: AssemblyPortProjectionBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, String> {
        checkpoint(cx)?;
        let (n, p) = (model.vibration.mode_count(), ports.len());
        if p == 0 || p > budget.maximum_ports {
            return Err("assembled harmonic ports require complete nonzero port capacity".into());
        }
        let equilibrium = model.vibration.equilibrium;
        let components = &equilibrium.displacements;
        let physical = components.iter().map(|v| v.len() as u128).sum::<u128>();
        let free = equilibrium.assembly.model.free as u128;
        let slave_terms = equilibrium
            .assembly
            .slave_couplings()
            .iter()
            .map(|s| s.master_terms.len() as u128)
            .sum::<u128>();
        let entries = n
            .checked_mul(p)
            .ok_or("assembled harmonic port map size overflow")?;
        let terms = p as u128 * (2 * physical + slave_terms + n as u128 * free + n as u128);
        let bytes = 8 * (p as u128 * physical + entries as u128 + free + n as u128)
            + 32 * (physical + p as u128)
            + 24 * ((p as u128 + 1) * components.len() as u128 + p as u128)
            + 1024;
        if terms > budget.maximum_projection_terms as u128
            || bytes > budget.maximum_payload_bytes as u128
        {
            return Err(
                "assembled harmonic physical ports exceed projection work or payload limits".into(),
            );
        }
        for port in &ports {
            checkpoint(cx)?;
            if port.component_force_shapes.len() != components.len()
                || port
                    .component_force_shapes
                    .iter()
                    .zip(components)
                    .any(|(row, shape)| {
                        row.len() != shape.len() || row.iter().any(|v| !v.is_finite())
                    })
            {
                return Err(
                    "assembled harmonic port shapes must match every finite physical component"
                        .into(),
                );
            }
        }
        let mut force_map = vec![0.0; entries];
        for (j, port) in ports.iter().enumerate() {
            checkpoint(cx)?;
            let source: Vec<_> = port
                .component_force_shapes
                .iter()
                .map(Vec::as_slice)
                .collect();
            let column = model.vibration.project_loads(&source, cx)?;
            if column.len() != n || column.iter().any(|v| !v.is_finite()) {
                return Err("assembled harmonic port projection is incomplete or nonfinite".into());
            }
            for (i, value) in column.into_iter().enumerate() {
                force_map[i * p + j] = value;
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            model,
            ports,
            force_map,
        })
    }

    /// Number of source attachments, preserving their original order.
    pub fn port_count(&self) -> usize {
        self.ports.len()
    }

    /// Original row-major B [retained coordinates x physical ports].
    pub fn force_map(&self) -> &[f64] {
        &self.force_map
    }

    /// Original immutable unit-effort force/moment shapes.
    pub fn physical_ports(&self) -> &[AssemblyHarmonicPort] {
        &self.ports
    }

    /// The exact accepted body and retained basis that own this physical map.
    /// Time-domain consumers use this identity, not just matching dimensions.
    pub fn vibration(&self) -> &'v super::AssemblyVibration<'body> {
        self.model.vibration
    }

    /// Load the complete source pencil with Z(omega) under exp(-i*omega*t).
    /// This matrix maps conjugate port flows to resisting efforts. It may have
    /// reactive, mutual and nonreciprocal terms. The common numerical Hermitian
    /// PSD admission is pointwise, not a broadband causality/passivity claim.
    pub fn prepare(
        &self,
        angular_frequency_rad_s: f64,
        impedance: &[C64],
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedAssemblyPortHarmonicResponse<'_, 'm, 'v, 'body>, String> {
        let response = PreparedHarmonicResponse::try_with_port_impedance(
            self.model.mass_matrix(),
            self.model.damping_matrix(),
            self.model.stiffness_matrix(),
            self.model.vibration.mode_count(),
            &self.force_map,
            self.port_count(),
            impedance,
            angular_frequency_rad_s,
            controls,
            port_controls,
            cx,
        )
        .map_err(|e| e.to_string())?;
        checkpoint(cx)?;
        Ok(PreparedAssemblyPortHarmonicResponse {
            ports: self,
            response: PreparedAssemblyHarmonicResponse {
                model: self.model,
                response,
                relative_tolerance: controls.maximum_relative_residual,
            },
            relative_tolerance: controls.maximum_relative_residual,
        })
    }
}

impl PreparedAssemblyPortHarmonicResponse<'_, '_, '_, '_> {
    /// Apply physical forces to the loaded system and audit its original ports.
    /// All reconstructed outputs are provisional until every check succeeds.
    pub fn solve(
        &self,
        component_loads: &[&[C64]],
        cx: &Cx<'_>,
    ) -> Result<AssemblyPortHarmonicResponse, String> {
        let assembly = self.response.solve(component_loads, cx)?;
        let load = assembly
            .modal
            .port_load
            .as_ref()
            .ok_or("assembled harmonic port response omitted its load")?;
        let p = self.ports.port_count();
        if load.flow.len() != p || load.effort.len() != p {
            return Err("assembled harmonic port response dimensions disagree".into());
        }
        let mut physical_port_flows = Vec::with_capacity(p);
        for (port, &expected) in self.ports.ports.iter().zip(&load.flow) {
            checkpoint(cx)?;
            let (mut flow, mut correction, mut scale) = (C64::ZERO, C64::ZERO, 0.0);
            for (row, velocity) in port.component_force_shapes.iter().zip(&assembly.velocities) {
                for (j, (&weight, &v)) in row.iter().zip(velocity).enumerate() {
                    if j % 1024 == 0 {
                        checkpoint(cx)?;
                    }
                    if weight == 0.0 {
                        continue;
                    }
                    add_complex(&mut flow, &mut correction, v.scale(weight));
                    scale += weight.abs() * v.abs();
                }
            }
            if !finite(flow)
                || !scale.is_finite()
                || (flow - expected).abs() > self.relative_tolerance * scale.max(f64::MIN_POSITIVE)
            {
                return Err("assembled harmonic physical/modal port flow disagrees".into());
            }
            physical_port_flows.push(flow);
        }
        let mut reaction_forces: Vec<Vec<_>> = assembly
            .velocities
            .iter()
            .map(|v| vec![C64::ZERO; v.len()])
            .collect();
        let (mut power, mut power_correction, mut scale, mut input_scale) = (0.0, 0.0, 0.0, 0.0);
        for (component, ((reaction, velocity), force)) in reaction_forces
            .iter_mut()
            .zip(&assembly.velocities)
            .zip(component_loads)
            .enumerate()
        {
            checkpoint(cx)?;
            for (j, value) in reaction.iter_mut().enumerate() {
                if j % 1024 == 0 {
                    checkpoint(cx)?;
                }
                let mut force_correction = C64::ZERO;
                let mut force_scale = 0.0;
                for (port, &effort) in self.ports.ports.iter().zip(&load.effort) {
                    let weight = port.component_force_shapes[component][j];
                    if weight == 0.0 {
                        continue;
                    }
                    add_complex(value, &mut force_correction, effort.scale(-weight));
                    force_scale += effort.abs() * weight.abs();
                }
                if !finite(*value) {
                    return Err("assembled harmonic physical reaction overflow".into());
                }
                let work = -0.5 * (value.conj() * velocity[j]).re;
                let adjusted = work - power_correction;
                let next = power + adjusted;
                power_correction = (next - power) - adjusted;
                power = next;
                scale += 0.5 * force_scale * velocity[j].abs();
                input_scale += 0.5 * force[j].abs() * velocity[j].abs();
            }
        }
        let balance =
            assembly.average_component_power_w - assembly.modal.average_viscous_power_w - power;
        let balance_scale = input_scale + assembly.modal.average_viscous_power_w.abs() + scale;
        if !power.is_finite()
            || !scale.is_finite()
            || !balance.is_finite()
            || !balance_scale.is_finite()
            || (power - load.average_absorbed_power_w).abs()
                > self.relative_tolerance * scale.max(f64::MIN_POSITIVE)
            || balance.abs() > self.relative_tolerance * balance_scale.max(f64::MIN_POSITIVE)
        {
            return Err("assembled harmonic physical port work or power balance disagrees".into());
        }
        checkpoint(cx)?;
        Ok(AssemblyPortHarmonicResponse {
            assembly,
            physical_port_flows,
            reaction_forces,
            average_physical_load_power_w: power,
            physical_power_balance_residual_w: balance,
        })
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

fn add_complex(sum: &mut C64, correction: &mut C64, value: C64) {
    let adjusted = value - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

#[cfg(test)]
mod tests;
