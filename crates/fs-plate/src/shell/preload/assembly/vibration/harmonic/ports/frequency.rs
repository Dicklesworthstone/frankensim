//! Exact physical frequency tangents of the complete loaded assembly.

use super::{PreparedAssemblyPortHarmonicResponse, add_real, checkpoint};
use crate::shell::preload::assembly::vibration::harmonic::AssemblyHarmonicFrequencyResponse;
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::HarmonicFrequencyControls;
use fs_modal::harmonic::resonance::HarmonicTransferSample;

/// Original-port output/audit limits, in addition to the shared modal controls.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyPortFrequencyBudget {
    /// Two physical reaction fields and two port-flow vectors, including vector
    /// headers: 32*(physical_dofs+ports)+48*components+256 bytes. Assembled
    /// motion, the immutable projection and shared modal scratch are separate.
    pub maximum_payload_bytes: usize,
    /// 4*ports*physical_dofs+12*physical_dofs+8*ports multiply/add observations.
    /// Includes complete primal/tangent port maps and original physical work.
    pub maximum_observation_terms: u64,
}

/// The original loaded motion and physical reactions, with all frequency slopes.
#[derive(Debug)]
pub struct AssemblyPortHarmonicFrequencyResponse {
    /// Original support/tie motion, source work and same-factor modal tangent.
    pub assembly: AssemblyHarmonicFrequencyResponse,
    /// W^T v in the original conjugate physical port units.
    pub physical_port_flows: Vec<C64>,
    /// W^T dv/dw, with the original geometry and attachments fixed.
    pub physical_port_flow_derivatives: Vec<C64>,
    /// Force/moment on the original structure, -W e, including supported DOFs.
    pub reaction_forces: Vec<Vec<C64>>,
    /// -W de/dw, retaining both Z_w flow and Z flow_w in the resisting effort.
    pub reaction_force_derivatives: Vec<Vec<C64>>,
    /// Original physical boundary work -Re(reaction^H v)/2 [W].
    pub average_physical_load_power_w: f64,
    /// Independently reconstructed derivative of boundary work [W per rad/s].
    pub average_physical_load_power_derivative: f64,
    /// Physical source work minus glue heat and physical load work [W].
    pub physical_power_balance_residual_w: f64,
    /// Derivative of the same original physical balance [W per rad/s].
    pub physical_power_balance_residual_derivative: f64,
}

impl PreparedAssemblyPortHarmonicResponse<'_, '_, '_, '_> {
    /// Differentiate the actual loaded physical response with this frequency's
    /// existing factor. Source force/moment slopes use original component DOFs;
    /// None explicitly holds source phasors fixed. The complete finite row-major
    /// dZ/dw is REQUIRED, including mutual and nonreciprocal entries. An explicit
    /// zero matrix declares a locally constant load. This slope has no PSD test.
    ///
    /// Full projected M/K/glue C, source phase and load dispersion are retained.
    /// The physical load owner supplies its own constitutive derivative on the
    /// exp(-i*w*t) axis; no differencing, second factor or fitted load is created.
    /// Physical flows, reactions, source/load work and both work slopes are
    /// reconstructed on the original support/tie map before publication. The
    /// load account is boundary work, not a material-heat/storage certificate.
    /// Limits are checked before physical allocation; refusal/cancellation
    /// preserves the accepted body and the reusable factor.
    pub fn solve_with_frequency_derivative(
        &self,
        component_loads: &[&[C64]],
        component_load_derivatives: Option<&[&[C64]]>,
        impedance_derivative: &[C64],
        controls: HarmonicFrequencyControls,
        budget: AssemblyPortFrequencyBudget,
        cx: &Cx<'_>,
    ) -> Result<AssemblyPortHarmonicFrequencyResponse, String> {
        checkpoint(cx)?;
        let components = &self.ports.model.vibration.equilibrium.displacements;
        let physical = components.iter().map(|v| v.len() as u128).sum::<u128>();
        let p = self.ports.port_count() as u128;
        let bytes = 32 * (physical + p) + 48 * components.len() as u128 + 256;
        let terms = 4 * p * physical + 12 * physical + 8 * p;
        if bytes > budget.maximum_payload_bytes as u128
            || terms > budget.maximum_observation_terms as u128
        {
            return Err(
                "assembled harmonic frequency port outputs exceed payload or work limits".into(),
            );
        }
        let assembly = self.response.solve_frequency(
            component_loads,
            component_load_derivatives,
            Some(impedance_derivative),
            controls,
            cx,
        )?;
        let load = assembly
            .modal
            .response
            .port_load
            .as_ref()
            .ok_or("assembled harmonic frequency response omitted its original load")?;
        let tangent = assembly
            .modal
            .port_load
            .as_ref()
            .ok_or("assembled harmonic frequency response omitted its load tangent")?;
        let physical_port_flows = self.physical_flows(&assembly.velocities, &load.flow, cx)?;
        let physical_port_flow_derivatives =
            self.physical_flows(&assembly.velocity_derivatives, &tangent.flow_derivative, cx)?;
        let (mut power, mut power_correction, mut power_scale, mut input_scale) =
            (0.0, 0.0, 0.0, 0.0);
        let (mut slope, mut slope_correction, mut slope_scale, mut input_slope_scale) =
            (0.0, 0.0, 0.0, 0.0);
        let reaction_forces =
            self.physical_reactions(&load.effort, cx, |c, j, reaction, envelope| {
                let (v, dv) = (
                    assembly.velocities[c][j],
                    assembly.velocity_derivatives[c][j],
                );
                add_real(
                    &mut power,
                    &mut power_correction,
                    -0.5 * (reaction.conj() * v).re,
                );
                add_real(
                    &mut slope,
                    &mut slope_correction,
                    -0.5 * (reaction.conj() * dv).re,
                );
                power_scale += 0.5 * envelope * v.abs();
                slope_scale += 0.5 * envelope * dv.abs();
                let source = component_loads[c][j].abs();
                let source_slope =
                    component_load_derivatives.map_or(0.0, |loads| loads[c][j].abs());
                input_scale += 0.5 * source * v.abs();
                input_slope_scale += 0.5 * (source * dv.abs() + source_slope * v.abs());
            })?;
        let reaction_force_derivatives = self.physical_reactions(
            &tangent.effort_derivative,
            cx,
            |c, j, reaction_slope, envelope| {
                let v = assembly.velocities[c][j];
                add_real(
                    &mut slope,
                    &mut slope_correction,
                    -0.5 * (reaction_slope.conj() * v).re,
                );
                slope_scale += 0.5 * envelope * v.abs();
            },
        )?;
        let viscous = assembly.modal.response.average_viscous_power_w;
        let viscous_slope = assembly.modal.average_viscous_power_derivative;
        let balance = assembly.average_component_power_w - viscous - power;
        let balance_slope = assembly.average_component_power_derivative - viscous_slope - slope;
        self.check_work(
            power,
            load.average_absorbed_power_w,
            power_scale,
            balance,
            input_scale + viscous.abs() + power_scale,
        )?;
        self.check_work(
            slope,
            tangent.average_absorbed_power_derivative,
            slope_scale,
            balance_slope,
            input_slope_scale + viscous_slope.abs() + slope_scale,
        )?;
        checkpoint(cx)?;
        Ok(AssemblyPortHarmonicFrequencyResponse {
            assembly,
            physical_port_flows,
            physical_port_flow_derivatives,
            reaction_forces,
            reaction_force_derivatives,
            average_physical_load_power_w: power,
            average_physical_load_power_derivative: slope,
            physical_power_balance_residual_w: balance,
            physical_power_balance_residual_derivative: balance_slope,
        })
    }
}

impl AssemblyPortHarmonicFrequencyResponse {
    /// A physical attachment's actual flow/slope, ready for local resonance or
    /// antiresonance refinement. A displacement maximum need not maximize flow.
    pub fn port_flow_transfer(&self, port: usize) -> Result<HarmonicTransferSample, String> {
        let amplitude = self
            .physical_port_flows
            .get(port)
            .ok_or("assembled harmonic receiver port is out of range")?;
        let angular_derivative = self
            .physical_port_flow_derivatives
            .get(port)
            .ok_or("assembled harmonic receiver port tangent is missing")?;
        Ok(HarmonicTransferSample {
            amplitude: *amplitude,
            angular_derivative: *angular_derivative,
        })
    }

    /// A load's resisting effort/slope, including the original dispersive law.
    /// Its peak/notch may differ from either structural displacement or flow.
    pub fn port_effort_transfer(&self, port: usize) -> Result<HarmonicTransferSample, String> {
        let amplitude = self
            .assembly
            .modal
            .response
            .port_load
            .as_ref()
            .and_then(|load| load.effort.get(port))
            .ok_or("assembled harmonic receiver effort is missing or out of range")?;
        let angular_derivative = self
            .assembly
            .modal
            .port_load
            .as_ref()
            .and_then(|load| load.effort_derivative.get(port))
            .ok_or("assembled harmonic receiver effort tangent is missing or out of range")?;
        Ok(HarmonicTransferSample {
            amplitude: *amplitude,
            angular_derivative: *angular_derivative,
        })
    }
}
