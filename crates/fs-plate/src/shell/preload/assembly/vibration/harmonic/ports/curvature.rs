//! Original supported/tied motion and work curvature of a dispersively loaded assembly.

use super::{
    AssemblyPortHarmonicFrequencyResponse, PreparedAssemblyPortHarmonicResponse, add_real,
    checkpoint,
};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{HarmonicPowerCurvatureSample, HarmonicPowerSample};
use fs_modal::harmonic::{
    HarmonicFrequencyCurvatureControls, HarmonicPortFrequencyCurvatureResponse,
    harmonic_mean_power_curvature_with_cancel,
};

/// Complete physical projection/recovery limits, additional to native modal limits.
#[derive(Clone, Copy, Debug)]
pub struct AssemblyPortCurvatureBudget {
    /// 256*d + 128*(n+f+p) + 512*c + 1024 numerical bytes, where d is the
    /// physical coordinate count, n the retained count, f free coordinates,
    /// p original ports and c components. Includes six motion fields, three
    /// reaction fields, source projection scratch and all three port flows.
    pub maximum_payload_bytes: usize,
    /// 24*(2*d+t+n*f+n) + 12*p*d + 256*d + 64*p multiply/add observations,
    /// where t counts original slave-map terms. Covers all source projections,
    /// six original support/tie expansions, port recovery and physical work.
    pub maximum_observation_terms: u64,
}

/// Original physical response, slope and curvature at one admitted frequency.
///
/// The accepted equilibrium, supports, ties, attachment and retained basis are
/// fixed. No equilibrium offset enters any velocity or perturbation field.
#[derive(Debug)]
pub struct AssemblyPortHarmonicCurvatureResponse {
    pub frequency_response: AssemblyPortHarmonicFrequencyResponse,
    pub displacement_second_derivatives: Vec<Vec<C64>>,
    pub velocity_second_derivatives: Vec<Vec<C64>>,
    /// Complete retained-coordinate second derivatives, without another basis.
    pub modal_displacement_second_derivative: Vec<C64>,
    pub modal_velocity_second_derivative: Vec<C64>,
    pub modal_port_curvature: HarmonicPortFrequencyCurvatureResponse,
    pub physical_port_flow_second_derivatives: Vec<C64>,
    pub reaction_force_second_derivatives: Vec<Vec<C64>>,
    /// Independently reconstructed original source work [W per (rad/s)^2].
    pub average_component_power_second_derivative: f64,
    /// Original symmetric glue work; external load work remains separate.
    pub average_viscous_power_second_derivative: f64,
    /// Independently reconstructed signed boundary work [W per (rad/s)^2].
    pub average_physical_load_power_second_derivative: f64,
    pub physical_power_balance_residual_second_derivative: f64,
    /// Native original-equation backward residual, not a forward error bound.
    pub relative_second_derivative_residual: f64,
}

impl PreparedAssemblyPortHarmonicResponse<'_, '_, '_, '_> {
    /// Reuse this loaded factor for the full second frequency derivative.
    ///
    /// Source derivatives remain in original component force/moment coordinates.
    /// None explicitly declares that source derivative zero. Both complete
    /// impedance derivatives are required on the structural exp(-i*w*t) axis;
    /// mutual entries and every source/load product-rule term are retained.
    /// Physical source work, reactions and port work pass the same original
    /// support/tie audit as the first-frequency response. This does not refit a
    /// load, infer heat from boundary work, or replace the response by a line shape.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_frequency_curvature(
        &self,
        component_loads: &[&[C64]],
        component_load_derivatives: Option<&[&[C64]]>,
        component_load_second_derivatives: Option<&[&[C64]]>,
        impedance_derivative: &[C64],
        impedance_second_derivative: &[C64],
        controls: HarmonicFrequencyCurvatureControls,
        budget: AssemblyPortCurvatureBudget,
        cx: &Cx<'_>,
    ) -> Result<AssemblyPortHarmonicCurvatureResponse, String> {
        checkpoint(cx)?;
        let equilibrium = self.ports.model.vibration.equilibrium;
        let components = &equilibrium.displacements;
        let d = components.iter().map(|v| v.len() as u128).sum::<u128>();
        let n = self.ports.model.vibration.mode_count() as u128;
        let f = equilibrium.assembly.model.free as u128;
        let p = self.ports.port_count() as u128;
        let t = equilibrium
            .assembly
            .slave_couplings()
            .iter()
            .map(|slave| slave.master_terms.len() as u128)
            .sum::<u128>();
        let bytes = 256 * d + 128 * (n + f + p) + 512 * components.len() as u128 + 1024;
        let terms = 24 * (2 * d + t + n * f + n) + 12 * p * d + 256 * d + 64 * p;
        if bytes > budget.maximum_payload_bytes as u128
            || terms > budget.maximum_observation_terms as u128
        {
            return Err(
                "assembled harmonic curvature exceeds physical payload or work limits".into(),
            );
        }
        let force = self.response.project_loads(component_loads, cx)?;
        let first_force = component_load_derivatives
            .map(|loads| self.response.project_loads(loads, cx))
            .transpose()?;
        let second_force = component_load_second_derivatives
            .map(|loads| self.response.project_loads(loads, cx))
            .transpose()?;
        let native = self
            .response
            .response
            .solve_with_frequency_curvature(
                &force,
                first_force.as_deref(),
                second_force.as_deref(),
                Some(impedance_derivative),
                Some(impedance_second_derivative),
                controls,
                cx,
            )
            .map_err(|error| error.to_string())?;
        let displacement_second_derivatives = self
            .response
            .expand(&native.displacement_second_derivative, cx)?;
        let velocity_second_derivatives = self
            .response
            .expand(&native.velocity_second_derivative, cx)?;
        let assembly = self.response.recover_frequency(
            native.frequency_response,
            component_loads,
            component_load_derivatives,
            cx,
        )?;
        let frequency_response =
            self.recover_frequency(assembly, component_loads, component_load_derivatives, cx)?;
        let port = native
            .port_load
            .ok_or("assembled harmonic curvature omitted its original load")?;
        let physical_port_flow_second_derivatives = self.physical_flows(
            &velocity_second_derivatives,
            &port.flow_second_derivative,
            cx,
        )?;
        let reaction_force_second_derivatives =
            self.physical_reactions(&port.effort_second_derivative, cx, |_, _, _, _| {})?;

        let assembly = &frequency_response.assembly;
        let zero = vec![C64::ZERO; components.iter().map(Vec::len).max().unwrap_or(0)];
        let (mut input, mut input_correction, mut input_scale) = (0.0, 0.0, 0.0);
        let (mut work, mut work_correction, mut work_scale) = (0.0, 0.0, 0.0);
        let load = assembly
            .modal
            .response
            .port_load
            .as_ref()
            .ok_or("assembled harmonic curvature omitted its original effort")?;
        let tangent = assembly
            .modal
            .port_load
            .as_ref()
            .ok_or("assembled harmonic curvature omitted its effort slope")?;
        for c in 0..components.len() {
            checkpoint(cx)?;
            let v = &assembly.velocities[c];
            let dv = &assembly.velocity_derivatives[c];
            let ddv = &velocity_second_derivatives[c];
            let df = component_load_derivatives.map_or(&zero[..v.len()], |loads| loads[c]);
            let ddf = component_load_second_derivatives.map_or(&zero[..v.len()], |loads| loads[c]);
            let source = harmonic_mean_power_curvature_with_cancel(
                v,
                dv,
                ddv,
                component_loads[c],
                df,
                ddf,
                || cx.checkpoint().is_err(),
            )
            .map_err(|error| error.to_string())?;
            add_real(&mut input, &mut input_correction, source);
            let boundary = harmonic_mean_power_curvature_with_cancel(
                v,
                dv,
                ddv,
                &frequency_response.reaction_forces[c],
                &frequency_response.reaction_force_derivatives[c],
                &reaction_force_second_derivatives[c],
                || cx.checkpoint().is_err(),
            )
            .map_err(|error| error.to_string())?;
            add_real(&mut work, &mut work_correction, -boundary);
            for j in 0..v.len() {
                if j % 1024 == 0 {
                    checkpoint(cx)?;
                }
                input_scale += 0.5
                    * (component_loads[c][j].abs() * ddv[j].abs()
                        + 2.0 * df[j].abs() * dv[j].abs()
                        + ddf[j].abs() * v[j].abs());
                // Retain the uncancelled port envelope in the audit scale.
                // Opposing physical attachments may have a tiny net reaction.
                for (index, attachment) in self.ports.ports.iter().enumerate() {
                    let weight = attachment.component_force_shapes[c][j].abs();
                    work_scale += 0.5
                        * weight
                        * (load.effort[index].abs() * ddv[j].abs()
                            + 2.0 * tangent.effort_derivative[index].abs() * dv[j].abs()
                            + port.effort_second_derivative[index].abs() * v[j].abs());
                }
            }
        }
        self.response.check_power(
            input,
            native.average_input_power_second_derivative,
            input_scale,
        )?;
        let viscous = native.average_viscous_power_second_derivative;
        let balance = input - viscous - work;
        self.check_work(
            work,
            port.average_absorbed_power_second_derivative,
            work_scale,
            balance,
            input_scale + viscous.abs() + work_scale,
        )?;
        checkpoint(cx)?;
        Ok(AssemblyPortHarmonicCurvatureResponse {
            frequency_response,
            displacement_second_derivatives,
            velocity_second_derivatives,
            modal_displacement_second_derivative: native.displacement_second_derivative,
            modal_velocity_second_derivative: native.velocity_second_derivative,
            modal_port_curvature: port,
            physical_port_flow_second_derivatives,
            reaction_force_second_derivatives,
            average_component_power_second_derivative: input,
            average_viscous_power_second_derivative: viscous,
            average_physical_load_power_second_derivative: work,
            physical_power_balance_residual_second_derivative: balance,
            relative_second_derivative_residual: native.relative_second_derivative_residual,
        })
    }
}

impl AssemblyPortHarmonicCurvatureResponse {
    /// Actual independently audited source work for shared power-peak queries.
    pub fn input_power_curvature_sample(&self) -> HarmonicPowerCurvatureSample {
        HarmonicPowerCurvatureSample {
            sample: HarmonicPowerSample {
                average_power_w: self.frequency_response.assembly.average_component_power_w,
                angular_derivative: self
                    .frequency_response
                    .assembly
                    .average_component_power_derivative,
            },
            angular_second_derivative: self.average_component_power_second_derivative,
        }
    }

    /// Actual independently audited signed load work, with no heat interpretation.
    pub fn port_power_curvature_sample(&self) -> HarmonicPowerCurvatureSample {
        HarmonicPowerCurvatureSample {
            sample: HarmonicPowerSample {
                average_power_w: self.frequency_response.average_physical_load_power_w,
                angular_derivative: self
                    .frequency_response
                    .average_physical_load_power_derivative,
            },
            angular_second_derivative: self.average_physical_load_power_second_derivative,
        }
    }

    /// Original symmetric glue loss, excluding external boundary work.
    pub fn viscous_power_curvature_sample(&self) -> HarmonicPowerCurvatureSample {
        HarmonicPowerCurvatureSample {
            sample: HarmonicPowerSample {
                average_power_w: self
                    .frequency_response
                    .assembly
                    .modal
                    .response
                    .average_viscous_power_w,
                angular_derivative: self
                    .frequency_response
                    .assembly
                    .modal
                    .average_viscous_power_derivative,
            },
            angular_second_derivative: self.average_viscous_power_second_derivative,
        }
    }
}
