//! The same nonlinear material valve between two original exact DG owners.
//!
//! Both duct inlets point away from the valve. For downstream flow Q and
//! independently supplied downstream body flow B, the tract flow is -(Q-B).
//! Each original owner supplies p = H + Z Q. Therefore the valve sees
//! Z = Z_up + Z_down, H = H_down, and an effective source
//! P_s + H_up + Z_up B. This is algebraic elimination of the two mean-port
//! laws, not another pressure solver or a delayed sample of tract feedback.
//! The physical source supplies only P_s (Q-B) dt. Each candidate retains its
//! original acoustic, branch, material-wall and terminal histories until both
//! port-work audits and the original mechanics have succeeded.

use super::{ApertureDg, ApertureDgLoad, DgApertureFrame, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dynamic::ApertureDrive;
use crate::bernoulli_aperture::paired::BladePropagator;
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::bernoulli_aperture::two_sided::{LoadedAperture, ValveTrial};
use crate::bernoulli_aperture::two_sided::supply::SupplyChannelFrame;
use crate::bernoulli_aperture::slit_channel::FixedChannel;
use crate::dg_port::{DgTerminalFrame, sum_terminal_work};
use fs_duct::dg::DgStepOutput;
use fs_material::gas::GasState;

/// Actual accepted upstream pressure, flow and original DG exchange.
#[derive(Debug, Clone, PartialEq)]
pub struct DgTractFrame {
    /// Aggregate original tract/branch propagation ledger and external ports.
    pub duct: DgStepOutput,
    /// Mean pressure at the tract inlet, excluding the ideal series source [Pa].
    pub inlet_pressure_pa: f64,
    /// Flow INTO the tract; the negative of aperture through-flow [m³/s].
    pub inlet_flow_m3_s: f64,
    /// Prescribed series pressure plus actual tract inlet pressure, less the
    /// explicitly attached fixed passage's pressure drop [Pa].
    pub upstream_face_pressure_pa: f64,
    /// Actual upstream-minus-downstream pressure across the valve [Pa].
    pub pressure_drop_pa: f64,
    /// Exact original DG inlet work [J], signed positive into the tract.
    pub inlet_work_j: f64,
    /// Every original tract terminal, in its admitted order.
    pub terminals: Vec<DgTerminalFrame>,
    /// Integrated tract inlet work minus mean pressure times volume [J].
    pub port_work_residual_j: f64,
}

/// Resolve the two cold inlet gases and SUM of original exact mean impedances.
/// The same small-signal volume flows cannot conserve mass across different
/// inlet densities. Both complete inlet gas states and clocks must match;
/// remote gas fields, geometry, walls and terminal loads remain independent.
/// Each source retains its original explicit DG, branch and terminal budgets.
pub fn series_inlet_state(
    upstream: &ApertureDgLoad,
    downstream: &ApertureDgLoad,
) -> Result<(GasState, f64), AcousticRealizeError> {
    let gas = *downstream.duct().inlet_gas();
    let impedance = upstream.impedance_pa_s_m3() + downstream.impedance_pa_s_m3();
    if upstream.has_driven_mouth() || downstream.has_driven_mouth()
        || upstream.state().step_count() != 0 || downstream.state().step_count() != 0
        || upstream.propagator().dt().to_bits() != downstream.propagator().dt().to_bits()
        || upstream.duct().inlet_gas() != &gas
        || !impedance.is_finite() || impedance <= 0.0
    {
        return Err(invalid(
            "DG tract requires two cold loads with one complete inlet gas and clock, finite summed impedance and no driven exterior mouth",
        ));
    }
    Ok((gas, impedance))
}

impl ApertureDg<LoadedAperture> {
    /// Bind unchanged original mechanics to a bore AND an upstream tract.
    /// Construct the valve with series_inlet_state's gas and impedance. All
    /// original signed modes, source loads, contact and material memory remain
    /// in that valve. Inlet and outlet observers continue to observe the bore.
    pub fn new_loaded_with_upstream(
        aperture: impl Into<LoadedAperture>,
        upstream: ApertureDgLoad,
        load: ApertureDgLoad,
    ) -> Result<Self, AcousticRealizeError> {
        let (gas, impedance) = series_inlet_state(&upstream, &load)?;
        let aperture = aperture.into();
        let (density, declared_impedance) = aperture.load();
        if aperture.accepted_steps() != 0
            || aperture.blade_propagator() != BladePropagator::ExactZoh
            || aperture.time_step_s().to_bits() != load.propagator().dt().to_bits()
            || density.to_bits() != gas.density.to_bits()
            || declared_impedance.to_bits() != impedance.to_bits()
            || matches!(&aperture, LoadedAperture::Paired(a) if a.gas().is_some_and(|g| g != &gas))
            || matches!(&aperture, LoadedAperture::Modal(a) if a.gas().is_some_and(|g| g != &gas))
        {
            return Err(invalid(
                "DG tract mechanics require cold ExactZoh state and the actual summed inlet impedance, complete gas and clock",
            ));
        }
        Ok(Self { aperture, load, upstream: Some(upstream), supply_channel: None })
    }

    pub(super) fn step_with_upstream(
        &mut self,
        drive: TubeDrive,
        forces_n: [f64; 2],
        pressures_pa: [f64; 2],
    ) -> Result<DgApertureFrame, AcousticRealizeError> {
        let upstream = self.upstream.as_mut().expect("selected upstream tract");
        let (mechanical, frame) = prepare_upstream_step(
            &self.aperture, upstream, &mut self.load, drive, forces_n, pressures_pa)?;
        self.aperture.accept(mechanical);
        Ok(frame)
    }
}

/// Original native series sample. Acoustic candidates publish only after every
/// pressure/work check; the mechanical trial stays private for the caller's
/// final publication. Finite backgrounds call this once on private load copies.
pub(super) fn prepare_upstream_step(
    source: &LoadedAperture,
    upstream: &mut ApertureDgLoad,
    load: &mut ApertureDgLoad,
    drive: TubeDrive,
    forces_n: [f64; 2],
    pressures_pa: [f64; 2],
) -> Result<(ValveTrial, DgApertureFrame), AcousticRealizeError> {
    let (mechanical, frame, _) = prepare_series_step(
        source, Some(upstream), load, None, drive, forces_n, pressures_pa,
    )?;
    Ok((mechanical, frame))
}

// Eliminate the original mean-port laws together. An optional fixed channel
// has dp_c = Z_c (Q-B) - H_c, so the valve sees Z_down + Z_up + Z_c and
// P_s + H_up + H_c + (Z_up + Z_c) B. No fluid or duct is stepped on lagged
// pressure. The source's private nonlinear trial remains its sole owner.
pub(super) fn prepare_series_step(
    source: &LoadedAperture,
    upstream: Option<&mut ApertureDgLoad>,
    load: &mut ApertureDgLoad,
    channel: Option<&FixedChannel>,
    drive: TubeDrive,
    forces_n: [f64; 2],
    pressures_pa: [f64; 2],
) -> Result<(ValveTrial, DgApertureFrame, Option<FixedChannel>), AcousticRealizeError> {
        if !drive.upstream_pressure_pa.is_finite() || !drive.body_flow_m3_s.is_finite() {
            return Err(invalid("DG tract source pressure and body flow must be finite"));
        }
        let up = upstream.as_deref().map(ApertureDgLoad::port_response).transpose()?;
        let down = load.port_response()?;
        let supply_load = channel.map(FixedChannel::pressure_load).transpose()?;
        let dt = load.propagator().dt();
        let effective_source = match up {
            Some(up) => drive.upstream_pressure_pa + up.history_pressure_pa
                + up.impedance_pa_s_m3 * drive.body_flow_m3_s,
            None => drive.upstream_pressure_pa,
        };
        let effective_source = match supply_load {
            Some(supply) => effective_source + supply.history_pressure_pa
                + supply.impedance_pa_s_m3 * drive.body_flow_m3_s,
            None => effective_source,
        };
        let mechanical = source.preview(ApertureDrive {
            upstream_pressure_pa: effective_source,
            incoming_pressure_pa: down.history_pressure_pa / 2.0,
            body_flow_m3_s: drive.body_flow_m3_s,
        }, forces_n, pressures_pa)?;
        let mut aperture = mechanical.frame();
        let spatial_pressure = mechanical.spatial_pressure();
        let q = aperture.bore_flow_m3_s;
        let through = q - drive.body_flow_m3_s;
        let expected_drop = effective_source - aperture.bore_pressure_pa;
        let up_id = upstream.as_deref().map(ApertureDgLoad::inlet);
        let down_id = load.inlet();

        // Both original candidates exist before either one can be accepted.
        // Dropping either trial leaves its complete material/acoustic state intact.
        let up_trial = upstream.map(|tract| tract.preview_step(-through)).transpose()?;
        let down_trial = load.preview_step(q)?;
        let supply_trial = channel.map(|channel| channel.preview_at_flow(through)).transpose()?;
        let supply = match (supply_load, &supply_trial) {
            (Some(relation), Some(channel)) => Some(SupplyChannelFrame {
                pressure_drop_pa: relation.pressure_drop_pa(through)?,
                channel: channel.last_step(),
            }),
            _ => None,
        };
        let up_output = up_trial.as_ref().map(|trial| trial.output());
        let down_output = down_trial.output();
        let pu = match (up_output, up_id) {
            (Some(output), Some(id)) => output.port_outputs.iter().find(|p| p.id == id)
                .ok_or_else(|| invalid("DG tract lost its original inlet"))?.conjugate_rate,
            _ => 0.0,
        };
        let pd = down_output.port_outputs.iter().find(|p| p.id == down_id)
            .ok_or_else(|| invalid("DG bore lost its original inlet"))?.conjugate_rate;
        let face_pressure = drive.upstream_pressure_pa + pu;
        let face_pressure = supply.map_or(face_pressure, |supply| face_pressure - supply.pressure_drop_pa);
        let pressure_drop = face_pressure - pd;
        let up_scale = up.map_or(0.0, |up| up.history_pressure_pa.abs()
            + (up.impedance_pa_s_m3 * through).abs() + pu.abs());
        let supply_scale = supply_load.map_or(0.0, |supply| supply.history_pressure_pa.abs()
            + (supply.impedance_pa_s_m3 * through).abs())
            + supply.map_or(0.0, |supply| supply.pressure_drop_pa.abs());
        let down_scale = down.history_pressure_pa.abs()
            + (down.impedance_pa_s_m3 * q).abs() + pd.abs();
        let pressure_scale = up_scale + down_scale + supply_scale + drive.upstream_pressure_pa.abs()
            + effective_source.abs() + aperture.bore_pressure_pa.abs();
        let up_work = up_trial.as_ref().map_or(0.0, |trial| trial.inlet_work_j());
        let down_work = down_trial.inlet_work_j();
        let up_defect = up_work - pu * (-through * dt);
        let down_defect = down_work - pd * (q * dt);
        let supply_defect = supply.map_or(0.0, |supply| supply.channel.pressure_work_j
            - supply.pressure_drop_pa * (through * dt));
        let work_scale = up_work.abs() + down_work.abs()
            + supply.map_or(0.0, |supply| supply.channel.pressure_work_j.abs())
            + ((up_scale + supply_scale) * through.abs() + down_scale * q.abs()) * dt;
        if ![pu, pd, face_pressure, pressure_drop, expected_drop, up_scale,
            down_scale, supply_scale, pressure_scale, up_defect, down_defect, supply_defect, work_scale]
            .iter().all(|x| x.is_finite())
            || up.is_some_and(|up| (pu - (up.history_pressure_pa - up.impedance_pa_s_m3 * through)).abs()
                > 1e-10 * up_scale.max(1e-300))
            || (pd - (down.history_pressure_pa + down.impedance_pa_s_m3 * q)).abs()
                > 1e-10 * down_scale.max(1e-300)
            || (pressure_drop - expected_drop).abs() > 1e-10 * pressure_scale.max(1e-300)
            || up_defect.abs() + down_defect.abs() + supply_defect.abs()
                > 1e-10 * work_scale.max(1e-300)
            || supply.is_some_and(|supply| (supply.channel.flow_m3_s - through).abs()
                > 1e-10 * (supply.channel.flow_m3_s.abs() + through.abs()).max(1e-300))
        {
            return Err(invalid("DG tract, supply passage, bore and valve disagree on same-step flow, pressure or work"));
        }

        // Map only observation coordinates back from the eliminated series
        // impedance. The original immutable trial still owns all mechanics.
        aperture.bore_pressure_pa = pd;
        aperture.outgoing_pressure_pa = down.history_pressure_pa / 2.0
            + down.impedance_pa_s_m3 * q;
        let stored_energy_j = aperture.stored_energy_j
            + up_output.map_or(0.0, |output| output.ledger.stored_energy_j)
            + down_output.ledger.stored_energy_j
            + supply.map_or(0.0, |supply| supply.channel.kinetic_energy_j);
        let storage_change_j = aperture.storage_change_j
            + up_output.map_or(0.0, |output| output.ledger.stored_energy_j)
            - up_trial.as_ref().map_or(0.0, |trial| trial.previous_stored_energy_j())
            + down_output.ledger.stored_energy_j - down_trial.previous_stored_energy_j()
            + supply.map_or(0.0, |supply| supply.channel.stored_delta_j);
        let inferred_total_dissipation_j = aperture.dissipated_energy_j
            + up_output.map_or(0.0, |output| output.ledger.cumulative_dissipated_j)
            - up_trial.as_ref().map_or(0.0, |trial| trial.previous_cumulative_dissipated_j())
            + down_output.ledger.cumulative_dissipated_j - down_trial.previous_cumulative_dissipated_j()
            + supply.map_or(0.0, |supply| supply.channel.viscous_dissipation_j
                + supply.channel.numerical_dissipation_j + supply.channel.exit_work_j);
        let upstream_work_j = drive.upstream_pressure_pa * (through * dt);
        let body_work_j = pd * (drive.body_flow_m3_s * dt);
        let terminal = down_trial.terminal_frame().copied();
        let terminals = down_trial.terminal_frames().to_vec();
        let up_terminals = up_trial.as_ref().map_or_else(Vec::new, |trial| trial.terminal_frames().to_vec());
        let terminal_work_j = sum_terminal_work(&terminals) + sum_terminal_work(&up_terminals);
        let balance = storage_change_j + inferred_total_dissipation_j + terminal_work_j
            - upstream_work_j - body_work_j - aperture.mechanical_work_j;
        if ![aperture.outgoing_pressure_pa, stored_energy_j, storage_change_j,
            inferred_total_dissipation_j, upstream_work_j, body_work_j,
            terminal_work_j, balance].iter().all(|x| x.is_finite())
        {
            return Err(invalid("DG tract candidate left the finite total energy/work domain"));
        }
        // No fallible operation follows the first publication.
        let up_duct = up_trial.map(|trial| trial.accept());
        let duct = down_trial.accept();
        Ok((mechanical, DgApertureFrame {
            aperture, spatial_pressure, duct,
            upstream: up_duct.map(|up_duct| DgTractFrame {
                duct: up_duct, inlet_pressure_pa: pu, inlet_flow_m3_s: -through,
                upstream_face_pressure_pa: face_pressure, pressure_drop_pa: pressure_drop,
                inlet_work_j: up_work, terminals: up_terminals, port_work_residual_j: up_defect,
            }),
            supply, terminal, terminals, stored_energy_j, storage_change_j,
            parameter_work_j: 0.0,
            inferred_total_dissipation_j, upstream_work_j, body_work_j,
            port_work_residual_j: up_defect + down_defect + supply_defect,
        }, supply_trial))
}

#[cfg(test)]
mod tests;
