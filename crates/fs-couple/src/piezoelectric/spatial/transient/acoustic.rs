//! Real incident/outgoing waves on the original piezoelectric midpoint solve.
//!
//! A retained uniform guide mode has U=b^T v and p=Z*U+2*p_in, Z=rho*c/A.
//! Its traction -b*p contributes 2*Z*b*b^T/dt to the SAME half-increment
//! mechanical operator. Source work, support reactions, charge, intrinsic
//! memories and optional physical shunts all share that coupled solve.
//!
//! These are real midpoint values: wave energies use dt*p^2/Z, without the
//! harmonic peak-phasor factor of one half. Outgoing acoustic energy leaves
//! the body; it is not ceramic heat. The caller owns the declared waveform
//! bandwidth and time refinement. A timestep alone cannot certify a spectrum,
//! and midpoint work closure does not certify temporal accuracy or a real guide.
//! Geometry, material coefficients and the retained uniform mode stay fixed.

use super::{
    BoundaryMap, ElectrodeCondition, Result, TetPiezoelectricControls, TetPiezoelectricError,
    TetPiezoelectricModel, TetPiezoelectricTransient, TetPiezoelectricTransientCandidate,
    invalid, load, poll, positive,
};
use fs_exec::Cx;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use std::collections::BTreeSet;

pub use super::super::harmonic::acoustic::{PlaneWaveAcousticPort, PlaneWavePortControls};

/// Original fluid owner and real incoming pressure on a prepared mesh trace.
#[derive(Clone, Copy)]
pub struct TransientAcousticPort<'a> {
    /// Exact original pressure/flow map, shared with harmonic acoustics.
    pub port: &'a PlaneWaveAcousticPort<'a>,
    /// Original immutable fluid/material/EOS owner, evaluated once per window.
    pub medium: &'a dyn FluidAcoustics,
    /// Real incident traveling-wave pressure at the interval midpoint [Pa].
    pub midpoint_incident_pressure_pa: f64,
    /// Caller-admitted upper frequency of BOTH the incident signal and actual
    /// interface motion [rad/s]. Must lie inside this guide's retained-mode
    /// domain. This declaration is not inferred or certified from timestep size.
    pub maximum_angular_frequency_rad_s: f64,
}

/// Waves and independent physical work for one complete midpoint window.
#[derive(Clone, Copy, Debug)]
pub struct TransientAcousticPortResponse {
    /// Coherent original density, adiabatic sound speed and bulk modulus.
    pub medium: FluidAcousticState,
    /// Full guide section, including any rigid part [m^2].
    pub section_area_m2: f64,
    /// Real midpoint volume velocity out of the solid into the guide [m^3/s].
    pub midpoint_volume_velocity_m3_s: f64,
    /// Total signed acoustic pressure acting on the solid [Pa], positive in compression.
    pub midpoint_pressure_pa: f64,
    /// Supplied incoming traveling-wave pressure [Pa].
    pub midpoint_incident_pressure_pa: f64,
    /// Computed outgoing traveling-wave pressure [Pa].
    pub midpoint_outgoing_pressure_pa: f64,
    /// Nonnegative incoming wave energy under midpoint quadrature [J].
    pub incident_energy_j: f64,
    /// Nonnegative outgoing wave energy under midpoint quadrature [J].
    pub outgoing_energy_j: f64,
    /// Signed pressure work through actual surface displacement [J].
    /// Negative values represent sound supplying the body.
    pub net_outward_work_j: f64,
    /// Traction work minus outgoing-minus-incoming wave energy [J].
    pub wave_energy_residual_j: f64,
}

impl<'model, 'mesh> TetPiezoelectricTransient<'model, 'mesh> {
    /// Advance acoustic sensing/actuation with the existing electrode topology.
    /// Complete native inertia, reciprocal coupling and intrinsic memories
    /// remain in the same solve. Stationary supports retain their accepted state.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_acoustic_ports_candidate(
        &self,
        dt_s: f64,
        acoustic_ports: &[TransientAcousticPort<'_>],
        endpoint_electrode_conditions: &[ElectrodeCondition],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricTransientCandidate<'model, 'mesh>> {
        self.advance_with_acoustic_loads(
            dt_s, None, endpoint_electrode_conditions, &[], &[], &[], acoustic_ports,
            midpoint_nodal_forces_n, endpoint_nodal_free_charge_c, controls, cx,
        )
        .map(|(candidate, _, _)| candidate)
    }

    /// The same guide loading with explicit velocities for every original
    /// support DOF. Actual reactions supply inertia, field and outgoing-wave
    /// work; prescribed support motion is never a separately fitted source.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_acoustic_support_motion_candidate(
        &self,
        dt_s: f64,
        acoustic_ports: &[TransientAcousticPort<'_>],
        endpoint_support_velocities_m_s: &[(usize, f64)],
        endpoint_electrode_conditions: &[ElectrodeCondition],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<TetPiezoelectricTransientCandidate<'model, 'mesh>> {
        self.advance_with_acoustic_loads(
            dt_s, Some(endpoint_support_velocities_m_s), endpoint_electrode_conditions,
            &[], &[], &[], acoustic_ports, midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c, controls, cx,
        )
        .map(|(candidate, _, _)| candidate)
    }

    /// Sound-to-electric harvesting with original physical resistors and
    /// held-current/endpoint-voltage drives. The load owner's complete external
    /// source-work audit includes signed acoustic outflow exactly once.
    /// Use None for stationary supports or Some with every original support DOF.
    /// Ceramic memory heat, resistor heat and outgoing waves stay distinct.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_acoustic_loads_candidate(
        &self,
        dt_s: f64,
        acoustic_ports: &[TransientAcousticPort<'_>],
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        drives: &[load::ElectrodeDrive],
        loads: &[load::ResistiveElectrodeLoad],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<(TetPiezoelectricTransientCandidate<'model, 'mesh>, load::LoadedPiezoelectricStep)> {
        self.advance_loaded_with_acoustic_candidate(
            dt_s, endpoint_support_velocities_m_s, acoustic_ports, drives, loads,
            midpoint_nodal_forces_n, endpoint_nodal_free_charge_c, controls, cx,
        )
    }
}

pub(super) struct PreparedTransientAcousticPort<'a> {
    port: &'a PlaneWaveAcousticPort<'a>,
    medium: FluidAcousticState,
    impedance: f64,
    incident: f64,
}

pub(super) fn prepare<'a>(
    model: &TetPiezoelectricModel<'_>,
    inputs: &[TransientAcousticPort<'a>],
    cx: &Cx<'_>,
) -> Result<Vec<PreparedTransientAcousticPort<'a>>> {
    if inputs.len() > model.tetrahedra.len().saturating_mul(4) {
        return Err(invalid("acoustic port count exceeds original exterior-face capacity"));
    }
    let mut faces = BTreeSet::new();
    let mut prepared = Vec::with_capacity(inputs.len());
    for input in inputs {
        poll(cx)?;
        if !core::ptr::eq(model, input.port.model)
            || !input.midpoint_incident_pressure_pa.is_finite()
            || !positive(input.maximum_angular_frequency_rad_s)
        {
            return Err(invalid("original acoustic trace, finite pressure and declared bandwidth required"));
        }
        for &face in &input.port.faces {
            poll(cx)?;
            if !faces.insert(face) {
                return Err(invalid("acoustic ports cannot load one face twice"));
            }
        }
        let medium = input.medium.acoustic_state().map_err(TetPiezoelectricError::AcousticMaterial)?;
        let impedance = medium.characteristic_impedance_pa_s_m() / input.port.section_area_m2();
        if !positive(impedance)
            || input.maximum_angular_frequency_rad_s / medium.sound_speed_m_s()
                > input.port.controls.maximum_wavenumber_rad_m
            || input.midpoint_incident_pressure_pa.abs() / medium.adiabatic_bulk_modulus_pa()
                > input.port.controls.maximum_relative_amplitude
        {
            return Err(invalid("transient acoustic port leaves its declared mode or amplitude domain"));
        }
        prepared.push(PreparedTransientAcousticPort {
            port: input.port, medium, impedance, incident: input.midpoint_incident_pressure_pa,
        });
    }
    Ok(prepared)
}

pub(super) fn stamp(
    ports: &[PreparedTransientAcousticPort<'_>],
    dt: f64,
    boundary: &mut BoundaryMap,
    entries: &mut Vec<(usize, usize, f64)>,
    cx: &Cx<'_>,
) -> Result<()> {
    let coefficient = 2.0 / dt;
    for port in ports {
        let weights = port.port.weights_m2();
        for (i, &(row, a)) in weights.iter().enumerate() {
            poll(cx)?;
            if let Some(r) = boundary.index[row] {
                boundary.rhs[r] -= 2.0 * (a * port.incident);
            }
            // Compute each reciprocal entry once, including the exact columns
            // of prescribed support motion in native affine elimination.
            for &(col, b) in &weights[i..] {
                let value = coefficient * (port.impedance * a) * b;
                if !value.is_finite() {
                    return Err(invalid("transient acoustic reaction is unrepresentable"));
                }
                entries.push((row, col, value));
                if row != col {
                    entries.push((col, row, value));
                }
            }
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn recover(
    ports: &[PreparedTransientAcousticPort<'_>],
    dt: f64,
    increment: &[f64],
    before_displacement: &[[f64; 3]],
    after_displacement: &[[f64; 3]],
    effective_force: &mut [[f64; 3]],
    energy_tolerance_j: f64,
    cx: &Cx<'_>,
) -> Result<Vec<TransientAcousticPortResponse>> {
    let mut responses = Vec::with_capacity(ports.len());
    for port in ports {
        poll(cx)?;
        let mut half_volume = 0.0;
        let mut actual_volume = 0.0;
        for &(dof, b) in port.port.weights_m2() {
            poll(cx)?;
            half_volume += b * increment[dof];
            actual_volume += b * (after_displacement[dof / 3][dof % 3]
                - before_displacement[dof / 3][dof % 3]);
        }
        let flow = (2.0 / dt) * half_volume;
        let outgoing = port.incident + port.impedance * flow;
        let pressure = outgoing + port.incident;
        let incident_energy_j = dt * port.incident * (port.incident / port.impedance);
        let outgoing_energy_j = dt * outgoing * (outgoing / port.impedance);
        let net_outward_work_j = pressure * actual_volume;
        // Factored wave difference avoids subtracting two nearly equal
        // positive energies, independently of actual displacement work.
        let wave_net = dt * (outgoing - port.incident) * (pressure / port.impedance);
        let wave_energy_residual_j = net_outward_work_j - wave_net;
        let area = port.port.section_area_m2();
        let amplitude = port.port.controls.maximum_relative_amplitude;
        if [flow, pressure, outgoing, incident_energy_j, outgoing_energy_j,
            net_outward_work_j, wave_energy_residual_j].iter().any(|v| !v.is_finite())
            || pressure.abs() / port.medium.adiabatic_bulk_modulus_pa() > amplitude
            || outgoing.abs() / port.medium.adiabatic_bulk_modulus_pa() > amplitude
            || flow.abs() / area / port.medium.sound_speed_m_s() > amplitude
        {
            return Err(invalid("solved transient acoustic waves leave their small-signal domain"));
        }
        if wave_energy_residual_j.abs() > energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "actual acoustic traction work differs from outgoing minus incoming wave energy",
            ));
        }
        for &(dof, b) in port.port.weights_m2() {
            poll(cx)?;
            let force = &mut effective_force[dof / 3][dof % 3];
            *force -= b * pressure;
            if !force.is_finite() {
                return Err(invalid("complete acoustic traction is unrepresentable"));
            }
        }
        responses.push(TransientAcousticPortResponse {
            medium: port.medium,
            section_area_m2: area,
            midpoint_volume_velocity_m3_s: flow,
            midpoint_pressure_pa: pressure,
            midpoint_incident_pressure_pa: port.incident,
            midpoint_outgoing_pressure_pa: outgoing,
            incident_energy_j, outgoing_energy_j, net_outward_work_j, wave_energy_residual_j,
        });
    }
    Ok(responses)
}
