//! The original modal/Maxwell/air equations joined to a complete DG generator.
//! Retaining DG energy coordinates avoids division by a bare duct resonance.

pub use super::super::dg::ModalDgError as DgHarmonicCouplingError;
use super::*;
use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::dg::{DgDuct, DgHarmonicResponse, PortId};

/// Continuous AC response of both original owners, in exp(-i omega t).
#[derive(Clone, Debug, PartialEq)]
pub struct DgLoadedModalHarmonicResponse {
    /// Original modal response. Its input power includes the opposing duct
    /// reaction; `authored_external_power_w` is the external power of the joint
    /// system. Its residual describes the complete coupled equations.
    pub modal: ModalHarmonicResponse,
    /// Actual DG energy coordinates, physical fields, port pressures and losses.
    pub duct: DgHarmonicResponse,
    /// Physical flow INTO the original duct inlet [m3/s].
    pub inlet_flow_m3_s: C64,
    /// Physical inlet pressure opposing that flow [Pa].
    pub inlet_pressure_pa: C64,
    /// Mean power of the authored generalized forces only [W].
    pub authored_external_power_w: f64,
    /// Mean storage in all retained modal, material, radiation and duct states [J].
    pub mean_stored_energy_j: f64,
    /// Authored power minus modal input power minus DG input power [W].
    pub interface_power_residual_w: f64,
    /// Combined physical loss minus authored external power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of the complete original equations.
    pub solve_relative_residual: f64,
}
impl DgLoadedModalHarmonicResponse {
    /// Original mechanical/material/radiation loss plus DG wall loss and
    /// actual matched-boundary absorption [W].
    pub fn total_dissipation_w(&self) -> f64 {
        self.modal.total_dissipation_w()
            + self.duct.wall_dissipation_w
            + self.duct.matched_absorption_w
    }
}

fn duct_error(error: fs_duct::DuctError, gate: &CancelGate) -> DgHarmonicCouplingError {
    if gate.is_requested() {
        ModalCouplingError::Cancelled.into()
    } else {
        DgHarmonicCouplingError::Acoustic(AcousticRealizeError::Duct(error))
    }
}

impl CoupledModalSystem {
    /// Continuous forced response with one complete source-bound DG duct.
    ///
    /// `flow_column` maps ORIGINAL modal velocities to inlet volume flow:
    /// Q=b^T v. The opposite generalized force is -b p. The duct retains all
    /// its original energy coordinates, including wall memory and its actual
    /// closed/open/matched distal boundary. Exactly one declared DG flow port
    /// must exist and equal `inlet`; omitted external ports are never closed
    /// implicitly. Branches and rational terminations require their own retained
    /// states and are not supplied by a bare DgDuct.
    ///
    /// The existing complex LU solves the complete matrix, without a standalone
    /// impedance inversion, fitted poles or a coupling timestep. Constitutive
    /// band and spatial-resolution obligations stay with their original owners.
    /// No accepted history changes. The dense allowance is
    /// d^3+d^2*(connections+material_arms+5), with
    /// d=original_mechanical_modes+air_poles+DG_states <=1024.
    #[allow(clippy::too_many_arguments)]
    pub fn harmonic_response_with_duct(
        &self,
        omega_rad_s: f64,
        force: &[C64],
        duct: &DgDuct,
        inlet: PortId,
        flow_column: &[f64],
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<DgLoadedModalHarmonicResponse, DgHarmonicCouplingError> {
        poll(Some(gate))?;
        let view = duct.continuous_operator();
        view.admit_frequency(omega_rad_s)
            .map_err(|e| duct_error(e, gate))?;
        let n = self.mode_count();
        if flow_column.len() != n
            || flow_column.iter().any(|b| !b.is_finite())
            || !flow_column.iter().any(|b| *b != 0.0)
            || view.ports().len() != 1
            || view.ports()[0].0 != inlet
        {
            return Err(invalid("continuous modal/DG coupling requires one original inlet and a complete nonzero volume projection").into());
        }
        let acoustic_dimension = view.dimension();
        let mut operator =
            self.harmonic_operator(omega_rad_s, force, acoustic_dimension, maximum_work, gate)?;
        let dimension = operator.dimension;
        let offset = operator.original_dimension;
        let source = view.generator_per_s();
        let pressure_column = &view.ports()[0].1;
        let s = C64::new(0.0, -omega_rad_s);
        for i in 0..acoustic_dimension {
            poll(Some(gate))?;
            for j in 0..acoustic_dimension {
                operator.values[(offset + i) * dimension + offset + j] =
                    C64::from_re(-source[i * acoustic_dimension + j]);
            }
            let diagonal = (offset + i) * dimension + offset + i;
            operator.values[diagonal] = operator.values[diagonal] + s;
            for (j, &b) in flow_column.iter().enumerate() {
                // D_modal q + b B^T x = F; (s I-L)x - s B b^T q = 0.
                let coupling = pressure_column[i] * b;
                operator.values[j * dimension + offset + i] = C64::from_re(coupling);
                operator.values[(offset + i) * dimension + j] = -s.scale(coupling);
            }
        }
        let (state, relative) = solve_harmonic(
            &operator.values,
            dimension,
            force,
            self.config.solve_relative_tolerance,
            gate,
        )?;
        let flow = s * projection(flow_column, &state[..n]);
        let mut cancelled = || gate.is_requested();
        let duct = view
            .observe_harmonic(
                omega_rad_s,
                &state[offset..],
                &[(inlet, flow)],
                maximum_work,
                self.config.solve_relative_tolerance,
                &mut cancelled,
            )
            .map_err(|e| duct_error(e, gate))?;
        let pressure = duct
            .ports
            .iter()
            .find(|p| p.id == inlet)
            .ok_or_else(|| invalid("continuous DG observation omitted its original inlet"))?
            .pressure_pa;
        let applied: Vec<_> = force
            .iter()
            .zip(flow_column)
            .map(|(&f, &b)| f - pressure.scale(b))
            .collect();
        if applied.iter().any(|&f| !complex_finite(f)) {
            return Err(invalid("continuous DG reaction force overflows").into());
        }
        let modal =
            self.observe_harmonic_state(omega_rad_s, &applied, &state[..offset], relative, gate)?;
        let authored_external_power_w = force
            .iter()
            .zip(&modal.velocity_m_sqrt_kg_per_s)
            .map(|(&f, &v)| 0.5 * (f * v.conj()).re)
            .sum::<f64>();
        let mean_stored_energy_j = modal.mean_stored_energy_j + duct.mean_stored_energy_j;
        let interface_power_residual_w =
            authored_external_power_w - modal.external_power_w - duct.input_power_w;
        let mut result = DgLoadedModalHarmonicResponse {
            modal,
            duct,
            inlet_flow_m3_s: flow,
            inlet_pressure_pa: pressure,
            authored_external_power_w,
            mean_stored_energy_j,
            interface_power_residual_w,
            power_residual_w: 0.0,
            solve_relative_residual: relative,
        };
        result.power_residual_w = result.total_dissipation_w() - authored_external_power_w;
        if [
            result.authored_external_power_w,
            result.mean_stored_energy_j,
            result.interface_power_residual_w,
            result.power_residual_w,
            result.total_dissipation_w(),
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("continuous modal/DG physical power or storage overflows").into());
        }
        poll(Some(gate))?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
