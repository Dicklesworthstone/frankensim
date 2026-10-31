//! Original player-side complete acoustic fields on the opposite NSF face.
//! The ideal series source obeys p_zone - p_field = rise and supplies
//! -rise * Q_out * dt. Neither acoustic assembly carries mean gas advection.

use super::*;
use crate::nsf_dg::NsfDgRadiationSolve;

pub(super) struct UpstreamField<'a, 'c> {
    pub(super) load: ThermoviscousDgLoad<'a, 'c>,
    reference: GasState,
}

/// Accepted upstream interface and complete acoustic storage. The original
/// field, molecular/wall states and every finite duct remain in
/// `upstream_load().system().last_frame()` without duplicating those fields.
#[derive(Clone, Debug)]
pub struct NsfThermoviscousDgUpstreamFrame {
    pub face: usize,
    pub incoming_characteristic_m_s: f64,
    /// Selected original upstream patch/face, already included in `faces`.
    pub field_exchange: NsfBoundaryExchange,
    /// Flow OUT of the nonlinear zone and into the upstream FE port [m³/s].
    pub inlet_flow_m3_s: f64,
    pub inlet_pressure_pa: f64,
    pub inlet_work_j: f64,
    /// Work of the explicit ideal series rise, -rise * Q_out * dt [J].
    /// A transverse finite supply has zero series rise and its own inventory.
    pub source_work_j: f64,
    /// This upstream FE plus every one of its original attached finite ducts.
    pub acoustic_stored_energy_j: f64,
    pub acoustic_balance_residual_j: f64,
    pub acoustic_balance_tolerance_j: f64,
    pub pressure_residual_pa: f64,
    pub projection_covariance_j: f64,
    /// Original resolved work plus series work minus upstream FE inlet work.
    pub acoustic_work_residual_j: f64,
    pub port_work_residual_j: f64,
    pub mass_projection_difference_kg: f64,
}

pub(super) struct UpstreamTrial<'s> {
    pub(super) trial: ThermoviscousDgLoadTrial<'s>,
    pub(super) frame: NsfThermoviscousDgUpstreamFrame,
}

impl<'a, 'c, L: Lattice3<27>> NsfThermoviscousDg<'a, 'c, L> {
    /// Bind a second complete field at the opposite original face or selected
    /// source patch. Its authored initial fields, material/thermal memory and
    /// all finite ducts remain intact. Both acoustic fields and the nonlinear
    /// zone must occupy distinct physical volumes supplied by the caller.
    ///
    /// One explicit nonlinear policy serves both faces and any downstream
    /// mouths. A later attachment must agree with that policy. Upstream mouth
    /// coordinates require their own source binding and are refused here.
    /// A transverse finite supply may coexist with this passive acoustic load.
    pub fn with_upstream_load(mut self, load: ThermoviscousDgLoad<'a, 'c>,
        solve: NsfDgRadiationSolve, gate: &CancelGate) -> Result<Self, NsfDgError>
    {
        poll(gate)?;
        self.admit_joint_solve(solve)?;
        let high = !self.high;
        let reference = self.zone.characteristic_reference_gas(high).map_err(NsfDgError::Flow)?;
        let traces = if self.source_patch().is_some() { self.zone.open_patch_traces_on(high) }
            else { self.zone.open_boundary_traces(high) }
            .ok_or(NsfDgError::Input("upstream complete field requires its original selected quadrature"))?;
        let area = traces.first().ok_or(NsfDgError::Input("upstream field interface is empty"))?
            .area_m2 * traces.len() as f64;
        let other_area = load.area_m2();
        if self.steps != 0 || self.zone.steps() != 0
            || self.load.system().accepted_steps() != 0 || self.upstream.is_some()
            || load.system().accepted_steps() != 0 || load.mouth_count() != 0
            || load.system().source().gas() != &reference
            || !area.is_finite() || area <= 0.0 || !other_area.is_finite() || other_area <= 0.0
            || (area - other_area).abs() > 32.0 * f64::EPSILON * area.max(other_area)
            || !same_clock(self.load.duration_s(), load.duration_s())
            || !same_clock(self.zone.law().dt_s * self.limits.substeps as f64, load.duration_s())
        {
            return Err(NsfDgError::Input(
                "upstream complete field requires one cold original owner, matching face gas/area/clock and no unbound mouths"));
        }
        poll(gate)?;
        self.upstream = Some(UpstreamField { load, reference });
        self.bind_joint_solve(solve);
        Ok(self)
    }

    pub fn upstream_load(&self) -> Option<&ThermoviscousDgLoad<'a, 'c>> {
        self.upstream.as_ref().map(|owner| &owner.load)
    }

    /// Ideal series pressure rise between the upstream field and nonlinear
    /// source. A finite material supply retains its exclusive own drive form.
    pub fn step_with_upstream_pressure_rise(&mut self, pressure_rise_pa: f64,
        gate: &CancelGate) -> Result<NsfThermoviscousDgFrame, NsfDgError>
    {
        self.step_with_drive(NsfDgDrive::UpstreamPressureRise(pressure_rise_pa), gate)
    }
}

impl UpstreamField<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn preview_step(&mut self, face: usize, exchange: NsfBoundaryExchange,
        incoming: f64, rise: f64, relation: DgPortResponse, limits: NsfDgLimits,
        gate: &CancelGate) -> Result<UpstreamTrial<'_>, NsfDgError>
    {
        let dt = self.load.duration_s();
        let trial = self.load.preview_step(exchange.normal_volume_flow_m3_s, gate)
            .map_err(NsfDgError::Thermoviscous)?;
        let pressure = trial.mean_pressure_pa();
        let work = trial.inlet_work_j();
        let audit = acoustic_interface::audit(exchange, pressure, work, rise, relation,
            self.reference, dt, limits)?;
        let frame = NsfThermoviscousDgUpstreamFrame {
            face, incoming_characteristic_m_s: incoming, field_exchange: exchange,
            inlet_flow_m3_s: exchange.normal_volume_flow_m3_s,
            inlet_pressure_pa: pressure, inlet_work_j: work, source_work_j: audit.source_work,
            acoustic_stored_energy_j: trial.frame().stored_energy_j,
            acoustic_balance_residual_j: trial.frame().balance_residual_j,
            acoustic_balance_tolerance_j: trial.frame().balance_tolerance_j,
            pressure_residual_pa: audit.pressure, projection_covariance_j: audit.covariance,
            acoustic_work_residual_j: audit.work, port_work_residual_j: audit.port,
            mass_projection_difference_kg: audit.mass,
        };
        Ok(UpstreamTrial { trial, frame })
    }
}
