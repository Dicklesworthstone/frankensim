//! Moving compliant stops without changing the string's material basis or span.
use super::*;
use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::geometry::ObstacleGeometryConfig;

/// Physical and numerical admission for prescribed longitudinal stop motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StopStationMotionLimits {
    /// Maximum mean station speed along the full source span [m/s].
    pub maximum_speed_m_s: f64,
    /// Maximum longitudinal distance over one mechanical sample [m]. The
    /// original normal root also requires positive endpoint/mean-map compliance.
    pub maximum_step_m: f64,
    /// Complete per-sample map/response work, as in `ObstacleGeometryConfig`.
    pub max_update_terms: usize,
}
impl StopStationMotionLimits {
    /// Positive finite physical limits and a nonzero numerical budget.
    pub fn validate(self) -> Result<(), RenderError> {
        if !self.maximum_speed_m_s.is_finite() || self.maximum_speed_m_s <= 0.0
            || !self.maximum_step_m.is_finite() || self.maximum_step_m <= 0.0
            || self.max_update_terms == 0 {
            return Err(error("stop migration needs positive finite speed/distance and bounded update work"));
        }
        Ok(())
    }
}

pub(in super::super) struct StationState {
    pub(in super::super) accepted: Vec<f64>,
    limits: StopStationMotionLimits,
}

impl MultiBowStringBank {
    /// Admit changing stop stations along the original supported string spans.
    ///
    /// Each step projects the source centreline map at both endpoints. The
    /// original normal contact supplies -mean_N R and reports the additional
    /// work R delta_N^T mean_q. Both transverse planes, torsion, other strings,
    /// body/cavity modes, causal gas history and geometric storage are retained.
    /// A receiver footprint and normal direction stay as authored. This is a
    /// prescribed moving compliant stop, not a solved finger or a pitch edit.
    ///
    /// Stops supplying bow normal loads move only with their matching tangent
    /// through the bow-station entry point. Independent stops can migrate in
    /// the same bank. Endpoint motion may refuse on the existing
    /// contact, energy, pressure, nonlinear-domain or compliance limits.
    pub fn with_stop_station_motion(mut self, limits: StopStationMotionLimits,
        gate: &CancelGate) -> Result<Self, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        limits.validate()?;
        if self.source.stops.is_empty() || self.stop_stations.is_some() || self.system.samples_rendered() != 0 {
            return Err(error("stop station motion requires one binding to existing stops before playback"));
        }
        self.system = self.system.with_obstacle_geometry_motion(ObstacleGeometryConfig {
            max_update_terms: limits.max_update_terms,
        }, gate).map_err(RenderError::Friction)?;
        self.stop_stations = Some(StationState {
            accepted: self.source.stops.iter().map(|s| s.location.station_fraction).collect(), limits,
        });
        Ok(self)
    }

    /// Accepted stations in original stop order; source descriptions are unchanged.
    #[must_use]
    pub fn stop_station_fractions(&self) -> Option<&[f64]> {
        self.stop_stations.as_ref().map(|s| s.accepted.as_slice())
    }

    /// Original physical migration and numerical work limits.
    #[must_use]
    pub fn stop_station_motion_limits(&self) -> Option<StopStationMotionLimits> {
        self.stop_stations.as_ref().map(|s| s.limits)
    }

    /// Advance all bows, stops and modal coordinates with prescribed endpoint
    /// stations and optional final normal offsets, both in original stop order.
    /// A missing normal-offset slice holds the accepted normal translations.
    /// Positions publish only with the complete physical step. Ordinary `step`
    /// and `step_with_stop_offsets` hold the last accepted stations.
    pub fn step_with_stop_stations(&mut self, external: &[f64], drives: &[FrictionDrive],
        stations: &[f64], offsets_m: Option<&[f64]>) -> Result<&MultiFrictionFrame, RenderError>
    {
        self.step_with_stop_stations_under_gate(external, drives, stations, offsets_m, &CancelGate::new_clock_free())
    }

    /// Same physical transaction with cancellation polling.
    pub fn step_with_stop_stations_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive],
        stations: &[f64], offsets_m: Option<&[f64]>, gate: &CancelGate) -> Result<&MultiFrictionFrame, RenderError>
    {
        let endpoints = self.prepare_stop_station_geometry(stations, gate)?;
        if self.bow_stations.is_some() {
            self.step_with_bow_geometry(external, drives, None, Some(&endpoints), offsets_m, gate)
        } else {
            self.system.step_with_obstacle_geometry_under_gate(external, drives, &endpoints, offsets_m, gate)
        }.map_err(RenderError::Friction)?;
        // No fallible operation follows publication; accepted source locations
        // and the generic owner's accepted maps always describe the same sample.
        self.stop_stations.as_mut().expect("admitted stop station state").accepted.copy_from_slice(stations);
        Ok(self.system.last_frame().expect("complete stop step was accepted"))
    }

    pub(in super::super) fn prepare_stop_station_geometry(&self, stations: &[f64],
        gate: &CancelGate) -> Result<Vec<ObstacleGeometry>, RenderError>
    {
        self.prepare_stop_station_geometry_with_bows(stations, None, gate)
    }

    pub(in super::super) fn prepare_stop_station_geometry_with_bows(&self, stations: &[f64],
        bow_stations: Option<&[f64]>, gate: &CancelGate) -> Result<Vec<ObstacleGeometry>, RenderError>
    {
        if gate.is_requested() { return Err(RenderError::Coupled(ModalCouplingError::Cancelled)); }
        let state = self.stop_stations.as_ref().ok_or_else(|| error("stop station motion has not been admitted"))?;
        if stations.len() != self.source.stops.len() {
            return Err(error("stop migration requires one station per original contact"));
        }
        for (i, ((&end, &start), stop)) in stations.iter().zip(&state.accepted).zip(&self.source.stops).enumerate() {
            let distance = (end-start)*self.cards[stop.location.string].length_m;
            let speed = distance/self.system.sample_period_s();
            if !end.is_finite() || !(0.0 < end && end < 1.0) || !distance.is_finite()
                || distance.abs() > state.limits.maximum_step_m || !speed.is_finite()
                || speed.abs() > state.limits.maximum_speed_m_s {
                return Err(error("stop migration exceeds its interior span, distance or speed limit"));
            }
            // Spatial identity matters even when a truncated modal map happens
            // to be equal at two distinct stations (e.g. a symmetric sine mode).
            if let Some(binding) = self.source.bow_contact_loads.iter().find(|b| b.stop == i) {
                match bow_stations {
                    Some(bows) if bows.get(binding.bow).is_none_or(|bow| bow.to_bits() != end.to_bits()) =>
                        return Err(error("moving normal and tangential contacts must occupy the same string station")),
                    None if end.to_bits() != start.to_bits() =>
                        return Err(error("a contact-loaded bow normal cannot migrate without its tangential contact")),
                    _ => {},
                }
            }
        }
        self.stop_geometry(&self.source.stops, Some(stations), state.limits.max_update_terms, gate)
    }
}
