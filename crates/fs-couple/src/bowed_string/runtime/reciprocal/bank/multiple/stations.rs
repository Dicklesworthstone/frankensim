//! Prescribed bow-position changes on the original complete string/body state.
//!
//! A contact's tangent is the source surface map at the interval midpoint.
//! It transforms material velocity and force, not a displacement coordinate:
//! slip is u*dt - B^T*delta_q, with no fictitious Bdot*q contribution. The
//! original simultaneous friction solver accounts for work and dissipation.
//! Changing a station never changes length, tension, eigenmodes or history.

use super::*;
use crate::render::RenderError;
use crate::render::schedule::force::coupled::contact::friction::multiple::{
    obstacles::geometry::ObstacleGeometry,
    surfaces::geometry::{FrictionEndpointGeometry, FrictionGeometry, FrictionGeometryConfig},
};

/// Physical and numerical admission of prescribed longitudinal bow motion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BowStationMotionLimits {
    /// Maximum mean speed along the original supported string span [m/s].
    pub maximum_speed_m_s: f64,
    /// Maximum travel along that span during one mechanical sample [m].
    pub maximum_step_m: f64,
    /// Per-sample bound on complete contact response reconstruction. It includes
    /// original material and radiation memory, as in `FrictionGeometryConfig`.
    pub max_update_terms: usize,
}
impl BowStationMotionLimits {
    /// Admit positive finite physical limits and a nonzero work budget.
    pub fn validate(self) -> Result<(), ModalFrictionError> {
        if !self.maximum_speed_m_s.is_finite()
            || self.maximum_speed_m_s <= 0.0
            || !self.maximum_step_m.is_finite()
            || self.maximum_step_m <= 0.0
            || self.max_update_terms == 0
        {
            return Err(error(
                "bow migration needs positive finite speed/distance and bounded update work",
            ));
        }
        Ok(())
    }
}

pub(super) struct BowStationState {
    pub(super) accepted: Vec<f64>,
    pub(super) limits: BowStationMotionLimits,
    pub(super) fiber_support: Option<ribbon::motion::FiberSupportState>,
}

fn error(what: &'static str) -> ModalFrictionError {
    ModalCouplingError::Invalid { what }.into()
}
fn source_error(failure: ReciprocalStringError) -> ModalFrictionError {
    match failure {
        ReciprocalStringError::Network(error) => error.into(),
        ReciprocalStringError::Friction(error) => error,
        // A station projection evaluates no new material law. Original source
        // coefficients are retained; an invalid projection is a geometry refusal.
        ReciprocalStringError::Material(_) => {
            error("original material-string surface map could not be projected")
        }
    }
}

impl MultiBowStringBank {
    /// Admit bow-position changes without replacing the physical realization.
    ///
    /// Locations are fractions of the SAME supported span, including its
    /// afterlength. Both transverse polarizations and any surface-torque map
    /// come from the original material source and loaded modal mass. Receiver
    /// footprints and fiber contact stations stay fixed in their own bases.
    /// A finite ribbon translates as a whole: width and cell order are retained.
    /// This models motion along the string, not travel along bow hair, support
    /// acceleration, a change of speaking length, or a frequency assignment.
    /// Separate `with_fiber_support_motion` can add hair travel and inertia for
    /// dry friction and source-projected normal fiber contacts.
    ///
    /// Bind once before playback, after any normal-load bindings. Each bound
    /// normal follows its bow's SAME endpoint automatically. Its existing
    /// contact law supplies the solved load and accounts for migration work;
    /// there is no independently prescribed normal-force waveform. Other stops
    /// hold or follow separately declared paths. Material shear layers still
    /// refuse because this does not transport their finite patch histories.
    /// This is a midpoint held-map discretization; refine its step size and the
    /// original physical basis. It does not promise allocation-free playback.
    pub fn with_bow_station_motion(
        mut self,
        limits: BowStationMotionLimits,
        gate: &CancelGate,
    ) -> Result<Self, RenderError> {
        if gate.is_requested() {
            return Err(RenderError::Coupled(ModalCouplingError::Cancelled));
        }
        limits.validate().map_err(RenderError::Friction)?;
        if self.bow_stations.is_some() || self.system.samples_rendered() != 0 {
            return Err(RenderError::Friction(error(
                "bow migration binds once before playback",
            )));
        }
        self.admit_bow_station_sources()
            .map_err(RenderError::Friction)?;
        self.system
            .validate_friction_geometry_config(FrictionGeometryConfig {
                max_update_terms: limits.max_update_terms,
            })
            .map_err(RenderError::Friction)?;
        if !self.source.bow_contact_loads.is_empty() && self.stop_stations.is_none() {
            self = self.with_stop_station_motion(stops::stations::StopStationMotionLimits {
                maximum_speed_m_s: limits.maximum_speed_m_s,
                maximum_step_m: limits.maximum_step_m,
                max_update_terms: limits.max_update_terms,
            }, gate)?;
        }
        self.bow_stations = Some(BowStationState {
            accepted: self
                .source
                .bows
                .iter()
                .map(|b| b.location.station_fraction)
                .collect(),
            limits,
            fiber_support: None,
        });
        Ok(self)
    }

    pub(super) fn admit_bow_station_sources(&self) -> Result<(), ModalFrictionError> {
        for port in 0..self.source.bows.len() {
            if self.system.material_shear_layer(port).is_some() {
                return Err(error(
                    "bow migration requires material shear-layer history transport",
                ));
            }
        }
        Ok(())
    }

    /// Last accepted physical endpoint of every bow, in original contact order.
    /// The original source and its initial stations remain available separately.
    #[must_use]
    pub fn bow_station_fractions(&self) -> Option<&[f64]> {
        Some(&self.bow_stations.as_ref()?.accepted)
    }

    /// Original physical and per-sample numerical migration limits.
    #[must_use]
    pub fn bow_station_motion_limits(&self) -> Option<BowStationMotionLimits> {
        Some(self.bow_stations.as_ref()?.limits)
    }

    /// Advance with final bow stations and optional original normal translations.
    /// One value per actual bow is required; ribbon cells must translate together.
    /// The accepted station, material/contact state, pressure and clock publish
    /// together. Ordinary `step` holds the final station, settling the midpoint
    /// map from the previous interval onto that actual endpoint.
    pub fn step_with_bow_stations(
        &mut self,
        external: &[f64],
        drives: &[FrictionDrive],
        stations: &[f64],
        stop_offsets_m: Option<&[f64]>,
    ) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        self.step_with_bow_stations_under_gate(
            external,
            drives,
            stations,
            stop_offsets_m,
            &CancelGate::new_clock_free(),
        )
    }

    /// Same complete transaction with cancellation polling before publication.
    pub fn step_with_bow_stations_under_gate(
        &mut self,
        external: &[f64],
        drives: &[FrictionDrive],
        stations: &[f64],
        stop_offsets_m: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        self.step_with_bow_geometry(external, drives, Some(stations), None, stop_offsets_m, gate)
    }

    /// Advance prescribed bow and stop locations on the SAME clock.
    /// Stop locations are in original stop order and use their existing source
    /// projection and actuator-work account. Explicit endpoints of contact-load
    /// normals must equal their associated bow endpoints. If stop locations are
    /// omitted, those normals follow their bows automatically and independent
    /// stops hold. Optional normal translations are absolute offsets.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_bow_and_stop_stations_under_gate(
        &mut self,
        external: &[f64],
        drives: &[FrictionDrive],
        bow_stations: &[f64],
        stop_stations: Option<&[f64]>,
        stop_offsets_m: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<&MultiFrictionFrame, RenderError> {
        let normal_geometry = stop_stations
            .map(|stations| self.prepare_stop_station_geometry_with_bows(stations, Some(bow_stations), gate))
            .transpose()?;
        self.step_with_bow_geometry(
            external,
            drives,
            Some(bow_stations),
            normal_geometry.as_deref(),
            stop_offsets_m,
            gate,
        )
        .map_err(RenderError::Friction)?;
        if let Some(stations) = stop_stations {
            self.stop_stations
                .as_mut()
                .expect("admitted stop motion")
                .accepted
                .copy_from_slice(stations);
        }
        Ok(self
            .system
            .last_frame()
            .expect("complete bow/stop step was accepted"))
    }

    pub(super) fn step_with_bow_geometry(
        &mut self,
        external: &[f64],
        drives: &[FrictionDrive],
        stations: Option<&[f64]>,
        normal_geometry: Option<&[ObstacleGeometry]>,
        stop_offsets_m: Option<&[f64]>,
        gate: &CancelGate,
    ) -> Result<&MultiFrictionFrame, ModalFrictionError> {
        if gate.is_requested() {
            return Err(ModalCouplingError::Cancelled.into());
        }
        self.admit_bow_station_sources()?;
        let state = self
            .bow_stations
            .as_ref()
            .ok_or_else(|| error("bow station motion has not been admitted"))?;
        let end = stations.unwrap_or(&state.accepted);
        if end.len() != self.source.bows.len() || drives.len() != self.source.bows.len() {
            return Err(error(
                "bow migration needs one station and drive per original contact",
            ));
        }
        let limits = state.limits;
        let support = state.fiber_support.as_ref()
            .map(|motion| motion.prepare(&self.source, &self.system, external, drives, gate))
            .transpose()?;
        for ((&end, &start), bow) in end.iter().zip(&state.accepted).zip(&self.source.bows) {
            let distance = (end - start) * self.cards[bow.location.string].length_m;
            let speed = distance / self.system.sample_period_s();
            if !interior(end)
                || !distance.is_finite()
                || distance.abs() > limits.maximum_step_m
                || !speed.is_finite()
                || speed.abs() > limits.maximum_speed_m_s
            {
                return Err(error(
                    "bow migration exceeds its interior span, distance or speed limit",
                ));
            }
        }
        if let Some(ribbon) = self.source.ribbon() {
            validate_ribbon_stations(ribbon, self.cards[ribbon.centre().string].length_m, end)?;
        }
        // A mapped normal is the SAME spatial patch as its tangent. A bow-only
        // command therefore moves that normal, retaining all independent stops
        // and the original normal law, receiver, material state and clearance.
        // Explicit stop endpoints have already been checked by the joint entry
        // point. Keep automatically generated locations staged until acceptance.
        let normal_stations = if normal_geometry.is_none() && (!self.source.bow_contact_loads.is_empty()
            || support.as_ref().is_some_and(|s| s.has_normal_contacts())) {
            let mut stations = self.stop_station_fractions()
                .ok_or_else(|| error("moving bow patch lost its admitted normal geometry"))?.to_vec();
            for binding in &self.source.bow_contact_loads {
                stations[binding.stop] = end[binding.bow];
            }
            Some(stations)
        } else { None };
        let automatic_normal_geometry = normal_stations.as_ref()
            .map(|stations| self.prepare_stop_station_geometry_with_bows(stations, Some(end), gate))
            .transpose().map_err(|failure| match failure {
                RenderError::Coupled(error) => error.into(),
                RenderError::Friction(error) => error,
                _ => error("joint bow/normal station projection exceeds its source or motion limits"),
            })?;
        let normal_geometry = normal_geometry.or(automatic_normal_geometry.as_deref());
        let fiber_normal_geometry = support.as_ref()
            .map(|s| s.normal_geometry(normal_geometry)).transpose()?.flatten();
        let normal_geometry = fiber_normal_geometry.as_deref().or(normal_geometry);
        // An unchanged interval uses the same source arithmetic as original
        // binding. After a move it also replaces the previous midpoint map by
        // the accepted endpoint, without altering vibration or material memory.
        let mut geometry = Vec::with_capacity(end.len());
        let mut endpoints = (0..end.len()).any(|port| self.system.interface(port).is_some())
            .then(|| Vec::with_capacity(end.len()));
        for (i, ((&end, &start), bow)) in end
            .iter()
            .zip(&state.accepted)
            .zip(&self.source.bows)
            .enumerate()
        {
            if gate.is_requested() {
                return Err(ModalCouplingError::Cancelled.into());
            }
            let station = if end.to_bits() == start.to_bits() {
                start
            } else {
                start + 0.5 * (end - start)
            };
            let string = bow.location.string;
            let shapes = self.source.bank.strings[string]
                .source
                .surface_shape(&self.cards[string], station)
                .map_err(source_error)?;
            if let Some(endpoints) = &mut endpoints {
                let left = |station| self.source.bank.strings[string].source
                    .surface_shape(&self.cards[string], station)
                    .map(|shapes| ModalAttachment { component: string, shapes })
                    .map_err(source_error);
                let initial = left(start)?;
                let final_attachment = left(end)?;
                endpoints.push(match &support {
                    Some(support) => support.endpoint_geometry(i, initial, final_attachment),
                    None => {
                        let right = self.system.port(i)
                            .ok_or_else(|| error("original bow port is absent"))?.right.clone();
                        FrictionEndpointGeometry {
                            start: FrictionGeometry { left: initial, right: right.clone() },
                            end: FrictionGeometry { left: final_attachment, right },
                            start_speed_m_s: drives[i].speed_m_s, end_speed_m_s: drives[i].speed_m_s,
                        }
                    },
                });
            }
            geometry.push(FrictionGeometry {
                left: ModalAttachment {
                    component: string,
                    shapes,
                },
                right: match &support {
                    Some(support) => Some(support.receiver(i)),
                    None => self.system.port(i)
                        .ok_or_else(|| error("original bow port is absent"))?.right.clone(),
                },
            });
        }
        let (external, drives) = support.as_ref()
            .map_or((external, drives), |s| (s.forces.as_slice(), s.drives.as_slice()));
        let config = FrictionGeometryConfig { max_update_terms: limits.max_update_terms };
        if let Some(endpoints) = &endpoints {
            self.system.step_with_friction_geometry_and_endpoints(external, drives, &geometry,
                endpoints, normal_geometry, stop_offsets_m, config, gate)?;
        } else {
            self.system.step_with_friction_geometry(external, drives, &geometry,
                normal_geometry, stop_offsets_m, config, gate)?;
        }
        let frame = self.system.last_frame().expect("accepted friction geometry step");
        if let Some(support) = support {
            self.bow_stations.as_mut().expect("admitted bow motion")
                .fiber_support.as_mut().expect("admitted fiber support").accept(support, &self.system, frame);
        }
        if let Some(end) = stations {
            self.bow_stations
                .as_mut()
                .expect("admitted bow motion")
                .accepted
                .copy_from_slice(end);
        }
        if let Some(end) = normal_stations {
            self.stop_stations.as_mut().expect("admitted normal patch motion")
                .accepted.copy_from_slice(&end);
        }
        Ok(frame)
    }
}

// A ribbon is one finite geometric interval, not independently steered point
// bows. Allow only floating-point reconstruction roundoff in its source spacing.
pub(super) fn validate_ribbon_stations(
    ribbon: &ribbon::UniformBowRibbon,
    length_m: f64,
    stations: &[f64],
) -> Result<(), ModalFrictionError> {
    if stations.len() != ribbon.stations().len() || stations.is_empty() {
        return Err(error(
            "moving ribbon needs one endpoint per original spatial cell",
        ));
    }
    let expected = ribbon::UniformBowRibbon::try_new(
        StringBowLocation {
            string: ribbon.centre().string,
            station_fraction: stations[0],
        },
        length_m,
        ribbon.width_m(),
        stations.len(),
    )
    .map_err(source_error)?;
    for (&actual, &expected) in stations.iter().zip(expected.stations()) {
        if !actual.is_finite()
            || (actual - expected).abs() > 8.0 * f64::EPSILON * actual.abs().max(expected.abs())
        {
            return Err(error(
                "all ribbon cells must translate together at the original width and order",
            ));
        }
    }
    Ok(())
}
