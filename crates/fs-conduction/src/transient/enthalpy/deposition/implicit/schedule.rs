//! Exact schedule boundaries over the same transactional thermal owner.
use super::*;

impl<G: Clone> DepositionSimulation<'_, G> {
    /// Advance through the prefix of `events` at or before `end_time_s`.
    ///
    /// All supplied event headers must be ordered by time and strictly increasing
    /// identity, and none may precede the accepted state. Equal-time events are
    /// allowed. Only the returned `report.events` prefix is consumed; pass the
    /// remaining suffix on the next call. No future cell/contact is allocated or
    /// included in heat transfer before its event. A birth exactly at the horizon
    /// is included, but receives no pre-birth cooling.
    ///
    /// Every substep is shortened at the next event. The entire call, including
    /// already traversed events, is atomic on cancellation, budget, domain, or
    /// balance refusal. A clone before the call is an in-memory restart point.
    /// Resuming at the same event/substep boundaries preserves operation order;
    /// splitting a timestep at a new observation time changes the discretization.
    pub fn run_until(
        &mut self,
        cx: &Cx<'_>,
        events: &[DepositionEvent<G>],
        end_time_s: f64,
        config: DepositionStepConfig,
    ) -> Result<DepositionReport, DepositionError> {
        validate_config(config)?;
        poll(cx, 0)?;
        if !end_time_s.is_finite() || end_time_s < self.time_s {
            return Err(DepositionError::Invalid("time must be finite and nondecreasing"));
        }
        if events.len() > self.budget.max_events {
            return Err(DepositionError::Budget("schedule headers"));
        }
        let mut previous_time = self.time_s;
        let mut previous_id = self.last_event;
        for (index, event) in events.iter().enumerate() {
            poll(cx, index)?;
            if !event.time_s.is_finite() || event.time_s < previous_time
                || previous_id.is_some_and(|id| event.id <= id) {
                return Err(DepositionError::Invalid("schedule must have nondecreasing times and increasing IDs"));
            }
            previous_time = event.time_s;
            previous_id = Some(event.id);
        }
        let mut trial = self.clone();
        let before = trial.stored_enthalpy_j();
        let mut report = DepositionReport::default();
        for event in events.iter().take_while(|e| e.time_s <= end_time_s) {
            trial.march(cx, event.time_s, config, &mut report)?;
            trial.activate(cx, event, &mut report)?;
        }
        trial.march(cx, end_time_s, config, &mut report)?;
        finish(&trial, before, &mut report, config.energy_tolerance_j)?;
        poll(cx, report.steps)?;
        *self = trial;
        Ok(report)
    }
}

use crate::transient::deposition as resolved;

impl DepositionSimulation<'_, resolved::CellGeometry> {
    /// Consume the existing explicit lane's resolved physical event vocabulary.
    /// This is an input adapter, not a slicer or a state conversion: both lanes
    /// retain their distinct time discretizations. Material IDs must name the
    /// same charts supplied at construction. Only the returned event prefix is
    /// consumed. `first_event_id` identifies its first element across retries.
    ///
    /// The caller explicitly supplies the large-enclosure radiation temperature;
    /// material emissivities are those in this implicit simulation. The complete
    /// original centreline geometry is retained in `geometry.detail`. A stable
    /// orthonormal frame is constructed from the nonzero centreline direction;
    /// its transverse axes are a convention, not measured anisotropy data.
    ///
    /// Finite positive conductances map to the existing declared K/W resistance
    /// type. Closed zero-conductance contacts and infinite bed reservoirs are
    /// deliberately refused, not silently dropped or replaced by free surface.
    /// Represent a finite bed/support with ordinary cells, or retain the explicit
    /// lane for reservoir/insulating-contact schedules. The richer native birth
    /// interface accepts caller-supplied material frames without this convention.
    #[allow(clippy::too_many_arguments)]
    pub fn run_resolved_until(
        &mut self,
        cx: &Cx<'_>,
        events: &[resolved::TimedEvent],
        first_event_id: EventId,
        end_time_s: f64,
        radiation_temperature_k: f64,
        config: DepositionStepConfig,
    ) -> Result<DepositionReport, DepositionError> {
        positive(radiation_temperature_k)?;
        validate_config(config)?;
        poll(cx, 0)?;
        if !end_time_s.is_finite() || end_time_s < self.time_s {
            return Err(DepositionError::Invalid("time must be finite and nondecreasing"));
        }
        if events.len() > self.budget.max_events {
            return Err(DepositionError::Budget("resolved schedule headers"));
        }
        let mut previous_time = self.time_s;
        let mut converted = Vec::new();
        for (index, event) in events.iter().enumerate() {
            poll(cx, index)?;
            if !event.time_s.is_finite() || event.time_s < previous_time {
                return Err(DepositionError::Invalid("resolved event times must be ordered"));
            }
            previous_time = event.time_s;
            if event.time_s > end_time_s { continue; }
            let offset = u64::try_from(index).map_err(|_| DepositionError::Budget("event identity"))?;
            let id = first_event_id.0.checked_add(offset)
                .ok_or(DepositionError::Budget("event identity"))?;
            let mut output = DepositionEvent { id: EventId(id), time_s: event.time_s,
                births: Vec::new(), contacts: Vec::new(), environments: Vec::new() };
            match &event.event {
                resolved::Event::Birth(cell) => {
                    let geometry = cell.geometry;
                    let mut position_m = [0.0; 3];
                    let mut direction = [0.0; 3];
                    for i in 0..3 {
                        position_m[i] = finite(0.5 * geometry.start_m[i] + 0.5 * geometry.end_m[i])?;
                        direction[i] = finite(geometry.end_m[i] - geometry.start_m[i])?;
                    }
                    output.births.push(MaterialBirth {
                        id: CellId(cell.id.0), body: BodyId(cell.body_id), segment: SegmentId(cell.segment_id),
                        material: MaterialId(cell.material_id.0), mass_kg: cell.mass_kg,
                        specific_enthalpy_j_kg: cell.specific_h_j_kg,
                        geometry: CellGeometry { position_m, material_axes: centreline_frame(direction)?,
                            volume_m3: geometry.volume_m3, surface_area_m2: geometry.surface_area_m2,
                            detail: geometry },
                        environment: resolved_environment(cell.boundary, radiation_temperature_k)?,
                    });
                }
                resolved::Event::Contact(link) => {
                    positive(link.conductance_w_k)?;
                    let resistance = crate::ThermalResistanceTerm::declared(
                        format!("resolved-interface-{}", link.id), 1.0 / link.conductance_w_k,
                        crate::ResistanceUncertainty::Unstated,
                        "caller-resolved conductance; inverse is K/W, not adhesion strength",
                    )?;
                    output.contacts.push(ContactChange::Set(ThermalContact {
                        id: InterfaceId(link.id), a: CellId(link.a.0), b: CellId(link.b.0),
                        area_m2: link.area_m2, resistance: SeriesThermalResistance::new(vec![resistance])?,
                    }));
                }
                resolved::Event::Separation(id) => output.contacts.push(ContactChange::Open(InterfaceId(*id))),
                resolved::Event::Boundary(id, boundary) => output.environments.push((CellId(id.0),
                    resolved_environment(*boundary, radiation_temperature_k)?)),
            }
            converted.push(output);
        }
        self.run_until(cx, &converted, end_time_s, config)
    }
}

fn resolved_environment(boundary: resolved::Boundary, radiation_k: f64) -> Result<Environment, DepositionError> {
    if boundary.reservoir.is_some() {
        return Err(DepositionError::Invalid("implicit adapter requires finite bed/support cells, not a reservoir"));
    }
    Ok(Environment { ambient_k: boundary.ambient_k, radiation_k,
        convection_w_m2_k: boundary.htc_w_m2_k, source_w: boundary.source_w })
}

fn centreline_frame(direction: [f64; 3]) -> Result<[[f64; 3]; 3], DepositionError> {
    let norm = direction[0].hypot(direction[1]).hypot(direction[2]);
    positive(norm)?;
    let x = direction.map(|v| v / norm);
    // Choose the least-aligned world axis to avoid an ill-conditioned cross product.
    let axis = if x[0].abs() <= x[1].abs() && x[0].abs() <= x[2].abs() { 0 }
        else if x[1].abs() <= x[2].abs() { 1 } else { 2 };
    let mut helper = [0.0; 3]; helper[axis] = 1.0;
    let y = cross(helper, x);
    let norm_y = y[0].hypot(y[1]).hypot(y[2]);
    positive(norm_y)?;
    let y = y.map(|v| v / norm_y);
    Ok([x, y, cross(x, y)])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
}
