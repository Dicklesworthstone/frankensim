//! Prescribed pump and inlet histories on the original fixed-mass flow graph.
use super::*;
use std::sync::Arc;

/// One continuous piecewise-linear pump/inlet sample.
///
/// The multiplier applies to every interior and boundary mass rate, so fixed
/// cell masses and an admitted stationary composition remain consistent.
/// Each knot must name the same subset of inlet boundary IDs. Unlisted inlets
/// retain their original enthalpy, in their own original material's datum.
#[derive(Clone, Debug, PartialEq)]
pub struct AdvectionScheduleKnot {
    pub time_s: f64,
    pub flow_multiplier: f64,
    pub inlet_enthalpies_j_kg: BTreeMap<u64, f64>,
}

/// Immutable, bounded, continuous prescribed flow and inlet-enthalpy history.
///
/// Values interpolate between strictly ordered knots and remain constant
/// outside their time range. Nonnegative common scaling preserves directions,
/// local mass continuity and frozen stationary species. It does not solve
/// momentum, changing composition, density or cell mass. Discontinuous jumps
/// require a separate event law and are not approximated by duplicate times.
#[derive(Clone, Debug, PartialEq)]
pub struct AdvectionSchedule {
    knots: Vec<AdvectionScheduleKnot>,
}

struct ScheduleSample<'a> {
    from: &'a AdvectionScheduleKnot,
    to: &'a AdvectionScheduleKnot,
    fraction: f64,
    multiplier: f64,
}

/// Original prescribed data, shared across cloned explicit candidates. Only
/// the current native flow buffers are resampled; physical observers and their
/// histories are never cloned or rebuilt by a source-epoch preparation.
#[derive(Clone)]
pub(super) struct ExplicitAdvectionSchedule {
    pub(super) schedule: Arc<AdvectionSchedule>,
    flows: Arc<[MassFlow]>,
    boundaries: Arc<[MassBoundary]>,
    admitted: bool,
}

impl ExplicitAdvectionSchedule {
    pub(super) fn new<O>(
        flow: &Advected<O>,
        schedule: Arc<AdvectionSchedule>,
    ) -> Result<Self, DepositionError> {
        schedule.inlets(&flow.boundaries, &mut || Ok(()))?;
        Ok(Self {
            schedule,
            flows: flow.flows.clone().into(),
            boundaries: flow.boundaries.clone().into(),
            admitted: false,
        })
    }

    pub(super) fn reset(&mut self) {
        self.admitted = false;
    }

    pub(super) fn prepare<O>(
        &mut self,
        flow: &mut Advected<O>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        if !self.admitted {
            // Validate original inlet laws even during a stopped interval.
            // begin() renews this admission against each actual input history.
            self.schedule
                .admit_sources(flow, &self.boundaries, thermal, checkpoint)?;
            self.admitted = true;
        }
        let sample = self.schedule.sample(thermal.time_s(), checkpoint)?;
        sample.apply(&self.flows, &self.boundaries, flow, checkpoint)
    }
}

impl AdvectionSchedule {
    pub fn new(
        knots: Vec<AdvectionScheduleKnot>,
        maximum_knots: usize,
        maximum_inlets: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if knots.is_empty() {
            return Err(invalid("advection schedule requires at least one knot"));
        }
        if knots.len() > maximum_knots || knots[0].inlet_enthalpies_j_kg.len() > maximum_inlets {
            return Err(DepositionError::Budget("advection schedule samples"));
        }
        for (i, knot) in knots.iter().enumerate() {
            checkpoint()?;
            if !knot.time_s.is_finite()
                || knot.time_s < 0.0
                || !knot.flow_multiplier.is_finite()
                || knot.flow_multiplier < 0.0
                || (i != 0 && knot.time_s <= knots[i - 1].time_s)
                || knot.inlet_enthalpies_j_kg.len() != knots[0].inlet_enthalpies_j_kg.len()
            {
                return Err(invalid(
                    "advection schedule needs ordered times, nonnegative flow and identical inlet IDs",
                ));
            }
            for (j, ((id, h), expected_id)) in knot
                .inlet_enthalpies_j_kg
                .iter()
                .zip(knots[0].inlet_enthalpies_j_kg.keys())
                .enumerate()
            {
                if j % 256 == 0 {
                    checkpoint()?;
                }
                if id != expected_id || !h.is_finite() {
                    return Err(invalid(
                        "scheduled inlets need identical IDs and finite enthalpies",
                    ));
                }
            }
        }
        checkpoint()?;
        Ok(Self { knots })
    }

    pub fn knots(&self) -> &[AdvectionScheduleKnot] {
        &self.knots
    }

    /// The existing time owner clips to this strictly future slope change.
    pub fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.knots
            .get(self.knots.partition_point(|knot| knot.time_s <= after_s))
            .map(|knot| knot.time_s)
    }

    /// Freeze the prescribed rates and inlet h at one physical time.
    ///
    /// Original schedule, inlet-law, stationary-composition and carrier data
    /// are admitted against `thermal`, then the SAME sampler used by native
    /// implicit stages supplies the snapshot. This advances no physical state
    /// or audit. The receiving native solve still owns continuity, actual cell
    /// calorics/liquid support, transport and acceptance. Explicit upwind
    /// entropy histories are not implicit endpoint snapshots.
    pub fn flow_at_time(
        &self,
        flow: &Advected<()>,
        thermal: &DepositionSimulation<'_>,
        time_s: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Advected<()>, DepositionError> {
        self.admit(flow, thermal, &mut checkpoint)?;
        self.flow_at(flow, time_s, &mut checkpoint)
            .map(|(flow, _)| flow)
    }

    pub(in crate::enthalpy) fn admit(
        &self,
        flow: &Advected<()>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        checkpoint()?;
        flow.admit_implicit_transport()?;
        self.admit_sources(flow, &flow.boundaries, thermal, checkpoint)
    }

    fn inlets<'a>(
        &self,
        boundaries: &'a [MassBoundary],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, &'a MassBoundary>, DepositionError> {
        let mut inlets = BTreeMap::new();
        for boundary in boundaries {
            checkpoint()?;
            if let MassBoundary::Inlet { id, .. } = *boundary {
                inlets.insert(id, boundary);
            }
        }
        for id in self.knots[0].inlet_enthalpies_j_kg.keys() {
            checkpoint()?;
            if !inlets.contains_key(id) {
                return Err(invalid(
                    "scheduled enthalpy must name an original inlet boundary",
                ));
            }
        }
        Ok(inlets)
    }

    fn admit_sources<O>(
        &self,
        flow: &Advected<O>,
        boundaries: &[MassBoundary],
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        for boundary in boundaries {
            checkpoint()?;
            if !thermal.cells().contains_key(&boundary.identity().1) {
                return Err(invalid("inactive scheduled mass-flow boundary cell"));
            }
        }
        let inlets = self.inlets(boundaries, checkpoint)?;
        // Keep the original composition/common-reference owner, including
        // during pump-off intervals. Never solve an underdetermined zero-flow
        // species graph or reinterpret inlet h in a receiving mixture's law.
        if let Some(species) = &flow.stationary_species {
            species.validate_laws(
                thermal,
                boundaries,
                &flow.selective_outlets,
                &flow.selective_outlet_laws,
            )?;
        }
        for knot in &self.knots {
            checkpoint()?;
            for (&id, boundary) in &inlets {
                checkpoint()?;
                let MassBoundary::Inlet {
                    cell,
                    material_id,
                    specific_h_j_kg,
                    ..
                } = **boundary
                else {
                    return Err(DepositionError::Arithmetic);
                };
                if material_id != thermal.cells()[&cell].material_id
                    && flow.stationary_species.is_none()
                {
                    return Err(invalid(
                        "scheduled inlet material/datum does not match receiving cell",
                    ));
                }
                let h = knot
                    .inlet_enthalpies_j_kg
                    .get(&id)
                    .copied()
                    .unwrap_or(specific_h_j_kg);
                let law = thermal
                    .materials
                    .get(&material_id)
                    .ok_or(invalid("unregistered scheduled inlet caloric law"))?;
                if law.storage_basis() != ThermalStorageBasis::Enthalpy {
                    return Err(invalid("scheduled inlet requires an enthalpy caloric law"));
                }
                let t = law.temperature_k(h)?;
                if !t.is_finite() || t <= 0.0 {
                    return Err(DepositionError::MaterialDomain);
                }
                if let Some(carrier) = flow.liquid_carrier {
                    let mut inlet = **boundary;
                    if let MassBoundary::Inlet {
                        specific_h_j_kg, ..
                    } = &mut inlet
                    {
                        *specific_h_j_kg = h;
                    }
                    carrier.admit_boundary(&inlet)?;
                }
            }
        }
        checkpoint()?;
        Ok(())
    }

    /// One scratch snapshot per implicit stage, reused for every Newton trial.
    /// The original flow's pending epoch and accepted audit are never touched.
    pub(in crate::enthalpy) fn flow_at(
        &self,
        flow: &Advected<()>,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Advected<()>, f64), DepositionError> {
        let sample = self.sample(time_s, checkpoint)?;
        let mut result = flow.clone();
        sample.apply(&flow.flows, &flow.boundaries, &mut result, checkpoint)?;
        Ok((result, sample.multiplier))
    }

    fn sample(
        &self,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ScheduleSample<'_>, DepositionError> {
        checkpoint()?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(invalid(
                "scheduled advection needs a finite nonnegative stage time",
            ));
        }
        let upper = self.knots.partition_point(|knot| knot.time_s <= time_s);
        let (a, b, fraction) = if upper == 0 {
            (&self.knots[0], &self.knots[0], 0.0)
        } else if upper == self.knots.len() {
            (&self.knots[upper - 1], &self.knots[upper - 1], 0.0)
        } else {
            let a = &self.knots[upper - 1];
            let b = &self.knots[upper];
            (a, b, finite((time_s - a.time_s) / (b.time_s - a.time_s))?)
        };
        let multiplier = interpolate(a.flow_multiplier, b.flow_multiplier, fraction)?;
        Ok(ScheduleSample {
            from: a,
            to: b,
            fraction,
            multiplier,
        })
    }
}

impl ScheduleSample<'_> {
    fn apply<O>(
        &self,
        original_flows: &[MassFlow],
        original_boundaries: &[MassBoundary],
        result: &mut Advected<O>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), DepositionError> {
        if result.flows.len() != original_flows.len() {
            return Err(invalid("scheduled flow topology changed"));
        }
        for (edge, original) in result.flows.iter_mut().zip(original_flows) {
            checkpoint()?;
            edge.kg_s = scale_rate(original.kg_s, self.multiplier)?;
        }
        // Retain capacity while stopped. A later sample rebuilds these same
        // original streams, rather than scaling the last sample cumulatively.
        result.boundaries.clear();
        for &original in original_boundaries {
            checkpoint()?;
            let mut boundary = original;
            match &mut boundary {
                MassBoundary::Inlet {
                    id,
                    kg_s,
                    specific_h_j_kg,
                    ..
                } => {
                    *kg_s = scale_rate(*kg_s, self.multiplier)?;
                    if let Some(&ha) = self.from.inlet_enthalpies_j_kg.get(&*id) {
                        *specific_h_j_kg =
                            interpolate(ha, self.to.inlet_enthalpies_j_kg[&*id], self.fraction)?;
                    }
                }
                MassBoundary::Outlet { kg_s, .. } => {
                    *kg_s = scale_rate(*kg_s, self.multiplier)?;
                }
            }
            if self.multiplier != 0.0 {
                result.boundaries.push(boundary);
            }
        }
        // A stopped sample carries no external mass. Original composition,
        // component references and carrier identity remain bound; nonzero
        // samples recheck native mobility and per-cell mass balance.
        checkpoint()?;
        Ok(())
    }
}

fn scale_rate(rate: f64, multiplier: f64) -> Result<f64, DepositionError> {
    let value = finite(rate * multiplier)?;
    if rate != 0.0 && multiplier != 0.0 && value == 0.0 {
        return Err(DepositionError::Arithmetic);
    }
    Ok(value)
}

fn interpolate(a: f64, b: f64, fraction: f64) -> Result<f64, DepositionError> {
    if a == b || fraction == 0.0 {
        return Ok(a);
    }
    if fraction == 1.0 {
        return Ok(b);
    }
    // Difference coordinates preserve small changes on a large common datum.
    // Opposite signs use the convex form to avoid overflowing b-a.
    finite(if a.is_sign_negative() == b.is_sign_negative() {
        a + fraction * (b - a)
    } else {
        (1.0 - fraction) * a + fraction * b
    })
}
