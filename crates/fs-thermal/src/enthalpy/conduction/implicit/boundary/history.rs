//! Continuous prescribed bath temperatures on the native conductive boundary.
use super::{ImplicitBoundaryExchangeState, ImplicitConductionChain, ImplicitConductionError};
use crate::enthalpy::backward_euler::PrescribedBathHeat;
use crate::enthalpy::conduction::implicit::{Endpoint, conductance_derivative, endpoint_distance};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, finite, positive};
use std::collections::BTreeMap;

/// A finite sample of an externally maintained bath's absolute temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BathTemperatureKnot {
    pub time_s: f64,
    pub temperature_k: f64,
}

/// Bounded, continuous, piecewise-linear prescribed bath temperature.
///
/// The explicit closed time range never extrapolates. Positive Kelvin samples,
/// strict time ordering and a caller-supplied sample budget are required.
/// Author constant tails explicitly and use native events for sharp jumps.
/// The same immutable schedule may be borrowed by any number of original
/// reservoir or ambient faces without duplicating its samples.
#[derive(Clone, Debug, PartialEq)]
pub struct BathTemperatureSchedule {
    knots: Vec<BathTemperatureKnot>,
}

impl BathTemperatureSchedule {
    pub fn new(
        knots: Vec<BathTemperatureKnot>,
        maximum_knots: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, DepositionError> {
        checkpoint()?;
        if knots.len() < 2 {
            return Err(DepositionError::InvalidInput(
                "bath-temperature history requires at least two knots",
            ));
        }
        if knots.len() > maximum_knots {
            return Err(DepositionError::Budget("bath-temperature history samples"));
        }
        for (i, knot) in knots.iter().enumerate() {
            checkpoint()?;
            if !knot.time_s.is_finite()
                || knot.time_s < 0.0
                || !positive(knot.temperature_k)
                || (i != 0 && knot.time_s <= knots[i - 1].time_s)
            {
                return Err(DepositionError::InvalidInput(
                    "bath-temperature history requires ordered times and positive finite Kelvin",
                ));
            }
        }
        checkpoint()?;
        Ok(Self { knots })
    }

    pub fn knots(&self) -> &[BathTemperatureKnot] {
        &self.knots
    }

    pub fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.knots
            .get(self.knots.partition_point(|knot| knot.time_s <= after_s))
            .map(|knot| knot.time_s)
    }

    /// Deterministic logarithmic lookup at the physical residual/stage time.
    pub fn temperature_at(&self, time_s: f64) -> Result<f64, DepositionError> {
        if !time_s.is_finite()
            || time_s < self.knots[0].time_s
            || time_s > self.knots[self.knots.len() - 1].time_s
        {
            return Err(DepositionError::InvalidInput(
                "bath-temperature query is outside the prescribed time domain",
            ));
        }
        let upper = self.knots.partition_point(|knot| knot.time_s <= time_s);
        let left = self.knots[upper - 1];
        if upper == self.knots.len() || time_s == left.time_s {
            return Ok(left.temperature_k);
        }
        let right = self.knots[upper];
        let fraction = (time_s - left.time_s) / (right.time_s - left.time_s);
        let temperature =
            finite((1.0 - fraction) * left.temperature_k + fraction * right.temperature_k)?;
        if !positive(temperature) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(temperature)
    }
}

#[derive(Clone, Copy)]
pub(in crate::enthalpy::conduction::implicit) struct BoundBathHistory<'a> {
    pub schedule: &'a BathTemperatureSchedule,
    authored_temperature_k: f64,
}

pub(in crate::enthalpy::conduction::implicit) fn temperature_at(
    bindings: &BTreeMap<CellId, BoundBathHistory<'_>>,
    id: CellId,
    authored_temperature_k: f64,
    time_s: f64,
) -> Result<f64, DepositionError> {
    match bindings.get(&id) {
        Some(history) => history.schedule.temperature_at(time_s),
        None => Ok(authored_temperature_k),
    }
}

fn single_bath(
    temperature_k: f64,
    cell_temperature_k: f64,
    conductance_w_k: f64,
    derivative_w_kg_j: f64,
) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
    ImplicitBoundaryExchangeState {
        baths: [
            PrescribedBathHeat {
                outward_power_w: finite(conductance_w_k * (cell_temperature_k - temperature_k))?,
                temperature_k,
            },
            PrescribedBathHeat {
                outward_power_w: 0.0,
                temperature_k,
            },
        ],
        power_enthalpy_derivative_w_kg_j: derivative_w_kg_j,
    }
    .validate()
}

impl<'a> ImplicitConductionChain<'a> {
    /// Prescribe continuous external bath temperatures on selected original faces.
    ///
    /// Native backward Euler, both SDIRK2 stages, observed runs and adaptive
    /// comparisons sample these immutable histories at their physical times.
    /// Every interval clips to the next knot. The existing boundary power and
    /// bath-entropy owner gives each face exactly one heat credit.
    ///
    /// Reservoirs retain their geometric conduction binding. Ambient histories
    /// support both geometric faces and native exposed-area Robin coefficients.
    /// Integrated, radiative and composed boundary laws receive the scheduled
    /// reference bath temperature; another bath declared by a composed law
    /// remains that law's responsibility. Original calorics, material-density
    /// geometry and all complete enthalpy tangents remain active.
    ///
    /// Authored boundaries are retained on accepted states. A native event
    /// changing a history-owned temperature refuses transactionally rather
    /// than silently replacing either owner. Other boundary event fields
    /// retain their existing admission rules.
    pub fn with_bath_temperature_histories(
        mut self,
        thermal: &DepositionSimulation<'_>,
        reservoirs: &[(CellId, &'a BathTemperatureSchedule)],
        ambient: &[(CellId, &'a BathTemperatureSchedule)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.reservoir_histories.is_empty()
            || !self.ambient_histories.is_empty()
            || (reservoirs.is_empty() && ambient.is_empty())
            || reservoirs.len() > self.conduction.reservoirs.len()
            || ambient.len() > self.ids.len()
        {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty bath-temperature history attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        for (bindings, selected, reservoir) in [
            (&mut self.reservoir_histories, reservoirs, true),
            (&mut self.ambient_histories, ambient, false),
        ] {
            for &(id, schedule) in selected {
                checkpoint()?;
                let cell = thermal
                    .cells
                    .get(&id)
                    .ok_or(ImplicitConductionError::Invalid(
                        "bath-temperature history requires an active original cell",
                    ))?;
                let authored_temperature_k = if reservoir {
                    cell.boundary
                        .reservoir
                        .ok_or(ImplicitConductionError::Invalid(
                            "bath-temperature history requires an original reservoir",
                        ))?
                        .temperature_k
                } else {
                    cell.boundary.ambient_k
                };
                schedule.temperature_at(thermal.time_s)?;
                if bindings
                    .insert(
                        id,
                        BoundBathHistory {
                            schedule,
                            authored_temperature_k,
                        },
                    )
                    .is_some()
                {
                    return Err(ImplicitConductionError::Invalid(
                        "bath-temperature history requires a unique original channel",
                    ));
                }
            }
        }
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }

    pub(in crate::enthalpy::conduction::implicit) fn admit_bath_histories(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        for (histories, reservoir) in [
            (&self.reservoir_histories, true),
            (&self.ambient_histories, false),
        ] {
            for (&id, history) in histories {
                checkpoint()?;
                let cell = thermal
                    .cells
                    .get(&id)
                    .ok_or(ImplicitConductionError::Invalid(
                        "bath-temperature history cell disappeared",
                    ))?;
                let actual = if reservoir {
                    cell.boundary
                        .reservoir
                        .ok_or(ImplicitConductionError::Invalid(
                            "bath-temperature history reservoir disappeared",
                        ))?
                        .temperature_k
                } else {
                    cell.boundary.ambient_k
                };
                if actual != history.authored_temperature_k {
                    return Err(ImplicitConductionError::Invalid(
                        "native boundary event conflicts with prescribed bath-temperature history",
                    ));
                }
                history.schedule.temperature_at(thermal.time_s)?;
            }
        }
        Ok(())
    }

    pub(in crate::enthalpy::conduction::implicit) fn apply_bath_histories(
        &self,
        thermal: &DepositionSimulation<'_>,
        endpoint: &mut Endpoint,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        // Only native unbound Robin channels need original exposed-area
        // reduction. Bound faces already retain that geometric admission.
        let exposed = if self
            .ambient_histories
            .keys()
            .any(|id| !self.conduction.ambient_faces.contains_key(id))
        {
            thermal.exposed_areas_m2(&self.ids, &mut *checkpoint)?
        } else {
            BTreeMap::new()
        };
        for (&id, history) in &self.reservoir_histories {
            checkpoint()?;
            if self.reservoir_exchanges.contains_key(&id) {
                // The original composed law was already sampled with this
                // actual bath temperature and owns both signed contributions.
                continue;
            }
            let index = self.index[&id];
            let state = endpoint.states[index];
            let bath = history.schedule.temperature_at(thermal.time_s)?;
            let binding = self.conduction.reservoirs[&id];
            let g = endpoint
                .boundary_values
                .reservoir_w_k
                .remove(&id)
                .ok_or(DepositionError::Arithmetic)?;
            let derivative = if let Some(boundary) = endpoint.reservoir_states.get(&id) {
                boundary.power_enthalpy_derivative_w_kg_j
            } else {
                let dg = conductance_derivative(
                    g,
                    binding.area_m2,
                    endpoint_distance(
                        &endpoint.layer_distances,
                        index,
                        binding.geometry.half_distance_m,
                    ),
                    state,
                )?;
                finite(g * state.temperature_slope_kg_k_j + dg * (state.temperature_k - bath))?
            };
            let exchange = single_bath(bath, state.temperature_k, g, derivative)?;
            endpoint
                .exchange_values
                .reservoir
                .insert(id, exchange.baths);
            endpoint.reservoir_exchange_states.insert(id, exchange);
        }
        for (&id, history) in &self.ambient_histories {
            checkpoint()?;
            if self.ambient_exchanges.contains_key(&id) {
                continue;
            }
            let index = self.index[&id];
            let state = endpoint.states[index];
            let bath = history.schedule.temperature_at(thermal.time_s)?;
            let (g, derivative) = if let Some(binding) = self.conduction.ambient_faces.get(&id) {
                let g = endpoint
                    .boundary_values
                    .ambient_w_k
                    .remove(&id)
                    .ok_or(DepositionError::Arithmetic)?;
                let derivative = if let Some(boundary) = endpoint.ambient_states.get(&id) {
                    boundary.power_enthalpy_derivative_w_kg_j
                } else {
                    let dg = conductance_derivative(
                        g,
                        binding.area_m2,
                        endpoint_distance(
                            &endpoint.layer_distances,
                            index,
                            binding.geometry.half_distance_m,
                        ),
                        state,
                    )?;
                    finite(g * state.temperature_slope_kg_k_j + dg * (state.temperature_k - bath))?
                };
                (g, derivative)
            } else {
                let g = finite(exposed[&id] * thermal.cells[&id].boundary.htc_w_m2_k)?;
                (g, finite(g * state.temperature_slope_kg_k_j)?)
            };
            let exchange = single_bath(bath, state.temperature_k, g, derivative)?;
            endpoint.exchange_values.ambient.insert(id, exchange.baths);
            endpoint.ambient_exchange_states.insert(id, exchange);
        }
        Ok(())
    }
}
