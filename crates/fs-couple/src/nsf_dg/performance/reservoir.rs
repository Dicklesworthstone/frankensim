//! Physical radiation of the original transverse escape aperture.
//! Each fluid face tile retains its own signed continuum flow. The explicit
//! rigid placement puts that actual surface in the existing common baffle;
//! no bore-flow proxy, uniform mouth piston or mass/rho substitution is used.

use super::{ApertureObservation, NsfDgPerformance, RenderError, sizing};
use crate::acoustic_realize::AcousticRealizeError;
use crate::pcm_wav::baffled::gas::GasPropagationReport;
use crate::pcm_wav::baffled::tiles::{BaffledTilePressure, FlowSurfaceSample, TilePressureLimits};
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::NsfReservoirFlow;
use fs_math::det;
use std::collections::BTreeMap;

/// Explicit rigid placement of the original reservoir face in the shared
/// infinite baffle, facing its +z half space. Tangential coordinates are the
/// cyclic grid axes `(axis+1)%3, (axis+2)%3`; the second is reversed on a low
/// face so that its normal points OUT of the fluid. Cell edge zero maps to
/// `origin_m`, then `angle_rad` rotates both coordinates in the baffle plane.
/// This declares a coplanar assembly, not a free-field instrument body.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirObservation {
    pub origin_m: [f64; 2],
    pub angle_rad: f64,
}

pub(super) struct ReservoirPressure {
    pub aperture: usize,
    pub input: ReservoirObservation,
    pub surface: Vec<FlowSurfaceSample>,
    pub receiver: BaffledTilePressure,
    pub gas: Option<GasPropagationReport>,
    pub cells: Vec<[usize; 3]>,
    pub flows: Vec<f64>,
}

impl ReservoirPressure {
    pub fn accumulate(
        &mut self,
        accepted: &[NsfReservoirFlow],
        weight: f64,
    ) -> Result<(), RenderError> {
        let selected = accepted.iter().filter(|flow| flow.aperture == self.aperture);
        if selected.clone().count() != self.cells.len()
            || selected.clone()
                .zip(&self.cells)
                .any(|(flow, cell)| flow.cell != *cell)
        {
            return Err(sizing(
                "resolved reservoir observation lost its original tile basis",
            ));
        }
        for (sum, flow) in self.flows.iter_mut().zip(selected) {
            *sum += weight * flow.outward_volume_flow_m3_s;
        }
        Ok(())
    }
}

// All apertures share one original grid spacing. Equal-size rotated squares
// overlap only if their centres lie within sqrt(2)*dx. Spatial bins restrict
// this cold geometry check to neighbours instead of all pairs of flow tiles.
fn disjoint_surfaces(
    existing: &[ReservoirPressure], surface: &[FlowSurfaceSample],
    input: ReservoirObservation, dx: f64, cancelled: &mut impl FnMut() -> bool,
) -> Result<(), RenderError> {
    if existing.is_empty() { return Ok(()); }
    let span = core::f64::consts::SQRT_2 * dx;
    let key = |position: [f64; 3]| -> Result<(i64, i64), RenderError> {
        let x = (position[0] / span).floor();
        let y = (position[1] / span).floor();
        if !x.is_finite() || !y.is_finite() || x.abs().max(y.abs()) >= (1_u64 << 52) as f64 {
            return Err(sizing("reservoir placement exceeds its representable tile geometry"));
        }
        Ok((x as i64, y as i64))
    };
    let mut bins: BTreeMap<(i64, i64), Vec<([f64; 3], [f64; 2])>> = BTreeMap::new();
    for observer in existing {
        let orientation = [det::cos(observer.input.angle_rad), det::sin(observer.input.angle_rad)];
        for tile in &observer.surface {
            if cancelled() { return Err(sizing("reservoir geometry admission cancelled")); }
            bins.entry(key(tile.position_m)?).or_default().push((tile.position_m, orientation));
        }
    }
    let orientation = [det::cos(input.angle_rad), det::sin(input.angle_rad)];
    for tile in surface {
        if cancelled() { return Err(sizing("reservoir geometry admission cancelled")); }
        let (x, y) = key(tile.position_m)?;
        for i in x - 1..=x + 1 {
            for j in y - 1..=y + 1 {
                if let Some(candidates) = bins.get(&(i, j)) {
                    for &(position, other) in candidates {
                        let delta = [tile.position_m[0] - position[0], tile.position_m[1] - position[1]];
                        let relative_cosine = orientation[0] * other[0] + orientation[1] * other[1];
                        let relative_sine = orientation[0] * other[1] - orientation[1] * other[0];
                        let extent = 0.5 * dx * (1.0 + relative_cosine.abs() + relative_sine.abs());
                        // Separating axes of both actual squares. Touching tile
                        // edges are allowed; roundoff cannot invent finite area.
                        let tolerance = 64.0 * f64::EPSILON * extent;
                        let overlaps = [orientation, other].into_iter().all(|axis| {
                            let along = (delta[0] * axis[0] + delta[1] * axis[1]).abs();
                            let across = (-delta[0] * axis[1] + delta[1] * axis[0]).abs();
                            along < extent - tolerance && across < extent - tolerance
                        });
                        if overlaps {
                            return Err(sizing("observed reservoir surfaces must be disjoint in the common baffle"));
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

pub(super) fn admit_storage<'a>(
    receivers: impl IntoIterator<Item = &'a BaffledTilePressure>,
) -> Result<(), RenderError> {
    let limits = TilePressureLimits::surface_reference();
    let (mut tiles, mut histories, mut kernels) = (0usize, 0usize, 0usize);
    for receiver in receivers {
        let (history, kernel) = receiver.retained_values();
        tiles = tiles.checked_add(receiver.input_count())
            .ok_or_else(|| sizing("resolved observation tile count overflow"))?;
        histories = histories.checked_add(history)
            .ok_or_else(|| sizing("resolved observation history count overflow"))?;
        kernels = kernels.checked_add(kernel)
            .ok_or_else(|| sizing("resolved observation kernel count overflow"))?;
    }
    if tiles > limits.maximum_tiles || histories > limits.maximum_history_values
        || kernels > limits.maximum_kernel_values
    {
        return Err(sizing("observed reservoir apertures exceed their shared tile/storage allowance"));
    }
    Ok(())
}

impl<L: Lattice3<27>> NsfDgPerformance<L> {
    /// Add the escape aperture's original per-tile flow to the same microphone
    /// as the bore outlets. Physical pressure is signed superposition with
    /// each original position and travel time. The numerical derivative and
    /// delay owner is shared; a subsequent `with_observer_gas` applies its law
    /// to every path and aligns only extra numerical padding.
    ///
    /// The source and all its reaction/history stay unchanged. This is one-way
    /// baffled observation of a prescribed-reservoir boundary, not a matched
    /// exterior load, transparent boundary, acoustic/hydrodynamic separation or
    /// free-field flute qualification. The ghost must be quiescent in the same
    /// reference gas as the computed exterior; gas mixing is not inferred.
    /// Each physical grid tile keeps its own signed flow and path in bounded
    /// groups of original receivers; midpoint quadrature must resolve the
    /// declared band to pi/4 per tile.
    pub fn with_reservoir_observation(
        self, input: ReservoirObservation,
    ) -> Result<Self, RenderError> {
        if self.system.zone().reservoir_aperture_count() != 1 {
            return Err(sizing("select the observed reservoir explicitly when multiple apertures exist"));
        }
        self.with_reservoir_observation_at(0, input)
    }

    /// Add one explicitly selected original escape surface. Repeated calls
    /// preserve disjoint source identities, tile flows and propagation paths;
    /// an original aperture may appear only once. All observed surfaces share
    /// the existing aggregate tile/storage allowance and exterior material.
    /// Independently prescribed breath tiles are not selected implicitly.
    /// This does not add reaction to apertures without an original joint load.
    pub fn with_reservoir_observation_at(
        mut self, aperture: usize, input: ReservoirObservation,
    ) -> Result<Self, RenderError> {
        if self.poisoned
            || self.system.accepted_steps() != 0
            || self.reservoirs.iter().any(|observer| observer.aperture == aperture)
            || self.observer_gas.is_some()
            || self.gate.is_requested()
            || !input
                .origin_m
                .iter()
                .chain([input.angle_rad].iter())
                .all(|x| x.is_finite())
        {
            return Err(sizing(
                "reservoir observation needs a unique cold finite placement before gas preparation",
            ));
        }
        if self.system.reservoir_radiation_at(aperture).is_some_and(|reaction|
            reaction.origin_m != input.origin_m || input.angle_rad != 0.0)
        {
            return Err(sizing(
                "reservoir observation must retain the original joint radiation placement",
            ));
        }
        let ApertureObservation::TubeBaffled(location) = self.observation else {
            return Err(sizing(
                "reservoir observation requires the original exterior microphone",
            ));
        };
        let zone = self.system.zone();
        let opening = zone.reservoir_aperture_at(aperture).ok_or_else(|| {
            sizing("reservoir observation requires an original transverse aperture")
        })?;
        let source =
            self.system.load().radiation_record().ok_or_else(|| {
                sizing("reservoir observation requires the bore's computed exterior")
            })?;
        let exterior = source.gas_record().state();
        let near = |a: f64, b: f64| (a - b).abs() <= 32.0 * f64::EPSILON * a.abs().max(b.abs());
        if self.system.reference_gas() != exterior
            || !near(opening.state.density_kg_m3, exterior.density)
            || !near(opening.state.temperature_k, exterior.temperature)
            || opening.state.velocity_m_s != [0.0; 3]
        {
            return Err(sizing(
                "reservoir radiation requires one quiescent exterior and source reference gas",
            ));
        }
        let accepted = zone.last_reservoir_flows_at(aperture);
        let limits = TilePressureLimits::surface_reference();
        let existing_tiles = self.reservoirs.iter().map(|observer| observer.surface.len()).sum::<usize>();
        if accepted.is_empty() || accepted.len() > limits.maximum_tiles.saturating_sub(existing_tiles) {
            return Err(sizing(
                "reservoir radiation exceeds its original per-tile input budget",
            ));
        }
        let dx = zone.law().dx_m;
        let phase = core::f64::consts::TAU * location.maximum_frequency_hz * (2.0_f64.sqrt() * dx)
            / exterior.sound_speed;
        if !phase.is_finite() || phase > core::f64::consts::FRAC_PI_4 {
            return Err(sizing(
                "reservoir grid tiles underresolve the declared acoustic band",
            ));
        }
        let a = (opening.axis + 1) % 3;
        let b = (opening.axis + 2) % 3;
        let sign = if opening.high { 1.0 } else { -1.0 };
        let cosine = det::cos(input.angle_rad);
        let sine = det::sin(input.angle_rad);
        let mut surface = Vec::with_capacity(accepted.len());
        for flow in accepted {
            if self.gate.is_requested() {
                return Err(sizing("reservoir observation cancelled"));
            }
            let local = [
                (flow.cell[a] as f64 + 0.5) * dx,
                sign * (flow.cell[b] as f64 + 0.5) * dx,
            ];
            let center = [
                input.origin_m[0] + cosine * local[0] - sine * local[1],
                input.origin_m[1] + sine * local[0] + cosine * local[1],
            ];
            for outlet in source.apertures() {
                let d = [
                    outlet.center_m[0] - center[0],
                    outlet.center_m[1] - center[1],
                ];
                let nearest = [
                    (cosine * d[0] + sine * d[1]).abs() - 0.5 * dx,
                    (-sine * d[0] + cosine * d[1]).abs() - 0.5 * dx,
                ];
                if nearest[0].max(0.0).hypot(nearest[1].max(0.0)) < outlet.radius_m {
                    return Err(sizing(
                        "reservoir tiles and bore radiation apertures must not overlap",
                    ));
                }
            }
            surface.push(FlowSurfaceSample {
                position_m: [center[0], center[1], 0.0],
                area_m2: dx * dx,
            });
        }
        disjoint_surfaces(&self.reservoirs, &surface, input, dx, &mut || self.gate.is_requested())?;
        let receiver = BaffledTilePressure::new_differentiating(
            &surface,
            self.sample_rate_hz,
            location.position_m,
            exterior.into(),
            limits,
            &mut || self.gate.is_requested(),
        )
        .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?;
        admit_storage(self.reservoirs.iter().map(|observer| &observer.receiver)
            .chain(std::iter::once(&receiver)))?;
        self.reservoirs.push(ReservoirPressure {
            aperture,
            input,
            surface,
            receiver,
            gas: None,
            cells: accepted.iter().map(|flow| flow.cell).collect(),
            flows: vec![0.0; accepted.len()],
        });
        Ok(self)
    }

    /// Number of explicitly observed original escape apertures.
    pub fn reservoir_observation_count(&self) -> usize { self.reservoirs.len() }

    /// Explicit placement for an original aperture identity.
    pub fn reservoir_observation_at(&self, aperture: usize) -> Option<ReservoirObservation> {
        self.reservoirs.iter().find(|observer| observer.aperture == aperture).map(|observer| observer.input)
    }
    /// Actual shared Rayleigh bank for one original aperture identity.
    pub fn reservoir_receiver_at(&self, aperture: usize) -> Option<&BaffledTilePressure> {
        self.reservoirs.iter().find(|observer| observer.aperture == aperture).map(|observer| &observer.receiver)
    }
    /// Same material law and common numerical padding on a selected aperture.
    pub fn reservoir_observer_gas_at(&self, aperture: usize) -> Option<&GasPropagationReport> {
        self.reservoirs.iter().find(|observer| observer.aperture == aperture)
            .and_then(|observer| observer.gas.as_ref())
    }
    /// Singular compatibility view; multiple selections require explicit identity.
    pub fn reservoir_observation(&self) -> Option<ReservoirObservation> {
        let [observer] = self.reservoirs.as_slice() else { return None; };
        Some(observer.input)
    }
    /// Singular compatibility view of the original observed aperture.
    pub fn reservoir_observation_aperture(&self) -> Option<usize> {
        let [observer] = self.reservoirs.as_slice() else { return None; };
        Some(observer.aperture)
    }
    /// Singular compatibility view of the actual receiver.
    pub fn reservoir_receiver(&self) -> Option<&BaffledTilePressure> {
        let [observer] = self.reservoirs.as_slice() else { return None; };
        Some(&observer.receiver)
    }
    /// Singular compatibility view of the original exterior propagation report.
    pub fn reservoir_observer_gas(&self) -> Option<&GasPropagationReport> {
        let [observer] = self.reservoirs.as_slice() else { return None; };
        observer.gas.as_ref()
    }
}
