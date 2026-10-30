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

    /// Observe one explicitly selected original escape surface. Independently
    /// prescribed breath tiles are retained by the fluid, not mixed into it.
    pub fn with_reservoir_observation_at(
        mut self, aperture: usize, input: ReservoirObservation,
    ) -> Result<Self, RenderError> {
        if self.poisoned
            || self.system.accepted_steps() != 0
            || self.reservoir.is_some()
            || self.observer_gas.is_some()
            || self.gate.is_requested()
            || !input
                .origin_m
                .iter()
                .chain([input.angle_rad].iter())
                .all(|x| x.is_finite())
        {
            return Err(sizing(
                "reservoir observation needs one cold finite placement before gas preparation",
            ));
        }
        if self.system.reservoir_radiation().is_some_and(|reaction|
            self.system.reservoir_radiation_aperture() != Some(aperture)
                || reaction.origin_m != input.origin_m || input.angle_rad != 0.0)
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
        if accepted.is_empty() || accepted.len() > limits.maximum_tiles {
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
        let receiver = BaffledTilePressure::new_differentiating(
            &surface,
            self.sample_rate_hz,
            location.position_m,
            exterior.into(),
            limits,
            &mut || self.gate.is_requested(),
        )
        .map_err(|error| RenderError::Voice(AcousticRealizeError::Nonlinear(error)))?;
        self.reservoir = Some(ReservoirPressure {
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

    /// Explicit placement of the original aperture, if selected.
    pub fn reservoir_observation(&self) -> Option<ReservoirObservation> {
        self.reservoir.as_ref().map(|observer| observer.input)
    }
    /// Stable original owner of the observed escape opening.
    pub fn reservoir_observation_aperture(&self) -> Option<usize> {
        self.reservoir.as_ref().map(|observer| observer.aperture)
    }
    /// Actual shared Rayleigh receiver bank and its separate physical paths.
    pub fn reservoir_receiver(&self) -> Option<&BaffledTilePressure> {
        self.reservoir.as_ref().map(|observer| &observer.receiver)
    }
    /// Same material propagation law, applied on the original escape paths.
    pub fn reservoir_observer_gas(&self) -> Option<&GasPropagationReport> {
        self.reservoir
            .as_ref()
            .and_then(|observer| observer.gas.as_ref())
    }
}
