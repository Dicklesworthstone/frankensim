//! Complete vibrating bodies at a common-velocity attachment.
//!
//! Composition retains the original constitutive owner and accepted histories.
//! Connecting an already moving body requires compatible attachment velocity;
//! an impulsive projection is a different physical operation and is not inferred.

use super::AcousticRealizeError;
use crate::thin_plate::{PlateBankImage, PlateBankObservation};
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, join_port, transformer};

pub(super) struct DiracBodyLoad {
    pub(super) system: PortHamiltonian,
    pub(super) state: Vec<f64>,
    /// Pressure/volume port connected to the body's outward surface.
    pub(super) face: usize,
    /// Optional independent fluid source in the original load's port order.
    pub(super) drive_port: Option<usize>,
}

pub(super) struct DiracBodyJoin {
    pub(super) system: DescriptorPortHamiltonian,
    pub(super) state: Vec<f64>,
    pub(super) observation: DiracBodyObservation,
    pub(super) fluid_port: Option<usize>,
}

pub(super) struct DiracBodyObservation {
    offset: usize,
    extent: usize,
    original: PlateBankObservation,
}

impl DiracBodyObservation {
    pub(super) fn pressure_between(
        &self,
        before: &[f64],
        after: &[f64],
        dt: f64,
        rho: f64,
        distance_m: f64,
    ) -> Result<f64, AcousticRealizeError> {
        if before.len() != self.extent || after.len() != self.extent {
            return Err(invalid(
                "Dirac observation requires the complete joined state",
            ));
        }
        let end = self.offset + self.original.state_dim();
        self.original.pressure_between(
            &before[self.offset..end],
            &after[self.offset..end],
            dt,
            rho,
            distance_m,
        )
    }
}

/// Original state order is `[waveguide, complete bank, optional load, lambda]`.
/// Every unused waveguide port remains first, in its original order. The bank's
/// independent pressure port stays open with zero applied pressure if no load is
/// supplied; it is never mislabeled as an external volume-flow source.
pub(super) fn join_complete_body(
    waveguide: PortHamiltonian,
    waveguide_state: Vec<f64>,
    image: PlateBankImage,
    load: Option<DiracBodyLoad>,
) -> Result<DiracBodyJoin, AcousticRealizeError> {
    let (body, body_state, observation) = image.into_parts();
    let (string_velocity, string_scale) = checked_port(&waveguide, &waveguide_state, 0)?;
    let (body_velocity, body_scale) = checked_port(&body, &body_state, 0)?;
    let scale = string_scale + body_scale;
    let allowance = 64.0 * f64::EPSILON * (waveguide.state_dim() + body.state_dim()) as f64 * scale;
    if !scale.is_finite()
        || !allowance.is_finite()
        || (string_velocity - body_velocity).abs() > allowance
    {
        return Err(invalid(
            "Dirac initial attachment velocities disagree; an impulse or compatible history is required",
        ));
    }
    let offset = waveguide.state_dim();
    let external_string_ports = waveguide.port_dim() - 1;
    let mut initial = waveguide_state;
    initial.extend_from_slice(&body_state);
    let (body, fluid_port) = if let Some(load) = load {
        checked_state(&load.system, &load.state)?;
        if load.face >= load.system.port_dim()
            || load
                .drive_port
                .is_some_and(|p| p >= load.system.port_dim() || p == load.face)
        {
            return Err(invalid(
                "Dirac fluid load has an invalid surface or source port",
            ));
        }
        let source = load
            .drive_port
            .map(|p| external_string_ports + p - usize::from(p > load.face));
        initial.extend_from_slice(&load.state);
        (
            transformer(body, load.system, 1, load.face, 1.0)
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
            source,
        )
    } else {
        (body, None)
    };
    let system = join_port(waveguide, body, 0, 0)
        .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
    // The multiplier carries no energy. Zero is an initial representative;
    // the selected integrator supplies its reaction (instantaneous endpoint
    // for exact quadratic flow, interval average for the descriptor step).
    initial.push(0.0);
    if initial.len() != system.state_dim() {
        return Err(invalid(
            "Dirac composition changed the complete state extent",
        ));
    }
    let observation = DiracBodyObservation {
        offset,
        extent: system.state_dim(),
        original: observation,
    };
    Ok(DiracBodyJoin {
        system,
        state: initial,
        observation,
        fluid_port,
    })
}

fn checked_state(system: &PortHamiltonian, state: &[f64]) -> Result<(), AcousticRealizeError> {
    if state.len() != system.state_dim()
        || state.iter().any(|x| !x.is_finite())
        || !system.hamiltonian(state).is_finite()
    {
        return Err(invalid(
            "Dirac composition requires complete finite physical history",
        ));
    }
    Ok(())
}

fn checked_port(
    system: &PortHamiltonian,
    state: &[f64],
    port: usize,
) -> Result<(f64, f64), AcousticRealizeError> {
    checked_state(system, state)?;
    if port >= system.port_dim() {
        return Err(invalid("Dirac attachment port is absent"));
    }
    let effort = system.effort(state);
    if effort.iter().any(|e| !e.is_finite()) {
        return Err(invalid("Dirac attachment effort is nonfinite"));
    }
    let (_, _, g) = system.structure();
    let mut velocity = 0.0;
    let mut scale = 0.0;
    for (i, &e) in effort.iter().enumerate() {
        let value = g[i * system.port_dim() + port] * e;
        velocity += value;
        scale += value.abs();
    }
    if !velocity.is_finite() || !scale.is_finite() {
        return Err(invalid("Dirac attachment velocity is nonfinite"));
    }
    Ok((velocity, scale))
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

#[cfg(test)]
mod tests;
