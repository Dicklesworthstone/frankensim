//! Unresolved vortex-tangle memory on the channel's interior dual cells.

use super::*;
use fs_material::superfluid::vinen::VinenTangle;
use fs_material::superfluid::{MutualFrictionStep, TwoFluidVelocity};

mod transport;
pub(super) use transport::TransportChange;
use transport::VortexTransportMemory;
pub use transport::{VortexTransportBoundaries, VortexTransportBoundary, VortexTransportControls};

#[derive(Debug)]
pub(super) struct VinenTangleMemory {
    model: VinenTangle,
    accepted: Vec<f64>,
    pending: Vec<f64>,
    transport: Option<VortexTransportMemory>,
}

impl VinenTangleMemory {
    fn new(
        model: VinenTangle,
        initial: &[f64],
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ChannelError> {
        let mut memory = Self {
            model,
            accepted: filled(initial.len())?,
            pending: filled(initial.len())?,
            transport: None,
        };
        for (i, &density) in initial.iter().enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            memory.model.line_energy_density_j_m3(density)?;
            memory.accepted[i] = density;
            memory.pending[i] = density;
        }
        memory.energy(false, volume, cancelled)?;
        Ok(memory)
    }

    pub(super) fn copy_pending(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        for (i, (pending, &accepted)) in self.pending.iter_mut().zip(&self.accepted).enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            *pending = accepted;
        }
        Ok(())
    }

    /// Each material call has its own bounded substep budget. Check cancellation
    /// before EVERY call rather than hiding up to 128 source solves in one tile.
    pub(super) fn advance(
        &mut self,
        face: usize,
        velocity: TwoFluidVelocity,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MutualFrictionStep, ChannelError> {
        poll(cancelled)?;
        let density = &mut self.pending[face - 1];
        let change = self.model.advance(velocity, *density, dt)?;
        *density = change.line_density_per_m2;
        Ok(MutualFrictionStep {
            velocity: change.velocity,
            normal_impulse_kg_m2_s: change.normal_impulse_kg_m2_s,
            superfluid_impulse_kg_m2_s: change.superfluid_impulse_kg_m2_s,
            dissipated_energy_j_m3: change.dissipated_energy_j_m3,
        })
    }

    pub(super) fn energy(
        &self,
        pending: bool,
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ChannelError> {
        let densities = if pending {
            &self.pending
        } else {
            &self.accepted
        };
        let mut total = 0.0;
        let mut correction = 0.0;
        for (i, &density) in densities.iter().enumerate() {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let energy = self.model.line_energy_density_j_m3(density)? * volume;
            accumulate(&mut total, &mut correction, energy);
        }
        if total.is_finite() && total >= 0.0 {
            Ok(total)
        } else {
            Err(ChannelError::Invalid(
                "vortex-line storage is not representable",
            ))
        }
    }

    pub(super) fn publish(&mut self) {
        core::mem::swap(&mut self.accepted, &mut self.pending);
    }
}

impl TwoFluidChannel {
    /// Attach a spatially varying local Vinen tangle before the first step.
    ///
    /// `interior_line_density_per_m2[i]` is the total line length per volume
    /// on interior face `i+1`, whose dual volume is `area*dx`. There must be
    /// exactly N-1 values; endpoint half-cell vortex storage belongs to a
    /// separately authored boundary model. Model and channel must have the
    /// same frozen material. Reconfiguration is allowed only before stepping.
    ///
    /// Growth, decay and section-mean mutual drag then require the explicit
    /// `EvolvingVinenTangle` selector on EVERY step and endpoint-response query.
    /// Other friction selectors refuse, preventing duplicate drag or silently
    /// frozen memory. The source includes its authored effective line energy
    /// in the channel balance and exports only the remaining dissipative heat.
    ///
    /// Optional authored line transport is configured separately and retained
    /// when reseeding before the first step. The isotropic closure has no
    /// orientation, resolved reconnections, or nucleation from an empty tangle.
    /// Normal-profile deviations remain unchanged by its section-mean force.
    /// Allocation occurs here only; cancellation/refusal preserves all accepted
    /// fields, profile and previously configured tangle memory.
    pub fn configure_vinen_tangle(
        &mut self,
        model: VinenTangle,
        interior_line_density_per_m2: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        poll(&mut cancelled)?;
        if self.time_s != 0.0
            || interior_line_density_per_m2.len() != self.cells() - 1
            || *model.material() != self.material
        {
            return Err(ChannelError::Invalid(
                "Vinen configuration needs N-1 interior densities and the same frozen material before stepping",
            ));
        }
        let mut memory = VinenTangleMemory::new(
            model,
            interior_line_density_per_m2,
            self.cell_volume,
            &mut cancelled,
        )?;
        let total = self.energy_with_profile(false, &mut cancelled)?
            + memory.energy(false, self.cell_volume, &mut cancelled)?;
        if !total.is_finite() {
            return Err(ChannelError::Invalid(
                "initial wave and vortex storage is not representable",
            ));
        }
        if let Some(transport) = self
            .vinen_tangle
            .as_ref()
            .and_then(|old| old.transport.as_ref())
        {
            transport.validate_boundary_densities(&model)?;
        }
        poll(&mut cancelled)?;
        if let Some(old) = &mut self.vinen_tangle {
            memory.transport = old.transport.take();
        }
        self.vinen_tangle = Some(memory);
        Ok(())
    }

    /// The authored local tangle law, if configured.
    pub fn vinen_tangle(&self) -> Option<&VinenTangle> {
        self.vinen_tangle.as_ref().map(|memory| &memory.model)
    }

    /// Accepted line length per volume [1/m²] at an interior face, numbered
    /// 1..N-1. Endpoints, out-of-range faces and unconfigured tangles return None.
    pub fn vortex_line_density_per_m2(&self, face: usize) -> Option<f64> {
        self.vinen_tangle
            .as_ref()?
            .accepted
            .get(face.checked_sub(1)?)
            .copied()
    }

    /// Authored effective vortex-line energy on the interior dual volumes [J].
    /// Separate from wave availability and included once in represented_energy_j.
    pub fn vortex_energy_j(&self) -> Result<f64, ChannelError> {
        self.vortex_energy_with(false, &mut || false)
    }

    pub(super) fn vortex_energy_with(
        &self,
        pending: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ChannelError> {
        match &self.vinen_tangle {
            Some(memory) => memory.energy(pending, self.cell_volume, cancelled),
            None => Ok(0.0),
        }
    }

    pub(super) fn transport_vortex_half_step(
        &mut self,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportChange, ChannelError> {
        match &mut self.vinen_tangle {
            Some(memory) => match &mut memory.transport {
                Some(transport) => transport.advance(
                    &mut memory.pending,
                    dt,
                    self.dx,
                    self.cell_volume,
                    cancelled,
                ),
                None => Ok(TransportChange::default()),
            },
            None => Ok(TransportChange::default()),
        }
    }
}
