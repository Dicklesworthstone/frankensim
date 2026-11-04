//! Material-owned pressure waves in the original branched network solver.
//!
//! One admitted fluid snapshot supplies every edge and material boundary at
//! the queried frequency. Actual density and sound speed set propagation phase
//! and impedance; thermoviscous calls retain the original viscosity,
//! conductivity, Cp and Cv. The source is never converted to a gas surrogate.
//!
//! Cylinders, exact rigid lossless cones, refinable lossy cones and authored
//! profiles use the existing material chain and stable scattering composer.
//! Junctions, cycles, parallel edges, supplied multiports, coherent pressure /
//! flow sources, mutual source impedances and complete coupled loads use the
//! original retained network equations. An open-circuit impedance pole need
//! not prevent a finite loaded response. Prepared application and source
//! statistics reuse the existing operator without another material query.
//!
//! The graph contains one uniform, stationary fluid, also filling its radiation
//! exterior and porous loads. A liquid free surface, mean flow and moving phase
//! interfaces are not inferred. Embedded gas tone-hole shorthand still refuses:
//! model those chimneys as explicit graph branches with their own boundaries.
//! Both compact junction-correction policies are equivalent on these admitted
//! hole-free edges. This is the original plane-wave/Webster graph rung, not a
//! claim of resolved transverse or three-dimensional junction fields.

use super::{
    CoupledSourceBank, CoupledTermination, DuctNetwork, MAX_NETWORK_CONDITION,
    NetworkSourceResponse, PortSource, PortTermination, PreparedNetworkSources, bad,
    finite_impedance, loaded,
};
use crate::compose::{
    ComposeOptions, MAX_COMPOSITION_BLOCKS, PreparedChain, TwoPort, compose_scattering,
};
use crate::fluid::{FluidDuctError, boundary_from_state};
use crate::{ConeSlices, DuctError, LossModel, SegmentWave};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

pub(super) fn admit(net: &DuctNetwork, omega: f64) -> Result<(), FluidDuctError> {
    net.validate()?;
    if !omega.is_finite() || omega <= 0.0 {
        return Err(bad("material network frequency must be positive and finite").into());
    }
    Ok(())
}

/// Prepare coherent physical sources from one original acoustic material.
/// Source laws and explicit scalar/coupled loads partition every external
/// port exactly once. An optional mutual-impedance source bank retains its
/// complete signed coupling in the same solve.
///
/// Rigid lossless cones are exact; profile step counts remain caller-authored.
/// The graph's original size/condition limits apply, and aggregate propagation
/// work is bounded by [`MAX_COMPOSITION_BLOCKS`] across all edges. This is the
/// existing synchronous dense solve. No open-circuit impedance is formed.
///
/// # Errors
/// Original material, geometry, propagation, source and load refusals. A porous
/// load needs [`prepare_sources_with_loss`] and its complete fluid snapshot,
/// even if the propagation loss selected there is [`LossModel::Lossless`].
#[allow(clippy::too_many_arguments)]
pub fn prepare_sources<M: FluidAcoustics + ?Sized>(
    net: &DuctNetwork,
    material: &M,
    omega: f64,
    sources: &[PortSource],
    bank: Option<&CoupledSourceBank>,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<PreparedNetworkSources, FluidDuctError> {
    admit(net, omega)?;
    let state = MaterialSnapshot::Acoustic(material.acoustic_state()?);
    prepare_snapshot(
        net,
        &state,
        omega,
        LossModel::Lossless,
        ConeSlices::Auto,
        sources,
        bank,
        loads,
        coupled,
    )
}

/// Prepare the same graph with the original material's wall diffusion and
/// porous-boundary thermodynamics. The cone refinement axis is explicit;
/// no geometry step or original material coefficient is silently capped.
/// One complete snapshot supplies both propagation and all material loads.
///
/// # Errors
/// As [`prepare_sources`], plus original transport/caloric, wall-diffusion and
/// cone-resolution refusals. Radiation retains its original model and range.
#[allow(clippy::too_many_arguments)]
pub fn prepare_sources_with_loss<M: ThermoviscousFluid + ?Sized>(
    net: &DuctNetwork,
    material: &M,
    omega: f64,
    loss: LossModel,
    cone_slices: ConeSlices,
    sources: &[PortSource],
    bank: Option<&CoupledSourceBank>,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<PreparedNetworkSources, FluidDuctError> {
    admit(net, omega)?;
    let state = MaterialSnapshot::Thermoviscous(material.thermoviscous_state()?);
    prepare_snapshot(
        net,
        &state,
        omega,
        loss,
        cone_slices,
        sources,
        bank,
        loads,
        coupled,
    )
}

/// Actual coherent port pressures, flows and work for the supplied material
/// and source amplitudes. The existing signed generator and passive source
/// loss accounting includes all interference and mutual-load terms.
/// Use [`prepare_sources`] to reuse this boundary for many source vectors.
///
/// # Errors
/// As [`prepare_sources`], including a true pole of the fully loaded graph.
#[allow(clippy::too_many_arguments)]
pub fn drive_sources<M: FluidAcoustics + ?Sized>(
    net: &DuctNetwork,
    material: &M,
    omega: f64,
    sources: &[PortSource],
    bank: Option<&CoupledSourceBank>,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkSourceResponse, FluidDuctError> {
    drive_default(prepare_sources(
        net, material, omega, sources, bank, loads, coupled,
    )?)
}

/// Actual material-driven port fields and work with the same original wall
/// and boundary laws as [`prepare_sources_with_loss`].
///
/// # Errors
/// Original material, geometry, propagation, source and load refusals.
#[allow(clippy::too_many_arguments)]
pub fn drive_sources_with_loss<M: ThermoviscousFluid + ?Sized>(
    net: &DuctNetwork,
    material: &M,
    omega: f64,
    loss: LossModel,
    cone_slices: ConeSlices,
    sources: &[PortSource],
    bank: Option<&CoupledSourceBank>,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<NetworkSourceResponse, FluidDuctError> {
    drive_default(prepare_sources_with_loss(
        net,
        material,
        omega,
        loss,
        cone_slices,
        sources,
        bank,
        loads,
        coupled,
    )?)
}

fn drive_default(
    prepared: PreparedNetworkSources,
) -> Result<NetworkSourceResponse, FluidDuctError> {
    let mut response = prepared.response_buffer();
    prepared.drive_into(prepared.default_excitations(), &mut response)?;
    Ok(response)
}

/// One immutable original query shared by graph, spatial and spectral owners.
#[derive(Clone, Copy)]
pub(super) enum MaterialSnapshot {
    Acoustic(FluidAcousticState),
    Thermoviscous(ThermoviscousFluidState),
}

impl MaterialSnapshot {
    pub(super) fn acoustic(self) -> FluidAcousticState {
        match self {
            Self::Acoustic(state) => state,
            Self::Thermoviscous(state) => state.acoustic(),
        }
    }

    fn thermoviscous(&self) -> Option<&ThermoviscousFluidState> {
        match self {
            Self::Acoustic(_) => None,
            Self::Thermoviscous(state) => Some(state),
        }
    }

    pub(super) fn wave(
        &self,
        radius: f64,
        omega: f64,
        loss: LossModel,
    ) -> Result<SegmentWave, DuctError> {
        match self {
            Self::Thermoviscous(state) => {
                crate::fluid::thermoviscous_wave_from_state(state, radius, omega, loss)
            }
            Self::Acoustic(state) if loss == LossModel::Lossless => {
                crate::fluid::wave_from_state(state, radius, omega)
            }
            Self::Acoustic(_) => Err(bad(
                "wall diffusion requires the original thermoviscous material",
            )),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn chain(
        &self,
        edge: &super::NetworkEdge,
        omega: f64,
        loss: LossModel,
        cones: ConeSlices,
        options: ComposeOptions,
        boundaries: Option<&mut Vec<usize>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedChain, DuctError> {
        crate::fluid::duct_chain_from_snapshot(
            &edge.duct,
            &self.acoustic(),
            self.thermoviscous(),
            omega,
            loss,
            cones,
            options,
            boundaries,
            cancelled,
        )
    }

    fn boundary(
        &self,
        omega: f64,
        load: &super::NetworkLoad,
    ) -> Result<Option<fs_math::c64::C64>, FluidDuctError> {
        Ok(
            boundary_from_state(&self.acoustic(), self.thermoviscous(), omega, load)?
                .map(finite_impedance)
                .transpose()?,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn edge_operator(
        &self,
        edge: &super::NetworkEdge,
        omega: f64,
        loss: LossModel,
        cones: ConeSlices,
        blocks: &mut usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(TwoPort, f64), DuctError> {
        let options = ComposeOptions {
            maximum_condition: MAX_NETWORK_CONDITION,
            maximum_blocks: MAX_COMPOSITION_BLOCKS - *blocks,
            ..ComposeOptions::default()
        };
        let chain = self.chain(edge, omega, loss, cones, options, None, &mut cancelled)?;
        *blocks += chain.blocks.len();
        compose_scattering(&chain.blocks, options, &mut cancelled)
    }
}

#[allow(clippy::too_many_arguments)]
fn prepare_snapshot(
    net: &DuctNetwork,
    state: &MaterialSnapshot,
    omega: f64,
    loss: LossModel,
    cone_slices: ConeSlices,
    sources: &[PortSource],
    bank: Option<&CoupledSourceBank>,
    loads: &[PortTermination],
    coupled: Option<&CoupledTermination>,
) -> Result<PreparedNetworkSources, FluidDuctError> {
    let mut blocks = 0usize;
    loaded::prepare_sources_with_operators(
        net,
        omega,
        sources,
        bank,
        loads,
        coupled,
        |edge| Ok(state.edge_operator(edge, omega, loss, cone_slices, &mut blocks, || false)?),
        |load, _node| state.boundary(omega, load),
    )
}

pub(super) fn prepare_spatial_snapshot(
    net: &DuctNetwork,
    state: &MaterialSnapshot,
    cones: ConeSlices,
    request: &super::field::NetworkFieldRequest<'_>,
    selected: &[usize],
    mut cancelled: impl FnMut() -> bool,
) -> Result<(PreparedNetworkSources, Vec<[fs_math::c64::C64; 4]>), FluidDuctError> {
    let mut blocks = 0usize;
    loaded::prepare_spatial_with_operators(
        net,
        request.omega,
        request.sources,
        request.source_bank,
        request.loads,
        request.coupled_load,
        selected,
        |edge| {
            Ok(state.edge_operator(
                edge,
                request.omega,
                request.loss,
                cones,
                &mut blocks,
                &mut cancelled,
            )?)
        },
        |load, _node| state.boundary(request.omega, load),
    )
}

#[cfg(test)]
mod tests;
