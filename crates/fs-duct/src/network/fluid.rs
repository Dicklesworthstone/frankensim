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
use crate::compose::{ComposeOptions, MAX_COMPOSITION_BLOCKS, compose_scattering};
use crate::fluid::{FluidDuctError, boundary_from_state};
use crate::{ConeSlices, LossModel, WallModel};
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

fn admit(net: &DuctNetwork, omega: f64) -> Result<(), FluidDuctError> {
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
    let state = material.acoustic_state()?;
    prepare_snapshot(
        net,
        &state,
        None,
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
    let state = material.thermoviscous_state()?;
    prepare_snapshot(
        net,
        &state.acoustic(),
        Some(&state),
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

#[allow(clippy::too_many_arguments)]
fn prepare_snapshot(
    net: &DuctNetwork,
    acoustic: &FluidAcousticState,
    thermoviscous: Option<&ThermoviscousFluidState>,
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
        |edge| {
            let options = ComposeOptions {
                maximum_condition: MAX_NETWORK_CONDITION,
                maximum_blocks: MAX_COMPOSITION_BLOCKS - blocks,
                ..ComposeOptions::default()
            };
            let chain = if let Some(state) = thermoviscous {
                crate::fluid::duct_chain_with_loss(
                    &edge.duct,
                    state,
                    omega,
                    loss,
                    WallModel::Rigid,
                    cone_slices,
                    options,
                    || false,
                )?
            } else {
                crate::fluid::duct_chain(
                    &edge.duct,
                    acoustic,
                    omega,
                    WallModel::Rigid,
                    cone_slices,
                    options,
                    || false,
                )?
            };
            blocks += chain.blocks.len();
            Ok(compose_scattering(&chain.blocks, options, &mut || false)?)
        },
        |load, _node| {
            Ok(boundary_from_state(acoustic, thermoviscous, omega, load)?
                .map(finite_impedance)
                .transpose()?)
        },
    )
}

#[cfg(test)]
mod tests;
