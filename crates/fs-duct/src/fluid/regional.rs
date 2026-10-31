//! Original fluid materials joined through the existing pressure/flow owner.
//!
//! Each authored axial segment has one original acoustic or thermoviscous
//! material. All blocks share one numerical power-wave reference, while their
//! actual impedances and diffusion use their own material. The existing stable
//! composer therefore enforces pressure and volume-flow continuity, including
//! reflection at a material interface; coefficients are never averaged.
//!
//! These are prescribed frozen regions at mechanical equilibrium, which the
//! caller must establish. Acoustic response alone does not certify background
//! pressure balance, miscibility or phase coexistence. There is no mean flow,
//! mixing, heat/species exchange, moving interface or capillary interface law.
//! Source material/query receipts remain with the caller. Interfaces coincide
//! with authored segment ends; split geometry explicitly to place a new one.

use super::{
    FluidDuctError, admit_blocks, admit_geometry, admit_load, bad, chain_from_waves_at_reference,
    poll, thermoviscous_wave_from_state, wave_from_state,
};
use crate::compose::{self, ChainResult, ComposeOptions, PreparedChain};
use crate::{ConeSlices, Duct, DuctError, LossModel, SegmentWave, WallModel};
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;

fn admit(
    duct: &Duct,
    count: usize,
    options: ComposeOptions,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DuctError> {
    admit_geometry(duct, options, cancelled)?;
    if count != duct.segments.len() {
        return Err(bad(
            "regional fluid duct needs one material per axial segment",
        ));
    }
    Ok(())
}

/// Lossless pressure/flow blocks with one material per authored axial segment.
/// Every source is evaluated once, after geometry/count/cancellation admission.
/// References to different source types may share a `&dyn FluidAcoustics` slice;
/// no conversion to an ideal gas or synthesized interface material is needed.
/// Cylinders, cones, profiles and local wall laws retain their original owners.
#[allow(clippy::too_many_arguments)]
pub fn duct_chain<M: FluidAcoustics>(
    duct: &Duct,
    materials: &[M],
    omega: f64,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedChain, FluidDuctError> {
    admit(duct, materials.len(), options, &mut cancelled)?;
    let mut states = Vec::with_capacity(materials.len());
    for (index, material) in materials.iter().enumerate() {
        poll(&mut cancelled, index)?;
        states.push(material.acoustic_state()?);
    }
    Ok(assemble(
        duct,
        &states,
        omega,
        wall,
        cone_slices,
        options,
        true,
        &mut cancelled,
        &mut |index, radius| wave_from_state(&states[index], radius, omega),
    )?)
}

/// Regional propagation and original local viscosity/heat-diffusion response.
/// Each source supplies one coherent acoustic/transport/Cp-Cv evaluation; all
/// profile nodes and cone slices in its segment use that same frozen state.
/// The same global block cap covers every region, slice and profile step.
#[allow(clippy::too_many_arguments)]
pub fn duct_chain_with_loss<M: ThermoviscousFluid>(
    duct: &Duct,
    materials: &[M],
    omega: f64,
    loss: LossModel,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedChain, FluidDuctError> {
    admit(duct, materials.len(), options, &mut cancelled)?;
    let mut states = Vec::with_capacity(materials.len());
    for (index, material) in materials.iter().enumerate() {
        poll(&mut cancelled, index)?;
        states.push(material.thermoviscous_state()?);
    }
    let acoustic: Vec<_> = states.iter().map(|state| state.acoustic()).collect();
    Ok(assemble(
        duct,
        &acoustic,
        omega,
        wall,
        cone_slices,
        options,
        loss == LossModel::Lossless,
        &mut cancelled,
        &mut |index, radius| thermoviscous_wave_from_state(&states[index], radius, omega, loss),
    )?)
}

#[allow(clippy::too_many_arguments)]
fn assemble(
    duct: &Duct,
    states: &[FluidAcousticState],
    omega: f64,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    lossless_bulk: bool,
    cancelled: &mut impl FnMut() -> bool,
    wave_at: &mut impl FnMut(usize, f64) -> Result<SegmentWave, DuctError>,
) -> Result<PreparedChain, DuctError> {
    // Both arrays are nonempty and equally sized after admission. This is a
    // numerical coordinate choice only, never a replacement physical material.
    let reference = wave_from_state(&states[0], duct.segments[0].inlet_radius(), omega)?
        .characteristic_impedance
        .re;
    let mut blocks = Vec::new();
    let mut min_shear_number = f64::INFINITY;
    for (index, (segment, state)) in duct.segments.iter().zip(states).enumerate() {
        poll(cancelled, index)?;
        admit_blocks(blocks.len(), 1, options.maximum_blocks)?;
        let remaining = ComposeOptions {
            maximum_blocks: options.maximum_blocks - blocks.len(),
            ..options
        };
        let local = chain_from_waves_at_reference(
            std::slice::from_ref(segment),
            state,
            omega,
            wall,
            cone_slices,
            remaining,
            lossless_bulk,
            reference,
            cancelled,
            &mut |radius| wave_at(index, radius),
        )?;
        min_shear_number = min_shear_number.min(local.min_shear_number);
        blocks.extend(local.blocks);
    }
    Ok(PreparedChain {
        blocks,
        reference_impedance: reference,
        min_shear_number,
    })
}

/// Terminated regional pressure and volume-flow response in exp(-i omega t).
/// The explicit passive load uses the actual terminal fluid/geometry owner's
/// Pa s/m3 convention. `None` is a rigid cap; zero is ideal pressure release.
/// The load may come from an independently admitted radiation/porous owner.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_load<M: FluidAcoustics>(
    duct: &Duct,
    materials: &[M],
    omega: f64,
    load: Option<C64>,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ChainResult, FluidDuctError> {
    admit_load(load)?;
    let chain = duct_chain(
        duct,
        materials,
        omega,
        wall,
        cone_slices,
        options,
        &mut cancelled,
    )?;
    Ok(compose::compose_chain_controlled(
        &chain.blocks,
        load,
        options,
        &mut cancelled,
    )?)
}

/// The same terminated regional response with each original fluid's boundary
/// diffusion. Local wall dissipation and terminal absorption remain separate
/// contributions in the same pressure/volume-flow solve.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_load_with_loss<M: ThermoviscousFluid>(
    duct: &Duct,
    materials: &[M],
    omega: f64,
    loss: LossModel,
    load: Option<C64>,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ChainResult, FluidDuctError> {
    admit_load(load)?;
    let chain = duct_chain_with_loss(
        duct,
        materials,
        omega,
        loss,
        wall,
        cone_slices,
        options,
        &mut cancelled,
    )?;
    Ok(compose::compose_chain_controlled(
        &chain.blocks,
        load,
        options,
        &mut cancelled,
    )?)
}
