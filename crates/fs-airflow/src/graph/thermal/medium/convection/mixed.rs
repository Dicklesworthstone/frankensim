//! Aiding mixed convection on the original material/correlation owners.

use super::{
    CorrelationId, CorrelationInputs, FluidTransport, FluidTransportState, Length,
    MaterialConvection, MaterialConvectionError, Temperature, finish,
};
pub use fs_material::fluid_transport::{AidingVerticalPlateFrame, AidingVerticalPlateSettings};

pub(super) fn inputs(
    transport: FluidTransportState,
    settings: AidingVerticalPlateSettings,
    expansion_per_k: f64,
    wall: Temperature,
    reference: Temperature,
    length: Length,
) -> Result<CorrelationInputs, MaterialConvectionError> {
    let groups = settings.groups(transport, expansion_per_k, wall, reference, length)?;
    Ok(CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl))
}

/// Evaluate one laminar aiding mixed-convection card from one material snapshot.
///
/// Length is plate height along the free stream. The same snapshot supplies Re,
/// Ra, Pr and the conductivity in h=Nu*k/L. Signed expansion, actual wall/bulk
/// temperatures and signed flow enforce aiding alignment; slowing the flow to
/// zero continuously reaches this card's natural-convection limit.
///
/// The source relation assumes constant properties. The source-cell and gas
/// adapters declare a mean-film frozen-property approximation and re-query that
/// state inside their existing thermal residual. Film averaging is an adapter
/// policy, not a property rule asserted by the original paper. The caller must
/// also admit property variation; the density bound alone does not establish it.
///
/// # Errors
/// Refuses opposing flow, invalid states, excessive declared density departure,
/// unsupported laminar/Prandtl ranges or vanishing forced AND buoyant driving.
/// A thermal owner can represent exact resting equilibrium as zero heat without
/// evaluating a zero/zero correlation.
pub fn evaluate_aiding_vertical_plate_from_material<F: FluidTransport + ?Sized>(
    fluid: &F,
    settings: AidingVerticalPlateSettings,
    expansion_per_k: f64,
    wall: Temperature,
    reference: Temperature,
    characteristic_length: Length,
) -> Result<MaterialConvection<'_, F>, MaterialConvectionError> {
    let transport = fluid.transport_state()?;
    let inputs = inputs(
        transport,
        settings,
        expansion_per_k,
        wall,
        reference,
        characteristic_length,
    )?;
    finish(
        CorrelationId::ManningQureshiAidingVerticalPlate,
        fluid,
        transport,
        characteristic_length,
        inputs,
    )
}
