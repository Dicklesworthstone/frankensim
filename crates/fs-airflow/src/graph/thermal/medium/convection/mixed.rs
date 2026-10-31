//! Aiding mixed convection on the original material/correlation owners.

use super::{
    CorrelationId, CorrelationInputs, FluidTransport, FluidTransportState, Length,
    MaterialConvection, MaterialConvectionError, MaterialError, Temperature, finish,
    natural_rayleigh,
};

/// Prescribed vertical free-stream motion and an explicit Boussinesq bound.
///
/// The Manning-Qureshi relation is for an unconfined isothermal vertical plate
/// in laminar AIDING flow. Positive speed is opposite gravity; negative speed
/// is along gravity. The actual signed expansion and wall/bulk contrast select
/// the buoyant direction. Opposing flow cannot borrow this aiding correlation.
/// The prescribed flow is an external input, not a finite pump-energy owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AidingVerticalPlateSettings {
    /// Signed vertical free-stream velocity [m/s], positive upward.
    pub speed_up_m_s: f64,
    /// Nonnegative gravity magnitude [m/s2]. Zero removes natural forcing.
    pub gravity_m_s2: f64,
    /// Caller-admitted upper bound on |beta * (T_wall - T_bulk)|, in [0, 1).
    /// This is an approximation policy, not a source-published accuracy bound.
    pub maximum_relative_density_departure: f64,
}

fn invalid(what: &str) -> MaterialConvectionError {
    MaterialError::Parameters { what: what.into() }.into()
}

impl AidingVerticalPlateSettings {
    pub(super) fn validate(self) -> Result<(), MaterialConvectionError> {
        if !self.speed_up_m_s.is_finite()
            || !self.gravity_m_s2.is_finite()
            || self.gravity_m_s2 < 0.0
            || !self.maximum_relative_density_departure.is_finite()
            || !(0.0..1.0).contains(&self.maximum_relative_density_departure)
        {
            return Err(invalid(
                "aiding mixed convection needs finite signed speed, nonnegative gravity and a density-departure bound in [0,1)",
            ));
        }
        Ok(())
    }
}

pub(super) fn inputs(
    transport: FluidTransportState,
    settings: AidingVerticalPlateSettings,
    expansion_per_k: f64,
    wall: Temperature,
    reference: Temperature,
    length: Length,
) -> Result<CorrelationInputs, MaterialConvectionError> {
    settings.validate()?;
    let rayleigh = natural_rayleigh(
        transport,
        settings.gravity_m_s2,
        expansion_per_k,
        wall,
        reference,
        length,
    )?;
    let departure = expansion_per_k * (wall.value() - reference.value());
    if !departure.is_finite() || departure.abs() > settings.maximum_relative_density_departure {
        return Err(invalid(
            "mixed-convection density departure exceeds its declared Boussinesq bound",
        ));
    }
    if rayleigh != 0.0
        && settings.speed_up_m_s != 0.0
        && rayleigh.is_sign_positive() != settings.speed_up_m_s.is_sign_positive()
    {
        return Err(invalid(
            "opposing vertical flow requires a different admitted law or resolved momentum/heat solve",
        ));
    }
    let reynolds = transport.reynolds(settings.speed_up_m_s.abs(), length.value())?;
    Ok(CorrelationInputs::mixed(
        reynolds,
        rayleigh.abs(),
        transport.prandtl(),
    ))
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
