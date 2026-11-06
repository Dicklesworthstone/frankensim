//! One exterior gas for radiation, with the bore gas retained in every chimney.
//!
//! Pressure and volume flow are continuous at each chimney mouth. Its finite
//! transmission line therefore receives the exterior's radiation impedance in
//! Pa s/m³ directly, without a density ratio or a second end correction. This
//! is frozen, pressure-matched acoustics: diffusion and mean flow are absent.

use crate::{
    Duct, DuctError, DuctResponse, LoadSpec, LossModel, OpenHoleFlow, Termination, WallModel,
    compose,
};
use fs_material::gas::GasState;
use fs_math::c64::C64;

/// Most open holes whose orifice terminals [`input_impedance_holes_and_ports`] exposes at once.
pub use crate::compose::MAX_ORIFICE_PORTS;

pub(crate) fn admit(interior: &GasState, exterior: &GasState) -> Result<(), DuctError> {
    crate::validate_gas(interior)?;
    crate::validate_gas(exterior)?;
    if interior.pressure != exterior.pressure {
        return Err(DuctError::BadParameter {
            what: "duct and exterior need a common static pressure; mean flow is not modeled",
        });
    }
    Ok(())
}

/// Radiation into the exterior with the local duct's higher-mode impedances
/// retained in the flanged Schur complement. Other termination models retain
/// their existing scope; closed/ideal ports acquire no radiating load.
pub(crate) fn termination_impedance(
    termination: Termination,
    interior: &GasState,
    exterior: &GasState,
    radius: f64,
    omega: f64,
) -> Result<(Option<C64>, f64), DuctError> {
    admit(interior, exterior)?;
    if interior != exterior && termination == Termination::FlangedOpen {
        return Ok((Some(crate::flanged_pipe_impedance_with_exterior(
            interior, exterior, radius, omega, crate::FLANGED_MODES,
        )?), omega / exterior.sound_speed * radius));
    }
    crate::termination_impedance(termination, exterior, radius, omega)
}

fn load(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    termination: Termination,
) -> Result<(LoadSpec, Option<f64>), DuctError> {
    admit(interior, exterior)?;
    if interior == exterior || matches!(termination, Termination::Closed | Termination::IdealOpen) {
        return Ok((LoadSpec::Analytic(termination), None));
    }
    let radius = duct
        .segments
        .last()
        .ok_or(DuctError::EmptyDuct)?
        .outlet_radius();
    let (impedance, _) = termination_impedance(termination, interior, exterior, radius, omega)?;
    Ok((
        LoadSpec::Fixed(impedance.expect("radiating termination")),
        Some(omega * radius / exterior.sound_speed),
    ))
}

/// Input response with independent exterior radiation at the bell and every
/// open/vented side hole. Closed chimney compliance, viscothermal propagation,
/// internal junction corrections and wall response keep the interior gas.
/// Identical media preserve the original arithmetic and validity checks.
///
/// # Errors
/// Incoherent media, unequal static pressure, or existing duct/wall/load refusals.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: WallModel<'_>,
) -> Result<DuctResponse, DuctError> {
    let (load, ka) = load(duct, interior, exterior, omega, termination)?;
    let mut response = compose::evaluate_in_exterior(
        duct, interior, exterior, omega, loss, load, wall,
    )?;
    if let Some(ka) = ka {
        response.mouth_ka = ka;
    }
    Ok(response)
}

/// The same solution with each side hole's outward volume-flow transfer.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_and_holes(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: WallModel<'_>,
) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
    let (load, ka) = load(duct, interior, exterior, omega, termination)?;
    let (mut response, holes) = compose::evaluate_in_exterior_and_holes(
        duct, interior, exterior, omega, loss, load, wall,
    )?;
    if let Some(ka) = ka {
        response.mouth_ka = ka;
    }
    Ok((response, holes))
}

/// [`input_impedance_and_holes`] plus the admittance of the linear network at the orifice terminals of
/// the first `ports` open holes (at most [`compose::MAX_ORIFICE_PORTS`]): row-major `m * m` with
/// `U_i = - sum_j Y_ij dp_j`; see `compose::open_hole_port_admittance`. A coupled orifice needs its
/// own network, so the cross-hole terms the one-port reflections of [`OpenHoleFlow`] do not carry live here.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_holes_and_ports(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: WallModel<'_>,
    ports: usize,
) -> Result<(DuctResponse, Vec<OpenHoleFlow>, Vec<C64>), DuctError> {
    let (load, ka) = load(duct, interior, exterior, omega, termination)?;
    let (mut response, holes, admittance) = compose::evaluate_in_exterior_holes_and_ports(
        duct, interior, exterior, omega, loss, load, wall, ports,
    )?;
    if let Some(ka) = ka {
        response.mouth_ka = ka;
    }
    Ok((response, holes, admittance))
}

/// Direct inlet reflection of the same two-medium chain. The characteristic
/// feedback remains representable at a lossless input-impedance pole.
#[allow(clippy::too_many_arguments)]
pub fn input_reflection(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    termination: Termination,
    wall: WallModel<'_>,
) -> Result<compose::InputReflection, DuctError> {
    admit(interior, exterior)?;
    let options = compose::ComposeOptions::default();
    let (chain, _) = compose::duct_chain_in_exterior(
        duct,
        interior,
        exterior,
        omega,
        loss,
        wall,
        options,
        true,
        || false,
    )?;
    let radius = duct
        .segments
        .last()
        .ok_or(DuctError::EmptyDuct)?
        .outlet_radius();
    let (load, _) = termination_impedance(termination, interior, exterior, radius, omega)?;
    compose::compose_chain_reflection(&chain.blocks, load, options)
}

/// Explicit terminal load with the same independently supplied exterior at
/// side holes. The terminal is a pressure/volume-flow port; it may connect to
/// another gas region or a source-bound baked radiation load.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_load(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    z_load: C64,
    wall: WallModel<'_>,
) -> Result<DuctResponse, DuctError> {
    admit(interior, exterior)?;
    if !z_load.re.is_finite() || !z_load.im.is_finite() {
        return Err(DuctError::BadParameter {
            what: "tabulated load must be finite",
        });
    }
    compose::evaluate_in_exterior(
        duct,
        interior,
        exterior,
        omega,
        loss,
        LoadSpec::Fixed(z_load),
        wall,
    )
}

/// Explicit terminal load and exterior, including outward side-hole flows.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_load_and_holes(
    duct: &Duct,
    interior: &GasState,
    exterior: &GasState,
    omega: f64,
    loss: LossModel,
    z_load: C64,
    wall: WallModel<'_>,
) -> Result<(DuctResponse, Vec<OpenHoleFlow>), DuctError> {
    admit(interior, exterior)?;
    if !z_load.re.is_finite() || !z_load.im.is_finite() {
        return Err(DuctError::BadParameter {
            what: "tabulated load must be finite",
        });
    }
    compose::evaluate_in_exterior_and_holes(
        duct,
        interior,
        exterior,
        omega,
        loss,
        LoadSpec::Fixed(z_load),
        wall,
    )
}

#[cfg(test)]
mod tests;
