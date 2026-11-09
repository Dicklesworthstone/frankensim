//! Live constitutive resistance with one independent finite support.

use super::super::recoil::{complete_motion, prepare_contact};
use super::super::{MaterialWallCaloric, MixtureCaloric, RecoilControl, RecoilUpdate};
use super::{
    Cell, MixtureGas, Path, Primitive, SlidingSupport, SourceCoupledDragControl,
    SourceCoupledDragReport, Time, WallContact, WallDrag, WallExchangeError, WallMotion, advance,
    nonnegative, poll, positive, work_error,
};

/// Complete original gas/support/thermal endpoint with source-integration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SourceCoupledRecoilUpdate {
    /// Original recoil, caloric storage, physical slip partition and balances.
    pub recoil: RecoilUpdate,
    /// Estimates on relative slip and total physical loss before its partition.
    pub source: SourceCoupledDragReport,
}

/// Integrate current-state friction against one independent finite support.
///
/// The original recoil owner gives every trial gas impulse its opposite support
/// impulse and actual mechanical work. Each source selection sees the original
/// gas primitive and current support velocity. Its returned drag must retain that
/// velocity. Density, volume, support mass and the authored `gas_heat_fraction`
/// stay fixed; a shared support requires a collective momentum solve.
///
/// During mechanics the gas source receives its actual share of slip heating;
/// wall heat is deposited once through the same original finite-storage owner
/// as frozen recoil. No tangent substitutes for a material energy integral. The
/// source controller measures total relative-slip loss before partition, so a
/// zero gas heat fraction still controls the heat delivered to the finite wall.
///
/// After mechanics, this entry performs the existing first-order thermal step
/// using the authored thermal contact, then admits the actual returned gas and
/// support through the original resolver again. Network callers set thermal
/// conductance to zero here before applying their separately refreshed thermal
/// and conduction owners. Sources refusing a predictor or final thermal endpoint,
/// work exhaustion and cancellation return neither accepted state.
#[allow(clippy::too_many_arguments)]
pub fn exchange_cell_with_source_recoil(
    gas: &MixtureGas,
    state: Cell<'_>,
    volume_m3: f64,
    duration: Time,
    support: SlidingSupport,
    gas_heat_fraction: f64,
    contact: WallContact,
    wall: Option<&MaterialWallCaloric<'_>>,
    control: RecoilControl,
    source_control: SourceCoupledDragControl,
    mut resolve: impl FnMut(
        Cell<'_>,
        SlidingSupport,
        Primitive,
        &mut dyn FnMut() -> bool,
    ) -> Result<WallDrag, WallExchangeError>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SourceCoupledRecoilUpdate, WallExchangeError> {
    poll(&mut cancelled)?;
    contact.validate()?;
    source_control.validate()?;
    if !positive(volume_m3) || !nonnegative(duration.value()) {
        return Err(WallExchangeError::Invalid(
            "source-coupled recoil volume or interval",
        ));
    }
    let initial_contact = prepare_contact(
        contact.thermal,
        wall,
        gas_heat_fraction,
        control,
        &mut cancelled,
    )?;
    let primitive = gas.primitive(state, &mut cancelled)?;
    let path = Path {
        gas,
        state,
        density: primitive.density,
        wall: WallMotion::Recoil {
            support,
            gas_heat_fraction,
        },
        volume: volume_m3,
    };
    let (stage, mut report) = advance(
        &path,
        duration.value(),
        source_control,
        &mut |cell, primitive, velocity_m_s, cancelled| {
            resolve(
                cell,
                SlidingSupport {
                    velocity_m_s,
                    ..support
                },
                primitive,
                cancelled,
            )
        },
        &mut cancelled,
    )?;
    let motion = path.motion(stage.momentum)?;
    let mut recoil = complete_motion(
        &MixtureCaloric {
            gas,
            partial_densities: state.partial_densities,
        },
        state.energy,
        volume_m3,
        duration,
        support,
        gas_heat_fraction,
        initial_contact,
        wall,
        control,
        motion,
        &mut cancelled,
    )?;
    // Unlike the mechanical-only prescribed entry, original thermal completion
    // can change gas temperature here. Charge and admit that actual endpoint.
    poll(&mut cancelled)?;
    if report.source_evaluations >= source_control.maximum_source_evaluations {
        return Err(work_error(report));
    }
    report.source_evaluations += 1;
    let actual = Cell {
        partial_densities: state.partial_densities,
        momentum: recoil.gas.momentum,
        energy: recoil.gas.energy,
    };
    let actual_primitive = gas.primitive(actual, &mut cancelled)?;
    let drag = resolve(actual, recoil.support, actual_primitive, &mut cancelled)?;
    drag.validate()?;
    if drag.moving_velocity().unwrap_or(0.0) != recoil.support.velocity_m_s {
        return Err(WallExchangeError::Invalid(
            "source changed recoil endpoint velocity",
        ));
    }
    recoil.gas.contact.drag = drag;
    poll(&mut cancelled)?;
    Ok(SourceCoupledRecoilUpdate {
        recoil,
        source: report,
    })
}

#[cfg(test)]
mod tests;
