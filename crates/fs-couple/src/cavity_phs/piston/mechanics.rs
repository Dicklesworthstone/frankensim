//! The original midpoint wall kinematics and pressure-impulse account.
//!
//! Gas and equilibrium-water owners supply their own reciprocal work pressure.
//! This shared kernel owns no fluid state, equation of state or heat closure.

/// Mechanical boundary retained by the original inertial-wall owners.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PistonSupport {
    /// The original finite-mass translational momentum equation.
    Free,
    /// An ideal support holds the accepted wall position exactly.
    Fixed,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct MechanicalTrial {
    pub(crate) displacement: f64,
    pub(crate) velocity: f64,
    pub(crate) residual: f64,
    pub(crate) support_impulse: f64,
}

/// Retain the original gas-piston arithmetic and admission order verbatim.
/// The storage owner validates its initial mass, geometry and physical clock.
#[allow(clippy::too_many_arguments)]
pub(crate) fn trial(
    mass_kg: f64,
    area_m2: f64,
    initial_volume_m3: f64,
    initial_velocity_m_s: f64,
    support: PistonSupport,
    volume_m3: f64,
    pressure_pa: f64,
    duration_s: f64,
    exterior_pressure_pa: f64,
    external_force_n: f64,
) -> Result<MechanicalTrial, &'static str> {
    let displacement = (volume_m3 - initial_volume_m3) / area_m2;
    let velocity = 2.0 * (displacement / duration_s) - initial_velocity_m_s;
    let force = area_m2 * (pressure_pa - exterior_pressure_pa) + external_force_n;
    let support_impulse = if support == PistonSupport::Fixed {
        if initial_velocity_m_s != 0.0 || volume_m3 != initial_volume_m3 {
            return Err("fixed support requires its original volume and zero velocity");
        }
        -duration_s * force
    } else {
        0.0
    };
    let residual =
        mass_kg * (velocity - initial_velocity_m_s) - duration_s * force - support_impulse;
    if !positive(pressure_pa)
        || !displacement.is_finite()
        || !velocity.is_finite()
        || !force.is_finite()
        || !support_impulse.is_finite()
        || !residual.is_finite()
        || !positive(volume_m3 / area_m2)
        || !(0.5 * mass_kg * velocity * velocity).is_finite()
    {
        return Err("piston candidate position, velocity or impulse is unrepresentable");
    }
    Ok(MechanicalTrial {
        displacement,
        velocity,
        residual,
        support_impulse,
    })
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
