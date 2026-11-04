//! Zero-duration momentum accounting for the original translating wall.

use fs_math::dd::Dd;

/// Physical endpoint ceilings and impulse-local numerical admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonImpulseControls {
    /// Maximum absolute accepted wall velocity [m/s].
    pub maximum_abs_velocity_m_s: f64,
    /// Maximum accepted wall kinetic energy [J]. Gas storage is unchanged.
    pub maximum_kinetic_energy_j: f64,
    /// Relative represented-momentum defect, scaled by the impulse, in [0, 1).
    pub momentum_relative_tolerance: f64,
    /// Absolute work-balance allowance [J], independent of gas/kinetic storage.
    pub energy_absolute_tolerance_j: f64,
    /// Relative work-balance allowance on this event's signed work, in [0, 1).
    pub energy_relative_tolerance: f64,
}

impl Default for PistonImpulseControls {
    fn default() -> Self {
        Self {
            maximum_abs_velocity_m_s: 1.0e6,
            maximum_kinetic_energy_j: 1.0e18,
            momentum_relative_tolerance: 1.0e-11,
            energy_absolute_tolerance_j: 1.0e-12,
            energy_relative_tolerance: 1.0e-10,
        }
    }
}

/// Accepted wall momentum and signed work at one fixed configuration.
/// Gas storage, volume, heat and finite-time port records are unchanged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PistonImpulseReceipt {
    /// Original and accepted wall velocities [m/s].
    pub velocity_m_s: [f64; 2],
    /// Requested external impulse, positive toward primary-gas expansion [N s].
    pub impulse_n_s: f64,
    /// Impulse ON the wall from its fixed support [N s]; zero for a free wall.
    pub support_impulse_n_s: f64,
    /// Mass times the represented velocity change [N s].
    pub represented_momentum_change_n_s: f64,
    /// Polarized kinetic-energy increment, including negative braking work [J].
    pub kinetic_energy_change_j: f64,
    /// Requested external impulse times the mean endpoint velocity [J].
    pub external_work_j: f64,
    /// Support impulse times mean velocity [J]; exactly zero for fixed support.
    pub support_work_j: f64,
    /// Represented momentum change minus external and support impulses [N s].
    pub momentum_residual_n_s: f64,
    /// Momentum defect relative to this event's net impulse.
    pub relative_momentum_residual: f64,
    /// Actual kinetic increment minus external and support work [J].
    pub energy_residual_j: f64,
}

/// Refusal of an instantaneous candidate; no accepted owner is mutated.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PistonImpulseError {
    /// Invalid physical data, controls or unrepresentable arithmetic.
    Invalid(&'static str),
    /// A physical endpoint ceiling was exceeded.
    Budget {
        what: &'static str,
        value: f64,
        limit: f64,
    },
    /// The represented wall momentum does not resolve the applied impulse.
    Momentum { relative: f64, tolerance: f64 },
    /// Actual kinetic change and conjugate impulse work disagree.
    Energy { residual_j: f64, tolerance_j: f64 },
    /// Cancellation was observed before publication.
    Cancelled,
}

impl core::fmt::Display for PistonImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "gas-piston impulse refused: {self:?}")
    }
}
impl std::error::Error for PistonImpulseError {}

fn finite(value: Dd) -> Result<Dd, PistonImpulseError> {
    if value.hi.is_finite() && value.lo.is_finite() {
        Ok(value)
    } else {
        Err(PistonImpulseError::Invalid(
            "wall impulse account is unrepresentable",
        ))
    }
}

pub(super) fn candidate(
    mass_kg: f64,
    before: f64,
    fixed: bool,
    impulse_n_s: f64,
    controls: PistonImpulseControls,
) -> Result<PistonImpulseReceipt, PistonImpulseError> {
    if !(mass_kg > 0.0 && mass_kg.is_finite())
        || !before.is_finite()
        || !impulse_n_s.is_finite()
        || (fixed && before != 0.0)
    {
        return Err(PistonImpulseError::Invalid(
            "positive finite wall mass and finite admitted velocity/impulse required",
        ));
    }
    if ![
        controls.maximum_abs_velocity_m_s,
        controls.maximum_kinetic_energy_j,
    ]
    .into_iter()
    .all(|x| x > 0.0 && x.is_finite())
        || ![
            controls.momentum_relative_tolerance,
            controls.energy_relative_tolerance,
        ]
        .into_iter()
        .all(|x| x.is_finite() && (0.0..1.0).contains(&x))
        || !controls.energy_absolute_tolerance_j.is_finite()
        || controls.energy_absolute_tolerance_j < 0.0
    {
        return Err(PistonImpulseError::Invalid(
            "finite positive endpoint ceilings and bounded nonnegative impulse tolerances required",
        ));
    }
    let after = if fixed {
        before
    } else {
        before + impulse_n_s / mass_kg
    };
    if !after.is_finite() {
        return Err(PistonImpulseError::Invalid(
            "wall velocity increment is unrepresentable",
        ));
    }
    if !fixed && impulse_n_s != 0.0 && after == before {
        return Err(PistonImpulseError::Invalid(
            "nonzero wall impulse was lost to velocity rounding",
        ));
    }
    let kinetic_after = (mass_kg * (0.5 * after)) * after;
    if !kinetic_after.is_finite() {
        return Err(PistonImpulseError::Invalid(
            "wall kinetic energy is unrepresentable",
        ));
    }
    for (what, value, limit) in [
        (
            "absolute wall velocity",
            after.abs(),
            controls.maximum_abs_velocity_m_s,
        ),
        (
            "wall kinetic energy",
            kinetic_after,
            controls.maximum_kinetic_energy_j,
        ),
    ] {
        if value > limit {
            return Err(PistonImpulseError::Budget { what, value, limit });
        }
    }
    let support = if fixed { -impulse_n_s } else { 0.0 };
    let requested = Dd::from_f64(impulse_n_s) + Dd::from_f64(support);
    let velocity_change = Dd::from_f64(after) - Dd::from_f64(before);
    let change = finite(Dd::from_f64(mass_kg) * velocity_change)?;
    let residual = finite(change - requested)?;
    let scale = change.to_f64().abs().max(requested.to_f64().abs());
    let relative = if scale == 0.0 {
        0.0
    } else {
        residual.to_f64().abs() / scale
    };
    if !relative.is_finite() || relative > controls.momentum_relative_tolerance {
        return Err(PistonImpulseError::Momentum {
            relative,
            tolerance: controls.momentum_relative_tolerance,
        });
    }
    let mean = Dd::from_f64(0.5 * before) + Dd::from_f64(0.5 * after);
    let kinetic = finite(change * mean)?;
    let external_work = finite(Dd::from_f64(impulse_n_s) * mean)?;
    let support_work = finite(Dd::from_f64(support) * mean)?;
    // Form the defect before multiplication, not by subtracting large works.
    let energy_residual = finite(residual * mean)?.to_f64();
    let tolerance = controls.energy_absolute_tolerance_j
        + controls.energy_relative_tolerance
            * kinetic
                .to_f64()
                .abs()
                .max(finite(external_work + support_work)?.to_f64().abs());
    if !tolerance.is_finite() || energy_residual.abs() > tolerance {
        return Err(PistonImpulseError::Energy {
            residual_j: energy_residual,
            tolerance_j: tolerance,
        });
    }
    Ok(PistonImpulseReceipt {
        velocity_m_s: [before, after],
        impulse_n_s,
        support_impulse_n_s: support,
        represented_momentum_change_n_s: change.to_f64(),
        kinetic_energy_change_j: kinetic.to_f64(),
        external_work_j: external_work.to_f64(),
        support_work_j: support_work.to_f64(),
        momentum_residual_n_s: residual.to_f64(),
        relative_momentum_residual: relative,
        energy_residual_j: energy_residual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn free_wall_kick_brake_and_reversal_close_momentum_and_work() {
        let controls = PistonImpulseControls::default();
        let kick = candidate(4.0, 0.0, false, 12.0, controls).unwrap();
        assert_eq!(kick.velocity_m_s, [0.0, 3.0]);
        assert_eq!(kick.represented_momentum_change_n_s, 12.0);
        assert_eq!(kick.kinetic_energy_change_j, 18.0);
        assert_eq!(kick.external_work_j, 18.0);
        assert_eq!(kick.support_impulse_n_s, 0.0);
        let brake = candidate(4.0, 3.0, false, -12.0, controls).unwrap();
        assert_eq!(brake.velocity_m_s, [3.0, 0.0]);
        assert_eq!(brake.kinetic_energy_change_j, -18.0);
        assert_eq!(brake.external_work_j, -18.0);
        assert_eq!(brake.energy_residual_j, 0.0);
    }

    #[test]
    fn fixed_support_receives_the_opposite_impulse_without_work() {
        for impulse in [-7.0, 0.0, 13.0] {
            let event =
                candidate(2.0, 0.0, true, impulse, PistonImpulseControls::default()).unwrap();
            assert_eq!(event.velocity_m_s, [0.0, 0.0]);
            assert_eq!(event.support_impulse_n_s, -impulse);
            assert_eq!(event.represented_momentum_change_n_s, 0.0);
            assert_eq!(event.kinetic_energy_change_j, 0.0);
            assert_eq!(event.external_work_j, 0.0);
            assert_eq!(event.support_work_j, 0.0);
            assert_eq!(event.momentum_residual_n_s, 0.0);
        }
    }

    #[test]
    fn nonrepresentable_momentum_and_endpoint_limits_refuse() {
        let controls = PistonImpulseControls::default();
        assert!(candidate(1.0, 1.0, false, 1.0e-30, controls).is_err());
        assert!(candidate(1.0, 1.0, true, -1.0, controls).is_err());
        assert!(candidate(0.0, 0.0, false, 1.0, controls).is_err());
        assert!(candidate(1.0, 0.0, false, f64::NAN, controls).is_err());
        let velocity = PistonImpulseControls {
            maximum_abs_velocity_m_s: 0.5,
            ..controls
        };
        assert!(matches!(
            candidate(1.0, 0.0, false, 1.0, velocity),
            Err(PistonImpulseError::Budget {
                what: "absolute wall velocity",
                ..
            })
        ));
        let energy = PistonImpulseControls {
            maximum_kinetic_energy_j: 0.25,
            ..controls
        };
        assert!(matches!(
            candidate(1.0, 0.0, false, 1.0, energy),
            Err(PistonImpulseError::Budget {
                what: "wall kinetic energy",
                ..
            })
        ));
    }

    #[test]
    fn subdivision_preserves_signed_work_at_the_same_configuration() {
        let controls = PistonImpulseControls::default();
        let whole = candidate(2.0, -1.0, false, 3.0, controls).unwrap();
        let first = candidate(2.0, -1.0, false, 1.0, controls).unwrap();
        let second = candidate(2.0, first.velocity_m_s[1], false, 2.0, controls).unwrap();
        assert_eq!(whole.velocity_m_s[1], second.velocity_m_s[1]);
        assert_eq!(
            whole.external_work_j,
            first.external_work_j + second.external_work_j
        );
        assert_eq!(whole.external_work_j, -0.75);
    }
}
