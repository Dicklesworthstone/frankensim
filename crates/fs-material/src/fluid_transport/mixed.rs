//! Dimensional admission for laminar aiding flow along a vertical plate.
//!
//! This owns the material-state and orientation checks, not a heat-transfer
//! correlation. Consumers pass the admitted groups to the original correlation
//! owner and retain that owner's domain, evidence and model discrepancy.

use super::{FluidTransportState, invalid};
use crate::MaterialError;
use fs_qty::{Length, Temperature};

/// Prescribed vertical free-stream motion and an explicit Boussinesq bound.
///
/// Positive speed is opposite gravity; negative speed is along gravity. The
/// actual signed expansion and wall/bulk contrast select the buoyant direction.
/// Opposing flow cannot borrow an aiding-flow correlation. Prescribed speed is
/// an external input, not a finite pump-energy owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AidingVerticalPlateSettings {
    /// Signed vertical free-stream velocity [m/s], positive upward.
    pub speed_up_m_s: f64,
    /// Nonnegative gravity magnitude [m/s²]. Zero removes natural forcing.
    pub gravity_m_s2: f64,
    /// Caller-admitted bound on `|beta * (T_wall - T_bulk)|`, in `[0, 1)`.
    /// This is an approximation policy, not a source-published accuracy bound.
    pub maximum_relative_density_departure: f64,
}

/// Groups from one admitted material snapshot and the full wall/bulk contrast.
///
/// Both drivers may be zero. This does not supply a stagnant-gas conductance:
/// a thermal owner may then certify zero heat only at exact equilibrium, or
/// defer to a correlation that explicitly admits the requested resting state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AidingVerticalPlateGroups {
    /// Reynolds number from axial speed magnitude and plate height.
    pub reynolds: f64,
    /// Nonnegative magnitude after signed buoyancy/flow alignment is checked.
    pub rayleigh: f64,
    /// Prandtl number of the same retained transport state.
    pub prandtl: f64,
}

impl AidingVerticalPlateSettings {
    /// Admit explicit finite forcing and a caller-selected density bound.
    pub fn validate(self) -> Result<(), MaterialError> {
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

    /// Form Re, |Ra| and Pr without changing the original property snapshot.
    ///
    /// The same snapshot supplies all groups. `wall` and `reference` are the
    /// actual wall/bulk temperatures, even when the snapshot is at mean film
    /// temperature. Signed expansion is explicit; no gas-only `1/T` law is
    /// substituted for a liquid or another source's isobaric derivative.
    ///
    /// This admits alignment and the declared Boussinesq density departure;
    /// it does not establish constant-property accuracy or a correlation's
    /// laminar/Prandtl domain. Those obligations remain with their owners.
    pub fn groups(
        self,
        transport: FluidTransportState,
        expansion_per_k: f64,
        wall: Temperature,
        reference: Temperature,
        length: Length,
    ) -> Result<AidingVerticalPlateGroups, MaterialError> {
        self.validate()?;
        if [wall.value(), reference.value()]
            .into_iter().any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(invalid("mixed convection needs positive finite absolute temperatures"));
        }
        let difference = wall.value() - reference.value();
        let rayleigh = transport.rayleigh(
            self.gravity_m_s2, expansion_per_k, difference, length.value(),
        )?;
        let departure = expansion_per_k * difference;
        if !departure.is_finite() || departure.abs() > self.maximum_relative_density_departure {
            return Err(invalid(
                "mixed-convection density departure exceeds its declared Boussinesq bound",
            ));
        }
        if rayleigh != 0.0
            && self.speed_up_m_s != 0.0
            && rayleigh.is_sign_positive() != self.speed_up_m_s.is_sign_positive()
        {
            return Err(invalid(
                "opposing vertical flow requires a different admitted law or resolved momentum/heat solve",
            ));
        }
        Ok(AidingVerticalPlateGroups {
            reynolds: transport.reynolds(self.speed_up_m_s.abs(), length.value())?,
            rayleigh: rayleigh.abs(),
            prandtl: transport.prandtl(),
        })
    }
}

/// Authored vertical plate axis in the flow/gravity world frame.
///
/// Nonzero gravity must be opposite `plate_up`; moving flow must be parallel
/// to that axis. Neither inclined gravity nor crossflow is projected into an
/// admitted vertical case. The signed material response is checked by
/// [`AidingVerticalPlateSettings::groups`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AidingVerticalPlateFrame {
    plate_up: [f64; 3],
    gravity_m_s2: f64,
}

// Only normalization roundoff is tolerated, not an empirical angular range.
const ALIGNMENT_ROUNDOFF: f64 = 64.0 * f64::EPSILON;

fn direction(vector: [f64; 3]) -> Result<([f64; 3], f64), MaterialError> {
    if vector.iter().any(|value| !value.is_finite()) {
        return Err(invalid("mixed-convection frame vectors must be finite"));
    }
    let scale = vector.iter().fold(0.0_f64, |a, value| a.max(value.abs()));
    if scale == 0.0 {
        return Err(invalid("mixed-convection direction must be nonzero"));
    }
    let scaled = vector.map(|value| value / scale);
    let length = fs_math::det::sqrt(scaled.iter().map(|value| value * value).sum());
    Ok((scaled.map(|value| value / length), scale * length))
}

fn aligned(actual: [f64; 3], expected: [f64; 3]) -> bool {
    actual.into_iter().zip(expected)
        .all(|(a, b)| (a - b).abs() <= ALIGNMENT_ROUNDOFF)
}

impl AidingVerticalPlateFrame {
    /// Normalize the plate's upward axis and admit the actual gravity vector.
    pub fn try_new(plate_up: [f64; 3], gravity_m_s2: [f64; 3]) -> Result<Self, MaterialError> {
        let (plate_up, _) = direction(plate_up)?;
        let gravity_m_s2 = if gravity_m_s2 == [0.0; 3] {
            0.0
        } else {
            let (down, magnitude) = direction(gravity_m_s2)?;
            if !magnitude.is_finite() || !aligned(down, plate_up.map(|value| -value)) {
                return Err(invalid(
                    "aiding vertical-plate gravity must be representable and opposite the authored plate-up axis",
                ));
            }
            magnitude
        };
        Ok(Self { plate_up, gravity_m_s2 })
    }

    /// Admit actual relative flow, preserving its signed axial speed.
    pub fn settings(
        self,
        relative_velocity_m_s: [f64; 3],
        maximum_relative_density_departure: f64,
    ) -> Result<AidingVerticalPlateSettings, MaterialError> {
        let speed_up_m_s = if relative_velocity_m_s == [0.0; 3] {
            0.0
        } else {
            let (flow, magnitude) = direction(relative_velocity_m_s)?;
            let projection: f64 = flow.into_iter().zip(self.plate_up)
                .map(|(a, b)| a * b).sum();
            let sign = if projection < 0.0 { -1.0 } else { 1.0 };
            if !magnitude.is_finite() || !aligned(flow, self.plate_up.map(|value| sign * value)) {
                return Err(invalid(
                    "aiding vertical-plate relative flow must be representable and parallel to the authored plate axis",
                ));
            }
            sign * magnitude
        };
        let settings = AidingVerticalPlateSettings {
            speed_up_m_s, gravity_m_s2: self.gravity_m_s2,
            maximum_relative_density_departure,
        };
        settings.validate()?;
        Ok(settings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings(speed_up_m_s: f64) -> AidingVerticalPlateSettings {
        AidingVerticalPlateSettings {
            speed_up_m_s, gravity_m_s2: 9.81, maximum_relative_density_departure: 0.01,
        }
    }

    #[test]
    fn g0_shared_groups_preserve_full_contrast_and_signed_expansion() {
        let fluid = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        for (wall, beta, speed) in [(310.0, 2e-4, 0.1), (290.0, 2e-4, -0.1),
            (310.0, -2e-4, -0.1), (290.0, -2e-4, 0.1)]
        {
            let groups = settings(speed).groups(fluid, beta, Temperature::new(wall),
                Temperature::new(300.0), Length::new(0.1)).unwrap();
            assert!((groups.reynolds / 10000.0 - 1.0).abs() < 1e-14);
            assert!((groups.rayleigh / 156_960_000.0 - 1.0).abs() < 1e-14);
            assert_eq!(groups.prandtl, 8.0);
            assert!(settings(-speed).groups(fluid, beta, Temperature::new(wall),
                Temperature::new(300.0), Length::new(0.1)).is_err());
        }
    }

    #[test]
    fn g0_pure_limits_do_not_discard_density_or_absolute_temperature_admission() {
        let fluid = FluidTransportState::try_new(1.0, 1000.0, 0.03, 1.8e-5).unwrap();
        let sample = |s: AidingVerticalPlateSettings, wall| s.groups(fluid, 0.002,
            Temperature::new(wall), Temperature::new(300.0), Length::new(0.1));
        let natural = sample(settings(0.0), 305.0).unwrap();
        assert_eq!(natural.reynolds, 0.0);
        assert!(natural.rayleigh > 0.0);
        let forced = sample(AidingVerticalPlateSettings { gravity_m_s2: 0.0, ..settings(-0.1) }, 305.0).unwrap();
        assert!(forced.reynolds > 0.0);
        assert_eq!(forced.rayleigh, 0.0);
        let resting = sample(settings(0.0), 300.0).unwrap();
        assert_eq!((resting.reynolds, resting.rayleigh), (0.0, 0.0));
        assert!(sample(AidingVerticalPlateSettings {
            maximum_relative_density_departure: 0.01_f64.next_down(), ..settings(0.1)
        }, 305.0).is_err());
        assert!(sample(AidingVerticalPlateSettings { gravity_m_s2: 0.0, ..settings(0.1) }, 310.0).is_err());
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(sample(settings(0.1), value).is_err());
            assert!(settings(0.1).groups(fluid, 0.002, Temperature::new(300.0),
                Temperature::new(value), Length::new(0.1)).is_err());
        }
        for bound in [-1.0, 1.0, f64::NAN, f64::INFINITY] {
            assert!(AidingVerticalPlateSettings {
                maximum_relative_density_departure: bound, ..settings(0.1)
            }.validate().is_err());
        }
    }

    #[test]
    fn g0_frame_refuses_oblique_flow_instead_of_projecting_it() {
        let frame = AidingVerticalPlateFrame::try_new([0.0, 0.0, 3.0], [0.0, 0.0, -9.81]).unwrap();
        assert_eq!(frame.settings([0.0, 0.0, -0.1], 0.01).unwrap().speed_up_m_s, -0.1);
        assert_eq!(frame.settings([0.0; 3], 0.01).unwrap().speed_up_m_s, 0.0);
        assert!(frame.settings([0.01, 0.0, 0.1], 0.01).is_err());
        assert!(frame.settings([f64::NAN, 0.0, 0.1], 0.01).is_err());
        assert!(AidingVerticalPlateFrame::try_new([0.0; 3], [0.0; 3]).is_err());
        assert!(AidingVerticalPlateFrame::try_new([0.0, 0.0, 1.0], [0.01, 0.0, -9.81]).is_err());
        let zero_g = AidingVerticalPlateFrame::try_new([1.0, 0.0, 0.0], [0.0; 3]).unwrap();
        assert_eq!(zero_g.settings([-0.1, 0.0, 0.0], 0.0).unwrap().gravity_m_s2, 0.0);
    }
}
