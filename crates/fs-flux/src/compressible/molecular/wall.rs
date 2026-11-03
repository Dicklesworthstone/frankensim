//! Authored impermeable walls shared by planar and volumetric molecular flow.
use super::WallThermalCondition;

/// Molecular behavior of one flat physical reflecting side.
///
/// `D` is the dimension of the original Cartesian velocity field. The wall
/// velocity is spatially constant on this side. Its normal component must be
/// zero: translating geometry requires a separate conservative moving-volume
/// owner. Slip walls retain their original stationary adiabatic behavior.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MolecularWall<const D: usize> {
    /// Stationary, adiabatic, impermeable slip wall.
    SlipAdiabatic,
    /// No-slip wall with prescribed tangential velocity and thermal condition.
    NoSlip {
        /// Wall velocity in the original global Cartesian frame [m/s].
        velocity: [f64; D],
        /// Adiabatic or prescribed-temperature thermal boundary.
        thermal: WallThermalCondition,
    },
}

impl<const D: usize> Default for MolecularWall<D> {
    fn default() -> Self {
        Self::SlipAdiabatic
    }
}

/// Authored molecular wall for the original two-dimensional Cartesian field.
pub type MolecularWall2d = MolecularWall<2>;

impl<const D: usize> MolecularWall<D> {
    /// Geometry owners independently require a physical reflecting side.
    pub(in crate::compressible) fn validate(
        self,
        normal_axis: usize,
        minimum_temperature: f64,
        maximum_temperature: f64,
    ) -> Result<(), &'static str> {
        if let Self::NoSlip { velocity, thermal } = self {
            if velocity.iter().any(|v| !v.is_finite())
                || velocity.get(normal_axis).is_none_or(|v| *v != 0.0)
            {
                return Err("wall velocity must be finite and tangential");
            }
            if let WallThermalCondition::Isothermal(t) = thermal {
                if !t.value().is_finite()
                    || !(minimum_temperature..=maximum_temperature).contains(&t.value())
                {
                    return Err("wall temperature outside caloric domain");
                }
            }
        }
        Ok(())
    }

    /// Original transport law sampling at the cell/wall mean temperature.
    pub(in crate::compressible) fn face_temperature(self, cell_temperature: f64) -> f64 {
        match self {
            Self::NoSlip {
                thermal: WallThermalCondition::Isothermal(t),
                ..
            } => 0.5 * cell_temperature + 0.5 * t.value(),
            _ => cell_temperature,
        }
    }

    /// No-slip velocity ghosts for the neighboring interior cross derivatives.
    /// Slip sides keep the original Euler owner's normal-reflection ghost.
    pub(in crate::compressible) fn velocity_ghost(
        self,
        cell_velocity: [f64; D],
    ) -> Option<[f64; D]> {
        match self {
            Self::SlipAdiabatic => None,
            Self::NoSlip { velocity, .. } => Some(core::array::from_fn(|axis| {
                2.0 * velocity[axis] - cell_velocity[axis]
            })),
        }
    }

    /// One-sided resolved Newtonian shear and its wall work plus Fourier heat.
    /// Coefficients come from the caller's original admitted transport source.
    /// Components are oriented along increasing Cartesian coordinates; the
    /// grid's boundary integration applies the outward sign exactly once.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible) fn shear_and_energy(
        self,
        normal_axis: usize,
        outward_sign: f64,
        cell_velocity: [f64; D],
        cell_temperature: f64,
        wall_distance: f64,
        viscosity: f64,
        conductivity: f64,
    ) -> ([f64; D], f64) {
        let Self::NoSlip { velocity, thermal } = self else {
            return ([0.0; D], 0.0);
        };
        let shear = core::array::from_fn(|axis| {
            if axis == normal_axis {
                0.0
            } else {
                viscosity * outward_sign * (velocity[axis] - cell_velocity[axis]) / wall_distance
            }
        });
        let mut work = 0.0;
        for axis in 0..D {
            if axis != normal_axis {
                work -= shear[axis] * velocity[axis];
            }
        }
        let heat = match thermal {
            WallThermalCondition::Adiabatic => 0.0,
            WallThermalCondition::Isothermal(t) => {
                -conductivity * outward_sign * (t.value() - cell_temperature) / wall_distance
            }
        };
        (shear, work + heat)
    }
}
