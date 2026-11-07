//! Authored impermeable walls shared by planar and volumetric molecular flow.
use super::WallThermalCondition;
use crate::compressible::stress::Slip;

/// Molecular behavior of one flat physical reflecting side.
///
/// `D` is the dimension of the original Cartesian velocity field. The wall
/// velocity is spatially constant on this side. Its normal component must be
/// zero: translating geometry requires a separate conservative moving-volume
/// owner. `SlipAdiabatic` imposes zero tangential traction, including the
/// normal-velocity cross derivative in its natural velocity trace.
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
    /// Finite Navier slip at a massless interface with authored friction heat.
    ///
    /// The complete tangential traction on the fluid satisfies
    /// `t = (mu/b) (u_wall-u_f)` for positive slip length `b`. Its symmetric
    /// Newtonian strain includes the tangential derivative of normal velocity.
    /// The compatible cell closure solves this trace together with adjacent
    /// wall traces; the elementary `h+b` series formula applies when that
    /// cross derivative vanishes. The interface dissipates
    /// `D = -t·(u_f-u_wall) >= 0` [W/m²]. A fraction `chi` returns
    /// to the fluid as heat; the complementary heat goes to the wall. Thus
    /// the inward fluid energy transfer is `t·u_wall - (1-chi)D + q_heat`.
    /// This equals actual fluid traction work `t·u_f + chi D + q_heat`.
    /// No separate volume heat source is added.
    ///
    /// Frozen-composition and evolving-species planar and spatial molecular
    /// owners support this condition, including their coupled and
    /// reactive steps. It is a continuum boundary law, not a rarefied-gas
    /// model: no mean free path, accommodation or temperature jump is inferred.
    NavierSlip {
        /// Wall velocity in the original global Cartesian frame [m/s].
        velocity: [f64; D],
        /// Adiabatic or prescribed-temperature conductive heat transfer.
        thermal: WallThermalCondition,
        /// Finite nonnegative slip length [m]; zero is exactly `NoSlip`.
        slip_length_m: f64,
        /// Authored fraction of interfacial friction heat deposited in the
        /// fluid, in `[0,1]`. `1` represents a thermally insulating wall that
        /// receives none of this heat; `0` delivers it all to the wall. The
        /// complementary wall heat is included in the boundary energy flux,
        /// but this prescribed-wall owner does not evolve wall temperature.
        fluid_heat_fraction: f64,
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
    /// The zero-length Navier limit is the same prescribed velocity trace as
    /// no slip, independently of its now-inactive friction heat partition.
    pub(in crate::compressible) fn no_slip_velocity(self) -> Option<[f64; D]> {
        match self {
            Self::NoSlip { velocity, .. }
            | Self::NavierSlip {
                velocity,
                slip_length_m: 0.0,
                ..
            } => Some(velocity),
            _ => None,
        }
    }

    /// Reflecting wall traces admitted by the cell-weighted mechanical energy
    /// identity. Construction independently validates the authored slip data.
    pub(in crate::compressible) fn supports_compatible_stress(self) -> bool {
        match self {
            Self::SlipAdiabatic | Self::NoSlip { .. } => true,
            Self::NavierSlip { slip_length_m, .. } => {
                slip_length_m.is_finite() && slip_length_m >= 0.0
            }
        }
    }

    /// Prescribed velocity, tangential slip law, and fluid friction-heat share
    /// for the shared compatible gradient and traction kernels. Normalize the
    /// zero-slip case before selecting its now-inactive heat fraction.
    pub(in crate::compressible) fn stress_boundary(self) -> ([f64; D], Slip, f64) {
        if let Some(velocity) = self.no_slip_velocity() {
            return (velocity, Slip::Length(0.0), 1.0);
        }
        match self {
            Self::SlipAdiabatic => ([0.0; D], Slip::Free, 1.0),
            Self::NavierSlip {
                velocity,
                slip_length_m,
                fluid_heat_fraction,
                ..
            } => (velocity, Slip::Length(slip_length_m), fluid_heat_fraction),
            Self::NoSlip { .. } => unreachable!(),
        }
    }

    /// Geometry owners independently require a physical reflecting side.
    pub(in crate::compressible) fn validate(
        self,
        normal_axis: usize,
        minimum_temperature: f64,
        maximum_temperature: f64,
    ) -> Result<(), &'static str> {
        if let Self::NavierSlip {
            slip_length_m,
            fluid_heat_fraction,
            ..
        } = self
            && (!slip_length_m.is_finite()
                || slip_length_m < 0.0
                || !fluid_heat_fraction.is_finite()
                || !(0.0..=1.0).contains(&fluid_heat_fraction))
        {
            return Err(
                "Navier slip requires a finite nonnegative length and heat fraction in [0,1]",
            );
        }
        if let Self::NoSlip { velocity, thermal }
        | Self::NavierSlip {
            velocity, thermal, ..
        } = self
        {
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
            }
            | Self::NavierSlip {
                thermal: WallThermalCondition::Isothermal(t),
                ..
            } => 0.5 * cell_temperature + 0.5 * t.value(),
            _ => cell_temperature,
        }
    }

    /// Fraction of a cell's tangential derivative retained at the physical
    /// face. The prescribed wall velocity and slip length are constant on a
    /// side, but a finite-slip fluid trace need not be tangentially uniform.
    pub(in crate::compressible) fn tangential_trace_fraction(self, wall_distance: f64) -> f64 {
        match self {
            Self::SlipAdiabatic => 1.0,
            Self::NoSlip { .. } => 0.0,
            Self::NavierSlip { slip_length_m, .. } => {
                if slip_length_m <= wall_distance {
                    let ratio = slip_length_m / wall_distance;
                    ratio / (1.0 + ratio)
                } else {
                    1.0 / (1.0 + wall_distance / slip_length_m)
                }
            }
        }
    }

    /// Native fallback mirror about the fluid trace used by its wall traction.
    /// Normal velocity always reflects about zero. Zero slip uses the original
    /// no-slip arithmetic; free-slip sides retain the Euler reflection ghost.
    /// Compatible owners instead solve the full symmetric-strain slip trace
    /// through `stress::velocity_gradient`, including cross derivatives.
    pub(in crate::compressible) fn velocity_ghost(
        self,
        normal_axis: usize,
        wall_distance: f64,
        cell_velocity: [f64; D],
    ) -> Option<[f64; D]> {
        match self {
            Self::SlipAdiabatic => None,
            Self::NoSlip { velocity, .. } => Some(core::array::from_fn(|axis| {
                2.0 * velocity[axis] - cell_velocity[axis]
            })),
            Self::NavierSlip {
                velocity,
                slip_length_m,
                ..
            } => {
                // Form the wall weight directly, retaining small wall
                // influence when the slip length is much larger than h.
                let wall_weight = if slip_length_m <= wall_distance {
                    1.0 / (1.0 + slip_length_m / wall_distance)
                } else {
                    let ratio = wall_distance / slip_length_m;
                    ratio / (1.0 + ratio)
                };
                Some(core::array::from_fn(|axis| {
                    if axis == normal_axis || slip_length_m == 0.0 {
                        2.0 * velocity[axis] - cell_velocity[axis]
                    } else {
                        // Relative form preserves exact co-moving translation.
                        cell_velocity[axis]
                            + 2.0 * wall_weight * (velocity[axis] - cell_velocity[axis])
                    }
                }))
            }
        }
    }

    /// One-sided resolved Newtonian shear and complete fluid energy flux.
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
        let (velocity, slip_length, fluid_heat_fraction) = match self {
            Self::SlipAdiabatic => return ([0.0; D], 0.0),
            Self::NoSlip { velocity, .. } => (velocity, 0.0, 1.0),
            Self::NavierSlip {
                velocity,
                slip_length_m,
                fluid_heat_fraction,
                ..
            } => (velocity, slip_length_m, fluid_heat_fraction),
        };
        let shear = core::array::from_fn(|axis| {
            if axis == normal_axis {
                0.0
            } else if slip_length == 0.0 {
                viscosity * outward_sign * (velocity[axis] - cell_velocity[axis]) / wall_distance
            } else {
                // Scale the series length before summing: h+b need not fit
                // even when both authored lengths and the traction do.
                let scale = wall_distance.max(slip_length);
                viscosity * outward_sign * ((velocity[axis] - cell_velocity[axis]) / scale)
                    / (wall_distance / scale + slip_length / scale)
            }
        });
        let mut work = 0.0;
        for axis in 0..D {
            if axis != normal_axis {
                work -= shear[axis] * velocity[axis];
            }
        }
        if slip_length > 0.0 && fluid_heat_fraction < 1.0 {
            let fraction = self.tangential_trace_fraction(wall_distance);
            let mut dissipation = 0.0;
            for axis in 0..D {
                if axis != normal_axis {
                    let relative_trace = fraction * (cell_velocity[axis] - velocity[axis]);
                    dissipation -= outward_sign * shear[axis] * relative_trace;
                }
            }
            // -tau·u_wall is wall mechanical work with the flux sign. The
            // complement below is actual slip heat delivered to the wall.
            // Keeping all heat in the fluid needs no extra energy source:
            // its dissipated kinetic energy is already internal energy.
            work += outward_sign * (1.0 - fluid_heat_fraction) * dissipation;
        }
        let heat = self.conductive_energy_flux(
            outward_sign,
            cell_temperature,
            wall_distance,
            conductivity,
        );
        (shear, work + heat)
    }

    /// Original Fourier heat flux before any mechanical work. Keeping this
    /// term separate avoids recovering finite heat by subtracting large work.
    pub(in crate::compressible) fn conductive_energy_flux(
        self,
        outward_sign: f64,
        cell_temperature: f64,
        wall_distance: f64,
        conductivity: f64,
    ) -> f64 {
        match self {
            Self::NoSlip {
                thermal: WallThermalCondition::Isothermal(t),
                ..
            }
            | Self::NavierSlip {
                thermal: WallThermalCondition::Isothermal(t),
                ..
            } => -conductivity * outward_sign * (t.value() - cell_temperature) / wall_distance,
            _ => 0.0,
        }
    }
}
