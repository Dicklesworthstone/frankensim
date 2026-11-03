//! SI construction and live reservoir targets on the original scaled flow owner.

use super::*;

#[cfg(test)]
#[path = "boundary/tests.rs"]
mod tests;

/// Prescribed reservoir data in the fixed world/grid frame, in SI units.
/// Open faces use the original native regularized reconstruction. They do not
/// describe finite reservoirs, material acoustic pressure, or pressure work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FaceBoundarySi {
    /// Both faces of this axis must select periodic wrapping.
    Periodic,
    /// Stationary halfway bounce-back wall, as admitted by immersed-body steps.
    StationaryWall,
    /// Prescribed material-fluid velocity [m/s] at an open face.
    Velocity {
        /// Cartesian components in the fixed world/grid frame [m/s].
        velocity_m_s: [f64; 3],
    },
    /// Gauge pressure perturbation [Pa] about the source's reference density.
    /// Uses the original numerical isothermal equation, not absolute gas pressure.
    Pressure {
        /// Signed incompressible perturbation about reference density [Pa].
        perturbation_pa: f64,
    },
}

impl FaceBoundarySi {
    /// A stationary wall with the original halfway bounce-back law.
    #[must_use]
    pub const fn stationary_wall() -> Self {
        Self::StationaryWall
    }
}

/// Explicit faces paired with the unchanged resolution/allocation controls.
/// All current material immersed constructors accept this or `ImmersedDomain`.
/// Accepted owners can update open target values through
/// `rebind_open_boundaries_si`; boundary kinds and topology remain fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImmersedBoundaryDomain {
    pub(super) domain: ImmersedDomain,
    faces: [FaceBoundarySi; 6],
}

impl ImmersedDomain {
    /// Override the legacy periodic/stationary selector with explicit SI faces.
    /// Order is x-min, x-max, y-min, y-max, z-min, z-max, as `Face3::ALL`.
    /// All other domain controls remain authoritative. Native face admission
    /// and SI conversion occur at construction before grid allocation.
    #[must_use]
    pub const fn with_boundaries(self, faces: [FaceBoundarySi; 6]) -> ImmersedBoundaryDomain {
        ImmersedBoundaryDomain {
            domain: self,
            faces,
        }
    }
}

impl From<ImmersedDomain> for ImmersedBoundaryDomain {
    fn from(domain: ImmersedDomain) -> Self {
        domain.with_boundaries(core::array::from_fn(|i| {
            if domain.periodic[i / 2] {
                FaceBoundarySi::Periodic
            } else {
                FaceBoundarySi::stationary_wall()
            }
        }))
    }
}

impl ImmersedBoundaryDomain {
    /// Authored SI faces in the same order as `Face3::ALL`.
    #[must_use]
    pub const fn faces(&self) -> &[FaceBoundarySi; 6] {
        &self.faces
    }

    pub(super) fn lattice_boundaries<F: FluidTransport + ?Sized>(
        &self,
        scaling: &MaterialScaling<'_, F>,
    ) -> Result<BoundarySpec3, Error> {
        lattice_boundaries(self.faces, scaling, self.domain.maximum_speed_m_s)
    }
}

/// One source-owned conversion for construction and live held targets.
pub(super) fn lattice_boundaries<F: FluidTransport + ?Sized>(
    faces_si: [FaceBoundarySi; 6],
    scaling: &MaterialScaling<'_, F>,
    maximum_speed_m_s: f64,
) -> Result<BoundarySpec3, Error> {
    let velocity = |value: [f64; 3]| -> Result<[f64; 3], Error> {
        let value = convert(value, 1.0)?;
        let speed = value[0].hypot(value[1]).hypot(value[2]);
        if speed > maximum_speed_m_s {
            return Err(Error::Invalid(
                "boundary velocity exceeds the source-scaled scenario speed bound",
            ));
        }
        Ok([
            scaling.velocity_lattice(value[0])?,
            scaling.velocity_lattice(value[1])?,
            scaling.velocity_lattice(value[2])?,
        ])
    };
    let mut faces = [FaceBoundary3::stationary_wall(); 6];
    for (to, from) in faces.iter_mut().zip(faces_si) {
        *to = match from {
            FaceBoundarySi::Periodic => FaceBoundary3::Periodic,
            FaceBoundarySi::StationaryWall => FaceBoundary3::stationary_wall(),
            FaceBoundarySi::Velocity { velocity_m_s } => FaceBoundary3::Velocity {
                velocity: velocity(velocity_m_s)?,
            },
            FaceBoundarySi::Pressure { perturbation_pa } => FaceBoundary3::Pressure {
                density: scaling.density_lattice_from_pressure_perturbation(perturbation_pa)?,
            },
        };
    }
    let boundaries = BoundarySpec3::new(faces);
    boundaries.try_validate().map_err(Error::Flow)?;
    Ok(boundaries)
}
