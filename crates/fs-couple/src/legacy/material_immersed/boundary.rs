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

/// Control of the prescribed trace's signed outward-normal mean velocity.
/// This constrains the velocity profile, not density-weighted mass transport.
/// The accepted interval's native `open` receipt measures actual exchange.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum InletFluxSi {
    /// Retain the complete sampled profile, including its bulk flow.
    Unconstrained,
    /// Shift every normal component to impose this plane mean [m/s].
    /// Inflow is negative on either face. Mixed wall/open seams are refused
    /// because their partial reconstruction cannot impose all velocity moments.
    MeanNormalVelocity(f64),
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

/// Convert a spatial SI trace on the same source-owned grid and fixed dx/dt.
/// Numerical admission and both reusable profile buffers belong to the native owner.
pub(super) fn set_velocity_profile_si<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    scaling: &MaterialScaling<'_, F>,
    maximum_speed_m_s: f64,
    face: Face3,
    flux: InletFluxSi,
    gate: &CancelGate,
    mut sample: impl FnMut([f64; 3]) -> Result<[f64; 3], BoundaryProfileError3>,
) -> Result<VelocityProfileReport3, Error> {
    if gate.is_requested() {
        return Err(Error::Profile(BoundaryProfileError3::Cancelled));
    }
    let flux = match flux {
        InletFluxSi::Unconstrained => InletFlux3::Unconstrained,
        InletFluxSi::MeanNormalVelocity(mean) => {
            InletFlux3::MeanNormal(scaling.velocity_lattice(mean)?)
        }
    };
    let maximum_speed_lu = scaling.velocity_lattice(maximum_speed_m_s)?;
    grid.set_velocity_profile_with_speed_limit(face, flux, maximum_speed_lu, gate, |cell| {
        let position = convert(cell.map(|i| i as f64 + 0.5), scaling.cell_size_m())
            .map_err(|error| BoundaryProfileError3::Source(error.to_string()))?;
        let value = sample(position)?;
        let component = |v| {
            scaling
                .velocity_lattice(v)
                .map_err(|error| BoundaryProfileError3::Source(error.to_string()))
        };
        Ok([component(value[0])?, component(value[1])?, component(value[2])?])
    })
    .map_err(Error::Profile)
}


/// Convert a held physical pressure-perturbation trace with the original
/// reference density and numerical isothermal pressure mapping.
pub(super) fn set_pressure_profile_si<F: FluidTransport + ?Sized>(
    grid: &mut BoundaryGrid3,
    scaling: &MaterialScaling<'_, F>,
    face: Face3,
    gate: &CancelGate,
    mut sample: impl FnMut([f64; 3]) -> Result<f64, BoundaryProfileError3>,
) -> Result<PressureProfileReport3, Error> {
    if gate.is_requested() {
        return Err(Error::Profile(BoundaryProfileError3::Cancelled));
    }
    grid.set_pressure_profile(face, gate, |cell| {
        let position = convert(cell.map(|i| i as f64 + 0.5), scaling.cell_size_m())
            .map_err(|error| BoundaryProfileError3::Source(error.to_string()))?;
        scaling
            .density_lattice_from_pressure_perturbation(sample(position)?)
            .map_err(|error| BoundaryProfileError3::Source(error.to_string()))
    })
    .map_err(Error::Profile)
}
