//! Per-face structural loads from original accepted Cartesian Euler steps.
use super::{Error, ROUNDING, WallPressureImpulse, WallTractionImpulse, dot, unit};
use fs_flux::compressible::mixture::{MixtureEuler2d, reactive::ReactingEuler2d};
use fs_flux::compressible::molecular::MolecularEuler2d;
use fs_flux::compressible::spatial::mixture::{MixtureEuler3d, reactive::ReactingEuler3d};
use fs_flux::compressible::spatial::molecular::MolecularEuler3d;
use fs_flux::compressible::spatial::{Boundaries3d, BoundarySide as SpatialSide, Euler3d};
use fs_flux::compressible::{
    AxisBoundaries, Boundaries2d, BoundaryCondition, BoundarySide, Euler2d,
};
use fs_gas::multicomponent::dissipation::TransportLaw;

/// Borrow only the original accepted interval and physical face transfer.
struct AcceptedWallFace<'a, const N: usize> {
    duration_s: f64,
    end_time_s: f64,
    area_m2: f64,
    species_masses: &'a [f64],
    momentum: [f64; N],
    energy: f64,
}

impl WallPressureImpulse {
    /// Consume one face of the actual latest accepted Cartesian wall step.
    /// `axes` embeds the x/y grid axes as orthonormal vectors in the structural
    /// frame. Face order and area remain the original grid's, so a nonuniform
    /// shock footprint can produce a moment instead of being spatially averaged.
    /// The caller assigns this explicit face index to the matching structural
    /// patch. No interpolation, geometric search or matching certificate is
    /// implied by equal area/normal admission.
    pub fn from_euler2d(
        flow: &Euler2d,
        side: BoundarySide,
        face: usize,
        axes: [[f64; 3]; 2],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let frame = planar_frame(flow.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no accepted Cartesian wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        Self::from_cartesian_face(
            AcceptedWallFace {
                duration_s: report.elapsed.value(),
                end_time_s: flow.time().value(),
                area_m2: flow.grid().face_area(side),
                species_masses: std::slice::from_ref(&impulse.mass),
                momentum: impulse.momentum,
                energy: impulse.energy,
            },
            frame,
            exterior_pressure_pa,
        )
    }

    /// Consume one original boundary face of the latest accepted 3-D Euler step.
    ///
    /// `axes` maps the original Cartesian x/y/z axes to orthonormal vectors in
    /// the structural frame. The source's physical face area and face order
    /// are retained. The caller explicitly matches each face to a structural
    /// patch; this is not an automatic geometric correspondence or remap.
    ///
    /// Both accepted RK stages contribute to this impulse. Reflecting faces
    /// alone are admitted: pressure cannot be extracted from transported
    /// momentum at an open boundary. No gas state is advanced or replaced.
    pub fn from_euler3d(
        flow: &Euler3d,
        side: SpatialSide,
        face: usize,
        axes: [[f64; 3]; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let frame = spatial_frame(flow.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no accepted Cartesian wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        Self::from_cartesian_face(
            AcceptedWallFace {
                duration_s: report.elapsed.value(),
                end_time_s: flow.time().value(),
                area_m2: flow.grid().face_area(side),
                species_masses: std::slice::from_ref(&impulse.mass),
                momentum: impulse.momentum,
                energy: impulse.energy,
            },
            frame,
            exterior_pressure_pa,
        )
    }

    /// Actual accepted composition-local pressure impulse on one planar face.
    /// Every species wall transfer must vanish individually; a zero total mass
    /// flux cannot hide counterflow. Geometry and source face order are retained.
    pub fn from_mixture_euler2d(
        flow: &MixtureEuler2d,
        side: BoundarySide,
        face: usize,
        axes: [[f64; 3]; 2],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let frame = planar_frame(flow.boundaries(), side, axes)?;
        let report = flow.last_step().ok_or(Error::Invalid(
            "no accepted Cartesian mixture wall interval",
        ))?;
        let impulse = flow
            .boundary_impulses(side)
            .nth(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        Self::from_cartesian_face(
            AcceptedWallFace {
                duration_s: report.elapsed.value(),
                end_time_s: flow.time().value(),
                area_m2: flow.grid().face_area(side),
                species_masses: impulse.species_masses,
                momentum: impulse.momentum,
                energy: impulse.energy,
            },
            frame,
            exterior_pressure_pa,
        )
    }

    /// Original accepted 3-D mixture wall pressure, including both RK stages.
    ///
    /// All species and three momenta share the original source Riemann solver;
    /// this adapter retains its per-face impulse, actual clock and physical
    /// area. `axes` is the explicit orthonormal structural-frame embedding.
    pub fn from_mixture_euler3d(
        flow: &MixtureEuler3d,
        side: SpatialSide,
        face: usize,
        axes: [[f64; 3]; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let frame = spatial_frame(flow.boundaries(), side, axes)?;
        let report = flow.last_step().ok_or(Error::Invalid(
            "no accepted Cartesian mixture wall interval",
        ))?;
        let impulse = flow
            .boundary_impulses(side)
            .nth(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        Self::from_cartesian_face(
            AcceptedWallFace {
                duration_s: report.elapsed.value(),
                end_time_s: flow.time().value(),
                area_m2: flow.grid().face_area(side),
                species_masses: impulse.species_masses,
                momentum: impulse.momentum,
                energy: impulse.energy,
            },
            frame,
            exterior_pressure_pa,
        )
    }

    /// Pressure load from the complete accepted chemistry/flow/chemistry step.
    ///
    /// Uses the outer accepted interval and its original per-face transfers.
    /// The middle transport report does not describe the final reacting state,
    /// and endpoint pressure times duration is not the accepted wall impulse.
    /// Only reflecting stationary adiabatic slip faces are admitted.
    pub fn from_reacting_euler3d(
        flow: &ReactingEuler3d,
        side: SpatialSide,
        face: usize,
        axes: [[f64; 3]; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let accepted = flow.flow();
        let frame = spatial_frame(accepted.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no complete reacting wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .nth(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        Self::from_cartesian_face(
            AcceptedWallFace {
                duration_s: report.elapsed.value(),
                end_time_s: accepted.time().value(),
                area_m2: accepted.grid().face_area(side),
                species_masses: impulse.species_masses,
                momentum: impulse.momentum,
                energy: impulse.energy,
            },
            frame,
            exterior_pressure_pa,
        )
    }

    /// Planar counterpart consuming the complete original reacting interval.
    /// Local chemistry changes the final calorics but carries no wall mass or
    /// energy flux. The outer report and its accepted face account determine
    /// the structural pressure-time integral, including the actual duration.
    /// An interval containing molecular stress requires the signed
    /// [`WallTractionImpulse::from_reacting_euler2d`] adapter instead.
    pub fn from_reacting_euler2d(
        flow: &ReactingEuler2d,
        side: BoundarySide,
        face: usize,
        axes: [[f64; 3]; 2],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let accepted = flow.flow();
        let frame = planar_frame(accepted.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no complete reacting wall interval"))?;
        if report.dissipation.is_some() {
            return Err(Error::Invalid(
                "molecular wall stress requires signed traction, not pressure",
            ));
        }
        let impulse = flow
            .boundary_impulses(side)
            .nth(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        Self::from_cartesian_face(
            AcceptedWallFace {
                duration_s: report.elapsed.value(),
                end_time_s: accepted.time().value(),
                area_m2: accepted.grid().face_area(side),
                species_masses: impulse.species_masses,
                momentum: impulse.momentum,
                energy: impulse.energy,
            },
            frame,
            exterior_pressure_pa,
        )
    }

    fn from_cartesian_face<const N: usize>(
        face: AcceptedWallFace<'_, N>,
        (axis, sign, normal): (usize, f64, [f64; 3]),
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        species_normal_only(face.species_masses, face.energy, face.momentum, axis)?;
        Self::new(
            face.duration_s,
            face.end_time_s,
            face.area_m2,
            normal,
            sign * face.momentum[axis],
            exterior_pressure_pa,
        )
    }
}

impl WallTractionImpulse {
    /// Full signed normal traction from the accepted planar reacting interval.
    ///
    /// When molecular transport is enabled, this retains BOTH molecular half
    /// steps and the wave stage, using the outer duration, clock and combined
    /// face account. The gas contribution is `p - tau_nn` and may be tensile.
    /// The same adapter admits inviscid intervals, where the viscous term is zero.
    /// Every species, tangent traction and wall-energy transfer must vanish.
    pub fn from_reacting_euler2d(
        flow: &ReactingEuler2d,
        side: BoundarySide,
        face: usize,
        axes: [[f64; 3]; 2],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let accepted = flow.flow();
        let (axis, sign, normal) = planar_frame(accepted.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no complete reacting wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .nth(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        species_normal_only(
            impulse.species_masses,
            impulse.energy,
            impulse.momentum,
            axis,
        )?;
        Self::new(
            report.elapsed.value(),
            accepted.time().value(),
            accepted.grid().face_area(side),
            normal,
            sign * impulse.momentum[axis],
            exterior_pressure_pa,
        )
    }

    /// Retain pressure and viscous normal stress from the complete accepted
    /// planar molecular/Euler split, including BOTH molecular half steps.
    ///
    /// The underlying Euler report describes only the middle advective stage;
    /// it is deliberately not used here. Stationary adiabatic slip walls may
    /// carry tensile total normal traction. The caller matches each original
    /// face to a structural patch in the orthonormal `axes` frame.
    pub fn from_molecular_euler2d<L: TransportLaw>(
        flow: &MolecularEuler2d<L>,
        side: BoundarySide,
        face: usize,
        axes: [[f64; 3]; 2],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let accepted = flow.flow();
        let (axis, sign, normal) = planar_frame(accepted.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no complete molecular wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        normal_only(impulse.mass, impulse.energy, impulse.momentum, axis)?;
        Self::new(
            report.elapsed.value(),
            accepted.time().value(),
            accepted.grid().face_area(side),
            normal,
            sign * impulse.momentum[axis],
            exterior_pressure_pa,
        )
    }

    /// Consume one original 3-D face of the full accepted molecular/Euler
    /// interval, preserving signed `p - tau_nn` and its exact duration.
    ///
    /// Uses the combined molecular boundary account, never an endpoint pressure
    /// or the underlying middle Euler stage. The caller supplies the original
    /// face-to-patch correspondence and orthonormal structural-frame embedding.
    pub fn from_molecular_euler3d<L: TransportLaw>(
        flow: &MolecularEuler3d<L>,
        side: SpatialSide,
        face: usize,
        axes: [[f64; 3]; 3],
        exterior_pressure_pa: f64,
    ) -> Result<Self, Error> {
        let accepted = flow.flow();
        let (axis, sign, normal) = spatial_frame(accepted.boundaries(), side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no complete molecular wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        normal_only(impulse.mass, impulse.energy, impulse.momentum, axis)?;
        Self::new(
            report.elapsed.value(),
            accepted.time().value(),
            accepted.grid().face_area(side),
            normal,
            sign * impulse.momentum[axis],
            exterior_pressure_pa,
        )
    }
}

fn planar_frame(
    boundaries: Boundaries2d,
    side: BoundarySide,
    axes: [[f64; 3]; 2],
) -> Result<(usize, f64, [f64; 3]), Error> {
    if !orthonormal(axes) {
        return Err(Error::Invalid(
            "Cartesian embedding requires orthonormal x/y axes",
        ));
    }
    let (axis, sign, pair, upper) = match side {
        BoundarySide::XMin => (0, -1.0, boundaries.x, false),
        BoundarySide::XMax => (0, 1.0, boundaries.x, true),
        BoundarySide::YMin => (1, -1.0, boundaries.y, false),
        BoundarySide::YMax => (1, 1.0, boundaries.y, true),
    };
    reflecting(pair, upper)?;
    Ok((axis, sign, axes[axis].map(|v| sign * v)))
}

fn spatial_frame(
    boundaries: Boundaries3d,
    side: SpatialSide,
    axes: [[f64; 3]; 3],
) -> Result<(usize, f64, [f64; 3]), Error> {
    if !orthonormal(axes) {
        return Err(Error::Invalid(
            "Cartesian embedding requires orthonormal x/y/z axes",
        ));
    }
    let axis = side.axis();
    let sign = side.sign();
    reflecting(boundaries.side(side), sign > 0.0)?;
    Ok((axis, sign, axes[axis].map(|v| sign * v)))
}

fn normal_only<const N: usize>(
    mass: f64,
    energy: f64,
    momentum: [f64; N],
    normal: usize,
) -> Result<(), Error> {
    if mass != 0.0
        || energy != 0.0
        || momentum
            .into_iter()
            .enumerate()
            .any(|(axis, value)| axis != normal && value != 0.0)
    {
        return Err(Error::Invalid(
            "stationary adiabatic slip wall transfers only normal momentum",
        ));
    }
    Ok(())
}

fn species_normal_only<const N: usize>(
    species_masses: &[f64],
    energy: f64,
    momentum: [f64; N],
    normal: usize,
) -> Result<(), Error> {
    if species_masses.iter().any(|&mass| mass != 0.0) {
        return Err(Error::Invalid(
            "stationary adiabatic wall transfers no individual species",
        ));
    }
    normal_only(0.0, energy, momentum, normal)
}

fn reflecting(pair: AxisBoundaries, upper: bool) -> Result<(), Error> {
    let condition = match pair {
        AxisBoundaries::Ends { lower, upper: high } => {
            if upper {
                high
            } else {
                lower
            }
        }
        AxisBoundaries::Periodic => {
            return Err(Error::Invalid("periodic faces have no wall load"));
        }
    };
    if condition != BoundaryCondition::Reflecting {
        return Err(Error::Invalid(
            "only a reflecting face supplies stationary support pressure",
        ));
    }
    Ok(())
}

fn orthonormal<const N: usize>(axes: [[f64; 3]; N]) -> bool {
    axes.iter().all(|&axis| unit(axis))
        && (0..N).all(|i| (0..i).all(|j| dot(axes[i], axes[j]).abs() <= ROUNDING))
}
