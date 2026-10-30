//! Per-face structural loads from original accepted Cartesian Euler steps.
use super::{Error, ROUNDING, WallPressureImpulse, WallTractionImpulse, dot, unit};
use fs_flux::compressible::molecular::MolecularEuler2d;
use fs_flux::compressible::spatial::molecular::MolecularEuler3d;
use fs_flux::compressible::spatial::{BoundarySide as SpatialSide, Euler3d};
use fs_flux::compressible::{AxisBoundaries, BoundaryCondition, BoundarySide, Euler2d};
use fs_gas::multicomponent::dissipation::TransportLaw;

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
        let (axis, sign, normal) = planar_frame(flow, side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no accepted Cartesian wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        normal_only(impulse.mass, impulse.energy, impulse.momentum, axis)?;
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.grid().face_area(side),
            normal,
            sign * impulse.momentum[axis],
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
        let (axis, sign, normal) = spatial_frame(flow, side, axes)?;
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no accepted Cartesian wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        normal_only(impulse.mass, impulse.energy, impulse.momentum, axis)?;
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.grid().face_area(side),
            normal,
            sign * impulse.momentum[axis],
            exterior_pressure_pa,
        )
    }
}

impl WallTractionImpulse {
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
        let (axis, sign, normal) = planar_frame(accepted, side, axes)?;
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
        let (axis, sign, normal) = spatial_frame(accepted, side, axes)?;
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
    flow: &Euler2d,
    side: BoundarySide,
    axes: [[f64; 3]; 2],
) -> Result<(usize, f64, [f64; 3]), Error> {
    if !orthonormal(axes) {
        return Err(Error::Invalid(
            "Cartesian embedding requires orthonormal x/y axes",
        ));
    }
    let (axis, sign, pair, upper) = match side {
        BoundarySide::XMin => (0, -1.0, flow.boundaries().x, false),
        BoundarySide::XMax => (0, 1.0, flow.boundaries().x, true),
        BoundarySide::YMin => (1, -1.0, flow.boundaries().y, false),
        BoundarySide::YMax => (1, 1.0, flow.boundaries().y, true),
    };
    reflecting(pair, upper)?;
    Ok((axis, sign, axes[axis].map(|v| sign * v)))
}

fn spatial_frame(
    flow: &Euler3d,
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
    reflecting(flow.boundaries().side(side), sign > 0.0)?;
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
