//! Per-face spatial loads from the original two-dimensional accepted Euler step.
use super::{Error, ROUNDING, WallPressureImpulse, dot, unit};
use fs_flux::compressible::{AxisBoundaries, BoundaryCondition, BoundarySide, Euler2d};

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
        if !unit(axes[0]) || !unit(axes[1]) || dot(axes[0], axes[1]).abs() > ROUNDING {
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
        let report = flow
            .last_step()
            .ok_or(Error::Invalid("no accepted Cartesian wall interval"))?;
        let impulse = flow
            .boundary_impulses(side)
            .get(face)
            .ok_or(Error::Invalid(
                "wall face index outside the original Cartesian side",
            ))?;
        if impulse.mass != 0.0 || impulse.energy != 0.0 || impulse.momentum[1 - axis] != 0.0 {
            return Err(Error::Invalid(
                "stationary inviscid wall transfers only normal momentum",
            ));
        }
        Self::new(
            report.elapsed.value(),
            flow.time().value(),
            flow.grid().face_area(side),
            axes[axis].map(|v| sign * v),
            sign * impulse.momentum[axis],
            exterior_pressure_pa,
        )
    }
}
