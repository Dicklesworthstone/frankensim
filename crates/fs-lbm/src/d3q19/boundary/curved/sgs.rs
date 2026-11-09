//! Local LES on the actual off-lattice, prescribed moving-wall geometry.

use super::super::sgs::{BoundarySgsError3, Neighbour};
use super::{
    BoundaryGrid3, CancelGate, CurvedWall3, CurvedWallStep3, checkpoint, invalid, low_mach,
};
use crate::d3q19::CollisionModel3;
use crate::sgs3::{SgsConfig, SgsReport, SgsWorkspace3};

#[derive(Clone, Copy)]
enum Sample {
    Point { distance: f64, velocity: [f64; 3] },
    Velocity,
    Pressure,
}

impl Sample {
    fn native(neighbour: Neighbour) -> Self {
        match neighbour {
            Neighbour::Fluid(velocity) => Self::Point {
                distance: 1.0,
                velocity,
            },
            Neighbour::Wall(velocity) => Self::Point {
                distance: 0.5,
                velocity,
            },
            Neighbour::Velocity => Self::Velocity,
            Neighbour::Pressure => Self::Pressure,
        }
    }
}

// Derivative at the center of the quadratic through (-a, low), (0, center)
// and (b, high). Distances belong to the coordinate direction, so an axial
// curved link uses its actual fluid-to-wall fraction. Difference form avoids
// subtracting separately weighted large velocities under a uniform boost.
fn derivative(center: f64, low: Sample, high: Sample, component: usize) -> f64 {
    match (low, high) {
        (Sample::Pressure, _) | (_, Sample::Pressure) => 0.0,
        (
            Sample::Point {
                distance: a,
                velocity: low,
            },
            Sample::Point {
                distance: b,
                velocity: high,
            },
        ) => {
            (b / (a + b)) * ((center - low[component]) / a)
                + (a / (a + b)) * ((high[component] - center) / b)
        }
        (Sample::Velocity, Sample::Point { distance, velocity }) => {
            (velocity[component] - center) / distance
        }
        (Sample::Point { distance, velocity }, Sample::Velocity) => {
            (center - velocity[component]) / distance
        }
        (Sample::Velocity, Sample::Velocity) => 0.0,
    }
}

impl BoundaryGrid3 {
    /// Prepare local Smagorinsky, WALE or Vreman on the accepted curved geometry.
    ///
    /// Supply one physical inertial-frame velocity per selected link, in
    /// `wall.links()` order. These frozen values also drive the combined step;
    /// no callback is resampled between the gradient and wall impulse.
    /// Held uniform/local body forces enter the original Guo-centered velocity.
    ///
    /// An axial selected link places its wall sample at the actual fraction.
    /// The opposite fluid or wall sample and the cell center give a quadratic
    /// nonuniform derivative, exact for affine and quadratic velocity profiles.
    /// Unselected axes retain the original periodic, halfway and open-face
    /// stencils. Diagonal intersections remain owned by the BFL population and
    /// load operator; they are not projected onto fictitious axial wall points.
    ///
    /// Reuses `sgs_workspace()` without allocating a dense distance/velocity
    /// field. One ordered sparse cursor visits selected links once. Collision
    /// rates use native tile/lane order and preserve all non-shear controls.
    /// This is resolved-gradient LES, with no wall law or subgrid roughness.
    pub fn prepare_curved_sgs(
        &self,
        wall: &CurvedWall3,
        wall_velocity: &[[f64; 3]],
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
    ) -> core::result::Result<SgsReport, BoundarySgsError3> {
        self.prepare_curved_sgs_sampled(wall, wall_velocity, gate, config, workspace, |_, _| {
            Ok(self.collision_model)
        })
    }

    // The actual intersection stencil supplies BOTH molecular rheology and
    // SGS. Source errors retain their original type; no halfway-wall rate or
    // independently sampled wall velocity enters the material composition.
    #[allow(clippy::too_many_arguments)]
    fn prepare_curved_sgs_sampled<E: From<BoundarySgsError3>>(
        &self,
        wall: &CurvedWall3,
        wall_velocity: &[[f64; 3]],
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
        sample: impl FnMut([usize; 3], [[f64; 3]; 3]) -> core::result::Result<CollisionModel3, E>,
    ) -> core::result::Result<SgsReport, E> {
        wall.validate(self, gate).map_err(BoundarySgsError3::from)?;
        if wall_velocity.len() != wall.links.len() {
            return Err(BoundarySgsError3::from(invalid(
                "curved SGS requires one velocity per selected link",
            ))
            .into());
        }
        for &velocity in wall_velocity {
            checkpoint(gate).map_err(BoundarySgsError3::from)?;
            if !low_mach(velocity) {
                return Err(BoundarySgsError3::from(invalid(
                    "curved wall velocity is non-finite or outside low Mach",
                ))
                .into());
            }
        }
        let mut cursor = 0;
        let report = self.prepare_sgs_with_gradient(
            gate,
            config,
            workspace,
            |p, velocity| {
                let (tile, lane) = self.addr(p[0], p[1], p[2]);
                let mut axial = [None; 6];
                while let Some(link) = wall.links.get(cursor) {
                    if (link.geometry.link.tile, link.geometry.link.lane) != (tile, lane) {
                        break;
                    }
                    let direction = link.geometry.link.direction;
                    if (1..=6).contains(&direction) {
                        axial[direction - 1] = Some(Sample::Point {
                            distance: link.geometry.fraction,
                            velocity: wall_velocity[cursor],
                        });
                    }
                    cursor += 1;
                }
                let mut gradient = self.sgs_gradient(p, velocity);
                let center = velocity[tile * super::TILE_CELLS + lane];
                for axis in 0..3 {
                    if axial[2 * axis].is_none() && axial[2 * axis + 1].is_none() {
                        continue;
                    }
                    let low = axial[2 * axis].unwrap_or_else(|| {
                        Sample::native(self.sgs_neighbour(p, axis, false, velocity))
                    });
                    let high = axial[2 * axis + 1].unwrap_or_else(|| {
                        Sample::native(self.sgs_neighbour(p, axis, true, velocity))
                    });
                    for component in 0..3 {
                        gradient[component][axis] =
                            derivative(center[component], low, high, component);
                    }
                }
                gradient
            },
            sample,
        )?;
        debug_assert_eq!(cursor, wall.links.len());
        Ok(report)
    }

    /// Advance one off-lattice moving-wall step with explicit local SGS.
    ///
    /// Both closures and actual wall populations use the same supplied wall
    /// velocities. The native curved owner retains BFL interpolation, periodic
    /// and open transport, mass change, body-force impulse, torque and work.
    /// The returned pair contains those actual loads and the pre-step SGS
    /// diagnostics; modeled dissipation is not deposited into a thermal field.
    ///
    /// Geometry remains fixed for this tick. All closure, Mach, collision, load
    /// and cancellation gates precede publication. On refusal populations and
    /// topology admission remain unchanged; workspace and post/next scratch may
    /// change. Zero coefficient preserves native curved-step population and
    /// load bits. No second fluid clone or per-step turbulence array is needed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_curved_walls_and_sgs(
        &mut self,
        wall: &CurvedWall3,
        wall_velocity: &[[f64; 3]],
        moment_origin: [f64; 3],
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
    ) -> core::result::Result<(CurvedWallStep3, SgsReport), BoundarySgsError3> {
        let molecular = self.collision_model;
        self.step_with_curved_walls_and_sgs_sampled(
            wall,
            wall_velocity,
            moment_origin,
            gate,
            config,
            workspace,
            |_, _| Ok(molecular),
        )
    }

    // Material consumers stage their original source epoch through this same
    // collision/curved-stream/load transaction, then publish it after success.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_with_curved_walls_and_sgs_sampled<E: From<BoundarySgsError3>>(
        &mut self,
        wall: &CurvedWall3,
        wall_velocity: &[[f64; 3]],
        moment_origin: [f64; 3],
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
        sample: impl FnMut([usize; 3], [[f64; 3]; 3]) -> core::result::Result<CollisionModel3, E>,
    ) -> core::result::Result<(CurvedWallStep3, SgsReport), E> {
        checkpoint(gate).map_err(BoundarySgsError3::from)?;
        if !moment_origin.into_iter().all(f64::is_finite) {
            return Err(BoundarySgsError3::from(invalid(
                "curved wall moment origin must be finite",
            ))
            .into());
        }
        let report =
            self.prepare_curved_sgs_sampled(wall, wall_velocity, gate, config, workspace, sample)?;
        self.collide_sgs_prepared(gate, workspace)?;
        let mut cursor = 0;
        let loads = self
            .stream_curved_collided(
                wall,
                |_| {
                    let velocity = wall_velocity[cursor];
                    cursor += 1;
                    velocity
                },
                moment_origin,
                gate,
            )
            .map_err(BoundarySgsError3::from)?;
        Ok((loads, report))
    }
}
