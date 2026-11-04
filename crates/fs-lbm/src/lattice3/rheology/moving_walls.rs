//! Source-resolved viscosity and prescribed Couette/channel walls in one step.

use super::*;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;

// A channel needs only its two prescribed plane velocities. Keeping these in
// fixed storage avoids per-cell wall maps, probe searches and step allocation.
// Native probe construction guarantees unique reflected links, so the exact
// link count proves that no moving or stationary partner face was omitted.
fn channel_walls<L: ViscousLattice3<Q>, const Q: usize>(
    grid: &Grid3<L, Q>,
    motion: &WallMotion3<'_, L, Q>,
    mut poll: impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<(usize, [[f64; 3]; 2]), GridError3<L::Error>> {
    poll()?;
    let mut normal = None;
    for (axis, &periodic) in grid.periodic.iter().enumerate() {
        if !periodic && normal.replace(axis).is_some() {
            return Err(GridError3::Input(
                "material channel requires one wall pair and two periodic axes",
            ));
        }
    }
    let normal = normal.ok_or(GridError3::Input(
        "material channel requires one wall pair and two periodic axes",
    ))?;
    for tile in grid.solid.chunks(64) {
        poll()?;
        if tile.iter().any(|&solid| solid) {
            return Err(GridError3::Input(
                "material channel requires an all-fluid rectangle",
            ));
        }
    }
    motion.admit(grid, &mut poll)?;
    let layers = L::E.iter().try_fold(0usize, |sum, direction| {
        sum.checked_add((direction[normal].unsigned_abs() as usize).min(grid.dims[normal]))
            .ok_or(GridError3::Input("material channel wall count overflow"))
    })?;
    let required = (grid.f.len() / grid.dims[normal])
        .checked_mul(layers)
        .ok_or(GridError3::Input("material channel wall count overflow"))?;
    if motion.probe.measured_links() != required {
        return Err(GridError3::Input(
            "material channel requires all native links of both walls",
        ));
    }
    let mut walls = [None; 2];
    for (i, (link, &velocity)) in motion.probe.links.iter().zip(motion.velocities).enumerate() {
        if i % 64 == 0 {
            poll()?;
        }
        let side = usize::from(L::E[link.incoming][normal] < 0);
        if let Some(expected) = walls[side] {
            if velocity != expected {
                return Err(GridError3::Input(
                    "material channel requires a uniform velocity on each wall",
                ));
            }
        } else {
            walls[side] = Some(velocity);
        }
    }
    poll()?;
    let [Some(lower), Some(upper)] = walls else {
        return Err(GridError3::Input(
            "material channel requires both native wall planes",
        ));
    };
    Ok((normal, [lower, upper]))
}

impl<L: ViscousLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Resolve material viscosity against the SAME prescribed halfway-wall
    /// velocities used by the original native collision/stream step.
    ///
    /// The all-fluid rectangular channel has one nonperiodic normal and two
    /// periodic tangents. `probe` includes every native link of both walls,
    /// including a stationary partner; `velocities` supplies one uniform,
    /// tangential lattice velocity per plane in that exact native link order.
    /// Geometry stays fixed. Neither normal motion nor partial faces is admitted.
    ///
    /// The existing quadratic wall stencil uses wall-relative velocities, so
    /// translation does not generate spurious shear or alter source viscosity.
    /// A one-cell gap uses the velocity difference across its two halfway walls.
    /// Source sampling, dx/dt/rho0 and collision-rate preservation follow
    /// [`Self::step_serial_with_material_rheology_workspace`].
    ///
    /// The ordinary rheology workspace is reused without additional buffers.
    /// The returned wall impulse/work measures this step's actual native
    /// stream; it is not a thermal-energy update. Source report, populations
    /// and step count publish only after both streaming and wall measurement
    /// succeed. Refusal preserves the previous accepted report and flow.
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_channel_material_rheology_and_moving_walls<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<(&'w RheologyStep3<L::Collision, F>, MovingWallMomentum3), RheologyError3<L::Error>>
    {
        let poll = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let walls = channel_walls(self, &WallMotion3 { probe, velocities }, poll)?;
        prepare_workspace_with_walls(self, reference, None, Some(walls), workspace, resolve, poll)?;
        let momentum = self.step_serial_with_moving_wall_momentum(
            cx,
            probe,
            velocities,
            Some(&workspace.trial.collisions),
        )?;
        Ok((workspace.publish(), momentum))
    }

    /// Pooled counterpart with deterministic serial material queries and the
    /// original native tile owner. Wall and rheology receipts describe one
    /// successful stream, and no fallible work follows native publication.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_channel_material_rheology_and_moving_walls<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<(&'w RheologyStep3<L::Collision, F>, MovingWallMomentum3), RheologyError3<L::Error>>
    {
        let poll = || gate_poll(gate);
        let walls = channel_walls(self, &WallMotion3 { probe, velocities }, poll)?;
        prepare_workspace_with_walls(self, reference, None, Some(walls), workspace, resolve, poll)?;
        let momentum = self.step_pooled_with_moving_wall_momentum(
            pool,
            gate,
            probe,
            velocities,
            Some(&workspace.trial.collisions),
        )?;
        Ok((workspace.publish(), momentum))
    }
}
