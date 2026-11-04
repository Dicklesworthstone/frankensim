//! One sourced viscosity field for simultaneous moving immersed bodies.

use super::*;
use crate::d3q19::MovingPartialSaturationBodiesExchange3;

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Advance an actual moving rigid assembly with one original material epoch.
    ///
    /// The existing incoming physical velocity gradient is evaluated once for
    /// the whole grid, independent of body count. Each non-solid cell supplies
    /// one original source response at `sqrt(2 D:D) / dt`. Its local viscosity
    /// controls fluid collision, every body's partial-saturation weight, and
    /// the simultaneous implicit rigid response. Bodies occupying disjoint
    /// subvolumes of one cell share its same total-occupancy denominator.
    ///
    /// Both native free half-drifts, all shape samples and all final pose,
    /// momentum, force-history and cancellation checks precede the one fluid
    /// publication. The one-body array follows the original single-body law.
    /// There is no sequential per-body fluid step or repeated material query.
    /// Native pair and wall/halo admission still apply; scratch may change on
    /// failure, while all accepted bodies, populations and force state remain.
    ///
    /// Inputs and momentum accounts use native lattice units. Reference rho0,
    /// dx and dt retain the same conversion and incompressible density rule as
    /// [`Self::step_moving_body_with_material_rheology`]. Sources are held at
    /// the incoming time level; temperature/composition remain caller-owned.
    /// This does not infer subgrid film viscosity, dry-body displaced-fluid
    /// inertia, thermal feedback, or an implicit endpoint constitutive solve.
    /// Each supplied solid inertia remains additional to retained fluid.
    #[allow(clippy::too_many_arguments)]
    pub fn step_moving_bodies_with_material_rheology<M: FluidTransport, const N: usize>(
        &self,
        grid: &mut BoundaryGrid3,
        samplers: &mut [BodyFrameSampler3; N],
        shapes: [&dyn BodyFrameShape3; N],
        bodies: &mut [MovingPartialSaturationBody3; N],
        external_impulses: [[f64; 6]; N],
        moment_origin: [f64; 3],
        workspace: &mut KinematicsWorkspace3,
        max_cells: usize,
        endpoint_force_scale: Option<f64>,
        gate: &CancelGate,
        resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<
        MaterialImmersedStep3<M, MovingPartialSaturationBodiesExchange3<N>>,
        MaterialImmersedError3,
    > {
        if N == 0 {
            return Err(MaterialImmersedError3::Invalid(
                "immersed material assembly must contain at least one body",
            ));
        }
        let prepared = self.prepare_immersed_rheology(grid, workspace, max_cells, gate, resolve)?;
        let momentum = grid.step_with_moving_partial_saturation_bodies_local_relaxation(
            samplers,
            shapes,
            bodies,
            external_impulses,
            moment_origin,
            gate,
            endpoint_force_scale,
            &prepared.relaxation,
        )?;
        Ok(MaterialImmersedStep3 {
            momentum,
            applied: prepared.applied,
        })
    }
}
