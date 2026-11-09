//! Original molecular rheology and SGS at the actual off-lattice wall geometry.

use super::{
    BoundaryGrid3, BoundaryRheologySgsWorkspace3, BoundarySgsReport3, CancelGate, FluidTransport,
    MaterialBoundarySgsError3, MaterialError, MaterialImmersedPoint3, MaterialScaling, SgsConfig,
};
use crate::d3q19::{CurvedWall3, CurvedWallStep3};
use crate::sgs3::SgsReport;

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Advance original curved-wall flow with source molecular rheology and SGS.
    ///
    /// Reuse [`Self::boundary_rheology_sgs_workspace`]. Each fluid cell is
    /// sampled once at the incoming physical deviatoric shear rate [1/s]. The
    /// SAME Guo-centered gradient supplies that rate and the full-gradient SGS
    /// kernel. Axial selected wall links use their actual intersection fractions
    /// and the supplied inertial-frame wall velocities, including unequal
    /// distances; all other axes retain the native boundary stencil.
    ///
    /// Supply one lattice-unit wall velocity per `wall.links()` entry. Those
    /// exact frozen values drive both gradient reconstruction and the original
    /// BFL population/impulse/work operator. No wall callback is resampled, and
    /// diagonal intersections are never replaced with invented axial points.
    /// The geometry and occupancy stay fixed for this tick.
    ///
    /// The original source response, physical scaling, reference-density
    /// requirement, resolved-shear molecular approximation and thermal limits
    /// are those of [`Self::step_boundary_with_material_rheology_and_sgs`].
    /// Actual local molecular viscosity sets the eddy/molecular ratio gate.
    /// Non-shear collision controls remain unchanged, including an independently
    /// authored bulk rate. Molecular and modeled dissipation remain separate
    /// instantaneous diagnostics; neither is silently deposited as heat.
    ///
    /// The complete returned native receipt retains curved and other-wall
    /// exchanges, BFL mass change, held body force, and open/periodic transport.
    /// `workspace.last_report()` keeps the SGS report and selected curved-wall
    /// exchange. Source samples and that report publish only after the native
    /// population swap succeeds. Any source, closure, wall, collision, load or
    /// cancellation refusal preserves both populations and the applied source
    /// epoch. Source callbacks must not publish external constitutive history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_curved_walls_with_material_rheology_and_sgs<M: FluidTransport>(
        &self,
        grid: &mut BoundaryGrid3,
        wall: &CurvedWall3,
        wall_velocity: &[[f64; 3]],
        moment_origin: [f64; 3],
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut BoundaryRheologySgsWorkspace3<M>,
        mut resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<(CurvedWallStep3, SgsReport), MaterialBoundarySgsError3> {
        self.admit_boundary_rheology_sgs(grid, gate, workspace)?;
        let pending = &mut workspace.pending;
        let model = workspace.model;
        let (momentum, sgs) = grid.step_with_curved_walls_and_sgs_sampled(
            wall,
            wall_velocity,
            moment_origin,
            gate,
            config,
            &mut workspace.flow,
            |cell, gradient| -> Result<_, MaterialBoundarySgsError3> {
                let (collision, sample) =
                    self.sample_boundary_rheology_sgs(model, cell, gradient, gate, &mut resolve)?;
                pending.push(sample);
                Ok(collision)
            },
        )?;
        workspace.accept(BoundarySgsReport3 {
            sgs,
            wall_momentum: Some(momentum.curved),
        });
        Ok((momentum, sgs))
    }
}
