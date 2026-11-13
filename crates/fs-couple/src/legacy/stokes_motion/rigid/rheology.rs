//! Constitutive dispatch through the original finite body/flow endpoint solve.

use super::{PlanarRigidControls, PlanarRigidError, PlanarRigidMotion, PlanarRigidStep};
use fs_flux::ale::rigid::{
    RigidAleCandidate, RigidAleError, RigidAleFlow, RigidAleMotion, RigidAleMotionTolerance,
};
use fs_flux::{MaterialFluxRheologyPoint, MaterialNavierStokesBoundary, NavierStokesControls};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

impl<F: FluidTransport + ?Sized> PlanarRigidMotion<'_, '_, '_, F> {
    /// Close all three body momentum equations with the original material law
    /// resolved at every candidate's actual fluid strain.
    ///
    /// resolve receives the endpoint WORLD centroid, objective strain invariant
    /// [1/s] and immutable reference transport state. Captured temperature and
    /// composition may vary spatially but remain fixed during this mechanical
    /// interval. Every response must retain the original reference density.
    /// Mechanical pressure is not a thermodynamic state.
    ///
    /// The native nonlinear solver admits convection AND constitutive stress at
    /// its returned velocity. The coupled response search repeats that solve at
    /// each body trial, including geometric and world-source changes. It never
    /// substitutes a viscosity frozen at the initial body or fluid rate.
    ///
    /// Original accepted material owners are available from
    /// result.flow.rheology_samples(); querying them does not re-evaluate the
    /// source. No source response, body state, fluid state or clock is partially
    /// published when a constitutive response or mechanical step refuses.
    /// The callback must be repeatable for identical input; its own external
    /// side effects cannot be rolled back by this owner.
    ///
    /// Native no-slip/slip, momentum/work, search-budget, rotation-increment and
    /// cancellation semantics are identical to step(). Viscoelastic memory,
    /// thermal evolution, variable density and a whole-system passivity theorem
    /// are outside this generalized-Newtonian mechanical step.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_rheology<G: FluidTransport>(
        &mut self,
        timestep_s: f64,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialNavierStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: [f64; 2],
        external_torque_n_m: f64,
        fluid_controls: NavierStokesControls,
        controls: PlanarRigidControls,
        resolve: impl FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PlanarRigidStep<G>, PlanarRigidError> {
        self.step_impl(
            timestep_s,
            force_world,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            external_torque_n_m,
            fluid_controls,
            controls,
            Rheology(resolve),
            cancelled,
        )
    }
}

// Static dispatch keeps both public constitutive paths on one mechanical loop.
// The returned candidate borrows only the native flow owner; material responses
// are retained by value, never as references to a previous response probe.
pub(super) trait FlowTrial<F: FluidTransport + ?Sized, G> {
    #[allow(clippy::too_many_arguments)]
    fn candidate<'o, 's, 'm, 'a>(
        &mut self,
        flow: &'o mut RigidAleFlow<'s, 'm, 'a, F>,
        timestep_s: f64,
        motion: RigidAleMotion,
        tolerance: RigidAleMotionTolerance,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries_world: &[MaterialNavierStokesBoundary<'_>],
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<RigidAleCandidate<'o, 's, 'm, 'a, F, G>, RigidAleError>;
}

pub(super) struct Newtonian;

impl<F: FluidTransport + ?Sized> FlowTrial<F, ()> for Newtonian {
    fn candidate<'o, 's, 'm, 'a>(
        &mut self,
        flow: &'o mut RigidAleFlow<'s, 'm, 'a, F>,
        timestep_s: f64,
        motion: RigidAleMotion,
        tolerance: RigidAleMotionTolerance,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries_world: &[MaterialNavierStokesBoundary<'_>],
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<RigidAleCandidate<'o, 's, 'm, 'a, F>, RigidAleError> {
        flow.candidate(
            timestep_s,
            motion,
            tolerance,
            force_world,
            boundaries_world,
            controls,
            absolute_velocity_tolerance_m_s,
            cancelled,
        )
    }
}

struct Rheology<R>(R);

impl<F, G, R> FlowTrial<F, G> for Rheology<R>
where
    F: FluidTransport + ?Sized,
    G: FluidTransport,
    R: FnMut(MaterialFluxRheologyPoint, &FluidTransportState) -> Result<G, MaterialError>,
{
    fn candidate<'o, 's, 'm, 'a>(
        &mut self,
        flow: &'o mut RigidAleFlow<'s, 'm, 'a, F>,
        timestep_s: f64,
        motion: RigidAleMotion,
        tolerance: RigidAleMotionTolerance,
        force_world: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries_world: &[MaterialNavierStokesBoundary<'_>],
        controls: NavierStokesControls,
        absolute_velocity_tolerance_m_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<RigidAleCandidate<'o, 's, 'm, 'a, F, G>, RigidAleError> {
        flow.candidate_with_rheology(
            timestep_s,
            motion,
            tolerance,
            force_world,
            boundaries_world,
            controls,
            absolute_velocity_tolerance_m_s,
            &mut self.0,
            cancelled,
        )
    }
}
