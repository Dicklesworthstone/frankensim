//! Original moving-wall exchange with cellwise source-resolved viscosity.

use super::*;
use crate::d3q19::freesurface3::{PlanarWallFieldControls3, PlanarWallVelocityField3};
use crate::material::free_surface::{
    MaterialPlanarWallFieldReport3, MaterialPlanarWallVelocityField3,
};

impl MaterialFreeSurfaceRheology3<'_, '_> {
    /// Drive a nonuniform tangential wall with the actual source viscosity.
    ///
    /// The accepted source viscosity at each wet cell supplies relaxation in
    /// the SAME native collision that reflects the prescribed wall field.
    /// Actual endpoint pressure and force-corrected strain then refresh the
    /// original constitutive law before either populations or source samples
    /// publish. Pressure-dependent laws keep their construction-time absolute
    /// datum; wall pressure is not substituted for interior liquid pressure.
    ///
    /// The returned SI face impulses, torques, couples and actuator work are
    /// measured from those actual reflected populations. This remains an
    /// isothermal material owner: no caloric data, viscous heating or hidden
    /// heat sink is introduced, and wall work is not a total-energy balance.
    /// Use the thermal rheology owner to evolve temperature and molecular heat.
    ///
    /// Geometry and field admission are the existing fully flooded, fixed
    /// planar-wall rules. Face velocities are held for one physical timestep.
    /// Normal motion, moving contact lines and simultaneous gas actuation are
    /// outside this entry. Local constitutive solves remain explicit in time;
    /// endpoint admission is not an implicit wall/fluid/material iteration.
    ///
    /// Field budgets precede the one native candidate clone. No unused
    /// thermal mass-route history is allocated. Invalid motion, late source
    /// or convergence refusal, unrepresentable SI receipts, excessive native
    /// speed and cancellation leave all accepted flow and source state intact.
    pub fn step_with_planar_wall_velocity_field_si(
        &mut self,
        motion: MaterialPlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<MaterialPlanarWallFieldReport3, Error> {
        poll(gate)?;
        let (velocities, origin) = self
            .flow
            .prepare_planar_wall_velocity_field_si(motion, controls, gate)?;
        let tau = self.local_relaxation(gate)?;
        let (lattice, exchange) = self
            .flow
            .lattice
            .planar_wall_field_candidate_with_relaxation(
                PlanarWallVelocityField3 {
                    wall: motion.wall,
                    moment_origin: origin,
                    tangential_velocity: &velocities,
                },
                controls,
                Some(&tau),
                gate,
            )?;
        let samples = self.sample_lattice(&lattice, gate)?;
        let report =
            self.flow
                .planar_wall_field_report_si(exchange, motion.moment_origin_m, gate)?;
        poll(gate)?;
        self.flow.lattice = lattice;
        self.samples = samples;
        Ok(report)
    }
}
