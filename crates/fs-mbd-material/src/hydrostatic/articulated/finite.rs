//! One conserved liquid inventory shared by every link of a mechanism.
use super::{ArticulatedHydrostaticError, PreparedArticulatedHydrostatics, link_wrench};
use crate::hydrostatic::{
    FiniteHydrostaticTank, FiniteLayeredHydrostaticTank, FiniteTankBodiesLoad, FiniteTankBody,
    HydrostaticError, MAX_TANK_BODIES, poll,
};
use fs_ga::{Se3, Twist, Wrench};
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::RigidBodyState;

/// Articulated pressure loading and ONE shared finite-liquid energy account.
/// Supply the tank's gravity once to native articulated dynamics for body weight.
#[derive(Clone, Debug, PartialEq)]
pub struct ArticulatedFiniteTankLoad {
    /// Native external body wrenches in authored-link axes, about link origins,
    /// in exactly the original model's order. They include buoyancy alone.
    pub external_body_wrenches: Vec<Wrench>,
    /// Original common surface, liquid potential and per-link mass moments.
    /// The surface potential is the complete liquid gravitational-energy change
    /// from the same inventories in an empty tank. Count it once, alongside
    /// body energy; there are no independently additive link tank potentials.
    pub tank_load: FiniteTankBodiesLoad,
}

impl PreparedArticulatedHydrostatics<'_> {
    /// Displace one homogeneous finite liquid with the WHOLE mechanism.
    ///
    /// Native forward kinematics supplies every link pose before the existing
    /// shared-tank solver evaluates A*h - sum(V_sub,link(h)) = V_liquid.
    /// Moving one link can therefore change the force on another. The same
    /// accepted cuts supply all wrenches and one liquid potential, including
    /// the cross-link terms in total displaced volume squared.
    ///
    /// This path supports at most MAX_TANK_BODIES links and retains the global
    /// cell budget. Envelopes must remain disjoint, above the stationary floor
    /// and strictly inside the side walls. Overflow, contact and a sealing
    /// piston refuse as in the original finite-tank owner. The tank is fixed:
    /// a floating mechanism does not turn it into a moving-container model.
    pub fn load_finite_tank<F: FluidTransport + ?Sized>(
        &self,
        tank: &FiniteHydrostaticTank<'_, F>,
        world_from_base: Se3,
        joint_position: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ArticulatedFiniteTankLoad, ArticulatedHydrostaticError> {
        self.load_shared_tank(
            world_from_base,
            joint_position,
            cancelled,
            |bodies, states, cancel| tank.load_bodies(bodies, states, cancel),
        )
    }

    /// Displace one finite stack of conserved, stably ordered immiscible liquids.
    ///
    /// Each cumulative inventory solves one common interface across ALL links;
    /// positive density jumps combine the original exact mass moments. The
    /// returned liquid potential covers the complete stack exactly once.
    /// Equal-density layers retain the native uniform-tank reduction. This
    /// does not add mixing, waves, fluid inertia, drag or added mass.
    pub fn load_finite_layered_tank(
        &self,
        tank: &FiniteLayeredHydrostaticTank<'_>,
        world_from_base: Se3,
        joint_position: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ArticulatedFiniteTankLoad, ArticulatedHydrostaticError> {
        self.load_shared_tank(
            world_from_base,
            joint_position,
            cancelled,
            |bodies, states, cancel| tank.load_bodies(bodies, states, cancel),
        )
    }

    fn load_shared_tank<C, L>(
        &self,
        world_from_base: Se3,
        joint_position: &[f64],
        mut cancelled: C,
        evaluate: L,
    ) -> Result<ArticulatedFiniteTankLoad, ArticulatedHydrostaticError>
    where
        C: FnMut() -> bool,
        L: FnOnce(
            &[FiniteTankBody<'_>],
            &[RigidBodyState],
            &mut C,
        ) -> Result<FiniteTankBodiesLoad, HydrostaticError>,
    {
        poll(&mut cancelled)?;
        if self.frames.len() > MAX_TANK_BODIES {
            return Err(HydrostaticError::Invalid(
                "shared articulated tank exceeds the native body limit",
            )
            .into());
        }
        let poses = self.link_poses(world_from_base, joint_position, &mut cancelled)?;
        let mut states = Vec::with_capacity(self.frames.len());
        let mut bodies = Vec::with_capacity(self.frames.len());
        for ((frame, envelope), pose) in self.frames.iter().zip(self.displacement).zip(poses) {
            poll(&mut cancelled)?;
            states.push(
                frame
                    .state_from_kinematics(pose, Twist::zero())
                    .map_err(HydrostaticError::from)?,
            );
            bodies.push(FiniteTankBody { envelope, frame });
        }
        // Exactly one shared solve; never call a single-body tank per link.
        let tank_load = evaluate(&bodies, &states, &mut cancelled)?;
        let mut external_body_wrenches = Vec::with_capacity(self.frames.len());
        for ((frame, &state), load) in self.frames.iter().zip(&states).zip(&tank_load.bodies) {
            poll(&mut cancelled)?;
            external_body_wrenches.push(link_wrench(frame, state, load.wrench())?);
        }
        poll(&mut cancelled)?;
        Ok(ArticulatedFiniteTankLoad {
            external_body_wrenches,
            tank_load,
        })
    }
}
