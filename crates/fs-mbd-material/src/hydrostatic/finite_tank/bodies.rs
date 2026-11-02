//! Several native rigid bodies displacing one conserved liquid inventory.
use super::{
    DisplacedVolume, DisplacementMesh, FiniteHydrostaticTank, FiniteTankSurface,
    PreparedDisplacement,
};
use crate::hydrostatic::{
    HydrostaticBuoyancy, HydrostaticError, MAX_DISPLACEMENT_CELLS, finite, poll,
};
use crate::principal::PrincipalFrame;
use fs_ga::Vec3;
use fs_material::fluid_transport::FluidTransport;
use fs_mbd::RigidBodyState;

mod step;
pub use step::FiniteTankBodiesStep;
pub(in crate::hydrostatic) use step::step_bodies_with;

/// Hard bound on independently moving envelopes in one finite tank.
pub const MAX_TANK_BODIES: usize = 64;

/// One immutable material frame and its fluid-excluding geometry.
/// Bodies and states use the same stable order throughout a motion history.
#[derive(Clone, Copy)]
pub struct FiniteTankBody<'a> {
    /// Actual closed displacement envelope in the original body frame.
    pub envelope: &'a DisplacementMesh,
    /// Existing material-derived mass, centroid and full-inertia frame.
    pub frame: &'a PrincipalFrame,
}

/// Forces on every body below one free surface, with one fluid-energy account.
/// For a layered tank the energy covers the complete immiscible liquid stack.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteTankBodiesLoad {
    /// Top free surface, total displaced volume and complete liquid storage.
    pub surface: FiniteTankSurface,
    /// Individual original-frame mass moments and native buoyancy wrenches.
    /// These omit body weight and do not duplicate the shared fluid potential.
    pub bodies: Vec<HydrostaticBuoyancy>,
}

/// Geometry admission and cuts shared by every cumulative liquid inventory.
pub(in crate::hydrostatic) struct PreparedTankBodies<'a> {
    envelopes: Vec<PreparedDisplacement<'a>>,
    whole_volume_m3: f64,
    pub(in crate::hydrostatic) cuts: Vec<DisplacedVolume>,
}

impl<F: FluidTransport + ?Sized> FiniteHydrostaticTank<'_, F> {
    /// Solve one finite inventory displaced by ALL supplied body poses.
    ///
    /// A*h - sum(V_sub,i) = V_liquid determines the common waterline. The
    /// original analytic clipping supplies each body's force and torque; the
    /// original scalar solve accounts for the shared gravitational potential
    /// rho*g*(sum(I_depth,i) - sum(V_sub,i)^2/(2*A)) exactly once.
    ///
    /// Work is bounded by 64 bodies and the original total displacement-cell
    /// budget per complete surface evaluation. Prepared frame maps and clipped
    /// results use O(bodies) storage, reused across root trials; no combined mesh
    /// is rebuilt and no source is queried. Dry/full envelope caches still apply.
    ///
    /// Envelopes must be physically disjoint. As for each input tetrahedral
    /// partition, numeric admission does not prove overlap or coverage. Contact,
    /// collision, fluid inertia and waves require other owners. Each envelope
    /// must remain above the tank floor and strictly inside its side walls.
    /// Failure or cancellation returns no partial set of forces.
    pub fn load_bodies(
        &self,
        bodies: &[FiniteTankBody<'_>],
        states: &[RigidBodyState],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FiniteTankBodiesLoad, HydrostaticError> {
        let mut prepared = self.prepare_bodies(bodies, states, &mut cancelled)?;
        let surface = self.fill_prepared_bodies(&mut prepared, &mut cancelled)?;
        // Accepted root samples remain in cuts: no extra clipping pass or
        // independent per-body surface solve is needed for the actual wrench.
        let mut loads = Vec::with_capacity(bodies.len());
        for ((body, &state), cut) in bodies.iter().zip(states).zip(prepared.cuts) {
            poll(&mut cancelled)?;
            let mass = finite(self.transport.density_kg_m3() * cut.volume_m3())?;
            if cut.volume_m3() > 0.0 && mass <= 0.0 {
                return Err(HydrostaticError::Invalid(
                    "displaced shared-tank mass is unrepresentable",
                ));
            }
            loads.push(HydrostaticBuoyancy::from_mass_moments(
                cut,
                mass,
                cut.centroid_body_m(),
                self.gravity,
                body.frame,
                state,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(FiniteTankBodiesLoad {
            surface,
            bodies: loads,
        })
    }

    pub(in crate::hydrostatic) fn prepare_bodies<'b>(
        &self,
        bodies: &[FiniteTankBody<'b>],
        states: &[RigidBodyState],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedTankBodies<'b>, HydrostaticError> {
        poll(cancelled)?;
        validate_bodies(bodies, states)?;
        let mut prepared = Vec::with_capacity(bodies.len());
        let mut whole_volume = 0.0;
        let mut lower = Vec3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut upper = Vec3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for (body, &state) in bodies.iter().zip(states) {
            poll(cancelled)?;
            let envelope = self.prepare_displacement(
                body.envelope,
                body.frame.world_from_body(state)?,
                cancelled,
            )?;
            whole_volume = finite(whole_volume + body.envelope.volume_m3())?;
            lower = Vec3::new(
                lower.x.min(envelope.lower.x),
                lower.y.min(envelope.lower.y),
                lower.z.min(envelope.lower.z),
            );
            upper = Vec3::new(
                upper.x.max(envelope.upper.x),
                upper.y.max(envelope.upper.y),
                upper.z.max(envelope.upper.z),
            );
            prepared.push(envelope);
        }
        self.check_open_area(lower, upper)?;
        Ok(PreparedTankBodies {
            envelopes: prepared,
            whole_volume_m3: whole_volume,
            cuts: vec![DisplacedVolume::dry(); bodies.len()],
        })
    }

    pub(in crate::hydrostatic) fn fill_prepared_bodies(
        &self,
        prepared: &mut PreparedTankBodies<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FiniteTankSurface, HydrostaticError> {
        self.solve_surface(prepared.whole_volume_m3, cancelled, |height, cancel| {
            let mut volume = 0.0;
            let mut depth = 0.0;
            for (envelope, cut) in prepared.envelopes.iter().zip(&mut prepared.cuts) {
                poll(cancel)?;
                *cut = envelope.at_height(height, cancel)?;
                volume = finite(volume + cut.volume_m3())?;
                depth = finite(depth + cut.depth_integral_m4)?;
            }
            Ok((volume, depth))
        })
    }
}

fn validate_bodies(
    bodies: &[FiniteTankBody<'_>],
    states: &[RigidBodyState],
) -> Result<(), HydrostaticError> {
    if bodies.is_empty() || bodies.len() > MAX_TANK_BODIES || bodies.len() != states.len() {
        return Err(HydrostaticError::Invalid(
            "shared tank requires 1..=64 matching bodies and states",
        ));
    }
    let mut cells = 0;
    for body in bodies {
        cells += body.envelope.cells().len();
        if cells > MAX_DISPLACEMENT_CELLS {
            return Err(HydrostaticError::Invalid(
                "shared tank total displacement-cell budget exceeded",
            ));
        }
    }
    Ok(())
}
