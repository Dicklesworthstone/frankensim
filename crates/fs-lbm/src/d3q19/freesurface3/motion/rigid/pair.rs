//! Two opposing rigid planes acting on one original liquid timestep.

use super::*;

/// One native fluid interval and both original rigid-wall exchanges.
///
/// The two entries preserve caller order and their own moment origins. Liquid
/// mass is counted once. Net impulse/work are sums of actual reflected-link
/// contributions, not endpoint pressures or a second elapsed fluid timestep.
/// The native momentum lattice has no total-energy equation. A thermal owner
/// transports its source enthalpy separately; boundary work alone is not a
/// claim of liquid total-energy conservation or resolved viscous heating.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarRigidWallPairReport3 {
    /// Original responses, reciprocal wrenches and work for both planes.
    pub walls: [PlanarRigidWallReport3; 2],
    /// Net actual raw population impulse delivered to liquid by these planes.
    pub fluid_impulse: [f64; 3],
    /// Exactly reciprocal net impulse on the pair of walls.
    pub wall_impulse: [f64; 3],
    /// Sum of both original boundary mechanical work contributions.
    pub fluid_work: f64,
    /// Exactly reciprocal net actuator work.
    pub wall_work: f64,
    /// The single liquid inventory before and after this candidate.
    pub mass_before: f64,
    pub mass_after: f64,
    /// Absolute relative native liquid mass defect.
    pub relative_mass_defect: f64,
    /// Maximum physical wet-cell speed at either endpoint.
    pub maximum_speed: f64,
    /// Maximum prescribed local wall speed over either finite plane.
    pub maximum_wall_speed: f64,
}

impl PlanarRigidWallPairReport3 {
    fn from_walls(walls: [PlanarRigidWallReport3; 2]) -> Result<Self, FreeSurfaceError3> {
        let fluid_impulse = core::array::from_fn(|d| {
            (Dd::from_f64(walls[0].fluid_impulse[d]) + Dd::from_f64(walls[1].fluid_impulse[d]))
                .to_f64()
        });
        let fluid_work =
            (Dd::from_f64(walls[0].fluid_work) + Dd::from_f64(walls[1].fluid_work)).to_f64();
        if !fluid_impulse
            .into_iter()
            .chain([fluid_work])
            .all(f64::is_finite)
        {
            return Err(refused(
                None,
                "rigid wall-pair impulse or work is unrepresentable",
            ));
        }
        Ok(Self {
            fluid_impulse,
            wall_impulse: fluid_impulse.map(|value| -value),
            fluid_work,
            wall_work: -fluid_work,
            mass_before: walls[0].mass_before,
            mass_after: walls[0].mass_after,
            relative_mass_defect: walls[0].relative_mass_defect,
            maximum_speed: walls[0].maximum_speed.max(walls[1].maximum_speed),
            maximum_wall_speed: walls[0].maximum_wall_speed.max(walls[1].maximum_wall_speed),
            walls,
        })
    }

    /// Sum both actual liquid angular impulses about a common fixed origin.
    ///
    /// Each original wall moment is shifted with its own linear impulse before
    /// summation. Adding moments at different origins directly would discard
    /// the orbital torque of separated planes. Raw numerical pressure remains
    /// in each original wrench; opposite complete-plane reference wrenches
    /// cancel only after this common-origin shift.
    pub fn fluid_angular_impulse_at(
        &self,
        origin: [f64; 3],
    ) -> Result<[f64; 3], FreeSurfaceError3> {
        if !origin.into_iter().all(f64::is_finite) {
            return Err(refused(None, "wall-pair moment origin must be finite"));
        }
        let mut moment = [Dd::ZERO; 3];
        for wall in self.walls {
            let arm = core::array::from_fn::<_, 3, _>(|d| {
                Dd::from_f64(wall.motion.origin[d]) - Dd::from_f64(origin[d])
            });
            for d in 0..3 {
                let a = (d + 1) % 3;
                let b = (d + 2) % 3;
                moment[d] = moment[d]
                    + Dd::from_f64(wall.fluid_angular_impulse[d])
                    + arm[a] * Dd::from_f64(wall.fluid_impulse[b])
                    - arm[b] * Dd::from_f64(wall.fluid_impulse[a]);
            }
        }
        let result = moment.map(Dd::to_f64);
        if !result.into_iter().all(f64::is_finite) {
            return Err(refused(
                None,
                "wall-pair angular impulse is unrepresentable",
            ));
        }
        Ok(result)
    }

    /// Exactly reciprocal wall angular impulse at the requested common origin.
    pub fn wall_angular_impulse_at(&self, origin: [f64; 3]) -> Result<[f64; 3], FreeSurfaceError3> {
        Ok(self.fluid_angular_impulse_at(origin)?.map(|value| -value))
    }
}

impl FreeSurface3 {
    /// Propose independent rigid motion of TWO opposing flooded domain planes.
    ///
    /// Both original planes must be on the same Cartesian axis and on opposite
    /// sides. Caller order is preserved. Every crossing link belongs to at
    /// most one selected plane; adjacent planes are refused rather than assigning
    /// an arbitrary corner reflection. Each plane admits its original tangential
    /// translation and axial spin, local absolute/relative Mach limits and fully
    /// flooded endpoints. The geometric planes do not move or tilt.
    ///
    /// There is one private liquid candidate, ONE native collision/streaming/
    /// conversion step, and two original reflected-link wrench/work responses.
    /// These same-input affine response blocks have zero instantaneous cross
    /// derivative: collision sees the accepted liquid, before either boundary
    /// reflection. Propagation through the actual fluid couples the walls over
    /// subsequent steps. A finite-body owner must solve both motions from that
    /// same accepted input and replay the pair together, never advance the
    /// liquid once for each body.
    ///
    /// Failure or cancellation leaves the caller's accepted liquid unchanged.
    /// Normal wall motion, remapping, partial wetting, moving contact lines,
    /// finite plate edges, thermal feedback and solid dynamics are not added.
    pub fn planar_rigid_wall_pair_candidate(
        &self,
        motions: [PlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarRigidWallPairReport3), FreeSurfaceError3> {
        let prepared =
            self.prepare_planar_rigid_wall_pair(motions, maximum_relative_mass_defect, gate)?;
        let mut candidate = self.clone();
        let report = candidate.step_prepared_planar_rigid_wall_pair(
            motions,
            maximum_relative_mass_defect,
            prepared,
            None,
            None,
            None,
            gate,
        )?;
        Ok((candidate, report))
    }

    fn prepare_planar_rigid_wall_pair(
        &self,
        motions: [PlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<[PreparedRigidWall3; 2], FreeSurfaceError3> {
        check_cancel(gate)?;
        if motions[0].wall.axis() != motions[1].wall.axis()
            || motions[0].wall.upper() == motions[1].wall.upper()
        {
            return Err(refused(
                None,
                "rigid wall pair requires opposite planes on one axis",
            ));
        }
        Ok([
            self.prepare_planar_rigid_wall_motion(motions[0], maximum_relative_mass_defect, gate)?,
            self.prepare_planar_rigid_wall_motion(motions[1], maximum_relative_mass_defect, gate)?,
        ])
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared_planar_rigid_wall_pair(
        &mut self,
        motions: [PlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        prepared: [PreparedRigidWall3; 2],
        transport: Option<&mut MassTransport3>,
        tau: Option<&BTreeMap<(u64, usize), f64>>,
        gas_pressure: Option<&GasPressure3<'_>>,
        gate: &CancelGate,
    ) -> Result<PlanarRigidWallPairReport3, FreeSurfaceError3> {
        let [mut first, mut second] = prepared;
        let mut exchanges = [first.exchange, second.exchange];
        self.step_inner_with_boundaries(
            gate,
            transport,
            tau,
            gas_pressure,
            Some(&mut exchanges),
            None,
        )?;
        [first.exchange, second.exchange] = exchanges;
        let walls = [
            self.finish_planar_rigid_wall_motion(
                motions[0],
                maximum_relative_mass_defect,
                first,
                gate,
            )?,
            self.finish_planar_rigid_wall_motion(
                motions[1],
                maximum_relative_mass_defect,
                second,
                gate,
            )?,
        ];
        let report = PlanarRigidWallPairReport3::from_walls(walls)?;
        check_cancel(gate)?;
        Ok(report)
    }

    /// The thermal owner supplies its sole private liquid candidate. Both
    /// rigid reflections and the cellwise relaxation/force field act during
    /// one original native tick, whose actual mass routes transport enthalpy.
    /// Optional gas loading uses that tick's original pressure reconstruction,
    /// separately from either wall's measured reflection impulse and work.
    /// No population clone, second collision or substitute mass flux is made.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_with_planar_rigid_wall_pair_transport(
        &mut self,
        motions: [PlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        maximum_routes: usize,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
        gate: &CancelGate,
    ) -> Result<(MassTransport3, PlanarRigidWallPairReport3), FreeSurfaceError3> {
        check_cancel(gate)?;
        let pressure = gas_pressure
            .map(|(bound, field)| GasPressure3::new(bound, field))
            .transpose()?;
        if let Some(forces) = forces {
            self.rebind_force_field(forces, gate)?;
        }
        let prepared =
            self.prepare_planar_rigid_wall_pair(motions, maximum_relative_mass_defect, gate)?;
        let mut trace = MassTransport3::new(self.mass_snapshot(), maximum_routes);
        let report = self.step_prepared_planar_rigid_wall_pair(
            motions,
            maximum_relative_mass_defect,
            prepared,
            Some(&mut trace),
            Some(tau),
            pressure.as_ref(),
            gate,
        )?;
        trace.after = self.mass_snapshot();
        check_cancel(gate)?;
        Ok((trace, report))
    }

    /// Recheck both actual endpoints after the thermal owner's force/wetting
    /// refresh. Keep the measured reflection impulses and work unchanged; only
    /// the shared inventory and endpoint admission are refreshed. A late second
    /// wall refusal leaves this report intact for the enclosing transaction.
    pub(crate) fn refresh_planar_rigid_wall_pair_endpoint(
        &self,
        report: &mut PlanarRigidWallPairReport3,
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        let mut walls = report.walls;
        for wall in &mut walls {
            self.refresh_planar_rigid_wall_endpoint(wall, maximum_relative_mass_defect, gate)?;
        }
        let refreshed = PlanarRigidWallPairReport3::from_walls(walls)?;
        check_cancel(gate)?;
        *report = refreshed;
        Ok(())
    }

    /// Publish the same fully admitted two-plane candidate in one transaction.
    pub fn step_with_planar_rigid_wall_pair(
        &mut self,
        motions: [PlanarRigidWallMotion3; 2],
        maximum_relative_mass_defect: f64,
        gate: &CancelGate,
    ) -> Result<PlanarRigidWallPairReport3, FreeSurfaceError3> {
        let (candidate, report) =
            self.planar_rigid_wall_pair_candidate(motions, maximum_relative_mass_defect, gate)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }
}
