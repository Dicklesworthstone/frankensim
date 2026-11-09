//! Nonuniform face-quadrature motion through the original wall reflection.

use super::super::pressure::{GasPressure3, GasPressureField3};
use super::super::transport::MassTransport3;
use super::*;
use fs_math::dd::Dd;

/// Tangential material velocity on every cell face of one flooded plane.
///
/// With `axis = wall.axis()`, the two components are along
/// `a = (axis + 1) % 3` and `b = (axis + 2) % 3`. Entries are in
/// `index_b * domain_cells[a] + index_a` order. Each entry is the velocity at the center
/// of that Cartesian wall face, held over one native tick. Every crossing
/// link arriving in its adjacent liquid cell uses that SAME value. Thus
/// opposite tangential diagonals cancel their added mass locally, even for
/// in-plane stretching or nonuniform shear.
///
/// This face-quadrature discretization differs from evaluating continuous
/// rigid motion at each link midpoint. Use `PlanarRigidWallMotion3` for that
/// existing rule. Normal displacement, cover/vacate, finite plate edges,
/// moving contact lines and a solid constitutive law are not supplied here.
#[derive(Clone, Copy, Debug)]
pub struct PlanarWallVelocityField3<'a> {
    /// Fixed, completely flooded domain plane.
    pub wall: PlanarWall3,
    /// Common origin for the actual link angular impulses [lattice length].
    pub moment_origin: [f64; 3],
    /// Tangential face velocities [lattice length/tick], in the order above.
    pub tangential_velocity: &'a [[f64; 2]],
}

/// Explicit storage/work ceilings and actual population mass acceptance.
#[derive(Clone, Copy, Debug)]
pub struct PlanarWallFieldControls3 {
    /// Complete domain cell ceiling, including inactive gas and solid cells.
    /// Checked before the original fluid clone or full endpoint scan.
    pub maximum_domain_cells: usize,
    /// Face-entry ceiling, bounding all added storage and boundary work.
    pub maximum_wall_faces: usize,
    /// Native liquid ledger mass change divided by initial liquid mass.
    pub maximum_relative_mass_defect: f64,
    /// Absolute sum of actual incoming-minus-outgoing wall populations on
    /// each face, divided by that face's positive outgoing population sum.
    /// This is measured after floating-point reflection, not assumed zero.
    pub maximum_relative_face_mass_defect: f64,
}

impl PlanarWallFieldControls3 {
    pub(crate) fn admit_shape(
        self,
        wall: PlanarWall3,
        cells: [usize; 3],
        entries: usize,
    ) -> Result<usize, FreeSurfaceError3> {
        if cells.contains(&0)
            || self.maximum_domain_cells == 0
            || self.maximum_wall_faces == 0
            || !self.maximum_relative_mass_defect.is_finite()
            || self.maximum_relative_mass_defect < 0.0
            || !self.maximum_relative_face_mass_defect.is_finite()
            || self.maximum_relative_face_mass_defect < 0.0
        {
            return Err(refused(
                None,
                "wall field requires finite nonnegative mass bounds and positive cell ceilings",
            ));
        }
        let count = cells.into_iter().try_fold(1usize, usize::checked_mul);
        if count.is_none_or(|n| n > self.maximum_domain_cells) {
            return Err(refused(
                None,
                "wall field exceeds its complete-domain cell ceiling",
            ));
        }
        let a = (wall.axis() + 1) % 3;
        let b = (wall.axis() + 2) % 3;
        let faces = cells[a]
            .checked_mul(cells[b])
            .filter(|n| *n <= self.maximum_wall_faces && *n == entries)
            .ok_or_else(|| refused(None, "wall field shape or face-storage ceiling is invalid"))?;
        Ok(faces)
    }
}

/// Actual reciprocal exchange on one Cartesian face, in lattice units.
///
/// Angular impulse uses the actual link midpoints. Transferring only the
/// force at `center` generally loses a transverse moment: also transfer the
/// reported couple at that face. Its angular velocity is zero in this
/// piecewise translational face kinematics, so it adds no wall work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarWallFaceExchange3 {
    /// Adjacent liquid cell; face order matches the supplied velocity field.
    pub liquid_cell: [usize; 3],
    /// Center of the Cartesian wall face [lattice length].
    pub center: [f64; 3],
    /// Actual imposed tangential Cartesian velocity.
    pub velocity: [f64; 3],
    /// Actual fluid impulse and its reciprocal wall impulse.
    pub fluid_impulse: [f64; 3],
    pub wall_impulse: [f64; 3],
    /// Actual link moments about the requested common origin.
    pub fluid_angular_impulse: [f64; 3],
    pub wall_angular_impulse: [f64; 3],
    /// Residual couple about the face center, retained separately from force.
    pub fluid_couple_at_face: [f64; 3],
    pub wall_couple_at_face: [f64; 3],
    /// Density-one numerical pressure impulse and its same-origin moment.
    pub reference_pressure_fluid_impulse: [f64; 3],
    pub reference_pressure_fluid_angular_impulse: [f64; 3],
    /// Diagonal derivative of the two tangential fluid impulses with respect
    /// to their face velocities. Off-diagonal and cross-face entries are zero
    /// for this original single-tick D3Q19 reflection from the SAME input.
    /// A changed field still requires native replay and endpoint admission.
    pub fluid_tangential_response: [f64; 2],
    /// Local prescribed velocity dotted with actual fluid/wall impulse.
    pub fluid_work: f64,
    pub wall_work: f64,
    /// Actual reflected incoming-minus-outgoing mass and relative defect.
    pub reflected_mass_defect: f64,
    pub relative_mass_defect: f64,
}

/// One native VOF candidate and the resolved deforming-wall exchange.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanarWallFieldReport3 {
    pub wall: PlanarWall3,
    pub moment_origin: [f64; 3],
    /// One entry per face in exactly the supplied field order.
    pub faces: Vec<PlanarWallFaceExchange3>,
    /// Sums of actual link exchanges, accumulated with double-double sums.
    pub fluid_impulse: [f64; 3],
    pub wall_impulse: [f64; 3],
    pub fluid_angular_impulse: [f64; 3],
    pub wall_angular_impulse: [f64; 3],
    pub reference_pressure_fluid_impulse: [f64; 3],
    pub reference_pressure_fluid_angular_impulse: [f64; 3],
    /// Mechanical boundary work only; no viscous heat or whole-fluid energy
    /// conservation is inferred from this isothermal population operator.
    pub fluid_work: f64,
    pub wall_work: f64,
    pub wetted_area: f64,
    pub mass_before: f64,
    pub mass_after: f64,
    pub relative_mass_defect: f64,
    pub maximum_relative_face_mass_defect: f64,
    pub maximum_speed: f64,
    pub maximum_wall_speed: f64,
}

struct FaceAccumulator {
    liquid_cell: [usize; 3],
    center: [f64; 3],
    velocity: [f64; 3],
    impulse: [Dd; 3],
    angular: [Dd; 3],
    couple: [Dd; 3],
    response: [Dd; 2],
    work: Dd,
    mass_defect: Dd,
    outgoing_mass: Dd,
}

pub(super) struct WallFieldLink3 {
    index: usize,
    arm: [f64; 3],
    local_arm: [f64; 3],
    pub(super) velocity: [f64; 3],
}

pub(super) struct WallFieldExchange3 {
    wall: PlanarWall3,
    cells: [usize; 3],
    origin: [f64; 3],
    faces: Vec<FaceAccumulator>,
    maximum_speed: f64,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

impl WallFieldExchange3 {
    fn new(
        motion: PlanarWallVelocityField3<'_>,
        cells: [usize; 3],
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<Self, FreeSurfaceError3> {
        let count = controls.admit_shape(motion.wall, cells, motion.tangential_velocity.len())?;
        if !motion.moment_origin.into_iter().all(f64::is_finite) {
            return Err(refused(None, "wall field moment origin must be finite"));
        }
        let axis = motion.wall.axis();
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let mut faces = Vec::new();
        faces
            .try_reserve_exact(count)
            .map_err(|_| refused(None, "wall field face allocation refused"))?;
        let mut maximum_speed = 0.0_f64;
        for (index, &tangent) in motion.tangential_velocity.iter().enumerate() {
            if index % 64 == 0 {
                check_cancel(gate)?;
            }
            let local_speed = tangent[0].hypot(tangent[1]);
            if !local_speed.is_finite() || local_speed / CS2_3.sqrt() >= MACH_LIMIT {
                return Err(refused(
                    None,
                    "wall field local velocity exceeds the lattice Mach domain",
                ));
            }
            maximum_speed = maximum_speed.max(local_speed);
            let mut liquid_cell = [0usize; 3];
            liquid_cell[axis] = if motion.wall.upper() {
                cells[axis] - 1
            } else {
                0
            };
            liquid_cell[a] = index % cells[a];
            liquid_cell[b] = index / cells[a];
            let mut center = liquid_cell.map(|n| n as f64 + 0.5);
            center[axis] = if motion.wall.upper() {
                cells[axis] as f64
            } else {
                0.0
            };
            let mut velocity = [0.0; 3];
            velocity[a] = tangent[0];
            velocity[b] = tangent[1];
            faces.push(FaceAccumulator {
                liquid_cell,
                center,
                velocity,
                impulse: [Dd::ZERO; 3],
                angular: [Dd::ZERO; 3],
                couple: [Dd::ZERO; 3],
                response: [Dd::ZERO; 2],
                work: Dd::ZERO,
                mass_defect: Dd::ZERO,
                outgoing_mass: Dd::ZERO,
            });
        }
        Ok(Self {
            wall: motion.wall,
            cells,
            origin: motion.moment_origin,
            faces,
            maximum_speed,
        })
    }

    pub(super) fn link(
        &self,
        source: [i64; 3],
        q: usize,
    ) -> Result<WallFieldLink3, FreeSurfaceError3> {
        let e = [E3[q].0 as i64, E3[q].1 as i64, E3[q].2 as i64];
        let target: [i64; 3] = core::array::from_fn(|d| source[d] + e[d]);
        if (0..3).any(|d| target[d] < 0 || target[d] >= self.cells[d] as i64) {
            return Err(refused(None, "wall field link has no adjacent liquid cell"));
        }
        let axis = self.wall.axis();
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let index = target[b] as usize * self.cells[a] + target[a] as usize;
        let face = &self.faces[index];
        if face.liquid_cell != target.map(|n| n as usize) {
            return Err(refused(
                None,
                "wall field link does not arrive on its declared plane",
            ));
        }
        let point: [f64; 3] = core::array::from_fn(|d| source[d] as f64 + 0.5 + 0.5 * e[d] as f64);
        let arm = core::array::from_fn(|d| point[d] - self.origin[d]);
        let local_arm = core::array::from_fn(|d| point[d] - face.center[d]);
        if !arm.into_iter().all(f64::is_finite) {
            return Err(refused(None, "wall field moment arm is unrepresentable"));
        }
        Ok(WallFieldLink3 {
            index,
            arm,
            local_arm,
            velocity: face.velocity,
        })
    }

    pub(super) fn record(
        &mut self,
        link: WallFieldLink3,
        direction: [f64; 3],
        incoming: f64,
        outgoing: f64,
        derivative: f64,
    ) -> Result<(), FreeSurfaceError3> {
        let face = &mut self.faces[link.index];
        // Both already-rounded native populations belong to the actual impulse.
        let strength = Dd::from_f64(incoming) + Dd::from_f64(outgoing);
        let moment = cross(link.arm, direction);
        let local_moment = cross(link.local_arm, direction);
        for d in 0..3 {
            let impulse = strength * Dd::from_f64(direction[d]);
            face.impulse[d] = face.impulse[d] + impulse;
            face.angular[d] = face.angular[d] + strength * Dd::from_f64(moment[d]);
            face.couple[d] = face.couple[d] + strength * Dd::from_f64(local_moment[d]);
            face.work = face.work + Dd::from_f64(link.velocity[d]) * impulse;
        }
        for j in 0..2 {
            let component = direction[(self.wall.axis() + 1 + j) % 3];
            face.response[j] =
                face.response[j] + Dd::from_f64(derivative) * Dd::from_f64(component * component);
        }
        face.mass_defect = face.mass_defect + (Dd::from_f64(incoming) - Dd::from_f64(outgoing));
        face.outgoing_mass = face.outgoing_mass + Dd::from_f64(outgoing);
        if !face
            .impulse
            .into_iter()
            .chain(face.angular)
            .chain(face.couple)
            .chain(face.response)
            .chain([face.work, face.mass_defect, face.outgoing_mass])
            .all(|x| x.to_f64().is_finite())
        {
            return Err(refused(
                Some(face.liquid_cell),
                "wall field exchange or response is unrepresentable",
            ));
        }
        Ok(())
    }

    fn finish(
        self,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<PlanarWallFieldReport3, FreeSurfaceError3> {
        let mut reports = Vec::new();
        reports
            .try_reserve_exact(self.faces.len())
            .map_err(|_| refused(None, "wall field report allocation refused"))?;
        let mut total_impulse = [Dd::ZERO; 3];
        let mut total_angular = [Dd::ZERO; 3];
        let mut reference_angular = [Dd::ZERO; 3];
        let mut total_work = Dd::ZERO;
        let mut maximum_mass_defect = 0.0_f64;
        let axis = self.wall.axis();
        let mut reference = [0.0; 3];
        reference[axis] = if self.wall.upper() { -CS2_3 } else { CS2_3 };
        let area = self.faces.len() as f64;
        for (index, face) in self.faces.into_iter().enumerate() {
            if index % 64 == 0 {
                check_cancel(gate)?;
            }
            let outgoing = face.outgoing_mass.to_f64();
            let mass_defect = face.mass_defect.to_f64();
            let relative_mass_defect = (mass_defect / outgoing).abs();
            if outgoing <= 0.0
                || !relative_mass_defect.is_finite()
                || relative_mass_defect > controls.maximum_relative_face_mass_defect
            {
                return Err(refused(
                    Some(face.liquid_cell),
                    "wall field violates its local reflected-mass bound",
                ));
            }
            maximum_mass_defect = maximum_mass_defect.max(relative_mass_defect);
            let reference_moment = cross(
                core::array::from_fn(|d| face.center[d] - self.origin[d]),
                reference,
            );
            if !reference_moment.into_iter().all(f64::is_finite) {
                return Err(refused(
                    Some(face.liquid_cell),
                    "wall field reference-pressure moment is unrepresentable",
                ));
            }
            for d in 0..3 {
                total_impulse[d] = total_impulse[d] + face.impulse[d];
                total_angular[d] = total_angular[d] + face.angular[d];
                reference_angular[d] = reference_angular[d] + Dd::from_f64(reference_moment[d]);
            }
            total_work = total_work + face.work;
            let fluid_impulse = face.impulse.map(Dd::to_f64);
            let fluid_angular_impulse = face.angular.map(Dd::to_f64);
            let fluid_couple_at_face = face.couple.map(Dd::to_f64);
            let fluid_work = face.work.to_f64();
            reports.push(PlanarWallFaceExchange3 {
                liquid_cell: face.liquid_cell,
                center: face.center,
                velocity: face.velocity,
                fluid_impulse,
                wall_impulse: fluid_impulse.map(|v| -v),
                fluid_angular_impulse,
                wall_angular_impulse: fluid_angular_impulse.map(|v| -v),
                fluid_couple_at_face,
                wall_couple_at_face: fluid_couple_at_face.map(|v| -v),
                reference_pressure_fluid_impulse: reference,
                reference_pressure_fluid_angular_impulse: reference_moment,
                fluid_tangential_response: face.response.map(Dd::to_f64),
                fluid_work,
                wall_work: -fluid_work,
                reflected_mass_defect: mass_defect,
                relative_mass_defect,
            });
        }
        if !total_impulse
            .into_iter()
            .chain(total_angular)
            .chain(reference_angular)
            .chain([total_work])
            .all(|x| x.to_f64().is_finite())
        {
            return Err(refused(
                None,
                "wall field total exchange is unrepresentable",
            ));
        }
        let fluid_impulse = total_impulse.map(Dd::to_f64);
        let fluid_angular_impulse = total_angular.map(Dd::to_f64);
        let fluid_work = total_work.to_f64();
        Ok(PlanarWallFieldReport3 {
            wall: self.wall,
            moment_origin: self.origin,
            faces: reports,
            fluid_impulse,
            wall_impulse: fluid_impulse.map(|v| -v),
            fluid_angular_impulse,
            wall_angular_impulse: fluid_angular_impulse.map(|v| -v),
            reference_pressure_fluid_impulse: reference.map(|v| v * area),
            reference_pressure_fluid_angular_impulse: reference_angular.map(Dd::to_f64),
            fluid_work,
            wall_work: -fluid_work,
            wetted_area: area,
            mass_before: 0.0,
            mass_after: 0.0,
            relative_mass_defect: 0.0,
            maximum_relative_face_mass_defect: maximum_mass_defect,
            maximum_speed: 0.0,
            maximum_wall_speed: self.maximum_speed,
        })
    }
}

impl FreeSurface3 {
    /// Propose one original native step with a nonuniform tangential wall field.
    ///
    /// Local face velocities and both liquid endpoints retain the native Mach,
    /// population-positivity and completely flooded-plane requirements. Added
    /// storage is linear in wall faces, and one original fluid candidate is
    /// cloned after count and field admission. The same reflected populations
    /// supply every mass, momentum, couple, response and work account.
    /// Refusal or cancellation returns no changed liquid.
    pub fn planar_wall_velocity_field_candidate(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarWallFieldReport3), FreeSurfaceError3> {
        self.planar_wall_field_candidate_with_relaxation(motion, controls, None, gate)
    }

    // The source-rheology owner supplies its accepted cellwise viscosity to
    // this same candidate. Isothermal flow needs no thermal mass-route trace.
    pub(crate) fn planar_wall_field_candidate_with_relaxation(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        tau: Option<&BTreeMap<(u64, usize), f64>>,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarWallFieldReport3), FreeSurfaceError3> {
        self.planar_wall_field_candidate_with_relaxation_and_bulk(motion, controls, tau, None, gate)
    }

    // A changed constitutive/SGS input may also change source-owned bulk
    // viscosity. Install that same source epoch before collision on the ONE
    // native candidate, without cloning an already-private flow a second time.
    pub(crate) fn planar_wall_field_candidate_with_relaxation_and_bulk(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        tau: Option<&BTreeMap<(u64, usize), f64>>,
        bulk: Option<BTreeMap<(u64, usize), f64>>,
        gate: &CancelGate,
    ) -> Result<(Self, PlanarWallFieldReport3), FreeSurfaceError3> {
        let (mut exchange, mass_before, before_speed) =
            self.prepare_planar_wall_velocity_field(motion, controls, gate)?;
        let mut candidate = self.clone();
        if let Some(bulk) = bulk {
            candidate.rebind_bulk_relaxations(bulk, gate)?;
        }
        candidate.step_inner_with_wall(gate, None, tau, None, Some(&mut exchange))?;
        let report = candidate.finish_planar_wall_velocity_field(
            motion,
            controls,
            exchange,
            mass_before,
            before_speed,
            gate,
        )?;
        Ok((candidate, report))
    }

    /// Thermal owner supplies the one private candidate and records the exact
    /// same collision/streaming/conversion routes as its stationary-wall path.
    /// Prescribed gas traction reconstructs only actual gas links during this
    /// same tick; it is never included in the reflected solid-wall wrench.
    /// No second flow solve or velocity-reconstructed companion flux occurs.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn step_with_planar_wall_field_transport(
        &mut self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        tau: &BTreeMap<(u64, usize), f64>,
        forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
        maximum_routes: usize,
        gas_pressure: Option<(f64, GasPressureField3<'_>)>,
        gate: &CancelGate,
    ) -> Result<(MassTransport3, PlanarWallFieldReport3), FreeSurfaceError3> {
        check_cancel(gate)?;
        let pressure = gas_pressure
            .map(|(bound, field)| GasPressure3::new(bound, field))
            .transpose()?;
        if let Some(forces) = forces {
            self.rebind_force_field(forces, gate)?;
        }
        let (mut exchange, mass_before, before_speed) =
            self.prepare_planar_wall_velocity_field(motion, controls, gate)?;
        let mut trace = MassTransport3::new(self.mass_snapshot(), maximum_routes);
        self.step_inner_with_wall(
            gate,
            Some(&mut trace),
            Some(tau),
            pressure.as_ref(),
            Some(&mut exchange),
        )?;
        trace.after = self.mass_snapshot();
        let report = self.finish_planar_wall_velocity_field(
            motion,
            controls,
            exchange,
            mass_before,
            before_speed,
            gate,
        )?;
        Ok((trace, report))
    }

    fn admit_planar_wall_velocity_field_state(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        gate: &CancelGate,
    ) -> Result<f64, FreeSurfaceError3> {
        let cells = self.domain_cells();
        let axis = motion.wall.axis();
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        let admit = |p: [usize; 3], liquid: [f64; 3]| {
            let velocity = motion.tangential_velocity[p[b] * cells[a] + p[a]];
            let mut relative = liquid;
            relative[a] -= velocity[0];
            relative[b] -= velocity[1];
            if !speed(relative).is_finite() || speed(relative) / CS2_3.sqrt() >= MACH_LIMIT {
                return Err(refused(
                    Some(p),
                    "local liquid/wall field velocity exceeds the lattice Mach domain",
                ));
            }
            Ok(())
        };
        self.admit_planar_wall_state_with(motion.wall, gate, admit)
            .map(|(_, maximum_speed)| maximum_speed)
    }

    fn prepare_planar_wall_velocity_field(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<(PlanarWallStep3, f64, f64), FreeSurfaceError3> {
        check_cancel(gate)?;
        let cells = self.domain_cells();
        let field = WallFieldExchange3::new(motion, cells, controls, gate)?;
        let before_speed = self.admit_planar_wall_velocity_field_state(motion, gate)?;
        let mass_before = self.ledger_mass();
        if !mass_before.is_finite() || mass_before <= 0.0 {
            return Err(refused(
                None,
                "wall field requires finite positive liquid mass",
            ));
        }
        let exchange = PlanarWallStep3 {
            wall: motion.wall,
            extent: cells[motion.wall.axis()] as i64,
            velocity: [0.0; 3],
            impulse: [0.0; 3],
            rigid: None,
            field: Some(field),
        };
        Ok((exchange, mass_before, before_speed))
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_planar_wall_velocity_field(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        exchange: PlanarWallStep3,
        mass_before: f64,
        before_speed: f64,
        gate: &CancelGate,
    ) -> Result<PlanarWallFieldReport3, FreeSurfaceError3> {
        let mut report = exchange
            .field
            .expect("field retained by original streaming")
            .finish(controls, gate)?;
        report.mass_before = mass_before;
        report.maximum_speed = before_speed;
        self.refresh_planar_wall_field_endpoint(motion, controls, &mut report, gate)?;
        Ok(report)
    }

    /// A coupled owner may refresh forces/geometry after the native tick. Its
    /// final candidate must still satisfy the same wall and ledger admission.
    pub(crate) fn refresh_planar_wall_field_endpoint(
        &self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        report: &mut PlanarWallFieldReport3,
        gate: &CancelGate,
    ) -> Result<(), FreeSurfaceError3> {
        let after_speed = self.admit_planar_wall_velocity_field_state(motion, gate)?;
        let mass_after = self.ledger_mass();
        let mass_before = report.mass_before;
        let relative_mass_defect = ((mass_after - mass_before) / mass_before).abs();
        if !mass_after.is_finite()
            || !relative_mass_defect.is_finite()
            || relative_mass_defect > controls.maximum_relative_mass_defect
        {
            return Err(refused(
                None,
                "wall field exceeds its native liquid mass-conservation tolerance",
            ));
        }
        report.mass_after = mass_after;
        report.relative_mass_defect = relative_mass_defect;
        report.maximum_speed = report.maximum_speed.max(after_speed);
        check_cancel(gate)?;
        Ok(())
    }

    /// Atomically publish the same fully admitted spatial wall-field candidate.
    pub fn step_with_planar_wall_velocity_field(
        &mut self,
        motion: PlanarWallVelocityField3<'_>,
        controls: PlanarWallFieldControls3,
        gate: &CancelGate,
    ) -> Result<PlanarWallFieldReport3, FreeSurfaceError3> {
        let (candidate, report) =
            self.planar_wall_velocity_field_candidate(motion, controls, gate)?;
        check_cancel(gate)?;
        *self = candidate;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
