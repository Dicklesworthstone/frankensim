//! A finite wall with two free tangent translations and one principal spin.
//!
//! The original material liquid measures a symmetric three-mode impedance in
//! one pilot. A mass-whitened 3-by-3 solve finds the common midpoint twist, and
//! one final native replay supplies the actual load to the original MBD owner.
//! No scalar coordinate is solved independently of the other two.
//!
//! The coupling origin is the free half-drift COM. For tangent forces and a
//! principal axial torque, the native MBD endpoint is also the two free
//! half-drifts around that impulse configuration. This choice retains orbital
//! as well as spin momentum and midpoint work. It does not establish the
//! temporal order of the complete liquid coupling. The flooded domain plane
//! stays fixed: finite plate edges, tilt, normal motion and heat are absent.

use super::{array, norm, normal_axis, vector};
use fs_exec::CancelGate;
use fs_la::factor::cholesky;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialPlanarRigidWallMotion3, MaterialPlanarRigidWallReport3,
};
use fs_math::dd::Dd;
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3,
    Wrench,
};

/// Physical and complete-grid work bounds for one coupled liquid timestep.
#[derive(Clone, Copy, Debug)]
pub struct PlanarWallControls {
    /// Complete domain cells across the pilot and accepted replay.
    pub maximum_trial_cells: usize,
    /// Original native mass-conservation allowance on both liquid candidates.
    pub maximum_relative_mass_defect: f64,
    /// Bound on the sum of the two free half-drift distances [m].
    pub maximum_translation_m: f64,
    /// Bound on the sum of the two free half-drift angular distances [rad].
    pub maximum_rotation_rad: f64,
    /// Maximum absolute mechanical wall work in the final liquid replay [J].
    pub maximum_absolute_wall_work_j: f64,
    /// Actual linear-impulse and translation-response allowance [kg m/s].
    pub momentum_absolute_kg_m_s: f64,
    /// Actual spin/orbital impulse and spin-response allowance [kg m²/s].
    pub angular_momentum_absolute_kg_m2_s: f64,
    /// Actual versus imposed midpoint translation allowance [m/s].
    pub midpoint_absolute_m_s: f64,
    /// Actual versus imposed midpoint angular-velocity allowance [rad/s].
    pub midpoint_absolute_rad_s: f64,
    /// Rigid kinetic/boundary/external work defect allowance [J].
    pub work_absolute_j: f64,
    /// Relative tolerance in [0,1), with no unit-sized floor or scaling by
    /// numerical normal-pressure forces that the guide cancels.
    pub relative_tolerance: f64,
}

/// Refusal leaves the original liquid, native body and physical clock intact.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarWallError {
    /// Unsupported physical state, limit or representability boundary.
    Invalid(&'static str),
    /// Original liquid source, scaling, speed or population refusal.
    Fluid(MaterialLbmError),
    /// Original native rigid state, inertia or integration refusal.
    Mechanics(DynamicsError),
    /// Pilot and final replay exceed the declared complete-grid work ceiling.
    WorkBudget {
        /// Two complete native liquid domain traversals.
        required_trial_cells: usize,
        /// Caller-supplied ceiling.
        maximum: usize,
    },
    /// Actual represented endpoints or the final affine replay failed admission.
    Balance {
        /// Actual represented linear increment minus supplied impulses.
        momentum_residual_kg_m_s: Vec3,
        /// Actual spatial spin increment minus supplied couples.
        spin_residual_kg_m2_s: Vec3,
        /// Full orbital-plus-spin increment about the fixed world origin.
        angular_momentum_residual_kg_m2_s: Vec3,
        /// Actual rigid midpoint minus imposed tangent velocity.
        midpoint_residual_m_s: Vec3,
        /// Actual rigid midpoint minus imposed axial angular velocity.
        midpoint_residual_rad_s: Vec3,
        /// Two tangent impulse residuals, then axial angular impulse residual.
        response_residual: [f64; 3],
        /// Kinetic increment minus mechanical boundary and external work.
        work_residual_j: f64,
    },
    /// Cancellation observed before any complete candidate escaped.
    Cancelled,
}

impl core::fmt::Display for PlanarWallError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite planar liquid wall refused: {self:?}")
    }
}
impl std::error::Error for PlanarWallError {}
impl From<MaterialLbmError> for PlanarWallError {
    fn from(value: MaterialLbmError) -> Self {
        Self::Fluid(value)
    }
}
impl From<DynamicsError> for PlanarWallError {
    fn from(value: DynamicsError) -> Self {
        Self::Mechanics(value)
    }
}

/// Original source-bound liquid and native finite wall retained for continuation.
#[derive(Clone)]
pub struct FreeSurfacePlanarWall<'a> {
    liquid: MaterialFreeSurface3<'a>,
    properties: MassProperties,
    body: RigidBodyState,
    wall: PlanarWall3,
    axis_body: Vec3,
    inertia_kg_m2: f64,
    domain_cells: usize,
    elapsed_time_s: f64,
    accepted_steps: u64,
}

/// Jointly admitted wall, liquid, guide reactions and physical work.
pub struct PlanarWallStep<'a> {
    /// Both accepted original owners; continue by stepping this motion.
    pub motion: FreeSurfacePlanarWall<'a>,
    /// Final original reflection. Its moments use the free half-drift COM.
    pub wall: MaterialPlanarRigidWallReport3,
    /// Native MBD endpoint; includes both the actual COM drift and attitude.
    pub mechanics: StepReceipt,
    /// Configuration used by both fluid candidates and every coupled impulse.
    pub coupling_origin_world_m: Vec3,
    /// Applied tangent force integrated over the native liquid timestep [kg m/s].
    pub external_linear_impulse_world: Vec3,
    /// Applied axial couple, about the coupling origin [kg m²/s].
    pub external_couple_impulse_world: Vec3,
    /// Normal guide impulse cancelling the raw wall force [kg m/s].
    /// Includes numerical pressure; it is not an absolute physical support load.
    pub raw_guide_impulse_world: Vec3,
    /// Transverse guide couple about the coupling origin [kg m²/s].
    pub raw_guide_angular_impulse_world: Vec3,
    /// Physical normal support impulse about the declared pressure reference.
    pub pressure_perturbation_guide_impulse_world: Vec3,
    /// Physical transverse guide couple about that same reference.
    pub pressure_perturbation_guide_angular_impulse_world: Vec3,
    /// Normal support impulse including an explicitly supplied pressure datum.
    pub absolute_guide_impulse_world: Option<Vec3>,
    /// Transverse support couple including that same absolute pressure datum.
    pub absolute_guide_angular_impulse_world: Option<Vec3>,
    /// Common tangent velocity imposed at the free half-drift COM.
    pub midpoint_velocity_m_s: Vec3,
    /// Common angular speed about the positive coordinate normal.
    pub midpoint_angular_speed_rad_s: f64,
    /// Actual native midpoint velocity minus the imposed wall velocity.
    pub midpoint_residual_m_s: Vec3,
    /// Actual native midpoint angular velocity minus imposed axial motion.
    pub midpoint_residual_rad_s: Vec3,
    /// Actual body linear-momentum increment minus every supplied impulse.
    pub momentum_residual_kg_m_s: Vec3,
    /// Actual spatial spin increment minus the wall, guide and applied couples.
    pub spin_residual_kg_m2_s: Vec3,
    /// Actual total orbital-plus-spin increment minus all corresponding
    /// impulses, about the fixed world origin, in kg m²/s.
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Actual final wall impulse minus the pilot's full affine prediction:
    /// tangent a, tangent b [kg m/s], then normal spin [kg m²/s].
    pub response_residual: [f64; 3],
    /// Actual native momentum-endpoint kinetic increment, in joules.
    pub kinetic_energy_change_j: f64,
    /// Applied force/couple impulses paired with actual native midpoint motion.
    pub external_work_j: f64,
    /// Normal guide force/transverse couple paired with the actual midpoint.
    pub raw_guide_work_j: f64,
    /// Kinetic change + liquid wall work - external and guide work [J].
    /// This numerical defect is never deposited as heat.
    pub mechanical_work_residual_j: f64,
    /// Exactly two original liquid candidates, with only one accepted tick.
    pub native_trials: usize,
    /// Sum of complete-domain cell counts in those two candidates.
    pub trial_cells: usize,
}

impl<'a> FreeSurfacePlanarWall<'a> {
    /// Bind one fixed flooded plane to the original finite rigid-body inertia.
    ///
    /// The normal is a principal inertia axis. Initial linear momentum must
    /// be tangential and spin must be axial to arithmetic roundoff; no state
    /// is projected. COM may be off the plane. The guide prevents normal
    /// translation and tilt, while both tangent translations and spin are free.
    pub fn new(
        liquid: MaterialFreeSurface3<'a>,
        properties: MassProperties,
        body: RigidBodyState,
        wall: PlanarWall3,
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PlanarWallError> {
        poll(gate)?;
        body.kinetic_energy(properties)?;
        let axis = normal_axis(wall);
        if array(body.linear_momentum_world())[axis] != 0.0 {
            return Err(invalid("planar wall requires tangential initial momentum"));
        }
        let axis_body = body
            .pose()
            .orientation()
            .rotate_world_to_body(axis_vector(axis));
        let unit_momentum = properties.angular_momentum_body_checked(axis_body)?;
        let inertia_kg_m2 = finite(axis_body.dot(unit_momentum))?;
        if inertia_kg_m2 <= 0.0
            || norm(unit_momentum.sub(axis_body.scale(inertia_kg_m2)))
                > 128.0 * f64::EPSILON * inertia_kg_m2
        {
            return Err(invalid(
                "wall normal must be a principal axis of the original inertia",
            ));
        }
        let spin = body.angular_momentum_body();
        if norm(spin.sub(axis_body.scale(axis_body.dot(spin)))) > 128.0 * f64::EPSILON * norm(spin)
        {
            return Err(invalid(
                "initial wall spin must lie on its fixed normal axis",
            ));
        }
        let domain_cells = liquid
            .lattice()
            .domain_cells()
            .into_iter()
            .try_fold(1_usize, |n, x| n.checked_mul(x))
            .ok_or_else(|| invalid("planar wall cell count overflowed"))?;
        if domain_cells == 0 || domain_cells > maximum_cells {
            return Err(invalid(
                "liquid domain exceeds the planar wall cell ceiling",
            ));
        }
        let result = Self {
            liquid,
            properties,
            body,
            wall,
            axis_body,
            inertia_kg_m2,
            domain_cells,
            elapsed_time_s: 0.0,
            accepted_steps: 0,
        };
        result.admit_endpoint(body)?;
        poll(gate)?;
        Ok(result)
    }

    /// Original accepted source-bound liquid and its population history.
    pub const fn liquid(&self) -> &MaterialFreeSurface3<'a> {
        &self.liquid
    }
    /// Original native body, including the actual COM and attitude.
    pub const fn body(&self) -> RigidBodyState {
        self.body
    }
    /// Original physical mass and inertia; no fluid-added surrogate replaces it.
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
    /// Physical time since binding; discarded fluid trials advance no time.
    pub const fn elapsed_time_s(&self) -> f64 {
        self.elapsed_time_s
    }
    /// Number of complete jointly admitted native timesteps.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Signed physical translation and normal spin of the accepted native body.
    pub fn twist(&self) -> Result<(Vec3, f64), PlanarWallError> {
        let q = self.body_twist(self.body)?;
        Ok((self.tangent(q), q[2]))
    }

    fn body_twist(&self, body: RigidBodyState) -> Result<[f64; 3], PlanarWallError> {
        let axis = normal_axis(self.wall);
        let p = array(body.linear_momentum_world());
        let omega = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        Ok([
            finite(p[(axis + 1) % 3] / self.properties.mass())?,
            finite(p[(axis + 2) % 3] / self.properties.mass())?,
            finite(omega.dot(self.axis_body))?,
        ])
    }

    fn tangent(&self, q: [f64; 3]) -> Vec3 {
        let axis = normal_axis(self.wall);
        let mut v = [0.0; 3];
        v[(axis + 1) % 3] = q[0];
        v[(axis + 2) % 3] = q[1];
        vector(v)
    }

    fn prescribed(&self, origin: Vec3, q: [f64; 3]) -> MaterialPlanarRigidWallMotion3 {
        MaterialPlanarRigidWallMotion3 {
            wall: self.wall,
            origin_m: array(origin),
            velocity_at_origin_m_s: array(self.tangent(q)),
            angular_speed_rad_s: q[2],
        }
    }

    fn admit_endpoint(&self, body: RigidBodyState) -> Result<(), PlanarWallError> {
        let bound = self
            .liquid
            .scaling()
            .velocity_si(self.liquid.scaling().flow().u_lattice)?;
        let actual = self
            .prescribed(body.pose().position_world(), self.body_twist(body)?)
            .maximum_speed_m_s(
                self.liquid.lattice().domain_cells(),
                self.liquid.scaling().cell_size_m(),
            )?;
        if actual > bound {
            return Err(invalid(
                "planar wall endpoint exceeds the original material speed domain",
            ));
        }
        Ok(())
    }

    /// Advance held tangential force [N] and axial couple [N m] with liquid.
    ///
    /// For q=(V_a,V_b,omega), the pilot measures J_wall(q0) and C=-dJ_wall/dq.
    /// The three simultaneous MIDPOINT changes satisfy
    /// (M+C/2) delta_q=(J_wall(q0)+J_external)/2, M=diag(m,m,I).
    /// The physical endpoint velocity jump is twice delta_q, up to represented
    /// native arithmetic. Final fluid replay supplies the actual MBD wrench;
    /// the solved twist is not substituted for the native endpoint state.
    #[allow(clippy::too_many_lines)] // One coupled candidate and its endpoint admission.
    pub fn step_candidate(
        &self,
        external_force_n: Vec3,
        external_torque_nm: f64,
        controls: PlanarWallControls,
        gate: &CancelGate,
    ) -> Result<PlanarWallStep<'a>, PlanarWallError> {
        poll(gate)?;
        controls.validate()?;
        let axis = normal_axis(self.wall);
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        if !external_force_n.is_finite() || array(external_force_n)[axis] != 0.0 {
            return Err(invalid(
                "planar wall applied force must be finite and tangential",
            ));
        }
        finite(external_torque_nm)?;
        let trial_cells = self
            .domain_cells
            .checked_mul(2)
            .ok_or_else(|| invalid("planar wall native trial work overflowed"))?;
        if trial_cells > controls.maximum_trial_cells {
            return Err(PlanarWallError::WorkBudget {
                required_trial_cells: trial_cells,
                maximum: controls.maximum_trial_cells,
            });
        }
        let dt = self.liquid.scaling().timestep_s();
        let q0 = self.body_twist(self.body)?;
        let v0 = self.tangent(q0);
        let initial_center = self.body.pose().position_world();
        let origin = checked_vector(initial_center.add(v0.scale(0.5 * dt)))?;
        if finite(norm(v0) * (0.5 * dt))? > controls.maximum_translation_m
            || finite(q0[2].abs() * (0.5 * dt))? > controls.maximum_rotation_rad
        {
            return Err(invalid(
                "planar wall first free half-drift exceeds its motion bound",
            ));
        }
        let external_linear = checked_vector(external_force_n.scale(dt))?;
        let external_spin = finite(external_torque_nm * dt)?;
        let external_couple = axis_vector(axis).scale(external_spin);
        let pilot = {
            let (candidate, report) = self.liquid.planar_rigid_wall_motion_candidate_si(
                self.prescribed(origin, q0),
                controls.maximum_relative_mass_defect,
                gate,
            )?;
            drop(candidate);
            report
        };
        poll(gate)?;
        let pilot_load = [
            pilot.wall_impulse_kg_m_s[a],
            pilot.wall_impulse_kg_m_s[b],
            pilot.wall_angular_impulse_kg_m2_s[axis],
        ];
        let external = [
            array(external_linear)[a],
            array(external_linear)[b],
            external_spin,
        ];
        let mass = [
            self.properties.mass(),
            self.properties.mass(),
            self.inertia_kg_m2,
        ];
        let roots = mass.map(f64::sqrt);
        let mut system = [0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                let c = pilot.fluid_response[i][j];
                let transposed = pilot.fluid_response[j][i];
                if !c.is_finite()
                    || !transposed.is_finite()
                    || (c - transposed).abs() > 64.0 * f64::EPSILON * c.abs().max(transposed.abs())
                    || (i == j && c < 0.0)
                {
                    return Err(invalid(
                        "native planar wall response is not a finite symmetric resistance",
                    ));
                }
                let symmetric = f64::midpoint(c, transposed);
                system[3 * i + j] = finite(
                    (0.5 * symmetric / roots[i]) / roots[j] + if i == j { 1.0 } else { 0.0 },
                )?;
            }
        }
        let mut change = [0.0; 3];
        for i in 0..3 {
            change[i] = finite(
                ((Dd::from_f64(pilot_load[i]) + Dd::from_f64(external[i])) * Dd::from_f64(0.5)
                    / Dd::from_f64(roots[i]))
                .to_f64(),
            )?;
        }
        let factor = cholesky(&system, 3).map_err(|_| {
            invalid("finite wall and liquid resistance is not representably positive definite")
        })?;
        factor.solve(&mut change);
        let mut imposed = [0.0; 3];
        for i in 0..3 {
            imposed[i] = finite(
                (Dd::from_f64(q0[i]) + Dd::from_f64(change[i]) / Dd::from_f64(roots[i])).to_f64(),
            )?;
        }
        let (liquid, wall) = self.liquid.planar_rigid_wall_motion_candidate_si(
            self.prescribed(origin, imposed),
            controls.maximum_relative_mass_defect,
            gate,
        )?;
        if wall.fluid_work_j.abs() > controls.maximum_absolute_wall_work_j {
            return Err(invalid(
                "planar wall work exceeds the declared exchange bound",
            ));
        }
        poll(gate)?;
        let raw_guide = normal_reaction(wall.wall_impulse_kg_m_s, axis, 1.0);
        let raw_guide_angular = transverse_reaction(wall.wall_angular_impulse_kg_m2_s, axis, 1.0);
        let linear_input = sum3(
            wall.wall_impulse_kg_m_s,
            array(raw_guide),
            array(external_linear),
        );
        let spin_input = sum3(
            wall.wall_angular_impulse_kg_m2_s,
            array(raw_guide_angular),
            array(external_couple),
        );
        let force = checked_vector(vector(
            linear_input.map(|j| (j / Dd::from_f64(dt)).to_f64()),
        ))?;
        let torque = finite((spin_input[axis] / Dd::from_f64(dt)).to_f64())?;
        let mechanics = RigidBodyIntegrator::new(Gravity::ZERO).step(
            self.body,
            self.properties,
            Wrench {
                force_world: force,
                torque_body: self.axis_body.scale(torque),
            },
            dt,
        )?;
        let body = mechanics.state_after;
        self.admit_endpoint(body)?;
        if array(body.linear_momentum_world())[axis] != 0.0
            || array(body.pose().position_world())[axis] != array(initial_center)[axis]
            || norm(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(self.axis_body)
                    .sub(axis_vector(axis)),
            ) > 256.0 * f64::EPSILON
        {
            return Err(invalid(
                "native wall step departed from its fixed planar guide",
            ));
        }
        let q1 = self.body_twist(body)?;
        if finite((0.5 * dt) * (norm(v0) + norm(self.tangent(q1))))?
            > controls.maximum_translation_m
            || finite((0.5 * dt) * (q0[2].abs() + q1[2].abs()))? > controls.maximum_rotation_rad
        {
            return Err(invalid(
                "actual native wall motion exceeds its two half-drift bounds",
            ));
        }
        let p0 = dd(array(self.body.linear_momentum_world()));
        let p1 = dd(array(body.linear_momentum_world()));
        let linear_increment = sub(p1, p0);
        let velocity0 =
            array(self.body.linear_momentum_world()).map(|p| p / self.properties.mass());
        let velocity1 = array(body.linear_momentum_world()).map(|p| p / self.properties.mass());
        let velocity_mid = average(velocity0, velocity1);
        let omega0 = self
            .properties
            .angular_velocity_body_checked(self.body.angular_momentum_body())?;
        let omega1 = self
            .properties
            .angular_velocity_body_checked(body.angular_momentum_body())?;
        let omega_world0 = self.body.pose().orientation().rotate_body_to_world(omega0);
        let omega_world1 = body.pose().orientation().rotate_body_to_world(omega1);
        let omega_mid = average(array(omega_world0), array(omega_world1));
        let spin0 = self
            .body
            .pose()
            .orientation()
            .rotate_body_to_world(self.body.angular_momentum_body());
        let spin1 = body
            .pose()
            .orientation()
            .rotate_body_to_world(body.angular_momentum_body());
        let spin_increment = sub(dd(array(spin1)), dd(array(spin0)));
        // Difference the actual endpoint positions/momenta before crossing:
        // delta(x cross p)=x0 cross delta_p + delta_x cross p1.
        let displacement = sub(
            dd(array(body.pose().position_world())),
            dd(array(initial_center)),
        );
        let orbital_increment = add(
            cross_dd(dd(array(initial_center)), linear_increment),
            cross_dd(displacement, p1),
        );
        let angular_increment = add(orbital_increment, spin_increment);
        let angular_input = add(cross_dd(dd(array(origin)), linear_input), spin_input);
        let momentum_residual = values(sub(linear_increment, linear_input));
        let spin_residual = values(sub(spin_increment, spin_input));
        let angular_residual = values(sub(angular_increment, angular_input));
        let midpoint_residual = values(sub(velocity_mid, dd(array(self.tangent(imposed)))));
        let midpoint_spin_residual = values(sub(
            omega_mid,
            dd(array(axis_vector(axis).scale(imposed[2]))),
        ));
        let final_load = [
            wall.wall_impulse_kg_m_s[a],
            wall.wall_impulse_kg_m_s[b],
            wall.wall_angular_impulse_kg_m2_s[axis],
        ];
        let response_residual = core::array::from_fn(|i| {
            let difference = Dd::from_f64(final_load[i]) - Dd::from_f64(pilot_load[i]);
            (0..3)
                .fold(difference, |sum, j| {
                    sum + Dd::from_f64(pilot.fluid_response[i][j])
                        * (Dd::from_f64(imposed[j]) - Dd::from_f64(q0[j]))
                })
                .to_f64()
        });
        let kinetic = dot_dd(velocity_mid, linear_increment)
            + dot_dd(
                average(array(omega0), array(omega1)),
                sub(
                    dd(array(body.angular_momentum_body())),
                    dd(array(self.body.angular_momentum_body())),
                ),
            );
        let external_work = dot_dd(velocity_mid, dd(array(external_linear)))
            + dot_dd(omega_mid, dd(array(external_couple)));
        let guide_work = dot_dd(velocity_mid, dd(array(raw_guide)))
            + dot_dd(omega_mid, dd(array(raw_guide_angular)));
        let work_residual = finite(
            (kinetic + Dd::from_f64(wall.fluid_work_j) - external_work - guide_work).to_f64(),
        )?;
        let kinetic_energy_change_j = finite(kinetic.to_f64())?;
        let external_work_j = finite(external_work.to_f64())?;
        let raw_guide_work_j = finite(guide_work.to_f64())?;
        let active_wall = sum3(wall.wall_impulse_kg_m_s, array(raw_guide), [0.0; 3]);
        let active_couple = sum3(
            wall.wall_angular_impulse_kg_m2_s,
            array(raw_guide_angular),
            [0.0; 3],
        );
        let linear_scale = norm(values(linear_increment))
            .max(norm(values(active_wall)))
            .max(norm(external_linear));
        let spin_scale = norm(values(spin_increment))
            .max(norm(values(active_couple)))
            .max(external_spin.abs());
        let angular_scale = norm(values(angular_increment))
            .max(norm(values(angular_input)))
            .max(norm(values(add(
                cross_dd(dd(array(origin)), active_wall),
                active_couple,
            ))))
            .max(norm(values(add(
                cross_dd(dd(array(origin)), dd(array(external_linear))),
                dd(array(external_couple)),
            ))));
        let speed_scale = norm(values(velocity_mid)).max(norm(self.tangent(imposed)));
        let spin_speed_scale = norm(values(omega_mid)).max(imposed[2].abs());
        let work_scale = kinetic_energy_change_j
            .abs()
            .max(wall.fluid_work_j.abs())
            .max(external_work_j.abs())
            .max(raw_guide_work_j.abs());
        let momentum_limit = allowance(
            controls.momentum_absolute_kg_m_s,
            controls.relative_tolerance,
            linear_scale,
        )?;
        let spin_limit = allowance(
            controls.angular_momentum_absolute_kg_m2_s,
            controls.relative_tolerance,
            spin_scale,
        )?;
        if !momentum_residual.is_finite()
            || !spin_residual.is_finite()
            || !angular_residual.is_finite()
            || !midpoint_residual.is_finite()
            || !midpoint_spin_residual.is_finite()
            || !response_residual.into_iter().all(f64::is_finite)
            || norm(momentum_residual) > momentum_limit
            || norm(spin_residual) > spin_limit
            || norm(angular_residual)
                > allowance(
                    controls.angular_momentum_absolute_kg_m2_s,
                    controls.relative_tolerance,
                    angular_scale,
                )?
            || response_residual[0].hypot(response_residual[1]) > momentum_limit
            || response_residual[2].abs() > spin_limit
            || norm(midpoint_residual)
                > allowance(
                    controls.midpoint_absolute_m_s,
                    controls.relative_tolerance,
                    speed_scale,
                )?
            || norm(midpoint_spin_residual)
                > allowance(
                    controls.midpoint_absolute_rad_s,
                    controls.relative_tolerance,
                    spin_speed_scale,
                )?
            || work_residual.abs()
                > allowance(
                    controls.work_absolute_j,
                    controls.relative_tolerance,
                    work_scale,
                )?
        {
            return Err(PlanarWallError::Balance {
                momentum_residual_kg_m_s: momentum_residual,
                spin_residual_kg_m2_s: spin_residual,
                angular_momentum_residual_kg_m2_s: angular_residual,
                midpoint_residual_m_s: midpoint_residual,
                midpoint_residual_rad_s: midpoint_spin_residual,
                response_residual,
                work_residual_j: work_residual,
            });
        }
        let time = finite(self.elapsed_time_s + dt)?;
        if time <= self.elapsed_time_s {
            return Err(invalid(
                "planar wall elapsed time cannot represent its next native tick",
            ));
        }
        let accepted_steps = self
            .accepted_steps
            .checked_add(1)
            .ok_or_else(|| invalid("planar wall step count overflowed"))?;
        let physical_guide = checked_vector(normal_reaction(
            wall.pressure_perturbation_wall_force_n,
            axis,
            dt,
        ))?;
        let physical_couple = checked_vector(transverse_reaction(
            wall.pressure_perturbation_wall_torque_nm,
            axis,
            dt,
        ))?;
        let absolute_guide = wall
            .absolute_wall_force_n
            .map(|f| checked_vector(normal_reaction(f, axis, dt)))
            .transpose()?;
        let absolute_couple = wall
            .absolute_wall_torque_nm
            .map(|m| checked_vector(transverse_reaction(m, axis, dt)))
            .transpose()?;
        poll(gate)?;
        Ok(PlanarWallStep {
            motion: Self {
                liquid,
                properties: self.properties,
                body,
                wall: self.wall,
                axis_body: self.axis_body,
                inertia_kg_m2: self.inertia_kg_m2,
                domain_cells: self.domain_cells,
                elapsed_time_s: time,
                accepted_steps,
            },
            wall,
            mechanics,
            coupling_origin_world_m: origin,
            external_linear_impulse_world: external_linear,
            external_couple_impulse_world: external_couple,
            raw_guide_impulse_world: raw_guide,
            raw_guide_angular_impulse_world: raw_guide_angular,
            pressure_perturbation_guide_impulse_world: physical_guide,
            pressure_perturbation_guide_angular_impulse_world: physical_couple,
            absolute_guide_impulse_world: absolute_guide,
            absolute_guide_angular_impulse_world: absolute_couple,
            midpoint_velocity_m_s: self.tangent(imposed),
            midpoint_angular_speed_rad_s: imposed[2],
            midpoint_residual_m_s: midpoint_residual,
            midpoint_residual_rad_s: midpoint_spin_residual,
            momentum_residual_kg_m_s: momentum_residual,
            spin_residual_kg_m2_s: spin_residual,
            angular_momentum_residual_kg_m2_s: angular_residual,
            response_residual,
            kinetic_energy_change_j,
            external_work_j,
            raw_guide_work_j,
            mechanical_work_residual_j: work_residual,
            native_trials: 2,
            trial_cells,
        })
    }
}

impl PlanarWallControls {
    fn validate(self) -> Result<(), PlanarWallError> {
        if [
            self.maximum_relative_mass_defect,
            self.maximum_absolute_wall_work_j,
            self.momentum_absolute_kg_m_s,
            self.angular_momentum_absolute_kg_m2_s,
            self.midpoint_absolute_m_s,
            self.midpoint_absolute_rad_s,
            self.work_absolute_j,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
            || [self.maximum_translation_m, self.maximum_rotation_rad]
                .into_iter()
                .any(|v| !v.is_finite() || v <= 0.0)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
        {
            return Err(invalid(
                "planar wall bounds need finite nonnegative allowances and positive motion limits",
            ));
        }
        Ok(())
    }
}

fn dd(a: [f64; 3]) -> [Dd; 3] {
    a.map(Dd::from_f64)
}
fn values(a: [Dd; 3]) -> Vec3 {
    vector(a.map(Dd::to_f64))
}
fn add(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn average(a: [f64; 3], b: [f64; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| (Dd::from_f64(a[i]) + Dd::from_f64(b[i])) * Dd::from_f64(0.5))
}
fn sum3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [Dd; 3] {
    core::array::from_fn(|i| Dd::from_f64(a[i]) + Dd::from_f64(b[i]) + Dd::from_f64(c[i]))
}
fn cross_dd(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot_dd(a: [Dd; 3], b: [Dd; 3]) -> Dd {
    (0..3).fold(Dd::ZERO, |sum, i| sum + a[i] * b[i])
}
fn axis_vector(axis: usize) -> Vec3 {
    vector(core::array::from_fn(|d| if d == axis { 1.0 } else { 0.0 }))
}
fn normal_reaction(force: [f64; 3], axis: usize, scale: f64) -> Vec3 {
    vector(core::array::from_fn(|d| {
        if d == axis { -force[d] * scale } else { 0.0 }
    }))
}
fn transverse_reaction(moment: [f64; 3], axis: usize, scale: f64) -> Vec3 {
    vector(core::array::from_fn(|d| {
        if d == axis { 0.0 } else { -moment[d] * scale }
    }))
}
fn checked_vector(value: Vec3) -> Result<Vec3, PlanarWallError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable planar wall vector"))
    }
}
fn allowance(absolute: f64, relative: f64, scale: f64) -> Result<f64, PlanarWallError> {
    finite(absolute + relative * scale)
}
fn invalid(reason: &'static str) -> PlanarWallError {
    PlanarWallError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64, PlanarWallError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable planar wall arithmetic"))
    }
}
fn poll(gate: &CancelGate) -> Result<(), PlanarWallError> {
    if gate.is_requested() {
        Err(PlanarWallError::Cancelled)
    } else {
        Ok(())
    }
}
