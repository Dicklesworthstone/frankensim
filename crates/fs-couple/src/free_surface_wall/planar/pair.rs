//! Two opposing finite walls on one original liquid collision and transport tick.
//!
//! Each wall retains two tangent translations and its principal normal spin.
//! The native explicit boundary map has two independent three-mode response
//! blocks at a fixed initial fluid state. One shared pilot measures them; the
//! original mass-whitened solves determine both midpoint twists, and one shared
//! replay advances the liquid once. Momentum communicates between the walls
//! through the evolving original populations over subsequent physical ticks.

/// Both finite walls on the original source-driven thermal liquid trajectory.
pub mod thermal;

use super::{
    PlanarWallCandidate, PlanarWallControls, PlanarWallError, PlanarWallLiquid, PlanarWallMotion,
    add, allowance, array, cross_dd, dd, finite, norm, normal_axis, sum3, values,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialPlanarRigidWallMotion3, MaterialPlanarRigidWallPairReport3,
};
use fs_math::dd::Dd;
use fs_mbd::{MassProperties, RigidBodyState, Vec3};

/// One mechanical owner shared by the original isothermal and thermal liquids.
/// The source-bound liquid is stored once; each wall retains only its mechanics.
#[derive(Clone)]
pub struct PlanarWallPairMotion<L> {
    liquid: L,
    walls: [PlanarWallMotion<()>; 2],
}

/// Original source-bound isothermal liquid and both finite wall mechanisms.
pub type FreeSurfacePlanarWallPair<'a> = PlanarWallPairMotion<MaterialFreeSurface3<'a>>;

/// Refusal publishes neither wall, liquid history nor physical clock.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarWallPairError {
    /// The selected planes or combined exchange are not admissible.
    Invalid(&'static str),
    /// Original per-wall mechanics, with the caller's wall index retained.
    Wall {
        index: usize,
        error: PlanarWallError,
    },
    /// Original joint liquid candidate refused before publication.
    Fluid(MaterialLbmError),
    /// Original common thermal source, transport or endpoint material refused.
    Thermal(FreeSurfaceEnthalpyError),
    /// Boundary exchange failed its combined physical balance.
    Balance {
        momentum_residual_kg_m_s: Vec3,
        angular_momentum_residual_kg_m2_s: Vec3,
        work_residual_j: f64,
    },
    /// Cancellation before complete publication.
    Cancelled,
}

impl core::fmt::Display for PlanarWallPairError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite opposing liquid walls refused: {self:?}")
    }
}
impl std::error::Error for PlanarWallPairError {}
impl From<MaterialLbmError> for PlanarWallPairError {
    fn from(error: MaterialLbmError) -> Self {
        Self::Fluid(error)
    }
}
impl From<FreeSurfaceEnthalpyError> for PlanarWallPairError {
    fn from(error: FreeSurfaceEnthalpyError) -> Self {
        Self::Thermal(error)
    }
}

/// Both accepted finite walls and exactly one accepted native liquid tick.
pub struct PlanarWallPairCandidate<L, R = ()> {
    /// Continue from this jointly accepted state.
    pub motion: PlanarWallPairMotion<L>,
    /// The one final native liquid replay, in the caller's wall order.
    pub fluid: MaterialPlanarRigidWallPairReport3,
    /// Original source/transport account from the final replay, never the pilot.
    pub liquid_report: R,
    /// Complete original per-wall mechanics, raw/physical guide reactions,
    /// spin/orbital impulses and work. Their mechanical records own no liquid.
    /// Each account's two trial counts refer to the SAME pilot and replay.
    pub walls: [PlanarWallCandidate<()>; 2],
    /// Actual combined body and liquid boundary exchange minus external and
    /// guide impulses [kg m/s]. This does not account for unrelated bulk forces.
    pub momentum_residual_kg_m_s: Vec3,
    /// The same exchange including both spins and orbital moments, about the
    /// fixed world origin [kg m²/s].
    pub angular_momentum_residual_kg_m2_s: Vec3,
    /// Sum of actual native body kinetic changes [J].
    pub kinetic_energy_change_j: f64,
    /// Sum of applied forces/couples paired with actual midpoint motion [J].
    pub external_work_j: f64,
    /// Sum of raw guide work on the actual native midpoints [J].
    pub raw_guide_work_j: f64,
    /// Body kinetic change plus native liquid boundary work minus applied and
    /// guide work [J]. Numerical defects are never converted into heat.
    pub mechanical_work_residual_j: f64,
    /// Exactly two complete shared liquid candidates, not two per wall.
    pub native_trials: usize,
    /// Complete-domain cell count across the shared pilot and replay.
    pub trial_cells: usize,
}

/// Complete original isothermal liquid/body pair candidate.
pub type PlanarWallPairStep<'a> = PlanarWallPairCandidate<MaterialFreeSurface3<'a>>;

impl<L: PlanarWallLiquid> PlanarWallPairMotion<L> {
    /// Bind distinct opposing planes of one axis to two original finite bodies.
    /// The original single-wall principal-axis, tangent-motion, speed, material
    /// and cell-domain admission applies independently to each body.
    pub fn new(
        liquid: L,
        properties: [MassProperties; 2],
        bodies: [RigidBodyState; 2],
        planes: [PlanarWall3; 2],
        maximum_cells: usize,
        gate: &CancelGate,
    ) -> Result<Self, PlanarWallPairError> {
        poll(gate)?;
        if planes[0] == planes[1] || normal_axis(planes[0]) != normal_axis(planes[1]) {
            return Err(PlanarWallPairError::Invalid(
                "finite walls require distinct opposing planes of one axis",
            ));
        }
        let prepare = |i: usize| {
            PlanarWallMotion::new(
                &liquid,
                properties[i],
                bodies[i],
                planes[i],
                maximum_cells,
                gate,
            )
            .map(|wall| wall.with_liquid(()))
            .map_err(|error| wall_error(i, error))
        };
        let walls = [prepare(0)?, prepare(1)?];
        poll(gate)?;
        Ok(Self { liquid, walls })
    }

    /// The only accepted liquid population history.
    pub const fn liquid(&self) -> &L {
        &self.liquid
    }
    /// Both original native body states, in constructor order.
    pub fn bodies(&self) -> [RigidBodyState; 2] {
        self.walls.each_ref().map(|wall| wall.body())
    }
    /// Both original finite mass and inertia values.
    pub fn properties(&self) -> [MassProperties; 2] {
        self.walls.each_ref().map(|wall| wall.properties())
    }
    /// Both fixed domain planes in constructor order.
    pub fn planes(&self) -> [PlanarWall3; 2] {
        self.walls.each_ref().map(|wall| wall.wall)
    }
    /// Shared physical time since binding.
    pub const fn elapsed_time_s(&self) -> f64 {
        self.walls[0].elapsed_time_s
    }
    /// Number of jointly accepted liquid and two-body intervals.
    pub const fn accepted_steps(&self) -> u64 {
        self.walls[0].accepted_steps
    }
    /// Both actual tangent velocities and signed normal angular speeds.
    pub fn twists(&self) -> Result<[(Vec3, f64); 2], PlanarWallPairError> {
        let twist = |i: usize| {
            self.walls[i]
                .with_liquid(&self.liquid)
                .twist()
                .map_err(|error| wall_error(i, error))
        };
        Ok([twist(0)?, twist(1)?])
    }

    // Both native liquid owners measure and replay this same mechanical solve.
    // The pilot and final candidate always begin at the unchanged original
    // liquid; only the final source/transport report is retained.
    fn step_candidate_with<R>(
        &self,
        external_forces_n: [Vec3; 2],
        external_torques_nm: [f64; 2],
        controls: PlanarWallControls,
        gate: &CancelGate,
        mut trial: impl FnMut(
            [MaterialPlanarRigidWallMotion3; 2],
            f64,
            &CancelGate,
        )
            -> Result<(L, MaterialPlanarRigidWallPairReport3, R), PlanarWallPairError>,
    ) -> Result<PlanarWallPairCandidate<L, R>, PlanarWallPairError> {
        poll(gate)?;
        let walls = self
            .walls
            .each_ref()
            .map(|wall| wall.with_liquid(&self.liquid));
        let prepare = |i: usize| {
            walls[i]
                .prepare_candidate(external_forces_n[i], external_torques_nm[i], controls, gate)
                .map_err(|error| wall_error(i, error))
        };
        let prepared = [prepare(0)?, prepare(1)?];
        let pilot = {
            let motions =
                core::array::from_fn(|i| walls[i].prescribed(prepared[i].origin, prepared[i].q0));
            let (candidate, report, _) =
                trial(motions, controls.maximum_relative_mass_defect, gate)?;
            drop(candidate);
            report
        };
        poll(gate)?;
        let midpoint = |i: usize| {
            walls[i]
                .solve_midpoint(&prepared[i], &pilot.walls[i])
                .map_err(|error| wall_error(i, error))
        };
        let imposed = [midpoint(0)?, midpoint(1)?];
        let motions = core::array::from_fn(|i| walls[i].prescribed(prepared[i].origin, imposed[i]));
        let (liquid, fluid, liquid_report) =
            trial(motions, controls.maximum_relative_mass_defect, gate)?;
        if fluid.walls[0].fluid_work_j.abs() + fluid.walls[1].fluid_work_j.abs()
            > controls.maximum_absolute_wall_work_j
        {
            return Err(PlanarWallPairError::Invalid(
                "combined absolute wall work exceeds its bound",
            ));
        }
        let finish = |i: usize| {
            walls[i]
                .finish_candidate(
                    &prepared[i],
                    &pilot.walls[i],
                    imposed[i],
                    (),
                    fluid.walls[i],
                    (),
                    controls,
                    gate,
                )
                .map_err(|error| wall_error(i, error))
        };
        let accounts = [finish(0)?, finish(1)?];
        // Each retained per-wall defect already uses actual endpoint differences
        // and polarized products. Add the original pair boundary closure, not
        // differences of large absolute kinetic or angular inventories.
        let mut linear = dd(fluid.fluid_impulse_kg_m_s);
        let mut angular = dd(fluid.fluid_angular_impulse_at_si([0.0; 3])?);
        let mut kinetic = Dd::ZERO;
        let mut external_work = Dd::ZERO;
        let mut guide_work = Dd::ZERO;
        let mut work = Dd::from_f64(fluid.fluid_work_j);
        let mut linear_scale = 0.0;
        let mut angular_scale = 0.0;
        let mut work_scale = 0.0;
        for account in &accounts {
            let origin = dd(array(account.coupling_origin_world_m));
            linear = add(
                linear,
                add(
                    dd(account.wall.wall_impulse_kg_m_s),
                    dd(array(account.momentum_residual_kg_m_s)),
                ),
            );
            angular = add(
                angular,
                add(
                    add(
                        cross_dd(origin, dd(account.wall.wall_impulse_kg_m_s)),
                        dd(account.wall.wall_angular_impulse_kg_m2_s),
                    ),
                    dd(array(account.angular_momentum_residual_kg_m2_s)),
                ),
            );
            kinetic = kinetic + Dd::from_f64(account.kinetic_energy_change_j);
            external_work = external_work + Dd::from_f64(account.external_work_j);
            guide_work = guide_work + Dd::from_f64(account.raw_guide_work_j);
            work = work + Dd::from_f64(account.mechanical_work_residual_j)
                - Dd::from_f64(account.wall.fluid_work_j);
            let active = sum3(
                account.wall.wall_impulse_kg_m_s,
                array(account.raw_guide_impulse_world),
                [0.0; 3],
            );
            let couple = sum3(
                account.wall.wall_angular_impulse_kg_m2_s,
                array(account.raw_guide_angular_impulse_world),
                [0.0; 3],
            );
            let force_scale = norm(values(active)).max(norm(account.external_linear_impulse_world));
            linear_scale += force_scale;
            let active_angular = add(cross_dd(origin, active), couple);
            let external_angular = add(
                cross_dd(origin, dd(array(account.external_linear_impulse_world))),
                dd(array(account.external_couple_impulse_world)),
            );
            angular_scale += norm(values(active_angular)).max(norm(values(external_angular)));
            work_scale += account
                .kinetic_energy_change_j
                .abs()
                .max(account.wall.fluid_work_j.abs())
                .max(account.external_work_j.abs())
                .max(account.raw_guide_work_j.abs());
        }
        let momentum_residual = values(linear);
        let angular_residual = values(angular);
        let work_residual = work.to_f64();
        let numeric = |result: Result<f64, PlanarWallError>| {
            result
                .map_err(|_| PlanarWallPairError::Invalid("unrepresentable combined wall account"))
        };
        let linear_limit = numeric(allowance(
            2.0 * controls.momentum_absolute_kg_m_s,
            controls.relative_tolerance,
            linear_scale,
        ))?;
        let angular_limit = numeric(allowance(
            2.0 * controls.angular_momentum_absolute_kg_m2_s,
            controls.relative_tolerance,
            angular_scale,
        ))?;
        let work_limit = numeric(allowance(
            2.0 * controls.work_absolute_j,
            controls.relative_tolerance,
            work_scale,
        ))?;
        if !momentum_residual.is_finite()
            || !angular_residual.is_finite()
            || !work_residual.is_finite()
            || norm(momentum_residual) > linear_limit
            || norm(angular_residual) > angular_limit
            || work_residual.abs() > work_limit
        {
            return Err(PlanarWallPairError::Balance {
                momentum_residual_kg_m_s: momentum_residual,
                angular_momentum_residual_kg_m2_s: angular_residual,
                work_residual_j: work_residual,
            });
        }
        let kinetic_energy_change_j = numeric(finite(kinetic.to_f64()))?;
        let external_work_j = numeric(finite(external_work.to_f64()))?;
        let raw_guide_work_j = numeric(finite(guide_work.to_f64()))?;
        poll(gate)?;
        Ok(PlanarWallPairCandidate {
            motion: Self {
                liquid,
                walls: accounts
                    .each_ref()
                    .map(|account| account.motion.with_liquid(())),
            },
            fluid,
            liquid_report,
            walls: accounts,
            momentum_residual_kg_m_s: momentum_residual,
            angular_momentum_residual_kg_m2_s: angular_residual,
            kinetic_energy_change_j,
            external_work_j,
            raw_guide_work_j,
            mechanical_work_residual_j: work_residual,
            native_trials: 2,
            trial_cells: prepared[0].trial_cells,
        })
    }
}

impl<'a> PlanarWallPairMotion<MaterialFreeSurface3<'a>> {
    /// Advance both finite walls under their held tangent forces and axial couples.
    ///
    /// The shared cell budget counts two full-grid trials. Translation/rotation
    /// bounds and physical tolerances apply to each original wall; the wall-work
    /// ceiling additionally bounds the sum of both absolute boundary works.
    /// Both planes must remain flooded in the original native model. The planes
    /// stay fixed: there is no normal displacement, tilt, changing coverage or
    /// finite plate edge. This is an isothermal liquid step; mechanical boundary
    /// work is not a fluid energy-conservation or heat-production claim.
    pub fn step_candidate(
        &self,
        external_forces_n: [Vec3; 2],
        external_torques_nm: [f64; 2],
        controls: PlanarWallControls,
        gate: &CancelGate,
    ) -> Result<PlanarWallPairStep<'a>, PlanarWallPairError> {
        self.step_candidate_with(
            external_forces_n,
            external_torques_nm,
            controls,
            gate,
            |motions, mass_tolerance, gate| {
                let (candidate, report) = self.liquid.planar_rigid_wall_pair_candidate_si(
                    motions,
                    mass_tolerance,
                    gate,
                )?;
                Ok((candidate, report, ()))
            },
        )
    }
}

fn wall_error(index: usize, error: PlanarWallError) -> PlanarWallPairError {
    if error == PlanarWallError::Cancelled {
        PlanarWallPairError::Cancelled
    } else {
        PlanarWallPairError::Wall { index, error }
    }
}
fn poll(gate: &CancelGate) -> Result<(), PlanarWallPairError> {
    if gate.is_requested() {
        Err(PlanarWallPairError::Cancelled)
    } else {
        Ok(())
    }
}
