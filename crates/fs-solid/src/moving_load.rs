//! Moving force/couple ports on the original six-DOF flat shell surrogate.
//!
//! P1 material-point motion and its exact transpose distribute each physical
//! wrench, including its reference offset. Traversing a patch changes the next
//! load stencil; it never adds the traversal velocity to material velocity.
//! The original spatial matrices and prepared average-acceleration Newmark
//! owner advance the motion. This retains the shell owner's linear, small-
//! rotation, estimate-only scope: angular momentum is about reference positions,
//! not a finite-rotation shell or an inferred contact law.

use crate::shell::{ShellAssembly, ShellError, ShellPlate};
use crate::shell_time::{DynamicState, PreparedNewmark, ShellTimeConfig, ShellTimeError};
use fs_la::LuWorkspace;

#[path = "moving_load/driven.rs"]
mod driven;
pub use driven::ShellSupportMotion;

#[path = "moving_load/payload.rs"]
mod payload;
pub use payload::ShellRigidPayload;

#[path = "moving_load/impulse.rs"]
mod impulse;
pub use impulse::{MovingShellImpulse, ShellPointImpulse};

#[path = "moving_load/inverse_inertia.rs"]
mod inverse_inertia;
pub use inverse_inertia::ShellPointInverseInertia;

#[path = "moving_load/normal_impact.rs"]
mod normal_impact;
pub use normal_impact::{
    PreparedShellNormalImpacts, ShellNormalImpactControls, ShellNormalImpactError,
    ShellNormalImpactReport, ShellNormalImpactSite, ShellObstacleMotion,
};

#[path = "moving_load/frictional_impact.rs"]
mod frictional_impact;
pub use frictional_impact::{
    NormalVelocityLaw, PreparedShellFrictionalImpacts, ShellFrictionalImpactControls,
    ShellFrictionalImpactError, ShellFrictionalImpactReport, ShellFrictionalImpactSite,
};

/// A reference material location on one original triangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellLoadPoint {
    /// Original triangle index. Triangle identity is retained at shared edges.
    pub triangle: usize,
    /// Nonnegative P1 weights summing to one; no extrapolation or normalization.
    pub barycentric: [f64; 3],
    /// Reference Cartesian arm from the midsurface point, in metres.
    /// For a surface load this may be the signed half-thickness normal.
    pub offset_m: [f64; 3],
}

/// One physical wrench. A distributed patch supplies its integrated quadrature
/// forces and couples as several entries; no area normalization is performed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellPointWrench {
    /// Reference material point and force lever arm.
    pub point: ShellLoadPoint,
    /// Cartesian force ON the shell, N.
    pub force_n: [f64; 3],
    /// Cartesian free couple ON the shell about that point, N m.
    pub couple_nm: [f64; 3],
}

/// Force/couple or linear/angular impulse in the fixed Cartesian frame.
/// Moments are about the same reference world origin as the shell nodes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ShellSpatialResultant {
    /// Force [N], or linear impulse/momentum [N s], as named by the owner field.
    pub linear: [f64; 3],
    /// Moment [N m], or angular impulse/momentum [N m s].
    pub angular: [f64; 3],
}

/// Material-point motion sampled with the SAME rows that transpose the load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellPointMotion {
    /// Reference location including the supplied arm, m.
    pub reference_position_m: [f64; 3],
    /// First-order displacement including director rotation of the arm, m.
    pub displacement_m: [f64; 3],
    /// Reference position plus first-order displacement, m.
    pub position_m: [f64; 3],
    /// Interpolated small Cartesian rotation, rad.
    pub rotation_rad: [f64; 3],
    /// Material velocity, without a moving-stencil convective contribution, m/s.
    pub velocity_m_s: [f64; 3],
    /// Interpolated Cartesian angular velocity, rad/s.
    pub angular_velocity_rad_s: [f64; 3],
    /// Unit direction of the first-order rotated reference triangle normal.
    pub normal: [f64; 3],
}

/// Full loads retain restrained rows; reduced loads enter the native time solve.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellLoadProjection {
    /// Node-major `[Fx,Fy,Fz,Mx,My,Mz]`, N and N m.
    pub full: Vec<f64>,
    /// The original assembly's free-DOF restriction of `full`.
    pub reduced: Vec<f64>,
    /// Independently reconstructed nodal force and reference-origin moment.
    pub resultant: ShellSpatialResultant,
}

/// No failed call publishes new motion or changes its accepted endpoint load.
#[derive(Debug)]
pub enum MovingLoadError {
    /// Invalid geometry location, load, state, timestamp or step budget.
    Invalid(&'static str),
    /// Original spatial shell owner refused.
    Spatial(ShellError),
    /// Original native time or inertia owner refused.
    Time(ShellTimeError),
    /// Cancellation observed before publication.
    Cancelled,
}
impl core::fmt::Display for MovingLoadError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moving shell load refused: {self:?}")
    }
}
impl std::error::Error for MovingLoadError {}
impl From<ShellTimeError> for MovingLoadError {
    fn from(error: ShellTimeError) -> Self {
        Self::Time(error)
    }
}

/// Original geometry and assembly bound once, without accepting unrelated
/// externally assembled operators. Work remains bounded by the shell budget.
#[derive(Debug)]
pub struct MovingShellModel {
    plate: ShellPlate,
    assembly: ShellAssembly,
    free_index: Vec<Option<usize>>,
    max_load_points: usize,
    payloads: Vec<ShellRigidPayload>,
}

impl MovingShellModel {
    /// Assemble one fixed shell and bound each load sample to at most 4096
    /// explicit quadrature points. Cancellation brackets the original bounded
    /// spatial assembly and is polled between point loads and native solves.
    pub fn new(
        plate: ShellPlate,
        max_load_points: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MovingLoadError> {
        poll(&mut cancelled)?;
        if max_load_points == 0 || max_load_points > 4096 {
            return Err(invalid("moving shell load bound must be in 1..=4096"));
        }
        let assembly = plate.assemble().map_err(MovingLoadError::Spatial)?;
        let mut free_index = vec![None; assembly.full_mass.dimension()];
        for (free, &full) in assembly.free_dofs.iter().enumerate() {
            free_index[full] = Some(free);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            plate,
            assembly,
            free_index,
            max_load_points,
            payloads: Vec::new(),
        })
    }

    /// Validated, immutable original spatial model.
    pub fn plate(&self) -> &ShellPlate {
        &self.plate
    }
    /// Native full/reduced operators, including any attached payload inertia,
    /// and the original fixed support membership.
    pub fn assembly(&self) -> &ShellAssembly {
        &self.assembly
    }

    fn point(&self, point: ShellLoadPoint) -> Result<([usize; 3], [f64; 3]), MovingLoadError> {
        let nodes = *self
            .plate
            .triangles
            .get(point.triangle)
            .ok_or_else(|| invalid("load triangle is outside the original shell"))?;
        if point
            .barycentric
            .iter()
            .any(|w| !w.is_finite() || *w < 0.0 || *w > 1.0)
            || (point.barycentric.iter().sum::<f64>() - 1.0).abs() > 8.0 * f64::EPSILON
            || point.offset_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite offset and convex unit-sum barycentric weights required",
            ));
        }
        let mut position = point.offset_m;
        for (&node, &weight) in nodes.iter().zip(&point.barycentric) {
            for (x, value) in position.iter_mut().zip(self.plate.nodes[node].position_m) {
                *x += weight * value;
            }
        }
        finite(&position)?;
        Ok((nodes, position))
    }

    /// Transfer physical force/couple by the transpose of material-point motion.
    /// Each nodal couple includes `offset cross force`. Restrained rows survive
    /// in `full`; no force or torque is dropped at a pinned node or shared edge.
    pub fn project(
        &self,
        loads: &[ShellPointWrench],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ShellLoadProjection, MovingLoadError> {
        self.project_components(
            loads.len(),
            loads
                .iter()
                .map(|load| (load.point, load.force_n, load.couple_nm)),
            &mut cancelled,
        )
    }

    // Forces and event impulses use the SAME material kinematics transpose.
    // The caller owns whether the components are N / N m or N s / N m s.
    fn project_components(
        &self,
        count: usize,
        components: impl IntoIterator<Item = (ShellLoadPoint, [f64; 3], [f64; 3])>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ShellLoadProjection, MovingLoadError> {
        poll(cancelled)?;
        if count > self.max_load_points {
            return Err(invalid("moving shell load point budget exceeded"));
        }
        let mut full = vec![0.0; self.free_index.len()];
        for (point, force, couple) in components {
            poll(cancelled)?;
            let (nodes, _) = self.point(point)?;
            finite(&force)?;
            finite(&couple)?;
            let moment = add(couple, cross(point.offset_m, force));
            finite(&moment)?;
            for (&node, &weight) in nodes.iter().zip(&point.barycentric) {
                for axis in 0..3 {
                    full[6 * node + axis] += weight * force[axis];
                    full[6 * node + 3 + axis] += weight * moment[axis];
                }
            }
        }
        finite(&full)?;
        let reduced = self.assembly.free_dofs.iter().map(|&i| full[i]).collect();
        let resultant = self.resultant(&full)?;
        poll(cancelled)?;
        Ok(ShellLoadProjection {
            full,
            reduced,
            resultant,
        })
    }

    /// Sample reduced displacement/velocity at one material point. Fixed support
    /// components are exactly zero. A different point in the next call does not
    /// turn path traversal into a material velocity or an extra work source.
    pub fn observe(
        &self,
        point: ShellLoadPoint,
        displacement: &[f64],
        velocity: &[f64],
    ) -> Result<ShellPointMotion, MovingLoadError> {
        self.observe_with_supports(point, displacement, velocity, None)
    }

    fn observe_with_supports(
        &self,
        point: ShellLoadPoint,
        displacement: &[f64],
        velocity: &[f64],
        supports: Option<&ShellSupportMotion>,
    ) -> Result<ShellPointMotion, MovingLoadError> {
        self.check_field(displacement)?;
        self.check_field(velocity)?;
        let (nodes, reference_position_m) = self.point(point)?;
        let q = self.material_components(
            point,
            nodes,
            displacement,
            supports.map(|s| &s.displacement_m),
        );
        let v = self.material_components(point, nodes, velocity, supports.map(|s| &s.velocity_m_s));
        let rotation_rad = [q[3], q[4], q[5]];
        let angular_velocity_rad_s = [v[3], v[4], v[5]];
        let displacement_m = [q[0], q[1], q[2]];
        let velocity_m_s = [v[0], v[1], v[2]];
        let [a, b, c] = nodes.map(|i| self.plate.nodes[i].position_m);
        let normal = unit(cross(sub(b, a), sub(c, a)))?;
        let normal = unit(add(normal, cross(rotation_rad, normal)))?;
        let position_m = add(reference_position_m, displacement_m);
        for value in [&q[..], &v, &position_m, &displacement_m, &velocity_m_s] {
            finite(value)?;
        }
        Ok(ShellPointMotion {
            reference_position_m,
            displacement_m,
            position_m,
            rotation_rad,
            velocity_m_s,
            angular_velocity_rad_s,
            normal,
        })
    }

    // Exact material-point kinematics paired with project_components.
    // Callers validate the original point, free field and prescribed values.
    fn material_components(
        &self,
        point: ShellLoadPoint,
        nodes: [usize; 3],
        field: &[f64],
        supports: Option<&[[f64; 3]; 3]>,
    ) -> [f64; 6] {
        let mut result = [0.0; 6];
        for (&node, &weight) in nodes.iter().zip(&point.barycentric) {
            for (axis, value) in result.iter_mut().enumerate() {
                let row = 6 * node + axis;
                if let Some(free) = self.free_index[row] {
                    *value += weight * field[free];
                } else if let Some(supports) = supports {
                    *value += weight * self.support_value(row, supports);
                }
            }
        }
        let offset = cross([result[3], result[4], result[5]], point.offset_m);
        for axis in 0..3 {
            result[axis] += offset[axis];
        }
        result
    }

    /// Start at t=0 with the original free-DOF displacement and velocity. The
    /// initial acceleration solves `M a = f - C v - K u`; it is never guessed
    /// zero. All subsequent steps reuse one native effective-mass factorization.
    /// Only average-acceleration parameters (beta=1/4, gamma=1/2) are admitted.
    pub fn start<'a>(
        &'a self,
        displacement: &[f64],
        velocity: &[f64],
        initial_loads: &[ShellPointWrench],
        config: ShellTimeConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellMotion<'a>, MovingLoadError> {
        self.start_with_supports(
            displacement,
            velocity,
            initial_loads,
            config,
            None,
            cancelled,
        )
    }

    fn start_with_supports<'a>(
        &'a self,
        displacement: &[f64],
        velocity: &[f64],
        initial_loads: &[ShellPointWrench],
        config: ShellTimeConfig,
        supports: Option<ShellSupportMotion>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellMotion<'a>, MovingLoadError> {
        poll(&mut cancelled)?;
        self.check_field(displacement)?;
        self.check_field(velocity)?;
        if let Some(motion) = &supports {
            self.check_supports(motion)?;
        }
        if config.beta != 0.25
            || config.gamma != 0.5
            || config.max_steps == 0
            || !config.duration_s.is_finite()
            || config.duration_s <= 0.0
        {
            return Err(invalid(
                "moving shell requires bounded average-acceleration Newmark",
            ));
        }
        let load = self.project(initial_loads, &mut cancelled)?;
        let a = &self.assembly;
        let n = a.mass.dimension();
        let rhs = self.acceleration_rhs(displacement, velocity, load.reduced, supports.as_ref())?;
        let mut acceleration = vec![0.0; n];
        let mut mass_solver = LuWorkspace::new(n).map_err(linear)?;
        mass_solver.factor(a.mass.values()).map_err(linear)?;
        mass_solver
            .solve_factored_into(&rhs, &mut acceleration)
            .map_err(linear)?;
        poll(&mut cancelled)?;
        let kinetic_energy_j = a.mass.quadratic_energy(velocity);
        let strain_energy_j = a.stiffness.quadratic_energy(displacement);
        let total_energy_j = kinetic_energy_j + strain_energy_j;
        finite(&[kinetic_energy_j, strain_energy_j, total_energy_j])?;
        let mut state = DynamicState {
            time_s: 0.0,
            step: 0,
            displacement: displacement.to_vec(),
            velocity: velocity.to_vec(),
            acceleration,
            kinetic_energy_j,
            strain_energy_j,
            total_energy_j,
        };
        if let Some(motion) = &supports {
            self.refresh_full_energy(&mut state, motion)?;
        }
        let reactions = self.reactions(&state, supports.as_ref(), &load.full)?;
        let integrator = PreparedNewmark::new(
            n,
            a.mass.values(),
            a.damping.as_ref().map(|c| c.values()),
            a.stiffness.values(),
            &config,
        )?;
        poll(&mut cancelled)?;
        Ok(MovingShellMotion {
            model: self,
            integrator,
            mass_solver,
            config,
            state,
            supports,
            full_load: load.full,
            reactions,
        })
    }

    // The original held-force equilibrium owns both initialization and the
    // post-impulse acceleration. Prescribed support acceleration stays held.
    fn acceleration_rhs(
        &self,
        displacement: &[f64],
        velocity: &[f64],
        load: Vec<f64>,
        supports: Option<&ShellSupportMotion>,
    ) -> Result<Vec<f64>, MovingLoadError> {
        let mut rhs = self.effective_load(load, supports)?;
        let elastic = self.assembly.stiffness.apply(displacement);
        let viscous = self.assembly.damping.as_ref().map(|c| c.apply(velocity));
        for (i, value) in rhs.iter_mut().enumerate() {
            *value -= elastic[i] + viscous.as_ref().map_or(0.0, |c| c[i]);
        }
        finite(&rhs)?;
        Ok(rhs)
    }

    fn check_field(&self, field: &[f64]) -> Result<(), MovingLoadError> {
        if field.len() != self.assembly.free_dofs.len() {
            return Err(invalid(
                "motion field must cover the original free shell DOFs",
            ));
        }
        finite(field)
    }
    fn expand(&self, field: &[f64]) -> Vec<f64> {
        self.free_index
            .iter()
            .map(|i| i.map_or(0.0, |i| field[i]))
            .collect()
    }
    fn resultant(&self, full: &[f64]) -> Result<ShellSpatialResultant, MovingLoadError> {
        self.resultant_by(|i| full[i])
    }
    fn resultant_by(
        &self,
        value: impl Fn(usize) -> f64,
    ) -> Result<ShellSpatialResultant, MovingLoadError> {
        let mut result = ShellSpatialResultant::default();
        for (index, node) in self.plate.nodes.iter().enumerate() {
            let row: [f64; 6] = std::array::from_fn(|i| value(6 * index + i));
            let force = [row[0], row[1], row[2]];
            result.linear = add(result.linear, force);
            result.angular = add(
                result.angular,
                add(cross(node.position_m, force), [row[3], row[4], row[5]]),
            );
        }
        finite(&result.linear)?;
        finite(&result.angular)?;
        Ok(result)
    }
    fn reactions(
        &self,
        state: &DynamicState,
        supports: Option<&ShellSupportMotion>,
        load: &[f64],
    ) -> Result<Vec<f64>, MovingLoadError> {
        let n = self.free_index.len();
        let mut reaction = vec![0.0; n];
        // Only fixed rows need a reaction. Keep the original full M/K/C
        // coupling, without three complete dense actions for a free body or
        // for the shell owner's nine restrained translational components.
        let row_action = |matrix: &[f64], row: usize, field: &[f64]| -> f64 {
            self.assembly
                .free_dofs
                .iter()
                .zip(field)
                .map(|(&column, value)| matrix[row * n + column] * value)
                .sum()
        };
        for i in 0..n {
            if self.free_index[i].is_none() {
                reaction[i] = row_action(self.assembly.full_mass.values(), i, &state.acceleration)
                    + row_action(
                        self.assembly.full_stiffness.values(),
                        i,
                        &state.displacement,
                    )
                    + self
                        .assembly
                        .full_damping
                        .as_ref()
                        .map_or(0.0, |c| row_action(c.values(), i, &state.velocity))
                    - load[i];
                if let Some(motion) = supports {
                    reaction[i] += self.support_row(
                        self.assembly.full_mass.values(),
                        i,
                        &motion.acceleration_m_s2,
                    ) + self.support_row(
                        self.assembly.full_stiffness.values(),
                        i,
                        &motion.displacement_m,
                    ) + self.assembly.full_damping.as_ref().map_or(0.0, |c| {
                        self.support_row(c.values(), i, &motion.velocity_m_s)
                    });
                }
            }
        }
        finite(&reaction)?;
        Ok(reaction)
    }
}

/// One accepted force interval, with actual support and damping exchanges.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingShellStep {
    /// Accepted endpoint time, s.
    pub time_s: f64,
    /// Applied load integrated by the native endpoint trapezoid, N s / N m s.
    pub applied_impulse: ShellSpatialResultant,
    /// Actual fixed or driven support reaction integrated by the same rule.
    pub support_impulse: ShellSpatialResultant,
    /// Signed impulse carried OUT of the shell by the declared damping matrix.
    /// Mass-proportional damping is a momentum sink, not an internal heat proof.
    pub damping_impulse: ShellSpatialResultant,
    /// Actual consistent-mass velocity change, about original reference nodes.
    pub momentum_change: ShellSpatialResultant,
    /// Change minus applied and support impulses, plus damping impulse.
    pub momentum_residual: ShellSpatialResultant,
    /// External endpoint-average force dotted with the native displacement increment, J.
    pub applied_work_j: f64,
    /// Endpoint-average support reaction dotted with prescribed motion, J.
    /// Zero for stationary supports. Positive work is supplied TO the shell.
    pub support_work_j: f64,
    /// `dt * v_mid^T C v_mid` in the original shell coordinates, J.
    pub damping_work_j: f64,
    /// Full mechanical-energy change plus damping work minus applied and support work, J.
    pub energy_residual_j: f64,
}

/// Continuing motion with immutable reference geometry/support membership and
/// cached native factors. Prescribed supports move only through explicit inputs.
pub struct MovingShellMotion<'a> {
    model: &'a MovingShellModel,
    integrator: PreparedNewmark<'a>,
    // The mass factorization already used at initialization also owns kicks.
    mass_solver: LuWorkspace,
    config: ShellTimeConfig,
    state: DynamicState,
    supports: Option<ShellSupportMotion>,
    full_load: Vec<f64>,
    reactions: Vec<f64>,
}
impl MovingShellMotion<'_> {
    /// Accepted native free-DOF motion. Energy fields include the complete
    /// prescribed and free motion for a driven-support trajectory.
    pub fn state(&self) -> &DynamicState {
        &self.state
    }
    /// Full endpoint support wrench ON the shell, zero on every free row.
    pub fn support_reactions(&self) -> &[f64] {
        &self.reactions
    }
    /// Current material point motion, including restrained components.
    pub fn observe(&self, point: ShellLoadPoint) -> Result<ShellPointMotion, MovingLoadError> {
        self.model.observe_with_supports(
            point,
            &self.state.displacement,
            &self.state.velocity,
            self.supports.as_ref(),
        )
    }

    /// Advance one sample with explicitly supplied endpoint load locations.
    /// `time_s` must equal the next native time exactly: stale, duplicate and
    /// future samples refuse. The old endpoint force is retained, so motion of
    /// the load integrates the correct force AND reference moment impulse.
    /// Cancellation/refusal retains the entire old state, load and reactions.
    /// A driven-support trajectory must instead use [`Self::advance_driven`].
    pub fn advance(
        &mut self,
        time_s: f64,
        loads: &[ShellPointWrench],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellStep, MovingLoadError> {
        self.advance_with_supports(time_s, loads, None, cancelled)
    }

    fn advance_with_supports(
        &mut self,
        time_s: f64,
        loads: &[ShellPointWrench],
        acceleration: Option<[[f64; 3]; 3]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingShellStep, MovingLoadError> {
        poll(&mut cancelled)?;
        let dt = self.config.dt_s;
        let expected = self.state.time_s + dt;
        if !time_s.is_finite()
            || time_s != expected
            || time_s <= self.state.time_s
            || time_s > self.config.duration_s
            || self.state.step >= self.config.max_steps
        {
            return Err(invalid(
                "load sample must be the next time within the admitted horizon",
            ));
        }
        let next_supports = match (&self.supports, acceleration) {
            (None, None) => None,
            (Some(old), Some(acceleration)) => Some(old.advance(dt, acceleration)?),
            _ => {
                return Err(invalid(
                    "driven support motion needs an explicit endpoint acceleration",
                ));
            }
        };
        let load = self.model.project(loads, &mut cancelled)?;
        let effective_load = self
            .model
            .effective_load(load.reduced, next_supports.as_ref())?;
        poll(&mut cancelled)?;
        let mut next = self.integrator.step(&self.state, &effective_load)?;
        if let Some(motion) = &next_supports {
            self.model.refresh_full_energy(&mut next, motion)?;
        }
        poll(&mut cancelled)?;
        let reactions = self
            .model
            .reactions(&next, next_supports.as_ref(), &load.full)?;
        let mean_velocity: Vec<_> = self
            .state
            .velocity
            .iter()
            .zip(&next.velocity)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let delta_velocity: Vec<_> = next
            .velocity
            .iter()
            .zip(&self.state.velocity)
            .map(|(a, b)| a - b)
            .collect();
        let mean_support_velocity = self.supports.zip(next_supports).map(|(old, next)| {
            std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    f64::midpoint(old.velocity_m_s[i][j], next.velocity_m_s[i][j])
                })
            })
        });
        let delta_support_velocity = self.supports.zip(next_supports).map(|(old, next)| {
            std::array::from_fn(|i| {
                std::array::from_fn(|j| next.velocity_m_s[i][j] - old.velocity_m_s[i][j])
            })
        });
        let mean_full = self
            .model
            .expand_with_supports(&mean_velocity, mean_support_velocity.as_ref());
        let damping = self
            .model
            .assembly
            .full_damping
            .as_ref()
            .map(|c| c.apply(&mean_full));
        // These are views of the accepted endpoints, not four additional
        // full-sized vectors allocated merely for impulse/work accounting.
        let applied = |i| dt * f64::midpoint(self.full_load[i], load.full[i]);
        let support = |i| dt * f64::midpoint(self.reactions[i], reactions[i]);
        let damping_increment = |i| dt * damping.as_ref().map_or(0.0, |c| c[i]);
        let change = self.model.assembly.full_mass.apply(
            &self
                .model
                .expand_with_supports(&delta_velocity, delta_support_velocity.as_ref()),
        );
        let momentum_change = self.model.resultant(&change)?;
        let applied_impulse = self.model.resultant_by(applied)?;
        let support_impulse = self.model.resultant_by(support)?;
        let damping_impulse = self.model.resultant_by(damping_increment)?;
        let momentum_residual = self
            .model
            .resultant_by(|i| change[i] - applied(i) - support(i) + damping_increment(i))?;
        let applied_work_j = mean_full
            .iter()
            .enumerate()
            .map(|(i, v)| applied(i) * v)
            .sum();
        let damping_work_j = mean_full
            .iter()
            .enumerate()
            .map(|(i, v)| damping_increment(i) * v)
            .sum();
        let support_work_j = if self.supports.is_some() {
            mean_full
                .iter()
                .enumerate()
                .map(|(i, v)| support(i) * v)
                .sum()
        } else {
            0.0
        };
        let energy_residual_j = next.total_energy_j - self.state.total_energy_j + damping_work_j
            - applied_work_j
            - support_work_j;
        finite(&[
            applied_work_j,
            support_work_j,
            damping_work_j,
            energy_residual_j,
        ])?;
        poll(&mut cancelled)?;
        self.state = next;
        self.supports = next_supports;
        self.full_load = load.full;
        self.reactions = reactions;
        Ok(MovingShellStep {
            time_s,
            applied_impulse,
            support_impulse,
            damping_impulse,
            momentum_change,
            momentum_residual,
            applied_work_j,
            support_work_j,
            damping_work_j,
            energy_residual_j,
        })
    }
}

fn invalid(message: &'static str) -> MovingLoadError {
    MovingLoadError::Invalid(message)
}
fn linear(source: fs_la::LuWorkspaceError) -> MovingLoadError {
    MovingLoadError::Time(ShellTimeError::LinearSolve { source })
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MovingLoadError> {
    if cancelled() {
        Err(MovingLoadError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(values: &[f64]) -> Result<(), MovingLoadError> {
    if values.iter().any(|v| !v.is_finite()) {
        Err(invalid("nonfinite shell load or motion"))
    } else {
        Ok(())
    }
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(value: [f64; 3]) -> Result<[f64; 3], MovingLoadError> {
    let norm = value[0].hypot(value[1]).hypot(value[2]);
    if !norm.is_finite() || norm == 0.0 {
        return Err(invalid("shell normal is unrepresentable"));
    }
    Ok(value.map(|x| x / norm))
}
