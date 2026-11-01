//! Fluid inertia attached to the original articulated body's moving links.
//!
//! A strip is constant in its link frame. The native articulated owner supplies
//! the link twist `V = J u` and acceleration bias `c = Jdot u`; its same map
//! serves dry and fluid inertia. The existing positive strip modes give a
//! constant local metric `G`, so the fluid contribution is `Jᵀ G J` and its
//! right-hand-side bias is `−Jᵀ (G c + ad*_V G V)`. This includes serial moving
//! hinges, sliding/helical joints, all control cross blocks and their analytic
//! configuration terms without differentiating assembled matrices.
//!
//! Only this adapter permutes native `[angular, linear, joints]` columns into
//! Flyer's `[linear, angular, controls]` order. Cross-surface fluid interaction,
//! deforming strips and density time derivatives remain outside this model.
//!
//! [`assemble_articulated_strips_in_uniform_flow`] extends the same modes to a
//! prescribed spatially uniform world flow. It retains the relative momentum,
//! rotating-frame current derivative and source work, rather than treating
//! ground speed as fluid-relative speed. Its storage diagnostics concern smooth
//! motion; existing impact receipts do not account for this affine momentum.

use super::{
    AeroGeneralizedLoads, MAX_CONTROL_COORDS, MAX_STRIPS, Strip, admit_strip, refuse, strip_modes,
};
use crate::Refusal;
use fs_ga::{Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::GeneralizedKinematics;

#[derive(Clone, Copy, Debug, PartialEq)]
struct LocalMode {
    inertia_per_density: f64,
    covector: Wrench,
}

/// A finite strip rigidly attached to one native articulated link.
#[derive(Clone, Debug, PartialEq)]
pub struct LinkedStrip {
    link: usize,
    modes: [LocalMode; 3],
}

impl LinkedStrip {
    /// Bind strip geometry expressed in the link's own coordinates.
    ///
    /// Control fields must be inactive: all control velocities come from the
    /// articulated model, including inherited motion of upstream joints. Link
    /// existence is checked against the shared kinematics at assembly time.
    ///
    /// # Errors
    /// Invalid strip geometry or any preassigned control velocity.
    pub fn new(link: usize, strip_in_link: Strip) -> Result<Self, Refusal> {
        admit_strip(&strip_in_link, 0)?;
        if strip_in_link.control_coord.is_some()
            || strip_in_link.control_gain != 0.0
            || strip_in_link.control_angular_gain != [0.0; 3]
        {
            return Err(refuse(
                "linked-strip-invalid",
                "linked strip has preassigned control velocity".into(),
                "supply inactive control fields; the native link map owns all joint motion",
            ));
        }
        let modes = strip_modes(1.0, &strip_in_link, 0)?.modes.map(|mode| {
            let j = mode.rigid;
            LocalMode {
                inertia_per_density: mode.inertia,
                covector: Wrench::new(Vec3::new(j[3], j[4], j[5]), Vec3::new(j[0], j[1], j[2])),
            }
        });
        Ok(Self { link, modes })
    }

    /// Index in the native model's root-first link order.
    #[must_use]
    pub fn link(&self) -> usize {
        self.link
    }
}

/// Prescribed spatially uniform flow in an inertial world frame.
///
/// Velocity and its time derivative must describe the same instant. Spatial
/// shear, fluid rotation and density derivatives are not represented. The
/// pressure-gradient force on displaced fluid (Froude–Krylov force) is separate
/// from added mass: author that force as an external link wrench when required.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct UniformFlow {
    /// World translational fluid velocity [m/s].
    pub velocity_world_mps: Vec3,
    /// Time derivative of that world velocity [m/s²].
    pub acceleration_world_mps2: Vec3,
}

#[derive(Clone, Debug, Default, PartialEq)]
struct FlowDiagnostics {
    energy: f64,
    momentum: [f64; 6 + MAX_CONTROL_COORDS],
    homogeneous_momentum: [f64; 6 + MAX_CONTROL_COORDS],
    storage_bias_rate: f64,
    bias_body_power: f64,
    work_scale: f64,
}

/// Fluid blocks and immutable smooth-motion diagnostics at one native state.
///
/// The added-mass metric is independent of uniform flow, but momentum is affine
/// in generalized velocity: `p_fluid = Σ Jᵀ G (V − U_local)`. Keeping the blocks
/// private binds reaction and power calculations to the same admitted state.
/// Clone [`Self::loads`] for the existing smooth load callback; altering that
/// clone cannot alter this snapshot or its diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformFlowAssembly {
    loads: AeroGeneralizedLoads,
    flow: UniformFlow,
    base_axes_world: [Vec3; 3],
    diagnostics: FlowDiagnostics,
}

/// Added-mass reaction and independently checked instantaneous work balance.
#[derive(Clone, Debug, PartialEq)]
pub struct UniformFlowReaction {
    /// Fluid force/torque/joint effort in Flyer `[linear, angular, joints]` order.
    pub generalized_force: Vec<f64>,
    /// Sum of added-mass forces on all links, in world axes [N].
    pub force_world: Vec3,
    /// Mechanical power delivered by fluid to the actual body motion [W].
    pub body_power_w: f64,
    /// Power supplied by the prescribed ambient flow, `U_world · F_world` [W].
    pub ambient_power_w: f64,
    /// Independently reconstructed disturbance kinetic-energy derivative [W].
    pub storage_rate_w: f64,
    /// `storage_rate + body_power − ambient_power` [W].
    pub balance_residual_w: f64,
    /// Arithmetic-resolution tolerance used to admit the work balance [W].
    pub balance_tolerance_w: f64,
}

impl UniformFlowAssembly {
    /// Complete original fluid load blocks, with applied-load fields zero.
    #[must_use]
    pub fn loads(&self) -> &AeroGeneralizedLoads {
        &self.loads
    }

    /// Consume the snapshot to supply the existing articulated smooth solver.
    #[must_use]
    pub fn into_loads(self) -> AeroGeneralizedLoads {
        self.loads
    }

    /// Disturbance kinetic energy `Σ (V − U_local)ᵀ G (V − U_local) / 2` [J].
    ///
    /// This is not the homogeneous `uᵀ M_added u / 2` used by existing impact
    /// receipts when the ambient flow is nonzero.
    #[must_use]
    pub fn disturbance_energy_j(&self) -> f64 {
        self.diagnostics.energy
    }

    /// Affine fluid momentum in Flyer `[linear, angular, joints]` order.
    #[must_use]
    pub fn generalized_fluid_momentum(&self) -> &[f64] {
        &self.diagnostics.momentum[..self.loads.q_added_bias.len()]
    }

    /// Evaluate the fluid reaction for the actual generalized acceleration.
    ///
    /// Supply body-coordinate velocity derivatives in Flyer order, followed by
    /// joint accelerations, at this snapshot's instant. The reaction is
    /// `Q_fluid = bias − M_added * acceleration`. Body power is `Q_fluid · u`,
    /// evaluated using the original link modes; storage rate is independently
    /// `p_fluid · acceleration + Σ inertia * relative_speed * relative_bias`.
    /// Their sum must equal the prescribed flow's source power. No pressure,
    /// buoyancy, drag, impact, or density-change work is included.
    ///
    /// # Errors
    /// Wrong acceleration shape, nonfinite input or result, or a work balance
    /// outside its arithmetic-resolution tolerance.
    pub fn reaction(&self, acceleration: &[f64]) -> Result<UniformFlowReaction, Refusal> {
        let n = self.loads.q_added_bias.len();
        if acceleration.len() != n || !acceleration.iter().all(|a| a.is_finite()) {
            return Err(refuse(
                "uniform-flow-acceleration-invalid",
                format!("expected {n} finite Flyer-order acceleration components"),
                "supply the actual acceleration at this assembly's native state",
            ));
        }
        let mut generalized_force = self.loads.q_added_bias.clone();
        let mut storage_rate_w = self.diagnostics.storage_bias_rate;
        let mut body_power_w = self.diagnostics.bias_body_power;
        let mut work_scale = self.diagnostics.work_scale;
        for (i, &a) in acceleration.iter().enumerate() {
            let storage = self.diagnostics.momentum[i] * a;
            let body = self.diagnostics.homogeneous_momentum[i] * a;
            storage_rate_w += storage;
            body_power_w -= body;
            work_scale += storage.abs() + body.abs();
        }
        for (i, force) in generalized_force.iter_mut().enumerate() {
            // The original bounded blocks are multiplied directly. There is no
            // second metric, map, factorization, or reconstructed velocity.
            let mut force_scale = force.abs();
            for (j, &a) in acceleration.iter().enumerate() {
                let mass = match (i < 6, j < 6) {
                    (true, true) => self.loads.m_added_rr[i][j],
                    (true, false) => self.loads.m_added_rc[j - 6][i],
                    (false, true) => self.loads.m_added_rc[i - 6][j],
                    (false, false) => self.loads.m_added_cc[i - 6][j - 6],
                };
                let term = mass * a;
                *force -= term;
                force_scale += term.abs();
            }
            if i < 3 {
                work_scale += self
                    .flow
                    .velocity_world_mps
                    .dot(self.base_axes_world[i])
                    .abs()
                    * force_scale;
            }
        }
        let force_world = self.base_axes_world[0].scale(generalized_force[0])
            + self.base_axes_world[1].scale(generalized_force[1])
            + self.base_axes_world[2].scale(generalized_force[2]);
        let ambient_power_w = self.flow.velocity_world_mps.dot(force_world);
        let balance_residual_w = storage_rate_w + body_power_w - ambient_power_w;
        let balance_tolerance_w = 2048.0 * f64::EPSILON * (1.0 + work_scale);
        if !generalized_force.iter().all(|v| v.is_finite())
            || !finite_vector(force_world)
            || ![
                storage_rate_w,
                body_power_w,
                ambient_power_w,
                balance_residual_w,
                balance_tolerance_w,
            ]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(refuse(
                "uniform-flow-reaction-invalid",
                "unrepresentable fluid reaction or work rate".into(),
                "use representable velocities, accelerations, geometry and density",
            ));
        }
        if balance_residual_w.abs() > balance_tolerance_w {
            return Err(refuse(
                "uniform-flow-work-balance",
                format!("work residual {balance_residual_w:e} W exceeds {balance_tolerance_w:e} W"),
                "use a representable native state and acceleration at the same instant",
            ));
        }
        Ok(UniformFlowReaction {
            generalized_force,
            force_world,
            body_power_w,
            ambient_power_w,
            storage_rate_w,
            balance_residual_w,
            balance_tolerance_w,
        })
    }
}

fn finite_vector(v: Vec3) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

/// Assemble the same linked-strip law in a prescribed uniform moving fluid.
///
/// In each link frame, `W = V − (0, Rᵀ U)` and
/// `c_relative = c + (0, omega × Rᵀ U − Rᵀ Udot)`. The complete bias is
/// `−Jᵀ (G c_relative + ad*_W G W)`: replacing only strip normal speed would
/// miss the rotating-current acceleration and translational transport. The
/// original rank modes, native kinematics and mass blocks are shared exactly.
///
/// Feed the returned blocks into the existing smooth articulated load callback
/// at each evaluated instant. Fluid density is held constant. Displaced-volume
/// pressure-gradient loads remain explicit external wrenches; this API does not
/// infer a displaced volume from strip geometry. Its affine momentum and work
/// diagnostics do not extend the energy claims of existing impact solvers.
///
/// # Errors
/// The original strip/density/map domain, nonfinite flow, or unrepresentable
/// transport, momentum, or disturbance energy.
pub fn assemble_articulated_strips_in_uniform_flow(
    rho_kg_m3: f64,
    strips: &[LinkedStrip],
    kinematics: &GeneralizedKinematics,
    flow: UniformFlow,
) -> Result<UniformFlowAssembly, Refusal> {
    if !finite_vector(flow.velocity_world_mps) || !finite_vector(flow.acceleration_world_mps2) {
        return Err(refuse(
            "uniform-flow-invalid",
            "world flow velocity and acceleration must be finite".into(),
            "supply a spatially uniform flow and its world time derivative",
        ));
    }
    let mut diagnostics = FlowDiagnostics::default();
    let loads = assemble(
        rho_kg_m3,
        strips,
        kinematics,
        Some(flow),
        Some(&mut diagnostics),
    )?;
    if ![
        diagnostics.energy,
        diagnostics.storage_bias_rate,
        diagnostics.bias_body_power,
        diagnostics.work_scale,
    ]
    .iter()
    .chain(&diagnostics.momentum)
    .chain(&diagnostics.homogeneous_momentum)
    .all(|v| v.is_finite())
    {
        return Err(refuse(
            "uniform-flow-invalid",
            "unrepresentable disturbance energy, fluid momentum or work rate".into(),
            "use representable flow, geometry, density and native motion",
        ));
    }
    // Native root translation columns map base translation into root-link
    // axes. This remains correct for a nonidentity root mounting transform.
    let root_rotation = kinematics.world_from_link(0).rotation();
    let root_columns = kinematics.body_jacobian(0);
    let mut base_axes_world = [Vec3::default(); 3];
    for i in 0..3 {
        base_axes_world[i] = root_rotation
            .rotate(root_columns[3 + i].linear)
            .map_err(|_| {
                refuse(
                    "uniform-flow-invalid",
                    "unrepresentable native world translation axes".into(),
                    "supply admitted finite native transforms",
                )
            })?;
    }
    Ok(UniformFlowAssembly {
        loads,
        flow,
        base_axes_world,
        diagnostics,
    })
}

fn native_column(flyer_column: usize) -> usize {
    if flyer_column < 6 {
        (flyer_column + 3) % 6
    } else {
        flyer_column
    }
}

/// Assemble fluid blocks and their full bias using already admitted kinematics.
///
/// Fixed/root links and arbitrarily nested scalar joints use the same strip
/// law. The output's applied-load fields are zero; its bias contains only fluid
/// inertia. Sum it with the complete dry generalized inertia and dry bias from
/// the same native state before solving; adding a constant body and independent
/// scalar control inertias would omit the moving dry linkage's cross terms.
///
/// The native map is prepared once and reused. Each strip's admitted local
/// rank modes are cached at construction, then scaled by density and projected
/// using one bounded stack row. There is no per-strip matrix or numerical
/// configuration derivative. Density is held constant in the bias calculation.
///
/// # Errors
/// Nonpositive/nonfinite density, control or strip caps, missing links, or
/// unrepresentable assembled mass/bias.
pub fn assemble_articulated_strips(
    rho_kg_m3: f64,
    strips: &[LinkedStrip],
    kinematics: &GeneralizedKinematics,
) -> Result<AeroGeneralizedLoads, Refusal> {
    assemble(rho_kg_m3, strips, kinematics, None, None)
}

fn assemble(
    rho_kg_m3: f64,
    strips: &[LinkedStrip],
    kinematics: &GeneralizedKinematics,
    flow: Option<UniformFlow>,
    mut diagnostics: Option<&mut FlowDiagnostics>,
) -> Result<AeroGeneralizedLoads, Refusal> {
    let nc = kinematics.dof_count();
    if nc > MAX_CONTROL_COORDS {
        return Err(refuse(
            "control-coords-outside-domain",
            format!("nc = {nc} exceeds the cap {MAX_CONTROL_COORDS}"),
            "use at most the admitted number of active dynamic coordinates",
        ));
    }
    if strips.len() > MAX_STRIPS {
        return Err(refuse(
            "strip-count-exceeded",
            format!("{} strips exceed {MAX_STRIPS}", strips.len()),
            "aggregate panels into strips",
        ));
    }
    if !rho_kg_m3.is_finite() || rho_kg_m3 <= 0.0 {
        return Err(refuse(
            "non-finite-input",
            "fluid density must be finite and positive".into(),
            "supply the admitted fluid density",
        ));
    }
    if strips.iter().any(|s| s.link >= kinematics.link_count()) {
        return Err(refuse(
            "linked-strip-invalid",
            "strip refers to a missing articulated link".into(),
            "bind strips to links in the same native articulated model",
        ));
    }

    let mut result = AeroGeneralizedLoads {
        q_rigid_nonaccel: [0.0; 6],
        q_control_nonaccel: vec![0.0; nc],
        m_added_rr: [[0.0; 6]; 6],
        m_added_rc: vec![[0.0; 6]; nc],
        m_added_cc: vec![vec![0.0; nc]; nc],
        q_added_bias: vec![0.0; 6 + nc],
    };
    // The still/zero-flow path keeps the original arithmetic, including signed
    // zeros, and skips all frame transforms. Diagnostics are optional.
    let moving_flow = flow.filter(|f| *f != UniformFlow::default());
    for strip in strips {
        let velocity = kinematics.body_twist(strip.link);
        let acceleration = kinematics.body_bias_acceleration(strip.link);
        let jacobian = kinematics.body_jacobian(strip.link);
        let (relative_velocity, relative_acceleration) = if let Some(flow) = moving_flow {
            let world_to_link = kinematics.world_from_link(strip.link).rotation().inverse();
            let local_flow = world_to_link.rotate(flow.velocity_world_mps).map_err(|_| {
                refuse(
                    "uniform-flow-invalid",
                    "unrepresentable local flow velocity".into(),
                    "supply representable world flow and native transforms",
                )
            })?;
            let local_acceleration =
                world_to_link
                    .rotate(flow.acceleration_world_mps2)
                    .map_err(|_| {
                        refuse(
                            "uniform-flow-invalid",
                            "unrepresentable local flow acceleration".into(),
                            "supply representable world flow and native transforms",
                        )
                    })?;
            (
                Twist::new(velocity.angular, velocity.linear - local_flow),
                acceleration.plus(Twist::new(
                    Vec3::default(),
                    velocity.angular.cross(local_flow) - local_acceleration,
                )),
            )
        } else {
            (velocity, acceleration)
        };
        for mode in strip.modes {
            let inertia = rho_kg_m3 * mode.inertia_per_density;
            if !inertia.is_finite() || inertia <= 0.0 {
                return Err(refuse(
                    "linked-strip-invalid",
                    "unrepresentable strip inertia at the supplied density".into(),
                    "use representable density and strip dimensions",
                ));
            }
            let covector = mode.covector;
            let speed = covector.pairing(relative_velocity);
            let bias_speed = covector.pairing(relative_acceleration);
            let mut row = [0.0; 6 + MAX_CONTROL_COORDS];
            for a in 0..6 + nc {
                let column = jacobian[native_column(a)];
                row[a] = covector.pairing(column);
                // <J_a, ad*_W p> = -<p, [W, J_a]>: use the original
                // Lie bracket to include angular AND linear transport once.
                result.q_added_bias[a] += inertia
                    * (speed * covector.pairing(relative_velocity.bracket(column))
                        - row[a] * bias_speed);
            }
            if let Some(d) = diagnostics.as_deref_mut() {
                let absolute_speed = covector.pairing(velocity);
                let storage = inertia * (speed * bias_speed);
                let transport =
                    inertia * (speed * covector.pairing(relative_velocity.bracket(velocity)));
                let body = inertia * (absolute_speed * bias_speed);
                d.energy += 0.5 * inertia * speed * speed;
                d.storage_bias_rate += storage;
                d.bias_body_power += transport - body;
                d.work_scale += storage.abs() + transport.abs() + body.abs();
                for a in 0..6 + nc {
                    d.momentum[a] += inertia * (speed * row[a]);
                    d.homogeneous_momentum[a] += inertia * (absolute_speed * row[a]);
                }
            }
            for a in 0..6 {
                for b in 0..6 {
                    result.m_added_rr[a][b] += inertia * (row[a] * row[b]);
                }
            }
            for a in 0..nc {
                for b in 0..6 {
                    result.m_added_rc[a][b] += inertia * (row[6 + a] * row[b]);
                }
                for b in 0..nc {
                    result.m_added_cc[a][b] += inertia * (row[6 + a] * row[6 + b]);
                }
            }
        }
    }
    if !result
        .m_added_rr
        .iter()
        .flatten()
        .chain(result.m_added_rc.iter().flatten())
        .chain(result.m_added_cc.iter().flatten())
        .chain(&result.q_added_bias)
        .all(|x| x.is_finite())
    {
        return Err(refuse(
            "linked-strip-invalid",
            "unrepresentable fluid inertia or transport bias".into(),
            "use representable geometry, density and joint motion",
        ));
    }
    Ok(result)
}
