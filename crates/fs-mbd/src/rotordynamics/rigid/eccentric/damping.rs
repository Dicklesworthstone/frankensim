//! Passive stationary-bearing relaxation on the original eccentric rigid body.

use super::{
    ElasticBearing, ElasticEccentricRotor, ElasticEccentricRotorError, ElasticEccentricRotorStep,
};
use crate::articulated::SpatialInertia;
use crate::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, BushingMatrixError, PassiveBushingMatrix,
};
use crate::rotordynamics::rigid::{MAX_ROTOR_STATIONS, RadialBearing};
use crate::rotordynamics::transient::{RotorTransientState, RotorTransientStep};
use crate::rotordynamics::{RotorBearingModel, RotorError};
use crate::{
    DynamicsError, MomentumImpulseReceipt, Pose, PreparedRigidPointMotion, RigidBodyIntegrator,
    RigidBodyState, RigidKineticEnergyChange, Vec3, Wrench,
};
use fs_ga::Se3;
use fs_math::{dd::Dd, det};

/// A finite-rotation eccentric body on original passive radial K/C bearings.
///
/// The elastic component retains its original native geometry and full inertia.
/// Each supplied damping matrix is admitted by the existing passive material
/// owner and retained exactly, including anisotropy and symmetric cross terms.
/// Bearing velocity is that of the actual material shaft point relative to its
/// stationary support. Its force and reciprocal COM moment act together.
/// Active/nonsymmetric damping, moving supports, nonlinear coefficient histories,
/// fluid films and a thermal destination are not inferred.
#[derive(Clone, Debug)]
pub struct DampedEccentricRotor {
    elastic: ElasticEccentricRotor,
    damping: Vec<[[f64; 2]; 2]>,
    source: RotorBearingModel,
    active: bool,
}

/// Explicit source resolution, native transfer and complete-step energy limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EccentricRotorDampingControl {
    /// Original full-mass midpoint impulse-equation relative tolerance in (0,1).
    pub relative_tolerance: f64,
    /// Maximum source-duration times a conservative damping-rate trace bound.
    /// Must lie in (0,2]; the upper cap excludes midpoint's negative modal
    /// amplification on unresolved stiff decay. This is not an accuracy estimate.
    pub maximum_modal_step: f64,
    /// Maximum componentwise actual-versus-constitutive force impulse defect [N s].
    pub linear_impulse_tolerance_ns: f64,
    /// Maximum componentwise actual-versus-constitutive COM moment defect [N m s].
    pub angular_impulse_tolerance_nms: f64,
    /// Maximum actual-versus-native-source endpoint COM velocity defect [m/s].
    pub linear_speed_tolerance_m_s: f64,
    /// Maximum actual-versus-native-source endpoint body angular-rate defect [rad/s].
    pub angular_speed_tolerance_rad_s: f64,
    /// Maximum source kinetic-plus-heat and impulse-work-plus-heat defects [J].
    pub source_energy_tolerance_j: f64,
    /// Maximum full-step storage-plus-heat-minus-external-work defect [J].
    /// Elastic trajectory truncation remains numerical error, never bearing heat.
    pub total_energy_tolerance_j: f64,
}

/// Actual elastic/damping wrench and separate instantaneous storage/loss.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DampedEccentricRotorLoads {
    /// Complete bearing wrench plus supplied shaft torque, excluding gravity.
    pub wrench: Wrench,
    /// Original elastic bearing storage [J].
    pub bearing_potential_j: f64,
    /// Original nonnegative sum of bearing velocity^T C velocity [W].
    pub dissipation_power_w: f64,
}

/// One fixed-pose damping source and its independently evaluated physical heat.
#[derive(Clone, Debug, PartialEq)]
pub struct EccentricBearingDampingStep {
    /// One original full-body momentum update; pose is unchanged.
    pub impulse: MomentumImpulseReceipt,
    /// Original six-coordinate solve in [world COM velocity, body angular rate].
    /// Its displacement is a source-work quadrature, not an extra body pose.
    /// Zero damping/duration skips the solve and returns None.
    pub native: Option<RotorTransientStep>,
    /// Source interval [s]; no physical pose or clock advances in this substep.
    pub duration_seconds: f64,
    /// Conservative original-inertia trace bound on the largest decay rate [1/s].
    /// Zero duration skips preparation and does not measure a physical zero rate.
    pub decay_rate_bound_per_s: f64,
    /// Midpoint heat from each original bearing, in supplied order [J].
    pub bearing_heat_j: Vec<f64>,
    /// Sum of the constitutive midpoint heats [J], not a kinetic-energy difference.
    pub dissipated_energy_j: f64,
    /// Actual original-body kinetic change plus independently evaluated heat [J].
    pub kinetic_heat_residual_j: f64,
    /// Original commanded-impulse work plus independently evaluated heat [J].
    pub impulse_heat_residual_j: f64,
    /// Actual momentum increment minus the law at actual midpoint rates [N s].
    pub linear_impulse_residual_world: Vec3,
    /// Same constitutive check on the original body-frame angular momentum [N m s].
    pub angular_impulse_residual_body: Vec3,
}

/// Symmetric damping/elastic/damping motion with complete original storage.
#[derive(Clone, Debug, PartialEq)]
pub struct DampedEccentricRotorStep {
    /// Sole physical endpoint after both damping sources and elastic motion.
    pub state_after: RigidBodyState,
    /// Sources prepared at the actual initial and final native poses.
    pub sources: [EccentricBearingDampingStep; 2],
    /// Original configuration-force Verlet receipt under held shaft torque.
    pub elastic: ElasticEccentricRotorStep,
    /// Complete initial-to-final native full-body kinetic change [J].
    pub kinetic_energy_change_j: f64,
    /// Sum of the independently evaluated source heats [J].
    pub dissipated_energy_j: f64,
    /// delta(K + gravity potential + bearing potential) + heat - shaft work [J].
    pub total_energy_residual_j: f64,
}

/// Refusal leaves the original body, coefficient law and geometry unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum EccentricRotorDampingError {
    /// Invalid numerical budget or unrepresentable derived quantity.
    Invalid(&'static str),
    /// Original elastic/full-body refusal.
    Mechanical(ElasticEccentricRotorError),
    /// Original passive damping admission refused.
    Material(BushingMatrixError),
    /// Original full-mass midpoint source refused.
    Native(RotorError),
    /// Damping is too stiff for the requested source interval; refine time.
    ModalStep { measure: f64, limit: f64 },
    /// Actual momentum, rate or work does not resolve the original source law.
    SourceBalance {
        linear_impulse_residual_ns: f64,
        angular_impulse_residual_nms: f64,
        linear_speed_residual_m_s: f64,
        angular_speed_residual_rad_s: f64,
        kinetic_heat_residual_j: f64,
        impulse_heat_residual_j: f64,
    },
    /// Complete source/storage balance exceeds the authored numerical allowance.
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
    /// Cancellation observed before complete publication.
    Cancelled,
}
impl From<DynamicsError> for EccentricRotorDampingError {
    fn from(e: DynamicsError) -> Self {
        Self::Mechanical(e.into())
    }
}
impl From<ElasticEccentricRotorError> for EccentricRotorDampingError {
    fn from(e: ElasticEccentricRotorError) -> Self {
        Self::Mechanical(e)
    }
}
impl From<BushingMatrixError> for EccentricRotorDampingError {
    fn from(e: BushingMatrixError) -> Self {
        Self::Material(e)
    }
}
impl From<RotorError> for EccentricRotorDampingError {
    fn from(e: RotorError) -> Self {
        Self::Native(e)
    }
}
impl core::fmt::Display for EccentricRotorDampingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "eccentric bearing damping refused: {self:?}")
    }
}
impl std::error::Error for EccentricRotorDampingError {}

impl DampedEccentricRotor {
    /// Prepare the original material K/C coefficients and full COM inertia.
    ///
    /// The reference transform and stations have the same meaning as for
    /// ElasticEccentricRotor. K is admitted by that elastic owner. C must pass
    /// native exact-symmetry/passive-matrix admission; no damping floor,
    /// coefficient projection, fictitious housing body or inferred mode is added.
    pub fn new(
        source: SpatialInertia,
        world_from_shaft: Se3,
        bearings: &[RadialBearing],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, EccentricRotorDampingError> {
        poll(&mut cancelled)?;
        if bearings.len() > MAX_ROTOR_STATIONS {
            return Err(EccentricRotorDampingError::Invalid(
                "too many rotor bearings",
            ));
        }
        let mut springs = Vec::with_capacity(bearings.len());
        let mut damping = Vec::with_capacity(bearings.len());
        let mut active = false;
        for bearing in bearings {
            poll(&mut cancelled)?;
            let c = bearing.damping_ns_per_m;
            let mut native = [[0.0; 6]; 6];
            for i in 0..2 {
                for j in 0..2 {
                    native[i][j] = c[i][j];
                }
            }
            PassiveBushingMatrix::try_new(
                native,
                BushingCoordinateScales {
                    translation_m: 1.0,
                    rotation_rad: 1.0,
                },
            )?;
            active |= c.into_iter().flatten().any(|x| x != 0.0);
            damping.push(c);
            springs.push(RadialBearing {
                damping_ns_per_m: [[0.0; 2]; 2],
                ..*bearing
            });
        }
        let elastic =
            ElasticEccentricRotor::new(source, world_from_shaft, &springs, &mut cancelled)?;
        let properties = elastic.mass_properties();
        let tensor = properties.inertia_tensor_body();
        let mut mass = [0.0; 36];
        for i in 0..3 {
            mass[6 * i + i] = properties.mass();
            for j in 0..3 {
                mass[6 * (i + 3) + j + 3] = tensor.m[3 * i + j];
            }
        }
        let source =
            RotorBearingModel::new(6, &mass, &[0.0; 36], &[0.0; 36], &[0.0; 36], &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            elastic,
            damping,
            source,
            active,
        })
    }

    /// Original elastic geometry/inertia component. Stepping it alone excludes C.
    #[must_use]
    pub const fn elastic_component(&self) -> &ElasticEccentricRotor {
        &self.elastic
    }

    /// Original unmodified radial damping coefficients for one supplied bearing.
    #[must_use]
    pub fn damping_coefficients(&self, bearing: usize) -> Option<[[f64; 2]; 2]> {
        self.damping.get(bearing).copied()
    }

    /// Actual velocity-dependent bearing load, including reciprocal shaft torque.
    pub fn loads(
        &self,
        state: RigidBodyState,
        shaft_torque_nm: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DampedEccentricRotorLoads, EccentricRotorDampingError> {
        let elastic = self.elastic.loads(state, shaft_torque_nm, &mut cancelled)?;
        let (force, heat) = self.source_law(
            state.pose(),
            rates(state, &self.elastic)?,
            1.0,
            &mut cancelled,
        )?;
        let f = elastic.wrench.force_world;
        let t = elastic.wrench.torque_body;
        let original = [f.x, f.y, f.z, t.x, t.y, t.z];
        let total: [f64; 6] =
            core::array::from_fn(|i| (Dd::from_f64(original[i]) + force[i]).to_f64());
        for x in total {
            finite(x)?;
        }
        let power = finite(heat.into_iter().fold(Dd::ZERO, |sum, x| sum + x).to_f64())?;
        poll(&mut cancelled)?;
        Ok(DampedEccentricRotorLoads {
            wrench: Wrench {
                force_world: vector(&total[..3]),
                torque_body: vector(&total[3..]),
            },
            bearing_potential_j: elastic.bearing_potential_j,
            dissipation_power_w: power,
        })
    }

    /// Original full-inertia instantaneous acceleration, with actual damping force.
    pub fn prepare_motion(
        &self,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        shaft_torque_nm: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRigidPointMotion, EccentricRotorDampingError> {
        let loads = self.loads(state, shaft_torque_nm, cancelled)?;
        Ok(integrator.prepare_point_motion(state, self.elastic.mass_properties(), loads.wrench)?)
    }

    /// One simultaneous passive damping source at the original fixed pose.
    ///
    /// With v=[v_COM_world, omega_body], each bearing has rows
    /// J=[e_world, arm_body cross R^-1 e_world]. Native midpoint solves
    /// (M+h D/2) delta_v = -h D v with D=sum J^T C J. One combined native
    /// impulse applies all reactions. Heat is independently h*sum(u_mid^T C u_mid)
    /// at actual returned body midpoint rates, never the measured kinetic loss.
    /// Zero duration/damping is an identity; no pose or clock advances here.
    pub fn relax_damping(
        &self,
        state: RigidBodyState,
        duration_seconds: f64,
        control: EccentricRotorDampingControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<EccentricBearingDampingStep, EccentricRotorDampingError> {
        poll(&mut cancelled)?;
        validate(control)?;
        let h = duration_seconds;
        if !h.is_finite() || h < 0.0 {
            return Err(EccentricRotorDampingError::Invalid(
                "finite nonnegative source duration required",
            ));
        }
        let before = rates(state, &self.elastic)?;
        let mut predicted = before;
        let mut native = None;
        let mut rate_bound = 0.0;
        if self.active && h != 0.0 {
            let (damping, bound) = self.assemble(state.pose(), &mut cancelled)?;
            rate_bound = bound;
            let measure = finite((Dd::from_f64(h) * Dd::from_f64(bound)).to_f64())?;
            if measure > control.maximum_modal_step {
                return Err(EccentricRotorDampingError::ModalStep {
                    measure,
                    limit: control.maximum_modal_step,
                });
            }
            // Reuse the once-admitted original mass, zero K and zero G. Only
            // the finite passive geometric congruence changes with pose.
            let mut model = self.source.clone();
            model.damping.copy_from_slice(&damping);
            let operating = model.at_speed(0.0, &mut cancelled)?;
            let step = operating
                .midpoint_stepper(h, control.relative_tolerance, &mut cancelled)?
                .step(
                    &RotorTransientState {
                        displacement: vec![0.0; 6],
                        velocity: before.to_vec(),
                    },
                    &[0.0; 6],
                    &mut cancelled,
                )?;
            predicted.copy_from_slice(&step.state_after.velocity);
            native = Some(step);
        }
        let predicted_mid = core::array::from_fn(|i| 0.5 * before[i] + 0.5 * predicted[i]);
        let (command, _) = self.source_law(state.pose(), predicted_mid, h, &mut cancelled)?;
        let command = command.map(|x| x.to_f64());
        for x in command {
            finite(x)?;
        }
        let impulse = state.apply_momentum_impulse(
            self.elastic.mass_properties(),
            vector(&command[..3]),
            vector(&command[3..]),
        )?;
        let after = rates(impulse.state_after, &self.elastic)?;
        let actual_mid = core::array::from_fn(|i| 0.5 * before[i] + 0.5 * after[i]);
        let (actual_law, heat) = self.source_law(state.pose(), actual_mid, h, &mut cancelled)?;
        let actual_before = momenta(state);
        let actual_after = momenta(impulse.state_after);
        let residual: [f64; 6] = core::array::from_fn(|i| {
            (Dd::from_f64(actual_after[i]) - Dd::from_f64(actual_before[i]) - actual_law[i])
                .to_f64()
        });
        for x in residual {
            finite(x)?;
        }
        let energy = heat.iter().copied().fold(Dd::ZERO, |sum, x| sum + x);
        let kinetic_residual =
            finite((Dd::from_f64(impulse.work.kinetic_energy_change_j) + energy).to_f64())?;
        let work_residual = finite((Dd::from_f64(impulse.work.impulse_work_j) + energy).to_f64())?;
        let speed_error: [f64; 6] = core::array::from_fn(|i| after[i] - predicted[i]);
        for x in speed_error {
            finite(x)?;
        }
        let (linear_impulse, angular_impulse, linear_speed, angular_speed) = (
            maximum(&residual[..3]),
            maximum(&residual[3..]),
            maximum(&speed_error[..3]),
            maximum(&speed_error[3..]),
        );
        if linear_impulse > control.linear_impulse_tolerance_ns
            || angular_impulse > control.angular_impulse_tolerance_nms
            || linear_speed > control.linear_speed_tolerance_m_s
            || angular_speed > control.angular_speed_tolerance_rad_s
            || kinetic_residual.abs() > control.source_energy_tolerance_j
            || work_residual.abs() > control.source_energy_tolerance_j
        {
            return Err(EccentricRotorDampingError::SourceBalance {
                linear_impulse_residual_ns: linear_impulse,
                angular_impulse_residual_nms: angular_impulse,
                linear_speed_residual_m_s: linear_speed,
                angular_speed_residual_rad_s: angular_speed,
                kinetic_heat_residual_j: kinetic_residual,
                impulse_heat_residual_j: work_residual,
            });
        }
        let result = EccentricBearingDampingStep {
            impulse,
            native,
            duration_seconds: h,
            decay_rate_bound_per_s: rate_bound,
            bearing_heat_j: heat.iter().map(|x| x.to_f64()).collect(),
            dissipated_energy_j: finite(energy.to_f64())?,
            kinetic_heat_residual_j: kinetic_residual,
            impulse_heat_residual_j: work_residual,
            linear_impulse_residual_world: vector(&residual[..3]),
            angular_impulse_residual_body: vector(&residual[3..]),
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Advance D(h/2), original elastic Verlet(h), D(h/2) at actual endpoint poses.
    ///
    /// Symmetric composition is second order for smooth fixed coefficients and
    /// resolved elastic/free rotation. Midpoint is not L-stable; both sources
    /// enforce their modal-resolution gate. Gyroscopic acceleration belongs to
    /// the native free-body motion, not a damping term. No velocity load is
    /// passed to configuration-only Verlet. Shaft/gravity work and physical
    /// bearing heat stay separate; a complete-step defect may require refinement.
    pub fn step(
        &self,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        shaft_torque_nm: f64,
        duration_seconds: f64,
        control: EccentricRotorDampingControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DampedEccentricRotorStep, EccentricRotorDampingError> {
        let half = 0.5 * duration_seconds;
        if !half.is_finite() || half <= 0.0 || 2.0 * half != duration_seconds {
            return Err(EccentricRotorDampingError::Invalid(
                "representable positive half interval required",
            ));
        }
        let first = self.relax_damping(state, half, control, &mut cancelled)?;
        let elastic = self.elastic.step(
            integrator,
            first.impulse.state_after,
            shaft_torque_nm,
            duration_seconds,
            &mut cancelled,
        )?;
        let last = self.relax_damping(elastic.native.state_after, half, control, &mut cancelled)?;
        let state_after = last.impulse.state_after;
        let mut kinetic = RigidKineticEnergyChange::new();
        kinetic.add_body(self.elastic.mass_properties(), state, state_after)?;
        let kinetic = kinetic.finish_retained()?;
        let heat = Dd::from_f64(first.dissipated_energy_j) + Dd::from_f64(last.dissipated_energy_j);
        let residual = finite(
            (kinetic
                + Dd::from_f64(elastic.native.work.gravitational_potential_energy_change_j)
                + Dd::from_f64(elastic.bearing_potential_change_j)
                + heat
                - Dd::from_f64(elastic.shaft_work_j))
            .to_f64(),
        )?;
        if residual.abs() > control.total_energy_tolerance_j {
            return Err(EccentricRotorDampingError::EnergyBalance {
                residual_j: residual,
                tolerance_j: control.total_energy_tolerance_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(DampedEccentricRotorStep {
            state_after,
            sources: [first, last],
            elastic,
            kinetic_energy_change_j: finite(kinetic.to_f64())?,
            dissipated_energy_j: finite(heat.to_f64())?,
            total_energy_residual_j: residual,
        })
    }

    fn jacobian(&self, bearing: &ElasticBearing, axes_body: [Vec3; 2]) -> [[f64; 6]; 2] {
        core::array::from_fn(|i| {
            let f = self.elastic.radial_axes_world[i];
            let t = bearing.arm_body.cross(axes_body[i]);
            [f.x, f.y, f.z, t.x, t.y, t.z]
        })
    }

    fn assemble(
        &self,
        pose: Pose,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; 36], f64), EccentricRotorDampingError> {
        let axes = self
            .elastic
            .radial_axes_world
            .map(|v| pose.orientation().rotate_world_to_body(v));
        let properties = self.elastic.mass_properties();
        let principal = properties.principal_axes_to_body();
        let moments = properties.principal_inertia_body();
        let mass_root = det::sqrt(properties.mass());
        let inertia_roots = [moments.x, moments.y, moments.z].map(det::sqrt);
        let mut matrix = [Dd::ZERO; 36];
        let mut trace = Dd::ZERO;
        for (bearing, c) in self.elastic.bearings.iter().zip(&self.damping) {
            poll(cancelled)?;
            let j = self.jacobian(bearing, axes);
            for row in 0..6 {
                for col in row..6 {
                    for a in 0..2 {
                        for b in 0..2 {
                            matrix[6 * row + col] = matrix[6 * row + col]
                                + Dd::from_f64(j[a][row])
                                    * Dd::from_f64(c[a][b])
                                    * Dd::from_f64(j[b][col]);
                        }
                    }
                }
            }
            // trace(M^-1 D) bounds every PSD decay eigenvalue. Evaluate it as
            // original C forms of mass-whitened Jacobian columns, using the
            // native principal axes only for the full inverse-inertia action.
            let whitened = j.map(|row| {
                let angular = principal.rotate_world_to_body(vector(&row[3..]));
                [
                    row[0] / mass_root,
                    row[1] / mass_root,
                    row[2] / mass_root,
                    angular.x / inertia_roots[0],
                    angular.y / inertia_roots[1],
                    angular.z / inertia_roots[2],
                ]
            });
            for k in 0..6 {
                trace = trace + quadratic(*c, [whitened[0][k], whitened[1][k]])?;
            }
        }
        let mut result = [0.0; 36];
        for i in 0..6 {
            for j in i..6 {
                let value = finite(matrix[6 * i + j].to_f64())?;
                result[6 * i + j] = value;
                result[6 * j + i] = value;
            }
        }
        // Leave a numerical margin for the native coordinate maps; this is a
        // resolution gate, not an interval spectral certificate.
        let bound = finite((trace * Dd::from_f64(1.0 + 256.0 * f64::EPSILON)).to_f64())?;
        Ok((result, bound))
    }

    fn source_law(
        &self,
        pose: Pose,
        velocity: [f64; 6],
        h: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([Dd; 6], Vec<Dd>), EccentricRotorDampingError> {
        let axes = self
            .elastic
            .radial_axes_world
            .map(|v| pose.orientation().rotate_world_to_body(v));
        let mut impulse = [Dd::ZERO; 6];
        let mut heat = Vec::with_capacity(self.damping.len());
        for (bearing, c) in self.elastic.bearings.iter().zip(&self.damping) {
            poll(cancelled)?;
            let j = self.jacobian(bearing, axes);
            let u: [f64; 2] = core::array::from_fn(|a| {
                (0..6)
                    .fold(Dd::ZERO, |sum, k| {
                        sum + Dd::from_f64(j[a][k]) * Dd::from_f64(velocity[k])
                    })
                    .to_f64()
            });
            let force: [Dd; 2] = core::array::from_fn(|a| {
                -(Dd::from_f64(c[a][0]) * Dd::from_f64(u[0])
                    + Dd::from_f64(c[a][1]) * Dd::from_f64(u[1]))
                    * Dd::from_f64(h)
            });
            for k in 0..6 {
                for a in 0..2 {
                    impulse[k] = impulse[k] + Dd::from_f64(j[a][k]) * force[a];
                }
            }
            let loss = quadratic(*c, u)? * Dd::from_f64(h);
            finite(loss.to_f64())?;
            heat.push(loss);
        }
        for x in impulse {
            finite(x.to_f64())?;
        }
        Ok((impulse, heat))
    }
}

fn quadratic(c: [[f64; 2]; 2], u: [f64; 2]) -> Result<Dd, EccentricRotorDampingError> {
    for x in u {
        finite(x)?;
    }
    let mut value = Dd::ZERO;
    for i in 0..2 {
        for j in 0..2 {
            value = value + Dd::from_f64(u[i]) * Dd::from_f64(c[i][j]) * Dd::from_f64(u[j]);
        }
    }
    finite(value.to_f64())?;
    if value.lt(Dd::ZERO) {
        return Err(EccentricRotorDampingError::Invalid(
            "unresolved negative original damping quadratic",
        ));
    }
    Ok(value)
}
fn rates(
    state: RigidBodyState,
    rotor: &ElasticEccentricRotor,
) -> Result<[f64; 6], EccentricRotorDampingError> {
    let v = state.center_of_mass_velocity_world(rotor.mass_properties())?;
    let w = rotor
        .mass_properties()
        .angular_velocity_body_checked(state.angular_momentum_body())?;
    Ok([v.x, v.y, v.z, w.x, w.y, w.z])
}
fn momenta(state: RigidBodyState) -> [f64; 6] {
    let p = state.linear_momentum_world();
    let l = state.angular_momentum_body();
    [p.x, p.y, p.z, l.x, l.y, l.z]
}
fn vector(v: &[f64]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn maximum(v: &[f64]) -> f64 {
    v.iter().fold(0.0_f64, |a, x| a.max(x.abs()))
}
fn finite(x: f64) -> Result<f64, EccentricRotorDampingError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(EccentricRotorDampingError::Invalid(
            "nonfinite damping input or response",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), EccentricRotorDampingError> {
    if cancelled() {
        Err(EccentricRotorDampingError::Cancelled)
    } else {
        Ok(())
    }
}
fn validate(c: EccentricRotorDampingControl) -> Result<(), EccentricRotorDampingError> {
    if !c.relative_tolerance.is_finite()
        || c.relative_tolerance <= 0.0
        || c.relative_tolerance >= 1.0
        || !c.maximum_modal_step.is_finite()
        || c.maximum_modal_step <= 0.0
        || c.maximum_modal_step > 2.0
        || [
            c.linear_impulse_tolerance_ns,
            c.angular_impulse_tolerance_nms,
            c.linear_speed_tolerance_m_s,
            c.angular_speed_tolerance_rad_s,
            c.source_energy_tolerance_j,
            c.total_energy_tolerance_j,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
    {
        return Err(EccentricRotorDampingError::Invalid(
            "valid source resolution and nonnegative dimensional allowances required",
        ));
    }
    Ok(())
}
