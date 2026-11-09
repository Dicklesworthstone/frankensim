//! Finite-inertia shaft acceleration with reciprocal elastic-bearing reactions.

use super::{MAX_ROTOR_STATIONS, RadialBearing};
use crate::articulated::SpatialInertia;
use crate::{
    DynamicsError, MassProperties, Pose, PreparedRigidPointMotion, RigidBodyIntegrator,
    RigidBodyState, StepReceipt, UnitQuaternion, Vec3, Wrench,
};
use fs_ga::Se3;
use fs_math::dd::Dd;

mod damping;
pub use damping::{
    DampedEccentricRotor, DampedEccentricRotorLoads, DampedEccentricRotorStep,
    EccentricBearingDampingStep, EccentricRotorDampingControl, EccentricRotorDampingError,
};

/// A finite-rotation rigid rotor whose actual centre of mass may miss its shaft.
///
/// Body +z is the shaft axis. The supplied spatial inertia is about the actual
/// COM; its COM offset locates that COM from the authored shaft reference frame.
/// Bearing forces act on material shaft points, with their complete reciprocal
/// moment about the COM. No prescribed centrifugal load or extra unbalance mass
/// is added. Asymmetric full inertia, tilt and spin reversal remain native rigid
/// motion, not a four-coordinate nominal-rotor linearization.
///
/// Bearings are fixed-world, bilateral radial springs. Their coefficients must
/// be symmetric positive semidefinite with zero damping. Clearance, fluid-film
/// cross coupling, velocity damping, flexible shafts and electrical source
/// limits require their own constitutive owners.
#[derive(Clone, Debug)]
pub struct ElasticEccentricRotor {
    source: SpatialInertia,
    properties: MassProperties,
    reference_pose: Pose,
    radial_axes_world: [Vec3; 2],
    bearings: Vec<ElasticBearing>,
}

#[derive(Clone, Copy, Debug)]
struct ElasticBearing {
    arm_body: Vec3,
    anchor_world: Vec3,
    stiffness: [[f64; 2]; 2],
}

/// Instantaneous native load and elastic storage at one rotor configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticEccentricRotorLoads {
    /// Bearing forces and their COM moments, plus the supplied shaft torque.
    /// Configured native gravity is excluded and must not be added again.
    pub wrench: Wrench,
    /// Sum of `displacement^T K displacement / 2` over the bearings [J].
    pub bearing_potential_j: f64,
}

/// Accepted native Verlet step with independently measured source and storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElasticEccentricRotorStep {
    /// Original native state, work quadrature and gravity diagnostics. Its applied
    /// work includes bearing work; that work is not an additional source here.
    pub native: StepReceipt,
    /// Original bearing storage before the step [J].
    pub bearing_potential_before_j: f64,
    /// Original bearing storage after the step [J].
    pub bearing_potential_after_j: f64,
    /// Polarized bearing storage increment, without subtracting stored totals [J].
    pub bearing_potential_change_j: f64,
    /// Held body-z drive torque times endpoint-trapezoidal angular travel [J].
    /// This second-order quadrature is independent of the storage change.
    pub shaft_work_j: f64,
    /// `delta(K + gravity potential + bearing potential) - shaft_work` [J].
    /// Signed numerical defect, never assigned to heat or motor efficiency.
    pub balance_residual_j: f64,
}

/// Refusal before any caller-owned state is advanced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ElasticEccentricRotorError {
    /// Invalid or unsupported bearing/source data.
    Invalid(&'static str),
    /// Original full-inertia rigid-body refusal.
    Native(DynamicsError),
    /// Cooperative cancellation observed before publishing a receipt.
    Cancelled,
}

impl From<DynamicsError> for ElasticEccentricRotorError {
    fn from(value: DynamicsError) -> Self {
        Self::Native(value)
    }
}

impl core::fmt::Display for ElasticEccentricRotorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(message) => f.write_str(message),
            Self::Native(error) => core::fmt::Display::fmt(error, f),
            Self::Cancelled => f.write_str("eccentric rotor operation cancelled"),
        }
    }
}
impl std::error::Error for ElasticEccentricRotorError {}

impl ElasticEccentricRotor {
    /// Prepare the original COM inertia and material bearing geometry once.
    ///
    /// `world_from_shaft` fixes the stationary bearing centres and radial axes.
    /// Stations are measured along +z from this authored reference origin.
    /// At [`Self::reference_pose`] all bearing displacements are zero, including
    /// when the actual COM has radial or axial offset. A subsequent native state
    /// always locates the COM, not the shaft origin. Zero bearings are admitted.
    pub fn new(
        source: SpatialInertia,
        world_from_shaft: Se3,
        bearings: &[RadialBearing],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ElasticEccentricRotorError> {
        poll(&mut cancelled)?;
        if bearings.len() > MAX_ROTOR_STATIONS {
            return Err(ElasticEccentricRotorError::Invalid(
                "too many rotor bearings",
            ));
        }
        let properties =
            MassProperties::from_inertia_tensor(source.mass(), Vec3::ZERO, source.inertia_com())?;
        let rotation = world_from_shaft.rotation();
        let q = rotation.as_quat();
        let orientation = UnitQuaternion::new(q.w, q.x, q.y, q.z)?;
        let origin = native(world_from_shaft.translation());
        let com = native(source.center_of_mass());
        let reference_pose = Pose::new(
            origin.add(orientation.rotate_body_to_world(com)),
            orientation,
        )?;
        let radial_axes_world = [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)]
            .map(|axis| orientation.rotate_body_to_world(axis));
        let mut prepared = Vec::with_capacity(bearings.len());
        for bearing in bearings {
            poll(&mut cancelled)?;
            validate_bearing(bearing)?;
            let station = Vec3::new(0.0, 0.0, bearing.station_m);
            let arm_body = station.sub(com);
            let anchor_world = origin.add(orientation.rotate_body_to_world(station));
            if !arm_body.is_finite() || !anchor_world.is_finite() {
                return Err(ElasticEccentricRotorError::Invalid(
                    "unrepresentable bearing geometry",
                ));
            }
            prepared.push(ElasticBearing {
                arm_body,
                anchor_world,
                stiffness: bearing.stiffness_n_per_m,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            properties,
            reference_pose,
            radial_axes_world,
            bearings: prepared,
        })
    }

    /// Original mass, actual COM offset and full centroidal tensor.
    #[must_use]
    pub const fn source_inertia(&self) -> SpatialInertia {
        self.source
    }

    /// Native COM properties retaining the original full body tensor.
    #[must_use]
    pub const fn mass_properties(&self) -> MassProperties {
        self.properties
    }

    /// COM pose corresponding to the undeformed authored shaft configuration.
    #[must_use]
    pub const fn reference_pose(&self) -> Pose {
        self.reference_pose
    }

    /// Body-fixed shaft point relative to the actual COM, in metres.
    pub fn shaft_arm_body(&self, station_m: f64) -> Result<Vec3, ElasticEccentricRotorError> {
        let arm = Vec3::new(0.0, 0.0, station_m).sub(native(self.source.center_of_mass()));
        if !arm.is_finite() {
            return Err(ElasticEccentricRotorError::Invalid("invalid shaft station"));
        }
        Ok(arm)
    }

    /// Evaluate springs and a finite held torque about the actual body +z shaft.
    ///
    /// Force and moment use the same material point: `tau_body = arm cross
    /// R^-1 force_world`. Thus an eccentric lateral bearing force accelerates
    /// the spin as well as the COM; spin speed is never prescribed here.
    pub fn loads(
        &self,
        state: RigidBodyState,
        shaft_torque_nm: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElasticEccentricRotorLoads, ElasticEccentricRotorError> {
        poll(&mut cancelled)?;
        finite(shaft_torque_nm)?;
        let mut sum = [Dd::from_f64(0.0); 7];
        sum[5] = Dd::from_f64(shaft_torque_nm);
        for bearing in &self.bearings {
            poll(&mut cancelled)?;
            let d = self.displacement(state.pose(), bearing)?;
            let radial_force = multiply(bearing.stiffness, d).map(|f| -f);
            let force = self.radial_axes_world[0]
                .scale(radial_force[0])
                .add(self.radial_axes_world[1].scale(radial_force[1]));
            let torque = bearing
                .arm_body
                .cross(state.pose().orientation().rotate_world_to_body(force));
            for (i, value) in [force.x, force.y, force.z, torque.x, torque.y, torque.z]
                .into_iter()
                .enumerate()
            {
                sum[i] = sum[i] + Dd::from_f64(finite(value)?);
            }
            sum[6] = sum[6]
                - (Dd::from_f64(d[0]) * Dd::from_f64(radial_force[0])
                    + Dd::from_f64(d[1]) * Dd::from_f64(radial_force[1]))
                    * Dd::from_f64(0.5);
        }
        let values = sum.map(|x| x.to_f64());
        for value in values {
            finite(value)?;
        }
        poll(&mut cancelled)?;
        Ok(ElasticEccentricRotorLoads {
            wrench: Wrench {
                force_world: Vec3::new(values[0], values[1], values[2]),
                torque_body: Vec3::new(values[3], values[4], values[5]),
            },
            bearing_potential_j: values[6],
        })
    }

    /// Native full-inertia instantaneous acceleration, reusable at shaft points.
    pub fn prepare_motion(
        &self,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        shaft_torque_nm: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedRigidPointMotion, ElasticEccentricRotorError> {
        let loads = self.loads(state, shaft_torque_nm, cancelled)?;
        Ok(integrator.prepare_point_motion(state, self.properties, loads.wrench)?)
    }

    /// Advance actual COM and angular momenta through native configuration Verlet.
    ///
    /// The same elastic law is evaluated at both endpoint poses. The held torque
    /// follows body +z and may brake or reverse the shaft. Gravity comes only
    /// from `integrator`. The caller must resolve bearing frequencies and rotor
    /// motion by timestep refinement; no exact energy or unconditional stability
    /// is claimed. Failure/cancellation returns no candidate state.
    pub fn step(
        &self,
        integrator: RigidBodyIntegrator,
        state: RigidBodyState,
        shaft_torque_nm: f64,
        duration_seconds: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElasticEccentricRotorStep, ElasticEccentricRotorError> {
        let before = self.loads(state, shaft_torque_nm, &mut cancelled)?;
        let (native, after) = integrator.step_verlet_with(
            state,
            self.properties,
            before.wrench,
            duration_seconds,
            |endpoint| {
                let after = self.loads(endpoint, shaft_torque_nm, &mut cancelled)?;
                Ok::<_, ElasticEccentricRotorError>((after.wrench, after))
            },
        )?;
        let mut change = Dd::from_f64(0.0);
        for bearing in &self.bearings {
            poll(&mut cancelled)?;
            let d0 = self.displacement(state.pose(), bearing)?;
            let d1 = self.displacement(native.state_after.pose(), bearing)?;
            // Polarize in retained arithmetic so a small storage change is not
            // lost by subtracting two large elastic energies.
            let delta =
                core::array::from_fn::<_, 2, _>(|i| Dd::from_f64(d1[i]) - Dd::from_f64(d0[i]));
            let mean = core::array::from_fn::<_, 2, _>(|i| {
                (Dd::from_f64(d1[i]) + Dd::from_f64(d0[i])) * Dd::from_f64(0.5)
            });
            for i in 0..2 {
                for j in 0..2 {
                    change = change + mean[i] * Dd::from_f64(bearing.stiffness[i][j]) * delta[j];
                }
            }
        }
        let w0 = self
            .properties
            .angular_velocity_body_checked(state.angular_momentum_body())?
            .z;
        let w1 = self
            .properties
            .angular_velocity_body_checked(native.state_after.angular_momentum_body())?
            .z;
        let work = Dd::from_f64(shaft_torque_nm)
            * Dd::from_f64(duration_seconds)
            * (Dd::from_f64(w0) + Dd::from_f64(w1))
            * Dd::from_f64(0.5);
        let residual = Dd::from_f64(native.work.mechanical_energy_change_j) + change - work;
        let result = ElasticEccentricRotorStep {
            native,
            bearing_potential_before_j: before.bearing_potential_j,
            bearing_potential_after_j: after.bearing_potential_j,
            bearing_potential_change_j: finite(change.to_f64())?,
            shaft_work_j: finite(work.to_f64())?,
            balance_residual_j: finite(residual.to_f64())?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn displacement(
        &self,
        pose: Pose,
        bearing: &ElasticBearing,
    ) -> Result<[f64; 2], ElasticEccentricRotorError> {
        let offset = pose
            .point_world_from_body(bearing.arm_body)?
            .sub(bearing.anchor_world);
        let result = self.radial_axes_world.map(|axis| axis.dot(offset));
        for value in result {
            finite(value)?;
        }
        Ok(result)
    }
}

fn validate_bearing(bearing: &RadialBearing) -> Result<(), ElasticEccentricRotorError> {
    finite(bearing.station_m)?;
    for value in bearing.stiffness_n_per_m.into_iter().flatten() {
        finite(value)?;
    }
    if bearing
        .damping_ns_per_m
        .into_iter()
        .flatten()
        .any(|x| x != 0.0)
    {
        return Err(ElasticEccentricRotorError::Invalid(
            "elastic rotor requires zero bearing damping",
        ));
    }
    let [[a, b], [c, d]] = bearing.stiffness_n_per_m;
    if b != c || a < 0.0 || d < 0.0 {
        return Err(ElasticEccentricRotorError::Invalid(
            "elastic stiffness must be symmetric positive semidefinite",
        ));
    }
    if b != 0.0 {
        let scale = a.max(d).max(b.abs());
        let ad = Dd::from_f64(a / scale) * Dd::from_f64(d / scale);
        let bb = Dd::from_f64(b / scale) * Dd::from_f64(b / scale);
        if bb.to_f64() == 0.0 || (ad - bb).to_f64() < 0.0 {
            return Err(ElasticEccentricRotorError::Invalid(
                "elastic stiffness is indefinite or unresolved",
            ));
        }
    }
    Ok(())
}
fn native(v: fs_ga::Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, v.z)
}
fn multiply(k: [[f64; 2]; 2], d: [f64; 2]) -> [f64; 2] {
    k.map(|row| {
        (Dd::from_f64(row[0]) * Dd::from_f64(d[0]) + Dd::from_f64(row[1]) * Dd::from_f64(d[1]))
            .to_f64()
    })
}
fn finite(value: f64) -> Result<f64, ElasticEccentricRotorError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ElasticEccentricRotorError::Invalid(
            "nonfinite eccentric rotor input or response",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ElasticEccentricRotorError> {
    if cancelled() {
        Err(ElasticEccentricRotorError::Cancelled)
    } else {
        Ok(())
    }
}
