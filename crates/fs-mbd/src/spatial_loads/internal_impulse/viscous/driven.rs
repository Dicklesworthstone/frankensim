//! Exact simultaneous internal actuation and viscous relaxation at fixed poses.
//!
//! Let B be the original transported wrench map and G = B^T M^-1 B its
//! six-coordinate relative mobility. For relative twist s and held internal
//! effort f, s' = G (f - C s). The conjugate coordinate x = G^-1 s has
//! H = x^T G x / 2 and x' = f - C G x. Thus the existing quadratic pHS owner
//! integrates the complete source, its actuator work and its independent heat.
//! The actual bodies receive B (x_end - x_start) through the original native
//! impulse owner. Their common momentum and complete inertia remain intact.
//!
//! This source avoids the terminal-slip bias of separate force kicks and
//! unforced damping. Geometry is held during this source interval; it does not
//! advance either pose or claim an exact changing-geometry trajectory.

use super::{
    InternalWrenchImpulse, InternalWrenchImpulseError,
    InternalWrenchImpulseReceipt, Matrix, PassiveBushingMatrix, PreparedViscousConnector,
    ViscousConnectorControls, ViscousConnectorError, apply_internal_wrench_impulse, finite,
    kinetic_product, matrix_vector, poll, rotate_wide, rounded, rounded_six, validate_controls,
};
use crate::{MassProperties, Pose, RigidBodyState, Vec3};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::dd::Dd;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, ResistanceWitness};

/// Accepted native momentum exchange, held-actuator work and physical heat.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use]
pub struct DrivenViscousConnectorReceipt {
    pub impulse: InternalWrenchImpulseReceipt,
    /// Original relative twist in A's material axes: [m/s; 3], [rad/s; 3].
    pub relative_twist_before: [f64; 6],
    pub relative_twist_after: [f64; 6],
    /// Integral of the actual source twist divided by its duration. At zero
    /// duration this is the instantaneous initial twist.
    pub average_relative_twist: [f64; 6],
    /// Native endpoint twist minus the simultaneous continuous-source solution.
    pub relative_twist_residual: [f64; 6],
    /// Signed integral of the original held effort paired with relative twist.
    /// Negative work is actuator energy recovery, never negative heat [J].
    pub actuator_work_j: f64,
    /// Independently integrated s^T C s, using the original passive form [J].
    pub dissipated_energy_j: f64,
    /// Actual native kinetic change + heat - actuator work [J].
    pub kinetic_balance_residual_j: f64,
    /// Original native net-impulse work + heat - actuator work [J].
    pub connector_balance_residual_j: f64,
}

/// No refused source changes either body or publishes a partial result.
#[derive(Clone, Debug, PartialEq)]
pub enum DrivenViscousConnectorError {
    Source(ViscousConnectorError),
    Flow(PhsError),
    Unresolved {
        linear_speed_residual_m_s: f64,
        angular_speed_residual_rad_s: f64,
        kinetic_balance_residual_j: f64,
        connector_balance_residual_j: f64,
    },
}
impl From<ViscousConnectorError> for DrivenViscousConnectorError {
    fn from(value: ViscousConnectorError) -> Self {
        Self::Source(value)
    }
}
impl From<InternalWrenchImpulseError> for DrivenViscousConnectorError {
    fn from(value: InternalWrenchImpulseError) -> Self {
        Self::Source(value.into())
    }
}
impl From<PhsError> for DrivenViscousConnectorError {
    fn from(value: PhsError) -> Self {
        Self::Flow(value)
    }
}
impl core::fmt::Display for DrivenViscousConnectorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven spatial viscous connector refused: {self:?}")
    }
}
impl std::error::Error for DrivenViscousConnectorError {}

/// Full finite-body spatial damper and held internal actuator on one trajectory.
///
/// The authored symmetric PSD C acts on the same transported relative twist
/// as [`PreparedViscousConnector`]. Both original masses, full inertias,
/// material frames and attachment arms supply the mobility. All six relative
/// coordinates remain, including undamped directions and translation/spin
/// cross-coupling. No inverse damping, modal truncation or timestep fit is used.
///
/// Preparation factors one six-coordinate mobility and prepares the existing
/// bounded pHS exponential/heat operators. Warm calls perform one six-coordinate
/// triangular solve, cached matrix operations and one original native impulse;
/// they allocate nothing and evaluate no exponentials. Only internal scratch
/// changes between calls. Reprepare after geometry, inertia, C or interval
/// changes; initial momenta and the held actuator effort may change freely.
pub struct PreparedDrivenViscousConnector {
    geometry: PreparedViscousConnector,
    mobility: Matrix,
    mobility_factor: Cholesky,
    flow: PreparedExactZoh<'static>,
}

impl PreparedDrivenViscousConnector {
    pub fn new(
        poses: [Pose; 2],
        properties: [MassProperties; 2],
        arms_body_m: [Vec3; 2],
        damping: PassiveBushingMatrix,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DrivenViscousConnectorError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(ViscousConnectorError::Invalid("finite nonnegative duration").into());
        }
        // Zero source duration obtains only the original native geometry; this
        // adapter reuses the pHS propagator instead of preparing a second spectrum.
        let geometry = PreparedViscousConnector::new(
            poses,
            properties,
            arms_body_m,
            damping,
            0.0,
            &mut cancelled,
        )?;
        let columns = geometry.kinetic_columns()?;
        let mut mobility = [[0.0; 6]; 6];
        for i in 0..6 {
            poll(&mut cancelled)?;
            for j in i..6 {
                let entry = finite(kinetic_product(columns[i], columns[j]).to_f64())?;
                mobility[i][j] = entry;
                mobility[j][i] = entry;
            }
        }
        let q: Vec<_> = mobility.iter().flatten().copied().collect();
        let mobility_factor = cholesky(&q, 6).map_err(|_| {
            ViscousConnectorError::Invalid("resolved positive finite-body relative mobility")
        })?;
        let c: Vec<_> = geometry.damping.native_coefficients().iter().flatten().copied().collect();
        let factor = geometry.damping.gram_factor().map_err(ViscousConnectorError::from)?;
        // The passive-form owner supplies columns W with C = W W^T. The pHS
        // source witness uses rows W^T, retaining exactly missing loss channels.
        let rows = (0..6).flat_map(|i| (0..6).map(move |j| factor[j][i])).collect();
        let witness = ResistanceWitness::try_from_source_rows(6, c.clone(), rows)?;
        let mut ports = vec![0.0; 36];
        for i in 0..6 {
            ports[6 * i + i] = 1.0;
        }
        let system = PortHamiltonian::new(
            6,
            6,
            vec![0.0; 36],
            c,
            ports,
            Box::new(QuadraticStorage::new(q, 6)?),
        )?
        .with_resistance_witness(witness)?;
        poll(&mut cancelled)?;
        let flow = PreparedExactZoh::new_owned(
            system,
            duration_s,
            LinearSplitBudget { max_modes: 6, max_omega_dt: 1.0 },
        )?;
        poll(&mut cancelled)?;
        Ok(Self { geometry, mobility, mobility_factor, flow })
    }

    pub fn duration_s(&self) -> f64 {
        self.flow.duration()
    }

    /// Apply simultaneous held internal effort and passive damping.
    ///
    /// `effort` is [force N; 3], [free couple N m; 3] on attachment B in A's
    /// prepared material axes. The transported opposite wrench acts on A.
    /// Both poses stay fixed during the source; this duration is not a drift
    /// clock. Native momentum/work admission, initial-coordinate reconstruction,
    /// actual endpoint twist and continuous work/heat all gate the result.
    pub fn advance(
        &mut self,
        states: [RigidBodyState; 2],
        effort: [f64; 6],
        controls: ViscousConnectorControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenViscousConnectorReceipt, DrivenViscousConnectorError> {
        poll(&mut cancelled)?;
        validate_controls(controls.impulse)?;
        if effort.iter().any(|x| !x.is_finite())
            || [controls.linear_speed_tolerance_m_s, controls.angular_speed_tolerance_rad_s]
                .iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(ViscousConnectorError::Invalid("finite effort and nonnegative twist tolerances").into());
        }
        if states.map(|body| body.pose()) != self.geometry.poses {
            return Err(ViscousConnectorError::Invalid("changed prepared body poses").into());
        }
        let before = self.geometry.relative_twist(states)?;
        let mut initial = rounded_six(before)?;
        self.mobility_factor.solve(&mut initial);
        let represented = matrix_vector(self.mobility, initial)?;
        let initial_residual = rounded_six(core::array::from_fn(|i| {
            Dd::from_f64(represented[i]) - before[i]
        }))?;
        admit_twist(initial_residual, controls, 0.0, 0.0)?;
        let mut next = [0.0; 6];
        let mut average = [0.0; 6];
        let record = self.flow.step_into(&initial, &effort, &mut next, &mut average)?;
        poll(&mut cancelled)?;
        let source: [_; 6] = core::array::from_fn(|i| Dd::from_f64(next[i]) - Dd::from_f64(initial[i]));
        let impulse = apply_internal_wrench_impulse(
            states,
            self.geometry.properties,
            InternalWrenchImpulse {
                arms_body_m: self.geometry.arms_body_m,
                impulse_on_second_world_ns: rounded(rotate_wide(
                    self.geometry.axes_world,
                    [source[0], source[1], source[2]],
                ))?,
                couple_on_second_world_nms: rounded(rotate_wide(
                    self.geometry.axes_world,
                    [source[3], source[4], source[5]],
                ))?,
            },
            controls.impulse,
            &mut cancelled,
        )?;
        let after = self.geometry.relative_twist(impulse.bodies.map(|body| body.state_after))?;
        let expected = matrix_vector(self.mobility, next)?;
        let residual = rounded_six(core::array::from_fn(|i| after[i] - Dd::from_f64(expected[i])))?;
        let net_work = Dd::from_f64(record.supplied) - Dd::from_f64(record.dissipated);
        let kinetic_balance = finite((Dd::from_f64(impulse.kinetic_energy_change_j) - net_work).to_f64())?;
        let connector_balance = finite((Dd::from_f64(impulse.connector_work_j) - net_work).to_f64())?;
        admit_twist(residual, controls, kinetic_balance, connector_balance)?;
        let result = DrivenViscousConnectorReceipt {
            impulse,
            relative_twist_before: rounded_six(before)?,
            relative_twist_after: rounded_six(after)?,
            average_relative_twist: average,
            relative_twist_residual: residual,
            actuator_work_j: record.supplied,
            dissipated_energy_j: record.dissipated,
            kinetic_balance_residual_j: kinetic_balance,
            connector_balance_residual_j: connector_balance,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn admit_twist(
    residual: [f64; 6],
    controls: ViscousConnectorControls,
    kinetic_balance_residual_j: f64,
    connector_balance_residual_j: f64,
) -> Result<(), DrivenViscousConnectorError> {
    let linear_speed_residual_m_s = residual[..3].iter().map(|x| x.abs()).fold(0.0, f64::max);
    let angular_speed_residual_rad_s = residual[3..].iter().map(|x| x.abs()).fold(0.0, f64::max);
    if linear_speed_residual_m_s > controls.linear_speed_tolerance_m_s
        || angular_speed_residual_rad_s > controls.angular_speed_tolerance_rad_s
        || kinetic_balance_residual_j.abs() > controls.impulse.work_tolerance_j
        || connector_balance_residual_j.abs() > controls.impulse.work_tolerance_j
    {
        return Err(DrivenViscousConnectorError::Unresolved {
            linear_speed_residual_m_s,
            angular_speed_residual_rad_s,
            kinetic_balance_residual_j,
            connector_balance_residual_j,
        });
    }
    Ok(())
}
