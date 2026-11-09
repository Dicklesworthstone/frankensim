//! Simultaneous held actuation and full spatial damping on a finite graph.
//!
//! The ORIGINAL native wrench columns are written in body kinetic coordinates:
//! y = [p/sqrt(m), I^(-1/2) h], H = y.y/2, and s = B^T y. Thus
//! y' = B f - B C B^T y. Every body, undamped direction, cycle, parallel link
//! and cross-coupled source remains in the existing exact quadratic pHS owner.
//! No inverse edge mobility, rank cutoff or sequence of pair relaxations occurs.
//! Original positive source rows integrate each link's heat along this SAME
//! trajectory. The net impulses then enter the native simultaneous wrench
//! owner once; its actual momentum, work and endpoint rates gate publication.

/// Complete material-frame Maxwell histories on the same native graph.
#[path = "driven/relaxation.rs"]
pub mod relaxation;

use super::super::{
    PreparedViscousConnector, ViscousConnectorControls, ViscousConnectorError,
    checked_vector_divide, finite, poll, rotate_wide, rounded, validate_controls,
};
use super::ViscousNetworkLink;
use crate::spatial_loads::internal_impulse::{
    InternalWrenchImpulse, InternalWrenchImpulseError,
    network::{
        InternalWrenchLink, InternalWrenchNetworkControls, InternalWrenchNetworkError, admit_limit,
        apply_internal_wrench_network, reserve,
    },
};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, Vec3};
use fs_math::{dd::Dd, det};
use fs_phs::{
    PhsError, PortHamiltonian, PreparedStepError, QuadraticStorage, ResistanceWitness,
    linear_exact::{
        DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        PreparedExactZoh, PreparedSourceLoss,
    },
};

/// Complete cold preparation limits for a simultaneous graph source.
///
/// The full body coordinate count is 6*bodies. Both that count doubled and
/// 6*(bodies+connectors)+1 must fit the original 256-coordinate integral owner.
/// The existing sparse passive network remains available for larger graphs.
/// Bounds cover the complete operator, all original source observers, scratch
/// and one native result, before allocation; no physical state is truncated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DrivenViscousNetworkLimits {
    pub maximum_bodies: usize,
    pub maximum_connectors: usize,
    /// Conservative logical peak payload, excluding allocator metadata [bytes].
    pub maximum_live_bytes: usize,
    /// Conservative dense preparation scalar multiply-add count, including all
    /// independently integrated connector heat channels.
    pub maximum_preparation_terms: u64,
}

/// Accepted original body bank and distinct source-work and heat accounts.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct DrivenViscousNetworkReceipt {
    pub states_after: Vec<RigidBodyState>,
    pub connector_actuator_work_j: Vec<f64>,
    pub connector_dissipation_j: Vec<f64>,
    /// Original link order, each in its first body's material frame.
    pub average_relative_twist: Vec<[f64; 6]>,
    pub source_duration_s: f64,
    pub actuator_work_j: f64,
    pub dissipated_energy_j: f64,
    pub kinetic_energy_change_j: f64,
    pub connector_work_j: f64,
    pub angular_origin_world_m: Vec3,
    pub linear_momentum_change_world_ns: Vec3,
    pub angular_momentum_change_world_nms: Vec3,
    pub kinetic_balance_residual_j: f64,
    pub connector_balance_residual_j: f64,
    pub heat_partition_residual_j: f64,
    pub actuator_partition_residual_j: f64,
    pub maximum_linear_speed_residual_m_s: f64,
    pub maximum_angular_speed_residual_rad_s: f64,
    pub planned_live_bytes: usize,
    pub planned_preparation_terms: u64,
}

/// Refusal publishes no intermediate body update or partial source account.
#[derive(Clone, Debug, PartialEq)]
pub enum DrivenViscousNetworkError {
    Invalid(&'static str),
    Network(InternalWrenchNetworkError),
    Connector {
        connector: usize,
        error: ViscousConnectorError,
    },
    Flow(PhsError),
    Unresolved {
        linear_speed_residual_m_s: f64,
        angular_speed_residual_rad_s: f64,
        kinetic_balance_residual_j: f64,
        connector_balance_residual_j: f64,
        heat_partition_residual_j: f64,
        actuator_partition_residual_j: f64,
    },
    Cancelled,
}
impl From<InternalWrenchNetworkError> for DrivenViscousNetworkError {
    fn from(error: InternalWrenchNetworkError) -> Self {
        match error {
            InternalWrenchNetworkError::Cancelled => Self::Cancelled,
            error => Self::Network(error),
        }
    }
}
impl From<InternalWrenchImpulseError> for DrivenViscousNetworkError {
    fn from(error: InternalWrenchImpulseError) -> Self {
        InternalWrenchNetworkError::from(error).into()
    }
}
impl From<DynamicsError> for DrivenViscousNetworkError {
    fn from(error: DynamicsError) -> Self {
        InternalWrenchNetworkError::from(error).into()
    }
}
impl From<PhsError> for DrivenViscousNetworkError {
    fn from(error: PhsError) -> Self {
        Self::Flow(error)
    }
}
impl From<PreparedStepError> for DrivenViscousNetworkError {
    fn from(error: PreparedStepError) -> Self {
        match error {
            PreparedStepError::Cancelled => Self::Cancelled,
            PreparedStepError::Solver(error) => Self::Flow(error),
        }
    }
}
impl core::fmt::Display for DrivenViscousNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven spatial viscous network refused: {self:?}")
    }
}
impl std::error::Error for DrivenViscousNetworkError {}

struct PreparedLink {
    bodies: [usize; 2],
    geometry: PreparedViscousConnector,
    loss: Option<PreparedSourceLoss>,
}
struct Scratch {
    initial: Vec<f64>,
    endpoint: Vec<f64>,
    input: Vec<f64>,
    average: Vec<f64>,
    heat_coordinates: Vec<f64>,
}

/// Fixed original graph, geometry, inertia, material laws and duration.
///
/// Preparation reuses the native finite-body wrench columns and existing exact
/// pHS/positive source-square integrals. Warm calls perform cached contractions
/// and one simultaneous native impulse network, with no exponential, eigensolve
/// or factorization. Returned body/link arrays and the original native event
/// allocate bounded O(bodies+connectors) storage; source scratch is retained.
///
/// Held efforts and momenta may change between calls. Reprepare after changing
/// a pose, mass, inertia, attachment, damping law or the source interval.
pub struct PreparedDrivenViscousNetwork {
    poses: Vec<Pose>,
    properties: Vec<MassProperties>,
    links: Vec<PreparedLink>,
    input_map: Vec<f64>,
    flow: Option<PreparedExactZoh<'static>>,
    scratch: Scratch,
    duration_s: f64,
    planned_live_bytes: usize,
    planned_preparation_terms: u64,
}

impl PreparedDrivenViscousNetwork {
    pub fn new(
        poses: &[Pose],
        properties: &[MassProperties],
        links: &[ViscousNetworkLink],
        duration_s: f64,
        limits: DrivenViscousNetworkLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DrivenViscousNetworkError> {
        poll(&mut cancelled)?;
        if poses.len() != properties.len() || !duration_s.is_finite() || duration_s < 0.0 {
            return Err(DrivenViscousNetworkError::Invalid(
                "matching original poses/inertias and finite nonnegative duration required",
            ));
        }
        admit_limit("bodies", poses.len() as u128, limits.maximum_bodies as u128)?;
        admit_limit(
            "connectors",
            links.len() as u128,
            limits.maximum_connectors as u128,
        )?;
        let (bytes, terms) = planned_costs(poses.len(), links.len(), duration_s)?;
        admit_limit("live bytes", bytes, limits.maximum_live_bytes as u128)?;
        admit_limit(
            "preparation terms",
            terms,
            limits.maximum_preparation_terms as u128,
        )?;
        for link in links {
            poll(&mut cancelled)?;
            let [a, b] = link.bodies;
            if a == b || a >= poses.len() || b >= poses.len() {
                return Err(DrivenViscousNetworkError::Invalid(
                    "every connector requires two distinct existing original bodies",
                ));
            }
        }
        let (n, m) = (6 * poses.len(), 6 * links.len());
        let GraphGeometry {
            poses: bound_poses,
            properties: bound_properties,
            links: mut prepared,
            input_map,
            rows,
            resistance,
        } = prepare_geometry(poses, properties, links, &mut cancelled)?;
        let flow = if duration_s == 0.0 || links.is_empty() {
            None
        } else {
            let r = resistance
                .into_iter()
                .map(|x| finite(x.to_f64()))
                .collect::<Result<Vec<_>, _>>()?;
            let generator: Vec<_> = r.iter().map(|x| -*x).collect();
            let witness = ResistanceWitness::try_from_source_rows(n, r.clone(), rows.clone())?;
            let mut q = vec![0.0; n * n];
            for i in 0..n {
                q[i * n + i] = 1.0;
            }
            let system = PortHamiltonian::new(
                n,
                m,
                vec![0.0; n * n],
                r,
                input_map.clone(),
                Box::new(QuadraticStorage::new(q, n)?),
            )?
            .with_resistance_witness(witness)?;
            poll(&mut cancelled)?;
            let flow = PreparedExactZoh::new_owned(
                system,
                duration_s,
                LinearSplitBudget {
                    max_modes: n,
                    max_omega_dt: 1.0,
                },
            )?;
            poll(&mut cancelled)?;
            let budget = DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: bytes as usize,
                max_multiply_terms: terms as u64,
            };
            for (e, link) in prepared.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                link.loss = Some(PreparedSourceLoss::new_controlled(
                    n,
                    m,
                    &generator,
                    &input_map,
                    &rows[6 * e * n..6 * (e + 1) * n],
                    duration_s,
                    budget,
                    &mut cancelled,
                )?);
            }
            Some(flow)
        };
        poll(&mut cancelled)?;
        Ok(Self {
            poses: bound_poses,
            properties: bound_properties,
            links: prepared,
            input_map,
            flow,
            scratch: Scratch {
                initial: vec![0.0; n],
                endpoint: vec![0.0; n],
                input: vec![0.0; m],
                average: vec![0.0; m],
                heat_coordinates: vec![0.0; 2 * n],
            },
            duration_s,
            planned_live_bytes: bytes as usize,
            planned_preparation_terms: terms as u64,
        })
    }

    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn planned_live_bytes(&self) -> usize {
        self.planned_live_bytes
    }
    pub const fn planned_preparation_terms(&self) -> u64 {
        self.planned_preparation_terms
    }

    /// Simultaneously advance every held internal effort and original damper.
    /// Each effort is [force N;3, free couple N m;3] on the second attachment,
    /// in the first body's prepared material axes. All inputs remain borrowed.
    pub fn advance(
        &mut self,
        states: &[RigidBodyState],
        efforts: &[[f64; 6]],
        controls: ViscousConnectorControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenViscousNetworkReceipt, DrivenViscousNetworkError> {
        poll(&mut cancelled)?;
        validate_controls(controls.impulse)?;
        if states.len() != self.poses.len()
            || states.iter().zip(&self.poses).any(|(s, p)| s.pose() != *p)
            || efforts.len() != self.links.len()
            || efforts.iter().flatten().any(|x| !x.is_finite())
            || [
                controls.linear_speed_tolerance_m_s,
                controls.angular_speed_tolerance_rad_s,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(DrivenViscousNetworkError::Invalid(
                "original poses, complete finite efforts and nonnegative twist tolerances required",
            ));
        }
        let m = 6 * self.links.len();
        for (i, state) in states.iter().enumerate() {
            poll(&mut cancelled)?;
            let p = checked_vector_divide(
                state.linear_momentum_world(),
                det::sqrt(self.properties[i].mass()),
            )
            .map_err(|_| {
                DrivenViscousNetworkError::Invalid(
                    "original body kinetic coordinate is unrepresentable",
                )
            })?;
            let h = self.properties[i].angular_mobility_factor(state.angular_momentum_body());
            let values = [p.x, p.y, p.z, h.x, h.y, h.z];
            if values.iter().any(|x| !x.is_finite()) {
                return Err(DrivenViscousNetworkError::Invalid(
                    "finite original kinetic coordinates required",
                ));
            }
            self.scratch.initial[6 * i..6 * (i + 1)].copy_from_slice(&values);
        }
        for (e, effort) in efforts.iter().enumerate() {
            self.scratch.input[6 * e..6 * (e + 1)].copy_from_slice(effort);
        }
        let mut linear_error = 0.0_f64;
        let mut angular_error = 0.0_f64;
        for (e, link) in self.links.iter().enumerate() {
            poll(&mut cancelled)?;
            let native = link
                .geometry
                .relative_twist(link.bodies.map(|i| states[i]))
                .map_err(|error| connector_error(e, error))?;
            for c in 0..6 {
                let expected = port_rate(&self.input_map, m, 6 * e + c, &self.scratch.initial);
                let error = finite((native[c] - expected).to_f64())?.abs();
                if c < 3 {
                    linear_error = linear_error.max(error);
                } else {
                    angular_error = angular_error.max(error);
                }
            }
        }
        admit(controls, linear_error, angular_error, [0.0; 4])?;
        let (flow_work, flow_heat) = if let Some(flow) = &mut self.flow {
            let record = flow.step_into(
                &self.scratch.initial,
                &self.scratch.input,
                &mut self.scratch.endpoint,
                &mut self.scratch.average,
            )?;
            (record.supplied, record.dissipated)
        } else {
            self.scratch.endpoint.copy_from_slice(&self.scratch.initial);
            for p in 0..m {
                self.scratch.average[p] =
                    finite(port_rate(&self.input_map, m, p, &self.scratch.initial).to_f64())?;
            }
            (0.0, 0.0)
        };
        poll(&mut cancelled)?;
        let mut commands = reserve(self.links.len())?;
        let mut connector_dissipation_j = reserve(self.links.len())?;
        let mut connector_actuator_work_j = reserve(self.links.len())?;
        let mut average_relative_twist = reserve(self.links.len())?;
        let mut heat = Dd::ZERO;
        let mut work = Dd::ZERO;
        for (e, link) in self.links.iter().enumerate() {
            poll(&mut cancelled)?;
            let average: [f64; 6] = core::array::from_fn(|c| self.scratch.average[6 * e + c]);
            let c = link.geometry.damping.native_coefficients();
            let impulse: [Dd; 6] = if self.duration_s == 0.0 {
                [Dd::ZERO; 6]
            } else {
                core::array::from_fn(|i| {
                    let mut viscous = Dd::ZERO;
                    for j in 0..6 {
                        viscous = viscous + Dd::from_f64(c[i][j]) * Dd::from_f64(average[j]);
                    }
                    (Dd::from_f64(efforts[e][i]) - viscous) * Dd::from_f64(self.duration_s)
                })
            };
            let loss = match &link.loss {
                Some(loss) => loss.evaluate_into_controlled(
                    &self.scratch.initial,
                    &self.scratch.input,
                    &mut self.scratch.heat_coordinates,
                    &mut cancelled,
                )?,
                None => 0.0,
            };
            let mut supplied = Dd::ZERO;
            if self.duration_s != 0.0 {
                for i in 0..6 {
                    supplied = supplied
                        + Dd::from_f64(efforts[e][i])
                            * Dd::from_f64(average[i])
                            * Dd::from_f64(self.duration_s);
                }
            }
            heat = heat + Dd::from_f64(loss);
            work = work + supplied;
            connector_dissipation_j.push(loss);
            connector_actuator_work_j.push(finite(supplied.to_f64())?);
            average_relative_twist.push(average);
            commands.push(InternalWrenchLink {
                bodies: link.bodies,
                impulse: InternalWrenchImpulse {
                    arms_body_m: link.geometry.arms_body_m,
                    impulse_on_second_world_ns: rounded(rotate_wide(
                        link.geometry.axes_world,
                        [impulse[0], impulse[1], impulse[2]],
                    ))?,
                    couple_on_second_world_nms: rounded(rotate_wide(
                        link.geometry.axes_world,
                        [impulse[3], impulse[4], impulse[5]],
                    ))?,
                },
            });
        }
        let native = apply_internal_wrench_network(
            states,
            &self.properties,
            &commands,
            InternalWrenchNetworkControls {
                balance: controls.impulse,
                maximum_bodies: states.len(),
                maximum_connectors: self.links.len(),
                maximum_live_bytes: self.planned_live_bytes,
            },
            &mut cancelled,
        )?;
        let mut states_after = reserve(states.len())?;
        states_after.extend(native.bodies.iter().map(|body| body.state_after));
        for (e, link) in self.links.iter().enumerate() {
            poll(&mut cancelled)?;
            let after = link
                .geometry
                .relative_twist(link.bodies.map(|i| states_after[i]))
                .map_err(|error| connector_error(e, error))?;
            for c in 0..6 {
                let expected = port_rate(&self.input_map, m, 6 * e + c, &self.scratch.endpoint);
                let error = finite((after[c] - expected).to_f64())?.abs();
                if c < 3 {
                    linear_error = linear_error.max(error);
                } else {
                    angular_error = angular_error.max(error);
                }
            }
        }
        let kinetic_balance =
            finite((Dd::from_f64(native.kinetic_energy_change_j) + heat - work).to_f64())?;
        let connector_balance =
            finite((Dd::from_f64(native.connector_work_j) + heat - work).to_f64())?;
        let heat_partition = finite((heat - Dd::from_f64(flow_heat)).to_f64())?;
        let actuator_partition = finite((work - Dd::from_f64(flow_work)).to_f64())?;
        admit(
            controls,
            linear_error,
            angular_error,
            [
                kinetic_balance,
                connector_balance,
                heat_partition,
                actuator_partition,
            ],
        )?;
        let result = DrivenViscousNetworkReceipt {
            states_after,
            connector_actuator_work_j,
            connector_dissipation_j,
            average_relative_twist,
            source_duration_s: self.duration_s,
            actuator_work_j: finite(work.to_f64())?,
            dissipated_energy_j: finite(heat.to_f64())?,
            kinetic_energy_change_j: native.kinetic_energy_change_j,
            connector_work_j: native.connector_work_j,
            angular_origin_world_m: controls.impulse.angular_origin_world_m,
            linear_momentum_change_world_ns: native.linear_momentum_change_world_ns,
            angular_momentum_change_world_nms: native.angular_momentum_change_world_nms,
            kinetic_balance_residual_j: kinetic_balance,
            connector_balance_residual_j: connector_balance,
            heat_partition_residual_j: heat_partition,
            actuator_partition_residual_j: actuator_partition,
            maximum_linear_speed_residual_m_s: linear_error,
            maximum_angular_speed_residual_rad_s: angular_error,
            planned_live_bytes: self.planned_live_bytes,
            planned_preparation_terms: self.planned_preparation_terms,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

fn port_rate(map: &[f64], ports: usize, port: usize, state: &[f64]) -> Dd {
    let mut result = Dd::ZERO;
    for (i, &x) in state.iter().enumerate() {
        result = result + Dd::from_f64(map[i * ports + port]) * Dd::from_f64(x);
    }
    result
}
fn connector_error(connector: usize, error: ViscousConnectorError) -> DrivenViscousNetworkError {
    match error {
        ViscousConnectorError::Impulse(InternalWrenchImpulseError::Cancelled) => {
            DrivenViscousNetworkError::Cancelled
        }
        error => DrivenViscousNetworkError::Connector { connector, error },
    }
}
fn admit(
    controls: ViscousConnectorControls,
    linear: f64,
    angular: f64,
    balance: [f64; 4],
) -> Result<(), DrivenViscousNetworkError> {
    if !linear.is_finite()
        || !angular.is_finite()
        || linear > controls.linear_speed_tolerance_m_s
        || angular > controls.angular_speed_tolerance_rad_s
        || balance
            .iter()
            .any(|x| !x.is_finite() || x.abs() > controls.impulse.work_tolerance_j)
    {
        return Err(DrivenViscousNetworkError::Unresolved {
            linear_speed_residual_m_s: linear,
            angular_speed_residual_rad_s: angular,
            kinetic_balance_residual_j: balance[0],
            connector_balance_residual_j: balance[1],
            heat_partition_residual_j: balance[2],
            actuator_partition_residual_j: balance[3],
        });
    }
    Ok(())
}

// Shared with the moving owner: its two preparations and retained receipts are
// admitted before either endpoint source allocates.
pub(in crate::spatial_loads) fn planned_costs(
    bodies: usize,
    connectors: usize,
    duration_s: f64,
) -> Result<(u128, u128), DrivenViscousNetworkError> {
    let n = 6 * bodies as u128;
    let m = 6 * connectors as u128;
    let d = n + m + 1;
    if d > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128
        || 2 * n > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128
    {
        return Err(DrivenViscousNetworkError::Invalid(
            "complete driven network exceeds the original dense coordinate ceiling",
        ));
    }
    // Conservative simultaneous envelope for the exact flow, complete source
    // rows, each retained source-square cache and the largest cold Gramian.
    // Native body/connector result banks and all source scratch are included.
    let dense = 8 * ((64 + 4 * m) * d * d + connectors as u128 * (8 * n * n + n * m));
    let native = bodies as u128
        * (core::mem::size_of::<Pose>()
            + core::mem::size_of::<MassProperties>()
            + 4 * core::mem::size_of::<RigidBodyState>()
            + 512) as u128
        + connectors as u128
            * (core::mem::size_of::<PreparedLink>()
                + core::mem::size_of::<InternalWrenchLink>()
                + 512) as u128;
    let terms = if duration_s == 0.0 || connectors == 0 {
        0
    } else {
        d.pow(3) * (224 * m.max(1) + 2048)
            + connectors as u128 * ((2 * n).pow(3) * 2048 + n * m)
            + 64 * n.pow(3)
            + 64 * connectors as u128 * n * n
    };
    Ok((dense + native, terms))
}

// Shared native graph geometry and original positive Kelvin rows. This stage
// performs no source evolution and does not assign a material memory model.
struct GraphGeometry {
    poses: Vec<Pose>,
    properties: Vec<MassProperties>,
    links: Vec<PreparedLink>,
    input_map: Vec<f64>,
    rows: Vec<f64>,
    resistance: Vec<Dd>,
}
fn prepare_geometry(
    poses: &[Pose],
    properties: &[MassProperties],
    links: &[ViscousNetworkLink],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GraphGeometry, DrivenViscousNetworkError> {
    let (n, m) = (6 * poses.len(), 6 * links.len());
    let mut bound_poses = reserve(poses.len())?;
    let mut bound_properties = reserve(properties.len())?;
    for (&pose, &properties) in poses.iter().zip(properties) {
        poll(&mut *cancelled)?;
        let zero = RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO)?;
        zero.kinetic_energy_change_to(properties, zero)?;
        bound_poses.push(pose);
        bound_properties.push(properties);
    }
    let mut prepared = reserve(links.len())?;
    let mut input_map = vec![0.0; n * m];
    let mut rows = vec![0.0; m * n];
    let mut resistance = vec![Dd::ZERO; n * n];
    for (e, link) in links.iter().enumerate() {
        poll(&mut *cancelled)?;
        let geometry = PreparedViscousConnector::new(
            link.bodies.map(|i| poses[i]),
            link.bodies.map(|i| properties[i]),
            link.arms_body_m,
            link.damping.clone(),
            0.0,
            &mut *cancelled,
        )
        .map_err(|error| connector_error(e, error))?;
        let columns = geometry
            .kinetic_columns()
            .map_err(|error| connector_error(e, error))?;
        for c in 0..6 {
            for side in 0..2 {
                for axis in 0..6 {
                    input_map[(6 * link.bodies[side] + axis) * m + 6 * e + c] =
                        finite(columns[c][6 * side + axis].to_f64())?;
                }
            }
        }
        let factor = link
            .damping
            .gram_factor()
            .map_err(|error| connector_error(e, error.into()))?;
        let c = link.damping.native_coefficients();
        // Original C supplies the generator; its native checked source
        // factor supplies independent positive heat. Neither substitutes
        // a diagonal or clipped law for the other.
        let mut cb = vec![0.0; 6 * n];
        for k in 0..6 {
            poll(&mut *cancelled)?;
            for i in 0..n {
                let mut value = Dd::ZERO;
                let mut source = Dd::ZERO;
                for j in 0..6 {
                    value =
                        value + Dd::from_f64(c[k][j]) * Dd::from_f64(input_map[i * m + 6 * e + j]);
                    source = source
                        + Dd::from_f64(factor[j][k]) * Dd::from_f64(input_map[i * m + 6 * e + j]);
                }
                cb[k * n + i] = finite(value.to_f64())?;
                rows[(6 * e + k) * n + i] = finite(source.to_f64())?;
            }
        }
        for i in 0..n {
            poll(&mut *cancelled)?;
            for j in i..n {
                let mut entry = Dd::ZERO;
                for k in 0..6 {
                    entry = entry
                        + Dd::from_f64(input_map[i * m + 6 * e + k]) * Dd::from_f64(cb[k * n + j]);
                }
                resistance[i * n + j] = resistance[i * n + j] + entry;
                resistance[j * n + i] = resistance[i * n + j];
            }
        }
        prepared.push(PreparedLink {
            bodies: link.bodies,
            geometry,
            loss: None,
        });
    }

    Ok(GraphGeometry {
        poses: bound_poses,
        properties: bound_properties,
        links: prepared,
        input_map,
        rows,
        resistance,
    })
}
