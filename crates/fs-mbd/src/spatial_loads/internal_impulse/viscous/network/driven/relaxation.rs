//! Generalized Maxwell histories on the complete native spatial graph.
//!
//! Each original branch supplies six material-frame Gram rows L and a positive
//! time tau. Its COMPLETE elastic memory a obeys a' = L s - a/tau,
//! wrench = L^T a, storage = |a|^2/2 and heat rate = |a|^2/tau.
//! The original body kinetic coordinates y obey
//! y' = B(f - C B^T y - sum L^T a), s = B^T y.
//! Identity storage, skew body-memory coupling and positive original resistance
//! enter the existing exact quadratic pHS owner. Every factor row and history,
//! including uncoupled rows and equal-time branches, survives.
//!
//! This is an explicitly authored corotated rate law on the native wrench
//! port. It supplies no equilibrium spring, rest geometry, fitted coefficient,
//! strain validity bound, contact constraint or measured-material claim.

use super::super::super::PassiveBushingMatrix;
use super::{
    DrivenViscousNetworkError, DrivenViscousNetworkLimits, GraphGeometry, PreparedLink,
    ViscousConnectorControls, ViscousNetworkLink, checked_vector_divide, connector_error, finite,
    poll, port_rate, prepare_geometry, rotate_wide, rounded, validate_controls,
};
use crate::spatial_loads::internal_impulse::{
    InternalWrenchImpulse, InternalWrenchImpulseError,
    network::{
        InternalWrenchLink, InternalWrenchNetworkControls, InternalWrenchNetworkError, admit_limit,
        apply_internal_wrench_network, reserve,
    },
};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, Vec3};
use fs_math::{
    dd::Dd,
    det,
    eft::{two_prod, two_sum},
};
use fs_phs::{
    PhsError, PortHamiltonian, PreparedStepError, QuadraticStorage, ResistanceWitness,
    linear_exact::{
        DenseExponentialBudget, LinearSplitBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        PreparedSourceLoss, prepared::PreparedQuadraticFlow, split_quadratic,
    },
};

/// One original PSD spatial stiffness factor and relaxation time.
/// L maps [translation m; rotation rad] to elastic energy coordinates [sqrt(J)].
/// Thus its actual stiffness is L^T L, with native work-conjugate SI units.
/// The factor, not a refactorization or eigenspace, defines the history basis.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialRelaxationBranch {
    gram_rows: [[f64; 6]; 6],
    relaxation_time_s: f64,
}
impl SpatialRelaxationBranch {
    pub fn from_gram_rows(
        gram_rows: [[f64; 6]; 6],
        relaxation_time_s: f64,
    ) -> Result<Self, RelaxationNetworkError> {
        let rate = relaxation_time_s.recip();
        if gram_rows.iter().flatten().any(|x| !x.is_finite())
            || !relaxation_time_s.is_finite()
            || relaxation_time_s <= 0.0
            || !rate.is_finite()
            || rate <= 0.0
        {
            return Err(RelaxationNetworkError::Invalid(
                "finite original Gram rows and positive representable relaxation rate required",
            ));
        }
        // Admit the actual authored stiffness without changing its factor.
        for i in 0..6 {
            for j in 0..6 {
                let mut value = Dd::ZERO;
                for row in &gram_rows {
                    value = value + Dd::from_f64(row[i]) * Dd::from_f64(row[j]);
                }
                if !value.to_f64().is_finite() {
                    return Err(RelaxationNetworkError::Invalid(
                        "original branch stiffness overflow",
                    ));
                }
            }
        }
        Ok(Self {
            gram_rows,
            relaxation_time_s,
        })
    }

    /// Use the existing strictly admitted spatial stiffness factor. Its native
    /// reconstruction rounding is retained; it never supplies a material fit.
    pub fn from_stiffness(
        stiffness: &PassiveBushingMatrix,
        relaxation_time_s: f64,
    ) -> Result<Self, RelaxationNetworkError> {
        let columns = stiffness.gram_factor().map_err(|_| {
            RelaxationNetworkError::Invalid("original spatial stiffness factor is unrepresentable")
        })?;
        Self::from_gram_rows(
            core::array::from_fn(|i| core::array::from_fn(|j| columns[j][i])),
            relaxation_time_s,
        )
    }
    pub const fn gram_rows(&self) -> &[[f64; 6]; 6] {
        &self.gram_rows
    }
    pub const fn relaxation_time_s(&self) -> f64 {
        self.relaxation_time_s
    }

    /// Initialize a physically specified elastic displacement in THIS factor
    /// basis. Directly authored memories remain allowed and are never projected.
    pub fn elastic_memory(
        &self,
        displacement: [f64; 6],
    ) -> Result<[f64; 6], RelaxationNetworkError> {
        if displacement.iter().any(|x| !x.is_finite()) {
            return Err(RelaxationNetworkError::Invalid(
                "finite elastic displacement required",
            ));
        }
        let mut memory = [0.0; 6];
        for (row, value) in self.gram_rows.iter().zip(&mut memory) {
            let mut sum = Dd::ZERO;
            for (&a, &x) in row.iter().zip(&displacement) {
                sum = sum + Dd::from_f64(a) * Dd::from_f64(x);
            }
            *value = finite(sum.to_f64())?;
        }
        Ok(memory)
    }
}

/// Original link, Kelvin viscosity and distinct Maxwell branch sources.
/// An empty branch list reproduces the original driven viscous graph.
#[derive(Clone, Debug, PartialEq)]
pub struct RelaxationNetworkLink {
    pub viscous: ViscousNetworkLink,
    pub branches: Vec<SpatialRelaxationBranch>,
}

/// Held viscous kinetics relative to the ORIGINAL reference laws.
///
/// One multiplier per Kelvin form and one per original Maxwell branch. Material
/// Gram rows, branch order and complete elastic histories keep their reference
/// identity; rate changes alter no stored energy and create no parameter work.
/// This low-level bank carries no empirical temperature law or validity claim.
/// Higher material owners evaluate their declared calibration and always scale
/// the reference law, never an already shifted result.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeldRelaxationKinetics<'a> {
    /// Finite nonnegative multipliers on original Kelvin viscosity.
    pub viscosity_scales: &'a [f64],
    /// Positive time ratios tau/tau_ref in original [link][branch] order.
    pub relaxation_time_scales: &'a [Vec<f64>],
}
impl HeldRelaxationKinetics<'_> {
    /// Admit every effective coefficient without allocating or preparing flow.
    /// Work is linear in the supplied graph and branch counts, with bounded
    /// six-coordinate matrix work. Original Kelvin Gram directions are retained
    /// by their existing scaling owner; memory factors are never refactored.
    pub fn validate(&self, links: &[RelaxationNetworkLink]) -> Result<(), RelaxationNetworkError> {
        self.validate_layout(links)?;
        for (e, link) in links.iter().enumerate() {
            self.damping(e, &link.viscous.damping)?;
            for (branch, source) in link.branches.iter().enumerate() {
                self.relaxation_time(e, branch, source.relaxation_time_s)?;
            }
        }
        Ok(())
    }
    fn validate_layout(
        &self,
        links: &[RelaxationNetworkLink],
    ) -> Result<(), RelaxationNetworkError> {
        if self.viscosity_scales.len() != links.len()
            || self.relaxation_time_scales.len() != links.len()
            || self
                .relaxation_time_scales
                .iter()
                .zip(links)
                .any(|(scales, link)| scales.len() != link.branches.len())
        {
            return Err(RelaxationNetworkError::Invalid(
                "held kinetics require complete original connector and branch banks",
            ));
        }
        Ok(())
    }
    fn damping(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
    ) -> Result<PassiveBushingMatrix, RelaxationNetworkError> {
        Ok(reference
            .scaled(self.viscosity_scales[connector])
            .map_err(|error| connector_error(connector, error.into()))?)
    }
    fn relaxation_time(
        &self,
        connector: usize,
        branch: usize,
        reference_s: f64,
    ) -> Result<f64, RelaxationNetworkError> {
        let scale = self.relaxation_time_scales[connector][branch];
        let time = reference_s * scale;
        let rate = time.recip();
        if !scale.is_finite()
            || scale <= 0.0
            || !time.is_finite()
            || time <= 0.0
            || !rate.is_finite()
            || rate <= 0.0
        {
            return Err(RelaxationNetworkError::Invalid(
                "held Maxwell time scales require positive representable times and rates",
            ));
        }
        Ok(time)
    }
}

/// Complete [link][original branch][six original factor rows] elastic histories.
/// A zero row may hold independently relaxing stored energy. Such a supplied
/// history is kept, never projected away or relabeled as mechanical work.
pub type RelaxationNetworkMemory = Vec<Vec<[f64; 6]>>;

/// Complete source bounds, including every branch coordinate and source cache.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RelaxationNetworkLimits {
    pub source: DrivenViscousNetworkLimits,
    pub maximum_branches: usize,
}

#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct RelaxationNetworkReceipt {
    pub states_after: Vec<RigidBodyState>,
    pub memory_after: RelaxationNetworkMemory,
    pub connector_actuator_work_j: Vec<f64>,
    pub connector_dissipation_j: Vec<f64>,
    pub average_relative_twist: Vec<[f64; 6]>,
    pub source_duration_s: f64,
    pub actuator_work_j: f64,
    pub dissipated_energy_j: f64,
    pub memory_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub connector_work_j: f64,
    pub angular_origin_world_m: Vec3,
    pub linear_momentum_change_world_ns: Vec3,
    pub angular_momentum_change_world_nms: Vec3,
    /// Native kinetic change + memory change + positive heat - actuator work.
    pub kinetic_memory_balance_residual_j: f64,
    /// Independent native impulse work + memory change + heat - drive work.
    pub connector_memory_balance_residual_j: f64,
    /// Full kinetic-coordinate source storage change + heat - drive work.
    pub source_balance_residual_j: f64,
    pub maximum_linear_speed_residual_m_s: f64,
    pub maximum_angular_speed_residual_rad_s: f64,
    pub planned_live_bytes: usize,
    pub planned_preparation_terms: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum RelaxationNetworkError {
    Invalid(&'static str),
    Native(DrivenViscousNetworkError),
    Unresolved {
        linear_speed_residual_m_s: f64,
        angular_speed_residual_rad_s: f64,
        kinetic_memory_balance_residual_j: f64,
        connector_memory_balance_residual_j: f64,
        source_balance_residual_j: f64,
    },
    Cancelled,
}
impl From<DrivenViscousNetworkError> for RelaxationNetworkError {
    fn from(error: DrivenViscousNetworkError) -> Self {
        match error {
            DrivenViscousNetworkError::Cancelled => Self::Cancelled,
            error => Self::Native(error),
        }
    }
}
impl From<InternalWrenchImpulseError> for RelaxationNetworkError {
    fn from(error: InternalWrenchImpulseError) -> Self {
        DrivenViscousNetworkError::from(error).into()
    }
}
impl From<InternalWrenchNetworkError> for RelaxationNetworkError {
    fn from(error: InternalWrenchNetworkError) -> Self {
        DrivenViscousNetworkError::from(error).into()
    }
}
impl From<DynamicsError> for RelaxationNetworkError {
    fn from(error: DynamicsError) -> Self {
        DrivenViscousNetworkError::from(error).into()
    }
}
impl From<PhsError> for RelaxationNetworkError {
    fn from(error: PhsError) -> Self {
        DrivenViscousNetworkError::from(error).into()
    }
}
impl From<PreparedStepError> for RelaxationNetworkError {
    fn from(error: PreparedStepError) -> Self {
        DrivenViscousNetworkError::from(error).into()
    }
}
impl core::fmt::Display for RelaxationNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial relaxation network refused: {self:?}")
    }
}
impl std::error::Error for RelaxationNetworkError {}

struct Scratch {
    initial: Vec<f64>,
    endpoint: Vec<f64>,
    input: Vec<f64>,
    augmented: Vec<f64>,
    balanced: Vec<f64>,
    integrals: Vec<f64>,
    heat_coordinates: Vec<f64>,
    null_initial: Vec<f64>,
    null_heat: Vec<Dd>,
}

// An original null-stress initial history is an exact invariant solution,
// even with nonzero graph drive. Its cross energy with the driven response
// vanishes because that response lies in range(L). Cache its scalar law once.
struct PreparedMemory {
    rows: [[f64; 6]; 6],
    relaxation_time_s: f64,
    decay: f64,
    heat_fraction: f64,
    separable: bool,
}
impl PreparedMemory {
    fn exactly_null(&self, memory: &[f64; 6]) -> bool {
        if !self.separable {
            return false;
        }
        // Identical conservative EFT witness to the native source-loss
        // stationarity rule. The floor only proves no product remainder could
        // underflow; failure keeps the complete original dense evolution.
        let floor = f64::MIN_POSITIVE / (f64::EPSILON * f64::EPSILON);
        for c in 0..6 {
            let mut sum = 0.0;
            for (row, &value) in self.rows.iter().zip(memory) {
                if row[c] == 0.0 || value == 0.0 {
                    continue;
                }
                let (product, product_error) = two_prod(row[c], value);
                if !product.is_finite() || product.abs() < floor || product_error != 0.0 {
                    return false;
                }
                let (next, sum_error) = two_sum(sum, product);
                if !next.is_finite() || sum_error != 0.0 {
                    return false;
                }
                sum = next;
            }
            if sum != 0.0 {
                return false;
            }
        }
        true
    }
}

/// Fixed original geometry, inertia, Kelvin forms, Maxwell factors and clock.
/// Cold work reuses the native wrench map, existing exact pHS flow and positive
/// source-square integrals. Warm calls use cached contractions and one native
/// simultaneous impulse application; no exponential, eigensolve or branch fit.
pub struct PreparedRelaxationNetwork {
    poses: Vec<Pose>,
    properties: Vec<MassProperties>,
    links: Vec<PreparedLink>,
    branch_counts: Vec<usize>,
    memories: Vec<PreparedMemory>,
    input_map: Vec<f64>,
    flow: Option<PreparedQuadraticFlow>,
    scratch: Scratch,
    duration_s: f64,
    planned_live_bytes: usize,
    planned_preparation_terms: u64,
}
impl PreparedRelaxationNetwork {
    pub fn new(
        poses: &[Pose],
        properties: &[MassProperties],
        links: &[RelaxationNetworkLink],
        duration_s: f64,
        limits: RelaxationNetworkLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RelaxationNetworkError> {
        Self::prepare(
            poses, properties, links, None, duration_s, limits, cancelled,
        )
    }

    /// Freeze explicitly supplied kinetics for this source interval.
    ///
    /// Only the original Kelvin viscosity and Maxwell times change. Full
    /// original material factors and history coordinates retain their identity.
    /// The existing exact source integrates its OWN held coefficients and heat;
    /// variation within an interval belongs to the composing time/thermal owner.
    /// No empirical law, entropy term or stiffness-change work is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_kinetics(
        poses: &[Pose],
        properties: &[MassProperties],
        links: &[RelaxationNetworkLink],
        kinetics: HeldRelaxationKinetics<'_>,
        duration_s: f64,
        limits: RelaxationNetworkLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RelaxationNetworkError> {
        Self::prepare(
            poses,
            properties,
            links,
            Some(kinetics),
            duration_s,
            limits,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        poses: &[Pose],
        properties: &[MassProperties],
        links: &[RelaxationNetworkLink],
        kinetics: Option<HeldRelaxationKinetics<'_>>,
        duration_s: f64,
        limits: RelaxationNetworkLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RelaxationNetworkError> {
        poll(&mut cancelled)?;
        if poses.len() != properties.len() || !duration_s.is_finite() || duration_s < 0.0 {
            return Err(RelaxationNetworkError::Invalid(
                "matching original poses/inertias and finite nonnegative duration required",
            ));
        }
        let branches = branch_count(links)?;
        admit_limit(
            "bodies",
            poses.len() as u128,
            limits.source.maximum_bodies as u128,
        )?;
        admit_limit(
            "connectors",
            links.len() as u128,
            limits.source.maximum_connectors as u128,
        )?;
        admit_limit(
            "branches",
            branches as u128,
            limits.maximum_branches as u128,
        )?;
        let (bytes, terms) = planned_costs(poses.len(), links.len(), branches, duration_s)?;
        admit_limit(
            "live bytes",
            bytes,
            limits.source.maximum_live_bytes as u128,
        )?;
        admit_limit(
            "preparation terms",
            terms,
            limits.source.maximum_preparation_terms as u128,
        )?;
        for link in links {
            poll(&mut cancelled)?;
            let [a, b] = link.viscous.bodies;
            if a == b || a >= poses.len() || b >= poses.len() {
                return Err(RelaxationNetworkError::Invalid(
                    "every connector requires two distinct existing original bodies",
                ));
            }
        }
        if let Some(kinetics) = kinetics {
            kinetics.validate_layout(links)?;
        }
        let mut viscous = reserve(links.len())?;
        let mut branch_counts = reserve(links.len())?;
        let mut memories = reserve(branches)?;
        for (e, link) in links.iter().enumerate() {
            poll(&mut cancelled)?;
            let mut effective = link.viscous.clone();
            if let Some(kinetics) = kinetics {
                effective.damping = kinetics.damping(e, &link.viscous.damping)?;
            }
            viscous.push(effective);
            branch_counts.push(link.branches.len());
            for (b, branch) in link.branches.iter().enumerate() {
                poll(&mut cancelled)?;
                let relaxation_time_s = match kinetics {
                    Some(kinetics) => kinetics.relaxation_time(e, b, branch.relaxation_time_s)?,
                    None => branch.relaxation_time_s,
                };
                let distance = duration_s * relaxation_time_s.recip();
                let separable = duration_s > 0.0 && distance > 0.0 && distance.is_finite();
                memories.push(PreparedMemory {
                    rows: branch.gram_rows,
                    relaxation_time_s,
                    decay: if separable { det::exp(-distance) } else { 1.0 },
                    heat_fraction: if separable {
                        -det::expm1(-2.0 * distance)
                    } else {
                        0.0
                    },
                    separable,
                });
            }
        }
        let GraphGeometry {
            poses: bound_poses,
            properties: bound_properties,
            links: mut prepared,
            input_map: body_map,
            rows: kelvin_rows,
            resistance: body_resistance,
        } = prepare_geometry(poses, properties, &viscous, &mut cancelled)?;
        let (body_n, n, m) = (
            6 * poses.len(),
            6 * (poses.len() + branches),
            6 * links.len(),
        );
        let mut input_map = vec![0.0; n * m];
        input_map[..body_n * m].copy_from_slice(&body_map);
        let flow = if duration_s == 0.0 || links.is_empty() {
            None
        } else {
            let mut q = vec![0.0; n * n];
            let mut j = vec![0.0; n * n];
            let mut r = vec![0.0; n * n];
            let mut rows = vec![0.0; (m + 6 * branches) * n];
            // The first m linear forms observe actual slip. The next m observe
            // the complete passive wrench directly, so a weak integrated memory
            // force is never recovered by subtracting two large endpoints.
            let mut observations = vec![vec![0.0; n]; 2 * m];
            for i in 0..n {
                q[i * n + i] = 1.0;
            }
            for i in 0..body_n {
                poll(&mut cancelled)?;
                for k in 0..body_n {
                    r[i * n + k] = finite(body_resistance[i * body_n + k].to_f64())?;
                }
                for p in 0..m {
                    observations[p][i] = body_map[i * m + p];
                }
            }
            for p in 0..m {
                rows[p * n..p * n + body_n]
                    .copy_from_slice(&kelvin_rows[p * body_n..(p + 1) * body_n]);
            }
            let mut memory = body_n;
            let mut row_offset = m;
            for (e, link) in links.iter().enumerate() {
                poll(&mut cancelled)?;
                let c = viscous[e].damping.native_coefficients();
                for component in 0..6 {
                    for i in 0..body_n {
                        let mut value = Dd::ZERO;
                        for k in 0..6 {
                            value = value
                                + Dd::from_f64(c[component][k])
                                    * Dd::from_f64(body_map[i * m + 6 * e + k]);
                        }
                        observations[m + 6 * e + component][i] = finite(value.to_f64())?;
                    }
                }
                for branch in &link.branches {
                    poll(&mut cancelled)?;
                    let rate = memories[(memory - body_n) / 6].relaxation_time_s.recip();
                    let root = det::sqrt(rate);
                    for a in 0..6 {
                        let index = memory + a;
                        r[index * n + index] = rate;
                        rows[(row_offset + a) * n + index] = root;
                        for i in 0..body_n {
                            let mut coupling = Dd::ZERO;
                            for c in 0..6 {
                                coupling = coupling
                                    + Dd::from_f64(body_map[i * m + 6 * e + c])
                                        * Dd::from_f64(branch.gram_rows[a][c]);
                            }
                            let value = finite(coupling.to_f64())?;
                            j[i * n + index] = -value;
                            j[index * n + i] = value;
                        }
                        for c in 0..6 {
                            observations[m + 6 * e + c][index] = branch.gram_rows[a][c];
                        }
                    }
                    memory += 6;
                    row_offset += 6;
                }
            }
            let generator: Vec<_> = j.iter().zip(&r).map(|(j, r)| j - r).collect();
            let witness = ResistanceWitness::try_from_source_rows(n, r.clone(), rows.clone())?;
            let system = PortHamiltonian::new(
                n,
                m,
                j,
                r,
                input_map.clone(),
                Box::new(QuadraticStorage::new(q, n)?),
            )?
            .with_resistance_witness(witness)?;
            poll(&mut cancelled)?;
            let split = split_quadratic(
                &system,
                LinearSplitBudget {
                    max_modes: n,
                    max_omega_dt: 1.0,
                },
            )?;
            let budget = DenseExponentialBudget {
                max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
                max_dense_bytes: bytes as usize,
                max_multiply_terms: terms as u64,
            };
            let flow = PreparedQuadraticFlow::new_controlled(
                &system,
                &split,
                duration_s,
                &[],
                &observations,
                budget,
                &mut cancelled,
            )?;
            let mut row_offset = m;
            for (e, link) in prepared.iter_mut().enumerate() {
                poll(&mut cancelled)?;
                let count = branch_counts[e] * 6;
                let mut source_rows = vec![0.0; (6 + count) * n];
                source_rows[..6 * n].copy_from_slice(&rows[6 * e * n..6 * (e + 1) * n]);
                source_rows[6 * n..]
                    .copy_from_slice(&rows[row_offset * n..(row_offset + count) * n]);
                link.loss = Some(PreparedSourceLoss::new_controlled(
                    n,
                    m,
                    &generator,
                    &input_map,
                    &source_rows,
                    duration_s,
                    budget,
                    &mut cancelled,
                )?);
                row_offset += count;
            }
            Some(flow)
        };
        poll(&mut cancelled)?;
        Ok(Self {
            poses: bound_poses,
            properties: bound_properties,
            links: prepared,
            branch_counts,
            memories,
            input_map,
            flow,
            scratch: Scratch {
                initial: vec![0.0; n],
                endpoint: vec![0.0; n],
                input: vec![0.0; m],
                augmented: vec![0.0; n + m + 1],
                balanced: vec![0.0; n],
                integrals: vec![0.0; 2 * m],
                heat_coordinates: vec![0.0; 2 * n],
                null_initial: vec![0.0; 6 * branches],
                null_heat: vec![Dd::ZERO; links.len()],
            },
            duration_s,
            planned_live_bytes: bytes as usize,
            planned_preparation_terms: terms as u64,
        })
    }

    /// The actual held Kelvin form used by dynamics, impulse and heat.
    pub fn effective_damping(&self, connector: usize) -> Option<&PassiveBushingMatrix> {
        self.links.get(connector).map(|link| &link.geometry.damping)
    }

    /// Actual held relaxation time in original connector/branch order.
    /// No allocation or numerical preparation occurs.
    pub fn effective_relaxation_time_s(&self, connector: usize, branch: usize) -> Option<f64> {
        if branch >= *self.branch_counts.get(connector)? {
            return None;
        }
        let offset: usize = self.branch_counts[..connector].iter().sum();
        Some(self.memories[offset + branch].relaxation_time_s)
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

    /// Simultaneous held actuation, Kelvin resistance and complete memory.
    /// Input banks stay borrowed on every success, refusal and cancellation.
    pub fn advance(
        &mut self,
        states: &[RigidBodyState],
        memory: &[Vec<[f64; 6]>],
        efforts: &[[f64; 6]],
        controls: ViscousConnectorControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RelaxationNetworkReceipt, RelaxationNetworkError> {
        poll(&mut cancelled)?;
        validate_controls(controls.impulse)?;
        if states.len() != self.poses.len()
            || states.iter().zip(&self.poses).any(|(s, p)| s.pose() != *p)
            || efforts.len() != self.links.len()
            || efforts.iter().flatten().any(|x| !x.is_finite())
            || memory.len() != self.branch_counts.len()
            || memory
                .iter()
                .zip(&self.branch_counts)
                .any(|(a, b)| a.len() != *b)
            || memory.iter().flatten().flatten().any(|x| !x.is_finite())
            || [
                controls.linear_speed_tolerance_m_s,
                controls.angular_speed_tolerance_rad_s,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(RelaxationNetworkError::Invalid(
                "original poses, complete finite memory/efforts and nonnegative rate tolerances required",
            ));
        }
        let (body_n, m) = (6 * states.len(), 6 * self.links.len());
        for (i, state) in states.iter().enumerate() {
            poll(&mut cancelled)?;
            let p = checked_vector_divide(
                state.linear_momentum_world(),
                det::sqrt(self.properties[i].mass()),
            )
            .map_err(|_| {
                RelaxationNetworkError::Invalid("original body kinetic coordinate overflow")
            })?;
            let h = self.properties[i].angular_mobility_factor(state.angular_momentum_body());
            let values = [p.x, p.y, p.z, h.x, h.y, h.z];
            if values.iter().any(|x| !x.is_finite()) {
                return Err(RelaxationNetworkError::Invalid(
                    "finite original kinetic coordinates required",
                ));
            }
            self.scratch.initial[6 * i..6 * (i + 1)].copy_from_slice(&values);
        }
        self.scratch.null_initial.fill(0.0);
        self.scratch.null_heat.fill(Dd::ZERO);
        let mut offset = body_n;
        let mut branch_index = 0;
        for (e, link) in memory.iter().enumerate() {
            poll(&mut cancelled)?;
            for branch in link {
                let prepared = &self.memories[branch_index];
                if prepared.exactly_null(branch) {
                    self.scratch.initial[offset..offset + 6].fill(0.0);
                    self.scratch.null_initial[offset - body_n..offset - body_n + 6]
                        .copy_from_slice(branch);
                    for &value in branch {
                        self.scratch.null_heat[e] = self.scratch.null_heat[e]
                            + Dd::from_f64(value)
                                * Dd::from_f64(value)
                                * Dd::from_f64(0.5 * prepared.heat_fraction);
                    }
                    finite(self.scratch.null_heat[e].to_f64())?;
                } else {
                    self.scratch.initial[offset..offset + 6].copy_from_slice(branch);
                }
                offset += 6;
                branch_index += 1;
            }
        }
        for (e, effort) in efforts.iter().enumerate() {
            self.scratch.input[6 * e..6 * (e + 1)].copy_from_slice(effort);
        }
        let (mut linear_error, mut angular_error) = (0.0_f64, 0.0_f64);
        self.rate_errors(
            states,
            &self.scratch.initial,
            &mut linear_error,
            &mut angular_error,
            &mut cancelled,
        )?;
        admit(controls, linear_error, angular_error, [0.0; 3])?;
        if let Some(flow) = &self.flow {
            flow.evaluate_into(
                &self.scratch.initial,
                &self.scratch.input,
                &mut self.scratch.augmented,
                &mut self.scratch.balanced,
                &mut self.scratch.endpoint,
                &mut self.scratch.integrals,
            )?;
        } else {
            self.scratch.endpoint.copy_from_slice(&self.scratch.initial);
            self.scratch.integrals.fill(0.0);
        }
        poll(&mut cancelled)?;
        let mut commands = reserve(self.links.len())?;
        let mut connector_dissipation_j = reserve(self.links.len())?;
        let mut connector_actuator_work_j = reserve(self.links.len())?;
        let mut average_relative_twist = reserve(self.links.len())?;
        let (mut heat, mut work) = (Dd::ZERO, Dd::ZERO);
        for (e, link) in self.links.iter().enumerate() {
            poll(&mut cancelled)?;
            let impulse: [Dd; 6] = if self.duration_s == 0.0 {
                [Dd::ZERO; 6]
            } else {
                core::array::from_fn(|c| {
                    Dd::from_f64(efforts[e][c]) * Dd::from_f64(self.duration_s)
                        - Dd::from_f64(self.scratch.integrals[m + 6 * e + c])
                })
            };
            let mut average = [0.0; 6];
            for (c, value) in average.iter_mut().enumerate() {
                *value = if self.duration_s == 0.0 {
                    finite(
                        port_rate(&self.input_map, m, 6 * e + c, &self.scratch.initial).to_f64(),
                    )?
                } else {
                    finite(self.scratch.integrals[6 * e + c] / self.duration_s)?
                };
            }
            let loss = match &link.loss {
                Some(loss) => loss.evaluate_into_controlled(
                    &self.scratch.initial,
                    &self.scratch.input,
                    &mut self.scratch.heat_coordinates,
                    &mut cancelled,
                )?,
                None => 0.0,
            };
            let loss = finite((Dd::from_f64(loss) + self.scratch.null_heat[e]).to_f64())?;
            let mut supplied = Dd::ZERO;
            if self.duration_s != 0.0 {
                for c in 0..6 {
                    supplied = supplied
                        + Dd::from_f64(efforts[e][c])
                            * Dd::from_f64(self.scratch.integrals[6 * e + c]);
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
        self.rate_errors(
            &states_after,
            &self.scratch.endpoint,
            &mut linear_error,
            &mut angular_error,
            &mut cancelled,
        )?;
        let mut memory_after = reserve(memory.len())?;
        let mut offset = body_n;
        for source in memory {
            poll(&mut cancelled)?;
            let mut destination = reserve(source.len())?;
            for _ in source {
                let mut branch = [0.0; 6];
                let prepared = &self.memories[(offset - body_n) / 6];
                for (c, value) in branch.iter_mut().enumerate() {
                    *value = finite(
                        (Dd::from_f64(self.scratch.endpoint[offset + c])
                            + Dd::from_f64(self.scratch.null_initial[offset - body_n + c])
                                * Dd::from_f64(prepared.decay))
                        .to_f64(),
                    )?;
                }
                destination.push(branch);
                offset += 6;
            }
            memory_after.push(destination);
        }
        let memory_change = memory_energy_change_retained(memory, &memory_after)?;
        // Use the actual complete original memory endpoint, including the
        // independent homogeneous solution, for this source-coordinate gate.
        let source_change = coordinate_energy_change(
            &self.scratch.initial[..body_n],
            &self.scratch.endpoint[..body_n],
        )? + memory_change;
        let balances = [
            finite(
                (Dd::from_f64(native.kinetic_energy_change_j) + memory_change + heat - work)
                    .to_f64(),
            )?,
            finite((Dd::from_f64(native.connector_work_j) + memory_change + heat - work).to_f64())?,
            finite((source_change + heat - work).to_f64())?,
        ];
        admit(controls, linear_error, angular_error, balances)?;
        let result = RelaxationNetworkReceipt {
            states_after,
            memory_after,
            connector_actuator_work_j,
            connector_dissipation_j,
            average_relative_twist,
            source_duration_s: self.duration_s,
            actuator_work_j: finite(work.to_f64())?,
            dissipated_energy_j: finite(heat.to_f64())?,
            memory_energy_change_j: finite(memory_change.to_f64())?,
            kinetic_energy_change_j: native.kinetic_energy_change_j,
            connector_work_j: native.connector_work_j,
            angular_origin_world_m: controls.impulse.angular_origin_world_m,
            linear_momentum_change_world_ns: native.linear_momentum_change_world_ns,
            angular_momentum_change_world_nms: native.angular_momentum_change_world_nms,
            kinetic_memory_balance_residual_j: balances[0],
            connector_memory_balance_residual_j: balances[1],
            source_balance_residual_j: balances[2],
            maximum_linear_speed_residual_m_s: linear_error,
            maximum_angular_speed_residual_rad_s: angular_error,
            planned_live_bytes: self.planned_live_bytes,
            planned_preparation_terms: self.planned_preparation_terms,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn rate_errors(
        &self,
        states: &[RigidBodyState],
        coordinates: &[f64],
        linear: &mut f64,
        angular: &mut f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RelaxationNetworkError> {
        let m = 6 * self.links.len();
        for (e, link) in self.links.iter().enumerate() {
            poll(&mut *cancelled)?;
            let native = link
                .geometry
                .relative_twist(link.bodies.map(|i| states[i]))
                .map_err(|error| connector_error(e, error))?;
            for c in 0..6 {
                let expected = port_rate(&self.input_map, m, 6 * e + c, coordinates);
                let error = finite((native[c] - expected).to_f64())?.abs();
                if c < 3 {
                    *linear = linear.max(error);
                } else {
                    *angular = angular.max(error);
                }
            }
        }
        Ok(())
    }
}

/// Original-to-final stored-energy change, retaining differences before
/// multiplication. It never subtracts rounded absolute stored energies.
pub fn relaxation_memory_energy_change(
    before: &[Vec<[f64; 6]>],
    after: &[Vec<[f64; 6]>],
) -> Result<f64, RelaxationNetworkError> {
    Ok(finite(
        memory_energy_change_retained(before, after)?.to_f64(),
    )?)
}
pub(in crate::spatial_loads) fn memory_energy_change_retained(
    before: &[Vec<[f64; 6]>],
    after: &[Vec<[f64; 6]>],
) -> Result<Dd, RelaxationNetworkError> {
    if before.len() != after.len() || before.iter().zip(after).any(|(a, b)| a.len() != b.len()) {
        return Err(RelaxationNetworkError::Invalid(
            "complete matching original memory banks required",
        ));
    }
    let mut energy = Dd::ZERO;
    for (a, b) in before.iter().flatten().zip(after.iter().flatten()) {
        energy = energy + coordinate_energy_change(a, b)?;
    }
    finite(energy.to_f64())?;
    Ok(energy)
}
fn coordinate_energy_change(before: &[f64], after: &[f64]) -> Result<Dd, RelaxationNetworkError> {
    let mut result = Dd::ZERO;
    for (&a, &b) in before.iter().zip(after) {
        if !a.is_finite() || !b.is_finite() {
            return Err(RelaxationNetworkError::Invalid(
                "finite original energy coordinates required",
            ));
        }
        let (a, b) = (Dd::from_f64(a), Dd::from_f64(b));
        result = result + (b - a) * (b + a) * Dd::from_f64(0.5);
    }
    finite(result.to_f64())?;
    Ok(result)
}
fn admit(
    controls: ViscousConnectorControls,
    linear: f64,
    angular: f64,
    balance: [f64; 3],
) -> Result<(), RelaxationNetworkError> {
    if !linear.is_finite()
        || !angular.is_finite()
        || linear > controls.linear_speed_tolerance_m_s
        || angular > controls.angular_speed_tolerance_rad_s
        || balance
            .iter()
            .any(|x| !x.is_finite() || x.abs() > controls.impulse.work_tolerance_j)
    {
        return Err(RelaxationNetworkError::Unresolved {
            linear_speed_residual_m_s: linear,
            angular_speed_residual_rad_s: angular,
            kinetic_memory_balance_residual_j: balance[0],
            connector_memory_balance_residual_j: balance[1],
            source_balance_residual_j: balance[2],
        });
    }
    Ok(())
}
pub(in crate::spatial_loads) fn branch_count(
    links: &[RelaxationNetworkLink],
) -> Result<usize, RelaxationNetworkError> {
    links.iter().try_fold(0usize, |sum, link| {
        sum.checked_add(link.branches.len())
            .ok_or(RelaxationNetworkError::Invalid(
                "complete branch extent overflow",
            ))
    })
}
pub(in crate::spatial_loads) fn planned_costs(
    bodies: usize,
    connectors: usize,
    branches: usize,
    duration_s: f64,
) -> Result<(u128, u128), RelaxationNetworkError> {
    let n = 6 * (bodies as u128 + branches as u128);
    let m = 6 * connectors as u128;
    let d = n + m + 1;
    if d > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128
        || 2 * n > MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as u128
    {
        return Err(RelaxationNetworkError::Invalid(
            "complete hereditary network exceeds original dense coordinate ceiling",
        ));
    }
    // Includes complete native geometry, all original memory/observer rows,
    // held-flow and positive-loss caches, largest cold Gramian and output banks.
    let dense = 8
        * ((80 + 8 * m) * d * d
            + connectors as u128 * (8 * n * n + n * m)
            + (m + 6 * branches as u128) * n);
    let native = bodies as u128
        * (core::mem::size_of::<Pose>()
            + core::mem::size_of::<MassProperties>()
            + 4 * core::mem::size_of::<RigidBodyState>()
            + 512) as u128
        + connectors as u128
            * (core::mem::size_of::<PreparedLink>()
                + core::mem::size_of::<InternalWrenchLink>()
                + 1024) as u128
        + branches as u128 * (1024 + core::mem::size_of::<PreparedMemory>() + 6 * 8) as u128
        + connectors as u128 * core::mem::size_of::<Dd>() as u128;
    let terms = if duration_s == 0.0 || connectors == 0 {
        0
    } else {
        d.pow(3) * (448 * m.max(1) + 2048)
            + connectors as u128 * ((2 * n).pow(3) * 2048 + n * m)
            + 64 * n.pow(3)
            + 128 * (connectors as u128 + branches as u128) * n * n
            + 4096 * branches as u128
    };
    Ok((dense + native, terms))
}
