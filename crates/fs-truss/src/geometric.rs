//! Objective, geometrically nonlinear elasticity of the original truss network.
//!
//! Current chord lengths and directions determine equilibrium. The authored
//! axial law remains linear in engineering strain: N = EA ((l-L)/L - e0),
//! U = EA L ((l-L)/L-e0)^2/2. Thus large rotations are exact, but large
//! material strain is not inferred from a Young modulus. The caller supplies
//! an elastic-strain ceiling. Stress is nominal N/A in the original area;
//! no current area, Poisson contraction, plasticity or member bending is guessed.
//!
//! Native Newton--Krylov owns fixed-load equilibrium; native Keller bordered
//! continuation owns folds, step control and resumable trajectories. Both use
//! the SAME analytical member tangent, including its actual N/l transverse
//! term. No private nonlinear solver or linearized buckling approximation is
//! substituted for the deformed force equation.

use crate::elastic::ElasticImposedState;
use crate::{GroundStructure, LayoutCase};
use fs_exec::Cx;
use fs_solid::SolidError;
use fs_solid::continuation::{self, ArcSettings, PathResidual, PathState};
use fs_solver::nonlinear::{
    NewtonKrylovConfig, NewtonKrylovState, NonlinearProblem, SolverRunProgress,
};
use fs_sparse::{Coo, Csr};
use std::cell::RefCell;

/// Fixed loads/misfit and prescribed support motion in original node/member order.
#[derive(Debug, Clone, Copy, Default)]
pub struct GeometricLoading<'a> {
    /// Fixed stress-free engineering strains and exact support translations.
    /// They stay fixed while the continuation factor changes.
    pub fixed: ElasticImposedState<'a>,
    /// Fixed physical dead forces [N], such as previously computed gravity.
    pub fixed_loads_n: Option<&'a [[f64; 2]]>,
    /// Additional exact support translation per load factor [m]. Free entries
    /// must be zero. Actual support motion is fixed + lambda * this field.
    pub support_motion_per_factor_m: Option<&'a [[f64; 2]]>,
}

/// Physical units of the dimensionless continuation coordinates and equations.
#[derive(Debug, Clone, Copy)]
pub struct GeometricScales {
    /// u_free = displacement_m * q; arclength uses q and lambda.
    pub displacement_m: f64,
    /// R = physical free force residual / force_n.
    pub force_n: f64,
}

/// Explicit material-domain and numerical-resource admission.
#[derive(Debug, Clone, Copy)]
pub struct GeometricLimits {
    /// Authored bound on |(l-L)/L - e0|, strictly between zero and one.
    /// This is a model restriction, not a measured material linearity limit.
    pub maximum_absolute_elastic_strain: f64,
    /// At most 100000 free translations; no unsupported coordinate is removed.
    pub maximum_free_dofs: usize,
    /// At most 8000000 staged tangent entries, bounding native COO sorting.
    pub maximum_staged_entries: usize,
    /// Conservative combined preparation/assembly/native-solver payload bound.
    /// Includes 64-column Krylov storage or the native small dense border,
    /// two cached assemblies and bounded path/iteration history.
    pub maximum_payload_bytes: usize,
}

/// Physical evaluation. A nonzero free residual means this is not equilibrium.
#[derive(Debug, Clone)]
pub struct GeometricResponse {
    pub load_factor: f64,
    pub displacements_m: Vec<[f64; 2]>,
    /// Original member order; removed members have zero outputs.
    pub member_engineering_strain: Vec<f64>,
    pub member_elastic_strain: Vec<f64>,
    /// Positive tension [N], along the CURRENT chord.
    pub member_force_n: Vec<f64>,
    /// N / original area [Pa], not an inferred Cauchy stress.
    pub member_nominal_stress_pa: Vec<f64>,
    /// Internal minus actual dead force on supported coordinates; zero on free.
    pub support_reactions_n: Vec<[f64; 2]>,
    /// Physical residual [N] in the owner's increasing free_dofs order.
    pub free_force_residual_n: Vec<f64>,
    pub strain_energy_j: f64,
    /// U - (fixed_force + lambda * reference_force) dot u at this factor.
    /// This is a fixed-parameter potential, not accumulated support/load work.
    pub potential_energy_j: f64,
}

struct Member {
    ends: [usize; 2],
    length: f64,
    area: f64,
    rigidity: f64,
    stiffness: f64,
    free_strain: f64,
    // Canonical row offsets of the 4x4 free-coordinate tangent contributions.
    slots: [usize; 16],
}

/// Immutable source geometry, material coefficients and loading family.
pub struct GeometricTruss<'a> {
    ground: &'a GroundStructure,
    case: &'a LayoutCase,
    members: Vec<Member>,
    fixed_load: Vec<[f64; 2]>,
    fixed_motion: Vec<[f64; 2]>,
    ramp_motion: Vec<[f64; 2]>,
    free: Vec<usize>,
    reduced: Vec<usize>,
    pattern: Csr,
    scales: GeometricScales,
    limits: GeometricLimits,
    payload_bytes: usize,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}
fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |sum, &value| sum.hypot(value))
}
fn finite(value: f64) -> Result<f64, SolidError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable geometric truss arithmetic"))
    }
}
fn scaled(value: f64, factor: f64) -> Result<f64, SolidError> {
    let result = finite(value * factor)?;
    if value != 0.0 && factor != 0.0 && result == 0.0 {
        Err(invalid("geometric truss scaling underflowed"))
    } else {
        Ok(result)
    }
}

struct Analysis {
    response: GeometricResponse,
    residual: Vec<f64>,
    tangent: Csr,
    load_column: Vec<f64>,
}

impl<'a> GeometricTruss<'a> {
    /// Admit a complete selected truss. The reference case supplies support
    /// flags and dead forces multiplied by lambda. Zero area removes a member;
    /// every original modulus and supplied input must nevertheless be valid.
    /// No support, mechanism or near-zero stiffness is silently regularized.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn try_new(
        ground: &'a GroundStructure,
        reference_case: &'a LayoutCase,
        areas_m2: &[f64],
        young_pa: &[f64],
        loading: GeometricLoading<'_>,
        scales: GeometricScales,
        limits: GeometricLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let nodes = ground.nodes().len();
        let count = ground.members().len();
        let ndof = nodes
            .checked_mul(2)
            .ok_or_else(|| invalid("truss size overflow"))?;
        if ndof > 800000
            || reference_case.node_count() != nodes
            || areas_m2.len() != count
            || young_pa.len() != count
            || loading
                .fixed
                .member_stress_free_strain
                .is_some_and(|s| s.len() != count)
            || [
                scales.displacement_m,
                scales.force_n,
                1.0 / scales.force_n,
                scales.displacement_m / scales.force_n,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !limits.maximum_absolute_elastic_strain.is_finite()
            || !(0.0..1.0).contains(&limits.maximum_absolute_elastic_strain)
            || limits.maximum_absolute_elastic_strain == 0.0
            || limits.maximum_free_dofs > 100000
            || limits.maximum_staged_entries > 8000000
        {
            return Err(invalid(
                "geometric truss shape, scales, strain or resource limits invalid",
            ));
        }
        let free_count = reference_case
            .supported()
            .iter()
            .flatten()
            .filter(|&&s| !s)
            .count();
        let staged = areas_m2
            .iter()
            .try_fold(0usize, |sum, &area| {
                sum.checked_add(if area > 0.0 { 16 } else { 0 })
            })
            .ok_or_else(|| invalid("geometric truss tangent size overflow"))?;
        let dense = if free_count < 1024 {
            32 * (free_count + 1).pow(2)
        } else {
            0
        };
        let payload = ndof
            .checked_mul(512)
            .and_then(|v| v.checked_add(count.checked_mul(256)?))
            .and_then(|v| v.checked_add(staged.checked_mul(128)?))
            .and_then(|v| v.checked_add(free_count.checked_mul(2048)?))
            .and_then(|v| v.checked_add(dense))
            // 64-column Hessenberg/rotations and bounded cycle histories.
            .and_then(|v| v.checked_add(65536))
            .ok_or_else(|| invalid("geometric truss payload overflow"))?;
        if free_count == 0
            || free_count > limits.maximum_free_dofs
            || staged > limits.maximum_staged_entries
            || payload > limits.maximum_payload_bytes
        {
            return Err(invalid(
                "geometric truss requires admitted nonempty free coordinates, tangent and payload",
            ));
        }
        for (field, supports_only) in [
            (loading.fixed_loads_n, false),
            (loading.fixed.prescribed_displacements_m, true),
            (loading.support_motion_per_factor_m, true),
        ] {
            if let Some(values) = field {
                if values.len() != nodes {
                    return Err(invalid(
                        "geometric truss nodal field has wrong original shape",
                    ));
                }
                for (node, value) in values.iter().enumerate() {
                    if node % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    for (axis, &v) in value.iter().enumerate() {
                        if !v.is_finite()
                            || (supports_only
                                && !reference_case.supported()[node][axis]
                                && v != 0.0)
                        {
                            return Err(invalid(
                                "geometric truss nodal fields must be finite; motion belongs only on supports",
                            ));
                        }
                    }
                }
            }
        }
        let mut members = Vec::with_capacity(count);
        for (index, &(a, b)) in ground.members().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            let area = areas_m2[index];
            let young = young_pa[index];
            let free_strain = loading
                .fixed
                .member_stress_free_strain
                .map_or(0.0, |s| s[index]);
            let length = ground.lengths()[index];
            let rigidity = young * area;
            let stiffness = rigidity / length;
            if !area.is_finite()
                || area < 0.0
                || !young.is_finite()
                || young <= 0.0
                || !free_strain.is_finite()
                || 1.0 + free_strain <= 0.0
                || !rigidity.is_finite()
                || !stiffness.is_finite()
                || (area > 0.0 && (rigidity <= 0.0 || stiffness <= 0.0))
                || !(length * (1.0 + free_strain)).is_finite()
            {
                return Err(invalid(
                    "geometric truss requires representable nonnegative areas, positive E/EA/L and natural lengths",
                ));
            }
            members.push(Member {
                ends: [a, b],
                length,
                area,
                rigidity,
                stiffness,
                free_strain,
                slots: [usize::MAX; 16],
            });
        }
        let mut free = Vec::with_capacity(free_count);
        let mut reduced = vec![usize::MAX; ndof];
        for (dof, entry) in reduced.iter_mut().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if !reference_case.supported()[dof / 2][dof % 2] {
                *entry = free.len();
                free.push(dof);
            }
        }
        // The graph never changes within this owner. Canonicalize its sparse
        // pattern once; warm assemblies only accumulate known row slots.
        let mut coo = Coo::new(free_count, free_count);
        for (index, member) in members.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if member.area == 0.0 {
                continue;
            }
            for i in 0..4 {
                let row = reduced[2 * member.ends[i / 2] + i % 2];
                if row == usize::MAX {
                    continue;
                }
                for j in 0..4 {
                    let column = reduced[2 * member.ends[j / 2] + j % 2];
                    if column != usize::MAX {
                        coo.push(row, column, 0.0);
                    }
                }
            }
        }
        checkpoint(cx)?;
        let pattern = coo.assemble();
        checkpoint(cx)?;
        for (index, member) in members.iter_mut().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if member.area == 0.0 {
                continue;
            }
            for i in 0..4 {
                let row = reduced[2 * member.ends[i / 2] + i % 2];
                if row == usize::MAX {
                    continue;
                }
                for j in 0..4 {
                    let column = reduced[2 * member.ends[j / 2] + j % 2];
                    if column != usize::MAX {
                        member.slots[4 * i + j] =
                            pattern.row(row).0.binary_search(&column).map_err(|_| {
                                invalid("geometric truss pattern lost a member contribution")
                            })?;
                    }
                }
            }
        }
        let prepared = Self {
            ground,
            case: reference_case,
            members,
            fixed_load: loading
                .fixed_loads_n
                .map_or_else(|| vec![[0.0; 2]; nodes], <[_]>::to_vec),
            fixed_motion: loading
                .fixed
                .prescribed_displacements_m
                .map_or_else(|| vec![[0.0; 2]; nodes], <[_]>::to_vec),
            ramp_motion: loading
                .support_motion_per_factor_m
                .map_or_else(|| vec![[0.0; 2]; nodes], <[_]>::to_vec),
            free,
            reduced,
            pattern,
            scales,
            limits,
            payload_bytes: payload,
        };
        checkpoint(cx)?;
        Ok(prepared)
    }

    /// Increasing reduced-to-original node-major xy map.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free
    }
    pub const fn scales(&self) -> GeometricScales {
        self.scales
    }

    fn displacement(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<[f64; 2]>, SolidError> {
        checkpoint(cx)?;
        if q.len() != self.free.len() || !lambda.is_finite() {
            return Err(invalid(
                "dimension-matched finite geometric truss state required",
            ));
        }
        let mut u = self.fixed_motion.clone();
        for (node, value) in u.iter_mut().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for (axis, component) in value.iter_mut().enumerate() {
                *component = finite(*component + finite(lambda * self.ramp_motion[node][axis])?)?;
            }
        }
        for (&dof, &value) in self.free.iter().zip(q) {
            u[dof / 2][dof % 2] = finite(value * self.scales.displacement_m)?;
            if value != 0.0 && u[dof / 2][dof % 2] == 0.0 {
                return Err(invalid(
                    "unrepresentable geometric truss displacement scale",
                ));
            }
        }
        Ok(u)
    }

    /// Pack physical displacement, enforcing the family's exact current supports.
    pub fn free_coordinates(
        &self,
        u: &[[f64; 2]],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        checkpoint(cx)?;
        if u.len() != self.ground.nodes().len() || !lambda.is_finite() {
            return Err(invalid("complete original-order displacement required"));
        }
        let mut q = Vec::with_capacity(self.free.len());
        for (node, value) in u.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for (axis, &v) in value.iter().enumerate() {
                finite(v)?;
                if self.case.supported()[node][axis] {
                    let target = finite(
                        self.fixed_motion[node][axis]
                            + finite(lambda * self.ramp_motion[node][axis])?,
                    )?;
                    if v != target {
                        return Err(invalid(
                            "geometric truss displacement violates prescribed support motion",
                        ));
                    }
                } else {
                    let value = finite(v / self.scales.displacement_m)?;
                    if v != 0.0 && value == 0.0 {
                        return Err(invalid("unrepresentable geometric truss coordinate scale"));
                    }
                    q.push(value);
                }
            }
        }
        Ok(q)
    }

    #[allow(clippy::too_many_lines)] // Single member pass owns force, energy and its exact linearization.
    fn analyze(&self, q: &[f64], lambda: f64, cx: &Cx<'_>) -> Result<Analysis, SolidError> {
        let u = self.displacement(q, lambda, cx)?;
        let n = self.free.len();
        let mut residual = vec![[0.0; 2]; u.len()];
        let mut load_column = vec![0.0; n];
        let mut tangent = self.pattern.clone();
        let mut response = GeometricResponse {
            load_factor: lambda,
            displacements_m: u,
            member_engineering_strain: vec![0.0; self.members.len()],
            member_elastic_strain: vec![0.0; self.members.len()],
            member_force_n: vec![0.0; self.members.len()],
            member_nominal_stress_pa: vec![0.0; self.members.len()],
            support_reactions_n: vec![[0.0; 2]; self.ground.nodes().len()],
            free_force_residual_n: Vec::with_capacity(n),
            strain_energy_j: 0.0,
            potential_energy_j: 0.0,
        };
        let mut dead_work = 0.0;
        for (node, value) in residual.iter_mut().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for (axis, entry) in value.iter_mut().enumerate() {
                let load = finite(
                    self.fixed_load[node][axis] + finite(lambda * self.case.loads()[node][axis])?,
                )?;
                *entry = -load;
                dead_work = finite(dead_work + load * response.displacements_m[node][axis])?;
                let row = self.reduced[2 * node + axis];
                if row != usize::MAX {
                    load_column[row] = self.case.loads()[node][axis];
                }
            }
        }
        for (index, member) in self.members.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if member.area == 0.0 {
                continue;
            }
            let [a, b] = member.ends;
            let u = &response.displacements_m;
            // Relative differences retain invariance under a common translation.
            let chord = [
                finite(
                    (self.ground.nodes()[b][0] - self.ground.nodes()[a][0]) + (u[b][0] - u[a][0]),
                )?,
                finite(
                    (self.ground.nodes()[b][1] - self.ground.nodes()[a][1]) + (u[b][1] - u[a][1]),
                )?,
            ];
            let length = chord[0].hypot(chord[1]);
            if !length.is_finite() || length <= 0.0 {
                return Err(invalid(
                    "retained geometric truss member collapsed or overflowed",
                ));
            }
            let direction = [chord[0] / length, chord[1] / length];
            let strain = finite((length - member.length) / member.length)?;
            let elastic = finite(strain - member.free_strain)?;
            if elastic.abs() > self.limits.maximum_absolute_elastic_strain {
                return Err(invalid(
                    "geometric truss elastic strain exceeds its authored material domain",
                ));
            }
            let force = finite(member.rigidity * elastic)?;
            let transverse = finite(force / length)?;
            response.member_engineering_strain[index] = strain;
            response.member_elastic_strain[index] = elastic;
            response.member_force_n[index] = force;
            response.member_nominal_stress_pa[index] = finite(force / member.area)?;
            response.strain_energy_j =
                finite(response.strain_energy_j + (0.5 * force * elastic) * member.length)?;
            for axis in 0..2 {
                let component = finite(force * direction[axis])?;
                residual[a][axis] = finite(residual[a][axis] - component)?;
                residual[b][axis] = finite(residual[b][axis] + component)?;
            }
            let mut local = [[0.0; 2]; 2];
            for i in 0..2 {
                for j in 0..2 {
                    let product = direction[i] * direction[j];
                    local[i][j] = finite(
                        member.stiffness * product
                            + transverse * (if i == j { 1.0 } else { 0.0 } - product),
                    )?;
                }
            }
            for (local_i, (end_i, sign_i)) in [(a, -1.0), (b, 1.0)].into_iter().enumerate() {
                for axis_i in 0..2 {
                    let row = self.reduced[2 * end_i + axis_i];
                    if row == usize::MAX {
                        continue;
                    }
                    for (local_j, (end_j, sign_j)) in [(a, -1.0), (b, 1.0)].into_iter().enumerate()
                    {
                        for axis_j in 0..2 {
                            let value = sign_i * sign_j * local[axis_i][axis_j];
                            let column = self.reduced[2 * end_j + axis_j];
                            if column != usize::MAX {
                                let offset =
                                    member.slots[4 * (2 * local_i + axis_i) + 2 * local_j + axis_j];
                                let entry = &mut tangent.row_values_mut(row)[offset];
                                *entry = finite(
                                    *entry
                                        + scaled(
                                            value,
                                            self.scales.displacement_m / self.scales.force_n,
                                        )?,
                                )?;
                            } else {
                                // -dR/dlambda = F_reference - K_free,support g.
                                load_column[row] = finite(
                                    load_column[row] - value * self.ramp_motion[end_j][axis_j],
                                )?;
                            }
                        }
                    }
                }
            }
        }
        response.potential_energy_j = finite(response.strain_energy_j - dead_work)?;
        let mut scaled_residual = Vec::with_capacity(n);
        for (node, value) in residual.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for (axis, &entry) in value.iter().enumerate() {
                if self.case.supported()[node][axis] {
                    response.support_reactions_n[node][axis] = entry;
                } else {
                    response.free_force_residual_n.push(entry);
                    scaled_residual.push(scaled(entry, 1.0 / self.scales.force_n)?);
                }
            }
        }
        for value in &mut load_column {
            *value = scaled(*value, 1.0 / self.scales.force_n)?;
        }
        checkpoint(cx)?;
        for row in 0..n {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            if tangent.row(row).1.iter().any(|value| !value.is_finite()) {
                return Err(invalid("assembled geometric truss tangent overflowed"));
            }
        }
        Ok(Analysis {
            response,
            residual: scaled_residual,
            tangent,
            load_column,
        })
    }

    /// Dimensionless residual; the physical residual remains in response().
    pub fn residual(&self, q: &[f64], lambda: f64, cx: &Cx<'_>) -> Result<Vec<f64>, SolidError> {
        Ok(self.analyze(q, lambda, cx)?.residual)
    }
    /// Exact derivative dR/dq at fixed lambda, including the deformed geometric term.
    pub fn tangent(&self, q: &[f64], lambda: f64, cx: &Cx<'_>) -> Result<Csr, SolidError> {
        Ok(self.analyze(q, lambda, cx)?.tangent)
    }
    /// Exact -dR/dlambda at fixed q, including current support-motion coupling.
    pub fn load_vector_at(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        Ok(self.analyze(q, lambda, cx)?.load_column)
    }
    pub fn response(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<GeometricResponse, SolidError> {
        Ok(self.analyze(q, lambda, cx)?.response)
    }

    /// Solve at one factor through the native Newton--Krylov owner. The seed
    /// must satisfy the exact prescribed supports. Residual tolerances act on
    /// the dimensionless force equation; no numerical residual implies stability.
    pub fn equilibrate(
        &self,
        displacement_m: &[[f64; 2]],
        lambda: f64,
        config: NewtonKrylovConfig,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<GeometricResponse, SolidError> {
        checkpoint(cx)?;
        if maximum_iterations > 10000
            || config.linear_restart > 64
            || config.max_linear_cycles > 320
        {
            return Err(invalid("geometric truss Newton work exceeds admission"));
        }
        self.admit_history(maximum_iterations, cx)?;
        let q = self.free_coordinates(displacement_m, lambda, cx)?;
        let adapter = Adapter {
            owner: self,
            cx,
            cache: RefCell::new(None),
        };
        let initial = adapter.with_analysis(&q, lambda, |a| norm(&a.residual))?;
        let fixed = FixedLoad {
            adapter: &adapter,
            lambda,
        };
        let mut state = NewtonKrylovState::new(&fixed, q, config)
            .map_err(|error| invalid(&format!("geometric truss Newton admission: {error}")))?;
        let report = state.run_cancellable(&fixed, maximum_iterations, cx);
        checkpoint(cx)?;
        let result = self.analyze(&state.x, lambda, cx)?;
        let target = config
            .absolute_tolerance
            .max(config.relative_tolerance * initial);
        if report.progress != SolverRunProgress::Complete
            || !report.report.converged
            || norm(&result.residual) > target
        {
            return Err(SolidError::SolveFailed {
                iters: report.report.iterations,
                rel_residual: norm(&result.residual) / initial.max(f64::MIN_POSITIVE),
            });
        }
        Ok(result.response)
    }

    fn admit_history(&self, entries: usize, cx: &Cx<'_>) -> Result<(), SolidError> {
        checkpoint(cx)?;
        if self
            .payload_bytes
            .checked_add(
                entries
                    .checked_mul(512)
                    .ok_or_else(|| invalid("geometric history overflow"))?,
            )
            .is_none_or(|v| v > self.limits.maximum_payload_bytes)
        {
            Err(invalid(
                "geometric truss solver history exceeds payload admission",
            ))
        } else {
            Ok(())
        }
    }

    /// Start only from an actually equilibrated physical configuration. Fixed
    /// thermal/preload baselines can first use equilibrate(). The native path
    /// stores dimensionless free displacements and the actual factor.
    pub fn start_from(
        &self,
        displacement_m: &[[f64; 2]],
        lambda: f64,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<PathState, SolidError> {
        checkpoint(cx)?;
        let mut state = PathState::start(self.free.len(), settings);
        state.u = self.free_coordinates(displacement_m, lambda, cx)?;
        state.lambda = lambda;
        state.trace = vec![(lambda, state.u.iter().map(|v| v.abs()).fold(0.0, f64::max))];
        self.advance(&state, settings, 0, cx)
    }

    /// Native pseudo-arclength across simple folds. A complete candidate is
    /// returned only after successful work; the previous checkpoint is unchanged
    /// by cancellation, domain refusal or failure of any requested step.
    /// Native event diagnostics are not certified bifurcation/stability claims.
    pub fn advance(
        &self,
        previous: &PathState,
        settings: &ArcSettings,
        steps: usize,
        cx: &Cx<'_>,
    ) -> Result<PathState, SolidError> {
        checkpoint(cx)?;
        if steps > 1000000 || settings.max_corrector_iters > 128 {
            return Err(invalid(
                "geometric truss continuation work exceeds admission",
            ));
        }
        if previous.u.len() != self.free.len() || previous.tangent.0.len() != self.free.len() {
            return Err(invalid(
                "geometric truss checkpoint has the wrong free-coordinate shape",
            ));
        }
        let histories = previous
            .trace
            .len()
            .checked_add(previous.events.len())
            .and_then(|n| n.checked_add(steps.checked_mul(2)?))
            .ok_or_else(|| invalid("geometric path history overflow"))?;
        self.admit_history(histories, cx)?;
        let adapter = Adapter {
            owner: self,
            cx,
            cache: RefCell::new(None),
        };
        let mut candidate = previous.clone();
        continuation::advance(&adapter, &mut candidate, settings, 0)?;
        if norm(&adapter.residual(&previous.u, previous.lambda)?) > settings.tol {
            return Err(invalid(
                "geometric truss path seed is not an equilibrium of its original force equation",
            ));
        }
        continuation::advance(&adapter, &mut candidate, settings, steps)?;
        checkpoint(cx)?;
        Ok(candidate)
    }
}

struct Cached {
    q: Vec<f64>,
    lambda: f64,
    analysis: Analysis,
}
struct Adapter<'a, 'g, 'c, 'cx> {
    owner: &'a GeometricTruss<'g>,
    cx: &'c Cx<'cx>,
    cache: RefCell<Option<Cached>>,
}
impl Adapter<'_, '_, '_, '_> {
    // Native callbacks at the same exact state share one physical assembly.
    fn with_analysis<T>(
        &self,
        q: &[f64],
        lambda: f64,
        read: impl FnOnce(&Analysis) -> T,
    ) -> Result<T, SolidError> {
        checkpoint(self.cx)?;
        let hit = self.cache.borrow().as_ref().is_some_and(|c| {
            c.lambda.to_bits() == lambda.to_bits()
                && c.q.len() == q.len()
                && c.q.iter().zip(q).all(|(a, b)| a.to_bits() == b.to_bits())
        });
        if !hit {
            let analysis = self.owner.analyze(q, lambda, self.cx)?;
            *self.cache.borrow_mut() = Some(Cached {
                q: q.to_vec(),
                lambda,
                analysis,
            });
        }
        Ok(read(
            &self
                .cache
                .borrow()
                .as_ref()
                .expect("exact analyzed state")
                .analysis,
        ))
    }
}
impl PathResidual for Adapter<'_, '_, '_, '_> {
    fn checkpoint(&self) -> Result<(), SolidError> {
        checkpoint(self.cx)
    }
    fn ndof(&self) -> usize {
        self.owner.free.len()
    }
    fn residual(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.with_analysis(q, lambda, |a| a.residual.clone())
    }
    fn tangent(&self, q: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        self.with_analysis(q, lambda, |a| a.tangent.clone())
    }
    fn load_vector(&self) -> Vec<f64> {
        self.owner
            .free
            .iter()
            .map(|&dof| self.owner.case.loads()[dof / 2][dof % 2] / self.owner.scales.force_n)
            .collect()
    }
    fn load_vector_at(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.with_analysis(q, lambda, |a| a.load_column.clone())
    }
}
struct FixedLoad<'a, 'b, 'g, 'c, 'cx> {
    adapter: &'a Adapter<'b, 'g, 'c, 'cx>,
    lambda: f64,
}
impl NonlinearProblem for FixedLoad<'_, '_, '_, '_, '_> {
    fn dimension(&self) -> usize {
        self.adapter.ndof()
    }
    fn residual(&self, q: &[f64], output: &mut [f64]) {
        if self
            .adapter
            .with_analysis(q, self.lambda, |a| output.copy_from_slice(&a.residual))
            .is_err()
        {
            output.fill(f64::NAN);
        }
    }
    fn jacobian_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        let result = self.adapter.with_analysis(q, self.lambda, |a| {
            for (i, value) in output.iter_mut().enumerate() {
                if i % 256 == 0 {
                    checkpoint(self.adapter.cx)?;
                }
                let (columns, entries) = a.tangent.row(i);
                *value = columns
                    .iter()
                    .zip(entries)
                    .map(|(&j, &k)| k * direction[j])
                    .sum();
            }
            Ok::<_, SolidError>(())
        });
        if !matches!(result, Ok(Ok(()))) {
            output.fill(f64::NAN);
        }
    }
    fn preconditioner_apply(
        &self,
        q: &[f64],
        _outer: usize,
        _inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        let result = self.adapter.with_analysis(q, self.lambda, |a| {
            for (i, value) in output.iter_mut().enumerate() {
                if i % 256 == 0 {
                    checkpoint(self.adapter.cx)?;
                }
                let diagonal = a.tangent.get(i, i);
                // Right Jacobi changes only the solve coordinates. A zero
                // diagonal keeps the original equation; it adds no spring.
                *value = if diagonal == 0.0 {
                    residual[i]
                } else {
                    residual[i] / diagonal
                };
            }
            Ok::<_, SolidError>(())
        });
        if !matches!(result, Ok(Ok(()))) {
            output.fill(f64::NAN);
        }
    }
}
