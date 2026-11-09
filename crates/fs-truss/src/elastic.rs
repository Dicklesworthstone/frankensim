//! Displacement-compatible elasticity on the native pin-jointed ground.
//!
//! Plastic layout forces need only satisfy nodal equilibrium. After selecting
//! areas, redundant members instead share load according to their actual EA/L.
//! This owner solves K u = f with K = B diag(EA/L) B^T, applies that operator
//! without assembling a global matrix, and recovers N = EA/L B^T u. All member
//! directions, supports and loads retain their original ground ordering.
//!
//! This is infinitesimal axial elasticity with fixed translational supports.
//! It contains no bending, prestress stiffness, buckling or plastic return map.
//! Acceptance measures equilibrium of the returned compatible displacement;
//! it does not prove uniqueness or exclude unloaded mechanisms. Zero-stiffness
//! unloaded coordinates receive zero displacement, never an invented spring.

use crate::{GroundStructure, LayoutCase, TrussConstructionError};
use fs_exec::Cx;
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, checked_cg};
use fs_sparse::precond::Precond;

/// Accuracy and bounded work for the native matrix-free elastic solve.
#[derive(Debug, Clone, Copy)]
pub struct ElasticControls {
    /// Relative Euclidean residual of actual member forces on free DOFs.
    pub relative_tolerance: f64,
    /// Total native CG iterations, including at most two defect corrections.
    pub maximum_iterations: usize,
}

impl Default for ElasticControls {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-10,
            maximum_iterations: 4096,
        }
    }
}

/// No accepted response is returned on malformed data or an unresolved load.
#[derive(Debug)]
pub enum ElasticError {
    /// Original geometry/material/control admission or cancellation.
    Construction(TrussConstructionError),
    /// A load acts on a free coordinate with no axial elastic stiffness.
    UnsupportedLoad {
        /// Original node index.
        node: usize,
        /// Global translation component, x = 0 and y = 1.
        axis: usize,
    },
    /// Native CG or the recovered-force residual did not converge.
    NotConverged {
        /// Measured relative Euclidean equilibrium error.
        relative_residual: f64,
        /// Total native CG iterations used.
        iterations: usize,
    },
}

impl From<TrussConstructionError> for ElasticError {
    fn from(error: TrussConstructionError) -> Self {
        Self::Construction(error)
    }
}
impl core::fmt::Display for ElasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "truss elastic response: {self:?}")
    }
}
impl std::error::Error for ElasticError {}

/// One compatible, numerically equilibrated linear response in coherent SI.
#[derive(Debug, Clone)]
pub struct ElasticResponse {
    /// Original node ordering; supported coordinates are exactly zero [m].
    pub displacements_m: Vec<[f64; 2]>,
    /// Original member ordering; zero-area members return zero strain/force/stress.
    pub member_strain: Vec<f64>,
    /// Positive tension, negative compression [N].
    pub member_force_n: Vec<f64>,
    /// Signed axial Cauchy stress to infinitesimal-strain accuracy [Pa].
    pub member_stress_pa: Vec<f64>,
    /// K u - f on supported coordinates, zero on free coordinates [N].
    /// Includes external loads applied directly to supports.
    pub support_reactions_n: Vec<[f64; 2]>,
    /// Sum of 1/2 N delta_length over actual members [J].
    pub strain_energy_j: f64,
    /// External nodal load dotted with compatible displacement [J].
    pub compliance_j: f64,
    /// Recomputed from recovered forces on ALL original free coordinates.
    pub relative_equilibrium_residual: f64,
    /// Total native CG iteration count.
    pub iterations: usize,
}

fn invalid(field: &'static str, requirement: &'static str) -> TrussConstructionError {
    TrussConstructionError::InvalidInput { field, requirement }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TrussConstructionError> {
    cx.checkpoint().map_err(|_| TrussConstructionError::Cancelled {
        stage: "truss elastic response",
    })
}

struct Member {
    ends: [usize; 4],
    direction: [f64; 4],
    stiffness: f64,
}

struct Operator<'a, 'cx> {
    members: Vec<Member>,
    dof_map: Vec<Option<usize>>,
    diagonal: Vec<f64>,
    cx: &'a Cx<'cx>,
}

impl LinearOp for Operator<'_, '_> {
    fn n(&self) -> usize {
        self.diagonal.len()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for (index, member) in self.members.iter().enumerate() {
            // The legacy action has no error channel. Nonfinite output stops
            // CG; the surrounding checkpoint returns the actual cancellation.
            if index % 256 == 0 && checkpoint(self.cx).is_err() {
                y.fill(f64::NAN);
                return;
            }
            if member.stiffness == 0.0 {
                continue;
            }
            let extension = (0..4)
                .map(|i| {
                    self.dof_map[member.ends[i]]
                        .map_or(0.0, |dof| member.direction[i] * x[dof])
                })
                .sum::<f64>();
            let force = member.stiffness * extension;
            for i in 0..4 {
                if let Some(dof) = self.dof_map[member.ends[i]] {
                    y[dof] += member.direction[i] * force;
                }
            }
        }
    }
}

impl Precond for Operator<'_, '_> {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for (i, ((value, &rhs), &diagonal)) in
            z.iter_mut().zip(r).zip(&self.diagonal).enumerate()
        {
            if i % 256 == 0 && checkpoint(self.cx).is_err() {
                *value = f64::NAN;
                return;
            }
            *value = rhs / diagonal;
        }
    }
}

/// Solve the actual selected areas using each member's explicit Young's modulus.
///
/// Arrays follow the complete original member order; area zero removes a
/// member, while every supplied modulus must remain finite and positive. The
/// same source geometry and case may therefore be reused after catalog sizing.
/// Storage and every operator action are O(nodes + members), with deterministic
/// member-order reductions. Member/recovery loops poll cancellation every 256
/// entries; native CG retains its own bounded iteration checkpoints.
///
/// The caller owns the validity of small strain and any stability/uniqueness
/// claim. An unloaded zero-stiffness coordinate receives zero motion. A loaded
/// one refuses; other loaded mechanisms refuse when true equilibrium fails.
/// No output is labeled converged from CG's recursive residual alone.
///
/// # Errors
/// Invalid arrays, material/controls, nonfinite arithmetic, unresolved loads
/// and cancellation return no accepted response.
#[allow(clippy::too_many_lines)] // Admission, native solve, and physical recovery.
pub fn solve(
    ground: &GroundStructure,
    case: &LayoutCase,
    areas_m2: &[f64],
    young_pa: &[f64],
    controls: ElasticControls,
    cx: &Cx<'_>,
) -> Result<ElasticResponse, ElasticError> {
    checkpoint(cx)?;
    let count = ground.members().len();
    for (field, actual) in [
        ("elastic areas", areas_m2.len()),
        ("elastic moduli", young_pa.len()),
    ] {
        if actual != count {
            return Err(TrussConstructionError::VectorLength {
                field,
                expected: count,
                actual,
            }
            .into());
        }
    }
    if case.node_count() != ground.nodes().len() {
        return Err(invalid("elastic load case", "must match original node ordering").into());
    }
    if !(controls.relative_tolerance > 0.0 && controls.relative_tolerance < 1.0)
        || controls.maximum_iterations > 1_000_000
    {
        return Err(invalid(
            "elastic controls",
            "require tolerance in (0,1) and at most 1000000 iterations",
        )
        .into());
    }
    let mut diagonal = vec![0.0; 2 * ground.nodes().len()];
    let mut members = Vec::with_capacity(count);
    for (index, &(a, b)) in ground.members().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        let area = areas_m2[index];
        let young = young_pa[index];
        if !area.is_finite() || area < 0.0 || !young.is_finite() || young <= 0.0 {
            return Err(invalid(
                "elastic material",
                "requires finite nonnegative areas and positive moduli",
            )
            .into());
        }
        let length = ground.lengths()[index];
        let stiffness = young * (area / length);
        if !stiffness.is_finite() || (area > 0.0 && stiffness <= 0.0) {
            return Err(invalid(
                "elastic axial stiffness",
                "EA/L must be positive and representable for present members",
            )
            .into());
        }
        let d = [
            (ground.nodes()[b][0] - ground.nodes()[a][0]) / length,
            (ground.nodes()[b][1] - ground.nodes()[a][1]) / length,
        ];
        let member = Member {
            ends: [2 * a, 2 * a + 1, 2 * b, 2 * b + 1],
            direction: [-d[0], -d[1], d[0], d[1]],
            stiffness,
        };
        for i in 0..4 {
            diagonal[member.ends[i]] += stiffness * member.direction[i].powi(2);
        }
        members.push(member);
    }
    let mut dof_map = vec![None; diagonal.len()];
    let mut free_diagonal = Vec::new();
    let mut rhs = Vec::new();
    for (node, support) in case.supported().iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for (axis, &fixed) in support.iter().enumerate() {
            let d = diagonal[2 * node + axis];
            if !d.is_finite() {
                return Err(invalid("elastic diagonal", "must remain finite").into());
            }
            if fixed {
                continue;
            }
            if d == 0.0 {
                if case.loads()[node][axis] != 0.0 {
                    return Err(ElasticError::UnsupportedLoad { node, axis });
                }
                continue;
            }
            dof_map[2 * node + axis] = Some(rhs.len());
            free_diagonal.push(d);
            rhs.push(case.loads()[node][axis]);
        }
    }
    let operator = Operator {
        members,
        dof_map,
        diagonal: free_diagonal,
        cx,
    };
    let solution = checked_cg(
        &operator,
        &operator,
        &rhs,
        CheckedCgConfig {
            tolerance: controls.relative_tolerance,
            max_iterations: controls.maximum_iterations,
            max_corrections: 2,
        },
        |_| checkpoint(cx),
    )
    .map_err(|error| match error {
        CheckedCgError::Interrupted(error) => ElasticError::Construction(error),
        CheckedCgError::InvalidInput(_) => {
            invalid("elastic solve", "native solver input refused").into()
        }
    })?;
    if !solution.report.converged_euclidean() {
        return Err(ElasticError::NotConverged {
            relative_residual: solution.report.rel_residual,
            iterations: solution.report.iters,
        });
    }
    let mut displacement = vec![[0.0; 2]; ground.nodes().len()];
    for (node, u) in displacement.iter_mut().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for (axis, value) in u.iter_mut().enumerate() {
            *value = operator.dof_map[2 * node + axis].map_or(0.0, |i| solution.x[i]);
        }
    }
    let mut forces = Vec::with_capacity(count);
    let mut strains = Vec::with_capacity(count);
    let mut stresses = Vec::with_capacity(count);
    let mut internal = vec![[0.0; 2]; displacement.len()];
    let mut energy = 0.0;
    for (index, member) in operator.members.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if member.stiffness == 0.0 {
            forces.push(0.0);
            strains.push(0.0);
            stresses.push(0.0);
            continue;
        }
        let extension = (0..4)
            .map(|i| {
                let slot = member.ends[i];
                member.direction[i] * displacement[slot / 2][slot % 2]
            })
            .sum::<f64>();
        let force = member.stiffness * extension;
        let strain = extension / ground.lengths()[index];
        let stress = young_pa[index] * strain;
        energy += 0.5 * force * extension;
        forces.push(force);
        strains.push(strain);
        stresses.push(stress);
        for i in 0..4 {
            let slot = member.ends[i];
            internal[slot / 2][slot % 2] += member.direction[i] * force;
        }
        if [force, strain, stress, energy].iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "elastic recovery",
                "member force, strain, stress and energy must remain finite",
            )
            .into());
        }
    }
    let mut reactions = vec![[0.0; 2]; displacement.len()];
    let mut residual = 0.0_f64;
    let mut load_norm = 0.0_f64;
    // Normalize first: overflowing absolute load norms cannot hide defects.
    let load_scale = rhs.iter().map(|v| v.abs()).fold(0.0, f64::max);
    let load_scale = if load_scale == 0.0 { 1.0 } else { load_scale };
    let mut compliance = 0.0;
    for node in 0..displacement.len() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        for axis in 0..2 {
            let load = case.loads()[node][axis];
            let defect = internal[node][axis] - load;
            if case.supported()[node][axis] {
                reactions[node][axis] = defect;
            } else {
                residual = residual.hypot(defect / load_scale);
                load_norm = load_norm.hypot(load / load_scale);
            }
            compliance += load * displacement[node][axis];
            if !defect.is_finite() || !compliance.is_finite() {
                return Err(invalid(
                    "elastic recovery",
                    "nodal reactions and compliance must remain finite",
                )
                .into());
            }
        }
    }
    let relative = if load_norm == 0.0 { residual } else { residual / load_norm };
    if !relative.is_finite() || relative >= controls.relative_tolerance {
        return Err(ElasticError::NotConverged {
            relative_residual: relative,
            iterations: solution.report.iters,
        });
    }
    checkpoint(cx)?;
    Ok(ElasticResponse {
        displacements_m: displacement,
        member_strain: strains,
        member_force_n: forces,
        member_stress_pa: stresses,
        support_reactions_n: reactions,
        strain_energy_j: energy,
        compliance_j: compliance,
        relative_equilibrium_residual: relative,
        iterations: solution.report.iters,
    })
}
