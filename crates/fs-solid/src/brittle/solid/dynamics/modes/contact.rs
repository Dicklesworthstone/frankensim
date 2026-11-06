//! Frozen brittle vibration with held, frictionless crack-face contact.
//!
//! Positive normal preloads select the closed-contact branch. Native paired
//! normals supply exact homogeneous displacement constraints, not springs.
//! Householder QR in fs-la and the existing fs-modal eigensolver act on the
//! original frozen stiffness and retained consistent reference mass. Reaction
//! increments have no inertia. Neither signs nor null modes are clipped.
//!
//! Preload is an explicit conditional assumption, NOT proof of equilibrium or
//! an endpoint force inferred from a sustained interval's mean force. Damage,
//! temperature and contact membership stay fixed. A caller claiming natural
//! frequencies must establish equilibrium independently. These modes do not
//! prove a finite amplitude that keeps contact loaded, resolve a switching
//! contact, include friction, geometric prestress stiffness or new collisions.

use super::{FrozenBrittleModesControls, FrozenBrittlePencil, assemble_frozen_pencil};
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::BrittleFaceImpactRow;
use crate::brittle::solid::fracture::ReleasedBrittleFace;
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::inertia::P1ReferenceMass;
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use fs_exec::Cx;
use fs_la::factor::qr;
use fs_material::brittle::{BrittleCrackBand, BrittleState};
use fs_modal::{ModalError, ModePair, eigh_gen_dense};
use fs_sparse::Csr;

/// This complete dense route is bounded before constraint reduction.
pub const MAX_CONTACT_MODE_DOFS: usize = 128;

/// Independent numerical admission gates; no value is a material stiffness.
#[derive(Clone, Copy, Debug)]
pub struct FrozenBrittleContactControls {
    pub assembly: FrozenBrittleModesControls,
    pub maximum_free_dofs: usize,
    pub maximum_contact_rows: usize,
    /// Conservative f64-slot allowance for dense scratch AND retained output.
    /// Includes full-node lifted shapes, even when almost every DOF is held.
    /// Original sparse assembly has its separate `assembly` budget.
    pub maximum_dense_entries: usize,
    pub gap_tolerance_m: f64,
    pub tangential_tolerance_m: f64,
    pub velocity_tolerance_mps: f64,
    /// A touching row must have preload strictly greater than this value.
    pub minimum_normal_preload_n: f64,
    /// Relative minimum QR diagonal after global constraint scaling, in (0,1).
    pub rank_tolerance: f64,
    /// Normwise backward-error gate for J times a basis vector or lifted mode.
    pub constraint_tolerance: f64,
    /// Normwise free-momentum residual and original-mass normalization gate.
    pub momentum_tolerance: f64,
    /// Native eigen-residual relative to the largest absolute eigenvalue.
    pub spectrum_tolerance: f64,
}

impl Default for FrozenBrittleContactControls {
    fn default() -> Self {
        Self {
            assembly: FrozenBrittleModesControls {
                maximum_dofs: 16_384,
                maximum_staged_entries: 1_048_576,
            },
            maximum_free_dofs: MAX_CONTACT_MODE_DOFS,
            maximum_contact_rows: MAX_CONTACT_MODE_DOFS,
            maximum_dense_entries: 2_097_152,
            gap_tolerance_m: 1e-10,
            tangential_tolerance_m: 1e-10,
            velocity_tolerance_mps: 1e-8,
            minimum_normal_preload_n: 1e-10,
            rank_tolerance: 1e-12,
            constraint_tolerance: 1e-10,
            momentum_tolerance: 1e-8,
            spectrum_tolerance: 1e-8,
        }
    }
}

/// Original paired trace and a physical work-conjugate force direction.
/// A positive multiplier applies -direction on A and +direction on B.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleModeConstraint {
    pub row: BrittleFaceImpactRow,
    pub direction: [f64; 3],
}

#[derive(Clone, Debug)]
pub struct FrozenBrittleContactMode {
    /// Native mass-normalized mode of the projected pencil, with signed value.
    pub coordinates: ModePair,
    /// FULL original node-major xyz shape; fixed coordinates are zero.
    pub displacement_shape: Vec<f64>,
    /// Constraint-order force increments. Multiply by the SAME modal coordinate
    /// as displacement_shape to obtain N. These are not equilibrium preloads.
    pub constraint_force_shape: Vec<f64>,
    pub constraint_relative_residual: f64,
    pub momentum_relative_residual: f64,
}

/// Original operators, admitted constraints and their complete projected modes.
#[derive(Clone, Debug)]
pub struct FrozenBrittleContactModes {
    pencil: FrozenBrittlePencil,
    constraints: Vec<BrittleModeConstraint>,
    basis: Vec<f64>,
    stiffness: Vec<f64>,
    mass: Vec<f64>,
    modes: Vec<FrozenBrittleContactMode>,
}

impl FrozenBrittleContactModes {
    pub const fn pencil(&self) -> &FrozenBrittlePencil {
        &self.pencil
    }
    /// Rows identically satisfied by the supports are omitted; their reaction
    /// increments cannot be separated from support reactions by free momentum.
    pub fn constraints(&self) -> &[BrittleModeConstraint] {
        &self.constraints
    }
    /// May be zero when independent contact rows constrain all free DOFs.
    /// In that case the projected operators and complete spectrum are empty.
    pub fn dimension(&self) -> usize {
        self.pencil.free_dofs().len() - self.constraints.len()
    }
    /// Row-major Z: free_dofs().len() rows by dimension() columns.
    pub fn basis(&self) -> &[f64] {
        &self.basis
    }
    pub fn projected_stiffness_n_m(&self) -> &[f64] {
        &self.stiffness
    }
    pub fn projected_mass_kg(&self) -> &[f64] {
        &self.mass
    }
    pub fn modes(&self) -> &[FrozenBrittleContactMode] {
        &self.modes
    }
}

#[derive(Debug)]
pub enum FrozenBrittleContactError {
    Solid(BrittleTetError),
    Modal(ModalError),
    Invalid(&'static str),
    Contact {
        face: u64,
        vertex: usize,
        reason: &'static str,
    },
    Residual {
        stage: &'static str,
        index: usize,
        relative: f64,
        tolerance: f64,
    },
}
impl From<BrittleTetError> for FrozenBrittleContactError {
    fn from(error: BrittleTetError) -> Self {
        Self::Solid(error)
    }
}
impl core::fmt::Display for FrozenBrittleContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frozen brittle contact modes refused: {self:?}")
    }
}
impl std::error::Error for FrozenBrittleContactError {}
type Error = FrozenBrittleContactError;

impl BrittleTetMotion<'_> {
    /// Observe conditional frictionless modes without changing native history.
    ///
    /// `normal_preload_n` follows retained face order and then vertex order.
    /// Tied vertices require zero preload. Strictly open pairs require exactly
    /// zero preload; touching independent pairs require strict positive preload,
    /// matched tangential geometry and negligible relative normal speed.
    /// Unloaded touching is a mode switch and refuses. Neither a supplied
    /// preload nor a previous interval's mean force establishes equilibrium.
    pub fn frozen_contact_modes(
        &self,
        normal_preload_n: &[[f64; 3]],
        fixed_dofs: &[usize],
        controls: FrozenBrittleContactControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittleContactModes, FrozenBrittleContactError> {
        assemble_contact_modes(
            self.solid.geometry,
            self.solid.budget,
            self.solid.materials.iter().map(|source| source.law()),
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.velocity_mps,
            &self.solid.released_faces,
            &self.mass,
            normal_preload_n,
            fixed_dofs,
            controls,
            cx,
        )
    }
}

/// Shared native material/mass path for mechanical and thermal owners.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_contact_modes<'m>(
    geometry: &TetLoadGeometry3<'_>,
    budget: TetAssemblyBudget,
    materials: impl ExactSizeIterator<Item = &'m BrittleCrackBand>,
    states: &[BrittleState<'_>],
    displacement_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    released_faces: &[ReleasedBrittleFace],
    mass: &P1ReferenceMass,
    normal_preload_n: &[[f64; 3]],
    fixed_dofs: &[usize],
    controls: FrozenBrittleContactControls,
    cx: &Cx<'_>,
) -> Result<FrozenBrittleContactModes, Error> {
    checkpoint(cx)?;
    let full = geometry
        .nodes_m()
        .len()
        .checked_mul(3)
        .ok_or(Error::Invalid("contact modal DOF overflow"))?;
    let free = full
        .checked_sub(fixed_dofs.len())
        .ok_or(Error::Invalid("too many fixed modal DOFs"))?;
    admit_counts(full, free, 0, controls)?;
    if displacement_m.len() != geometry.nodes_m().len()
        || velocity_mps.len() != geometry.nodes_m().len()
        || normal_preload_n.len() != released_faces.len()
    {
        return Err(Error::Invalid("contact modal source dimensions"));
    }
    let constraints = normal_constraints(
        displacement_m,
        velocity_mps,
        released_faces,
        normal_preload_n,
        controls,
        cx,
    )?;
    admit_counts(full, free, constraints.len(), controls)?;
    let pencil = assemble_frozen_pencil(
        geometry,
        budget,
        materials,
        states,
        displacement_m,
        mass,
        fixed_dofs,
        controls.assembly,
        cx,
    )?;
    project_contact_pencil(pencil, constraints, controls, cx)
}

fn admit_counts(
    full: usize,
    n: usize,
    m: usize,
    c: FrozenBrittleContactControls,
) -> Result<(), Error> {
    if full > c.assembly.maximum_dofs
        || n == 0
        || n > c.maximum_free_dofs
        || !(1..=MAX_CONTACT_MODE_DOFS).contains(&c.maximum_free_dofs)
        || !(1..=MAX_CONTACT_MODE_DOFS).contains(&c.maximum_contact_rows)
        || m > c.maximum_contact_rows
        || !c.rank_tolerance.is_finite()
        || !(0.0 < c.rank_tolerance && c.rank_tolerance < 1.0)
        || [
            c.gap_tolerance_m,
            c.tangential_tolerance_m,
            c.velocity_tolerance_mps,
            c.minimum_normal_preload_n,
            c.constraint_tolerance,
            c.momentum_tolerance,
            c.spectrum_tolerance,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(Error::Invalid("contact modal count or tolerance domain"));
    }
    // Includes native QR/eigensolver scratch, both projected operators, mode
    // coordinates and force shapes, and full-node lifted results (r <= n).
    // This conservative allocation admission is not an OS/RSS bound.
    let entries = n
        .checked_mul(n)
        .and_then(|x| x.checked_mul(48))
        .and_then(|x| {
            n.checked_mul(m)
                .and_then(|y| y.checked_mul(8))
                .and_then(|y| x.checked_add(y))
        })
        .and_then(|x| {
            full.checked_mul(n)
                .and_then(|y| y.checked_mul(2))
                .and_then(|y| x.checked_add(y))
        })
        .and_then(|x| {
            n.checked_add(m)
                .and_then(|y| y.checked_mul(32))
                .and_then(|y| x.checked_add(y))
        })
        .ok_or(Error::Invalid("contact modal dense count overflow"))?;
    if entries > c.maximum_dense_entries {
        return Err(Error::Invalid("contact modal dense storage budget"));
    }
    Ok(())
}

fn normal_constraints(
    displacement: &[[f64; 3]],
    velocity: &[[f64; 3]],
    faces: &[ReleasedBrittleFace],
    preload: &[[f64; 3]],
    c: FrozenBrittleContactControls,
    cx: &Cx<'_>,
) -> Result<Vec<BrittleModeConstraint>, Error> {
    let mut rows = Vec::new();
    for (face, forces) in faces.iter().zip(preload) {
        checkpoint(cx)?;
        for (vertex, &force) in forces.iter().enumerate() {
            let fail = |reason| Error::Contact {
                face: face.id,
                vertex,
                reason,
            };
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if !force.is_finite() || force < 0.0 {
                return Err(fail("finite nonnegative normal preload required"));
            }
            if a == b {
                if force != 0.0 {
                    return Err(fail("tied crack-front preload must be zero"));
                }
                continue;
            }
            let gap_vector = core::array::from_fn(|i| displacement[b][i] - displacement[a][i]);
            let gap = dot(face.normal, gap_vector);
            if !gap.is_finite() || gap < 0.0 {
                return Err(fail(
                    "contact modes cannot repair penetration or nonfinite geometry",
                ));
            }
            if gap > c.gap_tolerance_m.max(face.opening_tolerance_m) {
                if force != 0.0 {
                    return Err(fail("strictly open face cannot carry normal preload"));
                }
                continue;
            }
            if gap > c.gap_tolerance_m || force <= c.minimum_normal_preload_n {
                return Err(fail(
                    "touching contact needs an unambiguous positive preload",
                ));
            }
            let tangent: [f64; 3] = core::array::from_fn(|i| gap_vector[i] - gap * face.normal[i]);
            let speed = dot(
                face.normal,
                core::array::from_fn(|i| velocity[b][i] - velocity[a][i]),
            );
            if tangent.iter().any(|x| !x.is_finite())
                || tangent[0].hypot(tangent[1]).hypot(tangent[2]) > c.tangential_tolerance_m
                || !speed.is_finite()
                || speed.abs() > c.velocity_tolerance_mps
            {
                return Err(fail(
                    "loaded paired faces must be matched and normally resting",
                ));
            }
            if rows.len() == c.maximum_contact_rows {
                return Err(Error::Invalid(
                    "all loaded contact rows must fit the modal budget",
                ));
            }
            rows.push(BrittleModeConstraint {
                row: BrittleFaceImpactRow {
                    face_id: face.id,
                    vertex,
                    nodes: [a, b],
                },
                direction: face.normal,
            });
        }
    }
    Ok(rows)
}

/// Project physical trace directions separately from their contact admission.
/// Future sticking constraints can enter the SAME K/M projection, without
/// substituting tangent stiffness, contact mass or a new numerical solver.
pub(crate) fn project_contact_pencil(
    pencil: FrozenBrittlePencil,
    constraints: Vec<BrittleModeConstraint>,
    c: FrozenBrittleContactControls,
    cx: &Cx<'_>,
) -> Result<FrozenBrittleContactModes, Error> {
    checkpoint(cx)?;
    let free = pencil.free_dofs();
    let n = free.len();
    admit_counts(pencil.full_dofs(), n, constraints.len(), c)?;
    let mut rows = Vec::with_capacity(constraints.len());
    for constraint in constraints {
        checkpoint(cx)?;
        if constraint
            .row
            .nodes
            .iter()
            .any(|&node| node >= pencil.full_dofs() / 3)
            || constraint.direction.iter().any(|x| !x.is_finite())
        {
            return Err(Error::Invalid("invalid physical modal constraint"));
        }
        if (0..3).any(|axis| {
            constraint.direction[axis] != 0.0
                && constraint.row.nodes[0] != constraint.row.nodes[1]
                && constraint
                    .row
                    .nodes
                    .iter()
                    .any(|&node| free.binary_search(&(3 * node + axis)).is_ok())
        }) {
            rows.push(constraint);
        }
    }
    let m = rows.len();
    if m > n {
        return Err(Error::Invalid("unresolved contact constraint rank"));
    }
    let r = n - m;
    let mut transpose = vec![0.0; n * m];
    for (column, constraint) in rows.iter().enumerate() {
        checkpoint(cx)?;
        for (side, &node) in constraint.row.nodes.iter().enumerate() {
            for axis in 0..3 {
                if let Ok(row) = free.binary_search(&(3 * node + axis)) {
                    transpose[row * m + column] +=
                        (if side == 0 { -1.0 } else { 1.0 }) * constraint.direction[axis];
                }
            }
        }
    }
    let scale = transpose.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    let mut scaled = transpose.clone();
    if m > 0 {
        if !scale.is_finite() || scale == 0.0 {
            return Err(Error::Invalid("invalid contact constraint scale"));
        }
        for value in &mut scaled {
            *value /= scale;
        }
    }
    checkpoint(cx)?;
    let factor = qr(&scaled, n, m);
    for i in 0..m {
        if !factor.r(i, i).is_finite() || factor.r(i, i).abs() <= c.rank_tolerance {
            return Err(Error::Invalid(
                "unresolved contact constraint rank; no row was dropped",
            ));
        }
    }
    let mut basis = vec![0.0; n * r];
    for column in 0..r {
        checkpoint(cx)?;
        let mut vector = vec![0.0; n];
        vector[m + column] = 1.0;
        factor.apply_q(&mut vector);
        gate(
            "contact basis",
            column,
            constraint_residual(&transpose, n, m, &vector)?,
            c.constraint_tolerance,
        )?;
        for row in 0..n {
            basis[row * r + column] = vector[row];
        }
    }
    let stiffness = project(pencil.stiffness_n_m(), &basis, r, cx)?;
    let mass = project(pencil.mass_kg(), &basis, r, cx)?;
    checkpoint(cx)?;
    let native = eigh_gen_dense(&stiffness, &mass, r).map_err(Error::Modal)?;
    checkpoint(cx)?;
    if native.len() != r {
        return Err(Error::Invalid("incomplete contact modal basis"));
    }
    let spectrum_scale = native
        .iter()
        .map(|mode| mode.lambda.abs())
        .fold(0.0_f64, f64::max);
    let k_norm = operator_norm(pencil.stiffness_n_m())?;
    let m_norm = operator_norm(pencil.mass_kg())?;
    let jt_norm = (0..n)
        .map(|i| {
            transpose[i * m..(i + 1) * m]
                .iter()
                .map(|x| x.abs())
                .sum::<f64>()
        })
        .fold(0.0_f64, f64::max);
    let mut modes = Vec::with_capacity(r);
    for (index, coordinates) in native.into_iter().enumerate() {
        checkpoint(cx)?;
        if coordinates.phi.len() != r
            || !coordinates.lambda.is_finite()
            || !coordinates.residual.is_finite()
            || coordinates.residual < 0.0
            || !coordinates.interval.0.is_finite()
            || !coordinates.interval.1.is_finite()
            || coordinates.phi.iter().any(|x| !x.is_finite())
        {
            return Err(Error::Invalid("nonfinite native contact eigenpair"));
        }
        gate(
            "projected spectrum",
            index,
            relative(coordinates.residual, spectrum_scale)?,
            c.spectrum_tolerance,
        )?;
        let phi: Vec<_> = (0..n)
            .map(|i| {
                (0..r)
                    .map(|j| basis[i * r + j] * coordinates.phi[j])
                    .sum::<f64>()
            })
            .collect();
        let constraint_relative_residual = constraint_residual(&transpose, n, m, &phi)?;
        gate(
            "contact mode",
            index,
            constraint_relative_residual,
            c.constraint_tolerance,
        )?;
        let k_phi = action(pencil.stiffness_n_m(), &phi)?;
        let m_phi = action(pencil.mass_kg(), &phi)?;
        let mass_norm: f64 = phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum();
        gate(
            "original mass normalization",
            index,
            (mass_norm - 1.0).abs(),
            c.momentum_tolerance,
        )?;
        // Positive contact force opens the gap: K phi - lambda M phi = J^T df.
        let rhs: Vec<_> = k_phi
            .iter()
            .zip(&m_phi)
            .map(|(&k, &mass)| (-coordinates.lambda).mul_add(mass, k))
            .collect();
        vector_norm(&rhs)?;
        let mut constraint_force_shape = factor.solve_ls(&rhs);
        if m > 0 {
            for force in &mut constraint_force_shape {
                *force /= scale;
            }
        }
        let force_norm = vector_norm(&constraint_force_shape)?;
        let mut residual = 0.0_f64;
        for i in 0..n {
            let reaction: f64 = (0..m)
                .map(|j| transpose[i * m + j] * constraint_force_shape[j])
                .sum();
            if !reaction.is_finite() {
                return Err(Error::Invalid("nonfinite contact reaction action"));
            }
            residual = residual.max((rhs[i] - reaction).abs());
        }
        let denominator = (k_norm + coordinates.lambda.abs() * m_norm) * vector_norm(&phi)?
            + jt_norm * force_norm;
        let momentum_relative_residual = relative(residual, denominator)?;
        gate(
            "free momentum",
            index,
            momentum_relative_residual,
            c.momentum_tolerance,
        )?;
        let mut displacement_shape = vec![0.0; pencil.full_dofs()];
        for (&dof, &value) in free.iter().zip(&phi) {
            displacement_shape[dof] = value;
        }
        modes.push(FrozenBrittleContactMode {
            coordinates,
            displacement_shape,
            constraint_force_shape,
            constraint_relative_residual,
            momentum_relative_residual,
        });
    }
    checkpoint(cx)?;
    Ok(FrozenBrittleContactModes {
        pencil,
        constraints: rows,
        basis,
        stiffness,
        mass,
        modes,
    })
}

fn project(matrix: &Csr, basis: &[f64], r: usize, cx: &Cx<'_>) -> Result<Vec<f64>, Error> {
    let n = matrix.nrows();
    let mut action = vec![0.0; n * r];
    for i in 0..n {
        checkpoint(cx)?;
        let (cols, values) = matrix.row(i);
        for (&j, &value) in cols.iter().zip(values) {
            for k in 0..r {
                action[i * r + k] = value.mul_add(basis[j * r + k], action[i * r + k]);
            }
        }
    }
    let mut result = vec![0.0; r * r];
    for i in 0..r {
        checkpoint(cx)?;
        for j in 0..r {
            for k in 0..n {
                result[i * r + j] = basis[k * r + i].mul_add(action[k * r + j], result[i * r + j]);
            }
        }
    }
    if result.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("nonfinite contact projection"));
    }
    Ok(result)
}
fn action(matrix: &Csr, vector: &[f64]) -> Result<Vec<f64>, Error> {
    let result: Vec<_> = (0..matrix.nrows())
        .map(|i| {
            let (cols, values) = matrix.row(i);
            cols.iter()
                .zip(values)
                .map(|(&j, &value)| value * vector[j])
                .sum::<f64>()
        })
        .collect();
    vector_norm(&result)?;
    Ok(result)
}
fn operator_norm(matrix: &Csr) -> Result<f64, Error> {
    let rows: Vec<_> = (0..matrix.nrows())
        .map(|i| matrix.row(i).1.iter().map(|x| x.abs()).sum())
        .collect();
    vector_norm(&rows)
}
fn constraint_residual(
    transpose: &[f64],
    n: usize,
    m: usize,
    vector: &[f64],
) -> Result<f64, Error> {
    let mut residual = 0.0_f64;
    let mut norm = 0.0_f64;
    for row in 0..m {
        let value: f64 = (0..n).map(|i| transpose[i * m + row] * vector[i]).sum();
        let magnitude: f64 = (0..n).map(|i| transpose[i * m + row].abs()).sum();
        if !value.is_finite() || !magnitude.is_finite() {
            return Err(Error::Invalid("nonfinite contact constraint residual"));
        }
        residual = residual.max(value.abs());
        norm = norm.max(magnitude);
    }
    relative(residual, norm * vector_norm(vector)?)
}
fn vector_norm(values: &[f64]) -> Result<f64, Error> {
    if values.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("nonfinite contact mode arithmetic"));
    }
    Ok(values.iter().map(|x| x.abs()).fold(0.0_f64, f64::max))
}
fn relative(error: f64, scale: f64) -> Result<f64, Error> {
    if !error.is_finite() || !scale.is_finite() || error < 0.0 || scale < 0.0 {
        return Err(Error::Invalid("nonfinite contact residual scale"));
    }
    Ok(if scale == 0.0 {
        if error == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        error / scale
    })
}
fn gate(stage: &'static str, index: usize, relative: f64, tolerance: f64) -> Result<(), Error> {
    if !relative.is_finite() || relative > tolerance {
        return Err(Error::Residual {
            stage,
            index,
            relative,
            tolerance,
        });
    }
    Ok(())
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
