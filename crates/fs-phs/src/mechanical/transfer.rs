//! Instantaneous canonical changes with every original Maxwell history retained.

use super::LinearMechanicalSystem;
use crate::PhsError;
use crate::interconnection::accumulation::Sum;
use fs_la::factor::{FactorError, lu};
use fs_math::det;

/// Complete, fixed-frame coordinate correspondence; no modal ordering is used.
#[derive(Clone, Copy, Debug)]
pub enum MechanicalCoordinateMap<'a> {
    /// The same physical coordinates; displacement and momentum retain bits.
    Unchanged,
    /// Invertible row-major A, with q_target=A*q_previous. Canonical momenta
    /// transform by A^-T, independently of either body's physical inertia.
    TargetFromPrevious(&'a [f64]),
}

/// Cold allocation/conditioning and warm physical-event admission limits.
#[derive(Clone, Copy, Debug)]
pub struct MechanicalTransferControls {
    /// Additional cold-work bound; native construction already requires
    /// 2*coordinates+branches <= 256, hence coordinates <= 128.
    pub maximum_coordinates: usize,
    /// Conservative numeric payload of preparation and one workspace, including
    /// both retained original relaxed-K matrices. Caller input/output states
    /// (including the convenience method's returned state) and other pre-existing
    /// source storage are excluded. No count depends on rank or damping.
    pub maximum_payload_bytes: usize,
    /// Upper limit on ||A||_infinity*||A^-1||_infinity, at least one. Coordinates
    /// retain their authored units; admission never normalizes the physical map.
    pub maximum_condition_number: f64,
    /// Componentwise backward residual for the original inverse, displacement,
    /// momentum and branch-history reconstructions; finite in [0,1e-6].
    pub maximum_relative_mapping_residual: f64,
    /// Componentwise error in b_target^T*A=b_previous^T, relative to the sum
    /// of absolute source terms. Finite in [0,1e-6]. A changed physical strain
    /// law is not a coordinate change and needs its own material transport law.
    pub maximum_relative_strain_residual: f64,
    /// Maximum magnitude of the signed external Hamiltonian parameter work.
    pub maximum_absolute_parameter_work_j: f64,
    /// Absolute allowance on actual storage increment minus parameter work.
    /// This does not enlarge the work allowance or authorize history loss.
    pub energy_tolerance_j: f64,
}

/// All failures leave caller state, target output and both sources unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MechanicalTransferError {
    Invalid(&'static str),
    Source(PhsError),
    Factor(FactorError),
    Cancelled,
    Budget {
        requested: usize,
        maximum: usize,
    },
    Limit {
        quantity: &'static str,
        value: f64,
        maximum: f64,
    },
}
impl core::fmt::Display for MechanicalTransferError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "mechanical history transfer refused: {self:?}")
    }
}
impl std::error::Error for MechanicalTransferError {}
impl From<PhsError> for MechanicalTransferError {
    fn from(error: PhsError) -> Self {
        Self::Source(error)
    }
}
impl From<FactorError> for MechanicalTransferError {
    fn from(error: FactorError) -> Self {
        Self::Factor(error)
    }
}
type Result<T> = core::result::Result<T, MechanicalTransferError>;

/// Independent complete storage and physical parameter-work accounts [J].
#[derive(Clone, Copy, Debug)]
pub struct MechanicalTransferReport {
    pub previous_stored_energy_j: f64,
    pub target_stored_energy_j: f64,
    pub kinetic_parameter_work_j: f64,
    pub elastic_parameter_work_j: f64,
    pub memory_parameter_work_j: f64,
    /// Sum before rounding its separate observable channels. This is external
    /// parameter work, never force work, dissipative heat or projection loss.
    pub parameter_work_j: f64,
    /// Increment contracted from original source operators and actual stored
    /// states without subtracting their possibly large total Hamiltonians.
    pub storage_change_j: f64,
    pub energy_residual_j: f64,
    pub maximum_relative_state_residual: f64,
}

#[derive(Clone, Debug)]
pub struct MechanicalStateTransfer {
    pub state_after: Vec<f64>,
    pub report: MechanicalTransferReport,
}

/// A prepared map bound to the exact two immutable mechanical source owners.
/// All original coordinates and Maxwell arms remain. No eigenvectors, modal
/// identifiers, truncation, diagonal damping or second integrator are involved.
pub struct PreparedMechanicalTransfer<'source> {
    previous: &'source LinearMechanicalSystem,
    target: &'source LinearMechanicalSystem,
    // A and A^-T. Identity is kept structural, including signed zero state.
    maps: Option<(Vec<f64>, Vec<f64>)>,
    previous_branches: Vec<usize>,
    stiffness_change: Vec<f64>,
    inverse_mass_change: Vec<f64>,
    controls: MechanicalTransferControls,
    condition_number: f64,
    maximum_map_residual: f64,
    maximum_strain_residual: f64,
}

/// Reusable scratch intrinsically bound to one prepared source pair and map.
/// A warm transfer allocates nothing and performs no factorization or solve.
pub struct MechanicalTransferWorkspace<'plan, 'source> {
    plan: &'plan PreparedMechanicalTransfer<'source>,
    candidate: Vec<f64>,
}

impl LinearMechanicalSystem {
    /// Prepare a zero-duration parameter/basis change without erasing memory.
    ///
    /// The caller supplies an explicit bijection: entry j identifies which OLD
    /// physical Maxwell branch becomes target branch j. Stiffness and relaxation
    /// time may change, but equal times do not identify or merge histories. The
    /// source-binding caller establishes physical branch identity; this numerical
    /// owner checks the complete permutation and the original strain covectors.
    ///
    /// q_new=A*q_old and p_new=A^-T*p_old retain canonical motion. Changed mass
    /// therefore changes velocity. Each viscous strain survives, so z_new is
    /// sqrt(k_new)/sqrt(k_old)*z_old, never an unscaled copy or relaxed reset.
    /// We require b_new^T*A=b_old^T: arbitrary physical geometry/strain-law changes,
    /// branch births/deaths, remeshing and constraints are deliberately refused.
    /// Full-rank rotations through degenerate modes require no mode matching.
    ///
    /// The target's full M/K/C/P and relaxation times govern subsequent original
    /// pHS evolution. Rebuild propagation/observation owners from that target at
    /// the SAME accepted physical instant. This operation advances no time and
    /// supplies no force impulse or heat. The enclosing source/event owner must
    /// fund the reported signed parameter work and publish the complete update.
    ///
    /// Cold work is bounded O(n^3), using the existing partial-pivot LU. Warm
    /// work is O(n^2+n*branches). Cancellation is checked between bounded rows;
    /// the native bounded LU itself is a single noncancellable cold operation
    /// with n <= 128 under the existing complete 256-state source ceiling.
    pub fn prepare_state_transfer<'source>(
        &'source self,
        target: &'source Self,
        coordinates: MechanicalCoordinateMap<'_>,
        previous_branch_for_target: &[usize],
        controls: MechanicalTransferControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalTransfer<'source>> {
        poll(&mut cancelled)?;
        let n = self.coordinates;
        let branches = self.branches.len();
        if target.coordinates != n
            || target.branches.len() != branches
            || previous_branch_for_target.len() != branches
        {
            return Err(invalid(
                "complete equal-dimensional coordinates and branch bijection required",
            ));
        }
        admit_controls(controls)?;
        if n > controls.maximum_coordinates {
            return Err(MechanicalTransferError::Budget {
                requested: n,
                maximum: controls.maximum_coordinates,
            });
        }
        // Includes both original K arrays, retained maps/differences, LU and
        // congruence temporaries, branch permutation, and one candidate state.
        // Native construction already bounds n and the complete state by 256.
        let bytes = (16 * n * n + 8 * n + 4 * branches) * size_of::<f64>()
            + 2 * branches * size_of::<usize>();
        if bytes > controls.maximum_payload_bytes {
            return Err(MechanicalTransferError::Budget {
                requested: bytes,
                maximum: controls.maximum_payload_bytes,
            });
        }
        let mut seen = vec![false; branches];
        for &old in previous_branch_for_target {
            if old >= branches || seen[old] {
                return Err(invalid(
                    "Maxwell correspondence must preserve every original branch exactly once",
                ));
            }
            seen[old] = true;
        }
        let (maps, condition_number, maximum_map_residual) = match coordinates {
            MechanicalCoordinateMap::Unchanged => (None, 1.0, 0.0),
            MechanicalCoordinateMap::TargetFromPrevious(a) => {
                if a.len() != n * n || a.iter().any(|x| !x.is_finite()) {
                    return Err(invalid(
                        "complete finite full-rank physical coordinate map required",
                    ));
                }
                let factor = lu(a, n)?;
                poll(&mut cancelled)?;
                let mut inverse_transpose = vec![0.0; n * n];
                let mut column = vec![0.0; n];
                for j in 0..n {
                    poll(&mut cancelled)?;
                    column.fill(0.0);
                    column[j] = 1.0;
                    factor.solve_transpose(&mut column);
                    for i in 0..n {
                        inverse_transpose[i * n + j] = checked(column[i])?;
                    }
                }
                let mut maximum = 0.0_f64;
                let (mut a_norm, mut inverse_norm) = (0.0_f64, 0.0_f64);
                for i in 0..n {
                    poll(&mut cancelled)?;
                    let mut a_row = Sum::default();
                    let mut inverse_row = Sum::default();
                    for j in 0..n {
                        a_row
                            .product(a[i * n + j].abs(), 1.0)
                            .ok_or_else(arithmetic)?;
                        inverse_row
                            .product(inverse_transpose[j * n + i].abs(), 1.0)
                            .ok_or_else(arithmetic)?;
                        maximum = maximum.max(relative_residual(
                            (0..n).map(|k| (a[k * n + i], inverse_transpose[k * n + j])),
                            if i == j { 1.0 } else { 0.0 },
                        )?);
                    }
                    a_norm = a_norm.max(value(a_row)?);
                    inverse_norm = inverse_norm.max(value(inverse_row)?);
                }
                let condition = checked(a_norm * inverse_norm)?;
                limit(
                    "coordinate condition number",
                    condition,
                    controls.maximum_condition_number,
                )?;
                limit(
                    "coordinate inverse residual",
                    maximum,
                    controls.maximum_relative_mapping_residual,
                )?;
                (Some((a.to_vec(), inverse_transpose)), condition, maximum)
            }
        };
        let mut maximum_strain_residual = 0.0_f64;
        for (j, &old) in previous_branch_for_target.iter().enumerate() {
            poll(&mut cancelled)?;
            let before = &self.branches[old].projection;
            let after = &target.branches[j].projection;
            for i in 0..n {
                let residual = if let Some((a, _)) = &maps {
                    relative_residual((0..n).map(|k| (after[k], a[k * n + i])), before[i])?
                } else {
                    relative_residual([(after[i], 1.0)], before[i])?
                };
                maximum_strain_residual = maximum_strain_residual.max(residual);
            }
        }
        limit(
            "original branch strain correspondence",
            maximum_strain_residual,
            controls.maximum_relative_strain_residual,
        )?;
        let stiffness_change = pullback_difference(
            &self.relaxed_stiffness,
            &target.relaxed_stiffness,
            maps.as_ref().map(|(a, _)| a.as_slice()),
            n,
            &mut cancelled,
        )?;
        let inverse_mass_change = pullback_difference(
            &self.inverse_mass,
            &target.inverse_mass,
            maps.as_ref().map(|(_, b)| b.as_slice()),
            n,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedMechanicalTransfer {
            previous: self,
            target,
            maps,
            previous_branches: previous_branch_for_target.to_vec(),
            stiffness_change,
            inverse_mass_change,
            controls,
            condition_number,
            maximum_map_residual,
            maximum_strain_residual,
        })
    }
}

impl<'source> PreparedMechanicalTransfer<'source> {
    pub fn previous(&self) -> &'source LinearMechanicalSystem {
        self.previous
    }
    pub fn target(&self) -> &'source LinearMechanicalSystem {
        self.target
    }
    pub fn condition_number(&self) -> f64 {
        self.condition_number
    }
    pub fn maximum_map_residual(&self) -> f64 {
        self.maximum_map_residual
    }
    pub fn maximum_strain_residual(&self) -> f64 {
        self.maximum_strain_residual
    }
    pub fn workspace(&self) -> MechanicalTransferWorkspace<'_, 'source> {
        MechanicalTransferWorkspace {
            plan: self,
            candidate: vec![0.0; self.target.system.state_dim()],
        }
    }
    /// Allocating convenience; use the bound workspace for repeated transfers.
    pub fn transfer(
        &self,
        state: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MechanicalStateTransfer> {
        let mut workspace = self.workspace();
        let mut state_after = vec![0.0; self.target.system.state_dim()];
        let report = workspace.transfer_into(state, &mut state_after, cancelled)?;
        Ok(MechanicalStateTransfer {
            state_after,
            report,
        })
    }
}

impl MechanicalTransferWorkspace<'_, '_> {
    /// Publish complete target state only after motion, history and independent
    /// work/storage checks pass. Output is unchanged on every refusal, including
    /// late cancellation. Scratch may change and is reusable after refusal.
    pub fn transfer_into(
        &mut self,
        state: &[f64],
        target_state: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MechanicalTransferReport> {
        poll(&mut cancelled)?;
        let p = self.plan;
        p.previous.admit_state(state)?;
        if target_state.len() != self.candidate.len() {
            return Err(invalid("complete target mechanical state output required"));
        }
        let n = p.previous.coordinates;
        let (q, momentum) = (&state[..n], &state[n..2 * n]);
        if let Some((a, b)) = &p.maps {
            for i in 0..n {
                poll(&mut cancelled)?;
                self.candidate[i] = value(dot(&a[i * n..(i + 1) * n], q)?)?;
                self.candidate[n + i] = value(dot(&b[i * n..(i + 1) * n], momentum)?)?;
            }
        } else {
            self.candidate[..2 * n].copy_from_slice(&state[..2 * n]);
        }
        let mut maximum_relative_state_residual = 0.0_f64;
        for (j, &old) in p.previous_branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let before = det::sqrt(p.previous.branches[old].stiffness);
            let after = det::sqrt(p.target.branches[j].stiffness);
            let z = rescale(state[2 * n + old], after, before)?;
            self.candidate[2 * n + j] = z;
            maximum_relative_state_residual = maximum_relative_state_residual.max(
                relative_residual([(rescale(z, before, after)?, 1.0)], state[2 * n + old])?,
            );
        }
        let (next_q, next_p) = (&self.candidate[..n], &self.candidate[n..2 * n]);
        if let Some((a, b)) = &p.maps {
            for i in 0..n {
                poll(&mut cancelled)?;
                maximum_relative_state_residual = maximum_relative_state_residual
                    .max(relative_residual(
                        (0..n).map(|j| (b[j * n + i], next_q[j])),
                        q[i],
                    )?)
                    .max(relative_residual(
                        (0..n).map(|j| (a[j * n + i], next_p[j])),
                        momentum[i],
                    )?);
            }
        }
        limit(
            "canonical motion and memory residual",
            maximum_relative_state_residual,
            p.controls.maximum_relative_mapping_residual,
        )?;
        let elastic_work = quadratic(&p.stiffness_change, q, &mut cancelled)?;
        let kinetic_work = quadratic(&p.inverse_mass_change, momentum, &mut cancelled)?;
        let mut memory_work = Sum::default();
        let mut storage_change = quadratic_difference(
            &p.target.relaxed_stiffness,
            next_q,
            &p.previous.relaxed_stiffness,
            q,
            &mut cancelled,
        )?;
        storage_change
            .add_scaled(
                quadratic_difference(
                    &p.target.inverse_mass,
                    next_p,
                    &p.previous.inverse_mass,
                    momentum,
                    &mut cancelled,
                )?,
                1.0,
            )
            .ok_or_else(arithmetic)?;
        for (j, &old) in p.previous_branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let before = branch_extension(p.previous, old, state)?;
            let actual_after = branch_extension(p.target, j, &self.candidate)?;
            let (high, low) = before.components().ok_or_else(arithmetic)?;
            let old_root = det::sqrt(p.previous.branches[old].stiffness);
            let next_root = det::sqrt(p.target.branches[j].stiffness);
            let mut ideal_after = Sum::default();
            ideal_after
                .product(rescale(high, next_root, old_root)?, 1.0)
                .ok_or_else(arithmetic)?;
            // Carry the unrounded extension through its new material scale;
            // either unrepresentable word refuses instead of losing history.
            let tail = rescale(low, next_root, old_root)?;
            ideal_after.product(tail, 1.0).ok_or_else(arithmetic)?;
            memory_work
                .add_scaled(square_difference(ideal_after, before)?, 1.0)
                .ok_or_else(arithmetic)?;
            storage_change
                .add_scaled(square_difference(actual_after, before)?, 1.0)
                .ok_or_else(arithmetic)?;
        }
        let mut parameter_work = elastic_work;
        parameter_work
            .add_scaled(kinetic_work, 1.0)
            .ok_or_else(arithmetic)?;
        parameter_work
            .add_scaled(memory_work, 1.0)
            .ok_or_else(arithmetic)?;
        let parameter_work_j = value(parameter_work)?;
        limit(
            "external parameter work",
            parameter_work_j.abs(),
            p.controls.maximum_absolute_parameter_work_j,
        )?;
        let mut residual = storage_change;
        residual
            .add_scaled(parameter_work, -1.0)
            .ok_or_else(arithmetic)?;
        let energy_residual_j = value(residual)?;
        limit(
            "original storage and parameter work",
            energy_residual_j.abs(),
            p.controls.energy_tolerance_j,
        )?;
        let previous_stored_energy_j = checked(p.previous.system.hamiltonian(state))?;
        let target_stored_energy_j = checked(p.target.system.hamiltonian(&self.candidate))?;
        let report = MechanicalTransferReport {
            previous_stored_energy_j,
            target_stored_energy_j,
            kinetic_parameter_work_j: value(kinetic_work)?,
            elastic_parameter_work_j: value(elastic_work)?,
            memory_parameter_work_j: value(memory_work)?,
            parameter_work_j,
            storage_change_j: value(storage_change)?,
            energy_residual_j,
            maximum_relative_state_residual,
        };
        poll(&mut cancelled)?;
        target_state.copy_from_slice(&self.candidate);
        Ok(report)
    }
}

fn admit_controls(c: MechanicalTransferControls) -> Result<()> {
    if c.maximum_coordinates == 0
        || !c.maximum_condition_number.is_finite()
        || c.maximum_condition_number < 1.0
        || [
            c.maximum_relative_mapping_residual,
            c.maximum_relative_strain_residual,
        ]
        .iter()
        .any(|x| !x.is_finite() || !(0.0..=1e-6).contains(x))
        || [c.maximum_absolute_parameter_work_j, c.energy_tolerance_j]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(invalid(
            "finite bounded mechanical transfer controls required",
        ));
    }
    Ok(())
}
fn invalid(what: &'static str) -> MechanicalTransferError {
    MechanicalTransferError::Invalid(what)
}
fn arithmetic() -> MechanicalTransferError {
    invalid("unrepresentable mechanical transfer arithmetic")
}
fn checked(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(arithmetic())
    }
}
fn value(x: Sum) -> Result<f64> {
    x.value().ok_or_else(arithmetic)
}
fn dot(a: &[f64], b: &[f64]) -> Result<Sum> {
    Sum::dot(a, b).ok_or_else(arithmetic)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MechanicalTransferError::Cancelled)
    } else {
        Ok(())
    }
}
fn limit(quantity: &'static str, x: f64, maximum: f64) -> Result<()> {
    checked(x)?;
    if x > maximum {
        Err(MechanicalTransferError::Limit {
            quantity,
            value: x,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn relative_residual(terms: impl IntoIterator<Item = (f64, f64)>, expected: f64) -> Result<f64> {
    let mut residual = Sum::default();
    let mut scale = Sum::default();
    for (a, b) in terms {
        residual.product(a, b).ok_or_else(arithmetic)?;
        scale.product(a.abs(), b.abs()).ok_or_else(arithmetic)?;
    }
    residual.product(-expected, 1.0).ok_or_else(arithmetic)?;
    // max rather than a+b avoids overflowing a valid large history scale.
    let scale = value(scale)?.max(expected.abs());
    let r = value(residual)?.abs();
    if scale == 0.0 {
        Ok(0.0)
    } else {
        checked(r / scale)
    }
}

// Choose a finite evaluation order without normalizing coefficients or silently
// erasing a nonzero history. Equal roots retain every bit, including signed zero.
fn rescale(x: f64, numerator: f64, denominator: f64) -> Result<f64> {
    if x == 0.0 || numerator == denominator {
        return Ok(x);
    }
    for candidate in [
        x * (numerator / denominator),
        (x / denominator) * numerator,
        (x * numerator) / denominator,
    ] {
        if candidate.is_finite() && candidate != 0.0 {
            return Ok(candidate);
        }
    }
    Err(arithmetic())
}

fn pullback_difference(
    previous: &[f64],
    target: &[f64],
    map: Option<&[f64]>,
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let mut difference = vec![0.0; n * n];
    if let Some(a) = map {
        // Retain both words of target*A until A^T has been applied. A weak
        // physical coefficient must survive cancellation of large map terms.
        let mut intermediate = vec![Sum::default(); n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                for k in 0..n {
                    intermediate[i * n + j]
                        .product(target[i * n + k], a[k * n + j])
                        .ok_or_else(arithmetic)?;
                }
            }
        }
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                let mut entry = Sum::default();
                for k in 0..n {
                    entry
                        .add_scaled(intermediate[k * n + j], a[k * n + i])
                        .ok_or_else(arithmetic)?;
                }
                entry
                    .product(-previous[i * n + j], 1.0)
                    .ok_or_else(arithmetic)?;
                difference[i * n + j] = value(entry)?;
            }
        }
    } else {
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                difference[i * n + j] = checked(target[i * n + j] - previous[i * n + j])?;
            }
        }
    }
    Ok(difference)
}
fn quadratic(matrix: &[f64], x: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<Sum> {
    let n = x.len();
    let mut energy = Sum::default();
    for i in 0..n {
        poll(cancelled)?;
        energy
            .add_half_scaled(dot(&matrix[i * n..(i + 1) * n], x)?, x[i])
            .ok_or_else(arithmetic)?;
    }
    Ok(energy)
}
fn quadratic_difference(
    target: &[f64],
    after: &[f64],
    previous: &[f64],
    before: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Sum> {
    if after == before {
        // Fixed-coordinate parameter changes must not subtract complete
        // energies: even two words lose weak work below a large background's
        // product-rounding tail. Cancel identical source entries first, then
        // retain each changed coefficient's paired products before contraction.
        let n = before.len();
        let mut increment = Sum::default();
        for i in 0..n {
            poll(cancelled)?;
            let mut row = Sum::default();
            for j in 0..n {
                let index = i * n + j;
                if target[index] != previous[index] {
                    let mut change = Sum::default();
                    change
                        .product(target[index], before[j])
                        .ok_or_else(arithmetic)?;
                    change
                        .product(-previous[index], before[j])
                        .ok_or_else(arithmetic)?;
                    row.add_scaled(change, 1.0).ok_or_else(arithmetic)?;
                }
            }
            increment
                .add_half_scaled(row, before[i])
                .ok_or_else(arithmetic)?;
        }
        return Ok(increment);
    }
    let mut increment = quadratic(target, after, cancelled)?;
    increment
        .add_scaled(quadratic(previous, before, cancelled)?, -1.0)
        .ok_or_else(arithmetic)?;
    Ok(increment)
}
fn branch_extension(model: &LinearMechanicalSystem, branch: usize, state: &[f64]) -> Result<Sum> {
    let n = model.coordinates;
    let arm = &model.branches[branch];
    let root = det::sqrt(arm.stiffness);
    let mut extension = Sum::default();
    for i in 0..n {
        // These are the SAME represented couplings as the original storage.
        extension
            .product(root * arm.projection[i], state[i])
            .ok_or_else(arithmetic)?;
    }
    extension
        .product(-state[2 * n + branch], 1.0)
        .ok_or_else(arithmetic)?;
    Ok(extension)
}
fn square_difference(after: Sum, before: Sum) -> Result<Sum> {
    let mut difference = after;
    difference.add_scaled(before, -1.0).ok_or_else(arithmetic)?;
    let mut mean = Sum::default();
    mean.add_half_scaled(after, 1.0).ok_or_else(arithmetic)?;
    mean.add_half_scaled(before, 1.0).ok_or_else(arithmetic)?;
    let mut energy = Sum::default();
    energy
        .add_product(difference, mean)
        .ok_or_else(arithmetic)?;
    Ok(energy)
}
