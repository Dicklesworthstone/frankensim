//! Signed reference-geometry initial-stress stiffness on the original cut rules.
//!
//! G[(a,i),(b,j)] = delta_ij integral grad(N_a)^T sigma_0 grad(N_b).
//! Selected Nitsche supports add the symmetric geometric flux terms
//! -integral [N_a n^T sigma_0 grad(N_b) + N_b n^T sigma_0 grad(N_a)].
//! The original elastic penalty and ghost stabilization remain frozen in K_e;
//! neither creates an extra geometric penalty, artificial stress or inertia.
//!
//! This is the same conservative initial-stress approximation as body-fitted
//! linear elasticity, not a finite-strain dP/dF tangent or follower-load model.
//! A supplied state is physically admitted, but is not certified equilibrated.
//! G and K_e+G can be indefinite. No coercivity or stability claim is inferred.

use std::collections::{BTreeMap, btree_map::Entry};
use std::ops::ControlFlow;
use std::sync::Arc;

use fs_solver::op::LinearOp;

use super::{PreparedEigenstrain3, physical_response, poll};
use crate::elastic3::{Cell3, CutElasticity3, ElasticityError3, q1};

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

/// Immutable signed geometric stiffness at one original admitted thermal state.
///
/// The cached coefficients already contain the exact preparation-time stiffness
/// scales. No second scale or occupancy multiplier belongs in its action.
/// Rebuilding geometry/material/supports or changing scales invalidates its
/// combination with the current elastic operator; call `validate_prestress`.
/// An updated displacement or free-strain field requires a fresh preparation.
///
/// Scalar upper-triangular node pairs are shared by all three vector components.
/// Strong clamps are omitted and adaptive hanging constraints are condensed
/// once, so later actions neither reconstruct physical fields nor query sources.
pub struct PreparedPrestress3 {
    binding: Arc<()>,
    adaptive: bool,
    n: usize,
    fixed: Vec<bool>,
    scale_bits: Vec<u64>,
    pairs: Vec<[usize; 2]>,
    values: Vec<f64>,
}

impl core::fmt::Debug for PreparedPrestress3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("PreparedPrestress3")
            .field("displacement_dofs", &self.n)
            .field("scalar_entries", &self.values.len())
            .field("adaptive", &self.adaptive)
            .finish_non_exhaustive()
    }
}

impl PreparedPrestress3 {
    /// Original independent displacement dimension, including zero clamp rows.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Stored scalar upper-triangular entries before xyz replication.
    #[must_use]
    pub fn scalar_entries(&self) -> usize {
        self.values.len()
    }

    /// Original homogeneous node clamps; all three components are constrained.
    #[must_use]
    pub fn fixed(&self) -> &[bool] {
        &self.fixed
    }

    /// Read-only scalar upper-triangular pairs in deterministic cell/key order.
    ///
    /// Replicate each value over xyz and mirror off-diagonal entries. A pair can
    /// recur across original cells; an assembler must sum those contributions.
    /// This allocation-free view permits bounded sparse free-coordinate assembly
    /// without probing the operator on coordinate basis vectors.
    pub fn scalar_pairs(&self) -> impl ExactSizeIterator<Item = ([usize; 2], f64)> + '_ {
        self.pairs.iter().copied().zip(self.values.iter().copied())
    }

    fn admit_field(
        &self,
        field: &[f64],
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        poll(checkpoint)?;
        if field.len() != self.n {
            return Err(invalid("prestress field dimension mismatch"));
        }
        for chunk in field.chunks(192) {
            poll(checkpoint)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(invalid("nonfinite prestress field"));
            }
        }
        Ok(())
    }

    /// Apply the signed frozen geometric stiffness without allocating scratch.
    ///
    /// Directions are algebraic vectors; their arbitrary normalization does not
    /// undergo a physical strain-amplitude check. Input shapes and finite values
    /// are admitted before writing. Arithmetic failure or cancellation during
    /// accumulation can leave partial output, which the caller must discard.
    pub fn apply_into(
        &self,
        input: &[f64],
        output: &mut [f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<(), ElasticityError3> {
        self.admit_field(input, &mut checkpoint)?;
        if output.len() != self.n {
            return Err(invalid("prestress output dimension mismatch"));
        }
        output.fill(0.0);
        for (entry, ([a, b], value)) in self.scalar_pairs().enumerate() {
            if entry % 64 == 0 {
                poll(&mut checkpoint)?;
            }
            for component in 0..3 {
                output[3 * a + component] += value * input[3 * b + component];
                if a != b {
                    output[3 * b + component] += value * input[3 * a + component];
                }
            }
        }
        for chunk in output.chunks(192) {
            poll(&mut checkpoint)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(invalid("prestress action overflow"));
            }
        }
        poll(&mut checkpoint)
    }

    /// Signed diagonal of G. Negative entries are retained, not regularized.
    pub fn diagonal(
        &self,
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Vec<f64>, ElasticityError3> {
        poll(&mut checkpoint)?;
        let mut diagonal = vec![0.0; self.n];
        for (entry, ([a, b], value)) in self.scalar_pairs().enumerate() {
            if entry % 64 == 0 {
                poll(&mut checkpoint)?;
            }
            if a == b {
                for component in 0..3 {
                    diagonal[3 * a + component] += value;
                }
            }
        }
        for chunk in diagonal.chunks(192) {
            poll(&mut checkpoint)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(invalid("prestress diagonal overflow"));
            }
        }
        poll(&mut checkpoint)?;
        Ok(diagonal)
    }

    /// Exact cached bilinear work `left^T G right`, without a global action.
    /// No positivity or modal normalization is assumed for either direction.
    pub fn bilinear_work(
        &self,
        left: &[f64],
        right: &[f64],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<f64, ElasticityError3> {
        self.admit_field(left, &mut checkpoint)?;
        self.admit_field(right, &mut checkpoint)?;
        let mut work = 0.0;
        for (entry, ([a, b], value)) in self.scalar_pairs().enumerate() {
            if entry % 64 == 0 {
                poll(&mut checkpoint)?;
            }
            for component in 0..3 {
                work += value * left[3 * a + component] * right[3 * b + component];
                if a != b {
                    work += value * left[3 * b + component] * right[3 * a + component];
                }
            }
            if !work.is_finite() {
                return Err(invalid("prestress bilinear work overflow"));
            }
        }
        poll(&mut checkpoint)?;
        Ok(work)
    }

    pub(in crate::elastic3) fn validate_for(
        &self,
        raw: &CutElasticity3,
        adaptive: bool,
        n: usize,
    ) -> Result<(), ElasticityError3> {
        if !Arc::ptr_eq(&self.binding, &raw.material_body_binding)
            || self.adaptive != adaptive
            || self.n != n
            || self.scale_bits.len() != raw.scales.len()
            || self
                .scale_bits
                .iter()
                .zip(&raw.scales)
                .any(|(&bits, scale)| bits != scale.to_bits())
        {
            return Err(invalid(
                "prestress belongs to another owner, displacement space or stiffness state",
            ));
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::elastic3) fn prepare(
        raw: &CutElasticity3,
        rows: Option<&[Vec<(usize, f64)>]>,
        fixed: &[bool],
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        max_entries: usize,
        mut local_field: impl FnMut(&Cell3, &[f64]) -> [f64; 24],
        mut checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        poll(&mut checkpoint)?;
        raw.validate_eigenstrain(prepared)?;
        let n = fixed
            .len()
            .checked_mul(3)
            .ok_or(invalid("prestress displacement size overflow"))?;
        if u.len() != n
            || rows.is_some_and(|rows| rows.len() != raw.nodes.len())
            || (rows.is_none() && fixed.len() != raw.nodes.len())
            || prepared.cells.len() != raw.cells.len()
            || !maximum_kinematic_strain.is_finite()
            || maximum_kinematic_strain <= 0.0
            || maximum_kinematic_strain >= 1.0
        {
            return Err(invalid(
                "invalid prestress field, constraints or kinematic bound",
            ));
        }
        for chunk in u.chunks(192) {
            poll(&mut checkpoint)?;
            if chunk.iter().any(|value| !value.is_finite()) {
                return Err(invalid("nonfinite prestress displacement"));
            }
        }
        let mut points = 0usize;
        for (id, (cell, free)) in raw.cells.iter().zip(&prepared.cells).enumerate() {
            poll(&mut checkpoint)?;
            let selected = raw
                .embedded
                .as_ref()
                .map_or(0, |data| data.selected_points(id).len());
            if cell.rules.bulk().len() != free.bulk_strains.len()
                || selected != free.surface_strains.len()
            {
                return Err(invalid("prestress strain/quadrature alignment mismatch"));
            }
            points = points
                .checked_add(cell.rules.bulk().len())
                .and_then(|count| count.checked_add(selected))
                .ok_or(invalid("prestress point count overflow"))?;
            if points > max_points {
                return Err(invalid("prestress point allowance exhausted"));
            }
        }
        let allocation = || invalid("prestress allocation failed");
        let mut result = Self {
            binding: Arc::clone(&raw.material_body_binding),
            adaptive: rows.is_some(),
            n,
            fixed: Vec::new(),
            scale_bits: Vec::new(),
            pairs: Vec::new(),
            values: Vec::new(),
        };
        result
            .fixed
            .try_reserve(fixed.len())
            .map_err(|_| allocation())?;
        result.fixed.extend_from_slice(fixed);
        result
            .scale_bits
            .try_reserve(raw.scales.len())
            .map_err(|_| allocation())?;
        result
            .scale_bits
            .extend(raw.scales.iter().map(|scale| scale.to_bits()));
        for (id, (cell, free)) in raw.cells.iter().zip(&prepared.cells).enumerate() {
            poll(&mut checkpoint)?;
            let displacement = local_field(cell, u);
            let mut local = [[0.0; 8]; 8];
            for (&(position, weight), &strain) in cell.rules.bulk().iter().zip(&free.bulk_strains) {
                poll(&mut checkpoint)?;
                let stress = physical_response(
                    cell,
                    &displacement,
                    position,
                    strain,
                    free.mean_trace,
                    raw.scales[id],
                    maximum_kinematic_strain,
                )?
                .stress;
                let (_, gradients) = q1(cell.bounds, position);
                let flux = gradients.map(|gradient| stress_gradient(stress, gradient));
                for a in 0..8 {
                    for b in a..8 {
                        local[a][b] += weight * dot(gradients[a], flux[b]);
                        if !local[a][b].is_finite() {
                            return Err(invalid("prestress bulk integration overflow"));
                        }
                    }
                }
            }
            if let Some(support) = &raw.embedded {
                let selected = support.selected_points(id);
                if !selected.is_empty() {
                    let surface = cell
                        .rules
                        .surface()
                        .ok_or(invalid("prestress support quadrature is absent"))?;
                    for (&index, &strain) in selected.iter().zip(&free.surface_strains) {
                        poll(&mut checkpoint)?;
                        let point = surface
                            .points()
                            .get(index)
                            .ok_or(invalid("prestress support point is absent"))?;
                        let stress = physical_response(
                            cell,
                            &displacement,
                            point.position,
                            strain,
                            free.mean_trace,
                            raw.scales[id],
                            maximum_kinematic_strain,
                        )?
                        .stress;
                        let (shape, gradients) = q1(cell.bounds, point.position);
                        let flux = gradients
                            .map(|gradient| dot(point.normal, stress_gradient(stress, gradient)));
                        for a in 0..8 {
                            for b in a..8 {
                                local[a][b] -=
                                    point.weight * (shape[a] * flux[b] + shape[b] * flux[a]);
                                if !local[a][b].is_finite() {
                                    return Err(invalid("prestress boundary integration overflow"));
                                }
                            }
                        }
                    }
                }
            }
            for a in 0..8 {
                for b in a..8 {
                    local[b][a] = local[a][b];
                }
            }
            let mut block = BTreeMap::<[usize; 2], f64>::new();
            for a in 0..8 {
                let identity_a = [(cell.nodes[a], 1.0)];
                let row_a =
                    rows.map_or(identity_a.as_slice(), |rows| rows[cell.nodes[a]].as_slice());
                for b in 0..8 {
                    poll(&mut checkpoint)?;
                    let identity_b = [(cell.nodes[b], 1.0)];
                    let row_b =
                        rows.map_or(identity_b.as_slice(), |rows| rows[cell.nodes[b]].as_slice());
                    for &(i, wi) in row_a {
                        if i >= fixed.len() || !wi.is_finite() {
                            return Err(invalid("invalid prestress constraint row"));
                        }
                        if fixed[i] || wi == 0.0 {
                            continue;
                        }
                        for &(j, wj) in row_b {
                            poll(&mut checkpoint)?;
                            if j >= fixed.len() || !wj.is_finite() {
                                return Err(invalid("invalid prestress constraint row"));
                            }
                            if fixed[j] || wj == 0.0 || i > j {
                                continue;
                            }
                            let value = wi * local[a][b] * wj;
                            if !value.is_finite() {
                                return Err(invalid("constrained prestress overflow"));
                            }
                            if value == 0.0 {
                                continue;
                            }
                            let potential = result
                                .values
                                .len()
                                .checked_add(block.len())
                                .and_then(|count| count.checked_add(1))
                                .ok_or(invalid("prestress entry count overflow"))?;
                            let target = match block.entry([i, j]) {
                                Entry::Occupied(entry) => entry.into_mut(),
                                Entry::Vacant(entry) => {
                                    if potential > max_entries {
                                        return Err(invalid("prestress entry allowance exhausted"));
                                    }
                                    entry.insert(0.0)
                                }
                            };
                            *target += value;
                            if !target.is_finite() {
                                return Err(invalid("constrained prestress overflow"));
                            }
                        }
                    }
                }
            }
            result
                .pairs
                .try_reserve(block.len())
                .map_err(|_| allocation())?;
            result
                .values
                .try_reserve(block.len())
                .map_err(|_| allocation())?;
            for (pair, value) in block {
                poll(&mut checkpoint)?;
                if value != 0.0 {
                    result.pairs.push(pair);
                    result.values.push(value);
                }
            }
        }
        poll(&mut checkpoint)?;
        Ok(result)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

fn stress_gradient(stress: [f64; 6], gradient: [f64; 3]) -> [f64; 3] {
    [
        stress[0] * gradient[0] + stress[3] * gradient[1] + stress[5] * gradient[2],
        stress[3] * gradient[0] + stress[1] * gradient[1] + stress[4] * gradient[2],
        stress[5] * gradient[0] + stress[4] * gradient[1] + stress[2] * gradient[2],
    ]
}

impl CutElasticity3 {
    /// Freeze signed initial-stress stiffness at the supplied thermal state.
    ///
    /// Reuses original tensors, bulk rules, B-bar stress and selected Nitsche
    /// support points. `max_points` bounds bulk PLUS selected surface recovery;
    /// `max_entries` bounds retained scalar pairs and the current cell map before
    /// xyz replication. Physical domains are checked at those retained points.
    ///
    /// The original elastic Nitsche penalty and ghost coefficients remain
    /// frozen in K_e. No geometric penalty or ghost stress is added to G.
    /// The caller separately checks equilibrium before interpreting vibration
    /// or buckling about this state. No nonlinear stability claim is made.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_eigenstrain_prestress(
        &self,
        prepared: &PreparedEigenstrain3,
        u: &[f64],
        maximum_kinematic_strain: f64,
        max_points: usize,
        max_entries: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<PreparedPrestress3, ElasticityError3> {
        PreparedPrestress3::prepare(
            self,
            None,
            &self.fixed,
            prepared,
            u,
            maximum_kinematic_strain,
            max_points,
            max_entries,
            |cell, field| {
                std::array::from_fn(|i| {
                    let node = cell.nodes[i / 3];
                    if self.fixed[node] {
                        0.0
                    } else {
                        field[3 * node + i % 3]
                    }
                })
            },
            checkpoint,
        )
    }

    /// Require the same geometry/material/support space AND exact frozen scales.
    pub fn validate_prestress(
        &self,
        prestress: &PreparedPrestress3,
    ) -> Result<(), ElasticityError3> {
        prestress.validate_for(self, false, self.n())
    }
}
