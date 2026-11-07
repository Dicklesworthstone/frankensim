//! Native equilibrium support for the original sealed-cavity pressure owner.
//! Each gas contributes its own positive acoustic factor; no dense boundary
//! matrix or artificial coupling between independent inventories is introduced.

use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_math::dd::Dd;
use fs_solver::op::{CsrOp, LinearOp};

use super::{ClosedTetCavity, HyperTetError, TetCavityGasLoad, checkpoint, finite, invalid};

/// One stationary sealed-gas source bound to its original cavity geometry.
///
/// Every endpoint is reconstructed from `reference` through the same admitted
/// source isentrope. Other loads can share the native solid equilibrium solve;
/// this binding introduces neither a second gas model nor a separate solver.
#[derive(Debug, Clone, Copy)]
pub struct AdiabaticCavity<'a> {
    /// Prepared closed inner boundary of the complete reference solid.
    pub cavity: &'a ClosedTetCavity,
    /// Original gas inventory at the undeformed cavity volume.
    pub reference: &'a ClosedGasParcel,
    /// Fixed exterior absolute pressure [Pa], finite and nonnegative.
    pub exterior_pressure_pa: f64,
    /// Original source temperature bracket, entropy tolerance and work bound.
    pub controls: GasEntropySolve,
}

impl AdiabaticCavity<'_> {
    /// Resolve this same gas source on the admitted current solid geometry.
    pub fn evaluate(
        &self,
        displacement: &[f64],
        cx: &Cx<'_>,
    ) -> Result<TetCavityGasLoad, HyperTetError> {
        self.cavity.adiabatic_gas_load(
            displacement,
            self.reference,
            self.exterior_pressure_pa,
            self.controls,
            cx,
        )
    }

    /// Direct potential change between this binding's already admitted loads.
    /// The static owner retains these immutable endpoints from `evaluate`.
    /// Source calorics and cached compact cone geometry avoid subtracting two
    /// reference-relative potentials, re-solving entropy, or revisiting the
    /// full solid mesh. This uses the represented thermodynamic endpoints;
    /// their accuracy remains controlled by the declared entropy solve.
    pub(crate) fn potential_increment(
        &self,
        before: &TetCavityGasLoad,
        after: &TetCavityGasLoad,
        cx: &Cx<'_>,
    ) -> Result<Dd, HyperTetError> {
        checkpoint(cx)?;
        let stored = before
            .gas
            .stored_energy_change_to(&after.gas, || cx.checkpoint().is_err())
            .map_err(|error| match error {
                GasEnergyError::Cancelled => HyperTetError::Cancelled,
                other => invalid(&format!("cavity gas potential increment refused: {other}")),
            })?;
        let exterior_work = if self.exterior_pressure_pa == 0.0 {
            Dd::ZERO
        } else {
            Dd::from_f64(self.exterior_pressure_pa)
                * self
                    .cavity
                    .volume_increment(&before.positions_m, &after.positions_m, cx)?
        };
        let change = Dd::from_f64(stored) + exterior_work;
        finite(change.to_f64())?;
        checkpoint(cx)?;
        Ok(change)
    }
}

/// Free-coordinate factors of `sum_c w_c w_c^T` in original cavity order.
///
/// One cavity retains the original dense reduction, including zero entries and
/// its floating-point order. Multiple cavities retain only nonzero entries in
/// each free-coordinate factor. Their storage and Krylov application work scale
/// with the summed active boundary coordinates, not cavity count times all free
/// solid coordinates. Full source vectors and prescribed-coordinate work remain
/// with their original load owners; this is only a reduction for the linear solve.
#[derive(Debug)]
pub(crate) enum AcousticFactors {
    Empty,
    Single(Vec<f64>),
    Multiple {
        dimension: usize,
        offsets: Vec<usize>,
        rows: Vec<usize>,
        weights: Vec<f64>,
    },
}

impl AcousticFactors {
    pub const fn empty() -> Self {
        Self::Empty
    }

    /// Reduce original full-coordinate factors without merging gas inventories.
    /// The free chart must be unique and increasing, as in native elimination.
    /// Original source admission has already occurred; retained entries, common
    /// chart lengths and the continuation scale are checked again here.
    pub fn from_full_vectors<'a>(
        mut vectors: impl ExactSizeIterator<Item = &'a [f64]>,
        free: &[usize],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, HyperTetError> {
        checkpoint(cx)?;
        if !load.is_finite() || load < 0.0 {
            return Err(invalid(
                "cavity acoustic load scale must be finite and nonnegative",
            ));
        }
        let count = vectors.len();
        let Some(first) = vectors.next() else {
            return Ok(Self::Empty);
        };
        for (row, &dof) in free.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            if dof >= first.len() || (row > 0 && free[row - 1] >= dof) {
                return Err(invalid(
                    "cavity acoustic factors require the original ordered free chart",
                ));
            }
        }
        let amplitude = fs_math::det::sqrt(load);
        if count == 1 {
            let mut factor = Vec::with_capacity(free.len());
            for (row, &dof) in free.iter().enumerate() {
                if row % 256 == 0 {
                    checkpoint(cx)?;
                }
                factor.push(finite(amplitude * first[dof])?);
            }
            checkpoint(cx)?;
            return Ok(Self::Single(factor));
        }
        let mut offsets = Vec::new();
        let mut rows = Vec::new();
        let mut weights = Vec::new();
        offsets.push(0);
        for vector in std::iter::once(first).chain(vectors) {
            checkpoint(cx)?;
            if vector.len() != first.len() {
                return Err(invalid(
                    "cavity acoustic source vectors use different full charts",
                ));
            }
            for (row, &dof) in free.iter().enumerate() {
                if row % 256 == 0 {
                    checkpoint(cx)?;
                }
                let weight = finite(amplitude * vector[dof])?;
                if weight != 0.0 {
                    rows.push(row);
                    weights.push(weight);
                }
            }
            offsets.push(rows.len());
        }
        checkpoint(cx)?;
        Ok(Self::Multiple {
            dimension: free.len(),
            offsets,
            rows,
            weights,
        })
    }

    /// Add the complete acoustic diagonal to the caller's signed sparse
    /// diagonal before it forms an absolute-diagonal preconditioner.
    pub fn add_diagonal(&self, diagonal: &mut [f64], cx: &Cx<'_>) -> Result<(), HyperTetError> {
        checkpoint(cx)?;
        match self {
            Self::Empty => {}
            Self::Single(factor) => {
                if factor.len() != diagonal.len() {
                    return Err(invalid("cavity acoustic diagonal has the wrong free chart"));
                }
                for (row, (&weight, value)) in factor.iter().zip(diagonal).enumerate() {
                    if row % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    *value = finite(*value + finite(weight * weight)?)?;
                }
            }
            Self::Multiple {
                dimension,
                rows,
                weights,
                ..
            } => {
                if *dimension != diagonal.len() {
                    return Err(invalid("cavity acoustic diagonal has the wrong free chart"));
                }
                for (entry, (&row, &weight)) in rows.iter().zip(weights).enumerate() {
                    if entry % 256 == 0 {
                        checkpoint(cx)?;
                    }
                    diagonal[row] = finite(diagonal[row] + finite(weight * weight)?)?;
                }
            }
        }
        checkpoint(cx)?;
        Ok(())
    }

    // Allocation-free action, in fixed cavity and free-coordinate order. Keeping
    // projections separate avoids the cross terms from (sum w_c)(sum w_c)^T.
    fn add_action(&self, x: &[f64], y: &mut [f64]) {
        match self {
            Self::Empty => {}
            Self::Single(factor) => {
                let projection = factor
                    .iter()
                    .zip(x)
                    .fold(0.0, |sum, (&w, &v)| w.mul_add(v, sum));
                for (yi, &weight) in y.iter_mut().zip(factor) {
                    *yi = weight.mul_add(projection, *yi);
                }
            }
            Self::Multiple {
                offsets,
                rows,
                weights,
                ..
            } => {
                for bounds in offsets.windows(2) {
                    let range = bounds[0]..bounds[1];
                    let projection = rows[range.clone()]
                        .iter()
                        .zip(&weights[range.clone()])
                        .fold(0.0, |sum, (&row, &weight)| weight.mul_add(x[row], sum));
                    for (&row, &weight) in rows[range.clone()].iter().zip(&weights[range]) {
                        y[row] = weight.mul_add(projection, y[row]);
                    }
                }
            }
        }
    }
}

/// Original sparse operator plus every positive term in K_solid - sum J_gas.
/// Empty factors select the original general sparse operator exactly.
pub(crate) struct AcousticTangent<'a> {
    pub sparse: CsrOp,
    pub factors: &'a AcousticFactors,
}

impl LinearOp for AcousticTangent<'_> {
    fn n(&self) -> usize {
        self.sparse.n()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.sparse.apply(x, y);
        self.factors.add_action(x, y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.sparse.apply_transpose(x, y);
        self.factors.add_action(x, y);
    }
}

/// Translation-invariant bound on the dense acoustic deformation-force scale.
/// Three weighted component sums avoid an O(boundary DOFs squared) row walk.
pub(crate) fn acoustic_force_scale(
    factor: &[f64],
    displacement: &[f64],
    free: &[usize],
    cx: &Cx<'_>,
) -> Result<f64, HyperTetError> {
    if factor.len() != displacement.len() || displacement.len() < 3 {
        return Err(invalid(
            "cavity acoustic scale has the wrong displacement chart",
        ));
    }
    let Some(reference) = factor.iter().position(|&w| w != 0.0).map(|i| 3 * (i / 3)) else {
        return Ok(0.0);
    };
    let origin = &displacement[reference..reference + 3];
    let mut weight = [0.0; 3];
    let mut deviation = [0.0; 3];
    for (i, (&w, &u)) in factor.iter().zip(displacement).enumerate() {
        if i % 256 == 0 {
            checkpoint(cx)?;
        }
        let axis = i % 3;
        weight[axis] = finite(weight[axis] + w.abs())?;
        deviation[axis] = finite(w.abs().mul_add((u - origin[axis]).abs(), deviation[axis]))?;
    }
    let mut result = 0.0_f64;
    for (row, &dof) in free.iter().enumerate() {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        let node = 3 * (dof / 3);
        let mut bound = 0.0;
        for axis in 0..3 {
            bound = finite(
                bound
                    + deviation[axis]
                    + weight[axis] * (origin[axis] - displacement[node + axis]).abs(),
            )?;
        }
        result = finite(result.hypot(finite(factor[dof].abs() * bound)?))?;
    }
    checkpoint(cx)?;
    Ok(result)
}

/// Bound the full summed acoustic force by combining contributions on each
/// free row before taking its norm. Each cavity retains its own translation-
/// invariant anchor; distant chambers do not create a spurious lever arm.
/// Zero and one cavity retain the original normalization exactly.
pub(crate) fn acoustic_force_scale_many<'a>(
    mut vectors: impl ExactSizeIterator<Item = &'a [f64]>,
    displacement: &[f64],
    free: &[usize],
    cx: &Cx<'_>,
) -> Result<f64, HyperTetError> {
    checkpoint(cx)?;
    match vectors.len() {
        0 => return Ok(0.0),
        1 => {
            return acoustic_force_scale(
                vectors.next().expect("one original cavity factor"),
                displacement,
                free,
                cx,
            );
        }
        _ => {}
    }
    if displacement.len() < 3 || displacement.len() % 3 != 0 {
        return Err(invalid(
            "cavity acoustic scale has the wrong displacement chart",
        ));
    }
    for (row, &dof) in free.iter().enumerate() {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        if dof >= displacement.len() {
            return Err(invalid(
                "cavity acoustic scale has the wrong displacement chart",
            ));
        }
    }
    let mut bounds = vec![0.0; free.len()];
    for factor in vectors {
        checkpoint(cx)?;
        if factor.len() != displacement.len() {
            return Err(invalid(
                "cavity acoustic scale has the wrong displacement chart",
            ));
        }
        let mut reference = None;
        for (dof, &value) in factor.iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if value != 0.0 {
                reference = Some(3 * (dof / 3));
                break;
            }
        }
        let Some(reference) = reference else {
            continue;
        };
        let origin = &displacement[reference..reference + 3];
        let mut weight = [0.0; 3];
        let mut deviation = [0.0; 3];
        for (i, (&w, &u)) in factor.iter().zip(displacement).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            if w == 0.0 {
                continue;
            }
            let axis = i % 3;
            weight[axis] = finite(weight[axis] + w.abs())?;
            deviation[axis] = finite(w.abs().mul_add((u - origin[axis]).abs(), deviation[axis]))?;
        }
        for (row, (&dof, total)) in free.iter().zip(&mut bounds).enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            if factor[dof] == 0.0 {
                continue;
            }
            let node = 3 * (dof / 3);
            let mut bound = 0.0;
            for axis in 0..3 {
                bound = finite(
                    bound
                        + deviation[axis]
                        + weight[axis] * (origin[axis] - displacement[node + axis]).abs(),
                )?;
            }
            *total = finite(*total + finite(factor[dof].abs() * bound)?)?;
        }
    }
    let mut result = 0.0_f64;
    for (row, value) in bounds.into_iter().enumerate() {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        result = finite(result.hypot(value))?;
    }
    checkpoint(cx)?;
    Ok(result)
}

#[cfg(test)]
#[path = "../../../../tests/hyper3_cavity/gas_fixture.rs"]
mod gas_fixture;

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_sparse::Coo;

    fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let gate = CancelGate::new();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            f(
                &gate,
                &Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 31,
                        kernel_id: 59,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                ),
            );
        });
    }

    fn shell(cx: &Cx<'_>) -> (Vec<[f64; 3]>, ClosedTetCavity) {
        use crate::hyper3::pressure::ClosedTetPressureSurface;
        use crate::linear3::TetAssemblyBudget;
        use std::sync::Arc;
        let nodes = vec![
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
            [0.5, 0.5, 0.5],
            [1.5, 0.5, 0.5],
            [0.5, 1.5, 0.5],
            [0.5, 0.5, 1.5],
        ];
        let mut tets = Vec::new();
        for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            tets.extend([
                [a + 4, b + 4, c + 4, c],
                [a + 4, b + 4, b, c],
                [a + 4, a, b, c],
            ]);
        }
        let surface = Arc::new(
            ClosedTetPressureSurface::new(
                &nodes,
                &tets,
                TetAssemblyBudget::standard(),
                4 * tets.len(),
                cx,
            )
            .unwrap(),
        );
        let faces: Vec<_> = surface
            .faces()
            .iter()
            .enumerate()
            .filter_map(|(i, face)| face.iter().all(|&node| node >= 4).then_some(i))
            .collect();
        let cavity = ClosedTetCavity::new(surface, &faces, cx).unwrap();
        (nodes, cavity)
    }

    fn controls() -> GasEntropySolve {
        GasEntropySolve {
            minimum_temperature_k: 150.0,
            maximum_temperature_k: 1500.0,
            entropy_tolerance_j_kg_k: 1e-9,
            maximum_iterations: 128,
        }
    }

    #[test]
    fn g1_cavity_interval_uses_original_variable_calorics_and_matches_pressure_work() {
        with_cx(|_, cx| {
            let (nodes, cavity) = shell(cx);
            let mut direct = Vec::new();
            for datum in [0.0, 1e24] {
                let reference =
                    gas_fixture::gas(cavity.reference_volume_m3(), 100000.0, 0.001, datum);
                let source = AdiabaticCavity {
                    cavity: &cavity,
                    reference: &reference,
                    exterior_pressure_pa: 60000.0,
                    controls: controls(),
                };
                let at = |stretch: f64| {
                    source
                        .evaluate(
                            &nodes
                                .iter()
                                .flatten()
                                .map(|&x| (stretch - 1.0) * x)
                                .collect::<Vec<_>>(),
                            cx,
                        )
                        .unwrap()
                };
                let stretch = 0.94;
                let h = 2e-5;
                let (before, center, after) = (at(stretch - h), at(stretch), at(stretch + h));
                let change = source
                    .potential_increment(&before, &after, cx)
                    .unwrap()
                    .to_f64();
                let force_work = center
                    .force_n
                    .iter()
                    .zip(nodes.iter().flatten())
                    .map(|(force, direction)| force * direction)
                    .sum::<f64>();
                assert!((change / (2.0 * h) + force_work).abs() < 2e-7 * force_work.abs());
                let reverse = source
                    .potential_increment(&after, &before, cx)
                    .unwrap()
                    .to_f64();
                assert_eq!(reverse, -change);
                let t0 = before.gas.state().state().temperature;
                let t1 = after.gas.state().state().temperature;
                // Independent source integral, Cv/R = 2.5 + 0.001 T.
                let expected_u =
                    reference.mass_kg() * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
                        * (t1 - t0)
                        * (2.5 + 0.0005 * (t1 + t0));
                let expected_v = ((stretch + h).powi(3) - (stretch - h).powi(3)) / 6.0;
                let expected = expected_u + source.exterior_pressure_pa * expected_v;
                assert!((change - expected).abs() < 2e-9);
                if datum != 0.0 {
                    let rounded = reference.mass_kg()
                        * (after.gas.state().specific_internal_energy_j_kg()
                            - before.gas.state().specific_internal_energy_j_kg());
                    assert!((rounded - expected_u).abs() > 0.5 * expected_u.abs());
                }
                direct.push(change);
            }
            assert_eq!(direct[0], direct[1]);
        });
    }

    #[test]
    fn g1_cavity_interval_retains_work_below_a_large_reference_relative_potential() {
        with_cx(|gate, cx| {
            let (nodes, cavity) = shell(cx);
            let reference = gas_fixture::gas(cavity.reference_volume_m3(), 100000.0, 0.001, 1e24);
            let source = AdiabaticCavity {
                cavity: &cavity,
                reference: &reference,
                exterior_pressure_pa: 6.0 * 2.0_f64.powi(60),
                controls: controls(),
            };
            let stretch = 0.875;
            let before_u: Vec<_> = nodes
                .iter()
                .flatten()
                .map(|&x| (stretch - 1.0) * x)
                .collect();
            let mut after_u = before_u.clone();
            let epsilon = 2.0_f64.powi(-30);
            after_u[3 * 5 + 1] += epsilon;
            after_u[3 * 6] += epsilon;
            let before = source.evaluate(&before_u, cx).unwrap();
            let after = source.evaluate(&after_u, cx).unwrap();
            // Equal represented gas volumes retain the same admitted gas state.
            // The exterior load still performs -p_ext*s*epsilon^2/6 of work.
            assert_eq!(before.volume_m3, after.volume_m3);
            assert_eq!(before.gas, after.gas);
            assert_eq!(before.potential_change_j, after.potential_change_j);
            let change = source
                .potential_increment(&before, &after, cx)
                .unwrap()
                .to_f64();
            assert!((change + stretch).abs() < 3e-15);
            assert_eq!(
                source
                    .potential_increment(&after, &before, cx)
                    .unwrap()
                    .to_f64(),
                -change
            );
            gate.request();
            assert!(matches!(
                source.potential_increment(&before, &after, cx),
                Err(HyperTetError::Cancelled)
            ));
        });
    }

    #[test]
    fn cavity_acoustic_operator_preserves_the_full_action_and_transpose() {
        let mut sparse = Coo::new(3, 3);
        for (i, j, value) in [
            (0, 0, 3.0),
            (0, 1, -2.0),
            (1, 1, 4.0),
            (1, 2, 0.5),
            (2, 0, 1.0),
            (2, 2, -1.0),
        ] {
            sparse.push(i, j, value);
        }
        let matrix = sparse.assemble();
        let factor = [2.0, -3.0, 0.25];
        let factors = AcousticFactors::Single(factor.to_vec());
        let op = AcousticTangent {
            sparse: CsrOp::general(matrix.clone()),
            factors: &factors,
        };
        let x = [0.5, -2.0, 4.0];
        let y = [3.0, 0.25, -1.0];
        let mut ax = [0.0; 3];
        let mut aty = [0.0; 3];
        op.apply(&x, &mut ax);
        op.apply_transpose(&y, &mut aty);
        for i in 0..3 {
            let exact: f64 = (0..3)
                .map(|j| (matrix.get(i, j) + factor[i] * factor[j]) * x[j])
                .sum();
            assert!((ax[i] - exact).abs() < 1e-12);
        }
        let left: f64 = ax.iter().zip(y).map(|(a, b)| a * b).sum();
        let right: f64 = x.iter().zip(aty).map(|(a, b)| a * b).sum();
        assert!((left - right).abs() < 1e-12);
    }

    #[test]
    fn g0_multiple_cavity_factors_keep_independent_acoustic_energy_and_diagonal() {
        with_cx(|_, cx| {
            // Opposite factors still add stiffness. Merging them first would
            // erase both gas stores and invent cross-cavity terms.
            let source = [
                [2.0, 77.0, 0.0, -3.0, 0.0, 0.0, 0.5],
                [-2.0, 0.0, 0.0, 3.0, 11.0, 0.0, -0.5],
                [0.0, 0.0, 1.0, 0.0, 0.0, -4.0, 2.0],
            ];
            let free = [0, 2, 3, 6];
            let load = 0.25;
            let factors = AcousticFactors::from_full_vectors(
                source.iter().map(|v| v.as_slice()),
                &free,
                load,
                cx,
            )
            .unwrap();
            let AcousticFactors::Multiple { rows, weights, .. } = &factors else {
                panic!("several cavities must retain compact independent factors");
            };
            assert_eq!(weights.len(), 8);
            assert_eq!(rows.len(), weights.len());
            let matrix = [
                [3.0, -2.0, 0.0, 0.5],
                [0.0, 4.0, 1.0, 0.0],
                [2.0, 0.0, -8.0, 1.0],
                [0.0, 0.25, 0.0, -1.0],
            ];
            let mut sparse = Coo::new(4, 4);
            let mut expected = matrix;
            for i in 0..4 {
                for j in 0..4 {
                    sparse.push(i, j, matrix[i][j]);
                    for full in &source {
                        expected[i][j] += load * full[free[i]] * full[free[j]];
                    }
                }
            }
            let op = AcousticTangent {
                sparse: CsrOp::general(sparse.assemble()),
                factors: &factors,
            };
            let x = [0.5, -2.0, 4.0, -0.25];
            let mut action = [0.0; 4];
            let mut transpose = [0.0; 4];
            op.apply(&x, &mut action);
            op.apply_transpose(&x, &mut transpose);
            let mut diagonal: [f64; 4] = std::array::from_fn(|i| matrix[i][i]);
            factors.add_diagonal(&mut diagonal, cx).unwrap();
            for i in 0..4 {
                let exact: f64 = (0..4).map(|j| expected[i][j] * x[j]).sum();
                let exact_t: f64 = (0..4).map(|j| expected[j][i] * x[j]).sum();
                assert!((action[i] - exact).abs() < 2e-14);
                assert!((transpose[i] - exact_t).abs() < 2e-14);
                assert_eq!(diagonal[i], expected[i][i]);
            }
            assert!(
                diagonal[2] < 0.0,
                "the signed solid diagonal must survive until caller admission"
            );
        });
    }

    #[test]
    fn g0_empty_and_single_cavity_actions_preserve_the_original_bits() {
        with_cx(|_, cx| {
            let full = [0.0, -0.0, 0.25, -3.0, 1e-7, 7.0];
            let free = [0, 1, 2, 4, 5];
            let x = [-3.0, 0.5, 4.0, -1e8, 0.25];
            let mut entries = Coo::new(5, 5);
            for row in 0..5 {
                entries.push(row, row, row as f64 - 2.25);
                if row > 0 {
                    entries.push(row, row - 1, 0.3);
                }
            }
            let sparse = CsrOp::general(entries.assemble());
            let mut base = [0.0; 5];
            sparse.apply(&x, &mut base);
            let empty =
                AcousticFactors::from_full_vectors(std::iter::empty(), &free, 1.0, cx).unwrap();
            let op = AcousticTangent {
                sparse: CsrOp::general(sparse.matrix().clone()),
                factors: &empty,
            };
            let mut empty_action = [0.0; 5];
            op.apply(&x, &mut empty_action);
            assert_eq!(base.map(f64::to_bits), empty_action.map(f64::to_bits));
            for load in [0.0, 0.3, 1.0] {
                let factors = AcousticFactors::from_full_vectors(
                    std::iter::once(full.as_slice()),
                    &free,
                    load,
                    cx,
                )
                .unwrap();
                let original: Vec<_> = free
                    .iter()
                    .map(|&dof| fs_math::det::sqrt(load) * full[dof])
                    .collect();
                let projection = original
                    .iter()
                    .zip(x)
                    .fold(0.0, |sum, (&w, v)| w.mul_add(v, sum));
                let expected: [f64; 5] =
                    std::array::from_fn(|i| original[i].mul_add(projection, base[i]));
                let op = AcousticTangent {
                    sparse: CsrOp::general(sparse.matrix().clone()),
                    factors: &factors,
                };
                let mut action = [0.0; 5];
                op.apply(&x, &mut action);
                assert_eq!(action.map(f64::to_bits), expected.map(f64::to_bits));
            }
        });
    }

    #[test]
    fn g0_cavity_factor_reduction_checks_charts_scaling_and_cancellation() {
        with_cx(|gate, cx| {
            let full = [1.0, 2.0, 3.0];
            for free in [&[1, 1][..], &[2, 0][..], &[3][..]] {
                assert!(
                    AcousticFactors::from_full_vectors(
                        std::iter::once(full.as_slice()),
                        free,
                        1.0,
                        cx,
                    )
                    .is_err()
                );
            }
            for load in [-1.0, f64::INFINITY, f64::NAN] {
                assert!(
                    AcousticFactors::from_full_vectors(
                        std::iter::once(full.as_slice()),
                        &[0, 2],
                        load,
                        cx,
                    )
                    .is_err()
                );
            }
            assert!(
                AcousticFactors::from_full_vectors(
                    [full.as_slice(), &full[..2]].into_iter(),
                    &[0, 2],
                    1.0,
                    cx,
                )
                .is_err()
            );
            assert!(AcousticFactors::from_full_vectors(
                std::iter::once(&[f64::NAN][..]), &[0], 1.0, cx,
            ).is_err());
            let factor = AcousticFactors::from_full_vectors(
                std::iter::once(full.as_slice()),
                &[0, 2],
                1.0,
                cx,
            )
            .unwrap();
            assert!(factor.add_diagonal(&mut [0.0], cx).is_err());
            let huge =
                AcousticFactors::from_full_vectors(std::iter::once(&[f64::MAX][..]), &[0], 1.0, cx)
                    .unwrap();
            assert!(huge.add_diagonal(&mut [0.0], cx).is_err());
            // Cancellation while obtaining the second source discards the
            // first source's partial reduction without returning an operator.
            let mut count = 0;
            let cancelled = [full.as_slice(), full.as_slice()].into_iter().inspect(|_| {
                count += 1;
                if count == 2 {
                    gate.request();
                }
            });
            assert!(matches!(
                AcousticFactors::from_full_vectors(cancelled, &[0, 2], 1.0, cx),
                Err(HyperTetError::Cancelled)
            ));
            assert!(matches!(
                factor.add_diagonal(&mut [0.0; 2], cx),
                Err(HyperTetError::Cancelled)
            ));
        });
    }

    #[test]
    fn g0_multiple_cavity_force_scale_sums_rows_and_retains_translation_invariance() {
        with_cx(|gate, cx| {
            let first = [1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
            let second = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, -2.0, 0.0];
            let u = [1.0, 2.0, 3.0, 3.0, 2.0, 3.0, 7.0, 1.0, 0.0, 7.0, 4.0, 0.0];
            let free: Vec<_> = (0..u.len()).collect();
            let a = acoustic_force_scale(&first, &u, &free, cx).unwrap();
            let b = acoustic_force_scale(&second, &u, &free, cx).unwrap();
            let both = acoustic_force_scale_many(
                [first.as_slice(), second.as_slice()].into_iter(),
                &u,
                &free,
                cx,
            )
            .unwrap();
            assert!((both - a.hypot(b)).abs() < 2e-14);
            let repeated = acoustic_force_scale_many(
                [first.as_slice(), first.as_slice()].into_iter(),
                &u,
                &free,
                cx,
            )
            .unwrap();
            assert_eq!(repeated, 2.0 * a);
            let shifted = u.map(|value| value + 2.0_f64.powi(30));
            assert_eq!(
                both,
                acoustic_force_scale_many(
                    [first.as_slice(), second.as_slice()].into_iter(),
                    &shifted,
                    &free,
                    cx,
                )
                .unwrap()
            );
            assert_eq!(
                a.to_bits(),
                acoustic_force_scale_many(std::iter::once(first.as_slice()), &u, &free, cx,)
                    .unwrap()
                    .to_bits()
            );
            assert_eq!(
                acoustic_force_scale_many(std::iter::empty(), &u, &free, cx).unwrap(),
                0.0
            );
            assert!(
                acoustic_force_scale_many(
                    [first.as_slice(), &second[..9]].into_iter(),
                    &u,
                    &free,
                    cx,
                )
                .is_err()
            );
            gate.request();
            assert!(matches!(
                acoustic_force_scale_many(
                    [first.as_slice(), second.as_slice()].into_iter(),
                    &u,
                    &free,
                    cx,
                ),
                Err(HyperTetError::Cancelled)
            ));
        });
    }
}
