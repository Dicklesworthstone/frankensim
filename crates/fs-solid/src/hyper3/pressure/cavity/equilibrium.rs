//! Native equilibrium support for the original sealed-cavity pressure owner.
//! The acoustic stiffness remains one vector, never a dense boundary matrix.

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
    pub fn evaluate(&self, displacement: &[f64], cx: &Cx<'_>) -> Result<TetCavityGasLoad, HyperTetError> {
        self.cavity.adiabatic_gas_load(
            displacement, self.reference, self.exterior_pressure_pa, self.controls, cx,
        )
    }

    /// Direct potential change between this binding's already admitted loads.
    /// The static owner retains these immutable endpoints from `evaluate`.
    /// Source calorics and cached compact cone geometry avoid subtracting two
    /// reference-relative potentials, re-solving entropy, or revisiting the
    /// full solid mesh. This uses the represented thermodynamic endpoints;
    /// their accuracy remains controlled by the declared entropy solve.
    pub(crate) fn potential_increment(
        &self, before: &TetCavityGasLoad, after: &TetCavityGasLoad, cx: &Cx<'_>,
    ) -> Result<Dd, HyperTetError> {
        checkpoint(cx)?;
        let stored = before.gas.stored_energy_change_to(&after.gas, || cx.checkpoint().is_err())
            .map_err(|error| match error {
                GasEnergyError::Cancelled => HyperTetError::Cancelled,
                other => invalid(&format!("cavity gas potential increment refused: {other}")),
            })?;
        let exterior_work = if self.exterior_pressure_pa == 0.0 { Dd::ZERO } else {
            Dd::from_f64(self.exterior_pressure_pa)
                * self.cavity.volume_increment(&before.positions_m, &after.positions_m, cx)?
        };
        let change = Dd::from_f64(stored) + exterior_work;
        finite(change.to_f64())?;
        checkpoint(cx)?;
        Ok(change)
    }
}

/// Original sparse operator plus the positive acoustic term in K_solid - J_gas.
/// An empty factor selects the original general sparse operator exactly.
pub(crate) struct AcousticTangent<'a> {
    pub sparse: CsrOp,
    pub factor: &'a [f64],
}

impl AcousticTangent<'_> {
    fn add_acoustic(&self, x: &[f64], y: &mut [f64]) {
        if self.factor.is_empty() { return; }
        let projection = self.factor.iter().zip(x).fold(0.0, |s, (&w, &v)| w.mul_add(v, s));
        for (yi, &weight) in y.iter_mut().zip(self.factor) {
            *yi = weight.mul_add(projection, *yi);
        }
    }
}

impl LinearOp for AcousticTangent<'_> {
    fn n(&self) -> usize { self.sparse.n() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.sparse.apply(x, y);
        self.add_acoustic(x, y);
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.sparse.apply_transpose(x, y);
        self.add_acoustic(x, y);
    }
}

/// Translation-invariant bound on the dense acoustic deformation-force scale.
/// Three weighted component sums avoid an O(boundary DOFs squared) row walk.
pub(crate) fn acoustic_force_scale(
    factor: &[f64], displacement: &[f64], free: &[usize], cx: &Cx<'_>,
) -> Result<f64, HyperTetError> {
    if factor.len() != displacement.len() || displacement.len() < 3 {
        return Err(invalid("cavity acoustic scale has the wrong displacement chart"));
    }
    let Some(reference) = factor.iter().position(|&w| w != 0.0).map(|i| 3 * (i / 3)) else {
        return Ok(0.0);
    };
    let origin = &displacement[reference..reference + 3];
    let mut weight = [0.0; 3];
    let mut deviation = [0.0; 3];
    for (i, (&w, &u)) in factor.iter().zip(displacement).enumerate() {
        if i % 256 == 0 { checkpoint(cx)?; }
        let axis = i % 3;
        weight[axis] = finite(weight[axis] + w.abs())?;
        deviation[axis] = finite(w.abs().mul_add((u - origin[axis]).abs(), deviation[axis]))?;
    }
    let mut result = 0.0_f64;
    for (row, &dof) in free.iter().enumerate() {
        if row % 256 == 0 { checkpoint(cx)?; }
        let node = 3 * (dof / 3);
        let mut bound = 0.0;
        for axis in 0..3 {
            bound = finite(bound + deviation[axis]
                + weight[axis] * (origin[axis] - displacement[node + axis]).abs())?;
        }
        result = finite(result.hypot(finite(factor[dof].abs() * bound)?))?;
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
            f(&gate, &Cx::new(&gate, arena,
                StreamKey { seed: 31, kernel_id: 59, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic));
        });
    }

    fn shell(cx: &Cx<'_>) -> (Vec<[f64; 3]>, ClosedTetCavity) {
        use crate::hyper3::pressure::ClosedTetPressureSurface;
        use crate::linear3::TetAssemblyBudget;
        use std::sync::Arc;
        let nodes = vec![
            [0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0],
            [0.5, 0.5, 0.5], [1.5, 0.5, 0.5], [0.5, 1.5, 0.5], [0.5, 0.5, 1.5],
        ];
        let mut tets = Vec::new();
        for [a, b, c] in [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]] {
            tets.extend([[a + 4, b + 4, c + 4, c], [a + 4, b + 4, b, c], [a + 4, a, b, c]]);
        }
        let surface = Arc::new(ClosedTetPressureSurface::new(
            &nodes, &tets, TetAssemblyBudget::standard(), 4 * tets.len(), cx,
        ).unwrap());
        let faces: Vec<_> = surface.faces().iter().enumerate()
            .filter_map(|(i, face)| face.iter().all(|&node| node >= 4).then_some(i)).collect();
        let cavity = ClosedTetCavity::new(surface, &faces, cx).unwrap();
        (nodes, cavity)
    }

    fn controls() -> GasEntropySolve {
        GasEntropySolve { minimum_temperature_k: 150.0, maximum_temperature_k: 1500.0,
            entropy_tolerance_j_kg_k: 1e-9, maximum_iterations: 128 }
    }

    #[test]
    fn g1_cavity_interval_uses_original_variable_calorics_and_matches_pressure_work() {
        with_cx(|_, cx| {
            let (nodes, cavity) = shell(cx);
            let mut direct = Vec::new();
            for datum in [0.0, 1e24] {
                let reference = gas_fixture::gas(cavity.reference_volume_m3(), 100000.0, 0.001, datum);
                let source = AdiabaticCavity { cavity: &cavity, reference: &reference,
                    exterior_pressure_pa: 60000.0, controls: controls() };
                let at = |stretch: f64| source.evaluate(
                    &nodes.iter().flatten().map(|&x| (stretch - 1.0) * x).collect::<Vec<_>>(), cx,
                ).unwrap();
                let stretch = 0.94;
                let h = 2e-5;
                let (before, center, after) = (at(stretch - h), at(stretch), at(stretch + h));
                let change = source.potential_increment(&before, &after, cx).unwrap().to_f64();
                let force_work = center.force_n.iter().zip(nodes.iter().flatten())
                    .map(|(force, direction)| force * direction).sum::<f64>();
                assert!((change / (2.0 * h) + force_work).abs() < 2e-7 * force_work.abs());
                let reverse = source.potential_increment(&after, &before, cx).unwrap().to_f64();
                assert_eq!(reverse, -change);
                let t0 = before.gas.state().state().temperature;
                let t1 = after.gas.state().state().temperature;
                // Independent source integral, Cv/R = 2.5 + 0.001 T.
                let expected_u = reference.mass_kg()
                    * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028
                    * (t1 - t0) * (2.5 + 0.0005 * (t1 + t0));
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
            let source = AdiabaticCavity { cavity: &cavity, reference: &reference,
                exterior_pressure_pa: 6.0 * 2.0_f64.powi(60), controls: controls() };
            let stretch = 0.875;
            let before_u: Vec<_> = nodes.iter().flatten().map(|&x| (stretch - 1.0) * x).collect();
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
            let change = source.potential_increment(&before, &after, cx).unwrap().to_f64();
            assert!((change + stretch).abs() < 3e-15);
            assert_eq!(source.potential_increment(&after, &before, cx).unwrap().to_f64(), -change);
            gate.request();
            assert!(matches!(source.potential_increment(&before, &after, cx), Err(HyperTetError::Cancelled)));
        });
    }

    #[test]
    fn cavity_acoustic_operator_preserves_the_full_action_and_transpose() {
        let mut sparse = Coo::new(3, 3);
        for (i, j, value) in [(0, 0, 3.0), (0, 1, -2.0), (1, 1, 4.0),
            (1, 2, 0.5), (2, 0, 1.0), (2, 2, -1.0)] {
            sparse.push(i, j, value);
        }
        let matrix = sparse.assemble();
        let factor = [2.0, -3.0, 0.25];
        let op = AcousticTangent { sparse: CsrOp::general(matrix.clone()), factor: &factor };
        let x = [0.5, -2.0, 4.0];
        let y = [3.0, 0.25, -1.0];
        let mut ax = [0.0; 3];
        let mut aty = [0.0; 3];
        op.apply(&x, &mut ax);
        op.apply_transpose(&y, &mut aty);
        for i in 0..3 {
            let exact: f64 = (0..3).map(|j| (matrix.get(i, j) + factor[i]*factor[j])*x[j]).sum();
            assert!((ax[i] - exact).abs() < 1e-12);
        }
        let left: f64 = ax.iter().zip(y).map(|(a, b)| a*b).sum();
        let right: f64 = x.iter().zip(aty).map(|(a, b)| a*b).sum();
        assert!((left - right).abs() < 1e-12);
    }
}
