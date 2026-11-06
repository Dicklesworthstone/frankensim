//! Native equilibrium support for the original sealed-cavity pressure owner.
//! The acoustic stiffness remains one vector, never a dense boundary matrix.

use fs_exec::Cx;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_solver::op::{CsrOp, LinearOp};

use super::{ClosedTetCavity, HyperTetError, TetCavityGasLoad, checkpoint, finite, invalid};

#[derive(Clone, Copy)]
pub(crate) struct AdiabaticCavity<'a> {
    pub cavity: &'a ClosedTetCavity,
    pub reference: &'a ClosedGasParcel,
    pub exterior_pressure_pa: f64,
    pub controls: GasEntropySolve,
}

impl AdiabaticCavity<'_> {
    pub fn evaluate(&self, displacement: &[f64], cx: &Cx<'_>) -> Result<TetCavityGasLoad, HyperTetError> {
        self.cavity.adiabatic_gas_load(
            displacement, self.reference, self.exterior_pressure_pa, self.controls, cx,
        )
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
mod tests {
    use super::*;
    use fs_sparse::Coo;

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
