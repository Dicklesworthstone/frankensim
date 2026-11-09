//! Diagonal of the retained material action without global assembly.

use fs_exec::Cx;

use super::{
    Element, HyperTetError, PreparedHyperTetTangent, action_error, checkpoint, finish_action,
};

impl PreparedHyperTetTangent<'_> {
    /// Recover the original tangent diagonal for a matrix-free preconditioner.
    ///
    /// Retains every material and geometric term in the stored Piola Hessian.
    /// Values are signed: this observer neither clamps an unstable tangent nor
    /// replaces it with a positive material surrogate. All support rows remain.
    /// No constitutive evaluations, element matrices or allocations are needed.
    /// Scratch may change on refusal; output changes only on success.
    ///
    /// # Errors
    /// Incorrect buffer lengths, unrepresentable diagonal or cancellation.
    pub fn diagonal(
        &self,
        diagonal_n_m: &mut [f64],
        scratch_n_m: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), HyperTetError> {
        checkpoint(cx)?;
        if diagonal_n_m.len() != self.dofs || scratch_n_m.len() != self.dofs {
            return Err(action_error(
                "diagonal and scratch must have the full DOF count",
            ));
        }
        for chunk in scratch_n_m.chunks_mut(4096) {
            checkpoint(cx)?;
            chunk.fill(0.0);
        }
        for (element, tangent) in self.elements.iter().zip(&self.tangents) {
            checkpoint(cx)?;
            element_diagonal(element, tangent, scratch_n_m);
        }
        finish_action(diagonal_n_m, scratch_n_m, cx)
    }
}

fn element_diagonal(element: &Element, tangent: &[[f64; 9]; 9], output: &mut [f64]) {
    for (a, &node) in element.nodes.iter().enumerate() {
        // The forward action differentiates relative nodal increments. Its
        // first input column is therefore the negative sum of the other three
        // gradients, including their original floating-point summation order.
        let input_gradient = if a == 0 {
            core::array::from_fn(|j| {
                element.gradients[1..]
                    .iter()
                    .fold(0.0, |sum, gradient| (-1.0_f64).mul_add(gradient[j], sum))
            })
        } else {
            element.gradients[a]
        };
        for component in 0..3 {
            let force = (0..3).fold(0.0, |sum, j| {
                let piola = (0..3).fold(0.0, |value, k| {
                    tangent[3 * component + j][3 * component + k].mul_add(input_gradient[k], value)
                });
                piola.mul_add(element.gradients[a][j], sum)
            });
            let dof = 3 * node + component;
            output[dof] = element.volume.mul_add(force, output[dof]);
        }
    }
}
