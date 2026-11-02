//! Original P2 reference geometry with objective two-point material work.
//!
//! The endpoint Jacobian is generally nonsymmetric. Each local derivative is
//! retained in full; constrained rows remain available for support work. The
//! centered consistent mass is added before the single global CSR assembly.

use fs_exec::Cx;
use fs_sparse::Coo;

use super::{
    ConservativeTetError, ConservativeTetSolid, Evaluation, checkpoint, deformation, inertia_error,
    invalid, material_error,
};
use crate::inertia::MidpointInertia;

impl ConservativeTetSolid<'_> {
    pub(super) fn evaluate(
        &self,
        u: &[[f64; 3]],
        with_tangent: bool,
        inertia: &MidpointInertia<'_>,
        cx: &Cx<'_>,
    ) -> Result<Evaluation, ConservativeTetError> {
        checkpoint(cx)?;
        if u.len() != self.mesh.nodes_count {
            return Err(invalid(
                "conservative displacement must cover the admitted mesh",
            ));
        }
        for (node, value) in u.iter().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            if value.iter().any(|v| !v.is_finite()) {
                return Err(invalid("conservative displacement must be finite"));
            }
        }
        let nu = 3 * self.mesh.nodes_count;
        let mut forces_n = vec![0.0; nu];
        let mut matrix = with_tangent.then(|| Coo::new(nu, nu));
        let mut stored_j = 0.0;
        let mut material_work_j = 0.0;
        let mut min_det_f = f64::INFINITY;
        for (e, element) in self.mesh.elements.iter().enumerate() {
            checkpoint(cx)?;
            let dofs: [usize; 30] =
                core::array::from_fn(|d| 3 * element.connectivity[d / 3] + d % 3);
            let mut local_r = [0.0; 30];
            let mut local_k = [[0.0; 30]; 30];
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let f = deformation(element, point, u);
                let frozen = &self.points[e * self.mesh.points_per_element + q];
                let response = if with_tangent {
                    Some(
                        frozen
                            .response(&f)
                            .map_err(|error| material_error(e, q, error))?,
                    )
                } else {
                    None
                };
                let state = match &response {
                    Some(response) => response.state,
                    None => frozen
                        .state(&f)
                        .map_err(|error| material_error(e, q, error))?,
                };
                let determinant = f[0] * (f[4] * f[8] - f[5] * f[7])
                    - f[1] * (f[3] * f[8] - f[5] * f[6])
                    + f[2] * (f[3] * f[7] - f[4] * f[6]);
                min_det_f = min_det_f.min(determinant);
                stored_j = point.weight_m3.mul_add(state.energy_after_pa, stored_j);
                material_work_j = point.weight_m3.mul_add(state.work_pa, material_work_j);
                if !stored_j.is_finite()
                    || !material_work_j.is_finite()
                    || !determinant.is_finite()
                    || determinant <= 0.0
                {
                    return Err(invalid(
                        "conservative material storage, work or volume is unrepresentable",
                    ));
                }
                for a in 0..10 {
                    for i in 0..3 {
                        let force = (0..3).fold(0.0, |sum, j| {
                            state.first_piola_pa[3 * i + j].mul_add(point.gradients[a][j], sum)
                        });
                        local_r[3 * a + i] = point.weight_m3.mul_add(force, local_r[3 * a + i]);
                    }
                }
                if let Some(response) = response {
                    // Reuse each right shape contraction for every test shape:
                    // 9*30*3 + 30*30*3 products instead of 30*30*9.
                    let a_gradient: [[f64; 30]; 9] = core::array::from_fn(|ij| {
                        core::array::from_fn(|bk| {
                            (0..3).fold(0.0, |sum, l| {
                                response.endpoint_tangent_pa[ij][3 * (bk % 3) + l]
                                    .mul_add(point.gradients[bk / 3][l], sum)
                            })
                        })
                    });
                    for row in 0..30 {
                        let (a, i) = (row / 3, row % 3);
                        for col in 0..30 {
                            let value = (0..3).fold(0.0, |sum, j| {
                                point.gradients[a][j].mul_add(a_gradient[3 * i + j][col], sum)
                            });
                            local_k[row][col] = point.weight_m3.mul_add(value, local_k[row][col]);
                        }
                    }
                }
            }
            if local_r
                .iter()
                .chain(local_k.iter().flatten())
                .any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "conservative element force or tangent is unrepresentable",
                ));
            }
            // One scatter per element, retaining both off-diagonal directions.
            for row in 0..30 {
                forces_n[dofs[row]] += local_r[row];
                if let Some(coo) = matrix.as_mut() {
                    for col in 0..30 {
                        coo.push(dofs[row], dofs[col], local_k[row][col]);
                    }
                }
            }
        }
        inertia
            .add_force(u.as_flattened(), &mut forces_n, cx)
            .map_err(inertia_error)?;
        if let Some(coo) = matrix.as_mut() {
            inertia.add_tangent(coo, cx).map_err(inertia_error)?;
        }
        if forces_n.iter().any(|v| !v.is_finite()) || !min_det_f.is_finite() {
            return Err(invalid(
                "assembled conservative force or volume is unrepresentable",
            ));
        }
        checkpoint(cx)?;
        let tangent = matrix.map(|coo| coo.assemble());
        if let Some(k) = &tangent {
            for row in 0..k.nrows() {
                if row % 128 == 0 {
                    checkpoint(cx)?;
                }
                if k.row(row).1.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("assembled conservative tangent is unrepresentable"));
                }
            }
        }
        checkpoint(cx)?;
        Ok(Evaluation {
            forces_n,
            tangent,
            stored_j,
            material_work_j,
            min_det_f,
        })
    }
}
