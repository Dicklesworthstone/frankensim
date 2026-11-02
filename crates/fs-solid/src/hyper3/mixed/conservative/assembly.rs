//! Original P2 reference geometry with objective two-point material work.
//!
//! The endpoint Jacobian is generally nonsymmetric. Each local derivative is
//! retained in full; constrained rows remain available for support work. The
//! centered consistent mass is added before the single global CSR assembly.

use fs_exec::Cx;
use fs_sparse::Coo;

use super::{
    ConservativeTetError, ConservativeTetSolid, Evaluation, checkpoint, deformation, determinant,
    inertia_error, invalid, material_error,
};
use crate::inertia::MidpointInertia;

impl ConservativeTetSolid<'_> {
    pub(super) fn evaluate(
        &self,
        u: &[[f64; 3]],
        pressure_pa: &[f64],
        with_tangent: bool,
        inertia: &MidpointInertia<'_>,
        cx: &Cx<'_>,
    ) -> Result<Evaluation, ConservativeTetError> {
        checkpoint(cx)?;
        if u.len() != self.mesh.nodes_count
            || pressure_pa.len() != self.algorithmic_pressure_pa.len()
            || pressure_pa.iter().any(|p| !p.is_finite())
        {
            return Err(invalid(
                "finite conservative displacement and pressure must cover the admitted mesh",
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
        let np = pressure_pa.len();
        let local_dofs = if np == 0 { 30 } else { 34 };
        let volume_controls = self.definition.formulation.volume_controls();
        let mut forces_n = vec![0.0; nu];
        let mut pressure_residual_m3 = vec![0.0; np];
        let mut matrix = with_tangent.then(|| Coo::new(nu + np, nu + np));
        let mut stored_j = 0.0;
        let mut material_work_j = 0.0;
        let mut pressure_work_j = 0.0;
        let mut min_det_f = f64::INFINITY;
        for (e, element) in self.mesh.elements.iter().enumerate() {
            checkpoint(cx)?;
            let dofs: [usize; 34] = core::array::from_fn(|d| {
                if d < 30 {
                    3 * element.connectivity[d / 3] + d % 3
                } else {
                    nu + element.pressure_indices[d - 30]
                }
            });
            let mut local_r = [0.0; 34];
            let mut local_k = [[0.0; 34]; 34];
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
                let volume_response = match volume_controls {
                    Some((_, _, tolerance)) if with_tangent => Some(
                        frozen
                            .volume_response(&f, tolerance)
                            .map_err(|error| material_error(e, q, error))?,
                    ),
                    _ => None,
                };
                let volume_state = match (&volume_response, volume_controls) {
                    (Some(response), _) => Some(response.state),
                    (_, Some((_, _, tolerance))) => Some(
                        frozen
                            .volume_state(&f, tolerance)
                            .map_err(|error| material_error(e, q, error))?,
                    ),
                    _ => None,
                };
                let pressure = if np == 0 {
                    0.0
                } else {
                    (0..4).fold(0.0, |sum, a| {
                        point.pressure_shapes[a]
                            .mul_add(pressure_pa[element.pressure_indices[a]], sum)
                    })
                };
                let det_f = determinant(&f);
                min_det_f = min_det_f.min(det_f);
                stored_j = point.weight_m3.mul_add(state.energy_after_pa, stored_j);
                material_work_j = point.weight_m3.mul_add(state.work_pa, material_work_j);
                let stress: [f64; 9] = core::array::from_fn(|ij| {
                    volume_state
                        .as_ref()
                        .map_or(state.first_piola_pa[ij], |volume| {
                            pressure.mul_add(volume.discrete_cofactor[ij], state.first_piola_pa[ij])
                        })
                });
                if let Some(volume) = &volume_state {
                    pressure_work_j = point
                        .weight_m3
                        .mul_add(pressure * volume.volume_work, pressure_work_j);
                    for a in 0..4 {
                        local_r[30 + a] +=
                            point.weight_m3 * point.pressure_shapes[a] * (det_f - 1.0);
                    }
                }
                if !stored_j.is_finite()
                    || !material_work_j.is_finite()
                    || !pressure.is_finite()
                    || !pressure_work_j.is_finite()
                    || !det_f.is_finite()
                    || det_f <= 0.0
                {
                    return Err(invalid(
                        "conservative material storage, work or volume is unrepresentable",
                    ));
                }
                for a in 0..10 {
                    for i in 0..3 {
                        let force = (0..3).fold(0.0, |sum, j| {
                            stress[3 * i + j].mul_add(point.gradients[a][j], sum)
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
                                let kl = 3 * (bk % 3) + l;
                                let value = volume_response.as_ref().map_or(
                                    response.endpoint_tangent_pa[ij][kl],
                                    |volume| {
                                        pressure.mul_add(
                                            volume.endpoint_jacobian[ij][kl],
                                            response.endpoint_tangent_pa[ij][kl],
                                        )
                                    },
                                );
                                value.mul_add(point.gradients[bk / 3][l], sum)
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
                        if let Some(volume) = &volume_state {
                            let algorithmic = (0..3).fold(0.0, |s, j| {
                                volume.discrete_cofactor[3 * i + j]
                                    .mul_add(point.gradients[a][j], s)
                            });
                            let endpoint = (0..3).fold(0.0, |s, j| {
                                volume.endpoint_cofactor[3 * i + j]
                                    .mul_add(point.gradients[a][j], s)
                            });
                            for b in 0..4 {
                                let weight = point.weight_m3 * point.pressure_shapes[b];
                                // Momentum uses the work-conjugate volume gradient;
                                // the endpoint constraint has its own derivative.
                                local_k[row][30 + b] =
                                    weight.mul_add(algorithmic, local_k[row][30 + b]);
                                local_k[30 + b][row] =
                                    weight.mul_add(endpoint, local_k[30 + b][row]);
                            }
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
            for row in 0..local_dofs {
                if row < 30 {
                    forces_n[dofs[row]] += local_r[row];
                } else {
                    pressure_residual_m3[dofs[row] - nu] += local_r[row];
                }
                if let Some(coo) = matrix.as_mut() {
                    for col in 0..local_dofs {
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
        if forces_n
            .iter()
            .chain(&pressure_residual_m3)
            .any(|v| !v.is_finite())
            || !min_det_f.is_finite()
        {
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
            pressure_work_j,
            pressure_residual_m3,
            min_det_f,
        })
    }
}
