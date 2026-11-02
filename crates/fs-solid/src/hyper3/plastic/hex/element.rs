//! Original Hex8 geometry and variational finite mean dilatation.
//!
//! `Fbar_q = (Jbar/J_q)^(1/3) F_q`, with `Jbar = sum(w_q J_q)/V`.
//! The energy is the original point potential evaluated at these assumed
//! gradients. This owner differentiates the WHOLE element energy: material,
//! geometric and cross-quadrature terms. Full deviatoric quadrature remains;
//! pressure/dilatation condensation does not become one-point stiffness.
//!
//! Material domains apply to the assumed constitutive chart. Every physical
//! corner, centre and quadrature Jacobian is independently admitted, so a
//! positive average cannot hide local inversion. Sampled admission does not
//! certify global injectivity of a warped element.

use crate::linear3::hex::{HEX27_REFERENCE_NODES, hex8_point};
use crate::plastic3::SmallStrainTetError;
use fs_exec::Cx;
use fs_feec::highorder::quad1d::gauss_legendre;
use fs_material::plastic::finite::FiniteJ2State;
use fs_material::plastic::finite::thermal::FrozenThermalFiniteJ2;
use fs_math::det;

type Matrix24 = [[f64; 24]; 24];
type Gradients = [[f64; 3]; 8];

#[derive(Clone, Copy, Debug)]
pub(super) struct ReferencePoint {
    pub shapes: [f64; 8],
    pub gradients: Gradients,
    pub volume_m3: f64,
}

/// One original cell, prepared once and shared across accepted steps.
#[derive(Debug)]
pub(super) struct PreparedHex {
    pub nodes: [usize; 8],
    pub points: [ReferencePoint; 8],
    pub reference_volume_m3: f64,
    probe_gradients: [Gradients; 9],
    element_index: usize,
}

#[derive(Clone, Debug)]
pub(super) struct HexKinematics {
    pub physical_deformation: [[f64; 9]; 8],
    pub assumed_deformation: [[f64; 9]; 8],
    pub mean_volume_ratio: f64,
    jacobians: [f64; 8],
    scales: [f64; 8],
}

#[derive(Debug)]
pub(super) struct EvaluatedHex {
    pub local_force_n: [f64; 24],
    pub tangent_n_m: Option<Matrix24>,
    pub incremental_potential_j: f64,
}

fn invalid(what: impl Into<String>) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainTetError> {
    cx.checkpoint().map_err(|_| SmallStrainTetError::Cancelled)
}

// Match the exact derivative of translation-invariant relative displacement
// evaluation; the reference owner already supplies partition-of-unity gradients.
fn balanced(mut gradients: Gradients) -> Gradients {
    gradients[0] = core::array::from_fn(|i| -gradients[1..].iter().map(|g| g[i]).sum::<f64>());
    gradients
}

impl PreparedHex {
    pub fn new(
        nodes_m: &[[f64; 3]],
        connectivity: [usize; 8],
        element_index: usize,
        minimum_scaled_jacobian: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, SmallStrainTetError> {
        checkpoint(cx)?;
        if connectivity.iter().any(|&node| node >= nodes_m.len())
            || (0..8).any(|i| connectivity[..i].contains(&connectivity[i]))
        {
            return Err(invalid(format!(
                "finite hex {element_index} has missing or repeated nodes"
            )));
        }
        let positions = connectivity.map(|node| nodes_m[node]);
        if positions.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("finite hex coordinates must be finite"));
        }
        let (coordinate, weight) = gauss_legendre(2);
        let mut points = [ReferencePoint {
            shapes: [0.0; 8],
            gradients: [[0.0; 3]; 8],
            volume_m3: 0.0,
        }; 8];
        let mut reference_volume_m3 = 0.0;
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let q = 4 * i + 2 * j + k;
                    let point = hex8_point(
                        positions,
                        [coordinate[i], coordinate[j], coordinate[k]],
                        minimum_scaled_jacobian,
                        element_index,
                        q,
                    )
                    .map_err(|error| invalid(error.to_string()))?;
                    let volume_m3 = point.determinant_m3 * weight[i] * weight[j] * weight[k];
                    if !volume_m3.is_finite() || volume_m3 <= 0.0 {
                        return Err(invalid(
                            "finite hex quadrature volume is not positive and finite",
                        ));
                    }
                    reference_volume_m3 += volume_m3;
                    points[q] = ReferencePoint {
                        shapes: point.shapes,
                        gradients: balanced(point.gradients),
                        volume_m3,
                    };
                }
            }
        }
        let mut probe_gradients = [[[0.0; 3]; 8]; 9];
        for (probe, gradients) in probe_gradients.iter_mut().enumerate() {
            let location = if probe < 8 {
                HEX27_REFERENCE_NODES[probe].map(f64::from)
            } else {
                [0.0; 3]
            };
            *gradients = balanced(
                hex8_point(
                    positions,
                    location,
                    minimum_scaled_jacobian,
                    element_index,
                    8 + probe,
                )
                .map_err(|error| invalid(error.to_string()))?
                .gradients,
            );
        }
        if !reference_volume_m3.is_finite() || reference_volume_m3 <= 0.0 {
            return Err(invalid(
                "finite hex reference volume is not positive and finite",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            nodes: connectivity,
            points,
            reference_volume_m3,
            probe_gradients,
            element_index,
        })
    }

    pub fn kinematics(
        &self,
        local_u: &[f64; 24],
        mean_dilatation: bool,
        minimum_total_jacobian: f64,
        cx: &Cx<'_>,
    ) -> Result<HexKinematics, SmallStrainTetError> {
        checkpoint(cx)?;
        if local_u.iter().any(|x| !x.is_finite())
            || !minimum_total_jacobian.is_finite()
            || minimum_total_jacobian <= 0.0
        {
            return Err(invalid(
                "finite hex displacement and positive Jacobian bound required",
            ));
        }
        let admit = |f: &[f64; 9], point: usize| -> Result<f64, SmallStrainTetError> {
            let j = determinant(f);
            if f.iter().any(|x| !x.is_finite()) || !j.is_finite() || j < minimum_total_jacobian {
                return Err(invalid(format!(
                    "finite hex {} physical Jacobian at probe {point} is {j:e}, minimum {minimum_total_jacobian:e}",
                    self.element_index
                )));
            }
            Ok(j)
        };
        for (probe, gradients) in self.probe_gradients.iter().enumerate() {
            admit(&deformation(local_u, gradients), 8 + probe)?;
        }
        let mut physical_deformation = [[0.0; 9]; 8];
        let mut jacobians = [0.0; 8];
        let mut minimum_volume_ratio = f64::INFINITY;
        for q in 0..8 {
            physical_deformation[q] = deformation(local_u, &self.points[q].gradients);
            jacobians[q] = admit(&physical_deformation[q], q)?;
            minimum_volume_ratio = minimum_volume_ratio.min(jacobians[q]);
        }
        // A shifted positive weighted mean preserves a constant Jacobian
        // exactly, including the virgin reference state at arbitrarily stiff K.
        let mean_volume_ratio = minimum_volume_ratio
            + self
                .points
                .iter()
                .zip(jacobians)
                .fold(0.0, |sum, (point, j)| {
                    (point.volume_m3 / self.reference_volume_m3)
                        .mul_add(j - minimum_volume_ratio, sum)
                });
        if !mean_volume_ratio.is_finite() || mean_volume_ratio < minimum_total_jacobian {
            return Err(invalid("finite hex mean volume ratio is inadmissible"));
        }
        let mut assumed_deformation = physical_deformation;
        let mut scales = [1.0; 8];
        if mean_dilatation {
            let log_mean = det::ln(mean_volume_ratio);
            for q in 0..8 {
                scales[q] = det::exp((log_mean - det::ln(jacobians[q])) / 3.0);
                if !scales[q].is_finite() || scales[q] <= 0.0 {
                    return Err(invalid(
                        "finite hex mean-dilatation scale is not representable",
                    ));
                }
                assumed_deformation[q] = physical_deformation[q].map(|x| scales[q] * x);
                let assumed_j = determinant(&assumed_deformation[q]);
                if assumed_deformation[q].iter().any(|x| !x.is_finite())
                    || !assumed_j.is_finite()
                    || assumed_j <= 0.0
                {
                    return Err(invalid(
                        "finite hex assumed deformation is not representable",
                    ));
                }
            }
        }
        Ok(HexKinematics {
            physical_deformation,
            assumed_deformation,
            mean_volume_ratio,
            jacobians,
            scales,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        frozen: &FrozenThermalFiniteJ2,
        old: &[FiniteJ2State],
        local_u: &[f64; 24],
        mean_dilatation: bool,
        minimum_mechanical_jacobian: f64,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<EvaluatedHex, SmallStrainTetError> {
        if old.len() != 8 {
            return Err(invalid(
                "finite hex requires exactly eight accepted quadrature histories",
            ));
        }
        let theta = frozen.thermal_stretch();
        let state = self.kinematics(
            local_u,
            mean_dilatation,
            minimum_mechanical_jacobian * theta * theta * theta,
            cx,
        )?;
        let mut log_gradients = [[0.0; 24]; 8];
        let mut mean_log_gradient = [0.0; 24];
        let mut mean_log_hessian = [[0.0; 24]; 24];
        if mean_dilatation {
            for q in 0..8 {
                let inverse = inverse(&state.physical_deformation[q], state.jacobians[q])?;
                let gradient = &mut log_gradients[q];
                for a in 0..8 {
                    for i in 0..3 {
                        gradient[3 * a + i] = (0..3)
                            .map(|j| self.points[q].gradients[a][j] * inverse[3 * j + i])
                            .sum();
                    }
                }
                let fraction = (self.points[q].volume_m3 / self.reference_volume_m3)
                    * (state.jacobians[q] / state.mean_volume_ratio);
                for i in 0..24 {
                    mean_log_gradient[i] = fraction.mul_add(gradient[i], mean_log_gradient[i]);
                    if need_tangent {
                        for j in 0..24 {
                            let h = log_hessian(gradient, i, j);
                            mean_log_hessian[i][j] = fraction
                                .mul_add(gradient[i] * gradient[j] + h, mean_log_hessian[i][j]);
                        }
                    }
                }
            }
            if need_tangent {
                for i in 0..24 {
                    for j in 0..24 {
                        mean_log_hessian[i][j] -= mean_log_gradient[i] * mean_log_gradient[j];
                    }
                }
            }
        }
        let mut local_force_n = [0.0; 24];
        let mut tangent_n_m = need_tangent.then_some([[0.0; 24]; 24]);
        let mut incremental_potential_j = 0.0;
        for q in 0..8 {
            checkpoint(cx)?;
            let full = if need_tangent {
                Some(
                    frozen
                        .response(&state.assumed_deformation[q], &old[q])
                        .map_err(|error| {
                            invalid(format!(
                                "finite hex {} point {q}: {error}",
                                self.element_index
                            ))
                        })?,
                )
            } else {
                None
            };
            let (piola, potential) = if let Some(response) = &full {
                (response.first_piola_pa, response.incremental_potential_j_m3)
            } else {
                frozen
                    .stress_and_potential(&state.assumed_deformation[q], &old[q])
                    .map_err(|error| {
                        invalid(format!(
                            "finite hex {} point {q}: {error}",
                            self.element_index
                        ))
                    })?
            };
            let point = &self.points[q];
            let f = &state.physical_deformation[q];
            let scale = state.scales[q];
            let b: [f64; 24] = core::array::from_fn(|i| {
                if mean_dilatation {
                    (mean_log_gradient[i] - log_gradients[q][i]) / 3.0
                } else {
                    0.0
                }
            });
            let columns: [[f64; 9]; 24] = core::array::from_fn(|i| {
                core::array::from_fn(|k| {
                    let direct = if k / 3 == i % 3 {
                        point.gradients[i / 3][k % 3]
                    } else {
                        0.0
                    };
                    scale * b[i].mul_add(f[k], direct)
                })
            });
            incremental_potential_j = point.volume_m3.mul_add(potential, incremental_potential_j);
            for i in 0..24 {
                local_force_n[i] = point
                    .volume_m3
                    .mul_add(dot(&piola, &columns[i]), local_force_n[i]);
            }
            if let (Some(tangent), Some(response)) = (&mut tangent_n_m, &full) {
                // A-times-column once, then small dot products. No dense
                // 72-point-DOF Hessian or per-point matrix allocation is needed.
                let material_columns: [[f64; 9]; 24] = core::array::from_fn(|j| {
                    core::array::from_fn(|k| dot(&response.tangent_pa[k], &columns[j]))
                });
                let piola_f = dot(&piola, f);
                let piola_b: [f64; 24] = core::array::from_fn(|i| {
                    (0..3)
                        .map(|k| piola[3 * (i % 3) + k] * point.gradients[i / 3][k])
                        .sum()
                });
                for i in 0..24 {
                    for j in 0..24 {
                        let mut value = dot(&columns[i], &material_columns[j]);
                        if mean_dilatation {
                            let second_b = (mean_log_hessian[i][j]
                                - log_hessian(&log_gradients[q], i, j))
                                / 3.0;
                            value += scale
                                * (b[j] * piola_b[i]
                                    + b[i] * piola_b[j]
                                    + (b[i] * b[j] + second_b) * piola_f);
                        }
                        tangent[i][j] = point.volume_m3.mul_add(value, tangent[i][j]);
                    }
                }
            }
        }
        if !incremental_potential_j.is_finite()
            || local_force_n.iter().any(|x| !x.is_finite())
            || tangent_n_m
                .as_ref()
                .is_some_and(|matrix| matrix.iter().flatten().any(|x| !x.is_finite()))
        {
            return Err(invalid(
                "finite hex potential, force or full tangent is not representable",
            ));
        }
        checkpoint(cx)?;
        Ok(EvaluatedHex {
            local_force_n,
            tangent_n_m,
            incremental_potential_j,
        })
    }
}

fn log_hessian(gradient: &[f64; 24], i: usize, j: usize) -> f64 {
    -gradient[3 * (i / 3) + j % 3] * gradient[3 * (j / 3) + i % 3]
}

fn deformation(u: &[f64; 24], gradients: &Gradients) -> [f64; 9] {
    let mut f = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    for a in 1..8 {
        for i in 0..3 {
            let delta = u[3 * a + i] - u[i];
            for j in 0..3 {
                f[3 * i + j] = delta.mul_add(gradients[a][j], f[3 * i + j]);
            }
        }
    }
    f
}

fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (x, y)| x.mul_add(*y, sum))
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn inverse(f: &[f64; 9], jacobian: f64) -> Result<[f64; 9], SmallStrainTetError> {
    let matrix = [
        f[4] * f[8] - f[5] * f[7],
        f[2] * f[7] - f[1] * f[8],
        f[1] * f[5] - f[2] * f[4],
        f[5] * f[6] - f[3] * f[8],
        f[0] * f[8] - f[2] * f[6],
        f[2] * f[3] - f[0] * f[5],
        f[3] * f[7] - f[4] * f[6],
        f[1] * f[6] - f[0] * f[7],
        f[0] * f[4] - f[1] * f[3],
    ]
    .map(|x| x / jacobian);
    if matrix.iter().any(|x| !x.is_finite()) {
        return Err(invalid(
            "finite hex inverse deformation is not representable",
        ));
    }
    Ok(matrix)
}
