//! Conservative scalar weak form; original velocity and quadrature owners.
//!
//! SIP signs follow the symmetric consistency/adjoint-consistency form:
//! https://people.tamu.edu/~guermond/M661_FALL_2025/chap38.pdf, section 38.2.

use super::*;
use crate::bdm::{edge_gauss_pub, eval_basis};
use crate::ns::convection_edge_quadrature;

pub(in crate::ns) struct ScalarBasis {
    origin: [f64; 2],
    pub(super) gradient: [[f64; 2]; 3],
}
impl ScalarBasis {
    pub(in crate::ns) fn new(p: [[f64; 2]; 3], area: f64) -> Self {
        Self {
            origin: p[0],
            gradient: core::array::from_fn(|i| {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                [
                    (p[j][1] - p[k][1]) / (2.0 * area),
                    (p[k][0] - p[j][0]) / (2.0 * area),
                ]
            }),
        }
    }
    pub(in crate::ns) fn value(&self, p: [f64; 2]) -> [f64; 3] {
        let d = [p[0] - self.origin[0], p[1] - self.origin[1]];
        let b = dot(self.gradient[1], d);
        let c = dot(self.gradient[2], d);
        [1.0 - b - c, b, c]
    }
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn trace(samples: [f64; 2], fraction: f64) -> f64 {
    (1.0 - fraction) * samples[0] + fraction * samples[1]
}

impl<'m> FluxSystem<'m> {
    /// Bind a passive scalar to the original fixed BDM mesh and prescribed
    /// divergence-free velocity. `diffusivity_m2_s` is cellwise molecular or
    /// explicitly combined molecular+modeled diffusivity, finite and nonnegative.
    /// Boundary callbacks are sampled once at the original two Gauss points;
    /// the same affine trace supplies every sign-split inflow quadrature point.
    ///
    /// Interior diffusion uses harmonic alpha/normal-height weights and SIP;
    /// no orthogonal centroid connector or guessed diffusion distance is used.
    /// The original affine BDM velocity supplies conservative volume advection.
    /// A fixed mesh and constant storage density are required. No material,
    /// momentum, thermodynamic or scalar state is modified by preparation.
    pub fn prepare_scalar_transport(
        &self,
        velocity: &[f64],
        diffusivity_m2_s: &[f64],
        boundaries: &[ScalarBoundary<'_>],
        controls: ScalarTransportControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedScalarTransport<'m>> {
        poll(&mut cancelled)?;
        let mesh = self.mesh;
        let nc = mesh.tris.len();
        let n = nc
            .checked_mul(3)
            .ok_or(ScalarTransportError::Invalid("scalar dimension"))?;
        if nc == 0
            || diffusivity_m2_s.len() != nc
            || diffusivity_m2_s.iter().any(|a| !a.is_finite() || *a < 0.0)
            || !controls.penalty.is_finite()
            || controls.penalty < 12.0
            || !controls.relative_volume_divergence_tolerance.is_finite()
            || !(0.0..=1e-12).contains(&controls.relative_volume_divergence_tolerance)
            || mesh.areas.len() != nc
            || mesh.tri_edges.len() != nc
            || self.bases.len() != nc
        {
            return Err(ScalarTransportError::Invalid(
                "nonnegative cell diffusivities, P1 penalty >=12 and bounded continuity allowance required",
            ));
        }
        let rates = self.edge_volume_transport(velocity)?;
        poll(&mut cancelled)?;
        let mut owners = vec![None; mesh.edges.len()];
        for &boundary in boundaries {
            let edge = match boundary {
                ScalarBoundary::Value { edge, .. } | ScalarBoundary::Flux { edge, .. } => edge,
            };
            if edge >= owners.len() || !mesh.is_boundary(edge) || owners[edge].is_some() {
                return Err(ScalarTransportError::Invalid(
                    "one scalar owner per original boundary edge required",
                ));
            }
            owners[edge] = Some(boundary);
        }
        let mut bases = Vec::with_capacity(nc);
        for (cell, tri) in mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let area = mesh.areas[cell];
            if !area.is_finite() || area <= 0.0 {
                return Err(ScalarTransportError::Invalid(
                    "positive original triangle areas required",
                ));
            }
            let net = sum(mesh.tri_edges[cell]
                .iter()
                .map(|&(edge, sign)| sign * rates[edge].net_m2_s()));
            let streams = sum(mesh.tri_edges[cell].iter().map(|&(edge, _)| {
                rates[edge].along_normal_m2_s() + rates[edge].against_normal_m2_s()
            }));
            let relative = if streams == 0.0 {
                net.abs()
            } else {
                net.abs() / streams
            };
            if !relative.is_finite() || relative > controls.relative_volume_divergence_tolerance {
                return Err(ScalarTransportError::VolumeDivergence {
                    cell,
                    relative_residual: relative,
                });
            }
            let basis = ScalarBasis::new(tri.map(|vertex| mesh.verts[vertex]), area);
            if basis.gradient.iter().flatten().any(|v| !v.is_finite()) {
                return Err(ScalarTransportError::Arithmetic("scalar basis"));
            }
            bases.push(basis);
        }
        let mut matrix = Coo::new(n, n);
        let mut boundary_load = vec![0.0; n];
        let mut boundary_fluxes = Vec::with_capacity(boundaries.len());
        for (cell, tri) in mesh.tris.iter().enumerate() {
            poll(&mut cancelled)?;
            let basis = &bases[cell];
            let alpha = diffusivity_m2_s[cell];
            for i in 0..3 {
                for j in 0..3 {
                    matrix.push(
                        3 * cell + i,
                        3 * cell + j,
                        mesh.areas[cell] * alpha * dot(basis.gradient[i], basis.gradient[j]),
                    );
                }
            }
            for (point, weight) in tri_quad(tri.map(|vertex| mesh.verts[vertex]), mesh.areas[cell])
            {
                let phi = basis.value(point);
                let mut u = [0.0; 2];
                for local in 0..6 {
                    let dof = 2 * mesh.tri_edges[cell][local / 2].0 + local % 2;
                    let v = eval_basis(&self.bases[cell], local, point);
                    for axis in 0..2 {
                        u[axis] += velocity[dof] * v[axis];
                    }
                }
                for i in 0..3 {
                    for j in 0..3 {
                        matrix.push(
                            3 * cell + i,
                            3 * cell + j,
                            -weight * dot(u, basis.gradient[i]) * phi[j],
                        );
                    }
                }
            }
        }
        for (edge_id, edge) in mesh.edges.iter().enumerate() {
            poll(&mut cancelled)?;
            let (left, right) = edge.tris;
            let a = mesh.verts[edge.verts.0];
            let b = mesh.verts[edge.verts.1];
            let gauss = edge_gauss_pub(a, b);
            let moments = [velocity[2 * edge_id], velocity[2 * edge_id + 1]];
            let normal_rate = |fraction: f64| {
                moments[0] + (12.0 / 3.0_f64.sqrt()) * moments[1] * (fraction - 0.5)
            };
            let d_left = finite(
                diffusivity_m2_s[left] * (edge.len / (2.0 * mesh.areas[left])),
                "scalar face diffusivity",
            )?;
            if right != usize::MAX {
                let d_right = finite(
                    diffusivity_m2_s[right] * (edge.len / (2.0 * mesh.areas[right])),
                    "scalar face diffusivity",
                )?;
                let (weights, harmonic) = if d_left == 0.0 && d_right == 0.0 {
                    ([0.5, 0.5], 0.0)
                } else if d_left <= d_right {
                    let ratio = d_left / d_right;
                    (
                        [1.0 / (1.0 + ratio), ratio / (1.0 + ratio)],
                        d_left * (2.0 / (1.0 + ratio)),
                    )
                } else {
                    let ratio = d_right / d_left;
                    (
                        [ratio / (1.0 + ratio), 1.0 / (1.0 + ratio)],
                        d_right * (2.0 / (1.0 + ratio)),
                    )
                };
                // Form BOTH weights directly: subtracting one from the larger
                // loses the smaller side's finite flux at large contrast.
                let penalty = finite(controls.penalty * harmonic, "scalar SIP penalty")?;
                for (point, weight) in gauss {
                    let phi = [bases[left].value(point), bases[right].value(point)];
                    for (r, cell_r) in [left, right].into_iter().enumerate() {
                        let sign_r = if r == 0 { 1.0 } else { -1.0 };
                        for (c, cell_c) in [left, right].into_iter().enumerate() {
                            let sign_c = if c == 0 { 1.0 } else { -1.0 };
                            for i in 0..3 {
                                for j in 0..3 {
                                    let row_flux = weights[r]
                                        * diffusivity_m2_s[cell_r]
                                        * dot(bases[cell_r].gradient[i], edge.normal);
                                    let column_flux = weights[c]
                                        * diffusivity_m2_s[cell_c]
                                        * dot(bases[cell_c].gradient[j], edge.normal);
                                    matrix.push(
                                        3 * cell_r + i,
                                        3 * cell_c + j,
                                        weight
                                            * (-sign_r * phi[r][i] * column_flux
                                                - row_flux * sign_c * phi[c][j]
                                                + penalty
                                                    * sign_r
                                                    * sign_c
                                                    * phi[r][i]
                                                    * phi[c][j]),
                                    );
                                }
                            }
                        }
                    }
                }
                for (point, weight, fraction) in convection_edge_quadrature(a, b, moments) {
                    let normal = finite(normal_rate(fraction), "scalar normal velocity")?;
                    let donor = if normal >= 0.0 { left } else { right };
                    let donor_phi = bases[donor].value(point);
                    for (cell, sign) in [(left, 1.0), (right, -1.0)] {
                        let phi = bases[cell].value(point);
                        for i in 0..3 {
                            for j in 0..3 {
                                matrix.push(
                                    3 * cell + i,
                                    3 * donor + j,
                                    weight * sign * normal * phi[i] * donor_phi[j],
                                );
                            }
                        }
                    }
                }
                continue;
            }
            let owner = owners[edge_id].ok_or(ScalarTransportError::Invalid(
                "missing scalar boundary owner",
            ))?;
            let mut flux = BoundaryFlux {
                edge: edge_id,
                cell: left,
                advection: [0.0; 3],
                diffusion: [0.0; 3],
                advective_load: 0.0,
                diffusive_load: 0.0,
                reservoir_coefficients: [[0.0; 2]; 3],
                reservoir_values: [0.0; 2],
                prescribed_flux: false,
            };
            let (value, outward, inflow) = match owner {
                ScalarBoundary::Value { value, .. } => (Some(value), None, Some(value)),
                ScalarBoundary::Flux {
                    outward, inflow, ..
                } => (None, Some(outward), inflow),
            };
            let sample = |callback: &dyn Fn([f64; 2]) -> f64| -> Result<[f64; 2]> {
                Ok([
                    finite(callback(gauss[0].0), "scalar boundary trace")?,
                    finite(callback(gauss[1].0), "scalar boundary trace")?,
                ])
            };
            let value_samples = value.map(sample).transpose()?;
            let incoming = if rates[edge_id].against_normal_m2_s() > 0.0 {
                Some(if let Some(value) = value_samples {
                    value
                } else {
                    sample(inflow.ok_or(ScalarTransportError::MissingInflow { edge: edge_id })?)?
                })
            } else {
                None
            };
            let flux_samples = outward.map(sample).transpose()?;
            flux.reservoir_values = value_samples.or(incoming).unwrap_or([0.0; 2]);
            flux.prescribed_flux =
                flux_samples.is_some_and(|values| values.iter().any(|&v| v != 0.0));
            let penalty = finite(controls.penalty * d_left, "scalar boundary penalty")?;
            for (q, (point, weight)) in gauss.into_iter().enumerate() {
                let phi = bases[left].value(point);
                if let Some(values) = value_samples {
                    let g = values[q];
                    for i in 0..3 {
                        let row_flux =
                            diffusivity_m2_s[left] * dot(bases[left].gradient[i], edge.normal);
                        boundary_load[3 * left + i] += weight * (-row_flux + penalty * phi[i]) * g;
                        flux.reservoir_coefficients[i][q] +=
                            weight * (-row_flux + penalty * phi[i]);
                        flux.diffusion[i] += weight * (-row_flux + penalty * phi[i]);
                        for j in 0..3 {
                            let col_flux =
                                diffusivity_m2_s[left] * dot(bases[left].gradient[j], edge.normal);
                            matrix.push(
                                3 * left + i,
                                3 * left + j,
                                weight
                                    * (-phi[i] * col_flux - row_flux * phi[j]
                                        + penalty * phi[i] * phi[j]),
                            );
                        }
                    }
                    flux.diffusive_load -= weight * penalty * g;
                } else if let Some(values) = flux_samples {
                    for i in 0..3 {
                        boundary_load[3 * left + i] -= weight * values[q] * phi[i];
                    }
                    flux.diffusive_load += weight * values[q];
                }
            }
            for (point, weight, fraction) in convection_edge_quadrature(a, b, moments) {
                let normal = finite(normal_rate(fraction), "scalar normal velocity")?;
                let phi = bases[left].value(point);
                if normal >= 0.0 {
                    for i in 0..3 {
                        flux.advection[i] += weight * normal * phi[i];
                        for j in 0..3 {
                            matrix.push(
                                3 * left + i,
                                3 * left + j,
                                weight * normal * phi[i] * phi[j],
                            );
                        }
                    }
                } else {
                    let donor = finite(
                        trace(
                            incoming
                                .ok_or(ScalarTransportError::MissingInflow { edge: edge_id })?,
                            fraction,
                        ),
                        "scalar inflow trace",
                    )?;
                    for i in 0..3 {
                        boundary_load[3 * left + i] -= weight * normal * donor * phi[i];
                        let coefficient = -weight * normal * phi[i];
                        flux.reservoir_coefficients[i][0] += coefficient * (1.0 - fraction);
                        flux.reservoir_coefficients[i][1] += coefficient * fraction;
                    }
                    flux.advective_load += weight * normal * donor;
                }
            }
            boundary_fluxes.push(flux);
        }
        let operator = matrix.assemble();
        if boundary_load.iter().any(|v| !v.is_finite())
            || (0..n).any(|row| operator.row(row).1.iter().any(|v| !v.is_finite()))
        {
            return Err(ScalarTransportError::Arithmetic(
                "assembled scalar transport",
            ));
        }
        poll(&mut cancelled)?;
        Ok(PreparedScalarTransport {
            mesh,
            operator,
            boundary_load,
            boundary_fluxes,
            diffusivity_m2_s: diffusivity_m2_s.to_vec(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn high_contrast_interface_retains_both_one_sided_physical_fluxes() {
        let mesh = TriMesh::from_mesh2(&fs_solid::Mesh2::triangles(1.0, 1.0, 1, 1));
        let system = FluxSystem::new(&mesh);
        let zero = |_: [f64; 2]| 0.0;
        let boundaries: Vec<_> = mesh
            .edges
            .iter()
            .enumerate()
            .filter(|(edge, _)| mesh.is_boundary(*edge))
            .map(|(edge, _)| ScalarBoundary::Flux {
                edge,
                outward: &zero,
                inflow: None,
            })
            .collect();
        let (edge_id, edge) = mesh
            .edges
            .iter()
            .enumerate()
            .find(|(edge, _)| !mesh.is_boundary(*edge))
            .unwrap();
        for diffusivity in [[1e-18, 1.0], [1.0, 1e-18]] {
            let plan = system
                .prepare_scalar_transport(
                    &vec![0.0; system.n_u],
                    &diffusivity,
                    &boundaries,
                    ScalarTransportControls::default(),
                    || false,
                )
                .unwrap();
            let heights = [
                2.0 * mesh.areas[edge.tris.0] / edge.len,
                2.0 * mesh.areas[edge.tris.1] / edge.len,
            ];
            // Independent resistance-in-series conductance. The trial basis
            // opposite the interface has zero trace there. A constant test in
            // the neighbor removes all volume, penalty and adjoint terms.
            let expected = -edge.len
                / (heights[0] / diffusivity[edge.tris.0] + heights[1] / diffusivity[edge.tris.1]);
            for (test, trial) in [(edge.tris.0, edge.tris.1), (edge.tris.1, edge.tris.0)] {
                let opposite = mesh.tri_edges[trial]
                    .iter()
                    .position(|&(id, _)| id == edge_id)
                    .unwrap();
                let column = 3 * trial + opposite;
                let actual: f64 = (0..3)
                    .map(|i| {
                        let (columns, values) = plan.operator.row(3 * test + i);
                        columns
                            .iter()
                            .zip(values)
                            .filter(|(j, _)| **j == column)
                            .map(|(_, value)| *value)
                            .sum::<f64>()
                    })
                    .sum();
                assert!(
                    (actual / expected - 1.0).abs() < 2e-13,
                    "alpha={diffusivity:?}, test={test}, actual={actual:e}, expected={expected:e}"
                );
            }
        }
    }
}
