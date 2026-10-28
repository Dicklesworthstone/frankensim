//! Initial-stress vibration about an accepted quadratic thermal equilibrium.
//!
//! Keep the reference geometry, mass and homogeneous incremental supports.
//! Integrate K_g[(a,i),(b,j)] = delta_ij grad(N_a)^T sigma_0 grad(N_b)
//! at every original stiffness quadrature point. The actual quadratic strain
//! and stress enter this integral; element-mean stress would lose correlations
//! with the varying gradients, especially on curved isoparametric cells.
//!
//! This extends the same conservative initial-stress approximation as
//! crate::prestress. It is not a finite-strain tangent, follower-load model,
//! stability certificate, mode transfer or post-buckling continuation.

use super::super::{Point, mandel_b};
use super::{QuadraticThermalResponse, elastic_response_at};
use crate::linear3::{
    QuadraticThermalStrainField, TetElasticAssembly, TetElasticError, TetMaterialField, hash_usize,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_sparse::{Coo, Csr};

/// One solved equilibrium and its immutable reference vibration/stability pencil.
///
/// Tensile stress stiffens and compression softens the relevant motions.
/// The tangent is symmetric, but need not be positive definite. Downstream
/// spectral consumers must retain signed eigenvalues; nonpositive eigenvalues
/// cannot be converted into real vibration frequencies.
#[derive(Clone, Debug)]
pub struct QuadraticPrestressedAssembly {
    pencil: TetElasticAssembly,
    elastic_stiffness: Csr,
    geometric_stiffness: Csr,
    response: QuadraticThermalResponse,
    identity: ContentHash,
}

impl QuadraticPrestressedAssembly {
    /// Existing (K_e + K_g, M, free_dofs) carrier for the shared modal solver.
    #[must_use]
    pub const fn pencil(&self) -> &TetElasticAssembly {
        &self.pencil
    }

    /// Original elastic stiffness without initial-stress correction [N/m].
    #[must_use]
    pub const fn elastic_stiffness(&self) -> &Csr {
        &self.elastic_stiffness
    }

    /// Signed initial-stress contribution in the same free-DOF order [N/m].
    #[must_use]
    pub const fn geometric_stiffness(&self) -> &Csr {
        &self.geometric_stiffness
    }

    /// Actual equilibrium used at every quadrature point of the tangent.
    #[must_use]
    pub const fn response(&self) -> &QuadraticThermalResponse {
        &self.response
    }

    /// Binds the equilibrium, integration rule and complete signed pencil.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Private entry point: each public problem solves and assembles on the same
/// borrowed geometry/material/thermal inputs before handing ownership here.
/// A caller cannot inject an unrelated response or a hand-authored stress.
#[allow(clippy::too_many_arguments)]
pub(in crate::linear3) fn assemble<const N: usize, const D: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    thermal: QuadraticThermalStrainField<'_>,
    mut pencil: TetElasticAssembly,
    response: QuadraticThermalResponse,
    stiffness_rule: &[(R, f64)],
    evaluate: impl Fn([[f64; 3]; N], R, usize, usize) -> Result<Point<N>, TetElasticError>,
    maximum_absolute_mandel_strain: f64,
    cx: &Cx<'_>,
) -> Result<QuadraticPrestressedAssembly, TetElasticError> {
    checkpoint(cx)?;
    let n = pencil.free_dofs.len();
    if n == 0 {
        return Err(TetElasticError::NoFreeDofs);
    }
    let full_dofs = nodes_m
        .len()
        .checked_mul(3)
        .ok_or(TetElasticError::SizeOverflow)?;
    let mut reduced = vec![None; full_dofs];
    for (row, &dof) in pencil.free_dofs.iter().enumerate() {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        reduced[dof] = Some(row);
    }
    let mut geometric = Coo::new(n, n);
    for (element, connectivity) in elements.iter().enumerate() {
        checkpoint(cx)?;
        let nodes = connectivity.map(|node| nodes_m[node]);
        let material = materials.at(element);
        let mut local = [[0.0; N]; N];
        for (point, &(reference, weight)) in stiffness_rule.iter().enumerate() {
            checkpoint(cx)?;
            let p = evaluate(nodes, reference, element, point)?;
            let b = mandel_b::<N, D>(p.gradients);
            let (_, stress) = elastic_response_at(
                &b,
                connectivity,
                response.displacement_m(),
                &material.stiffness_mandel_pa,
                &thermal.at(element, point).free_strain_mandel,
                maximum_absolute_mandel_strain,
            )?;
            let h = core::f64::consts::FRAC_1_SQRT_2;
            let sigma = [
                [stress[0], h * stress[3], h * stress[5]],
                [h * stress[3], stress[1], h * stress[4]],
                [h * stress[5], h * stress[4], stress[2]],
            ];
            let dv = weight * p.determinant;
            if stress.iter().any(|value| !value.is_finite()) || !dv.is_finite() || dv <= 0.0 {
                return Err(invalid(
                    "quadratic initial-stress integration is nonfinite or nonpositive",
                ));
            }
            for a in 0..N {
                for b in a..N {
                    let mut contraction = 0.0;
                    for i in 0..3 {
                        for j in 0..3 {
                            contraction += p.gradients[a][i] * sigma[i][j] * p.gradients[b][j];
                        }
                    }
                    local[a][b] += dv * contraction;
                    if !local[a][b].is_finite() {
                        return Err(invalid(
                            "quadratic initial-stress element stiffness is nonfinite",
                        ));
                    }
                }
            }
        }
        for a in 0..N {
            for b in a..N {
                let value = local[a][b];
                if value == 0.0 {
                    continue;
                }
                for axis in 0..3 {
                    if let (Some(row), Some(col)) = (
                        reduced[3 * connectivity[a] + axis],
                        reduced[3 * connectivity[b] + axis],
                    ) {
                        geometric.push(row, col, value);
                        if a != b {
                            geometric.push(col, row, value);
                        }
                    }
                }
            }
        }
    }
    checkpoint(cx)?;
    let geometric_stiffness = geometric.assemble();
    checkpoint(cx)?;
    let elastic_stiffness = pencil.stiffness.clone();
    let mut tangent = Coo::new(n, n);
    for row in 0..n {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        for matrix in [&elastic_stiffness, &geometric_stiffness] {
            let (columns, values) = matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                if !value.is_finite() {
                    return Err(invalid("quadratic initial-stress pencil is nonfinite"));
                }
                tangent.push(row, column, value);
            }
        }
    }
    checkpoint(cx)?;
    pencil.stiffness = tangent.assemble();
    let mut identity =
        DomainHasher::new("org.frankensim.fs-solid.quadratic-initial-stress-pencil.v1");
    identity.update(response.identity().as_bytes());
    hash_usize(&mut identity, n)?;
    for row in 0..n {
        if row % 256 == 0 {
            checkpoint(cx)?;
        }
        hash_usize(&mut identity, pencil.free_dofs[row])?;
        for matrix in [&pencil.stiffness, &pencil.mass, &geometric_stiffness] {
            let (columns, values) = matrix.row(row);
            hash_usize(&mut identity, columns.len())?;
            for (&column, &value) in columns.iter().zip(values) {
                if !value.is_finite() {
                    return Err(invalid(
                        "quadratic assembled initial-stress pencil is unrepresentable",
                    ));
                }
                hash_usize(&mut identity, column)?;
                identity.update(&value.to_bits().to_le_bytes());
            }
        }
    }
    checkpoint(cx)?;
    Ok(QuadraticPrestressedAssembly {
        pencil,
        elastic_stiffness,
        geometric_stiffness,
        response,
        identity: identity.finalize(),
    })
}

fn invalid(what: &'static str) -> TetElasticError {
    TetElasticError::InvalidStressObservation { what }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TetElasticError> {
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)
}
