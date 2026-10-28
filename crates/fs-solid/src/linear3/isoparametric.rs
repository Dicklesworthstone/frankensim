//! The TET10 integration/assembly path, shared with quadratic hexahedra.
//! Local upper triangles, Mandel arithmetic and COO insertion retain the
//! original TET10 operation order. Geometry and reference rules are inputs.
use super::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetMaterialField, TetThermalLoad, free_dof_map, hash_usize,
};
use fs_blake3::DomainHasher;
use fs_exec::Cx;
use fs_sparse::Coo;

pub(super) mod thermal_field;

pub(super) struct Point<const N: usize> {
    pub shapes: [f64; N],
    pub gradients: [[f64; 3]; N],
    pub determinant: f64,
    pub scaled_jacobian: f64,
}

pub(super) fn assemble<const N: usize, const D: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    material_samples: Option<QuadraticThermalStrainField<'_>>,
    fixed_dofs: &[usize],
    budget: TetAssemblyBudget,
    stiffness_rule: &[(R, f64)],
    mass_rule: &[(R, f64)],
    evaluate: impl Fn([[f64; 3]; N], R, usize, usize) -> Result<Point<N>, TetElasticError>,
    cx: &Cx<'_>,
) -> Result<TetElasticAssembly, TetElasticError> {
    debug_assert_eq!(D, 3 * N);
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    if let Some(field) = material_samples {
        field.validate(elements.len(), stiffness_rule.len())?;
    }
    let (free_dofs, reduced) = free_dof_map(nodes_m.len(), fixed_dofs, budget.maximum_free_dofs)?;
    let n = free_dofs.len();
    let mut stiffness = Coo::new(n, n);
    let mut mass = Coo::new(n, n);
    let mut total_mass_kg = 0.0;
    let mut minimum_scaled_jacobian = f64::INFINITY;
    let mut element_volumes_m3 = Vec::with_capacity(elements.len());
    for (element, connectivity) in elements.iter().enumerate() {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let nodes = connectivity.map(|i| nodes_m[i]);
        let material = materials.at(element);
        let mut ke = [[0.0; D]; D];
        let mut me = [[0.0; D]; D];
        let mut volume = 0.0;
        for (point, &(reference, weight)) in stiffness_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let p = evaluate(nodes, reference, element, point)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(p.scaled_jacobian);
            let stiffness_material = match &material_samples {
                Some(field) => field.material_at(material, element, point)?,
                None => material,
            };
            let b = mandel_b::<N, D>(p.gradients);
            let mut cb = [[0.0; D]; 6];
            for i in 0..6 {
                for j in 0..D {
                    for k in 0..6 {
                        cb[i][j] += stiffness_material.stiffness_mandel_pa[i][k] * b[k][j];
                    }
                }
            }
            let dv = weight * p.determinant;
            for i in 0..D {
                for j in i..D {
                    let entry: f64 = (0..6).map(|k| b[k][i] * cb[k][j]).sum();
                    ke[i][j] += dv * entry;
                }
            }
        }
        for (point, &(reference, weight)) in mass_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let p = evaluate(nodes, reference, element, stiffness_rule.len() + point)?;
            minimum_scaled_jacobian = minimum_scaled_jacobian.min(p.scaled_jacobian);
            let dv = weight * p.determinant;
            volume += dv;
            for a in 0..N {
                for b in a..N {
                    let v = material.density_kg_m3 * dv * p.shapes[a] * p.shapes[b];
                    for c in 0..3 {
                        me[3 * a + c][3 * b + c] += v;
                    }
                }
            }
        }
        for i in 0..D {
            for j in i..D {
                if !ke[i][j].is_finite() || !me[i][j].is_finite() {
                    return Err(TetElasticError::InvalidMaterial {
                        what: "nonfinite quadratic element matrix",
                    });
                }
                ke[j][i] = ke[i][j];
                me[j][i] = me[i][j];
            }
        }
        for i in 0..D {
            let Some(row) = reduced[3 * connectivity[i / 3] + i % 3] else {
                continue;
            };
            for j in 0..D {
                let Some(col) = reduced[3 * connectivity[j / 3] + j % 3] else {
                    continue;
                };
                if ke[i][j] != 0.0 {
                    stiffness.push(row, col, ke[i][j]);
                }
                if me[i][j] != 0.0 {
                    mass.push(row, col, me[i][j]);
                }
            }
        }
        element_volumes_m3.push(volume);
        total_mass_kg += material.density_kg_m3 * volume;
    }
    if !(total_mass_kg.is_finite() && total_mass_kg > 0.0) {
        return Err(TetElasticError::InvalidMaterial {
            what: "nonfinite quadratic total mass",
        });
    }
    let stiffness = stiffness.assemble();
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    let mass = mass.assemble();
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    if (0..n).any(|row| {
        stiffness
            .row(row)
            .1
            .iter()
            .chain(mass.row(row).1)
            .any(|v| !v.is_finite())
    }) {
        return Err(TetElasticError::InvalidMaterial {
            what: "nonfinite assembled quadratic matrix",
        });
    }
    Ok(TetElasticAssembly {
        stiffness,
        mass,
        free_dofs,
        total_mass_kg,
        minimum_scaled_jacobian,
        element_volumes_m3,
    })
}
pub(super) fn mandel_b<const N: usize, const D: usize>(gradients: [[f64; 3]; N]) -> [[f64; D]; 6] {
    let mut b = [[0.0; D]; 6];
    for (a, g) in gradients.iter().enumerate() {
        let c = 3 * a;
        let s = core::f64::consts::FRAC_1_SQRT_2;
        b[0][c] = g[0];
        b[1][c + 1] = g[1];
        b[2][c + 2] = g[2];
        b[3][c] = s * g[1];
        b[3][c + 1] = s * g[0];
        b[4][c + 1] = s * g[2];
        b[4][c + 2] = s * g[1];
        b[5][c] = s * g[2];
        b[5][c + 2] = s * g[0];
    }
    b
}

/// The same reference geometry and Mandel weak form as elastic assembly:
/// `f_th = integral B^T C epsilon_free dV`, hence `K u = f_ext + f_th`.
/// The reference mass is constant per element. Local stiffness and integrated
/// eigenstrain may both vary over the SAME stiffness quadrature. Full rows
/// include constrained DOFs so downstream support reactions are not discarded.
pub(super) fn thermal_load<const N: usize, const D: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    thermal_strains: QuadraticThermalStrainField<'_>,
    fixed_dofs: &[usize],
    budget: TetAssemblyBudget,
    stiffness_rule: &[(R, f64)],
    mass_rule: &[(R, f64)],
    evaluate: impl Fn([[f64; 3]; N], R, usize, usize) -> Result<Point<N>, TetElasticError>,
    cx: &Cx<'_>,
) -> Result<TetThermalLoad, TetElasticError> {
    debug_assert_eq!(D, 3 * N);
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    thermal_strains.validate(elements.len(), stiffness_rule.len())?;
    let (free_dofs, _) = free_dof_map(nodes_m.len(), fixed_dofs, budget.maximum_free_dofs)?;
    let full_dofs = nodes_m
        .len()
        .checked_mul(3)
        .ok_or(TetElasticError::SizeOverflow)?;
    let mut full_equivalent_force_n = vec![0.0; full_dofs];
    let mut identity = DomainHasher::new("org.frankensim.fs-solid.isoparametric-thermal-load.v1");
    for count in [
        N,
        nodes_m.len(),
        elements.len(),
        free_dofs.len(),
        stiffness_rule.len(),
        mass_rule.len(),
    ] {
        hash_usize(&mut identity, count)?;
    }
    for (index, coordinates) in nodes_m.iter().enumerate() {
        if index % 256 == 0 {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        }
        for coordinate in coordinates {
            identity.update(&coordinate.to_bits().to_le_bytes());
        }
    }
    for (index, &dof) in free_dofs.iter().enumerate() {
        if index % 256 == 0 {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        }
        hash_usize(&mut identity, dof)?;
    }
    for (element, connectivity) in elements.iter().enumerate() {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let reference_material = materials.at(element);
        let material = thermal_strains.material_at(reference_material, element, 0)?;
        let thermal = thermal_strains.at(element, 0);
        if thermal.elastic_material_state_identity != material.material_state_identity {
            return Err(TetElasticError::ThermalElasticStateMismatch { element });
        }
        for &node in connectivity {
            hash_usize(&mut identity, node)?;
        }
        identity.update(&material.material_state_identity);
        identity.update(thermal.identity.as_bytes());
        if !thermal_strains.is_elementwise() {
            // Distinguish pointwise input while retaining the original identity
            // and arithmetic for all element-constant callers.
            identity.update(b"quadrature-thermal-strain-v1");
        }
        for coefficient in material.stiffness_mandel_pa.iter().flatten() {
            identity.update(&coefficient.to_bits().to_le_bytes());
        }
        let mut free_stress = [0.0; 6];
        for (row, stress) in free_stress.iter_mut().enumerate() {
            for column in 0..6 {
                *stress +=
                    material.stiffness_mandel_pa[row][column] * thermal.free_strain_mandel[column];
            }
        }
        if free_stress.iter().any(|value| !value.is_finite()) {
            return Err(TetElasticError::NonFiniteThermalLoad {
                element,
                local_dof: None,
            });
        }
        let nodes = connectivity.map(|node| nodes_m[node]);
        let mut local_force = [0.0; D];
        for (point, &(reference, weight)) in stiffness_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            if !thermal_strains.is_elementwise() {
                let sample = thermal_strains.at(element, point);
                let material = thermal_strains.material_at(reference_material, element, point)?;
                if sample.elastic_material_state_identity != material.material_state_identity {
                    return Err(TetElasticError::ThermalElasticStateMismatch { element });
                }
                identity.update(sample.identity.as_bytes());
                if thermal_strains.has_local_materials() {
                    identity.update(b"quadrature-thermomechanical-material-v1");
                    identity.update(&material.material_state_identity);
                    for coefficient in material.stiffness_mandel_pa.iter().flatten() {
                        identity.update(&coefficient.to_bits().to_le_bytes());
                    }
                }
                free_stress.fill(0.0);
                for (row, stress) in free_stress.iter_mut().enumerate() {
                    for column in 0..6 {
                        *stress += material.stiffness_mandel_pa[row][column]
                            * sample.free_strain_mandel[column];
                    }
                }
                if free_stress.iter().any(|value| !value.is_finite()) {
                    return Err(TetElasticError::NonFiniteThermalLoad {
                        element,
                        local_dof: None,
                    });
                }
            }
            let p = evaluate(nodes, reference, element, point)?;
            let b = mandel_b::<N, D>(p.gradients);
            let dv = weight * p.determinant;
            identity.update(&weight.to_bits().to_le_bytes());
            identity.update(&p.determinant.to_bits().to_le_bytes());
            for gradient in p.gradients.iter().flatten() {
                identity.update(&gradient.to_bits().to_le_bytes());
            }
            for (dof, force) in local_force.iter_mut().enumerate() {
                let density: f64 = (0..6)
                    .map(|component| b[component][dof] * free_stress[component])
                    .sum();
                *force += dv * density;
            }
        }
        // Preserve the elastic owner's geometry-admission domain even where
        // its mass rule samples points absent from the stiffness/load rule.
        for (point, &(reference, _)) in mass_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            evaluate(nodes, reference, element, stiffness_rule.len() + point)?;
        }
        for (local_dof, &force) in local_force.iter().enumerate() {
            let full_dof = 3 * connectivity[local_dof / 3] + local_dof % 3;
            full_equivalent_force_n[full_dof] += force;
            if !force.is_finite() || !full_equivalent_force_n[full_dof].is_finite() {
                return Err(TetElasticError::NonFiniteThermalLoad {
                    element,
                    local_dof: Some(local_dof),
                });
            }
        }
    }
    let mut reduced_force_n = Vec::with_capacity(free_dofs.len());
    for (index, &dof) in free_dofs.iter().enumerate() {
        if index % 256 == 0 {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        }
        reduced_force_n.push(full_equivalent_force_n[dof]);
    }
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    Ok(TetThermalLoad {
        full_equivalent_force_n,
        reduced_force_n,
        free_dofs,
        identity: identity.finalize(),
    })
}

// Quadratic thermal equilibrium and integrated stress/energy/support response.
pub(super) mod thermal_response;

// Public P1 support-motion methods reuse the private linear3 assembly/solve owners.
mod prescribed;

// Full-row P1 dynamic response reuses the same canonical element kernels.
mod response;

// Material-scale derivatives reuse the canonical P1 tiles.
mod material_gradient;

// Stress observations transpose the same P1 B and constitutive tensors.
mod stress_gradient;
