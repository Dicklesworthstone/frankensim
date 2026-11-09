//! Shared physical point response for observations and initial-stress stiffness.

use fs_material::tensor::{contract, from_mandel, to_mandel};

use super::{material_domain, material_stress, project_free};
use crate::elastic3::strain_domain::FrobeniusStrainBall;
use crate::elastic3::{Cell3, ElasticityError3, q1};

/// All returned arrays use physical tensor shear, not Mandel coordinates.
pub(super) struct PointResponse3 {
    pub(super) strain: [f64; 6],
    pub(super) elastic_strain: [f64; 6],
    pub(super) constitutive_elastic_strain: [f64; 6],
    pub(super) reference_stress: [f64; 6],
    pub(super) stress: [f64; 6],
    pub(super) energy_density: f64,
}

pub(super) fn physical_response(
    cell: &Cell3,
    local: &[f64; 24],
    position: [f64; 3],
    free_strain: [f64; 6],
    mean_free_trace: f64,
    scale: f64,
    maximum_kinematic_strain: f64,
) -> Result<PointResponse3, ElasticityError3> {
    let (_, gradients) = q1(cell.bounds, position);
    let strain = crate::elastic3::stress::mandel_strain(&gradients, local);
    let assumed = cell.project_mandel_strain(strain, local);
    for total in [strain, assumed] {
        if FrobeniusStrainBall::new(total, maximum_kinematic_strain).is_none() {
            return Err(ElasticityError3::Invalid(
                "total strain exceeds small-deformation budget; finite-deformation model required",
            ));
        }
    }
    let eigen = to_mandel(&free_strain);
    let projected_free = project_free(cell, eigen, mean_free_trace);
    let elastic = std::array::from_fn(|i| strain[i] - eigen[i]);
    let constitutive = std::array::from_fn(|i| assumed[i] - projected_free[i]);
    material_domain(cell, elastic)?;
    material_domain(cell, constitutive)?;
    // Subtract strains BEFORE the material action: nearly compatible thermal
    // expansion must not subtract two large independently rounded stresses.
    let reference_stress = material_stress(cell, constitutive);
    let stress = reference_stress.map(|s| scale * s);
    let constitutive_elastic_strain = from_mandel(&constitutive);
    let energy_density = 0.5 * contract(&constitutive_elastic_strain, &stress);
    if reference_stress
        .iter()
        .chain(&stress)
        .any(|v| !v.is_finite())
        || !energy_density.is_finite()
        || energy_density < 0.0
    {
        return Err(ElasticityError3::Invalid(
            "unrepresentable eigenstrain response",
        ));
    }
    Ok(PointResponse3 {
        strain: from_mandel(&strain),
        elastic_strain: from_mandel(&elastic),
        constitutive_elastic_strain,
        reference_stress,
        stress,
        energy_density,
    })
}
