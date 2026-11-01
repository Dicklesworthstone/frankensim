//! Native anisotropic adaptive modal regression fixtures.
use super::*;

use fs_ascent::projected_al::ProjectedAlOptions;
use fs_cutfem::elastic3::ElasticityError3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_topopt::sdf3::modal::design::Sdf3EigenfrequencyDesignOptions3;
use fs_topopt::sdf3_goal::GoalRefinementError3;

fn build(
    tree: &Octree3,
    density: f64,
    checkpoint: &mut dyn FnMut() -> ControlFlow<()>,
) -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
    let material = material(density, 1.0);
    let mut callback = |_| checkpoint();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 { depth: 1, ..Default::default() },
        &mut callback,
    )
    .map_err(ElasticityError3::from)?;
    let operator = AdaptiveElasticity3::build_with_material_field(
        domain(),
        tree,
        &Top,
        ElasticMaterialField3::Uniform(&material),
        &|p| p[0] == 0.0 || p[1] == 0.0 || p[2] == 0.0,
        Default::default(),
        &mut quadrature,
    )?;
    Ok(AdaptiveSolveSpace3::jacobi(operator, 10_000_000))
}

fn frequency() -> Sdf3EigenfrequencyOptions3 {
    Sdf3EigenfrequencyOptions3 {
        solver: MatrixFreeEigenOptions {
            count: 2,
            oversampling: 6,
            max_iterations: 192,
            relative_tolerance: 1e-8,
            linear_tolerance: 1e-11,
            seed: 17,
        },
        aggregation_beta: 30.0,
    }
}

fn design_settings(reference: f64) -> Sdf3EigenfrequencyDesignOptions3 {
    Sdf3EigenfrequencyDesignOptions3 {
        frequency: frequency(),
        eigenvalue_reference: reference,
        density_floor: 0.3,
        optimizer: ProjectedAlOptions {
            max_spectral_step: 0.01,
            ..Default::default()
        },
        ..Default::default()
    }
}

mod assessment;

mod refit;

mod campaign;
