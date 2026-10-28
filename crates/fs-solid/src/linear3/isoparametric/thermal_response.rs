//! Quadratic thermal equilibrium and work-conjugate response on one frozen mesh.
//! Numerical iteration remains in fs-solver; geometry, integration and material
//! admission remain in the existing isoparametric owners.
use super::super::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticError, TetMaterialField,
    TetStaticSolveConfig, add_nodal_forces_to_thermal_load, euclidean_norm, hash_usize,
    validate_rigid_mode_constraints,
};
use super::{Point, assemble, mandel_b, thermal_load};
use crate::linear::Jacobi;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_solver::krylov::CgState;
use fs_solver::op::{CsrOp, LinearOp};

/// Integrated response of one quadratic element in mesh-frame Mandel coordinates.
///
/// Strain and stress are volume-weighted averages over the stiffness rule.
/// Quadratic strain varies within each element: these are neither pointwise
/// stresses nor peak-stress estimates. Energy integrates the actual quadratic
/// field before averaging, so bending energy is retained even when mean strain
/// vanishes. All values use the undeformed reference geometry.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticElementThermalResponse {
    mean_elastic_strain_mandel: [f64; 6],
    mean_stress_mandel_pa: [f64; 6],
    elastic_energy_j: f64,
    quadrature_volume_m3: f64,
}

impl QuadraticElementThermalResponse {
    /// Volume average of compatible strain minus the thermal eigenstrain.
    #[must_use]
    pub const fn mean_elastic_strain_mandel(&self) -> &[f64; 6] {
        &self.mean_elastic_strain_mandel
    }

    /// Volume average of actual Cauchy stress C(Bu - epsilon_free) [Pa].
    #[must_use]
    pub const fn mean_stress_mandel_pa(&self) -> &[f64; 6] {
        &self.mean_stress_mandel_pa
    }

    /// Integral of one half elastic strain contracted with stress [J].
    #[must_use]
    pub const fn elastic_energy_j(&self) -> f64 {
        self.elastic_energy_j
    }

    /// Reference volume integrated by the stiffness/response rule [m^3].
    #[must_use]
    pub const fn quadrature_volume_m3(&self) -> f64 {
        self.quadrature_volume_m3
    }
}

/// Accepted thermal equilibrium and its physical response, published atomically.
///
/// Fixed displacement components are zero. Internal force integrates the actual
/// stress. The nodal force balance subtracts applied dead nodal loads, giving the
/// support force ON the solid at fixed DOFs and the remaining equilibrium residual
/// at free DOFs. Residual entries are retained, not overwritten with zero.
/// Identity binds thermal and mechanical loads, actual material coefficients,
/// solve/strain controls and the full result.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticThermalResponse {
    displacement_m: Vec<[f64; 3]>,
    elements: Vec<QuadraticElementThermalResponse>,
    nodal_internal_force_n: Vec<[f64; 3]>,
    nodal_force_balance_n: Vec<[f64; 3]>,
    elastic_energy_j: f64,
    iterations: usize,
    true_relative_residual: f64,
    mesh_frame: ContentHash,
    thermal_load_identity: ContentHash,
    load_identity: ContentHash,
    identity: ContentHash,
}

impl QuadraticThermalResponse {
    /// Compatible displacement on every geometry node [m].
    #[must_use]
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }

    /// Integrated element observations in connectivity order.
    #[must_use]
    pub fn elements(&self) -> &[QuadraticElementThermalResponse] {
        &self.elements
    }

    /// Full integral B^T sigma dV, including support rows [N].
    #[must_use]
    pub fn nodal_internal_force_n(&self) -> &[[f64; 3]] {
        &self.nodal_internal_force_n
    }

    /// Internal force minus applied dead nodal force [N]. At fixed DOFs this
    /// is the support reaction ON the solid; at free DOFs it is the unmodified
    /// equilibrium residual. Prescribed displacement components remain zero.
    #[must_use]
    pub fn nodal_force_balance_n(&self) -> &[[f64; 3]] {
        &self.nodal_force_balance_n
    }

    /// Total stored elastic energy in the thermal equilibrium [J].
    #[must_use]
    pub const fn elastic_energy_j(&self) -> f64 {
        self.elastic_energy_j
    }

    /// Iterations used by the existing reduced SPD solver.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }

    /// Independently recomputed norm(f_thermal + f_external - K u) divided
    /// by the norm of the complete reduced load.
    #[must_use]
    pub const fn true_relative_residual(&self) -> f64 {
        self.true_relative_residual
    }

    /// Caller-declared frame of mesh coordinates and Mandel components.
    #[must_use]
    pub const fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }

    /// Original thermal-load and geometry/support/quadrature identity.
    /// This excludes any separately supplied dead nodal forces.
    #[must_use]
    pub const fn thermal_load_identity(&self) -> ContentHash {
        self.thermal_load_identity
    }

    /// Complete load identity, including dead nodal forces when supplied.
    #[must_use]
    pub const fn load_identity(&self) -> ContentHash {
        self.load_identity
    }

    /// Identity binding the complete admitted response.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

#[allow(clippy::too_many_arguments)]
pub(in crate::linear3) fn solve<const N: usize, const D: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    thermal: QuadraticThermalStrainField<'_>,
    nodal_forces_n: Option<&[(usize, f64)]>,
    fixed_dofs: &[usize],
    budget: TetAssemblyBudget,
    stiffness_rule: &[(R, f64)],
    mass_rule: &[(R, f64)],
    evaluate: impl Fn([[f64; 3]; N], R, usize, usize) -> Result<Point<N>, TetElasticError>,
    config: TetStaticSolveConfig,
    mesh_frame: ContentHash,
    maximum_absolute_mandel_strain: f64,
    cx: &Cx<'_>,
) -> Result<QuadraticThermalResponse, TetElasticError> {
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    if !(config.relative_residual_tolerance.is_finite()
        && config.relative_residual_tolerance > 0.0
        && config.relative_residual_tolerance < 1.0)
        || config.maximum_iterations == 0
    {
        return Err(TetElasticError::InvalidStaticSolveConfig {
            what: "relative tolerance must lie in (0,1) and iteration limit must be positive",
        });
    }
    if mesh_frame == ContentHash([0; 32])
        || !(maximum_absolute_mandel_strain.is_finite() && maximum_absolute_mandel_strain > 0.0)
    {
        return Err(TetElasticError::InvalidStressObservation {
            what: "mesh frame must be nonzero and strain ceiling must be finite and positive",
        });
    }
    let load = thermal_load::<N, D, R>(
        nodes_m,
        elements,
        materials,
        thermal,
        fixed_dofs,
        budget,
        stiffness_rule,
        mass_rule,
        &evaluate,
        cx,
    )?;
    let thermal_load_identity = load.identity;
    let load = match nodal_forces_n {
        Some(forces) => add_nodal_forces_to_thermal_load(load, nodes_m.len(), forces, cx)?,
        None => load,
    };
    let (reduced, iterations, true_relative_residual) = if load.free_dofs.is_empty() {
        // A fully restrained thermal body has no linear solve, but still has
        // physical stress, stored energy and nonzero support forces.
        (Vec::new(), 0, 0.0)
    } else {
        validate_rigid_mode_constraints(nodes_m, elements, fixed_dofs)?;
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let assembly = assemble::<N, D, R>(
            nodes_m,
            elements,
            materials,
            Some(thermal),
            fixed_dofs,
            budget,
            stiffness_rule,
            mass_rule,
            &evaluate,
            cx,
        )?;
        let preconditioner = Jacobi::new(&assembly.stiffness);
        let operator = CsrOp::symmetric(assembly.stiffness);
        let mut state = CgState::new(&operator, &preconditioner, &load.reduced_force_n);
        while state.iters < config.maximum_iterations
            && state.rel_residual() >= config.relative_residual_tolerance
        {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let remaining = config.maximum_iterations - state.iters;
            let _ = state.run(
                &operator,
                &preconditioner,
                config.relative_residual_tolerance,
                remaining.min(64),
            );
        }
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let mut applied = vec![0.0; state.x.len()];
        operator.apply(&state.x, &mut applied);
        let residual = euclidean_norm(
            load.reduced_force_n
                .iter()
                .zip(&applied)
                .map(|(force, ku)| force - ku),
        );
        let force = euclidean_norm(load.reduced_force_n.iter().copied());
        let relative = residual / force.max(f64::MIN_POSITIVE);
        if !relative.is_finite() || relative > config.relative_residual_tolerance {
            return Err(TetElasticError::StaticSolveFailed {
                iterations: state.iters,
                true_relative_residual: relative,
            });
        }
        (state.x, state.iters, relative)
    };
    let mut displacement_m = vec![[0.0; 3]; nodes_m.len()];
    for (index, (&dof, &value)) in load.free_dofs.iter().zip(&reduced).enumerate() {
        if index % 256 == 0 {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        }
        displacement_m[dof / 3][dof % 3] = value;
    }

    let mut response = Vec::with_capacity(elements.len());
    let mut nodal_internal_force_n = vec![[0.0; 3]; nodes_m.len()];
    let mut elastic_energy_j = 0.0;
    let mut identity = DomainHasher::new("org.frankensim.fs-solid.quadratic-thermal-response.v1");
    identity.update(load.identity.as_bytes());
    identity.update(mesh_frame.as_bytes());
    identity.update(&maximum_absolute_mandel_strain.to_bits().to_le_bytes());
    identity.update(&config.relative_residual_tolerance.to_bits().to_le_bytes());
    hash_usize(&mut identity, config.maximum_iterations)?;
    hash_usize(&mut identity, iterations)?;
    identity.update(&true_relative_residual.to_bits().to_le_bytes());

    for (element, connectivity) in elements.iter().enumerate() {
        cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        let nodes = connectivity.map(|node| nodes_m[node]);
        let material = materials.at(element);
        identity.update(&material.density_kg_m3.to_bits().to_le_bytes());
        let mut mean_elastic_strain_mandel = [0.0; 6];
        let mut mean_stress_mandel_pa = [0.0; 6];
        let mut volume = 0.0;
        let mut energy = 0.0;
        for (point, &(reference, weight)) in stiffness_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let p = evaluate(nodes, reference, element, point)?;
            let material = thermal.material_at(material, element, point)?;
            let b = mandel_b::<N, D>(p.gradients);
            let (elastic, stress) = elastic_response_at(
                &b,
                connectivity,
                &displacement_m,
                &material.stiffness_mandel_pa,
                &thermal.at(element, point).free_strain_mandel,
                maximum_absolute_mandel_strain,
            )?;
            let dv = weight * p.determinant;
            let point_energy =
                (0..6).fold(0.0, |sum, i| elastic[i].mul_add(stress[i], sum)) * (0.5 * dv);
            if stress.iter().any(|v| !v.is_finite())
                || !point_energy.is_finite()
                || point_energy < 0.0
                || !dv.is_finite()
                || dv <= 0.0
            {
                return Err(TetElasticError::InvalidStressObservation {
                    what: "quadratic stress, energy or integration volume is inadmissible",
                });
            }
            volume += dv;
            energy += point_energy;
            for i in 0..6 {
                mean_elastic_strain_mandel[i] += dv * elastic[i];
                mean_stress_mandel_pa[i] += dv * stress[i];
            }
            for d in 0..D {
                let force = (0..6).fold(0.0, |sum, i| b[i][d].mul_add(stress[i], sum)) * dv;
                let accumulated = &mut nodal_internal_force_n[connectivity[d / 3]][d % 3];
                *accumulated += force;
                if !accumulated.is_finite() {
                    return Err(TetElasticError::InvalidStressObservation {
                        what: "quadratic nodal internal force is nonfinite",
                    });
                }
            }
        }
        for i in 0..6 {
            mean_elastic_strain_mandel[i] /= volume;
            mean_stress_mandel_pa[i] /= volume;
        }
        elastic_energy_j += energy;
        if !volume.is_finite()
            || volume <= 0.0
            || !elastic_energy_j.is_finite()
            || mean_elastic_strain_mandel
                .iter()
                .chain(&mean_stress_mandel_pa)
                .any(|v| !v.is_finite())
        {
            return Err(TetElasticError::InvalidStressObservation {
                what: "quadratic element response is unrepresentable",
            });
        }
        for value in mean_elastic_strain_mandel
            .iter()
            .chain(&mean_stress_mandel_pa)
            .chain([&energy, &volume])
        {
            identity.update(&value.to_bits().to_le_bytes());
        }
        response.push(QuadraticElementThermalResponse {
            mean_elastic_strain_mandel,
            mean_stress_mandel_pa,
            elastic_energy_j: energy,
            quadrature_volume_m3: volume,
        });
    }
    let mut nodal_force_balance_n = nodal_internal_force_n.clone();
    if let Some(forces) = nodal_forces_n {
        for (index, &(dof, applied)) in forces.iter().enumerate() {
            if index % 256 == 0 {
                cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            }
            let balance = &mut nodal_force_balance_n[dof / 3][dof % 3];
            *balance -= applied;
            if !balance.is_finite() {
                return Err(TetElasticError::InvalidStressObservation {
                    what: "quadratic support reaction or force residual is nonfinite",
                });
            }
        }
    }
    for (node, (displacement, force)) in displacement_m
        .iter()
        .zip(&nodal_internal_force_n)
        .enumerate()
    {
        if node % 256 == 0 {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
        }
        for value in displacement.iter().chain(force) {
            identity.update(&value.to_bits().to_le_bytes());
        }
    }
    identity.update(&elastic_energy_j.to_bits().to_le_bytes());
    cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
    Ok(QuadraticThermalResponse {
        displacement_m,
        elements: response,
        nodal_internal_force_n,
        nodal_force_balance_n,
        elastic_energy_j,
        iterations,
        true_relative_residual,
        mesh_frame,
        thermal_load_identity,
        load_identity: load.identity,
        identity: identity.finalize(),
    })
}

// Initial-stress vibration uses the same admitted quadrature response.
pub(in crate::linear3) mod prestress;

/// Shared kinematics and constitutive contraction for observation and tangent.
/// Keep the operation order identical to the original response recovery.
fn elastic_response_at<const N: usize, const D: usize>(
    b: &[[f64; D]; 6],
    connectivity: &[usize; N],
    displacement_m: &[[f64; 3]],
    stiffness_mandel_pa: &[[f64; 6]; 6],
    free_strain_mandel: &[f64; 6],
    maximum_absolute_mandel_strain: f64,
) -> Result<([f64; 6], [f64; 6]), TetElasticError> {
    let total: [f64; 6] = core::array::from_fn(|i| {
        (0..D).fold(0.0, |sum, d| {
            b[i][d].mul_add(displacement_m[connectivity[d / 3]][d % 3], sum)
        })
    });
    let elastic: [f64; 6] = core::array::from_fn(|i| total[i] - free_strain_mandel[i]);
    if total
        .iter()
        .chain(&elastic)
        .any(|v| !v.is_finite() || v.abs() > maximum_absolute_mandel_strain)
    {
        return Err(TetElasticError::InvalidStressObservation {
            what: "quadrature total or elastic strain exceeds the finite small-strain ceiling",
        });
    }
    let stress: [f64; 6] = core::array::from_fn(|i| {
        (0..6).fold(0.0, |sum, j| {
            stiffness_mandel_pa[i][j].mul_add(elastic[j], sum)
        })
    });
    Ok((elastic, stress))
}
