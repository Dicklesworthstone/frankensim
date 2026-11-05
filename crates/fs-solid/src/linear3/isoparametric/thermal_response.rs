//! Isoparametric thermal equilibrium and work-conjugate response on one frozen mesh.
//! Numerical iteration remains in fs-solver; geometry, integration and material
//! admission remain in the existing isoparametric owners.
use super::super::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticError, TetMaterialField,
    TetStaticSolveConfig, add_nodal_forces_to_thermal_load, hash_usize, relative_mandel_force,
    relative_mandel_strain, validate_rigid_mode_constraints,
};
use super::{
    Point, StrainProjection, assemble_with_formulation, mandel_b, thermal_load_with_formulation,
};
use crate::linear::{Formulation, Jacobi};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_solver::checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg};
use fs_solver::op::CsrOp;

/// Actual small-strain Cauchy stress at one stiffness quadrature point.
///
/// Points remain in the native integration-rule order. They are discrete
/// samples, not nodal extrapolations or a continuum peak-stress bound.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticStressPoint {
    stress_mandel_pa: [f64; 6],
    material_state_identity: ContentHash,
    integration_weight_m3: f64,
}

impl QuadraticStressPoint {
    /// Mesh-frame Mandel [xx, yy, zz, xy, yz, zx], with sqrt(2) shear entries.
    #[must_use]
    pub const fn stress_mandel_pa(&self) -> &[f64; 6] {
        &self.stress_mandel_pa
    }

    /// Actual local elastic source, including temperature-dependent updates.
    #[must_use]
    pub const fn material_state_identity(&self) -> ContentHash {
        self.material_state_identity
    }

    /// Reference-geometry integration weight, including the Jacobian [m^3].
    #[must_use]
    pub const fn integration_weight_m3(&self) -> f64 {
        self.integration_weight_m3
    }
}

/// Integrated response of one isoparametric element in mesh-frame Mandel coordinates.
///
/// Mean strain and stress are volume-weighted averages over the stiffness rule.
/// These means are neither pointwise stresses nor peak-stress estimates. Energy
/// integrates the actual compatible or B-bar assumed field before averaging,
/// so bending energy is retained even when mean strain vanishes. All values
/// use the undeformed reference geometry and the parent response's formulation.
/// The original point stresses are retained separately for strength assessments;
/// averaging opposite bending stresses must not erase a tensile hotspot.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadraticElementThermalResponse {
    mean_elastic_strain_mandel: [f64; 6],
    mean_stress_mandel_pa: [f64; 6],
    elastic_energy_j: f64,
    quadrature_volume_m3: f64,
    stress_points: Vec<QuadraticStressPoint>,
}

impl QuadraticElementThermalResponse {
    /// Volume average of the selected strain field minus thermal eigenstrain.
    /// B-bar reports its assumed strain, not the raw compatible strain.
    #[must_use]
    pub const fn mean_elastic_strain_mandel(&self) -> &[f64; 6] {
        &self.mean_elastic_strain_mandel
    }

    /// Volume average of C(Bu - epsilon_free) [Pa], with the same compatible
    /// or B-bar columns as stiffness, source load and internal force.
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

    /// Original stress samples in the stiffness quadrature's native order.
    /// B-bar samples use the same assumed strain as stiffness and forces.
    #[must_use]
    pub fn stress_points(&self) -> &[QuadraticStressPoint] {
        &self.stress_points
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
    formulation: Formulation,
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
    /// Strain formulation used consistently in stiffness, loads and recovery.
    #[must_use]
    pub const fn formulation(&self) -> Formulation {
        self.formulation
    }

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
    solve_with_formulation::<N, D, R>(
        nodes_m,
        elements,
        materials,
        thermal,
        nodal_forces_n,
        fixed_dofs,
        budget,
        Formulation::Standard,
        stiffness_rule,
        mass_rule,
        evaluate,
        config,
        mesh_frame,
        maximum_absolute_mandel_strain,
        cx,
    )
}

#[allow(clippy::too_many_arguments)]
pub(in crate::linear3) fn solve_with_formulation<const N: usize, const D: usize, R: Copy>(
    nodes_m: &[[f64; 3]],
    elements: &[[usize; N]],
    materials: TetMaterialField<'_>,
    thermal: QuadraticThermalStrainField<'_>,
    nodal_forces_n: Option<&[(usize, f64)]>,
    fixed_dofs: &[usize],
    budget: TetAssemblyBudget,
    formulation: Formulation,
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
    let load = thermal_load_with_formulation::<N, D, R>(
        nodes_m,
        elements,
        materials,
        thermal,
        fixed_dofs,
        budget,
        formulation,
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
        let assembly = assemble_with_formulation::<N, D, R>(
            nodes_m,
            elements,
            materials,
            Some(thermal),
            fixed_dofs,
            budget,
            formulation,
            stiffness_rule,
            mass_rule,
            &evaluate,
            |_, _| None,
            cx,
        )?;
        let preconditioner = Jacobi::new(&assembly.stiffness);
        let operator = CsrOp::symmetric(assembly.stiffness);
        let solved = checked_cg(
            &operator,
            &preconditioner,
            &load.reduced_force_n,
            CheckedCgConfig {
                tolerance: config.relative_residual_tolerance,
                max_iterations: config.maximum_iterations,
                max_corrections: 2,
            },
            |_| cx.checkpoint().map_err(|_| TetElasticError::Cancelled),
        )
        .map_err(|error| match error {
            CheckedCgError::InvalidInput(what) => {
                TetElasticError::InvalidStaticSolveConfig { what }
            }
            CheckedCgError::Interrupted(error) => error,
        })?;
        let relative = solved.report.rel_residual;
        if !relative.is_finite() || relative > config.relative_residual_tolerance {
            return Err(TetElasticError::StaticSolveFailed {
                iterations: solved.report.iters,
                true_relative_residual: relative,
            });
        }
        (solved.x, solved.report.iters, relative)
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
    let mut identity = DomainHasher::new("org.frankensim.fs-solid.quadratic-thermal-response.v2");
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
        let projection =
            StrainProjection::prepare(nodes, element, formulation, stiffness_rule, &evaluate, cx)?;
        identity.update(&material.density_kg_m3.to_bits().to_le_bytes());
        let mut mean_elastic_strain_mandel = [0.0; 6];
        let mut mean_stress_mandel_pa = [0.0; 6];
        let mut volume = 0.0;
        let mut energy = 0.0;
        let mut stress_points = Vec::with_capacity(stiffness_rule.len());
        for (point, &(reference, weight)) in stiffness_rule.iter().enumerate() {
            cx.checkpoint().map_err(|_| TetElasticError::Cancelled)?;
            let compatible_point;
            let p = if let Some(prepared) = projection.point(point) {
                prepared
            } else {
                compatible_point = evaluate(nodes, reference, element, point)?;
                &compatible_point
            };
            let material = thermal.material_at(material, element, point)?;
            let free_strain = &thermal.at(element, point).free_strain_mandel;
            if formulation == Formulation::BBar {
                // A projection must not hide an inadmissible compatible
                // displacement field. Bound it independently of assumed strain.
                elastic_strain_at(
                    &mandel_b::<N, D>(p.gradients),
                    connectivity,
                    &displacement_m,
                    free_strain,
                    maximum_absolute_mandel_strain,
                )?;
            }
            let b = projection.columns::<D>(p);
            let (elastic, stress) = elastic_response_at(
                &b,
                connectivity,
                &displacement_m,
                &material.stiffness_mandel_pa,
                free_strain,
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
            identity.update(&material.material_state_identity);
            for value in stress.iter().chain(core::iter::once(&dv)) {
                identity.update(&value.to_bits().to_le_bytes());
            }
            stress_points.push(QuadraticStressPoint {
                stress_mandel_pa: stress,
                material_state_identity: ContentHash(material.material_state_identity),
                integration_weight_m3: dv,
            });
            for i in 0..6 {
                mean_elastic_strain_mandel[i] += dv * elastic[i];
                mean_stress_mandel_pa[i] += dv * stress[i];
            }
            let force = relative_mandel_force(&b, &stress, dv);
            for d in 0..D {
                let accumulated = &mut nodal_internal_force_n[connectivity[d / 3]][d % 3];
                *accumulated += force[d];
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
            stress_points,
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
        formulation,
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
    let elastic = elastic_strain_at(
        b,
        connectivity,
        displacement_m,
        free_strain_mandel,
        maximum_absolute_mandel_strain,
    )?;
    let stress: [f64; 6] = core::array::from_fn(|i| {
        (0..6).fold(0.0, |sum, j| {
            stiffness_mandel_pa[i][j].mul_add(elastic[j], sum)
        })
    });
    Ok((elastic, stress))
}

fn elastic_strain_at<const N: usize, const D: usize>(
    b: &[[f64; D]; 6],
    connectivity: &[usize; N],
    displacement_m: &[[f64; 3]],
    free_strain_mandel: &[f64; 6],
    maximum_absolute_mandel_strain: f64,
) -> Result<[f64; 6], TetElasticError>{
    debug_assert_eq!(D, 3 * N);
    let total = relative_mandel_strain(
        N,
        |d| core::array::from_fn(|i| b[i][d]),
        |d| displacement_m[connectivity[d / 3]][d % 3],
    );
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
    Ok(elastic)
}

#[cfg(test)]
mod strain_tests {
    use super::*;
    use crate::linear3::{hex, tet10};

    fn geometry([x, y, z]: [f64; 3]) -> [f64; 3] {
        [
            0.25 + 1.25 * x + 0.125 * y,
            -0.5 + 0.25 * x + 0.875 * y + 0.125 * z,
            0.75 - 0.125 * x + 1.5 * z,
        ]
    }

    fn check<const N: usize, const D: usize>(nodes: [[f64; 3]; N], gradients: [[f64; 3]; N]) {
        let b = mandel_b::<N, D>(gradients);
        let connectivity: [usize; N] = core::array::from_fn(|i| i);
        // Exactly representable nodal data keeps the test independent of input
        // quantization when a common, much larger displacement is added.
        let u = nodes.map(|[x, y, z]| {
            [
                x / 1024.0 + y / 2048.0 - z / 4096.0,
                -x / 4096.0 + y / 2048.0 + z / 1024.0,
                x / 2048.0 - y / 4096.0 + z / 2048.0,
            ]
        });
        let shift = [2.0_f64.powi(35), -2.0_f64.powi(36), 2.0_f64.powi(34)];
        let translated = u.map(|x| core::array::from_fn(|c| x[c] + shift[c]));
        let free = [1.0, -0.5, 0.25, -0.75, 0.5, 0.25].map(|x| x / 8192.0);
        let c = core::array::from_fn(|i| {
            core::array::from_fn(|j| if i == j { 30.0 * (i + 1) as f64 } else { 0.0 })
        });
        let (elastic, stress) =
            elastic_response_at(&b, &connectivity, &u, &c, &free, 0.01).unwrap();
        assert_eq!(
            (elastic, stress),
            elastic_response_at(&b, &connectivity, &translated, &c, &free, 0.01).unwrap()
        );
        assert_eq!(
            elastic_strain_at(&b, &connectivity, &[shift; N], &[0.0; 6], 1e-20).unwrap(),
            [0.0; 6]
        );
        let s = core::f64::consts::FRAC_1_SQRT_2;
        let total = [
            1.0 / 1024.0,
            1.0 / 2048.0,
            1.0 / 2048.0,
            s / 4096.0,
            3.0 * s / 4096.0,
            s / 4096.0,
        ];
        for i in 0..6 {
            assert!((elastic[i] - (total[i] - free[i])).abs() < 2e-16);
        }
        let force = relative_mandel_force(&b, &stress, 1.7);
        let work: f64 = force
            .iter()
            .enumerate()
            .map(|(d, f)| f * u[d / 3][d % 3])
            .sum();
        let exact_work: f64 = total.iter().zip(stress).map(|(e, s)| 1.7 * e * s).sum();
        assert!((work - exact_work).abs() < 2e-16);
        for axis in 0..3 {
            assert!(force.iter().skip(axis).step_by(3).sum::<f64>().abs() < 2e-15);
        }
    }

    #[test]
    fn g3_all_isoparametric_families_preserve_translated_affine_thermal_strain() {
        let vertices = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
        .map(geometry);
        let tet_nodes = core::array::from_fn(|a| {
            if a < 4 {
                vertices[a]
            } else {
                let [i, j] = tet10::TET10_EDGES[a - 4];
                core::array::from_fn(|c| 0.5 * (vertices[i][c] + vertices[j][c]))
            }
        });
        let tet_point = tet10::tet10_point(tet_nodes, [0.1, 0.2, 0.3, 0.4], 1e-12, 0, 0).unwrap();
        check::<10, 30>(tet_nodes, tet_point.gradients);
        let hex8_nodes =
            core::array::from_fn(|a| geometry(hex::HEX27_REFERENCE_NODES[a].map(f64::from)));
        let hex20_nodes = hex::HEX20_REFERENCE_NODES.map(|x| geometry(x.map(f64::from)));
        let hex27_nodes = hex::HEX27_REFERENCE_NODES.map(|x| geometry(x.map(f64::from)));
        for p in [[0.17, -0.23, 0.41], [-0.51, 0.31, -0.11]] {
            check::<8, 24>(
                hex8_nodes,
                hex::hex8_point(hex8_nodes, p, 1e-12, 0, 0)
                    .unwrap()
                    .gradients,
            );
            check::<20, 60>(
                hex20_nodes,
                hex::hex20_point(hex20_nodes, p, 1e-12, 0, 0)
                    .unwrap()
                    .gradients,
            );
            check::<27, 81>(
                hex27_nodes,
                hex::hex27_point(hex27_nodes, p, 1e-12, 0, 0)
                    .unwrap()
                    .gradients,
            );
        }
    }
}
