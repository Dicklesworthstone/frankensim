//! Independent current-carrying 3D magnetic oracles and nontrivial topology.
//! The source cards are manufactured, not measured engineering materials.
mod common;
use common::{material, near};
use fs_em::{
    planar::permanent::MagneticEnergyReference,
    spatial::{
        SpatialMagneticError, SpatialMagneticMaterial, SpatialRecoilMaterial,
        vector_potential::{
            SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
            VectorPotentialFace,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::{
    hodge::hodge_diagonal_barycentric,
    weighted::WeightedAssemblyLimits,
    whitney::{ElementGeometry, try_element_geometry},
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_solver::CheckedCgConfig;
use std::collections::BTreeSet;
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 4000,
        max_dofs: 20000,
        max_triplets: 150000,
    }
}
fn controls() -> VectorMagneticControls {
    VectorMagneticControls {
        field: CheckedCgConfig {
            tolerance: 1e-11,
            max_iterations: 2500,
            max_corrections: 2,
        },
        gauge: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 2500,
            max_corrections: 2,
        },
        maximum_gauge_product_terms: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 3e-8,
    }
}
fn grid(
    nx: usize,
    ny: usize,
    nz: usize,
    spacing: f64,
    hole: bool,
) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut p = Vec::new();
    for k in 0..=nz {
        for j in 0..=ny {
            for i in 0..=nx {
                p.push([i as f64 * spacing, j as f64 * spacing, k as f64 * spacing]);
            }
        }
    }
    let id = |i, j, k| ((k * (ny + 1) + j) * (nx + 1) + i) as u32;
    let mut t = Vec::new();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                if hole && i == 1 && j == 1 {
                    continue;
                }
                let (a, b, c, d) = (
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i, j + 1, k),
                    id(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i, j + 1, k + 1),
                    id(i + 1, j + 1, k + 1),
                );
                t.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    (p, t)
}
fn cube(n: usize, mu: f64) -> SpatialVectorMagnetostatics {
    let (p, t) = grid(n, n, n, 1.0 / n as f64, false);
    let m = vec![SpatialMagneticMaterial::Linear(material(mu)); t.len()];
    SpatialVectorMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap()
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn trace(problem: &SpatialVectorMagnetostatics, a: &[f64]) -> VectorMagneticBoundary {
    VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: problem.face_edges(face).unwrap().map(|e| a[e]),
            })
            .collect(),
    }
}
fn uniform_trace(problem: &SpatialVectorMagnetostatics, b: [f64; 3], gauge: bool) -> Vec<f64> {
    let psi = |v: usize| {
        if gauge {
            ((v * 7) % 11) as f64 * 0.17
        } else {
            0.0
        }
    };
    problem
        .mesh()
        .edges
        .iter()
        .map(|&[u, v]| {
            let (p, q) = (
                problem.positions_m()[u as usize],
                problem.positions_m()[v as usize],
            );
            let mid = std::array::from_fn(|j| 0.5 * (p[j] + q[j]));
            let delta = std::array::from_fn(|j| q[j] - p[j]);
            0.5 * dot(cross(b, mid), delta) + psi(v as usize) - psi(u as usize)
        })
        .collect()
}
fn assert_physical(field: &fs_em::spatial::vector_potential::SpatialVectorField<'_>) {
    assert!(field.maximum_residual_a < controls().current_tolerance_a);
    assert!(field.source_compatibility_defect_a < controls().current_tolerance_a);
    assert!(field.gauge_flux_defect_wb < controls().flux_tolerance_wb);
    assert!(field.maximum_divergence_flux_wb < controls().flux_tolerance_wb);
    assert!(field.maximum_recovered_flux_defect_wb < controls().flux_tolerance_wb);
    assert!(field.work_residual_j.abs() < controls().energy_tolerance_j);
    for report in [
        &field.field_report,
        &field.source_report,
        &field.gauge_report,
    ]
    .into_iter()
    .flatten()
    {
        assert!(report.converged_euclidean());
    }
}
#[test]
fn exact_uniform_field_and_discrete_gradient_trace_have_the_same_physics() {
    let mu = 2.0;
    let b = [0.2, -0.3, 0.4];
    let problem = cube(2, mu);
    let j = vec![[0.0; 3]; problem.mesh().tets.len()];
    let a = uniform_trace(&problem, b, false);
    let shifted = uniform_trace(&problem, b, true);
    let field = problem
        .solve(&j, &trace(&problem, &a), controls(), &mut || false)
        .unwrap();
    let gauged = problem
        .solve(&j, &trace(&problem, &shifted), controls(), &mut || false)
        .unwrap();
    for solved in [&field, &gauged] {
        assert_physical(solved);
        for actual in solved.flux_density_t() {
            for k in 0..3 {
                near(actual[k], b[k]);
            }
        }
        near(solved.stored_energy_j(), 0.5 * dot(b, b) / mu);
        assert!(matches!(
            solved.energy_reference(),
            MagneticEnergyReference::ZeroFlux
        ));
        for (face, &flux) in problem.mesh().faces.iter().zip(solved.face_flux_wb()) {
            let [p, q, r] = face.map(|v| problem.positions_m()[v as usize]);
            let u = std::array::from_fn(|k| q[k] - p[k]);
            let v = std::array::from_fn(|k| r[k] - p[k]);
            near(flux, 0.5 * dot(b, cross(u, v)));
        }
    }
    for (a, b) in field.face_flux_wb().iter().zip(gauged.face_flux_wb()) {
        near(*a, *b);
    }
}
#[test]
fn constant_current_quadratic_potential_has_first_order_magnetic_field_convergence() {
    let mu = 2.0;
    let current = 3.0;
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let problem = cube(n, mu);
        // Exact integral of A_z=-mu*J*x^2/2 along each physical edge.
        let a: Vec<_> = problem
            .mesh()
            .edges
            .iter()
            .map(|&[u, v]| {
                let (p, q) = (
                    problem.positions_m()[u as usize],
                    problem.positions_m()[v as usize],
                );
                let dx = q[0] - p[0];
                -0.5 * mu * current * (q[2] - p[2]) * (p[0] * p[0] + p[0] * dx + dx * dx / 3.0)
            })
            .collect();
        let j = vec![[0.0, 0.0, current]; problem.mesh().tets.len()];
        let field = problem
            .solve(&j, &trace(&problem, &a), controls(), &mut || false)
            .unwrap();
        assert_physical(&field);
        let volume = 1.0 / (6.0 * (n * n * n) as f64);
        let mut squared = 0.0;
        for (tet, b) in problem.mesh().tets.iter().zip(field.flux_density_t()) {
            let e = tet.map(|v| {
                [
                    b[0],
                    b[1] - mu * current * problem.positions_m()[v as usize][0],
                    b[2],
                ]
            });
            let sum: [f64; 3] = std::array::from_fn(|k| e.iter().map(|e| e[k]).sum());
            squared += volume / 20.0 * (dot(sum, sum) + e.iter().map(|&e| dot(e, e)).sum::<f64>());
        }
        errors.push(squared.sqrt());
    }
    for pair in errors.windows(2) {
        assert!(
            (1.7..2.3).contains(&(pair[0] / pair[1])),
            "first-order B refinement: {errors:?}"
        );
    }
    assert!(
        errors[2] < 0.6,
        "resolved current-carrying field error: {errors:?}"
    );
}
fn recoil_source() -> ResolvedRecoilMaterial {
    let values = [
        (
            RECOIL_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            2.0,
        ),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, 0.7),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            -10.0,
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            10.0,
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            10.0,
        ),
    ];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured 3D vector recoil".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic magnet".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    ResolvedRecoilMaterial::try_from_material(
        resolve_material_state_point(
            &card,
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
    .unwrap()
}
#[test]
fn oriented_recoil_natural_field_preserves_shifted_energy_and_operating_limits() {
    let (p, t) = grid(2, 2, 2, 0.5, false);
    let direction = [0.6, 0.0, 0.8];
    let magnet = SpatialRecoilMaterial::new(recoil_source(), direction).unwrap();
    let m = vec![SpatialMagneticMaterial::Recoil(magnet); t.len()];
    let problem = SpatialVectorMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap();
    let j = vec![[0.0; 3]; problem.mesh().tets.len()];
    let field = problem
        .solve(
            &j,
            &VectorMagneticBoundary::default(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_physical(&field);
    for (b, h) in field
        .flux_density_t()
        .iter()
        .zip(field.field_strength_a_m())
    {
        for k in 0..3 {
            near(b[k], 0.7 * direction[k]);
            near(h[k], 0.0);
        }
    }
    near(field.stored_energy_j(), 0.0);
    assert!(field.recoil_operating_points().iter().all(Option::is_some));
    assert!(matches!(
        field.energy_reference(),
        MagneticEnergyReference::FixedRemanence
    ));
    let outside = uniform_trace(&problem, direction.map(|x| x * 50.0), false);
    assert!(matches!(
        problem.solve(&j, &trace(&problem, &outside), controls(), &mut || false),
        Err(SpatialMagneticError::Recoil { .. })
    ));
}
fn current_from_edge_field(
    problem: &SpatialVectorMagnetostatics,
    geometry: &ElementGeometry,
    field: &[f64],
) -> Vec<[f64; 3]> {
    problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .map(|(cell, tet)| {
            let mut current = [0.0; 3];
            for a in 0..4 {
                for b in a + 1..4 {
                    let (a, b) = if tet[a] < tet[b] { (a, b) } else { (b, a) };
                    let edge = problem
                        .mesh()
                        .edges
                        .binary_search(&[tet[a], tet[b]])
                        .unwrap();
                    let curl = cross(geometry.grads[cell][a], geometry.grads[cell][b]);
                    for k in 0..3 {
                        current[k] += 2.0 * field[edge] * curl[k];
                    }
                }
            }
            current
        })
        .collect()
}
#[test]
fn handled_domain_checks_harmonic_source_compatibility_and_complete_coexact_gauge() {
    let (p, t) = grid(3, 3, 1, 1.0, true);
    let m = vec![SpatialMagneticMaterial::Linear(material(2.0)); t.len()];
    let problem = SpatialVectorMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap();
    let geo = try_element_geometry(problem.mesh(), problem.positions_m()).unwrap();
    let exterior: BTreeSet<_> = problem
        .exterior_faces()
        .into_iter()
        .flat_map(|f| problem.face_edges(f).unwrap())
        .collect();
    let edge = (0..problem.mesh().edges.len())
        .find(|e| !exterior.contains(e))
        .unwrap();
    let mut source = vec![0.0; problem.mesh().edges.len()];
    source[edge] = 0.3;
    let current = current_from_edge_field(&problem, &geo, &source);
    let field = problem
        .solve(
            &current,
            &VectorMagneticBoundary::default(),
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_physical(&field);
    assert!(field.stored_energy_j() > 0.0);
    // Wrapped angle increments form a closed, non-exact cochain around the
    // hole. Test the gauge against that harmonic class without random bases.
    let harmonic: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .map(|&[u, v]| {
            let angle = |v: u32| {
                let p = problem.positions_m()[v as usize];
                (p[1] - 1.5).atan2(p[0] - 1.5)
            };
            let mut delta = angle(v) - angle(u);
            if delta > std::f64::consts::PI {
                delta -= std::f64::consts::TAU;
            }
            if delta < -std::f64::consts::PI {
                delta += std::f64::consts::TAU;
            }
            delta
        })
        .collect();
    for row in problem.mesh().d1().rows {
        assert!(
            row.iter()
                .map(|&(e, s)| f64::from(s) * harmonic[e])
                .sum::<f64>()
                .abs()
                < 1e-12
        );
    }
    let metric = hodge_diagonal_barycentric(problem.mesh(), problem.positions_m(), &geo, 1);
    near(
        metric
            .iter()
            .zip(&harmonic)
            .zip(field.potential_wb())
            .map(|((m, h), a)| m * h * a)
            .sum(),
        0.0,
    );
    let mut gauge_div = vec![0.0; problem.mesh().vertex_count];
    for (edge, &[u, v]) in problem.mesh().edges.iter().enumerate() {
        let value = metric[edge] * field.potential_wb()[edge];
        gauge_div[u as usize] -= value;
        gauge_div[v as usize] += value;
    }
    assert!(gauge_div.iter().all(|v| v.abs() < 1e-9));
    // Stream function is constant on each boundary, but different on inner
    // and outer boundaries: J is divergence-free and circulates through the
    // handle. Natural n×H=0 cannot support that linked current.
    let psi: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| {
            if (p[0] == 1.0 || p[0] == 2.0) && (p[1] == 1.0 || p[1] == 2.0) {
                1.0
            } else {
                0.0
            }
        })
        .collect();
    let incompatible: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .map(|(cell, tet)| {
            let grad: [f64; 3] = std::array::from_fn(|k| {
                (0..4)
                    .map(|a| psi[tet[a] as usize] * geo.grads[cell][a][k])
                    .sum()
            });
            [grad[1], -grad[0], 0.0]
        })
        .collect();
    let mut charge = vec![0.0; problem.mesh().vertex_count];
    for (cell, tet) in problem.mesh().tets.iter().enumerate() {
        for a in 0..4 {
            charge[tet[a] as usize] +=
                geo.vol_signed[cell].abs() * dot(incompatible[cell], geo.grads[cell][a]);
        }
    }
    assert!(charge.iter().all(|q| q.abs() < 1e-12));
    assert!(matches!(
        problem.solve(
            &incompatible,
            &VectorMagneticBoundary::default(),
            controls(),
            &mut || false
        ),
        Err(SpatialMagneticError::Invalid(
            "impressed current contains a gradient or harmonic component incompatible with the magnetic boundary"
        ))
    ));
}
#[test]
fn unsupported_current_traces_work_budgets_and_cancellation_publish_no_field() {
    let problem = cube(2, 2.0);
    let zeros = vec![[0.0; 3]; problem.mesh().tets.len()];
    let unclosed = vec![[0.0, 0.0, 1.0]; problem.mesh().tets.len()];
    assert!(matches!(
        problem.solve(
            &unclosed,
            &VectorMagneticBoundary::default(),
            controls(),
            &mut || false
        ),
        Err(SpatialMagneticError::Invalid(
            "impressed current violates interface or natural-exterior charge continuity"
        ))
    ));
    let trace = trace(&problem, &uniform_trace(&problem, [0.2, 0.3, 0.4], false));
    let mut duplicate = trace.clone();
    duplicate.potential.push(duplicate.potential[0]);
    assert!(
        problem
            .solve(&zeros, &duplicate, controls(), &mut || false)
            .is_err()
    );
    let mut limited = controls();
    limited.maximum_gauge_product_terms = 0;
    assert!(
        problem
            .solve(&zeros, &trace, limited, &mut || false)
            .is_err()
    );
    let mut limited = controls();
    limited.field.max_iterations = 0;
    assert!(matches!(
        problem.solve(&zeros, &trace, limited, &mut || false),
        Err(SpatialMagneticError::Nonconvergence(_))
    ));
    let mut calls = 0;
    problem
        .solve(&zeros, &trace, controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut checks = 0;
    assert!(matches!(
        problem.solve(&zeros, &trace, controls(), &mut || {
            checks += 1;
            checks > calls * 2 / 3
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    let field = problem
        .solve(&zeros, &trace, controls(), &mut || false)
        .unwrap();
    assert_physical(&field);
}

#[test]
fn vector_maxwell_wrench_uses_the_closed_current_free_exterior_collar() {
    let problem = cube(3, 2.0);
    let selected: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(cell, tet)| {
            let center: [f64; 3] = std::array::from_fn(|k| {
                tet.iter()
                    .map(|&v| problem.positions_m()[v as usize][k] / 4.0)
                    .sum()
            });
            center
                .iter()
                .all(|&x| x > 1.0 / 3.0 && x < 2.0 / 3.0)
                .then_some(cell)
        })
        .collect();
    assert_eq!(selected.len(), 6);
    let boundary = trace(&problem, &uniform_trace(&problem, [0.2, -0.3, 0.4], false));
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let field = problem
        .solve(&zero, &boundary, controls(), &mut || false)
        .unwrap();
    let limits = fs_em::spatial::force::SpatialMaxwellForceLimits {
        max_cells: 200,
        max_faces: 1000,
    };
    let wrench = field
        .maxwell_wrench(
            &selected,
            &fs_ga::Se3::identity(),
            [0.1, -0.7, 0.3],
            limits,
            &mut || false,
        )
        .unwrap();
    for component in [
        wrench.force.x,
        wrench.force.y,
        wrench.force.z,
        wrench.torque.x,
        wrench.torque.y,
        wrench.torque.z,
    ] {
        assert!(component.abs() < 1e-9);
    }
    let current = vec![[0.0, 0.0, 1.0]; problem.mesh().tets.len()];
    let current_field = problem
        .solve(&current, &boundary, controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        current_field.maxwell_wrench(
            &selected,
            &fs_ga::Se3::identity(),
            [0.0; 3],
            limits,
            &mut || false,
        ),
        Err(SpatialMagneticError::Invalid(
            "Maxwell surface requires a current-free exterior collar"
        ))
    ));
}
