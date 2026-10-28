//! Analytical current-free 3D fields, magnet interfaces, and mesh refinement.
//! Manufactured material sources test the discretization, not material quality.
mod common;
use common::{material, near};
use fs_em::{
    planar::permanent::MagneticEnergyReference,
    spatial::{
        MagneticFluxFace, MagneticPotentialFace, SpatialMagneticBoundary, SpatialMagneticControls,
        SpatialMagneticError, SpatialMagneticMaterial, SpatialMagnetostatics,
        SpatialRecoilMaterial,
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::{
    weighted::WeightedAssemblyLimits,
    whitney::{ElementGeometryError, try_element_geometry},
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
use fs_rep_mesh::TetComplex;
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 4000,
        max_dofs: 20000,
        max_triplets: 150000,
    }
}
fn controls() -> SpatialMagneticControls {
    SpatialMagneticControls {
        cg: common::config(),
        flux_tolerance_wb: 1e-10,
        energy_tolerance_j: 2e-6,
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn grid(xs: &[f64], ys: &[f64], zs: &[f64]) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let mut p = Vec::new();
    for &z in zs {
        for &y in ys {
            for &x in xs {
                p.push([x, y, z]);
            }
        }
    }
    let id = |i, j, k| ((k * ys.len() + j) * xs.len() + i) as u32;
    let mut t = Vec::new();
    for k in 0..zs.len() - 1 {
        for j in 0..ys.len() - 1 {
            for i in 0..xs.len() - 1 {
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
fn cube(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let q: Vec<_> = (0..=n).map(|i| i as f64 / n as f64).collect();
    grid(&q, &q, &q)
}
fn linear(n: usize, mu: f64) -> SpatialMagnetostatics {
    let (p, t) = cube(n);
    let m = vec![SpatialMagneticMaterial::Linear(material(mu)); t.len()];
    SpatialMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap()
}
fn potentials(
    problem: &SpatialMagnetostatics,
    phi: impl Fn([f64; 3]) -> f64,
) -> SpatialMagneticBoundary {
    SpatialMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| MagneticPotentialFace {
                face,
                potential_a: problem.mesh().faces[face]
                    .map(|v| phi(problem.positions_m()[v as usize])),
            })
            .collect(),
        flux: vec![],
    }
}
fn box_normal(problem: &SpatialMagnetostatics, face: usize) -> [f64; 3] {
    let p = problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize]);
    for axis in 0..3 {
        for bound in [0.0, 1.0] {
            if p.iter().all(|p| p[axis] == bound) {
                let mut n = [0.0; 3];
                n[axis] = if bound == 0.0 { -1.0 } else { 1.0 };
                return n;
            }
        }
    }
    panic!("expected unit-box exterior face")
}
fn recoil(mu: f64, br: f64, limit: f64) -> ResolvedRecoilMaterial {
    let values = [
        (RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, br),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            -limit,
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            limit,
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            limit,
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
                    source: "manufactured 3D recoil layer".into(),
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
fn uniform_vector_field_matches_dirichlet_neumann_and_large_gauge_shift() {
    let mu = 1.3e-6;
    let h = [200.0, -150.0, 70.0];
    let problem = linear(2, mu);
    let original = problem
        .solve(
            &potentials(&problem, |p| -dot(h, p)),
            controls(),
            &mut || false,
        )
        .unwrap();
    let offset = 2.0_f64.powi(36);
    let shifted = problem
        .solve(
            &potentials(&problem, |p| offset - dot(h, p)),
            controls(),
            &mut || false,
        )
        .unwrap();
    let flux = SpatialMagneticBoundary {
        potential: vec![],
        flux: problem
            .exterior_faces()
            .into_iter()
            .map(|face| {
                let n = box_normal(&problem, face);
                for (actual, expected) in problem.exterior_normal(face).unwrap().into_iter().zip(n)
                {
                    near(actual, expected);
                }
                MagneticFluxFace {
                    face,
                    outward_flux_density_t: mu * dot(h, n),
                }
            })
            .collect(),
    };
    let neumann = problem.solve(&flux, controls(), &mut || false).unwrap();
    for field in [&original, &shifted, &neumann] {
        for (actual_h, actual_b) in field
            .field_strength_a_m()
            .iter()
            .zip(field.flux_density_t())
        {
            for axis in 0..3 {
                near(actual_h[axis], h[axis]);
                near(actual_b[axis], mu * h[axis]);
            }
        }
        near(field.stored_energy_j(), 0.5 * mu * dot(h, h));
        assert!(field.maximum_weak_residual_wb() < controls().flux_tolerance_wb);
        assert!(
            field
                .interface_flux_jumps_wb()
                .iter()
                .all(|(_, jump)| jump.abs() < 1e-11)
        );
        assert!(matches!(
            field.energy_reference(),
            MagneticEnergyReference::ZeroFlux
        ));
    }
    assert_eq!(neumann.gauge_vertices(), &[0]);
    for (i, p) in problem.positions_m().iter().enumerate() {
        near(neumann.potential_a()[i], -dot(h, *p));
    }
    for (a, b) in original
        .relative_potential_a()
        .iter()
        .zip(shifted.relative_potential_a())
    {
        near(*a, *b);
    }
}
#[test]
fn recoil_air_layers_match_analytic_flux_demagnetizing_field_and_energy() {
    let mu_m = 2e-6;
    let mu_air = 1e-6;
    let br = 0.9;
    let lm = 0.4;
    let la = 0.6;
    let expected_b = (br * lm / mu_m) / (lm / mu_m + la / mu_air);
    for direction in [1.0, -1.0] {
        let (p, t) = grid(&[0.0, lm, 1.0], &[0.0, 0.5, 1.0], &[0.0, 0.5, 1.0]);
        let magnet =
            SpatialRecoilMaterial::new(recoil(mu_m, br, 1e6), [direction, 0.0, 0.0]).unwrap();
        let left: Vec<_> = t
            .iter()
            .map(|tet| tet.iter().map(|&v| p[v as usize][0] / 4.0).sum::<f64>() < lm)
            .collect();
        let materials = left
            .iter()
            .map(|&left| {
                if left {
                    SpatialMagneticMaterial::Recoil(magnet.clone())
                } else {
                    SpatialMagneticMaterial::Linear(material(mu_air))
                }
            })
            .collect();
        let problem = SpatialMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap();
        let boundary = SpatialMagneticBoundary {
            potential: problem
                .exterior_faces()
                .into_iter()
                .filter(|&face| {
                    let x =
                        problem.mesh().faces[face].map(|v| problem.positions_m()[v as usize][0]);
                    x.iter().all(|x| *x == 0.0) || x.iter().all(|x| *x == 1.0)
                })
                .map(|face| MagneticPotentialFace {
                    face,
                    potential_a: [0.0; 3],
                })
                .collect(),
            flux: vec![],
        };
        let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
        for (i, &left) in left.iter().enumerate() {
            let expected_h = direction
                * if left {
                    (expected_b - br) / mu_m
                } else {
                    expected_b / mu_air
                };
            near(field.field_strength_a_m()[i][0], expected_h);
            near(field.flux_density_t()[i][0], direction * expected_b);
            for axis in 1..3 {
                assert!(field.field_strength_a_m()[i][axis].abs() < 2e-6);
                near(field.flux_density_t()[i][axis], 0.0);
            }
            assert_eq!(field.recoil_operating_points()[i].is_some(), left);
        }
        near(
            field.stored_energy_j(),
            0.5 * ((expected_b - br).powi(2) * lm / mu_m + expected_b.powi(2) * la / mu_air),
        );
        assert!(matches!(
            field.energy_reference(),
            MagneticEnergyReference::FixedRemanence
        ));
        assert!(
            field
                .interface_flux_jumps_wb()
                .iter()
                .all(|(_, jump)| jump.abs() < 1e-10)
        );
        let mut end_flux = [0.0; 2];
        for (p, &flux) in problem
            .positions_m()
            .iter()
            .zip(field.nodal_boundary_flux_wb())
        {
            if p[0] == 0.0 {
                end_flux[0] += flux;
            } else if p[0] == 1.0 {
                end_flux[1] += flux;
            }
        }
        near(end_flux[0], -direction * expected_b);
        near(end_flux[1], direction * expected_b);
    }
}
#[test]
fn harmonic_potential_has_first_order_field_convergence_in_three_dimensions() {
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let problem = linear(n, 1.7e-6);
        let field = problem
            .solve(
                &potentials(&problem, |p| p[0] * p[0] - p[1] * p[1]),
                controls(),
                &mut || false,
            )
            .unwrap();
        let volume = 1.0 / (6.0 * (n * n * n) as f64);
        let mut squared = 0.0;
        for (tet, h) in problem.mesh().tets.iter().zip(field.field_strength_a_m()) {
            let e = tet.map(|v| {
                let p = problem.positions_m()[v as usize];
                [h[0] + 2.0 * p[0], h[1] - 2.0 * p[1], h[2]]
            });
            let sum: [f64; 3] = std::array::from_fn(|j| e.iter().map(|e| e[j]).sum());
            // Exact tetrahedral quadratic integral of the affine field error.
            squared += volume / 20.0 * (dot(sum, sum) + e.iter().map(|&e| dot(e, e)).sum::<f64>());
        }
        errors.push(squared.sqrt());
    }
    for pair in errors.windows(2) {
        let rate = pair[0] / pair[1];
        assert!(
            (1.9..2.1).contains(&rate),
            "first-order H convergence: {errors:?}"
        );
    }
    assert!(errors[2] < 0.12, "finest physical field error: {errors:?}");
}
#[test]
fn disconnected_domains_get_independent_gauges_and_require_local_flux_balance() {
    let (mut p, mut t) = cube(1);
    let (q, u) = cube(1);
    let offset = p.len() as u32;
    p.extend(q.into_iter().map(|p| [p[0] + 2.0, p[1], p[2]]));
    t.extend(u.into_iter().map(|tet| tet.map(|v| v + offset)));
    let m = vec![SpatialMagneticMaterial::Linear(material(1e-6)); t.len()];
    let problem = SpatialMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap();
    let zero = problem
        .solve(&SpatialMagneticBoundary::default(), controls(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(zero.gauge_vertices(), &[0, 8]);
    near(zero.stored_energy_j(), 0.0);
    let a = problem
        .exterior_faces()
        .into_iter()
        .find(|&f| problem.mesh().faces[f].iter().all(|&v| v < 8))
        .unwrap();
    let b = problem
        .exterior_faces()
        .into_iter()
        .find(|&f| problem.mesh().faces[f].iter().all(|&v| v >= 8))
        .unwrap();
    let incompatible = SpatialMagneticBoundary {
        potential: vec![],
        flux: vec![
            MagneticFluxFace {
                face: a,
                outward_flux_density_t: 0.1,
            },
            MagneticFluxFace {
                face: b,
                outward_flux_density_t: -0.1,
            },
        ],
    };
    assert!(matches!(
        problem.solve(&incompatible, controls(), &mut || false),
        Err(SpatialMagneticError::Invalid(
            "net outward magnetic flux must vanish in a pure-flux component"
        ))
    ));
}
#[test]
fn checked_geometry_and_volume_admission_refuse_invalid_external_cells() {
    let good = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let tetra = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let mut flat = good.clone();
    flat[3] = [0.5, 0.5, 0.0];
    assert!(matches!(
        try_element_geometry(&tetra, &flat),
        Err(ElementGeometryError::Cell(0))
    ));
    let mut nan = good.clone();
    nan[2][1] = f64::NAN;
    assert!(matches!(
        try_element_geometry(&tetra, &nan),
        Err(ElementGeometryError::Positions)
    ));
    assert!(matches!(
        try_element_geometry(&tetra, &good[..3]),
        Err(ElementGeometryError::Positions)
    ));
    let constructor = |p: Vec<[f64; 3]>, t: Vec<[u32; 4]>| {
        let m = vec![SpatialMagneticMaterial::Linear(material(1e-6)); t.len()];
        SpatialMagnetostatics::new(p, t, m, limits(), &mut || false)
    };
    assert!(constructor(flat, vec![[0, 1, 2, 3]]).is_err());
    assert!(constructor(good.clone(), vec![[0, 1, 2, 3], [3, 2, 1, 0]]).is_err());
    assert!(constructor(good.clone(), vec![[0, 1, 1, 3]]).is_err());
    assert!(constructor(good.clone(), vec![[0, 1, 2, 9]]).is_err());
    let mut same_side = good;
    same_side.push([0.2, 0.2, 2.0]);
    assert!(matches!(
        constructor(same_side, vec![[0, 1, 2, 3], [0, 1, 2, 4]]),
        Err(SpatialMagneticError::Invalid(
            "neighbor tetrahedra must lie on opposite face sides"
        ))
    ));
}
#[test]
fn boundary_traces_must_be_exterior_unique_and_consistent() {
    let problem = linear(1, 1e-6);
    let exterior = problem.exterior_faces();
    let interior = (0..problem.mesh().faces.len())
        .find(|&f| !exterior.contains(&f))
        .unwrap();
    let trace = MagneticPotentialFace {
        face: exterior[0],
        potential_a: [0.0; 3],
    };
    let cases = [
        SpatialMagneticBoundary {
            potential: vec![MagneticPotentialFace {
                face: interior,
                potential_a: [0.0; 3],
            }],
            flux: vec![],
        },
        SpatialMagneticBoundary {
            potential: vec![trace, trace],
            flux: vec![],
        },
        SpatialMagneticBoundary {
            potential: vec![trace],
            flux: vec![MagneticFluxFace {
                face: trace.face,
                outward_flux_density_t: 0.0,
            }],
        },
        SpatialMagneticBoundary {
            potential: vec![],
            flux: vec![MagneticFluxFace {
                face: exterior[0],
                outward_flux_density_t: 0.1,
            }],
        },
    ];
    for boundary in cases {
        assert!(problem.solve(&boundary, controls(), &mut || false).is_err());
    }
    let shared = problem.mesh().faces[trace.face][0];
    let neighbor = exterior
        .into_iter()
        .find(|&f| f != trace.face && problem.mesh().faces[f].contains(&shared))
        .unwrap();
    let conflicting = SpatialMagneticBoundary {
        potential: vec![
            trace,
            MagneticPotentialFace {
                face: neighbor,
                potential_a: problem.mesh().faces[neighbor]
                    .map(|v| if v == shared { 1.0 } else { 0.0 }),
            },
        ],
        flux: vec![],
    };
    assert!(matches!(
        problem.solve(&conflicting, controls(), &mut || false),
        Err(SpatialMagneticError::Invalid(
            "inconsistent potential at shared trace vertex"
        ))
    ));
}
#[test]
fn recoil_limits_work_budgets_and_late_cancellation_refuse_partial_fields() {
    let source = recoil(1e-6, 0.8, 100.0);
    assert!(SpatialRecoilMaterial::new(source.clone(), [1.0, 1.0, 0.0]).is_err());
    let magnet = SpatialRecoilMaterial::new(source, [0.0, 0.0, 1.0]).unwrap();
    let (p, t) = cube(1);
    let m = vec![SpatialMagneticMaterial::Recoil(magnet); t.len()];
    let problem = SpatialMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap();
    assert!(matches!(
        problem.solve(
            &potentials(&problem, |p| -200.0 * p[2]),
            controls(),
            &mut || false
        ),
        Err(SpatialMagneticError::Recoil { .. })
    ));
    assert!(matches!(
        problem.solve(
            &potentials(&problem, |p| -200.0 * p[0]),
            controls(),
            &mut || false
        ),
        Err(SpatialMagneticError::Recoil { .. })
    ));
    let problem = linear(3, 1e-6);
    let boundary = potentials(&problem, |p| -p[0] - 2.0 * p[1] - 3.0 * p[2]);
    let mut exhausted = controls();
    exhausted.cg.max_iterations = 0;
    assert!(matches!(
        problem.solve(&boundary, exhausted, &mut || false),
        Err(SpatialMagneticError::Nonconvergence(_))
    ));
    let mut calls = 0;
    problem
        .solve(&boundary, controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut checks = 0;
    assert!(matches!(
        problem.solve(&boundary, controls(), &mut || {
            checks += 1;
            checks > calls * 2 / 3
        }),
        Err(SpatialMagneticError::Cancelled)
    ));
    let accepted = problem.solve(&boundary, controls(), &mut || false).unwrap();
    near(accepted.stored_energy_j(), 7e-6);
    let (p, t) = cube(1);
    let m = vec![SpatialMagneticMaterial::Linear(material(1e-6)); t.len()];
    assert!(matches!(
        SpatialMagnetostatics::new(p, t, m, limits(), &mut || true),
        Err(SpatialMagneticError::Cancelled)
    ));
}
