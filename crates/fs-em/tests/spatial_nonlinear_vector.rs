//! G1/G3/G4 current-driven nonlinear 3D fields from original source curves.
//! Manufactured sources verify physics/discretization, not measured materials.
mod common;
use common::{material, near};
use fs_em::spatial::{
    SpatialMagneticError, SpatialMagneticMaterial, SpatialRecoilMaterial,
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        nonlinear::{
            NonlinearVectorControls, NonlinearVectorError, SpatialNonlinearVectorField,
            SpatialNonlinearVectorMagnetostatics,
        },
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
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
use fs_qty::QuantitySpec;
use fs_solver::CheckedCgConfig;
fn controls() -> NonlinearVectorControls {
    NonlinearVectorControls {
        field: VectorMagneticControls {
            field: CheckedCgConfig {
                tolerance: 2e-11,
                max_iterations: 3000,
                max_corrections: 2,
            },
            gauge: CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 3000,
                max_corrections: 2,
            },
            maximum_gauge_product_terms: 2_000_000,
            current_tolerance_a: 2e-8,
            flux_tolerance_wb: 3e-9,
            energy_tolerance_j: 8e-8,
        },
        relative_residual_tolerance: 2e-10,
        maximum_newton_iterations: 35,
        maximum_backtracks: 35,
    }
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn homogeneous(n: usize, nonlinear: bool) -> SpatialNonlinearVectorMagnetostatics {
    let (p, t) = cube(n);
    let m = if nonlinear {
        SpatialNonlinearMaterial::Reversible(law())
    } else {
        SpatialNonlinearMaterial::Linear(material(2.0))
    };
    SpatialNonlinearVectorMagnetostatics::new(p, t.clone(), vec![m; t.len()], limits(), &mut || {
        false
    })
    .unwrap()
}
fn trace(p: &SpatialNonlinearVectorMagnetostatics, a: &[f64]) -> VectorMagneticBoundary {
    VectorMagneticBoundary {
        potential: p
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: p.face_edges(face).unwrap().map(|e| a[e]),
            })
            .collect(),
    }
}
fn uniform(p: &SpatialNonlinearVectorMagnetostatics, b: [f64; 3], shift: bool) -> Vec<f64> {
    p.mesh()
        .edges
        .iter()
        .map(|&[i, j]| {
            let (x, y) = (p.positions_m()[i as usize], p.positions_m()[j as usize]);
            let mid = std::array::from_fn(|k| 0.5 * (x[k] + y[k]));
            let delta = std::array::from_fn(|k| y[k] - x[k]);
            let psi = |v: u32| {
                if shift {
                    ((v * 7) % 11) as f64 * 0.17
                } else {
                    0.0
                }
            };
            0.5 * dot(cross(b, mid), delta) + psi(j) - psi(i)
        })
        .collect()
}
fn physical(f: &SpatialNonlinearVectorField<'_>) {
    assert!(f.report.maximum_residual_a <= controls().field.current_tolerance_a);
    assert!(f.report.source_compatibility_defect_a <= controls().field.current_tolerance_a);
    assert!(f.report.maximum_gauge_flux_defect_wb <= controls().field.flux_tolerance_wb);
    assert!(f.maximum_divergence_flux_wb <= controls().field.flux_tolerance_wb);
    assert!(f.maximum_recovered_flux_defect_wb <= controls().field.flux_tolerance_wb);
    assert!(f.work_residual_j.abs() <= controls().field.energy_tolerance_j);
    assert!(f.legendre_defect_j <= controls().field.energy_tolerance_j);
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 4000,
        max_dofs: 20000,
        max_triplets: 150000,
    }
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn law() -> ReversibleMagneticMaterial {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS),
            value: PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots: vec![(0.0, 0.0), (1.0, 0.5), (2.0, 2.5), (4.0, 10.5)],
                dims: FIELD_STRENGTH_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "B",
                QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
                0.0,
                4.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured 3D saturation law".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic soft magnet".into(),
            phase: "solid".into(),
            process: "reversible fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap();
    ReversibleMagneticMaterial::resolve(
        &card,
        &query,
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
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
fn linear_reduction_matches_original_current_owner_and_gradient_trace_shift() {
    let p = homogeneous(3, false);
    let original = SpatialVectorMagnetostatics::new(
        p.positions_m().to_vec(),
        p.mesh().tets.clone(),
        vec![SpatialMagneticMaterial::Linear(material(2.0)); p.mesh().tets.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let a: Vec<_> = p
        .mesh()
        .edges
        .iter()
        .map(|&[i, j]| {
            let (x, y) = (p.positions_m()[i as usize], p.positions_m()[j as usize]);
            let dx = y[0] - x[0];
            -(y[2] - x[2]) * (x[0] * x[0] + x[0] * dx + dx * dx / 3.0)
        })
        .collect();
    let j = vec![[0.0, 0.0, 1.0]; p.mesh().tets.len()];
    let field = p
        .solve(&j, &trace(&p, &a), &a, controls(), &mut || false)
        .unwrap();
    let old = original
        .solve(&j, &trace(&p, &a), controls().field, &mut || false)
        .unwrap();
    physical(&field);
    for (a, b) in field
        .flux_density_t()
        .iter()
        .flatten()
        .zip(old.flux_density_t().iter().flatten())
    {
        assert!((a - b).abs() < 2e-8);
    }
    near(field.stored_energy_j(), old.stored_energy_j());
    let shifted: Vec<_> = p
        .mesh()
        .edges
        .iter()
        .zip(&a)
        .map(|(&[i, j], &a)| a + ((j * 7) % 11) as f64 * 0.17 - ((i * 7) % 11) as f64 * 0.17)
        .collect();
    let other = p
        .solve(&j, &trace(&p, &shifted), &shifted, controls(), &mut || {
            false
        })
        .unwrap();
    physical(&other);
    for (a, b) in field.face_flux_wb().iter().zip(other.face_flux_wb()) {
        assert!((a - b).abs() < 2e-9);
    }
}

#[test]
fn saturated_current_quadratic_vector_potential_has_first_order_field_convergence() {
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let xs: Vec<_> = (0..=n).map(|i| 1.0 + i as f64 / n as f64).collect();
        let ys: Vec<_> = (0..=n).map(|i| i as f64 / n as f64).collect();
        let (positions, tets) = grid(&xs, &ys, &ys);
        let p = SpatialNonlinearVectorMagnetostatics::new(
            positions,
            tets.clone(),
            vec![SpatialNonlinearMaterial::Reversible(law()); tets.len()],
            limits(),
            &mut || false,
        )
        .unwrap();
        // H_y=x, J_z=dH_y/dx=1; on this saturated source segment B_y=x/2+3/4.
        // A_z=-x^2/4-3x/4; these are exact physical edge integrals.
        let a: Vec<_> = p
            .mesh()
            .edges
            .iter()
            .map(|&[i, j]| {
                let (x, y) = (p.positions_m()[i as usize], p.positions_m()[j as usize]);
                let dx = y[0] - x[0];
                -(y[2] - x[2])
                    * (0.25 * (x[0] * x[0] + x[0] * dx + dx * dx / 3.0) + 0.75 * (x[0] + 0.5 * dx))
            })
            .collect();
        let j = vec![[0.0, 0.0, 1.0]; p.mesh().tets.len()];
        let field = p
            .solve(&j, &trace(&p, &a), &a, controls(), &mut || false)
            .unwrap();
        physical(&field);
        assert!(field.report.newton_iterations > 0);
        let volume = 1.0 / (6.0 * (n * n * n) as f64);
        let mut b_error = 0.0;
        let mut h_error = 0.0;
        let u = (5.0_f64 + 3.0 * 5.0_f64.sqrt()) / 20.0;
        let v = (5.0_f64 - 5.0_f64.sqrt()) / 20.0;
        for (cell, tet) in p.mesh().tets.iter().enumerate() {
            for corner in 0..4 {
                let x = (0..4)
                    .map(|k| if k == corner { u } else { v } * p.positions_m()[tet[k] as usize][0])
                    .sum::<f64>();
                let b = [0.0, 0.5 * x + 0.75, 0.0];
                let h = [0.0, x, 0.0];
                b_error += 0.25
                    * volume
                    * (0..3)
                        .map(|k| (field.flux_density_t()[cell][k] - b[k]).powi(2))
                        .sum::<f64>();
                h_error += 0.25
                    * volume
                    * (0..3)
                        .map(|k| (field.field_strength_a_m()[cell][k] - h[k]).powi(2))
                        .sum::<f64>();
            }
        }
        errors.push((b_error.sqrt(), h_error.sqrt()));
        assert!(
            field
                .flux_density_t()
                .iter()
                .all(|b| dot(*b, *b).sqrt() > 1.0)
        );
        assert!(
            (field.stored_energy_j() - 0.5 * (field.stored_energy_j() + field.coenergy_j())).abs()
                > 0.05
        );
    }
    for pair in errors.windows(2) {
        for (a, b) in [(pair[0].0, pair[1].0), (pair[0].1, pair[1].1)] {
            let order = (a / b).log2();
            assert!((0.8..1.2).contains(&order), "{errors:?}: order={order}");
        }
    }
}

#[test]
fn uniform_rotated_saturation_retains_source_energy_coenergy_and_derivative() {
    let p = homogeneous(2, true);
    let b = [0.9, 1.2, 0.0];
    let j = vec![[0.0; 3]; p.mesh().tets.len()];
    let a = uniform(&p, b, true);
    let f = p
        .solve(&j, &trace(&p, &a), &a, controls(), &mut || false)
        .unwrap();
    physical(&f);
    let response = law().evaluate(1.5).unwrap();
    near(f.stored_energy_j(), response.energy_density_j_m3);
    near(f.coenergy_j(), response.coenergy_density_j_m3);
    for actual in f.flux_density_t() {
        for k in 0..3 {
            near(actual[k], b[k]);
        }
    }
    let eps = 1e-5;
    for k in 0..3 {
        let mut bp = b;
        bp[k] += eps;
        let mut bm = b;
        bm[k] -= eps;
        let ap = uniform(&p, bp, false);
        let am = uniform(&p, bm, false);
        let fp = p
            .solve(&j, &trace(&p, &ap), &ap, controls(), &mut || false)
            .unwrap();
        let fm = p
            .solve(&j, &trace(&p, &am), &am, controls(), &mut || false)
            .unwrap();
        let h = 1.5 * b[k] / 1.5; // H magnitude=1.5 at B=1.5 on the supplied segment.
        assert!(((fp.stored_energy_j() - fm.stored_energy_j()) / (2.0 * eps) - h).abs() < 2e-8);
    }
}

#[test]
fn recoil_and_saturated_layers_share_tangential_h_and_shifted_energy() {
    let (p, t) = grid(&[0.0, 1.0, 2.0, 3.0], &[0.0, 0.5, 1.0], &[0.0, 0.5, 1.0]);
    let m: Vec<_> = t
        .iter()
        .map(|tet| {
            let x = tet.iter().map(|&v| p[v as usize][0]).sum::<f64>() / 4.0;
            if x < 1.0 {
                SpatialNonlinearMaterial::Recoil(
                    SpatialRecoilMaterial::new(recoil(1.0, 5.0, 10.0), [0.0, 0.0, 1.0]).unwrap(),
                )
            } else if x < 2.0 {
                SpatialNonlinearMaterial::Reversible(law())
            } else {
                SpatialNonlinearMaterial::Linear(material(1.0))
            }
        })
        .collect();
    let p = SpatialNonlinearVectorMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap();
    let ay = |x: f64| {
        if x <= 1.0 {
            5.8 * x
        } else if x <= 2.0 {
            5.8 + 1.15 * (x - 1.0)
        } else {
            6.95 + 0.8 * (x - 2.0)
        }
    };
    let a: Vec<_> = p
        .mesh()
        .edges
        .iter()
        .map(|&[i, j]| {
            let (x, y) = (p.positions_m()[i as usize], p.positions_m()[j as usize]);
            (y[1] - x[1]) * 0.5 * (ay(x[0]) + ay(y[0]))
        })
        .collect();
    let j = vec![[0.0; 3]; p.mesh().tets.len()];
    let f = p
        .solve(&j, &trace(&p, &a), &a, controls(), &mut || false)
        .unwrap();
    physical(&f);
    near(f.stored_energy_j(), 0.9875);
    near(f.coenergy_j(), 5.2125);
    for h in f.field_strength_a_m() {
        near(h[0], 0.0);
        near(h[1], 0.0);
        near(h[2], 0.8);
    }
    assert_eq!(
        f.recoil_operating_points()
            .iter()
            .filter(|v| v.is_some())
            .count(),
        24
    );
}

#[test]
fn real_source_energy_iteration_and_cancellation_failures_return_no_field() {
    let p = homogeneous(2, true);
    let zero = vec![0.0; p.mesh().edges.len()];
    let j = vec![[0.0, 0.0, 0.05]; p.mesh().tets.len()];
    let boundary = trace(&p, &zero);
    let mut c = controls();
    c.maximum_newton_iterations = 0;
    assert!(matches!(
        p.solve(&j, &boundary, &zero, c, &mut || false),
        Err(NonlinearVectorError::Nonconvergence(_))
    ));
    c = controls();
    c.field.field.max_iterations = 0;
    c.field.field.max_corrections = 0;
    assert!(matches!(
        p.solve(&j, &boundary, &zero, c, &mut || false),
        Err(NonlinearVectorError::Spatial(
            SpatialMagneticError::Nonconvergence(_)
        ))
    ));
    assert!(
        p.solve(
            &j,
            &VectorMagneticBoundary::default(),
            &zero,
            controls(),
            &mut || false
        )
        .is_err()
    );
    let out = uniform(&p, [0.0, 0.0, 8.0], false);
    assert!(matches!(
        p.solve(
            &vec![[0.0; 3]; j.len()],
            &trace(&p, &out),
            &out,
            controls(),
            &mut || false
        ),
        Err(NonlinearVectorError::FluxRange { .. })
    ));
    // A loose current allowance must not convert an unsolved interior field
    // into fictitious boundary supply: all prescribed edges are zero.
    let noisy: Vec<_> = (0..zero.len())
        .map(|i| 0.001 * ((i * 7) % 11) as f64)
        .collect();
    c = controls();
    c.field.current_tolerance_a = 10.0;
    c.field.energy_tolerance_j = 1e-12;
    c.maximum_newton_iterations = 0;
    assert!(matches!(
        p.solve(&vec![[0.0; 3]; j.len()], &boundary, &noisy, c, &mut || {
            false
        }),
        Err(NonlinearVectorError::Spatial(
            SpatialMagneticError::Invalid(_)
        ))
    ));
    let accepted = p
        .solve(&j, &boundary, &zero, controls(), &mut || false)
        .unwrap();
    physical(&accepted);
    let mut polls = 0;
    assert!(
        p.solve(&j, &boundary, &zero, controls(), &mut || {
            polls += 1;
            polls > 120
        })
        .is_err()
    );
    let again = p
        .solve(&j, &boundary, &zero, controls(), &mut || false)
        .unwrap();
    assert_eq!(accepted.potential_wb(), again.potential_wb());
}

#[test]
fn nonlinear_vector_field_uses_shared_static_force_admission() {
    use fs_em::spatial::force::SpatialMaxwellForceLimits;
    let p = homogeneous(4, false);
    let a = uniform(&p, [0.2, 0.3, 0.4], false);
    let j = vec![[0.0; 3]; p.mesh().tets.len()];
    let f = p
        .solve(&j, &trace(&p, &a), &a, controls(), &mut || false)
        .unwrap();
    let selected: Vec<_> = p
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(i, t)| {
            let mid: [f64; 3] = std::array::from_fn(|k| {
                t.iter()
                    .map(|&v| p.positions_m()[v as usize][k])
                    .sum::<f64>()
                    / 4.0
            });
            mid.iter().all(|&v| v > 0.25 && v < 0.75).then_some(i)
        })
        .collect();
    let limits = SpatialMaxwellForceLimits {
        max_cells: 4000,
        max_faces: 20000,
    };
    let wrench = f
        .maxwell_wrench(
            &selected,
            &fs_ga::Se3::identity(),
            [0.0; 3],
            limits,
            &mut || false,
        )
        .unwrap();
    assert!(wrench.to_array().iter().all(|v| v.abs() < 1e-10));
    let nonlinear = homogeneous(4, true);
    let a = uniform(&nonlinear, [0.2, 0.3, 0.4], false);
    let f = nonlinear
        .solve(&j, &trace(&nonlinear, &a), &a, controls(), &mut || false)
        .unwrap();
    assert!(
        f.maxwell_wrench(
            &selected,
            &fs_ga::Se3::identity(),
            [0.0; 3],
            limits,
            &mut || false
        )
        .is_err()
    );
}
