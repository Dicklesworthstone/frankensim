//! G1/G3/G4: exact sourced saturation, dual energy, recoil bias and 3D refinement.
//! Synthetic material cards verify the equations, not measured steel or magnets.
mod common;
use common::{material, near};
use fs_em::spatial::{
    MagneticFluxFace, MagneticPotentialFace, SpatialMagneticBoundary, SpatialMagneticControls,
    SpatialMagneticError, SpatialMagneticMaterial, SpatialMagnetostatics, SpatialRecoilMaterial,
    nonlinear::{
        SpatialNonlinearControls, SpatialNonlinearError, SpatialNonlinearMagnetostatics,
        SpatialNonlinearMaterial,
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
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 4000,
        max_dofs: 20000,
        max_triplets: 150000,
    }
}
fn controls() -> SpatialNonlinearControls {
    SpatialNonlinearControls {
        field: SpatialMagneticControls {
            cg: common::config(),
            flux_tolerance_wb: 2e-10,
            energy_tolerance_j: 2e-8,
        },
        relative_residual_tolerance: 1e-10,
        maximum_newton_iterations: 35,
        maximum_backtracks: 35,
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

fn potentials(
    problem: &SpatialNonlinearMagnetostatics,
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
fn ends(
    problem: &SpatialNonlinearMagnetostatics,
    axis: usize,
    length: f64,
    drop: f64,
) -> SpatialMagneticBoundary {
    SpatialMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .filter(|&f| {
                let p = problem.mesh().faces[f].map(|v| problem.positions_m()[v as usize][axis]);
                p.iter().all(|v| *v == 0.0) || p.iter().all(|v| *v == length)
            })
            .map(|face| MagneticPotentialFace {
                face,
                potential_a: problem.mesh().faces[face]
                    .map(|v| -drop * problem.positions_m()[v as usize][axis] / length),
            })
            .collect(),
        flux: vec![],
    }
}
fn homogeneous(n: usize) -> SpatialNonlinearMagnetostatics {
    let (p, t) = cube(n);
    let m = vec![SpatialNonlinearMaterial::Reversible(law()); t.len()];
    SpatialNonlinearMagnetostatics::new(p, t, m, limits(), &mut || false).unwrap()
}
fn uniform_flux(problem: &SpatialNonlinearMagnetostatics, b: [f64; 3]) -> SpatialMagneticBoundary {
    SpatialMagneticBoundary {
        potential: vec![],
        flux: problem
            .exterior_faces()
            .into_iter()
            .map(|face| MagneticFluxFace {
                face,
                outward_flux_density_t: dot(b, problem.exterior_normal(face).unwrap()),
            })
            .collect(),
    }
}

#[test]
fn original_linear_solution_and_flux_gauge_are_preserved() {
    let (p, t) = cube(2);
    let mu = 0.7;
    let source = material(mu);
    let nonlinear = SpatialNonlinearMagnetostatics::new(
        p.clone(),
        t.clone(),
        vec![SpatialNonlinearMaterial::Linear(source.clone()); t.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let linear = SpatialMagnetostatics::new(
        p,
        t.clone(),
        vec![SpatialMagneticMaterial::Linear(source); t.len()],
        limits(),
        &mut || false,
    )
    .unwrap();
    let h = [0.25, -0.5, 0.75];
    let bc = potentials(&nonlinear, |p| -dot(h, p));
    let expected = linear.solve(&bc, controls().field, &mut || false).unwrap();
    let initial = vec![0.0; nonlinear.positions_m().len()];
    let actual = nonlinear
        .solve(&bc, &initial, controls(), &mut || false)
        .unwrap();
    assert!(actual.report().newton_iterations > 0);
    let offset = 2.0_f64.powi(36);
    let shifted = nonlinear
        .solve(
            &potentials(&nonlinear, |p| offset - dot(h, p)),
            &vec![offset; initial.len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    let floating = nonlinear
        .solve(
            &uniform_flux(&nonlinear, h.map(|v| v * mu)),
            &vec![offset; initial.len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    for field in [&actual, &shifted, &floating] {
        for (a, b) in field.flux_density_t().iter().zip(expected.flux_density_t()) {
            for j in 0..3 {
                near(a[j], b[j]);
            }
        }
        near(field.stored_energy_j(), expected.stored_energy_j());
        near(field.coenergy_j(), expected.stored_energy_j());
        assert!(field.legendre_residual_j() < 1e-12);
        assert!(field.energy_pairing_residual_j().abs() < 1e-12);
    }
    assert_eq!(floating.gauge_vertices(), &[0]);
}

#[test]
fn saturated_layers_match_independent_flux_energy_and_rotated_solution() {
    for rotate in [false, true] {
        let (mut p, t) = grid(
            &[0.0, 0.5, 1.0, 1.5, 2.0],
            &[0.0, 0.5, 1.0],
            &[0.0, 0.5, 1.0],
        );
        let core: Vec<_> = t
            .iter()
            .map(|tet| tet.iter().map(|&v| p[v as usize][0] / 4.0).sum::<f64>() < 1.0)
            .collect();
        if rotate {
            for v in &mut p {
                *v = [v[2], v[0], v[1]];
            }
        }
        let axis = usize::from(rotate);
        let original = law();
        let materials = core
            .iter()
            .map(|&v| {
                if v {
                    SpatialNonlinearMaterial::Reversible(original.clone())
                } else {
                    SpatialNonlinearMaterial::Linear(material(1.0))
                }
            })
            .collect();
        let problem =
            SpatialNonlinearMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap();
        let initial: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| -2.0 * p[axis])
            .collect();
        let field = problem
            .solve(
                &ends(&problem, axis, 2.0, 4.0),
                &initial,
                controls(),
                &mut || false,
            )
            .unwrap();
        // H_core(B)=2B-1.5 and H_air(B)=B: H_core+H_air=4.
        let b = 11.0 / 6.0;
        let h_core = 13.0 / 6.0;
        let e_core = 0.25 + 0.5 * (b - 1.0) + (b - 1.0) * (b - 1.0);
        let energy = e_core + 0.5 * b * b;
        for (cell, &is_core) in core.iter().enumerate() {
            for j in 0..3 {
                near(
                    field.flux_density_t()[cell][j],
                    if j == axis { b } else { 0.0 },
                );
                near(
                    field.field_strength_a_m()[cell][j],
                    if j == axis {
                        if is_core { h_core } else { b }
                    } else {
                        0.0
                    },
                );
            }
        }
        near(field.stored_energy_j(), energy);
        near(field.coenergy_j(), 4.0 * b - energy);
        assert!(
            (field.stored_energy_j() - 2.0 * b).abs() > 0.1,
            "saturated energy is not B.H/2"
        );
        assert!(
            field
                .interface_flux_jumps_wb()
                .iter()
                .all(|(_, v)| v.abs() < 1e-10)
        );
        assert!(field.report().newton_iterations > 0);
        let left: f64 = problem
            .positions_m()
            .iter()
            .zip(field.nodal_boundary_flux_wb())
            .filter(|(p, _)| p[axis] == 0.0)
            .map(|(_, v)| v)
            .sum();
        near(left, -b);
    }
}

#[test]
fn stationary_recoil_bias_saturates_the_adjacent_core_with_correct_dual_energy() {
    let (p, t) = grid(&[0.0, 1.0, 2.0, 3.0], &[0.0, 0.5, 1.0], &[0.0, 0.5, 1.0]);
    let regions: Vec<_> = t
        .iter()
        .map(|tet| {
            tet.iter()
                .map(|&v| p[v as usize][0] / 4.0)
                .sum::<f64>()
                .floor() as usize
        })
        .collect();
    let magnet = SpatialRecoilMaterial::new(recoil(1.0, 5.0, 10.0), [1.0, 0.0, 0.0]).unwrap();
    let curve = law();
    let materials = regions
        .iter()
        .map(|r| match r {
            0 => SpatialNonlinearMaterial::Recoil(magnet.clone()),
            1 => SpatialNonlinearMaterial::Reversible(curve.clone()),
            _ => SpatialNonlinearMaterial::Linear(material(1.0)),
        })
        .collect();
    let problem =
        SpatialNonlinearMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap();
    let field = problem
        .solve(
            &ends(&problem, 0, 3.0, 0.0),
            &vec![0.0; problem.positions_m().len()],
            controls(),
            &mut || false,
        )
        .unwrap();
    // (B-5) + (2B-1.5) + B = 0 with B above the saturation knee.
    let b = 1.625;
    for (cell, &region) in regions.iter().enumerate() {
        near(field.flux_density_t()[cell][0], b);
        near(
            field.field_strength_a_m()[cell][0],
            match region {
                0 => -3.375,
                1 => 1.75,
                _ => 1.625,
            },
        );
        assert_eq!(field.recoil_operating_points()[cell].is_some(), region == 0);
    }
    near(field.stored_energy_j(), 7.96875);
    near(field.coenergy_j(), -7.96875);
    assert!(field.energy_pairing_residual_j().abs() < 1e-11);
    assert!(
        field
            .interface_flux_jumps_wb()
            .iter()
            .all(|(_, jump)| jump.abs() < 1e-10)
    );
}

fn radial_potential(p: [f64; 3]) -> f64 {
    let r = dot(p, p).sqrt();
    if r >= 6.0_f64.sqrt() {
        3.0 / r
    } else {
        12.0 / r + 1.5 * r - 3.0 * 6.0_f64.sqrt()
    }
}
fn radial_fields(p: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let r = dot(p, p).sqrt();
    let b = 6.0 / (r * r);
    let h = if b <= 1.0 { 0.5 * b } else { 2.0 * b - 1.5 };
    (p.map(|v| h * v / r), p.map(|v| b * v / r))
}
#[test]
fn nonlinear_current_free_radial_field_converges_under_mesh_refinement() {
    // B=6 r/|r|^3 is divergence-free away from the origin. H(B) is radial
    // and curl-free; the explicit piecewise potential crosses the source knee.
    let mut errors = Vec::new();
    for n in [2, 4, 8] {
        let q: Vec<_> = (0..=n).map(|i| 1.0 + i as f64 / n as f64).collect();
        let (p, t) = grid(&q, &q, &q);
        let materials = vec![SpatialNonlinearMaterial::Reversible(law()); t.len()];
        let problem =
            SpatialNonlinearMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap();
        let initial: Vec<_> = problem
            .positions_m()
            .iter()
            .copied()
            .map(radial_potential)
            .collect();
        let field = problem
            .solve(
                &potentials(&problem, radial_potential),
                &initial,
                controls(),
                &mut || false,
            )
            .unwrap();
        let alpha = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
        let beta = (5.0 - 5.0_f64.sqrt()) / 20.0;
        let volume = 1.0 / (6.0 * (n * n * n) as f64);
        let mut error = [0.0; 2];
        for (cell, tet) in problem.mesh().tets.iter().enumerate() {
            for large in 0..4 {
                let mut p = [0.0; 3];
                for (i, &v) in tet.iter().enumerate() {
                    let w = if i == large { alpha } else { beta };
                    for j in 0..3 {
                        p[j] += w * problem.positions_m()[v as usize][j];
                    }
                }
                let (h, b) = radial_fields(p);
                for j in 0..3 {
                    error[0] +=
                        0.25 * volume * (field.field_strength_a_m()[cell][j] - h[j]).powi(2);
                    error[1] += 0.25 * volume * (field.flux_density_t()[cell][j] - b[j]).powi(2);
                }
            }
        }
        errors.push(error.map(f64::sqrt));
        assert!(field.report().maximum_weak_residual_wb <= controls().field.flux_tolerance_wb);
    }
    for pair in errors.windows(2) {
        for j in 0..2 {
            let order = (pair[0][j] / pair[1][j]).log2();
            assert!(
                (0.8..1.2).contains(&order),
                "field{j} order={order}, errors={errors:?}"
            );
        }
    }
}

#[test]
fn inverse_law_tangent_and_legendre_derivatives_hold_in_three_dimensions() {
    let problem = homogeneous(1);
    let direction = [0.6, 0.0, 0.8];
    // The saturated segment has B=(H+1.5)/2 and energy
    // 1/4 +(B-1)/2 +(B-1)^2; dW*/dH=B.
    let evaluate = |h: f64| {
        let bc = potentials(&problem, |p| -h * dot(direction, p));
        let initial: Vec<_> = problem
            .positions_m()
            .iter()
            .map(|p| -h * dot(direction, *p))
            .collect();
        problem
            .solve(&bc, &initial, controls(), &mut || false)
            .unwrap()
    };
    let h = 1.75;
    let delta = 1e-5;
    let center = evaluate(h);
    let low = evaluate(h - delta);
    let high = evaluate(h + delta);
    let b = 1.625;
    near(dot(center.flux_density_t()[0], direction), b);
    near((high.coenergy_j() - low.coenergy_j()) / (2.0 * delta), b);
    near(
        (dot(high.flux_density_t()[0], direction) - dot(low.flux_density_t()[0], direction))
            / (2.0 * delta),
        0.5,
    );
    let e = 0.25 + 0.5 * (b - 1.0) + (b - 1.0) * (b - 1.0);
    near(center.stored_energy_j(), e);
    near(center.coenergy_j(), b * h - e);
    near(center.stored_energy_j() + center.coenergy_j(), b * h);
}

#[test]
fn source_domain_flux_work_and_cancellation_refuse_without_a_partial_field() {
    let problem = homogeneous(2);
    let initial = vec![0.0; problem.positions_m().len()];
    let bc = uniform_flux(&problem, [1.6, 0.0, 0.0]);
    // A loose flux solve may not convert its interior residual into external
    // supply: W+W*=.08 J, authored trace work=-.2 J for this trial.
    let loose_bc = uniform_flux(&problem, [1.0, 0.0, 0.0]);
    let loose_initial: Vec<_> = problem.positions_m().iter().map(|p| -0.2 * p[0]).collect();
    let mut loose = controls();
    loose.field.flux_tolerance_wb = 10.0;
    loose.field.energy_tolerance_j = 1e-9;
    loose.maximum_newton_iterations = 0;
    assert!(matches!(
        problem.solve(&loose_bc, &loose_initial, loose, &mut || false),
        Err(SpatialNonlinearError::Spatial(
            SpatialMagneticError::Invalid(_)
        ))
    ));
    let mut zero = controls();
    zero.maximum_newton_iterations = 0;
    assert!(matches!(
        problem.solve(&bc, &initial, zero, &mut || false),
        Err(SpatialNonlinearError::Nonconvergence(_))
    ));
    let mut krylov = controls();
    krylov.field.cg.max_iterations = 0;
    assert!(matches!(
        problem.solve(&bc, &initial, krylov, &mut || false),
        Err(SpatialNonlinearError::Spatial(
            SpatialMagneticError::Nonconvergence(_)
        ))
    ));
    let bad = potentials(&problem, |p| -11.0 * p[0]);
    let bad_initial: Vec<_> = problem.positions_m().iter().map(|p| -11.0 * p[0]).collect();
    assert!(matches!(
        problem.solve(&bad, &bad_initial, controls(), &mut || false),
        Err(SpatialNonlinearError::FieldRange { .. })
    ));
    let mut flux = bc.clone();
    flux.flux[0].outward_flux_density_t += 1.0;
    assert!(matches!(
        problem.solve(&flux, &initial, controls(), &mut || false),
        Err(SpatialNonlinearError::Spatial(
            SpatialMagneticError::Invalid(_)
        ))
    ));
    assert!(matches!(
        problem.solve(&bc, &initial, controls(), &mut || true),
        Err(SpatialNonlinearError::Spatial(
            SpatialMagneticError::Cancelled
        ))
    ));
    let mut checkpoints = 0;
    let accepted = problem
        .solve(&bc, &initial, controls(), &mut || {
            checkpoints += 1;
            false
        })
        .unwrap();
    let calls = checkpoints;
    let mut seen = 0;
    assert!(matches!(
        problem.solve(&bc, &initial, controls(), &mut || {
            seen += 1;
            seen >= calls - 1
        }),
        Err(SpatialNonlinearError::Spatial(
            SpatialMagneticError::Cancelled
        ))
    ));
    for h in accepted.field_strength_a_m() {
        near(h[0], 1.7);
    }
    let again = problem
        .solve(&bc, &initial, controls(), &mut || false)
        .unwrap();
    assert_eq!(accepted.potential_a(), again.potential_a());
    assert_eq!(accepted.flux_density_t(), again.flux_density_t());
}
