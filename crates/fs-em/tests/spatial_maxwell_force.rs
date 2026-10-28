//! Maxwell force and torque against an independent uniformly magnetized cuboid.
//! Manufactured material sources test the discretization, not material quality.
mod common;
use common::material;
use fs_em::spatial::force::SpatialMaxwellForceLimits;
use fs_em::spatial::nonlinear::{
    SpatialNonlinearControls, SpatialNonlinearMagnetostatics, SpatialNonlinearMaterial,
};
use fs_em::spatial::{
    MagneticPotentialFace, SpatialMagneticBoundary, SpatialMagneticControls, SpatialMagneticError,
    SpatialMagneticMaterial, SpatialMagnetostatics, SpatialRecoilMaterial,
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_ga::{Se3, So3, So3Tangent, Vec3, Wrench};
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
use fs_qty::{QuantitySpec, Temperature};

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
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 30000,
        max_dofs: 100000,
        max_triplets: 1_100_000,
    }
}
fn force_limits() -> SpatialMaxwellForceLimits {
    SpatialMaxwellForceLimits {
        max_cells: 30000,
        max_faces: 100000,
    }
}
fn controls() -> SpatialMagneticControls {
    SpatialMagneticControls {
        cg: fs_solver::CheckedCgConfig {
            tolerance: 2e-11,
            max_iterations: 1000,
            max_corrections: 2,
        },
        flux_tolerance_wb: 1e-9,
        energy_tolerance_j: 1e-7,
    }
}
fn mesh(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 4]>) {
    let q: Vec<_> = (0..=n).map(|i| -2.0 + 4.0 * i as f64 / n as f64).collect();
    grid(&q, &q, &q)
}
fn center(positions: &[[f64; 3]], tet: [u32; 4]) -> [f64; 3] {
    std::array::from_fn(|axis| tet.iter().map(|&v| positions[v as usize][axis] / 4.0).sum())
}
fn region(problem: &SpatialMagnetostatics, radius: f64) -> Vec<usize> {
    problem
        .mesh()
        .tets
        .iter()
        .enumerate()
        .filter_map(|(i, &tet)| {
            center(problem.positions_m(), tet)
                .iter()
                .all(|x| x.abs() < radius)
                .then_some(i)
        })
        .collect()
}
fn vector_error(a: Vec3, b: Vec3) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
fn assert_wrench(actual: Wrench, expected: Wrench, tolerance: f64) {
    assert!(
        vector_error(actual.force, expected.force) <= tolerance,
        "force {:?} != {:?}",
        actual.force,
        expected.force
    );
    assert!(
        vector_error(actual.torque, expected.torque) <= tolerance,
        "torque {:?} != {:?}",
        actual.torque,
        expected.torque
    );
}
fn linear(n: usize, mu: f64) -> SpatialMagnetostatics {
    let (p, t) = mesh(n);
    let materials = vec![SpatialMagneticMaterial::Linear(material(mu)); t.len()];
    SpatialMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap()
}

#[test]
fn uniform_field_has_no_net_force_or_torque_on_a_closed_interior_volume() {
    let problem = linear(4, 2.0);
    let field = problem
        .solve(
            &potentials(&problem, |p| -2.0 * p[0] + p[1] - 3.0 * p[2]),
            controls(),
            &mut || false,
        )
        .unwrap();
    let selected = region(&problem, 1.0);
    let first = field
        .maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert_wrench(first, Wrench::default(), 2e-9);
    let mut reordered = selected.clone();
    reordered.reverse();
    let second = field
        .maxwell_wrench(
            &reordered,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(first, second);
    let placement = Se3::from_parts(
        So3::exp(So3Tangent {
            angular: Vec3::new(0.3, -0.2, 0.4),
        })
        .unwrap(),
        Vec3::new(4.0, -3.0, 2.0),
    )
    .unwrap();
    let placed = field
        .maxwell_wrench(
            &selected,
            &placement,
            [2.0, -1.0, 0.5],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert_wrench(placed, Wrench::default(), 2e-8);
}

// Independent scalar potential of a UNIT uniformly magnetized cube [-.5,.5]^3.
// Its only magnetic-charge faces are z=+/- .5. The exterior potential is
// (1/4pi) integral_surface (M.n)/|r-r'| dS. An 8x8 Gauss rule is smooth and
// rapidly convergent on the remote box boundary; it does not use FE gradients,
// Hodge entries, face normals, or Maxwell force quadrature from production.
fn cube_magnet_potential(p: [f64; 3]) -> f64 {
    const Q: [f64; 8] = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.525532409916329,
        -0.1834346424956498,
        0.1834346424956498,
        0.525532409916329,
        0.7966664774136267,
        0.9602898564975363,
    ];
    const W: [f64; 8] = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.362683783378362,
        0.362683783378362,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ];
    let mut potential = 0.0;
    for i in 0..8 {
        for j in 0..8 {
            let x = p[0] - 0.5 * Q[i];
            let y = p[1] - 0.5 * Q[j];
            let top = x.hypot(y).hypot(p[2] - 0.5);
            let bottom = x.hypot(y).hypot(p[2] + 0.5);
            potential += 0.25 * W[i] * W[j] * (top.recip() - bottom.recip());
        }
    }
    potential / (4.0 * core::f64::consts::PI)
}
fn magnet(n: usize) -> SpatialMagnetostatics {
    let (p, t) = mesh(n);
    let air = material(1.0);
    let magnet = SpatialRecoilMaterial::new(recoil(1.0, 1.0, 20.0), [0.0, 0.0, 1.0]).unwrap();
    let materials = t
        .iter()
        .map(|&tet| {
            if center(&p, tet).iter().all(|x| x.abs() < 0.5) {
                SpatialMagneticMaterial::Recoil(magnet.clone())
            } else {
                SpatialMagneticMaterial::Linear(air.clone())
            }
        })
        .collect();
    SpatialMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap()
}

#[test]
fn remanent_cuboid_in_gradient_field_recovers_independent_nonzero_force_and_torque() {
    let expected = Wrench::new(Vec3::new(0.0, 0.4, 0.0), Vec3::new(0.7, 0.0, 0.0));
    let mut previous = None;
    for n in [8, 16] {
        let problem = magnet(n);
        // H_ext=(.4+.7*z,0,.7*x) is curl-free and divergence-free.
        // With mu=1 and Br=(0,0,1), the unit-volume magnet has m=(0,0,1).
        // Independently F=(m.grad)B_ext=(.7,0,0), tau=m x B_ext(0)=(0,.4,0).
        let boundary = potentials(&problem, |p| {
            cube_magnet_potential(p) - 0.4 * p[0] - 0.7 * p[0] * p[2]
        });
        let field = problem.solve(&boundary, controls(), &mut || false).unwrap();
        let selected = region(&problem, 1.0);
        let actual = field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let error = vector_error(actual.force, expected.force)
            + vector_error(actual.torque, expected.torque);
        eprintln!(
            "n={n} force={:?} torque={:?} error={error:e}",
            actual.force, actual.torque
        );
        let outer = field
            .maxwell_wrench(
                &region(&problem, 1.5),
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        eprintln!(
            "n={n} outer force={:?} outer torque={:?}",
            outer.force, outer.torque
        );
        if let Some((coarse, coarse_outer, old_error)) = previous {
            let ratio: f64 = old_error / error;
            assert!(
                ratio > 1.8 && ratio < 2.2,
                "expected first-order exterior P1 stress convergence, ratio {ratio}"
            );
            // A TEST-ONLY first-order extrapolation checks the independently
            // known continuum limit. The public API returns the raw surface
            // integral, with no extrapolation or certified-error claim.
            let extrapolate = |fine: Wrench, coarse: Wrench| {
                let f = Vec3::new(
                    2.0 * fine.force.x - coarse.force.x,
                    2.0 * fine.force.y - coarse.force.y,
                    2.0 * fine.force.z - coarse.force.z,
                );
                let t = Vec3::new(
                    2.0 * fine.torque.x - coarse.torque.x,
                    2.0 * fine.torque.y - coarse.torque.y,
                    2.0 * fine.torque.z - coarse.torque.z,
                );
                Wrench::new(t, f)
            };
            assert_wrench(extrapolate(actual, coarse), expected, 0.02);
            assert_wrench(extrapolate(outer, coarse_outer), expected, 0.02);
        }
        previous = Some((actual, outer, error));
        assert!(actual.force.x > 0.0 && actual.torque.y > 0.0);
        // Compare origin/frame behavior on the SAME nonzero discrete load.
        let rotation = So3::exp(So3Tangent {
            angular: Vec3::new(0.3, -0.2, 0.4),
        })
        .unwrap();
        let placement = Se3::from_parts(rotation, Vec3::new(4.0, -3.0, 2.0)).unwrap();
        let origin = Vec3::new(0.4, -0.7, 0.2);
        let transformed = actual.transform_by(&placement);
        let shift = origin.cross(transformed.force);
        let expected_frame = Wrench::new(
            Vec3::new(
                transformed.torque.x - shift.x,
                transformed.torque.y - shift.y,
                transformed.torque.z - shift.z,
            ),
            transformed.force,
        );
        let moved = field
            .maxwell_wrench(
                &selected,
                &placement,
                [origin.x, origin.y, origin.z],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        assert_wrench(moved, expected_frame, 2e-11);
        let shifted = field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [1.0, 2.0, 3.0],
                force_limits(),
                &mut || false,
            )
            .unwrap();
        let shift = Vec3::new(1.0, 2.0, 3.0).cross(actual.force);
        assert_wrench(
            shifted,
            Wrench::new(
                Vec3::new(
                    actual.torque.x - shift.x,
                    actual.torque.y - shift.y,
                    actual.torque.z - shift.z,
                ),
                actual.force,
            ),
            2e-11,
        );
        // The magnet surface itself is not an admitted homogeneous exterior.
        assert!(matches!(
            field.maxwell_wrench(
                &region(&problem, 0.5),
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false
            ),
            Err(SpatialMagneticError::Invalid(_))
        ));
    }
}

#[test]
fn surface_refuses_different_original_sources_even_with_identical_permeability() {
    let base = linear(4, 1.0);
    let selected = region(&base, 1.0);
    let card = common::source(
        PropertyValue::Scalar {
            value: 1.0,
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
        ValidityDomain::unconstrained(),
    );
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
        .unwrap();
    let distinct = common::resolve(&card, &point);
    let mut sources = base.materials().to_vec();
    sources[selected[0]] = SpatialMagneticMaterial::Linear(distinct);
    let problem = SpatialMagnetostatics::new(
        base.positions_m().to_vec(),
        base.mesh().tets.clone(),
        sources,
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = problem
        .solve(&potentials(&problem, |p| -p[0]), controls(), &mut || false)
        .unwrap();
    assert!(matches!(
        field.maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false
        ),
        Err(SpatialMagneticError::Invalid(_))
    ));
}

#[test]
fn invalid_open_boundary_regions_budgets_and_cancellation_publish_no_wrench() {
    let problem = linear(4, 1.0);
    let field = problem
        .solve(
            &potentials(&problem, |p| -p[0] - 0.2 * p[1]),
            controls(),
            &mut || false,
        )
        .unwrap();
    let selected = region(&problem, 1.0);
    let mut repeated = selected.clone();
    repeated.push(selected[0]);
    for cells in [vec![], vec![0], vec![problem.mesh().tets.len()], repeated] {
        assert!(matches!(
            field.maxwell_wrench(
                &cells,
                &Se3::identity(),
                [0.0; 3],
                force_limits(),
                &mut || false
            ),
            Err(SpatialMagneticError::Invalid(_))
        ));
    }
    let limit = SpatialMaxwellForceLimits {
        max_cells: problem.mesh().tets.len() - 1,
        max_faces: 100000,
    };
    assert!(
        field
            .maxwell_wrench(&selected, &Se3::identity(), [0.0; 3], limit, &mut || false)
            .is_err()
    );
    assert!(
        field
            .maxwell_wrench(
                &selected,
                &Se3::identity(),
                [f64::NAN, 0.0, 0.0],
                force_limits(),
                &mut || false
            )
            .is_err()
    );
    let mut polls = 0;
    let cancel_after = selected.len() + problem.mesh().faces.len() / 2;
    assert!(matches!(
        field.maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || {
                polls += 1;
                polls > cancel_after
            }
        ),
        Err(SpatialMagneticError::Cancelled)
    ));
    let original = field
        .maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert_wrench(original, Wrench::default(), 1e-9);
}

fn saturation_law() -> ReversibleMagneticMaterial {
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
                source: "manufactured saturation for 3D exterior force".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic soft magnet".into(),
            phase: "solid".into(),
            process: "reversible force fixture".into(),
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

#[test]
fn saturated_core_uses_linear_exterior_traction_and_refuses_its_material_interface() {
    let q = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0];
    let (p, t) = grid(&q, &q, &q);
    let law = saturation_law();
    let air = material(1.0);
    let selected: Vec<_> = t
        .iter()
        .enumerate()
        .filter_map(|(cell, &tet)| {
            center(&p, tet)
                .iter()
                .all(|x| x.abs() < 1.0)
                .then_some(cell)
        })
        .collect();
    let core: Vec<_> = t
        .iter()
        .enumerate()
        .filter_map(|(cell, &tet)| {
            center(&p, tet)
                .iter()
                .all(|x| x.abs() < 0.5)
                .then_some(cell)
        })
        .collect();
    assert_eq!(core.len(), 6);
    let materials = t
        .iter()
        .enumerate()
        .map(|(cell, _)| {
            if core.contains(&cell) {
                SpatialNonlinearMaterial::Reversible(law.clone())
            } else {
                SpatialNonlinearMaterial::Linear(air.clone())
            }
        })
        .collect();
    let problem =
        SpatialNonlinearMagnetostatics::new(p, t, materials, limits(), &mut || false).unwrap();
    // A centered cubic core in a uniform diagonal applied field is invariant
    // under inversion and coordinate permutations. Its net force and torque
    // therefore vanish even after the original nonlinear core law saturates.
    let boundary = SpatialMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| MagneticPotentialFace {
                face,
                potential_a: problem.mesh().faces[face]
                    .map(|v| -problem.positions_m()[v as usize].iter().sum::<f64>()),
            })
            .collect(),
        flux: vec![],
    };
    let settings = SpatialNonlinearControls {
        field: controls(),
        relative_residual_tolerance: 1e-10,
        maximum_newton_iterations: 35,
        maximum_backtracks: 35,
    };
    let initial: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| -p.iter().sum::<f64>())
        .collect();
    let field = problem
        .solve(&boundary, &initial, settings, &mut || false)
        .unwrap();
    for &cell in &core {
        let b = field.flux_density_t()[cell];
        let magnitude = b[0].hypot(b[1]).hypot(b[2]);
        assert!(
            magnitude > 1.0 && magnitude < 4.0,
            "core must enter an admitted saturated source segment, B={magnitude}"
        );
    }
    let wrench = field
        .maxwell_wrench(
            &selected,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        )
        .unwrap();
    assert_wrench(wrench, Wrench::default(), 2e-8);
    // Sampling the core/air interface would require a different constitutive
    // stress law; using exterior air alone cannot authorize that surface.
    assert!(matches!(
        field.maxwell_wrench(
            &core,
            &Se3::identity(),
            [0.0; 3],
            force_limits(),
            &mut || false,
        ),
        Err(SpatialMagneticError::Invalid(
            "3D Maxwell surface requires induced-linear material on both sides"
        ))
    ));
}

#[path = "spatial_maxwell_force/eggshell.rs"]
mod eggshell;
