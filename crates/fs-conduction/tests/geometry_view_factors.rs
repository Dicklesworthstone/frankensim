//! G0/G1 visibility and heat-balance evidence; G4 refusal and G5 replay checks.
//! These finite geometry fixtures establish estimates, not certified bounds.

mod support;

use fs_conduction::fixtures::box_grid;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::radiation::view_factors::{
    DiffuseViewFactorConfig, ViewFactorProjection, estimate_diffuse_view_factors,
};
use fs_conduction::{
    ConductionError, EMISSIVITY_DIMS, GrayDiffuseEnclosure, RadiationSurface,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ViewFactorEvidence,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;
use support::{with_cancelled_cx, with_cx};

fn emissivity() -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: 0.8,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "G1 geometry-view-factor fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("fixture emissivity claim");
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture-alloy".to_string(),
            phase: "solid".to_string(),
            process: "diffuse-gray".to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("fixture card");
    SurfaceEmissivity::from_card(
        "wall-finish",
        &card,
        350.0,
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("fixture emissivity")
}

fn config() -> DiffuseViewFactorConfig {
    DiffuseViewFactorConfig {
        seed: 0xC0FF_EE42,
        samples_per_replicate: 1024,
        replicates: 4,
        max_triangles: 128,
        max_ray_triangle_tests: 16_000_000,
    }
}

fn projection() -> ViewFactorProjection {
    ViewFactorProjection {
        max_iterations: 1024,
        row_sum_tolerance: 1.0e-12,
        max_abs_adjustment: 0.03,
    }
}

/// Six separate solids surround a 2 m by 1 m by 1 m void. All normals are
/// extracted from real tetrahedra; the selected outward normals point into it.
fn cavity(blocked: bool) -> (ConductionMesh, Vec<RadiationSurface>) {
    cavity_geometry(blocked, false)
}

fn cavity_geometry(blocked: bool, open: bool) -> (ConductionMesh, Vec<RadiationSurface>) {
    let size = [2.0, 1.0, 1.0];
    let mut positions = Vec::new();
    let mut tets = Vec::new();
    let mut append_box = |origin: [f64; 3], extent: [f64; 3]| {
        let (complex, points) = box_grid([1, 1, 1], extent);
        let offset = u32::try_from(positions.len()).expect("small fixture");
        positions.extend(
            points
                .into_iter()
                .map(|p| [p[0] + origin[0], p[1] + origin[1], p[2] + origin[2]]),
        );
        tets.extend(complex.tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    };
    for axis in 0..3 {
        for side in 0..2 {
            if open && axis == 0 && side == 1 {
                continue;
            }
            let mut origin = [0.0; 3];
            origin[axis] = if side == 0 { -0.125 } else { size[axis] };
            let mut extent = size;
            extent[axis] = 0.125;
            append_box(origin, extent);
        }
    }
    if blocked {
        append_box([0.9375, 0.0, 0.0], [0.125, 1.0, 1.0]);
    }
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions)
        .expect("six panels around a rectangular cavity");
    let finish = emissivity();
    let mut surfaces = Vec::new();
    for (axis, &length) in size.iter().enumerate() {
        for side in 0..2 {
            if open && axis == 0 && side == 1 {
                continue;
            }
            let coordinate = if side == 0 { 0.0 } else { length };
            let normal = if side == 0 { 1.0 } else { -1.0 };
            surfaces.push(
                RadiationSurface::new(
                    &mesh,
                    format!("wall-{axis}-{side}"),
                    |face| face.centroid[axis] == coordinate && face.outward_normal[axis] == normal,
                    finish.clone(),
                )
                .expect("interior wall"),
            );
        }
    }
    (mesh, surfaces)
}

#[test]
fn g1_finite_parallel_square_factor_and_unequal_area_energy_closure() {
    let (mesh, surfaces) = cavity(false);
    with_cx(|cx| {
        let estimate = estimate_diffuse_view_factors(
            cx,
            &mesh,
            &surfaces,
            DiffuseViewFactorConfig {
                samples_per_replicate: 4096,
                ..config()
            },
        )
        .expect("geometry estimate");
        assert_eq!(estimate.areas_m2(), &[1.0, 1.0, 2.0, 2.0, 2.0, 2.0]);
        assert_eq!(estimate.samples_per_surface(), 16_384);
        assert!(estimate.ray_triangle_tests() > 0);
        assert!(estimate.ray_triangle_tests() <= config().max_ray_triangle_tests);

        // Closed-form parallel-rectangle integral with X = Y = width/gap = 1/2:
        // 2/(pi X Y) [log sqrt((1+X²)(1+Y²)/(1+X²+Y²))
        // + X sqrt(1+Y²) atan(X/sqrt(1+Y²)) + (X <-> Y)
        // - X atan(X) - Y atan(Y)]. Independent of the ray estimator.
        const REFERENCE: f64 = 0.068_589_588_818_552_6;
        for (i, j) in [(0, 1), (1, 0)] {
            assert!(
                (estimate.factors()[i][j] - REFERENCE).abs() < 0.006,
                "finite plate estimate {} differs from {REFERENCE}",
                estimate.factors()[i][j]
            );
        }
        for i in 0..surfaces.len() {
            assert_eq!(
                estimate.hit_counts()[i].iter().sum::<u64>(),
                estimate.samples_per_surface()
            );
            assert_eq!(estimate.escaped_counts()[i], 0);
            assert_eq!(estimate.unassigned_counts()[i], 0);
            assert_eq!(estimate.backface_counts()[i], 0);
            assert_eq!(estimate.factors()[i][i], 0.0);
            for j in 0..surfaces.len() {
                assert_eq!(
                    estimate.factors()[i][j],
                    estimate.hit_counts()[i][j] as f64 / estimate.samples_per_surface() as f64
                );
                assert!(estimate.standard_errors()[i][j].is_finite());
                assert!(estimate.standard_errors()[i][j] >= 0.0);
            }
        }
        assert!(
            estimate
                .project_closed(
                    cx,
                    ViewFactorProjection {
                        max_abs_adjustment: 1.0e-15,
                        ..projection()
                    }
                )
                .is_err(),
            "a reciprocal repair cannot silently exceed its adjustment budget"
        );
        let projected = estimate
            .project_closed(cx, projection())
            .expect("closed projection");
        assert!(projected.max_abs_adjustment() <= projection().max_abs_adjustment);
        assert!(projected.iterations() <= projection().max_iterations);
        let matrix = projected.matrix();
        assert!(matches!(
            matrix.evidence(),
            ViewFactorEvidence::GeometryQmc { .. }
        ));
        for i in 0..surfaces.len() {
            assert!((matrix.row_sums()[i] - 1.0).abs() <= 2.0e-12);
            for j in 0..surfaces.len() {
                let exchange = matrix.areas_m2()[i] * matrix.factors()[i][j];
                let reverse = matrix.areas_m2()[j] * matrix.factors()[j][i];
                assert!((exchange - reverse).abs() <= 1.0e-12);
            }
        }
        let enclosure = GrayDiffuseEnclosure::new(surfaces, matrix.clone()).expect("enclosure");
        let report = enclosure
            .solve(cx, &[450.0, 300.0, 350.0, 325.0, 375.0, 400.0])
            .expect("radiosity from computed geometry");
        assert!(report.net_outward_heat_w[0] > 0.0);
        assert!(report.net_outward_heat_w[1] < 0.0);
        assert!(report.relative_energy_closure() < 1.0e-10);
    });
}

#[test]
fn g5_replay_and_surface_reorder_preserve_named_random_streams() {
    let (mesh, mut surfaces) = cavity(false);
    with_cx(|cx| {
        let first = estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()).unwrap();
        let replay = estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()).unwrap();
        assert_eq!(first.identity(), replay.identity());
        assert_eq!(first.hit_counts(), replay.hit_counts());
        assert_eq!(first.standard_errors(), replay.standard_errors());
        let changed_seed = estimate_diffuse_view_factors(
            cx,
            &mesh,
            &surfaces,
            DiffuseViewFactorConfig {
                seed: config().seed + 1,
                ..config()
            },
        )
        .unwrap();
        assert_ne!(first.identity(), changed_seed.identity());
        assert_ne!(first.hit_counts(), changed_seed.hit_counts());
        surfaces.reverse();
        let reversed = estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()).unwrap();
        let n = surfaces.len();
        for i in 0..n {
            assert_eq!(first.names()[i], reversed.names()[n - 1 - i]);
            for j in 0..n {
                assert_eq!(
                    first.hit_counts()[i][j],
                    reversed.hit_counts()[n - 1 - i][n - 1 - j]
                );
                assert_eq!(
                    first.standard_errors()[i][j],
                    reversed.standard_errors()[n - 1 - i][n - 1 - j]
                );
            }
        }
    });
}

#[test]
fn unassigned_faces_occlude_targets_and_cannot_be_projected_closed() {
    let (mesh, surfaces) = cavity(false);
    with_cx(|cx| {
        let partial = estimate_diffuse_view_factors(cx, &mesh, &surfaces[..2], config()).unwrap();
        assert!(partial.hit_counts()[0][1] > 0);
        assert!(partial.unassigned_counts()[0] > 0);
        assert_eq!(partial.escaped_counts(), &[0, 0]);
        assert!(partial.project_closed(cx, projection()).is_err());
        assert!(
            partial.project_open(cx, projection()).is_err(),
            "opaque unnamed walls are not a black external reservoir"
        );
        let (blocked_mesh, blocked_surfaces) = cavity(true);
        let blocked =
            estimate_diffuse_view_factors(cx, &blocked_mesh, &blocked_surfaces[..2], config())
                .unwrap();
        for i in 0..2 {
            assert_eq!(blocked.hit_counts()[i].iter().sum::<u64>(), 0);
            assert_eq!(
                blocked.unassigned_counts()[i],
                blocked.samples_per_surface()
            );
            assert_eq!(blocked.escaped_counts()[i], 0);
            assert_eq!(blocked.backface_counts()[i], 0);
        }
        assert!(blocked.project_closed(cx, projection()).is_err());
        assert!(blocked.project_open(cx, projection()).is_err());
    });
}

#[test]
fn outward_open_surfaces_report_escape_without_renormalizing() {
    let (complex, points) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, points).unwrap();
    let surfaces = [0.0, 1.0].map(|x| {
        RadiationSurface::new(
            &mesh,
            format!("outer-{x}"),
            |face| face.centroid[0] == x,
            emissivity(),
        )
        .unwrap()
    });
    with_cx(|cx| {
        let estimate = estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()).unwrap();
        for i in 0..2 {
            assert_eq!(estimate.escaped_counts()[i], estimate.samples_per_surface());
            assert_eq!(estimate.hit_counts()[i], [0, 0]);
            assert_eq!(estimate.factors()[i], [0.0, 0.0]);
            assert_eq!(estimate.unassigned_counts()[i], 0);
            assert_eq!(estimate.backface_counts()[i], 0);
        }
        assert!(estimate.project_closed(cx, projection()).is_err());
        let open = estimate
            .project_open(cx, projection())
            .expect("escaped rays reach the bath");
        assert_eq!(open.surface_factors(), &[vec![0.0, 0.0], vec![0.0, 0.0]]);
        assert_eq!(open.ambient_factors(), &[1.0, 1.0]);
        assert_eq!(open.max_abs_adjustment(), 0.0);
    });
}

#[test]
fn g1_open_cavity_retains_geometry_escape_and_closes_surface_plus_bath_heat() {
    use fs_conduction::radiation::open::OpenDiffuseEnclosure;
    let (mesh, surfaces) = cavity_geometry(false, true);
    with_cx(|cx| {
        let estimate = estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()).unwrap();
        assert!(estimate.escaped_counts().iter().all(|&n| n > 0));
        assert!(estimate.unassigned_counts().iter().all(|&n| n == 0));
        assert!(estimate.project_closed(cx, projection()).is_err());
        let open = estimate
            .project_open(cx, projection())
            .expect("open geometry projection");
        for (i, &ambient) in open.ambient_factors().iter().enumerate() {
            let raw_escape =
                estimate.escaped_counts()[i] as f64 / estimate.samples_per_surface() as f64;
            assert!(ambient > 0.0 && ambient < 1.0);
            assert!((ambient - raw_escape).abs() <= open.max_abs_adjustment());
            assert!((open.surface_factors()[i].iter().sum::<f64>() + ambient - 1.0).abs() < 2e-15);
            for j in 0..surfaces.len() {
                let forward = surfaces[i].area_m2() * open.surface_factors()[i][j];
                let reverse = surfaces[j].area_m2() * open.surface_factors()[j][i];
                assert!((forward - reverse).abs() < 1e-12);
            }
        }
        assert!(
            matches!(open.evidence(), ViewFactorEvidence::GeometryQmc { estimate: hash, .. }
            if *hash == estimate.identity())
        );
        let enclosure = OpenDiffuseEnclosure::new(
            cx,
            surfaces.clone(),
            open.surface_factors().to_vec(),
            open.evidence().clone(),
            open.tolerance(),
        )
        .unwrap();
        let report = enclosure
            .solve(cx, &vec![400.0; surfaces.len()], 300.0)
            .unwrap();
        // Uniform gray walls surrounding a cooler opening must all cool.
        assert!(report.surface_heat_w().iter().all(|&power| power > 0.0));
        assert!(report.ambient_heat_w() < 0.0);
        assert!(
            (report.surface_heat_w().iter().sum::<f64>() + report.ambient_heat_w()).abs() < 1e-9
        );
        let repeated = estimate.project_open(cx, projection()).unwrap();
        assert_eq!(repeated, open);
        assert!(matches!(
            with_cancelled_cx(|cancelled| estimate.project_open(cancelled, projection())),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn g1_single_exposed_patch_samples_the_complete_surrounding_reservoir() {
    let (complex, points) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, points).unwrap();
    let surface = RadiationSurface::new(
        &mesh,
        "exposed-face",
        |face| face.centroid[0] == 0.0,
        emissivity(),
    )
    .unwrap();
    with_cx(|cx| {
        let estimate = estimate_diffuse_view_factors(cx, &mesh, &[surface], config()).unwrap();
        let open = estimate.project_open(cx, projection()).unwrap();
        assert_eq!(open.surface_factors(), &[vec![0.0]]);
        assert_eq!(open.ambient_factors(), &[1.0]);
    });
}

#[test]
fn g4_cancellation_and_geometry_work_budgets_refuse() {
    let (mesh, surfaces) = cavity(false);
    let cancelled =
        with_cancelled_cx(|cx| estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()));
    assert!(matches!(cancelled, Err(ConductionError::Cancelled { .. })));
    with_cx(|cx| {
        for limits in [
            DiffuseViewFactorConfig {
                max_triangles: mesh.boundary().len() - 1,
                ..config()
            },
            DiffuseViewFactorConfig {
                max_ray_triangle_tests: 1,
                ..config()
            },
            DiffuseViewFactorConfig {
                samples_per_replicate: 3,
                ..config()
            },
            DiffuseViewFactorConfig {
                replicates: 1,
                ..config()
            },
        ] {
            assert!(estimate_diffuse_view_factors(cx, &mesh, &surfaces, limits).is_err());
        }
        assert!(estimate_diffuse_view_factors(cx, &mesh, &[], config()).is_err());
        assert!(
            estimate_diffuse_view_factors(
                cx,
                &mesh,
                &[surfaces[0].clone(), surfaces[0].clone()],
                config()
            )
            .is_err()
        );
        let estimate = estimate_diffuse_view_factors(cx, &mesh, &surfaces, config()).unwrap();
        let cancelled =
            with_cancelled_cx(|cancelled_cx| estimate.project_closed(cancelled_cx, projection()));
        assert!(matches!(cancelled, Err(ConductionError::Cancelled { .. })));
    });
}

#[test]
fn two_patch_full_visibility_projects_within_the_unit_interval() {
    let (mesh, walls) = cavity(false);
    let remaining_faces: std::collections::BTreeSet<_> = walls
        .iter()
        .skip(1)
        .flat_map(|wall| {
            wall.face_slots()
                .iter()
                .map(|&slot| mesh.boundary()[slot].face)
        })
        .collect();
    let remaining = RadiationSurface::new(
        &mesh,
        "remaining-walls",
        |face| remaining_faces.contains(&face.face),
        emissivity(),
    )
    .unwrap();
    let surfaces = vec![walls[0].clone(), remaining];
    with_cx(|cx| {
        let mut overestimated_reverse_seen = false;
        for seed in 0..16 {
            let estimate = estimate_diffuse_view_factors(
                cx,
                &mesh,
                &surfaces,
                DiffuseViewFactorConfig {
                    seed: config().seed + seed,
                    samples_per_replicate: 128,
                    replicates: 2,
                    ..config()
                },
            )
            .unwrap();
            assert_eq!(estimate.areas_m2(), &[1.0, 9.0]);
            assert_eq!(estimate.factors()[0], [0.0, 1.0]);
            overestimated_reverse_seen |= estimate.factors()[1][0] > 1.0 / 9.0;
            assert!(
                estimate
                    .project_closed(
                        cx,
                        ViewFactorProjection {
                            max_abs_adjustment: 0.0,
                            ..projection()
                        }
                    )
                    .is_err(),
                "the non-dyadic 1/9 reciprocal factor requires an actual correction"
            );
            let projected = estimate
                .project_closed(
                    cx,
                    ViewFactorProjection {
                        max_abs_adjustment: 0.1,
                        ..projection()
                    },
                )
                .expect("full visibility must not round above the admissible factor range");
            let matrix = projected.matrix();
            assert!((matrix.factors()[0][1] - 1.0).abs() <= 1.0e-12);
            assert!(
                matrix
                    .factors()
                    .iter()
                    .flatten()
                    .all(|f| (0.0..=1.0).contains(f))
            );
            assert!(
                matrix
                    .row_sums()
                    .iter()
                    .all(|sum| (sum - 1.0).abs() <= 1.0e-12)
            );
            assert!(matrix.max_reciprocity_residual() <= 1.0e-12);
        }
        assert!(
            overestimated_reverse_seen,
            "this fixed-seed fixture must exercise the pre-normalization overshoot"
        );
    });
}
