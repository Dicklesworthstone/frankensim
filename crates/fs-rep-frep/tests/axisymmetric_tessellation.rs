//! Axisymmetric tessellation observations and representation conversion (bead frankensim-b8bxd.4).
//!
//! Verifies:
//! - Analytic sagitta bounds over azimuthal and meridional features;
//! - Strict compliance with error budgets (Hausdorff upper bound <= requested budget);
//! - Monotone mesh refinement under tightened error budgets;
//! - Topological invariants: watertightness, outward orientation, Euler characteristic chi = 2;
//! - Feature tracking: every triangle preserves its generating meridian feature ID;
//! - Strongly-typed domain artifacts: distinct `AxisymmetricRenderMesh` vs `AxisymmetricCollisionMesh`;
//! - Structured refusals: non-positive budget, infeasible budget exceeding caps, empty/degenerate chart;
//! - Cooperative cancellation at checkpoints;
//! - Bit-identical deterministic replay.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::SagittaEnclosure;
use fs_rep_frep::{
    AxisymmetricChart, AxisymmetricCollisionMesh, AxisymmetricRenderMesh,
    AxisymmetricTessellationConfig, AxisymmetricTessellationError, SquatDiscEdgeTreatment,
    TessellationPurpose, tessellate_axisymmetric,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xB8B_4001,
                kernel_id: 42,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

#[test]
fn axt_001_sharp_disc_tessellation_satisfies_budget() {
    with_cx(|cx| {
        let outer_radius = 1.0;
        let thickness = 0.4;
        let disc =
            AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
                .expect("sharp disc");

        let budget = 0.05;
        let config = AxisymmetricTessellationConfig::new(budget, TessellationPurpose::Rendering)
            .expect("valid config");

        let mesh = tessellate_axisymmetric(&disc, config, cx).expect("tessellation");

        assert!(mesh.receipt.total_hausdorff_bound <= budget);
        assert!(mesh.receipt.is_watertight);
        assert!(mesh.receipt.is_outward_oriented);
        assert_eq!(mesh.receipt.euler_characteristic, 2);
        assert_eq!(mesh.receipt.purpose, TessellationPurpose::Rendering);

        // Positions and normals are finite
        for p in &mesh.positions {
            assert!(p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
        }
        for n in &mesh.normals {
            assert!(n.x.is_finite() && n.y.is_finite() && n.z.is_finite());
            let len = (n.x * n.x + n.y * n.y + n.z * n.z).sqrt();
            assert!((len - 1.0).abs() < 1e-6);
        }

        // Triangle count and feature mapping
        assert_eq!(mesh.triangles.len(), mesh.triangle_features.len());
        assert!(!mesh.triangles.is_empty());
    });
}

#[test]
fn axt_002_filleted_disc_refines_arc_features_by_budget() {
    with_cx(|cx| {
        let outer_radius = 1.0;
        let thickness = 0.4;
        let fillet = 0.1;
        let disc = AxisymmetricChart::squat_disc(
            outer_radius,
            thickness,
            SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
        )
        .expect("filleted disc");

        let coarse_budget = 0.05;
        let fine_budget = 0.005;

        let coarse_config =
            AxisymmetricTessellationConfig::new(coarse_budget, TessellationPurpose::Collision)
                .expect("coarse config");
        let fine_config =
            AxisymmetricTessellationConfig::new(fine_budget, TessellationPurpose::Collision)
                .expect("fine config");

        let coarse_mesh =
            tessellate_axisymmetric(&disc, coarse_config, cx).expect("coarse tessellation");
        let fine_mesh = tessellate_axisymmetric(&disc, fine_config, cx).expect("fine tessellation");

        // Fine mesh must have strictly more triangles and smaller bounds
        assert!(fine_mesh.triangles.len() > coarse_mesh.triangles.len());
        assert!(fine_mesh.positions.len() > coarse_mesh.positions.len());
        assert!(fine_mesh.receipt.total_hausdorff_bound <= fine_budget);
        assert!(coarse_mesh.receipt.total_hausdorff_bound <= coarse_budget);
        assert!(
            fine_mesh.receipt.total_hausdorff_bound < coarse_mesh.receipt.total_hausdorff_bound
        );

        // Topology is preserved
        assert!(coarse_mesh.receipt.is_watertight);
        assert!(fine_mesh.receipt.is_watertight);
        assert_eq!(coarse_mesh.receipt.euler_characteristic, 2);
        assert_eq!(fine_mesh.receipt.euler_characteristic, 2);
    });
}

#[test]
fn axt_003_distinct_typed_rendering_and_collision_wrappers() {
    with_cx(|cx| {
        let outer_radius = 0.8;
        let thickness = 0.3;
        let disc =
            AxisymmetricChart::squat_disc(outer_radius, thickness, SquatDiscEdgeTreatment::Sharp)
                .expect("sharp disc");

        let config_render =
            AxisymmetricTessellationConfig::new(0.02, TessellationPurpose::Rendering)
                .expect("render config");
        let config_collision =
            AxisymmetricTessellationConfig::new(0.02, TessellationPurpose::Collision)
                .expect("collision config");

        let mesh_render = tessellate_axisymmetric(&disc, config_render, cx).expect("render mesh");
        let mesh_collision =
            tessellate_axisymmetric(&disc, config_collision, cx).expect("collision mesh");

        let render_artifact = AxisymmetricRenderMesh(mesh_render);
        let collision_artifact = AxisymmetricCollisionMesh(mesh_collision);

        assert_eq!(
            render_artifact.0.receipt.purpose,
            TessellationPurpose::Rendering
        );
        assert_eq!(
            collision_artifact.0.receipt.purpose,
            TessellationPurpose::Collision
        );
    });
}

#[test]
fn axt_004_refusals_on_invalid_and_infeasible_budgets() {
    // Non-positive budget
    let result_neg = AxisymmetricTessellationConfig::new(-0.01, TessellationPurpose::Rendering);
    assert!(matches!(
        result_neg,
        Err(AxisymmetricTessellationError::InvalidBudget { .. })
    ));

    let result_zero = AxisymmetricTessellationConfig::new(0.0, TessellationPurpose::Rendering);
    assert!(matches!(
        result_zero,
        Err(AxisymmetricTessellationError::InvalidBudget { .. })
    ));

    let result_nan = AxisymmetricTessellationConfig::new(f64::NAN, TessellationPurpose::Rendering);
    assert!(matches!(
        result_nan,
        Err(AxisymmetricTessellationError::InvalidBudget { .. })
    ));
    with_cx(|cx| {
        let disc = AxisymmetricChart::squat_disc(
            2.0,
            2.0,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
        )
        .unwrap();
        let mut config =
            AxisymmetricTessellationConfig::new(1e-4, TessellationPurpose::Rendering).unwrap();
        assert!(tessellate_axisymmetric(&disc, config, cx).is_ok());
        config.min_arc_subdivisions = 1;
        config.max_arc_subdivisions = 1;
        let error = tessellate_axisymmetric(&disc, config, cx).unwrap_err();
        let AxisymmetricTessellationError::BudgetInfeasible {
            requested,
            min_achievable,
        } = error
        else {
            panic!("wrong capped-arc refusal: {error:?}");
        };
        let arc_error = 0.5 * (1.0 - core::f64::consts::FRAC_PI_4.cos());
        assert_eq!(requested.to_bits(), 1e-4_f64.to_bits());
        assert!(min_achievable >= arc_error && min_achievable > requested);
    });
}

#[test]
fn axt_005_sagitta_enclosure_analytic_checks() {
    let radius = 2.0;
    let sectors = 16;
    let arc_r = 0.5;
    let arc_sweep = core::f64::consts::FRAC_PI_2;
    let arc_subdivs = 4;

    let sagitta = SagittaEnclosure::compute(radius, sectors, arc_r, arc_sweep, arc_subdivs);

    // Exact azimuthal sagitta: 2.0 * (1 - cos(pi / 16))
    let expected_az = radius * (1.0 - (core::f64::consts::PI / 16.0).cos());
    assert!((sagitta.azimuthal_sagitta - expected_az).abs() < 1e-12);

    // Exact arc sagitta: 0.5 * (1 - cos(pi / 16))
    let expected_arc = arc_r * (1.0 - (arc_sweep / (2.0 * 4.0)).cos());
    assert!((sagitta.meridian_sagitta - expected_arc).abs() < 1e-12);

    let expected_total = expected_az + expected_arc;
    assert!((sagitta.total_hausdorff_bound - expected_total).abs() < 1e-12);
}

#[test]
fn axt_006_deterministic_bit_identical_replay() {
    let disc = AxisymmetricChart::squat_disc(
        1.0,
        0.4,
        SquatDiscEdgeTreatment::CircularFillet { radius: 0.08 },
    )
    .expect("filleted disc");

    let config = AxisymmetricTessellationConfig::new(0.01, TessellationPurpose::Rendering).unwrap();

    let mesh1 = with_cx(|cx| tessellate_axisymmetric(&disc, config, cx).unwrap());
    let mesh2 = with_cx(|cx| tessellate_axisymmetric(&disc, config, cx).unwrap());

    assert_eq!(mesh1.positions.len(), mesh2.positions.len());
    assert_eq!(mesh1.triangles.len(), mesh2.triangles.len());
    assert_eq!(mesh1.triangle_features, mesh2.triangle_features);

    for (p1, p2) in mesh1.positions.iter().zip(&mesh2.positions) {
        assert_eq!(p1.x.to_bits(), p2.x.to_bits());
        assert_eq!(p1.y.to_bits(), p2.y.to_bits());
        assert_eq!(p1.z.to_bits(), p2.z.to_bits());
    }
}

#[test]
fn axt_007_fixed_resolution_budget_boundary_and_adjacent_values() {
    with_cx(|cx| {
        let disc = AxisymmetricChart::squat_disc(
            2.0,
            2.0,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
        )
        .unwrap();
        let theta_error = 2.0 * (1.0 - (core::f64::consts::PI / 100.0).cos());
        let arc_error = 0.5 * (1.0 - core::f64::consts::FRAC_PI_4.cos());
        let bound = theta_error + arc_error;
        let mut config = AxisymmetricTessellationConfig {
            max_hausdorff_error: bound,
            min_azimuthal_sectors: 100,
            max_azimuthal_sectors: 100,
            min_arc_subdivisions: 1,
            max_arc_subdivisions: 1,
            purpose: TessellationPurpose::Rendering,
        };
        for budget in [bound, bound.next_up()] {
            config.max_hausdorff_error = budget;
            let mesh = tessellate_axisymmetric(&disc, config, cx).unwrap();
            assert_eq!(
                mesh.receipt.total_hausdorff_bound.to_bits(),
                bound.to_bits()
            );
            assert!(mesh.receipt.total_hausdorff_bound <= budget);
        }
        config.max_hausdorff_error = bound.next_down();
        assert!(matches!(tessellate_axisymmetric(&disc, config, cx),
            Err(AxisymmetricTessellationError::BudgetInfeasible { min_achievable, .. })
                if min_achievable.to_bits() == bound.to_bits()));

        // A sharp profile has no meridional error: equality to the evaluated
        // azimuthal bound must not be rejected as a zero arc allowance.
        let sharp = AxisymmetricChart::squat_disc(2.0, 2.0, SquatDiscEdgeTreatment::Sharp).unwrap();
        for budget in [theta_error, theta_error.next_up()] {
            config.max_hausdorff_error = budget;
            let mesh = tessellate_axisymmetric(&sharp, config, cx).unwrap();
            assert_eq!(
                mesh.receipt.total_hausdorff_bound.to_bits(),
                theta_error.to_bits()
            );
        }
        config.max_hausdorff_error = theta_error.next_down();
        assert!(matches!(tessellate_axisymmetric(&sharp, config, cx),
            Err(AxisymmetricTessellationError::BudgetInfeasible { min_achievable, .. })
                if min_achievable.to_bits() == theta_error.to_bits()));
    });
}

#[test]
fn axt_008_malformed_limits_and_cancelled_first_error_order() {
    let disc = AxisymmetricChart::squat_disc(1.0, 0.4, SquatDiscEdgeTreatment::Sharp).unwrap();
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 8,
                kernel_id: 42,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let base =
            AxisymmetricTessellationConfig::new(0.05, TessellationPurpose::Rendering).unwrap();
        assert!(matches!(
            tessellate_axisymmetric(&disc, base, &cx),
            Err(AxisymmetricTessellationError::Cancelled)
        ));
        for config in [
            AxisymmetricTessellationConfig {
                min_azimuthal_sectors: 0,
                ..base
            },
            AxisymmetricTessellationConfig {
                min_azimuthal_sectors: 2,
                ..base
            },
            AxisymmetricTessellationConfig {
                max_azimuthal_sectors: 7,
                ..base
            },
            AxisymmetricTessellationConfig {
                min_arc_subdivisions: 0,
                ..base
            },
            AxisymmetricTessellationConfig {
                max_arc_subdivisions: 1,
                ..base
            },
        ] {
            assert!(matches!(
                tessellate_axisymmetric(&disc, config, &cx),
                Err(AxisymmetricTessellationError::InvalidResolutionLimits { .. })
            ));
        }
        for budget in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let config = AxisymmetricTessellationConfig {
                max_hausdorff_error: budget,
                min_arc_subdivisions: 0,
                ..base
            };
            assert!(matches!(
                tessellate_axisymmetric(&disc, config, &cx),
                Err(AxisymmetricTessellationError::InvalidBudget { .. })
            ));
        }
    });
}

#[test]
fn axt_009_provenance_binds_each_semantic_input_and_is_stable() {
    with_cx(|cx| {
        let disc = AxisymmetricChart::squat_disc(1.0, 0.4, SquatDiscEdgeTreatment::Sharp).unwrap();
        let base =
            AxisymmetricTessellationConfig::new(0.05, TessellationPurpose::Rendering).unwrap();
        let original = tessellate_axisymmetric(&disc, base, cx).unwrap();
        let replay = tessellate_axisymmetric(&disc, base, cx).unwrap();
        // Independent 232-byte documented-format fixture (provenance_reference.py).
        assert_eq!(original.receipt.provenance.0, 12_347_313_292_366_907_400);
        assert_eq!(original.receipt.provenance, replay.receipt.provenance);
        for config in [
            AxisymmetricTessellationConfig {
                max_hausdorff_error: 0.06,
                ..base
            },
            AxisymmetricTessellationConfig {
                min_azimuthal_sectors: 9,
                ..base
            },
            AxisymmetricTessellationConfig {
                max_azimuthal_sectors: 4095,
                ..base
            },
            AxisymmetricTessellationConfig {
                min_arc_subdivisions: 3,
                ..base
            },
            AxisymmetricTessellationConfig {
                max_arc_subdivisions: 1023,
                ..base
            },
            AxisymmetricTessellationConfig {
                purpose: TessellationPurpose::Collision,
                ..base
            },
        ] {
            let mesh = tessellate_axisymmetric(&disc, config, cx).unwrap();
            assert_ne!(
                mesh.receipt.provenance, original.receipt.provenance,
                "{config:?}"
            );
            if config.max_hausdorff_error.to_bits() == base.max_hausdorff_error.to_bits()
                && config.min_azimuthal_sectors == base.min_azimuthal_sectors
            {
                assert_eq!(mesh.positions, original.positions);
                assert_eq!(mesh.normals, original.normals);
                assert_eq!(mesh.triangles, original.triangles);
                assert_eq!(mesh.triangle_features, original.triangle_features);
                assert_eq!(
                    mesh.receipt.total_hausdorff_bound.to_bits(),
                    original.receipt.total_hausdorff_bound.to_bits()
                );
            }
        }
        for changed in [
            AxisymmetricChart::squat_disc(1.1, 0.4, SquatDiscEdgeTreatment::Sharp).unwrap(),
            AxisymmetricChart::squat_disc(1.0, 0.5, SquatDiscEdgeTreatment::Sharp).unwrap(),
            AxisymmetricChart::squat_disc(
                1.0,
                0.4,
                SquatDiscEdgeTreatment::CircularFillet { radius: 0.08 },
            )
            .unwrap(),
        ] {
            assert_ne!(
                tessellate_axisymmetric(&changed, base, cx)
                    .unwrap()
                    .receipt
                    .provenance,
                original.receipt.provenance
            );
        }
    });
}

#[test]
fn axt_010_cylinder_faces_point_outward_and_name_their_generating_feature() {
    with_cx(|cx| {
        let disc = AxisymmetricChart::squat_disc(2.0, 2.0, SquatDiscEdgeTreatment::Sharp).unwrap();
        let config =
            AxisymmetricTessellationConfig::new(0.2, TessellationPurpose::Rendering).unwrap();
        let mesh = tessellate_axisymmetric(&disc, config, cx).unwrap();
        assert!(mesh.receipt.is_outward_oriented);
        assert_eq!((mesh.positions.len(), mesh.triangles.len()), (22, 40));
        let mut volume = 0.0;
        let mut seen = [false; 3];
        for (triangle, &feature) in mesh.triangles.iter().zip(&mesh.triangle_features) {
            let [a, b, c] = triangle.map(|i| mesh.positions[i as usize]);
            let nx = (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y);
            let ny = (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z);
            let nz = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
            let expected = if [a, b, c]
                .iter()
                .all(|p| p.z.to_bits() == (-1.0_f64).to_bits())
            {
                assert!(nz < 0.0);
                0
            } else if [a, b, c].iter().all(|p| p.z.to_bits() == 1.0_f64.to_bits()) {
                assert!(nz > 0.0);
                2
            } else {
                assert!(nx * (a.x + b.x + c.x) + ny * (a.y + b.y + c.y) > 0.0);
                1
            };
            assert_eq!(feature, expected);
            seen[expected] = true;
            volume += (a.x * (b.y * c.z - b.z * c.y)
                + a.y * (b.z * c.x - b.x * c.z)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0;
        }
        assert!(seen.iter().all(|x| *x));
        // Independent regular-decagon prism volume, opposite the old native loss.
        let expected = 40.0 * (core::f64::consts::TAU / 10.0).sin();
        assert!(volume > 0.0 && (volume - expected).abs() < 1e-12);
    });
}

#[test]
fn axt_011_fillets_and_bore_have_outward_faces_and_correct_feature_ids() {
    with_cx(|cx| {
        // Known six-feature meridians: bottom, lower fillet, outer wall,
        // upper fillet, top, then axis closure or bore. Include three sectors
        // and one arc chord to expose Cartesian-centroid radial contraction.
        for (bore, fillet, coarse) in [(0.0, 0.5, false), (0.0, 0.05, true), (0.5, 0.5, false)] {
            let chart = if bore == 0.0 {
                AxisymmetricChart::squat_disc(
                    2.0,
                    2.0,
                    SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
                )
                .unwrap()
            } else {
                AxisymmetricChart::annular_disc_outer_fillets(2.0, bore, 2.0, fillet).unwrap()
            };
            let mut config = AxisymmetricTessellationConfig::new(
                if coarse { 2.0 } else { 0.05 },
                TessellationPurpose::Rendering,
            )
            .unwrap();
            if coarse {
                config.min_azimuthal_sectors = 3;
                config.max_azimuthal_sectors = 3;
                config.min_arc_subdivisions = 1;
                config.max_arc_subdivisions = 1;
            }
            let mesh = tessellate_axisymmetric(&chart, config, cx).unwrap();
            assert_eq!(mesh.triangles.len(), mesh.triangle_features.len());
            let mut seen = [false; 6];
            for (triangle, &feature) in mesh.triangles.iter().zip(&mesh.triangle_features) {
                let points = triangle.map(|i| mesh.positions[i as usize]);
                let [a, b, c] = points;
                let radii = points.map(|p| p.x.hypot(p.y));
                let nx = (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y);
                let ny = (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z);
                let nz = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
                let x = (a.x + b.x + c.x) / 3.0;
                let y = (a.y + b.y + c.y) / 3.0;
                let z = (a.z + b.z + c.z) / 3.0;
                let radial_alignment = (nx * x + ny * y) / x.hypot(y);
                let expected = if bore > 0.0 && radii.iter().all(|r| (*r - bore).abs() < 1e-12) {
                    assert!(radial_alignment < 0.0);
                    5
                } else if points.iter().all(|p| (p.z + 1.0).abs() < 1e-12) {
                    assert!(nz < 0.0);
                    0
                } else if points.iter().all(|p| (p.z - 1.0).abs() < 1e-12) {
                    assert!(nz > 0.0);
                    4
                } else if points.iter().all(|p| p.z <= -1.0 + fillet + 1e-12) {
                    let radial = radii.iter().sum::<f64>() / 3.0 - (2.0 - fillet);
                    assert!(radial_alignment * radial + nz * (z + 1.0 - fillet) > 0.0);
                    1
                } else if points.iter().all(|p| p.z >= 1.0 - fillet - 1e-12) {
                    let radial = radii.iter().sum::<f64>() / 3.0 - (2.0 - fillet);
                    assert!(radial_alignment * radial + nz * (z - 1.0 + fillet) > 0.0);
                    3
                } else {
                    assert!(radii.iter().all(|r| (*r - 2.0).abs() < 1e-12));
                    assert!(radial_alignment > 0.0);
                    2
                };
                assert_eq!(
                    feature, expected,
                    "bore={bore}, fillet={fillet}, coarse={coarse}"
                );
                seen[expected] = true;
            }
            assert!(seen[..5].iter().all(|s| *s));
            assert_eq!(seen[5], bore > 0.0);
            if bore == 0.0 {
                assert!(mesh.receipt.is_watertight && mesh.receipt.is_outward_oriented);
                assert_eq!(mesh.receipt.euler_characteristic, 2);
            }
            // Bore seam closure is an independent existing topology concern;
            // per-face direction does not assert that its mesh is watertight.
        }
    });
}

#[test]
fn axt_012_composed_bound_covers_a_true_fillet_point_beyond_mesh_support_plane() {
    with_cx(|cx| {
        let chart = AxisymmetricChart::squat_disc(
            2.0,
            2.0,
            SquatDiscEdgeTreatment::CircularFillet { radius: 0.5 },
        )
        .unwrap();
        let config = AxisymmetricTessellationConfig {
            max_hausdorff_error: 0.01,
            min_azimuthal_sectors: 45,
            max_azimuthal_sectors: 45,
            min_arc_subdivisions: 6,
            max_arc_subdivisions: 6,
            purpose: TessellationPurpose::Rendering,
        };
        let mesh = tessellate_axisymmetric(&chart, config, cx).unwrap();
        // This facet is retained in the actual native counterexample. The
        // complete mesh lies below its outward plane; hence signed distance
        // beyond that plane is a lower bound on distance to every mesh face.
        assert!(mesh.triangles.contains(&[226, 227, 272]));
        let [a, b, c] = [226, 227, 272].map(|i| mesh.positions[i]);
        let nx = (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y);
        let ny = (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z);
        let nz = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        let length = nx.hypot(ny).hypot(nz);
        let plane =
            |p: fs_geom::Point3| (nx * (p.x - a.x) + ny * (p.y - a.y) + nz * (p.z - a.z)) / length;
        assert!(mesh.positions.iter().all(|p| plane(*p) <= 1e-12));
        let phi = -core::f64::consts::PI / 24.0;
        let theta = core::f64::consts::PI / 45.0;
        let radius = 1.5 + 0.5 * phi.cos();
        let point = fs_geom::Point3::new(
            radius * theta.cos(),
            radius * theta.sin(),
            -0.5 + 0.5 * phi.sin(),
        );
        let distance_lower_observation = plane(point) - 1e-12;
        let receipt = &mesh.receipt;
        let obsolete_hypot = receipt
            .azimuthal_sagitta_bound
            .hypot(receipt.meridian_sagitta_bound);
        assert!(distance_lower_observation > obsolete_hypot);
        assert!(receipt.total_hausdorff_bound >= distance_lower_observation);
        assert!(receipt.total_hausdorff_bound <= config.max_hausdorff_error);
    });
}
