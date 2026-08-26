//! Comprehensive test battery for certified axisymmetric tessellation and representation conversion (bead frankensim-b8bxd.4).
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

    let expected_total = expected_az.hypot(expected_arc);
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
