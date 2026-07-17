//! G3 metamorphic wiring against the PV-skeleton physics kernel
//! (bead frankensim-epic-gauntlet-6nb.4).
//!
//! `EdgeLaw` is the skeleton's one-source-of-truth physics: both the primal
//! operator and the adjoint contraction consume it, so a frame or unit bug
//! here poisons everything downstream. Two canonical relations are wired
//! directly against it, and a seeded violator for each proves the harness
//! actually catches the bug class it exists for.

use fs_propcheck::Stream;
use fs_propcheck::metamorphic::{
    RelationCase, Tolerance, check_relation, rigid_motion, unit_rescaling,
};
use fs_vskeleton::model::EdgeLaw;

const HOLE_CENTER: (f64, f64) = (0.5, 0.5);

fn law(radius: f64, width: f64) -> EdgeLaw {
    EdgeLaw {
        eps: 1e-3,
        width,
        radius,
    }
}

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-vskeleton/metamorphic\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{detail}\"}}"
    );
}

/// Frame invariance: the ersatz density is a function of DISTANCE from the
/// hole center, so rotating the sample point about that center must leave
/// rho (and its radius sensitivity) unchanged to rotation round-off.
#[test]
fn g3_edge_law_density_is_frame_invariant_under_rotation_about_the_hole() {
    // Input: (radius, dx, dy) — a sample offset from the hole center.
    // Transform: rotation angle about the center.
    let relation = rigid_motion(
        "pv-skeleton-rho-rotation",
        Tolerance::AbsoluteRelative {
            max_abs: 1e-12,
            max_relative: 1e-12,
        },
        |input: &(f64, f64, f64), angle: &f64| {
            let (radius, dx, dy) = *input;
            let (sin, cos) = angle.sin_cos();
            (radius, dx * cos - dy * sin, dx * sin + dy * cos)
        },
        |base: &f64, transformed: &f64, _angle: &f64, tolerance: Tolerance| {
            tolerance.evaluate_scalar(*base, *transformed)
        },
    );
    let generate = |stream: &mut Stream| {
        let radius = stream.f64_in(0.05, 0.45);
        let dx = stream.f64_in(-0.5, 0.5);
        let dy = stream.f64_in(-0.5, 0.5);
        let angle = stream.f64_in(0.0, core::f64::consts::TAU);
        RelationCase::new((radius, dx, dy), angle)
    };
    let rho = |input: &(f64, f64, f64)| {
        let (radius, dx, dy) = *input;
        law(radius, 0.03).rho(HOLE_CENTER.0 + dx, HOLE_CENTER.1 + dy)
    };
    let drho = |input: &(f64, f64, f64)| {
        let (radius, dx, dy) = *input;
        law(radius, 0.03).d_rho_d_radius(HOLE_CENTER.0 + dx, HOLE_CENTER.1 + dy)
    };
    check_relation(
        "pv-skeleton-rho",
        0x6E_B4_01,
        512,
        generate,
        &rho,
        &relation,
    );
    check_relation(
        "pv-skeleton-drho-dr",
        0x6E_B4_02,
        512,
        generate,
        &drho,
        &relation,
    );
    verdict(
        "frame-invariance",
        "rho and d(rho)/dr invariant under 1024 random rotations about the hole center",
    );
}

/// Unit rescaling THROUGH RUNTIME PATHS: scale every length (offset, hole
/// radius, Heaviside width) by one factor. rho is dimensionless and must be
/// unchanged; d(rho)/dr carries dimension 1/length and must scale exactly by
/// the inverse factor.
#[test]
fn g3_edge_law_units_rescale_coherently() {
    type UnitInput = ((f64, f64), f64, f64); // ((radius, width), dx, dy)

    let dimensionless = unit_rescaling(
        "pv-skeleton-rho-unit-rescale",
        Tolerance::AbsoluteRelative {
            max_abs: 1e-12,
            max_relative: 1e-12,
        },
        |input: &UnitInput, scale: &f64| {
            let ((radius, width), dx, dy) = *input;
            ((radius * scale, width * scale), dx * scale, dy * scale)
        },
        |base: &f64, transformed: &f64, _scale: &f64, tolerance: Tolerance| {
            tolerance.evaluate_scalar(*base, *transformed)
        },
    );
    let dimensioned = unit_rescaling(
        "pv-skeleton-drho-dr-unit-rescale",
        Tolerance::AbsoluteRelative {
            max_abs: 1e-12,
            max_relative: 1e-12,
        },
        |input: &UnitInput, scale: &f64| {
            let ((radius, width), dx, dy) = *input;
            ((radius * scale, width * scale), dx * scale, dy * scale)
        },
        // d(rho)/dr has dimension 1/length: candidate * scale must equal base.
        |base: &f64, transformed: &f64, scale: &f64, tolerance: Tolerance| {
            tolerance.evaluate_scalar(*base, *transformed * *scale)
        },
    );
    let generate = |stream: &mut Stream| {
        let radius = stream.f64_in(0.05, 0.45);
        let width = stream.f64_in(0.01, 0.1);
        let dx = stream.f64_in(-0.5, 0.5);
        let dy = stream.f64_in(-0.5, 0.5);
        // Millimetres to kilometres relative to the base unit.
        let scale = stream.f64_in(1e-3, 1e3);
        RelationCase::new(((radius, width), dx, dy), scale)
    };
    // The scaled evaluation places the hole center at the scaled origin
    // offset, exercising the full length pipeline rather than spec constants.
    let rho = |input: &UnitInput| {
        let ((radius, width), dx, dy) = *input;
        let law = law(radius, width);
        // signed_gap uses the fixed (0.5, 0.5) center; sample in centered
        // coordinates so every length in play is scaled by the transform.
        law.rho(HOLE_CENTER.0 + dx, HOLE_CENTER.1 + dy)
    };
    let drho = |input: &UnitInput| {
        let ((radius, width), dx, dy) = *input;
        law(radius, width).d_rho_d_radius(HOLE_CENTER.0 + dx, HOLE_CENTER.1 + dy)
    };
    check_relation(
        "pv-skeleton-rho",
        0x6E_B4_03,
        512,
        generate,
        &rho,
        &dimensionless,
    );
    check_relation(
        "pv-skeleton-drho-dr",
        0x6E_B4_04,
        512,
        generate,
        &drho,
        &dimensioned,
    );
    verdict(
        "unit-rescaling",
        "rho dimensionless-invariant and d(rho)/dr scales as 1/length across 1024 rescalings",
    );
}

/// Seeded violators: the harness must CATCH the exact bug classes it exists
/// for — an anisotropic distance (frame bug) and a hardcoded absolute width
/// (unit bug, the classic silent meter assumption).
#[test]
fn g3_seeded_frame_and_unit_violations_are_caught() {
    // Frame violator: y-weighted distance breaks rotational symmetry.
    let frame_relation = rigid_motion(
        "pv-skeleton-rho-rotation",
        Tolerance::AbsoluteRelative {
            max_abs: 1e-12,
            max_relative: 1e-12,
        },
        |input: &(f64, f64, f64), angle: &f64| {
            let (radius, dx, dy) = *input;
            let (sin, cos) = angle.sin_cos();
            (radius, dx * cos - dy * sin, dx * sin + dy * cos)
        },
        |base: &f64, transformed: &f64, _angle: &f64, tolerance: Tolerance| {
            tolerance.evaluate_scalar(*base, *transformed)
        },
    );
    let anisotropic_rho = |input: &(f64, f64, f64)| {
        let (radius, dx, dy) = *input;
        let gap = (dx * dx + 2.0 * dy * dy).sqrt() - radius;
        1e-3 + (1.0 - 1e-3) * 0.5 * (1.0 + (gap / 0.03).tanh())
    };
    let frame_caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        check_relation(
            "pv-skeleton-anisotropic-rho",
            0x6E_B4_05,
            512,
            |stream: &mut Stream| {
                let radius = stream.f64_in(0.05, 0.45);
                let dx = stream.f64_in(-0.5, 0.5);
                let dy = stream.f64_in(-0.5, 0.5);
                let angle = stream.f64_in(0.5, 1.5);
                RelationCase::new((radius, dx, dy), angle)
            },
            &anisotropic_rho,
            &frame_relation,
        );
    }));
    assert!(
        frame_caught.is_err(),
        "the anisotropic-distance frame bug must be caught by the rotation relation"
    );

    // Unit violator: the Heaviside width stays in absolute units while every
    // declared length scales — the hardcoded meter assumption.
    let unit_relation = unit_rescaling(
        "pv-skeleton-rho-unit-rescale",
        Tolerance::AbsoluteRelative {
            max_abs: 1e-12,
            max_relative: 1e-12,
        },
        |input: &((f64, f64), f64, f64), scale: &f64| {
            let ((radius, _width), dx, dy) = *input;
            // The transform scales radius and offsets, but the violator's
            // width is baked into the operator below and never scales.
            ((radius * scale, 0.03), dx * scale, dy * scale)
        },
        |base: &f64, transformed: &f64, _scale: &f64, tolerance: Tolerance| {
            tolerance.evaluate_scalar(*base, *transformed)
        },
    );
    let fixed_width_rho = |input: &((f64, f64), f64, f64)| {
        let ((radius, _width), dx, dy) = *input;
        let gap = (dx * dx + dy * dy).sqrt() - radius;
        1e-3 + (1.0 - 1e-3) * 0.5 * (1.0 + (gap / 0.03).tanh())
    };
    let unit_caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        check_relation(
            "pv-skeleton-fixed-width-rho",
            0x6E_B4_06,
            512,
            |stream: &mut Stream| {
                let radius = stream.f64_in(0.05, 0.45);
                let dx = stream.f64_in(-0.5, 0.5);
                let dy = stream.f64_in(-0.5, 0.5);
                let scale = stream.f64_in(2.0, 100.0);
                RelationCase::new(((radius, 0.03), dx, dy), scale)
            },
            &fixed_width_rho,
            &unit_relation,
        );
    }));
    assert!(
        unit_caught.is_err(),
        "the hardcoded-width unit bug must be caught by the rescaling relation"
    );
    verdict(
        "seeded-violations",
        "anisotropic-distance frame bug and hardcoded-width unit bug both refused",
    );
}
