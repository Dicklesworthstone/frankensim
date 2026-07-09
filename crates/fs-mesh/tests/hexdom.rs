//! Hex-dominant conformance (the wqd.18 bead; runs under
//! `frontier-hexmesh`). Acceptance: MBO smoothing decreases the SH9
//! Dirichlet energy monotonically with boundary alignment held (G0);
//! singularity structures valid (a smooth field is singularity-free, a
//! seeded twist is detected, deterministically); hex extraction meets
//! quality targets on the box and the polycube fallback engages with
//! documented decisions; failure routes to IGA/CutFEM with the honest
//! diagnostic; the accuracy-per-DOF harness reports honestly.

#![cfg(feature = "frontier-hexmesh")]

use fs_mesh::hexdom::{
    FrameField, Quat, Route, accuracy_per_dof, cube_group, dirichlet_energy, extract_hex_dominant,
    mbo_step, singular_edges,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-mesh/hexdom\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn axis_quat(axis: usize, angle: f64) -> Quat {
    let (s, c) = ((angle / 2.0).sin(), (angle / 2.0).cos());
    let mut q = Quat {
        w: c,
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    match axis {
        0 => q.x = s,
        1 => q.y = s,
        _ => q.z = s,
    }
    q
}

/// A field with interior twist noise (boundary identity).
fn noisy_field(n: usize) -> FrameField {
    let mut frames = vec![Quat::identity(); n * n * n];
    for k in 1..n - 1 {
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let a = (k * n + j) * n + i;
                let angle = 0.35 * (((i * 7 + j * 13 + k * 29) % 11) as f64 / 11.0 - 0.5);
                frames[a] = axis_quat((i + j + k) % 3, angle);
            }
        }
    }
    FrameField {
        dims: [n, n, n],
        frames,
    }
}

#[test]
fn hd_001_cube_group_and_mbo_monotone() {
    // The matching set is exactly the 24 cube rotations.
    let group = cube_group();
    assert_eq!(group.len(), 24, "the cube group has 24 elements");
    // MBO smoothing decreases the SH9 Dirichlet energy monotonically
    // on the noisy field, and boundary frames stay pinned (aligned).
    let mut field = noisy_field(5);
    let mut energies = vec![dirichlet_energy(&field)];
    for step in 0..3 {
        mbo_step(&mut field, 0x5eed + step);
        energies.push(dirichlet_energy(&field));
    }
    println!(
        "{{\"metric\":\"mbo\",\"energies\":{:?}}}",
        energies
            .iter()
            .map(|e| (e * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    for w in energies.windows(2) {
        assert!(
            w[1] <= w[0] + 1e-9,
            "the SH9 energy decreases per MBO step: {energies:?}"
        );
    }
    assert!(
        energies.last().expect("last") < &(0.5 * energies[0]),
        "smoothing genuinely smooths"
    );
    // Boundary pinned (normal-aligned on the box by construction).
    let n = 5;
    for j in 0..n {
        for i in 0..n {
            let a = (0 * n + j) * n + i;
            assert_eq!(field.frames[a], Quat::identity(), "boundary stays aligned");
        }
    }
    verdict(
        "hd-001",
        "24-element matching group; MBO halves the SH9 Dirichlet energy monotonically \
         in 3 steps with boundary alignment pinned (G0)",
    );
}

#[test]
fn hd_002_singularity_detection_and_determinism() {
    // A smooth (identity) field is singularity-free.
    let smooth = FrameField {
        dims: [4, 4, 4],
        frames: vec![Quat::identity(); 64],
    };
    assert!(
        singular_edges(&smooth).is_empty(),
        "the aligned field has an empty singular set"
    );
    // A QUARTER-TURN WINDING around one plaquette column: frames at
    // the four surrounding cells rotate 0/30/60/90 degrees about z, so
    // each pairwise matching snaps to identity except the closing gap
    // (90 -> 0 snaps to the 90-degree cube rotation) — the holonomy is
    // non-identity and the column is singular. (A first draft used
    // isolated 45-degree cells: every matching snapped consistently
    // and NOTHING was singular — winding, not local twist, is what a
    // singularity is.)
    let n = 4;
    let mut frames = vec![Quat::identity(); n * n * n];
    for k in 0..n {
        let cells = [(1usize, 1usize), (2, 1), (2, 2), (1, 2)];
        for (e, &(ci, cj)) in cells.iter().enumerate() {
            let a = std::f64::consts::FRAC_PI_6 * e as f64; // 0, 30, 60, 90 deg
            frames[(k * n + cj) * n + ci] = axis_quat(2, a);
        }
    }
    let twisted = FrameField {
        dims: [n, n, n],
        frames,
    };
    let singular = singular_edges(&twisted);
    println!(
        "{{\"metric\":\"singularities\",\"count\":{},\"sample\":{:?}}}",
        singular.len(),
        singular.first()
    );
    assert!(
        !singular.is_empty(),
        "the quarter-turn winding column IS detected as singular"
    );
    // Determinism: identical runs give identical singular sets.
    assert_eq!(singular, singular_edges(&twisted), "bitwise deterministic");
    verdict(
        "hd-002",
        "aligned fields are singularity-free; the seeded twist column is detected; the \
         singular set replays bitwise",
    );
}

#[test]
fn hd_003_box_extraction_full_hex_quality() {
    let field = FrameField {
        dims: [8, 8, 8],
        frames: vec![Quat::identity(); 512],
    };
    let mesh = extract_hex_dominant(&field, &|_, _, _| true, 0.8);
    println!(
        "{{\"metric\":\"box-hex\",\"hexes\":{},\"transitions\":{},\"fraction\":{:.3},\
         \"min_sj\":{:.3},\"route\":\"{:?}\"}}",
        mesh.hexes, mesh.transitions, mesh.hex_fraction, mesh.min_scaled_jacobian, mesh.route
    );
    assert_eq!(
        mesh.route,
        Route::FrameField,
        "integrable: the direct route"
    );
    assert!(
        (mesh.hex_fraction - 1.0).abs() < 1e-12,
        "the box is pure hex"
    );
    assert!(
        (mesh.min_scaled_jacobian - 1.0).abs() < 1e-12,
        "axis-aligned lattice hexes are perfect"
    );
    verdict(
        "hd-003",
        "the box extracts 100% hexes at scaled Jacobian exactly 1.0 through the \
         frame-field route",
    );
}

#[test]
fn hd_004_polycube_fallback_and_honest_refusal() {
    // A stepped (staircase) domain through the twisted field: the
    // singular set forces the POLYCUBE fallback, whose decision is
    // documented in the diagnostic.
    let n = 6;
    let mut frames = vec![Quat::identity(); n * n * n];
    for k in 0..n {
        frames[(k * n + 2) * n + 2] = axis_quat(2, std::f64::consts::FRAC_PI_4);
        frames[(k * n + 3) * n + 3] = axis_quat(0, std::f64::consts::FRAC_PI_4);
    }
    let field = FrameField {
        dims: [n, n, n],
        frames,
    };
    let stepped = |x: f64, y: f64, _z: f64| x + y < 1.4;
    let mesh = extract_hex_dominant(&field, &stepped, 0.5);
    println!(
        "{{\"metric\":\"fallback\",\"route\":\"{:?}\",\"fraction\":{:.3},\
         \"diagnostic\":\"{}\"}}",
        mesh.route, mesh.hex_fraction, mesh.diagnostic
    );
    assert_eq!(mesh.route, Route::Polycube, "the fallback engaged");
    assert!(
        mesh.diagnostic.contains("singular") || mesh.diagnostic.contains("energy"),
        "the decision is documented: {}",
        mesh.diagnostic
    );
    assert!(mesh.hex_fraction >= 0.5, "hex-dominance target met");
    // A hostile thin-shell domain misses any reasonable hex floor:
    // the HONEST refusal routes to IGA/CutFEM.
    let shell = |x: f64, _y: f64, _z: f64| (x - 0.5).abs() < 0.08;
    let refused = extract_hex_dominant(&field, &shell, 0.9);
    assert_eq!(refused.route, Route::RoutedToAlternatives);
    assert!(
        refused.diagnostic.contains("IGA") && refused.diagnostic.contains("CutFEM"),
        "the diagnostic routes to the honest alternatives: {}",
        refused.diagnostic
    );
    verdict(
        "hd-004",
        "the twist-seeded field engages the polycube fallback with its decision \
         documented; the hostile thin shell refuses and routes to IGA/CutFEM by name",
    );
}

#[test]
fn hd_005_accuracy_per_dof_reported_honestly() {
    // The comparison harness: trilinear hex vs linear tet at matched
    // nodes on a smooth field. The REPORT is the deliverable — if hex
    // loses, the numbers say so (here trilinear should win on the
    // smooth fixture, but the assert demands only the honest record).
    let (hex_err, tet_err) = accuracy_per_dof(12);
    let winner = if hex_err < tet_err { "hex" } else { "tet" };
    println!(
        "{{\"metric\":\"accuracy-per-dof\",\"hex_l2\":{hex_err:.3e},\"tet_l2\":{tet_err:.3e},\
         \"winner\":\"{winner}\"}}"
    );
    assert!(hex_err.is_finite() && tet_err.is_finite() && hex_err > 0.0 && tet_err > 0.0);
    verdict(
        "hd-005",
        "accuracy-per-DOF measured and reported for both element classes — the honest \
         comparison the doctrine demands, whichever way it falls",
    );
}
