use super::*;

#[test]
fn resolved_geometry_keeps_source_features_and_material_authority() {
    for published in [false, true] {
        for (divisions, step) in [(8, 0.025), (16, 0.0125), (32, 0.02), (48, 0.025)] {
            let preset = build_resolved(divisions, step, published).unwrap();
            super::super::board_geometry::BoardGeometry::read(&preset.geometry).unwrap();
            assert_eq!(preset.geometry.lines().filter(|s| s.starts_with("bridge,")).count(), 88);
            assert!(preset.geometry.contains("source,mixed,"));
            assert!(preset.geometry.contains("support,clamped"));
            assert!(preset.geometry.contains("pretension,0"));
            assert!(preset.geometry.contains("damping,0.015"));
            for rib in 1..=17 {
                assert!(preset.obj.contains(&format!("g sugar_pine_rib_{rib}\n")));
            }
            assert_eq!(preset.geometry, build_resolved(divisions, step, published).unwrap().geometry);
        }
    }
}

fn nodes(geometry: &str) -> Vec<[f64; 2]> {
    geometry.lines().filter(|l| l.starts_with("node,")).map(|l| {
        let f = l.split(',').skip(2).map(|v| v.parse::<f64>().unwrap()).collect::<Vec<_>>();
        [f[0], f[1]]
    }).collect()
}

/// Drawing sweep station t = y + SLOPE*x of an SI node.
fn sweep_station(p: [f64; 2]) -> f64 {
    let raster = [p[0] / scale() + 55.0, 970.0 - p[1] / scale()];
    raster[1] + SLOPE * raster[0]
}

#[test]
fn unrefined_resolution_differs_from_the_default_only_by_replaced_regular_nodes() {
    let old = build(8).unwrap();
    let resolved = build_resolved(8, 1.0, false).unwrap();
    let mesh = ["#", "node,", "triangle,", "stiffener,", "bridge,", "fixed,"];
    let rows = |g: &str| g.lines().filter(|l| !mesh.iter().any(|m| l.starts_with(m)))
        .map(str::to_owned).collect::<Vec<_>>();
    assert_eq!(rows(&resolved.geometry), rows(&old.geometry));
    let (kept, all) = (nodes(&resolved.geometry), nodes(&old.geometry));
    assert!(kept.len() < all.len(), "the default d8 board has crossings beside regular nodes");
    // Every resolved node is a default node, bit for bit, and the clamped
    // boundary keeps every default clamped node that survives.
    assert!(kept.iter().all(|p| all.contains(p)));
    let fixed = |g: &str, xy: &[[f64; 2]]| g.lines().filter_map(|l| l.strip_prefix("fixed,"))
        .map(|v| xy[v.parse::<usize>().unwrap()]).collect::<Vec<_>>();
    let (fixed_kept, fixed_all) = (fixed(&resolved.geometry, &kept), fixed(&old.geometry, &all));
    assert!(fixed_kept.iter().all(|p| fixed_all.contains(p)));
    assert!(fixed_all.iter().all(|p| !kept.contains(p) || fixed_kept.contains(p)));
    // Every dropped node had a structural crossing within a quarter of its
    // row's regular spacing; that crossing is a retained node of the same row.
    for p in all.iter().filter(|p| !kept.contains(p)) {
        let t = sweep_station(*p);
        let row = all.iter().filter(|q| (sweep_station(**q) - t).abs() < 1e-6).collect::<Vec<_>>();
        let x = |q: &[f64; 2]| q[0] / scale() + 55.0;
        let (lo, hi) = row.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), q| (lo.min(x(q)), hi.max(x(q))));
        let near = row.iter().filter(|q| kept.contains(q))
            .map(|q| (x(q) - x(p)).abs()).fold(f64::INFINITY, f64::min);
        assert!(near < 0.25 * (hi - lo) / 8.0, "node at station {t} dropped without a nearby crossing");
    }
    assert_eq!(resolved.geometry.lines().filter(|l| l.starts_with("bridge,")).count(), 88);
    assert!(build(33).is_err());
    assert!(build_with_rt0425_contacts(33).is_err());
}

/// frankensim-52edk: these exported resolutions once failed the shift-0
/// factorization because a bridge crossing beside a regular rib-row node left
/// rib segments 0.04-0.13 mm long (12EI/L^3 dwarfed the panel's pivots).
#[test]
fn resolved_rib_rows_keep_no_sliver_beam_segments() {
    for (divisions, step) in [(20, 0.016), (20, 0.015), (22, 0.0155), (12, 0.03)] {
        let preset = build_resolved(divisions, step, false).unwrap();
        let xy = nodes(&preset.geometry);
        let shortest = preset.geometry.lines().filter(|l| l.starts_with("stiffener,")).map(|l| {
            let ends = l.split(',').rev().take(2).map(|v| v.parse::<usize>().unwrap()).collect::<Vec<_>>();
            let (a, b) = (xy[ends[0]], xy[ends[1]]);
            (a[0] - b[0]).hypot(a[1] - b[1])
        }).fold(f64::INFINITY, f64::min);
        // 1.479 mm is the drawing's own closest pair of feature stations
        // (a long-bridge vertex 0.8 raster units from rib 3), present at every resolution.
        assert!(shortest > 1.4e-3, "{divisions}/{step}: {shortest} m stiffener segment");
    }
}

#[test]
fn formerly_refused_resolutions_factor_at_shift_zero() {
    for (divisions, step) in [(20, 0.016), (20, 0.015), (22, 0.0155)] {
        let preset = build_resolved(divisions, step, false).unwrap();
        let board = super::super::board_geometry::BoardGeometry::read(&preset.geometry).unwrap();
        board.prepare(&[57], 150.0).unwrap_or_else(|e| panic!("{divisions}/{step}: {e}"));
    }
}

#[test]
fn resolution_requires_explicit_finite_physical_inputs() {
    for count in [0, 3, 65, usize::MAX] {
        assert!(build_resolved(count, 0.02, true).is_err());
    }
    for step in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MIN_POSITIVE] {
        assert!(build_resolved(16, step, true).is_err());
    }
}
