//! End-to-end battery: a learned neural implicit whose Lipschitz bound, safe
//! rendering, and bounded-non-empty-region topology are all PROVEN.

use fs_evidence::Color;
use fs_neuroshape_e2e::{blob_sdf_net, run_campaign};

#[test]
fn the_neural_shape_topology_is_certified() {
    let net = blob_sdf_net();
    let report = run_campaign(&net, 2.5, 0.3);
    // a finite certified Lipschitz bound underwrites everything.
    assert!(
        report.lipschitz.is_finite() && report.lipschitz > 0.0,
        "L {}",
        report.lipschitz
    );
    // sound sphere tracing: the origin is negative and the safe step is a
    // positive, finite, non-tunneling distance.
    assert!(report.origin_value < 0.0, "origin {}", report.origin_value);
    assert!(report.safe_radius > 0.0 && report.safe_radius.is_finite());
    // the safe radius must UNDER-estimate the distance to the NEAREST surface
    // point — the actual no-tunnel soundness guarantee (a Lipschitz theorem).
    assert!(
        report.safe_radius < report.nearest_surface_radius,
        "safe {} !< nearest surface {}",
        report.safe_radius,
        report.nearest_surface_radius
    );
    assert!(report.nearest_surface_radius <= report.max_crossing_radius);
    // TOPOLOGY, PROVEN: a certified-inside interior enclosed by a CLOSED,
    // fully-certified boundary frame ⇒ a non-empty, bounded region.
    assert!(
        report.certified_inside,
        "inside interval {:?}",
        report.inside_interval
    );
    assert_eq!(report.boundary_segments, 4);
    assert_eq!(report.boundary_certified, report.boundary_segments);
    assert!(report.bounded);
    assert!(matches!(report.topology_color, Color::Verified { .. }));
    // Morse cross-check: one interior minimum.
    assert!(report.single_minimum);
    // the visualization localizes a closed surface, all inside the ring.
    assert!(report.surface_crossings > 0);
    assert!(
        report.max_crossing_radius < 2.5,
        "surface escaped the ring: {}",
        report.max_crossing_radius
    );
    println!(
        "{{\"campaign\":\"neuroshapecert\",\"L\":{:.3},\"origin\":{:.3},\"safe_radius\":{:.3},\
         \"inside\":[{:.3},{:.3}],\"boundary\":{}/{},\"single_min\":{},\"crossings\":{},\
         \"max_r\":{:.3}}}",
        report.lipschitz,
        report.origin_value,
        report.safe_radius,
        report.inside_interval.0,
        report.inside_interval.1,
        report.boundary_certified,
        report.boundary_segments,
        report.single_minimum,
        report.surface_crossings,
        report.max_crossing_radius,
    );
}

#[test]
fn an_open_ring_yields_no_topology_certificate() {
    // too small a box: its boundary frame overlaps the surface → not certified.
    let net = blob_sdf_net();
    let report = run_campaign(&net, 0.3, 0.3);
    assert!(!report.bounded || !report.certified_inside);
    assert!(matches!(report.topology_color, Color::Estimated { .. }));
}

#[test]
fn the_campaign_is_deterministic() {
    let net = blob_sdf_net();
    let a = run_campaign(&net, 2.5, 0.3);
    let b = run_campaign(&net, 2.5, 0.3);
    assert_eq!(a.lipschitz.to_bits(), b.lipschitz.to_bits());
    assert_eq!(a.surface_crossings, b.surface_crossings);
    assert_eq!(a.safe_radius.to_bits(), b.safe_radius.to_bits());
}
