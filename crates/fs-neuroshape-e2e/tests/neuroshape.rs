//! End-to-end battery: a learned neural implicit whose Lipschitz bound, safe
//! rendering, and existence of an enclosed negative component are PROVEN.
//! The global component count deliberately remains inexact.

use fs_evidence::Color;
use fs_neuroshape_e2e::{ComponentCountEvidence, blob_sdf_net, run_campaign};

#[test]
fn an_enclosed_component_lower_bound_is_certified() {
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
    // fully-certified boundary frame implies at least one enclosed component.
    // It does not exclude additional components inside or outside the frame.
    assert!(
        report.certified_inside,
        "inside interval {:?}",
        report.inside_interval
    );
    assert_eq!(report.boundary_segments, 4);
    assert_eq!(report.boundary_certified, report.boundary_segments);
    assert!(report.boundary_frame_certified);
    assert!(matches!(
        &report.component_enclosure_color,
        Color::Verified { .. }
    ));
    assert_eq!(report.component_count_evidence.lower_bound(), 1);
    assert_eq!(report.component_count_evidence.exact_count(), None);
    let ComponentCountEvidence::LowerBound(enclosed) = &report.component_count_evidence else {
        panic!("closed frame must produce a typed lower-bound witness");
    };
    assert_eq!(
        enclosed.central_box_half_width().to_bits(),
        0.3_f64.to_bits()
    );
    let enclosed_interval = enclosed.central_box_interval();
    assert_eq!(
        enclosed_interval.0.to_bits(),
        report.inside_interval.0.to_bits()
    );
    assert_eq!(
        enclosed_interval.1.to_bits(),
        report.inside_interval.1.to_bits()
    );
    assert_eq!(
        enclosed.boundary_frame_outer_half_width().to_bits(),
        2.5_f64.to_bits()
    );
    assert_eq!(
        enclosed.boundary_frame_inner_half_width().to_bits(),
        2.1_f64.to_bits()
    );
    assert!(
        enclosed
            .boundary_strip_intervals()
            .iter()
            .all(|(lo, hi)| lo.is_finite() && hi.is_finite() && *lo > 0.0 && lo <= hi)
    );
    // The origin Hessian is positive definite, but no zero-gradient witness is
    // present. This curvature check must never promote the lower bound to an
    // exact count or even claim a critical point.
    assert!(report.origin_hessian_positive_definite);
    // the visualization localizes a closed surface, all inside the ring.
    assert!(report.surface_crossings > 0);
    assert!(
        report.max_crossing_radius < 2.5,
        "surface escaped the ring: {}",
        report.max_crossing_radius
    );
    println!(
        "{{\"campaign\":\"neuroshapecert\",\"L\":{:.3},\"origin\":{:.3},\"safe_radius\":{:.3},\
         \"inside\":[{:.3},{:.3}],\"boundary\":{}/{},\"component_count_lower_bound\":{},\
         \"exact_component_count\":null,\"origin_hessian_positive_definite\":{},\
         \"crossings\":{},\"max_r\":{:.3}}}",
        report.lipschitz,
        report.origin_value,
        report.safe_radius,
        report.inside_interval.0,
        report.inside_interval.1,
        report.boundary_certified,
        report.boundary_segments,
        report.component_count_evidence.lower_bound(),
        report.origin_hessian_positive_definite,
        report.surface_crossings,
        report.max_crossing_radius,
    );
}

#[test]
fn an_open_ring_yields_no_topology_certificate() {
    // too small a box: its boundary frame overlaps the surface → not certified.
    let net = blob_sdf_net();
    let report = run_campaign(&net, 0.3, 0.3);
    assert!(!report.boundary_frame_certified || !report.certified_inside);
    assert!(matches!(
        &report.component_enclosure_color,
        Color::Estimated { .. }
    ));
    assert!(matches!(
        &report.component_count_evidence,
        ComponentCountEvidence::Unknown
    ));
    assert_eq!(report.component_count_evidence.lower_bound(), 0);
    assert_eq!(report.component_count_evidence.exact_count(), None);
}

#[test]
fn the_campaign_is_deterministic() {
    let net = blob_sdf_net();
    let a = run_campaign(&net, 2.5, 0.3);
    let b = run_campaign(&net, 2.5, 0.3);
    assert_eq!(a.lipschitz.to_bits(), b.lipschitz.to_bits());
    assert_eq!(a.surface_crossings, b.surface_crossings);
    assert_eq!(a.safe_radius.to_bits(), b.safe_radius.to_bits());
    assert_eq!(a.component_count_evidence, b.component_count_evidence);
}
