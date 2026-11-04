#![allow(dead_code)]
#[path = "../source_point/common.rs"]
pub mod common;

use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_contact::normal_patch::dmt::{DmtEnvelope, SourceDmtSphereContact};
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
use fs_contact::normal_patch::{
    FiniteGapChartSamplingAuthority, FiniteGapResponseCurve, FiniteGapResponseCurveAuthority,
    FiniteGapResponseNode,
};

pub fn jkr() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    let normal = common::bind(&a, &b, None, 0.05).normal_state().clone();
    SourceJkrSphereContact::new(
        normal,
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap()
}

pub fn point(chi: Option<f64>) -> SourcePointContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
    let mut envelope = common::envelope();
    envelope.maximum_indentation_m = 1e-5;
    SourcePointContact::new(
        common::bind(&a, &b, chi, 0.0),
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.01,
        },
        envelope,
    )
    .unwrap()
}

pub fn dmt() -> SourceDmtSphereContact {
    SourceDmtSphereContact::new(
        jkr().normal_state().clone(),
        0.01,
        1e-5,
        DmtEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}

/// Manufactured retained curve; no claim that these nodes came from a BEM solve.
pub fn curve(length: f64, force: f64) -> FiniteGapResponseCurve {
    let identity = fs_blake3::hash_domain("test/mixed-normal-derivative", b"manufactured curve");
    let mut nodes: Vec<FiniteGapResponseNode> = Vec::new();
    for (x, y) in [(0.0, 0.0), (0.2, 0.5), (0.5, 1.0), (1.0, 2.0)] {
        let approach_m = x * length;
        let normal_force_n = y * force;
        let energy = nodes.last().map_or(0.0, |old| {
            old.curve_reversible_energy_j
                + (0.5 * old.normal_force_n + 0.5 * normal_force_n) * (approach_m - old.approach_m)
        });
        nodes.push(FiniteGapResponseNode {
            approach_m,
            contact_identity: identity,
            normal_force_n,
            peak_pressure_pa: 0.0,
            pressure_centroid_m: [0.0; 2],
            equivalent_pressure_semiaxes_m: [0.0; 2],
            reference_reversible_energy_j: energy,
            curve_reversible_energy_j: energy,
            energy_residual_j: 0.0,
        });
    }
    FiniteGapResponseCurve {
        identity,
        gap_source_identity: identity,
        gap_source_authority: FiniteGapChartSamplingAuthority::Estimate,
        response_authority: FiniteGapResponseCurveAuthority::Estimate,
        nodes,
        maximum_absolute_energy_residual_j: 0.0,
    }
}

pub fn close(actual: f64, expected: f64, scale: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * scale,
        "{actual:e} != {expected:e}; scaled error {:e}",
        (actual - expected).abs() / scale
    );
}
