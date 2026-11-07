//! G0/G1/G3: original JKR area, tensile rolling memory and independent work.
#[path = "source_point/common.rs"]
mod common;
use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::rolling::*;
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_tribo::InputAuthority;

fn normal(label: &str) -> SourceJkrSphereContact {
    let a = common::solid(label, 1000.0, 1e7);
    let b = common::solid("rolling-B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
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
fn controls() -> JkrRollingControls {
    JkrRollingControls {
        critical_rolling_displacement_m: 1e-8,
        maximum_critical_displacement_to_contact_radius: 0.01,
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn original_jkr_zero_and_tensile_loads_have_nonzero_rolling_stiffness() {
    let normal = normal("rolling-A");
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let r = normal.effective_radius_m();
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let a0 = (4.5 * PI * w * r * r / e).cbrt();
    let fc = 1.5 * PI * w * r;
    let mut stiffness = Vec::new();
    for load in [0.0, -0.5 * fc, -fc] {
        let response = normal.evaluate_force_controlled(load).unwrap();
        let source =
            SourceJkrRollingContact::new(normal.clone(), response.indentation_m, controls())
                .unwrap();
        let factor = (response.contact_radius_m / a0).powf(1.5);
        let k = 4.0 * fc / r * factor;
        stiffness.push(source.law().stiffness_n_per_m());
        close(source.law().stiffness_n_per_m(), k);
        assert!(
            (source.normal_response().normal_force_n - response.normal_force_n).abs() < 1e-12 * fc
        );
        assert_eq!(
            source.normal_response().normal_state_identity(),
            normal.normal_state().identity()
        );
        assert_eq!(
            source.normal_response().authority,
            InputAuthority::SyntheticFixture
        );
        let angle = 0.25 * controls().critical_rolling_displacement_m / r;
        let step = source.advance(&source.initial_state(), angle).unwrap();
        close(step.rolling().endpoint_moment_n_m, -k * r * r * angle);
        close(
            step.rolling().storage_after_j,
            0.5 * k * (r * angle).powi(2),
        );
        close(
            step.rolling().applied_work_j,
            step.rolling().storage_after_j,
        );
        assert_eq!(step.rolling().irreversible_loss_j, 0.0);
        assert!(step.rolling().endpoint_moment_n_m < 0.0);
    }
    close(stiffness[2] / stiffness[0], 0.5);
}

#[test]
fn source_owned_yield_reversal_and_partition_preserve_distinct_storage_and_loss() {
    let normal = normal("rolling-A");
    let response = normal.evaluate_force_controlled(0.0).unwrap();
    let source = SourceJkrRollingContact::new(normal, response.indentation_m, controls()).unwrap();
    let theta = source.law().critical_rolling_angle_rad();
    let initial = source.initial_state();
    let full = source.advance(&initial, 3.0 * theta).unwrap();
    let k = source.law().stiffness_n_per_m();
    let xi = controls().critical_rolling_displacement_m;
    close(full.rolling().irreversible_loss_j, 2.0 * k * xi * xi);
    close(full.rolling().storage_after_j, 0.5 * k * xi * xi);
    close(
        -full.rolling().mean_moment_n_m * 3.0 * theta,
        full.rolling().applied_work_j,
    );
    let mut state = initial;
    let (mut work, mut loss) = (0.0, 0.0);
    for increment in [0.25, 0.75, 0.5, 1.5] {
        let step = source.advance(&state, increment * theta).unwrap();
        work += step.rolling().applied_work_j;
        loss += step.rolling().irreversible_loss_j;
        state = step.state();
    }
    close(state.elastic_shift_m(), full.state().elastic_shift_m());
    close(work, full.rolling().applied_work_j);
    close(loss, full.rolling().irreversible_loss_j);
    let unload = source.advance(&state, -theta).unwrap();
    assert_eq!(unload.rolling().irreversible_loss_j, 0.0);
    assert!(unload.rolling().applied_work_j < 0.0);
    close(
        source.stored_energy_j(&state).unwrap(),
        full.rolling().storage_after_j,
    );
    let reverse = source.advance(&state, -3.0 * theta).unwrap();
    close(reverse.state().elastic_shift_m(), -xi);
    close(reverse.rolling().irreversible_loss_j, k * xi * xi);
}

#[test]
fn changed_normal_state_or_source_cannot_discard_or_rebind_rolling_memory() {
    let original = normal("rolling-A");
    let indentation = original
        .evaluate_force_controlled(0.0)
        .unwrap()
        .indentation_m;
    let source = SourceJkrRollingContact::new(original.clone(), indentation, controls()).unwrap();
    let loaded = source
        .advance(
            &source.initial_state(),
            0.5 * source.law().critical_rolling_angle_rad(),
        )
        .unwrap()
        .state();
    let saved = loaded;
    let changed = [
        SourceJkrRollingContact::new(original.clone(), indentation * 1.1, controls()).unwrap(),
        SourceJkrRollingContact::new(normal("other-source-A"), indentation, controls()).unwrap(),
        SourceJkrRollingContact::new(
            original,
            indentation,
            JkrRollingControls {
                critical_rolling_displacement_m: 2e-8,
                ..controls()
            },
        )
        .unwrap(),
    ];
    for other in changed {
        assert!(matches!(
            other.advance(&loaded, 0.0),
            Err(SourceJkrRollingError::StateMismatch)
        ));
        assert!(other.stored_energy_j(&loaded).is_err());
    }
    assert!(source.advance(&loaded, f64::NAN).is_err());
    assert_eq!(loaded, saved);
    let no_op = source.advance(&loaded, 0.0).unwrap();
    assert_eq!(no_op.state(), loaded);
    assert_eq!(no_op.rolling().applied_work_j, 0.0);
    assert_eq!(no_op.rolling().irreversible_loss_j, 0.0);
}
