//! G1/G3/G4: original-source rough adhesion, independent sums and branch memory.
#[path = "source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrEnvelope, JkrJumpKind, SourceJkrSphereContact, asperity::*,
};
use fs_exec::CancelGate;
use fs_tribo::InputAuthority;

fn source() -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, 1e7);
    let b = common::solid("B", 2000.0, 2e7);
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
fn envelope(h: f64) -> JkrAsperityEnvelope {
    JkrAsperityEnvelope {
        minimum_separation_m: -4.0 * h,
        maximum_separation_m: 10.0 * h,
        nominal_area_m2: 1e-3,
        maximum_contact_area_fraction: 0.05,
    }
}
fn bin(source: &SourceJkrSphereContact, height_m: f64, multiplicity: f64) -> JkrAsperity<'_> {
    JkrAsperity {
        source,
        height_m,
        multiplicity,
    }
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn single_asperity_retains_original_area_force_work_and_all_branch_events() {
    let law = source();
    let h = law.indentation_scale_m();
    let z = 2.0 * h;
    let gate = CancelGate::new();
    let population =
        SourceJkrAsperityContact::new([bin(&law, z, 1.0)], envelope(h), &gate).unwrap();
    let mut state = population.detached_state(z + 4.0 * h, &gate).unwrap();
    let mut original = law.detached_state(z - (z + 4.0 * h)).unwrap();
    for gap in [
        z + h,
        z,
        z - 2.0 * h,
        z + h,
        z + 4.0 * h,
        z + 5.0 * h,
        z - h,
    ] {
        let trial = population.advance_monotone(&state, gap, &gate).unwrap();
        let expected = law.advance_monotone(original, z - gap).unwrap();
        assert_eq!(
            trial.response().normal_force_n,
            expected.response.normal_force_n
        );
        assert_eq!(
            trial.response().contact_area_m2,
            expected.response.contact_area_m2
        );
        assert_eq!(
            trial.response().reversible_energy_j,
            expected.response.reversible_energy_j
        );
        assert_eq!(trial.prescribed_work_j, expected.prescribed_work_j);
        assert_eq!(trial.released_energy_j, expected.irreversible_work_j);
        assert_eq!(trial.state().local_states()[0], expected.response.state());
        assert_eq!(trial.response().authority, InputAuthority::SyntheticFixture);
        assert_eq!(
            trial.jumps[0].map(|j| j.kind),
            expected.jump.map(|j| j.kind)
        );
        if let (Some(jump), Some(original_jump)) = (trial.jumps[0], expected.jump) {
            assert_eq!(jump.separation_m, z - original_jump.indentation_m);
            assert_eq!(jump.released_energy_j, original_jump.released_energy_j);
        }
        state = trial.state();
        original = expected.response.state();
    }
}

#[test]
fn rough_summits_attach_and_detach_at_distinct_physical_heights() {
    let law = source();
    let h = law.indentation_scale_m();
    let gate = CancelGate::new();
    let population = SourceJkrAsperityContact::new(
        [
            bin(&law, -2.0 * h, 1.0),
            bin(&law, 0.0, 2.0),
            bin(&law, 2.0 * h, 3.0),
        ],
        envelope(h),
        &gate,
    )
    .unwrap();
    let initial = population.detached_state(7.0 * h, &gate).unwrap();
    let first = population.advance_monotone(&initial, h, &gate).unwrap();
    assert_eq!(
        first.state().local_states().map(|s| s.branch()),
        [
            JkrBranch::Detached,
            JkrBranch::Detached,
            JkrBranch::Attached
        ]
    );
    assert_eq!(first.response().attached_multiplicity, 3.0);
    assert!(first.response().normal_force_n < 0.0);
    assert!(first.response().contact_area_m2 > 0.0);
    assert_eq!(first.jumps[2].unwrap().separation_m, 2.0 * h);
    assert_eq!(first.jumps[2].unwrap().kind, JkrJumpKind::Attach);
    let partial_unload = population
        .advance_monotone(&first.state(), 2.5 * h, &gate)
        .unwrap();
    assert_eq!(partial_unload.response().attached_multiplicity, 3.0);
    assert!(partial_unload.jumps.iter().all(Option::is_none));
    // The middle summit sees negative approach throughout these two phases;
    // its neighbour's adhesive branch must not nucleate a bond there.
    assert_eq!(
        partial_unload.state().local_states()[1].branch(),
        JkrBranch::Detached
    );

    let loaded = population
        .advance_monotone(&initial, -3.0 * h, &gate)
        .unwrap();
    assert_eq!(loaded.response().attached_multiplicity, 6.0);
    let mut state = loaded.state();
    let mut work = loaded.prescribed_work_j;
    let mut release = loaded.released_energy_j;
    for (step_index, gap) in [2.0 * h, 4.0 * h, 6.0 * h].into_iter().enumerate() {
        let step = population.advance_monotone(&state, gap, &gate).unwrap();
        let event = step.jumps[step_index].unwrap();
        assert_eq!(event.kind, JkrJumpKind::Detach);
        close(
            event.separation_m,
            [-2.0 * h, 0.0, 2.0 * h][step_index] - law.minimum_attached_indentation_m(),
            h,
            2e-15,
        );
        work += step.prescribed_work_j;
        release += step.released_energy_j;
        state = step.state();
    }
    assert_eq!(state.response().attached_multiplicity, 0.0);
    assert_eq!(state.response().reversible_energy_j, 0.0);
    assert!(release > 0.0);
    close(work, release, release, 2e-13);
    // A prescribed monotone overshoot contains the same per-asperity jumps and work.
    let direct = population
        .advance_monotone(&loaded.state(), 6.0 * h, &gate)
        .unwrap();
    assert_eq!(direct.state(), state);
    close(
        loaded.prescribed_work_j + direct.prescribed_work_j,
        work,
        release,
        2e-13,
    );
    close(
        loaded.released_energy_j + direct.released_energy_j,
        release,
        release,
        2e-13,
    );
}

#[test]
fn distinct_curvatures_sum_dimensional_jkr_force_area_energy_and_tangent() {
    let law = source();
    let second =
        SourceJkrSphereContact::new(law.normal_state().clone(), 0.02, 1e-9, law.envelope())
            .unwrap();
    let h = law.indentation_scale_m();
    let gate = CancelGate::new();
    let population = SourceJkrAsperityContact::new(
        [bin(&law, 6.0 * h, 1.0), bin(&second, -h, 5.0)],
        envelope(h),
        &gate,
    )
    .unwrap();
    let state = population.attached_state(0.0, &gate).unwrap();
    let mut expected_force = 0.0;
    let mut expected_area = 0.0;
    let mut expected_repulsive = 0.0;
    let mut expected_energy = 0.0;
    let mut expected_tangent = 0.0;
    for summit in population.asperities() {
        let a = summit
            .source
            .evaluate_attached(summit.height_m)
            .unwrap()
            .contact_radius_m;
        let e = summit.source.normal_state().elastic().reduced_modulus_pa();
        let w = summit.source.normal_state().adhesion_energy_j_per_m2();
        let r = summit.source.effective_radius_m();
        let c = (2.0 * PI * w / e).sqrt();
        let repulsive = 4.0 * e * a.powi(3) / (3.0 * r);
        let force = repulsive - (8.0 * PI * w * e * a.powi(3)).sqrt();
        let energy = 8.0 * e * a.powi(5) / (15.0 * r * r) - 4.0 * e * c * a.powf(3.5) / (3.0 * r)
            + PI * w * a * a;
        let tangent =
            (4.0 * e * a * a / r - 3.0 * e * c * a.sqrt()) / (2.0 * a / r - c / (2.0 * a.sqrt()));
        expected_force += summit.multiplicity * force;
        expected_area += summit.multiplicity * PI * a * a;
        expected_repulsive += summit.multiplicity * repulsive;
        expected_energy += summit.multiplicity * energy;
        expected_tangent += summit.multiplicity * tangent;
    }
    let response = state.response();
    close(
        response.normal_force_n,
        expected_force,
        expected_repulsive,
        1e-13,
    );
    close(
        response.contact_area_m2,
        expected_area,
        expected_area,
        1e-13,
    );
    close(
        response.reversible_energy_j,
        expected_energy,
        expected_repulsive * h,
        1e-12,
    );
    close(
        response.tangent_n_per_m.unwrap(),
        expected_tangent,
        expected_repulsive / h,
        1e-13,
    );
    assert!(law.evaluate_attached(6.0 * h).unwrap().normal_force_n > 0.0);
    assert!(second.evaluate_attached(-h).unwrap().normal_force_n < 0.0);
    assert!(response.normal_force_n < 0.0);
    let friction = response.friction_contact_state().unwrap();
    close(friction.area_m2(), expected_area, expected_area, 1e-13);
    close(
        friction.repulsive_normal_force_n(),
        expected_repulsive,
        expected_repulsive,
        1e-13,
    );
    let dx = 1e-4 * h;
    let left = population.attached_state(-dx, &gate).unwrap();
    let right = population.attached_state(dx, &gate).unwrap();
    close(
        -(right.response().reversible_energy_j - left.response().reversible_energy_j) / (2.0 * dx),
        response.normal_force_n,
        expected_repulsive,
        2e-8,
    );
    close(
        -(right.response().normal_force_n - left.response().normal_force_n) / (2.0 * dx),
        response.tangent_n_per_m.unwrap(),
        expected_repulsive / h,
        2e-8,
    );
}

#[test]
fn bin_partition_and_height_reference_preserve_population_physics() {
    let law = source();
    let h = law.indentation_scale_m();
    let gate = CancelGate::new();
    let original = SourceJkrAsperityContact::new([bin(&law, h, 5.0)], envelope(h), &gate).unwrap();
    let split =
        SourceJkrAsperityContact::new([bin(&law, h, 2.0), bin(&law, h, 3.0)], envelope(h), &gate)
            .unwrap();
    let shift = 32.0 * h;
    let mut shifted_envelope = envelope(h);
    shifted_envelope.minimum_separation_m += shift;
    shifted_envelope.maximum_separation_m += shift;
    let shifted =
        SourceJkrAsperityContact::new([bin(&law, h + shift, 5.0)], shifted_envelope, &gate)
            .unwrap();
    let a = original.detached_state(7.0 * h, &gate).unwrap();
    let b = split.detached_state(7.0 * h, &gate).unwrap();
    let c = shifted.detached_state(7.0 * h + shift, &gate).unwrap();
    let x = original.advance_monotone(&a, -h, &gate).unwrap();
    let y = split.advance_monotone(&b, -h, &gate).unwrap();
    let z = shifted.advance_monotone(&c, -h + shift, &gate).unwrap();
    for response in [y.response(), z.response()] {
        close(
            response.normal_force_n,
            x.response().normal_force_n,
            x.response().repulsive_normal_force_n,
            2e-13,
        );
        close(
            response.contact_area_m2,
            x.response().contact_area_m2,
            x.response().contact_area_m2,
            2e-13,
        );
    }
    close(
        y.prescribed_work_j,
        x.prescribed_work_j,
        x.released_energy_j,
        2e-13,
    );
    close(
        z.prescribed_work_j,
        x.prescribed_work_j,
        x.released_energy_j,
        2e-13,
    );
    assert_eq!(x.released_energy_j, z.released_energy_j);
}

#[test]
fn refusal_and_cancellation_preserve_accepted_branches_and_original_source() {
    let law = source();
    let h = law.indentation_scale_m();
    let gate = CancelGate::new();
    let population =
        SourceJkrAsperityContact::new([bin(&law, 0.0, 1.0)], envelope(h), &gate).unwrap();
    let initial = population.detached_state(5.0 * h, &gate).unwrap();
    let expected = population.advance_monotone(&initial, -h, &gate).unwrap();
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        population.advance_monotone(&initial, -h, &cancelled),
        Err(JkrAsperityError::Cancelled)
    );
    assert_eq!(
        population.advance_monotone(&initial, 11.0 * h, &gate),
        Err(JkrAsperityError::OutsideSeparation)
    );
    let other = SourceJkrAsperityContact::new([bin(&law, h, 1.0)], envelope(h), &gate).unwrap();
    assert_eq!(
        other.advance_monotone(&initial, -h, &gate),
        Err(JkrAsperityError::StateMismatch)
    );
    assert_eq!(initial.response().attached_multiplicity, 0.0);
    assert_eq!(
        population.advance_monotone(&initial, -h, &gate).unwrap(),
        expected
    );
    let at_touch = population.advance_monotone(&initial, 0.0, &gate).unwrap();
    assert_eq!(at_touch.response().tangent_n_per_m, None);
    assert_eq!(at_touch.jumps[0].unwrap().kind, JkrJumpKind::Attach);

    let empty = SourceJkrAsperityContact::<0>::new([], envelope(h), &gate);
    assert!(matches!(empty, Err(JkrAsperityError::Invalid(_))));
    let oversized = SourceJkrAsperityContact::new(
        [bin(&law, 0.0, 1.0); MAX_JKR_ASPERITY_BINS + 1],
        envelope(h),
        &gate,
    );
    assert!(matches!(oversized, Err(JkrAsperityError::Invalid(_))));
    assert!(matches!(
        SourceJkrAsperityContact::new([bin(&law, 0.0, 0.0)], envelope(h), &gate),
        Err(JkrAsperityError::Invalid(_))
    ));
    let mut too_small = envelope(h);
    too_small.nominal_area_m2 = 1e-12;
    assert!(matches!(
        SourceJkrAsperityContact::new([bin(&law, 0.0, 1.0)], too_small, &gate),
        Err(JkrAsperityError::ContactAreaFraction { .. })
    ));
    let mut too_deep = envelope(h);
    too_deep.minimum_separation_m = -1.0;
    assert!(matches!(
        SourceJkrAsperityContact::new([bin(&law, 0.0, 1.0)], too_deep, &gate),
        Err(JkrAsperityError::Source { index: 0, .. })
    ));
}
