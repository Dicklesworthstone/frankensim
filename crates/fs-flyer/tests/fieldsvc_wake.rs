//! Physical near/far wake fields through the original native field service.

use fs_flyer::fieldsvc::{
    BoundSystem, C_BOUND_CIRCULATION, C_GROUND_IMAGES, C_MEAN_ATMO, C_PHYSICAL_WAKE,
    FieldSourceStateV1, GridSpec, PreparedWakeField, ambient_atmosphere_v1, claim_total_flow,
    sample_field,
};
use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, hybrid_kinematics};
use fs_vpm::filament3d::{FilamentWake, WakeRateCertificate};

fn fixture() -> (FilamentWake, FarField) {
    let mut near = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 120.0,
            n_stations: 2,
            max_rows: 12,
        },
        vec![[0.0, -1.0, 2.0], [0.0, 0.0, 2.1], [0.0, 1.0, 2.0]],
    )
    .unwrap();
    for i in 0..6 {
        near.shed(
            &[0.8 + 0.1 * f64::from(i), 0.6 - 0.02 * f64::from(i)],
            [-0.3, 0.02, -0.01],
        )
        .unwrap();
        near.rows.last_mut().unwrap().core2_m2 = 0.01 + 0.002 * f64::from(i);
    }
    let far = FarField {
        mode: WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.015 },
        shed_dt_s: 1.0 / 120.0,
        cells: Vec::new(),
    };
    (near, far)
}

fn state() -> FieldSourceStateV1 {
    FieldSourceStateV1 {
        tick: 6,
        source_state_digest: "physical-wake-fixture".into(),
        atmosphere: ambient_atmosphere_v1(1903, 8.0).unwrap(),
        // Deliberately unrelated legacy source: the prepared connected wake
        // must own its bound line instead of adding this second horseshoe.
        bound: Some(BoundSystem {
            gamma_m2ps: 1000.0,
            tip_left_m: [0.0, -4.0, 3.0],
            tip_right_m: [0.0, 4.0, 3.0],
            trail_m: 20.0,
            core_m: 0.2,
        }),
        images_active: true,
    }
}

fn grid(p: [f64; 3]) -> GridSpec {
    GridSpec {
        origin_m: p,
        dx_m: 0.1,
        nx: 1,
        ny: 1,
        nz: 1,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn physical_wake_survives_aggregation_and_owns_its_bound_line() {
    let (mut near, far) = fixture();
    let whole = PreparedWakeField::new(&near, &far, false, 1000, || false).unwrap();
    let far = FarField::aggregate(&mut near, 4, 2, far.mode, far.shed_dt_s).unwrap();
    let partitioned = PreparedWakeField::new(&near, &far, false, 1000, || false).unwrap();
    let state = state();
    let grid = GridSpec {
        nx: 3,
        ny: 2,
        ..grid([-0.73, 0.24, 1.72])
    };
    let a = whole
        .sample(&state, &grid, C_PHYSICAL_WAKE, 100000, || false)
        .unwrap();
    let b = partitioned
        .sample(&state, &grid, C_PHYSICAL_WAKE, 100000, || false)
        .unwrap();
    for n in 0..grid.n_points() {
        let direct = hybrid_kinematics(&near, &far, grid.point(n)).unwrap();
        assert!(a.validity_mask[n] && b.validity_mask[n]);
        for i in 0..3 {
            close(a.u[n][i], b.u[n][i], 2e-13);
            close(b.u[n][i], direct.velocity[i], 2e-13);
            close(a.omega[n][i], b.omega[n][i], 2e-13);
            for j in 0..3 {
                close(a.grad_u[n][i][j], b.grad_u[n][i][j], 2e-13);
                close(b.grad_u[n][i][j], direct.gradient[i][j], 2e-13);
            }
        }
    }
    assert_eq!(b.meta.source_modes & C_BOUND_CIRCULATION, 0);
    assert_ne!(b.meta.source_modes & C_PHYSICAL_WAKE, 0);
    claim_total_flow(&b).unwrap();
    assert_eq!(
        partitioned
            .sample(
                &state,
                &grid,
                C_PHYSICAL_WAKE | C_BOUND_CIRCULATION,
                100000,
                || false
            )
            .unwrap_err()
            .code,
        "component-unsupported"
    );
    let ambient = partitioned
        .sample(&state, &grid, C_MEAN_ATMO, 0, || false)
        .unwrap();
    assert_eq!(
        claim_total_flow(&ambient).unwrap_err().code,
        "forbidden-claim-total-flow"
    );
    let original = sample_field(&state, &grid, C_MEAN_ATMO).unwrap();
    assert_eq!(ambient.u, original.u);
    assert_eq!(ambient.grad_u, original.grad_u);
}

#[test]
fn complete_wake_images_obey_plane_condition_and_axial_curl_reflection() {
    let (mut near, far) = fixture();
    let far = FarField::aggregate(&mut near, 4, 2, far.mode, far.shed_dt_s).unwrap();
    let field = PreparedWakeField::new(&near, &far, true, 1000, || false).unwrap();
    let state = state();
    let grid = grid([-0.77, 0.3, 0.0]);
    let physical = field
        .sample(&state, &grid, C_PHYSICAL_WAKE, 100000, || false)
        .unwrap();
    let image = field
        .sample(&state, &grid, C_GROUND_IMAGES, 100000, || false)
        .unwrap();
    let total = field
        .sample(
            &state,
            &grid,
            C_PHYSICAL_WAKE | C_GROUND_IMAGES,
            100000,
            || false,
        )
        .unwrap();
    for i in 0..3 {
        let sign = if i == 2 { -1.0 } else { 1.0 };
        close(image.u[0][i], sign * physical.u[0][i], 2e-14);
        close(image.omega[0][i], -sign * physical.omega[0][i], 2e-14);
    }
    close(total.u[0][2], 0.0, 2e-14);
    close(total.grad_u[0][2][0], 0.0, 2e-14);
    close(total.grad_u[0][2][1], 0.0, 2e-14);
    close(total.div_analytic[0], total.div_finite_difference[0], 2e-7);
    assert!(total.validity_mask[0]);
    assert!(
        claim_total_flow(&physical).is_err(),
        "active image flow cannot be omitted from total flow"
    );
    claim_total_flow(&total).unwrap();
}

#[test]
fn independent_velocity_probes_recover_full_wake_gradient() {
    let (mut near, far) = fixture();
    let far = FarField::aggregate(&mut near, 4, 2, far.mode, far.shed_dt_s).unwrap();
    let field = PreparedWakeField::new(&near, &far, true, 1000, || false).unwrap();
    let state = state();
    let p = [-0.55, 0.21, 1.6];
    let mask = C_PHYSICAL_WAKE | C_GROUND_IMAGES | C_MEAN_ATMO;
    let center = field
        .sample(&state, &grid(p), mask, 100000, || false)
        .unwrap();
    for axis in 0..3 {
        let h = 1e-5;
        let mut plus = p;
        let mut minus = p;
        plus[axis] += h;
        minus[axis] -= h;
        let high = field
            .sample(&state, &grid(plus), mask, 100000, || false)
            .unwrap();
        let low = field
            .sample(&state, &grid(minus), mask, 100000, || false)
            .unwrap();
        for component in 0..3 {
            close(
                (high.u[0][component] - low.u[0][component]) / (2.0 * h),
                center.grad_u[0][component][axis],
                2e-8,
            );
        }
    }
}

#[test]
fn finite_core_axis_keeps_vorticity_and_zero_core_is_invalid() {
    let (near, mut far) = fixture();
    // Use an exactly represented endpoint; a rounded midpoint of a tilted
    // segment need not lie exactly on its zero-radius mathematical axis.
    let p = near.line_nodes[0];
    let field = PreparedWakeField::new(&near, &far, false, 1000, || false).unwrap();
    let sample = field
        .sample(&state(), &grid(p), C_PHYSICAL_WAKE, 100000, || false)
        .unwrap();
    assert!(sample.validity_mask[0]);
    assert!(sample.singularity_core_mask[0]);
    assert!(sample.omega[0].iter().map(|w| w * w).sum::<f64>() > 1.0);
    let mut zero = near.clone();
    for row in &mut zero.rows {
        row.core2_m2 = 0.0;
    }
    far.mode = WakeCoreEvolutionMode::Frozen;
    let field = PreparedWakeField::new(&zero, &far, false, 1000, || false).unwrap();
    let sample = field
        .sample(&state(), &grid(p), C_PHYSICAL_WAKE, 100000, || false)
        .unwrap();
    assert!(!sample.validity_mask[0]);
    assert!(sample.normalized_divergence(0).is_none());
}

#[test]
fn admission_cancellation_and_snapshot_identity_preserve_sources() {
    let (mut near, mut far) = fixture();
    let before = near.digest();
    let field = PreparedWakeField::new(&near, &far, true, 1000, || false).unwrap();
    let count = field.segment_count();
    assert_eq!(
        PreparedWakeField::new(&near, &far, true, count - 1, || false)
            .unwrap_err()
            .code,
        "field-wake-budget"
    );
    let state = state();
    let grid = grid([0.3, 0.2, 1.6]);
    let visits = 16 * count as u64;
    let mask = C_PHYSICAL_WAKE | C_GROUND_IMAGES;
    assert_eq!(
        field
            .sample(&state, &grid, mask, visits - 1, || false)
            .unwrap_err()
            .code,
        "field-wake-budget"
    );
    let first = field.sample(&state, &grid, mask, visits, || false).unwrap();
    let mut polls = 0;
    assert_eq!(
        field
            .sample(&state, &grid, mask, visits, || {
                polls += 1;
                polls == 6
            })
            .unwrap_err()
            .code,
        "field-cancelled"
    );
    let replay = field.sample(&state, &grid, mask, visits, || false).unwrap();
    assert_eq!(first.digest(), replay.digest());
    assert_eq!(near.digest(), before);
    near.rows[0].core2_m2 *= 2.0;
    let changed = PreparedWakeField::new(&near, &far, true, 1000, || false).unwrap();
    let changed = changed
        .sample(&state, &grid, mask, 100000, || false)
        .unwrap();
    assert_ne!(
        first.meta.field_source_snapshot_id,
        changed.meta.field_source_snapshot_id
    );
    assert_eq!(
        first.digest(),
        field
            .sample(&state, &grid, mask, visits, || false)
            .unwrap()
            .digest()
    );
    far.shed_dt_s *= 2.0;
    assert_eq!(
        PreparedWakeField::new(&near, &far, true, 1000, || false)
            .unwrap_err()
            .code,
        "field-wake-invalid"
    );
    let overflow = GridSpec {
        nx: usize::MAX,
        ny: 2,
        ..grid
    };
    assert!(
        field
            .sample(&state, &overflow, mask, u64::MAX, || false)
            .is_err()
    );
}

#[test]
fn overflowing_vortex_invariants_never_mark_a_sample_valid() {
    let (mut near, far) = fixture();
    for row in &mut near.rows {
        for gamma in &mut row.gamma {
            *gamma *= 1e160;
        }
    }
    let field = PreparedWakeField::new(&near, &far, false, 1000, || false).unwrap();
    let sample = field
        .sample(
            &state(),
            &grid([-0.55, 0.21, 1.6]),
            C_PHYSICAL_WAKE,
            100000,
            || false,
        )
        .unwrap();
    assert!(sample.u[0].iter().all(|v| v.is_finite()));
    assert!(sample.grad_u[0].iter().flatten().all(|v| v.is_finite()));
    assert!(!sample.validity_mask[0]);
    assert!(sample.normalized_divergence(0).is_none());
}
