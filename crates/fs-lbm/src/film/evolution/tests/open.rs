//! Independent open-film balances, including both physical half-cell profiles.
use super::*;

fn ends(left: f64, right: f64, shear: f64) -> ThinFilmOpenEnds {
    ThinFilmOpenEnds {
        left: ThinFilmOpenEnd {
            pressure_pa: left,
            surface_shear_pa: shear,
        },
        right: ThinFilmOpenEnd {
            pressure_pa: right,
            surface_shear_pa: shear,
        },
    }
}
fn opening_from_endpoint(h: &[f64], q: &[f64], dt: f64) -> Vec<f64> {
    (0..h.len())
        .map(|i| h[i] + dt * (q[i + 1] - q[i]) / 0.5)
        .collect()
}
// Direct central differences with natural zero-slope ghosts, independently
// evaluated at the authored endpoint rather than inferred from solver reports.
fn pressures(h: &[f64]) -> Vec<f64> {
    (0..h.len())
        .map(|i| {
            let left = h[i.saturating_sub(1)];
            let right = h[(i + 1).min(h.len() - 1)];
            3.0 * h[i] - 0.04 * (left - 2.0 * h[i] + right) / (0.5 * 0.5)
        })
        .collect()
}

#[test]
fn uniform_open_throughflow_preserves_height_and_counts_both_half_cell_powers() {
    let h: f64 = 0.1;
    let (n, dx, dt, mu) = (4, 0.5, 0.4, 0.5);
    for sign in [-1.0, 1.0] {
        let (acceleration, shear) = (0.3 * sign, 0.03 * sign);
        let body = 2.0 * acceleration;
        let q = (body * h.powi(3) / 3.0 + shear * h * h / 2.0) / mu;
        let surface = (body * h * h / 2.0 + shear * h) / mu;
        let dissipation =
            (body * body * h.powi(3) / 3.0 + body * shear * h * h + shear * shear * h) / mu;
        let initial = vec![h; n];
        let ends = ends(3.0 * h, 3.0 * h, shear);
        let mut line = ThinFilmLine::newtonian(
            material(0.04),
            geometry(ThinFilmBoundary::Open, 1.5),
            initial.clone(),
            0.0,
        )
        .unwrap();
        assert_eq!(line.face_count(), n - 1);
        assert_eq!(line.profile_source_count(), n + 1);
        let report = line
            .advance_open(
                dt,
                acceleration,
                &vec![shear; n - 1],
                ends,
                limits(),
                controls(),
                &initial,
                || false,
            )
            .unwrap();
        let open = report.open_ends.unwrap();
        assert_eq!(open.ends, ends);
        for profile in [open.left, open.right]
            .iter()
            .chain(report.faces.iter().map(|f| &f.profile))
        {
            near(profile.thickness_m, h, 2e-12);
            near(profile.flow_per_width_m2_s, q, 2e-13);
            near(profile.surface_velocity_m_s, surface, 2e-13);
            near(profile.viscous_dissipation_w_m2, dissipation, 2e-13);
        }
        for &actual in line.heights_m() {
            near(actual, h, 2e-12);
        }
        near(
            report.viscous_heat_j_per_m,
            dt * n as f64 * dx * dissipation,
            2e-12,
        );
        near(
            report.body_work_j_per_m,
            dt * n as f64 * dx * body * q,
            2e-12,
        );
        near(
            report.surface_work_j_per_m,
            dt * n as f64 * dx * shear * surface,
            2e-12,
        );
        near(report.boundary_work_j_per_m, 0.0, 2e-13);
        near(report.volume_change_m2, 0.0, 1e-11);
        near(report.volume_balance_residual_m2, 0.0, 1e-11);
        near(report.energy_residual_j_per_m, 0.0, 1e-11);
    }
}

#[test]
fn manufactured_open_endpoints_recover_signed_volume_capillary_pressure_and_reservoir_work() {
    let target = [0.10, 0.11, 0.09, 0.12];
    let pressure = pressures(&target);
    let (dx, dt, mu) = (0.5, 1.0, 0.5);
    let mut filling = false;
    let mut draining = false;
    let mut exporting_work = false;
    for (acceleration, left, right) in [(0.2, 0.44, 0.20), (-0.2, 0.14, 0.48), (2.0, 0.10, 0.80)] {
        let body = 2.0 * acceleration;
        let forces = [
            body - (pressure[0] - left) / (dx / 2.0),
            body - (pressure[1] - pressure[0]) / dx,
            body - (pressure[2] - pressure[1]) / dx,
            body - (pressure[3] - pressure[2]) / dx,
            body - (right - pressure[3]) / (dx / 2.0),
        ];
        let heights: [f64; 5] = std::array::from_fn(|i| {
            if i == 0 {
                target[0]
            } else if i == 4 {
                target[3]
            } else if forces[i] >= 0.0 {
                target[i - 1]
            } else {
                target[i]
            }
        });
        let q: [f64; 5] = std::array::from_fn(|i| forces[i] * heights[i].powi(3) / (3.0 * mu));
        let old = opening_from_endpoint(&target, &q, dt);
        let mut line = ThinFilmLine::newtonian(
            material(0.04),
            geometry(ThinFilmBoundary::Open, 1.5),
            old.clone(),
            0.0,
        )
        .unwrap();
        let report = line
            .advance_open(
                dt,
                acceleration,
                &[0.0; 3],
                ends(left, right, 0.0),
                limits(),
                controls(),
                &old,
                || false,
            )
            .unwrap();
        for (&actual, &expected) in line.heights_m().iter().zip(&target) {
            near(actual, expected, 3e-11);
        }
        let open = report.open_ends.unwrap();
        near(open.left.flow_per_width_m2_s, q[0], 3e-11);
        near(open.right.flow_per_width_m2_s, q[4], 3e-11);
        let transported = dt * (q[0] - q[4]);
        let work = dt * (left * q[0] - right * q[4]);
        let heat = dt
            * dx
            * (0..5)
                .map(|i| {
                    let weight = if i == 0 || i == 4 { 0.5 } else { 1.0 };
                    weight * forces[i] * forces[i] * heights[i].powi(3) / (3.0 * mu)
                })
                .sum::<f64>();
        let body_work = dt * dx * body * (0.5 * q[0] + q[1] + q[2] + q[3] + 0.5 * q[4]);
        let delta: Vec<_> = target.iter().zip(&old).map(|(&a, &b)| a - b).collect();
        let change =
            energy(&target, 2.0, 1.5, 0.04, dx, false) - energy(&old, 2.0, 1.5, 0.04, dx, false);
        let defect = energy(&delta, 2.0, 1.5, 0.04, dx, false);
        near(report.volume_change_m2, transported, 1e-11);
        near(dx * delta.iter().sum::<f64>(), transported, 2e-15);
        near(report.volume_balance_residual_m2, 0.0, 1e-11);
        near(report.boundary_work_j_per_m, work, 1e-11);
        near(report.viscous_heat_j_per_m, heat, 1e-11);
        near(report.body_work_j_per_m, body_work, 1e-11);
        near(report.stored_energy_change_j_per_m, change, 1e-11);
        near(report.backward_euler_energy_defect_j_per_m, defect, 1e-11);
        near(change + defect + heat - body_work - work, 0.0, 3e-15);
        near(report.energy_residual_j_per_m, 0.0, 1e-11);
        filling |= transported > 0.0;
        draining |= transported < 0.0;
        exporting_work |= work < 0.0;
    }
    assert!(filling && draining && exporting_work);
}

#[test]
fn original_carreau_open_profiles_keep_separate_end_and_interior_temperature_sources() {
    let law = source(0.4);
    let temperatures = [330.0, 300.0, 310.0, 320.0, 340.0];
    let target = [0.1; 4];
    let (dt, body, left, right) = (0.02, 0.6, 0.4, 0.2);
    let forces = [
        body - (0.3 - left) / 0.25,
        body,
        body,
        body,
        body - (right - 0.3) / 0.25,
    ];
    let reference: [[f64; 3]; 5] =
        std::array::from_fn(|i| original_profile(law, temperatures[i], 0.1, forces[i], 0.0));
    let q: [f64; 5] = std::array::from_fn(|i| reference[i][0]);
    let old = opening_from_endpoint(&target, &q, dt);
    let mut line = ThinFilmLine::carreau(
        law,
        &temperatures,
        2.0,
        0.04,
        20.0,
        geometry(ThinFilmBoundary::Open, 1.5),
        old.clone(),
        0.0,
    )
    .unwrap();
    assert_eq!(line.profile_source_count(), 5);
    let report = line
        .advance_open(
            dt,
            body / 2.0,
            &[0.0; 3],
            ends(left, right, 0.0),
            limits(),
            controls(),
            &old,
            || false,
        )
        .unwrap();
    for &h in line.heights_m() {
        near(h, 0.1, 3e-11);
    }
    let open = report.open_ends.unwrap();
    let profiles = [
        open.left,
        report.faces[0].profile,
        report.faces[1].profile,
        report.faces[2].profile,
        open.right,
    ];
    for (actual, expected) in profiles.iter().zip(&reference) {
        near(
            actual.flow_per_width_m2_s,
            expected[0],
            2e-8 * expected[0].abs(),
        );
        near(
            actual.surface_velocity_m_s,
            expected[1],
            2e-8 * expected[1].abs(),
        );
        near(
            actual.viscous_dissipation_w_m2,
            expected[2],
            2e-8 * expected[2].abs(),
        );
    }
    let heat = dt
        * 0.5
        * (0.5 * reference[0][2]
            + reference[1][2]
            + reference[2][2]
            + reference[3][2]
            + 0.5 * reference[4][2]);
    near(report.viscous_heat_j_per_m, heat, 1e-11);
    near(
        report.boundary_work_j_per_m,
        dt * (left * q[0] - right * q[4]),
        1e-11,
    );
    near(report.volume_change_m2, dt * (q[0] - q[4]), 1e-11);
    near(report.energy_residual_j_per_m, 0.0, 1e-11);
}

#[test]
fn explicit_open_end_admission_and_late_cancellation_preserve_the_accepted_inventory() {
    let initial = vec![0.10, 0.12, 0.09, 0.11];
    let make = || {
        ThinFilmLine::newtonian(
            material(0.04),
            geometry(ThinFilmBoundary::Open, 1.5),
            initial.clone(),
            0.0,
        )
        .unwrap()
    };
    let ends = ends(0.4, 0.2, 0.03);
    let mut victim = make();
    assert!(
        victim
            .advance(0.4, 0.3, &[0.03; 3], limits(), controls(), &initial, || {
                false
            })
            .is_err()
    );
    let invalid = ThinFilmOpenEnds {
        right: ThinFilmOpenEnd {
            pressure_pa: f64::NAN,
            ..ends.right
        },
        ..ends
    };
    assert!(
        victim
            .advance_open(
                0.4,
                0.3,
                &[0.03; 3],
                invalid,
                limits(),
                controls(),
                &initial,
                || false
            )
            .is_err()
    );
    assert!(
        victim
            .advance_open(
                0.4,
                0.3,
                &[0.03; 3],
                ends,
                limits(),
                ThinFilmControls {
                    maximum_profile_evaluations: 1,
                    ..controls()
                },
                &initial,
                || false
            )
            .is_err()
    );
    let mut closed = ThinFilmLine::newtonian(
        material(0.04),
        geometry(ThinFilmBoundary::Closed, 1.5),
        initial.clone(),
        0.0,
    )
    .unwrap();
    assert!(
        closed
            .advance_open(
                0.4,
                0.3,
                &[0.03; 3],
                ends,
                limits(),
                controls(),
                &initial,
                || false
            )
            .is_err()
    );
    assert!(
        ThinFilmLine::carreau(
            source(0.4),
            &[300.0; 3],
            2.0,
            0.04,
            20.0,
            geometry(ThinFilmBoundary::Open, 1.5),
            initial.clone(),
            0.0
        )
        .is_err()
    );
    assert!(
        ThinFilmLine::carreau(
            source(0.4),
            &[300.0, 300.0, 300.0, 300.0, 361.0],
            2.0,
            0.04,
            20.0,
            geometry(ThinFilmBoundary::Open, 1.5),
            initial.clone(),
            0.0
        )
        .is_err()
    );
    assert_eq!(victim.heights_m(), initial);
    assert_eq!(victim.accepted_time_s(), 0.0);
    let mut reference = make();
    let mut polls = 0;
    let expected = reference
        .advance_open(
            0.4,
            0.3,
            &[0.03; 3],
            ends,
            limits(),
            controls(),
            &initial,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut calls = 0;
    assert!(
        victim
            .advance_open(
                0.4,
                0.3,
                &[0.03; 3],
                ends,
                limits(),
                controls(),
                &initial,
                || {
                    calls += 1;
                    calls == polls
                }
            )
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(victim.heights_m(), initial);
    assert_eq!(victim.accepted_time_s(), 0.0);
    let actual = victim
        .advance_open(
            0.4,
            0.3,
            &[0.03; 3],
            ends,
            limits(),
            controls(),
            &initial,
            || false,
        )
        .unwrap();
    assert_eq!(victim.heights_m(), reference.heights_m());
    assert_eq!(actual.open_ends, expected.open_ends);
    assert_eq!(
        actual.boundary_work_j_per_m.to_bits(),
        expected.boundary_work_j_per_m.to_bits()
    );
    assert_eq!(
        actual.volume_balance_residual_m2.to_bits(),
        expected.volume_balance_residual_m2.to_bits()
    );
    assert_eq!(
        actual.energy_residual_j_per_m.to_bits(),
        expected.energy_residual_j_per_m.to_bits()
    );
}
