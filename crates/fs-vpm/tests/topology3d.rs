//! G0/G2/G3: evaluate the physical edges of the wake, not a tautological
//! circulation formula. Independent closed rings and nodal incidence are
//! the oracles for sign, shared-edge cancellation, and near/far conversion.
use fs_vpm::coarsen3d::{coarsen_oldest, station_impulse};
use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode, cell_from_segs};
use fs_vpm::filament3d::{
    FilamentWake, WakeRateCertificate, WakeSegment, segment_velocity_with_core,
};
use fs_vpm::images3d::{ImagePlane, ImageSystem, grounded_hybrid_velocity};
use std::collections::BTreeMap;

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

fn near(a: [f64; 3], b: [f64; 3]) {
    let error = norm(std::array::from_fn(|i| a[i] - b[i]));
    assert!(error < 3e-12 * (1.0 + norm(b)), "{a:?} != {b:?}");
}

fn wake(rows: usize, unsteady: bool) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 100.0,
            n_stations: 3,
            max_rows: 64,
        },
        (0..4).map(|i| [0.0, f64::from(i) - 1.5, 1.0]).collect(),
    )
    .unwrap();
    for i in 0..rows {
        let t = i as f64;
        let gamma = if unsteady {
            [1.0 + 0.125 * t, 0.5 - 0.0625 * t, -0.5 + 0.125 * t]
        } else {
            [1.0; 3]
        };
        wake.shed(&gamma, [-0.125, 0.0, 0.0]).unwrap();
    }
    for row in &mut wake.rows {
        row.core2_m2 = 0.04;
    }
    wake
}

fn ring(p: [f64; 3], corners: [[f64; 3]; 4], gamma: f64, core2: f64) -> [f64; 3] {
    let mut velocity = [0.0; 3];
    for edge in 0..4 {
        let value = segment_velocity_with_core(p, corners[edge], corners[(edge + 1) % 4], core2);
        for i in 0..3 {
            velocity[i] += gamma * value[i];
        }
    }
    velocity
}

fn segments(wake: &FilamentWake, growth: f64) -> Vec<WakeSegment> {
    let mut out = Vec::new();
    wake.visit_segments_with_core_growth(growth, |s| out.push(s));
    out
}

/// Assemble incidence independently from the actual emitted source endpoints.
fn assert_closed(segments: &[WakeSegment]) {
    let mut balance = BTreeMap::<[u64; 3], (f64, f64)>::new();
    for s in segments {
        for (point, sign) in [(s.start, -1.0), (s.end, 1.0)] {
            let (net, scale) = balance.entry(point.map(f64::to_bits)).or_default();
            *net += sign * s.circulation;
            *scale += s.circulation.abs();
        }
    }
    for (point, (net, scale)) in balance {
        assert!(
            net.abs() <= 8.0 * f64::EPSILON * scale,
            "unbalanced vortex endpoint {point:?}: {net}"
        );
    }
}

fn impulse(segments: &[WakeSegment]) -> [f64; 3] {
    let mut value = [0.0; 3];
    for s in segments {
        let (a, b) = (s.start, s.end);
        let cross = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        for i in 0..3 {
            value[i] += 0.5 * s.circulation * cross[i];
        }
    }
    value
}

#[test]
fn g2_steady_wake_equals_one_outer_ring_at_every_shedding_resolution() {
    let length = 1.0;
    for rows in [3, 9, 17] {
        let mut w = wake(rows, false);
        for (r, row) in w.rows.iter_mut().enumerate() {
            let x = -length * (rows - r - 1) as f64 / (rows - 1) as f64;
            for point in &mut row.nodes {
                point[0] = x;
            }
        }
        let outer = [
            w.line_nodes[0],
            w.line_nodes[3],
            w.rows[0].nodes[3],
            w.rows[0].nodes[0],
        ];
        for p in [[-0.3, 0.0, 1.2], [0.4, 0.2, 1.1], [-0.8, 1.7, 0.9]] {
            near(w.induced_velocity(p), ring(p, outer, 1.0, 0.04));
        }
        assert_closed(&segments(&w, 0.0));
    }
}

#[test]
fn g2_unsteady_sheet_matches_independent_oriented_panel_rings() {
    let mut w = wake(9, true);
    for row in &mut w.rows {
        for (i, point) in row.nodes.iter_mut().enumerate() {
            point[0] += 0.02 * (i * i) as f64;
            point[2] += 0.03 * i as f64;
        }
    }
    for (i, point) in w.line_nodes.iter_mut().enumerate() {
        point[0] += 0.02 * (i * i) as f64;
        point[2] += 0.03 * i as f64;
    }
    for p in [[-0.3, 0.2, 1.4], [0.4, -0.7, 0.8], [-2.0, 0.5, 1.0]] {
        let mut reference = [0.0; 3];
        for (r, row) in w.rows.iter().enumerate() {
            let newer = w.rows.get(r + 1).map_or(&w.line_nodes, |r| &r.nodes);
            for s in 0..3 {
                let v = ring(
                    p,
                    [newer[s], newer[s + 1], row.nodes[s + 1], row.nodes[s]],
                    row.gamma[s],
                    0.04,
                );
                for i in 0..3 {
                    reference[i] += v[i];
                }
            }
        }
        near(w.induced_velocity(p), reference);
    }
    assert_closed(&segments(&w, 0.0));
    for row in 0..=w.rows.len() {
        for edge in 0..=3 {
            assert_eq!(w.node_circulation_residual(row, edge), Some(0.0));
        }
    }
    assert_eq!(w.node_circulation_residual(w.rows.len() + 1, 0), None);
    assert_eq!(w.node_circulation_residual(0, 4), None);
}

#[test]
fn g0_steady_interior_vorticity_cancels_before_distinct_cores_are_evaluated() {
    let mut w = wake(12, false);
    for (r, row) in w.rows.iter_mut().enumerate() {
        row.core2_m2 = 0.001 * r as f64;
    }
    let sources = segments(&w, 0.02);
    let spans: Vec<_> = sources.iter().filter(|s| s.start[1] != s.end[1]).collect();
    assert_eq!(spans.len(), 6, "only the starting and closing spans remain");
    for s in spans {
        assert!(s.start[0] == w.rows[0].nodes[0][0] || s.start[0] == 0.0);
    }
    assert_closed(&sources);
}

#[test]
fn g0_physical_source_impulse_matches_the_coarsening_invariant() {
    let mut w = wake(16, true);
    let before = impulse(&segments(&w, 0.0));
    let expected_z: f64 = (0..3).map(|s| station_impulse(&w, s)).sum();
    near(before, [0.0, 0.0, expected_z]);
    coarsen_oldest(&mut w, 4).unwrap();
    let after = segments(&w, 0.0);
    assert_closed(&after);
    near(impulse(&after), before);
}

#[test]
fn g3_repeated_near_far_partition_preserves_physical_edges_and_the_viscous_field() {
    let full = wake(12, true);
    let growth = 0.02;
    let mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 };
    let mut near_wake = full.clone();
    let mut far = FarField::aggregate(&mut near_wake, 4, 2, mode, 0.01).unwrap();
    let second = FarField::aggregate(&mut near_wake, 2, 1, mode, 0.01).unwrap();
    far.cells.extend(second.cells);
    assert_eq!(
        near_wake.downstream_gamma.as_ref(),
        Some(&full.rows[5].gamma)
    );
    let mut combined = segments(&near_wake, growth);
    for cell in &far.cells {
        for ((start, end, circulation), core2_m2) in cell.segs.iter().zip(&cell.source_core2_m2) {
            combined.push(WakeSegment {
                start: *start,
                end: *end,
                circulation: *circulation,
                core2_m2: *core2_m2,
            });
        }
    }
    assert_closed(&combined);
    for p in [[-0.3, 0.2, 1.2], [0.4, -0.7, 0.8], [-2.0, 0.5, 1.0]] {
        let expected = full.induced_velocity_with_core_growth(p, growth);
        let mut exact = near_wake.induced_velocity_with_core_growth(p, growth);
        for cell in &far.cells {
            let v = cell.exact(p);
            for i in 0..3 {
                exact[i] += v[i];
            }
        }
        near(exact, expected);
    }
    let images = ImageSystem::build(&near_wake, &far, ImagePlane { z0_m: 0.0 }).unwrap();
    for p in [[-0.3, 0.2, 0.0], [0.4, -0.7, 0.0], [-2.0, 0.5, 0.0]] {
        let v = grounded_hybrid_velocity(&near_wake, &far, &images, p);
        assert!(v[2].abs() < 1e-12 * (1.0 + norm(v)));
    }
}

#[test]
fn g0_attached_interface_cannot_be_coarsened_away() {
    let mut w = wake(12, true);
    let _far = FarField::aggregate(&mut w, 4, 2, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
    let before = w.clone();
    assert_eq!(
        coarsen_oldest(&mut w, 2).unwrap_err().code,
        "coarsen-interface-attached"
    );
    assert_eq!(w, before);
}

#[test]
fn g0_identity_tracks_closing_line_and_partition_circulation() {
    let w = wake(8, true);
    let baseline = w.digest();
    let mut displaced = w.clone();
    displaced.line_nodes[0][0] += 0.05;
    assert_ne!(baseline, displaced.digest());
    let mut connected = w.clone();
    connected.downstream_gamma = Some(vec![0.25, 0.5, 0.75]);
    assert_ne!(baseline, connected.digest());
    let p = [-0.9, 0.0, 1.1];
    assert_ne!(w.induced_velocity(p), connected.induced_velocity(p));
}

#[test]
fn g2_wake_jacobian_differentiates_the_actual_mixed_core_field() {
    let mut w = wake(9, true);
    for (i, row) in w.rows.iter_mut().enumerate() {
        row.core2_m2 = 0.01 * (i + 1) as f64;
    }
    let growth = 0.02;
    let p = [-0.5, 0.3, 1.2];
    let k = w.induced_kinematics_with_core_growth(p, growth).unwrap();
    near(k.velocity, w.induced_velocity_with_core_growth(p, growth));
    let h = 2e-5;
    for axis in 0..3 {
        let mut left = p;
        let mut right = p;
        left[axis] -= h;
        right[axis] += h;
        let vl = w.induced_velocity_with_core_growth(left, growth);
        let vr = w.induced_velocity_with_core_growth(right, growth);
        for component in 0..3 {
            let fd = (vr[component] - vl[component]) / (2.0 * h);
            assert!((k.gradient[component][axis] - fd).abs() < 2e-8 * (1.0 + fd.abs()));
        }
    }
    assert!(k.divergence().abs() < 1e-13);
    assert!(
        norm(k.vorticity()) > 0.01,
        "finite-core vorticity is resolved"
    );
}

#[test]
fn g2_closed_inviscid_sheet_is_irrotational_off_filaments_and_refuses_the_axis() {
    let mut w = wake(9, true);
    for row in &mut w.rows {
        row.core2_m2 = 0.0;
    }
    for p in [[-0.5, 0.3, 1.2], [-1.2, -0.8, 0.9], [2.0, 3.0, 4.0]] {
        let k = w.induced_kinematics_with_core_growth(p, 0.0).unwrap();
        near(k.velocity, w.induced_velocity(p));
        assert!(
            norm(k.vorticity()) < 2e-12,
            "closed-sheet curl {:?}",
            k.vorticity()
        );
        assert!(k.divergence().abs() < 1e-13);
    }
    let on_start = [w.rows[0].nodes[0][0], 0.0, 1.0];
    assert!(
        w.induced_kinematics_with_core_growth(on_start, 0.0)
            .is_none()
    );
    for row in &mut w.rows {
        row.core2_m2 = 0.04;
    }
    let regularized = w
        .induced_kinematics_with_core_growth(on_start, 0.0)
        .unwrap();
    assert!(norm(regularized.vorticity()) > 0.1);
    assert!(
        w.induced_kinematics_with_core_growth([f64::NAN; 3], 0.0)
            .is_none()
    );
    assert!(
        w.induced_kinematics_with_core_growth([0.0; 3], -0.1)
            .is_none()
    );
}

#[test]
fn g3_uniform_translation_preserves_exact_and_multipole_velocity() {
    let mut w = wake(12, true);
    let mut far = FarField::aggregate(&mut w, 4, 2, WakeCoreEvolutionMode::Frozen, 0.01).unwrap();
    let before = far.clone();
    far.translate([0.0; 3]).unwrap();
    assert_eq!(before.digest(), far.digest());
    let delta = [1.25, -0.5, 0.25];
    far.translate(delta).unwrap();
    for (old, translated) in before.cells.iter().zip(&far.cells) {
        assert_eq!(old.w_sum, translated.w_sum);
        assert_eq!(old.dip, translated.dip);
        assert_eq!(old.quad, translated.quad);
        assert_eq!(old.core_moment, translated.core_moment);
        assert_eq!(old.source_core2_m2, translated.source_core2_m2);
        for p in [[-0.5, 0.1, 1.2], [12.0, 9.0, 8.0]] {
            let shifted = std::array::from_fn(|i| p[i] + delta[i]);
            near(translated.exact(shifted), old.exact(p));
            near(translated.eval(shifted), old.eval(p));
        }
    }
    assert_ne!(before.digest(), far.digest());
}

#[test]
fn g3_prescribed_shedding_keeps_near_far_connected_and_matches_unsplit_evolution() {
    let mut full = wake(12, true);
    let mut near_wake = full.clone();
    let growth = 0.02;
    let mut far = FarField::aggregate(
        &mut near_wake,
        4,
        2,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 },
        0.01,
    )
    .unwrap();
    for tick in 0..5 {
        let gamma = [1.0 + 0.5 * f64::from(tick), -0.75, 0.25];
        let displacement = [-0.125, 0.015625, 0.03125];
        full.shed(&gamma, displacement).unwrap();
        far.shed_prescribed(&mut near_wake, &gamma, displacement)
            .unwrap();
        let mut physical = segments(&near_wake, growth);
        for cell in &far.cells {
            for ((a, b, gamma), core2) in cell.segs.iter().zip(&cell.source_core2_m2) {
                physical.push(WakeSegment {
                    start: *a,
                    end: *b,
                    circulation: *gamma,
                    core2_m2: *core2 + cell.diffusion_offset_m2,
                });
            }
        }
        assert_closed(&physical);
        for p in [[-0.3, 0.2, 1.2], [0.4, -0.7, 0.8], [-2.0, 0.5, 1.0]] {
            let mut exact = near_wake.induced_velocity_with_core_growth(p, growth);
            for cell in &far.cells {
                let v = cell.exact(p);
                for i in 0..3 {
                    exact[i] += v[i];
                }
            }
            near(exact, full.induced_velocity_with_core_growth(p, growth));
        }
        assert_eq!(near_wake.rows.last().unwrap().age_ticks, 0.0);
    }
}

#[test]
fn g0_translation_and_coupled_shedding_refusals_leave_both_states_unchanged() {
    let mut w = wake(12, true);
    let mut far = FarField::aggregate(
        &mut w,
        4,
        2,
        WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.5 },
        0.01,
    )
    .unwrap();
    w.cert.max_rows = w.rows.len();
    let before_w = w.clone();
    let before_far = far.clone();
    assert_eq!(
        far.shed_prescribed(&mut w, &[1.0; 3], [-0.125, 0.0, 0.0])
            .unwrap_err()
            .code,
        "filament-rows-exhausted"
    );
    assert_eq!(w, before_w);
    assert_eq!(far, before_far);
    w.cert.max_rows += 1;
    let before_w = w.clone();
    assert!(
        far.shed_prescribed(&mut w, &[1.0; 3], [f64::NAN; 3])
            .is_err()
    );
    assert_eq!(w, before_w);
    assert_eq!(far, before_far);
    // One cell could translate, while a later cell would overflow: no partial
    // far update, and the near wake's age/geometry must not advance either.
    let huge = 0.75 * f64::MAX;
    far.cells.push(cell_from_segs(
        vec![([huge, -0.5, 1.0], [huge, 0.5, 1.0], 1.0)],
        0.04,
    ));
    let before_far = far.clone();
    assert_eq!(
        far.translate([huge, 0.0, 0.0]).unwrap_err().code,
        "farfield-translation-invalid"
    );
    assert_eq!(far, before_far);
    assert!(
        far.shed_prescribed(&mut w, &[1.0; 3], [huge, 0.0, 0.0])
            .is_err()
    );
    assert_eq!(w, before_w);
    assert_eq!(far, before_far);
    // Near geometry has its own complete preflight, even without a far field.
    w.rows.last_mut().unwrap().nodes[0][0] = huge;
    let before_w = w.clone();
    assert_eq!(
        w.shed(&[1.0; 3], [huge, 0.0, 0.0]).unwrap_err().code,
        "filament-convection-invalid"
    );
    assert_eq!(w, before_w);
}
