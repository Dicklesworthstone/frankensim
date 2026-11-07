use super::*;
use fs_solid::dynamics::sustained::events::{
    CollinearCoulombContact, CoulombPhaseMode, SolidCoulombEventStep,
};

fn fixture(
    h: f64,
    vx: f64,
    vs: f64,
    support_a: f64,
    spectator: bool,
    cx: &Cx<'_>,
) -> SupportedSpectralSolid {
    let n = if spectator { 5 } else { 4 };
    let mut k = Coo::new(n, n);
    let mut m = Coo::new(n, n);
    for i in 0..n {
        m.push(i, i, if i == 3 { 2.0 } else { 1.0 });
    }
    // Original consistent cross mass couples free x motion to the actuator.
    m.push(0, 3, 0.25);
    m.push(3, 0, 0.25);
    if spectator {
        // A real elastic mode driven by the accelerating support. It does not
        // participate in the rigid friction trace, but its whole native source
        // and contribution to the actuator reaction must remain present.
        for (i, j, value) in [(3, 3, 4.0), (4, 4, 4.0), (3, 4, -4.0), (4, 3, -4.0)] {
            k.push(i, j, value);
        }
    }
    let assembly = TetElasticAssembly {
        stiffness: k.assemble(),
        mass: m.assemble(),
        free_dofs: (0..n).collect(),
        total_mass_kg: 1.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    };
    let nf = n - 1;
    let mut basis = vec![0.0; nf * nf];
    for i in 0..nf {
        basis[i * nf + i] = 1.0;
    }
    let mut lambda = vec![0.0; nf];
    let (mut q, mut v) = (vec![0.0; n], vec![0.0; n]);
    v[0] = vx;
    v[3] = vs;
    if spectator {
        lambda[3] = 4.0;
        q[4] = 0.1;
        v[4] = 0.2;
    }
    SupportedSpectralSolid::prepare_with_basis(
        &assembly,
        &[3],
        &lambda,
        &basis,
        h,
        &q,
        &v,
        &[support_a],
        native_controls(),
        cx,
    )
    .unwrap()
}

fn weights(n: usize) -> [Vec<f64>; 3] {
    let mut rows: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; n]);
    rows[0][1] = 1.0;
    rows[1][0] = 1.0;
    rows[1][3] = -1.0;
    rows[2][2] = 1.0;
    rows
}
fn law(rows: &[Vec<f64>; 3], h: f64) -> CollinearCoulombContact<'_> {
    CollinearCoulombContact {
        contact: SolidSustainedContact {
            id: 91,
            velocity_weights: [&rows[0], &rows[1], &rows[2]],
            gap_m: 0.0,
            prescribed_surface_velocity_mps: [0.0; 3],
            friction_coefficient: 0.6,
            scales: NaturalMapScales {
                normal_n_per_m: 2.0 / (h * h),
                tangent_n_per_m: 2.0 / (h * h),
            },
        },
        kinetic_coefficient: 0.3,
        tangent_axis: [1.0, 0.0],
    }
}
fn balance<P>(r: &SolidCoulombEventStep<P, SupportedSolidStep>) {
    near(r.energy_balance_error_j, 0.0, 1e-9);
    near(r.native.energy.balance_error_j, 0.0, 1e-9);
    near(r.native.energy.cumulative_balance_error_j, 0.0, 1e-9);
    near(
        r.native.energy.applied_work_j,
        r.external_work_j
            + r.prescribed_surface_work_j
            + r.normal_contact_work_j
            + r.frictional_work_j,
        1e-9,
    );
    assert_eq!(
        r.support_work_j.to_bits(),
        r.native.energy.support_work_j.to_bits()
    );
    assert!(r.frictional_work_j <= 1e-10);
}

#[test]
fn supported_stop_recaptures_while_actuator_supplies_original_cross_mass_work() {
    with_cx(|cx, _| {
        for boost in [0.0, 5.0] {
            let mut s = fixture(1.0, 2.0 + boost, 1.0 + boost, 0.5, false, cx);
            let rows = weights(4);
            let r = s
                .advance_collinear_coulomb_events(
                    &[0.5],
                    &law(&rows, 1.0),
                    [0.0; 3],
                    &[0.0, -20.0, 0.0, 0.0],
                    controls(),
                    cx,
                    |path| {
                        assert_eq!(path.before.time_s(), 0.0);
                        assert_eq!(path.before.support_acceleration_mps2(), [0.5]);
                        assert_eq!(path.phases.len(), 2);
                        Ok::<_, ()>("fixed material point and driven support")
                    },
                )
                .unwrap();
            let stop = 8.0 / 53.0;
            let travel = |t: f64| (1.0 + boost) * t + 0.25 * t * t;
            let work = -6.53125 * travel(stop) + 1.75 * (travel(1.0) - travel(stop));
            near(r.stopping_time_offset_s.unwrap(), stop, 1e-14);
            near(r.phases[0].force_n[1], -6.0, 1e-8);
            near(r.phases[1].force_n[1], 0.625, 1e-8);
            assert_eq!(r.phases[1].mode, CoulombPhaseMode::Sticking);
            near(s.velocity_mps()[0], 1.5 + boost, 1e-9);
            near(s.velocity_mps()[3], 1.5 + boost, 1e-12);
            near(s.displacement_m()[0], travel(1.0) + 0.5 * stop, 1e-10);
            near(s.displacement_m()[3], travel(1.0), 1e-12);
            near(r.native.support_reaction_n[0], 1.75, 1e-9);
            near(r.support_work_j, work, 1e-9);
            near(r.frictional_dissipation_j, 3.0 * stop, 1e-9);
            near(r.momentum_change_kg_m_s[0], -0.375, 1e-9);
            near(r.momentum_change_kg_m_s[3], 0.875, 1e-9);
            near(
                r.native.energy.kinetic_after_j - r.native.energy.kinetic_before_j,
                0.4375 + 0.5 * boost,
                1e-9,
            );
            assert!(r.support_work_j > r.frictional_dissipation_j);
            assert_eq!(r.native.step, 1);
            assert_eq!(r.native.time_s, 1.0);
            balance(&r);
        }
    });
}

#[test]
fn supported_reversal_integrates_two_reactions_and_gross_friction_loss() {
    with_cx(|cx, _| {
        let h = 0.5;
        let mut s = fixture(h, 2.0, 1.0, 0.5, false, cx);
        let rows = weights(4);
        let r = s
            .advance_collinear_coulomb_events(
                &[0.5],
                &law(&rows, h),
                [0.0; 3],
                &[-20.0, -20.0, 0.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        let stop = 1.0 / 26.625;
        let remaining = h - stop;
        let forward = 0.5 * stop;
        let backward = 0.5 * 14.625 * remaining * remaining;
        let support_travel = |t: f64| t + 0.25 * t * t;
        let work =
            -11.53125 * support_travel(stop) + 3.46875 * (support_travel(h) - support_travel(stop));
        near(r.stopping_time_offset_s.unwrap(), stop, 1e-14);
        assert_eq!(r.phases.len(), 2);
        assert!(
            r.phases
                .iter()
                .all(|phase| phase.mode == CoulombPhaseMode::Sliding)
        );
        near(r.phases[0].force_n[1], -6.0, 1e-8);
        near(r.phases[1].force_n[1], 6.0, 1e-8);
        near(
            s.velocity_mps()[0],
            1.0 + 0.5 * h - 14.625 * remaining,
            1e-9,
        );
        near(
            s.displacement_m()[0],
            support_travel(h) + forward - backward,
            1e-10,
        );
        near(r.sliding_distance_m, forward + backward, 1e-10);
        near(r.frictional_dissipation_j, 6.0 * (forward + backward), 1e-9);
        near(r.native.support_reaction_n[0], 3.46875, 1e-9);
        near(r.support_work_j, work, 1e-9);
        near(r.external_work_j, -20.0 * s.displacement_m()[0], 1e-9);
        assert!(
            r.response_evaluations > r.phases.iter().map(|p| p.solve.response_evaluations).sum()
        );
        balance(&r);
    });
}

#[test]
fn supported_static_capture_transmits_acceleration_without_frictional_loss() {
    with_cx(|cx, _| {
        let mut s = fixture(1.0, 1.0, 1.0, 0.5, false, cx);
        let rows = weights(4);
        let r = s
            .advance_collinear_coulomb_events(
                &[0.5],
                &law(&rows, 1.0),
                [20.0, 0.625, 0.0],
                &[0.0, -20.0, 0.0, 0.0],
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        assert_eq!(r.phases.len(), 1);
        assert_eq!(r.phases[0].mode, CoulombPhaseMode::Sticking);
        near(r.phases[0].force_n[1], 0.625, 1e-9);
        near(s.displacement_m()[0], 1.25, 1e-10);
        near(s.velocity_mps()[0], 1.5, 1e-9);
        near(r.support_work_j, 2.1875, 1e-9);
        near(r.frictional_dissipation_j, 0.0, 1e-10);
        balance(&r);
    });
}

#[test]
fn supported_event_preserves_driven_elastic_mode_and_native_modal_continuation() {
    with_cx(|cx, _| {
        let (mut s, mut free) = (
            fixture(1.0, 2.0, 1.0, 0.5, true, cx),
            fixture(1.0, 2.0, 1.0, 0.5, true, cx),
        );
        let rows = weights(5);
        let load = [0.0, -20.0, 0.0, 0.0, 0.0];
        free.advance_affine_supports(&[0.5], &load, &load, cx)
            .unwrap();
        let r = s
            .advance_collinear_coulomb_events(
                &[0.5],
                &law(&rows, 1.0),
                [0.0; 3],
                &load,
                controls(),
                cx,
                |_| Ok::<_, ()>(()),
            )
            .unwrap();
        assert_eq!(
            s.displacement_m()[4].to_bits(),
            free.displacement_m()[4].to_bits()
        );
        assert_eq!(
            s.velocity_mps()[4].to_bits(),
            free.velocity_mps()[4].to_bits()
        );
        assert_eq!(
            s.stored_energy_j().to_bits(),
            free.stored_energy_j().to_bits()
        );
        balance(&r);
        let held = [0.625, 0.0, 0.0, -0.625, 0.0];
        let next = s.advance_affine_supports(&[0.5], &held, &held, cx).unwrap();
        free.advance_affine_supports(&[0.5], &load, &load, cx)
            .unwrap();
        near(s.displacement_m()[0], 3.0 + 4.0 / 53.0, 1e-10);
        near(s.velocity_mps()[0], 2.0, 1e-9);
        assert_eq!(
            s.displacement_m()[4].to_bits(),
            free.displacement_m()[4].to_bits()
        );
        assert_eq!(
            s.velocity_mps()[4].to_bits(),
            free.velocity_mps()[4].to_bits()
        );
        near(next.energy.cumulative_balance_error_j, 0.0, 1e-9);
    });
}

#[test]
fn supported_event_acceleration_budget_path_and_cancel_refusals_preserve_both_owners() {
    with_cx(|cx, gate| {
        let mut s = fixture(1.0, 2.0, 1.0, 0.5, false, cx);
        let rows = weights(4);
        let law = law(&rows, 1.0);
        let load = [0.0, -20.0, 0.0, 0.0];
        let before = bits(&s);
        let r = s.advance_collinear_coulomb_events(
            &[0.75],
            &law,
            [0.0; 3],
            &load,
            controls(),
            cx,
            |_| -> Result<(), ()> { panic!("changing support acceleration") },
        );
        assert!(matches!(
            r,
            Err(SolidSustainedError::Native(SupportedSolidError::Invalid(_)))
        ));
        assert_eq!(bits(&s), before);
        let mut elastic = fixture(1.0, 2.0, 1.0, 0.5, true, cx);
        let mut nonrigid_rows = weights(5);
        nonrigid_rows[1][4] = 0.125;
        let elastic_before = bits(&elastic);
        let result = elastic.advance_collinear_coulomb_events(
            &[0.5],
            &self::law(&nonrigid_rows, 1.0),
            [0.0; 3],
            &[0.0, -20.0, 0.0, 0.0, 0.0],
            controls(),
            cx,
            |_| -> Result<(), ()> { panic!("nonrigid supported contact") },
        );
        assert!(matches!(result, Err(SolidSustainedError::Native(_))));
        assert_eq!(bits(&elastic), elastic_before);
        let mut reference = fixture(1.0, 2.0, 1.0, 0.5, false, cx);
        let accepted = reference
            .advance_collinear_coulomb_events(&[0.5], &law, [0.0; 3], &load, controls(), cx, |_| {
                Ok::<_, ()>(())
            })
            .unwrap();
        let mut limited = controls();
        limited.solve.maximum_response_evaluations =
            accepted.phases[0].solve.response_evaluations + 1;
        let r = s.advance_collinear_coulomb_events(
            &[0.5],
            &law,
            [0.0; 3],
            &load,
            limited,
            cx,
            |_| -> Result<(), ()> { panic!("event response budget") },
        );
        assert!(matches!(r, Err(SolidSustainedError::Budget { .. })));
        assert_eq!(bits(&s), before);
        let r = s.advance_collinear_coulomb_events(
            &[0.5],
            &law,
            [0.0; 3],
            &load,
            controls(),
            cx,
            |_| Err::<(), _>("changed feature"),
        );
        assert!(matches!(
            r,
            Err(SolidSustainedError::Path("changed feature"))
        ));
        assert_eq!(bits(&s), before);
        let r = s.advance_collinear_coulomb_events(
            &[0.5],
            &law,
            [0.0; 3],
            &load,
            controls(),
            cx,
            |_| {
                gate.request();
                Ok::<_, ()>(())
            },
        );
        assert!(matches!(r, Err(SolidSustainedError::Cancelled)));
        assert_eq!(bits(&s), before);
        with_cx(|fresh, _| {
            // A normal native retry observes the actual unmodified modal
            // predecessor after all refused event phases and support trials.
            let mut original = fixture(1.0, 2.0, 1.0, 0.5, false, fresh);
            s.advance_affine_supports(&[0.5], &load, &load, fresh)
                .unwrap();
            original
                .advance_affine_supports(&[0.5], &load, &load, fresh)
                .unwrap();
            assert_eq!(bits(&s), bits(&original));
        });
    });
}
