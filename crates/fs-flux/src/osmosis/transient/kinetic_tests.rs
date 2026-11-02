use super::*;

fn bath(site: FilmLangmuirSite, ratio: f64, bulk_c: f64, film_c: f64) -> TransientOsmoticPair {
    TransientOsmoticPair::new(
        [bulk(&[bulk_c], 1e-6, 0.0), bulk(&[0.0], 1e-6, 0.0)],
        1e-4,
        membrane(0.0, 1.0, 0.0, 1),
        Some(
            TransientFilm::new(1e-3, 1, vec![1e-8], vec![film_c])
                .unwrap()
                .with_linear_sorption(vec![ratio])
                .unwrap()
                .with_langmuir_sorption(vec![Some(site)])
                .unwrap(),
        ),
        None,
        1 << 20,
    )
    .unwrap()
}
fn phases(pair: &TransientOsmoticPair) -> Vec<u64> {
    let mut result = snapshot(pair);
    let k = pair.kinetic.as_ref().unwrap();
    result.extend(
        k.state
            .equilibrium
            .iter()
            .chain(&k.state.trapped)
            .map(|v| v.to_bits()),
    );
    result
}

#[test]
fn kinetic_finite_bath_matches_implicit_quadratic_and_separates_all_three_pools() {
    let q = 1e-5;
    let ka = 0.02;
    let kd = 0.1;
    let dt = 10.0;
    let mut pair = bath(
        FilmLangmuirSite::new(q, ka, kd, 0.0).unwrap(),
        3.0,
        100.0,
        0.0,
    );
    let before = pair.total_amount(0);
    // Eliminate the finite bulk from the ORIGINAL two mobile balances. The
    // remaining Langmuir equation is quadratic in z=R*V*c, not a fitted oracle.
    let capacity = 3e-7;
    let h = dt * 2e-9; // area*D/(half width)*dt
    let beta = h / (1.0 + h / 1e-6);
    let alpha = 1.0 + beta / capacity;
    let target = beta * 100.0;
    let a = dt * ka / capacity;
    let d = 1.0 + dt * kd;
    let qa = alpha * a;
    let qb = alpha * d + a * q - target * a;
    let qc = -target * d;
    let expected_z = -2.0 * qc / (qb + (qb * qb - 4.0 * qa * qc).sqrt());
    let expected_bound = a * q * expected_z / (d + a * expected_z);
    let report = pair
        .step(dt, TransientStepControls::default(), || false)
        .unwrap();
    let mobile = pair.film_mobile_solute_moles(0, 0, 0).unwrap();
    let bound = pair.film_trapped_solute_moles(0, 0, 0).unwrap();
    close(mobile * 3.0, expected_z, 1e-10, 1e-18);
    close(bound, expected_bound, 1e-10, 1e-18);
    close(
        pair.film_sorbed_solute_moles(0, 0, 0).unwrap(),
        2.0 * mobile + bound,
        1e-14,
        1e-22,
    );
    close(
        pair.film_solute_moles(0, 0).unwrap()[0],
        3.0 * mobile + bound,
        1e-14,
        1e-22,
    );
    close(pair.total_amount(0), before, 1e-12, 1e-18);
    assert!(report.maximum_sorption_residual_mol < 1e-15);
    assert_eq!(report.maximum_linear_residual_mol, 0.0);
}

#[test]
fn finite_sites_saturate_and_preloaded_sites_release_with_analytical_time_refinement() {
    let q = 1e-5;
    let mut saturating = bath(
        FilmLangmuirSite::new(q, 0.1, 0.0, 0.0).unwrap(),
        1.0,
        100.0,
        0.0,
    );
    let before = saturating.total_amount(0);
    for _ in 0..20 {
        saturating
            .step(100.0, TransientStepControls::default(), || false)
            .unwrap();
    }
    let loaded = saturating.film_trapped_solute_moles(0, 0, 0).unwrap();
    assert!(loaded > 0.999 * q && loaded <= q);
    close(saturating.total_amount(0), before, 1e-10, 1e-16);

    let make = || {
        bath(
            FilmLangmuirSite::new(q, 0.0, 0.2, q).unwrap(),
            2.0,
            0.0,
            0.0,
        )
    };
    let duration: f64 = 5.0;
    let expected = q * (-0.2 * duration).exp();
    let mut errors = Vec::new();
    for tolerance in [1e-2, 1e-4] {
        let mut release = make();
        let mut workspace = release.adaptive_workspace(1 << 20).unwrap();
        let report = release
            .advance(
                duration,
                TransientAdvanceControls {
                    maximum_step_s: duration,
                    relative_tolerance: tolerance,
                    absolute_amount_tolerance_mol: 1e-12 * tolerance,
                    absolute_face_concentration_tolerance_mol_m3: 1e6,
                    ..Default::default()
                },
                &mut workspace,
                || false,
            )
            .unwrap();
        let trapped = release.film_trapped_solute_moles(0, 0, 0).unwrap();
        errors.push((trapped - expected).abs());
        assert!(
            report.rejected_steps > 0,
            "internal pool exchange must trigger time refinement"
        );
        assert!(release.film_mobile_solute_moles(0, 0, 0).unwrap() > 0.0);
        close(release.total_amount(0), q, 1e-10, 1e-16);
        assert!(report.solute_balance_residual_mol[0].abs() < 1e-15);
    }
    assert!(errors[1] < 0.2 * errors[0]);
    assert!(errors[1] < 0.006 * q);
}

#[test]
fn different_site_loadings_at_mobile_equilibrium_create_no_osmotic_force() {
    let make_film = |q, ratio| {
        film(1, 25.0)
            .with_linear_sorption(vec![ratio])
            .unwrap()
            .with_langmuir_sorption(vec![Some(
                FilmLangmuirSite::new(q, 0.02, 0.5, 0.5 * q).unwrap(),
            )])
            .unwrap()
    };
    let mut pair = TransientOsmoticPair::new(
        [bulk(&[25.0], 1e-6, 0.0), bulk(&[25.0], 1e-6, 0.0)],
        1e-4,
        membrane(2e-11, 1.0, 1e-6, 1),
        Some(make_film(1e-5, 2.0)),
        Some(make_film(4e-7, 5.0)),
        1 << 20,
    )
    .unwrap();
    let before = phases(&pair);
    let report = pair
        .step(1e4, TransientStepControls::default(), || false)
        .unwrap();
    assert_eq!(report.solvent_transfer_m3, 0.0);
    for side in 0..2 {
        close(
            pair.film_mobile_concentration_mol_m3(side, 0, 0).unwrap(),
            25.0,
            1e-14,
            1e-14,
        );
    }
    assert_eq!(&phases(&pair)[1..], &before[1..]);
    assert_eq!(
        pair.film_langmuir_site(0, 0, 0).unwrap().capacity_mol(),
        1e-5
    );

    // A finite preloaded pool many orders larger than the mobile inventory may
    // round the total to its bound part; the original mobile pool must survive.
    let mut trace = bath(
        FilmLangmuirSite::new(1.0, 0.0, 0.0, 1.0).unwrap(),
        1.0,
        1e-15,
        1e-15,
    );
    let mobile = trace.film_mobile_solute_moles(0, 0, 0).unwrap();
    assert_eq!(trace.film_solute_moles(0, 0).unwrap()[0], 1.0);
    trace
        .step(1.0, TransientStepControls::default(), || false)
        .unwrap();
    assert_eq!(trace.film_mobile_solute_moles(0, 0, 0), Some(mobile));
    assert!(mobile > 0.0);
}

#[test]
fn kinetic_root_work_and_late_adaptive_cancellation_restore_both_pools() {
    let make = || {
        let mut pair = bath(
            FilmLangmuirSite::new(1e-5, 0.02, 0.1, 0.0).unwrap(),
            3.0,
            100.0,
            0.0,
        );
        // Enable actual water/root trials on the same finite storage problem.
        pair.membrane = fs_material::osmosis::polarization::PolarizedMembraneTransport::new(
            membrane(2e-11, 1.0, 0.0, 1),
            vec![Some(
                FilmPermeance::from_diffusion_layer(1e-8, 0.5e-3).unwrap(),
            )],
            vec![None],
        )
        .unwrap();
        pair
    };
    let mut pair = make();
    let before = phases(&pair);
    assert!(
        pair.step(
            10.0,
            TransientStepControls {
                maximum_sorption_iterations: 1,
                ..Default::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(phases(&pair), before);
    assert!(
        pair.step(
            10.0,
            TransientStepControls {
                maximum_root_evaluations: 2,
                ..Default::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(phases(&pair), before);
    let controls = TransientAdvanceControls {
        maximum_step_s: 1.0,
        relative_tolerance: 1e-3,
        ..Default::default()
    };
    let mut workspace = pair.adaptive_workspace(1 << 20).unwrap();
    let mut polls = 0;
    let report = pair
        .advance(1.0, controls, &mut workspace, || {
            polls += 1;
            false
        })
        .unwrap();
    assert_ne!(report.solvent_transfer_m3, 0.0);
    assert_eq!(
        pair.film_langmuir_site(0, 0, 0).unwrap().capacity_mol(),
        1e-5
    );
    let mut rollback = make();
    let mut work = rollback.adaptive_workspace(1 << 20).unwrap();
    let mut seen = 0;
    let result = rollback.advance(1.0, controls, &mut work, || {
        seen += 1;
        seen == polls
    });
    assert!(matches!(
        result,
        Err(TransientOsmosisError::Network(
            OsmoticNetworkError::Cancelled
        ))
    ));
    assert_eq!(phases(&rollback), before);
}

#[test]
fn finite_site_admission_and_inactive_definitions_preserve_legacy_storage() {
    assert!(FilmLangmuirSite::new(0.0, 1.0, 1.0, 0.0).is_err());
    assert!(FilmLangmuirSite::new(1.0, -1.0, 1.0, 0.0).is_err());
    assert!(FilmLangmuirSite::new(1.0, 1.0, 1.0, 2.0).is_err());
    assert!(film(2, 0.0).with_langmuir_sorption(vec![None]).is_err());
    let plain = film(2, 50.0);
    let inactive = plain.clone().with_langmuir_sorption(vec![None; 2]).unwrap();
    let prepare = |left| {
        TransientOsmoticPair::new(
            [bulk(&[50.0], 1e-6, 0.0), bulk(&[100.0], 1e-6, 0.0)],
            1e-4,
            membrane(2e-11, 1.0, 0.0, 1),
            Some(left),
            None,
            1 << 20,
        )
        .unwrap()
    };
    let mut original = prepare(plain);
    let mut equivalent = prepare(inactive);
    assert_eq!(original.workspace_bytes(), equivalent.workspace_bytes());
    original
        .step(1.0, TransientStepControls::default(), || false)
        .unwrap();
    equivalent
        .step(1.0, TransientStepControls::default(), || false)
        .unwrap();
    assert_eq!(snapshot(&original), snapshot(&equivalent));
}

#[test]
fn sub_ulp_bound_uptake_refuses_accuracy_instead_of_erasing_site_kinetics() {
    let mut pair = bath(
        FilmLangmuirSite::new(2e12, 1e-17, 0.0, 1e12).unwrap(),
        1.0,
        0.1,
        0.1,
    );
    let before = phases(&pair);
    // The shared site's broad capacity-scaled roundoff bound admits this
    // representable gross uptake, although the stored bound endpoint cannot
    // resolve it. Coupled storage alone would then falsely see no reaction.
    let step = fs_material::adsorption::competitive::sites::scalar::backward_euler_site(
        1e12, 2e12, 1e-18, 0.0, 1.0,
    )
    .unwrap();
    assert_eq!(step.occupied_moles, 1e12);
    assert!(step.gross_adsorption_mol > 1e-7);
    assert!(step.maximum_balance_residual_mol > 1e-7);
    let result = pair.step(1.0, TransientStepControls::default(), || false);
    assert!(matches!(
        result,
        Err(TransientOsmosisError::SorptionConvergence { .. })
    ));
    assert_eq!(phases(&pair), before);
}
