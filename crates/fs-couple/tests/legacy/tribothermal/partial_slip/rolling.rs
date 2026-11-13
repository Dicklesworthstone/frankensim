//! Same original rolling and loaded-contact laws in one finite-temperature solve.
//! All coefficients are authored numerical references, not measured calibration.
use super::*;
use fs_couple::tribothermal::partial_slip::rolling::*;
use fs_tribo::rolling_loss::*;
use fs_tribo::{ApplicabilityRange, InterfaceMedium, InterfaceSystemRef};

fn rolling_law(kind: usize, upper_k: f64) -> RollingLossLaw {
    let domain = RollingLossApplicability::new(
        ApplicabilityRange::new(200.0, upper_k).unwrap(),
        ApplicabilityRange::new(0.0, 100.0).unwrap(),
    )
    .unwrap();
    match kind {
        0 => RollingLossLaw::CoulombContour(
            CoulombContourCard::new(
                "declared contour",
                "synthetic rolling card",
                InputAuthority::SyntheticFixture,
                0.2,
                domain,
            )
            .unwrap(),
        ),
        1 => RollingLossLaw::ViscousContour(
            ViscousContourCard::new(
                "declared viscous contour",
                "synthetic rolling card",
                InputAuthority::SyntheticFixture,
                2.0,
                domain,
            )
            .unwrap(),
        ),
        _ => RollingLossLaw::HystereticRollingMoment(
            HystereticRollingCard::new(
                "declared rolling moment",
                "synthetic rolling card",
                InputAuthority::SyntheticFixture,
                0.002,
                0.25,
                domain,
            )
            .unwrap(),
        ),
    }
}
fn rolling<'a>(
    partial: CaloricPartialSlip<'a>,
    kind: usize,
    upper_k: f64,
    coordinate: &str,
    temperature: RollingTemperatureSurface,
) -> CaloricRollingContact<'a> {
    let p = partial.checkpoint().normal_patch();
    let i = partial.checkpoint().interface();
    let patch = RollingPatchReceipt::new(
        p.patch_id(),
        p.normal_model_id(),
        p.source_id(),
        InputAuthority::SyntheticFixture,
        p.normal_load_n(),
        core::f64::consts::PI * p.semi_axis_longitudinal_m() * p.semi_axis_lateral_m(),
        PatchCurvature::Principal {
            first_per_m: 10.0,
            second_per_m: 10.0,
        },
    )
    .unwrap();
    let interface = InterfaceSystemRef::new(
        i.ordered_interface_id(),
        i.history_id(),
        i.source_id(),
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .unwrap();
    let owner = RollingWorkOwnership::new(
        p.patch_id(),
        "original rolling interval",
        coordinate,
        if kind < 2 {
            RollingLossChannel::ContourDeformation
        } else {
            RollingLossChannel::RollingHysteresis
        },
    )
    .unwrap();
    CaloricRollingContact::new(
        partial,
        patch,
        interface,
        rolling_law(kind, upper_k),
        RollingLossState::zero(),
        None,
        owner,
        temperature,
        HeatPartition::new(0.6, 0.4, 0.0).unwrap(),
    )
    .unwrap()
}
fn rolling_motion() -> RollingMotion {
    RollingMotion {
        contour_speed_mps: 0.0,
        rolling_rate_rad_s: 40.0,
        excitation_frequency_hz: 10.0,
    }
}

#[test]
fn rolling_heats_loaded_contact_without_fabricating_elastic_heat() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let partial = contact(&thermal, &curves, model(1000.0, slopes()), [0.1, 0.0, 0.0]);
    let initial_u = partial.thermodynamics().internal_energy_j;
    let initial_psi = partial.thermodynamics().free_energy_j;
    let old = rolling(
        partial,
        2,
        600.0,
        "rolling",
        RollingTemperatureSurface::SurfaceA,
    );
    let next = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.0),
            &work("one physical interval"),
            rolling_motion(),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // N*l*eta*|omega|*dt = 10*.002*.25*40 = .2 J, split .6/.4.
    near(next.report().rolling().dissipation.total_heat_j, 0.2);
    near(next.report().total_mechanical_work_j, 0.2);
    near(next.report().external_heat_j, 0.0);
    near(next.contact().partial_slip().temperature_k()[0], 300.012);
    near(next.contact().partial_slip().temperature_k()[1], 300.004);
    near(next.report().partial_slip.source.temperature_k, 300.012);
    near(
        next.contact()
            .partial_slip()
            .thermodynamics()
            .internal_energy_j,
        initial_u,
    );
    assert!(next.contact().partial_slip().thermodynamics().free_energy_j < initial_psi);
    assert_eq!(next.report().partial_slip.physical_heat_j, [0.0; 2]);
    assert_eq!(next.report().partial_slip.reversible_heat_j, [0.0; 2]);
    assert_eq!(
        next.report()
            .rolling()
            .generalized_work
            .ownership
            .interval_id(),
        "one physical interval"
    );
    assert_eq!(next.contact().rolling_state().accepted_steps(), 1);
    assert_eq!(old.rolling_state(), &RollingLossState::zero());
    assert!(old.rolling_checkpoint().is_none());

    // A true continuation preserves native history, including reversal.
    let mut reversed = rolling_motion();
    reversed.rolling_rate_rad_s = -40.0;
    let (thermal, old, _) = next.into_parts();
    let next = old
        .step_candidate(
            &thermal,
            2.0,
            frame(),
            motion(0.0),
            &work("second interval"),
            reversed,
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    near(next.contact().rolling_state().cumulative_heat_j(), 0.4);
    assert_eq!(next.contact().rolling_state().accepted_steps(), 2);
    near(next.report().rolling().wrench.rolling_moment_nm, 0.005);
    assert!(next.report().total_energy_residual_j.abs() < 1e-8);
}

#[test]
fn all_original_loss_channels_match_independent_unequal_temperature_conduction() {
    let curves = curves();
    for (kind, expected_heat) in [(0, 6.0), (1, 18.0), (2, 0.02)] {
        let thermal = initial(&curves, [320.0, 280.0], 1.0);
        let partial = contact(&thermal, &curves, model(1000.0, zero_slopes()), [0.0; 3]);
        let old = rolling(
            partial,
            kind,
            600.0,
            "separate rolling",
            RollingTemperatureSurface::SurfaceB,
        );
        let next = old
            .step_candidate(
                &thermal,
                1.0,
                frame(),
                motion(0.0),
                &work("conduction"),
                RollingMotion {
                    contour_speed_mps: -3.0,
                    rolling_rate_rad_s: -4.0,
                    excitation_frequency_hz: 10.0,
                },
                controls(),
                &[0.0; 2],
                |_| {
                    Ok::<_, Infallible>(StepCoefficients {
                        sources_w: BTreeMap::from([(DESTINATIONS[1], 2.0)]),
                        conductances_w_k: BTreeMap::new(),
                        payload: "heater",
                    })
                },
                || false,
            )
            .unwrap();
        // Independently invert [Ca+G,-G;-G,Cb+G] at dt=1.
        let a = 10.0 * 320.0 + 0.6 * expected_heat;
        let b = 20.0 * 280.0 + 0.4 * expected_heat + 2.0;
        let expected = [(21.0 * a + b) / 230.0, (a + 11.0 * b) / 230.0];
        for i in 0..2 {
            near(
                next.contact().partial_slip().temperature_k()[i],
                expected[i],
            );
        }
        near(next.report().partial_slip.source.temperature_k, expected[1]);
        near(
            next.report().rolling().dissipation.total_heat_j,
            expected_heat,
        );
        near(next.report().external_heat_j, 2.0);
        near(next.report().total_mechanical_work_j, expected_heat);
        near(
            next.report().partial_slip.thermal.stored_energy_change_j,
            expected_heat + 2.0,
        );
        assert_eq!(*next.report().source(), "heater");
        assert!(next.report().total_energy_residual_j.abs() < 1e-8);
    }
}

#[test]
fn rolling_heat_changes_pinned_friction_at_the_same_solved_temperature() {
    let (a, ca) = common::material(1.0, 0.05, 0.0);
    let (b, cb) = common::material(1.0, 0.1, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = source(&cards, [(200.0, 0.1), (600.0, 0.5)], 600.0);
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let partial = loaded(
        &thermal,
        &curves,
        data(&cards, &curves, &card, &point),
        0.2,
        [[0.0; 2]; 2],
        1.0,
    );
    let old = rolling(
        partial,
        2,
        600.0,
        "rolling",
        RollingTemperatureSurface::SurfaceA,
    );
    let next = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(1.0),
            &work("loaded slip and roll"),
            rolling_motion(),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // D=10*mu-.292*mu^2 and Ta=300+5*D+.6*.2/.05.
    // Hence mu=.2024+.005*D, with NO duplicated rolling heat credit.
    let slip_heat = root(0.0, 4.0, |d| {
        let mu = 0.2024 + 0.005 * d;
        d - 10.0 * mu + 0.292 * mu * mu
    });
    let mu = 0.2024 + 0.005 * slip_heat;
    near(
        next.report().partial_slip.contact.dissipation.heat_j,
        slip_heat,
    );
    near(next.report().rolling().dissipation.total_heat_j, 0.2);
    near(
        next.contact().partial_slip().temperature_k()[0],
        300.0 + 5.0 * slip_heat + 2.4,
    );
    near(
        next.contact().partial_slip().temperature_k()[1],
        300.0 + 7.5 * slip_heat + 0.8,
    );
    near(
        next.report().partial_slip.contact.tangent_force_n[0],
        -10.0 * mu,
    );
    near(
        next.report()
            .partial_slip
            .material
            .as_ref()
            .unwrap()
            .speed_mps,
        1.0,
    );
    near(next.report().external_heat_j, 0.0);
    assert!(next.report().total_energy_residual_j.abs() < 1e-8);
    assert_eq!(next.contact().rolling_state().accepted_steps(), 1);
}

#[test]
fn torsional_microslip_and_rolling_have_separate_work_and_heat() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let partial = contact(&thermal, &curves, model(0.3, zero_slopes()), [0.0; 3]);
    let old = rolling(
        partial,
        2,
        600.0,
        "rolling",
        RollingTemperatureSurface::SurfaceA,
    );
    let mut spin = motion(0.0);
    spin.torsional_spin_rad_per_s = 0.1;
    let next = old
        .step_candidate(
            &thermal,
            1.0,
            frame(),
            spin,
            &work("torsion and rolling"),
            rolling_motion(),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    // Independent two-Jenkins twist return: ci=N*mu*r*(.15,.85), ki=(7.5,2.5).
    let capacity = 10.0 * 0.3 * 0.005_f64.sqrt();
    let mut physical = 0.0;
    let mut stored = 0.0;
    for (fraction, stiffness) in [(0.15, 7.5), (0.85, 2.5)] {
        let c = fraction * capacity;
        physical += c * (0.1 - c / stiffness);
        stored += 0.5 * c * c / stiffness;
    }
    near(
        next.report().partial_slip.contact.dissipation.heat_j,
        physical,
    );
    near(next.report().partial_slip.contact.stored_energy_j, stored);
    near(
        next.report()
            .partial_slip
            .contact
            .generalized_work
            .numerical_dissipation_j,
        stored,
    );
    near(next.report().total_mechanical_work_j, capacity * 0.1 + 0.2);
    near(next.report().rolling().wrench.spin_moment_nm, 0.0);
    for (i, cp) in [10.0, 20.0].into_iter().enumerate() {
        let slip_fraction = [0.25, 0.75][i];
        let rolling_fraction = [0.6, 0.4][i];
        near(
            next.contact().partial_slip().temperature_k()[i],
            300.0 + (slip_fraction * physical + rolling_fraction * 0.2) / cp,
        );
    }
}

#[test]
fn late_original_rolling_domain_refusal_and_work_overlap_leave_both_owners_unchanged() {
    let curves = curves();
    let thermal = initial(&curves, [300.0; 2], 0.0);
    let make = || {
        contact(
            &thermal,
            &curves,
            model(1000.0, zero_slopes()),
            [0.1, 0.0, 0.0],
        )
    };
    let old = rolling(
        make(),
        2,
        300.005,
        "rolling",
        RollingTemperatureSurface::SurfaceA,
    );
    let before = old.partial_slip().checkpoint().clone();
    let attempt = || {
        old.step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.0),
            &work("late domain"),
            rolling_motion(),
            controls(),
            &[0.0; 2],
            other,
            || false,
        )
    };
    assert!(matches!(
        attempt(),
        Err(CaloricRollingError::Contact(CaloricContactError::Source(
            RollingEndpointError::Native(RollingLossError::OutsideApplicability { .. })
        )))
    ));
    assert_eq!(old.partial_slip().checkpoint(), &before);
    assert_eq!(old.rolling_state(), &RollingLossState::zero());
    assert_eq!(thermal.time_s(), 0.0);
    assert!(attempt().is_err());

    let collision = rolling(
        make(),
        2,
        600.0,
        "twist",
        RollingTemperatureSurface::SurfaceA,
    );
    assert!(matches!(
        collision.step_candidate(
            &thermal,
            1.0,
            frame(),
            motion(0.0),
            &work("same coordinate"),
            rolling_motion(),
            controls(),
            &[0.0; 2],
            other,
            || false,
        ),
        Err(CaloricRollingError::Invalid(_))
    ));

    let admitted = rolling(
        make(),
        2,
        600.0,
        "rolling",
        RollingTemperatureSurface::SurfaceA,
    );
    assert!(
        admitted
            .step_candidate(
                &thermal,
                1.0,
                frame(),
                motion(0.0),
                &work("cancel"),
                rolling_motion(),
                controls(),
                &[0.0; 2],
                other,
                || true,
            )
            .is_err()
    );
    let run = || {
        admitted
            .step_candidate(
                &thermal,
                1.0,
                frame(),
                motion(0.0),
                &work("replay"),
                rolling_motion(),
                controls(),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap()
    };
    let a = run();
    let b = run();
    assert_eq!(a.contact().rolling_state(), b.contact().rolling_state());
    assert_eq!(
        a.contact().partial_slip().checkpoint(),
        b.contact().partial_slip().checkpoint()
    );
    assert_eq!(
        a.contact().partial_slip().temperature_k(),
        b.contact().partial_slip().temperature_k()
    );
}
