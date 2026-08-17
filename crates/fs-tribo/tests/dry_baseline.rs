//! G0/G1/G3 public conformance for the dependency-independent dry baseline.
//!
//! Every coefficient and thermal property here is synthetic. These tests check
//! equations, refusal, and replay; they do not validate a material system.

use fs_tribo::{
    ApplicabilityRange, ArchardLaw, ContactFrame, DryFrictionApplicability, DryFrictionState,
    DryInterfaceSystemCard, FlashTemperatureEstimate, FlashTemperatureInput,
    FlashTemperatureUnknown, FrictionLaw, FrictionRegime, HeatPartition, HertzCylinderPlane,
    HertzSpherePlane, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    SurfaceThermalProperties, TangentialSlip, TriboError, WearState, WorkLedger,
    flash_temperature_candidate,
};

fn interface() -> InterfaceSystemRef {
    InterfaceSystemRef::new(
        "fixture/steel-coated-a->ceramic-b",
        "fixture/dry-history-v1",
        "fixture/synthetic-dry-card-v1",
        InputAuthority::SyntheticFixture,
        InterfaceMedium::Dry,
    )
    .expect("synthetic dry identity")
}

fn slip(speed_mps: f64) -> TangentialSlip {
    let frame = ContactFrame::new([0.0, 0.0, 1.0]).expect("unit normal");
    TangentialSlip::new(&frame, [speed_mps, 0.0, 0.0]).expect("tangent velocity")
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1.0e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g0_card_requires_ordered_domain_and_matching_slip_speed() {
    let card = DryInterfaceSystemCard::new(
        interface(),
        FrictionLaw::Coulomb {
            static_mu: 0.6,
            kinetic_mu: 0.4,
        },
        DryFrictionApplicability::new(
            ApplicabilityRange::new(280.0, 320.0).unwrap(),
            ApplicabilityRange::new(1.0e5, 5.0e6).unwrap(),
            ApplicabilityRange::new(0.0, 2.0).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let state = DryFrictionState::new(300.0, 2.0e6, 1.5).unwrap();
    let receipt = card.query(state, 10.0, slip(1.5)).unwrap();
    assert_eq!(receipt.response().regime, FrictionRegime::Sliding);
    assert_eq!(
        receipt.response().provenance().authority(),
        InputAuthority::SyntheticFixture
    );
    assert!(matches!(
        card.query(
            DryFrictionState::new(321.0, 2.0e6, 1.5).unwrap(),
            10.0,
            slip(1.5)
        ),
        Err(TriboError::OutsideApplicability {
            field: "temperature_kelvin",
            ..
        })
    ));
    assert!(matches!(
        card.query(
            DryFrictionState::new(300.0, 2.0e6, 1.0).unwrap(),
            10.0,
            slip(1.5)
        ),
        Err(TriboError::InvalidInput {
            field: "state.slip_speed_mps"
        })
    ));
}

#[test]
fn g1_hertz_and_block_on_incline_match_closed_forms() {
    let sphere = HertzSpherePlane {
        effective_radius: 0.02,
        reduced_modulus: 200.0e9,
    }
    .response(&interface(), 2.0e-6)
    .unwrap();
    close(sphere.contact_radius_m, 2.0e-4);
    close(sphere.normal_force_n, 106.666_666_666_666_67);

    let cylinder = HertzCylinderPlane {
        effective_radius: 0.02,
        reduced_modulus: 200.0e9,
    }
    .response(&interface(), 1_000.0)
    .unwrap();
    close(cylinder.half_width_m, 1.128_379_167_095_512_6e-5);
    close(cylinder.peak_pressure_pa, 56_418_958.354_775_63);

    let static_mu = 0.5;
    let law = FrictionLaw::Coulomb {
        static_mu,
        kinetic_mu: 0.4,
    };
    let mass_kg = 3.0;
    let gravity_m_per_s2 = 9.806_65;
    let theta = static_mu.atan(); // det-ok: test-only incline oracle
    let normal_force = mass_kg * gravity_m_per_s2 * theta.cos(); // det-ok: test-only incline oracle
    let downslope_force = mass_kg * gravity_m_per_s2 * theta.sin(); // det-ok: test-only incline oracle
    let response = law.evaluate(&interface(), normal_force, slip(0.0)).unwrap();
    close(response.static_limit, downslope_force);
}

#[test]
fn g3_reversal_scaling_and_replay_are_deterministic() {
    let law = FrictionLaw::Stribeck {
        static_mu: 0.7,
        kinetic_mu: 0.3,
        characteristic_speed: 1.0,
        viscous_per_speed: 0.02,
    };
    let forward = law.evaluate(&interface(), 10.0, slip(0.5)).unwrap();
    let reverse = law.evaluate(&interface(), 10.0, slip(-0.5)).unwrap();
    let scaled = law.evaluate(&interface(), 20.0, slip(0.5)).unwrap();
    close(reverse.traction_n()[0], -forward.traction_n()[0]);
    close(reverse.dissipated_power_w(), forward.dissipated_power_w());
    close(scaled.traction_n()[0], 2.0 * forward.traction_n()[0]);

    let law = ArchardLaw {
        wear_coefficient: 2.0e-6,
    };
    let mut first = WearState::default();
    let mut second = WearState::default();
    for (force, distance) in [(100.0, 0.2), (200.0, 0.1), (50.0, 0.4)] {
        law.advance(&interface(), &mut first, force, distance, 2.0e9)
            .unwrap();
        law.advance(&interface(), &mut second, force, distance, 2.0e9)
            .unwrap();
    }
    assert_eq!(first, second);
}

#[test]
fn g0_heat_partition_ledger_and_flash_candidate_are_closed_or_unknown() {
    let partition = HeatPartition::new(0.5, 0.25, 0.25).unwrap();
    let step = fs_tribo::DissipationStep::from_power(100.0, 0.2, partition).unwrap();
    let mut ledger = WorkLedger::default();
    ledger.record(step).unwrap();
    close(ledger.dissipated_work_j(), 20.0);
    close(ledger.surface_a_heat_j(), 10.0);
    close(ledger.surface_b_heat_j(), 5.0);
    close(ledger.other_work_j(), 5.0);

    let surface_a = SurfaceThermalProperties::new(100.0, 1.0e-5).unwrap();
    let surface_b = SurfaceThermalProperties::new(50.0, 4.0e-6).unwrap();
    let input = FlashTemperatureInput::new(
        100.0,
        partition,
        0.01,
        0.002,
        0.5,
        Some(surface_a),
        Some(surface_b),
    )
    .unwrap();
    let estimate = flash_temperature_candidate(&interface(), input).unwrap();
    assert!(matches!(estimate, FlashTemperatureEstimate::Candidate(_)));
    let FlashTemperatureEstimate::Candidate(candidate) = estimate else {
        return;
    };
    close(candidate.dwell_time_s(), 0.004);
    close(candidate.surface_a_heat_flux_w_per_m2(), 5_000.0);
    close(candidate.surface_b_heat_flux_w_per_m2(), 2_500.0);
    assert!(candidate.surface_a_rise_k() > 0.0);
    assert!(candidate.surface_b_rise_k() > 0.0);
    assert_eq!(
        candidate.provenance().authority(),
        InputAuthority::SyntheticFixture
    );

    let missing_b =
        FlashTemperatureInput::new(100.0, partition, 0.01, 0.002, 0.5, Some(surface_a), None)
            .unwrap();
    assert_eq!(
        flash_temperature_candidate(&interface(), missing_b).unwrap(),
        FlashTemperatureEstimate::Unknown(
            FlashTemperatureUnknown::MissingSurfaceBThermalProperties
        )
    );
}

#[test]
fn g3_flash_candidate_scales_with_power() {
    let partition = HeatPartition::new(0.5, 0.5, 0.0).unwrap();
    let thermal = Some(SurfaceThermalProperties::new(60.0, 1.2e-5).unwrap());
    let candidate = |power| {
        let input =
            FlashTemperatureInput::new(power, partition, 0.01, 0.004, 1.0, thermal, thermal)
                .unwrap();
        let estimate = flash_temperature_candidate(&interface(), input).unwrap();
        assert!(matches!(estimate, FlashTemperatureEstimate::Candidate(_)));
        let FlashTemperatureEstimate::Candidate(value) = estimate else {
            return None;
        };
        Some(value)
    };
    let one = candidate(25.0).unwrap();
    let two = candidate(50.0).unwrap();
    close(two.surface_a_rise_k(), 2.0 * one.surface_a_rise_k());
    close(two.surface_b_rise_k(), 2.0 * one.surface_b_rise_k());
}
