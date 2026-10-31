//! Source-bound wetting coupled to the existing NASA-9 gas energy owner.
use super::super::contact_line::{limits, sourced_wetting};
use super::*;
use fs_flux::capillary::ContactLineMotion;
use fs_flux::capillary::material::contact_line::ContactLineTubeLimits;
use fs_flux::capillary::material::sealed::advance::SealedCapillaryAdvanceControls;

#[test]
fn pinned_gas_heats_at_fixed_volume_then_thermally_depins() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let original = gas(volume(g), 3.5);
    let mut owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap();
    let young_pressure = 2.0 * 0.03 * 1.0_f64.cos() / g.radius_m;
    // cp/R=3.5 gives cv/R=2.5. Each fixed-volume pulse raises pressure by
    // 40 Pa: the first two endpoints stay inside the receding threshold.
    let heat = 2.5 * original.volume_m3() * 40.0;
    let input = SealedCapillaryInput {
        inlet_pressure_pa: 100000.0 - young_pressure,
        gas_heat_power_w: heat,
        ..boundary()
    };
    for pulse in 1..=2 {
        let candidate = owner.step_candidate(input, controls(), || false).unwrap();
        let step = candidate.step();
        assert_eq!(step.contact_line_motion, Some(ContactLineMotion::Pinned));
        assert_eq!(step.transferred_liquid_volume_m3, 0.0);
        assert_eq!(step.wet_length_m, g.wet_length_m);
        assert_eq!(step.contact_line_viscous_dissipation_j, 0.0);
        assert_eq!(step.pinning_dissipation_j, 0.0);
        assert_eq!(step.liquid_viscous_dissipation_j, 0.0);
        assert_eq!(candidate.endpoint().gas().volume_m3(), original.volume_m3());
        near(
            candidate.endpoint().gas().state().state().pressure,
            100000.0 + 40.0 * f64::from(pulse),
            1e-6,
        );
        near(step.gas.stored_energy_change_j, heat, 1e-13);
        assert!(step.energy_balance_residual_j.abs() < 1e-13);
        assert!(step.apparent_contact_angle_radians.unwrap() < 1.0);
        owner = candidate.into_parts().0;
    }
    let candidate = owner.step_candidate(input, controls(), || false).unwrap();
    assert_eq!(
        candidate.step().contact_line_motion,
        Some(ContactLineMotion::Receding)
    );
    assert!(candidate.step().transferred_liquid_volume_m3 < 0.0);
    assert!(candidate.step().contact_line_viscous_dissipation_j > 0.0);
    assert!(candidate.step().pinning_dissipation_j > 0.0);
    assert!(candidate.endpoint().gas().state().state().temperature > 300.0);
    assert_eq!(candidate.endpoint().gas().mass_kg(), original.mass_kg());
    assert_eq!(candidate.endpoint().gas().model(), original.model());
    assert!(std::ptr::eq(
        candidate.endpoint().tube().wetting(),
        &wetting
    ));
}

#[test]
fn original_gas_energy_closes_with_separate_moving_line_and_pinning_losses() {
    let (liquid, _) = fixture(1000.0, 0.1, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let original = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap();
    let area = PI * g.radius_m.powi(2);
    let capillary_scale = 2.0 * 0.03 / g.radius_m;
    for (displacement, threshold_angle, motion) in [
        (0.0002, 1.2_f64, ContactLineMotion::Advancing),
        (-0.0002, 0.8_f64, ContactLineMotion::Receding),
    ] {
        let heat = 1e-6;
        let gravity = 9.81;
        let transferred = area * displacement;
        let midpoint = g.wet_length_m + 0.5 * displacement;
        let v0 = original.volume_m3();
        let v1 = v0 - transferred;
        // Independent ideal-gas midpoint relation for the actual NASA-9 fixture.
        let pm = (2.0 * 100000.0 * v0 + 0.4 * heat) / (1.4 * v1 + 0.6 * v0);
        let threshold = capillary_scale * (1.0_f64.cos() - threshold_angle.cos());
        let pin = pm - capillary_scale * 1.0_f64.cos()
            + 1000.0 * gravity * midpoint
            + 8.0 * 0.1 * midpoint * displacement / g.radius_m.powi(2)
            + 2.0 * 2.0 * displacement / g.radius_m
            + threshold;
        let candidate = owner
            .step_candidate(
                SealedCapillaryInput {
                    inlet_pressure_pa: pin,
                    gas_heat_power_w: heat,
                    opposing_gravity_m_s2: gravity,
                    ..boundary()
                },
                controls(),
                || false,
            )
            .unwrap();
        let step = candidate.step();
        near(step.wet_length_m, g.wet_length_m + displacement, 2e-12);
        near(step.midpoint_gas_pressure_pa, pm, 1e-6);
        assert_eq!(step.contact_line_motion, Some(motion));
        near(
            step.contact_line_viscous_dissipation_j,
            2.0 * PI * g.radius_m * 2.0 * displacement.powi(2),
            1e-16,
        );
        near(step.pinning_dissipation_j, threshold * transferred, 1e-16);
        near(
            step.liquid_viscous_dissipation_j,
            8.0 * PI * 0.1 * midpoint * displacement.powi(2),
            1e-16,
        );
        near(
            step.wall_surface_energy_change_j,
            -capillary_scale * 1.0_f64.cos() * transferred,
            1e-16,
        );
        near(
            step.apparent_contact_angle_radians.unwrap().cos(),
            threshold_angle.cos() - 2.0 * displacement / 0.03,
            1e-10,
        );
        near(
            step.gas.stored_energy_change_j,
            heat + pm * transferred,
            1e-13,
        );
        near(
            step.inlet_pressure_work_j + heat,
            step.gas.stored_energy_change_j
                + step.wall_surface_energy_change_j
                + step.gravitational_energy_change_j
                + step.liquid_viscous_dissipation_j
                + step.contact_line_viscous_dissipation_j
                + step.pinning_dissipation_j,
            1e-13,
        );
        assert!(step.energy_balance_residual_j.abs() < 1e-13);
        assert!(std::ptr::eq(
            candidate
                .endpoint()
                .tube()
                .contact_line_source()
                .unwrap()
                .wetting(),
            &wetting
        ));
        assert_eq!(
            candidate
                .endpoint()
                .gas()
                .point()
                .axes()
                .get("specimen-lot"),
            Some(&42.0)
        );
        assert_eq!(owner.gas(), &original);
        assert_eq!(owner.geometry(), g);
    }
}

#[test]
fn contact_line_limit_and_cancel_refusals_leave_sealed_sources_unchanged() {
    let (liquid, _) = fixture(1000.0, 0.1, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let original_liquid = liquid.clone();
    let original_wetting = wetting.clone();
    let g = geometry();
    let original = gas(volume(g), 3.5);
    // This accepted reference has mean speed below 1e-4 m/s, while its
    // independently calculated initial speed exceeds that bound. Admission
    // must check the physical endpoint state, not only the mean displacement.
    let reference = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap()
        .step_candidate(boundary(), controls(), || false)
        .unwrap();
    let mean_speed = (reference.step().wet_length_m - g.wet_length_m) / boundary().duration_s;
    let initial_speed = (2.0 * 0.03 * 1.2_f64.cos() / g.radius_m)
        / (8.0 * 0.1 * g.wet_length_m / g.radius_m.powi(2) + 2.0 * 2.0 / g.radius_m);
    assert!(mean_speed > 0.0 && mean_speed < 1e-4);
    assert!(initial_speed > 1e-4);
    for domain in [
        ContactLineTubeLimits {
            maximum_speed_m_s: 1e-4,
            ..limits()
        },
        ContactLineTubeLimits {
            maximum_speed_m_s: 1e-12,
            ..limits()
        },
        ContactLineTubeLimits {
            maximum_angle_departure_radians: 0.01,
            ..limits()
        },
    ] {
        let owner = MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .with_contact_line(domain)
            .unwrap()
            .seal(original.clone(), g)
            .unwrap();
        assert!(matches!(
            owner.step_candidate(boundary(), controls(), || false),
            Err(SealedCapillaryError::Invalid(_))
        ));
        assert!(matches!(
            owner.step_candidate(boundary(), controls(), || true),
            Err(SealedCapillaryError::Cancelled)
        ));
        assert_eq!(owner.gas(), &original);
        assert_eq!(owner.geometry(), g);
    }
    assert_eq!(liquid, original_liquid);
    assert_eq!(wetting, original_wetting);
}

#[test]
fn heated_pinning_and_advancing_overshoot_refuse_an_inconsistent_endpoint_branch() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let original = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap();
    let scale = 2.0 * 0.03 / g.radius_m;
    let young_pressure = scale * 1.0_f64.cos();
    // A strong heat pulse selects receding motion at the midpoint although
    // the original state is already advancing. One net displacement must not
    // erase the intervening reversal and its positive friction losses.
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                gas_heat_power_w: 2.5 * original.volume_m3() * 1000.0,
                ..boundary()
            },
            controls(),
            || false,
        ),
        Err(SealedCapillaryError::StepTooLarge(_))
    ));
    // The hypothetical pinned midpoint needs only 60 Pa of reaction, while
    // the endpoint needs 120 Pa, outside the 93.84 Pa receding threshold.
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                inlet_pressure_pa: 100000.0 - young_pressure,
                gas_heat_power_w: 2.5 * original.volume_m3() * 120.0,
                ..boundary()
            },
            controls(),
            || false
        ),
        Err(SealedCapillaryError::StepTooLarge(_))
    ));

    // A moving midpoint can also pass its static threshold without reversing
    // the Young drive. That still hides a stick/slip event and must subdivide.
    let displacement = 0.0002;
    let heat = 1e-6;
    let gravity = 9.81;
    let midpoint = g.wet_length_m + 0.5 * displacement;
    let v0 = original.volume_m3();
    let v1 = v0 - PI * g.radius_m.powi(2) * displacement;
    let pm = (2.0 * 100000.0 * v0 + 0.4 * heat) / (1.4 * v1 + 0.6 * v0);
    let threshold = scale * (1.0_f64.cos() - 1.2_f64.cos());
    let pin = pm - young_pressure
        + 1000.0 * gravity * midpoint
        + 8.0 * 0.05 * midpoint * displacement / g.radius_m.powi(2)
        + 2.0 * 2.0 * displacement / g.radius_m
        + threshold;
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                inlet_pressure_pa: pin,
                gas_heat_power_w: heat,
                opposing_gravity_m_s2: gravity,
                ..boundary()
            },
            controls(),
            || false
        ),
        Err(SealedCapillaryError::StepTooLarge(_))
    ));
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), g);
}

#[test]
fn adaptive_interval_resolves_thermal_depinning_and_sums_only_accepted_losses() {
    let (liquid, _) = fixture(1000.0, 0.05, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let original = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap();
    let scale = 2.0 * 0.03 / g.radius_m;
    let input = SealedCapillaryInput {
        inlet_pressure_pa: 100000.0 - scale * 1.0_f64.cos(),
        duration_s: 4.0,
        gas_heat_power_w: 2.5 * original.volume_m3() * 40.0,
        ..boundary()
    };
    let adaptive = SealedCapillaryAdvanceControls {
        initial_step_s: 4.0,
        minimum_step_s: 1e-8,
        maximum_step_s: 4.0,
        maximum_attempts: 2000,
        relative_tolerance: 1e-6,
        length_absolute_tolerance_m: 1e-10,
        temperature_absolute_tolerance_k: 1e-6,
        pressure_absolute_tolerance_pa: 1e-3,
    };
    let result = owner
        .advance_candidate(input, controls(), adaptive, || false)
        .unwrap();
    let step = result.step();
    assert_eq!(step.elapsed_time_s, input.duration_s);
    assert!(step.accepted_intervals > 1);
    assert!(step.rejected_intervals > 0);
    assert_eq!(
        step.attempts,
        step.accepted_intervals + step.rejected_intervals
    );
    assert!(step.step_evaluations <= 3 * step.attempts);
    assert!(step.maximum_normalized_local_error <= 1.0);
    assert!(step.wet_length_m < g.wet_length_m);
    assert!(step.contact_line_viscous_dissipation_j > 0.0);
    assert!(step.pinning_dissipation_j > 0.0);
    let receding_threshold = scale * (0.8_f64.cos() - 1.0_f64.cos());
    near(
        step.pinning_dissipation_j,
        -receding_threshold * step.transferred_liquid_volume_m3,
        1e-18,
    );
    near(
        step.gas_supplied_heat_j,
        input.gas_heat_power_w * input.duration_s,
        1e-19,
    );
    near(
        step.transferred_liquid_volume_m3,
        PI * g.radius_m.powi(2) * (step.wet_length_m - g.wet_length_m),
        1e-20,
    );
    near(
        step.inlet_pressure_work_j + step.gas_supplied_heat_j,
        step.gas_stored_energy_change_j
            + step.wall_surface_energy_change_j
            + step.gravitational_energy_change_j
            + step.liquid_viscous_dissipation_j
            + step.contact_line_viscous_dissipation_j
            + step.pinning_dissipation_j,
        1e-13,
    );
    assert!(step.energy_balance_residual_j.abs() < 1e-13);
    assert_eq!(result.endpoint().gas().mass_kg(), original.mass_kg());
    assert_eq!(result.endpoint().gas().model(), original.model());
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), g);
}
