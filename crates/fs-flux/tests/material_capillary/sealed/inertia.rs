//! Inertial liquid motion uses the original NASA-9 gas and sourced wetting.
use super::super::contact_line::{limits, sourced_wetting};
use super::*;
use fs_flux::capillary::ContactLineMotion;
use fs_flux::capillary::inertia::InertialCapillaryEnd;
use fs_flux::capillary::material::inertia::CapillaryInertiaParameters;
use fs_flux::capillary::material::sealed::advance::SealedCapillaryAdvanceControls;
use fs_flux::capillary::material::sealed::inertia::SealedCapillaryInertiaAccuracy;

fn parameters(g: SealedCapillaryGeometry) -> CapillaryInertiaParameters {
    CapillaryInertiaParameters {
        radius_m: g.radius_m,
        tube_length_m: g.tube_length_m,
        liquid_added_length_m: g.radius_m,
        displaced_added_length_m: 0.0,
        liquid_entry_loss_coefficient: 1.25,
        liquid_exit_loss_coefficient: 0.75,
        displaced_entry_loss_coefficient: 0.0,
        displaced_exit_loss_coefficient: 0.0,
        maximum_velocity_m_s: 0.01,
    }
}

fn accuracy() -> SealedCapillaryInertiaAccuracy {
    SealedCapillaryInertiaAccuracy {
        velocity_absolute_tolerance_m_s: 1e-10,
        stopping_time_absolute_tolerance_s: 1e-10,
    }
}

#[test]
fn manufactured_motion_retains_original_gas_and_every_mechanical_account() {
    let (liquid, _) = fixture(1000.0, 0.001, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let params = parameters(g);
    let original = gas(volume(g), 3.5);
    let area = PI * g.radius_m.powi(2);
    let capillary_scale = 2.0 * 0.03 / g.radius_m;
    let young_pressure = capillary_scale * 1.0_f64.cos();
    for (u0, u1, threshold_angle, loss_coefficient, branch) in [
        (
            0.0004_f64,
            0.0008_f64,
            1.2_f64,
            1.25_f64,
            ContactLineMotion::Advancing,
        ),
        (
            -0.0004_f64,
            -0.0008_f64,
            0.8_f64,
            0.75_f64,
            ContactLineMotion::Receding,
        ),
    ] {
        let owner = MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .with_contact_line(limits())
            .unwrap()
            .seal(original.clone(), g)
            .unwrap()
            .with_inertia(params, u0, accuracy())
            .unwrap();
        let elapsed = 0.002;
        let heat_power = 2e-7;
        let heat = elapsed * heat_power;
        let gravity = 9.81;
        let mean_velocity = 0.5 * (u0 + u1);
        let displacement = elapsed * mean_velocity;
        let final_length = g.wet_length_m + displacement;
        let midpoint_length = g.wet_length_m + 0.5 * displacement;
        let transferred = area * displacement;
        let v0 = original.volume_m3();
        let v1 = v0 - transferred;
        // Independent perfect-gas relation for this actual calorically perfect
        // NASA-9 fixture. Production still evaluates its original material law.
        let denominator = 1.4 * v1 + 0.6 * v0;
        let pm = (2.0 * 100000.0 * v0 + 0.4 * heat) / denominator;
        let mass0 = area * 1000.0 * (g.wet_length_m + params.liquid_added_length_m);
        let mass1 = area * 1000.0 * (final_length + params.liquid_added_length_m);
        let inertial_force = 0.5 * (mass0 + mass1) * (u1 - u0) / elapsed
            + 0.25 * area * 1000.0 * (u1 * u1 + u0 * u0);
        let pin_pressure = capillary_scale * (1.0_f64.cos() - threshold_angle.cos());
        let reservoir_force =
            0.5 * area * 1000.0 * loss_coefficient * mean_velocity.abs() * mean_velocity;
        let inlet = pm - young_pressure
            + 1000.0 * gravity * midpoint_length
            + pin_pressure
            + (inertial_force
                + 8.0 * PI * 0.001 * midpoint_length * mean_velocity
                + 2.0 * PI * g.radius_m * 2.0 * mean_velocity
                + reservoir_force)
                / area;
        let candidate = owner
            .step_candidate(
                SealedCapillaryInput {
                    inlet_pressure_pa: inlet,
                    duration_s: elapsed,
                    gas_heat_power_w: heat_power,
                    opposing_gravity_m_s2: gravity,
                },
                controls(),
                || false,
            )
            .unwrap();
        let step = candidate.step();
        near(step.wet_length_m, final_length, 2e-12);
        near(step.velocity_m_s.unwrap(), u1, 1e-10);
        assert_eq!(step.initial_velocity_m_s, Some(u0));
        assert_eq!(step.inertial_end, Some(InertialCapillaryEnd::TimeLimit));
        assert_eq!(step.elapsed_time_s, elapsed);
        assert_eq!(step.contact_line_motion, Some(branch));
        near(step.transferred_liquid_volume_m3, transferred, 1e-20);
        near(step.transferred_liquid_mass_kg, 1000.0 * transferred, 1e-17);
        near(step.midpoint_gas_pressure_pa, pm, 1e-6);
        near(step.gas.work_by_gas_j, -pm * transferred, 1e-13);
        near(
            step.gas.stored_energy_change_j,
            heat + pm * transferred,
            1e-13,
        );
        near(
            step.gas.pressure_final_volume_derivative_pa_m3,
            -1.4 * pm / denominator,
            100.0,
        );
        near(
            step.gas.pressure_heat_derivative_pa_j,
            0.4 / denominator,
            0.001,
        );
        near(
            step.kinetic_energy_change_j,
            0.5 * (mass1 * u1 * u1 - mass0 * u0 * u0),
            1e-19,
        );
        near(
            step.liquid_viscous_dissipation_j,
            8.0 * PI * 0.001 * midpoint_length * mean_velocity.powi(2) * elapsed,
            1e-19,
        );
        near(
            step.contact_line_viscous_dissipation_j,
            2.0 * PI * g.radius_m * 2.0 * mean_velocity.powi(2) * elapsed,
            1e-19,
        );
        near(
            step.pinning_dissipation_j,
            pin_pressure * transferred,
            1e-18,
        );
        near(
            step.reservoir_dissipation_j,
            0.5 * area * 1000.0 * loss_coefficient * mean_velocity.abs().powi(3) * elapsed,
            1e-22,
        );
        near(
            step.wall_surface_energy_change_j,
            -young_pressure * transferred,
            1e-18,
        );
        near(
            step.gravitational_energy_change_j,
            1000.0 * gravity * midpoint_length * transferred,
            1e-18,
        );
        assert!(step.liquid_viscous_dissipation_j > 0.0);
        assert!(step.contact_line_viscous_dissipation_j > 0.0);
        assert!(step.pinning_dissipation_j > 0.0);
        assert!(step.reservoir_dissipation_j > 0.0);
        near(
            step.inlet_pressure_work_j + heat,
            step.gas.stored_energy_change_j
                + step.kinetic_energy_change_j
                + step.wall_surface_energy_change_j
                + step.gravitational_energy_change_j
                + step.liquid_viscous_dissipation_j
                + step.contact_line_viscous_dissipation_j
                + step.pinning_dissipation_j
                + step.reservoir_dissipation_j,
            1e-13,
        );
        assert!(step.energy_balance_residual_j.abs() < 1e-13);
        assert_eq!(candidate.endpoint().inertia_parameters(), Some(&params));
        assert_eq!(candidate.endpoint().inertia_accuracy(), Some(accuracy()));
        assert_eq!(candidate.endpoint().gas().mass_kg(), original.mass_kg());
        assert_eq!(candidate.endpoint().gas().model(), original.model());
        assert_eq!(
            candidate
                .endpoint()
                .gas()
                .point()
                .axes()
                .get("specimen-lot"),
            Some(&42.0)
        );
        assert!(std::ptr::eq(candidate.endpoint().tube().liquid(), &liquid));
        assert!(std::ptr::eq(
            candidate.endpoint().tube().wetting(),
            &wetting
        ));
        assert_eq!(owner.gas(), &original);
        assert_eq!(owner.geometry(), g);
        assert_eq!(owner.velocity_m_s(), Some(u0));
    }
}

#[test]
fn manufactured_stop_uses_actual_elapsed_heat_and_returns_the_original_gas_candidate() {
    let (liquid, _) = fixture(1000.0, 0.001, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let params = parameters(g);
    let original = gas(volume(g), 3.5);
    let u0: f64 = 0.001;
    let stop_time: f64 = 1e-4;
    let requested_time = 1.5 * stop_time;
    let heat_power = 1e-7;
    let mean_velocity = 0.5 * u0;
    let displacement = stop_time * mean_velocity;
    let area = PI * g.radius_m.powi(2);
    let final_length = g.wet_length_m + displacement;
    let midpoint_length = g.wet_length_m + 0.5 * displacement;
    let transferred = area * displacement;
    let v0 = original.volume_m3();
    let pm = (2.0 * 100000.0 * v0 + 0.4 * heat_power * stop_time)
        / (1.4 * (v0 - transferred) + 0.6 * v0);
    let mass0 = area * 1000.0 * (g.wet_length_m + params.liquid_added_length_m);
    let inlet = pm - 2.0 * 0.03 * 1.2_f64.cos() / g.radius_m
        + (-mass0 * u0 / stop_time
            + 8.0 * PI * 0.001 * midpoint_length * mean_velocity
            + 2.0 * PI * g.radius_m * 2.0 * mean_velocity
            + 0.5 * area * 1000.0 * params.liquid_entry_loss_coefficient * mean_velocity.powi(2))
            / area;
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap()
        .with_inertia(params, u0, accuracy())
        .unwrap();
    let candidate = owner
        .step_candidate(
            SealedCapillaryInput {
                inlet_pressure_pa: inlet,
                duration_s: requested_time,
                gas_heat_power_w: heat_power,
                opposing_gravity_m_s2: 0.0,
            },
            controls(),
            || false,
        )
        .unwrap();
    let step = candidate.step();
    assert_eq!(step.inertial_end, Some(InertialCapillaryEnd::Stopped));
    assert_eq!(step.velocity_m_s, Some(0.0));
    assert_eq!(candidate.endpoint().velocity_m_s(), Some(0.0));
    near(step.elapsed_time_s, stop_time, 1e-10);
    assert!(step.elapsed_time_s < requested_time);
    near(step.wet_length_m, final_length, 1e-12);
    assert_eq!(step.gas.supplied_heat_j, heat_power * step.elapsed_time_s);
    assert!(step.gas.supplied_heat_j < heat_power * requested_time);
    near(step.kinetic_energy_change_j, -0.5 * mass0 * u0 * u0, 1e-19);
    near(
        step.gas.stored_energy_change_j,
        heat_power * step.elapsed_time_s
            + step.midpoint_gas_pressure_pa * step.transferred_liquid_volume_m3,
        1e-13,
    );
    let direct = original
        .midpoint_heat_volume(
            candidate.endpoint().gas().volume_m3(),
            step.gas.supplied_heat_j,
            controls().gas,
            || false,
        )
        .unwrap();
    assert_eq!(candidate.endpoint().gas(), direct.endpoint());
    assert_eq!(&step.gas, direct.step());
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), g);
    assert_eq!(owner.velocity_m_s(), Some(u0));
}

#[test]
fn adaptive_coasting_inside_static_band_stops_once_then_completes_the_heated_interval() {
    let (liquid, _) = fixture(1000.0, 0.001, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let params = parameters(g);
    let original = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap()
        .with_inertia(params, 0.001, accuracy())
        .unwrap();
    let input = SealedCapillaryInput {
        inlet_pressure_pa: 100000.0 - 2.0 * 0.03 * 1.0_f64.cos() / g.radius_m + 50.0,
        duration_s: 0.001,
        gas_heat_power_w: 1e-7,
        opposing_gravity_m_s2: 0.0,
    };
    let candidate = owner
        .advance_candidate(
            input,
            controls(),
            SealedCapillaryAdvanceControls {
                initial_step_s: input.duration_s,
                minimum_step_s: 1e-10,
                maximum_step_s: input.duration_s,
                maximum_attempts: 4000,
                relative_tolerance: 1e-6,
                length_absolute_tolerance_m: 1e-11,
                temperature_absolute_tolerance_k: 1e-6,
                pressure_absolute_tolerance_pa: 1e-3,
            },
            || false,
        )
        .unwrap();
    let step = candidate.step();
    assert_eq!(step.elapsed_time_s, input.duration_s);
    assert_eq!(step.velocity_m_s, Some(0.0));
    assert_eq!(candidate.endpoint().velocity_m_s(), Some(0.0));
    assert_eq!(step.stopping_events, 1);
    assert!(step.wet_length_m > g.wet_length_m);
    assert!(step.rejected_intervals > 0);
    assert!(step.maximum_normalized_local_error <= 1.0);
    assert!(step.step_evaluations <= 5 * step.attempts);
    near(
        step.gas_supplied_heat_j,
        input.gas_heat_power_w * input.duration_s,
        1e-20,
    );
    near(
        step.gas_stored_energy_change_j + step.gas_work_by_gas_j,
        step.gas_supplied_heat_j,
        1e-13,
    );
    assert!(step.kinetic_energy_change_j < 0.0);
    assert!(step.pinning_dissipation_j > 0.0);
    assert!(step.reservoir_dissipation_j > 0.0);
    assert!(step.energy_balance_residual_j.abs() < 1e-13);
    let pinned = candidate
        .endpoint()
        .step_candidate(input, controls(), || false)
        .unwrap();
    assert_eq!(
        pinned.step().contact_line_motion,
        Some(ContactLineMotion::Pinned)
    );
    assert_eq!(pinned.step().wet_length_m, step.wet_length_m);
    assert_eq!(pinned.step().velocity_m_s, Some(0.0));
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), g);
    assert_eq!(owner.velocity_m_s(), Some(0.001));
}

#[test]
fn invalid_bindings_source_refusals_and_cancellation_leave_the_owner_unchanged() {
    let (liquid, _) = fixture(1000.0, 0.001, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let params = parameters(g);
    let original = gas(volume(g), 3.5);
    let plain = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap();
    for invalid in [
        CapillaryInertiaParameters {
            radius_m: 2.0 * g.radius_m,
            ..params
        },
        CapillaryInertiaParameters {
            displaced_added_length_m: 1e-4,
            ..params
        },
        CapillaryInertiaParameters {
            displaced_entry_loss_coefficient: 1.0,
            ..params
        },
        CapillaryInertiaParameters {
            displaced_exit_loss_coefficient: 1.0,
            ..params
        },
        CapillaryInertiaParameters {
            liquid_entry_loss_coefficient: -1.0,
            ..params
        },
    ] {
        assert!(matches!(
            plain.clone().with_inertia(invalid, 0.0, accuracy()),
            Err(SealedCapillaryError::Invalid(_))
        ));
    }
    assert!(matches!(
        plain.clone().with_inertia(
            CapillaryInertiaParameters {
                maximum_velocity_m_s: 0.1,
                ..params
            },
            0.02,
            accuracy()
        ),
        Err(SealedCapillaryError::Invalid(_))
    ));
    assert_eq!(plain.velocity_m_s(), None);
    let owner = plain.with_inertia(params, 0.001, accuracy()).unwrap();
    let input = SealedCapillaryInput {
        duration_s: 1e-4,
        ..boundary()
    };
    assert!(matches!(
        owner.step_candidate(input, controls(), || true),
        Err(SealedCapillaryError::Cancelled)
    ));
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                gas_heat_power_w: 1.0,
                ..input
            },
            SealedCapillarySolve {
                gas: GasEnergySolve {
                    maximum_temperature_k: 300.0001,
                    ..controls().gas
                },
                ..controls()
            },
            || false,
        ),
        Err(SealedCapillaryError::Gas(_))
    ));
    assert!(matches!(
        owner.step_candidate(
            SealedCapillaryInput {
                duration_s: 1000.0,
                ..input
            },
            controls(),
            || false
        ),
        Err(SealedCapillaryError::StepTooLarge(_))
    ));
    let dry = SealedCapillaryGeometry {
        wet_length_m: 0.0,
        ..g
    };
    assert!(matches!(
        MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .seal(gas(volume(dry), 3.5), dry)
            .unwrap()
            .with_inertia(
                CapillaryInertiaParameters {
                    liquid_added_length_m: 0.0,
                    ..parameters(dry)
                },
                0.0,
                accuracy()
            ),
        Err(SealedCapillaryError::Invalid(_))
    ));
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), g);
    assert_eq!(owner.velocity_m_s(), Some(0.001));
}

#[test]
fn heating_cannot_hide_or_reverse_the_original_depinning_branch() {
    let (liquid, _) = fixture(1000.0, 0.001, 0.03, 1.0);
    let wetting = sourced_wetting(&liquid, 1.0, 2.0, 0.8, 1.2);
    let g = geometry();
    let original = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(original.clone(), g)
        .unwrap()
        .with_inertia(parameters(g), 0.0, accuracy())
        .unwrap();
    let scale = 2.0 * 0.03 / g.radius_m;
    let advancing_threshold = scale * (1.0_f64.cos() - 1.2_f64.cos());
    // Initial drive lies 20 Pa above advancing onset. A 60 Pa final
    // pressure rise places both midpoint and endpoint in the static band.
    // A 500 Pa rise instead selects a receding midpoint branch. Neither
    // trial may erase the original advancing excursion from rest.
    for pressure_rise in [60.0, 500.0] {
        let input = SealedCapillaryInput {
            inlet_pressure_pa: 100000.0 - scale * 1.0_f64.cos() + advancing_threshold + 20.0,
            duration_s: 0.001,
            gas_heat_power_w: 2.5 * original.volume_m3() * pressure_rise / 0.001,
            opposing_gravity_m_s2: 0.0,
        };
        assert!(matches!(
            owner.step_candidate(input, controls(), || false),
            Err(SealedCapillaryError::StepTooLarge(_))
        ));
    }
    assert_eq!(owner.gas(), &original);
    assert_eq!(owner.geometry(), g);
    assert_eq!(owner.velocity_m_s(), Some(0.0));
}
