use super::super::{CircularMeniscusInput, TwoFluidCircularMeniscusInput};
use super::*;

fn input() -> InertialCapillaryInput {
    InertialCapillaryInput {
        meniscus: InclinedTwoFluidCircularMeniscusInput {
            meniscus: TwoFluidCircularMeniscusInput {
                meniscus: CircularMeniscusInput {
                    dynamic_viscosity_pa_s: 0.0,
                    radius_m: 1.0 / PI.sqrt(),
                    tube_length_m: 10.0,
                    wet_length_m: 0.0,
                    surface_tension_n_m: 0.0,
                    contact_angle_cosine: 0.0,
                    applied_pressure_pa: 1.0,
                    duration_s: 0.01,
                },
                displaced_dynamic_viscosity_pa_s: 0.0,
            },
            density_kg_m3: 1.0,
            displaced_density_kg_m3: 0.0,
            opposing_gravity_m_s2: 0.0,
        },
        initial_velocity_m_s: 0.0,
        added_mass_kg: 1.0,
        advancing_reservoir_loss_kg_m: 0.5,
        receding_reservoir_loss_kg_m: 0.5,
        maximum_velocity_m_s: 10.0,
        contact_line: None,
    }
}
fn controls() -> InertialCapillarySolve {
    InertialCapillarySolve {
        maximum_iterations: 128,
        pressure_absolute_tolerance_pa: 1e-11,
        velocity_absolute_tolerance_m_s: 1e-12,
        time_absolute_tolerance_s: 1e-12,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + 2e-10 * expected.abs(),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn restart(input: &mut InertialCapillaryInput, step: InertialCapillaryStep) {
    input.meniscus.meniscus.meniscus.wet_length_m = step.wet_length_m;
    input.initial_velocity_m_s = step.velocity_m_s;
}

#[test]
fn inertial_dry_start_has_finite_acceleration_and_second_order_bosanquet_limit() {
    // In dimensionless A=rho=ell=drive=1, the corrected reservoir convention
    // gives x=sqrt(1+t²)-1, u=t/sqrt(1+t²), D_res=x²/[2(x+1)].
    let end_time = 3.0_f64.sqrt();
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let mut i = input();
        i.meniscus.meniscus.meniscus.duration_s = end_time / f64::from(n);
        let mut loss = 0.0;
        for _ in 0..n {
            let step = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
            assert_eq!(step.end, InertialCapillaryEnd::TimeLimit);
            assert!(step.velocity_m_s.is_finite() && step.velocity_m_s > 0.0);
            assert!(step.energy_balance_residual_j.abs() < 1e-12);
            loss += step.reservoir_dissipation_j;
            restart(&mut i, step);
        }
        let x = i.meniscus.meniscus.meniscus.wet_length_m;
        errors.push((x - 1.0).abs());
        if n == 128 {
            near(x, 1.0, 1.8e-5);
            near(i.initial_velocity_m_s, 3.0_f64.sqrt() / 2.0, 9e-7);
            near(loss, 0.25, 2e-5);
        }
    }
    assert!(errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9);
    let mut short = input();
    short.meniscus.meniscus.meniscus.duration_s = 1e-5;
    let step = step_inertial_circular_meniscus(short, controls(), || false).unwrap();
    near(step.velocity_m_s / 1e-5, 1.0, 2e-9);
    near(step.wet_length_m / (0.5e-10), 1.0, 2e-9);
}

#[test]
fn inviscid_recession_keeps_velocity_and_exports_exact_lost_column_energy() {
    let mut i = input();
    i.meniscus.meniscus.meniscus.wet_length_m = 2.0;
    i.meniscus.meniscus.meniscus.applied_pressure_pa = 0.0;
    i.meniscus.meniscus.meniscus.duration_s = 0.2;
    i.initial_velocity_m_s = -0.2;
    let s = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
    near(s.velocity_m_s, -0.2, 1e-13);
    near(s.wet_length_m, 1.96, 1e-13);
    near(s.kinetic_energy_change_j, -0.0008, 1e-14);
    near(s.reservoir_dissipation_j, 0.0008, 1e-14);
    assert_eq!(s.motion, ContactLineMotion::Receding);
    assert!(s.energy_balance_residual_j.abs() < 1e-14);
}

#[test]
fn manufactured_both_directions_close_independent_phase_line_pin_and_mass_accounts() {
    for sign in [-1.0, 1.0] {
        let mut i = input();
        i.meniscus.meniscus.meniscus.radius_m = 0.2;
        i.meniscus.meniscus.meniscus.tube_length_m = 1.0;
        i.meniscus.meniscus.meniscus.wet_length_m = 0.3;
        i.meniscus.meniscus.meniscus.dynamic_viscosity_pa_s = 0.04;
        i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = 0.08;
        i.meniscus.meniscus.meniscus.surface_tension_n_m = 0.1;
        i.meniscus.meniscus.meniscus.contact_angle_cosine = 1.0_f64.cos();
        i.meniscus.density_kg_m3 = 1.5;
        i.meniscus.displaced_density_kg_m3 = 0.7;
        i.meniscus.opposing_gravity_m_s2 = 0.3;
        i.added_mass_kg = 0.01;
        i.advancing_reservoir_loss_kg_m = 0.3;
        i.receding_reservoir_loss_kg_m = 0.2;
        i.initial_velocity_m_s = sign * 0.02;
        let u1 = sign * 0.024;
        let w = 0.5 * (i.initial_velocity_m_s + u1);
        let h = 0.05;
        i.meniscus.meniscus.meniscus.duration_s = h;
        i.contact_line = Some(InertialContactLine {
            friction_pa_s: 0.04,
            advancing_angle_radians: 1.2,
            receding_angle_radians: 0.8,
            maximum_angle_departure_radians: 0.5,
        });
        let area = PI * 0.2_f64.powi(2);
        let x0 = 0.3;
        let x1 = x0 + h * w;
        let xm = 0.5 * (x0 + x1);
        let m0 = area * (1.5 * x0 + 0.7 * (1.0 - x0)) + 0.01;
        let m1 = area * (1.5 * x1 + 0.7 * (1.0 - x1)) + 0.01;
        let theta = if sign > 0.0 { 1.2_f64 } else { 0.8_f64 };
        let threshold = (2.0 * 0.1 / 0.2) * (1.0_f64.cos() - theta.cos());
        let c = if sign > 0.0 { 0.3 } else { 0.2 };
        let bulk = 8.0 * PI * (0.04 * xm + 0.08 * (1.0 - xm));
        let line = 2.0 * PI * 0.2 * 0.04;
        let inertia = 0.5 * (m0 + m1) * (u1 - i.initial_velocity_m_s) / h
            + 0.25 * area * (1.5 - 0.7) * (u1 * u1 + i.initial_velocity_m_s.powi(2));
        i.meniscus.meniscus.meniscus.applied_pressure_pa =
            (inertia + (bulk + line) * w + c * w.abs() * w) / area + threshold
                - (2.0 * 0.1 / 0.2) * 1.0_f64.cos()
                + 0.3 * (1.5 * xm + 0.7 * (1.0 - xm));
        let s = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
        let volume = area * (x1 - x0);
        near(s.wet_length_m, x1, 2e-13);
        near(s.velocity_m_s, u1, 2e-12);
        near(s.effective_mass_kg, m1, 1e-15);
        near(
            s.kinetic_energy_change_j,
            0.5 * m1 * u1 * u1 - 0.5 * m0 * i.initial_velocity_m_s.powi(2),
            1e-15,
        );
        near(
            s.liquid_viscous_dissipation_j,
            8.0 * PI * 0.04 * xm * w * w * h,
            1e-15,
        );
        near(
            s.displaced_viscous_dissipation_j,
            8.0 * PI * 0.08 * (1.0 - xm) * w * w * h,
            1e-15,
        );
        near(
            s.contact_line_viscous_dissipation_j,
            line * w * w * h,
            1e-15,
        );
        near(s.pinning_dissipation_j, threshold * volume, 1e-15);
        near(s.reservoir_dissipation_j, c * w.abs().powi(3) * h, 1e-15);
        near(
            s.gravitational_storage_change_j,
            0.3 * (1.5 - 0.7) * xm * volume,
            1e-15,
        );
        near(s.outlet_potential_energy_j, 0.7 * 0.3 * volume, 1e-15);
        near(
            s.wall_surface_energy_change_j,
            -(2.0 * 0.1 / 0.2) * 1.0_f64.cos() * volume,
            1e-15,
        );
        assert!(s.energy_balance_residual_j.abs() < 1e-13);
    }
}

#[test]
fn constant_mass_resistance_matches_the_exact_midpoint_map_and_refines() {
    let mut i = input();
    i.meniscus.meniscus.meniscus.tube_length_m = 1.0;
    i.meniscus.meniscus.meniscus.wet_length_m = 0.1;
    i.meniscus.displaced_density_kg_m3 = 1.0;
    i.added_mass_kg = 0.0;
    i.meniscus.meniscus.meniscus.dynamic_viscosity_pa_s = 1.0 / (8.0 * PI);
    i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = 1.0 / (8.0 * PI);
    i.advancing_reservoir_loss_kg_m = 0.0;
    i.receding_reservoir_loss_kg_m = 0.0;
    i.initial_velocity_m_s = 0.2;
    i.meniscus.meniscus.meniscus.duration_s = 0.1;
    let s = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
    let u1 = ((1.0 - 0.05) * 0.2 + 0.1) / (1.0 + 0.05);
    near(s.velocity_m_s, u1, 1e-13);
    near(s.wet_length_m, 0.1 + 0.05 * (0.2 + u1), 1e-13);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut j = i;
        j.meniscus.meniscus.meniscus.duration_s = 0.5 / f64::from(n);
        for _ in 0..n {
            let s = step_inertial_circular_meniscus(j, controls(), || false).unwrap();
            restart(&mut j, s);
        }
        let exact = 1.0 + (0.2 - 1.0) * (-0.5_f64).exp();
        errors.push((j.initial_velocity_m_s - exact).abs());
    }
    assert!(errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9);
}

#[test]
fn physical_stop_publishes_zero_velocity_then_the_same_state_can_pin() {
    let mut i = input();
    i.meniscus.meniscus.meniscus.tube_length_m = 1.0;
    i.meniscus.meniscus.meniscus.wet_length_m = 0.2;
    i.meniscus.displaced_density_kg_m3 = 1.0;
    i.added_mass_kg = 0.0;
    i.advancing_reservoir_loss_kg_m = 0.0;
    i.receding_reservoir_loss_kg_m = 0.0;
    i.meniscus.meniscus.meniscus.surface_tension_n_m = 0.1;
    i.meniscus.meniscus.meniscus.contact_angle_cosine = 1.0_f64.cos();
    i.contact_line = Some(InertialContactLine {
        friction_pa_s: 0.0,
        advancing_angle_radians: 1.2,
        receding_angle_radians: 0.8,
        maximum_angle_departure_radians: 0.5,
    });
    let scale = 2.0 * 0.1 / i.meniscus.meniscus.meniscus.radius_m;
    let threshold = scale * (1.0_f64.cos() - 1.2_f64.cos());
    i.meniscus.meniscus.meniscus.applied_pressure_pa = 0.03 - scale * 1.0_f64.cos();
    i.initial_velocity_m_s = 0.02;
    i.meniscus.meniscus.meniscus.duration_s = 1.0;
    let expected_time = 0.02 / (threshold - 0.03);
    let s = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
    assert_eq!(s.end, InertialCapillaryEnd::Stopped);
    assert_eq!(s.velocity_m_s, 0.0);
    near(s.elapsed_time_s, expected_time, 1e-12);
    near(s.wet_length_m, 0.2 + 0.01 * expected_time, 1e-13);
    assert!(s.pinning_dissipation_j > 0.0);
    restart(&mut i, s);
    i.meniscus.meniscus.meniscus.duration_s = 1.0 - s.elapsed_time_s;
    let pinned = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
    assert_eq!(pinned.motion, ContactLineMotion::Pinned);
    assert_eq!(pinned.end, InertialCapillaryEnd::TimeLimit);
    assert_eq!(pinned.transferred_volume_m3, 0.0);
    assert_eq!(pinned.total_dissipation_j, 0.0);
}

#[test]
fn viscous_coasting_never_becomes_a_fabricated_finite_stop() {
    let mut i = input();
    i.meniscus.meniscus.meniscus.tube_length_m = 1.0;
    i.meniscus.meniscus.meniscus.wet_length_m = 0.2;
    i.meniscus.meniscus.meniscus.applied_pressure_pa = 0.0;
    i.meniscus.meniscus.meniscus.dynamic_viscosity_pa_s = 1.0 / (8.0 * PI);
    i.meniscus.meniscus.displaced_dynamic_viscosity_pa_s = 1.0 / (8.0 * PI);
    i.meniscus.displaced_density_kg_m3 = 1.0;
    i.added_mass_kg = 0.0;
    i.advancing_reservoir_loss_kg_m = 0.0;
    i.receding_reservoir_loss_kg_m = 0.0;
    i.initial_velocity_m_s = 0.2;
    i.meniscus.meniscus.meniscus.duration_s = 3.0;
    assert!(matches!(
        step_inertial_circular_meniscus(i, controls(), || false),
        Err(InertialCapillaryError::StepTooLarge(_))
    ));
}

#[test]
fn stopping_uses_trial_elapsed_pressure_and_retains_exact_accepted_payload() {
    let mut i = input();
    i.meniscus.meniscus.meniscus.tube_length_m = 1.0;
    i.meniscus.meniscus.meniscus.wet_length_m = 0.2;
    i.meniscus.meniscus.meniscus.applied_pressure_pa = 0.0;
    i.meniscus.displaced_density_kg_m3 = 1.0;
    i.added_mass_kg = 0.0;
    i.advancing_reservoir_loss_kg_m = 0.0;
    i.receding_reservoir_loss_kg_m = 0.0;
    i.initial_velocity_m_s = 0.4;
    i.meniscus.meniscus.meniscus.duration_s = 1.0;
    let s = solve_with_pressure(
        i,
        controls(),
        |x, t, _| {
            Ok::<_, ()>(Pressure {
                downstream_pressure_pa: 1.0 + t,
                endpoint_pressure_pa: 1.0 + 2.0 * t,
                derivative_wet_length_pa_m: 0.0,
                derivative_time_pa_s: 1.0,
                payload: (x, t),
            })
        },
        &mut || false,
    )
    .unwrap();
    let t = (2.6_f64.sqrt() - 1.0) / 2.0;
    assert_eq!(s.step.end, InertialCapillaryEnd::Stopped);
    near(s.step.elapsed_time_s, t, 1e-12);
    near(s.pressure.payload.1, t, 1e-12);
    assert_eq!(s.pressure.payload.0, s.step.wet_length_m);
    near(s.step.kinetic_energy_change_j, -0.08, 1e-13);
    near(s.step.pressure_work_j, -0.08, 1e-13);
}

#[test]
fn invalid_geometry_limits_cancellation_and_callback_work_are_bounded() {
    let mut i = input();
    i.added_mass_kg = 0.0;
    assert!(matches!(
        step_inertial_circular_meniscus(i, controls(), || false),
        Err(InertialCapillaryError::Invalid(_))
    ));
    i = input();
    assert_eq!(
        step_inertial_circular_meniscus(i, controls(), || true),
        Err(InertialCapillaryError::Cancelled)
    );
    let mut c = controls();
    c.maximum_iterations = 1;
    assert!(matches!(
        step_inertial_circular_meniscus(i, c, || false),
        Err(InertialCapillaryError::NotConverged { iterations: 1, .. })
    ));
    let original = i;
    let result = solve_with_pressure(
        i,
        controls(),
        |_, _, _| Err::<Pressure<()>, _>("source refusal"),
        &mut || false,
    );
    assert!(matches!(
        result,
        Err(InertialSolveError::Pressure("source refusal"))
    ));
    assert_eq!(i, original);
    i.meniscus.meniscus.meniscus.wet_length_m = 9.99;
    i.initial_velocity_m_s = 1.0;
    i.meniscus.meniscus.meniscus.duration_s = 1.0;
    assert!(matches!(
        step_inertial_circular_meniscus(i, controls(), || false),
        Err(InertialCapillaryError::StepTooLarge(_))
    ));
    i = input();
    i.initial_velocity_m_s = 0.2;
    i.meniscus.meniscus.meniscus.duration_s = 0.0;
    let identity = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
    assert_eq!(identity.wet_length_m, 0.0);
    assert_eq!(identity.velocity_m_s, 0.2);
    assert_eq!(identity.kinetic_energy_change_j, 0.0);
    assert_eq!(identity.total_dissipation_j, 0.0);
}

#[test]
fn equal_wetting_angles_survive_the_young_cosine_roundtrip() {
    for theta in [1.0_f64, PI / 2.0] {
        let mut i = input();
        i.meniscus.meniscus.meniscus.wet_length_m = 0.2;
        i.meniscus.meniscus.meniscus.surface_tension_n_m = 1.0;
        i.meniscus.meniscus.meniscus.contact_angle_cosine = theta.cos();
        let radius = i.meniscus.meniscus.meniscus.radius_m;
        i.meniscus.meniscus.meniscus.applied_pressure_pa = 1.0 - 2.0 * theta.cos() / radius;
        i.contact_line = Some(InertialContactLine {
            friction_pa_s: 0.1,
            advancing_angle_radians: theta,
            receding_angle_radians: theta,
            maximum_angle_departure_radians: 0.3,
        });
        let s = step_inertial_circular_meniscus(i, controls(), || false).unwrap();
        assert_eq!(s.pinning_dissipation_j, 0.0);
        assert!(s.contact_line_viscous_dissipation_j > 0.0);
    }
    let mut i = input();
    i.meniscus.meniscus.meniscus.surface_tension_n_m = 1.0;
    i.meniscus.meniscus.meniscus.contact_angle_cosine = 1.0_f64.cos();
    i.contact_line = Some(InertialContactLine {
        friction_pa_s: 0.1,
        advancing_angle_radians: 1.0 + 1e-9,
        receding_angle_radians: 1.0 + 1e-10,
        maximum_angle_departure_radians: 0.3,
    });
    assert!(matches!(
        step_inertial_circular_meniscus(i, controls(), || false),
        Err(InertialCapillaryError::Invalid(_))
    ));
}
