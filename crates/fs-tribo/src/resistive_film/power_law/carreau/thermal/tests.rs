use super::*;
use crate::resistive_film::carreau::{CarreauFilmControls, WallVelocityPorts};
use crate::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

pub(super) const WIDTH: f64 = 0.01;
pub(super) const LENGTH: f64 = 0.02;
pub(super) const GAP: f64 = 2e-4;
const MU: f64 = 0.1;
const ACTIVATION: f64 = 12_000.0;
const GAS_CONSTANT: f64 = 8.314_462_618_153_24;

pub(super) fn near(a: f64, b: f64, abs: f64, rel: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= abs + rel * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

pub(super) fn controls() -> ThermalFilmControls {
    ThermalFilmControls {
        temperature_tolerance_k: 1e-10,
        absolute_energy_tolerance_j: 2e-11,
        relative_energy_tolerance: 2e-10,
        maximum_temperature_change_k: 30.0,
        ..ThermalFilmControls::default()
    }
}

pub(super) fn insulated() -> ThermalFilmBoundary {
    ThermalFilmBoundary {
        prescribed_power_w: 0.0,
        bath_conductance_w_k: 0.0,
        bath_temperature_k: 300.0,
    }
}

pub(super) fn film(
    temperature: f64,
    nonlinear: bool,
    rate_shift: bool,
    supply: f64,
) -> CarreauFilm {
    let gap = GapPort {
        reference_m: GAP,
        closure: vec![1.0, 0.0, 0.0],
    };
    let geometry = super::super::ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 1e-4,
            gap: gap.clone(),
        }],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: WIDTH,
                length_m: LENGTH,
                gap
            };
            2
        ],
        3,
        MU,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 1e-3,
            maximum_pressure_pa: 1e7,
        },
    )
    .unwrap();
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            MU,
            if nonlinear { 0.002 } else { MU },
            0.001,
            0.45,
            2.0,
            [0.0, 1e8],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: ACTIVATION,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        if rate_shift {
            ViscosityShiftMode::VerticalAndRate
        } else {
            ViscosityShiftMode::Vertical
        },
    );
    let passage = LaminarCarreauSlit::new(law, temperature, WIDTH, GAP, LENGTH, 1e6).unwrap();
    let walls = vec![
        WallVelocityPorts {
            lower: vec![0.0, -1.0, 0.0],
            upper: vec![0.0, 0.0, -1.0],
        },
        WallVelocityPorts {
            lower: vec![0.0, 1.0, 0.0],
            upper: vec![0.0, 0.0, 1.0],
        },
    ];
    let mut numerical = CarreauFilmControls::default();
    numerical.pressure.absolute_flow_tolerance_m3_s = 1e-22;
    numerical.pressure.relative_flow_tolerance = 1e-10;
    numerical.pressure.absolute_power_tolerance_w = 1e-18;
    numerical.pressure.relative_power_tolerance = 1e-9;
    numerical.passage.local.relative_constitutive_tolerance = 1e-13;
    numerical.passage.local.relative_quadrature_tolerance = 1e-10;
    CarreauFilm::new(
        &geometry,
        vec![passage; 2],
        &[Some(supply), Some(-supply)],
        walls,
        numerical,
    )
    .unwrap()
}

fn body(
    storage: &TabulatedHeatStorage,
    temperature: f64,
    nonlinear: bool,
    rate_shift: bool,
    supply: f64,
) -> ThermalCarreauFilm<'_> {
    ThermalCarreauFilm::new(
        film(temperature, nonlinear, rate_shift, supply),
        storage,
        &[0.0; 3],
        4.0,
    )
    .unwrap()
}

pub(super) fn viscosity(t: f64, rate: f64, nonlinear: bool, rate_shift: bool) -> f64 {
    let shift = (ACTIVATION / GAS_CONSTANT * (1.0 / t - 1.0 / 300.0)).exp();
    if !nonlinear {
        return MU * shift;
    }
    let reduced = rate * if rate_shift { shift } else { 1.0 };
    shift * (0.002 + (MU - 0.002) * (1.0 + (0.001 * reduced).powi(2)).powf((0.45 - 1.0) / 2.0))
}

fn couette_power(t: f64, relative: f64, rate_shift: bool) -> f64 {
    2.0 * WIDTH
        * LENGTH
        * viscosity(t, relative.abs() / GAP, true, rate_shift)
        * relative
        * relative
        / GAP
}

/// Independent piecewise-affine caloric integral, used only as a test oracle.
fn heat_integral(knots: &[[f64; 2]], initial: f64, final_t: f64) -> f64 {
    let mut sum = 0.0;
    for pair in knots.windows(2) {
        let low = initial.max(pair[0][0]);
        let high = final_t.min(pair[1][0]);
        if low < high {
            let slope = (pair[1][1] - pair[0][1]) / (pair[1][0] - pair[0][0]);
            let left = pair[0][1] + slope * (low - pair[0][0]);
            let right = pair[0][1] + slope * (high - pair[0][0]);
            sum += 0.5 * (left + right) * (high - low);
        }
    }
    sum
}

pub(super) fn independent_root(mut low: f64, mut high: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(low) <= 0.0 && residual(high) >= 0.0);
    for _ in 0..80 {
        let mid = low + 0.5 * (high - low);
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    low + 0.5 * (high - low)
}

#[test]
fn g1_original_capacity_and_midpoint_cooling_match_exact_linear_balance() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.25], [360.0, 0.25]]).unwrap();
    let old = body(&storage, 310.0, true, true, 0.0);
    let boundary = ThermalFilmBoundary {
        prescribed_power_w: 0.05,
        bath_conductance_w_k: 0.04,
        bath_temperature_k: 290.0,
    };
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [CarreauFilmChannelReport::default(); 2];
    let step = old
        .advance_into(
            &[0.0; 3],
            1.0,
            boundary,
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    let expected = 310.0 + (0.05 - 0.04 * (310.0 - 290.0)) / (0.25 + 0.5 * 0.04);
    near(step.next.temperature_k(), expected, 2e-10, 1e-13);
    near(
        step.report.stored_energy_change_j,
        0.25 * (expected - 310.0),
        5e-11,
        1e-12,
    );
    assert_eq!(step.report.dissipated_heat_j, 0.0);
    near(
        step.report.prescribed_heat_j + step.report.bath_heat_j,
        step.report.stored_energy_change_j,
        5e-11,
        1e-12,
    );
    assert_eq!(old.temperature_k(), 310.0);
    assert_eq!(step.next.time_s(), 5.0);
    assert!(std::ptr::eq(step.next.storage(), &storage));
}

#[test]
fn g1_nonlinear_self_heating_crosses_caloric_knots_and_preserves_both_source_shifts() {
    let knots = [[280.0, 0.02], [305.0, 0.08], [330.0, 0.06], [360.0, 0.10]];
    let storage = TabulatedHeatStorage::try_new(&knots).unwrap();
    let boundary = ThermalFilmBoundary {
        prescribed_power_w: 0.001,
        bath_conductance_w_k: 0.0005,
        bath_temperature_k: 295.0,
    };
    let velocity = [0.0, -0.25, 0.25];
    let duration = 24.0;
    let mut results = [0.0; 2];
    for (index, rate_shift) in [false, true].into_iter().enumerate() {
        let old = body(&storage, 300.0, true, rate_shift, 0.0);
        let expected = independent_root(300.0, 330.0, |t| {
            let mid = 0.5 * (300.0 + t);
            heat_integral(&knots, 300.0, t)
                - duration * (couette_power(mid, 0.5, rate_shift) + 0.001 - 0.0005 * (mid - 295.0))
        });
        assert!(expected > 305.0);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [CarreauFilmChannelReport::default(); 2];
        let step = old
            .advance_into(
                &velocity,
                duration,
                boundary,
                controls(),
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .unwrap();
        near(step.next.temperature_k(), expected, 2e-8, 1e-12);
        let mid = step.report.midpoint_temperature_k;
        near(
            step.report.dissipated_heat_j,
            duration * couette_power(mid, 0.5, rate_shift),
            1e-11,
            1e-10,
        );
        near(
            step.report.mechanical_work_j,
            step.report.dissipated_heat_j,
            1e-11,
            1e-9,
        );
        near(step.report.reservoir_work_j, 0.0, 1e-20, 0.0);
        near(
            step.report.thermal_residual_j,
            0.0,
            step.report.thermal_residual_allowance_j,
            0.0,
        );
        near(
            step.report.combined_energy_residual_j,
            step.report.thermal_residual_j - step.report.mechanical_work_residual_j,
            1e-20,
            0.0,
        );
        assert!(step.report.dissipated_heat_j < duration * couette_power(300.0, 0.5, rate_shift));
        for (k, passage) in step.next.film().passages().iter().enumerate() {
            assert_eq!(passage.temperature_k(), step.next.temperature_k());
            assert_eq!(passage.law(), old.film().passages()[k].law());
            near(
                passage.zero_shear_viscosity_pa_s(),
                viscosity(expected, 0.0, true, rate_shift),
                1e-11,
                1e-9,
            );
            near(
                channels[k].dissipated_power_w,
                0.5 * couette_power(mid, 0.5, rate_shift),
                1e-12,
                1e-9,
            );
            near(channels[k].volume_flow_m3_s, 0.0, 1e-21, 0.0);
        }
        near(
            force[1] * velocity[1] + force[2] * velocity[2],
            couette_power(mid, 0.5, rate_shift),
            1e-11,
            1e-9,
        );
        results[index] = step.next.temperature_k();
    }
    assert!((results[0] - results[1]).abs() > 0.01);
}

#[test]
fn g1_pressure_supply_work_heats_the_returning_inventory_and_changes_flow() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.002], [360.0, 0.002]]).unwrap();
    let supply = 2000.0;
    let duration = 80.0;
    let old = body(&storage, 300.0, false, false, supply);
    let power = |t| {
        2.0 * WIDTH * GAP.powi(3) * supply * supply
            / (12.0 * viscosity(t, 0.0, false, false) * LENGTH)
    };
    let expected = independent_root(300.0, 320.0, |t| {
        0.002 * (t - 300.0) - duration * power(0.5 * (300.0 + t))
    });
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [CarreauFilmChannelReport::default(); 2];
    let step = old
        .advance_into(
            &[0.0; 3],
            duration,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    near(step.next.temperature_k(), expected, 2e-8, 1e-12);
    assert_eq!(step.report.mechanical_work_j, 0.0);
    near(
        step.report.reservoir_work_j,
        step.report.dissipated_heat_j,
        1e-12,
        1e-11,
    );
    near(
        step.report.stored_energy_change_j,
        step.report.dissipated_heat_j,
        5e-11,
        1e-10,
    );
    assert!(step.report.dissipated_heat_j > duration * power(300.0));
    near(
        channels[0].volume_flow_m3_s + channels[1].volume_flow_m3_s,
        0.0,
        1e-22,
        0.0,
    );
    let expected_flow = WIDTH * GAP.powi(3) * supply
        / (12.0 * viscosity(step.report.midpoint_temperature_k, 0.0, false, false) * LENGTH);
    near(channels[0].volume_flow_m3_s, -expected_flow, 1e-22, 1e-11);
    near(channels[1].volume_flow_m3_s, expected_flow, 1e-22, 1e-11);
}

#[test]
fn g3_thermal_path_has_second_order_refinement_for_original_nonlinear_shear() {
    let capacity = 0.08;
    let storage = TabulatedHeatStorage::try_new(&[[280.0, capacity], [360.0, capacity]]).unwrap();
    let final_time = 24.0;
    let derivative = |t| couette_power(t, 0.5, true) / capacity;
    let mut reference = 300.0;
    let dt = final_time / 4096.0;
    for _ in 0..4096 {
        let a = derivative(reference);
        let b = derivative(reference + 0.5 * dt * a);
        let c = derivative(reference + 0.5 * dt * b);
        let d = derivative(reference + dt * c);
        reference += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let integrate = |steps| {
        let mut state = body(&storage, 300.0, true, true, 0.0);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [CarreauFilmChannelReport::default(); 2];
        for _ in 0..steps {
            state = state
                .advance_into(
                    &[0.0, -0.25, 0.25],
                    final_time / steps as f64,
                    insulated(),
                    controls(),
                    &mut force,
                    &mut pressure,
                    &mut channels,
                )
                .unwrap()
                .next;
        }
        state.temperature_k()
    };
    let coarse = (integrate(1) - reference).abs();
    let fine = (integrate(2) - reference).abs();
    assert!(
        fine > 1e-8 && coarse / fine > 3.5 && coarse / fine < 4.5,
        "{coarse} / {fine}"
    );
}

#[test]
fn g4_state_domains_fixed_geometry_and_local_work_limits_refuse_atomically() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.002], [350.0, 0.002]]).unwrap();
    let old = body(&storage, 349.9, true, false, 0.0);
    let mut force = [17.0; 3];
    let mut pressure = [19.0];
    let sentinel = CarreauFilmChannelReport {
        volume_flow_m3_s: 23.0,
        ..Default::default()
    };
    let mut channels = [sentinel; 2];
    for (velocity, duration, numerical) in [
        ([1e-6, 0.0, 0.0], 1.0, controls()),
        ([0.0, -0.25, 0.25], 100.0, controls()),
        (
            [0.0, -0.25, 0.25],
            0.001,
            ThermalFilmControls {
                maximum_evaluations: 2,
                ..controls()
            },
        ),
        (
            [0.0, -0.25, 0.25],
            1.0,
            ThermalFilmControls {
                maximum_temperature_change_k: 1e-6,
                ..controls()
            },
        ),
    ] {
        assert!(
            old.advance_into(
                &velocity,
                duration,
                insulated(),
                numerical,
                &mut force,
                &mut pressure,
                &mut channels
            )
            .is_err()
        );
        assert_eq!(force, [17.0; 3]);
        assert_eq!(pressure, [19.0]);
        assert_eq!(channels, [sentinel; 2]);
        assert_eq!(old.time_s(), 4.0);
        assert_eq!(old.temperature_k(), 349.9);
    }
    let original = film(300.0, true, false, 0.0);
    let mut sources = original.passages().to_vec();
    sources[1] =
        LaminarCarreauSlit::new(*sources[1].law(), 301.0, WIDTH, GAP, LENGTH, 1e6).unwrap();
    assert!(
        ThermalCarreauFilm::new(
            original.with_passages(sources).unwrap(),
            &storage,
            &[0.0; 3],
            0.0
        )
        .is_err()
    );
    assert!(
        ThermalCarreauFilm::new(film(355.0, true, false, 0.0), &storage, &[0.0; 3], 0.0).is_err()
    );
}

#[test]
fn g4_stiff_cooling_cannot_cross_equilibrium_and_refinement_recovers() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.25], [360.0, 0.25]]).unwrap();
    let old = body(&storage, 310.0, true, false, 0.0);
    let boundary = ThermalFilmBoundary {
        prescribed_power_w: 0.0,
        bath_conductance_w_k: 1.0,
        bath_temperature_k: 300.0,
    };
    let mut force = [13.0; 3];
    let mut pressure = [17.0];
    let mut channels = [CarreauFilmChannelReport::default(); 2];
    let error = old
        .advance_into(
            &[0.0; 3],
            1.0,
            boundary,
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap_err();
    assert_eq!(
        error,
        FilmError("Carreau thermal step crosses a stationary temperature; refine duration")
    );
    assert_eq!(force, [13.0; 3]);
    assert_eq!(pressure, [17.0]);
    let step = old
        .advance_into(
            &[0.0; 3],
            0.1,
            boundary,
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert!(step.next.temperature_k() > 300.0 && step.next.temperature_k() < 310.0);
}

#[test]
fn g4_late_cancellation_preserves_mechanics_caloric_clock_and_exact_retry() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.08], [360.0, 0.08]]).unwrap();
    let old = body(&storage, 300.0, true, true, 0.0);
    let velocity = [0.0, -0.25, 0.25];
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [CarreauFilmChannelReport::default(); 2];
    let mut calls = 0;
    let clean = old
        .advance_with_check_into(
            &velocity,
            1.0,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let expected = (force, pressure, channels);
    force = [29.0; 3];
    pressure = [31.0];
    channels = [CarreauFilmChannelReport {
        pressure_flow_m3_s: 37.0,
        ..Default::default()
    }; 2];
    let before = channels;
    let mut at = 0;
    let error = old
        .advance_with_check_into(
            &velocity,
            1.0,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
            || {
                at += 1;
                if at == calls {
                    Err(FilmError("thermal final cancellation"))
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
    assert_eq!(error, FilmError("thermal final cancellation"));
    assert_eq!(force, [29.0; 3]);
    assert_eq!(pressure, [31.0]);
    assert_eq!(channels, before);
    assert_eq!(old.time_s(), 4.0);
    assert_eq!(old.temperature_k(), 300.0);
    let retry = old
        .advance_into(
            &velocity,
            1.0,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert_eq!((force, pressure, channels), expected);
    assert_eq!(
        retry.next.temperature_k().to_bits(),
        clean.next.temperature_k().to_bits()
    );
    assert_eq!(
        retry.report.dissipated_heat_j.to_bits(),
        clean.report.dissipated_heat_j.to_bits()
    );
}

#[test]
fn g1_moving_newtonian_squeeze_couples_midgap_pressure_heat_and_volume() {
    let capacity = 2e-5;
    let area = 1e-4;
    let velocity = [2e-4, 0.0, 0.0];
    let duration = 0.4;
    let storage = TabulatedHeatStorage::try_new(&[[280.0, capacity], [360.0, capacity]]).unwrap();
    let old = body(&storage, 300.0, false, false, 0.0);
    let midpoint_gap = GAP - 0.5 * duration * velocity[0];
    let outgoing = 0.5 * area * velocity[0];
    let conductance = |temperature| {
        WIDTH * midpoint_gap.powi(3) / (12.0 * viscosity(temperature, 0.0, false, false) * LENGTH)
    };
    let expected = independent_root(300.0, 330.0, |temperature| {
        capacity * (temperature - 300.0)
            - duration * 2.0 * outgoing * outgoing / conductance(0.5 * (300.0 + temperature))
    });
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = [CarreauFilmChannelReport::default(); 2];
    let step = old
        .advance_moving_into(
            &velocity,
            duration,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    let midpoint_pressure = outgoing / conductance(0.5 * (300.0 + expected));
    near(step.next.temperature_k(), expected, 2e-7, 1e-12);
    near(pressure[0], midpoint_pressure, 1e-6, 1e-10);
    near(force[0], area * midpoint_pressure, 1e-10, 1e-10);
    for channel in channels {
        near(channel.volume_flow_m3_s, outgoing, 1e-22, 1e-11);
    }
    near(
        step.report.dissipated_heat_j,
        duration * 2.0 * outgoing * midpoint_pressure,
        2e-12,
        1e-9,
    );
    near(
        step.report.mechanical_work_j,
        duration * force[0] * velocity[0],
        1e-16,
        1e-11,
    );
    near(
        step.report.stored_energy_change_j,
        step.report.mechanical_work_j,
        step.report.thermal_residual_allowance_j,
        1e-9,
    );
    near(step.report.reservoir_work_j, 0.0, 1e-20, 0.0);
    near(
        step.report.initial_film_cell_volume_m3,
        area * GAP,
        1e-23,
        1e-13,
    );
    near(
        step.report.final_film_cell_volume_m3,
        area * (GAP - duration * velocity[0]),
        1e-23,
        1e-13,
    );
    near(
        step.report.film_cell_volume_change_m3,
        -area * duration * velocity[0],
        1e-23,
        1e-13,
    );
    near(
        step.report.reservoir_outflow_volume_m3,
        area * duration * velocity[0],
        1e-23,
        1e-11,
    );
    near(step.report.volume_balance_residual_m3, 0.0, 1e-21, 0.0);
    assert_eq!(
        step.next.configuration(),
        &[duration * velocity[0], 0.0, 0.0]
    );
    assert_eq!(old.configuration(), &[0.0; 3]);
    assert_eq!(old.temperature_k(), 300.0);
    assert_eq!(step.next.time_s(), 4.0 + duration);
    assert!(std::ptr::eq(step.next.storage(), &storage));
    // The accepted owner is prepared at endpoint temperature and can continue
    // at the actual narrower endpoint, not at the midpoint pressure state.
    step.next
        .film()
        .evaluate_into(
            step.next.configuration(),
            &velocity,
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    let final_gap = GAP - duration * velocity[0];
    let final_conductance =
        WIDTH * final_gap.powi(3) / (12.0 * viscosity(expected, 0.0, false, false) * LENGTH);
    near(pressure[0], outgoing / final_conductance, 2e-6, 1e-10);
}

#[test]
fn g2_nonlinear_squeeze_sliding_heat_retains_original_thermal_laws_and_normal_work() {
    let capacity = 0.02;
    let storage = TabulatedHeatStorage::try_new(&[[280.0, capacity], [360.0, capacity]]).unwrap();
    let velocity = [2e-5, -0.25, 0.25];
    let duration = 2.0;
    let mut temperatures = [0.0; 2];
    for (index, rate_shift) in [false, true].into_iter().enumerate() {
        let old = body(&storage, 300.0, true, rate_shift, 0.0);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [CarreauFilmChannelReport::default(); 2];
        let step = old
            .advance_moving_into(
                &velocity,
                duration,
                insulated(),
                controls(),
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .unwrap();
        assert!(step.next.temperature_k() > 300.0);
        assert!(pressure[0] > 0.0 && force[0] > 0.0);
        near(
            channels[0].volume_flow_m3_s + channels[1].volume_flow_m3_s,
            1e-4 * velocity[0],
            1e-21,
            1e-9,
        );
        near(
            step.report.mechanical_work_j,
            duration * dot(&force, &velocity),
            1e-14,
            1e-10,
        );
        near(
            step.report.dissipated_heat_j,
            duration * channels.iter().map(|c| c.dissipated_power_w).sum::<f64>(),
            1e-14,
            1e-10,
        );
        near(
            step.report.stored_energy_change_j,
            step.report.dissipated_heat_j,
            step.report.thermal_residual_allowance_j,
            1e-9,
        );
        near(step.report.combined_energy_residual_j, 0.0, 1e-9, 0.0);
        near(step.report.volume_balance_residual_m3, 0.0, 1e-20, 0.0);
        for (position, rate) in step.next.configuration().iter().zip(velocity) {
            near(*position, duration * rate, 1e-20, 1e-14);
        }
        for (source, original) in step
            .next
            .film()
            .passages()
            .iter()
            .zip(old.film().passages())
        {
            assert_eq!(source.law(), original.law());
            assert_eq!(source.temperature_k(), step.next.temperature_k());
        }
        temperatures[index] = step.next.temperature_k();
    }
    assert!((temperatures[0] - temperatures[1]).abs() > 1e-4);
}

#[test]
fn g3_moving_squeeze_has_second_order_geometry_and_thermal_refinement() {
    let capacity = 2e-5;
    let area = 1e-4;
    let rate = 2e-4;
    let duration = 0.4;
    let storage = TabulatedHeatStorage::try_new(&[[280.0, capacity], [360.0, capacity]]).unwrap();
    // Independent Newtonian Reynolds reduction with the original Arrhenius
    // viscosity: D=6*mu(T)*L*(A*v)^2/(width*h(t)^3). Time, geometry and
    // temperature evolve together in this separate fourth-order reference.
    let derivative = |time: f64, temperature: f64| {
        let gap = GAP - rate * time;
        6.0 * viscosity(temperature, 0.0, false, false) * LENGTH * (area * rate).powi(2)
            / (WIDTH * gap.powi(3) * capacity)
    };
    let mut reference = 300.0;
    let dt = duration / 4096.0;
    for index in 0..4096 {
        let time = f64::from(index) * dt;
        let a = derivative(time, reference);
        let b = derivative(time + 0.5 * dt, reference + 0.5 * dt * a);
        let c = derivative(time + 0.5 * dt, reference + 0.5 * dt * b);
        let d = derivative(time + dt, reference + dt * c);
        reference += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let integrate = |steps| {
        let mut state = body(&storage, 300.0, false, false, 0.0);
        let mut force = [0.0; 3];
        let mut pressure = [0.0];
        let mut channels = [CarreauFilmChannelReport::default(); 2];
        let mut heat = 0.0;
        let mut displacement = 0.0;
        for _ in 0..steps {
            let step = state
                .advance_moving_into(
                    &[rate, 0.0, 0.0],
                    duration / f64::from(steps),
                    insulated(),
                    controls(),
                    &mut force,
                    &mut pressure,
                    &mut channels,
                )
                .unwrap();
            heat += step.report.dissipated_heat_j;
            displacement += step.report.reservoir_outflow_volume_m3;
            state = step.next;
        }
        near(
            capacity * (state.temperature_k() - 300.0),
            heat,
            1e-10,
            1e-9,
        );
        near(displacement, area * rate * duration, 1e-22, 1e-10);
        near(state.configuration()[0], rate * duration, 1e-20, 1e-13);
        (state.temperature_k() - reference).abs()
    };
    let errors = [integrate(1), integrate(2), integrate(4)];
    assert!(
        errors[2] > 1e-6
            && errors[0] / errors[1] > 3.0
            && errors[1] / errors[2] > 3.6
            && errors[1] / errors[2] < 4.5,
        "{errors:?}"
    );
}

#[test]
fn g4_moving_thermal_endpoint_admission_budget_and_cancellation_are_atomic() {
    let storage = TabulatedHeatStorage::try_new(&[[280.0, 0.08], [360.0, 0.08]]).unwrap();
    let old = body(&storage, 300.0, false, false, 0.0);
    let sentinel = CarreauFilmChannelReport {
        volume_flow_m3_s: 23.0,
        ..Default::default()
    };
    let mut force = [17.0; 3];
    let mut pressure = [19.0];
    let mut channels = [sentinel; 2];
    for (duration, numerical) in [
        (1.0, controls()),  // midpoint valid; endpoint closes the volume cell
        (0.98, controls()), // endpoint gap valid; endpoint pressure exceeds support
        (
            0.1,
            ThermalFilmControls {
                maximum_evaluations: 4,
                ..controls()
            },
        ),
    ] {
        assert!(
            old.advance_moving_into(
                &[2e-4, 0.0, 0.0],
                duration,
                insulated(),
                numerical,
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .is_err()
        );
        assert_eq!(force, [17.0; 3]);
        assert_eq!(pressure, [19.0]);
        assert_eq!(channels, [sentinel; 2]);
        assert_eq!(old.configuration(), &[0.0; 3]);
        assert_eq!(old.temperature_k(), 300.0);
        assert_eq!(old.time_s(), 4.0);
    }
    // A wide mid-gap cannot hide an unsupported starting squeeze state.
    let narrow = ThermalCarreauFilm::new(
        film(300.0, false, false, 0.0),
        &storage,
        &[GAP - 5e-6, 0.0, 0.0],
        4.0,
    )
    .unwrap();
    assert!(
        narrow
            .advance_moving_into(
                &[-2e-4, 0.0, 0.0],
                0.4,
                insulated(),
                controls(),
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .is_err()
    );
    assert_eq!(force, [17.0; 3]);
    assert_eq!(pressure, [19.0]);
    assert_eq!(channels, [sentinel; 2]);

    let velocity = [2e-5, -0.25, 0.25];
    let mut calls = 0;
    let clean = old
        .advance_moving_with_check_into(
            &velocity,
            0.2,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let expected = (force, pressure, channels);
    force = [17.0; 3];
    pressure = [19.0];
    channels = [sentinel; 2];
    let mut at = 0;
    let error = old
        .advance_moving_with_check_into(
            &velocity,
            0.2,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
            || {
                at += 1;
                if at == calls {
                    Err(FilmError("moving thermal final cancellation"))
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
    assert_eq!(error, FilmError("moving thermal final cancellation"));
    assert_eq!(
        (force, pressure, channels),
        ([17.0; 3], [19.0], [sentinel; 2])
    );
    assert_eq!(old.configuration(), &[0.0; 3]);
    assert_eq!(old.time_s(), 4.0);
    let retry = old
        .advance_moving_into(
            &velocity,
            0.2,
            insulated(),
            controls(),
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert_eq!((force, pressure, channels), expected);
    assert_eq!(retry.next.configuration(), clean.next.configuration());
    assert_eq!(
        retry.next.temperature_k().to_bits(),
        clean.next.temperature_k().to_bits()
    );
    assert_eq!(
        retry.report.dissipated_heat_j.to_bits(),
        clean.report.dissipated_heat_j.to_bits()
    );
}
