//! G3: original per-line Maxwell loads, one caloric solve, and atomic feedback.

mod incident;

use super::*;
use fs_couple::layered_wave_heating::{
    BroadbandThermoOpticalLimits, BroadbandThermoOpticalWaveHeating, ThermoOpticalSpectralLine,
};

const LIMITS: BroadbandThermoOpticalLimits = BroadbandThermoOpticalLimits {
    maximum_lines: 8,
    maximum_cell_line_pairs: 64,
    maximum_cell_factor_pairs: 128,
};

fn line(
    wavelength_m: f64,
    incidence_cosine: f64,
    polarization: Polarization,
) -> ThermoOpticalSpectralLine {
    ThermoOpticalSpectralLine {
        angular_frequency_rad_s: core::f64::consts::TAU * LIGHT_SPEED / wavelength_m,
        incidence_cosine,
        polarization,
    }
}

fn relative(a: f64, b: f64, tolerance: f64) {
    near(a, b, tolerance * a.abs().max(b.abs()).max(1e-30));
}

#[test]
fn one_spectral_line_exactly_matches_the_original_temperature_refreshed_step() {
    let card = infrared_card(true, 360.0, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, true, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let native = body(&[&card], &[4e-6], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[290.0, 325.0]);
    let sources = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let spectral_line = line(11e-6, 0.9, Polarization::Tm);
    let mono = ThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        spectral_line.angular_frequency_rad_s,
        spectral_line.incidence_cosine,
        spectral_line.polarization,
        &sources,
        || false,
    )
    .unwrap();
    let broadband = BroadbandThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        &[spectral_line],
        &sources,
        LIMITS,
        || false,
    )
    .unwrap();
    let factors = [
        [C64::new(1e6, 0.0), C64::new(2e5, -1e5)],
        [C64::new(0.0, 3e5), C64::new(-1e5, 2e5)],
    ];
    let drive = LayeredWaveDrive::Coherence {
        factors: &factors,
        budget: COHERENCE,
    };
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let controls = LayeredThermalSolveControls::default();
    let original = mono
        .advance_to(
            &initial,
            2e-10,
            drive,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    let result = broadband
        .advance_to(
            &initial,
            2e-10,
            &[drive],
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(result.line_power, vec![original.step.optical_power]);
    assert_eq!(result.step.optical_power, original.step.optical_power);
    assert_eq!(result.step.optical, original.step.optical);
    assert_eq!(result.step.cell_energy_j, original.step.cell_energy_j);
    assert_eq!(result.step.thermal.cells(), original.step.thermal.cells());
    assert_eq!(
        result.applied_optical_states,
        original.applied_optical_states
    );
    assert_eq!(
        result.endpoint_optical_states,
        original.endpoint_optical_states
    );
}

#[test]
fn different_wavelengths_add_cell_watts_once_and_ignore_only_cross_line_phase() {
    let card = infrared_card(false, 360.0, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, false, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let width = 1.5e-6;
    let native = body(&[&card], &[width], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[290.0, 330.0]);
    let before = initial.cells().clone();
    let temperatures: Vec<_> = bound
        .cell_ids()
        .iter()
        .map(|&id| initial.temperature_k(id).unwrap())
        .collect();
    let refreshed: Vec<_> = temperatures
        .iter()
        .map(|&t| source.at_query_point(&card, &at(t, 100_000.0)).unwrap())
        .collect();
    let lines = [
        line(8.5e-6, 0.92, Polarization::Te),
        line(13e-6, 0.7, Polarization::Tm),
    ];
    let inputs = [
        [C64::new(8e5, 0.0), C64::new(4e5, 2e5)],
        [C64::new(1e6, 2e5), C64::new(-2e5, 4e5)],
    ];
    let drives = inputs.map(LayeredWaveDrive::Coherent);
    let sources = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    let owner = BroadbandThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        &lines,
        &sources,
        LIMITS,
        || false,
    )
    .unwrap();
    let mut expected_power = [0.0; 2];
    let mut original_line_power = Vec::new();
    for (line, drive) in lines.iter().zip(drives) {
        let responses: Vec<_> = refreshed
            .iter()
            .map(|state| state.wave_response(line.angular_frequency_rad_s).unwrap())
            .collect();
        let reference = LayeredWaveHeating::try_new_cellwise(
            &bound,
            &outside,
            DielectricResponse::Lossless(&outside),
            line.angular_frequency_rad_s,
            &responses,
            || false,
        )
        .unwrap();
        let mut powers = [0.0; 2];
        let power = reference
            .prepare(line.incidence_cosine, line.polarization, || false)
            .unwrap()
            .cell_power_into(drive, &mut powers, || false)
            .unwrap();
        original_line_power.push(power);
        for (sum, value) in expected_power.iter_mut().zip(powers) {
            *sum += value;
        }
    }
    let dt = 2e-10;
    let controls = LayeredThermalSolveControls {
        cell_tolerance_j: 1e-11,
        total_tolerance_j: 1e-11,
        ..LayeredThermalSolveControls::default()
    };
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let result = owner
        .advance_to(
            &initial,
            dt,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(result.line_power, original_line_power);
    assert_eq!(result.applied_optical_states, refreshed);
    // Independent constant-capacity two-cell backward-Euler solution. Each
    // line must contribute during the same dt, while conduction runs once.
    let capacity = RHO * AREA * (width / 2.0) * CP;
    let g_dt = AREA * K / (width / 2.0) * dt;
    let mean = 0.5 * (temperatures[0] + temperatures[1])
        + 0.5 * dt * (expected_power[0] + expected_power[1]) / capacity;
    let contrast = (0.5 * capacity * (temperatures[0] - temperatures[1])
        + 0.5 * dt * (expected_power[0] - expected_power[1]))
        / (capacity + 2.0 * g_dt);
    for (i, expected) in [mean + contrast, mean - contrast].into_iter().enumerate() {
        let id = bound.cell_ids()[i];
        let t = result.step.thermal.temperature_k(id).unwrap();
        near(t, expected, 3e-9);
        relative(
            result.step.cell_energy_j[i].1,
            expected_power[i] * dt,
            2e-13,
        );
        assert_eq!(result.step.cell_energy_j[i].0, id);
        assert_eq!(
            result.endpoint_optical_states[i],
            source.at_query_point(&card, &at(t, 100_000.0)).unwrap()
        );
    }
    near(
        result.step.heat_report.energy.stored_energy_change_j,
        result.step.optical.absorbed_j,
        3e-11,
    );
    near(result.step.optical.residual_j, 0.0, 3e-11);
    for side in 0..2 {
        relative(
            result.step.optical.incoming_j[side],
            original_line_power
                .iter()
                .map(|p| p.incoming_w_m2[side])
                .sum::<f64>()
                * AREA
                * dt,
            2e-13,
        );
        relative(
            result.step.optical.outgoing_j[side],
            original_line_power
                .iter()
                .map(|p| p.outgoing_w_m2[side])
                .sum::<f64>()
                * AREA
                * dt,
            2e-13,
        );
    }
    // Rotate both faces of only the second line: its global phase cannot
    // interfere with another wavelength in this declared incoherent model.
    let phase = C64::new(0.0, 1.0);
    let rotated = [
        drives[0],
        LayeredWaveDrive::Coherent(inputs[1].map(|e| e * phase)),
    ];
    let phase_result = owner
        .advance_to(
            &initial,
            dt,
            &rotated,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    for (a, b) in result
        .step
        .cell_energy_j
        .iter()
        .zip(&phase_result.step.cell_energy_j)
    {
        relative(a.1, b.1, 2e-12);
    }
    // Relative phase between faces of ONE line remains coherent and changes
    // the finite film's standing-wave loss; incoherence is only between lines.
    let changed = [
        drives[0],
        LayeredWaveDrive::Coherent([inputs[1][0], -inputs[1][1]]),
    ];
    let changed_result = owner
        .advance_to(
            &initial,
            dt,
            &changed,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert!(
        result
            .step
            .cell_energy_j
            .iter()
            .zip(&changed_result.step.cell_energy_j)
            .any(|(a, b)| (a.1 - b.1).abs() > 1e-5 * a.1.abs())
    );
    let next = owner
        .advance_to(
            &result.step.thermal,
            2.0 * dt,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || false,
        )
        .unwrap();
    assert_eq!(next.applied_optical_states, result.endpoint_optical_states);
    assert_ne!(next.applied_optical_states, result.applied_optical_states);
    assert_eq!(initial.cells(), &before);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn spectrum_and_factor_limits_and_early_mid_late_cancellation_are_atomic() {
    let card = infrared_card(true, 360.0, false);
    let source =
        SpectralThermoOpticalSource::try_new(infrared_source(&card, true, true), MAX_GAP_M)
            .unwrap();
    let outside = exterior();
    let native = body(&[&card], &[4e-6], &[2], AREA);
    let laws = native.thermal_laws(false).unwrap();
    let bound = native.bind(&laws, ids()).unwrap();
    let initial = initial_at(&bound, &[300.0, 300.0]);
    let before = initial.cells().clone();
    let lines = [
        line(8.5e-6, 0.9, Polarization::Te),
        line(13e-6, 0.8, Polarization::Tm),
    ];
    let sources = [ThermoOpticalLayer {
        source: &source,
        card: &card,
    }];
    for limits in [
        BroadbandThermoOpticalLimits {
            maximum_lines: 1,
            ..LIMITS
        },
        BroadbandThermoOpticalLimits {
            maximum_cell_line_pairs: 3,
            ..LIMITS
        },
    ] {
        assert!(
            BroadbandThermoOpticalWaveHeating::try_new(
                &bound,
                &outside,
                DielectricResponse::Lossless(&outside),
                &lines,
                &sources,
                limits,
                || false,
            )
            .is_err()
        );
    }
    assert!(
        BroadbandThermoOpticalWaveHeating::try_new(
            &bound,
            &outside,
            DielectricResponse::Lossless(&outside),
            &[],
            &sources,
            LIMITS,
            || false,
        )
        .is_err()
    );
    let owner = BroadbandThermoOpticalWaveHeating::try_new(
        &bound,
        &outside,
        DielectricResponse::Lossless(&outside),
        &lines,
        &sources,
        BroadbandThermoOpticalLimits {
            maximum_cell_factor_pairs: 4,
            ..LIMITS
        },
        || false,
    )
    .unwrap();
    let input = [C64::new(1e6, 0.0), C64::new(2e5, -1e5)];
    let drives = [LayeredWaveDrive::Coherent(input); 2];
    let factors = [input; 2];
    let overflow = [
        drives[0],
        LayeredWaveDrive::Coherence {
            factors: &factors,
            budget: COHERENCE,
        },
    ];
    let own_limit = [
        drives[0],
        LayeredWaveDrive::Coherence {
            factors: &factors[..1],
            budget: IlluminationBudget {
                maximum_modes: 0,
                maximum_layer_mode_pairs: 2,
            },
        },
    ];
    let mut workspace = LayeredThermalWorkspace::new(2).unwrap();
    let controls = LayeredThermalSolveControls::default();
    for invalid in [&drives[..1], &overflow[..], &own_limit[..]] {
        assert!(
            owner
                .advance_to(
                    &initial,
                    2e-10,
                    invalid,
                    conduction(),
                    controls,
                    &mut workspace,
                    || false,
                )
                .is_err()
        );
    }
    let mut total_polls = 0;
    owner
        .advance_to(
            &initial,
            2e-10,
            &drives,
            conduction(),
            controls,
            &mut workspace,
            || {
                total_polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            owner.advance_to(
                &initial,
                2e-10,
                &drives,
                conduction(),
                controls,
                &mut workspace,
                || {
                    polls += 1;
                    polls == cutoff
                },
            ),
            Err(LayeredWaveHeatingError::Cancelled)
        ));
        assert_eq!(initial.cells(), &before);
        assert_eq!(initial.time_s(), 0.0);
    }
}

#[test]
fn every_line_keeps_its_original_band_and_endpoint_temperature_and_gap_support() {
    for (dielectric, upper, widen_gap) in [(true, 300.01, false), (false, 360.0, true)] {
        let card = infrared_card(dielectric, upper, widen_gap);
        let source = SpectralThermoOpticalSource::try_new(
            infrared_source(&card, dielectric, true),
            MAX_GAP_M,
        )
        .unwrap();
        let outside = exterior();
        let native = body(&[&card], &[2e-6], &[1], AREA);
        let laws = native.thermal_laws(false).unwrap();
        let bound = native.bind(&laws, ids()).unwrap();
        let initial = initial_at(&bound, &[300.0]);
        let before = initial.cells().clone();
        let lines = [
            line(11e-6, 1.0, Polarization::Te),
            line(8.5e-6, 1.0, Polarization::Te),
        ];
        let sources = [ThermoOpticalLayer {
            source: &source,
            card: &card,
        }];
        let out_of_band = [lines[0], line(550e-9, 1.0, Polarization::Te)];
        assert!(
            BroadbandThermoOpticalWaveHeating::try_new(
                &bound,
                &outside,
                DielectricResponse::Lossless(&outside),
                &out_of_band,
                &sources,
                LIMITS,
                || false,
            )
            .is_err()
        );
        let owner = BroadbandThermoOpticalWaveHeating::try_new(
            &bound,
            &outside,
            DielectricResponse::Lossless(&outside),
            &lines,
            &sources,
            LIMITS,
            || false,
        )
        .unwrap();
        let drive = LayeredWaveDrive::Coherent([C64::new(1.5e6, 0.0), C64::ZERO]);
        let mut workspace = LayeredThermalWorkspace::new(1).unwrap();
        let controls = LayeredThermalSolveControls::default();
        // Even the first line alone admits a native thermal candidate beyond
        // the optical endpoint support. The second nonnegative load cannot
        // repair that source failure or authorize publishing a partial state.
        let responses = [source.wave_response(IR_OMEGA).unwrap()];
        let frozen = LayeredWaveHeating::try_new(
            &bound,
            &outside,
            DielectricResponse::Lossless(&outside),
            IR_OMEGA,
            &responses,
            || false,
        )
        .unwrap();
        let candidate = frozen
            .prepare(1.0, Polarization::Te, || false)
            .unwrap()
            .advance_to(
                &initial,
                1e-9,
                drive,
                conduction(),
                controls,
                &mut workspace,
                || false,
            )
            .unwrap();
        let t = candidate
            .thermal
            .temperature_k(bound.cell_ids()[0])
            .unwrap();
        assert!(t > 300.05 && t < 360.0);
        let result = owner.advance_to(
            &initial,
            1e-9,
            &[drive; 2],
            conduction(),
            controls,
            &mut workspace,
            || false,
        );
        if widen_gap {
            assert!(matches!(
                result,
                Err(LayeredWaveHeatingError::Electromagnetic(_))
            ));
        } else {
            assert!(matches!(result, Err(LayeredWaveHeatingError::Material(_))));
        }
        assert_eq!(initial.cells(), &before);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(source.resolved().state_coordinate("T"), Some(300.0));
    }
}
