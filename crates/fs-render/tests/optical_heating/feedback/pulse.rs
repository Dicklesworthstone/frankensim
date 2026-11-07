use super::*;
use fs_render::volumes::heating::feedback::nonlinear::{
    TemperatureTwoPhotonTable, ThermalTwoPhotonCell, TwoPhotonFeedback,
};
use fs_render::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};
use fs_thermal::enthalpy::DepositionObserver;

fn segment(start: f64, end: f64, a: f64, b: f64) -> SpectralPowerSegment<1> {
    SpectralPowerSegment {
        start_s: start,
        end_s: end,
        start_power_w: [a],
        end_power_w: [b],
    }
}

fn profile(segments: &[SpectralPowerSegment<1>]) -> SpectralPowerSchedule<1> {
    SpectralPowerSchedule::try_new(segments, 8).unwrap()
}

fn beta(bounds: (f64, f64), values: (f64, f64), ceiling: f64) -> TemperatureTwoPhotonTable {
    TemperatureTwoPhotonTable::try_new(
        ContentHash([107; 32]),
        &[(bounds.0, values.0), (bounds.1, values.1)],
        (400.0, 800.0),
        ceiling,
        2,
    )
    .unwrap()
}

#[allow(clippy::too_many_arguments)]
fn pulsed(
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'_>,
    peak: f64,
    linear: &[TemperatureAbsorptionTable],
    nonlinear: &[TemperatureTwoPhotonTable],
    schedule: SpectralPowerSchedule<1>,
    work: OpticalFeedbackBudget,
) -> Result<TemperatureOpticalHeating<1>, OpticalFeedbackError> {
    let density = vec![1.0; linear.len()];
    let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(1.0).unwrap());
    let cells: Vec<_> = linear
        .iter()
        .enumerate()
        .map(|(i, _)| binding(i, (i as u64 + 1) * 11, i))
        .collect();
    let nonlinear_cells: Vec<_> = linear
        .iter()
        .enumerate()
        .map(|(i, _)| ThermalTwoPhotonCell {
            grid_cell: i,
            two_photon_table: i,
        })
        .collect();
    TemperatureOpticalHeating::try_new_two_photon(
        cx,
        &medium,
        beam(peak, linear.len()),
        linear,
        &cells,
        TwoPhotonFeedback {
            beam_area_m2: 1.0,
            tables: nonlinear,
            cells: &nonlinear_cells,
        },
        thermal,
        work,
        1.0e-11,
        || Ok(()),
    )?
    .with_power_schedule(thermal, schedule, 8, || Ok(()))
}

#[test]
fn g1_spectral_pulses_preserve_exact_fluence_native_events_dark_gaps_and_restart() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![TimedEvent {
                time_s: 3.5,
                event: Event::Boundary(CellId(11), boundary(2.0)),
            }],
        );
        let density = [1.0];
        let medium = grid(cx, &density, MassAbsorptionSpectrum::constant(1.0).unwrap());
        let source = SpectralPowerSchedule::try_new(
            &[
                SpectralPowerSegment {
                    start_s: 1.0,
                    end_s: 2.0,
                    start_power_w: [0.0, 6.0],
                    end_power_w: [4.0, 2.0],
                },
                SpectralPowerSegment {
                    start_s: 3.0,
                    end_s: 4.0,
                    start_power_w: [3.0, 0.0],
                    end_power_w: [3.0, 8.0],
                },
            ],
            2,
        )
        .unwrap();
        let optics = TemperatureOpticalHeating::try_new(
            cx,
            &medium,
            IncidentBeam {
                origin: Point3::new(-1.0, 0.5, 0.5),
                unit_direction: Vec3::new(1.0, 0.0, 0.0),
                max_distance_m: 3.0,
                wavelengths_nm: [500.0, 700.0],
                incident_power_w: [10.0, 20.0],
            },
            &[table(300.0, 400.0, 0.2, 0.2)],
            &[binding(0, 11, 0)],
            &original,
            budget(),
            1.0e-11,
            || Ok(()),
        )
        .unwrap()
        .with_power_schedule(&original, source, 2, || Ok(()))
        .unwrap();
        assert!(optics.requires_step_interval());
        assert_eq!(optics.power().incident_w, [0.0; 2]);
        let (whole, accepted, audit) = original
            .advance_observed(&optics, 5.0, config(10.0), || Ok(()))
            .unwrap();
        assert_eq!(audit.substeps, 6);
        assert_eq!(audit.events, 1);
        assert_eq!(accepted.source_evaluations(), audit.substeps);
        close(accepted.energy().incident_j[0], 5.0);
        close(accepted.energy().incident_j[1], 8.0);
        close(
            accepted.energy().thermal_input_j,
            13.0 * -(-0.2_f64).exp_m1(),
        );
        close(
            audit.external_heat_j,
            accepted.energy().thermal_input_j + 3.0,
        );
        close(whole.energy_j(), audit.external_heat_j);
        assert_eq!(accepted.power().incident_w, [0.0; 2]);
        let (middle, mid_source, _) = original
            .advance_observed(&optics, 3.0, config(10.0), || Ok(()))
            .unwrap();
        let (split, split_source, _) = middle
            .advance_observed(&mid_source, 5.0, config(10.0), || Ok(()))
            .unwrap();
        assert_eq!(split.cells(), whole.cells());
        assert_eq!(split_source.energy(), accepted.energy());
    });
}

#[test]
fn g1_two_photon_ramps_integrate_actual_intensity_and_preserve_time_subdivision() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let optics = pulsed(
            cx,
            &original,
            4.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta((300.0, 400.0), (0.5, 0.5), 4.0)],
            profile(&[segment(1.0, 3.0, 0.0, 4.0)]),
            budget(),
        )
        .unwrap();
        // Independent integral: P(t) traverses 0..4 W over two seconds,
        // P_out=P/(1+0.5P), and integral P_out dt=2*(2-ln(3)).
        let transmitted = 2.0 * (2.0 - 3.0_f64.ln());
        let absorbed = 4.0 - transmitted;
        // Applying beta to the mean incident 2 W would incorrectly give 2 J.
        assert!(absorbed > 2.1);
        for dt in [10.0, 0.25, 0.05] {
            let (next, accepted, audit) = original
                .advance_observed(&optics, 4.0, config(dt), || Ok(()))
                .unwrap();
            close(accepted.energy().incident_j[0], 4.0);
            close(accepted.energy().transmitted_j[0], transmitted);
            close(accepted.energy().absorbed_j[0], absorbed);
            close(audit.external_heat_j, absorbed);
            close(
                next.temperature_k(CellId(11)).unwrap(),
                300.0 + absorbed / 4.0,
            );
        }
        let (next, accepted, report) = optics
            .advance_refined(
                &original,
                4.0,
                FeedbackRefinement {
                    enthalpy_absolute_tolerance_j_kg: 1.0e-9,
                    temperature_absolute_tolerance_k: 1.0e-9,
                    optical_energy_absolute_tolerance_j: 1.0e-9,
                    ..refinement()
                },
                || Ok(()),
            )
            .unwrap();
        assert!(report.roundoff_limited);
        assert_eq!(report.observed_order, None);
        close(next.energy_j(), absorbed);
        close(accepted.energy().thermal_input_j, absorbed);
    });
}

#[test]
fn g3_affine_pulses_resolve_heterogeneous_nonlinear_shielding_cell_by_cell() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let optics = pulsed(
            cx,
            &original,
            9.0,
            &[table(300.0, 400.0, 0.1, 0.1), table(300.0, 400.0, 0.7, 0.7)],
            &[
                beta((300.0, 400.0), (0.15, 0.15), 10.0),
                beta((300.0, 400.0), (0.4, 0.4), 10.0),
            ],
            profile(&[segment(0.0, 1.0, 1.0, 9.0), segment(1.0, 2.0, 9.0, 1.0)]),
            budget(),
        )
        .unwrap();
        let (next, accepted, _) = original
            .advance_observed(&optics, 2.0, config(2.0), || Ok(()))
            .unwrap();
        let mut reference = [0.0; 3];
        let intervals = 4096;
        for sample in 0..=intervals {
            let mut power = 1.0 + 8.0 * f64::from(sample) / f64::from(intervals);
            let weight = if sample == 0 || sample == intervals {
                1.0
            } else if sample % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for (i, (alpha, beta)) in [(0.1_f64, 0.15), (0.7, 0.4)].into_iter().enumerate() {
                let remaining =
                    power * (-alpha).exp() / (1.0 + beta * power * -(-alpha).exp_m1() / alpha);
                reference[i] += weight * (power - remaining);
                power = remaining;
            }
            reference[2] += weight * power;
        }
        for value in &mut reference {
            *value *= 2.0 / (3.0 * f64::from(intervals));
        }
        close(next.cells()[&CellId(11)].specific_h_j_kg, reference[0]);
        close(next.cells()[&CellId(22)].specific_h_j_kg, reference[1]);
        close(accepted.energy().transmitted_j[0], reference[2]);
        close(accepted.energy().incident_j[0], 10.0);
        close(
            accepted.energy().thermal_input_j + accepted.energy().transmitted_j[0],
            10.0,
        );
    });
}

#[test]
fn g1_pulsed_temperature_feedback_converges_against_independent_runge_kutta() {
    with_cx(|cx| {
        let original = thermal(vec![cell(11, 1.0, 1.0, 0.0)], &SENSIBLE, vec![]);
        let optics = pulsed(
            cx,
            &original,
            6.0,
            &[table(300.0, 400.0, 0.0, 0.0)],
            &[beta((300.0, 400.0), (0.2, 10.2), 10.0)],
            profile(&[segment(0.0, 1.0, 2.0, 6.0)]),
            OpticalFeedbackBudget {
                max_source_evaluations: 4000,
                ..budget()
            },
        )
        .unwrap();
        let rhs = |time: f64, rise: f64| {
            let power = 2.0 + 4.0 * time;
            let beta = 0.2 + 0.1 * rise;
            (power - power / (1.0 + beta * power)) / 4.0
        };
        let mut rise = 0.0;
        let dt = 1.0 / 8192.0;
        for step in 0..8192 {
            let time = f64::from(step) * dt;
            let k1 = rhs(time, rise);
            let k2 = rhs(time + 0.5 * dt, rise + 0.5 * dt * k1);
            let k3 = rhs(time + 0.5 * dt, rise + 0.5 * dt * k2);
            let k4 = rhs(time + dt, rise + dt * k3);
            rise += dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        }
        let mut errors = Vec::new();
        for dt in [0.125, 0.0625, 0.03125] {
            let (next, _, _) = original
                .advance_observed(&optics, 1.0, config(dt), || Ok(()))
                .unwrap();
            errors.push((next.temperature_k(CellId(11)).unwrap() - 300.0 - rise).abs());
        }
        assert!((1.8..2.2).contains(&(errors[0] / errors[1])));
        assert!((1.8..2.2).contains(&(errors[1] / errors[2])));
        let (next, accepted, convergence) = optics
            .advance_refined(
                &original,
                1.0,
                FeedbackRefinement {
                    initial_max_dt_s: 0.25,
                    minimum_max_dt_s: 0.0001,
                    maximum_halvings: 10,
                    temperature_absolute_tolerance_k: 1.0e-4,
                    enthalpy_absolute_tolerance_j_kg: 4.0e-4,
                    optical_energy_absolute_tolerance_j: 4.0e-4,
                    ..refinement()
                },
                || Ok(()),
            )
            .unwrap();
        assert!((next.temperature_k(CellId(11)).unwrap() - 300.0 - rise).abs() < 1.1e-4);
        assert!(convergence.normalized_difference <= 1.0);
        close(accepted.energy().incident_j[0], 4.0);
    });
}

#[test]
fn g4_jump_endpoint_domains_work_limits_and_cancellation_are_transactional() {
    with_cx(|cx| {
        let original = thermal(
            vec![cell(11, 1.0, 1.0, 0.0), cell(22, 1.0, 1.0, 0.0)],
            &SENSIBLE,
            vec![],
        );
        let linear = [table(300.0, 400.0, 0.0, 0.0), table(300.0, 400.0, 0.0, 0.0)];
        let nonlinear = [
            beta((300.0, 300.2), (1.0, 0.0), 10.0),
            beta((300.0, 400.0), (0.0, 0.0), 1.0),
        ];
        let source = profile(&[segment(0.0, 1.0, 1.0, 1.0), segment(1.0, 2.0, 8.0, 8.0)]);
        let optics = pulsed(cx, &original, 8.0, &linear, &nonlinear, source, budget()).unwrap();
        // The low pulse heats/clears the front. At the final requested endpoint,
        // the new 8 W level overexposes the rear source even without another step.
        assert!(matches!(
            original.advance_observed(&optics, 1.0, config(1.0), || Ok(())),
            Err(OpticalFeedbackError::IrradianceDomain { grid_cell: 1, .. })
        ));
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(optics.energy().absorbed_j, [0.0]);
        // A pulse that actually ends here is valid: future dark power does not
        // inherit the original beam ceiling as a fictitious exposure.
        let ending = pulsed(
            cx,
            &original,
            8.0,
            &linear,
            &nonlinear,
            profile(&[segment(0.0, 1.0, 1.0, 1.0)]),
            budget(),
        )
        .unwrap();
        let (next, _, _) = original
            .advance_observed(&ending, 1.0, config(1.0), || Ok(()))
            .unwrap();
        close(next.temperature_k(CellId(11)).unwrap(), 300.125);
        let bounded = pulsed(
            cx,
            &original,
            8.0,
            &linear,
            &nonlinear,
            profile(&[segment(0.0, 1.0, 1.0, 1.0)]),
            OpticalFeedbackBudget {
                max_source_evaluations: 1,
                ..budget()
            },
        )
        .unwrap();
        assert!(matches!(
            original.advance_observed(&bounded, 1.0, config(1.0), || Ok(())),
            Err(OpticalFeedbackError::Budget)
        ));
        assert!(matches!(
            original.advance_observed(&ending, 1.0, config(1.0), || Err(
                DepositionError::Cancelled
            )),
            Err(OpticalFeedbackError::Thermal(DepositionError::Cancelled))
        ));
        assert!(
            pulsed(
                cx,
                &original,
                8.0,
                &linear,
                &nonlinear,
                profile(&[segment(0.0, 1.0, 9.0, 9.0)]),
                budget()
            )
            .is_err()
        );
        assert!(
            ending
                .clone()
                .with_power_schedule(
                    &original,
                    profile(&[segment(0.0, 1.0, 1.0, 1.0)]),
                    0,
                    || Ok(())
                )
                .is_err()
        );
        let mut candidate = optics.clone();
        assert!(
            candidate
                .step_sources(&original, 0.0, 1.5, 1.5, &mut || Ok(()))
                .is_err()
        );
        assert_eq!(original.time_s(), 0.0);
        assert_eq!(ending.energy().thermal_input_j, 0.0);
    });
}
