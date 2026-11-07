//! G1/G3/G4: physical-time forcing, source phase inventory and atomic refusal.
use super::*;
use crate::material::compression::MAX_COMPRESSION_SUBSTEPS;
use crate::material::compression::thermal::adaptive::{
    AdaptiveMaterialThermalCompressionError as AdaptiveError,
    AdaptiveMaterialThermalCompressionSolve,
};
use crate::material::compression::thermal::schedule::{
    ConductanceKnot, ConductanceSchedule, HeatBathKnot, HeatBathSchedule,
};

fn bath_knot(time: f64, temperature: f64, conductance: f64) -> HeatBathKnot {
    HeatBathKnot {
        elapsed_time: Time::new(time),
        bath: HeatBath::new(
            ThermalConductance::new(conductance),
            Temperature::new(temperature),
        )
        .unwrap(),
    }
}

fn conductance_knot(time: f64, conductance: f64) -> ConductanceKnot {
    ConductanceKnot {
        elapsed_time: Time::new(time),
        conductance: ThermalConductance::new(conductance),
    }
}

fn adaptation(tolerance: f64) -> AdaptiveMaterialThermalCompressionSolve {
    AdaptiveMaterialThermalCompressionSolve {
        volume: path(1).volume,
        absolute_temperature_tolerance_k: tolerance,
        absolute_amount_tolerance_mol: 1e-9,
        absolute_work_tolerance_j: tolerance,
        absolute_heat_tolerance_j: tolerance,
        absolute_body_enthalpy_tolerance_j: tolerance,
        relative_tolerance: 0.0,
        initial_step_s: 0.125,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.125,
        maximum_intervals: 4096,
        maximum_evaluations: 12288,
    }
}

fn bath_program() -> [HeatBathKnot; 4] {
    [
        bath_knot(0.0, 300.0, 0.1),
        bath_knot(0.03125, 800.0, 1.5),
        bath_knot(0.078125, 250.0, 0.2),
        bath_knot(0.125, 650.0, 0.8),
    ]
}

// Independent coupled ODE integration: no native split, kinetics, source
// inversion, schedule sampling or adaptive controller enters these references.
fn rk4<const N: usize>(
    times: &[f64],
    mut state: [f64; N],
    rate: impl Fn(f64, [f64; N]) -> [f64; N],
) -> [f64; N] {
    for pair in times.windows(2) {
        let dt = (pair[1] - pair[0]) / 4096.0;
        for index in 0..4096 {
            let time = pair[0] + index as f64 * dt;
            let shift = |change: [f64; N], fraction: f64| {
                std::array::from_fn(|i| state[i] + fraction * dt * change[i])
            };
            let a = rate(time, state);
            let b = rate(time + 0.5 * dt, shift(a, 0.5));
            let c = rate(time + 0.5 * dt, shift(b, 0.5));
            let d = rate(time + dt, shift(c, 1.0));
            for i in 0..N {
                state[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
    }
    state
}

fn bath_reference_program(final_volume: f64) -> [f64; 4] {
    let rows: [[f64; 3]; 4] = [
        [0.0, 300.0, 0.1],
        [0.03125, 800.0, 1.5],
        [0.078125, 250.0, 0.2],
        [0.125, 650.0, 0.8],
    ];
    let logarithmic_rate = (0.01 / final_volume).ln() / 0.125;
    rk4(
        &[0.0, 0.03125, 0.078125, 0.125],
        [400.0, 0.0, 0.0, 0.0],
        |time, state| {
            let upper = rows
                .partition_point(|row| row[0] < time)
                .clamp(1, rows.len() - 1);
            let a = rows[upper - 1];
            let b = rows[upper];
            let fraction = (time - a[0]) / (b[0] - a[0]);
            let bath = a[1] + fraction * (b[1] - a[1]);
            let conductance = a[2] + fraction * (b[2] - a[2]);
            let heat = conductance * (bath - state[0]);
            let work = 0.01 * R * state[0] * logarithmic_rate;
            [
                (work + heat + 0.01 * 12_000.0 * 4.0 * (-4.0 * time).exp()) / (0.03 * R),
                work,
                heat,
                -heat / bath,
            ]
        },
    )
}

#[test]
fn bath_temperature_and_conductance_ramps_refine_compression_expansion_and_rigid_heat() {
    let knots = bath_program();
    let schedule = HeatBathSchedule::new(&knots, || false).unwrap();
    let (bank, mechanism, original) = fixture(-12_000.0, 2000.0);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    for final_volume in [0.005, 0.02, 0.01] {
        let expected = bath_reference_program(final_volume);
        let mut errors = Vec::new();
        for substeps in [8, 16, 32] {
            let mut parcel = original.clone();
            let step = reactor
                .advance_with_heat_bath_schedule(
                    &mut parcel,
                    Volume::new(final_volume),
                    schedule,
                    path(substeps),
                    || false,
                )
                .unwrap();
            errors.push((parcel.state().state().temperature - expected[0]).abs());
            let stored = parcel.mass_kg()
                * (parcel.state().specific_internal_energy_j_kg()
                    - original.state().specific_internal_energy_j_kg());
            close(step.internal_energy_change_j, stored, 5e-6);
            close(
                stored,
                step.boundary_work_on_gas_j + step.heat_input_j + step.energy_residual_j,
                5e-6,
            );
            close(step.energy_residual_j, 0.0, 2e-6);
            assert!(step.entropy_production_j_k >= -1e-7);
            assert_eq!(step.substeps, substeps); // All forcing knots are grid points on these rungs.
            assert_eq!(step.elapsed_time, Time::new(0.125));
            assert_eq!(parcel.volume_m3().to_bits(), final_volume.to_bits());
            assert_eq!(parcel.mass_kg().to_bits(), original.mass_kg().to_bits());
            assert_eq!(parcel.point().axes()["sample_id"], 7.0);
            close(
                parcel
                    .model()
                    .thermodynamics()
                    .mole_composition()
                    .fractions()[1],
                1.0 - (-0.5_f64).exp(),
                5e-7,
            );
            parcel.state().state().validate_coherence().unwrap();
            if substeps == 32 {
                close(parcel.state().state().temperature, expected[0], 0.21);
                close(step.boundary_work_on_gas_j, expected[1], 0.009);
                close(step.heat_input_j, expected[2], 0.046);
                // Entropy must integrate the actual varying bath temperature,
                // not divide the complete heat by its final temperature.
                close(step.thermal_partner_entropy_change_j_k, expected[3], 4.1e-5);
                assert!(
                    (step.thermal_partner_entropy_change_j_k + step.heat_input_j / 650.0).abs()
                        > 1e-3
                );
            }
        }
        assert!(
            errors.windows(2).all(|pair| pair[0] / pair[1] > 3.9),
            "{errors:?}"
        );
    }
}

#[test]
fn heat_only_adaptation_resolves_a_varying_bath_at_exactly_fixed_gas_volume() {
    let knots = bath_program();
    let schedule = HeatBathSchedule::new(&knots, || false).unwrap();
    let expected = bath_reference_program(0.01);
    let (bank, mechanism, original) = fixture(-12_000.0, 2000.0);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut errors = Vec::new();
    for tolerance in [0.01, 0.0001] {
        let mut solve = adaptation(1e6);
        solve.absolute_amount_tolerance_mol = 1e6;
        solve.absolute_heat_tolerance_j = tolerance;
        let mut parcel = original.clone();
        let report = reactor
            .advance_with_heat_bath_schedule_adaptive(
                &mut parcel,
                Volume::new(0.01),
                schedule,
                solve,
                || false,
            )
            .unwrap();
        let physical = report.physical;
        errors.push((parcel.state().state().temperature - expected[0]).abs());
        assert!(report.accepted_intervals > 3 && report.rejected_intervals > 0);
        assert_eq!(
            report.evaluations,
            3 * (report.accepted_intervals + report.rejected_intervals)
        );
        assert_eq!(physical.substeps, 2 * report.accepted_intervals);
        assert!(report.total_chemical_accepted_steps > physical.chemical_accepted_steps);
        assert!(
            report
                .maximum_accepted_error_ratios
                .into_iter()
                .all(|error| error <= 1.0)
        );
        assert_eq!(report.maximum_accepted_error_ratios[4], 0.0);
        assert!(
            report.maximum_accepted_error_ratios[..3]
                .iter()
                .all(|&error| error < 1e-5)
        );
        assert_eq!(parcel.volume_m3().to_bits(), original.volume_m3().to_bits());
        assert_eq!(physical.elapsed_time, schedule.duration());
        assert_eq!(physical.boundary_work_on_gas_j, 0.0);
        close(
            physical.internal_energy_change_j,
            physical.heat_input_j,
            2e-5,
        );
        if tolerance == 0.0001 {
            close(parcel.state().state().temperature, expected[0], 0.007);
            close(physical.heat_input_j, expected[2], 0.0018);
        }
    }
    assert!(errors[0] / errors[1] > 15.0, "{errors:?}");
}

fn latent_curve(bank: &ThermochemicalGasModel, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, temperature_k, liquid_mass_fraction, bulk_density_kg_m3| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: datum + h,
        temperature_k,
        liquid_mass_fraction,
        bulk_density_kg_m3,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        bank.thermodynamics().components()[0].regions()[0]
            .card()
            .content_hash(),
        vec![
            knot(0.0, 200.0, 0.0, 1100.0),
            knot(20_000.0, 400.0, 0.0, 1000.0),
            knot(120_000.0, 400.0, 1.0, 800.0),
            knot(200_000.0, 1200.0, 1.0, 700.0),
        ],
    )
    .unwrap()
}

fn body_reference_program() -> [f64; 3] {
    // Variable gas Cv/R=3+0.004T; the original body remains inside its 400 K
    // latent interval. Chemical formation energy is zero in this fixture.
    rk4(&[0.0, 0.046875, 0.125], [400.0, 0.0, 0.0], |time, state| {
        let conductance = if time <= 0.046875 {
            0.1 + 1.9 * time / 0.046875
        } else {
            2.0 - 1.5 * (time - 0.046875) / (0.125 - 0.046875)
        };
        let work = 0.01 * R * state[0] * 2.0_f64.ln() / 0.125;
        let heat = conductance * (400.0 - state[0]);
        [
            (work + heat) / (0.01 * R * (3.0 + 0.004 * state[0])),
            work,
            heat,
        ]
    })
}

#[test]
fn programmed_body_conductance_preserves_original_latent_inventory_and_variable_gas_cv() {
    let knots = [
        conductance_knot(0.0, 0.1),
        conductance_knot(0.046875, 2.0),
        conductance_knot(0.125, 0.5),
    ];
    let schedule = ConductanceSchedule::new(&knots, || false).unwrap();
    let expected = body_reference_program();
    for datum in [0.0, 1e8] {
        let (bank, mechanism, original) = variable_fixture(0.0);
        let curve = latent_curve(&bank, datum);
        let body = FiniteHeatBody::new(
            &curve,
            Mass::new(0.001),
            Pressure::new(2e5),
            datum + 70_000.0,
            ThermalConductance::new(0.1),
        )
        .unwrap();
        let retuned = body.with_conductance(ThermalConductance::new(0.7)).unwrap();
        assert_eq!(retuned.state(), body.state());
        assert_eq!(retuned.mass(), body.mass());
        assert_eq!(retuned.pressure(), body.pressure());
        assert!(std::ptr::eq(retuned.curve(), body.curve()));
        let mut reactor =
            MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
        let mut errors = Vec::new();
        for tolerance in [0.01, 0.00002] {
            let mut solve = adaptation(1e6);
            solve.absolute_amount_tolerance_mol = 1e6;
            solve.absolute_body_enthalpy_tolerance_j = tolerance;
            let mut parcel = original.clone();
            let report = reactor
                .advance_with_finite_body_schedule_adaptive(
                    &mut parcel,
                    Volume::new(0.005),
                    body,
                    schedule,
                    solve,
                    || false,
                )
                .unwrap();
            let physical = report.physical;
            errors.push((parcel.state().state().temperature - expected[0]).abs());
            assert_eq!(physical.body.state().temperature_k(), 400.0);
            assert_eq!(
                physical.body.conductance().value().to_bits(),
                0.5_f64.to_bits()
            );
            assert_eq!(physical.body.mass(), body.mass());
            assert_eq!(physical.body.pressure(), body.pressure());
            assert!(std::ptr::eq(physical.body.curve(), &curve));
            assert!(physical.body.state().liquid_mass_fraction() > 0.5);
            assert!(report.rejected_intervals > 0);
            assert!(
                report.maximum_accepted_error_ratios[..4]
                    .iter()
                    .all(|&e| e < 1e-5)
            );
            close(
                physical.body_enthalpy_change_j,
                -physical.gas.heat_input_j,
                2e-5,
            );
            close(
                physical.body_internal_energy_change_j + physical.body_boundary_work_j,
                physical.body_enthalpy_change_j,
                2e-5,
            );
            close(
                physical.body_entropy_change_j_k,
                physical.body_enthalpy_change_j / 400.0,
                5e-8,
            );
            close(physical.combined_heat_residual_j, 0.0, 2e-5);
            close(physical.combined_energy_residual_j, 0.0, 2e-5);
            if tolerance == 0.00002 {
                close(parcel.state().state().temperature, expected[0], 0.001);
                close(physical.gas.boundary_work_on_gas_j, expected[1], 0.0001);
                close(physical.gas.heat_input_j, expected[2], 0.0002);
                close(
                    physical.body.state().liquid_mass_fraction(),
                    0.5 - expected[2] / 100.0,
                    2e-6,
                );
            }
        }
        assert!(errors[0] / errors[1] > 20.0, "{errors:?}");
        // Fixed forcing uses the same actual source/body state and terminal G.
        let mut parcel = original.clone();
        let physical = reactor
            .advance_with_finite_body_schedule(
                &mut parcel,
                Volume::new(0.005),
                body,
                schedule,
                path(64),
                || false,
            )
            .unwrap();
        close(parcel.state().state().temperature, expected[0], 0.002);
        assert_eq!(
            physical.body.conductance().value().to_bits(),
            0.5_f64.to_bits()
        );
    }
}

#[test]
fn constant_and_zero_time_programs_preserve_the_original_fixed_and_adaptive_paths() {
    let (bank, mechanism, original) = fixture(-12_000.0, 2000.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 800.0, 0.2);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let knots = [
        bath_knot(0.0, 700.0, 0.2),
        bath_knot(0.037, 700.0, 0.2),
        bath_knot(0.1, 700.0, 0.2),
    ];
    let bath = HeatBathSchedule::new(&knots, || false).unwrap();
    let coefficients = [
        conductance_knot(0.0, 0.2),
        conductance_knot(0.037, 0.2),
        conductance_knot(0.1, 0.2),
    ];
    let program = ConductanceSchedule::new(&coefficients, || false).unwrap();
    let mut a = original.clone();
    let mut b = original.clone();
    let native = reactor
        .advance_with_heat_bath(
            &mut a,
            Volume::new(0.005),
            Time::new(0.1),
            knots[0].bath,
            path(8),
            || false,
        )
        .unwrap();
    let scheduled = reactor
        .advance_with_heat_bath_schedule(&mut b, Volume::new(0.005), bath, path(8), || false)
        .unwrap();
    assert_eq!((a, native), (b, scheduled));
    let mut a = original.clone();
    let mut b = original.clone();
    let native = reactor
        .advance_with_finite_body(
            &mut a,
            Volume::new(0.005),
            Time::new(0.1),
            body,
            path(8),
            || false,
        )
        .unwrap();
    let scheduled = reactor
        .advance_with_finite_body_schedule(
            &mut b,
            Volume::new(0.005),
            body,
            program,
            path(8),
            || false,
        )
        .unwrap();
    assert_eq!((a, native), (b, scheduled));
    let mut a = original.clone();
    let mut b = original.clone();
    let native = reactor
        .advance_with_heat_bath_adaptive(
            &mut a,
            Volume::new(0.01),
            Time::new(0.1),
            knots[0].bath,
            adaptation(0.01),
            || false,
        )
        .unwrap();
    let scheduled = reactor
        .advance_with_heat_bath_schedule_adaptive(
            &mut b,
            Volume::new(0.01),
            bath,
            adaptation(0.01),
            || false,
        )
        .unwrap();
    assert_eq!(scheduled.evaluations, 1);
    assert_eq!((a, native), (b, scheduled));
    let mut a = original.clone();
    let mut b = original.clone();
    let native = reactor
        .advance_with_finite_body_adaptive(
            &mut a,
            Volume::new(0.01),
            Time::new(0.1),
            body,
            adaptation(0.01),
            || false,
        )
        .unwrap();
    let scheduled = reactor
        .advance_with_finite_body_schedule_adaptive(
            &mut b,
            Volume::new(0.01),
            body,
            program,
            adaptation(0.01),
            || false,
        )
        .unwrap();
    assert_eq!((a, native), (b, scheduled));
    let zero = ConductanceSchedule::new(&coefficients[..1], || false).unwrap();
    let mut parcel = original.clone();
    let step = reactor
        .advance_with_finite_body_schedule(
            &mut parcel,
            Volume::new(0.01),
            body,
            zero,
            path(4096),
            || false,
        )
        .unwrap();
    assert_eq!(parcel, original);
    assert_eq!(step.body, body);
    assert_eq!(step.gas.heat_input_j, 0.0);
    assert_eq!(step.gas.boundary_work_on_gas_j, 0.0);
    let coefficients = [conductance_knot(0.0, 0.0), conductance_knot(0.001, -0.0)];
    let program = ConductanceSchedule::new(&coefficients, || false).unwrap();
    let body = body.with_conductance(ThermalConductance::new(0.0)).unwrap();
    let step = reactor
        .advance_with_finite_body_schedule(
            &mut parcel,
            Volume::new(0.01),
            body,
            program,
            path(1),
            || false,
        )
        .unwrap();
    assert_eq!(
        step.body.conductance().value().to_bits(),
        (-0.0f64).to_bits()
    );
}

#[test]
fn forcing_knots_are_exact_global_work_boundaries_including_rigid_strokes() {
    let knots = [
        bath_knot(0.0, 400.0, 0.2),
        bath_knot(0.03125, 800.0, 0.3),
        bath_knot(0.0625, 300.0, 0.4),
        bath_knot(0.09375, 600.0, 0.5),
    ];
    let schedule = HeatBathSchedule::new(&knots, || false).unwrap();
    let (bank, mechanism, original) = fixture(-12_000.0, 2000.0);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut parcel = original.clone();
    let fixed = reactor
        .advance_with_heat_bath_schedule(&mut parcel, Volume::new(0.01), schedule, path(1), || {
            false
        })
        .unwrap();
    assert_eq!(fixed.substeps, 3); // Mandatory knots split the requested one-interval grid.
    let mut loose = adaptation(1e6);
    loose.absolute_amount_tolerance_mol = 1e6;
    loose.minimum_step_s = 0.05; // Exact forcing remainders may be shorter.
    loose.maximum_intervals = 3;
    loose.maximum_evaluations = 9;
    let mut parcel = original.clone();
    let report = reactor
        .advance_with_heat_bath_schedule_adaptive(
            &mut parcel,
            Volume::new(0.01),
            schedule,
            loose,
            || false,
        )
        .unwrap();
    assert_eq!(
        (
            report.accepted_intervals,
            report.rejected_intervals,
            report.evaluations
        ),
        (3, 0, 9)
    );
    assert_eq!(report.physical.substeps, 6);
    assert_eq!(report.minimum_accepted_step_s, 0.03125);
    assert_eq!(report.maximum_accepted_step_s, 0.03125);
    for (intervals, evaluations) in [(2, 9), (3, 8)] {
        let mut solve = loose;
        solve.maximum_intervals = intervals;
        solve.maximum_evaluations = evaluations;
        let mut parcel = original.clone();
        assert!(matches!(
            reactor.advance_with_heat_bath_schedule_adaptive(
                &mut parcel,
                Volume::new(0.01),
                schedule,
                solve,
                || false,
            ),
            Err(AdaptiveError::WorkLimit {
                elapsed_time_s: 0.0625,
                accepted_intervals: 2,
                rejected_intervals: 0,
                evaluations: 6
            })
        ));
        assert_eq!(parcel, original);
    }
    let mut tight = adaptation(1e-12);
    tight.maximum_evaluations = 3;
    let mut parcel = original.clone();
    assert!(matches!(
        reactor.advance_with_heat_bath_schedule_adaptive(
            &mut parcel,
            Volume::new(0.01),
            schedule,
            tight,
            || false,
        ),
        Err(AdaptiveError::WorkLimit {
            elapsed_time_s: 0.0,
            accepted_intervals: 0,
            rejected_intervals: 1,
            evaluations: 3
        })
    ));
    assert_eq!(parcel, original);
    tight.minimum_step_s = tight.initial_step_s;
    tight.maximum_evaluations = 100;
    assert!(matches!(
        reactor.advance_with_heat_bath_schedule_adaptive(
            &mut parcel,
            Volume::new(0.01),
            schedule,
            tight,
            || false,
        ),
        Err(AdaptiveError::Accuracy {
            elapsed_time_s: 0.0,
            step_s: 0.03125,
            ..
        })
    ));
    assert_eq!(parcel, original);
    let crowded = [knots[0], bath_knot(0.03, 600.0, 0.3), knots[3]];
    let crowded = HeatBathSchedule::new(&crowded, || false).unwrap();
    assert!(matches!(
        reactor.advance_with_heat_bath_schedule(
            &mut parcel,
            Volume::new(0.01),
            crowded,
            path(4096),
            || false
        ),
        Err(MaterialCompressionError::Invalid(_))
    ));
    assert_eq!(parcel, original);
}

#[test]
fn invalid_time_and_conductance_programs_refuse_before_publishing_physical_state() {
    assert!(HeatBathSchedule::new(&[], || false).is_err());
    for knots in [
        vec![conductance_knot(0.01, 1.0)],
        vec![conductance_knot(0.0, 1.0), conductance_knot(0.0, 2.0)],
        vec![conductance_knot(0.0, 1.0), conductance_knot(f64::NAN, 2.0)],
        vec![conductance_knot(0.0, -1.0)],
        vec![conductance_knot(0.0, f64::INFINITY)],
    ] {
        assert!(ConductanceSchedule::new(&knots, || false).is_err());
    }
    let maximum: Vec<_> = (0..=MAX_COMPRESSION_SUBSTEPS)
        .map(|i| conductance_knot(i as f64, 1.0))
        .collect();
    assert!(ConductanceSchedule::new(&maximum, || false).is_ok());
    let mut excessive = maximum.clone();
    excessive.push(conductance_knot((MAX_COMPRESSION_SUBSTEPS + 1) as f64, 1.0));
    assert!(ConductanceSchedule::new(&excessive, || false).is_err());
    assert!(matches!(
        ConductanceSchedule::new(&maximum, || true),
        Err(MaterialCompressionError::Cancelled)
    ));
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 800.0, 0.2);
    assert!(
        body.with_conductance(ThermalConductance::new(-1.0))
            .is_err()
    );
    assert!(
        body.with_conductance(ThermalConductance::new(f64::NAN))
            .is_err()
    );
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mismatch = [conductance_knot(0.0, 0.1), conductance_knot(0.1, 1.0)];
    let mismatch = ConductanceSchedule::new(&mismatch, || false).unwrap();
    let mut parcel = original.clone();
    assert!(
        reactor
            .advance_with_finite_body_schedule(
                &mut parcel,
                Volume::new(0.005),
                body,
                mismatch,
                path(8),
                || false
            )
            .is_err()
    );
    assert_eq!(parcel, original);
    let tiny = f64::from_bits(1);
    let knots = [bath_knot(0.0, 400.0, 0.2), bath_knot(tiny, 500.0, 0.3)];
    let schedule = HeatBathSchedule::new(&knots, || false).unwrap();
    let mut solve = adaptation(0.1);
    solve.minimum_step_s = tiny;
    solve.initial_step_s = tiny;
    solve.maximum_step_s = tiny;
    assert!(matches!(
        reactor.advance_with_heat_bath_schedule_adaptive(
            &mut parcel,
            Volume::new(0.01),
            schedule,
            solve,
            || false
        ),
        Err(AdaptiveError::TimeResolution { .. })
    ));
    assert!(matches!(
        reactor.advance_with_heat_bath_schedule(
            &mut parcel,
            Volume::new(0.01),
            schedule,
            path(1),
            || false
        ),
        Err(MaterialCompressionError::Invalid(_))
    ));
    assert_eq!(parcel, original);
}

#[test]
fn late_body_program_cancellation_and_source_failure_roll_back_then_replay() {
    let knots = [
        conductance_knot(0.0, 0.1),
        conductance_knot(0.046875, 2.0),
        conductance_knot(0.125, 0.5),
    ];
    let schedule = ConductanceSchedule::new(&knots, || false).unwrap();
    let (bank, mechanism, original) = variable_fixture(0.0);
    let curve = latent_curve(&bank, 0.0);
    let body = FiniteHeatBody::new(
        &curve,
        Mass::new(0.001),
        Pressure::new(2e5),
        70_000.0,
        ThermalConductance::new(0.1),
    )
    .unwrap();
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut solve = adaptation(1e6);
    solve.absolute_amount_tolerance_mol = 1e6;
    solve.absolute_body_enthalpy_tolerance_j = 0.00002;
    let mut accepted = original.clone();
    let mut polls = 0;
    let expected = reactor
        .advance_with_finite_body_schedule_adaptive(
            &mut accepted,
            Volume::new(0.005),
            body,
            schedule,
            solve,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(expected.rejected_intervals > 0);
    for stop in [polls / 2, polls] {
        let mut parcel = original.clone();
        let mut calls = 0;
        assert!(matches!(
            reactor.advance_with_finite_body_schedule_adaptive(
                &mut parcel,
                Volume::new(0.005),
                body,
                schedule,
                solve,
                || {
                    calls += 1;
                    calls >= stop
                }
            ),
            Err(AdaptiveError::Cancelled)
        ));
        assert_eq!(parcel, original);
        assert_eq!(body.state().specific_enthalpy_j_kg(), 70_000.0);
        assert_eq!(body.conductance().value(), 0.1);
        let retry = reactor
            .advance_with_finite_body_schedule_adaptive(
                &mut parcel,
                Volume::new(0.005),
                body,
                schedule,
                solve,
                || false,
            )
            .unwrap();
        assert_eq!((parcel, retry), (accepted.clone(), expected));
    }
    let knots = [
        conductance_knot(0.0, 0.01),
        conductance_knot(0.025, 0.04),
        conductance_knot(0.1, 0.02),
    ];
    let schedule = ConductanceSchedule::new(&knots, || false).unwrap();
    let prefix_schedule = ConductanceSchedule::new(&knots[..2], || false).unwrap();
    let (bank, mechanism, original) = fixture(0.0, 450.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 300.0, 0.01);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let prefix_volume = 0.01 * 0.5_f64.powf(0.25);
    let mut prefix = original.clone();
    let expected = reactor
        .advance_with_finite_body_schedule(
            &mut prefix,
            Volume::new(prefix_volume),
            body,
            prefix_schedule,
            path(2),
            || false,
        )
        .unwrap();
    assert!(prefix.state().state().temperature > 400.0);
    assert!(expected.body.state().temperature_k() > body.state().temperature_k());
    let mut parcel = original.clone();
    assert!(matches!(
        reactor.advance_with_finite_body_schedule(
            &mut parcel,
            Volume::new(0.005),
            body,
            schedule,
            path(8),
            || false
        ),
        Err(MaterialCompressionError::Kinetics(_) | MaterialCompressionError::Volume(_))
    ));
    assert_eq!(parcel, original);
    let retry = reactor
        .advance_with_finite_body_schedule(
            &mut parcel,
            Volume::new(prefix_volume),
            body,
            prefix_schedule,
            path(2),
            || false,
        )
        .unwrap();
    assert_eq!((parcel, retry), (prefix, expected));
}
