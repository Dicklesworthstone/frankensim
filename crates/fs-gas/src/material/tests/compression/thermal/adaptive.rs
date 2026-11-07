//! Local error channels and atomic finite-body state in reacting compression.
use super::*;
use crate::material::compression::thermal::adaptive::{
    AdaptiveMaterialThermalCompressionError as Error,
    AdaptiveMaterialThermalCompressionSolve,
};

fn adaptive(tolerance: f64) -> AdaptiveMaterialThermalCompressionSolve {
    AdaptiveMaterialThermalCompressionSolve {
        volume: path(1).volume,
        absolute_temperature_tolerance_k: tolerance,
        absolute_amount_tolerance_mol: 1e-10,
        absolute_work_tolerance_j: 0.25 * tolerance,
        absolute_heat_tolerance_j: 0.25 * tolerance,
        absolute_body_enthalpy_tolerance_j: 0.25 * tolerance,
        relative_tolerance: 0.0,
        initial_step_s: 0.15,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.15,
        maximum_intervals: 4096,
        maximum_evaluations: 12288,
    }
}

fn body_only(tolerance: f64) -> AdaptiveMaterialThermalCompressionSolve {
    let mut solve = adaptive(1e6);
    solve.absolute_amount_tolerance_mol = 1e6;
    solve.absolute_work_tolerance_j = 1e6;
    solve.absolute_heat_tolerance_j = 1e6;
    solve.absolute_body_enthalpy_tolerance_j = tolerance;
    solve.initial_step_s = 0.1;
    solve.maximum_step_s = 0.1;
    solve
}

fn melting_curve(bank: &ThermochemicalGasModel, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    let knot = |h, t, fraction, rho| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: datum + h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: rho,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        bank.thermodynamics().components()[0].regions()[0].card().content_hash(),
        vec![
            knot(0.0, 200.0, 0.0, 1100.0),
            knot(20_000.0, 400.0, 0.0, 1000.0),
            knot(120_000.0, 400.0, 1.0, 800.0),
            knot(200_000.0, 1200.0, 1.0, 700.0),
        ],
    ).unwrap()
}

fn melting_body(curve: &EquilibriumEnthalpyPhaseCurve, datum: f64) -> FiniteHeatBody<'_> {
    FiniteHeatBody::new(
        curve, Mass::new(0.001), Pressure::new(2e5),
        datum + 20_000.0, ThermalConductance::new(1.0),
    ).unwrap()
}

#[test]
fn bath_temperature_work_and_heat_refine_against_the_independent_coupled_solution() {
    for (final_volume, boundary) in [(0.005, 300.0), (0.02, 700.0)] {
        let expected = bath_reference(final_volume, 0.15, boundary, 0.2, 12_000.0);
        let (bank, mechanism, initial) = fixture(-12_000.0, 2000.0);
        let mut reactor =
            MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
        let workspace = reactor.workspace_bytes();
        let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(boundary)).unwrap();
        let mut errors = Vec::new();
        let mut intervals = Vec::new();
        for tolerance in [0.04, 0.0004] {
            let mut parcel = initial.clone();
            let mut controls = adaptive(tolerance);
            if final_volume > 0.01 {
                // In this row only transferred heat drives outer refinement.
                controls.absolute_temperature_tolerance_k = 1e6;
                controls.absolute_amount_tolerance_mol = 1e6;
                controls.absolute_work_tolerance_j = 1e6;
            }
            let report = reactor.advance_with_heat_bath_adaptive(
                &mut parcel, Volume::new(final_volume), Time::new(0.15),
                bath, controls, || false,
            ).unwrap();
            let physical = report.physical;
            errors.push((parcel.state().state().temperature - expected[0]).abs());
            intervals.push(report.accepted_intervals);
            assert!(report.rejected_intervals > 0);
            assert_eq!(physical.substeps, 2 * report.accepted_intervals);
            assert_eq!(report.evaluations, 3 * (report.accepted_intervals + report.rejected_intervals));
            assert!(report.total_chemical_accepted_steps > physical.chemical_accepted_steps);
            assert!(report.total_chemical_rejected_steps >= physical.chemical_rejected_steps);
            assert!(report.maximum_accepted_error_ratios.into_iter().all(|value| value <= 1.0));
            assert_eq!(report.maximum_accepted_error_ratios[4], 0.0);
            if final_volume > 0.01 {
                assert!(report.maximum_accepted_error_ratios[..3].iter().all(|value| *value < 1e-3));
                assert!(report.maximum_accepted_error_ratios[3] > 0.0);
            }
            assert_eq!(physical.elapsed_time, Time::new(0.15));
            assert_eq!(parcel.volume_m3().to_bits(), final_volume.to_bits());
            assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
            let energy = parcel.mass_kg()
                * (parcel.state().specific_internal_energy_j_kg()
                    - initial.state().specific_internal_energy_j_kg());
            close(physical.internal_energy_change_j, energy, 3e-6);
            close(
                energy,
                physical.boundary_work_on_gas_j + physical.heat_input_j + physical.energy_residual_j,
                3e-6,
            );
            close(physical.energy_residual_j, 0.0, 3e-6);
            close(
                physical.thermal_partner_entropy_change_j_k,
                -physical.heat_input_j / boundary, 3e-9,
            );
            assert!(physical.entropy_production_j_k >= -1e-7);
            assert_eq!(physical.heat_input_j.is_sign_positive(), boundary > 400.0);
            parcel.state().state().validate_coherence().unwrap();
            if tolerance == 0.0004 {
                close(parcel.state().state().temperature, expected[0], 0.003);
                close(physical.boundary_work_on_gas_j, expected[1], 0.002);
                close(physical.heat_input_j, expected[2], 0.002);
            }
        }
        assert!(intervals[1] > intervals[0]);
        assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
        assert_eq!(reactor.workspace_bytes(), workspace);
    }
}

#[test]
fn actual_body_enthalpy_alone_resolves_latent_inventory_and_ignores_energy_datums() {
    let expected = bath_reference(0.005, 0.1, 400.0, 1.0, 0.0);
    let expected_fraction = -expected[2] / 100.0; // 0.001 kg * 100 kJ/kg latent interval.
    for datum in [0.0, 1e8] {
        let (bank, mechanism, initial) = fixture(0.0, 2000.0);
        let curve = melting_curve(&bank, datum);
        let original_body = melting_body(&curve, datum);
        let mut reactor =
            MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
        let mut errors = Vec::new();
        let mut intervals = Vec::new();
        for tolerance in [0.02, 0.00002] {
            let mut parcel = initial.clone();
            let report = reactor.advance_with_finite_body_adaptive(
                &mut parcel, Volume::new(0.005), Time::new(0.1),
                original_body, body_only(tolerance), || false,
            ).unwrap();
            let physical = report.physical;
            let fraction = physical.body.state().liquid_mass_fraction();
            errors.push((fraction - expected_fraction).abs());
            intervals.push(report.accepted_intervals);
            // Coarse and fine temperatures are both exactly 400 K while H and
            // phase change. The only active outer comparison is actual body H.
            assert_eq!(physical.body.state().temperature_k(), 400.0);
            assert!(fraction > 0.0 && fraction < 1.0);
            assert!(report.maximum_accepted_error_ratios[..4].iter().all(|value| *value < 1e-3));
            assert!(report.maximum_accepted_error_ratios[4] > 0.0);
            assert!(report.maximum_accepted_error_ratios[4] <= 1.0);
            assert!(report.rejected_intervals > 0);
            assert_eq!(physical.gas.substeps, 2 * report.accepted_intervals);
            assert_eq!(report.evaluations, 3 * (report.accepted_intervals + report.rejected_intervals));
            assert_eq!(physical.body.mass(), original_body.mass());
            assert_eq!(physical.body.pressure(), original_body.pressure());
            assert!(std::ptr::eq(physical.body.curve(), &curve));
            assert!(physical.body_boundary_work_j > 0.0);
            close(physical.body_enthalpy_change_j, -physical.gas.heat_input_j, 3e-6);
            close(
                physical.body_internal_energy_change_j + physical.body_boundary_work_j,
                physical.body_enthalpy_change_j, 3e-6,
            );
            let gas_u = parcel.mass_kg()
                * (parcel.state().specific_internal_energy_j_kg()
                    - initial.state().specific_internal_energy_j_kg());
            close(
                gas_u + physical.body_internal_energy_change_j + physical.body_boundary_work_j,
                physical.gas.boundary_work_on_gas_j, 3e-6,
            );
            close(physical.combined_heat_residual_j, 0.0, 3e-6);
            close(physical.combined_energy_residual_j, 0.0, 3e-6);
            if tolerance == 0.00002 {
                close(fraction, expected_fraction, 3e-6);
                close(parcel.state().state().temperature, expected[0], 0.003);
            }
        }
        assert!(intervals[1] > intervals[0]);
        assert!(errors[1] < 0.1 * errors[0], "{errors:?}");
        assert_eq!(original_body.state().liquid_mass_fraction(), 0.0);
        assert_eq!(original_body.state().specific_enthalpy_j_kg(), datum + 20_000.0);
    }
}

#[test]
fn collision_species_control_is_independent_of_work_and_thermal_allowances() {
    let (bank, original, initial) = fixture(0.0, 2000.0);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("thermal-compression-collision").unwrap(),
        vec![(id("A"), 2)], vec![(id("A"), 1), (id("B"), 1)],
        ArrheniusRateV1::new(
            Temperature::new(400.0), ReactionProgressRateQuantityV1::new(4.0),
            0.0, MolarEnergyQuantityV1::new(0.0),
        ).unwrap(),
    ).unwrap().irreversible();
    let mechanism = IdealGasKineticsV1::new(
        original.thermodynamics().clone(), "thermal-compression-collision",
        1, vec![reaction], 0.0,
    ).unwrap();
    // N_A'=-4 N_A^2/V(t), independent of temperature despite finite cooling.
    // Stop before half conversion, preserving increasing chemical mixing entropy.
    let integral = 0.1 * (100.0 - 200.0) / (0.5_f64).ln();
    let exact_a = (100.0 + 4.0 * integral).recip();
    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(300.0)).unwrap();
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut errors = Vec::new();
    let mut intervals = Vec::new();
    for amount_tolerance in [1e-5, 1e-8] {
        let mut controls = adaptive(1e6);
        controls.absolute_amount_tolerance_mol = amount_tolerance;
        let mut parcel = initial.clone();
        let report = reactor.advance_with_heat_bath_adaptive(
            &mut parcel, Volume::new(0.005), Time::new(0.1), bath, controls, || false,
        ).unwrap();
        let a = 0.01 * parcel.model().thermodynamics().mole_composition().fractions()[0];
        errors.push((a - exact_a).abs());
        intervals.push(report.accepted_intervals);
        assert!(report.rejected_intervals > 0);
        assert!(report.maximum_accepted_error_ratios[1] > 0.0);
        assert!(report.maximum_accepted_error_ratios[1] <= 1.0);
        for channel in [0, 2, 3, 4] {
            assert!(report.maximum_accepted_error_ratios[channel] < 1e-3);
        }
        assert!(report.physical.heat_input_j < 0.0);
        assert_eq!(parcel.mass_kg().to_bits(), initial.mass_kg().to_bits());
        assert_eq!(parcel.volume_m3(), 0.005);
    }
    assert!(intervals[1] > intervals[0]);
    assert!(errors[1] < 0.15 * errors[0], "{errors:?}");
    assert!(errors[1] < 2e-7, "{errors:?}");
}

#[test]
fn adaptive_finite_body_keeps_original_variable_gas_and_body_heat_capacities() {
    for (final_volume, boundary) in [(0.005, 800.0), (0.02, 300.0)] {
        let (bank, mechanism, mut parcel) = variable_fixture(24_000.0);
        let initial = parcel.clone();
        let curve = body_curve(&bank);
        let body = finite_body(&curve, boundary, 2.0);
        let expected = finite_body_reference(final_volume, boundary);
        let mut reactor =
            MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
        let report = reactor.advance_with_finite_body_adaptive(
            &mut parcel, Volume::new(final_volume), Time::new(0.1),
            body, adaptive(0.001), || false,
        ).unwrap();
        let physical = report.physical;
        close(parcel.state().state().temperature, expected[0], 0.02);
        close(physical.body.state().temperature_k(), expected[1], 0.02);
        close(physical.gas.boundary_work_on_gas_j, expected[2], 0.01);
        close(physical.gas.heat_input_j, expected[3], 0.01);
        let gas_u = parcel.mass_kg()
            * (parcel.state().specific_internal_energy_j_kg()
                - initial.state().specific_internal_energy_j_kg());
        close(
            gas_u + physical.body_internal_energy_change_j + physical.body_boundary_work_j,
            physical.gas.boundary_work_on_gas_j, 3e-6,
        );
        close(physical.gas.energy_residual_j, 0.0, 3e-6);
        close(physical.combined_energy_residual_j, 0.0, 3e-6);
        assert!(physical.body_boundary_work_j * physical.body_enthalpy_change_j > 0.0);
        assert!(report.maximum_accepted_error_ratios.into_iter().all(|value| value <= 1.0));
        assert!(std::ptr::eq(physical.body.curve(), &curve));
        assert_eq!(parcel.model().transport_materials(), reactor.bank().transport_materials());
        parcel.state().state().validate_coherence().unwrap();
    }
}

#[test]
fn fixed_volume_no_op_accuracy_and_global_work_refusals_preserve_original_endpoints() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let curve = melting_curve(&bank, 0.0);
    let body = melting_body(&curve, 0.0);
    let warm_curve = body_curve(&bank);
    let warm_body = finite_body(&warm_curve, 800.0, 2.0);
    let bath = HeatBath::new(ThermalConductance::new(0.2), Temperature::new(300.0)).unwrap();
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut native = original.clone();
    let reference = reactor.advance_with_heat_bath(
        &mut native, Volume::new(0.01), Time::new(0.05), bath, path(99), || false,
    ).unwrap();
    let mut parcel = original.clone();
    let report = reactor.advance_with_heat_bath_adaptive(
        &mut parcel, Volume::new(0.01), Time::new(0.05), bath, adaptive(1e-6), || false,
    ).unwrap();
    assert_eq!(parcel, native);
    assert_eq!(report.physical, reference);
    assert_eq!((report.accepted_intervals, report.evaluations), (1, 1));
    assert_eq!(report.maximum_accepted_error_ratios, [0.0; 5]);
    let mut native = original.clone();
    let reference = reactor.advance_with_finite_body(
        &mut native, Volume::new(0.01), Time::new(0.05), warm_body, path(99), || false,
    ).unwrap();
    let mut parcel = original.clone();
    let report = reactor.advance_with_finite_body_adaptive(
        &mut parcel, Volume::new(0.01), Time::new(0.05), warm_body, adaptive(1e-6), || false,
    ).unwrap();
    assert_eq!(parcel, native);
    assert_eq!(report.physical, reference);
    assert_eq!((report.accepted_intervals, report.evaluations), (1, 1));
    let mut parcel = original.clone();
    let report = reactor.advance_with_finite_body_adaptive(
        &mut parcel, Volume::new(0.01), Time::new(0.0), body, adaptive(1e-6), || false,
    ).unwrap();
    assert_eq!(parcel, original);
    assert_eq!(report.physical.body, body);
    assert_eq!((report.accepted_intervals, report.evaluations), (0, 1));

    let mut limited = body_only(1e6);
    limited.initial_step_s = 0.025;
    limited.maximum_step_s = 0.025;
    for (intervals, evaluations) in [(1, 12), (4, 3)] {
        limited.maximum_intervals = intervals;
        limited.maximum_evaluations = evaluations;
        let mut parcel = original.clone();
        assert!(matches!(
            reactor.advance_with_finite_body_adaptive(
                &mut parcel, Volume::new(0.005), Time::new(0.1), body, limited, || false,
            ),
            Err(Error::WorkLimit {
                elapsed_time_s, accepted_intervals: 1, rejected_intervals: 0, evaluations: 3,
            }) if elapsed_time_s == 0.025
        ));
        assert_eq!(parcel, original);
    }
    let mut rejected = body_only(1e-8);
    rejected.maximum_evaluations = 3;
    let mut parcel = original.clone();
    assert!(matches!(
        reactor.advance_with_finite_body_adaptive(
            &mut parcel, Volume::new(0.005), Time::new(0.1), body, rejected, || false,
        ),
        Err(Error::WorkLimit {
            elapsed_time_s, accepted_intervals: 0, rejected_intervals: 1, evaluations: 3,
        }) if elapsed_time_s == 0.0
    ));
    assert_eq!(parcel, original);
    rejected.maximum_evaluations = 12;
    rejected.minimum_step_s = 0.1;
    assert!(matches!(
        reactor.advance_with_finite_body_adaptive(
            &mut parcel, Volume::new(0.005), Time::new(0.1), body, rejected, || false,
        ),
        Err(Error::Accuracy { elapsed_time_s, step_s, error_ratios })
            if elapsed_time_s == 0.0 && step_s == 0.1 && error_ratios[4] > 1.0 && error_ratios[..4].iter().all(|value| *value < 1.0)
    ));
    assert_eq!(parcel, original);

    let mut unresolved = body_only(1.0);
    let tiny = f64::from_bits(1);
    unresolved.minimum_step_s = tiny;
    unresolved.initial_step_s = tiny;
    unresolved.maximum_step_s = tiny;
    assert!(matches!(
        reactor.advance_with_finite_body_adaptive(
            &mut parcel, Volume::new(0.005), Time::new(tiny), body, unresolved, || false,
        ),
        Err(Error::TimeResolution { elapsed_time_s, target_time_s })
            if elapsed_time_s == 0.0 && target_time_s == tiny
    ));
    assert_eq!(parcel, original);
    for field in 0..3 {
        let mut invalid = adaptive(0.01);
        match field {
            0 => invalid.absolute_body_enthalpy_tolerance_j = 0.0,
            1 => invalid.absolute_heat_tolerance_j = f64::NAN,
            _ => invalid.minimum_step_s = 1.0,
        }
        assert!(matches!(
            reactor.advance_with_heat_bath_adaptive(
                &mut parcel, Volume::new(0.005), Time::new(0.1), bath, invalid, || false,
            ),
            Err(Error::Invalid(_))
        ));
        assert_eq!(parcel, original);
    }
}

#[test]
fn rejected_body_history_cancellation_and_late_source_failure_publish_neither_endpoint() {
    let (bank, mechanism, original) = fixture(0.0, 2000.0);
    let curve = melting_curve(&bank, 0.0);
    let body = melting_body(&curve, 0.0);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let solve = body_only(0.00002);
    let mut accepted = original.clone();
    let mut polls = 0;
    let expected = reactor.advance_with_finite_body_adaptive(
        &mut accepted, Volume::new(0.005), Time::new(0.1), body, solve, || {
            polls += 1;
            false
        },
    ).unwrap();
    assert!(expected.rejected_intervals > 0);
    for stop in [1, polls / 2, polls] {
        let mut parcel = original.clone();
        let mut calls = 0;
        assert_eq!(
            reactor.advance_with_finite_body_adaptive(
                &mut parcel, Volume::new(0.005), Time::new(0.1), body, solve, || {
                    calls += 1;
                    calls >= stop
                },
            ).unwrap_err(),
            Error::Cancelled,
        );
        assert_eq!(parcel, original);
        assert_eq!(body.state().liquid_mass_fraction(), 0.0);
        let report = reactor.advance_with_finite_body_adaptive(
            &mut parcel, Volume::new(0.005), Time::new(0.1), body, solve, || false,
        ).unwrap();
        assert_eq!(parcel, accepted);
        assert_eq!(report, expected);
    }

    let (bank, mechanism, original) = fixture(0.0, 450.0);
    let curve = body_curve(&bank);
    let body = finite_body(&curve, 300.0, 0.05);
    let mut reactor =
        MaterialThermalCompressionReactor::new(bank, mechanism, kinetics(), || false).unwrap();
    let mut solve = body_only(1e6);
    solve.initial_step_s = 0.0125;
    solve.maximum_step_s = 0.0125;
    let first_volume = 0.01 * (-(2.0_f64).ln() / 8.0).exp();
    let mut prefix = original.clone();
    reactor.advance_with_finite_body_adaptive(
        &mut prefix, Volume::new(first_volume), Time::new(0.0125), body, solve, || false,
    ).unwrap();
    assert!(prefix.state().state().temperature > 400.0);
    let mut parcel = original.clone();
    assert!(matches!(
        reactor.advance_with_finite_body_adaptive(
            &mut parcel, Volume::new(0.005), Time::new(0.1), body, solve, || false,
        ),
        Err(Error::Physical(MaterialCompressionError::Volume(_)))
            | Err(Error::Physical(MaterialCompressionError::Kinetics(_)))
    ));
    assert_eq!(parcel, original);
    reactor.advance_with_finite_body_adaptive(
        &mut parcel, Volume::new(first_volume), Time::new(0.0125), body, solve, || false,
    ).unwrap();
    assert_eq!(parcel, prefix);
}
