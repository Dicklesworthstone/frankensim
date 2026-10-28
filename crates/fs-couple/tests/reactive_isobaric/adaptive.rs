//! Independent reactor trajectories through the existing source fixtures.
use super::*;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::isobaric::advance::{
    IsobaricAdvanceError, IsobaricAdvanceSolve,
};

fn accuracy() -> IsobaricAdvanceSolve {
    IsobaricAdvanceSolve {
        absolute_amount_mol: 2e-5,
        absolute_temperature_k: 0.15,
        absolute_volume_m3: 4e-7,
        relative_tolerance: 0.0,
        maximum_substeps: 256,
        maximum_endpoint_evaluations: 511,
        maximum_chemical_evaluations: 131072,
    }
}

#[test]
fn adaptive_stiff_and_transient_chains_match_independent_exact_kinetics() {
    for k in [20.0, 10000.0] {
        let channels = channels(&[("ab", [-1, 1, 0], k), ("bc", [0, -1, 1], k)],
            [1; 3], [0.0; 3]);
        let law = implicit(GasReaction::combine(channels.clone()).unwrap());
        let initial = gas(&law, [1.0, 0.0, 0.0], 400.0);
        let result = law.isobaric_adaptive_candidate(&initial, 1.0, 0.0,
            controls(), accuracy(), || false).unwrap();
        let n = initial.mass_kg() / 0.028;
        let alpha = k * R * 400.0 / 100000.0;
        let slow = (-alpha).exp();
        let fast = (-3.0 * alpha).exp();
        let expected = [n * (1.0 / 3.0 + 0.5 * slow + fast / 6.0),
            n * (1.0 - fast) / 3.0, n * (1.0 / 3.0 - 0.5 * slow + fast / 6.0)];
        for (row, expected) in result.report().species.iter().zip(expected) {
            near(row.final_mol, expected, 1.1 * accuracy().absolute_amount_mol);
            assert!(row.final_mol > 0.0, "reachable absent intermediates must be produced");
        }
        near(result.gas().state().state().temperature, 400.0, 1e-7);
        near(result.gas().volume_m3(), initial.volume_m3(), 1e-12);
        assert_eq!(result.gas().mass_kg(), initial.mass_kg());
        assert_eq!(result.gas().state().state().pressure, 100000.0);
        let report = result.report();
        assert!(report.steps.len().is_power_of_two() && report.steps.len() >= 2);
        assert_eq!(report.refused_meshes, 0);
        assert_eq!(report.endpoint_evaluations, 2 * report.steps.len() - 1);
        assert!(report.chemical_evaluations > report.steps.iter().map(|step|
            step.chemistry.boundary.reaction.implicit_mechanism.as_ref().unwrap().evaluations).sum());
        assert!(report.estimated_error_ratio <= 1.0);
        near(report.reaction_extents_mol[0].1, n - report.species[0].final_mol,
            controls().species_tolerance_mol);
        near(report.reaction_extents_mol[1].1, report.species[2].final_mol,
            controls().species_tolerance_mol);
        let reversed = implicit(GasReaction::combine(channels.into_iter().rev().collect()).unwrap());
        assert_eq!(result, reversed.isobaric_adaptive_candidate(&initial, 1.0, 0.0,
            controls(), accuracy(), || false).unwrap());
    }
}

// Independent scalar ODE: n_A=n_A0-xi, n_B=xi, N constant;
// T(t)=[H0+release*xi+power*t]/(N Cp). Activities are x_i at p=p0,
// and the reverse factor is exp[-release/(R T)] from these authored isomers.
fn reference_extent(n: f64, duration: f64, heat: f64, release: f64, k: f64) -> f64 {
    let rate = |time: f64, extent: f64| {
        let temperature = 400.0 + (release * extent + heat * time / duration) / (n * 3.5 * R);
        k * R * temperature / 100000.0
            * (0.75 * n - extent - (-release / (R * temperature)).exp() * extent)
    };
    let dt = duration / 8192.0;
    let mut extent = 0.0;
    for i in 0..8192 {
        let time = i as f64 * dt;
        let a = rate(time, extent);
        let b = rate(time + 0.5 * dt, extent + 0.5 * dt * a);
        let c = rate(time + 0.5 * dt, extent + 0.5 * dt * b);
        let d = rate(time + dt, extent + dt * c);
        extent += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    extent
}

#[test]
fn adaptive_exothermic_heating_and_cooling_close_enthalpy_work_and_heat_entropy() {
    let (release, k, duration) = (10000.0, 100.0, 0.1);
    let law = implicit(channels(&[("ab", [-1, 1, 0], k)], [1; 3],
        [0.0, -release, 0.0]).remove(0));
    let initial = gas(&law, [0.75, 0.0, 0.25], 400.0);
    let n = initial.mass_kg() / 0.028;
    let energy_tolerance = initial.mass_kg() * controls().gas.energy_tolerance_j_kg;
    for heat in [-5.0, 0.0, 5.0] {
        let result = law.isobaric_adaptive_candidate(&initial, duration, heat,
            controls(), accuracy(), || false).unwrap();
        let report = result.report();
        let extent = report.reaction_extents_mol[0].1;
        let reference = reference_extent(n, duration, heat, release, k);
        near(extent, reference, 1.1 * accuracy().absolute_amount_mol);
        let expected_t = 400.0 + (release * extent + heat) / (n * 3.5 * R);
        let reference_t = 400.0 + (release * reference + heat) / (n * 3.5 * R);
        near(result.gas().state().state().temperature, expected_t, 1e-7);
        near(result.gas().state().state().temperature, reference_t,
            1.1 * accuracy().absolute_temperature_k);
        near(result.gas().volume_m3(), n * R * expected_t / 100000.0, 1e-12);
        assert_eq!(result.gas().state().state().pressure, 100000.0);
        assert_eq!(result.gas().point().axes().get("specimen-lot"), Some(&42.0));
        near(report.enthalpy_change_j, heat, energy_tolerance);
        near(report.stored_energy_change_j + report.work_by_gas_j, heat, energy_tolerance);
        near(report.work_by_gas_j, 100000.0 * (result.gas().volume_m3() - initial.volume_m3()), 0.0);
        near(initial.mass_kg() * (result.gas().state().specific_enthalpy_j_kg()
            - initial.state().specific_enthalpy_j_kg()), heat, energy_tolerance);
        assert!(report.maximum_species_residual_mol <= controls().species_tolerance_mol);
        assert!(report.entropy_residual_j_k >= -controls().entropy_tolerance_j_k);
        assert_eq!(report.steps.last().unwrap().end_s, duration);
        // Replaying ONLY accepted original endpoints reproduces the complete
        // candidate. Coarse and refused paths cannot contribute physical heat.
        let mut replay = initial.clone();
        let (mut heat_sum, mut entropy_sum, mut start) = (0.0, 0.0, 0.0);
        for interval in &report.steps {
            assert_eq!(interval.start_s, start);
            let q = interval.chemistry.boundary.boundary.heat_j;
            let step = law.isobaric_implicit_candidate(&replay, interval.end_s - interval.start_s,
                q, controls(), || false).unwrap();
            assert_eq!(step.step(), &interval.chemistry);
            heat_sum += q;
            entropy_sum += q / step.gas().state().state().temperature;
            replay = step.into_parts().0;
            start = interval.end_s;
        }
        assert_eq!(replay, *result.gas());
        near(heat_sum, heat, 1e-14);
        near(report.endpoint_heat_entropy_j_k, entropy_sum, 1e-14);
        if heat != 0.0 {
            assert!((report.endpoint_heat_entropy_j_k
                - heat / result.gas().state().state().temperature).abs() > 1e-5,
                "heat entropy must retain the actual temperature history");
        }
    }
}

#[test]
fn adaptive_mole_changing_reaction_resolves_actual_expanding_volume() {
    let release = 70000.0;
    let law = implicit(channels(&[("dimer", [-2, 1, 0], 10.0)], [1, 2, 1],
        [0.0, -release, 0.0]).remove(0));
    let initial = gas(&law, [0.75, 0.0, 0.25], 400.0);
    // Leave amount and temperature estimates loose, so actual EOS volume
    // controls refinement independently of the other physical observables.
    let accuracy = IsobaricAdvanceSolve { absolute_amount_mol: 1.0,
        absolute_temperature_k: 1e6, absolute_volume_m3: 1e-8, ..accuracy() };
    let result = law.isobaric_adaptive_candidate(&initial, 0.1, 0.0,
        controls(), accuracy, || false).unwrap();
    assert!(result.report().steps.len() > 2);
    let n0 = initial.mass_kg() / 0.028;
    let extent = result.report().reaction_extents_mol[0].1;
    let temperature = (n0 * 3.5 * R * 400.0 + release * extent) / ((n0 - extent) * 3.5 * R);
    near(result.gas().state().state().temperature, temperature, 1e-7);
    near(result.gas().volume_m3(), (n0 - extent) * R * temperature / 100000.0, 1e-12);
    for (row, expected) in result.report().species.iter()
        .zip([0.75 * n0 - 2.0 * extent, extent, 0.25 * n0])
    { near(row.final_mol, expected, controls().species_tolerance_mol); }
    near(result.report().stored_energy_change_j + result.report().work_by_gas_j,
        0.0, initial.mass_kg() * controls().gas.energy_tolerance_j_kg);
}

#[test]
fn refused_attempts_spend_exact_chemical_budget_and_only_fine_history_is_retained() {
    let law = implicit(channels(&[("ab", [-1, 1, 0], 20.0)], [1; 3], [0.0; 3]).remove(0));
    let initial = gas(&law, [1.0, 0.0, 0.0], 400.0);
    let original = initial.clone();
    let single = law.isobaric_implicit_candidate(&initial, 1.0, 0.0, controls(), || false).unwrap();
    let spent = single.step().boundary.reaction.implicit_mechanism.as_ref().unwrap().evaluations;
    let strict = GasReactionSolve { maximum_consumed_fraction: 0.1, ..controls() };
    assert!(matches!(law.isobaric_implicit_candidate(&initial, 1.0, 0.0, strict, || false),
        Err(GasReactionError::ConsumptionLimit { .. })));
    let limited = IsobaricAdvanceSolve { maximum_chemical_evaluations: spent, ..accuracy() };
    assert_eq!(law.isobaric_adaptive_candidate(&initial, 1.0, 0.0, strict, limited, || false),
        Err(IsobaricAdvanceError::ChemicalEvaluationLimit {
            evaluations: spent, maximum: spent, endpoint_evaluations: 1,
        }));
    let loose = IsobaricAdvanceSolve { absolute_amount_mol: 0.1, ..accuracy() };
    let result = law.isobaric_adaptive_candidate(&initial, 1.0, 0.0, strict, loose, || false).unwrap();
    assert!(result.report().refused_meshes > 0);
    assert!(result.report().endpoint_evaluations > result.report().steps.len());
    assert!(result.report().steps.iter().all(|step|
        step.chemistry.boundary.reaction.maximum_consumed_fraction <= 0.1));
    // The unrelated absent species stays exactly absent, not at a floor.
    assert_eq!(result.report().species[2].final_mol, 0.0);
    assert_eq!(result.gas().model().thermodynamics().components().len(), 2);
    assert_eq!(initial, original);
}

#[test]
fn temporal_and_work_limits_and_late_cancellation_never_publish_partial_state() {
    let law = implicit(channels(&[("ab", [-1, 1, 0], 20.0)], [1; 3], [0.0; 3]).remove(0));
    let initial = gas(&law, [1.0, 0.0, 0.0], 400.0);
    let original = initial.clone();
    let unresolved = IsobaricAdvanceSolve { absolute_amount_mol: 1e-20,
        maximum_substeps: 2, ..accuracy() };
    assert!(matches!(law.isobaric_adaptive_candidate(&initial, 1.0, 0.0,
        controls(), unresolved, || false), Err(IsobaricAdvanceError::RefinementLimit {
            substeps: 2, endpoint_evaluations: 3, last_error_ratio: Some(ratio), ..
        }) if ratio > 1.0));
    let endpoint_cap = IsobaricAdvanceSolve { maximum_endpoint_evaluations: 1, ..accuracy() };
    assert!(matches!(law.isobaric_adaptive_candidate(&initial, 1.0, 0.0,
        controls(), endpoint_cap, || false), Err(IsobaricAdvanceError::EndpointEvaluationLimit {
            evaluations: 1, maximum: 1,
        })));
    let mut polls = 0;
    law.isobaric_adaptive_candidate(&initial, 0.1, 1.0, controls(), accuracy(), || {
        polls += 1; false
    }).unwrap();
    let mut replay_polls = 0;
    assert_eq!(law.isobaric_adaptive_candidate(&initial, 0.1, 1.0, controls(), accuracy(), || {
        replay_polls += 1; replay_polls == polls
    }), Err(IsobaricAdvanceError::Cancelled));
    assert!(law.isobaric_adaptive_candidate(&initial, 0.1, 1e9,
        controls(), accuracy(), || false).is_err());
    for invalid in [IsobaricAdvanceSolve { absolute_volume_m3: 0.0, ..accuracy() },
        IsobaricAdvanceSolve { maximum_substeps: 3, ..accuracy() },
        IsobaricAdvanceSolve { maximum_chemical_evaluations: 0, ..accuracy() }]
    { assert!(law.isobaric_adaptive_candidate(&initial, 0.1, 0.0,
        controls(), invalid, || false).is_err()); }
    assert_eq!(initial, original);
}

#[test]
fn equilibrated_reaction_retains_original_heating_cooling_and_common_reference_cancellation() {
    let mut reference = Vec::new();
    for formation in [0.0, 1e12, -1e12] {
        // Identical isomers at equilibrium: every rate/extent stays exactly zero.
        // A common formation reference changes no physical energy difference.
        let law = implicit(channels(&[("ab", [-1, 1, 0], 20.0)],
            [1; 3], [formation; 3]).remove(0));
        let initial = gas(&law, [0.5, 0.5, 0.0], 400.0);
        let saved = initial.clone();
        for (case, heat) in [-0.000123, 0.000123].into_iter().enumerate() {
            let mut direct = initial.clone();
            let original = direct.add_heat_at_constant_pressure(heat, controls().gas, || false).unwrap();
            let endpoint = law.isobaric_implicit_candidate(&initial, 0.1, heat,
                controls(), || false).unwrap();
            assert_eq!(endpoint.gas(), &direct);
            assert_eq!(endpoint.step().enthalpy_change_j, original.enthalpy_change_j);
            assert_eq!(endpoint.step().boundary.reaction.stored_energy_change_j,
                original.stored_energy_change_j);
            assert_eq!(endpoint.step().boundary.reaction.entropy_change_j_k,
                original.entropy_change_j_k);
            assert_eq!(endpoint.step().boundary.reaction.extent_mol, 0.0);

            let adaptive = law.isobaric_adaptive_candidate(&initial, 0.1, heat,
                controls(), accuracy(), || false).unwrap();
            let mut replay = initial.clone();
            for interval in &adaptive.report().steps {
                let q = interval.chemistry.boundary.boundary.heat_j;
                let direct_step = replay.add_heat_at_constant_pressure(q, controls().gas, || false).unwrap();
                assert_eq!(interval.chemistry.enthalpy_change_j, direct_step.enthalpy_change_j);
                assert_eq!(interval.chemistry.boundary.reaction.stored_energy_change_j,
                    direct_step.stored_energy_change_j);
                assert_eq!(interval.chemistry.boundary.reaction.entropy_change_j_k,
                    direct_step.entropy_change_j_k);
                assert_eq!(interval.chemistry.boundary.reaction.extent_mol, 0.0);
            }
            assert_eq!(adaptive.gas(), &replay);
            near(adaptive.report().enthalpy_change_j, heat,
                initial.mass_kg() * controls().gas.energy_tolerance_j_kg);
            let observed = [adaptive.gas().state().state().temperature,
                adaptive.gas().volume_m3(), adaptive.report().enthalpy_change_j,
                adaptive.report().stored_energy_change_j, adaptive.report().entropy_change_j_k];
            if formation == 0.0 { reference.push(observed); }
            else { assert_eq!(observed, reference[case]); }

            // The shared admission also serves fixed-volume chemistry. Its
            // pressure changes, so entropy must retain Cp ln(T1/T0)-R ln(p1/p0).
            let mut fixed = initial.clone();
            let fixed_heat = fixed.add_heat(heat, controls().gas, || false).unwrap();
            let reaction = law.heated_implicit_candidate(&initial, 0.1, heat,
                controls(), || false).unwrap();
            assert_eq!(reaction.gas(), &fixed);
            assert_eq!(reaction.step().reaction.stored_energy_change_j, fixed_heat.stored_energy_change_j);
            let expected_entropy = (initial.mass_kg() / 0.028) * 2.5 * R
                * (fixed.state().state().temperature / 400.0).ln();
            near(reaction.step().reaction.entropy_change_j_k, expected_entropy, 1e-12);
            assert_ne!(fixed.state().state().pressure, initial.state().state().pressure);
            assert_eq!(initial, saved);
        }
    }
}

#[test]
fn frozen_moving_boundary_uses_actual_pressure_entropy_and_refuses_free_compression() {
    use fs_couple::cavity_phs::piston::reactive::boundary::GasReactionBoundary;
    let law = implicit(channels(&[("ab", [-1, 1, 0], 20.0)],
        [1; 3], [0.0; 3]).remove(0));
    let initial = gas(&law, [0.5, 0.5, 0.0], 400.0);
    let expanded = law.boundary_implicit_candidate(&initial, 0.1,
        GasReactionBoundary { final_volume_m3: 2.0 * initial.volume_m3(),
            heat_j: 0.0, work_by_gas_j: 0.0 }, controls(), || false).unwrap();
    near(expanded.gas().state().state().temperature, 400.0, 1e-10);
    near(expanded.gas().state().state().pressure, 50000.0, 1e-8);
    near(expanded.step().reaction.entropy_change_j_k,
        (initial.mass_kg() / 0.028) * R * 2.0_f64.ln(), 1e-12);
    assert!(matches!(law.boundary_implicit_candidate(&initial, 0.1,
        GasReactionBoundary { final_volume_m3: 0.5 * initial.volume_m3(),
            heat_j: 0.0, work_by_gas_j: 0.0 }, controls(), || false),
        Err(GasReactionError::EntropyDecrease { .. })));
}
