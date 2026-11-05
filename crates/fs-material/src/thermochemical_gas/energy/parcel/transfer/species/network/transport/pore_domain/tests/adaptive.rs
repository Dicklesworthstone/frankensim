//! G1/G3: native dusty-gas breakthrough against independent network eigenmodes.
use super::super::adaptive::{PoreAdvanceError, PoreTimeControls};
use super::*;

fn accuracy(species_kg: f64) -> PoreTimeControls {
    PoreTimeControls {
        gas_species_tolerance_kg: species_kg,
        gas_temperature_tolerance_k: 1e-6,
        relative_tolerance: 0.0,
    }
}

#[test]
fn tighter_accuracy_resolves_three_volume_breakthrough_and_retains_whole_window_transfers() {
    let (domain, throats, card, context) = fixture();
    let pairs = [binding(&card, &context)];
    let faces = throats
        .each_ref()
        .map(|m| face(m, THROAT_AREA, &pairs, &context));
    let duration = 4.0;
    let solve = controls(1.0);
    let rate = equivalent_diffusivity(&throats[0]) * THROAT_AREA
        / (throats[0].thickness_m() * domain.pore_gas(0).unwrap().volume_m3());
    // The equal-mass binary DGM limit has independent eigenvalues 0, k and 3k.
    let (slow, fast) = ((-rate * duration).exp(), (-3.0 * rate * duration).exp());
    let exact = [
        0.4 + 0.3 * slow + 0.1 * fast,
        0.4 - 0.2 * fast,
        0.4 - 0.3 * slow + 0.1 * fast,
    ];
    let loose = domain
        .step_candidate_adaptive_to(
            duration,
            &faces,
            &[energy(); 3],
            solve,
            accuracy(1e-6),
            || false,
        )
        .unwrap();
    let tight = domain
        .step_candidate_adaptive_to(
            duration,
            &faces,
            &[energy(); 3],
            solve,
            accuracy(1e-9),
            || false,
        )
        .unwrap();
    let error = |candidate: &super::super::adaptive::PoreAdvanceCandidate| {
        exact
            .iter()
            .enumerate()
            .map(|(id, expected)| (x(&candidate.domain().parcels()[&(id as u64)]) - expected).abs())
            .fold(0.0_f64, f64::max)
    };
    assert!(error(&tight) < error(&loose) / 8.0);
    assert!(
        error(&tight) > 0.0,
        "local control is not exact time integration"
    );
    assert!(tight.report().rejected_intervals > 0);
    assert!(tight.report().evaluations > loose.report().evaluations);
    for candidate in [&loose, &tight] {
        let report = candidate.report();
        let physical = &report.physical;
        assert_eq!(candidate.domain().time_s(), duration);
        assert_eq!(report.interval_s, [0.0, duration]);
        assert_eq!(physical.duration_s, duration);
        assert_eq!(physical.accepted_substeps, 2 * report.accepted_intervals);
        assert_eq!(physical.attempts, report.evaluations);
        assert!(
            report
                .maximum_accepted_error_ratios
                .iter()
                .all(|&ratio| ratio <= 1.0)
        );
        assert!(
            physical.maximum_species_mass_residual_kg <= solve.transport.exchange.mass_tolerance_kg
        );
        assert!(
            physical.total_mass_residual_kg.abs() <= solve.transport.exchange.mass_tolerance_kg
        );
        assert!(
            physical.total_energy_residual_j.abs() <= solve.transport.exchange.energy_tolerance_j
        );
        assert!(physical.entropy_change_j_k > 0.0);
        assert_eq!(
            physical.neglected_throat_volume_m3,
            loose.report().physical.neglected_throat_volume_m3
        );
        let a = SpeciesId::new("A").unwrap();
        let b = SpeciesId::new("B").unwrap();
        for face in 0..2 {
            let forward = physical
                .faces
                .iter()
                .find(|row| row.face == face && row.species == a)
                .unwrap();
            let reverse = physical
                .faces
                .iter()
                .find(|row| row.face == face && row.species == b)
                .unwrap();
            assert!(forward.mass_left_to_right_kg > 0.0);
            assert_eq!(forward.mass_right_to_left_kg, 0.0);
            assert!(reverse.mass_right_to_left_kg > 0.0);
            assert_eq!(reverse.mass_left_to_right_kg, 0.0);
            let inventory_change = if face == 0 {
                a_mass(domain.left_reservoir()) - a_mass(candidate.domain().left_reservoir())
            } else {
                a_mass(candidate.domain().right_reservoir()) - a_mass(domain.right_reservoir())
            };
            close(forward.mass_left_to_right_kg, inventory_change, 3e-14);
            close(reverse.mass_right_to_left_kg, inventory_change, 3e-14);
        }
    }
    assert_eq!(domain.time_s(), 0.0);
}

#[test]
fn native_refusals_count_as_probes_and_late_work_or_cancellation_preserves_every_pore() {
    let (mut domain, throats, card, context) = fixture();
    let original = domain.clone();
    let pairs = [binding(&card, &context)];
    let faces = throats
        .each_ref()
        .map(|m| face(m, THROAT_AREA, &pairs, &context));
    let mut limited = controls(0.5);
    limited.maximum_attempts = 3;
    let failure =
        domain.advance_adaptive_to(1.0, &faces, &[energy(); 3], limited, accuracy(1e-6), || {
            false
        });
    assert!(matches!(
        failure,
        Err(PoreAdvanceError::WorkLimit {
            time_s: 0.5,
            accepted_substeps: 2,
            evaluations: 3,
        })
    ));
    assert_eq!(domain, original);

    let mut withdrawal = controls(0.25);
    withdrawal.transport.exchange.maximum_species_fraction = 0.002;
    let mut polls = 0;
    let mut accepted = domain.clone();
    let report = accepted
        .advance_adaptive_to(
            0.25,
            &faces,
            &[energy(); 3],
            withdrawal,
            accuracy(1e-6),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(accepted.time_s(), 0.25);
    assert!(report.rejected_intervals > 0);
    assert!(report.evaluations > 3 * report.accepted_intervals);
    assert!(report.physical.maximum_accepted_step_s < withdrawal.maximum_step_s);
    assert_eq!(report.physical.attempts, report.evaluations);

    let mut count = 0;
    let cancelled = domain.advance_adaptive_to(
        0.25,
        &faces,
        &[energy(); 3],
        withdrawal,
        accuracy(1e-6),
        || {
            count += 1;
            count == polls
        },
    );
    assert_eq!(cancelled, Err(PoreAdvanceError::Cancelled));
    assert_eq!(domain, original);
}
