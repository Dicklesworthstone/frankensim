//! Source-bound caloric preparation through the ordinary instrument frontend.
use super::*;

#[test]
fn native_thermal_cells_keep_caloric_sources_and_exact_initial_temperature() {
    scoped(|cx, gate| {
        for fibers in [false, true] {
            let profile = Profile::parse(profile_text(fibers).as_bytes()).unwrap();
            let ready = ready(cx, gate, &profile, fibers, 20.0, 7, 48_000, None);
            let prepared = profile.prepare(&ready, gate).unwrap();
            let network = prepared.network(gate).unwrap();
            for (cell, source) in prepared.profile.cells.iter().zip(&prepared.sources) {
                let initial = cell.temperatures[1];
                let state = source.curve.state_at_specific_enthalpy(0.0).unwrap();
                assert_eq!(state.temperature_k().to_bits(), initial.to_bits());
                assert_eq!(network.temperature_k(cell.id).unwrap().to_bits(), initial.to_bits());
                assert_eq!(state.material_card_identity(), source.card.content_hash());
                let lower = source.source_curves[0].receipts();
                let upper = source.source_curves[1].receipts();
                let retained = source.curve.source_receipts();
                assert_eq!(retained.len(), lower.len() + upper.len() - 2);
                assert_eq!(&retained[..lower.len()], lower);
                assert_eq!(&retained[lower.len()..], &upper[2..]);
                for knot in source.curve.knots() {
                    assert!((knot.specific_enthalpy_j_kg - 20.0 * (knot.temperature_k - initial)).abs() < 1e-11);
                }
                let mass = network.cells()[&cell.id].mass_kg;
                let heated = source.curve.advance_isobaric_heat(state, mass, 100_000.0,
                    mass * 20.0 * 2.0).unwrap();
                assert!((heated.state.temperature_k() - initial - 2.0).abs() < 1e-11);
                // The fixed-density instrument model still has no volume work.
                assert_eq!(heated.boundary_work_j, 0.0);
            }
        }
    });
}

#[test]
fn retaining_caloric_sources_does_not_change_fixed_density_pressure_pcm() {
    scoped(|cx, gate| {
        let profile = Profile::parse(profile_text(false).as_bytes()).unwrap();
        let mut original = ready(cx, gate, &profile, false, 20.0, 7, 48_000, None);
        let mut sourced = ready(cx, gate, &profile, false, 20.0, 7, 48_000, None);
        let mut legacy = profile.prepare(&original, gate).unwrap();
        let prepared = profile.prepare(&sourced, gate).unwrap();
        // Reconstruct the previous anonymous-chart path ONLY in this regression.
        // No anonymous fallback is added to the production source admission.
        for source in &mut legacy.sources {
            source.curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
                source.card.content_hash(), SolidLiquidPhase::Solid,
                source.curve.knots().to_vec(),
            ).unwrap();
            assert!(source.curve.source_receipts().is_empty());
        }
        let (old_wav, old_summary) = legacy.write(&mut original,
            Cursor::new(Vec::new()), None, gate).unwrap();
        let (new_wav, new_summary) = prepared.write(&mut sourced,
            Cursor::new(Vec::new()), None, gate).unwrap();
        assert_eq!(old_summary.samples, new_summary.samples);
        assert_eq!(old_wav.into_inner(), new_wav.into_inner());
    });
}

#[path = "gas_convection_tests.rs"]
mod gas_convection;
