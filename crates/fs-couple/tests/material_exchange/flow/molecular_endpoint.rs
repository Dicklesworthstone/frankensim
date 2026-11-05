//! Actual molecular endpoints must retain their original sources after feedback.
use super::*;
use fs_gas::multicomponent::dissipation::DissipationError;

fn limited_molecular(
    original: &Molecular,
    maximum_temperature_k: f64,
    maximum_pressure_pa: f64,
) -> ParcelMolecularTransport<ConstantTransport> {
    let law = ConstantTransport::new(
        original.face().gas().source().clone(),
        "bounded endpoint diffusion",
        1,
        TransportDomain {
            minimum_temperature_k: 200.0,
            maximum_temperature_k,
            minimum_pressure_pa: 1000.0,
            maximum_pressure_pa,
        },
        original.face().law().law().coefficients(),
        vec![0.0, DIFFUSIVITY, DIFFUSIVITY, 0.0],
    )
    .unwrap();
    ParcelMolecularTransport::new(
        PlanarTransport::new(original.face().gas().clone(), law, || false).unwrap(),
    )
}

#[test]
fn zero_flow_faces_refresh_cross_sources_at_the_complete_feedback_endpoint() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let state = MaterialExchangeState::from_parts(
                thermal,
                circuit,
                molecular_gases(&bank, [(0.5, 400.0, 1e5); 2], true),
            );
            let mut drives = drives(state.chambers());
            drives[1].heat_j = 0.02;
            let (next, report) = state
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives,
                    &[],
                    &[],
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    controls(),
                    &[0.0; 2],
                )
                .unwrap();
            // The declared connection had no initial transport. Qualification
            // must still see the subsequent prescribed heat and pressure rise.
            assert_eq!(report.flow.links[0].mass_left_to_right_kg, vec![0.0; 2]);
            assert_eq!(report.flow.links[0].heat_correction_left_to_right_j, 0.0);
            let endpoint = next.chambers()[&20].gas().state().state();
            assert!(endpoint.temperature > 400.01);
            assert!(endpoint.pressure > 1e5);
            let resolved = molecular.face().law().last_resolved();
            close(resolved.temperature().value(), endpoint.temperature, 2e-8);
            close(resolved.pressure().value(), endpoint.pressure, 1e-7);
            close(
                resolved.heats()[0].value(),
                20000.0 + 100.0 * endpoint.temperature,
                2e-6,
            );
            assert_eq!(resolved.mole_fractions(), &[0.5, 0.5]);
            assert_eq!(state.time_s(), 0.0);
        });
    });
}

#[test]
fn final_temperature_or_pressure_outside_the_molecular_source_rolls_back_every_owner() {
    for (maximum_temperature_k, maximum_pressure_pa) in [(400.001, 1e7), (1800.0, 100000.1)] {
        fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
            with_cx(|cx, _| {
                let (bank, original) = molecular_fixture(MASS, 1.0, 0.0);
                let mut molecular =
                    limited_molecular(&original, maximum_temperature_k, maximum_pressure_pa);
                let prepared = assembly.prepare(cx, &thermal).unwrap();
                let mut state = MaterialExchangeState::from_parts(
                    thermal,
                    circuit,
                    molecular_gases(&bank, [(0.5, 400.0, 1e5); 2], true),
                );
                let mut drives = drives(state.chambers());
                drives[1].heat_j = 0.01;
                let parcels = state
                    .chambers()
                    .iter()
                    .map(|(&id, gas)| (id, gas.gas().clone()))
                    .collect();
                let endpoints = drives
                    .iter()
                    .map(|drive| (drive.chamber, drive.piston_solve.gas))
                    .collect();
                let native = molecular
                    .step_candidate(
                        &parcels,
                        &[molecular_link()],
                        GasTransportPaths::default(),
                        DT,
                        &endpoints,
                        molecular_controls(),
                        || false,
                    )
                    .unwrap();
                assert_eq!(native.step().links[0].mass_left_to_right_kg, vec![0.0; 2]);
                // Both constituent steps are otherwise admitted: only the
                // molecular source's final validity is narrower than the gas EOS.
                let (feedback, _) = state
                    .step_candidate(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &[],
                        controls(),
                        &[0.0; 2],
                    )
                    .unwrap();
                let endpoint = feedback.chambers()[&20].gas().state().state();
                assert!(
                    endpoint.temperature > maximum_temperature_k
                        || endpoint.pressure > maximum_pressure_pa
                );
                let old_gases = state.chambers().clone();
                let old_cells = state.thermal().cells().clone();
                let old_circuit = state.circuit().clone();
                assert!(matches!(
                    state.advance_with_molecular_transport(
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &[],
                        &mut molecular,
                        &[molecular_link()],
                        GasTransportPaths::default(),
                        molecular_controls(),
                        controls(),
                        &[0.0; 2],
                    ),
                    Err(MaterialFlowError::Molecular(ParcelMolecularError::Face(
                        DissipationError::Invalid(_)
                    )))
                ));
                assert_eq!(state.chambers(), &old_gases);
                assert_eq!(state.thermal().cells(), &old_cells);
                assert_eq!(state.circuit(), &old_circuit);
                assert_eq!(state.time_s(), 0.0);
                for (&id, gas) in state.chambers() {
                    assert_eq!(
                        gas.last_window_transfers(),
                        old_gases[&id].last_window_transfers()
                    );
                    assert_eq!(
                        gas.last_window_duration_s(),
                        old_gases[&id].last_window_duration_s()
                    );
                }
            });
        });
    }
}

#[test]
fn molecular_qualification_excludes_unconnected_species_banks() {
    fixture([0.2, 0.2], 150.0, |assembly, thermal, circuit| {
        with_cx(|cx, _| {
            let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let mut gases = molecular_gases(&bank, [(0.5, 400.0, 1e5); 2], true);
            let unrelated = source_gas("C", MASS, 1e5).with_fixed_support().unwrap();
            let unrelated_mass = unrelated.gas().mass_kg();
            gases.insert(30, unrelated);
            let state = MaterialExchangeState::from_parts(thermal, circuit, gases);
            let mut settings = controls();
            settings.maximum_chambers = 3;
            let (next, _) = state
                .molecular_transport_step_candidate(
                    &prepared,
                    cx,
                    DT,
                    &[],
                    &drives(state.chambers()),
                    &[],
                    &[],
                    &mut molecular,
                    &[molecular_link()],
                    GasTransportPaths::default(),
                    molecular_controls(),
                    settings,
                    &[0.0; 2],
                )
                .unwrap();
            let gas = next.chambers()[&30].gas();
            assert_eq!(
                gas.model().thermodynamics().components()[0].species(),
                &species_id("C")
            );
            close(gas.mass_kg(), unrelated_mass, 1e-16);
            assert_eq!(next.time_s(), DT);
        });
    });
}

#[test]
fn native_qualification_preserves_the_original_species_identity_and_cancellation() {
    let (bank, mut molecular) = molecular_fixture(MASS, 1.0, 0.0);
    let gases = molecular_gases(&bank, [(0.5, 400.0, 1e5); 2], true);
    molecular
        .qualify_parcel(gases[&10].gas(), molecular_controls().face, || false)
        .unwrap();
    let unsupported = source_gas("C", MASS, 1e5);
    assert!(matches!(
        molecular.qualify_parcel(unsupported.gas(), molecular_controls().face, || false),
        Err(ParcelMolecularError::SourceMismatch(species)) if species == species_id("C")
    ));
    assert!(matches!(
        molecular.qualify_parcel(gases[&10].gas(), molecular_controls().face, || true),
        Err(ParcelMolecularError::Cancelled)
    ));
}
