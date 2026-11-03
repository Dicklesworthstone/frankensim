//! Source-pinned film and sphere laws on the original multi-chamber owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    AidingVerticalPlateSettings,
    finite_gas::{
        FiniteGasConvectionSample,
        film::{AidingVerticalPlateGasSurface, FilmGasHeatLaw, FilmGasSurface},
        sphere::ForcedSphereGasSurface,
    },
};
use fs_couple::cavity_phs::piston::GasPiston;
use fs_couple::cavity_phs::piston::thermal::boundary::ports::network::material::state::flow::MaterialGasLawFlowError;
use fs_material::fluid_transport::FluidTransport;
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};

fn gas(fraction: f64, temperature: f64, pressure: f64) -> GasPiston {
    GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(fraction, true),
            source::point(temperature, pressure),
            0.001,
        )
        .unwrap(),
        1.0,
        0.01,
        0.0,
    )
    .unwrap()
}

fn film(model: &ThermochemicalGasModel, law: FilmGasHeatLaw) -> GasHeatLaw<'_> {
    GasHeatLaw::FilmConvection(FilmGasSurface { model, law })
}

fn sphere(model: &ThermochemicalGasModel, area_m2: f64) -> GasHeatLaw<'_> {
    GasHeatLaw::ForcedSphereConvection(ForcedSphereGasSurface {
        model,
        area_m2,
        diameter: Length::new(0.1),
        speed_m_s: 1.0,
    })
}

fn natural(area_m2: f64) -> FilmGasHeatLaw {
    FilmGasHeatLaw::Natural(NaturalGasSurface {
        geometry: fs_convection::NaturalConvectionGeometry::HorizontalCylinder,
        area_m2,
        characteristic_length: Length::new(0.1),
        gravity_m_s2: 9.81,
    })
}

fn aiding(area_m2: f64) -> FilmGasHeatLaw {
    FilmGasHeatLaw::AidingVerticalPlate(AidingVerticalPlateGasSurface {
        area_m2,
        characteristic_length: Length::new(0.1),
        settings: AidingVerticalPlateSettings {
            speed_up_m_s: 0.25,
            gravity_m_s2: 9.81,
            maximum_relative_density_departure: 0.1,
        },
    })
}

fn patch<'a>(id: u64, chamber: u64, cell: CellId, law: GasHeatLaw<'a>) -> ExchangeContact<'a> {
    ExchangeContact {
        contact: Contact { id, chamber, cell },
        law,
    }
}

fn check_source(sample: &GasHeatSample) -> &FiniteGasConvectionSample {
    let (interface, property) = match sample {
        GasHeatSample::FilmConvection(sample) => {
            let property = sample.property_gas.as_ref().unwrap();
            near(
                property.state().temperature,
                0.5 * (sample.interface.body.temperature_k()
                    + sample.interface.gas.state().temperature),
                1e-12,
            );
            assert_eq!(
                sample.interface.evaluation.as_ref().unwrap().transport,
                property.transport_state().unwrap()
            );
            (&sample.interface, property)
        }
        GasHeatSample::ForcedSphereConvection(sample) => {
            let property = sample.surface_gas.as_ref().unwrap();
            assert_eq!(
                property.state().temperature,
                sample.interface.body.temperature_k()
            );
            let transport = sample.interface.gas.transport_state().unwrap();
            let evaluated = sample.interface.evaluation.as_ref().unwrap();
            assert_eq!(evaluated.transport, transport);
            near(
                evaluated.nusselt.groups()["mu_bulk_over_surface"],
                transport.dynamic_viscosity_pa_s()
                    / property.transport_state().unwrap().dynamic_viscosity_pa_s(),
                1e-13,
            );
            (&sample.interface, property)
        }
        _ => panic!("source-pinned patch lost its concrete material evaluation"),
    };
    assert_eq!(property.state().pressure, interface.gas.state().pressure);
    assert_eq!(property.source_packs(), interface.gas.source_packs());
    assert_eq!(
        property.thermodynamics().receipt().mole_fraction_bits(),
        interface
            .gas
            .thermodynamics()
            .receipt()
            .mole_fraction_bits()
    );
    assert_eq!(
        property.components().len(),
        interface.gas.components().len()
    );
    for (query, bulk) in property.components().iter().zip(interface.gas.components()) {
        assert_eq!(query.card_identity(), bulk.card_identity());
        for (query, original) in query.properties().iter().zip(bulk.properties()) {
            assert_eq!(
                query.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
    }
    interface
}

#[test]
fn sourced_network_laws_close_sensible_latent_and_two_chamber_energy_together() {
    fixture_at([260000.0, 180000.0], |assembly, thermal, circuit| {
        scope(|cx, _| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            // The bank's declared fractions deliberately differ from BOTH chambers.
            let bank = source::model_transport(0.25, true);
            let gas = BTreeMap::from([
                (7, gas(1.0, 600.0, 100000.0)),
                (8, gas(0.75, 480.0, 100000.0)),
            ]);
            let saved = gas.clone();
            let chambers = gas
                .iter()
                .map(|(&id, piston)| Chamber {
                    id,
                    piston,
                    input: input(),
                    piston_solve: controls(),
                    exchange: exchange(),
                })
                .collect::<Vec<_>>();
            let patches = [
                patch(1, 7, COLD, sphere(&bank, 0.02)),
                patch(
                    2,
                    8,
                    HOT,
                    film(
                        &bank,
                        FilmGasHeatLaw::Forced(super::super::ForcedGasSurface {
                            area_m2: 0.005,
                            ..surface()
                        }),
                    ),
                ),
                patch(3, 8, HOT, film(&bank, natural(0.005))),
                patch(4, 8, HOT, film(&bank, aiding(0.01))),
            ];
            let controls = NetworkSolve {
                maximum_chambers: 2,
                maximum_contacts: 4,
                maximum_boundary_solves: 2000,
                ..network_controls()
            };
            let next = coupled::step_candidate(
                &prepared,
                cx,
                &thermal,
                DT,
                &circuit,
                &[],
                &chambers,
                &[],
                &patches,
                controls,
                &[0.0; 2],
            )
            .unwrap();
            let mut cell_heat = BTreeMap::<CellId, f64>::new();
            let mut balance = next.thermal().energy_j() - thermal.energy_j();
            for chamber in next.chambers() {
                let mut heat = 0.0;
                for patch in &chamber.surfaces {
                    let sample = check_source(&patch.source.1);
                    assert_eq!(sample.gas, chamber.boundary.piston.midpoint);
                    assert_eq!(
                        sample.body.specific_enthalpy_j_kg(),
                        next.thermal().cells()[&patch.source.0.cell].specific_h_j_kg
                    );
                    near(
                        patch.heat_j,
                        DT * patch.conductance_w_k
                            * (sample.body.temperature_k() - sample.gas.state().temperature),
                        1e-10,
                    );
                    near(
                        patch.conductance_w_k,
                        sample.area_m2
                            * sample
                                .evaluation
                                .as_ref()
                                .unwrap()
                                .coefficient
                                .value
                                .value(),
                        1e-13,
                    );
                    *cell_heat.entry(patch.source.0.cell).or_default() += patch.heat_j;
                    heat += patch.heat_j;
                    if patch.source.0.cell == COLD {
                        assert_eq!(sample.body.temperature_k(), 450.0);
                        assert!(sample.body.liquid_mass_fraction() > 0.5);
                    }
                }
                let before = &gas[&chamber.id];
                let after = &chamber.piston;
                let gas_energy = before.gas().mass_kg()
                    * (after.gas().state().specific_internal_energy_j_kg()
                        - before.gas().state().specific_internal_energy_j_kg());
                let kinetic = 0.5 * before.mass_kg() * after.velocity_m_s().powi(2);
                let outside = 100000.0 * (after.gas().volume_m3() - before.gas().volume_m3());
                near(gas_energy + kinetic + outside, heat, 1e-7);
                balance += gas_energy + kinetic + outside;
            }
            assert!(cell_heat[&HOT] > 0.0 && cell_heat[&COLD] < 0.0);
            for cell in [HOT, COLD] {
                near(
                    MASS * (next.thermal().cells()[&cell].specific_h_j_kg
                        - thermal.cells()[&cell].specific_h_j_kg),
                    -cell_heat[&cell],
                    1e-7,
                );
            }
            near(balance, 0.0, 1e-6);
            assert!(next.step().total_energy_residual_j.abs() <= controls.total_energy_tolerance_j);
            assert_eq!(next.circuit().charges_c(), circuit.charges_c());
            assert_eq!(gas, saved);
            assert_eq!(thermal.time_s(), 0.0);
        })
    });
}

fn transported_patches(bank: &ThermochemicalGasModel) -> [ExchangeContact<'_>; 3] {
    [
        patch(1, 7, HOT, film(bank, natural(0.01))),
        patch(2, 7, HOT, film(bank, aiding(0.01))),
        patch(3, 7, COLD, sphere(bank, 0.02)),
    ]
}

#[test]
fn transported_composition_reaches_sources_and_late_refusals_roll_back_all_owners() {
    fixture_at([120000.0, 50000.0], |assembly, thermal, circuit| {
        scope(|cx, gate| {
            let card = diffusion_card(2e-5);
            let context = QueryPoint::new();
            let pairs = [BinaryDiffusionMaterial {
                species: [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
                material: &card,
                context: &context,
                selection: MaterialPropertySelection::SingleClaimOnly,
            }];
            let faces = [diffusion_face(&pairs)];
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            // B-rich synthetic transport keeps the sphere's actual Pr above 0.71;
            // the bank's unrelated initial fractions must never select that value.
            let gases = BTreeMap::from([
                (7, gas(0.25, 400.0, 100000.0)),
                (8, gas(0.75, 400.0, 100000.0)),
                (9, gas(0.25, 400.0, 60000.0)),
            ]);
            let mut state = MaterialExchangeState::from_parts(thermal, circuit, gases);
            let saved_cells = state.thermal().cells().clone();
            let saved_circuit = state.circuit().clone();
            let saved_gases = state.chambers().clone();
            let initial_thermal_energy = state.thermal().energy_j();
            let drives = transport_drives(&state);
            let controls = NetworkSolve {
                maximum_chambers: 3,
                ..network_controls()
            };
            let bank = source::model_transport(0.6, true);
            let patches = transported_patches(&bank);
            let wrong_bank = source::model_transport(0.6, false);
            let mut unsupported_ratio = patches;
            unsupported_ratio[0].law = sphere(&bank, 0.01); // Gas is colder than HOT.
            for invalid in [transported_patches(&wrong_bank), unsupported_ratio] {
                assert!(matches!(
                    coupled::advance_with_transport(
                        &mut state,
                        &prepared,
                        cx,
                        DT,
                        &[],
                        &drives,
                        &[],
                        &invalid,
                        &[nozzle()],
                        &faces,
                        transport_controls(),
                        controls,
                        &[0.0; 2]
                    ),
                    Err(MaterialGasLawFlowError::Exchange(_))
                ));
                assert_eq!(state.time_s(), 0.0);
                assert_eq!(state.thermal().cells(), &saved_cells);
                assert_eq!(state.circuit(), &saved_circuit);
                assert_eq!(state.chambers(), &saved_gases);
            }
            // Independently stage the original conservative transport, then reuse
            // the original no-flow heat path. Both paths must publish identical owners.
            let parcels = saved_gases
                .iter()
                .map(|(&id, g)| (id, g.gas().clone()))
                .collect();
            let endpoints = drives
                .iter()
                .map(|d| (d.chamber, d.piston_solve.gas))
                .collect();
            let (moved, flow) = shared::step_candidate(
                &parcels,
                &[nozzle()],
                &faces,
                DT,
                &endpoints,
                transport_controls(),
                || false,
            )
            .unwrap()
            .into_parts();
            let staged = MaterialExchangeState::from_parts(
                state.thermal().clone(),
                state.circuit().clone(),
                moved
                    .iter()
                    .map(|(&id, gas)| (id, GasPiston::new(gas.clone(), 1.0, 0.01, 0.0).unwrap()))
                    .collect(),
            );
            let empty_flow = GasNetworkSolve {
                maximum_parcels: 3,
                maximum_links: 0,
                maximum_withdrawn_fraction: 0.05,
                mass_tolerance_kg: 1e-13,
                energy_tolerance_j: 1e-7,
            };
            let (expected, _) = coupled::flow_step_candidate(
                &staged,
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &patches,
                &[],
                empty_flow,
                controls,
                &[0.0; 2],
            )
            .unwrap();
            let report = coupled::advance_with_transport(
                &mut state,
                &prepared,
                cx,
                DT,
                &[],
                &drives,
                &[],
                &patches,
                &[nozzle()],
                &faces,
                transport_controls(),
                controls,
                &[0.0; 2],
            )
            .unwrap();
            assert_eq!(report.flow, flow);
            assert_ne!(moved[&7].model(), saved_gases[&7].gas().model());
            assert_eq!(state.thermal().cells(), expected.thermal().cells());
            assert_eq!(state.circuit(), expected.circuit());
            for (&id, gas) in state.chambers() {
                assert_eq!(gas.gas(), expected.chambers()[&id].gas());
                assert_eq!(gas.velocity_m_s(), expected.chambers()[&id].velocity_m_s());
                assert_eq!(gas.position_m(), expected.chambers()[&id].position_m());
            }
            let mut cell_heat = BTreeMap::<CellId, f64>::new();
            for chamber in &report.exchange.report.source.chambers {
                for patch in &chamber.surfaces {
                    let sample = check_source(&patch.source.1);
                    assert_eq!(sample.gas, chamber.boundary.piston.midpoint);
                    assert_eq!(
                        sample.gas.thermodynamics().receipt().mole_fraction_bits(),
                        moved[&chamber.id]
                            .state()
                            .thermodynamics()
                            .receipt()
                            .mole_fraction_bits()
                    );
                    assert_ne!(
                        sample.gas.thermodynamics().receipt().mole_fraction_bits(),
                        saved_gases[&chamber.id]
                            .gas()
                            .state()
                            .thermodynamics()
                            .receipt()
                            .mole_fraction_bits()
                    );
                    near(
                        patch.heat_j,
                        DT * patch.conductance_w_k
                            * (sample.body.temperature_k() - sample.gas.state().temperature),
                        1e-10,
                    );
                    *cell_heat.entry(patch.source.0.cell).or_default() += patch.heat_j;
                }
            }
            for cell in [HOT, COLD] {
                near(
                    MASS * (state.thermal().cells()[&cell].specific_h_j_kg
                        - saved_cells[&cell].specific_h_j_kg),
                    -cell_heat[&cell],
                    1e-7,
                );
            }
            let mut balance = state.thermal().energy_j() - initial_thermal_energy
                + (state.circuit().charges_c()[0].powi(2) - saved_circuit.charges_c()[0].powi(2))
                    / (2.0 * 0.1);
            for (&id, before) in &saved_gases {
                let after = &state.chambers()[&id];
                balance += after.gas().mass_kg()
                    * after.gas().state().specific_internal_energy_j_kg()
                    - before.gas().mass_kg() * before.gas().state().specific_internal_energy_j_kg()
                    + 0.5 * before.mass_kg() * after.velocity_m_s().powi(2)
                    + before.gas().state().state().pressure
                        * (after.gas().volume_m3() - before.gas().volume_m3());
            }
            near(balance, 0.0, 1e-6);
            assert!(
                report.total_mass_residual_kg.abs()
                    <= transport_controls().exchange.mass_tolerance_kg
            );
            assert!(report.total_energy_residual_j.abs() <= controls.total_energy_tolerance_j);
            let accepted_cells = state.thermal().cells().clone();
            let accepted_gases = state.chambers().clone();
            let accepted_circuit = state.circuit().clone();
            gate.request();
            assert!(matches!(
                coupled::advance_with_transport(
                    &mut state,
                    &prepared,
                    cx,
                    2.0 * DT,
                    &[],
                    &drives,
                    &[],
                    &patches,
                    &[nozzle()],
                    &faces,
                    transport_controls(),
                    controls,
                    &[0.0; 2]
                ),
                Err(MaterialGasLawFlowError::Cancelled)
            ));
            assert_eq!(state.time_s(), DT);
            assert_eq!(state.thermal().cells(), &accepted_cells);
            assert_eq!(state.circuit(), &accepted_circuit);
            assert_eq!(state.chambers(), &accepted_gases);
        })
    });
}
