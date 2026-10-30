//! The full heat assembly shares one native timeline with finite fluid storage.
use super::*;
use fs_couple::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::{
        ContactThermalCoefficients,
        heat::{HeatAssembly, HeatNetworkError, PreparedHeatNetwork, observe_heat_assembly},
    },
};
use fs_thermal::enthalpy::DepositionObserver;
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};

// A complete implicit provider remains valid for implicit drivers, but cannot
// enter an explicit timeline without declaring its physical source-rate bound.
// Neither this provider nor the real coolant assembly needs to be Clone.
struct ImplicitOnly<'a, A: ?Sized>(&'a A);
impl<A: HeatAssembly + ?Sized> HeatAssembly for ImplicitOnly<'_, A> {
    type Sample = A::Sample;
    type Error = A::Error;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.0.prepared_heat_network()
    }
    fn resolve_heat(
        &self,
        cx: &fs_exec::Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.0.resolve_heat(cx, thermal)
    }
}

const EVENT_S: f64 = 0.2;
const END_S: f64 = 0.4;

#[derive(Clone, Debug, Default, PartialEq)]
struct Trace {
    // Start, duration, old wall/fluid, accepted wall/fluid.
    steps: Vec<[f64; 6]>,
    events: Vec<f64>,
}
impl DepositionObserver for Trace {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.events.push(thermal.time_s());
        Ok(())
    }
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        self.steps.push([
            start,
            dt,
            old[&BODY],
            old[&FLUID],
            thermal.temperature_k(BODY)?,
            thermal.temperature_k(FLUID)?,
        ]);
        Ok(())
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-8,
    }
}
fn initial<'a>(
    body: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
    fluid_cp: f64,
) -> DepositionSimulation<'a> {
    let boundary = |source_w| Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    };
    let mut events = [(BODY, 0, CP * 30.0), (FLUID, 1, fluid_cp * 10.0)]
        .map(|(id, material, specific_h_j_kg)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: material,
                body_id: material,
                material_id: MaterialId(material),
                mass_kg: MASS,
                specific_h_j_kg,
                geometry: CellGeometry {
                    start_m: [0.0, material as f64, 0.0],
                    end_m: [0.1, material as f64, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: if id == BODY { 2.0 * SURFACE } else { SURFACE },
                },
                boundary: boundary(0.0),
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: BODY,
            b: FLUID,
            area_m2: SURFACE,
            conductance_w_k: 0.0,
        }),
    });
    // The native scheduled source becomes active after the exact event boundary.
    events.push(TimedEvent {
        time_s: EVENT_S,
        event: Event::Boundary(BODY, boundary(1.0)),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), body.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 4,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, config(0.5), || Ok(()))
    .unwrap()
    .0
}

// Independent flat-plate formula and explicit updates, including combined
// conductance restriction BEFORE clipping to the actual event/end time.
fn coefficient(cp: f64, t: f64) -> f64 {
    let k = 0.5 + 0.002 * (t - 300.0);
    let mu = 0.001 * (1.0 - 0.01 * (t - 300.0));
    0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * cp / k).cbrt() * k / 0.1
}
fn reference(cp: f64, max_dt: f64) -> (Trace, f64) {
    let mut trace = Trace::default();
    let (mut time, mut wall, mut fluid) = (0.0, 310.0, 290.0);
    let mut reservoir_heat = 0.0;
    let wall_capacity = MASS * CP;
    let fluid_capacity = MASS * cp;
    while time < END_S {
        let g = SURFACE * coefficient(cp, 0.5 * (wall + fluid));
        let reservoir_g = SURFACE * coefficient(cp, 0.5 * (wall + 300.0));
        let native = max_dt.min(0.8 / (g / wall_capacity).max(g / fluid_capacity));
        let combined = 1.0 / (1.0 / native + reservoir_g / wall_capacity / 0.8);
        let stop = if time < EVENT_S { EVENT_S } else { END_S };
        let dt = combined.min(stop - time);
        let next_time = if dt == stop - time { stop } else { time + dt };
        let transfer = g * (wall - fluid);
        let reservoir = reservoir_g * (300.0 - wall);
        let authored = if time < EVENT_S { 0.0 } else { 1.0 };
        let next_wall = wall + dt * (reservoir + authored - transfer) / wall_capacity;
        let next_fluid = fluid + dt * transfer / fluid_capacity;
        trace
            .steps
            .push([time, dt, wall, fluid, next_wall, next_fluid]);
        reservoir_heat += dt * reservoir;
        time = next_time;
        wall = next_wall;
        fluid = next_fluid;
        if time == EVENT_S {
            trace.events.push(time);
        }
    }
    (trace, reservoir_heat)
}

#[test]
fn g1_g3_complete_material_heat_and_both_fluid_boundaries_share_native_substeps() {
    let mut outcomes = Vec::new();
    for cp in [4000.0, 8000.0] {
        let (body_card, body_curve) = body_material();
        let (fluid_card, fluid_curve) = fluid_material(1.0, cp);
        let body_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial(&body_curve, &fluid_curve, cp);
        let saved = thermal.cells().clone();
        let saved_links = thermal.links().clone();
        let fluid = liquid_law(&fluid_card);
        let cells = [
            MaterialCell {
                cell: BODY,
                material_id: MaterialId(0),
                curve: body_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &body_card,
                    context: &body_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
            MaterialCell {
                cell: FLUID,
                material_id: MaterialId(1),
                curve: fluid_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &fluid_card,
                    context: &fluid_context,
                    temperature_axis: "T",
                    liquid_fraction_axis: None,
                },
                conductivity: None,
            },
        ];
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let links = [binding()];
        let reservoirs = [SourceConvectiveCell {
            law: SourceCellConvectionLaw::Forced {
                correlation: CorrelationId::FlatPlateLaminarAverage,
                speed_m_s: 1.0,
                settings: ForcedConvectionSettings::default(),
            },
            ..film(BODY)
        }];
        let initial_links = evaluate_source_link_convection(
            &thermal,
            &links,
            1,
            |_, _, t| fluid.at_temperature(t),
            || false,
        )
        .unwrap();
        let gate = fs_exec::CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 179,
                    kernel_id: 29,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat
                .prepare_with_link_conductances(&cx, &thermal, initial_links.conductances_w_k())
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &reservoirs,
                maximum_cells: 1,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let observer = network
                .observe_with_links(
                    &cx,
                    &links,
                    1,
                    |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
                        fluid.at_temperature(t)
                    },
                    Trace::default(),
                    1e-8,
                )
                .unwrap();
            let (next, accepted, audit) = thermal
                .advance_observed(&observer, END_S, config(0.5), || Ok(()))
                .unwrap();
            // The common assembly observer must follow the same independent
            // finite-link + prescribed-reservoir stability restriction and
            // event sequence as the original specialized coolant observer.
            let assembly = network.with_links(
                &links,
                1,
                |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
                    fluid.at_temperature(t)
                },
            );
            let common = observe_heat_assembly(&assembly, &cx, Trace::default(), 1e-8).unwrap();
            let (common_next, common_history, common_audit) = thermal
                .advance_observed(&common, END_S, config(0.5), || Ok(()))
                .unwrap();
            assert_eq!(common_next.cells(), next.cells());
            assert_eq!(common_history.second, accepted.second);
            assert_eq!(common_audit.substeps, audit.substeps);
            near(common_audit.external_heat_j, audit.external_heat_j, 1e-12);
            near(common_history.first.closed_radiation_net_j(), 0.0, 0.0);
            let common_endpoint = common_history.first.sample().unwrap();
            let common_link = &common_endpoint.finite_convection.as_ref().unwrap().samples()[0];
            let specialized_link = &accepted.first.sample().unwrap()
                .finite_convection.as_ref().unwrap().samples()[0];
            near(common_link.wall_temperature_k, specialized_link.wall_temperature_k, 0.0);
            near(common_link.fluid_temperature_k, specialized_link.fluid_temperature_k, 0.0);
            assert_eq!(
                common_link.material.material().identity(),
                specialized_link.material.material().identity(),
            );
            assert!(common.first.sample().is_none() && common.second.steps.is_empty());
            let (common_middle, common_part, first_audit) = thermal
                .advance_observed(&common, EVENT_S, config(0.5), || Ok(()))
                .unwrap();
            let (common_resumed, common_part, last_audit) = common_middle
                .advance_observed(&common_part, END_S, config(0.5), || Ok(()))
                .unwrap();
            assert_eq!(common_resumed.cells(), common_next.cells());
            assert_eq!(common_part.second, common_history.second);
            near(first_audit.external_heat_j + last_audit.external_heat_j,
                common_audit.external_heat_j, 1e-10);

            let undeclared = ImplicitOnly(&assembly);
            // Resolving the implicit physical assembly itself remains valid.
            undeclared.resolve_heat(&cx, &thermal).unwrap();
            let refused = observe_heat_assembly(&undeclared, &cx, Trace::default(), 1e-8).unwrap();
            assert!(matches!(
                thermal.advance_observed(&refused, END_S, config(0.5), || Ok(())),
                Err(ObserverCouplingError::First(CircuitEnthalpyError::Source(
                    ConvectionHeatError::Heat(HeatNetworkError::Invalid(
                        "heat assembly has not declared its explicit source-rate bounds",
                    )),
                )))
            ));
            assert!(refused.first.sample().is_none() && refused.second.steps.is_empty());
            assert_eq!(thermal.cells(), &saved);

            let (expected, reservoir_heat) = reference(cp, 0.5);
            assert_eq!(accepted.second.steps.len(), expected.steps.len());
            for (actual, expected) in accepted.second.steps.iter().zip(&expected.steps) {
                for (&actual, &expected) in actual.iter().zip(expected) {
                    near(actual, expected, 1e-9);
                }
            }
            assert_eq!(accepted.second.events, vec![EVENT_S]);
            assert!(audit.substeps > 1);
            near(
                accepted.first.external_convection_heat_j(),
                reservoir_heat,
                1e-8,
            );
            near(accepted.first.closed_radiation_net_j(), 0.0, 0.0);
            near(
                audit.external_heat_j,
                reservoir_heat + END_S - EVENT_S,
                1e-8,
            );
            near(
                next.energy_j() - thermal.energy_j(),
                audit.external_heat_j,
                1e-8,
            );
            let endpoint = accepted.first.sample().unwrap();
            let applied = accepted.first.applied_sample().unwrap();
            let last = accepted.second.steps.last().unwrap();
            let finite = &endpoint.finite_convection.as_ref().unwrap().samples()[0];
            let starting = &applied.finite_convection.as_ref().unwrap().samples()[0];
            near(starting.wall_temperature_k, last[2], 1e-12);
            near(starting.fluid_temperature_k, last[3], 1e-12);
            near(finite.wall_temperature_k, last[4], 1e-12);
            near(finite.fluid_temperature_k, last[5], 1e-12);
            near(
                finite.property_temperature_k,
                0.5 * (last[4] + last[5]),
                1e-12,
            );
            near(
                endpoint.convection.convection().samples()[0].wall_temperature_k,
                last[4],
                1e-12,
            );
            assert_eq!(
                finite.material.material().card_identity(),
                fluid_card.content_hash()
            );
            for property in finite.material.material().properties() {
                let original = fluid
                    .reference_state()
                    .material()
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    original.answer().receipt.selected
                );
            }
            for sample in &endpoint.heat.material.cells {
                near(
                    sample.phase.temperature_k(),
                    next.temperature_k(sample.cell).unwrap(),
                    1e-12,
                );
            }
            assert_eq!(next.links(), &saved_links);
            assert_eq!(next.cells()[&BODY].boundary.source_w, 1.0);
            assert_eq!(next.cells()[&BODY].boundary.htc_w_m2_k, 0.0);
            assert_eq!(thermal.cells(), &saved);
            assert!(observer.first.sample().is_none() && observer.first.applied_sample().is_none());
            assert!(observer.second.steps.is_empty());

            // Split at the real event boundary: native histories and the
            // generated-coefficient sequence continue identically.
            let (middle, history, _) = thermal
                .advance_observed(&observer, EVENT_S, config(0.5), || Ok(()))
                .unwrap();
            assert!(history.first.sample().is_none());
            assert!(history.first.applied_sample().is_some());
            let first_heat = history.first.external_convection_heat_j();
            let (resumed, history, _) = middle
                .advance_observed(&history, END_S, config(0.5), || Ok(()))
                .unwrap();
            assert_eq!(resumed.cells(), next.cells());
            assert_eq!(history.second, accepted.second);
            near(
                first_heat + history.first.external_convection_heat_j(),
                reservoir_heat,
                1e-8,
            );
            outcomes.push([
                next.temperature_k(BODY).unwrap(),
                next.temperature_k(FLUID).unwrap(),
            ]);

            assert!(matches!(
                thermal.advance_observed_limited(&observer, END_S, config(0.5), 1, || Ok(()),),
                Err(ObserverCouplingError::Thermal(DepositionError::Budget(_)))
            ));
            let duplicated = CoupledObservers::new(observer.first.clone(), observer.first.clone());
            assert!(matches!(
                thermal.advance_observed(&duplicated, END_S, config(0.5), || Ok(()),),
                Err(ObserverCouplingError::ConflictingConductance(LINK))
            ));

            // Electrical storage/work may not vanish through a thermal view.
            let resistors = [CellResistor {
                destination: ResistorCell {
                    branch: 1,
                    cell: BODY,
                },
                selection: MaterialPropertySelection::SingleClaimOnly,
                length_m: 0.1,
                area_m2: 1e-5,
            }];
            let powered_bulk = MaterialNetwork {
                cells: &cells,
                resistors: &resistors,
                links: &[],
            };
            let powered_electrical = ElectricalContactNetwork {
                bulk: &powered_bulk,
                contacts: &[],
                electrical_bases: &[(1, ElectricalPropertyBasis::Resistivity)],
            };
            let powered = MaterialHeatNetwork {
                electrical: &powered_electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let powered = powered
                .prepare_with_link_conductances(&cx, &thermal, initial_links.conductances_w_k())
                .unwrap();
            let rejected = SourceConvectionHeatNetwork {
                heat: &powered,
                cells: &reservoirs,
                maximum_cells: 1,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            assert!(matches!(
                rejected.observe_with_links(
                    &cx,
                    &links,
                    1,
                    |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| fluid
                        .at_temperature(t),
                    (),
                    1e-8,
                ),
                Err(ConvectionHeatError::Heat(HeatNetworkError::Invalid(_)))
            ));

            // A new endpoint's source refusal discards both thermal bodies
            // and the independently composed accepted history.
            let fallible = SourceConvectionHeatNetwork {
                heat: &prepared,
                cells: &reservoirs,
                maximum_cells: 1,
                resolve_fluid: |_: &SourceConvectiveCell,
                                _: &DepositionSimulation<'_>,
                                t: Temperature| {
                    fluid.at_temperature(t).map_err(|_| "reservoir source")
                },
            };
            let refused = fallible
                .observe_with_links(
                    &cx,
                    &links,
                    1,
                    |_: &SourceConvectiveLink, state: &DepositionSimulation<'_>, t: Temperature| {
                        if state.time_s() > 0.0 {
                            return Err("finite source refused");
                        }
                        fluid.at_temperature(t).map_err(|_| "finite source")
                    },
                    Trace::default(),
                    1e-8,
                )
                .unwrap();
            assert!(matches!(
                thermal.advance_observed(&refused, END_S, config(0.5), || Ok(()),),
                Err(ObserverCouplingError::First(
                    ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                        source: "finite source refused",
                        ..
                    })
                ))
            ));
            assert!(refused.first.sample().is_none() && refused.first.applied_sample().is_none());
            assert!(refused.second.steps.is_empty());

            let mut polls = 0;
            thermal
                .advance_observed(&observer, END_S, config(0.5), || {
                    polls += 1;
                    Ok(())
                })
                .unwrap();
            let mut count = 0;
            assert!(
                thermal
                    .advance_observed(&observer, END_S, config(0.5), || {
                        count += 1;
                        if count == polls {
                            Err(DepositionError::Cancelled)
                        } else {
                            Ok(())
                        }
                    })
                    .is_err()
            );
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(thermal.links(), &saved_links);
            assert!(observer.second.steps.is_empty());
        });
    }
    assert!(outcomes[1][0] < outcomes[0][0] && outcomes[1][1] < outcomes[0][1]);
}
