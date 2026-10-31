//! Native finite coolant and disjoint reservoirs through nested complete owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::{
    heat::composition::SourceConvectionHeatLayer,
    links::{SourceConvectiveLink, evaluate_source_link_convection},
};
use fs_couple::circuit_enthalpy::contacts::{
    ContactThermalCoefficients,
    heat::{HeatAssembly, HeatNetworkError, PreparedHeatNetwork, observe_heat_assembly},
};
use fs_exec::Cx;
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};

const FLUID: CellId = CellId(72);
const LINK: u64 = 17;
const FLUID_CP: f64 = 4000.0;

// Diagnostic counting changes neither a source result nor physical history.
struct Counted<'a, A: ?Sized> {
    inner: &'a A,
    resolves: Cell<usize>,
    bounds: Cell<usize>,
}
impl<A: HeatAssembly + ?Sized> HeatAssembly for Counted<'_, A> {
    type Sample = A::Sample;
    type Error = A::Error;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.inner.prepared_heat_network()
    }
    fn admit_additional_reservoir_convection(
        &self,
        cx: &Cx<'_>,
        cells: &BTreeSet<CellId>,
    ) -> Result<(), CircuitEnthalpyError<Self::Error>> {
        self.inner.admit_additional_reservoir_convection(cx, cells)
    }
    fn explicit_source_conductance_bounds_w_k(
        &self,
        sample: &Self::Sample,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, CircuitEnthalpyError<Self::Error>> {
        self.bounds.set(self.bounds.get() + 1);
        self.inner
            .explicit_source_conductance_bounds_w_k(sample, checkpoint)
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.resolves.set(self.resolves.get() + 1);
        self.inner.resolve_heat(cx, thermal)
    }
}

// An opaque complete owner may not be assumed free of reservoir convection.
struct Undeclared<'a, A: ?Sized>(&'a A);
impl<A: HeatAssembly + ?Sized> HeatAssembly for Undeclared<'_, A> {
    type Sample = A::Sample;
    type Error = A::Error;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.0.prepared_heat_network()
    }
    fn resolve_heat(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        self.0.resolve_heat(cx, thermal)
    }
}

fn coolant() -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let original = liquid_card(0.0, 0.01, 320.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == SPECIFIC_HEAT_CAPACITY_PROPERTY {
            claim.value = PropertyValue::Scalar {
                value: FLUID_CP,
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &fluid_point(280.0),
            upper: &fluid_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-6,
    )
    .unwrap();
    (card, curve)
}
fn fluid_point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
        .with("p", 100000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap()
}
fn forced() -> SourceCellConvectionLaw {
    SourceCellConvectionLaw::Forced {
        correlation: CorrelationId::FlatPlateLaminarAverage,
        speed_m_s: 1.0,
        settings: ForcedConvectionSettings::default(),
    }
}
fn reservoir(cell: CellId) -> SourceConvectiveCell {
    SourceConvectiveCell {
        law: forced(),
        ..film(cell)
    }
}
fn conductance(film: f64) -> f64 {
    let k = 0.5 + 0.002 * (film - 300.0);
    let mu = 0.001 * (1.0 - 0.01 * (film - 300.0));
    SURFACE * 0.664 * (1000.0 * 0.1 / mu).sqrt() * (mu * FLUID_CP / k).cbrt() * k / 0.1
}
fn initial<'a>(
    body: &'a ResolvedHeatCapacityCurve,
    fluid: &'a ResolvedHeatCapacityCurve,
) -> DepositionSimulation<'a> {
    let mut events = [(BODY, 0, CP * 30.0), (FLUID, 1, FLUID_CP * 10.0)]
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
                    surface_area_m2: 2.0 * SURFACE,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
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
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), body.curve() as &dyn ThermalLaw),
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, step(), || Ok(()))
    .unwrap()
    .0
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 0.005,
        energy_tolerance_j: 1e-8,
    }
}

#[test]
fn nested_complete_sources_share_one_native_timeline_and_retain_both_error_types() {
    let (body_card, body_curve) = body_material();
    let (fluid_card, fluid_curve) = coolant();
    let body_context = body_point(310.0);
    let fluid_context = fluid_point(290.0);
    let thermal = initial(&body_curve, &fluid_curve);
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
    let links = [SourceConvectiveLink {
        link: LINK,
        fluid_cell: FLUID,
        characteristic_length: Length::new(0.1),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: forced(),
    }];
    let initial_links = evaluate_source_link_convection(
        &thermal,
        &links,
        1,
        |_, _, t| fluid.at_temperature(t),
        || false,
    )
    .unwrap();
    let inner_cells = [reservoir(FLUID)];
    let outer_cells = [reservoir(BODY)];
    let both_cells = [outer_cells[0], inner_cells[0]];
    let refuse_inner = Cell::new(false);
    let refuse_outer = Cell::new(false);
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 189,
                kernel_id: 91,
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
            cells: &inner_cells,
            maximum_cells: 1,
            resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t| {
                fluid.at_temperature(t).map_err(|_| "inner reservoir")
            },
        };
        let finite = network.with_links(&links, 1, |_, _, t| {
            if refuse_inner.get() {
                return Err("original finite source refused");
            }
            fluid.at_temperature(t).map_err(|_| "finite source")
        });
        let counted = Counted {
            inner: &finite,
            resolves: Cell::new(0),
            bounds: Cell::new(0),
        };
        let layer = SourceConvectionHeatLayer {
            inner: &counted,
            cells: &outer_cells,
            maximum_cells: 1,
            resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t| {
                if refuse_outer.get() {
                    return Err(7_u8);
                }
                fluid.at_temperature(t).map_err(|_| 8_u8)
            },
        };
        let resolved = layer.resolve_heat(&cx, &thermal).unwrap();
        assert_eq!(counted.resolves.get(), 1);
        near(
            resolved.bulk.conductances_w_k[&LINK],
            conductance(300.0),
            1e-10,
        );
        near(
            resolved.bulk.sources_w[&BODY],
            -10.0 * conductance(305.0),
            1e-10,
        );
        near(
            resolved.bulk.sources_w[&FLUID],
            10.0 * conductance(295.0),
            1e-10,
        );
        let bounds = layer
            .explicit_source_conductance_bounds_w_k(&resolved.bulk.source, &mut || Ok(()))
            .unwrap();
        assert_eq!(counted.resolves.get(), 1); // Bounds reuse the retained sample.
        assert_eq!(counted.bounds.get(), 1);
        near(bounds[&BODY], conductance(305.0), 1e-10);
        near(bounds[&FLUID], conductance(295.0), 1e-10);
        let inner = &resolved.bulk.source.heat;
        assert!(inner.finite_convection.is_some());
        assert_eq!(inner.convection.materials()[0].cell, FLUID);
        assert_eq!(resolved.bulk.source.convection.materials()[0].cell, BODY);
        assert!(std::ptr::eq(layer.prepared_heat_network(), &prepared));

        // Match the original all-in-one physics path on the exact native
        // explicit timeline, then close its independent global heat account.
        let complete = SourceConvectionHeatNetwork {
            heat: &prepared,
            cells: &both_cells,
            maximum_cells: 2,
            resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t| {
                fluid.at_temperature(t).map_err(|_| "flat reservoir")
            },
        };
        let complete = complete.with_links(&links, 1, |_, _, t| {
            fluid.at_temperature(t).map_err(|_| "flat finite")
        });
        let observed = observe_heat_assembly(&layer, &cx, (), 1e-8).unwrap();
        let flat = observe_heat_assembly(&complete, &cx, (), 1e-8).unwrap();
        let (next, accepted, audit) = thermal
            .advance_observed(&observed, DT, step(), || Ok(()))
            .unwrap();
        let (expected, _, expected_audit) = thermal
            .advance_observed(&flat, DT, step(), || Ok(()))
            .unwrap();
        assert_eq!(next.cells(), expected.cells());
        assert_eq!(audit.substeps, expected_audit.substeps);
        near(audit.external_heat_j, expected_audit.external_heat_j, 1e-12);
        near(
            next.energy_j() - thermal.energy_j(),
            audit.external_heat_j,
            1e-8,
        );
        let endpoint = accepted.first.sample().unwrap();
        let finite_sample = &endpoint.heat.finite_convection.as_ref().unwrap().samples()[0];
        near(
            finite_sample.wall_temperature_k,
            next.temperature_k(BODY).unwrap(),
            1e-12,
        );
        near(
            finite_sample.fluid_temperature_k,
            next.temperature_k(FLUID).unwrap(),
            1e-12,
        );
        assert_eq!(
            finite_sample.material.material().card_identity(),
            fluid_card.content_hash()
        );
        for sample in endpoint
            .convection
            .materials()
            .iter()
            .chain(endpoint.heat.convection.materials())
        {
            assert_eq!(
                sample.material.material().card_identity(),
                fluid_card.content_hash()
            );
            for property in sample.material.material().properties() {
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
        }
        let (continued, resumed, later) = next
            .advance_observed(&accepted, 2.0 * DT, step(), || Ok(()))
            .unwrap();
        near(
            continued.energy_j() - thermal.energy_j(),
            audit.external_heat_j + later.external_heat_j,
            1e-8,
        );
        assert!(
            resumed
                .first
                .sample()
                .unwrap()
                .heat
                .finite_convection
                .is_some()
        );
        assert_eq!(continued.links(), &saved_links);

        // Zero flux still owns an exposed boundary; the outer layer cannot
        // reinterpret it as unclaimed area. Refusal precedes inner evaluation.
        let equilibrium_cells = [SourceConvectiveCell {
            reference_temperature: Temperature::new(310.0),
            ..outer_cells[0]
        }];
        let equilibrium = SourceConvectionHeatLayer {
            inner: &counted,
            cells: &equilibrium_cells,
            maximum_cells: 1,
            resolve_fluid: &layer.resolve_fluid,
        };
        let eq = equilibrium.resolve_heat(&cx, &thermal).unwrap();
        near(
            eq.bulk.source.convection.convection().samples()[0].heat_into_cell_w,
            0.0,
            0.0,
        );
        let duplicate = SourceConvectionHeatLayer {
            inner: &equilibrium,
            cells: &outer_cells,
            maximum_cells: 1,
            resolve_fluid: &layer.resolve_fluid,
        };
        let calls = counted.resolves.get();
        assert!(duplicate.resolve_heat(&cx, &thermal).is_err());
        assert_eq!(counted.resolves.get(), calls);
        let unknown = Undeclared(&counted);
        let undeclared = SourceConvectionHeatLayer {
            inner: &unknown,
            cells: &outer_cells,
            maximum_cells: 1,
            resolve_fluid: &layer.resolve_fluid,
        };
        assert!(undeclared.resolve_heat(&cx, &thermal).is_err());
        assert_eq!(counted.resolves.get(), calls);

        refuse_inner.set(true);
        assert!(matches!(
            layer.resolve_heat(&cx, &thermal),
            Err(CircuitEnthalpyError::Source(ConvectionHeatError::Heat(
                ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                    source: "original finite source refused",
                    ..
                })
            )))
        ));
        refuse_inner.set(false);
        refuse_outer.set(true);
        assert!(matches!(
            layer.resolve_heat(&cx, &thermal),
            Err(CircuitEnthalpyError::Source(
                ConvectionHeatError::Convection(SourceCellConvectionError::Source {
                    source: 7,
                    ..
                })
            ))
        ));
        assert!(
            thermal
                .advance_observed(&observed, DT, step(), || Ok(()))
                .is_err()
        );
        assert!(observed.first.sample().is_none());
        assert_eq!(thermal.cells(), &saved);
        assert_eq!(thermal.links(), &saved_links);
        refuse_outer.set(false);
        gate.request();
        assert!(matches!(
            layer.resolve_heat(&cx, &thermal),
            Err(CircuitEnthalpyError::Cancelled)
        ));
    });
}
