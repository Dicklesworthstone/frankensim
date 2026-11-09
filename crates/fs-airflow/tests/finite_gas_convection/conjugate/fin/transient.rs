//! Original source calorics and films share finite base/fin/fluid inventories.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    ForcedConvectionSettings,
    deposition::{
        fin::{
            FinPlacement, TransientFinBudget, TransientFinCell, TransientFinInput,
            TransientFinTopology,
        },
        source::{
            ConvectionPropertyTemperature, SourceCellConvectionError, SourceCellConvectionLaw,
            fin::{
                PreparedSourceTransientFin, SourceTransientFinError, SourceTransientFinObserver,
                TransientFinSourceControls,
            },
            links::SourceConvectiveLink,
        },
    },
};
use fs_convection::CorrelationId;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergySolve, enthalpy::GasEnthalpyLaw},
};
use fs_qty::Length;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};

const BASE: CellId = CellId(1);
const FLUID: CellId = CellId(2);
const FIN: [CellId; 2] = [CellId(3), CellId(4)];
const SOLID_MATERIAL: MaterialId = MaterialId(10);
const FLUID_MATERIAL: MaterialId = MaterialId(11);

fn boundary() -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: 0.0,
    }
}
fn thermal_controls(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn source_controls() -> TransientFinSourceControls {
    TransientFinSourceControls {
        maximum_links: 5,
        maximum_transverse_biot: 0.1,
        layer: config(),
    }
}
fn topology(initial_k: f64, mass_scale: f64, exposed_tip: bool) -> TransientFinTopology {
    let geometry = StraightFinGeometry::try_new(0.04, 0.00002, 0.042).unwrap();
    let inventory = FIN.map(|id| TransientFinCell {
        id,
        material_id: SOLID_MATERIAL,
        mass_kg: mass_scale * 0.0004,
        specific_h_j_kg: 1000.0 * (initial_k - 200.0),
        boundary: boundary(),
    });
    TransientFinTopology::new(
        TransientFinInput {
            geometry,
            axial_edges_m: &[0.0, 0.02, 0.04],
            inventory: &inventory,
            placement: FinPlacement {
                root_m: [0.0; 3],
                direction: [1.0, 0.0, 0.0],
                body_id: 1,
                segment_id: 1,
            },
            base_cell: BASE,
            fluid_cell: FLUID,
            root_link: 10,
            axial_links: &[11],
            side_links: &[20, 21],
            tip_link: exposed_tip.then_some(22),
            base_half_distance_m: 0.001,
            root_contact_resistance_m2_k_w: 0.002,
        },
        TransientFinBudget {
            maximum_cells: 2,
            maximum_links: 5,
        },
        || Ok(()),
    )
    .unwrap()
}
fn gas_law(model: &ThermochemicalGasModel) -> GasEnthalpyLaw {
    GasEnthalpyLaw::from_mole(
        model,
        model.thermodynamics().declared_composition().clone(),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 440.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 96,
        },
    )
    .unwrap()
}
fn state<'m>(
    topology: &TransientFinTopology,
    solid: &'m EquilibriumEnthalpyPhaseCurve,
    fluid: &'m GasEnthalpyLaw,
) -> DepositionSimulation<'m> {
    state_with_events(topology, solid, fluid, &[])
}
fn state_with_events<'m>(
    topology: &TransientFinTopology,
    solid: &'m EquilibriumEnthalpyPhaseCurve,
    fluid: &'m GasEnthalpyLaw,
    events: &[TimedEvent],
) -> DepositionSimulation<'m> {
    let mut schedule = [
        (BASE, SOLID_MATERIAL, 200000.0),
        (
            FLUID,
            FLUID_MATERIAL,
            fluid.enthalpy_at_temperature(300.0).unwrap(),
        ),
    ]
    .into_iter()
    .map(|(id, material_id, h)| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(DepositedCell {
            id,
            segment_id: id.0,
            body_id: id.0,
            material_id,
            mass_kg: 0.02,
            specific_h_j_kg: h,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [0.001, 0.0, 0.0],
                volume_m3: 0.00002,
                surface_area_m2: 0.02,
            },
            boundary: boundary(),
        }),
    })
    .collect::<Vec<_>>();
    schedule.extend(topology.events_at(0.0, || Ok(())).unwrap());
    schedule.extend_from_slice(events);
    DepositionSimulation::new(
        schedule,
        &[
            (SOLID_MATERIAL, solid as &dyn ThermalLaw),
            (FLUID_MATERIAL, fluid as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 4,
            max_links: 5,
            max_events: 16,
            max_substeps: 20000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_controls(0.01), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn source_fin_rejects_occupied_or_changed_insulated_tip_at_the_final_event() {
    let card = solid(false, 100.0);
    let caloric = caloric(&card);
    let material = layer(&card, &caloric, [250.0, 440.0]);
    let model = source::model_transport(1.0, false);
    let fluid = gas_law(&model);
    let topology = topology(300.0, 1.0, false);
    let extra_tip = Event::Contact(fs_thermal::enthalpy::ThermalLink {
        id: 30,
        a: FIN[1],
        b: BASE,
        area_m2: topology.geometry().cross_section_m2(),
        conductance_w_k: 0.01,
    });
    let occupied = state_with_events(
        &topology,
        &caloric,
        &fluid,
        &[TimedEvent {
            time_s: 0.0,
            event: extra_tip.clone(),
        }],
    );
    assert!(matches!(
        PreparedSourceTransientFin::new(
            &occupied,
            &topology,
            &material,
            &material,
            &films(&topology),
            source_controls(),
            || Ok(()),
        ),
        Err(SourceTransientFinError::Thermal(
            DepositionError::InvalidInput(_)
        ))
    ));
    let bad_boundary = Event::Boundary(
        FIN[1],
        Boundary {
            htc_w_m2_k: 1.0,
            ..boundary()
        },
    );
    for event in [extra_tip, bad_boundary] {
        let initial = state_with_events(
            &topology,
            &caloric,
            &fluid,
            &[TimedEvent {
                time_s: 0.0001,
                event,
            }],
        );
        let saved = initial.cells().clone();
        let prepared = PreparedSourceTransientFin::new(
            &initial,
            &topology,
            &material,
            &material,
            &films(&topology),
            source_controls(),
            || Ok(()),
        )
        .unwrap();
        let observer = SourceTransientFinObserver::new(
            &prepared,
            |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
                fluid.model().evaluate(&source::point(t.value(), 100000.0))
            },
        );
        assert!(matches!(
            initial.advance_observed(&observer, 0.0001, thermal_controls(0.0001), || Ok(())),
            Err(SourceTransientFinError::Thermal(
                DepositionError::InvalidInput(_)
            ))
        ));
        assert_eq!(initial.cells(), &saved);
        assert_eq!(initial.time_s(), 0.0);
        assert!(observer.sample().is_none());
    }
}
fn films(topology: &TransientFinTopology) -> Vec<SourceConvectiveLink> {
    topology
        .side_links()
        .iter()
        .copied()
        .chain(topology.tip_link())
        .map(|link| SourceConvectiveLink {
            link,
            fluid_cell: FLUID,
            characteristic_length: Length::new(0.04),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: SourceCellConvectionLaw::Forced {
                correlation: CorrelationId::FlatPlateLaminarAverage,
                speed_m_s: 0.5,
                settings: ForcedConvectionSettings::default(),
            },
        })
        .collect()
}

#[test]
fn source_fin_stores_startup_heat_on_original_native_inventory_and_retains_endpoint_films() {
    let card = solid(true, 100.0);
    let caloric = caloric(&card);
    let material = layer(&card, &caloric, [250.0, 440.0]);
    let model = source::model_transport(1.0, true);
    let fluid = gas_law(&model);
    let mut first_rises = Vec::new();
    for mass_scale in [1.0, 2.0] {
        let topology = topology(300.0, mass_scale, true);
        let thermal = state(&topology, &caloric, &fluid);
        let saved = thermal.cells().clone();
        let bindings = films(&topology);
        let prepared = PreparedSourceTransientFin::new(
            &thermal,
            &topology,
            &material,
            &material,
            &bindings,
            source_controls(),
            || Ok(()),
        )
        .unwrap();
        let resolver = |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
            fluid.model().evaluate(&source::point(t.value(), 100000.0))
        };
        let initial = prepared.resolve(&thermal, resolver, || Ok(())).unwrap();
        assert!(initial.payload.root_heat_into_fin_w > 0.0);
        near(initial.payload.heat_into_fluid_w, 0.0, 1e-12);
        near(
            initial.payload.net_power_into_fin_w,
            initial.payload.root_heat_into_fin_w,
            1e-12,
        );
        let observer = SourceTransientFinObserver::new(&prepared, resolver);
        let dt = 0.0001;
        let (next, accepted, audit) = thermal
            .advance_observed(&observer, dt, thermal_controls(dt), || Ok(()))
            .unwrap();
        near(next.temperature_k(FLUID).unwrap(), 300.0, 1e-9);
        let stored = FIN
            .iter()
            .map(|id| {
                let a = &thermal.cells()[id];
                let b = &next.cells()[id];
                a.mass_kg * (b.specific_h_j_kg - a.specific_h_j_kg)
            })
            .sum::<f64>();
        near(stored, dt * initial.payload.root_heat_into_fin_w, 1e-10);
        assert!(next.temperature_k(BASE).unwrap() < 400.0);
        assert!(next.temperature_k(FIN[0]).unwrap() > 300.0);
        near(next.energy_j(), thermal.energy_j(), 1e-7);
        near(audit.external_heat_j, 0.0, 1e-12);
        let endpoint = accepted.sample().unwrap();
        assert!(endpoint.heat_into_fluid_w > 0.0);
        for sample in endpoint.sides.samples() {
            near(
                sample.wall_temperature_k,
                next.temperature_k(sample.wall_cell).unwrap(),
                1e-10,
            );
            near(
                sample.property_temperature_k,
                0.5 * (sample.wall_temperature_k + next.temperature_k(FLUID).unwrap()),
                1e-10,
            );
            assert_eq!(
                sample.material.source_packs(),
                fluid
                    .model()
                    .evaluate(&source::point(300.0, 100000.0))
                    .unwrap()
                    .source_packs()
            );
        }
        for sample in &endpoint.conduction {
            assert!(Arc::ptr_eq(&sample.source_lines[0], material.source_line()));
            assert!(Arc::ptr_eq(&sample.source_lines[1], material.source_line()));
        }
        assert_eq!(thermal.cells(), &saved);
        assert!(observer.sample().is_none() && observer.applied_sample().is_none());
        first_rises.push(next.temperature_k(FIN[0]).unwrap() - 300.0);
    }
    near(first_rises[0], 2.0 * first_rises[1], 1e-10);
}

#[test]
fn source_fin_tip_uses_actual_surface_and_axial_half_cell_with_root_contact_once() {
    let card = solid(false, 100.0); // k=60 W/(m K), independently known.
    let caloric = caloric(&card);
    let material = layer(&card, &caloric, [250.0, 440.0]);
    let model = source::model_transport(1.0, false);
    let fluid = gas_law(&model);
    let topology = topology(380.0, 1.0, true);
    let thermal = state(&topology, &caloric, &fluid);
    let prepared = PreparedSourceTransientFin::new(
        &thermal,
        &topology,
        &material,
        &material,
        &films(&topology),
        source_controls(),
        || Ok(()),
    )
    .unwrap();
    let result = prepared
        .resolve(
            &thermal,
            |_, _, t| fluid.model().evaluate(&source::point(t.value(), 100000.0)),
            || Ok(()),
        )
        .unwrap();
    let area = topology.geometry().cross_section_m2();
    let root_g = area / (0.001 / 60.0 + 0.002 + 0.01 / 60.0);
    near(result.payload.root_heat_into_fin_w, root_g * 20.0, 1e-10);
    let tip = &result.payload.tip.as_ref().unwrap().samples()[0];
    let surface = tip.layer.surface_temperature_k;
    let bare_g = tip.layer.film.conductance_w_k;
    let total_g = 1.0 / (0.01 / (60.0 * area) + 1.0 / bare_g);
    assert!(surface > 300.0 && surface < 380.0);
    near(tip.layer.film.wall_temperature_k, surface, 1e-12);
    near(tip.layer.power_to_bath_w, bare_g * (surface - 300.0), 1e-9);
    near(
        tip.layer.solid_power_to_bath_w,
        60.0 * area / 0.01 * (380.0 - surface),
        1e-9,
    );
    near(result.conductances_w_k[&22], total_g, 1e-10);
    near(tip.applied_heat_into_fluid_w, total_g * 80.0, 1e-8);
    near(
        result.payload.heat_into_fluid_w,
        result
            .payload
            .sides
            .samples()
            .iter()
            .map(|s| s.heat_into_fluid_w)
            .sum::<f64>()
            + tip.applied_heat_into_fluid_w,
        1e-10,
    );
    assert!(result.sources_w.is_empty());
    assert_eq!(result.conductances_w_k.len(), 5);
}

#[test]
fn source_fin_refuses_caloric_substitution_biot_missing_films_and_late_source_without_publication()
{
    let card = solid(false, 100.0);
    let caloric = caloric(&card);
    let material = layer(&card, &caloric, [250.0, 440.0]);
    let model = source::model_transport(1.0, false);
    let fluid = gas_law(&model);
    let topology = topology(300.0, 1.0, true);
    let thermal = state(&topology, &caloric, &fluid);
    let saved = thermal.cells().clone();
    let bindings = films(&topology);
    assert!(
        PreparedSourceTransientFin::new(
            &thermal,
            &topology,
            &material,
            &material,
            &bindings[..2],
            source_controls(),
            || Ok(())
        )
        .is_err()
    );
    let other_caloric = caloric.clone();
    let substituted = layer(&card, &other_caloric, [250.0, 440.0]);
    assert!(
        PreparedSourceTransientFin::new(
            &thermal,
            &topology,
            &material,
            &substituted,
            &bindings,
            source_controls(),
            || Ok(())
        )
        .is_err()
    );
    let strict = PreparedSourceTransientFin::new(
        &thermal,
        &topology,
        &material,
        &material,
        &bindings,
        TransientFinSourceControls {
            maximum_transverse_biot: 0.0,
            ..source_controls()
        },
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        strict.resolve(
            &thermal,
            |_, _, t| fluid.model().evaluate(&source::point(t.value(), 100000.0)),
            || Ok(())
        ),
        Err(SourceTransientFinError::TransverseBiot { .. })
    ));
    let prepared = PreparedSourceTransientFin::new(
        &thermal,
        &topology,
        &material,
        &material,
        &bindings,
        source_controls(),
        || Ok(()),
    )
    .unwrap();
    let observer = SourceTransientFinObserver::new(
        &prepared,
        |_: &SourceConvectiveLink, thermal: &DepositionSimulation<'_>, t: Temperature| {
            if thermal.time_s() > 0.0 {
                Err("late original source refusal")
            } else {
                Ok(fluid
                    .model()
                    .evaluate(&source::point(t.value(), 100000.0))
                    .unwrap())
            }
        },
    );
    assert!(matches!(
        thermal.advance_observed(&observer, 0.0001, thermal_controls(0.0001), || Ok(())),
        Err(SourceTransientFinError::Convection(
            SourceCellConvectionError::Source {
                source: "late original source refusal",
                ..
            }
        ))
    ));
    assert_eq!(thermal.cells(), &saved);
    assert!(observer.sample().is_none() && observer.applied_sample().is_none());
    let resolver = |_: &SourceConvectiveLink, _: &DepositionSimulation<'_>, t: Temperature| {
        fluid.model().evaluate(&source::point(t.value(), 100000.0))
    };
    let observer = SourceTransientFinObserver::new(&prepared, resolver);
    let mut polls = 0;
    let expected = thermal
        .advance_observed(&observer, 0.0001, thermal_controls(0.0001), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut at = 0;
    assert!(
        thermal
            .advance_observed(&observer, 0.0001, thermal_controls(0.0001), || {
                at += 1;
                if at == polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(thermal.cells(), &saved);
    let retried = thermal
        .advance_observed(&observer, 0.0001, thermal_controls(0.0001), || Ok(()))
        .unwrap();
    assert_eq!(retried.0.cells(), expected.0.cells());
}
