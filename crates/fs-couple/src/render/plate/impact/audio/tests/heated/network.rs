//! Spatial heat transport remains the original owner while its motion is heard.
use super::*;
use crate::render::plate::impact::audio::network::ThermalNetworkSource;
use crate::render::plate::impact::material::MaterialContact;
use crate::render::plate::impact::material::heating::network::ReservoirCell;
use fs_contact::normal_patch::source::SourcePointContact;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    ThermalLink, TimedEvent,
};
use std::collections::BTreeMap;

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 1.0 / f64::from(RATE) / 4.0,
        energy_tolerance_j: 1e-13,
    }
}

fn bindings() -> Vec<ReservoirCell> {
    [("liquid", 1), ("sink", 2)]
        .into_iter()
        .map(|(name, id)| ReservoirCell {
            reservoir: name.into(),
            cell: CellId(id),
        })
        .collect()
}

fn frozen_contact(
    _: &MaterialContact,
    _: &[CaloricReservoir<'_>],
) -> Result<Option<SourcePointContact>, ImpactError> {
    Ok(None)
}

fn heated_pair<'m>(
    source: &ResolvedLiquidState,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
) -> HeatedImpactSystem<'m> {
    let stores = ["liquid", "sink"]
        .into_iter()
        .map(|name| {
            CaloricReservoir::new(
                name.into(),
                curve,
                curve.state_at_specific_enthalpy(0.0).unwrap(),
                0.001,
                1e-6,
                0.0006,
            )
            .unwrap()
        })
        .collect();
    HeatedImpactSystem::new_with_films(
        mechanical(source),
        stores,
        vec![],
        vec![SqueezeFilmHeatRoute {
            film: 0,
            to_a: "liquid".into(),
            to_b: "sink".into(),
            fraction_to_a: 1.0,
        }],
        1e-13,
    )
    .unwrap()
}

fn liquid_pair<'m>(
    card: &'m MaterialCard,
    curve: &'m EquilibriumEnthalpyPhaseCurve,
) -> LiquidFilmImpactSystem<'m> {
    let law = law(card);
    LiquidFilmImpactSystem::new(
        heated_pair(law.reference_state(), curve),
        vec![LiquidFilmBinding {
            film: 0,
            reservoir: "liquid".into(),
            law,
            maximum_relative_viscosity_change: 0.1,
        }],
    )
    .unwrap()
}

fn thermal(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for id in [1, 2] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: 0.001,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64 + 0.01, 0.0, 0.0],
                    volume_m3: 1e-6,
                    surface_area_m2: 0.0006,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 1e-6,
            conductance_w_k: 1.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn rendered_heated_network_matches_direct_steps_and_original_audits() {
    let card = card(0.01);
    let law = law(&card);
    let curve = chart(&card, 0.0002);
    let radiation = artifact();
    let mut direct = heated_pair(law.reference_state(), &curve);
    let mut direct_thermal = thermal(&curve);
    let routes = bindings();
    let gate = CancelGate::new_clock_free();
    let source = ThermalNetworkSource::new(
        heated_pair(law.reference_state(), &curve),
        thermal(&curve),
        (),
        bindings(),
        config(),
        frozen_contact,
    );
    let mut renderer = observer(source, &radiation);
    let mut peak = 0.0_f64;
    for _ in 0..256 {
        let expected = direct
            .step_with_thermal_network(
                &mut direct_thermal,
                &mut (),
                &routes,
                config(),
                &[0.0; 2],
                &gate,
                frozen_contact,
            )
            .unwrap();
        let mut pressure = [0.0];
        renderer.block(&mut pressure).unwrap();
        peak = peak.max(pressure[0].abs());
        let actual = renderer.mechanics();
        assert_eq!(
            actual.source().mechanics().state(),
            direct.mechanics().state()
        );
        assert_eq!(actual.source().reservoirs(), direct.reservoirs());
        assert_eq!(actual.thermal().cells(), direct_thermal.cells());
        assert_eq!(actual.thermal().time_s(), direct_thermal.time_s());
        let frame = renderer.last_mechanical_frame().unwrap();
        assert_eq!(frame.network, expected.network);
        assert_eq!(frame.external_heat_j, expected.external_heat_j);
        assert_eq!(
            frame.impact.thermal.impact.deposited_heat_j,
            expected.impact.thermal.impact.deposited_heat_j
        );
        assert_eq!(
            frame.impact.thermal.impact.mechanical.sample,
            renderer.samples_rendered()
        );
        assert!(frame.network.substeps >= 4);
    }
    assert!(peak > 1e-6);
    assert!(
        renderer
            .mechanics()
            .thermal()
            .temperature_k(CellId(2))
            .unwrap()
            > 300.0
    );
    assert_eq!(
        renderer.mechanics().source().reservoirs()[1].received_heat_j(),
        0.0
    );
}

#[test]
fn spatial_conduction_changes_live_film_viscosity_motion_and_pressure() {
    let card = card(0.001);
    let law = law(&card);
    let curve = chart(&card, 0.0002);
    let radiation = artifact();
    let source = ThermalNetworkSource::new(
        liquid_pair(&card, &curve),
        thermal(&curve),
        (),
        bindings(),
        config(),
        frozen_contact,
    );
    let mut cooled = observer(source, &radiation);
    let mut adiabatic = observer(coupled(&card, &curve), &radiation);
    let mut difference = 0.0_f64;
    for _ in 0..512 {
        let (mut a, mut b) = ([0.0], [0.0]);
        cooled.block(&mut a).unwrap();
        adiabatic.block(&mut b).unwrap();
        difference = difference.max((a[0] - b[0]).abs());
        let retained = cooled.mechanics();
        let live = retained.source();
        let stores = live.system().reservoirs();
        let expected = law
            .at_temperature(Temperature::new(stores[0].state().temperature_k()))
            .unwrap();
        assert_eq!(live.sources()[0], expected);
        assert_eq!(
            live.system()
                .mechanics()
                .squeeze_film(0)
                .unwrap()
                .viscosity_pa_s(),
            Some(expected.dynamic_viscosity_pa_s())
        );
        for (index, id) in [1, 2].into_iter().enumerate() {
            assert_eq!(
                stores[index].state().specific_enthalpy_j_kg(),
                retained.thermal().cells()[&CellId(id)].specific_h_j_kg
            );
        }
        assert!(
            cooled
                .last_mechanical_frame()
                .unwrap()
                .impact
                .thermal
                .impact
                .combined_energy_residual_j
                .abs()
                < 1e-12
        );
    }
    let live = cooled.mechanics().source();
    assert!(live.system().reservoirs()[1].state().temperature_k() > 300.0);
    assert_eq!(live.system().reservoirs()[1].received_heat_j(), 0.0);
    assert!(
        live.sources()[0].dynamic_viscosity_pa_s()
            > adiabatic.mechanics().sources()[0].dynamic_viscosity_pa_s()
    );
    assert!(
        difference > 1e-6,
        "conduction-dependent viscosity must reach the pressure observer"
    );
}

#[derive(Clone, Default)]
struct RefuseSecondSubstep {
    steps: usize,
}

impl DepositionObserver for RefuseSecondSubstep {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.steps += 1;
        if self.steps == 2 {
            Err(DepositionError::MaterialDomain)
        } else {
            Ok(())
        }
    }
}

#[test]
fn late_network_refusal_preserves_all_retained_owners_and_output() {
    let card = card(0.001);
    let curve = chart(&card, 0.0002);
    let radiation = artifact();
    let live = liquid_pair(&card, &curve);
    let heat = thermal(&curve);
    let motion = live.state().to_vec();
    let stores = live.system().reservoirs().to_vec();
    let sources = live.sources().to_vec();
    let viscosity = live
        .system()
        .mechanics()
        .squeeze_film(0)
        .unwrap()
        .viscosity_pa_s();
    let cells = heat.cells().clone();
    let source = ThermalNetworkSource::new(
        live,
        heat,
        RefuseSecondSubstep::default(),
        bindings(),
        config(),
        frozen_contact,
    );
    let mut renderer = observer(source, &radiation);
    let mut output = [42.0];
    let error = renderer.block(&mut output).unwrap_err();
    assert!(error.to_string().contains("MaterialDomain"), "{error}");
    let retained = renderer.mechanics();
    assert_eq!(retained.samples_rendered(), 0);
    assert_eq!(retained.source().state(), motion);
    assert_eq!(retained.source().system().reservoirs(), stores);
    assert_eq!(retained.source().sources(), sources);
    assert_eq!(
        retained
            .source()
            .system()
            .mechanics()
            .squeeze_film(0)
            .unwrap()
            .viscosity_pa_s(),
        viscosity
    );
    assert_eq!(retained.thermal().cells(), &cells);
    assert_eq!(retained.thermal().time_s(), 0.0);
    assert_eq!(retained.observer().steps, 0);
    assert_eq!(renderer.samples_rendered(), 0);
    assert!(renderer.last_mechanical_frame().is_none());
    assert_eq!(output, [42.0]);
    assert!(matches!(
        renderer.block(&mut output),
        Err(RenderError::Poisoned)
    ));
}
