//! G1/G3: sourced radiation constrains the original finite-cell timeline.
use super::*;
use std::collections::BTreeMap;

const MAX_DT_S: f64 = 0.2;
const TARGET_S: f64 = 0.4;
const SIGMA: f64 = 5.670_374_419e-8;
const IDS: [CellId; 2] = [CellId(30), CellId(10)];

#[derive(Clone, Default)]
struct Trace {
    steps: Vec<(f64, f64, [f64; 2])>,
}
impl DepositionObserver for Trace {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        check()?;
        let temperatures = [
            thermal.temperature_k(IDS[0])?,
            thermal.temperature_k(IDS[1])?,
        ];
        // This fixture's declared step ceiling approaches equilibrium without
        // crossing it. This is not a general temporal-accuracy claim.
        assert!(temperatures[0] <= old[&IDS[0]] + 1e-10);
        assert!(temperatures[1] >= old[&IDS[1]] - 1e-10);
        assert!(temperatures[0] >= 350.0 - 1e-10);
        assert!(temperatures[1] <= 350.0 + 1e-10);
        self.steps.push((start_s, duration_s, temperatures));
        Ok(())
    }
}
fn with_radiating_cells<R>(
    epsilon: f64,
    max_substeps: usize,
    run: impl FnOnce(&Cx<'_>, &MaterialHeatNetwork<'_>, DepositionSimulation<'_>) -> R,
) -> R {
    let cards = [card(epsilon, epsilon), card(epsilon, epsilon)];
    let curves = [curve(&cards[0], 100.0), curve(&cards[1], 100.0)];
    let context = context();
    let cells = [0, 1].map(|i| MaterialCell {
        cell: IDS[i],
        material_id: MaterialId(i as u64),
        curve: &curves[i],
        point: binding(&cards[i], &context),
        conductivity: None,
    });
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
    let surface_cells = [0, 1].map(|i| RadiatingCell {
        surface: i,
        cell: IDS[i],
        material_id: MaterialId(i as u64),
        curve: &curves[i],
        material: cells[i].point,
        selection: MaterialPropertySelection::SingleClaimOnly,
    });
    let enclosure = enclosure([
        emissivity(&cards[0], &curves[0], -20.0),
        emissivity(&cards[1], &curves[1], -80.0),
    ]);
    let radiation = MaterialRadiation {
        enclosure: &enclosure,
        cells: &surface_cells,
        power_tolerance_w: 1e-8,
        radiosity_tolerance_w_m2: 1e-8,
    };
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: Some(&radiation),
    };
    let events = [0, 1]
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: [-20.0, -80.0][i],
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 2.0,
                },
                boundary: Boundary {
                    ambient_k: 350.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 0,
            max_events: 2,
            max_substeps,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: MAX_DT_S,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 73,
                kernel_id: 53,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &assembly, thermal)
    })
}
// Independent equal-area plate formula and scalar explicit-Euler recurrence.
// Each original solid cell has m*Cp=1 J/K. No production radiation/timestep
// helper supplies this reference.
fn plate_reference(epsilon: f64) -> Vec<(f64, f64, [f64; 2])> {
    let mut temperatures = [380.0_f64, 320.0_f64];
    let mut time = 0.0;
    let mut steps = Vec::new();
    while time < TARGET_S {
        let maximum = temperatures[0].max(temperatures[1]);
        let bound = 4.0 * SIGMA * epsilon * maximum.powi(3);
        let proposed = MAX_DT_S / (1.0 + MAX_DT_S * bound / 0.8);
        let remaining = TARGET_S - time;
        let dt = proposed.min(remaining);
        let [hot, cold] = temperatures;
        let outward =
            SIGMA / (2.0 / epsilon - 1.0) * (hot - cold) * (hot + cold) * (hot * hot + cold * cold);
        temperatures = [hot - dt * outward, cold + dt * outward];
        steps.push((time, dt, temperatures));
        time = if dt == remaining { TARGET_S } else { time + dt };
        assert!(steps.len() < 100);
    }
    steps
}
#[test]
fn material_emissivity_sets_native_substeps_and_preserves_closed_finite_heat() {
    let run = |epsilon| {
        with_radiating_cells(epsilon, 100, |cx, assembly, thermal| {
            let selected = assembly.prepare(cx, &thermal).unwrap();
            let initial = selected.resolve(cx, &thermal).unwrap();
            let radiation = initial.bulk.source.radiation.as_ref().unwrap();
            let bounds = radiation.conductance_bounds_w_k(|| Ok(())).unwrap();
            for id in IDS {
                near(bounds[&id], 4.0 * SIGMA * epsilon * 380.0_f64.powi(3));
            }
            let original = thermal.cells().clone();
            let observer = selected.observe(cx, Trace::default(), 1e-8).unwrap();
            let (next, observed, audit) = thermal
                .advance_observed(
                    &observer,
                    TARGET_S,
                    StepConfig {
                        max_dt_s: MAX_DT_S,
                        energy_tolerance_j: 1e-9,
                    },
                    || Ok(()),
                )
                .unwrap();
            let reference = plate_reference(epsilon);
            assert_eq!(audit.substeps, reference.len());
            assert_eq!(observed.second.steps.len(), reference.len());
            for (actual, expected) in observed.second.steps.iter().zip(&reference) {
                near(actual.0, expected.0);
                near(actual.1, expected.1);
                for i in 0..2 {
                    near(actual.2[i], expected.2[i]);
                }
            }
            near(next.energy_j(), thermal.energy_j());
            near(audit.stored_energy_change_j, 0.0);
            near(
                audit.external_heat_j - observed.first.closed_radiation_net_j(),
                0.0,
            );
            let sample = observed.first.sample().unwrap().radiation.as_ref().unwrap();
            for surface in sample.surfaces() {
                near(
                    surface.phase.temperature_k(),
                    next.temperature_k(surface.cell).unwrap(),
                );
                let original_surface = &radiation.surfaces()[surface.surface];
                assert_eq!(
                    surface.material.card_identity(),
                    original_surface.material.card_identity()
                );
                for (actual, original) in surface
                    .material
                    .properties()
                    .iter()
                    .zip(original_surface.material.properties())
                {
                    assert_eq!(
                        actual.answer().receipt.selected,
                        original.answer().receipt.selected
                    );
                }
            }
            assert_eq!(next.time_s(), TARGET_S);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(thermal.cells(), &original);
            assert!(observer.first.sample().is_none());
            assert!(observer.second.steps.is_empty());
            (audit.substeps, next.temperature_k(IDS[0]).unwrap())
        })
    };
    let reflective = run(0.2);
    let emissive = run(0.9);
    assert_eq!((reflective.0, emissive.0), (4, 7));
    assert!(emissive.1 < reflective.1 - 10.0);
}
#[test]
fn radiation_substep_budget_refuses_without_publishing_heat_or_observations() {
    with_radiating_cells(0.9, 1, |cx, assembly, thermal| {
        let selected = assembly.prepare(cx, &thermal).unwrap();
        let observer = selected.observe(cx, Trace::default(), 1e-8).unwrap();
        let original = thermal.cells().clone();
        let failed = thermal.advance_observed(
            &observer,
            TARGET_S,
            StepConfig {
                max_dt_s: MAX_DT_S,
                energy_tolerance_j: 1e-9,
            },
            || Ok(()),
        );
        assert!(matches!(
            failed,
            Err(fs_thermal::observers::ObserverCouplingError::Thermal(
                DepositionError::Budget("thermal substeps")
            ))
        ));
        assert_eq!(thermal.cells(), &original);
        assert_eq!(thermal.time_s(), 0.0);
        assert!(observer.first.sample().is_none());
        assert_eq!(observer.first.closed_radiation_net_j(), 0.0);
        assert!(observer.second.steps.is_empty());
    });
}
