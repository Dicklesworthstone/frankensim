//! Nonlinear original NASA calorics with reciprocal finite-pipe conduction.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::conduction::GraphConductiveLink;
use fs_conduction::transient::deposition::{DepositionError, backward_euler::BackwardEulerStep};

fn linked(law: &GasEnthalpyLaw, temperatures: [f64; 6]) -> DepositionSimulation<'_> {
    let initial = thermal(law, temperatures, None);
    let mut schedule: Vec<_> = initial
        .cells()
        .values()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    schedule.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 30,
            a: CellId(2),
            b: CellId(3),
            area_m2: AREA,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 16,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn stagnant_gas_conduction_keeps_original_nasa_storage_pressure_and_composition() {
    context(|cx| {
        let law = gas_law();
        let initial = linked(&law, [300.0, 300.0, 340.0, 280.0, 300.0, 300.0]);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let owner = GasGraphSources::try_new(
            cx,
            &law,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            source_limits(),
            0.1,
        )
        .unwrap();
        let mass = initial.cells()[&CellId(2)].mass_kg;
        let dt = 10_000.0;
        // The fixture has constant original species transport. Its nonlinear
        // NASA enthalpy has the independent analytic inverse in t_from_h.
        let k = law
            .model()
            .evaluate(&point(300.0))
            .unwrap()
            .transport_state()
            .unwrap()
            .thermal_conductivity_w_m_k();
        let conductance = k * AREA / 0.5;
        let mut lower = 0.0;
        let mut upper = 0.5 * mass * (h(340.0) - h(280.0));
        for _ in 0..80 {
            let energy = 0.5 * (lower + upper);
            let hot = t_from_h(h(340.0) - energy / mass);
            let cold = t_from_h(h(280.0) + energy / mass);
            if energy > dt * conductance * (hot - cold) {
                upper = energy;
            } else {
                lower = energy;
            }
        }
        let exchanged = 0.5 * (lower + upper);
        let step = BackwardEulerStep::new(&initial, dt, || Ok(())).unwrap();
        let trial = step
            .evaluate(
                &[0.0, 0.0, -exchanged, exchanged, 0.0, 0.0],
                |candidate| {
                    let resolved = owner.resolve(cx, candidate, &[0.0; 3]).unwrap();
                    let conduction = owner
                        .evaluate_conduction_links(
                            cx,
                            candidate,
                            &resolved.sources,
                            &[GraphConductiveLink {
                                link: 30,
                                branch: 0,
                            }],
                            16,
                        )
                        .unwrap();
                    let mut map = BTreeMap::new();
                    conduction.apply_to(cx, &mut map).unwrap();
                    Ok::<_, DepositionError>(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: map,
                        payload: (resolved.sources, conduction),
                    })
                },
                || Ok(()),
            )
            .unwrap();
        let (next, (sources, conduction), audit) = trial.accept(1e-7, 1e-7, || Ok(())).unwrap();
        near(
            next.temperature_k(CellId(2)).unwrap(),
            t_from_h(h(340.0) - exchanged / mass),
        );
        near(
            next.temperature_k(CellId(3)).unwrap(),
            t_from_h(h(280.0) + exchanged / mass),
        );
        near(dt * conduction.samples()[0].heat_a_to_b_w, exchanged);
        near(conduction.samples()[0].conductance_w_k, conductance);
        assert!(std::ptr::eq(
            &sources.sections()[0][0],
            &conduction.source_sections()[0][0]
        ));
        for (section, id) in [CellId(2), CellId(3)].into_iter().enumerate() {
            let temperature = next.temperature_k(id).unwrap();
            let original = &conduction.source_sections()[0][section].momentum.medium;
            assert_eq!(
                original,
                &law.model().evaluate(&point(temperature)).unwrap()
            );
            near(original.state().pressure, PRESSURE);
            near(next.cells()[&id].mass_kg, initial.cells()[&id].mass_kg);
        }
        assert_eq!(owner.source().mole_composition().fractions(), &[0.75, 0.25]);
        assert!(audit.residual_j.abs() < 1e-7);
        near(audit.external_heat_j, 0.0);
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 340.0);
        near(initial.temperature_k(CellId(3)).unwrap(), 280.0);
    });
}
