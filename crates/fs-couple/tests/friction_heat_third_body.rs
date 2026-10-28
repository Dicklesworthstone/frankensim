//! Finite-cell third-body physics, with analytical synthetic material fixtures.
use fs_couple::friction_heat::tangential_heat_candidate;
use fs_couple::friction_heat::third_body::*;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::{FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use std::cell::Cell;
use std::collections::BTreeMap;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static CALORIC: Caloric = Caloric;
const LAYER: CellId = CellId(2);
const LINKS: [u64; 2] = [7, 8];
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn tolerances() -> ThirdBodyHeatTolerances {
    ThirdBodyHeatTolerances {
        cell_j: 1e-12,
        thermal_j: 1e-12,
        layer_j: 1e-12,
        total_j: 1e-12,
    }
}
fn network(
    capacity: [f64; 3],
    conductivity: [f64; 3],
    resistance: [f64; 2],
    layer_source: f64,
) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..3 {
        let x = [-1.0, 1.0, 0.0][i];
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: capacity[i],
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [x, 0.0, 0.0],
                    end_m: [x + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: if i == 2 { layer_source } else { 0.0 },
                },
            }),
        });
    }
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: CellId(i as u64),
                b: LAYER,
                area_m2: 1.0,
                conductance_w_k: series_contact_conductance_w_k(
                    1.0,
                    [0.5; 2],
                    [conductivity[i], conductivity[2]],
                    resistance[i],
                )
                .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 4,
            max_events: 8,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.01,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn source<P>(work: f64, k: [f64; 3], resistance: [f64; 2], payload: P) -> ThirdBodyHeatInput<P> {
    ThirdBodyHeatInput {
        dissipated_j: work,
        conductivity_w_m_k: k,
        area_resistance_m2_k_w: resistance,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload,
        },
    }
}
fn input(dt: f64, du: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 10.0,
        slip_increment_m: [du, 0.0],
        dt_s: dt,
        contact_temperature_k: None,
        nominal_pressure_pa: Some(10.0),
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
    }
}
fn endpoint(capacity: [f64; 3], g: [f64; 2], dt: f64, work: f64) -> [f64; 3] {
    let weight: [f64; 2] = std::array::from_fn(|i| dt * g[i] / (capacity[i] + dt * g[i]));
    let rise = work / (capacity[2] + capacity[0] * weight[0] + capacity[1] * weight[1]);
    [
        capacity[0] * weight[0] * rise,
        capacity[1] * weight[1] * rise,
        capacity[2] * rise,
    ]
}

#[test]
fn g0_finite_layer_capacity_and_contact_resistance_close_original_work_once() {
    let capacity = [2.0, 3.0, 0.25];
    let k = [4.0, 1.0, 2.0];
    let resistance = [0.2, 0.8];
    let thermal = network(capacity, k, resistance, 0.3);
    let original = thermal.cells().clone();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, k, resistance).unwrap();
    assert_eq!(route.body_cells(), [CellId(0), CellId(1)]);
    assert_eq!(route.distances_m(), [[0.5; 2]; 2]);
    let g = LINKS.map(|id| thermal.links()[&id].conductance_w_k);
    let step = ThirdBodyHeatStep::new(&thermal, 0.1, &route, || Ok(())).unwrap();
    let energy = endpoint(capacity, g, 0.1, 0.4 + 0.03);
    let trial = step
        .evaluate(
            &energy,
            |_, _| Ok::<_, DepositionError>(source(0.4, k, resistance, 17)),
            || Ok(()),
        )
        .unwrap();
    let (accepted, payload, report) = trial.accept(tolerances(), || Ok(())).unwrap();
    assert_eq!(payload, 17);
    close(report.layer_stored_change_j, energy[2], 1e-14);
    close(report.layer_external_heat_j, 0.03, 1e-14);
    close(report.external_heat_j, 0.03, 1e-14);
    close(
        report.layer_stored_change_j + report.body_heat_j.iter().sum::<f64>(),
        0.43,
        1e-13,
    );
    assert!(report.layer_stored_change_j > 0.0 && report.layer_stored_change_j < 0.4);
    for i in 0..2 {
        close(report.body_heat_j[i], energy[i], 1e-13);
    }
    close(accepted.energy_j(), 0.43, 1e-13);
    assert_eq!(thermal.cells(), &original);
    assert_eq!(thermal.time_s(), 0.0);
    // Original contact coefficients remain authored topology, not a trial side effect.
    assert_eq!(thermal.links(), accepted.links());
}

#[test]
fn g1_same_endpoint_temperature_changes_the_original_tangential_law() {
    let capacity = [2.0, 3.0, 0.25];
    let k = [4.0, 1.0, 2.0];
    let resistance = [0.2, 0.8];
    let thermal = network(capacity, k, resistance, 0.0);
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, k, resistance).unwrap();
    let g = LINKS.map(|id| thermal.links()[&id].conductance_w_k);
    let dt = 0.1;
    let beta = 0.1;
    let unit = endpoint(capacity, g, dt, 1.0);
    // N=10, v=1: W = (0.8-beta*theta)*10*dt, theta=W*unit_layer/C_layer.
    let work = 0.8 * 10.0 * dt / (1.0 + beta * 10.0 * dt * unit[2] / capacity[2]);
    let energy = endpoint(capacity, g, dt, work);
    let cold = FrictionLaw::Coulomb {
        static_mu: 0.9,
        kinetic_mu: 0.8,
    };
    let previous = cold.initial_state();
    let step = ThirdBodyHeatStep::new(&thermal, dt, &route, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &energy,
            |_, t| {
                let law = FrictionLaw::Coulomb {
                    static_mu: 0.9,
                    kinetic_mu: 0.8 - beta * (t - 300.0),
                };
                let mechanical =
                    tangential_heat_candidate(&law, &previous, input(dt, dt), t, 1e-12)
                        .map_err(|_| DepositionError::MaterialDomain)?;
                Ok::<_, DepositionError>(source(mechanical.dissipated_j, k, resistance, mechanical))
            },
            || Ok(()),
        )
        .unwrap();
    let (next, mechanical, report) = trial.accept(tolerances(), || Ok(())).unwrap();
    close(report.dissipated_j, work, 1e-13);
    close(mechanical.traction_n[0], work / dt, 1e-12);
    close(
        report.layer_temperature_k,
        300.0 + energy[2] / capacity[2],
        1e-13,
    );
    close(
        next.temperature_k(LAYER).unwrap(),
        report.layer_temperature_k,
        1e-13,
    );
    assert!(report.dissipated_j < 0.8);
    assert_eq!(previous, cold.initial_state());
}

#[test]
fn g0_loaded_shear_layer_preserves_elastic_storage_and_heats_only_plastic_work() {
    let thermal = network([1.0; 3], [1.0; 3], [0.0; 2], 0.0);
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, [1.0; 3], [0.0; 2]).unwrap();
    let law_at = |t| ThirdBodyLayerContact::new(100.0, 1.0, 0.5, t, 1.0).unwrap();
    let cold = law_at(300.0);
    let previous = cold.initial_state();
    let base = tangential_heat_candidate(&cold, &previous, input(0.1, 0.1), 300.0, 1e-12).unwrap();
    assert!(base.stored_change_j > 0.0 && base.dissipated_j > 0.0);
    let energy = endpoint([1.0; 3], [1.0; 2], 0.1, base.dissipated_j);
    let step = ThirdBodyHeatStep::new(&thermal, 0.1, &route, || Ok(())).unwrap();
    let trial = step
        .evaluate(
            &energy,
            |_, t| {
                let law = law_at(t);
                let retained = law
                    .continue_state(&previous)
                    .map_err(|_| DepositionError::MaterialDomain)?;
                let mechanical =
                    tangential_heat_candidate(&law, &retained, input(0.1, 0.1), t, 1e-12)
                        .map_err(|_| DepositionError::MaterialDomain)?;
                Ok::<_, DepositionError>(source(
                    mechanical.dissipated_j,
                    [1.0; 3],
                    [0.0; 2],
                    mechanical,
                ))
            },
            || Ok(()),
        )
        .unwrap();
    let (next, mechanical, report) = trial.accept(tolerances(), || Ok(())).unwrap();
    close(next.energy_j(), base.dissipated_j, 1e-12);
    close(
        report.dissipated_j + mechanical.stored_change_j,
        mechanical.traction_n[0] * 0.1,
        1e-12,
    );
    assert!(next.energy_j() < mechanical.traction_n[0] * 0.1);
    assert_eq!(previous, cold.initial_state());
}

#[test]
fn g0_insulated_finite_layer_stores_heat_and_refusals_publish_nothing() {
    let thermal = network([1.0; 3], [0.0; 3], [0.0; 2], 0.0);
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, [0.0; 3], [0.0; 2]).unwrap();
    let step = ThirdBodyHeatStep::new(&thermal, 0.1, &route, || Ok(())).unwrap();
    let accepted = step
        .evaluate(
            &[0.0, 0.0, 0.5],
            |_, _| Ok::<_, DepositionError>(source(0.5, [0.0; 3], [0.0; 2], ())),
            || Ok(()),
        )
        .unwrap()
        .accept(tolerances(), || Ok(()))
        .unwrap();
    close(accepted.2.layer_stored_change_j, 0.5, 1e-14);
    assert_eq!(accepted.2.body_heat_j, [0.0; 2]);
    let bad = step
        .evaluate(
            &[0.0; 3],
            |_, _| Ok::<_, DepositionError>(source(0.5, [0.0; 3], [0.0; 2], ())),
            || Ok(()),
        )
        .unwrap();
    assert!(bad.accept(tolerances(), || Ok(())).is_err());
    let duplicate = step.evaluate(
        &[0.0; 3],
        |_, _| {
            let mut s = source(0.0, [0.0; 3], [0.0; 2], ());
            s.transport.conductances_w_k.insert(LINKS[0], 1.0);
            Ok::<_, DepositionError>(s)
        },
        || Ok(()),
    );
    assert!(duplicate.is_err());
    let polls = Cell::new(0);
    let cancelled = step.evaluate(
        &[0.0, 0.0, 0.5],
        |_, _| Ok::<_, DepositionError>(source(0.5, [0.0; 3], [0.0; 2], ())),
        || {
            polls.set(polls.get() + 1);
            if polls.get() >= 5 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(cancelled, Err(DepositionError::Cancelled)));
    assert_eq!(thermal.energy_j(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    let extra = thermal
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: 9,
                expected: None,
                replacement: Some(ContactPatch::symmetric(ThermalLink {
                    id: 9,
                    a: LAYER,
                    b: CellId(0),
                    area_m2: 0.1,
                    conductance_w_k: 0.0,
                })),
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(ThirdBodyHeatStep::new(&extra, 0.1, &route, || Ok(())).is_err());
    assert!(ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, [1.0; 3], [0.0; 2]).is_err());
}

#[test]
fn g1_finite_layer_transient_converges_on_three_timestep_rungs() {
    // Equal C=1 J/K and G=1 W/K, P=1 W. The layer/body difference has
    // d'=1-3d, and the energy mean rises t/3. This is an independent exact ODE.
    let exact = 1.0 / 3.0 + 2.0 / 9.0 * (1.0 - (-3.0_f64).exp());
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut thermal = network([1.0; 3], [1.0; 3], [0.0; 2], 0.0);
        let route =
            ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, [1.0; 3], [0.0; 2]).unwrap();
        let dt = 1.0 / n as f64;
        let mut difference = 0.0;
        for i in 0..n {
            difference = (difference + dt) / (1.0 + 3.0 * dt);
            let mean = (i + 1) as f64 * dt / 3.0;
            let rise = [
                mean - difference / 3.0,
                mean - difference / 3.0,
                mean + 2.0 * difference / 3.0,
            ];
            let energies: [f64; 3] = std::array::from_fn(|j| {
                rise[j] - thermal.cells()[&CellId(j as u64)].specific_h_j_kg
            });
            let step =
                ThirdBodyHeatStep::new(&thermal, (i + 1) as f64 * dt, &route, || Ok(())).unwrap();
            thermal = step
                .evaluate(
                    &energies,
                    |_, _| Ok::<_, DepositionError>(source(dt, [1.0; 3], [0.0; 2], ())),
                    || Ok(()),
                )
                .unwrap()
                .accept(tolerances(), || Ok(()))
                .unwrap()
                .0;
        }
        close(thermal.energy_j(), 1.0, 1e-12);
        let error = (thermal.cells()[&LAYER].specific_h_j_kg - exact).abs();
        eprintln!("third-body steps={n} layer_rise_error_k={error:.8e}");
        errors.push(error);
    }
    assert!(errors[0] / errors[1] > 1.9 && errors[1] / errors[2] > 1.9);
    assert!(errors[2] < 0.002);
}
