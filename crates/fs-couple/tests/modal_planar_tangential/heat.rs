//! Same-step thermal feedback through the actual matrix modal/enthalpy owners.
use super::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::planar::heat::{
    ModalFrictionHeatSolve, PlanarTangentialHeatError, TangentialHeatSource,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{DryTangentialState, TangentialStep, TangentialStepInput, TriboError};
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 1.0 / 48000.0;
const CAPACITIES: [f64; 2] = [1e-6, 3e-6];
const K: [f64; 2] = [0.1, 0.2];
const MASSES: [[f64; 2]; 2] = [[0.002, 0.008], [0.01, 0.02]];
const ROTATION: [[f64; 2]; 2] = [[0.8, -0.6], [0.6, 0.8]];

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + h;
        if temperature.is_finite() && temperature > 0.0 {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static CALORIC: Caloric = Caloric;

fn thermal(temperatures: [f64; 2]) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..2 {
        let x = i as f64 - 1.0;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: CAPACITIES[i],
                specific_h_j_kg: temperatures[i] - 300.0,
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
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(0),
            b: CellId(1),
            area_m2: 1.0,
            conductance_w_k: series_conductance_w_k(1.0, [0.5; 2], K).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn free_system<L: TangentialLaw>(law: L, initial_speed: f64) -> PlanarTangentialModalSystem<L> {
    let models = MASSES
        .iter()
        .enumerate()
        .map(|(body, masses)| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                48000,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    };
                    2
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&std::array::from_fn::<_, 2, _>(|mode| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: if body == 0 && mode == 0 {
                        initial_speed * masses[0].sqrt()
                    } else {
                        0.0
                    },
                }))
                .unwrap();
            model
        })
        .collect();
    let network = CoupledModalSystem::new(
        models,
        vec![],
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 4,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
        },
        &CancelGate::new(),
    )
    .unwrap();
    let attachment = |body| {
        std::array::from_fn(|axis| ModalAttachment {
            component: body,
            shapes: (0..2)
                .map(|mode| ROTATION[axis][mode] / MASSES[body][mode].sqrt())
                .collect(),
        })
    };
    PlanarTangentialModalSystem::new(
        network,
        attachment(0),
        Some(attachment(1)),
        law,
        config(),
        matrix_config(),
        &CancelGate::new(),
    )
    .unwrap()
}

fn controls() -> ModalFrictionHeatSolve {
    ModalFrictionHeatSolve {
        maximum_cells: 2,
        maximum_trials: 256,
        maximum_newton_iterations: 16,
        difference_step_j: 1e-8,
        temperature_scale_j_per_k: 1e-6,
        heat: FrictionHeatTolerances {
            cell_j: 1e-12,
            thermal_j: 1e-12,
            interface_j: 1e-12,
            total_j: 1e-12,
        },
        combined_energy_tolerance_j: 1e-12,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-10,
            linear_restart: 3,
            ..Default::default()
        },
    }
}

fn source(
    _: &DepositionSimulation<'_>,
    temperature: f64,
) -> Result<TangentialHeatSource<f64>, DepositionError> {
    Ok(TangentialHeatSource {
        conductivity_w_m_k: K,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: temperature,
        },
    })
}

fn input(speed: [f64; 2]) -> PlanarTangentialModalDrive {
    PlanarTangentialModalDrive {
        contact_temperature_k: None,
        ..drive(speed)
    }
}

fn admittances() -> [f64; 2] {
    std::array::from_fn(|i| 2.0 * K[i] * CAPACITIES[i] / (CAPACITIES[i] + DT * 2.0 * K[i]))
}

// Synthetic thermal dependence only: mechanics remains the original dry law.
// These coefficients define an analytic reference, not material calibration.
struct SofteningCoulomb;
impl TangentialLaw for SofteningCoulomb {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        "test-planar-affine-thermal-coulomb"
    }
    fn initial_state(&self) -> Self::State {
        FrictionLaw::Coulomb {
            static_mu: 1.0,
            kinetic_mu: 0.8,
        }
        .initial_state()
    }
    fn stored_energy_j(&self, _: &Self::State) -> f64 {
        0.0
    }
    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let temperature = input
            .contact_temperature_k
            .ok_or(TriboError::InvalidInput {
                field: "temperature",
            })?;
        FrictionLaw::Coulomb {
            static_mu: 1.0,
            kinetic_mu: 0.8 - 0.05 * (temperature - 300.0),
        }
        .advance(state, input)
    }
}

#[test]
fn solved_temperature_changes_both_force_axes_and_both_body_recoils() {
    let mut owner = free_system(SofteningCoulomb, 1.0);
    let mut cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let before = owner.total_energy_j().unwrap();
    let compliance = owner.compliance_m_per_n();
    assert!(compliance[0][1].abs() > 0.2 * compliance[0][0]);
    let result = owner
        .step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 4],
            input([0.0; 2]),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    // The physical eigenaxis is rotated relative to both supplied tangents.
    // H=mu*(dt-c*mu)=dt*sum(Y_i)*rise, mu=0.8-0.05*rise.
    let c = 0.5 * DT * DT * (MASSES[0][0].recip() + MASSES[1][0].recip());
    let y = admittances();
    let a = c * 0.05 * 0.05;
    let b = DT * (y[0] + y[1] + 0.05) - 2.0 * c * 0.8 * 0.05;
    let d = 0.8 * DT - c * 0.8 * 0.8;
    let rise = 2.0 * d / (b + (b * b + 4.0 * a * d).sqrt());
    let mu = 0.8 - 0.05 * rise;
    close(result.heat.interface_temperature_k, 300.0 + rise, 2e-8);
    assert!(mu < 0.7);
    assert!(result.trials > 3);
    for axis in 0..2 {
        close(
            result.mechanical.traction_n[axis],
            -ROTATION[axis][0] * mu,
            2e-9,
        );
        close(
            result.mechanical.slip_distance_m[axis],
            -ROTATION[axis][0] * (DT - c * mu),
            1e-13,
        );
        for body in 0..2 {
            let initial = if body == 0 { 1.0 } else { 0.0 };
            let force = if body == 0 { -mu } else { mu };
            let actual: f64 = (0..2)
                .map(|mode| {
                    ROTATION[axis][mode]
                        * owner.network().components()[body].states()[mode].velocity_m_sqrt_kg_per_s
                        / MASSES[body][mode].sqrt()
                })
                .sum();
            close(
                actual,
                ROTATION[axis][0] * (initial + force * DT / MASSES[body][0]),
                2e-10,
            );
        }
    }
    for body in 0..2 {
        close(result.heat.body_heat_j[body], DT * y[body] * rise, 1e-12);
    }
    close(
        cells.energy_j(),
        result.mechanical.friction_dissipation_j,
        1e-12,
    );
    close(
        owner.total_energy_j().unwrap() + cells.energy_j(),
        before,
        1e-12,
    );
    assert_eq!(owner.state().contact_temperature_k(), Some(result.source));
    assert!(result.combined_energy_residual_j.abs() < 1e-12);
}

#[test]
fn transverse_elastic_storage_and_return_work_never_become_heat() {
    let mut owner = system(layer(100.0), true, IDENTITY);
    let mut direct = system(layer(100.0), true, IDENTITY);
    let mut cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let mut work = 0.0;
    let mut network_loss = 0.0;
    let mut released = false;
    for speed in [[0.01, 0.02], [0.02, -0.01], [-0.02, -0.02], [-0.01, 0.005]] {
        let result = owner
            .step_with_friction_heat(
                &mut cells,
                &route,
                &[0.0; 4],
                input(speed),
                controls(),
                &[0.0; 2],
                300.0,
                source,
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(
            result.mechanical,
            *direct.step(&[0.0; 4], drive(speed)).unwrap()
        );
        assert_eq!(owner.state(), direct.state());
        assert_eq!(states(&owner), states(&direct));
        let elastic = owner.state().elastic_slip_m();
        assert!(elastic[1].abs() > 1e-10);
        close(
            result.mechanical.contact_energy_j,
            1e6 * (elastic[0] * elastic[0] + elastic[1] * elastic[1]),
            1e-15,
        );
        assert_eq!(result.mechanical.friction_dissipation_j, 0.0);
        assert_eq!(result.heat.body_heat_j, [0.0; 2]);
        assert_eq!(cells.energy_j(), 0.0);
        work += result.mechanical.drive_work_j;
        network_loss += result.mechanical.network_dissipation_j;
        released |= result.mechanical.contact_stored_change_j < -1e-10;
        close(
            owner.total_energy_j().unwrap() + network_loss + cells.energy_j(),
            work,
            1e-12,
        );
    }
    assert!(released);
}

#[test]
fn hot_to_cold_conduction_survives_exact_zero_slip_in_both_axes() {
    let mut owner = free_system(
        FrictionLaw::Coulomb {
            static_mu: 1.0,
            kinetic_mu: 0.8,
        },
        0.0,
    );
    let mut cells = thermal([310.0, 290.0]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let before = cells.energy_j();
    let result = owner
        .step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 4],
            input([0.0; 2]),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let y = admittances();
    close(
        result.heat.interface_temperature_k,
        (y[0] * 310.0 + y[1] * 290.0) / (y[0] + y[1]),
        2e-8,
    );
    assert_eq!(result.mechanical.slip_distance_m, [0.0; 2]);
    assert_eq!(result.mechanical.traction_n, [0.0; 2]);
    assert_eq!(result.mechanical.friction_dissipation_j, 0.0);
    assert!(result.heat.body_heat_j[0] < 0.0);
    assert!(result.heat.body_heat_j[1] > 0.0);
    close(cells.energy_j(), before, 1e-12);
}

#[test]
fn loaded_vector_history_and_heat_survive_late_refusal_cancellation_and_retry() {
    let mut owner = system(layer(100.0), true, IDENTITY);
    let mut reference = system(layer(100.0), true, IDENTITY);
    let mut cells = thermal([300.0; 2]);
    let mut reference_cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    for (body, thermal) in [
        (&mut owner, &mut cells),
        (&mut reference, &mut reference_cells),
    ] {
        body.step_with_friction_heat(
            thermal,
            &route,
            &[0.0; 4],
            input([0.01, 0.02]),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    }
    let history = owner.state().clone();
    assert!(history.elastic_slip_m().iter().all(|x| x.abs() > 1e-10));
    let bodies = states(&owner);
    let last = owner.last_frame().copied();
    let caloric = cells.cells().clone();
    let next = input([-0.01, 0.005]);
    let mut exhausted = controls();
    exhausted.maximum_trials = 1;
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 4],
            next,
            exhausted,
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        ),
        Err(PlanarTangentialHeatError::TrialBudget)
    ));
    let calls = Cell::new(0usize);
    let refuse = |thermal: &DepositionSimulation<'_>, temperature| {
        calls.set(calls.get() + 1);
        if calls.get() == 2 {
            Err(DepositionError::MaterialDomain)
        } else {
            source(thermal, temperature)
        }
    };
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 4],
            next,
            controls(),
            &[0.0; 2],
            300.0,
            refuse,
            &CancelGate::new(),
        ),
        Err(PlanarTangentialHeatError::Source(
            DepositionError::MaterialDomain
        ))
    ));
    let gate = CancelGate::new();
    let calls = Cell::new(0usize);
    let cancel = |thermal: &DepositionSimulation<'_>, temperature| {
        calls.set(calls.get() + 1);
        if calls.get() == 2 {
            gate.request();
        }
        source(thermal, temperature)
    };
    assert!(
        owner
            .step_with_friction_heat(
                &mut cells,
                &route,
                &[0.0; 4],
                next,
                controls(),
                &[0.0; 2],
                300.0,
                cancel,
                &gate,
            )
            .is_err()
    );
    assert_eq!(owner.state(), &history);
    assert_eq!(states(&owner), bodies);
    assert_eq!(owner.last_frame().copied(), last);
    assert_eq!(owner.network().samples_rendered(), 1);
    assert_eq!(cells.cells(), &caloric);
    assert_eq!(cells.time_s(), DT);
    let actual = owner
        .step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 4],
            next,
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let expected = reference
        .step_with_friction_heat(
            &mut reference_cells,
            &route,
            &[0.0; 4],
            next,
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(actual.mechanical, expected.mechanical);
    assert_eq!(owner.state(), reference.state());
    assert_eq!(states(&owner), states(&reference));
    assert_eq!(cells.cells(), reference_cells.cells());
}

#[test]
fn frozen_elastic_source_refuses_heating_without_resetting_or_rebinding() {
    let mut owner = system(layer(0.01), false, IDENTITY);
    let mut cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let history = owner.state().clone();
    let law = owner.law().clone();
    assert!(matches!(
        owner.step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 4],
            input([0.2, 0.1]),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        ),
        Err(PlanarTangentialHeatError::Mechanical(_))
    ));
    assert_eq!(owner.state(), &history);
    assert_eq!(owner.law(), &law);
    assert_eq!(owner.network().samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
    assert_eq!(cells.time_s(), 0.0);
    assert_eq!(cells.energy_j(), 0.0);
}
