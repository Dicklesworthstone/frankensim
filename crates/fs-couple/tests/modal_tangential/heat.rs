//! Physical contact memory and existing caloric cells share one transaction.
use super::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::heat::{
    ModalFrictionHeatError, ModalFrictionHeatSolve, TangentialHeatSource,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::{DryTangentialState, TangentialStep, TangentialStepInput, TriboError};
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 1.0 / 48000.0;
const MASSES: [f64; 2] = [0.002, 0.01];
const CAPACITIES: [f64; 2] = [1e-6, 3e-6];
const K: [f64; 2] = [0.1, 0.2];

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

fn free_system<L: TangentialLaw>(law: L, initial_speed: f64) -> TangentialModalSystem<L> {
    let models = MASSES
        .iter()
        .enumerate()
        .map(|(i, mass)| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                48000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: if i == 0 {
                        initial_speed * mass.sqrt()
                    } else {
                        0.0
                    },
                }])
                .unwrap();
            model
        })
        .collect();
    let network = CoupledModalSystem::new(models, vec![], budgets(), &CancelGate::new()).unwrap();
    TangentialModalSystem::new(
        network,
        attachment(0),
        Some(attachment(1)),
        law,
        config(),
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
    t: f64,
) -> Result<TangentialHeatSource<f64>, DepositionError> {
    Ok(TangentialHeatSource {
        conductivity_w_m_k: K,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: t,
        },
    })
}
fn input(speed: f64, normal: f64) -> TangentialModalDrive {
    TangentialModalDrive {
        speed_m_s: speed,
        normal_force_n: normal,
        contact_temperature_k: None,
        nominal_pressure_pa: None,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn admittances() -> [f64; 2] {
    std::array::from_fn(|i| 2.0 * K[i] * CAPACITIES[i] / (CAPACITIES[i] + DT * 2.0 * K[i]))
}

#[test]
fn mindlin_microslip_heats_bodies_once_and_preserves_reversible_storage() {
    let mut body = free_system(mindlin(), 0.0);
    let mut reference = free_system(mindlin(), 0.0);
    let mut cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let mut work = 0.0;
    let mut loss = 0.0;
    let mut released = false;
    let mut temperature = 300.0;
    for speed in [0.04, 0.04, 0.04, -0.01, -0.04, -0.04, 0.04, 0.04] {
        let before = cells.energy_j();
        let result = body
            .step_with_friction_heat(
                &mut cells,
                &route,
                &[0.0; 2],
                input(speed, 10.0),
                controls(),
                &[0.0; 2],
                temperature,
                source,
                &CancelGate::new(),
            )
            .unwrap();
        let direct = *reference
            .step(
                &[0.0; 2],
                TangentialModalDrive {
                    contact_temperature_k: Some(result.heat.interface_temperature_k),
                    ..input(speed, 10.0)
                },
            )
            .unwrap();
        assert_eq!(result.mechanical, direct);
        assert_eq!(body.state(), reference.state());
        assert_eq!(states(&body), states(&reference));
        assert_eq!(result.source, result.heat.interface_temperature_k);
        close(
            cells.energy_j() - before,
            result.mechanical.friction_dissipation_j,
            1e-12,
        );
        close(
            result.heat.body_heat_j.iter().sum(),
            result.mechanical.friction_dissipation_j,
            1e-12,
        );
        assert!(result.combined_energy_residual_j.abs() < 1e-12);
        work += result.mechanical.drive_work_j;
        loss += result.mechanical.friction_dissipation_j;
        released |= result.mechanical.contact_stored_change_j < -1e-10;
        temperature = result.heat.interface_temperature_k;
        close(
            body.total_energy_j().unwrap() + cells.energy_j(),
            work,
            1e-12,
        );
    }
    assert!(loss > 0.0);
    assert!(body.law().stored_energy_j(body.state()) > 0.0);
    assert!(
        released,
        "a reversal must return stored shear work to the bodies"
    );
    assert!(cells.temperature_k(CellId(0)).unwrap() > 300.0);
    assert!(cells.temperature_k(CellId(1)).unwrap() > 300.0);
}

// A synthetic temperature law delegates mechanics to the original Coulomb
// owner. Its coefficients establish an independent algebraic oracle, not data.
struct SofteningCoulomb;
impl TangentialLaw for SofteningCoulomb {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        "test-affine-thermal-coulomb"
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
        x: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let t = x.contact_temperature_k.ok_or(TriboError::InvalidInput {
            field: "temperature",
        })?;
        FrictionLaw::Coulomb {
            static_mu: 1.0,
            kinetic_mu: 0.8 - 0.05 * (t - 300.0),
        }
        .advance(state, x)
    }
}

#[test]
fn same_step_temperature_changes_friction_and_both_body_momenta() {
    let mut body = free_system(SofteningCoulomb, 1.0);
    let mut cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let before = body.total_energy_j().unwrap();
    let result = body
        .step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 2],
            input(0.0, 1.0),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    // Two free masses and two finite BE cells: H=mu*(dt-c*mu),
    // H=dt*sum(Y_i)*(T-300), mu=0.8-0.05*(T-300).
    let c = 0.5 * DT * DT * (MASSES[0].recip() + MASSES[1].recip());
    let y = admittances();
    let a = c * 0.05 * 0.05;
    let b = DT * (y[0] + y[1] + 0.05) - 2.0 * c * 0.8 * 0.05;
    let d = 0.8 * DT - c * 0.8 * 0.8;
    let rise = 2.0 * d / (b + (b * b + 4.0 * a * d).sqrt());
    let mu = 0.8 - 0.05 * rise;
    close(result.heat.interface_temperature_k, 300.0 + rise, 2e-8);
    close(result.mechanical.traction_n, -mu, 2e-9);
    assert!(mu < 0.7);
    assert!(result.trials > 3);
    for i in 0..2 {
        let force = if i == 0 { -mu } else { mu };
        let initial = if i == 0 { 1.0 } else { 0.0 };
        let state = body.network().components()[i].states()[0];
        close(
            state.velocity_m_sqrt_kg_per_s / MASSES[i].sqrt(),
            initial + force * DT / MASSES[i],
            2e-10,
        );
        close(result.heat.body_heat_j[i], DT * y[i] * rise, 1e-12);
    }
    close(
        body.total_energy_j().unwrap() + cells.energy_j(),
        before,
        1e-12,
    );
    assert_eq!(body.state().contact_temperature_k(), Some(result.source));
}

#[test]
fn ordinary_conduction_survives_zero_friction_without_artificial_partition() {
    let mut body = free_system(
        FrictionLaw::Coulomb {
            static_mu: 0.0,
            kinetic_mu: 0.0,
        },
        0.0,
    );
    let mut cells = thermal([310.0, 290.0]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let before = cells.energy_j();
    let result = body
        .step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 2],
            input(0.0, 1.0),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let y = admittances();
    let expected = (y[0] * 310.0 + y[1] * 290.0) / (y[0] + y[1]);
    close(result.heat.interface_temperature_k, expected, 2e-8);
    assert_eq!(result.mechanical.friction_dissipation_j, 0.0);
    assert!(result.heat.body_heat_j[0] < 0.0);
    assert!(result.heat.body_heat_j[1] > 0.0);
    close(cells.energy_j(), before, 1e-12);
}

#[test]
fn loaded_history_and_thermal_state_survive_late_refusal_and_cancellation() {
    let mut body = free_system(mindlin(), 0.0);
    let mut reference = free_system(mindlin(), 0.0);
    let mut cells = thermal([300.0; 2]);
    let mut reference_cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    for (s, t) in [
        (&mut body, &mut cells),
        (&mut reference, &mut reference_cells),
    ] {
        s.step_with_friction_heat(
            t,
            &route,
            &[0.0; 2],
            input(0.04, 10.0),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    }
    let history = body.state().clone();
    let before = states(&body);
    let frame = body.last_frame().copied();
    let heat = cells.cells().clone();
    let input = input(-0.01, 10.0);
    let mut exhausted = controls();
    exhausted.maximum_trials = 1;
    assert!(matches!(
        body.step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 2],
            input,
            exhausted,
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new()
        ),
        Err(ModalFrictionHeatError::TrialBudget)
    ));
    let calls = Cell::new(0usize);
    let refuse = |t: &DepositionSimulation<'_>, temperature| {
        calls.set(calls.get() + 1);
        if calls.get() == 3 {
            Err(DepositionError::MaterialDomain)
        } else {
            source(t, temperature)
        }
    };
    assert!(matches!(
        body.step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 2],
            input,
            controls(),
            &[0.0; 2],
            300.0,
            refuse,
            &CancelGate::new()
        ),
        Err(ModalFrictionHeatError::Source(
            DepositionError::MaterialDomain
        ))
    ));
    let gate = CancelGate::new();
    let calls = Cell::new(0usize);
    let cancel = |t: &DepositionSimulation<'_>, temperature| {
        calls.set(calls.get() + 1);
        if calls.get() == 2 {
            gate.request();
        }
        source(t, temperature)
    };
    assert!(
        body.step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 2],
            input,
            controls(),
            &[0.0; 2],
            300.0,
            cancel,
            &gate
        )
        .is_err()
    );
    assert_eq!(body.state(), &history);
    assert_eq!(states(&body), before);
    assert_eq!(body.last_frame().copied(), frame);
    assert_eq!(cells.cells(), &heat);
    assert_eq!(cells.time_s(), DT);
    let actual = body
        .step_with_friction_heat(
            &mut cells,
            &route,
            &[0.0; 2],
            input,
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
            &[0.0; 2],
            input,
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(actual.mechanical, expected.mechanical);
    assert_eq!(body.state(), reference.state());
    assert_eq!(states(&body), states(&reference));
    assert_eq!(cells.cells(), reference_cells.cells());
}

#[test]
fn frozen_temperature_law_is_not_replaced_or_reset_to_force_thermal_acceptance() {
    let mut body = free_system(layer(0.01), 0.0);
    let mut cells = thermal([300.0; 2]);
    let route = FrictionHeatRoute::from_link(&cells, 7, K).unwrap();
    let original = body.state().clone();
    let law = body.law().clone();
    let result = body.step_with_friction_heat(
        &mut cells,
        &route,
        &[0.0; 2],
        input(0.2, 1.0),
        controls(),
        &[0.0; 2],
        300.0,
        source,
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(ModalFrictionHeatError::Mechanical(_))));
    assert_eq!(body.state(), &original);
    assert_eq!(body.law(), &law);
    assert_eq!(body.network().samples_rendered(), 0);
    assert_eq!(cells.time_s(), 0.0);
    assert_eq!(cells.energy_j(), 0.0);
    assert!(body.last_frame().is_none());
}
