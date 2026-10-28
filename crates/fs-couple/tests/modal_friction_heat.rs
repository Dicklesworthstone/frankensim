// G1/G3: actual two-body motion and a finite thermal network accept together.
// Caloric/conductivity coefficients are analytical fixtures, not material data.
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::heat::{
    ModalFrictionHeatError, ModalFrictionHeatSolve, ModalHeatSource,
};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort, ModalFrictionRegime,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::FrictionLaw;
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 1.0 / 48000.0;
const MASSES: [f64; 2] = [0.002, 0.01];
const CAPACITIES: [f64; 2] = [1e-6, 3e-6];
const K: [f64; 2] = [0.1, 0.2];
const COLD_MU: f64 = 0.8;
const MU_SLOPE: f64 = 0.05;

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

fn law(t: f64) -> FrictionLaw {
    FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: COLD_MU - MU_SLOPE * (t - 300.0),
    }
}
fn mechanics() -> FrictionModalSystem {
    mechanics_with_law(law(300.0))
}
fn mechanics_with_law(law: FrictionLaw) -> FrictionModalSystem {
    let images = MASSES
        .iter()
        .enumerate()
        .map(|(i, &mass)| {
            let mut m = ModalAcousticTimeModel::try_new_with_free_coordinates(
                48000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            m.restore_states(&[ModalAcousticState {
                displacement_m_sqrt_kg: 0.0,
                velocity_m_sqrt_kg_per_s: if i == 0 { mass.sqrt() } else { 0.0 },
            }])
            .unwrap();
            m
        })
        .collect();
    let network = CoupledModalSystem::new(
        images,
        vec![],
        ModalCouplingConfig {
            max_modes: 2,
            max_connections: 1,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 100.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-9,
        },
        &CancelGate::new(),
    )
    .unwrap();
    let attachment = |i: usize| ModalAttachment {
        component: i,
        shapes: vec![MASSES[i].sqrt().recip()],
    };
    FrictionModalSystem::new(
        network,
        ModalFrictionPort {
            left: attachment(0),
            right: Some(attachment(1)),
            law,
        },
        ModalFrictionConfig {
            max_iterations: 80,
            maximum_force_n: 100.0,
            maximum_slip_speed_m_s: 10.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-11,
            sticking_tolerance_m: 1e-12,
        },
        &CancelGate::new(),
    )
    .unwrap()
}
fn thermal() -> DepositionSimulation<'static> {
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
fn source(_: &DepositionSimulation<'_>, t: f64) -> Result<ModalHeatSource<f64>, DepositionError> {
    Ok(ModalHeatSource {
        law: law(t),
        conductivity_w_m_k: K,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: t,
        },
    })
}
fn drive() -> FrictionDrive {
    FrictionDrive {
        speed_m_s: 0.0,
        normal_force_n: 1.0,
    }
}
fn states(s: &FrictionModalSystem) -> Vec<ModalAcousticState> {
    s.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

#[test]
fn g1_temperature_feedback_matches_independent_two_body_and_caloric_solution() {
    let mut mechanical = mechanics();
    let mut thermal = thermal();
    let route = FrictionHeatRoute::from_link(&thermal, 7, K).unwrap();
    let before = mechanical.total_energy_j().unwrap();
    let result = mechanical
        .step_with_friction_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    // Independent reduction of two free masses and two backward-Euler cells:
    // H = mu*(dt - compliance*mu) = dt*(Y_a+Y_b)*(T_interface-300).
    let compliance = 0.5 * DT * DT * (MASSES[0].recip() + MASSES[1].recip());
    let y = std::array::from_fn::<_, 2, _>(|i| {
        2.0 * K[i] * CAPACITIES[i] / (CAPACITIES[i] + DT * 2.0 * K[i])
    });
    let a = compliance * MU_SLOPE * MU_SLOPE;
    let b = DT * (y[0] + y[1] + MU_SLOPE) - 2.0 * compliance * COLD_MU * MU_SLOPE;
    let c = COLD_MU * DT - compliance * COLD_MU * COLD_MU;
    let rise = 2.0 * c / (b + (b * b + 4.0 * a * c).sqrt());
    let mu = COLD_MU - MU_SLOPE * rise;
    close(result.heat.interface_temperature_k, 300.0 + rise, 2e-8);
    close(result.mechanical.traction_n, -mu, 2e-9);
    assert_eq!(result.mechanical.regime, ModalFrictionRegime::Sliding);
    assert_eq!(mechanical.port().law, law(result.source));
    assert_eq!(result.source, result.heat.interface_temperature_k);
    assert!(
        mu < COLD_MU - 0.1,
        "same-step thermal softening must affect motion"
    );
    for i in 0..2 {
        let force = if i == 0 { -mu } else { mu };
        let initial_v = if i == 0 { 1.0 } else { 0.0 };
        let state = mechanical.components()[i].states()[0];
        close(
            state.velocity_m_sqrt_kg_per_s / MASSES[i].sqrt(),
            initial_v + force * DT / MASSES[i],
            2e-10,
        );
        close(
            state.displacement_m_sqrt_kg / MASSES[i].sqrt(),
            initial_v * DT + 0.5 * force * DT * DT / MASSES[i],
            2e-14,
        );
        close(result.heat.body_heat_j[i], DT * y[i] * rise, 1e-12);
    }
    close(
        thermal.energy_j(),
        result.mechanical.friction_dissipation_j,
        1e-12,
    );
    close(
        mechanical.total_energy_j().unwrap() + thermal.energy_j(),
        before,
        1e-12,
    );
    assert!(result.combined_energy_residual_j.abs() < 1e-12);
    assert_eq!(mechanical.samples_rendered(), 1);
    assert_eq!(thermal.time_s(), DT);
}

#[test]
fn g3_stribeck_load_line_uses_each_trial_temperature_and_keeps_the_original_branch_rule() {
    let stribeck = |t: f64| {
        let factor = 1.0 - 0.04 * (t - 300.0);
        FrictionLaw::Stribeck {
            static_mu: 0.9 * factor,
            kinetic_mu: 0.6 * factor,
            characteristic_speed: 0.3,
            viscous_per_speed: 0.0,
        }
    };
    let mut mechanical = mechanics_with_law(stribeck(300.0));
    let mut thermal = thermal();
    let route = FrictionHeatRoute::from_link(&thermal, 7, K).unwrap();
    let resolve = |endpoint: &DepositionSimulation<'_>, t| {
        let mut data = source(endpoint, t)?;
        data.law = stribeck(t);
        Ok::<_, DepositionError>(data)
    };
    let result = mechanical
        .step_with_friction_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 2],
            300.0,
            resolve,
            &CancelGate::new(),
        )
        .unwrap();
    assert!(
        result.trials > 3,
        "temperature must change over multiple original mechanical trials"
    );
    let accepted_law = stribeck(result.heat.interface_temperature_k);
    let speed = (result.mechanical.slip_distance_m / DT).abs();
    close(
        -result.mechanical.traction_n,
        accepted_law.kinetic_coefficient(speed).unwrap(),
        1e-10,
    );
    assert!(result.mechanical.traction_n.abs() < 0.55);
    // A fresh original mechanical owner at the solved temperature must choose
    // the identical root and trajectory. Reusing a cold load line fails this.
    let mut independent = mechanics_with_law(accepted_law);
    let expected = *independent.step(&[0.0; 2], drive()).unwrap();
    assert_eq!(result.mechanical, expected);
    assert_eq!(states(&mechanical), states(&independent));
    assert_eq!(mechanical.port().law, independent.port().law);
}

#[test]
fn g3_refused_trials_restore_motion_law_and_heat_then_retry_identically() {
    let mut mechanical = mechanics();
    let mut thermal = thermal();
    let route = FrictionHeatRoute::from_link(&thermal, 7, K).unwrap();
    let before = states(&mechanical);
    let before_law = mechanical.port().law.clone();
    let mut exhausted = controls();
    exhausted.maximum_trials = 1;
    assert!(matches!(
        mechanical.step_with_friction_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            exhausted,
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new()
        ),
        Err(ModalFrictionHeatError::TrialBudget)
    ));
    let calls = Cell::new(0usize);
    let refuse = |thermal: &DepositionSimulation<'_>, t| {
        calls.set(calls.get() + 1);
        if calls.get() == 3 {
            Err(DepositionError::MaterialDomain)
        } else {
            source(thermal, t)
        }
    };
    assert!(matches!(
        mechanical.step_with_friction_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
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
    assert!(calls.get() >= 3);
    assert_eq!(states(&mechanical), before);
    assert_eq!(mechanical.port().law, before_law);
    assert_eq!(mechanical.samples_rendered(), 0);
    assert!(mechanical.last_frame().is_none());
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.energy_j(), 0.0);
    let actual = mechanical
        .step_with_friction_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let mut direct = mechanics();
    let mut direct_heat = self::thermal();
    let expected = direct
        .step_with_friction_heat(
            &mut direct_heat,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(actual.mechanical, expected.mechanical);
    assert_eq!(states(&mechanical), states(&direct));
    assert_eq!(thermal.cells(), direct_heat.cells());
    assert_eq!(
        actual.heat.interface_temperature_k,
        expected.heat.interface_temperature_k
    );
}

#[test]
fn g3_cancellation_after_a_staged_trial_and_misaligned_clock_publish_nothing() {
    let mut mechanical = mechanics();
    let mut thermal = thermal();
    let route = FrictionHeatRoute::from_link(&thermal, 7, K).unwrap();
    let before = states(&mechanical);
    let gate = CancelGate::new();
    let calls = Cell::new(0usize);
    let cancel = |thermal: &DepositionSimulation<'_>, t| {
        calls.set(calls.get() + 1);
        if calls.get() == 2 {
            gate.request();
        }
        source(thermal, t)
    };
    assert!(
        mechanical
            .step_with_friction_heat(
                &mut thermal,
                &route,
                &[0.0; 2],
                drive(),
                controls(),
                &[0.0; 2],
                300.0,
                cancel,
                &gate
            )
            .is_err()
    );
    assert!(calls.get() >= 2);
    assert_eq!(states(&mechanical), before);
    assert_eq!(mechanical.port().law, law(300.0));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.energy_j(), 0.0);
    mechanical.step(&[0.0; 2], drive()).unwrap();
    let after = states(&mechanical);
    assert!(matches!(
        mechanical.step_with_friction_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 2],
            300.0,
            source,
            &CancelGate::new()
        ),
        Err(ModalFrictionHeatError::Invalid(_))
    ));
    assert_eq!(states(&mechanical), after);
    assert_eq!(mechanical.samples_rendered(), 1);
    assert_eq!(thermal.time_s(), 0.0);
}
