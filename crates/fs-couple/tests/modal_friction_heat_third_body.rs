//! Native free-body friction coupled to finite layer storage and two body cells.
//! All material coefficients are analytical synthetic fixtures.
use fs_couple::friction_heat::third_body::{ThirdBodyHeatRoute, ThirdBodyHeatTolerances};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::heat::ModalFrictionHeatError;
use fs_couple::render::schedule::force::coupled::contact::friction::heat::third_body::*;
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionDrive, FrictionModalSystem, ModalFrictionConfig, ModalFrictionPort,
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
use std::collections::BTreeMap;

const DT: f64 = 1.0 / 48000.0;
const MASS: [f64; 2] = [0.002, 0.01];
const CAPACITY: [f64; 3] = [1e-6, 3e-6, 1e-6];
const K: [f64; 3] = [0.1, 0.2, 0.3];
const R: [f64; 2] = [0.2, 0.4];
const LAYER: CellId = CellId(2);
const LINKS: [u64; 2] = [7, 8];
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:.16e} != {b:.16e}");
}
fn law(t: f64) -> FrictionLaw {
    FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 0.8 - 0.05 * (t - 300.0),
    }
}
fn mechanics() -> FrictionModalSystem {
    let images = MASS
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
        shapes: vec![MASS[i].sqrt().recip()],
    };
    FrictionModalSystem::new(
        network,
        ModalFrictionPort {
            left: attachment(0),
            right: Some(attachment(1)),
            law: law(300.0),
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
fn thermal() -> DepositionSimulation<'static> {
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
                mass_kg: CAPACITY[i],
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
    for i in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: LINKS[i],
                a: CellId(i as u64),
                b: LAYER,
                area_m2: 1.0,
                conductance_w_k: series_contact_conductance_w_k(1.0, [0.5; 2], [K[i], K[2]], R[i])
                    .unwrap(),
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
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
fn controls() -> ModalThirdBodyHeatSolve {
    ModalThirdBodyHeatSolve {
        maximum_cells: 3,
        maximum_trials: 256,
        maximum_newton_iterations: 20,
        difference_step_j: 1e-8,
        heat: ThirdBodyHeatTolerances {
            cell_j: 1e-12,
            thermal_j: 1e-12,
            layer_j: 1e-12,
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
) -> Result<ModalThirdBodyHeatSource<f64>, DepositionError> {
    Ok(ModalThirdBodyHeatSource {
        law: law(t),
        conductivity_w_m_k: K,
        area_resistance_m2_k_w: R,
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

#[test]
fn g1_finite_layer_feedback_matches_independent_motion_and_caloric_balance() {
    let mut mechanical = mechanics();
    let mut thermal = thermal();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let before = mechanical.total_energy_j().unwrap();
    let g = LINKS.map(|id| thermal.links()[&id].conductance_w_k);
    let weights: [f64; 2] = std::array::from_fn(|i| DT * g[i] / (CAPACITY[i] + DT * g[i]));
    let effective = CAPACITY[2] + CAPACITY[0] * weights[0] + CAPACITY[1] * weights[1];
    let compliance = 0.5 * DT * DT * (1.0 / MASS[0] + 1.0 / MASS[1]);
    // Independent scalar balance: F=.8-.05*theta,
    // effective*theta=F*(dt-compliance*F), N=1.
    let mut lo = 0.0;
    let mut hi = 0.8;
    for _ in 0..70 {
        let f = 0.5 * (lo + hi);
        let residual = (0.8 - f) * effective - 0.05 * f * (DT - compliance * f);
        if residual > 0.0 {
            lo = f;
        } else {
            hi = f;
        }
    }
    let force = 0.5 * (lo + hi);
    let rise = (0.8 - force) / 0.05;
    let result = mechanical
        .step_with_third_body_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    close(result.heat.layer_temperature_k, 300.0 + rise, 2e-8);
    close(result.mechanical.traction_n, -force, 2e-9);
    close(result.heat.layer_stored_change_j, CAPACITY[2] * rise, 1e-12);
    assert!(result.heat.layer_stored_change_j > 0.0);
    assert_eq!(result.source, result.heat.layer_temperature_k);
    assert_eq!(mechanical.port().law, law(result.source));
    for i in 0..2 {
        close(
            thermal.temperature_k(CellId(i as u64)).unwrap(),
            300.0 + weights[i] * rise,
            2e-8,
        );
        close(
            result.heat.body_heat_j[i],
            CAPACITY[i] * weights[i] * rise,
            1e-12,
        );
        let v = mechanical.components()[i].states()[0].velocity_m_sqrt_kg_per_s / MASS[i].sqrt();
        close(
            v,
            if i == 0 {
                1.0 - force * DT / MASS[i]
            } else {
                force * DT / MASS[i]
            },
            2e-10,
        );
    }
    close(
        mechanical.total_energy_j().unwrap() + thermal.energy_j(),
        before,
        1e-12,
    );
    close(
        thermal.energy_j(),
        result.mechanical.friction_dissipation_j,
        1e-12,
    );
    assert!(result.heat.layer_energy_residual_j.abs() < 1e-12);
    assert!(result.combined_energy_residual_j.abs() < 1e-12);
    assert_eq!(mechanical.samples_rendered(), 1);
    assert_eq!(thermal.time_s(), DT);
}

#[test]
fn g3_failed_layer_solve_retries_without_motion_heat_or_branch_history_changes() {
    let mut mechanical = mechanics();
    let mut thermal = thermal();
    let route = ThirdBodyHeatRoute::from_links(&thermal, LAYER, LINKS, K, R).unwrap();
    let before = states(&mechanical);
    let mut exhausted = controls();
    exhausted.maximum_trials = 1;
    let failure = mechanical.step_with_third_body_heat(
        &mut thermal,
        &route,
        &[0.0; 2],
        drive(),
        exhausted,
        &[0.0; 3],
        source,
        &CancelGate::new(),
    );
    assert!(matches!(failure, Err(ModalFrictionHeatError::TrialBudget)));
    assert_eq!(states(&mechanical), before);
    assert_eq!(mechanical.samples_rendered(), 0);
    assert_eq!(mechanical.port().law, law(300.0));
    assert_eq!(thermal.energy_j(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    let stopped = CancelGate::new();
    stopped.request();
    assert!(
        mechanical
            .step_with_third_body_heat(
                &mut thermal,
                &route,
                &[0.0; 2],
                drive(),
                controls(),
                &[0.0; 3],
                source,
                &stopped
            )
            .is_err()
    );
    assert_eq!(states(&mechanical), before);
    let retried = mechanical
        .step_with_third_body_heat(
            &mut thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    let mut fresh = mechanics();
    let mut fresh_thermal = self::thermal();
    let direct = fresh
        .step_with_third_body_heat(
            &mut fresh_thermal,
            &route,
            &[0.0; 2],
            drive(),
            controls(),
            &[0.0; 3],
            source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(states(&mechanical), states(&fresh));
    assert_eq!(thermal.cells(), fresh_thermal.cells());
    assert_eq!(retried.mechanical, direct.mechanical);
    assert_eq!(
        retried.heat.layer_temperature_k.to_bits(),
        direct.heat.layer_temperature_k.to_bits()
    );
}
