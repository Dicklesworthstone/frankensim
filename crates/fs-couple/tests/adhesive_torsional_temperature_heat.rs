//! Finite temperatures, original thermal claims and loaded changing-area JKR.
//! All material/transport coefficients are manufactured qualification data.
#[path = "../../fs-contact/tests/adhesive_temperature/common.rs"]
mod common;

#[path = "adhesive_torsional_temperature_heat/combined.rs"]
mod combined;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::mixed_mode::torsion::motion::JkrTorsionalMotionBudget;
use fs_contact::normal_patch::adhesive::mixed_mode::{JkrMixedModeEnvelope, SourceJkrMixedMode};
use fs_contact::normal_patch::adhesive::temperature::AdhesiveTemperatures;
use fs_couple::friction_heat::adhesive::torsional_temperature::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::adhesive::jkr_torsional::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::TangentFrame;
use std::cell::Cell;
use std::collections::BTreeMap;

const DT: f64 = 0.001;
const TWIST: f64 = 0.01;
const HALF_LENGTH: f64 = 0.005;
const MASS: [f64; 2] = [0.1, 0.2];
const CAPACITY: [f64; 2] = [0.05, 0.08];
const K: [f64; 2] = [100.0, 40.0];
const BODY_IDS: [u64; 2] = [73, 29];
const LINK: u64 = 19;

fn source(data: &common::Fixture) -> SourceJkrMixedMode {
    SourceJkrMixedMode::new(
        data.jkr(),
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}
fn native_controls() -> RigidJkrTorsionalControls {
    RigidJkrTorsionalControls {
        motion: JkrTorsionalMotionBudget {
            maximum_step_s: 0.1,
            maximum_coordinate_increment: [0.001, 0.1],
            maximum_rate: [100.0, 100.0],
            maximum_effort: [1e4, 1e4],
            maximum_iterations: 64,
            maximum_evaluations: 2048,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance: [1e-16, 1e-16],
            coordinate_relative_tolerance: 1e-12,
            rate_absolute_tolerance: [1e-13, 1e-13],
            rate_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-13,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: 2,
        geometry_tolerance_m: 1e-12,
        displacement_tolerance_m: 1e-12,
        twist_tolerance_rad: 1e-12,
        velocity_tolerance_m_s: 1e-10,
        angular_velocity_tolerance_rad_s: 1e-10,
        linear_momentum_tolerance_ns: 1e-12,
        angular_momentum_tolerance_nms: 1e-12,
        energy_tolerance_j: 1e-11,
    }
}
fn controls() -> JkrTorsionalTemperatureHeatControls {
    JkrTorsionalTemperatureHeatControls {
        thermal: JkrTorsionalThermalControls {
            maximum_cells: 2,
            maximum_interfaces: 1,
            maximum_unknowns: 3,
            maximum_trials: 256,
            maximum_newton_iterations: 16,
            difference_step_j: 1e-7,
            temperature_scale_j_per_k: 0.05,
            heat: FrictionHeatTolerances {
                cell_j: 2e-9,
                thermal_j: 2e-9,
                interface_j: 2e-9,
                total_j: 2e-9,
            },
            combined_energy_tolerance_j: 4e-9,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-11,
                linear_restart: 3,
                forcing_minimum: 1e-6,
                ..Default::default()
            },
            audit_entropy: true,
        },
        geometry_tolerance_m: 1e-12,
        maximum_fractional_area_change: 0.8,
        maximum_absolute_parameter_work_j: 1.0,
    }
}

struct Pair {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
}
impl Pair {
    fn new() -> Self {
        let properties = [
            MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.3)).unwrap(),
            MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.4, 0.4, 0.6)).unwrap(),
        ];
        let bodies = std::array::from_fn(|i| {
            RigidBodyState::new(
                Pose::new(
                    Vec3::new(0.0, 0.0, if i == 0 { -0.5 } else { 0.5 }),
                    UnitQuaternion::IDENTITY,
                )
                .unwrap(),
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .unwrap()
        });
        Self { bodies, properties }
    }
    fn prepared<'a>(&self, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrTorsionalPair<'a> {
        PreparedRigidJkrTorsionalPair::new(
            &self.bodies,
            &self.properties,
            RigidJkrTorsionalContact {
                binding: RigidContactBinding {
                    body_a: 0,
                    body_b: 1,
                    arm_a_body: Vec3::new(0.0, 0.0, 0.5),
                    arm_b_body: Vec3::new(0.0, 0.0, -0.5),
                },
                frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
                reference_approach_m: 0.0,
                source,
            },
            native_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap()
    }
    fn hold_phase(
        &mut self,
        source: &SourceJkrMixedMode,
        previous: Option<RigidJkrTorsionalState>,
    ) -> RigidJkrTorsionalState {
        let native = self.prepared(source);
        let gate = CancelGate::new_clock_free();
        let state = previous.unwrap_or_else(|| native.initialize(TWIST, &gate).unwrap());
        let response = state.motion_state().contact();
        let wrenches = [1.0, -1.0].map(|sign| Wrench {
            force_world: Vec3::new(0.0, 0.0, sign * response.normal_force_n),
            torque_body: Vec3::new(0.0, 0.0, sign * response.torsional_moment_n_m),
        });
        let step = native
            .advance_with_wrenches(state, wrenches, DT, &gate)
            .unwrap();
        self.bodies = step.bodies_after;
        step.state_after
    }
}

fn curves(source: &SourceJkrMixedMode) -> [EquilibriumEnthalpyPhaseCurve; 2] {
    let e = source.normal_source().normal_state().elastic();
    let ids = [
        e.surface_a_state().resolved().card_identity(),
        e.surface_b_state().resolved().card_identity(),
    ];
    std::array::from_fn(|i| {
        EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ids[i],
            SolidLiquidPhase::Solid,
            [280.0, 320.0]
                .map(|t| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: (t - 280.0) * CAPACITY[i] / MASS[i],
                    temperature_k: t,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 1000.0,
                })
                .to_vec(),
        )
        .unwrap()
    })
}
fn finite_thermal<'t>(
    curves: &'t [EquilibriumEnthalpyPhaseCurve; 2],
    area: f64,
) -> DepositionSimulation<'t> {
    let mut events: Vec<_> = (0..2)
        .map(|i| {
            let volume = MASS[i] / 1000.0;
            let face = volume / (2.0 * HALF_LENGTH);
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(i as u64),
                    segment_id: i as u64,
                    body_id: BODY_IDS[i],
                    material_id: MaterialId(i as u64),
                    mass_kg: MASS[i],
                    specific_h_j_kg: 20.0 * CAPACITY[i] / MASS[i],
                    geometry: CellGeometry {
                        start_m: [0.0, 0.0, if i == 0 { -2.0 * HALF_LENGTH } else { 0.0 }],
                        end_m: [0.0, 0.0, if i == 0 { 0.0 } else { 2.0 * HALF_LENGTH }],
                        volume_m3: volume,
                        surface_area_m2: 2.0 * face + 8.0 * face.sqrt() * HALF_LENGTH,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(0),
            b: CellId(1),
            area_m2: area,
            conductance_w_k: series_conductance_w_k(area, [HALF_LENGTH; 2], K).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &curves[0]), (MaterialId(1), &curves[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 4,
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
fn binding<'a>(
    data: &'a common::Fixture,
    curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    route: &'a FrictionHeatRoute,
) -> JkrTorsionalTemperatureHeatBinding<'a> {
    JkrTorsionalTemperatureHeatBinding {
        body_ids: BODY_IDS,
        material_cards: data.cards(),
        caloric: [&curves[0], &curves[1]],
        route,
        parameter_work_source: "prescribed source-coefficient controller",
    }
}
fn resolve(
    q: [f64; 2],
    t: AdhesiveTemperatures,
) -> Result<JkrTorsionalTemperatureHeatSource<AdhesiveTemperatures>, &'static str> {
    Ok(JkrTorsionalTemperatureHeatSource {
        conductivity_w_m_k: K,
        transport: StepCoefficients {
            sources_w: BTreeMap::from([(CellId(0), q[0] / DT), (CellId(1), q[1] / DT)]),
            conductances_w_k: BTreeMap::new(),
            payload: t,
        },
    })
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-8 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

// Independent dimensional source equations at held delta=0:
// E*a^3/R^2 + 8*Gbar*a*phi^2 = 2*pi*W. No production contact solver is used.
fn independent_source(t: AdhesiveTemperatures, twist: f64) -> (f64, [f64; 3]) {
    let e_a = 1e6 * (1.0 - 0.01 * (t.surface_a_k - 300.0));
    let e_b = 2e6 * (1.0 - 0.01 * (t.surface_b_k - 300.0));
    let e = 1.0 / ((1.0 - 0.25_f64.powi(2)) / e_a + (1.0 - 0.25_f64.powi(2)) / e_b);
    let g = 1.0 / (2.5 / e_a + 2.5 / e_b);
    let w = 1.0 - 0.03 * (t.interface_k - 300.0);
    let mut low = 0.0;
    let mut high = (2.0 * PI * w / e).cbrt();
    for _ in 0..100 {
        let a = 0.5 * (low + high);
        if e * a.powi(3) + 8.0 * g * a * twist * twist > 2.0 * PI * w {
            high = a;
        } else {
            low = a;
        }
    }
    let a = 0.5 * (low + high);
    (
        a,
        [
            e * a.powi(5) / 5.0,
            (8.0 / 3.0) * g * a.powi(3) * twist * twist,
            -PI * w * a * a,
        ],
    )
}
fn temperatures_at_area(area: f64, q: [f64; 2]) -> AdhesiveTemperatures {
    let h = DT * series_conductance_w_k(area, [HALF_LENGTH; 2], K).unwrap();
    let d = CAPACITY[0] * CAPACITY[1] + h * (CAPACITY[0] + CAPACITY[1]);
    let a = 300.0 + ((CAPACITY[1] + h) * q[0] + h * q[1]) / d;
    let b = 300.0 + ((CAPACITY[0] + h) * q[1] + h * q[0]) / d;
    AdhesiveTemperatures {
        surface_a_k: a,
        surface_b_k: b,
        interface_k: (K[0] * a + K[1] * b) / (K[0] + K[1]),
    }
}
fn independent_feedback(q: [f64; 2]) -> (AdhesiveTemperatures, f64, [f64; 3]) {
    let original = independent_source(
        AdhesiveTemperatures {
            surface_a_k: 300.0,
            surface_b_k: 300.0,
            interface_k: 300.0,
        },
        TWIST,
    )
    .0;
    let f = |a: f64| a - independent_source(temperatures_at_area(PI * a * a, q), TWIST).0;
    let mut low = 0.2 * original;
    let mut high = 2.0 * original;
    assert!(f(low) < 0.0 && f(high) > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if f(mid) > 0.0 {
            high = mid;
        } else {
            low = mid;
        }
    }
    let a = 0.5 * (low + high);
    let t = temperatures_at_area(PI * a * a, q);
    let storage = independent_source(t, TWIST).1;
    (t, a, storage)
}

#[test]
fn implicit_finite_temperatures_update_original_loaded_area_without_charging_work_as_heat() {
    let data = common::Fixture::new();
    let source = source(&data);
    let curves = curves(&source);
    let mut pair = Pair::new();
    let state = pair.hold_phase(&source, None);
    let native = pair.prepared(&source);
    let before = state.motion_state().contact();
    assert!(state.coordinates()[0].abs() < 1e-13);
    let thermal = finite_thermal(&curves, before.contact_area_m2);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let q = [0.4, -0.08];
    let candidate = propose_jkr_torsional_temperature_heat(
        &native,
        state,
        &thermal,
        binding(&data, &curves, &route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| resolve(q, t),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let report = candidate.report();
    let (expected_t, expected_a, storage) = independent_feedback(q);
    close(
        report.temperatures.surface_a_k,
        expected_t.surface_a_k,
        2e-7,
    );
    close(
        report.temperatures.surface_b_k,
        expected_t.surface_b_k,
        2e-7,
    );
    close(
        report.temperatures.interface_k,
        expected_t.interface_k,
        2e-7,
    );
    close(
        candidate.route().contact_area_m2(),
        PI * expected_a * expected_a,
        2e-12,
    );
    let after = report.temperature.after;
    close(after.contact_radius_m, expected_a, 2e-11);
    close(after.normal_elastic_energy_j, storage[0], 1e-11);
    close(after.torsional_elastic_energy_j, storage[1], 1e-11);
    close(after.surface_energy_j, storage[2], 1e-11);
    close(
        report.temperature.parameter_work_j,
        storage.iter().sum::<f64>() - before.reversible_energy_j,
        2e-11,
    );
    assert!(report.temperature.parameter_work_j.abs() > 1e-7);
    assert!(after.contact_area_m2 < before.contact_area_m2);
    assert!(report.temperature.torsional_elastic_change_j != 0.0);
    assert_eq!(
        report.parameter_work_source,
        "prescribed source-coefficient controller"
    );
    close(
        report.heat.thermal.stored_energy_change_j,
        q.iter().sum(),
        2e-9,
    );
    assert_eq!(report.heat.dissipated_j, 0.0);
    assert_eq!(report.heat.entropy.unwrap().friction_generation_j_k, 0.0);
    close(report.combined_energy_residual_j, 0.0, 4e-9);
    assert_eq!(
        candidate.state().coordinates().map(f64::to_bits),
        state.coordinates().map(f64::to_bits)
    );
    assert_eq!(
        candidate.state().rates().map(f64::to_bits),
        state.rates().map(f64::to_bits)
    );
    assert_eq!(candidate.state().steps(), state.steps());
    assert_eq!(candidate.state().time_s(), state.time_s());
    assert_eq!(candidate.thermal().time_s(), state.time_s());
    assert_eq!(candidate.bodies(), &pair.bodies);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(report.source, report.temperatures);
    let refreshed = candidate.source().normal_source().normal_state();
    assert_eq!(
        refreshed
            .elastic()
            .surface_a_state()
            .resolved()
            .state_coordinate("T"),
        Some(report.temperatures.surface_a_k)
    );
    assert_eq!(
        refreshed
            .elastic()
            .surface_b_state()
            .resolved()
            .state_coordinate("T"),
        Some(report.temperatures.surface_b_k)
    );
    assert_eq!(
        refreshed.interface_state().state_coordinate("T"),
        Some(report.temperatures.interface_k)
    );
    assert_eq!(refreshed.interface_state().state_coordinate("p"), Some(2e5));
    let patch = candidate.thermal().contact_patch(LINK).unwrap();
    close(patch.link.area_m2, candidate.route().contact_area_m2(), 0.0);
    close(
        patch.link.conductance_w_k,
        series_conductance_w_k(patch.link.area_m2, [HALF_LENGTH; 2], K).unwrap(),
        1e-12,
    );
    // The coefficient-temperature successor continues the original native clock
    // and preload; another source phase never initializes a new elastic twist.
    let (successor, bodies, next_state, accepted, _, _) = candidate.into_parts();
    pair.bodies = bodies;
    let later = pair.hold_phase(&successor, Some(next_state));
    assert_eq!(later.steps(), state.steps() + 1);
    assert_eq!(later.time_s(), 2.0 * DT);
    close(later.coordinates()[1], TWIST, 1e-12);
    let next_native = pair.prepared(&successor);
    // Install the actual next mechanical endpoint's coverage before its thermal
    // source phase. Even a holding solve can retain a tiny admitted motion.
    let old_patch = accepted.contact_patch(LINK).unwrap();
    let area = later.motion_state().contact().contact_area_m2;
    let mut next_patch = old_patch;
    next_patch.link.area_m2 = area;
    next_patch.link.conductance_w_k = series_conductance_w_k(area, [HALF_LENGTH; 2], K).unwrap();
    next_patch.covered_area_m2 = [area; 2];
    let accepted = accepted
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: LINK,
                expected: Some(old_patch),
                replacement: Some(next_patch),
            }],
            || Ok(()),
        )
        .unwrap();
    let next_route = FrictionHeatRoute::from_link(&accepted, LINK, K).unwrap();
    let cooling = propose_jkr_torsional_temperature_heat(
        &next_native,
        later,
        &accepted,
        binding(&data, &curves, &next_route),
        controls(),
        &[0.0; 2],
        expected_t.interface_k,
        |_, t, _| resolve([-0.3, -0.02], t),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(cooling.report().temperatures.surface_a_k < expected_t.surface_a_k);
    assert!(cooling.report().temperature.parameter_work_j < 0.0);
    assert_eq!(cooling.report().heat.dissipated_j, 0.0);
    close(
        cooling.report().heat.thermal.stored_energy_change_j,
        -0.32,
        2e-9,
    );
    assert!(cooling.state().coordinates()[1] > 0.0);
}

#[test]
fn source_failure_budget_and_final_trial_cancellation_keep_all_inputs_retryable() {
    let data = common::Fixture::new();
    let source = source(&data);
    let curves = curves(&source);
    let mut pair = Pair::new();
    let state = pair.hold_phase(&source, None);
    let native = pair.prepared(&source);
    let thermal = finite_thermal(&curves, state.motion_state().contact().contact_area_m2);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let original_h: Vec<_> = thermal
        .cells()
        .values()
        .map(|c| c.specific_h_j_kg)
        .collect();
    let calls = Cell::new(0usize);
    let pilot = propose_jkr_torsional_temperature_heat(
        &native,
        state,
        &thermal,
        binding(&data, &curves, &route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| {
            calls.set(calls.get() + 1);
            resolve([0.4, -0.08], t)
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let final_call = calls.get();
    assert_eq!(final_call, pilot.report().thermal_trials);
    assert!(final_call > 2);
    let gate = CancelGate::new_clock_free();
    calls.set(0);
    let cancelled = propose_jkr_torsional_temperature_heat(
        &native,
        state,
        &thermal,
        binding(&data, &curves, &route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| {
            calls.set(calls.get() + 1);
            if calls.get() == final_call {
                gate.request();
            }
            resolve([0.4, -0.08], t)
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), final_call);
    let mut too_little = controls();
    too_little.thermal.maximum_trials = 1;
    assert!(
        propose_jkr_torsional_temperature_heat(
            &native,
            state,
            &thermal,
            binding(&data, &curves, &route),
            too_little,
            &[0.0; 2],
            300.0,
            |_, t, _| resolve([0.4, -0.08], t),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut no_work = controls();
    no_work.maximum_absolute_parameter_work_j = 1e-15;
    assert!(
        propose_jkr_torsional_temperature_heat(
            &native,
            state,
            &thermal,
            binding(&data, &curves, &route),
            no_work,
            &[0.0; 2],
            300.0,
            |_, t, _| resolve([0.4, -0.08], t),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    // A coefficient query cannot extrapolate beyond the original source range,
    // even when the caller proposes a positive finite interface temperature.
    assert!(
        propose_jkr_torsional_temperature_heat(
            &native,
            state,
            &thermal,
            binding(&data, &curves, &route),
            controls(),
            &[0.0; 2],
            321.0,
            |_, t, _| resolve([0.0; 2], t),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        thermal
            .cells()
            .values()
            .map(|c| c.specific_h_j_kg)
            .collect::<Vec<_>>(),
        original_h
    );
    assert_eq!(native.contact().source, &source);
    let retry = propose_jkr_torsional_temperature_heat(
        &native,
        state,
        &thermal,
        binding(&data, &curves, &route),
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| resolve([0.4, -0.08], t),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(retry.state(), pilot.state());
    assert_eq!(retry.source(), pilot.source());
    assert_eq!(retry.bodies(), &pair.bodies);
    assert_eq!(retry.report().temperature, pilot.report().temperature);
}

#[test]
fn original_caloric_cards_area_work_port_and_clock_are_required() {
    let data = common::Fixture::new();
    let source = source(&data);
    let curves = curves(&source);
    let mut pair = Pair::new();
    let initial = pair
        .prepared(&source)
        .initialize(TWIST, &CancelGate::new_clock_free())
        .unwrap();
    let state = pair.hold_phase(&source, Some(initial));
    let native = pair.prepared(&source);
    let thermal = finite_thermal(&curves, state.motion_state().contact().contact_area_m2);
    let route = FrictionHeatRoute::from_link(&thermal, LINK, K).unwrap();
    let invoke = |b, s| {
        propose_jkr_torsional_temperature_heat(
            &native,
            s,
            &thermal,
            b,
            controls(),
            &[0.0; 2],
            300.0,
            |_, t, _| resolve([0.0; 2], t),
            &CancelGate::new_clock_free(),
        )
    };
    let original = binding(&data, &curves, &route);
    let mut wrong = original;
    wrong.caloric.swap(0, 1);
    assert!(invoke(wrong, state).is_err());
    wrong = original;
    wrong.material_cards.swap(0, 1);
    assert!(invoke(wrong, state).is_err());
    wrong = original;
    wrong.parameter_work_source = " ";
    assert!(invoke(wrong, state).is_err());
    wrong = original;
    wrong.body_ids.swap(0, 1);
    assert!(invoke(wrong, state).is_err());
    assert!(invoke(original, initial).is_err());
    // Another owner supplied actual heat at the accepted clock, but the native
    // source was not continued. A zero-supply step must not silently repair it.
    let (stale, _) = thermal
        .apply_heat_j(
            &BTreeMap::from([(CellId(0), CAPACITY[0])]),
            1e-12,
            1e-12,
            || Ok(()),
        )
        .unwrap();
    let attempted = Cell::new(false);
    let refusal = propose_jkr_torsional_temperature_heat(
        &native,
        state,
        &stale,
        original,
        controls(),
        &[0.0; 2],
        300.0,
        |_, t, _| {
            attempted.set(true);
            resolve([0.0; 2], t)
        },
        &CancelGate::new_clock_free(),
    );
    assert!(matches!(
        refusal,
        Err(JkrTorsionalTemperatureHeatError::Invalid(
            "incoming bulk source temperatures differ from their accepted caloric cells"
        ))
    ));
    assert!(!attempted.get());
    close(stale.temperature_k(CellId(0)).unwrap(), 301.0, 1e-12);
    let narrow = finite_thermal(
        &curves,
        0.5 * state.motion_state().contact().contact_area_m2,
    );
    let narrow_route = FrictionHeatRoute::from_link(&narrow, LINK, K).unwrap();
    assert!(
        propose_jkr_torsional_temperature_heat(
            &native,
            state,
            &narrow,
            binding(&data, &curves, &narrow_route),
            controls(),
            &[0.0; 2],
            300.0,
            |_, t, _| resolve([0.0; 2], t),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}
