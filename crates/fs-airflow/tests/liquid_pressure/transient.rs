//! Independent finite-volume and inertance references, using real source cards.
use super::*;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::transient::{
    FiniteLiquidLoop, FiniteLoopControls, FiniteLoopError, FiniteLoopLimits, PipeCell,
};
use fs_conduction::transient::deposition::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent,
    backward_euler::{BoundaryStepCoefficients, StepCoefficients},
};
use fs_solver::NewtonKrylovConfig;
use std::collections::BTreeMap;
use std::convert::Infallible;

const PIPE_AREA: f64 = 1e-4;
const CELL_MASS: f64 = 0.1;
const INERTANCE: f64 = 4e7;
const RESISTANCE: f64 = 4e7;
const CP: f64 = 2000.0;

fn vertices() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
    ]
}
fn pipe_sections() -> Vec<PipeCell> {
    let resistance = MaterialFlowResistance::try_new(
        5e9,
        0.0,
        Area::new(PIPE_AREA),
        Length::new(0.01),
        (0.0, 2300.0),
    )
    .unwrap();
    (0..4)
        .map(|i| PipeCell {
            cell: CellId(i),
            resistance,
        })
        .collect()
}
fn finite_limits() -> FiniteLoopLimits {
    FiniteLoopLimits {
        maximum_cells: 8,
        maximum_relative_density_departure: 0.08,
    }
}
fn finite_controls() -> FiniteLoopControls {
    FiniteLoopControls {
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_newton_iterations: 32,
        maximum_residual_evaluations: 2048,
        energy_scale_j: 1.0,
        flow_scale_m3_s: 1e-5,
        finite_difference_step: 1e-3,
        cell_tolerance_j: 1e-6,
        total_thermal_tolerance_j: 1e-6,
        momentum_tolerance_pa_s: 1e-5,
        mechanical_energy_tolerance_j: 1e-9,
        maximum_neglected_power_w: 1.0,
        maximum_neglected_energy_j: 1.0,
        maximum_drag_step: 2.0,
    }
}
fn heat(powers: [f64; 4]) -> BoundaryStepCoefficients<()> {
    StepCoefficients {
        sources_w: powers
            .into_iter()
            .enumerate()
            .map(|(i, p)| (CellId(i as u64), p))
            .collect(),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    }
    .into()
}
fn temperatures(thermal: &DepositionSimulation<'_>) -> [f64; 4] {
    core::array::from_fn(|i| thermal.temperature_k(CellId(i as u64)).unwrap())
}
fn with_finite<R>(
    card: &MaterialCard,
    temperatures: [f64; 4],
    run: impl FnOnce(
        &Cx<'_>,
        &LiquidTemperatureLaw<'_>,
        &ResolvedHeatCapacityCurve,
        DepositionSimulation<'_>,
        &CancelGate,
    ) -> R,
) -> R {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let cp = PropertyKey::with_quantity(
        "specific_heat_capacity",
        QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
    );
    let rho = PropertyKey::with_quantity(
        "density",
        QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
    );
    let source = ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &cp,
            density: &rho,
            temperature_axis: "T",
            lower: &point(250.0),
            upper: &point(450.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        0.5,
    )
    .unwrap();
    let path = vertices();
    let schedule = temperatures
        .into_iter()
        .enumerate()
        .map(|(i, t)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(1),
                mass_kg: CELL_MASS,
                specific_h_j_kg: CP * (t - 250.0),
                geometry: CellGeometry {
                    start_m: path[i],
                    end_m: path[i + 1],
                    volume_m3: PIPE_AREA,
                    surface_area_m2: 0.04,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    let thermal = DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), source.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links: 16,
            max_events: 8,
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
    .0;
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 743,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &law, &source, thermal, &gate)
    })
}

// Direct cyclic four-volume backward-Euler solution. No production transport,
// caloric inverse, pressure root or nonlinear solver is used by this reference.
fn finite_volume(old: [f64; 4], powers: [f64; 4], average_flow: f64, dt: f64) -> [f64; 4] {
    let c = 1000.0 * average_flow.abs() * dt / CELL_MASS;
    let ratio = c / (1.0 + c);
    let b: [f64; 4] =
        core::array::from_fn(|i| (old[i] + powers[i] * dt / (CELL_MASS * CP)) / (1.0 + c));
    core::array::from_fn(|i| {
        let mut sum = 0.0;
        let mut weight = 1.0;
        for k in 0..4 {
            let upstream = if average_flow >= 0.0 {
                (i + 4 - k) % 4
            } else {
                (i + k) % 4
            };
            sum += weight * b[upstream];
            weight *= ratio;
        }
        sum / (1.0 - ratio.powi(4))
    })
}

#[test]
fn finite_pipe_heats_from_rest_and_closes_buoyancy_momentum_and_enthalpy() {
    with_finite(
        &buoyant_card(),
        [300.0; 4],
        |cx, law, source, mut thermal, _| {
            let mut owner = FiniteLiquidLoop::try_new(
                cx,
                law,
                source,
                &thermal,
                &pipe_sections(),
                SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap(),
                0.0,
                finite_limits(),
            )
            .unwrap();
            assert_eq!(owner.inertance_pa_s2_m3(), INERTANCE);
            let stationary = owner
                .advance(
                    cx,
                    &mut thermal,
                    0.25,
                    0.0,
                    |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                    finite_controls(),
                    &[0.0; 4],
                )
                .unwrap();
            assert_eq!(stationary.hydraulics.average_volume_flow_m3_s, 0.0);
            assert_eq!(owner.volume_flow_m3_s(), 0.0);
            assert_eq!(temperatures(&thermal), [300.0; 4]);

            let powers = [100.0, 200.0, 2000.0, 300.0];
            let dt = 0.5;
            let residual = |q: f64| {
                let t = finite_volume([300.0; 4], powers, q, dt);
                let head = 4.905 * (t[2] - t[0]);
                q - dt * head / (2.0 * INERTANCE + dt * RESISTANCE)
            };
            let (mut lower, mut upper) = (0.0, 1e-5);
            assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lower + upper);
                if residual(mid) > 0.0 {
                    upper = mid;
                } else {
                    lower = mid;
                }
            }
            let q = 0.5 * (lower + upper);
            let expected = finite_volume([300.0; 4], powers, q, dt);
            let report = owner
                .advance(
                    cx,
                    &mut thermal,
                    0.75,
                    0.0,
                    |_, _| Ok::<_, Infallible>(heat(powers)),
                    finite_controls(),
                    &[0.0; 4],
                )
                .unwrap();
            assert!((report.hydraulics.average_volume_flow_m3_s - q).abs() < 1e-13);
            assert!((owner.volume_flow_m3_s() - 2.0 * q).abs() < 2e-13);
            for (actual, expected) in temperatures(&thermal).into_iter().zip(expected) {
                assert!((actual - expected).abs() < 1e-8);
            }
            assert!(report.hydraulics.buoyancy_work_j > 0.0);
            assert!(report.hydraulics.sources[2].density_kg_m3() < 1000.0);
            assert!(
                report
                    .hydraulics
                    .losses
                    .iter()
                    .all(|loss| loss.transport.density_kg_m3() == 1000.0)
            );
            assert!((report.thermal.stored_energy_change_j - 1300.0).abs() < 1e-6);
            assert!((report.thermal.external_heat_j - 1300.0).abs() < 1e-10);
            assert!(report.thermal.residual_j.abs() < 1e-6);
            assert!(report.hydraulics.mechanical_balance_residual_j.abs() < 1e-9);
            assert!(
                report.hydraulics.neglected_energy_bound_j >= report.hydraulics.hydraulic_loss_j
            );
            assert_eq!(owner.accepted_time_s(), thermal.time_s());
        },
    );
}

#[test]
fn actual_midpoint_transport_survives_endpoint_reversal_and_retains_pipe_energy() {
    let old = [320.0, 310.0, 300.0, 290.0];
    with_finite(&card(false, 1.0), old, |cx, law, source, mut thermal, _| {
        let mut owner = FiniteLiquidLoop::try_new(
            cx,
            law,
            source,
            &thermal,
            &pipe_sections(),
            SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap(),
            1e-5,
            finite_limits(),
        )
        .unwrap();
        let dt = 0.5;
        let drive = -1200.0;
        let q = (2.0 * INERTANCE * 1e-5 + dt * drive) / (2.0 * INERTANCE + dt * RESISTANCE);
        assert!(q > 0.0 && 2.0 * q - 1e-5 < 0.0);
        let expected = finite_volume(old, [0.0; 4], q, dt);
        let report = owner
            .advance(
                cx,
                &mut thermal,
                dt,
                drive,
                |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                finite_controls(),
                &[0.0; 4],
            )
            .unwrap();
        assert!((report.hydraulics.average_volume_flow_m3_s - q).abs() < 1e-14);
        assert!((owner.volume_flow_m3_s() - (2.0 * q - 1e-5)).abs() < 1e-14);
        assert!(owner.volume_flow_m3_s() < 0.0);
        for (actual, expected) in temperatures(&thermal).into_iter().zip(expected) {
            assert!((actual - expected).abs() < 1e-8);
        }
        assert!(report.thermal.stored_energy_change_j.abs() < 1e-6);
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert!(report.hydraulics.mechanical_balance_residual_j.abs() < 1e-9);

        let start = owner.volume_flow_m3_s();
        let reverse_average = 2.0 * INERTANCE * start / (2.0 * INERTANCE + dt * RESISTANCE);
        let next_expected = finite_volume(expected, [0.0; 4], reverse_average, dt);
        let reverse = owner
            .advance(
                cx,
                &mut thermal,
                1.0,
                0.0,
                |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                finite_controls(),
                &[0.0; 4],
            )
            .unwrap();
        assert!(reverse.hydraulics.average_volume_flow_m3_s < 0.0);
        for (actual, expected) in temperatures(&thermal).into_iter().zip(next_expected) {
            assert!((actual - expected).abs() < 1e-8);
        }
        assert!(owner.kinetic_energy_j() < 0.5 * INERTANCE * start * start);
        assert!(reverse.thermal.residual_j.abs() < 1e-6);
    });
}

#[test]
fn original_source_viscosity_controls_native_inertial_acceleration() {
    for viscosity_scale in [1.0, 4.0] {
        with_finite(
            &card(false, viscosity_scale),
            [300.0; 4],
            |cx, law, source, mut thermal, _| {
                let mut owner = FiniteLiquidLoop::try_new(
                    cx,
                    law,
                    source,
                    &thermal,
                    &pipe_sections(),
                    SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap(),
                    0.0,
                    finite_limits(),
                )
                .unwrap();
                let q = 0.2 * 400.0 / (2.0 * INERTANCE + 0.2 * RESISTANCE * viscosity_scale);
                let report = owner
                    .advance(
                        cx,
                        &mut thermal,
                        0.2,
                        400.0,
                        |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                        finite_controls(),
                        &[0.0; 4],
                    )
                    .unwrap();
                assert!((report.hydraulics.average_volume_flow_m3_s - q).abs() < 1e-14);
                assert!((owner.volume_flow_m3_s() - 2.0 * q).abs() < 1e-14);
                assert_eq!(temperatures(&thermal), [300.0; 4]);
            },
        );
    }
}

#[test]
fn failures_and_explicit_approximations_leave_both_native_histories_unchanged() {
    with_finite(
        &buoyant_card(),
        [300.0; 4],
        |cx, law, source, mut thermal, gate| {
            let sections = pipe_sections();
            let gravity = || SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap();
            let foreign = source.clone();
            assert!(matches!(
                FiniteLiquidLoop::try_new(
                    cx,
                    law,
                    &foreign,
                    &thermal,
                    &sections,
                    gravity(),
                    1e-5,
                    finite_limits()
                ),
                Err(FiniteLoopError::Invalid(_))
            ));
            let mut owner = FiniteLiquidLoop::try_new(
                cx,
                law,
                source,
                &thermal,
                &sections,
                gravity(),
                1e-5,
                finite_limits(),
            )
            .unwrap();
            let before = thermal.cells().clone();
            let flow = owner.volume_flow_m3_s();

            let mut controls = finite_controls();
            controls.maximum_neglected_power_w = 0.0;
            controls.maximum_neglected_energy_j = 0.0;
            assert!(matches!(
                owner.advance(
                    cx,
                    &mut thermal,
                    0.2,
                    0.0,
                    |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                    controls,
                    &[0.0; 4]
                ),
                Err(FiniteLoopError::NeglectedWork { .. })
            ));
            let failed_source = owner.advance(
                cx,
                &mut thermal,
                0.2,
                0.0,
                |_, _| Err::<BoundaryStepCoefficients<()>, _>("unsupported original wall source"),
                finite_controls(),
                &[0.0; 4],
            );
            assert!(matches!(
                failed_source,
                Err(FiniteLoopError::Heat("unsupported original wall source"))
            ));
            controls = finite_controls();
            controls.maximum_residual_evaluations = 1;
            assert!(matches!(
                owner.advance(
                    cx,
                    &mut thermal,
                    0.2,
                    0.0,
                    |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                    controls,
                    &[0.0; 4]
                ),
                Err(FiniteLoopError::WorkBudget {
                    performed: 1,
                    maximum: 1
                })
            ));
            assert_eq!(thermal.cells(), &before);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(owner.volume_flow_m3_s(), flow);
            assert_eq!(owner.accepted_time_s(), 0.0);

            gate.request();
            assert!(matches!(
                owner.advance(
                    cx,
                    &mut thermal,
                    0.2,
                    0.0,
                    |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                    finite_controls(),
                    &[0.0; 4]
                ),
                Err(FiniteLoopError::Cancelled)
            ));
            assert_eq!(thermal.cells(), &before);
            assert_eq!(owner.volume_flow_m3_s(), flow);
        },
    );
    with_finite(
        &buoyant_card(),
        [301.0; 4],
        |cx, law, source, thermal, _| {
            assert!(matches!(
                FiniteLiquidLoop::try_new(
                    cx,
                    law,
                    source,
                    &thermal,
                    &pipe_sections(),
                    SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap(),
                    0.0,
                    FiniteLoopLimits {
                        maximum_relative_density_departure: 1e-4,
                        ..finite_limits()
                    }
                ),
                Err(FiniteLoopError::DensityDeparture { .. })
            ));
        },
    );
}

#[test]
fn symmetric_reversal_cannot_hide_neglected_work_in_a_zero_average_flow() {
    with_finite(
        &card(false, 1.0),
        [300.0; 4],
        |cx, law, source, mut thermal, _| {
            let initial_flow = 1e-5;
            let dt = 0.5;
            // Independent midpoint identity:
            // Qbar=(2*L*Qold+dt*drive)/(2*L+dt*R)=0, Qnew=-Qold.
            let drive = -2.0 * INERTANCE * initial_flow / dt;
            let mut owner = FiniteLiquidLoop::try_new(
                cx,
                law,
                source,
                &thermal,
                &pipe_sections(),
                SeriesGravity::try_new([0.0, 0.0, -9.81], vertices()).unwrap(),
                initial_flow,
                finite_limits(),
            )
            .unwrap();
            let before = thermal.cells().clone();
            let expected_power_bound = RESISTANCE * initial_flow * initial_flow;
            let expected_energy_bound = dt * expected_power_bound;
            let mut small_cap = finite_controls();
            small_cap.maximum_neglected_energy_j = 1e-4;
            let error = owner
                .advance(
                    cx,
                    &mut thermal,
                    dt,
                    drive,
                    |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                    small_cap,
                    &[0.0; 4],
                )
                .unwrap_err();
            match error {
                FiniteLoopError::NeglectedWork { power_w, energy_j } => {
                    assert!((power_w - expected_power_bound).abs() < 1e-12);
                    assert!((energy_j - expected_energy_bound).abs() < 1e-12);
                    assert!(energy_j > small_cap.maximum_neglected_energy_j);
                }
                other => panic!("symmetric reversal must retain its omitted work: {other:?}"),
            }
            assert_eq!(thermal.cells(), &before);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(owner.accepted_time_s(), 0.0);
            assert_eq!(owner.volume_flow_m3_s(), initial_flow);

            let report = owner
                .advance(
                    cx,
                    &mut thermal,
                    dt,
                    drive,
                    |_, _| Ok::<_, Infallible>(heat([0.0; 4])),
                    finite_controls(),
                    &[0.0; 4],
                )
                .unwrap();
            let hydraulic = &report.hydraulics;
            assert!(hydraulic.average_volume_flow_m3_s.abs() < 1e-14);
            assert!((hydraulic.final_volume_flow_m3_s + initial_flow).abs() < 1e-14);
            assert!(hydraulic.hydraulic_loss_j.abs() < 1e-20);
            assert!(hydraulic.momentum.dissipated.abs() < 1e-20);
            assert!((hydraulic.neglected_work_speed_envelope_m3_s - initial_flow).abs() < 1e-14);
            assert!(
                (hydraulic.neglected_hydraulic_energy_bound_j - expected_energy_bound).abs()
                    < 1e-12
            );
            assert!((hydraulic.neglected_energy_bound_j - expected_energy_bound).abs() < 1e-12);
            assert!(
                hydraulic
                    .envelope_losses
                    .iter()
                    .all(|loss| loss.reynolds > 0.0)
            );
            assert_eq!(temperatures(&thermal), [300.0; 4]);
        },
    );
}
