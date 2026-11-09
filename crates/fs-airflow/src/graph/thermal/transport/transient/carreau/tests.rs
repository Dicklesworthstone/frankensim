//! Independent serial finite-volume enthalpy and original local-viscosity checks.

mod series;

use super::*;
use crate::graph::rheology::CarreauPipeBranch;
use crate::graph::tests::with_cx;
use fs_material::fluid_transport::rheology::pipe::LaminarCarreauPipe;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::{Pressure, VolumetricFlowRate};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, Event, MaterialId,
    StepConfig, ThermalLaw, TimedEvent,
};
use std::collections::BTreeMap;

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.5;
const RHO: f64 = 1000.0;
const CP: f64 = 500.0;
const IDS: [CellId; 5] = [CellId(0), CellId(20), CellId(1), CellId(10), CellId(2)];

struct Caloric {
    datum: f64,
    slope: f64,
}
impl Caloric {
    fn sensible(&self, t: f64) -> f64 {
        let x = t - 300.0;
        CP * x + 0.5 * self.slope * x * x
    }
    fn temperature(&self, h: f64) -> f64 {
        300.0 + 2.0 * h / (CP + (CP * CP + 2.0 * self.slope * h).sqrt())
    }
}
impl ThermalLaw for Caloric {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut id = [37; 32];
        id[..8].copy_from_slice(&self.datum.to_bits().to_le_bytes());
        id[8..16].copy_from_slice(&self.slope.to_bits().to_le_bytes());
        Some(id)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = self.temperature(h - self.datum);
        if t.is_finite() && (280.0..=360.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / (CP - 20.0 * self.slope)
    }
}

fn graph(shear_thinning: bool) -> CarreauPipeGraph {
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            0.2,
            if shear_thinning { 0.04 } else { 0.2 },
            0.02,
            0.7,
            1.5,
            [0.0, 1e7],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 9000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    let pipe = LaminarCarreauPipe::new(law, 300.0, RADIUS, LENGTH, 1e5).unwrap();
    CarreauPipeGraph::new(
        3,
        vec![
            CarreauPipeBranch {
                name: "left".into(),
                from: 0,
                to: 1,
                pipe,
            },
            // Positive physical flow is deliberately opposite this declaration.
            CarreauPipeBranch {
                name: "right".into(),
                from: 2,
                to: 1,
                pipe,
            },
        ],
    )
    .unwrap()
}

fn inventory<'a>(
    law: &'a Caloric,
    temperatures: [f64; 5],
    pipe_volume_scale: f64,
) -> DepositionSimulation<'a> {
    let pipe_volume = core::f64::consts::PI * RADIUS * RADIUS * LENGTH;
    let masses = [
        0.04,
        RHO * pipe_volume * pipe_volume_scale,
        0.08,
        RHO * pipe_volume,
        0.03,
    ];
    let endpoints = [
        ([-0.1, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ([0.0; 3], [0.5, 0.0, 0.0]),
        ([0.5, 0.0, 0.0], [0.5, 0.1, 0.0]),
        ([1.0, 0.0, 0.0], [0.5, 0.0, 0.0]),
        ([1.0, 0.0, 0.0], [1.1, 0.0, 0.0]),
    ];
    let events = (0..5)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: masses[i],
                specific_h_j_kg: law.datum + law.sensible(temperatures[i]),
                geometry: CellGeometry {
                    start_m: endpoints[i].0,
                    end_m: endpoints[i].1,
                    volume_m3: masses[i] / RHO,
                    surface_area_m2: 1.0,
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
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 10,
            max_substeps: 100,
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

fn binding(
    cx: &Cx<'_>,
    state: &DepositionSimulation<'_>,
    shear_thinning: bool,
) -> RegionalCarreauTransport {
    RegionalCarreauTransport::new(
        cx,
        graph(shear_thinning),
        Density::new(RHO),
        &[CellId(0), CellId(1), CellId(2)],
        &[CellId(20), CellId(10)],
        state,
        8,
    )
    .unwrap()
}
fn hydraulics() -> CarreauGraphControls {
    CarreauGraphControls {
        absolute_flow_tolerance: VolumetricFlowRate::new(1e-20),
        relative_flow_tolerance: 2e-15,
        ..CarreauGraphControls::default()
    }
}
fn controls() -> ImplicitTransportControls {
    ImplicitTransportControls {
        solver: NewtonKrylovConfig {
            // Scaled stopping tolerance stays below the independently gated
            // cell/total joules and admits the large-datum coordinate quantum.
            absolute_tolerance: 2e-9,
            relative_tolerance: 0.0,
            ..NewtonKrylovConfig::default()
        },
        maximum_cells: 8,
        maximum_newton_iterations: 30,
        maximum_endpoint_evaluations: 1024,
        energy_scale_j: 1000.0,
        finite_difference_step: 1e-5,
        cell_tolerance_j: 2e-6,
        total_tolerance_j: 5e-6,
    }
}
fn pressures(dp: f64) -> [FixedPressure; 2] {
    [
        FixedPressure {
            node: 0,
            pressure: Pressure::new(dp),
        },
        FixedPressure {
            node: 2,
            pressure: Pressure::new(0.0),
        },
    ]
}
fn exchange(
    law: &Caloric,
    inlet_k: f64,
    flow: &CarreauGraphSolution,
) -> CarreauTransportExchange<()> {
    let q = flow.branches[0].profile.flow().volumetric_flow_m3_s;
    let (inlet, outlet) = if q >= 0.0 {
        (CellId(0), CellId(2))
    } else {
        (CellId(2), CellId(0))
    };
    CarreauTransportExchange {
        mass_boundaries: if q == 0.0 {
            vec![]
        } else {
            vec![
                MassBoundary::Inlet {
                    id: 10,
                    cell: inlet,
                    kg_s: RHO * q.abs(),
                    material_id: MaterialId(1),
                    specific_h_j_kg: law.datum + law.sensible(inlet_k),
                },
                MassBoundary::Outlet {
                    id: 11,
                    cell: outlet,
                    kg_s: RHO * q.abs(),
                },
            ]
        },
        heat: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        }
        .into(),
    }
}
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

// Analytical Newtonian-limit source, evaluated independently of production
// constitutive queries. Both pipes retain separate actual temperatures.
fn resistance(t: f64) -> f64 {
    let mu = 0.2 * (9000.0 / 8.31446261815324 * (1.0 / t - 1.0 / 300.0)).exp();
    8.0 * mu * LENGTH / (core::f64::consts::PI * RADIUS.powi(4))
}
fn bisect(mut low: f64, mut high: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(low) <= 0.0 && f(high) >= 0.0);
    for _ in 0..90 {
        let middle = 0.5 * (low + high);
        if f(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    0.5 * (low + high)
}
fn oracle(
    law: &Caloric,
    old: &DepositionSimulation<'_>,
    dt: f64,
    dp: f64,
    inlet: f64,
) -> (f64, [f64; 5]) {
    let reversed = dp < 0.0;
    let path = if reversed {
        [4, 3, 2, 1, 0]
    } else {
        [0, 1, 2, 3, 4]
    };
    let march = |q: f64| {
        let mut next = [0.0; 5];
        let mut upstream_h = law.sensible(inlet);
        for i in path {
            let cell = &old.cells()[&IDS[i]];
            let old_h = cell.specific_h_j_kg - law.datum;
            let mass = cell.mass_kg;
            let t = if i == 1 || i == 3 {
                bisect(280.0, 360.0, |t| {
                    let h = law.sensible(t);
                    mass * (h - old_h) - dt * (RHO * q * (upstream_h - h) + resistance(t) * q * q)
                })
            } else {
                law.temperature((mass * old_h + dt * RHO * q * upstream_h) / (mass + dt * RHO * q))
            };
            next[i] = t;
            upstream_h = law.sensible(t);
        }
        next
    };
    let q = bisect(0.0, 0.002, |q| {
        let t = march(q);
        q * (resistance(t[1]) + resistance(t[3])) - dp.abs()
    });
    (q.copysign(dp), march(q))
}

fn check_step(
    law: &Caloric,
    old: &DepositionSimulation<'_>,
    next: &DepositionSimulation<'_>,
    report: &ImplicitTransportReport<CarreauTransportEndpoint<()>>,
    dt: f64,
    dp: f64,
    inlet: f64,
) {
    let (q, expected) = oracle(law, old, dt, dp, inlet);
    near(
        report.heat.hydraulics.branches[0]
            .profile
            .flow()
            .volumetric_flow_m3_s,
        q,
        1e-12,
        2e-8,
    );
    let mut stored = 0.0;
    for (i, id) in IDS.into_iter().enumerate() {
        near(next.temperature_k(id).unwrap(), expected[i], 3e-6, 0.0);
        stored += old.cells()[&id].mass_kg
            * (next.cells()[&id].specific_h_j_kg - old.cells()[&id].specific_h_j_kg);
    }
    for (i, slot) in [1, 3].into_iter().enumerate() {
        near(
            report.heat.branch_viscous_heat_w[i],
            resistance(expected[slot]) * q * q,
            1e-8,
            3e-8,
        );
        near(
            report.heat.hydraulics.branches[i]
                .profile
                .pipe()
                .temperature_k(),
            next.temperature_k(IDS[slot]).unwrap(),
            0.0,
            0.0,
        );
    }
    let outlet = if dp > 0.0 { 4 } else { 0 };
    let streams = RHO * q.abs() * (law.sensible(inlet) - law.sensible(expected[outlet]));
    near(report.heat.external_stream_enthalpy_w, streams, 1e-5, 3e-8);
    near(stored, dt * (streams + dp * q), 1e-5, 3e-8);
    near(report.thermal.stored_energy_change_j, stored, 1e-6, 1e-10);
    near(report.thermal.external_heat_j, stored, 5e-6, 1e-10);
    assert!(report.maximum_cell_residual_j < controls().cell_tolerance_j);
}

#[test]
fn g1_regional_carreaus_follow_independent_pipe_inventory_flow_and_heat_root() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 0.0,
            slope: 0.0,
        };
        let old = inventory(&law, [300.0, 310.0, 305.0, 315.0, 298.0], 1.0);
        let owner = binding(cx, &old, false);
        let (next, report) = owner
            .advance(
                cx,
                &old,
                0.5,
                &pressures(20000.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, flow| Ok::<_, ()>(exchange(&law, 340.0, flow)),
            )
            .unwrap();
        check_step(&law, &old, &next, &report, 0.5, 20000.0, 340.0);
        assert!(
            (next.temperature_k(CellId(20)).unwrap() - next.temperature_k(CellId(10)).unwrap())
                .abs()
                > 1.0
        );
        assert_eq!(old.time_s(), 0.0);
    });
}

#[test]
fn g1_reversal_transports_retained_nonconstant_cp_enthalpy_and_separates_loss_heat() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 1e8,
            slope: 4.0,
        };
        let old = inventory(&law, [300.0, 310.0, 305.0, 315.0, 298.0], 1.0);
        let owner = binding(cx, &old, false);
        let (first, report) = owner
            .advance(
                cx,
                &old,
                0.5,
                &pressures(20000.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, flow| Ok::<_, ()>(exchange(&law, 340.0, flow)),
            )
            .unwrap();
        check_step(&law, &old, &first, &report, 0.5, 20000.0, 340.0);
        let (reversed, report) = owner
            .advance(
                cx,
                &first,
                0.75,
                &pressures(-20000.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, flow| Ok::<_, ()>(exchange(&law, 290.0, flow)),
            )
            .unwrap();
        check_step(&law, &first, &reversed, &report, 0.25, -20000.0, 290.0);
        assert!(
            report.heat.hydraulics.branches[0]
                .profile
                .flow()
                .volumetric_flow_m3_s
                < 0.0
        );
        assert!(
            report.heat.hydraulics.branches[1]
                .profile
                .flow()
                .volumetric_flow_m3_s
                > 0.0
        );
        assert!(
            report
                .heat
                .branch_viscous_heat_w
                .iter()
                .all(|&power| power > 0.0)
        );
        assert!(report.heat.external_stream_enthalpy_w < 0.0);
        assert_eq!(first.time_s(), 0.5);
    });
}

#[test]
fn g3_original_shear_thinning_feedback_uses_each_transported_pipe_temperature() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 0.0,
            slope: 2.0,
        };
        let old = inventory(&law, [300.0; 5], 1.0);
        let owner = binding(cx, &old, true);
        let frozen = owner
            .graph()
            .solve(&pressures(2000.0), hydraulics(), cx)
            .unwrap();
        let (next, report) = owner
            .advance(
                cx,
                &old,
                3.0,
                &pressures(2000.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, flow| Ok::<_, ()>(exchange(&law, 345.0, flow)),
            )
            .unwrap();
        assert!(
            report.heat.hydraulics.branches[0]
                .profile
                .flow()
                .volumetric_flow_m3_s
                > 1.05 * frozen.branches[0].profile.flow().volumetric_flow_m3_s,
            "flow must respond to transported original local viscosity"
        );
        let temperatures = [
            next.temperature_k(CellId(20)).unwrap(),
            next.temperature_k(CellId(10)).unwrap(),
        ];
        assert!(temperatures[0] > temperatures[1] + 1.0);
        for (branch, t) in report.heat.hydraulics.branches.iter().zip(temperatures) {
            let profile = &branch.profile;
            assert_eq!(profile.pipe().temperature_k(), t);
            assert_eq!(profile.pipe().law(), owner.graph().branches()[0].pipe.law());
            let source = profile
                .pipe()
                .law()
                .evaluate(t, profile.flow().wall_shear_rate_per_s.abs())
                .unwrap();
            let sample = profile.sample(RADIUS).unwrap();
            near(
                sample.dynamic_viscosity_pa_s,
                source.dynamic_viscosity_pa_s(),
                1e-12,
                1e-10,
            );
        }
    });
}

#[test]
fn g1_stationary_heat_keeps_inventory_and_uses_native_caloric_energy() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 0.0,
            slope: 4.0,
        };
        let old = inventory(&law, [300.0; 5], 1.0);
        let owner = binding(cx, &old, true);
        let (next, report) = owner
            .advance(
                cx,
                &old,
                0.5,
                &pressures(0.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, flow| {
                    let mut response = exchange(&law, 340.0, flow);
                    response.heat.bulk.sources_w.insert(CellId(20), 40.0);
                    Ok::<_, ()>(response)
                },
            )
            .unwrap();
        let cell = &old.cells()[&CellId(20)];
        near(
            next.cells()[&CellId(20)].specific_h_j_kg - cell.specific_h_j_kg,
            20.0 / cell.mass_kg,
            1e-7,
            0.0,
        );
        for id in [CellId(0), CellId(1), CellId(2), CellId(10)] {
            assert_eq!(
                next.cells()[&id].specific_h_j_kg,
                old.cells()[&id].specific_h_j_kg
            );
        }
        assert_eq!(report.heat.branch_viscous_heat_w, vec![0.0; 2]);
        assert_eq!(report.heat.external_stream_enthalpy_w, 0.0);
        near(report.thermal.stored_energy_change_j, 20.0, 1e-7, 0.0);
        // The pre-existing prescribed-flow entry shares the extracted native
        // driver and must retain exactly the same zero-flow heat trajectory.
        let (prescribed, prescribed_report) = owner
            .inventory
            .advance_implicit_with_heat(cx, &old, 0.5, controls(), &[0.0; 5], |_, _| {
                Ok::<_, ()>(
                    StepCoefficients {
                        sources_w: BTreeMap::from([(CellId(20), 40.0)]),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    }
                    .into(),
                )
            })
            .unwrap();
        assert_eq!(prescribed.cells(), next.cells());
        assert_eq!(prescribed.time_s(), next.time_s());
        assert_eq!(prescribed_report.thermal, report.thermal);
    });
}

#[test]
fn g0_native_inventory_stream_source_and_work_refusals_preserve_inputs() {
    with_cx(|cx| {
        let law = Caloric {
            datum: 0.0,
            slope: 0.0,
        };
        let old = inventory(&law, [300.0; 5], 1.0);
        let owner = binding(cx, &old, false);
        let before = old.cells().clone();
        let bad = inventory(&law, [300.0; 5], 1.01);
        assert!(
            RegionalCarreauTransport::new(
                cx,
                graph(false),
                Density::new(RHO),
                &[CellId(0), CellId(1), CellId(2)],
                &[CellId(20), CellId(10)],
                &bad,
                8
            )
            .is_err()
        );
        let other = Caloric {
            datum: 1.0,
            slope: 0.0,
        };
        let replacement = inventory(&other, [300.0; 5], 1.0);
        assert!(
            owner
                .advance(
                    cx,
                    &replacement,
                    0.5,
                    &pressures(20000.0),
                    &[],
                    hydraulics(),
                    controls(),
                    &[0.0; 5],
                    |_, _, flow| Ok::<_, ()>(exchange(&other, 340.0, flow))
                )
                .is_err()
        );
        assert!(
            owner
                .advance(
                    cx,
                    &old,
                    0.5,
                    &pressures(20000.0),
                    &[],
                    hydraulics(),
                    controls(),
                    &[0.0; 5],
                    |_, _, flow| {
                        let mut response = exchange(&law, 340.0, flow);
                        response.mass_boundaries.clear();
                        Ok::<_, ()>(response)
                    }
                )
                .is_err()
        );
        let short = ImplicitTransportControls {
            maximum_endpoint_evaluations: 1,
            ..controls()
        };
        assert!(matches!(
            owner.advance(
                cx,
                &old,
                0.5,
                &pressures(20000.0),
                &[],
                hydraulics(),
                short,
                &[0.0; 5],
                |_, _, flow| Ok::<_, ()>(exchange(&law, 340.0, flow))
            ),
            Err(ImplicitTransportError::WorkBudget { .. })
        ));
        assert!(matches!(
            owner.advance(
                cx,
                &old,
                0.5,
                &pressures(20000.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, _| Err::<CarreauTransportExchange<()>, _>("source refused")
            ),
            Err(ImplicitTransportError::Resolve(
                CarreauTransportError::Exchange("source refused")
            ))
        ));
        assert_eq!(old.cells(), &before);
        assert_eq!(old.time_s(), 0.0);
        owner
            .advance(
                cx,
                &old,
                0.5,
                &pressures(20000.0),
                &[],
                hydraulics(),
                controls(),
                &[0.0; 5],
                |_, _, flow| Ok::<_, ()>(exchange(&law, 340.0, flow)),
            )
            .unwrap();
    });
}

#[test]
fn g4_cancellation_after_real_flow_resolution_publishes_no_thermal_endpoint() {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 711,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let law = Caloric {
            datum: 0.0,
            slope: 0.0,
        };
        let old = inventory(&law, [300.0; 5], 1.0);
        let owner = binding(&cx, &old, false);
        let before = old.cells().clone();
        assert!(
            owner
                .advance(
                    &cx,
                    &old,
                    0.5,
                    &pressures(20000.0),
                    &[],
                    hydraulics(),
                    controls(),
                    &[0.0; 5],
                    |_, _, flow| {
                        gate.request();
                        Ok::<_, ()>(exchange(&law, 340.0, flow))
                    }
                )
                .is_err()
        );
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(old.cells(), &before);
    });
}
