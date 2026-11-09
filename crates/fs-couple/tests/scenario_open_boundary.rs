//! G0/G3: one control-time transaction reaches both native boundary operators.

#[path = "scenario_open_boundary/thermal.rs"]
mod thermal;

use fs_couple::scenario_open_boundary::lattice3;
use fs_couple::scenario_open_boundary::{
    BoundarySample3, ScenarioBoundary3, ScenarioBoundaryError, rebind_open_faces3,
};
use fs_couple::scenario_pressure::{
    PressureBoundary3, PressureReference, ScenarioPressureError, lattice3 as pressure_lattice3,
};
use fs_exec::{CancelGate, Cx, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, CollisionModel3, Face3, FaceBoundary3, InletFlux3,
};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
use fs_lbm::lattice3::{D3q19, EquilibriumLattice3, Grid3, GridBudget3, GridError3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, Payload, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, FrameTree, LinearChirp, Physics, TimeSignal, WORLD,
};
use std::sync::atomic::{AtomicUsize, Ordering};

const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn velocity(value: [f64; 3]) -> BoundaryCondition {
    velocity_source(SampleSource::Fixed(
        value.map(|v| QtyAny::new(v, VELOCITY)).to_vec(),
    ))
}

fn velocity_source(source: SampleSource<Vec<QtyAny>>) -> BoundaryCondition {
    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(VELOCITY),
        PayloadId::new("basis/world-cartesian").unwrap(),
        WORLD,
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .unwrap();
    BoundaryCondition {
        region: "inlet".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::Dirichlet,
        compatibility: None,
        frame: WORLD.0,
        value: Some(BcValue::Typed(Payload::Vector(
            VectorPayload::new(meta, source).unwrap(),
        ))),
    }
}

fn pressure(value: f64) -> BoundaryCondition {
    BoundaryCondition {
        region: "outlet".into(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        compatibility: None,
        frame: WORLD.0,
        value: Some(BcValue::Uniform(QtyAny::new(value, PRESSURE))),
    }
}

fn controls<'a>(
    inlet: &'a BoundaryCondition,
    outlet: &'a BoundaryCondition,
) -> [Option<ScenarioBoundary3<'a>>; 6] {
    [
        None,
        None,
        None,
        None,
        Some(ScenarioBoundary3::Velocity(inlet)),
        Some(ScenarioBoundary3::Pressure(PressureBoundary3 {
            boundary: outlet,
            reference: PressureReference::Gauge,
        })),
    ]
}

fn boundaries() -> BoundarySpec3 {
    BoundarySpec3::new([
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::Periodic,
        FaceBoundary3::Velocity { velocity: [0.0; 3] },
        FaceBoundary3::Pressure { density: 1.0 },
    ])
}

fn fluid() -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 4200.0, 0.6, 0.001).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 5e-14, "{actual} != {expected}");
}

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            fs_exec::StreamKey {
                seed: 0xB0A0,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn generic_grid(cx: &Cx<'_>, tau: f64, maximum_mach: f64) -> OpenBoundaryGrid3<D3q27, 27> {
    grid_with_collision(cx,
        CollisionConfig27 {
            model: CollisionModel27::Bgk {
                tau,
                equilibrium: Equilibrium27::Factorized,
            },
            maximum_mach,
        },
    )
}

fn grid_with_collision<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
) -> OpenBoundaryGrid3<L, Q> {
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        [4, 4, 5],
        collision,
        [2e-5, -1e-5, 3e-5],
        [true, true, false],
        GridBudget3 { maximum_cells: 80, maximum_mach: 0.2 },
        |_| false,
        |_| (1.0, [0.0; 3]),
    )
    .unwrap();
    OpenBoundaryGrid3::try_new(grid, cx, boundaries()).unwrap()
}

#[test]
fn mixed_controls_reach_native_velocity_and_pressure_on_the_same_step() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let mut flow = BoundaryGrid3::new(
        8, 8, 8, scaling.flow().tau, [1e-5, -2e-5, 3e-5], boundaries(),
    );
    let inlet = velocity([0.2, -0.1, 0.4]);
    let mut outlet = pressure(0.0);
    outlet.value = Some(BcValue::Signal(TimeSignal::Ramp {
        t_start: 0.0,
        t_end: 1.0,
        from: QtyAny::new(100.0, PRESSURE),
        to: QtyAny::new(300.0, PRESSURE),
    }));
    let before = (flow.populations(3, 3, 0), flow.populations(3, 3, 7));
    let force = flow.force_density(3, 3, 3);
    let samples = rebind_open_faces3(
        &mut flow, controls(&inlet, &outlet), &FrameTree::new(), 0.5,
        &scaling, &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(samples[..4].iter().all(Option::is_none));
    let Some(BoundarySample3::Velocity(v)) = samples[4] else { panic!("missing velocity"); };
    let Some(BoundarySample3::Pressure(p)) = samples[5] else { panic!("missing pressure"); };
    assert_eq!(v.time_s, p.time_s);
    assert_eq!(p.time_s, 0.5);
    assert_eq!(p.perturbation_pa, 200.0);
    let expected_density = 1.0 + 200.0 / ((1.0 / 3.0) * 1000.0 * 20.0 * 20.0);
    near(p.density_lattice, expected_density);
    assert_eq!((flow.populations(3, 3, 0), flow.populations(3, 3, 7)), before);
    assert_eq!(flow.force_density(3, 3, 3), force);
    flow.step();
    for y in 0..8 {
        for x in 0..8 {
            let actual = flow.velocity(x, y, 0);
            for (value, expected) in actual.into_iter().zip([0.01, -0.005, 0.02]) {
                near(value, expected);
            }
            near(flow.density(x, y, 7), expected_density);
        }
    }
}

#[test]
fn later_pressure_and_native_kind_refusals_leave_the_inlet_unchanged() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity([0.0, 0.0, 0.4]);
    let outlet = pressure(f64::NAN);
    let gate = CancelGate::new_clock_free();
    let frames = FrameTree::new();
    let mut flow = BoundaryGrid3::new(8, 8, 8, scaling.flow().tau, [0.0; 3], boundaries());
    let before = flow.boundaries();
    let populations = flow.populations(3, 3, 0);
    assert!(matches!(
        rebind_open_faces3(&mut flow, controls(&inlet, &outlet), &frames, 0.5, &scaling, &gate),
        Err(ScenarioBoundaryError::Pressure(_))
    ));
    assert_eq!(flow.boundaries(), before);
    assert_eq!(flow.populations(3, 3, 0), populations);

    let mut paired = Face3::ALL.map(|face| boundaries().face(face));
    paired[Face3::ZMax as usize] = FaceBoundary3::Velocity { velocity: [0.0; 3] };
    let mut wrong_kind = BoundaryGrid3::new(
        8, 8, 8, scaling.flow().tau, [0.0; 3], BoundarySpec3::new(paired),
    );
    let before = wrong_kind.boundaries();
    assert!(matches!(
        rebind_open_faces3(
            &mut wrong_kind, controls(&inlet, &pressure(200.0)), &frames, 0.5, &scaling, &gate,
        ),
        Err(ScenarioBoundaryError::Boundary(_))
    ));
    assert_eq!(wrong_kind.boundaries(), before);

    gate.request();
    assert!(rebind_open_faces3(
        &mut flow, controls(&inlet, &pressure(200.0)), &frames, 0.5, &scaling, &gate,
    ).is_err());
    assert_eq!(flow.boundaries(), boundaries());
    assert_eq!(flow.populations(3, 3, 0), populations);
}

#[test]
fn mixed_d3q27_controls_keep_the_native_step_and_collision_envelope() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity([0.2, -0.1, 0.4]);
    let outlet = pressure(200.0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = generic_grid(cx, scaling.flow().tau, 0.2);
        let before = flow.populations().to_vec();
        let force = flow.force_density();
        let samples = lattice3::rebind_open_faces3(
            &mut flow, controls(&inlet, &outlet), &frames, 0.5, &scaling, &gate,
        ).unwrap();
        assert_eq!(flow.populations(), before);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(flow.force_density(), force);
        let Some(BoundarySample3::Pressure(p)) = samples[5] else { panic!("missing pressure"); };
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        assert_eq!(flow.completed_steps(), 1);
        for y in 0..4 {
            for x in 0..4 {
                let inlet_state = flow.macroscopic(y * 4 + x).unwrap();
                for (value, expected) in inlet_state.1.into_iter().zip([0.01, -0.005, 0.02]) {
                    near(value, expected);
                }
                near(flow.macroscopic((4 * 4 + y) * 4 + x).unwrap().0, p.density_lattice);
            }
        }

        let mut restricted = generic_grid(cx, scaling.flow().tau, 0.04);
        let before = restricted.boundaries();
        let populations = restricted.populations().to_vec();
        let too_fast = velocity([0.0, 0.0, 0.8]);
        assert!(matches!(
            lattice3::rebind_open_faces3(
                &mut restricted, controls(&too_fast, &outlet), &frames, 0.5, &scaling, &gate,
            ),
            Err(lattice3::Error::Boundary(GridError3::Mach { maximum, .. })) if maximum == 0.04
        ));
        assert_eq!(restricted.boundaries(), before);
        assert_eq!(restricted.populations(), populations);
        assert_eq!(restricted.completed_steps(), 0);
    });
}

#[test]
fn mixed_profile_refusal_keeps_both_controls_at_the_same_time() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity([0.0, 0.0, 0.4]);
    let dormant = pressure(0.0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = generic_grid(cx, scaling.flow().tau, 0.2);
        flow.set_pressure_profile(Face3::ZMax, &gate, |point| {
            Ok(1.002 + 0.0001 * point[0] as f64)
        }).unwrap();
        let profile = flow.pressure_profile(Face3::ZMax).unwrap().to_vec();
        let before = flow.boundaries();
        let populations = flow.populations().to_vec();
        let force = flow.force_density();
        assert!(matches!(
            lattice3::rebind_open_faces3(
                &mut flow, controls(&inlet, &dormant), &frames, 0.5, &scaling, &gate,
            ),
            Err(lattice3::Error::Scenario(ScenarioBoundaryError::Pressure(_)))
        ));
        assert_eq!(flow.boundaries(), before);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(flow.force_density(), force);
        assert_eq!(flow.pressure_profile(Face3::ZMax).unwrap(), profile);

        let samples = lattice3::rebind_open_faces3(
            &mut flow,
            [None, None, None, None, Some(ScenarioBoundary3::Velocity(&inlet)), None],
            &frames, 0.5, &scaling, &gate,
        ).unwrap();
        assert!(samples[5].is_none());
        assert_eq!(flow.pressure_profile(Face3::ZMax).unwrap(), profile);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        for y in 0..4 {
            for x in 0..4 {
                near(flow.macroscopic(y * 4 + x).unwrap().1[2], 0.02);
                near(
                    flow.macroscopic((4 * 4 + y) * 4 + x).unwrap().0,
                    1.002 + 0.0001 * x as f64,
                );
            }
        }

        flow.clear_pressure_profile(Face3::ZMax);
        let samples = lattice3::rebind_open_faces3(
            &mut flow, controls(&inlet, &dormant), &frames, 1.0, &scaling, &gate,
        ).unwrap();
        assert!(matches!(samples[5], Some(BoundarySample3::Pressure(_))));
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        assert_eq!(flow.completed_steps(), 2);
        for y in 0..4 {
            for x in 0..4 {
                near(flow.macroscopic((4 * 4 + y) * 4 + x).unwrap().0, 1.0);
            }
        }
    });
}

fn exercise_scenario_steps<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity_source(
        SampleSource::table(
            vec![
                QtyAny::new(0.0, Dims([0, 0, 1, 0, 0, 0])),
                QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
            ],
            [[0.1, -0.05, 0.2], [0.3, -0.15, 0.6]]
                .into_iter()
                .map(|v| v.map(|value| QtyAny::new(value, VELOCITY)).to_vec())
                .collect(),
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        )
        .unwrap(),
    );
    let mut outlet = pressure(0.0);
    outlet.value = Some(BcValue::Signal(TimeSignal::Ramp {
        t_start: 0.0,
        t_end: 1.0,
        from: QtyAny::new(100.0, PRESSURE),
        to: QtyAny::new(300.0, PRESSURE),
    }));
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut serial = grid_with_collision::<L, Q>(cx, collision);
        let mut parallel = grid_with_collision::<L, Q>(cx, collision);
        let mut reference = grid_with_collision::<L, Q>(cx, collision);
        let pool = TilePool::for_host(2, 0xB0A0);
        let origin = [0.5, -0.25, 1.0];
        for (step, time_s) in [0.25, 0.75].into_iter().enumerate() {
            let previous_targets = serial.boundaries();
            let samples = lattice3::rebind_open_faces3(
                &mut reference,
                controls(&inlet, &outlet),
                &frames,
                time_s,
                &scaling,
                &gate,
            )
            .unwrap();
            let transport = reference.step_serial(cx, None, origin).unwrap();
            let accepted = lattice3::step_serial(
                &mut serial,
                controls(&inlet, &outlet),
                &frames,
                time_s,
                &scaling,
                cx,
                None,
                origin,
            )
            .unwrap();
            let pooled = lattice3::step_pooled(
                &mut parallel,
                controls(&inlet, &outlet),
                &frames,
                time_s,
                &scaling,
                &pool,
                &gate,
                None,
                origin,
            )
            .unwrap();
            assert_eq!(accepted.samples, samples);
            assert_eq!(accepted.transport, transport);
            assert_eq!(pooled, accepted);
            assert_eq!(serial.populations(), reference.populations());
            assert_eq!(parallel.populations(), serial.populations());
            assert_ne!(serial.boundaries(), previous_targets);
            assert_eq!(serial.boundaries(), reference.boundaries());
            assert_eq!(parallel.boundaries(), serial.boundaries());
            assert_eq!(serial.completed_steps(), step + 1);
            assert_eq!(parallel.completed_steps(), step + 1);

            let expected_velocity = [
                (0.1 + 0.2 * time_s) / 20.0,
                (-0.05 - 0.1 * time_s) / 20.0,
                (0.2 + 0.4 * time_s) / 20.0,
            ];
            let expected_density =
                1.0 + (100.0 + 200.0 * time_s) / ((1.0 / 3.0) * 1000.0 * 20.0 * 20.0);
            for cell in 0..16 {
                for (actual, expected) in serial
                    .macroscopic(cell)
                    .unwrap()
                    .1
                    .into_iter()
                    .zip(expected_velocity)
                {
                    near(actual, expected);
                }
                near(serial.macroscopic(64 + cell).unwrap().0, expected_density);
            }
        }
    });
}

#[test]
fn changing_scenario_histories_advance_both_lattices_and_execution_paths() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let tau = scaling.flow().tau;
    exercise_scenario_steps::<D3q19, 19>(CollisionModel3::Bgk { tau });
    exercise_scenario_steps::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn unresolved_pressure_cycles_refuse_before_publication_and_finer_scaling_resolves_them() {
    let material = fluid();
    let dt = 0.001;
    let scaling = plan_scaling_from_material(&material, 0.02, dt, 1.0).unwrap();
    let inlet = velocity([0.0, 0.0, 0.4]);
    let mut outlet = pressure(0.0);
    outlet.value = Some(BcValue::Signal(TimeSignal::Chirp(
        LinearChirp::new(
            0.0,
            1.0,
            1.0 / dt,
            1.0 / dt,
            0.0,
            QtyAny::new(100.0, PRESSURE),
            QtyAny::new(100.0, PRESSURE),
        )
        .unwrap(),
    )));
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    let pool = TilePool::for_host(2, 0xB0A0);
    with_cx(&gate, |cx| {
        let mut flow = generic_grid(cx, scaling.flow().tau, 0.2);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        let targets = flow.boundaries();
        let populations = flow.populations().to_vec();
        let force = flow.force_density();
        // Point sampling at this clock aliases every complete cycle to its peak.
        for time in [0.0, dt, 2.0 * dt] {
            let sample = fs_couple::scenario_pressure::sample_pressure3(
                PressureBoundary3 {
                    boundary: &outlet,
                    reference: PressureReference::Gauge,
                },
                &frames,
                time,
                &scaling,
                &mut || std::ops::ControlFlow::Continue(()),
            )
            .unwrap();
            near(sample.supplied_pa, 200.0);
        }
        for pooled in [false, true] {
            let result = if pooled {
                lattice3::step_pooled(
                    &mut flow,
                    controls(&inlet, &outlet),
                    &frames,
                    dt,
                    &scaling,
                    &pool,
                    &gate,
                    None,
                    [0.0; 3],
                )
            } else {
                lattice3::step_serial(
                    &mut flow,
                    controls(&inlet, &outlet),
                    &frames,
                    dt,
                    &scaling,
                    cx,
                    None,
                    [0.0; 3],
                )
            };
            let error = result.unwrap_err();
            assert!(error.to_string().contains("smaller timestep"));
            assert!(matches!(error,
                lattice3::Error::Scenario(ScenarioBoundaryError::FrequencyResolution {
                    face: Face3::ZMax, time_s, timestep_s,
                }) if time_s == dt && timestep_s == dt
            ));
            assert_eq!(flow.boundaries(), targets);
            assert_eq!(flow.populations(), populations);
            assert_eq!(flow.force_density(), force);
            assert_eq!(flow.completed_steps(), 1);
        }
        // Explicit point rebinding still has no implied integration interval.
        lattice3::rebind_open_faces3(
            &mut flow,
            controls(&inlet, &outlet),
            &frames,
            dt,
            &scaling,
            &gate,
        )
        .unwrap();
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 1);

        let fine = plan_scaling_from_material(&material, 0.02, dt / 32.0, 1.0).unwrap();
        let mut serial = generic_grid(cx, fine.flow().tau, 0.2);
        let mut parallel = generic_grid(cx, fine.flow().tau, 0.2);
        let sample_time = fine.timestep_s() * 0.5;
        let accepted = lattice3::step_serial(
            &mut serial,
            controls(&inlet, &outlet),
            &frames,
            sample_time,
            &fine,
            cx,
            None,
            [0.0; 3],
        )
        .unwrap();
        let pooled = lattice3::step_pooled(
            &mut parallel,
            controls(&inlet, &outlet),
            &frames,
            sample_time,
            &fine,
            &pool,
            &gate,
            None,
            [0.0; 3],
        )
        .unwrap();
        assert_eq!(pooled, accepted);
        assert_eq!(parallel.populations(), serial.populations());
        assert_eq!(serial.completed_steps(), 1);
        assert_eq!(parallel.completed_steps(), 1);
        let Some(BoundarySample3::Pressure(sample)) = accepted.samples[5] else {
            panic!("missing pressure");
        };
        assert!(sample.supplied_pa < 200.0);
    });
}

#[test]
fn refused_native_step_keeps_scenario_controls_and_the_accepted_fluid_prefix() {
    struct CancelAfterStream {
        pool: TilePool,
        calls: AtomicUsize,
    }
    impl KernelRunner for CancelAfterStream {
        fn workers(&self) -> usize {
            self.pool.workers()
        }
        fn run_with_gate<K: TileKernel>(
            &self,
            kernel: &K,
            gate: &CancelGate,
        ) -> (Result<K::Out, RunError>, RunReport) {
            let result = self.pool.run_with_gate(kernel, gate);
            if self.calls.fetch_add(1, Ordering::AcqRel) == 1 {
                gate.request();
            }
            result
        }
    }
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity([0.0, 0.0, 0.4]);
    let outlet = pressure(200.0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = generic_grid(cx, scaling.flow().tau, 0.2);
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        let targets = flow.boundaries();
        let populations = flow.populations().to_vec();
        let mut laws = vec![
            CollisionConfig27 {
                model: CollisionModel27::Bgk {
                    tau: scaling.flow().tau,
                    equilibrium: Equilibrium27::Factorized,
                },
                maximum_mach: 0.2,
            };
            80
        ];
        // The accepted near-rest state admits this law. The changed inlet
        // fails only after collision and the earlier reconstructed face cells.
        laws[15].maximum_mach = 0.01;
        assert!(matches!(
            lattice3::step_serial(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                Some(&laws),
                [0.0; 3],
            ),
            Err(lattice3::Error::Boundary(GridError3::Mach { cell: 15, .. }))
        ));
        assert_eq!(flow.boundaries(), targets);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 1);

        let runner = CancelAfterStream {
            pool: TilePool::for_host(2, 0xB0A0),
            calls: AtomicUsize::new(0),
        };
        assert!(matches!(
            lattice3::step_pooled(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                &runner,
                &CancelGate::new_clock_free(),
                None,
                [0.0; 3],
            ),
            Err(lattice3::Error::Boundary(GridError3::Cancelled))
        ));
        assert_eq!(runner.calls.load(Ordering::Acquire), 2);
        assert_eq!(flow.boundaries(), targets);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 1);

        let mut reference = generic_grid(cx, scaling.flow().tau, 0.2);
        reference.step_serial(cx, None, [0.0; 3]).unwrap();
        let expected_samples = lattice3::rebind_open_faces3(
            &mut reference,
            controls(&inlet, &outlet),
            &frames,
            0.5,
            &scaling,
            &gate,
        )
        .unwrap();
        let expected_transport = reference.step_serial(cx, None, [0.0; 3]).unwrap();
        let retry = lattice3::step_serial(
            &mut flow,
            controls(&inlet, &outlet),
            &frames,
            0.5,
            &scaling,
            cx,
            None,
            [0.0; 3],
        )
        .unwrap();
        assert_eq!(retry.samples, expected_samples);
        assert_eq!(retry.transport, expected_transport);
        assert_eq!(flow.populations(), reference.populations());
        assert_eq!(flow.boundaries(), reference.boundaries());
        assert_eq!(flow.completed_steps(), 2);
    });
}

#[test]
fn selected_pressure_history_cannot_be_hidden_by_an_installed_profile() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity([0.0, 0.0, 0.4]);
    // This maps to the dormant uniform density already held by the grid.
    let outlet = pressure(0.0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = generic_grid(cx, scaling.flow().tau, 0.2);
        flow.set_pressure_profile(Face3::ZMax, &gate, |p| Ok(1.01 + 0.001 * p[0] as f64))
            .unwrap();
        let profile = flow.pressure_profile(Face3::ZMax).unwrap().to_vec();
        let targets = flow.boundaries();
        let populations = flow.populations().to_vec();
        let input = PressureBoundary3 {
            boundary: &outlet,
            reference: PressureReference::Gauge,
        };
        assert!(matches!(
            pressure_lattice3::rebind_pressure_face3(
                &mut flow, Face3::ZMax, input, &frames, 0.5, &scaling, &gate,
            ),
            Err(pressure_lattice3::Error::Scenario(ScenarioPressureError::Boundary(_)))
        ));
        let mut pressures = [None; 6];
        pressures[Face3::ZMax as usize] = Some(input);
        assert!(matches!(
            pressure_lattice3::rebind_pressure_faces3(
                &mut flow, pressures, &frames, 0.5, &scaling, &gate,
            ),
            Err(pressure_lattice3::Error::Scenario(ScenarioPressureError::Boundary(_)))
        ));
        assert!(matches!(
            lattice3::rebind_open_faces3(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                &gate,
            ),
            Err(lattice3::Error::Scenario(ScenarioBoundaryError::Pressure(
                _
            )))
        ));
        assert!(matches!(
            lattice3::step_serial(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                None,
                [0.0; 3],
            ),
            Err(lattice3::Error::Scenario(ScenarioBoundaryError::Pressure(
                _
            )))
        ));
        assert!(matches!(
            lattice3::step_pooled(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                &TilePool::for_host(2, 0xB0A0),
                &gate,
                None,
                [0.0; 3],
            ),
            Err(lattice3::Error::Scenario(ScenarioBoundaryError::Pressure(
                _
            )))
        ));
        assert_eq!(flow.boundaries(), targets);
        assert_eq!(flow.populations(), populations);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(flow.pressure_profile(Face3::ZMax), Some(profile.as_slice()));

        let mut inlet_only = controls(&inlet, &outlet);
        inlet_only[Face3::ZMax as usize] = None;
        let accepted = lattice3::step_serial(
            &mut flow, inlet_only, &frames, 0.5, &scaling, cx, None, [0.0; 3],
        )
        .unwrap();
        assert_eq!(accepted.samples[5], None);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(flow.pressure_profile(Face3::ZMax), Some(profile.as_slice()));
        for (cell, density) in profile.into_iter().enumerate() {
            near(flow.macroscopic(64 + cell).unwrap().0, density);
        }
    });
}

fn force_epoch(epoch: usize) -> Vec<[f64; 3]> {
    let sign = if epoch % 2 == 0 { 1.0 } else { -1.0 };
    (0..80)
        .map(|cell| {
            [
                sign * 2e-5 * (1 + cell % 3) as f64,
                -sign * 1e-5 * (1 + cell % 5) as f64,
                sign * 3e-5 * (1 + cell % 7) as f64,
            ]
        })
        .collect()
}

fn forced_grid<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
) -> ForcedOpenBoundaryGrid3<L, Q> {
    ForcedOpenBoundaryGrid3::try_new_preserving_velocity(
        grid_with_collision(cx, collision),
        cx,
        force_epoch(0),
    )
    .unwrap()
}

// Independent physical moments, including each cell's accepted half force.
fn forced_totals<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &ForcedOpenBoundaryGrid3<L, Q>,
) -> (f64, [f64; 3]) {
    let mass = flow.populations().iter().flat_map(|f| f.iter()).sum();
    let momentum = std::array::from_fn(|axis| {
        flow.populations()
            .iter()
            .zip(flow.forces())
            .map(|(f, force)| {
                f.iter()
                    .zip(L::E)
                    .map(|(value, direction)| value * f64::from(direction[axis]))
                    .sum::<f64>()
                    + 0.5 * force[axis]
            })
            .sum()
    });
    (mass, momentum)
}

fn forced_targets(
    inlet: &BoundaryCondition,
    outlet: &BoundaryCondition,
    frames: &FrameTree,
    time_s: f64,
    scaling: &fs_lbm::material::MaterialScaling<'_, FluidTransportState>,
) -> BoundarySpec3 {
    // Keep the pre-existing source/scaling arithmetic identical so population
    // equality isolates native force/reservoir composition, not decimal ulps.
    let velocity =
        fs_couple::scenario_velocity::sample_velocity3(inlet, frames, time_s, scaling, &mut || {
            std::ops::ControlFlow::Continue(())
        })
        .unwrap();
    let pressure = fs_couple::scenario_pressure::sample_pressure3(
        PressureBoundary3 {
            boundary: outlet,
            reference: PressureReference::Gauge,
        },
        frames,
        time_s,
        scaling,
        &mut || std::ops::ControlFlow::Continue(()),
    )
    .unwrap();
    let mut targets = Face3::ALL.map(|face| boundaries().face(face));
    targets[4] = FaceBoundary3::Velocity {
        velocity: velocity.lattice,
    };
    targets[5] = FaceBoundary3::Pressure {
        density: pressure.density_lattice,
    };
    BoundarySpec3::new(targets)
}

fn exercise_forced_scenario<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity_source(
        SampleSource::table(
            [0.0, 1.0]
                .map(|t| QtyAny::new(t, Dims([0, 0, 1, 0, 0, 0])))
                .to_vec(),
            [[0.1, -0.05, 0.2], [0.3, -0.15, 0.6]]
                .into_iter()
                .map(|v| v.map(|value| QtyAny::new(value, VELOCITY)).to_vec())
                .collect(),
            TableInterpolation::Linear,
            OutsideDomainPolicy::Refuse,
        )
        .unwrap(),
    );
    let mut outlet = pressure(0.0);
    outlet.value = Some(BcValue::Signal(TimeSignal::Ramp {
        t_start: 0.0,
        t_end: 1.0,
        from: QtyAny::new(100.0, PRESSURE),
        to: QtyAny::new(300.0, PRESSURE),
    }));
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut actual = forced_grid::<L, Q>(cx, collision);
        let mut reference = forced_grid::<L, Q>(cx, collision);
        for (step, time_s) in [0.25, 0.5, 0.75].into_iter().enumerate() {
            let opening = actual.forces().to_vec();
            let endpoint = if step == 0 {
                opening.clone()
            } else {
                force_epoch(step)
            };
            let before = forced_totals(&actual);
            reference
                .rebind_open_boundaries(
                    forced_targets(&inlet, &outlet, &frames, time_s, &scaling),
                    &gate,
                )
                .unwrap();
            let expected = if step == 0 {
                reference.step_serial(cx, None, [0.0; 3]).unwrap()
            } else {
                reference
                    .advance_serial(cx, None, endpoint.clone(), [0.0; 3])
                    .unwrap()
            };
            let accepted = if step == 0 {
                lattice3::step_forced_serial(
                    &mut actual,
                    controls(&inlet, &outlet),
                    &frames,
                    time_s,
                    &scaling,
                    cx,
                    None,
                    [0.0; 3],
                )
                .unwrap()
            } else {
                lattice3::advance_forced_serial(
                    &mut actual,
                    controls(&inlet, &outlet),
                    &frames,
                    time_s,
                    &scaling,
                    cx,
                    None,
                    endpoint.clone(),
                    [0.0; 3],
                )
                .unwrap()
            };
            assert_eq!(accepted.transport, expected);
            assert_eq!(actual.populations(), reference.populations());
            assert_eq!(actual.forces(), endpoint);
            assert_eq!(actual.boundaries(), reference.boundaries());
            assert_eq!(actual.completed_steps(), step + 1);
            let after = forced_totals(&actual);
            assert!((after.0 - before.0 - accepted.transport.mass_change).abs() < 3e-12);
            for axis in 0..3 {
                let body: f64 = opening
                    .iter()
                    .zip(&endpoint)
                    .map(|(a, b)| 0.5 * a[axis] + 0.5 * b[axis])
                    .sum();
                assert!(
                    (after.1[axis] - before.1[axis] - body - accepted.transport.impulse[axis])
                        .abs()
                        < 3e-13
                );
            }
            // Separate SI oracle: dx/dt=20 m/s, rho=1000 kg/m^3, cs^2=1/3.
            let velocity = [
                (0.1 + 0.2 * time_s) / 20.0,
                (-0.05 - 0.1 * time_s) / 20.0,
                (0.2 + 0.4 * time_s) / 20.0,
            ];
            let density = 1.0 + (100.0 + 200.0 * time_s) / ((1.0 / 3.0) * 1000.0 * 20.0 * 20.0);
            for cell in 0..16 {
                for (value, target) in actual
                    .macroscopic(cell)
                    .unwrap()
                    .1
                    .into_iter()
                    .zip(velocity)
                {
                    near(value, target);
                }
                near(actual.macroscopic(64 + cell).unwrap().0, density);
            }
        }
    });
}

#[test]
fn typed_reservoirs_share_actual_force_epochs_and_impulse_on_both_lattices() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let tau = scaling.flow().tau;
    exercise_forced_scenario::<D3q19, 19>(CollisionModel3::Bgk { tau });
    exercise_forced_scenario::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn late_forced_step_refusal_preserves_sources_reservoirs_and_retry() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let law = CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: scaling.flow().tau,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    };
    let inlet = velocity([0.0, 0.0, 0.4]);
    let outlet = pressure(200.0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut actual = forced_grid::<D3q27, 27>(cx, law);
        let mut reference = forced_grid::<D3q27, 27>(cx, law);
        actual.step_serial(cx, None, [0.0; 3]).unwrap();
        reference.step_serial(cx, None, [0.0; 3]).unwrap();
        let previous = (
            actual.populations().to_vec(),
            actual.forces().to_vec(),
            actual.boundaries(),
            actual.last_wall_exchange(),
        );
        let mut laws = vec![law; 80];
        laws[15].maximum_mach = 0.01;
        // Opening fluid admits the law; the new inlet fails at its last cell
        // after earlier cells have already reconstructed with endpoint forces.
        assert!(matches!(
            lattice3::advance_forced_serial(
                &mut actual,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                Some(&laws),
                force_epoch(1),
                [0.0; 3],
            ),
            Err(lattice3::Error::Boundary(GridError3::Mach { cell: 15, .. }))
        ));
        assert_eq!(actual.populations(), previous.0);
        assert_eq!(actual.forces(), previous.1);
        assert_eq!(actual.boundaries(), previous.2);
        assert_eq!(actual.last_wall_exchange(), previous.3);
        assert_eq!(actual.completed_steps(), 1);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_cx(&cancelled, |stopped| {
            assert!(
                lattice3::step_forced_serial(
                    &mut actual,
                    controls(&inlet, &outlet),
                    &frames,
                    0.5,
                    &scaling,
                    stopped,
                    None,
                    [0.0; 3],
                )
                .is_err()
            )
        });
        assert_eq!(actual.populations(), previous.0);
        assert_eq!(actual.forces(), previous.1);
        assert_eq!(actual.boundaries(), previous.2);
        assert_eq!(actual.completed_steps(), 1);

        reference
            .rebind_open_boundaries(
                forced_targets(&inlet, &outlet, &frames, 0.5, &scaling),
                &gate,
            )
            .unwrap();
        let expected = reference
            .advance_serial(cx, None, force_epoch(1), [0.0; 3])
            .unwrap();
        let retry = lattice3::advance_forced_serial(
            &mut actual,
            controls(&inlet, &outlet),
            &frames,
            0.5,
            &scaling,
            cx,
            None,
            force_epoch(1),
            [0.0; 3],
        )
        .unwrap();
        assert_eq!(retry.transport, expected);
        assert_eq!(actual.populations(), reference.populations());
        assert_eq!(actual.forces(), reference.forces());
        assert_eq!(actual.boundaries(), reference.boundaries());
        assert_eq!(actual.completed_steps(), 2);
    });
}

#[test]
fn generic_velocity_profiles_keep_ownership_across_scenario_adapters() {
    let material = fluid();
    let scaling = plan_scaling_from_material(&material, 0.02, 0.001, 1.0).unwrap();
    let inlet = velocity([0.0; 3]); // Identical to the dormant uniform target.
    let outlet = pressure(100.0);
    let frames = FrameTree::new();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut flow = generic_grid(cx, scaling.flow().tau, 0.2);
        flow.set_velocity_profile(Face3::ZMin, InletFlux3::Unconstrained, &gate, |p| {
            Ok([0.002, -0.001, 0.01 + 0.001 * p[0] as f64])
        })
        .unwrap();
        let profile = flow.velocity_profile(Face3::ZMin).unwrap().to_vec();
        let before = (flow.populations().to_vec(), flow.boundaries());
        assert!(
            fs_couple::scenario_velocity::lattice3::rebind_velocity_face3(
                &mut flow,
                Face3::ZMin,
                &inlet,
                &frames,
                0.5,
                &scaling,
                &gate,
            )
            .is_err()
        );
        assert!(
            lattice3::rebind_open_faces3(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                &gate
            )
            .is_err()
        );
        assert!(
            lattice3::step_serial(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                None,
                [0.0; 3]
            )
            .is_err()
        );
        assert!(
            lattice3::step_pooled(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                &TilePool::for_host(1, 0xB0A0),
                &gate,
                None,
                [0.0; 3]
            )
            .is_err()
        );
        assert_eq!(flow.populations(), before.0);
        assert_eq!(flow.boundaries(), before.1);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(flow.velocity_profile(Face3::ZMin), Some(profile.as_slice()));

        let mut flow =
            ForcedOpenBoundaryGrid3::try_new_preserving_velocity(flow, cx, force_epoch(0)).unwrap();
        let before = (
            flow.populations().to_vec(),
            flow.forces().to_vec(),
            flow.boundaries(),
        );
        assert!(
            lattice3::advance_forced_serial(
                &mut flow,
                controls(&inlet, &outlet),
                &frames,
                0.5,
                &scaling,
                cx,
                None,
                force_epoch(1),
                [0.0; 3]
            )
            .is_err()
        );
        assert_eq!(flow.populations(), before.0);
        assert_eq!(flow.forces(), before.1);
        assert_eq!(flow.boundaries(), before.2);
        let mut outlet_only = controls(&inlet, &outlet);
        outlet_only[4] = None;
        let accepted = lattice3::step_forced_serial(
            &mut flow,
            outlet_only,
            &frames,
            0.5,
            &scaling,
            cx,
            None,
            [0.0; 3],
        )
        .unwrap();
        assert_eq!(accepted.samples[4], None);
        assert_eq!(flow.velocity_profile(Face3::ZMin), Some(profile.as_slice()));
        assert_eq!(flow.forces(), before.1);
        for (cell, target) in profile.into_iter().enumerate() {
            for (actual, expected) in flow.macroscopic(cell).unwrap().1.into_iter().zip(target) {
                near(actual, expected);
            }
        }
    });
}
