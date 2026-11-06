//! G0/G3: one control-time transaction reaches both native boundary operators.

use fs_couple::scenario_open_boundary::{
    BoundarySample3, ScenarioBoundary3, ScenarioBoundaryError, rebind_open_faces3,
};
use fs_couple::scenario_open_boundary::lattice3;
use fs_couple::scenario_pressure::{PressureBoundary3, PressureReference};
use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, Face3, FaceBoundary3};
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use fs_lbm::lattice3::open_boundaries::OpenBoundaryGrid3;
use fs_lbm::lattice3::{Grid3, GridBudget3, GridError3};
use fs_lbm::material::plan_scaling_from_material;
use fs_material::fluid_transport::FluidTransportState;
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, Payload, PayloadId, PayloadMeta, QuantityContract, ReferenceSemantics,
    SampleSource, VectorPayload,
};
use fs_scenario::{BcKind, BcValue, BoundaryCondition, FrameTree, Physics, TimeSignal, WORLD};

const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn velocity(value: [f64; 3]) -> BoundaryCondition {
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
            VectorPayload::new(
                meta,
                SampleSource::Fixed(value.map(|v| QtyAny::new(v, VELOCITY)).to_vec()),
            )
            .unwrap(),
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
    let grid = Grid3::<D3q27, 27>::try_new_equilibrium(
        [4, 4, 5],
        CollisionConfig27 {
            model: CollisionModel27::Bgk {
                tau,
                equilibrium: Equilibrium27::Factorized,
            },
            maximum_mach,
        },
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
