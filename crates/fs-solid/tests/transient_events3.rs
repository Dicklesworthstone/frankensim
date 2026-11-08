//! G1/G3/G4 automatic load knots, exact finite jumps and retained trajectories.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{FrameTree, TimeSignal, WORLD};
use fs_scenario::frame::{Frame, FrameId, FrameMotion};
use fs_scenario::signal::Interp;
use fs_scenario::payload::{OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta,
    QuantityContract, ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload};
use fs_solid::linear3::*;
use fs_solid::loads3::{SurfaceLoad3, TetLoadGeometry3, VolumeLoad3};
use fs_solid::loads3::history::{SurfaceHistory3, TetLoadHistory3, VolumeHistory3, ACCELERATION_DIMS, TRACTION_DIMS};
use fs_solid::loads3::history::dynamics::events::EventLoadHistory3;
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::transient::events::MaterialRebindConfig3;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{MaterialPropertySelection, resolve_isotropic_elastic_state_point};
use fs_qty::{Density, Pressure};

const NODES: [[f64;3];4] = [[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
const TETS: [[usize;4];1] = [[0,1,2,3]];
const ZERO: [[f64;3];4] = [[0.;3];4];
fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free(); if cancelled { gate.request(); }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(&gate,arena,
        StreamKey { seed:83,kernel_id:83,tile:0,iteration:0 },Budget::INFINITE,ExecMode::Deterministic)))
}
fn config() -> TransientConfig3 { TransientConfig3 { maximum_displacement_m:1.0,..Default::default() } }
fn problem(material: &TetElasticMaterial) -> TetLinearElasticProblem<'_> {
    TetLinearElasticProblem { nodes_m:&NODES,tetrahedra:&TETS,materials:TetMaterialField::Uniform(material),
        fixed_dofs:&[],budget:TetAssemblyBudget::default() }
}
fn payload(times: &[f64], values: &[[f64;3]], interp: TableInterpolation,
    outside: OutsideDomainPolicy, frame: FrameId) -> VectorPayload {
    let meta = PayloadMeta::new(QuantityContract::Dimensions(ACCELERATION_DIMS),
        PayloadId::new(if frame == WORLD { "basis/world-cartesian" } else { "basis/frame-cartesian" }).unwrap(),
        frame,OrientationParity::Even,ReferenceSemantics::Continuous).unwrap();
    VectorPayload::new(meta,SampleSource::table(
        times.iter().map(|&t| QtyAny::new(t,Dims([0,0,1,0,0,0]))).collect(),
        values.iter().map(|v| v.map(|x| QtyAny::new(x,ACCELERATION_DIMS)).to_vec()).collect(),interp,outside).unwrap()).unwrap()
}
#[test]
fn chirped_frame_and_pressure_sweep_knots_are_retained() {
    with_cx(false, |cx| {
        let chirp = |start, end, dims| TimeSignal::Chirp(
            fs_scenario::signal::LinearChirp::new(start, end, 1.0, 4.0, 0.2,
                QtyAny::new(0.0, dims), QtyAny::new(0.01, dims)).unwrap());
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1), parent: WORLD, name: "chirped-drive".into(),
            motion: FrameMotion::Translating {
                direction: [1.0, 0.0, 0.0],
                displacement: chirp(0.25, 0.75, Dims([1, 0, 0, 0, 0, 0])),
            },
        });
        let pressure = chirp(0.1, 0.9, TRACTION_DIMS);
        let surfaces = [SurfaceHistory3::Pressure { faces: &[0], value: &pressure }];
        let history = TetLoadHistory3 {
            surfaces: &surfaces, volume: VolumeHistory3::Fixed(VolumeLoad3::None),
            nodal_forces_n: &[], frames: &frames,
        };
        let schedule = EventLoadHistory3::new(history, 16, cx).unwrap();
        assert_eq!(schedule.knots_s(), &[0.1, 0.25, 0.75, 0.9]);
    });
}

#[test]
fn g1_chirped_pressure_cannot_alias_whole_cycles_in_native_fixed_or_adaptive_steps() {
    use fs_solid::loads3::history::dynamics::events::adaptive::TransientTimeControl3;

    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.25, [83; 32]).unwrap();
        let p = problem(&material);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let frequency = 4.0;
        let omega = core::f64::consts::TAU * frequency;
        let amplitude = 0.02;
        let pressure = TimeSignal::Chirp(fs_scenario::LinearChirp::new(
            0.0, 1.0, frequency, frequency, 0.0,
            QtyAny::new(0.0, TRACTION_DIMS), QtyAny::new(amplitude, TRACTION_DIMS),
        ).unwrap());
        let surfaces = [SurfaceHistory3::Pressure { faces: &[0], value: &pressure }];
        let history = TetLoadHistory3 {
            surfaces: &surfaces, volume: VolumeHistory3::Fixed(VolumeLoad3::None),
            nodal_forces_n: &[], frames: &frames,
        };
        let schedule = EventLoadHistory3::new(history, 16, cx).unwrap();
        // The unconstrained body's mass is rho*V = 6*(1/6) = 1 kg.
        // Internal elastic forces cancel in its center-of-mass equation.
        let face = geometry.faces()[0];
        let acceleration = face.outward_normal.map(|n| -amplitude * face.area_m2 * n);
        let control = TransientTimeControl3 {
            maximum_step_s: 1.0, minimum_step_s: 1e-10, ..Default::default()
        };
        for adaptive in [false, true] {
            let mut run = schedule.begin_transient(
                &geometry, &p, None, 0.0, &ZERO, &ZERO, config(), cx,
            ).unwrap();
            // The original full and half steps at 0, .25, .5 all sampled
            // cosine peaks, incorrectly treating two periods as a constant load.
            while run.run().state().time_s < 0.5 {
                let start = run.run().state().time_s;
                let end = if adaptive {
                    run.advance_adaptive_to(0.5, 1.0, control, cx).unwrap().second.end_time_s
                } else {
                    run.advance_to(0.5, 1.0, cx).unwrap().end_time_s
                };
                assert!((end - start) * frequency <= fs_scenario::DEFAULT_MAXIMUM_CYCLES_PER_STEP);
                let motion = run.run().nodal_motion(cx).unwrap();
                for axis in 0..3 {
                    let q = motion[0].iter().map(|v| v[axis]).sum::<f64>() / 4.0;
                    let v = motion[1].iter().map(|v| v[axis]).sum::<f64>() / 4.0;
                    let expected_q = acceleration[axis] * (1.0 - (omega * end).cos()) / (omega * omega);
                    let expected_v = acceleration[axis] * (omega * end).sin() / omega;
                    // Sixteen steps/cycle bounds this fixture's coarse Newmark
                    // phase error; adaptive error control may refine further.
                    assert!((q - expected_q).abs() < 0.06 * acceleration[axis].abs() / (omega * omega) + 1e-10,
                        "center displacement {q:e} != {expected_q:e} at {end:e}");
                    assert!((v - expected_v).abs() < 0.02 * acceleration[axis].abs() / omega + 1e-10,
                        "center velocity {v:e} != {expected_v:e} at {end:e}");
                }
            }
            assert!(run.run().state().step >= 32);
        }
    });
}

#[test]
fn g0_chirped_frame_frequency_limits_the_actual_event_step_and_refusal_is_atomic() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.25, [83; 32]).unwrap();
        let p = problem(&material);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1), parent: WORLD, name: "oscillating load frame".into(),
            motion: FrameMotion::Tilt {
                axis: [0.0, 0.0, 1.0], center: fs_ga::Vec3::new(0.0, 0.0, 0.0),
                angle: TimeSignal::Chirp(fs_scenario::LinearChirp::new(
                    0.0, 1.0, 32.0, 32.0, 0.0,
                    QtyAny::dimensionless(0.0), QtyAny::dimensionless(0.1),
                ).unwrap()),
            },
        });
        let acceleration = payload(&[0.0, 1.0], &[[0.01, 0.0, 0.0]; 2],
            TableInterpolation::Linear, OutsideDomainPolicy::Clamp, FrameId(1));
        let history = TetLoadHistory3 {
            surfaces: &[], volume: VolumeHistory3::Acceleration(&acceleration),
            nodal_forces_n: &[], frames: &frames,
        };
        let schedule = EventLoadHistory3::new(history, 16, cx).unwrap();
        let mut run = schedule.begin_transient(
            &geometry, &p, None, 0.0, &ZERO, &ZERO, config(), cx,
        ).unwrap();
        let step = run.advance_to(0.5, 0.5, cx).unwrap();
        assert!(step.end_time_s > 0.0 && step.end_time_s <= 1.0 / (16.0 * 32.0));
        let saved = run.run().state().clone();
        with_cx(true, |cancelled| {
            assert!(run.advance_to(0.5, 0.5, cancelled).is_err());
        });
        assert_eq!(run.run().state(), &saved);
    });
}

#[test]
fn g1_unused_chirped_frames_keep_knots_without_throttling_active_loads() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.25, [83; 32]).unwrap();
        let p = problem(&material);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(7), parent: WORLD, name: "unused fast tilt".into(),
            motion: FrameMotion::Tilt {
                axis: [0.0, 0.0, 1.0], center: fs_ga::Vec3::new(0.0, 0.0, 0.0),
                angle: TimeSignal::Chirp(fs_scenario::LinearChirp::new(
                    0.125, 0.75, 1e9, 1e9, 0.0,
                    QtyAny::dimensionless(0.0), QtyAny::dimensionless(0.1),
                ).unwrap()),
            },
        });
        let acceleration = payload(&[0.0, 1.0], &[[0.01, 0.0, 0.0]; 2],
            TableInterpolation::Linear, OutsideDomainPolicy::Clamp, WORLD);
        for driven_pressure in [false, true] {
            let pressure = if driven_pressure {
                TimeSignal::Chirp(fs_scenario::LinearChirp::new(
                    0.0, 1.0, 4.0, 4.0, 0.0,
                    QtyAny::new(0.0, TRACTION_DIMS), QtyAny::new(0.001, TRACTION_DIMS),
                ).unwrap())
            } else {
                TimeSignal::Constant(QtyAny::new(0.0, TRACTION_DIMS))
            };
            let surfaces = [SurfaceHistory3::Pressure { faces: &[0], value: &pressure }];
            let history = TetLoadHistory3 {
                surfaces: &surfaces, volume: VolumeHistory3::Acceleration(&acceleration),
                nodal_forces_n: &[], frames: &frames,
            };
            let schedule = EventLoadHistory3::new(history, 32, cx).unwrap();
            assert_eq!(schedule.knots_s(), &[0.0, 0.125, 0.75, 1.0]);
            let mut run = schedule.begin_transient(
                &geometry, &p, None, 0.0, &ZERO, &ZERO, config(), cx,
            ).unwrap();
            let step = run.advance_to(0.5, 0.5, cx).unwrap();
            if driven_pressure {
                // The directly used pressure still resolves 4 Hz, not the
                // unrelated GHz frame and not just its own authored knots.
                assert!(step.end_time_s > 1.0 / 128.0);
                assert!(step.end_time_s <= 1.0 / 64.0);
            } else {
                assert_eq!(step.end_time_s, 0.125);
                assert_eq!(run.advance_to(0.5, 0.5, cx).unwrap().end_time_s, 0.5);
            }
        }
    });
}

#[test]
fn g1_used_ancestor_tilt_resolves_angular_travel_and_refuses_a_coarse_clock() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.25, [83; 32]).unwrap();
        let p = problem(&material);
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frequency = 2.0;
        let amplitude = 100.0;
        let angle = TimeSignal::Chirp(fs_scenario::LinearChirp::new(
            0.0, 1.0, frequency, frequency, core::f64::consts::FRAC_PI_2,
            QtyAny::dimensionless(0.0), QtyAny::dimensionless(amplitude),
        ).unwrap());
        let mut frames = FrameTree::new();
        // The payload names a fixed child; its dynamic parent is deliberately
        // stored later so frame ancestry must follow identities, not row order.
        frames.add(Frame {
            id: FrameId(2), parent: FrameId(1), name: "load coordinates".into(),
            motion: FrameMotion::Fixed {
                orientation: fs_ga::Quat::identity(), translation: fs_ga::Vec3::new(0.0, 0.0, 0.0),
            },
        });
        frames.add(Frame {
            id: FrameId(1), parent: WORLD, name: "large angular excursion".into(),
            motion: FrameMotion::Tilt {
                axis: [0.0, 0.0, 1.0], center: fs_ga::Vec3::new(0.0, 0.0, 0.0), angle: angle.clone(),
            },
        });
        let acceleration = payload(&[0.0, 1.0], &[[0.01, 0.0, 0.0]; 2],
            TableInterpolation::Linear, OutsideDomainPolicy::Clamp, FrameId(2));
        let history = TetLoadHistory3 {
            surfaces: &[], volume: VolumeHistory3::Acceleration(&acceleration),
            nodal_forces_n: &[], frames: &frames,
        };
        let schedule = EventLoadHistory3::new(history, 16, cx).unwrap();
        let mut run = schedule.begin_transient(
            &geometry, &p, None, 0.0, &ZERO, &ZERO, config(), cx,
        ).unwrap();
        let step = run.advance_to(0.25, 0.25, cx).unwrap();
        let allowance = fs_scenario::DEFAULT_MAXIMUM_CYCLES_PER_STEP;
        assert!(step.end_time_s > 0.0 && step.end_time_s * frequency * amplitude <= allowance);
        let travel = (angle.eval(step.end_time_s).unwrap().value - angle.eval(0.0).unwrap().value).abs();
        assert!(travel <= core::f64::consts::TAU * allowance);

        // A representable source-phase step alone is too coarse for this
        // tilt's angular speed on a remote clock. Refuse without publishing.
        let remote = 2.0_f64.powi(44);
        let mut run = schedule.begin_transient(
            &geometry, &p, None, remote, &ZERO, &ZERO, config(), cx,
        ).unwrap();
        let saved = run.run().state().clone();
        assert!(run.advance_to(remote + 1.0, 1.0, cx).is_err());
        assert_eq!(run.run().state(), &saved);
    });
}

#[test]
fn translating_frame_history_knots_are_retained_by_the_native_event_owner() {
    with_cx(false, |cx| {
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1), parent: WORLD, name: "linear-drive".into(),
            motion: FrameMotion::Translating {
                direction: [1.0, 0.0, 0.0],
                displacement: TimeSignal::Table {
                    times: vec![0.0, 0.25, 0.75], values: vec![0.0, 0.1, 0.0],
                    dims: Dims([1, 0, 0, 0, 0, 0]), interp: Interp::Hold,
                },
            },
        });
        let history = TetLoadHistory3 {
            surfaces: &[], volume: VolumeHistory3::Fixed(VolumeLoad3::None),
            nodal_forces_n: &[], frames: &frames,
        };
        let schedule = EventLoadHistory3::new(history, 16, cx).unwrap();
        assert_eq!(schedule.knots_s(), &[0.0, 0.25, 0.75]);
    });
}
#[test]
fn g1_short_held_pulse_has_no_pre_response_and_exact_momentum() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let acceleration=payload(&[0.,0.25,0.375,1.],&[[0.;3],[4.,0.,0.],[0.;3],[0.;3]],
            TableInterpolation::StepLeft,OutsideDomainPolicy::Clamp,WORLD);
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,64,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let first=run.advance_to(0.75,1.,cx).unwrap();
        assert_eq!(first.end_time_s,0.25); assert!(first.jump.is_some());
        assert_eq!(run.run().nodal_motion(cx).unwrap()[0],ZERO);
        assert_eq!(run.run().nodal_motion(cx).unwrap()[1],ZERO);
        let second=run.advance_to(0.75,1.,cx).unwrap(); assert_eq!(second.end_time_s,0.375);
        assert!(second.jump.is_some()); assert!((second.interval.mechanical_work_j-0.125).abs()<1e-12);
        let third=run.advance_to(0.75,1.,cx).unwrap(); assert_eq!(third.end_time_s,0.75);
        let motion=run.run().nodal_motion(cx).unwrap();
        for node in 0..4 {
            assert!((motion[0][node][0]-0.21875).abs()<1e-10);
            assert!((motion[1][node][0]-0.5).abs()<1e-10);
        }
        assert_eq!(run.run().state().step,3);
    });
}
#[test]
fn g3_held_pressure_uses_the_exact_left_load_then_updates_reactions() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let faces=[0usize];
        let pressure=TimeSignal::Table { times:vec![0.,0.125,0.25],values:vec![0.,0.02,0.],dims:TRACTION_DIMS,interp:Interp::Hold };
        let surfaces=[SurfaceHistory3::Pressure { faces:&faces,value:&pressure }];
        let history=TetLoadHistory3 { surfaces:&surfaces,volume:VolumeHistory3::Fixed(VolumeLoad3::None),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,32,cx).unwrap();
        let empty=geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let loaded=geometry.assemble(&[SurfaceLoad3::Pressure { faces:&faces,pascals:0.02 }],VolumeLoad3::None,&[],cx).unwrap();
        let mut manual=empty.begin_transient(&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        manual.advance_to_load_jump(&empty,&loaded,0.125,cx).unwrap(); run.advance_to(0.5,0.5,cx).unwrap();
        assert_eq!(run.run().state(),manual.state());
        manual.advance_to_load_jump(&loaded,&empty,0.125,cx).unwrap(); run.advance_to(0.5,0.5,cx).unwrap();
        assert_eq!(run.run().state(),manual.state());
    });
}
#[test]
fn g1_scheduled_frame_tilt_changes_direction_without_moving_the_clock() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let mut frames=FrameTree::new();
        frames.frames.push(Frame { id:FrameId(1),name:"tilting load frame".into(),parent:WORLD,
            motion:FrameMotion::Tilt { axis:[0.,0.,1.],center:fs_ga::Vec3::new(10.,20.,30.),
                angle:TimeSignal::Table { times:vec![0.,0.25,0.5],values:vec![0.,core::f64::consts::FRAC_PI_2,core::f64::consts::FRAC_PI_2],dims:Dims::NONE,interp:Interp::Hold } } });
        let acceleration=payload(&[0.,0.125],&[[0.2,0.,0.];2],TableInterpolation::Linear,OutsideDomainPolicy::Clamp,FrameId(1));
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,64,cx).unwrap();
        assert_eq!(schedule.knots_s(),&[0.,0.125,0.25,0.5]);
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        while run.run().state().time_s<0.5 { run.advance_to(0.5,1.,cx).unwrap(); }
        let motion=run.run().nodal_motion(cx).unwrap();
        for node in 0..4 {
            for c in 0..2 { assert!((motion[1][node][c]-0.05).abs()<1e-10); }
            assert!((motion[0][node][0]-0.01875).abs()<1e-10);
            assert!((motion[0][node][1]-0.00625).abs()<1e-10);
        }
    });
}
#[test]
fn g4_right_limit_failure_rolls_back_the_interval_and_checkpoint_replays() {
    with_cx(false,|cx| {
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap(); let frames=FrameTree::new();
        let acceleration=payload(&[0.,0.25,0.5],&[[0.;3],[f64::MAX,0.,0.],[0.;3]],TableInterpolation::StepLeft,OutsideDomainPolicy::Clamp,WORLD);
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,32,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let before=run.run().state().clone();
        assert!(run.advance_to(1.,1.,cx).is_err()); assert_eq!(run.run().state(),&before);
        with_cx(true,|cancelled| { assert!(run.advance_to(0.125,1.,cancelled).is_err()); });
        assert_eq!(run.run().state(),&before);
        let mut fork=run.clone();
        assert_eq!(run.advance_to(0.125,1.,cx).unwrap(),fork.advance_to(0.125,1.,cx).unwrap());
        assert_eq!(run.run().state(),fork.run().state());
        assert!(run.advance_to(1.,1.,cx).is_err()); assert_eq!(run.run().state(),fork.run().state());
    });
}
#[test]
fn g0_schedule_budget_periodic_tables_and_malformed_scalars_refuse() {
    with_cx(false,|cx| {
        let frames=FrameTree::new();
        let values=[[0.;3];2];
        let periodic=payload(&[0.,1.],&values,TableInterpolation::Linear,OutsideDomainPolicy::Periodic,WORLD);
        let base=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&periodic),nodal_forces_n:&[],frames:&frames };
        assert!(EventLoadHistory3::new(base,32,cx).is_err());
        let finite=payload(&[0.,1.],&values,TableInterpolation::Linear,OutsideDomainPolicy::Clamp,WORLD);
        let base=TetLoadHistory3 { volume:VolumeHistory3::Acceleration(&finite),..base };
        for cap in [0,7,1_048_577] { assert!(EventLoadHistory3::new(base,cap,cx).is_err()); }
        assert!(EventLoadHistory3::new(base,8,cx).is_ok());
        let malformed=TimeSignal::Table { times:vec![0.,1.],values:vec![0.],dims:TRACTION_DIMS,interp:Interp::Hold };
        let surfaces=[SurfaceHistory3::Pressure { faces:&[0],value:&malformed }];
        assert!(EventLoadHistory3::new(TetLoadHistory3 { surfaces:&surfaces,..base },32,cx).is_err());
    });
}
#[test]
fn g0_adjacent_knots_remain_distinct_and_continuous_knots_need_no_jump() {
    with_cx(false,|cx| {
        let frames=FrameTree::new();
        let a=1.0_f64; let b=a.next_up();
        let acceleration=payload(&[0.,a,b],&[[0.;3];3],TableInterpolation::Linear,OutsideDomainPolicy::Clamp,WORLD);
        let history=TetLoadHistory3 { surfaces:&[],volume:VolumeHistory3::Acceleration(&acceleration),nodal_forces_n:&[],frames:&frames };
        let schedule=EventLoadHistory3::new(history,12,cx).unwrap(); assert_eq!(schedule.knots_s(),&[0.,a,b]);
        let material=TetElasticMaterial::try_new(6.,60.,0.25,[83;32]).unwrap(); let p=problem(&material);
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let mut run=schedule.begin_transient(&geometry,&p,None,0.,&ZERO,&ZERO,config(),cx).unwrap();
        let first=run.advance_to(2.,2.,cx).unwrap(); assert_eq!(first.end_time_s,a); assert!(first.jump.is_none());
        let second=run.advance_to(2.,2.,cx).unwrap(); assert_eq!(second.end_time_s,b); assert!(second.jump.is_none());
        let before=run.run().state().clone();
        for cap in [0.,f64::NAN,1e-300] { assert!(run.advance_to(2.,cap,cx).is_err()); assert_eq!(run.run().state(),&before); }
    });
}

fn rebind_controls() -> MaterialRebindConfig3 {
    MaterialRebindConfig3 { maximum_absolute_work_j:1.0,
        absolute_momentum_tolerance_kg_m_s:1e-13, maximum_strain:0.1 }
}
fn source_material(density: f64, young: f64) -> TetElasticMaterial {
    let mut claims=ClaimSet::new();
    for (name,dims,value) in [("density",Density::DIMS,density),
        ("young_modulus",Pressure::DIMS,young),("poisson_ratio",Dims::NONE,0.0)] {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained(), interpolation:InterpolationPolicy::ConstantWithinValidity,
            uncertainty:UncertaintyModel::Unstated, observations:vec![],
            provenance:Provenance {source:"synthetic native material rebind".into(),
                license:"test-only".into(),artifact:None} }).unwrap();
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"solid".into(),
        process:"rebind regression".into(),revision:0},claims,vec![]).unwrap();
    let point=QueryPoint::new().with("T",300.0).unwrap();
    let resolved=resolve_isotropic_elastic_state_point(&card,&point,MaterialPropertySelection::SingleClaimOnly).unwrap();
    TetElasticMaterial::from_resolved_elastic_state(&resolved)
}

#[test]
fn g1_source_material_rebind_preserves_momentum_and_changes_native_motion_and_stress() {
    with_cx(false,|cx| {
        let original=source_material(6.0,60.0); let replacement=source_material(12.0,180.0);
        let fixed:Vec<_>=(0..12).filter(|d| *d!=3).collect();
        let before_problem=TetLinearElasticProblem {fixed_dofs:&fixed,..problem(&original)};
        let after_problem=TetLinearElasticProblem {fixed_dofs:&fixed,..problem(&replacement)};
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,before_problem.budget,4,cx).unwrap();
        let old_load=geometry.assemble(&[],VolumeLoad3::None,&[(3,0.1)],cx).unwrap();
        let new_load=geometry.assemble(&[],VolumeLoad3::None,&[(3,0.25)],cx).unwrap();
        let frame=fs_blake3::ContentHash([84;32]);
        let strain=TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            replacement.material_state_identity,frame,0.1).unwrap();
        let mut u=ZERO; u[1][0]=0.01;
        let mut v=ZERO; v[1][0]=0.2;
        let cfg=TransientConfig3 {mass_damping_per_s:0.4,stiffness_damping_s:0.02,..config()};
        let run=old_load.begin_transient(&before_problem,None,2.0,&u,&v,cfg,cx).unwrap();
        let saved=run.state().clone();
        let (mut replaced,report)=run.rebind_material_candidate(&after_problem,
            Some(TetThermalStrainField::Uniform(&strain)),&new_load,frame,rebind_controls(),cx).unwrap();
        // Independent one-DOF mechanics: (m,k) changes (.1,10) -> (.2,30).
        // Momentum .02 is unchanged; new C=.4*.2+.02*30=.68.
        assert_eq!(replaced.state().displacement,vec![0.01]);
        assert!((replaced.state().velocity[0]-0.1).abs()<1e-14);
        assert!((replaced.state().acceleration[0]+0.44).abs()<1e-12);
        assert_eq!(replaced.state().time_s,2.0); assert_eq!(replaced.state().step,0);
        assert!((report.kinetic_work_j+0.001).abs()<1e-14);
        assert!((report.elastic_work_j-0.000715).abs()<1e-14);
        assert!((report.supplied_work_j+0.000285).abs()<1e-14);
        assert!((report.reference_mass_kg[0]-1.0).abs()<1e-14);
        assert!((report.reference_mass_kg[1]-2.0).abs()<1e-14);
        assert!(report.reduced_momentum_residual_kg_m_s<1e-13);
        assert_eq!(report.material_state_identities,vec![[original.material_state_identity,
            replacement.material_state_identity]]);
        assert_eq!(report.thermal_state_identities,vec![[None,Some(strain.identity())]]);
        let response=replaced.recover_response(frame,0.1,cx).unwrap();
        assert!((response.stress.elements()[0].stress_mandel_pa()[0]-1.62).abs()<1e-12);
        assert!(response.free_force_residual_n<1e-12);
        assert_eq!(run.state(),&saved);

        // The refreshed K/M, finite force and re-equilibrated a launch the next interval.
        replaced.advance(&new_load,0.02,cx).unwrap();
        let qp=0.01+0.02*0.1-0.0001*0.44;
        let vp=0.1-0.01*0.44;
        let acceleration=(0.25+0.03-30.0*qp-0.68*vp)/(0.2+0.01*0.68+0.0001*30.0);
        assert!((replaced.state().displacement[0]-(qp+0.0001*acceleration)).abs()<1e-13);
        let retained=replaced.state().clone();
        let zero_work=MaterialRebindConfig3 {maximum_absolute_work_j:0.0,..rebind_controls()};
        let (same,unchanged)=replaced.rebind_material_candidate(&after_problem,
            Some(TetThermalStrainField::Uniform(&strain)),&new_load,frame,zero_work,cx).unwrap();
        assert_eq!(same.state().displacement,retained.displacement);
        assert_eq!(same.state().velocity,retained.velocity);
        assert_eq!(same.state().time_s,retained.time_s); assert_eq!(same.state().step,1);
        assert_eq!(unchanged.supplied_work_j,0.0);
        assert_eq!(unchanged.step,1);

        // A late work refusal, after the mass/acceleration/response solves,
        // cannot replace the accepted source, history or motion.
        assert!(replaced.rebind_material_candidate(&before_problem,None,&new_load,
            frame,zero_work,cx).is_err());
        let (back,reverse)=replaced.rebind_material_candidate(&before_problem,None,&new_load,
            frame,rebind_controls(),cx).unwrap();
        assert_eq!(back.state().time_s,retained.time_s); assert_eq!(back.state().step,1);
        assert!((back.state().velocity[0]-2.0*retained.velocity[0]).abs()<1e-13);
        assert!((reverse.kinetic_energy_j[1]-2.0*reverse.kinetic_energy_j[0]).abs()<1e-13);
        let changed_supports=problem(&original);
        assert!(replaced.rebind_material_candidate(&changed_supports,None,&new_load,
            frame,rebind_controls(),cx).is_err());
        with_cx(true,|cancelled| {
            assert!(replaced.rebind_material_candidate(&before_problem,None,&new_load,
                frame,rebind_controls(),cancelled).is_err());
        });
        assert_eq!(replaced.state(),&retained);
        assert!(std::ptr::eq(replaced.material_problem(),&after_problem));
    });
}

#[test]
fn g1_clamped_material_rebind_accounts_eigenstrain_storage_without_motion() {
    with_cx(false,|cx| {
        let original=source_material(6.0,60.0); let replacement=source_material(12.0,180.0);
        let fixed:Vec<_>=(0..12).collect();
        let before_problem=TetLinearElasticProblem {fixed_dofs:&fixed,..problem(&original)};
        let after_problem=TetLinearElasticProblem {fixed_dofs:&fixed,..problem(&replacement)};
        let geometry=TetLoadGeometry3::new(&NODES,&TETS,before_problem.budget,4,cx).unwrap();
        let load=geometry.assemble(&[],VolumeLoad3::None,&[],cx).unwrap();
        let run=load.begin_transient(&before_problem,None,0.0,&ZERO,&ZERO,config(),cx).unwrap();
        let frame=fs_blake3::ContentHash([84;32]);
        let strain=TetThermalStrainState::try_new(310.0,300.0,[0.001,0.0,0.0,0.0,0.0,0.0],
            replacement.material_state_identity,frame,0.1).unwrap();
        let (next,report)=run.rebind_material_candidate(&after_problem,
            Some(TetThermalStrainField::Uniform(&strain)),&load,frame,rebind_controls(),cx).unwrap();
        assert_eq!(report.kinetic_work_j,0.0);
        assert!((report.elastic_work_j-15e-6).abs()<1e-15);
        assert!((report.supplied_work_j-15e-6).abs()<1e-15);
        assert_eq!(report.reduced_momentum_residual_kg_m_s,0.0);
        assert_eq!(next.nodal_motion(cx).unwrap(),[ZERO.to_vec(),ZERO.to_vec(),ZERO.to_vec()]);
        let limited=MaterialRebindConfig3 {maximum_absolute_work_j:10e-6,..rebind_controls()};
        assert!(run.rebind_material_candidate(&after_problem,Some(TetThermalStrainField::Uniform(&strain)),
            &load,frame,limited,cx).is_err());
        let excessive=TetThermalStrainState::try_new(310.0,300.0,[0.2,0.0,0.0,0.0,0.0,0.0],
            replacement.material_state_identity,frame,0.3).unwrap();
        assert!(run.rebind_material_candidate(&after_problem,Some(TetThermalStrainField::Uniform(&excessive)),
            &load,frame,rebind_controls(),cx).is_err());
        assert!(run.thermal_strain_field().is_none());
        assert_eq!(run.state().time_s,0.0); assert_eq!(run.state().step,0);
    });
}

#[test]
fn g1_g3_material_work_survives_unchanged_consistent_mass_inventory() {
    with_cx(false, |cx| {
        // Two disconnected P1 elements each retain one free x component.
        // Their consistent reduced masses are 1 and 2^54 kg. Every velocity
        // stays at or below 1 m/s; the extreme mass is a numerical stress case.
        let background = 2.0_f64.powi(54);
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3],
            [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.0, 0.0, 1.0]];
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let old_materials = [source_material(60.0, 60.0),
            source_material(60.0 * background, 60.0)];
        let new_materials = [source_material(120.0, 60.0),
            source_material(60.0 * background, 60.0)];
        let fixed: Vec<_> = (0..24).filter(|&d| d != 3 && d != 15).collect();
        let old_problem = TetLinearElasticProblem { nodes_m: &nodes, tetrahedra: &tets,
            materials: TetMaterialField::PerElement(&old_materials), fixed_dofs: &fixed,
            budget: Default::default() };
        let new_problem = TetLinearElasticProblem {
            materials: TetMaterialField::PerElement(&new_materials), ..old_problem };
        let geometry = TetLoadGeometry3::new(&nodes, &tets, old_problem.budget, 8, cx).unwrap();
        let load = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let zero = [[0.0; 3]; 8];
        let mut velocity = zero;
        velocity[1][0] = 1.0;
        velocity[5][0] = 1.0;
        let run = load.begin_transient(&old_problem, None, 3.0, &zero, &velocity,
            config(), cx).unwrap();
        let accepted = run.state().clone();
        let frame = fs_blake3::ContentHash([85; 32]);
        let (next, report) = run.rebind_material_candidate(&new_problem, None, &load,
            frame, rebind_controls(), cx).unwrap();
        assert_eq!(next.state().velocity, vec![0.5, 1.0]);
        assert_eq!(report.kinetic_energy_j[0], report.kinetic_energy_j[1]);
        assert_eq!(report.kinetic_work_j, -0.25);
        assert_eq!(report.retained_kinetic_work_j(), fs_math::dd::Dd::from_f64(-0.25));
        assert_eq!(report.supplied_work_j, -0.25);
        assert_eq!(report.reduced_momentum_residual_kg_m_s, 0.0);
        let zero_work = MaterialRebindConfig3 { maximum_absolute_work_j: 0.0,
            ..rebind_controls() };
        assert!(run.rebind_material_candidate(&new_problem, None, &load,
            frame, zero_work, cx).is_err());
        assert_eq!(run.state(), &accepted);
        let (restored, reverse) = next.rebind_material_candidate(&old_problem, None,
            &load, frame, rebind_controls(), cx).unwrap();
        assert_eq!(restored.state().velocity, accepted.velocity);
        assert_eq!(reverse.retained_supplied_work_j(), fs_math::dd::Dd::from_f64(0.25));
        assert_eq!(restored.state().time_s, 3.0);
        assert_eq!(restored.state().step, 0);
    });
}
