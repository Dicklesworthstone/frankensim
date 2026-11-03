//! G1/G3/G4: typed displacement sources drive actual 3-D inertial supports.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::Vec3;
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, OutsideDomainPolicy, PayloadId, PayloadMeta, QuantityContract,
    ReferenceSemantics, SampleSource, TableInterpolation, VectorPayload,
};
use fs_scenario::{Frame, FrameId, FrameMotion, FrameTree, WORLD};
use fs_solid::linear3::*;
use fs_solid::loads3::history::supports::{DISPLACEMENT_DIMS, SupportMotion3};
use fs_solid::loads3::history::{TRACTION_DIMS, TetLoadHistory3, VolumeHistory3};
use fs_solid::loads3::history::ACCELERATION_DIMS;
use fs_solid::loads3::history::dynamics::events::EventLoadHistory3;
use fs_solid::loads3::history::dynamics::events::adaptive::{AdaptiveHistoryError3, TransientTimeControl3};
use fs_solid::loads3::transient::TransientConfig3;
use fs_solid::loads3::{TetLoadGeometry3, VolumeLoad3};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ZERO: [[f64; 3]; 4] = [[0.0; 3]; 4];
const ALL: [usize; 4] = [0, 1, 2, 3];
const FIXED: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 75,
                kernel_id: 75,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn source(
    dims: Dims,
    frame: FrameId,
    a: [f64; 3],
    b: [f64; 3],
    outside: OutsideDomainPolicy,
) -> VectorPayload {
    let basis = if frame == WORLD {
        "basis/world-cartesian"
    } else {
        "basis/frame-cartesian"
    };
    VectorPayload::new(
        PayloadMeta::new(
            QuantityContract::Dimensions(dims),
            PayloadId::new(basis).unwrap(),
            frame,
            OrientationParity::Even,
            ReferenceSemantics::Continuous,
        )
        .unwrap(),
        SampleSource::table(
            vec![
                QtyAny::new(0.0, Dims([0, 0, 1, 0, 0, 0])),
                QtyAny::new(1.0, Dims([0, 0, 1, 0, 0, 0])),
            ],
            vec![
                a.map(|v| QtyAny::new(v, dims)).to_vec(),
                b.map(|v| QtyAny::new(v, dims)).to_vec(),
            ],
            TableInterpolation::Linear,
            outside,
        )
        .unwrap(),
    )
    .unwrap()
}

fn problem<'a>(
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
) -> TetLinearElasticProblem<'a> {
    TetLinearElasticProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::default(),
    }
}

fn history(frames: &FrameTree) -> TetLoadHistory3<'_> {
    TetLoadHistory3 {
        surfaces: &[],
        volume: VolumeHistory3::Fixed(VolumeLoad3::None),
        nodal_forces_n: &[],
        frames,
    }
}

fn table(dims: Dims, times: &[f64], x: &[f64], interpolation: TableInterpolation) -> VectorPayload {
    VectorPayload::new(
        PayloadMeta::new(QuantityContract::Dimensions(dims),
            PayloadId::new("basis/world-cartesian").unwrap(), WORLD,
            OrientationParity::Even, ReferenceSemantics::Continuous).unwrap(),
        SampleSource::table(
            times.iter().map(|&t| QtyAny::new(t, Dims([0,0,1,0,0,0]))).collect(),
            x.iter().map(|&x| [x,0.,0.].map(|v| QtyAny::new(v,dims)).to_vec()).collect(),
            interpolation, OutsideDomainPolicy::Clamp).unwrap()).unwrap()
}

#[test]
fn g1_support_knots_and_force_pulses_share_one_inertial_clock_and_work_balance() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.,60.,0.,[75;32]).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let p = problem(&material, &fixed);
        let g = TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let frames = FrameTree::new();
        let rate = 1.0/1024.0;
        let pulse = 1.0/128.0;
        let motion = table(DISPLACEMENT_DIMS, &[0.,0.0625,0.5],
            &[0.,rate*0.0625,rate*0.5], TableInterpolation::Linear);
        let acceleration = table(ACCELERATION_DIMS, &[0.,0.125,0.25,0.5],
            &[0.,pulse,0.,0.], TableInterpolation::StepLeft);
        let rows = [SupportMotion3 { nodes:&[0], components:[true,false,false], displacement:&motion }];
        let h = TetLoadHistory3 { volume:VolumeHistory3::Acceleration(&acceleration), ..history(&frames) };
        let schedule = EventLoadHistory3::new_with_supports(h,&rows,28,cx).unwrap();
        assert_eq!(schedule.knots_s(), &[0.,0.0625,0.125,0.25,0.5]);
        assert!(EventLoadHistory3::new_with_supports(h,&rows,27,cx).is_err());
        assert!(schedule.begin_transient(&g,&p,None,0.,&ZERO,&ZERO,Default::default(),cx).is_err());
        let mut velocity = ZERO; velocity[0][0] = rate;
        let mut errors = Vec::new();
        for maximum_step in [1.0/32.0,1.0/64.0,1.0/128.0] {
            let mut run = schedule.begin_transient_with_supports(&g,&p,None,0.,&ZERO,
                &velocity,&ZERO,Default::default(),cx).unwrap();
            let initial_energy = run.run().state().total_energy_j;
            let mut work = 0.; let mut jumps = Vec::new(); let mut visited = Vec::new();
            while run.run().state().time_s < 0.5 {
                let report = run.advance_to(0.5,maximum_step,cx).unwrap();
                assert!(report.interval.balance_error_j.abs() < 1e-11);
                work += report.interval.mechanical_work_j + report.interval.support_work_j;
                visited.push(report.end_time_s);
                if let Some(jump) = report.jump {
                    jumps.push(jump.time_s);
                    // Consistent tetra mass gives f_f=(rho V/4)*a=.25*a,
                    // M_ff=.1, hence the exact acceleration jump is 2.5*a.
                    assert!((jump.acceleration_change_l2_m_s2-2.5*pulse).abs() < 1e-10);
                }
            }
            assert_eq!(jumps,vec![0.125,0.25]); assert!(visited.contains(&0.0625));
            assert!((run.run().state().total_energy_j-initial_energy-work).abs() < 1e-11);
            // Original equation: q''+100q=100*rate*t+2.5*a(t).
            // This independent continuum-in-time response includes both pulse edges.
            let t = 0.5_f64;
            let expected = rate*(t-(10.*t).sin()/10.)
                + 2.5*pulse/100.*((1.-(10.*(t-0.125)).cos())-(1.-(10.*(t-0.25)).cos()));
            errors.push((run.run().state().displacement[0]-expected).abs());
            let full = run.run().nodal_motion(cx).unwrap();
            assert_eq!(full[0][0][0],rate*t);
            assert_eq!(full[1][0][0],rate);
        }
        assert!(errors[1] < errors[0]/3.5 && errors[2] < errors[1]/3.5, "{errors:?}");
        assert!(errors[2] < 2e-6, "{errors:?}");
    });
}

#[test]
fn g4_support_discontinuities_and_cancellation_leave_the_accepted_event_run_unchanged() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.,60.,0.,[75;32]).unwrap();
        let p = problem(&material,&FIXED);
        let g = TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let frames = FrameTree::new();
        let motion = table(DISPLACEMENT_DIMS,&[0.,0.25,0.5],&[0.,0.001,0.001],TableInterpolation::StepLeft);
        let rows = [SupportMotion3 { nodes:&ALL,components:[true;3],displacement:&motion }];
        let schedule = EventLoadHistory3::new_with_supports(history(&frames),&rows,12,cx).unwrap();
        let mut run = schedule.begin_transient_with_supports(&g,&p,None,0.,&ZERO,&ZERO,&ZERO,
            Default::default(),cx).unwrap();
        let before = run.run().state().clone();
        let full = run.run().nodal_motion(cx).unwrap();
        assert!(run.advance_to(0.5,1.,cx).is_err());
        with_cx(true, |cancelled| { assert!(run.advance_to(0.125,1.,cancelled).is_err()); });
        assert_eq!(run.run().state(),&before);
        assert_eq!(run.run().nodal_motion(cx).unwrap(),full);
        let mut fork = run.clone();
        assert_eq!(run.advance_to(0.125,1.,cx).unwrap(),fork.advance_to(0.125,1.,cx).unwrap());
        assert!(run.advance_to(0.5,1.,cx).is_err());
        assert_eq!(run.run().state(),fork.run().state());
        assert_eq!(run.run().nodal_motion(cx).unwrap(),fork.run().nodal_motion(cx).unwrap());

        // A continuous vector in a discontinuously tilted frame still jumps in
        // the physical world. The support path must use the same exact left
        // frame as the force path, rather than sampling both at the right side.
        let mut frames = FrameTree::new();
        frames.add(Frame { id:FrameId(1),name:"held support frame".into(),parent:WORLD,
            motion:FrameMotion::Tilt { axis:[0.,0.,1.],center:Vec3::new(0.,0.,0.),
                angle:fs_scenario::TimeSignal::Table { times:vec![0.,0.25],
                    values:vec![0.,core::f64::consts::FRAC_PI_2],dims:Dims::NONE,
                    interp:fs_scenario::signal::Interp::Hold } } });
        let motion = source(DISPLACEMENT_DIMS,FrameId(1),[0.001,0.,0.],
            [0.001,0.,0.],OutsideDomainPolicy::Clamp);
        let rows = [SupportMotion3 { nodes:&ALL,components:[true;3],displacement:&motion }];
        let schedule = EventLoadHistory3::new_with_supports(history(&frames),&rows,12,cx).unwrap();
        let mut run = schedule.begin_transient_with_supports(&g,&p,None,0.,&[[0.001,0.,0.];4],
            &ZERO,&ZERO,Default::default(),cx).unwrap();
        let full = run.run().nodal_motion(cx).unwrap();
        assert!(run.advance_to(0.5,1.,cx).is_err());
        assert_eq!(run.run().state().time_s,0.);
        assert_eq!(run.run().nodal_motion(cx).unwrap(),full);
    });
}

#[test]
fn g3_adaptive_motion_error_includes_prescribed_nodes_and_retains_original_half_steps() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.,60.,0.,[75;32]).unwrap();
        let p = problem(&material,&FIXED);
        let g = TetLoadGeometry3::new(&NODES,&TETS,p.budget,4,cx).unwrap();
        let mut frames = FrameTree::new();
        frames.add(Frame { id:FrameId(1),name:"rotating prescribed translation".into(),parent:WORLD,
            motion:FrameMotion::Rotating { axis:[0.,0.,1.], center:Vec3::new(0.,0.,0.),
                rate:QtyAny::new(2.,Dims([0,0,-1,0,0,0])) } });
        let amplitude = 1.0/1024.0;
        let motion = source(DISPLACEMENT_DIMS,FrameId(1),[amplitude,0.,0.],
            [amplitude,0.,0.],OutsideDomainPolicy::Clamp);
        let rows = [SupportMotion3 { nodes:&ALL,components:[true;3],displacement:&motion }];
        let schedule = EventLoadHistory3::new_with_supports(history(&frames),&rows,8,cx).unwrap();
        let mut run = schedule.begin_transient_with_supports(&g,&p,None,0.,
            &[[amplitude,0.,0.];4],&[[0.,2.*amplitude,0.];4],&[[-4.*amplitude,0.,0.];4],
            Default::default(),cx).unwrap();
        assert!(run.run().state().displacement.is_empty());
        let before = run.run().state().clone();
        let full = run.run().nodal_motion(cx).unwrap();
        let strict = TransientTimeControl3 { absolute_displacement_m:1e-12,
            absolute_velocity_m_s:1e-12,absolute_acceleration_m_s2:1e-12,
            maximum_step_s:0.125,maximum_attempts:1,relative:0.,..Default::default() };
        assert!(matches!(run.advance_adaptive_to(0.125,0.125,strict,cx),Err(AdaptiveHistoryError3::Rejected { .. })));
        assert_eq!(run.run().state(),&before); assert_eq!(run.run().nodal_motion(cx).unwrap(),full);
        let mut coarse = run.clone();
        coarse.advance_to(0.125,0.125,cx).unwrap();
        let mut manual = run.clone();
        let first = manual.advance_to(0.0625,0.0625,cx).unwrap();
        let second = manual.advance_to(0.125,0.0625,cx).unwrap();
        let coarse_motion = coarse.run().nodal_motion(cx).unwrap();
        let fine_motion = manual.run().nodal_motion(cx).unwrap();
        let gap = (fine_motion[2][0][1]-coarse_motion[2][0][1]).abs();
        let exact_acceleration = -4.*amplitude*0.25_f64.sin();
        let actual_error = (fine_motion[2][0][1]-exact_acceleration).abs();
        // The cubic term of sin(2t) makes Newmark's prescribed acceleration
        // fail the free-state second-order Richardson assumption. In this
        // analytic example /3 would UNDERestimate the fine derivative error.
        assert!(actual_error > gap/3. && actual_error < gap);
        let consistency = TransientTimeControl3 { absolute_displacement_m:0.01,
            absolute_velocity_m_s:0.01,absolute_acceleration_m_s2:gap/2.,..strict };
        assert!(matches!(run.advance_adaptive_to(0.125,0.125,consistency,cx),Err(AdaptiveHistoryError3::Rejected { .. })));
        assert_eq!(run.run().nodal_motion(cx).unwrap(),full);
        let loose = TransientTimeControl3 { absolute_displacement_m:0.01,
            absolute_velocity_m_s:0.01,absolute_acceleration_m_s2:0.01,..strict };
        let accepted = run.advance_adaptive_to(0.125,0.125,loose,cx).unwrap();
        assert!(accepted.error_ratio > 0. && accepted.error_ratio < 1.);
        assert!((accepted.error_ratio*loose.absolute_acceleration_m_s2-gap).abs() < 1e-12);
        assert_eq!(accepted.first,first); assert_eq!(accepted.second,second);
        assert_eq!(run.run().nodal_motion(cx).unwrap(),manual.run().nodal_motion(cx).unwrap());
    });
}

#[test]
fn g1_uniform_translation_keeps_support_velocity_and_accounts_for_damping_work() {
    with_cx(false, |cx| {
        // Total mass = rho*V = 1 kg. A rigid translation has zero elastic
        // energy; prescribed supports supply the eta*M*v damping force.
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let p = problem(&material, &FIXED);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let h = history(&frames);
        let rate = 1.0 / 1024.0;
        let motion = source(
            DISPLACEMENT_DIMS,
            WORLD,
            [0.0; 3],
            [rate, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let rows = [SupportMotion3 {
            nodes: &ALL,
            components: [true; 3],
            displacement: &motion,
        }];
        let config = TransientConfig3 {
            mass_damping_per_s: 0.25,
            ..Default::default()
        };
        let mut run = h
            .begin_transient_with_support_history(
                &g,
                &p,
                None,
                &rows,
                0.0,
                &ZERO,
                &[[rate, 0.0, 0.0]; 4],
                &ZERO,
                config,
                cx,
            )
            .unwrap();
        let mut work = 0.0;
        for _ in 0..8 {
            let step = run
                .advance_history_with_supports(&h, &rows, 1.0 / 32.0, cx)
                .unwrap();
            assert!(step.balance_error_j.abs() < 1e-12);
            assert!((step.support_work_j - step.damping_loss_j).abs() < 1e-12);
            work += step.support_work_j;
        }
        let state = run.nodal_motion(cx).unwrap();
        for n in 0..4 {
            assert!((state[0][n][0] - rate * 0.25).abs() < 1e-13);
            assert!((state[1][n][0] - rate).abs() < 1e-13);
            assert!(state[2][n][0].abs() < 1e-12);
        }
        assert!((run.state().kinetic_energy_j - 0.5 * rate * rate).abs() < 1e-12);
        assert!((work - 0.25 * rate * rate * 0.25).abs() < 1e-12);
    });
}

#[test]
fn g3_typed_force_and_support_samples_match_manual_coupled_endpoints() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let fixed: Vec<_> = FIXED.into_iter().filter(|d| *d != 3).collect();
        let p = problem(&material, &fixed);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let h = TetLoadHistory3 {
            nodal_forces_n: &[(3, 0.01)],
            ..history(&frames)
        };
        let rate = 1.0 / 1024.0;
        let motion = source(
            DISPLACEMENT_DIMS,
            WORLD,
            [0.0; 3],
            [rate, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let rows = [SupportMotion3 {
            nodes: &[0],
            components: [true, false, false],
            displacement: &motion,
        }];
        let mut velocity = ZERO;
        velocity[0][0] = rate;
        let load = g
            .assemble(&[], VolumeLoad3::None, &[(3, 0.01)], cx)
            .unwrap();
        let config = TransientConfig3::default();
        let mut scheduled = h
            .begin_transient_with_support_history(
                &g, &p, None, &rows, 0.0, &ZERO, &velocity, &ZERO, config, cx,
            )
            .unwrap();
        let mut direct = load
            .begin_transient_with_supports(&p, None, 0.0, &ZERO, &velocity, &ZERO, config, cx)
            .unwrap();
        for _ in 0..16 {
            let time = scheduled.state().time_s + 1.0 / 128.0;
            let a = scheduled
                .advance_history_with_supports(&h, &rows, 1.0 / 128.0, cx)
                .unwrap();
            let b = direct
                .advance_with_supports(&load, &[(0, rate * time)], 1.0 / 128.0, cx)
                .unwrap();
            assert_eq!(a, b);
            assert_eq!(scheduled.state(), direct.state());
            assert_eq!(
                scheduled.nodal_motion(cx).unwrap(),
                direct.nodal_motion(cx).unwrap()
            );
        }
        assert!(scheduled.state().velocity[0].abs() > 1e-6);
    });
}

#[test]
fn g3_rotating_support_keeps_physical_frame_time_beyond_clamped_value_table() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let p = problem(&material, &FIXED);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let mut frames = FrameTree::new();
        frames.add(Frame {
            id: FrameId(1),
            name: "moving-support-vector".into(),
            parent: WORLD,
            motion: FrameMotion::Rotating {
                axis: [0.0, 0.0, 1.0],
                center: Vec3::new(1e6, -2e6, 0.0),
                rate: QtyAny::new(2.0, Dims([0, 0, -1, 0, 0, 0])),
            },
        });
        let h = history(&frames);
        let amplitude = 1.0 / 1024.0;
        let motion = source(
            DISPLACEMENT_DIMS,
            FrameId(1),
            [0.0; 3],
            [amplitude, 0.0, 0.0],
            OutsideDomainPolicy::Clamp,
        );
        let rows = [SupportMotion3 {
            nodes: &ALL,
            components: [true; 3],
            displacement: &motion,
        }];
        // Use the exact admitted endpoint for q0; the history adapter forbids
        // mixing an independently rounded initialization with the sampled BC.
        let start = 1.25;
        let sample = h.at_with_supports(&g, &p, &rows, start, 0.01, cx).unwrap();
        let mut displacement = ZERO;
        for &(d, u) in sample.prescribed_dofs_m() {
            displacement[d / 3][d % 3] = u;
        }
        let angle = 2.0_f64 * start;
        let velocity = [[
            -2.0 * amplitude * angle.sin(),
            2.0 * amplitude * angle.cos(),
            0.0,
        ]; 4];
        let acceleration = displacement.map(|row| row.map(|u| -4.0 * u));
        let mut run = h
            .begin_transient_with_support_history(
                &g,
                &p,
                None,
                &rows,
                start,
                &displacement,
                &velocity,
                &acceleration,
                Default::default(),
                cx,
            )
            .unwrap();
        run.advance_history_with_supports(&h, &rows, 1.0 / 128.0, cx)
            .unwrap();
        let endpoint = run.nodal_motion(cx).unwrap();
        let angle = 2.0 * run.state().time_s;
        for q in &endpoint[0] {
            assert!((q[0] - amplitude * angle.cos()).abs() < 1e-12);
            assert!((q[1] - amplitude * angle.sin()).abs() < 1e-12);
        }
        assert!((endpoint[0][0][0] - amplitude * 2.0_f64.cos()).abs() > 1e-4);
    });
}

#[test]
fn g4_invalid_histories_initial_mismatch_and_cancellation_publish_no_partial_motion() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(6.0, 60.0, 0.0, [75; 32]).unwrap();
        let p = problem(&material, &FIXED);
        let g = TetLoadGeometry3::new(&NODES, &TETS, p.budget, 4, cx).unwrap();
        let frames = FrameTree::new();
        let h = history(&frames);
        let motion = source(
            DISPLACEMENT_DIMS,
            WORLD,
            [0.0; 3],
            [0.001, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let wrong = source(
            TRACTION_DIMS,
            WORLD,
            [0.0; 3],
            [0.001, 0.0, 0.0],
            OutsideDomainPolicy::Refuse,
        );
        let valid = SupportMotion3 {
            nodes: &ALL,
            components: [true; 3],
            displacement: &motion,
        };
        assert!(
            h.begin_transient_with_support_history(
                &g,
                &p,
                None,
                &[valid],
                0.0,
                &[[0.001, 0.0, 0.0]; 4],
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .is_err()
        );
        let mut run = h
            .begin_transient_with_support_history(
                &g,
                &p,
                None,
                &[valid],
                0.0,
                &ZERO,
                &ZERO,
                &ZERO,
                Default::default(),
                cx,
            )
            .unwrap();
        let before = run.state().clone();
        let full_before = run.nodal_motion(cx).unwrap();
        for bad in [
            vec![SupportMotion3 {
                displacement: &wrong,
                ..valid
            }],
            vec![valid, valid],
            vec![SupportMotion3 {
                nodes: &[4],
                ..valid
            }],
        ] {
            assert!(
                run.advance_history_with_supports(&h, &bad, 0.01, cx)
                    .is_err()
            );
        }
        assert!(
            run.advance_history_with_supports(&h, &[valid], 2.0, cx)
                .is_err()
        );
        for step in [0.0, -0.01, f64::NAN, f64::INFINITY] {
            assert!(
                run.advance_history_with_supports(&h, &[valid], step, cx)
                    .is_err()
            );
        }
        with_cx(true, |cancelled| {
            assert!(
                run.advance_history_with_supports(&h, &[valid], 0.01, cancelled)
                    .is_err()
            );
        });
        assert_eq!(run.state(), &before);
        assert_eq!(run.nodal_motion(cx).unwrap(), full_before);
        let mut fork = run.clone();
        assert_eq!(
            run.advance_history_with_supports(&h, &[valid], 0.01, cx)
                .unwrap(),
            fork.advance_history_with_supports(&h, &[valid], 0.01, cx)
                .unwrap()
        );
        assert_eq!(
            run.nodal_motion(cx).unwrap(),
            fork.nodal_motion(cx).unwrap()
        );
    });
}
