//! Local time accuracy, actual half-step energy and atomic adaptive publication.
mod support;

use fs_conduction::bc::history::{ThermalBoundaryHistory};
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryMethod,
    ThermalHistoryProblem, ThermalHistoryRun};
use fs_conduction::bc::history::march::events::adaptive::AdaptiveHistoryConfig;
use fs_conduction::transient::VolumetricHeatCapacity;
use fs_conduction::transient::backward_euler::{NonlinearStepConfig, StepConfig};
use fs_conduction::{ConductionError, ConductionMesh, ConductivityModel, ConductivityTable, LinearConfig,
    ScalarField, ThermalBc, ThermalBoundaryBuilder, VOLUMETRIC_SOURCE_DIMS};
use fs_exec::Cx;
use fs_qty::QtyAny;
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{DEFAULT_MAXIMUM_CYCLES_PER_STEP, Interp, LinearChirp, TimeSignal};
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    ConductionMesh::new(TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]).unwrap()
}
fn controls() -> StepConfig {
    StepConfig { linear: LinearConfig { tolerance: 1e-12, ..LinearConfig::default() },
        energy_tolerance_j: 1e-7 }
}
fn adaptive(tolerance: f64) -> AdaptiveHistoryConfig {
    AdaptiveHistoryConfig { absolute_temperature_k: tolerance, relative_temperature: 0.0,
        min_step_s: 1e-7, max_step_s: 0.25, max_attempts: 24 }
}
fn cold<'m>(cx: &Cx<'_>, mesh: &'m ConductionMesh, time: f64) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(mesh)
        .region("cold", |f| (f.centroid.iter().sum::<f64>() - 1.0).abs() < 1e-12,
            ThermalBc::dirichlet(300.0).unwrap()).unwrap().adiabatic_remainder(),
        &Environment::earth_lab(), Vec::new(), time).unwrap()
}
fn problem<'h, 'm>(history: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    source: HeatSourceHistory) -> ThermalHistoryProblem<'h, 'm>
{
    ThermalHistoryProblem { boundary: history, material, source, element_materials: None, interfaces: None }
}
fn cooling<'h, 'm>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>, material: &'h ConductivityModel,
    budget: usize, method: ThermalHistoryMethod) -> ThermalHistoryRun<'h, 'm>
{
    ThermalHistoryRun::new_with_method(cx, problem(history, material,
        HeatSourceHistory::Fixed(ScalarField::Uniform(0.0))),
        VolumetricHeatCapacity::declared(15.0).unwrap(), &[400., 300., 300., 300.], budget, method).unwrap()
}

#[test]
fn g1_adaptive_cooling_publishes_the_original_two_half_steps_not_an_extrapolation() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    with_cx(|cx| {
        let history = cold(cx, &mesh, 0.0);
        for method in [ThermalHistoryMethod::Linear, ThermalHistoryMethod::Nonlinear(NonlinearStepConfig::default())] {
            let mut run = cooling(cx, &history, &material, 2, method);
            let step = run.advance_adaptive_step(cx, 0.25, 0.25, controls(), adaptive(0.02)).unwrap();
            assert!(step.rejected_attempts > 0);
            assert_eq!(step.endpoint_solves, 3 * (step.rejected_attempts + 1));
            assert!(step.normalized_error <= 1.0);
            // Independent one-free-node equation: C00=15/24, K00=20/2,
            // so d(T0-300)/dt=-16*(T0-300). The other three nodes are fixed.
            let first = 300.0 + 100.0 / (1.0 + 16.0 * step.midpoint_time_s);
            let second = 300.0 + (first - 300.0) / (1.0 + 16.0 * (step.end_time_s - step.midpoint_time_s));
            assert!((step.first_half.solution().temperature[0] - first).abs() < 1e-9);
            assert!((step.second_half.solution().temperature[0] - second).abs() < 1e-9);
            for half in [&step.first_half, &step.second_half] {
                assert!(half.solution().energy_residual_j.abs() <= controls().energy_tolerance_j);
                assert_eq!(&half.solution().temperature[1..], &[300.; 3]);
            }
            assert_eq!(run.temperature(), step.second_half.solution().temperature.as_slice());
            assert_eq!(run.time_s(), step.end_time_s);
            assert_eq!(run.accepted_steps(), 2);
        }
    });
}

#[test]
fn g1_tighter_local_tolerance_improves_the_independent_transient_solution() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    let curved = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 15.0), (500.0, 40.0)]).unwrap());
    with_cx(|cx| {
        let history = cold(cx, &mesh, 0.0);
        for (material, method, nonlinear) in [(&material, ThermalHistoryMethod::Linear, false),
            (&curved, ThermalHistoryMethod::Nonlinear(NonlinearStepConfig::default()), true)] {
        let mut errors = Vec::new();
        let mut counts = Vec::new();
        for tolerance in [0.2, 0.002] {
            let mut run = cooling(cx, &history, material, 10_000, method);
            let mut proposal = 0.25;
            while run.time_s() < 0.25 {
                let step = run.advance_adaptive_step(cx, 0.25, proposal, controls(), adaptive(tolerance)).unwrap();
                assert!(step.normalized_error <= 1.0);
                proposal = step.suggested_step_s;
            }
            // With k=0.1*T_mean-10, x=T0-300 obeys x'=-16*x-0.02*x^2.
            let decay = (-4.0_f64).exp();
            let exact = 300.0 + 100.0 * decay / if nonlinear { 1.0 + 0.125 * (1.0 - decay) } else { 1.0 };
            errors.push((run.temperature()[0] - exact).abs());
            counts.push(run.accepted_steps());
        }
        assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
        assert!(counts[1] > counts[0]);
        }
    });
}

#[test]
fn g1_adaptation_does_not_drop_the_interval_ending_at_a_heater_switch() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.0).unwrap();
        let source = HeatSourceHistory::Uniform(TimeSignal::Table { times: vec![0.0, 0.125],
            values: vec![120.0, 0.0], dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold });
        let mut run = ThermalHistoryRun::new(cx, problem(&history, &material, source),
            VolumetricHeatCapacity::declared(30.0).unwrap(), &[300.; 4], 4).unwrap();
        let step = run.advance_adaptive_step(cx, 0.25, 0.25, controls(), adaptive(1e-6)).unwrap();
        assert_eq!(step.end_time_s, 0.125);
        assert_eq!(step.rejected_attempts, 0);
        assert!(run.temperature().iter().all(|t| (t - 300.5).abs() < 1e-9));
        let work = step.first_half.solution().source_w * step.midpoint_time_s
            + step.second_half.solution().source_w * (step.end_time_s - step.midpoint_time_s);
        assert!((work - 2.5).abs() < 1e-10); // 120 * (tet volume 1/6) * 0.125.
        run.advance_adaptive_step(cx, 0.25, 0.25, controls(), adaptive(1e-6)).unwrap();
        assert!(run.temperature().iter().all(|t| (t - 300.5).abs() < 1e-9));
    });
}

#[test]
fn g1_frequency_resolved_history_prevents_whole_cycle_heat_aliasing() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    with_cx(|cx| {
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.0).unwrap();
        let signal = TimeSignal::Chirp(LinearChirp::new(0.0, 1.0, 8.0, 8.0, 0.0,
            QtyAny::new(120.0, VOLUMETRIC_SOURCE_DIMS),
            QtyAny::new(120.0, VOLUMETRIC_SOURCE_DIMS)).unwrap());
        let make_run = |budget| ThermalHistoryRun::new(cx,
            problem(&history, &material, HeatSourceHistory::Uniform(signal.clone())),
            VolumetricHeatCapacity::declared(30.0).unwrap(), &[300.0; 4], budget).unwrap();

        // Two complete source cycles are invisible to both the full interval
        // and its midpoint: all three endpoint loads equal 240 W/m^3. The
        // deliberately unaligned API therefore gives agreeing WRONG histories.
        let mut coarse = make_run(1);
        coarse.advance_to(cx, 0.25, controls()).unwrap();
        let mut fine = make_run(2);
        fine.advance_to(cx, 0.125, controls()).unwrap();
        fine.advance_to(cx, 0.25, controls()).unwrap();
        for (&a, &b) in coarse.temperature().iter().zip(fine.temperature()) {
            assert!((a - b).abs() < 1e-9);
            assert!((a - 302.0).abs() < 1e-9);
        }

        // Source-aware event selection resolves the same physical history
        // before the local temperature-error comparison. A loose tolerance
        // isolates that frequency guard from ordinary adaptive refinement.
        let mut resolved = make_run(128);
        let mut work = 0.0;
        let mut pairs = 0;
        while resolved.time_s() < 0.25 {
            let step = resolved.advance_adaptive_step(
                cx, 0.25, 0.25, controls(), adaptive(1.0)).unwrap();
            assert!(8.0 * (step.end_time_s - step.start_time_s)
                <= DEFAULT_MAXIMUM_CYCLES_PER_STEP);
            assert_eq!(step.rejected_attempts, 0);
            work += step.first_half.solution().source_w
                * (step.midpoint_time_s - step.start_time_s)
                + step.second_half.solution().source_w
                * (step.end_time_s - step.midpoint_time_s);
            pairs += 1;
        }
        assert!(pairs >= 32);
        // Integral q dt = 120*0.25; the cosine contributes zero. The original
        // tetrahedron has volume 1/6, so supplied heat is 5 J and delta T is 1 K.
        assert!((work - 5.0).abs() < 1e-8);
        assert!(resolved.temperature().iter().all(|&t| (t - 301.0).abs() < 1e-8));
        assert_eq!(resolved.accepted_steps(), 2 * pairs);
    });
}

#[test]
fn g0_rejection_cancellation_and_budget_exhaustion_keep_the_entire_accepted_state() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    let history = with_cx(|cx| cold(cx, &mesh, 0.0));
    let mut run = with_cx(|cx| cooling(cx, &history, &material, 2, ThermalHistoryMethod::Linear));
    with_cancelled_cx(|cx| assert!(matches!(run.advance_adaptive_step(cx, 0.25, 0.25,
        controls(), adaptive(0.02)), Err(ConductionError::Cancelled { .. }))));
    with_cx(|cx| {
        let impossible = AdaptiveHistoryConfig { max_attempts: 1, ..adaptive(1e-20) };
        assert!(run.advance_adaptive_step(cx, 0.25, 0.25, controls(), impossible).is_err());
        assert_eq!(run.temperature(), &[400., 300., 300., 300.]);
        assert_eq!(run.time_s(), 0.0);
        assert_eq!(run.accepted_steps(), 0);
        let first = run.advance_adaptive_step(cx, 0.25, 0.25, controls(), adaptive(0.02)).unwrap();
        let accepted = run.temperature().to_vec();
        assert!(run.advance_adaptive_step(cx, 0.25, first.suggested_step_s, controls(), adaptive(0.02)).is_err());
        assert_eq!(run.temperature(), accepted.as_slice());
        assert_eq!(run.time_s(), first.end_time_s);
        assert_eq!(run.accepted_steps(), 2);
        run.extend_step_budget(4).unwrap();
        run.advance_adaptive_step(cx, 0.25, first.suggested_step_s, controls(), adaptive(0.02)).unwrap();
        let mut straight = cooling(cx, &history, &material, 4, ThermalHistoryMethod::Linear);
        let first = straight.advance_adaptive_step(cx, 0.25, 0.25, controls(), adaptive(0.02)).unwrap();
        straight.advance_adaptive_step(cx, 0.25, first.suggested_step_s, controls(), adaptive(0.02)).unwrap();
        assert_eq!(run.temperature(), straight.temperature());
        assert_eq!(run.time_s(), straight.time_s());
        assert_eq!(run.accepted_steps(), 4);
    });
}

#[test]
fn g0_no_midpoint_or_invalid_time_policy_never_silently_accepts_an_unchecked_step() {
    let mesh = mesh();
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    with_cx(|cx| {
        let history = cold(cx, &mesh, 1.0);
        let mut run = cooling(cx, &history, &material, 2, ThermalHistoryMethod::Linear);
        assert!(run.advance_adaptive_step(cx, 1.0_f64.next_up(), 0.25, controls(), adaptive(0.02)).is_err());
        for config in [AdaptiveHistoryConfig { max_attempts: 0, ..adaptive(0.02) },
            AdaptiveHistoryConfig { absolute_temperature_k: f64::NAN, ..adaptive(0.02) },
            AdaptiveHistoryConfig { relative_temperature: 1.0, ..adaptive(0.02) }] {
            assert!(run.advance_adaptive_step(cx, 1.25, 0.25, controls(), config).is_err());
        }
        assert_eq!(run.time_s(), 1.0);
        assert_eq!(run.accepted_steps(), 0);
        assert_eq!(run.temperature(), &[400., 300., 300., 300.]);
    });
}
