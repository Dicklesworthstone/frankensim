//! Source-bound reduced motion and physical transfer, using an independent
//! one-free-DOF tetrahedron oracle. Manufactured data, not measured acoustics.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_modal::transfer::ModalTransferControls;
use fs_solid::linear3::{
    TetElasticMaterial, TetMaterialField, TetStaticSolveConfig, TetThermalStrainField,
    TetThermalStrainState,
};
use fs_solid::reduce::thermal::motion::ThermalValveModalMotion;
use fs_solid::reduce::thermal::{
    ThermalPressureLoadedValve, ThermalValveControls, reduce_pressure_loaded_valve,
};
use fs_solid::reduce::{OrificeSpec, ValveCardRequest, ValveElementOrder};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const FRAME: ContentHash = ContentHash([41; 32]);
const RATE: u32 = 1000;

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(gate, arena, StreamKey { seed: 1, kernel_id: 712, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic))
    })
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new_clock_free(), f)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}

fn valve(young: f64, density: f64, expansion: f64, free: usize, frame: ContentHash,
    cx: &Cx<'_>) -> ThermalPressureLoadedValve
{
    let material = TetElasticMaterial::try_new(density, young, 0.25, [7; 32]).unwrap();
    let thermal = TetThermalStrainState::try_new(
        300.0 + 1000.0 * expansion, 300.0, [expansion, expansion, expansion, 0.0, 0.0, 0.0],
        [7; 32], ContentHash([9; 32]), 0.1,
    ).unwrap();
    let fixed = (0..12).filter(|&dof| dof != free).collect::<Vec<_>>();
    let mut pressure = [[0.0; 3]; 4];
    pressure[free / 3][free % 3] = 0.03;
    reduce_pressure_loaded_valve(
        &ValveCardRequest {
            nodes_m: &NODES,
            tetrahedra: &TETS,
            element_order: ValveElementOrder::Tet4,
            fixed_dofs: &fixed,
            material: TetMaterialField::Uniform(&material),
            orifice: OrificeSpec {
                face_nodes: &[1, 2, 3], opening_axis: [0.0, 0.0, 1.0],
                width_axis: [0.0, 1.0, 0.0], opposing_plane_offset_m: -0.1,
            },
            window_hz: (0.1, 100.0), retained_compliance_floor: 0.999,
            loss_factor: 0.0, source_id: "manufactured/one-tet-thermal-motion",
        },
        TetThermalStrainField::Uniform(&thermal), &[], &pressure,
        ThermalValveControls { solve: TetStaticSolveConfig::default(), mesh_frame: frame,
            maximum_absolute_mandel_strain: 0.1, maximum_nodal_displacement_m: 0.02 },
        cx,
    ).unwrap()
}

fn observations(valve: &ThermalPressureLoadedValve) -> Vec<C64> {
    // Explicit observer p = 2 Pa s/m * physical velocity of the free node.
    valve.valve().source_modes().iter().map(|mode| C64::new(2.0 * mode.phi[0], 0.0)).collect()
}
fn controls() -> ModalTransferControls {
    ModalTransferControls {
        maximum_dofs: 12, maximum_modes: 12, maximum_factor_bytes: 1 << 20,
        maximum_mass_orthogonality_error: 1e-9, maximum_relative_eigen_residual: 1e-9,
        maximum_displacement_residual_m_sqrt_kg: 1e-9,
        maximum_kinetic_residual_j: 1e-12, maximum_elastic_residual_j: 1e-12,
        maximum_abs_parameter_work_j: 1.0, energy_tolerance_j: 1e-10,
    }
}
fn motion(cx: &Cx<'_>) -> ThermalValveModalMotion {
    let source = valve(1000.0, 6.0, 0.0, 11, FRAME, cx);
    assert_eq!(source.valve().source_modes().len(), 1);
    let observation = observations(&source);
    ThermalValveModalMotion::new(source, RATE,
        &[ModalAcousticState { displacement_m_sqrt_kg: 3e-4, velocity_m_sqrt_kg_per_s: 4e-3 }],
        &observation, ModalAcousticTimeBudget::audible_reference(), cx).unwrap()
}

#[test]
fn pressure_port_matches_independent_tetrahedron_oscillator_and_replays() {
    with_cx(|cx| {
        let mut run = motion(cx);
        let mut replay = run.clone();
        let source = run.valve().equilibrium().pencil();
        // Unit tet volume 1/6; consistent corner mass rho*V/10 = 0.1.
        // C_zz = 1.2 E and grad N_3 = e_z, so K_zz = E/5 = 200.
        close(source.mass.get(0, 0), 0.1, 1e-14);
        close(source.stiffness.get(0, 0), 200.0, 1e-11);
        let initial = run.nodal_motion(cx).unwrap();
        let mut u = initial.displacement_m[3][2];
        let mut v = initial.velocity_m_s[3][2];
        let omega = (2000.0_f64).sqrt();
        let h = 1.0 / f64::from(RATE);
        for step in 1..=16 {
            let pressure = if step <= 8 { 0.4 } else { 0.0 };
            let equilibrium = 0.03 * pressure / 200.0;
            let (sin, cos) = (omega * h).sin_cos();
            let next_u = equilibrium + (u - equilibrium) * cos + v * sin / omega;
            let next_v = -(u - equilibrium) * omega * sin + v * cos;
            let frame = run.step_pressure(pressure, cx).unwrap();
            let repeat = replay.step_pressure(pressure, cx).unwrap();
            let nodes = run.nodal_motion(cx).unwrap();
            close(nodes.displacement_m[3][2], next_u, 2e-13);
            close(nodes.velocity_m_s[3][2], next_v, 2e-12);
            close(frame.swept_volume_change_m3, 0.03 * (next_u - u), 1e-15);
            close(frame.modal.input_work_j, pressure * frame.swept_volume_change_m3, 1e-15);
            close(frame.modal.total_modal_energy_j, 0.05 * next_v.powi(2) + 100.0 * next_u.powi(2), 1e-13);
            close(frame.modal.observer_pressure_pa, 2.0 * next_v, 2e-12);
            assert_eq!(frame.step, step);
            assert_eq!(frame.time_s.to_bits(), (step as f64 * h).to_bits());
            assert_eq!(frame.modal, repeat.modal);
            assert_eq!(run.model().states(), replay.model().states());
            (u, v) = (next_u, next_v);
        }
    });
}

#[test]
fn thermal_source_update_preserves_physical_motion_clock_and_explicit_storage_change() {
    with_cx(|cx| {
        let mut old = motion(cx);
        old.step_pressure(0.4, cx).unwrap();
        let before = old.nodal_motion(cx).unwrap();
        let hot = valve(1600.0, 6.0, 0.001, 11, FRAME, cx);
        close(hot.equilibrium().displacement_m()[3][2], 5e-3 / 3.0, 1e-12);
        close(hot.equilibrium().pencil().stiffness.get(0, 0), 320.0, 1e-9);
        let observation = observations(&hot);
        let (mut updated, audit) = old.rebind_candidate(hot, &observation, controls(), cx).unwrap();
        let after = updated.nodal_motion(cx).unwrap();
        for (a, b) in before.displacement_m.iter().flatten().zip(after.displacement_m.iter().flatten()) {
            close(*a, *b, 1e-12);
        }
        for (a, b) in before.velocity_m_s.iter().flatten().zip(after.velocity_m_s.iter().flatten()) {
            close(*a, *b, 1e-12);
        }
        assert_eq!(updated.accepted_steps(), old.accepted_steps());
        assert_eq!(updated.time_s().to_bits(), old.time_s().to_bits());
        let u = before.displacement_m[3][2];
        let expected_work = 160.0 * (u - 5e-3 / 3.0).powi(2) - 100.0 * u.powi(2);
        close(audit.parameter_work_j, expected_work, 1e-12);
        assert!(audit.parameter_work_j.abs() > 1e-8);
        assert!(audit.kinetic_residual_j < 1e-20 && audit.elastic_residual_j < 1e-20);
        assert!(audit.projection_energy_residual_j.abs() < 1e-12);
        let cold_frame = old.step_pressure(0.4, cx).unwrap();
        let hot_frame = updated.step_pressure(0.4, cx).unwrap();
        assert_eq!(cold_frame.time_s, hot_frame.time_s);
        assert_ne!(cold_frame.modal.observer_pressure_pa, hot_frame.modal.observer_pressure_pa);
        assert_ne!(old.nodal_motion(cx).unwrap(), updated.nodal_motion(cx).unwrap());
    });
}

#[test]
fn inadmissible_source_observer_budget_and_cancellation_preserve_the_accepted_prefix() {
    let mut old = with_cx(motion);
    with_cx(|cx| {
        old.step_pressure(0.4, cx).unwrap();
        let state = old.model().states().to_vec();
        let identity = old.valve().identity();
        let time = old.time_s().to_bits();
        for (density, free, frame) in [(12.0, 11, FRAME), (6.0, 8, FRAME), (6.0, 11, ContentHash([42; 32]))] {
            let next = valve(1600.0, density, 0.0, free, frame, cx);
            let observation = observations(&next);
            assert!(old.rebind_candidate(next, &observation, controls(), cx).is_err());
        }
        let next = valve(1600.0, 6.0, 0.001, 11, FRAME, cx);
        let observation = observations(&next);
        let no_work = ModalTransferControls { maximum_abs_parameter_work_j: 0.0, ..controls() };
        assert!(old.rebind_candidate(next.clone(), &observation, no_work, cx).is_err());
        assert!(old.rebind_candidate(next.clone(), &[], controls(), cx).is_err());
        assert!(old.rebind_candidate(next, &[C64::new(f64::NAN, 0.0)], controls(), cx).is_err());
        assert!(old.step_pressure(f64::NAN, cx).is_err());
        assert_eq!(old.model().states(), state);
        assert_eq!(old.valve().identity(), identity);
        assert_eq!(old.time_s().to_bits(), time);
        assert_eq!(old.accepted_steps(), 1);
    });
    let state = old.model().states().to_vec();
    let next = with_cx(|cx| valve(1600.0, 6.0, 0.001, 11, FRAME, cx));
    let observation = observations(&next);
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert!(old.rebind_candidate(next, &observation, controls(), cx).is_err());
        assert!(old.step_pressure(0.4, cx).is_err());
    });
    assert_eq!(old.model().states(), state);
    assert_eq!(old.accepted_steps(), 1);
}
