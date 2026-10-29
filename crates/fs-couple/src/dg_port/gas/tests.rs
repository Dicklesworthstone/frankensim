use super::*;
use fs_duct::dg::{DgMesh, DgSidePort, MAX_DG_PROPAGATOR_DOFS};
use fs_duct::{Duct, Segment};
use fs_matdb::QueryPoint;
use fs_material::gas_field::PrimitiveGasState;
use fs_material::gas_species::input::ThermochemicalPhase;
use fs_material::thermochemical_gas::packs::from_packs;
use fs_math::c64::C64;
use fs_thermochem::{Composition, CompositionBasis};
use fs_vfit::mimo::MimoModel;
use fs_vfit::mimo::load::MimoLoadLimits;
use fs_vfit::model::PoleTerm;

#[path = "../../render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const DT: f64 = 1.0 / 96_000.0;
const INLET: PortId = PortId(7);
const OUTLET: PortId = PortId(3);
const SIDE: PortId = PortId(11);

fn field(temperature: f64, pressure: f64, length: f64) -> AxialGasField {
    let source = gas_support::packs("live-port-gas", 0.028, 3.5, 200.0);
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap();
    let model = from_packs(&[source.binding()], composition.clone(), None).unwrap();
    AxialGasField::new(
        vec![0.0, length],
        pressure,
        vec![
            PrimitiveGasState {
                temperature_k: temperature,
                composition
            };
            2
        ],
        model,
        QueryPoint::new(),
        ThermochemicalPhase::FrozenGas,
    )
    .unwrap()
}

fn duct(field: &AxialGasField, elements: usize, inlet: PortId, outlet: PortId) -> DgDuct {
    DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: field.length_m(),
            }],
        },
        field,
        DgMesh { elements, order: 2 },
        DgBoundary::Port(inlet),
        DgBoundary::Port(outlet),
    )
    .unwrap()
}

fn load(field: &AxialGasField) -> DgPortLoad {
    let main = duct(field, 2, INLET, OUTLET)
        .with_side_ports(&[DgSidePort {
            id: SIDE,
            downstream_element: 1,
        }])
        .unwrap();
    // The branch is its own explicitly supplied, warmer material owner.
    let branch = duct(
        &self::field(320.0, 101_325.0, 0.012),
        1,
        PortId(0),
        PortId(1),
    );
    let model = MimoModel {
        ports: 2,
        poles: vec![PoleTerm::Real {
            pole: -400.0,
            residue: 0.0,
        }],
        residues: vec![vec![
            C64::from_re(2e9),
            C64::from_re(5e8),
            C64::from_re(5e8),
            C64::from_re(1.8e9),
        ]],
        d: vec![2e7, 5e6, 5e6, 3e7],
        e: vec![0.0; 4],
        reciprocal: true,
    };
    let terminal = PreparedMimoLoad::new(
        &model,
        DT,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 129,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    DgPortLoad::new_with_terminal_branches(
        main,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        vec![SIDE, OUTLET],
        terminal,
        vec![DgBranch {
            port: SIDE,
            duct: branch,
        }],
    )
    .unwrap()
}

fn slit(port: PortId) -> DgTerminalChannelInput {
    DgTerminalChannelInput {
        port,
        kind: DgTerminalChannelKind::Slit {
            width_m: 0.002,
            height_m: 4e-5,
            length_m: 0.001,
            transverse_cells: 8,
            maximum_mach: 0.2,
        },
    }
}

fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 3e-9 * scale.max(1e-30), "{a:e} != {b:e}");
}

#[test]
fn g3_cold_gas_refresh_rebuilds_branch_and_terminal_law_like_native_assembly() {
    let before = field(300.0, 101_325.0, 0.05);
    let after = field(450.0, 101_325.0, 0.05);
    let mut refreshed = load(&before);
    let report = refreshed
        .rebind_gas_field_controlled(&after, DgGasRefresh::Lossless, 0.0, &mut || false)
        .unwrap();
    let mut assembled = load(&after);
    assert_eq!(report.parameter_work_j, 0.0);
    assert_eq!(refreshed.duct().gas_field(), Some(&after));
    assert_eq!(refreshed.state().step_count(), 0);
    near(
        refreshed.impedance_pa_s_m3(),
        assembled.impedance_pa_s_m3(),
        assembled.impedance_pa_s_m3(),
    );
    for sample in 0..20 {
        let flow = if sample < 9 { 1e-9 } else { -2e-10 };
        let a = refreshed.preview_step(flow).unwrap().accept();
        let b = assembled.preview_step(flow).unwrap().accept();
        let scale = b
            .port_outputs
            .iter()
            .map(|port| port.conjugate_rate.abs())
            .fold(1e-20, f64::max);
        for (a, b) in a.port_outputs.iter().zip(&b.port_outputs) {
            near(a.conjugate_rate, b.conjugate_rate, scale);
        }
        near(
            refreshed.stored_energy_j(),
            assembled.stored_energy_j(),
            assembled.stored_energy_j(),
        );
    }
}

#[test]
fn g1_warm_bore_refresh_preserves_branch_gap_and_causal_terminal_histories() {
    let before = field(300.0, 101_325.0, 0.05);
    let after = field(450.0, 101_325.0, 0.05);
    let mut live = load(&before)
        .with_terminal_channels(vec![slit(SIDE)], 8)
        .unwrap();
    for _ in 0..30 {
        live.preview_step(1e-9).unwrap().accept();
    }
    let accepted = live.clone();
    assert!(
        accepted
            .terminal()
            .unwrap()
            .pole_states()
            .iter()
            .any(|z| *z != C64::ZERO)
    );
    assert!(accepted.terminal_channels()[0].stored_energy_j() > 0.0);
    let report = live
        .rebind_gas_field_controlled(&after, DgGasRefresh::Lossless, 1.0, &mut || false)
        .unwrap();
    assert_eq!(live.state().step_count(), accepted.state().step_count());
    assert_eq!(live.branches()[0].state(), accepted.branches()[0].state());
    assert_eq!(
        live.branches()[0].duct().gas_field(),
        accepted.branches()[0].duct().gas_field()
    );
    assert_eq!(
        live.terminal().unwrap().pole_states(),
        accepted.terminal().unwrap().pole_states()
    );
    assert_eq!(live.last_terminal_frames(), accepted.last_terminal_frames());
    let gap = &live.terminal_channels()[0];
    let old_gap = &accepted.terminal_channels()[0];
    assert_eq!(gap.gas(), old_gap.gas());
    assert_eq!(
        gap.channel().unwrap().velocities_m_s(),
        old_gap.channel().unwrap().velocities_m_s()
    );
    assert_eq!(
        gap.cumulative_pressure_work_j(),
        old_gap.cumulative_pressure_work_j()
    );
    near(
        live.stored_energy_j() - accepted.stored_energy_j(),
        report.parameter_work_j,
        live.stored_energy_j() + accepted.stored_energy_j(),
    );
    assert_ne!(live.impedance_pa_s_m3(), accepted.impedance_pa_s_m3());
    let mut unchanged = accepted;
    let mut changed_signal = 0.0;
    for _ in 0..20 {
        let a = live.preview_step(0.0).unwrap().accept();
        let b = unchanged.preview_step(0.0).unwrap().accept();
        changed_signal +=
            (a.port_outputs[0].conjugate_rate - b.port_outputs[0].conjugate_rate).abs();
        let frame = &live.last_branch_frames()[0];
        near(
            frame.bore_port.conjugate_rate,
            frame.branch_inlet.conjugate_rate,
            frame.bore_port.conjugate_rate.abs() + frame.branch_inlet.conjugate_rate.abs(),
        );
    }
    assert!(changed_signal > 1e-8);
}

#[test]
fn g4_incompatible_branch_pressure_and_terminal_gas_refuse_without_resetting_history() {
    let before = field(300.0, 101_325.0, 0.05);
    for (port, replacement) in [
        (SIDE, field(450.0, 120_000.0, 0.05)),
        (OUTLET, field(450.0, 101_325.0, 0.05)),
    ] {
        let mut live = load(&before)
            .with_terminal_channels(vec![slit(port)], 8)
            .unwrap();
        for _ in 0..10 {
            live.preview_step(1e-9).unwrap().accept();
        }
        let accepted = live.clone();
        assert!(
            live.rebind_gas_field_controlled(
                &replacement,
                DgGasRefresh::Lossless,
                1.0,
                &mut || false,
            )
            .is_err()
        );
        assert_eq!(live.state(), accepted.state());
        assert_eq!(
            live.duct().operator_identity(),
            accepted.duct().operator_identity()
        );
        assert_eq!(live.branches()[0].state(), accepted.branches()[0].state());
        assert_eq!(
            live.port_response().unwrap(),
            accepted.port_response().unwrap()
        );
        assert_eq!(
            live.terminal().unwrap().pole_states(),
            accepted.terminal().unwrap().pole_states()
        );
        assert_eq!(
            live.terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s(),
            accepted.terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s()
        );
        let mut retry = accepted;
        assert_eq!(
            live.preview_step(-1e-9).unwrap().accept(),
            retry.preview_step(-1e-9).unwrap().accept()
        );
    }
}
