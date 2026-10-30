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
fn g4_incompatible_branch_pressure_and_unfunded_gas_refuse_without_resetting_history() {
    let before = field(300.0, 101_325.0, 0.05);
    for (port, replacement, allowance) in [
        (SIDE, field(450.0, 120_000.0, 0.05), 1.0),
        (OUTLET, field(450.0, 101_325.0, 0.05), 0.0),
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
                allowance,
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

// Actual FV cell volumes, independently integrated from declared geometry.
fn volumes(kind: DgTerminalChannelKind) -> Vec<f64> {
    match kind {
        DgTerminalChannelKind::Slit {
            width_m,
            height_m,
            length_m,
            transverse_cells,
            ..
        } => vec![width_m * height_m * length_m / transverse_cells as f64; transverse_cells],
        DgTerminalChannelKind::Circular {
            radius_m,
            length_m,
            radial_cells,
            ..
        } => {
            let dr = radius_m / radial_cells as f64;
            (0..radial_cells)
                .map(|i| {
                    core::f64::consts::PI * length_m * dr * dr * ((i + 1).pow(2) - i.pow(2)) as f64
                })
                .collect()
        }
        DgTerminalChannelKind::Closed => Vec::new(),
    }
}

#[test]
fn g1_live_terminal_gas_preserves_cell_momenta_and_closes_aggregate_work() {
    let before = field(300.0, 101_325.0, 0.05);
    let after = field(450.0, 101_325.0, 0.05);
    let round = DgTerminalChannelInput {
        port: OUTLET,
        kind: DgTerminalChannelKind::Circular {
            radius_m: 8e-5,
            length_m: 0.002,
            radial_cells: 8,
            maximum_mach: 0.2,
        },
    };
    for outlet in [slit(OUTLET), round] {
        let mut live = load(&before)
            .with_terminal_channels(vec![slit(SIDE), outlet], 16)
            .unwrap();
        for _ in 0..30 {
            live.preview_step(1e-9).unwrap().accept();
        }
        let accepted = live.clone();
        let old = &accepted.terminal_channels()[1];
        assert!(old.stored_energy_j() > 0.0);
        let old_velocities = old.channel().unwrap().velocities_m_s();
        let weights = volumes(outlet.kind);
        let old_energy: f64 = weights
            .iter()
            .zip(old_velocities)
            .map(|(volume, velocity)| 0.5 * old.gas().density * volume * velocity * velocity)
            .sum();
        near(old_energy, old.stored_energy_j(), old_energy);
        let report = live
            .rebind_gas_field_controlled(&after, DgGasRefresh::Lossless, 1.0, &mut || false)
            .unwrap();
        let channel = &live.terminal_channels()[1];
        let new_velocities = channel.channel().unwrap().velocities_m_s();
        let new_energy: f64 = weights
            .iter()
            .zip(new_velocities)
            .map(|(volume, velocity)| 0.5 * channel.gas().density * volume * velocity * velocity)
            .sum();
        for ((volume, old_speed), new_speed) in
            weights.iter().zip(old_velocities).zip(new_velocities)
        {
            let old_momentum = old.gas().density * volume * old_speed;
            near(
                channel.gas().density * volume * new_speed,
                old_momentum,
                old_momentum.abs(),
            );
        }
        assert_eq!(channel.gas(), after.uniform_state().unwrap().state());
        assert_eq!(
            channel.channel().unwrap().last_step(),
            old.channel().unwrap().last_step()
        );
        assert_eq!(
            channel.cumulative_pressure_work_j(),
            old.cumulative_pressure_work_j()
        );
        assert_eq!(
            channel.cumulative_viscous_dissipation_j(),
            old.cumulative_viscous_dissipation_j()
        );
        assert_eq!(
            channel.cumulative_numerical_dissipation_j(),
            old.cumulative_numerical_dissipation_j()
        );
        assert_eq!(
            live.last_terminal_channel_frames(),
            accepted.last_terminal_channel_frames()
        );
        assert_eq!(
            live.terminal().unwrap().pole_states(),
            accepted.terminal().unwrap().pole_states()
        );
        assert_eq!(live.branches()[0].state(), accepted.branches()[0].state());
        let branch_channel = &live.terminal_channels()[0];
        assert_eq!(branch_channel.gas(), accepted.terminal_channels()[0].gas());
        assert_eq!(
            branch_channel.channel().unwrap().velocities_m_s(),
            accepted.terminal_channels()[0]
                .channel()
                .unwrap()
                .velocities_m_s()
        );
        assert_eq!(
            report.terminal_channels[0]
                .transfer
                .unwrap()
                .parameter_work_j,
            0.0
        );
        let transfer = report.terminal_channels[1].transfer.unwrap();
        near(transfer.kinetic_energy_before_j, old_energy, old_energy);
        near(transfer.kinetic_energy_after_j, new_energy, new_energy);
        near(
            transfer.parameter_work_j,
            new_energy - old_energy,
            new_energy + old_energy,
        );
        assert_eq!(
            report.parameter_work_j,
            report.duct.parameter_work_j + transfer.parameter_work_j
        );
        assert_eq!(
            report.absolute_parameter_work_j,
            report.duct.parameter_work_j.abs() + transfer.parameter_work_j.abs()
        );
        near(
            report.energy_after_j - report.energy_before_j,
            report.parameter_work_j,
            report.energy_after_j + report.energy_before_j,
        );
        assert_eq!(live.cumulative_parameter_work_j(), report.parameter_work_j);
        for _ in 0..12 {
            let old_storage = live.stored_energy_j();
            let old_loss = live.cumulative_dissipated_j();
            let trial = live.preview_step(-2e-10).unwrap();
            let inlet_work = trial.inlet_work_j();
            let exterior_work = sum_terminal_work(trial.terminal_frames());
            let output = trial.accept();
            let storage = live.stored_energy_j() - old_storage;
            let loss = live.cumulative_dissipated_j() - old_loss;
            near(
                storage + loss + exterior_work,
                inlet_work,
                storage.abs() + loss.abs() + exterior_work.abs() + inlet_work.abs(),
            );
            assert_eq!(
                output.ledger.cumulative_parameter_work_j,
                report.parameter_work_j
            );
            near(
                output.ledger.energy_balance_residual_j,
                0.0,
                output.ledger.stored_energy_j
                    + output.ledger.cumulative_port_work_j.abs()
                    + output.ledger.cumulative_dissipated_j
                    + report.absolute_parameter_work_j,
            );
        }
    }
}

#[test]
fn g4_channel_work_budget_and_final_cancel_leave_all_gas_owners_unchanged() {
    let before = field(300.0, 101_325.0, 0.05);
    let after = field(450.0, 101_325.0, 0.05);
    let mut accepted = load(&before)
        .with_terminal_channels(vec![slit(OUTLET)], 8)
        .unwrap();
    for _ in 0..30 {
        accepted.preview_step(1e-9).unwrap().accept();
    }
    let mut clean = accepted.clone();
    let mut polls = 0;
    let report = clean
        .rebind_gas_field_controlled(&after, DgGasRefresh::Lossless, 1.0, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let channel_work = report.terminal_channels[0]
        .transfer
        .unwrap()
        .parameter_work_j;
    assert!(channel_work != 0.0);
    let allowance = report.duct.parameter_work_j.abs().max(channel_work.abs());
    assert!(report.absolute_parameter_work_j > allowance);
    for cancel_at_end in [false, true] {
        let mut live = accepted.clone();
        let mut count = 0;
        assert!(
            live.rebind_gas_field_controlled(
                &after,
                DgGasRefresh::Lossless,
                if cancel_at_end { 1.0 } else { allowance },
                &mut || {
                    count += 1;
                    cancel_at_end && count == polls
                }
            )
            .is_err()
        );
        if cancel_at_end {
            assert_eq!(count, polls);
        }
        assert_eq!(live.state(), accepted.state());
        assert_eq!(
            live.port_response().unwrap(),
            accepted.port_response().unwrap()
        );
        assert_eq!(live.duct().gas_field(), accepted.duct().gas_field());
        assert_eq!(
            live.terminal().unwrap().pole_states(),
            accepted.terminal().unwrap().pole_states()
        );
        assert_eq!(
            live.terminal_channels()[0].gas(),
            accepted.terminal_channels()[0].gas()
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
        assert_eq!(
            live.terminal_channels()[0].cumulative_parameter_work_j(),
            0.0
        );
        assert_eq!(
            live.rebind_gas_field_controlled(&after, DgGasRefresh::Lossless, 1.0, &mut || false)
                .unwrap(),
            report
        );
        let mut reference = clean.clone();
        assert_eq!(
            live.preview_step(0.0).unwrap().accept(),
            reference.preview_step(0.0).unwrap().accept()
        );
    }
}

#[test]
fn g3_closed_terminal_rebinds_source_without_fabricating_fluid_storage() {
    let before = field(300.0, 101_325.0, 0.05);
    let after = field(450.0, 101_325.0, 0.05);
    let mut live = load(&before)
        .with_terminal_channels(
            vec![DgTerminalChannelInput {
                port: OUTLET,
                kind: DgTerminalChannelKind::Closed,
            }],
            0,
        )
        .unwrap();
    for _ in 0..8 {
        live.preview_step(1e-9).unwrap().accept();
    }
    let report = live
        .rebind_gas_field_controlled(&after, DgGasRefresh::Lossless, 1.0, &mut || false)
        .unwrap();
    assert_eq!(report.terminal_channels[0].transfer, None);
    assert_eq!(live.terminal_channels()[0].stored_energy_j(), 0.0);
    assert_eq!(
        live.terminal_channels()[0].cumulative_parameter_work_j(),
        0.0
    );
    assert_eq!(
        live.terminal_channels()[0].gas(),
        after.uniform_state().unwrap().state()
    );
    live.preview_step(0.0).unwrap().accept();
    assert_eq!(
        live.last_terminal_channel_frames()[0].outward_flow_m3_s,
        0.0
    );
}

mod molecular;

