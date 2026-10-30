//! Regional source changes reach the continuing original moving valve.
use super::*;
use fs_couple::dg_port::{DgBranch, DgBranchGasInput};
use fs_duct::dg::DgSidePort;

const BRANCH: PortId = PortId(8);
const BRANCH_LENGTH: f64 = 0.015;

fn branched(main: &AxialGasField, branch: &AxialGasField) -> DgPortLoad {
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: main.length_m(),
            }],
        },
        main,
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::Port(PortId(1)),
    )
    .unwrap()
    .with_side_ports(&[DgSidePort {
        id: BRANCH,
        downstream_element: 1,
    }])
    .unwrap();
    let branch_duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: branch.length_m(),
            }],
        },
        branch,
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(20)),
        DgBoundary::Port(PortId(21)),
    )
    .unwrap();
    let terminal = PreparedMimoLoad::new(
        &MimoModel {
            ports: 2,
            poles: vec![PoleTerm::Real {
                pole: -2000.0,
                residue: 0.0,
            }],
            residues: vec![vec![
                C64::from_re(2e11),
                C64::ZERO,
                C64::ZERO,
                C64::from_re(1e11),
            ]],
            d: vec![1e8, 0.0, 0.0, 2e8],
            e: vec![0.0; 4],
            reciprocal: true,
        },
        DT,
        MimoLoadLimits {
            band_rad_s: [1.0, core::f64::consts::PI / DT],
            sampled_passivity_points: 257,
            max_states: 8,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap();
    DgPortLoad::new_with_terminal_branches(
        duct,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        vec![PortId(1), BRANCH],
        terminal,
        vec![DgBranch {
            port: BRANCH,
            duct: branch_duct,
        }],
    )
    .unwrap()
    .with_terminal_channels(
        vec![DgTerminalChannelInput {
            port: BRANCH,
            kind: DgTerminalChannelKind::Circular {
                radius_m: 4e-5,
                length_m: 0.001,
                radial_cells: 8,
                maximum_mach: 0.2,
            },
        }],
        8,
    )
    .unwrap()
}

fn input(field: &AxialGasField) -> [DgBranchGasInput<'_>; 1] {
    [DgBranchGasInput {
        port: BRANCH,
        field,
        refresh: DgGasRefresh::Lossless,
    }]
}

fn unchanged(load: &DgPortLoad, old: &DgPortLoad) {
    assert_eq!(load.state(), old.state());
    assert_eq!(
        load.duct().operator_identity(),
        old.duct().operator_identity()
    );
    assert_eq!(load.duct().gas_field(), old.duct().gas_field());
    assert_eq!(load.branches()[0].state(), old.branches()[0].state());
    assert_eq!(
        load.branches()[0].duct().gas_field(),
        old.branches()[0].duct().gas_field()
    );
    assert_eq!(
        load.terminal().unwrap().pole_states(),
        old.terminal().unwrap().pole_states()
    );
    assert_eq!(load.last_branch_frames(), old.last_branch_frames());
    assert_eq!(load.last_terminal_frames(), old.last_terminal_frames());
    assert_eq!(
        load.terminal_channels()[0].gas(),
        old.terminal_channels()[0].gas()
    );
    assert_eq!(
        load.terminal_channels()[0]
            .channel()
            .unwrap()
            .velocities_m_s(),
        old.terminal_channels()[0]
            .channel()
            .unwrap()
            .velocities_m_s()
    );
    assert_eq!(
        load.cumulative_parameter_work_j(),
        old.cumulative_parameter_work_j()
    );
    assert_eq!(load.impedance_pa_s_m3(), old.impedance_pa_s_m3());
}

#[test]
fn g1_branch_only_gas_change_reaches_the_scalar_source_without_resetting_motion() {
    let main = field(300.0, 1.0, 101325.0, LENGTH);
    let branch = field(350.0, 0.5, 101325.0, BRANCH_LENGTH);
    let next = field(460.0, 0.2, 101325.0, BRANCH_LENGTH);
    let build = || {
        let load = branched(&main, &branch);
        ApertureDg::new(aperture(&load, load.impedance_pa_s_m3()), load).unwrap()
    };
    let (mut changed, mut frozen) = (build(), build());
    for sample in 0..24 {
        assert_eq!(
            changed.step(drive(sample)).unwrap(),
            frozen.step(drive(sample)).unwrap()
        );
    }
    let old = changed.load().clone();
    let motion = changed.aperture().state();
    let mechanical_energy = changed.aperture().stored_energy_j();
    let viscosity = changed.aperture().viscous_slit();
    let report = changed
        .rebind_regional_gas_fields_controlled(None, &input(&next), 1.0, &mut || false)
        .unwrap();
    assert!(report.main.is_none());
    assert_eq!(report.branches[0].port, BRANCH);
    assert_eq!(changed.duct_state(), old.state());
    assert_eq!(
        changed.load().duct().operator_identity(),
        old.duct().operator_identity()
    );
    assert_eq!(
        changed.load().terminal().unwrap().pole_states(),
        old.terminal().unwrap().pole_states()
    );
    assert_eq!(
        changed.load().last_branch_frames(),
        old.last_branch_frames()
    );
    assert_eq!(
        changed.load().last_terminal_frames(),
        old.last_terminal_frames()
    );
    assert_eq!(changed.load().branches()[0].duct().gas_field(), Some(&next));
    assert_eq!(changed.aperture().state(), motion);
    assert_eq!(changed.aperture().stored_energy_j(), mechanical_energy);
    assert_eq!(changed.aperture().viscous_slit(), viscosity);
    assert_eq!(
        changed.aperture().spec().density_kg_m3,
        old.duct().inlet_gas().density
    );
    assert_ne!(changed.load().impedance_pa_s_m3(), old.impedance_pa_s_m3());
    assert_eq!(
        changed.aperture().spec().impedance_pa_s_m3,
        changed.load().impedance_pa_s_m3()
    );
    close(
        report.parameter_work_j,
        changed.load().stored_energy_j() - old.stored_energy_j(),
        changed.load().stored_energy_j() + old.stored_energy_j(),
    );
    let mut difference = 0.0;
    for sample in 24..32 {
        let frame = changed.step(drive(sample)).unwrap();
        let reference = frozen.step(drive(sample)).unwrap();
        difference += (frame.aperture.bore_flow_m3_s - reference.aperture.bore_flow_m3_s).abs();
        close(
            frame.balance_residual_j(),
            0.0,
            frame.stored_energy_j
                + frame.upstream_work_j.abs()
                + frame.inferred_total_dissipation_j,
        );
    }
    assert!(
        difference > 0.0,
        "regional physical load must reach the next source solve"
    );
}

#[test]
fn g1_regional_loading_preserves_both_original_blades_and_their_source_viscosity() {
    let main = field(300.0, 1.0, 101325.0, LENGTH);
    let branch = field(350.0, 0.5, 101325.0, BRANCH_LENGTH);
    let next = field(460.0, 0.2, 101325.0, BRANCH_LENGTH);
    let load = branched(&main, &branch);
    let mut model =
        LoadedApertureDg::new_loaded(paired_source(&load, load.impedance_pa_s_m3()), load).unwrap();
    for sample in 0..24 {
        model.step(drive(sample)).unwrap();
    }
    let motion = *pair(&model).blade_states();
    let blades = pair(&model).spec().blades;
    let energy = pair(&model).stored_energy_j();
    let viscosity = pair(&model).viscous_slit();
    let steps = pair(&model).accepted_steps();
    assert_ne!(motion[0], motion[1]);
    model
        .rebind_regional_gas_fields_controlled(None, &input(&next), 1.0, &mut || false)
        .unwrap();
    assert_eq!(pair(&model).blade_states(), &motion);
    assert_eq!(pair(&model).spec().blades, blades);
    assert_eq!(pair(&model).stored_energy_j(), energy);
    assert_eq!(pair(&model).viscous_slit(), viscosity);
    assert_eq!(pair(&model).accepted_steps(), steps);
    assert_eq!(
        pair(&model).spec().impedance_pa_s_m3,
        model.load().impedance_pa_s_m3()
    );
    let frame = model.step(drive(24)).unwrap();
    let pressure_drop = drive(24).upstream_pressure_pa - frame.aperture.bore_pressure_pa;
    for i in 0..2 {
        // The actual new pressure drives each independent original oscillator.
        let b = blades[i];
        let omega = (b.stiffness_n_m / b.mass_kg).sqrt();
        let equilibrium = -b.pressure_area_m2 * pressure_drop / b.stiffness_n_m;
        let (sin, cos) = (omega * DT).sin_cos();
        let displacement = equilibrium
            + (motion[i].displacement_m - equilibrium) * cos
            + motion[i].velocity_m_s / omega * sin;
        let velocity =
            -(motion[i].displacement_m - equilibrium) * omega * sin + motion[i].velocity_m_s * cos;
        let amplitude = motion[i].displacement_m.abs()
            + equilibrium.abs()
            + motion[i].velocity_m_s.abs() / omega;
        close(
            pair(&model).blade_states()[i].displacement_m,
            displacement,
            amplitude,
        );
        close(
            pair(&model).blade_states()[i].velocity_m_s,
            velocity,
            omega * amplitude,
        );
    }
}

#[test]
fn g4_series_regional_refresh_keeps_one_source_epoch_and_retries_exactly() {
    let old = [
        field(300.0, 1.0, 101325.0, 0.035),
        field(300.0, 1.0, 101325.0, LENGTH),
    ];
    let branches = [
        field(350.0, 0.5, 101325.0, BRANCH_LENGTH),
        field(330.0, 0.3, 101325.0, BRANCH_LENGTH),
    ];
    let next = [
        field(420.0, 0.2, 125000.0, 0.035),
        field(420.0, 0.2, 125000.0, LENGTH),
    ];
    let next_branches = [
        field(440.0, 0.6, 125000.0, BRANCH_LENGTH),
        field(380.0, 0.8, 125000.0, BRANCH_LENGTH),
    ];
    let build = || {
        let up = branched(&old[0], &branches[0]);
        let down = branched(&old[1], &branches[1]);
        let (_, impedance) = series_inlet_state(&up, &down).unwrap();
        LoadedApertureDg::new_loaded_with_upstream(paired_source(&down, impedance), up, down)
            .unwrap()
    };
    let (mut changed, mut clean) = (build(), build());
    for sample in 0..24 {
        assert_eq!(
            changed.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
    let loads = [
        changed.upstream_load().unwrap().clone(),
        changed.load().clone(),
    ];
    let motion = *pair(&changed).blade_states();
    let viscosity = pair(&changed).viscous_slit();
    let source_spec = *pair(&changed).spec();
    let branch_inputs = [input(&next_branches[0]), input(&next_branches[1])];
    let main = [
        Some((&next[0], DgGasRefresh::Lossless)),
        Some((&next[1], DgGasRefresh::Lossless)),
    ];
    let mut polls = 0;
    let reports = clean
        .rebind_series_regional_gas_fields_controlled(
            main,
            [&branch_inputs[0], &branch_inputs[1]],
            1.0,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    // Each tract's own proposal fits, but the enclosing allowance is smaller
    // than their sum: only the shared source transaction can reject this case.
    let [a, b] = reports
        .each_ref()
        .map(|r| r.as_ref().unwrap().absolute_parameter_work_j);
    assert!(a > 0.0 && b > 0.0);
    let budget = a.max(b) + 0.5 * a.min(b);
    assert!(
        changed
            .rebind_series_regional_gas_fields_controlled(
                main,
                [&branch_inputs[0], &branch_inputs[1]],
                budget,
                &mut || false,
            )
            .is_err()
    );
    // The upstream region update succeeds privately, but the unchanged inlet
    // on the other side cannot participate in the original volume-flow valve.
    assert!(
        changed
            .rebind_series_regional_gas_fields_controlled(
                [main[0], None],
                [&branch_inputs[0], &[]],
                1.0,
                &mut || false,
            )
            .is_err()
    );
    let mut calls = 0;
    assert!(
        changed
            .rebind_series_regional_gas_fields_controlled(
                main,
                [&branch_inputs[0], &branch_inputs[1]],
                1.0,
                &mut || {
                    calls += 1;
                    calls == polls
                },
            )
            .is_err()
    );
    assert_eq!(calls, polls);
    unchanged(changed.upstream_load().unwrap(), &loads[0]);
    unchanged(changed.load(), &loads[1]);
    assert_eq!(pair(&changed).blade_states(), &motion);
    assert_eq!(pair(&changed).viscous_slit(), viscosity);
    assert_eq!(pair(&changed).spec(), &source_spec);
    assert_eq!(
        changed
            .rebind_series_regional_gas_fields_controlled(
                main,
                [&branch_inputs[0], &branch_inputs[1]],
                1.0,
                &mut || false,
            )
            .unwrap(),
        reports
    );
    let gas = changed.load().duct().inlet_gas();
    assert_eq!(changed.upstream_load().unwrap().duct().inlet_gas(), gas);
    assert_eq!(pair(&changed).gas(), Some(gas));
    assert_eq!(pair(&changed).viscous_slit().unwrap().gas(), Some(gas));
    assert_eq!(
        pair(&changed).spec().impedance_pa_s_m3,
        changed.load().impedance_pa_s_m3() + changed.upstream_load().unwrap().impedance_pa_s_m3()
    );
    for sample in 24..32 {
        assert_eq!(
            changed.step(drive(sample)).unwrap(),
            clean.step(drive(sample)).unwrap()
        );
    }
    // An explicitly unselected tract remains exactly unchanged, while a remote
    // branch on its peer can still alter the valve's summed mean impedance.
    let before_down = changed.load().clone();
    let branch_only = field(470.0, 0.4, 125000.0, BRANCH_LENGTH);
    let partial = changed
        .rebind_series_regional_gas_fields_controlled(
            [None; 2],
            [&input(&branch_only), &[]],
            1.0,
            &mut || false,
        )
        .unwrap();
    assert!(partial[0].is_some());
    assert!(partial[1].is_none());
    unchanged(changed.load(), &before_down);
    assert_eq!(
        pair(&changed).spec().impedance_pa_s_m3,
        changed.load().impedance_pa_s_m3() + changed.upstream_load().unwrap().impedance_pa_s_m3()
    );
}
