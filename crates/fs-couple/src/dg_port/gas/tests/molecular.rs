//! Same complete port transaction, with the original molecular gas histories.
use super::*;
use fs_duct::dg::DgWallResolution;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn molecular_load(source: &AxialGasField) -> DgPortLoad {
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: SpeciesId::new("live-port-gas").unwrap(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                theta_k: 900.0,
                degeneracy: 1,
            }]),
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 3.0,
                t_ref_k: 300.0,
                exponent: 1.2,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }],
    }])
    .unwrap();
    // This is the same explicit spatial-wall route used by dg_bore gas_field.
    // Its original quadrature is retained even for a uniform cylinder.
    let main = duct(source, 2, INLET, OUTLET)
        .with_side_ports(&[DgSidePort {
            id: SIDE,
            downstream_element: 1,
        }])
        .unwrap()
        .with_spatial_viscothermal(DgWallResolution {
            poles_per_field: 8,
            quadrature_points: 3,
            band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
            maximum_relative_error: 0.02,
        })
        .unwrap()
        .with_molecular_relaxation_field(law)
        .unwrap();
    let template = load(source);
    DgPortLoad::new_with_terminal_branches(
        main,
        DT,
        1e-11,
        MAX_DG_PROPAGATOR_DOFS,
        vec![SIDE, OUTLET],
        template.terminal().unwrap().clone(),
        vec![DgBranch {
            port: SIDE,
            duct: template.branches()[0].duct().clone(),
        }],
    )
    .unwrap()
    .with_terminal_channels(vec![slit(SIDE)], 8)
    .unwrap()
}

#[test]
fn molecular_refresh_keeps_the_complete_branch_and_mutual_terminal_transaction() {
    let before = field(300.0, 101_325.0, 0.05);
    let after = field(420.0, 101_325.0, 0.05);
    let mut live = molecular_load(&before);
    for _ in 0..20 {
        live.preview_step(1e-9).unwrap().accept();
    }
    let accepted = live.clone();
    assert!(
        accepted
            .duct()
            .molecular_dissipated_power_w(accepted.state())
            .unwrap()
            > 0.0
    );
    assert!(
        accepted
            .terminal()
            .unwrap()
            .pole_states()
            .iter()
            .any(|z| *z != C64::ZERO)
    );
    let mut original_duct = accepted.duct().clone();
    let mut original_state = accepted.state().clone();
    let native = original_duct
        .rebind_uniform_molecular_gas_field(&mut original_state, &after, 0.02, 1.0, &mut || false)
        .unwrap();
    let report = live
        .rebind_gas_field_controlled(
            &after,
            DgGasRefresh::UniformMolecularDiffusion {
                maximum_relative_error: 0.02,
            },
            1.0,
            &mut || false,
        )
        .unwrap();
    assert_eq!(report.duct, native);
    assert_eq!(live.state(), &original_state);
    assert_eq!(
        live.duct().operator_identity(),
        original_duct.operator_identity()
    );
    assert_eq!(live.state().step_count(), accepted.state().step_count());
    assert_eq!(live.branches()[0].state(), accepted.branches()[0].state());
    assert_eq!(
        live.terminal().unwrap().pole_states(),
        accepted.terminal().unwrap().pole_states()
    );
    assert_eq!(live.last_terminal_frames(), accepted.last_terminal_frames());
    let channel = &live.terminal_channels()[0];
    let old_channel = &accepted.terminal_channels()[0];
    assert_eq!(channel.gas(), old_channel.gas());
    assert_eq!(
        channel.channel().unwrap().velocities_m_s(),
        old_channel.channel().unwrap().velocities_m_s()
    );
    near(
        report.energy_after_j - report.energy_before_j,
        report.parameter_work_j,
        report.energy_after_j + report.energy_before_j,
    );
    let mut reference = accepted;
    let mut difference = 0.0;
    for _ in 0..12 {
        let updated = live.preview_step(0.0).unwrap().accept();
        let unchanged = reference.preview_step(0.0).unwrap().accept();
        difference += (updated.port_outputs[0].conjugate_rate
            - unchanged.port_outputs[0].conjugate_rate)
            .abs();
        assert_eq!(
            updated.ledger.cumulative_parameter_work_j,
            report.parameter_work_j
        );
        near(
            updated.ledger.energy_balance_residual_j,
            0.0,
            updated.ledger.stored_energy_j
                + updated.ledger.cumulative_port_work_j.abs()
                + updated.ledger.cumulative_dissipated_j
                + report.absolute_parameter_work_j,
        );
    }
    assert!(difference > 1e-8);
}
