//! The accepted SGS field closes heat transport as well as nonlinear stress.

use super::*;
use super::super::super::super::super::momentum::coupled_controls;
use crate::conjugate::FvSgsConfig;
use crate::sgs3::{SgsModel, SgsUse};

fn model_config(model: SgsModel) -> FvSgsConfig {
    FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 1.0,
        usage: SgsUse::Exploration,
    }
}

/// Eliminate the original three molecular-heating caloric rows at native
/// faces [0,b,-b,0]. Only fluid-fluid conductivity receives rho*Cp*nu_t/Pr.
/// Half-cell series resistance gives t=dt*G=ke*kc/(ke+kc), not an arithmetic
/// conductivity average. The optional frozen speed supplies a wrong-closure
/// comparison without reusing any production thermal or momentum solver.
fn thermal_offsets(b: f64, cs: f64, pr: f64, frozen: Option<f64>) -> [f64; 3] {
    let delta_k = 2.0 * cs * cs * 2.0_f64.sqrt() * frozen.unwrap_or(b) / pr;
    let ke = 2.0 + delta_k;
    let kc = 2.0 + 2.0 * delta_k;
    let t = ke * kc / (ke + kc);
    let alpha = 2.0 + t - (31.0 / 6.0) * b * b;
    let beta = 2.0 + 2.0 * t - (62.0 / 3.0) * b * b;
    let edge = 0.25 * b * b + (17.0 / 48.0) * b.powi(3);
    let center = b * b + (17.0 / 6.0) * b.powi(3);
    let determinant = alpha * beta - 2.0 * t * t;
    [
        (edge * beta + t * center) / determinant,
        (alpha * center + 2.0 * t * edge) / determinant,
        t,
    ]
}

fn reference(cs: f64, pr: f64, frozen: Option<f64>) -> (f64, [f64; 3], f64) {
    let e_sgs = (8.0 / 3.0) * cs * cs * 2.0_f64.sqrt();
    let (mut low, mut high) = (0.0, 0.3);
    for _ in 0..64 {
        let b = 0.5 * (low + high);
        let [e, c, _] = thermal_offsets(b, cs, pr, frozen);
        let le = 0.5 + (31.0 / 3.0) * e + (17.0 / 24.0) * b;
        let lc = 0.5 + (31.0 / 3.0) * c + (17.0 / 12.0) * b;
        if b * (1.0 + le / 4.0 + lc / 2.0 + 1.25 * e_sgs * b) > 0.5 {
            high = b;
        } else {
            low = b;
        }
    }
    let b = 0.5 * (low + high);
    let [e, c, t] = thermal_offsets(b, cs, pr, frozen);
    (b, [1.0 + e, 1.0 + c, 1.0 + e], t)
}

#[test]
fn turbulent_links_use_the_accepted_rate_and_prandtl_with_molecular_solid_interfaces() {
    let (material, melting) = thermal_material_with_rate(true);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    let cs = 0.3;
    let quadratic = 1.25 * (17.0 / 24.0 + (8.0 / 3.0) * cs * cs * 2.0_f64.sqrt());
    // Original-T nonlinear predictor, before any conductive feedback.
    let predictor = 1.0 / (1.375 + (1.375_f64.powi(2) + 2.0 * quadratic).sqrt());
    let mut previous: Option<[f64; 3]> = None;
    for pr in [0.2, 1.0] {
        let (b, expected_h, t) = reference(cs, pr, None);
        let (_, frozen_h, _) = reference(cs, pr, Some(predictor));
        assert!((expected_h[1] - frozen_h[1]).abs() > 1e-5);
        if let Some(stronger_mixing) = previous {
            assert!(expected_h[1] > stronger_mixing[1]);
            assert!(expected_h[0] < stronger_mixing[0]);
        }
        previous = Some(expected_h);
        // The solid link keeps the ORIGINAL k=2 on both half-cells.
        let expected_solid_h = (2.0 * 1.125 + expected_h[2]) / 3.0;
        let wall_heat = expected_solid_h - expected_h[2];
        for finite_wall in [false, true] {
            let n = if finite_wall { 4 } else { 3 };
            let mut domain = VoxelDomain::new(n, 1, 1, 1.0).unwrap();
            let mut mass = vec![2.0; n];
            let mut phases = vec![Some(caloric.state_at_specific_enthalpy(1.0).unwrap()); n];
            let mut authored = vec![0.0; n];
            let mut solids = vec![];
            if finite_wall {
                domain.set(3, 0, 0, Voxel::Solid(0));
                mass[3] = 0.0;
                phases[3] = None;
                authored[2] = -wall_heat;
                solids.push(ThermalVolumeSolidCell {
                    cell: 3,
                    conductivity: &conductivity,
                    mass_kg: 2.0,
                    phase: caloric.state_at_specific_enthalpy(1.125).unwrap(),
                    net_enthalpy_increment_j: 0.0,
                });
            }
            let field = FlowField::quiescent(&domain);
            let mut momentum = vec![[0.0; 3]; n];
            momentum[0][0] = 2.0;
            momentum[2][0] = -2.0;
            let forces = vec![[0.0; 3]; n];
            let result = field
                .advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs_transport(
                    &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &solids, &[],
                    ThermalVolumeMomentumInput {
                        volume: input(&mass, &phases, &authored),
                        momentum_kg_m_s: &momentum,
                        acceleration_m_s2: [0.0; 3],
                        force_density_n_m3: &forces,
                        boundary_velocity_m_s: [None; 6],
                    },
                    [FvBoundary::Symmetry; 6], model_config(SgsModel::Smagorinsky { cs }), pr,
                    coupled_controls(0.5), molecular_controls(), thermal_solve_controls(),
                    feedback_controls(), &mut ImplicitThermalVolumeWorkspace::new(n).unwrap(),
                    &CancelGate::new_clock_free(),
                    |_, phase, rate| joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate),
                )
                .unwrap();
            let accepted = &result.evolution;
            let viscous = &accepted.evolution.evolution.viscosity;
            let volume = &accepted.evolution.evolution.evolution.thermal.volume;
            let work = &accepted.feedback;
            assert_eq!(work.turbulent_prandtl_number, Some(pr));
            assert_eq!(work.used_cell_eddy_viscosity_m2_s.len(), n);
            assert_eq!(work.cell_eddy_viscosity_m2_s.len(), n);
            let (mut total_heat, mut modeled) = (0.0, 0.0);
            for cell in 0..3 {
                let g = if cell == 1 { 2.0 * b } else { b };
                let rate = 2.0 * g / 3.0_f64.sqrt();
                let mu = 0.125 + 4.0 * (expected_h[cell] - 1.0) + g / 4.0;
                let bulk = 1.0 / 3.0 + 5.0 * (expected_h[cell] - 1.0) + 3.0 * g / 8.0;
                let heat = 0.5 * ((4.0 / 3.0) * mu + bulk) * g * g;
                let transfer = (4.0 / 3.0) * cs * cs * 2.0_f64.sqrt() * g.powi(3);
                total_heat += heat;
                modeled += transfer;
                near(result.endpoint_shear_rate_per_s[cell], rate);
                near(work.cell_eddy_viscosity_m2_s[cell], cs * cs * 2.0_f64.sqrt() * g);
                assert!(work.used_cell_eddy_viscosity_m2_s[cell].is_finite());
                assert!(work.used_cell_eddy_viscosity_m2_s[cell] >= 0.0);
                let source = accepted.endpoint_sources[cell].as_ref().unwrap();
                near(source.phase_state().specific_enthalpy_j_kg(), expected_h[cell]);
                near(source.transport_state().unwrap().dynamic_viscosity_pa_s(), mu);
                near(source.material().property(K_KEY).unwrap().value_si(), 2.0);
                near(accepted.bulk_viscosity_pa_s[cell], bulk);
                for property in source.material().properties() {
                    for (axis, expected) in [("T", 300.0 + expected_h[cell]), ("rate", rate)] {
                        near(property.answer().receipt.query_point.iter()
                            .find(|(name, _)| name == axis).unwrap().1, expected);
                    }
                }
                near(viscous.cell_molecular_heat_j[cell], heat);
                near(accepted.used_molecular_heat_j[cell], heat);
                near(viscous.cell_sgs_transfer_j[cell], transfer);
                near(volume.phase_states[cell].unwrap().specific_enthalpy_j_kg(), expected_h[cell]);
                let mut conductive = if cell == 1 {
                    2.0 * t * (expected_h[0] - expected_h[1])
                } else {
                    t * (expected_h[1] - expected_h[0])
                };
                if finite_wall && cell == 2 { conductive += wall_heat; }
                near(accepted.evolution.conductive_enthalpy_increment_j[cell], conductive);
                near(2.0 * (expected_h[cell] - 1.0), authored[cell] + heat + conductive);
            }
            near(viscous.molecular_heat_j, total_heat);
            near(viscous.sgs_transfer_j, modeled);
            near(viscous.native_kinetic_energy_change_j + total_heat + modeled
                + viscous.backward_euler_loss_j - viscous.wall_work_j, 0.0);
            let report = &accepted.evolution.report;
            near(report.conduction.solve.energy.external_heat_j, total_heat + authored.iter().sum::<f64>());
            near(volume.report.material_enthalpy_change_j + report.solid_enthalpy_change_j,
                total_heat + authored.iter().sum::<f64>());
            if finite_wall {
                near(accepted.evolution.solid_endpoints[0].phase.specific_enthalpy_j_kg(), expected_solid_h);
                near(report.solid_to_fluid_heat_j, wall_heat);
                near(report.solid_enthalpy_change_j, -wall_heat);
                assert_eq!(solids[0].phase.specific_enthalpy_j_kg(), 1.125);
                assert_eq!(work.used_cell_eddy_viscosity_m2_s[3], 0.0);
                assert_eq!(work.cell_eddy_viscosity_m2_s[3], 0.0);
            }
            // Compare physical link work, not coefficient equality: a zero
            // temperature difference legitimately makes a closure mismatch
            // thermally inactive. Here T=300+h, so each exact delta_h is delta_T.
            let dt_conductance = |eddy: &[f64], left: usize, right: usize| {
                let kl = 2.0 + 2.0 * eddy[left] / pr;
                let kr = 2.0 + 2.0 * eddy[right] / pr;
                kl * kr / (kl + kr)
            };
            let mut transported_defect = 0.0;
            for (left, right) in [(0, 1), (1, 2)] {
                let delta_t = accepted.endpoint_sources[right].as_ref().unwrap()
                    .phase_state().specific_enthalpy_j_kg()
                    - accepted.endpoint_sources[left].as_ref().unwrap()
                        .phase_state().specific_enthalpy_j_kg();
                let used = dt_conductance(&work.used_cell_eddy_viscosity_m2_s, left, right);
                let actual = dt_conductance(&work.cell_eddy_viscosity_m2_s, left, right);
                transported_defect += 2.0 * ((used - actual) * delta_t).abs();
            }
            assert!(transported_defect <= work.absolute_sgs_heat_coupling_defect_j + 2e-14);
            assert!(work.trials > 1);
            assert!(work.absolute_sgs_heat_coupling_defect_j.is_finite());
            assert!(work.absolute_sgs_heat_coupling_defect_j >= 0.0);
            assert!(work.absolute_heat_coupling_defect_j + work.absolute_sgs_heat_coupling_defect_j
                <= feedback_controls().maximum_heat_defect_j);
            assert_eq!(result.rheology_source_updates, work.molecular_sweeps + work.molecular_solves);
            assert_eq!(work.source_evaluations, 3 * result.rheology_source_updates);
            assert_eq!(field, FlowField::quiescent(&domain));
            assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
        }
    }
}

#[test]
fn disabled_controls_and_transport_refusals_preserve_original_workspace_retry() {
    let (material, melting) = thermal_material_with_rate(true);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 3];
    let phases = [Some(caloric.state_at_specific_enthalpy(1.0).unwrap()); 3];
    let momentum = [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0]];
    let data = ThermalVolumeMomentumInput {
        volume: input(&mass, &phases, &[0.0; 3]),
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]; 3],
        boundary_velocity_m_s: [None; 6],
    };
    let calls = Cell::new(0);
    let run = |model, pr, viscosity, cancel: bool,
               workspace: &mut ImplicitThermalVolumeWorkspace, gate: &CancelGate| {
        field.advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs_transport(
            &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &[], &[], data,
            [FvBoundary::Symmetry; 6], model_config(model), pr, coupled_controls(0.5), viscosity,
            thermal_solve_controls(), feedback_controls(), workspace, gate,
            |_, phase, rate| {
                calls.set(calls.get() + 1);
                if cancel && phase.specific_enthalpy_j_kg() > 1.0 {
                    gate.request();
                    return Err(refuse("thermal source refused after cancellation"));
                }
                joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate)
            },
        )
    };
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(3).unwrap();
    let original = field
        .advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths(
            &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &[], &[], data,
            [FvBoundary::Symmetry; 6], coupled_controls(0.5), molecular_controls(),
            thermal_solve_controls(), feedback_controls(), &mut workspace, &gate,
            |_, phase, rate| joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate),
        ).unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let disabled = run(model, 0.4, molecular_controls(), false, &mut workspace, &gate).unwrap();
        assert_eq!(disabled, original);
        assert_eq!(disabled.evolution.feedback.turbulent_prandtl_number, None);
        assert!(disabled.evolution.feedback.used_cell_eddy_viscosity_m2_s.is_empty());
        assert!(disabled.evolution.feedback.cell_eddy_viscosity_m2_s.is_empty());
    }
    let active = SgsModel::Smagorinsky { cs: 0.3 };
    for model in [active, SgsModel::Smagorinsky { cs: 0.0 }] {
        for invalid in [0.0, -0.1, f64::NAN, f64::INFINITY] {
            calls.set(0);
            assert!(matches!(run(model, invalid, molecular_controls(), false, &mut workspace, &gate),
                Err(ChtError::InvalidInput { .. })));
            assert_eq!(calls.get(), 0);
        }
    }
    calls.set(0);
    let expected = run(active, 0.4, molecular_controls(), false, &mut workspace, &gate).unwrap();
    let total = expected.evolution.feedback.molecular_sweeps;
    assert!(total > expected.evolution.evolution.evolution.viscosity.sweeps + 1);
    assert!(run(active, 0.4, MaterialViscosityControls { max_sweeps: total - 1, ..molecular_controls() },
        false, &mut workspace, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    calls.set(0);
    assert!(matches!(run(active, 0.4, molecular_controls(), true, &mut workspace, &cancelled),
        Err(ChtError::Cancelled)));
    assert!(calls.get() > 3);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 3]);
    assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
    assert_eq!(momentum, [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0]]);
    calls.set(0);
    let retry = run(active, 0.4, molecular_controls(), false, &mut workspace, &gate).unwrap();
    assert_eq!(retry, expected);
    assert_eq!(calls.get(), retry.evolution.feedback.source_evaluations);
}
