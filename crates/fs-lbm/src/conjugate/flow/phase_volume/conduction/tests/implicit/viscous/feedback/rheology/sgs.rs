//! One implicit thermal/native rheology endpoint with separate modeled transfer.

use super::*;
use super::super::super::super::super::momentum::coupled_controls;
use crate::conjugate::FvSgsConfig;
use crate::sgs3::{SgsModel, SgsUse};

fn sgs(model: SgsModel) -> FvSgsConfig {
    FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 1.0,
        usage: SgsUse::Exploration,
    }
}

/// For native faces [0,b,-b,0], let e,c be edge/central h-1. The real
/// source gives lambda=4mu/3+zeta=1/2+(31/3)*(h-1)+(17/24)*|du/dx|.
/// Eliminating the two original molecular-only heat rows gives this 2x2
/// system; SGS does not enter it. dt=1/2, k=2, Cp=1, M=2 in every cell.
fn thermal_offsets(b: f64) -> [f64; 2] {
    let alpha = 3.0 - (31.0 / 6.0) * b * b;
    let beta = 4.0 - (62.0 / 3.0) * b * b;
    let edge = 0.25 * b * b + (17.0 / 48.0) * b.powi(3);
    let center = b * b + (17.0 / 6.0) * b.powi(3);
    let determinant = alpha * beta - 2.0;
    [
        (edge * beta + center) / determinant,
        (alpha * center + 2.0 * edge) / determinant,
    ]
}

fn reference(cs: f64) -> (f64, [f64; 3]) {
    // The native original inertia closes b*(1+lambda_e/4+lambda_c/2
    // +5*E*b/4)=1/2, E=(8/3)*Cs^2*sqrt(2). A bracketed scalar root
    // is independent of the production nested thermal/momentum iteration.
    let e_sgs = (8.0 / 3.0) * cs * cs * 2.0_f64.sqrt();
    let (mut low, mut high) = (0.0, 0.25);
    for _ in 0..64 {
        let b = 0.5 * (low + high);
        let [e, c] = thermal_offsets(b);
        let le = 0.5 + (31.0 / 3.0) * e + (17.0 / 24.0) * b;
        let lc = 0.5 + (31.0 / 3.0) * c + (17.0 / 12.0) * b;
        if b * (1.0 + le / 4.0 + lc / 2.0 + 1.25 * e_sgs * b) > 0.5 {
            high = b;
        } else {
            low = b;
        }
    }
    let b = 0.5 * (low + high);
    let [e, c] = thermal_offsets(b);
    (b, [1.0 + e, 1.0 + c, 1.0 + e])
}

#[test]
fn temperature_rate_and_sgs_close_against_original_fluid_and_finite_wall_storage() {
    let (material, melting) = thermal_material_with_rate(true);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    let cs = 0.3;
    let (b, expected_h) = reference(cs);
    assert!((b - 0.25).abs() > 0.006);
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
            // Balance the independently derived finite-wall supply so both
            // geometries retain the same three-fluid scalar reference.
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
        let calls = Cell::new(0);
        let result = field
            .advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs(
                &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &solids, &[],
                ThermalVolumeMomentumInput {
                    volume: input(&mass, &phases, &authored),
                    momentum_kg_m_s: &momentum,
                    acceleration_m_s2: [0.0; 3],
                    force_density_n_m3: &forces,
                    boundary_velocity_m_s: [None; 6],
                },
                [FvBoundary::Symmetry; 6], sgs(SgsModel::Smagorinsky { cs }),
                coupled_controls(0.5), molecular_controls(), thermal_solve_controls(),
                feedback_controls(), &mut ImplicitThermalVolumeWorkspace::new(n).unwrap(),
                &CancelGate::new_clock_free(),
                |cell, phase, rate| {
                    assert!(cell < 3);
                    calls.set(calls.get() + 1);
                    joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate)
                },
            )
            .unwrap();
        let accepted = &result.evolution;
        let viscous = &accepted.evolution.evolution.viscosity;
        let volume = &accepted.evolution.evolution.evolution.thermal.volume;
        let (mut total_heat, mut total_transfer) = (0.0, 0.0);
        for cell in 0..3 {
            let gradient = if cell == 1 { 2.0 * b } else { b };
            let rate = 2.0 * gradient / 3.0_f64.sqrt();
            let mu = 0.125 + 4.0 * (expected_h[cell] - 1.0) + gradient / 4.0;
            let bulk = 1.0 / 3.0 + 5.0 * (expected_h[cell] - 1.0) + 3.0 * gradient / 8.0;
            let mu_t = 2.0 * cs * cs * 2.0_f64.sqrt() * gradient;
            let heat = 0.5 * ((4.0 / 3.0) * mu + bulk) * gradient * gradient;
            let transfer = 0.5 * (4.0 / 3.0) * mu_t * gradient * gradient;
            total_heat += heat;
            total_transfer += transfer;
            near(result.endpoint_shear_rate_per_s[cell], rate);
            let source = accepted.endpoint_sources[cell].as_ref().unwrap();
            near(source.phase_state().specific_enthalpy_j_kg(), expected_h[cell]);
            near(source.transport_state().unwrap().dynamic_viscosity_pa_s(), mu);
            near(accepted.bulk_viscosity_pa_s[cell], bulk);
            for property in source.material().properties() {
                for (axis, expected) in [("T", 300.0 + expected_h[cell]), ("rate", rate)] {
                    let actual = property.answer().receipt.query_point.iter()
                        .find(|(name, _)| name == axis).unwrap().1;
                    near(actual, expected);
                }
            }
            near(viscous.cell_molecular_heat_j[cell], heat);
            near(accepted.used_molecular_heat_j[cell], heat);
            near(viscous.cell_sgs_transfer_j[cell], transfer);
            near(volume.phase_states[cell].unwrap().specific_enthalpy_j_kg(), expected_h[cell]);
            near(volume.transport.mass_kg[cell], 2.0);
            near(
                accepted.evolution.conductive_enthalpy_increment_j[cell],
                2.0 * (expected_h[cell] - 1.0) - authored[cell] - heat,
            );
        }
        near(viscous.molecular_heat_j, total_heat);
        near(viscous.sgs_transfer_j, total_transfer);
        assert!(total_transfer > 0.02);
        near(viscous.native_kinetic_energy_change_j, 2.0 * b * b - 0.5);
        near(viscous.backward_euler_loss_j, 2.0 * (b - 0.5).powi(2));
        near(
            viscous.native_kinetic_energy_change_j + total_heat + total_transfer
                + viscous.backward_euler_loss_j - viscous.wall_work_j,
            0.0,
        );
        let report = &accepted.evolution.report;
        near(
            report.conduction.solve.energy.external_heat_j,
            total_heat + authored.iter().sum::<f64>(),
        );
        near(
            volume.report.material_enthalpy_change_j + report.solid_enthalpy_change_j,
            total_heat + authored.iter().sum::<f64>(),
        );
        if finite_wall {
            near(
                accepted.evolution.solid_endpoints[0].phase.specific_enthalpy_j_kg(),
                expected_solid_h,
            );
            near(report.solid_to_fluid_heat_j, wall_heat);
            near(report.solid_enthalpy_change_j, -wall_heat);
            assert_eq!(solids[0].phase.specific_enthalpy_j_kg(), 1.125);
            assert!(accepted.endpoint_sources[3].is_none());
            assert_eq!(result.endpoint_shear_rate_per_s[3], 0.0);
            assert_eq!(viscous.cell_sgs_transfer_j[3], 0.0);
        }
        let work = &accepted.feedback;
        assert!(work.trials > 1);
        assert_eq!(work.molecular_solves, work.trials + 1);
        assert_eq!(result.rheology_source_updates, work.molecular_sweeps + work.molecular_solves);
        assert_eq!(calls.get(), 3 * result.rheology_source_updates);
        assert_eq!(work.source_evaluations, calls.get());
        assert!(work.molecular_sweeps <= molecular_controls().max_sweeps);
        assert!(work.absolute_heat_coupling_defect_j <= feedback_controls().maximum_heat_defect_j);
        assert_eq!(field, FlowField::quiescent(&domain));
        assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
        assert_eq!(momentum[0], [2.0, 0.0, 0.0]);
    }
}

#[test]
fn all_disabled_models_retain_the_exact_original_implicit_source_and_heat_history() {
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
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(3).unwrap();
    let resolve = |_: usize, phase: &EquilibriumPhaseState, rate| {
        joint_response(&material, &law, phase.specific_enthalpy_j_kg(), rate)
    };
    let original = field
        .advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths(
            &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &[], &[], data,
            [FvBoundary::Symmetry; 6], coupled_controls(0.5), molecular_controls(),
            thermal_solve_controls(), feedback_controls(), &mut workspace, &gate, resolve,
        ).unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let result = field
            .advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs(
                &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &[], &[], data,
                [FvBoundary::Symmetry; 6], FvSgsConfig { maximum_eddy_to_molecular: 0.0, ..sgs(model) },
                coupled_controls(0.5), molecular_controls(), thermal_solve_controls(), feedback_controls(),
                &mut workspace, &gate, resolve,
            ).unwrap();
        assert_eq!(result, original);
        assert!(result.evolution.evolution.evolution.viscosity.cell_sgs_transfer_j.is_empty());
    }
}

#[test]
fn cooled_molecular_caps_whole_call_work_and_late_cancellation_preserve_retry() {
    let (material, melting) = thermal_material_with_rate(true);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let conductivity = thermal_conductivity(&material, caloric);
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 3];
    let phases = [Some(caloric.state_at_specific_enthalpy(1.0).unwrap()); 3];
    let momentum = [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0]];
    let selected = sgs(SgsModel::Smagorinsky { cs: 0.3 });
    let calls = Cell::new(0);
    let changed_temperature = Cell::new(false);
    #[derive(Clone, Copy, PartialEq)]
    enum Fault {
        None,
        StalePhase,
        Cancel,
    }
    let run = |selected, viscosity, phases: &[Option<EquilibriumPhaseState>], authored: &[f64], fault,
               workspace: &mut ImplicitThermalVolumeWorkspace, gate: &CancelGate| {
        let original_h = phases[0].unwrap().specific_enthalpy_j_kg();
        field.advance_isotropic_implicit_conjugate_thermorheological_thermal_volume_momentum_with_baths_and_sgs(
            &domain, None, &conductivity, ThermalVolumeBoundaries::default(), &[], &[],
            ThermalVolumeMomentumInput {
                volume: input(&mass, phases, authored), momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3], force_density_n_m3: &[[0.0; 3]; 3],
                boundary_velocity_m_s: [None; 6],
            },
            [FvBoundary::Symmetry; 6], selected, coupled_controls(0.5), viscosity,
            thermal_solve_controls(), feedback_controls(), workspace, gate,
            |_, phase, rate| {
                calls.set(calls.get() + 1);
                let mut h = phase.specific_enthalpy_j_kg();
                if (h - original_h).abs() > 1e-10 {
                    changed_temperature.set(true);
                    if fault == Fault::Cancel {
                        gate.request();
                        return Err(refuse("source refusal after requested cancellation"));
                    }
                    if fault == Fault::StalePhase {
                        h += 0.001;
                    }
                }
                joint_response(&material, &law, h, rate)
            },
        )
    };
    let mut workspace = ImplicitThermalVolumeWorkspace::new(3).unwrap();
    let gate = CancelGate::new_clock_free();
    let expected = run(
        selected, molecular_controls(), &phases, &[0.0; 3], Fault::None, &mut workspace, &gate,
    ).unwrap();
    let total = expected.evolution.feedback.molecular_sweeps;
    assert!(total > expected.evolution.evolution.evolution.viscosity.sweeps + 1);
    assert!(run(
        selected, MaterialViscosityControls { max_sweeps: total - 1, ..molecular_controls() },
        &phases, &[0.0; 3], Fault::None, &mut workspace, &gate,
    ).is_err());
    assert!(run(
        selected, molecular_controls(), &phases, &[0.0; 3], Fault::StalePhase,
        &mut workspace, &gate,
    ).is_err());
    let cancelled = CancelGate::new_clock_free();
    calls.set(0);
    assert!(matches!(run(
        selected, molecular_controls(), &phases, &[0.0; 3], Fault::Cancel,
        &mut workspace, &cancelled,
    ), Err(ChtError::Cancelled)));
    assert!(calls.get() > 3);

    // Manufacture a colder supported endpoint h=[1.01,1.02,1.01] from
    // original h=1.07. The native momentum root fixes the molecular heat;
    // authored heat closes the original three caloric rows independently.
    let cold = [1.01, 1.02, 1.01];
    let original_cold = [Some(caloric.state_at_specific_enthalpy(1.07).unwrap()); 3];
    let e_sgs = (8.0 / 3.0) * 0.3_f64.powi(2) * 2.0_f64.sqrt();
    let quadratic = 1.25 * (17.0 / 24.0 + e_sgs);
    let linear = 1.375 + (31.0 / 3.0) * (0.01 / 4.0 + 0.02 / 2.0);
    let b = 1.0 / (linear + (linear * linear + 2.0 * quadratic).sqrt());
    let cooling: [f64; 3] = core::array::from_fn(|cell| {
        let g = if cell == 1 { 2.0 * b } else { b };
        let lambda = 0.5 + (31.0 / 3.0) * (cold[cell] - 1.0) + (17.0 / 24.0) * g;
        let conduction: f64 = (0..3).filter(|&other| cell.abs_diff(other) == 1)
            .map(|other| cold[other] - cold[cell]).sum();
        2.0 * (cold[cell] - 1.07) - 0.5 * lambda * g * g - conduction
    });
    let cooled = run(
        selected, molecular_controls(), &original_cold, &cooling, Fault::None,
        &mut workspace, &gate,
    ).unwrap();
    for (cell, expected_h) in cold.into_iter().enumerate() {
        near(
            cooled.evolution.endpoint_sources[cell].as_ref().unwrap()
                .phase_state().specific_enthalpy_j_kg(),
            expected_h,
        );
    }
    // Initial maximum nu_t/nu=.38863884, cooled endpoint=.40364420.
    // A cap at .395 must therefore reject only after refreshed thermal mu.
    let initial_ratio = 2.0 * 0.3_f64.powi(2) * 2.0_f64.sqrt() / 0.655;
    let final_ratio = 2.0 * 0.3_f64.powi(2) * 2.0_f64.sqrt() * 2.0 * b / (0.205 + b / 2.0);
    assert!(initial_ratio < 0.395 && final_ratio > 0.395);
    calls.set(0);
    changed_temperature.set(false);
    assert!(matches!(run(
        FvSgsConfig { maximum_eddy_to_molecular: 0.395, ..selected },
        molecular_controls(), &original_cold, &cooling, Fault::None, &mut workspace, &gate,
    ), Err(ChtError::InvalidInput { field: "simple.sgs", .. })));
    assert!(changed_temperature.get());
    assert!(calls.get() > 3);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 1.0);
    assert_eq!(original_cold[0].unwrap().specific_enthalpy_j_kg(), 1.07);
    assert_eq!(momentum, [[2.0, 0.0, 0.0], [0.0; 3], [-2.0, 0.0, 0.0]]);
    calls.set(0);
    let retry = run(
        selected, molecular_controls(), &phases, &[0.0; 3], Fault::None, &mut workspace, &gate,
    ).unwrap();
    assert_eq!(retry, expected);
    assert_eq!(calls.get(), retry.evolution.feedback.source_evaluations);
}
