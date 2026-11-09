//! Nonlinear original material sources and SGS on one native viscous endpoint.

use super::*;
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

fn coefficients(rheology: Rheology) -> [f64; 4] {
    match rheology {
        Rheology::Thickening => [0.5, 0.5, 0.0, 0.0],
        Rheology::Thinning => [0.875, -0.25, 0.0, 0.0],
        Rheology::ShearAndBulk => [0.25, 0.25, 0.25, 0.5],
    }
}

/// Native faces [0,a,0], rho=2, dt=1: gamma=2a/sqrt(3), mu=m0+m1*a,
/// zeta=b0+b1*a. The declared Smagorinsky invariant is sqrt(2)*a, so
/// mu_t=2*Cs^2*sqrt(2)*a. Dividing M(a-1)+2*(4(mu+mu_t)/3+zeta)*a=0
/// by the original face mass M=2 gives A*a^2+B*a-1=0.
fn endpoint_speed(rheology: Rheology, cs: f64) -> f64 {
    let [m0, m1, b0, b1] = coefficients(rheology);
    let a = (4.0 / 3.0) * m1 + b1 + (8.0 / 3.0) * cs * cs * 2.0_f64.sqrt();
    let b = 1.0 + (4.0 / 3.0) * m0 + b0;
    2.0 / (b + (b * b + 4.0 * a).sqrt())
}

#[test]
fn original_nonlinear_sources_and_sgs_close_on_one_physical_endpoint() {
    for rheology in [Rheology::Thinning, Rheology::Thickening, Rheology::ShearAndBulk] {
        let (material, melting) = fixture(rheology);
        let caloric = melting.curve();
        let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
        let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
        let field = FlowField::quiescent(&domain);
        let phase = caloric.state_at_specific_enthalpy(1.0).unwrap();
        let mass = [2.0; 2];
        let phases = [Some(phase); 2];
        let momentum = [[2.0, 0.0, 0.0]; 2];
        let authored_heat = [0.125, -0.0625];
        let cs = 0.35;
        let mut queries = vec![];
        let result = field
            .advance_piecewise_thermal_volume_momentum_with_endpoint_rheology_and_sgs(
                &domain,
                None,
                caloric,
                input(&mass, &phases, &momentum, &authored_heat, &[[0.0; 3]; 2]),
                [FvBoundary::Symmetry; 6],
                sgs(SgsModel::Smagorinsky { cs }),
                controls(),
                viscosity_controls(),
                &CancelGate::new_clock_free(),
                |cell, supplied_phase, rate| {
                    assert_eq!(supplied_phase, &phase);
                    queries.push((cell, rate));
                    resolve(
                        &material, &law, supplied_phase.specific_enthalpy_j_kg(), rate, rheology,
                    )
                },
            )
            .unwrap();
        let a = endpoint_speed(rheology, cs);
        assert!((a - 0.5).abs() > 0.02); // Original no-SGS root for every source.
        let [m0, m1, b0, b1] = coefficients(rheology);
        let mu = m0 + m1 * a;
        let bulk = b0 + b1 * a;
        let mu_t = 2.0 * cs * cs * 2.0_f64.sqrt() * a;
        let rate = 2.0 * a / 3.0_f64.sqrt();
        let heat = ((4.0 / 3.0) * mu + bulk) * a * a;
        let transfer = (4.0 / 3.0) * mu_t * a * a;
        let viscous = &result.evolution.viscosity;
        let volume = &result.evolution.evolution.thermal.volume;
        assert!(viscous.sweeps > 1);
        assert_eq!(result.source_updates, viscous.sweeps + 1);
        assert_eq!(result.source_evaluations, 2 * result.source_updates);
        assert_eq!(queries.len(), result.source_evaluations);
        for pair in queries.chunks_exact(2) {
            assert_eq!((pair[0].0, pair[1].0), (0, 1));
            close(pair[0].1, pair[1].1);
        }
        for cell in 0..2 {
            // Neither molecular weights nor SGS coefficients may contaminate
            // unit-quadrature rate observation, initially or after each solve.
            close(queries[cell].1, 2.0 / 3.0_f64.sqrt());
            close(queries[queries.len() - 2 + cell].1, rate);
            close(result.endpoint_shear_rate_per_s[cell], rate);
            let source = result.endpoint_sources[cell].as_ref().unwrap();
            assert_eq!(source.phase_state(), &phase);
            close(source.transport_state().unwrap().dynamic_viscosity_pa_s(), mu);
            close(source.material().property(NU_KEY).unwrap().value_si(), mu / 2.0);
            for property in source.material().properties() {
                let coordinates = &property.answer().receipt.query_point;
                close(
                    coordinates.iter().find(|(axis, _)| axis == "rate").unwrap().1,
                    rate,
                );
                assert!(coordinates.contains(&("T".into(), 301.0)));
            }
            close(result.bulk_viscosity_pa_s[cell], bulk);
            close(viscous.cell_molecular_heat_j[cell], heat);
            close(viscous.cell_sgs_transfer_j[cell], transfer);
            // The original caloric path receives authored heat plus ONLY the
            // physical molecular loss, not the separate modeled transfer.
            let expected_h = 1.0 + (authored_heat[cell] + heat) / mass[cell];
            close(
                volume.phase_states[cell].unwrap().specific_enthalpy_j_kg(),
                expected_h,
            );
            close(volume.transport.mass_kg[cell], mass[cell]);
        }
        close(viscous.face_impulse_kg_m_s[0], 2.0 * (a - 1.0));
        close(viscous.native_kinetic_energy_change_j, a * a - 1.0);
        close(viscous.backward_euler_loss_j, (a - 1.0).powi(2));
        close(viscous.molecular_heat_j, 2.0 * heat);
        close(viscous.sgs_transfer_j, 2.0 * transfer);
        assert!(viscous.sgs_transfer_j > 0.05);
        close(
            viscous.native_kinetic_energy_change_j + 2.0 * (heat + transfer)
                + viscous.backward_euler_loss_j - viscous.wall_work_j,
            0.0,
        );
        close(
            viscous.cell_kinetic_energy_change_j + 2.0 * (heat + transfer)
                + viscous.backward_euler_loss_j + viscous.numerical_scatter_loss_j
                - viscous.wall_work_j - viscous.representation_work_j,
            viscous.energy_defect_j,
        );
        close(viscous.energy_defect_j, 0.0);
        close(
            volume.report.material_enthalpy_change_j,
            authored_heat.iter().sum::<f64>() + 2.0 * heat,
        );
        assert!(viscous.relative_residual <= viscosity_controls().relative_tolerance);
        assert_eq!(field, FlowField::quiescent(&domain));
        assert_eq!(phases, [Some(phase); 2]);
        assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
    }
}

#[test]
fn zero_closures_are_the_exact_original_source_and_caloric_step() {
    let (material, melting) = fixture(Rheology::ShearAndBulk);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let phase = caloric.state_at_specific_enthalpy(1.0).unwrap();
    let mass = [2.0; 2];
    let phases = [Some(phase); 2];
    let momentum = [[2.0, 0.0, 0.0]; 2];
    let gate = CancelGate::new_clock_free();
    let source = |_: usize, phase: &EquilibriumPhaseState, rate| {
        resolve(
            &material, &law, phase.specific_enthalpy_j_kg(), rate, Rheology::ShearAndBulk,
        )
    };
    let original = field
        .advance_piecewise_thermal_volume_momentum_with_endpoint_rheology(
            &domain,
            None,
            caloric,
            input(&mass, &phases, &momentum, &[0.0; 2], &[[0.0; 3]; 2]),
            [FvBoundary::Symmetry; 6],
            controls(),
            viscosity_controls(),
            &gate,
            source,
        )
        .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let result = field
            .advance_piecewise_thermal_volume_momentum_with_endpoint_rheology_and_sgs(
                &domain,
                None,
                caloric,
                input(&mass, &phases, &momentum, &[0.0; 2], &[[0.0; 3]; 2]),
                [FvBoundary::Symmetry; 6],
                FvSgsConfig { maximum_eddy_to_molecular: 0.0, ..sgs(model) },
                controls(),
                viscosity_controls(),
                &gate,
                source,
            )
            .unwrap();
        assert_eq!(result, original);
        assert!(result.evolution.viscosity.cell_sgs_transfer_j.is_empty());
        assert_eq!(result.evolution.viscosity.sgs_transfer_j, 0.0);
    }
}

#[test]
fn source_cap_work_and_cancel_refusals_leave_one_original_retry() {
    let (material, melting) = fixture(Rheology::Thickening);
    let caloric = melting.curve();
    let law = MushyPhaseLaw::try_from_melting_curve(&melting, 1e-8).unwrap();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let phase = caloric.state_at_specific_enthalpy(1.0).unwrap();
    let mass = [2.0; 2];
    let phases = [Some(phase); 2];
    let momentum = [[2.0, 0.0, 0.0]; 2];
    let selected = sgs(SgsModel::Smagorinsky { cs: 0.35 });
    #[derive(Clone, Copy, PartialEq)]
    enum Fault {
        None,
        Source,
        StalePhase,
        Cancel,
    }
    let run = |selected, viscosity, fault, calls: &Cell<usize>, gate: &CancelGate| {
        field.advance_piecewise_thermal_volume_momentum_with_endpoint_rheology_and_sgs(
            &domain,
            None,
            caloric,
            input(&mass, &phases, &momentum, &[0.0; 2], &[[0.0; 3]; 2]),
            [FvBoundary::Symmetry; 6],
            selected,
            controls(),
            viscosity,
            gate,
            |_, supplied_phase, rate| {
                calls.set(calls.get() + 1);
                assert_eq!(supplied_phase, &phase);
                if calls.get() > 2 && matches!(fault, Fault::Source | Fault::Cancel) {
                    if fault == Fault::Cancel {
                        gate.request();
                    }
                    return Err(ChtError::InvalidInput {
                        field: "test.endpoint_sgs_source",
                        reason: "late source refusal".into(),
                    });
                }
                let h = supplied_phase.specific_enthalpy_j_kg()
                    + if calls.get() > 2 && fault == Fault::StalePhase { 0.001 } else { 0.0 };
                resolve(&material, &law, h, rate, Rheology::Thickening)
            },
        )
    };
    let gate = CancelGate::new_clock_free();
    let calls = Cell::new(0);
    let expected = run(selected, viscosity_controls(), Fault::None, &calls, &gate).unwrap();
    for fault in [Fault::Source, Fault::StalePhase] {
        calls.set(0);
        assert!(run(selected, viscosity_controls(), fault, &calls, &gate).is_err());
        assert_eq!(calls.get(), 3);
    }
    calls.set(0);
    assert!(matches!(
        run(
            FvSgsConfig { maximum_eddy_to_molecular: 0.0, ..selected },
            viscosity_controls(), Fault::None, &calls, &gate,
        ),
        Err(ChtError::InvalidInput { field: "simple.sgs", .. })
    ));
    assert_eq!(calls.get(), 2); // Original molecular field precedes the SGS cap.
    calls.set(0);
    assert!(matches!(
        run(
            selected,
            MaterialViscosityControls { max_sweeps: 1, ..viscosity_controls() },
            Fault::None, &calls, &gate,
        ),
        Err(ChtError::SolverNotConverged { .. })
    ));
    assert_eq!(calls.get(), 4);
    for late in [false, true] {
        let cancelled = CancelGate::new_clock_free();
        if !late {
            cancelled.request();
        }
        calls.set(0);
        assert!(matches!(
            run(selected, viscosity_controls(), Fault::Cancel, &calls, &cancelled),
            Err(ChtError::Cancelled)
        ));
        assert_eq!(calls.get(), if late { 3 } else { 0 });
    }
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(phases, [Some(phase); 2]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]; 2]);
    calls.set(0);
    assert_eq!(
        run(selected, viscosity_controls(), Fault::None, &calls, &gate).unwrap(),
        expected,
    );
    assert_eq!(calls.get(), expected.source_evaluations);
}
