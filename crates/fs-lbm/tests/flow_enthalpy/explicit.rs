//! Actual solved momentum with native conservative second-order h transport.
use super::*;
use fs_conduction::transient::deposition::advection::reconstruction::{
    AxialCellStencil, LimitedEnthalpyReconstruction,
};
use fs_conduction::transient::deposition::backward_euler::{
    BoundaryStepCoefficients, StepCoefficients,
};
use fs_lbm::conjugate::enthalpy::advection::explicit::{
    AdaptiveSsprk2Controls, AdaptiveSsprk2Error, Ssprk2Controls,
};

fn controls(dt: f64) -> Ssprk2Controls {
    Ssprk2Controls {
        step: StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: 1e-8,
        },
        maximum_stage_evaluations: 100_000,
    }
}

fn adaptive(dt: f64) -> AdaptiveSsprk2Controls {
    AdaptiveSsprk2Controls {
        step: controls(dt),
        initial_dt_s: dt,
        minimum_dt_s: 1e-9,
        maximum_trials: 10_000,
        enthalpy_absolute_tolerance_j_kg: 1e-4,
        temperature_absolute_tolerance_k: Some(1e-9),
        relative_tolerance: 0.0,
    }
}

#[test]
fn solved_mass_flux_has_second_order_latent_throughflow_accuracy() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(2.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let gate = CancelGate::new();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let law = curve(2.0, 100.0);
    let old = initial_profile(&domain, 2.0, &law, |_| 20.0);
    let ids = [Some(CellId(0)), Some(CellId(1)), Some(CellId(2))];
    let mut binding = input(&ids, &law);
    binding.inlet_specific_h_j_kg = [Some(100.0); 6];
    let owner = solved.advect_enthalpy(&old, binding, (), &gate).unwrap();
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let dt = 4.0 / count as f64;
        let (next, accepted, report) = owner
            .advance_ssprk2_to(&old, 4.0, controls(dt / (1.0 - dt * 0.25 / 0.8)), &gate)
            .unwrap();
        let mut squared = 0.0;
        // Exact semidiscrete upwind chain: h_j=100-80 exp(-z) sum(z^k/k!).
        // Here z=(rho*u/m)*t=1, independently of the sensible Cp.
        for (i, series) in [1.0, 2.0, 2.5].into_iter().enumerate() {
            let exact = 100.0 - 80.0 * (-1.0_f64).exp() * series;
            squared += (next.cells()[&CellId(i as u64)].specific_h_j_kg - exact).powi(2);
            near(next.temperature_k(CellId(i as u64)).unwrap(), 310.0);
        }
        errors.push(squared.sqrt());
        near(accepted.advection().audit().inlet_mass_kg, 2.0);
        near(accepted.advection().audit().outlet_mass_kg, 2.0);
        near(
            report.energy.external_heat_j,
            accepted.advection().audit().external_enthalpy_j,
        );
        near(report.energy.residual_j, 0.0);
        assert!(core::ptr::eq(accepted.material(), &source));
        assert!(core::ptr::eq(accepted.field(), &solved.flow().field));
        assert_eq!(accepted.transport(), solved.transport());
        assert_eq!(accepted.advection().flows(), owner.advection().flows());
        assert_eq!(report.stage_evaluations, 2 * report.energy.substeps);
    }
    for pair in errors.windows(2) {
        assert!(
            (3.5..4.6).contains(&(pair[0] / pair[1])),
            "errors {errors:?}"
        );
    }
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(owner.advection().audit().substeps, 0);
}

#[test]
fn adaptive_solved_flow_retains_reconstruction_and_fine_only_latent_ledgers() {
    let n = 24;
    let domain = VoxelDomain::new(n, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let law = curve(1.0, 100.0);
    let profile = |c| if (9..15).contains(&c) { 120.0 } else { 20.0 };
    let old = initial_profile(&domain, 1.0, &law, profile);
    let ids: Vec<_> = (0..n).map(|i| Some(CellId(i as u64))).collect();
    let stencils: Vec<_> = (1..n - 1)
        .map(|i| AxialCellStencil {
            previous: CellId(i as u64 - 1),
            cell: CellId(i as u64),
            next: CellId(i as u64 + 1),
        })
        .collect();
    let plan = LimitedEnthalpyReconstruction::try_new(&old, &stencils, n, || Ok(())).unwrap();
    let gate = CancelGate::new();
    for forward in [true, false] {
        let mut config = throughflow();
        if !forward {
            config.faces[0] = FvBoundary::Outlet;
            config.faces[1] = FvBoundary::Inlet {
                velocity: [-0.25, 0.0, 0.0],
            };
        }
        let solved = fluid.solve_flow(&domain, &config, &gate).unwrap();
        let mut binding = input(&ids, &law);
        binding.inlet_specific_h_j_kg = [Some(20.0); 6];
        let donor = solved.advect_enthalpy(&old, binding, (), &gate).unwrap();
        let limited = donor.clone().with_reconstruction(plan.clone()).unwrap();
        let mut errors = [0.0; 2];
        for (scheme, owner) in [&donor, &limited].into_iter().enumerate() {
            let (next, accepted, report) = owner
                .advance_ssprk2_adaptive_to(&old, 12.0, adaptive(2.0), &gate)
                .unwrap();
            for i in 0..n {
                let upstream = if forward {
                    i.checked_sub(3)
                } else {
                    Some(i + 3)
                };
                let exact = upstream.map_or(20.0, profile);
                let h = next.cells()[&CellId(i as u64)].specific_h_j_kg;
                assert!((20.0 - 1e-10..=120.0 + 1e-10).contains(&h));
                near(next.temperature_k(CellId(i as u64)).unwrap(), 310.0);
                errors[scheme] += (h - exact).abs();
            }
            if scheme == 1 {
                assert!(core::ptr::eq(
                    accepted.reconstruction().unwrap(),
                    limited.reconstruction().unwrap()
                ));
            }
            assert!(report.accuracy_rejections > 0);
            assert!(report.maximum_accepted_normalized_error <= 1.0);
            assert_eq!(report.work.energy.substeps, 2 * report.accepted_intervals);
            assert_eq!(
                accepted.advection().audit().substeps,
                report.work.energy.substeps
            );
            assert_eq!(accepted.advection().flows(), owner.advection().flows());
            near(accepted.advection().audit().inlet_mass_kg, 3.0);
            near(accepted.advection().audit().outlet_mass_kg, 3.0);
            near(
                report.work.energy.external_heat_j,
                accepted.advection().audit().external_enthalpy_j,
            );
            near(report.work.energy.residual_j, 0.0);
        }
        assert!(
            errors[1] < 0.9 * errors[0],
            "forward={forward}, errors {errors:?}"
        );
    }
}

#[test]
fn stage_source_heat_is_separate_from_solved_stream_enthalpy() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(2.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let gate = CancelGate::new();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let law = curve(2.0, 100.0);
    let old = initial_profile(&domain, 2.0, &law, |_| 20.0);
    let ids = [Some(CellId(0)), Some(CellId(1)), Some(CellId(2))];
    let mut binding = input(&ids, &law);
    binding.inlet_specific_h_j_kg = [Some(80.0); 6];
    let owner = solved.advect_enthalpy(&old, binding, (), &gate).unwrap();
    let source = |state: &DepositionSimulation<'_>| -> Result<BoundaryStepCoefficients<()>, DepositionError> {
        Ok(StepCoefficients { sources_w: BTreeMap::from([(CellId(1), 2.0 * state.time_s())]),
            conductances_w_k: BTreeMap::new(), payload: () }.into())
    };
    let (_, fixed, fixed_report) = owner
        .advance_ssprk2_with_coefficients(&old, 2.0, controls(0.125), source, &gate)
        .unwrap();
    near(
        fixed_report.energy.external_heat_j - fixed.advection().audit().external_enthalpy_j,
        4.0,
    );
    let (next, accepted, report) = owner
        .advance_ssprk2_adaptive_with_coefficients(&old, 2.0, adaptive(2.0), source, &gate)
        .unwrap();
    near(
        report.work.energy.external_heat_j - accepted.advection().audit().external_enthalpy_j,
        4.0,
    );
    near(accepted.advection().audit().inlet_mass_kg, 1.0);
    near(accepted.advection().audit().outlet_mass_kg, 1.0);
    near(
        next.energy_j() - old.energy_j(),
        report.work.energy.external_heat_j,
    );
    near(report.work.energy.residual_j, 0.0);
}

#[test]
fn changed_map_law_budget_and_cancellation_publish_no_flow_candidate() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(2.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let gate = CancelGate::new();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let law = curve(2.0, 100.0);
    let old = initial_profile(&domain, 2.0, &law, |_| 20.0);
    let ids = [Some(CellId(0)), Some(CellId(1)), Some(CellId(2))];
    let owner = solved
        .advect_enthalpy(&old, input(&ids, &law), (), &gate)
        .unwrap();
    let wrong_mass = initial_profile(&domain, 1.0, &law, |_| 20.0);
    assert!(matches!(
        owner.advance_ssprk2_to(&wrong_mass, 0.0, controls(1.0), &gate),
        Err(DepositionError::InvalidInput(_))
    ));
    let foreign = curve(2.0, 100.0);
    let wrong_law = initial_profile(&domain, 2.0, &foreign, |_| 20.0);
    assert!(matches!(
        owner.advance_ssprk2_adaptive_to(&wrong_law, 0.0, adaptive(1.0), &gate),
        Err(AdaptiveSsprk2Error::Deposition(
            DepositionError::InvalidInput(_)
        ))
    ));
    let mut bounded = adaptive(1.0);
    bounded.step.maximum_stage_evaluations = 4;
    assert!(matches!(
        owner.advance_ssprk2_adaptive_to(&old, 1.0, bounded, &gate),
        Err(AdaptiveSsprk2Error::Deposition(DepositionError::Budget(_)))
    ));
    let mut queries = 0;
    assert!(matches!(
        owner.advance_ssprk2_adaptive_with_coefficients(
            &old,
            1.0,
            adaptive(1.0),
            |_| {
                queries += 1;
                if queries == 5 {
                    gate.request();
                }
                Ok(StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                }
                .into())
            },
            &gate
        ),
        Err(AdaptiveSsprk2Error::Deposition(DepositionError::Cancelled))
    ));
    assert_eq!(queries, 5);
    assert_eq!(owner.advection().audit().substeps, 0);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, 20.0);
}
