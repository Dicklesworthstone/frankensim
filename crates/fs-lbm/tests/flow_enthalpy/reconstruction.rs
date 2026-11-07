//! The material-flow bridge uses the original h reconstruction and CFL owner.

use super::*;
use fs_conduction::transient::deposition::advection::reconstruction::{
    AxialCellStencil, LimitedEnthalpyReconstruction,
};

fn plan(thermal: &DepositionSimulation<'_>, n: usize) -> LimitedEnthalpyReconstruction {
    let stencils: Vec<_> = (1..n - 1)
        .map(|c| AxialCellStencil {
            previous: CellId(c as u64 - 1),
            cell: CellId(c as u64),
            next: CellId(c as u64 + 1),
        })
        .collect();
    LimitedEnthalpyReconstruction::try_new(thermal, &stencils, n, || Ok(())).unwrap()
}

#[test]
fn solved_flow_reconstructs_latent_fronts_conservatively_in_both_directions() {
    // Passive equal-density latent inventory is a transported caloric tracer.
    // This does not give the stationary solid fraction a mobility law.
    let n = 24;
    let domain = VoxelDomain::new(n, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let law = curve(1.0, 100.0);
    let initial_h = |c| if (9..15).contains(&c) { 120.0 } else { 20.0 };
    let old = initial_profile(&domain, 1.0, &law, initial_h);
    let ids: Vec<_> = (0..n).map(|c| Some(CellId(c as u64))).collect();
    let gate = CancelGate::new();
    let mut binding = input(&ids, &law);
    binding.inlet_specific_h_j_kg = [Some(20.0); 6];
    let dt = 0.125;
    let duration = 12.0; // At 0.25 m/s the exact profile moves three cells.
    let decay = 0.25;
    // Both schemes use the SAME physical Euler timestep after their native
    // combined limits; the error comparison isolates spatial reconstruction.
    let controls = |multiplier| StepConfig {
        max_dt_s: dt / (1.0 - multiplier * dt * decay / 0.8),
        energy_tolerance_j: 1e-8,
    };
    for forward in [true, false] {
        let mut flow_config = throughflow();
        if !forward {
            flow_config.faces[0] = FvBoundary::Outlet;
            flow_config.faces[1] = FvBoundary::Inlet {
                velocity: [-0.25, 0.0, 0.0],
            };
        }
        let solved = fluid.solve_flow(&domain, &flow_config, &gate).unwrap();
        let donor = solved.advect_enthalpy(&old, binding, (), &gate).unwrap();
        let reconstructed = donor.clone().with_reconstruction(plan(&old, n)).unwrap();
        assert!(donor.reconstruction().is_none());
        assert_eq!(reconstructed.reconstruction().unwrap().len(), n - 2);
        let (low, low_flow, _) = old
            .advance_observed(&donor, duration, controls(1.0), || Ok(()))
            .unwrap();
        let (high, high_flow, energy) = old
            .advance_observed(&reconstructed, duration, controls(2.0), || Ok(()))
            .unwrap();
        near(low_flow.advection().audit().maximum_courant, dt * decay);
        near(high_flow.advection().audit().maximum_courant, dt * decay);
        let mut low_error = 0.0;
        let mut high_error = 0.0;
        for c in 0..n {
            let upstream = if forward {
                c.checked_sub(3)
            } else {
                Some(c + 3)
            };
            let exact = upstream.map_or(20.0, initial_h);
            let low_h = low.cells()[&CellId(c as u64)].specific_h_j_kg;
            let high_h = high.cells()[&CellId(c as u64)].specific_h_j_kg;
            assert!((20.0 - 1e-10..=120.0 + 1e-10).contains(&high_h));
            near(high.temperature_k(CellId(c as u64)).unwrap(), 310.0);
            low_error += (low_h - exact).abs();
            high_error += (high_h - exact).abs();
        }
        assert!(
            high_error < 0.9 * low_error,
            "forward={forward}: reconstructed={high_error}, donor={low_error}"
        );
        near(high.energy_j() - old.energy_j(), energy.external_heat_j);
        near(
            energy.external_heat_j,
            high_flow.advection().audit().external_enthalpy_j,
        );
        near(high.mass_kg(), old.mass_kg());
        assert_eq!(high_flow.advection().flows(), donor.advection().flows());
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(reconstructed.advection().audit().substeps, 0);
    }
}

#[test]
fn reconstruction_preserves_empty_plan_arithmetic_and_failed_candidate_history() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let source = FluidTransportState::try_new(1.0, 2.0, 0.01, 0.1).unwrap();
    let fluid = MaterialChtFluid::try_new(&source).unwrap();
    let law = curve(1.0, 100.0);
    let old = initial(&domain, 1.0, &law);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let solved = fluid.solve_flow(&domain, &throughflow(), &gate).unwrap();
    let donor = solved
        .advect_enthalpy(&old, input(&ids, &law), Heater::default(), &gate)
        .unwrap();
    let empty = LimitedEnthalpyReconstruction::try_new(&old, &[], 0, || Ok(())).unwrap();
    let wrapped = donor.clone().with_reconstruction(empty).unwrap();
    let (expected, expected_flow, expected_energy) = old
        .advance_observed(&donor, 1.0, config(), || Ok(()))
        .unwrap();
    let (actual, actual_flow, actual_energy) = old
        .advance_observed(&wrapped, 1.0, config(), || Ok(()))
        .unwrap();
    assert_eq!(actual.cells(), expected.cells());
    assert_eq!(actual_energy, expected_energy);
    assert_eq!(
        actual_flow.advection().audit(),
        expected_flow.advection().audit()
    );
    assert_eq!(
        actual_flow.advection().inner(),
        expected_flow.advection().inner()
    );

    let history = Heater {
        reject_step: Some(2),
        ..Heater::default()
    };
    let limited = solved
        .advect_enthalpy(&old, input(&ids, &law), history.clone(), &gate)
        .unwrap()
        .with_reconstruction(plan(&old, 4))
        .unwrap();
    assert!(
        old.advance_observed(&limited, 1.0, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(limited.advection().inner(), &history);
    assert_eq!(limited.advection().audit().substeps, 0);
    assert_eq!(old.time_s(), 0.0);
    assert_eq!(old.energy_j(), 0.0);
}
