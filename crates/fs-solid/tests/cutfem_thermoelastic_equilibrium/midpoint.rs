//! Changing native source endpoints retain actual inertia and midpoint work.
use super::*;

use fs_solid::cutfront::dynamics3::{CutSecondOrderProblem3, CutThermalMidpointSettings3};
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig, norm2};
use fs_time::galpha::initialization::{InitialSolveConfig, second_order_acceleration};
use fs_time::galpha::{ImplicitSolveConfig, SecondOrderProblem, SecondOrderState, TimeSolveError};

struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}

fn step_settings(n: usize) -> CutThermalMidpointSettings3 {
    CutThermalMidpointSettings3 {
        solve: ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-11,
                linear_restart: n,
                max_linear_cycles: 8,
                forcing_maximum: 1e-3,
                ..Default::default()
            },
            max_newton_iterations: 8,
        },
        maximum_steps: 4,
        maximum_workspace_components: 2_000_000,
        maximum_force_residual_n: 2e-8,
    }
}

#[test]
fn changing_temperature_uses_both_native_force_endpoints_and_original_mass() {
    let law = law(0.05);
    let cold = frozen(&law, 300.0);
    let hot = |bounds: HexCell| Ok(frozen(
        &law, 325.0 + 12.0 * bounds.lo()[0] + 8.0 * bounds.lo()[2],
    ));
    let body0 = |_| [0.01, -0.02, 0.01];
    let body1 = |_| [-0.015, 0.01, -0.005];
    for adapted in [false, true] {
        let mut owner = build(adapted, &cold, false);
        let native_n = owner.operator().n();
        let cells = native!(&owner, cells);
        let scales: Vec<_> = (0..cells).map(|i| 0.6 + 0.01 * (i % 9) as f64).collect();
        native!(&mut owner, set_scales, &scales).unwrap();
        let occupancy = vec![0.73; cells];
        let old_mass = native!(&owner, prepare_material_mass, POINTS, proceed).unwrap();
        let old_strain = native!(&owner, prepare_thermoelastic, POINTS, proceed).unwrap();
        let original_inertia = old_mass.diagonal(&occupancy, proceed).unwrap();
        let old_force = native!(&owner, body_load, &body0, proceed).unwrap();
        let next_force = native!(&owner, body_load, &body1, proceed).unwrap();
        let native_q: Vec<_> = native!(&owner, nodes).iter()
            .flat_map(|p| [1e-4*p[0], -2e-4*p[0], 0.5e-4*p[0]]).collect();
        let native_v: Vec<_> = native!(&owner, nodes).iter()
            .flat_map(|p| [2e-4*p[0], 1e-4*p[0], -3e-4*p[0]]).collect();
        let gate = CancelGate::new_clock_free();
        let (before, old_resistance, reduced_old_resistance, old_inertia) = {
            let original = CutSecondOrderProblem3::new_with_eigenstrain(
                owner.source(), &old_mass, &occupancy, &old_strain,
                0.08, POINTS, native_n, &gate,
            ).unwrap();
            let q = original.restrict_force(&native_q).unwrap();
            let v = original.restrict_force(&native_v).unwrap();
            let external = original.restrict_force(&old_force).unwrap();
            let initial = second_order_acceleration(
                &original, &q, &v, &external, &Identity,
                InitialSolveConfig { restart: original.dimension(), max_cycles: 8, tolerance: 1e-12 },
                2_000_000, &mut || gate.is_requested(),
            ).unwrap();
            let mut resistance = vec![0.0; original.dimension()];
            original.internal_force(&q, &mut resistance);
            let full = original.expand_field(&resistance).unwrap();
            let mut inertia = resistance.clone();
            original.mass_apply(&initial.value, &mut inertia);
            (SecondOrderState::new(0.25, &q, &v, &initial.value), full, resistance, inertia)
        };
        let unchanged = before.clone();
        native!(&mut owner, refresh_material_field,
            ElasticMaterialField3::ThermoelasticByCell(&hot), Default::default(), proceed,
        ).unwrap();
        let mass = native!(&owner, prepare_material_mass, POINTS, proceed).unwrap();
        let strain = native!(&owner, prepare_thermoelastic, POINTS, proceed).unwrap();
        assert_eq!(mass.diagonal(&occupancy, proceed).unwrap(), original_inertia);
        let current = CutSecondOrderProblem3::new_with_eigenstrain(
            owner.source(), &mass, &occupancy, &strain,
            0.08, POINTS, native_n, &gate,
        ).unwrap();
        let controls = step_settings(current.dimension());
        let h = 0.001;
        let candidate = current.thermal_midpoint_candidate(
            &before, &old_resistance, &old_force, &next_force, h, controls,
        ).unwrap();
        assert_eq!(before, unchanged);
        assert_eq!(candidate.state.steps, 1);
        assert_eq!(candidate.state.history.len(), 1);
        assert_eq!(candidate.state.history[0], candidate.telemetry);
        close(candidate.state.t, before.t + h, 1e-15);
        assert!(candidate.initial_residual_norm_n <= controls.maximum_force_residual_n);
        assert!(candidate.endpoint_residual_norm_n <= controls.maximum_force_residual_n);

        // Evaluate the endpoint through the original independent M, K and
        // source-owned free force, not through the private centered wrapper.
        let q1 = current.expand_field(&candidate.state.q).unwrap();
        let a1 = current.expand_field(&candidate.state.a).unwrap();
        let elastic = action(&owner, &q1);
        let free = native!(&owner, eigenstrain_load, &strain, proceed).unwrap();
        let mut inertia = vec![0.0; native_n];
        mass.apply_into(&occupancy, &a1, &mut inertia, proceed).unwrap();
        let f0 = current.restrict_force(&old_force).unwrap();
        let f1 = current.restrict_force(&next_force).unwrap();
        let mut midpoint_residual = Vec::new();
        let mut endpoint_residual = Vec::new();
        for (i, &dof) in current.free_dofs().iter().enumerate() {
            let resistance1 = elastic[dof] - free[dof];
            close(candidate.native_resistance_n[dof], resistance1, 2e-12);
            endpoint_residual.push(inertia[dof] + resistance1 - f1[i]);
            midpoint_residual.push(0.5 * (old_inertia[i] + inertia[dof]
                + reduced_old_resistance[i] + resistance1 - f0[i] - f1[i]));
            close(candidate.state.q[i] - before.q[i],
                0.5*h*(before.v[i] + candidate.state.v[i]), 2e-12);
            close(candidate.state.v[i] - before.v[i],
                0.5*h*(before.a[i] + candidate.state.a[i]), 2e-12);
        }
        assert!(norm2(&endpoint_residual) <= controls.maximum_force_residual_n);
        assert!(norm2(&midpoint_residual) <= controls.maximum_force_residual_n);

        // The common frozen-new-source shortcut omits half the actual thermal
        // endpoint force change. This case has enough C(T)/free-strain change
        // to make that defect observable, independent of the solver tolerance.
        let midpoint: Vec<_> = before.q.iter().zip(&candidate.state.q)
            .map(|(old, new)| 0.5*(old+new)).collect();
        let mut frozen_new = vec![0.0; current.dimension()];
        current.internal_force(&midpoint, &mut frozen_new);
        let next_resistance = current.restrict_force(&candidate.native_resistance_n).unwrap();
        let omitted: Vec<_> = frozen_new.iter().zip(&reduced_old_resistance)
            .zip(&next_resistance)
            .map(|((shortcut, old), new)| shortcut - 0.5*(old+new)).collect();
        assert!(norm2(&omitted) > 1e-4);

        // Workspace/step/time and inconsistent-acceleration refusals leave the
        // same accepted checkpoint untouched, even after a successful trial.
        let workspace = CutThermalMidpointSettings3 { maximum_workspace_components: 0, ..controls };
        assert!(current.thermal_midpoint_candidate(
            &before, &old_resistance, &old_force, &next_force, h, workspace,
        ).is_err());
        let exhausted = CutThermalMidpointSettings3 { maximum_steps: 0, ..controls };
        assert!(current.thermal_midpoint_candidate(
            &before, &old_resistance, &old_force, &next_force, h, exhausted,
        ).is_err());
        assert!(matches!(current.thermal_midpoint_candidate(
            &before, &old_resistance, &old_force, &next_force, 0.0, controls,
        ), Err(TimeSolveError::InvalidTimeAdvance { .. })));
        let mut inconsistent = before.clone();
        inconsistent.a.fill(0.0);
        assert!(current.thermal_midpoint_candidate(
            &inconsistent, &old_resistance, &old_force, &next_force, h, controls,
        ).is_err());
        assert_eq!(before, unchanged);
    }
}
