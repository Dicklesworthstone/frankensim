//! Accepted heated equilibrium drives the original consistent-mass trajectory.
use super::*;

use fs_solid::cutfront::dynamics3::CutSecondOrderProblem3;
use fs_solver::{FlexiblePreconditioner, NewtonKrylovConfig};
use fs_time::galpha::initialization::{InitialSolveConfig, second_order_acceleration};
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderProblem, SecondOrderState,
};

struct Identity;
impl FlexiblePreconditioner for Identity {
    fn apply(&self, _: usize, residual: &[f64], output: &mut [f64]) {
        output.copy_from_slice(residual);
    }
}

#[test]
fn accepted_heated_state_releases_load_through_original_mass_and_signed_tangent() {
    let law = law(0.05);
    let cold = frozen(&law, 300.0);
    let hot = |bounds: HexCell| Ok(frozen(
        &law, 325.0 + 12.0 * bounds.lo()[0] + 8.0 * bounds.lo()[2],
    ));
    let body = |_| [0.0, -0.02, 0.01];
    for adapted in [false, true] {
        // Dynamics holds homogeneous supports; embedded motion has its own
        // support-work requirements and is intentionally absent from this case.
        let mut model = build(adapted, &cold, false);
        let native_n = model.operator().n();
        let cells = native!(&model, cells);
        let scales = vec![0.7; cells];
        native!(&mut model, set_scales, &scales).unwrap();
        let occupancy = vec![1.0; cells];
        let cold_mass = native!(&model, prepare_material_mass, POINTS, proceed).unwrap();
        let original_inertia = cold_mass.diagonal(&occupancy, proceed).unwrap();
        let controls = settings(native_n);
        let accepted = native!(&mut model, refresh_thermoelastic_equilibrium,
            ElasticMaterialField3::ThermoelasticByCell(&hot),
            ReferenceLoad3::body(&body), None, None, controls, |_| proceed(),
        ).unwrap();
        assert_eq!(accepted.mass.diagonal(&occupancy, proceed).unwrap(), original_inertia);

        let gate = CancelGate::new_clock_free();
        let trajectory = CutSecondOrderProblem3::new_with_eigenstrain(
            model.source(), &accepted.mass, &occupancy, &accepted.eigenstrain,
            controls.maximum_kinematic_strain, POINTS, native_n, &gate,
        ).unwrap().with_reference_prestress(
            accepted.solution.coefficients(), POINTS, POINTS,
        ).unwrap();
        let reference = trajectory.restrict_force(accepted.solution.coefficients()).unwrap();
        let n = trajectory.dimension();
        let zero = vec![0.0; n];
        let external = native!(&model, body_load, &body, proceed).unwrap();
        let held_force = trajectory.restrict_force(&external).unwrap();
        let mut resistance = zero.clone();
        trajectory.internal_force(&reference, &mut resistance);
        same(&resistance, &held_force, 3e-8);

        // The trajectory's inseparable reference produces the same signed G as
        // the accepted equilibrium packet, including the adaptive reduction.
        let geometric = accepted.prestress.as_ref().unwrap();
        assert!(geometric.scalar_pairs().any(|(_, value)| value.abs() > 1e-10));
        let direction: Vec<_> = native!(&model, nodes).iter()
            .flat_map(|position| [0.0, position[0], 0.0]).collect();
        let mut from_packet = vec![0.0; native_n];
        let mut from_trajectory = from_packet.clone();
        geometric.apply_into(&direction, &mut from_packet, proceed).unwrap();
        trajectory.reference_prestress().unwrap()
            .apply_into(&direction, &mut from_trajectory, proceed).unwrap();
        assert_eq!(from_packet, from_trajectory);

        // Release the physical external load. Thermal free strain remains held
        // inside native resistance; adding its force again would fail this gate.
        let initial = second_order_acceleration(
            &trajectory, &reference, &zero, &zero, &Identity,
            InitialSolveConfig { restart: n, max_cycles: 8, tolerance: 1e-12 },
            2_000_000, &mut || gate.is_requested(),
        ).unwrap();
        let acceleration_work: f64 = held_force.iter().zip(&initial.value)
            .map(|(force, acceleration)| force * acceleration).sum();
        assert!(acceleration_work < -1e-10);
        let native_acceleration = trajectory.expand_field(&initial.value).unwrap();
        let mut inertia_force = vec![0.0; native_n];
        accepted.mass.apply_into(
            &occupancy, &native_acceleration, &mut inertia_force, proceed,
        ).unwrap();
        for &dof in trajectory.free_dofs() {
            close(inertia_force[dof], -external[dof], 3e-8);
        }

        let dt = 0.001;
        let method = OperatorGeneralizedAlpha::new(n, dt, 1.0, ImplicitSolveConfig {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-10,
                relative_tolerance: 1e-10,
                linear_restart: n,
                max_linear_cycles: 8,
                forcing_maximum: 1e-3,
                ..Default::default()
            },
            max_newton_iterations: 8,
        });
        let mut state = SecondOrderState::new(0.0, &reference, &zero, &initial.value);
        trajectory.advance(&method, &mut state, &vec![0.0; native_n]).unwrap();
        assert_eq!(state.steps, 1);
        close(state.t, dt, 1e-15);
        assert!(state.v.iter().any(|value| value.abs() > 1e-8));
        assert!(state.q.iter().zip(&reference)
            .any(|(value, reference)| (value - reference).abs() > 1e-12));

        // At midpoint parameters and a consistent initial acceleration, the
        // linear model also satisfies the endpoint balance. Evaluate it using
        // only the original returned mass, native K/free load and accepted G.
        let displacement = trajectory.expand_field(&state.q).unwrap();
        let acceleration = trajectory.expand_field(&state.a).unwrap();
        accepted.mass.apply_into(
            &occupancy, &acceleration, &mut inertia_force, proceed,
        ).unwrap();
        let elastic = action(&model, &displacement);
        let free = native!(&model, eigenstrain_load, &accepted.eigenstrain, proceed).unwrap();
        let increment: Vec<_> = displacement.iter().zip(accepted.solution.coefficients())
            .map(|(value, reference)| value - reference).collect();
        geometric.apply_into(&increment, &mut from_packet, proceed).unwrap();
        for &dof in trajectory.free_dofs() {
            close(inertia_force[dof] + elastic[dof] - free[dof] + from_packet[dof], 0.0, 4e-8);
        }
    }
}
