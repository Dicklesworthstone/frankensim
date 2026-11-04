//! Source continuation in the actual midpoint trajectory, with an independent
//! scalar Maxwell/Kelvin solution on the unit tetrahedron's node-1 x DOF.

use super::*;

const MASS: f64 = 1000.0 / 60.0;
const H: f64 = 0.02;

fn dynamic_config(steps: usize) -> ViscoDynamicConfig {
    ViscoDynamicConfig {
        maximum_steps: steps,
        maximum_free_dofs: 1,
        solve: solve_config(),
        absolute_energy_tolerance_j: 1e-14,
        relative_energy_tolerance: 2e-10,
    }
}

#[derive(Clone, Copy)]
struct ScalarMotion {
    u: f64,
    v: f64,
    q: f64,
    free: f64,
}

impl ScalarMotion {
    fn storage(self, e: f64, branch: f64) -> f64 {
        0.5 * VOLUME * (e * (self.u - self.free).powi(2) + branch * self.q * self.q)
    }
    // q' + q/tau = (u-free)'; exact ramp integrals give the physical work
    // mean. Combining that expression with midpoint momentum leaves one scalar
    // linear equation, without any native element/material/solver calls.
    fn advance(&mut self, coefficients: [f64; 4], force: f64, free: f64) -> [f64; 4] {
        let [e, branch, tau, viscosity] = coefficients;
        let decay = (-H / tau).exp();
        let phi = tau / H * (1.0 - decay);
        let tangent = 0.5 * e + branch * tau / H * (1.0 - phi) + viscosity / H;
        let delta_free = free - self.free;
        let old_strain = self.u - self.free;
        let hold_mean = e * old_strain + branch * phi * self.q - tangent * delta_free;
        let w = (MASS * self.v + 0.5 * H * (force - VOLUME * hold_mean))
            / (MASS + 0.5 * H * H * VOLUME * tangent);
        let increment = H * w;
        let delta = increment - delta_free;
        let mean_stress = e * (old_strain + 0.5 * delta)
            + branch * (phi * self.q + tau / H * (1.0 - phi) * delta)
            + viscosity * delta / H;
        let stationary = tau * delta / H;
        let transient = self.q - stationary;
        let heat = VOLUME
            * (branch / tau
                * (stationary * stationary * H
                    + 2.0 * stationary * transient * tau * (1.0 - decay)
                    + 0.5 * transient * transient * tau * (1.0 - decay * decay))
                + viscosity * delta * delta / H);
        self.u += increment;
        self.v = 2.0 * w - self.v;
        self.q = decay * self.q + phi * delta;
        self.free = free;
        let stress = e * (self.u - free) + branch * self.q + viscosity * delta / H;
        [
            VOLUME * stress,
            heat,
            -VOLUME * mean_stress * delta_free,
            force * increment,
        ]
    }
}

#[test]
fn live_source_change_preserves_motion_memory_clock_budget_and_parameter_work() {
    with_cx(|cx, gate| {
        let card = card(true);
        let source = source(&card);
        let original = tet(&source, &[3], cx);
        let warm = original
            .at_query_points(&[&card], &[point(325.0)], cx)
            .unwrap();
        let next = warm.at_query_points(&[&card], &[point(350.0)], cx).unwrap();
        let configuration = dynamic_config(6);
        let old_operator = warm.dynamic_operator(H, configuration, cx).unwrap();
        let next_operator = next.dynamic_operator(H, configuration, cx).unwrap();
        close(warm.equilibrium().mass.get(0, 0), MASS);
        assert_eq!(
            warm.equilibrium().mass.row(0),
            next.equilibrium().mass.row(0)
        );
        assert!(
            (old_operator.effective_mass().get(0, 0) - next_operator.effective_mass().get(0, 0))
                .abs()
                > 1e-4
        );

        let mut history = original.state(cx).unwrap();
        let mut displacement = [0.0; 12];
        displacement[3] = 0.002;
        history.step_internal_force(&displacement, H, cx).unwrap();
        let initial_change = history.rebind_materials(&warm, cx).unwrap();
        assert!(initial_change.parameter_work_j < 0.0);
        let preloaded_parameter_work = initial_change.history.parameter_work_j();
        let preloaded_heat = initial_change.history.energy().dissipated_j;
        let mut velocity = [0.0; 12];
        velocity[3] = 0.01;
        let mut motion = old_operator
            .start(initial_change.history, &velocity, cx)
            .unwrap();
        assert_eq!(motion.parameter_work_j(), 0.0);
        assert_eq!(
            motion.material_history().parameter_work_j(),
            preloaded_parameter_work
        );
        let initial_total = motion.kinetic_energy_j() + motion.material_history().energy().stored_j;
        let mut oracle = ScalarMotion {
            u: displacement[3],
            v: velocity[3],
            q: phi(H, 0.07) * displacement[3],
            free: 0.0,
        };
        let mut external = 0.0;
        let mut heat = 0.0;
        for load in [0.03, -0.01, 0.02] {
            let expected = oracle.advance([900.0, 450.0, 0.105, 6.0], load, 0.0);
            let mut force = [0.0; 12];
            force[3] = load;
            let step = motion.advance(&force, cx).unwrap();
            close(step.constitutive.full_force_n[3], expected[0]);
            close(step.constitutive.step_dissipated_j, expected[1]);
            close(motion.material_history().displacement_m()[3], oracle.u);
            close(motion.velocity_mps()[3], oracle.v);
            external += expected[3];
            heat += expected[1];
        }
        let before = motion.clone();
        let expected_work = 0.5
            * VOLUME
            * ((600.0 - 900.0) * oracle.u * oracle.u + (300.0 - 450.0) * oracle.q * oracle.q);
        let changed = motion.rebind_materials(&next_operator, cx).unwrap();
        close(changed.parameter_work_j, expected_work);
        let mut continued = changed.trajectory;
        close(continued.parameter_work_j(), expected_work);
        close(
            continued.material_history().parameter_work_j(),
            preloaded_parameter_work + expected_work,
        );
        assert_eq!(continued.time_s(), before.time_s());
        assert_eq!(continued.steps(), 3);
        assert_eq!(continued.velocity_mps(), before.velocity_mps());
        assert_eq!(continued.kinetic_energy_j(), before.kinetic_energy_j());
        assert_eq!(
            continued.material_history().displacement_m(),
            before.material_history().displacement_m()
        );
        assert_eq!(continued.external_work_j(), before.external_work_j());
        assert_eq!(continued.eigenstrain_work_j(), before.eigenstrain_work_j());
        assert_eq!(
            continued.material_history().energy().work_j,
            before.material_history().energy().work_j
        );
        assert_eq!(
            continued.material_history().energy().dissipated_j,
            before.material_history().energy().dissipated_j
        );
        close(
            continued.material_history().energy().stored_j,
            oracle.storage(600.0, 300.0),
        );

        for load in [-0.015, 0.0, 0.025] {
            let expected = oracle.advance([600.0, 300.0, 0.14, 4.0], load, 0.0);
            let mut force = [0.0; 12];
            force[3] = load;
            let step = continued.advance(&force, cx).unwrap();
            close(step.constitutive.full_force_n[3], expected[0]);
            close(step.constitutive.step_dissipated_j, expected[1]);
            close(continued.material_history().displacement_m()[3], oracle.u);
            close(continued.velocity_mps()[3], oracle.v);
            external += expected[3];
            heat += expected[1];
            close(continued.external_work_j(), external);
            close(
                continued.material_history().energy().dissipated_j - preloaded_heat,
                heat,
            );
            close(
                continued.kinetic_energy_j()
                    + continued.material_history().energy().stored_j
                    + heat,
                initial_total + external + expected_work,
            );
        }
        assert_eq!(continued.steps(), configuration.maximum_steps);
        close(continued.time_s(), 6.0 * H);
        assert_eq!(continued.eigenstrain_work_j(), 0.0);
        let final_energy = continued.material_history().energy();
        assert!(continued.advance(&[0.0; 12], cx).is_err());
        assert_eq!(continued.material_history().energy(), final_energy);
        assert_eq!(continued.steps(), 6);

        // A rebuild cannot restart the original clock or enlarge any controls.
        let enlarged = next
            .dynamic_operator(
                H,
                ViscoDynamicConfig {
                    maximum_steps: 7,
                    ..configuration
                },
                cx,
            )
            .unwrap();
        assert!(continued.rebind_materials(&enlarged, cx).is_err());
        let changed_duration = next.dynamic_operator(2.0 * H, configuration, cx).unwrap();
        assert!(continued.rebind_materials(&changed_duration, cx).is_err());
        let foreign = tet(&source, &[3], cx);
        let foreign_operator = foreign.dynamic_operator(H, configuration, cx).unwrap();
        assert!(continued.rebind_materials(&foreign_operator, cx).is_err());
        let tight = warm.at_query_points(&[&card], &[point(400.0)], cx).unwrap();
        let tight_operator = tight.dynamic_operator(H, configuration, cx).unwrap();
        assert!(continued.rebind_materials(&tight_operator, cx).is_err());
        // Same-source no-op refresh remains exact even after consuming all time
        // intervals, and still cannot create a seventh accepted interval.
        let unchanged = continued.rebind_materials(&next_operator, cx).unwrap();
        assert_eq!(unchanged.parameter_work_j, 0.0);
        let mut terminal = unchanged.trajectory;
        assert_eq!(terminal.time_s(), continued.time_s());
        assert_eq!(terminal.steps(), 6);
        assert!(terminal.advance(&[0.0; 12], cx).is_err());
        gate.request();
        assert!(matches!(
            continued.rebind_materials(&next_operator, cx),
            Err(ViscoSourceError::Native(ViscoElementError::Elastic(
                TetElasticError::Cancelled
            )))
        ));
        assert_eq!(continued.material_history().energy(), final_energy);
        assert_eq!(continued.velocity_mps(), terminal.velocity_mps());
        assert_eq!(continued.external_work_j(), terminal.external_work_j());
        assert_eq!(continued.steps(), 6);
    });
}

fn scalar_thermal(solid: &ViscoElements, t: f64) -> TetThermalStrainState {
    TetThermalStrainState::try_new(
        t,
        300.0,
        [1e-5 * (t - 300.0), 0.0, 0.0, 0.0, 0.0, 0.0],
        solid.reference_material(0).unwrap().material_state_identity,
        ContentHash([74; 32]),
        0.03,
    )
    .unwrap()
}

#[test]
fn thermal_source_refresh_keeps_existing_work_and_drives_the_same_native_trajectory() {
    with_cx(|cx, _| {
        let card = card(true);
        let source = source(&card);
        let original = tet(&source, &[3], cx);
        let next = original
            .at_query_points(&[&card], &[point(350.0)], cx)
            .unwrap();
        let configuration = dynamic_config(3);
        let initial_operator = original.dynamic_operator(H, configuration, cx).unwrap();
        let next_operator = next.dynamic_operator(H, configuration, cx).unwrap();
        let mut motion = initial_operator
            .start(original.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        let mut oracle = ScalarMotion {
            u: 0.0,
            v: 0.0,
            q: 0.0,
            free: 0.0,
        };
        let original_thermal = scalar_thermal(&original, 350.0);
        let replacement = scalar_thermal(&next, 350.0);
        let mut force = [0.0; 12];
        force[3] = 0.02;
        let first = oracle.advance(
            [1200.0, 600.0, 0.07, 8.0],
            force[3],
            original_thermal.free_strain_mandel()[0],
        );
        let response = motion
            .advance_with_thermal_strain(
                &force,
                fs_solid::linear3::TetThermalStrainField::Uniform(&original_thermal),
                cx,
            )
            .unwrap();
        close(response.constitutive.full_force_n[3], first[0]);
        close(response.constitutive.step_eigenstrain_work_j, first[2]);
        assert!(motion.eigenstrain_work_j().abs() > 1e-6);
        assert!(motion.rebind_materials(&next_operator, cx).is_err());
        let wrong_thermal = scalar_thermal(&next, 351.0);
        assert!(
            motion
                .rebind_materials_with_thermal_strain(
                    &next_operator,
                    fs_solid::linear3::TetThermalStrainField::Uniform(&wrong_thermal),
                    cx
                )
                .is_err()
        );
        let old_energy = motion.material_history().energy();
        let before = motion.clone();
        let change = motion
            .rebind_materials_with_thermal_strain(
                &next_operator,
                fs_solid::linear3::TetThermalStrainField::Uniform(&replacement),
                cx,
            )
            .unwrap();
        let work = -0.5
            * VOLUME
            * (600.0 * (oracle.u - oracle.free).powi(2) + 300.0 * oracle.q * oracle.q);
        close(change.parameter_work_j, work);
        let mut continued = change.trajectory;
        assert_eq!(continued.time_s(), before.time_s());
        assert_eq!(continued.steps(), 1);
        assert_eq!(continued.velocity_mps(), before.velocity_mps());
        assert_eq!(continued.kinetic_energy_j(), before.kinetic_energy_j());
        assert_eq!(continued.external_work_j(), before.external_work_j());
        assert_eq!(continued.eigenstrain_work_j(), before.eigenstrain_work_j());
        assert_eq!(
            continued.material_history().energy().dissipated_j,
            old_energy.dissipated_j
        );
        assert_eq!(
            continued.material_history().thermal_strains().unwrap(),
            &[replacement]
        );
        let next_thermal = scalar_thermal(&next, 360.0);
        let second = oracle.advance(
            [600.0, 300.0, 0.14, 4.0],
            force[3],
            next_thermal.free_strain_mandel()[0],
        );
        let response = continued
            .advance_with_thermal_strain(
                &force,
                fs_solid::linear3::TetThermalStrainField::Uniform(&next_thermal),
                cx,
            )
            .unwrap();
        close(response.constitutive.full_force_n[3], second[0]);
        close(response.constitutive.step_eigenstrain_work_j, second[2]);
        close(continued.material_history().displacement_m()[3], oracle.u);
        close(continued.velocity_mps()[3], oracle.v);
        let third = oracle.advance([600.0, 300.0, 0.14, 4.0], force[3], oracle.free);
        let response = continued.advance(&force, cx).unwrap();
        assert_eq!(response.constitutive.step_eigenstrain_work_j, 0.0);
        close(response.constitutive.full_force_n[3], third[0]);
        close(continued.material_history().displacement_m()[3], oracle.u);
        close(continued.velocity_mps()[3], oracle.v);
        close(continued.external_work_j(), first[3] + second[3] + third[3]);
        close(continued.eigenstrain_work_j(), first[2] + second[2]);
        close(
            continued.material_history().energy().dissipated_j,
            first[1] + second[1] + third[1],
        );
        close(continued.parameter_work_j(), work);
        close(
            continued.kinetic_energy_j()
                + continued.material_history().energy().stored_j
                + continued.material_history().energy().dissipated_j,
            continued.external_work_j() + continued.eigenstrain_work_j() + work,
        );
        assert_eq!(
            continued.material_history().thermal_strains().unwrap(),
            &[next_thermal]
        );
        assert_eq!(continued.steps(), 3);
    });
}
