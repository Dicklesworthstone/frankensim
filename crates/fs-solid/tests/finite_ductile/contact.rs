//! Native rigid tooling with analytical compression and physical-rate oracles.
use super::*;
use fs_contact::sustained::solve::SolveControls;
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_solid::hyper3::plastic::ductile::contact::{
    FiniteDuctilePlaneContactError, PlaneContactControls, RigidPlaneTool, TetPlaneContact,
};
use fs_solid::hyper3::plastic::ductile::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileTetMotion,
};
use fs_solver::NewtonKrylovConfig;

const V: f64 = 1.0 / 6.0;

fn contact_controls() -> PlaneContactControls {
    PlaneContactControls {
        solve: SolveControls {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-7,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            difference_step_n: 0.01,
            force_residual_tolerance_n: 1e-6,
            ..SolveControls::default()
        },
        maximum_normal_work_defect_j: 2e-6,
        work_balance_tolerance_j: 2e-6,
        ..PlaneContactControls::default()
    }
}
fn material_controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 2e-7,
        ..controls()
    }
}
fn moving_controls() -> FiniteDuctileMotionControls {
    FiniteDuctileMotionControls {
        material: material_controls(),
        maximum_numerical_kinetic_loss_j: 1.0,
    }
}
fn tool(axis: usize, before: f64, after: f64) -> RigidPlaneTool {
    let mut normal = [0.0; 3];
    let mut point_before_m = [0.0; 3];
    let mut translation_m = [0.0; 3];
    normal[axis] = -1.0;
    point_before_m[axis] = before;
    translation_m[axis] = after - before;
    RigidPlaneTool {
        normal,
        point_before_m,
        translation_m,
    }
}
fn row(initial_force_n: f64) -> TetPlaneContact {
    TetPlaneContact {
        node: 1,
        plane: 0,
        initial_force_n,
        scale_n_per_m: 1e4,
    }
}
fn supports(axis: usize) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 3 + axis)
        .map(|dof| (dof, 0.0))
        .collect()
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}

#[test]
fn compressive_tools_keep_bulk_resistance_progressive_damage_rotation_and_release() {
    with_cx(|_, cx| {
        for axis in [0, 1] {
            let nodes = if axis == 0 {
                NODES
            } else {
                NODES.map(|v| [v[2], v[0], v[1]])
            };
            let geometry =
                TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
            let mut normal = [0.0; 3];
            normal[axis] = 1.0;
            let bands = [TetCrackBandGeometry::from_element(&geometry, 0, normal).unwrap()];
            let material = model(bands[0].band_width_m(), 1.0);
            let materials = [&material];
            let solid = FiniteDuctileTetSolid::try_new(
                &geometry,
                &materials,
                &bands,
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let a = -0.03_f64;
            let end = a.exp();
            let alpha = (1000.0 * -a - 3.0) / 1500.0;
            let integrity = (-bands[0].band_width_m() * 3.0 * (alpha - ONSET) / GP).exp();
            // Uniaxial strain, compressive volume retained: tau_11=K*a-2*g.
            let expected_piola = (1000.0 * a - 2.0 * integrity) / end;
            let formed = solid
                .advance_with_plane_contact(
                    &supports(axis),
                    &[],
                    [0.0; 3],
                    conditions(),
                    &[tool(axis, 1.0, end)],
                    &[row(0.0)],
                    IncrementSettings::default(),
                    material_controls(),
                    contact_controls(),
                    cx,
                )
                .unwrap();
            close(formed.step.solid.displacement_m()[1][axis], end - 1.0, 3e-9);
            close(
                formed.step.solid.element_states()[0].native_state().alpha,
                alpha,
                3e-9,
            );
            close(
                formed.step.solid.element_states()[0].integrity(),
                integrity,
                3e-8,
            );
            close(
                formed.contact.solve.forces_n[0][0],
                -V * expected_piola,
                2e-5,
            );
            close(
                formed.contact.force_on_tool_n[0][axis],
                -V * expected_piola,
                2e-5,
            );
            close(
                formed.contact.body_work_j,
                V * expected_piola * (end - 1.0),
                2e-6,
            );
            close(formed.contact.tool_work_j, formed.contact.body_work_j, 2e-6);
            close(formed.contact.work_balance_error_j, 0.0, 2e-6);
            assert!(formed.step.energy.plastic_dissipation_j > 0.0);
            assert_eq!(formed.step.energy.damage_release_j, 0.0);
            assert!(formed.contact.minimum_gap_m >= 0.0);
            assert!(formed.contact.path_windows > 0);
            assert_eq!(solid.element_states()[0].damage(), 0.0);
            let held = formed
                .step
                .solid
                .advance_with_plane_contact(
                    &supports(axis),
                    &[],
                    [0.0; 3],
                    conditions(),
                    &[tool(axis, end, end)],
                    &[row(formed.contact.solve.forces_n[0][0])],
                    IncrementSettings::default(),
                    material_controls(),
                    contact_controls(),
                    cx,
                )
                .unwrap();
            close(
                held.step.solid.element_states()[0].native_state().alpha,
                alpha,
                3e-9,
            );
            close(
                held.step.solid.element_states()[0].integrity(),
                integrity,
                3e-8,
            );
            close(held.step.energy.plastic_dissipation_j, 0.0, 2e-8);
            // Release is compared with the exact original free endpoint. Lateral
            // constraints can cause reverse yielding, which must remain native.
            let free = held
                .step
                .solid
                .advance(
                    &supports(axis),
                    &[],
                    [0.0; 3],
                    conditions(),
                    IncrementSettings::default(),
                    material_controls(),
                    cx,
                )
                .unwrap();
            let released = held
                .step
                .solid
                .advance_with_plane_contact(
                    &supports(axis),
                    &[],
                    [0.0; 3],
                    conditions(),
                    &[tool(axis, end, 1.05)],
                    &[row(0.0)],
                    IncrementSettings::default(),
                    material_controls(),
                    contact_controls(),
                    cx,
                )
                .unwrap();
            close(released.contact.solve.forces_n[0][0], 0.0, 1e-7);
            assert!(released.contact.solve.paths_m[0][0] > 0.04);
            close(
                released.step.solid.displacement_m()[1][axis],
                free.solid.displacement_m()[1][axis],
                3e-9,
            );
            assert!(released.step.solid.displacement_m()[1][axis] < -1e-5);
            close(
                released.step.solid.element_states()[0].damage(),
                free.solid.element_states()[0].damage(),
                3e-9,
            );
        }
    });
}

#[test]
fn moving_rate_dependent_tool_keeps_contact_kinetic_and_damage_work_distinct() {
    with_cx(|_, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let viscosity = 150.0;
        let material = model(bands[0].band_width_m(), 1.0)
            .with_overstress(FiniteOverstressParameters {
                viscosity_pa_s: viscosity,
                stress_scale_pa: 3.0,
                exponent: 1.0,
            })
            .unwrap();
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let motion = FiniteDuctileTetMotion::try_new(solid, &[12.0], &[[0.0; 3]; 4], cx).unwrap();
        let before = motion.solid().element_states().to_vec();
        let a = -0.03_f64;
        let end = a.exp();
        let dt = conditions().duration_s;
        let alpha = (1000.0 * -a - 3.0) / (1500.0 + viscosity / dt);
        let q = 3.0 + viscosity * alpha / dt;
        let integrity = (-bands[0].band_width_m() * 3.0 * (alpha - ONSET) / GP).exp();
        let piola = (1000.0 * a - 2.0 * q * integrity / 3.0) / end;
        let u = end - 1.0;
        let velocity = u / dt;
        // Exact P1 consistent diagonal rho*V/10=.2 kg at the one moving DOF.
        let inertia_force = 0.2 * u / (dt * dt);
        let force = 0.5 - V * piola - inertia_force;
        let tools = [tool(0, 1.01, end)];
        let mut strict = contact_controls();
        strict.maximum_normal_work_defect_j = 0.0;
        let refused = motion.advance_with_plane_contact(
            &supports(0),
            &[(3, 0.25)],
            [6.0, 0.0, 0.0],
            conditions(),
            &tools,
            &[row(0.0)],
            IncrementSettings::default(),
            moving_controls(),
            strict,
            cx,
        );
        assert!(matches!(
            refused,
            Err(FiniteDuctilePlaneContactError::Energy {
                quantity: "absolute normal endpoint work defect",
                ..
            })
        ));
        assert_eq!(motion.solid().element_states(), before);
        assert_eq!(motion.velocity_mps(), &[[0.0; 3]; 4]);
        let mut admitted = contact_controls();
        admitted.maximum_normal_work_defect_j = 1.0;
        let formed = motion
            .advance_with_plane_contact(
                &supports(0),
                &[(3, 0.25)],
                [6.0, 0.0, 0.0],
                conditions(),
                &tools,
                &[row(0.0)],
                IncrementSettings::default(),
                moving_controls(),
                admitted,
                cx,
            )
            .unwrap();
        close(
            formed.step.motion.solid().element_states()[0]
                .native_state()
                .alpha,
            alpha,
            3e-9,
        );
        close(
            formed.step.motion.solid().element_states()[0].integrity(),
            integrity,
            3e-8,
        );
        close(formed.contact.solve.forces_n[0][0], force, 2e-5);
        close(formed.contact.normal_work_defect_j, force * 0.01, 2e-6);
        close(formed.contact.dead_load_work_j, 0.5 * u, 2e-8);
        close(
            formed.step.inertia.kinetic_after_j,
            0.1 * velocity * velocity,
            2e-8,
        );
        close(
            formed.step.inertia.numerical_loss_j,
            formed.step.inertia.kinetic_after_j,
            2e-8,
        );
        close(formed.contact.work_balance_error_j, 0.0, 2e-6);
        close(formed.contact.load_decomposition_error_j, 0.0, 2e-6);
        assert!(formed.step.material_energy.plastic_dissipation_j > 0.0);
        assert_eq!(formed.step.material_energy.damage_release_j, 0.0);
        let held = formed
            .step
            .motion
            .advance_with_plane_contact(
                &supports(0),
                &[],
                [0.0; 3],
                conditions(),
                &[tool(0, end, end)],
                &[row(formed.contact.solve.forces_n[0][0])],
                IncrementSettings::default(),
                moving_controls(),
                contact_controls(),
                cx,
            )
            .unwrap();
        assert!(held.step.material_energy.damage_release_j > 0.0);
        let alpha2 = alpha + (q - 3.0) / (1500.0 + viscosity / dt);
        close(
            held.step.motion.solid().element_states()[0]
                .native_state()
                .alpha,
            alpha2,
            3e-9,
        );
        assert!(
            held.step.motion.solid().element_states()[0].damage()
                > formed.step.motion.solid().element_states()[0].damage()
        );
        close(held.step.inertia.kinetic_after_j, 0.0, 2e-8);
        close(
            held.step.inertia.kinetic_change_j + held.step.inertia.numerical_loss_j,
            0.0,
            2e-8,
        );
        close(held.step.motion.reference_mass_kg(), 2.0, 1e-12);
        close(held.contact.work_balance_error_j, 0.0, 2e-6);
    });
}

#[test]
fn omitted_penetrating_vertices_damage_limits_and_cancellation_publish_no_endpoint() {
    with_cx(|gate, cx| {
        let geometry =
            TetLoadGeometry3::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let bands = [TetCrackBandGeometry::from_element(&geometry, 0, [1.0, 0.0, 0.0]).unwrap()];
        let material = model(bands[0].band_width_m(), 1e-4);
        let materials = [&material];
        let solid = FiniteDuctileTetSolid::try_new(
            &geometry,
            &materials,
            &bands,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let omitted = solid.advance_with_plane_contact(
            &supports(1),
            &[],
            [0.0; 3],
            conditions(),
            &[tool(1, 1.0, 0.9)],
            &[row(0.0)],
            IncrementSettings::default(),
            material_controls(),
            contact_controls(),
            cx,
        );
        assert!(matches!(
            omitted,
            Err(FiniteDuctilePlaneContactError::Penetration {
                node: 2,
                at_end: true,
                ..
            })
        ));
        let before = solid.element_states().to_vec();
        assert!(
            solid
                .advance_with_plane_contact(
                    &supports(0),
                    &[],
                    [0.0; 3],
                    conditions(),
                    &[tool(0, 1.0, (-0.03_f64).exp())],
                    &[row(0.0)],
                    IncrementSettings::default(),
                    material_controls(),
                    contact_controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(solid.element_states(), before);
        assert_eq!(solid.displacement_m(), &[[0.0; 3]; 4]);
        gate.request();
        assert!(matches!(
            solid.advance_with_plane_contact(
                &supports(0),
                &[],
                [0.0; 3],
                conditions(),
                &[tool(0, 1.0, 1.0)],
                &[row(0.0)],
                IncrementSettings::default(),
                material_controls(),
                contact_controls(),
                cx,
            ),
            Err(FiniteDuctilePlaneContactError::Cancelled)
        ));
        assert_eq!(solid.element_states(), before);
    });
}
