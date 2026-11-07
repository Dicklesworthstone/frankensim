//! Independent dynamic cavity/mass/work oracles on the original admitted fixtures.

use super::*;
use fs_solid::hyper3::plastic::porous::dynamics::{
    FinitePlasticMotionControls, FinitePorousCavitiesMotionStep, FinitePorousTetMotion,
};

const RHO: f64 = 6.0;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

/// Independent reference volume and linear-tet basis gradients.
fn reference(nodes: &[[f64; 3]], tet: [usize; 4]) -> (f64, [[f64; 3]; 4]) {
    let a = sub(nodes[tet[1]], nodes[tet[0]]);
    let b = sub(nodes[tet[2]], nodes[tet[0]]);
    let c = sub(nodes[tet[3]], nodes[tet[0]]);
    let determinant = dot(a, cross(b, c));
    assert!(determinant != 0.0);
    let g1 = cross(b, c).map(|x| x / determinant);
    let g2 = cross(c, a).map(|x| x / determinant);
    let g3 = cross(a, b).map(|x| x / determinant);
    let g0 = std::array::from_fn(|i| -g1[i] - g2[i] - g3[i]);
    (determinant.abs() / 6.0, [g0, g1, g2, g3])
}

fn motion<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a FinitePorous,
    velocity: &[[f64; 3]],
    cx: &Cx<'_>,
) -> FinitePorousTetMotion<'a> {
    FinitePorousTetMotion::try_new(
        solid(nodes, tets, material, cx),
        &vec![RHO; tets.len()],
        velocity,
        cx,
    )
    .unwrap()
}

fn motion_controls() -> FinitePlasticMotionControls {
    FinitePlasticMotionControls {
        material: controls(),
        maximum_numerical_kinetic_loss_j: 20.0,
    }
}

fn settings() -> HyperTetSettings {
    HyperTetSettings {
        relative_tolerance: 1e-10,
        max_newton_iterations: 80,
        ..HyperTetSettings::default()
    }
}

/// Exact constant-reference-density P1 mass action, M_ab=rho*V*(1+delta_ab)/20.
fn mass_action(
    nodes: &[[f64; 3]],
    tets: &[[usize; 4]],
    acceleration: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut force = vec![[0.0; 3]; nodes.len()];
    for &tet in tets {
        let (volume, _) = reference(nodes, tet);
        let sum: [f64; 3] = std::array::from_fn(|i| tet.iter().map(|&n| acceleration[n][i]).sum());
        for node in tet {
            for i in 0..3 {
                force[node][i] += RHO * volume / 20.0 * (sum[i] + acceleration[node][i]);
            }
        }
    }
    force
}

fn kinetic(nodes: &[[f64; 3]], tets: &[[usize; 4]], velocity: &[[f64; 3]]) -> f64 {
    tets.iter()
        .map(|&tet| {
            let (volume, _) = reference(nodes, tet);
            RHO * volume / 40.0
                * (0..3)
                    .map(|i| {
                        tet.iter().map(|&n| velocity[n][i].powi(2)).sum::<f64>()
                            + tet.iter().map(|&n| velocity[n][i]).sum::<f64>().powi(2)
                    })
                    .sum::<f64>()
        })
        .sum()
}

fn material_force(nodes: &[[f64; 3]], tets: &[[usize; 4]], nominal: f64) -> Vec<[f64; 3]> {
    let mut force = vec![[0.0; 3]; nodes.len()];
    for &tet in tets {
        let (volume, gradient) = reference(nodes, tet);
        for (local, node) in tet.into_iter().enumerate() {
            for i in 0..3 {
                force[node][i] += volume * nominal * gradient[local][i];
            }
        }
    }
    force
}

fn complete_balance(result: &FinitePorousCavitiesMotionStep<'_>, tolerance: f64) {
    let step = &result.step;
    let energy = &step.material_energy;
    let held = step.pressure.as_ref().map_or(0.0, |work| work.work_j);
    let held_correction = step
        .pressure
        .as_ref()
        .map_or(0.0, |work| work.endpoint_work_correction_j);
    let exterior = result
        .cavity_energy
        .iter()
        .map(|work| work.exterior_work_j)
        .sum::<f64>();
    let gas_change = result
        .cavity_energy
        .iter()
        .map(|work| work.gas_stored_energy_change_j)
        .sum::<f64>();
    near(
        energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j + held + exterior,
        energy.stored_after_j - energy.stored_before_j
            + gas_change
            + energy.plastic_dissipation_j
            + step.inertia.kinetic_change_j
            + step.inertia.numerical_loss_j
            + result.external_work_defect_j
            + step.balance_error_j,
        tolerance,
    );
    near(
        result.external_work_defect_j,
        energy.endpoint_work_defect_j
            + held_correction
            + result
                .cavity_energy
                .iter()
                .map(|work| work.endpoint_work_correction_j)
                .sum::<f64>(),
        tolerance,
    );
    near(step.balance_error_j, 0.0, tolerance);
}

#[test]
fn free_cavity_compaction_solves_gas_stiffness_and_consistent_momentum() {
    with_cx(|cx, _| {
        let (mut nodes, tets) = shell();
        // Both outer and inner tetrahedra have centroid (3/4,3/4,3/4).
        // Centering eliminates rigid acceleration from the manufactured crush.
        for x in &mut nodes {
            for coordinate in x {
                *coordinate -= 0.75;
            }
        }
        let surface = surface(&nodes, &tets, cx);
        let cavity = cavity(&surface, 4, cx);
        let material = material();
        let initial = motion(&nodes, &tets, &material, &vec![[0.0; 3]; nodes.len()], cx);
        let stretch = 0.94_f64;
        let dt = 0.25;
        let (alpha, plastic_j, tau) = hydrostatic(stretch);
        let pressure = -tau / stretch.powi(3);
        let exterior = 0.5;
        // At the target inner gas and outer pressure balance the analytic
        // material stress. Independently integrated M*a supplies the nonzero
        // acceleration; omitting inertia does not solve this manufactured load.
        let gas = gas_fixture::gas(
            cavity.reference_volume_m3(),
            (pressure + exterior) * stretch.powf(4.2),
            0.0,
            1e24,
        );
        let original = gas.clone();
        let source = [binding(&cavity, &gas, exterior)];
        let face_pressure: Vec<_> = (0..surface.faces().len())
            .map(|face| {
                if cavity.solid_face_indices().contains(&face) {
                    0.0
                } else {
                    pressure
                }
            })
            .collect();
        let acceleration: Vec<_> = nodes
            .iter()
            .map(|x| x.map(|v| (stretch - 1.0) * v / (dt * dt)))
            .collect();
        let force = mass_action(&nodes, &tets, &acceleration);
        let loads: Vec<_> = force.iter().flatten().copied().enumerate().collect();
        let result = initial
            .advance_with_cavities(
                &[],
                &loads,
                [0.0; 3],
                dt,
                &source,
                Some((surface.as_ref(), TetPressure::PerFace(&face_pressure))),
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert!(result.step.equilibrium.reactions_n.is_empty());
        assert!(result.step.equilibrium.steps[0].linear_iterations > 0);
        for (node, x) in nodes.iter().enumerate() {
            for i in 0..3 {
                near(
                    result.step.motion.solid().displacement_m()[node][i],
                    (stretch - 1.0) * x[i],
                    8e-7,
                );
                near(
                    result.step.motion.velocity_mps()[node][i],
                    (stretch - 1.0) * x[i] / dt,
                    3e-6,
                );
            }
        }
        for (state, &j) in result
            .step
            .motion
            .solid()
            .element_states()
            .iter()
            .zip(result.step.motion.solid().element_plastic_jacobians())
        {
            near(state.alpha, alpha, 8e-7);
            near(j, plastic_j, 8e-7);
        }
        for &density in &result.step.motion.current_density_kg_m3(cx).unwrap() {
            near(density, RHO / stretch.powi(3), 2e-5);
        }
        let material_volume = 13.0 / 3.0;
        near(
            result.step.motion.reference_mass_kg(),
            RHO * material_volume,
            1e-12,
        );
        near(
            result.step.material_energy.stored_after_j,
            material_volume * (0.5 * tau * tau / K + 0.5 * H * alpha * alpha),
            3e-6,
        );
        near(
            result.step.material_energy.plastic_dissipation_j,
            material_volume * Y * alpha,
            3e-6,
        );
        near(
            result.step.material_energy.internal_endpoint_work_j,
            3.0 * material_volume * tau / stretch * (stretch - 1.0),
            3e-6,
        );
        // The independent shell polar second moment about its centroid is
        // integral |X|^2 dV = (9/2)*(81/80) - (1/6)*(9/80) = 363/80.
        let expected_kinetic = 0.5 * RHO * (363.0 / 80.0) * (stretch - 1.0).powi(2) / (dt * dt);
        near(result.step.inertia.kinetic_change_j, expected_kinetic, 4e-6);
        near(result.step.inertia.numerical_loss_j, expected_kinetic, 4e-6);
        near(
            result.step.inertia.inertial_work_j,
            2.0 * expected_kinetic,
            6e-6,
        );
        near(
            result.step.material_energy.dead_load_endpoint_work_j,
            2.0 * expected_kinetic,
            6e-6,
        );
        assert!(result.step.inertia.kinetic_after_j > 0.7);
        let gas_temperature = 300.0 * stretch.powf(-1.2);
        near(
            result.cavities[0].gas.state().state().temperature,
            gas_temperature,
            8e-5,
        );
        near(result.cavities[0].net_pressure_pa, pressure, 8e-6);
        near(
            result.cavities[0].volume_m3,
            cavity.reference_volume_m3() * stretch.powi(3),
            5e-7,
        );
        let work = result.cavity_energy[0];
        // The 1e24 J/kg datum makes absolute stored-energy subtraction unusable.
        near(
            work.gas_stored_energy_change_j,
            gas_energy(&gas, 0.0, 300.0, gas_temperature),
            2e-6,
        );
        near(work.gas_work_j, -work.gas_stored_energy_change_j, 0.0);
        near(
            work.exterior_work_j,
            exterior * cavity.reference_volume_m3() * (1.0 - stretch.powi(3)),
            2e-7,
        );
        near(
            work.endpoint_work_j,
            3.0 * pressure * cavity.reference_volume_m3() * stretch.powi(2) * (stretch - 1.0),
            2e-6,
        );
        near(
            result.step.pressure.unwrap().work_j,
            pressure * 4.5 * (1.0 - stretch.powi(3)),
            6e-6,
        );
        complete_balance(&result, 6e-6);
        assert_eq!(
            result.cavities[0].gas.mass_kg().to_bits(),
            gas.mass_kg().to_bits()
        );
        assert_eq!(result.cavities[0].gas.model(), gas.model());
        assert_eq!(gas, original);
        assert!(
            initial
                .solid()
                .element_states()
                .iter()
                .all(|state| state.alpha == 0.0)
        );
        assert_eq!(initial.velocity_mps(), vec![[0.0; 3]; nodes.len()]);
    });
}

#[test]
fn retained_cavities_release_with_independent_support_work_and_atomic_kinetic_gates() {
    with_cx(|cx, gate| {
        let (nodes, tets) = connected_shells();
        let surface = surface(&nodes, &tets, cx);
        let cavities = [cavity(&surface, 4, cx), cavity(&surface, 9, cx)];
        let slopes = [0.0, 0.002];
        let pressures = [5.0, 7.0];
        let exteriors = [0.7, 1.1];
        let gases = [
            gas_fixture::gas(
                cavities[0].reference_volume_m3(),
                pressures[0],
                slopes[0],
                0.0,
            ),
            gas_fixture::gas(
                cavities[1].reference_volume_m3(),
                pressures[1],
                slopes[1],
                1e24,
            ),
        ];
        let original_gases = gases.clone();
        // Deliberately reverse source order; connected solid nodes do not mix gases.
        let sources = [
            binding(&cavities[1], &gases[1], exteriors[1]),
            binding(&cavities[0], &gases[0], exteriors[0]),
        ];
        let material = material();
        let initial_velocity = [0.2, -0.1, 0.05];
        let initial = motion(
            &nodes,
            &tets,
            &material,
            &vec![initial_velocity; nodes.len()],
            cx,
        );
        let peak = 0.92_f64;
        let first_dt = 0.3;
        let (alpha, plastic_j, _) = hydrostatic(peak);
        let first = initial
            .advance_with_cavities(
                &affine(&nodes, peak),
                &[],
                [0.0; 3],
                first_dt,
                &sources,
                None,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let first_velocity: Vec<_> = nodes
            .iter()
            .map(|x| x.map(|v| (peak - 1.0) * v / first_dt))
            .collect();
        let first_difference: Vec<[f64; 3]> = first_velocity
            .iter()
            .map(|v| std::array::from_fn(|i| v[i] - initial_velocity[i]))
            .collect();
        near(
            first.step.inertia.kinetic_before_j,
            kinetic(&nodes, &tets, &vec![initial_velocity; nodes.len()]),
            2e-8,
        );
        near(
            first.step.inertia.kinetic_after_j,
            kinetic(&nodes, &tets, &first_velocity),
            2e-8,
        );
        near(
            first.step.inertia.numerical_loss_j,
            kinetic(&nodes, &tets, &first_difference),
            2e-8,
        );
        near(
            first.step.material_energy.plastic_dissipation_j,
            (26.0 / 3.0) * Y * alpha,
            2e-6,
        );
        complete_balance(&first, 3e-6);

        let accepted = &first.step.motion;
        let original_history = accepted.solid().element_states().to_vec();
        let original_u = accepted.solid().displacement_m().to_vec();
        let original_velocity = accepted.velocity_mps().to_vec();
        let original_storage = accepted.solid().stored_energy_j();
        // Stop halfway to the stress-free plastic stretch: elastic release
        // retains compaction even while support motion reverses the velocity.
        let stretch = (peak * plastic_j.cbrt()).sqrt();
        let dt = 0.2;
        let supports = affine(&nodes, stretch);
        let velocity: Vec<_> = nodes
            .iter()
            .map(|x| x.map(|v| (stretch - peak) * v / dt))
            .collect();
        let difference: Vec<[f64; 3]> = velocity
            .iter()
            .zip(&first_velocity)
            .map(|(after, before)| sub(*after, *before))
            .collect();
        let acceleration: Vec<_> = difference.iter().map(|v| v.map(|x| x / dt)).collect();
        let second = accepted
            .advance_with_cavities(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &sources,
                None,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        for (before, after) in original_history
            .iter()
            .zip(second.step.motion.solid().element_states())
        {
            near(after.alpha, before.alpha, 2e-8);
        }
        for &j in second.step.motion.solid().element_plastic_jacobians() {
            near(j, plastic_j, 2e-8);
        }
        near(second.step.material_energy.plastic_dissipation_j, 0.0, 2e-8);
        near(
            second.step.inertia.kinetic_before_j,
            kinetic(&nodes, &tets, &first_velocity),
            2e-8,
        );
        near(
            second.step.inertia.kinetic_after_j,
            kinetic(&nodes, &tets, &velocity),
            2e-8,
        );
        near(
            second.step.inertia.kinetic_change_j,
            kinetic(&nodes, &tets, &velocity) - kinetic(&nodes, &tets, &first_velocity),
            3e-8,
        );
        near(
            second.step.inertia.numerical_loss_j,
            kinetic(&nodes, &tets, &difference),
            3e-8,
        );
        assert!(second.step.inertia.numerical_loss_j > 1.0);

        let tau = K * (stretch.powi(3) / plastic_j).ln();
        let mut reactions = material_force(&nodes, &tets, tau / stretch);
        let inertia = mass_action(&nodes, &tets, &acceleration);
        for (force, added) in reactions.iter_mut().zip(inertia) {
            for i in 0..3 {
                force[i] += added[i];
            }
        }
        for (output, source) in [1_usize, 0].into_iter().enumerate() {
            let t1 = temperature(slopes[source], peak.powi(3));
            let t2 = temperature(slopes[source], stretch.powi(3));
            let net_pressure =
                pressures[source] * (t2 / 300.0) / stretch.powi(3) - exteriors[source];
            let volume = cavities[source].reference_volume_m3();
            let work = second.cavity_energy[output];
            near(
                second.cavities[output].gas.state().state().temperature,
                t2,
                3e-7,
            );
            near(second.cavities[output].net_pressure_pa, net_pressure, 3e-8);
            near(
                second.cavities[output].volume_m3,
                volume * stretch.powi(3),
                2e-14,
            );
            near(
                work.gas_stored_energy_change_j,
                gas_energy(&gases[source], slopes[source], t1, t2),
                3e-9,
            );
            near(
                first.cavity_energy[output].gas_stored_energy_change_j
                    + work.gas_stored_energy_change_j,
                gas_energy(&gases[source], slopes[source], 300.0, t2),
                4e-9,
            );
            near(
                work.exterior_work_j,
                -exteriors[source] * volume * (stretch.powi(3) - peak.powi(3)),
                2e-12,
            );
            near(
                work.endpoint_work_j,
                3.0 * net_pressure * volume * stretch.powi(2) * (stretch - peak),
                2e-9,
            );
            near(work.work_j, work.gas_work_j + work.exterior_work_j, 1e-13);
            assert_eq!(
                second.cavities[output].gas.mass_kg().to_bits(),
                gases[source].mass_kg().to_bits()
            );
            assert_eq!(second.cavities[output].gas.model(), gases[source].model());
            // Independent gas-volume derivative on the four original inner nodes.
            let first_node = if source == 0 { 4 } else { 9 };
            let tet = [first_node, first_node + 1, first_node + 2, first_node + 3];
            let (reference_volume, gradient) = reference(&nodes, tet);
            for (local, node) in tet.into_iter().enumerate() {
                for i in 0..3 {
                    reactions[node][i] -=
                        net_pressure * reference_volume * stretch.powi(2) * gradient[local][i];
                }
            }
        }
        assert_eq!(second.step.equilibrium.reactions_n.len(), 3 * nodes.len());
        for &(dof, force) in &second.step.equilibrium.reactions_n {
            near(force, reactions[dof / 3][dof % 3], 4e-6);
        }
        let expected_support_work: f64 = reactions
            .iter()
            .zip(&nodes)
            .map(|(force, x)| dot(*force, x.map(|v| (stretch - peak) * v)))
            .sum();
        near(
            second.step.material_energy.support_endpoint_work_j,
            expected_support_work,
            3e-6,
        );
        complete_balance(&second, 3e-6);
        near(
            second.step.motion.reference_mass_kg(),
            RHO * (26.0 / 3.0),
            1e-12,
        );

        let mut strict = motion_controls();
        strict.maximum_numerical_kinetic_loss_j = 0.0;
        let refused = accepted
            .advance_with_cavities(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &sources,
                None,
                settings(),
                strict,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePorousTetError::Solid(FinitePlasticTetError::EnergyGate {
                quantity: "porous backward-Euler numerical kinetic loss",
                ..
            })
        ));
        let retry = accepted
            .advance_with_cavities(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &sources,
                None,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.step.equilibrium, second.step.equilibrium);
        assert_eq!(retry.step.material_energy, second.step.material_energy);
        assert_eq!(retry.step.inertia, second.step.inertia);
        assert_eq!(retry.cavity_energy, second.cavity_energy);
        near(
            retry.external_work_defect_j,
            second.external_work_defect_j,
            0.0,
        );

        // The empty-cavity overload preserves the original pressure dynamics.
        let zero = accepted
            .advance_with_pressure(
                &supports,
                &[],
                [0.0; 3],
                dt,
                surface.as_ref(),
                0.0,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        let empty = accepted
            .advance_with_cavities::<ClosedGasParcel>(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &[],
                Some((surface.as_ref(), TetPressure::Uniform(0.0))),
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap();
        assert!(empty.cavities.is_empty() && empty.cavity_energy.is_empty());
        assert_eq!(empty.step.equilibrium, zero.equilibrium);
        assert_eq!(empty.step.material_energy, zero.material_energy);
        assert_eq!(empty.step.inertia, zero.inertia);
        assert_eq!(empty.step.pressure, zero.pressure);
        near(
            empty.external_work_defect_j,
            zero.pressure.unwrap().external_work_defect_j,
            0.0,
        );

        // A deformed gas endpoint cannot become the undeformed-volume reference.
        assert!(
            accepted
                .advance_with_cavities(
                    &supports,
                    &[],
                    [0.0; 3],
                    dt,
                    &[binding(&cavities[0], &first.cavities[1].gas, exteriors[0])],
                    None,
                    settings(),
                    motion_controls(),
                    cx,
                )
                .is_err()
        );
        gate.request();
        let cancelled = accepted
            .advance_with_cavities(
                &supports,
                &[],
                [0.0; 3],
                dt,
                &sources,
                None,
                settings(),
                motion_controls(),
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            cancelled,
            FinitePorousTetError::Solid(FinitePlasticTetError::Native(HyperTetError::Cancelled))
        ));
        assert_eq!(accepted.solid().element_states(), original_history);
        assert_eq!(accepted.solid().displacement_m(), original_u);
        assert_eq!(accepted.velocity_mps(), original_velocity);
        near(accepted.solid().stored_energy_j(), original_storage, 0.0);
        near(accepted.reference_mass_kg(), RHO * (26.0 / 3.0), 1e-12);
        assert_eq!(gases, original_gases);
    });
}
