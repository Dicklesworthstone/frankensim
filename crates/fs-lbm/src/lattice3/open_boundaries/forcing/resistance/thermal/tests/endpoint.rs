//! Actual evolving-law velocity in thermal source kicks and joint rollback.

use super::*;

fn source_epochs<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let cells = 24;
        let mut law = DarcyForchheimer3::from(LinearResistance3 {
            principal: [0.24; 3],
            skeleton_velocity: [0.01, 0.005, -0.002],
            body_force: [0.0003, -0.0002, 0.0001],
            ..LinearResistance3::default()
        });
        let mut owner = ResistanceOpenGrid3::try_new_preserving_velocity(
            native::<L, Q>(cx, [4, 3, 2], collision, BoundarySpec3::periodic(), |_| {
                [0.04, -0.02, 0.03]
            }),
            cx,
            vec![law; cells],
            [0.0; 3],
        )
        .unwrap();
        let mut heat = owner
            .prepare_thermal_field(cx, &vec![0.7; cells], 0.83, periodic_heat(), limits())
            .unwrap();
        let tensor = DiffusivityTensor::from_lattice([
            [0.08, 0.015, 0.0],
            [0.015, 0.1, -0.01],
            [0.0, -0.01, 0.09],
        ])
        .unwrap();
        let increment = 0.004;
        for tick in 1..=4 {
            let u0 = owner.flow().macroscopic(0).unwrap().1;
            let force0 = owner.flow().forces()[0];
            let next = DarcyForchheimer3::from(LinearResistance3 {
                principal: [0.11 + 0.02 * tick as f64; 3],
                skeleton_velocity: [-0.003 * tick as f64, 0.004, 0.002 * tick as f64],
                body_force: [-0.0002, 0.0001 * tick as f64, 0.0004],
                ..LinearResistance3::default()
            });
            // Uniform physical trapezoidal momentum, independent of native
            // pull/reconstruction and its nonlinear constitutive solver.
            let expected: [f64; 3] = core::array::from_fn(|a| {
                (RHO * u0[a]
                    + 0.5 * force0[a]
                    + 0.5 * next.linear.body_force[a]
                    + 0.5 * next.linear.principal[a] * next.linear.skeleton_velocity[a])
                    / (RHO + 0.5 * next.linear.principal[a])
            });
            let r = if tick % 2 == 0 {
                owner
                    .advance_serial_with_thermal_tensor_diffusion(
                        &mut heat,
                        cx,
                        None,
                        vec![next; cells],
                        1.0,
                        &vec![increment; cells],
                        DiffusivityField::Uniform(&tensor),
                    )
                    .unwrap()
            } else {
                owner
                    .advance_serial_with_thermal(
                        &mut heat,
                        cx,
                        None,
                        vec![next; cells],
                        1.0,
                        &vec![increment; cells],
                        None,
                    )
                    .unwrap()
            };
            for cell in 0..cells {
                let (rho, u) = owner.flow().macroscopic(cell).unwrap();
                near(rho, RHO, 3e-13);
                for a in 0..3 {
                    near(u[a], expected[a], 3e-13);
                    near(
                        owner.flow().forces()[cell][a],
                        next.linear.body_force[a]
                            + next.linear.principal[a]
                                * (next.linear.skeleton_velocity[a] - expected[a]),
                        3e-13,
                    );
                    if tick == 1 {
                        let (plus, minus) = crate::thermal::advection::D3Q7.direction_pair(a);
                        let g = heat.populations()[cell];
                        near(
                            g[plus] - g[minus],
                            (0.7 + 0.5 * increment) * u0[a] + 0.5 * increment * expected[a],
                            4e-15,
                        );
                    }
                }
                near(
                    heat.temperature(cell).unwrap(),
                    0.7 + tick as f64 * increment,
                    3e-13,
                );
            }
            let initial_power: f64 = (0..3)
                .map(|a| {
                    law.linear.principal[a] * (u0[a] - law.linear.skeleton_velocity[a]).powi(2)
                })
                .sum();
            let final_power: f64 = (0..3)
                .map(|a| {
                    next.linear.principal[a]
                        * (expected[a] - next.linear.skeleton_velocity[a]).powi(2)
                })
                .sum();
            near(
                r.momentum.resistance.endpoint_relative_dissipation,
                0.5 * cells as f64 * (initial_power + final_power),
                2e-13,
            );
            // Explicit supplied heat alone enters thermal storage; temporal
            // work and drag remain independently reported mechanical quantities.
            near(r.heat.source, cells as f64 * increment, 2e-15);
            near(r.heat.stored_change, cells as f64 * increment, 3e-13);
            near(r.heat.defect, 0.0, 3e-13);
            assert_eq!(owner.laws(), vec![next; cells]);
            assert_eq!(owner.last_endpoint_exchange(), Some(r.momentum));
            assert_eq!(owner.last_exchange(), None);
            assert_eq!(owner.flow().completed_steps(), heat.completed_steps());
            law = next;
        }
    });
}

#[test]
fn evolving_porous_velocity_drives_both_thermal_source_epochs_and_keeps_work_separate() {
    source_epochs::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    source_epochs::<D3q27, 27>(bgk27());
}

fn changed(laws: &[DarcyForchheimer3]) -> Vec<DarcyForchheimer3> {
    laws.iter()
        .map(|law| {
            let mut next = *law;
            next.linear.principal = next.linear.principal.map(|v| 0.7 * v);
            next.linear.body_force = [-0.00003, 0.00001, 0.00002];
            next.linear.skeleton_velocity = [0.007, -0.002, 0.003];
            next.quadratic *= 0.6;
            next
        })
        .collect()
}

#[test]
fn late_thermal_work_or_cancel_refusal_keeps_both_fields_and_laws_retryable() {
    with_cx(|cx| {
        let source = vec![0.001; 48];
        let (mut expected, mut expected_heat) = channel(cx);
        let laws = changed(expected.laws());
        let mut polls = 0;
        let expected_receipt = expected
            .advance_with_thermal_poll(
                &mut expected_heat,
                None,
                laws.clone(),
                1.0,
                &source,
                None,
                None,
                &mut || {
                    polls += 1;
                    Ok(())
                },
            )
            .unwrap();
        assert!(
            expected_receipt
                .momentum
                .resistance
                .absolute_relative_quadrature_remainder
                > 0.0
        );
        for failure in 0..3 {
            let (mut owner, mut heat) = channel(cx);
            let before = accepted(&owner, &heat);
            let pending_before = owner.flow.pending_forces.clone();
            let excessive_heat: Vec<_> = (0..48)
                .map(|cell| 1.8 * (1.0 - heat.temperature(cell).unwrap()))
                .collect();
            let mut visits = 0;
            let result = owner.advance_with_thermal_poll(
                &mut heat,
                None,
                laws.clone(),
                if failure == 1 { 0.0 } else { 1.0 },
                if failure == 0 {
                    &excessive_heat
                } else {
                    &source
                },
                None,
                None,
                &mut || {
                    visits += 1;
                    if failure == 2 && visits == polls {
                        Err(GridError3::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            );
            match failure {
                0 => assert!(matches!(
                    result,
                    Err(OpenThermalError3::Thermal(
                        ThermalOpenError3::Temperature { .. }
                    ))
                )),
                1 => assert!(matches!(
                    result,
                    Err(OpenThermalError3::Flow(GridError3::Input(_)))
                )),
                _ => assert!(matches!(
                    result,
                    Err(OpenThermalError3::Flow(GridError3::Cancelled))
                )),
            }
            assert_ne!(
                owner.flow.pending_forces, pending_before,
                "momentum was staged before refusal"
            );
            assert_eq!(accepted(&owner, &heat), before);
            let retry = owner
                .advance_serial_with_thermal(&mut heat, cx, None, laws.clone(), 1.0, &source, None)
                .unwrap();
            assert_eq!(retry, expected_receipt);
            assert_eq!(accepted(&owner, &heat), accepted(&expected, &expected_heat));
            let mut storage = retry.heat.stored_change;
            let mut supplied = retry.heat.source;
            for reservoir in retry.heat.reservoirs {
                supplied += reservoir.heat_in;
            }
            for wall in retry.heat.walls.into_iter().flatten() {
                let solid = wall.solid.unwrap();
                storage += solid.stored_change;
                supplied += solid.source;
            }
            near(storage, supplied, 2e-12);
        }
    });
}

#[test]
fn unchanged_law_thermal_transaction_has_original_population_bits() {
    with_cx(|cx| {
        let (mut held, mut held_heat) = channel(cx);
        let (mut evolving, mut evolving_heat) = channel(cx);
        let laws = held.laws().to_vec();
        let source = vec![0.001; 48];
        let a = held
            .step_serial_with_thermal(&mut held_heat, cx, None, &source, None)
            .unwrap();
        let b = evolving
            .advance_serial_with_thermal(&mut evolving_heat, cx, None, laws, 1.0, &source, None)
            .unwrap();
        assert_eq!(held.flow().populations(), evolving.flow().populations());
        assert_eq!(held.flow().forces(), evolving.flow().forces());
        assert_eq!(held.drag, evolving.drag);
        assert_eq!(held_heat.populations(), evolving_heat.populations());
        assert_eq!(a.heat, b.heat);
        assert_eq!(a.momentum.open, b.momentum.open);
        assert_eq!(a.momentum.walls, b.momentum.resistance.walls);
        assert_eq!(
            a.momentum.fluid_drag_work,
            b.momentum.resistance.fluid_drag_work.impulse_work
        );
        for face in 0..6 {
            assert_eq!(held_heat.wall_states(face), evolving_heat.wall_states(face));
        }
    });
}
