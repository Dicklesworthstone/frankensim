use super::*;

fn advance_background<L: EquilibriumLattice3<Q>, const Q: usize>(
    coupled: &mut BoussinesqOpenGrid3<L, Q>,
    background_at_end: &mut Vec<[f64; 3]>,
    heat: &[f64],
    diffusion: Option<DiffusivityField<'_, 3>>,
) -> BoussinesqOpenStep3 {
    coupled
        .advance_with_diffusion_poll(
            Some(background_at_end),
            None,
            heat,
            None,
            diffusion,
            [0.0; 3],
            &mut || Ok(()),
        )
        .unwrap()
}

fn check_affine_mixed_convection<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
) {
    let mut coupled = periodic::<L, Q>(collision, 0.5);
    let cells = coupled.flow().populations().len();
    let source = 0.002;
    let ramp = [3e-6, -2e-6, 1e-6];
    let mut endpoint = vec![[0.0; 3]; cells];
    let thermal_slope = GRAVITY.map(|g| -REFERENCE_RHO * BETA * g);
    for step in 1..=20 {
        let n = step as f64;
        let previous_background = coupled.background_forces().to_vec();
        let authored = core::array::from_fn(|axis| BACKGROUND[axis] + n * ramp[axis]);
        endpoint.fill(authored);
        let report = advance_background(&mut coupled, &mut endpoint, &vec![source; cells], None);
        assert_eq!(endpoint, previous_background);
        assert!(coupled.background_forces().iter().all(|b| *b == authored));
        for cell in 0..cells {
            let (rho, velocity) = coupled.flow().macroscopic(cell).unwrap();
            near(rho, RHO, 2e-12);
            near(
                coupled.thermal().temperature(cell).unwrap(),
                0.5 + n * source,
                2e-13,
            );
            for axis in 0..3 {
                // Exact integral for an affine mechanical force plus buoyancy
                // from uniform heating, using physical momentum rho*u.
                let expected = INITIAL_VELOCITY[axis]
                    + (n * BACKGROUND[axis]
                        + 0.5 * n * n * (ramp[axis] + thermal_slope[axis] * source))
                        / RHO;
                near(velocity[axis], expected, 2e-12);
                near(
                    coupled.flow().forces()[cell][axis],
                    authored[axis] + thermal_slope[axis] * n * source,
                    2e-15,
                );
            }
        }
        for axis in 0..3 {
            let mechanical = cells as f64 * (BACKGROUND[axis] + (n - 0.5) * ramp[axis]);
            let buoyant = cells as f64 * thermal_slope[axis] * source * (n - 0.5);
            near(report.background_impulse[axis], mechanical, 2e-14);
            near(report.buoyancy_impulse[axis], buoyant, 2e-13);
            near(report.body_impulse[axis], mechanical + buoyant, 2e-13);
        }
        near(report.heat.combined_defect, 0.0, 2e-12);
        assert!(report.maximum_force_defect <= 1e-12);
    }
}

#[test]
fn affine_background_and_heating_integrate_their_distinct_impulses_on_both_lattices() {
    check_affine_mixed_convection::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    check_affine_mixed_convection::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

#[test]
fn identical_endpoint_background_preserves_held_scalar_and_tensor_steps_bitwise() {
    let mut explicit = channel(BETA);
    let mut held = channel(BETA);
    let coefficient =
        DiffusivityTensor::from_lattice([[0.09, 0.035, 0.0], [0.035, 0.08, 0.0], [0.0, 0.0, 0.06]])
            .unwrap();
    let field = vec![coefficient; 72];
    let heat: Vec<_> = (0..72).map(|cell| 0.0002 * (1 + cell % 5) as f64).collect();
    let mut endpoint = vec![BACKGROUND; 72];
    for step in 0..6 {
        let expected = if step % 2 == 0 {
            advance(&mut held, &heat)
        } else {
            tensor_advance(&mut held, &heat, DiffusivityField::PerCell(&field))
        };
        let actual = advance_background(
            &mut explicit,
            &mut endpoint,
            &heat,
            if step % 2 == 0 {
                None
            } else {
                Some(DiffusivityField::PerCell(&field))
            },
        );
        assert_eq!(actual, expected);
        assert_eq!(accepted(&explicit), accepted(&held));
        assert_eq!(endpoint, vec![BACKGROUND; 72]);
    }
}

fn spatial_background(step: f64) -> Vec<[f64; 3]> {
    (0..72)
        .map(|cell| {
            [
                BACKGROUND[0] + step * 1e-6 * (1 + cell % 6) as f64,
                BACKGROUND[1] - step * 0.5e-6 * (1 + (cell / 6) % 4) as f64,
                BACKGROUND[2] + step * 0.25e-6 * (1 + cell / 24) as f64,
            ]
        })
        .collect()
}

#[test]
fn changing_spatial_background_retains_tensor_buoyancy_and_conducting_wall_energy() {
    let mut coupled = channel(BETA);
    let coefficient =
        DiffusivityTensor::from_lattice([[0.09, 0.035, 0.0], [0.035, 0.08, 0.0], [0.0, 0.0, 0.06]])
            .unwrap();
    let heat: Vec<_> = (0..72).map(|cell| 0.0002 * (1 + cell % 5) as f64).collect();
    let mut endpoint = vec![BACKGROUND; 72];
    for step in 1..=8 {
        let energy_before = total_energy(&coupled);
        let opening = coupled.background_forces().to_vec();
        let opening_temperature: Vec<_> = (0..72)
            .map(|cell| coupled.thermal().temperature(cell).unwrap())
            .collect();
        let authored = spatial_background(step as f64);
        endpoint.copy_from_slice(&authored);
        let report = advance_background(
            &mut coupled,
            &mut endpoint,
            &heat,
            Some(DiffusivityField::Uniform(&coefficient)),
        );
        assert_eq!(endpoint, opening);
        assert_eq!(coupled.background_forces(), authored);
        let reservoir_heat: f64 = report.heat.reservoirs.iter().map(|r| r.heat_in).sum();
        near(
            total_energy(&coupled) - energy_before,
            heat.iter().sum::<f64>() + reservoir_heat,
            2e-11,
        );
        near(report.heat.combined_defect, 0.0, 2e-12);
        for axis in 0..3 {
            let expected_background: f64 = (0..72)
                .map(|cell| 0.5 * opening[cell][axis] + 0.5 * authored[cell][axis])
                .sum();
            let expected_buoyancy: f64 = (0..72)
                .map(|cell| {
                    -REFERENCE_RHO
                        * BETA
                        * (0.5
                            * (opening_temperature[cell]
                                + coupled.thermal().temperature(cell).unwrap())
                            - 0.5)
                        * GRAVITY[axis]
                })
                .sum();
            near(report.background_impulse[axis], expected_background, 2e-14);
            near(report.buoyancy_impulse[axis], expected_buoyancy, 2e-13);
            near(
                report.body_impulse[axis],
                expected_background + expected_buoyancy,
                2e-13,
            );
            for cell in 0..72 {
                let temperature = coupled.thermal().temperature(cell).unwrap();
                let expected = authored[cell][axis]
                    - REFERENCE_RHO * BETA * (temperature - 0.5) * GRAVITY[axis];
                near(coupled.flow().forces()[cell][axis], expected, 2e-15);
            }
        }
    }
}

fn force_bits(values: &[[f64; 3]]) -> Vec<[u64; 3]> {
    values.iter().map(|v| v.map(f64::to_bits)).collect()
}

#[test]
fn varying_background_refusal_and_final_cancellation_leave_input_and_joint_epoch_intact() {
    let mut actual = channel(BETA);
    let mut reference = channel(BETA);
    let coefficient =
        DiffusivityTensor::from_lattice([[0.09, 0.035, 0.0], [0.035, 0.08, 0.0], [0.0, 0.0, 0.06]])
            .unwrap();
    let field = vec![coefficient; 72];
    let heat = vec![0.0004; 72];
    advance(&mut actual, &heat);
    advance(&mut reference, &heat);
    let before = accepted(&actual);

    let mut malformed = spatial_background(1.0);
    malformed.pop();
    let short_before = force_bits(&malformed);
    let refused = actual.advance_with_diffusion_poll(
        Some(&mut malformed),
        None,
        &heat,
        None,
        None,
        [0.0; 3],
        &mut || Ok(()),
    );
    assert!(matches!(refused, Err(BoussinesqOpenError3::Invalid(_))));
    assert_eq!(force_bits(&malformed), short_before);
    assert_eq!(accepted(&actual), before);

    malformed = spatial_background(1.0);
    malformed[35][1] = f64::NAN;
    let nan_before = force_bits(&malformed);
    let refused = actual.advance_with_diffusion_poll(
        Some(&mut malformed),
        None,
        &heat,
        None,
        None,
        [0.0; 3],
        &mut || Ok(()),
    );
    assert!(matches!(refused, Err(BoussinesqOpenError3::Invalid(_))));
    assert_eq!(force_bits(&malformed), nan_before);
    assert_eq!(accepted(&actual), before);

    let mut endpoint = spatial_background(1.0);
    let endpoint_before = endpoint.clone();
    let refused = actual.advance_with_diffusion_poll(
        Some(&mut endpoint),
        None,
        &vec![2.2; 72],
        None,
        None,
        [0.0; 3],
        &mut || Ok(()),
    );
    assert!(
        matches!(
            refused,
            Err(BoussinesqOpenError3::Law {
                source: BoussinesqLawError3::DensityDeparture { .. },
                ..
            }) | Err(BoussinesqOpenError3::BoundaryLaw {
                source: BoussinesqLawError3::DensityDeparture { .. },
                ..
            })
        ),
        "{refused:?}"
    );
    assert_eq!(endpoint, endpoint_before);
    assert_eq!(accepted(&actual), before);

    let mut reference_endpoint = endpoint.clone();
    let mut polls = 0;
    let expected = reference
        .advance_with_diffusion_poll(
            Some(&mut reference_endpoint),
            None,
            &heat,
            None,
            Some(DiffusivityField::PerCell(&field)),
            [0.0; 3],
            &mut || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    let mut calls = 0;
    let refused = actual.advance_with_diffusion_poll(
        Some(&mut endpoint),
        None,
        &heat,
        None,
        Some(DiffusivityField::PerCell(&field)),
        [0.0; 3],
        &mut || {
            calls += 1;
            if calls == polls {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(BoussinesqOpenError3::Flow(GridError3::Cancelled))
    ));
    assert_eq!(calls, polls);
    assert_eq!(endpoint, endpoint_before);
    assert_eq!(accepted(&actual), before);
    let retried = advance_background(
        &mut actual,
        &mut endpoint,
        &heat,
        Some(DiffusivityField::PerCell(&field)),
    );
    assert_eq!(retried, expected);
    assert_eq!(endpoint, reference_endpoint);
    assert_eq!(accepted(&actual), accepted(&reference));
}
