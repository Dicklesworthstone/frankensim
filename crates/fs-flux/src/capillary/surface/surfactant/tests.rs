use super::*;

fn laws() -> [SurfactantLaw; 2] {
    [
        SurfactantLaw::IdealDilute {
            clean_tension_n_m: 1.4,
            thermal_energy_j_mol: 0.7,
            reference_concentration_mol_m2: 0.3,
        },
        SurfactantLaw::Langmuir {
            clean_tension_n_m: 1.4,
            thermal_energy_j_mol: 0.7,
            saturation_concentration_mol_m2: 1.0,
        },
    ]
}

fn triangle() -> [Vector; 3] {
    [[0.0; 3], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}

// An independent endpoint expression using the platform log as an oracle.
fn energy(law: SurfactantLaw, area: f64, n: f64) -> f64 {
    let (clean, rt, reference) = law.parameters();
    if n == 0.0 {
        return clean * area;
    }
    let gamma = n / area;
    match law {
        SurfactantLaw::IdealDilute { .. } => {
            clean * area + rt * n * ((gamma / reference).ln() - 1.0)
        }
        SurfactantLaw::Langmuir { .. } => {
            let theta = gamma / reference;
            area * (clean
                + rt * reference * (theta * theta.ln() + (1.0 - theta) * (-theta).ln_1p()))
        }
    }
}

#[test]
fn g1_equations_of_state_follow_the_free_energy_area_and_amount_derivatives() {
    for law in laws() {
        let face = SurfactantFace {
            law,
            amount_mol: 0.17,
        };
        for area in [0.45, 1.0, 2.7] {
            let state = face.at_area(area).unwrap();
            close(
                state.free_energy_j,
                energy(law, area, face.amount_mol),
                2e-15,
            );
            close(state.concentration_mol_m2 * area, face.amount_mol, 3e-17);
            let da = area * 1e-5;
            let derivative = (energy(law, area - 2.0 * da, face.amount_mol)
                - 8.0 * energy(law, area - da, face.amount_mol)
                + 8.0 * energy(law, area + da, face.amount_mol)
                - energy(law, area + 2.0 * da, face.amount_mol))
                / (12.0 * da);
            close(state.tension_n_m, derivative, 4e-11);
            let dn = 1e-6;
            let chemical = (energy(law, area, face.amount_mol + dn)
                - energy(law, area, face.amount_mol - dn))
                / (2.0 * dn);
            close(state.chemical_potential_j_mol.unwrap(), chemical, 5e-10);
            let slope = (face.at_area(area + da).unwrap().tension_n_m
                - face.at_area(area - da).unwrap().tension_n_m)
                / (2.0 * da);
            close(state.dilatational_elasticity_n_m, area * slope, 3e-10);
        }
    }
}

#[test]
fn g0_material_advection_closes_free_energy_work_and_momentum_without_coefficient_supply() {
    let before = triangle();
    let after = [[0.7, -0.2, 0.4], [1.8, 0.3, 0.7], [0.3, 1.2, 0.1]];
    for law in laws() {
        let face = SurfactantFace {
            law,
            amount_mol: 0.17,
        };
        let response = triangle_surfactant_step(before, after, face).unwrap();
        let mut work = 0.0;
        let mut force = [0.0; 3];
        let mut torque = [0.0; 3];
        for i in 0..3 {
            let midpoint = core::array::from_fn(|d| 0.5 * before[i][d] + 0.5 * after[i][d]);
            let moment = cross(midpoint, response.forces_n[i]);
            work += dot(response.forces_n[i], sub(after[i], before[i]));
            for d in 0..3 {
                force[d] += response.forces_n[i][d];
                torque[d] += moment[d];
            }
        }
        let independent_change = energy(law, geometry(after).unwrap().area, face.amount_mol)
            - energy(law, geometry(before).unwrap().area, face.amount_mol);
        close(response.energy.energy_change_j, independent_change, 7e-16);
        close(work + independent_change, 0.0, 9e-16);
        assert_eq!(response.energy.coefficient_supply_j, 0.0);
        for d in 0..3 {
            close(force[d], 0.0, 2e-16);
            close(torque[d], 0.0, 3e-16);
        }
        let reversed = triangle_surfactant_step(
            [before[0], before[2], before[1]],
            [after[0], after[2], after[1]],
            face,
        )
        .unwrap();
        for (i, j) in [0, 2, 1].into_iter().enumerate() {
            for d in 0..3 {
                close(response.forces_n[i][d], reversed.forces_n[j][d], 2e-15);
            }
        }
    }
}

#[test]
fn g1_exact_endpoint_tangent_retains_dilution_at_equal_and_finitely_changed_area() {
    let before = triangle();
    let direction = [[0.1, 0.3, -0.2], [-0.4, 0.1, 0.5], [0.2, -0.3, 0.4]];
    for law in laws() {
        let face = SurfactantFace {
            law,
            amount_mol: 0.17,
        };
        for scale in [0.7, 1.0, 1.0 + 1e-11, 1.8] {
            let after = before.map(|p| p.map(|x| scale * x));
            let tangent = triangle_surfactant_step_tangent(before, after, face, direction).unwrap();
            let epsilon = 2e-5;
            let force = |factor: f64| {
                triangle_surfactant_step(
                    before,
                    core::array::from_fn(|i| {
                        core::array::from_fn(|d| after[i][d] + factor * epsilon * direction[i][d])
                    }),
                    face,
                )
                .unwrap()
                .forces_n
            };
            let minus2 = force(-2.0);
            let minus = force(-1.0);
            let plus = force(1.0);
            let plus2 = force(2.0);
            for i in 0..3 {
                for d in 0..3 {
                    let reference = (minus2[i][d] - 8.0 * minus[i][d] + 8.0 * plus[i][d]
                        - plus2[i][d])
                        / (12.0 * epsilon);
                    close(tangent[i][d], reference, 4e-10);
                }
            }
        }
        let state = face.at_area(1.0).unwrap();
        let (_, derivative) = area_secant(face, 1.0, 1.0, state, state).unwrap();
        close(derivative, 0.5 * state.dilatational_elasticity_n_m, 0.0);
        assert!(derivative > 0.0);
    }
}

#[test]
fn g3_vanishing_coverage_preserves_the_ideal_limit_and_zero_inventory_recovers_clean_forces() {
    let ideal = SurfactantLaw::IdealDilute {
        clean_tension_n_m: 1.4,
        thermal_energy_j_mol: 0.7,
        reference_concentration_mol_m2: 1e20,
    };
    let langmuir = SurfactantLaw::Langmuir {
        clean_tension_n_m: 1.4,
        thermal_energy_j_mol: 0.7,
        saturation_concentration_mol_m2: 1e20,
    };
    for area in [0.45, 1.0, 1.0 + 2e-13, 2.7] {
        let evaluate = |law| {
            let face = SurfactantFace {
                law,
                amount_mol: 0.17,
            };
            area_secant(
                face,
                1.0,
                area,
                face.at_area(1.0).unwrap(),
                face.at_area(area).unwrap(),
            )
            .unwrap()
        };
        let (a, da) = evaluate(ideal);
        let (b, db) = evaluate(langmuir);
        close(a, b, 3e-16);
        close(da, db, 3e-16);
    }
    let before = triangle();
    let after = [[0.1, -0.1, 0.2], [2.2, 0.3, 0.1], [-0.2, 0.9, 0.4]];
    let direction = [[0.1, 0.2, 0.3], [-0.2, 0.1, 0.4], [0.3, -0.4, 0.2]];
    for law in laws() {
        let face = SurfactantFace {
            law,
            amount_mol: 0.0,
        };
        assert_eq!(
            triangle_surfactant_step(before, after, face).unwrap(),
            triangle_surface_tension_step(before, after, 1.4, 1.4).unwrap()
        );
        assert_eq!(
            triangle_surfactant_step_tangent(before, after, face, direction).unwrap(),
            triangle_surface_tension_step_tangent(before, after, 1.4, 1.4, direction).unwrap()
        );
        assert_eq!(face.at_area(1.0).unwrap().chemical_potential_j_mol, None);
    }
}

#[test]
fn g1_nonuniform_inventory_generates_the_assembled_marangoni_line_force() {
    let q = [[0.0; 3], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
    let faces = [[0, 1, 2], [0, 2, 3]];
    let amounts = [0.05, 0.15];
    let law = laws()[0];
    let response = assemble_surfactant_step(
        &q,
        &q,
        &faces,
        SurfactantField::Uniform {
            law,
            amounts_mol: &amounts,
        },
    )
    .unwrap();
    let a = 1.4 - 0.7 * amounts[0] / 0.5;
    let b = 1.4 - 0.7 * amounts[1] / 0.5;
    let reference = [
        [a / 2.0, b / 2.0, 0.0],
        [-a / 2.0, a / 2.0, 0.0],
        [-b / 2.0, -a / 2.0, 0.0],
        [b / 2.0, -b / 2.0, 0.0],
    ];
    for (force, expected) in response.forces_n.iter().zip(reference) {
        for d in 0..3 {
            close(force[d], expected[d], 2e-16);
        }
    }
    let per_face = amounts.map(|amount_mol| SurfactantFace { law, amount_mol });
    let heterogeneous =
        assemble_surfactant_step(&q, &q, &faces, SurfactantField::PerFace(&per_face)).unwrap();
    assert_eq!(response, heterogeneous);
    close(response.energy.energy_change_j, 0.0, 0.0);
}

#[test]
fn g0_amount_exchange_energy_is_stable_at_empty_endpoints_and_reversible() {
    for law in laws() {
        for (before, after) in [(0.0, 0.17), (0.17, 0.0), (0.17, 0.17 + 1e-13), (0.17, 0.29)] {
            let face = SurfactantFace {
                law,
                amount_mol: before,
            };
            let change = face.energy_change_at_fixed_area(1.0, after).unwrap();
            close(
                change,
                energy(law, 1.0, after) - energy(law, 1.0, before),
                4e-16,
            );
            let reverse = SurfactantFace {
                law,
                amount_mol: after,
            }
            .energy_change_at_fixed_area(1.0, before)
            .unwrap();
            assert_eq!(change, -reverse);
        }
    }
}

#[test]
fn invalid_states_refuse_without_clipping_tension_or_packing() {
    for law in laws() {
        for n in [-1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                SurfactantFace { law, amount_mol: n }.at_area(1.0),
                Err(SurfactantError::InvalidAmount)
            );
        }
        for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                SurfactantFace {
                    law,
                    amount_mol: 0.1
                }
                .at_area(area),
                Err(SurfactantError::InvalidArea)
            );
        }
        assert!(
            SurfactantFace {
                law,
                amount_mol: 3.0
            }
            .at_area(1.0)
            .is_err()
        );
    }
    let negative = SurfactantFace {
        law: laws()[0],
        amount_mol: 2.01,
    };
    assert_eq!(negative.at_area(1.0), Err(SurfactantError::NegativeTension));
    let packed = SurfactantFace {
        law: laws()[1],
        amount_mol: 1.0,
    };
    assert_eq!(packed.at_area(1.0), Err(SurfactantError::PackingLimit));
    let bad = SurfactantLaw::IdealDilute {
        clean_tension_n_m: 1.0,
        thermal_energy_j_mol: 0.0,
        reference_concentration_mol_m2: 1.0,
    };
    assert_eq!(
        SurfactantFace {
            law: bad,
            amount_mol: 0.0
        }
        .at_area(1.0),
        Err(SurfactantError::InvalidLaw)
    );
    assert!(
        assemble_surfactant_step(
            &triangle(),
            &triangle(),
            &[[0, 1, 3]],
            SurfactantField::Uniform {
                law: laws()[0],
                amounts_mol: &[0.1]
            }
        )
        .is_err()
    );
    assert!(
        assemble_surfactant_step(
            &triangle(),
            &triangle(),
            &[[0, 1, 2]],
            SurfactantField::Uniform {
                law: laws()[0],
                amounts_mol: &[]
            }
        )
        .is_err()
    );
}
