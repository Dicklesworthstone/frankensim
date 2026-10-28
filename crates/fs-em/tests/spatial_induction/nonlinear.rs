//! G1/G3/G4: original-source saturation, charge elimination and mechanical work.
//! Closed-form constants below come from exact integration on the unit cube;
//! neither reference trajectories nor dissipation use production matrices.
use super::*;
use fs_em::spatial::{
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        induction::nonlinear::{
            SpatialNonlinearInduction, SpatialNonlinearInductionControls,
            SpatialNonlinearInductionError,
        },
        nonlinear::{
            NonlinearVectorControls, SpatialNonlinearVectorField,
            SpatialNonlinearVectorMagnetostatics,
        },
    },
};
use fs_material::field_storage::nonlinear::{
    REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial,
};
use fs_qty::QuantitySpec;

fn saturation(sigma: f64) -> (SpatialNonlinearMaterial, ResolvedConductorState) {
    let mut claims = ClaimSet::new();
    for (name, dims, value, interpolation, validity) in [
        (
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            FIELD_STRENGTH_DIMS,
            PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots: vec![(0.0, 0.0), (1.0, 0.5), (2.0, 2.5), (4.0, 10.5)],
                dims: FIELD_STRENGTH_DIMS,
            },
            InterpolationPolicy::LinearInside,
            ValidityDomain::unconstrained().with_quantity(
                "B",
                QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
                0.0,
                4.0,
            ),
        ),
        (
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            PropertyValue::Scalar {
                value: sigma,
                dims: ELECTRICAL_CONDUCTIVITY_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            ValidityDomain::unconstrained(),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured three-dimensional nonlinear induction".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic soft conductor".into(),
            phase: "solid".into(),
            process: "reversible fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap();
    let law = ReversibleMagneticMaterial::resolve(
        &card,
        &query,
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &query,
        &[
            ScalarPropertyRequirement::try_new(
                REVERSIBLE_FIELD_STRENGTH_PROPERTY,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                ELECTRICAL_CONDUCTIVITY_PROPERTY,
                ELECTRICAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        SpatialNonlinearMaterial::Reversible(law),
        ResolvedConductorState::try_from_material(material).unwrap(),
    )
}
fn model(n: usize, material: SpatialNonlinearMaterial) -> SpatialNonlinearVectorMagnetostatics {
    let (positions, tets) = grid(n);
    let count = tets.len();
    SpatialNonlinearVectorMagnetostatics::new(
        positions,
        tets,
        vec![material; count],
        assembly(),
        &mut || false,
    )
    .unwrap()
}
fn nonlinear_controls() -> SpatialNonlinearInductionControls {
    SpatialNonlinearInductionControls {
        induction: controls(),
        maximum_newton_iterations: 40,
        maximum_backtracks: 35,
    }
}
fn boundary(p: &SpatialNonlinearVectorMagnetostatics, a: &[f64]) -> VectorMagneticBoundary {
    VectorMagneticBoundary {
        potential: p
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: p.face_edges(face).unwrap().map(|e| a[e]),
            })
            .collect(),
    }
}
fn solve<'m>(
    p: &'m SpatialNonlinearVectorMagnetostatics,
    a: &[f64],
    source: &[[f64; 3]],
) -> SpatialNonlinearVectorField<'m> {
    p.solve(
        source,
        &boundary(p, a),
        a,
        NonlinearVectorControls {
            field: static_controls(),
            relative_residual_tolerance: 2e-10,
            maximum_newton_iterations: 40,
            maximum_backtracks: 35,
        },
        &mut || false,
    )
    .unwrap()
}
fn mode(
    p: &SpatialNonlinearVectorMagnetostatics,
    amplitude: f64,
    current: f64,
) -> SpatialNonlinearVectorField<'_> {
    let mut a = vec![0.0; p.mesh().edges.len()];
    a[p.mesh().edges.binary_search(&[0, 7]).unwrap()] = amplitude;
    solve(p, &a, &vec![[current; 3]; 6])
}
fn history<'m>(
    field: &SpatialNonlinearVectorField<'m>,
    conductor: &ResolvedConductorState,
) -> SpatialNonlinearInduction<'m> {
    SpatialNonlinearInduction::new(
        field,
        vec![Some(conductor.clone()); field.problem().mesh().tets.len()],
        preparation(),
        2e-9,
        &mut || false,
    )
    .unwrap()
}
fn energy(a: f64) -> f64 {
    let b = 2.0 * a;
    if b <= 1.0 {
        0.25 * b * b
    } else if b <= 2.0 {
        0.25 + 0.5 * (b - 1.0) + (b - 1.0).powi(2)
    } else {
        1.75 + 2.5 * (b - 2.0) + 2.0 * (b - 2.0).powi(2)
    }
}

#[test]
fn linear_limit_obeys_implicit_decay_with_separate_numerical_loss() {
    let (material, conductor) = source(2.0, Some(5.0), None);
    let SpatialMagneticMaterial::Linear(material) = material else {
        unreachable!()
    };
    let p = model(1, SpatialNonlinearMaterial::Linear(material));
    let initial = mode(&p, 1.0, 4.0);
    let state = history(&initial, &conductor.unwrap());
    let (m, k, dt): (f64, f64, f64) = (3.0 * 5.0 / 40.0, 4.0 / 2.0, 0.04);
    let expected = m / (m + k * dt);
    let candidate = state
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            dt,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    let edge = p.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(candidate.state().potential_wb()[edge], expected, 2e-10);
    let step = candidate.step();
    near(step.heat_j, m * (expected - 1.0).powi(2) / dt, 2e-10);
    near(
        step.numerical_dissipation_j,
        0.5 * k * (expected - 1.0).powi(2),
        2e-10,
    );
    near(
        step.stored_change_j,
        0.5 * k * (expected.powi(2) - 1.0),
        2e-10,
    );
    near(step.impressed_work_j, 0.0, 1e-15);
    near(step.energy_residual_j, 0.0, 2e-10);
    assert!(step.maximum_charge_residual_a < 2e-8);
    assert!(step.maximum_ampere_residual_a < 2e-8);
    assert_eq!(state.time_s(), 0.0);
    near(candidate.state().time_s(), dt, 1e-15);
}

#[test]
fn fixed_recoil_retains_original_energy_reference_and_operating_bounds() {
    let (material, conductor) = source(2.0, Some(5.0), Some(0.2));
    let SpatialMagneticMaterial::Recoil(material) = material else {
        unreachable!()
    };
    let p = model(1, SpatialNonlinearMaterial::Recoil(material));
    let initial = mode(&p, 1.0, 4.0);
    let state = history(&initial, &conductor.unwrap());
    let (m, k, dt): (f64, f64, f64) = (0.375, 2.0, 0.04);
    let expected = m / (m + k * dt);
    let candidate = state
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            dt,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    let edge = p.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(candidate.state().potential_wb()[edge], expected, 2e-10);
    let step = candidate.step();
    // Constant Br has zero work against an interior-edge curl, but its source
    // reference 0.5*|Br|^2/mu remains in absolute storage throughout the step.
    near(step.stored_before_j, 1.01, 2e-10);
    near(step.stored_after_j, expected.powi(2) + 0.01, 2e-10);
    near(step.stored_change_j, expected.powi(2) - 1.0, 2e-10);
    near(step.heat_j, m * (expected - 1.0).powi(2) / dt, 2e-10);
    near(
        step.numerical_dissipation_j,
        (expected - 1.0).powi(2),
        2e-10,
    );
    near(step.energy_residual_j, 0.0, 2e-10);
    near(
        candidate
            .state()
            .field(&mut || false)
            .unwrap()
            .stored_energy_j(),
        expected.powi(2) + 0.01,
        2e-10,
    );
    assert!(
        state
            .step_candidate(
                &[[1e6; 3]; 6],
                &[[0.0; 3]; 8],
                1.0,
                nonlinear_controls(),
                &mut || false,
            )
            .is_err()
    );
    near(state.potential_wb()[edge], 1.0, 2e-10);
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn saturated_mode_matches_source_exact_solution_and_first_order_refinement() {
    // |curl(W)|=2 in every tetrahedron, total volume=1. Thus f(a)=2 H(2a),
    // j=J_amplitude/2 and charge-eliminated mass M=3*sigma/40.
    let (material, conductor) = saturation(5.0);
    let p = model(1, material);
    let initial = mode(&p, 0.8, 6.8);
    let state = history(&initial, &conductor);
    let edge = p.mesh().edges.binary_search(&[0, 7]).unwrap();
    let (m, dt): (f64, f64) = (0.375, 0.02);
    let expected = (m * 0.8 + dt * (1.0 + 3.0)) / (m + 8.0 * dt);
    let candidate = state
        .step_candidate(
            &[[2.0; 3]; 6],
            &[[0.0; 3]; 8],
            dt,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    let step = candidate.step();
    near(candidate.state().potential_wb()[edge], expected, 2e-10);
    near(step.stored_before_j, energy(0.8), 2e-10);
    near(step.stored_after_j, energy(expected), 2e-10);
    near(step.heat_j, m * (expected - 0.8).powi(2) / dt, 2e-10);
    near(
        step.numerical_dissipation_j,
        4.0 * (expected - 0.8).powi(2),
        2e-10,
    );
    near(step.impressed_work_j, expected - 0.8, 2e-10);
    near(step.energy_residual_j, 0.0, 2e-10);
    for h in &step.endpoint_field_strength_a_m {
        near(
            h.iter().fold(0.0_f64, |n, v| n.hypot(*v)),
            4.0 * expected - 1.5,
            2e-9,
        );
    }
    // m*a'+8*a-3=0 remains entirely on the second source segment.
    // Backward Euler must converge at order one, with every step independently
    // closing physical heat and numerical dissipation.
    let total: f64 = 0.04;
    let exact = 0.375 + (0.8 - 0.375) * (-8.0 * total / m).exp();
    let mut errors = vec![];
    for steps in [16, 32, 64] {
        let mut state = history(&initial, &conductor);
        for _ in 0..steps {
            let candidate = state
                .step_candidate(
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    total / steps as f64,
                    nonlinear_controls(),
                    &mut || false,
                )
                .unwrap();
            assert!(candidate.step().numerical_dissipation_j >= 0.0);
            near(candidate.step().energy_residual_j, 0.0, 2e-10);
            state = candidate.into_parts().0;
        }
        errors.push((state.potential_wb()[edge] - exact).abs());
    }
    for pair in errors.windows(2) {
        let order = (pair[0] / pair[1]).log2();
        assert!((0.9..1.1).contains(&order), "{errors:?}, order={order}");
    }
}

#[test]
fn crossing_saturation_knots_preserves_exact_energy_without_false_heat() {
    let (material, conductor) = saturation(5.0);
    let p = model(1, material);
    let initial = mode(&p, 1.1, 13.2);
    let state = history(&initial, &conductor);
    let candidate = state
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.5,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    let edge = p.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(candidate.state().potential_wb()[edge], 0.3, 2e-10);
    let step = candidate.step();
    near(step.stored_before_j, 2.33, 2e-10);
    near(step.stored_after_j, 0.09, 2e-10);
    near(step.stored_change_j, -2.24, 2e-10);
    near(step.heat_j, 0.48, 2e-10);
    near(step.numerical_dissipation_j, 1.76, 2e-10);
    near(step.energy_residual_j, 0.0, 2e-10);
    assert!(step.report.newton_iterations >= 2);
    near(step.cells.iter().map(|c| c.heat_j).sum(), 0.48, 2e-10);
    near(
        step.magnetic_cells
            .iter()
            .map(|c| c.numerical_dissipation_j)
            .sum(),
        1.76,
        2e-10,
    );
}

#[test]
fn nonlinear_radial_motion_matches_charge_and_lorentz_work_oracles() {
    let (material, conductor) = saturation(5.0);
    let p = model(1, material);
    let initial = mode(&p, 0.8, 6.8);
    let state = history(&initial, &conductor);
    let (alpha, m, dt): (f64, f64, f64) = (0.6, 0.375, 0.03);
    let velocity: Vec<_> = p
        .positions_m()
        .iter()
        .map(|x| x.map(|v| alpha * v))
        .collect();
    // Projection of v cross curl(W) is -2*alpha*(W-grad(q)). The exact
    // nonlinear endpoint equation is m*a'+(8+2*alpha*m)*a-3=j.
    let expected = (m * 0.8 + dt * (3.4 + 3.0)) / (m + dt * (8.0 + 2.0 * alpha * m));
    let rate = (expected - 0.8) / dt;
    let conductor_rate = rate + 2.0 * alpha * expected;
    let candidate = state
        .step_candidate(
            &[[6.8; 3]; 6],
            &velocity,
            dt,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    let edge = p.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(candidate.state().potential_wb()[edge], expected, 2e-10);
    let step = candidate.step();
    near(step.heat_j, m * conductor_rate.powi(2) * dt, 2e-10);
    near(
        step.mechanical_work_j,
        -2.0 * alpha * expected * m * conductor_rate * dt,
        2e-10,
    );
    near(step.impressed_work_j, 3.4 * (expected - 0.8), 2e-10);
    near(step.stored_change_j, energy(expected) - energy(0.8), 2e-10);
    near(
        step.numerical_dissipation_j,
        4.0 * (expected - 0.8).powi(2),
        2e-10,
    );
    near(step.energy_residual_j, 0.0, 2e-10);
    assert!(step.heat_j > 0.0 && step.mechanical_work_j.abs() > 1e-4);
    let q = |v| {
        if v == 0 {
            0.0
        } else if v == 7 {
            0.5
        } else {
            0.25
        }
    };
    let r = |v| {
        if v == 0 {
            0.0
        } else if v == 7 {
            1.0
        } else {
            -0.5
        }
    };
    for cell in &step.cells {
        let tet = p.mesh().tets[cell.cell];
        for i in 0..4 {
            let phi = -rate * (q(tet[i]) - q(tet[0])) + alpha * expected * (r(tet[i]) - r(tet[0]));
            near(cell.potential_v[i] - cell.potential_v[0], phi, 2e-9);
        }
    }
}

#[test]
fn rotated_saturated_translation_is_cancelled_by_charge_in_both_gauges() {
    let (material, conductor) = saturation(5.0);
    let p = model(2, material);
    let b = [0.8, -0.4, 1.2];
    let velocity = vec![[0.7, -0.3, 0.2]; p.mesh().vertex_count];
    let zero = vec![[0.0; 3]; p.mesh().tets.len()];
    let gradient = cross(velocity[0], b);
    let mut outcomes = vec![];
    for gauge in [false, true] {
        let psi = |v: u32| {
            if gauge {
                ((v * 7) % 11) as f64 * 0.17
            } else {
                0.0
            }
        };
        let a: Vec<_> = p
            .mesh()
            .edges
            .iter()
            .map(|&[u, v]| {
                let x = p.positions_m()[u as usize];
                let y = p.positions_m()[v as usize];
                let mid = std::array::from_fn(|j| 0.5 * (x[j] + y[j]));
                let vector = cross(b, mid);
                0.5 * (0..3).map(|j| vector[j] * (y[j] - x[j])).sum::<f64>() + psi(v) - psi(u)
            })
            .collect();
        let initial = solve(&p, &a, &zero);
        let state = history(&initial, &conductor);
        let candidate = state
            .step_candidate(&zero, &velocity, 0.1, nonlinear_controls(), &mut || false)
            .unwrap();
        let step = candidate.step();
        near(step.heat_j, 0.0, 1e-17);
        near(step.numerical_dissipation_j, 0.0, 1e-17);
        near(step.stored_change_j, 0.0, 2e-9);
        for field in &step.endpoint_flux_density_t {
            for j in 0..3 {
                near(field[j], b[j], 2e-9);
            }
        }
        for cell in &step.cells {
            let vertices = p.mesh().tets[cell.cell];
            for i in 0..4 {
                let x = p.positions_m()[vertices[i] as usize];
                let y = p.positions_m()[vertices[0] as usize];
                let phi = (0..3).map(|j| gradient[j] * (x[j] - y[j])).sum::<f64>();
                near(cell.potential_v[i] - cell.potential_v[0], phi, 2e-9);
                for j in 0..3 {
                    near(cell.current_density_a_m2[i][j], 0.0, 2e-8);
                }
            }
        }
        outcomes.push(candidate.state().field(&mut || false).unwrap());
    }
    for (a, b) in outcomes[0]
        .flux_density_t()
        .iter()
        .zip(outcomes[1].flux_density_t())
    {
        for j in 0..3 {
            near(a[j], b[j], 2e-9);
        }
    }
}

#[test]
fn source_bounds_real_work_exhaustion_and_late_cancellation_leave_history_atomic() {
    let (material, conductor) = saturation(5.0);
    let p = model(1, material);
    let initial = mode(&p, 0.8, 6.8);
    let state = history(&initial, &conductor);
    let old = state.potential_wb().to_vec();
    let (_, other) = saturation(6.0);
    assert!(
        SpatialNonlinearInduction::new(
            &initial,
            vec![Some(other); 6],
            preparation(),
            2e-9,
            &mut || false
        )
        .is_err()
    );
    let mut small = preparation();
    small.maximum_assembly_terms = 1;
    assert!(
        SpatialNonlinearInduction::new(
            &initial,
            vec![Some(conductor); 6],
            small,
            2e-9,
            &mut || false
        )
        .is_err()
    );
    let mut no_newton = nonlinear_controls();
    no_newton.maximum_newton_iterations = 0;
    assert!(matches!(
        state.step_candidate(&[[0.0; 3]; 6], &[[0.0; 3]; 8], 0.1, no_newton, &mut || {
            false
        }),
        Err(SpatialNonlinearInductionError::Nonconvergence(_))
    ));
    let mut overflow = nonlinear_controls();
    overflow.maximum_backtracks = usize::MAX;
    assert!(
        state
            .step_candidate(&[[0.0; 3]; 6], &[[0.0; 3]; 8], 0.1, overflow, &mut || false,)
            .is_err()
    );
    let mut no_krylov = nonlinear_controls();
    no_krylov.induction.restart = 1;
    no_krylov.induction.maximum_cycles = 1;
    assert!(
        state
            .step_candidate(&[[0.0; 3]; 6], &[[0.0; 3]; 8], 0.1, no_krylov, &mut || {
                false
            })
            .is_err()
    );
    assert!(
        state
            .step_candidate(
                &[[1e6; 3]; 6],
                &[[0.0; 3]; 8],
                1.0,
                nonlinear_controls(),
                &mut || false
            )
            .is_err()
    );
    let mut small_pe = nonlinear_controls();
    small_pe.induction.maximum_cell_magnetic_peclet = 0.01;
    assert!(
        state
            .step_candidate(
                &[[0.0; 3]; 6],
                &[[1.0, 0.0, 0.0]; 8],
                0.1,
                small_pe,
                &mut || false
            )
            .is_err()
    );
    let mut calls = 0;
    let success = state
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            nonlinear_controls(),
            &mut || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 50);
    let mut counter = 0;
    assert!(matches!(
        state.step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            nonlinear_controls(),
            &mut || {
                counter += 1;
                counter == calls - 1
            }
        ),
        Err(SpatialNonlinearInductionError::Spatial(
            SpatialMagneticError::Cancelled
        ))
    ));
    assert_eq!(state.potential_wb(), old);
    assert_eq!(state.time_s(), 0.0);
    let repeat = state
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            0.1,
            nonlinear_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        success.state().potential_wb(),
        repeat.state().potential_wb()
    );
    assert_eq!(
        success.step().heat_j.to_bits(),
        repeat.step().heat_j.to_bits()
    );
    assert_eq!(
        success.step().numerical_dissipation_j.to_bits(),
        repeat.step().numerical_dissipation_j.to_bits()
    );
}
