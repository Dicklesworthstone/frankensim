use super::*;

fn refreshed(
    state: &PlanarEddy<'_>,
    card: &MaterialCard,
    t: f64,
) -> (
    Vec<fs_material::state_point::ResolvedMaterialStatePoint>,
    Vec<EddyConductor>,
) {
    let materials: Vec<_> = state
        .problem()
        .materials()
        .iter()
        .map(|m| m.at_query_point(card, &point(t)).unwrap())
        .collect();
    let sources = state
        .conductors()
        .iter()
        .map(|c| EddyConductor {
            triangle: c.triangle,
            region: c.region,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                materials[c.triangle].clone(),
                c.conductor.source_basis(),
            )
            .unwrap(),
        })
        .collect();
    (materials, sources)
}

#[test]
fn changed_permeability_uses_both_endpoint_efforts_and_named_parameter_work() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let depth = 0.3;
    let dt = 0.1;
    for slope in [-0.0002, 0.0, 0.001] {
        let card = material::affine_magnetic_card(1.0 / 18.0, 0.0, 1.0, slope, basis);
        let problem = square(&card, basis, depth);
        let initial = problem
            .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
            .unwrap();
        let state = PlanarEddy::new(
            &initial,
            conductors(&problem, basis),
            &[],
            &boundary(0.0),
            0.0,
            controls(),
            &mut || false,
        )
        .unwrap();
        let (materials, sources) = refreshed(&state, &card, 500.0);
        let next = state
            .step_with_materials_candidate(
                dt,
                materials,
                sources,
                &[0.0; 4],
                limits(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let c = depth;
        let k0 = 4.0 * depth;
        let k1 = 4.0 * depth / (1.0 + 200.0 * slope);
        let a = (c / dt - k0 / 2.0) / (c / dt + k1 / 2.0);
        near(next.state().potential_wb_m()[4], a);
        near(next.step().parameter_work_j, 0.5 * a * (k1 - k0));
        near(next.step().stored_after_j, 0.5 * k1 * a * a);
        near(next.step().heat_j, c * (a - 1.0) * (a - 1.0) / dt);
        near(next.step().energy_residual_j, 0.0);
        near(
            next.step()
                .magnetic_cells
                .iter()
                .map(|c| c.parameter_work_j)
                .sum(),
            next.step().parameter_work_j,
        );
        // The accepted candidate OWNS the rebuilt operator; a second interval
        // continues that exact source/field history instead of resetting A.
        let (retained, _) = next.into_parts();
        let (materials, sources) = refreshed(&retained, &card, 600.0);
        let second = retained
            .step_with_materials_candidate(
                2.0 * dt,
                materials,
                sources,
                &[0.0; 4],
                limits(),
                controls(),
                &mut || false,
            )
            .unwrap();
        let k2 = 4.0 * depth / (1.0 + 300.0 * slope);
        near(
            second.state().potential_wb_m()[4],
            a * (c / dt - k1 / 2.0) / (c / dt + k2 / 2.0),
        );
        near(second.step().energy_residual_j, 0.0);
        near(state.potential_wb_m()[4], 1.0);
        assert_eq!(state.time_s(), 0.0);
    }
}

#[test]
fn insulating_algebraic_field_rebalances_without_a_fake_fixed_a_reset() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let source = material::affine_magnetic_card(1.0, 0.0, 1.0, 0.001, basis);
    let problem = square(&source, basis, 0.4);
    // Balanced current and a nonuniform fixed boundary produce an actual field
    // with an interior algebraic coordinate (C=0 everywhere).
    let current = [12.0, -4.0, -4.0, -4.0];
    let bc = vec![(0, 0.0), (1, 1.0), (2, 0.0), (3, 0.0)];
    let initial = problem
        .solve(&current, &bc, linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new_with_impressed_current(
        &initial,
        vec![],
        &[0, 1, 2, 3],
        &bc,
        &current,
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    // Spatially unequal temperatures change the harmonic interior trace.
    let temperatures = [400.0, 500.0, 600.0, 700.0];
    let materials: Vec<_> = problem
        .materials()
        .iter()
        .zip(temperatures)
        .map(|(m, t)| m.at_query_point(&source, &point(t)).unwrap())
        .collect();
    let reference = PlanarMagnetostatics::new(
        problem.mesh().clone(),
        materials.clone(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let equilibrium = reference
        .solve(&current, &bc, linear(), &mut || false)
        .unwrap();
    let next = state
        .step_with_materials_candidate(
            0.1,
            materials,
            vec![],
            &current,
            limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (a, b) in next
        .state()
        .potential_wb_m()
        .iter()
        .zip(equilibrium.potential_wb_m())
    {
        near(*a, *b);
    }
    assert!((next.state().potential_wb_m()[4] - state.potential_wb_m()[4]).abs() > 1e-4);
    near(next.step().heat_j, 0.0);
    near(next.step().energy_residual_j, 0.0);
}

#[test]
fn thermal_refresh_refusals_preserve_original_field_and_sources() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = material::affine_magnetic_card(1.0 / 18.0, 0.0, 1.0, 0.001, basis);
    let problem = square(&card, basis, 0.3);
    let initial = problem
        .solve(&[12.0; 4], &boundary(0.0), linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        conductors(&problem, basis),
        &[],
        &boundary(0.0),
        0.0,
        controls(),
        &mut || false,
    )
    .unwrap();
    let (materials, sources) = refreshed(&state, &card, 400.0);
    assert!(
        state
            .step_to_candidate(0.1, sources.clone(), &[0.0; 4], controls(), &mut || false)
            .is_err()
    );
    let mut tiny = limits();
    tiny.max_cells = 3;
    assert!(
        state
            .step_with_materials_candidate(
                0.1,
                materials.clone(),
                sources.clone(),
                &[0.0; 4],
                tiny,
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        state
            .step_with_materials_candidate(
                0.1,
                materials.clone(),
                sources.clone(),
                &[0.0; 4],
                limits(),
                controls(),
                &mut || true
            )
            .is_err()
    );
    let foreign = material::card(1.0 / 18.0, 0.0, 1.1, basis);
    assert!(
        state
            .step_with_materials_candidate(
                0.1,
                vec![bundle(&foreign, 400.0, basis); 4],
                sources,
                &[0.0; 4],
                limits(),
                controls(),
                &mut || false
            )
            .is_err()
    );
    assert_eq!(state.time_s(), 0.0);
    near(state.potential_wb_m()[4], 1.0);
    assert_eq!(state.problem().permeability_h_m(), &[1.0; 4]);
}
