//! Analytic 3D dielectric memory, total-current continuity and physical heat.
//! Synthetic material sources; the parent target owns the cube mesh fixture.

use super::{controls, cube, limits, near, near_vector};
use fs_em::{
    electrostatic::{ElectrostaticError, dispersion::DielectricFieldState},
    spatial::electrostatic::dispersion::SpatialDebyeDielectric,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
#[cfg(feature = "resolved-leaky-dielectric")]
use fs_material::field_storage::dielectric::leaky::ResolvedLeakyDielectric;
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{DebyePoleProperties, bulk::ResolvedDebyeMaterial},
    },
    state_point::MaterialPropertySelection,
};
use fs_qty::Time;

fn material_card(
    epsilon: f64,
    poles: &[(f64, f64)],
    sigma: Option<f64>,
) -> (MaterialCard, Vec<DebyePoleProperties>) {
    let mut claims = ClaimSet::new();
    let mut values = vec![(
        "epsilon_inf".to_owned(),
        ELECTRIC_PERMITTIVITY_DIMS,
        epsilon,
    )];
    let names = poles
        .iter()
        .enumerate()
        .map(|(i, &(delta, tau))| {
            let permittivity_increment = format!("delta_{i}");
            let relaxation_time = format!("tau_{i}");
            values.push((
                permittivity_increment.clone(),
                ELECTRIC_PERMITTIVITY_DIMS,
                delta,
            ));
            values.push((relaxation_time.clone(), Time::DIMS, tau));
            DebyePoleProperties {
                permittivity_increment,
                relaxation_time,
            }
        })
        .collect();
    if let Some(sigma) = sigma {
        values.push(("sigma_dc".into(), ELECTRICAL_CONDUCTIVITY_DIMS, sigma));
    }
    for (name, dims, value) in values {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic tetrahedral dielectric regression".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (card, names)
}
fn law(epsilon: f64, poles: &[(f64, f64)]) -> ResolvedDebyeMaterial {
    let (card, names) = material_card(epsilon, poles, None);
    ResolvedDebyeMaterial::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
    )
    .unwrap()
}
#[cfg(feature = "resolved-leaky-dielectric")]
fn leaky_law(epsilon: f64, poles: &[(f64, f64)], sigma: f64) -> ResolvedLeakyDielectric {
    let (card, names) = material_card(epsilon, poles, Some(sigma));
    ResolvedLeakyDielectric::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
        "sigma_dc",
    )
    .unwrap()
}
fn affine_positions() -> Vec<[f64; 3]> {
    vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 3.0],
    ]
}
fn affine_voltages() -> Vec<f64> {
    vec![3.0, 2.0, 7.0, 1.5]
}
fn affine(tet: [u32; 4]) -> SpatialDebyeDielectric {
    SpatialDebyeDielectric::new(
        affine_positions(),
        vec![tet],
        vec![law(2e-11, &[(3e-11, 0.2)])],
        affine_voltages(),
        vec![vec![[0.0; 3]]],
        &[0, 1, 2, 3],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap()
}
fn held(field: &SpatialDebyeDielectric) -> Vec<(usize, f64)> {
    field.voltage_v().iter().copied().enumerate().collect()
}
fn end_values(positions: &[[f64; 3]], ground_only: bool) -> Vec<(usize, f64)> {
    positions
        .iter()
        .enumerate()
        .filter_map(|(v, p)| (p[0] == 0.0 || (!ground_only && p[0] == 1.0)).then_some((v, p[0])))
        .collect()
}

#[test]
fn oblique_vector_memory_uses_physical_volume_exact_weak_moments_and_original_orientation() {
    let e = [1.0, -2.0, 0.5];
    let square = 5.25;
    for tet in [[0, 1, 2, 3], [1, 0, 2, 3]] {
        let initial = affine(tet);
        near_vector(initial.electric_field_v_m()[0], e);
        near(initial.cell_volumes_m3()[0], 1.0);
        near(initial.stored_energy_j(), 1e-11 * square);
        let prepared = initial.prepare_step(0.1, &mut || false).unwrap();
        near(prepared.time_step_s(), 0.1);
        let candidate = initial
            .step_candidate(&prepared, &held(&initial), controls(), &mut || false)
            .unwrap();
        let next = candidate.field();
        let step = candidate.step();
        near_vector(next.electric_field_v_m()[0], e);
        near_vector(next.polarization_c_m2()[0][0], e.map(|v| v * 1.2e-11));
        near_vector(next.displacement_c_m2()[0], e.map(|v| v * 3.2e-11));
        near_vector(
            step.polarization_current_density_a_m2[0],
            e.map(|v| v * 1.2e-10),
        );
        for (&actual, factor) in
            next.nodal_charge_c()
                .iter()
                .zip([1.0 / 6.0, -1.0, 1.0, -1.0 / 6.0])
        {
            near(actual, 3.2e-11 * factor);
        }
        near(step.stored_change_j, 0.24e-11 * square);
        near(step.polarization_heat_j, 0.96e-11 * square);
        near(step.supplied_work_j, 1.2e-11 * square);
        assert_eq!(step.joule_heat_j, 0.0);
        assert!(step.energy_residual_j.abs() < 1e-24);
        assert!(step.boundary_charge_residual_c.abs() < 1e-24);
        assert!(initial.same_problem(next));
        assert_eq!(initial.time_s(), 0.0);
        fn receiving_state<F: DielectricFieldState>(field: &F) {
            near(field.cell_volumes_m3().iter().sum(), 1.0);
            assert_eq!(field.conductivity_s_m(), &[0.0]);
            assert_eq!(field.conductivity_properties(), &[None]);
        }
        receiving_state(next);
    }
}

#[test]
fn heterogeneous_tetrahedra_redistribute_voltage_with_second_order_dielectric_relaxation() {
    fn run(steps: usize) -> f64 {
        let (positions, tets) = cube(&[0.0, 0.5, 1.0], 1, 1);
        let boundary = end_values(&positions, false);
        let vertices: Vec<_> = boundary.iter().map(|&(v, _)| v).collect();
        let voltage = positions
            .iter()
            .map(|p| {
                if p[0] <= 0.5 {
                    4.0 * p[0] / 3.0
                } else {
                    2.0 / 3.0 + 2.0 * (p[0] - 0.5) / 3.0
                }
            })
            .collect();
        let left = law(2e-11, &[(6e-11, 0.1)]);
        let right = law(4e-11, &[]);
        let materials: Vec<_> = (0..12)
            .map(|cell| {
                if cell < 6 {
                    left.clone()
                } else {
                    right.clone()
                }
            })
            .collect();
        let polarization = materials
            .iter()
            .map(|material| vec![[0.0; 3]; material.poles().len()])
            .collect();
        let mut state = SpatialDebyeDielectric::new(
            positions,
            tets,
            materials,
            voltage,
            polarization,
            &vertices,
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
        let prepared = state
            .prepare_step(0.1 / steps as f64, &mut || false)
            .unwrap();
        for _ in 0..steps {
            let (next, step) = state
                .step_candidate(&prepared, &boundary, controls(), &mut || false)
                .unwrap()
                .into_parts();
            assert!(step.maximum_free_charge_residual_c < 1e-23);
            assert!(step.energy_residual_j.abs() < 1e-23);
            assert!(step.cell_heat_j.iter().all(|q| *q >= 0.0));
            for field in &next.displacement_c_m2()[1..] {
                near_vector(*field, next.displacement_c_m2()[0]);
            }
            state = next;
        }
        state.electric_field_v_m()[0][0]
    }
    // Equal layer widths: tau_eff = tau (epsilon_inf,left + epsilon_right)
    // / (epsilon_static,left + epsilon_right) = 0.05 s.
    let exact = -2.0 / 3.0 - (2.0 / 3.0) * (-2.0_f64).exp();
    let coarse = (run(20) - exact).abs();
    let fine = (run(40) - exact).abs();
    assert!(fine < 5e-5);
    assert!(coarse > 3.8 * fine && coarse < 4.2 * fine);
}

#[test]
fn isolated_precharge_relaxes_polarization_without_free_charge_or_source_work() {
    let (positions, tets) = cube(&[0.0, 1.0], 1, 1);
    let boundary = end_values(&positions, true);
    let vertices: Vec<_> = boundary.iter().map(|&(v, _)| v).collect();
    let voltage = positions.iter().map(|p| p[0]).collect();
    let mut state = SpatialDebyeDielectric::new(
        positions,
        tets,
        vec![law(2e-11, &[(6e-11, 0.2)]); 6],
        voltage,
        vec![vec![[0.0; 3]]; 6],
        &vertices,
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let original_charge = state.nodal_charge_c().to_vec();
    let initial_energy = state.stored_energy_j();
    let prepared = state.prepare_step(0.005, &mut || false).unwrap();
    let mut heat = 0.0;
    for _ in 0..40 {
        let (next, step) = state
            .step_candidate(&prepared, &boundary, controls(), &mut || false)
            .unwrap()
            .into_parts();
        assert!(step.supplied_work_j.abs() < 1e-25);
        assert!(step.energy_residual_j.abs() < 1e-24);
        heat += step.heat_j;
        state = next;
    }
    let recurrence = (0.095_f64 / 0.105).powi(40);
    let expected = -0.25 - 0.75 * recurrence;
    for &field in state.electric_field_v_m() {
        near_vector(field, [expected, 0.0, 0.0]);
    }
    for (&charge, &initial) in state.nodal_charge_c().iter().zip(&original_charge) {
        assert!((charge - initial).abs() < 1e-23);
    }
    near(state.stored_energy_j() + heat, initial_energy);
}

#[test]
fn disconnected_components_retain_resolved_fields_below_the_common_voltage_ulp() {
    fn initial(biases: [f64; 2]) -> (SpatialDebyeDielectric, Vec<(usize, f64)>) {
        let (mut positions, mut tets) = cube(&[0.0, 1.0], 1, 1);
        let count = positions.len();
        let shifted: Vec<_> = positions.iter().map(|p| [p[0], p[1] + 3.0, p[2]]).collect();
        let extra: Vec<_> = tets
            .iter()
            .map(|tet| tet.map(|v| v + count as u32))
            .collect();
        positions.extend(shifted);
        tets.extend(extra);
        let voltage: Vec<_> = (0..positions.len())
            .map(|v| biases[usize::from(v >= count)])
            .collect();
        let boundary: Vec<_> = positions
            .iter()
            .enumerate()
            .filter_map(|(v, p)| (p[0] == 0.0).then_some((v, voltage[v])))
            .collect();
        let vertices: Vec<_> = boundary.iter().map(|&(v, _)| v).collect();
        let polarization = (0..12)
            .map(|cell| vec![[if cell < 6 { 2e-11 } else { -3e-11 }, 0.0, 0.0]])
            .collect();
        let state = SpatialDebyeDielectric::new(
            positions,
            tets,
            vec![law(2e-11, &[(6e-11, 0.2)]); 12],
            voltage,
            polarization,
            &vertices,
            0.0,
            limits(),
            controls(),
            &mut || false,
        )
        .unwrap();
        (state, boundary)
    }
    let (zero, zero_boundary) = initial([0.0, 0.0]);
    let (biased, biased_boundary) = initial([1e20, -2e20]);
    let zero_step = zero.prepare_step(0.05, &mut || false).unwrap();
    let biased_step = biased.prepare_step(0.05, &mut || false).unwrap();
    let a = zero
        .step_candidate(&zero_step, &zero_boundary, controls(), &mut || false)
        .unwrap();
    let b = biased
        .step_candidate(&biased_step, &biased_boundary, controls(), &mut || false)
        .unwrap();
    for (&actual, &expected) in b
        .field()
        .electric_field_v_m()
        .iter()
        .zip(a.field().electric_field_v_m())
    {
        near_vector(actual, expected);
    }
    assert!(b.field().electric_field_v_m()[0][0].abs() > 0.1);
    assert_eq!(b.field().voltage_v(), biased.voltage_v());
    near(b.field().stored_energy_j(), a.field().stored_energy_j());
    near(b.step().heat_j, a.step().heat_j);
    assert!(b.step().supplied_work_j.abs() < 1e-25);
}

#[test]
fn prepared_native_hodge_cancellation_foreign_operators_and_trace_changes_are_transactional() {
    let state = affine([0, 1, 2, 3]);
    let mut polls = 0;
    let prepared = state
        .prepare_step(0.1, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    for stop in 1..=polls {
        let mut seen = 0;
        assert!(matches!(
            state.prepare_step(0.1, &mut || {
                seen += 1;
                seen == stop
            }),
            Err(ElectrostaticError::Cancelled)
        ));
    }
    let unrelated = affine([0, 1, 2, 3]);
    assert!(!state.same_problem(&unrelated));
    assert!(matches!(
        unrelated.step_candidate(&prepared, &held(&unrelated), controls(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        state.step_candidate(&prepared, &held(&state)[..3], controls(), &mut || false),
        Err(ElectrostaticError::Invalid(_))
    ));
    assert!(matches!(
        state.step_candidate(&prepared, &held(&state), controls(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
    let replay = state
        .step_candidate(&prepared, &held(&state), controls(), &mut || false)
        .unwrap();
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.polarization_c_m2(), &[vec![[0.0; 3]]]);
    near(replay.field().time_s(), 0.1);
}

#[test]
fn point_edge_duplicate_and_unanchored_component_voltage_traces_are_refused() {
    for vertices in [vec![0], vec![0, 1], vec![0, 1, 1, 2]] {
        assert!(matches!(
            SpatialDebyeDielectric::new(
                affine_positions(),
                vec![[0, 1, 2, 3]],
                vec![law(2e-11, &[])],
                affine_voltages(),
                vec![vec![]],
                &vertices,
                0.0,
                limits(),
                controls(),
                &mut || false
            ),
            Err(ElectrostaticError::Invalid(_))
        ));
    }
    let mut positions = affine_positions();
    positions.extend(affine_positions().iter().map(|p| [p[0] + 5.0, p[1], p[2]]));
    assert!(matches!(
        SpatialDebyeDielectric::new(
            positions,
            vec![[0, 1, 2, 3], [4, 5, 6, 7]],
            vec![law(2e-11, &[]); 2],
            vec![0.0; 8],
            vec![vec![]; 2],
            &[0, 1, 2],
            0.0,
            limits(),
            controls(),
            &mut || false
        ),
        Err(ElectrostaticError::Invalid(_))
    ));
}

#[cfg(feature = "resolved-leaky-dielectric")]
#[test]
fn oblique_leakage_adds_source_current_and_separate_positive_joule_heat() {
    let state = SpatialDebyeDielectric::new_leaky(
        affine_positions(),
        vec![[0, 1, 2, 3]],
        vec![leaky_law(2e-11, &[(3e-11, 0.2)], 4e-11)],
        affine_voltages(),
        vec![vec![[0.0; 3]]],
        &[0, 1, 2, 3],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let prepared = state.prepare_step(0.1, &mut || false).unwrap();
    let candidate = state
        .step_candidate(&prepared, &held(&state), controls(), &mut || false)
        .unwrap();
    let step = candidate.step();
    near_vector(
        step.conduction_current_density_a_m2[0],
        [4e-11, -8e-11, 2e-11],
    );
    near(step.cell_polarization_heat_j[0], 0.96e-11 * 5.25);
    near(step.cell_joule_heat_j[0], 0.4e-11 * 5.25);
    near(step.cell_heat_j[0], 1.36e-11 * 5.25);
    near(step.supplied_work_j, 1.6e-11 * 5.25);
    assert!(step.energy_residual_j.abs() < 1e-24);
    assert_eq!(
        candidate.field().conductivity_properties(),
        &[Some("sigma_dc".into())]
    );
    for vertex in 0..4 {
        near(
            step.nodal_boundary_supplied_charge_c[vertex].unwrap(),
            step.nodal_boundary_charge_increment_c[vertex].unwrap()
                + step.nodal_conduction_outflow_c[vertex],
        );
    }
}

#[cfg(feature = "resolved-leaky-dielectric")]
#[test]
fn isolated_ohmic_precharge_decays_at_maxwell_time_and_storage_becomes_joule_heat() {
    let (positions, tets) = cube(&[0.0, 1.0], 1, 1);
    let boundary = end_values(&positions, true);
    let vertices: Vec<_> = boundary.iter().map(|&(v, _)| v).collect();
    let voltage = positions.iter().map(|p| p[0]).collect();
    let mut state = SpatialDebyeDielectric::new_leaky(
        positions,
        tets,
        vec![leaky_law(2e-11, &[], 4e-11); 6],
        voltage,
        vec![vec![]; 6],
        &vertices,
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let initial_energy = state.stored_energy_j();
    let prepared = state.prepare_step(0.02, &mut || false).unwrap();
    let mut heat = 0.0;
    for _ in 0..40 {
        let (next, step) = state
            .step_candidate(&prepared, &boundary, controls(), &mut || false)
            .unwrap()
            .into_parts();
        assert!(step.supplied_work_j.abs() < 1e-25);
        assert_eq!(step.polarization_heat_j, 0.0);
        assert!(step.joule_heat_j > 0.0);
        assert!(step.maximum_free_charge_residual_c < 1e-23);
        assert!(step.energy_residual_j.abs() < 1e-24);
        heat += step.joule_heat_j;
        state = next;
    }
    // tau = epsilon/sigma = 0.5 s. Midpoint gives this exact discrete decay.
    let expected = -(0.98_f64 / 1.02).powi(40);
    for &field in state.electric_field_v_m() {
        near_vector(field, [expected, 0.0, 0.0]);
    }
    near(state.stored_energy_j() + heat, initial_energy);
}
