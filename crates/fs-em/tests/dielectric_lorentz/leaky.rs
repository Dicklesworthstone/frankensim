//! G1/G3/G4: original charge continuity with simultaneous inertial polarization.

use super::*;
use fs_material::field_storage::dielectric::lorentz::leaky::ResolvedLeakyLorentzMaterial;

fn leaky_law(poles: &[(f64, f64, f64)], sigma: f64) -> ResolvedLeakyLorentzMaterial {
    let (card, names) = source(poles, Some(sigma));
    ResolvedLeakyLorentzMaterial::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names,
        "sigma_dc",
    )
    .unwrap()
}

fn field(
    material: ResolvedLeakyLorentzMaterial,
    state: LorentzState<2>,
    traces: &[usize],
) -> PlanarLorentzDielectric {
    PlanarLorentzDielectric::new_leaky(
        layered(0.5, 1.0, 0.2),
        vec![material; 4],
        vec![0.0, 1.0, 2.0, 0.0, 1.0, 2.0],
        vec![vec![state]; 4],
        traces,
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn isolated_leaky_resonance_matches_original_three_state_odes_and_continuous_roots() {
    // E'=-5E-j, p'=j, j'=12E-10p-3j, with p=P/epsilon_inf,
    // j=J/epsilon_inf. Its exact roots are -2 and -3 +/- 4i.
    // The independent closed-form field below uses E(0)=1, E'(0)=-5,
    // E''(0)=13; no condensed algorithm enters this reference.
    let duration: f64 = 0.7;
    let exact_ratio = 8.0 / 17.0 * (-2.0 * duration).exp()
        + (-3.0 * duration).exp()
            * (9.0 / 17.0 * (4.0 * duration).cos() - 21.0 / 34.0 * (4.0 * duration).sin());
    let mut errors = Vec::new();
    for intervals in [24, 48, 96] {
        let dt = duration / intervals as f64;
        let sigma = 5.0 * EPS;
        let delta = 1.2 * EPS;
        let mut accepted = field(
            leaky_law(&[(delta, 10.0_f64.sqrt(), 3.0)], sigma),
            LorentzState::default(),
            &[0, 3],
        );
        let e0 = accepted.electric_field_v_m()[0][0];
        let energy0 = accepted.stored_energy_j();
        let charge0 = accepted.nodal_charge_c().to_vec();
        let prepared = accepted.prepare_step(dt, &mut || false).unwrap();
        let (mut heat, mut polarization_heat, mut joule_heat) = (0.0, 0.0, 0.0);
        for _ in 0..intervals {
            let (next, report) = accepted
                .step_candidate(&prepared, &[(0, 0.0), (3, 0.0)], controls(), &mut || false)
                .unwrap()
                .into_parts();
            for cell in 0..4 {
                let before = accepted.polarization_states()[cell][0];
                let after = next.polarization_states()[cell][0];
                let e = 0.5
                    * (accepted.electric_field_v_m()[cell][0] + next.electric_field_v_m()[cell][0]);
                let p = 0.5 * (before.polarization_c_m2[0] + after.polarization_c_m2[0]);
                let j = 0.5 * (before.current_a_m2[0] + after.current_a_m2[0]);
                let dp = after.polarization_c_m2[0] - before.polarization_c_m2[0];
                let dj = after.current_a_m2[0] - before.current_a_m2[0];
                assert!((dp - dt * j).abs() < 1e-23);
                assert!((dj - dt * (10.0 * (delta * e - p) - 3.0 * j)).abs() < 1e-22);
                let de =
                    next.electric_field_v_m()[cell][0] - accepted.electric_field_v_m()[cell][0];
                assert!((EPS * de + dp + dt * sigma * e).abs() < 1e-23);
                near(
                    report.cell_joule_heat_j[cell],
                    accepted.cell_volumes_m3()[cell] * dt * sigma * e * e,
                    2e-10,
                );
            }
            for vertex in [1, 2, 4, 5] {
                let balance = next.nodal_charge_c()[vertex] - accepted.nodal_charge_c()[vertex]
                    + report.nodal_conduction_outflow_c[vertex];
                assert!(balance.abs() < 1e-23);
            }
            assert_eq!(report.supplied_work_j, 0.0);
            assert!(report.energy_residual_j.abs() < 1e-23);
            heat += report.heat_j;
            polarization_heat += report.polarization_heat_j;
            joule_heat += report.joule_heat_j;
            accepted = next;
        }
        assert!(polarization_heat > 0.0 && joule_heat > 0.0);
        near(accepted.stored_energy_j() + heat, energy0, 2e-10);
        assert!((accepted.nodal_charge_c()[2] - charge0[2]).abs() > 1e-13);
        errors.push((accepted.electric_field_v_m()[0][0] / e0 - exact_ratio).abs());
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.7..4.3).contains(&ratio), "second-order ratio {ratio}");
    }
}

#[test]
fn equilibrium_dc_current_heats_without_changing_resonant_storage() {
    let sigma = 4e-10;
    let accepted = field(
        leaky_law(&[(DELTA, OMEGA, 2.0)], sigma),
        LorentzState {
            polarization_c_m2: [-2.0 * DELTA, 0.0],
            current_a_m2: [0.0; 2],
        },
        &[0, 3, 2, 5],
    );
    assert_eq!(accepted.conductivity_s_m(), &[sigma; 4]);
    assert!(
        accepted
            .conductivity_properties()
            .iter()
            .all(|name| name.as_deref() == Some("sigma_dc"))
    );
    let dt = 0.2;
    let prepared = accepted.prepare_step(dt, &mut || false).unwrap();
    let proposal = accepted
        .step_candidate(&prepared, &ends(2.0), controls(), &mut || false)
        .unwrap();
    let report = proposal.step();
    let expected = dt * sigma * 4.0 * 0.2;
    near(report.joule_heat_j, expected, 2e-12);
    near(report.supplied_conduction_work_j, expected, 2e-12);
    assert!(report.polarization_heat_j < 1e-34);
    assert!(report.stored_change_j.abs() < 1e-23);
    near(
        proposal.field().stored_energy_j(),
        accepted.stored_energy_j(),
        2e-12,
    );
    for current in &report.conduction_current_density_a_m2 {
        near(current[0], -2.0 * sigma, 2e-12);
        assert!(current[1].abs() < 1e-23);
    }
    assert!(report.energy_residual_j.abs() < 1e-23);
}

#[test]
fn revolved_and_oblique_tetrahedral_current_use_original_physical_volumes() {
    let sigma = 4e-10;
    let dt = 0.03;
    let angle = 1.3;
    let mesh = TriComplex2::from_indexed_triangles(
        "conducting resonant meridian",
        vec![[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]],
        vec![[0, 1, 2], [0, 2, 3]],
        Metric2::axisymmetric(angle).unwrap(),
    )
    .unwrap();
    let axial = AxisymmetricLorentzDielectric::new_leaky(
        mesh,
        vec![leaky_law(&[(DELTA, OMEGA, 2.0)], sigma); 2],
        vec![0.0, 0.0, 2.0, 2.0],
        vec![vec![LorentzState::default()]; 2],
        &[0, 1, 2, 3],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let prepared = axial.prepare_step(dt, &mut || false).unwrap();
    let proposal = axial
        .step_candidate(
            &prepared,
            &[(0, 0.0), (1, 0.0), (2, 2.0), (3, 2.0)],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(
        proposal.step().joule_heat_j,
        dt * sigma * 4.0 * 1.5 * angle,
        2e-12,
    );
    assert!(proposal.step().polarization_heat_j > 0.0);
    assert!(proposal.step().energy_residual_j.abs() < 1e-23);

    let spatial = SpatialLorentzDielectric::new_leaky(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ],
        vec![[0, 1, 2, 3]],
        vec![leaky_law(&[(DELTA, OMEGA, 2.0)], sigma)],
        vec![3.0, 2.0, 7.0, 1.5],
        vec![vec![LorentzState::default()]],
        &[0, 1, 2, 3],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let prepared = spatial.prepare_step(dt, &mut || false).unwrap();
    let traces: Vec<_> = spatial.voltage_v().iter().copied().enumerate().collect();
    let proposal = spatial
        .step_candidate(&prepared, &traces, controls(), &mut || false)
        .unwrap();
    near(spatial.cell_volumes_m3()[0], 1.0, 1e-13);
    near(proposal.step().joule_heat_j, dt * sigma * 5.25, 2e-12);
    for (axis, e) in [1.0, -2.0, 0.5].into_iter().enumerate() {
        near(
            proposal.step().conduction_current_density_a_m2[0][axis],
            sigma * e,
            2e-12,
        );
    }
    assert!(proposal.step().polarization_heat_j > 0.0);
    assert!(proposal.step().energy_residual_j.abs() < 1e-23);
}

#[test]
fn sourced_zero_conduction_matches_insulation_and_refusals_preserve_inertial_history() {
    let material = leaky_law(&[(DELTA, OMEGA, 2.0)], 0.0);
    let original = field(material.clone(), LorentzState::default(), &[0, 3, 2, 5]);
    let mut leaking = original.clone();
    let mut insulating = PlanarLorentzDielectric::new(
        original.mesh().clone(),
        vec![material.dielectric().clone(); 4],
        original.voltage_v().to_vec(),
        original.polarization_states().to_vec(),
        &[0, 3, 2, 5],
        0.0,
        limits(),
        controls(),
        &mut || false,
    )
    .unwrap();
    let a = leaking.prepare_step(0.03, &mut || false).unwrap();
    let b = insulating.prepare_step(0.03, &mut || false).unwrap();
    assert!(
        leaking
            .step_candidate(&b, &ends(2.0), controls(), &mut || false)
            .is_err()
    );
    assert!(matches!(
        leaking.step_candidate(&a, &ends(2.0), controls(), &mut || true),
        Err(ElectrostaticError::Cancelled)
    ));
    assert_eq!(leaking.time_s(), 0.0);
    assert_eq!(
        leaking.polarization_states(),
        original.polarization_states()
    );
    for step in 0..8 {
        let endpoint = ends(2.0 * (0.3 * step as f64).cos());
        let (next_a, report_a) = leaking
            .step_candidate(&a, &endpoint, controls(), &mut || false)
            .unwrap()
            .into_parts();
        let (next_b, report_b) = insulating
            .step_candidate(&b, &endpoint, controls(), &mut || false)
            .unwrap()
            .into_parts();
        assert_eq!(next_a.polarization_states(), next_b.polarization_states());
        assert_eq!(next_a.voltage_v(), next_b.voltage_v());
        assert_eq!(next_a.nodal_charge_c(), next_b.nodal_charge_c());
        assert_eq!(report_a.joule_heat_j, 0.0);
        assert_eq!(report_a.polarization_heat_j, report_b.polarization_heat_j);
        assert_eq!(report_a.supplied_work_j, report_b.supplied_work_j);
        leaking = next_a;
        insulating = next_b;
    }
}
