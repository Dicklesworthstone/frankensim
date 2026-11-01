//! Source-bound graph resistance must change the native trajectory without
//! restarting the accepted field, terminal current, or other circuit storage.
use super::*;
#[path = "../winding_ohmic/common.rs"]
mod resistor_source;
use fs_em::winding::nonlinear::resistors::GraphResistorMaterial;
use fs_material::conductor::ElectricalPropertyBasis;

fn graph_sources(card: &MaterialCard, temperature: f64) -> Vec<GraphResistorMaterial> {
    let state = ResolvedConductorState::try_from_material_with_basis(
        resistor_source::bundle(card, temperature, ElectricalPropertyBasis::Resistivity),
        ElectricalPropertyBasis::Resistivity,
    )
    .unwrap();
    // Deliberately unordered: branch indices are physical identity.
    [(4, 3.0), (1, 3.5)]
        .into_iter()
        .map(|(branch, length)| GraphResistorMaterial {
            branch,
            conductor: state.clone().into_uniform_conductor(length, 1.0).unwrap(),
        })
        .collect()
}

#[test]
fn graph_resistivity_refresh_retains_live_inductive_history_and_changes_native_response() {
    let (magnetic, conductor) = material(false, 18.0);
    let problem = cube(magnetic);
    let profiles = vec![profiles(&problem).remove(0)];
    let end = 0.03;
    let inductance = 2.0 * solve(&problem, &profiles[0]).stored_energy_j() + end;
    let winding_r = resistance(&profiles[0], 18.0, 0.6);
    let card = resistor_source::card(0.2, 0.002, 2.0, ElectricalPropertyBasis::Resistivity);
    let mut start = initial(vec![0.0], vec![0.0]);
    start.remaining_fluxes_wb = vec![0.08];
    start.charges_c = vec![0.06];
    let unbound = SpatialNonlinearWindingCircuit::new(
        history(&problem, &current(&profiles, &[0.0]), vec![None; 48]),
        CircuitGraph {
            node_count: 4,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: 99.0 }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
                (3, 0, Branch::Inductor { henries: 0.4 }),
                (1, 3, Branch::Resistor { ohms: 99.0 }),
                (3, 0, Branch::Capacitor { farads: 0.3 }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![2],
            profiles_per_m2: profiles.clone(),
            end_inductance_h: vec![end],
        },
        sources(&profiles, &conductor, 0.6),
        start,
        config(),
        &mut || false,
    )
    .unwrap();
    let bound = unbound
        .with_graph_resistors(graph_sources(&card, 300.0), &mut || false)
        .unwrap();
    assert!(unbound.graph_resistor_materials().is_empty());
    let accepted = bound
        .step_to_candidate(0.04, &[[0.0; 3]; 27], &[1.3], controls(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let old_a = accepted.field().potential_wb().to_vec();
    let old_state = accepted.state().to_vec();
    let old_current = accepted.currents_a()[0];
    let old_linkage = accepted.flux_linkages_wb().to_vec();
    let old_phi = old_state[accepted.dae().flux_index[0]];
    let old_q = old_state[accepted.dae().charge_index[0]];
    assert!(old_current > 0.0 && old_a.iter().any(|&a| a != 0.0));
    let mut polls = 0;
    let hot = accepted
        .with_graph_resistors(graph_sources(&card, 500.0), &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert_eq!(hot.time_s().to_bits(), accepted.time_s().to_bits());
    assert_eq!(hot.field().potential_wb(), old_a);
    assert_eq!(hot.state(), old_state);
    assert_eq!(hot.currents_a(), accepted.currents_a());
    assert_eq!(hot.flux_linkages_wb(), old_linkage);
    assert_eq!(
        hot.dae().system.hamiltonian(hot.state()).to_bits(),
        accepted
            .dae()
            .system
            .hamiltonian(accepted.state())
            .to_bits()
    );

    let mut outcomes = vec![];
    let dt = 0.04;
    for (owner, temperature) in [(&accepted, 300.0), (&hot, 500.0)] {
        // Source curve and uniform geometry independently give these two R.
        let rho = 0.2 * (1.0 + 0.002 * (temperature - 300.0));
        let series_r = 3.5 * rho;
        let lc_r = 3.0 * rho;
        let expected_current =
            (1.3 + inductance * old_current / dt) / (winding_r + series_r + inductance / dt);
        // The extra native LC mode remains trapezoidal and is driven from the
        // ideal source, so it cannot change the winding's scalar BE equation.
        let (a, b, d) = (-dt / 0.6, dt / 0.8, 1.0 + dt / (0.6 * lc_r));
        let rhs = [
            old_phi + dt * old_q / 0.6,
            old_q + dt * 1.3 / lc_r - dt * old_q / (0.6 * lc_r) - dt * old_phi / 0.8,
        ];
        let determinant = d - a * b;
        let phi = (d * rhs[0] - a * rhs[1]) / determinant;
        let q = (rhs[1] - b * rhs[0]) / determinant;
        let lc_voltage = (q + old_q) / 0.6;
        let series_heat = dt * series_r * expected_current.powi(2);
        let lc_heat = dt * (1.3 - lc_voltage).powi(2) / lc_r;
        let candidate = owner
            .step_to_candidate(0.08, &[[0.0; 3]; 27], &[1.3], controls(), &mut || false)
            .unwrap();
        let next = candidate.circuit();
        let step = candidate.step();
        near(next.currents_a()[0], expected_current, 3e-7);
        near(
            next.flux_linkages_wb()[0],
            inductance * expected_current,
            3e-7,
        );
        near(next.state()[next.dae().flux_index[0]], phi, 3e-7);
        near(next.state()[next.dae().charge_index[0]], q, 3e-7);
        near(
            step.terminal_voltage_v[0],
            1.3 - series_r * expected_current,
            3e-7,
        );
        near(
            step.graph_resistors
                .iter()
                .find(|r| r.branch == 1)
                .unwrap()
                .heat_j,
            series_heat,
            3e-7,
        );
        near(
            step.graph_resistors
                .iter()
                .find(|r| r.branch == 4)
                .unwrap()
                .heat_j,
            lc_heat,
            3e-7,
        );
        near(
            step.heat_j,
            dt * winding_r * expected_current.powi(2) + series_heat + lc_heat,
            3e-7,
        );
        near(
            step.supplied_work_j,
            dt * 1.3 * (expected_current + (1.3 - lc_voltage) / lc_r),
            3e-7,
        );
        near(
            step.stored_before_j,
            0.5 * (inductance * old_current.powi(2) + old_phi.powi(2) / 0.4 + old_q.powi(2) / 0.3),
            3e-7,
        );
        near(
            step.stored_after_j,
            0.5 * (inductance * expected_current.powi(2) + phi.powi(2) / 0.4 + q.powi(2) / 0.3),
            3e-7,
        );
        near(
            step.numerical_dissipation_j,
            0.5 * inductance * (expected_current - old_current).powi(2),
            3e-7,
        );
        near(step.parameter_work_j, 0.0, 1e-14);
        near(step.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert_eq!(
            step.graph_materials.as_slice(),
            owner.graph_resistor_materials()
        );
        assert_eq!(
            step.graph_materials
                .iter()
                .map(|s| s.branch)
                .collect::<Vec<_>>(),
            [1, 4]
        );
        for source in &step.graph_materials {
            assert_eq!(
                source.conductor.material().card_identity(),
                card.content_hash()
            );
            assert_eq!(
                source.conductor.material().query().axes().get("T"),
                Some(&temperature)
            );
            assert_eq!(
                source.conductor.source_basis(),
                ElectricalPropertyBasis::Resistivity
            );
        }
        outcomes.push(next.currents_a()[0]);
    }
    assert!(outcomes[1] < outcomes[0]);

    let mut counter = 0;
    assert!(matches!(
        accepted.with_graph_resistors(graph_sources(&card, 500.0), &mut || {
            counter += 1;
            counter == polls
        }),
        Err(SpatialNonlinearWindingError::Cancelled)
    ));
    let mut wrong = graph_sources(&card, 500.0);
    wrong[0].branch = 2;
    assert!(matches!(
        accepted.with_graph_resistors(wrong, &mut || false),
        Err(SpatialNonlinearWindingError::Invalid(_))
    ));
    assert_eq!(accepted.field().potential_wb(), old_a);
    assert_eq!(accepted.state(), old_state);
    assert_eq!(accepted.flux_linkages_wb(), old_linkage);
    assert_eq!(accepted.time_s(), 0.04);
    assert!(
        accepted
            .graph_resistor_materials()
            .iter()
            .all(|source| source.conductor.material().query().axes().get("T") == Some(&300.0))
    );
}
