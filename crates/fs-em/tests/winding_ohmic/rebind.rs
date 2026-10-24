//! G1/G3 canonical material replacement in an actual mutually coupled circuit.
use super::*;
use fs_em::winding::{circuit::WindingCircuit, ohmic::WindingMaterialRebindControls};
use fs_phs::circuit::{Branch, CircuitGraph};

fn rebind_controls() -> WindingMaterialRebindControls {
    WindingMaterialRebindControls {
        maximum_absolute_work_j: 1.0,
        energy_tolerance_j: 1e-10,
    }
}

fn custom_circuit<'a>(
    problem: &'a fs_em::planar::PlanarMagnetostatics,
    profiles: Vec<Vec<f64>>,
    boundary: &[usize],
    graph: CircuitGraph,
    branches: Vec<usize>,
) -> WindingCircuit<'a> {
    PlanarWindingBasis::new(
        problem,
        profiles,
        boundary,
        vec![0.002, 0.003],
        winding_config(),
        linear(),
        &mut || false,
    )
    .unwrap()
    .into_circuit(graph, branches, &mut || false)
    .unwrap()
}

#[test]
fn a_live_material_swap_preserves_flux_and_charge_then_drives_new_current_and_heat() {
    let resistivity = ElectricalPropertyBasis::Resistivity;
    let conductivity = ElectricalPropertyBasis::Conductivity;
    let original_cards = [
        card(2e-5, 0.0, 2e-4, resistivity),
        card(3e-5, 0.0, 2e-4, resistivity),
    ];
    let replacement_cards = [
        card(4e-5, 0.0, 4e-4, conductivity),
        card(1.5e-5, 0.0, 1e-4, conductivity),
    ];
    let (old_problem, profiles, boundary) = problem(&original_cards, resistivity, 0.5);
    let (new_problem, new_profiles, new_boundary) = problem(&replacement_cards, conductivity, 0.5);
    let old_magnetic = winding(&old_problem, profiles, &boundary);
    let new_magnetic = winding(&new_problem, new_profiles, &new_boundary);
    let initial = OhmicWindingCircuit::new(
        &old_magnetic,
        sources(&old_magnetic, resistivity),
        0.0,
        &[1.0, 0.25],
        &[0.2],
        &[1.4],
        &mut || false,
    )
    .unwrap();
    let (accepted, _) = initial
        .step_to_candidate(
            0.01,
            initial.conductors().to_vec(),
            &[1.4],
            controls(),
            &mut || false,
        )
        .unwrap()
        .into_parts();
    let original_state = accepted.state().to_vec();
    let (rebound, report) = accepted
        .rebind_material_candidate(
            &new_magnetic,
            sources(&new_magnetic, conductivity),
            &[1.4],
            rebind_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(rebound.time_s(), 0.01);
    assert_eq!(report.time_s, 0.01);
    let nd = old_magnetic.dae().system.differential_dim();
    assert_eq!(
        rebound.state()[..nd]
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>(),
        accepted.state()[..nd]
            .iter()
            .map(|v| v.to_bits())
            .collect::<Vec<_>>()
    );
    let lambda = [0, 1].map(|i| original_state[old_magnetic.flux_coordinates()[i]]);
    let l = new_magnetic.basis().inductance_h();
    assert!(l[1].abs() > 1e-6);
    let expected_current = two_by_two([l[0], l[1], l[2], l[3]], lambda);
    let old_current = accepted.currents_a().unwrap();
    for i in 0..2 {
        near(report.currents_a[1][i], expected_current[i]);
    }
    assert!((old_current[0] - expected_current[0]).abs() > 1e-3);
    let expected_work = 0.5
        * (lambda[0] * (expected_current[0] - old_current[0])
            + lambda[1] * (expected_current[1] - old_current[1]));
    near(report.supplied_work_j, expected_work);
    near(
        report.stored_energy_j[1] - report.stored_energy_j[0],
        expected_work,
    );
    assert!(report.energy_residual_j.abs() < rebind_controls().energy_tolerance_j);
    near(rebound.resistance_ohm()[0], 1.6);
    near(rebound.resistance_ohm()[1], 0.6);
    for (triangle, identity) in report.material_state_identities.iter().enumerate() {
        assert_eq!(identity[0], old_problem.materials()[triangle].identity().0);
        assert_eq!(identity[1], new_problem.materials()[triangle].identity().0);
        assert_ne!(identity[0], identity[1]);
    }
    for source in rebound.conductors() {
        assert_eq!(source.conductor.source_basis(), conductivity);
    }

    // The next original midpoint step consumes the new FULL inductance block
    // and spatial conductor resistance; the original capacitor is retained.
    let dt = 0.005;
    let candidate = rebound
        .step_to_candidate(
            0.015,
            rebound.conductors().to_vec(),
            &[1.4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let r = [
        1.0 + rebound.resistance_ohm()[0],
        3.0 + rebound.resistance_ohm()[1],
    ];
    let expected_next = two_by_two(
        [l[0] + dt * r[0] / 2.0, l[1], l[2], l[3] + dt * r[1] / 2.0],
        [
            lambda[0] - dt * r[0] * expected_current[0] / 2.0 + dt * 1.4,
            lambda[1] - dt * r[1] * expected_current[1] / 2.0,
        ],
    );
    for (got, want) in candidate
        .circuit()
        .currents_a()
        .unwrap()
        .iter()
        .zip(expected_next)
    {
        near(*got, want);
    }
    for cell in &candidate.step().cells {
        assert_eq!(
            cell.source.conductor.material().card_identity(),
            replacement_cards[cell.winding].content_hash()
        );
        let current = f64::midpoint(expected_current[cell.winding], expected_next[cell.winding]);
        let j =
            new_magnetic.basis().profiles_per_m2()[cell.winding][cell.source.triangle] * current;
        near(
            cell.heat_j,
            dt * cell.source.conductor.resistivity_ohm_m() * j * j * cell.volume_m3,
        );
    }
    near(
        candidate.step().stored_after_j + candidate.step().heat_j,
        report.stored_energy_j[0] + report.supplied_work_j + candidate.step().native.supplied,
    );
    let (reversed, reverse) = rebound
        .rebind_material_candidate(
            &old_magnetic,
            accepted.conductors().to_vec(),
            &[1.4],
            rebind_controls(),
            &mut || false,
        )
        .unwrap();
    near(reverse.supplied_work_j, -report.supplied_work_j);
    for (a, b) in reversed.currents_a().unwrap().iter().zip(&old_current) {
        near(*a, *b);
    }
    assert_eq!(accepted.state(), original_state);
    assert_eq!(accepted.time_s(), 0.01);
}

#[test]
fn resistance_only_replacement_has_zero_work_and_changes_subsequent_dissipation() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [card(2e-5, 0.0, 2e-4, basis), card(3e-5, 0.0, 2e-4, basis)];
    let new_cards = [card(4e-5, 0.0, 2e-4, basis), cards[1].clone()];
    let (p0, profiles, boundary) = problem(&cards, basis, 0.5);
    let (p1, profiles1, boundary1) = problem(&new_cards, basis, 0.5);
    let m0 = winding(&p0, profiles, &boundary);
    let m1 = winding(&p1, profiles1, &boundary1);
    let original = OhmicWindingCircuit::new(
        &m0,
        sources(&m0, basis),
        0.0,
        &[1.0, 0.0],
        &[0.2],
        &[1.4],
        &mut || false,
    )
    .unwrap();
    let mut admission = rebind_controls();
    admission.maximum_absolute_work_j = 0.0;
    let (rebound, report) = original
        .rebind_material_candidate(&m1, sources(&m1, basis), &[1.4], admission, &mut || false)
        .unwrap();
    assert_eq!(report.supplied_work_j, 0.0);
    assert_eq!(report.stored_energy_j[0], report.stored_energy_j[1]);
    assert_eq!(report.currents_a[0], report.currents_a[1]);
    let a = original
        .step_to_candidate(
            0.01,
            original.conductors().to_vec(),
            &[1.4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let b = rebound
        .step_to_candidate(
            0.01,
            rebound.conductors().to_vec(),
            &[1.4],
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(b.circuit().currents_a().unwrap()[0] < a.circuit().currents_a().unwrap()[0]);
    assert!(
        b.step()
            .cells
            .iter()
            .filter(|c| c.winding == 0)
            .map(|c| c.heat_j)
            .sum::<f64>()
            > a.step()
                .cells
                .iter()
                .filter(|c| c.winding == 0)
                .map(|c| c.heat_j)
                .sum::<f64>()
    );
}

#[test]
fn changed_geometry_boundary_winding_graph_work_and_late_cancellation_refuse_atomically() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [card(2e-5, 0.0, 2e-4, basis), card(3e-5, 0.0, 2e-4, basis)];
    let changed = [card(4e-5, 0.0, 4e-4, basis), card(3e-5, 0.0, 4e-4, basis)];
    let (p0, profiles, boundary) = problem(&cards, basis, 0.5);
    let (p1, profiles1, boundary1) = problem(&changed, basis, 0.5);
    let m0 = winding(&p0, profiles.clone(), &boundary);
    let m1 = winding(&p1, profiles1, &boundary1);
    let circuit = OhmicWindingCircuit::new(
        &m0,
        sources(&m0, basis),
        0.0,
        &[1.0, 0.25],
        &[0.2],
        &[1.4],
        &mut || false,
    )
    .unwrap();
    let original = circuit.state().to_vec();
    let mut denied = rebind_controls();
    denied.maximum_absolute_work_j = 0.0;
    assert!(
        circuit
            .rebind_material_candidate(&m1, sources(&m1, basis), &[1.4], denied, &mut || false)
            .is_err()
    );
    let (deep, dp, db) = problem(&cards, basis, 1.0);
    let deeper = winding(&deep, dp, &db);
    let bottom: Vec<_> = (0..=8).collect();
    let different_boundary = winding(&p0, profiles.clone(), &bottom);
    let different_turns = winding(
        &p0,
        profiles
            .iter()
            .map(|p| p.iter().map(|j| j * 2.0).collect())
            .collect(),
        &boundary,
    );
    let mut changed_graph = graph();
    changed_graph.branches[0].2 = Branch::Capacitor { farads: 0.8 };
    let different_graph =
        custom_circuit(&p0, profiles.clone(), &boundary, changed_graph, vec![4, 2]);
    for replacement in [
        &deeper,
        &different_boundary,
        &different_turns,
        &different_graph,
    ] {
        assert!(matches!(
            circuit.rebind_material_candidate(
                replacement,
                sources(replacement, basis),
                &[1.4],
                rebind_controls(),
                &mut || false
            ),
            Err(OhmicWindingError::Invalid(_))
        ));
    }
    let mut polls = 0;
    circuit
        .rebind_material_candidate(
            &m1,
            sources(&m1, basis),
            &[1.4],
            rebind_controls(),
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        circuit.rebind_material_candidate(
            &m1,
            sources(&m1, basis),
            &[1.4],
            rebind_controls(),
            &mut || {
                count += 1;
                count == polls
            }
        ),
        Err(OhmicWindingError::Cancelled)
    ));
    assert_eq!(circuit.state(), original);
    assert_eq!(circuit.time_s(), 0.0);
}

#[test]
fn incompatible_ideal_current_constraints_cannot_create_a_flux_impulse() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [card(2e-5, 0.0, 2e-4, basis), card(3e-5, 0.0, 2e-4, basis)];
    let changed = [card(2e-5, 0.0, 4e-4, basis), card(3e-5, 0.0, 4e-4, basis)];
    let (p0, profiles, boundary) = problem(&cards, basis, 0.5);
    let (p1, profiles1, boundary1) = problem(&changed, basis, 0.5);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (2, 0, Branch::Inductor { henries: 1.0 }),
            (1, 0, Branch::CurrentSource { port: 0 }),
            (2, 0, Branch::CurrentSource { port: 1 }),
        ],
        transformers: vec![],
    };
    let m0 = custom_circuit(&p0, profiles, &boundary, graph.clone(), vec![0, 1]);
    let m1 = custom_circuit(&p1, profiles1, &boundary1, graph, vec![0, 1]);
    let circuit = OhmicWindingCircuit::new(
        &m0,
        sources(&m0, basis),
        0.0,
        &[1.0, 0.25],
        &[],
        &[1.0, 0.25],
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        circuit.rebind_material_candidate(
            &m1,
            sources(&m1, basis),
            &[1.0, 0.25],
            rebind_controls(),
            &mut || false
        ),
        Err(OhmicWindingError::Circuit(_))
    ));
    assert_eq!(circuit.time_s(), 0.0);
}
