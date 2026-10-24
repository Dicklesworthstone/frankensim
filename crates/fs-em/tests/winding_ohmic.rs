#[path = "winding_ohmic/common.rs"]
mod common;
#[path = "winding_ohmic/rebind.rs"]
mod rebind;
use common::*;
use fs_em::winding::{PlanarWindingBasis, ohmic::{OhmicWindingCircuit, OhmicWindingError}};
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductorState};

#[test]
fn material_resistance_changes_actual_mutual_circuit_currents_and_complete_heat() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let mut outcomes = Vec::new();
    for scale in [1.0, 2.0] {
        let cards = [card(2e-5 * scale, 0.0, 2e-4, basis), card(3e-5, 0.0, 2e-4, basis)];
        let (problem, profiles, boundary) = problem(&cards, basis, 0.5);
        let magnetic = winding(&problem, profiles, &boundary);
        let circuit = OhmicWindingCircuit::new(&magnetic, sources(&magnetic, basis),
            0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).unwrap();
        near(circuit.resistance_ohm()[0], 0.8 * scale);
        near(circuit.resistance_ohm()[1], 1.2);
        assert_eq!(magnetic.flux_coordinates(), &[1, 0]);
        let original = circuit.state().to_vec();
        let dt = 0.01;
        let candidate = circuit.step_to_candidate(dt, circuit.conductors().to_vec(), &[1.4],
            controls(), &mut || false).unwrap();
        let l = magnetic.basis().inductance_h();
        assert!(l[1].abs() > 1e-6);
        let r = candidate.circuit().resistance_ohm();
        let expected = two_by_two([
            l[0] + dt * (1.0 + r[0]) / 2.0, l[1],
            l[2], l[3] + dt * (3.0 + r[1]) / 2.0,
        ], [l[0] - dt * (1.0 + r[0]) / 2.0 + dt * 1.4, l[2]]);
        let next_i = candidate.circuit().currents_a().unwrap();
        for i in 0..2 { near(next_i[i], expected[i]); }
        let step = candidate.step();
        let mut coil = [0.0; 2];
        let mut volumes = [0.0; 2];
        for loss in &step.cells {
            let w = loss.winding;
            let i = f64::midpoint(if w == 0 { 1.0 } else { 0.0 }, expected[w]);
            let j = magnetic.basis().profiles_per_m2()[w][loss.source.triangle] * i;
            near(loss.current_density_a_m2, j);
            near(loss.electric_field_v_m, loss.source.conductor.resistivity_ohm_m() * j);
            near(loss.heat_j, dt * loss.electric_field_v_m * j * loss.volume_m3);
            assert_eq!(loss.source.conductor.material().card_identity(), cards[w].content_hash());
            coil[w] += loss.heat_j;
            volumes[w] += loss.volume_m3;
        }
        near(volumes[0], 0.25); near(volumes[1], 0.25);
        let i0 = f64::midpoint(1.0, expected[0]);
        let i1 = expected[1] / 2.0;
        near(coil[0], dt * r[0] * i0 * i0);
        near(coil[1], dt * r[1] * i1 * i1);
        assert_eq!(step.graph_resistors.iter().map(|h| h.branch).collect::<Vec<_>>(), vec![1, 3, 5]);
        near(step.graph_resistors[0].heat_j, dt * (1.4 - i0) * (1.4 - i0));
        near(step.graph_resistors[2].heat_j, dt * 3.0 * i1 * i1);
        let q1 = 0.2 * (1.0 - dt / (2.0 * 2.0 * 0.4)) / (1.0 + dt / (2.0 * 2.0 * 0.4));
        near(candidate.circuit().state()[magnetic.dae().charge_index[0]], q1);
        near(step.graph_resistors[1].heat_j, (0.2 * 0.2 - q1 * q1) / (2.0 * 0.4));
        near(step.heat_j, step.native.dissipated);
        near(step.stored_after_j - step.stored_before_j + step.heat_j, step.native.supplied);
        assert!(step.energy_residual_j.abs() < controls().energy_tolerance_j);
        assert_eq!(circuit.state(), original);
        assert_eq!(circuit.time_s(), 0.0);
        let (accepted, _) = candidate.into_parts();
        let continuation = accepted.step_to_candidate(0.02, accepted.conductors().to_vec(), &[1.4],
            controls(), &mut || false).unwrap();
        near(continuation.step().stored_before_j, magnetic.dae().system.hamiltonian(accepted.state()));
        assert_eq!(continuation.circuit().time_s(), 0.02);
        outcomes.push((next_i[0], coil[0]));
    }
    assert!(outcomes[1].0 < outcomes[0].0);
    assert!(outcomes[1].1 > outcomes[0].1);
}

#[test]
fn electrical_basis_and_actual_extrusion_are_data_not_solver_switches() {
    let mut reference: Option<(Vec<f64>, f64)> = None;
    for basis in [ElectricalPropertyBasis::Resistivity, ElectricalPropertyBasis::Conductivity] {
        let cards = [card(2e-5, 0.0, 2e-4, basis), card(3e-5, 0.0, 2e-4, basis)];
        for depth in [0.5, 1.0] {
            let (problem, profiles, boundary) = problem(&cards, basis, depth);
            let magnetic = winding(&problem, profiles, &boundary);
            let circuit = OhmicWindingCircuit::new(&magnetic, sources(&magnetic, basis),
                0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).unwrap();
            near(circuit.resistance_ohm()[0], 1.6 * depth);
            near(circuit.resistance_ohm()[1], 2.4 * depth);
            let next = circuit.step_to_candidate(0.01, circuit.conductors().to_vec(), &[1.4],
                controls(), &mut || false).unwrap();
            for cell in &next.step().cells { assert_eq!(cell.source.conductor.source_basis(), basis); }
            if depth == 0.5 {
                let values = (next.circuit().currents_a().unwrap(), next.step().heat_j);
                if let Some((currents, heat)) = &reference {
                    for (a, b) in values.0.iter().zip(currents) { near(*a, *b); }
                    near(values.1, *heat);
                } else { reference = Some(values); }
            }
        }
    }
}

#[test]
fn live_temperature_refresh_retains_history_source_receipts_and_fixed_magnetic_storage() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [card(2e-5, 0.002, 2e-4, basis), card(3e-5, 0.002, 2e-4, basis)];
    let (problem, profiles, boundary) = problem(&cards, basis, 0.5);
    let magnetic = winding(&problem, profiles, &boundary);
    let circuit = OhmicWindingCircuit::new(&magnetic, sources(&magnetic, basis),
        0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).unwrap();
    let frozen = circuit.step_to_candidate(0.01, circuit.conductors().to_vec(), &[1.4],
        controls(), &mut || false).unwrap();
    let hot = circuit.step_to_candidate(0.01, refreshed(&circuit, &cards, [400.0, 350.0]), &[1.4],
        controls(), &mut || false).unwrap();
    near(hot.circuit().resistance_ohm()[0], 0.96);
    near(hot.circuit().resistance_ohm()[1], 1.32);
    near(hot.step().stored_before_j, frozen.step().stored_before_j);
    assert!(hot.circuit().currents_a().unwrap()[0] < frozen.circuit().currents_a().unwrap()[0]);
    for (old, new) in circuit.conductors().iter().zip(hot.circuit().conductors()) {
        let a = &old.conductor.source_property().answer().receipt;
        let b = &new.conductor.source_property().answer().receipt;
        assert_eq!(a.selected, b.selected);
        assert_eq!(a.policy, "single-claim-only");
        assert_eq!(b.policy, "pinned-claim");
    }
    let foreign = [card(4e-5, 0.002, 2e-4, basis), cards[1].clone()];
    let mut invalid = circuit.conductors().to_vec();
    invalid[0].conductor = ResolvedConductorState::try_from_material_with_basis(bundle(&foreign[0], 300.0, basis), basis).unwrap();
    assert!(matches!(circuit.step_to_candidate(0.01, invalid, &[1.4], controls(), &mut || false),
        Err(OhmicWindingError::Invalid(_))));
    assert_eq!(circuit.time_s(), 0.0);
}

#[test]
fn a_temperature_dependent_mu_cannot_hide_parameter_work_between_equal_endpoints() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [magnetic_curve_card(2e-5, 2e-4, basis), card(3e-5, 0.001, 2e-4, basis)];
    let (problem, profiles, boundary) = problem(&cards, basis, 0.5);
    let magnetic = winding(&problem, profiles, &boundary);
    let circuit = OhmicWindingCircuit::new(&magnetic, sources(&magnetic, basis),
        0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).unwrap();
    // The selected mu curve returns to exactly the original value at 2000 K.
    let candidates = refreshed(&circuit, &cards, [2000.0, 300.0]);
    assert_eq!(candidates[0].conductor.material().property("magnetic_permeability").unwrap().value_si(),
        circuit.conductors()[0].conductor.material().property("magnetic_permeability").unwrap().value_si());
    assert!(matches!(circuit.step_to_candidate(0.01, candidates, &[1.4], controls(), &mut || false),
        Err(OhmicWindingError::Invalid(_))));
}

#[test]
fn fixed_pressure_sources_drive_resistance_and_field_while_temperature_can_refresh() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [
        pressure_electrical_magnetic_card(2e-5, 2e-4, basis),
        pressure_electrical_magnetic_card(3e-5, 2e-4, basis),
    ];
    let mut outcomes = Vec::new();
    for pressure in [1e5, 2e5] {
        let scale = pressure / 1e5;
        let context = pressure_point(300.0, pressure);
        let (problem, profiles, boundary) = problem_at(&cards, basis, 0.5, &context);
        let magnetic = winding(&problem, profiles, &boundary);
        let circuit = OhmicWindingCircuit::new(&magnetic, sources(&magnetic, basis),
            0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).unwrap();
        near(circuit.resistance_ohm()[0], 0.8 * scale);
        near(circuit.resistance_ohm()[1], 1.2 * scale);
        let at = |temperature, pressure| circuit.conductors().iter().zip(circuit.conductor_windings())
            .map(|(source, &winding)| {
                let material = source.conductor.material().at_query_point(
                    &cards[winding], &pressure_point(temperature, pressure)).unwrap();
                fs_em::winding::ohmic::WindingConductor { triangle: source.triangle,
                    conductor: ResolvedConductorState::try_from_material_with_basis(material, basis).unwrap() }
            }).collect::<Vec<_>>();
        let hot = circuit.step_to_candidate(0.01, at(400.0, pressure), &[1.4],
            controls(), &mut || false).unwrap();
        near(hot.circuit().resistance_ohm()[0], 0.8 * scale);
        near(hot.circuit().resistance_ohm()[1], 1.2 * scale);
        for source in hot.circuit().conductors() {
            let mu = source.conductor.material().property("magnetic_permeability").unwrap();
            assert_eq!(mu.source_axes()[0].0, "pressure");
            assert_eq!(source.conductor.source_property().source_axes()[0].0, "pressure");
            near(mu.value_si(), 2e-4 * scale);
        }
        near(hot.step().stored_after_j - hot.step().stored_before_j + hot.step().heat_j,
            hot.step().native.supplied);
        assert!(matches!(circuit.step_to_candidate(0.01, at(400.0, 1.5e5), &[1.4],
            controls(), &mut || false), Err(OhmicWindingError::Invalid(_))));
        outcomes.push((circuit.resistance_ohm()[0], magnetic.basis().inductance_h()[0] - 0.002));
        assert_eq!(circuit.time_s(), 0.0);
    }
    near(outcomes[1].0, 2.0 * outcomes[0].0);
    near(outcomes[1].1, 2.0 * outcomes[0].1);
}

#[test]
fn incomplete_overlapping_foreign_and_cancelled_proposals_never_publish() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [card(2e-5, 0.0, 2e-4, basis), card(3e-5, 0.0, 2e-4, basis)];
    let (problem, profiles, boundary) = problem(&cards, basis, 0.5);
    let magnetic = winding(&problem, profiles.clone(), &boundary);
    let circuit = OhmicWindingCircuit::new(&magnetic, sources(&magnetic, basis),
        0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).unwrap();
    let original = circuit.state().to_vec();
    let mut incomplete = circuit.conductors().to_vec(); incomplete.pop();
    assert!(circuit.step_to_candidate(0.01, incomplete, &[1.4], controls(), &mut || false).is_err());
    let mut duplicate = circuit.conductors().to_vec(); duplicate[1] = duplicate[0].clone();
    assert!(circuit.step_to_candidate(0.01, duplicate, &[1.4], controls(), &mut || false).is_err());
    let mut bad_triangle = circuit.conductors().to_vec(); bad_triangle[0].triangle = usize::MAX;
    assert!(circuit.step_to_candidate(0.01, bad_triangle, &[1.4], controls(), &mut || false).is_err());
    let overlapping = vec![profiles[0].clone(), profiles[0].iter().zip(&profiles[1]).map(|(a, b)| a + b).collect()];
    let overlap = PlanarWindingBasis::new(&problem, overlapping, &boundary, vec![0.002, 0.003],
        winding_config(), linear(), &mut || false).unwrap().into_circuit(graph(), vec![4, 2], &mut || false).unwrap();
    assert!(OhmicWindingCircuit::new(&overlap, sources(&overlap, basis),
        0.0, &[1.0, 0.0], &[0.2], &[1.4], &mut || false).is_err());
    let mut polls = 0;
    let replay = circuit.step_to_candidate(0.01, circuit.conductors().to_vec(), &[1.4], controls(),
        &mut || { polls += 1; false }).unwrap();
    let mut count = 0;
    assert!(matches!(circuit.step_to_candidate(0.01, circuit.conductors().to_vec(), &[1.4], controls(),
        &mut || { count += 1; count == polls }), Err(OhmicWindingError::Cancelled)));
    let again = circuit.step_to_candidate(0.01, circuit.conductors().to_vec(), &[1.4], controls(), &mut || false).unwrap();
    assert_eq!(again.circuit().state(), replay.circuit().state());
    assert_eq!(again.step().heat_j, replay.step().heat_j);
    assert_eq!(circuit.state(), original);
    assert_eq!(circuit.time_s(), 0.0);
}
