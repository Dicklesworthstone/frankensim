//! Global voltage compatibility, circulating current and finite inventories.
use super::*;

fn ionic() -> ZeroCurrentDiffusion {
    ZeroCurrentDiffusion::new(vec![1, -1], 1e-15, 1e-16).unwrap()
}

fn selective_law(model: &Arc<VirialSolutionModel>, plus: f64, minus: f64) -> SolutionMaxwellStefan {
    SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        vec![
            (SpeciesId::new("solute-0").unwrap(), MolarMass::new(0.04)),
            (SpeciesId::new("solute-1").unwrap(), MolarMass::new(0.06)),
        ],
        vec![0.0, plus, minus, plus, 0.0, 5e-10, minus, 5e-10, 0.0],
    )
    .unwrap()
}

pub(super) fn loop_network(controls: ElectroneutralDiffusionControls) -> OsmoticNetwork {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let links = [
        (0, 1, 1e-4, 1e-9, 3e-9),
        (0, 2, 2e-4, 3e-9, 1e-9),
        (1, 2, 1.5e-4, 2e-9, 1e-9),
    ]
    .into_iter()
    .map(|(a, b, area, plus, minus)| {
        SolutionDiffusionLink::new(a, b, area, 0.001, selective_law(&model, plus, minus)).unwrap()
    })
    .collect();
    OsmoticNetwork::with_solution_diffusion(
        vec!["solute-0".into(), "solute-1".into()],
        [100.0, 40.0, 15.0].map(|c| cell(&model, &[c, c])).into(),
        links,
        physical(),
        1 << 22,
    )
    .unwrap()
    .with_electroneutral_diffusion(ionic(), controls, 1 << 22)
    .unwrap()
}

#[test]
fn cyclic_diffusion_has_one_compatible_potential_and_nonzero_circulating_current() {
    let mut system = loop_network(ElectroneutralDiffusionControls::default());
    let mut response = Vec::new();
    for link in &mut system.diffusion {
        let mut flux = [0.0; 2];
        let DiffusionMaterial::Constant(law) = &link.material else {
            unreachable!()
        };
        let sample = law
            .current_response_into(
                &system.cells[link.left].solution,
                &system.cells[link.right].solution,
                link.distance_m,
                physical().constitutive,
                link.electrical.as_deref().unwrap(),
                &mut link.workspace,
                &mut flux,
                || false,
            )
            .unwrap();
        response.push((
            sample.chemical_charge_flux_mol_m2_s(),
            sample.electric_conductance_mol_m2_s(),
            link.area_m2,
        ));
    }
    // Independent circuit identity: sum of chemical diffusion voltages around
    // the loop divided by sum of physical edge resistances. Face order is
    // 01,02,12, so the 02 contribution has the opposite loop orientation.
    let voltage = response[0].0 / response[0].1 - response[1].0 / response[1].1
        + response[2].0 / response[2].1;
    let resistance: f64 = response.iter().map(|&(_, g, area)| 1.0 / (area * g)).sum();
    let circulating = voltage / resistance;
    assert!(
        voltage.abs() > 0.1,
        "independent edge zero-current voltages are incompatible"
    );
    assert!(circulating.abs() > 1e-10);
    let mut electrical = system.electroneutral_diffusion.take().unwrap();
    let mut evaluations = 0;
    electrical
        .sweep(
            &mut system.cells,
            &mut system.diffusion,
            0.01,
            physical(),
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    assert_eq!(electrical.potentials[0], 0.0);
    for (edge, &(chemical, conductance, area)) in system.diffusion.iter().zip(&response) {
        let difference = electrical.potentials[edge.right] - electrical.potentials[edge.left];
        let actual = area * (chemical - conductance * difference);
        let expected = if edge.endpoints() == [0, 2] {
            -circulating
        } else {
            circulating
        };
        assert!(
            (actual - expected).abs() < 2e-17,
            "{actual:e} vs {expected:e}"
        );
    }
    assert_eq!(evaluations, 6);
    assert_eq!(electrical.report().face_evaluations, evaluations);
    assert!(electrical.report().maximum_node_charge_rate_mol_s <= 1e-16);
    assert!((electrical.report().maximum_edge_charge_rate_mol_s - circulating.abs()).abs() < 2e-17);
    for cell in &system.cells {
        ionic()
            .admit_inventories(cell.solution.solute_moles())
            .unwrap();
    }
}

#[test]
fn cyclic_ionic_evolution_conserves_each_species_and_global_free_energy() {
    let mut system = loop_network(ElectroneutralDiffusionControls::default());
    assert!(system.zero_current_diffusion().is_none());
    let before = bits(&system);
    let report = system.advance(2.0, numerical(1e-6), || false).unwrap();
    assert_eq!(system.time_s(), 2.0);
    assert_ne!(bits(&system), before);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    assert!(
        report.solute_balance_residual_mol[..2]
            .iter()
            .all(|n| n.abs() <= 1e-16)
    );
    let electrical = report.electroneutral_diffusion.unwrap();
    assert!(electrical.maximum_edge_charge_rate_mol_s > 1e-10);
    assert!(electrical.maximum_node_charge_rate_mol_s <= 1e-16);
    assert_eq!(electrical.face_evaluations, report.residual_evaluations);
    assert!(electrical.potential_solves >= 3 * report.accepted_steps);
    for cell in system.compartments() {
        ionic()
            .admit_inventories(cell.solution().solute_moles())
            .unwrap();
    }
}

#[test]
fn isolated_and_carrier_free_components_have_independent_voltage_gauges() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let mut system = network(
        &model,
        &[
            &[100.0, 100.0],
            &[40.0, 40.0],
            &[0.0, 0.0],
            &[0.0, 0.0],
            &[60.0, 60.0],
        ],
        &[(0, 1, 1e-4), (2, 3, 1e-4)],
    )
    .with_electroneutral_diffusion(ionic(), ElectroneutralDiffusionControls::default(), 1 << 22)
    .unwrap();
    let absent = system.cells[2..]
        .iter()
        .map(|c| c.solution.solute_moles().to_vec())
        .collect::<Vec<_>>();
    let report = system.advance(1.0, numerical(1e-6), || false).unwrap();
    assert!(
        report
            .electroneutral_diffusion
            .unwrap()
            .maximum_edge_charge_rate_mol_s
            <= 1e-16
    );
    for (cell, saved) in system.cells[2..].iter().zip(absent) {
        assert_eq!(cell.solution.solute_moles(), saved);
    }
    for &node in &[0usize, 2, 3, 4] {
        assert_eq!(
            system.electroneutral_diffusion.as_ref().unwrap().potentials[node],
            0.0
        );
    }
}

#[test]
fn global_work_exhaustion_and_late_cancellation_publish_no_partial_state() {
    let mut controls = ElectroneutralDiffusionControls::default();
    // One complete coarse solve/replay, followed by exhaustion during the
    // fine replay. Work is not reset between the three adaptive probes.
    controls.maximum_face_evaluations = 11;
    let mut bounded = loop_network(controls);
    let saved = bits(&bounded);
    assert!(matches!(
        bounded.advance(0.1, numerical(1e-4), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&bounded), saved);
    let mut reference = loop_network(ElectroneutralDiffusionControls::default());
    let mut polls = 0;
    reference
        .advance(0.1, numerical(1e-4), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 10);
    let mut cancelled = loop_network(ElectroneutralDiffusionControls::default());
    let saved = bits(&cancelled);
    let mut seen = 0;
    assert!(matches!(
        cancelled.advance(0.1, numerical(1e-4), || {
            seen += 1;
            seen == polls - 1
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&cancelled), saved);
}

#[test]
fn global_mode_refuses_incompatible_implicit_coordinates_and_checks_storage() {
    let system = loop_network(ElectroneutralDiffusionControls::default());
    assert!(
        system
            .clone()
            .with_implicit_salt_diffusion(ImplicitSaltDiffusionControls::default(), 1 << 22)
            .is_err()
    );
    assert!(
        system
            .clone()
            .with_implicit_multicomponent_diffusion(
                ImplicitMulticomponentDiffusionControls::default(),
                1 << 22
            )
            .is_err()
    );
    assert!(
        system
            .clone()
            .with_zero_current_diffusion(ionic(), 1 << 22)
            .is_err()
    );
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let original = network(&model, &[&[100.0, 100.0], &[40.0, 40.0]], &[(0, 1, 1e-4)]);
    let budget = original.workspace_bytes();
    assert!(matches!(
        original.with_electroneutral_diffusion(
            ionic(),
            ElectroneutralDiffusionControls::default(),
            budget
        ),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
}

#[test]
fn cyclic_source_transport_uses_actual_neutral_face_coefficients_after_evolution() {
    use crate::osmosis::coupled::CoupledSolutionControls;
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let sourced = super::source::multicomponent_material(&model);
    let system = OsmoticNetwork::with_electroneutral_solution_transport(
        vec!["solute-0".into(), "solute-1".into()],
        [100.0, 40.0, 15.0].map(|c| cell(&model, &[c, c])).into(),
        Vec::new(),
        vec![
            SolutionDiffusionLink::new(0, 1, 0.001, 0.001, selective_law(&model, 1e-9, 3e-9))
                .unwrap(),
            SolutionDiffusionLink::from_material_sources(0, 2, 0.001, 0.001, sourced).unwrap(),
            SolutionDiffusionLink::new(1, 2, 0.001, 0.001, selective_law(&model, 2e-9, 1e-9))
                .unwrap(),
        ],
        CoupledSolutionControls {
            diffusion: physical(),
            water_balance_tolerance_m3: 1e-18,
        },
        ionic(),
        ElectroneutralDiffusionControls::default(),
        1 << 22,
    )
    .unwrap();
    for mut system in [
        system.clone(),
        system
            .with_implicit_electroneutral_diffusion(
                ImplicitElectroneutralDiffusionControls::default(),
                1 << 22,
            )
            .unwrap(),
    ] {
        let report = system.advance(2.0, numerical(1e-5), || false).unwrap();
        assert!(
            report
                .electroneutral_diffusion
                .unwrap()
                .maximum_edge_charge_rate_mol_s
                > 0.0
        );
        assert!(
            report
                .coupled_potential_change
                .unwrap()
                .total_potential_change_j
                < 0.0
        );
        let link = &system.solution_diffusion_links()[1];
        let (left, right) = (
            system.compartments()[0].solution(),
            system.compartments()[2].solution(),
        );
        let (a, b) = (
            left.solute_moles()[0] / VOLUME,
            right.solute_moles()[0] / VOLUME,
        );
        let coefficient = 1e-9 * (1.0 + (a - b) / (a / b).ln() / 100.0);
        let original = link.material_sources().unwrap();
        let resolved = original
            .resolve_zero_current_inventories(
                [VOLUME, VOLUME],
                [left.solute_moles(), right.solute_moles()],
                &ionic(),
                || false,
            )
            .unwrap();
        assert!(
            (coefficient - 1.6e-9).abs() > 1e-11,
            "reference-state coefficients cannot represent the actual face"
        );
        for pair in resolved.pairs() {
            assert!((pair.face.properties()[0].value_si() - coefficient).abs() < 1e-21);
            assert!(
                original
                    .sources()
                    .iter()
                    .any(|source| source.selected_claim()
                        == pair.face.properties()[0].answer().receipt.selected)
            );
        }
        for cell in system.compartments() {
            ionic()
                .admit_inventories(cell.solution().solute_moles())
                .unwrap();
        }
    }
}
