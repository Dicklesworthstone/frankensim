//! G1/G3: native elastic history, objective tensile onset, and fracture replay.
//! These synthetic coupons check implementation, not a calibrated crack path.

use super::*;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_solid::brittle::solid::initiation::{
    TensileElasticTetSolid, TensileTetOnset, TensileTetProbe,
};

fn initiation_settings() -> TensileInitiationSettings {
    TensileInitiationSettings {
        maximum_strain: 0.03,
        relative_stress_tolerance: 1e-8,
        fraction_tolerance: 1e-10,
        maximum_bisections: 64,
    }
}

/// Transverse normal coordinates remain free under the native elastic solve.
fn prescribed_y(strain: f64) -> Vec<(usize, f64)> {
    (0..12)
        .filter(|&dof| dof != 3 && dof != 11)
        .map(|dof| (dof, if dof == 7 { 0.02 * strain } else { 0.0 }))
        .collect()
}

fn elastic<'a>(probe: TensileTetProbe<'a>) -> TensileElasticTetSolid<'a> {
    match probe {
        TensileTetProbe::Elastic(solid) => solid,
        TensileTetProbe::Onset(_) => panic!("expected an entirely elastic endpoint"),
    }
}

fn onset<'a>(probe: TensileTetProbe<'a>) -> TensileTetOnset<'a> {
    match probe {
        TensileTetProbe::Onset(onset) => onset,
        TensileTetProbe::Elastic(_) => panic!("expected a tensile crossing"),
    }
}

#[test]
fn discovered_plane_preserves_elastic_preload_and_drives_native_tensile_softening() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "discovered y plane");
        let materials = [material(&source, &geometry, 0)];
        assert_eq!(materials[0].law().normal(), [1.0, 0.0, 0.0]);
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let limit_strain = STRENGTH / YOUNG;
        let preload = 0.4 * limit_strain;
        let target = 0.002;
        let accepted = elastic(
            virgin
                .probe(
                    &prescribed_y(preload),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        close(
            accepted.solution().element_stress[0][1],
            YOUNG * preload,
            2e-10,
        );
        close(
            accepted.stored_energy_j(),
            0.5 * YOUNG * preload * preload * geometry.element_volumes_m3()[0],
            2e-10,
        );
        let snapshot = accepted.solution().clone();
        let prediction = onset(
            accepted
                .probe(
                    &prescribed_y(target),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        assert_eq!(prediction.initial_solution(), &snapshot);
        assert_eq!(prediction.first_candidate_element(), 0);
        let fraction = (limit_strain - preload) / (target - preload);
        let bracket = prediction.first_fraction_bracket();
        assert!(bracket[0] <= fraction + 1e-9 && bracket[1] >= fraction - 1e-9);
        assert!(bracket[1] - bracket[0] <= initiation_settings().fraction_tolerance);
        let event = prediction.events()[0].as_ref().unwrap();
        let normal = event.native().normal();
        assert!(normal[0].abs() < 1e-12 && normal[2].abs() < 1e-12);
        close(normal[1].abs(), 1.0, 2e-12);
        assert_eq!(event.source().material(), materials[0].material());
        close(
            prediction.elastic_predictor().element_stress[0][1],
            YOUNG * target,
            2e-10,
        );

        let bands = prediction.resolve_bands(cx).unwrap();
        let directional =
            TetCrackBandGeometry::from_element(&geometry, 0, [0.0, 1.0, 0.0]).unwrap();
        close(
            bands[0].law().band_width_m(),
            directional.band_width_m(),
            2e-12,
        );
        let handoff = prediction.handoff(&bands, 1e-9, cx).unwrap();
        assert_eq!(
            handoff.solid.solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            handoff.solid.solution().element_strain,
            snapshot.element_strain
        );
        close(handoff.elastic_storage_j, accepted.stored_energy_j(), 2e-10);
        close(
            handoff.solid.stored_energy_j(),
            accepted.stored_energy_j(),
            2e-10,
        );
        assert!(handoff.storage_defect_j.abs() < 1e-12);
        assert!(handoff.maximum_stress_defect_pa < 1e-6 * STRENGTH);
        assert_eq!(handoff.solid.fracture_dissipation_j(), 0.0);
        assert_eq!(handoff.solid.solution().element_states[0].damage(), 0.0);
        assert!(handoff.solid.solution().increments.is_empty());
        assert!(std::ptr::eq(
            handoff.solid.materials()[0].law(),
            bands[0].law()
        ));

        let fractured = handoff
            .solid
            .advance(
                &prescribed_y(target),
                &[],
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        let failure_strain = 2.0 * FRACTURE_WORK / (STRENGTH * directional.band_width_m());
        let expected = STRENGTH * (failure_strain - target) / (failure_strain - limit_strain);
        close(
            fractured.solid.solution().element_stress[0][1],
            expected,
            2e-8,
        );
        assert!(fractured.solid.solution().element_states[0].damage() > 0.0);
        assert!(fractured.energy.fracture_dissipation_j > 0.0);
        close(
            fractured.energy.stored_before_j,
            accepted.stored_energy_j(),
            2e-10,
        );
        close(
            fractured.energy.material_work_j,
            fractured.solid.stored_energy_j() - accepted.stored_energy_j()
                + fractured.energy.fracture_dissipation_j,
            2e-9,
        );
        assert_eq!(accepted.solution(), &snapshot);
    });
}

#[test]
fn source_plane_count_and_cancellation_refusals_preserve_elastic_history_for_retry() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "onset rollback");
        let materials = [material(&source, &geometry, 0)];
        assert!(
            TensileElasticTetSolid::try_new(&geometry, &[], TetAssemblyBudget::standard(), cx)
                .is_err()
        );
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let accepted = elastic(
            virgin
                .probe(
                    &prescribed_y(0.25 * STRENGTH / YOUNG),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let snapshot = accepted.solution().clone();
        let storage = accepted.stored_energy_j();
        let prediction = onset(
            accepted
                .probe(
                    &prescribed_y(0.002),
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        let bands = prediction.resolve_bands(cx).unwrap();
        let foreign_card = card(STRENGTH, "same numbers, different source");
        let foreign = [resolve_brittle_crack_band(
            &foreign_card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            bands[0].law().settings(),
        )
        .unwrap()];
        assert_eq!(foreign[0].law().parameters(), bands[0].law().parameters());
        assert!(prediction.handoff(&foreign, 1e-9, cx).is_err());
        assert!(prediction.handoff(&materials, 1e-9, cx).is_err());
        assert!(prediction.handoff(&[], 1e-9, cx).is_err());
        with_cx(true, |cancelled| {
            assert!(prediction.resolve_bands(cancelled).is_err());
            assert!(prediction.handoff(&bands, 1e-9, cancelled).is_err());
            assert!(
                accepted
                    .probe(
                        &prescribed_y(0.002),
                        &[],
                        IncrementSettings::default(),
                        initiation_settings(),
                        None,
                        cancelled,
                    )
                    .is_err()
            );
        });
        assert_eq!(accepted.solution(), &snapshot);
        assert_eq!(accepted.stored_energy_j(), storage);
        let first = prediction.handoff(&bands, 1e-9, cx).unwrap();
        let repeated = prediction.handoff(&bands, 1e-9, cx).unwrap();
        assert_eq!(first.solid.solution(), repeated.solid.solution());
        assert_eq!(first.storage_defect_j, repeated.storage_defect_j);
        assert_eq!(first.solid.fracture_dissipation_j(), 0.0);
    });
}

#[test]
fn repeated_principal_tension_requires_a_physical_hint_instead_of_the_source_plane() {
    with_cx(false, |cx| {
        let nodes = points(0.03);
        let geometry =
            TetLoadGeometry3::new(&nodes, &TETS, TetAssemblyBudget::standard(), 4, cx).unwrap();
        let source = card(STRENGTH, "hydrostatic onset");
        let materials = [material(&source, &geometry, 0)];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let strain = 2.0 * (1.0 - 2.0 * POISSON) * STRENGTH / YOUNG;
        // One shear coordinate remains free; its zero solution is determined
        // by the native equilibrium equations.
        let dilation: Vec<_> = (0..12)
            .filter(|&dof| dof != 4)
            .map(|dof| {
                let length = match dof {
                    3 => 0.03,
                    7 => 0.02,
                    11 => 0.01,
                    _ => 0.0,
                };
                (dof, length * strain)
            })
            .collect();
        assert!(
            virgin
                .probe(
                    &dilation,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .is_err()
        );
        let hints = [Some([0.0, 1.0, 0.0])];
        let prediction = onset(
            virgin
                .probe(
                    &dilation,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    Some(&hints),
                    cx,
                )
                .unwrap(),
        );
        let event = prediction.events()[0].as_ref().unwrap();
        assert_eq!(event.native().principal_multiplicity(), 3);
        let normal = event.native().normal();
        assert!(normal[0].abs() < 1e-12 && normal[2].abs() < 1e-12);
        close(normal[1].abs(), 1.0, 2e-12);
        close(prediction.first_fraction_bracket()[1], 0.5, 2e-8);
        let bands = prediction.resolve_bands(cx).unwrap();
        let handoff = prediction.handoff(&bands, 1e-9, cx).unwrap();
        assert_eq!(handoff.solid.stored_energy_j(), 0.0);
        assert_eq!(handoff.solid.fracture_dissipation_j(), 0.0);
        assert_eq!(virgin.stored_energy_j(), 0.0);
    });
}

#[test]
fn independent_sources_keep_nontriggered_elements_elastic_through_fracture_replay() {
    with_cx(false, |cx| {
        let mut nodes = points(0.03).to_vec();
        nodes.extend(points(0.03).map(|mut x| {
            x[0] += 0.1;
            x
        }));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let geometry =
            TetLoadGeometry3::new(&nodes, &tets, TetAssemblyBudget::standard(), 8, cx).unwrap();
        let weak = card(STRENGTH, "weak onset");
        let strong = card(2.0 * STRENGTH, "strong onset");
        let materials = [
            material(&weak, &geometry, 0),
            material(&strong, &geometry, 1),
        ];
        let virgin = TensileElasticTetSolid::try_new(
            &geometry,
            &materials,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let targets: Vec<_> = prescribed_y(1.5 * STRENGTH / YOUNG)
            .into_iter()
            .flat_map(|(dof, value)| [(dof, value), (dof + 12, value)])
            .collect();
        let prediction = onset(
            virgin
                .probe(
                    &targets,
                    &[],
                    IncrementSettings::default(),
                    initiation_settings(),
                    None,
                    cx,
                )
                .unwrap(),
        );
        assert!(prediction.events()[0].is_some());
        assert!(prediction.events()[1].is_none());
        assert_eq!(prediction.first_candidate_element(), 0);
        let bands = prediction.resolve_bands(cx).unwrap();
        assert_eq!(bands[0].material(), materials[0].material());
        assert_eq!(bands[1].material(), materials[1].material());
        assert_eq!(bands[1].law().normal(), materials[1].law().normal());
        let handoff = prediction.handoff(&bands, 1e-9, cx).unwrap();
        let step = handoff
            .solid
            .advance(&targets, &[], IncrementSettings::default(), controls(), cx)
            .unwrap();
        assert!(step.solid.solution().element_states[0].damage() > 0.0);
        assert_eq!(step.solid.solution().element_states[1].damage(), 0.0);
        assert!(step.solid.solution().element_stress[0][1] < STRENGTH);
        close(
            step.solid.solution().element_stress[1][1],
            1.5 * STRENGTH,
            2e-8,
        );
        assert_eq!(virgin.stored_energy_j(), 0.0);
    });
}
