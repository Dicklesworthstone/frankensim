//! G0/G3/G4: physical thermal packets for continuing non-equilibrium motion.
use super::*;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_material::tensor::from_mandel;
use fs_cutfem::elastic3::thermal::{
    ThermoelasticEquilibriumSettings3, ThermoelasticStateSettings3,
};

fn settings(physical_dofs: usize) -> ThermoelasticStateSettings3 {
    let equilibrium = ThermoelasticEquilibriumSettings3 {
        refresh: MaterialRefreshLimits3::default(),
        maximum_native_dofs: physical_dofs,
        relative_residual_tolerance: 1e-11,
        maximum_iterations: 0,
        poll_iterations: 16,
        maximum_free_residual_n: 0.0,
        maximum_kinematic_strain: TOTAL_STRAIN,
        maximum_points: POINTS,
        maximum_mass_terms: POINTS,
        maximum_prestress_entries: Some(POINTS),
        moment_origin_m: [0.0; 3],
    };
    let state = ThermoelasticStateSettings3::from(equilibrium);
    assert_eq!(state.maximum_native_dofs, physical_dofs);
    assert_eq!(state.maximum_kinematic_strain, TOTAL_STRAIN);
    assert_eq!(state.maximum_points, POINTS);
    assert_eq!(state.maximum_mass_terms, POINTS);
    assert_eq!(state.maximum_prestress_entries, Some(POINTS));
    state
}

macro_rules! physical_packet {
    ($op:ident, $physical_dofs:expr) => {{
        let controls = settings($physical_dofs);
        let scales: Vec<_> = (0..$op.cells())
            .map(|i| 0.6 + 0.25 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        let canonical = probe($op.nodes());
        let mut supplied = canonical.clone();
        for (node, fixed) in $op.fixed().iter().enumerate() {
            if *fixed {
                supplied[3 * node..3 * node + 3].copy_from_slice(&[17.0, -23.0, 31.0]);
            }
        }
        let unchanged_input = supplied.clone();
        let before_action = action(&$op, &canonical);
        let original_mass = $op.prepare_material_mass(POINTS, poll).unwrap();
        let original_entries: Vec<_> = original_mass.reference_scalar_entries().collect();
        let occupancy: Vec<_> = (0..$op.cells())
            .map(|i| 0.7 + 0.2 * i as f64 / $op.cells() as f64).collect();
        let original_mass_action = original_mass.apply(&occupancy, &canonical, poll).unwrap();

        // These inherited equilibrium controls allow zero solve iterations and
        // zero static force defect. A physical packet must not run a solve or
        // impose static balance on a moving configuration.
        let packet = $op.prepare_thermoelastic_state(&supplied, controls, poll).unwrap();
        assert_eq!(supplied, unchanged_input);
        assert_eq!(packet.displacement, canonical);
        for (node, fixed) in $op.fixed().iter().enumerate() {
            if *fixed {
                assert!(packet.displacement[3 * node..3 * node + 3]
                    .iter().all(|value| value.to_bits() == 0.0_f64.to_bits()));
            }
        }
        let full_free = $op.eigenstrain_equilibrium_load(
            &packet.eigenstrain, ReferenceLoad3::default(), None, poll,
        ).unwrap();
        let reference = $op.eigenstrain_equilibrium_response(
            &packet.eigenstrain, &canonical, &full_free, TOTAL_STRAIN,
            POINTS, [0.0; 3], poll,
        ).unwrap();
        assert!(reference.free_residual_norm > 1e-5);
        assert!(reference.residual.chunks_exact(3).zip($op.fixed())
            .any(|(force, fixed)| *fixed && force.iter().any(|f| f.abs() > 1e-8)));
        close_vectors(&packet.resistance, &reference.residual);
        $op.validate_eigenstrain(&packet.eigenstrain).unwrap();
        $op.validate_material_mass(&packet.mass).unwrap();
        assert_eq!(packet.mass.reference_scalar_entries().collect::<Vec<_>>(), original_entries);
        assert_eq!(packet.mass.apply(&occupancy, &canonical, poll).unwrap(), original_mass_action);
        $op.validate_prestress(packet.prestress.as_ref().unwrap()).unwrap();
        let direction: Vec<_> = (0..$op.n()).map(|i| (0.31 * i as f64).sin()).collect();
        let mut old_g = vec![0.0; $op.n()];
        packet.prestress.as_ref().unwrap().apply_into(&direction, &mut old_g, poll).unwrap();

        let changed: Vec<_> = canonical.iter().map(|u| 0.65 * u).collect();
        let expected_g = $op.prepare_eigenstrain_prestress(
            &packet.eigenstrain, &changed, TOTAL_STRAIN, POINTS, POINTS, poll,
        ).unwrap();
        let mut expected_g_action = vec![0.0; $op.n()];
        expected_g.apply_into(&direction, &mut expected_g_action, poll).unwrap();
        let recovered = $op.recover_thermoelastic_state(
            &changed, packet.eigenstrain, packet.mass, controls, poll,
        ).unwrap();
        assert_eq!(recovered.displacement, changed);
        let full_free = $op.eigenstrain_equilibrium_load(
            &recovered.eigenstrain, ReferenceLoad3::default(), None, poll,
        ).unwrap();
        let reference = $op.eigenstrain_equilibrium_response(
            &recovered.eigenstrain, &changed, &full_free,
            TOTAL_STRAIN, POINTS, [0.0; 3], poll,
        ).unwrap();
        close_vectors(&recovered.resistance, &reference.residual);
        assert_eq!(recovered.mass.reference_scalar_entries().collect::<Vec<_>>(), original_entries);
        assert_eq!(recovered.mass.apply(&occupancy, &canonical, poll).unwrap(), original_mass_action);
        $op.validate_eigenstrain(&recovered.eigenstrain).unwrap();
        $op.validate_material_mass(&recovered.mass).unwrap();
        $op.validate_prestress(recovered.prestress.as_ref().unwrap()).unwrap();
        let mut new_g = vec![0.0; $op.n()];
        recovered.prestress.as_ref().unwrap().apply_into(&direction, &mut new_g, poll).unwrap();
        close_vectors(&new_g, &expected_g_action);
        assert!(new_g.iter().zip(old_g).any(|(a, b)| (a - b).abs() > 1e-8));

        // Optional G is only storage policy; complete physical force recovery
        // remains available with the SAME moved source and mass preparations.
        let without_g = $op.recover_thermoelastic_state(
            &changed, recovered.eigenstrain, recovered.mass,
            ThermoelasticStateSettings3 { maximum_prestress_entries: None, ..controls }, poll,
        ).unwrap();
        assert!(without_g.prestress.is_none());
        close_vectors(&without_g.resistance, &reference.residual);
        assert_eq!(without_g.mass.reference_scalar_entries().collect::<Vec<_>>(), original_entries);
        assert_eq!(action(&$op, &canonical), before_action);
        assert_eq!($op.scales(), scales.as_slice());
    }};
}

#[test]
fn raw_physical_packets_retain_full_thermal_resistance_without_solving_static_equilibrium() {
    let law = law();
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let mut op = cartesian(ElasticMaterialField3::ThermoelasticByCell(&hot), formulation);
        physical_packet!(op, op.n());
    }
}

#[test]
fn adaptive_physical_packets_preserve_hanging_coordinates_consistent_mass_and_signed_g() {
    let law = law();
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let tree = tree();
    let mut op = adaptive(
        ElasticMaterialField3::ThermoelasticByCell(&hot), &tree, ElasticityFormulation3::BBar,
    );
    assert!(op.physical_nodes().len() > op.nodes().len());
    let physical_nodes = op.physical_nodes().to_vec();
    let field = probe(op.nodes());
    let reconstructed = op.physical_displacements(&field).unwrap();
    let physical_dofs = 3 * op.physical_nodes().len();
    physical_packet!(op, physical_dofs);
    assert_eq!(op.physical_nodes(), physical_nodes.as_slice());
    assert_eq!(op.physical_displacements(&field).unwrap(), reconstructed);
    assert!(op.prepare_thermoelastic_state(
        &field, ThermoelasticStateSettings3 {
            maximum_native_dofs: physical_dofs - 1, ..settings(physical_dofs)
        }, poll,
    ).is_err());
}

macro_rules! refuse_wrong_owner {
    ($op:ident, $foreign:ident, $changed:ident, $physical_dofs:expr) => {{
        let controls = settings($physical_dofs);
        let u = probe($op.nodes());
        let before = action(&$op, &u);
        assert!( $op.recover_thermoelastic_state(
            &u, $foreign.prepare_thermoelastic(POINTS, poll).unwrap(),
            $op.prepare_material_mass(POINTS, poll).unwrap(), controls, poll,
        ).is_err());
        assert!( $op.recover_thermoelastic_state(
            &u, $op.prepare_thermoelastic(POINTS, poll).unwrap(),
            $foreign.prepare_material_mass(POINTS, poll).unwrap(), controls, poll,
        ).is_err());
        assert_eq!(action(&$op, &u), before);
        // A valid token alone cannot turn a caller-authored free strain into
        // this source's integrated thermal expansion.
        let other_free = $op.prepare_eigenstrain(&|_, _| Ok([0.0; 6]), POINTS, poll).unwrap();
        $op.validate_eigenstrain(&other_free).unwrap();
        assert!($op.recover_thermoelastic_state(
            &u, other_free, $op.prepare_material_mass(POINTS, poll).unwrap(), controls, poll,
        ).is_err());
        let changed_surface = Cell::new(false);
        let surface_only = $op.prepare_eigenstrain(&|cell, position| {
            let source = $op.cell_thermoelastic_state(cell).unwrap().source();
            let mut free = from_mandel(&source.free_strain_mandel());
            if (position[2] - 0.73).abs() < 1e-12 {
                free[0] += 0.00001;
                changed_surface.set(true);
            }
            Ok(free)
        }, POINTS, poll).unwrap();
        assert!(changed_surface.get());
        assert!($op.recover_thermoelastic_state(
            &u, surface_only, $op.prepare_material_mass(POINTS, poll).unwrap(), controls, poll,
        ).is_err());
        let stale_source = $op.prepare_thermoelastic(POINTS, poll).unwrap();
        let stale_mass = $op.prepare_material_mass(POINTS, poll).unwrap();
        $op.refresh_material_field(
            ElasticMaterialField3::ThermoelasticByCell(&$changed),
            MaterialRefreshLimits3::default(), poll,
        ).unwrap();
        let protected = $op.prepare_thermoelastic(POINTS, poll).unwrap();
        let before = action(&$op, &u);
        assert!( $op.recover_thermoelastic_state(
            &u, stale_source, $op.prepare_material_mass(POINTS, poll).unwrap(), controls, poll,
        ).is_err());
        assert!( $op.recover_thermoelastic_state(
            &u, $op.prepare_thermoelastic(POINTS, poll).unwrap(), stale_mass, controls, poll,
        ).is_err());
        $op.validate_eigenstrain(&protected).unwrap();
        assert_eq!(action(&$op, &u), before);
    }};
}

#[test]
fn physical_recovery_rejects_each_foreign_and_stale_preparation_in_both_spaces() {
    let law = law();
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let changed = |bounds| Ok(frozen(&law, temperature(bounds) + 1.0));
    let formulation = ElasticityFormulation3::BBar;
    let mut op = cartesian(ElasticMaterialField3::ThermoelasticByCell(&hot), formulation);
    let foreign = cartesian(ElasticMaterialField3::ThermoelasticByCell(&hot), formulation);
    refuse_wrong_owner!(op, foreign, changed, op.n());
    let tree = tree();
    let mut op = adaptive(ElasticMaterialField3::ThermoelasticByCell(&hot), &tree, formulation);
    let foreign = adaptive(ElasticMaterialField3::ThermoelasticByCell(&hot), &tree, formulation);
    refuse_wrong_owner!(op, foreign, changed, 3 * op.physical_nodes().len());
}

#[test]
fn state_packet_bounds_physical_refusal_and_late_cancellation_leave_native_owner_unchanged() {
    let law = law();
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let op = cartesian(
        ElasticMaterialField3::ThermoelasticByCell(&hot), ElasticityFormulation3::BBar,
    );
    let u = probe(op.nodes());
    let before = action(&op, &u);
    let protected = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let controls = settings(op.n());
    for insufficient in [
        ThermoelasticStateSettings3 { maximum_native_dofs: op.n() - 1, ..controls },
        ThermoelasticStateSettings3 { maximum_points: 0, ..controls },
        ThermoelasticStateSettings3 { maximum_mass_terms: 0, ..controls },
        ThermoelasticStateSettings3 { maximum_prestress_entries: Some(0), ..controls },
    ] {
        assert!(op.prepare_thermoelastic_state(&u, insufficient, poll).is_err());
        assert!(op.recover_thermoelastic_state(
            &u, op.prepare_thermoelastic(POINTS, poll).unwrap(),
            op.prepare_material_mass(POINTS, poll).unwrap(), insufficient, poll,
        ).is_err());
    }
    let outside: Vec<_> = u.iter().map(|u| 100.0 * u).collect();
    assert!(op.prepare_thermoelastic_state(
        &outside, ThermoelasticStateSettings3 { maximum_prestress_entries: None, ..controls }, poll,
    ).is_err());
    assert!(op.prepare_thermoelastic_state(&u[..u.len() - 1], controls, poll).is_err());
    let mut calls = 0;
    op.prepare_thermoelastic_state(&u, controls, || { calls += 1; poll() }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut at = 0;
        assert!(matches!(op.prepare_thermoelastic_state(&u, controls, || {
            at += 1;
            if at == stop { ControlFlow::Break(()) } else { poll() }
        }), Err(ElasticityError3::Cancelled)));
    }
    let mut calls = 0;
    op.recover_thermoelastic_state(
        &u, op.prepare_thermoelastic(POINTS, poll).unwrap(),
        op.prepare_material_mass(POINTS, poll).unwrap(), controls,
        || { calls += 1; poll() },
    ).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut at = 0;
        assert!(matches!(op.recover_thermoelastic_state(
            &u, op.prepare_thermoelastic(POINTS, poll).unwrap(),
            op.prepare_material_mass(POINTS, poll).unwrap(), controls,
            || { at += 1; if at == stop { ControlFlow::Break(()) } else { poll() } },
        ), Err(ElasticityError3::Cancelled)));
    }
    op.validate_eigenstrain(&protected).unwrap();
    assert_eq!(action(&op, &u), before);
}
