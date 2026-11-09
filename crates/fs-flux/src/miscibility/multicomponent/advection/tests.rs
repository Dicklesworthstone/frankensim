//! Independent parcel chronology, counterflow and fixed-volume conservation.

use super::*;
use crate::miscibility::multicomponent::MulticomponentMiscibilityError;
use crate::miscibility::multicomponent::tests::{
    face, field, gibbs_energy, matrix, near, network, uniform_matrix,
};

fn opposing(network: &MulticomponentMiscibilityNetwork, rate: f64) -> MulticomponentMiscibilityAdvection {
    network.prepare_advection(vec![MulticomponentMiscibilityAdvectionFace {
        forward_m3_s: rate, reverse_m3_s: rate,
    }], 0.0, || false).unwrap()
}

fn circulation(network: &MulticomponentMiscibilityNetwork, rate: f64) -> MulticomponentMiscibilityAdvection {
    network.prepare_advection(vec![
        MulticomponentMiscibilityAdvectionFace { forward_m3_s: rate, reverse_m3_s: 0.0 },
        MulticomponentMiscibilityAdvectionFace { forward_m3_s: 0.0, reverse_m3_s: rate },
        MulticomponentMiscibilityAdvectionFace { forward_m3_s: rate, reverse_m3_s: 0.0 },
    ], 0.0, || false).unwrap()
}

#[test]
fn a_three_species_circulation_uses_only_the_original_synchronous_donors() {
    let chi = matrix(3, &[(0, 1, 0.2), (0, 2, -0.3), (1, 2, 0.1)]);
    // Canonical faces become [0,1], [0,2], [1,2], retaining the cycle
    // 0->1->2->0 rather than cascading a newly received parcel onward.
    let network = network(&[0.002, 0.028, 0.044], chi.clone(), &[0.001; 3],
        vec![face(2, 1), face(2, 0), face(1, 0)], vec![uniform_matrix(3, 1e-4)]);
    let fractions = [[0.8, 0.15, 0.05], [0.2, 0.5, 0.3], [0.1, 0.2, 0.7]];
    let initial = field(&network, &[&fractions[0], &fractions[1], &fractions[2]]);
    let flow = circulation(&network, 0.0005);
    let candidate = initial.advect_candidate(&flow, 0.5, MulticomponentMiscibilityAdvectionControls::default(), || false).unwrap();
    assert_eq!(candidate.step().substeps, 1);
    near(candidate.step().maximum_courant, 0.25, 2e-15);
    for (cell, donor) in [2, 0, 1].into_iter().enumerate() {
        for species in 0..3 {
            let expected = 0.75 * fractions[cell][species] + 0.25 * fractions[donor][species];
            near(candidate.field().cell_component_moles(cell).unwrap()[species], expected, 2e-15);
        }
    }
    for (face, (donor, sign)) in [(0, 1.0), (2, -1.0), (1, 1.0)].into_iter().enumerate() {
        for species in 0..3 {
            near(candidate.step().transferred_mol[face][species], sign * 0.25 * fractions[donor][species], 2e-15);
        }
    }
    let change = gibbs_energy(candidate.field(), 3, &chi) - gibbs_energy(&initial, 3, &chi);
    assert!(change < 0.0);
    near(candidate.step().free_energy_change_j, change, 4e-9);
    assert_eq!(candidate.field().elapsed_s(), 0.5);
    assert_eq!(initial.elapsed_s(), 0.0);
}

#[test]
fn zero_net_volume_counterflow_mixes_every_species_through_bounded_substeps() {
    let network = network(&[0.002, 0.028, 0.044], vec![0.0; 9], &[0.001, 0.002],
        vec![face(0, 1)], vec![uniform_matrix(3, 1e-4)]);
    let left = [0.2, 0.3, 0.5];
    let right = [0.6, 0.15, 0.25];
    let initial = field(&network, &[&left, &right]);
    let flow = opposing(&network, 0.0005);
    let candidate = initial.advect_candidate(&flow, 2.0,
        MulticomponentMiscibilityAdvectionControls { maximum_courant: 0.25, ..Default::default() },
        || false,
    ).unwrap();
    assert_eq!(candidate.step().substeps, 4);
    near(candidate.step().maximum_courant, 0.25, 2e-15);
    // Each 0.5-s substep multiplies the concentration contrast by
    // 1 - Q*dt*(1/V_left+1/V_right) = 0.625. The volume-weighted mean stays fixed.
    for species in 0..3 {
        let mean = (left[species] + 2.0 * right[species]) / 3.0;
        let contrast = (right[species] - left[species]) * 0.625_f64.powi(4);
        let expected = [mean - 2.0 * contrast / 3.0, mean + contrast / 3.0];
        for cell in 0..2 {
            near(candidate.field().cell_mole_fractions(cell).unwrap()[species], expected[cell], 3e-15);
            let delta = candidate.field().cell_component_moles(cell).unwrap()[species]
                - initial.cell_component_moles(cell).unwrap()[species];
            let sign = if cell == 0 { -1.0 } else { 1.0 };
            near(delta, sign * candidate.step().transferred_mol[0][species], 2e-14);
        }
        near(candidate.field().component_moles()[species] + candidate.field().component_moles()[3 + species],
            left[species] + 2.0 * right[species], 2e-14);
    }
    for cell in 0..2 {
        near(candidate.field().cell_component_moles(cell).unwrap().iter().sum(), [1.0, 2.0][cell], 2e-14);
    }
    assert!(candidate.step().transferred_mol[0].iter().any(|q| q.abs() > 0.1));
    assert!(candidate.step().free_energy_change_j < 0.0);
}

#[test]
fn uniform_composition_survives_unequal_volumes_and_masses_with_real_throughflow() {
    let masses = [0.002, 0.028, 0.044];
    let volumes = [0.001, 0.002, 0.003];
    let chi = matrix(3, &[(0, 1, 0.2), (0, 2, -0.3), (1, 2, 0.1)]);
    let network = network(&masses, chi, &volumes,
        vec![face(0, 1), face(0, 2), face(1, 2)], vec![uniform_matrix(3, 1e-4)]);
    let x = [0.2, 0.3, 0.5];
    let initial = field(&network, &[&x, &x, &x]);
    let flow = circulation(&network, 0.0005);
    let candidate = initial.advect_candidate(&flow, 2.0,
        MulticomponentMiscibilityAdvectionControls { maximum_courant: 0.2, ..Default::default() },
        || false,
    ).unwrap();
    assert_eq!(candidate.step().substeps, 5);
    for cell in 0..3 {
        for species in 0..3 {
            near(candidate.field().cell_component_moles(cell).unwrap()[species], initial.cell_component_moles(cell).unwrap()[species], 2e-14);
            near(candidate.field().cell_mole_fractions(cell).unwrap()[species], x[species], 2e-14);
        }
        let initial_mass: f64 = (0..3).map(|s| initial.cell_component_moles(cell).unwrap()[s] * masses[s]).sum();
        let final_mass: f64 = (0..3).map(|s| candidate.field().cell_component_moles(cell).unwrap()[s] * masses[s]).sum();
        near(final_mass, initial_mass, 2e-15);
        near(candidate.field().cell_component_moles(cell).unwrap().iter().sum(), 1000.0 * volumes[cell], 2e-14);
    }
    // One mol of mixture crosses each edge in two seconds, although each
    // cell's accepted inventory is unchanged. These are physical receipts.
    for (face, sign) in [1.0, -1.0, 1.0].into_iter().enumerate() {
        for species in 0..3 {
            near(candidate.step().transferred_mol[face][species], sign * x[species], 2e-14);
        }
    }
    near(candidate.step().free_energy_change_j, 0.0, 4e-9);
}

#[test]
fn work_divergence_foreign_flow_and_last_poll_refusals_leave_the_owner_unchanged() {
    let network = network(&[0.002, 0.028, 0.044], vec![0.0; 9], &[0.001; 2],
        vec![face(0, 1)], vec![uniform_matrix(3, 1e-4)]);
    let initial = field(&network, &[&[0.2, 0.3, 0.5], &[0.6, 0.15, 0.25]]);
    let saved = initial.clone();
    let flow = opposing(&network, 0.0005);
    assert!(matches!(initial.advect_candidate(&flow, 2.0,
        MulticomponentMiscibilityAdvectionControls { maximum_courant: 0.1, maximum_substeps: 1, ..Default::default() }, || false),
        Err(MulticomponentMiscibilityError::AdvectionWorkLimit { .. })));
    assert!(matches!(network.prepare_advection(vec![MulticomponentMiscibilityAdvectionFace {
        forward_m3_s: 0.0005, reverse_m3_s: 0.0,
    }], 0.0, || false), Err(MulticomponentMiscibilityError::VolumeDivergence { .. })));
    let foreign_network = crate::miscibility::multicomponent::tests::network(
        &[0.002, 0.028, 0.044], uniform_matrix(3, 0.2), &[0.001; 2],
        vec![face(0, 1)], vec![uniform_matrix(3, 1e-4)]);
    let foreign = opposing(&foreign_network, 0.0005);
    assert!(initial.advect_candidate(&foreign, 0.2, MulticomponentMiscibilityAdvectionControls::default(), || false).is_err());
    let mut polls = 0usize;
    let expected = initial.advect_candidate(&flow, 0.2, MulticomponentMiscibilityAdvectionControls::default(), || { polls += 1; false }).unwrap();
    assert!(polls > 3);
    let mut seen = 0usize;
    assert_eq!(initial.advect_candidate(&flow, 0.2, MulticomponentMiscibilityAdvectionControls::default(), || { seen += 1; seen == polls }).unwrap_err(), MulticomponentMiscibilityError::Cancelled);
    assert_eq!(initial, saved);
    assert_eq!(initial.advect_candidate(&flow, 0.2, MulticomponentMiscibilityAdvectionControls::default(), || false).unwrap(), expected);
}
