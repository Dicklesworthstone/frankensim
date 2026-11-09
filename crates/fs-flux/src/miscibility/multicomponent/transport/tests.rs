//! Same-interval reductions, original free-energy accounting and coupled chronology.

use super::*;
use crate::miscibility::multicomponent::advection::{
    MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionControls,
    MulticomponentMiscibilityAdvectionFace,
};
use crate::miscibility::multicomponent::tests::{
    face, field, gibbs_energy, matrix, near, network, uniform_matrix,
};
use crate::miscibility::multicomponent::{
    MulticomponentMiscibilityError, MulticomponentMiscibilityNetwork,
};

fn opposing(network: &MulticomponentMiscibilityNetwork, rate: f64) -> MulticomponentMiscibilityAdvection {
    network.prepare_advection(vec![MulticomponentMiscibilityAdvectionFace {
        forward_m3_s: rate, reverse_m3_s: rate,
    }], 0.0, || false).unwrap()
}

#[test]
fn zero_flow_is_original_diffusion_and_a_uniform_advected_endpoint_has_no_diffusive_stage() {
    let network = network(&[0.002, 0.028, 0.044], vec![0.0; 9], &[0.001; 2],
        vec![face(0, 1)], vec![uniform_matrix(3, 1e-4)]);
    let initial = field(&network, &[&[0.2, 0.3, 0.5], &[0.6, 0.1, 0.3]]);
    let controls = MulticomponentMiscibilitySolve::default();
    let advection = MulticomponentMiscibilityAdvectionControls::default();
    let rest = opposing(&network, 0.0);
    let diffusion = initial.step_candidate(0.03, controls, || false).unwrap();
    let combined = initial.transport_candidate(&rest, 0.03, advection, controls, || false).unwrap();
    assert_eq!(combined.field(), diffusion.field());
    assert_eq!(&combined.step().diffusion, diffusion.step());
    assert_eq!(combined.step().advection.transferred_mol, vec![vec![0.0; 3]]);
    assert_eq!(combined.field().elapsed_s(), 0.03);

    // An equal half-volume exchange makes both cells exactly uniform. The
    // same positive original D matrix then produces zero chemical transport.
    let stirring = opposing(&network, 0.001);
    let stirred = initial.advect_candidate(&stirring, 0.5, advection, || false).unwrap();
    let combined = initial.transport_candidate(&stirring, 0.5, advection, controls, || false).unwrap();
    assert_eq!(combined.field(), stirred.field());
    assert_eq!(&combined.step().advection, stirred.step());
    assert_eq!(combined.step().diffusion.transferred_mol, vec![vec![0.0; 3]]);
    assert_eq!(combined.field().elapsed_s(), 0.5);
    near(combined.step().stage_energy_residual_j, 0.0, 2e-10);
}

#[test]
fn coupled_transport_uses_the_advected_mixture_and_one_clock_with_complete_receipts() {
    let chi = matrix(3, &[(0, 1, 0.2), (0, 2, -0.3), (1, 2, 0.1)]);
    let network = network(&[0.002, 0.028, 0.044], chi.clone(), &[0.001, 0.002],
        vec![face(0, 1)], vec![matrix(3, &[(0, 1, 1e-6), (0, 2, 1e-4), (1, 2, 1e-4)])]);
    let initial = field(&network, &[&[0.3, 0.3, 0.4], &[0.31, 0.1, 0.59]]);
    let saved = initial.clone();
    let flow = opposing(&network, 0.0005);
    let duration = 0.3;
    let controls = MulticomponentMiscibilitySolve::default();
    let advection = MulticomponentMiscibilityAdvectionControls::default();
    let advected = initial.advect_candidate(&flow, duration, advection, || false).unwrap();
    let composed = advected.field().step_candidate(duration, controls, || false).unwrap();
    let independent_diffusion = initial.step_candidate(duration, controls, || false).unwrap();
    let mut polls = 0usize;
    let combined = initial.transport_candidate(&flow, duration, advection, controls, || { polls += 1; false }).unwrap();
    assert_eq!(combined.field().component_moles(), composed.field().component_moles());
    assert_eq!(combined.field().mole_fractions(), composed.field().mole_fractions());
    assert_eq!(&combined.step().advection, advected.step());
    assert_eq!(&combined.step().diffusion, composed.step());
    assert_eq!(combined.field().elapsed_s(), duration);
    assert_eq!(composed.field().elapsed_s(), 2.0 * duration);
    let mut independent_difference = 0.0_f64;
    for cell in 0..2 {
        for species in 0..3 {
            let index = 3 * cell + species;
            let independent = advected.field().component_moles()[index]
                + independent_diffusion.field().component_moles()[index]
                - initial.component_moles()[index];
            independent_difference = independent_difference.max((combined.field().component_moles()[index] - independent).abs());
            let transferred = combined.step().transferred_mol[0][species];
            let sign = if cell == 0 { -1.0 } else { 1.0 };
            near(combined.field().component_moles()[index] - initial.component_moles()[index], sign * transferred, 2e-14);
        }
        near(combined.field().cell_component_moles(cell).unwrap().iter().sum(), [1.0, 2.0][cell], 2e-14);
    }
    assert!(independent_difference > 1e-5);
    for species in 0..3 {
        near(combined.step().transferred_mol[0][species],
            combined.step().advection.transferred_mol[0][species] + combined.step().diffusion.transferred_mol[0][species], 2e-15);
        near(combined.step().component_residual_mol[species], 0.0, 2e-14);
    }
    let actual_energy = gibbs_energy(combined.field(), 3, &chi) - gibbs_energy(&initial, 3, &chi);
    assert!(actual_energy < 0.0);
    near(combined.step().free_energy_change_j, actual_energy, 4e-9);
    near(combined.step().free_energy_change_j,
        combined.step().advection.free_energy_change_j + combined.step().diffusion.free_energy_change_j, 2e-10);
    near(combined.step().stage_energy_residual_j, 0.0, 2e-10);
    let continuation = combined.field().transport_candidate(&flow, 0.2, advection, controls, || false).unwrap();
    assert_eq!(continuation.field().elapsed_s(), 0.5);
    assert_eq!(initial.transport_candidate(&flow, duration, advection, controls, || false).unwrap(), combined);
    let mut seen = 0usize;
    assert_eq!(initial.transport_candidate(&flow, duration, advection, controls, || { seen += 1; seen == polls }).unwrap_err(), MulticomponentMiscibilityError::Cancelled);
    assert_eq!(initial, saved);
}
