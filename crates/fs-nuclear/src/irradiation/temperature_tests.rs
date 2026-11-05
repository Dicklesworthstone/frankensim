//! Manufactured temperature laws exercise original irradiation and activation
//! owners. Values are analytical fixtures, not an evaluated isotope library.
use super::*;
use crate::cross_sections::{CrossSectionData, CrossSectionLimits};
use crate::decay::activation::{CaptureActivationBudget, CaptureActivationState, CaptureProduct};
use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
use fs_la::decay::DecayLimits;
use fs_thermochem::decay::DecayLawV1;

fn channels(deposited_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j: 0.5 * deposited_j,
        escaped_j: 0.25 * deposited_j,
    }
}
fn budget() -> IrradiationBudget {
    IrradiationBudget {
        max_neutron_groups: 4,
        source: FissionBudget {
            max_decay_groups: 4,
            energy_tolerance_j: 1e-8,
            fuel_tolerance_nuclei: 1e-8,
        },
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-10, "{a:?} != {b:?}");
}

fn table(
    groups: Arc<NeutronGroups>,
    nuclide: NuclideId,
    reaction: NeutronReaction,
    values: Vec<Vec<f64>>,
) -> Arc<CrossSectionTable> {
    Arc::new(
        CrossSectionTable::new(
            CrossSectionData {
                data_id: format!("manufactured-temperature-{reaction:?}"),
                nuclide,
                reaction,
                groups,
                temperatures_k: vec![300.0, 600.0],
                cross_sections_m2: values,
            },
            CrossSectionLimits {
                max_temperatures: 4,
                max_neutron_groups: 4,
                max_values: 16,
            },
        )
        .unwrap(),
    )
}

fn data() -> TabulatedIrradiationData {
    let groups =
        Arc::new(NeutronGroups::new("manufactured-T-grid".into(), vec![0.0, 1.0, 2.0]).unwrap());
    let nuclide = NuclideId::new(92, 235, 0).unwrap();
    TabulatedIrradiationData {
        data_id: "manufactured-T-parent".into(),
        fission: table(
            groups.clone(),
            nuclide,
            NeutronReaction::Fission,
            vec![vec![0.1, 0.2], vec![0.4, 0.5]],
        ),
        capture: table(
            groups,
            nuclide,
            NeutronReaction::Capture,
            vec![vec![0.6, 0.3], vec![0.2, 0.1]],
        ),
        temperature_min_k: 300.0,
        temperature_max_k: 600.0,
        fission_model: Arc::new(
            FissionModel::new(
                "manufactured-T-prompt-and-tail".into(),
                channels(3.0),
                vec![crate::DecayHeatGroup {
                    decay_constant_s_inv: 0.7,
                    energy_per_fission: channels(1.0),
                }],
            )
            .unwrap(),
        ),
        capture_energy_per_event: channels(2.0),
    }
}

fn fixed_at(data: &TabulatedIrradiationData, temperature: f64) -> Arc<IrradiationModel> {
    Arc::new(
        IrradiationModel::new(IrradiationData {
            data_id: "independent-frozen-coefficients".into(),
            nuclide: data.fission.nuclide(),
            groups: data.fission.groups().clone(),
            cross_sections: (0..data.fission.groups().len())
                .map(|g| GroupCrossSections {
                    fission_m2: data.fission.cross_section_at_group(temperature, g).unwrap(),
                    capture_m2: data.capture.cross_section_at_group(temperature, g).unwrap(),
                })
                .collect(),
            temperature_min_k: 300.0,
            temperature_max_k: 600.0,
            fission_model: data.fission_model.clone(),
            capture_energy_per_event: data.capture_energy_per_event,
        })
        .unwrap(),
    )
}

#[test]
fn temperature_changes_actual_fuel_competition_and_the_original_decay_heat_convolution() {
    let data = data();
    let model = Arc::new(IrradiationModel::new_tabulated(data.clone()).unwrap());
    let flux = GroupFlux::new(model.groups().clone(), vec![1.0, 2.0]).unwrap();
    assert_eq!(model.cross_section_reference_temperature_k(), Some(300.0));
    assert!(model.temperature_tables().is_some());
    let mut counts = Vec::new();
    for (temperature, kf, kc) in [(300.0, 0.5, 1.2), (450.0, 0.95, 0.8), (600.0, 1.4, 0.4)] {
        let rates = model
            .rate_constants(&flux, temperature, budget(), || Ok(()))
            .unwrap();
        close(rates.fission_s_inv, kf);
        close(rates.capture_s_inv, kc);
        let initial = IrradiationState::new(model.clone(), 100.0, vec![2.0]).unwrap();
        let (next, step) = initial
            .trial(0.5, &flux, temperature, budget(), || Ok(()))
            .unwrap();
        let removed = 100.0 * -(-(kf + kc) * 0.5_f64).exp_m1();
        close(step.fissions, removed * kf / (kf + kc));
        close(step.captures, removed * kc / (kf + kc));
        close(
            next.remaining_parent_nuclei(),
            100.0 * (-(kf + kc) * 0.5_f64).exp(),
        );
        let reference =
            IrradiationState::new(fixed_at(&data, temperature), 100.0, vec![2.0]).unwrap();
        let (reference, expected) = reference
            .trial(0.5, &flux, temperature, budget(), || Ok(()))
            .unwrap();
        close(step.released.deposited_j, expected.released.deposited_j);
        close(step.released.radiation_j, expected.released.radiation_j);
        close(next.delayed_energy_j()[0], reference.delayed_energy_j()[0]);
        counts.push((step.fissions, step.captures));
    }
    assert!(
        counts
            .windows(2)
            .all(|x| x[1].0 > x[0].0 && x[1].1 < x[0].1)
    );
}

fn activated(model: Arc<IrradiationModel>) -> CaptureActivationState {
    let parent = IrradiationState::new(model, 100.0, vec![2.0]).unwrap();
    let products = Arc::new(
        DecayModel::new(
            "manufactured-activation-T".into(),
            vec![
                DecayNuclide {
                    id: "capture-product".into(),
                    law: DecayLawV1::from_decay_constant(fs_qty::Frequency::new(0.4)).unwrap(),
                },
                DecayNuclide {
                    id: "stable-daughter".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: channels(1.5),
            }],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    CaptureActivationState::new(
        parent,
        DecayInventory::new(products, vec![0.0, 0.0]).unwrap(),
        CaptureProduct {
            node: 0,
            nuclide: NuclideId::new(92, 236, 0).unwrap(),
            data_id: "manufactured-capture-binding".into(),
        },
    )
    .unwrap()
}

#[test]
fn existing_capture_activation_uses_the_temperature_selected_parent_history() {
    let data = data();
    let model = Arc::new(IrradiationModel::new_tabulated(data.clone()).unwrap());
    let flux = GroupFlux::new(model.groups().clone(), vec![1.0, 2.0]).unwrap();
    let budget = CaptureActivationBudget {
        irradiation: budget(),
        decay: DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
    };
    let initial = activated(model);
    for temperature in [300.0, 450.0, 600.0] {
        let (candidate, step) = initial
            .trial(0.5, &flux, temperature, budget, || Ok(()))
            .unwrap();
        let (reference, expected) = activated(fixed_at(&data, temperature))
            .trial(0.5, &flux, temperature, budget, || Ok(()))
            .unwrap();
        close(step.irradiation.captures, expected.irradiation.captures);
        close(step.released.deposited_j, expected.released.deposited_j);
        for (a, b) in candidate
            .products()
            .populations()
            .iter()
            .zip(reference.products().populations())
        {
            close(*a, *b);
        }
    }
}

#[test]
fn tabulated_parent_rejects_identity_domain_and_budget_mismatches_before_advance() {
    let original = data();
    let mut wrong = original.clone();
    wrong.capture = original.fission.clone();
    assert!(IrradiationModel::new_tabulated(wrong).is_err());
    let mut wrong = original.clone();
    wrong.capture = table(
        original.capture.groups().clone(),
        NuclideId::new(92, 238, 0).unwrap(),
        NeutronReaction::Capture,
        vec![vec![0.2, 0.2], vec![0.3, 0.3]],
    );
    assert!(IrradiationModel::new_tabulated(wrong).is_err());
    let mut wrong = original.clone();
    wrong.temperature_max_k = 601.0;
    assert!(IrradiationModel::new_tabulated(wrong).is_err());
    let mut wrong = original.clone();
    let other = Arc::new(
        NeutronGroups::new("different-grid-identity".into(), vec![0.0, 1.0, 2.0]).unwrap(),
    );
    wrong.capture = table(
        other,
        original.capture.nuclide(),
        NeutronReaction::Capture,
        vec![vec![0.2, 0.2], vec![0.3, 0.3]],
    );
    assert!(IrradiationModel::new_tabulated(wrong).is_err());
    let model = Arc::new(IrradiationModel::new_tabulated(original).unwrap());
    let flux = GroupFlux::new(model.groups().clone(), vec![1.0, 2.0]).unwrap();
    let mut state = IrradiationState::new(model, 100.0, vec![2.0]).unwrap();
    let before = state.clone();
    assert!(
        state
            .advance(1.0, &flux, 299.0, budget(), || Ok(()))
            .is_err()
    );
    let mut small = budget();
    small.max_neutron_groups = 1;
    assert!(matches!(
        state.advance(1.0, &flux, 400.0, small, || Ok(())),
        Err(FissionError::Budget(_))
    ));
    let mut polls = 0;
    let failure = state.advance(1.0, &flux, 400.0, budget(), || {
        polls += 1;
        if polls == 4 {
            Err(FissionError::Budget("caller-stop"))
        } else {
            Ok(())
        }
    });
    assert_eq!(failure, Err(FissionError::Budget("caller-stop")));
    assert_eq!(state, before);
}
