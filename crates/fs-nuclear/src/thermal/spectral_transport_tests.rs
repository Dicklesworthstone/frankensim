//! G1/G3: spectral radiation follows original burnup events and residence.
//! The shared manufactured evaluation has one-second parent/daughter half-lives;
//! its 100/200 eV lines are equation fixtures, not isotope calibration data.

use super::*;
use crate::burnup::transport::{
    BurnupTransportBudget, BurnupTransportLimits, BurnupTransportState, NuclideTransfer,
    NuclideTransport,
};
use crate::burnup::{
    BurnupBudget, BurnupInventory, BurnupLimits, BurnupModel, BurnupRate, BurnupReaction,
};
use crate::endf::decay::model::EvaluatedDecayAssembly;
use crate::endf::decay::model::emissions::{DecayEmissionLimits, PreparedDecayEmission};
use crate::endf::decay::spectra::DecayRadiation;
use crate::radiation::RadiationDeposit;
use crate::radiation::decay::{DecayRadiationClosure, PreparedDecayRadiation, tests::gamma_chain};
use crate::radiation::grouped::{
    GroupedRadiationResponse, GroupedRadiationResponseLimits, GroupedRadiationSourceResponse,
    RadiationGroupResponse,
};

const EV_J: f64 = 1.602_176_634e-19;
const INITIAL_NUCLEI: f64 = 1e18;

fn partition(assembly: Arc<EvaluatedDecayAssembly>) -> Arc<PreparedDecayRadiation> {
    let limits = DecayEmissionLimits::default();
    let spectrum = Arc::new(
        PreparedDecayEmission::new(
            assembly,
            DecayRadiation::Gamma,
            vec![0.0, 150.0 * EV_J, 300.0 * EV_J],
            limits,
            || Ok(()),
        )
        .unwrap(),
    );
    Arc::new(
        PreparedDecayRadiation::new(
            spectrum,
            DecayRadiationClosure {
                absolute_j_per_decay: 0.0,
                relative: 1e-12,
            },
            limits,
            || Ok(()),
        )
        .unwrap(),
    )
}

fn inventory_model(assembly: &EvaluatedDecayAssembly, capture: bool) -> Arc<BurnupModel> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured spectral burnup grid".into(), vec![0.0, 1.0]).unwrap(),
    );
    let mut data = assembly
        .burnup_data(groups, 200.0, 1000.0, || Ok(()))
        .unwrap();
    // Binding must follow the original physical channels, not a reaction index.
    data.reactions.reverse();
    if capture {
        data.reactions.insert(
            0,
            BurnupReaction {
                data_id: "independent prompt capture radiation".into(),
                parent: 0,
                rate: BurnupRate::Capture {
                    cross_sections_m2: vec![core::f64::consts::LN_2],
                },
                products: vec![],
                energy_per_event: EnergyChannels {
                    deposited_j: 0.0,
                    radiation_j: 7.0 * EV_J,
                    escaped_j: 0.0,
                },
            },
        );
    }
    Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap())
}

fn source_budget() -> BurnupHeatingBudget {
    BurnupHeatingBudget {
        source: BurnupBudget {
            limits: BurnupLimits::default(),
            nuclei_tolerance: 1e6,
            energy_tolerance_j: 1e-8,
        },
        max_composition_feedback_terms: 0,
    }
}

fn response(
    partition: &PreparedDecayRadiation,
    receiver: u64,
    rows: &[(u64, [f64; 2])],
) -> Arc<GroupedRadiationResponse> {
    Arc::new(
        GroupedRadiationResponse::new(
            "manufactured location-and-energy deposition fractions".into(),
            partition.radiation(),
            partition.group_edges_j().to_vec(),
            rows.iter()
                .map(|&(source_id, fractions)| GroupedRadiationSourceResponse {
                    source_id,
                    groups: fractions
                        .into_iter()
                        .map(|fraction| RadiationGroupResponse {
                            deposits: vec![RadiationDeposit {
                                receiver_id: receiver,
                                fraction,
                            }],
                            escape_fraction: 1.0 - fraction,
                        })
                        .collect(),
                })
                .collect(),
            GroupedRadiationResponseLimits {
                max_sources: 2,
                max_groups: 2,
                max_receivers: 1,
                max_source_groups: 4,
                max_transfers: 4,
            },
        )
        .unwrap(),
    )
}

#[test]
fn g1_competing_capture_keeps_prompt_radiation_outside_evaluated_decay_groups() {
    let assembly = gamma_chain(100.0, 200.0, 0.5, 0.25);
    let partition = partition(Arc::clone(&assembly));
    let model = inventory_model(&assembly, true);
    let flux = GroupFlux::new(Arc::clone(model.groups()), vec![1.0]).unwrap();
    let source = BurnupCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventory: BurnupInventory::new(model, vec![INITIAL_NUCLEI, 0.0, 0.0]).unwrap(),
        fluxes: NeutronFluxSchedule::new(flux, vec![]).unwrap(),
        composition_feedback: vec![],
    };
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating = BurnupHeating::new(&initial, vec![source], source_budget())
        .unwrap()
        .with_grouped_decay_radiation(
            &initial,
            response(&partition, 2, &[(1, [0.75, 0.25])]),
            vec![DecayRadiationEmitter {
                cell_id: CellId(1),
                partition,
            }],
            vec![RadiationReceiver {
                cell_id: CellId(2),
                material_id: MaterialId(2),
            }],
            || Ok(()),
        )
        .unwrap();
    let (heated, accepted, thermal_audit) = heating
        .advance_to(&initial, 1.0, controls(2.0), || Ok(()))
        .unwrap();

    // Parent loss rate is 2 ln(2); half its removals are captures. Daughter
    // decay rate is ln(2), giving P=N/4, D=N/4 and N/8 daughter decays at t=1.
    let parent_decays = INITIAL_NUCLEI * 0.375;
    let daughter_decays = INITIAL_NUCLEI * 0.125;
    let prompt_radiation = parent_decays * 7.0 * EV_J;
    let low = parent_decays * 50.0 * EV_J;
    let high = daughter_decays * 50.0 * EV_J;
    let local = parent_decays * 50.0 * EV_J + daughter_decays * 150.0 * EV_J;
    let received = 0.75 * low + 0.25 * high;
    let account = accepted.audit();
    close(account.captures, parent_decays, 1e5);
    close(account.decays, parent_decays + daughter_decays, 1e5);
    close(account.released.deposited_j, local, 1e-9);
    close(
        account.released.radiation_j,
        low + high + prompt_radiation,
        1e-9,
    );
    close(account.radiation_deposited_j, received, 1e-9);
    close(account.radiation_ungrouped_j, prompt_radiation, 1e-9);
    close(
        account.radiation_escaped_j,
        0.25 * low + 0.75 * high + prompt_radiation,
        1e-9,
    );
    close(account.thermal_credit_j, local + received, 1e-9);
    close(thermal_audit.external_heat_j, local + received, 1e-9);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + local / 8.0,
        1e-9,
    );
    close(
        heated.temperature_k(CellId(2)).unwrap(),
        300.0 + received / 6.0,
        1e-9,
    );
    let population = accepted.source(CellId(1)).unwrap().inventory.populations();
    for (&actual, fraction) in population.iter().zip([0.25, 0.25, 0.125]) {
        close(actual, INITIAL_NUCLEI * fraction, 1e5);
    }
}

#[test]
fn g1_g3_movement_keeps_spectral_deposition_at_the_actual_decay_location() {
    let assembly = gamma_chain(100.0, 200.0, 0.5, 0.25);
    let partition = partition(Arc::clone(&assembly));
    let model = inventory_model(&assembly, false);
    let limits = BurnupTransportLimits::default();
    let movement_rate = 0.8;
    let state = Arc::new(
        BurnupTransportState::new(
            vec![
                BurnupInventory::new(Arc::clone(&model), vec![INITIAL_NUCLEI, 0.0, 0.0]).unwrap(),
                BurnupInventory::new(Arc::clone(&model), vec![0.0; 3]).unwrap(),
            ],
            vec![0.01, 0.01],
            NuclideTransport {
                transfers: vec![NuclideTransfer {
                    from_compartment: 0,
                    to_compartment: 1,
                    rates_per_s: vec![movement_rate; 3],
                }],
                outflows: vec![],
            },
            limits,
            || Ok(()),
        )
        .unwrap(),
    );
    let fluxes = (0..2)
        .map(|_| {
            NeutronFluxSchedule::new(
                GroupFlux::new(Arc::clone(model.groups()), vec![0.0]).unwrap(),
                vec![],
            )
            .unwrap()
        })
        .collect();
    // Compartment order deliberately opposes thermal identity order.
    let sources = BurnupTransportCellSource::bind(
        state,
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
        fluxes,
    )
    .unwrap();
    let initial = thermal(
        vec![
            cell(20, 1, 2.0, 1200.0, 0.0),
            cell(10, 2, 3.0, 600.0, 0.0),
            cell(30, 1, 1.0, 1200.0, 0.0),
        ],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating = BurnupTransportHeating::new(
        &initial,
        sources,
        BurnupTransportBudget {
            limits,
            nuclei_tolerance: 1e6,
            energy_tolerance_j: 1e-8,
        },
    )
    .unwrap()
    .with_grouped_decay_radiation(
        &initial,
        response(&partition, 30, &[(20, [0.2, 0.8]), (10, [0.9, 0.1])]),
        [10, 20]
            .into_iter()
            .map(|id| DecayRadiationEmitter {
                cell_id: CellId(id),
                partition: Arc::clone(&partition),
            })
            .collect(),
        vec![RadiationReceiver {
            cell_id: CellId(30),
            material_id: MaterialId(1),
        }],
        || Ok(()),
    )
    .unwrap();

    let duration = 2.0;
    let (heated, accepted, thermal_audit) = heating
        .advance_to(&initial, duration, controls(duration), || Ok(()))
        .unwrap();
    let lambda = core::f64::consts::LN_2;
    let total_rate = lambda + movement_rate;
    let remaining = (-lambda * duration).exp();
    let remains_left = (-total_rate * duration).exp();
    let parent_left = INITIAL_NUCLEI * lambda / total_rate * (1.0 - remains_left);
    let daughter_left = INITIAL_NUCLEI
        * (lambda / total_rate).powi(2)
        * (1.0 - (1.0 + total_rate * duration) * remains_left);
    let parent_total = INITIAL_NUCLEI * (1.0 - remaining);
    let daughter_total = INITIAL_NUCLEI * (1.0 - (1.0 + lambda * duration) * remaining);
    let parent_right = parent_total - parent_left;
    let daughter_right = daughter_total - daughter_left;
    let local_left = (50.0 * parent_left + 150.0 * daughter_left) * EV_J;
    let local_right = (50.0 * parent_right + 150.0 * daughter_right) * EV_J;
    let received = 50.0
        * EV_J
        * (0.2 * parent_left + 0.8 * daughter_left + 0.9 * parent_right + 0.1 * daughter_right);
    let outward = 50.0 * EV_J * (parent_total + daughter_total);
    let account = accepted.audit();
    close(
        heated.temperature_k(CellId(20)).unwrap(),
        300.0 + local_left / 8.0,
        1e-9,
    );
    close(
        heated.temperature_k(CellId(10)).unwrap(),
        300.0 + local_right / 6.0,
        1e-9,
    );
    close(
        heated.temperature_k(CellId(30)).unwrap(),
        300.0 + received / 4.0,
        1e-9,
    );
    close(account.decays, parent_total + daughter_total, 1e5);
    close(account.radiation_deposited_j, received, 1e-9);
    close(account.radiation_escaped_j, outward - received, 1e-9);
    close(account.radiation_ungrouped_j, 0.0, 1e-10);
    close(
        thermal_audit.external_heat_j,
        local_left + local_right + received,
        1e-9,
    );
    assert!(Arc::ptr_eq(
        accepted.source(CellId(20)).unwrap().state(),
        accepted.source(CellId(10)).unwrap().state(),
    ));

    let (subdivided, divided_source, _) = heating
        .advance_to(&initial, duration, controls(0.25), || Ok(()))
        .unwrap();
    for id in [10, 20, 30] {
        close(
            subdivided.temperature_k(CellId(id)).unwrap(),
            heated.temperature_k(CellId(id)).unwrap(),
            1e-9,
        );
    }
    close(divided_source.audit().radiation_deposited_j, received, 1e-9);
    close(
        divided_source.audit().radiation_escaped_j,
        outward - received,
        1e-9,
    );
    assert!(divided_source.audit().substeps > account.substeps);
}
