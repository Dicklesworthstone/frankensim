//! Exercise the original nozzle, source merge and caloric reconstruction owners.
use super::*;
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use crate::thermochemical_gas::ThermochemicalGasModel;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_thermochem::{Composition, CompositionBasis};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn parcel(name: &str, molar_mass: f64, pressure: f64, slope: f64) -> ClosedGasParcel {
    reference_parcel(name, molar_mass, pressure, slope, 0.0)
}
fn reference_parcel(name: &str, molar_mass: f64, pressure: f64, slope: f64, reference: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference(name, molar_mass, 3.5, slope, reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport(name, 1.7e-5, 0.025)], None).unwrap();
    ClosedGasParcel::new(model, point(500.0, pressure).with("cell-lot", pressure).unwrap(), 0.001).unwrap()
}
fn links() -> [GasNetworkLink; 2] {
    let energy = GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 };
    [(7, [30, 10]), (11, [10, 90])].map(|(id, parcels)| GasNetworkLink {
        id, parcels, nozzle: GasNozzle { area_m2: 1e-7, discharge_coefficient: 0.8 },
        solve: GasNozzleStepSolve { nozzle: energy, transfer: GasMassTransferSolve {
            donor: energy, receiver: energy, mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 },
            maximum_donor_fraction: 0.05 },
    })
}
fn controls() -> GasNetworkSolve {
    GasNetworkSolve { maximum_parcels: 3, maximum_links: 3, maximum_withdrawn_fraction: 0.05,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-7 }
}
fn initial() -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([(30, parcel("A", 0.028, 300000.0, 0.0)),
        (10, parcel("B", 0.004, 200000.0, 0.001)), (90, parcel("C", 0.040, 100000.0, 0.0))])
}

#[test]
fn g3_sequential_network_energy_balance_is_invariant_to_independent_species_datums() {
    let ordinary = initial();
    let shifted = BTreeMap::from([
        (30, reference_parcel("A", 0.028, 300000.0, 0.0, 1e24)),
        (10, reference_parcel("B", 0.004, 200000.0, 0.001, -2e24)),
        (90, reference_parcel("C", 0.040, 100000.0, 0.0, 3e24)),
    ]);
    let mut edges = links().to_vec();
    edges.push(GasNetworkLink { id: 20, parcels: [30, 90], ..edges[0] });
    let expected = step_candidate(&ordinary, &edges, 1e-4, controls(), || false).unwrap();
    let actual = step_candidate(&shifted, &edges, 1e-4, controls(), || false).unwrap();
    let mut physical_balance = 0.0;
    let mut rounded_absolute_balance = 0.0;
    for (&id, new) in actual.parcels() {
        let original = &ordinary[&id];
        let reference = &expected.parcels()[&id];
        assert_eq!(new.state().state(), reference.state().state());
        assert_eq!(new.mass_kg(), reference.mass_kg());
        physical_balance += reference.mass_kg() * reference.state().specific_internal_energy_j_kg()
            - original.mass_kg() * original.state().specific_internal_energy_j_kg();
        let old = &shifted[&id];
        rounded_absolute_balance += new.mass_kg()
            * (new.state().specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg())
            + (new.mass_kg() - old.mass_kg()) * old.state().specific_internal_energy_j_kg();
    }
    assert!(physical_balance.abs() <= controls().energy_tolerance_j);
    assert!(rounded_absolute_balance.abs() > controls().energy_tolerance_j,
        "absolute formation-energy subtraction cannot resolve this network balance");
    assert_eq!(actual.step().total_energy_residual_j, expected.step().total_energy_residual_j);
    assert!(actual.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(actual.step().total_mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    assert!(actual.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    for (a, b) in actual.step().links.iter().zip(&expected.step().links) {
        assert_eq!(a.flow.mass_flow_kg_s, b.flow.mass_flow_kg_s);
        assert_eq!(a.flow.velocity_m_s, b.flow.velocity_m_s);
        assert_eq!(a.transfer.transferred_mass_kg, b.transfer.transferred_mass_kg);
    }
    edges.reverse();
    assert_eq!(actual, step_candidate(&shifted, &edges, 1e-4, controls(), || false).unwrap());
}
#[test]
fn one_link_is_the_exact_original_pair_and_zero_links_preserve_all_parcels() {
    let old = initial(); let link = links()[0]; let dt = 1e-4;
    let (pair, flow) = old[&30].through_nozzle_to(&old[&10], link.nozzle, dt, link.solve, || false).unwrap();
    let result = step_candidate(&old, &[link], dt, controls(), || false).unwrap();
    assert_eq!(&result.parcels()[&30], pair.donor());
    assert_eq!(&result.parcels()[&10], pair.receiver());
    assert_eq!(&result.step().links[0].transfer, pair.step());
    assert_eq!(result.step().links[0].flow, flow);
    assert_eq!(result.parcels()[&90], old[&90]);
    let idle = step_candidate(&old, &[], dt, controls(), || false).unwrap();
    assert_eq!(idle.parcels(), &old);
    assert_eq!(idle.step().total_energy_residual_j, 0.0);
    assert_eq!(idle.step().maximum_species_mass_residual_kg, 0.0);
}
#[test]
fn branching_and_transit_keep_sources_species_energy_and_canonical_edge_order() {
    let old = initial(); let mut edges = links().to_vec();
    // A branch and a transit: parcel 10 receives and supplies on the same interval.
    edges.push(GasNetworkLink { id: 20, parcels: [30, 90], ..edges[0] });
    let result = step_candidate(&old, &edges, 1e-4, controls(), || false).unwrap();
    edges.reverse();
    assert_eq!(result, step_candidate(&old, &edges, 1e-4, controls(), || false).unwrap());
    assert!(result.parcels()[&90].state().components().len() == 3);
    assert_ne!(result.parcels()[&90].state().state().sound_speed, old[&90].state().state().sound_speed);
    for (&id, next) in result.parcels() {
        assert_eq!(next.volume_m3(), old[&id].volume_m3());
        assert_eq!(next.point().axes().get("cell-lot"), old[&id].point().axes().get("cell-lot"));
        for source in next.model().thermodynamics().components() {
            let original = old.values().find_map(|p| p.model().thermodynamics().components().iter()
                .find(|s| s.species() == source.species())).unwrap();
            assert_eq!(source, original);
        }
    }
    assert!(result.step().total_mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    assert!(result.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
    assert!(result.step().total_energy_residual_j.abs() <= controls().energy_tolerance_j);
}
#[test]
fn parallel_paths_cannot_hide_aggregate_withdrawal_and_reverse_flow_uses_actual_donor() {
    let old = initial(); let first = links()[0];
    let edges = [first, GasNetworkLink { id: 12, parcels: first.parcels, ..first }];
    let one = step_candidate(&old, &[first], 1e-4, controls(), || false).unwrap();
    let bound = 1.5 * one.step().links[0].transfer.transferred_mass_kg / old[&30].mass_kg();
    assert!(matches!(step_candidate(&old, &edges, 1e-4,
        GasNetworkSolve { maximum_withdrawn_fraction: bound, ..controls() }, || false),
        Err(GasNetworkError::Withdrawal { id: 30, .. })));
    let reversed = GasNetworkLink { parcels: [10,30], ..first };
    let next = step_candidate(&old, &[reversed], 1e-4, controls(), || false).unwrap();
    assert_eq!(next.parcels(), one.parcels());
    assert_eq!(next.step().links[0].donor, 30);
}
#[test]
fn late_source_failure_cancellation_bad_maps_and_closed_paths_are_atomic() {
    let old = initial(); let saved = old.clone();
    let mut polls = 0;
    step_candidate(&old, &links(), 1e-4, controls(), || { polls += 1; false }).unwrap();
    let mut seen = 0;
    assert!(matches!(step_candidate(&old, &links(), 1e-4, controls(), || { seen += 1; seen == polls }),
        Err(GasNetworkError::Cancelled)));
    let link = links()[0];
    for edges in [vec![link, link], vec![GasNetworkLink { parcels: [30,30], ..link }],
        vec![GasNetworkLink { parcels: [30,1234], ..link }]] {
        assert!(step_candidate(&old, &edges, 1e-4, controls(), || false).is_err());
    }
    let mut conflict = old.clone();
    conflict.insert(90, parcel("A", 0.028, 100000.0, 0.002));
    let snapshot = conflict.clone();
    // Edge 7 succeeds locally, then edge 11 refuses the conflicting A source.
    assert!(matches!(step_candidate(&conflict, &links(), 1e-4, controls(), || false),
        Err(GasNetworkError::Link { id: 11, .. })));
    assert_eq!(conflict, snapshot);
    let closed = GasNetworkLink { nozzle: GasNozzle { area_m2: 0.0, ..link.nozzle }, ..link };
    assert_eq!(step_candidate(&old, &[closed], 1e-4, controls(), || false).unwrap().parcels(), &old);
    assert_eq!(old, saved);
}
