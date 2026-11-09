//! G1/G3/G4: original thermal beam capture, finite-wall storage and exact vacuum.

mod trajectory;
use super::super::super::tests::{controls, energy, gas, interface, near};
use super::super::thermal::tests::material;
use super::*;
use crate::thermal::transport::volumetric::timeline::FixedVolumeThermalLaw;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        channel::{ClausingControls, CylindricalMolecularGeometry, SpeciesThermalAccommodation},
    },
    tests::{point, species_with_energy_reference, transport_in_domain},
};
use fs_matdb::QueryPoint;
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const AREA: f64 = 1e-6;
const DT: f64 = 0.001;
const M: f64 = 0.028;

fn channel(datum: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(datum),
        CylindricalMolecularGeometry::try_new(
            1e-8,
            2e-8,
            1e-4,
            "synthetic fixed diffuse first-fill passage",
        )
        .unwrap(),
    )
    .unwrap()
}
fn empty(volume: f64) -> EffusionReceiver {
    EffusionReceiver::Evacuated(EvacuatedGasVolume::new(volume, QueryPoint::new()).unwrap())
}
fn link(channel: &DiffuseMolecularChannel) -> FiniteChannelFillLink<'_> {
    FiniteChannelFillLink {
        id: 7,
        donor: 0,
        receiver: 1,
        wall: 9,
        channel,
        area_m2: AREA,
    }
}
fn initial<'a>(candidate: &'a FiniteChannelFillCandidate<'_, '_>) -> &'a FiniteChannelFirstFill {
    match &candidate.step().transfer {
        FiniteChannelFillTransfer::Initial(first) => first,
        FiniteChannelFillTransfer::Exchange(_) => panic!("expected first-fill account"),
    }
}

#[test]
fn first_fill_temperature_follows_transmitted_wall_energy_and_finite_storage() {
    let source = material(1e6, 1000.0, [200.0, 2000.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = empty(0.001);
    let path = channel(0.0);
    let incident = AREA * DT * 100_000.0 / (2.0 * core::f64::consts::PI * M * R * 400.0).sqrt();
    let tau = path.geometry().transmission_probability();
    let direct = path.geometry().direct_probability();
    for tw in [250.0, 900.0] {
        let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, tw).unwrap();
        let candidate = step_candidate(
            &donor,
            &receiver,
            &wall,
            link(&path),
            DT,
            energy(),
            energy(),
            controls(),
            || false,
        )
        .unwrap();
        let filled = candidate.receiver().gas().unwrap();
        let n = tau * incident;
        let expected_t = (4.0 / 3.0) * (direct * 400.0 + (tau - direct) * tw) / tau;
        near(filled.mass_kg(), n * M, 1e-19);
        near(filled.state().state().temperature, expected_t, 3e-8);
        near(
            filled.state().state().pressure,
            n * R * expected_t / 0.001,
            1e-8,
        );
        let wall_heat = -incident * (1.0 - direct) * 2.0 * R * (tw - 400.0);
        near(
            candidate.step().balance.wall_energy_change_j,
            wall_heat,
            1e-9,
        );
        near(
            candidate.step().balance.gas_energy_change_j,
            -wall_heat,
            1e-9,
        );
        assert!(candidate.step().balance.total_entropy_change_j_k > 0.0);
        assert!(
            candidate.step().balance.energy_residual_j.abs()
                <= controls().exchange.energy_tolerance_j
        );
        let report = initial(&candidate);
        assert_eq!(
            report.captured_species_kg[&SpeciesId::new("A").unwrap()],
            report.transport.gas.outlets[0].escaped_species[0].mass_kg
        );
        near(
            report.reference_relative_energy_j,
            n * R * 400.0 / 2.0 + incident * (tau - direct) * 2.0 * R * (tw - 400.0),
            1e-12,
        );
        assert!(matches!(receiver, EffusionReceiver::Evacuated(_)));
        near(wall.temperature_k(), tw, 1e-10);
    }
}

#[test]
fn the_filled_receiver_uses_the_original_reciprocal_finite_wall_owner() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 700.0).unwrap();
    let donor = gas(400.0, 100_000.0, 0.0);
    let path = channel(0.0);
    let first = step_candidate(
        &donor,
        &empty(0.001),
        &wall,
        link(&path),
        DT,
        energy(),
        energy(),
        controls(),
        || false,
    )
    .unwrap();
    let second = step_candidate(
        first.donor(),
        first.receiver(),
        first.wall(),
        link(&path),
        DT,
        energy(),
        energy(),
        controls(),
        || false,
    )
    .unwrap();
    let original = BTreeMap::from([
        (0, first.donor().clone()),
        (1, first.receiver().gas().unwrap().clone()),
    ]);
    let native = thermal::step_candidate(
        &original,
        &BTreeMap::from([(9, first.wall().clone())]),
        GasTransportPaths::default(),
        &[FiniteChannelNetworkLink {
            id: 7,
            parcels: [0, 1],
            channel: &path,
            area_m2: AREA,
            wall: 9,
        }],
        DT,
        &BTreeMap::from([(0, energy()), (1, energy())]),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(second.donor(), &native.parcels()[&0]);
    assert_eq!(second.receiver().gas().unwrap(), &native.parcels()[&1]);
    assert_eq!(
        second.wall().temperature_k(),
        native.walls()[&9].temperature_k()
    );
    match &second.step().transfer {
        FiniteChannelFillTransfer::Exchange(report) => {
            assert_eq!(report, native.step());
            assert!(
                report.channels[0].transmitted.species[0]
                    .right_to_left
                    .mass_rate_kg_s
                    > 0.0
            );
        }
        _ => panic!("a real receiving gas must supply its return beam"),
    }
}

#[test]
fn closed_or_elastic_channels_preserve_empty_and_inactive_wall_semantics() {
    let source = material(1e6, 0.0, [200.0, 2000.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap();
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = empty(0.001);
    let path = channel(0.0);
    let closed = step_candidate(
        &donor,
        &receiver,
        &wall,
        FiniteChannelFillLink {
            area_m2: 0.0,
            ..link(&path)
        },
        DT,
        energy(),
        energy(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(closed.receiver(), &receiver);
    assert_eq!(closed.donor(), &donor);
    assert_eq!(closed.wall().temperature_k(), wall.temperature_k());
    assert_eq!(initial(&closed).captured_mass_kg, 0.0);
    assert!(initial(&closed).fill.is_none());
    let elastic = path
        .with_thermal_accommodation(
            0.0,
            ClausingControls::default(),
            "synthetic elastic diffuse surface",
            || false,
        )
        .unwrap();
    let filled = step_candidate(
        &donor,
        &receiver,
        &wall,
        link(&elastic),
        DT,
        energy(),
        energy(),
        controls(),
        || false,
    )
    .unwrap();
    near(
        filled.receiver().gas().unwrap().state().state().temperature,
        400.0 * 4.0 / 3.0,
        3e-8,
    );
    assert_eq!(filled.wall().temperature_k(), wall.temperature_k());
    assert_eq!(filled.step().balance.wall_energy_change_j, 0.0);
    assert!(
        initial(&filled).transport.gas.outlets[0]
            .wall
            .reflected_species[0]
            .mass_kg
            > 0.0
    );
}

#[test]
fn different_species_accommodation_controls_first_fill_energy_and_fractionation() {
    let rows = [("A", 0.028), ("B", 0.004)];
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter()
                .map(|&(id, m)| species_with_energy_reference(id, m, 2.5, 0.0, 0.0))
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
        )
        .unwrap(),
        rows.iter()
            .map(|&(id, _)| transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1e-8, 1e7]))
            .collect(),
        None,
    )
    .unwrap();
    let donor = ClosedGasParcel::new(model.clone(), point(400.0, 200.0), 0.001).unwrap();
    let mouth = EffusiveInterface::try_new(
        model,
        rows.iter()
            .map(|&(id, _)| {
                EffusiveSpeciesTransmission::try_new(
                    SpeciesId::new(id).unwrap(),
                    1.0,
                    "synthetic open channel mouth",
                )
                .unwrap()
            })
            .collect(),
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1e-8, 1e7),
        },
        EffusiveRegime::try_new(1e-8, 1e-11, 1e-4, "synthetic collisionless support").unwrap(),
    )
    .unwrap();
    let path = DiffuseMolecularChannel::try_new(
        mouth,
        CylindricalMolecularGeometry::try_new(1e-8, 2e-8, 1e-4, "synthetic mixture channel")
            .unwrap(),
    )
    .unwrap()
    .with_species_thermal_accommodation(
        vec![
            SpeciesThermalAccommodation::try_new(
                SpeciesId::new("A").unwrap(),
                0.0,
                "elastic species A",
            )
            .unwrap(),
            SpeciesThermalAccommodation::try_new(
                SpeciesId::new("B").unwrap(),
                1.0,
                "thermal species B",
            )
            .unwrap(),
        ],
        ClausingControls::default(),
        || false,
    )
    .unwrap();
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap();
    let result = step_candidate(
        &donor,
        &empty(0.001),
        &wall,
        link(&path),
        DT,
        energy(),
        energy(),
        controls(),
        || false,
    )
    .unwrap();
    let tau = path.geometry().transmission_probability();
    let direct = path.geometry().direct_probability();
    let common = AREA * DT * 200.0 * tau / (2.0 * core::f64::consts::PI * R * 400.0).sqrt();
    let na = common * 0.25 / 0.028_f64.sqrt();
    let nb = common * 0.75 / 0.004_f64.sqrt();
    let tb = (direct * 400.0 + (tau - direct) * 900.0) / tau;
    let expected_t = (4.0 / 3.0) * (na * 400.0 + nb * tb) / (na + nb);
    near(
        result.receiver().gas().unwrap().state().state().temperature,
        expected_t,
        3e-8,
    );
    near(
        initial(&result).captured_species_kg[&SpeciesId::new("A").unwrap()],
        na * 0.028,
        1e-21,
    );
    near(
        initial(&result).captured_species_kg[&SpeciesId::new("B").unwrap()],
        nb * 0.004,
        1e-21,
    );
    assert!(
        result.step().balance.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j
    );
}

#[test]
fn formation_references_do_not_change_the_first_fill_or_finite_wall() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap();
    let mut reference: Option<[f64; 3]> = None;
    for datum in [0.0, 1e20] {
        let path = channel(datum);
        let result = step_candidate(
            &gas(400.0, 100_000.0, datum),
            &empty(0.001),
            &wall,
            link(&path),
            DT,
            energy(),
            energy(),
            controls(),
            || false,
        )
        .unwrap();
        let values = [
            result.donor().state().state().temperature,
            result.receiver().gas().unwrap().state().state().temperature,
            result.wall().temperature_k(),
        ];
        if let Some(reference) = reference {
            for (value, expected) in values.into_iter().zip(reference) {
                near(value, expected, 1e-8);
            }
        } else {
            reference = Some(values);
        }
    }
}

#[test]
fn late_source_failure_inventory_limits_and_cancellation_are_atomic() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap();
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = empty(0.001);
    let path = channel(0.0);
    let original = donor.clone();
    assert!(
        step_candidate(
            &donor,
            &empty(1e-12),
            &wall,
            link(&path),
            DT,
            energy(),
            energy(),
            controls(),
            || false
        )
        .is_err()
    );
    assert!(matches!(
        step_candidate(
            &donor,
            &receiver,
            &wall,
            link(&path),
            100.0,
            energy(),
            energy(),
            controls(),
            || false
        ),
        Err(Error::Physical(_))
    ));
    let mut count = 0;
    step_candidate(
        &donor,
        &receiver,
        &wall,
        link(&path),
        DT,
        energy(),
        energy(),
        controls(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        step_candidate(
            &donor,
            &receiver,
            &wall,
            link(&path),
            DT,
            energy(),
            energy(),
            controls(),
            || {
                seen += 1;
                seen == count
            }
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(donor, original);
    assert!(matches!(receiver, EffusionReceiver::Evacuated(_)));
    assert_eq!(wall.temperature_k(), 900.0);
}
