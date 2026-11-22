//! Stream compatibility retains the native material catalog, not just calorics.

use super::*;

fn catalog_words(model: &ThermochemicalGasModel) -> Vec<u64> {
    let mut words = Vec::new();
    model.catalog_identity_words(&mut words);
    words
}

#[test]
fn equal_caloric_streams_refuse_other_transport_selection_and_water_sources() {
    use MaterialPropertySelection::{PreferObservationBacked, SingleClaimOnly};

    let reference = parcel(0.0, 0.0);
    for foreign in [
        parcel_of(0.0, 0.0, 3.0e-5, SingleClaimOnly, false),
        parcel_of(0.0, 0.0, 2.0e-5, PreferObservationBacked, false),
        parcel_of(0.0, 0.0, 2.0e-5, SingleClaimOnly, true),
    ] {
        assert_eq!(
            reference.model().thermodynamics(),
            foreign.model().thermodynamics()
        );
        assert_ne!(
            catalog_words(reference.model()),
            catalog_words(foreign.model())
        );
        let original = streaming(reference.clone());
        let changed = streaming(foreign);
        assert_eq!(original.ports()[2].width(), changed.ports()[2].width());
        assert!(matches!(
            connect_streams(original, changed),
            Err(AdmissionError::StreamBasis { .. })
        ));
    }
}

#[test]
fn stream_catalog_retains_pack_bindings_even_with_identical_physical_cards() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    use fs_material::gas_species::input::packs::{SourcePackBytes, read_source_parcel};

    let sources = [SourcePackBytes {
        species: include_bytes!("../../../../examples/system/n2/species.pack"),
        nasa9: include_bytes!("../../../../examples/system/n2/nasa9.pack"),
        transport: include_bytes!("../../../../examples/system/n2/transport.fsmcdpk"),
    }];
    let reference = ClosedGasParcel::from_isochoric(
        read_source_parcel(
            "thermochemical 300 90000 mole 1 N2 1",
            &sources,
            None,
            ThermochemicalPhase::FrozenGas,
            V0,
        )
        .unwrap(),
    );
    let model = reference.model();
    let stripped = ThermochemicalGasModel::new(
        model.thermodynamics().clone(),
        model.transport_materials().to_vec(),
        model.water_species(),
    )
    .unwrap();
    assert_eq!(model.thermodynamics(), stripped.thermodynamics());
    assert_eq!(model.transport_materials(), stripped.transport_materials());
    assert!(!model.source_packs().is_empty());
    assert!(stripped.source_packs().is_empty());
    assert_ne!(catalog_words(model), catalog_words(&stripped));
    let stripped = ClosedGasParcel::new(stripped, reference.point().clone(), V0).unwrap();
    let cavity = |state| {
        GasCavity::new(
            state,
            GasEnergySolve {
                maximum_temperature_k: 900.0,
                ..controls()
            },
        )
        .unwrap()
        .streaming(1.0e-15)
        .unwrap()
    };
    assert!(matches!(
        connect_streams(cavity(reference), cavity(stripped)),
        Err(AdmissionError::StreamBasis { .. })
    ));
}

#[test]
fn catalog_identity_excludes_composition_and_local_material_axes() {
    let reference = mixture(0.75, 1.0e5);
    let other = mixture(0.25, 1.5e5);
    let other = ClosedGasParcel::new(
        other.model().clone(),
        other.point().clone().with("material-lot", 7.0).unwrap(),
        1.5 * V0,
    )
    .unwrap();
    assert_ne!(reference.model(), other.model());
    assert_ne!(reference.point(), other.point());
    assert_eq!(
        catalog_words(reference.model()),
        catalog_words(other.model())
    );
    connect_streams(streaming(reference), streaming(other)).unwrap();
}
