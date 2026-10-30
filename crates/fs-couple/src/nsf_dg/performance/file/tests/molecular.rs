//! Species-aware exterior propagation on the original resolved-flow render.
//! Reservoir data is synthetic: these are G0/G3 integration regressions, not
//! validation of molecular coefficients or a measured instrument's sound.
use super::*;
use crate::air_path::GasAbsorptionModel;
use crate::pcm_wav::baffled::gas::GasPropagationLimits;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

const OBSERVER: &str = "observer_gas molecular 4096 1000000 50000000 declared 2 N2 1 vibration 1 1000 1 power 10 300 0 250 350 O2 0\n";

fn exterior_source() -> String {
    include_str!("../../../../../examples/resolved-dg.performance")
        .replace("output 2048 20", "output 256 20")
        .replace("hold 128 -0.1", "hold 64 -0.1")
        .replace("hold 256 0.1", "hold 128 0.1")
        .replace("hold 2048 0", "hold 256 0")
        .replace(
            "dg_exterior 300 100000 0",
            "dg_exterior mixture 300 100000 mole 2 N2 0.8 O2 0.2",
        )
        .replace("observer_gas stokes-kirchhoff 4096 1000000 50000000\n", "")
}

fn observed(source: &str, row: &str) -> String {
    source.replace("boundary 3", &format!("{row}boundary 3"))
}

#[test]
fn molecular_exterior_matches_shared_owner_and_changes_only_observed_pressure() {
    let source = exterior_source();
    let file = read(&observed(&source, OBSERVER), 17);
    let exterior = file.gas_records()[1].clone();
    let report = file.renderer().observer_gas().unwrap();
    assert_eq!(report.law, GasAbsorptionModel::MolecularRelaxation);
    assert_eq!(&report.gas, exterior.state());
    assert_ne!(&report.gas, file.renderer().system().reference_gas());
    // The resolved volume deliberately has only the legacy dry-air record.
    // Accepting molecular propagation must use the independently declared
    // exterior species, never borrow or infer the volume's composition.
    assert!(matches!(file.gas_records()[0], GasRecord::Legacy(_)));
    assert!(!report.filters.is_empty());

    // Independently construct the same constitutive law, then attach it through
    // the existing typed owner. The file reader must yield identical pressure.
    let law = RelaxationAbsorption::new(vec![
        SpeciesRelaxation {
            species: SpeciesId::new("N2").unwrap(),
            processes: vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1000.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 10.0,
                    t_ref_k: 300.0,
                    exponent: 0.0,
                    t_min_k: 250.0,
                    t_max_k: 350.0,
                },
            }],
        },
        SpeciesRelaxation {
            species: SpeciesId::new("O2").unwrap(),
            processes: vec![],
        },
    ])
    .unwrap();
    let medium = fs_material::relaxation::input::medium_from_record(&exterior).unwrap();
    let propagation = GasPropagation {
        law: GasAbsorptionModel::MolecularRelaxation,
        relaxation: Some(law.evaluate(&medium).unwrap()),
        limits: GasPropagationLimits {
            max_unique_paths: 4096,
            max_filter_taps: 1_000_000,
            max_kernel_terms: 50_000_000,
        },
    };
    let mut direct = read(&source, 31)
        .into_renderer()
        .with_observer_gas(&propagation)
        .unwrap();
    let (pressure, completed) = render(file, 17);
    let mut reference = vec![0.0; pressure.len()];
    for block in reference.chunks_mut(31) {
        direct.block(block).unwrap();
    }
    assert_eq!(
        pressure, reference,
        "file and typed source paths must agree across callbacks"
    );
    assert!(
        pressure.iter().any(|p| p.abs() > 1e-12),
        "receiver must hear the actual source"
    );

    let lossless = observed(&source, "observer_gas lossless 4096 1000000 50000000\n");
    let (comparison, unabsorbed) = render(read(&lossless, 17), 17);
    assert!(
        pressure
            .iter()
            .zip(&comparison)
            .any(|(a, b)| (a - b).abs() > 1e-12),
        "molecular absorption and dispersion must reach the emitted pressure"
    );
    assert_eq!(
        completed.system().zone().state_hash(),
        unabsorbed.system().zone().state_hash()
    );
    assert_eq!(
        completed.system().load().state(),
        unabsorbed.system().load().state()
    );
    assert_eq!(
        completed.system().load().state(),
        direct.system().load().state()
    );
}

#[test]
fn exterior_molecular_source_mismatches_refuse_before_physical_preparation() {
    let source = observed(&exterior_source(), OBSERVER);
    for refused in [
        source.replace("N2 0.8 O2 0.2", "N2 0.8 He 0.2"),
        source.replace(
            "dg_exterior mixture 300 100000 mole 2 N2 0.8 O2 0.2",
            "dg_exterior 300 100000 0",
        ),
        source.replace("power 10 300 0 250 350", "power 10 300 0 310 350"),
        source.replace(
            "declared 2 N2 1 vibration 1 1000 1 power 10 300 0 250 350 O2 0",
            "iso9613-terrestrial-air",
        ),
    ] {
        assert!(
            matches!(
                ResolvedDgFilePerformance::from_bytes(
                    refused.as_bytes(),
                    17,
                    Arc::new(CancelGate::new_clock_free())
                ),
                Err(ResolvedDgInputError::Gas(GasInputError::Material(_)))
            ),
            "mismatched exterior species or constitutive domain must refuse at source admission"
        );
    }
}
