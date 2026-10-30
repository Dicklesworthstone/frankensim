use super::super::{ObservationPolicy, ReedPerformance};
use super::*;

// All reservoir coefficients in this example are synthetic causal test data.
// The fixture makes no measured nitrogen/oxygen or instrument-fidelity claim.
const CLI_INPUT: &str =
    include_str!("../../../../../examples/reed-molecular-synthetic.performance");
// Keep targeted library checks inexpensive. The shipped example uses the
// existing CLI's 48 kHz output clock; these explicit data edits retain its
// duration and pressure schedule at 8 kHz with the corresponding wall band.
static INPUT: std::sync::LazyLock<String> = std::sync::LazyLock::new(|| {
    CLI_INPUT
        .replace("audio 48000 9606", "audio 8000 1601")
        .replace("diffusive 32 0.0001 160000", "diffusive 16 0.0001 30000")
});
const AMBIENT: &str = "ambient mixture 293.15 101325 mole 2 N2 0.7 O2 0.3";

fn render(mut source: ReedPerformance, block: usize) -> Vec<f64> {
    let mut pressure = vec![0.0; source.info().samples as usize];
    for part in pressure.chunks_mut(block) {
        source.pressure_block(part).unwrap();
    }
    assert_eq!(source.remaining_samples(), 0);
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
fn reservoir_data_alone_reaches_physical_pressure_with_exact_block_replay() {
    let faster = INPUT.replace("power 2.0265", "power 0.20265");
    let baseline = ReedPerformance::from_bytes(INPUT.as_bytes(), 37).unwrap();
    let changed = ReedPerformance::from_bytes(faster.as_bytes(), 37).unwrap();
    let a = baseline.molecular_driving_point().unwrap();
    let b = changed.molecular_driving_point().unwrap();
    assert_eq!(a.geometry(), b.geometry());
    assert_eq!(a.inlet_gas(), baseline.ambient_gas());
    assert_eq!(a.inlet_gas(), b.inlet_gas());
    assert_eq!(a.one_way_wavefront_s(), b.one_way_wavefront_s());
    assert_ne!(
        a.material().state().processes(),
        b.material().state().processes()
    );
    assert_eq!(baseline.outlet_receiver().unwrap().gas, *a.exterior_gas());
    assert_eq!(
        baseline.info().compiled_controls,
        changed.info().compiled_controls
    );
    let pa = render(baseline, 37);
    let pb = render(changed, 37);
    let norm: f64 = pa.iter().map(|p| p * p).sum();
    let difference: f64 = pa.iter().zip(&pb).map(|(p, q)| (p - q).powi(2)).sum();
    assert!(
        norm > 0.0 && difference > 1e-8 * norm,
        "reservoir data must change physical Pa: norm={norm}, difference={difference}"
    );
    assert_eq!(
        pa,
        render(
            ReedPerformance::from_bytes(INPUT.as_bytes(), 113).unwrap(),
            113
        )
    );
}

#[test]
fn the_declared_exterior_is_shared_by_molecular_loading_and_the_receiver() {
    let input = INPUT.replace(
        "baffled_outlet 0 0 1 2 8 800",
        "baffled_outlet 0 0 1 2 8 800\nexterior mixture 320 101325 mole 2 N2 0.7 O2 0.3",
    );
    let source = ReedPerformance::from_bytes_with_policy(
        input.as_bytes(),
        37,
        ObservationPolicy::Acceptance,
    )
    .unwrap();
    let point = source.molecular_driving_point().unwrap();
    assert_eq!(point.exterior_gas(), &source.outlet_receiver().unwrap().gas);
    assert_eq!(point.inlet_gas(), source.ambient_gas());
    assert_ne!(point.inlet_gas(), point.exterior_gas());
    assert_eq!(point.exterior_gas().temperature, 320.0);
}

#[test]
fn unsupported_sources_and_missing_species_refuse_instead_of_dropping_the_law() {
    let segment_gas = INPUT.replace(
        "cylinder 0.001 0.15",
        "cylinder 0.001 0.15 gas mixture 293.15 101325 mole 2 N2 0.7 O2 0.3",
    );
    assert!(matches!(ReedPerformance::from_bytes(segment_gas.as_bytes(), 37),
        Err(ReedPerformanceError::Input { what, .. }) if what.contains("homogeneous ambient species")));
    let field = crate::driving_point::regional::field::tests::field_with_viscosity(
        280.0,
        [1.33e-5, 4.15e-5],
    );
    assert!(
        matches!(ReedPerformance::from_bytes_with_gas_field(INPUT.as_bytes(), 37,
        ObservationPolicy::Diagnostic, None, field, 8),
        Err(ReedPerformanceError::Input { what, .. }) if what.contains("explicit continuous"))
    );
    for input in [
        INPUT.replace(AMBIENT, "ambient 293.15 101325 0"),
        INPUT.replace("O2 1 vibration", "CO2 1 vibration"),
        INPUT.replace(
            "power 2.0265 293.15 0 200 500",
            "power 2.0265 293.15 0 310 500",
        ),
    ] {
        let error = ReedPerformance::from_bytes(input.as_bytes(), 37)
            .err()
            .unwrap();
        assert!(
            matches!(error, ReedPerformanceError::Gas { line: 10, .. }),
            "{error}"
        );
        assert!(std::error::Error::source(&error).is_some());
    }
}

#[test]
fn molecular_resolution_and_source_bands_are_explicit_and_bounded() {
    for (before, after) in [
        ("molecular_resolution 8 4096", "molecular_resolution 0 4096"),
        (
            "molecular_resolution 8 4096",
            "molecular_resolution 65537 4096",
        ),
        (
            "molecular_resolution 8 4096",
            "molecular_resolution 8 131073",
        ),
        ("molecular_resolution 8 4096", "molecular_resolution 8 0"),
        (
            "molecular_resolution 8 4096",
            "molecular_resolution 8 4096 extra",
        ),
        (
            "molecular_wall diffusive 16 0.0001 30000 0.05",
            "molecular_wall unknown",
        ),
        (
            "molecular_wall diffusive 16 0.0001 30000 0.05",
            "molecular_wall diffusive 0 0.0001 30000 0.05",
        ),
        (
            "molecular_wall diffusive 16 0.0001 30000 0.05",
            "molecular_wall diffusive 65 0.0001 30000 0.05",
        ),
        (
            "molecular_wall diffusive 16 0.0001 30000 0.05",
            "molecular_wall none extra",
        ),
        (
            "molecular declared 2",
            "molecular declared 99999999999999999999",
        ),
    ] {
        let input = INPUT.replace(before, after);
        assert!(
            ReedPerformance::from_bytes(input.as_bytes(), 37).is_err(),
            "accepted {after}"
        );
    }
    // A source band that excludes the positive near-DC feedback probe must
    // refuse; no low-frequency equilibrium replacement is allowed.
    let band = INPUT.replace("diffusive 16 0.0001", "diffusive 16 1");
    assert!(matches!(
        ReedPerformance::from_bytes(band.as_bytes(), 37),
        Err(ReedPerformanceError::Model(_))
    ));
    let existing = include_str!("../../../../../examples/reed-duct.performance");
    assert!(
        ReedPerformance::from_bytes(existing.as_bytes(), 37)
            .unwrap()
            .molecular_driving_point()
            .is_none()
    );
}

