//! Native input and source feedback for the original spatial material law.
use super::*;
use fs_material::gas_species::input::{ThermochemicalPhase, ThermochemicalSources};

#[path = "../../../../../render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const EXAMPLE: &str =
    include_str!("../../../../../../examples/flue-dg-field-molecular.performance");
const FIELD_LAW: &str =
    "dg_molecular declared 2 N2 0 CO2 1 vibration 1 1000 1 power 10 293.15 0 273.15 323.15\n";

fn internal_source() -> String {
    input()
        .replace(
            "gas mixture 293.15 101325 mole 2 N2 0.8 O2 0.2",
            "gas thermochemical 308.15 101325 mole 2 N2 0.96 CO2 0.04",
        )
        .replace(
            "dg_tube 0.05 0.001 open 1 2 8 50 2000 0.01 1e-11 256",
            concat!(
                "dg_bore 1 open 2 2 12 50 2000 0.01 1e-11 512 4 gas_field 3 mole 2 N2 CO2\n",
                "dg_segment cone 0.001 0.002 0.05\n",
                "gas_station 0 308.15 0.96 0.04\n",
                "gas_station 0.025 300.15 0.98 0.02\n",
                "gas_station 0.05 293.15 0.99 0.01"
            ),
        )
        .replace("\nflue ", &format!("\n{FIELD_LAW}flue "))
}

#[test]
fn spatial_molecular_flue_keeps_local_gas_feedback_and_block_identity() {
    // Manufactured NASA-9/transport packs and reservoir coefficients test the
    // exact data path; these are not measured molecular relaxation constants.
    let a = gas_support::packs("N2", 0.028, 3.5, 200.0);
    let b = gas_support::packs("CO2", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let gate = CancelGate::new_clock_free();
    let build = |text: &str, block| {
        FlueDgFilePerformance::from_bytes_with_gas_sources(
            text.as_bytes(),
            block,
            &gate,
            Some(&sources),
        )
        .unwrap()
    };
    let text = internal_source();
    let source = build(&text, 512);
    let duct = source.renderer().voice().dg().unwrap().load().duct();
    assert!(duct.molecular_relaxation_field().is_some());
    assert!(duct.molecular_relaxation().is_none());
    assert_eq!(duct.total_dofs(), 300);
    assert_eq!(duct.gas_field().unwrap().stations_m(), &[0.0, 0.025, 0.05]);
    assert_ne!(duct.inlet_gas(), duct.outlet_gas());
    let mut whole = source.into_renderer();
    let mut partitioned = build(&text, 37).into_renderer();
    let (mut reference, mut other) = (vec![0.0; 512], vec![0.0; 512]);
    whole.block(&mut reference).unwrap();
    for block in other.chunks_mut(37) {
        partitioned.block(block).unwrap();
    }
    assert_eq!(reference, other);
    assert_eq!(
        whole.voice().dg().unwrap().load().state(),
        partitioned.voice().dg().unwrap().load().state()
    );
    for changed in [
        text.replace(FIELD_LAW, ""),
        text.replace("power 10 ", "power 1 "),
        text.replace(
            "gas_station 0.05 293.15 0.99 0.01",
            "gas_station 0.05 313.15 0.99 0.01",
        ),
        text.replace(
            "gas_station 0.05 293.15 0.99 0.01",
            "gas_station 0.05 293.15 0.9 0.1",
        ),
    ] {
        let mut changed = build(&changed, 512).into_renderer();
        assert_eq!(
            whole.voice().dg().unwrap().load().duct().inlet_gas(),
            changed.voice().dg().unwrap().load().duct().inlet_gas()
        );
        changed.block(&mut other).unwrap();
        assert!(
            reference
                .iter()
                .zip(&other)
                .any(|(a, b)| (a - b).abs() > 1e-12)
        );
    }
    // Exercise the bundled exterior-pressure source with the exact same source
    // registry. Its atmosphere remains independent of the internal field law.
    let external = build(EXAMPLE, 37);
    assert!(
        external
            .renderer()
            .voice()
            .dg()
            .unwrap()
            .load()
            .duct()
            .molecular_relaxation_field()
            .is_some()
    );
    assert_eq!(
        external.renderer().observer_gas().unwrap().law,
        crate::air_path::GasAbsorptionModel::StokesKirchhoff
    );
}

#[test]
fn spatial_molecular_native_admission_checks_remote_sources_and_complete_state_budget() {
    let a = gas_support::packs("N2", 0.028, 3.5, 200.0);
    let b = gas_support::packs("CO2", 0.044, 4.5, 200.0);
    let bindings = [a.binding(), b.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let gate = CancelGate::new_clock_free();
    let text = internal_source();
    for refused in [
        text.replace("1e-11 512 4", "1e-11 299 4"),
        text.replace("gas_station 0.05 293.15", "gas_station 0.05 333.15"),
        text.replace(FIELD_LAW.trim(), "dg_molecular declared 1 N2 0"),
        text.replace(FIELD_LAW.trim(), "dg_molecular iso9613-terrestrial-air"),
    ] {
        assert!(
            FlueDgFilePerformance::from_bytes_with_gas_sources(
                refused.as_bytes(),
                37,
                &gate,
                Some(&sources),
            )
            .is_err(),
            "{refused}"
        );
    }
    // Exact capacity is admitted; no nodal-count estimate may replace the four
    // original thermal quadrature points. A missing registry stays an error.
    assert!(
        FlueDgFilePerformance::from_bytes_with_gas_sources(
            text.replace("1e-11 512 4", "1e-11 300 4").as_bytes(),
            37,
            &gate,
            Some(&sources),
        )
        .is_ok()
    );
    assert!(FlueDgFilePerformance::from_bytes(text.as_bytes(), 37, &gate).is_err());
}
