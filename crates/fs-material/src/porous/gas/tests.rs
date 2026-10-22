use super::*;
use crate::gas::GasSpec;
use crate::gas_species::{GasSpecies, resolve_gas_mixture};
use fs_qty::semantic::{Composition, CompositionBasis};

fn air(temperature: f64, pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, pressure).unwrap()
}

#[test]
fn characterization_state_keeps_original_response_and_source() {
    let gas = air(293.15, 101_325.0);
    let mut card = PorousCard::melamine_uf_nguyen2024();
    // Analytic fixture: express the same specified intrinsic permeability in
    // this exactly evaluated gas. This is not a new measured calibration.
    card.flow_resistivity = gas.dynamic_viscosity / 31.47e-10;
    card.provenance = "analytic gas-rescaling fixture with explicitly fixed pore geometry";
    let law = PorousGasLaw::try_new(card.clone(), gas.dynamic_viscosity).unwrap();
    assert_eq!(law.reference_card().provenance, card.provenance);
    assert_eq!(law.reference_card().flow_resistivity, card.flow_resistivity);
    assert!((law.viscous_permeability_m2() / 31.47e-10 - 1.0).abs() < 1e-15);
    for frequency in [100.0, 600.0, 1700.0, 4000.0] {
        let omega = core::f64::consts::TAU * frequency;
        let original = card.surface_impedance_rigid(
            PorousAmbient::from_gas(&gas).unwrap(), omega, 0.05,
        ).unwrap();
        assert_eq!(law.surface_impedance_rigid(&gas, omega, 0.05).unwrap(), original);
    }
}

#[test]
fn gas_substitution_refreshes_drag_and_all_ambient_coefficients() {
    let card = PorousCard::melamine_uf_nguyen2024();
    let law = PorousGasLaw::try_new(card.clone(), 1.81e-5).unwrap();
    let composition = Composition::new(
        CompositionBasis::MoleFraction, vec![0.6, 0.2, 0.2],
    ).unwrap();
    let mixture = resolve_gas_mixture(
        &[GasSpecies::Nitrogen, GasSpecies::Oxygen, GasSpecies::CarbonDioxide],
        &composition, 293.15, 101_325.0,
    ).unwrap();
    let states = [
        air(293.15, 101_325.0), air(323.15, 101_325.0),
        air(293.15, 200_000.0), mixture,
    ];
    let omega = core::f64::consts::TAU * 600.0;
    let baseline = law.surface_impedance_rigid(&states[0], omega, 0.05).unwrap();
    for (index, gas) in states.iter().enumerate() {
        let ambient = PorousAmbient::from_gas(gas).unwrap();
        assert_eq!(ambient.rho0, gas.density);
        assert_eq!(ambient.c0, gas.sound_speed);
        assert_eq!(ambient.eta, gas.dynamic_viscosity);
        assert_eq!(ambient.gamma, gas.gamma);
        assert_eq!(ambient.p0, gas.pressure);
        assert_eq!(ambient.prandtl, gas.prandtl);
        let mut resolved = card.clone();
        resolved.flow_resistivity = gas.dynamic_viscosity / law.viscous_permeability_m2();
        let expected = resolved.surface_impedance_rigid(ambient, omega, 0.05).unwrap();
        let actual = law.surface_impedance_rigid(gas, omega, 0.05).unwrap();
        assert!((actual - expected).abs() / expected.abs() < 1e-14);
        if index != 0 {
            assert!((actual - baseline).abs() / baseline.abs() > 1e-3);
        }
    }
    // A hot gas cannot retain the characterization's drag coefficient.
    let hot = states[1];
    let frozen = card.surface_impedance_rigid(
        PorousAmbient::from_gas(&hot).unwrap(), omega, 0.05,
    ).unwrap();
    let actual = law.surface_impedance_rigid(&hot, omega, 0.05).unwrap();
    assert!((actual - frozen).abs() / actual.abs() > 1e-3);
    assert_eq!(law.reference_card().flow_resistivity, card.flow_resistivity);
}

#[test]
fn missing_characterization_invalid_gas_and_out_of_support_queries_refuse() {
    let card = PorousCard::melamine_uf_nguyen2024();
    for viscosity in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(PorousGasLaw::try_new(card.clone(), viscosity).is_err());
    }
    let mut anonymous = card.clone();
    anonymous.provenance = " ";
    assert!(PorousGasLaw::try_new(anonymous, 1.81e-5).is_err());
    let law = PorousGasLaw::try_new(card, 1.81e-5).unwrap();
    let gas = air(293.15, 101_325.0);
    for frequency in [0.0, 99.0, 4001.0, f64::NAN, f64::INFINITY] {
        assert!(law.surface_impedance_rigid(&gas, core::f64::consts::TAU * frequency, 0.05).is_err());
    }
    for thickness in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(law.surface_impedance_rigid(&gas, core::f64::consts::TAU * 600.0, thickness).is_err());
    }
    let mut stale = gas;
    stale.prandtl *= 2.0;
    assert!(law.surface_impedance_rigid(&stale, core::f64::consts::TAU * 600.0, 0.05).is_err());
}
