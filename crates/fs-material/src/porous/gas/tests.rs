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
        let original = card
            .surface_impedance_rigid(PorousAmbient::from_gas(&gas).unwrap(), omega, 0.05)
            .unwrap();
        assert_eq!(
            law.surface_impedance_rigid(&gas, omega, 0.05).unwrap(),
            original
        );
    }
}

#[test]
fn gas_substitution_refreshes_drag_and_all_ambient_coefficients() {
    let card = PorousCard::melamine_uf_nguyen2024();
    let law = PorousGasLaw::try_new(card.clone(), 1.81e-5).unwrap();
    let composition =
        Composition::new(CompositionBasis::MoleFraction, vec![0.6, 0.2, 0.2]).unwrap();
    let mixture = resolve_gas_mixture(
        &[
            GasSpecies::Nitrogen,
            GasSpecies::Oxygen,
            GasSpecies::CarbonDioxide,
        ],
        &composition,
        293.15,
        101_325.0,
    )
    .unwrap();
    let states = [
        air(293.15, 101_325.0),
        air(323.15, 101_325.0),
        air(293.15, 200_000.0),
        mixture,
    ];
    let omega = core::f64::consts::TAU * 600.0;
    let baseline = law
        .surface_impedance_rigid(&states[0], omega, 0.05)
        .unwrap();
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
        let expected = resolved
            .surface_impedance_rigid(ambient, omega, 0.05)
            .unwrap();
        let actual = law.surface_impedance_rigid(gas, omega, 0.05).unwrap();
        assert!((actual - expected).abs() / expected.abs() < 1e-14);
        if index != 0 {
            assert!((actual - baseline).abs() / baseline.abs() > 1e-3);
        }
    }
    // A hot gas cannot retain the characterization's drag coefficient.
    let hot = states[1];
    let frozen = card
        .surface_impedance_rigid(PorousAmbient::from_gas(&hot).unwrap(), omega, 0.05)
        .unwrap();
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
        assert!(
            law.surface_impedance_rigid(&gas, core::f64::consts::TAU * frequency, 0.05)
                .is_err()
        );
    }
    for thickness in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(
            law.surface_impedance_rigid(&gas, core::f64::consts::TAU * 600.0, thickness)
                .is_err()
        );
    }
    let mut stale = gas;
    stale.prandtl *= 2.0;
    assert!(
        law.surface_impedance_rigid(&stale, core::f64::consts::TAU * 600.0, 0.05)
            .is_err()
    );
}

fn liquid(viscosity: f64, conductivity: f64, ratio: f64) -> ThermoviscousFluidState {
    // Manufactured single-phase liquid, not a measured water calibration.
    ThermoviscousFluidState::from_heat_capacity_ratio(
        crate::fluid_acoustics::FluidAcousticState::try_new(1000.0, 1000.0).unwrap(),
        crate::fluid_transport::FluidTransportState::try_new(
            1000.0,
            4200.0,
            conductivity,
            viscosity,
        )
        .unwrap(),
        ratio,
    )
    .unwrap()
}

fn numerical_skeleton(lafarge: bool) -> PorousFluidLaw {
    // Deliberately broad frequency support for asymptotic numerical references;
    // it does not extend the measured foam cards' source bands.
    PorousFluidLaw::try_new(
        PorousCard {
            name: "manufactured-rigid-pores",
            porosity: 0.45,
            flow_resistivity: 1.0e-3 / 2.0e-10,
            tortuosity: 1.4,
            viscous_length_m: 100.0e-6,
            thermal_length_m: 200.0e-6,
            thermal_permeability_m2: lafarge.then_some(3.0e-10),
            provenance: "manufactured single-fluid rigid-frame limit fixture; no measurement claim",
            valid_hz: (1.0e-9, 1.0e16),
        },
        1.0e-3,
    )
    .unwrap()
}

#[test]
fn gas_and_general_fluid_views_share_the_original_constitutive_kernel() {
    let law = PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap();
    for gas in [air(293.15, 101_325.0), air(323.15, 200_000.0)] {
        for frequency in [100.0, 600.0, 1700.0, 4000.0] {
            let omega = core::f64::consts::TAU * frequency;
            let reference = law.surface_impedance_rigid(&gas, omega, 0.05).unwrap();
            let response = law.effective_fluid(&gas, omega).unwrap();
            let generic = response.surface_impedance_rigid(0.05).unwrap();
            // The legacy gas route retains gamma*p; the generic route retains
            // rho*c*c. The gas EOS makes them equal up to rounding.
            assert!((generic - reference).abs() / reference.abs() < 2.0e-12);
            assert_eq!(*response.fluid(), gas.thermoviscous_state().unwrap());
            assert_eq!(response.angular_frequency_rad_s(), omega);
            assert_eq!(
                generic,
                law.surface_impedance_rigid_fluid(&gas, omega, 0.05)
                    .unwrap()
            );
        }
    }
}

#[test]
fn liquid_low_and_high_frequency_limits_keep_its_actual_bulk_modulus() {
    let fluid = liquid(1.0e-3, 0.6, 1.2);
    for lafarge in [false, true] {
        let law = numerical_skeleton(lafarge);
        let card = law.reference_card();
        let low = law.effective_fluid(&fluid, 1.0e-6).unwrap();
        let high = law.effective_fluid(&fluid, 1.0e14).unwrap();
        let bulk = fluid.acoustic().adiabatic_bulk_modulus_pa();
        let isothermal = bulk / (fluid.heat_capacity_ratio() * card.porosity);
        let adiabatic = bulk / card.porosity;
        // G1: isothermal and adiabatic thermodynamic limits, and the Darcy
        // low-frequency loss versus the tortuosity-controlled high limit.
        assert!((low.dynamic_bulk_modulus_pa().re / isothermal - 1.0).abs() < 1.0e-12);
        assert!(low.dynamic_bulk_modulus_pa().im / isothermal < 1.0e-8);
        assert!((high.dynamic_bulk_modulus_pa().re / adiabatic - 1.0).abs() < 1.0e-6);
        assert!(
            ((-low.dynamic_density_kg_m3().im * 1.0e-6) / card.flow_resistivity - 1.0).abs()
                < 1.0e-12
        );
        let inertial_density = card.tortuosity * fluid.acoustic().density_kg_m3() / card.porosity;
        assert!((high.dynamic_density_kg_m3().re / inertial_density - 1.0).abs() < 2.0e-6);
        // A thin rigid-backed layer tends to the liquid's isothermal
        // compliance. There is no atmospheric pressure in this reference.
        let depth = 1.0e-5;
        let surface = low.surface_impedance_rigid(depth).unwrap();
        assert!((-surface.im * 1.0e-6 * depth / isothermal - 1.0).abs() < 1.0e-10);
        assert!(surface.re >= 0.0);
    }
}

#[test]
fn cp_equals_cv_removes_only_thermal_relaxation_exactly() {
    let law = numerical_skeleton(true);
    let first = liquid(1.0e-3, 0.1, 1.0);
    let second = liquid(1.0e-3, 10.0, 1.0);
    let more_viscous = liquid(4.0e-3, 0.1, 1.0);
    for omega in [1.0e-6, 10.0, 1.0e5, 1.0e14] {
        let a = law.effective_fluid(&first, omega).unwrap();
        let b = law.effective_fluid(&second, omega).unwrap();
        let c = law.effective_fluid(&more_viscous, omega).unwrap();
        let exact_bulk = C64::new(
            first.acoustic().adiabatic_bulk_modulus_pa() / law.reference_card().porosity,
            0.0,
        );
        assert_eq!(a.dynamic_bulk_modulus_pa(), exact_bulk);
        assert_eq!(a.dynamic_bulk_modulus_pa(), b.dynamic_bulk_modulus_pa());
        assert_eq!(a.dynamic_bulk_modulus_pa(), c.dynamic_bulk_modulus_pa());
        assert_eq!(a.characteristic(), b.characteristic());
        assert_eq!(
            a.surface_impedance_rigid(0.01).unwrap(),
            b.surface_impedance_rigid(0.01).unwrap()
        );
        assert_ne!(a.dynamic_density_kg_m3(), c.dynamic_density_kg_m3());
        assert!(a.dynamic_density_kg_m3().im < 0.0);
    }
}

#[test]
fn fluid_substitution_refreshes_drag_and_thermal_response_without_mutating_pores() {
    let law = numerical_skeleton(true);
    let original = law.reference_card().clone();
    let first = liquid(1.0e-3, 0.1, 1.2);
    let second = liquid(2.0e-3, 0.8, 1.2);
    let omega = 1000.0;
    let a = law.effective_fluid(&first, omega).unwrap();
    let b = law.effective_fluid(&second, omega).unwrap();
    assert_eq!(
        b.flow_resistivity_pa_s_m2(),
        second.transport().dynamic_viscosity_pa_s() / law.viscous_permeability_m2()
    );
    assert!((b.flow_resistivity_pa_s_m2() / a.flow_resistivity_pa_s_m2() - 2.0).abs() < 1.0e-15);
    assert!(
        (a.dynamic_density_kg_m3() - b.dynamic_density_kg_m3()).abs()
            / a.dynamic_density_kg_m3().abs()
            > 0.1
    );
    assert!(
        (a.dynamic_bulk_modulus_pa() - b.dynamic_bulk_modulus_pa()).abs()
            / a.dynamic_bulk_modulus_pa().abs()
            > 1.0e-3
    );
    assert_ne!(
        a.surface_impedance_rigid(0.01).unwrap(),
        b.surface_impedance_rigid(0.01).unwrap()
    );
    assert_eq!(law.reference_card(), &original);
    assert_eq!(law.effective_fluid(&first, omega).unwrap(), a);
}

#[test]
fn generic_fluid_refuses_missing_source_transport_and_original_support_limits() {
    struct MissingCv;
    impl crate::fluid_acoustics::FluidAcoustics for MissingCv {
        fn acoustic_state(
            &self,
        ) -> Result<crate::fluid_acoustics::FluidAcousticState, MaterialError> {
            panic!("the porous law must query the complete original fluid owner once")
        }
    }
    impl crate::fluid_transport::FluidTransport for MissingCv {
        fn transport_state(
            &self,
        ) -> Result<crate::fluid_transport::FluidTransportState, MaterialError> {
            panic!("the porous law cannot construct a separate incomplete transport state")
        }
    }
    impl ThermoviscousFluid for MissingCv {
        fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
            Err(MaterialError::State {
                what: "original liquid has no Cv claim".into(),
            })
        }
    }
    let law = PorousFluidLaw::try_new(PorousCard::melamine_uf_nguyen2024(), 1.81e-5).unwrap();
    let omega = core::f64::consts::TAU * 600.0;
    assert_eq!(
        law.effective_fluid(&MissingCv, omega).unwrap_err(),
        MaterialError::State {
            what: "original liquid has no Cv claim".into(),
        }
    );
    let valid = liquid(1.0e-3, 0.6, 1.2);
    for frequency in [0.0, 99.0, 4001.0, f64::NAN, f64::INFINITY] {
        assert!(
            law.effective_fluid(&valid, core::f64::consts::TAU * frequency)
                .is_err()
        );
    }
    let sample = law.effective_fluid(&valid, omega).unwrap();
    for thickness in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(sample.surface_impedance_rigid(thickness).is_err());
    }
    let mut stale = air(293.15, 101_325.0);
    stale.dynamic_viscosity *= 2.0;
    assert!(law.effective_fluid(&stale, omega).is_err());
    assert_eq!(law.effective_fluid(&valid, omega).unwrap(), sample);
}
