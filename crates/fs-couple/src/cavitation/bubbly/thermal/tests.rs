use super::*;
use fs_flux::cavitation::inception::EquilibriumControls;
use fs_flux::cavitation::{
    BubbleGas, BubbleLimits, BubblePopulation, BubblePopulationClass, KellerMiksisLimits,
};
use fs_material::gas::GasSpec;

fn liquid() -> CavitationState {
    CavitationState::try_new(300.0, 1000.0, 0.001, 0.072, 2000.0, 1500.0).unwrap()
}
fn acoustic_limits() -> BubblyLiquidLimits {
    BubblyLiquidLimits {
        maximum_void_fraction: 0.001,
        maximum_wavenumber_radius: 0.1,
        maximum_wavenumber_spacing: 1.0,
        maximum_fractional_radius: 0.01,
        maximum_neighbor_scattered_pressure_ratio: 0.1,
        maximum_multiple_scattering_parameter: 0.1,
    }
}
fn thermal_limits() -> ThermalBubbleLimits {
    ThermalBubbleLimits {
        maximum_gas_wavenumber_radius: 0.1,
        maximum_gas_liquid_density_ratio: 0.01,
        maximum_vapor_pressure_fraction: 0.1,
    }
}
fn class(radius: f64, fraction: f64) -> BubblePopulationClass {
    let bubble = KellerMiksisBubble::try_new(
        liquid(),
        BubbleGas::equilibrium(liquid(), radius, 100_000.0, 1.0).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.25 * radius,
                maximum_radius_m: 3.0 * radius,
                maximum_wall_mach_number: 0.01,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap();
    BubblePopulationClass::try_new(
        bubble,
        fraction / ((4.0 * core::f64::consts::PI / 3.0) * radius.powi(3)),
    )
    .unwrap()
}
fn gas_for(class: &BubblePopulationClass) -> Option<GasState> {
    class.linearization().map(|linear| {
        let pressure = class
            .bubble()
            .gas()
            .pressure_at_radius(linear.equilibrium_radius_m)
            .unwrap();
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, pressure).unwrap()
    })
}
fn gases(population: &BubblePopulation) -> Vec<Option<GasState>> {
    population.classes().iter().map(gas_for).collect()
}
fn distribution() -> BubblePopulation {
    BubblePopulation::try_new(liquid(), vec![class(5e-5, 1e-5), class(2e-4, 2e-6)], 2).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-30),
        "{a:e} != {b:e}"
    );
}
fn near_complex(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-30),
        "{a:?} != {b:?}"
    );
}

// Independent spherical heat-diffusion eigenmodes, not the native coth kernel.
// Absolute F-tail <= 2*x/(pi^4*N^3), below 1e-11 over these physical fixtures.
fn modal_exponent(gamma: f64, x: f64) -> C64 {
    let mut correction = C64::ZERO;
    for n in 1..=8192 {
        let a = (core::f64::consts::PI * f64::from(n)).powi(2);
        correction = correction + C64::new(0.0, 6.0 * x) / C64::new(a * a, -a * x);
    }
    C64::from_re(gamma) / (C64::ONE + (C64::ONE + correction).scale(gamma - 1.0))
}

#[test]
fn thermal_population_matches_independent_diffusion_dispersion_and_pressure_volume_work() {
    let material = liquid();
    let original =
        PolydisperseBubblyLiquid::new(&material, distribution(), acoustic_limits()).unwrap();
    let thermal = ThermalPolydisperseBubblyLiquid::new(
        &original,
        &gases(original.population()),
        thermal_limits(),
        &mut || false,
    )
    .unwrap();
    assert!(core::ptr::eq(thermal.original(), &original));
    let local_pressure = 0.02;
    for omega in [1000.0_f64, 5000.0, 15000.0] {
        let response = thermal.response(omega, 0.1).unwrap();
        let mut k_squared = C64::from_re((omega / 1500.0).powi(2));
        let mut thermal_power = 0.0;
        let mut viscous_power = 0.0;
        for member in original.population().classes() {
            let r = member.linearization().unwrap().equilibrium_radius_m;
            let gas = gas_for(member).unwrap();
            let diffusivity = gas.thermal_conductivity / gas.density / gas.specific_heat_cp;
            let kappa = modal_exponent(gas.gamma, omega * r * r / diffusivity);
            let stiffness = (kappa.scale(3.0 * gas.pressure) - C64::from_re(2.0 * 0.072 / r))
                .scale(1.0 / (1000.0 * r * r));
            let tau = r / 1500.0;
            let viscous = 4.0 * 0.001 / (1000.0 * r * r);
            let drive = C64::new(1.0, -omega * tau);
            let denominator = drive * stiffness
                - C64::new((1.0 + viscous * tau) * omega * omega, viscous * omega);
            let chi = -(drive / denominator).scale(1.0 / (1000.0 * r * r));
            let n = member.number_density_per_m3();
            k_squared = k_squared
                + chi.scale(-4.0 * core::f64::consts::PI * n * 1000.0 * r.powi(3) * omega * omega);

            // Direct harmonic p*dV power, independent of the loss-coefficient
            // formula. Positive compression work heats the gas each cycle.
            let displacement = chi.scale(local_pressure);
            let gas_pressure = kappa * displacement.scale(-3.0 * gas.pressure);
            let volume_velocity =
                C64::new(0.0, -omega * 4.0 * core::f64::consts::PI * r.powi(3)) * displacement;
            thermal_power += -0.5 * n * (gas_pressure * volume_velocity.conj()).re;
            // Integrate the native instantaneous Newtonian dissipation over
            // one period, without using the implementation's mean coefficient.
            let mut cycle = 0.0;
            for sample in 0..128 {
                let phase = core::f64::consts::TAU * f64::from(sample) / 128.0;
                let velocity =
                    omega * r * (displacement.re * phase.sin() + displacement.im * phase.cos());
                cycle += 16.0 * core::f64::consts::PI * 0.001 * r * velocity * velocity / 128.0;
            }
            viscous_power += n * cycle;
        }
        near_complex(
            response.acoustic().wavenumber_per_m(),
            k_squared.sqrt(),
            2e-10,
        );
        let heat = response.mean_absorption_w_m3(local_pressure).unwrap();
        near(heat.gas_thermal_w_m3, thermal_power, 2e-9);
        near(heat.liquid_viscous_w_m3, viscous_power, 2e-10);
        assert!(heat.gas_thermal_w_m3 > 0.0 && heat.liquid_viscous_w_m3 > 0.0);
        let wave = response.acoustic().segment_wave(0.005).unwrap();
        near_complex(
            wave.specific_impedance,
            C64::from_re(1000.0 * omega) / k_squared.sqrt(),
            2e-10,
        );
        assert_ne!(
            response.acoustic().wavenumber_per_m(),
            original.response(omega, 0.1).unwrap().wavenumber_per_m()
        );
        assert!(response.mean_absorption_w_m3(0.10001).is_err());
        assert!(response.mean_absorption_w_m3(f64::NAN).is_err());
        assert_eq!(
            response.mean_absorption_w_m3(0.0).unwrap(),
            BubblyAbsorption {
                gas_thermal_w_m3: 0.0,
                liquid_thermal_w_m3: 0.0,
                liquid_viscous_w_m3: 0.0
            }
        );
        let doubled = response.mean_absorption_w_m3(2.0 * local_pressure).unwrap();
        near(doubled.gas_thermal_w_m3, 4.0 * heat.gas_thermal_w_m3, 2e-15);
    }
}

#[test]
fn thermal_single_class_and_split_histogram_preserve_the_same_population_response() {
    let material = liquid();
    let member = class(1e-4, 1e-5);
    let mono =
        DiluteBubblyLiquid::from_population_class(&material, member, acoustic_limits()).unwrap();
    let one = PolydisperseBubblyLiquid::new(
        &material,
        BubblePopulation::try_new(material, vec![member], 1).unwrap(),
        acoustic_limits(),
    )
    .unwrap();
    let half =
        BubblePopulationClass::try_new(*member.bubble(), member.number_density_per_m3() / 2.0)
            .unwrap();
    let split = PolydisperseBubblyLiquid::new(
        &material,
        BubblePopulation::try_new(material, vec![half, half], 2).unwrap(),
        acoustic_limits(),
    )
    .unwrap();
    let a = ThermalDiluteBubblyLiquid::new(&mono, gas_for(&member), thermal_limits())
        .unwrap()
        .response(10000.0, 0.1)
        .unwrap();
    let b = ThermalPolydisperseBubblyLiquid::new(
        &one,
        &gases(one.population()),
        thermal_limits(),
        &mut || false,
    )
    .unwrap()
    .response(10000.0, 0.1)
    .unwrap();
    let c = ThermalPolydisperseBubblyLiquid::new(
        &split,
        &gases(split.population()),
        thermal_limits(),
        &mut || false,
    )
    .unwrap()
    .response(10000.0, 0.1)
    .unwrap();
    near_complex(
        a.acoustic().wavenumber_per_m(),
        b.acoustic().wavenumber_per_m(),
        5e-15,
    );
    near_complex(
        c.acoustic().wavenumber_per_m(),
        b.acoustic().wavenumber_per_m(),
        5e-15,
    );
    for value in [a.diagnostics(), c.diagnostics()] {
        near(
            value.thermal_loss_per_pressure_squared_w_m3_pa2,
            b.diagnostics().thermal_loss_per_pressure_squared_w_m3_pa2,
            5e-15,
        );
        near(
            value.viscous_loss_per_pressure_squared_w_m3_pa2,
            b.diagnostics().viscous_loss_per_pressure_squared_w_m3_pa2,
            5e-15,
        );
    }
    near(
        c.acoustic().diagnostics().multiple_scattering_parameter,
        b.acoustic().diagnostics().multiple_scattering_parameter,
        5e-15,
    );
}

#[test]
fn pressure_continuation_requires_current_gas_transport_and_retains_the_original_inventory() {
    let material = liquid();
    let original =
        PolydisperseBubblyLiquid::new(&material, distribution(), acoustic_limits()).unwrap();
    let old_gases = gases(original.population());
    let controls = EquilibriumControls {
        relative_radius_tolerance: 1e-13,
        absolute_radius_tolerance_m: 1e-19,
        absolute_pressure_tolerance_pa: 1e-8,
        maximum_iterations: 100,
    };
    let (continued, _) = original
        .at_equilibrium_pressure(40_000.0, &controls, &mut || false)
        .unwrap();
    assert!(matches!(
        ThermalPolydisperseBubblyLiquid::new(&continued, &old_gases, thermal_limits(), &mut || {
            false
        },),
        Err(BubblyLiquidError::ThermalClass { class_index: 0, .. })
    ));
    let fresh = gases(continued.population());
    let prepared =
        ThermalPolydisperseBubblyLiquid::new(&continued, &fresh, thermal_limits(), &mut || false)
            .unwrap();
    prepared.response(5000.0, 0.01).unwrap();
    for ((before, after), (old, new)) in original
        .population()
        .classes()
        .iter()
        .zip(continued.population().classes())
        .zip(old_gases.iter().zip(&fresh))
    {
        assert_eq!(before.bubble(), after.bubble());
        assert_eq!(
            before.number_density_per_m3(),
            after.number_density_per_m3()
        );
        let old = old.unwrap();
        let new = new.unwrap();
        assert!(new.density < old.density);
        assert_eq!(new.thermal_conductivity, old.thermal_conductivity);
        assert_eq!(new.specific_heat_cp, old.specific_heat_cp);
    }
    assert!(core::ptr::eq(
        prepared.original().source(),
        original.source()
    ));
}

#[test]
fn empty_thermal_population_is_exact_liquid_and_failed_queries_leave_accepted_owners_usable() {
    let material = liquid();
    let empty_class = class(1e-4, 0.0);
    let mono = DiluteBubblyLiquid::from_population_class(&material, empty_class, acoustic_limits())
        .unwrap();
    let empty = PolydisperseBubblyLiquid::new(
        &material,
        BubblePopulation::try_new(material, vec![empty_class], 1).unwrap(),
        acoustic_limits(),
    )
    .unwrap();
    let native = empty.response(1000.0, 1e5).unwrap();
    let a = ThermalDiluteBubblyLiquid::new(&mono, None, thermal_limits())
        .unwrap()
        .response(1000.0, 1e5)
        .unwrap();
    let b = ThermalPolydisperseBubblyLiquid::new(&empty, &[None], thermal_limits(), &mut || false)
        .unwrap()
        .response(1000.0, 1e5)
        .unwrap();
    assert_eq!(b.acoustic(), &native);
    assert_eq!(a.acoustic().wavenumber_per_m(), native.wavenumber_per_m());
    assert_eq!(
        a.mean_absorption_w_m3(1e5).unwrap(),
        b.mean_absorption_w_m3(1e5).unwrap()
    );
    assert_eq!(
        b.diagnostics().thermal_loss_per_pressure_squared_w_m3_pa2,
        0.0
    );

    let original =
        PolydisperseBubblyLiquid::new(&material, distribution(), acoustic_limits()).unwrap();
    let states = gases(original.population());
    let limits = ThermalBubbleLimits {
        maximum_gas_wavenumber_radius: 0.005,
        ..thermal_limits()
    };
    let mut preparation_polls = 0;
    let thermal = ThermalPolydisperseBubblyLiquid::new(&original, &states, limits, &mut || {
        preparation_polls += 1;
        false
    })
    .unwrap();
    let mut response_polls = 0;
    let accepted = thermal
        .response_controlled(1000.0, 0.01, &mut || {
            response_polls += 1;
            false
        })
        .unwrap();
    for stop in 1..=preparation_polls {
        let mut count = 0;
        assert!(
            ThermalPolydisperseBubblyLiquid::new(&original, &states, limits, &mut || {
                count += 1;
                count == stop
            })
            .is_err()
        );
    }
    for stop in 1..=response_polls {
        let mut count = 0;
        assert!(
            thermal
                .response_controlled(1000.0, 0.01, &mut || {
                    count += 1;
                    count == stop
                })
                .is_err()
        );
    }
    assert!(thermal.response(20000.0, 0.01).is_err()); // gas, not liquid, acoustic compactness
    assert!(thermal.response(1000.0, 10000.0).is_err()); // original linear-amplitude bound
    assert_eq!(thermal.response(1000.0, 0.01).unwrap(), accepted);
    assert!(
        ThermalPolydisperseBubblyLiquid::new(&original, &states[..1], limits, &mut || false)
            .is_err()
    );
}

#[derive(Clone, Copy)]
struct ConductingTestLiquid {
    material: CavitationState,
    thermal: ThermalFluidTransportState,
}

impl CavitationProperties for ConductingTestLiquid {
    fn cavitation_state(&self) -> Result<CavitationState, MaterialError> {
        Ok(self.material)
    }
}

impl fs_material::fluid_transport::FluidTransport for ConductingTestLiquid {
    fn transport_state(
        &self,
    ) -> Result<fs_material::fluid_transport::FluidTransportState, MaterialError> {
        Ok(self.thermal.transport())
    }
}

impl ThermalFluidTransport for ConductingTestLiquid {
    fn thermal_transport_state(&self) -> Result<ThermalFluidTransportState, MaterialError> {
        Ok(self.thermal)
    }
}

fn conducting_material(conductivity: f64) -> ConductingTestLiquid {
    ConductingTestLiquid {
        material: liquid(),
        thermal: ThermalFluidTransportState::try_new(
            fs_material::fluid_transport::FluidTransportState::try_new(
                1000.0,
                4200.0,
                conductivity,
                0.001,
            )
            .unwrap(),
            300.0,
            None,
        )
        .unwrap(),
    }
}

fn liquid_limits() -> ConductingLiquidPopulationLimits {
    ConductingLiquidPopulationLimits {
        maximum_penetration_to_spacing: 0.1,
        maximum_fractional_interface_temperature: 0.01,
    }
}

#[test]
fn conducting_liquid_population_preserves_class_weights_and_closes_both_heat_channels() {
    let material = conducting_material(0.02);
    let member = class(1e-4, 1e-5);
    let half =
        BubblePopulationClass::try_new(*member.bubble(), 0.5 * member.number_density_per_m3())
            .unwrap();
    let mono =
        DiluteBubblyLiquid::from_population_class(&material, member, acoustic_limits()).unwrap();
    let split = PolydisperseBubblyLiquid::new(
        &material,
        BubblePopulation::try_new(liquid(), vec![half, half], 2).unwrap(),
        acoustic_limits(),
    )
    .unwrap();
    let monothermal = ThermalDiluteBubblyLiquid::new_with_conducting_liquid(
        &mono,
        gas_for(&member),
        thermal_limits(),
        liquid_limits(),
    )
    .unwrap();
    let population = ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid(
        &split,
        &gases(split.population()),
        thermal_limits(),
        liquid_limits(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(monothermal.liquid_thermal_state(), Some(material.thermal));
    assert!(core::ptr::eq(population.original().source(), &material));
    assert_eq!(population.thermal_bubbles().len(), 2);
    let omega = 10_000.0;
    let amplitude = 0.1;
    let single = monothermal.response(omega, amplitude).unwrap();
    let many = population.response(omega, amplitude).unwrap();
    near_complex(
        single.acoustic().wavenumber_per_m(),
        many.acoustic().wavenumber_per_m(),
        1e-14,
    );
    let heat = many.mean_absorption_w_m3(amplitude).unwrap();
    assert!(heat.gas_thermal_w_m3 > 0.0);
    assert!(heat.liquid_thermal_w_m3 > 0.0);
    assert!(heat.liquid_viscous_w_m3 > 0.0);
    near(
        heat.gas_thermal_w_m3 + heat.liquid_thermal_w_m3,
        many.diagnostics()
            .thermal_loss_per_pressure_squared_w_m3_pa2
            * amplitude
            * amplitude,
        5e-13,
    );
    let single_heat = single.mean_absorption_w_m3(amplitude).unwrap();
    near(single_heat.gas_thermal_w_m3, heat.gas_thermal_w_m3, 1e-14);
    near(
        single_heat.liquid_thermal_w_m3,
        heat.liquid_thermal_w_m3,
        1e-14,
    );
    let half_heat = many.mean_absorption_w_m3(amplitude * 0.5).unwrap();
    near(
        half_heat.liquid_thermal_w_m3,
        0.25 * heat.liquid_thermal_w_m3,
        1e-14,
    );
    assert!(many.diagnostics().maximum_fractional_interface_temperature > 0.0);

    // Independent gas eigenmode sum, then interface heat-flux continuity and
    // the original mass balance. No call to the native thermal exponent.
    let gas = gas_for(&member).unwrap();
    let r = member.linearization().unwrap().equilibrium_radius_m;
    let x = omega * r * r * gas.density * gas.specific_heat_cp / gas.thermal_conductivity;
    let fixed_wall = modal_exponent(gas.gamma, x);
    let f = (C64::from_re(gas.gamma) / fixed_wall - C64::ONE).scale(1.0 / (gas.gamma - 1.0));
    let yg = C64::new(0.0, -x / 3.0) * f.scale(gas.thermal_conductivity);
    let a = (0.5 * omega * r * r / material.thermal.transport().thermal_diffusivity_m2_s()).sqrt();
    let yl = C64::new(1.0 + a, -a).scale(material.thermal.transport().thermal_conductivity_w_m_k());
    let kappa = C64::from_re(gas.gamma) / (C64::ONE + (f * yl / (yg + yl)).scale(gas.gamma - 1.0));
    let stiffness =
        kappa.scale(3.0 * gas.pressure) - C64::new(2.0 * 0.072 / r, 4.0 * 0.001 * omega);
    let inertia = C64::from_re(1000.0 * r * r * omega * omega) / C64::new(1.0, -omega * r / 1500.0);
    let chi = -(stiffness - inertia).recip();
    let expected_squared = C64::from_re((omega / 1500.0).powi(2))
        - chi.scale(
            4.0 * core::f64::consts::PI
                * member.number_density_per_m3()
                * 1000.0
                * r.powi(3)
                * omega
                * omega,
        );
    let wave = many.acoustic().wavenumber_per_m();
    near_complex(wave * wave, expected_squared, 3e-9);
    let imposed = ThermalDiluteBubblyLiquid::new(&mono, Some(gas), thermal_limits())
        .unwrap()
        .response(omega, amplitude)
        .unwrap();
    assert!(
        (single.acoustic().wavenumber_per_m() - imposed.acoustic().wavenumber_per_m()).abs() > 1e-8
    );
}

#[test]
fn conducting_liquid_uses_original_if97_water_transport() {
    use fs_material::cavitation::SaturatedWaterCavitation;
    use fs_thermochem::water::if97::If97SaturationPair;
    let water = SaturatedWaterCavitation::try_new(
        If97SaturationPair::at_temperature(fs_qty::Temperature::new(293.15)).unwrap(),
        0.072,
    )
    .unwrap();
    let state = water.cavitation_state().unwrap();
    let template = class(1e-4, 1e-5);
    let bubble = KellerMiksisBubble::try_new(
        state,
        BubbleGas::equilibrium(state, 1e-4, 101_325.0, 1.0).unwrap(),
        template.bubble().compressibility_limits(),
    )
    .unwrap();
    let source = DiluteBubblyLiquid::new(&water, bubble, 1e-5, acoustic_limits()).unwrap();
    let gas = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        293.15,
        bubble.gas().pressure_at_radius(1e-4).unwrap(),
    )
    .unwrap();
    let thermal = ThermalDiluteBubblyLiquid::new_with_conducting_liquid(
        &source,
        Some(gas),
        thermal_limits(),
        liquid_limits(),
    )
    .unwrap();
    assert!(core::ptr::eq(thermal.original().source(), &water));
    assert_eq!(
        thermal.liquid_thermal_state().unwrap(),
        water.liquid_transport().thermal_transport_state().unwrap()
    );
    let response = thermal.response(10_000.0, 0.1).unwrap();
    assert!(
        response
            .mean_absorption_w_m3(0.1)
            .unwrap()
            .liquid_thermal_w_m3
            > 0.0
    );
}

#[test]
fn conducting_liquid_refuses_stale_sources_overlapping_layers_and_large_temperature_excursions() {
    let material = conducting_material(0.02);
    let original =
        PolydisperseBubblyLiquid::new(&material, distribution(), acoustic_limits()).unwrap();
    let gas = gases(original.population());
    let thermal = ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid(
        &original,
        &gas,
        thermal_limits(),
        liquid_limits(),
        &mut || false,
    )
    .unwrap();
    let accepted = thermal.response(10_000.0, 0.1).unwrap();
    assert!(matches!(
        thermal.response(1e-8, 0.1),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "liquid thermal penetration / mean spacing",
            ..
        })
    ));
    let strict = ConductingLiquidPopulationLimits {
        maximum_fractional_interface_temperature: accepted
            .diagnostics()
            .maximum_fractional_interface_temperature
            * 0.5,
        ..liquid_limits()
    };
    let bounded = ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid(
        &original,
        &gas,
        thermal_limits(),
        strict,
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        bounded.response(10_000.0, 0.1),
        Err(BubblyLiquidError::OutsideModel {
            quantity: "fractional interface temperature",
            ..
        })
    ));
    assert!(
        ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid(
            &original,
            &gas,
            thermal_limits(),
            liquid_limits(),
            &mut || true,
        )
        .is_err()
    );
    assert_eq!(thermal.response(10_000.0, 0.1).unwrap(), accepted);

    let mut stale = material;
    stale.thermal =
        ThermalFluidTransportState::try_new(material.thermal.transport(), 301.0, None).unwrap();
    let stale_population =
        PolydisperseBubblyLiquid::new(&stale, distribution(), acoustic_limits()).unwrap();
    assert!(matches!(
        ThermalPolydisperseBubblyLiquid::new_with_conducting_liquid(
            &stale_population,
            &gas,
            thermal_limits(),
            liquid_limits(),
            &mut || false,
        ),
        Err(BubblyLiquidError::MaterialMismatch)
    ));
}
