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
