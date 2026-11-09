use super::*;
use crate::cavitation::{BubbleGas, BubbleLimits, KellerMiksisLimits};
use fs_material::{cavitation::CavitationState, gas::GasSpec};

const RADIUS: f64 = 50.0e-6;

fn limits() -> ThermalBubbleLimits {
    ThermalBubbleLimits {
        maximum_gas_wavenumber_radius: 0.1,
        maximum_gas_liquid_density_ratio: 0.01,
        maximum_vapor_pressure_fraction: 0.05,
    }
}

fn bubble(kappa: f64, viscosity: f64) -> KellerMiksisBubble {
    let liquid =
        CavitationState::try_new(293.15, 1000.0, viscosity, 0.072, 2000.0, 1500.0).unwrap();
    KellerMiksisBubble::try_new(
        liquid,
        BubbleGas::equilibrium(liquid, RADIUS, 101_325.0, kappa).unwrap(),
        KellerMiksisLimits {
            bubble: BubbleLimits {
                minimum_radius_m: 0.1 * RADIUS,
                maximum_radius_m: 3.0 * RADIUS,
                maximum_wall_mach_number: 0.1,
                maximum_vapor_density_ratio: 0.01,
            },
            maximum_pressure_ratio: 0.01,
        },
    )
    .unwrap()
}

fn gas_at(bubble: KellerMiksisBubble, radius: f64) -> GasState {
    GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        bubble.liquid().temperature_k(),
        bubble.gas().pressure_at_radius(radius).unwrap(),
    )
    .unwrap()
}

fn model() -> IsothermalWallBubble {
    let bubble = bubble(1.0, 0.001);
    IsothermalWallBubble::new(bubble, RADIUS, gas_at(bubble, RADIUS), limits()).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()),
        "{actual:.17e} != {expected:.17e}"
    );
}

// Independent eigenfunctions of the spherical Dirichlet heat equation:
// F = 1 + 6 sum_n i*x/[a_n(a_n-i*x)], a_n=pi^2*n^2.
// The omitted tail is bounded by 2*x/(pi^4*N^3). No coth or exponential is used.
fn heat_eigenmodes(gamma: f64, x: f64) -> (C64, f64) {
    const N: usize = 32768;
    let mut sum = C64::ZERO;
    let mut compensation = C64::ZERO;
    for n in 1..=N {
        let a = PI * PI * (n as f64) * (n as f64);
        let denominator = a * a + x * x;
        let term = C64::new(-x * x / (a * denominator), x / denominator);
        let adjusted = term - compensation;
        let next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    let f = C64::ONE + sum.scale(6.0);
    let kappa = C64::from_re(gamma) / (C64::ONE + f.scale(gamma - 1.0));
    let error = gamma * (gamma - 1.0) * 2.0 * x / (PI.powi(4) * (N as f64).powi(3));
    (kappa, error)
}

#[test]
fn heat_diffusion_matches_independent_spherical_eigenmodes() {
    for gamma in [1.01, 1.4, 5.0 / 3.0] {
        for x in [
            1e-6, 0.01, 0.25, 0.25000001, 1.0, 10.0, 100.0, 800.0, 1000.0, 1e4,
        ] {
            let actual = effective_polytropic_exponent(gamma, x);
            let (expected, tail) = heat_eigenmodes(gamma, x);
            assert!(
                (actual - expected).abs() <= tail * 1.01 + 3e-14,
                "gamma={gamma}, x={x}: {actual:?} != {expected:?}, tail={tail}"
            );
            assert!(actual.re >= 1.0 && actual.re <= gamma);
            assert!(actual.im < 0.0);
        }
    }
}

#[test]
fn isothermal_adiabatic_and_branch_limits_preserve_small_loss() {
    for gamma in [1.01, 1.4, 5.0 / 3.0] {
        assert_eq!(effective_polytropic_exponent(gamma, 0.0), C64::ONE);
        for x in [1e-100, 1e-16, 1e-8, 1e-5] {
            let actual = effective_polytropic_exponent(gamma, x);
            close(actual.im, -(gamma - 1.0) * x / (15.0 * gamma), 2e-11);
            close(actual.re, 1.0, 2e-11);
        }
        for x in [1e16, 1e100, 1e308] {
            let actual = effective_polytropic_exponent(gamma, x);
            close(actual.re, gamma, 2e-8);
            // Large-frequency loss is O(x^-1/2), not exactly zero.
            close(
                actual.im,
                -3.0 * gamma * (gamma - 1.0) / (2.0_f64.sqrt() * x.sqrt()),
                5e-8,
            );
        }
        for x in [0.25_f64, 800.0] {
            let left = effective_polytropic_exponent(gamma, x);
            let right = effective_polytropic_exponent(gamma, x.next_up());
            assert!(
                (left - right).abs() < 3e-14,
                "branch at {x}: {left:?}, {right:?}"
            );
        }
    }
}

#[test]
fn current_equilibrium_retains_inventory_and_native_dc_compliance() {
    let bubble = bubble(1.0, 0.001);
    for radius in [RADIUS, 1.1 * RADIUS] {
        let gas = gas_at(bubble, radius);
        let thermal = IsothermalWallBubble::new(bubble, radius, gas, limits()).unwrap();
        assert_eq!(thermal.bubble(), bubble);
        assert_eq!(thermal.gas(), gas);
        assert_eq!(
            bubble.gas().pressure_at_radius(radius).unwrap(),
            gas.pressure
        );
        assert_eq!(
            thermal.linearization(),
            bubble.linearization_at_radius(radius).unwrap()
        );
        assert_eq!(thermal.limits(), limits());
        assert_eq!(
            thermal.thermal_diffusivity_m2_s(),
            gas.transport_state().unwrap().thermal_diffusivity_m2_s()
        );
        let dc = thermal.response(0.0).unwrap();
        assert_eq!(dc.effective_polytropic_exponent, C64::ONE);
        assert_eq!(
            dc.pressure_transfer_s2,
            C64::from_re(1.0 / thermal.linearization().stiffness_coefficient_per_s2)
        );
        assert_eq!(dc.thermal_diffusion_parameter, 0.0);
        assert_eq!(dc.gas_wavenumber_radius, 0.0);
        assert_eq!(dc.thermal_power_per_fractional_radius_squared_w, 0.0);
    }
    assert!(matches!(
        IsothermalWallBubble::new(bubble, 1.1 * RADIUS, gas_at(bubble, RADIUS), limits()),
        Err(ThermalBubbleError::StateMismatch(
            "gas equilibrium partial pressure"
        ))
    ));
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(bubble.gas().pressure_at_radius(radius).is_err());
    }
    let tiny = BubbleGas::try_new(1.0, 1.0, 1.0).unwrap();
    assert!(matches!(
        tiny.pressure_at_radius(1e200),
        Err(CavitationError::NonFiniteDerived(_))
    ));
    let empty = BubbleGas::try_new(1.0, 0.0, 1.0).unwrap();
    assert_eq!(empty.pressure_at_radius(1e200).unwrap(), 0.0);
}

#[test]
fn dynamic_operator_retains_native_inertia_and_pressure_derivative() {
    let thermal = model();
    let linear = thermal.linearization();
    let gas = thermal.gas();
    let liquid = thermal.bubble().liquid();
    let inertia = liquid.liquid_density_kg_m3() * RADIUS * RADIUS;
    let viscous = 4.0 * liquid.dynamic_viscosity_pa_s() / inertia;
    for omega in [1e-7, 1.0, 1e4, 1e5, 6e5] {
        let response = thermal.response(omega).unwrap();
        let kappa = response.effective_polytropic_exponent;
        let stiffness = C64::new(
            (3.0 * gas.pressure * kappa.re - 2.0 * liquid.surface_tension_n_m() / RADIUS) / inertia,
            3.0 * gas.pressure * kappa.im / inertia,
        );
        let numerator = C64::new(1.0, -omega * linear.pressure_derivative_time_s);
        let denominator = numerator * stiffness
            - C64::new(linear.inertia_multiplier * omega * omega, viscous * omega);
        let expected = numerator / denominator;
        assert!((response.pressure_transfer_s2 - expected).abs() < 3e-14 * expected.abs());
        assert!(response.pressure_transfer_s2.im > 0.0);
        assert!(response.thermal_power_per_fractional_radius_squared_w > 0.0);
        close(
            response.thermal_diffusion_parameter,
            omega * RADIUS * RADIUS / thermal.thermal_diffusivity_m2_s(),
            5e-16,
        );
    }
}

#[test]
fn thermal_power_equals_independent_cycle_average_of_gas_work() {
    let thermal = model();
    let omega = 2e5;
    let response = thermal.response(omega).unwrap();
    let kappa = response.effective_polytropic_exponent;
    let volume = 4.0 * PI * RADIUS.powi(3) / 3.0;
    let mut work = 0.0;
    for i in 0..1024 {
        let (sin, cos) = (2.0 * PI * (i as f64) / 1024.0).sin_cos();
        let pressure = -3.0 * thermal.gas().pressure * (kappa.re * cos + kappa.im * sin);
        let volume_rate = -3.0 * volume * omega * sin;
        work -= pressure * volume_rate / 1024.0;
    }
    close(
        response.thermal_power_per_fractional_radius_squared_w,
        work,
        3e-14,
    );
}

#[test]
fn material_inventory_stability_and_thermal_domain_refusals() {
    let native = bubble(1.0, 0.001);
    let gas = gas_at(native, RADIUS);
    let mut invalid = limits();
    invalid.maximum_gas_wavenumber_radius = 0.10001;
    assert!(invalid.validate().is_err());
    let mut invalid = limits();
    invalid.maximum_gas_liquid_density_ratio = 1e-4;
    assert!(matches!(
        IsothermalWallBubble::new(native, RADIUS, gas, invalid),
        Err(ThermalBubbleError::OutsideModel { .. })
    ));
    let mut invalid = limits();
    invalid.maximum_vapor_pressure_fraction = 0.001;
    assert!(matches!(
        IsothermalWallBubble::new(native, RADIUS, gas, invalid),
        Err(ThermalBubbleError::OutsideModel { .. })
    ));
    let mut incoherent = gas;
    incoherent.thermal_conductivity *= 2.0;
    assert!(matches!(
        IsothermalWallBubble::new(native, RADIUS, incoherent, limits()),
        Err(ThermalBubbleError::Material(_))
    ));
    let wrong_temperature = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        gas.temperature + 1.0,
        gas.pressure,
    )
    .unwrap();
    assert!(matches!(
        IsothermalWallBubble::new(native, RADIUS, wrong_temperature, limits()),
        Err(ThermalBubbleError::StateMismatch(_))
    ));
    let mut wet = gas;
    wet.water_mole_fraction = 0.01;
    assert!(matches!(
        IsothermalWallBubble::new(native, RADIUS, wet, limits()),
        Err(ThermalBubbleError::InvalidInput(_))
    ));
    let adiabatic = bubble(1.4, 0.001);
    assert!(matches!(
        IsothermalWallBubble::new(adiabatic, RADIUS, gas_at(adiabatic, RADIUS), limits()),
        Err(ThermalBubbleError::InvalidInput(_))
    ));
    let unstable = KellerMiksisBubble::try_new(
        native.liquid(),
        BubbleGas::try_new(RADIUS, 100.0, 1.0).unwrap(),
        native.compressibility_limits(),
    )
    .unwrap();
    assert!(matches!(
        IsothermalWallBubble::new(unstable, RADIUS, gas_at(unstable, RADIUS), limits()),
        Err(ThermalBubbleError::Bubble(_))
    ));
    let soft_liquid = CavitationState::try_new(293.15, 1000.0, 0.001, 0.072, 2000.0, 10.0).unwrap();
    let soft =
        KellerMiksisBubble::try_new(soft_liquid, native.gas(), native.compressibility_limits())
            .unwrap();
    assert!(matches!(
        IsothermalWallBubble::new(soft, 1.1 * RADIUS, gas_at(soft, 1.1 * RADIUS), limits()),
        Err(ThermalBubbleError::Bubble(
            CavitationError::PressureRatioExceeded { .. }
        ))
    ));
    let model = model();
    assert!(model.response(-1.0).is_err());
    assert!(model.response(f64::NAN).is_err());
    let too_fast = 1.01 * limits().maximum_gas_wavenumber_radius * gas.sound_speed / RADIUS;
    assert!(matches!(
        model.response(too_fast),
        Err(ThermalBubbleError::OutsideModel {
            quantity: "gas compactness",
            ..
        })
    ));
}

fn liquid_thermal(conductivity: f64) -> ThermalFluidTransportState {
    ThermalFluidTransportState::try_new(
        fs_material::fluid_transport::FluidTransportState::try_new(
            1000.0,
            4200.0,
            conductivity,
            0.001,
        )
        .unwrap(),
        293.15,
        None,
    )
    .unwrap()
}

#[test]
fn conducting_liquid_preserves_state_and_dc_compliance() {
    let original = model();
    let liquid = liquid_thermal(0.6);
    let thermal = ConductingLiquidBubble::new(original, liquid).unwrap();
    assert_eq!(thermal.base(), &original);
    assert_eq!(thermal.liquid(), liquid);
    let dc = thermal.response(0.0).unwrap();
    assert_eq!(dc.bubble, original.response(0.0).unwrap());
    assert_eq!(dc.interface_temperature_per_fractional_radius_k, C64::ZERO);
    assert_eq!(dc.outward_heat_per_fractional_radius_w, C64::ZERO);
    assert_eq!(dc.gas_thermal_power_per_fractional_radius_squared_w, 0.0);
    assert_eq!(dc.liquid_thermal_power_per_fractional_radius_squared_w, 0.0);
    assert_eq!(dc.liquid_thermal_diffusion_parameter, 0.0);
    for (temperature, density, viscosity) in [
        (294.15, 1000.0, 0.001),
        (293.15, 999.0, 0.001),
        (293.15, 1000.0, 0.002),
    ] {
        let wrong = ThermalFluidTransportState::try_new(
            fs_material::fluid_transport::FluidTransportState::try_new(
                density, 4200.0, 0.6, viscosity,
            )
            .unwrap(),
            temperature,
            None,
        )
        .unwrap();
        assert!(matches!(
            ConductingLiquidBubble::new(original, wrong),
            Err(ThermalBubbleError::StateMismatch(
                "surrounding liquid thermal state"
            ))
        ));
    }
    assert!(thermal.response(-1.0).is_err());
    assert!(thermal.response(f64::NAN).is_err());
    assert!(thermal.response(f64::INFINITY).is_err());
    assert!(matches!(
        thermal.response(
            1.01 * limits().maximum_gas_wavenumber_radius * original.gas().sound_speed / RADIUS
        ),
        Err(ThermalBubbleError::OutsideModel {
            quantity: "gas compactness",
            ..
        })
    ));
}

#[test]
fn finite_liquid_recovers_isothermal_wall_and_insulated_gas_limits() {
    let original = model();
    let conducting = ConductingLiquidBubble::new(original, liquid_thermal(1e20)).unwrap();
    let insulating = ConductingLiquidBubble::new(original, liquid_thermal(1e-30)).unwrap();
    for omega in [1.0, 1e3, 1e5, 6e5] {
        let fixed = original.response(omega).unwrap();
        let wet = conducting.response(omega).unwrap();
        assert!(
            (wet.bubble.effective_polytropic_exponent - fixed.effective_polytropic_exponent).abs()
                < 2e-14
        );
        assert!(
            (wet.bubble.pressure_transfer_s2 - fixed.pressure_transfer_s2).abs()
                < 2e-14 * fixed.pressure_transfer_s2.abs()
        );
        assert!(wet.interface_temperature_per_fractional_radius_k.abs() < 1e-15);
        close(
            wet.gas_thermal_power_per_fractional_radius_squared_w,
            fixed.thermal_power_per_fractional_radius_squared_w,
            2e-14,
        );
        let insulated = insulating.response(omega).unwrap();
        close(
            insulated.bubble.effective_polytropic_exponent.re,
            original.gas().gamma,
            1e-10,
        );
        assert!(insulated.bubble.effective_polytropic_exponent.im <= 0.0);
        let expected_temperature = -3.0 * original.gas().temperature * (original.gas().gamma - 1.0);
        close(
            insulated.interface_temperature_per_fractional_radius_k.re,
            expected_temperature,
            2e-10,
        );
        assert!(insulated.outward_heat_per_fractional_radius_w.abs() < 1e-11);
    }
}

#[test]
fn gas_and_liquid_entropy_losses_close_original_pressure_volume_work() {
    for conductivity in [1e-12, 1e-4, 0.02, 0.6, 1e8] {
        let thermal = ConductingLiquidBubble::new(model(), liquid_thermal(conductivity)).unwrap();
        for omega in [1e-96, 1e-12, 1e-4, 1.0, 1e3, 1e5, 6e5] {
            let response = thermal.response(omega).unwrap();
            let gas_loss = response.gas_thermal_power_per_fractional_radius_squared_w;
            let liquid_loss = response.liquid_thermal_power_per_fractional_radius_squared_w;
            assert!(gas_loss > 0.0 && liquid_loss > 0.0);
            assert!(response.bubble.effective_polytropic_exponent.im < 0.0);
            close(
                gas_loss + liquid_loss,
                response
                    .bubble
                    .thermal_power_per_fractional_radius_squared_w,
                3e-13,
            );
            // Complex heat exchange cancels at the interface, without using
            // the native harmonic stiffness or mechanical work expression.
            let liquid_a = (0.5 * omega * RADIUS * RADIUS
                / thermal.liquid().transport().thermal_diffusivity_m2_s())
            .sqrt();
            let heat_from_liquid = (C64::new(1.0 + liquid_a, -liquid_a)
                * response.interface_temperature_per_fractional_radius_k)
                .scale(4.0 * PI * RADIUS * conductivity);
            assert!(
                (response.outward_heat_per_fractional_radius_w - heat_from_liquid).abs()
                    <= 5e-14 * heat_from_liquid.abs()
            );
        }
    }
}

fn complex_exp(z: C64) -> C64 {
    let (sin, cos) = z.im.sin_cos();
    C64::new(cos, sin).scale(z.re.exp())
}

fn simpson_integral(f: impl Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    const N: usize = 4096;
    let h = (hi - lo) / N as f64;
    let mut total = f(lo) + f(hi);
    for n in 1..N {
        total += (if n % 2 == 0 { 2.0 } else { 4.0 }) * f(lo + n as f64 * h);
    }
    total * h / 3.0
}

#[test]
fn independent_radial_heat_fields_reproduce_each_entropy_integral() {
    let thermal = ConductingLiquidBubble::new(model(), liquid_thermal(0.02)).unwrap();
    let gas = thermal.base().gas();
    for omega in [10.0, 1e3, 1e5] {
        let response = thermal.response(omega).unwrap();
        let x = omega * RADIUS * RADIUS / thermal.base().thermal_diffusivity_m2_s();
        let z = C64::new((0.5 * x).sqrt(), -(0.5 * x).sqrt());
        let sinh_z = (complex_exp(z) - complex_exp(-z)).scale(0.5);
        // Directly integrate k|grad T|^2/(2T0) over the gas volume. This
        // uses neither z*coth(z)-1 nor the production admittance/heat average.
        let gradient_integral = simpson_integral(
            |s| {
                if s == 0.0 {
                    return 0.0;
                }
                let plus = complex_exp(z.scale(s));
                let minus = complex_exp(z.scale(-s));
                let sinh = (plus - minus).scale(0.5);
                let cosh = (plus + minus).scale(0.5);
                let derivative = (z * cosh - sinh.scale(1.0 / s)) / sinh_z;
                derivative.norm_sq()
            },
            0.0,
            1.0,
        );
        let adiabatic = response
            .bubble
            .effective_polytropic_exponent
            .scale(-3.0 * gas.temperature * (gas.gamma - 1.0) / gas.gamma);
        let drop = adiabatic - response.interface_temperature_per_fractional_radius_k;
        let gas_loss = 2.0 * PI * RADIUS * gas.thermal_conductivity / gas.temperature
            * drop.norm_sq()
            * gradient_integral;
        close(
            gas_loss,
            response.gas_thermal_power_per_fractional_radius_squared_w,
            3e-9,
        );
        // Exterior coordinate t=Re(z_l)*(r/R-1) resolves the diffusion
        // layer independently of its thickness; the omitted tail is <1e-17.
        let a = (0.5 * response.liquid_thermal_diffusion_parameter).sqrt();
        let exterior_integral = simpson_integral(
            |t| {
                let s = 1.0 + t / a;
                C64::new(a + 1.0 / s, -a).norm_sq() * (-2.0 * t).exp() / a
            },
            0.0,
            24.0,
        );
        let liquid_loss = 2.0 * PI * RADIUS * 0.02 / gas.temperature
            * response
                .interface_temperature_per_fractional_radius_k
                .norm_sq()
            * exterior_integral;
        close(
            liquid_loss,
            response.liquid_thermal_power_per_fractional_radius_squared_w,
            3e-8,
        );
    }
}

#[test]
fn finite_liquid_changes_real_resonance_and_thermal_loss_in_original_km_operator() {
    let original = model();
    let thermal = ConductingLiquidBubble::new(original, liquid_thermal(1e-5)).unwrap();
    let omega = 2e5;
    let response = thermal.response(omega).unwrap();
    let fixed = original.response(omega).unwrap();
    assert!(
        (response.bubble.effective_polytropic_exponent - fixed.effective_polytropic_exponent).abs()
            > 0.01
    );
    assert!(
        (response.bubble.pressure_transfer_s2 - fixed.pressure_transfer_s2).abs()
            > 0.01 * fixed.pressure_transfer_s2.abs()
    );
    let liquid = original.bubble().liquid();
    let inertia = liquid.liquid_density_kg_m3() * RADIUS * RADIUS;
    let viscous = 4.0 * liquid.dynamic_viscosity_pa_s() / inertia;
    let kappa = response.bubble.effective_polytropic_exponent;
    let stiffness = (kappa.scale(3.0 * original.gas().pressure)
        - C64::from_re(2.0 * liquid.surface_tension_n_m() / RADIUS))
    .scale(1.0 / inertia);
    let numerator = C64::new(
        1.0,
        -omega * original.linearization().pressure_derivative_time_s,
    );
    let denominator = numerator * stiffness
        - C64::new(
            original.linearization().inertia_multiplier * omega * omega,
            viscous * omega,
        );
    let expected = numerator / denominator;
    assert!((response.bubble.pressure_transfer_s2 - expected).abs() < 3e-14 * expected.abs());
    let volume = 4.0 * PI * RADIUS.powi(3) / 3.0;
    let mut work = 0.0;
    for i in 0..1024 {
        let (sin, cos) = (2.0 * PI * i as f64 / 1024.0).sin_cos();
        let pressure = -3.0 * original.gas().pressure * (kappa.re * cos + kappa.im * sin);
        let volume_rate = -3.0 * volume * omega * sin;
        work -= pressure * volume_rate / 1024.0;
    }
    close(
        response.gas_thermal_power_per_fractional_radius_squared_w
            + response.liquid_thermal_power_per_fractional_radius_squared_w,
        work,
        3e-13,
    );
}
