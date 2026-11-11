#![cfg(feature = "plane-wave")]
//! G1/G3: original collisional carriers in native Maxwell half-space fields.
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;

#[path = "plane_wave_plasma/magnetized.rs"]
mod magnetized;

#[path = "plane_wave_plasma/layers.rs"]
mod layer_fields;

#[path = "plane_wave_plasma/one_sided.rs"]
mod one_sided;

#[path = "plane_wave_plasma/depth_heat.rs"]
mod depth_heat;

#[path = "plane_wave_plasma/dispersion.rs"]
mod dispersion;

use fs_em::plane_wave::{
    DielectricResponse, IonizedGasWaveResponse, PlanarDielectricInterface, Polarization,
    ReflectionError,
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    ionized_gas::{IonizedGasConductivity, IonizedGasError},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const EPSILON: f64 = 1e-11;
const MU: f64 = 1e-6;

fn gas() -> IonizedGasConductivity {
    let state = fixture::frozen_state(0.001, 6000.0, 12000.0, 1.0);
    IonizedGasConductivity::from_state(
        state.clone(),
        fixture::frozen_rates(&state),
        fixture::limits(),
        || false,
    )
    .unwrap()
}

fn background_at(th: f64, te: f64, pressure: f64, include_mu: bool) -> ResolvedMaterialStatePoint {
    background_with_properties(th, te, pressure, include_mu, EPSILON, None)
}

fn background_with_properties(
    th: f64,
    te: f64,
    pressure: f64,
    include_mu: bool,
    epsilon: f64,
    conductivity: Option<f64>,
) -> ResolvedMaterialStatePoint {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let pressure_type =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let validity = ValidityDomain::unconstrained()
        .with_quantity("T_h", temperature, 2000.0, 20000.0)
        .with_quantity("T_e", temperature, 2000.0, 20000.0)
        .with_quantity("p", pressure_type, 1.0, 1e8);
    let point = QueryPoint::new()
        .with_quantity("T_h", temperature, th)
        .unwrap()
        .with_quantity("T_e", temperature, te)
        .unwrap()
        .with_quantity("p", pressure_type, pressure)
        .unwrap();
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        (
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            epsilon,
        ),
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            MU,
        ),
    ]
    .into_iter()
    .take(if include_mu { 2 } else { 1 })
    .chain(conductivity.map(|value| {
        (
            fs_material::conductor::ELECTRICAL_CONDUCTIVITY_PROPERTY,
            fs_material::conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
            value,
        )
    }))
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic linear material fixture; optional Ohmic conductivity is separately authored".into(),
                    license: "analytic test fixture".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "explicit bound-charge background".into(),
            phase: "gas".into(),
            process: "analytic Maxwell fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn background(gas: &IonizedGasConductivity) -> ResolvedMaterialStatePoint {
    let state = gas.state();
    background_at(
        state.heavy_temperature().value(),
        state.electron_temperature().value(),
        state.pressure().value(),
        true,
    )
}
fn complex_near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-100),
        "{a:?} != {b:?}"
    );
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-100),
        "{a:e} != {b:e}"
    );
}
fn cross(k: [C64; 3], v: [C64; 3]) -> [C64; 3] {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let l = (i + 2) % 3;
        k[j] * v[l] - k[l] * v[j]
    })
}

#[test]
fn g1_plasma_dispersion_and_cutoff_reuse_the_original_native_interface() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let weight: f64 = gas
        .carriers()
        .iter()
        .map(|carrier| {
            let charge = carrier.charge_number() as f64 * fixture::E;
            carrier.number_density_per_m3() * charge * charge / carrier.mass_kg()
        })
        .sum();
    let plasma_frequency = (weight / EPSILON).sqrt();
    let mut reflection = Vec::new();
    for ratio in [0.25_f64, 1.0, 4.0] {
        let omega = ratio * plasma_frequency;
        let law =
            IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
        assert!(std::ptr::eq(law.transport().source(), &gas));
        assert!(std::ptr::eq(law.material(), &background));
        let sigma = gas.carriers().iter().fold(C64::ZERO, |sum, carrier| {
            let charge = carrier.charge_number() as f64 * fixture::E;
            sum + C64::from_re(
                carrier.number_density_per_m3() * charge * charge / carrier.mass_kg(),
            ) / C64::new(carrier.relaxation().rate_per_s(), omega)
        });
        complex_near(law.conductivity_s_m(), sigma, 3e-13);
        let epsilon = C64::new(EPSILON + sigma.im / omega, -sigma.re / omega);
        complex_near(law.effective_permittivity_f_m(), epsilon, 3e-12);
        let interface =
            PlanarDielectricInterface::try_new(&background, DielectricResponse::IonizedGas(law))
                .unwrap();
        let normal = (epsilon.scale(1.0 / EPSILON)).sqrt();
        let expected_r = (C64::ONE - normal) / (C64::ONE + normal);
        for polarization in [Polarization::Te, Polarization::Tm] {
            let result = interface.response(omega, 1.0, polarization).unwrap();
            complex_near(result.reflected_tangential_e, expected_r, 3e-12);
            complex_near(
                result.transmitted_conjugate_h,
                normal.scale(2.0) / (C64::ONE + normal),
                3e-12,
            );
            near(
                result.reflected_power_fraction + result.transmitted_power_fraction,
                1.0,
                3e-13,
            );
            assert!(result.transmitted_normal_wavenumber_rad_m.re > 0.0);
            assert!(result.transmitted_normal_wavenumber_rad_m.im < 0.0);
        }
        reflection.push(
            interface
                .response(omega, 1.0, Polarization::Te)
                .unwrap()
                .reflected_power_fraction,
        );
        if ratio < 1.0 {
            assert!(epsilon.re < 0.0);
        }
        if ratio > 1.0 {
            assert!(epsilon.re > 0.0);
        }
    }
    assert!(reflection[0] > 0.95);
    assert!(reflection[2] < 0.01);
}

#[test]
fn g1_vector_fields_close_maxwell_poynting_and_original_carrier_heat_at_depth() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let amplitude = C64::new(0.003, -0.004);
    for omega in [1e9_f64, 1e12, 1e14] {
        let law =
            IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
        let interface =
            PlanarDielectricInterface::try_new(&background, DielectricResponse::IonizedGas(law))
                .unwrap();
        for cosine in [1.0_f64, 0.5, 0.03] {
            for polarization in [Polarization::Te, Polarization::Tm] {
                let wave = interface
                    .prepare_ionized_gas_wave(cosine, polarization, amplitude, || false)
                    .unwrap();
                assert!(std::ptr::eq(wave.surface_response().source, &gas));
                assert!(wave.power_balance_relative_residual() < 2e-12);
                let kz = wave
                    .interface_response()
                    .transmitted_normal_wavenumber_rad_m;
                let kx = omega * (EPSILON * MU).sqrt() * ((1.0 - cosine) * (1.0 + cosine)).sqrt();
                let k = [C64::from_re(kx), C64::ZERO, kz];
                let entrance = wave.field(0.0).unwrap();
                for optical_depth in [0.0_f64, 0.1, 2.0] {
                    let depth = optical_depth / -kz.im;
                    let field = wave.field(depth).unwrap();
                    let local = transport
                        .response(omega, field.electric_v_m, || false)
                        .unwrap();
                    let k_cross_e = cross(k, field.electric_v_m);
                    let k_cross_h = cross(k, field.magnetic_a_m);
                    for axis in 0..3 {
                        complex_near(
                            k_cross_e[axis],
                            field.magnetic_a_m[axis].scale(omega * MU),
                            2e-11,
                        );
                        complex_near(
                            k_cross_h[axis],
                            -(law.effective_permittivity_f_m() * field.electric_v_m[axis])
                                .scale(omega),
                            2e-11,
                        );
                        complex_near(
                            field.current_density_a_m2[axis],
                            local.current_density_a_m2[axis],
                            3e-12,
                        );
                    }
                    near(
                        field.mean_heating_w_m3(),
                        local.mean_heating_w_m3.total_power_density_w_m3(),
                        3e-12,
                    );
                    near(
                        field.mean_drift_kinetic_energy_j_m3(),
                        local.mean_drift_kinetic_energy_j_m3,
                        3e-12,
                    );
                    near(
                        field.reactive_carrier_power_var_m3(),
                        2.0 * omega * field.mean_drift_kinetic_energy_j_m3(),
                        3e-12,
                    );
                    near(
                        (-2.0 * kz.im) * field.complex_poynting_w_m2[2].re,
                        field.mean_heating_w_m3(),
                        3e-12,
                    );
                    near(
                        field.absorbed_power_to_depth_w_m2 + field.complex_poynting_w_m2[2].re,
                        entrance.complex_poynting_w_m2[2].re,
                        3e-12,
                    );
                    assert!(field.mean_background_field_energy_j_m3() > 0.0);
                    for i in 0..gas.carriers().len() {
                        for (actual, expected) in field
                            .carrier_velocity_m_s(i)
                            .unwrap()
                            .iter()
                            .zip(local.carrier_velocity_m_s[i])
                        {
                            complex_near(*actual, expected, 3e-12);
                        }
                    }
                    if polarization == Polarization::Tm && cosine < 1.0 {
                        assert!(field.electric_v_m[2].abs() > 0.0);
                    }
                }
                let opaque = wave.field(900.0 / -kz.im).unwrap();
                assert_eq!(opaque.electric_v_m, [C64::ZERO; 3]);
                assert_eq!(opaque.mean_heating_w_m3(), 0.0);
                near(
                    opaque.absorbed_power_to_depth_w_m2,
                    entrance.complex_poynting_w_m2[2].re,
                    3e-13,
                );
            }
        }
    }
}

#[test]
fn g1_propagation_scales_match_the_original_complex_wave_number_and_fields() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    for omega in [1e9_f64, 1e12, 1e14] {
        let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
        let interface = PlanarDielectricInterface::try_new(
            &background, DielectricResponse::IonizedGas(law),
        ).unwrap();
        for cosine in [1.0_f64, 0.5, 0.03] {
            let wave = interface.prepare_ionized_gas_wave(
                cosine, Polarization::Te, C64::new(0.003, -0.004), || false,
            ).unwrap();
            let kz = wave.interface_response().transmitted_normal_wavenumber_rad_m;
            let scales = wave.propagation_scales();
            let wavelength = scales.normal_phase_wavelength_m.unwrap();
            let phase_speed = scales.normal_phase_velocity_m_s.unwrap();
            let amplitude_depth = scales.amplitude_skin_depth_m.unwrap();
            let intensity_depth = scales.intensity_skin_depth_m.unwrap();
            near(wavelength * kz.re, core::f64::consts::TAU, 3e-15);
            near(phase_speed * kz.re, omega, 3e-15);
            near(amplitude_depth * -kz.im, 1.0, 3e-15);
            near(intensity_depth * -2.0 * kz.im, 1.0, 3e-15);
            near(
                scales.phase_cycles_per_amplitude_efold.unwrap(),
                kz.re / (-core::f64::consts::TAU * kz.im),
                3e-15,
            );

            let entrance = wave.field(0.0).unwrap();
            let one_cycle = wave.field(wavelength).unwrap();
            let axis = entrance.electric_v_m.iter().position(|e| e.abs() > 0.0).unwrap();
            let ratio = one_cycle.electric_v_m[axis] / entrance.electric_v_m[axis];
            // Over one normal wavelength the outgoing factor exp(-i kz z) is
            // exp(kz.im lambda) exp(-i 2 pi): real. A nonzero imaginary part is
            // rounding alone, so it is judged against |ratio|, not against 0.
            // The bound counts the operations that build the ratio, in units
            // of u = 2^-53 times |ratio|:
            // - lambda = fl(TAU / kz.re) moves the exact DD phase kz.re lambda
            //   from 2 pi by at most u TAU + |TAU - 2 pi| < 8.5u;
            // - sin_cos_dd's cosine and sine round to f64 within u/2;
            // - the two products by the attenuation add u;
            // - the complex product with the surface field (Kahan's FMA form,
            //   componentwise 2u, Jeannerod et al., Math. Comp. 86, 2017) 2u;
            // - Smith's quotient as fs-math's C64 evaluates it with FMAs,
            //   under 8u to first order (3u in the denominator, 2u per
            //   numerator component plus 2u from its r-products, u per final
            //   division);
            // - det::exp against the libm exp of the same rounded argument,
            //   its declared 3 ulp budget, at most 6u.
            // That is 26u; 32u covers the second-order terms.
            complex_near(
                ratio,
                C64::new((kz.im * wavelength).exp(), 0.0),
                32.0 * (f64::EPSILON / 2.0),
            );

            let one_intensity = wave.field(intensity_depth).unwrap();
            near(
                one_intensity.complex_poynting_w_m2[2].re
                    / entrance.complex_poynting_w_m2[2].re,
                (-1.0_f64).exp(),
                3e-12,
            );
        }
    }
}

#[test]
fn g3_coulomb_drag_and_signed_thermal_partition_survive_wave_propagation() {
    let original = gas();
    let gas = original
        .clone()
        .with_coulomb_drag(fixture::coulomb_model(original.state()), || false)
        .unwrap();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let omega = 1e9_f64;
    let law = IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).unwrap();
    let interface =
        PlanarDielectricInterface::try_new(&background, DielectricResponse::IonizedGas(law))
            .unwrap();
    let wave = interface
        .prepare_ionized_gas_wave(0.3, Polarization::Tm, C64::new(0.01, 0.02), || false)
        .unwrap();
    let kz = wave
        .interface_response()
        .transmitted_normal_wavenumber_rad_m;
    let field = wave.field(0.5 / -kz.im).unwrap();
    let local = transport
        .response(omega, field.electric_v_m, || false)
        .unwrap();
    near(
        field.mean_electron_heating_w_m3(),
        local.mean_heating_w_m3.electron_power_density_w_m3(),
        3e-11,
    );
    near(
        field.mean_heavy_heating_w_m3(),
        local.mean_heating_w_m3.heavy_power_density_w_m3(),
        3e-11,
    );
    near(
        field.mean_exchange_correction_to_electrons_w_m3(),
        local
            .mean_heating_w_m3
            .thermal_exchange_correction_to_electrons_w_m3(),
        3e-11,
    );
    assert!(field.mean_exchange_correction_to_electrons_w_m3() > 0.0);
    near(
        field.mean_electron_heating_w_m3() + field.mean_heavy_heating_w_m3(),
        field.mean_heating_w_m3(),
        3e-12,
    );
}

#[test]
fn g0_original_conditions_frequency_and_actual_surface_drift_are_required() {
    let gas = gas();
    let transport = gas.prepare_harmonic(|| false).unwrap();
    let background = background(&gas);
    let p = gas.state().pressure().value();
    for wrong in [
        background_at(6001.0, 12000.0, p, true),
        background_at(6000.0, 12001.0, p, true),
        background_at(6000.0, 12000.0, p * 1.01, true),
        background_at(6000.0, 12000.0, p, false),
    ] {
        assert!(IonizedGasWaveResponse::try_new(&transport, &wrong, 1e9, || false).is_err());
    }
    for omega in [0.0_f64, -1.0, f64::INFINITY, f64::NAN] {
        assert!(IonizedGasWaveResponse::try_new(&transport, &background, omega, || false).is_err());
    }
    assert!(matches!(
        IonizedGasWaveResponse::try_new(&transport, &background, 1e9, || true),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
    let law = IonizedGasWaveResponse::try_new(&transport, &background, 1e9, || false).unwrap();
    let interface =
        PlanarDielectricInterface::try_new(&background, DielectricResponse::IonizedGas(law))
            .unwrap();
    assert!(interface.response(2e9, 1.0, Polarization::Te).is_err());
    assert!(
        interface
            .prepare_ionized_gas_wave(1.0, Polarization::Te, C64::from_re(1e12), || false)
            .is_err()
    );
    assert!(matches!(
        interface.prepare_ionized_gas_wave(0.5, Polarization::Tm, C64::ONE, || true),
        Err(ReflectionError::IonizedGas(IonizedGasError::Cancelled))
    ));
    let silence = interface
        .prepare_ionized_gas_wave(0.5, Polarization::Tm, C64::ZERO, || false)
        .unwrap();
    let point = silence.field(1.0).unwrap();
    assert_eq!(point.mean_heating_w_m3(), 0.0);
    assert_eq!(point.current_density_a_m2, [C64::ZERO; 3]);
    assert!(point.carrier_velocity_m_s(gas.carriers().len()).is_none());
    for depth in [-1.0_f64, f64::INFINITY, f64::NAN] {
        assert!(silence.field(depth).is_err());
    }
    // Generic films/beams cannot bypass the whole-field weak-drift admission.
    assert!(
        PlanarDielectricLayers::try_new(
            &background,
            DielectricResponse::Lossless(&background),
            1e9,
            &[DielectricLayer {
                response: DielectricResponse::IonizedGas(law),
                thickness_m: 1.0
            }]
        )
        .is_err()
    );
}
