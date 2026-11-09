//! Coupled pore diffusion and molecular calorics: analytic limits and energy.
//! Reservoir frequencies and pore scaling are manufactured numerical fixtures,
//! not measured CO2/air attenuation or a new qualification of any foam.
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::layers::{
    AngularSpectrumBudget, IncidentPlaneWave, Layer, LayerAbsorptionBudget,
    RigidBackedLayers, TransmittingLayers,
};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    RelaxationState, RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;

fn relaxing_gas(modes: &[(f64, f64)]) -> RelaxationState {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 101_325.0).unwrap();
    let species = SpeciesId::new("synthetic-reservoir-gas").unwrap();
    let medium = RelaxingMedium::new(gas, vec![(species.clone(), 1.0)]).unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species,
        processes: modes
            .iter()
            .map(|&(theta_k, tau)| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: gas.pressure * tau,
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 400.0,
                },
            })
            .collect(),
    }])
    .unwrap()
    .evaluate(&medium)
    .unwrap()
}

fn pores(lafarge: bool, scale: f64) -> PorousFluidLaw {
    PorousFluidLaw::try_new(
        PorousCard {
            name: "synthetic-relaxing-pores",
            porosity: 0.45,
            flow_resistivity: 1.81e-5 / (2.0e-10 * scale * scale),
            tortuosity: 1.4,
            viscous_length_m: 100.0e-6 * scale,
            thermal_length_m: 200.0e-6 * scale,
            thermal_permeability_m2: lafarge.then_some(3.0e-10 * scale * scale),
            provenance: "manufactured coupled caloric limit fixture; no measurement claim",
            valid_hz: (1.0e-9, 1.0e16),
        },
        1.81e-5,
    )
    .unwrap()
}

fn relative(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:?} != {expected:?} at relative tolerance {tolerance}"
    );
}

#[test]
fn absent_and_zero_strength_reservoirs_preserve_original_response_exactly() {
    let absent = relaxing_gas(&[]);
    // A thermally unpopulated vibrational level has exactly zero capacity.
    let zero = relaxing_gas(&[(1.0e8, 2.0e-4)]);
    assert_eq!(zero.processes().len(), 1);
    assert_eq!(zero.processes()[0].heat_capacity_j_kg_k, 0.0);
    for lafarge in [false, true] {
        let law = pores(lafarge, 1.0);
        for relaxation in [&absent, &zero] {
            for omega in [1.0e-6, 50.0, 5.0e3, 5.0e5, 1.0e14] {
                assert_eq!(
                    law.effective_relaxing(relaxation, omega).unwrap(),
                    law.effective_fluid(relaxation.gas(), omega).unwrap()
                );
            }
        }
    }
}

#[test]
fn isothermal_pores_suppress_molecular_compression_loss_at_its_peak() {
    // Slow synthetic relaxation keeps omega*tau=1 while pore-wall diffusion
    // equilibrates temperature. Independent free-gas attenuation would survive
    // this limit incorrectly.
    let relaxation = relaxing_gas(&[(1000.0, 1.0e7)]);
    let gas = relaxation.gas();
    let omega = 1.0e-7;
    let free = relaxation.caloric_response(omega).unwrap().bulk_modulus_pa;
    assert!(-free.im / free.re > 0.02);
    for lafarge in [false, true] {
        let law = pores(lafarge, 1.0);
        let response = law.effective_relaxing(&relaxation, omega).unwrap();
        let isothermal = C64::from_re(gas.pressure / law.reference_card().porosity);
        relative(response.dynamic_bulk_modulus_pa(), isothermal, 2.0e-11);
        relative(
            response.dynamic_bulk_modulus_pa(),
            law.effective_fluid(gas, omega).unwrap().dynamic_bulk_modulus_pa(),
            2.0e-13,
        );
    }
}

#[test]
fn wide_pores_recover_molecular_dispersion_and_high_frequency_freezing() {
    let relaxation = relaxing_gas(&[(1000.0, 2.0e-4)]);
    for lafarge in [false, true] {
        // Geometrical scaling is an asymptotic algebra test, not a physical
        // claim for a finite specimen with kilometre-scale pore dimensions.
        let wide = pores(lafarge, 1.0e8);
        let phi = wide.reference_card().porosity;
        for omega in [50.0, 5000.0, 5.0e5] {
            let expected = relaxation
                .caloric_response(omega)
                .unwrap()
                .bulk_modulus_pa
                .conj()
                .scale(1.0 / phi);
            relative(
                wide.effective_relaxing(&relaxation, omega)
                    .unwrap()
                    .dynamic_bulk_modulus_pa(),
                expected,
                5.0e-8,
            );
        }
        let law = pores(lafarge, 1.0);
        let gas = relaxation.gas();
        let frozen = gas.pressure
            * (1.0 + gas.specific_gas_constant / relaxation.frozen_cv_j_kg_k())
            / phi;
        relative(
            law.effective_relaxing(&relaxation, 1.0e16)
                .unwrap()
                .dynamic_bulk_modulus_pa(),
            C64::from_re(frozen),
            1.0e-6,
        );
    }
}

#[test]
fn finite_pore_crossover_recovers_equilibrated_and_frozen_caloric_endpoints() {
    // Both sources have the same populated reservoir. Only its time differs:
    // instantaneous exchange reproduces equilibrium Cp, while frozen exchange
    // changes BOTH the adiabatic modulus and the thermal penetration depth.
    let equilibrated = relaxing_gas(&[(1000.0, 1.0e-18)]);
    let frozen = relaxing_gas(&[(1000.0, 1.0e18)]);
    assert!(frozen.processes()[0].heat_capacity_j_kg_k > 0.0);
    let gas = frozen.gas();
    let frozen_cv = frozen.frozen_cv_j_kg_k();
    let frozen_cp = frozen_cv + gas.specific_gas_constant;
    let frozen_fluid = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(gas.density, frozen.frozen_sound_speed_m_s()).unwrap(),
        FluidTransportState::try_new(
            gas.density,
            frozen_cp,
            gas.thermal_conductivity,
            gas.dynamic_viscosity,
        )
        .unwrap(),
        frozen_cp / frozen_cv,
    )
    .unwrap();

    for lafarge in [false, true] {
        let law = pores(lafarge, 1.0);
        let card = law.reference_card();
        let permeability = card.thermal_permeability_m2.unwrap_or(
            card.porosity * card.thermal_length_m * card.thermal_length_m / 8.0,
        );
        let crossover = gas.thermal_conductivity * card.porosity
            / (gas.density * gas.specific_heat_cp * permeability);
        for ratio in [0.1, 1.0, 10.0] {
            let omega = ratio * crossover;
            // Independent owner: the existing real-caloric general-fluid
            // kernel supplies the full finite-wall response at each endpoint.
            let equilibrium_reference = law.effective_fluid(gas, omega).unwrap();
            let frozen_reference = law.effective_fluid(&frozen_fluid, omega).unwrap();
            let equilibrium_response = law.effective_relaxing(&equilibrated, omega).unwrap();
            let frozen_response = law.effective_relaxing(&frozen, omega).unwrap();
            relative(
                equilibrium_response.dynamic_bulk_modulus_pa(),
                equilibrium_reference.dynamic_bulk_modulus_pa(),
                2.0e-12,
            );
            relative(
                frozen_response.dynamic_bulk_modulus_pa(),
                frozen_reference.dynamic_bulk_modulus_pa(),
                2.0e-12,
            );
            assert_eq!(
                frozen_response.dynamic_density_kg_m3(),
                frozen_reference.dynamic_density_kg_m3()
            );
            if ratio == 1.0 {
                assert!(
                    (frozen_reference.dynamic_bulk_modulus_pa()
                        - equilibrium_reference.dynamic_bulk_modulus_pa())
                    .abs()
                        > 1.0e-3 * equilibrium_reference.dynamic_bulk_modulus_pa().abs()
                );
            }
        }
    }
}

#[test]
fn coupled_reservoirs_preserve_passivity_and_the_same_viscous_drag() {
    let relaxation = relaxing_gas(&[(1000.0, 2.0e-4), (2000.0, 5.0e-3)]);
    for lafarge in [false, true] {
        let law = pores(lafarge, 1.0);
        for exponent in -6..=14 {
            let omega = 10.0_f64.powi(exponent);
            let response = law.effective_relaxing(&relaxation, omega).unwrap();
            let original = law.effective_fluid(relaxation.gas(), omega).unwrap();
            assert_eq!(
                response.dynamic_density_kg_m3(),
                original.dynamic_density_kg_m3()
            );
            assert_eq!(
                response.flow_resistivity_pa_s_m2(),
                original.flow_resistivity_pa_s_m2()
            );
            let bulk = response.dynamic_bulk_modulus_pa();
            let (impedance, wavenumber) = response.characteristic();
            assert!(bulk.re > 0.0 && bulk.im >= 0.0);
            assert!(bulk.recip().im <= 0.0);
            assert!(impedance.re > 0.0 && wavenumber.re > 0.0 && wavenumber.im <= 0.0);
            for thickness in [1.0e-6, 0.05, 100.0] {
                assert!(response.surface_impedance_rigid(thickness).unwrap().re >= 0.0);
            }
        }
        let coupled = law.effective_relaxing(&relaxation, 5000.0).unwrap();
        let original = law.effective_fluid(relaxation.gas(), 5000.0).unwrap();
        assert!(
            (coupled.dynamic_bulk_modulus_pa() - original.dynamic_bulk_modulus_pa()).abs()
                > 1.0e-4 * original.dynamic_bulk_modulus_pa().abs()
        );
    }
}

#[test]
fn relaxing_layers_retain_sources_and_balance_oblique_and_evanescent_power() {
    let relaxation = relaxing_gas(&[(1000.0, 2.0e-4)]);
    let gas = relaxation.gas();
    let receiver =
        FluidAcousticState::try_new(0.8 * gas.density, 1.4 * gas.sound_speed).unwrap();
    for lafarge in [false, true] {
        let law = pores(lafarge, 1.0);
        let layers = [
            Layer::RelaxingPorous {
                law: &law,
                relaxation: &relaxation,
                thickness_m: 0.03,
            },
            Layer::FluidGap { thickness_m: 0.02 },
            Layer::RelaxingPorous {
                law: &law,
                relaxation: &relaxation,
                thickness_m: 0.015,
            },
        ];
        let stack = TransmittingLayers::new(gas, &receiver, 5000.0, &layers).unwrap();
        let prepared = stack.prepared_layers();
        assert_eq!(
            prepared.material_response(0),
            Some(&law.effective_relaxing(&relaxation, 5000.0).unwrap())
        );
        match prepared.layer(0).unwrap() {
            Layer::RelaxingPorous { relaxation: retained, .. } => {
                assert!(core::ptr::eq(retained, &relaxation));
            }
            _ => panic!("the original molecular source was discarded"),
        }
        for cosine in [1.0, 0.85, 0.5, 0.05] {
            for account in [
                prepared
                    .layer_absorption(
                        cosine, 3.0, LayerAbsorptionBudget { max_layers: 3 }, || false,
                    )
                    .unwrap(),
                stack
                    .layer_absorption(
                        cosine, 3.0, LayerAbsorptionBudget { max_layers: 3 }, || false,
                    )
                    .unwrap(),
            ] {
                assert!(account.absorbed_fraction > 0.0);
                assert!((account.reflected_fraction + account.transmitted_fraction
                    + account.absorbed_fraction - 1.0).abs() < 2.0e-11);
                assert!(account.balance_residual_w.abs() < 6.0e-11);
                assert_eq!(account.layers()[1].absorbed_fraction(), 0.0);
                for index in [0, 2] {
                    assert!(account.layers()[index].viscous_absorbed_fraction >= 0.0);
                    assert!(account.layers()[index].thermal_absorbed_fraction > 0.0);
                }
                if cosine < 0.6 {
                    assert_eq!(account.transmitted_fraction, 0.0);
                }
            }
        }
        let different =
            GasState::try_new(&GasSpec::dry_air_ussa1976(), 310.0, gas.pressure).unwrap();
        assert!(RigidBackedLayers::new(&different, 5000.0, &layers).is_err());
        assert!(RigidBackedLayers::new(gas, 1.0e-12, &layers).is_err());
    }
}

#[test]
fn local_molecular_and_pore_heat_equals_intensity_loss_and_layer_absorption() {
    let relaxation = relaxing_gas(&[(1000.0, 2.0e-4)]);
    let gas = relaxation.gas();
    let law = pores(false, 1.0);
    let thickness = 0.045;
    let stack = TransmittingLayers::new(
        gas,
        gas,
        5000.0,
        &[Layer::RelaxingPorous {
            law: &law,
            relaxation: &relaxation,
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let spectrum = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.6, 0.0, 0.8],
                pressure_pa: C64::new(2.0, 1.0),
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let interior = spectrum.prepare_layer_field(0, &mut || false).unwrap();
    let sample = |depth| {
        interior.field_with_heat([0.0, 0.0, depth], &mut || false).unwrap()
    };
    let middle = sample(0.01);
    let heat = middle.viscous_heating_w_m3 + middle.thermal_heating_w_m3;
    let h = 1.0e-7;
    let divergence = (sample(0.01 - h).field.active_intensity_w_m2[2]
        - sample(0.01 + h).field.active_intensity_w_m2[2])
        / (2.0 * h);
    assert!(middle.thermal_heating_w_m3 > 0.0);
    assert!((divergence - heat).abs() < 1.0e-6 * heat);

    // Independent spatial quadrature of local loss, without using a power
    // residual or the production analytic layer integral as the heat source.
    let intervals = 256;
    let dz = thickness / f64::from(intervals);
    let mut integral = 0.0;
    for index in 0..=intervals {
        let q = sample(f64::from(index) * dz);
        let weight = if index == 0 || index == intervals {
            1.0
        } else if index % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * (q.viscous_heating_w_m3 + q.thermal_heating_w_m3);
    }
    integral *= dz / 3.0;
    let incident_flux = 5.0 * 0.8 / (2.0 * gas.characteristic_impedance);
    let account = stack
        .layer_absorption(
            0.8, incident_flux, LayerAbsorptionBudget { max_layers: 1 }, || false,
        )
        .unwrap();
    assert!((integral - account.absorbed_power_w).abs()
        < 1.0e-7 * account.absorbed_power_w);
}

// Independent conventional matrix for a moderate free-gas layer. The bulk
// modulus is reconstructed from the original reservoir sum, without using the
// production caloric response, equivalent-fluid sample or layer propagator.
fn free_layer_reference(
    relaxation: &RelaxationState,
    incident: &GasState,
    receiver: &GasState,
    omega: f64,
    thickness: f64,
    cosine: f64,
) -> (C64, C64, C64, C64) {
    let gas = relaxation.gas();
    let cv = relaxation.processes().iter().fold(
        C64::from_re(relaxation.frozen_cv_j_kg_k()),
        |sum, process| {
            sum + C64::from_re(process.heat_capacity_j_kg_k)
                / C64::new(1.0, omega * process.relaxation_time_s)
        },
    );
    let bulk = ((cv + C64::from_re(gas.specific_gas_constant)) / cv).scale(gas.pressure);
    let kt = omega / incident.sound_speed * (1.0 - cosine * cosine).sqrt();
    let k = (C64::from_re(gas.density) / bulk).sqrt().scale(omega);
    let kn = (k * k - C64::from_re(kt * kt)).sqrt();
    let zn = C64::from_re(omega * gas.density) / kn;
    let zn1 = incident.characteristic_impedance / cosine;
    let kn3 = ((omega / receiver.sound_speed).powi(2) - kt * kt).sqrt();
    let zn3 = omega * receiver.density / kn3;
    let phase = kn.scale(thickness);
    let sin = C64::new(
        phase.re.sin() * phase.im.cosh(),
        phase.re.cos() * phase.im.sinh(),
    );
    let cos = C64::new(
        phase.re.cos() * phase.im.cosh(),
        -phase.re.sin() * phase.im.sinh(),
    );
    let i_sin = C64::new(-sin.im, sin.re);
    let pressure = cos + i_sin * zn.scale(1.0 / zn3);
    let velocity = i_sin * C64::from_re(zn1) / zn + cos.scale(zn1 / zn3);
    (
        (pressure - velocity) / (pressure + velocity),
        C64::from_re(2.0) / (pressure + velocity),
        bulk,
        kn,
    )
}

#[test]
fn free_relaxing_layer_refracts_with_original_calorics_and_no_fictitious_drag() {
    let relaxation = relaxing_gas(&[(1000.0, 2.0e-4)]);
    let incident =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 270.0, 101_325.0).unwrap();
    let receiver =
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 290.0, 120_000.0).unwrap();
    let thickness = 0.035;
    for omega in [50.0, 5000.0, 5.0e5] {
        let stack = TransmittingLayers::new(
            &incident,
            &receiver,
            omega,
            &[Layer::RelaxingFluid {
                relaxation: &relaxation,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let response = stack.prepared_layers().material_response(0).unwrap();
        assert_eq!(response.flow_resistivity_pa_s_m2(), 0.0);
        assert_eq!(
            response.dynamic_density_kg_m3(),
            C64::from_re(relaxation.gas().density)
        );
        relative(
            response.characteristic().1,
            relaxation.relaxing_wavenumber(omega).unwrap().conj(),
            4e-15,
        );
        assert!(matches!(stack.prepared_layers().layer(0),
            Some(Layer::RelaxingFluid { relaxation: original, .. }) if core::ptr::eq(original, &relaxation)));
        for cosine in [1.0, 0.8, 0.4] {
            let actual = stack.scattering(cosine).unwrap();
            let (reflection, transmission, bulk, _) = free_layer_reference(
                &relaxation,
                &incident,
                &receiver,
                omega,
                thickness,
                cosine,
            );
            relative(actual.pressure_reflection, reflection, 2e-12);
            relative(actual.pressure_transmission, transmission, 2e-12);
            relative(response.dynamic_bulk_modulus_pa(), bulk, 3e-15);
            let heat = stack
                .layer_absorption(
                    cosine,
                    3.0,
                    LayerAbsorptionBudget { max_layers: 1 },
                    || false,
                )
                .unwrap();
            assert_eq!(heat.layers()[0].viscous_absorbed_fraction, 0.0);
            assert!(heat.layers()[0].thermal_absorbed_fraction > 0.0);
            assert!(
                (heat.reflected_fraction + heat.transmitted_fraction + heat.absorbed_fraction - 1.0)
                    .abs()
                    < 2e-13
            );
        }
    }
}

#[test]
fn free_relaxing_layer_interior_heat_closes_flux_and_independent_spatial_integral() {
    let relaxation = relaxing_gas(&[(1000.0, 2.0e-4)]);
    let gas = relaxation.gas();
    let omega = 5000.0;
    let thickness = 0.045;
    let cosine = 0.8;
    let stack = TransmittingLayers::new(
        gas,
        gas,
        omega,
        &[Layer::RelaxingFluid {
            relaxation: &relaxation,
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let spectrum = stack
        .prepare_angular_spectrum(
            &[IncidentPlaneWave {
                direction: [0.6, 0.0, cosine],
                pressure_pa: C64::ONE,
            }],
            AngularSpectrumBudget { max_components: 1 },
            &mut || false,
        )
        .unwrap();
    let interior = spectrum.prepare_layer_field(0, &mut || false).unwrap();
    let (_, transmitted, bulk, kn) =
        free_layer_reference(&relaxation, gas, gas, omega, thickness, cosine);
    let zn = C64::from_re(omega * gas.density) / kn;
    let reference = gas.characteristic_impedance / cosine;
    let intervals = 1024;
    let dz = thickness / f64::from(intervals);
    let mut integrated_heat = 0.0;
    for index in 0..=intervals {
        let depth = f64::from(index) * dz;
        let phase = kn.scale(thickness - depth);
        let sin = C64::new(
            phase.re.sin() * phase.im.cosh(),
            phase.re.cos() * phase.im.sinh(),
        );
        let cos = C64::new(
            phase.re.cos() * phase.im.cosh(),
            -phase.re.sin() * phase.im.sinh(),
        );
        let pressure =
            transmitted * (cos + C64::new(-sin.im, sin.re) * zn.scale(1.0 / reference));
        let expected_heat = -0.5 * omega * bulk.recip().im * pressure.norm_sq();
        let actual = interior
            .field_with_heat([0.0, 0.0, depth], &mut || false)
            .unwrap();
        relative(actual.field.pressure_pa, pressure, 4e-14);
        assert_eq!(actual.viscous_heating_w_m3, 0.0);
        assert!((actual.thermal_heating_w_m3 / expected_heat - 1.0).abs() < 8e-14);
        let weight = if index == 0 || index == intervals {
            1.0
        } else if index % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integrated_heat += weight * expected_heat;
    }
    integrated_heat *= dz / 3.0;
    let entrance = interior.field([0.0, 0.0, 0.0], &mut || false).unwrap();
    let exit = interior
        .field([0.0, 0.0, thickness], &mut || false)
        .unwrap();
    let flux_drop = entrance.active_intensity_w_m2[2] - exit.active_intensity_w_m2[2];
    let account = stack
        .layer_absorption(
            cosine,
            1.0,
            LayerAbsorptionBudget { max_layers: 1 },
            || false,
        )
        .unwrap();
    assert!((flux_drop / integrated_heat - 1.0).abs() < 3e-12);
    assert!(
        (account.absorbed_fraction / (2.0 * reference * integrated_heat) - 1.0).abs() < 3e-12
    );
}

#[test]
fn absent_free_reservoirs_preserve_lossless_refraction_and_critical_incidence() {
    let absent = relaxing_gas(&[]);
    let zero = relaxing_gas(&[(1e8, 2e-4)]);
    let incident = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1.0, 0.8 * absent.gas().sound_speed).unwrap(),
        FluidTransportState::try_new(1.0, 1000.0, 0.02, 1e-5).unwrap(),
        1.4,
    )
    .unwrap();
    for relaxation in [&absent, &zero] {
        let gas = relaxation.gas();
        let original = TransmittingLayers::new(
            &incident,
            &incident,
            5000.0,
            &[Layer::FluidLayer {
                fluid: gas,
                thickness_m: 0.04,
            }],
        )
        .unwrap();
        let new = TransmittingLayers::new(
            &incident,
            &incident,
            5000.0,
            &[Layer::RelaxingFluid {
                relaxation,
                thickness_m: 0.04,
            }],
        )
        .unwrap();
        assert!(new.prepared_layers().material_response(0).is_none());
        for cosine in [1.0, 0.8, 0.6, 0.4, 1e-8] {
            let a = new.scattering(cosine).unwrap();
            let b = original.scattering(cosine).unwrap();
            assert_eq!(a.pressure_reflection, b.pressure_reflection);
            assert_eq!(a.pressure_transmission, b.pressure_transmission);
            assert_eq!(a.absorbed_fraction, 0.0);
        }
    }
}

#[test]
fn free_relaxing_layer_split_and_thick_limits_preserve_passive_transport() {
    let relaxation = relaxing_gas(&[(1000.0, 2e-4)]);
    let gas = relaxation.gas();
    for thickness in [0.0, 1e-9, 0.04, 1e4] {
        let whole = TransmittingLayers::new(
            gas,
            gas,
            5000.0,
            &[Layer::RelaxingFluid {
                relaxation: &relaxation,
                thickness_m: thickness,
            }],
        )
        .unwrap();
        let split = TransmittingLayers::new(
            gas,
            gas,
            5000.0,
            &[Layer::RelaxingFluid {
                relaxation: &relaxation,
                thickness_m: thickness * 0.5,
            }; 2],
        )
        .unwrap();
        for cosine in [1.0, 0.5, 0.1] {
            let a = whole.scattering(cosine).unwrap();
            let b = split.scattering(cosine).unwrap();
            assert!((a.pressure_reflection - b.pressure_reflection).abs() < 2e-13);
            assert!((a.pressure_transmission - b.pressure_transmission).abs() < 2e-13);
            let heat = whole
                .layer_absorption(
                    cosine,
                    1.0,
                    LayerAbsorptionBudget { max_layers: 1 },
                    || false,
                )
                .unwrap();
            assert_eq!(heat.layers()[0].viscous_absorbed_fraction, 0.0);
            if thickness == 0.0 {
                assert_eq!(heat.absorbed_fraction, 0.0);
                assert_eq!(a.transmitted_fraction, 1.0);
            } else {
                assert!(heat.absorbed_fraction > 0.0);
            }
            if thickness == 1e4 {
                assert_eq!(a.transmitted_fraction, 0.0);
            }
        }
        assert!(
            whole
                .layer_absorption(
                    0.8,
                    1.0,
                    LayerAbsorptionBudget { max_layers: 1 },
                    || true,
                )
                .is_err()
        );
    }
}
