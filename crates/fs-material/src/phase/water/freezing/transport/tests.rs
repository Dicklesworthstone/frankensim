use super::*;
use crate::phase::water::freezing::FreezingWaterControls;
use fs_thermochem::water::ice::IceMeltingSolve;

fn chart(p: f64, minimum: f64, maximum: f64) -> FreezingWater {
    FreezingWater::new(
        Pressure::new(p),
        FreezingWaterControls {
            minimum_temperature_k: minimum,
            maximum_temperature_k: maximum,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            maximum_gibbs_mismatch_j_kg: 30.0,
            melting: IceMeltingSolve {
                pressure_tolerance_pa: 1e-5,
                temperature_tolerance_k: 1e-10,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap()
}

#[test]
fn pure_phase_transport_retains_original_calorics_without_absent_phase_extrapolation() {
    let chart = chart(NORMAL_PRESSURE_PA, 243.15, 303.15);
    let transport =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Series, || false).unwrap();
    for t in [243.15, 253.15, 273.15, 280.0, 303.15] {
        let original = chart
            .state_at_temperature(Temperature::new(t), || false)
            .unwrap();
        let state = transport
            .state_at_material_state(original, || false)
            .unwrap();
        assert_eq!(*state.material_state(), original);
        let h = original.specific_enthalpy_j_kg();
        let inverse = transport.state_at_specific_enthalpy(h, || false).unwrap();
        assert!((inverse.conductivity_w_m_k() - state.conductivity_w_m_k()).abs() < 1e-11);
        assert!((inverse.material_state().specific_enthalpy_j_kg() - h).abs() < 1e-7);
        assert_eq!(
            transport.conductivity_w_m_k(h).unwrap(),
            inverse.conductivity_w_m_k()
        );
        if original.liquid_mass_fraction() == 0.0 {
            assert!(state.liquid_source().is_none());
            assert_eq!(
                *state.ice_source().unwrap(),
                ice_at_temperature(Temperature::new(t)).unwrap()
            );
        } else {
            assert!(state.ice_source().is_none());
            assert_eq!(
                *state.liquid_source().unwrap(),
                liquid_at_pressure_temperature(
                    Pressure::new(NORMAL_PRESSURE_PA),
                    Temperature::new(t)
                )
                .unwrap()
            );
        }
    }
}

#[test]
fn fusion_conductivity_uses_physical_phase_volumes_and_declared_layer_orientation() {
    let chart = chart(NORMAL_PRESSURE_PA, 243.15, 303.15);
    let series =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Series, || false).unwrap();
    let parallel =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Parallel, || false)
            .unwrap();
    let t = chart.melting().temperature();
    let ki = ice_at_temperature(t).unwrap().conductivity().value();
    let kl = liquid_at_pressure_temperature(Pressure::new(NORMAL_PRESSURE_PA), t)
        .unwrap()
        .conductivity()
        .value();
    for x in [0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0] {
        let material = chart.state_at_liquid_fraction(x, || false).unwrap();
        let a = series.state_at_material_state(material, || false).unwrap();
        let b = parallel
            .state_at_material_state(material, || false)
            .unwrap();
        assert_eq!(a.material_state(), b.material_state());
        assert_eq!(a.material_state().temperature_k(), t.value());
        let ice_volume = (1.0 - x) * chart.melting().ice().specific_volume().value();
        let liquid_volume = x * chart.melting().liquid().specific_volume().value();
        let total_volume = ice_volume + liquid_volume;
        // Unit-area layers: length = volume, R_total = sum(length/k).
        let layer_k = total_volume / (ice_volume / ki + liquid_volume / kl);
        let parallel_k = (ki * ice_volume + kl * liquid_volume) / total_volume;
        assert!((a.conductivity_w_m_k() - layer_k).abs() < 2e-15);
        assert!((b.conductivity_w_m_k() - parallel_k).abs() < 2e-15);
        assert!(a.conductivity_w_m_k() >= kl - 2e-15);
        assert!(b.conductivity_w_m_k() <= ki + 2e-15);
        assert!(a.conductivity_w_m_k() <= b.conductivity_w_m_k() + 2e-15);
        if x == 0.0 {
            assert_eq!(a.conductivity_w_m_k(), ki);
            assert!(a.liquid_source().is_none());
        } else if x == 1.0 {
            assert_eq!(a.conductivity_w_m_k(), kl);
            assert!(a.ice_source().is_none());
        } else {
            assert!(a.ice_source().is_some() && a.liquid_source().is_some());
        }
        if x == 0.5 {
            assert!((material.liquid_volume_fraction() - x).abs() > 0.02);
            let incorrect_mass_average = 1.0 / ((1.0 - x) / ki + x / kl);
            assert!((a.conductivity_w_m_k() - incorrect_mass_average).abs() > 0.02);
        }
    }
}

#[test]
fn transport_domain_binding_and_cancellation_refuse_without_altering_source_state() {
    let original = chart(NORMAL_PRESSURE_PA, 243.15, 303.15);
    let law = FreezingWaterConductivity::new(original, PhaseConductivityMixture::Series, || false)
        .unwrap();
    let other = chart(NORMAL_PRESSURE_PA, 243.15, 313.15);
    let other_state = other.state_at_liquid_fraction(0.5, || false).unwrap();
    assert_eq!(
        law.state_at_material_state(other_state, || false),
        Err(FreezingConductivityError::ChartMismatch)
    );
    for p in [10_000.0, 100_000.0, 110_000.0] {
        assert!(matches!(
            FreezingWaterConductivity::new(
                chart(p, 243.15, 303.15),
                PhaseConductivityMixture::Series,
                || false
            ),
            Err(FreezingConductivityError::PressureOutsideModel { .. })
        ));
    }
    assert!(matches!(
        FreezingWaterConductivity::new(
            chart(NORMAL_PRESSURE_PA, 233.15, 303.15),
            PhaseConductivityMixture::Series,
            || false
        ),
        Err(FreezingConductivityError::Source(_))
    ));
    let state = original.state_at_liquid_fraction(0.5, || false).unwrap();
    assert_eq!(
        law.state_at_material_state(state, || true),
        Err(FreezingConductivityError::Cancelled)
    );
    let mut checks = 0;
    assert_eq!(
        law.state_at_material_state(state, || {
            checks += 1;
            checks == 2
        }),
        Err(FreezingConductivityError::Cancelled)
    );
    assert_eq!(
        law.state_at_specific_enthalpy(state.specific_enthalpy_j_kg(), || true),
        Err(FreezingConductivityError::Cancelled)
    );
    assert_eq!(
        FreezingWaterConductivity::new(original, PhaseConductivityMixture::Series, || true),
        Err(FreezingConductivityError::Cancelled)
    );
    for h in [
        original.enthalpy_bounds_j_kg()[0] - 1.0,
        original.enthalpy_bounds_j_kg()[1] + 1.0,
        f64::NAN,
    ] {
        assert!(law.state_at_specific_enthalpy(h, || false).is_err());
        assert_eq!(
            law.conductivity_w_m_k(h),
            Err(DepositionError::MaterialDomain)
        );
    }
    assert_eq!(law.chart(), &original);
    assert_eq!(state.liquid_mass_fraction(), 0.5);
}

#[test]
fn sensible_enthalpy_tangents_follow_original_cp_and_transport_derivatives() {
    let chart = chart(NORMAL_PRESSURE_PA, 243.15, 303.15);
    let transport =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Series, || false).unwrap();
    for t in [244.0, 253.15, 270.0, 274.0, 290.0, 302.0] {
        let at_t = |t| {
            let material = chart
                .state_at_temperature(Temperature::new(t), || false)
                .unwrap();
            transport
                .state_at_material_state(material, || false)
                .unwrap()
        };
        let state = at_t(t);
        assert_eq!(state.temperature_k(), t);
        let cp = state
            .material_state()
            .ice_source()
            .map(|source| source.cp().value())
            .or_else(|| {
                state
                    .material_state()
                    .liquid_source()
                    .map(|source| source.cp().value())
            })
            .unwrap();
        let dtdh = state.temperature_enthalpy_derivative_kg_k_j();
        let dkdh = state.conductivity_enthalpy_derivative_kg_w_m_k_j();
        assert_eq!(dtdh, 1.0 / cp);
        let transport_derivative = state
            .ice_source()
            .or(state.liquid_source())
            .unwrap()
            .temperature_derivative()
            .value();
        assert_eq!(dkdh, transport_derivative * dtdh);
        let lo = at_t(t - 0.002);
        let hi = at_t(t + 0.002);
        let dh = hi.material_state().specific_enthalpy_j_kg()
            - lo.material_state().specific_enthalpy_j_kg();
        let numerical_temperature = (hi.temperature_k() - lo.temperature_k()) / dh;
        let numerical_conductivity = (hi.conductivity_w_m_k() - lo.conductivity_w_m_k()) / dh;
        assert!((dtdh / numerical_temperature - 1.0).abs() < 2e-8);
        assert!((dkdh / numerical_conductivity - 1.0).abs() < 2e-8);
        assert!(dtdh > 0.0);
        if t < chart.melting().temperature().value() {
            assert!(dkdh < 0.0);
        } else {
            assert!(dkdh > 0.0);
        }
    }
}

#[test]
fn latent_tangents_include_phase_volume_change_and_choose_plateau_side_at_corners() {
    let chart = chart(NORMAL_PRESSURE_PA, 243.15, 303.15);
    let latent = chart.melting().latent_enthalpy_j_kg();
    for mixture in [
        PhaseConductivityMixture::Series,
        PhaseConductivityMixture::Parallel,
    ] {
        let law = FreezingWaterConductivity::new(chart, mixture, || false).unwrap();
        let at_x = |x| {
            law.state_at_material_state(
                chart.state_at_liquid_fraction(x, || false).unwrap(),
                || false,
            )
            .unwrap()
        };
        for x in [0.0, 0.1, 0.5, 0.9, 1.0] {
            let state = at_x(x);
            assert_eq!(state.temperature_enthalpy_derivative_kg_k_j(), 0.0);
            assert_eq!(state.temperature_k(), chart.melting().temperature().value());
            let delta = 1e-6;
            let (lo, hi) = if x == 0.0 {
                (state, at_x(delta))
            } else if x == 1.0 {
                (at_x(1.0 - delta), state)
            } else {
                (at_x(x - delta), at_x(x + delta))
            };
            let mass_fraction_change = hi.material_state().liquid_mass_fraction()
                - lo.material_state().liquid_mass_fraction();
            let numerical = (hi.conductivity_w_m_k() - lo.conductivity_w_m_k())
                / (latent * mass_fraction_change);
            let analytic = state.conductivity_enthalpy_derivative_kg_w_m_k_j();
            assert!(analytic < 0.0);
            assert!(
                (analytic / numerical - 1.0).abs() < 4e-6,
                "{mixture:?} at x={x}: analytic={analytic}, numerical={numerical}"
            );
        }
        // The ice endpoint selects the latent-side slope; the adjacent admitted
        // sensible ice has its original finite positive reciprocal Cp instead.
        let ice_h = chart.melting().ice().specific_enthalpy().value();
        let sensible = law
            .state_at_specific_enthalpy(ice_h - 1.0, || false)
            .unwrap();
        let endpoint = law.state_at_specific_enthalpy(ice_h, || false).unwrap();
        assert!(sensible.temperature_enthalpy_derivative_kg_k_j() > 0.0);
        assert_eq!(endpoint.temperature_enthalpy_derivative_kg_k_j(), 0.0);
        let liquid_h = chart.melting().liquid().specific_enthalpy().value();
        let endpoint = law.state_at_specific_enthalpy(liquid_h, || false).unwrap();
        let sensible = law
            .state_at_specific_enthalpy(liquid_h + 1.0, || false)
            .unwrap();
        assert_eq!(endpoint.temperature_enthalpy_derivative_kg_k_j(), 0.0);
        assert!(endpoint.conductivity_enthalpy_derivative_kg_w_m_k_j() < 0.0);
        assert!(sensible.temperature_enthalpy_derivative_kg_k_j() > 0.0);
        assert!(sensible.conductivity_enthalpy_derivative_kg_w_m_k_j() > 0.0);
    }
}

#[test]
fn implicit_binding_retains_original_objects_domain_and_fused_source_state() {
    use crate::phase::water::freezing::thermal::FreezingThermalControls;
    use fs_thermochem::water::ice::bounds::IceIhHeatCapacityBoundSolve;
    use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;

    let chart = chart(NORMAL_PRESSURE_PA, 243.15, 303.15);
    let caloric = FreezingWaterThermalLaw::new(
        chart,
        FreezingThermalControls {
            ice_heat_capacity: IceIhHeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
            liquid_heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )
    .unwrap();
    let transport =
        FreezingWaterConductivity::new(chart, PhaseConductivityMixture::Series, || false).unwrap();
    let binding = FreezingWaterImplicitLaw::new(&caloric, &transport, || false).unwrap();
    assert!(core::ptr::eq(binding.caloric(), &caloric));
    assert!(core::ptr::eq(binding.transport(), &transport));
    assert!(core::ptr::addr_eq(
        binding.caloric_law(),
        &caloric as &dyn ThermalLaw
    ));
    assert!(core::ptr::addr_eq(
        binding.conductivity_law(),
        &transport as &dyn ConductivityLaw
    ));
    assert_eq!(binding.enthalpy_bounds_j_kg(), chart.enthalpy_bounds_j_kg());
    let pair = chart.melting();
    for h in [
        chart.enthalpy_bounds_j_kg()[0],
        pair.ice().specific_enthalpy().value(),
        pair.ice().specific_enthalpy().value() + 0.37 * pair.latent_enthalpy_j_kg(),
        pair.liquid().specific_enthalpy().value(),
        chart.enthalpy_bounds_j_kg()[1],
    ] {
        let direct = transport.state_at_specific_enthalpy(h, || false).unwrap();
        let fused = binding.state_at_specific_enthalpy(h).unwrap();
        assert_eq!(fused.temperature_k, direct.temperature_k());
        assert_eq!(
            fused.temperature_slope_kg_k_j,
            direct.temperature_enthalpy_derivative_kg_k_j()
        );
        assert_eq!(fused.conductivity_w_m_k, direct.conductivity_w_m_k());
        assert_eq!(
            fused.conductivity_slope_kg_w_m_k_j,
            direct.conductivity_enthalpy_derivative_kg_w_m_k_j()
        );
        assert_eq!(fused.temperature_k, caloric.temperature_k(h).unwrap());
    }
    for h in [
        f64::NAN,
        chart.enthalpy_bounds_j_kg()[0] - 1.0,
        chart.enthalpy_bounds_j_kg()[1] + 1.0,
    ] {
        assert!(matches!(
            binding.state_at_specific_enthalpy(h),
            Err(DepositionError::MaterialDomain)
        ));
    }
    let other = FreezingWaterConductivity::new(
        super::tests::chart(NORMAL_PRESSURE_PA, 243.15, 313.15),
        PhaseConductivityMixture::Series,
        || false,
    )
    .unwrap();
    assert!(matches!(
        FreezingWaterImplicitLaw::new(&caloric, &other, || false),
        Err(FreezingConductivityError::ChartMismatch)
    ));
    assert!(matches!(
        FreezingWaterImplicitLaw::new(&caloric, &transport, || true),
        Err(FreezingConductivityError::Cancelled)
    ));
    let mut checkpoints = 0;
    assert!(matches!(
        FreezingWaterImplicitLaw::new(&caloric, &transport, || {
            checkpoints += 1;
            checkpoints == 2
        }),
        Err(FreezingConductivityError::Cancelled)
    ));
}
