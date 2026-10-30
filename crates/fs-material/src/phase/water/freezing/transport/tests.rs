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
