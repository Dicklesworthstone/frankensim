use super::*;
use crate::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::regions::If97FluidRegion;

fn chart(pressure_pa: f64, lo: f64, hi: f64) -> RegionalWater {
    RegionalWater::new(
        Pressure::new(pressure_pa),
        RegionalWaterControls {
            minimum_temperature_k: lo,
            maximum_temperature_k: hi,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite()
            && (actual - expected).abs() <= relative * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g3_homogeneous_regions_keep_original_roots_for_both_coexistence_policies() {
    for (p, lo, hi, t, region) in [
        (1e5, 300.0, 500.0, 320.0, If97FluidRegion::Region1),
        (1e5, 300.0, 500.0, 450.0, If97FluidRegion::Region2),
        (20e6, 625.0, 645.0, 630.0, If97FluidRegion::Region3),
        (20e6, 625.0, 645.0, 642.0, If97FluidRegion::Region3),
        (25e6, 640.0, 660.0, 650.0, If97FluidRegion::Region3),
    ] {
        let law = chart(p, lo, hi);
        let thermal = law
            .state_at_temperature(Temperature::new(t), || false)
            .unwrap();
        let source = thermal.single_phase_source().unwrap();
        assert_eq!(source.region(), region);
        let native = source.acoustic_state().unwrap();
        let mut previous = None;
        for closure in [
            CoexistenceAcousticClosure::FrozenPhaseMasses,
            CoexistenceAcousticClosure::HomogeneousEquilibrium,
        ] {
            let view = RegionalWaterAcoustics::try_new(thermal, closure).unwrap();
            assert_eq!(*view.source(), thermal);
            assert_eq!(view.coexistence_closure(), closure);
            assert_eq!(
                view.response(),
                &RegionalWaterAcousticResponse::SinglePhase(source)
            );
            assert_eq!(
                view.state().density_kg_m3().to_bits(),
                thermal.bulk_density_kg_m3().to_bits()
            );
            assert_eq!(
                view.state().adiabatic_bulk_modulus_pa(),
                native.adiabatic_bulk_modulus_pa()
            );
            near(view.state().sound_speed_m_s(), native.sound_speed_m_s(), 4e-15);
            if native.density_kg_m3().to_bits() == thermal.bulk_density_kg_m3().to_bits() {
                assert_eq!(view.state(), native);
            }
            if let Some(before) = previous {
                assert_eq!(view.state(), before);
            }
            previous = Some(view.state());
            let consumer: &dyn FluidAcoustics = &view;
            assert_eq!(consumer.acoustic_state().unwrap(), view.state());
        }
    }
}

#[test]
fn g1_latent_heating_reaches_acoustics_with_one_density_and_an_explicit_closure() {
    for (p, lo, hi) in [(1e5, 300.0, 500.0), (20e6, 625.0, 645.0)] {
        let law = chart(p, lo, hi);
        let initial = law.state_at_quality(0.2, || false).unwrap();
        let target = law.state_at_quality(0.6, || false).unwrap();
        let mass = 0.3;
        let heat = mass * (target.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
        let step = law.advance_heat(initial, mass, heat, 2e-5, || false).unwrap();
        let thermal = step.state;
        assert_eq!(thermal.temperature_k(), initial.temperature_k());
        assert!(thermal.bulk_density_kg_m3() < initial.bulk_density_kg_m3());
        near(thermal.vapor_quality().unwrap(), 0.6, 2e-14);
        let frozen = RegionalWaterAcoustics::try_new(
            thermal,
            CoexistenceAcousticClosure::FrozenPhaseMasses,
        )
        .unwrap();
        let equilibrium = RegionalWaterAcoustics::try_new(
            thermal,
            CoexistenceAcousticClosure::HomogeneousEquilibrium,
        )
        .unwrap();
        for view in [&frozen, &equilibrium] {
            assert_eq!(*view.source(), thermal);
            assert_eq!(view.source().law(), &law);
            assert_eq!(
                view.state().density_kg_m3().to_bits(),
                thermal.bulk_density_kg_m3().to_bits()
            );
            near(view.state().density_kg_m3(), mass / step.volumes_m3[1], 2e-14);
            near(
                view.state().sound_speed_m_s().powi(2),
                view.state().adiabatic_bulk_modulus_pa() / thermal.bulk_density_kg_m3(),
                2e-14,
            );
            near(
                view.state().characteristic_impedance_pa_s_m(),
                thermal.bulk_density_kg_m3() * view.state().sound_speed_m_s(),
                2e-14,
            );
        }
        let pair = *thermal.law().saturation().unwrap();
        let x = thermal.vapor_quality().unwrap();
        let a = pair.liquid().acoustic_state().unwrap();
        let b = pair.vapor().acoustic_state().unwrap();
        // Independent Wood relation from phase volumes per kg, never alpha=x.
        let vl = (1.0 - x) / a.density_kg_m3();
        let vv = x / b.density_kg_m3();
        let compliance = (vl * a.adiabatic_compressibility_per_pa()
            + vv * b.adiabatic_compressibility_per_pa())
            / (vl + vv);
        near(frozen.state().adiabatic_bulk_modulus_pa(), 1.0 / compliance, 2e-14);
        let RegionalWaterAcousticResponse::FrozenPhaseMasses(original_frozen) = frozen.response()
        else {
            panic!("an interior frozen state must retain its mixture owner");
        };
        assert_eq!(original_frozen.dispersed_mass_fraction(), x);
        assert_eq!(
            frozen.density_representation_difference_kg_m3(),
            thermal.bulk_density_kg_m3() - original_frozen.state().density_kg_m3()
        );
        let RegionalWaterAcousticResponse::HomogeneousEquilibrium(original_equilibrium) =
            equilibrium.response()
        else {
            panic!("an interior equilibrium state must retain its tangent owner");
        };
        assert_eq!(*original_equilibrium.source(), pair);
        assert_eq!(original_equilibrium.vapor_mass_fraction(), x);
        assert_eq!(
            original_equilibrium.tangent().source_phase_pressure_residuals_pa,
            thermal.source_pressure_residuals_pa()
        );
        assert_eq!(
            equilibrium.state().adiabatic_bulk_modulus_pa(),
            original_equilibrium.state().adiabatic_bulk_modulus_pa()
        );
        assert_eq!(
            equilibrium.density_representation_difference_kg_m3(),
            thermal.bulk_density_kg_m3() - original_equilibrium.state().density_kg_m3()
        );
        assert!(equilibrium.state().sound_speed_m_s() < frozen.state().sound_speed_m_s());
        assert_eq!(
            frozen.state().density_kg_m3().to_bits(),
            equilibrium.state().density_kg_m3().to_bits()
        );
    }
}

#[test]
fn g0_complete_condensation_and_boiling_have_explicit_one_sided_endpoint_semantics() {
    for (p, lo, hi) in [(1e5, 300.0, 500.0), (20e6, 625.0, 645.0)] {
        let law = chart(p, lo, hi);
        let initial = law.state_at_quality(0.25, || false).unwrap();
        for quality in [0.0, 1.0] {
            let endpoint = law.state_at_quality(quality, || false).unwrap();
            let heat = endpoint.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg();
            let thermal = law.advance_heat(initial, 1.0, heat, 2e-5, || false).unwrap().state;
            assert_eq!(thermal.vapor_quality(), Some(quality));
            let pair = *law.saturation().unwrap();
            let phase = if quality == 0.0 { pair.liquid() } else { pair.vapor() };
            let native = phase.acoustic_state().unwrap();
            let frozen = RegionalWaterAcoustics::try_new(
                thermal,
                CoexistenceAcousticClosure::FrozenPhaseMasses,
            )
            .unwrap();
            assert_eq!(*frozen.source(), thermal);
            assert_eq!(
                frozen.response(),
                &RegionalWaterAcousticResponse::SinglePhase(phase)
            );
            assert_eq!(
                frozen.state().density_kg_m3().to_bits(),
                thermal.bulk_density_kg_m3().to_bits()
            );
            assert_eq!(
                frozen.state().adiabatic_bulk_modulus_pa(),
                native.adiabatic_bulk_modulus_pa()
            );
            near(frozen.state().sound_speed_m_s(), native.sound_speed_m_s(), 4e-15);
            if native.density_kg_m3().to_bits() == thermal.bulk_density_kg_m3().to_bits() {
                assert_eq!(frozen.state(), native);
            }
            assert!(
                RegionalWaterAcoustics::try_new(
                    thermal,
                    CoexistenceAcousticClosure::HomogeneousEquilibrium,
                )
                .is_err()
            );
        }
    }
}
