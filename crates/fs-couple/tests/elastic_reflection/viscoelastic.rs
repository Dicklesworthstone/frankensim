//! Original lossy coating work and the actual mechanical backing remain distinct.

use super::close;
use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::{ElasticSupportError, OscillatorBackedViscoelasticLayers};
use fs_material::elastic::acoustic::viscoelastic::{
    ViscoelasticSlabError, ViscoelasticSlabSource,
    layers::{ViscoelasticBacking, ViscoelasticLayer, ViscoelasticLayerBudget},
};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_phs::{WallPin, wall_specific_impedance};

fn pin(resistance: f64) -> WallPin {
    // At omega=20, the native -i*omega impedance is r-7i; the coating's
    // +i*omega convention therefore requires r+7i.
    WallPin {
        surface_density: 1.0,
        stiffness_per_area: 260.0,
        resistance,
    }
}

fn budget(layers: usize) -> ViscoelasticLayerBudget {
    ViscoelasticLayerBudget {
        maximum_layers: layers,
        maximum_total_branches: 2 * layers,
    }
}

fn source<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
    density: f64,
    id: u8,
) -> ViscoelasticSlabSource<'a> {
    ViscoelasticSlabSource {
        bulk,
        shear,
        density_kg_m3: density,
        frequency_band_rad_s: (10.0, 100.0),
        maximum_absolute_strain: 0.1,
        source_identity: ContentHash([id; 32]),
    }
}

fn cc(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_absent_stack_uses_original_pin_phase_displacement_and_storage() {
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let backing = pin(13.0);
    let omega = 20.0;
    let owner = OscillatorBackedViscoelasticLayers::prepare(&[], &fluid, backing, omega, budget(0))
        .unwrap();
    assert_eq!(owner.backing(), &backing);
    assert_eq!(owner.omega_rad_s(), omega);
    assert_eq!(owner.layers().layer_count(), 0);
    assert!(core::ptr::eq(
        owner.layers().incident_source(),
        &fluid as &dyn FluidAcoustics
    ));
    let z = C64::new(13.0, 7.0);
    assert_eq!(wall_specific_impedance(&backing, omega).unwrap().conj(), z);
    assert_eq!(owner.layers().backing(), ViscoelasticBacking::Impedance(z));
    let input = C64::new(1.0, 0.2);
    let area = 0.7;
    let response = owner.drive(input, area, &mut || false).unwrap();
    let expected_r = (z - C64::from_re(20.0)) / (z + C64::from_re(20.0));
    let expected_v = input.scale(2.0) / (z + C64::from_re(20.0));
    // For exp(+i omega t), v=i*omega*x, so x=-i*v/omega.
    let expected_x = C64::new(expected_v.im / omega, -expected_v.re / omega);
    cc(
        response.layers.reflected_pressure_pa(),
        input * expected_r,
        3e-14,
    );
    cc(
        response.layers.support_normal_velocity_m_s(),
        expected_v,
        3e-15,
    );
    cc(
        response.layers.support_compressive_traction_pa(),
        z * expected_v,
        3e-14,
    );
    cc(response.support_displacement_m, expected_x, 3e-16);
    let kinetic = 0.25 * backing.surface_density * expected_v.norm_sq();
    let elastic = 0.25 * backing.stiffness_per_area * expected_x.norm_sq();
    close(response.mean_kinetic_energy_j_m2, kinetic, 3e-16);
    close(response.mean_elastic_energy_j_m2, elastic, 3e-16);
    close(
        response.mean_support_stored_energy_j_m2,
        kinetic + elastic,
        3e-16,
    );
    close(
        response.mean_support_stored_energy_j,
        area * (kinetic + elastic),
        3e-16,
    );
    let power = response.layers.power();
    close(
        power.incident_power_w,
        0.5 * area * input.norm_sq() / 20.0,
        3e-15,
    );
    close(
        power.backing_absorbed_power_w,
        0.5 * area * backing.resistance * expected_v.norm_sq(),
        3e-15,
    );
    assert_eq!(power.bulk_absorbed_power_w, 0.0);
    assert_eq!(power.shear_absorbed_power_w, 0.0);
    close(
        power.incident_power_w,
        power.reflected_power_w + power.backing_absorbed_power_w,
        3e-15,
    );
    assert!(response.layers.field(0, 0.0, &mut || false).is_err());
}

#[test]
fn g1_original_heterogeneous_layers_and_backing_match_independent_heat_integrals() {
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let bulk = [
        GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap(),
        GeneralizedMaxwell::new(130.0, vec![(15.0, 0.02)]).unwrap(),
    ];
    let shear = [
        GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap(),
        GeneralizedMaxwell::new(45.0, vec![(27.0, 0.03)]).unwrap(),
    ];
    let layers = [
        ViscoelasticLayer {
            source: source(&bulk[0], &shear[0], 4.0, 3),
            thickness_m: 0.2,
        },
        ViscoelasticLayer {
            source: source(&bulk[1], &shear[1], 7.0, 8),
            thickness_m: 0.3,
        },
    ];
    for count in [1, 2] {
        let owner = OscillatorBackedViscoelasticLayers::prepare(
            &layers[..count],
            &fluid,
            pin(13.0),
            20.0,
            budget(count),
        )
        .unwrap();
        assert!(core::ptr::eq(
            owner.layers().incident_source(),
            &fluid as &dyn FluidAcoustics
        ));
        let response = owner.drive(C64::ONE, 1.0, &mut || false).unwrap();
        let power = response.layers.power();
        let mut summed_bulk = 0.0;
        let mut summed_shear = 0.0;
        for i in 0..count {
            let retained = owner.layers().layer(i).unwrap();
            assert!(core::ptr::eq(retained.source.bulk, &bulk[i]));
            assert!(core::ptr::eq(retained.source.shear, &shear[i]));
            assert_eq!(
                retained.source.source_identity,
                layers[i].source.source_identity
            );
            let local = response.layers.layer_absorption(i, &mut || false).unwrap();
            assert_eq!(local.source_identity, layers[i].source.source_identity);
            assert!(local.bulk_absorbed_power_w > 0.0 && local.shear_absorbed_power_w > 0.0);
            assert!(local.maximum_absolute_strain_bound <= retained.source.maximum_absolute_strain);
            summed_bulk += local.bulk_absorbed_power_w;
            summed_shear += local.shear_absorbed_power_w;
        }
        close(power.bulk_absorbed_power_w, summed_bulk, 3e-15);
        close(power.shear_absorbed_power_w, summed_shear, 3e-15);
        let rear = response
            .layers
            .field(count - 1, layers[count - 1].thickness_m, &mut || false)
            .unwrap();
        cc(
            rear.normal_velocity_m_s,
            response.layers.support_normal_velocity_m_s(),
            3e-15,
        );
        cc(
            rear.compressive_traction_pa,
            C64::new(13.0, 7.0) * rear.normal_velocity_m_s,
            3e-14,
        );
        close(
            power.backing_absorbed_power_w,
            6.5 * rear.normal_velocity_m_s.norm_sq(),
            3e-14,
        );
        assert!(response.mean_support_stored_energy_j_m2 > 0.0);
        close(
            power.incident_power_w,
            power.reflected_power_w + summed_bulk + summed_shear + power.backing_absorbed_power_w,
            3e-14,
        );
        if count == 2 {
            // Independent 80-digit physical transfer matrices and direct positive
            // bulk/shear strain-work quadrature, without a reference-port model.
            cc(
                response.layers.reflected_pressure_pa(),
                C64::new(0.21661124142492808, -0.4698365313672166),
                3e-13,
            );
            let velocity = C64::new(-0.01860854261520848, -0.03580350559133733);
            cc(rear.normal_velocity_m_s, velocity, 3e-14);
            close(power.reflected_power_w, 0.006691669902970649, 3e-14);
            close(power.backing_absorbed_power_w, 0.010583097660791197, 3e-14);
            let expected = [
                [0.0033753594392632683, 0.0023278340960436333],
                [0.0004967436606334001, 0.0015252952402978522],
            ];
            for (i, expected) in expected.into_iter().enumerate() {
                let actual = response.layers.layer_absorption(i, &mut || false).unwrap();
                close(actual.bulk_absorbed_power_w, expected[0], 3e-14);
                close(actual.shear_absorbed_power_w, expected[1], 3e-14);
            }
            let left = response.layers.field(0, 0.2, &mut || false).unwrap();
            let right = response.layers.field(1, 0.0, &mut || false).unwrap();
            cc(
                left.compressive_traction_pa,
                right.compressive_traction_pa,
                3e-13,
            );
            cc(left.normal_velocity_m_s, right.normal_velocity_m_s, 3e-14);
            let kinetic = 0.25 * velocity.norm_sq();
            let elastic = 0.25 * 260.0 * velocity.norm_sq() / (20.0 * 20.0);
            close(response.mean_kinetic_energy_j_m2, kinetic, 3e-15);
            close(response.mean_elastic_energy_j_m2, elastic, 3e-15);
            close(
                response.mean_support_stored_energy_j,
                kinetic + elastic,
                3e-15,
            );
        }
    }
}

#[test]
fn g3_area_scales_total_work_while_reactive_backing_keeps_storage_without_heat() {
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let layer = ViscoelasticLayer {
        source: source(&bulk, &shear, 4.0, 12),
        thickness_m: 0.2,
    };
    let owner =
        OscillatorBackedViscoelasticLayers::prepare(&[layer], &fluid, pin(0.0), 20.0, budget(1))
            .unwrap();
    let input = C64::new(0.7, 0.2);
    let small = owner.drive(input, 0.25, &mut || false).unwrap();
    let large = owner.drive(input, 2.0, &mut || false).unwrap();
    assert_eq!(
        small.layers.support_normal_velocity_m_s(),
        large.layers.support_normal_velocity_m_s()
    );
    assert_eq!(small.support_displacement_m, large.support_displacement_m);
    assert_eq!(
        small.mean_support_stored_energy_j_m2,
        large.mean_support_stored_energy_j_m2
    );
    assert!(small.mean_support_stored_energy_j_m2 > 0.0);
    close(
        large.mean_support_stored_energy_j,
        8.0 * small.mean_support_stored_energy_j,
        3e-15,
    );
    let a = small.layers.power();
    let b = large.layers.power();
    assert_eq!(a.backing_absorbed_power_w, 0.0);
    assert_eq!(b.backing_absorbed_power_w, 0.0);
    assert!(a.bulk_absorbed_power_w > 0.0 && a.shear_absorbed_power_w > 0.0);
    for (small, large) in [
        (a.incident_power_w, b.incident_power_w),
        (a.reflected_power_w, b.reflected_power_w),
        (a.bulk_absorbed_power_w, b.bulk_absorbed_power_w),
        (a.shear_absorbed_power_w, b.shear_absorbed_power_w),
    ] {
        close(large, 8.0 * small, 3e-14);
    }
    // Removing constitutive loss as well recovers a unit-power reflector.
    // Neither a standing wave nor the reactive pin's stored energy is heat.
    let bulk = GeneralizedMaxwell::new(60.0, vec![]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![]).unwrap();
    let layer = ViscoelasticLayer {
        source: source(&bulk, &shear, 4.0, 13),
        thickness_m: 0.2,
    };
    let lossless =
        OscillatorBackedViscoelasticLayers::prepare(&[layer], &fluid, pin(0.0), 20.0, budget(1))
            .unwrap();
    let response = lossless.drive(input, 0.7, &mut || false).unwrap();
    let power = response.layers.power();
    assert_eq!(power.bulk_absorbed_power_w, 0.0);
    assert_eq!(power.shear_absorbed_power_w, 0.0);
    assert_eq!(power.backing_absorbed_power_w, 0.0);
    close(power.incident_power_w, power.reflected_power_w, 3e-14);
    assert!(response.mean_support_stored_energy_j > 0.0);
}

#[test]
fn g0_g4_original_source_domains_pin_admission_and_cancellation_remain_binding() {
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let original = source(&bulk, &shear, 4.0, 21);
    let layer = ViscoelasticLayer {
        source: original,
        thickness_m: 0.2,
    };
    assert!(matches!(
        OscillatorBackedViscoelasticLayers::prepare(&[layer], &fluid, pin(13.0), 5.0, budget(1)),
        Err(ElasticSupportError::Viscoelastic(
            ViscoelasticSlabError::OutsideSourceSupport
        ))
    ));
    assert!(matches!(
        OscillatorBackedViscoelasticLayers::prepare(&[layer], &fluid, pin(13.0), 20.0, budget(0)),
        Err(ElasticSupportError::Viscoelastic(
            ViscoelasticSlabError::Budget
        ))
    ));
    let mut missing_identity = layer;
    missing_identity.source.source_identity = ContentHash([0; 32]);
    assert!(
        OscillatorBackedViscoelasticLayers::prepare(
            &[missing_identity],
            &fluid,
            pin(13.0),
            20.0,
            budget(1)
        )
        .is_err()
    );
    let mut limited = layer;
    limited.source.maximum_absolute_strain = 1e-12;
    let owner =
        OscillatorBackedViscoelasticLayers::prepare(&[limited], &fluid, pin(13.0), 20.0, budget(1))
            .unwrap();
    assert!(matches!(
        owner.drive(C64::ONE, 1.0, &mut || false),
        Err(ElasticSupportError::Viscoelastic(
            ViscoelasticSlabError::OutsideSourceSupport
        ))
    ));
    assert!(owner.drive(C64::ZERO, 1.0, &mut || false).is_ok());
    assert!(matches!(
        OscillatorBackedViscoelasticLayers::prepare_controlled(
            &[layer],
            &fluid,
            pin(13.0),
            20.0,
            budget(1),
            &mut || true
        ),
        Err(ElasticSupportError::Cancelled)
    ));
    let owner =
        OscillatorBackedViscoelasticLayers::prepare(&[layer], &fluid, pin(13.0), 20.0, budget(1))
            .unwrap();
    assert!(matches!(
        owner.drive(C64::ONE, 1.0, &mut || true),
        Err(ElasticSupportError::Viscoelastic(
            ViscoelasticSlabError::Cancelled
        ))
    ));
    let good = owner.drive(C64::ONE, 1.0, &mut || false).unwrap();
    assert!(core::ptr::eq(good.layers.prepared(), owner.layers()));
    assert!(core::ptr::eq(
        owner.layers().layer(0).unwrap().source.bulk,
        &bulk
    ));
    assert!(good.layers.field(0, 0.1, &mut || false).is_ok());
    for area in [0.0, -1.0, f64::NAN] {
        assert!(owner.drive(C64::ONE, area, &mut || false).is_err());
    }
    let invalid_pin = WallPin {
        surface_density: 0.0,
        ..pin(13.0)
    };
    assert!(
        OscillatorBackedViscoelasticLayers::prepare(&[layer], &fluid, invalid_pin, 20.0, budget(1))
            .is_err()
    );
}
