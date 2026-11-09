//! G1/G3: native elastic reflection and original causal support mechanics.

use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::{ElasticSupportError, OscillatorSupportedSurface};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_math::c64::C64;
use fs_phs::{WallPin, wall_specific_impedance};

#[path = "elastic_reflection/coating.rs"]
mod coating;
#[path = "elastic_reflection/normal.rs"]
mod normal;
#[path = "elastic_reflection/oriented.rs"]
mod oriented;
#[path = "elastic_reflection/transient.rs"]
mod transient;
#[path = "elastic_reflection/viscoelastic.rs"]
mod viscoelastic;
#[path = "elastic_reflection/viscoelastic_heating.rs"]
mod viscoelastic_heating;

fn solid(scale: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        2500.0 * scale,
        2500.0 * (4000.0_f64.powi(2) - 4.0 / 3.0 * 2000.0_f64.powi(2)) * scale,
        2500.0 * 2000.0_f64.powi(2) * scale,
        ContentHash([19; 32]),
        0.002,
    )
    .unwrap()
}

fn pins(resistance: f64) -> [WallPin; 2] {
    [
        WallPin {
            surface_density: 1e5,
            stiffness_per_area: 1e11,
            resistance,
        },
        WallPin {
            surface_density: 2e4,
            stiffness_per_area: 2e10,
            resistance: 0.5 * resistance,
        },
    ]
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn g1_resonance_absorbs_matched_waves_and_retains_real_mechanical_storage() {
    let state = solid(1.0);
    let [normal, tangent] = pins(1e7);
    let source = [C64::from_re(0.001), C64::ZERO, C64::ZERO];
    for omega in [100.0, 1000.0, 10_000.0] {
        let owner = OscillatorSupportedSurface::prepare(
            &state,
            normal,
            tangent,
            omega,
            ElasticMode::P,
            1.0,
        )
        .unwrap();
        assert!(core::ptr::eq(owner.state(), &state));
        assert_eq!(*owner.normal_support(), normal);
        assert_eq!(*owner.tangential_support(), tangent);
        assert_eq!(owner.omega_rad_s(), omega);
        let response = owner.drive(source).unwrap();
        let z = C64::new(
            normal.resistance,
            normal.stiffness_per_area / omega - omega * normal.surface_density,
        );
        let expected = (z - C64::from_re(1e7)) / (z + C64::from_re(1e7));
        let r = response.reflection.reflected_velocity_m_s[0] / source[0];
        assert!((r - expected).abs() < 3e-14);
        let fraction = response.reflection.absorbed_normal_power_w_m2
            / response.reflection.total_incoming_power_w_m2;
        if omega == 1000.0 {
            close(fraction, 1.0, 2e-14);
            close(response.mean_kinetic_energy_j_m2[2], 0.025, 2e-15);
            close(response.mean_elastic_energy_j_m2[2], 0.025, 2e-15);
            close(response.mean_stored_energy_j_m2, 0.05, 3e-15);
        } else {
            assert!(fraction < 0.001);
        }
        close(
            response.reflection.total_reflected_power_w_m2
                + response.reflection.absorbed_normal_power_w_m2,
            response.reflection.total_incoming_power_w_m2,
            2e-12,
        );
    }
}

#[test]
fn g1_lossless_resonance_has_storage_without_absorption() {
    let state = solid(1.0);
    let [normal, tangent] = pins(0.0);
    let owner =
        OscillatorSupportedSurface::prepare(&state, normal, tangent, 1000.0, ElasticMode::P, 1.0)
            .unwrap();
    let response = owner
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    assert_eq!(response.reflection.absorbed_normal_power_w_m2, 0.0);
    close(
        response.reflection.reflected_velocity_m_s[0].re,
        -0.001,
        1e-17,
    );
    close(response.mean_kinetic_energy_j_m2[2], 0.1, 3e-15);
    close(response.mean_elastic_energy_j_m2[2], 0.1, 3e-15);
    close(response.support_displacement_m[2].im, 2e-6, 2e-20);
}

#[test]
fn g3_coherent_mode_conversion_preserves_each_supports_work_and_storage() {
    let state = solid(1.0);
    let [normal, tangent] = pins(3e6);
    let omega = 1700.0;
    let owner =
        OscillatorSupportedSurface::prepare(&state, normal, tangent, omega, ElasticMode::P, 0.7)
            .unwrap();
    let input = [
        C64::new(0.001, 0.0003),
        C64::new(-0.0007, 0.0002),
        C64::new(0.0004, -0.0001),
    ];
    let response = owner.drive(input).unwrap();
    let mut expected_loss = 0.0;
    for axis in 0..3 {
        let pin = if axis == 2 { normal } else { tangent };
        let v = response.reflection.boundary_velocity_m_s[axis];
        let x = response.support_displacement_m[axis];
        let traction = response.reflection.boundary_traction_pa[axis];
        let z = wall_specific_impedance(&pin, omega).unwrap();
        assert!((traction + z * v).abs() < 1e-8);
        assert!((x * C64::new(0.0, -omega) - v).abs() < 1e-17);
        let kinetic = 0.25 * pin.surface_density * v.norm_sq();
        let elastic = 0.25 * pin.stiffness_per_area * x.norm_sq();
        close(response.mean_kinetic_energy_j_m2[axis], kinetic, 2e-15);
        close(response.mean_elastic_energy_j_m2[axis], elastic, 2e-15);
        // The reactive traction stores energy instead of becoming heat.
        close(
            0.5 * (traction * v.conj()).im,
            2.0 * omega * (kinetic - elastic),
            2e-12,
        );
        expected_loss += 0.5 * pin.resistance * v.norm_sq();
    }
    close(
        response.reflection.absorbed_normal_power_w_m2,
        expected_loss,
        2e-13,
    );
    let replay = owner.drive(input).unwrap();
    assert_eq!(
        replay.mean_stored_energy_j_m2,
        response.mean_stored_energy_j_m2
    );
    assert_eq!(
        replay.reflection.reflected_velocity_m_s,
        response.reflection.reflected_velocity_m_s
    );
}

#[test]
fn g3_original_solid_impedance_and_pin_mass_change_the_reflected_wave() {
    let first = solid(1.0);
    let second = solid(2.0);
    let [normal, tangent] = pins(1e7);
    let input = [C64::from_re(0.001), C64::ZERO, C64::ZERO];
    let original =
        OscillatorSupportedSurface::prepare(&first, normal, tangent, 1000.0, ElasticMode::P, 1.0)
            .unwrap()
            .drive(input)
            .unwrap();
    let changed =
        OscillatorSupportedSurface::prepare(&second, normal, tangent, 1000.0, ElasticMode::P, 1.0)
            .unwrap()
            .drive(input)
            .unwrap();
    assert!(original.reflection.total_reflected_power_w_m2 < 1e-25);
    assert!(changed.reflection.total_reflected_power_w_m2 > 1.0);
    let heavier = WallPin {
        surface_density: 2.0 * normal.surface_density,
        ..normal
    };
    let shifted =
        OscillatorSupportedSurface::prepare(&first, heavier, tangent, 1000.0, ElasticMode::P, 1.0)
            .unwrap()
            .drive(input)
            .unwrap();
    assert!(shifted.reflection.total_reflected_power_w_m2 > 4.0);
}

#[test]
fn g0_g4_original_pin_domains_and_cancellation_refuse_transactionally() {
    let state = solid(1.0);
    let [normal, tangent] = pins(1e7);
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            OscillatorSupportedSurface::prepare(
                &state,
                normal,
                tangent,
                omega,
                ElasticMode::P,
                1.0
            )
            .is_err()
        );
    }
    for invalid in [
        WallPin {
            resistance: -1.0,
            ..normal
        },
        WallPin {
            surface_density: 0.0,
            ..normal
        },
        WallPin {
            stiffness_per_area: f64::INFINITY,
            ..normal
        },
    ] {
        assert!(
            OscillatorSupportedSurface::prepare(
                &state,
                invalid,
                tangent,
                1000.0,
                ElasticMode::P,
                1.0
            )
            .is_err()
        );
    }
    assert!(matches!(
        OscillatorSupportedSurface::prepare_controlled(
            &state,
            normal,
            tangent,
            1000.0,
            ElasticMode::P,
            1.0,
            &mut || true
        ),
        Err(ElasticSupportError::Cancelled)
    ));
    let owner =
        OscillatorSupportedSurface::prepare(&state, normal, tangent, 1000.0, ElasticMode::P, 1.0)
            .unwrap();
    assert_eq!(
        owner.drive([C64::ZERO; 3]).unwrap().mean_stored_energy_j_m2,
        0.0
    );
    assert!(
        owner
            .drive([C64::from_re(1e9), C64::ZERO, C64::ZERO])
            .is_err()
    );
}
