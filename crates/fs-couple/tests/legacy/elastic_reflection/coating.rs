//! G1/G3/G4: resolved coating interference and the original rear support work.

use super::{close, pins, solid};
use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::{
    ElasticSupportError, OrientedOscillatorSupport, OrientedOscillatorSupportedSurface,
    OscillatorBackedElasticCoating,
};
use fs_ga::{So3, So3Tangent, Vec3};
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_math::c64::C64;
use fs_phs::WallPin;

fn material(rho: f64, cp: f64, cs: f64, limit: f64) -> ElasticAcousticState {
    ElasticAcousticState::from_bulk_and_shear(
        rho,
        rho * (cp * cp - 4.0 / 3.0 * cs * cs),
        rho * cs * cs,
        ContentHash([23; 32]),
        limit,
    )
    .unwrap()
}

fn cc(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}

fn phase(x: f64) -> C64 {
    C64::new(x.cos(), x.sin())
}

fn rotation() -> So3 {
    let x = So3::exp(So3Tangent::new(Vec3::new(0.4, 0.0, 0.0))).unwrap();
    let y = So3::exp(So3Tangent::new(Vec3::new(0.0, -0.2, 0.0))).unwrap();
    let z = So3::exp(So3Tangent::new(Vec3::new(0.0, 0.0, 0.3))).unwrap();
    z.compose(y).unwrap().compose(x).unwrap()
}

fn support() -> OrientedOscillatorSupport {
    let [normal, tangent] = pins(3e6);
    OrientedOscillatorSupport::new(
        [
            tangent,
            WallPin {
                surface_density: 3e4,
                stiffness_per_area: 6e10,
                resistance: 2e6,
            },
            normal,
        ],
        rotation(),
    )
}

#[test]
fn g3_absent_coating_recovers_direct_oriented_support_exactly() {
    let state = solid(1.0);
    let film = material(1800.0, 3000.0, 1500.0, 0.001);
    let support = support();
    let frequency = 1700.0 / core::f64::consts::TAU;
    let omega = core::f64::consts::TAU * frequency;
    let input = [
        C64::new(0.001, 0.0002),
        C64::new(-0.0004, 0.0003),
        C64::new(0.0006, -0.0001),
    ];
    let direct =
        OrientedOscillatorSupportedSurface::prepare(&state, support, omega, ElasticMode::P, 0.7)
            .unwrap();
    let owner = OscillatorBackedElasticCoating::prepare(
        &state,
        &film,
        0.0,
        support,
        frequency,
        ElasticMode::P,
        0.7,
    )
    .unwrap();
    let original = direct.drive(input).unwrap();
    let coated = owner.drive(input).unwrap();
    assert!(core::ptr::eq(owner.coating().incident_state(), &state));
    assert!(core::ptr::eq(owner.coating().coating_state(), &film));
    assert_eq!(owner.support(), &support);
    assert_eq!(owner.omega_rad_s(), omega);
    assert_eq!(owner.coating().frequency_hz(), frequency);
    assert!(coated.coating.coating_faces.is_none());
    assert!(coated.coating.field_at(0.0).is_err());
    assert_eq!(
        coated.coating.reflection.reflected_velocity_m_s,
        original.reflection.reflected_velocity_m_s
    );
    assert_eq!(
        coated.principal_velocity_m_s,
        original.principal_velocity_m_s
    );
    assert_eq!(
        coated.support_displacement_m,
        original.support_displacement_m
    );
    assert_eq!(
        coated.mean_support_stored_energy_j_m2,
        original.mean_stored_energy_j_m2
    );
}

#[test]
fn g1_homogeneous_extension_preserves_all_modal_phases_and_rear_storage() {
    let state = solid(1.0);
    let support = support();
    let frequency = 350.0;
    let omega = core::f64::consts::TAU * frequency;
    let thickness = 0.4;
    let input = [
        C64::new(0.001, 0.0002),
        C64::new(-0.0004, 0.0003),
        C64::new(0.0006, -0.0001),
    ];
    let propagation = [4000.0, 2000.0, 2000.0].map(|c| phase(omega * thickness / c));
    let direct =
        OrientedOscillatorSupportedSurface::prepare(&state, support, omega, ElasticMode::P, 1.0)
            .unwrap();
    let original = direct
        .drive(core::array::from_fn(|m| input[m] * propagation[m]))
        .unwrap();
    let owner = OscillatorBackedElasticCoating::prepare(
        &state,
        &state,
        thickness,
        support,
        frequency,
        ElasticMode::P,
        1.0,
    )
    .unwrap();
    let response = owner.drive(input).unwrap();
    for m in 0..3 {
        cc(
            response.coating.reflection.reflected_velocity_m_s[m],
            original.reflection.reflected_velocity_m_s[m] * propagation[m],
            3e-16,
        );
        cc(
            response.principal_velocity_m_s[m],
            original.principal_velocity_m_s[m],
            3e-16,
        );
    }
    close(
        response.mean_support_stored_energy_j_m2,
        original.mean_stored_energy_j_m2,
        3e-13,
    );
    let rear = response.coating.field_at(thickness).unwrap();
    for axis in 0..3 {
        cc(
            rear.velocity_m_s[axis],
            response.coating.support_velocity_m_s[axis],
            3e-16,
        );
    }
}

#[test]
fn g1_quarter_wave_coating_matches_impedance_and_uses_rear_motion() {
    let state = solid(1.0);
    let film = material(2000.0, 3000.0, 1500.0, 0.001);
    let frequency = 100_000.0;
    let omega = core::f64::consts::TAU * frequency;
    let pin = WallPin {
        surface_density: 1.0,
        stiffness_per_area: omega * omega,
        resistance: 3e6,
    };
    let support = OrientedOscillatorSupport::new([pin; 3], So3::identity());
    let owner = OscillatorBackedElasticCoating::prepare(
        &state,
        &film,
        3000.0 / (4.0 * frequency),
        support,
        frequency,
        ElasticMode::P,
        1.0,
    )
    .unwrap();
    let response = owner
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    // A quarter-wave elastic slab transforms Z_back to Z_film²/Z_back.
    let effective = 6e6 * 6e6 / 3e6;
    let expected_r = (effective - 1e7) / (effective + 1e7);
    cc(
        response.coating.reflection.reflected_velocity_m_s[0],
        C64::from_re(expected_r * 0.001),
        3e-16,
    );
    let rear_velocity = C64::new(0.0, 2.0 * 0.001 * (1.0 - expected_r));
    cc(
        response.coating.support_velocity_m_s[2],
        rear_velocity,
        3e-16,
    );
    close(
        response.coating.reflection.absorbed_normal_power_w_m2,
        5.0 * (1.0 - expected_r * expected_r),
        3e-11,
    );
    let kinetic = 0.25 * rear_velocity.norm_sq();
    close(response.mean_kinetic_energy_j_m2[2], kinetic, 3e-18);
    close(response.mean_elastic_energy_j_m2[2], kinetic, 3e-18);
    // Front motion is half the magnitude: using it for pin storage gives a
    // factor-of-four error that a mere reflection-coefficient check would miss.
    let front = response.coating.reflection.boundary_velocity_m_s[2];
    close(rear_velocity.norm_sq(), 4.0 * front.norm_sq(), 3e-18);
}

#[test]
fn g1_coated_tensor_support_matches_independent_continuum_solution() {
    let state = solid(1.0);
    let film = material(1800.0, 3000.0, 1500.0, 0.001);
    let owner = OscillatorBackedElasticCoating::prepare(
        &state,
        &film,
        0.4,
        support(),
        1700.0 / core::f64::consts::TAU,
        ElasticMode::P,
        1.0,
    )
    .unwrap();
    let response = owner
        .drive([
            C64::from_re(0.001),
            C64::new(0.0004, -0.0001),
            C64::new(-0.0002, 0.0003),
        ])
        .unwrap();
    // Independent 70-digit continuum transfer exponential, with all three
    // source impedances and the original nonaligned principal wall mechanics.
    let reflected = [
        C64::new(0.0008635191696483835, 0.0004855286540046558),
        C64::new(-0.0002601286855269892, -0.0003078010899039823),
        C64::new(0.0001493027010564802, 0.0002838794726613362),
    ];
    for (actual, expected) in response
        .coating
        .reflection
        .reflected_velocity_m_s
        .into_iter()
        .zip(reflected)
    {
        cc(actual, expected, 3e-15);
    }
    close(
        response.coating.reflection.total_reflected_power_w_m2,
        5.570235389852326,
        3e-10,
    );
    close(
        response.coating.reflection.absorbed_normal_power_w_m2,
        0.17976461014767364,
        3e-11,
    );
    close(
        response.mean_support_stored_energy_j_m2,
        0.0024276631014983144,
        3e-13,
    );
    for depth in [0.0, 0.07, 0.21, 0.4] {
        let field = response.coating.field_at(depth).unwrap();
        close(field.normal_flux_w_m2, 0.17976461014767364, 3e-10);
    }
}

#[test]
fn g0_g4_coating_source_domain_and_cancellation_are_retained() {
    let state = solid(1.0);
    let film = material(1800.0, 3000.0, 1500.0, 1e-12);
    let support = support();
    let owner = OscillatorBackedElasticCoating::prepare(
        &state,
        &film,
        0.4,
        support,
        350.0,
        ElasticMode::P,
        0.7,
    )
    .unwrap();
    assert!(
        owner
            .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
            .is_err()
    );
    assert!(owner.drive([C64::ZERO; 3]).is_ok());
    assert!(
        owner
            .drive([C64::from_re(f64::NAN), C64::ZERO, C64::ZERO])
            .is_err()
    );
    for frequency in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            OscillatorBackedElasticCoating::prepare(
                &state,
                &film,
                0.4,
                support,
                frequency,
                ElasticMode::P,
                1.0
            )
            .is_err()
        );
    }
    assert!(matches!(
        OscillatorBackedElasticCoating::prepare_controlled(
            &state,
            &film,
            0.4,
            support,
            350.0,
            ElasticMode::P,
            1.0,
            &mut || true
        ),
        Err(ElasticSupportError::Cancelled)
    ));
}
