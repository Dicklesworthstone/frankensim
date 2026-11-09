use super::{close, pins, solid};
use fs_couple::elastic_reflection::{
    ElasticSupportError, OrientedOscillatorSupport, OrientedOscillatorSupportedSurface,
    OscillatorSupportedSurface,
};
use fs_ga::{So3, So3Tangent, Vec3};
use fs_material::elastic::acoustic::interface::ElasticMode;
use fs_math::c64::C64;
use fs_phs::{WallPin, wall_specific_impedance};

fn rotation(x: f64, y: f64, z: f64) -> So3 {
    So3::exp(So3Tangent::new(Vec3::new(x, y, z))).unwrap()
}

fn complex_close(a: C64, b: C64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:?} != {b:?}");
}

#[test]
fn g1_aligned_supports_recover_original_two_pin_reflection_and_storage() {
    let state = solid(1.0);
    let [normal, tangent] = pins(3e6);
    let support = OrientedOscillatorSupport::new([tangent, tangent, normal], So3::identity());
    let input = [
        C64::new(0.001, 0.0002),
        C64::new(-0.0004, 0.0003),
        C64::from_re(0.0006),
    ];
    for omega in [1000.0, 1700.0] {
        let original = OscillatorSupportedSurface::prepare(
            &state,
            normal,
            tangent,
            omega,
            ElasticMode::P,
            0.7,
        )
        .unwrap()
        .drive(input)
        .unwrap();
        let owner = OrientedOscillatorSupportedSurface::prepare(
            &state,
            support,
            omega,
            ElasticMode::P,
            0.7,
        )
        .unwrap();
        assert!(core::ptr::eq(owner.state(), &state));
        assert_eq!(owner.support(), &support);
        assert_eq!(owner.omega_rad_s(), omega);
        let oriented = owner.drive(input).unwrap();
        for axis in 0..3 {
            complex_close(
                oriented.reflection.reflected_velocity_m_s[axis],
                original.reflection.reflected_velocity_m_s[axis],
                2e-16,
            );
            complex_close(
                oriented.support_displacement_m[axis],
                original.support_displacement_m[axis],
                2e-18,
            );
            complex_close(
                oriented.principal_velocity_m_s[axis],
                original.reflection.boundary_velocity_m_s[axis],
                2e-16,
            );
            close(
                oriented.mean_kinetic_energy_j_m2[axis],
                original.mean_kinetic_energy_j_m2[axis],
                2e-14,
            );
            close(
                oriented.mean_elastic_energy_j_m2[axis],
                original.mean_elastic_energy_j_m2[axis],
                2e-14,
            );
        }
        close(
            oriented.reflection.absorbed_normal_power_w_m2,
            original.reflection.absorbed_normal_power_w_m2,
            2e-12,
        );
        close(
            oriented.mean_stored_energy_j_m2,
            original.mean_stored_energy_j_m2,
            3e-14,
        );
    }
}

#[test]
fn g1_rotated_matched_and_lossless_pins_convert_shear_and_split_energy() {
    let state = solid(1.0);
    let [normal, matched] = pins(1e7);
    let free = WallPin {
        resistance: 0.0,
        ..matched
    };
    let support = OrientedOscillatorSupport::new(
        [matched, free, normal],
        rotation(0.0, 0.0, core::f64::consts::FRAC_PI_4),
    );
    let owner =
        OrientedOscillatorSupportedSurface::prepare(&state, support, 1000.0, ElasticMode::Sh, 1.0)
            .unwrap();
    let response = owner
        .drive([C64::ZERO, C64::ZERO, C64::from_re(0.001)])
        .unwrap();
    // In the surface plane the reflection map is the projector onto the
    // lossless principal axis (-1,1)/sqrt(2); the matched axis reflects zero.
    complex_close(
        response.reflection.reflected_velocity_m_s[0],
        C64::ZERO,
        1e-17,
    );
    complex_close(
        response.reflection.reflected_velocity_m_s[1],
        C64::from_re(-0.0005),
        2e-17,
    );
    complex_close(
        response.reflection.reflected_velocity_m_s[2],
        C64::from_re(0.0005),
        2e-17,
    );
    close(response.reflection.absorbed_normal_power_w_m2, 1.25, 3e-13);
    close(response.reflection.total_reflected_power_w_m2, 1.25, 3e-13);
    close(response.mean_kinetic_energy_j_m2[0], 0.0025, 3e-16);
    close(response.mean_kinetic_energy_j_m2[1], 0.01, 3e-16);
    close(response.mean_stored_energy_j_m2, 0.025, 1e-15);
    // Drive only the authored lossless principal direction. Rounded Cartesian
    // resistance must not introduce a weak extra dissipative channel here.
    let r = support.principal_to_surface().matrix().m;
    let lossless = owner
        .drive([
            C64::ZERO,
            C64::from_re(0.001 * r[1]),
            C64::from_re(0.001 * r[4]),
        ])
        .unwrap();
    assert!(
        lossless.reflection.absorbed_normal_power_w_m2
            < 1e-25 * lossless.reflection.total_incoming_power_w_m2
    );
}

#[test]
fn g1_tilted_normal_support_converts_p_to_sh_against_independent_elimination() {
    let state = solid(1.0);
    let [normal, tangent] = pins(3e6);
    let support =
        OrientedOscillatorSupport::new([tangent, tangent, normal], rotation(0.4, 0.0, 0.0));
    let omega = 1700.0;
    let r = support.principal_to_surface().matrix().m;
    let zt = wall_specific_impedance(&tangent, omega).unwrap();
    let zn = wall_specific_impedance(&normal, omega).unwrap();
    let zyy = zt.scale(r[4] * r[4]) + zn.scale(r[5] * r[5]);
    let zzz = zt.scale(r[7] * r[7]) + zn.scale(r[8] * r[8]);
    let zyz = zt.scale(r[4] * r[7]) + zn.scale(r[5] * r[8]);
    let shear = zyy + C64::from_re(5e6);
    let effective_normal = zzz - zyz * zyz / shear;
    let reflected_p =
        (effective_normal - C64::from_re(1e7)) / (effective_normal + C64::from_re(1e7));
    let reflected_sh = -(zyz / shear) * (C64::ONE - reflected_p);
    let response =
        OrientedOscillatorSupportedSurface::prepare(&state, support, omega, ElasticMode::P, 1.0)
            .unwrap()
            .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
            .unwrap();
    complex_close(
        response.reflection.reflected_velocity_m_s[0],
        reflected_p.scale(0.001),
        3e-16,
    );
    complex_close(
        response.reflection.reflected_velocity_m_s[1],
        C64::ZERO,
        1e-17,
    );
    complex_close(
        response.reflection.reflected_velocity_m_s[2],
        reflected_sh.scale(0.001),
        3e-16,
    );
    assert!(response.reflection.reflected_velocity_m_s[2].abs() > 1e-5);
}

#[test]
fn g3_coherent_oriented_pins_preserve_original_work_and_each_storage_potential() {
    let state = solid(1.0);
    let [normal, tangent] = pins(3e6);
    let principal = [
        tangent,
        WallPin {
            surface_density: 4e4,
            stiffness_per_area: 3e10,
            resistance: 2e6,
        },
        normal,
    ];
    let support = OrientedOscillatorSupport::new(principal, rotation(0.2, -0.4, 0.6));
    let omega = 1700.0;
    let owner =
        OrientedOscillatorSupportedSurface::prepare(&state, support, omega, ElasticMode::P, 0.7)
            .unwrap();
    let input = [
        C64::new(0.001, 0.0003),
        C64::new(-0.0007, 0.0002),
        C64::new(0.0004, -0.0001),
    ];
    let response = owner.drive(input).unwrap();
    let rotation = support.principal_to_surface().matrix().m;
    let mut loss = 0.0;
    let mut stored = 0.0;
    for axis in 0..3 {
        let v = response.principal_velocity_m_s[axis];
        let x = response.principal_displacement_m[axis];
        let mut principal_traction = C64::ZERO;
        let mut principal_velocity = C64::ZERO;
        for row in 0..3 {
            principal_traction = principal_traction
                + response.reflection.boundary_traction_pa[row].scale(rotation[row * 3 + axis]);
            principal_velocity = principal_velocity
                + response.reflection.boundary_velocity_m_s[row].scale(rotation[row * 3 + axis]);
        }
        complex_close(v, principal_velocity, 1e-17);
        complex_close(x * C64::new(0.0, -omega), v, 1e-17);
        complex_close(
            principal_traction + wall_specific_impedance(&principal[axis], omega).unwrap() * v,
            C64::ZERO,
            2e-8,
        );
        let kinetic = 0.25 * principal[axis].surface_density * v.norm_sq();
        let elastic = 0.25 * principal[axis].stiffness_per_area * x.norm_sq();
        close(response.mean_kinetic_energy_j_m2[axis], kinetic, 2e-14);
        close(response.mean_elastic_energy_j_m2[axis], elastic, 2e-14);
        close(
            0.5 * (principal_traction * v.conj()).im,
            2.0 * omega * (kinetic - elastic),
            3e-11,
        );
        loss += 0.5 * principal[axis].resistance * v.norm_sq();
        stored += kinetic + elastic;
        complex_close(
            response.support_displacement_m[axis] * C64::new(0.0, -omega),
            response.reflection.boundary_velocity_m_s[axis],
            1e-17,
        );
    }
    close(response.reflection.absorbed_normal_power_w_m2, loss, 2e-12);
    close(response.mean_stored_energy_j_m2, stored, 3e-14);
    close(
        response.reflection.total_incoming_power_w_m2,
        response.reflection.total_reflected_power_w_m2 + loss,
        3e-12,
    );
    let replay = owner.drive(input).unwrap();
    assert_eq!(
        replay.principal_velocity_m_s,
        response.principal_velocity_m_s
    );
    assert_eq!(
        replay.mean_stored_energy_j_m2,
        response.mean_stored_energy_j_m2
    );
}

#[test]
fn g0_oriented_support_keeps_lossless_resonance_original_domains_and_cancellation() {
    let state = solid(1.0);
    let [normal, tangent] = pins(0.0);
    let support =
        OrientedOscillatorSupport::new([tangent, tangent, normal], rotation(0.2, 0.4, 0.6));
    let owner =
        OrientedOscillatorSupportedSurface::prepare(&state, support, 1000.0, ElasticMode::P, 1.0)
            .unwrap();
    let response = owner
        .drive([C64::from_re(0.001), C64::ZERO, C64::ZERO])
        .unwrap();
    assert_eq!(response.reflection.absorbed_normal_power_w_m2, 0.0);
    complex_close(
        response.reflection.reflected_velocity_m_s[0],
        C64::from_re(-0.001),
        1e-17,
    );
    assert!(response.mean_stored_energy_j_m2 > 0.0);
    assert_eq!(
        owner.drive([C64::ZERO; 3]).unwrap().mean_stored_energy_j_m2,
        0.0
    );
    for omega in [0.0, -1.0, f64::NAN] {
        assert!(
            OrientedOscillatorSupportedSurface::prepare(
                &state,
                support,
                omega,
                ElasticMode::P,
                1.0
            )
            .is_err()
        );
    }
    let invalid = OrientedOscillatorSupport::new(
        [
            WallPin {
                resistance: -1.0,
                ..tangent
            },
            tangent,
            normal,
        ],
        So3::identity(),
    );
    assert!(
        OrientedOscillatorSupportedSurface::prepare(&state, invalid, 1000.0, ElasticMode::P, 1.0)
            .is_err()
    );
    assert!(
        owner
            .drive([C64::from_re(1e9), C64::ZERO, C64::ZERO])
            .is_err()
    );
    assert!(matches!(
        OrientedOscillatorSupportedSurface::prepare_controlled(
            &state,
            support,
            1000.0,
            ElasticMode::P,
            1.0,
            &mut || true
        ),
        Err(ElasticSupportError::Cancelled)
    ));
}
