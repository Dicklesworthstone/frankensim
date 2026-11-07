use super::super::super::ResolvedFrequency;
use super::*;

const MU: f64 = 1.256_637_062_12e-6;
const C0: f64 = 299_792_458.0;
fn data() -> kernel::Input {
    kernel::Input {
        crystal: ResolvedFrequency {
            omega: C0,
            ni: 1.0,
            no: 1.5,
            ne: 1.9,
            mui: MU,
            mut_: MU,
        },
        exit_index: 1.2,
        exit_mu_h_m: MU,
        thickness_m: 0.7,
        axis: [1.0, 2.0, 3.0],
    }
}
fn input(p: f64, axis: [f64; 2], e: [C64; 2]) -> UniaxialIncidentMode {
    UniaxialIncidentMode {
        tangential_phase_index: p,
        incidence_plane_axis: axis,
        electric_v_m: e,
    }
}
fn close(a: C64, b: C64, t: f64) {
    assert!((a - b).abs() < t, "{a:?} != {b:?}");
}

#[test]
fn global_azimuth_rotation_preserves_both_slab_faces_and_internal_fields() {
    let d = data();
    let e = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let modes = prepare_modes(d, &[input(-0.6, [1.0, 0.0], e)], &mut || false).unwrap();
    let rotated = prepare_modes(
        kernel::Input {
            axis: [-2.0, 1.0, 3.0],
            ..d
        },
        &[input(-0.6, [0.0, 1.0], [-e[1], e[0]])],
        &mut || false,
    )
    .unwrap();
    for (region, z) in [
        (SlabRegion::TotalIncidentSide, -0.4),
        (SlabRegion::Crystal, 0.3),
        (SlabRegion::Transmitted, 1.1),
    ] {
        let original =
            field_at(&modes, d.thickness_m, region, [0.2, -0.5, z], &mut || false).unwrap();
        let turned = field_at(&rotated, d.thickness_m, region, [0.5, 0.2, z], &mut || {
            false
        })
        .unwrap();
        for (a, b) in [
            (original.electric_v_m, turned.electric_v_m),
            (original.magnetic_a_m, turned.magnetic_a_m),
        ] {
            close(b[0], -a[1], 4e-15);
            close(b[1], a[0], 4e-15);
            close(b[2], a[2], 4e-15);
        }
    }
    let several = prepare_modes(
        d,
        &[
            input(-0.6, [3.0, 4.0], e),
            input(0.4, [-4.0, 3.0], [C64::new(0.1, 0.4), C64::ONE]),
        ],
        &mut || false,
    )
    .unwrap();
    for (z, left, right) in [
        (0.0, SlabRegion::TotalIncidentSide, SlabRegion::Crystal),
        (d.thickness_m, SlabRegion::Crystal, SlabRegion::Transmitted),
    ] {
        let a = field_at(&several, d.thickness_m, left, [0.3, -0.2, z], &mut || false).unwrap();
        let b = field_at(&several, d.thickness_m, right, [0.3, -0.2, z], &mut || {
            false
        })
        .unwrap();
        for j in 0..2 {
            close(a.electric_v_m[j], b.electric_v_m[j], 5e-15);
            close(a.magnetic_a_m[j], b.magnetic_a_m[j], 5e-17);
        }
        assert!((a.poynting_w_m2[2] - b.poynting_w_m2[2]).abs() < 5e-17);
    }
}

#[test]
fn prepared_multidirection_slab_keeps_homogeneous_bright_and_dark_fringes() {
    let mut d = data();
    d.crystal.no = 1.0;
    d.crystal.ne = 1.0;
    d.exit_index = 1.0;
    let modes = prepare_modes(
        d,
        &[
            input(0.6, [1.0, 0.0], [C64::ZERO, C64::ONE]),
            input(-0.6, [1.0, 0.0], [C64::ZERO, C64::ONE]),
        ],
        &mut || false,
    )
    .unwrap();
    for (region, z) in [
        (SlabRegion::TotalIncidentSide, -0.4),
        (SlabRegion::Crystal, 0.3),
        (SlabRegion::Transmitted, 1.1),
    ] {
        let bright = field_at(&modes, d.thickness_m, region, [0.0, 0.0, z], &mut || false).unwrap();
        close(
            bright.electric_v_m[1],
            C64::new((0.8 * z).cos(), -(0.8 * z).sin()).scale(2.0),
            4e-15,
        );
        assert!((bright.poynting_w_m2[2] - 2.0 * 0.8 / (MU * C0)).abs() < 6e-17);
        let dark = field_at(
            &modes,
            d.thickness_m,
            region,
            [core::f64::consts::PI / 1.2, 0.0, z],
            &mut || false,
        )
        .unwrap();
        assert!(dark.electric_v_m[1].abs() < 1e-15);
        assert!(dark.poynting_w_m2[2].abs() < 1e-30);
        assert!((dark.magnetic_a_m[2].abs() - 1.2 / (MU * C0)).abs() < 6e-17);
    }
}

#[test]
fn slab_preparation_limits_and_cancelled_field_leave_reusable_source() {
    let limits = UniaxialSpectrumLimits {
        maximum_modes: 2,
        maximum_workspace_bytes: usize::MAX,
    };
    assert!(admit(3, limits).is_err());
    assert!(
        admit(
            2,
            UniaxialSpectrumLimits {
                maximum_workspace_bytes: 1,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        admit(
            usize::MAX,
            UniaxialSpectrumLimits {
                maximum_modes: usize::MAX,
                ..limits
            }
        )
        .is_err()
    );
    let d = data();
    let modes = prepare_modes(
        d,
        &[
            input(0.6, [1.0, 0.0], [C64::ONE, C64::ZERO]),
            input(10.0, [1.0, 0.0], [C64::ZERO; 2]),
        ],
        &mut || false,
    )
    .unwrap();
    assert_eq!(modes.len(), 1);
    let mut calls = 0;
    assert!(matches!(
        field_at(
            &modes,
            d.thickness_m,
            SlabRegion::Crystal,
            [0.0, 0.0, 0.3],
            &mut || {
                calls += 1;
                calls == 2
            }
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    let value = field_at(
        &modes,
        d.thickness_m,
        SlabRegion::Crystal,
        [0.0, 0.0, 0.3],
        &mut || false,
    )
    .unwrap();
    assert!(value.poynting_w_m2[2] > 0.0);
    assert!(
        field_at(
            &[],
            d.thickness_m,
            SlabRegion::Crystal,
            [0.0, 0.0, 0.9],
            &mut || false
        )
        .is_err()
    );
}
