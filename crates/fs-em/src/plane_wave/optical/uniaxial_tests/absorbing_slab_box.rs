//! G1/G3 physical-box controls using the original finite-film sources.
//! Native execution remains a separate check; these are authored fixtures.
use super::*;
use crate::plane_wave::uniaxial::oriented::{
    OrientedUniaxialInterface,
    absorbing::{
        AbsorbingCrystalBox, AbsorbingCrystalBoxLimits, CrystalIlluminationError,
        slab::{AbsorbingSlabCoherentIllumination, AbsorbingSlabIllumination},
    },
    slab::{OrientedUniaxialSlab, SlabRegion},
};

fn box_limits() -> AbsorbingCrystalBoxLimits {
    AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: 256,
        maximum_workspace_bytes: 1 << 20,
    }
}

// Independent tensor Simpson integration of original point E/H and tensor work.
// No analytic modal-box helper or difference of boundary powers supplies heat.
fn point_quadrature(
    source: &AbsorbingSlabCoherentIllumination<'_, '_, '_, '_>,
    region: AbsorbingCrystalBox,
    n: usize,
) -> (f64, [f64; 6]) {
    assert!(n > 0 && n % 2 == 0);
    let origin = region.origin_m();
    let extent = region.extent_m();
    let w = |i| if i == 0 || i == n { 1.0 } else if i % 2 == 0 { 2.0 } else { 4.0 };
    let mut heat = 0.0;
    for i in 0..=n {
        for j in 0..=n {
            for k in 0..=n {
                let index = [i, j, k];
                let point = core::array::from_fn(|axis| {
                    origin[axis] + extent[axis] * (index[axis] as f64 / n as f64)
                });
                heat += w(i) * w(j) * w(k)
                    * source.field_at(SlabRegion::Crystal, point, || false)
                        .unwrap().absorbed_power_density_w_m3;
            }
        }
    }
    heat *= region.volume_m3() / (27.0 * (n as f64).powi(3));
    let mut faces = [0.0; 6];
    for axis in 0..3 {
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        for upper in 0..2 {
            let mut sum = 0.0;
            for i in 0..=n {
                for j in 0..=n {
                    let mut point = origin;
                    point[axis] += upper as f64 * extent[axis];
                    point[u] += extent[u] * (i as f64 / n as f64);
                    point[v] += extent[v] * (j as f64 / n as f64);
                    sum += w(i) * w(j)
                        * source.poynting_w_m2(SlabRegion::Crystal, point, || false).unwrap()[axis];
                }
            }
            faces[2 * axis + upper] = sum * extent[u] * extent[v] / (9.0 * (n as f64).powi(2))
                * if upper == 0 { -1.0 } else { 1.0 };
        }
    }
    (heat, faces)
}

#[test]
fn opaque_two_sided_film_preserves_actual_boundary_planes_and_lateral_work() {
    let outside = incident(1.0, MU);
    let exit = incident(1.3, MU);
    let (o, e) = axes([1.5; 2], [1.9; 2], [0.35, 0.5], Some(MU), [400.0, 800.0], 2);
    let law = UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap();
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside), law, [0.0, 0.0, 1.0],
    ).unwrap();
    let wavelength = 550e-9;
    let thickness = 400.0 * wavelength;
    let film = OrientedUniaxialSlab::try_new(
        &entrance, DielectricResponse::Lossless(&exit), thickness,
    ).unwrap();
    let fields = [film.prepare_absorbing_field_from_both_sides_at_tangential_index(
        frequency(550.0), 0.35,
        [C64::ZERO, C64::ONE], [C64::ZERO, C64::new(0.7, 0.4)],
    ).unwrap()];
    let source = AbsorbingSlabCoherentIllumination::try_new(&film, &fields, 1, || false).unwrap();
    let extent = [0.3 * wavelength, 0.4 * wavelength, thickness];
    let region = AbsorbingCrystalBox::try_new([0.13 * wavelength, -0.2 * wavelength, 0.0], extent).unwrap();
    let result = source.box_absorption(region, box_limits(), || false).unwrap();
    let area = extent[0] * extent[1];
    let heat = fields[0].slab_absorption(area, 0.0, thickness).unwrap().absorbed_power_w();
    near(result.absorbed_power_w(), heat, 4e-11 * heat);
    let response = fields[0].response();
    assert!(-response.forward[0].normal_wavenumber_rad_m.im * thickness > 750.0);
    // TE Maxwell identity: Sx = kx*q / (mu*omega²*epsilon_loss).
    // It holds for the complete forward/backward standing wave at every z.
    let ratio = response.tangential_wavenumber_rad_m
        / (MU * response.omega_rad_s() * response.omega_rad_s()
            * response.principal_loss_permittivities_f_m()[0]);
    let lateral = heat * ratio / extent[0];
    let lower = fields[0].field_at(SlabRegion::Crystal, [0.0; 3]).unwrap().field.poynting_w_m2[2];
    let upper = fields[0].field_at(SlabRegion::Crystal, [0.0, 0.0, thickness]).unwrap().field.poynting_w_m2[2];
    let expected = [-lateral, lateral, 0.0, 0.0, -area * lower, area * upper];
    let scale = expected.iter().map(|v| v.abs()).sum::<f64>();
    for (actual, expected) in result.outward_face_power_w().into_iter().zip(expected) {
        near(actual, expected, 5e-11 * scale);
    }
    assert!(expected[4] < 0.0 && expected[5] < 0.0);
    near(result.residual_w(), 0.0, 5e-11 * scale);
    assert!(core::ptr::eq(result.source(), &film));
    assert!(core::ptr::eq(result.illumination(), &source));
    assert_eq!(result.mode_pair_count(), 16);
}

#[test]
fn tilted_crossed_two_port_beams_match_original_local_work_and_six_face_partition() {
    let outside = incident(1.0, MU);
    let exit = incident(1.25, MU);
    let (o, e) = axes([1.5; 2], [1.9; 2], [0.08, 0.2], Some(MU), [400.0, 800.0], 3);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap(), [0.3, -0.4, 0.8],
    ).unwrap();
    let wavelength = 550e-9;
    let film = OrientedUniaxialSlab::try_new(
        &entrance, DielectricResponse::Lossless(&exit), wavelength,
    ).unwrap();
    let fields = [
        film.prepare_absorbing_field_from_both_sides_at_tangential_index_in_plane(
            frequency(550.0), 0.25, [1.0, 0.0],
            [C64::new(0.7, 0.2), C64::new(0.3, -0.1)],
            [C64::new(0.2, -0.3), C64::new(0.1, 0.4)],
        ).unwrap(),
        film.prepare_absorbing_field_from_both_sides_at_tangential_index_in_plane(
            frequency(550.0), 0.7, [0.6, 0.8],
            [C64::new(0.1, 0.6), C64::new(0.5, 0.2)],
            [C64::new(-0.2, 0.1), C64::new(0.3, -0.2)],
        ).unwrap(),
    ];
    let source = AbsorbingSlabCoherentIllumination::try_new(&film, &fields, 2, || false).unwrap();
    let region = AbsorbingCrystalBox::try_new(
        [0.1 * wavelength, -0.2 * wavelength, 0.2 * wavelength],
        [0.29 * wavelength, 0.21 * wavelength, 0.15 * wavelength],
    ).unwrap();
    let result = source.box_absorption(region, box_limits(), || false).unwrap();
    let (heat, expected_faces) = point_quadrature(&source, region, 24);
    near(result.absorbed_power_w(), heat, 3e-6 * heat);
    let scale = expected_faces.iter().map(|v| v.abs()).sum::<f64>();
    for (actual, expected) in result.outward_face_power_w().into_iter().zip(expected_faces) {
        near(actual, expected, 3e-6 * scale);
    }
    near(result.residual_w(), 0.0, 8e-11 * scale);
    assert!(expected_faces[..4].iter().map(|v| v.abs()).sum::<f64>() > 1e-4 * scale);
    let origin = region.origin_m();
    let extent = region.extent_m();
    let mut sum_heat = 0.0;
    let mut sum_faces = [0.0; 6];
    for index in 0..8 {
        let side = [index & 1, (index >> 1) & 1, (index >> 2) & 1];
        let half = extent.map(|value| value * 0.5);
        let start = core::array::from_fn(|axis| origin[axis] + half[axis] * side[axis] as f64);
        let part = source.box_absorption(
            AbsorbingCrystalBox::try_new(start, half).unwrap(), box_limits(), || false,
        ).unwrap();
        sum_heat += part.absorbed_power_w();
        for axis in 0..3 {
            sum_faces[2 * axis + side[axis]] += part.outward_face_power_w()[2 * axis + side[axis]];
        }
    }
    near(sum_heat, result.absorbed_power_w(), 5e-11 * result.absorbed_power_w());
    for (actual, expected) in sum_faces.into_iter().zip(result.outward_face_power_w()) {
        near(actual, expected, 5e-11 * scale);
    }
}

#[test]
fn short_standing_wave_node_retains_positive_cubic_heat() {
    let outside = incident(1.0, MU);
    let (o, e) = axes([1.4; 2], [2.1; 2], [0.12, 0.42], Some(MU), [400.0, 800.0], 1);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap(), [0.0, 0.0, 1.0],
    ).unwrap();
    let wavelength = 550e-9;
    let thickness = 0.8 * wavelength;
    let film = OrientedUniaxialSlab::try_new(
        &entrance, DielectricResponse::Lossless(&outside), thickness,
    ).unwrap();
    let fields = [film.prepare_absorbing_field_from_both_sides_at_tangential_index(
        frequency(550.0), 0.3, [C64::ZERO, C64::ONE], [C64::ZERO, C64::new(-1.0, 0.0)],
    ).unwrap()];
    let source = AbsorbingSlabCoherentIllumination::try_new(&film, &fields, 1, || false).unwrap();
    let area = (0.27 * wavelength) * (0.31 * wavelength);
    // Symmetric ports with opposite drives give an electric node at d/2.
    // Positive short-box field integration retains the |z-d/2|² loss.
    for fraction in [1e-4, 1e-7] {
        let mut powers = [0.0; 2];
        for (index, shrink) in [1.0, 0.5].into_iter().enumerate() {
            let width = wavelength * fraction * shrink;
            let start = thickness * 0.5 - width * 0.5;
            let region = AbsorbingCrystalBox::try_new(
                [-0.1 * wavelength, 0.2 * wavelength, start],
                [0.27 * wavelength, 0.31 * wavelength, width],
            ).unwrap();
            let result = source.box_absorption(region, box_limits(), || false).unwrap();
            let reference = fields[0].slab_absorption(area, start, start + width).unwrap();
            assert!(result.absorbed_power_w() > 0.0);
            near(result.absorbed_power_w(), reference.absorbed_power_w(),
                2e-7 * reference.absorbed_power_w());
            powers[index] = result.absorbed_power_w();
        }
        near(powers[0] / powers[1], 8.0, 2e-5);
    }
}

#[test]
fn spectral_group_boxes_bound_total_work_peak_storage_and_exact_film_support() {
    let outside = incident(1.0, MU);
    let exit = incident(1.2, MU);
    let (o, e) = axes([1.3, 1.8], [1.7, 2.3], [0.1, 0.3], Some(MU), [400.0, 800.0], 4);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&o, &e, SPAN).unwrap(), [0.0, 0.0, 1.0],
    ).unwrap();
    let thickness = 0.9e-6;
    let film = OrientedUniaxialSlab::try_new(&entrance, DielectricResponse::Lossless(&exit), thickness).unwrap();
    let angular = [0.2, -0.4].map(|p| film.prepare_absorbing_field_from_both_sides_at_tangential_index(
        frequency(500.0), p, [C64::ZERO, C64::ONE], [C64::ZERO, C64::new(0.3, 0.4)],
    ).unwrap());
    let red = [film.prepare_absorbing_field_from_both_sides_at_tangential_index(
        frequency(700.0), 0.6, [C64::new(0.4, -0.2), C64::ONE], [C64::ZERO; 2],
    ).unwrap()];
    let groups = [
        AbsorbingSlabCoherentIllumination::try_new(&film, &angular, 2, || false).unwrap(),
        AbsorbingSlabCoherentIllumination::try_new(&film, &red, 1, || false).unwrap(),
    ];
    let bank = AbsorbingSlabIllumination::try_new(&film, &groups, 2, 3, || false).unwrap();
    let one = AbsorbingSlabIllumination::try_new(&film, &groups[..1], 1, 2, || false).unwrap();
    let region = AbsorbingCrystalBox::try_new([0.13e-6, -0.17e-6, 0.1e-6], [0.3e-6, 0.4e-6, 0.6e-6]).unwrap();
    let result = bank.box_absorption(region, box_limits(), || false).unwrap();
    let first = groups[0].box_absorption(region, box_limits(), || false).unwrap();
    let second = groups[1].box_absorption(region, box_limits(), || false).unwrap();
    near(result.absorbed_power_w(), first.absorbed_power_w() + second.absorbed_power_w(),
        3e-13 * result.absorbed_power_w());
    let scale = result.outward_face_power_w().iter().map(|v| v.abs()).sum::<f64>();
    for face in 0..6 {
        near(result.outward_face_power_w()[face],
            first.outward_face_power_w()[face] + second.outward_face_power_w()[face], 3e-13 * scale);
    }
    assert!(core::ptr::eq(result.illumination(), &bank));
    assert!(core::ptr::eq(result.source(), &film));
    assert_eq!(result.region(), region);
    assert_eq!(result.volume_m3(), region.volume_m3());
    assert_eq!(groups[0].mode_pair_count(), 64);
    assert_eq!(groups[1].mode_pair_count(), 16);
    assert_eq!(bank.mode_pair_count(), 80);
    assert_eq!(result.mode_pair_count(), 80);
    assert_eq!(result.workspace_bytes(), one.box_absorption(region, box_limits(), || false).unwrap().workspace_bytes());
    let exact = AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: result.mode_pair_count(),
        maximum_workspace_bytes: result.workspace_bytes(),
    };
    let repeated = bank.box_absorption(region, exact, || false).unwrap();
    assert_eq!(repeated.absorbed_power_w(), result.absorbed_power_w());
    assert_eq!(repeated.outward_face_power_w(), result.outward_face_power_w());
    assert!(bank.box_absorption(region, AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: 79, ..exact
    }, || false).is_err());
    assert!(bank.box_absorption(region, AbsorbingCrystalBoxLimits {
        maximum_workspace_bytes: exact.maximum_workspace_bytes - 1, ..exact
    }, || false).is_err());
    let coherent_exact = AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: first.mode_pair_count(), maximum_workspace_bytes: first.workspace_bytes(),
    };
    assert!(groups[0].box_absorption(region, coherent_exact, || false).is_ok());
    assert!(groups[0].box_absorption(region, AbsorbingCrystalBoxLimits {
        maximum_workspace_bytes: coherent_exact.maximum_workspace_bytes - 1, ..coherent_exact
    }, || false).is_err());
    let mut count = 0;
    bank.box_absorption(region, exact, || { count += 1; false }).unwrap();
    for cutoff in [1, count / 2, count] {
        let mut calls = 0;
        assert!(matches!(bank.box_absorption(region, exact, || {
            calls += 1;
            calls >= cutoff
        }), Err(CrystalIlluminationError::Cancelled)));
    }
    let empty_group = AbsorbingSlabCoherentIllumination::try_new(&film, &[], 0, || false).unwrap();
    let empty_bank = AbsorbingSlabIllumination::try_new(&film, &[], 0, 0, || false).unwrap();
    let tiny = thickness * f64::EPSILON * 0.125;
    assert_eq!(thickness + tiny, thickness);
    let outside_support = AbsorbingCrystalBox::try_new([0.0, 0.0, thickness], [1.0, 1.0, tiny]).unwrap();
    assert!(empty_group.box_absorption(outside_support, box_limits(), || false).is_err());
    assert!(empty_bank.box_absorption(outside_support, box_limits(), || false).is_err());
    let dark = empty_bank.box_absorption(region, AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: 0, ..box_limits()
    }, || false).unwrap();
    assert_eq!(dark.absorbed_power_w(), 0.0);
    assert_eq!(dark.outward_face_power_w(), [0.0; 6]);
}
