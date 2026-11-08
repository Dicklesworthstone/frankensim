use super::*;
use crate::plane_wave::uniaxial::oriented::{
    Incidence,
    absorbing::{kernel::Input, slab::kernel::solve},
    slab::SlabRegion,
};

const C0: f64 = 299_792_458.0;
const MU: f64 = 1.256_637_062_12e-6;

fn input(indices: [C64; 2]) -> Input {
    Input {
        omega: C0,
        ni: 1.2,
        mui: MU,
        mut_: MU,
        axis: [0.4, 0.7, 0.6],
        indices,
        incidence: Incidence::TangentialIndex(0.85),
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
        loss_f_m: indices.map(|n| -2.0 * n.re * n.im / (MU * C0 * C0)),
    }
}

fn relative(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs(),
        "{a:.17e} != {b:.17e}; relative tolerance {tolerance}"
    );
}

fn quadrature(r: &AbsorbingSlabResponse, area: f64, start: f64, end: f64) -> f64 {
    let n = 800;
    let step = (end - start) / n as f64;
    let heat = |z| {
        r.field_at(SlabRegion::Crystal, [0.0, 0.0, z])
            .unwrap()
            .absorbed_power_density_w_m3
    };
    let mut integral = heat(start) + heat(end);
    for j in 1..n {
        integral += (if j % 2 == 0 { 2.0 } else { 4.0 }) * heat(start + j as f64 * step);
    }
    area * integral * step / 3.0
}

#[test]
fn four_wave_volume_heat_matches_local_quadrature_flux_and_subdivision() {
    let r = solve(
        input([C64::new(1.5, -0.08), C64::new(2.0, -0.23)]),
        1.0,
        MU,
        2.3,
    )
    .unwrap();
    let area = 3.7;
    for (start, end) in [(0.0, 2.3), (0.17, 0.89), (0.9, 2.3), (1.1, 1.10000001)] {
        let actual = integrated_heat(&r, area, start, end).unwrap();
        relative(actual, quadrature(&r, area, start, end), 2e-10);
        assert!(actual > 0.0);
        if end - start > 0.01 {
            let flux = |z| {
                r.field_at(SlabRegion::Crystal, [0.0, 0.0, z])
                    .unwrap()
                    .field
                    .poynting_w_m2[2]
            };
            relative(actual, area * (flux(start) - flux(end)), 2e-13);
        }
    }
    let whole = integrated_heat(&r, area, 0.0, 2.3).unwrap();
    let pieces: f64 = [(0.0, 0.35), (0.35, 1.1), (1.1, 2.3)]
        .into_iter()
        .map(|(a, b)| integrated_heat(&r, area, a, b).unwrap())
        .sum();
    relative(pieces, whole, 3e-14);
    relative(
        whole,
        area * r.incident_normal_power_w_m2 * r.absorbed_power_fraction,
        2e-13,
    );
}

#[test]
fn original_weak_loss_remains_positive_below_port_subtraction_resolution() {
    let r = solve(
        input([C64::new(1.5, -1e-22), C64::new(1.9, -3e-22)]),
        1.0,
        MU,
        0.8,
    )
    .unwrap();
    let actual = integrated_heat(&r, 1.0, 0.0, 0.8).unwrap();
    assert!(actual > 0.0);
    assert!(actual < 1e-18 * r.incident_normal_power_w_m2);
    relative(actual, quadrature(&r, 1.0, 0.0, 0.8), 2e-11);
}

#[test]
fn lossless_dominant_polarization_does_not_erase_finite_orthogonal_heat() {
    let mut d = input([C64::new(1.5, -0.1), C64::from_re(1.9)]);
    d.axis = [1.0, 0.0, 0.0];
    d.incidence = Incidence::Cosine(1.0);
    d.electric = [C64::ZERO, C64::ONE];
    let reference = solve(d, 1.0, MU, 0.8).unwrap();
    let expected = integrated_heat(&reference, 1.0, 0.0, 0.8).unwrap();
    // Ex is lossless and 300 orders larger than absorbing Ey. Linear field
    // scaling and compensating area leave the deposited watts unchanged.
    d.electric = [C64::from_re(1e150), C64::from_re(1e-150)];
    let mixed = solve(d, 1.0, MU, 0.8).unwrap();
    let actual = integrated_heat(&mixed, 1e300, 0.0, 0.8).unwrap();
    assert!(actual > 0.0);
    relative(actual, expected, 3e-13);
}

#[test]
fn counterpropagating_node_retains_cubic_thin_volume_heat() {
    // Prescribed field fixture for the integral itself: the two waves cancel
    // exactly at z=1/2 and |E|²=4 sin²(z-1/2). The dielectric coefficient is
    // explicit; Maxwell/source consistency is exercised by the tests above.
    let mut r = solve(
        input([C64::new(1.5, -0.08), C64::new(2.0, -0.23)]),
        1.0,
        MU,
        1.0,
    )
    .unwrap();
    r.axis = [0.0, 0.0, 1.0];
    r.loss_f_m = [1e-11, 0.0];
    for w in r.forward.iter_mut().chain(&mut r.backward) {
        w.electric_v_m = [C64::ZERO; 3];
    }
    r.forward[0].normal_wavenumber_rad_m = C64::ONE;
    r.backward[0].normal_wavenumber_rad_m = -C64::ONE;
    r.forward[0].electric_v_m[0] = C64::ONE;
    r.backward[0].electric_v_m[0] = -C64::ONE;
    for requested_width in [1e-4, 1e-8, 1e-16] {
        let end = 0.5 + requested_width;
        let width = end - 0.5;
        let expected = (0.5 * r.omega * r.loss_f_m[0])
            * (4.0 / 3.0)
            * width
            * width
            * width
            * (1.0 - width * width / 5.0);
        let actual = integrated_heat(&r, 1.0, 0.5, end).unwrap();
        relative(actual, expected, 2e-13);
    }
}

#[test]
fn opaque_and_exactly_lossless_limits_need_no_growing_factors() {
    let r = solve(
        input([C64::new(1.5, -0.1), C64::new(2.0, -0.3)]),
        1.0,
        MU,
        1e5,
    )
    .unwrap();
    let power = integrated_heat(&r, 2.0, 0.0, r.thickness_m).unwrap();
    relative(
        power,
        2.0 * r.incident_normal_power_w_m2 * r.absorbed_power_fraction,
        3e-13,
    );
    assert_eq!(integrated_heat(&r, 2.0, 9e4, 1e5).unwrap(), 0.0);
    let transparent = solve(input([C64::from_re(1.5), C64::from_re(1.9)]), 1.0, MU, 1.0).unwrap();
    assert_eq!(integrated_heat(&transparent, 2.0, 0.0, 1.0).unwrap(), 0.0);
}

#[test]
fn volume_domain_is_checked_even_for_zero_loss() {
    let r = solve(input([C64::from_re(1.5), C64::from_re(1.9)]), 1.0, MU, 1.0).unwrap();
    for (area, start, end) in [
        (0.0, 0.0, 1.0),
        (1.0, -0.1, 0.5),
        (1.0, 0.5, 0.5),
        (1.0, 0.0, 1.1),
        (f64::INFINITY, 0.0, 1.0),
        (1.0, 0.0, f64::NAN),
    ] {
        assert!(integrated_heat(&r, area, start, end).is_err());
    }
    assert_eq!(admit_volume(&r, 3.0, 0.2, 0.7).unwrap(), 3.0 * (0.7 - 0.2));
}
