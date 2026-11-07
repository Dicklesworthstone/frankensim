//! Cross-direction interference is retained until the final physical sum.
use super::super::super::super::kernel::UniaxialWaveKind;
use super::super::super::slab::SlabWave;
use super::*;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

fn fixture() -> Frequency {
    let mu = 1.256_637_062_12e-6;
    let mut entrance = 0.0;
    let layers = [0.6, 0.4]
        .into_iter()
        .map(|thickness| {
            let exit = entrance + thickness;
            let layer = super::super::ResolvedLayer {
                no: 1.0,
                ne: 1.0,
                mu,
                axis: [1.0, 0.0, 0.0],
                thickness,
                entrance,
                exit,
            };
            entrance = exit;
            layer
        })
        .collect();
    Frequency {
        omega: super::super::C0,
        ni: 1.0,
        nt: 1.0,
        mui: mu,
        mut_: mu,
        layers,
    }
}

#[test]
fn solved_nearby_directions_keep_relative_phase_at_large_coordinates() {
    let data = fixture();
    let p = [0.3, f64::from_bits(0.3_f64.to_bits() + 1)];
    let input = core::array::from_fn::<_, 2, _>(|j| UniaxialIncidentMode {
        tangential_phase_index: p[j],
        incidence_plane_axis: [1.0, 0.0],
        electric_v_m: [C64::ZERO, if j == 0 { C64::ONE } else { -C64::ONE }],
    });
    let modes = prepare_modes(&data, &input, &mut || false).unwrap();
    let delta = modes[1].response.tangential_wavenumber_rad_m
        - modes[0].response.tangential_wavenumber_rad_m;
    assert!(delta > 0.0);
    let admittance = -modes[0].response.incident_magnetic_a_m[0].re
        + modes[1].response.incident_magnetic_a_m[0].re;
    for x in [-1e16, 1e16] {
        let at = field_at(&modes, StackRegion::Incident, [x, 0.0, 0.0], &mut || false).unwrap();
        let sine = (0.5 * delta * x).sin();
        close(at.electric_v_m[1].abs(), 2.0 * sine.abs(), 3e-15);
        close(at.poynting_w_m2[2] / (admittance * sine * sine), 1.0, 3e-14);
    }
}

#[test]
fn all_directions_share_one_layer_sum_before_rounding_weak_fields() {
    let data = fixture();
    let mut response = kernel::solve(
        &data,
        Incidence::TangentialIndex(0.0),
        [C64::ONE, C64::ZERO],
        [1.0, 0.0],
        &mut || false,
    )
    .unwrap();
    let zero = SlabWave {
        wave_kind: UniaxialWaveKind::Propagating,
        normal_wavenumber_rad_m: C64::ONE,
        electric_v_m: [C64::ZERO; 3],
        magnetic_a_m: [C64::ZERO; 3],
    };
    response.layers[0].forward = [zero; 2];
    response.layers[0].backward = [zero; 2];
    response.tangential_wavenumber_rad_m = 1e20;
    let k = [1.0, f64::from_bits(1.0_f64.to_bits() + 1)];
    let modes = core::array::from_fn::<_, 2, _>(|j| {
        let mut response = response.clone();
        let a = if j == 0 { C64::ONE } else { -C64::ONE };
        response.layers[0].forward[0] = SlabWave {
            normal_wavenumber_rad_m: C64::from_re(k[j]),
            electric_v_m: [a, C64::ZERO, C64::ZERO],
            magnetic_a_m: [C64::ZERO, a, C64::ZERO],
            ..zero
        };
        Mode {
            axis: [1.0, 0.0],
            response,
        }
    });
    let depth = 1e-4;
    let at = field_at(
        &modes,
        StackRegion::Layer(0),
        [1.0, 0.0, depth],
        &mut || false,
    )
    .unwrap();
    let amplitude = (k[1] - k[0]) * depth;
    close(at.electric_v_m[0].abs() / amplitude, 1.0, 3e-14);
    close(
        at.poynting_w_m2[2] / (0.5 * amplitude * amplitude),
        1.0,
        3e-14,
    );
}
