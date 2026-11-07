use super::*;

fn layers(depths: &[f64], indices: &[f64]) -> Vec<AbsorbingLayer> {
    depths
        .iter()
        .zip(indices)
        .map(|(&depth, &index)| {
            AbsorbingLayer::new(1.0, depth, [100.0, 3000.0], "beam absorption")
                .unwrap()
                .with_refractive_index(index, "beam index")
                .unwrap()
        })
        .collect()
}

fn beam(depths: &[f64], indices: &[f64], mu: f64, s: f64) -> SlabBeamResponse {
    SlabBeamResponse::new(layers(depths, indices), mu, s, || false).unwrap()
}

#[test]
fn g1_polarized_oblique_beams_match_independent_boundary_equations() {
    let depths = [0.2, 0.9, 0.4];
    let indices = [1.4, 2.1, 1.2];
    let h = [11.0, 5.0];
    for mu in [1.0_f64, 0.7, 0.2] {
        for s in [0.0, 0.3, 1.0] {
            let response = beam(&depths, &indices, mu, s);
            let p = ((1.0 - mu) * (1.0 + mu)).sqrt();
            let (qs, js) = heterogeneous::ray(&depths, &indices, p, 0, &[0.0; 3], h);
            let (qp, jp) = heterogeneous::ray(&depths, &indices, p, 1, &[0.0; 3], h);
            let mut heat = [0.0; 3];
            let report = response.exchange(h, &mut heat).unwrap();
            for i in 0..3 {
                close(heat[i], s * qs[i] + (1.0 - s) * qp[i], 2e-13);
            }
            for i in 0..2 {
                close(
                    report.reflected_w_m2[i] + report.transmitted_w_m2[i],
                    s * js[i] + (1.0 - s) * jp[i],
                    2e-13,
                );
            }
            close(report.absorbed_w_m2, heat.iter().sum(), 2e-14);
            assert!(report.balance_residual_w_m2.abs() < 2e-13);
        }
    }
}

#[test]
fn g1_beer_brewster_and_clear_slab_limits_resolve_polarization() {
    let matched = beam(&[0.2, 0.4], &[1.0, 1.0], 0.3, 0.37);
    let transmission = (-2.0_f64).exp();
    assert_eq!(matched.reflectances(), [0.0; 2]);
    close(matched.transmittance(), transmission, 2e-14);
    let mut heat = [0.0; 2];
    matched.exchange([7.0, 0.0], &mut heat).unwrap();
    close(heat[0], 7.0 * -(-2.0_f64 / 3.0).exp_m1(), 2e-14);
    close(heat.iter().sum(), 7.0 * (1.0 - transmission), 2e-14);

    let clear = beam(&[0.0], &[1.5], 1.0, 0.5);
    close(clear.reflectances()[0], 2.0 * 0.04 / 1.04, 2e-14);
    close(clear.transmittance(), 0.96 / 1.04, 2e-14);
    assert_eq!(clear.layer_absorptance(), &[[0.0; 2]]);
    let brewster = (1.0_f64 + 1.5 * 1.5).sqrt().recip();
    let p = beam(&[0.0], &[1.5], brewster, 0.0);
    let s = beam(&[0.0], &[1.5], brewster, 1.0);
    assert!(p.reflectances()[0] < 1e-29);
    assert!(s.reflectances()[0] > 0.2);
    close(p.transmittance(), 1.0, 2e-14);
}

#[test]
fn g1_directional_integration_recovers_original_diffuse_coefficients() {
    let slab = AbsorbingSlab::new(
        layers(&[0.2, 0.9, 0.4], &[1.4, 2.1, 1.2]),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let mut r = [0.0; 2];
    let mut t = 0.0;
    let mut a = [[0.0; 2]; 3];
    let panels = 64;
    for panel in 0..panels {
        for &(node, weight) in &GL16 {
            for sign in [-1.0, 1.0] {
                let mu = (panel as f64 + 0.5 + 0.5 * sign * node) / panels as f64;
                let q = mu * weight / panels as f64;
                let ray = slab.beam_response(mu, 0.5, || false).unwrap();
                for side in 0..2 {
                    r[side] += q * ray.reflectances()[side];
                    for i in 0..3 {
                        a[i][side] += q * ray.layer_absorptance()[i][side];
                    }
                }
                t += q * ray.transmittance();
            }
        }
    }
    close(t, slab.transmittance(), 2e-11);
    for side in 0..2 {
        close(r[side], slab.reflectances()[side], 2e-11);
        for i in 0..3 {
            close(a[i][side], slab.face_coupling()[i][side], 2e-11);
        }
    }
}

#[test]
fn g1_directional_reversal_subdivision_and_native_thermal_superposition() {
    let slab = AbsorbingSlab::new(
        layers(&[0.2, 0.9], &[1.4, 2.1]),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap();
    let original = slab.beam_response(0.4, 0.7, || false).unwrap();
    let reversed = beam(&[0.9, 0.2], &[2.1, 1.4], 0.4, 0.7);
    let split = beam(&[0.1, 0.1, 0.9], &[1.4, 1.4, 2.1], 0.4, 0.7);
    close(original.transmittance(), reversed.transmittance(), 2e-14);
    close(original.transmittance(), split.transmittance(), 2e-14);
    for side in 0..2 {
        close(
            original.reflectances()[side],
            reversed.reflectances()[1 - side],
            2e-14,
        );
        close(
            original.reflectances()[side],
            split.reflectances()[side],
            2e-14,
        );
        close(
            original.layer_absorptance()[0][side],
            split.layer_absorptance()[0][side] + split.layer_absorptance()[1][side],
            2e-14,
        );
    }
    let mut directed = [0.0; 2];
    original.exchange([13.0, 7.0], &mut directed).unwrap();
    let mut diffuse = [0.0; 2];
    slab.layer_net_fluxes([11.0, 5.0], &[3.0, 7.0], &mut diffuse)
        .unwrap();
    let mut total = [0.0; 2];
    slab.layer_net_fluxes_with_beam([11.0, 5.0], &[3.0, 7.0], &original, [13.0, 7.0], &mut total)
        .unwrap();
    for i in 0..2 {
        close(total[i], diffuse[i] + directed[i], 2e-14);
    }
    slab.layer_net_fluxes_with_beam([31.0; 2], &[31.0; 2], &original, [13.0, 7.0], &mut total)
        .unwrap();
    assert_eq!(total, directed);
}

#[test]
fn g1_weak_absorption_and_grazing_beam_escape_remain_positive() {
    let weak = beam(&[1e-20], &[1.5], 0.6, 0.5);
    let doubled = beam(&[2e-20], &[1.5], 0.6, 0.5);
    assert!(weak.layer_absorptance()[0][0] > 0.0);
    close(
        doubled.layer_absorptance()[0][0] / weak.layer_absorptance()[0][0],
        2.0,
        2e-14,
    );
    let mut heat = [0.0];
    let report = weak.exchange([1e20, 0.0], &mut heat).unwrap();
    assert!(heat[0] > 1.0);
    assert_eq!(heat[0], report.absorbed_w_m2);
    let grazing = beam(&[0.0], &[1.5], 1e-20, 0.5);
    assert_eq!(grazing.reflectances(), [1.0; 2]);
    assert!(grazing.transmittance() > 0.0);
    let matched = beam(&[0.0], &[1.0], f64::from_bits(1), 0.5);
    assert_eq!(matched.transmittance(), 1.0);
    assert_eq!(matched.reflectances(), [0.0; 2]);
}

#[test]
fn g0_beam_admission_cancellation_and_heat_publication_are_transactional() {
    for (mu, s) in [
        (0.0, 0.5),
        (1.1, 0.5),
        (f64::NAN, 0.5),
        (0.5, -0.1),
        (0.5, f64::NAN),
    ] {
        assert!(SlabBeamResponse::new(layers(&[0.1], &[1.0]), mu, s, || false).is_err());
    }
    assert!(SlabBeamResponse::new(vec![], 1.0, 0.5, || false).is_err());
    assert!(SlabBeamResponse::new(layers(&[0.1; 65], &[1.0; 65]), 1.0, 0.5, || false).is_err());
    assert!(SlabBeamResponse::new(layers(&[0.1], &[1.0]), 1.0, 0.5, || true).is_err());
    let response = beam(&[100.0], &[1.0], 1.0, 0.5);
    let mut heat = [17.0];
    assert!(response.exchange([-1.0, 0.0], &mut heat).is_err());
    assert_eq!(heat, [17.0]);
    assert!(response.exchange([f64::MAX; 2], &mut heat).is_err());
    assert_eq!(heat, [17.0]);
    let slab = slab(&[0.1]);
    assert!(
        slab.layer_net_fluxes_with_beam([0.0; 2], &[0.0], &response, [1.0, 0.0], &mut heat)
            .is_err()
    );
    assert_eq!(heat, [17.0]);
}
