//! Diagonal square orthotropy stays planar without becoming isotropic.
use super::*;

fn square_ply(angle: f64) -> Ply {
    Ply {
        isotropic: false,
        angle,
        e: [[70e9, -30e6]; 2],
        nu: [0.28, 0.0004],
        g: [10e9, -5e6],
        alpha: [[2e-5, 1e-7]; 2],
        ..ply()
    }
}

// Independent uniaxial compliance at 45 degrees:
// 1/E45=(1-nu)/(2E)+1/(4G). Differentiate that rational expression,
// not the native rotated stiffness, Schur reduction or thermal response.
fn analytic(p: Ply, t: f64, strain: [f64; 2]) -> [f64; 8] {
    let dt = t - 300.0;
    let [er, et] = p.e[0];
    let [gr, gt] = p.g;
    let [nr, nt] = p.nu;
    let e = er + et * dt;
    let g = gr + gt * dt;
    let nu = nr + nt * dt;
    let numerator = 4.0 * e * g;
    let numerator_t = 4.0 * (et * g + e * gt);
    let numerator_tt = 8.0 * et * gt;
    let denominator = e + 2.0 * g * (1.0 - nu);
    let denominator_t = et + 2.0 * gt * (1.0 - nu) - 2.0 * g * nt;
    let denominator_tt = -4.0 * gt * nt;
    let axial = numerator / denominator;
    let axial_t = (numerator_t - axial * denominator_t) / denominator;
    let axial_tt =
        (numerator_tt - axial * denominator_tt - 2.0 * axial_t * denominator_t) / denominator;
    let [ar, at] = p.alpha[0];
    let alpha = ar + at * dt;
    let free = (ar + 0.5 * at * dt) * dt;
    let area = WIDTH * p.h;
    let inertia = WIDTH * p.h.powi(3) / 12.0;
    let ex = strain[0] - free;
    let square = area * ex * ex + inertia * strain[1] * strain[1];
    let psi = 0.5 * axial * square;
    let entropy = axial * alpha * area * ex - 0.5 * axial_t * square;
    let capacity = -t
        * (0.5 * axial_tt * square - 2.0 * alpha * axial_t * area * ex
            + axial * alpha * alpha * area
            - at * axial * area * ex);
    [
        axial * area * ex,
        axial * inertia * strain[1],
        psi,
        entropy,
        psi + t * entropy,
        capacity,
        area * (axial_t * ex - axial * alpha),
        axial_t * inertia * strain[1],
    ]
}

#[test]
fn square_diagonals_use_transformed_axial_stiffness_entropy_and_capacity() {
    let strain = [0.003, 0.06];
    for angle in [core::f64::consts::FRAC_PI_4, -core::f64::consts::FRAC_PI_4] {
        let p = square_ply(angle);
        let source = source(&[p], ReferenceSurface::Midsurface);
        let laws = laws(&source);
        for t in [285.0, 332.0, 355.0] {
            let prepared = prepare(&source, &laws, &[t], 0.0);
            let state = prepared.layer_state(0, strain[0], strain[1]).unwrap();
            let expected = analytic(p, t, strain);
            let actual = [
                state.axial_force_n,
                state.bending_moment_n_m,
                state.free_energy_j_m,
                state.entropy_j_m_k,
                state.internal_energy_j_m,
                state.heat_capacity_j_m_k,
                state.generalized_force_temperature_derivative[0],
                state.generalized_force_temperature_derivative[1],
            ];
            for (a, b) in actual.into_iter().zip(expected) {
                close(a, b, 1e-10, 3e-12);
            }
            let original_e = p.e[0][0] + p.e[0][1] * (t - 300.0);
            let axial = prepared.beam().tangent()[0][0] / (WIDTH * p.h);
            assert!((axial - original_e).abs() > 0.25 * original_e);
            assert!(
                prepared.layer_minimum_heat_capacity_j_m_k(0).unwrap() <= state.heat_capacity_j_m_k
            );
            close(
                state.free_energy_j_m,
                prepared
                    .beam()
                    .evaluate(strain[0], strain[1], || false)
                    .unwrap()
                    .strain_energy_j_m,
                1e-11,
                2e-13,
            );
        }
    }
}

#[test]
fn square_layer_finite_exchange_reuses_the_original_rotated_potential() {
    let p = square_ply(core::f64::consts::FRAC_PI_4);
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let (t0, t1) = (312.0, 345.0);
    let before = [0.0025, -0.04];
    let after = [0.0015, 0.08];
    let old = prepare(&source, &laws, &[t0], 0.0);
    let new = prepare(&source, &laws, &[t1], 0.0);
    let exchange = old.layer_exchange(&new, 0, before, after).unwrap();
    let a = analytic(p, t0, before);
    let b = analytic(p, t1, after);
    let work =
        0.5 * (a[0] + b[0]) * (after[0] - before[0]) + 0.5 * (a[1] + b[1]) * (after[1] - before[1]);
    for (actual, expected) in [
        (exchange.mechanical_work_j_m, work),
        (exchange.free_energy_change_j_m, b[2] - a[2]),
        (exchange.entropy_change_j_m_k, b[3] - a[3]),
        (exchange.internal_energy_change_j_m, b[4] - a[4]),
        (exchange.reversible_heat_j_m, work - (b[4] - a[4])),
    ] {
        close(actual, expected, 1e-10, 3e-12);
    }
    assert!(exchange.free_energy_balance_residual_j_m.abs() < 1e-10);
}

#[test]
fn symmetry_requires_original_profiles_not_current_or_roundoff_close_values() {
    for asymmetry in 0..3 {
        let mut p = square_ply(core::f64::consts::FRAC_PI_4);
        match asymmetry {
            0 => p.e[1][1] = 20e6,
            1 => p.alpha[1][1] = -2e-7,
            _ => p.e[1][0] = p.e[0][0].next_up(),
        }
        let source = source(&[p], ReferenceSurface::Midsurface);
        let laws = laws(&source);
        let current = source.at_temperatures(&[300.0], || false).unwrap();
        // The current mechanical state looks aligned (including its admitted
        // rotation roundoff); that does not establish the thermal source family.
        let beam =
            PreparedLaminateBeam::try_new(Arc::clone(current.prepared()), WIDTH, 0.0, || false)
                .unwrap();
        assert!(
            ResolvedLaminateBeamThermoelastic::try_new(&laws, &current, &beam, || false,).is_err()
        );
    }
    let p = square_ply(core::f64::consts::FRAC_PI_6);
    let source = source(&[p], ReferenceSurface::Midsurface);
    let current = source.at_temperatures(&[330.0], || false).unwrap();
    assert!(
        PreparedLaminateBeam::try_new(Arc::clone(current.prepared()), WIDTH, 0.0, || false,)
            .is_err()
    );
}
