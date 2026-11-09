use super::*;
use fs_material::laminate::thermal::beam::thermoelastic::LaminateBeamThermoelasticExchange;

fn values(exchange: LaminateBeamThermoelasticExchange) -> [f64; 5] {
    [
        exchange.mechanical_work_j_m,
        exchange.free_energy_change_j_m,
        exchange.entropy_change_j_m_k,
        exchange.internal_energy_change_j_m,
        exchange.reversible_heat_j_m,
    ]
}

#[test]
fn finite_uniform_layer_exchange_matches_the_uniaxial_analytic_potential() {
    let mut p = ply();
    p.e[0][1] = 0.0;
    p.nu[1] = 0.0;
    p.alpha[0][1] = 0.0;
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let old = prepare(&source, &laws, &[315.0], 0.0);
    let new = prepare(&source, &laws, &[345.0], 0.0);
    let before = [0.002, 0.03];
    let after = [0.003, -0.02];
    let exchange = old.layer_exchange(&new, 0, before, after).unwrap();
    let (e, alpha, area, inertia) = (
        p.e[0][0],
        p.alpha[0][0],
        WIDTH * p.h,
        WIDTH * p.h.powi(3) / 12.0,
    );
    let x0 = before[0] - alpha * 15.0;
    let x1 = after[0] - alpha * 45.0;
    let work = 0.5 * e * area * (x0 + x1) * (after[0] - before[0])
        + 0.5 * e * inertia * (before[1] + after[1]) * (after[1] - before[1]);
    let dpsi = 0.5 * e * area * (x1 + x0) * (x1 - x0)
        + 0.5 * e * inertia * (after[1] + before[1]) * (after[1] - before[1]);
    let ds = e * area * alpha * (x1 - x0);
    let du = dpsi + e * area * alpha * (345.0 * x1 - 315.0 * x0);
    for (actual, expected) in values(exchange)
        .into_iter()
        .zip([work, dpsi, ds, du, work - du])
    {
        close(actual, expected, 1e-11, 3e-13);
    }
    assert!(exchange.reversible_heat_j_m < 0.0);
    let reverse = new.layer_exchange(&old, 0, after, before).unwrap();
    for (forward, backward) in values(exchange).into_iter().zip(values(reverse)) {
        close(forward, -backward, 1e-12, 2e-13);
    }
    assert!(exchange.free_energy_balance_residual_j_m.abs() < 1e-11);
}

#[test]
fn bimetal_exchange_keeps_self_stress_and_original_layer_temperatures() {
    let first = ply();
    let mut second = ply();
    second.h = 0.004;
    second.isotropic = false;
    second.angle = core::f64::consts::FRAC_PI_2;
    second.e = [[105e9, 25e6], [55e9, -15e6]];
    second.alpha = [[4e-5, -1e-7], [5e-5, 1e-7]];
    let source = source(&[first, second], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let old = prepare(&source, &laws, &[300.0, 300.0], 0.0);
    let new = prepare(&source, &laws, &[315.0, 345.0], 0.0);
    let before = old.beam().free_generalized_strain();
    let after = new.beam().free_generalized_strain();
    let state = new.state(after[0], after[1], || false).unwrap();
    assert!(state.free_energy_j_m > 0.0);
    close(state.axial_force_n, 0.0, 1e-8, 0.0);
    close(state.bending_moment_n_m, 0.0, 1e-10, 0.0);
    let offset = 0.0015;
    let shifted_old = prepare(&source, &laws, &[300.0, 300.0], offset);
    let shifted_new = prepare(&source, &laws, &[315.0, 345.0], offset);
    let mut total_work = 0.0;
    let mut total_heat = 0.0;
    for layer in 0..2 {
        let exchange = old.layer_exchange(&new, layer, before, after).unwrap();
        let start = old.layer_state(layer, before[0], before[1]).unwrap();
        let end = new.layer_state(layer, after[0], after[1]).unwrap();
        for (actual, expected) in [
            (
                exchange.free_energy_change_j_m,
                end.free_energy_j_m - start.free_energy_j_m,
            ),
            (
                exchange.entropy_change_j_m_k,
                end.entropy_j_m_k - start.entropy_j_m_k,
            ),
            (
                exchange.internal_energy_change_j_m,
                end.internal_energy_j_m - start.internal_energy_j_m,
            ),
            (
                exchange.mechanical_work_j_m,
                0.5 * (start.axial_force_n + end.axial_force_n) * (after[0] - before[0])
                    + 0.5
                        * (start.bending_moment_n_m + end.bending_moment_n_m)
                        * (after[1] - before[1]),
            ),
        ] {
            close(actual, expected, 1e-11, 2e-12);
        }
        let shifted = shifted_old
            .layer_exchange(
                &shifted_new,
                layer,
                [before[0] - offset * before[1], before[1]],
                [after[0] - offset * after[1], after[1]],
            )
            .unwrap();
        for (a, b) in values(exchange).into_iter().zip(values(shifted)) {
            close(a, b, 1e-11, 2e-12);
        }
        total_work += exchange.mechanical_work_j_m;
        total_heat += exchange.reversible_heat_j_m;
    }
    close(total_work, 0.0, 1e-11, 0.0);
    close(total_heat, -state.internal_energy_j_m, 1e-10, 2e-12);
}

#[test]
fn tiny_isothermal_increment_preserves_work_below_stored_energy_roundoff() {
    let mut p = ply();
    p.e[0][1] = 0.0;
    p.nu = [0.0; 2];
    p.alpha[0] = [0.0; 2];
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let state = prepare(&source, &laws, &[320.0], 0.0);
    let eps: f64 = 0.01;
    let after = eps.next_up();
    let exchange = state
        .layer_exchange(&state, 0, [eps, 0.0], [after, 0.0])
        .unwrap();
    let expected = WIDTH * p.h * p.e[0][0] * (0.5 * eps + 0.5 * after) * (after - eps);
    assert!(expected > 0.0);
    close(exchange.free_energy_change_j_m, expected, 0.0, 2e-14);
    close(exchange.mechanical_work_j_m, expected, 0.0, 2e-14);
    close(exchange.internal_energy_change_j_m, expected, 0.0, 2e-14);
    close(exchange.reversible_heat_j_m, 0.0, 1e-28, 0.0);
    assert_eq!(
        state
            .layer_exchange(&state, 0, [eps, 0.0], [eps, 0.0])
            .unwrap(),
        LaminateBeamThermoelasticExchange::default(),
    );
}

#[test]
fn tiny_increment_under_thermal_prestress_retains_entropy_and_reversible_heat() {
    let mut p = ply();
    p.e[0][1] = 0.0;
    p.nu = [0.28, 0.0];
    p.alpha[0] = [1e-4, 0.0];
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let t = 350.0;
    let state = prepare(&source, &laws, &[t], 0.0);
    let eps: f64 = 0.0001;
    let after = eps.next_up();
    let free = laws[0].at_temperature(t).unwrap().free_strain()[0];
    // The generalized increment is real but subtraction of thermal strain
    // erases it in the endpoint elastic-coordinate representation.
    assert_eq!((eps - free).to_bits(), (after - free).to_bits());
    let exchange = state
        .layer_exchange(&state, 0, [eps, 0.0], [after, 0.0])
        .unwrap();
    let ea = WIDTH * p.h * p.e[0][0];
    let work = ea * (eps - free) * (after - eps);
    let entropy = ea * p.alpha[0][0] * (after - eps);
    close(exchange.mechanical_work_j_m, work, 0.0, 2e-13);
    close(exchange.free_energy_change_j_m, work, 0.0, 2e-13);
    close(exchange.entropy_change_j_m_k, entropy, 0.0, 2e-13);
    close(
        exchange.internal_energy_change_j_m,
        work + t * entropy,
        0.0,
        2e-13,
    );
    close(exchange.reversible_heat_j_m, -t * entropy, 0.0, 2e-13);
    assert!(exchange.reversible_heat_j_m < 0.0);
}

#[test]
fn finite_exchange_refines_to_the_same_relaxed_thermoelastic_heat_rate() {
    let mut p = ply();
    p.isotropic = false;
    p.angle = core::f64::consts::FRAC_PI_2;
    let source = source(&[p], ReferenceSurface::Midsurface);
    let laws = laws(&source);
    let t = 322.0;
    let before = [0.003, 0.08];
    let rates = [2e-4, -0.03];
    let td = 4.0;
    let old = prepare(&source, &laws, &[t], 0.001);
    let at = old.layer_state(0, before[0], before[1]).unwrap();
    let rate = t
        * (at.generalized_force_temperature_derivative[0] * rates[0]
            + at.generalized_force_temperature_derivative[1] * rates[1])
        - at.heat_capacity_j_m_k * td;
    let mut errors = Vec::new();
    for dt in [1e-2, 5e-3, 2.5e-3] {
        let new = prepare(&source, &laws, &[t + dt * td], 0.001);
        let after = [before[0] + dt * rates[0], before[1] + dt * rates[1]];
        let heat = old
            .layer_exchange(&new, 0, before, after)
            .unwrap()
            .reversible_heat_j_m;
        errors.push((heat / dt - rate).abs());
    }
    assert!(errors[1] < 0.6 * errors[0]);
    assert!(errors[2] < 0.6 * errors[1]);
}

#[test]
fn exchange_refuses_foreign_sources_geometry_intervals_and_invalid_states() {
    let original = source(&[ply()], ReferenceSurface::Midsurface);
    let original_laws = laws(&original);
    let old = prepare(&original, &original_laws, &[320.0], 0.0);
    let end = prepare(&original.clone(), &original_laws, &[330.0], 0.0);
    let strain = [0.002, 0.0];
    assert!(old.layer_exchange(&end, 0, strain, strain).is_ok());
    let foreign = source(&[ply()], ReferenceSurface::Midsurface);
    let foreign = prepare(&foreign, &laws(&foreign), &[330.0], 0.0);
    assert!(old.layer_exchange(&foreign, 0, strain, strain).is_err());
    let shifted = prepare(&original, &original_laws, &[330.0], 0.001);
    assert!(old.layer_exchange(&shifted, 0, strain, strain).is_err());
    let changed_laws =
        vec![AffineLaminaThermoelastic::new(&original.layers()[0], [290.0, 350.0], 0.02).unwrap()];
    let changed = prepare(&original, &changed_laws, &[330.0], 0.0);
    assert!(old.layer_exchange(&changed, 0, strain, strain).is_err());
    let wider = PreparedLaminateBeam::try_new(
        Arc::clone(end.source().prepared()),
        2.0 * WIDTH,
        0.0,
        || false,
    )
    .unwrap();
    let wider =
        ResolvedLaminateBeamThermoelastic::try_new(&original_laws, end.source(), &wider, || false)
            .unwrap();
    assert!(old.layer_exchange(&wider, 0, strain, strain).is_err());
    assert!(old.layer_exchange(&end, 1, strain, strain).is_err());
    assert!(old.layer_exchange(&end, 0, strain, [0.03, 0.0]).is_err());
    assert!(
        old.layer_exchange(&end, 0, [f64::NAN, 0.0], strain)
            .is_err()
    );
}
