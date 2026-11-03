use super::*;

fn observe(
    film: &ShearedFilm,
    q: &[f64; 3],
    v: &[f64; 3],
) -> ([f64; 3], f64, Vec<ShearChannelReport>, Vec<f64>, f64) {
    let mut force = [0.0; 3];
    let mut pressure = [0.0];
    let mut channels = vec![ShearChannelReport::default(); 2];
    let mut work = vec![0.0; film.layer_count()];
    let report = film
        .evaluate_with_layers_into(q, v, &mut force, &mut pressure, &mut channels, &mut work)
        .unwrap();
    (
        force,
        pressure[0],
        channels,
        work,
        report.dissipated_power_w,
    )
}

#[test]
fn viscosity_refresh_changes_mobility_entrainment_drag_and_layer_heat_together() {
    // G1: independent two-layer Couette integrals. Oppositely oriented equal
    // passages give zero pressure, so no pressure work can hide a shear error.
    let original = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H, H]), vec![layers(); 2]).unwrap(),
        vec![walls(1.0), walls(-1.0)],
    )
    .unwrap();
    let mu = [0.025, 0.005];
    let next = original
        .with_layer_viscosities(&[mu[0], mu[1], mu[0], mu[1]])
        .unwrap();
    let v = [0.0, 0.013, -0.02];
    let (force, pressure, channels, work, power) = observe(&next, &[0.0; 3], &v);
    let a = 0.5 / mu[0] + 0.5 / mu[1];
    let first = 0.125 / mu[0] + 0.375 / mu[1];
    let second = (0.125 / mu[0] + 0.875 / mu[1]) / 3.0;
    let neutral = first / a;
    let mobility = second - first * first / a;
    let shear = W * L / (H * a);
    let relative = v[2] - v[1];
    let entrainment = W * H * (neutral * v[1] + (1.0 - neutral) * v[2]);
    near(pressure, 0.0, 1e-12);
    near(
        next.pressure_operator().channel_viscosities_pa_s()[0],
        1.0 / (12.0 * mobility),
        1e-15,
    );
    near(channels[0].entrained_flow_m3_s, entrainment, 1e-20);
    near(channels[1].entrained_flow_m3_s, -entrainment, 1e-20);
    near(force[1], -2.0 * shear * relative, 1e-16);
    near(force[2], 2.0 * shear * relative, 1e-16);
    near(power, 2.0 * shear * relative * relative, 1e-17);
    for (index, &mu) in mu.iter().enumerate() {
        let expected = W * L * relative * relative * 0.5 / (H * a * a * mu);
        near(work[index], expected, 1e-17);
        near(work[2 + index], expected, 1e-17);
    }
    let old = observe(&original, &[0.0; 3], &v);
    assert_ne!(old.0, force);
    assert_ne!(
        old.2[0].entrained_flow_m3_s,
        channels[0].entrained_flow_m3_s
    );
    assert_ne!(old.3, work);
    near(dot(&force, &v), power, 1e-17);
    assert_eq!(
        original.stratified_profiles().unwrap()[0].layers()[0].dynamic_viscosity_pa_s,
        0.01
    );
}

#[test]
fn viscosity_refresh_preserves_barus_graph_walls_and_complete_tangent() {
    let graph = film(&[H, 0.8 * H]).with_barus_viscosity(2e-5).unwrap();
    let original = ShearedFilm::from_stratified(
        StratifiedFilm::new(&graph, vec![layers(); 2]).unwrap(),
        vec![walls(1.0), walls(-1.0)],
    )
    .unwrap();
    let next = original
        .with_layer_viscosities(&[0.02, 0.005, 0.03, 0.012])
        .unwrap();
    assert!(std::sync::Arc::ptr_eq(
        &original.film.cells,
        &next.film.cells
    ));
    assert!(std::sync::Arc::ptr_eq(
        &original.film.channels,
        &next.film.channels
    ));
    assert_eq!(next.wall_ports(), original.wall_ports());
    assert_eq!(next.film.pressure_viscosity_coefficient_per_pa(), 2e-5);
    for profile in next.stratified_profiles().unwrap() {
        assert_eq!(profile.layers()[0].upper_fraction, 0.5);
        assert_eq!(profile.layers()[1].upper_fraction, 1.0);
    }
    let q = [0.01 * H, 0.0, 0.0];
    let v = [0.0001, 0.013, -0.02];
    let dq = [1e-5, 0.0, 0.0];
    let dv = [0.0003, -0.01, 0.006];
    let base = observe(&next, &q, &v);
    near(dot(&base.0, &v), base.4, 1e-15);
    let mut tangent = [0.0; 3];
    next.tangent_into(&q, &v, &dq, &dv, &mut tangent).unwrap();
    let epsilon = 1e-5;
    let sample = |sign: f64| {
        observe(
            &next,
            &core::array::from_fn(|i| q[i] + sign * epsilon * dq[i]),
            &core::array::from_fn(|i| v[i] + sign * epsilon * dv[i]),
        )
        .0
    };
    let plus = sample(1.0);
    let minus = sample(-1.0);
    for i in 0..3 {
        let expected = (plus[i] - minus[i]) / (2.0 * epsilon);
        near(tangent[i], expected, 2e-8 * (1.0 + expected.abs()));
    }
}

#[test]
fn uniform_layer_limit_and_unchanged_viscosity_preserve_native_mechanics() {
    let uniform = ShearedFilm::new(film(&[H, H]), vec![walls(1.0), walls(-1.0)]).unwrap();
    let layered = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H, H]), vec![layers(); 2]).unwrap(),
        vec![walls(1.0), walls(-1.0)],
    )
    .unwrap();
    let q = [0.0; 3];
    let v = [0.0001, 0.013, -0.02];
    let old = observe(&uniform, &q, &v);
    assert_eq!(
        observe(&uniform.with_layer_viscosities(&[MU; 2]).unwrap(), &q, &v),
        old
    );
    let next = layered.with_layer_viscosities(&[MU; 4]).unwrap();
    let expanded = observe(&next, &q, &v);
    assert_eq!(
        (&expanded.0, expanded.1, &expanded.2, expanded.4),
        (&old.0, old.1, &old.2, old.4)
    );
    for (channel, range) in [0..2, 2..4].into_iter().enumerate() {
        near(
            expanded.3[range].iter().sum(),
            old.2[channel].dissipated_power_w(),
            1e-17,
        );
    }
    assert_eq!(
        observe(&next.with_layer_viscosities(&[MU; 4]).unwrap(), &q, &v),
        expanded
    );
    let doubled = uniform.with_layer_viscosities(&[2.0 * MU; 2]).unwrap();
    let expected = ShearedFilm::new(
        film(&[H, H]).with_viscosity(2.0 * MU).unwrap(),
        vec![walls(1.0), walls(-1.0)],
    )
    .unwrap();
    assert_eq!(observe(&doubled, &q, &v), observe(&expected, &q, &v));
}

#[test]
fn malformed_viscosity_updates_leave_the_original_operator_usable() {
    let original = ShearedFilm::from_stratified(
        StratifiedFilm::new(&film(&[H, H]), vec![layers(); 2]).unwrap(),
        vec![walls(1.0), walls(-1.0)],
    )
    .unwrap();
    let q = [0.0; 3];
    let v = [0.0001, 0.013, -0.02];
    let before = observe(&original, &q, &v);
    for values in [
        vec![],
        vec![MU; 3],
        vec![MU; 5],
        vec![MU, 0.0, MU, MU],
        vec![MU, -MU, MU, MU],
        vec![MU, f64::NAN, MU, MU],
        vec![MU, f64::INFINITY, MU, MU],
    ] {
        assert!(original.with_layer_viscosities(&values).is_err());
        assert_eq!(observe(&original, &q, &v), before);
    }
}
