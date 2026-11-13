//! G0/G1/G3/G4: finite inventories, independent kinetic moments and full windows.
use super::*;
use crate::molecular_exchange::network::trajectory::MolecularNetworkStreamMoments;
use crate::molecular_exchange::network::trajectory::{
    MolecularExchangeAdaptiveControls, MolecularNetworkAdaptiveControls,
    MolecularNetworkTrajectoryCandidate, MolecularNetworkTrajectoryError,
    advance_candidate,
};

fn adaptive(h: f64) -> MolecularNetworkAdaptiveControls {
    MolecularNetworkAdaptiveControls {
        states: MolecularExchangeAdaptiveControls {
            initial_step_s: h,
            minimum_step_s: 1e-7,
            maximum_step_s: h,
            maximum_attempts: 20_000,
            maximum_accepted_intervals: 10_000,
            relative_tolerance: 0.0,
            temperature_tolerance_k: 1e-3,
            pressure_tolerance_pa: 0.05,
            species_mass_tolerance_kg: 1e-10,
            velocity_tolerance_m_s: 1e-3,
        },
        stream_momentum_tolerance_n_s: 1e-7,
        stream_energy_tolerance_j: 1e-4,
        window_mass_tolerance_kg: 1e-10,
        window_momentum_tolerance_n_s: 1e-9,
        window_energy_tolerance_j: 1e-3,
    }
}
fn fixed(h: f64) -> MolecularNetworkAdaptiveControls {
    let mut controls = adaptive(h);
    controls.states.temperature_tolerance_k = 1e6;
    controls.states.pressure_tolerance_pa = 1e12;
    controls.states.species_mass_tolerance_kg = 1.0;
    controls.states.velocity_tolerance_m_s = 1e6;
    controls.stream_momentum_tolerance_n_s = 1e6;
    controls.stream_energy_tolerance_j = 1e12;
    controls
}
fn evolve(
    gases: &BTreeMap<u64, MovingGasParcel>,
    links: &[MolecularNetworkLink<'_>],
    duration: f64,
    physical: MolecularNetworkControls,
    adaptive: MolecularNetworkAdaptiveControls,
) -> MolecularNetworkTrajectoryCandidate {
    advance_candidate(gases, links, duration, &solves(gases), physical, adaptive, || false).unwrap()
}
fn ports<'a>(law: &'a EffusiveDriftInterface, boost: [f64; 3]) -> [MolecularNetworkLink<'a>; 3] {
    let mut links = [link(12, [2, 4], law), link(7, [4, 8], law), link(22, [8, 2], law)];
    links[1].surface.normal_left_to_right = [0.0, 1.0, 0.0];
    links[2].surface.normal_left_to_right = [0.0, 0.0, 1.0];
    for (link, velocity) in links.iter_mut().zip([[2.0, -4.0, 1.0], [-3.0, 5.0, 2.0], [1.0, 3.0, -2.0]]) {
        link.surface.velocity_world_m_s = core::array::from_fn(|d| velocity[d] + boost[d]);
    }
    links
}

#[test]
fn new_species_cross_the_network_on_accepted_states_and_history_replays_exactly() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let only_a = model.with_composition(
        Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap()).unwrap();
    let only_b = model.with_composition(
        Composition::new(CompositionBasis::MoleFraction, vec![0.0, 1.0]).unwrap()).unwrap();
    let gases = BTreeMap::from([
        (2, gas(only_a.clone(), T, P, [0.0; 3])),
        (4, gas(only_a, 500.0, 1.2 * P, [0.0, 10.0, 0.0])),
        (8, gas(only_b, 600.0, 1.3 * P, [5.0, -10.0, 0.0])),
    ]);
    let saved = gases.clone();
    let mut links = [link(1, [2, 4], &law), link(2, [4, 8], &law)];
    links[1].surface.normal_left_to_right = [0.0, 1.0, 0.0];
    links[1].surface.velocity_world_m_s = [1.0, 2.0, 3.0];
    let duration = 0.25;
    let physical = budgets([7.0, -2.0, 4.0]);
    let result = evolve(&gases, &links, duration, physical, adaptive(duration));
    let a = SpeciesId::new("A").unwrap();
    let b = SpeciesId::new("B").unwrap();
    assert!(result.step().accepted_intervals > 1);
    assert!(result.step().rejected_intervals > 0);
    assert_eq!(result.step().accepted_steps.len(), 2 * result.step().accepted_intervals);
    assert_eq!(result.step().attempted_intervals,
        result.step().accepted_intervals + result.step().rejected_intervals);
    assert_eq!(result.step().duration_s, duration);
    assert!(result.step().maximum_accepted_error <= 1.0);
    let first = &result.step().accepted_steps[0].links[0];
    assert_eq!(first.id, 1);
    assert!(!first.exchange.species.iter().any(|row| row.species == b));
    assert!(result.step().accepted_steps[1].links[0].exchange.species.iter()
        .any(|row| row.species == b && row.mass_kg[1] > 0.0));
    assert!(result.gases()[&2].parcel().species_mass_kg(&b).unwrap() > 0.0);
    for species in [&a, &b] {
        near(result.gases().values().map(|g| g.parcel().species_mass_kg(species).unwrap()).sum(),
            gases.values().map(|g| g.parcel().species_mass_kg(species).unwrap()).sum(), 1e-13);
    }
    let mut replay = gases.clone();
    for step in &result.step().accepted_steps {
        let native = network::step_candidate(&replay, &links, step.duration_s,
            &solves(&replay), physical, || false).unwrap();
        assert_eq!(native.step(), step);
        replay = native.into_parts().0;
    }
    assert_eq!(result.gases(), &replay);
    for total in &result.step().links {
        for side in 0..2 {
            let gross: f64 = result.step().accepted_steps.iter()
                .flat_map(|step| step.links.iter().filter(|link| link.id == total.id))
                .flat_map(|link| &link.exchange.species).map(|row| row.mass_kg[side]).sum();
            near(total.gross[side].mass_kg, gross, 1e-16);
        }
    }
    near(result.gases().values().map(|g| internal(g) + g.kinetic_energy_j()).sum::<f64>()
        - gases.values().map(|g| internal(g) + g.kinetic_energy_j()).sum::<f64>(), 0.0, 1e-5);
    links.reverse();
    assert_eq!(result, evolve(&gases, &links, duration, physical, adaptive(duration)));
    assert_eq!(gases, saved);
}

type OdeState = [[f64; 5]; 3]; // m, P_x, P_y, P_z, U+K in the laboratory frame.

fn ode_state(gases: &BTreeMap<u64, MovingGasParcel>) -> OdeState {
    [2, 4, 8].map(|id| {
        let g = &gases[&id];
        let p = g.momentum_kg_m_s();
        [g.parcel().mass_kg(), p[0], p[1], p[2], internal(g) + g.kinetic_energy_j()]
    })
}
fn ode_rates(state: OdeState, links: &[MolecularNetworkLink<'_>]) -> OdeState {
    let mut rates = [[0.0; 5]; 3];
    let ids = [2, 4, 8];
    for link in links {
        let endpoints = link.gases.map(|id| ids.iter().position(|&x| x == id).unwrap());
        for side in 0..2 {
            let donor = endpoints[side];
            let receiver = endpoints[1 - side];
            let y = state[donor];
            let velocity: [f64; 3] = core::array::from_fn(|d| y[1 + d] / y[0]);
            let kinetic = 0.5 * y[0] * velocity.iter().map(|v| v * v).sum::<f64>();
            let temperature = (y[4] - kinetic) / (y[0] * 1.5 * R / M);
            let sigma = (R * temperature / M).sqrt();
            let n = link.surface.normal_left_to_right.map(|v| if side == 0 { v } else { -v });
            let drift: [f64; 3] = core::array::from_fn(|d| velocity[d] - link.surface.velocity_world_m_s[d]);
            let normal = (0..3).map(|d| drift[d] * n[d]).sum::<f64>();
            let tangent: [f64; 3] = core::array::from_fn(|d| drift[d] - normal * n[d]);
            let z = normal / sigma;
            let phi = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
            let cdf = 0.5 * (1.0 + fs_math::det::erf(z / 2.0_f64.sqrt()));
            // Independent positive-half Gaussian raw moments E[w^k 1(w>0)].
            let j1 = sigma * phi + normal * cdf;
            let j2 = (normal * normal + sigma * sigma) * cdf + normal * sigma * phi;
            let j3 = (normal.powi(3) + 3.0 * normal * sigma * sigma) * cdf
                + (normal * normal * sigma + 2.0 * sigma.powi(3)) * phi;
            let weight = y[0] / VOLUME * link.surface.area_m2;
            let mass = weight * j1;
            let momentum: [f64; 3] = core::array::from_fn(|d| weight * (j2 * n[d] + j1 * tangent[d]));
            // Monatomic molecules have only translational thermal energy.
            let mut energy = 0.5 * weight * (j3 + j1 *
                (tangent.iter().map(|v| v * v).sum::<f64>() + 2.0 * sigma * sigma));
            let mut world = [mass, 0.0, 0.0, 0.0, 0.0];
            for d in 0..3 {
                let v = link.surface.velocity_world_m_s[d];
                world[1 + d] = momentum[d] + v * mass;
                energy += v * momentum[d] + 0.5 * mass * v * v;
            }
            world[4] = energy;
            for k in 0..5 {
                rates[donor][k] -= world[k];
                rates[receiver][k] += world[k];
            }
        }
    }
    rates
}
fn rk4(mut state: OdeState, links: &[MolecularNetworkLink<'_>], duration: f64, count: usize) -> OdeState {
    let h = duration / count as f64;
    let shifted = |y: OdeState, k: OdeState, fraction: f64| core::array::from_fn(|i|
        core::array::from_fn(|j| y[i][j] + fraction * h * k[i][j]));
    for _ in 0..count {
        let a = ode_rates(state, links);
        let b = ode_rates(shifted(state, a, 0.5), links);
        let c = ode_rates(shifted(state, b, 0.5), links);
        let d = ode_rates(shifted(state, c, 1.0), links);
        for i in 0..3 { for j in 0..5 {
            state[i][j] += h * (a[i][j] + 2.0 * b[i][j] + 2.0 * c[i][j] + d[i][j]) / 6.0;
        } }
    }
    state
}
fn ode_distance(actual: OdeState, reference: OdeState, original: OdeState) -> f64 {
    (0..3).map(|i| {
        let mass = original[i][0];
        (actual[i][0] - reference[i][0]).abs() / mass
            + (1..4).map(|d| (actual[i][d] - reference[i][d]).abs()
                / (mass * (R * T / M).sqrt())).sum::<f64>()
            + (actual[i][4] - reference[i][4]).abs() / (mass * 1.5 * R * T / M)
    }).sum()
}

#[test]
fn network_refinement_converges_to_independent_continuous_gaussian_moment_equations() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = three(&model, [0.0; 3]);
    let links = ports(&law, [0.0; 3]);
    let initial = ode_state(&gases);
    let duration = 0.8;
    let reference = rk4(initial, &links, duration, 1024);
    assert!(ode_distance(rk4(initial, &links, duration, 512), reference, initial) < 1e-10);
    let errors = [0.4, 0.2, 0.1].map(|h| {
        let result = evolve(&gases, &links, duration, budgets([7.0, -2.0, 4.0]), fixed(h));
        ode_distance(ode_state(result.gases()), reference, initial)
    });
    assert!(errors[1] < 0.65 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.65 * errors[1], "{errors:?}");
    let mut accuracy = adaptive(0.4);
    accuracy.states.temperature_tolerance_k = 1e-5;
    accuracy.states.velocity_tolerance_m_s = 1e-5;
    accuracy.states.species_mass_tolerance_kg = 1e-11;
    accuracy.stream_momentum_tolerance_n_s = 1e-8;
    accuracy.stream_energy_tolerance_j = 1e-5;
    let result = evolve(&gases, &links, duration, budgets([7.0, -2.0, 4.0]), accuracy);
    assert!(ode_distance(ode_state(result.gases()), reference, initial) < 0.5 * errors[2]);
    assert!(result.step().maximum_accepted_error <= 1.0);
    // This reference evolves only the aperture operator; no remaining wall
    // pressure forces or isolated-vessel equilibrium are implied.
}

#[test]
fn gross_stream_moments_can_drive_refinement_independently_of_loose_gas_state_tolerances() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = three(&model, [0.0; 3]);
    let links = ports(&law, [0.0; 3]);
    let duration = 0.5;
    let loose = evolve(&gases, &links, duration, budgets([0.0; 3]), fixed(duration));
    let mut strict = fixed(duration);
    strict.stream_momentum_tolerance_n_s = 1e-9;
    strict.stream_energy_tolerance_j = 1e-6;
    let refined = evolve(&gases, &links, duration, budgets([0.0; 3]), strict);
    assert_eq!(loose.step().accepted_intervals, 1);
    assert!(refined.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(refined.step().rejected_intervals > 0);
    assert!(refined.step().maximum_accepted_error <= 1.0);
    assert!(refined.step().links.iter().all(|l| l.gross.iter().all(|g| g.mass_kg > 0.0)));
}

#[test]
fn common_frame_streams_and_thermal_states_survive_formation_shifts_and_observer_boosts() {
    let run = |boost: [f64; 3], datum| {
        let model = catalog(MIXTURE, datum);
        let law = interface(&model, 2000.0);
        let gases = three(&model, boost);
        let links = ports(&law, boost);
        evolve(&gases, &links, 0.2,
            budgets(core::array::from_fn(|d| [7.0, -2.0, 4.0][d] + boost[d])), fixed(0.05))
    };
    let rest = run([0.0; 3], 0.0);
    let boost = [500.0, -700.0, 300.0];
    let moved = run(boost, 1e12);
    assert_eq!(rest.step().accepted_intervals, moved.step().accepted_intervals);
    for (&id, gas) in rest.gases() {
        let shifted = &moved.gases()[&id];
        near(temperature(shifted), temperature(gas), 2e-7);
        near(shifted.parcel().mass_kg(), gas.parcel().mass_kg(), 1e-14);
        for d in 0..3 {
            near(shifted.velocity_m_s()[d] - boost[d], gas.velocity_m_s()[d], 2e-9);
        }
    }
    for (a, b) in rest.step().links.iter().zip(&moved.step().links) {
        for (a, b) in a.gross.iter().chain([&a.net_left_to_right])
            .zip(b.gross.iter().chain([&b.net_left_to_right])) {
            near(a.mass_kg, b.mass_kg, 1e-14);
            for d in 0..3 { near(a.momentum_n_s[d], b.momentum_n_s[d], 2e-10); }
            near(a.reference_energy_j, b.reference_energy_j, 2e-7);
        }
    }
}

#[test]
fn fixed_window_budgets_native_limits_and_late_cancellation_publish_nothing() {
    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = three(&model, [0.0; 3]);
    let saved = gases.clone();
    let links = ports(&law, [0.0; 3]);
    let duration = 0.2;
    let controls = fixed(0.025);
    let physical = budgets([7.0, -2.0, 4.0]);
    let mut polls = 0;
    let baseline = advance_candidate(&gases, &links, duration, &solves(&gases), physical, controls, || {
        polls += 1;
        false
    }).unwrap();
    let spent = baseline.step().debit.energy_j.max(baseline.step().debit.laboratory_energy_j);
    assert!(spent > 0.0);
    let mut strict = controls;
    strict.window_energy_tolerance_j = 0.5 * spent;
    assert!(matches!(advance_candidate(&gases, &links, duration, &solves(&gases), physical,
        strict, || false), Err(MolecularNetworkTrajectoryError::WindowBalance { .. })));
    let mut limited = controls;
    limited.states.maximum_accepted_intervals = 1;
    assert!(matches!(advance_candidate(&gases, &links, duration, &solves(&gases), physical,
        limited, || false), Err(MolecularNetworkTrajectoryError::WorkBudget { accepted_intervals: 1, elapsed_s, .. }) if elapsed_s > 0.0));
    let mut minimum = adaptive(duration);
    minimum.states.minimum_step_s = duration;
    minimum.states.temperature_tolerance_k = 1e-14;
    minimum.states.species_mass_tolerance_kg = 1e-20;
    minimum.states.velocity_tolerance_m_s = 1e-14;
    assert!(matches!(advance_candidate(&gases, &links, duration, &solves(&gases), physical,
        minimum, || false), Err(MolecularNetworkTrajectoryError::MinimumStep { elapsed_s: 0.0, .. })));
    let mut remaining = polls;
    assert!(matches!(advance_candidate(&gases, &links, duration, &solves(&gases), physical,
        controls, || { remaining -= 1; remaining == 0 }),
        Err(MolecularNetworkTrajectoryError::Cancelled { elapsed_s }) if elapsed_s == duration));
    let restrictive = interface(&model, 1.0);
    assert!(matches!(advance_candidate(&gases, &ports(&restrictive, [0.0; 3]), duration,
        &solves(&gases), physical, controls, || false),
        Err(MolecularNetworkTrajectoryError::Step { elapsed_s: 0.0, error: MolecularNetworkError::Link { .. }, .. })));
    assert_eq!(gases, saved);
}

#[test]
fn shared_withdrawal_retries_and_closed_apertures_keep_inactive_sources_inactive() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = three(&model, [0.0; 3]);
    let links = ports(&law, [0.0; 3]);
    let mut physical = budgets([0.0; 3]);
    physical.maximum_duration_s = 0.015;
    physical.maximum_species_fraction = 1e-4;
    let accuracy = fixed(0.04);
    let result = evolve(&gases, &links, 0.04, physical, accuracy);
    assert!(result.step().rejected_intervals > 0);
    assert!(result.step().accepted_steps.iter().all(|s| s.duration_s <= 0.0075));
    let mut work = accuracy;
    work.states.maximum_attempts = 1;
    assert!(matches!(advance_candidate(&gases, &links, 0.04, &solves(&gases), physical,
        work, || false), Err(MolecularNetworkTrajectoryError::WorkBudget {
            elapsed_s: 0.0, attempted_intervals: 1, accepted_intervals: 0 })));

    let unsupported = interface(&model, 0.0);
    let closed_gases = BTreeMap::from([
        (2, gas(model.clone(), 1900.0, P, [5000.0, 0.0, 0.0])),
        (4, gas(model, 1900.0, P, [-5000.0, 0.0, 0.0])),
    ]);
    let closed = [MolecularNetworkLink {
        id: 1, gases: [2, 4], interface: &unsupported,
        surface: MolecularExchangeSurface {
            area_m2: 0.0, normal_left_to_right: [f64::NAN; 3], velocity_world_m_s: [f64::NAN; 3],
        },
    }];
    let result = evolve(&closed_gases, &closed, 0.2, budgets([0.0; 3]), fixed(0.05));
    assert_eq!(result.gases(), &closed_gases);
    assert_eq!(result.step().residual, Default::default());
    assert_eq!(result.step().debit, Default::default());
    assert_eq!(
        result.step().links[0].gross,
        [MolecularNetworkStreamMoments::default(); 2]
    );
    let empty = evolve(&closed_gases, &[], 0.2, budgets([0.0; 3]), fixed(0.05));
    assert_eq!(empty.gases(), &closed_gases);
    assert!(empty.step().links.is_empty());
}
