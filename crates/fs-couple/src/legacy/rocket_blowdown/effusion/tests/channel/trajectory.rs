//! G1/G3/G4: adaptive channel/wall/recoil evolution, physical stream errors,
//! nonlinear conduction and fixed-window conservation with atomic refusal.

use super::*;
use crate::rocket_blowdown::effusion::channel::trajectory::{
    EffusionChannelAdaptiveSolve, EffusionChannelTrajectoryCandidate,
    EffusionChannelTrajectoryError,
};
use crate::rocket_blowdown::effusion::trajectory::{EffusionAdaptiveSolve, EffusionTrajectoryError};
use fs_material::thermal::transport::volumetric::timeline::FixedVolumeConductionLink;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::transport::effusion::channel::thermal::{
    FiniteChannelError, conduction::ChannelWallConduction,
};
use fs_thermal::enthalpy::{CellId, ThermalLink};
use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1, mixture::FrozenIdealGasMixtureModelV1};

fn adaptive(scale: f64) -> EffusionChannelAdaptiveSolve {
    EffusionChannelAdaptiveSolve {
        motion: EffusionAdaptiveSolve {
            initial_step_s: 0.1,
            minimum_step_s: 1e-6,
            maximum_step_s: 0.1,
            maximum_attempts: 4096,
            maximum_accepted_intervals: 2048,
            relative_tolerance: 1e-3 * scale,
            temperature_tolerance_k: 0.05 * scale,
            pressure_tolerance_pa: 5.0 * scale,
            species_mass_tolerance_kg: 1e-9 * scale,
            position_tolerance_m: 1e-5 * scale,
            velocity_tolerance_m_s: 1e-4 * scale,
        },
        wall_temperature_tolerance_k: 0.01 * scale,
        wall_energy_tolerance_j: 0.01 * scale,
        emitted_momentum_tolerance_n_s: 1e-7 * scale,
        emitted_energy_tolerance_j: 1e-4 * scale,
        window_mass_tolerance_kg: 1e-11,
        window_momentum_tolerance_n_s: 1e-9,
        window_energy_tolerance_j: 1e-4,
    }
}

fn check_window(before: &GasRocketState, duration: f64, result: &EffusionChannelTrajectoryCandidate<'_, '_>) {
    let step = result.step();
    let after = result.endpoint();
    assert_eq!(step.duration_s, duration);
    assert_eq!(step.accepted_steps.len(), 2 * step.accepted_intervals);
    assert_eq!(step.attempted_intervals, step.accepted_intervals + step.rejected_intervals);
    assert!(step.maximum_accepted_error <= 1.0);
    near(after.reservoir().mass_kg() - before.reservoir().mass_kg() + step.escaped_mass_kg, 0.0, 1e-11);
    vector_near(after.vehicle().velocity_m_s.scale(mass(after))
        .sub(before.vehicle().velocity_m_s.scale(mass(before)))
        .add(step.escaped_momentum_world_kg_m_s).sub(step.external_impulse_n_s), Vec3::ZERO, 2e-9);
    near(internal(after) - internal(before) + kinetic(after) - kinetic(before)
        + step.wall_storage_change_j + step.escaped_energy_world_j - step.external_work_j, 0.0, 1e-4);
    assert!(step.debit.energy_j <= adaptive(1.0).window_energy_tolerance_j);
    assert!(step.debit.mass_kg <= adaptive(1.0).window_mass_tolerance_kg);
    assert!(step.maximum_species_mass_residual_kg <= adaptive(1.0).window_mass_tolerance_kg);
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
    assert_eq!(after.vehicle().propellant_mass_kg, after.reservoir().mass_kg());
}

fn rk4<const N: usize>(mut state: [f64; N], duration: f64, derivative: impl Fn([f64; N]) -> [f64; N]) -> [f64; N] {
    let dt = duration / 8192.0;
    for _ in 0..8192 {
        let a = derivative(state);
        let b = derivative(core::array::from_fn(|i| state[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| state[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| state[i] + dt * c[i]));
        for i in 0..N { state[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0; }
    }
    state
}

// Independent continuous monatomic balances, using only resolved geometric
// coefficients. No production heat/effusion/time/mechanics step enters this
// reference. The wall has the authored nonlinear volumetric Cv = 1e6 + 1000 T.
fn continuous(path: &DiffuseChannelMomentumInterface, duration: f64) -> [f64; 7] {
    let tau = path.channel().geometry().transmission_probability();
    let heat = path.channel().thermal_accommodation().unwrap();
    let ht = heat.transmitted_heat_probability();
    let hr = heat.reflected_heat_probability();
    let id = path.channel().mouth().gas_catalog().thermodynamics().components()[0].species();
    let angular = path.axial_accommodation(id).unwrap();
    let retained = angular.retained_transmitted_axial_cosine();
    let reset = angular.reset_transmitted_axial_cosine();
    rk4([P0 * VOLUME / (R * T0), T0, 900.0, 0.0, 0.0, 0.0, 0.0], duration, |y| {
        let [n, t, w, _, v, _, _] = y;
        let incident = AREA * n / VOLUME * (R * t / (2.0 * core::f64::consts::PI * 0.028)).sqrt();
        let dn = -tau * incident;
        let mass_out = -0.028 * dn;
        let heat_return = 2.0 * R * incident * hr * (w - t);
        let heat_wall = -2.0 * R * incident * (ht + hr) * (w - t);
        let speed = |temperature: f64| (core::f64::consts::PI * R * temperature / (2.0 * 0.028)).sqrt();
        let impulse = 1.5 * incident * 0.028 * (retained * speed(t) + reset * speed(w));
        [dn, (-0.5 * R * incident * tau * t + heat_return) / (1.5 * R * n),
            heat_wall / (1e-6 * (1e6 + 1000.0 * w)), v,
            -impulse / (0.01 + 0.028 * n),
            mass_out * v + impulse,
            2.0 * R * incident * (tau * t + ht * (w - t)) + v * impulse + 0.5 * mass_out * v * v]
    })
}

#[test]
fn partial_accommodation_gas_wall_and_recoil_converge_to_independent_continuous_balances() {
    let before = initial(PURE, Vec3::ZERO);
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0).unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([(9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap())]);
    let path = channel(&before, 2.0, 0.35);
    let ports = [channel_outlet(&path, Vec3::new(1.0, 0.0, 0.0))];
    let duration = 0.4;
    let run = |dt| {
        let mut temporal = adaptive(1e8);
        temporal.motion.initial_step_s = dt;
        temporal.motion.maximum_step_s = dt;
        before.advance_effusion_channels_adaptive(environment(), &walls, &ports,
            duration, endpoint(), controls(), temporal, || false).unwrap()
    };
    let coarse = run(0.1);
    let fine = run(0.0125);
    let oracle = continuous(&path, duration);
    let discrepancy = |candidate: &EffusionChannelTrajectoryCandidate<'_, '_>| {
        let state = candidate.endpoint();
        [state.reservoir().mass_kg() / 0.028, state.reservoir().state().state().temperature,
            candidate.walls()[&9].temperature_k(), state.vehicle().position_m.x,
            state.vehicle().velocity_m_s.x, candidate.step().escaped_momentum_world_kg_m_s.x,
            candidate.step().escaped_energy_world_j]
            .into_iter().zip(oracle).zip([0.03, 400.0, 900.0, 0.01, 0.1, 0.002, 2.0])
            .map(|((actual, expected), scale)| (actual - expected).abs() / scale).sum::<f64>()
    };
    assert!(discrepancy(&fine) < 0.3 * discrepancy(&coarse));
    assert!(fine.step().accepted_intervals > coarse.step().accepted_intervals);
    assert!(fine.walls()[&9].temperature_k() < 900.0);
    assert!(fine.endpoint().vehicle().velocity_m_s.x < 0.0);
    check_window(&before, duration, &fine);
}

#[test]
fn opposing_shared_wall_channels_control_each_beam_and_replay_accepted_steps_exactly() {
    let before = initial(MIXTURE, Vec3::ZERO);
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0).unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([(9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap())]);
    let path = channel(&before, 2.0, 0.35);
    let ports = [channel_outlet(&path, Vec3::new(1.0, 0.0, 0.0)),
        EffusionChannelOutlet { id: 3, ..channel_outlet(&path, Vec3::new(-1.0, 0.0, 0.0)) }];
    let mut loose = adaptive(1e12);
    loose.motion.relative_tolerance = 0.0;
    let mut tight = loose;
    tight.emitted_momentum_tolerance_n_s = 1e-9;
    let run = |outlets: &[EffusionChannelOutlet<'_>], temporal| before
        .advance_effusion_channels_adaptive(environment(), &walls, outlets, 0.2,
            endpoint(), controls(), temporal, || false).unwrap();
    let coarse = run(&ports, loose);
    let fine = run(&ports, tight);
    let reordered = run(&[ports[1], ports[0]], tight);
    assert!(fine.step().accepted_intervals > coarse.step().accepted_intervals);
    assert_eq!(fine.step(), reordered.step());
    assert_eq!(fine.endpoint(), reordered.endpoint());
    vector_near(fine.endpoint().vehicle().velocity_m_s, Vec3::ZERO, 1e-14);
    assert_eq!(fine.step().escaped_momentum_world_kg_m_s, Vec3::ZERO);
    let mut replay = before.clone();
    let mut replay_walls = walls.clone();
    for report in &fine.step().accepted_steps {
        let candidate = replay.advance_effusion_channels(environment(), &replay_walls, &ports,
            report.duration_s, endpoint(), controls(), || false).unwrap();
        assert_eq!(candidate.step(), report);
        (replay, replay_walls, _) = candidate.into_parts();
    }
    assert_eq!(&replay, fine.endpoint());
    assert_eq!(replay_walls[&9].specific_energy_j_kg(), fine.walls()[&9].specific_energy_j_kg());
    let fractions = fine.endpoint().reservoir().model().thermodynamics().mole_composition().fractions();
    assert!(fractions[1] < 0.5); // Lighter B preferentially escapes both paths.
    check_window(&before, 0.2, &fine);
}

fn connection<'l, 'p>(law: &'l FixedVolumeThermalLaw<'p>) -> ChannelWallConduction<'l, 'p> {
    ChannelWallConduction::Cellwise(FixedVolumeConductionLink::try_aligned(
        &ThermalLink { id: 17, a: CellId(9), b: CellId(12), area_m2: 0.01, conductance_w_k: 0.0 },
        [law, law], [0.005; 2], 0.0).unwrap())
}

#[test]
fn conduction_only_walls_refine_and_retry_overshoot_without_adding_mass_thrust_or_external_heat() {
    let before = initial(PURE, Vec3::new(7.0, -1.0, 0.0));
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0).unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([
        (9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 400.0).unwrap()),
        (12, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1000.0).unwrap()),
    ]);
    let links = [connection(&law)];
    assert!(matches!(before.advance_effusion_channels_with_conduction(environment(), &walls, &[],
        &links, 2.0, endpoint(), controls(), || false),
        Err(EffusionChannelError::Channel(FiniteChannelError::WallOvershoot { .. }))));
    let run = |tolerance| {
        let mut temporal = adaptive(1e12);
        temporal.motion.relative_tolerance = 0.0;
        temporal.motion.initial_step_s = 2.0;
        temporal.motion.maximum_step_s = 2.0;
        temporal.wall_temperature_tolerance_k = tolerance;
        temporal.wall_energy_tolerance_j = tolerance;
        before.advance_effusion_channels_with_conduction_adaptive(environment(), &walls, &[],
            &links, 2.0, endpoint(), controls(), temporal, || false).unwrap()
    };
    let coarse = run(0.1);
    let fine = run(0.001);
    let oracle = rk4([400.0, 1000.0], 2.0, |[a, b]| {
        let power = b - a; // k=1, area=0.01, distance=0.01 => G=1 W/K.
        [power / (1.0 + 0.001 * a), -power / (1.0 + 0.001 * b)]
    });
    let error = |candidate: &EffusionChannelTrajectoryCandidate<'_, '_>| {
        (candidate.walls()[&9].temperature_k() - oracle[0]).abs()
            + (candidate.walls()[&12].temperature_k() - oracle[1]).abs()
    };
    assert!(error(&fine) < 0.3 * error(&coarse));
    assert!(fine.step().rejected_intervals > 0);
    assert!(fine.step().accepted_intervals > coarse.step().accepted_intervals);
    assert_eq!(fine.endpoint().reservoir(), before.reservoir());
    vector_near(fine.endpoint().vehicle().velocity_m_s, before.vehicle().velocity_m_s, 1e-13);
    assert_eq!(fine.step().escaped_mass_kg, 0.0);
    assert_eq!(fine.step().escaped_energy_world_j, 0.0);
    assert_eq!(fine.step().external_work_j, 0.0);
    near(fine.step().wall_storage_change_j, 0.0, 1e-8);
    assert!(fine.step().wall_storage_changes_j[&9] > 0.0);
    assert!(fine.step().wall_storage_changes_j[&12] < 0.0);
    check_window(&before, 2.0, &fine);
}

fn shifted_model(datum: f64) -> ThermochemicalGasModel {
    let model = source::load(PURE, "effusive recoil analytic fixture").unwrap();
    let species = model.thermodynamics().components().iter().map(|species| {
        Nasa9StandardStateModelV1::new(species.species().clone(), species.molar_mass(),
            species.convention().clone(), species.regions().iter().map(|region| {
                let mut card = region.card().clone();
                card.parameters.get_mut("a7").unwrap().value += datum;
                Nasa9RegionV1::from_card(card).unwrap()
            }).collect()).unwrap()
    }).collect();
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(species,
        model.thermodynamics().mole_composition().clone()).unwrap(),
        model.transport_materials().to_vec(), None).unwrap()
}

#[test]
fn accepted_export_energy_uses_fixed_caloric_references_and_common_inertial_transforms() {
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0).unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([(9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap())]);
    let mut temporal = adaptive(1e8);
    temporal.motion.initial_step_s = 0.025;
    temporal.motion.maximum_step_s = 0.025;
    let run = |datum, boost| {
        let gas = ClosedGasParcel::new(shifted_model(datum), source::point(T0, P0).unwrap(), VOLUME).unwrap();
        let before = GasRocketState::new(gas, 0.01, Vec3::ZERO, boost, FRAME).unwrap();
        let path = channel(&before, 2.0, 0.35);
        let result = before.advance_effusion_channels_adaptive(environment(), &walls,
            &[channel_outlet(&path, Vec3::new(1.0, 0.0, 0.0))], DT,
            endpoint(), controls(), temporal, || false).unwrap();
        (before, result)
    };
    let (before, a) = run(0.0, Vec3::ZERO);
    let (_, shifted) = run(1e12, Vec3::ZERO);
    let boost = Vec3::new(11.0, -3.0, 7.0);
    let (moving, b) = run(0.0, boost);
    near(shifted.step().escaped_reference_energy_world_j, a.step().escaped_reference_energy_world_j, 1e-8);
    near(shifted.endpoint().reservoir().state().state().temperature,
        a.endpoint().reservoir().state().state().temperature, 1e-7);
    assert_eq!(a.step().accepted_intervals, b.step().accepted_intervals);
    assert_eq!(a.endpoint().reservoir(), b.endpoint().reservoir());
    vector_near(b.step().escaped_momentum_world_kg_m_s,
        a.step().escaped_momentum_world_kg_m_s.add(boost.scale(a.step().escaped_mass_kg)), 1e-12);
    let correction = boost.dot(a.step().escaped_momentum_world_kg_m_s)
        + 0.5 * boost.norm_squared() * a.step().escaped_mass_kg;
    near(b.step().escaped_reference_energy_world_j,
        a.step().escaped_reference_energy_world_j + correction, 1e-10);
    vector_near(b.endpoint().vehicle().position_m.sub(a.endpoint().vehicle().position_m), boost.scale(DT), 1e-12);
    check_window(&before, DT, &a);
    check_window(&moving, DT, &b);
}

#[test]
fn fixed_window_budgets_and_typed_physical_temporal_refusals_leave_all_inputs_unchanged() {
    let before = initial(PURE, Vec3::ZERO);
    let saved = before.clone();
    let material = wall_material();
    let law = FixedVolumeThermalLaw::try_new(&material, 1000.0, 200.0).unwrap().with_sensible_entropy();
    let walls = BTreeMap::from([(9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap())]);
    let original_wall = walls[&9].specific_energy_j_kg();
    let path = channel(&before, 2.0, 0.35);
    let ports = [channel_outlet(&path, Vec3::new(1.0, 0.0, 0.0))];
    let mut gross = controls();
    gross.transport.exchange.maximum_species_fraction = 2e-4;
    let retried = before.advance_effusion_channels_adaptive(environment(), &walls, &ports,
        DT, endpoint(), gross, adaptive(1e8), || false).unwrap();
    assert!(retried.step().rejected_intervals > 0);
    let mut minimal = adaptive(1e-8);
    minimal.motion.minimum_step_s = minimal.motion.initial_step_s;
    assert!(matches!(before.advance_effusion_channels_adaptive(environment(), &walls, &ports,
        DT, endpoint(), controls(), minimal, || false),
        Err(EffusionChannelTrajectoryError::Control(EffusionTrajectoryError::MinimumStep { .. }))));
    let mut bounded = adaptive(1e-8);
    bounded.motion.maximum_attempts = 1;
    assert!(matches!(before.advance_effusion_channels_adaptive(environment(), &walls, &ports,
        DT, endpoint(), controls(), bounded, || false),
        Err(EffusionChannelTrajectoryError::Control(EffusionTrajectoryError::WorkBudget { .. }))));
    let unsupported = BTreeMap::from([(9, FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap())]);
    assert!(matches!(before.advance_effusion_channels_adaptive(environment(), &unsupported, &ports,
        DT, endpoint(), controls(), adaptive(1.0), || false),
        Err(EffusionChannelTrajectoryError::Step { duration_s, .. }) if duration_s == DT));
    let mut checks = 0;
    assert!(before.advance_effusion_channels_adaptive(environment(), &walls, &ports,
        DT, endpoint(), controls(), adaptive(1.0), || { checks += 1; checks >= 80 }).is_err());
    let mut perturbed = initial(PURE, Vec3::new(1.0, 0.0, 0.0));
    perturbed.vehicle.propellant_mass_kg = perturbed.vehicle.propellant_mass_kg.next_up();
    let mut exact = adaptive(1e8);
    exact.window_mass_tolerance_kg = 0.0;
    assert!(matches!(perturbed.advance_effusion_channels_adaptive(environment(), &walls, &[],
        DT, endpoint(), controls(), exact, || false),
        Err(EffusionChannelTrajectoryError::WindowBalance { debit, .. }) if debit.mass_kg > 0.0));
    let closed = before.advance_effusion_channels_adaptive(environment(), &unsupported,
        &[EffusionChannelOutlet { area_m2: 0.0, ..ports[0] }], DT,
        endpoint(), controls(), adaptive(1.0), || false).unwrap();
    assert_eq!(closed.endpoint(), &before);
    assert_eq!(closed.walls()[&9].specific_energy_j_kg(), unsupported[&9].specific_energy_j_kg());
    assert_eq!(before, saved);
    assert_eq!(walls[&9].specific_energy_j_kg(), original_wall);
}
