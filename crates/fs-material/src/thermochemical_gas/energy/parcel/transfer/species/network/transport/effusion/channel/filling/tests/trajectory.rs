//! Adaptive vacuum startup, real receiver backflow and original finite storage.
//! Synthetic source laws are numerical fixtures, not empirical validation.

use super::*;
use super::super::trajectory::{
    EffusionAdvanceControls, FiniteChannelFillAdvanceCandidate,
    FiniteChannelFillAdvanceError, advance_candidate as advance,
};
use super::super::super::thermal::trajectory::FiniteChannelAdvanceError;

fn accuracy(scale: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: 0.05,
        minimum_step_s: 1e-6,
        maximum_step_s: 0.05,
        maximum_attempts: 16_384,
        maximum_accepted_intervals: 4096,
        relative_tolerance: 1e-3 * scale,
        temperature_tolerance_k: 0.05 * scale,
        pressure_tolerance_pa: 2.0 * scale,
        species_mass_tolerance_kg: 1e-10 * scale,
        stream_energy_tolerance_j: 5e-5 * scale,
    }
}
fn vacuum(volume: f64) -> EvacuatedGasVolume {
    EvacuatedGasVolume::new(volume, QueryPoint::new()).unwrap()
}
fn temperatures(candidate: &FiniteChannelFillAdvanceCandidate<'_, '_>) -> [f64; 3] {
    [
        candidate.donor().state().state().temperature,
        candidate.receiver().gas().unwrap().state().state().temperature,
        candidate.wall().temperature_k(),
    ]
}
fn closed_balance(candidate: &FiniteChannelFillAdvanceCandidate<'_, '_>) {
    let report = candidate.step();
    let b = report.balance;
    assert!(b.total_mass_residual_kg.abs() <= controls().exchange.mass_tolerance_kg);
    assert!(b.maximum_species_mass_residual_kg <= controls().exchange.mass_tolerance_kg);
    assert!(b.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    assert!(b.total_entropy_change_j_k >= -controls().entropy_tolerance_j_k);
    assert!(report.maximum_accepted_error <= 1.0);
    if report.first_reciprocal.is_some() {
        assert_eq!(2 + report.continuation.len(), 2 * report.accepted_intervals);
        assert!(report.attempted_steps >= 3 * report.accepted_intervals);
    }
}

// Independent continuous monatomic mass/internal-energy balances. The receiver
// starts with exactly n=U=0; its nonexistent T is never evaluated. Incoming and
// outgoing beams use the same frozen geometric probabilities as the test input,
// but this oracle uses no library gas step, wall inversion or trajectory driver.
fn continuous(path: &DiffuseMolecularChannel, duration: f64) -> [f64; 3] {
    let tau = path.geometry().transmission_probability();
    let hit_transmission = path.geometry().wall_transmitted_probability();
    let cv = 1.5 * R;
    let derivative = |y: [f64; 5]| {
        let [nl, ul, nr, ur, tw] = y;
        let tl = ul / (cv * nl);
        let left = AREA * nl / 0.001
            * (R * tl / (2.0 * core::f64::consts::PI * M)).sqrt();
        let right = if nr > 0.0 {
            let tr = ur / (cv * nr);
            let rate = AREA * nr / 0.001
                * (R * tr / (2.0 * core::f64::consts::PI * M)).sqrt();
            [rate, rate * tr]
        } else {
            [0.0; 2]
        };
        let ql = 2.0 * R
            * ((1.0 - tau) * left * (tw - tl)
                + hit_transmission * (tw * right[0] - right[1]));
        let qr = 2.0 * R
            * ((1.0 - tau) * (tw * right[0] - right[1])
                + hit_transmission * left * (tw - tl));
        let amount = tau * (left - right[0]);
        let energy = 2.0 * R * tau * (left * tl - right[1]);
        [
            -amount,
            -energy + ql,
            amount,
            energy + qr,
            -(ql + qr) / (1e-6 * (1e6 + 700.0 * tw)),
        ]
    };
    let n0 = 100_000.0 * 0.001 / (R * 400.0);
    let mut y = [n0, n0 * cv * 400.0, 0.0, 0.0, 900.0];
    let dt = duration / 8192.0;
    for _ in 0..8192 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..5 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    [y[1] / (cv * y[0]), y[3] / (cv * y[2]), y[4]]
}

#[test]
fn tighter_thermal_filling_converges_to_continuous_backflow_and_finite_wall_feedback() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap();
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = vacuum(0.001);
    let path = channel(0.0);
    let duration = 0.4;
    let run = |scale| advance(
        &donor, &receiver, &wall, link(&path), duration,
        energy(), energy(), controls(), accuracy(scale), || false,
    ).unwrap();
    let loose = run(1.0);
    let tight = run(0.04);
    let reference = continuous(&path, duration);
    let discrepancy = |candidate: &FiniteChannelFillAdvanceCandidate<'_, '_>| {
        temperatures(candidate).into_iter().zip(reference)
            .map(|(a, b)| (a - b).abs()).sum::<f64>()
    };
    assert!(discrepancy(&tight) < 0.4 * discrepancy(&loose));
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().rejected_intervals > 0);
    assert_eq!(tight.step().duration_s, duration);
    assert!(tight.wall().temperature_k() < wall.temperature_k());
    assert!(tight.donor().state().state().temperature > 400.0);
    assert!(tight.step().first_reciprocal.as_ref().unwrap().channels[0]
        .transmitted.species[0].right_to_left.mass_rate_kg_s > 0.0);
    assert!(tight.step().balance.total_entropy_change_j_k > 0.0);
    closed_balance(&tight);
    assert_eq!(wall.temperature_k(), 900.0);
    assert_eq!(donor.state().state().temperature, 400.0);
}

fn mixture(datum: f64) -> (ClosedGasParcel, DiffuseMolecularChannel) {
    let rows = [("A", 0.028, 0.0), ("B", 0.004, 1.0)];
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter().map(|&(id, m, _)| {
                species_with_energy_reference(id, m, 2.5, 0.0, datum)
            }).collect(),
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
        ).unwrap(),
        rows.iter().map(|&(id, _, _)| {
            transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1e-8, 1e7])
        }).collect(),
        None,
    ).unwrap();
    let donor = ClosedGasParcel::new(model.clone(), point(400.0, 5000.0), 0.001).unwrap();
    let mouth = EffusiveInterface::try_new(
        model,
        rows.iter().map(|&(id, _, _)| EffusiveSpeciesTransmission::try_new(
            SpeciesId::new(id).unwrap(), 1.0, "synthetic reciprocal mouth",
        ).unwrap()).collect(),
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1e-8, 1e7),
        },
        EffusiveRegime::try_new(1e-8, 1e-11, 1e-4, "synthetic collisionless support").unwrap(),
    ).unwrap();
    let path = DiffuseMolecularChannel::try_new(
        mouth,
        CylindricalMolecularGeometry::try_new(
            1e-8, 2e-8, 1e-4, "synthetic two-species thermal passage",
        ).unwrap(),
    ).unwrap().with_species_thermal_accommodation(
        rows.iter().map(|&(id, _, coefficient)| SpeciesThermalAccommodation::try_new(
            SpeciesId::new(id).unwrap(), coefficient, "synthetic species-surface qualification",
        ).unwrap()).collect(),
        ClausingControls::default(), || false,
    ).unwrap();
    (donor, path)
}

#[test]
fn species_accommodation_histories_replay_exactly_without_formation_datum_error() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap();
    let receiver = vacuum(0.001);
    let duration = 0.08;
    let mut original: Option<([f64; 3], usize)> = None;
    for datum in [0.0, 1e20] {
        let (donor, path) = mixture(datum);
        let result = advance(
            &donor, &receiver, &wall, link(&path), duration,
            energy(), energy(), controls(), accuracy(1.0), || false,
        ).unwrap();
        if let Some((values, intervals)) = original {
            for (actual, expected) in temperatures(&result).into_iter().zip(values) {
                near(actual, expected, 1e-7);
            }
            assert_eq!(result.step().accepted_intervals, intervals);
        } else {
            original = Some((temperatures(&result), result.step().accepted_intervals));
        }
        let report = result.step();
        let first = step_candidate(
            &donor, &EffusionReceiver::Evacuated(receiver.clone()), &wall,
            link(&path), report.first_fill.duration_s,
            energy(), energy(), controls(), || false,
        ).unwrap();
        assert_eq!(first.step(), &report.first_fill);
        let (mut left, mut right, mut current_wall, _) = first.into_parts();
        let mut elapsed = Dd::from_f64(report.first_fill.duration_s);
        for accepted in report.first_reciprocal.iter().chain(&report.continuation) {
            let next = step_candidate(
                &left, &right, &current_wall, link(&path), accepted.duration_s,
                energy(), energy(), controls(), || false,
            ).unwrap();
            assert!(matches!(&next.step().transfer,
                FiniteChannelFillTransfer::Exchange(actual) if actual == accepted));
            (left, right, current_wall, _) = next.into_parts();
            elapsed = elapsed + Dd::from_f64(accepted.duration_s);
        }
        assert_eq!(&left, result.donor());
        assert_eq!(&right, result.receiver());
        assert_eq!(current_wall.specific_energy_j_kg(), result.wall().specific_energy_j_kg());
        near(elapsed.to_f64(), duration, 1e-12);
        let gas = result.receiver().gas().unwrap();
        let fractions = gas.model().thermodynamics().mole_composition().fractions();
        assert!(fractions[1] > 0.75, "lighter B preferentially fills the receiver");
        closed_balance(&result);
    }
}

#[test]
fn closed_and_elastic_limits_and_pressure_limited_startup_preserve_exact_vacuum() {
    let source = material(1e6, 0.0, [200.0, 2000.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap();
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = vacuum(0.001);
    let path = channel(0.0);
    let closed = advance(
        &donor, &receiver, &wall,
        FiniteChannelFillLink { area_m2: 0.0, ..link(&path) }, 1.0,
        energy(), energy(), controls(), accuracy(1.0), || false,
    ).unwrap();
    assert_eq!(closed.donor(), &donor);
    assert_eq!(closed.receiver(), &EffusionReceiver::Evacuated(receiver.clone()));
    assert_eq!(closed.wall().specific_energy_j_kg(), wall.specific_energy_j_kg());
    assert_eq!(closed.step().attempted_steps, 1);
    assert_eq!(closed.step().first_fill.duration_s, 1.0);
    assert!(closed.step().first_reciprocal.is_none());
    assert!(closed.step().continuation.is_empty());
    let elastic = path.with_thermal_accommodation(
        0.0, ClausingControls::default(), "synthetic elastic diffuse wall", || false,
    ).unwrap();
    let filled = advance(
        &donor, &receiver, &wall, link(&elastic), 0.025,
        energy(), energy(), controls(), accuracy(1.0), || false,
    ).unwrap();
    assert!(filled.receiver().gas().is_some());
    assert_eq!(filled.wall().specific_energy_j_kg(), wall.specific_energy_j_kg());
    assert_eq!(filled.step().balance.wall_energy_change_j, 0.0);
    closed_balance(&filled);

    let small = vacuum(1e-9);
    let mut limited = accuracy(1.0);
    limited.minimum_step_s = 0.001;
    assert!(matches!(advance(
        &donor, &small, &wall, link(&elastic), 0.01,
        energy(), energy(), controls(), limited, || false,
    ), Err(FiniteChannelFillAdvanceError::InitialPressureLimit {
        maximum_step_s, minimum_step_s,
    }) if maximum_step_s < minimum_step_s));
    limited.minimum_step_s = 1e-8;
    let admitted = advance(
        &donor, &small, &wall, link(&elastic), 1e-4,
        energy(), energy(), controls(), limited, || false,
    ).unwrap();
    let transmitted_rate = AREA * 100_000.0 * elastic.geometry().transmission_probability()
        / (2.0 * core::f64::consts::PI * M * R * 400.0).sqrt();
    let pressure_ceiling = 100_000.0 * small.volume_m3()
        / (transmitted_rate * R * energy().maximum_temperature_k);
    assert!(2.0 * admitted.step().first_fill.duration_s <= pressure_ceiling);
    assert_eq!(admitted.receiver().gas().unwrap().volume_m3(), small.volume_m3());
    closed_balance(&admitted);
}

#[test]
fn source_time_work_and_cancellation_refusals_publish_none_of_the_three_endpoints() {
    let source = material(1e6, 700.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap().with_sensible_entropy();
    let wall = FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap();
    let donor = gas(400.0, 100_000.0, 0.0);
    let receiver = vacuum(0.001);
    let original = donor.clone();
    let original_empty = receiver.clone();
    let path = channel(0.0);
    let mut strict = accuracy(1e-8);
    strict.minimum_step_s = strict.initial_step_s;
    assert!(matches!(advance(
        &donor, &receiver, &wall, link(&path), 0.4,
        energy(), energy(), controls(), strict, || false,
    ), Err(FiniteChannelFillAdvanceError::Advance(
        FiniteChannelAdvanceError::MinimumStep { elapsed_s: 0.0, .. }
    ))));
    let work = EffusionAdvanceControls { maximum_attempts: 3, ..accuracy(1.0) };
    assert!(matches!(advance(
        &donor, &receiver, &wall, link(&path), 0.4,
        energy(), energy(), controls(), work, || false,
    ), Err(FiniteChannelFillAdvanceError::Advance(
        FiniteChannelAdvanceError::WorkBudget { attempted_steps: 3, .. }
    ))));
    assert!(matches!(advance(
        &donor, &receiver, &wall, link(&path), 0.1,
        energy(), GasEnergySolve { maximum_temperature_k: 300.0, ..energy() },
        controls(), accuracy(1.0), || false,
    ), Err(FiniteChannelFillAdvanceError::Fill {
        duration_s, error: FiniteChannelFillError::Receiver(_),
    }) if duration_s == accuracy(1.0).initial_step_s));

    let mut polls = 0;
    assert!(matches!(advance(
        &donor, &receiver, &wall, link(&path), 0.02,
        energy(), energy(), controls(), accuracy(1.0),
        || { polls += 1; polls == 3 },
    ), Err(FiniteChannelFillAdvanceError::Cancelled)));
    let mut total = 0;
    advance(
        &donor, &receiver, &wall, link(&path), 0.02,
        energy(), energy(), controls(), accuracy(1.0),
        || { total += 1; false },
    ).unwrap();
    let mut seen = 0;
    assert!(matches!(advance(
        &donor, &receiver, &wall, link(&path), 0.02,
        energy(), energy(), controls(), accuracy(1.0),
        || { seen += 1; seen == total },
    ), Err(FiniteChannelFillAdvanceError::Cancelled)));
    assert_eq!(donor, original);
    assert_eq!(receiver, original_empty);
    assert_eq!(wall.temperature_k(), 900.0);
}
