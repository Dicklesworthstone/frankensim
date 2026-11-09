//! Shared local optical/kinetic closure for one or many simultaneous paths.

use super::*;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Channel {
    pub incident_rate: f64,
    pub energy: f64,
    pub depth: [f64; 2],
    pub passive_depth: f64,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ChannelStep {
    pub incident: f64,
    pub transmitted: f64,
    pub absorbed: f64,
    /// Original target counts before spectral energy partition.
    pub absorbed_by_state: [f64; 2],
    pub passive: f64,
    pub transmitted_rate: f64,
    pub residual: f64,
}

pub(crate) struct SharedTrial {
    pub kinetics: SaturableSingletTripletStep,
    pub channels: Vec<ChannelStep>,
    pub mean_ground: f64,
    pub passive_heat_j: f64,
}

/// A Beer field at an independently supplied common ground inventory. Kinetics
/// are deliberately absent: tracing a beam must never advance its own copy of
/// a shared emitter inventory.
pub(crate) fn transmitted_rate(
    channel: Channel,
    mean_ground: f64,
) -> Result<f64, SingletTripletError> {
    Ok(attenuate(channel.incident_rate, depth(channel, mean_ground)?))
}

fn depth(channel: Channel, mean_ground: f64) -> Result<f64, SingletTripletError> {
    finite(
        mean_ground * channel.depth.iter().sum::<f64>() + channel.passive_depth,
        "pump optical depth",
    )
}

pub(crate) fn consistent(trial: &SharedTrial, guess: f64, tolerance: f64) -> bool {
    let close = (trial.mean_ground - guess).abs()
        <= 0.0625 * tolerance * trial.mean_ground.max(guess);
    close && photons_close(trial, 0.25 * tolerance)
}

pub(crate) fn photons_close(trial: &SharedTrial, tolerance: f64) -> bool {
    trial.channels.iter().all(|p| {
        admit_photons(
            [p.incident],
            [p.transmitted],
            [p.absorbed],
            [p.passive],
            [p.residual],
            tolerance,
        )
        .is_ok()
    })
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn solve(
    model: SingletTripletDecay,
    populations: [f64; 3],
    channels: &[Channel],
    duration: f64,
    limits: PumpColumnLimits,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SharedTrial, SingletTripletError> {
    let mut total = Sum::default();
    for value in populations {
        total.add(value)?;
    }
    let total = total.value()?;
    let active = duration > 0.0
        && total > 0.0
        && channels
            .iter()
            .any(|c| c.incident_rate > 0.0 && c.depth.iter().any(|&x| x > 0.0));
    if !active {
        return evaluate(
            model, populations, channels, duration, populations[0], limits, evaluations, cancelled,
        );
    }
    let mut low = 0.0;
    let mut high = total;
    let mut f_low = 0.0;
    let mut f_high = 0.0;
    for iteration in 0..limits.maximum_root_iterations {
        poll(cancelled)?;
        let guess = if iteration == 0 {
            low
        } else if iteration == 1 {
            high
        } else if iteration % 4 == 0 {
            low + 0.5 * (high - low)
        } else {
            let scale = f64::max(f_low, -f_high);
            let weight = (f_low / scale) / (f_low / scale - f_high / scale);
            let candidate = low + weight * (high - low);
            if candidate > low && candidate < high {
                candidate
            } else {
                low + 0.5 * (high - low)
            }
        };
        let trial = evaluate(
            model, populations, channels, duration, guess, limits, evaluations, cancelled,
        )?;
        let residual = finite(trial.mean_ground - guess, "mean ground residual")?;
        if consistent(&trial, guess, limits.kinetic.relative_tolerance) {
            return Ok(trial);
        }
        if residual > 0.0 {
            low = guess;
            f_low = residual;
        } else {
            high = guess;
            f_high = residual;
        }
    }
    cap(
        "pump column ground-residence iterations",
        limits.maximum_root_iterations.saturating_add(1),
        limits.maximum_root_iterations,
    )?;
    unreachable!("exceeded root iteration cap")
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn evaluate(
    model: SingletTripletDecay,
    populations: [f64; 3],
    channels: &[Channel],
    duration: f64,
    mean_ground: f64,
    limits: PumpColumnLimits,
    evaluations: &mut usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SharedTrial, SingletTripletError> {
    *evaluations = evaluations
        .checked_add(1)
        .ok_or(SingletTripletError::Unrepresentable("pump cell evaluation counter"))?;
    cap("pump column cell evaluations", *evaluations, limits.maximum_cell_evaluations)?;
    poll(cancelled)?;
    let mut hazards = Vec::with_capacity(channels.len());
    let mut packets = Vec::with_capacity(channels.len());
    let mut rates = [Sum::default(); 2];
    for &channel in channels {
        poll(cancelled)?;
        let depth = depth(channel, mean_ground)?;
        let mean_transmission = if depth == 0.0 {
            1.0
        } else {
            decayed_fraction(depth) / depth
        };
        let incident = positive_product(channel.incident_rate, duration, "incident photons")?;
        let transmitted_rate = attenuate(channel.incident_rate, depth);
        let transmitted = attenuate(incident, depth);
        if transmitted_rate == 0.0 && transmitted > 0.0 {
            return Err(SingletTripletError::Unrepresentable("transmitted photon rate"));
        }
        let passive = product3(incident, channel.passive_depth, mean_transmission)?;
        let mut hazard = [0.0; 2];
        for target in 0..2 {
            hazard[target] = product3(
                channel.incident_rate, channel.depth[target], mean_transmission,
            )?;
            rates[target].add(hazard[target])?;
        }
        hazards.push(hazard);
        packets.push(ChannelStep {
            incident, transmitted, absorbed: 0.0, absorbed_by_state: [0.0; 2], passive,
            transmitted_rate, residual: 0.0,
        });
    }
    let rates = values(rates)?;
    // Original channel quanta are retained until after actual event splitting.
    // Averaging energies first erases weak but positive pump-relaxation heat.
    let plan = SaturableSingletTripletPlan::prepare_with_population_integrals(
        model, rates, model.excitation_energies_j(), duration,
        limits.kinetic, &mut *cancelled,
    )?;
    let mut kinetics = plan.advance_with_cancel(populations, &mut *cancelled)?;
    let residence = plan.integrate_populations_with_cancel(populations, &mut *cancelled)?;
    let mean_ground = if duration == 0.0 {
        populations[0]
    } else {
        finite(residence[0] / duration, "mean ground population")?
    };
    let mut pump_energy = [Sum::default(); 2];
    let mut pump_heat = [Sum::default(); 2];
    let mut passive_heat = Sum::default();
    for ((channel, hazard), packet) in channels.iter().zip(hazards).zip(&mut packets) {
        let mut events = Sum::default();
        for target in 0..2 {
            if rates[target] > 0.0 {
                let weight = hazard[target] / rates[target];
                if hazard[target] > 0.0 && weight == 0.0 {
                    return Err(SingletTripletError::Unrepresentable("spectral pump branching"));
                }
                let count = positive_product(
                    kinetics.pump_events[target], weight, "spectral pump events",
                )?;
                packet.absorbed_by_state[target] = count;
                events.add(count)?;
                pump_energy[target].add(positive_product(
                    count, channel.energy, "spectral pump energy",
                )?)?;
                let excess = if hazard[target] > 0.0 {
                    channel.energy - model.excitation_energies_j()[target]
                } else {
                    0.0
                };
                pump_heat[target].add(positive_product(
                    count, excess, "spectral pump excess heat",
                )?)?;
            }
        }
        packet.absorbed = events.value()?;
        packet.residual = photon_balance(
            [packet.incident], [packet.transmitted], [packet.absorbed], [packet.passive],
        )?[0];
        passive_heat.add(positive_product(
            packet.passive, channel.energy, "passive pump energy",
        )?)?;
    }
    kinetics.absorbed_pump_energy_j = values(pump_energy)?;
    kinetics.pump_relaxation_heat_j = values(pump_heat)?;
    kinetics.residual_j = kinetic_energy_residual(
        model, populations, &kinetics, limits.kinetic.relative_tolerance,
    )?;
    Ok(SharedTrial {
        kinetics, channels: packets, mean_ground, passive_heat_j: passive_heat.value()?,
    })
}
