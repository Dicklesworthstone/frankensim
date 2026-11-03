//! Conservative finite thermal ports on the original wall inventories.
//!
//! Each source owner supplies the two current stage states and owns ordered
//! composition and publication. Frozen and evolving-mixture networks use the
//! same original scalar equation; no capacity or latent chart is reconstructed.

use super::*;

#[cfg(feature = "material")]
mod conduction;

/// One stationary, passive thermal port between two original wall inventories.
/// Order is significant at finite timestep: links execute in the supplied
/// order. Parallel interfaces and cycles are admitted without merging them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedWallLink {
    /// Distinct indices in the bound shared-wall array, oriented a then b.
    pub walls: [usize; 2],
    /// Complete authored conductance [W/K], finite and nonnegative. Geometry
    /// or temperature-dependent conductivity is not inferred from this value.
    pub conductance_w_per_k: f64,
}

/// Original states and physical accounts of one accepted link stage.
/// Later links may change these walls again; these temperatures are deliberately
/// retained rather than relabelled as a simultaneous whole-network endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedWallLinkReport {
    /// Actual source-stage duration, used once by this link stage.
    pub elapsed: Time,
    /// Original states entering this particular composition stage.
    pub initial_state: [WallThermal; 2],
    /// Accepted original states leaving this particular composition stage.
    pub endpoint_state: [WallThermal; 2],
    /// Initial temperatures in the link's a,b order [K].
    pub initial_temperature_k: [f64; 2],
    /// Endpoint temperatures in the same order [K].
    pub endpoint_temperature_k: [f64; 2],
    /// Actual represented energy received by wall b [J]; negative reverses it.
    pub heat_a_to_b_j: f64,
    /// Actual original storage changes of wall a and b [J].
    pub wall_energy_change_j: [f64; 2],
    /// Sum of the two actual storage changes [J].
    pub energy_residual_j: f64,
    /// Q - dt G (Ta_new - Tb_new) at this stage [J].
    pub flux_residual_j: f64,
    /// Safeguarded scalar iterations, excluding bracket queries.
    pub iterations: usize,
}

pub(crate) struct WallLinks {
    pub(crate) ports: Vec<SharedWallLink>,
    pub(crate) candidate: Vec<SharedWallLinkReport>,
    pub(crate) accepted: Vec<SharedWallLinkReport>,
    #[cfg(feature = "material")]
    pub(crate) conduction: Option<Vec<Option<conduction::LinkConduction>>>,
}

// Heat leaving this original wall before either constitutive endpoint. An
// authored constant-capacity wall uses the same passive initial temperature
// hull as the existing gas/wall exchange. Source h/u retains its latent range.
fn heat_limits(
    group: &Group<'_>,
    state: WallThermal,
    other_temperature: f64,
    work: &Work,
) -> Result<[f64; 2], WallExchangeError> {
    work.charge()?;
    if let Some(source) = &group.accepted.storage.heat_storage {
        return source.heat_limits(
            WallContact {
                drag: WallDrag::None,
                thermal: state,
            },
            other_temperature,
        );
    }
    let WallThermal::Lumped {
        temperature_k,
        heat_capacity_j_per_k,
        ..
    } = state
    else {
        return Err(WallExchangeError::Invalid("wall link finite storage"));
    };
    let heat = heat_capacity_j_per_k * (temperature_k - other_temperature);
    if !heat.is_finite() {
        return Err(WallExchangeError::Invalid("wall link storage bound"));
    }
    Ok([heat.min(0.0), heat.max(0.0)])
}

pub(crate) fn exchange(
    [a, b]: [&Group<'_>; 2],
    links: &WallLinks,
    index: usize,
    duration: Time,
    control: WallControl,
    work: &Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SharedWallLinkReport, WallExchangeError> {
    poll(cancelled)?;
    let link = links.ports[index];
    let states = [a.candidate, b.candidate];
    let initial = [
        a.after_state(states[0], 0.0, work)?,
        b.after_state(states[1], 0.0, work)?,
    ];
    let temperatures = initial.each_ref().map(|p| p.temperature);
    let conductance = link.conductance_w_per_k;
    #[cfg(feature = "material")]
    let prepared = links
        .conduction
        .as_ref()
        .and_then(|bindings| bindings[index].as_ref())
        .map(|binding| binding.prepare([a, b], states, conductance, work, cancelled))
        .transpose()?;
    #[cfg(feature = "material")]
    let conductance = prepared
        .as_ref()
        .map_or(conductance, |binding| binding.initial);
    let integrated = duration.value() * conductance;
    let initial_flux = integrated * (temperatures[0] - temperatures[1]);
    if !integrated.is_finite() || !initial_flux.is_finite() {
        return Err(WallExchangeError::Invalid(
            "wall link integrated conductance",
        ));
    }
    let mut receipt = SharedWallLinkReport {
        elapsed: duration,
        initial_state: states,
        endpoint_state: states,
        initial_temperature_k: temperatures,
        endpoint_temperature_k: temperatures,
        heat_a_to_b_j: 0.0,
        wall_energy_change_j: [0.0; 2],
        energy_residual_j: 0.0,
        flux_residual_j: -initial_flux,
        iterations: 0,
    };
    if initial_flux.abs() <= control.tolerance(initial_flux.abs()) {
        poll(cancelled)?;
        return Ok(receipt);
    }
    let la = heat_limits(a, states[0], temperatures[1], work)?;
    let lb = heat_limits(b, states[1], temperatures[0], work)?;
    let flux_bound = initial_flux;
    #[cfg(feature = "material")]
    let flux_bound = if prepared.as_ref().is_some_and(|binding| binding.implicit) {
        // Increasing k can make the actual transfer exceed the frozen initial
        // flux. Only original storage support and passive temperature bounds
        // constrain this nonlinear path; the fixed-G path retains its bound.
        if initial_flux > 0.0 {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        }
    } else {
        flux_bound
    };
    let (mut low, mut high) = if initial_flux > 0.0 {
        (0.0, flux_bound.min(la[1]).min(-lb[0]))
    } else {
        (flux_bound.max(la[0]).max(-lb[1]), 0.0)
    };
    if low < 0.0 {
        low = low.next_up();
    }
    if high > 0.0 {
        high = high.next_down();
    }
    #[cfg(feature = "material")]
    if prepared.as_ref().is_some_and(|binding| binding.implicit) {
        if initial_flux > 0.0 {
            high = (-a.bound(
                states[0],
                -high,
                temperatures[1],
                true,
                control,
                work,
                cancelled,
            )?)
            .min(b.bound(
                states[1],
                high,
                temperatures[0],
                false,
                control,
                work,
                cancelled,
            )?);
        } else {
            low = (-a.bound(
                states[0],
                -low,
                temperatures[1],
                false,
                control,
                work,
                cancelled,
            )?)
            .max(b.bound(
                states[1],
                low,
                temperatures[0],
                true,
                control,
                work,
                cancelled,
            )?);
        }
    }
    if !low.is_finite() || !high.is_finite() || low >= high {
        return Err(WallExchangeError::Invalid(
            "wall link supported heat interval",
        ));
    }
    let (heat, iterations) = solve_scalar(
        low,
        high,
        0.0,
        control,
        |heat, checkpoint| {
            if checkpoint() {
                return Err(WallExchangeError::Cancelled);
            }
            let wa = a.after_state(states[0], -heat, work)?;
            let wb = b.after_state(states[1], heat, work)?;
            #[cfg(feature = "material")]
            let integrated = match &prepared {
                Some(binding) => {
                    duration.value()
                        * binding.at(
                            [a, b],
                            [wa.thermal, wb.thermal],
                            link.conductance_w_per_k,
                            work,
                            checkpoint,
                        )?
                }
                None => integrated,
            };
            let flux = integrated * (wa.temperature - wb.temperature);
            let residual = heat - flux;
            let slope = wa
                .slope
                .zip(wb.slope)
                .map(|(a, b)| 1.0 + integrated * (a + b));
            #[cfg(feature = "material")]
            let slope = if prepared.as_ref().is_some_and(|binding| binding.implicit) {
                None
            } else {
                slope
            };
            if !residual.is_finite() || slope.is_some_and(|s| !positive(s)) {
                return Err(WallExchangeError::Invalid("wall link heat equation"));
            }
            Ok((residual, slope, heat.abs().max(flux.abs())))
        },
        cancelled,
    )?;
    let endpoint = [
        a.after_state(states[0], -heat, work)?,
        b.after_state(states[1], heat, work)?,
    ];
    let actual_heat = endpoint[1].change;
    let energy_residual = endpoint[0].change + actual_heat;
    #[cfg(feature = "material")]
    let integrated = match &prepared {
        Some(binding) => {
            let thermal = endpoint.each_ref().map(|point| point.thermal);
            binding.admit_lagged_endpoint([a, b], thermal, work, cancelled)?;
            duration.value()
                * binding.at([a, b], thermal, link.conductance_w_per_k, work, cancelled)?
        }
        None => integrated,
    };
    let flux = integrated * (endpoint[0].temperature - endpoint[1].temperature);
    let flux_residual = actual_heat - flux;
    let scale = endpoint[0]
        .change
        .abs()
        .max(actual_heat.abs())
        .max(flux.abs());
    if !energy_residual.is_finite()
        || !flux_residual.is_finite()
        || energy_residual.abs() > control.tolerance(scale)
        || flux_residual.abs() > control.tolerance(scale)
    {
        return Err(WallExchangeError::Balance {
            energy_residual_j: energy_residual,
            thermal_residual_j: flux_residual,
        });
    }
    receipt.endpoint_state = endpoint.each_ref().map(|p| p.thermal);
    receipt.endpoint_temperature_k = endpoint.each_ref().map(|p| p.temperature);
    receipt.heat_a_to_b_j = actual_heat;
    receipt.wall_energy_change_j = endpoint.each_ref().map(|p| p.change);
    receipt.energy_residual_j = energy_residual;
    receipt.flux_residual_j = flux_residual;
    receipt.iterations = iterations;
    poll(cancelled)?;
    Ok(receipt)
}
