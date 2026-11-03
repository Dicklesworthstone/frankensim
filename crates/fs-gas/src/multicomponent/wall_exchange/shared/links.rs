//! Conservative finite thermal ports on the original wall inventories.
//!
//! Each source owner supplies the two current stage states and owns ordered
//! composition and publication. Frozen and evolving-mixture networks use the
//! same original scalar equation; no capacity or latent chart is reconstructed.

use super::*;

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
    link: SharedWallLink,
    duration: Time,
    control: WallControl,
    work: &Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SharedWallLinkReport, WallExchangeError> {
    poll(cancelled)?;
    let states = [a.candidate, b.candidate];
    let initial = [
        a.after_state(states[0], 0.0, work)?,
        b.after_state(states[1], 0.0, work)?,
    ];
    let temperatures = initial.each_ref().map(|p| p.temperature);
    let integrated = duration.value() * link.conductance_w_per_k;
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
    let (mut low, mut high) = if initial_flux > 0.0 {
        (0.0, initial_flux.min(la[1]).min(-lb[0]))
    } else {
        (initial_flux.max(la[0]).max(-lb[1]), 0.0)
    };
    if low < 0.0 {
        low = low.next_up();
    }
    if high > 0.0 {
        high = high.next_down();
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
            let flux = integrated * (wa.temperature - wb.temperature);
            let residual = heat - flux;
            let slope = wa
                .slope
                .zip(wb.slope)
                .map(|(a, b)| 1.0 + integrated * (a + b));
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
