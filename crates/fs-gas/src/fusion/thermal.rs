//! Shared original-source thermal kernels for the admitted fusion integrators.
//!
//! The caller owns the plasma/source admission and the thermal-pool accounts.
//! These allocation-free kernels return coefficients or transferred energy;
//! they do not advance time, change composition or publish a material state.

use super::{DtFusionError, DtFusionTransport, poll};
use fs_math::det::{expm1, sqrt};
use fs_qty::{Amount, Temperature, Volume};
use fs_thermochem::coulomb::CoulombErrorV1;
use fs_thermochem::plasma_radiation::BremsstrahlungErrorV1;

#[derive(Clone, Copy)]
pub(super) struct ThermalSources {
    pub(super) temperatures: [f64; 2],
    pub(super) conductance_w_k: f64,
    pub(super) radiation_w_per_sqrt_k: f64,
}

/// Refresh the selected thermal laws at the actual admitted populations and
/// temperatures. The owning integrator checks their original-bank identity.
pub(super) fn evaluate(
    transport: &DtFusionTransport,
    amounts: &[Amount],
    volume: Volume,
    temperatures: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ThermalSources, DtFusionError> {
    let [th, te] = temperatures.map(Temperature::new);
    let conductance_w_k = if let Some(model) = &transport.coulomb {
        model
            .evaluate_point(amounts, volume, th, te, &mut *cancelled)
            .map_err(coulomb_error)?
            .electron_heavy_thermal_conductance_w_m3_k()
            * volume.value()
    } else {
        0.0
    };
    let radiation_w_per_sqrt_k = if let Some(model) = &transport.bremsstrahlung {
        model
            .evaluate_point(amounts, volume, th, te, &mut *cancelled)
            .map_err(radiation_error)?
            .total_power_over_sqrt_temperature_w_k_sqrt()
    } else {
        0.0
    };
    if [conductance_w_k, radiation_w_per_sqrt_k]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(DtFusionError::Arithmetic(
            "finite thermal source coefficients",
        ));
    }
    poll(cancelled)?;
    Ok(ThermalSources {
        temperatures,
        conductance_w_k,
        radiation_w_per_sqrt_k,
    })
}

/// Signed electron-to-heavy energy for admitted positive temperatures and
/// heat capacities [J/K]. With frozen G [W/K], the temperature difference
/// decays as exp[-G(1/Ch+1/Ce)dt]; the reciprocal update cannot overshoot
/// equilibrium even for an exposure much larger than one.
pub(super) fn exchange_energy(
    temperatures: [f64; 2],
    capacities: [f64; 2],
    dt: f64,
    conductance: f64,
) -> Result<f64, DtFusionError> {
    if conductance == 0.0 || dt == 0.0 {
        return Ok(0.0);
    }
    let [th, te] = temperatures;
    let [ch, ce] = capacities;
    let exposure = (conductance / ch + conductance / ce) * dt;
    if exposure.is_nan() || exposure < 0.0 {
        return Err(DtFusionError::Arithmetic("Coulomb thermal exposure"));
    }
    // ch*ce/(ch+ce), without multiplying two extensive capacities.
    let reduced_capacity = ch.min(ce) / (1.0 + ch.min(ce) / ch.max(ce));
    let transfer = reduced_capacity * (te - th) * (-expm1(-exposure));
    if !transfer.is_finite() {
        return Err(DtFusionError::Arithmetic(
            "reciprocal Coulomb energy transfer",
        ));
    }
    Ok(transfer)
}

/// Escaping photon energy [J] for an admitted electron pool, freezing its
/// instantaneous coefficient C=P/sqrt(Te). Then sqrt(Te)'=-C/(2Ce).
/// Avoid subtraction of nearly equal squared temperatures; crossing zero
/// refuses the step instead of clipping the pool or extending a source domain.
pub(super) fn radiated_energy(
    te: f64,
    ce: f64,
    dt: f64,
    coefficient: f64,
) -> Result<f64, DtFusionError> {
    if coefficient == 0.0 || dt == 0.0 {
        return Ok(0.0);
    }
    let root = sqrt(te);
    let decrement = (0.5 * coefficient / ce) * dt;
    if !decrement.is_finite() || decrement >= root {
        return Err(DtFusionError::OutsideDomain);
    }
    let loss = (coefficient * dt) * (root - 0.5 * decrement);
    if !loss.is_finite() || loss < 0.0 {
        return Err(DtFusionError::Arithmetic(
            "escaping free-free photon energy",
        ));
    }
    Ok(loss)
}

fn coulomb_error(error: CoulombErrorV1) -> DtFusionError {
    match error {
        CoulombErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::Coulomb(error),
    }
}

fn radiation_error(error: BremsstrahlungErrorV1) -> DtFusionError {
    match error {
        BremsstrahlungErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::Radiation(error),
    }
}
