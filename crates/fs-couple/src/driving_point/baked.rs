//! Full-band mouth loads on the existing time-domain line realization.
//! Positive bins come only from the table. Geometry, gas and the DC limit are
//! explicit; the original Bessel propagation and FIR/ladder remain the owners.
//! The bore/wall law must also admit the existing 1e-3 rad/s quasi-static
//! probe; mouth DC metadata cannot authorize out-of-domain wall material data.

use super::{
    DrivingPointError, LineDc, LineLadder, LineRealization, PASSIVITY_TOL, PoleState,
    build_impedance_response_filter, bulk_delay_pole_state,
    characteristic_line_response_with_poles, duct_axial_length, impedance_reflectance,
};
use fs_duct::{Duct, DuctError, DuctResponse, LossModel, TabulatedLoad, WallModel};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_vfit::discretize::DelayedFilter;

pub(super) fn matched_refusal() -> DrivingPointError {
    DrivingPointError::Invalid {
        what: "Matched is a multimodal continuation; scalar lines require an analytic or referenced tabulated mouth",
    }
}

pub(super) fn validate_band(
    table: &TabulatedLoad,
    rate: u32,
    n_fft: usize,
) -> Result<(), DrivingPointError> {
    if rate == 0 {
        return Err(DrivingPointError::Invalid {
            what: "positive sample rate required",
        });
    }
    let dt = 1.0 / f64::from(rate);
    for bin in [1, n_fft / 2] {
        table
            .z_at(core::f64::consts::TAU * bin as f64 / (n_fft as f64 * dt))
            .map_err(DrivingPointError::Duct)?;
    }
    Ok(())
}

pub(super) fn response(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    load: C64,
    wall: WallModel<'_>,
) -> Result<DuctResponse, DuctError> {
    fs_duct::input_impedance_load_with_wall(duct, gas, omega, LossModel::Bessel, load, wall)
}

pub(super) fn reflection(
    duct: &Duct,
    gas: &GasState,
    omega: f64,
    load: C64,
    wall: WallModel<'_>,
    zc: f64,
) -> Result<C64, DuctError> {
    let r = match response(duct, gas, omega, load, wall) {
        Ok(value) => impedance_reflectance(&value, zc),
        Err(DuctError::Singular) => fs_duct::compose::input_reflection_load_with_wall(
            duct,
            gas,
            omega,
            LossModel::Bessel,
            load,
            wall,
        )?
        .referenced_to(zc)?,
        Err(error) => return Err(error),
    };
    if !r.re.is_finite() || !r.im.is_finite() || r.abs() > 1.0 + PASSIVITY_TOL {
        return Err(DuctError::BadParameter {
            what: "tabulated characteristic response is not finite and passive",
        });
    }
    Ok(r)
}

pub(super) fn characteristic(
    duct: &Duct,
    gas: &GasState,
    table: &TabulatedLoad,
    rate: u32,
    n: usize,
    zc: f64,
    wall: Option<&WallPin>,
    n_fft: Option<usize>,
    ladder: Option<&LineLadder>,
) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
    let length = duct_axial_length(duct);
    let radius = duct
        .segments
        .last()
        .ok_or(DrivingPointError::Duct(DuctError::EmptyDuct))?
        .outlet_radius();
    let mouth = table
        .validate_mouth(length, radius, gas)
        .map_err(DrivingPointError::Duct)?;
    let default_ladder = LineLadder::default();
    let ladder = ladder.unwrap_or(&default_ladder);
    let band_n = match n_fft {
        Some(n) => n.clamp(256, 8192).next_power_of_two(),
        None => super::validate_ladder(ladder)?,
    };
    // Admit the complete requested ladder before constructing any filter.
    validate_band(table, rate, band_n)?;
    if !zc.is_finite() || zc <= 0.0 {
        return Err(DrivingPointError::Invalid {
            what: "finite positive port impedance required",
        });
    }
    let dt = 1.0 / f64::from(rate);
    let delay = 2.0 * length / gas.sound_speed / dt;
    if !(delay >= 2.0 && delay < n as f64 - 2.0) {
        return Err(DrivingPointError::Invalid {
            what: "round-trip delay does not fit the realized history",
        });
    }
    let wall = WallModel::from_pin(wall);
    // Same quasi-static bore limit as the analytic ladder; the mouth load is
    // its DECLARED DC value, never an out-of-band table query or fitted splice.
    let dc = reflection(
        duct,
        gas,
        1.0e-3,
        C64::new(mouth.dc_impedance, 0.0),
        wall,
        zc,
    )
    .map_err(DrivingPointError::Duct)?
    .re;
    let poles: Vec<PoleState> = vec![bulk_delay_pole_state(delay)];
    let sample = |omega| reflection(duct, gas, omega, table.z_at(omega)?, wall, zc);
    if n_fft.is_none() {
        super::realize_characteristic_response(
            delay,
            LineDc::Explicit(dc),
            rate,
            zc,
            poles,
            ladder,
            sample,
        )
    } else {
        characteristic_line_response_with_poles(
            delay,
            LineDc::Explicit(dc),
            rate,
            n,
            zc,
            n_fft,
            poles,
            sample,
        )
    }
}

pub(super) fn impedance(
    duct: &Duct,
    gas: &GasState,
    table: &TabulatedLoad,
    rate: u32,
    n: usize,
    wall: Option<&WallPin>,
) -> Result<DelayedFilter, DrivingPointError> {
    let radius = duct
        .segments
        .last()
        .ok_or(DrivingPointError::Duct(DuctError::EmptyDuct))?
        .outlet_radius();
    let mouth = table
        .validate_mouth(duct_axial_length(duct), radius, gas)
        .map_err(DrivingPointError::Duct)?;
    let n_fft = n.clamp(256, 8192).next_power_of_two();
    validate_band(table, rate, n_fft)?;
    let wall = WallModel::from_pin(wall);
    let dc = response(duct, gas, 1.0e-3, C64::new(mouth.dc_impedance, 0.0), wall)
        .map_err(DrivingPointError::Duct)?
        .p_mouth_over_u_in
        .re;
    build_impedance_response_filter(1.0 / f64::from(rate), n_fft, dc, 0.0, &[], |omega| {
        response(duct, gas, omega, table.z_at(omega)?, wall).map(|value| value.p_mouth_over_u_in)
    })
}

#[cfg(test)]
mod tests;
