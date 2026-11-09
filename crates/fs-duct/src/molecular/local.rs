//! Shared local caloric/wall law for the original plane and modal propagators.
//! The circular thermal measure acts on the original dynamic Cp; a separate
//! equilibrium thermal shunt would count the same heat exchange twice.

use super::{MolecularDuctError, bad, finite};
use crate::{DuctError, LossModel, SegmentWave};
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};
use fs_phs::ZkDiffusive;

pub(crate) fn caloric(
    state: &RelaxationState,
    omega: f64,
) -> Result<(C64, C64), MolecularDuctError> {
    let response = state.caloric_response(omega)?;
    let k = response.bulk_modulus_pa.scale(1.0 / state.gas().density)
        .sqrt().recip().scale(omega);
    if !finite(k) || k.re <= 0.0 || k.im < 0.0 {
        return Err(bad("molecular duct material wavenumber is unrepresentable").into());
    }
    Ok((response.cp_j_kg_k, k))
}

pub(crate) fn wave(
    state: &RelaxationState,
    radius: f64,
    x_m: f64,
    wall: Option<&ZkDiffusive>,
    omega: f64,
    cp: C64,
    relaxing_k: C64,
) -> Result<SegmentWave, DuctError> {
    let gas = state.gas();
    let area = core::f64::consts::PI * radius * radius;
    let Some(model) = wall else {
        if state.processes().is_empty() {
            return crate::segment_wave(gas, radius, omega, LossModel::Lossless);
        }
        let k = relaxing_k;
        let z = C64::from_re(omega * gas.density) / k;
        return Ok(SegmentWave {
            wavenumber: k,
            characteristic_impedance: z.scale(1.0 / area),
            specific_impedance: z,
            shear_number: f64::INFINITY,
        });
    };
    let (series, equilibrium_shunt) = model.response_at(omega)
        .map_err(|reason| DuctError::WallDiffusionUnresolved {
            reason, radius_m: radius, x_m,
        })?;
    let z_plus = series.scale(gas.density / area);
    let y_plus = if state.processes().is_empty() {
        equilibrium_shunt.scale(area / (gas.density * gas.sound_speed * gas.sound_speed))
    } else {
        let mut temperature = C64::ZERO;
        for term in model.shunt_terms() {
            let weight = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
            // s/(s cp+d) = 1/(cp+i d/omega), avoiding large s*cp.
            let denominator =
                cp + C64::new(0.0, gas.specific_heat_cp * (term.rate_per_s / omega));
            temperature = temperature + denominator.recip().scale(weight);
        }
        let kappa = (C64::ONE - temperature.scale(gas.specific_gas_constant))
            .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature));
        (C64::new(0.0, -omega) * kappa.scale(area)).conj()
    };
    let (k, zc) = crate::helmholtz_wave_from_telegraph(z_plus, y_plus);
    if !finite(k) || k.re <= 0.0 || k.im < 0.0 || !finite(zc) || zc.re <= 0.0 {
        return Err(bad("molecular duct wave is not finite, forward and passive"));
    }
    Ok(SegmentWave {
        wavenumber: k,
        characteristic_impedance: zc,
        specific_impedance: zc.scale(area),
        shear_number: radius * det::sqrt(gas.density * omega / gas.dynamic_viscosity),
    })
}
