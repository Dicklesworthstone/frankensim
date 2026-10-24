//! One caloric account for the original strip and spatial beam quadrature.
//!
//! The solid owner has already advanced both material points and integrated
//! their extensive free-energy/work/loss receipt. Bind those exact histories
//! to the original caloric source; do not advance a point again or infer a
//! second constitutive law from the section's aggregate stiffness.

use std::collections::BTreeMap;
use fs_material::thermovisco::{Energy, thermal::CaloricPoint};
use fs_solid::deposition::Layer;
use fs_thermal::enthalpy::{CellId, DepositedCell, DepositionError};
use super::{LayerThermalBinding, StripEnthalpyError};

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct LayerThermalExchange {
    pub(super) mechanical_internal_energy_change_j: f64,
    pub(super) mechanical_entropy_energy_change_j: f64,
    pub(super) reversible_heat_j: f64,
    pub(super) dissipated_j: f64,
    pub(super) total_material_heat_j: f64,
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() { Ok(value) } else { Err(DepositionError::Arithmetic) }
}

/// Evaluate the same two equal Gauss weights as `LayeredStrip`. The native
/// `energy` is already integrated over this layer's full reference volume.
/// `CaloricPoint` owns admission of the affine modulus and constant free strain,
/// and its getters include the actual retained branch history at each point.
pub(super) fn layer_exchange<E>(
    old: &Layer<'_>, next: &Layer<'_>, old_heat: &DepositedCell,
    next_heat: &DepositedCell, binding: &LayerThermalBinding<'_>,
    length_m: f64, energy: Energy,
) -> Result<LayerThermalExchange, StripEnthalpyError<E>> {
    if old.geometry != next.geometry || old.birth_time_s != next.birth_time_s
        || old_heat.id != binding.cell || next_heat.id != binding.cell
        || old_heat.material_id != next_heat.material_id
        || old_heat.mass_kg != next_heat.mass_kg
        || old_heat.geometry.volume_m3 != next_heat.geometry.volume_m3
        || old_heat.geometry.surface_area_m2 != next_heat.geometry.surface_area_m2
        || !length_m.is_finite() || length_m <= 0.0 {
        return Err(StripEnthalpyError::Invalid("material heat requires unchanged reference layer and thermal geometry"));
    }
    let volume = finite(old.geometry.area_m2 * length_m)?;
    if !volume.is_finite() || volume <= 0.0
        || (volume-old_heat.geometry.volume_m3).abs()
            > 64.0*f64::EPSILON*volume.abs().max(old_heat.geometry.volume_m3.abs()) {
        return Err(StripEnthalpyError::Invalid("material heat uses the original layer reference volume"));
    }
    if !energy.dissipated.is_finite() || energy.dissipated < 0.0
        || !energy.modulus_work.is_finite() || energy.thermal_strain_work != 0.0 {
        return Err(StripEnthalpyError::Invalid("invalid material loss or unsupported thermal eigenstrain work"));
    }
    let mut internal_change = 0.0;
    let mut entropy_energy_change = 0.0;
    for (before, after) in old.material_points().iter().zip(next.material_points()) {
        if !std::ptr::eq(before.law(), binding.material.law())
            || !std::ptr::eq(after.law(), binding.material.law()) {
            return Err(StripEnthalpyError::Invalid("thermal exchange must retain the original mechanical source"));
        }
        let bind = |point, heat: &DepositedCell| CaloricPoint::try_new(point,
            binding.caloric.curve(), heat.specific_h_j_kg, 0.5*heat.mass_kg,
            0.5*volume, 0.5*heat.geometry.surface_area_m2)
            .map_err(StripEnthalpyError::CaloricAdmission);
        let before = bind(before.clone(), old_heat)?;
        let after = bind(after.clone(), next_heat)?;
        internal_change = finite(internal_change + finite(after.mechanical_internal_energy_j()
            - before.mechanical_internal_energy_j())?)?;
        let old_ts = finite(before.mechanical().temperature_k()*before.mechanical_entropy_j_k())?;
        let next_ts = finite(after.mechanical().temperature_k()*after.mechanical_entropy_j_k())?;
        entropy_energy_change = finite(entropy_energy_change + finite(next_ts-old_ts)?)?;
    }
    // Native identity: Δpsi + D + N = W + Cpsi. Hence Qrev=-Cpsi-Δ(T*s)
    // closes ΔU + ΔH + N = W + Qexternal without thermalizing numerical loss.
    let reversible_heat = finite(-energy.modulus_work-entropy_energy_change)?;
    Ok(LayerThermalExchange {
        mechanical_internal_energy_change_j: internal_change,
        mechanical_entropy_energy_change_j: entropy_energy_change,
        reversible_heat_j: reversible_heat,
        dissipated_j: energy.dissipated,
        total_material_heat_j: finite(energy.dissipated+reversible_heat)?,
    })
}

/// Insert signed physical material heat into the actual original cell source.
/// Return its representable energy increment, including rounding in source
/// addition. A small contribution may not silently disappear beside other heat.
pub(super) fn deposit_heat_j<E>(
    sources: &mut BTreeMap<CellId, f64>, cell: CellId,
    heat_j: f64, dt: f64, tolerance_j: f64,
) -> Result<f64, StripEnthalpyError<E>> {
    if !dt.is_finite() || dt <= 0.0 || !tolerance_j.is_finite() || tolerance_j < 0.0 {
        return Err(StripEnthalpyError::Invalid("finite positive duration and nonnegative heat tolerance required"));
    }
    let power = finite(heat_j/dt)?;
    if heat_j != 0.0 && power == 0.0 {
        return Err(StripEnthalpyError::Invalid("material heat vanished in power conversion"));
    }
    let other = sources.entry(cell).or_insert(0.0);
    let merged = finite(*other+power)?;
    let increment = finite(merged-*other)?;
    let delivered = finite(dt*increment)?;
    if heat_j != 0.0 && (delivered == 0.0 || (delivered > 0.0) != (heat_j > 0.0)) {
        return Err(StripEnthalpyError::Invalid("material heat vanished or changed sign in the receiving cell"));
    }
    if finite(delivered-heat_j)?.abs() > tolerance_j {
        return Err(StripEnthalpyError::Invalid("actual receiving-cell source increment does not preserve material heat"));
    }
    *other = merged;
    Ok(delivered)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_heat_uses_the_actual_source_increment_and_refuses_erasure() {
        let cell = CellId(7);
        let mut sources = BTreeMap::from([(cell,2.0)]);
        let delivered = deposit_heat_j::<core::convert::Infallible>(
            &mut sources,cell,-0.25,0.2,1e-14).unwrap();
        assert_eq!(delivered,-0.25);
        assert_eq!(sources[&cell],0.75);
        let mut overwhelming = BTreeMap::from([(cell,1e30)]);
        let before = overwhelming.clone();
        assert!(deposit_heat_j::<core::convert::Infallible>(
            &mut overwhelming,cell,-1e-4,0.2,1e-9).is_err());
        assert_eq!(overwhelming,before);
    }
}
