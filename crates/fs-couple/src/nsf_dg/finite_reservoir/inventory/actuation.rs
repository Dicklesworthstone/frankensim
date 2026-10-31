//! Source-owned implicit midpoint in (specific internal energy, volume).
//! NASA9 delegates to its existing parcel owner. The named mixture has constant
//! Cv, so the same conservative midpoint equation has an analytic inverse.
//! These are tentative closed-mass states, not a new breath or heat-transfer law.

use super::{Inventory, NsfDgError, NsfFiniteReservoirActuationFrame,
    NsfFiniteReservoirLimits, material, poll};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_material::thermochemical_gas::energy::parcel::midpoint::GasHeatVolumeStep;

pub(super) fn candidate(
    old: &Inventory, volume_m3: f64, heat_j: f64,
    limits: &NsfFiniteReservoirLimits, gate: &CancelGate,
) -> Result<(Inventory, NsfFiniteReservoirActuationFrame, GasState), NsfDgError> {
    poll(gate)?;
    let initial_volume_m3 = old.volume_m3();
    let mass = old.mass_kg();
    let vm = 0.5 * initial_volume_m3 + 0.5 * volume_m3;
    let dv = volume_m3 - initial_volume_m3;
    if ![initial_volume_m3, volume_m3, vm, mass].into_iter()
        .all(|v| v.is_finite() && v > 0.0) || !heat_j.is_finite()
    {
        return Err(NsfDgError::Input("finite reservoir actuation requires positive finite volumes/mass and finite heat"));
    }
    let unchanged = dv == 0.0 && heat_j == 0.0;
    let (next, midpoint, original_step) = match old {
            Inventory::Thermochemical { gas, energy_solve, energy_reference_j_kg } => {
                let candidate = gas.midpoint_heat_volume(volume_m3, heat_j, *energy_solve,
                    || gate.is_requested()).map_err(material)?;
                let midpoint = *candidate.midpoint().state();
                let (gas, step) = candidate.into_parts();
                (Inventory::Thermochemical { gas, energy_solve: *energy_solve,
                    energy_reference_j_kg: *energy_reference_j_kg }, midpoint, Some(step))
            }
            Inventory::NamedMixture { .. } if unchanged => {
                // Keep the admitted source and passive storage bits. NASA's
                // unchanged branch above also retains its original receipt.
                (old.clone(), *old.state(), None)
            }
            Inventory::NamedMixture { source, gas, temperature_window_k, .. } => {
                let reference = source.reference_state();
                let cv = reference.specific_heat_cp - reference.specific_gas_constant;
                let r = reference.specific_gas_constant;
                let q = heat_j / mass;
                let a = 0.5 * (dv / vm) * r;
                let denominator = cv + a;
                if !q.is_finite() || (heat_j != 0.0 && q == 0.0)
                    || !denominator.is_finite() || denominator <= 0.0
                {
                    return Err(NsfDgError::Input("finite reservoir named heat/work increment is unrepresentable"));
                }
                // Cv*(Tm-T0) + pm*dV/(2m) = Q/(2m), with pm=m*R*Tm/Vm.
                // Reconstruct the endpoint from du(T0,T1)=2*du(T0,Tm).
                let delta_mid = (0.5 * q - a * gas.temperature) / denominator;
                let tm = gas.temperature + delta_mid;
                let t1 = gas.temperature + 2.0 * delta_mid;
                let midpoint = source.at(tm, mass / vm * r * tm).map_err(NsfDgError::Flow)?;
                poll(gate)?;
                let endpoint = source.at(t1, mass / volume_m3 * r * t1).map_err(NsfDgError::Flow)?;
                (Inventory::NamedMixture { source: source.clone(), gas: endpoint,
                    volume_m3, mass_kg: mass, temperature_window_k: *temperature_window_k }, midpoint, None)
            }
    };
    if !unchanged { admit_path(old, &midpoint, &next, gate)?; }

    // The retained mass is exact; independently audit EOS coherence at both
    // material states. Never silently replace conserved mass with rounded rho*V.
    let midpoint_mass_residual = midpoint.density * vm - mass;
    let mass_residual_kg = next.state().density * volume_m3 - mass;
    if next.mass_kg().to_bits() != mass.to_bits()
        || ![midpoint_mass_residual, mass_residual_kg].into_iter()
            .all(|v| v.is_finite() && v.abs() <= limits.mass_absolute_tolerance_kg)
    {
        return Err(NsfDgError::Input("finite reservoir actuation did not preserve closed mass and EOS density"));
    }
    let du = old.specific_energy_change(old.state().temperature, next.state().temperature,
        "finite reservoir actuation endpoint energy is invalid")?;
    let du_mid = old.specific_energy_change(old.state().temperature, midpoint.temperature,
        "finite reservoir actuation midpoint energy is invalid")?;
    let stored = mass * du;
    let work = midpoint.pressure * dv;
    let requested = heat_j - work;
    let residual = (stored + work) - heat_j;
    let stage_residual = mass * (2.0 * du_mid - du);
    let allowance = limits.energy_absolute_tolerance_j
        + limits.relative_tolerance * (stored.abs() + requested.abs());
    let stage_allowance = limits.energy_absolute_tolerance_j
        + limits.relative_tolerance * ((2.0 * mass * du_mid).abs() + stored.abs());
    if ![stored, work, requested, residual, stage_residual, allowance, stage_allowance]
        .into_iter().all(f64::is_finite)
        || (dv != 0.0 && work == 0.0)
        || (requested != 0.0 && (stored == 0.0 || stored.signum() != requested.signum()))
        || residual.abs() > allowance || stage_residual.abs() > stage_allowance
    {
        return Err(NsfDgError::Input("finite reservoir actuation failed original-law midpoint heat/work balance"));
    }
    let step = if let Some(step) = original_step {
        // Preserve every published NASA receipt bit, including its original
        // entropy increment and derivatives; the checks above add no debit.
        step
    } else {
        let reference = old.state();
        let cv = reference.specific_heat_cp - reference.specific_gas_constant;
        let r = reference.specific_gas_constant;
        let denominator = cv + 0.5 * (dv / vm) * r;
        // Frozen ideal-mixture entropy difference; no formation/reference
        // entropy or changing-composition term enters this closed-mass step.
        let entropy = mass * (cv * ((next.state().temperature - reference.temperature)
            / reference.temperature).ln_1p() + r * (dv / initial_volume_m3).ln_1p());
        // Differentiate the same discrete pm=N/B, B=Cv*Vm+R*dV/2.
        let dpdv = -0.5 * ((cv + r) / denominator) * (midpoint.pressure / vm);
        let dpdq = 0.5 * (r / denominator) / vm;
        if ![denominator, entropy, dpdv, dpdq].into_iter().all(f64::is_finite)
            || denominator <= 0.0 || dpdv >= 0.0 || dpdq <= 0.0
        {
            return Err(NsfDgError::Input("finite reservoir named entropy or discrete pressure derivatives are unrepresentable"));
        }
        GasHeatVolumeStep {
            initial_volume_m3, final_volume_m3: volume_m3, supplied_heat_j: heat_j,
            work_by_gas_j: work, stored_energy_change_j: stored, energy_residual_j: residual,
            entropy_change_j_k: entropy, midpoint_energy_residual_j_kg: 2.0 * du_mid - du,
            pressure_final_volume_derivative_pa_m3: dpdv, pressure_heat_derivative_pa_j: dpdq,
            iterations: 0,
        }
    };
    poll(gate)?;
    Ok((next, NsfFiniteReservoirActuationFrame { midpoint_pressure_pa: midpoint.pressure, step }, midpoint))
}

/// Bound the conservative segment between the old, midpoint and final states.
/// The named source's constant Cv makes T affine and p a ratio of affine
/// functions in (u,V), so their extrema are at the endpoints. Tmin/Pmax also
/// conservatively bounds its original monotone Buck phase admission. NASA9 may
/// have variable Cv: use the wider independent positive-T/volume enclosure.
/// Every query still goes through the retained original material owner.
fn admit_path(
    old: &Inventory, midpoint: &GasState, next: &Inventory, gate: &CancelGate,
) -> Result<(), NsfDgError> {
    let first = old.state();
    let last = next.state();
    let t_low = first.temperature.min(midpoint.temperature).min(last.temperature);
    let t_high = first.temperature.max(midpoint.temperature).max(last.temperature);
    let window = old.temperature_window_k();
    if ![t_low, t_high].into_iter().all(|t| t.is_finite() && t >= window[0] && t <= window[1]) {
        return Err(NsfDgError::Input("finite reservoir actuation leaves the NSF temperature window; refine the interval"));
    }
    poll(gate)?;
    match old {
        Inventory::NamedMixture { source, .. } => {
            let p_low = first.pressure.min(midpoint.pressure).min(last.pressure);
            let p_high = first.pressure.max(midpoint.pressure).max(last.pressure);
            source.at(t_low, p_high).map_err(NsfDgError::Flow)?;
            source.at(t_high, p_low).map_err(NsfDgError::Flow)?;
        }
        Inventory::Thermochemical { gas, .. } => {
            let v_low = old.volume_m3().min(next.volume_m3());
            let v_high = old.volume_m3().max(next.volume_m3());
            let r = first.specific_gas_constant;
            let p_low = (((old.mass_kg() / v_high).next_down() * r).next_down() * t_low).next_down();
            let p_high = (((old.mass_kg() / v_low).next_up() * r).next_up() * t_high).next_up();
            for t in [t_low, t_high] {
                for p in [p_low, p_high] {
                    poll(gate)?;
                    let point = gas.point().clone().with("temperature", t)
                        .and_then(|point| point.with("pressure", p))
                        .map_err(|_| NsfDgError::Input("finite reservoir actuation domain coordinates are invalid"))?;
                    gas.model().evaluate(&point).map_err(|error| material(error.into()))?;
                }
            }
        }
    }
    poll(gate)
}
