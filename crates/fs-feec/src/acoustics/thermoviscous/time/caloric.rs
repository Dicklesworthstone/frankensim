//! Exact caloric memory, without eliminating the resolved thermal field.
//!
//! With theta = T'/T0 and eta_j = T'_j/T0, each original reservoir obeys
//! rho C_j T0 M eta'_j = rho C_j T0/tau_j M (theta - eta_j).
//! Eliminating eta_j at s = i omega gives the ORIGINAL
//! Cp(s) = Cp_frozen + sum C_j/(1+s tau_j). Its positive stored energy is
//! rho C_j T0 eta_j^T M eta_j/2; exchange loss is
//! rho C_j T0/tau_j (theta-eta_j)^T M (theta-eta_j).
//!
//! M is the original FULL hierarchical P2 mass. Only theta is restricted by
//! the translational thermal boundary. Reservoirs have no spatial conduction
//! and no separate Dirichlet condition: a molecule at an isothermal wall can
//! still start out of internal equilibrium and release its stored energy.

use super::*;

/// One original material heat reservoir; not a fit or a bulk viscosity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermoviscousCaloricReservoir {
    /// Nonnegative mass-specific relaxing capacity [J/(kg K)].
    pub heat_capacity_j_kg_k: f64,
    /// Positive molecular heat-exchange time [s].
    pub relaxation_time_s: f64,
}

/// Explicit bounded constitutive data, independent of the smaller state cap.
const MAX_CALORIC_RESERVOIRS: usize = 32;

pub(super) fn layout(
    model: &ThermoviscousFe<'_>,
    frozen_cp: f64,
    reservoirs: &[ThermoviscousCaloricReservoir],
) -> Result<(usize, Vec<Option<usize>>), ThermoviscousTimeError> {
    let minimum_cp = if model.ideal_gas {
        model.medium.gas_constant_j_kg_k()
    } else {
        model.thermodynamics.minimum_cp(model.medium)
    };
    if reservoirs.len() > MAX_CALORIC_RESERVOIRS
        || !frozen_cp.is_finite()
        || frozen_cp <= minimum_cp
    {
        return Err(invalid(
            "bounded original caloric spectrum and stable frozen storage required",
        ));
    }
    let mut equilibrium_cp = frozen_cp;
    let mut mixed = model.n;
    let mut offsets = Vec::with_capacity(reservoirs.len());
    for reservoir in reservoirs {
        let c = reservoir.heat_capacity_j_kg_k;
        let tau = reservoir.relaxation_time_s;
        if !c.is_finite() || c < 0.0 || !tau.is_finite() || tau <= 0.0 {
            return Err(invalid(
                "finite nonnegative reservoir capacity and positive time required",
            ));
        }
        equilibrium_cp += c;
        offsets.push(if c > 0.0 {
            let start = mixed;
            mixed = mixed
                .checked_add(model.space.ndof)
                .ok_or_else(|| invalid("complete molecular state dimension overflow"))?;
            Some(start)
        } else {
            None
        });
    }
    let tolerance = 128.0 * f64::EPSILON * (reservoirs.len() + 1) as f64 * model.medium.cp_j_kg_k;
    if !equilibrium_cp.is_finite() || (equilibrium_cp - model.medium.cp_j_kg_k).abs() > tolerance {
        return Err(invalid(
            "original equilibrium Cp must equal frozen Cp plus every reservoir",
        ));
    }
    Ok((mixed, offsets))
}

pub(super) fn assemble(
    model: &ThermoviscousFe<'_>,
    frozen_cp: f64,
    reservoirs: &[ThermoviscousCaloricReservoir],
    offsets: &[Option<usize>],
    mixed: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<f64>, Vec<f64>), ThermoviscousTimeError> {
    let base = model.n;
    let mut mass = vec![0.0; mixed * mixed];
    let mut exchange = vec![0.0; mixed * mixed];
    for i in 0..base {
        poll(cancelled)?;
        mass[i * mixed..i * mixed + base].copy_from_slice(&model.mass[i * base..(i + 1) * base]);
    }
    if frozen_cp != model.medium.cp_j_kg_k {
        for i in model.temperature_start..base {
            for j in model.temperature_start..base {
                // Same original thermal mass as harmonic complex-Cp assembly.
                mass[i * mixed + j] = model.mass[i * base + j] / model.medium.cp_j_kg_k * frozen_cp;
            }
        }
    }
    if offsets.iter().all(Option::is_none) {
        return Ok((mass, exchange));
    }
    let ndof = model.space.ndof;
    let mut full_mass = vec![0.0; ndof * ndof];
    let quad = duffy_quadrature(4);
    for t in 0..model.space.complex.tets.len() {
        poll(cancelled)?;
        let dofs = model.space.element_dofs(t);
        let jac = 6.0 * model.geometry.vol_signed[t].abs();
        for &(bary, w) in &quad {
            let weight = jac * w;
            let values: Vec<_> = dofs.iter().map(|(_, f)| f.eval(bary, 2)).collect();
            for (i, &(gi, _)) in dofs.iter().enumerate() {
                for (j, &(gj, _)) in dofs.iter().enumerate() {
                    full_mass[gi * ndof + gj] += weight * values[i] * values[j];
                }
            }
        }
    }
    for (reservoir, &offset) in reservoirs.iter().zip(offsets) {
        let Some(start) = offset else { continue };
        let capacity = model.medium.density_kg_m3
            * model.medium.temperature_k
            * reservoir.heat_capacity_j_kg_k;
        let conductance = capacity / reservoir.relaxation_time_s;
        if !capacity.is_finite()
            || capacity <= 0.0
            || !conductance.is_finite()
            || conductance <= 0.0
        {
            return Err(invalid(
                "original molecular storage or exchange is unrepresentable",
            ));
        }
        for i in 0..ndof {
            poll(cancelled)?;
            let row = start + i;
            let thermal_i = model.temperature_map[i].map(|k| model.temperature_start + k);
            for j in 0..ndof {
                let col = start + j;
                let m = full_mass[i * ndof + j];
                let loss = conductance * m;
                mass[row * mixed + col] = capacity * m;
                exchange[row * mixed + col] += loss;
                if let Some(ti) = thermal_i {
                    exchange[ti * mixed + col] -= loss;
                    if let Some(tj) = model.temperature_map[j] {
                        exchange[ti * mixed + model.temperature_start + tj] += loss;
                    }
                }
                if let Some(tj) = model.temperature_map[j] {
                    exchange[row * mixed + model.temperature_start + tj] -= loss;
                }
            }
        }
    }
    if mass.iter().chain(&exchange).any(|value| !value.is_finite()) {
        return Err(invalid("nonfinite original caloric field assembly"));
    }
    Ok((mass, exchange))
}
