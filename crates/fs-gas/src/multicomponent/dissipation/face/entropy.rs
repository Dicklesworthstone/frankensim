//! Source-conjugate axial molecular faces for the opt-in entropy-admitted path.
//!
//! This is a two-point discretization of the ORIGINAL ideal-mixture calorics
//! and Maxwell--Stefan/Onsager law, not a constant-gamma replacement. For
//! beta=1/T, use the divided source enthalpy
//!   h_i^* = Delta(h_i/T - s_i^0) / Delta beta.
//! It is formed from sensible source intervals about the sampled temperature,
//! before formation references can cancel in floating-point arithmetic. Thus
//! a_i=Delta(mu_i/T)-h_i^* Delta beta equals R Delta ln(p_i) up to the retained
//! interval/reduction roundoff. At equal temperatures h_i^*=h_i(T).
//!
//! With z_i=j_i/M_i and a common positive face composition x_i, the balanced
//! force is d_i=x_i(a_i-M_i lambda)/(R distance), lambda=sum(x_i a_i)/Mbar.
//! Authored heats H_i add -x_i M_i H_i Delta beta/(R distance), paired with
//! reduced heat q=-k Delta T/distance+sum(H_i j_i). Arithmetic velocity in
//! stress work gives Delta(u/T)-ubar Delta beta=betabar Delta u.
//!
//! Consequently the paired semidiscrete entropy rate per area is
//!   -R distance sum(z_i d_i/x_i)
//!   + k Delta T^2/(T_left T_right distance)
//!   + longitudinal_viscosity betabar Delta u^2/distance.
//! Original Maxwell--Stefan row/mass residuals and floating-point reductions
//! remain explicit limits of this identity. This is not a fully discrete
//! theorem; the caller must still admit the actual finite-step state.
//!
//! A species may be absent on BOTH traces. One-sided zeros are refused before
//! a source query: no logarithmic floor, artificial trace inventory or fallback
//! to the engineering face is introduced. Source sampling remains the original
//! width-weighted T/p/x, and the original material/kinetic law is evaluated once.
//!
//! Continuum force and Onsager convention: Van-Brunt, Farrell and Monroe,
//! https://arxiv.org/abs/2109.05082, Eq. 2.8 and Section 3. The discrete divided
//! caloric/force pairing above is derived here for the admitted source model.

use super::*;
use crate::multicomponent::wall::log_ratio;
use fs_math::dd::Dd;

#[cfg(test)]
mod tests;

/// Original axial inputs with entropy-conjugate interior energy and diffusion.
/// wl and wr are physical cell widths; neither trace nor caller state changes.
#[allow(clippy::too_many_arguments)]
pub(in crate::multicomponent::dissipation) fn interior<L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    workspace: &mut MaxwellStefanWorkspaceV1,
    left: Trace<'_>,
    right: Trace<'_>,
    wl: f64,
    wr: f64,
    mass_ratio: f64,
    binary: &mut [f64],
    out: &mut [f64],
    controls: MaxwellStefanSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    poll(cancelled)?;
    let n = gas.species_count();
    if out.len() != n + 2 || !positive(wl) || !positive(wr) {
        return Err(DissipationError::Invalid(
            "entropy-conjugate molecular face dimensions or distances",
        ));
    }
    let pl = left.primitive;
    let pr = right.primitive;
    let distance = 0.5 * wl + 0.5 * wr;
    let right_weight = (0.5 * wl) / distance;
    let left_weight = 1.0 - right_weight;
    let t = if pl.temperature == pr.temperature {
        pl.temperature
    } else {
        left_weight * pl.temperature + right_weight * pr.temperature
    };
    let p = if pl.pressure == pr.pressure {
        pl.pressure
    } else {
        left_weight * pl.pressure + right_weight * pr.pressure
    };
    let delta_t = pr.temperature - pl.temperature;
    let grad_t = delta_t / distance;
    let grad_u = (pr.velocity - pl.velocity) / distance;
    let delta_beta = reciprocal_difference(pl.temperature, pr.temperature)?;
    if !positive(distance)
        || !positive(t)
        || !positive(p)
        || !grad_t.is_finite()
        || !grad_u.is_finite()
    {
        return Err(DissipationError::Invalid(
            "unrepresentable entropy-conjugate molecular face state",
        ));
    }

    let mut xl = [0.0; MAX_SPECIES];
    let mut xr = [0.0; MAX_SPECIES];
    let mut x = [0.0; MAX_SPECIES];
    fractions(gas, left.densities, &mut xl[..n])?;
    fractions(gas, right.densities, &mut xr[..n])?;
    let mass_ratio = active_mass_ratio(diffusion, &xl[..n], &xr[..n], mass_ratio);
    let mut contrast =
        (t / pl.temperature.min(pr.temperature)).max(p / pl.pressure.min(pr.pressure));
    for k in 0..n {
        poll(cancelled)?;
        if (left.densities[k] == 0.0) != (right.densities[k] == 0.0) {
            return Err(DissipationError::Invalid(
                "entropy-conjugate diffusion requires common positive species support",
            ));
        }
        x[k] = left_weight * xl[k] + right_weight * xr[k];
        if left.densities[k] > 0.0 {
            if !positive(xl[k]) || !positive(xr[k]) || !positive(x[k]) {
                return Err(DissipationError::Invalid(
                    "entropy-conjugate active mole fraction is unrepresentable",
                ));
            }
            contrast = contrast.max(x[k] / xl[k].min(xr[k]));
        }
    }
    let pivot = (0..n)
        .max_by(|&a, &b| x[a].total_cmp(&x[b]))
        .expect("nonempty species bank");
    close_simplex(&mut x[..n], pivot)?;
    let masses = diffusion.molar_masses();
    let mean_mass = sum((0..n).map(|k| x[k] * masses[k].value()));
    if !positive(mean_mass) || !positive(contrast) {
        return Err(DissipationError::Invalid(
            "entropy-conjugate mixture scaling overflow",
        ));
    }

    let coefficients = law.evaluate(
        Temperature::new(t),
        Pressure::new(p),
        &x[..n],
        binary,
        &mut || cancelled(),
    )?;
    coefficients.validate()?;
    let maximum_diffusion = maximum_active_diffusivity(binary, &x[..n])?;
    let mut heats = [0.0; MAX_SPECIES];
    let coupled = law.heats_of_transport(
        Temperature::new(t),
        Pressure::new(p),
        &x[..n],
        &mut heats[..n],
        &mut || cancelled(),
    )?;
    let mut thermal_force_scale = 0.0_f64;
    let mut maximum_heat = 0.0_f64;
    if coupled {
        if !heats[..n].iter().all(|h| h.is_finite()) {
            return Err(DissipationError::Invalid(
                "non-finite source heat of transport",
            ));
        }
        let reference = heats[pivot];
        let mean_heat =
            sum((0..n).map(|k| (x[k] * (masses[k].value() / mean_mass)) * (heats[k] - reference)));
        for k in 0..n {
            heats[k] = (heats[k] - reference) - mean_heat;
            if x[k] > 0.0 {
                let scale =
                    masses[k].value() * ((heats[k] / R) / pl.temperature.min(pr.temperature));
                thermal_force_scale = thermal_force_scale.max(scale.abs());
                maximum_heat = maximum_heat.max(heats[k].abs());
            }
        }
        if !mean_heat.is_finite()
            || !heats[..n].iter().all(|h| h.is_finite())
            || !thermal_force_scale.is_finite()
        {
            return Err(DissipationError::Invalid(
                "heat-of-transport scaling overflow",
            ));
        }
    }

    let mut enthalpies = [Dd::ZERO; MAX_SPECIES];
    let mut efforts = [0.0; MAX_SPECIES];
    let log_t = log_ratio(pr.temperature, pl.temperature);
    for (k, model) in gas.source().models().iter().enumerate() {
        poll(cancelled)?;
        if x[k] == 0.0 {
            continue;
        }
        let sampled = Temperature::new(t);
        let initial = Temperature::new(pl.temperature);
        let final_t = Temperature::new(pr.temperature);
        let h_sample = Dd::from_f64(
            model
                .molar_internal_energy_change(
                    Temperature::new(gas.domain().minimum_temperature),
                    sampled,
                )
                .map_err(FlowError::from)?
                .value(),
        ) + Dd::from_f64(R) * Dd::from_f64(t);
        let relative = if pl.temperature == pr.temperature {
            Dd::ZERO
        } else {
            Dd::from_f64(
                model
                    .molar_enthalpy_change(sampled, final_t)
                    .map_err(FlowError::from)?
                    .value(),
            ) / Dd::from_f64(pr.temperature)
                - Dd::from_f64(
                    model
                        .molar_enthalpy_change(sampled, initial)
                        .map_err(FlowError::from)?
                        .value(),
                ) / Dd::from_f64(pl.temperature)
                - Dd::from_f64(
                    model
                        .molar_entropy_change(initial, final_t)
                        .map_err(FlowError::from)?
                        .value(),
                )
        };
        let h = if delta_beta == 0.0 {
            h_sample
        } else {
            h_sample + relative / Dd::from_f64(delta_beta)
        };
        checked(h)?;
        enthalpies[k] = h;
        // Retain the actual finite-precision caloric remainder, rather than
        // replacing it by zero on the strength of the exact-arithmetic identity.
        let remainder = relative - (h - h_sample) * Dd::from_f64(delta_beta);
        efforts[k] = checked(
            remainder
                + Dd::from_f64(R)
                    * (Dd::from_f64(log_ratio(right.densities[k], left.densities[k]))
                        + Dd::from_f64(log_t)),
        )?;
    }

    // Center the molar efforts on one species BEFORE reduction. This is the
    // explicit barycentric pressure projection, not a post-solve flux repair.
    // The pivot effort is exactly zero even for a one-species gas.
    let reference = efforts[pivot];
    for k in 0..n {
        efforts[k] = if k == pivot || x[k] == 0.0 {
            0.0
        } else {
            checked(
                Dd::from_f64(efforts[k])
                    - Dd::from_f64(masses[k].value()) / Dd::from_f64(masses[pivot].value())
                        * Dd::from_f64(reference),
            )?
        };
    }
    let mut mean_effort = Dd::ZERO;
    for k in 0..n {
        mean_effort = mean_effort + Dd::from_f64(x[k]) * Dd::from_f64(efforts[k]);
    }
    mean_effort = mean_effort / Dd::from_f64(mean_mass);
    let mut forces = [0.0; MAX_SPECIES];
    for k in 0..n {
        if x[k] == 0.0 {
            continue;
        }
        let mut effort = Dd::from_f64(efforts[k]) - Dd::from_f64(masses[k].value()) * mean_effort;
        if coupled {
            effort = effort
                - Dd::from_f64(masses[k].value())
                    * Dd::from_f64(heats[k])
                    * Dd::from_f64(delta_beta);
        }
        forces[k] =
            checked(Dd::from_f64(x[k]) * effort / Dd::from_f64(R) / Dd::from_f64(distance))?;
    }
    let evaluated = diffusion.evaluate_driving_forces(
        &x[..n],
        &forces[..n],
        AmountConcentration::new((p / R) / t),
        binary,
        controls,
        workspace,
        &mut *cancelled,
    )?;
    out[..n].copy_from_slice(evaluated.mass_fluxes());
    let mut report = TransportFaceReport {
        diffusivity_m2_s: 0.0,
        dynamic_viscosity_pa_s: coefficients.dynamic_viscosity.value(),
        thermal_conductivity_w_m_k: coefficients.thermal_conductivity.value(),
        maximum_gradient_residual_per_m: evaluated.maximum_gradient_residual_per_m(),
        mass_frame_residual_kg_m2_s: evaluated.mass_frame_residual_kg_m2_s(),
        // The material law samples this temperature; species energy uses the
        // divided source enthalpies, which need not share one effective T.
        enthalpy_temperature_k: t,
        reduced_heat_flux_w_m2: 0.0,
    };
    let longitudinal = coefficients.longitudinal_viscosity()?;
    let stress = longitudinal * grad_u;
    out[n] = -stress;
    let mut heat =
        Dd::ZERO - Dd::from_f64(coefficients.thermal_conductivity.value()) * Dd::from_f64(grad_t);
    let mut energy = Dd::ZERO;
    for k in 0..n {
        if out[k] != 0.0 {
            energy =
                energy + enthalpies[k] * Dd::from_f64(out[k]) / Dd::from_f64(masses[k].value());
            if coupled {
                heat = heat + Dd::from_f64(heats[k]) * Dd::from_f64(out[k]);
            }
        }
    }
    report.reduced_heat_flux_w_m2 = checked(heat)?;
    let work_velocity = 0.5 * pl.velocity + 0.5 * pr.velocity;
    out[n + 1] = checked(energy + heat - Dd::from_f64(work_velocity) * Dd::from_f64(stress))?;
    if !out.iter().all(|v| v.is_finite()) {
        return Err(DissipationError::Invalid(
            "non-finite entropy-conjugate molecular flux",
        ));
    }

    // Logarithmic chemical forces require trace-composition/temperature
    // contrast in the existing engineering time-step scale. This is not an
    // eigenvalue, positivity or invariant-domain bound for the nonlinear flux.
    let gamma = 1.0
        + (pl.pressure / pl.temperature / left.cv_density)
            .max(pr.pressure / pr.temperature / right.cv_density);
    let molecular_scale = maximum_diffusion * mass_ratio * gamma * contrast;
    let thermal_scale =
        coefficients.thermal_conductivity.value() / left.cv_density.min(right.cv_density);
    report.diffusivity_m2_s = molecular_scale
        .max(longitudinal / pl.density.min(pr.density))
        .max(thermal_scale);
    if thermal_force_scale > 0.0 {
        let storage_scale = (maximum_heat / pl.temperature.min(pr.temperature))
            * (pl.density.max(pr.density) / left.cv_density.min(right.cv_density));
        report.diffusivity_m2_s = report.diffusivity_m2_s.max(
            molecular_scale * (1.0 + thermal_force_scale) * (1.0 + storage_scale) + thermal_scale,
        );
    }
    if !report.diffusivity_m2_s.is_finite() || report.diffusivity_m2_s < 0.0 {
        return Err(DissipationError::Invalid(
            "entropy-conjugate face scale overflow",
        ));
    }
    poll(cancelled)?;
    Ok(report)
}

/// Conjugate cell-to-bath layer over the original physical half-cell distance.
///
/// The bath is an infinite stationary source of energy and species. Use its
/// actual entropy effort on the SAME conservative flux: E_ref/T_b minus
/// sum(mu_i,ref/T_b * j_i/M_i). Relative to bath partial enthalpies its heat is
/// E_ref-sum(h_i,ref(T_b)*j_i/M_i), including the enthalpy-layer correction and
/// bath share of viscous conversion. The returned reduced_heat_flux retains
/// the interior Fourier/Dufour convention; it is NOT that bath heat.
///
/// Symmetric source sampling and arithmetic stress work give the original
/// two-trace entropy pairing while neither bath state nor source is advanced.
/// Reusing the kernel also preserves common-support refusal before source
/// evaluation and avoids another constitutive or caloric query for accounting.
#[allow(clippy::too_many_arguments)]
pub(in crate::multicomponent::dissipation) fn reservoir<L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    diffusion: &MaxwellStefanV1,
    workspace: &mut MaxwellStefanWorkspaceV1,
    left: Trace<'_>,
    right: Trace<'_>,
    distance: f64,
    boundary_on_left: bool,
    mass_ratio: f64,
    binary: &mut [f64],
    out: &mut [f64],
    controls: MaxwellStefanSolveV1,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    let bath = if boundary_on_left { &left } else { &right };
    if bath.primitive.velocity != 0.0 {
        return Err(DissipationError::Invalid(
            "entropy-conjugate molecular reservoir is not stationary",
        ));
    }
    interior(
        gas, law, diffusion, workspace, left, right, distance, distance,
        mass_ratio, binary, out, controls, cancelled,
    )
}

fn reciprocal_difference(left: f64, right: f64) -> Result<f64, DissipationError> {
    if left == right {
        return Ok(0.0);
    }
    let (lo, hi, sign) = if right > left {
        (left, right, -1.0)
    } else {
        (right, left, 1.0)
    };
    let difference = ((hi - lo) / hi) / lo;
    if !positive(difference) {
        return Err(DissipationError::Invalid(
            "unrepresentable reciprocal temperature difference",
        ));
    }
    Ok(sign * difference)
}

fn checked(value: Dd) -> Result<f64, DissipationError> {
    let result = value.to_f64();
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(DissipationError::Invalid(
            "non-finite entropy-conjugate source effort",
        ))
    }
}
