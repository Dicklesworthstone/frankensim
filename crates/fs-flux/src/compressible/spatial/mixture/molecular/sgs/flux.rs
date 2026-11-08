//! Shared conservative modeled species, heat and stress at an admitted face.
//!
//! Interior and complete-reservoir faces share the same common-species
//! mass-frame closure and original caloric references. Constitutive Cp uses the
//! midpoint temperature; the source receipt selects the partial-enthalpy
//! temperature, and the caller supplies the physical mechanical work velocity.

use super::*;

pub(super) struct Face<'a> {
    pub(super) partial_densities: [&'a [f64]; 2],
    pub(super) primitive: [Primitive; 2],
    pub(super) cv_density: [f64; 2],
    /// Component-first physical tensor du_i/dx_j.
    pub(super) gradient: [[f64; 3]; 3],
    pub(super) axis: usize,
    pub(super) distance: f64,
    pub(super) work_velocity: [f64; 3],
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add(
    gas: &MixtureGas,
    face: Face<'_>,
    mu: f64,
    dynamic_heat: Option<(f64, f64)>,
    dynamic_species: Option<f64>,
    prandtl: Option<f64>,
    schmidt: Option<f64>,
    mut receipt: TransportFaceReport,
    controls: MaxwellStefanSolveV1,
    out: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<TransportFaceReport, DissipationError> {
    let Face {
        partial_densities,
        primitive,
        cv_density,
        gradient,
        axis,
        distance,
        work_velocity,
    } = face;
    let n = gas.species_count();
    let [pl, pr] = primitive;
    let temperature = receipt.enthalpy_temperature_k;
    let constitutive_temperature = 0.5 * pl.temperature + 0.5 * pr.temperature;
    let mut y = [0.0; MAX_SPECIES];
    let mut species_flux = [0.0; MAX_SPECIES];
    let mass_diffusion =
        dynamic_species.unwrap_or_else(|| mu / schmidt.expect("fixed-Schmidt species closure"));
    for k in 0..n {
        let yl = partial_densities[0][k] / pl.density;
        let yr = partial_densities[1][k] / pr.density;
        y[k] = 0.5 * yl + 0.5 * yr;
        species_flux[k] = -mass_diffusion * ((yr - yl) / distance);
    }
    let pivot = (0..n)
        .max_by(|&a, &b| y[a].total_cmp(&y[b]))
        .expect("admitted species bank");
    // Close a present species, so exactly absent species stay absent. With
    // common Schmidt this changes only summation roundoff, not the model.
    y[pivot] = 1.0 - sum((0..n).filter(|&k| k != pivot).map(|k| y[k]));
    species_flux[pivot] = -sum((0..n).filter(|&k| k != pivot).map(|k| species_flux[k]));
    let mut cp = [0.0; MAX_SPECIES];
    let mut enthalpy = [0.0; MAX_SPECIES];
    for (k, model) in gas.source().models().iter().enumerate() {
        if k % 16 == 0 {
            poll(cancelled)?;
        }
        if y[k] == 0.0 {
            continue;
        }
        let inverse_mass = model.molar_mass().value().recip();
        cp[k] = y[k]
            * inverse_mass
            * model
                .evaluate_properties(Temperature::new(constitutive_temperature))
                .map_err(FlowError::from)?
                .cp()
                .value();
        // Evaluate changes within each original caloric source. Subtracting
        // two large formation enthalpies would destroy gauge invariance.
        if mass_diffusion != 0.0 {
            enthalpy[k] = inverse_mass
                * (model
                    .molar_internal_energy_change(
                        Temperature::new(gas.domain().minimum_temperature),
                        Temperature::new(temperature),
                    )
                    .map_err(FlowError::from)?
                    .value()
                    + R * temperature);
        }
    }
    let cp = sum(cp[..n].iter().copied());
    let conductivity = if let Some((weight, maximum_ratio)) = dynamic_heat {
        let conductivity = weight * cp;
        if conductivity != 0.0 {
            // The original face reports the same source evaluation that
            // supplied its Fourier law; stateful laws are not queried twice.
            let ratio = conductivity / receipt.thermal_conductivity_w_m_k;
            if !ratio.is_finite() || ratio > maximum_ratio {
                return Err(DissipationError::Invalid(
                    "mixture subgrid heat/molecular conductivity ratio exceeded",
                ));
            }
        }
        conductivity
    } else {
        mu * (cp / prandtl.expect("fixed-Prandtl heat closure"))
    };
    let heat = -conductivity * ((pr.temperature - pl.temperature) / distance);
    let nonmechanical_energy = heat + sum((0..n).map(|k| enthalpy[k] * species_flux[k]));
    let mut energy = nonmechanical_energy;
    if mass_diffusion != 0.0 {
        for k in 0..n {
            out[k] += species_flux[k];
        }
    }
    if mu != 0.0 {
        let [a, b] = tangents(axis);
        let divergence = gradient[axis][axis] + gradient[a][a] + gradient[b][b];
        for (frame, component) in [axis, a, b].into_iter().enumerate() {
            let normal = gradient[component][axis];
            let transpose = if component == axis {
                normal - (2.0 / 3.0) * divergence
            } else {
                gradient[axis][component]
            };
            let stress = mu * (normal + transpose);
            out[n + frame] -= stress;
            energy -= work_velocity[component] * stress;
        }
    }
    out[n + 3] += energy;
    receipt.reduced_heat_flux_w_m2 += heat;
    receipt.nonmechanical_energy_flux_w_m2 += nonmechanical_energy;
    receipt.mass_frame_residual_kg_m2_s = sum(out[..n].iter().copied());
    let frame_scale = sum(out[..n].iter().map(|j| j.abs()));
    if receipt.mass_frame_residual_kg_m2_s.abs()
        > controls.absolute_mass_flux_tolerance_kg_m2_s + controls.relative_tolerance * frame_scale
    {
        return Err(DissipationError::Invalid(
            "mixture subgrid total mass-frame residual",
        ));
    }
    let cv_min = cv_density[0].min(cv_density[1]);
    let rho_min = pl.density.min(pr.density);
    let gamma = 1.0
        + (pl.pressure / (pl.temperature * cv_density[0]))
            .max(pr.pressure / (pr.temperature * cv_density[1]));
    let enthalpy_span = (0..n)
        .filter(|&k| y[k] > 0.0)
        .map(|k| (enthalpy[k] - enthalpy[pivot]).abs())
        .fold(0.0, f64::max);
    let species_scale = (mass_diffusion / rho_min)
        * gamma
        * (1.0 + (enthalpy_span / temperature) * (pl.density.max(pr.density) / cv_min));
    // Sum molecular and modeled bounds: max(D_m,D_t) alone underestimates
    // their combined diffusion. Include enthalpy/species coupling and the
    // full 3-D stress row sum, then let the owner apply all inverse lengths.
    let scalar_scale = conductivity / cv_min + species_scale;
    receipt.diffusivity_m2_s += ((20.0 / 3.0) * mu / rho_min).max(scalar_scale);
    receipt.nonmechanical_diffusivity_m2_s += scalar_scale;
    if !positive(cp)
        || !conductivity.is_finite()
        || !mass_diffusion.is_finite()
        || !out.iter().all(|v| v.is_finite())
        || !receipt.reduced_heat_flux_w_m2.is_finite()
        || !receipt.nonmechanical_energy_flux_w_m2.is_finite()
        || !receipt.nonmechanical_diffusivity_m2_s.is_finite()
        || !receipt.diffusivity_m2_s.is_finite()
        || !receipt.mass_frame_residual_kg_m2_s.is_finite()
        || !frame_scale.is_finite()
        || !gamma.is_finite()
        || !enthalpy_span.is_finite()
        || !species_scale.is_finite()
    {
        return Err(DissipationError::Invalid(
            "unrepresentable mixture subgrid face transport",
        ));
    }
    poll(cancelled)?;
    Ok(receipt)
}
