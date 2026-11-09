//! Actual chemical-potential/energy pairing on the existing mass-frame solver.
//! The PR residual face enthalpy is a symmetric endpoint approximation, not an
//! exact finite-temperature Gibbs--Helmholtz integral. The same chosen h* enters
//! both species effort and donor-to-face energy correction, so their finite
//! entropy pairing remains exact independently of that spatial approximation.
use super::super::force;
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanErrorV1 as Error, MaxwellStefanSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};
use fs_thermochem::{Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

pub(super) struct Evaluation<'a> {
    pub kinetic: force::Evaluation<'a>,
    pub centered_heats: Vec<f64>,
    pub face_residual_enthalpies: Vec<f64>,
    pub fourier: f64,
    pub reciprocal: f64,
    pub correction: f64,
    pub additional: f64,
    pub entropy_from_sources: f64,
    pub entropy_from_dissipation: f64,
}

fn finite(v: Dd) -> Result<f64, Error> {
    let x = v.to_f64();
    if [v.hi, v.lo, x].iter().all(|x| x.is_finite()) {
        Ok(x)
    } else {
        Err(Error::Arithmetic(
            "nonisothermal real-gas reduction overflowed",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

/// All arrays are on the original MS/NASA species axis. log_fugacity_change is
/// Delta ln(x_i p phi_i) at ACTUAL endpoint T,p,composition; residual_enthalpy
/// comes from original fixed-T,p partial PR calorics at those actual densities.
/// No source data, ideal density, perturbed root or coefficient is invented.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate<'a>(
    model: &MaxwellStefanV1,
    calorics: &[Nasa9StandardStateModelV1],
    x: &[f64],
    temperatures: [f64; 2],
    log_fugacity_change: &[f64],
    residual_enthalpy: [&[f64]; 2],
    concentration: f64,
    binary: &[f64],
    heats: &[f64],
    conductivity: f64,
    distance: f64,
    controls: MaxwellStefanSolveV1,
    entropy_tolerance: f64,
    work: &'a mut MaxwellStefanWorkspaceV1,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Evaluation<'a>, Error> {
    poll(&mut cancelled)?;
    let n = model.species().len();
    if n == 0
        || n > 32
        || x.len() != n
        || calorics.len() != n
        || heats.len() != n
        || log_fugacity_change.len() != n
        || residual_enthalpy.iter().any(|h| h.len() != n)
        || !temperatures.iter().all(|t| t.is_finite() && *t > 0.0)
        || !distance.is_finite()
        || distance <= 0.0
        || !conductivity.is_finite()
        || conductivity < 0.0
        || heats
            .iter()
            .chain(log_fugacity_change)
            .chain(residual_enthalpy[0])
            .chain(residual_enthalpy[1])
            .any(|v| !v.is_finite())
        || calorics.iter().enumerate().any(|(i, s)| {
            s.species() != &model.species()[i] || s.molar_mass() != model.molar_masses()[i]
        })
    {
        return Err(Error::Invalid(
            "bounded original thermal species and physical data required",
        ));
    }
    let [tl, tr] = temperatures;
    let t = 0.5 * tl + 0.5 * tr;
    let delta_t = tr - tl;
    let beta = if tl == tr {
        0.0
    } else {
        ((tl - tr) / tl.max(tr)) / tl.min(tr)
    };
    if !beta.is_finite() || (tl != tr && beta == 0.0) {
        return Err(Error::Arithmetic(
            "reciprocal temperature difference is unrepresentable",
        ));
    }
    let mut masses = [0.0; 32];
    for (i, mass) in model.molar_masses().iter().enumerate() {
        masses[i] = mass.value();
    }
    let pivot = (0..n)
        .max_by(|&i, &j| masses[i].total_cmp(&masses[j]))
        .unwrap();
    let mut mean_mass = Dd::ZERO;
    let mut mean_heat = Dd::ZERO;
    for i in 0..n {
        mean_mass = mean_mass + Dd::from_f64(x[i]) * Dd::from_f64(masses[i]);
        mean_heat = mean_heat
            + Dd::from_f64(x[i])
                * Dd::from_f64(masses[i])
                * (Dd::from_f64(heats[i]) - Dd::from_f64(heats[pivot]));
    }
    if finite(mean_mass)? <= 0.0 {
        return Err(Error::Invalid(
            "positive original thermal mean mass required",
        ));
    }
    mean_heat = mean_heat / mean_mass;
    let mut centered_heats = Vec::with_capacity(n);
    let mut face_residual_enthalpies = Vec::with_capacity(n);
    let mut residual_faces = [Dd::ZERO; 32];
    let mut relative = [Dd::ZERO; 32];
    let mut divided = [Dd::ZERO; 32];
    let mut gradients = [0.0; 32];
    for i in 0..n {
        poll(&mut cancelled)?;
        let h = finite(Dd::from_f64(heats[i]) - Dd::from_f64(heats[pivot]) - mean_heat)?;
        centered_heats.push(h);
        let hr = Dd::from_f64(0.5)
            * (Dd::from_f64(residual_enthalpy[0][i]) + Dd::from_f64(residual_enthalpy[1][i]));
        face_residual_enthalpies.push(finite(hr)?);
        residual_faces[i] = hr;
        let ideal = if beta == 0.0 {
            Dd::ZERO
        } else {
            Dd::from_f64(
                calorics[i]
                    .molar_enthalpy_change(Temperature::new(t), Temperature::new(tr))
                    .map_err(|_| Error::Invalid("original thermal enthalpy increment refused"))?
                    .value(),
            ) / Dd::from_f64(tr)
                - Dd::from_f64(
                    calorics[i]
                        .molar_enthalpy_change(Temperature::new(t), Temperature::new(tl))
                        .map_err(|_| Error::Invalid("original thermal enthalpy increment refused"))?
                        .value(),
                ) / Dd::from_f64(tl)
                - Dd::from_f64(
                    calorics[i]
                        .molar_entropy_change(Temperature::new(tl), Temperature::new(tr))
                        .map_err(|_| Error::Invalid("original thermal entropy increment refused"))?
                        .value(),
                )
        };
        let correction = if beta == 0.0 {
            Dd::ZERO
        } else {
            ideal / Dd::from_f64(beta)
        };
        finite(correction)?;
        relative[i] = ideal;
        divided[i] = correction;
        let effort = ideal - correction * Dd::from_f64(beta)
            + Dd::from_f64(R) * Dd::from_f64(log_fugacity_change[i])
            - (hr + Dd::from_f64(masses[i]) * Dd::from_f64(h)) * Dd::from_f64(beta);
        gradients[i] = finite(effort / Dd::from_f64(R) / Dd::from_f64(distance))?;
        if effort.to_f64() != 0.0 && gradients[i] == 0.0 {
            return Err(Error::Arithmetic("real-gas thermal effort underflowed"));
        }
    }
    let kinetic = force::evaluate(
        model,
        x,
        &gradients[..n],
        concentration,
        binary,
        controls,
        entropy_tolerance,
        work,
        &mut cancelled,
    )?;
    let fourier =
        Dd::ZERO - Dd::from_f64(conductivity) * Dd::from_f64(delta_t) / Dd::from_f64(distance);
    let mut reciprocal = Dd::ZERO;
    let mut correction = Dd::ZERO;
    let mut donor_gradients = [0.0; 32];
    for i in 0..n {
        poll(&mut cancelled)?;
        let molar = kinetic.diffusion.molar_fluxes()[i];
        let mass = kinetic.diffusion.mass_fluxes()[i];
        let side = usize::from(mass < 0.0);
        let ideal_donor_to_mid = Dd::from_f64(
            calorics[i]
                .molar_enthalpy_change(Temperature::new(temperatures[side]), Temperature::new(t))
                .map_err(|_| Error::Invalid("original donor caloric increment refused"))?
                .value(),
        );
        let donor_to_face = ideal_donor_to_mid + divided[i] + residual_faces[i]
            - Dd::from_f64(residual_enthalpy[side][i]);
        correction = correction + Dd::from_f64(molar) * donor_to_face;
        reciprocal = reciprocal + Dd::from_f64(mass) * Dd::from_f64(centered_heats[i]);
        // Actual Delta(mu/T) minus actual donor h_i Delta(beta), with the
        // arbitrary NASA h_i(Tmid) datum cancelled BEFORE floating arithmetic.
        donor_gradients[i] = finite(
            (relative[i]
                + Dd::from_f64(R) * Dd::from_f64(log_fugacity_change[i])
                + (ideal_donor_to_mid - Dd::from_f64(residual_enthalpy[side][i]))
                    * Dd::from_f64(beta))
                / Dd::from_f64(R)
                / Dd::from_f64(distance),
        )?;
    }
    let mut ignored_forces = [0.0; 32];
    let donor_reduced = force::project_with_gauge(
        x,
        &masses[..n],
        &donor_gradients[..n],
        &mut ignored_forces[..n],
    )
    .map_err(Error::Invalid)?;
    // Round the COMPLETE physical additional flux once, then audit exactly
    // that exposed flux. Readding rounded diagnostic components could lose a
    // small residual when Fourier, reciprocal and caloric terms oppose.
    let additional = finite(fourier + reciprocal + correction)?;
    let mut source_power = Dd::from_f64(additional) * Dd::from_f64(beta);
    for (i, &molar) in kinetic.diffusion.molar_fluxes().iter().enumerate() {
        poll(&mut cancelled)?;
        source_power = source_power
            - Dd::from_f64(molar) * donor_reduced[i] * Dd::from_f64(R) * Dd::from_f64(distance);
    }
    let dissipative =
        Dd::from_f64(kinetic.drag_entropy) * Dd::from_f64(distance) + fourier * Dd::from_f64(beta);
    let source_power = finite(source_power)?;
    let dissipative = finite(dissipative)?;
    let allowance = entropy_tolerance * distance
        + controls.relative_tolerance * source_power.abs().max(dissipative.abs());
    if !allowance.is_finite()
        || dissipative < 0.0
        || source_power < -allowance
        || (source_power - dissipative).abs() > allowance
    {
        return Err(Error::Arithmetic(
            "actual PR thermal chemical/energy power did not close dissipation",
        ));
    }
    poll(&mut cancelled)?;
    Ok(Evaluation {
        kinetic,
        centered_heats,
        face_residual_enthalpies,
        fourier: finite(fourier)?,
        reciprocal: finite(reciprocal)?,
        correction: finite(correction)?,
        additional,
        entropy_from_sources: source_power,
        entropy_from_dissipation: dissipative,
    })
}

#[cfg(test)]
#[path = "constitutive/tests.rs"]
mod tests;
