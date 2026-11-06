//! Photon-conserving ionization of a rigid, uniformly mixed gas column.
//!
//! A collimated monochromatic beam crosses area A and volume V. Its absorbed
//! photon rate is A*F*(1-exp(-tau)), where tau=sigma*NA*n_A/A. Mixing maintains
//! one neutral density through the column as it ionizes. For fixed incident
//! number flux F and cross section sigma, the optical-depth equation has the
//! exact solution expm1(tau(t))=expm1(tau(0))*exp(-sigma*F*t).
//! This removes a chemical time-step restriction and accounts for the changing
//! opacity without consuming more photons or neutrals than are available.
//!
//! Each absorbed photon creates one ion and one electron. Its energy enters
//! original-source chemical storage once; the excess enters the electron
//! sensible pool, whose heat capacity changes with the electron population.
//! Exact zero initial electrons are supported. Heavy temperature is unchanged
//! in the stationary-ion approximation. Source reference offsets cancel before
//! energy subtraction, including arbitrary conserved elemental/charge offsets.
//!
//! This is an exact update of a declared zero-dimensional mixing model, not
//! a resolved ionization front. The radiation is quasi-static, photoelectrons
//! thermalize immediately, and recoil, scattering, secondary ionization,
//! excited states and radiative recombination are absent. Other reaction and
//! thermal owners can consume the returned common gas state. No radiation
//! entropy model is supplied, so no closed-parcel entropy inequality is claimed.

use fs_math::det::{exp, expm1, ln};
use fs_qty::{Amount, Area, Energy, Temperature, Time};
use fs_thermochem::electron_collision::AVOGADRO_PER_MOL as NA;
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationErrorV1, PhotonNumberFluxV1,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

/// Constant illumination of the whole declared gas column.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoionizationExposure {
    /// Photon energy and incident number flux, not an energy flux or density.
    pub beam: PhotonNumberFluxV1,
    /// Area normal to the beam [m²]; column length is the gas volume divided by A.
    pub illuminated_area: Area,
    /// Physical exposure duration [s].
    pub duration: Time,
}

/// Numerical publication allowances; these are not measured model errors.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotoionizationControls {
    /// Charge residual scaled by actual carriers, strictly in (0,1).
    pub relative_charge_tolerance: f64,
    /// Absolute electron sensible plus chemical minus absorbed energy allowance [J].
    pub energy_tolerance_j: f64,
}

/// Original-source, domain, arithmetic, or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum PhotoionizationError {
    /// Invalid exposure, control, or mismatched original species bank.
    Invalid(&'static str),
    /// Bound-free cross-section evaluation refused.
    Source(PhotoionizationErrorV1),
    /// A finite positive physical increment cannot be represented.
    Arithmetic(&'static str),
    /// Published source energy does not close within the requested allowance.
    EnergyBalance { residual_j: f64 },
    /// Common charged-state or original species temperature admission refused.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation before publishing a complete endpoint.
    Cancelled,
}
impl core::fmt::Display for PhotoionizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "monochromatic gas photoionization: {self:?}")
    }
}
impl std::error::Error for PhotoionizationError {}

/// Photon, population and energy accounts for the same physical events.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PhotoionizationBalance {
    /// Incident photons crossing the declared area over the interval.
    pub incident_photons: f64,
    /// One absorbed photon for each new ion/electron pair.
    pub absorbed_photons: f64,
    /// Incident photons that leave the column.
    pub transmitted_photons: f64,
    /// Original-source ionization formation energy stored [J].
    pub chemical_storage_change_j: f64,
    /// Actual electron sensible-energy increment, including population change [J].
    pub electron_sensible_energy_change_j: f64,
    /// Absorbed photon energy [J], before chemical/sensible partition.
    pub absorbed_energy_j: f64,
    /// Electron sensible plus chemical storage minus absorbed energy [J].
    pub energy_residual_j: f64,
}

/// Source-retained endpoint ready for existing collisional, conductivity and optics owners.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotoionizationResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicPhotoionizationV1,
    exposure: PhotoionizationExposure,
    controls: PhotoionizationControls,
    extent: Amount,
    optical_depths: [f64; 2],
    balance: PhotoionizationBalance,
}
impl PhotoionizationResult {
    /// Accepted common gas state; no equilibrium projection was performed.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 { &self.state }
    /// Original immutable initial state.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 { &self.initial }
    /// Original cross section, species bank and level-energy source.
    pub const fn source(&self) -> &AtomicPhotoionizationV1 { &self.source }
    /// Applied constant illumination.
    pub const fn exposure(&self) -> PhotoionizationExposure { self.exposure }
    /// Numerical publication allowances.
    pub const fn controls(&self) -> PhotoionizationControls { self.controls }
    /// Nonnegative number of ionization events in moles.
    pub const fn reaction_extent(&self) -> Amount { self.extent }
    /// Initial and final neutral bound-free optical depths at the beam energy.
    pub const fn optical_depths(&self) -> [f64; 2] { self.optical_depths }
    /// Photon and energy accounts.
    pub const fn balance(&self) -> PhotoionizationBalance { self.balance }
    /// Energy supplied by absorbed photons [J].
    pub const fn absorbed_energy(&self) -> Energy { Energy::new(self.balance.absorbed_energy_j) }
}

/// Evolve the monochromatic, perfectly mixed column analytically.
///
/// All incoming photons are accounted for even with zero absorption. A zero
/// flux/duration is an identity, but invalid source energy or bank binding still
/// refuses. Unsupported final electron temperatures refuse transactionally;
/// no temperature clipping, artificial electron seed or missing opacity is used.
pub fn advance_monochromatic_photoionization(
    initial: &FrozenIonizedGasStateV1,
    source: &AtomicPhotoionizationV1,
    exposure: PhotoionizationExposure,
    controls: PhotoionizationControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PhotoionizationResult, PhotoionizationError> {
    poll(&mut cancelled)?;
    let channel = source.channel();
    let area = exposure.illuminated_area.value();
    let time = exposure.duration.value();
    let flux = exposure.beam.number_flux_per_m2_s;
    let photon_energy = exposure.beam.photon_energy.value();
    if !positive(area) || !nonnegative(time) || !nonnegative(flux)
        || !positive(photon_energy) || !positive(controls.energy_tolerance_j)
        || !positive(controls.relative_charge_tolerance)
        || controls.relative_charge_tolerance >= 1.0
    {
        return Err(PhotoionizationError::Invalid("finite exposure and numerical controls"));
    }
    if initial.thermodynamics() != channel.thermodynamics()
        || initial.electron_index() != Some(channel.electron_index())
        || initial.relative_charge_residual() > controls.relative_charge_tolerance
    {
        return Err(PhotoionizationError::Invalid("state and original atomic source disagree"));
    }
    let sigma = source.cross_section(exposure.beam.photon_energy)
        .map_err(PhotoionizationError::Source)?.value();
    let [a, i, e] = [channel.neutral_index(), channel.ion_index(), channel.electron_index()];
    let neutral = initial.amounts()[a].value();
    let electrons = initial.amounts()[e].value();
    let fluence = flux * time;
    let photons = fluence * area;
    let tau = ((sigma * NA) * neutral) / area;
    let dose = sigma * fluence;
    if [fluence, photons, tau, dose].iter().any(|x| !nonnegative(*x)) {
        return Err(PhotoionizationError::Arithmetic("photon exposure or optical depth"));
    }
    let (extent, remaining, tau_final) = if sigma == 0.0 || neutral == 0.0
        || flux == 0.0 || time == 0.0
    {
        (0.0, neutral, tau)
    } else {
        if tau == 0.0 || dose == 0.0 {
            return Err(PhotoionizationError::Arithmetic("positive optical depth or photon dose underflow"));
        }
        let (fraction, remaining_fraction, final_depth) = depletion(tau, dose)?;
        let lost = neutral * fraction;
        let retained = neutral * remaining_fraction;
        if lost == 0.0 || !nonnegative(lost) || !nonnegative(retained) {
            return Err(PhotoionizationError::Arithmetic("positive ionization extent"));
        }
        (lost, retained, final_depth)
    };
    let absorbed = extent * NA;
    // The exact solution obeys both inventory and photon bounds. Only tolerate
    // last-bit account roundoff; never rescale the species to repair a violation.
    let photon_allowance = 64.0 * f64::EPSILON * photons.max(absorbed);
    if !nonnegative(absorbed) || absorbed > photons + photon_allowance {
        return Err(PhotoionizationError::Arithmetic("absorbed photon budget"));
    }
    let transmitted = if absorbed > photons { 0.0 } else { photons - absorbed };
    let mut amounts = initial.amounts().to_vec();
    amounts[a] = Amount::new(remaining);
    amounts[i] = Amount::new(initial.amounts()[i].value() + extent);
    amounts[e] = Amount::new(electrons + extent);
    let chi = channel.ionization_energy_j_mol();
    let absorbed_energy = absorbed * photon_energy;
    let chemical = extent * chi;
    let cv = 1.5 * R;
    let initial_te = initial.electron_temperature().value();
    let excess = (photon_energy - chi / NA) * NA;
    let next_te = if extent == 0.0 {
        initial_te
    } else {
        if !nonnegative(excess) {
            return Err(PhotoionizationError::Arithmetic("photon energy below original ionization gap"));
        }
        let born_temperature = excess / cv;
        let total = electrons + extent;
        if extent <= electrons {
            initial_te + (extent / total) * (born_temperature - initial_te)
        } else {
            (electrons / total) * initial_te + (extent / total) * born_temperature
        }
    };
    let sensible = cv * (electrons * (next_te - initial_te) + extent * next_te);
    // Formation reference offsets cancel before this physical energy difference.
    let residual = sensible + (chemical - absorbed_energy);
    if [absorbed_energy, chemical, sensible, residual].iter().any(|x| !x.is_finite()) {
        return Err(PhotoionizationError::Arithmetic("photoelectron energy account"));
    }
    if residual.abs() > controls.energy_tolerance_j {
        return Err(PhotoionizationError::EnergyBalance { residual_j: residual });
    }
    poll(&mut cancelled)?;
    let state = FrozenIonizedGasStateV1::new(
        channel.thermodynamics(), &amounts, initial.volume(),
        initial.heavy_temperature(), Temperature::new(next_te),
        controls.relative_charge_tolerance, &mut cancelled,
    ).map_err(|error| match error {
        FrozenIonizedGasStateErrorV1::Cancelled => PhotoionizationError::Cancelled,
        error => PhotoionizationError::State(error),
    })?;
    poll(&mut cancelled)?;
    Ok(PhotoionizationResult {
        initial: initial.clone(), state, source: source.clone(), exposure, controls,
        extent: Amount::new(extent), optical_depths: [tau, tau_final],
        balance: PhotoionizationBalance {
            incident_photons: photons, absorbed_photons: absorbed,
            transmitted_photons: transmitted, chemical_storage_change_j: chemical,
            electron_sensible_energy_change_j: sensible,
            absorbed_energy_j: absorbed_energy, energy_residual_j: residual,
        },
    })
}

// ln(1+x) for this module's bounded arguments, preserving the rounding of 1+x.
// The correction is essential when thin optical depths would otherwise vanish.
fn log1p(x: f64) -> f64 {
    let y = 1.0 + x;
    if y == 1.0 { x } else { ln(y) * (x / (y - 1.0)) }
}

// Return absorbed fraction, surviving fraction and surviving optical depth.
// Calculate the small side directly instead of subtracting it from one.
fn depletion(tau: f64, dose: f64) -> Result<(f64, f64, f64), PhotoionizationError> {
    let absorbed_opacity = -expm1(-tau);
    let exposed = -expm1(-dose);
    let product = absorbed_opacity * exposed;
    let candidate = if product < 0.5 {
        // Divide before multiplying: tau*dose may underflow even when the
        // absorbed fraction (approximately dose in a thin gas) is representable.
        let correction = if product == 0.0 { 1.0 } else { -log1p(-product) / product };
        (absorbed_opacity / tau) * exposed * correction
    } else {
        let small = tau.min(dose);
        (small - log1p(exp(-(tau - dose).abs()) * (-expm1(-small)))) / tau
    };
    let (fraction, remaining_fraction, final_depth) = if candidate <= 0.5 {
        let remaining = 1.0 - candidate;
        (candidate, remaining, tau * remaining)
    } else {
        // Each term is nonnegative. Compute the surviving small side directly;
        // the absorbed candidate can otherwise round a little above one.
        let final_depth = if tau >= dose {
            let gap = tau - dose;
            gap + log1p(exp(-gap) * (-expm1(-dose)))
        } else {
            log1p(exp(-(dose - tau)) * (-expm1(-tau)))
        };
        let remaining = final_depth / tau;
        (1.0 - remaining, remaining, final_depth)
    };
    if !nonnegative(fraction) || fraction > 1.0 || !nonnegative(remaining_fraction)
        || remaining_fraction > 1.0 || !nonnegative(final_depth)
    {
        return Err(PhotoionizationError::Arithmetic("analytic photon depletion"));
    }
    Ok((fraction, remaining_fraction, final_depth))
}

fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn nonnegative(x: f64) -> bool { x.is_finite() && x >= 0.0 }
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), PhotoionizationError> {
    if cancelled() { Err(PhotoionizationError::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_depletion_matches_independent_high_precision_optical_depth_oracles() {
        // Independent 100-decimal-digit fixed-point exp/log evaluation.
        for (tau, dose, absorbed, neutral) in [
            (1e-12, 1e-12, 9.99999999999e-13, 0.999999999999),
            (1e-12, 1.0, 0.6321205588284414, 0.3678794411715586),
            (1e-12, 50.0, 1.0, 1.9287498479648822e-22),
            (1.0, 1.0, 0.51011987435525, 0.48988012564475),
            (50.0, 50.0, 0.9861370563888011, 0.013862943611198906),
            (1e6, 1e-12, 1e-18, 1.0),
            (1e6, 999999.0, 0.9999986867383125, 1.3132616875182228e-6),
            (1e6, 1e6, 0.9999993068528195, 6.931471805599453e-7),
            (1e6, 1000001.0, 0.9999996867383125, 3.1326168751822283e-7),
            (1e-200, 1e-200, 1e-200, 1.0),
        ] {
            let (a, n, _) = depletion(tau, dose).unwrap();
            assert!((a - absorbed).abs() <= 5e-13 * absorbed, "{tau} {dose}: {a}");
            assert!((n - neutral).abs() <= 5e-13 * neutral, "{tau} {dose}: {n}");
            assert!((a + n - 1.0).abs() <= 2e-15);
        }
    }
}
