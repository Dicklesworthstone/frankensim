//! Simultaneous discrete-spectrum photoionization of one uniformly mixed cell.
//!
//! Each group retains its own cross section, photon inventory and excess
//! photoelectron energy. For neutral amount n, its loss is
//! `A*F_g*(1-exp(-sigma_g*NA*n/A))/NA`. Averaging photon energies or cross
//! sections would change both shielding and heating as this opacity evolves.
//!
//! Symmetric splitting composes the parent's exact single-group subflows.
//! Full/two-half-step differences, divided by three for this second-order
//! method, control estimated local population, temperature AND group-photon
//! errors. Positive substeps preserve populations and photon budgets even in
//! stiff depletion; there is no nonlinear solve. Cross sections are prepared
//! once, and every trial reuses bounded buffers without allocating.
//! See Strang (1968), <https://doi.org/10.1137/0705041> for the composition.
//!
//! Groups are discrete beams along the same column, with constant incident
//! fluxes over the exposure. This is neither a spectral quadrature nor a
//! resolved front. The parent's fixed geometry, immediate thermalization,
//! stationary heavy targets and absence of other absorbers/recombination apply.
//! Splitting error estimates are not global or empirical error certificates.

pub mod column;

use super::*;
use fs_math::dd::Dd;
use fs_thermochem::kinetics::atomic_ionization::photoionization::MAX_PHOTON_GROUPS_V1;
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};

/// Source interval, local integration controls and whole-exposure allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralPhotoionizationSolve {
    /// Supported electron temperature lower bound [K].
    pub minimum_electron_temperature_k: f64,
    /// Supported electron temperature upper bound [K].
    pub maximum_electron_temperature_k: f64,
    /// Estimated local relative amount, temperature and per-group photon error.
    pub relative_tolerance: f64,
    /// Absolute amount error scale, including photon accounts divided by NA [mol].
    pub absolute_amount_tolerance_mol: f64,
    /// Absolute local electron temperature error scale [K].
    pub temperature_tolerance_k: f64,
    /// Final original-source chemical plus sensible minus absorbed energy [J].
    pub energy_tolerance_j: f64,
    /// Carrier-scaled charge admission, in (0,1).
    pub relative_charge_tolerance: f64,
    /// Per-group photon conservation allowance relative to incoming photons.
    pub relative_photon_tolerance: f64,
    /// First proposed step [s].
    pub initial_step_s: f64,
    /// Smallest retry [s]; a shorter final remainder is allowed.
    pub minimum_step_s: f64,
    /// Largest proposed step [s].
    pub maximum_step_s: f64,
    /// Maximum accepted plus rejected attempts.
    pub maximum_steps: usize,
    /// Maximum admitted groups, capped by the original L1 source policy.
    pub maximum_groups: usize,
}
impl SpectralPhotoionizationSolve {
    /// Explicit source interval with numerical defaults, not measured accuracy.
    pub fn new(minimum: Temperature, maximum: Temperature) -> Self {
        Self {
            minimum_electron_temperature_k: minimum.value(),
            maximum_electron_temperature_k: maximum.value(),
            relative_tolerance: 1e-5,
            absolute_amount_tolerance_mol: 1e-14,
            temperature_tolerance_k: 1e-6,
            energy_tolerance_j: 1e-7,
            relative_charge_tolerance: 1e-10,
            relative_photon_tolerance: 1e-12,
            initial_step_s: 1e-6,
            minimum_step_s: 1e-24,
            maximum_step_s: 1.0,
            maximum_steps: 100_000,
            maximum_groups: MAX_PHOTON_GROUPS_V1,
        }
    }
    fn validate(self) -> Result<(), SpectralPhotoionizationError> {
        if [
            self.minimum_electron_temperature_k,
            self.maximum_electron_temperature_k,
            self.relative_tolerance,
            self.absolute_amount_tolerance_mol,
            self.temperature_tolerance_k,
            self.energy_tolerance_j,
            self.relative_charge_tolerance,
            self.relative_photon_tolerance,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .iter()
        .any(|x| !positive(*x))
            || self.minimum_electron_temperature_k >= self.maximum_electron_temperature_k
            || self.relative_tolerance > 0.1
            || self.relative_charge_tolerance >= 1.0
            || self.relative_photon_tolerance >= 1.0
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=1_000_000).contains(&self.maximum_steps)
            || !(1..=MAX_PHOTON_GROUPS_V1).contains(&self.maximum_groups)
        {
            return Err(SpectralPhotoionizationError::Invalid(
                "spectral integration controls",
            ));
        }
        Ok(())
    }
}

/// Atomic refusal: no caller-owned source or gas state is modified.
#[derive(Clone, Debug, PartialEq)]
pub enum SpectralPhotoionizationError {
    /// Invalid controls, exposure or original bank binding.
    Invalid(&'static str),
    /// Explicit group budget exceeded before allocation.
    GroupLimit { requested: usize, maximum: usize },
    /// Source energy, cross section or numerical subflow refused in this group.
    Group {
        index: usize,
        error: PhotoionizationError,
    },
    /// Original electron caloric support does not span the declared interval.
    Admission(KineticsErrorV1),
    /// Electron temperature left the declared source interval.
    OutsideDomain,
    /// Step count, minimum step or representable time increment exhausted.
    StepBudget,
    /// Required physical value or bounded result allocation is not representable.
    Arithmetic(&'static str),
    /// Absorbed plus transmitted photons failed this group's incoming account.
    PhotonBalance { group: usize, residual_photons: f64 },
    /// Surviving neutrals plus absorbed photon amount failed the original inventory.
    AmountBalance { residual_mol: f64 },
    /// Actual chemical plus electron sensible energy failed absorbed energy closure.
    EnergyBalance { residual_j: f64 },
    /// Original source state or temperature publication refused.
    State(FrozenIonizedGasStateErrorV1),
    /// Cancellation observed before complete publication.
    Cancelled,
}
impl core::fmt::Display for SpectralPhotoionizationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spectral gas photoionization: {self:?}")
    }
}
impl std::error::Error for SpectralPhotoionizationError {}

/// Original-source state and frequency-resolved photon/energy accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralPhotoionizationResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    source: AtomicPhotoionizationV1,
    groups: Vec<PhotonNumberFluxV1>,
    group_balances: Vec<PhotoionizationBalance>,
    area: Area,
    elapsed: Time,
    controls: SpectralPhotoionizationSolve,
    extent: Amount,
    balance: PhotoionizationBalance,
    accepted_steps: usize,
    rejected_steps: usize,
    maximum_normalized_local_error: f64,
}
impl SpectralPhotoionizationResult {
    /// Actual endpoint for existing collision, electrical and optical consumers.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Original source-bound input.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Original cross section, atomic levels and thermodynamic bank.
    pub const fn source(&self) -> &AtomicPhotoionizationV1 {
        &self.source
    }
    /// Original discrete energy/flux groups in their declared deterministic order.
    pub fn groups(&self) -> &[PhotonNumberFluxV1] {
        &self.groups
    }
    /// Per-group budgets in the same order. Group sensible changes are each
    /// group's deposited excess energy; aggregate sensible change is measured
    /// from the actual final electron temperature and population.
    pub fn group_balances(&self) -> &[PhotoionizationBalance] {
        &self.group_balances
    }
    /// Common area normal to the incident beams [m²].
    pub const fn illuminated_area(&self) -> Area {
        self.area
    }
    /// Completed exposure [s].
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    /// Applied source and integration controls.
    pub const fn controls(&self) -> SpectralPhotoionizationSolve {
        self.controls
    }
    /// One ion/electron pair per absorbed photon [mol].
    pub const fn reaction_extent(&self) -> Amount {
        self.extent
    }
    /// Summed photons and original-state gas-energy account.
    pub const fn balance(&self) -> PhotoionizationBalance {
        self.balance
    }
    /// Complete accepted steps, each using two symmetric half steps.
    pub const fn accepted_steps(&self) -> usize {
        self.accepted_steps
    }
    /// Rejected attempts counted against the same work budget.
    pub const fn rejected_steps(&self) -> usize {
        self.rejected_steps
    }
    /// Largest accepted local full/two-half difference divided by its allowance.
    pub const fn maximum_normalized_local_error(&self) -> f64 {
        self.maximum_normalized_local_error
    }
}

#[derive(Clone, Copy)]
struct Group {
    opacity_per_mol: f64,
    dose_rate: f64,
    photon_rate: f64,
    excess_j_mol: f64,
    incident_photons: f64,
}
#[derive(Clone, Copy)]
struct Point {
    neutral: f64,
    extent: Dd,
    deposited: Dd,
}
impl Point {
    fn temperature(self, initial: &FrozenIonizedGasStateV1, electron: usize) -> f64 {
        let old = initial.electron_temperature().value();
        let extent = self.extent.to_f64();
        if extent == 0.0 {
            return old;
        }
        let electrons = initial.amounts()[electron].value();
        let total = electrons + extent;
        let born = (self.deposited.to_f64() / extent) / (1.5 * R);
        if extent <= electrons {
            old + (extent / total) * (born - old)
        } else {
            (electrons / total) * old + (extent / total) * born
        }
    }
}

/// Evolve all discrete photon groups simultaneously, with evolving band opacity.
///
/// Each group weight is a bandwidth-integrated number flux. Empty, dark and
/// subthreshold spectra preserve the input while retaining transmitted photons.
/// All source energies are admitted even for dark groups. Exact zero initial
/// electrons are supported. The returned temperature is never clipped or
/// projected to equilibrium, and no mean-energy replacement spectrum is used.
pub fn advance_spectral_photoionization(
    initial: &FrozenIonizedGasStateV1,
    source: &AtomicPhotoionizationV1,
    groups: &[PhotonNumberFluxV1],
    area: Area,
    duration: Time,
    controls: SpectralPhotoionizationSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SpectralPhotoionizationResult, SpectralPhotoionizationError> {
    check_cancel(&mut cancelled)?;
    controls.validate()?;
    if groups.len() > controls.maximum_groups {
        return Err(SpectralPhotoionizationError::GroupLimit {
            requested: groups.len(),
            maximum: controls.maximum_groups,
        });
    }
    let channel = source.channel();
    let [a, i, e] = [
        channel.neutral_index(),
        channel.ion_index(),
        channel.electron_index(),
    ];
    if !positive(area.value())
        || !nonnegative(duration.value())
        || initial.thermodynamics() != channel.thermodynamics()
        || initial.electron_index() != Some(e)
        || initial.relative_charge_residual() > controls.relative_charge_tolerance
    {
        return Err(SpectralPhotoionizationError::Invalid(
            "exposure or original source bank",
        ));
    }
    check_temperature(initial.electron_temperature().value(), controls)?;
    admit_caloric_interval(
        &initial.thermodynamics().models()[e],
        Temperature::new(controls.minimum_electron_temperature_k),
        Temperature::new(controls.maximum_electron_temperature_k),
        &mut cancelled,
    )
    .map_err(|error| match error {
        KineticsErrorV1::Cancelled => SpectralPhotoionizationError::Cancelled,
        error => SpectralPhotoionizationError::Admission(error),
    })?;
    let mut prepared = Vec::new();
    prepared
        .try_reserve_exact(groups.len())
        .map_err(|_| SpectralPhotoionizationError::Arithmetic("group allocation"))?;
    let cell_controls = PhotoionizationControls {
        relative_charge_tolerance: controls.relative_charge_tolerance,
        energy_tolerance_j: controls.energy_tolerance_j,
    };
    let chi = channel.ionization_energy_j_mol();
    for (index, &beam) in groups.iter().enumerate() {
        check_cancel(&mut cancelled)?;
        let (_, incident) = validate_exposure(
            PhotoionizationExposure {
                beam,
                illuminated_area: area,
                duration,
            },
            cell_controls,
        )
        .map_err(|error| SpectralPhotoionizationError::Group { index, error })?;
        let sigma = source
            .cross_section(beam.photon_energy)
            .map_err(|error| SpectralPhotoionizationError::Group {
                index,
                error: PhotoionizationError::Source(error),
            })?
            .value();
        let group = Group {
            opacity_per_mol: (sigma * NA) / area.value(),
            dose_rate: sigma * beam.number_flux_per_m2_s,
            photon_rate: beam.number_flux_per_m2_s * area.value(),
            excess_j_mol: if sigma == 0.0 {
                0.0
            } else {
                (beam.photon_energy.value() - chi / NA) * NA
            },
            incident_photons: incident,
        };
        if [
            group.opacity_per_mol,
            group.dose_rate,
            group.photon_rate,
            group.excess_j_mol,
        ]
        .iter()
        .any(|x| !nonnegative(*x))
        {
            return Err(SpectralPhotoionizationError::Arithmetic(
                "prepared photon group",
            ));
        }
        prepared.push(group);
    }
    let mut full_accounts = account_buffer(groups.len())?;
    let mut half_accounts = account_buffer(groups.len())?;
    let mut accounts = account_buffer(groups.len())?;
    let mut point = Point {
        neutral: initial.amounts()[a].value(),
        extent: Dd::ZERO,
        deposited: Dd::ZERO,
    };
    let mut time = 0.0;
    let mut step = controls.initial_step_s;
    let (mut accepted, mut rejected, mut maximum_error) = (0, 0, 0.0_f64);
    if point.neutral == 0.0 || prepared.iter().all(|group| group.dose_rate == 0.0) {
        for (account, group) in accounts.iter_mut().zip(&prepared) {
            account[1] = Dd::from_f64(group.incident_photons);
        }
        time = duration.value();
    }
    while time < duration.value() {
        check_cancel(&mut cancelled)?;
        if accepted + rejected >= controls.maximum_steps {
            return Err(SpectralPhotoionizationError::StepBudget);
        }
        let remaining = duration.value() - time;
        let dt = step.min(remaining);
        let end = if dt == remaining {
            duration.value()
        } else {
            time + dt
        };
        if end <= time || dt * 0.5 == 0.0 {
            return Err(SpectralPhotoionizationError::StepBudget);
        }
        let dt = end - time;
        full_accounts.fill([Dd::ZERO; 2]);
        half_accounts.fill([Dd::ZERO; 2]);
        let mut full = point;
        split_step(&mut full, dt, &prepared, &mut full_accounts, &mut cancelled)?;
        let mut half = point;
        split_step(
            &mut half,
            dt * 0.5,
            &prepared,
            &mut half_accounts,
            &mut cancelled,
        )?;
        let middle_temperature = half.temperature(initial, e);
        split_step(
            &mut half,
            dt * 0.5,
            &prepared,
            &mut half_accounts,
            &mut cancelled,
        )?;
        let amount_error = |x: f64, y: f64| {
            (x - y).abs()
                / (3.0
                    * (controls.absolute_amount_tolerance_mol
                        + controls.relative_tolerance * x.abs().max(y.abs())))
        };
        let t_full = full.temperature(initial, e);
        let t_half = half.temperature(initial, e);
        let mut error = amount_error(full.neutral, half.neutral)
            .max(amount_error(full.extent.to_f64(), half.extent.to_f64()))
            .max(
                (t_full - t_half).abs()
                    / (3.0
                        * (controls.temperature_tolerance_k
                            + controls.relative_tolerance * t_full.abs().max(t_half.abs()))),
            );
        for (f, h) in full_accounts.iter().zip(&half_accounts) {
            for k in 0..2 {
                error = error.max(amount_error(f[k].to_f64() / NA, h[k].to_f64() / NA));
            }
        }
        if !error.is_finite() || !positive(t_full) || !positive(t_half) {
            return Err(SpectralPhotoionizationError::Arithmetic(
                "spectral step error or temperature",
            ));
        }
        if error <= 1.0 {
            check_temperature(middle_temperature, controls)?;
            check_temperature(t_half, controls)?;
            for (sum, local) in accounts.iter_mut().zip(&half_accounts) {
                for k in 0..2 {
                    sum[k] = sum[k] + local[k];
                }
            }
            point = half;
            time = end;
            accepted += 1;
            maximum_error = maximum_error.max(error);
            let factor = if error == 0.0 {
                2.0
            } else {
                (0.85 * exp(-ln(error) / 3.0)).clamp(0.25, 2.0)
            };
            step = (dt * factor).clamp(controls.minimum_step_s, controls.maximum_step_s);
        } else {
            rejected += 1;
            step = dt * (0.8 * exp(-ln(error) / 3.0)).clamp(0.1, 0.5);
            if step < controls.minimum_step_s {
                return Err(SpectralPhotoionizationError::StepBudget);
            }
        }
    }
    let mut group_balances = Vec::new();
    group_balances
        .try_reserve_exact(groups.len())
        .map_err(|_| SpectralPhotoionizationError::Arithmetic("group result allocation"))?;
    let mut totals = [Dd::ZERO; 6];
    for (index, ((group, prepared), account)) in
        groups.iter().zip(&prepared).zip(&accounts).enumerate()
    {
        let absorbed = account[0].to_f64();
        let transmitted = account[1].to_f64();
        let photon_residual =
            (account[0] + account[1] - Dd::from_f64(prepared.incident_photons)).to_f64();
        if !photon_residual.is_finite()
            || photon_residual.abs()
                > controls.relative_photon_tolerance * prepared.incident_photons
        {
            return Err(SpectralPhotoionizationError::PhotonBalance {
                group: index,
                residual_photons: photon_residual,
            });
        }
        let chemical = (absorbed / NA) * chi;
        let sensible = (absorbed / NA) * prepared.excess_j_mol;
        let energy = absorbed * group.photon_energy.value();
        let values = [
            prepared.incident_photons,
            absorbed,
            transmitted,
            chemical,
            sensible,
            energy,
        ];
        if values.iter().any(|x| !nonnegative(*x)) {
            return Err(SpectralPhotoionizationError::Arithmetic(
                "group photon or energy balance",
            ));
        }
        for (sum, value) in totals.iter_mut().zip(values) {
            *sum = *sum + Dd::from_f64(value);
        }
        group_balances.push(PhotoionizationBalance {
            incident_photons: values[0],
            absorbed_photons: absorbed,
            transmitted_photons: transmitted,
            chemical_storage_change_j: chemical,
            electron_sensible_energy_change_j: sensible,
            absorbed_energy_j: energy,
            energy_residual_j: sensible + (chemical - energy),
        });
    }
    let extent = point.extent.to_f64();
    let amount_residual = (Dd::from_f64(point.neutral) + point.extent
        - Dd::from_f64(initial.amounts()[a].value()))
    .to_f64();
    if !amount_residual.is_finite()
        || amount_residual.abs() > controls.absolute_amount_tolerance_mol
    {
        return Err(SpectralPhotoionizationError::AmountBalance {
            residual_mol: amount_residual,
        });
    }
    let next_te = point.temperature(initial, e);
    let sensible = 1.5
        * R
        * (initial.amounts()[e].value() * (next_te - initial.electron_temperature().value())
            + extent * next_te);
    let residual = (Dd::from_f64(sensible) + totals[3] - totals[5]).to_f64();
    if totals.iter().any(|x| !x.to_f64().is_finite()) || !residual.is_finite() {
        return Err(SpectralPhotoionizationError::Arithmetic(
            "total photon or energy balance",
        ));
    }
    if residual.abs() > controls.energy_tolerance_j {
        return Err(SpectralPhotoionizationError::EnergyBalance {
            residual_j: residual,
        });
    }
    let mut amounts = initial.amounts().to_vec();
    amounts[a] = Amount::new(point.neutral);
    amounts[i] = Amount::new(initial.amounts()[i].value() + extent);
    amounts[e] = Amount::new(initial.amounts()[e].value() + extent);
    let state = FrozenIonizedGasStateV1::new(
        channel.thermodynamics(),
        &amounts,
        initial.volume(),
        initial.heavy_temperature(),
        Temperature::new(next_te),
        controls.relative_charge_tolerance,
        &mut cancelled,
    )
    .map_err(|error| match error {
        FrozenIonizedGasStateErrorV1::Cancelled => SpectralPhotoionizationError::Cancelled,
        error => SpectralPhotoionizationError::State(error),
    })?;
    check_cancel(&mut cancelled)?;
    Ok(SpectralPhotoionizationResult {
        initial: initial.clone(),
        state,
        source: source.clone(),
        groups: groups.to_vec(),
        group_balances,
        area,
        elapsed: duration,
        controls,
        extent: Amount::new(extent),
        balance: PhotoionizationBalance {
            incident_photons: totals[0].to_f64(),
            absorbed_photons: totals[1].to_f64(),
            transmitted_photons: totals[2].to_f64(),
            chemical_storage_change_j: totals[3].to_f64(),
            electron_sensible_energy_change_j: sensible,
            absorbed_energy_j: totals[5].to_f64(),
            energy_residual_j: residual,
        },
        accepted_steps: accepted,
        rejected_steps: rejected,
        maximum_normalized_local_error: maximum_error,
    })
}

fn split_step(
    point: &mut Point,
    dt: f64,
    groups: &[Group],
    accounts: &mut [[Dd; 2]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SpectralPhotoionizationError> {
    // Palindromic group order. Merge the two identical central half-steps.
    for index in 0..groups.len() {
        let h = if index + 1 == groups.len() {
            dt
        } else {
            dt * 0.5
        };
        group_step(
            point,
            index,
            groups[index],
            h,
            &mut accounts[index],
            cancelled,
        )?;
    }
    for index in (0..groups.len().saturating_sub(1)).rev() {
        group_step(
            point,
            index,
            groups[index],
            dt * 0.5,
            &mut accounts[index],
            cancelled,
        )?;
    }
    Ok(())
}
fn group_step(
    point: &mut Point,
    index: usize,
    group: Group,
    dt: f64,
    account: &mut [Dd; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SpectralPhotoionizationError> {
    check_cancel(cancelled)?;
    let photons = group.photon_rate * dt;
    let tau = group.opacity_per_mol * point.neutral;
    let dose = group.dose_rate * dt;
    if !nonnegative(photons) || !nonnegative(tau) || !nonnegative(dose) {
        return Err(SpectralPhotoionizationError::Arithmetic(
            "group substep exposure",
        ));
    }
    if group.dose_rate == 0.0 || point.neutral == 0.0 {
        account[1] = account[1] + Dd::from_f64(photons);
        return Ok(());
    }
    if tau == 0.0 || dose == 0.0 {
        return Err(SpectralPhotoionizationError::Arithmetic(
            "positive group optical depth or dose underflow",
        ));
    }
    let (lost_fraction, retained_fraction, _) = depletion(tau, dose)
        .map_err(|error| SpectralPhotoionizationError::Group { index, error })?;
    let lost = point.neutral * lost_fraction;
    let retained = point.neutral * retained_fraction;
    let transmitted = transmitted_photons(photons, tau, dose);
    if !positive(lost) || !nonnegative(retained) || !nonnegative(transmitted) {
        return Err(SpectralPhotoionizationError::Arithmetic(
            "group photon depletion",
        ));
    }
    point.neutral = retained;
    point.extent = point.extent + Dd::from_f64(lost);
    point.deposited = point.deposited + Dd::from_f64(lost * group.excess_j_mol);
    account[0] = account[0] + Dd::from_f64(lost * NA);
    account[1] = account[1] + Dd::from_f64(transmitted);
    Ok(())
}
fn account_buffer(count: usize) -> Result<Vec<[Dd; 2]>, SpectralPhotoionizationError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| SpectralPhotoionizationError::Arithmetic("spectral work allocation"))?;
    result.resize(count, [Dd::ZERO; 2]);
    Ok(result)
}
fn check_temperature(
    t: f64,
    controls: SpectralPhotoionizationSolve,
) -> Result<(), SpectralPhotoionizationError> {
    if t.is_finite()
        && (controls.minimum_electron_temperature_k..=controls.maximum_electron_temperature_k)
            .contains(&t)
    {
        Ok(())
    } else {
        Err(SpectralPhotoionizationError::OutsideDomain)
    }
}
fn check_cancel(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpectralPhotoionizationError> {
    if cancelled() {
        Err(SpectralPhotoionizationError::Cancelled)
    } else {
        Ok(())
    }
}
