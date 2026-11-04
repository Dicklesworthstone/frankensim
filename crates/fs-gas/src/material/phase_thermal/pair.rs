//! Two finite original-law inventories joined by one passive thermal contact.

use super::{
    CaloricReport, Error, IsobaricMixtureEquilibriumParcel, IsobaricMixtureHeatStep,
    IsochoricMixtureEquilibriumParcel, MixtureEquilibriumHeatStep, PhaseEnergyVolumeControls,
    PhaseEnthalpyControls, ThermalParcel, Time, WallControl, Work, poll, positive, solve_scalar,
};

/// Bounded heat exchange between two finite phase-equilibrium inventories.
///
/// Each side retains its own original species, phase policy and caloric
/// constraint. Positive heat enters the left parcel and leaves the right one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixturePairThermalControl<L = PhaseEnergyVolumeControls, R = PhaseEnergyVolumeControls> {
    /// Signed heat into the left parcel [J]. Increasing finite endpoints must
    /// contain zero. The passive initial-flux bound further limits this bracket;
    /// every remaining endpoint must be supported by BOTH original sources.
    pub heat_bracket_j: [f64; 2],
    /// Original left source bounds, phase policy and per-trial work limits.
    pub left_equilibrium: L,
    /// Original right source bounds, phase policy and per-trial work limits.
    pub right_equilibrium: R,
    /// Physical heat/energy residual allowance and scalar iteration limit.
    pub thermal: WallControl,
    /// Maximum admitted negative combined source entropy change [J/K].
    /// This is a numerical allowance, never a physical entropy source.
    pub entropy_tolerance_j_k: f64,
    /// Total paired heat-coordinate queries, including bracket endpoints.
    pub maximum_heat_evaluations: usize,
    /// Total original PT selections across BOTH sides and all heat trials.
    pub maximum_pt_evaluations: usize,
}

impl<L, R> MixturePairThermalControl<L, R> {
    /// Declare the heat interval and independent original caloric controls.
    #[must_use]
    pub fn new(heat_bracket_j: [f64; 2], left_equilibrium: L, right_equilibrium: R) -> Self {
        Self {
            heat_bracket_j,
            left_equilibrium,
            right_equilibrium,
            thermal: WallControl {
                absolute_energy_tolerance_j: 1e-5,
                relative_energy_tolerance: 1e-9,
                maximum_iterations: 96,
            },
            entropy_tolerance_j_k: 1e-7,
            maximum_heat_evaluations: 132,
            maximum_pt_evaluations: 65_536,
        }
    }

    fn validate_common(&self) -> Result<(), Error> {
        self.thermal.validate()?;
        let [lo, hi] = self.heat_bracket_j;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo > 0.0
            || hi < 0.0
            || lo >= hi
            || !self.entropy_tolerance_j_k.is_finite()
            || self.entropy_tolerance_j_k < 0.0
            || !(1..=132).contains(&self.maximum_heat_evaluations)
            || !(1..=1_048_576).contains(&self.maximum_pt_evaluations)
        {
            return Err(Error::Invalid("mixture-pair thermal controls"));
        }
        Ok(())
    }
}

impl MixturePairThermalControl {
    /// Validate both U,V controls and the common physical/work allowances.
    pub fn validate(self) -> Result<(), Error> {
        self.left_equilibrium.validate()?;
        self.right_equilibrium.validate()?;
        self.validate_common()
    }
}

/// Two independently prescribed constant-pressure phase inventories.
pub type IsobaricMixturePairThermalControl =
    MixturePairThermalControl<PhaseEnthalpyControls, PhaseEnthalpyControls>;

impl IsobaricMixturePairThermalControl {
    /// Validate both PH controls and the common physical/work allowances.
    pub fn validate(self) -> Result<(), Error> {
        self.left_equilibrium.validate()?;
        self.right_equilibrium.validate()?;
        self.validate_common()
    }
}

/// A rigid left inventory and an independently prescribed-pressure right one.
pub type MixedMixturePairThermalControl =
    MixturePairThermalControl<PhaseEnergyVolumeControls, PhaseEnthalpyControls>;

impl MixedMixturePairThermalControl {
    /// Validate the U,V/PH controls and common physical/work allowances.
    pub fn validate(self) -> Result<(), Error> {
        self.left_equilibrium.validate()?;
        self.right_equilibrium.validate()?;
        self.validate_common()
    }
}

/// Accepted heat exchange, original source endpoints and complete work accounts.
/// Arrays are ordered `[left, right]`; neither parcel exchanges mass or species.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixturePairThermalStep<L = MixtureEquilibriumHeatStep, R = MixtureEquilibriumHeatStep> {
    /// Physical interval supplied to the backward-Euler contact law.
    pub elapsed: Time,
    /// Actual source heat into the left parcel [J]: Delta U or Delta H.
    pub heat_to_left_j: f64,
    /// Actual source heat into the right parcel [J], opposite to the left heat
    /// within the independently checked physical energy allowance.
    pub heat_to_right_j: f64,
    /// Independently source-integrated internal-energy changes [J].
    pub internal_energy_change_j: [f64; 2],
    /// Work exported by each parcel [J]: zero at fixed V, actual p Delta V at
    /// prescribed p. The two mechanical loads are separate from heat exchange.
    pub boundary_work_by_parcels_j: [f64; 2],
    /// Combined Delta U plus both outward boundary works [J], evaluated against
    /// the exactly opposed trial heat to avoid subtracting absolute energies.
    pub energy_residual_j: f64,
    /// Left source heat minus dt G times the endpoint temperature gap [J].
    /// The corresponding right equation is independently checked as well.
    pub thermal_residual_j: f64,
    /// Original reference-invariant entropy changes [J/K].
    pub entropy_change_j_k: [f64; 2],
    /// Combined source entropy change [J/K]. The isolated pair has no external
    /// heat or mass entropy transfer; quasistatic boundary work carries none.
    pub entropy_production_j_k: f64,
    /// Original accepted left source receipt; counts cover that trial only.
    pub left_equilibrium: L,
    /// Original accepted right source receipt; counts cover that trial only.
    pub right_equilibrium: R,
    /// Initial parcel temperatures [K].
    pub initial_temperature_k: [f64; 2],
    /// Accepted parcel temperatures [K].
    pub endpoint_temperature_k: [f64; 2],
    /// Safeguarded scalar iterations; zero for an exact resting contact.
    pub iterations: usize,
    /// Total paired heat-coordinate queries across this call.
    pub heat_evaluations: usize,
    /// Total original PT selections across both sources and all trials.
    pub pt_evaluations: usize,
    /// Total original pressure-root samples; PH participants contribute zero.
    pub pressure_evaluations: usize,
    /// Sum of original source EOS ceilings, not an observed evaluation count.
    pub eos_evaluation_upper_bound: usize,
}

/// Two constant-pressure source receipts and their joint thermal/work account.
pub type IsobaricMixturePairThermalStep =
    MixturePairThermalStep<IsobaricMixtureHeatStep, IsobaricMixtureHeatStep>;

/// A rigid left receipt and a constant-pressure right receipt in one account.
pub type MixedMixturePairThermalStep =
    MixturePairThermalStep<MixtureEquilibriumHeatStep, IsobaricMixtureHeatStep>;

/// Exchange heat between two finite rigid phase-equilibrium mixtures.
///
/// The original U,V owners determine both endpoint temperatures, pressures and
/// phase amounts. One bounded scalar solve closes
/// `Q = dt G [T_right(-Q) - T_left(Q)]`; no effective Cp/Cv or extra latent source
/// is introduced. `G` is the complete nonnegative interface conductance [W/K],
/// held over this interval, and positive Q heats the left inventory.
///
/// Both original energies, equal/opposite source heats, endpoint heat equations
/// and combined original entropy change are checked before either state is
/// published. Source refusal, work exhaustion and cancellation preserve both
/// inputs. Zero duration/conductance or equal initial temperatures preserve all
/// state bits with no phase evaluations.
///
/// This is first-order lumped equilibrium heat exchange: it needs time refinement
/// and uniform-inventory justification. It resolves neither interfacial phase
/// kinetics nor spatial temperature fields. Original finite phase evidence is
/// retained without a new global-equilibrium or uniqueness claim.
pub fn exchange_isochoric_mixture_pair(
    left: &mut IsochoricMixtureEquilibriumParcel,
    right: &mut IsochoricMixtureEquilibriumParcel,
    conductance_w_per_k: f64,
    duration: Time,
    controls: MixturePairThermalControl,
    cancelled: impl FnMut() -> bool,
) -> Result<MixturePairThermalStep, Error> {
    exchange(
        left,
        right,
        conductance_w_per_k,
        duration,
        controls,
        cancelled,
    )
}

/// Exchange heat between two finite mixtures at their own fixed pressures.
///
/// Uses the same implicit contact as [`exchange_isochoric_mixture_pair`], with
/// native PH owners on both sides. Each source heat is its original Delta H;
/// the combined first law independently includes both Delta U values and both
/// actual outward p Delta V works. Pressures, species and masses remain fixed.
/// Each load is quasistatic; piston inertia and changing contact geometry are
/// outside this held-conductance interval. All pair acceptance/refusal rules
/// of the rigid interface apply unchanged.
pub fn exchange_isobaric_mixture_pair(
    left: &mut IsobaricMixtureEquilibriumParcel,
    right: &mut IsobaricMixtureEquilibriumParcel,
    conductance_w_per_k: f64,
    duration: Time,
    controls: IsobaricMixturePairThermalControl,
    cancelled: impl FnMut() -> bool,
) -> Result<IsobaricMixturePairThermalStep, Error> {
    exchange(
        left,
        right,
        conductance_w_per_k,
        duration,
        controls,
        cancelled,
    )
}

/// Exchange heat between a rigid left mixture and a constant-pressure right one.
///
/// The same heat enters native U on the left and leaves native H on the right;
/// only the right mechanical load does boundary work. Positive Q heats the rigid
/// inventory, negative Q heats the expanding inventory. Both source laws and all
/// physical/work/cancellation guarantees of [`exchange_isochoric_mixture_pair`]
/// are retained, including latent redistribution and original entropy changes.
pub fn exchange_isochoric_isobaric_mixture_pair(
    left: &mut IsochoricMixtureEquilibriumParcel,
    right: &mut IsobaricMixtureEquilibriumParcel,
    conductance_w_per_k: f64,
    duration: Time,
    controls: MixedMixturePairThermalControl,
    cancelled: impl FnMut() -> bool,
) -> Result<MixedMixturePairThermalStep, Error> {
    exchange(
        left,
        right,
        conductance_w_per_k,
        duration,
        controls,
        cancelled,
    )
}

struct Trial<L: ThermalParcel, R: ThermalParcel> {
    heat: f64,
    left: L,
    right: R,
    left_equilibrium: CaloricReport<L::Report>,
    right_equilibrium: CaloricReport<R::Report>,
}

fn resolve<P: ThermalParcel>(
    parcel: &mut P,
    heat: f64,
    controls: P::Control,
    maximum_pt: usize,
    work: &mut Work,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<CaloricReport<P::Report>, Error> {
    let remaining = maximum_pt
        .checked_sub(work.pt)
        .ok_or_else(|| work.error())?;
    if heat != 0.0 && remaining == 0 {
        return Err(work.error());
    }
    let source = parcel
        .resolve(heat, controls, remaining, cancelled)
        .map_err(|e| work.source_error(e))?;
    work.add(source.pt, source.pressure, source.eos)?;
    if work.pt > maximum_pt {
        return Err(work.error());
    }
    Ok(source)
}

fn trial<L: ThermalParcel, R: ThermalParcel>(
    left: &L,
    right: &R,
    heat: f64,
    controls: MixturePairThermalControl<L::Control, R::Control>,
    work: &mut Work,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<Trial<L, R>, Error> {
    if work.heat == controls.maximum_heat_evaluations {
        return Err(work.error());
    }
    work.heat += 1;
    let (mut left, mut right) = (left.clone(), right.clone());
    let left_equilibrium = resolve(
        &mut left,
        heat,
        controls.left_equilibrium,
        controls.maximum_pt_evaluations,
        work,
        cancelled,
    )?;
    // Charge the left solve before admitting the right one. It receives only
    // the work still available to the complete pair, never a second full budget.
    let right_equilibrium = resolve(
        &mut right,
        if heat == 0.0 { 0.0 } else { -heat },
        controls.right_equilibrium,
        controls.maximum_pt_evaluations,
        work,
        cancelled,
    )?;
    Ok(Trial {
        heat,
        left,
        right,
        left_equilibrium,
        right_equilibrium,
    })
}

fn exchange<L: ThermalParcel, R: ThermalParcel>(
    left: &mut L,
    right: &mut R,
    conductance: f64,
    duration: Time,
    mut controls: MixturePairThermalControl<L::Control, R::Control>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MixturePairThermalStep<L::Report, R::Report>, Error> {
    poll(&mut cancelled)?;
    L::validate(controls.left_equilibrium)?;
    R::validate(controls.right_equilibrium)?;
    controls.validate_common()?;
    if !conductance.is_finite()
        || conductance < 0.0
        || !duration.value().is_finite()
        || duration.value() < 0.0
    {
        return Err(Error::Invalid(
            "nonnegative finite mixture-pair conductance and duration required",
        ));
    }
    let initial_temperature_k = [left.temperature(), right.temperature()];
    let coefficient = duration.value() * conductance;
    let initial_gap = initial_temperature_k[1] - initial_temperature_k[0];
    let initial_flux = coefficient * initial_gap;
    if !coefficient.is_finite()
        || !initial_flux.is_finite()
        || (coefficient == 0.0 && duration.value() != 0.0 && conductance != 0.0)
        || (initial_flux == 0.0 && coefficient != 0.0 && initial_gap != 0.0)
    {
        return Err(Error::Invalid(
            "mixture-pair thermal flux is unrepresentable",
        ));
    }
    let mut work = Work::default();
    let mut accepted = None;
    let (heat, iterations) = if initial_flux == 0.0 {
        accepted = Some(trial(
            left,
            right,
            0.0,
            controls,
            &mut work,
            &mut cancelled,
        )?);
        (0.0, 0)
    } else {
        controls.thermal.absolute_energy_tolerance_j = controls
            .thermal
            .absolute_energy_tolerance_j
            .min(0.125 * initial_flux.abs());
        // Reserve an independent eighth of the absolute physical allowance for
        // each native source, so caloric error cannot consume the contact budget.
        let source_tolerance_j = 0.125 * controls.thermal.absolute_energy_tolerance_j;
        L::refine(
            &mut controls.left_equilibrium,
            source_tolerance_j / left.mass(),
        );
        R::refine(
            &mut controls.right_equilibrium,
            source_tolerance_j / right.mass(),
        );
        L::validate(controls.left_equilibrium)?;
        R::validate(controls.right_equilibrium)?;
        controls.validate_common()?;
        let lo = controls.heat_bracket_j[0].max(initial_flux.min(0.0));
        let hi = controls.heat_bracket_j[1].min(initial_flux.max(0.0));
        if lo >= hi {
            return Err(Error::Invalid("mixture-pair thermal heat bracket is empty"));
        }
        let mut previous: Option<(f64, f64)> = None;
        let mut equation = |q: f64, cancel: &mut dyn FnMut() -> bool| {
            let candidate = trial(left, right, q, controls, &mut work, cancel)?;
            let flux = coefficient * (candidate.right.temperature() - candidate.left.temperature());
            let residual = q - flux;
            if !residual.is_finite() {
                return Err(Error::Invalid("mixture-pair thermal trial flux"));
            }
            // A positive numerical secant accelerates the common safeguarded
            // root; it has no heat-capacity or physical-tangent authority.
            let slope = previous.and_then(|(old_q, old_r)| {
                let slope = (residual - old_r) / (q - old_q);
                positive(slope).then_some(slope)
            });
            previous = Some((q, residual));
            accepted = Some(candidate);
            Ok((residual, slope, q.abs().max(flux.abs())))
        };
        let scalar = WallControl {
            absolute_energy_tolerance_j: 0.5 * controls.thermal.absolute_energy_tolerance_j,
            relative_energy_tolerance: 0.5 * controls.thermal.relative_energy_tolerance,
            ..controls.thermal
        };
        scalar.validate()?;
        solve_scalar(lo, hi, 0.0, scalar, &mut equation, &mut cancelled)?
    };
    let accepted = accepted.ok_or(Error::Invalid("mixture-pair thermal endpoint missing"))?;
    if accepted.heat.to_bits() != heat.to_bits() {
        return Err(Error::Invalid(
            "mixture-pair thermal endpoint disagrees with scalar root",
        ));
    }
    let l = accepted.left_equilibrium;
    let r = accepted.right_equilibrium;
    let endpoint_temperature_k = [accepted.left.temperature(), accepted.right.temperature()];
    let flux = coefficient * (endpoint_temperature_k[1] - endpoint_temperature_k[0]);
    let heat_balance = l.heat_change_j + r.heat_change_j;
    let energy_residual = ((l.internal_energy_change_j - heat) + l.boundary_work_j)
        + ((r.internal_energy_change_j + heat) + r.boundary_work_j);
    let thermal_residual = l.heat_change_j - flux;
    let right_thermal_residual = -r.heat_change_j - flux;
    let scale = l
        .heat_change_j
        .abs()
        .max(r.heat_change_j.abs())
        .max(flux.abs());
    let tolerance = controls.thermal.tolerance(scale);
    if !heat_balance.is_finite() || heat_balance.abs() > tolerance {
        return Err(Error::Invalid(
            "mixture-pair original heats are not equal and opposite",
        ));
    }
    if !energy_residual.is_finite()
        || !thermal_residual.is_finite()
        || !right_thermal_residual.is_finite()
        || energy_residual.abs() > tolerance
        || thermal_residual.abs() > tolerance
        || right_thermal_residual.abs() > tolerance
    {
        return Err(Error::Balance {
            energy_residual_j: energy_residual,
            thermal_residual_j: if right_thermal_residual.abs() > thermal_residual.abs() {
                right_thermal_residual
            } else {
                thermal_residual
            },
        });
    }
    let entropy_production = l.entropy_change_j_k + r.entropy_change_j_k;
    if !entropy_production.is_finite() || entropy_production < -controls.entropy_tolerance_j_k {
        return Err(Error::Invalid(
            "mixture-pair original entropy balance refused",
        ));
    }
    let report = MixturePairThermalStep {
        elapsed: duration,
        heat_to_left_j: l.heat_change_j,
        heat_to_right_j: r.heat_change_j,
        internal_energy_change_j: [l.internal_energy_change_j, r.internal_energy_change_j],
        boundary_work_by_parcels_j: [l.boundary_work_j, r.boundary_work_j],
        energy_residual_j: energy_residual,
        thermal_residual_j: thermal_residual,
        entropy_change_j_k: [l.entropy_change_j_k, r.entropy_change_j_k],
        entropy_production_j_k: entropy_production,
        left_equilibrium: l.source,
        right_equilibrium: r.source,
        initial_temperature_k,
        endpoint_temperature_k,
        iterations,
        heat_evaluations: work.heat,
        pt_evaluations: work.pt,
        pressure_evaluations: work.pressure,
        eos_evaluation_upper_bound: work.eos,
    };
    // Both admitted source-bearing candidates are retained from the last scalar
    // query. There is no second flash or fallible operation after publication.
    poll(&mut cancelled)?;
    *left = accepted.left;
    *right = accepted.right;
    Ok(report)
}
