//! Reciprocal solid-vapor mass, recoil and heat on original finite stores.
//!
//! The original shifted-Maxwellian interface supplies complete species beams,
//! including reflected carrier pressure and heat. One frozen-rate interval
//! changes both actual bulk momenta and subtracts their actual kinetic changes
//! before reconstructing either original caloric store. Latent heat belongs to
//! the source enthalpies; it is not an additional heat debit or receiver loss.
//!
//! This is a translation-only exchange stage on a uniform pure-solid patch and
//! a homogeneous, fixed-volume gas. Patch area, normal and reference geometry
//! are held over the interval. No pose, angular momentum, moving-volume work,
//! Knudsen layer, spatial gas flow or whole-vessel pressure balance is inferred.
//! The caller must supply those owners when needed; this stage's momentum and
//! energy closure is not a whole-assembly passivity or entropy claim. Source
//! rates are first order in time, with explicit inventory/temperature limits.

use crate::liquid_gas::evaporation::momentum::{
    MovingGasParcel,
    boundary::{
        MovingGasBoundaryError, MovingGasBoundarySolve, MovingGasBoundaryStep, kinetic_change,
    },
};
use fs_material::sublimation::gas::{
    GasSurfaceError,
    chamber::GasChamberSurface,
    kinetic::{KineticInterfaceRate, MolecularInterfaceVelocity, evaluate_molecular_interface},
    momentum::{MaxwellReflection, MolecularSurfaceFrame},
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError, parcel::transfer::species::boundary::SpeciesBoundaryRelativeTransfer,
    },
};
use fs_math::dd::Dd;
use fs_thermal::sublimation::{
    SolidSurfaceState, SublimationError, SurfaceKinetics, VaporReservoir, surface_flux,
};

/// The original solid mass/enthalpy with independently retained bulk momentum.
/// This local source state has no pose, rotational inertia or implied volume.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MovingSolidSurface {
    solid: SolidSurfaceState,
    momentum_kg_m_s: [f64; 3],
}

impl MovingSolidSurface {
    pub fn new(
        solid: SolidSurfaceState,
        velocity_m_s: [f64; 3],
    ) -> Result<Self, KineticSublimationError> {
        let momentum = velocity_m_s.map(|v| solid.mass_kg * v);
        if !velocity_m_s.into_iter().all(f64::is_finite)
            || (0..3).any(|d| velocity_m_s[d] != 0.0 && momentum[d] == 0.0)
        {
            return Err(Error::Invalid(
                "finite representable original solid velocity required",
            ));
        }
        Self::from_momentum(solid, momentum)
    }

    /// Retain original represented momenta without a velocity round trip.
    pub fn from_momentum(
        solid: SolidSurfaceState,
        momentum_kg_m_s: [f64; 3],
    ) -> Result<Self, KineticSublimationError> {
        if !positive(solid.mass_kg)
            || !solid.specific_enthalpy_j_kg.is_finite()
            || !nonnegative(solid.area_m2)
            || !momentum_kg_m_s.into_iter().all(f64::is_finite)
        {
            return Err(Error::Invalid(
                "finite pure-solid storage and momentum required",
            ));
        }
        let state = Self {
            solid,
            momentum_kg_m_s,
        };
        if !state.velocity_m_s().into_iter().all(f64::is_finite)
            || !state.kinetic_energy_j().is_finite()
        {
            return Err(Error::Invalid(
                "solid velocity or kinetic storage is unrepresentable",
            ));
        }
        Ok(state)
    }

    pub const fn solid(&self) -> SolidSurfaceState {
        self.solid
    }
    pub const fn momentum_kg_m_s(&self) -> [f64; 3] {
        self.momentum_kg_m_s
    }
    pub fn velocity_m_s(&self) -> [f64; 3] {
        self.momentum_kg_m_s.map(|p| p / self.solid.mass_kg)
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.momentum_kg_m_s
            .into_iter()
            .zip(self.velocity_m_s())
            .fold(Dd::ZERO, |sum, (p, v)| {
                sum + Dd::from_f64(0.5 * p) * Dd::from_f64(v)
            })
            .to_f64()
    }
}

/// Original gas controls and physical limits for one frozen-rate interval.
#[derive(Clone, Copy, Debug)]
pub struct KineticSublimationControls {
    pub gas: MovingGasBoundarySolve,
    pub maximum_duration_s: f64,
    /// Both gross phase streams are bounded by their ORIGINAL donors, in (0,1).
    /// Simultaneous incoming mass cannot fund outgoing mass in this interval.
    pub maximum_gross_inventory_fraction: f64,
    pub maximum_temperature_change_k: f64,
    /// Absolute agreement of original vapor minus solid enthalpy with L(T).
    pub latent_compatibility_tolerance_j_kg: f64,
}

#[derive(Debug)]
pub enum KineticSublimationError {
    Invalid(&'static str),
    Surface(GasSurfaceError),
    Gas(GasEnergyError),
    MovingGas(MovingGasBoundaryError),
    InventoryLimit {
        maximum_duration_s: f64,
    },
    Duration {
        maximum_duration_s: f64,
    },
    TemperatureIncrement {
        maximum_duration_s: f64,
    },
    Balance {
        account: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}
impl core::fmt::Display for KineticSublimationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite solid/gas recoil refused: {self:?}")
    }
}
impl std::error::Error for KineticSublimationError {}
impl From<GasEnergyError> for KineticSublimationError {
    fn from(value: GasEnergyError) -> Self {
        if value == GasEnergyError::Cancelled {
            Self::Cancelled
        } else {
            Self::Gas(value)
        }
    }
}
impl From<GasSurfaceError> for KineticSublimationError {
    fn from(value: GasSurfaceError) -> Self {
        if matches!(value, GasSurfaceError::Surface(SublimationError::Cancelled)) {
            Self::Cancelled
        } else {
            Self::Surface(value)
        }
    }
}
impl From<SublimationError> for KineticSublimationError {
    fn from(value: SublimationError) -> Self {
        GasSurfaceError::from(value).into()
    }
}
impl From<MovingGasBoundaryError> for KineticSublimationError {
    fn from(value: MovingGasBoundaryError) -> Self {
        if matches!(
            value,
            MovingGasBoundaryError::Cancelled
                | MovingGasBoundaryError::Gas(GasEnergyError::Cancelled)
        ) {
            Self::Cancelled
        } else {
            Self::MovingGas(value)
        }
    }
}

/// Source-relative endpoint balances in coherent SI units.
/// Thermal and kinetic increments use actual represented states. Formation
/// reference transport uses requested dm, while represented mass roundoff is
/// independently bounded by the original total and species mass tolerances.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticSublimationStep {
    pub duration_s: f64,
    /// Original immutable molecular rates evaluated once at the initial state.
    pub rate: KineticInterfaceRate,
    pub outgoing_mass_kg: f64,
    pub incoming_mass_kg: f64,
    pub transferred_mass_kg: f64,
    pub supplied_heat_j: f64,
    pub surface_temperature_k: [f64; 2],
    pub gas_temperature_k: [f64; 2],
    /// Actual remaining solid sensible storage, m1*(h1-h0).
    pub solid_thermal_increment_j: f64,
    /// Requested transferred-mass reference energy plus actual thermal storage.
    /// Reporting only near large datums; represented mass rounding is bounded
    /// separately and is not multiplied into this source-relative energy test.
    pub solid_energy_change_j: f64,
    /// Actual solid kinetic change in the initial-solid-velocity frame.
    pub solid_kinetic_energy_change_j: f64,
    pub solid_laboratory_kinetic_energy_change_j: f64,
    pub solid_momentum_change_n_s: [f64; 3],
    pub gas: MovingGasBoundaryStep,
    pub stream_mass_rounding_kg: f64,
    pub mass_residual_kg: f64,
    pub maximum_species_mass_residual_kg: f64,
    pub momentum_residual_n_s: [f64; 3],
    /// Joint internal + bulk kinetic closure in the translating energy frame.
    pub energy_residual_j: f64,
    /// Same account in the retained momenta's laboratory frame.
    pub laboratory_energy_residual_j: f64,
    /// Positive retreat, negative growth; geometry itself remains caller-owned.
    pub recession_m: f64,
}

/// Both endpoints publish together; errors never mutate an accepted owner.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticSublimationCandidate {
    solid: MovingSolidSurface,
    gas: MovingGasParcel,
    step: KineticSublimationStep,
}
impl KineticSublimationCandidate {
    pub const fn solid(&self) -> &MovingSolidSurface {
        &self.solid
    }
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn step(&self) -> &KineticSublimationStep {
        &self.step
    }
    pub fn into_parts(self) -> (MovingSolidSurface, MovingGasParcel, KineticSublimationStep) {
        (self.solid, self.gas, self.step)
    }
}

type Error = KineticSublimationError;
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: Dd) -> Result<f64, Error> {
    if x.hi.is_finite() && x.lo.is_finite() && x.to_f64().is_finite() {
        Ok(x.to_f64())
    } else {
        Err(Error::Invalid(
            "molecular mass, energy or momentum is unrepresentable",
        ))
    }
}
fn checked_balance(account: &'static str, residual: Dd, allowance: f64) -> Result<f64, Error> {
    let residual = finite(residual)?;
    if !positive(allowance) || residual.abs() > allowance {
        Err(Error::Balance {
            account,
            residual,
            allowance,
        })
    } else {
        Ok(residual)
    }
}
fn admit_enthalpy(h: f64, range: [f64; 2]) -> Result<(), Error> {
    if !range.into_iter().all(f64::is_finite)
        || range[0] >= range[1]
        || !h.is_finite()
        || h < range[0]
        || h > range[1]
    {
        Err(Error::Invalid(
            "enthalpy lies outside the declared original pure-solid branch",
        ))
    } else {
        Ok(())
    }
}
fn inventory_bound(dt: f64, gross: f64, mass: f64, fraction: f64) -> Option<f64> {
    if gross == 0.0 || (mass > 0.0 && gross / mass <= fraction) {
        None
    } else {
        Some((dt * (mass * fraction / gross)).next_down().max(0.0))
    }
}

/// Advance the original pure solid and finite moving gas under reciprocal beams.
///
/// The source's `solid` must equal the supplied motion state's original store.
/// `pure_solid_enthalpy_range_j_kg` is an explicit inclusive branch of that
/// ORIGINAL chart; a melting plateau must not be admitted as a pure solid.
/// `frame` uses the momenta's inertial axes and represents a translation-only
/// resultant through the declared patch centroid. Angular/geometry evolution
/// requires a separate owner and is not silently inferred from this state.
///
/// External `supplied_heat_j` enters only the original solid thermal store.
/// All phase/reflected stream energy is integrated in the inertial frame moving
/// at the INITIAL solid velocity. Actual kinetic changes in that frame are
/// removed before caloric inversion, preserving a small heat balance beside a
/// large common drift. Both laboratory and translated balances must pass.
/// Frozen rates require physical time refinement, even when energy closes.
#[allow(clippy::too_many_arguments)]
pub fn advance_solid_with_moving_gas(
    surface: &GasChamberSurface<'_>,
    solid: &MovingSolidSurface,
    gas: &MovingGasParcel,
    catalog: &ThermochemicalGasModel,
    frame: MolecularSurfaceFrame,
    pure_solid_enthalpy_range_j_kg: [f64; 2],
    scattering: Option<&MaxwellReflection>,
    duration_s: f64,
    supplied_heat_j: f64,
    controls: KineticSublimationControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<KineticSublimationCandidate, Error> {
    poll(&mut cancelled)?;
    if surface.solid != solid.solid()
        || !nonnegative(duration_s)
        || !supplied_heat_j.is_finite()
        || !positive(controls.maximum_duration_s)
        || !positive(controls.maximum_temperature_change_k)
        || !positive(controls.maximum_gross_inventory_fraction)
        || controls.maximum_gross_inventory_fraction >= 1.0
        || !nonnegative(controls.latent_compatibility_tolerance_j_kg)
        || !positive(controls.gas.gas.mass_tolerance_kg)
        || !positive(controls.gas.gas.energy_tolerance_j)
        || !positive(controls.gas.momentum_absolute_n_s)
        || !nonnegative(controls.gas.momentum_relative_tolerance)
        || controls.gas.momentum_relative_tolerance >= 1.0
    {
        return Err(Error::Invalid(
            "coherent original storage, duration and positive budgets required",
        ));
    }
    if duration_s > controls.maximum_duration_s {
        return Err(Error::Duration {
            maximum_duration_s: controls.maximum_duration_s,
        });
    }
    let old = solid.solid();
    admit_enthalpy(old.specific_enthalpy_j_kg, pure_solid_enthalpy_range_j_kg)?;
    let energy_frame = solid.velocity_m_s();
    let rate = evaluate_molecular_interface(
        surface,
        None,
        frame,
        frame.centroid_m,
        MolecularInterfaceVelocity {
            gas_velocity_m_s: gas.velocity_m_s(),
            surface_velocity_m_s: energy_frame,
        },
        scattering,
        gas.parcel(),
        catalog,
        controls.latent_compatibility_tolerance_j_kg,
        &mut cancelled,
    )?;
    let exposure = finite(Dd::from_f64(old.area_m2) * Dd::from_f64(duration_s))?;
    if old.area_m2 != 0.0 && duration_s != 0.0 && exposure == 0.0 {
        return Err(Error::Invalid("surface exposure underflowed"));
    }
    let before_species = gas.parcel().species_inventory_kg()?;
    let mut transfers = Vec::with_capacity(rate.species.len());
    let (mut outgoing, mut incoming, mut mass, mut relative_energy) =
        (Dd::ZERO, Dd::ZERO, Dd::ZERO, Dd::ZERO);
    let mut relative_impulse = [Dd::ZERO; 3];
    let mut stream_rounding = Dd::ZERO;
    let mut maximum_duration: Option<f64> = None;
    for row in &rate.species {
        poll(&mut cancelled)?;
        let emitted = finite(Dd::from_f64(row.emitted_mass_flux_kg_m2_s) * Dd::from_f64(exposure))?;
        let captured =
            finite(Dd::from_f64(row.sticking_mass_flux_kg_m2_s) * Dd::from_f64(exposure))?;
        let dm = finite(Dd::from_f64(row.net_mass_flux_into_gas_kg_m2_s) * Dd::from_f64(exposure))?;
        let heat = finite(
            Dd::from_f64(row.relative_energy_flux_beyond_gas_enthalpy_w_m2)
                * Dd::from_f64(exposure),
        )?;
        let source_mass = before_species.get(&row.species).copied().unwrap_or(0.0);
        for limit in [
            inventory_bound(
                duration_s,
                emitted,
                old.mass_kg,
                controls.maximum_gross_inventory_fraction,
            ),
            inventory_bound(
                duration_s,
                captured,
                source_mass,
                controls.maximum_gross_inventory_fraction,
            ),
        ]
        .into_iter()
        .flatten()
        {
            maximum_duration = Some(maximum_duration.map_or(limit, |old| old.min(limit)));
        }
        let rounding = (Dd::from_f64(emitted) - Dd::from_f64(captured)) - Dd::from_f64(dm);
        checked_balance(
            "molecular gross/net mass",
            rounding,
            controls.gas.gas.mass_tolerance_kg,
        )?;
        stream_rounding = stream_rounding + rounding;
        outgoing = outgoing + Dd::from_f64(emitted);
        incoming = incoming + Dd::from_f64(captured);
        mass = mass + Dd::from_f64(dm);
        relative_energy = relative_energy + Dd::from_f64(heat);
        for d in 0..3 {
            relative_impulse[d] = relative_impulse[d]
                + Dd::from_f64(row.relative_momentum_flux_into_gas_pa[d]) * Dd::from_f64(exposure);
        }
        transfers.push(SpeciesBoundaryRelativeTransfer {
            species: row.species.clone(),
            mass_into_gas_kg: dm,
            relative_energy_into_gas_j: heat,
        });
    }
    if let Some(maximum_duration_s) = maximum_duration {
        return Err(Error::InventoryLimit { maximum_duration_s });
    }
    let dm = finite(mass)?;
    let next_mass = old.mass_kg - dm;
    if !positive(next_mass) || (dm != 0.0 && next_mass == old.mass_kg) {
        return Err(Error::Invalid(
            "solid exhaustion or unrepresentable original mass change",
        ));
    }
    let solid_mass_change = Dd::from_f64(next_mass) - Dd::from_f64(old.mass_kg);
    checked_balance(
        "solid mass",
        solid_mass_change + mass,
        controls.gas.gas.mass_tolerance_kg,
    )?;
    let mut impulse = [0.0; 3];
    let mut next_momentum = solid.momentum_kg_m_s();
    let mut solid_momentum_change = [0.0; 3];
    for d in 0..3 {
        impulse[d] = finite(relative_impulse[d] + Dd::from_f64(energy_frame[d]) * mass)?;
        next_momentum[d] -= impulse[d];
        let change = Dd::from_f64(next_momentum[d]) - Dd::from_f64(solid.momentum_kg_m_s()[d]);
        solid_momentum_change[d] = finite(change)?;
        let allowance = controls.gas.momentum_absolute_n_s
            + controls.gas.momentum_relative_tolerance * impulse[d].abs();
        checked_balance(
            "solid impulse",
            change + Dd::from_f64(impulse[d]),
            allowance,
        )?;
        if impulse[d] != 0.0 && solid_momentum_change[d] == 0.0 {
            return Err(Error::Invalid(
                "molecular impulse vanished at original solid momentum",
            ));
        }
    }
    let solid_kinetic = kinetic_change(
        old.mass_kg,
        next_mass,
        solid.momentum_kg_m_s(),
        next_momentum,
        energy_frame,
    )?;
    let solid_lab_kinetic = kinetic_change(
        old.mass_kg,
        next_mass,
        solid.momentum_kg_m_s(),
        next_momentum,
        [0.0; 3],
    )?;
    let old_phase = surface.law.properties(rate.surface_temperature_k)?;
    let delta_h = catalog.species_specific_enthalpy_change_j_kg(
        &surface.binding.species,
        rate.gas_temperature_k,
        rate.surface_temperature_k,
    )?;
    // h_s+offset-h_v(T_gas) = Delta h_v(T_gas,T_s)-L-source_alignment.
    // Cancel formation references BEFORE multiplying by the transported mass.
    let thermal_target = mass
        * (Dd::from_f64(delta_h)
            - Dd::from_f64(old_phase.latent_heat_j_kg)
            - Dd::from_f64(rate.source_alignment_residual_j_kg))
        - relative_energy
        - Dd::from_f64(solid_kinetic)
        + Dd::from_f64(supplied_heat_j);
    let dh = finite(thermal_target / Dd::from_f64(next_mass))?;
    let next_h = old.specific_enthalpy_j_kg + dh;
    if dh != 0.0 && next_h == old.specific_enthalpy_j_kg {
        return Err(Error::Invalid(
            "solid caloric increment vanished in its original chart",
        ));
    }
    admit_enthalpy(next_h, pure_solid_enthalpy_range_j_kg)?;
    let next_temperature = (surface.temperature_from_enthalpy)(next_h)?;
    // Reuse the original thermal admission (including phase range and finite
    // properties) without another molecular species loop or gas inversion.
    surface_flux(
        surface.law,
        next_temperature,
        SurfaceKinetics {
            accommodation_coefficient: 0.0,
            mass_transfer_coefficient_m_s: None,
        },
        VaporReservoir::Prescribed {
            partial_pressure_pa: 0.0,
            temperature_k: next_temperature,
        },
    )?;
    let next_phase = surface.law.properties(next_temperature)?;
    if (next_phase.molar_mass_kg_mol / old_phase.molar_mass_kg_mol - 1.0).abs()
        > 64.0 * f64::EPSILON
    {
        return Err(Error::Invalid(
            "phase law changed its original molecular identity",
        ));
    }
    let vapor_h =
        catalog.species_specific_enthalpy_j_kg(&surface.binding.species, next_temperature)?;
    let alignment = (vapor_h - (next_h + surface.binding.solid_enthalpy_offset_j_kg))
        - next_phase.latent_heat_j_kg;
    if !alignment.is_finite() || alignment.abs() > controls.latent_compatibility_tolerance_j_kg {
        return Err(GasSurfaceError::CaloricMismatch {
            temperature_k: next_temperature,
            residual_j_kg: alignment,
            tolerance_j_kg: controls.latent_compatibility_tolerance_j_kg,
        }
        .into());
    }
    if (next_temperature - rate.surface_temperature_k).abs() > controls.maximum_temperature_change_k
    {
        return Err(Error::TemperatureIncrement {
            maximum_duration_s: 0.5 * duration_s,
        });
    }
    let thermal_change =
        Dd::from_f64(next_mass) * (Dd::from_f64(next_h) - Dd::from_f64(old.specific_enthalpy_j_kg));
    let solid_energy_residual = checked_balance(
        "solid caloric energy",
        thermal_change - thermal_target,
        controls.gas.gas.energy_tolerance_j,
    )?;
    let next_solid = MovingSolidSurface::from_momentum(
        SolidSurfaceState {
            mass_kg: next_mass,
            specific_enthalpy_j_kg: next_h,
            area_m2: old.area_m2,
        },
        next_momentum,
    )?;
    poll(&mut cancelled)?;
    let gas_candidate = gas.exchange_species_momentum_relative(
        catalog,
        &transfers,
        impulse,
        energy_frame,
        controls.gas,
        &mut cancelled,
    )?;
    let (next_gas, _, gas_step) = gas_candidate.into_parts();
    let next_gas_temperature = next_gas.parcel().state().state().temperature;
    if (next_gas_temperature - rate.gas_temperature_k).abs() > controls.maximum_temperature_change_k
    {
        return Err(Error::TemperatureIncrement {
            maximum_duration_s: 0.5 * duration_s,
        });
    }
    let gas_mass_change =
        Dd::from_f64(next_gas.parcel().mass_kg()) - Dd::from_f64(gas.parcel().mass_kg());
    let total_mass = solid_mass_change + gas_mass_change;
    let mass_residual = checked_balance(
        "joint solid/gas mass",
        total_mass,
        controls.gas.gas.mass_tolerance_kg,
    )?;
    let after_species = next_gas.parcel().species_inventory_kg()?;
    let mut max_species_residual: f64 = 0.0;
    for species in catalog.thermodynamics().components() {
        poll(&mut cancelled)?;
        let id = species.species();
        let mut residual = Dd::from_f64(after_species.get(id).copied().unwrap_or(0.0))
            - Dd::from_f64(before_species.get(id).copied().unwrap_or(0.0));
        if id == &surface.binding.species {
            residual = residual + solid_mass_change;
        }
        let value = checked_balance(
            "original individual species mass",
            residual,
            controls.gas.gas.mass_tolerance_kg,
        )?;
        max_species_residual = max_species_residual.max(value.abs());
    }
    let mut momentum_residual = [0.0; 3];
    for d in 0..3 {
        let change = (Dd::from_f64(next_momentum[d]) - Dd::from_f64(solid.momentum_kg_m_s()[d]))
            + (Dd::from_f64(next_gas.momentum_kg_m_s()[d])
                - Dd::from_f64(gas.momentum_kg_m_s()[d]));
        momentum_residual[d] = checked_balance(
            "joint momentum",
            change,
            controls.gas.momentum_absolute_n_s
                + controls.gas.momentum_relative_tolerance * impulse[d].abs(),
        )?;
    }
    let energy_residual = checked_balance(
        "joint translated total energy",
        Dd::from_f64(solid_energy_residual) + Dd::from_f64(gas_step.total_energy_residual_j),
        controls.gas.gas.energy_tolerance_j,
    )?;
    let mut laboratory_residual = Dd::from_f64(energy_residual);
    for d in 0..3 {
        laboratory_residual = laboratory_residual
            + Dd::from_f64(energy_frame[d]) * Dd::from_f64(momentum_residual[d])
            - (Dd::from_f64(0.5 * energy_frame[d]) * Dd::from_f64(energy_frame[d])) * total_mass;
    }
    let laboratory_energy_residual = checked_balance(
        "joint laboratory total energy",
        laboratory_residual,
        controls.gas.gas.energy_tolerance_j,
    )?;
    let recession = if old.area_m2 == 0.0 {
        0.0
    } else {
        dm / next_phase.solid_density_kg_m3 / old.area_m2
    };
    if !recession.is_finite() {
        return Err(Error::Invalid("surface recession is unrepresentable"));
    }
    let solid_energy_change = thermal_change
        - mass
            * (Dd::from_f64(old.specific_enthalpy_j_kg)
                + Dd::from_f64(surface.binding.solid_enthalpy_offset_j_kg));
    let step = KineticSublimationStep {
        duration_s,
        outgoing_mass_kg: finite(outgoing)?,
        incoming_mass_kg: finite(incoming)?,
        transferred_mass_kg: dm,
        supplied_heat_j,
        surface_temperature_k: [rate.surface_temperature_k, next_temperature],
        gas_temperature_k: [rate.gas_temperature_k, next_gas_temperature],
        solid_thermal_increment_j: finite(thermal_change)?,
        solid_energy_change_j: finite(solid_energy_change)?,
        solid_kinetic_energy_change_j: solid_kinetic,
        solid_laboratory_kinetic_energy_change_j: solid_lab_kinetic,
        solid_momentum_change_n_s: solid_momentum_change,
        gas: gas_step,
        rate,
        stream_mass_rounding_kg: finite(stream_rounding)?,
        mass_residual_kg: mass_residual,
        maximum_species_mass_residual_kg: max_species_residual,
        momentum_residual_n_s: momentum_residual,
        energy_residual_j: energy_residual,
        laboratory_energy_residual_j: laboratory_energy_residual,
        recession_m: recession,
    };
    poll(&mut cancelled)?;
    Ok(KineticSublimationCandidate {
        solid: next_solid,
        gas: next_gas,
        step,
    })
}
