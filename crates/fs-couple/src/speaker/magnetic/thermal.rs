//! Material voice-coil heating and power compression on the native enthalpy owner.
//!
//! Electrical power dissipated in R(T) heats the finite coil; the refreshed
//! resistance changes current, force and sound on the SAME accepted step. This
//! is the thermal feedback described by Klippel's voice-coil measurements:
//! <https://www.klippel.de/know-how/measurements/power-handling-heat-dissipation/voice-coil-temperature.html>.
//! No copper coefficient, heat capacity, cooling time or temperature-dependent
//! magnet strength is invented. Caloric, resistivity/conductivity and optional
//! thermal-link data come from the original MaterialNetwork and its cards.
//!
//! The coil is one uniform, wholly solid, fixed-density cell. Other network
//! cells, authored reservoirs and sourced links may represent the former,
//! magnet and enclosure. This does not infer motion-dependent convection or
//! skin/eddy losses. Magnetic profiles and mechanical properties remain fixed.
//! Suspension, compact radiation and panel losses leave this represented
//! system through an explicit mechanical/acoustic sink; they are NOT coil heat.

use super::{MagneticTsDriver, PanelMode, SpeakerError, TsCard, VoiceCoilMagnetics};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients,
    DriverAccount, ThermalDriver, cell_driver_candidate,
    material::{
        MaterialNetwork, MaterialNetworkError, MaterialNetworkSample, PreparedMaterialNetwork,
    },
};
use fs_material::conductor::{ElectricalPropertyBasis, ResolvedConductor};
use fs_material::gas::GasState;
use fs_phs::{StepRecord, discrete_gradient};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap, sync::Arc};

/// Original material or magnetic speaker refusal, without fallback constants.
#[derive(Debug)]
pub enum SpeakerThermalError {
    Material(MaterialNetworkError),
    Speaker(SpeakerError),
}
impl core::fmt::Display for SpeakerThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "magnetic speaker/thermal source: {self:?}")
    }
}
impl std::error::Error for SpeakerThermalError {}
type Error = CircuitEnthalpyError<SpeakerThermalError>;

fn material(error: MaterialNetworkError) -> Error {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => {
            CircuitEnthalpyError::Cancelled
        }
        other => CircuitEnthalpyError::Source(SpeakerThermalError::Material(other)),
    }
}
fn speaker(error: SpeakerError) -> Error {
    CircuitEnthalpyError::Source(SpeakerThermalError::Speaker(error))
}
fn invalid(what: &'static str) -> Error {
    CircuitEnthalpyError::Invalid(what)
}

/// Mechanical inputs independent of the source-defined electrical coil.
///
/// The non-coil mass includes the cone, former and other moving parts, but
/// excludes both the coil and the native acoustic added mass. Suspension,
/// damping and piston geometry remain explicit model inputs; coil composition
/// does not imply any of them or a permanent-magnet constitutive law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeakerMechanics {
    /// Finite nonnegative moving mass excluding the whole electrical coil [kg].
    pub non_coil_moving_mass_kg: f64,
    /// Suspension compliance [m/N].
    pub cms_m_per_n: f64,
    /// Mechanical loss [N s/m], not electrical heating.
    pub rms_n_s_m: f64,
    /// Effective piston area [m²].
    pub sd_m2: f64,
    /// Linear suspension excursion limit [m].
    pub x_max_m: f64,
}

// Both new-specimen construction and continuing-state binding admit the same
// whole-cell geometry, original source pins and frozen solid mass inventory.
fn prepare_coil<'n, 'a>(
    thermal: &DepositionSimulation<'_>,
    network: &'n MaterialNetwork<'a>,
    basis: ElectricalPropertyBasis,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<
    (
        PreparedMaterialNetwork<'n, 'a>,
        CellId,
        ResolvedConductor,
        f64,
    ),
    Error,
> {
    if cancelled() {
        return Err(CircuitEnthalpyError::Cancelled);
    }
    if network.resistors.len() != 1 || network.resistors[0].destination.branch != 0 {
        return Err(invalid("one whole-cell coil binding at branch 0 required"));
    }
    let prepared = network
        .prepare_with_bases(thermal, &[(0, basis)], &mut *cancelled)
        .map_err(material)?;
    let coil_cell = network.resistors[0].destination.cell;
    let cell = &thermal.cells()[&coil_cell];
    let binding = network
        .cells
        .iter()
        .find(|b| b.cell == coil_cell)
        .ok_or_else(|| invalid("coil lacks its original caloric binding"))?;
    let density = cell.mass_kg / cell.geometry.volume_m3;
    if !density.is_finite()
        || density <= 0.0
        || binding.curve.knots().iter().any(|k| {
            k.liquid_mass_fraction != 0.0
                || (k.bulk_density_kg_m3 - density).abs() > 1e-10 * density
        })
    {
        return Err(invalid(
            "coil requires fixed solid density and unchanged matter",
        ));
    }
    let mut initial = prepared
        .resolve(thermal, &mut *cancelled)
        .map_err(material)?;
    let coil = initial.resistors.remove(0).conductor;
    if cancelled() {
        return Err(CircuitEnthalpyError::Cancelled);
    }
    Ok((prepared, coil_cell, coil, cell.mass_kg))
}

/// One native electrical/mechanical work record with physically separated losses.
#[derive(Clone, Debug)]
pub struct MagneticSpeakerThermalStep {
    /// Original descriptor step, including ALL electrical/mechanical dissipation.
    pub native: StepRecord,
    /// Same compact 1-m pressure observer as the original magnetic driver [Pa].
    pub pressure_pa: f64,
    /// Original source and geometry at the accepted coil temperature.
    pub coil: ResolvedConductor,
    /// Coil I^2 R dt [J], deposited exactly once into its receiving cell.
    pub coil_heat_j: f64,
    /// Suspension, compact radiation and panel dissipation [J]. This explicit
    /// external sink is excluded from thermal cell heating.
    pub mechanical_acoustic_loss_j: f64,
}

/// Magnetic speaker retaining one immutable, prepared material-network binding.
///
/// All source selections stay pinned across accepted steps and Newton trials.
/// Electrical geometry and caloric storage describe the SAME whole coil cell.
/// The ordinary MagneticTsDriver remains available for fixed-temperature work.
pub struct ThermalMagneticTsDriver<'n, 'a> {
    driver: MagneticTsDriver,
    prepared: Arc<PreparedMaterialNetwork<'n, 'a>>,
    coil_cell: CellId,
    coil: ResolvedConductor,
    coil_mass_kg: f64,
}

impl<'n, 'a> ThermalMagneticTsDriver<'n, 'a> {
    /// Construct a material-defined coil and a new mechanical state at rest.
    ///
    /// Initial Re and coil inertia come from the SAME prepared material network,
    /// caloric cell and whole-cell conductor geometry. Total dry moving mass is
    /// that original cell mass plus `non_coil_moving_mass_kg`; the native driver
    /// separately adds the supplied gas's acoustic load. Bl and Le at zero are
    /// taken directly from the original magnetic profile, with no duplicated
    /// card parameters to become stale under substitution.
    ///
    /// Motion starts at rest at the supplied thermal clock. This creates a new
    /// specimen, not a momentum-conserving replacement of a moving coil. Use
    /// [`Self::step_candidate`] to continue its physical history, or [`Self::bind`]
    /// to attach an already-running native driver without resetting its state.
    /// The original solid-density, source-domain and excursion limits apply;
    /// heating updates R(T), not coil geometry, mass or magnetic profiles.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material(
        mechanics: SpeakerMechanics,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        gas: GasState,
        profile: VoiceCoilMagnetics,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        basis: ElectricalPropertyBasis,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::from_material_with_media(
            mechanics,
            box_volume_m3,
            panel,
            gas,
            box_volume_m3.map(|_| gas),
            profile,
            thermal,
            network,
            basis,
            cancelled,
        )
    }

    /// Resolve the same whole material coil with distinct exterior/cavity media.
    ///
    /// Gas-region roles and static-pressure admission are exactly those of
    /// `MagneticTsDriver::new_with_media`. Coil heating retains both media and
    /// their original mechanical operators; it does not silently heat either
    /// gas or replace the sealed fill with the exterior. Any finite heat exchange
    /// remains an explicit original thermal-network binding.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_with_media(
        mechanics: SpeakerMechanics,
        box_volume_m3: Option<f64>,
        panel: Option<PanelMode>,
        exterior_gas: GasState,
        cavity_gas: Option<GasState>,
        profile: VoiceCoilMagnetics,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        basis: ElectricalPropertyBasis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        if !mechanics.non_coil_moving_mass_kg.is_finite() || mechanics.non_coil_moving_mass_kg < 0.0
        {
            return Err(invalid(
                "non-coil moving mass must be finite and nonnegative",
            ));
        }
        let (prepared, coil_cell, coil, coil_mass_kg) =
            prepare_coil(thermal, network, basis, &mut cancelled)?;
        let origin = profile.sample(0.0).map_err(speaker)?;
        let card = TsCard {
            re_ohm: coil.resistance_ohm(),
            le_h: origin.inductance_h,
            bl: origin.force_factor_tm,
            mms_kg: mechanics.non_coil_moving_mass_kg + coil_mass_kg,
            cms_m_per_n: mechanics.cms_m_per_n,
            rms_n_s_m: mechanics.rms_n_s_m,
            sd_m2: mechanics.sd_m2,
            x_max_m: mechanics.x_max_m,
        };
        let mut driver = MagneticTsDriver::new_with_media(
            card,
            box_volume_m3,
            panel,
            exterior_gas,
            cavity_gas,
            profile,
        )
        .map_err(speaker)?;
        driver.time_s = thermal.time_s();
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(Self {
            driver,
            prepared: Arc::new(prepared),
            coil_cell,
            coil,
            coil_mass_kg,
        })
    }

    /// Bind a driver to actual thermal/material history at its current clock.
    ///
    /// `network` must contain exactly one CellResistor, with branch 0 denoting
    /// this driver's sole electrical coil coordinate. Its length times area
    /// must cover the whole receiving cell, as checked by MaterialNetwork.
    /// Other cells are passive heat recipients, not additional electrical coils.
    /// The sourced initial resistance must agree with the accepted driver's Re.
    /// The coil chart must be wholly solid and constant-density; coil mass may
    /// not exceed the card's cone-plus-coil moving mass. No melting, expansion
    /// or undeclared mass change is silently applied to fixed speaker geometry.
    pub fn bind(
        driver: MagneticTsDriver,
        thermal: &DepositionSimulation<'_>,
        network: &'n MaterialNetwork<'a>,
        basis: ElectricalPropertyBasis,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        if driver.time_s().to_bits() != thermal.time_s().to_bits() {
            return Err(invalid("matching speaker/thermal clocks required"));
        }
        let (prepared, coil_cell, coil, coil_mass_kg) =
            prepare_coil(thermal, network, basis, &mut cancelled)?;
        if coil_mass_kg > driver.card.mms_kg * (1.0 + 1e-10) {
            return Err(invalid(
                "coil mass must fit within the original moving assembly",
            ));
        }
        if (coil.resistance_ohm() - driver.coil_resistance_ohm()).abs()
            > 1e-12 * driver.coil_resistance_ohm()
        {
            return Err(invalid(
                "initial sourced coil resistance must match the accepted magnetic driver",
            ));
        }
        if cancelled() {
            return Err(CircuitEnthalpyError::Cancelled);
        }
        Ok(Self {
            driver,
            prepared: Arc::new(prepared),
            coil_cell,
            coil,
            coil_mass_kg,
        })
    }

    /// Accepted original magnetic/electromechanical state.
    pub const fn driver(&self) -> &MagneticTsDriver {
        &self.driver
    }

    /// Accepted source-bound conductor at the actual coil temperature.
    pub const fn coil(&self) -> &ResolvedConductor {
        &self.coil
    }

    /// Actual receiving thermal cell; only electrical coil loss is deposited.
    pub const fn coil_cell(&self) -> CellId {
        self.coil_cell
    }

    /// Original whole-cell coil mass, retained through thermal continuation [kg].
    pub const fn coil_mass_kg(&self) -> f64 {
        self.coil_mass_kg
    }

    /// Coil plus explicitly authored other moving mass, excluding gas load [kg].
    pub const fn moving_assembly_mass_kg(&self) -> f64 {
        self.driver.card.mms_kg
    }

    /// Advance voltage, magnetic motion, R(T), finite heat storage and cooling
    /// jointly through the EXISTING cell enthalpy/Newton solver.
    ///
    /// Every trial starts from the same accepted linkage and mechanical state.
    /// Coil loss uses the native discrete-gradient current and endpoint R(T).
    /// Thermal backward Euler makes the general coupled scheme first order;
    /// resolve accuracy by timestep refinement. No coefficient-lag fallback is
    /// permitted. Finite-difference stencils and physical endpoints must stay
    /// inside the original sources and the admitted magnetic excursion range.
    ///
    /// Inputs remain unchanged on source refusal, excursion failure, exhausted
    /// work, cancellation or energy failure. Process start events first and
    /// split intervals at interior thermal/source events. Publish both returned
    /// histories together; the returned thermal state already includes Joule heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        voltage_v: f64,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MagneticSpeakerEnthalpyCandidate<'t, 'n, 'a>, Error> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(CircuitEnthalpyError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        if self.driver.time_s().to_bits() != thermal.time_s().to_bits()
            || thermal
                .cells()
                .get(&self.coil_cell)
                .is_none_or(|cell| cell.mass_kg.to_bits() != self.coil_mass_kg.to_bits())
            || !target_s.is_finite()
            || target_s <= thermal.time_s()
            || !voltage_v.is_finite()
            || [
                controls.circuit.energy_tolerance_j,
                controls.circuit.power_tolerance_w,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "matching clocks and original coil mass, positive interval and finite voltage/work controls required",
            ));
        }
        let initial = self
            .prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        if initial.resistors[0].conductor != self.coil {
            return Err(invalid(
                "accepted coil source must match the actual thermal endpoint",
            ));
        }
        let dt = target_s - thermal.time_s();
        let (thermal, driver, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |state, _remaining| {
                check()?;
                let mut coefficients = self
                    .prepared
                    .resolve(state, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let coil = coefficients.resistors.remove(0).conductor;
                let mut driver = self
                    .driver
                    .with_coil_resistance(coil.resistance_ohm())
                    .map_err(speaker)?;
                let (native, pressure_pa) = driver.step(voltage_v, dt).map_err(speaker)?;
                driver.time_s = target_s;
                // Original diagonal dissipation, evaluated at the SAME effort
                // used by the native step. This is a partition, not a second
                // force law or an estimate using endpoint physical current.
                let effort = discrete_gradient(driver.phs.storage(), &self.driver.x, &native.x);
                let r = driver.phs.structure().1;
                let n = effort.len();
                let coil_power = r[0] * effort[0] * effort[0];
                let mechanical_power = (1..n)
                    .map(|i| r[i * n + i] * effort[i] * effort[i])
                    .sum::<f64>();
                let coil_heat_j = dt * coil_power;
                let mechanical_acoustic_loss_j = dt * mechanical_power;
                if !coil_heat_j.is_finite()
                    || !mechanical_acoustic_loss_j.is_finite()
                    || coil_heat_j < 0.0
                    || mechanical_acoustic_loss_j < 0.0
                    || native.balance_residual().abs() > controls.circuit.energy_tolerance_j
                    || (native.dissipated - coil_heat_j - mechanical_acoustic_loss_j).abs()
                        > controls.circuit.energy_tolerance_j
                    || (native.dissipated / dt - coil_power - mechanical_power).abs()
                        > controls.circuit.power_tolerance_w
                {
                    return Err(invalid(
                        "native speaker work and separate coil/mechanical losses do not close",
                    ));
                }
                check()?;
                let step = MagneticSpeakerThermalStep {
                    native,
                    pressure_pa,
                    coil: coil.clone(),
                    coil_heat_j,
                    mechanical_acoustic_loss_j,
                };
                let next = Self {
                    driver,
                    prepared: Arc::clone(&self.prepared),
                    coil_cell: self.coil_cell,
                    coil,
                    coil_mass_kg: self.coil_mass_kg,
                };
                Ok(DrivenCoefficients {
                    circuit: SpeakerCandidate { driver: next, step },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: coefficients.source,
                    circuit_trials: 1,
                })
            },
            |candidate, dt, checkpoint| {
                checkpoint()?;
                let power = candidate.step.coil_heat_j / dt;
                if !power.is_finite() {
                    return Err(invalid("coil heat power is unrepresentable"));
                }
                Ok(BTreeMap::from([(self.coil_cell, power)]))
            },
            || (*cancel.borrow_mut())(),
        )?;
        check()?;
        Ok(MagneticSpeakerEnthalpyCandidate {
            thermal,
            driver,
            report,
        })
    }
}

struct SpeakerCandidate<'n, 'a> {
    driver: ThermalMagneticTsDriver<'n, 'a>,
    step: MagneticSpeakerThermalStep,
}
impl<'n, 'a> ThermalDriver for SpeakerCandidate<'n, 'a> {
    type Model = ThermalMagneticTsDriver<'n, 'a>;
    type Step = MagneticSpeakerThermalStep;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.step.coil_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native.delta_h,
            // The original driver emits this mechanical/acoustic energy to
            // unmodeled receivers. Subtract it from work entering the coupled
            // magnetic-plus-coil-thermal store, never deposit it into copper.
            supplied_j: self.step.native.supplied - self.step.mechanical_acoustic_loss_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.driver, self.step)
    }
}

/// Complete jointly admitted speaker and thermal endpoints. Neither input moved.
pub struct MagneticSpeakerEnthalpyCandidate<'t, 'n, 'a> {
    thermal: DepositionSimulation<'t>,
    driver: ThermalMagneticTsDriver<'n, 'a>,
    report: CircuitEnthalpyReport<MaterialNetworkSample, MagneticSpeakerThermalStep>,
}
impl<'t, 'n, 'a> MagneticSpeakerEnthalpyCandidate<'t, 'n, 'a> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn driver(&self) -> &ThermalMagneticTsDriver<'n, 'a> {
        &self.driver
    }
    pub const fn report(
        &self,
    ) -> &CircuitEnthalpyReport<MaterialNetworkSample, MagneticSpeakerThermalStep> {
        &self.report
    }
    /// Publish BOTH histories; the coil heat has already been deposited once.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalMagneticTsDriver<'n, 'a>,
        CircuitEnthalpyReport<MaterialNetworkSample, MagneticSpeakerThermalStep>,
    ) {
        (self.thermal, self.driver, self.report)
    }
}

#[cfg(test)]
mod tests;
