//! Finite wall heat replenishment and freely recoiling molecular escape.
//!
//! One first-order Lie-split interval applies the original implicit wall/gas
//! heat transaction, then the original effusion/recoil transaction. Both stages
//! represent the SAME physical duration, not two elapsed intervals. The heated
//! gas supplies the emitted species, momentum and energy; its subsequent cooling
//! feeds the next interval's wall exchange. Time refinement remains necessary.
//!
//! The wall is part of the vessel's existing dry mass. Its thermal storage is
//! added once, without adding its mass or translational kinetic energy again.
//! Constant-density, fixed-volume charts make the wall's enthalpy increment an
//! internal-energy increment; expanding walls require a separate work owner.
//! This does not infer a rarefied heat-transfer coefficient, wall accommodation,
//! pressure integrity, rotation or spatial thermal equilibration. Original
//! source, Biot, caloric, finite-inventory and mechanical limits remain active.

/// Finite-wall heat with original bidirectional and drifting molecular ports.
pub mod exchange;
/// Error-controlled physical windows with finite wall, gas and vehicle feedback.
pub mod trajectory;

use super::{
    EffusionOutlet, EffusionRocketEnvironment, EffusionRocketError, EffusionRocketSolve,
    EffusionRocketStep, GasRocketState, MAX_SNAPSHOT_TRANSFERS, Sum,
};
use crate::{
    cavity_phs::piston::thermal::HeatExchangeSolve,
    finite_gas_thermal::{
        ConvectionCoefficient, ConvectionError, FiniteGasThermalEnvironment, FiniteGasThermalReport,
    },
};
use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_exec::Cx;
use fs_material::{
    phase::{EquilibriumPhaseState, UniformEnthalpyStep, UniformEnthalpyStepInput},
    thermochemical_gas::{
        ResolvedThermochemicalGasState,
        energy::{
            GasEnergySolve,
            parcel::{ClosedGasParcel, transfer::species::effusion::EffusiveInterface},
        },
    },
};

/// Original stage controls plus one noncancelling budget for the whole proposal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffusionWallSolve {
    /// Used by both heat inversion and the subsequent finite-gas reconstruction.
    pub gas: GasEnergySolve,
    pub exchange: HeatExchangeSolve,
    pub maximum_gas_temperature_increment_k: f64,
    pub recoil: EffusionRocketSolve,
    /// Sum of absolute stage energy defects plus native rounding reserves [J].
    /// This is a balance limit, not a temporal-error or source-accuracy estimate.
    pub combined_energy_tolerance_j: f64,
}

/// Original owners' refusals remain distinct; neither endpoint is published.
#[derive(Debug)]
pub enum EffusionWallError<E = core::convert::Infallible> {
    Invalid(&'static str),
    Thermal(ConvectionError<E>),
    Effusion(EffusionRocketError),
    EnergyBalance {
        residual_j: f64,
        absolute_residual_sum_j: f64,
        uncertainty_bound_j: f64,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for EffusionWallError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-wall effusive recoil: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for EffusionWallError<E> {}
impl<E> From<EffusionRocketError> for EffusionWallError<E> {
    fn from(error: EffusionRocketError) -> Self {
        Self::Effusion(error)
    }
}

/// Both actual owner reports for a single physical interval.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallStep<S = ()> {
    pub duration_s: f64,
    /// Contains the accepted wall and the heated gas BEFORE molecular escape.
    pub thermal: UniformEnthalpyStep<FiniteGasThermalReport<S>>,
    pub recoil: EffusionRocketStep,
    /// Body internal source + gas source + enclosing radiation [J]. Reciprocal
    /// wall/gas heat cancels. Mechanical work remains in `recoil.exhaust`.
    pub external_thermal_heat_j: f64,
    /// Nominal ΔU_gas + ΔU_wall + ΔK_vehicle + E_escaped - Q_external - W_external.
    /// Composes original reference-free storage/transport accounts; native
    /// mass-chart/gravity uncertainty is reserved separately, never called heat.
    pub combined_energy_residual_j: f64,
    pub absolute_energy_residual_sum_j: f64,
    /// Initial native mass remapping and unresolved gravity mass-history work.
    pub energy_uncertainty_bound_j: f64,
}

/// Publish the finite wall and the gas/vehicle together after acceptance.
#[derive(Clone, Debug, PartialEq)]
pub struct EffusionWallCandidate<S = ()> {
    endpoint: GasRocketState,
    step: EffusionWallStep<S>,
}
impl<S> EffusionWallCandidate<S> {
    #[must_use]
    pub const fn endpoint(&self) -> &GasRocketState {
        &self.endpoint
    }
    #[must_use]
    pub const fn wall(&self) -> &EquilibriumPhaseState {
        &self.step.thermal.state
    }
    #[must_use]
    pub const fn step(&self) -> &EffusionWallStep<S> {
        &self.step
    }
    #[must_use]
    pub fn into_parts(self) -> (GasRocketState, EquilibriumPhaseState, EffusionWallStep<S>) {
        (self.endpoint, self.step.thermal.state, self.step)
    }
}

// The vacuum and bidirectional wrappers share these exact finite storage and
// caloric owners. Only their original particle-exchange transaction differs.
fn admit_wall<E>(
    original: &GasRocketState,
    cx: &Cx<'_>,
    environment: EffusionRocketEnvironment,
    wall: UniformEnthalpyStepInput<'_>,
    controls: EffusionWallSolve,
    port_count: usize,
) -> Result<(), EffusionWallError<E>> {
    cx.checkpoint().map_err(|_| EffusionWallError::Cancelled)?;
    let density = wall.initial.bulk_density_kg_m3();
    let material_volume = wall.mass_kg / density;
    if !wall.mass_kg.is_finite()
        || wall.mass_kg <= 0.0
        || wall.mass_kg > original.vehicle.dry_mass_kg
        || !material_volume.is_finite()
        || material_volume <= 0.0
        || !wall.volume_m3.is_finite()
        || wall.volume_m3 <= 0.0
        || (wall.volume_m3 - material_volume).abs()
            > 64.0 * f64::EPSILON * wall.volume_m3.max(material_volume)
        || !wall.duration_s.is_finite()
        || wall.duration_s <= 0.0
        || !controls.combined_energy_tolerance_j.is_finite()
        || controls.combined_energy_tolerance_j < 0.0
        || environment.inertial_frame != original.inertial_frame
        || port_count > MAX_SNAPSHOT_TRANSFERS
    {
        return Err(EffusionWallError::Invalid(
            "finite wall must fit dry mass and fixed material volume; require positive duration, matching frame and finite energy budget",
        ));
    }
    for knot in wall.curve.knots() {
        cx.checkpoint().map_err(|_| EffusionWallError::Cancelled)?;
        if knot.bulk_density_kg_m3 != density {
            return Err(EffusionWallError::Invalid(
                "fixed-volume wall requires constant-density calorics; expansion needs its work owner",
            ));
        }
    }
    Ok(())
}

fn admit_original_gas<'a, E>(
    original: &GasRocketState,
    cx: &Cx<'_>,
    apertures: impl IntoIterator<Item = (u64, f64, &'a EffusiveInterface)>,
) -> Result<(), EffusionWallError<E>> {
    // Heating must not hide an initially unqualified active aperture. Closed
    // paths retain the original no-source-query behavior.
    for (id, area_m2, interface) in apertures {
        cx.checkpoint().map_err(|_| EffusionWallError::Cancelled)?;
        if !area_m2.is_finite() || area_m2 < 0.0 {
            return Err(EffusionWallError::Invalid(
                "finite nonnegative outlet area required",
            ));
        }
        if area_m2 > 0.0 {
            interface
                .admit_parcel(&original.reservoir, || cx.checkpoint().is_err())
                .map_err(|error| EffusionRocketError::Source { outlet: id, error })?;
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn exchange_heat<S, E>(
    original: &GasRocketState,
    cx: &Cx<'_>,
    wall: UniformEnthalpyStepInput<'_>,
    body: &LumpedThermalEnvironment,
    gas_heat_j: f64,
    controls: EffusionWallSolve,
    resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConvectionCoefficient<S>, E>,
) -> Result<
    (
        GasRocketState,
        UniformEnthalpyStep<FiniteGasThermalReport<S>>,
    ),
    EffusionWallError<E>,
> {
    let gas = original.reservoir.clone().into_isochoric();
    let thermal = FiniteGasThermalEnvironment {
        gas: &gas,
        body,
        gas_heat_j,
        gas_solve: controls.gas,
        exchange: controls.exchange,
        maximum_gas_temperature_increment_k: controls.maximum_gas_temperature_increment_k,
    }
    .advance_with_convection(cx, wall, resolve)
    .map_err(EffusionWallError::Thermal)?;
    let mut heated = original.clone();
    heated.reservoir = ClosedGasParcel::from_isochoric(thermal.report.next_gas.clone());
    if heated.reservoir.mass_kg().to_bits() != original.reservoir.mass_kg().to_bits()
        || heated.reservoir.volume_m3().to_bits() != original.reservoir.volume_m3().to_bits()
    {
        return Err(EffusionWallError::Invalid(
            "thermal stage changed fixed gas inventory or volume",
        ));
    }
    Ok((heated, thermal))
}

struct WallEnergyAccount {
    external_thermal_heat_j: f64,
    combined_energy_residual_j: f64,
    absolute_energy_residual_sum_j: f64,
    energy_uncertainty_bound_j: f64,
}

fn energy_account<S, E>(
    thermal: &UniformEnthalpyStep<FiniteGasThermalReport<S>>,
    wall: UniformEnthalpyStepInput<'_>,
    gas_heat_j: f64,
    native_energy_residual_j: f64,
    energy_uncertainty_bound_j: f64,
    tolerance_j: f64,
) -> Result<WallEnergyAccount, EffusionWallError<E>> {
    let mut energy = Sum::default();
    energy.add(thermal.report.combined_energy_residual_j)?;
    energy.add(native_energy_residual_j)?;
    let combined_energy_residual_j = energy.total();
    let absolute_energy_residual_sum_j =
        thermal.report.combined_energy_residual_j.abs() + native_energy_residual_j.abs();
    if !absolute_energy_residual_sum_j.is_finite()
        || !energy_uncertainty_bound_j.is_finite()
        || absolute_energy_residual_sum_j > tolerance_j
        || energy_uncertainty_bound_j > tolerance_j - absolute_energy_residual_sum_j
    {
        return Err(EffusionWallError::EnergyBalance {
            residual_j: combined_energy_residual_j,
            absolute_residual_sum_j: absolute_energy_residual_sum_j,
            uncertainty_bound_j: energy_uncertainty_bound_j,
        });
    }
    let mut external_heat = Sum::default();
    for heat in [
        wall.internal_heat_j,
        gas_heat_j,
        thermal.report.radiation_heat_j,
    ] {
        external_heat.add(heat)?;
    }
    Ok(WallEnergyAccount {
        external_thermal_heat_j: external_heat.total(),
        combined_energy_residual_j,
        absolute_energy_residual_sum_j,
        energy_uncertainty_bound_j,
    })
}

impl GasRocketState {
    /// Exchange heat with a finite material wall and emit from the heated gas.
    ///
    /// `wall.mass_kg` is already included in this vessel's dry mass. The body
    /// geometry and all wall-chart densities must agree with that invariant
    /// mass. The original body owner supplies implicit heat storage, radiation
    /// and Biot admission; `body.temperature_k` is replaced by the trial gas
    /// temperature. The caller supplies a physically qualified whole-surface h.
    /// Prescribed `gas_heat_j` excludes the reciprocal wall transfer.
    ///
    /// Heating then escape is first-order operator splitting. Both stages use
    /// `wall.duration_s`; this does not assert that heat precedes every physical
    /// molecular event. Refine the interval for coupled cooling/heat feedback.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_with_wall(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        outlets: &[EffusionOutlet<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        controls: EffusionWallSolve,
    ) -> Result<EffusionWallCandidate, EffusionWallError> {
        self.advance_effusion_with_wall_convection(
            cx,
            environment,
            outlets,
            wall,
            body,
            gas_heat_j,
            controls,
            |_| {
                Ok(ConvectionCoefficient {
                    coefficient_w_per_m2_k: body.convection_w_per_m2_k,
                    source: (),
                })
            },
        )
    }

    /// Resolve a repeatable convection law against each original thermal trial.
    ///
    /// The callback sees the actual heated gas's caloric/transport/composition
    /// state and retains its concrete source evaluation. Its validity, including
    /// rarefaction, remains the caller's responsibility. It cannot depend on an
    /// unresolved wall endpoint or assume that no rejected trials occur.
    /// Callback side effects are not rolled back. Physical inputs are immutable.
    ///
    /// Original and heated gas states must satisfy every open aperture source.
    /// Original owners admit their final states too. The combined balance adds
    /// stage defects without subtracting absolute formation energies; opposite
    /// errors cannot cancel their way through the global budget. External heat
    /// carries no separately modeled photon or other mechanical impulse here.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_effusion_with_wall_convection<S, E>(
        &self,
        cx: &Cx<'_>,
        environment: EffusionRocketEnvironment,
        outlets: &[EffusionOutlet<'_>],
        wall: UniformEnthalpyStepInput<'_>,
        body: &LumpedThermalEnvironment,
        gas_heat_j: f64,
        controls: EffusionWallSolve,
        resolve: impl FnMut(&ResolvedThermochemicalGasState) -> Result<ConvectionCoefficient<S>, E>,
    ) -> Result<EffusionWallCandidate<S>, EffusionWallError<E>> {
        admit_wall(self, cx, environment, wall, controls, outlets.len())?;
        admit_original_gas(
            self,
            cx,
            outlets
                .iter()
                .map(|outlet| (outlet.id, outlet.area_m2, outlet.interface.interface())),
        )?;
        let (heated, thermal) = exchange_heat(self, cx, wall, body, gas_heat_j, controls, resolve)?;
        let (endpoint, recoil) = heated
            .advance_effusion_with_environment(
                environment,
                outlets,
                wall.duration_s,
                controls.gas,
                controls.recoil,
                || cx.checkpoint().is_err(),
            )?
            .into_parts();
        let account = energy_account(
            &thermal,
            wall,
            gas_heat_j,
            recoil.energy_residual_j,
            recoil.initial_inventory_kinetic_correction_j.abs() + recoil.gravity_work_uncertainty_j,
            controls.combined_energy_tolerance_j,
        )?;
        cx.checkpoint().map_err(|_| EffusionWallError::Cancelled)?;
        Ok(EffusionWallCandidate {
            endpoint,
            step: EffusionWallStep {
                duration_s: wall.duration_s,
                thermal,
                recoil,
                external_thermal_heat_j: account.external_thermal_heat_j,
                combined_energy_residual_j: account.combined_energy_residual_j,
                absolute_energy_residual_sum_j: account.absolute_energy_residual_sum_j,
                energy_uncertainty_bound_j: account.energy_uncertainty_bound_j,
            },
        })
    }
}
