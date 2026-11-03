//! Distributed phase transfer on the retained finite moving gas.
//!
//! All original donors are sampled before the single native transaction.
//! Condensation withdraws the old gas drift; evaporation leaves at the old
//! local liquid drift. Gross withdrawals use the old gas inventory, then all
//! retained and incoming gas mixes once. This is the existing frozen-donor,
//! first-order mean-drift closure, not a resolved vapor jet or angular gas model.

use std::collections::BTreeMap;

use fs_material::evaporation::gas::chamber::heat::GasLiquidChamberHeatSolve;
use fs_math::dd::Dd;

use super::super::multiple::{
    DistributedLiquidGasReport, DistributedLiquidGasSolve, LiquidGasHeatSurface,
    advance_native_with_momentum_and_flow,
};
use super::*;

/// Closed mass-transfer momentum and kinetic/receiver-heat account BEFORE flow.
/// Pressure traction, wall reactions and subsequent flow dissipation retain
/// their original owners and are not included in this source-only balance.
#[derive(Clone, Debug)]
pub struct DistributedPhaseMomentumReport {
    pub duration_s: f64,
    pub evaporated_mass_kg: f64,
    pub condensed_mass_kg: f64,
    pub gas_momentum_before_kg_m_s: [f64; 3],
    pub gas_momentum_after_kg_m_s: [f64; 3],
    pub gas_velocity_before_m_s: [f64; 3],
    pub gas_velocity_after_m_s: [f64; 3],
    pub liquid_kinetic_energy_change_j: f64,
    pub gas_kinetic_energy_change_j: f64,
    pub liquid_mixing_heat_j: f64,
    pub gas_mixing_heat_j: f64,
    pub gas_heat: Option<GasHeatStep>,
    pub momentum_residual_kg_m_s: [f64; 3],
    pub absolute_momentum_residual_kg_m_s: f64,
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub struct MovingGasDistributedEvaporationReport {
    /// The original combined source/native-flow account includes the mechanical
    /// source defect below, and subtracts each native mixing heat exactly once.
    pub native: DistributedLiquidGasReport,
    pub momentum: DistributedPhaseMomentumReport,
}

pub struct MovingGasDistributedEvaporationCandidate<'s, 'm, 'catalog> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: MovingGasDistributedEvaporationReport,
}

impl<'s, 'm, 'catalog> MovingGasDistributedEvaporationCandidate<'s, 'm, 'catalog> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &MovingGasDistributedEvaporationReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        MovingGasParcel,
        &'catalog ThermochemicalGasModel,
        MovingGasDistributedEvaporationReport,
    ) {
        (self.liquid, self.gas, self.catalog, self.report)
    }
}

/// Joint phase exchange from distinct native cells sharing ONE moving gas.
/// The original chamber selects every represented mass before its caloric
/// solve. One native material refresh includes inflows, comoving evaporation
/// and equilibrium zero-mass patches. Optional pressure uses the final mixed
/// gas EOS in the ONE subsequent original flow step. No accepted input changes
/// on any late source refusal, balance failure or cancellation.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_moving_gas<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasEvaporationSurface<'_>],
    controls: DistributedLiquidGasSolve,
    maximum_density_excursion: Option<f64>,
    gate: &CancelGate,
) -> Result<MovingGasDistributedEvaporationCandidate<'s, 'm, 'catalog>, EvaporatingLiquidGasError> {
    advance(
        liquid,
        gas,
        catalog,
        surfaces.iter().copied().map(|s| (s, None)),
        controls,
        None,
        maximum_density_excursion,
        gate,
    )
}

/// The same joint donor/mixing transaction after the original reciprocal
/// sensible-heat solve. Heat, phase and flow share one retained physical
/// interval; receiver mixing is a separate kinetic-to-internal conversion.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_moving_gas_and_heat<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasHeatSurface<'_>],
    controls: DistributedLiquidGasSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    maximum_density_excursion: Option<f64>,
    gate: &CancelGate,
) -> Result<MovingGasDistributedEvaporationCandidate<'s, 'm, 'catalog>, EvaporatingLiquidGasError> {
    advance(
        liquid,
        gas,
        catalog,
        surfaces.iter().map(|s| {
            (
                s.surface,
                Some((s.heat_transfer_coefficient_w_m2_k, s.heat_controls)),
            )
        }),
        controls,
        Some(heat_controls),
        maximum_density_excursion,
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance<'s, 'm, 'catalog, 'surface>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: impl ExactSizeIterator<
        Item = (
            LiquidGasEvaporationSurface<'surface>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    controls: DistributedLiquidGasSolve,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    maximum_density_excursion: Option<f64>,
    gate: &CancelGate,
) -> Result<MovingGasDistributedEvaporationCandidate<'s, 'm, 'catalog>, EvaporatingLiquidGasError> {
    let (candidate, motion) = advance_native_with_momentum_and_flow(
        liquid,
        gas.parcel(),
        catalog,
        surfaces,
        controls,
        heat_controls,
        maximum_density_excursion,
        Some(gas.momentum_kg_m_s()),
        gate,
        |liquid, next_gas| {
            super::super::advance_flow(
                liquid,
                gas.parcel(),
                next_gas,
                maximum_density_excursion,
                gate,
            )
        },
    )?;
    let momentum = motion.ok_or(EvaporatingLiquidGasError::Invalid(
        "distributed moving gas omitted its momentum account",
    ))?;
    let (liquid, parcel, catalog, native) = candidate.into_parts();
    let gas = MovingGasParcel::from_momentum(parcel, momentum.gas_momentum_after_kg_m_s)?;
    poll(gate)?;
    Ok(MovingGasDistributedEvaporationCandidate {
        liquid,
        gas,
        catalog,
        report: MovingGasDistributedEvaporationReport { native, momentum },
    })
}

pub(in crate::liquid_gas::evaporation) fn exchange(
    liquid: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
    original_gas: &ClosedGasParcel,
    gas: &mut ClosedGasParcel,
    original_momentum: [f64; 3],
    transfers: &[FreeSurfaceMassExchange],
    controls: DistributedLiquidGasSolve,
    gate: &CancelGate,
) -> Result<
    (
        FreeSurfaceMassExchangeReport,
        DistributedPhaseMomentumReport,
    ),
    EvaporatingLiquidGasError,
> {
    type Error = EvaporatingLiquidGasError;
    poll(gate)?;
    let gas_velocity = velocity(original_gas.mass_kg(), original_momentum)?;
    let mut incoming = Vec::new();
    let mut comoving = Vec::new();
    let mut donors = BTreeMap::new();
    for transfer in transfers {
        poll(gate)?;
        let donor = if transfer.mass_kg > 0.0 {
            incoming.push(FreeSurfaceMassInflow {
                exchange: *transfer,
                velocity_m_s: gas_velocity,
            });
            gas_velocity
        } else {
            comoving.push(*transfer);
            let [x, y, z] = transfer.coordinates;
            liquid
                .flow()
                .sample_at_si(x, y, z)
                .map_err(FreeSurfaceEnthalpyError::from)?
                .ok_or(Error::Invalid("phase donor has no native liquid velocity"))?
                .velocity_m_s
        };
        donors.insert(transfer.coordinates, donor);
    }
    let native =
        liquid.mix_interface_inflow_with_exchange(&incoming, &comoving, controls.native, gate)?;
    let samples: BTreeMap<_, _> = native
        .exchange
        .samples
        .iter()
        .map(|sample| (sample.coordinates, sample))
        .collect();
    let mut streams = Vec::with_capacity(donors.len());
    let mut evaporated = Dd::ZERO;
    let mut condensed = Dd::ZERO;
    // Native reports have tile order. Reduce gas sources in physical-coordinate
    // order, just like the original chamber, independently of caller ordering.
    for (coordinates, donor) in donors {
        let sample = samples
            .get(&coordinates)
            .ok_or(Error::Invalid("native phase batch omitted a donor"))?;
        let mass = -sample.mass_change_kg;
        if mass > 0.0 {
            evaporated = evaporated + Dd::from_f64(mass);
        }
        if mass < 0.0 {
            condensed = condensed - Dd::from_f64(mass);
        }
        streams.push(kinematics::GasStream {
            mass_kg: mass,
            velocity_m_s: donor,
        });
    }
    let motion = kinematics::advance_many(
        original_gas.mass_kg(),
        gas.mass_kg(),
        original_momentum,
        &streams,
    )
    .ok_or(Error::Invalid(
        "distributed gas recoil or mixing is unrepresentable",
    ))?;
    let gas_heat = if motion.mixing_heat_j != 0.0 {
        Some(
            gas.add_heat(
                motion.mixing_heat_j,
                controls.chamber.surface.gas.gas,
                || gate.is_requested(),
            )
            .map_err(gas_error)?,
        )
    } else {
        None
    };
    let residual: [f64; 3] = core::array::from_fn(|d| {
        (Dd::from_f64(native.exchange.momentum_change_kg_m_s[d]) + Dd::from_f64(motion.momentum[d])
            - Dd::from_f64(original_momentum[d]))
        .to_f64()
    });
    let absolute_momentum = residual.iter().map(|p| p.abs()).sum::<f64>();
    let energy = (Dd::from_f64(native.kinetic_energy_change_j)
        + Dd::from_f64(motion.kinetic_change_j)
        + Dd::from_f64(native.mixing_heat_j)
        + Dd::from_f64(gas_heat.as_ref().map_or(0.0, |h| h.stored_energy_change_j)))
    .to_f64();
    if !absolute_momentum.is_finite()
        || absolute_momentum > controls.native.momentum_tolerance_kg_m_s
    {
        return Err(Error::Invalid(
            "distributed phase momentum exceeded its native budget",
        ));
    }
    if !energy.is_finite() || energy.abs() > controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: energy,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok((
        native.exchange,
        DistributedPhaseMomentumReport {
            duration_s: liquid.flow().scaling().timestep_s(),
            evaporated_mass_kg: evaporated.to_f64(),
            condensed_mass_kg: condensed.to_f64(),
            gas_momentum_before_kg_m_s: original_momentum,
            gas_momentum_after_kg_m_s: motion.momentum,
            gas_velocity_before_m_s: gas_velocity,
            gas_velocity_after_m_s: motion.velocity,
            liquid_kinetic_energy_change_j: native.kinetic_energy_change_j,
            gas_kinetic_energy_change_j: motion.kinetic_change_j,
            liquid_mixing_heat_j: native.mixing_heat_j,
            gas_mixing_heat_j: motion.mixing_heat_j,
            gas_heat,
            momentum_residual_kg_m_s: residual,
            absolute_momentum_residual_kg_m_s: absolute_momentum,
            energy_residual_j: energy,
        },
    ))
}
