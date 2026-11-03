//! Phase transfer and the original finite-atmosphere thermal/flow timeline.
//!
//! The phase owner selects represented native mass and updates the same liquid
//! caloric field. Its post-phase gas then enters the existing finite-gas,
//! finite-wall conduction owner without reconstruction or a second flow step.

use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat,
    conduction::{
        FreeSurfaceBodyContact, FreeSurfaceBodyInterface, FreeSurfaceConductionControls,
        FreeSurfaceConductionFace, FreeSurfaceGasFace, FreeSurfaceRadiationFace,
        FreeSurfaceThermalBody,
    },
    viscous::FreeSurfaceViscousConductionReport,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

use super::{
    CancelGate, ClosedGasParcel, EvaporatingLiquidGasError, LiquidGasEvaporationSurface,
    MaterialFreeSurfaceEnthalpy3, ThermochemicalGasModel,
    multiple::{DistributedLiquidGasReport, DistributedLiquidGasSolve, advance_native_with_flow},
    poll,
};

type Error = EvaporatingLiquidGasError;

/// Original phase and thermal accounts over the SAME physical interval.
#[derive(Clone, Debug)]
pub struct ConductingLiquidGasReport {
    /// Source/native phase account followed by the one full thermal/flow step.
    /// Its exchange gas endpoint precedes conduction, and its flow defect
    /// already counts the native thermal/flow ledger once.
    pub phase: DistributedLiquidGasReport,
    /// Original liquid, finite solid, finite gas and optional electrical or
    /// molecular heat receipts. Its flow is exactly the one in `phase`.
    pub thermal: FreeSurfaceViscousConductionReport,
    /// [Post-phase, post-thermal] temperatures of the same original gas [K].
    pub thermal_gas_temperature_k: [f64; 2],
    /// Phase/native defect plus ONLY the original finite-gas caloric defect.
    /// Adding the full combined thermal defect would count native flow twice.
    pub energy_residual_j: f64,
    /// Absolute phase/native defects plus the absolute finite-gas defect.
    pub absolute_energy_defect_j: f64,
}

/// Liquid, original finite gas and finite solids publish as one endpoint.
pub struct ConductingLiquidGasCandidate<'s, 'm, 'catalog, 'b> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    bodies: Vec<FreeSurfaceThermalBody<'b>>,
    report: ConductingLiquidGasReport,
}

impl<'s, 'm, 'catalog, 'b> ConductingLiquidGasCandidate<'s, 'm, 'catalog, 'b> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub fn bodies(&self) -> &[FreeSurfaceThermalBody<'b>] {
        &self.bodies
    }
    pub const fn report(&self) -> &ConductingLiquidGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        Vec<FreeSurfaceThermalBody<'b>>,
        ConductingLiquidGasReport,
    ) {
        (
            self.liquid,
            self.gas,
            self.catalog,
            self.bodies,
            self.report,
        )
    }
}

/// Exchange every phase patch with one finite gas, then conduct heat and
/// advance native flow ONCE over the same retained physical timestep.
///
/// A single surface is a one-element slice. The original distributed phase
/// owner admits native mass, one volatile/datum binding and the common gas
/// donor before the original thermal owner sees those accepted phase states.
/// There is no separate authored sensible-film coefficient: `exposed` and
/// solid gas faces use the original source-backed thermal conductances.
///
/// The full thermal owner includes liquid conduction, finite solid storage and
/// contacts, solid gas films/radiation, attached electrical current, and optional
/// molecular heating. Every original source/geometry/work limit still applies.
/// Its gas conversion is exact: no new atmosphere or composition is inferred.
/// During thermal exchange mass, composition and physical gas volume are fixed.
///
/// Optional pressure uses that SAME post-thermal gas EOS for the one native
/// momentum/VOF step. Its pressure receipt begins AFTER phase transfer, unlike
/// the standalone phase/pressure route. An explicit fixed pressure reference
/// and stationary solid geometry remain required by the original pressure owner.
///
/// This is a first-order phase/thermal/flow split, not a simultaneous implicit
/// phase/heat solve; refine the original timestep to assess its accuracy. The
/// original well-mixed fixed-volume atmosphere and phase geometry assumptions
/// remain explicit. Solved surface-stress/contact histories, dry-out, interface
/// overflow and unsupported source changes retain their original refusals.
/// The thermal balance is not a total mechanical-energy or compression-work
/// claim. Any refusal or cancellation discards liquid, gas and solids together.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_gas_and_conduction<'s, 'm, 'catalog, 'b>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasEvaporationSurface<'_>],
    bodies: &[FreeSurfaceThermalBody<'b>],
    faces: &[FreeSurfaceConductionFace],
    heat: &[FreeSurfaceCellHeat],
    contacts: &[FreeSurfaceBodyContact],
    interfaces: &[FreeSurfaceBodyInterface<'_>],
    exposed: &[FreeSurfaceGasFace],
    phase_controls: DistributedLiquidGasSolve,
    thermal_controls: FreeSurfaceConductionControls,
    gas_controls: GasEnergySolve,
    maximum_density_excursion: Option<f64>,
    with_viscous_heating: bool,
    gate: &CancelGate,
) -> Result<ConductingLiquidGasCandidate<'s, 'm, 'catalog, 'b>, Error> {
    advance_liquid_surfaces_with_gas_and_radiative_conduction(
        liquid,
        gas,
        catalog,
        surfaces,
        bodies,
        faces,
        heat,
        contacts,
        interfaces,
        exposed,
        &[],
        phase_controls,
        thermal_controls,
        gas_controls,
        maximum_density_excursion,
        with_viscous_heating,
        gate,
    )
}

/// Add original resolved-cell liquid radiation to the same phase/thermal step.
///
/// Phase and thermal radiation use the native liquid caloric temperature. The
/// original gray-radiation law requires a resolved liquid thermal boundary layer
/// or an explicitly justified lumped surface value; this entry does not invent
/// a distinct interface temperature or a latent sink in the solved-surface law.
/// Its radiation is external supply and its gas-film transfer is internal to the
/// retained gas/liquid inventory. Both enter the original native thermal owner
/// exactly once, including combined thermal stability and source admission.
///
/// All phase/conduction source, fixed-volume, optional pressure, atomic state
/// and first-order splitting limits of
/// [`advance_liquid_surfaces_with_gas_and_conduction`] remain in force. Bound
/// solved-surface stress and wetting histories retain their original refusal.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_gas_and_radiative_conduction<'s, 'm, 'catalog, 'b>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasEvaporationSurface<'_>],
    bodies: &[FreeSurfaceThermalBody<'b>],
    faces: &[FreeSurfaceConductionFace],
    heat: &[FreeSurfaceCellHeat],
    contacts: &[FreeSurfaceBodyContact],
    interfaces: &[FreeSurfaceBodyInterface<'_>],
    exposed: &[FreeSurfaceGasFace],
    radiation: &[FreeSurfaceRadiationFace],
    phase_controls: DistributedLiquidGasSolve,
    thermal_controls: FreeSurfaceConductionControls,
    gas_controls: GasEnergySolve,
    maximum_density_excursion: Option<f64>,
    with_viscous_heating: bool,
    gate: &CancelGate,
) -> Result<ConductingLiquidGasCandidate<'s, 'm, 'catalog, 'b>, Error> {
    poll(gate)?;
    let mut next_bodies = bodies.to_vec();
    let mut thermal_result = None;
    let candidate = advance_native_with_flow(
        liquid,
        gas,
        catalog,
        surfaces.iter().copied().map(|surface| (surface, None)),
        phase_controls,
        None,
        maximum_density_excursion,
        gate,
        |next_liquid, phase_gas| {
            if phase_gas.volume_m3() != gas.volume_m3() {
                return Err(Error::Invalid(
                    "phase exchange changed the original fixed gas volume",
                ));
            }
            // ClosedGasParcel wraps this exact isochoric owner. Keep the
            // post-phase parcel until the entire thermal candidate succeeds.
            let mut thermal_gas = phase_gas.clone().into_isochoric();
            let report = next_liquid.step_with_radiative_finite_gas_environment(
                faces,
                heat,
                thermal_controls,
                &mut next_bodies,
                contacts,
                interfaces,
                &mut thermal_gas,
                exposed,
                radiation,
                gas_controls,
                maximum_density_excursion,
                with_viscous_heating,
                gate,
            )?;
            if thermal_gas.mass_kg() != phase_gas.mass_kg()
                || thermal_gas.volume_m3() != phase_gas.volume_m3()
                || thermal_gas.model() != phase_gas.model()
            {
                return Err(Error::Invalid(
                    "thermal exchange changed original gas mass, volume or composition sources",
                ));
            }
            let finite = report.conduction.finite_gas.as_ref().ok_or(Error::Invalid(
                "native conduction omitted its original finite-gas account",
            ))?;
            let flow = report.conduction.flow;
            let pressure = finite.pressure;
            let temperatures = [
                phase_gas.state().state().temperature,
                thermal_gas.state().state().temperature,
            ];
            *phase_gas = ClosedGasParcel::from_isochoric(thermal_gas);
            thermal_result = Some((report, temperatures));
            Ok((flow, pressure))
        },
    )?;
    let (liquid, gas, catalog, phase) = candidate.into_parts();
    let (thermal, thermal_gas_temperature_k) = thermal_result.ok_or(Error::Invalid(
        "phase transaction omitted its native thermal/flow step",
    ))?;
    let finite = thermal
        .conduction
        .finite_gas
        .as_ref()
        .ok_or(Error::Invalid(
            "native conduction omitted its original finite-gas account",
        ))?;
    // The phase account already includes conduction.flow.balance_defect_j.
    // Its gas endpoint ends before thermal exchange, so only the subsequent
    // original finite-gas storage defect belongs on this enclosing boundary.
    let energy_residual_j = phase.energy_residual_j + finite.gas_balance_defect_j;
    let absolute_energy_defect_j =
        phase.absolute_energy_defect_j + finite.gas_balance_defect_j.abs();
    if ![energy_residual_j, absolute_energy_defect_j]
        .into_iter()
        .all(f64::is_finite)
    {
        return Err(Error::Invalid(
            "combined phase/conduction energy account is unrepresentable",
        ));
    }
    let worst = energy_residual_j.abs().max(absolute_energy_defect_j);
    if worst > phase_controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: worst,
            tolerance_j: phase_controls.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok(ConductingLiquidGasCandidate {
        liquid,
        gas,
        catalog,
        bodies: next_bodies,
        report: ConductingLiquidGasReport {
            phase,
            thermal,
            thermal_gas_temperature_k,
            energy_residual_j,
            absolute_energy_defect_j,
        },
    })
}
