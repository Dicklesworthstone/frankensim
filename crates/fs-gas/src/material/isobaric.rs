//! Finite-rate chemistry at prescribed pressure on the shared material parcel.
//!
//! The original reactor advances species and reservoir heat together with H as
//! its energy constraint. The material owner resolves that composition and heat
//! at the original pressure, including the actual expansion work. Formation
//! enthalpy is already in those two source models; it is never added as heat.

use super::{
    ClosedGasParcel, MaterialKineticsError, MaterialKineticsSolve, ThermochemicalGasModel,
    admit_bank, admit_sources, compensated_sum, poll, positive,
};
use crate::chemistry::heat::{HeatBath, material::FiniteHeatBody};
use crate::chemistry::isobaric::{IsobaricKineticsView, IsobaricReactor};
use fs_material::phase::isobaric::IsobaricHeatStep;
use fs_material::thermochemical_gas::energy::parcel::isobaric::composition::GasIsobaricCompositionStep;
use fs_qty::{Amount, Pressure, Temperature, Time};
use fs_thermochem::kinetics::IdealGasKineticsV1;

type Error = MaterialKineticsError;

/// Cross-owner balances at the actual accepted material endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialIsobaricBalance {
    /// Reconstructed feed mass minus the exactly preserved parcel mass [kg].
    pub feed_mass_residual_kg: f64,
    /// Integrated chemical mass minus the preserved parcel mass [kg].
    pub chemical_mass_residual_kg: f64,
    /// Actual final material density times volume minus parcel mass [kg].
    pub density_mass_residual_kg: f64,
    /// Original chemical feed H minus initial material H [J].
    pub initial_enthalpy_residual_j: f64,
    /// Actual final material H minus the chemical endpoint H [J].
    pub freeze_enthalpy_residual_j: f64,
    /// Actual final material U minus the chemical endpoint U [J].
    pub freeze_internal_energy_residual_j: f64,
    /// Material p Delta V minus the chemical trajectory's p Delta V [J].
    pub freeze_boundary_work_residual_j: f64,
    /// Chemical initial volume relative to the original material volume.
    pub relative_initial_volume_correction: f64,
    /// Actual material volume relative to the chemical endpoint volume.
    pub relative_volume_correction: f64,
    /// Actual material temperature relative to the chemical endpoint temperature.
    pub relative_temperature_correction: f64,
    /// Largest absolute fraction change when freezing the exact simplex.
    pub maximum_mole_fraction_correction: f64,
}

/// Accepted common material state and the original prepared chemistry report.
///
/// The chemical view borrows the reactor until the next advance. Its existing
/// `to_owned` method can retain an independent mechanism/receipt snapshot.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialIsobaricStep<'a> {
    pub chemistry: IsobaricKineticsView<'a>,
    /// Actual formation-inclusive Delta H, Delta U and p Delta V, including
    /// both independent first-law residuals from the original material owner.
    pub material: GasIsobaricCompositionStep,
    pub balance: MaterialIsobaricBalance,
}

/// Simultaneously accepted gas and finite-body endpoints.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialIsobaricBodyStep<'a, 'm> {
    pub gas: MaterialIsobaricStep<'a>,
    /// Continue this same original caloric source at the next interval.
    pub body: FiniteHeatBody<'m>,
    /// Original sensible/latent heat, density change and body expansion work.
    pub body_report: IsobaricHeatStep,
    pub body_entropy_change_j_k: f64,
    /// Actual Delta Hgas + Delta Hbody [J].
    pub combined_enthalpy_residual_j: f64,
    /// Actual Delta Ugas + Delta Ubody + Wgas + Wbody [J].
    pub combined_energy_residual_j: f64,
}

/// Reusable constant-pressure reactor for source-bearing material gas parcels.
///
/// The original kinetic preparation, Newton storage and dense feed buffer are
/// reused across intervals and independent parcels. Full transport/acoustic
/// material resolution happens at the endpoint; it is not an allocation-free
/// path. The bank must cover every possible product with the SAME NASA-9 laws.
///
/// A quasistatic external load holds the parcel's current pressure. This owner
/// does not resolve that load's inertia, spatial gradients, streams or a valve.
/// The original kinetics still determines composition, heat and entropy. No
/// equilibrium projection or constant-Cp temperature update is introduced.
pub struct MaterialIsobaricReactor {
    bank: ThermochemicalGasModel,
    reactor: IsobaricReactor,
    controls: MaterialKineticsSolve,
    amounts: Vec<Amount>,
}

impl MaterialIsobaricReactor {
    pub fn new(
        bank: ThermochemicalGasModel,
        mechanism: IdealGasKineticsV1,
        controls: MaterialKineticsSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        controls.material_energy.validate()?;
        if !positive(controls.mass_tolerance_kg)
            || !positive(controls.relative_state_tolerance)
            || controls.relative_state_tolerance >= 1.0
        {
            return Err(Error::Invalid(
                "positive finite mass and relative state tolerances below one required",
            ));
        }
        admit_bank(&bank, &mechanism)?;
        let amounts = vec![Amount::new(0.0); mechanism.thermodynamics().models().len()];
        let reactor = IsobaricReactor::new(mechanism, controls.chemistry, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            bank,
            reactor,
            controls,
            amounts,
        })
    }

    pub const fn bank(&self) -> &ThermochemicalGasModel {
        &self.bank
    }
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }
    pub const fn controls(&self) -> MaterialKineticsSolve {
        self.controls
    }
    /// Reusable numerical storage; excludes immutable source banks and reports.
    pub fn workspace_bytes(&self) -> usize {
        self.reactor.workspace_bytes() + self.amounts.capacity() * core::mem::size_of::<Amount>()
    }

    /// Advance with H conserved and outward p Delta V taken from actual volume.
    /// Mass and pressure are preserved exactly. Zero duration preserves the
    /// complete parcel after source, control and cancellation admission.
    pub fn advance_adiabatic<'a>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialIsobaricStep<'a>, Error> {
        let feed_mass = self.prepare(parcel, &mut cancelled)?;
        let state = parcel.state().state();
        let chemistry = self.reactor.advance_adiabatic(
            &self.amounts,
            Pressure::new(state.pressure),
            Temperature::new(state.temperature),
            duration,
            &mut cancelled,
        )?;
        let (candidate, step) = reconcile(
            parcel,
            &self.bank,
            self.controls,
            feed_mass,
            chemistry,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(step)
    }

    /// Integrate chemistry and G(Tbath - Tgas) in the same implicit stages.
    /// The reservoir heat enters H once; its counterpart entropy stays in the
    /// original chemistry report. This does not split chemistry from heating.
    pub fn advance_with_heat_bath<'a>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        bath: HeatBath,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialIsobaricStep<'a>, Error> {
        let feed_mass = self.prepare(parcel, &mut cancelled)?;
        let state = parcel.state().state();
        let chemistry = self.reactor.advance_with_heat_bath(
            &self.amounts,
            Pressure::new(state.pressure),
            Temperature::new(state.temperature),
            duration,
            bath,
            &mut cancelled,
        )?;
        let (candidate, step) = reconcile(
            parcel,
            &self.bank,
            self.controls,
            feed_mass,
            chemistry,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(step)
    }

    /// Integrate chemistry and a finite original caloric body simultaneously.
    /// Its latent plateaus and density-dependent expansion stay in the existing
    /// body owner. Both endpoints remain unpublished if ANY final check fails.
    pub fn advance_with_finite_body<'a, 'm>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        body: FiniteHeatBody<'m>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialIsobaricBodyStep<'a, 'm>, Error> {
        let feed_mass = self.prepare(parcel, &mut cancelled)?;
        let state = parcel.state().state();
        let coupled = self.reactor.advance_with_finite_body(
            &self.amounts,
            Pressure::new(state.pressure),
            Temperature::new(state.temperature),
            duration,
            body,
            &mut cancelled,
        )?;
        let body_entropy_change_j_k = coupled.body_entropy_change_j_k();
        let (chemistry, body, body_report) = coupled.into_parts();
        let (candidate, gas) = reconcile(
            parcel,
            &self.bank,
            self.controls,
            feed_mass,
            chemistry,
            &mut cancelled,
        )?;
        let combined_enthalpy_residual_j =
            gas.material.enthalpy_change_j + body_report.enthalpy_change_j;
        let combined_energy_residual_j = compensated_sum(
            [
                gas.material.energy.stored_energy_change_j,
                body_report.internal_energy_change_j,
                gas.material.energy.work_by_gas_j,
                body_report.boundary_work_j,
            ]
            .into_iter(),
        );
        let tolerance = parcel.mass_kg() * self.controls.material_energy.energy_tolerance_j_kg
            + self.controls.chemistry.energy_tolerance_j;
        if !combined_enthalpy_residual_j.is_finite()
            || !combined_energy_residual_j.is_finite()
            || combined_enthalpy_residual_j.abs() > tolerance
            || combined_energy_residual_j.abs() > tolerance
        {
            return Err(Error::Invalid(
                "material gas and finite body do not close combined heat/work",
            ));
        }
        let step = MaterialIsobaricBodyStep {
            gas,
            body,
            body_report,
            body_entropy_change_j_k,
            combined_enthalpy_residual_j,
            combined_energy_residual_j,
        };
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(step)
    }

    fn prepare(
        &mut self,
        parcel: &ClosedGasParcel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        admit_sources(parcel, &self.bank, self.reactor.mechanism())?;
        let amount = parcel.mass_kg()
            / parcel
                .state()
                .thermodynamics()
                .properties()
                .molar_mass()
                .value();
        if !positive(amount) {
            return Err(Error::Invalid("unrepresentable material feed amount"));
        }
        self.amounts.fill(Amount::new(0.0));
        let original = parcel.model().thermodynamics();
        let models = self.bank.thermodynamics().components();
        for (source, &fraction) in original
            .components()
            .iter()
            .zip(original.mole_composition().fractions())
        {
            // Both original owners use canonical SpeciesId order. Reusing the
            // dense feed avoids constructing/cloning a sparse species map.
            let index = models
                .binary_search_by(|model| model.species().cmp(source.species()))
                .expect("admitted active source belongs to the complete bank");
            let value = amount * fraction;
            if !positive(value) {
                return Err(Error::Invalid(
                    "positive feed species amount is unrepresentable",
                ));
            }
            self.amounts[index] = Amount::new(value);
        }
        let feed_mass = compensated_sum(
            models
                .iter()
                .zip(&self.amounts)
                .map(|(model, amount)| model.molar_mass().value() * amount.value()),
        );
        if !feed_mass.is_finite()
            || (feed_mass - parcel.mass_kg()).abs() > self.controls.mass_tolerance_kg
        {
            return Err(Error::Invalid(
                "reconstructed feed does not preserve parcel mass",
            ));
        }
        poll(cancelled)?;
        Ok(feed_mass)
    }
}

fn reconcile<'a>(
    parcel: &ClosedGasParcel,
    bank: &ThermochemicalGasModel,
    controls: MaterialKineticsSolve,
    feed_mass: f64,
    chemistry: IsobaricKineticsView<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ClosedGasParcel, MaterialIsobaricStep<'a>), Error> {
    poll(cancelled)?;
    let mass = parcel.mass_kg();
    let tolerance = mass * controls.material_energy.energy_tolerance_j_kg
        + controls.chemistry.energy_tolerance_j;
    if !positive(tolerance) {
        return Err(Error::Invalid(
            "unrepresentable material/chemical energy tolerance",
        ));
    }
    let products = if chemistry.elapsed_time().value() == 0.0 {
        parcel.model().clone()
    } else {
        bank.with_mole_amounts(chemistry.amounts())?
    };
    let mut candidate = parcel.clone();
    let material = candidate.replace_model_with_heat_at_constant_pressure(
        products,
        chemistry.balance().heat_input_j,
        controls.material_energy,
        &mut *cancelled,
    )?;
    let state = candidate.state();
    let total = compensated_sum(chemistry.amounts().iter().map(|amount| amount.value()));
    if !positive(total) {
        return Err(Error::Invalid("unrepresentable final chemical amount"));
    }
    let product_models = candidate.model().thermodynamics().components();
    let product_fractions = candidate
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let correction = chemistry
        .mechanism()
        .thermodynamics()
        .models()
        .iter()
        .zip(chemistry.amounts())
        .map(|(model, amount)| {
            let fraction = product_models
                .binary_search_by(|p| p.species().cmp(model.species()))
                .map_or(0.0, |index| product_fractions[index]);
            (fraction - amount.value() / total).abs()
        })
        .fold(0.0_f64, f64::max);
    let balance = MaterialIsobaricBalance {
        feed_mass_residual_kg: feed_mass - mass,
        chemical_mass_residual_kg: chemistry.mass().value() - mass,
        density_mass_residual_kg: state.state().density * candidate.volume_m3() - mass,
        initial_enthalpy_residual_j: chemistry.target_enthalpy().value()
            - mass * parcel.state().specific_enthalpy_j_kg(),
        freeze_enthalpy_residual_j: mass * state.specific_enthalpy_j_kg()
            - chemistry.enthalpy().value(),
        freeze_internal_energy_residual_j: mass * state.specific_internal_energy_j_kg()
            - chemistry.internal_energy().value(),
        freeze_boundary_work_residual_j: material.energy.work_by_gas_j
            - chemistry.boundary_work().value(),
        relative_initial_volume_correction: chemistry.initial_volume().value() / parcel.volume_m3()
            - 1.0,
        relative_volume_correction: candidate.volume_m3() / chemistry.volume().value() - 1.0,
        relative_temperature_correction: state.state().temperature
            / chemistry.temperature().value()
            - 1.0,
        maximum_mole_fraction_correction: correction,
    };
    if candidate.mass_kg().to_bits() != mass.to_bits()
        || state.state().pressure.to_bits() != parcel.state().state().pressure.to_bits()
        || [
            balance.feed_mass_residual_kg,
            balance.chemical_mass_residual_kg,
            balance.density_mass_residual_kg,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > controls.mass_tolerance_kg)
        || [
            balance.initial_enthalpy_residual_j,
            balance.freeze_enthalpy_residual_j,
            balance.freeze_internal_energy_residual_j,
            balance.freeze_boundary_work_residual_j,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > tolerance)
        || [
            balance.relative_initial_volume_correction,
            balance.relative_volume_correction,
            balance.relative_temperature_correction,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > controls.relative_state_tolerance)
        || !correction.is_finite()
        || correction > 64.0 * f64::EPSILON * chemistry.amounts().len() as f64
    {
        return Err(Error::IsobaricBalance(balance));
    }
    poll(cancelled)?;
    Ok((
        candidate,
        MaterialIsobaricStep {
            chemistry,
            material,
            balance,
        },
    ))
}
