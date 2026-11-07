//! Simultaneous reaction and heat exchange in the original rigid material parcel.
//!
//! The prepared chemistry owner integrates the original species Cv and reaction
//! formation energy. The material owner resolves the products and signed heat
//! together at fixed mass and volume, including transport and acoustic response.
//! No intermediate adiabatic composition or separately heated feed is required.

use super::{
    ClosedGasParcel, MaterialKineticsError, MaterialKineticsSolve, ThermochemicalGasModel,
    admit_bank, admit_sources, compensated_sum, poll, positive,
};
use crate::chemistry::heat::{HeatBath, material::FiniteHeatBody};
use crate::chemistry::{IsochoricKineticsView, IsochoricReactor};
use fs_material::phase::isobaric::IsobaricHeatStep;
use fs_material::thermochemical_gas::energy::parcel::composition::GasCompositionStep;
use fs_qty::{Amount, Temperature, Time, Volume};
use fs_thermochem::kinetics::IdealGasKineticsV1;

type Error = MaterialKineticsError;

/// Cross-owner checks on the actual accepted material state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialIsochoricBalance {
    pub feed_mass_residual_kg: f64,
    pub chemical_mass_residual_kg: f64,
    pub density_mass_residual_kg: f64,
    /// Original chemical feed U minus original material U [J].
    pub initial_internal_energy_residual_j: f64,
    /// Actual final material U minus the chemical endpoint U [J].
    pub freeze_internal_energy_residual_j: f64,
    pub relative_temperature_correction: f64,
    pub relative_pressure_correction: f64,
    pub maximum_mole_fraction_correction: f64,
}

/// Accepted material endpoint and original prepared chemistry report.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialIsochoricStep<'a> {
    /// Borrows the prepared reactor; `to_owned` retains an independent snapshot.
    pub chemistry: IsochoricKineticsView<'a>,
    /// Actual Delta U, signed heat and the material first-law residual.
    pub material: GasCompositionStep,
    pub balance: MaterialIsochoricBalance,
}

/// Gas and finite-body states accepted in one transaction.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialIsochoricBodyStep<'a, 'm> {
    pub gas: MaterialIsochoricStep<'a>,
    pub body: FiniteHeatBody<'m>,
    /// Original body Delta H, Delta U and pressure-volume work.
    pub body_report: IsobaricHeatStep,
    pub body_entropy_change_j_k: f64,
    /// Delta Ugas + Delta Hbody: the gas volume is fixed [J].
    pub combined_heat_residual_j: f64,
    /// Delta Ugas + Delta Ubody + Wbody [J].
    pub combined_energy_residual_j: f64,
}

/// Reusable fixed-volume reactor bound to the complete original material bank.
///
/// Chemical stages, source calorics and the dense feed buffer are reused. Final
/// material resolution still allocates. Both neutral species sources and active
/// feed transport cards must match exactly; no equilibrium projection, nominal
/// heat capacity, additional reaction heat or pressure work is introduced.
pub struct MaterialIsochoricReactor {
    bank: ThermochemicalGasModel,
    reactor: IsochoricReactor,
    controls: MaterialKineticsSolve,
    amounts: Vec<Amount>,
}

impl MaterialIsochoricReactor {
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
        let reactor = IsochoricReactor::new(mechanism, controls.chemistry, &mut cancelled)?;
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
    /// Numerical workspace only; immutable banks and endpoint reports excluded.
    pub fn workspace_bytes(&self) -> usize {
        self.reactor.workspace_bytes() + self.amounts.capacity() * core::mem::size_of::<Amount>()
    }

    /// Conserve original total internal energy at exactly preserved mass and V.
    pub fn advance_adiabatic<'a>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialIsochoricStep<'a>, Error> {
        self.advance(parcel, duration, None, &mut cancelled)
    }

    /// Integrate signed G(Tbath-Tgas) and reaction in the same implicit stages.
    pub fn advance_with_heat_bath<'a>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        bath: HeatBath,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialIsochoricStep<'a>, Error> {
        self.advance(parcel, duration, Some(bath), &mut cancelled)
    }

    fn advance<'a>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        bath: Option<HeatBath>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MaterialIsochoricStep<'a>, Error> {
        let feed_mass = self.prepare(parcel, cancelled)?;
        let volume = Volume::new(parcel.volume_m3());
        let temperature = Temperature::new(parcel.state().state().temperature);
        let chemistry = if let Some(bath) = bath {
            self.reactor.advance_with_heat_bath(
                &self.amounts,
                volume,
                temperature,
                duration,
                bath,
                &mut *cancelled,
            )?
        } else {
            self.reactor.advance_adiabatic(
                &self.amounts,
                volume,
                temperature,
                duration,
                &mut *cancelled,
            )?
        };
        let (candidate, step) = reconcile(
            parcel,
            &self.bank,
            self.controls,
            feed_mass,
            chemistry,
            cancelled,
        )?;
        poll(cancelled)?;
        *parcel = candidate;
        Ok(step)
    }

    /// Exchange heat with an original finite Cp/enthalpy/phase source while
    /// reacting. The gas stays rigid; body expansion work remains in its own
    /// source report. Any source, balance or cancellation refusal publishes
    /// neither endpoint. Zero duration preserves both original states exactly.
    pub fn advance_with_finite_body<'a, 'm>(
        &'a mut self,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        body: FiniteHeatBody<'m>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialIsochoricBodyStep<'a, 'm>, Error> {
        let feed_mass = self.prepare(parcel, &mut cancelled)?;
        let coupled = self.reactor.advance_with_finite_body(
            &self.amounts,
            Volume::new(parcel.volume_m3()),
            Temperature::new(parcel.state().state().temperature),
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
        let combined_heat_residual_j =
            gas.material.stored_energy_change_j + body_report.enthalpy_change_j;
        let combined_energy_residual_j = compensated_sum(
            [
                gas.material.stored_energy_change_j,
                body_report.internal_energy_change_j,
                body_report.boundary_work_j,
            ]
            .into_iter(),
        );
        let tolerance = parcel.mass_kg() * self.controls.material_energy.energy_tolerance_j_kg
            + self.controls.chemistry.energy_tolerance_j;
        if !positive(tolerance)
            || [combined_heat_residual_j, combined_energy_residual_j]
                .iter()
                .any(|value| !value.is_finite() || value.abs() > tolerance)
        {
            return Err(Error::Invalid(
                "rigid material gas and finite body do not close heat/work",
            ));
        }
        let step = MaterialIsochoricBodyStep {
            gas,
            body,
            body_report,
            body_entropy_change_j_k,
            combined_heat_residual_j,
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
        let original = parcel.model().thermodynamics();
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
        let models = self.bank.thermodynamics().components();
        for (source, &fraction) in original
            .components()
            .iter()
            .zip(original.mole_composition().fractions())
        {
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
    chemistry: IsochoricKineticsView<'a>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ClosedGasParcel, MaterialIsochoricStep<'a>), Error> {
    poll(cancelled)?;
    let mass = parcel.mass_kg();
    let tolerance = mass * controls.material_energy.energy_tolerance_j_kg
        + controls.chemistry.energy_tolerance_j;
    if !positive(tolerance) {
        return Err(Error::Invalid(
            "unrepresentable material/chemical energy tolerance",
        ));
    }
    let mut candidate = parcel.clone();
    let material = if chemistry.elapsed_time().value() == 0.0 {
        // Admission already checked the original source. Do not re-solve a
        // physically unchanged state and introduce an inversion roundoff step.
        GasCompositionStep {
            initial: parcel.state().clone(),
            final_state: parcel.state().clone(),
            mass_kg: mass,
            initial_volume_m3: parcel.volume_m3(),
            volume_m3: parcel.volume_m3(),
            supplied_heat_j: 0.0,
            work_by_gas_j: 0.0,
            stored_energy_change_j: 0.0,
            energy_residual_j: 0.0,
            iterations: 0,
        }
    } else {
        candidate.replace_model_with_heat(
            bank.with_mole_amounts(chemistry.amounts())?,
            chemistry.balance().heat_input_j,
            controls.material_energy,
            &mut *cancelled,
        )?
    };
    let state = candidate.state();
    let total = compensated_sum(chemistry.amounts().iter().map(|amount| amount.value()));
    if !positive(total) {
        return Err(Error::Invalid("unrepresentable final chemical amount"));
    }
    let models = candidate.model().thermodynamics().components();
    let fractions = candidate
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
            let fraction = models
                .binary_search_by(|p| p.species().cmp(model.species()))
                .map_or(0.0, |index| fractions[index]);
            (fraction - amount.value() / total).abs()
        })
        .fold(0.0_f64, f64::max);
    let balance = MaterialIsochoricBalance {
        feed_mass_residual_kg: feed_mass - mass,
        chemical_mass_residual_kg: chemistry.mass().value() - mass,
        density_mass_residual_kg: state.state().density * candidate.volume_m3() - mass,
        initial_internal_energy_residual_j: chemistry.target_internal_energy().value()
            - mass * parcel.state().specific_internal_energy_j_kg(),
        freeze_internal_energy_residual_j: mass * state.specific_internal_energy_j_kg()
            - chemistry.internal_energy().value(),
        relative_temperature_correction: state.state().temperature
            / chemistry.temperature().value()
            - 1.0,
        relative_pressure_correction: state.state().pressure / chemistry.pressure().value() - 1.0,
        maximum_mole_fraction_correction: correction,
    };
    if candidate.mass_kg().to_bits() != mass.to_bits()
        || candidate.volume_m3().to_bits() != parcel.volume_m3().to_bits()
        || chemistry.volume().value().to_bits() != parcel.volume_m3().to_bits()
        || [
            balance.feed_mass_residual_kg,
            balance.chemical_mass_residual_kg,
            balance.density_mass_residual_kg,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > controls.mass_tolerance_kg)
        || [
            balance.initial_internal_energy_residual_j,
            balance.freeze_internal_energy_residual_j,
            material.energy_residual_j,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > tolerance)
        || [
            balance.relative_temperature_correction,
            balance.relative_pressure_correction,
        ]
        .iter()
        .any(|value| !value.is_finite() || value.abs() > controls.relative_state_tolerance)
        || !correction.is_finite()
        || correction > 64.0 * f64::EPSILON * chemistry.amounts().len() as f64
    {
        return Err(Error::IsochoricBalance(balance));
    }
    poll(cancelled)?;
    Ok((
        candidate,
        MaterialIsochoricStep {
            chemistry,
            material,
            balance,
        },
    ))
}
