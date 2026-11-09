//! Source-owned ideal-gas calorics on the native constituent diffusion owner.
//!
//! The immutable source bank fixes canonical species, NASA-9 formation
//! references and transport sources. Every proposed composition and total h
//! resolve together through the original material inversion. Pressure remains
//! the explicitly supplied query coordinate: source density is returned for a
//! subsequent flow/pressure solve, never substituted for the fixed carrier.
//! This is an isobaric constitutive closure, not compressible pressure work.

use super::{
    CancelGate, ChtError, MAX_MATERIAL_SPECIES, MaterialAdvectionStep, SpeciesDiffusionControls,
    SpeciesEnthalpyDiffusionStep, Sum, VoxelDomain, poll, poll_every,
};
use fs_matdb::QueryPoint;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::enthalpy::GasEnthalpyLaw;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_qty::chemistry::SpeciesId;
use fs_qty::{Amount, Temperature};

mod adaptive;
mod advection;
pub use advection::{
    GasSpeciesAdvectionInlet, GasSpeciesAdvectionInput, GasSpeciesVolumeSource,
    ThermochemicalSpeciesAdvectionControls, ThermochemicalSpeciesAdvectionReport,
    ThermochemicalSpeciesAdvectionStep,
};
mod boundary;
pub use boundary::{BoundaryThermochemicalSpeciesDiffusionStep, GasSpeciesEnthalpyBoundarySource};
#[cfg(feature = "enthalpy-coupling")]
mod conduction;
#[cfg(feature = "enthalpy-coupling")]
pub use conduction::{
    ThermochemicalConductionAdvance, ThermochemicalConductionControls,
    ThermochemicalConductionScheme, ThermochemicalSpeciesConductionReport,
    ThermochemicalSpeciesConductionStep,
    ThermochemicalConjugateConductionStep, ThermochemicalSolidConductionCell,
    ThermochemicalConductivityModel,
};
mod diffusivity;
pub use diffusivity::GasSpeciesDiffusivity;

fn invalid(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.thermochemical_species",
        reason: reason.into(),
    }
}

fn source_error(error: impl core::fmt::Display) -> ChtError {
    invalid(format!("original gas source refused: {error}"))
}

fn caloric_error(error: GasEnergyError) -> ChtError {
    match error {
        GasEnergyError::Cancelled => ChtError::Cancelled,
        other => source_error(other),
    }
}

/// One immutable full species bank, pressure/source context and diffusion law.
///
/// The explicit identity array must match the source bank's canonical order.
/// A species may be absent from individual cells without disappearing from the
/// bank. No coefficient is selected by a guessed species name or index.
#[derive(Clone, Copy, Debug)]
pub struct GasSpeciesEnthalpySource<'a, const N: usize> {
    source: &'a ThermochemicalGasModel,
    species: &'a [SpeciesId; N],
    context: &'a QueryPoint,
    diffusivity: GasSpeciesDiffusivity<'a>,
    controls: GasEnergySolve,
}

/// Compact source-bound state prepared from one complete mass/h candidate.
///
/// Original source identities remain borrowed from the full bank; numerical
/// coefficients come from its complete resolved material state. Per-cell
/// evaluation receipts are not duplicated in this field. Source evaluation
/// scratch and the immutable bank are additional to its owned array payload.
#[derive(Clone, Debug)]
pub struct GasSpeciesEnthalpyField<'a, const N: usize> {
    source: &'a ThermochemicalGasModel,
    species: &'a [SpeciesId; N],
    context: &'a QueryPoint,
    temperatures_k: Vec<f64>,
    diffusivity_m2_s: Vec<f64>,
    transport: Vec<Option<FluidTransportState>>,
    caloric_residuals_j_kg: Vec<f64>,
    absolute_caloric_residual_j: f64,
}

impl<'a, const N: usize> GasSpeciesEnthalpyField<'a, N> {
    /// Full original bank, including sources for currently absent species.
    pub fn source(&self) -> &'a ThermochemicalGasModel {
        self.source
    }
    /// Canonical identity attached to every constituent array slot.
    pub fn species(&self) -> &'a [SpeciesId; N] {
        self.species
    }
    /// Original pressure and auxiliary material-source coordinates.
    pub fn context(&self) -> &'a QueryPoint {
        self.context
    }
    /// Reconstructed fluid temperatures; solid entries are zero.
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    /// Supplied or original source-resolved common diffusivity; solids are zero.
    pub fn diffusivity_m2_s(&self) -> &[f64] {
        &self.diffusivity_m2_s
    }
    /// Cp, viscosity, conductivity and EOS density at the SAME resolved state.
    /// Density does not replace the native carrier mass divided by its volume.
    pub fn transport(&self, cell: usize) -> Option<&FluidTransportState> {
        self.transport.get(cell).and_then(Option::as_ref)
    }
    /// Actual source h(T,Y,p) minus supplied total h; solid entries are zero.
    pub fn caloric_residuals_j_kg(&self) -> &[f64] {
        &self.caloric_residuals_j_kg
    }
    /// Sum of carrier mass times absolute inversion residual, in joules.
    pub fn absolute_caloric_residual_j(&self) -> f64 {
        self.absolute_caloric_residual_j
    }
}

impl<'a, const N: usize> GasSpeciesEnthalpySource<'a, N> {
    /// Bind an explicitly supplied common diffusivity field in native cell order.
    /// This convenience entry neither infers molecular D nor selects a Schmidt number.
    pub fn new(
        source: &'a ThermochemicalGasModel,
        species: &'a [SpeciesId; N],
        context: &'a QueryPoint,
        diffusivity_m2_s: &'a [f64],
        controls: GasEnergySolve,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        Self::new_with_diffusivity(
            source,
            species,
            context,
            GasSpeciesDiffusivity::Frozen(diffusivity_m2_s),
            controls,
            gate,
        )
    }

    /// Bind a full original gas bank and an explicit diffusion constitutive owner.
    ///
    /// The original GasEnthalpyLaw admission checks a continuous increasing
    /// NASA-9 caloric branch for EVERY bank species once, including dormant
    /// species. Raw discontinuous region seams are not silently crossed.
    pub fn new_with_diffusivity(
        source: &'a ThermochemicalGasModel,
        species: &'a [SpeciesId; N],
        context: &'a QueryPoint,
        diffusivity: GasSpeciesDiffusivity<'a>,
        controls: GasEnergySolve,
        gate: &CancelGate,
    ) -> Result<Self, ChtError> {
        poll(gate)?;
        let components = source.thermodynamics().components();
        if N == 0
            || N > MAX_MATERIAL_SPECIES
            || components.len() != N
            || species
                .iter()
                .zip(components)
                .any(|(id, model)| id != model.species())
        {
            return Err(invalid(
                "explicit species identities must cover the full source bank in canonical order",
            ));
        }
        diffusivity.validate_species(species)?;
        GasEnthalpyLaw::from_mole(
            source,
            source.thermodynamics().mole_composition().clone(),
            context.clone(),
            controls,
        )
        .map_err(caloric_error)?;
        poll(gate)?;
        Ok(Self {
            source,
            species,
            context,
            diffusivity,
            controls,
        })
    }

    /// Full bank retained by this closure, including dormant species sources.
    pub fn source(&self) -> &'a ThermochemicalGasModel {
        self.source
    }
    /// Conservative bound on owned prepared arrays and their descriptors.
    /// Immutable source banks and temporary original material inversion scratch
    /// are additional; no per-cell source card clones survive preparation.
    pub fn maximum_prepared_bytes(&self, cells: usize) -> Result<usize, ChtError> {
        cells
            .checked_mul(
                3 * core::mem::size_of::<f64>()
                    + core::mem::size_of::<Option<FluidTransportState>>(),
            )
            .and_then(|bytes| {
                bytes.checked_add(core::mem::size_of::<GasSpeciesEnthalpyField<'a, N>>())
            })
            .ok_or_else(|| invalid("prepared gas field storage extent overflows"))
    }

    /// Reconstruct temperature and transport from the SAME species masses and h.
    ///
    /// Only machine-roundoff partition discrepancies are admitted before
    /// converting actual m_k/M_k amounts through the original mixture owner.
    /// A generous transport partition budget never authorizes material-state
    /// renormalization. Positive species that lose their molar amount refuse.
    /// All inputs and the accepted carrier remain immutable on any failure.
    pub fn prepare(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<GasSpeciesEnthalpyField<'a, N>, ChtError> {
        poll(gate)?;
        carrier.validate_diffusion_carrier(domain, gate)?;
        let count = domain.cell_count();
        self.maximum_prepared_bytes(count)?;
        if species_mass_kg.len() != count || specific_enthalpy_j_kg.len() != count {
            return Err(invalid(
                "gas inventories must cover the complete native diffusion domain",
            ));
        }
        if let GasSpeciesDiffusivity::Frozen(values) = self.diffusivity
            && values.len() != count
        {
            return Err(invalid(
                "supplied gas diffusivity must cover the complete native domain",
            ));
        }
        let mut temperatures_k = Vec::with_capacity(count);
        let mut diffusivity_m2_s = Vec::with_capacity(count);
        let mut transport = Vec::with_capacity(count);
        let mut caloric_residuals_j_kg = Vec::with_capacity(count);
        let mut absolute_residual = Sum::default();
        for cell in 0..count {
            poll_every(gate, cell)?;
            let masses = &species_mass_kg[cell];
            let h = specific_enthalpy_j_kg[cell];
            if !h.is_finite() || masses.iter().any(|m| !m.is_finite() || *m < 0.0) {
                return Err(invalid(
                    "gas constituent inventories and total enthalpy must be finite",
                ));
            }
            if !domain.is_fluid(cell) {
                if *masses != [0.0; N] || h != 0.0 {
                    return Err(invalid("solid cells cannot carry gas species or enthalpy"));
                }
                if let GasSpeciesDiffusivity::Frozen(values) = self.diffusivity
                    && values[cell] != 0.0
                {
                    return Err(invalid("solid cells cannot carry gas diffusivity"));
                }
                temperatures_k.push(0.0);
                diffusivity_m2_s.push(0.0);
                transport.push(None);
                caloric_residuals_j_kg.push(0.0);
                continue;
            }
            let carrier_mass = carrier.mass_kg[cell];
            let mut partition = Sum::default();
            partition.add(-carrier_mass);
            let mut amounts = [Amount::new(0.0); N];
            for (k, model) in self.source.thermodynamics().components().iter().enumerate() {
                partition.add(masses[k]);
                let amount = masses[k] / model.molar_mass().value();
                if !amount.is_finite() || (masses[k] > 0.0 && amount == 0.0) {
                    return Err(invalid(
                        "positive gas constituent molar amount is unrepresentable",
                    ));
                }
                amounts[k] = Amount::new(amount);
            }
            let relative_partition = partition.value.abs() / carrier_mass;
            if !relative_partition.is_finite()
                || relative_partition > 64.0 * N as f64 * f64::EPSILON
            {
                return Err(invalid(
                    "gas species must partition carrier mass to machine roundoff",
                ));
            }
            let model = self
                .source
                .with_mole_amounts(&amounts)
                .map_err(source_error)?;
            let solution = model
                .reconstruct_pressure_enthalpy(self.context, h, self.controls, || {
                    poll(gate).is_err()
                })
                .map_err(caloric_error)?;
            let state = solution.state();
            let temperature = state.state().temperature;
            let diffusivity =
                self.diffusivity
                    .sample(cell, temperature, state.state().pressure, gate)?;
            if !diffusivity.is_finite() || diffusivity < 0.0 {
                return Err(invalid("source diffusivity must be finite and nonnegative"));
            }
            let residual = solution.residual_j_kg();
            let energy_residual = carrier_mass * residual.abs();
            if !energy_residual.is_finite() {
                return Err(invalid("gas caloric residual is unrepresentable in joules"));
            }
            absolute_residual.add(energy_residual);
            temperatures_k.push(temperature);
            diffusivity_m2_s.push(diffusivity);
            transport.push(Some(state.transport_state().map_err(source_error)?));
            caloric_residuals_j_kg.push(residual);
        }
        if !absolute_residual.value.is_finite() {
            return Err(invalid("gas caloric residual sum overflowed"));
        }
        poll(gate)?;
        Ok(GasSpeciesEnthalpyField {
            source: self.source,
            species: self.species,
            context: self.context,
            temperatures_k,
            diffusivity_m2_s,
            transport,
            caloric_residuals_j_kg,
            absolute_caloric_residual_j: absolute_residual.value,
        })
    }

    /// Original partial SPECIFIC enthalpies at the geometric face temperature.
    ///
    /// Linear center-to-face interpolation uses each cell's physical width, so
    /// graded cells do not turn this into an unweighted temperature average.
    /// Ideal-gas partial h_k depends on T; no mixture Cp or composition-averaged
    /// enthalpy substitutes for the original species law.
    pub fn face_partials(
        &self,
        field: &GasSpeciesEnthalpyField<'a, N>,
        domain: &VoxelDomain,
        cells: [usize; 2],
        axis: usize,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        poll(gate)?;
        if !std::ptr::eq(field.source, self.source)
            || !std::ptr::eq(field.context, self.context)
            || field.species != self.species
            || axis >= 3
            || field.temperatures_k.len() != domain.cell_count()
            || cells.iter().any(|&c| c >= domain.cell_count())
        {
            return Err(invalid(
                "gas face query does not match its prepared source/domain",
            ));
        }
        let [left, right] = cells;
        if !domain.is_fluid(left)
            || !domain.is_fluid(right)
            || domain.neighbor(left, 2 * axis + 1) != Some(right)
        {
            return Err(invalid(
                "gas partial enthalpy requires an ordered adjacent fluid face",
            ));
        }
        let wl = domain.widths(left)[axis];
        let wr = domain.widths(right)[axis];
        let scale = wl.max(wr);
        let right_weight = (wl / scale) / (wl / scale + wr / scale);
        let tl = field.temperatures_k[left];
        let tr = field.temperatures_k[right];
        let temperature = right_weight.mul_add(tr - tl, tl);
        self.partials_at_temperature(temperature, gate)
    }

    // Interior and prescribed physical boundary faces share the original
    // canonical source evaluation. No donor-temperature or mixture-h closure.
    fn partials_at_temperature(
        &self,
        temperature: f64,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        poll(gate)?;
        if !temperature.is_finite()
            || temperature < self.controls.minimum_temperature_k
            || temperature > self.controls.maximum_temperature_k
        {
            return Err(invalid(
                "gas face temperature leaves the admitted source branch",
            ));
        }
        let mut partial_h = [0.0; N];
        for (k, model) in self.source.thermodynamics().components().iter().enumerate() {
            let h = model
                .evaluate_properties(Temperature::new(temperature))
                .map_err(source_error)?
                .h()
                .value()
                / model.molar_mass().value();
            if !h.is_finite() {
                return Err(invalid(
                    "source partial specific enthalpy is unrepresentable",
                ));
            }
            partial_h[k] = h;
        }
        poll(gate)?;
        Ok(partial_h)
    }
}

/// Source closure defects remain separate from the original conservative audit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermochemicalSpeciesDiffusionReport {
    pub initial_caloric_defect_j: f64,
    pub endpoint_caloric_defect_j: f64,
    /// Original transport defect plus both absolute source-inversion defects.
    pub absolute_enthalpy_defect_j: f64,
}

/// Native species/enthalpy endpoint and its refreshed original gas coefficients.
#[derive(Clone, Debug)]
pub struct ThermochemicalSpeciesDiffusionStep<'a, const N: usize> {
    pub step: SpeciesEnthalpyDiffusionStep<N>,
    pub state: GasSpeciesEnthalpyField<'a, N>,
    pub report: ThermochemicalSpeciesDiffusionReport,
}

impl MaterialAdvectionStep {
    /// Advance one frozen-coefficient interval and resolve its original gas state.
    ///
    /// Face h_k and D come from the initial source state for this interval.
    /// The endpoint re-evaluates temperature, composition and all transport
    /// properties together; repeated/adaptive callers must refresh the next
    /// interval from that state. No second heat deposit, reaction heat, pressure
    /// update or coupled second-order temporal accuracy is implied.
    #[allow(clippy::too_many_arguments)]
    pub fn diffuse_thermochemical_species<'a, const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        duration_s: f64,
        source: &GasSpeciesEnthalpySource<'a, N>,
        mut controls: SpeciesDiffusionControls<N>,
        maximum_enthalpy_defect_j: f64,
        gate: &CancelGate,
    ) -> Result<ThermochemicalSpeciesDiffusionStep<'a, N>, ChtError> {
        poll(gate)?;
        if !maximum_enthalpy_defect_j.is_finite() || maximum_enthalpy_defect_j < 0.0 {
            return Err(invalid(
                "finite nonnegative gas enthalpy defect budget required",
            ));
        }
        let prepared_bytes = source.maximum_prepared_bytes(domain.cell_count())?;
        controls.maximum_storage_bytes = controls
            .maximum_storage_bytes
            .checked_sub(prepared_bytes)
            .ok_or_else(|| invalid("gas preparation exceeds the diffusion storage budget"))?;
        let initial =
            source.prepare(domain, self, species_mass_kg, specific_enthalpy_j_kg, gate)?;
        let initial_caloric_defect_j = initial.absolute_caloric_residual_j();
        let remaining = maximum_enthalpy_defect_j - initial_caloric_defect_j;
        if remaining < 0.0 {
            return Err(invalid(
                "initial gas caloric closure exceeds the enthalpy budget",
            ));
        }
        let step = self.diffuse_species_with_enthalpy(
            domain,
            species_mass_kg,
            specific_enthalpy_j_kg,
            initial.diffusivity_m2_s(),
            duration_s,
            |cells, axis| source.face_partials(&initial, domain, cells, axis, gate),
            controls,
            remaining,
            gate,
        )?;
        drop(initial);
        let state = source.prepare(
            domain,
            self,
            &step.species.species_mass_kg,
            &step.specific_enthalpy_j_kg,
            gate,
        )?;
        let endpoint_caloric_defect_j = state.absolute_caloric_residual_j();
        let mut defect = Sum::default();
        for value in [
            step.report.absolute_enthalpy_defect_j,
            initial_caloric_defect_j,
            endpoint_caloric_defect_j,
        ] {
            defect.add(value);
        }
        if !defect.value.is_finite() || defect.value > maximum_enthalpy_defect_j {
            return Err(invalid(
                "gas source closure and conserved transport exceed the enthalpy budget",
            ));
        }
        poll(gate)?;
        Ok(ThermochemicalSpeciesDiffusionStep {
            step,
            state,
            report: ThermochemicalSpeciesDiffusionReport {
                initial_caloric_defect_j,
                endpoint_caloric_defect_j,
                absolute_enthalpy_defect_j: defect.value,
            },
        })
    }
}
