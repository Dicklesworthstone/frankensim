//! Source-bound gas streams on the original projected material/species routes.
//!
//! Incoming density, total specific enthalpy and composition come from one
//! original gas evaluation. The native material and species owners still own
//! every transferred mass. Endpoint temperature is reconstructed from their
//! accepted inventories; no stream heat or pressure-flow work is deposited twice.

use super::{
    Amount, CancelGate, ChtError, GasSpeciesDiffusivity, GasSpeciesEnthalpyField,
    GasSpeciesEnthalpySource, Sum, VoxelDomain, caloric_error, invalid, poll, poll_every,
    source_error,
};
use crate::conjugate::{
    AdvectedMaterial, FlowFace, FlowField, MaterialAdvectionControls, MaterialAdvectionInput,
    MaterialAdvectionStep, MaterialSpeciesControls, MaterialSpeciesInput, MaterialSpeciesStep,
    MaterialVolumeSource,
};
use fs_matdb::QueryPoint;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::QuantitySpec;
use fs_qty::chemistry::SpeciesId;
use fs_qty::semantic::{Composition, CompositionBasis, QuantityKind, SemanticType, ValueForm};

/// One immutable incoming state in the full bank's canonical species order.
///
/// Only the original source constructor can create this handle. Its material
/// state and fractions cannot be edited independently. Exact absent species
/// remain absent; the bank still retains their original laws for later mixing.
#[derive(Clone, Copy, Debug)]
pub struct GasSpeciesAdvectionInlet<'a, const N: usize> {
    source: &'a ThermochemicalGasModel,
    species: &'a [SpeciesId; N],
    context: &'a QueryPoint,
    mass_fractions: [f64; N],
    temperature_k: f64,
    material: AdvectedMaterial,
}

impl<const N: usize> GasSpeciesAdvectionInlet<'_, N> {
    pub fn mass_fractions(&self) -> &[f64; N] {
        &self.mass_fractions
    }
    pub fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    /// The SAME source density and total specific enthalpy used by native inflow.
    pub fn advected_material(&self) -> AdvectedMaterial {
        self.material
    }
}

/// Native volume-source convention with one coherently bound incoming gas.
#[derive(Clone, Copy, Debug, Default)]
pub enum GasSpeciesVolumeSource<'a, const N: usize> {
    #[default]
    None,
    /// Authored signed volume rate [m³/s], positive outward and negative inward.
    /// This changes the native pressure constraint without any local gas mass,
    /// species or enthalpy exchange. Actual face transfers retain their native
    /// donors. It neither infers thermal expansion nor imposes EOS consistency.
    Dilatation { volume_rate_m3_s: f64 },
    Inflow {
        volume_rate_m3_s: f64,
        inlet: GasSpeciesAdvectionInlet<'a, N>,
    },
    /// The original implicit donor supplies composition, density and h.
    Outflow { volume_rate_m3_s: f64 },
}

/// Conserved original inventories and explicitly located incoming states.
#[derive(Clone, Copy, Debug)]
pub struct GasSpeciesAdvectionInput<'i, 'a, const N: usize> {
    pub mass_kg: &'i [f64],
    pub species_mass_kg: &'i [[f64; N]],
    /// Total specific enthalpy in the SAME source/formation convention [J/kg].
    /// Retain accepted h and constituent masses rather than rebuilding them from T.
    pub specific_enthalpy_j_kg: &'i [f64],
    /// One native local volume source per cell; solids must contain None.
    pub volume_sources: &'i [GasSpeciesVolumeSource<'a, N>],
    /// Incoming/backflow states in original x-, x+, y-, y+, z-, z+ order.
    /// Every actually incoming face needs its own state. Outflow uses its donor.
    pub boundary_inflow: [Option<GasSpeciesAdvectionInlet<'a, N>>; 6],
}

/// Original solver budgets plus bounded source closure and adapter payload.
#[derive(Clone, Copy, Debug)]
pub struct ThermochemicalSpeciesAdvectionControls<const N: usize> {
    pub material: MaterialAdvectionControls,
    pub species: MaterialSpeciesControls<N>,
    /// Initial plus endpoint cell caloric inversions. Each uses the source's
    /// original bounded GasEnergySolve. Prebuilt inlet evaluations are separate
    /// single source operations; this transaction reuses their compact caches.
    pub maximum_source_evaluations: usize,
    /// Additional lowering arrays, prepared source field and adapter descriptors.
    /// Native material cells/connections and native species scratch/output retain
    /// their own controls above. Original material evaluation scratch is additional.
    pub maximum_additional_storage_bytes: usize,
    /// Native advective energy defect plus initial and endpoint caloric residuals.
    /// Absolute joules; no formation-energy scale loosens this allowance.
    pub maximum_enthalpy_defect_j: f64,
}

/// Independent source closure alongside the unchanged native transport accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermochemicalSpeciesAdvectionReport {
    pub initial_caloric_defect_j: f64,
    pub endpoint_caloric_defect_j: f64,
    pub absolute_enthalpy_defect_j: f64,
    pub source_evaluations: usize,
}

/// One complete projected carrier, constituent candidate and endpoint source image.
#[derive(Clone, Debug)]
pub struct ThermochemicalSpeciesAdvectionStep<'a, const N: usize> {
    pub carrier: MaterialAdvectionStep,
    pub species: MaterialSpeciesStep<N>,
    /// Reuse the original reconstructed temperature/transport without another
    /// caloric inversion. Changing returned inventories requires fresh preparation.
    pub state: GasSpeciesEnthalpyField<'a, N>,
    pub report: ThermochemicalSpeciesAdvectionReport,
}

impl<'a, const N: usize> GasSpeciesEnthalpySource<'a, N> {
    /// Evaluate a physical incoming temperature and authored canonical mass fractions.
    ///
    /// Pressure and all other source coordinates come from this source context.
    /// The original composition owner requires an exact unit sum; no authored
    /// fraction is clipped, seeded or normalized. One remixture and one complete
    /// material evaluation provide rho, h and the same species partition.
    pub fn advection_inlet(
        &self,
        mass_fractions: [f64; N],
        temperature_k: f64,
        gate: &CancelGate,
    ) -> Result<GasSpeciesAdvectionInlet<'a, N>, ChtError> {
        poll(gate)?;
        self.admit_advection_temperature(temperature_k)?;
        let composition = Composition::new(CompositionBasis::MassFraction, mass_fractions.to_vec())
            .map_err(source_error)?;
        let model = self
            .source
            .with_composition(composition)
            .map_err(source_error)?;
        poll(gate)?;
        let point = self
            .context
            .clone()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                temperature_k,
            )
            .map_err(source_error)?;
        let state = model.evaluate(&point).map_err(source_error)?;
        poll(gate)?;
        Ok(GasSpeciesAdvectionInlet {
            source: self.source,
            species: self.species,
            context: self.context,
            mass_fractions,
            temperature_k,
            material: AdvectedMaterial {
                density_kg_m3: state.state().density,
                specific_enthalpy_j_kg: state.specific_enthalpy_j_kg(),
            },
        })
    }

    fn admit_advection_temperature(&self, temperature_k: f64) -> Result<(), ChtError> {
        if !temperature_k.is_finite()
            || temperature_k < self.controls.minimum_temperature_k
            || temperature_k > self.controls.maximum_temperature_k
        {
            return Err(invalid(
                "gas inlet temperature leaves its admitted original caloric branch",
            ));
        }
        Ok(())
    }

    fn admit_advection_inlet(
        &self,
        inlet: &GasSpeciesAdvectionInlet<'_, N>,
    ) -> Result<(), ChtError> {
        if !std::ptr::eq(inlet.source, self.source)
            || !std::ptr::eq(inlet.context, self.context)
            || inlet.species != self.species
        {
            return Err(invalid(
                "gas advection inlet belongs to a different original source or context",
            ));
        }
        self.admit_advection_temperature(inlet.temperature_k)
    }

    // Initial inventories precede the native carrier transaction, so they
    // cannot use prepare(), which deliberately requires its accepted operator.
    // Reuse the SAME original remixture and pressure-h inversion without a
    // synthetic/no-op flow solve or a second numerical caloric law.
    fn initial_advection_caloric_defect(
        &self,
        domain: &VoxelDomain,
        input: &GasSpeciesAdvectionInput<'_, '_, N>,
        gate: &CancelGate,
    ) -> Result<f64, ChtError> {
        let mut defect = Sum::default();
        for cell in 0..domain.cell_count() {
            poll_every(gate, cell)?;
            if !domain.is_fluid(cell) {
                continue;
            }
            let mass = input.mass_kg[cell];
            let masses = input.species_mass_kg[cell];
            let mut partition = Sum::default();
            partition.add(-mass);
            let mut amounts = [Amount::new(0.0); N];
            for (k, species) in self.source.thermodynamics().components().iter().enumerate() {
                let value = masses[k];
                let fraction = value / mass;
                let amount = value / species.molar_mass().value();
                if !value.is_finite()
                    || value < 0.0
                    || !fraction.is_finite()
                    || !amount.is_finite()
                    || (value > 0.0 && (fraction == 0.0 || amount == 0.0))
                {
                    return Err(invalid(
                        "positive initial gas constituent inventory is unrepresentable",
                    ));
                }
                partition.add(value);
                amounts[k] = Amount::new(amount);
            }
            let relative_partition = partition.value.abs() / mass;
            if !relative_partition.is_finite()
                || relative_partition > 64.0 * N as f64 * f64::EPSILON
            {
                return Err(invalid(
                    "initial gas species must partition carrier mass to machine roundoff",
                ));
            }
            let model = self
                .source
                .with_mole_amounts(&amounts)
                .map_err(source_error)?;
            let solved = model
                .reconstruct_pressure_enthalpy(
                    self.context,
                    input.specific_enthalpy_j_kg[cell],
                    self.controls,
                    || poll(gate).is_err(),
                )
                .map_err(caloric_error)?;
            let residual = mass * solved.residual_j_kg().abs();
            if !residual.is_finite() {
                return Err(invalid(
                    "initial gas advection caloric residual is unrepresentable",
                ));
            }
            defect.add(residual);
        }
        if !defect.value.is_finite() {
            return Err(invalid(
                "initial gas advection caloric residual sum overflowed",
            ));
        }
        poll(gate)?;
        Ok(defect.value)
    }
}

impl FlowField {
    /// Project and advect one source-bound gas over the original native routes.
    ///
    /// Initial M, m_k and total h are authoritative conserved stores. Each inlet
    /// supplies source-consistent rho/h/Y at its declared physical temperature;
    /// local extraction and boundary outflow use the native implicit donor.
    /// Material transport runs once, then species use its retained mass routes
    /// and operator. The final source inversion sees those SAME inventories.
    /// Authored Dilatation enters the volume constraint without any local
    /// material/species source; it still drives actual conservative face fluxes.
    ///
    /// No reaction, Fourier heat, diffusion, latent/formation deposit or pressure
    /// heating is inferred. Stream h already contains pressure-flow work; native
    /// projection work stays mechanical. Source pressure is explicitly fixed.
    /// Endpoint EOS density is diagnostic: this operator split does not impose
    /// an additional EOS volume/pressure closure on the accepted carrier.
    ///
    /// Inputs and the original flow remain unchanged on any refusal/cancellation.
    /// Native work/mass/partition checks and the joint absolute caloric/transport
    /// budget must all pass before any candidate is returned.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn project_and_advect_thermochemical_species<'a, const N: usize>(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: GasSpeciesAdvectionInput<'_, 'a, N>,
        source: &GasSpeciesEnthalpySource<'a, N>,
        controls: ThermochemicalSpeciesAdvectionControls<N>,
        gate: &CancelGate,
    ) -> Result<ThermochemicalSpeciesAdvectionStep<'a, N>, ChtError> {
        poll(gate)?;
        let cells = domain.cell_count();
        if cells > controls.material.max_cells
            || input.mass_kg.len() != cells
            || input.species_mass_kg.len() != cells
            || input.specific_enthalpy_j_kg.len() != cells
            || input.volume_sources.len() != cells
        {
            return Err(invalid(
                "gas advection inputs must cover the admitted native domain",
            ));
        }
        if !controls.maximum_enthalpy_defect_j.is_finite()
            || controls.maximum_enthalpy_defect_j < 0.0
        {
            return Err(invalid(
                "gas advection requires a finite nonnegative absolute joule budget",
            ));
        }
        if !controls.material.maximum_enthalpy_defect_j.is_finite()
            || controls.material.maximum_enthalpy_defect_j < 0.0
        {
            return Err(invalid(
                "native gas advection requires a finite nonnegative joule budget",
            ));
        }
        if let GasSpeciesDiffusivity::Frozen(values) = source.diffusivity
            && values.len() != cells
        {
            return Err(invalid(
                "gas advection endpoint diffusivity must cover the native domain",
            ));
        }
        let mut fluid_cells = 0usize;
        for cell in 0..cells {
            poll_every(gate, cell)?;
            let mass = input.mass_kg[cell];
            let h = input.specific_enthalpy_j_kg[cell];
            if domain.is_fluid(cell) {
                if !mass.is_finite() || mass <= 0.0 || !h.is_finite() {
                    return Err(invalid(
                        "initial gas carrier mass must be positive and enthalpy finite",
                    ));
                }
                fluid_cells += 1;
            } else if mass != 0.0
                || h != 0.0
                || input.species_mass_kg[cell] != [0.0; N]
                || !matches!(input.volume_sources[cell], GasSpeciesVolumeSource::None)
            {
                return Err(invalid(
                    "solid cells cannot carry gas inventories or volume sources",
                ));
            }
        }
        let source_evaluations = fluid_cells
            .checked_mul(2)
            .ok_or_else(|| invalid("gas advection source count overflows"))?;
        if source_evaluations > controls.maximum_source_evaluations {
            return Err(invalid("gas advection source evaluation budget exhausted"));
        }
        let extra_bytes = cells
            .checked_mul(
                core::mem::size_of::<MaterialVolumeSource>()
                    + core::mem::size_of::<Option<[f64; N]>>(),
            )
            .and_then(|n| n.checked_add(source.maximum_prepared_bytes(cells).ok()?))
            .and_then(|n| {
                n.checked_add(
                    core::mem::size_of::<ThermochemicalSpeciesAdvectionStep<'a, N>>()
                        + core::mem::size_of::<GasSpeciesAdvectionInput<'_, 'a, N>>()
                        + core::mem::size_of::<[Amount; N]>()
                        + 6 * core::mem::size_of::<Option<AdvectedMaterial>>()
                        + 6 * core::mem::size_of::<Option<[f64; N]>>()
                        + 1024,
                )
            })
            .filter(|&n| isize::try_from(n).is_ok())
            .ok_or_else(|| invalid("gas advection additional storage extent overflows"))?;
        if extra_bytes > controls.maximum_additional_storage_bytes {
            return Err(invalid("gas advection additional storage budget exhausted"));
        }
        for inlet in input.boundary_inflow.iter().flatten() {
            source.admit_advection_inlet(inlet)?;
        }
        let mut volume_sources = Vec::with_capacity(cells);
        let mut source_fractions = Vec::with_capacity(cells);
        for (cell, entry) in input.volume_sources.iter().enumerate() {
            poll_every(gate, cell)?;
            match *entry {
                GasSpeciesVolumeSource::None => {
                    volume_sources.push(MaterialVolumeSource::None);
                    source_fractions.push(None);
                }
                GasSpeciesVolumeSource::Dilatation { volume_rate_m3_s } => {
                    volume_sources.push(MaterialVolumeSource::Dilatation { volume_rate_m3_s });
                    source_fractions.push(None);
                }
                GasSpeciesVolumeSource::Outflow { volume_rate_m3_s } => {
                    volume_sources.push(MaterialVolumeSource::Outflow { volume_rate_m3_s });
                    source_fractions.push(None);
                }
                GasSpeciesVolumeSource::Inflow {
                    volume_rate_m3_s,
                    inlet,
                } => {
                    source.admit_advection_inlet(&inlet)?;
                    volume_sources.push(MaterialVolumeSource::Inflow {
                        volume_rate_m3_s,
                        material: inlet.material,
                    });
                    source_fractions.push(Some(inlet.mass_fractions));
                }
            }
        }
        let initial_caloric_defect_j =
            source.initial_advection_caloric_defect(domain, &input, gate)?;
        let remaining = controls.maximum_enthalpy_defect_j - initial_caloric_defect_j;
        if remaining < 0.0 {
            return Err(invalid(
                "initial gas advection caloric closure exceeds the joule budget",
            ));
        }
        let mut material_controls = controls.material;
        material_controls.maximum_enthalpy_defect_j =
            material_controls.maximum_enthalpy_defect_j.min(remaining);
        let carrier = self.project_and_advect_material(
            domain,
            faces,
            MaterialAdvectionInput {
                mass_kg: input.mass_kg,
                specific_enthalpy_j_kg: input.specific_enthalpy_j_kg,
                volume_sources: &volume_sources,
                boundary_inflow: input.boundary_inflow.map(|entry| entry.map(|s| s.material)),
            },
            material_controls,
            gate,
        )?;
        let species = carrier.advect_species(
            MaterialSpeciesInput {
                species_mass_kg: input.species_mass_kg,
                source_mass_fractions: &source_fractions,
                boundary_mass_fractions: input
                    .boundary_inflow
                    .map(|entry| entry.map(|s| s.mass_fractions)),
            },
            controls.species,
            gate,
        )?;
        drop(volume_sources);
        drop(source_fractions);
        let state = source.prepare(
            domain,
            &carrier,
            &species.species_mass_kg,
            &carrier.specific_enthalpy_j_kg,
            gate,
        )?;
        let endpoint_caloric_defect_j = state.absolute_caloric_residual_j();
        let mut defect = Sum::default();
        for value in [
            carrier.report.absolute_enthalpy_defect_j,
            initial_caloric_defect_j,
            endpoint_caloric_defect_j,
        ] {
            defect.add(value);
        }
        if !defect.value.is_finite() || defect.value > controls.maximum_enthalpy_defect_j {
            return Err(invalid(
                "gas advection transport and source closure exceed the joule budget",
            ));
        }
        poll(gate)?;
        Ok(ThermochemicalSpeciesAdvectionStep {
            carrier,
            species,
            state,
            report: ThermochemicalSpeciesAdvectionReport {
                initial_caloric_defect_j,
                endpoint_caloric_defect_j,
                absolute_enthalpy_defect_j: defect.value,
                source_evaluations,
            },
        })
    }
}
