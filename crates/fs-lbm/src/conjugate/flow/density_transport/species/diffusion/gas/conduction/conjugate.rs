//! Explicit solid inventories on the same original gas/Fourier heat equation.

use super::{
    CancelGate, ChtError, ConductivityTensor, GasConductionCaloric, GasConductionQueryBudget,
    GasSpeciesEnthalpyField, GasSpeciesEnthalpySource, MaterialAdvectionStep, SolidMaterial,
    TensorCaloricLaw, TensorCaloricStorage, TensorConductivityLaw, ThermalSetup,
    ThermochemicalConductionAdvance, ThermochemicalConductionControls,
    ThermochemicalSpeciesConductionReport, VoxelDomain, admit_interval, admitted_defect,
    advance_native, conduction_report, gas_query_defect, invalid, poll, poll_every, reserved,
};

/// Original caloric and conductivity owners for one stationary solid voxel.
///
/// Bindings must be strictly increasing by cell and cover every solid voxel
/// exactly once. They cannot replace a gas cell. Mass is an explicit conserved
/// reference inventory; neither caloric density nor voxel volume recreates it.
/// Caloric and conductivity queries retain their original bounds, phase
/// restrictions, tensor directions, interval integrals and reference checks.
#[derive(Clone, Copy)]
pub struct ThermochemicalSolidConductionCell<'s> {
    pub cell: usize,
    pub reference_mass_kg: f64,
    /// Previous energy in the original law's declared storage basis, J/kg.
    /// This can be enthalpy, internal energy or an original sensible chart.
    pub specific_energy_j_kg: f64,
    pub caloric: &'s dyn TensorCaloricLaw,
    pub conductivity: &'s dyn TensorConductivityLaw,
    /// Required exactly when ThermalSetup prescribes this solid cell's T.
    /// A latent plateau does not have a unique temperature-to-energy inverse.
    /// The native owner checks this energy against the SAME original law and
    /// authored temperature. Both SDIRK stages use this prescribed state.
    pub prescribed_specific_energy_j_kg: Option<f64>,
}

/// Complete heterogeneous heat and a separately reusable original gas state.
#[derive(Clone, Debug)]
pub struct ThermochemicalConjugateConductionStep<'a, const N: usize> {
    /// Untouched native endpoint(s), face heat and whole-interval ledger.
    /// All-cell specific energies retain each gas/solid law's declared basis;
    /// solid temperatures and latent energies remain physical here.
    pub thermal: ThermochemicalConductionAdvance,
    /// Gas h on fluid cells and exact zero on solids, ready for the next gas
    /// transport or conduction interval. This is not the mixed storage field.
    pub gas_specific_enthalpy_j_kg: Vec<f64>,
    /// Original gas preparation from the accepted fluid h and composition.
    /// Solid temperatures/D are zero and solid transport entries are None.
    pub state: GasSpeciesEnthalpyField<'a, N>,
    /// Native energy closure covers BOTH gas and solid storage. Caloric query
    /// counts and source inversion defects refer to the original gas owner.
    pub report: ThermochemicalSpeciesConductionReport,
}

fn admit_solids(
    domain: &VoxelDomain,
    bindings: &[ThermochemicalSolidConductionCell<'_>],
    gate: &CancelGate,
) -> Result<usize, ChtError> {
    let mut next = 0usize;
    let mut prescribed = 0usize;
    for cell in 0..domain.cell_count() {
        poll_every(gate, cell)?;
        if domain.is_fluid(cell) {
            if bindings
                .get(next)
                .is_some_and(|binding| binding.cell == cell)
            {
                return Err(invalid("a solid caloric binding cannot replace a gas cell"));
            }
            continue;
        }
        let binding = bindings
            .get(next)
            .filter(|binding| binding.cell == cell)
            .ok_or_else(|| {
                invalid("solid caloric bindings must cover every solid in increasing cell order")
            })?;
        if !binding.reference_mass_kg.is_finite()
            || binding.reference_mass_kg <= 0.0
            || !binding.specific_energy_j_kg.is_finite()
            || binding
                .prescribed_specific_energy_j_kg
                .is_some_and(|energy| !energy.is_finite())
        {
            return Err(invalid(
                "solid caloric inventory needs positive finite mass and finite original energy",
            ));
        }
        prescribed += usize::from(binding.prescribed_specific_energy_j_kg.is_some());
        next += 1;
    }
    if next != bindings.len() {
        return Err(invalid(
            "solid caloric bindings contain a duplicate, fluid or out-of-domain cell",
        ));
    }
    Ok(prescribed)
}

impl MaterialAdvectionStep {
    /// Conduct gas and explicit solid inventories through ONE native heat solve.
    ///
    /// Fluid cells use the unchanged original gas T(h,Y), dT/dh=1/Cp, accepted
    /// carrier mass and source-resolved opening isotropic k. Each solid supplies
    /// its own original caloric law, fixed reference mass, previous energy and
    /// tensor conductivity law. Source-backed solid k(T) retains its actual
    /// interval integral; gas k is held for this interval just as in the
    /// gas-only entry. SDIRK2 order applies to the held-gas-k problem, not to an
    /// inferred evolving gas conductivity.
    ///
    /// The native storage coordinate advances m*(e-e0), including zero dT/de
    /// on a sharp latent plateau. The same physical interface heat enters the
    /// gas and solid equations once with opposite signs. No second latent,
    /// formation or contact-heat deposit is applied. Solid occupancy and mass
    /// stay fixed even if an explicitly supplied chart crosses melting; wall
    /// motion, expansion and species transfer are not inferred.
    ///
    /// The native SolidMaterial table validates occupancy and material-pair
    /// contacts; the supplied laws own conductivity and storage. Original
    /// face/source contacts, sources, sinks and thermal boundary rules pass
    /// through the existing tensor owner. All its support and work admissions
    /// remain active. Gas fixed-cell h comes from the original gas forward
    /// query; fixed solids require an explicit prescribed energy, including
    /// any selected latent fraction. The same setup is used at both stages.
    ///
    /// Gas species, input h and frozen D remain zero in solid slots. The result
    /// returns the entire mixed thermal endpoint AND a separately masked gas
    /// h vector; only the latter is passed to original gas preparation. EOS
    /// density never replaces gas carrier mass or the explicit solid mass.
    ///
    /// maximum_caloric_queries counts gas queries only. Solid caloric and
    /// conductivity calls retain the native shared Newton/backtrack/Krylov
    /// ceilings and their own original law admissions. Additional owned fields
    /// use fallible reservations; no separate adapter byte ceiling is claimed.
    /// Any invalid binding, source, solver, heat-budget or cancellation refusal
    /// returns no partial gas or solid endpoint and mutates no input inventory.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn conduct_thermochemical_species_with_solids<'a, const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        gas_specific_enthalpy_j_kg: &[f64],
        solid_materials: &[SolidMaterial],
        solid_cells: &[ThermochemicalSolidConductionCell<'_>],
        setup: &ThermalSetup,
        start_time_s: f64,
        duration_s: f64,
        source: &GasSpeciesEnthalpySource<'a, N>,
        controls: ThermochemicalConductionControls,
        gate: &CancelGate,
    ) -> Result<ThermochemicalConjugateConductionStep<'a, N>, ChtError> {
        poll(gate)?;
        admit_interval(start_time_s, duration_s, &controls)?;
        let declared_prescribed_solids = admit_solids(domain, solid_cells, gate)?;
        domain.check_materials(solid_materials.len())?;
        let cells = domain.cell_count();
        let fluid_cells = cells - solid_cells.len();
        let budget = GasConductionQueryBudget::new(controls.maximum_caloric_queries, gate)?;
        let initial = source.prepare(
            domain,
            self,
            species_mass_kg,
            gas_specific_enthalpy_j_kg,
            gate,
        )?;
        let initial_defect = admitted_defect(
            initial.absolute_caloric_residual_j(),
            controls.maximum_enthalpy_defect_j,
        )?;
        let mut calorics = reserved(fluid_cells)?;
        let mut tensors = reserved(fluid_cells)?;
        let mut mass = reserved(cells)?;
        let mut previous_energy = reserved(cells)?;
        let mut solid = 0usize;
        for cell in 0..cells {
            poll_every(gate, cell)?;
            if domain.is_fluid(cell) {
                calorics.push(GasConductionCaloric::new(
                    source,
                    species_mass_kg[cell],
                    self.mass_kg[cell],
                    &budget,
                )?);
                let k = initial
                    .transport(cell)
                    .ok_or_else(|| invalid("conjugate heat lost its original fluid transport"))?
                    .thermal_conductivity_w_m_k();
                tensors.push(ConductivityTensor::new([
                    [k, 0.0, 0.0],
                    [0.0, k, 0.0],
                    [0.0, 0.0, k],
                ])?);
                mass.push(self.mass_kg[cell]);
                previous_energy.push(gas_specific_enthalpy_j_kg[cell]);
            } else {
                let binding = &solid_cells[solid];
                mass.push(binding.reference_mass_kg);
                previous_energy.push(binding.specific_energy_j_kg);
                solid += 1;
            }
        }
        drop(initial);

        let mut prescribed = reserved(setup.fixed_temperature.len())?;
        let mut fixed_solids = 0usize;
        for (index, &(cell, temperature)) in setup.fixed_temperature.iter().enumerate() {
            poll_every(gate, index)?;
            if cell >= cells {
                return Err(invalid(
                    "prescribed conjugate temperature cell is outside its domain",
                ));
            }
            let energy = match solid_cells.binary_search_by_key(&cell, |binding| binding.cell) {
                Ok(index) => {
                    fixed_solids += 1;
                    solid_cells[index]
                        .prescribed_specific_energy_j_kg
                        .ok_or_else(|| {
                            invalid("a fixed solid temperature requires explicit original energy")
                        })?
                }
                Err(solids_before) => {
                    calorics[cell - solids_before].specific_enthalpy_at_temperature(temperature)?
                }
            };
            prescribed.push((cell, energy));
        }
        if fixed_solids != declared_prescribed_solids {
            return Err(invalid(
                "prescribed solid energies must exactly cover the fixed solid temperatures",
            ));
        }
        // Original geometry admission also rejects duplicate fixed-cell owners;
        // with unique entries the count above excludes every unused solid e.
        let mut thermal_laws: Vec<&dyn TensorCaloricLaw> = reserved(cells)?;
        let mut conductivity_laws: Vec<&dyn TensorConductivityLaw> = reserved(cells)?;
        let (mut fluid, mut solid) = (0usize, 0usize);
        for cell in 0..cells {
            poll_every(gate, cell)?;
            if domain.is_fluid(cell) {
                thermal_laws.push(&calorics[fluid]);
                conductivity_laws.push(&tensors[fluid]);
                fluid += 1;
            } else {
                thermal_laws.push(solid_cells[solid].caloric);
                conductivity_laws.push(solid_cells[solid].conductivity);
                solid += 1;
            }
        }
        let thermal = advance_native(
            domain,
            solid_materials,
            setup,
            &conductivity_laws,
            TensorCaloricStorage {
                laws: &thermal_laws,
                reference_mass_kg: &mass,
                previous_specific_energy_j_kg: &previous_energy,
            },
            &prescribed,
            start_time_s,
            duration_s,
            &controls,
            gate,
        )?;
        let maximum_query_caloric_defect_j =
            gas_query_defect(&calorics, controls.maximum_enthalpy_defect_j, gate)?;
        let caloric_queries = budget.queries();
        drop(thermal_laws);
        drop(conductivity_laws);
        drop(calorics);
        drop(tensors);
        drop(mass);
        drop(prescribed);

        // Reuse the initial combined-energy allocation for the gas-only output;
        // the native result retains every physical solid energy unchanged.
        for (cell, energy) in previous_energy.iter_mut().enumerate() {
            poll_every(gate, cell)?;
            *energy = if domain.is_fluid(cell) {
                thermal.endpoint().specific_energy_j_kg[cell]
            } else {
                0.0
            };
        }
        let gas_specific_enthalpy_j_kg = previous_energy;
        let state = source.prepare(
            domain,
            self,
            species_mass_kg,
            &gas_specific_enthalpy_j_kg,
            gate,
        )?;
        let report = conduction_report(
            &thermal,
            initial_defect,
            state.absolute_caloric_residual_j(),
            maximum_query_caloric_defect_j,
            caloric_queries,
            &controls,
        )?;
        poll(gate)?;
        Ok(ThermochemicalConjugateConductionStep {
            thermal,
            gas_specific_enthalpy_j_kg,
            state,
            report,
        })
    }
}
