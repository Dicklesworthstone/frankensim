//! Original mixture enthalpy on the existing implicit voxel Fourier owner.
//!
//! Composition and carrier inventory are fixed for this heat interval. The
//! conductivity follows an explicit held or original-integral choice, while
//! the original nonlinear caloric inverse drives the native heat solve.

use super::{
    CancelGate, ChtError, GasSpeciesEnthalpyField, GasSpeciesEnthalpySource, MaterialAdvectionStep,
    Sum, VoxelDomain, invalid, poll, poll_every,
};
use crate::conjugate::{
    ConductivityTensor, RadiationConfig, SolidMaterial, TensorCaloricEndpoint, TensorCaloricLaw,
    TensorCaloricRadiationReport, TensorCaloricSdirk2Step, TensorCaloricStep, TensorCaloricStorage,
    TensorConductivityLaw, TensorNonlinearConfig, ThermalSetup, TransientRecord,
    advance_caloric_tensor_conduction, advance_caloric_tensor_conduction_radiating,
    advance_caloric_tensor_conduction_radiating_sdirk2, advance_caloric_tensor_conduction_sdirk2,
};

mod conductivity;
pub use conductivity::ThermochemicalConductivityModel;
use conductivity::{GasConductionConductivity, GasConductivityBudget};

mod conjugate;
pub use conjugate::{ThermochemicalConjugateConductionStep, ThermochemicalSolidConductionCell};
mod source;
use source::{GasConductionCaloric, GasConductionQueryBudget};

/// Native caloric time formula with the explicitly selected conductivity law.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermochemicalConductionScheme {
    BackwardEuler,
    Sdirk2,
}

/// Original heat work limits plus independent caloric work and absolute joules.
#[derive(Clone, Copy, Debug)]
pub struct ThermochemicalConductionControls {
    pub scheme: ThermochemicalConductionScheme,
    /// Original source integration or an explicitly held opening gas k.
    pub conductivity: ThermochemicalConductivityModel,
    /// The native Newton, Krylov and backtracking ceilings apply to the entire
    /// heat interval, including both stages when SDIRK2 is selected.
    pub thermal: TensorNonlinearConfig,
    /// Shared original gas-source caloric query allowance across every fluid
    /// cell and native trial. Borrowed solid queries retain native work limits;
    /// the gas source inverse also retains its own bounded root iterations.
    pub maximum_caloric_queries: usize,
    /// Native interval energy closure plus initial and final source inversion
    /// defects must fit this absolute allowance. The sum of each cell's largest
    /// trial caloric residual must independently fit it as well.
    pub maximum_enthalpy_defect_j: f64,
}

/// Unmodified native heat solution, including actual face fluxes and receipts.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermochemicalConductionAdvance {
    BackwardEuler(TensorCaloricStep),
    Sdirk2(TensorCaloricSdirk2Step),
    RadiatingBackwardEuler {
        step: TensorCaloricStep,
        radiation: TensorCaloricRadiationReport,
    },
    RadiatingSdirk2 {
        step: TensorCaloricSdirk2Step,
        radiation: [TensorCaloricRadiationReport; 2],
    },
}

impl ThermochemicalConductionAdvance {
    pub fn endpoint(&self) -> &TensorCaloricEndpoint {
        match self {
            Self::BackwardEuler(step) => &step.endpoint,
            Self::Sdirk2(step) => &step.stages[1],
            Self::RadiatingBackwardEuler { step, .. } => &step.endpoint,
            Self::RadiatingSdirk2 { step, .. } => &step.stages[1],
        }
    }

    /// Complete physical interval; SDIRK stage heat is counted by its owner.
    pub fn record(&self) -> &TransientRecord {
        match self {
            Self::BackwardEuler(step) => &step.record,
            Self::Sdirk2(step) => &step.record,
            Self::RadiatingBackwardEuler { step, .. } => &step.record,
            Self::RadiatingSdirk2 { step, .. } => &step.record,
        }
    }

    /// Original radiation iteration evidence; empty for nonradiating heat.
    pub fn radiation_reports(&self) -> &[TensorCaloricRadiationReport] {
        match self {
            Self::BackwardEuler(_) | Self::Sdirk2(_) => &[],
            Self::RadiatingBackwardEuler { radiation, .. } => core::slice::from_ref(radiation),
            Self::RadiatingSdirk2 { radiation, .. } => radiation,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermochemicalSpeciesConductionReport {
    pub caloric_queries: usize,
    /// Actual source conductivity point queries across all cells and trials.
    /// Zero when the opening conductivity is held.
    pub conductivity_evaluations: usize,
    /// Largest estimated absolute interval-mean k error [W/(m K)], including
    /// discarded trials. This is a quadrature estimate, not a rigorous bound
    /// or an additional physical heat/energy defect.
    pub maximum_estimated_conductivity_error_w_m_k: f64,
    pub initial_caloric_defect_j: f64,
    pub endpoint_caloric_defect_j: f64,
    /// Sum over cells of their largest actual root residual times carrier mass.
    /// Includes discarded Newton states; it is a source-accuracy admission,
    /// not another physical heat deposit or accumulated interval energy.
    pub maximum_query_caloric_defect_j: f64,
    /// Absolute closure of the original native Fourier interval ledger.
    /// Local cell/face equations also retain the original thermal controls.
    pub native_energy_defect_j: f64,
    pub absolute_enthalpy_defect_j: f64,
}

#[derive(Clone, Debug)]
pub struct ThermochemicalSpeciesConductionStep<'a, const N: usize> {
    pub thermal: ThermochemicalConductionAdvance,
    /// Reusable original gas image at the returned composition and total h.
    /// Its transport includes refreshed k and D for the next physical stage.
    pub state: GasSpeciesEnthalpyField<'a, N>,
    pub report: ThermochemicalSpeciesConductionReport,
}

impl<const N: usize> ThermochemicalSpeciesConductionStep<'_, N> {
    pub fn specific_enthalpy_j_kg(&self) -> &[f64] {
        &self.thermal.endpoint().specific_energy_j_kg
    }
}

fn reserved<T>(count: usize) -> Result<Vec<T>, ChtError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| invalid("gas conduction preparation allocation refused"))?;
    Ok(values)
}

fn admitted_defect(value: f64, allowance: f64) -> Result<f64, ChtError> {
    if !value.is_finite() || value < 0.0 || value > allowance {
        Err(invalid(
            "gas conduction exceeds its absolute caloric/energy allowance",
        ))
    } else {
        Ok(value)
    }
}

fn admit_interval(
    start_time_s: f64,
    duration_s: f64,
    controls: &ThermochemicalConductionControls,
) -> Result<(), ChtError> {
    if !controls.maximum_enthalpy_defect_j.is_finite()
        || controls.maximum_enthalpy_defect_j < 0.0
        || !start_time_s.is_finite()
        || start_time_s < 0.0
        || !duration_s.is_finite()
        || duration_s <= 0.0
        || !(start_time_s + duration_s).is_finite()
        || start_time_s + duration_s <= start_time_s
    {
        return Err(invalid(
            "gas conduction needs resolved physical time and absolute joule controls",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn advance_native(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed: &[(usize, f64)],
    start_time_s: f64,
    duration_s: f64,
    controls: &ThermochemicalConductionControls,
    radiation: Option<&RadiationConfig>,
    gate: &CancelGate,
) -> Result<ThermochemicalConductionAdvance, ChtError> {
    if let Some(radiation) = radiation {
        return match controls.scheme {
            ThermochemicalConductionScheme::BackwardEuler => {
                let (step, radiation) = advance_caloric_tensor_conduction_radiating(
                    domain,
                    solids,
                    setup,
                    conductivity_laws,
                    storage,
                    prescribed,
                    start_time_s,
                    duration_s,
                    &controls.thermal,
                    radiation,
                    gate,
                )?;
                Ok(ThermochemicalConductionAdvance::RadiatingBackwardEuler { step, radiation })
            }
            ThermochemicalConductionScheme::Sdirk2 => {
                let (step, radiation) = advance_caloric_tensor_conduction_radiating_sdirk2(
                    domain,
                    solids,
                    [setup, setup],
                    conductivity_laws,
                    storage,
                    [prescribed, prescribed],
                    start_time_s,
                    duration_s,
                    &controls.thermal,
                    radiation,
                    gate,
                )?;
                Ok(ThermochemicalConductionAdvance::RadiatingSdirk2 { step, radiation })
            }
        };
    }
    match controls.scheme {
        ThermochemicalConductionScheme::BackwardEuler => Ok(
            ThermochemicalConductionAdvance::BackwardEuler(advance_caloric_tensor_conduction(
                domain,
                solids,
                setup,
                conductivity_laws,
                storage,
                prescribed,
                start_time_s,
                duration_s,
                &controls.thermal,
                gate,
            )?),
        ),
        ThermochemicalConductionScheme::Sdirk2 => Ok(ThermochemicalConductionAdvance::Sdirk2(
            advance_caloric_tensor_conduction_sdirk2(
                domain,
                solids,
                [setup, setup],
                conductivity_laws,
                storage,
                [prescribed, prescribed],
                start_time_s,
                duration_s,
                &controls.thermal,
                gate,
            )?,
        )),
    }
}

fn gas_query_defect(
    calorics: &[GasConductionCaloric<'_, '_>],
    allowance: f64,
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    let mut defect = Sum::default();
    for (cell, law) in calorics.iter().enumerate() {
        poll_every(gate, cell)?;
        defect.add(law.maximum_query_caloric_defect_j());
    }
    admitted_defect(defect.value, allowance)
}

fn conduction_report(
    thermal: &ThermochemicalConductionAdvance,
    initial_caloric_defect_j: f64,
    endpoint_caloric_defect_j: f64,
    maximum_query_caloric_defect_j: f64,
    caloric_queries: usize,
    conductivity_budget: &GasConductivityBudget<'_>,
    controls: &ThermochemicalConductionControls,
) -> Result<ThermochemicalSpeciesConductionReport, ChtError> {
    let endpoint_caloric_defect_j = admitted_defect(
        endpoint_caloric_defect_j,
        controls.maximum_enthalpy_defect_j,
    )?;
    let native_energy_defect_j = thermal.record().closure_j.abs();
    let mut defect = Sum::default();
    for value in [
        native_energy_defect_j,
        initial_caloric_defect_j,
        endpoint_caloric_defect_j,
    ] {
        defect.add(value);
    }
    let absolute_enthalpy_defect_j =
        admitted_defect(defect.value, controls.maximum_enthalpy_defect_j)?;
    Ok(ThermochemicalSpeciesConductionReport {
        caloric_queries,
        conductivity_evaluations: conductivity_budget.evaluations(),
        maximum_estimated_conductivity_error_w_m_k: conductivity_budget
            .maximum_estimated_error_w_m_k(),
        initial_caloric_defect_j,
        endpoint_caloric_defect_j,
        maximum_query_caloric_defect_j,
        native_energy_defect_j,
        absolute_enthalpy_defect_j,
    })
}

impl MaterialAdvectionStep {
    /// Conduct Fourier heat through the original gas mixture caloric laws.
    ///
    /// The existing native tensor-caloric owner advances m*h with immutable
    /// carrier mass, exact original T(h,Y), and original dT/dh=1/Cp. The opening
    /// prepared gas image supplies the original opening conductivity. The
    /// selected model either holds it or integrates the same original WMS
    /// k(T,Y) through every crossed pure-source interval. Cp is never replaced
    /// by a constant, and zero molecular species D does
    /// not disable heat conduction. All returned heat is already included in h.
    ///
    /// Composition and source context remain fixed. HeldAtOpening retains
    /// first-order conductivity-freezing error, including with SDIRK2.
    /// SourceIntegral prepares the original species mu/k source lines once
    /// and integrates their WMS mixture within explicit quadrature tolerances.
    /// Its original endpoint k supplies the native Kirchhoff tangent. Smooth
    /// SDIRK2 order requires resolving the remaining source and heat tolerances;
    /// unresolved source corners require temporal refinement. The integral
    /// error is an estimate, not a certified enclosure.
    /// The same authored ThermalSetup is used at both SDIRK2 stage times.
    /// Its original flux/temperature/film/source/contact/sink admissions apply.
    /// Fixed cell temperatures derive their unique h from the SAME gas law.
    ///
    /// Every voxel must be fluid. Solid storage requires its own positive
    /// reference mass and caloric law through the native heterogeneous API.
    /// The pressure-held enthalpy convention matches the species owner; EOS
    /// density remains diagnostic and never replaces the accepted carrier mass.
    /// This is a heat component, not a sealed rigid-volume u/Cv pressure update.
    ///
    /// Each cell's immutable caloric model is constructed once, outside native
    /// Newton queries. The native heat owner retains its spatial allocations;
    /// this adapter adds fallibly reserved source/law/tensor arrays and one
    /// source model per cell. Original source tables are not rebuilt per query.
    /// Any source, solver, work, energy or cancellation refusal returns no
    /// partial state and leaves the carrier and supplied inventories unchanged.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn conduct_thermochemical_species<'a, const N: usize>(
        &self,
        domain: &VoxelDomain,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        setup: &ThermalSetup,
        start_time_s: f64,
        duration_s: f64,
        source: &GasSpeciesEnthalpySource<'a, N>,
        controls: ThermochemicalConductionControls,
        gate: &CancelGate,
    ) -> Result<ThermochemicalSpeciesConductionStep<'a, N>, ChtError> {
        poll(gate)?;
        admit_interval(start_time_s, duration_s, &controls)?;
        let cells = domain.cell_count();
        for cell in 0..cells {
            poll_every(gate, cell)?;
            if !domain.is_fluid(cell) {
                return Err(invalid(
                    "gas conduction requires a fluid caloric owner in every voxel",
                ));
            }
        }
        let budget = GasConductionQueryBudget::new(controls.maximum_caloric_queries, gate)?;
        let conductivity_budget = GasConductivityBudget::new(controls.conductivity, gate)?;
        // The original preparation admits shapes, canonical species partition,
        // carrier topology, caloric branch and complete source transport.
        let initial =
            source.prepare(domain, self, species_mass_kg, specific_enthalpy_j_kg, gate)?;
        let initial_defect = admitted_defect(
            initial.absolute_caloric_residual_j(),
            controls.maximum_enthalpy_defect_j,
        )?;
        let mut calorics = reserved(cells)?;
        let mut tensors = reserved(cells)?;
        for cell in 0..cells {
            poll_every(gate, cell)?;
            calorics.push(GasConductionCaloric::new(
                source,
                species_mass_kg[cell],
                self.mass_kg[cell],
                &budget,
            )?);
            let k = initial
                .transport(cell)
                .ok_or_else(|| invalid("gas conduction lost its original fluid transport"))?
                .thermal_conductivity_w_m_k();
            tensors.push(GasConductionConductivity::new(
                &calorics[cell],
                [
                    source.controls.minimum_temperature_k,
                    source.controls.maximum_temperature_k,
                ],
                k,
                &conductivity_budget,
            )?);
        }
        drop(initial);
        let mut prescribed = reserved(setup.fixed_temperature.len())?;
        for (index, &(cell, temperature)) in setup.fixed_temperature.iter().enumerate() {
            poll_every(gate, index)?;
            let law = calorics
                .get(cell)
                .ok_or_else(|| invalid("prescribed gas temperature cell is outside its domain"))?;
            prescribed.push((cell, law.specific_enthalpy_at_temperature(temperature)?));
        }
        let mut thermal_laws: Vec<&dyn TensorCaloricLaw> = reserved(cells)?;
        let mut conductivity_laws: Vec<&dyn TensorConductivityLaw> = reserved(cells)?;
        for cell in 0..cells {
            poll_every(gate, cell)?;
            thermal_laws.push(&calorics[cell]);
            conductivity_laws.push(&tensors[cell]);
        }
        let storage = TensorCaloricStorage {
            laws: &thermal_laws,
            reference_mass_kg: &self.mass_kg,
            previous_specific_energy_j_kg: specific_enthalpy_j_kg,
        };
        let thermal = advance_native(
            domain,
            &[],
            setup,
            &conductivity_laws,
            storage,
            &prescribed,
            start_time_s,
            duration_s,
            &controls,
            None,
            gate,
        )?;
        let maximum_query_caloric_defect_j =
            gas_query_defect(&calorics, controls.maximum_enthalpy_defect_j, gate)?;
        let caloric_queries = budget.queries();
        // Release the per-cell source models and all borrowed law views before
        // preparing the reusable endpoint image. Native heat owns its own output.
        drop(thermal_laws);
        drop(conductivity_laws);
        drop(calorics);
        drop(tensors);
        drop(prescribed);
        let state = source.prepare(
            domain,
            self,
            species_mass_kg,
            &thermal.endpoint().specific_energy_j_kg,
            gate,
        )?;
        let report = conduction_report(
            &thermal,
            initial_defect,
            state.absolute_caloric_residual_j(),
            maximum_query_caloric_defect_j,
            caloric_queries,
            &conductivity_budget,
            &controls,
        )?;
        poll(gate)?;
        Ok(ThermochemicalSpeciesConductionStep {
            thermal,
            state,
            report,
        })
    }
}
