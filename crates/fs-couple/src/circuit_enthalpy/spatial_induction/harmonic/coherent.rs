//! Exact coherent AC heat on the original finite-cell enthalpy residual.
//!
//! The field owner integrates the actual thermal interval once, retaining
//! carrier and beat terms and the original spatial current distribution. Its
//! cell heat is then reused through the original nonlinear thermal solve.
//! Every candidate still resolves the complete material heat assembly.
//!
//! This receiver requires stationary geometry and electrical/magnetic source
//! laws independent of the evolving temperature. It neither averages carriers
//! nor adds the original cycle means. Only the caloric energy account closes
//! here: supplied electrical work, magnetic storage, switch-on transients,
//! magnetic entropy and electromagnetic feedback remain with their owners.

use super::{SpatialHarmonicThermalSample, field::HarmonicFieldMaterials};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, CircuitEnthalpySolve, DrivenCoefficients, DriverAccount, ThermalDriver,
    admit_driver_unknowns, cell_driver_candidate,
    contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance},
    finite,
    magnetic::{self, MagneticThermalError, MagneticThermalNetwork, check, invalid, material},
    material::MaterialNetworkError,
    spatial_induction::{
        SpatialMagneticMaterialSample, SpatialMagneticRegion, SpatialMagneticThermalNetwork,
    },
};
use fs_em::spatial::vector_potential::induction::harmonic::coherent_exposure::{
    SpatialHarmonicInductionExposure, SpatialHarmonicInductionExposureError,
    SpatialHarmonicInductionExposureResult,
};
use fs_exec::Cx;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::conductor::spatial::SpatialConductorMaterial;
use fs_math::dd::Dd;
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{CellId, DepositionSimulation, EnergyAudit};
use std::{cell::RefCell, collections::BTreeMap};

/// Preserve the original material, observation and heat-provider refusals.
#[derive(Debug)]
pub enum CoherentInductionHeatError<E> {
    Material(MagneticThermalError),
    Observation(SpatialHarmonicInductionExposureError),
    Heat(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for CoherentInductionHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coherent 3D induction heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for CoherentInductionHeatError<E> {}
type Error<E> = CircuitEnthalpyError<CoherentInductionHeatError<E>>;

fn native<E>(error: magnetic::Error) -> Error<E> {
    error.map_source(CoherentInductionHeatError::Material)
}
fn heat<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(CoherentInductionHeatError::Heat(other.into())),
    }
}
fn observation<E>(error: SpatialHarmonicInductionExposureError) -> Error<E> {
    match error {
        SpatialHarmonicInductionExposureError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(CoherentInductionHeatError::Observation(other)),
    }
}

/// Original physical observation and finite caloric balance at acceptance.
///
/// `caloric_energy_residual_j` is thermal storage change minus deposited Joule
/// heat and other external heat. It makes no electromagnetic supply/storage
/// claim. The borrowed source and original cell records remain inspectable.
pub struct CoherentInductionHeatReport<'s, 'm, C: SpatialConductorMaterial, S> {
    pub exposure: SpatialHarmonicInductionExposureResult<'s, 'm, C>,
    /// Native audit; external heat includes the deposited induction heat.
    pub thermal: EnergyAudit,
    pub source: SpatialHarmonicThermalSample<S>,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    /// Complete thermal/material evaluations in Newton and final acceptance.
    pub thermal_trials: usize,
    pub deposited_heat_j: f64,
    /// Excludes deposited induction heat and numerical closed-radiation net.
    pub other_external_heat_j: f64,
    pub caloric_energy_residual_j: f64,
}

/// Unpublished thermal endpoint. The prepared electromagnetic bank is unchanged.
pub struct CoherentInductionHeatCandidate<'t, 's, 'm, C: SpatialConductorMaterial, S> {
    thermal: DepositionSimulation<'t>,
    report: CoherentInductionHeatReport<'s, 'm, C, S>,
}
impl<'t, 's, 'm, C: SpatialConductorMaterial, S> CoherentInductionHeatCandidate<'t, 's, 'm, C, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn report(&self) -> &CoherentInductionHeatReport<'s, 'm, C, S> {
        &self.report
    }
    /// All reported induction heat is already included in this endpoint.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        CoherentInductionHeatReport<'s, 'm, C, S>,
    ) {
        (self.thermal, self.report)
    }
}

// Reuse the original deposition/Newton owner with a pure heat-transfer account.
// This is energy entering the CALORIC receiver, not an inferred electrical
// source work. The public report deliberately exposes only that caloric scope.
struct HeatTransfer(f64);
impl ThermalDriver for HeatTransfer {
    type Model = ();
    type Step = ();
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.0,
            numerical_dissipation_j: 0.0,
            stored_change_j: 0.0,
            supplied_j: self.0,
        }
    }
    fn into_parts(self) -> ((), ()) {
        ((), ())
    }
}

/// Deposit the original coherent field's exact interval heat into finite cells.
///
/// The interval is `[thermal.time_s(), target_s]` on the AC bank's authored
/// phase clock. All conducting tetrahedra require whole-volume solid caloric
/// receivers. Scalar and oriented conductivity, permeability, source selection
/// and nonthermal coordinates remain unchanged; temperature-dependent electrical
/// or magnetic laws refuse even at coincident endpoint values. Thermal laws,
/// contacts, layer conduction, radiation and finite coolant still resolve at
/// each actual proposed endpoint through the supplied `HeatAssembly`.
///
/// The electromagnetic observation runs once, outside Newton. Thermal rates
/// retain the original backward-Euler discretization and require time refinement.
/// No magnetic initial condition, mechanical work, source work or numerical loss
/// is deposited as heat. Process start events first and split at source/topology
/// discontinuities. A refusal changes neither accepted owner.
#[allow(clippy::too_many_arguments)]
pub fn heat_assembly_candidate<
    't,
    's,
    'm,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    exposure: &SpatialHarmonicInductionExposure<'s, 'm, C>,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<CoherentInductionHeatCandidate<'t, 's, 'm, C, A::Sample>, Error<A::Error>> {
    let cancel = RefCell::new(|| cx.checkpoint().is_err());
    check(&cancel).map_err(native)?;
    admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
    let dt = finite(target_s - thermal.time_s())?;
    let window = HarmonicWindow::new(thermal.time_s(), dt).map_err(|_| {
        native(invalid(
            "coherent induction heat needs a positive representable thermal interval",
        ))
    })?;
    for vertices in exposure.prescribed_velocity_m_s().chunks(256) {
        check(&cancel).map_err(native)?;
        if vertices.iter().flatten().any(|v| *v != 0.0) {
            return Err(native(invalid(
                "fixed coherent caloric receivers require stationary conductors",
            )));
        }
    }
    let prepared = provider.prepared_heat_network();
    let bulk = prepared.thermal_material_network().map_err(heat)?;
    prepared
        .admit_field_independent_boundaries(cx, &[])
        .map_err(heat)?;
    if bulk.cells.len() != thermal.cells().len() {
        return Err(native(invalid(
            "complete coherent induction thermal bindings required",
        )));
    }
    let layout = SpatialMagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: regions,
        volume_tolerance_m3,
    };
    let regions = layout.regions();
    let network = MagneticThermalNetwork {
        cells: bulk.cells,
        links: bulk.links,
        magnetic_regions: &regions,
        volume_tolerance_m3,
    };
    let problem = exposure.source().problem();
    let field_materials =
        HarmonicFieldMaterials::new(problem, exposure.source().conductors(), &cancel)
            .map_err(native)?;
    let routes = field_materials
        .routes(&network, thermal, &cancel)
        .map_err(native)?;
    let bindings = bulk
        .cells
        .iter()
        .map(|binding| (binding.cell, binding))
        .collect::<BTreeMap<_, _>>();

    let resolve = |state: &DepositionSimulation<'_>| {
        check(&cancel).map_err(native)?;
        let coefficients = provider
            .resolve_heat(cx, state)
            .map_err(|error| error.map_source(CoherentInductionHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(native(invalid(
                "coherent heating cannot discard graph electrical work",
            )));
        }
        let (refreshed, magnetic) = field_materials
            .resolve(
                &network,
                state,
                &coefficients.bulk.source.as_ref().material,
                &routes,
                &cancel,
            )
            .map_err(native)?;
        for (cell, source) in refreshed {
            check(&cancel).map_err(native)?;
            let original = exposure
                .source()
                .conductors()
                .get(cell)
                .and_then(Option::as_ref)
                .ok_or_else(|| native(invalid("coherent heat lost its original conductor")))?;
            if !original.has_same_electrical_source(&source)
                || original.conductivity() != source.conductivity()
            {
                return Err(native(invalid(
                    "coherent heat cannot change its prepared conductivity",
                )));
            }
            let receiver = routes
                .conductors
                .get(&cell)
                .ok_or_else(|| native(invalid("coherent heat lost its caloric route")))?;
            let binding = bindings
                .get(receiver)
                .ok_or_else(|| native(invalid("coherent heat lost its caloric binding")))?;
            for key in source.source_keys() {
                let property = source.material().property_by_key(key).ok_or_else(|| {
                    native(invalid("coherent heat lost its original electrical source"))
                })?;
                if property
                    .source_axes()
                    .iter()
                    .any(|(axis, _)| axis == binding.point.temperature_axis)
                {
                    return Err(native(invalid(
                        "coherent electrical sources must be independent of evolving temperature",
                    )));
                }
                binding
                    .point
                    .admit_static(property)
                    .map_err(|e| native(material(MaterialNetworkError::Material(e))))?;
            }
        }
        Ok((coefficients.bulk, magnetic))
    };
    // Pin the accepted correspondence before observing or constructing a trial.
    resolve(thermal)?;
    let observed = exposure
        .observe(window, &mut || (*cancel.borrow_mut())())
        .map_err(observation)?;
    if observed.cells().len() != routes.conductors.len() {
        return Err(native(invalid(
            "coherent heat omitted an original conductor",
        )));
    }
    let mut sums = BTreeMap::<CellId, Dd>::new();
    for record in observed.cells() {
        check(&cancel).map_err(native)?;
        let cell = *routes
            .conductors
            .get(&record.cell)
            .ok_or_else(|| native(invalid("coherent heat lost its original caloric receiver")))?;
        let sum = sums.entry(cell).or_insert(Dd::ZERO);
        *sum = *sum + Dd::from_f64(record.mean_ohmic_power_w);
    }
    let powers = sums
        .into_iter()
        .map(|(id, value)| Ok((id, finite(value.to_f64())?)))
        .collect::<Result<BTreeMap<_, _>, Error<A::Error>>>()?;
    let (next, (), report) = cell_driver_candidate(
        thermal,
        target_s,
        controls,
        initial_energy_increments_j,
        |state, _remaining| {
            let (coefficients, magnetic) = resolve(state)?;
            Ok(DrivenCoefficients {
                circuit: HeatTransfer(observed.ohmic_heat_j),
                sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k,
                circuit_trials: 1,
                source: SpatialHarmonicThermalSample {
                    network: coefficients.source,
                    magnetic: magnetic
                        .into_iter()
                        .map(|s| SpatialMagneticMaterialSample {
                            tetrahedron: s.triangle,
                            cell: s.cell,
                            material: s.material,
                        })
                        .collect(),
                },
            })
        },
        |_: &HeatTransfer, _, checkpoint| {
            checkpoint()?;
            Ok(powers.clone())
        },
        || (*cancel.borrow_mut())(),
    )?;
    let sample = report.source.network.as_ref();
    closed_balance(
        dt,
        report.total_energy_residual_j,
        sample,
        controls.total_energy_tolerance_j,
    )
    .map_err(heat)?;
    let radiation_net = finite(dt * sample.radiation_net_power_w())?;
    let other_external_heat_j = finite(report.non_electrical_heat_j - radiation_net)?;
    let caloric_energy_residual_j = finite(report.total_energy_residual_j + radiation_net)?;
    check(&cancel).map_err(native)?;
    Ok(CoherentInductionHeatCandidate {
        thermal: next,
        report: CoherentInductionHeatReport {
            exposure: observed,
            thermal: report.thermal,
            source: report.source,
            nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j,
            thermal_trials: report.circuit_trials,
            deposited_heat_j: report.deposited_joule_j,
            other_external_heat_j,
            caloric_energy_residual_j,
        },
    })
}
