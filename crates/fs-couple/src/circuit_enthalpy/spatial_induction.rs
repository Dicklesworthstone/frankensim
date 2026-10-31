//! Three-dimensional induction, source temperature and finite caloric energy.
//!
//! This composes the original coupled Whitney A–phi field with the original
//! finite-cell enthalpy/Newton owner. Every endpoint temperature proposal starts
//! from identical accepted magnetic and thermal histories. The original material
//! source supplies sigma(T), affine mu(T), field entropy and internal energy.
//! Actual local Ohmic heat and signed reversible heat have separate accounts;
//! no fitted core loss is added. Geometry and conducting occupancy stay fixed.

/// Source-pinned saturation/recoil and local conductivity-temperature feedback.
pub mod nonlinear;
// Complete sourced heat transport uses these same retained field residuals.
mod heat;
pub use heat::{
    SpatialInductionHeatError, heat_assembly_candidate, nonlinear_heat_assembly_candidate,
    nonlinear_thermal_heat_assembly_candidate,
};

use super::{
    CircuitEnthalpyError, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients,
    DrivenCoefficients, DriverAccount, ThermalDriver, cell_driver_candidate, finite,
    magnetic::{
        Error, MagneticRegion, MagneticRoutes, MagneticThermalError, MagneticThermalNetwork,
        caloric::{MagneticEnergyCell, MagneticThermalLaw, PreparedMagneticCaloric},
        check, invalid, material, same_material,
    },
    material::{ConductionLink, MaterialCell, MaterialNetwork, MaterialNetworkSample},
};
use fs_em::spatial::{
    SpatialMagneticError, SpatialMagneticMaterial,
    vector_potential::induction::{
        SpatialInduction, SpatialInductionCandidate, SpatialInductionCell,
        SpatialInductionControls, SpatialInductionStep,
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    field_storage::{FieldStorageError, admit_magnetic_temperature_refresh},
    state_point::ResolvedMaterialStatePoint,
};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

pub use fs_material::field_storage::thermal::{
    FieldThermalState, magnetic::ThermalMagneticMaterial,
};

/// Exactly one correspondence per original tetrahedron. All conductors require
/// a real receiving cell; fixed external regions are thermally unowned.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpatialMagneticRegion {
    Cell { tetrahedron: usize, cell: CellId },
    Fixed { tetrahedron: usize },
}

/// Original temperature-dependent source and admitted B envelope of one tet.
#[derive(Clone, Debug, PartialEq)]
pub struct SpatialMagneticThermalLaw<L = ThermalMagneticMaterial> {
    pub tetrahedron: usize,
    pub law: L,
    pub maximum_flux_density_t: f64,
}

/// Data-only original caloric/card/source and fixed-volume correspondence.
pub struct SpatialMagneticThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub magnetic_regions: &'a [SpatialMagneticRegion],
    pub volume_tolerance_m3: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpatialMagneticMaterialSample {
    pub tetrahedron: usize,
    pub cell: Option<CellId>,
    pub material: ResolvedMaterialStatePoint,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpatialMagneticCaloricSample {
    pub tetrahedron: usize,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    pub flux_density_t: [f64; 2],
    pub state: [FieldThermalState; 2],
    pub parameter_work_j: f64,
    pub reversible_heat_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SpatialInductionThermalSource<S = MaterialNetworkSample> {
    pub network: S,
    pub magnetic: Vec<SpatialMagneticMaterialSample>,
    pub fields: Vec<SpatialMagneticCaloricSample>,
    /// Includes fixed exterior regions, for which U=psi.
    pub field_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}

pub struct SpatialInductionThermalReport<D = SpatialInductionStep, S = MaterialNetworkSample> {
    pub coupled: CircuitEnthalpyReport<SpatialInductionThermalSource<S>, D>,
    /// Authored thermal sources, conduction/reservoir transfer and boundary heat.
    /// Internal reversible field heat is excluded from this external account.
    pub external_heat_j: f64,
    /// Delta U_field + Delta H_caloric + W_mechanical + numerical loss
    /// - W_impressed - Q_external. Numerical loss is zero for the midpoint law.
    pub total_energy_residual_j: f64,
}

pub struct SpatialInductionThermalCandidate<'t, 'm, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    field: SpatialInduction<'m>,
    report: SpatialInductionThermalReport<SpatialInductionStep, S>,
}
impl<'t, 'm, S> SpatialInductionThermalCandidate<'t, 'm, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &SpatialInduction<'m> {
        &self.field
    }
    pub fn report(&self) -> &SpatialInductionThermalReport<SpatialInductionStep, S> {
        &self.report
    }
    /// Publish BOTH histories together; all returned Ohmic/reversible heat has
    /// already entered the thermal state and must not be deposited again.
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        SpatialInduction<'m>,
        SpatialInductionThermalReport<SpatialInductionStep, S>,
    ) {
        (self.thermal, self.field, self.report)
    }
}

// The shared thermal residual sees only this actual field candidate. Parameter
// work is separately returned to the same caloric cells by the entropy owner.
impl<'m> ThermalDriver for SpatialInductionCandidate<'m> {
    type Model = SpatialInduction<'m>;
    type Step = SpatialInductionStep;
    fn account(&self) -> DriverAccount<'_> {
        let s = self.step();
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: s.heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: s.stored_after_j - s.stored_before_j - s.parameter_work_j,
            supplied_j: s.impressed_work_j - s.mechanical_work_j,
        }
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        SpatialInductionCandidate::into_parts(self)
    }
}

// Preserve original magnetic/material refusals while retaining each complete
// heat assembly's own typed source error in the shared residual.
fn native_error<E: From<MagneticThermalError>>(error: Error) -> CircuitEnthalpyError<E> {
    error.map_source(E::from)
}

fn induction(error: SpatialMagneticError) -> Error {
    match error {
        SpatialMagneticError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(
            super::magnetic::MagneticThermalError::SpatialInduction(other),
        ),
    }
}

impl SpatialMagneticThermalNetwork<'_> {
    fn regions(&self) -> Vec<MagneticRegion> {
        self.magnetic_regions
            .iter()
            .map(|r| match *r {
                SpatialMagneticRegion::Cell { tetrahedron, cell } => MagneticRegion::Cell {
                    triangle: tetrahedron,
                    cell,
                },
                SpatialMagneticRegion::Fixed { tetrahedron } => MagneticRegion::Fixed {
                    triangle: tetrahedron,
                },
            })
            .collect()
    }

    /// Joint source-resolved sigma(T), mu(T), A–phi and finite enthalpy step.
    ///
    /// All original conducting cells need whole-volume solid caloric receivers.
    /// Every thermally mapped magnetic cell needs a source-pinned derivative law
    /// and a B envelope whose combined caloric/field capacity stays positive.
    /// Fixed regions require temperature-independent permeability. The field
    /// entropy owner deposits Q_rev = -P - delta(T*S_field), separately from
    /// the exact native sigma*E_conductor^2 integrals. The physical energy audit
    /// uses U_field=psi+T*S_field, not a duplicated free-energy store.
    ///
    /// Prescribed velocity may transfer mechanical work to/from an external
    /// receiver; it does not move the mesh or integrate mechanics. Sources are
    /// linear magnetic laws with scalar conductivity. Endpoint coefficients and
    /// backward-Euler enthalpy give first-order coupled trajectories in general.
    /// Recoil, saturation, hysteresis, thermal expansion and moving interfaces
    /// require their own thermodynamic/formulation extensions.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInductionThermalCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        let material_network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = material_network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_resolved_candidate(
            thermal,
            target_s,
            field,
            endpoint_current_a_m2,
            velocity_m_s,
            bindings,
            controls,
            field_controls,
            assembly,
            initial_energy_increments_j,
            |trial| {
                prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)
            },
            |sample| sample,
            &cancel,
        )
    }

    // Bulk-only and complete sourced heat transport enter the SAME original
    // material/field/caloric residual. The source payload remains inspectable.
    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, 'm, S, E: From<MagneticThermalError>>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &SpatialInduction<'m>,
        endpoint_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        bindings: &[SpatialMagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: SpatialInductionControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        resolve_sources: impl Fn(
            &DepositionSimulation<'_>,
        ) -> Result<CircuitThermalCoefficients<S>, CircuitEnthalpyError<E>>,
        sample: impl Fn(&S) -> &MaterialNetworkSample,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<SpatialInductionThermalCandidate<'t, 'm, S>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(native_error)?;
        super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(native_error(invalid(
                "matching spatial field/thermal clocks and complete cells required",
            )));
        }
        let regions = self.regions();
        let network = MagneticThermalNetwork {
            cells: self.cells,
            links: self.links,
            magnetic_regions: &regions,
            volume_tolerance_m3: self.volume_tolerance_m3,
        };
        let problem = field.problem();
        let originals = problem
            .materials()
            .iter()
            .map(|m| match m {
                SpatialMagneticMaterial::Linear(m) => Ok(m),
                SpatialMagneticMaterial::Recoil(_) => Err(invalid(
                    "spatial caloric coupling requires induced-linear magnetic materials",
                )),
            })
            .collect::<Result<Vec<_>, Error>>()
            .map_err(native_error)?;
        let volumes = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(&cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("original tetrahedron has no physical volume"))
            })
            .collect::<Result<Vec<_>, Error>>()
            .map_err(native_error)?;
        let conductors = field
            .conductors()
            .iter()
            .enumerate()
            .filter_map(|(cell, source)| source.as_ref().map(|s| (cell, s)))
            .collect::<BTreeMap<_, _>>();
        let conducting_cells = conductors.keys().copied().collect::<Vec<_>>();
        let routes = network
            .routes_for_volumes(thermal, &volumes, &conducting_cells, &cancel)
            .map_err(native_error)?;
        let laws = bindings
            .iter()
            .map(|b| MagneticThermalLaw {
                triangle: b.tetrahedron,
                law: b.law.clone(),
                maximum_flux_density_t: b.maximum_flux_density_t,
            })
            .collect::<Vec<_>>();
        let b = field
            .flux_density_t(&mut || (*cancel.borrow_mut())())
            .map_err(induction)
            .map_err(native_error)?;
        let caloric = PreparedMagneticCaloric::new_for_volumes(
            &network, thermal, &volumes, &originals, &b, &laws, &routes, &cancel,
        )
        .map_err(native_error)?;
        let admit = |_: usize,
                     before: &ResolvedMaterialStatePoint,
                     after: &ResolvedMaterialStatePoint,
                     _: &MaterialCell<'_>| {
            admit_magnetic_temperature_refresh(before, after).map_err(|e| match e {
                FieldStorageError::Invalid(reason) => invalid(reason),
                other => CircuitEnthalpyError::Source(
                    super::magnetic::MagneticThermalError::Magnetic(other),
                ),
            })
        };
        let initial = resolve_sources(thermal)?;
        let (accepted, _) = network.materials_at_using(
            thermal,
            &originals,
            &conductors,
            sample(&initial.source),
            &routes,
            &cancel,
            admit,
        )
        .map_err(native_error)?;
        if accepted.len() != conductors.len()
            || accepted.iter().any(|(cell, source)| {
                conductors.get(cell).is_none_or(|old| {
                    old.source_basis() != source.source_basis()
                        || old.source_key() != source.source_key()
                        || !same_material(old.material(), source.material())
                })
            })
        {
            return Err(native_error(invalid(
                "accepted spatial conductor source differs from its actual thermal endpoint",
            )));
        }
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _| {
                check(&cancel).map_err(native_error)?;
                let mut coefficients = resolve_sources(trial)?;
                let (sources, magnetic) = network.materials_at_using(
                    trial,
                    &originals,
                    &conductors,
                    sample(&coefficients.source),
                    &routes,
                    &cancel,
                    admit,
                )
                .map_err(native_error)?;
                caloric
                    .validate_sources(trial, &magnetic, &cancel)
                    .map_err(native_error)?;
                let mut refreshed = vec![None; originals.len()];
                for (cell, source) in sources {
                    refreshed[cell] = Some(source);
                }
                let proposal = field
                    .step_with_materials_candidate(
                        endpoint_current_a_m2,
                        velocity_m_s,
                        dt,
                        magnetic
                            .iter()
                            .map(|s| SpatialMagneticMaterial::Linear(s.material.clone()))
                            .collect(),
                        refreshed,
                        assembly,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(induction)
                    .map_err(native_error)?;
                let exchange = caloric.deposit(
                    trial,
                    proposal
                        .step()
                        .magnetic_cells
                        .iter()
                        .map(|cell| MagneticEnergyCell {
                            triangle: cell.cell,
                            flux_density_t: cell.flux_density_t,
                            free_energy_j: cell.free_energy_j,
                            parameter_work_j: cell.parameter_work_j,
                        }),
                    &mut coefficients.sources_w,
                    dt,
                    controls.total_energy_tolerance_j,
                    &cancel,
                )
                .map_err(native_error)?;
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    circuit_trials: 1,
                    source: SpatialInductionThermalSource {
                        network: coefficients.source,
                        magnetic: magnetic
                            .into_iter()
                            .map(|s| SpatialMagneticMaterialSample {
                                tetrahedron: s.triangle,
                                cell: s.cell,
                                material: s.material,
                            })
                            .collect(),
                        fields: exchange
                            .fields
                            .into_iter()
                            .map(|s| SpatialMagneticCaloricSample {
                                tetrahedron: s.triangle,
                                cell: s.cell,
                                temperature_k: s.temperature_k,
                                flux_density_t: s.flux_density_t,
                                state: s.state,
                                parameter_work_j: s.parameter_work_j,
                                reversible_heat_j: s.reversible_heat_j,
                            })
                            .collect(),
                        field_internal_energy_j: exchange.field_internal_energy_j,
                        reversible_heat_j: exchange.reversible_heat_j,
                    },
                })
            },
            |proposal: &SpatialInductionCandidate<'_>, dt, checkpoint| {
                deposit_cell_heat(&proposal.step().cells, &routes, dt, checkpoint)
                    .map_err(native_error)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let external_heat_j =
            finite(report.non_electrical_heat_j - report.source.reversible_heat_j)?;
        let total_energy_residual_j = finite(
            report.source.field_internal_energy_j[1] - report.source.field_internal_energy_j[0]
                + report.thermal.stored_energy_change_j
                + report.electrical.mechanical_work_j
                - report.electrical.impressed_work_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(native_error(invalid(
                "spatial field internal energy, finite enthalpy and external work do not close",
            )));
        }
        check(&cancel).map_err(native_error)?;
        Ok(SpatialInductionThermalCandidate {
            thermal,
            field,
            report: SpatialInductionThermalReport {
                coupled: report,
                external_heat_j,
                total_energy_residual_j,
            },
        })
    }
}

// Linear and nonlinear induction recover the SAME original-cell current and
// heat type. Share its exact receiving-cell account rather than reintegrating J.
fn deposit_cell_heat(
    losses: &[SpatialInductionCell],
    routes: &MagneticRoutes,
    dt: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<BTreeMap<CellId, f64>, Error> {
    if losses.len() != routes.conductors.len() {
        return Err(invalid("spatial heat omitted original conducting cells"));
    }
    let mut powers = BTreeMap::new();
    for loss in losses {
        checkpoint()?;
        let id = *routes
            .conductors
            .get(&loss.cell)
            .ok_or_else(|| invalid("spatial eddy heat lacks its original thermal receiver"))?;
        let power = finite(loss.heat_j / dt)?;
        if loss.heat_j < 0.0 || loss.heat_j != 0.0 && power == 0.0 {
            return Err(invalid(
                "spatial eddy heat is negative or lost in power conversion",
            ));
        }
        let previous = *powers.get(&id).unwrap_or(&0.0);
        let combined = finite(previous + power)?;
        if power != 0.0 && combined == previous {
            return Err(invalid(
                "spatial eddy heat disappeared in its thermal receiver",
            ));
        }
        powers.insert(id, combined);
    }
    Ok(powers)
}
