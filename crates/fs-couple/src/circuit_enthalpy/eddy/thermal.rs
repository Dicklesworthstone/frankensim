//! Source-resolved mu(T), induction and caloric enthalpy on one accepted clock.
//! All trials restart from original A and H. The native magnetic endpoint-effort
//! step exposes parameter work P; this owner deposits Q_rev=−P−delta(T*S_field)
//! beside (and separately from) its actual spatial Ohmic heat. Consequently the
//! independent audit uses U_field=psi+T*S and original caloric energy exactly once.

use super::super::{
    CircuitEnthalpyReport, CircuitEnthalpySolve, DrivenCoefficients, cell_driver_candidate, finite,
    magnetic::{
        Error, MagneticRoutes, MagneticThermalError, MagneticThermalNetwork, check, invalid,
        material, same_material,
    },
    material::{MaterialNetwork, MaterialNetworkError},
};
use super::{EddyThermalSample, deposit_heat, induction};
use fs_em::planar::eddy::{EddyConductor, EddyControls, EddyStep, PlanarEddy, PlanarEddyCandidate};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::field_storage::{
    FieldStorageError, MAGNETIC_PERMEABILITY_PROPERTY,
    thermal::{FieldThermalState, magnetic::ThermalMagneticMaterial},
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// One original magnetic triangle's complete source and admitted field envelope.
/// Its receiving cell comes from the shared MagneticThermalNetwork region map.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticThermalLaw {
    pub triangle: usize,
    pub law: ThermalMagneticMaterial,
    pub maximum_flux_density_t: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticCaloricSample {
    pub triangle: usize,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    pub flux_density_t: [f64; 2],
    pub state: [FieldThermalState; 2],
    pub parameter_work_j: f64,
    pub reversible_heat_j: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct EddyCaloricSource {
    pub material: EddyThermalSample,
    pub fields: Vec<MagneticCaloricSample>,
    /// Includes temperature-independent fixed magnetic regions as psi=U.
    pub field_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}
pub struct EddyCaloricReport {
    /// Native shared thermal residual plus native spatial induction report.
    /// The driver's work account removes its explicit parameter term; the
    /// independent physical-U audit below closes its thermal destination.
    pub coupled: CircuitEnthalpyReport<EddyCaloricSource, EddyStep>,
    pub external_heat_j: f64,
    pub total_energy_residual_j: f64,
}
pub struct EddyCaloricCandidate<'t, 'm> {
    thermal: DepositionSimulation<'t>,
    field: PlanarEddy<'m>,
    report: EddyCaloricReport,
}
impl<'t, 'm> EddyCaloricCandidate<'t, 'm> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &PlanarEddy<'m> {
        &self.field
    }
    pub fn report(&self) -> &EddyCaloricReport {
        &self.report
    }
    /// Both accepted histories must publish together; heat is already deposited.
    pub fn into_parts(self) -> (DepositionSimulation<'t>, PlanarEddy<'m>, EddyCaloricReport) {
        (self.thermal, self.field, self.report)
    }
}
fn field_error(e: FieldStorageError) -> Error {
    super::super::CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e))
}
fn magnitude(b: [f64; 2]) -> f64 {
    b[0].hypot(b[1])
}
struct Initial {
    cell: CellId,
    temperature_k: f64,
    b_t: f64,
    volume_m3: f64,
    state: FieldThermalState,
}

fn admit(
    network: &MagneticThermalNetwork<'_>,
    thermal: &DepositionSimulation<'_>,
    field: &PlanarEddy<'_>,
    bindings: &[MagneticThermalLaw],
    routes: &MagneticRoutes,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<BTreeMap<usize, Initial>, Error> {
    let native = field
        .field(&mut || (*cancel.borrow_mut())())
        .map_err(induction)?;
    let expected = routes.magnetic.values().filter(|c| c.is_some()).count();
    if bindings.len() != expected {
        return Err(invalid(
            "every thermally mapped magnetic triangle needs exactly one derivative law",
        ));
    }
    let mut initial = BTreeMap::new();
    let mut capacity = BTreeMap::<CellId, f64>::new();
    let mut ranges = BTreeMap::<CellId, [f64; 2]>::new();
    for binding in bindings {
        check(cancel)?;
        let triangle = binding.triangle;
        let id = routes
            .magnetic
            .get(&triangle)
            .copied()
            .flatten()
            .ok_or_else(|| invalid("magnetic thermal law requires an original receiving cell"))?;
        if initial.contains_key(&triangle)
            || !binding.maximum_flux_density_t.is_finite()
            || binding.maximum_flux_density_t < 0.0
        {
            return Err(invalid(
                "unique original magnetic triangles and finite nonnegative B envelopes required",
            ));
        }
        let source = network
            .cells
            .iter()
            .find(|s| s.cell == id)
            .ok_or_else(|| invalid("magnetic law lacks its caloric material binding"))?;
        let cell = &thermal.cells()[&id];
        let phase = source
            .curve
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|e| material(MaterialNetworkError::Phase(e)))?;
        let query = source
            .point
            .at(phase)
            .map_err(|e| material(MaterialNetworkError::Material(e)))?;
        let at = binding
            .law
            .at_temperature(phase.temperature_k())
            .map_err(field_error)?;
        if source.point.material.content_hash() != binding.law.card().content_hash()
            || source.point.temperature_axis != binding.law.temperature_axis()
            || query != at.material().query()
            || !same_material(at.material(), &field.problem().materials()[triangle])
        {
            return Err(invalid(
                "magnetic field, derivative law and actual caloric context must use the same original source bundle",
            ));
        }
        let volume = field
            .problem()
            .mesh()
            .face_measure(triangle)
            .ok_or_else(|| invalid("missing original magnetic volume"))?;
        let b = magnitude(native.magnetic_flux_density_t()[triangle]);
        if !b.is_finite() || b > binding.maximum_flux_density_t {
            return Err(invalid(
                "accepted magnetic field exceeds its thermal capacity envelope",
            ));
        }
        let range = binding.law.temperature_range_k();
        let knots = source.curve.knots();
        if range[0] < knots[0].temperature_k || range[1] > knots[knots.len() - 1].temperature_k {
            return Err(invalid(
                "the complete magnetic temperature range requires caloric support",
            ));
        }
        let span = ranges.entry(id).or_insert(range);
        span[0] = span[0].max(range[0]);
        span[1] = span[1].min(range[1]);
        let bound = binding
            .law
            .minimum_heat_capacity_j_k(binding.maximum_flux_density_t, volume)
            .map_err(field_error)?;
        let total = capacity.entry(id).or_insert(0.0);
        *total = finite(*total + bound)?;
        initial.insert(
            triangle,
            Initial {
                cell: id,
                temperature_k: phase.temperature_k(),
                b_t: b,
                volume_m3: volume,
                state: at.state(b, volume).map_err(field_error)?,
            },
        );
    }
    for (id, range) in ranges {
        check(cancel)?;
        if range[0] >= range[1] {
            return Err(invalid(
                "a magnetic cell requires a common nonempty temperature interval",
            ));
        }
        let source = network
            .cells
            .iter()
            .find(|s| s.cell == id)
            .expect("admitted binding");
        let mut minimum = f64::INFINITY;
        for pair in source.curve.knots().windows(2) {
            check(cancel)?;
            let dt = pair[1].temperature_k - pair[0].temperature_k;
            if dt <= 0.0 {
                return Err(invalid(
                    "magnetic caloric chart requires strictly increasing temperature",
                ));
            }
            if pair[1].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                let value = finite(
                    thermal.cells()[&id].mass_kg
                        * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                        / dt,
                )?;
                minimum = minimum.min(value);
            }
        }
        if !minimum.is_finite() || finite(minimum + capacity[&id])? <= 0.0 {
            return Err(invalid(
                "caloric plus ALL magnetic field heat capacities must stay positive over the declared envelope",
            ));
        }
    }
    // A fixed external region has no caloric receiver. Only temperature-INDEPENDENT
    // permeability is admitted there, so its changing field has no hidden T*S.
    for (&triangle, &cell) in &routes.magnetic {
        check(cancel)?;
        if cell.is_some() {
            continue;
        }
        let m = &field.problem().materials()[triangle];
        let property = m
            .property(MAGNETIC_PERMEABILITY_PROPERTY)
            .ok_or_else(|| invalid("fixed magnetic region lacks permeability"))?;
        if property.source_axes().iter().any(|(name, _)| {
            name == "T"
                || m.axis_quantities().get(name)
                    == Some(&QuantitySpec::dimensional(Temperature::DIMS))
                || m.axis_quantities()
                    .get(name)
                    .copied()
                    .and_then(QuantitySpec::semantic_type)
                    .is_some_and(|q| {
                        q.kind() == QuantityKind::AbsoluteTemperature
                            && q.form() == ValueForm::Static
                    })
        }) {
            return Err(invalid(
                "temperature-dependent magnetic regions require a thermal receiver and derivative law",
            ));
        }
    }
    Ok(initial)
}

impl MagneticThermalNetwork<'_> {
    /// Evolve local mu(T), sigma(T), A and caloric enthalpy together using the
    /// existing spatial field CG and existing enthalpy Newton residual. All
    /// changing magnetic regions (including insulators) need explicit derivative
    /// laws and whole-volume caloric cells. Fixed external regions need mu with
    /// no temperature dependence. Density, geometry and source choices stay fixed.
    ///
    /// The B envelope is enforced at BOTH accepted endpoints and the complete
    /// temperature interval has positive combined caloric+field capacity. This
    /// is first-order thermal coupling; time refinement remains necessary and
    /// no unconditional entropy/passivity or physical-validation claim is made.
    #[allow(clippy::too_many_arguments)]
    pub fn step_eddy_thermal_candidate<'t, 'm>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarEddy<'m>,
        impressed_jz_a_m2: &[f64],
        bindings: &[MagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: EddyControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<EddyCaloricCandidate<'t, 'm>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
        {
            return Err(invalid(
                "matching field/thermal clocks and complete original cells required",
            ));
        }
        let network = MaterialNetwork {
            cells: self.cells,
            resistors: &[],
            links: self.links,
        };
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let triangles: Vec<_> = field.conductors().iter().map(|s| s.triangle).collect();
        let routes = self.routes(thermal, field.problem(), &triangles, &cancel)?;
        let initial = admit(self, thermal, field, bindings, &routes, &cancel)?;
        let laws: BTreeMap<_, _> = bindings.iter().map(|b| (b.triangle, b)).collect();
        let conductors: BTreeMap<_, _> = field
            .conductors()
            .iter()
            .map(|s| (s.triangle, &s.conductor))
            .collect();
        let regions: BTreeMap<_, _> = field
            .conductors()
            .iter()
            .map(|s| (s.triangle, s.region))
            .collect();
        let initial_sample = prepared
            .resolve(thermal, || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (accepted_sources, _) = self.materials_at_thermal(
            thermal,
            field.problem(),
            &conductors,
            &initial_sample.source,
            &routes,
            &cancel,
        )?;
        if accepted_sources.len() != field.conductors().len()
            || accepted_sources
                .iter()
                .zip(field.conductors())
                .any(|((triangle, c), s)| {
                    *triangle != s.triangle
                        || c.source_basis() != s.conductor.source_basis()
                        || !same_material(c.material(), s.conductor.material())
                })
        {
            return Err(invalid(
                "accepted conductor sources must match actual thermal state",
            ));
        }
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, _| {
                check(&cancel)?;
                let mut coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (sources, magnetic) = self.materials_at_thermal(
                    trial,
                    field.problem(),
                    &conductors,
                    &coefficients.source,
                    &routes,
                    &cancel,
                )?;
                for s in &magnetic {
                    check(&cancel)?;
                    if let Some(cell) = s.cell {
                        let at = laws[&s.triangle]
                            .law
                            .at_temperature(trial.temperature_k(cell)?)
                            .map_err(field_error)?;
                        if !same_material(&s.material, at.material()) {
                            return Err(invalid(
                                "native magnetic refresh and derivative law selected different source bundles",
                            ));
                        }
                    }
                }
                let proposal = field
                    .step_with_materials_candidate(
                        target_s,
                        magnetic.iter().map(|s| s.material.clone()).collect(),
                        sources
                            .into_iter()
                            .map(|(triangle, conductor)| EddyConductor {
                                triangle,
                                region: regions[&triangle],
                                conductor,
                            })
                            .collect(),
                        impressed_jz_a_m2,
                        assembly,
                        field_controls,
                        &mut || (*cancel.borrow_mut())(),
                    )
                    .map_err(induction)?;
                let mut fields = Vec::with_capacity(bindings.len());
                let mut internal = [0.0; 2];
                let mut reversible_heat_j = 0.0;
                for cell in &proposal.step().magnetic_cells {
                    check(&cancel)?;
                    if let Some(before) = initial.get(&cell.triangle) {
                        let binding = laws[&cell.triangle];
                        let b = magnitude(cell.flux_density_t[1]);
                        if !b.is_finite() || b > binding.maximum_flux_density_t {
                            return Err(invalid(
                                "proposed magnetic field exceeds its thermal capacity envelope",
                            ));
                        }
                        let at = binding
                            .law
                            .at_temperature(trial.temperature_k(before.cell)?)
                            .map_err(field_error)?;
                        let after = at.state(b, before.volume_m3).map_err(field_error)?;
                        if (after.free_energy_j - cell.free_energy_j[1]).abs()
                            > controls.total_energy_tolerance_j
                            || (before.state.free_energy_j - cell.free_energy_j[0]).abs()
                                > controls.total_energy_tolerance_j
                        {
                            return Err(invalid(
                                "native local magnetic energy and sourced thermodynamic law disagree",
                            ));
                        }
                        let delta_ts = finite(
                            at.temperature_k() * after.entropy_j_k
                                - before.temperature_k * before.state.entropy_j_k,
                        )?;
                        let heat = finite(-cell.parameter_work_j - delta_ts)?;
                        let power = finite(heat / dt)?;
                        if heat != 0.0 && power == 0.0 {
                            return Err(invalid(
                                "magnetic reversible heat disappeared during power conversion",
                            ));
                        }
                        let destination = coefficients.sources_w.entry(before.cell).or_insert(0.0);
                        *destination = finite(*destination + power)?;
                        reversible_heat_j = finite(reversible_heat_j + power * dt)?;
                        internal[0] = finite(internal[0] + before.state.internal_energy_j)?;
                        internal[1] = finite(internal[1] + after.internal_energy_j)?;
                        fields.push(MagneticCaloricSample {
                            triangle: cell.triangle,
                            cell: before.cell,
                            temperature_k: [before.temperature_k, at.temperature_k()],
                            flux_density_t: [before.b_t, b],
                            state: [before.state, after],
                            parameter_work_j: cell.parameter_work_j,
                            reversible_heat_j: heat,
                        });
                    } else {
                        if cell.parameter_work_j != 0.0 {
                            return Err(invalid(
                                "fixed region acquired unowned magnetic parameter work",
                            ));
                        }
                        for (target, value) in internal.iter_mut().zip(cell.free_energy_j) {
                            *target = finite(*target + value)?;
                        }
                    }
                }
                Ok(DrivenCoefficients {
                    circuit: proposal,
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    circuit_trials: 1,
                    source: EddyCaloricSource {
                        material: EddyThermalSample {
                            network: coefficients.source,
                            magnetic,
                        },
                        fields,
                        field_internal_energy_j: internal,
                        reversible_heat_j,
                    },
                })
            },
            |proposal: &PlanarEddyCandidate<'_>, dt, checkpoint| {
                deposit_heat(proposal, &routes, dt, checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let external_heat_j =
            finite(report.non_electrical_heat_j - report.source.reversible_heat_j)?;
        let total_energy_residual_j = finite(
            report.source.field_internal_energy_j[1] - report.source.field_internal_energy_j[0]
                + report.thermal.stored_energy_change_j
                - report.electrical.supplied_work_j
                - external_heat_j,
        )?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(invalid(
                "physical magnetic and caloric internal energy do not close external supply",
            ));
        }
        check(&cancel)?;
        Ok(EddyCaloricCandidate {
            thermal,
            field,
            report: EddyCaloricReport {
                coupled: report,
                external_heat_j,
                total_energy_residual_j,
            },
        })
    }
}
