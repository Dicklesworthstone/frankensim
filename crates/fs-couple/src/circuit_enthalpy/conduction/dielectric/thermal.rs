//! Native spatial charge, affine epsilon(T) and caloric enthalpy at one endpoint.
//!
//! The original field receipt supplies psi and signed parameter work P. Its
//! original derivative law supplies S and U=psi+TS; Qrev=-P-delta(TS) joins the
//! same thermal residual separately from original Ohmic heat. The enclosing
//! audit uses field INTERNAL energy, not a second credit for electric storage.

use super::*;
use fs_material::field_storage::{FieldStorageError, thermal::FieldThermalState};

#[derive(Clone, Debug, PartialEq)]
pub struct DielectricCaloricSample {
    pub triangle: usize,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    pub displacement_c_m2: [[f64; 2]; 2],
    pub state: [FieldThermalState; 2],
    pub parameter_work_j: f64,
    pub reversible_heat_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DielectricCaloricSource {
    pub material: MaterialNetworkSample,
    pub fields: Vec<DielectricCaloricSample>,
    pub field_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}

pub struct DielectricCaloricReport {
    /// Native Helmholtz/parameter-work account and original enthalpy residual.
    pub coupled: CircuitEnthalpyReport<DielectricCaloricSource, DielectricRelaxationStep>,
    /// Original prescribed/environmental heat, excluding internal Qrev/Joule.
    pub external_heat_j: f64,
    pub total_energy_residual_j: f64,
}

pub struct DielectricCaloricCandidate<'t, C = ResolvedConductorState> {
    thermal: DepositionSimulation<'t>,
    field: PlanarLeakyDielectric<C>,
    report: DielectricCaloricReport,
}
impl<'t, C> DielectricCaloricCandidate<'t, C> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn field(&self) -> &PlanarLeakyDielectric<C> { &self.field }
    pub fn report(&self) -> &DielectricCaloricReport { &self.report }
    /// Both histories publish together; all thermal transfers are already paid.
    pub fn into_parts(self) -> (DepositionSimulation<'t>, PlanarLeakyDielectric<C>, DielectricCaloricReport) {
        (self.thermal, self.field, self.report)
    }
}

struct Initial {
    cell: CellId,
    temperature_k: f64,
    volume_m3: f64,
    displacement_c_m2: [f64; 2],
    state: FieldThermalState,
}
fn source_error(error: FieldStorageError) -> Error {
    dielectric(DielectricRelaxationError::Electric(fs_em::electrostatic::ElectrostaticError::Material(error)))
}

impl ConductionThermalNetwork<'_> {
    /// Close the original native charge and heat residual with source-backed
    /// scalar/oriented epsilon(T), actual sigma(T), and local field entropy.
    /// `maximum_displacement_c_m2` is an ordered |D| envelope per original
    /// triangle, enforced at both accepted endpoints. All triangles sharing a
    /// receiver contribute to its combined-capacity admission.
    ///
    /// Mapped caloric charts must be fully solid at constant density, within
    /// complete affine source intervals. Geometry and material axes stay fixed.
    /// The existing constant-epsilon API separately retains its phase support.
    /// This first-order thermal coupling is not an unconditional entropy or
    /// time-accuracy certificate; dispersion, phase polarization and moving
    /// interfaces require their own constitutive/free-energy laws.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dielectric_thermal_candidate<'t, C: PlanarConductorMaterial>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarLeakyDielectric<C>,
        endpoint_voltage_v: &[(usize, f64)],
        maximum_displacement_c_m2: &[f64],
        assembly: WeightedAssemblyLimits,
        field_controls: DielectricRelaxationControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<DielectricCaloricCandidate<'t, C>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let laws = field.thermal_dielectric_laws()
            .ok_or_else(|| invalid("native dielectric must retain its original thermal derivative laws"))?;
        if field.time_s().to_bits() != thermal.time_s().to_bits()
            || self.cells.len() != thermal.cells().len()
            || laws.len() != maximum_displacement_c_m2.len()
        {
            return Err(invalid("matching dielectric/thermal clocks, complete cells and one displacement envelope per original triangle required"));
        }
        let network = MaterialNetwork { cells: self.cells, resistors: &[], links: self.links };
        let prepared = network.prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let problem = field.conduction();
        let routes = self.admit_routes(thermal, problem, &cancel)?;
        let sample = prepared.resolve(thermal, || (*cancel.borrow_mut())()).map_err(material)?;
        let sources = self.materials_at(problem, &sample.source, &routes, &cancel)?;
        if sources.iter().zip(problem.conductors()).any(|(a, b)| {
            !a.has_same_electrical_source(b) || a.frame() != b.frame()
                || !same_material(a.material(), b.material())
        }) {
            return Err(invalid("accepted dielectric sources must match the actual caloric state"));
        }
        let displacement = field.displacement_c_m2(&mut || (*cancel.borrow_mut())()).map_err(dielectric)?;
        let bindings: BTreeMap<_, _> = self.cells.iter().map(|s| (s.cell, s)).collect();
        let mut initial = Vec::with_capacity(laws.len());
        let mut capacities = BTreeMap::<CellId, f64>::new();
        let mut ranges = BTreeMap::<CellId, [f64; 2]>::new();
        for (triangle, law) in laws.iter().enumerate() {
            check(&cancel)?;
            let cell = routes[triangle];
            let binding = bindings[&cell];
            let temperature_k = thermal.temperature_k(cell)?;
            let resolved = law.at_material(problem.conductors()[triangle].material()).map_err(source_error)?;
            if binding.point.material.content_hash() != law.card().content_hash()
                || binding.point.temperature_axis != law.temperature_axis()
                || resolved.temperature_k().to_bits() != temperature_k.to_bits()
            {
                return Err(invalid("dielectric, derivative law and caloric receiver must share the same original material and temperature coordinate"));
            }
            let range = law.temperature_range_k();
            let knots = binding.curve.knots();
            if knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
                || range[0] < knots[0].temperature_k
                || range[1] > knots[knots.len()-1].temperature_k
            {
                return Err(invalid("dielectric thermal work requires a fully solid caloric chart covering its complete source interval"));
            }
            let maximum = maximum_displacement_c_m2[triangle];
            let d = displacement[triangle];
            if !maximum.is_finite() || maximum < 0.0 || !d[0].hypot(d[1]).is_finite()
                || d[0].hypot(d[1]) > maximum
            {
                return Err(invalid("accepted dielectric displacement exceeds its finite nonnegative capacity envelope"));
            }
            let volume = problem.mesh().face_measure(triangle)
                .ok_or_else(|| invalid("dielectric triangle lacks original physical volume"))?;
            let state = resolved.state(d, volume).map_err(source_error)?;
            if finite(state.free_energy_j - field.cell_energy_j()[triangle])?.abs()
                > controls.total_energy_tolerance_j
            {
                return Err(invalid("original native dielectric storage disagrees with its thermodynamic law"));
            }
            let contribution = law.minimum_heat_capacity_j_k(maximum, volume).map_err(source_error)?;
            let capacity = capacities.entry(cell).or_insert(0.0);
            *capacity = finite(*capacity + contribution)?;
            let common = ranges.entry(cell).or_insert(range);
            common[0] = common[0].max(range[0]);
            common[1] = common[1].min(range[1]);
            initial.push(Initial { cell, temperature_k, volume_m3: volume, displacement_c_m2: d, state });
        }
        for (&cell, &range) in &ranges {
            check(&cancel)?;
            if range[0] >= range[1] {
                return Err(invalid("dielectric cell requires a common nonempty thermal source interval"));
            }
            let mut minimum = f64::INFINITY;
            for pair in bindings[&cell].curve.knots().windows(2) {
                check(&cancel)?;
                let dt = pair[1].temperature_k - pair[0].temperature_k;
                if dt <= 0.0 {
                    return Err(invalid("dielectric caloric chart requires strictly increasing temperature"));
                }
                if pair[1].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                    minimum = minimum.min(finite(thermal.cells()[&cell].mass_kg
                        * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg) / dt)?);
                }
            }
            if !minimum.is_finite() || finite(minimum + capacities[&cell])? <= 0.0 {
                return Err(invalid("caloric plus all dielectric field heat capacities must stay positive over the declared envelope"));
            }
        }
        let dt = target_s - thermal.time_s();
        let (thermal, field, report) = cell_driver_candidate(
            thermal, target_s, controls, initial_energy_increments_j,
            |trial, _| {
                check(&cancel)?;
                let mut coefficients = prepared.resolve(trial, || (*cancel.borrow_mut())()).map_err(material)?;
                let sources = self.materials_at(problem, &coefficients.source, &routes, &cancel)?;
                let proposal = field.step_to_candidate(target_s, sources, endpoint_voltage_v,
                    assembly, field_controls, &mut || (*cancel.borrow_mut())()).map_err(dielectric)?;
                if proposal.step().electric_cells.len() != initial.len() {
                    return Err(invalid("native dielectric receipt lost original triangle correspondence"));
                }
                let mut fields = Vec::with_capacity(initial.len());
                let mut field_internal_energy_j = [0.0; 2];
                let mut reversible_heat_j = 0.0;
                for (triangle, before) in initial.iter().enumerate() {
                    check(&cancel)?;
                    let receipt = &proposal.step().electric_cells[triangle];
                    let d = receipt.displacement_c_m2[1];
                    if receipt.triangle != triangle || !d[0].hypot(d[1]).is_finite()
                        || d[0].hypot(d[1]) > maximum_displacement_c_m2[triangle]
                    {
                        return Err(invalid("proposed dielectric displacement exceeds its original capacity envelope"));
                    }
                    let at = laws[triangle].at_material(proposal.field().conduction().conductors()[triangle].material())
                        .map_err(source_error)?;
                    if at.temperature_k().to_bits() != trial.temperature_k(before.cell)?.to_bits() {
                        return Err(invalid("dielectric endpoint and caloric trial temperatures disagree"));
                    }
                    let after = at.state(d, before.volume_m3).map_err(source_error)?;
                    if finite(after.free_energy_j - receipt.free_energy_j[1])?.abs() > controls.total_energy_tolerance_j
                        || finite(before.state.free_energy_j - receipt.free_energy_j[0])?.abs() > controls.total_energy_tolerance_j
                    {
                        return Err(invalid("native dielectric endpoint energy and original derivative law disagree"));
                    }
                    let delta_ts = finite(at.temperature_k() * after.entropy_j_k
                        - before.temperature_k * before.state.entropy_j_k)?;
                    let heat = finite(-receipt.parameter_work_j - delta_ts)?;
                    let power = finite(heat / dt)?;
                    if heat != 0.0 && power == 0.0 {
                        return Err(invalid("dielectric reversible heat disappeared during power conversion"));
                    }
                    let source = coefficients.sources_w.entry(before.cell).or_insert(0.0);
                    let next = finite(*source + power)?;
                    if power != 0.0 && next == *source {
                        return Err(invalid("dielectric reversible heat disappeared in its receiving cell"));
                    }
                    *source = next;
                    reversible_heat_j = finite(reversible_heat_j + power * dt)?;
                    field_internal_energy_j[0] = finite(field_internal_energy_j[0] + before.state.internal_energy_j)?;
                    field_internal_energy_j[1] = finite(field_internal_energy_j[1] + after.internal_energy_j)?;
                    fields.push(DielectricCaloricSample {
                        triangle, cell: before.cell, temperature_k: [before.temperature_k, at.temperature_k()],
                        displacement_c_m2: [before.displacement_c_m2, d], state: [before.state, after],
                        parameter_work_j: receipt.parameter_work_j, reversible_heat_j: heat,
                    });
                }
                Ok(DrivenCoefficients { circuit: proposal,
                    sources_w: coefficients.sources_w, conductances_w_k: coefficients.conductances_w_k,
                    circuit_trials: 1, source: DielectricCaloricSource {
                        material: coefficients.source, fields, field_internal_energy_j, reversible_heat_j,
                    },
                })
            },
            |proposal: &DielectricRelaxationCandidate<C>, _dt, checkpoint| {
                deposit_heat(proposal, &routes, checkpoint)
            },
            || (*cancel.borrow_mut())(),
        )?;
        let external_heat_j = finite(report.non_electrical_heat_j - report.source.reversible_heat_j)?;
        let total_energy_residual_j = finite(report.source.field_internal_energy_j[1]
            - report.source.field_internal_energy_j[0] + report.thermal.stored_energy_change_j
            - report.electrical.supplied_work_j - external_heat_j)?;
        if total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(invalid("physical dielectric and caloric internal energy do not close external supply"));
        }
        check(&cancel)?;
        Ok(DielectricCaloricCandidate { thermal, field,
            report: DielectricCaloricReport { coupled: report, external_heat_j, total_energy_residual_j },
        })
    }
}
