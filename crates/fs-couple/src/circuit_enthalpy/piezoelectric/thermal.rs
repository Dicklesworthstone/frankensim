//! Ceramic temperature feeds the original charge/strain/circuit trajectory.
//! Every enthalpy trial begins with the same accepted rod/circuit coordinates.
//! The native coefficient update exposes P; the source law supplies S and U.
//! Signed Qrev=-P-delta(TS) enters the ceramic caloric cell alongside separately
//! routed native resistor heat. Final acceptance closes physical U, kinetic,
//! other circuit storage and caloric enthalpy against actual external supply.

use super::*;

// Complete sourced heat transport shares the same ceramic/circuit residual.
pub mod heat;

use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use fs_material::{
    conductor::ElectricalPropertyBasis,
    field_storage::thermal::FieldThermalState,
    piezoelectric::{LongitudinalPiezoelectric, PiezoelectricError,
        thermal::ThermalLongitudinalPiezoelectric},
};

/// One whole uniform ceramic body in the original finite thermal network.
pub struct PiezoelectricThermalBinding {
    pub cell: CellId,
    pub law: ThermalLongitudinalPiezoelectric,
    /// Both accepted |Q/A| endpoints must satisfy this bound. Strain uses the
    /// original rod's small-signal bound; neither is a measured strength claim.
    pub maximum_abs_displacement_c_m2: f64,
    pub volume_tolerance_m3: f64,
}

#[derive(Debug)]
pub enum PiezoelectricThermalError<E = crate::circuit_enthalpy::contacts::heat::HeatNetworkError> {
    Network(MaterialNetworkError),
    Heat(E),
    Material(PiezoelectricError),
    Invalid(&'static str),
}
impl<E: core::fmt::Debug> core::fmt::Display for PiezoelectricThermalError<E> {
    fn fmt(&self,f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"piezoelectric ceramic/enthalpy: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for PiezoelectricThermalError<E> {}
type Error<E = crate::circuit_enthalpy::contacts::heat::HeatNetworkError> =
    CircuitEnthalpyError<PiezoelectricThermalError<E>>;
fn invalid<E>(message: &'static str) -> Error<E> {
    CircuitEnthalpyError::Source(PiezoelectricThermalError::Invalid(message))
}
fn source<E>(error: PiezoelectricError) -> Error<E> {
    CircuitEnthalpyError::Source(PiezoelectricThermalError::Material(error))
}
fn network<E>(error: MaterialNetworkError) -> Error<E> {
    match error {
        MaterialNetworkError::Thermal(DepositionError::Cancelled) => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(PiezoelectricThermalError::Network(other)),
    }
}
fn check<E>(cancel: &RefCell<impl FnMut() -> bool>) -> Result<(),Error<E>> {
    if (*cancel.borrow_mut())() { Err(CircuitEnthalpyError::Cancelled) } else { Ok(()) }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PiezoelectricThermalSource<S = MaterialNetworkSample> {
    pub network: S,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    /// Original [axial strain, axial D] coordinates at the accepted endpoints.
    pub coordinates: [[f64; 2]; 2],
    pub state: [FieldThermalState; 2],
    pub material: LongitudinalPiezoelectric,
    pub reversible_heat_j: f64,
}
pub struct PiezoelectricThermalReport<S = MaterialNetworkSample> {
    pub coupled: PiezoelectricEnthalpyReport<PiezoelectricThermalSource<S>>,
    /// Authored heat, excluding reversible ceramic exchange and closed radiation.
    pub external_heat_j: f64,
    /// Retained numerical net of the closed enclosure, never external supply.
    pub closed_radiation_net_j: f64,
    pub total_energy_residual_j: f64,
}
pub struct PiezoelectricThermalCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    model: PiezoelectricCircuit,
    report: PiezoelectricThermalReport<S>,
}
impl<'t, S> PiezoelectricThermalCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub fn model(&self) -> &PiezoelectricCircuit { &self.model }
    pub fn report(&self) -> &PiezoelectricThermalReport<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>,PiezoelectricCircuit,PiezoelectricThermalReport<S>) {
        (self.thermal,self.model,self.report)
    }
}

/// Solve ceramic coefficients, rod motion, circuit current, resistor heating
/// and finite material temperatures on the existing enthalpy/Newton residual.
/// The ceramic has a fully solid constant-density chart covering its complete
/// admitted affine source interval and whole reference volume. Its caloric
/// capacity plus the complete coupled-field bound must remain positive.
///
/// Other cells retain the existing thermal/network behavior. This endpoint-T
/// split is first order in general, not unconditional entropy/passivity or a
/// time-error certificate. The native small-signal rod geometry, field bounds,
/// source constraints and circuit limits remain authoritative. Thermal strain,
/// pyroelectric offsets, phase polarization and irreversible ceramic loss are
/// not inferred. Refusal/cancellation publishes none of the three histories.
#[allow(clippy::too_many_arguments)]
pub fn step_candidate<'t>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
    model: &PiezoelectricCircuit, electrical_ports: &[f64], force_n: f64,
    ceramic: &PiezoelectricThermalBinding, electrical_bases: &[(usize,ElectricalPropertyBasis)],
    controls: CircuitEnthalpySolve, guess_j: &[f64], cancelled: impl FnMut() -> bool,
) -> Result<PiezoelectricThermalCandidate<'t>,Error> {
    let cancel=RefCell::new(cancelled); check(&cancel)?;
    super::super::admit_driver_unknowns(thermal,controls,guess_j)?;
    let prepared=material.prepare_with_bases(thermal,electrical_bases,||(*cancel.borrow_mut())()).map_err(network)?;
    resolved_candidate(material,thermal,target_s,model,electrical_ports,force_n,
        ceramic,controls,guess_j,
        |trial| prepared.resolve(trial,||(*cancel.borrow_mut())()).map_err(network),
        |_| 0.0,||(*cancel.borrow_mut())())
}

// Both the original bulk network and the complete heat assembly enter this
// same driver. Only coefficient/source resolution differs; accepted histories,
// ceramic reversible heat and native mechanical/electrical stepping do not.
#[allow(clippy::too_many_arguments)]
fn resolved_candidate<'t, S, E>(
    material: &MaterialNetwork<'_>, thermal: &DepositionSimulation<'t>, target_s: f64,
    model: &PiezoelectricCircuit, electrical_ports: &[f64], force_n: f64,
    ceramic: &PiezoelectricThermalBinding, controls: CircuitEnthalpySolve, guess_j: &[f64],
    resolve: impl Fn(&DepositionSimulation<'t>) -> Result<CircuitThermalCoefficients<S>,Error<E>>,
    radiation_net_w: impl Fn(&S) -> f64, cancelled: impl FnMut() -> bool,
) -> Result<PiezoelectricThermalCandidate<'t, S>,Error<E>> {
    let cancel=RefCell::new(cancelled); check(&cancel)?;
    super::super::admit_driver_unknowns(thermal,controls,guess_j)?;
    let cell=thermal.cells().get(&ceramic.cell).ok_or_else(||invalid("ceramic caloric receiver is inactive"))?;
    let binding=material.cells.iter().find(|s|s.cell==ceramic.cell)
        .ok_or_else(||invalid("ceramic requires its original material/caloric binding"))?;
    let phase=binding.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
        .map_err(|e|network(MaterialNetworkError::Phase(e)))?;
    let point=binding.point.at(phase).map_err(|e|network(MaterialNetworkError::Material(e)))?;
    let initial=ceramic.law.at_material(model.rod().material()).map_err(source)?;
    if binding.point.material.content_hash()!=ceramic.law.card().content_hash()
        || binding.point.temperature_axis!=ceramic.law.temperature_axis()
        || point!=initial.material().material().query()
        || phase.temperature_k().to_bits()!=initial.temperature_k().to_bits() {
        return Err(invalid("ceramic dynamics and caloric receiver require the same original source context and temperature"));
    }
    let volume=finite(model.rod().length_m()*model.rod().area_m2())?;
    if !ceramic.volume_tolerance_m3.is_finite() || ceramic.volume_tolerance_m3<0.0
        || (volume-cell.geometry.volume_m3).abs()>ceramic.volume_tolerance_m3 {
        return Err(invalid("ceramic thermal receiver must cover its whole reference volume"));
    }
    let range=ceramic.law.temperature_range_k(); let knots=binding.curve.knots();
    if range[0]<knots[0].temperature_k || range[1]>knots[knots.len()-1].temperature_k {
        return Err(invalid("ceramic caloric chart must cover its complete thermal source interval"));
    }
    for knot in knots {
        check(&cancel)?;
        if knot.liquid_mass_fraction!=0.0
            || (finite(cell.mass_kg/knot.bulk_density_kg_m3)?-volume).abs()>ceramic.volume_tolerance_m3
            || (finite(cell.mass_kg/model.rod().material().density_kg_m3())?-volume).abs()>ceramic.volume_tolerance_m3 {
            return Err(invalid("thermal piezoelectric motion requires a fully solid chart and matching constant reference density"));
        }
    }
    let x0=[model.extension_m()/model.rod().length_m(),model.charge_c()/model.rod().area_m2()];
    let bound=ceramic.law.minimum_heat_capacity_j_k(model.rod().limits().max_abs_strain,
        ceramic.maximum_abs_displacement_c_m2,volume).map_err(source)?;
    if x0[1].abs()>ceramic.maximum_abs_displacement_c_m2 {
        return Err(invalid("accepted ceramic charge exceeds its thermal capacity envelope"));
    }
    let mut minimum=f64::INFINITY;
    for pair in knots.windows(2) {
        check(&cancel)?; let dt=pair[1].temperature_k-pair[0].temperature_k;
        if dt<=0.0 { return Err(invalid("ceramic caloric chart requires strictly increasing temperature")); }
        if pair[1].temperature_k>=range[0] && pair[0].temperature_k<=range[1] {
            minimum=minimum.min(finite(cell.mass_kg*(pair[1].specific_enthalpy_j_kg-pair[0].specific_enthalpy_j_kg)/dt)?);
        }
    }
    if !minimum.is_finite() || finite(minimum+bound)?<=0.0 {
        return Err(invalid("caloric plus complete piezoelectric field capacity must remain positive over the admitted coordinate envelope"));
    }
    let before=initial.state(x0[0],x0[1],volume).map_err(source)?;
    let kinetic_before=finite(0.5*model.momentum_kg_m_s()*model.velocity_m_s())?;
    let destinations:Vec<_>=material.resistors.iter().map(|r|r.destination).collect();
    let dt=target_s-thermal.time_s();
    let (thermal,model,report)=coupled_driver_candidate(thermal,target_s,CircuitEnthalpyInput {
        circuit:model.circuit(),ports:electrical_ports,destinations:&destinations,
    },&[],controls,guess_j,|trial,_remaining| {
        check(&cancel)?;
        let mut coefficients=resolve(trial)?;
        let at=ceramic.law.at_temperature(trial.temperature_k(ceramic.cell)?).map_err(source)?;
        let candidate=model.step_with_thermal_material_candidate(&ceramic.law,at.temperature_k(),
            coefficients.resistors,electrical_ports,force_n,dt,controls.circuit,||(*cancel.borrow_mut())())?;
        let next=candidate.model();
        let x1=[next.extension_m()/next.rod().length_m(),next.charge_c()/next.rod().area_m2()];
        if x1[1].abs()>ceramic.maximum_abs_displacement_c_m2 {
            return Err(invalid("proposed ceramic charge exceeds its thermal capacity envelope"));
        }
        let after=at.state(x1[0],x1[1],volume).map_err(source)?;
        let frozen=at.state(x0[0],x0[1],volume).map_err(source)?;
        let parameter=finite(frozen.free_energy_j-before.free_energy_j)?;
        let kinetic_after=finite(0.5*next.momentum_kg_m_s()*next.velocity_m_s())?;
        if finite(parameter-candidate.step().parameter_work_j)?.abs()>controls.total_energy_tolerance_j
            || finite(before.free_energy_j+kinetic_before-candidate.step().piezoelectric_energy_j[0])?.abs()>controls.total_energy_tolerance_j
            || finite(after.free_energy_j+kinetic_after-candidate.step().piezoelectric_energy_j[1])?.abs()>controls.total_energy_tolerance_j {
            return Err(invalid("native piezoelectric storage/parameter work disagrees with its thermodynamic law"));
        }
        let heat=finite(-parameter-(at.temperature_k()*after.entropy_j_k-initial.temperature_k()*before.entropy_j_k))?;
        let power=finite(heat/dt)?;
        if heat!=0.0 && power==0.0 { return Err(invalid("ceramic reversible heat underflows power conversion")); }
        let destination=coefficients.sources_w.entry(ceramic.cell).or_insert(0.0);
        *destination=finite(*destination+power)?;
        Ok(DrivenCoefficients { circuit:candidate,sources_w:coefficients.sources_w,
            conductances_w_k:coefficients.conductances_w_k,circuit_trials:1,
            source:PiezoelectricThermalSource { network:coefficients.source,cell:ceramic.cell,
                temperature_k:[initial.temperature_k(),at.temperature_k()],coordinates:[x0,x1],
                state:[before,after],material:at.material().clone(),reversible_heat_j:finite(power*dt)?,
            },
        })
    },||(*cancel.borrow_mut())())?;
    let closed_radiation_net_j=finite(dt*radiation_net_w(&report.source.network))?;
    let external_heat_j=finite(report.non_electrical_heat_j-report.source.reversible_heat_j
        -closed_radiation_net_j)?;
    let kinetic_after=finite(0.5*model.momentum_kg_m_s()*model.velocity_m_s())?;
    let total_energy_residual_j=finite(report.source.state[1].internal_energy_j-before.internal_energy_j
        +kinetic_after-kinetic_before
        +report.electrical.other_electrical_energy_j[1]-report.electrical.other_electrical_energy_j[0]
        +report.thermal.stored_energy_change_j-report.electrical.electrical_supply_j
        -report.electrical.mechanical_supply_j-external_heat_j)?;
    if total_energy_residual_j.abs()>controls.total_energy_tolerance_j {
        return Err(invalid("physical ceramic, kinetic, circuit and caloric storage do not close external supply"));
    }
    check(&cancel)?;
    Ok(PiezoelectricThermalCandidate { thermal,model,
        report:PiezoelectricThermalReport { coupled:report,external_heat_j,closed_radiation_net_j,
            total_energy_residual_j } })
}
