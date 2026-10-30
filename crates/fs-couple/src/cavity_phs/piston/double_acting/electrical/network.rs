//! Electrical force/back-emf AND finite material heat on one shared mechanism.
//! Reuses the original circuit/enthalpy residual, paired thermal boundary and
//! electrical-motion fixed point. No resistor loss remains an unowned heat sink.

use super::*;
use crate::cavity_phs::piston::thermal::ThermalPistonError;
use crate::cavity_phs::piston::thermal::boundary::double_acting::{DoubleBoundaryCandidate, DoubleBoundaryStep};
use crate::cavity_phs::piston::thermal::boundary::double_acting::network::DoubleNetworkDrive;
use crate::cavity_phs::piston::thermal::boundary::network::{NetworkGasSolve, NetworkGasStep};
use crate::cavity_phs::piston::thermal::boundary::network::material::{GasConductionGeometry, resolve_contact};
use crate::circuit_enthalpy::{self, CircuitEnthalpyCandidate, CircuitEnthalpyError, CircuitEnthalpyInput,
    CircuitEnthalpyReport, DrivenCoefficients};
use crate::circuit_enthalpy::material::{MaterialNetwork, MaterialNetworkError, MaterialNetworkSample};
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_thermal::enthalpy::{DepositionSimulation, DepositionError};
use std::collections::BTreeMap;

/// Two gas contacts and one reciprocal transducer over the same physical step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectromagneticNetworkDrive {
    /// Original paired thermal-boundary, gas and mechanical inputs/controls.
    /// Its additional applied force EXCLUDES the electromagnetic force.
    pub gas: DoubleNetworkDrive,
    /// One exclusive voltage source, reserved with a zero placeholder input.
    pub transducer: TranslationalCircuitPort,
    /// Original coupled velocity controls. Electrical tolerances must equal the
    /// network's electrical controls; the remaining global budget caps each call.
    pub motion: ElectricalMotionSolve,
    /// Explicit starting guess for EVERY velocity solve, never reused trial history.
    pub initial_midpoint_velocity_m_s: f64,
    /// Original material gas/cell conductivity patches [primary, opposing].
    pub geometry: [GasConductionGeometry; 2],
}

/// Material failures and fully typed motion/thermal failures remain distinct.
#[derive(Debug)]
pub enum DrivenSourceError {
    /// Original cell/resistor/conduction material provider refused a state.
    Material(MaterialNetworkError),
    /// Original circuit, velocity, gas or thermal-boundary solve refused a state.
    Motion(ElectricalMotionError<ThermalPistonError<MaterialNetworkError>>),
}
impl core::fmt::Display for DrivenSourceError {
    fn fmt(&self, f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"driven material network: {self:?}")}
}
impl std::error::Error for DrivenSourceError {}
/// Failed proposals leave circuit, cells, both gases and momentum untouched.
pub type DrivenNetworkError = CircuitEnthalpyError<DrivenSourceError>;

/// Original material and gas/motion proposal at the admitted thermal endpoint.
pub struct DrivenNetworkSource {
    /// Original cells, phase, resistivity and conductivity source evaluations.
    pub material: MaterialNetworkSample,
    /// Exact paired gas proposal and both original interface source evaluations.
    pub boundary: DoubleBoundaryCandidate<ResolvedMaterialStatePoint>,
    /// Reciprocal native electrical/mechanical port work.
    pub motion: ElectricalMotionStep,
    /// Independent native thermal rate before inserting Joule or gas transfers [J].
    pub external_heat_j: f64,
}
/// All finite owners proposed together, using the original network format.
pub struct DrivenNetworkCandidate<'m> {
    network: CircuitEnthalpyCandidate<'m, DrivenNetworkSource>,
    coupled: NetworkGasStep,
}
/// Reports from the original electrical, thermal and gas owners.
pub struct DrivenNetworkReport {
    /// Original network report with the original material source payload.
    pub network: CircuitEnthalpyReport<MaterialNetworkSample>,
    /// Two thermal contacts, both gases and one kinetic-energy account.
    pub boundary: DoubleBoundaryStep,
    /// Original conductivity evaluations, not reconstructed material constants.
    pub interfaces: [ResolvedMaterialStatePoint; 2],
    /// Native force/back-emf/current and reciprocal work.
    pub motion: ElectricalMotionStep,
    /// Independent whole-system energy and complete nested work counts.
    pub coupled: NetworkGasStep,
}
impl<'m> DrivenNetworkCandidate<'m> {
    /// Complete original electrical/thermal proposal and matching gas source.
    pub const fn network(&self)->&CircuitEnthalpyCandidate<'m,DrivenNetworkSource> {&self.network}
    /// Exact tentative shared mechanism in this network's final source evaluation.
    pub fn piston(&self)->&DoubleActingGasPiston {self.network.report().source.boundary.piston()}
    /// Independent combined energy, with internal transfers excluded from supply.
    pub const fn coupled(&self)->&NetworkGasStep {&self.coupled}
    /// Publish all owners together. Neither Joule nor boundary heat is reapplied.
    pub fn into_parts(self)->(DepositionSimulation<'m>,MaterialCircuit,DoubleActingGasPiston,DrivenNetworkReport) {
        let (thermal,circuit,r)=self.network.into_parts();
        let DrivenNetworkSource {material,boundary,motion,..}=r.source;
        let (piston,interfaces,boundary)=boundary.into_parts();
        let network=CircuitEnthalpyReport {source:material,electrical:r.electrical,thermal:r.thermal,
            nonlinear:r.nonlinear,maximum_cell_residual_j:r.maximum_cell_residual_j,circuit_trials:r.circuit_trials,
            deposited_joule_j:r.deposited_joule_j,non_electrical_heat_j:r.non_electrical_heat_j,
            total_energy_residual_j:r.total_energy_residual_j};
        (thermal,circuit,piston,DrivenNetworkReport {network,boundary,interfaces,motion,coupled:self.coupled})
    }
}
fn material_error(e:MaterialNetworkError)->DrivenNetworkError {
    match e {
        MaterialNetworkError::Thermal(DepositionError::Cancelled)=>CircuitEnthalpyError::Cancelled,
        other=>CircuitEnthalpyError::Source(DrivenSourceError::Material(other)),
    }
}
fn motion_error(e:ElectricalMotionError<ThermalPistonError<MaterialNetworkError>>)->DrivenNetworkError {
    match e {
        ElectricalMotionError::Cancelled | ElectricalMotionError::Mechanical(ThermalPistonError::Cancelled)=>CircuitEnthalpyError::Cancelled,
        other=>CircuitEnthalpyError::Source(DrivenSourceError::Motion(other)),
    }
}

impl DoubleActingGasPiston {
    /// Couple direct electromagnetic work and both gas heat transfers to the
    /// existing finite circuit/material network. For each thermal trial the
    /// SAME endpoint resistance supplies circuit current, Joule heat and force.
    /// Velocity back emf and both heat laws close before the network residual
    /// sees their paired proposal. Every nested trial starts from original history.
    ///
    /// All resistor heat is deposited by the original circuit/enthalpy owner.
    /// Contact heats debit their named finite cells once, including a shared cell.
    /// Native external supply is evaluated before inserting internal transfers.
    /// Accepted output contains circuit, cells, both gases and the common motion.
    /// Fixed thermal geometry/C/L/topology, frozen linear field coupling, uniform
    /// closed gases; endpoint thermal rates imply generally first-order accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn electromagnetic_network_candidate<'m>(
        &self,binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
        circuit:&MaterialCircuit,inputs:&[f64],drive:ElectromagneticNetworkDrive,
        controls:NetworkGasSolve,guess_j:&[f64],cancelled:impl FnMut()->bool,
    )->Result<DrivenNetworkCandidate<'m>,DrivenNetworkError> {
        self.electromagnetic_network_candidate_with_midpoint_resistance(binding,thermal,target_s,
            circuit,inputs,drive,controls,guess_j,0.0,cancelled)
    }

    /// Shared acoustic/mechanical composition: retain the same endpoint
    /// resistance, heat and velocity solves, with a native affine load on motion.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::cavity_phs::piston) fn electromagnetic_network_candidate_with_midpoint_resistance<'m>(
        &self,binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
        circuit:&MaterialCircuit,inputs:&[f64],drive:ElectromagneticNetworkDrive,
        controls:NetworkGasSolve,guess_j:&[f64],resistance_n_s_m:f64,cancelled:impl FnMut()->bool,
    )->Result<DrivenNetworkCandidate<'m>,DrivenNetworkError> {
        let cancel=RefCell::new(cancelled);
        let check=||if (*cancel.borrow_mut())(){Err(DepositionError::Cancelled)}else{Ok(())};
        check()?;
        let dt=target_s-thermal.time_s();
        if !positive(dt) || drive.gas.input.duration_s!=dt || controls.maximum_boundary_solves==0
            || !positive(controls.total_energy_tolerance_j) || drive.motion.circuit!=controls.network.circuit
            || binding.resistors.len()>crate::material_circuit::MAX_MATERIAL_CIRCUIT_BRANCHES
            || !resistance_n_s_m.is_finite() || resistance_n_s_m<0.0 {
            return Err(CircuitEnthalpyError::Invalid("one clock, one electrical policy and bounded coupled work required"));
        }
        let mut coverage=BTreeMap::new();
        for (cell,patch) in drive.gas.cells.into_iter().zip(drive.geometry) {
            check()?;
            if !patch.area_m2.is_finite() || patch.area_m2<0.0 || !positive(patch.gas_distance_m) || !positive(patch.cell_distance_m) {
                return Err(CircuitEnthalpyError::Invalid("invalid material gas-contact geometry"));
            }
            let area=coverage.entry(cell).or_insert(0.0);*area+=patch.area_m2;
            if !area.is_finite() || *area>thermal.exposed_area_m2(cell)? {
                return Err(CircuitEnthalpyError::Invalid("paired patches exceed one cell's exposed surface"));
            }
        }
        let destinations:Vec<_>=binding.resistors.iter().map(|r|r.destination).collect();
        let calls=Cell::new(0usize);let heat_trials=Cell::new(0usize);let gas_trials=Cell::new(0usize);
        let candidate=circuit_enthalpy::driven_candidate(thermal,target_s,
            CircuitEnthalpyInput {circuit,ports:inputs,destinations:&destinations},controls.network,guess_j,
            |state,remaining| {
                let mut coefficients=binding.resolve(state,||(*cancel.borrow_mut())()).map_err(material_error)?;
                let external_heat_j=dt*state.external_step_power_w(&coefficients.sources_w,&coefficients.conductances_w_k,check)?;
                let temperatures=[state.temperature_k(drive.gas.cells[0])?,state.temperature_k(drive.gas.cells[1])?];
                let motion_controls=ElectricalMotionSolve {maximum_trials:drive.motion.maximum_trials.min(remaining),..drive.motion};
                let coupled=super::solve(ElectricalMotionInput {circuit,resistors:&coefficients.resistors,inputs,port:drive.transducer},
                    dt,motion_controls,drive.initial_midpoint_velocity_m_s,&cancel,|force| {
                        if calls.get()>=controls.maximum_boundary_solves {
                            return Err(ThermalPistonError::Invalid("global electromagnetic boundary budget exhausted"));
                        }
                        calls.set(calls.get()+1);
                        let total_force=drive.gas.input.external_force_n+force;
                        if !total_force.is_finite() || (force!=0.0 && total_force==drive.gas.input.external_force_n) {
                            return Err(ThermalPistonError::Invalid("electromagnetic force vanished in boundary load"));
                        }
                        let boundary=self.thermal_boundaries_candidate_with_midpoint_resistance(DoublePistonInput {external_force_n:total_force,..drive.gas.input},
                            drive.gas.piston,drive.gas.exchange,temperatures,resistance_n_s_m,|gases|Ok([
                                resolve_contact(binding,drive.gas.cells[0],state,&coefficients.source,&gases[0],drive.geometry[0])?,
                                resolve_contact(binding,drive.gas.cells[1],state,&coefficients.source,&gases[1],drive.geometry[1])?,
                            ]),||(*cancel.borrow_mut())())?;
                        heat_trials.set(heat_trials.get().checked_add(boundary.step().joint_trials)
                            .ok_or(ThermalPistonError::Invalid("boundary work count overflow"))?);
                        gas_trials.set(gas_trials.get().checked_add(boundary.step().material_trials)
                            .ok_or(ThermalPistonError::Invalid("material work count overflow"))?);
                        let p=&boundary.step().piston;
                        Ok(MotionTrial {displacement_m:p.displacement_m,
                            midpoint_velocity_m_s:f64::midpoint(p.velocity_m_s[0],p.velocity_m_s[1]),state:boundary})
                    }).map_err(motion_error)?;
                let mut debits=BTreeMap::new();
                for (cell,q) in drive.gas.cells.into_iter().zip(coupled.mechanical.step().exchanged_heat_j) {
                    let sum=debits.entry(cell).or_insert(0.0);*sum+=q;
                }
                let tolerance=drive.gas.exchange.faces[0].energy_tolerance_j.min(drive.gas.exchange.faces[1].energy_tolerance_j);
                for (cell,q) in debits {
                    check()?;
                    let watts=q/dt;let source=coefficients.sources_w.entry(cell).or_insert(0.0);
                    let previous=*source;let next=previous-watts;let actual=(previous-next)*dt;
                    if !q.is_finite() || !watts.is_finite() || !next.is_finite() || !actual.is_finite()
                        || (q!=0.0 && (watts==0.0 || next==previous)) || (actual-q).abs()>tolerance {
                        return Err(CircuitEnthalpyError::Invalid("actual finite-cell heat debit disagrees with gas transfer"));
                    }
                    *source=next;
                }
                let circuit_trials=coupled.step.trials;
                Ok(DrivenCoefficients {circuit:coupled.circuit,sources_w:coefficients.sources_w,
                    conductances_w_k:coefficients.conductances_w_k,circuit_trials,
                    source:DrivenNetworkSource {material:coefficients.source,boundary:coupled.mechanical,motion:coupled.step,external_heat_j}})
            },||(*cancel.borrow_mut())())?;
        let r=candidate.report();let source=&r.source;let p=&source.boundary.step().piston;
        // The additional load is outside this native circuit/thermal boundary.
        // Keep its ACTUAL mechanical work so a surrounding acoustic account can
        // cancel it. The zero-resistance path retains its original arithmetic.
        let mechanical_work=if resistance_n_s_m==0.0 {
            drive.gas.input.external_force_n*p.displacement_m
        } else {
            p.external_force_work_j-source.motion.mechanical_work_j
        };
        let total=(r.electrical.stored_after_j-r.electrical.stored_before_j)+r.thermal.stored_energy_change_j
            +p.gas[0].stored_energy_change_j+p.gas[1].stored_energy_change_j+p.kinetic_energy_change_j
            -source.motion.external_electrical_work_j-source.external_heat_j-drive.gas.input.heat_j[0]-drive.gas.input.heat_j[1]
            -mechanical_work;
        if !source.external_heat_j.is_finite() || !total.is_finite() || total.abs()>controls.total_energy_tolerance_j {
            return Err(CircuitEnthalpyError::Invalid("electrical, finite thermal, gas and mechanical storage do not close external supply"));
        }
        let coupled=NetworkGasStep {external_network_heat_j:source.external_heat_j,total_energy_residual_j:total,
            boundary_solves:calls.get(),exchange_trials:heat_trials.get(),material_trials:gas_trials.get()};
        check()?;
        Ok(DrivenNetworkCandidate {network:candidate,coupled})
    }
}

use crate::cavity_phs::piston::material::MaterialDoubleActingGasPiston;

/// Complete finite network with the original material-derived mechanical inertia.
pub struct MaterialDrivenNetworkCandidate<'m> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    piston: MaterialDoubleActingGasPiston,
    report: DrivenNetworkReport,
}
impl<'m> MaterialDrivenNetworkCandidate<'m> {
    /// Attach the unchanged material body to the complete proposed dynamics.
    /// Keeping the owner intact preserves its accepted port clocks and history.
    pub(in crate::cavity_phs::piston) fn from_native(
        original:&MaterialDoubleActingGasPiston,candidate:DrivenNetworkCandidate<'m>,
    )->Result<Self,DrivenNetworkError> {
        let (thermal,circuit,next,report)=candidate.into_parts();
        let piston=original.bind_piston_candidate(next)
            .map_err(|e|motion_error(ElectricalMotionError::Mechanical(ThermalPistonError::from(e))))?;
        Ok(Self {thermal,circuit,piston,report})
    }
    /// Original finite thermal timeline at the proposed endpoint.
    pub const fn thermal(&self)->&DepositionSimulation<'m> {&self.thermal}
    /// Original electrical graph and accepted-candidate charges/fluxes.
    pub const fn circuit(&self)->&MaterialCircuit {&self.circuit}
    /// Both gas histories, shared momentum and canonical material-derived inertia.
    pub const fn piston(&self)->&MaterialDoubleActingGasPiston {&self.piston}
    /// All physical work and original material/interface evaluations.
    pub const fn report(&self)->&DrivenNetworkReport {&self.report}
    /// Publish all endpoints together; no second circuit step or heat deposit.
    pub fn into_parts(self)->(DepositionSimulation<'m>,MaterialCircuit,MaterialDoubleActingGasPiston,DrivenNetworkReport) {
        (self.thermal,self.circuit,self.piston,self.report)
    }
}
impl MaterialDoubleActingGasPiston {
    /// Preserve regional density/geometry sources through the same electromagnetic,
    /// thermal and gas solve. The rigid specimen supplies the only inertial mass.
    /// Its thermal state remains frozen; network cells may represent other bodies.
    #[allow(clippy::too_many_arguments)]
    pub fn electromagnetic_network_candidate<'m>(
        &self,binding:&MaterialNetwork<'_>,thermal:&DepositionSimulation<'m>,target_s:f64,
        circuit:&MaterialCircuit,inputs:&[f64],drive:ElectromagneticNetworkDrive,
        controls:NetworkGasSolve,guess_j:&[f64],mut cancelled:impl FnMut()->bool,
    )->Result<MaterialDrivenNetworkCandidate<'m>,DrivenNetworkError> {
        let candidate=self.piston().electromagnetic_network_candidate(binding,thermal,target_s,circuit,inputs,
            drive,controls,guess_j,&mut cancelled)?;
        let candidate=MaterialDrivenNetworkCandidate::from_native(self,candidate)?;
        if cancelled(){return Err(CircuitEnthalpyError::Cancelled);}
        Ok(candidate)
    }
}

#[cfg(test)]
pub(in crate::cavity_phs::piston::double_acting) mod tests;
