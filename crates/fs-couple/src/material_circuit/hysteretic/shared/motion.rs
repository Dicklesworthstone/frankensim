//! Reciprocal native body motion of one source-owned hysteretic magnetic gap.
//!
//! One joint interface solve closes all winding Faraday equations and the gap
//! displacement. Every trial starts from the same accepted core memory, total
//! linkages and body momentum. The original circuit consumes the material's
//! path-mean currents; the native mechanical descriptor consumes its path-mean
//! gap force. Their internal work cancels, including threshold crossings.
//! Magnetic hysteresis is never presented as a path-independent Storage law.
//!
//! A source-derived guided body and optional frozen axial support share one
//! translational coordinate. Common-flux reversals require separate intervals;
//! three instantaneous rate probes reject detectable turns but do not certify
//! unresolved oscillations. No fringing, contact impact, changing leakage,
//! vector core rotation or thermal evolution is inferred by this owner.

use super::{
    HystereticCircuitError, HystereticCircuitSolve, NativeWindingPorts, Result,
    SharedHystereticCircuit, SharedHystereticStep, finite, instantaneous_voltages, invalid,
    solve_interface,
};
use crate::material_circuit::{
    MaterialCircuitStep, MaterialResistor, motion::support::AxialSupport,
};
use fs_matdb::QueryPoint;
use fs_mbd_material::ResolvedBody;
use fs_phs::{DescriptorPortHamiltonian, PhsError, StepRecord, Storage, step_descriptor};
use fs_solver::NewtonReport;
use std::cell::{Cell, RefCell};

/// Numerical interface scales and independent mechanical acceptance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SharedHystereticGapMotionSolve {
    /// Original electrical/magnetic policy. Each counted joint trial performs
    /// one circuit step and one fixed two-coordinate mechanical descriptor step.
    pub circuit: HystereticCircuitSolve,
    /// Mechanical Jacobian perturbation [m], paired with difference_step_wb.
    /// The interface uses their ratio to scale the displacement coordinate.
    pub difference_step_m: f64,
    pub displacement_tolerance_m: f64,
    /// Independent total and mechanical work allowance [J].
    pub energy_tolerance_j: f64,
}

/// One common magnetic history and source-derived guided body momentum.
#[derive(Clone, Debug, PartialEq)]
pub struct SharedHystereticGapMotion {
    circuit: SharedHystereticCircuit,
    body: ResolvedBody,
    support: Option<AxialSupport>,
    momentum_kg_m_s: f64,
}

/// Native substeps plus the independently closed physical work account.
pub struct SharedHystereticGapMotionStep {
    /// Native graph receipt, including its PRIVATE winding-current ports.
    pub circuit: MaterialCircuitStep,
    /// One common memory, gap/leakage storage and exactly one pinning heat.
    pub magnetic: SharedHystereticStep,
    /// Native mechanical state is [gap, momentum]; its force ports are
    /// [external force, magnetic force], with work-conjugate velocity outputs.
    pub mechanical: StepRecord,
    pub terminal_voltages_v: Vec<f64>,
    pub terminal_voltage_residuals_v: Vec<f64>,
    /// Actual common-flux rate numerator sum(r V/D) - R' Phi velocity.
    pub drive_rate_start_mid_end_a_s: [f64; 3],
    pub gap_m: [f64; 2],
    /// Material gap endpoint minus the native mechanical proposal [m].
    pub displacement_residual_m: f64,
    pub kinetic_energy_j: [f64; 2],
    pub elastic_energy_j: [f64; 2],
    pub electrical_supply_j: f64,
    pub mechanical_supply_j: f64,
    pub field_work_on_body_j: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    /// Original resistor heat plus one core pinning-heat receipt.
    pub heat_j: f64,
    pub energy_residual_j: f64,
    pub circuit_trials: usize,
    pub mechanical_trials: usize,
    pub nonlinear: NewtonReport,
}

pub struct SharedHystereticGapMotionCandidate {
    model: SharedHystereticGapMotion,
    step: SharedHystereticGapMotionStep,
}
impl SharedHystereticGapMotionCandidate {
    pub const fn model(&self) -> &SharedHystereticGapMotion {
        &self.model
    }
    pub const fn step(&self) -> &SharedHystereticGapMotionStep {
        &self.step
    }
    pub fn into_parts(self) -> (SharedHystereticGapMotion, SharedHystereticGapMotionStep) {
        (self.model, self.step)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn phs(error: PhsError) -> HystereticCircuitError {
    fs_phs::circuit::CircuitError::from(error).into()
}

impl SharedHystereticGapMotion {
    /// Bind an explicitly prepared shared core/gap circuit and native body.
    /// Positive velocity opens the gap. Mass is the original regional body
    /// inertia plus the optional support's consistent modal mass.
    /// Establish this complete assembly before advancing its physical history.
    pub fn new(
        circuit: SharedHystereticCircuit,
        body: ResolvedBody,
        support: Option<AxialSupport>,
        velocity_m_s: f64,
    ) -> Result<Self> {
        if !velocity_m_s.is_finite()
            || circuit.core.gap().is_none()
            || support
                .as_ref()
                .is_some_and(AxialSupport::is_thermally_live)
        {
            return Err(invalid(
                "moving hysteretic core requires its original gap, finite velocity and frozen support",
            ));
        }
        let held = circuit.core.hold_candidate(&circuit.state)?;
        let gap = held
            .state
            .gap_m()
            .ok_or_else(|| invalid("shared gap state is missing"))?;
        if let Some(support) = &support {
            support.force_and_energy(gap)?;
        }
        let mass =
            body.inertia().mass() + support.as_ref().map_or(0.0, AxialSupport::modal_mass_kg);
        let momentum = mass * velocity_m_s;
        let kinetic = 0.5 * momentum * (momentum / mass);
        if !positive(mass)
            || !positive(mass.recip())
            || !momentum.is_finite()
            || velocity_m_s != 0.0 && momentum == 0.0
            || !kinetic.is_finite()
            || momentum != 0.0 && kinetic == 0.0
        {
            return Err(invalid(
                "source-derived body mass or shared-gap momentum is unrepresentable",
            ));
        }
        Ok(Self {
            circuit,
            body,
            support,
            momentum_kg_m_s: momentum,
        })
    }
    pub const fn circuit(&self) -> &SharedHystereticCircuit {
        &self.circuit
    }
    pub const fn body(&self) -> &ResolvedBody {
        &self.body
    }
    pub const fn support(&self) -> Option<&AxialSupport> {
        self.support.as_ref()
    }
    pub fn gap_m(&self) -> f64 {
        self.circuit.state.gap_m().expect("admitted shared gap")
    }
    pub fn generalized_mass_kg(&self) -> f64 {
        self.body.inertia().mass()
            + self
                .support
                .as_ref()
                .map_or(0.0, AxialSupport::modal_mass_kg)
    }
    pub const fn momentum_kg_m_s(&self) -> f64 {
        self.momentum_kg_m_s
    }
    pub fn velocity_m_s(&self) -> f64 {
        self.momentum_kg_m_s / self.generalized_mass_kg()
    }
    pub const fn time_s(&self) -> f64 {
        self.circuit.time_s
    }
    /// Refresh the same constant core fit, retaining body/gap and every history.
    /// This does not evolve thermal inventory or replace the fixed gap medium.
    pub fn with_core_source_query(&self, query: &QueryPoint) -> Result<Self> {
        Ok(Self {
            circuit: self.circuit.with_core_source_query(query)?,
            ..self.clone()
        })
    }

    /// Advance the common electrical/mechanical interval transactionally.
    /// Inputs and external force are held over this interval. Geometry limits,
    /// incompatible native clamps, unresolved flux turns or work defects refuse.
    pub fn step_candidate(
        &self,
        target_s: f64,
        resistors: Vec<MaterialResistor>,
        external_inputs: &[f64],
        external_force_n: f64,
        controls: SharedHystereticGapMotionSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<SharedHystereticGapMotionCandidate> {
        let cancel = RefCell::new(cancelled);
        let check = || {
            if (*cancel.borrow_mut())() {
                Err(HystereticCircuitError::Cancelled)
            } else {
                Ok(())
            }
        };
        check()?;
        let electrical = controls.circuit;
        electrical.newton.validate()?;
        let dt = target_s - self.time_s();
        let gap_scale = controls.difference_step_m / electrical.difference_step_wb;
        if !target_s.is_finite()
            || !positive(dt)
            || !external_force_n.is_finite()
            || electrical.maximum_circuit_trials < 4
            || electrical.maximum_newton_iterations == 0
            || [
                electrical.difference_step_wb,
                electrical.terminal_voltage_tolerance_v,
                electrical.energy_tolerance_j,
                electrical.circuit.energy_tolerance_j,
                electrical.circuit.power_tolerance_w,
                controls.difference_step_m,
                controls.displacement_tolerance_m,
                controls.energy_tolerance_j,
                gap_scale,
                gap_scale.recip(),
            ]
            .iter()
            .any(|&v| !positive(v))
            || external_inputs.len() != self.circuit.external_ports
            || external_inputs.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "positive bounded moving-gap scales, interval and complete finite inputs required",
            ));
        }
        let ports = NativeWindingPorts::new(&self.circuit, resistors, external_inputs)?;
        let gap0 = self.gap_m();
        let mass = self.generalized_mass_kg();
        let law = MechanicalStorage {
            mass,
            support: self.support.clone(),
        };
        let mechanics = DescriptorPortHamiltonian::new(
            2,
            2,
            2,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0; 4],
            vec![0.0, 0.0, 1.0, 1.0],
            Box::new(law.clone()),
        )
        .map_err(phs)?;
        let initial_mechanical = [gap0, self.momentum_kg_m_s];
        let trials = Cell::new(0usize);
        let mechanical_trials = Cell::new(0usize);
        let reserve = || {
            check()?;
            if trials.get() >= electrical.maximum_circuit_trials {
                return Err(HystereticCircuitError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            Ok(())
        };
        let count = self.circuit.branch_indices.len();
        let evaluate =
            |delta: &[f64]| -> Result<(SharedHystereticStep, StepRecord, StepRecord, Vec<f64>)> {
                reserve()?;
                if delta.len() != count + 1 {
                    return Err(invalid("moving-gap trial has the wrong dimension"));
                }
                let linkages = self
                    .circuit
                    .state
                    .flux_linkages_wb()
                    .iter()
                    .zip(delta)
                    .map(|(&old, &change)| {
                        let next = finite(old + change)?;
                        if change != 0.0 && next == old {
                            return Err(invalid(
                                "winding increment is below representable resolution",
                            ));
                        }
                        Ok(next)
                    })
                    .collect::<Result<Vec<_>>>()?;
                let displacement = finite(delta[count] * gap_scale)?;
                let gap1 = finite(gap0 + displacement)?;
                if displacement != 0.0 && gap1 == gap0 {
                    return Err(invalid("gap increment is below representable resolution"));
                }
                let magnetic = self.circuit.core.candidate_linkages_and_gap(
                    &self.circuit.state,
                    &linkages,
                    gap1,
                )?;
                let gap = magnetic
                    .gap
                    .as_ref()
                    .ok_or_else(|| invalid("moving material candidate lost the common gap"))?;
                let native = ports.step(external_inputs, &magnetic.mean_currents_a, dt)?;
                check()?;
                mechanical_trials.set(mechanical_trials.get() + 1);
                let mechanical = step_descriptor(
                    &mechanics,
                    &initial_mechanical,
                    &[external_force_n, gap.mean_force_n],
                    dt,
                )
                .map_err(phs)?;
                let mut residuals = linkages
                    .iter()
                    .zip(self.circuit.state.flux_linkages_wb())
                    .enumerate()
                    .map(|(i, (&after, &before))| {
                        finite((after - before) - dt * native.y[self.circuit.external_ports + i])
                    })
                    .collect::<Result<Vec<_>>>()?;
                residuals.push(finite((gap1 - mechanical.x[0]) / gap_scale)?);
                check()?;
                Ok((magnetic, native, mechanical, residuals))
            };
        let (solution, nonlinear) = solve_interface(
            count + 1,
            electrical,
            |delta| evaluate(delta).map(|(_, _, _, residuals)| residuals),
            &mut || (*cancel.borrow_mut())(),
        )?;
        check()?;
        let (magnetic, native, mechanical, residuals) = evaluate(&solution)?;
        let gap1 = magnetic
            .state
            .gap_m()
            .ok_or_else(|| invalid("accepted magnetic gap is missing"))?;
        let gap = magnetic
            .gap
            .as_ref()
            .ok_or_else(|| invalid("accepted magnetic gap work is missing"))?;
        let displacement_residual_m = finite(gap1 - mechanical.x[0])?;
        let terminal_voltages_v = native.y[self.circuit.external_ports..].to_vec();
        let terminal_voltage_residuals_v = residuals[..count]
            .iter()
            .map(|&residual| finite(residual / dt))
            .collect::<Result<Vec<_>>>()?;
        if displacement_residual_m.abs() > controls.displacement_tolerance_m
            || terminal_voltage_residuals_v
                .iter()
                .any(|v| v.abs() > electrical.terminal_voltage_tolerance_v)
        {
            return Err(HystereticCircuitError::Balance(
                "native winding voltages or mechanical displacement disagree with the magnetic chord",
            ));
        }
        let next_fluxes: Vec<_> = ports.dae.flux_index.iter().map(|&i| native.x[i]).collect();
        let next_charges: Vec<_> = ports
            .dae
            .charge_index
            .iter()
            .map(|&i| native.x[i])
            .collect();
        let mid_fluxes: Vec<_> = self
            .circuit
            .circuit
            .fluxes_wb
            .iter()
            .zip(&next_fluxes)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let mid_charges: Vec<_> = self
            .circuit
            .circuit
            .charges_c
            .iter()
            .zip(&next_charges)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let mid_linkages: Vec<_> = self
            .circuit
            .state
            .flux_linkages_wb()
            .iter()
            .zip(magnetic.state.flux_linkages_wb())
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let middle = self.circuit.core.candidate_linkages_and_gap(
            &self.circuit.state,
            &mid_linkages,
            f64::midpoint(gap0, gap1),
        )?;
        let reluctance_rate = self
            .circuit
            .core
            .gap()
            .expect("admitted gap")
            .reluctance_per_m();
        let mut rates = [0.0; 3];
        let mut rate_tolerance = 0.0_f64;
        for (index, (fluxes, charges, currents, flux, momentum)) in [
            (
                &self.circuit.circuit.fluxes_wb,
                &self.circuit.circuit.charges_c,
                self.circuit.state.currents_a(),
                self.circuit.reference_flux_linkage_wb,
                self.momentum_kg_m_s,
            ),
            (
                &mid_fluxes,
                &mid_charges,
                middle.state.currents_a(),
                middle.core.flux_linkage_wb,
                f64::midpoint(self.momentum_kg_m_s, mechanical.x[1]),
            ),
            (
                &next_fluxes,
                &next_charges,
                magnetic.state.currents_a(),
                magnetic.core.flux_linkage_wb,
                mechanical.x[1],
            ),
        ]
        .into_iter()
        .enumerate()
        {
            reserve()?;
            let voltages = instantaneous_voltages(
                &ports.dae,
                &ports.graph,
                &self.circuit.branch_indices,
                fluxes,
                charges,
                external_inputs,
                currents,
            )?;
            let mut tolerance =
                finite((reluctance_rate * flux).abs() * controls.displacement_tolerance_m / dt)?;
            for (&voltage, winding) in voltages.iter().zip(self.circuit.core.windings()) {
                rates[index] =
                    finite(rates[index] + winding.turns_ratio * (voltage / winding.leakage_h))?;
                tolerance = finite(
                    tolerance
                        + winding.turns_ratio.abs() / winding.leakage_h
                            * electrical.terminal_voltage_tolerance_v,
                )?;
            }
            rates[index] =
                finite((-reluctance_rate * flux).mul_add(momentum / mass, rates[index]))?;
            rate_tolerance = rate_tolerance.max(tolerance);
        }
        let flux_change =
            finite(magnetic.core.flux_linkage_wb - self.circuit.reference_flux_linkage_wb)?;
        if rates.iter().any(|&rate| {
            if flux_change == 0.0 {
                rate.abs() > rate_tolerance
            } else {
                rate * flux_change.signum() < -rate_tolerance
            }
        }) {
            return Err(HystereticCircuitError::SharedMagneticReversal {
                reference_linkage_change_wb: flux_change,
                drive_rate_start_mid_end_a_s: rates,
            });
        }
        let kinetic_energy_j = [
            law.kinetic(self.momentum_kg_m_s),
            law.kinetic(mechanical.x[1]),
        ];
        let support_before = law.support(gap0)?;
        let support_after = law.support(gap1)?;
        let elastic_energy_j = [support_before.1, support_after.1];
        let kinetic_change = finite(
            0.5 * (mechanical.x[1] - self.momentum_kg_m_s)
                * ((mechanical.x[1] + self.momentum_kg_m_s) / mass),
        )?;
        let elastic_change =
            finite(-f64::midpoint(support_before.0, support_after.0) * (gap1 - gap0))?;
        let mechanical_change = finite(kinetic_change + elastic_change)?;
        let electrical_supply_j = finite(
            dt * external_inputs
                .iter()
                .zip(&native.y)
                .map(|(&u, &y)| u * y)
                .sum::<f64>(),
        )?;
        let mechanical_supply_j = finite(external_force_n * (gap1 - gap0))?;
        let field_work_on_body_j = gap.mechanical_work_j;
        let native_supplied = native.supplied;
        let stored_change_j =
            finite(native.delta_h + magnetic.stored_change_j + mechanical_change)?;
        let before_electrical = ports.stored_before_j;
        let (circuit, circuit_step) =
            ports.finish(&self.circuit, native, dt, electrical.circuit, &mut || {
                (*cancel.borrow_mut())()
            })?;
        let heat_j = finite(circuit_step.heat_j + magnetic.core.core_heat_j)?;
        let energy_residual_j =
            finite(stored_change_j + heat_j - electrical_supply_j - mechanical_supply_j)?;
        if energy_residual_j.abs() > controls.energy_tolerance_j
            || (electrical_supply_j - native_supplied - magnetic.magnetic_work_j).abs()
                > electrical.energy_tolerance_j
            || (mechanical_change - mechanical_supply_j - field_work_on_body_j).abs()
                > controls.energy_tolerance_j
            || (mechanical.supplied - mechanical_supply_j - field_work_on_body_j).abs()
                > controls.energy_tolerance_j
            || mechanical.balance_residual().abs() > controls.energy_tolerance_j
            || mechanical.dissipated != 0.0
        {
            return Err(HystereticCircuitError::Balance(
                "shared hysteretic electrical, mechanical, leakage/core/gap storage and physical heat do not close",
            ));
        }
        let model = Self {
            circuit: SharedHystereticCircuit {
                circuit,
                core: self.circuit.core.clone(),
                state: magnetic.state.clone(),
                branch_indices: self.circuit.branch_indices.clone(),
                external_ports: self.circuit.external_ports,
                reference_flux_linkage_wb: magnetic.core.flux_linkage_wb,
                time_s: target_s,
            },
            body: self.body.clone(),
            support: self.support.clone(),
            momentum_kg_m_s: mechanical.x[1],
        };
        let step = SharedHystereticGapMotionStep {
            stored_before_j: finite(
                before_electrical
                    + self.circuit.state.stored_energy_j()
                    + kinetic_energy_j[0]
                    + elastic_energy_j[0],
            )?,
            stored_after_j: finite(
                circuit_step.stored_after_j
                    + magnetic.stored_energy_j
                    + kinetic_energy_j[1]
                    + elastic_energy_j[1],
            )?,
            circuit: circuit_step,
            magnetic,
            mechanical,
            terminal_voltages_v,
            terminal_voltage_residuals_v,
            drive_rate_start_mid_end_a_s: rates,
            gap_m: [gap0, gap1],
            displacement_residual_m,
            kinetic_energy_j,
            elastic_energy_j,
            electrical_supply_j,
            mechanical_supply_j,
            field_work_on_body_j,
            stored_change_j,
            heat_j,
            energy_residual_j,
            circuit_trials: trials.get(),
            mechanical_trials: mechanical_trials.get(),
            nonlinear,
        };
        check()?;
        Ok(SharedHystereticGapMotionCandidate { model, step })
    }
}

#[derive(Clone)]
struct MechanicalStorage {
    mass: f64,
    support: Option<AxialSupport>,
}
impl MechanicalStorage {
    fn kinetic(&self, p: f64) -> f64 {
        0.5 * p * (p / self.mass)
    }
    fn support(&self, gap: f64) -> Result<(f64, f64)> {
        match &self.support {
            Some(support) => Ok(support.force_and_energy(gap)?),
            None => Ok((0.0, 0.0)),
        }
    }
}
impl Storage for MechanicalStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        self.support(x[0])
            .map_or(f64::NAN, |(_, energy)| energy + self.kinetic(x[1]))
    }
    fn gradient(&self, x: &[f64], output: &mut [f64]) {
        match self.support(x[0]) {
            Ok((force, _)) => {
                output[0] = -force;
                output[1] = x[1] / self.mass;
            }
            Err(_) => output.fill(f64::NAN),
        }
    }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        let (Ok((fa, _)), Ok((fb, _)), Ok((fm, _))) =
            (self.support(a[0]), self.support(b[0]), self.support(mid[0]))
        else {
            return Some(f64::NAN);
        };
        Some(
            0.5 * ((a[1] - mid[1]) + (b[1] - mid[1])) * ((b[1] - a[1]) / self.mass)
                - 0.5 * ((fa - fm) + (fb - fm)) * (b[0] - a[0]),
        )
    }
}
