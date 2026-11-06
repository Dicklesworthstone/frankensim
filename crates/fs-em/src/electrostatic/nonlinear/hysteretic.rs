//! Retained planar polarization on the original nonlinear Gauss/electrode solve.
//!
//! The frozen-history voltage potential is Phi = E.D - W - Q, whose gradient
//! is D and Hessian is the consistent dD/dE. Physical W and switching Q remain
//! separate from this Newton objective. Every Newton trial uses the SAME old
//! polarization; only a completely accepted field becomes the next history.
//!
//! This is a fixed-geometry, isotropic, quasistatic sequence. Each finite vector
//! increment uses endpoint electrical work, with an explicitly separate
//! nonnegative numerical remainder. Resolve reversals/rotations and refine the
//! path; neither that remainder nor a reversible Maxwell-force formula supplies
//! physical heat or switching strain. No rate, leakage, temperature evolution,
//! breakdown, continuum error bound or measured ferroelectric accuracy is implied.

use super::*;
use fs_material::field_storage::dielectric::hysteresis::vector::{
    ResolvedVectorDielectricHysteresis, VectorDielectricHysteresisControls,
    VectorDielectricHysteresisError, VectorDielectricHysteresisState,
    VectorDielectricHysteresisStep,
};

fn history_error(cell: usize, error: VectorDielectricHysteresisError) -> Error {
    match error {
        VectorDielectricHysteresisError::Cancelled => ElectrostaticError::Cancelled.into(),
        error => Error::Hysteresis { cell, error },
    }
}

/// One original material selection and immutable accepted polarization state.
/// This is an incremental constitutive boundary, not a reversible dielectric.
#[derive(Clone, Debug, PartialEq)]
pub struct FrozenHystereticDielectric {
    source: ResolvedVectorDielectricHysteresis,
    history: VectorDielectricHysteresisState,
    controls: VectorDielectricHysteresisControls,
}
impl FrozenHystereticDielectric {
    pub fn new(
        source: ResolvedVectorDielectricHysteresis,
        controls: VectorDielectricHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let history = source.law().initial_state();
        Self::from_history(source, history, controls, cancelled)
    }

    /// Explicit prepared memory must belong to this fit and lie in the mesh
    /// plane. A standalone frozen solve leaves ownership of prior global field
    /// consistency with its caller; PlanarHystereticElectrostatics owns both.
    pub fn from_history(
        source: ResolvedVectorDielectricHysteresis,
        history: VectorDielectricHysteresisState,
        controls: VectorDielectricHysteresisControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let result = Self {
            source,
            history,
            controls,
        };
        result.admit_planar(0, cancelled)?;
        Ok(result)
    }
    pub fn source(&self) -> &ResolvedVectorDielectricHysteresis {
        &self.source
    }
    pub fn history(&self) -> &VectorDielectricHysteresisState {
        &self.history
    }

    pub(super) fn admit_planar(
        &self,
        cell: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), Error> {
        poll(cancelled)?;
        if self.history.field_v_m()[2] != 0.0
            || self.history.displacement_c_m2()[2] != 0.0
            || self
                .history
                .polarizations_c_m2()
                .iter()
                .any(|p| p[2] != 0.0)
        {
            return Err(invalid(
                "planar history cannot hide out-of-plane polarization or displacement",
            ));
        }
        self.source
            .law()
            .candidate_field(
                &self.history,
                self.history.field_v_m(),
                self.controls,
                cancelled,
            )
            .map_err(|error| history_error(cell, error))?;
        Ok(())
    }
    pub(super) fn evaluate(
        &self,
        cell: usize,
        electric: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<VectorDielectricHysteresisStep, Error> {
        let magnitude = finite(electric[0].hypot(electric[1]))?;
        let maximum = self.source.law().field_limit_v_m();
        if magnitude > maximum {
            return Err(Error::FieldRange {
                cell,
                observed_v_m: magnitude,
                maximum_v_m: maximum,
            });
        }
        self.source
            .law()
            .candidate_field(
                &self.history,
                [electric[0], electric[1], 0.0],
                self.controls,
                cancelled,
            )
            .map_err(|error| history_error(cell, error))
    }
    pub(super) fn coenergy(
        &self,
        cell: usize,
        step: &VectorDielectricHysteresisStep,
    ) -> Result<f64, Error> {
        step.incremental_coenergy_j_m3()
            .map_err(|error| history_error(cell, error))
    }
    pub(super) fn coenergy_change(
        &self,
        cell: usize,
        before: &VectorDielectricHysteresisStep,
        after: &VectorDielectricHysteresisStep,
    ) -> Result<f64, Error> {
        self.source
            .law()
            .incremental_coenergy_difference_j_m3(&self.history, &before.state, &after.state)
            .map_err(|error| history_error(cell, error))
    }
}

/// One physical cell's discrete switching account, already volume scaled.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricSwitchingCell {
    pub triangle: usize,
    pub volume_m3: f64,
    pub stored_energy_change_j: f64,
    pub switching_heat_j: f64,
    pub endpoint_work_j: f64,
    /// Endpoint discretization loss, never material heating.
    pub numerical_dissipation_j: f64,
    pub energy_residual_j: f64,
}

/// Independent field and electrode/volume-charge work accounts of one increment.
pub struct HystereticElectrostaticStep {
    pub cells: Vec<DielectricSwitchingCell>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub stored_change_j: f64,
    pub switching_heat_j: f64,
    pub numerical_dissipation_j: f64,
    pub endpoint_field_work_j: f64,
    /// sum(V_new delta Q_electrode) + integral(V_new delta rho_free).
    /// This is endpoint quadrature, not an exact unresolved-path work integral.
    pub endpoint_source_work_j: f64,
    pub field_source_work_residual_j: f64,
    pub energy_residual_j: f64,
    pub charge_balance_residual_c: f64,
    pub maximum_free_charge_residual_c: f64,
    pub maximum_floating_charge_residual_c: f64,
    pub nonlinear: NonlinearElectrostaticReport,
}

/// Accepted whole-mesh memory. Geometry, gradients and electrode layout are
/// prepared once and shared; only physical field/memory arrays advance.
#[derive(Clone)]
pub struct PlanarHystereticElectrostatics {
    problem: PlanarNonlinearElectrostatics,
    layout: Arc<electrodes::Layout>,
    relative_voltage_v: Vec<f64>,
    reference_voltage_v: f64,
    terminal_voltage_v: Vec<f64>,
    terminal_charge_c: Vec<f64>,
    free_charge_density_c_m3: Vec<f64>,
    stored_energy_j: f64,
}

pub struct HystereticElectrostaticCandidate {
    field: PlanarHystereticElectrostatics,
    step: HystereticElectrostaticStep,
}
impl HystereticElectrostaticCandidate {
    pub fn field(&self) -> &PlanarHystereticElectrostatics {
        &self.field
    }
    pub fn step(&self) -> &HystereticElectrostaticStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarHystereticElectrostatics, HystereticElectrostaticStep) {
        (self.field, self.step)
    }
}

impl PlanarHystereticElectrostatics {
    /// Start at the compatible virgin zero-field/zero-charge state. Source
    /// roles are explicit in every cell; zero switching branches give ordinary
    /// linear dielectric cells. Prepare poled specimens through accepted steps.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedVectorDielectricHysteresis>,
        terminals: &[&[usize]],
        local: VectorDielectricHysteresisControls,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if terminals.is_empty() {
            return Err(invalid("original electrode traces are required"));
        }
        let mut laws = Vec::with_capacity(materials.len());
        for source in materials {
            poll(cancelled)?;
            laws.push(DielectricMaterial::Hysteretic(
                FrozenHystereticDielectric::new(source, local, cancelled)?,
            ));
        }
        let problem = PlanarNonlinearElectrostatics::new(mesh, laws, limits, cancelled)?;
        let relative_voltage_v = vec![0.0; problem.mesh.vertices().len()];
        let prepared = problem.prepare_electrodes(terminals, cancelled)?;
        prepared.coordinates(
            &vec![ElectrodeControl::Voltage(0.0); terminals.len()],
            &relative_voltage_v,
            cancelled,
        )?;
        let layout = prepared.layout.clone();
        let cells = problem.elements.len();
        Ok(Self {
            problem,
            layout,
            relative_voltage_v,
            reference_voltage_v: 0.0,
            terminal_voltage_v: vec![0.0; terminals.len()],
            terminal_charge_c: vec![0.0; terminals.len()],
            free_charge_density_c_m3: vec![0.0; cells],
            stored_energy_j: 0.0,
        })
    }
    pub fn problem(&self) -> &PlanarNonlinearElectrostatics {
        &self.problem
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        &self.relative_voltage_v
    }
    pub fn reference_voltage_v(&self) -> f64 {
        self.reference_voltage_v
    }
    pub fn terminal_voltage_v(&self) -> &[f64] {
        &self.terminal_voltage_v
    }
    pub fn terminal_charge_c(&self) -> &[f64] {
        &self.terminal_charge_c
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        &self.free_charge_density_c_m3
    }
    pub fn cell_history(&self, triangle: usize) -> Option<&VectorDielectricHysteresisState> {
        match self.problem.materials.get(triangle) {
            Some(DielectricMaterial::Hysteretic(law)) => Some(law.history()),
            _ => None,
        }
    }

    /// Solve mixed voltage/total-charge electrodes with the original matrix-free
    /// Newton/checked-CG owner. Refusal leaves every accepted cell unchanged.
    /// All source fits remain fixed; this increment does not advance heat.
    pub fn solve_candidate(
        &self,
        free_charge_density_c_m3: &[f64],
        boundary: &[ElectrodeControl],
        controls: NonlinearElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HystereticElectrostaticCandidate, Error> {
        poll(cancelled)?;
        let prepared = NonlinearElectrodes {
            problem: &self.problem,
            layout: self.layout.clone(),
        };
        let reference = boundary
            .iter()
            .position(|c| matches!(c, ElectrodeControl::Voltage(_)))
            .filter(|&index| index < self.layout.terminals.len())
            .ok_or_else(|| invalid("a voltage-controlled reference is required"))?;
        // Retain small fields even if an authored common voltage bias is huge.
        let shift = self.relative_voltage_v[self.layout.terminals[reference][0]];
        let initial: Vec<_> = self
            .relative_voltage_v
            .iter()
            .map(|v| finite(*v - shift))
            .collect::<Result<_, _>>()?;
        let solved = prepared.solve(
            free_charge_density_c_m3,
            boundary,
            &initial,
            controls,
            cancelled,
        )?;
        let mut cells = Vec::with_capacity(self.problem.elements.len());
        let mut next = self.problem.clone();
        let mut totals = [0.0; 4];
        let mut correction = [0.0; 4];
        for (cell, element) in self.problem.elements.iter().enumerate() {
            poll(cancelled)?;
            let response = solved.evaluation.switching[cell].as_ref().ok_or_else(|| {
                invalid("retained dielectric cell lost its native switching response")
            })?;
            let terms = [
                response.energy_change_j_m3,
                response.irreversible_heat_j_m3,
                response.numerical_dissipation_j_m3,
                response.endpoint_work_j_m3,
            ];
            let mut values = [0.0; 4];
            for i in 0..4 {
                values[i] = finite(element.volume_m3 * terms[i])?;
                if terms[i] != 0.0 && values[i] == 0.0 {
                    return Err(invalid("dielectric cell transfer underflowed"));
                }
                add(&mut totals[i], &mut correction[i], values[i])?;
            }
            let residual = finite(values[3] - values[0] - values[1] - values[2])?;
            if residual.abs() > controls.energy_tolerance_j {
                return Err(invalid(
                    "individual dielectric switching account did not close",
                ));
            }
            cells.push(DielectricSwitchingCell {
                triangle: cell,
                volume_m3: element.volume_m3,
                stored_energy_change_j: values[0],
                switching_heat_j: values[1],
                numerical_dissipation_j: values[2],
                endpoint_work_j: values[3],
                energy_residual_j: residual,
            });
            let DielectricMaterial::Hysteretic(law) = &mut next.materials[cell] else {
                return Err(invalid("retained dielectric source correspondence changed"));
            };
            law.history = response.state.clone();
        }
        for i in 0..4 {
            totals[i] = finite(totals[i] + correction[i])?;
        }
        // Independent external pairing; no constitutive work is substituted.
        let (mut supplied, mut cs) = (0.0, 0.0);
        for (t, trace) in self.layout.terminals.iter().enumerate() {
            poll(cancelled)?;
            add(
                &mut supplied,
                &mut cs,
                solved.relative_voltage_v[trace[0]]
                    * (solved.terminal_charge_c[t] - self.terminal_charge_c[t]),
            )?;
        }
        for (cell, element) in self.problem.elements.iter().enumerate() {
            poll(cancelled)?;
            let charge = finite(
                element.volume_m3
                    * (free_charge_density_c_m3[cell] - self.free_charge_density_c_m3[cell]),
            )?;
            let mean = element
                .vertices
                .iter()
                .map(|&v| solved.relative_voltage_v[v] / 3.0)
                .sum::<f64>();
            add(&mut supplied, &mut cs, charge * mean)?;
        }
        supplied = finite(supplied + cs)?;
        let work_residual = finite(supplied - totals[3])?;
        let energy_residual = finite(supplied - totals[0] - totals[1] - totals[2])?;
        if work_residual.abs() > controls.energy_tolerance_j
            || energy_residual.abs() > controls.energy_tolerance_j
            || (solved.switching_heat_j() - totals[1]).abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "electrode/volume work, retained storage, switching heat and numerical loss do not close",
            ));
        }
        let step = HystereticElectrostaticStep {
            cells,
            stored_before_j: self.stored_energy_j,
            stored_after_j: solved.energy_j(),
            stored_change_j: totals[0],
            switching_heat_j: totals[1],
            numerical_dissipation_j: totals[2],
            endpoint_field_work_j: totals[3],
            endpoint_source_work_j: supplied,
            field_source_work_residual_j: work_residual,
            energy_residual_j: energy_residual,
            charge_balance_residual_c: solved.charge_balance_residual_c,
            maximum_free_charge_residual_c: solved.maximum_free_charge_residual_c,
            maximum_floating_charge_residual_c: solved.maximum_floating_charge_residual_c,
            nonlinear: solved.report.clone(),
        };
        let field = Self {
            problem: next,
            layout: self.layout.clone(),
            relative_voltage_v: solved.relative_voltage_v,
            reference_voltage_v: solved.reference_voltage_v,
            terminal_voltage_v: solved.terminal_voltage_v,
            terminal_charge_c: solved.terminal_charge_c,
            free_charge_density_c_m3: free_charge_density_c_m3.to_vec(),
            stored_energy_j: step.stored_after_j,
        };
        poll(cancelled)?;
        Ok(HystereticElectrostaticCandidate { field, step })
    }
}
