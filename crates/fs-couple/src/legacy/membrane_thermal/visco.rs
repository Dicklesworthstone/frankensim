//! Finite-strain Maxwell heat coupled to the original membrane enthalpy clock.
//!
//! The native membrane owns deformation, momentum and every branch history.
//! The existing thermal composition owns interval selection, events, transport,
//! original source/inventory admission and conservative heat handoff. Only the
//! original per-triangle relaxation increments become thermal sources; material
//! endpoint defects, pressure quadrature error and numerical kinetic damping
//! remain separate. A failed interval publishes neither mechanical nor thermal
//! state, including when the user's observer refuses after heat was proposed.
//!
//! The baseline entry point holds mechanical coefficients and natural references.
//! The explicit fully-solid temperature interval admits that approximation;
//! `temperature` separately supplies actual-temperature kinetic feedback through
//! the original shift owner. Native prescribed-reference events remain separate,
//! with their signed external work.

/// Actual-temperature kinetic feedback using the original relaxation spectrum.
pub mod temperature;

use std::borrow::Cow;
use std::collections::BTreeSet;

use fs_exec::CancelGate;
use fs_material::sheet::MembraneLaw;
use fs_solid::membrane::visco::FiniteViscoMembraneEnergy;
use fs_solid::membrane::visco::dynamics::{
    FiniteViscoMembraneMotion, FiniteViscoMotionControls, InertialEnergy,
};
use fs_solid::membrane::{MembraneError, MembranePatch, MembraneSolution, MembraneSolveSettings};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation};

use super::{
    ElementThermalCell, MechanicalInterval, MembraneThermalConfig, MembraneThermalError,
    MembraneThermalReport, Result, Sum, ThermalMotion, advance_native, check, number,
};

/// Independent native mechanical controls on the shared physical timeline.
#[derive(Clone, Copy, Debug)]
pub struct ViscoMembraneThermalControls {
    pub solve: MembraneSolveSettings,
    pub motion: FiniteViscoMotionControls,
    /// Upper bound before the thermal owner applies its own CFL/event limits.
    pub maximum_step_s: f64,
    /// Held signed pressure jump on the current surface, in Pa. Nonzero pressure
    /// retains the native complete prescribed-perimeter admission.
    pub pressure_pa: f64,
}

// One work identity for both frozen and temperature-shifted native motion.
fn mechanical_interval<E>(
    energy: &FiniteViscoMembraneEnergy,
    inertia: &InertialEnergy,
    equilibrium: &MembraneSolution,
    time_s: f64,
    dt_s: f64,
    prescribed: &[(usize, f64)],
) -> Result<MechanicalInterval, E> {
    // Held loads use displacement work; pressure uses native path work. Its
    // endpoint quadrature difference remains numerical error, never heat.
    let mut work = Sum::default();
    work.add(energy.dead_load_endpoint_work_j)?;
    work.add(energy.support_endpoint_work_j)?;
    work.add(energy.pressure_path_work_j)?;
    let pressure_defect = number(energy.pressure_endpoint_work_j - energy.pressure_path_work_j)?;
    let mut defect = Sum::default();
    let mut absolute_defect = Sum::default();
    for term in [
        -energy.endpoint_work_defect_j,
        -inertia.numerical_loss_j,
        -inertia.balance_error_j,
        -energy.equilibrium_work_residual_j,
        pressure_defect,
    ] {
        defect.add(term)?;
        absolute_defect.add(term.abs())?;
    }
    let mut support = [Sum::default(); 3];
    for &(dof, _) in prescribed {
        support[dof % 3].add(number(equilibrium.residual_n[dof / 3][dof % 3] * dt_s)?)?;
    }
    Ok(MechanicalInterval {
        time_s,
        dissipated_j: energy.dissipated_j,
        mechanical_work_j: work.value()?,
        energy_defect_j: defect.value()?,
        absolute_energy_defect_j: absolute_defect.value()?,
        numerical_kinetic_loss_j: inertia.numerical_loss_j,
        support_impulse_n_s: [
            support[0].value()?,
            support[1].value()?,
            support[2].value()?,
        ],
    })
}

/// Whole-call proposal with retained native memory, original inventory and observer.
pub struct ViscoMembraneThermalCandidate<'patch, 'thermal, O> {
    motion: FiniteViscoMembraneMotion<'patch>,
    thermal: DepositionSimulation<'thermal>,
    observer: O,
    report: MembraneThermalReport,
}

impl<'patch, 'thermal, O> ViscoMembraneThermalCandidate<'patch, 'thermal, O> {
    pub fn motion(&self) -> &FiniteViscoMembraneMotion<'patch> {
        &self.motion
    }
    pub fn thermal(&self) -> &DepositionSimulation<'thermal> {
        &self.thermal
    }
    pub fn observer(&self) -> &O {
        &self.observer
    }
    pub fn report(&self) -> &MembraneThermalReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        FiniteViscoMembraneMotion<'patch>,
        DepositionSimulation<'thermal>,
        O,
        MembraneThermalReport,
    ) {
        (self.motion, self.thermal, self.observer, self.report)
    }
}

#[derive(Clone)]
struct Adapter<'call, 'patch, 'boundary> {
    motion: Cow<'call, FiniteViscoMembraneMotion<'patch>>,
    prescribed: &'boundary [(usize, f64)],
    controls: ViscoMembraneThermalControls,
    kinetic_j: f64,
    element_heat_j: Vec<f64>,
    source_temperature_k: Vec<Option<f64>>,
}

impl ThermalMotion for Adapter<'_, '_, '_> {
    fn patch(&self) -> &MembranePatch {
        self.motion.membrane().patch()
    }
    fn time_s(&self) -> f64 {
        self.motion.membrane().time_s()
    }
    fn strain_energy_j(&self) -> f64 {
        self.motion.membrane().stored_energy_j()
    }
    fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_j
    }
    fn max_step_s(&self) -> f64 {
        self.controls.maximum_step_s
    }
    fn source_temperature_k(&self, element: usize) -> Option<f64> {
        self.source_temperature_k[element]
    }
    fn element_dissipated_j(&self) -> &[f64] {
        &self.element_heat_j
    }

    fn advance<E>(
        &mut self,
        dt_s: f64,
        dead_load_n: &[[f64; 3]],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<MechanicalInterval, E> {
        let step = self
            .motion
            .advance_with_pressure(
                dt_s,
                dead_load_n,
                self.controls.pressure_pa,
                self.prescribed,
                self.controls.solve,
                self.controls.motion,
                cancelled,
            )
            .map_err(MembraneThermalError::Viscoelastic)?;
        let interval = mechanical_interval(
            &step.material_energy,
            &step.inertia,
            &step.equilibrium,
            step.motion.membrane().time_s(),
            dt_s,
            self.prescribed,
        )?;
        self.kinetic_j = step.inertia.kinetic_after_j;
        self.element_heat_j = step.element_dissipated_j;
        self.motion = Cow::Owned(step.motion);
        Ok(interval)
    }
}

/// Advance original finite-strain viscoelastic motion and enthalpy together.
///
/// The thermal owner chooses each actual interval, including event/source cuts;
/// the material return map and momentum equation use exactly that interval.
/// Native per-triangle physical heat is retained before cumulative rounding,
/// then handed to the existing caloric receiver as its average interval power.
/// Receiver source, mass, reference volume, solid phase and temperature-range
/// admission are identical to the existing Kelvin--Voigt thermal composition.
///
/// Prescribed coordinates must equal their accepted displacements and remain
/// stationary throughout this call. Dead loads and pressure are held. This API
/// does not reinterpret a final displacement as a loading ramp across thermal
/// event cuts. A separate native reference event can change prescribed natural
/// strain, with its reference work reported before this physical interval.
///
/// All inputs remain usable after cancellation or a late observer/work refusal.
/// Backward-Euler kinetic damping and signed material/pressure integration errors
/// are bounded independently of physical heat, and their absolute contributions
/// cannot cancel before the whole-call numerical-error gate.
#[allow(clippy::too_many_arguments)]
pub fn advance_candidate<'patch, 'thermal, O: DepositionObserver>(
    motion: &FiniteViscoMembraneMotion<'patch>,
    thermal: &DepositionSimulation<'thermal>,
    observer: &O,
    target_s: f64,
    dead_load_n: &[[f64; 3]],
    prescribed: &[(usize, f64)],
    cells: &[ElementThermalCell<'_>],
    mechanics: ViscoMembraneThermalControls,
    config: MembraneThermalConfig,
    gate: &CancelGate,
) -> Result<ViscoMembraneThermalCandidate<'patch, 'thermal, O>, O::Error> {
    let (kinetic_j, source_temperature_k) = admit_motion(motion, prescribed, mechanics, gate)?;
    let adapter = Adapter {
        motion: Cow::Borrowed(motion),
        prescribed,
        controls: mechanics,
        kinetic_j,
        element_heat_j: vec![0.0; motion.membrane().patch().sheets().len()],
        source_temperature_k,
    };
    let candidate = advance_native(
        &adapter,
        thermal,
        observer,
        target_s,
        dead_load_n,
        cells,
        config,
        gate,
    )?;
    Ok(ViscoMembraneThermalCandidate {
        motion: candidate.motion.motion.into_owned(),
        thermal: candidate.thermal,
        observer: candidate.observer,
        report: candidate.report,
    })
}

fn admit_motion<E>(
    motion: &FiniteViscoMembraneMotion<'_>,
    prescribed: &[(usize, f64)],
    mechanics: ViscoMembraneThermalControls,
    gate: &CancelGate,
) -> Result<(f64, Vec<Option<f64>>), E> {
    check(gate)?;
    if !mechanics.maximum_step_s.is_finite()
        || mechanics.maximum_step_s <= 0.0
        || !mechanics.pressure_pa.is_finite()
    {
        return Err(MembraneThermalError::Invalid(
            "positive finite membrane timestep ceiling and finite held pressure required",
        ));
    }
    let patch = motion.membrane().patch();
    let u = motion.membrane().displacement_m();
    let mut held = BTreeSet::new();
    for &(dof, value) in prescribed {
        check(gate)?;
        if dof >= 3 * u.len()
            || !value.is_finite()
            || !held.insert(dof)
            || value != u[dof / 3][dof % 3]
            || motion.velocity_mps()[dof / 3][dof % 3] != 0.0
        {
            return Err(MembraneThermalError::Invalid(
                "unique original coordinates held at accepted displacement with zero support velocity required",
            ));
        }
    }
    let momentum = motion
        .nodal_momentum_kg_m_s(|| gate.is_requested())
        .map_err(MembraneThermalError::Viscoelastic)?;
    let mut kinetic = Sum::default();
    for (v, p) in motion.velocity_mps().iter().zip(&momentum) {
        check(gate)?;
        for axis in 0..3 {
            kinetic.add(0.5 * v[axis] * p[axis])?;
        }
    }
    let kinetic_j = kinetic.value()?;
    if kinetic_j < 0.0 {
        return Err(MembraneThermalError::Invalid(
            "negative original membrane kinetic energy",
        ));
    }
    let mut source_temperature_k = Vec::with_capacity(patch.sheets().len());
    for (element, sheet) in patch.sheets().iter().enumerate() {
        check(gate)?;
        let law = sheet.membrane_law().map_err(|source| {
            MembraneThermalError::Membrane(MembraneError::Material { element, source })
        })?;
        let MembraneLaw::Hyperelastic(law) = law else {
            return Err(MembraneThermalError::Invalid(
                "original hyperelastic membrane source required",
            ));
        };
        source_temperature_k.push(law.model().context().axes().get("T").copied());
    }
    Ok((kinetic_j, source_temperature_k))
}
