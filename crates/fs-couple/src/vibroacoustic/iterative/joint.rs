//! Simultaneous original structural and outgoing acoustic equations.
//!
//! A p - B(i omega rho T v) = 0,  Z_s v + T^T W p = F,
//! Z_s = C + i(K/omega - omega M), using peak exp(-i omega t).
//! Neither the bare mechanical pencil nor an acoustic impedance is inverted.
//! Thus actual radiation can regularize a bare structural resonance. A truly
//! unresolved loaded mode still refuses through the original equation checks.

use super::{
    ComplexSum, StructuralRadiationError, StructuralRadiationPencil, balance, checkpoint, filled,
    finite, invalid, mean_work, positive, tolerance,
};
use fs_bem::{
    helmholtz::{
        Formulation, IterativeError, IterativeOptions, IterativeRadiation, Medium,
        RadiationSolution, ResourceReceipt, bc::Admittance,
    },
    panel3d::SpherePanels,
};
use fs_exec::{CancelGate, Cx};
use fs_math::{c64::C64, dd::Dd, det};
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedSupportedHarmonicResponse,
};
use fs_solver::{
    krylov::{GmresState, SolveReport},
    op::LinearOp,
};
use std::{
    cell::{Cell, RefCell},
    time::Instant,
};

/// Separate original mechanical, trace, solver-payload and physical-work limits.
#[derive(Clone, Copy, Debug)]
pub struct JointStructuralRadiationControls {
    /// Original mechanical admission; its residual ceiling also gates the final
    /// force-block defect in this owner's physical diagonal metric. That joint
    /// defect is a block norm, not the native standalone coefficientwise error.
    pub mechanical: HarmonicResponseControls,
    pub maximum_trace_entries: usize,
    /// Added real GMRES basis/Hessenberg, candidate and native action scratch.
    /// Original mechanical/FMM payload, borrowed source data and allocator
    /// metadata are additional and retain their own original admission.
    pub maximum_payload_bytes: usize,
    pub maximum_relative_power_defect: f64,
}

#[derive(Debug)]
pub enum JointStructuralRadiationError {
    Source(StructuralRadiationError),
    NotConverged { report: SolveReport },
    MechanicalResidual { measured: f64, maximum: f64 },
}
impl core::fmt::Display for JointStructuralRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint structural radiation: {self:?}")
    }
}
impl std::error::Error for JointStructuralRadiationError {}
impl From<StructuralRadiationError> for JointStructuralRadiationError {
    fn from(value: StructuralRadiationError) -> Self {
        Self::Source(value)
    }
}
impl From<IterativeError> for JointStructuralRadiationError {
    fn from(value: IterativeError) -> Self {
        Self::Source(value.into())
    }
}
impl From<HarmonicResponseError> for JointStructuralRadiationError {
    fn from(value: HarmonicResponseError) -> Self {
        Self::Source(value.into())
    }
}
type Result<T> = core::result::Result<T, JointStructuralRadiationError>;

/// One immutable original pencil, trace and prepared Neumann FMM at one omega.
/// The native all-prescribed mechanical owner evaluates H*q and physical heat
/// WITHOUT factoring H. Those trial coordinates are Krylov unknowns, not added
/// physical supports. Every coordinate is solved in the simultaneous system.
pub struct PreparedJointStructuralRadiation<'a> {
    surface: &'a SpherePanels,
    traces: &'a [Vec<f64>],
    mechanical: PreparedSupportedHarmonicResponse,
    acoustic: IterativeRadiation<'a>,
    omega: f64,
    scales: Vec<f64>,
    pressure_scale: f64,
    zero_force: Vec<C64>,
    options: IterativeOptions,
    controls: JointStructuralRadiationControls,
    added_payload_bytes: usize,
}

impl<'a> PreparedJointStructuralRadiation<'a> {
    /// Admit the actual full M/C/K and original closed-surface trace before
    /// preparing the existing FMM owner. Gyroscopic and dark coordinates stay.
    ///
    /// Positive diagonal changes of solver coordinates use original M/K/C and
    /// rho*c*integral(T_j^2), with pressure scale sqrt(rho*c/area). They add no
    /// stiffness, damping, mass, force or radiation term to either equation.
    /// The acoustic right preconditioner is the original bounded near-block LU.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &'a SpherePanels,
        traces: &'a [Vec<f64>],
        pencil: StructuralRadiationPencil<'_>,
        omega: f64,
        medium: Medium,
        options: &IterativeOptions,
        controls: JointStructuralRadiationControls,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        check(gate, cx)?;
        let n = pencil.coordinates;
        let panels = surface.areas().len();
        let count = n
            .checked_mul(panels)
            .ok_or(invalid("joint trace extent overflow"))?;
        let total = n
            .checked_add(panels)
            .ok_or(invalid("joint coordinate extent overflow"))?;
        let restart = options.restart;
        let bytes = total
            .checked_mul(16)
            .and_then(|v| restart.checked_add(20).and_then(|r| v.checked_mul(r)))
            .and_then(|v| {
                restart
                    .checked_mul(restart.checked_add(4)?)
                    .and_then(|r| r.checked_mul(8))
                    .and_then(|r| v.checked_add(r))
            })
            .and_then(|v| n.checked_mul(512).and_then(|m| v.checked_add(m)))
            .and_then(|v| panels.checked_mul(160).and_then(|p| v.checked_add(p)))
            // Restart histories retained by the state, trial and native report.
            .and_then(|v| options.max_iters.checked_add(1)?.checked_mul(48)?.checked_add(v))
            .and_then(|v| v.checked_add(4096))
            .ok_or(invalid("joint solver payload extent overflow"))?;
        if n == 0
            || panels == 0
            || traces.len() != n
            || count > controls.maximum_trace_entries
            || bytes > controls.maximum_payload_bytes
            || !positive(omega)
            || !positive(medium.density)
            || !positive(medium.sound_speed)
            || !positive(medium.density * medium.sound_speed)
            || !tolerance(controls.maximum_relative_power_defect)
        {
            return Err(
                invalid("joint source, trace or numerical-payload admission failed").into(),
            );
        }
        for row in traces {
            check(gate, cx)?;
            if row.len() != panels || row.iter().any(|v| !v.is_finite()) {
                return Err(
                    invalid("one finite original trace per complete coordinate required").into(),
                );
            }
        }
        let mechanical = PreparedSupportedHarmonicResponse::try_new(
            pencil.mass,
            pencil.damping,
            pencil.stiffness,
            n,
            &(0..n).collect::<Vec<_>>(),
            omega,
            controls.mechanical,
            cx,
        )?;
        let acoustic = IterativeRadiation::prepare_bc(
            surface,
            omega / medium.sound_speed,
            medium,
            Formulation::BurtonMiller,
            options,
            Admittance::Rigid,
            gate,
        )?;
        let z = medium.density * medium.sound_speed;
        let mut area = Dd::ZERO;
        for &a in surface.areas() {
            area = area + Dd::from_f64(a);
        }
        let pressure_scale = det::sqrt(z) / det::sqrt(area.to_f64());
        let mut scales = filled(n, 0.0)?;
        for (j, scale) in scales.iter_mut().enumerate() {
            check(gate, cx)?;
            let mut trace_work = Dd::ZERO;
            for (&t, &a) in traces[j].iter().zip(surface.areas()) {
                trace_work = trace_work + Dd::from_f64(t) * Dd::from_f64(t) * Dd::from_f64(a);
            }
            let diagonal = Dd::from_f64(omega) * Dd::from_f64(pencil.mass[j * n + j])
                + Dd::from_f64(pencil.stiffness[j * n + j].abs()) / Dd::from_f64(omega)
                + Dd::from_f64(pencil.damping[j * n + j].abs())
                + Dd::from_f64(z) * trace_work;
            *scale = det::sqrt(diagonal.to_f64());
        }
        if !positive(pressure_scale) || scales.iter().any(|&v| !positive(v)) {
            return Err(invalid("joint physical coordinate scaling is unrepresentable").into());
        }
        check(gate, cx)?;
        Ok(Self {
            surface,
            traces,
            mechanical,
            acoustic,
            omega,
            scales,
            pressure_scale,
            zero_force: filled(n, C64::ZERO)?,
            options: options.clone(),
            controls,
            added_payload_bytes: bytes,
        })
    }

    pub const fn surface(&self) -> &'a SpherePanels {
        self.surface
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn coordinates(&self) -> usize {
        self.traces.len()
    }
    pub const fn added_payload_bytes(&self) -> usize {
        self.added_payload_bytes
    }

    pub fn solve(
        &self,
        force: &[C64],
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<JointStructuralRadiationResponse> {
        self.solve_controlled(force, gate, cx, |_| false)
    }

    /// One original force excites the complete simultaneous system. Cancel is
    /// polled before every composed matvec; Cx and native CancelGate retain
    /// their own polls inside mechanical/FMM actions. A refusal changes no
    /// owner or source and publishes no partial physical field.
    pub fn solve_controlled(
        &self,
        force: &[C64],
        gate: &CancelGate,
        cx: &Cx<'_>,
        mut cancelled: impl FnMut(usize) -> bool,
    ) -> Result<JointStructuralRadiationResponse> {
        check(gate, cx)?;
        let started = Instant::now();
        let n = self.coordinates();
        let panels = self.surface.areas().len();
        let total = n + panels;
        if force.len() != n || force.iter().any(|&v| !finite(v)) {
            return Err(
                invalid("one finite original force per complete coordinate required").into(),
            );
        }
        let mut rhs = filled(2 * total, 0.0)?;
        let mut amplitude = 0.0_f64;
        for (j, (&f, &d)) in force.iter().zip(&self.scales).enumerate() {
            rhs[panels + j] = f.re / d;
            rhs[total + panels + j] = f.im / d;
            amplitude = amplitude
                .max(rhs[panels + j].abs())
                .max(rhs[total + panels + j].abs());
        }
        if !amplitude.is_finite() || rhs.iter().any(|v| !v.is_finite()) {
            return Err(invalid("joint force scaling is unrepresentable").into());
        }
        let amplitude = if amplitude == 0.0 { 1.0 } else { amplitude };
        for v in &mut rhs {
            *v /= amplitude;
        }
        let op = JointOp {
            owner: self,
            gate,
            cx,
            error: RefCell::new(None),
            calls: Cell::new(0),
            cancelled: RefCell::new(&mut cancelled),
        };
        let restart = self.options.restart.min(rhs.len());
        let mut state = GmresState::new(&rhs, restart);
        // Leave numerical margin for independent unscaled block/work checks.
        // These gates keep the original tolerances; no physical term changes.
        let target = 0.125
            * self
                .options
                .tol
                .min(self.controls.mechanical.maximum_relative_residual)
                .min(self.controls.maximum_relative_power_defect);
        if !positive(target) {
            return Err(invalid("joint stopping tolerance underflowed").into());
        }
        loop {
            check(gate, cx)?;
            if state.iters >= self.options.max_iters {
                break;
            }
            let mut trial = state.clone();
            trial.restart = restart.min(self.options.max_iters - state.iters);
            let report = trial.run(&op, &rhs, target, 1, false);
            if let Some(error) = op.error.borrow_mut().take() {
                return Err(error);
            }
            trial.restart = restart;
            state = trial;
            if report.converged_euclidean() || !report.rel_residual.is_finite() {
                break;
            }
        }
        let report = SolveReport::from_claim(
            state.iters,
            state.residual_claim(),
            target,
            state.history.clone(),
        );
        if !report.converged_euclidean() || state.x.iter().any(|v| !v.is_finite()) {
            return Err(JointStructuralRadiationError::NotConverged { report });
        }
        let (mut pressure, mut velocity) = self.physical_coordinates(&state.x, gate)?;
        for v in pressure.iter_mut().chain(&mut velocity) {
            *v = v.scale(amplitude);
        }
        if pressure.iter().chain(&velocity).any(|&v| !finite(v)) {
            return Err(invalid("joint physical solution is unrepresentable").into());
        }
        let boundary_velocity = self.trace(&velocity)?;
        let (field, boundary_defect) =
            self.acoustic
                .finish_coupled_traces(pressure, boundary_velocity, gate)?;
        let reaction = self.traction(&field.pressure)?;
        let displacement = velocity
            .iter()
            .map(|v| C64::new(-v.im / self.omega, v.re / self.omega))
            .collect::<Vec<_>>();
        let motion = self.mechanical.solve(&self.zero_force, &displacement, cx)?;
        let mechanical_defect = self.force_defect(&motion.support_reaction, &reaction, force)?;
        if !mechanical_defect.is_finite()
            || mechanical_defect > self.controls.mechanical.maximum_relative_residual
        {
            return Err(JointStructuralRadiationError::MechanicalResidual {
                measured: mechanical_defect,
                maximum: self.controls.mechanical.maximum_relative_residual,
            });
        }
        let external_work = mean_work(force, &velocity)?;
        let interface_work = mean_work(&reaction, &velocity)?;
        let heat = motion.average_viscous_power_w;
        let radiation = field.power_diagnostics.surface_power.re;
        if field.radiated_power_roundoff_interval.1 < 0.0 {
            return Err(StructuralRadiationError::NonPassiveRadiation {
                signed_power_w: radiation,
                roundoff_interval_w: field.radiated_power_roundoff_interval,
            }
            .into());
        }
        let interface_defect = balance(interface_work, radiation, 0.0);
        let total_defect = balance(external_work, heat, radiation);
        if !interface_defect.is_finite()
            || !total_defect.is_finite()
            || interface_defect > self.controls.maximum_relative_power_defect
            || total_defect > self.controls.maximum_relative_power_defect
        {
            return Err(StructuralRadiationError::CoupledResidual {
                velocity: boundary_defect,
                interface_work: interface_defect,
                total_work: total_defect,
            }
            .into());
        }
        check(gate, cx)?;
        let receipt = self.acoustic.coupled_resource_receipt(
            op.calls.get() + 1,
            state.iters,
            started.elapsed().as_secs_f64(),
        )?;
        Ok(JointStructuralRadiationResponse {
            displacement: motion.displacement,
            velocity,
            field,
            report,
            receipt,
            opposing_pressure_force: reaction,
            external_input_power_w: external_work,
            structural_dissipation_w: heat,
            structure_to_fluid_power_w: interface_work,
            signed_radiation_power_w: radiation,
            relative_boundary_equation_defect: boundary_defect,
            relative_mechanical_equation_defect: mechanical_defect,
            relative_interface_work_defect: interface_defect,
            relative_total_work_defect: total_defect,
        })
    }

    fn physical_coordinates(&self, x: &[f64], gate: &CancelGate) -> Result<(Vec<C64>, Vec<C64>)> {
        let panels = self.surface.areas().len();
        let total = panels + self.coordinates();
        let pressure = (0..panels)
            .map(|i| C64::new(x[i], x[total + i]))
            .collect::<Vec<_>>();
        let pressure = self
            .acoustic
            .coupled_pressure_precondition(&pressure, gate)?
            .into_iter()
            .map(|p| p.scale(self.pressure_scale))
            .collect();
        let velocity = (0..self.coordinates())
            .map(|j| {
                C64::new(
                    x[panels + j] / self.scales[j],
                    x[total + panels + j] / self.scales[j],
                )
            })
            .collect();
        Ok((pressure, velocity))
    }
    fn force_defect(&self, required: &[C64], reaction: &[C64], force: &[C64]) -> Result<f64> {
        // A zero-force dark row can contain cancelling mechanical terms. Its
        // tiny rounded Hq is not a meaningful denominator by itself. Check the
        // complete original force block in the same physical diagonal metric
        // as the solve; under an original coordinate-unit congruence this norm
        // is unchanged. This is explicitly NOT a coefficientwise backward error.
        let magnitude = required
            .iter()
            .chain(reaction)
            .chain(force)
            .fold(0.0_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
        if magnitude == 0.0 {
            return Ok(0.0);
        }
        let mut a = filled(2 * self.coordinates(), 0.0)?;
        let mut b = a.clone();
        let mut f = a.clone();
        let mut residual = a.clone();
        for j in 0..self.coordinates() {
            let d = self.scales[j];
            for (part, (x, (y, z))) in [required[j].re, required[j].im]
                .into_iter()
                .zip(
                    [reaction[j].re, reaction[j].im]
                        .into_iter()
                        .zip([force[j].re, force[j].im]),
                )
                .enumerate()
            {
                let i = 2 * j + part;
                a[i] = (x / magnitude) / d;
                b[i] = (y / magnitude) / d;
                f[i] = (z / magnitude) / d;
                residual[i] =
                    (Dd::from_f64(a[i]) + Dd::from_f64(b[i]) - Dd::from_f64(f[i])).to_f64();
            }
        }
        let defect = fs_solver::norm2(&residual)
            / (fs_solver::norm2(&a) + fs_solver::norm2(&b) + fs_solver::norm2(&f));
        if !defect.is_finite() {
            return Err(invalid("joint original force block residual is unrepresentable").into());
        }
        Ok(defect)
    }
    fn apply(&self, x: &[f64], gate: &CancelGate, cx: &Cx<'_>) -> Result<Vec<f64>> {
        check(gate, cx)?;
        let panels = self.surface.areas().len();
        let total = panels + self.coordinates();
        let (pressure, velocity) = self.physical_coordinates(x, gate)?;
        let a = self.acoustic.apply(&pressure, gate)?;
        let b = self
            .acoustic
            .coupled_velocity_image(&self.trace(&velocity)?, gate)?;
        let displacement = velocity
            .iter()
            .map(|v| C64::new(-v.im / self.omega, v.re / self.omega))
            .collect::<Vec<_>>();
        let motion = self.mechanical.solve(&self.zero_force, &displacement, cx)?;
        let reaction = self.traction(&pressure)?;
        let mut output = filled(2 * total, 0.0)?;
        for i in 0..panels {
            output[i] =
                (Dd::from_f64(a[i].re) - Dd::from_f64(b[i].re)).to_f64() / self.pressure_scale;
            output[total + i] =
                (Dd::from_f64(a[i].im) - Dd::from_f64(b[i].im)).to_f64() / self.pressure_scale;
        }
        for j in 0..self.coordinates() {
            output[panels + j] = (Dd::from_f64(motion.support_reaction[j].re)
                + Dd::from_f64(reaction[j].re))
            .to_f64()
                / self.scales[j];
            output[total + panels + j] = (Dd::from_f64(motion.support_reaction[j].im)
                + Dd::from_f64(reaction[j].im))
            .to_f64()
                / self.scales[j];
        }
        if output.iter().any(|v| !v.is_finite()) {
            return Err(invalid("joint original operator overflow").into());
        }
        check(gate, cx)?;
        Ok(output)
    }
    fn trace(&self, velocity: &[C64]) -> Result<Vec<C64>> {
        let mut result = filled(self.surface.areas().len(), C64::ZERO)?;
        for (i, value) in result.iter_mut().enumerate() {
            let mut sum = ComplexSum::default();
            for (&v, row) in velocity.iter().zip(self.traces) {
                sum.scaled(v, row[i]);
            }
            *value = sum
                .finish()
                .map_err(|e| StructuralRadiationError::Boundary(IterativeError::Helmholtz(e)))?;
        }
        Ok(result)
    }
    fn traction(&self, pressure: &[C64]) -> Result<Vec<C64>> {
        let mut result = filled(self.coordinates(), C64::ZERO)?;
        for (value, row) in result.iter_mut().zip(self.traces) {
            let mut sum = ComplexSum::default();
            for ((&p, &t), &a) in pressure.iter().zip(row).zip(self.surface.areas()) {
                sum.scaled_product(p, t, a);
            }
            *value = sum
                .finish()
                .map_err(|e| StructuralRadiationError::Boundary(IterativeError::Helmholtz(e)))?;
        }
        Ok(result)
    }
}

/// Native field plus original-coordinate motion and independent physical work.
pub struct JointStructuralRadiationResponse {
    pub displacement: Vec<C64>,
    pub velocity: Vec<C64>,
    pub field: RadiationSolution,
    /// True Euclidean residual of the explicitly scaled simultaneous system.
    /// Original acoustic and mechanical block defects are reported separately.
    pub report: SolveReport,
    pub receipt: ResourceReceipt,
    pub opposing_pressure_force: Vec<C64>,
    pub external_input_power_w: f64,
    pub structural_dissipation_w: f64,
    pub structure_to_fluid_power_w: f64,
    pub signed_radiation_power_w: f64,
    pub relative_boundary_equation_defect: f64,
    /// Original force-block defect in the solver's physical diagonal metric:
    /// ||D^-1(Z_s v+reaction-F)||/(||D^-1 Z_s v||+||D^-1 reaction||+||D^-1 F||).
    /// Unlike the native mechanical owner's residual this is not coefficientwise.
    pub relative_mechanical_equation_defect: f64,
    pub relative_interface_work_defect: f64,
    pub relative_total_work_defect: f64,
}

struct JointOp<'a, 's, 'c> {
    owner: &'a PreparedJointStructuralRadiation<'s>,
    gate: &'a CancelGate,
    cx: &'a Cx<'c>,
    error: RefCell<Option<JointStructuralRadiationError>>,
    calls: Cell<usize>,
    cancelled: RefCell<&'a mut dyn FnMut(usize) -> bool>,
}
impl LinearOp for JointOp<'_, '_, '_> {
    fn n(&self) -> usize {
        2 * (self.owner.coordinates() + self.owner.surface.areas().len())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        if x.iter().any(|v| !v.is_finite()) {
            // A Krylov breakdown belongs to the GMRES report, rather than
            // being mislabeled as invalid physical source data by an owner.
            y.fill(f64::NAN);
            return;
        }
        if self.error.borrow().is_some() {
            y.fill(f64::NAN);
            return;
        }
        let result = if (self.cancelled.borrow_mut())(self.calls.get()) {
            Err(StructuralRadiationError::Cancelled.into())
        } else {
            self.calls.set(self.calls.get() + 1);
            self.owner.apply(x, self.gate, self.cx)
        };
        match result {
            Ok(v) => y.copy_from_slice(&v),
            Err(e) => {
                *self.error.borrow_mut() = Some(e);
                y.fill(f64::NAN);
            }
        }
    }
    fn apply_transpose(&self, _: &[f64], y: &mut [f64]) {
        *self.error.borrow_mut() =
            Some(invalid("joint radiation transpose is not implemented").into());
        y.fill(f64::NAN);
    }
}
fn check(gate: &CancelGate, cx: &Cx<'_>) -> Result<()> {
    checkpoint(cx)?;
    if gate.is_requested() {
        Err(StructuralRadiationError::Cancelled.into())
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests;
