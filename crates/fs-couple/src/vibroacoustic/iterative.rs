//! Full original-coordinate structural reaction through native FMM/Robin acoustics.
//!
//! Peak exp(-i omega t) fields obey H q = F - T^T W p and v_n = T (-i omega q).
//! T contains the original real normal-velocity traces and W the original panel
//! areas. With the native structural mobility Y = -i omega H^{-1}, the exterior
//! equation uses the EXACT nonlocal normalized admittance rho*c*T*Y*T^T*W.
//! The same map transfers pressure traction and velocity; no modal, diagonal,
//! local-wall or independent-radiator replacement is introduced.
//!
//! This mobility route needs an invertible original bare mechanical pencil.
//! A singular lossless bare resonance refuses; the existing jointly loaded
//! dense owners handle cases where radiation alone resolves that singularity.
//! Frozen M/C/K includes full consistent inertia, viscous and gyroscopic terms.
//! No frequency-dependent material sample is silently made into such a pencil.

use fs_bem::helmholtz::{
    Formulation, HelmholtzError, IterativeError, IterativeOptions, IterativeRadiation, Medium,
    RadiationSolution, ResourceReceipt,
    bc::{Admittance, BoundaryOperator},
};
use fs_bem::panel3d::SpherePanels;
use fs_exec::{CancelGate, Cx};
use fs_math::{c64::C64, dd::Dd};
use fs_modal::harmonic::{
    HarmonicResponse, HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_solver::krylov::SolveReport;

/// Original real row-major physical matrices, including every off-diagonal.
#[derive(Clone, Copy, Debug)]
pub struct StructuralRadiationPencil<'a> {
    pub coordinates: usize,
    pub mass: &'a [f64],
    pub damping: &'a [f64],
    pub stiffness: &'a [f64],
}

/// Explicit structural work, trace storage and coupled physical residual gates.
/// FMM, quadrature, near-block and GMRES limits remain the original options.
#[derive(Clone, Copy, Debug)]
pub struct StructuralRadiationControls {
    pub mechanical: HarmonicResponseControls,
    /// Original coordinates times acoustic panels, including zero/dark traces.
    pub maximum_trace_entries: usize,
    /// Original full structural mobility entries. No acoustic N-by-N array.
    pub maximum_mobility_entries: usize,
    /// Additional owned numerical payload beyond the original mechanical and
    /// FMM owners and borrowed traces. Includes preparation and one solve.
    pub maximum_payload_bytes: usize,
    pub maximum_relative_velocity_defect: f64,
    pub maximum_relative_power_defect: f64,
}

#[derive(Debug)]
pub enum StructuralRadiationError {
    Invalid(&'static str),
    Mechanical(HarmonicResponseError),
    Boundary(IterativeError),
    Cancelled,
    /// The native outward-rounded dot-product enclosure is strictly negative.
    /// This is a pointwise discrete refusal, not a spatial error certificate.
    NonPassiveRadiation {
        signed_power_w: f64,
        roundoff_interval_w: (f64, f64),
    },
    CoupledResidual {
        velocity: f64,
        interface_work: f64,
        total_work: f64,
    },
}
impl core::fmt::Display for StructuralRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "structural exterior radiation: {self:?}")
    }
}
impl std::error::Error for StructuralRadiationError {}
impl From<HarmonicResponseError> for StructuralRadiationError {
    fn from(value: HarmonicResponseError) -> Self {
        if matches!(value, HarmonicResponseError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Mechanical(value)
        }
    }
}
impl From<IterativeError> for StructuralRadiationError {
    fn from(value: IterativeError) -> Self {
        if matches!(value, IterativeError::Cancelled { .. }) {
            Self::Cancelled
        } else {
            Self::Boundary(value)
        }
    }
}

/// Source-bound immutable mechanical mobility and original surface traces.
/// Preparing this owner performs structural unit-force solves ONLY, using the
/// existing mechanical factor. It performs no acoustic column solves.
pub struct StructuralRadiation<'a> {
    surface: &'a SpherePanels,
    traces: &'a [Vec<f64>],
    mechanical: PreparedHarmonicResponse,
    mobility: Vec<C64>,
    omega: f64,
    medium: Medium,
    controls: StructuralRadiationControls,
    fingerprint: u64,
}

impl<'a> StructuralRadiation<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &'a SpherePanels,
        traces: &'a [Vec<f64>],
        pencil: StructuralRadiationPencil<'_>,
        angular_frequency_rad_s: f64,
        medium: Medium,
        controls: StructuralRadiationControls,
        cx: &Cx<'_>,
    ) -> Result<Self, StructuralRadiationError> {
        checkpoint(cx)?;
        let n = pencil.coordinates;
        let panels = surface.areas().len();
        let entries = n
            .checked_mul(n)
            .ok_or(invalid("mobility extent overflow"))?;
        let trace_entries = n
            .checked_mul(panels)
            .ok_or(invalid("trace extent overflow"))?;
        let payload = entries
            .checked_mul(16)
            .and_then(|v| n.checked_mul(192).and_then(|m| v.checked_add(m)))
            .and_then(|v| panels.checked_mul(160).and_then(|p| v.checked_add(p)))
            .and_then(|v| v.checked_add(4096))
            .ok_or(invalid("structural radiation payload overflow"))?;
        if n == 0
            || panels == 0
            || traces.len() != n
            || trace_entries > controls.maximum_trace_entries
            || entries > controls.maximum_mobility_entries
            || payload > controls.maximum_payload_bytes
        {
            return Err(invalid(
                "original trace, mobility or payload limit exceeded",
            ));
        }
        if !positive(angular_frequency_rad_s)
            || !positive(medium.density)
            || !positive(medium.sound_speed)
            || !positive(medium.density * medium.sound_speed)
            || !positive(angular_frequency_rad_s / medium.sound_speed)
            || !tolerance(controls.maximum_relative_velocity_defect)
            || !tolerance(controls.maximum_relative_power_defect)
        {
            return Err(invalid(
                "positive finite frequency/medium and residual gates in (0,1) required",
            ));
        }
        for row in traces {
            checkpoint(cx)?;
            if row.len() != panels || row.iter().any(|v| !v.is_finite()) {
                return Err(invalid(
                    "one finite original normal-velocity trace per coordinate required",
                ));
            }
        }
        if surface.areas().iter().any(|a| !positive(*a)) {
            return Err(invalid("positive finite original panel areas required"));
        }
        let mechanical = PreparedHarmonicResponse::try_new(
            pencil.mass,
            pencil.damping,
            pencil.stiffness,
            n,
            angular_frequency_rad_s,
            controls.mechanical,
            cx,
        )?;
        let mut mobility = filled(entries, C64::ZERO)?;
        let mut force = filled(n, C64::ZERO)?;
        for column in 0..n {
            checkpoint(cx)?;
            force[column] = C64::ONE;
            let response = mechanical.solve(&force, cx)?;
            force[column] = C64::ZERO;
            for row in 0..n {
                mobility[row * n + column] = response.velocity[row];
            }
        }
        let mut fingerprint = 0xcbf29ce484222325_u64;
        let mut word = |value: u64| {
            fingerprint = (fingerprint ^ value).wrapping_mul(0x100000001b3);
        };
        for value in [
            n as u64,
            panels as u64,
            angular_frequency_rad_s.to_bits(),
            medium.density.to_bits(),
            medium.sound_speed.to_bits(),
        ] {
            word(value);
        }
        for row in traces {
            for &v in row {
                word(v.to_bits());
            }
        }
        for &area in surface.areas() {
            word(area.to_bits());
        }
        for value in &mobility {
            word(value.re.to_bits());
            word(value.im.to_bits());
        }
        checkpoint(cx)?;
        Ok(Self {
            surface,
            traces,
            mechanical,
            mobility,
            omega: angular_frequency_rad_s,
            medium,
            controls,
            fingerprint,
        })
    }

    pub const fn surface(&self) -> &'a SpherePanels {
        self.surface
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn medium(&self) -> Medium {
        self.medium
    }
    pub fn coordinates(&self) -> usize {
        self.traces.len()
    }
    pub fn normal_velocity_traces(&self) -> &'a [Vec<f64>] {
        self.traces
    }
    /// Original row-major velocity/generalized-force mobility at this frequency.
    /// Nonreciprocal gyroscopic entries remain; this is not a broadband law.
    pub fn mobility(&self) -> &[C64] {
        &self.mobility
    }

    /// Prepare the existing resonance-safe acoustic FMM/GMRES owner once.
    /// Full Triangle and directional options are passed through unchanged.
    /// Preparation polls Cx around the original synchronous FMM stages; the
    /// separate native CancelGate remains active inside those stages.
    pub fn prepare<'s>(
        &'s self,
        options: &IterativeOptions,
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralRadiation<'s, 'a>, StructuralRadiationError> {
        checkpoint(cx)?;
        let operator = IterativeRadiation::prepare_bc(
            self.surface,
            self.omega / self.medium.sound_speed,
            self.medium,
            Formulation::BurtonMiller,
            options,
            Admittance::NonLocal {
                operator: self,
                declared_active: false,
            },
            gate,
        )?;
        checkpoint(cx)?;
        Ok(PreparedStructuralRadiation {
            source: self,
            operator,
            maximum_iterations: options.max_iters,
        })
    }

    fn trace(&self, velocity: &[C64]) -> Result<Vec<C64>, StructuralRadiationError> {
        let mut out = filled(self.surface.areas().len(), C64::ZERO)?;
        self.trace_into(velocity, &mut out)
            .map_err(|e| StructuralRadiationError::Boundary(IterativeError::Helmholtz(e)))?;
        Ok(out)
    }

    fn trace_into(&self, velocity: &[C64], out: &mut [C64]) -> Result<(), HelmholtzError> {
        if velocity.len() != self.coordinates() || out.len() != self.surface.areas().len() {
            return Err(boundary_invalid("normal-velocity trace dimension mismatch"));
        }
        for (panel, value) in out.iter_mut().enumerate() {
            let mut sum = ComplexSum::default();
            for (row, &v) in self.traces.iter().zip(velocity) {
                sum.scaled(v, row[panel]);
            }
            *value = sum.finish()?;
        }
        Ok(())
    }

    fn traction(&self, pressure: &[C64]) -> Result<Vec<C64>, HelmholtzError> {
        if pressure.len() != self.surface.areas().len() || pressure.iter().any(|v| !finite(*v)) {
            return Err(boundary_invalid(
                "one finite pressure per original panel required",
            ));
        }
        let mut force = boundary_filled(self.coordinates(), C64::ZERO)?;
        for (value, row) in force.iter_mut().zip(self.traces) {
            let mut sum = ComplexSum::default();
            for ((&p, &shape), &area) in pressure.iter().zip(row).zip(self.surface.areas()) {
                sum.scaled_product(p, shape, area);
            }
            *value = sum.finish()?;
        }
        Ok(force)
    }
}

impl BoundaryOperator for StructuralRadiation<'_> {
    fn len(&self) -> usize {
        self.surface.areas().len()
    }
    fn fingerprint(&self) -> u64 {
        self.fingerprint
    }
    fn apply(&self, pressure: &[C64], output: &mut [C64]) -> Result<(), HelmholtzError> {
        let force = self.traction(pressure)?;
        let n = self.coordinates();
        let mut velocity = boundary_filled(n, C64::ZERO)?;
        for (row, value) in velocity.iter_mut().enumerate() {
            let mut sum = ComplexSum::default();
            for (column, &f) in force.iter().enumerate() {
                sum.product(self.mobility[row * n + column], f);
            }
            *value = sum.finish()?;
        }
        self.trace_into(&velocity, output)?;
        let z = self.medium.density * self.medium.sound_speed;
        for value in output {
            *value = value.scale(z);
            if !finite(*value) {
                return Err(boundary_invalid(
                    "normalized structural admittance overflow",
                ));
            }
        }
        Ok(())
    }
}

/// Reusable acoustically loaded operator. Each force requires one native
/// acoustic solve and two original mechanical responses through the same factor; no acoustic
/// impedance matrix or acoustic unit-coordinate columns are constructed.
pub struct PreparedStructuralRadiation<'s, 'a> {
    source: &'s StructuralRadiation<'a>,
    operator: IterativeRadiation<'s>,
    maximum_iterations: usize,
}
impl<'s, 'a> PreparedStructuralRadiation<'s, 'a> {
    pub const fn source(&self) -> &'s StructuralRadiation<'a> {
        self.source
    }

    /// Solve full original force/torque phasors. Pressure opposes outward
    /// motion; positive reported interface work flows from structure to fluid.
    /// The native mechanical response's input work uses NET force after that
    /// reaction, while external_input_power_w below uses the original drive.
    pub fn solve(
        &self,
        external_force: &[C64],
        gate: &CancelGate,
        cx: &Cx<'_>,
    ) -> Result<StructuralRadiationResponse, StructuralRadiationError> {
        checkpoint(cx)?;
        if gate.is_requested() {
            return Err(StructuralRadiationError::Cancelled);
        }
        let original = self.source.mechanical.solve(external_force, cx)?;
        let drive = self.source.trace(&original.velocity)?;
        let mut state = self.operator.start(&drive, gate)?;
        self.operator.advance_controlled(
            &mut state,
            self.maximum_iterations.saturating_add(1),
            gate,
            &mut |_| cx.checkpoint().is_err(),
        )?;
        let (field, report, receipt) = self.operator.finish(&state)?;
        checkpoint(cx)?;
        let reaction = self
            .source
            .traction(&field.pressure)
            .map_err(|e| StructuralRadiationError::Boundary(IterativeError::Helmholtz(e)))?;
        let mut net_force = filled(external_force.len(), C64::ZERO)?;
        for ((net, &drive), &load) in net_force.iter_mut().zip(external_force).zip(&reaction) {
            *net = C64::new(
                (Dd::from_f64(drive.re) - Dd::from_f64(load.re)).to_f64(),
                (Dd::from_f64(drive.im) - Dd::from_f64(load.im)).to_f64(),
            );
        }
        let motion = self.source.mechanical.solve(&net_force, cx)?;
        let reconstructed = self.source.trace(&motion.velocity)?;
        // Form a dimensionless mismatch before norms: a physical unit change
        // must not turn an overflowing denominator into an artificial zero.
        let amplitude = reconstructed
            .iter()
            .chain(&field.velocity)
            .fold(0.0_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
        let mut mismatch = 0.0_f64;
        let mut velocity_scale = 0.0_f64;
        if amplitude > 0.0 {
            for (&a, &b) in reconstructed.iter().zip(&field.velocity) {
                let a = C64::new(a.re / amplitude, a.im / amplitude);
                let b = C64::new(b.re / amplitude, b.im / amplitude);
                mismatch = mismatch.hypot((a - b).abs());
                velocity_scale = velocity_scale.hypot(a.abs()).hypot(b.abs());
            }
        }
        let velocity_defect = relative(mismatch, velocity_scale);
        let external_work = mean_work(external_force, &motion.velocity)?;
        let interface_work = mean_work(&reaction, &motion.velocity)?;
        let radiation_work = field.power_diagnostics.surface_power.re;
        let structural_heat = motion.average_viscous_power_w;
        if field.radiated_power_roundoff_interval.1 < 0.0 {
            return Err(StructuralRadiationError::NonPassiveRadiation {
                signed_power_w: radiation_work,
                roundoff_interval_w: field.radiated_power_roundoff_interval,
            });
        }
        let interface_defect = balance(interface_work, radiation_work, 0.0);
        let total_defect = balance(external_work, structural_heat, radiation_work);
        if !velocity_defect.is_finite()
            || !interface_defect.is_finite()
            || !total_defect.is_finite()
            || velocity_defect > self.source.controls.maximum_relative_velocity_defect
            || interface_defect > self.source.controls.maximum_relative_power_defect
            || total_defect > self.source.controls.maximum_relative_power_defect
        {
            return Err(StructuralRadiationError::CoupledResidual {
                velocity: velocity_defect,
                interface_work: interface_defect,
                total_work: total_defect,
            });
        }
        checkpoint(cx)?;
        if gate.is_requested() {
            return Err(StructuralRadiationError::Cancelled);
        }
        Ok(StructuralRadiationResponse {
            motion,
            field,
            report,
            receipt,
            opposing_pressure_force: reaction,
            external_input_power_w: external_work,
            structural_dissipation_w: structural_heat,
            structure_to_fluid_power_w: interface_work,
            signed_radiation_power_w: radiation_work,
            relative_velocity_defect: velocity_defect,
            relative_interface_work_defect: interface_defect,
            relative_total_work_defect: total_defect,
        })
    }
}

pub struct StructuralRadiationResponse {
    pub motion: HarmonicResponse,
    /// Same bound pressure/velocity field consumed by native BEM observers.
    pub field: RadiationSolution,
    pub report: SolveReport,
    pub receipt: ResourceReceipt,
    pub opposing_pressure_force: Vec<C64>,
    pub external_input_power_w: f64,
    pub structural_dissipation_w: f64,
    pub structure_to_fluid_power_w: f64,
    pub signed_radiation_power_w: f64,
    pub relative_velocity_defect: f64,
    pub relative_interface_work_defect: f64,
    pub relative_total_work_defect: f64,
}

#[derive(Default)]
struct ComplexSum {
    re: Dd,
    im: Dd,
}
impl ComplexSum {
    fn scaled(&mut self, z: C64, s: f64) {
        self.re = self.re + Dd::from_f64(z.re) * Dd::from_f64(s);
        self.im = self.im + Dd::from_f64(z.im) * Dd::from_f64(s);
    }
    fn scaled_product(&mut self, z: C64, a: f64, b: f64) {
        self.re = self.re + Dd::from_f64(z.re) * Dd::from_f64(a) * Dd::from_f64(b);
        self.im = self.im + Dd::from_f64(z.im) * Dd::from_f64(a) * Dd::from_f64(b);
    }
    fn product(&mut self, a: C64, b: C64) {
        self.re = self.re + Dd::from_f64(a.re) * Dd::from_f64(b.re)
            - Dd::from_f64(a.im) * Dd::from_f64(b.im);
        self.im = self.im
            + Dd::from_f64(a.re) * Dd::from_f64(b.im)
            + Dd::from_f64(a.im) * Dd::from_f64(b.re);
    }
    fn finish(self) -> Result<C64, HelmholtzError> {
        let value = C64::new(self.re.to_f64(), self.im.to_f64());
        if finite(value) {
            Ok(value)
        } else {
            Err(boundary_invalid(
                "physical pressure/velocity contraction overflow",
            ))
        }
    }
}
fn mean_work(a: &[C64], b: &[C64]) -> Result<f64, StructuralRadiationError> {
    let mut sum = Dd::ZERO;
    for (&a, &b) in a.iter().zip(b) {
        sum =
            sum + Dd::from_f64(a.re) * Dd::from_f64(b.re) + Dd::from_f64(a.im) * Dd::from_f64(b.im);
    }
    let value = (sum * Dd::from_f64(0.5)).to_f64();
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("physical structural work overflow"))
    }
}
fn relative(error: f64, scale: f64) -> f64 {
    if scale == 0.0 {
        if error == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        error / scale
    }
}
fn balance(a: f64, b: f64, c: f64) -> f64 {
    let scale = a.abs().max(b.abs()).max(c.abs());
    if scale == 0.0 {
        return 0.0;
    }
    let [a, b, c] = [a / scale, b / scale, c / scale];
    relative(
        (Dd::from_f64(a) - Dd::from_f64(b) - Dd::from_f64(c))
            .to_f64()
            .abs(),
        a.abs() + b.abs() + c.abs(),
    )
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn tolerance(v: f64) -> bool {
    positive(v) && v < 1.0
}
fn invalid(s: &'static str) -> StructuralRadiationError {
    StructuralRadiationError::Invalid(s)
}
fn boundary_invalid(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), StructuralRadiationError> {
    cx.checkpoint()
        .map_err(|_| StructuralRadiationError::Cancelled)
}
fn filled<T: Clone>(length: usize, value: T) -> Result<Vec<T>, StructuralRadiationError> {
    let mut out = Vec::new();
    out.try_reserve_exact(length)
        .map_err(|_| invalid("structural radiation allocation failed"))?;
    out.resize(length, value);
    Ok(out)
}
fn boundary_filled<T: Clone>(length: usize, value: T) -> Result<Vec<T>, HelmholtzError> {
    let mut out = Vec::new();
    out.try_reserve_exact(length)
        .map_err(|_| boundary_invalid("nonlocal structural allocation failed"))?;
    out.resize(length, value);
    Ok(out)
}

#[cfg(test)]
mod tests;
