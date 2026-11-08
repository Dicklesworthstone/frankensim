//! Incident sound and material dispersion in the original shell coordinates.
//!
//! The native screen supplies the blocked force f_b = -G[p_b] and the complete
//! motion impedance Z. The original tensor-Prony shell solves its material
//! pencil with load Z and force f_external + f_b. The total fluid reaction is
//! -G([p_b] + [p_motion]); the incident field itself has no pressure jump across
//! the zero-offset sheet. Both wet faces are therefore counted exactly once.
//!
//! All fields use peak exp(-i omega t) phasors. Molecular compression changes
//! k(omega), while the physical Neumann drive remains i*rho*omega*v. Signed
//! exterior work includes gas absorption and radiation, and can be negative
//! when sound drives the shell. It is not a far-field power or cross-section.
//! No modes, structural material branches, physical coordinates or mutual
//! reaction terms are truncated. This is a stationary harmonic composition,
//! not a causal realization, moving geometry or measured scattering claim.

use super::{
    PreparedOpenTensorShellRadiation, TensorShellRadiationError, admit, checkpoint, reserved,
};
use fs_bem::helmholtz::{bc::IncidentField, molecular::MolecularExterior};
use fs_bem::near_field::{Options as ReceiverOptions, VectorObservation};
use fs_bem::reciprocal_radiation::open_screen::{
    MolecularScreenIllumination, MolecularScreenIlluminationLimits,
    PreparedMolecularScreenIlluminationReceivers,
};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicPortControls, HarmonicResponse, HarmonicResponseControls};
use fs_plate::shell::viscosity::harmonic::PreparedTensorPronyShell;

/// Complete original-coordinate preparation and per-response resource limits.
#[derive(Clone, Copy, Debug)]
pub struct OpenTensorShellIlluminationControls {
    pub mechanical: HarmonicResponseControls,
    pub ports: HarmonicPortControls,
    /// One native boundary factor for all motion and incident source columns.
    /// Includes actual source clearance, quadrature work and parallel scratch.
    pub boundary: MolecularScreenIlluminationLimits,
    /// Additional logical live payload for effective force, native result,
    /// reaction and coherent jump scratch: 320*q + 64*V + 4096 bytes.
    /// Original geometry, boundary bank, mechanical factor/strain recovery,
    /// receiver rows and allocator metadata have their separate native limits.
    pub maximum_response_bytes: usize,
    /// Per-force additional work score: 96*q*q + 64*q*V + 512*q + 4096.
    /// Native mechanical solution/strain recovery retains its own admission.
    pub maximum_response_work: u64,
    /// Relative work/traction defects, finite and strictly between zero and one.
    pub power_relative_tolerance: f64,
}

/// One material-shell factor and one source-bound, two-sided boundary bank.
pub struct PreparedOpenTensorShellIllumination<'p, 'geometry, 'shell, 'incident> {
    radiation: MolecularScreenIllumination<'p, 'geometry, 'incident>,
    mechanical: PreparedTensorPronyShell<'shell>,
    tolerance: f64,
}

impl<'geometry, 'shell> PreparedOpenTensorShellRadiation<'geometry, 'shell> {
    /// Prepare the full material shell under one or more independent incident
    /// fields. Each field may itself contain coherent source terms. The exact
    /// original source gas must match the geometry's equilibrium rho/c, and
    /// both the mechanical material band and acoustic preparation band apply.
    ///
    /// An explicit zero-amplitude incident field selects force-only molecular
    /// radiation. Repeated external force vectors reuse both factors; no new
    /// fit, material evaluation, boundary solve or source sampling is needed.
    pub fn prepare_molecular_illumination<'p, 'incident>(
        &'p self,
        exterior: &MolecularExterior,
        omega_rad_s: f64,
        incidents: &[IncidentField<'incident>],
        controls: OpenTensorShellIlluminationControls,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedOpenTensorShellIllumination<'p, 'geometry, 'shell, 'incident>,
        TensorShellRadiationError,
    > {
        checkpoint(cx)?;
        let gas = exterior.source().gas();
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || omega_rad_s > self.highest_omega_rad_s
            || gas.density.to_bits() != self.source.medium.density.to_bits()
            || gas.sound_speed.to_bits() != self.source.medium.sound_speed.to_bits()
            || !controls.power_relative_tolerance.is_finite()
            || controls.power_relative_tolerance <= 0.0
            || controls.power_relative_tolerance >= 1.0
        {
            return Err(TensorShellRadiationError::Invalid(
                "illuminated material shell requires matched gas, admitted frequency and power tolerance",
            ));
        }
        let q = self.source.velocity_rows.len();
        let (n, v) = (q as u128, self.source.surface.vertices().len() as u128);
        admit(
            "illuminated shell coordinates",
            n,
            controls.mechanical.maximum_coordinates as u128,
        )?;
        admit(
            "illuminated shell ports",
            n,
            controls.ports.maximum_ports as u128,
        )?;
        admit(
            "illuminated shell response bytes",
            320 * n + 64 * v + 4096,
            controls.maximum_response_bytes as u128,
        )?;
        admit(
            "illuminated shell response work",
            96 * n * n + 64 * n * v + 512 * n + 4096,
            controls.maximum_response_work as u128,
        )?;
        let mut rows = reserved(q, "illuminated shell source rows")?;
        rows.extend(self.source.velocity_rows.iter().map(Vec::as_slice));
        let radiation = self.screen.solve_molecular_illumination(
            exterior,
            omega_rad_s,
            &rows,
            incidents,
            controls.boundary,
            || cx.checkpoint().is_err(),
        )?;
        let mechanical = self
            .source
            .shell
            .prepare_harmonic_response_with_port_impedance(
                omega_rad_s,
                &self.source.identity_ports,
                q,
                radiation.motion().impedance(),
                controls.mechanical,
                controls.ports,
                cx,
            )?;
        checkpoint(cx)?;
        Ok(PreparedOpenTensorShellIllumination {
            radiation,
            mechanical,
            tolerance: controls.power_relative_tolerance,
        })
    }
}

/// Original physical powers [W]; the total exterior term is signed.
#[derive(Clone, Copy, Debug)]
pub struct OpenTensorShellIlluminationPower {
    pub external_mechanical_input_w: f64,
    pub blocked_acoustic_input_w: f64,
    /// Original tensor-Prony and Kelvin heat, excluding exterior molecular heat.
    pub intrinsic_material_heat_w: f64,
    /// Work of motion-only radiation, including molecular heat in the gas.
    pub motion_exterior_work_w: f64,
    /// Direct work of the coherent TOTAL pressure-jump reaction. Can be negative.
    pub total_exterior_work_w: f64,
    /// Sum of absolute force/velocity products and intrinsic material heat.
    /// This retains reactive work when net active input is close to zero.
    pub apparent_power_scale_w: f64,
    /// Maximum defect in external = heat + total exterior, effective input =
    /// heat + motion exterior, and the independently observed native input.
    pub relative_balance: f64,
    pub relative_motion_work_defect: f64,
    /// Maximum componentwise defect in total reaction = blocked force - Z*v.
    pub relative_traction_defect: f64,
}

/// Motion and field share the same immutable source owner and incident case.
pub struct OpenTensorShellIlluminationResponse<'r, 'p, 'geometry, 'incident> {
    radiation: &'r MolecularScreenIllumination<'p, 'geometry, 'incident>,
    case: usize,
    motion: HarmonicResponse,
    total_fluid_reaction: Vec<C64>,
    power: OpenTensorShellIlluminationPower,
}

impl<'p, 'geometry, 'shell, 'incident>
    PreparedOpenTensorShellIllumination<'p, 'geometry, 'shell, 'incident>
{
    pub fn radiation(&self) -> &MolecularScreenIllumination<'p, 'geometry, 'incident> {
        &self.radiation
    }
    pub fn shell_preparation(&self) -> &PreparedTensorPronyShell<'shell> {
        &self.mechanical
    }
    pub fn impedance(&self) -> &[C64] {
        self.radiation.motion().impedance()
    }

    /// Solve one incident case plus a complete original physical force/torque
    /// vector. The native material owner admits strain over the full harmonic
    /// cycle, and source/total reaction work is independently checked before
    /// the response escapes. A zero external force allows sound-driven motion.
    pub fn solve(
        &self,
        case: usize,
        external_force: &[C64],
        cx: &Cx<'_>,
    ) -> Result<
        OpenTensorShellIlluminationResponse<'_, 'p, 'geometry, 'incident>,
        TensorShellRadiationError,
    > {
        checkpoint(cx)?;
        let blocked =
            self.radiation
                .blocked_forces()
                .get(case)
                .ok_or(TensorShellRadiationError::Invalid(
                    "illuminated shell incident case out of range",
                ))?;
        if external_force.len() != blocked.len() || external_force.iter().any(|z| !finite(*z)) {
            return Err(TensorShellRadiationError::Invalid(
                "illuminated shell requires a complete finite physical force vector",
            ));
        }
        let mut effective = reserved(blocked.len(), "illuminated shell effective force")?;
        for (&external, &acoustic) in external_force.iter().zip(blocked) {
            checkpoint(cx)?;
            let force = external + acoustic;
            if !finite(force) {
                return Err(TensorShellRadiationError::Invalid(
                    "illuminated shell effective force overflow",
                ));
            }
            effective.push(force);
        }
        let motion = self.mechanical.solve(&effective, cx)?;
        let total = self
            .radiation
            .total_fluid_reaction(case, &motion.velocity, || cx.checkpoint().is_err())?;
        let port = motion
            .port_load
            .as_ref()
            .ok_or(TensorShellRadiationError::Invalid(
                "illuminated shell omitted original physical port response",
            ))?;
        let motion_exterior = self
            .radiation
            .motion()
            .boundary_work_w(&motion.velocity, || cx.checkpoint().is_err())?;
        let mut external_work = Sum::default();
        let mut blocked_work = Sum::default();
        let mut total_work = Sum::default();
        let mut apparent = Sum::default();
        let mut traction_defect = 0.0_f64;
        for i in 0..blocked.len() {
            checkpoint(cx)?;
            let velocity = motion.velocity[i];
            let (external, acoustic, effort, reaction) =
                (external_force[i], blocked[i], port.effort[i], total[i]);
            external_work.add(0.5 * (external * velocity.conj()).re);
            blocked_work.add(0.5 * (acoustic * velocity.conj()).re);
            total_work.add(-0.5 * (reaction * velocity.conj()).re);
            apparent.add(
                0.5 * (external.abs() + acoustic.abs() + effort.abs() + reaction.abs())
                    * velocity.abs(),
            );
            let scale = acoustic.abs() + effort.abs() + reaction.abs();
            let defect = (reaction - (acoustic - effort)).abs();
            traction_defect = traction_defect.max(relative(defect, scale));
        }
        apparent.add(motion.average_viscous_power_w.abs());
        let (external_work, blocked_work, total_work, scale) = (
            external_work.value(),
            blocked_work.value(),
            total_work.value(),
            apparent.value(),
        );
        let heat = motion.average_viscous_power_w;
        let coupled = relative((external_work - heat - total_work).abs(), scale)
            .max(relative(
                (external_work + blocked_work - heat - motion_exterior).abs(),
                scale,
            ))
            .max(relative(
                (motion.average_input_power_w - external_work - blocked_work).abs(),
                scale,
            ));
        let motion_work_defect = relative(
            (port.average_absorbed_power_w - motion_exterior).abs(),
            scale,
        );
        if ![
            external_work,
            blocked_work,
            total_work,
            scale,
            heat,
            motion_exterior,
            coupled,
            motion_work_defect,
            traction_defect,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(TensorShellRadiationError::Invalid(
                "illuminated shell work is not representable",
            ));
        }
        if coupled > self.tolerance
            || motion_work_defect > self.tolerance
            || traction_defect > self.tolerance
        {
            return Err(TensorShellRadiationError::IlluminationBalance {
                coupled,
                motion_work: motion_work_defect,
                traction: traction_defect,
                tolerance: self.tolerance,
            });
        }
        let power = OpenTensorShellIlluminationPower {
            external_mechanical_input_w: external_work,
            blocked_acoustic_input_w: blocked_work,
            intrinsic_material_heat_w: heat,
            motion_exterior_work_w: motion_exterior,
            total_exterior_work_w: total_work,
            apparent_power_scale_w: scale,
            relative_balance: coupled,
            relative_motion_work_defect: motion_work_defect,
            relative_traction_defect: traction_defect,
        };
        checkpoint(cx)?;
        Ok(OpenTensorShellIlluminationResponse {
            radiation: &self.radiation,
            case,
            motion,
            total_fluid_reaction: total,
            power,
        })
    }

    /// Prepare the same total-pressure/Euler-velocity rows for repeated forces.
    /// Analytic incident fields are cached once and added exactly once, after
    /// blocked and motion jumps are combined coherently. Limits are the native
    /// row limits, including this wrapper's additional source reference.
    pub fn prepare_receivers(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: ReceiverOptions,
        maximum_row_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<
        PreparedOpenTensorShellIlluminationReceivers<'_, 'p, 'geometry, 'incident>,
        TensorShellRadiationError,
    > {
        checkpoint(cx)?;
        let inner = maximum_row_bytes
            .checked_sub(core::mem::size_of::<&MolecularScreenIllumination<'_, '_, '_>>())
            .ok_or(TensorShellRadiationError::Invalid(
                "illuminated shell receiver bytes are insufficient",
            ))?;
        let rows = self.radiation.prepare_receivers(
            points_m,
            minimum_clearance_m,
            options,
            inner,
            || cx.checkpoint().is_err(),
        )?;
        checkpoint(cx)?;
        Ok(PreparedOpenTensorShellIlluminationReceivers {
            source: &self.radiation,
            rows,
        })
    }
}

impl OpenTensorShellIlluminationResponse<'_, '_, '_, '_> {
    pub fn motion(&self) -> &HarmonicResponse {
        &self.motion
    }
    pub fn blocked_acoustic_force(&self) -> &[C64] {
        &self.radiation.blocked_forces()[self.case]
    }
    pub fn total_fluid_reaction(&self) -> &[C64] {
        &self.total_fluid_reaction
    }
    pub const fn power(&self) -> &OpenTensorShellIlluminationPower {
        &self.power
    }
    pub fn total_pressure_jump(&self, cx: &Cx<'_>) -> Result<Vec<C64>, TensorShellRadiationError> {
        checkpoint(cx)?;
        Ok(self
            .radiation
            .total_pressure_jump(self.case, &self.motion.velocity, || {
                cx.checkpoint().is_err()
            })?)
    }
}

/// Source-bound receiver bank for actual incident/force-driven shell responses.
pub struct PreparedOpenTensorShellIlluminationReceivers<'r, 'p, 'geometry, 'incident> {
    source: &'r MolecularScreenIllumination<'p, 'geometry, 'incident>,
    rows: PreparedMolecularScreenIlluminationReceivers<'r, 'p, 'geometry, 'incident>,
}
impl<'p, 'geometry, 'incident>
    PreparedOpenTensorShellIlluminationReceivers<'_, 'p, 'geometry, 'incident>
{
    pub fn points_m(&self) -> &[[f64; 3]] {
        self.rows.points_m()
    }
    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }
    /// A response from another boundary preparation refuses even if its sizes
    /// and frequencies agree. The total field includes acoustic interference.
    pub fn evaluate(
        &self,
        response: &OpenTensorShellIlluminationResponse<'_, 'p, 'geometry, 'incident>,
        cx: &Cx<'_>,
    ) -> Result<VectorObservation, TensorShellRadiationError> {
        checkpoint(cx)?;
        if !core::ptr::eq(self.source, response.radiation) {
            return Err(TensorShellRadiationError::Invalid(
                "illuminated shell receiver and response sources differ",
            ));
        }
        Ok(self
            .rows
            .evaluate(response.case, &response.motion.velocity, || {
                cx.checkpoint().is_err()
            })?)
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn relative(defect: f64, scale: f64) -> f64 {
    if !defect.is_finite() || !scale.is_finite() || scale < 0.0 {
        return f64::INFINITY;
    }
    if scale > 0.0 {
        defect / scale
    } else if defect == 0.0 {
        0.0
    } else {
        f64::INFINITY
    }
}
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.value + value;
        self.correction += if self.value.abs() >= value.abs() {
            (self.value - next) + value
        } else {
            (value - next) + self.value
        };
        self.value = next;
    }
    fn value(self) -> f64 {
        self.value + self.correction
    }
}
