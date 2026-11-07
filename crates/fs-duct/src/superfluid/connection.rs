//! Reciprocal serial connections, including an ideal entropy-blocking superleak.
//!
//! The right end of the first channel joins the left end of the second. Both
//! retain their original material, wave, wall and vortex owners. The connection
//! owns the omitted endpoint half-volume inertia. Eliminating their two native
//! endpoint responses requires only a constant-size symmetric midpoint solve.
//!
//! This is a linear, inviscid plug-flow junction with the SAME equilibrium
//! material on both sides. No entrance loss, junction wall drag, vortex-line
//! transfer, finite-temperature relinearization or nonlinear advective momentum
//! is inferred. An area change or superleak constraint has a stationary support
//! reaction, so energy conservation does not imply unsupported axial momentum
//! conservation. Interior vortex reservoirs remain separate explicit ports.
//!
//! Landau's currents and d(mu)=dp/rho-s*dT give the power identity
//! `mu'*J_mass + T'*J_entropy = p'*U + rho_s*s*T'*W`, with
//! `U=A*v` and `W=A*(v_n-v_s)`. See Landau, Phys. Rev. 60, 356 (1941),
//! <https://doi.org/10.1103/PhysRev.60.356>, and Soulaine et al., CPC 187
//! (2015), <https://doi.org/10.1016/j.cpc.2014.10.006>.

use super::*;
use fs_math::dd::Dd;

/// Which component currents the stationary connection admits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConnectionKind {
    /// Both components pass, conserving their area-integrated currents.
    Open,
    /// Normal flow and transported entropy are blocked. Only superfluid passes.
    /// This is an ideal constraint, not a porous-flow or critical-velocity law.
    Superleak,
}

/// Original channel sources and the two unconnected outer ports for one step.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConnectionDrive<'a> {
    /// Left is the first channel's left end; right is the second's right end.
    /// Both use their channels' +x orientation, as in `BoundaryVelocities`.
    pub outer: BoundaryVelocities,
    /// Original component body forces on each channel, including the connected
    /// endpoint half volumes. Outer prescribed half volumes remain external.
    pub forces: [BodyForce<'a>; 2],
    /// Original interior friction closures, in channel order.
    pub friction: [MutualFriction; 2],
}

/// One atomic interval of the native channels and their connection inertia.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConnectionStep {
    /// Unmodified native reports, in first/second channel order.
    pub channels: [ChannelStep; 2],
    /// Instantaneous barycentric/relative volume currents, positive first to second.
    pub flow_before_m3_s: [f64; 2],
    /// Interval midpoint currents imposed on the two connected ports.
    pub flow_midpoint_m3_s: [f64; 2],
    /// Instantaneous currents after the interval; channel endpoint samples keep
    /// the imposed MIDPOINT velocities, while this owner keeps endpoint inertia.
    pub flow_after_m3_s: [f64; 2],
    /// Kinetic storage of the two connected endpoint half volumes.
    pub junction_energy_before_j: f64,
    /// Kinetic storage after the interval, in J.
    pub junction_energy_after_j: f64,
    /// Pressure/temperature port work entering the junction from both channels.
    pub junction_wave_work_j: f64,
    /// Authored component-force work on the connected half volumes.
    pub junction_body_work_j: f64,
    /// Wave work entering the assembly through its two unconnected ends.
    pub outer_wave_work_j: f64,
    /// Junction storage change minus its port and body work.
    pub junction_balance_defect_j: f64,
    /// Native wave/vortex/junction storage change plus exported heat, minus
    /// external wave/body work and the original vortex-reservoir energy ports.
    pub combined_energy_defect_j: f64,
    /// Sum of channel mass imports minus the two external mass imports.
    pub mass_balance_defect_kg: f64,
    /// Corresponding FIRST-ORDER entropy-current residual. Frictional entropy
    /// production belongs to the caloric receiver of the separate heat arrays.
    pub entropy_balance_defect_j_k: f64,
}

/// Finite-inertia connection bound to two channel geometries and one material.
///
/// Preparation is allocation-free. Every wave and coupling gate passes before
/// either channel, the junction clock or its instantaneous flow is published.
#[derive(Clone, Debug)]
pub struct TwoFluidConnection {
    kind: ConnectionKind,
    material: TwoFluidState,
    lengths_m: [f64; 2],
    areas_m2: [f64; 2],
    cells: [usize; 2],
    inertance: [f64; 2],
    flow_m3_s: [f64; 2],
    last_midpoint_m3_s: [f64; 2],
    time_s: f64,
}

impl TwoFluidConnection {
    /// Attach already compatible currents without projecting away stored energy.
    /// An ideal superleak additionally requires zero initial normal current.
    pub fn try_new(
        first: &TwoFluidChannel,
        second: &TwoFluidChannel,
        kind: ConnectionKind,
    ) -> Result<Self, ChannelError> {
        if first.material != second.material || first.time_s != second.time_s {
            return Err(ChannelError::Invalid(
                "connected channels require identical reference material and clocks",
            ));
        }
        let first_flow = endpoint_flow(first, first.cells())?;
        let second_flow = endpoint_flow(second, 0)?;
        if !same_flow(first_flow, second_flow) {
            return Err(ChannelError::Invalid(
                "initial connected component volume currents must already agree",
            ));
        }
        let mut flow = [0, 1].map(|i| 0.5 * first_flow[i] + 0.5 * second_flow[i]);
        if kind == ConnectionKind::Superleak {
            let constrained = [-first.material.superfluid_fraction() * flow[1], flow[1]];
            if !same_flow(flow, constrained) {
                return Err(ChannelError::Invalid(
                    "a superleak cannot discard an initial normal-fluid current",
                ));
            }
            flow = constrained;
        }
        let geometric = 0.5 * first.dx / first.area_m2 + 0.5 * second.dx / second.area_m2;
        let rho = first.material.density_kg_m3();
        let reduced = rho * first.material.superfluid_fraction() * first.material.normal_fraction();
        let inertance = [rho * geometric, reduced * geometric];
        if !positive(geometric) || !inertance.into_iter().all(positive) {
            return Err(ChannelError::Invalid(
                "junction inertia is not representable",
            ));
        }
        let connection = Self {
            kind,
            material: first.material,
            lengths_m: [first.length_m, second.length_m],
            areas_m2: [first.area_m2, second.area_m2],
            cells: [first.cells(), second.cells()],
            inertance,
            flow_m3_s: flow,
            last_midpoint_m3_s: flow,
            time_s: first.time_s,
        };
        connection.kinetic_energy_j()?;
        Ok(connection)
    }

    /// Authored current constraint.
    pub const fn kind(&self) -> ConnectionKind {
        self.kind
    }

    /// Accepted time shared by the connection and both channels.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Instantaneous `[A*v, A*(v_n-v_s)]`, positive first to second, in m³/s.
    pub const fn flow_m3_s(&self) -> [f64; 2] {
        self.flow_m3_s
    }

    /// Physical kinetic energy excluded from the two native channel interiors.
    pub fn kinetic_energy_j(&self) -> Result<f64, ChannelError> {
        kinetic(self.inertance, self.flow_m3_s)
    }

    /// Advance both native channels and the finite-inertia connection together.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &mut self,
        first: &mut TwoFluidChannel,
        second: &mut TwoFluidChannel,
        dt_s: f64,
        drive: ConnectionDrive<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConnectionStep, ChannelError> {
        match self.step_with_heat(first, second, dt_s, drive, cancelled, |_, _, _| {
            Ok::<(), ChannelError>(())
        }) {
            Ok((report, ())) => Ok(report),
            Err(HeatStepError::Channel(error) | HeatStepError::Heat(error)) => Err(error),
        }
    }

    /// Prepare a thermal candidate using BOTH original cell-heat arrays.
    ///
    /// The callback must not mutate accepted external state. Its refusal, a
    /// second-channel refusal or final cancellation publishes NONE of the three
    /// owners. After success its returned candidate may be installed infallibly.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_heat<T, E>(
        &mut self,
        first: &mut TwoFluidChannel,
        second: &mut TwoFluidChannel,
        dt_s: f64,
        drive: ConnectionDrive<'_>,
        mut cancelled: impl FnMut() -> bool,
        prepare_heat: impl FnOnce(&ConnectionStep, &[f64], &[f64]) -> Result<T, E>,
    ) -> Result<(ConnectionStep, T), HeatStepError<E>> {
        let report = self
            .stage(first, second, dt_s, drive, &mut cancelled)
            .map_err(HeatStepError::Channel)?;
        let candidate = prepare_heat(
            &report,
            &first.pending_dissipation_heat_j,
            &second.pending_dissipation_heat_j,
        )
        .map_err(HeatStepError::Heat)?;
        poll(&mut cancelled).map_err(HeatStepError::Channel)?;
        first.publish_step(report.channels[0]);
        second.publish_step(report.channels[1]);
        self.flow_m3_s = report.flow_after_m3_s;
        self.last_midpoint_m3_s = report.flow_midpoint_m3_s;
        self.time_s = report.channels[0].time_s;
        Ok((report, candidate))
    }

    fn validate_binding(
        &self,
        first: &TwoFluidChannel,
        second: &TwoFluidChannel,
    ) -> Result<(), ChannelError> {
        for (i, channel) in [first, second].into_iter().enumerate() {
            if channel.material != self.material
                || channel.time_s != self.time_s
                || channel.length_m != self.lengths_m[i]
                || channel.area_m2 != self.areas_m2[i]
                || channel.cells() != self.cells[i]
            {
                return Err(ChannelError::Invalid(
                    "junction material, geometry and channel clocks must remain bound",
                ));
            }
        }
        if !same_flow(
            endpoint_flow(first, first.cells())?,
            self.last_midpoint_m3_s,
        ) || !same_flow(endpoint_flow(second, 0)?, self.last_midpoint_m3_s)
        {
            return Err(ChannelError::Invalid(
                "connected endpoint state changed outside its junction owner",
            ));
        }
        Ok(())
    }

    fn stage(
        &self,
        first: &mut TwoFluidChannel,
        second: &mut TwoFluidChannel,
        dt: f64,
        drive: ConnectionDrive<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConnectionStep, ChannelError> {
        poll(cancelled)?;
        self.validate_binding(first, second)?;
        if !positive(dt) || dt > first.maximum_step_s().min(second.maximum_step_s()) {
            return Err(ChannelError::Invalid(
                "junction step must satisfy both original channel Courant limits",
            ));
        }
        drive.forces[0].validate(first.cells() + 1, &self.material, cancelled)?;
        drive.forces[1].validate(second.cells() + 1, &self.material, cancelled)?;
        let source = half_volume_force(first, second, drive.forces);
        if !source.into_iter().all(f64::is_finite) {
            return Err(ChannelError::Invalid(
                "junction body force is not representable",
            ));
        }
        let source_correction = self.correction(0.5 * dt, [[0.0; 2]; 2], source)?;
        let mut reference = [0, 1].map(|i| self.flow_m3_s[i] + source_correction[i]);
        if self.kind == ConnectionKind::Superleak {
            reference[0] = -self.material.superfluid_fraction() * reference[1];
        }
        let outer = [
            outward_flow(first, drive.outer.left, -1.0),
            outward_flow(second, drive.outer.right, 1.0),
        ];
        let reference_ports = [[outer[0], reference], [reference.map(|q| -q), outer[1]]];
        let a = first.thermomechanical_ends_response(
            dt,
            reference_ports[0],
            drive.forces[0],
            drive.friction[0],
            &mut *cancelled,
        )?;
        let b = second.thermomechanical_ends_response(
            dt,
            reference_ports[1],
            drive.forces[1],
            drive.friction[1],
            &mut *cancelled,
        )?;
        let difference = [0, 1].map(|i| a.reference_effort_pa[1][i] - b.reference_effort_pa[0][i]);
        let resistance = [0, 1]
            .map(|i| [0, 1].map(|j| a.impedance_pa_s_m3[2 + i][2 + j] + b.impedance_pa_s_m3[i][j]));
        let correction = self.correction(0.5 * dt, resistance, difference)?;
        let mut midpoint = [0, 1].map(|i| reference[i] + correction[i]);
        if self.kind == ConnectionKind::Superleak {
            midpoint[0] = -self.material.superfluid_fraction() * midpoint[1];
        }
        let mut after = [0, 1].map(|i| self.flow_m3_s[i] + 2.0 * (midpoint[i] - self.flow_m3_s[i]));
        if self.kind == ConnectionKind::Superleak {
            // Store one constrained state coordinate, preventing roundoff drift
            // into normal current; this does not project an incompatible input.
            after[0] = -self.material.superfluid_fraction() * after[1];
        }
        for (i, channel) in [&*first, &*second].into_iter().enumerate() {
            let v = after[0] / channel.area_m2;
            let w = after[1] / channel.area_m2;
            let speeds = [
                v + self.material.superfluid_fraction() * w,
                v - self.material.normal_fraction() * w,
            ];
            if speeds
                .into_iter()
                .any(|speed| !speed.is_finite() || speed.abs() > channel.profile_speed_limit())
            {
                return Err(ChannelError::Domain {
                    field: "junction endpoint velocity",
                    index: i,
                });
            }
        }
        let ports = [
            first.thermomechanical_boundary_velocities([outer[0], midpoint])?,
            second.thermomechanical_boundary_velocities([midpoint.map(|q| -q), outer[1]])?,
        ];
        let first_report =
            first.stage_step(dt, ports[0], drive.friction[0], drive.forces[0], cancelled)?;
        let second_report =
            second.stage_step(dt, ports[1], drive.friction[1], drive.forces[1], cancelled)?;
        let ea = first.thermomechanical_midpoint_efforts()?;
        let eb = second.thermomechanical_midpoint_efforts()?;
        let work = dt * dot([ea[1][0] - eb[0][0], ea[1][1] - eb[0][1]], midpoint);
        let body_work = dt * dot(source, midpoint);
        let outer_work = -dt * (dot(ea[0], outer[0]) + dot(eb[1], outer[1]));
        let before_energy = self.kinetic_energy_j()?;
        let after_energy = kinetic(self.inertance, after)?;
        // Resolve a small increment without subtracting large kinetic stores.
        let change = sum([0, 1].map(|i| {
            (self.inertance[i] * (after[i] - self.flow_m3_s[i]))
                * (0.5 * after[i] + 0.5 * self.flow_m3_s[i])
        }));
        let junction_defect = sum([change, -work, -body_work]);
        let port_defect = sum([
            first_report.boundary_wave_work_j,
            second_report.boundary_wave_work_j,
            work,
            -outer_work,
        ]);
        let combined = sum([
            first_report.wave_balance_defect_j,
            second_report.wave_balance_defect_j,
            junction_defect,
            port_defect,
        ]);
        let scale = before_energy
            + after_energy
            + work.abs()
            + body_work.abs()
            + outer_work.abs()
            + first_report.wave_energy_before_j
            + first_report.wave_energy_after_j
            + second_report.wave_energy_before_j
            + second_report.wave_energy_after_j
            + first_report.dissipated_heat_j()
            + second_report.dissipated_heat_j()
            + first_report.vortex_energy_before_j
            + first_report.vortex_energy_after_j
            + second_report.vortex_energy_before_j
            + second_report.vortex_energy_after_j;
        let tolerance = 131072.0 * f64::EPSILON * scale;
        if ![
            work,
            body_work,
            outer_work,
            change,
            junction_defect,
            port_defect,
            combined,
            tolerance,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(ChannelError::Invalid(
                "junction work balance is not representable",
            ));
        }
        for defect in [junction_defect, port_defect, combined] {
            if defect.abs() > tolerance {
                return Err(ChannelError::Balance {
                    defect_j: defect,
                    tolerance_j: tolerance,
                });
            }
        }
        let rho = self.material.density_kg_m3();
        let s = self.material.entropy_j_kg_k();
        let fs = self.material.superfluid_fraction();
        let external_mass = -dt * rho * (outer[0][0] + outer[1][0]);
        let external_entropy =
            -dt * rho * s * (outer[0][0] + fs * outer[0][1] + outer[1][0] + fs * outer[1][1]);
        let mass_defect = sum([
            first_report.boundary_mass_kg,
            second_report.boundary_mass_kg,
            -external_mass,
        ]);
        let entropy_defect = sum([
            first_report.boundary_entropy_j_k,
            second_report.boundary_entropy_j_k,
            -external_entropy,
        ]);
        // Include both terms in the scale when the normal entropy current
        // vanishes through cancellation of bulk and relative currents.
        let current_scale = outer
            .into_iter()
            .flatten()
            .chain(midpoint)
            .map(f64::abs)
            .sum::<f64>();
        let mass_scale = dt * rho * current_scale;
        let entropy_scale = mass_scale * s;
        if ![mass_defect, entropy_defect, mass_scale, entropy_scale]
            .into_iter()
            .all(f64::is_finite)
            || mass_defect.abs() > 8192.0 * f64::EPSILON * mass_scale
            || entropy_defect.abs() > 8192.0 * f64::EPSILON * entropy_scale
        {
            return Err(ChannelError::Invalid(
                "junction mass or entropy transfer is unresolved",
            ));
        }
        poll(cancelled)?;
        Ok(ConnectionStep {
            channels: [first_report, second_report],
            flow_before_m3_s: self.flow_m3_s,
            flow_midpoint_m3_s: midpoint,
            flow_after_m3_s: after,
            junction_energy_before_j: before_energy,
            junction_energy_after_j: after_energy,
            junction_wave_work_j: work,
            junction_body_work_j: body_work,
            outer_wave_work_j: outer_work,
            junction_balance_defect_j: junction_defect,
            combined_energy_defect_j: combined,
            mass_balance_defect_kg: mass_defect,
            entropy_balance_defect_j_k: entropy_defect,
        })
    }

    fn correction(
        &self,
        half_dt: f64,
        resistance: [[f64; 2]; 2],
        effort: [f64; 2],
    ) -> Result<[f64; 2], ChannelError> {
        let answer = if self.kind == ConnectionKind::Superleak {
            let direction = [self.material.superfluid_fraction(), -1.0];
            let inertia = self.inertance[0] * direction[0] * direction[0] + self.inertance[1];
            let load = dot(direction, effort);
            let response = sum([
                direction[0] * direction[0] * resistance[0][0],
                -2.0 * direction[0] * resistance[0][1],
                resistance[1][1],
            ]);
            let denominator = inertia + half_dt * response;
            if !positive(denominator) || !response.is_finite() || response < 0.0 {
                return Err(ChannelError::Invalid(
                    "superleak inertial response is unresolved",
                ));
            }
            let change = half_dt * load / denominator;
            direction.map(|d| d * change)
        } else {
            // Inertia-normalized Cholesky preserves the symmetric response.
            let root = self.inertance.map(f64::sqrt);
            let diagonal = [0, 1].map(|i| 1.0 + half_dt * (resistance[i][i] / self.inertance[i]));
            let cross = half_dt * (resistance[0][1] / root[0]) / root[1];
            let rhs = [0, 1].map(|i| half_dt * (effort[i] / root[i]));
            let l00 = diagonal[0].sqrt();
            let l10 = cross / l00;
            let schur = diagonal[1] - l10 * l10;
            if !positive(l00) || !positive(schur) {
                return Err(ChannelError::Invalid(
                    "junction inertial response is unresolved",
                ));
            }
            let y0 = rhs[0] / l00;
            let x1 = (rhs[1] - l10 * y0) / schur;
            [(y0 - l10 * x1) / l00 / root[0], x1 / root[1]]
        };
        if answer.into_iter().all(f64::is_finite) {
            Ok(answer)
        } else {
            Err(ChannelError::Invalid(
                "junction flow correction is not representable",
            ))
        }
    }
}

fn endpoint_flow(channel: &TwoFluidChannel, face: usize) -> Result<[f64; 2], ChannelError> {
    let sample = channel
        .face(face)
        .ok_or(ChannelError::Invalid("missing junction endpoint"))?;
    let flow = [
        sample.barycentric_velocity_m_s,
        sample.counterflow_velocity_m_s,
    ]
    .map(|v| channel.area_m2 * v);
    if flow.into_iter().all(f64::is_finite) {
        Ok(flow)
    } else {
        Err(ChannelError::Invalid(
            "junction endpoint current is not representable",
        ))
    }
}

fn outward_flow(channel: &TwoFluidChannel, end: EndVelocities, sign: f64) -> [f64; 2] {
    [
        channel.material.normal_fraction() * end.normal_m_s
            + channel.material.superfluid_fraction() * end.superfluid_m_s,
        end.normal_m_s - end.superfluid_m_s,
    ]
    .map(|v| sign * channel.area_m2 * v)
}

fn half_volume_force(
    first: &TwoFluidChannel,
    second: &TwoFluidChannel,
    force: [BodyForce<'_>; 2],
) -> [f64; 2] {
    let a = force[0].components(first.cells(), &first.material);
    let b = force[1].components(0, &second.material);
    let rho = first.material.density_kg_m3();
    let reduced = rho * first.material.superfluid_fraction() * first.material.normal_fraction();
    [
        rho * (0.5 * first.dx * a.0 + 0.5 * second.dx * b.0),
        reduced * (0.5 * first.dx * a.1 + 0.5 * second.dx * b.1),
    ]
}

fn same_flow(a: [f64; 2], b: [f64; 2]) -> bool {
    let scale = a.into_iter().chain(b).map(f64::abs).fold(0.0, f64::max);
    [0, 1]
        .into_iter()
        .all(|i| (a[i] - b[i]).abs() <= 1024.0 * f64::EPSILON * scale)
}

fn kinetic(inertance: [f64; 2], flow: [f64; 2]) -> Result<f64, ChannelError> {
    let energy = sum([0, 1].map(|i| (0.5 * inertance[i] * flow[i]) * flow[i]));
    if energy.is_finite() && energy >= 0.0 {
        Ok(energy)
    } else {
        Err(ChannelError::Invalid(
            "junction kinetic storage is not representable",
        ))
    }
}

fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

fn sum<const N: usize>(values: [f64; N]) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |total, value| total + Dd::from_f64(value))
        .to_f64()
}
