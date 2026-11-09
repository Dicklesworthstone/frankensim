//! Landau--Khalatnikov bulk transport on the original channel fields.
//!
//! Kinetic whitening diagonalizes the symmetric two-component viscous form.
//! Fourier conduction is rank one in thermodynamic energy coordinates. Each
//! scalar diffusion solve is implicit midpoint, with a bounded diffusion number
//! to avoid the unresolved sign-alternating stiff tail of Crank--Nicolson.
//! The spatial incidence used in the solve also evaluates positive dissipation
//! and end traction. No difference of large kinetic stores defines the heat.

use super::*;
use fs_material::superfluid::transport::LongitudinalTransport;

/// Accuracy/work budget for stationary closed, thermally insulated transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClosedTransportControls {
    /// Maximum `diffusivity * substep / dx²`, in (0, 0.5]. This bound keeps
    /// every resolved scalar diffusion mode nonoscillatory under midpoint.
    pub max_diffusive_courant: f64,
    /// Maximum number of diffusion substeps in EACH half of the wave split.
    pub max_substeps: usize,
}

impl Default for ClosedTransportControls {
    fn default() -> Self {
        Self {
            max_diffusive_courant: 0.5,
            max_substeps: 1024,
        }
    }
}

#[derive(Debug)]
struct DiffusionWorkspace {
    old: Vec<f64>,
    next: Vec<f64>,
    solve: SolveWorkspace,
    insulated: bool,
}

impl DiffusionWorkspace {
    fn new(n: usize, insulated: bool) -> Result<Self, ChannelError> {
        Ok(Self {
            old: filled(n)?,
            next: filled(n)?,
            solve: SolveWorkspace::new(n)?,
            insulated,
        })
    }

    fn advance(
        &mut self,
        r: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, ChannelError> {
        let n = self.old.len();
        let work = &mut self.solve;
        if work.factor_step != Some((r, r)) {
            work.factor_step = None;
            for i in 0..n {
                if i % 128 == 0 {
                    poll(cancelled)?;
                }
                let neighbors = if self.insulated {
                    usize::from(i > 0) + usize::from(i + 1 < n)
                } else {
                    2
                };
                work.diagonal[i] = 1.0 + neighbors as f64 * r;
                if i + 1 < n {
                    work.off_diagonal[i] = -r;
                }
            }
            work.factor.factor(
                &work.off_diagonal,
                &work.diagonal,
                &work.off_diagonal,
                &mut *cancelled,
            )?;
            work.factor_step = Some((r, r));
        }
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let center = self.old[i];
            let left = if i > 0 {
                self.old[i - 1] - center
            } else if self.insulated {
                0.0
            } else {
                -center
            };
            let right = if i + 1 < n {
                self.old[i + 1] - center
            } else if self.insulated {
                0.0
            } else {
                -center
            };
            work.rhs[i] = center + r * (left + right);
        }
        work.factor
            .solve_factored_into(&work.rhs, &mut self.next, &mut *cancelled)?;
        let mut error: f64 = 0.0;
        let mut scale: f64 = 0.0;
        for i in 0..n {
            if i % 128 == 0 {
                poll(cancelled)?;
            }
            let mut value = work.diagonal[i] * self.next[i];
            let mut norm = value.abs() + work.rhs[i].abs();
            if i > 0 {
                value -= r * self.next[i - 1];
                norm += (r * self.next[i - 1]).abs();
            }
            if i + 1 < n {
                value -= r * self.next[i + 1];
                norm += (r * self.next[i + 1]).abs();
            }
            error = error.max((value - work.rhs[i]).abs());
            scale = scale.max(norm);
        }
        let relative = if scale == 0.0 { error } else { error / scale };
        if !relative.is_finite() || relative > 4096.0 * f64::EPSILON {
            return Err(ChannelError::Invalid(
                "longitudinal diffusion solve has an unresolved residual",
            ));
        }
        Ok(relative)
    }

    fn midpoint(&self, i: usize) -> f64 {
        0.5 * self.old[i] + 0.5 * self.next[i]
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct TransportChange {
    pub conductive_loss_j: f64,
    pub boundary_impulse_n_s: [f64; 2],
    pub substeps: usize,
    pub residual: f64,
}

#[derive(Debug)]
pub(super) struct TransportMemory {
    model: LongitudinalTransport,
    controls: ClosedTransportControls,
    rows: [[f64; 2]; 2],
    diffusivity: [f64; 2],
    momentum_scale: [f64; 2],
    thermal_row: [f64; 2],
    thermal_diffusivity: f64,
    velocity: [DiffusionWorkspace; 2],
    thermal: DiffusionWorkspace,
    pub wave_start: [Vec<f64>; 2],
    pub pending_heat_j: Vec<f64>,
    pub last_heat_j: Vec<f64>,
    total_heat_j: f64,
    total_conductive_loss_j: f64,
}

impl TransportMemory {
    fn new(
        model: LongitudinalTransport,
        controls: ClosedTransportControls,
        modes: Modes,
        n: usize,
    ) -> Result<Self, ChannelError> {
        let state = model.material();
        let matrix = model.kinetic_diffusivity_m2_s();
        let scale = matrix[0][0].max(matrix[1][1]).max(matrix[0][1].abs());
        let mut q = [1.0, 0.0];
        if scale > 0.0 {
            let difference = matrix[0][0] / scale - matrix[1][1] / scale;
            let off = matrix[0][1] / scale;
            q = if off == 0.0 {
                if difference >= 0.0 {
                    [1.0, 0.0]
                } else {
                    [0.0, 1.0]
                }
            } else {
                let gap = difference.hypot(2.0 * off);
                if difference >= 0.0 {
                    let c = (0.5 * (1.0 + difference / gap)).sqrt();
                    [c, off / (gap * c)]
                } else {
                    let s = (0.5 * (1.0 - difference / gap)).sqrt().copysign(off);
                    [off / (gap * s), s]
                }
            };
            let norm = q[0].hypot(q[1]);
            q = q.map(|x| x / norm);
        }
        let kinetic_rows = [q, [-q[1], q[0]]];
        let rho_root = state.density_kg_m3().sqrt();
        let relative_root = state.reduced_density_kg_m3().sqrt();
        let rows = kinetic_rows.map(|q| modes.face_to_modes(q[0] / rho_root, q[1] / relative_root));
        let mut diffusivity = [0.0; 2];
        for i in 0..2 {
            let v = kinetic_rows[i][0] / rho_root;
            let w = kinetic_rows[i][1] / relative_root;
            // Evaluate the material's positive squares instead of subtracting
            // nearly equal eigenvalue products near a rank-one viscosity.
            diffusivity[i] = model
                .response(
                    v + state.superfluid_fraction() * w,
                    -state.superfluid_density_kg_m3() * w,
                    0.0,
                )?
                .viscous_heat_w_m3;
        }
        let momentum_scale =
            rows.map(|row| state.density_kg_m3() * modes.face(state, row).barycentric_velocity_m_s);
        let temperature = [
            modes.cell(state, [1.0, 0.0]).temperature_perturbation_k,
            modes.cell(state, [0.0, 1.0]).temperature_perturbation_k,
        ];
        let norm = temperature[0].hypot(temperature[1]);
        let thermal_row = temperature.map(|v| v / norm);
        let thermal_diffusivity = model.thermal_diffusivity_m2_s();
        if !positive(norm)
            || rows.into_iter().flatten().any(|v| !v.is_finite())
            || momentum_scale.into_iter().any(|v| !v.is_finite())
            || diffusivity
                .into_iter()
                .chain([thermal_diffusivity])
                .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(ChannelError::Invalid(
                "bulk transport modes are not representable",
            ));
        }
        Ok(Self {
            model,
            controls,
            rows,
            diffusivity,
            momentum_scale,
            thermal_row,
            thermal_diffusivity,
            velocity: [
                DiffusionWorkspace::new(n - 1, false)?,
                DiffusionWorkspace::new(n - 1, false)?,
            ],
            thermal: DiffusionWorkspace::new(n, true)?,
            wave_start: [filled(n)?, filled(n)?],
            pending_heat_j: filled(n)?,
            last_heat_j: filled(n)?,
            total_heat_j: 0.0,
            total_conductive_loss_j: 0.0,
        })
    }

    pub fn publish(&mut self, report: ChannelStep) {
        core::mem::swap(&mut self.last_heat_j, &mut self.pending_heat_j);
        self.total_heat_j += report.longitudinal_viscous_heat_j;
        self.total_conductive_loss_j += report.conductive_availability_loss_j;
    }

    pub fn admit_totals(&self, heat: f64, loss: f64) -> Result<(), ChannelError> {
        if [
            heat,
            loss,
            self.total_heat_j + heat,
            self.total_conductive_loss_j + loss,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v >= 0.0)
        {
            Ok(())
        } else {
            Err(ChannelError::Invalid(
                "bulk transport account is not representable",
            ))
        }
    }

    fn advance(
        &mut self,
        fields: &mut [Field; 2],
        dt: f64,
        dx: f64,
        area: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportChange, ChannelError> {
        let max_rate = self.diffusivity[0]
            .max(self.diffusivity[1])
            .max(self.thermal_diffusivity);
        if max_rate == 0.0 {
            return Ok(TransportChange::default());
        }
        let diffusion_number = (max_rate / dx) * (dt / dx);
        let count = (diffusion_number / self.controls.max_diffusive_courant)
            .ceil()
            .max(1.0);
        if !diffusion_number.is_finite() || count > self.controls.max_substeps as f64 {
            return Err(ChannelError::Invalid(
                "bulk transport diffusion substep budget exceeded",
            ));
        }
        let count = count as usize;
        let h = dt / count as f64;
        if !positive(h) {
            return Err(ChannelError::Invalid(
                "bulk transport substep is not representable",
            ));
        }
        let n = fields[0].cells.len();
        let mut change = TransportChange {
            substeps: count,
            ..TransportChange::default()
        };
        let mut loss_correction = 0.0;
        let mut impulse_correction = [0.0; 2];
        for _ in 0..count {
            poll(cancelled)?;
            for mode in 0..2 {
                let rate = self.diffusivity[mode];
                if rate == 0.0 {
                    continue;
                }
                let row = self.rows[mode];
                let work = &mut self.velocity[mode];
                for i in 1..n {
                    if i % 128 == 0 {
                        poll(cancelled)?;
                    }
                    work.old[i - 1] = row[0] * fields[0].faces[i] + row[1] * fields[1].faces[i];
                }
                let r = 0.5 * (rate / dx) * (h / dx);
                change.residual = change.residual.max(work.advance(r, cancelled)?);
                for i in 0..n {
                    if i % 128 == 0 {
                        poll(cancelled)?;
                    }
                    let left = if i == 0 { 0.0 } else { work.midpoint(i - 1) };
                    let right = if i + 1 == n { 0.0 } else { work.midpoint(i) };
                    let difference = right - left;
                    let heat = ((h * rate / dx) * difference) * difference * area;
                    if !heat.is_finite() || heat < 0.0 {
                        return Err(ChannelError::Invalid(
                            "longitudinal viscous heat is not representable",
                        ));
                    }
                    self.pending_heat_j[i] += heat;
                    if i == 0 || i + 1 == n {
                        let end = usize::from(i != 0);
                        let sign = if end == 0 { -1.0 } else { 1.0 };
                        let impulse =
                            sign * area * (h * rate / dx) * difference * self.momentum_scale[mode];
                        accumulate(
                            &mut change.boundary_impulse_n_s[end],
                            &mut impulse_correction[end],
                            impulse,
                        );
                    }
                }
                // An incremental update retains the nullspace and gives an
                // exact no-op when the represented modal change rounds to zero.
                for i in 1..n {
                    if i % 128 == 0 {
                        poll(cancelled)?;
                    }
                    let increment = work.next[i - 1] - work.old[i - 1];
                    fields[0].faces[i] += row[0] * increment;
                    fields[1].faces[i] += row[1] * increment;
                }
            }
            if self.thermal_diffusivity > 0.0 {
                let row = self.thermal_row;
                let work = &mut self.thermal;
                for i in 0..n {
                    if i % 128 == 0 {
                        poll(cancelled)?;
                    }
                    work.old[i] = row[0] * fields[0].cells[i] + row[1] * fields[1].cells[i];
                }
                let r = 0.5 * (self.thermal_diffusivity / dx) * (h / dx);
                change.residual = change.residual.max(work.advance(r, cancelled)?);
                for i in 0..n {
                    if i % 128 == 0 {
                        poll(cancelled)?;
                    }
                    if i > 0 {
                        let difference = work.midpoint(i) - work.midpoint(i - 1);
                        let loss =
                            ((h * self.thermal_diffusivity / dx) * difference) * difference * area;
                        if !loss.is_finite() || loss < 0.0 {
                            return Err(ChannelError::Invalid(
                                "conductive availability loss is not representable",
                            ));
                        }
                        accumulate(&mut change.conductive_loss_j, &mut loss_correction, loss);
                    }
                    let increment = work.next[i] - work.old[i];
                    fields[0].cells[i] += row[0] * increment;
                    fields[1].cells[i] += row[1] * increment;
                }
            }
        }
        if change
            .boundary_impulse_n_s
            .into_iter()
            .any(|v| !v.is_finite())
            || !change.conductive_loss_j.is_finite()
        {
            return Err(ChannelError::Invalid(
                "bulk transport transfer is not representable",
            ));
        }
        Ok(change)
    }
}

impl TwoFluidChannel {
    /// Attach authored bulk transport before stepping. The source must retain
    /// exactly this material state. Both component velocities must already be
    /// zero at both ends, and every subsequent step must keep them closed.
    ///
    /// Temperature has insulated ends. Longitudinal traction is carried by the
    /// stationary supports and separately reported; it does zero boundary work.
    /// Scalar normal-wall drag, forces and local/transported Vinen tangles may
    /// coexist. Transverse normal-profile modes and moving/thermal endpoint
    /// response queries refuse because their additional transport is not modeled.
    ///
    /// The second-order symmetric transport/wave split is not an exact stiff
    /// driven-equilibrium solver. Refine the wave step as well as the declared
    /// diffusion number when splitting error matters. Configuration allocates
    /// O(cells) memory; subsequent steps reuse all factors and buffers.
    pub fn configure_closed_transport(
        &mut self,
        model: LongitudinalTransport,
        controls: ClosedTransportControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        poll(&mut cancelled)?;
        let n = self.cells();
        if self.time_s != 0.0
            || *model.material() != self.material
            || self.normal_profile.is_some()
            || !positive(controls.max_diffusive_courant)
            || controls.max_diffusive_courant > 0.5
            || controls.max_substeps == 0
            || self
                .accepted
                .iter()
                .any(|f| f.faces[0] != 0.0 || f.faces[n] != 0.0)
        {
            return Err(ChannelError::Invalid(
                "bulk transport requires the same material, stationary closed ends, no transverse profile and valid controls before stepping",
            ));
        }
        let memory = TransportMemory::new(model, controls, self.modes, n)?;
        poll(&mut cancelled)?;
        self.transport = Some(memory);
        Ok(())
    }

    /// Original authored bulk-transport material and coefficients, if enabled.
    pub fn longitudinal_transport(&self) -> Option<&LongitudinalTransport> {
        self.transport.as_ref().map(|t| &t.model)
    }

    /// Physical bulk viscous heat from the last accepted interval, by cell.
    pub fn last_longitudinal_viscous_heat_j(&self) -> Option<&[f64]> {
        self.transport.as_ref().map(|t| t.last_heat_j.as_slice())
    }

    /// Cumulative physical bulk viscous heat [J], excluding conduction.
    pub fn cumulative_longitudinal_viscous_heat_j(&self) -> f64 {
        self.transport.as_ref().map_or(0.0, |t| t.total_heat_j)
    }

    /// Cumulative loss of wave availability from Fourier conduction [J].
    /// Divide by the frozen reference temperature for produced entropy [J/K].
    pub fn cumulative_conductive_availability_loss_j(&self) -> f64 {
        self.transport
            .as_ref()
            .map_or(0.0, |t| t.total_conductive_loss_j)
    }

    pub(super) fn transport_half_step(
        &mut self,
        dt: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportChange, ChannelError> {
        match &mut self.transport {
            Some(t) => t.advance(&mut self.pending, dt, self.dx, self.area_m2, cancelled),
            None => Ok(TransportChange::default()),
        }
    }

    pub(super) fn wave_initial_cells(&self, mode: usize) -> &[f64] {
        match &self.transport {
            Some(t) => &t.wave_start[mode],
            None => &self.accepted[mode].cells,
        }
    }
}
