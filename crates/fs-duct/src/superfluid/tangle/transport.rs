//! Conservative transport of unresolved vortex length on interior dual cells.
//!
//! The authored law is `L_t + (u_L L)_x = 0`; `u_L` is independent of both
//! component velocities. Limited linear upwinding and SSPRK2 preserve positivity
//! under an outgoing Courant bound of 1/2. One interface transfer is shared by
//! its two neighbours. No clipping, diffusion coefficient, stretching, or
//! filament geometry is inferred. Constant effective line energy travels with L.

use super::{ChannelError, TwoFluidChannel, VinenTangle, accumulate, filled, poll};

/// Accuracy and work controls, independent of the local Vinen source budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VortexTransportControls {
    /// Maximum sum of the outgoing Courant numbers of a cell, in `(0, 0.5]`.
    /// Refinement controls temporal error; this is not an a posteriori bound.
    pub maximum_outgoing_courant: f64,
    /// Maximum SSPRK2 substeps per transport HALF-step, in `1..=4096`.
    pub maximum_substeps: usize,
}

/// Authored condition at one end of the interior dual-cell transport domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum VortexTransportBoundary {
    /// Zero line flux. The supplied velocity at this interface must be zero.
    Reflecting,
    /// This density is used only on inflow. Outflow uses the adjacent channel
    /// density. Transferred effective line energy is an external port, not heat.
    Reservoir { incoming_line_density_per_m2: f64 },
}

/// Left and right conditions at `dx/2` and `length-dx/2`, respectively.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VortexTransportBoundaries {
    pub left: VortexTransportBoundary,
    pub right: VortexTransportBoundary,
}

impl VortexTransportBoundaries {
    pub const CLOSED: Self = Self {
        left: VortexTransportBoundary::Reflecting,
        right: VortexTransportBoundary::Reflecting,
    };
}

#[derive(Debug)]
pub(super) struct VortexTransportMemory {
    velocities: Vec<f64>,
    boundaries: VortexTransportBoundaries,
    controls: VortexTransportControls,
    maximum_outgoing_speed: f64,
    active: bool,
    stage: Vec<f64>,
    slopes: Vec<f64>,
    transfer: Vec<f64>,
}

#[derive(Clone, Copy, Debug, Default)]
pub(in crate::superfluid) struct TransportChange {
    pub boundary_line_m: [f64; 2],
    pub balance_defect_m: f64,
    pub substeps: usize,
}

impl TransportChange {
    pub(in crate::superfluid) fn add(&mut self, other: Self) -> Result<(), ChannelError> {
        for end in 0..2 {
            self.boundary_line_m[end] += other.boundary_line_m[end];
        }
        self.balance_defect_m += other.balance_defect_m;
        self.substeps += other.substeps; // At most twice the admitted 4096 budget.
        if self.boundary_line_m.into_iter().all(f64::is_finite) && self.balance_defect_m.is_finite()
        {
            Ok(())
        } else {
            Err(ChannelError::Invalid(
                "vortex transport totals are not representable",
            ))
        }
    }
}

impl VortexTransportMemory {
    fn new(
        velocities: &[f64],
        boundaries: VortexTransportBoundaries,
        controls: VortexTransportControls,
        model: &VinenTangle,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ChannelError> {
        if !controls.maximum_outgoing_courant.is_finite()
            || controls.maximum_outgoing_courant <= 0.0
            || controls.maximum_outgoing_courant > 0.5
            || !(1..=4096).contains(&controls.maximum_substeps)
        {
            return Err(ChannelError::Invalid(
                "vortex transport needs outgoing Courant in (0,0.5] and 1..=4096 substeps",
            ));
        }
        let cells = velocities.len() - 1; // Channel admission supplies N >= 2.
        let mut memory = Self {
            velocities: filled(velocities.len())?,
            boundaries,
            controls,
            maximum_outgoing_speed: 0.0,
            active: false,
            stage: filled(cells)?,
            slopes: filled(cells)?,
            transfer: filled(velocities.len())?,
        };
        memory.validate_boundary_densities(model)?;
        for (edge, &velocity) in velocities.iter().enumerate() {
            if edge % 128 == 0 {
                poll(cancelled)?;
            }
            if !velocity.is_finite() {
                return Err(ChannelError::Domain {
                    field: "vortex drift velocity",
                    index: edge,
                });
            }
            memory.velocities[edge] = velocity;
            memory.active |= velocity != 0.0;
        }
        for (boundary, edge) in [(boundaries.left, 0), (boundaries.right, cells)] {
            if boundary == VortexTransportBoundary::Reflecting && velocities[edge] != 0.0 {
                return Err(ChannelError::Invalid(
                    "reflecting vortex boundaries require zero authored endpoint drift",
                ));
            }
        }
        for cell in 0..cells {
            if cell % 128 == 0 {
                poll(cancelled)?;
            }
            let outgoing = (-velocities[cell]).max(0.0) + velocities[cell + 1].max(0.0);
            if !outgoing.is_finite() {
                return Err(ChannelError::Invalid(
                    "vortex outgoing speed is not representable",
                ));
            }
            memory.maximum_outgoing_speed = memory.maximum_outgoing_speed.max(outgoing);
        }
        Ok(memory)
    }

    pub(super) fn validate_boundary_densities(
        &self,
        model: &VinenTangle,
    ) -> Result<(), ChannelError> {
        for boundary in [self.boundaries.left, self.boundaries.right] {
            if let VortexTransportBoundary::Reservoir {
                incoming_line_density_per_m2,
            } = boundary
            {
                model.line_energy_density_j_m3(incoming_line_density_per_m2)?;
            }
        }
        Ok(())
    }

    pub(super) fn advance(
        &mut self,
        density: &mut [f64],
        dt: f64,
        dx: f64,
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportChange, ChannelError> {
        poll(cancelled)?;
        if !self.active {
            return Ok(TransportChange::default());
        }
        let courant = (dt / dx) * self.maximum_outgoing_speed;
        let required = (courant / self.controls.maximum_outgoing_courant)
            .ceil()
            .max(1.0);
        if !courant.is_finite() || required > self.controls.maximum_substeps as f64 {
            return Err(ChannelError::Invalid(
                "vortex transport substep budget exhausted",
            ));
        }
        let substeps = required as usize;
        let ratio = (dt / substeps as f64) / dx;
        if !ratio.is_finite() || ratio < 0.0 {
            return Err(ChannelError::Invalid(
                "vortex transport interval is not representable",
            ));
        }
        let before = line_length(density, volume, cancelled)?;
        let cells = density.len();
        let mut boundary_line_m = [0.0; 2];
        let mut correction = [0.0; 2];
        for _ in 0..substeps {
            poll(cancelled)?;
            transfers(
                density,
                &self.velocities,
                self.boundaries,
                ratio,
                &mut self.slopes,
                &mut self.transfer,
                cancelled,
            )?;
            let first_ports = [self.transfer[0], -self.transfer[cells]];
            for cell in 0..cells {
                if cell % 128 == 0 {
                    poll(cancelled)?;
                }
                self.stage[cell] =
                    euler_cell(density[cell], self.transfer[cell], self.transfer[cell + 1])?;
            }
            transfers(
                &self.stage,
                &self.velocities,
                self.boundaries,
                ratio,
                &mut self.slopes,
                &mut self.transfer,
                cancelled,
            )?;
            for cell in 0..cells {
                if cell % 128 == 0 {
                    poll(cancelled)?;
                }
                let second = euler_cell(
                    self.stage[cell],
                    self.transfer[cell],
                    self.transfer[cell + 1],
                )?;
                density[cell] += 0.5 * (second - density[cell]);
            }
            for (end, second) in [self.transfer[0], -self.transfer[cells]]
                .into_iter()
                .enumerate()
            {
                let line = (first_ports[end] + 0.5 * (second - first_ports[end])) * volume;
                accumulate(&mut boundary_line_m[end], &mut correction[end], line);
            }
        }
        let after = line_length(density, volume, cancelled)?;
        let defect = (after - before) - boundary_line_m[0] - boundary_line_m[1];
        let tolerance = 32768.0
            * f64::EPSILON
            * (before + after + boundary_line_m[0].abs() + boundary_line_m[1].abs());
        if ![defect, tolerance, boundary_line_m[0], boundary_line_m[1]]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(ChannelError::Invalid(
                "vortex line transfer is not representable",
            ));
        }
        if defect.abs() > tolerance {
            return Err(ChannelError::VortexTransportBalance {
                defect_m: defect,
                tolerance_m: tolerance,
            });
        }
        Ok(TransportChange {
            boundary_line_m,
            balance_defect_m: defect,
            substeps,
        })
    }
}

fn line_length(
    density: &[f64],
    volume: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, ChannelError> {
    let mut total = 0.0;
    let mut correction = 0.0;
    for (cell, &value) in density.iter().enumerate() {
        if cell % 128 == 0 {
            poll(cancelled)?;
        }
        if !value.is_finite() || value < 0.0 {
            return Err(ChannelError::Domain {
                field: "transported vortex density",
                index: cell + 1,
            });
        }
        accumulate(&mut total, &mut correction, value * volume);
    }
    if total.is_finite() && total >= 0.0 {
        Ok(total)
    } else {
        Err(ChannelError::Invalid(
            "total vortex length is not representable",
        ))
    }
}

fn incoming(boundary: VortexTransportBoundary) -> f64 {
    match boundary {
        VortexTransportBoundary::Reflecting => 0.0, // Admitted endpoint drift is zero.
        VortexTransportBoundary::Reservoir {
            incoming_line_density_per_m2,
        } => incoming_line_density_per_m2,
    }
}

fn transfers(
    density: &[f64],
    velocities: &[f64],
    boundaries: VortexTransportBoundaries,
    ratio: f64,
    slopes: &mut [f64],
    transfer: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ChannelError> {
    let cells = density.len();
    slopes[0] = 0.0;
    slopes[cells - 1] = 0.0;
    for cell in 1..cells - 1 {
        if cell % 128 == 0 {
            poll(cancelled)?;
        }
        let left = density[cell] - density[cell - 1];
        let right = density[cell + 1] - density[cell];
        slopes[cell] = if left > 0.0 && right > 0.0 {
            left.min(right)
        } else if left < 0.0 && right < 0.0 {
            left.max(right)
        } else {
            0.0
        };
    }
    for (edge, &velocity) in velocities.iter().enumerate() {
        if edge % 128 == 0 {
            poll(cancelled)?;
        }
        let donor = if velocity > 0.0 {
            if edge == 0 {
                incoming(boundaries.left)
            } else {
                density[edge - 1] + 0.5 * slopes[edge - 1]
            }
        } else if velocity < 0.0 {
            if edge == cells {
                incoming(boundaries.right)
            } else {
                density[edge] - 0.5 * slopes[edge]
            }
        } else {
            0.0
        };
        let value = (ratio * velocity) * donor;
        if !value.is_finite() || !donor.is_finite() || donor < 0.0 {
            return Err(ChannelError::Invalid(
                "vortex interface transfer is not representable",
            ));
        }
        transfer[edge] = value;
    }
    Ok(())
}

fn euler_cell(value: f64, left: f64, right: f64) -> Result<f64, ChannelError> {
    // Each reconstructed donor is at most 1.5*its cell average. The outgoing
    // Courant sum <= 0.5 therefore leaves >= 0.25*value before adding inflow.
    // Keep positive inflow separate from outflow; never clip a negative result.
    let remainder = (value - (-left).max(0.0)) - right.max(0.0);
    let next = (remainder + left.max(0.0)) + (-right).max(0.0);
    if remainder >= 0.0 && next.is_finite() && next >= 0.0 {
        Ok(next)
    } else {
        Err(ChannelError::Invalid(
            "vortex transport positivity or finite arithmetic failed",
        ))
    }
}

impl TwoFluidChannel {
    /// Configure or replace an authored drift field for subsequent intervals.
    /// Requires a configured Vinen tangle. `interface_velocity_m_s` has N values
    /// in +x orientation at `x=(i+1/2)*dx`, bounding the N-1 stored dual cells.
    /// Endpoint half-cell vortex storage belongs to the boundary owner.
    ///
    /// This law transports line length and its constant effective line energy;
    /// it is NOT inferred from normal flow, superflow or vortex self-induction.
    /// Drift and reservoir values are held fixed over each channel step.
    /// Symmetric transport/source/wave splitting preserves the affine endpoint
    /// responses. Reservoir energy enters through separate reported ports;
    /// closed transport creates neither heat nor additional stored energy.
    ///
    /// Configuration allocates reusable O(N) workspace. It can change between
    /// accepted steps without resetting density. Refusal/cancellation leaves
    /// the previous configuration and all accepted fields intact.
    pub fn configure_vortex_line_transport(
        &mut self,
        interface_velocity_m_s: &[f64],
        boundaries: VortexTransportBoundaries,
        controls: VortexTransportControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), ChannelError> {
        poll(&mut cancelled)?;
        if interface_velocity_m_s.len() != self.cells() {
            return Err(ChannelError::Invalid(
                "vortex transport requires N interface velocities",
            ));
        }
        let memory = self.vinen_tangle.as_mut().ok_or(ChannelError::Invalid(
            "vortex transport requires a configured Vinen tangle",
        ))?;
        let transport = VortexTransportMemory::new(
            interface_velocity_m_s,
            boundaries,
            controls,
            &memory.model,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        memory.transport = Some(transport);
        Ok(())
    }

    /// Authored drift [m/s] at dual interface 0..N-1, if configured.
    pub fn vortex_transport_velocity_m_s(&self, interface: usize) -> Option<f64> {
        self.vinen_tangle
            .as_ref()?
            .transport
            .as_ref()?
            .velocities
            .get(interface)
            .copied()
    }

    pub fn vortex_transport_boundaries(&self) -> Option<VortexTransportBoundaries> {
        Some(self.vinen_tangle.as_ref()?.transport.as_ref()?.boundaries)
    }

    pub fn vortex_transport_controls(&self) -> Option<VortexTransportControls> {
        Some(self.vinen_tangle.as_ref()?.transport.as_ref()?.controls)
    }
}
