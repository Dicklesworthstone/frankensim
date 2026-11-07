//! Adaptive pressure and analytic Euler velocity of the original screen jump.
use super::*;
use crate::near_field::{
    self,
    linear::open::{OpenRows, required_row_bytes},
};

/// Reusable pressure/gradient rows bound to one immutable, real-frequency
/// screen solution. Both sides and the reactive near field use the original
/// P1 double layer; no far-field or plane-wave velocity approximation is made.
pub struct PreparedScreenReceivers<'b, 'p, 's> {
    source: &'b JumpSolution<'p, 's>,
    points: Vec<Point>,
    rows: OpenRows,
}

/// The same receiver rows for a complete velocity-port bank. Rates combine
/// the original solved pressure jumps before observation, preserving coherent
/// interference across every source and both sides of the screen.
pub struct PreparedScreenPortReceivers<'b, 'p, 's> {
    source: &'b ScreenVelocityPortSolution<'p, 's>,
    receivers: PreparedScreenReceivers<'b, 'p, 's>,
}

impl<'p, 's> JumpSolution<'p, 's> {
    /// Prepare pressure and analytic spatial-gradient quadrature once for all
    /// subsequent observations. Each gradient component has its own adaptive
    /// acceptance gate; pressure convergence alone does not admit velocity.
    /// Points on either face or an edge refuse at the positive `separation`
    /// in metres. The original mesh, frequency and material remain borrowed.
    ///
    /// `maximum_row_bytes` covers retained moments, receiver positions and the
    /// returned owner. It excludes triangle/normal preparation scratch
    /// (96 bytes/facet), the original solution, the 96-byte/receiver returned
    /// field payload, allocator metadata and recursive integration stacks.
    /// Quadrature discrepancies describe observation of the solved jump, not
    /// boundary-mesh or material error. Complex-frequency solutions refuse.
    pub fn prepare_receivers<'b>(
        &'b self,
        points: &[Point],
        separation: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedScreenReceivers<'b, 'p, 's>, HelmholtzError> {
        poll(&mut cancelled)?;
        options.validate()?;
        if !self.real_frequency {
            return Err(bad("open adaptive receivers require real frequency"));
        }
        let panels = self.prepared.facets.len();
        let additional_bytes = points
            .len()
            .checked_mul(core::mem::size_of::<Point>())
            .and_then(|bytes| {
                bytes.checked_add(
                    core::mem::size_of::<PreparedScreenReceivers<'_, '_, '_>>()
                        - core::mem::size_of::<OpenRows>(),
                )
            })
            .ok_or_else(|| bad("open screen receiver payload overflows"))?;
        let required = required_row_bytes(panels, points.len())
            .and_then(|bytes| bytes.checked_add(additional_bytes));
        if points.len() > near_field::MAX_RECEIVERS
            || required.is_none_or(|bytes| bytes > maximum_row_bytes)
            || panels
                .checked_mul(points.len())
                .and_then(|count| count.checked_mul(80))
                .is_none_or(|count| count > options.maximum_kernel_evaluations)
        {
            return Err(bad("open screen receivers exceed their row or work limit"));
        }
        let mut triangles = Vec::new();
        let mut normals = Vec::new();
        triangles
            .try_reserve_exact(panels)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "open screen receiver triangles",
            })?;
        normals
            .try_reserve_exact(panels)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "open screen receiver normals",
            })?;
        for (panel, facet) in self.prepared.facets.iter().enumerate() {
            if panel % 32 == 0 {
                poll(&mut cancelled)?;
            }
            triangles.push(facet.triangle);
            normals.push(facet.normal);
        }
        let rows = OpenRows::new(
            &triangles,
            &normals,
            DampedWaveNumber {
                k: self.k,
                phase_origin_m: 0.0,
            },
            self.k.re * self.medium.sound_speed * self.medium.density,
            points,
            separation,
            options,
            maximum_row_bytes - additional_bytes,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedScreenReceivers {
            source: self,
            points: points.to_vec(),
            rows,
        })
    }
}

impl<'p, 's> ScreenVelocityPortSolution<'p, 's> {
    /// Reuse one set of adaptive pressure/gradient rows for every coherent
    /// source-rate vector at this frequency. The limits and exclusions are
    /// those of [`JumpSolution::prepare_receivers`], with the returned port
    /// owner also charged to `maximum_row_bytes`. Evaluation additionally
    /// allocates one coherent nodal jump, as in [`Self::pressure_jump`].
    pub fn prepare_receivers<'b>(
        &'b self,
        points: &[Point],
        separation: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedScreenPortReceivers<'b, 'p, 's>, HelmholtzError> {
        poll(&mut cancelled)?;
        let additional_bytes = core::mem::size_of::<PreparedScreenPortReceivers<'_, '_, '_>>()
            - core::mem::size_of::<PreparedScreenReceivers<'_, '_, '_>>();
        let inner_bytes = maximum_row_bytes
            .checked_sub(additional_bytes)
            .ok_or_else(|| bad("open screen port receiver payload exceeds its limit"))?;
        let receivers = self.solutions()[0].prepare_receivers(
            points,
            separation,
            options,
            inner_bytes,
            &mut cancelled,
        )?;
        Ok(PreparedScreenPortReceivers {
            source: self,
            receivers,
        })
    }
}

impl PreparedScreenReceivers<'_, '_, '_> {
    /// Original receiver coordinates in metres, in their original order.
    pub fn points_m(&self) -> &[Point] {
        &self.points
    }

    /// Preparation work only; evaluation repeats no kernels or boundary solve.
    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }

    /// Pressure [Pa] and Cartesian fluid particle velocity [m/s], with
    /// separate propagated quadrature discrepancies. For `exp(-i omega t)`,
    /// Euler momentum gives `v = grad(p)/(i rho omega)`, with `omega = k*c`.
    /// These fields can be passed directly to [`near_field::FirstOrder`].
    pub fn evaluate(
        &self,
        cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        self.evaluate_jump(&self.source.jump, cancelled)
    }

    fn evaluate_jump(
        &self,
        jump: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        if jump.len() != self.source.prepared.surface.vertices.len() {
            return Err(bad(
                "open screen receiver jump extent differs from its source",
            ));
        }
        self.rows.evaluate(
            |panel| self.source.prepared.facets[panel].ids.map(|v| jump[v]),
            &mut cancelled,
        )
    }
}

impl PreparedScreenPortReceivers<'_, '_, '_> {
    /// Original receiver coordinates in metres, in their original order.
    pub fn points_m(&self) -> &[Point] {
        self.receivers.points_m()
    }

    /// Preparation work only; changing rates repeats no quadrature or solve.
    pub fn kernel_evaluations(&self) -> usize {
        self.receivers.kernel_evaluations()
    }

    /// Observe the coherent field of every supplied source rate, in the
    /// original port order and units. Complex amplitudes combine before the
    /// pressure and Euler-velocity contractions, never as independent powers.
    pub fn evaluate(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        let jump = self.source.pressure_jump(rates, &mut cancelled)?;
        self.receivers.evaluate_jump(&jump, cancelled)
    }
}
