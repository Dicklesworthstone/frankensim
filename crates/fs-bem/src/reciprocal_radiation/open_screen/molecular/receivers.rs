//! Adaptive pressure and Euler velocity of the original coherent jump field.
use super::*;
use crate::near_field::{
    self,
    linear::open::{OpenRows, required_row_bytes},
};

/// Reusable pressure/gradient rows bound by borrow to one immutable molecular
/// screen solution. Every rate vector combines its original unit-source jumps
/// before a single contraction, preserving coherent interference on both sides.
pub struct PreparedMolecularScreenReceivers<'b, 'p, 's> {
    source: &'b MolecularScreenVelocityPortSolution<'p, 's>,
    points: Vec<Point>,
    rows: OpenRows,
}

impl<'p, 's> MolecularScreenVelocityPortSolution<'p, 's> {
    /// Integrate the original retained P1 geometry once for all subsequent
    /// coherent rate vectors. Pressure and each gradient component have their
    /// own adaptive acceptance gate. A receiver on either face or an edge
    /// refuses at the declared positive separation.
    ///
    /// `maximum_row_bytes` covers the retained moments, points and owner. It
    /// excludes bounded triangle/normal preparation scratch (96 bytes/facet),
    /// the original solved source, the coherent nodal jump returned during
    /// evaluation and the 96-byte/receiver observation payload, allocator
    /// metadata and recursive integration stack. Quadrature estimates describe
    /// integration of the solved jump only, not mesh or material uncertainty.
    pub fn prepare_receivers<'b>(
        &'b self,
        points: &[Point],
        separation: f64,
        options: near_field::Options,
        maximum_row_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularScreenReceivers<'b, 'p, 's>, HelmholtzError> {
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        options.validate()?;
        let panels = self.prepared.facets.len();
        let additional_bytes = points
            .len()
            .checked_mul(core::mem::size_of::<Point>())
            .and_then(|bytes| {
                bytes.checked_add(
                    core::mem::size_of::<PreparedMolecularScreenReceivers<'_, '_, '_>>()
                        - core::mem::size_of::<OpenRows>(),
                )
            })
            .ok_or_else(|| bad("molecular screen receiver payload overflows"))?;
        let required = required_row_bytes(panels, points.len())
            .and_then(|bytes| bytes.checked_add(additional_bytes));
        if points.len() > near_field::MAX_RECEIVERS
            || required.is_none_or(|bytes| bytes > maximum_row_bytes)
            || panels
                .checked_mul(points.len())
                .and_then(|count| count.checked_mul(80))
                .is_none_or(|count| count > options.maximum_kernel_evaluations)
        {
            return Err(bad(
                "molecular screen receivers exceed their row or work limit",
            ));
        }
        let mut triangles = Vec::new();
        let mut normals = Vec::new();
        triangles
            .try_reserve_exact(panels)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "molecular screen receiver triangles",
            })?;
        normals
            .try_reserve_exact(panels)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "molecular screen receiver normals",
            })?;
        for (panel, facet) in self.prepared.facets.iter().enumerate() {
            if panel % 32 == 0 && cancelled() {
                return Err(HelmholtzError::Cancelled);
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
            self.source.source().gas().density * self.omega,
            points,
            separation,
            options,
            maximum_row_bytes - additional_bytes,
            &mut cancelled,
        )?;
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(PreparedMolecularScreenReceivers {
            source: self,
            points: points.to_vec(),
            rows,
        })
    }
}

impl PreparedMolecularScreenReceivers<'_, '_, '_> {
    /// Original receiver coordinates in metres, in their original order.
    pub fn points_m(&self) -> &[Point] {
        &self.points
    }

    /// Preparation work only. Applying another rate vector repeats no kernel
    /// integration or boundary solve.
    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }

    /// Coherent near pressure [Pa] and Euler particle velocity [m/s] in the
    /// original frame. The same real material density and physical omega
    /// divide the analytic gradient; complex `k*c_eq` is never substituted.
    pub fn evaluate(
        &self,
        rates: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        let jump = self.source.pressure_jump(rates, &mut cancelled)?;
        self.evaluate_jump(&jump, cancelled)
    }

    /// Internal matched jump observation for an incident-plus-motion owner.
    /// The immutable geometry/wave rows remain those of this source.
    pub(super) fn evaluate_jump(
        &self,
        jump: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<near_field::VectorObservation, HelmholtzError> {
        if jump.len() != self.source.prepared.surface.vertices.len() {
            return Err(bad(
                "molecular screen receiver jump extent differs from its source",
            ));
        }
        self.rows.evaluate(
            |panel| {
                self.source.prepared.facets[panel]
                    .ids
                    .map(|vertex| jump[vertex])
            },
            &mut cancelled,
        )
    }
}
