//! Mean-dependent outward flux independent of any convective control port.
//! The complete temperature correction is G (dq/dmean) W. Each retained
//! left factor is the exact P1 flux injection G_i; its right factor is one
//! row of (dq/dmean) W. Full-node factors retain prescribed-value columns.

use std::collections::BTreeSet;

use crate::{ConductionMesh, ThermalBoundary};
use super::{ConductionError, Cx, add, checked, invalid, poll, vector, zeros};

#[derive(Debug, Clone)]
struct FluxFactor {
    left: Vec<f64>,
    right: Vec<f64>,
}

/// Bounded full-node feedback from uniform fluxes on exact boundary regions.
///
/// A flux q_i is positive outward, in W/m2, and may depend on every selected
/// region's area-mean temperature. The row-major derivative matrix contains
/// dq_i/dmean_j in W/(m2 K). Robin controls are selected independently: zero
/// convective regions are valid when a transient storage operator supplies the
/// temperature anchor. This type supplies no primal solve or physical-law
/// admission; its caller binds the same mesh, trace partition and state to the
/// checked endpoint and its derivatives. Geometry and this matrix stay fixed.
#[derive(Debug, Clone)]
pub struct MeanFluxFeedback {
    vertices: usize,
    regions: Vec<String>,
    factors: Vec<FluxFactor>,
    entries: usize,
}

impl MeanFluxFeedback {
    /// Build G_i from every face in the named boundary region, retaining the
    /// caller's region order. No secondary geometric selector or supplied area
    /// can change the trace. At most 64 regions and 2*nodes*regions retained
    /// scalars must fit max_feedback_entries before any factor is allocated.
    ///
    /// # Errors
    /// Refuses empty, duplicate or unknown regions, empty/nonfinite traces,
    /// malformed derivative matrices, insufficient factor storage, nonfinite
    /// arithmetic, or cancellation. A zero matrix still validates the traces.
    pub fn new(
        cx: &Cx<'_>,
        mesh: &ConductionMesh,
        boundary: &ThermalBoundary,
        regions: &[&str],
        flux_derivatives_w_m2_k: &[f64],
        max_feedback_entries: usize,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let n = mesh.vertex_count();
        let p = regions.len();
        let entries = n.checked_mul(2).and_then(|value| value.checked_mul(p));
        if n == 0 || !(1..=64).contains(&p)
            || entries.is_none_or(|value| value > max_feedback_entries)
            || entries.and_then(|value| value.checked_mul(std::mem::size_of::<f64>()))
                .is_none_or(|value| isize::try_from(value).is_err())
        {
            return Err(invalid("mean-flux feedback exceeds its two-vector allowance or 1..=64 regions"));
        }
        vector(cx, flux_derivatives_w_m2_k, p * p)?;
        let mut seen = BTreeSet::new();
        let mut selected = Vec::with_capacity(p);
        for &name in regions {
            poll(cx, selected.len())?;
            if !seen.insert(name) {
                return Err(invalid("mean-flux feedback repeats a boundary region"));
            }
            selected.push(boundary.region_names().iter().position(|found| found == name)
                .ok_or_else(|| invalid("mean-flux feedback names an unknown boundary region"))?);
        }
        let mut lefts = (0..p).map(|_| zeros(n)).collect::<Result<Vec<_>, _>>()?;
        let mut areas = vec![0.0; p];
        for (slot, face) in mesh.boundary().iter().enumerate() {
            if slot % 512 == 0 { poll(cx, slot)?; }
            let Some(region) = boundary.region_for(slot)
                .and_then(|index| selected.iter().position(|&found| found == index))
            else { continue; };
            if !face.area.is_finite() || face.area <= 0.0 {
                return Err(invalid("mean-flux feedback requires finite positive face areas"));
            }
            add(&mut areas[region], face.area)?;
            for &vertex in &face.vertices {
                add(&mut lefts[region][vertex as usize], face.area / 3.0)?;
            }
        }
        if areas.iter().any(|&area| !area.is_finite() || area <= 0.0) {
            return Err(invalid("mean-flux feedback has an empty boundary trace"));
        }
        let mut rights = Vec::with_capacity(p);
        for row in 0..p {
            poll(cx, row)?;
            let mut right = zeros(n)?;
            for column in 0..p {
                let slope = flux_derivatives_w_m2_k[row * p + column];
                if slope == 0.0 { continue; }
                for (vertex, (value, &injection)) in right.iter_mut().zip(&lefts[column]).enumerate() {
                    if vertex % 512 == 0 { poll(cx, vertex)?; }
                    let mean_weight = checked(injection / areas[column])?;
                    *value = checked(slope.mul_add(mean_weight, *value))?;
                }
            }
            rights.push(right);
        }
        let factors = lefts.into_iter().zip(rights)
            .map(|(left, right)| FluxFactor { left, right }).collect();
        poll(cx, n)?;
        Ok(Self {
            vertices: n,
            regions: regions.iter().map(|name| (*name).to_owned()).collect(),
            factors,
            entries: entries.expect("checked factor count"),
        })
    }

    /// Full reference-mesh vertex count used by every factor and contraction.
    #[must_use]
    pub const fn vertex_count(&self) -> usize { self.vertices }

    /// Number of independently named flux regions.
    #[must_use]
    pub fn region_count(&self) -> usize { self.regions.len() }

    /// Boundary names in the flux-matrix and pullback order.
    #[must_use]
    pub fn region_names(&self) -> &[String] { &self.regions }

    /// Retained scalar factor entries, excluding allocation metadata.
    #[must_use]
    pub const fn factor_entries(&self) -> usize { self.entries }

    /// Return G^T lambda in region order, with a POSITIVE flux-injection sign.
    /// An outward-flux parameter enters the residual positively, so its
    /// objective partial uses the negative of this contraction. The caller
    /// supplies its physical multiplier and timestep convention exactly once;
    /// this method neither masks prescribed rows nor introduces another dt.
    pub fn flux_pullback(
        &self, cx: &Cx<'_>, nodal_load_adjoint: &[f64],
    ) -> Result<Vec<f64>, ConductionError> {
        vector(cx, nodal_load_adjoint, self.vertices)?;
        let mut result = Vec::with_capacity(self.factors.len());
        for factor in &self.factors {
            let mut value = 0.0;
            for (vertex, (&injection, &weight)) in factor.left.iter().zip(nodal_load_adjoint).enumerate() {
                if vertex % 512 == 0 { poll(cx, vertex)?; }
                value = checked(injection.mul_add(weight, value))?;
            }
            result.push(value);
        }
        poll(cx, result.len())?;
        Ok(result)
    }

    pub(in crate::adjoint::robin::nonlinear) fn check_vertices(
        &self, expected: usize,
    ) -> Result<(), ConductionError> {
        if self.vertices == expected { Ok(()) } else {
            Err(ConductionError::FieldLength {
                field: "mean-flux feedback binding", expected, found: self.vertices,
            })
        }
    }

    // The existing FGMRES operator already owns the free map. Direct indexing
    // avoids duplicating all retained full-node factors in reduced coordinates.
    pub(in crate::adjoint::robin::nonlinear) fn add_free_action(
        &self, free: &[usize], x: &[f64], y: &mut [f64], transposed: bool,
    ) {
        for factor in &self.factors {
            let (left, right) = if transposed { (&factor.right, &factor.left) }
                else { (&factor.left, &factor.right) };
            let projection: f64 = free.iter().zip(x).map(|(&vertex, &value)| right[vertex] * value).sum();
            for (&vertex, value) in free.iter().zip(y.iter_mut()) {
                *value += left[vertex] * projection;
            }
        }
    }

    pub(in crate::adjoint::robin::nonlinear) fn add_prescribed_pullback(
        &self, cx: &Cx<'_>, nodal_load_adjoint: &[f64], fixed: &[bool], result: &mut [f64],
    ) -> Result<(), ConductionError> {
        if fixed.len() != self.vertices || result.len() != self.vertices {
            return Err(invalid("mean-flux prescribed lift lost its full-node layout"));
        }
        let bars = self.flux_pullback(cx, nodal_load_adjoint)?;
        for (factor, bar) in self.factors.iter().zip(bars) {
            for (vertex, value) in result.iter_mut().enumerate() {
                if vertex % 512 == 0 { poll(cx, vertex)?; }
                if fixed[vertex] {
                    *value = checked(*value - checked(bar * factor.right[vertex])?)?;
                }
            }
        }
        poll(cx, self.vertices)
    }
}
