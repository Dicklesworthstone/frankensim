//! Original connected-component pressure bases, with every compression mode retained.
use super::{CavityModes, VibroError, invalid};

/// Original sealed compartments in one pressure coordinate system.
///
/// All uniform zero modes come first in compartment order. Higher modes then
/// follow frequency order with stable compartment/local-index tie breaking.
/// Interface samples concatenate in compartment order; each original pressure
/// trace has zero support on every other compartment's moving-wall samples.
/// Structural shapes/quadrature supplied to assemble_coupling must use this
/// SAME concatenated interface and the original outward cavity normals.
/// Gas density, equilibrium sound speed and scalar loss must agree exactly.
#[derive(Debug, Clone)]
pub struct CavityCompartmentBasis {
    pub(super) modes: CavityModes,
    pub(super) mode_indices: Vec<Vec<usize>>,
    interface_ranges: Vec<core::ops::Range<usize>>,
}
impl CavityCompartmentBasis {
    /// maximum_trace_entries bounds the complete zero-padded pressure trace
    /// payload in f64 entries, including zeros on other compartments.
    pub fn try_new(
        compartments: &[CavityModes],
        maximum_trace_entries: usize,
    ) -> Result<Self, VibroError> {
        Self::try_new_controlled(compartments, maximum_trace_entries, &mut || false)
    }

    /// The same bounded construction, polling between original compartments
    /// and output modes. A Cx caller may pass || cx.checkpoint().is_err().
    pub fn try_new_controlled(
        compartments: &[CavityModes],
        maximum_trace_entries: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, VibroError> {
        let poll = |cancelled: &mut dyn FnMut() -> bool| -> Result<(), VibroError> {
            if cancelled() {
                Err(VibroError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll(cancelled)?;
        if maximum_trace_entries == 0 {
            return Err(invalid(
                "combined cavity trace payload needs a positive explicit entry bound",
            ));
        }

        if compartments.len() < 2 || compartments.len() > 64 {
            return Err(invalid(
                "connected acoustic basis requires 2..=64 original compartments",
            ));
        }
        let reference = &compartments[0];
        if !reference.rho0.is_finite()
            || reference.rho0 <= 0.0
            || !reference.c0.is_finite()
            || reference.c0 <= 0.0
            || !reference.loss_factor.is_finite()
            || reference.loss_factor < 0.0
        {
            return Err(invalid(
                "compartments require one positive finite reference medium and loss",
            ));
        }
        let mut count = 0usize;
        let mut points = 0usize;
        let mut ranges = Vec::with_capacity(compartments.len());
        let mut order = Vec::new();
        let mut higher = Vec::new();
        for (c, basis) in compartments.iter().enumerate() {
            poll(cancelled)?;
            let n = basis.omegas.len();
            if n == 0
                || basis.lambdas.len() != n
                || basis.interface.len() != n
                || basis.interface[0].is_empty()
                || basis.omegas[0] != 0.0
                || basis.omegas[1..]
                    .iter()
                    .any(|w| !w.is_finite() || *w <= 0.0)
                || basis.lambdas.iter().any(|v| !v.is_finite() || *v <= 0.0)
                || basis.rho0 != reference.rho0
                || basis.c0 != reference.c0
                || basis.loss_factor != reference.loss_factor
            {
                return Err(invalid(
                    "each sealed compartment needs its original uniform mode, norms and common medium",
                ));
            }
            let width = basis.interface[0].len();
            let uniform = basis.interface[0][0];
            if !uniform.is_finite()
                || uniform == 0.0
                || basis.interface[0].iter().any(|p| *p != uniform)
                || basis
                    .interface
                    .iter()
                    .any(|row| row.len() != width || row.iter().any(|p| !p.is_finite()))
            {
                return Err(invalid(
                    "compartment pressure traces must retain a nonzero uniform mode and finite local interface",
                ));
            }
            count = count
                .checked_add(n)
                .ok_or_else(|| invalid("combined cavity mode count overflow"))?;
            let end = points
                .checked_add(width)
                .ok_or_else(|| invalid("combined cavity interface count overflow"))?;
            let payload = count
                .checked_mul(end)
                .ok_or_else(|| invalid("combined cavity trace storage overflow"))?;
            if payload > maximum_trace_entries {
                return Err(invalid(
                    "combined cavity pressure traces exceed the declared payload bound",
                ));
            }
            higher
                .try_reserve_exact(n - 1)
                .map_err(|_| invalid("combined cavity ordering allocation failed"))?;
            ranges.push(points..end);
            points = end;
            order.push((c, 0));
            higher.extend((1..n).map(|q| (c, q)));
        }
        count
            .checked_mul(points)
            .ok_or_else(|| invalid("combined cavity trace storage overflow"))?;
        poll(cancelled)?;
        higher.sort_unstable_by(|&(c, q), &(d, r)| {
            compartments[c].omegas[q]
                .total_cmp(&compartments[d].omegas[r])
                .then(c.cmp(&d))
                .then(q.cmp(&r))
        });
        order
            .try_reserve_exact(higher.len())
            .map_err(|_| invalid("combined cavity order allocation failed"))?;
        order.extend(higher);
        let mut indices = Vec::with_capacity(compartments.len());
        for basis in compartments {
            let mut map = Vec::new();
            map.try_reserve_exact(basis.omegas.len())
                .map_err(|_| invalid("combined cavity index allocation failed"))?;
            map.resize(basis.omegas.len(), 0);
            indices.push(map);
        }
        let mut modes = CavityModes {
            omegas: Vec::new(),
            lambdas: Vec::new(),
            interface: Vec::new(),
            loss_factor: reference.loss_factor,
            rho0: reference.rho0,
            c0: reference.c0,
        };
        modes
            .omegas
            .try_reserve_exact(count)
            .map_err(|_| invalid("combined cavity frequencies allocation failed"))?;
        modes
            .lambdas
            .try_reserve_exact(count)
            .map_err(|_| invalid("combined cavity norms allocation failed"))?;
        modes
            .interface
            .try_reserve_exact(count)
            .map_err(|_| invalid("combined cavity traces allocation failed"))?;
        for (global, (c, q)) in order.into_iter().enumerate() {
            poll(cancelled)?;
            indices[c][q] = global;
            modes.omegas.push(compartments[c].omegas[q]);
            modes.lambdas.push(compartments[c].lambdas[q]);
            let mut trace = Vec::new();
            trace
                .try_reserve_exact(points)
                .map_err(|_| invalid("combined cavity trace payload allocation failed"))?;
            trace.resize(points, 0.0);
            trace[ranges[c].clone()].copy_from_slice(&compartments[c].interface[q]);
            modes.interface.push(trace);
        }
        poll(cancelled)?;
        Ok(Self {
            modes,
            mode_indices: indices,
            interface_ranges: ranges,
        })
    }
    pub fn modes(&self) -> &CavityModes {
        &self.modes
    }
    pub fn compartment_count(&self) -> usize {
        self.mode_indices.len()
    }
    /// Local pressure-mode order mapped to the combined original coordinates.
    pub fn mode_indices(&self, compartment: usize) -> Option<&[usize]> {
        self.mode_indices.get(compartment).map(Vec::as_slice)
    }
    pub fn interface_ranges(&self) -> &[core::ops::Range<usize>] {
        &self.interface_ranges
    }
}
