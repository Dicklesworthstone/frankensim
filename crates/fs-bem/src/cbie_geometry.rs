//! Frequency-independent triangle integrals for repeated plain-CBIE solves.
use super::*;

/// Immutable geometry preparation for a sweep of plain-CBIE radiation solves.
///
/// The exact static single/double layers depend only on the source triangles
/// and target centroids. Compute them once, retaining two `f64` values per
/// panel pair; every frequency still integrates the original dynamic
/// remainders and builds and factors its own Helmholtz operator. No frequency
/// interpolation, source-mode reduction, or stored quadrature-distance bank
/// is involved. The borrowed surface cannot be changed or cross-wired.
///
/// This type does not establish that plain CBIE avoids interior resonances.
/// The caller must select that formulation for its physical boundary/band.
/// The existing uncached single/batch functions remain appropriate for a
/// one-shot solve or a caller that does not want the additional storage.
pub struct PreparedCbieGeometry<'a> {
    surface: &'a SpherePanels,
    static_layers: Vec<(f64, f64)>,
}

impl<'a> PreparedCbieGeometry<'a> {
    /// Prepare geometry after preflighting the highest planned wavenumber.
    /// `max_cache_bytes` bounds the logical static-table storage, exactly
    /// `16 * panels * panels` bytes; it does not include the dense per-frequency
    /// matrices/LU or allocator overhead. No table is allocated or integrated
    /// until the resolution and byte-budget checks pass.
    ///
    /// # Errors
    /// Missing triangles, invalid wavenumber, the original dense/resolution
    /// cap, insufficient cache budget, or failure to reserve bounded storage.
    pub fn new(
        surface: &'a SpherePanels,
        highest_wavenumber: f64,
        max_cache_bytes: usize,
    ) -> Result<Self, HelmholtzError> {
        Self::new_with_cancel(surface, highest_wavenumber, max_cache_bytes, || false)
    }

    /// The same preparation with cancellation before allocation and between
    /// target rows. A refusal publishes no partially prepared object. The
    /// callback does not interrupt a subsequent frequency assembly or LU.
    ///
    /// # Errors
    /// The same admissions as [`Self::new`], or [`HelmholtzError::Cancelled`].
    pub fn new_with_cancel(
        surface: &'a SpherePanels,
        highest_wavenumber: f64,
        max_cache_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, HelmholtzError> {
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        admit_resolution(surface, highest_wavenumber)?;
        let triangles = surface.triangles().ok_or(HelmholtzError::BadParameter {
            what: "CBIE geometry preparation requires retained source triangles",
        })?;
        let n = surface.centroids().len();
        // The dense panel cap was checked above, so these products fit usize.
        let entries = n * n;
        let required_bytes = entries * core::mem::size_of::<(f64, f64)>();
        if required_bytes > max_cache_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes,
                maximum_bytes: max_cache_bytes,
            });
        }
        let mut static_layers = Vec::new();
        static_layers
            .try_reserve_exact(entries)
            .map_err(|_| HelmholtzError::AllocationFailed {
                what: "CBIE static geometry table",
            })?;
        for (i, &target) in surface.centroids().iter().enumerate() {
            if cancelled() {
                return Err(HelmholtzError::Cancelled);
            }
            for (j, &triangle) in triangles.iter().enumerate() {
                // Preserve the self-term's original normal calculation too,
                // so the cached and uncached operation graphs agree bitwise.
                let layers = if i == j {
                    (triangle_self_static(target, triangle), 0.0)
                } else {
                    triangle_static_influence(target, triangle, surface.normals()[j])
                };
                static_layers.push(layers);
            }
        }
        if cancelled() {
            return Err(HelmholtzError::Cancelled);
        }
        Ok(Self {
            surface,
            static_layers,
        })
    }

    /// Logical bytes retained by the two dense static-layer coefficient tables.
    #[must_use]
    pub fn cache_bytes(&self) -> usize {
        self.static_layers.len() * core::mem::size_of::<(f64, f64)>()
    }

    /// Solve one frequency with the prepared geometry and a fresh shared LU.
    /// Returned solutions preserve the exact input order and semantics of
    /// [`solve_radiation_batch`] with [`Formulation::PlainCbie`]. The geometry
    /// is independent of frequency and medium; their original admissions run
    /// on every call, including any frequency above the preparation preflight.
    ///
    /// # Errors
    /// The same parameter, field-count/shape, resolution and factorization
    /// refusals as the ordinary plain-CBIE batch entry point.
    pub fn solve_batch(
        &self,
        k: f64,
        medium: Medium,
        velocity_fields: &[&[C64]],
    ) -> Result<Vec<RadiationSolution>, HelmholtzError> {
        let operator = RadiationOperator::prepare_with_static(
            self.surface,
            k,
            medium,
            velocity_fields,
            Formulation::PlainCbie,
            Some(&self.static_layers),
        )?;
        Ok(velocity_fields
            .iter()
            .map(|&velocity| operator.solve(velocity))
            .collect())
    }
}
