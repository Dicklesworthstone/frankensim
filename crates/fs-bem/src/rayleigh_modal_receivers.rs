//! Finite-triangle observation of the original coherent Rayleigh cell flows.
//! All spatial integration and gradient gates belong to `near_field`.
use super::*;
use crate::near_field::{self, RayleighSingleLayer};

/// Physical clearance and cold integration/retained-storage admission.
#[derive(Clone, Copy, Debug)]
pub struct TriangleReceiverOptions {
    /// Actual distance from each receiver to every finite triangle [m].
    pub minimum_clearance_m: f64,
    /// Shared pressure/gradient quadrature controls over the complete array.
    pub integration: near_field::Options,
    /// Retained rows, point/clearance arrays and prepared owner. Borrowed
    /// source geometry, returned fields and bounded integration stack excluded.
    pub maximum_row_bytes: usize,
}

/// Co-located peak fields on the original exp(-i omega t) source clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriangleReceiverField {
    pub pressure_pa: C64,
    pub pressure_gradient_pa_m: [C64; 3],
    pub particle_velocity_m_s: [C64; 3],
    /// Propagated coarse/fine integration discrepancies, not continuum,
    /// structural mesh, material or relative dark-field error bounds.
    pub pressure_quadrature_error_estimate_pa: f64,
    pub gradient_quadrature_error_estimate_pa_m: [f64; 3],
    pub velocity_quadrature_error_estimate_m_s: [f64; 3],
}

/// Reusable finite-cell rows for one surface and wavenumber. Every evaluation
/// contracts original DD cell flows; it repeats no quadrature or boundary solve.
pub struct PreparedTriangleReceivers<'s> {
    surface: &'s PreparedSurface,
    k: f64,
    rows: Vec<RayleighSingleLayer>,
    points: Vec<[f64; 3]>,
    clearances: Vec<f64>,
    kernel_evaluations: usize,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), HelmholtzError> {
    if cancelled() { Err(HelmholtzError::Cancelled) } else { Ok(()) }
}

impl<'s> PreparedCoherentRadiation<'s> {
    /// Complete retained receiver payload, including the owner and point
    /// diagnostics, before allocation. The original 64-receiver cap applies.
    pub fn triangle_receiver_row_bytes(&self, receivers: usize) -> Result<usize, HelmholtzError> {
        if receivers > near_field::MAX_RECEIVERS {
            return Err(invalid("too many finite-cell Rayleigh receivers"));
        }
        self.surface.count.checked_mul(receivers)
            .and_then(|n| n.checked_mul(core::mem::size_of::<RayleighSingleLayer>()))
            .and_then(|n| n.checked_add(receivers * (core::mem::size_of::<[f64; 3]>() + core::mem::size_of::<f64>())))
            .and_then(|n| n.checked_add(core::mem::size_of::<PreparedTriangleReceivers<'_>>()))
            .ok_or_else(|| invalid("finite-cell Rayleigh receiver payload overflow"))
    }

    /// Integrate finite triangles instead of replacing nearby cells by points.
    /// Each triangle is in z=0 and matches the corresponding retained cell's
    /// centroid and area to floating-point roundoff. Its prescribed velocity
    /// is uniform, with the SAME complete coherent cell volume flow. Triangle
    /// winding is immaterial: the radiating half-space is explicitly z>=0.
    ///
    /// The existing adaptive P0 Green owner integrates pressure and all three
    /// analytic derivatives with their own discrepancy/separation/phase gates.
    /// The caller supplies the physical tessellation; overlapping patches add
    /// their prescribed velocities. No finite-cell shape is inferred from area,
    /// and no receiver is clamped
    /// away from a source. Touching cells, unresolved clearance, or exhausted
    /// work/refinement refuse. Empty receiver lists are admitted.
    ///
    /// Rows may observe another coherent motion of this same surface and
    /// wavenumber. This is finite-cell prescribed radiation, not a new loading
    /// solve or finite-baffle scattering. Its spatial model differs from the
    /// existing area-point far field by finite-cell quadrature error; refine
    /// the source cells to compare those two images.
    pub fn prepare_triangle_receivers(
        &self,
        triangles: &[[[f64; 3]; 3]],
        points: &[[f64; 3]],
        options: TriangleReceiverOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedTriangleReceivers<'s>, HelmholtzError> {
        poll(&mut cancelled)?;
        options.integration.validate()?;
        let bytes = self.triangle_receiver_row_bytes(points.len())?;
        if triangles.len() != self.surface.count
            || bytes > options.maximum_row_bytes
            || triangles.len().checked_mul(points.len()).and_then(|n| n.checked_mul(80))
                .is_none_or(|n| n > options.integration.maximum_kernel_evaluations)
            || !options.minimum_clearance_m.is_finite() || options.minimum_clearance_m <= 0.0 {
            return Err(invalid("invalid finite-cell Rayleigh geometry, clearance or budget"));
        }
        let mut coordinate_scale = 1.0_f64;
        for (cell, triangle) in triangles.iter().enumerate() {
            if cell % 32 == 0 { poll(&mut cancelled)?; }
            for vertex in triangle {
                if vertex.iter().any(|x| !x.is_finite() || x.abs() > 1e9) || vertex[2] != 0.0 {
                    return Err(invalid("Rayleigh receiver triangles must be finite in the original baffle plane"));
                }
                coordinate_scale = coordinate_scale.max(vertex[0].abs()).max(vertex[1].abs());
            }
            let a = [triangle[1][0] - triangle[0][0], triangle[1][1] - triangle[0][1]];
            let b = [triangle[2][0] - triangle[0][0], triangle[2][1] - triangle[0][1]];
            let area = 0.5 * (a[0] * b[1] - a[1] * b[0]).abs();
            if !area.is_finite() || area <= 0.0
                || (area - self.surface.areas[cell]).abs() > 128.0 * f64::EPSILON * area {
                return Err(invalid("Rayleigh receiver triangle area does not match its source cell"));
            }
            for axis in 0..2 {
                let centroid = triangle[0][axis] + a[axis] / 3.0 + b[axis] / 3.0;
                if (centroid - self.surface.positions[cell][axis]).abs() > 64.0 * f64::EPSILON * coordinate_scale {
                    return Err(invalid("Rayleigh receiver triangle centroid does not match its source cell"));
                }
            }
        }
        for point in points {
            poll(&mut cancelled)?;
            if point.iter().any(|x| !x.is_finite() || x.abs() > 1e9) || point[2] < 0.0 {
                return Err(invalid("Rayleigh receivers require finite positions in the radiating half-space"));
            }
            for x in point { coordinate_scale = coordinate_scale.max(x.abs()); }
        }
        if options.minimum_clearance_m <= 64.0 * f64::EPSILON * coordinate_scale {
            return Err(invalid("Rayleigh receiver clearance is below coordinate resolution"));
        }
        let mut clearances = reserved(points.len(), "finite-cell Rayleigh clearances")?;
        for &point in points {
            let mut clearance = f64::INFINITY;
            for (cell, &triangle) in triangles.iter().enumerate() {
                if cell % 32 == 0 { poll(&mut cancelled)?; }
                clearance = clearance.min(near_field::distance(point, triangle));
            }
            if !clearance.is_finite() || clearance <= options.minimum_clearance_m {
                return Err(invalid("Rayleigh receiver touches or is too close to a finite cell"));
            }
            clearances.push(clearance);
        }
        let mut retained_points = reserved(points.len(), "finite-cell Rayleigh positions")?;
        retained_points.extend_from_slice(points);
        let (rows, kernel_evaluations) = near_field::rayleigh_single_layer_rows(
            self.k, triangles, points, options.integration, &mut cancelled)?;
        Ok(PreparedTriangleReceivers {
            surface: self.surface, k: self.k, rows, points: retained_points, clearances, kernel_evaluations,
        })
    }
}

fn multiply(value: [Dd; 2], kernel: C64) -> [Dd; 2] {
    let (re, im) = (Dd::from_f64(kernel.re), Dd::from_f64(kernel.im));
    [value[0] * re - value[1] * im, value[0] * im + value[1] * re]
}

impl PreparedTriangleReceivers<'_> {
    pub fn points_m(&self) -> &[[f64; 3]] { &self.points }
    pub fn clearances_m(&self) -> &[f64] { &self.clearances }
    /// Actual original-kernel evaluations, including rejected coarse tiles.
    pub fn kernel_evaluations(&self) -> usize { self.kernel_evaluations }

    /// Observe any prepared coherent motion of the same surface and wavenumber.
    /// Different source density may reuse the Green rows; its original pressure
    /// scale is applied here. Original cell flows combine before field rounding.
    /// Velocity is -grad(potential), retaining its regular low-frequency limit;
    /// no division by a small omega or plane-wave pressure/velocity substitution.
    /// Return no partial observation on cancellation, mismatch or overflow.
    pub fn evaluate(
        &self,
        source: &PreparedCoherentRadiation<'_>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<TriangleReceiverField>, HelmholtzError> {
        poll(&mut cancelled)?;
        if !core::ptr::eq(self.surface, source.surface) || self.k.to_bits() != source.k.to_bits() {
            return Err(invalid("Rayleigh receiver rows and coherent source differ"));
        }
        let mut fields = reserved(self.points.len(), "finite-cell Rayleigh observations")?;
        let pressure = |value: [Dd; 2]| {
            let scale = Dd::from_f64(source.pressure_scale);
            let tau = Dd::from_f64(2.0 * PI);
            C64::new((value[1] * scale * tau).to_f64(), (-value[0] * scale * tau).to_f64())
        };
        for row in self.rows.chunks_exact(self.surface.count) {
            poll(&mut cancelled)?;
            let mut potential = [[Dd::ZERO; 2]; 4];
            let mut error = [Dd::ZERO; 4];
            for (cell, (kernel, &flow)) in row.iter().zip(&source.flow).enumerate() {
                if cell % 32 == 0 { poll(&mut cancelled)?; }
                let magnitude = Dd::from_f64(det::hypot(flow[0].to_f64(), flow[1].to_f64()));
                for (component, (value, discrepancy)) in [kernel.value].into_iter().chain(kernel.gradient)
                    .zip([kernel.error_value].into_iter().chain(kernel.error_gradient)).enumerate() {
                    let term = multiply(flow, value);
                    for axis in 0..2 { potential[component][axis] = potential[component][axis] + term[axis]; }
                    error[component] = error[component] + magnitude * Dd::from_f64(discrepancy);
                }
            }
            let pressure_error = |value: Dd| (value * Dd::from_f64(source.pressure_scale) * Dd::from_f64(2.0 * PI)).to_f64();
            let field = TriangleReceiverField {
                pressure_pa: pressure(potential[0]),
                pressure_gradient_pa_m: std::array::from_fn(|axis| pressure(potential[axis + 1])),
                particle_velocity_m_s: std::array::from_fn(|axis| C64::new(
                    -potential[axis + 1][0].to_f64(), -potential[axis + 1][1].to_f64())),
                pressure_quadrature_error_estimate_pa: pressure_error(error[0]),
                gradient_quadrature_error_estimate_pa_m: std::array::from_fn(|axis| pressure_error(error[axis + 1])),
                velocity_quadrature_error_estimate_m_s: std::array::from_fn(|axis| error[axis + 1].to_f64()),
            };
            if [field.pressure_pa].into_iter().chain(field.pressure_gradient_pa_m).chain(field.particle_velocity_m_s)
                .any(|value| !value.re.is_finite() || !value.im.is_finite())
                || [field.pressure_quadrature_error_estimate_pa].into_iter().chain(field.gradient_quadrature_error_estimate_pa_m)
                    .chain(field.velocity_quadrature_error_estimate_m_s).any(|e| !e.is_finite()) {
                return Err(invalid("finite-cell Rayleigh observation exceeds finite numerical range"));
            }
            fields.push(field);
        }
        poll(&mut cancelled)?;
        Ok(fields)
    }
}

#[cfg(test)]
#[path = "rayleigh_modal_receiver_tests.rs"]
mod tests;
