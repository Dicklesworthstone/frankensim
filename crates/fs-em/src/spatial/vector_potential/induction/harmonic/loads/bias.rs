//! Fundamental-frequency forces from an original, stationary magnetic bias.

use super::{
    CellLoad, SpatialConductorMaterial, SpatialHarmonicInductionResponse, SpatialLorentzLoadError,
    SpatialLorentzLoadLimits, SpatialLorentzLoads, integrate,
};
use crate::spatial::{
    cross,
    vector_potential::{SpatialVectorField, induction::SpatialInduction},
};
use fs_math::c64::C64;

/// Peak conduction-current force at the AC carrier in a static magnetic bias.
///
/// Add these coefficients to the AC response's separate mean and twice-carrier
/// loads. At phase `theta`, the extra force is `real*cos(theta) +
/// imaginary*sin(theta)`, using the native `exp(-i omega t)` convention.
/// There is no half-factor: only one factor in `J_ac cross B_dc` oscillates.
#[derive(Debug)]
pub struct SpatialHarmonicBiasLorentzLoads {
    omega_rad_s: f64,
    real: SpatialLorentzLoads,
    imaginary: SpatialLorentzLoads,
}

impl SpatialHarmonicBiasLorentzLoads {
    /// Original electromagnetic angular frequency [rad/s].
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Peak cosine coefficient in the original mesh coordinates.
    pub fn real(&self) -> &SpatialLorentzLoads {
        &self.real
    }

    /// Peak sine coefficient in the original mesh coordinates.
    pub fn imaginary(&self) -> &SpatialLorentzLoads {
        &self.imaginary
    }
}

// Share the original static-field admission with coherent finite-window
// forcing. The original response remains the geometry/material/motion owner.
pub(crate) fn admit_static_bias<C: SpatialConductorMaterial>(
    source: &SpatialInduction<'_, C>,
    velocity_m_s: &[[f64; 3]],
    bias: &SpatialVectorField<'_>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), SpatialLorentzLoadError> {
    use SpatialLorentzLoadError::{Cancelled, Invalid};
    if cancelled() {
        return Err(Cancelled);
    }
    if !std::ptr::eq(source.problem(), bias.problem()) {
        return Err(Invalid(
            "static bias requires one original source and a bounded nonempty selection",
        ));
    }
    for velocity in velocity_m_s {
        if cancelled() {
            return Err(Cancelled);
        }
        if velocity.iter().any(|v| *v != 0.0) {
            return Err(Invalid(
                "static-bias force requires stationary conductors without motional DC current",
            ));
        }
    }
    for (dc, ac) in bias
        .prescribed_edges()
        .iter()
        .zip(&source.prepared.fixed_edges)
    {
        if cancelled() {
            return Err(Cancelled);
        }
        if dc != ac {
            return Err(Invalid("static and harmonic magnetic trace masks differ"));
        }
    }
    Ok(())
}

pub(crate) fn admit_static_bias_cell(
    bias: &SpatialVectorField<'_>,
    cell: usize,
) -> Result<(), SpatialLorentzLoadError> {
    if bias.impressed_current_density_a_m2()[cell]
        .iter()
        .any(|j| *j != 0.0)
    {
        return Err(SpatialLorentzLoadError::Invalid(
            "static-bias selection contains DC current requiring the full DC/AC force",
        ));
    }
    Ok(())
}

impl<C: SpatialConductorMaterial> SpatialHarmonicInductionResponse<'_, '_, C> {
    /// Retain the fundamental-frequency force induced by a static magnetic bias.
    ///
    /// The accepted magnetostatic field must use the SAME original material,
    /// geometry and prescribed-trace mask as this response. Its trace values
    /// and gauge may differ: linear superposition retains the authored static
    /// field and the homogeneous AC trace independently. No field is reset or
    /// solved again, and the original AC skin effect and phase are retained.
    ///
    /// This owner admits stationary conductors and selected cells without
    /// impressed DC current. Thus their settled conduction current has no DC
    /// part and the added force is exactly `J_ac cross B_dc`. A moving or
    /// DC-driven conductor also needs `J_dc cross B_ac` from a coupled DC
    /// conduction solution; it is refused rather than silently omitting it.
    /// Static sources outside the selected cells may generate the bias.
    ///
    /// Selection names distinct original conducting tetrahedra. Two consistent
    /// P1 nodal-force arrays preserve wrench and virtual power with no vertex
    /// lumping. Work is O(E + V + S log C + S log S), storage O(V + S), for E
    /// original edges, V vertices, S selected cells and C conducting cells.
    /// `limits.max_vertices` bounds each retained nodal array. Refusal or
    /// cancellation returns no partial loads and advances no physical state.
    ///
    /// Linear magnetic sources, prescribed static bias and one-way mechanical
    /// forcing are the scope. This adds no magnetic saturation, recoil-envelope
    /// admission, DC terminal model, impressed-current force or motion feedback.
    pub fn static_bias_lorentz_loads(
        &self,
        bias: &SpatialVectorField<'_>,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicBiasLorentzLoads, SpatialLorentzLoadError> {
        use SpatialLorentzLoadError::{Cancelled, Invalid};
        if cancelled() {
            return Err(Cancelled);
        }
        let problem = self.problem();
        if !std::ptr::eq(problem, bias.problem())
            || selected_cells.is_empty()
            || selected_cells.len() > limits.max_cells
            || problem.positions_m().len() > limits.max_vertices
        {
            return Err(Invalid(
                "static bias requires one original source and a bounded nonempty selection",
            ));
        }
        admit_static_bias(self.source, &self.velocity_m_s, bias, cancelled)?;
        let mut selected = Vec::with_capacity(selected_cells.len());
        for &cell in selected_cells {
            if cancelled() {
                return Err(Cancelled);
            }
            let index = self
                .cells()
                .binary_search_by_key(&cell, |record| record.cell)
                .map_err(|_| Invalid("static-bias selection requires original conducting cells"))?;
            admit_static_bias_cell(bias, cell)?;
            let record = &self.cells()[index];
            selected.push((
                cell,
                problem.mesh().tets[cell],
                record.volume_m3,
                force_density(record.current_density_a_m2, bias.flux_density_t()[cell]),
            ));
        }
        selected.sort_unstable_by_key(|entry| entry.0);
        if selected.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(Invalid("static-bias selection contains a repeated cell"));
        }
        let mut coefficient = |index| {
            integrate(
                problem.positions_m(),
                selected
                    .iter()
                    .map(|(_, vertices, volume, density)| CellLoad {
                        vertices: *vertices,
                        volume_m3: *volume,
                        force_density_n_m3: &density[index],
                    }),
                limits,
                cancelled,
            )
        };
        let real = coefficient(0)?;
        let imaginary = coefficient(1)?;
        if cancelled() {
            return Err(Cancelled);
        }
        Ok(SpatialHarmonicBiasLorentzLoads {
            omega_rad_s: self.omega_rad_s(),
            real,
            imaginary,
        })
    }
}

fn force_density(current: [[C64; 3]; 4], bias: [f64; 3]) -> [[[f64; 3]; 4]; 2] {
    [
        current.map(|j| cross(j.map(|v| v.re), bias)),
        current.map(|j| cross(j.map(|v| v.im), bias)),
    ]
}

#[cfg(test)]
#[path = "bias/tests.rs"]
mod tests;
