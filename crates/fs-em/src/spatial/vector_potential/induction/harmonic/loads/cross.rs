//! Exact two-carrier mixing of the retained linear induction fields.

use super::{
    CellLoad, SpatialConductorMaterial, SpatialHarmonicInductionResponse, SpatialLorentzLoadError,
    SpatialLorentzLoadLimits, SpatialLorentzLoads, integrate,
};
use fs_math::c64::C64;

/// Cross terms in `(J_a + J_b) cross (B_a + B_b)` on one original body.
///
/// The self terms are supplied separately by each response's `lorentz_loads`.
/// Add this interaction once per unordered carrier pair. The difference and
/// sum coefficients are PEAK amplitudes under `exp(-i omega t)`, so each real
/// component multiplies cosine and each imaginary component positive sine.
/// At equal carrier frequencies the difference-real coefficient is the full
/// additional DC force; its imaginary part is exactly zero. No independence
/// or cycle averaging is imposed on the two sources.
#[derive(Debug)]
pub struct SpatialHarmonicLorentzCrossLoads {
    difference_omega_rad_s: f64,
    sum_omega_rad_s: f64,
    difference_real: SpatialLorentzLoads,
    difference_imaginary: SpatialLorentzLoads,
    sum_real: SpatialLorentzLoads,
    sum_imaginary: SpatialLorentzLoads,
}

impl SpatialHarmonicLorentzCrossLoads {
    /// Nonnegative angular beat frequency `abs(omega_a - omega_b)` [rad/s].
    pub fn difference_omega_rad_s(&self) -> f64 {
        self.difference_omega_rad_s
    }

    /// Angular sum frequency `omega_a + omega_b` [rad/s].
    pub fn sum_omega_rad_s(&self) -> f64 {
        self.sum_omega_rad_s
    }

    /// Cosine coefficient at the difference frequency; DC when it is zero.
    pub fn difference_real(&self) -> &SpatialLorentzLoads {
        &self.difference_real
    }

    /// Sine coefficient at the difference frequency; zero for equal carriers.
    pub fn difference_imaginary(&self) -> &SpatialLorentzLoads {
        &self.difference_imaginary
    }

    /// Cosine coefficient at the sum frequency.
    pub fn sum_real(&self) -> &SpatialLorentzLoads {
        &self.sum_real
    }

    /// Sine coefficient at the sum frequency.
    pub fn sum_imaginary(&self) -> &SpatialLorentzLoads {
        &self.sum_imaginary
    }
}

impl<C: SpatialConductorMaterial> SpatialHarmonicInductionResponse<'_, '_, C> {
    /// Resolve both sum- and difference-frequency conduction-current forces.
    ///
    /// Both responses must borrow the SAME original induction source and use
    /// the same prescribed steady nodal velocity. Their separate, accepted
    /// linear field solves retain each carrier's skin effect and phase lag.
    /// The source phasors must share a physical time origin. Distinct sources,
    /// material states or motion operators cannot be silently superposed.
    ///
    /// With `omega_a >= omega_b`, the retained peak force coefficients are
    /// `F_sum = (J_a cross B_b + J_b cross B_a)/2` and
    /// `F_difference = (J_a cross conj(B_b) + conj(J_b) cross B_a)/2`.
    /// These follow directly by expanding the two real fields BEFORE taking
    /// their cross product. Thus 1f/3f excitation also drives 2f/4f motion;
    /// close tones produce a slow beat instead of independent averaged loads.
    ///
    /// Every coefficient reuses consistent P1 mechanical integration, keeping
    /// nodal forces, torque and virtual work. Work is O(V + S log C + S log S),
    /// storage O(V + S), for V original vertices, S selected cells and C
    /// conducting cells. Four real nodal arrays are retained; `limits` bounds
    /// each array and the selected-cell workspace. No new field solve, time
    /// sampling or quadratic storage in the mesh dimension is needed.
    ///
    /// This transfers prescribed linear AC conduction-current forces only.
    /// DC bias mixing, saturation, magnetostriction, impressed winding forces,
    /// interface traction and mechanical feedback remain outside this owner.
    /// Cancellation or refusal publishes no partial load and changes no state.
    pub fn cross_lorentz_loads(
        &self,
        other: &Self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicLorentzCrossLoads, SpatialLorentzLoadError> {
        use SpatialLorentzLoadError::{Cancelled, Invalid};
        if cancelled() {
            return Err(Cancelled);
        }
        let problem = self.problem();
        if !std::ptr::eq(self.source, other.source)
            || selected_cells.is_empty()
            || selected_cells.len() > limits.max_cells
            || problem.positions_m().len() > limits.max_vertices
        {
            return Err(Invalid(
                "cross-carrier loads require one original source and a nonempty bounded selection",
            ));
        }
        for (a, b) in self.velocity_m_s.iter().zip(&other.velocity_m_s) {
            if cancelled() {
                return Err(Cancelled);
            }
            if a != b {
                return Err(Invalid(
                    "cross-carrier loads require the same prescribed steady motion",
                ));
            }
        }
        // Canonical frequency order makes the beat frequency nonnegative and
        // ensures swapping the caller's two responses retains the same phase.
        let (a, b) = if self.omega_rad_s() >= other.omega_rad_s() {
            (self, other)
        } else {
            (other, self)
        };
        let difference_omega_rad_s = a.omega_rad_s() - b.omega_rad_s();
        let sum_omega_rad_s = a.omega_rad_s() + b.omega_rad_s();
        if !sum_omega_rad_s.is_finite() {
            return Err(Invalid("unrepresentable Lorentz sum frequency"));
        }
        let mut selected = Vec::with_capacity(selected_cells.len());
        for &cell in selected_cells {
            if cancelled() {
                return Err(Cancelled);
            }
            let index = a
                .cells()
                .binary_search_by_key(&cell, |record| record.cell)
                .map_err(|_| Invalid("cross-carrier loads require original conducting cells"))?;
            let ca = &a.cells()[index];
            let cb = &b.cells()[index];
            let ba = a.flux_density_t()[cell];
            let bb = b.flux_density_t()[cell];
            let mut density = [[[0.0; 3]; 4]; 4];
            for vertex in 0..4 {
                let ja = ca.current_density_a_m2[vertex];
                let jb = cb.current_density_a_m2[vertex];
                let sum_a = cross(ja, bb);
                let sum_b = cross(jb, ba);
                let difference_a = cross(ja, bb.map(C64::conj));
                let difference_b = cross(jb.map(C64::conj), ba);
                for axis in 0..3 {
                    // Half-sums avoid overflow from adding two large finite
                    // terms when their peak cross coefficient is representable.
                    density[0][vertex][axis] =
                        f64::midpoint(difference_a[axis].re, difference_b[axis].re);
                    density[1][vertex][axis] = if difference_omega_rad_s == 0.0 {
                        0.0
                    } else {
                        f64::midpoint(difference_a[axis].im, difference_b[axis].im)
                    };
                    density[2][vertex][axis] = f64::midpoint(sum_a[axis].re, sum_b[axis].re);
                    density[3][vertex][axis] = f64::midpoint(sum_a[axis].im, sum_b[axis].im);
                }
            }
            selected.push((cell, problem.mesh().tets[cell], ca.volume_m3, density));
        }
        selected.sort_unstable_by_key(|entry| entry.0);
        if selected.windows(2).any(|pair| pair[0].0 == pair[1].0) {
            return Err(Invalid(
                "cross-carrier load selection contains a repeated cell",
            ));
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
        let difference_real = coefficient(0)?;
        let difference_imaginary = coefficient(1)?;
        let sum_real = coefficient(2)?;
        let sum_imaginary = coefficient(3)?;
        if cancelled() {
            return Err(Cancelled);
        }
        Ok(SpatialHarmonicLorentzCrossLoads {
            difference_omega_rad_s,
            sum_omega_rad_s,
            difference_real,
            difference_imaginary,
            sum_real,
            sum_imaginary,
        })
    }
}

fn cross(a: [C64; 3], b: [C64; 3]) -> [C64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
