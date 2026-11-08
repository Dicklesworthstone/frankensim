//! Full-turn azimuthal Fourier components of the original meridian FE problem.
//!
//! For psi(r,phi,z)=psi_m(r,z) exp(i*m*phi) and exp(+i*omega*t), Euler gives
//! v=-(psi_r, i*m*psi/r, psi_z)/rho. The added stiffness is
//! K_phi=2*pi*m² integral(phi_i phi_j/(rho*r) dr dz). See the independently
//! separated Helmholtz equation in oomph-lib's original formulation:
//! <https://oomph-lib.github.io/oomph-lib/doc/fourier_decomposed_helmholtz/sphere_scattering/html/index.html>.
//!
//! Each signed order is one COMPLETE complex component, with the same 2*pi*r
//! volume and boundary-work normalization as m=0. Real cos/sin angular shapes
//! require their explicit coherent +m/-m amplitudes; no extra degeneracy or
//! factor of two is inferred. Geometry, materials and boundary coefficients
//! must be independent of phi. The supplied normal-velocity trace is its m-th
//! Fourier amplitude. There is no coupling between different orders here.
//!
//! Nonzero orders enforce zero axis pressure by coordinate restriction. This
//! is the finite-energy weak condition. P1 fields do not promise the exact
//! pointwise O(r^|m|) regularity or a unique Cartesian gradient at axis vertices;
//! such vector probes refuse, and physical h refinement remains necessary.

use super::*;
use crate::highorder::quad1d::gauss_legendre;
use fs_math::det;

/// Explicit integer-order admission; angular evaluation has fixed work.
pub const MAX_AZIMUTHAL_ORDER: u32 = 1024;
/// Maximum requested coarse quadrature; the retained rule uses twice this.
pub const MAX_FOURIER_QUADRATURE_POINTS: usize = 32;

/// Original triangle integration, independent of the meridian h refinement.
#[derive(Clone, Copy, Debug)]
pub struct FourierQuadrature {
    /// Native Gauss--Legendre points in the remaining collapsed coordinate,
    /// in 2..=32. Assembly retains the 2*q result and compares it with q.
    pub points: usize,
    /// Largest relative q-to-2q change of any nonzero angular element entry,
    /// in (0,1). This convergence diagnostic is not an interval error bound,
    /// and excludes floating-point, spatial and geometry error.
    pub maximum_relative_difference: f64,
}

pub(super) struct FourierForm {
    order: i32,
    pub(super) stiffness: Csr,
    axis: Vec<usize>,
    maximum_relative_difference: f64,
}

impl<'m> AxisymmetricPressureFe<'m> {
    /// Resolve one signed circumferential order in the ORIGINAL pressure space.
    ///
    /// All original cells and non-axis pressure coordinates remain. m=0 calls
    /// the original constructor without angular work; its quadrature argument
    /// is unused. Other orders add their angular stiffness to the existing
    /// sparse matrix, and preparation removes only the physical axis and
    /// authored pressure-release coordinates before the native harmonic solve.
    ///
    /// Angular assembly reuses the original triplet/cell/DOF limits. It adds
    /// one sparse matrix with the original sparsity and O(DOF) axis storage.
    /// Work is O(cells*q), plus the original cancellable sparse assembly. Native
    /// quadrature generation is a bounded tile of at most 64 points; all cell,
    /// quadrature-point, assembly and merge loops poll the original Cx.
    pub fn new_fourier(
        mesh: &'m TriComplex2,
        materials: &[AxisymmetricFluid],
        order: i32,
        quadrature: FourierQuadrature,
        limits: WeightedAssemblyLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, AxisymmetricError> {
        poll(cx)?;
        if order == 0 {
            return Self::new(mesh, materials, limits, cx);
        }
        if order.unsigned_abs() > MAX_AZIMUTHAL_ORDER
            || !(2..=MAX_FOURIER_QUADRATURE_POINTS).contains(&quadrature.points)
            || !quadrature.maximum_relative_difference.is_finite()
            || quadrature.maximum_relative_difference <= 0.0
            || quadrature.maximum_relative_difference >= 1.0
        {
            return Err(AxisymmetricError::Invalid(
                "Fourier order or radial quadrature admission",
            ));
        }
        let mut model = Self::new(mesh, materials, limits, cx)?;
        let coarse = gauss_legendre(quadrature.points);
        poll(cx)?;
        let fine = gauss_legendre(2 * quadrature.points);
        poll(cx)?;
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(9 * mesh.faces().len())
            .map_err(|_| WeightedError::Allocation)?;
        let mut maximum_relative_difference = 0.0_f64;
        for (cell, &face) in mesh.faces().iter().enumerate() {
            poll(cx)?;
            let points = face.map(|v| mesh.vertices()[v as usize]);
            let low = angular_element(points, materials[cell].density_kg_m3, order, &coarse, cx)?;
            let high = angular_element(points, materials[cell].density_kg_m3, order, &fine, cx)?;
            for i in 0..3 {
                for j in 0..3 {
                    let difference = if high[i][j] == 0.0 {
                        0.0
                    } else {
                        (high[i][j] - low[i][j]).abs() / high[i][j]
                    };
                    maximum_relative_difference = maximum_relative_difference.max(difference);
                    if !difference.is_finite()
                        || difference > quadrature.maximum_relative_difference
                    {
                        return Err(AxisymmetricError::Invalid(
                            "unresolved angular stiffness quadrature",
                        ));
                    }
                    entries.push((face[i] as usize, face[j] as usize, high[i][j]));
                }
            }
        }
        let stiffness = assemble_triplets(model.ndof(), entries, &mut || cx.checkpoint().is_err())?;
        // Both native assemblies retained every original cell pair, including
        // zero entries. Their identical canonical sparsity permits this O(nnz)
        // coefficient addition without another assembly or dense allocation.
        for row in 0..model.ndof() {
            poll(cx)?;
            let (columns, values) = stiffness.row(row);
            if columns != model.stiffness.row(row).0 {
                return Err(AxisymmetricError::Invalid(
                    "original Fourier stiffness sparsity changed",
                ));
            }
            for (target, value) in model.stiffness.row_values_mut(row).iter_mut().zip(values) {
                poll(cx)?;
                *target += value;
                if !target.is_finite() {
                    return Err(AxisymmetricError::Invalid(
                        "unrepresentable total Fourier stiffness",
                    ));
                }
            }
        }
        let mut axis = Vec::new();
        for (vertex, point) in mesh.vertices().iter().enumerate() {
            poll(cx)?;
            if point[0] == 0.0 {
                axis.push(vertex);
            }
        }
        model.fourier = Some(FourierForm {
            order,
            stiffness,
            axis,
            maximum_relative_difference,
        });
        poll(cx)?;
        Ok(model)
    }

    #[must_use]
    pub fn azimuthal_order(&self) -> i32 {
        self.fourier.as_ref().map_or(0, |f| f.order)
    }

    /// Exact zero-pressure constraints required by finite angular energy.
    /// The public full mass/stiffness forms act on vectors obeying these
    /// constraints. The original m=0 leaves all axis coordinates live.
    #[must_use]
    pub fn axis_pressure_constraints(&self) -> &[usize] {
        self.fourier.as_ref().map_or(&[], |f| f.axis.as_slice())
    }

    /// Angular contribution already included in [`Self::stiffness`].
    #[must_use]
    pub fn azimuthal_stiffness(&self) -> Option<&Csr> {
        self.fourier.as_ref().map(|f| &f.stiffness)
    }

    #[must_use]
    pub fn maximum_angular_quadrature_difference(&self) -> f64 {
        self.fourier
            .as_ref()
            .map_or(0.0, |f| f.maximum_relative_difference)
    }
}

/// Exact collapse in s followed by native positive quadrature in t:
/// lambda=(1-s,s*(1-t),s*t), r=a*(1-s)+s*b(t), dA=|det|*s ds dt.
/// The minimum-radius vertex is the collapse tip. At a=0, its constrained
/// shape never enters; the retained block has integral s² ds / b=1/(3b).
/// No singular axis basis, small-radius substitute or logarithmic cutoff occurs.
fn angular_element(
    points: [[f64; 2]; 3],
    density: f64,
    order: i32,
    rule: &(Vec<f64>, Vec<f64>),
    cx: &Cx<'_>,
) -> Result<[[f64; 3]; 3], AxisymmetricError> {
    let mut permutation = [0, 1, 2];
    permutation.sort_by(|&i, &j| {
        points[i][0]
            .total_cmp(&points[j][0])
            .then(points[i][1].total_cmp(&points[j][1]))
    });
    let radius = permutation.map(|i| points[i][0]);
    let scale_radius = radius[2];
    let scaled = radius.map(|r| r / scale_radius);
    let determinant = ((points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
        - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1]))
        .abs();
    let multiplier =
        (determinant / scale_radius) / density * core::f64::consts::TAU * f64::from(order).powi(2);
    if !multiplier.is_finite()
        || multiplier <= 0.0
        || radius
            .iter()
            .zip(scaled)
            .any(|(&r, q)| !q.is_finite() || (r > 0.0 && q == 0.0))
    {
        return Err(AxisymmetricError::Invalid(
            "unrepresentable angular radial geometry",
        ));
    }
    let mut matrix = [[Dd::ZERO; 3]; 3];
    for (&node, &weight) in rule.0.iter().zip(&rule.1) {
        poll(cx)?;
        let t = 0.5 * (1.0 + node);
        let b = (1.0 - t) * scaled[1] + t * scaled[2];
        let [a00, a01, a11] = collapsed_moments(scaled[0], b)?;
        let values = [
            [a00, a01 * (1.0 - t), a01 * t],
            [
                a01 * (1.0 - t),
                a11 * (1.0 - t).powi(2),
                a11 * (1.0 - t) * t,
            ],
            [a01 * t, a11 * (1.0 - t) * t, a11 * t * t],
        ];
        for i in 0..3 {
            for j in i..3 {
                if radius[i] == 0.0 || radius[j] == 0.0 {
                    continue;
                }
                let value = values[i][j] * (0.5 * weight);
                if !value.is_finite() || value <= 0.0 {
                    return Err(AxisymmetricError::Invalid(
                        "unrepresentable angular quadrature term",
                    ));
                }
                matrix[i][j] = matrix[i][j] + Dd::from_f64(value);
            }
        }
    }
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in i..3 {
            let value = (matrix[i][j] * Dd::from_f64(multiplier)).to_f64();
            if !value.is_finite() || (radius[i] > 0.0 && radius[j] > 0.0 && value <= 0.0) {
                return Err(AxisymmetricError::Invalid(
                    "unrepresentable angular element",
                ));
            }
            result[permutation[i]][permutation[j]] = value;
            result[permutation[j]][permutation[i]] = value;
        }
    }
    Ok(result)
}

// Integrals of [s(1-s)^2, s²(1-s), s³]/[a+(b-a)s] over [0,1].
// A short convergent series avoids cancellation near equal radii; otherwise
// polynomial division gives J_n=(1/n-a*J_(n-1))/(b-a). No fitted coefficients.
fn collapsed_moments(a: f64, b: f64) -> Result<[f64; 3], AxisymmetricError> {
    if a == 0.0 {
        return Ok([0.0, 0.0, (1.0 / b) / 3.0]);
    }
    let delta = (b - a) / a;
    let result = if delta.abs() <= 0.25 {
        let mut sum = [Dd::ZERO; 3];
        let mut power = 1.0;
        for n in 0..28 {
            let n = f64::from(n);
            let moments = [
                2.0 / ((n + 2.0) * (n + 3.0) * (n + 4.0)),
                1.0 / ((n + 3.0) * (n + 4.0)),
                1.0 / (n + 4.0),
            ];
            for i in 0..3 {
                sum[i] = sum[i] + Dd::from_f64(power) * Dd::from_f64(moments[i]);
            }
            power *= -delta;
        }
        sum.map(|v| v.to_f64() / a)
    } else {
        let width = Dd::from_f64(b - a);
        let a_dd = Dd::from_f64(a);
        let ratio = b / a;
        let log_ratio = if ratio.is_finite() {
            det::ln(ratio)
        } else {
            det::ln(b) - det::ln(a)
        };
        let j0 = Dd::from_f64(log_ratio) / width;
        let j1 = (Dd::from_f64(1.0) - a_dd * j0) / width;
        let j2 = (Dd::from_f64(0.5) - a_dd * j1) / width;
        let j3 = (Dd::from_f64(1.0 / 3.0) - a_dd * j2) / width;
        [
            (j1 - j2 - j2 + j3).to_f64(),
            (j2 - j3).to_f64(),
            j3.to_f64(),
        ]
    };
    if result.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(AxisymmetricError::Invalid(
            "unrepresentable collapsed angular moment",
        ));
    }
    Ok(result)
}

/// Physical three-dimensional point for one complete signed Fourier component.
#[derive(Clone, Copy, Debug)]
pub struct FourierAcousticPoint {
    pub position_xyz_m: [f64; 3],
    pub pressure_pa: C64,
    pub pressure_potential_pa_s: C64,
    pub particle_velocity_r_phi_z_m_s: [C64; 3],
    pub particle_velocity_xyz_m_s: [C64; 3],
    /// Peak-phasor active intensity Re(p conj(v))/2 [W/m²].
    pub active_intensity_xyz_w_m2: [f64; 3],
}

impl AxisymmetricHarmonicField<'_, '_> {
    /// Reconstruct the actual three-dimensional phase, Euler velocity and
    /// intensity. Angle is measured from +x towards +y about the +z axis.
    /// Nonzero-order vector samples exactly on the axis refuse because a P1
    /// cell does not supply a unique regular Cartesian limit. m=0 keeps the
    /// original finite cell-side radial-gradient interpretation there.
    pub fn probe_fourier(
        &self,
        cell: usize,
        barycentric: [f64; 3],
        azimuth_rad: f64,
    ) -> Result<FourierAcousticPoint, AxisymmetricError> {
        let base = self.probe(cell, barycentric)?;
        let order = self.model.azimuthal_order();
        let phase_angle = f64::from(order) * azimuth_rad;
        let r = base.position_r_z_m[0];
        if !azimuth_rad.is_finite() || !phase_angle.is_finite() || (order != 0 && r == 0.0) {
            return Err(AxisymmetricError::Invalid(
                "finite azimuth and off-axis Fourier vector probe required",
            ));
        }
        let (sin, cos) = det::sin_cos(azimuth_rad);
        let (phase_sin, phase_cos) = det::sin_cos(phase_angle);
        let phase = C64::new(phase_cos, phase_sin);
        let azimuthal = if order == 0 {
            C64::ZERO
        } else {
            let psi = base.pressure_potential_pa_s;
            C64::new(psi.re / r, psi.im / r)
                * C64::new(
                    0.0,
                    -f64::from(order) / self.model.materials[cell].density_kg_m3,
                )
        };
        let velocity = [
            base.particle_velocity_r_z_m_s[0],
            azimuthal,
            base.particle_velocity_r_z_m_s[1],
        ]
        .map(|v| v * phase);
        let cartesian = [
            velocity[0].scale(cos) - velocity[1].scale(sin),
            velocity[0].scale(sin) + velocity[1].scale(cos),
            velocity[2],
        ];
        let pressure_pa = base.pressure_pa * phase;
        let pressure_potential_pa_s = base.pressure_potential_pa_s * phase;
        let intensity = cartesian.map(|v| 0.5 * (pressure_pa * v.conj()).re);
        let position = [r * cos, r * sin, base.position_r_z_m[1]];
        if !finite(pressure_pa)
            || !finite(pressure_potential_pa_s)
            || velocity.iter().chain(&cartesian).any(|v| !finite(*v))
            || intensity.iter().chain(&position).any(|v| !v.is_finite())
        {
            return Err(AxisymmetricError::Invalid(
                "unrepresentable Fourier acoustic point",
            ));
        }
        Ok(FourierAcousticPoint {
            position_xyz_m: position,
            pressure_pa,
            pressure_potential_pa_s,
            particle_velocity_r_phi_z_m_s: velocity,
            particle_velocity_xyz_m_s: cartesian,
            active_intensity_xyz_w_m2: intensity,
        })
    }
}

#[cfg(test)]
mod tests;
