//! Signed circular mouth modes through the original retarded Rayleigh owner.
//!
//! With the duct's area-orthonormal basis, v(r,t)=sum U_n(t) psi_n(r/R)/S.
//! Higher U_n have zero net flow but radiate through their spatial phase and
//! distance variation. They are neither independent pistons nor silent ports.
//! Gauss-Legendre radial and periodic angular quadrature sample that original
//! trace; no discrete mean subtraction or basis renormalization is performed.
//!
//! This is a stationary circular opening in an infinite rigid baffle. It is an
//! observation of prescribed flow, not a new reaction law. A coupled consumer
//! must supply the same aperture, exterior gas and full flanged load used by
//! its source. An unflanged/BEM load cannot acquire this exterior implicitly.

use super::{BaffledPressure, CircularOutletReceiver, RayleighMedium, SurfaceSample};
use fs_duct::modal::{AxisymmetricMouthBasis, ModeCount};
use fs_feec::highorder::quad1d::gauss_legendre;
use fs_material::gas::GasState;
use fs_math::det;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

/// Total cold quadrature-work ceiling, independent of the source FIR budget.
pub const MAX_MODAL_KERNEL_TERMS: usize = 1 << 27;

/// Caller-declared spatial convergence and complete construction-work budget.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalOutletOptions {
    /// Final half-to-full radial AND angular kernel changes must be below
    /// this value for every mode. Changes use the unit-flow radiation scale
    /// rho/(2 pi d_min), not a possibly vanishing observed modal pressure.
    /// This measures quadrature change, not a continuum error enclosure.
    pub relative_kernel_tolerance: f64,
    /// Sum over all five axis rungs of (67 * points * modes + radial * modes²).
    /// The 67 terms include the original linear and 65-tap sinc preparations.
    pub maximum_kernel_terms: usize,
}

/// Convergence of the actual retarded kernels supplied to the running receiver.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalOutletReport {
    /// Original retained basis count, including the uniform mode.
    pub modes: usize,
    /// Quarter, half and full radial node counts; angles are held at full.
    pub radial_nodes: [usize; 3],
    /// Quarter, half and full angle counts; radial nodes are held at full.
    pub angular_points: [usize; 3],
    /// Per-mode kernel L1 changes: quarter-to-half, then half-to-full.
    pub radial_changes: [Vec<f64>; 2],
    /// Same two changes on the independent angular axis.
    pub angular_changes: [Vec<f64>; 2],
    /// Infinity norm of the fine discrete modal Gram matrix minus identity.
    pub basis_gram_defect: f64,
    /// Complete five-rung work admitted before allocation.
    pub kernel_terms: usize,
    /// Original requested tolerance, retained without loosening on failure.
    pub tolerance: f64,
}

impl BaffledPressure {
    /// Observe ALL signed mouth-flow coordinates at the declared finite point.
    /// Radial/angle counts are multiples of four in 4..=128 and 32..=512.
    /// Here `radial_rings` means Gauss-Legendre radial nodes. Both axes are
    /// checked on quarter/half/full rungs, while the other stays at full.
    ///
    /// The final kernel is the existing bandlimited derivative and causal
    /// travel-time image, with its original clock/history limits and reported
    /// numerical latency. `step_trace` takes complete interleaved modal flows.
    /// A kernel L1 difference bounds the change for arbitrary bounded modal
    /// acceleration histories; its normalization does not divide by a null.
    /// This is not a bound on omitted modes, nonlinear bandwidth, continuum
    /// quadrature error or finite-baffle/body scattering.
    ///
    /// # Errors
    /// Original gas/basis/receiver refusals, invalid declared budgets, spatial
    /// underresolution, unconverged modal kernels, or cancellation between rungs.
    #[allow(clippy::too_many_arguments)]
    pub fn circular_modal_outlet_multirate_differentiating(
        radius_m: f64, modes: ModeCount, mechanics_rate: u32, output_rate: u32,
        spec: CircularOutletReceiver, exterior: &GasState, options: ModalOutletOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, ModalOutletReport), String> {
        exterior.validate_coherence().map_err(|e| e.to_string())?;
        let basis = AxisymmetricMouthBasis::new(modes).map_err(|e| e.to_string())?;
        let (nr, na, n) = (spec.radial_rings, spec.angular_points, modes.0);
        let tolerance = options.relative_kernel_tolerance;
        if !(4..=128).contains(&nr) || nr % 4 != 0
            || !(32..=512).contains(&na) || na % 4 != 0
            || !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0
            || !(1..=MAX_MODAL_KERNEL_TERMS).contains(&options.maximum_kernel_terms)
            || !radius_m.is_finite() || radius_m <= 0.0
            || !spec.maximum_frequency_hz.is_finite() || spec.maximum_frequency_hz <= 0.0
            || spec.maximum_frequency_hz > 0.45 * f64::from(output_rate)
            || spec.position_m.iter().any(|v| !v.is_finite()) || spec.position_m[2] < 0.05
        {
            return Err("modal outlet needs bounded quarter/half/full quadrature, positive tolerance and admitted band/clearance".into());
        }
        let rungs = [(nr, na), (nr/2, na), (nr/4, na), (nr, na/2), (nr, na/4)];
        let kernel_terms = rungs.iter().map(|&(r,a)| 67*r*a*n + r*n*n).sum::<usize>();
        if kernel_terms > options.maximum_kernel_terms {
            return Err("complete modal outlet quadrature exceeds its declared kernel-work budget".into());
        }
        let medium = RayleighMedium::from(exterior);
        let minimum_distance = det::hypot(spec.position_m[2],
            det::hypot(spec.position_m[0], spec.position_m[1]).max(radius_m) - radius_m);
        let scale = exterior.density/(TAU*minimum_distance);
        if !scale.is_finite() || scale <= 0.0 {
            return Err("modal outlet radiation scale is unrepresentable".into());
        }
        let mut prepare = |r: usize, a: usize| -> Result<(Self, f64, f64), String> {
            if cancelled() { return Err("modal outlet preparation cancelled".into()); }
            let (surface, gram, spacing) = surface(radius_m, &basis, r, a)?;
            let receiver = Self::new_multirate_differentiating(&surface, n,
                mechanics_rate, output_rate, spec.position_m, medium)?;
            Ok((receiver, gram, spacing))
        };
        let (fine, basis_gram_defect, spacing) = prepare(nr, na)?;
        let phase = TAU*spec.maximum_frequency_hz*radius_m/exterior.sound_speed
            * (spacing + TAU/na as f64);
        let radial_phase = basis.radial_eigenvalues()[n-1]*spacing;
        if !phase.is_finite() || phase > FRAC_PI_4 || radial_phase > FRAC_PI_2
            || !basis_gram_defect.is_finite() || basis_gram_defect > tolerance
        {
            return Err("modal outlet quadrature underresolves the original basis or declared radiation band; refine it".into());
        }
        let mut changes = |radial: bool| -> Result<[Vec<f64>; 2], String> {
            let (half, _, _) = prepare(if radial {nr/2} else {nr}, if radial {na} else {na/2})?;
            let last = kernel_change(&fine, &half, scale);
            let (quarter, _, _) = prepare(if radial {nr/4} else {nr}, if radial {na} else {na/4})?;
            let first = kernel_change(&half, &quarter, scale);
            if last.iter().any(|v| !v.is_finite() || *v > tolerance)
                || first.iter().any(|v| !v.is_finite())
            {
                return Err(format!("modal outlet {} quadrature did not converge: final per-mode kernel changes {last:?}",
                    if radial {"radial"} else {"angular"}));
            }
            Ok([first, last])
        };
        let radial_changes = changes(true)?;
        let angular_changes = changes(false)?;
        if cancelled() { return Err("modal outlet preparation cancelled".into()); }
        Ok((fine, ModalOutletReport { modes: n, radial_nodes: [nr/4,nr/2,nr],
            angular_points: [na/4,na/2,na], radial_changes, angular_changes,
            basis_gram_defect, kernel_terms, tolerance }))
    }
}

fn surface(radius: f64, basis: &AxisymmetricMouthBasis, radial: usize, angular: usize)
    -> Result<(Vec<SurfaceSample>, f64, f64), String>
{
    let area = PI*radius*radius;
    if !area.is_finite() || area <= 0.0 || !area.recip().is_finite() {
        return Err("modal outlet area is unrepresentable".into());
    }
    let (nodes, weights) = gauss_legendre(radial);
    let n = basis.modes().0;
    let mut gram = vec![0.0; n*n];
    let mut samples = Vec::with_capacity(radial*angular);
    let (mut previous, mut spacing) = (0.0_f64, 0.0_f64);
    for (&node, &weight) in nodes.iter().zip(&weights) {
        let eta = 0.5*(node+1.0);
        spacing = spacing.max(eta-previous);
        previous = eta;
        let shapes = (0..n).map(|m| basis.shape(m, eta).map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        for i in 0..n { for j in 0..n { gram[i*n+j] += weight*eta*shapes[i]*shapes[j]; } }
        let velocity: Vec<_> = shapes.iter().map(|v| v/area).collect();
        for angle in 0..angular {
            let theta = TAU*(angle as f64+0.5)/angular as f64;
            samples.push(SurfaceSample { position_m: [radius*eta*det::cos(theta), radius*eta*det::sin(theta), 0.0],
                area_m2: area*eta*weight/angular as f64, mode_shape: velocity.clone() });
        }
    }
    spacing = spacing.max(1.0-previous);
    let defect = (0..n).map(|i| (0..n).map(|j|
        (gram[i*n+j]-if i == j {1.0} else {0.0}).abs()).sum::<f64>()).fold(0.0, f64::max);
    Ok((samples, defect, spacing))
}

// Both kernels use the same physical origin and numerical delay. Missing rows
// are zero. Sum absolute differences before normalizing, so signed cancellation
// cannot hide an unresolved lobe or a quiet higher-mode contribution.
fn kernel_change(a: &BaffledPressure, b: &BaffledPressure, scale: f64) -> Vec<f64> {
    let mut change = vec![0.0; a.modes];
    let (mut i, mut j) = (0, 0);
    while i < a.kernels.len() || j < b.kernels.len() {
        let da = a.kernels.get(i).map_or(usize::MAX, |row| row.0);
        let db = b.kernels.get(j).map_or(usize::MAX, |row| row.0);
        for (m, error) in change.iter_mut().enumerate() {
            let left = if da <= db {a.kernels[i].1[m]} else {0.0};
            let right = if db <= da {b.kernels[j].1[m]} else {0.0};
            *error += (left-right).abs()/scale;
        }
        if da <= db { i += 1; }
        if db <= da { j += 1; }
    }
    change
}

#[cfg(test)]
#[path = "baffled_pressure_modal_tests.rs"]
mod tests;
