//! Planar liquid-liquid interfaces from the SAME square-gradient free energy.
//!
//! For this symmetric pair the coexistence common tangent is horizontal.
//! With `Delta g = g(x)-g(x_alpha)`, the first integral is
//! `kappa (dx/dz)^2/2 = c Delta g`, and the tension is
//! `integral[sqrt(2 kappa c Delta g)] dx` between the original binodals.
//! These are the Cahn-Hilliard planar-interface relations, with our explicit
//! `kappa/2` convention (J. Chem. Phys. 28, 258, DOI:10.1063/1.1744102).
//! No measured surface-tension provenance or sharp-interface force is invented.

use super::super::{BinaryCoexistenceOptionsV1, BinaryCoexistenceV1, BinarySolutionErrorV1};
use super::BinarySquareGradientV1;
use fs_math::det::sqrt;

/// Independent bounds for the original binodal solve and interface quadrature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryPlanarInterfaceOptionsV1 {
    pub coexistence: BinaryCoexistenceOptionsV1,
    /// Positive relative quadrature tolerance; not a rigorous error certificate.
    pub relative_tolerance: f64,
    /// Absolute quadrature tolerance [N/m], added to the relative allowance.
    pub absolute_tolerance_n_m: f64,
    /// Total integrand samples, including the two analytic endpoints; 5..=65537.
    pub maximum_integrand_evaluations: usize,
}
impl Default for BinaryPlanarInterfaceOptionsV1 {
    fn default() -> Self {
        Self {
            coexistence: BinaryCoexistenceOptionsV1::default(),
            relative_tolerance: 1e-8,
            absolute_tolerance_n_m: 0.0,
            maximum_integrand_evaluations: 4097,
        }
    }
}
impl BinaryPlanarInterfaceOptionsV1 {
    fn validate(self) -> Result<(), BinarySolutionErrorV1> {
        if !self.relative_tolerance.is_finite()
            || !(64.0 * f64::EPSILON..=1e-2).contains(&self.relative_tolerance)
            || !self.absolute_tolerance_n_m.is_finite()
            || self.absolute_tolerance_n_m < 0.0
            || !(5..=65537).contains(&self.maximum_integrand_evaluations)
        {
            return Err(BinarySolutionErrorV1::Invalid("planar interface controls"));
        }
        Ok(())
    }
}

/// An equilibrium planar interface and bounded numerical receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryPlanarInterfaceV1 {
    model: BinarySquareGradientV1,
    coexistence: BinaryCoexistenceV1,
    options: BinaryPlanarInterfaceOptionsV1,
    surface_tension_n_m: f64,
    central_slope_width_m: f64,
    center_excess_gibbs_j_per_mol: f64,
    estimated_absolute_error_n_m: f64,
    integrand_evaluations: usize,
}
impl BinaryPlanarInterfaceV1 {
    /// Original pair, temperature, site concentration and independently supplied kappa.
    pub const fn model(&self) -> &BinarySquareGradientV1 {
        &self.model
    }
    pub const fn coexistence(&self) -> &BinaryCoexistenceV1 {
        &self.coexistence
    }
    pub const fn options(&self) -> BinaryPlanarInterfaceOptionsV1 {
        self.options
    }
    /// Interfacial excess free energy per area [N/m = J/m^2].
    pub const fn surface_tension_n_m(&self) -> f64 {
        self.surface_tension_n_m
    }
    /// Composition gap divided by the maximum (central) profile slope [m].
    /// This finite characteristic width is NOT the full profile extent: the
    /// exact binodals are approached only in infinite spatial tails.
    pub const fn central_slope_width_m(&self) -> f64 {
        self.central_slope_width_m
    }
    pub const fn center_excess_gibbs_j_per_mol(&self) -> f64 {
        self.center_excess_gibbs_j_per_mol
    }
    /// Summed adaptive-Simpson truncation estimates [N/m]. Excludes floating
    /// point roundoff and binodal uncertainty; see coexistence.log_odds_bracket.
    pub const fn estimated_absolute_error_n_m(&self) -> f64 {
        self.estimated_absolute_error_n_m
    }
    pub const fn integrand_evaluations(&self) -> usize {
        self.integrand_evaluations
    }
}

impl BinarySquareGradientV1 {
    /// Compute a planar equilibrium interface using the original coexistence
    /// states and homogeneous mixing law. Miscible and critical states return
    /// `None`; unresolved positive barriers, quadrature work exhaustion and
    /// cancellation refuse. Negative excess energy is never clipped to zero.
    pub fn equilibrium_interface(
        &self,
        options: BinaryPlanarInterfaceOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<BinaryPlanarInterfaceV1>, BinarySolutionErrorV1> {
        options.validate()?;
        let Some(coexistence) = self.thermodynamics().coexistence(
            self.temperature(),
            options.coexistence,
            &mut cancelled,
        )?
        else {
            return Ok(None);
        };
        let lower = coexistence.lower_a().mole_fractions();
        let half_gap = 0.5 - lower[0];
        let gap = 2.0 * half_gap;
        let barrier = positive(
            self.mixing_free_energy_change_j_per_mol(lower, [0.5, 0.5])?,
            "unresolved planar interface barrier",
        )?;
        // Scale before quadrature: its dimensionless shape stays O(1) as the
        // critical barrier vanishes. Separate square roots avoid c*Delta-g
        // overflow when both the final tension and width are representable.
        let root_kappa = sqrt(self.gradient_coefficient_j_m());
        let root_c = sqrt(self.molar_concentration_mol_m3());
        let root_barrier = sqrt(barrier);
        let scale = positive(
            gap * root_kappa * root_c * root_barrier * std::f64::consts::SQRT_2,
            "unrepresentable planar interface tension scale",
        )?;
        let width = positive(
            gap * root_kappa / root_c / root_barrier / std::f64::consts::SQRT_2,
            "unrepresentable planar interface width",
        )?;
        let mut samples = 0;
        let mut integrand = |t: f64| -> Result<f64, BinarySolutionErrorV1> {
            if cancelled() {
                return Err(BinarySolutionErrorV1::Cancelled);
            }
            if samples == options.maximum_integrand_evaluations {
                return Err(BinarySolutionErrorV1::IterationLimit);
            }
            samples += 1;
            // Exact boundary values from the definition of Delta g, not floors.
            if t == 0.0 || t == 1.0 {
                return Ok(t);
            }
            let x = lower[0] + half_gap * t;
            if x <= lower[0] || x >= 0.5 {
                return Err(BinarySolutionErrorV1::Unrepresentable(
                    "unresolved interior interface composition",
                ));
            }
            let excess = positive(
                self.mixing_free_energy_change_j_per_mol(lower, [x, 1.0 - x])?,
                "nonpositive interior common-tangent excess energy",
            )?;
            positive(sqrt(excess / barrier), "unresolved interface integrand")
        };
        let (integral, error) = integrate(
            &mut integrand,
            options.relative_tolerance,
            options.absolute_tolerance_n_m / scale,
        )?;
        if cancelled() {
            return Err(BinarySolutionErrorV1::Cancelled);
        }
        Ok(Some(BinaryPlanarInterfaceV1 {
            model: self.clone(),
            coexistence,
            options,
            surface_tension_n_m: positive(scale * integral, "planar interface tension")?,
            central_slope_width_m: width,
            center_excess_gibbs_j_per_mol: barrier,
            estimated_absolute_error_n_m: scale * error,
            integrand_evaluations: samples,
        }))
    }
}

fn positive(value: f64, label: &'static str) -> Result<f64, BinarySolutionErrorV1> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(BinarySolutionErrorV1::Unrepresentable(label))
    }
}

#[derive(Clone, Copy)]
struct Panel {
    a: f64,
    b: f64,
    f: [f64; 3],
    simpson: f64,
}
impl Panel {
    fn new(a: f64, b: f64, f: [f64; 3]) -> Self {
        Self {
            a,
            b,
            f,
            simpson: (b - a) * (f[0] + 4.0 * f[1] + f[2]) / 6.0,
        }
    }
}

fn integrate(
    sample: &mut impl FnMut(f64) -> Result<f64, BinarySolutionErrorV1>,
    relative: f64,
    absolute: f64,
) -> Result<(f64, f64), BinarySolutionErrorV1> {
    let mut stack = Vec::with_capacity(64);
    stack.push(Panel::new(
        0.0,
        1.0,
        [sample(0.0)?, sample(0.5)?, sample(1.0)?],
    ));
    // The common-tangent excess increases from the lower binodal to the
    // center. Therefore half the mid-sample is a conservative integral scale.
    // Distribute one global budget by panel width: a vanishing endpoint must
    // not demand arbitrarily small *local* relative errors from cancellation.
    let budget = absolute + relative * (0.5 * stack[0].f[1]);
    let (mut integral, mut error) = (0.0, 0.0);
    while let Some(panel) = stack.pop() {
        let mid = panel.a + 0.5 * (panel.b - panel.a);
        let left_mid = panel.a + 0.5 * (mid - panel.a);
        let right_mid = mid + 0.5 * (panel.b - mid);
        if left_mid <= panel.a || right_mid >= panel.b {
            return Err(BinarySolutionErrorV1::Unrepresentable(
                "interface quadrature resolution",
            ));
        }
        let left = Panel::new(panel.a, mid, [panel.f[0], sample(left_mid)?, panel.f[1]]);
        let right = Panel::new(mid, panel.b, [panel.f[1], sample(right_mid)?, panel.f[2]]);
        let refined = left.simpson + right.simpson;
        let correction = (refined - panel.simpson) / 15.0;
        let estimate = positive(refined + correction, "interface quadrature estimate")?;
        let local_error = correction.abs();
        if local_error <= budget * (panel.b - panel.a) {
            integral += estimate;
            error += local_error;
        } else {
            stack.push(right);
            stack.push(left);
        }
    }
    if error > absolute + relative * integral {
        return Err(BinarySolutionErrorV1::Unrepresentable(
            "unresolved global interface quadrature accuracy",
        ));
    }
    Ok((integral, error))
}

#[cfg(test)]
mod tests;
