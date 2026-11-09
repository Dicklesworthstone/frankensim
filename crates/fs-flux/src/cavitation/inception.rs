//! Static cavitation inception with the bubble's original gas inventory.
//!
//! At rest the spherical models share `p_eq(R) = p_v + p_g(R) - 2 sigma/R`.
//! Its negative logarithmic slope, `-R dp_eq/dR = 3 kappa p_g - 2 sigma/R`,
//! distinguishes the mechanically stable and unstable equilibrium branches.
//! For positive gas, positive surface tension and `3 kappa > 1`, their fold
//! is the Blake threshold. Below it the nucleus has no static equilibrium.
//! Continuation also admits monotone gas-only and vapor-only limits. For
//! `3 kappa < 1`, the fold is a maximum and the stable branch has larger radii;
//! it is not a Blake minimum. A flat, marginal curve has neither strict branch.
//!
//! These operations retain the admitted liquid and polytropic gas law; they
//! never reconstruct gas with `BubbleGas::equilibrium` at a changed pressure.
//! The caller must choose the appropriate thermal limit (isothermal gas for a
//! slow pressure change). This is not a frequency-dependent acoustic threshold,
//! homogeneous nucleation model, or gas-diffusion/thermal-equilibrium model.
//!
//! Reference: C. E. Brennen, *Cavitation and Bubble Dynamics*, §2.5,
//! equations 2.41–2.46, <https://media.library.caltech.edu/CaltechBOOK:1995.001/chap2.htm>.

use super::{
    BubbleState, CavitationError, SphericalBubbleModel, finite, nonnegative, poll, positive,
    validate_state,
};

/// Pressure that supports an existing nucleus at rest, using its original gas.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StaticBubbleResponse {
    pub radius_m: f64,
    /// Can be negative when the liquid is under tension.
    pub equilibrium_absolute_pressure_pa: f64,
    pub gas_partial_pressure_pa: f64,
    /// `-R dp_eq/dR`: positive is mechanically stable, negative is unstable.
    pub restoring_pressure_pa: f64,
}

/// Fold of the static pressure-radius curve.
/// The public threshold query admits its radius against the model's domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BlakeThreshold {
    pub critical_radius_m: f64,
    /// Minimum imposed pressure admitting static equilibrium, possibly tensile.
    pub critical_absolute_pressure_pa: f64,
    pub gas_partial_pressure_pa: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EquilibriumBranch {
    Stable,
    Unstable,
}

/// Finite bisection budget and simultaneous radius/pressure residual targets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumControls {
    pub relative_radius_tolerance: f64,
    pub absolute_radius_tolerance_m: f64,
    pub absolute_pressure_tolerance_pa: f64,
    pub maximum_iterations: usize,
}

impl EquilibriumControls {
    pub(super) fn validate(&self) -> Result<(), CavitationError> {
        if !nonnegative(self.relative_radius_tolerance)
            || self.relative_radius_tolerance > 1.0
            || !positive(self.absolute_radius_tolerance_m)
            || !nonnegative(self.absolute_pressure_tolerance_pa)
            || self.maximum_iterations == 0
        {
            return Err(CavitationError::InvalidInput("static equilibrium controls"));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleEquilibrium {
    pub response: StaticBubbleResponse,
    /// Numerical sign bracket, not an interval-arithmetic certificate.
    pub radius_bracket_m: [f64; 2],
    pub pressure_residual_pa: f64,
    pub iterations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum InceptionError {
    Input(CavitationError),
    /// Positive gas, positive surface tension and `3 kappa > 1` are required.
    NoFiniteBlakeMinimum,
    /// The requested branch has no root in the admitted radius interval.
    EquilibriumOutsideDomain {
        absolute_pressure_pa: f64,
        branch: EquilibriumBranch,
        radius_bracket_m: [f64; 2],
    },
    ToleranceNotMet {
        iterations: usize,
    },
}

impl std::fmt::Display for InceptionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "cavitation inception: {self:?}")
    }
}

impl std::error::Error for InceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Input(error) => Some(error),
            _ => None,
        }
    }
}

impl From<CavitationError> for InceptionError {
    fn from(error: CavitationError) -> Self {
        Self::Input(error)
    }
}

/// Evaluate static pressure and restoring slope without changing the source.
/// Also admits ideal gas-free/zero-surface-tension cases that have no Blake fold.
pub fn static_response(
    model: &impl SphericalBubbleModel,
    radius_m: f64,
) -> Result<StaticBubbleResponse, CavitationError> {
    let limits = model.limits();
    limits.validate()?;
    let state = BubbleState::try_new(0.0, radius_m, 0.0)?;
    validate_state(&state, model.liquid(), limits)?;
    let gas = model.gas();
    let gas_pressure = gas.pressure(radius_m)?;
    if gas.reference_pressure_pa() > 0.0 && !positive(gas_pressure) {
        return Err(CavitationError::NonFiniteDerived("static gas pressure"));
    }
    let capillary = finite(
        2.0 * (model.liquid().surface_tension_n_m() / radius_m),
        "static capillary pressure",
    )?;
    Ok(StaticBubbleResponse {
        radius_m,
        equilibrium_absolute_pressure_pa: finite(
            gas_pressure + model.liquid().saturation_pressure_pa() - capillary,
            "static equilibrium pressure",
        )?,
        gas_partial_pressure_pa: gas_pressure,
        restoring_pressure_pa: finite(
            (3.0 * gas.polytropic_exponent()) * gas_pressure - capillary,
            "static restoring pressure",
        )?,
    })
}

/// Determine loss of static stability from the retained gas law.
/// A critical radius outside the native model's limits is refused explicitly.
pub fn blake_threshold(
    model: &impl SphericalBubbleModel,
) -> Result<BlakeThreshold, InceptionError> {
    let threshold = analytic_blake_threshold(model)?;
    let limits = model.limits();
    limits.validate()?;
    validate_state(
        &BubbleState::try_new(0.0, threshold.critical_radius_m, 0.0)?,
        model.liquid(),
        limits,
    )?;
    Ok(threshold)
}

// The fold may lie outside the native radius window while an equilibrium on
// one of its branches lies inside. Keep analytic geometry separate from the
// public threshold query's native-domain admission.
fn analytic_blake_threshold(
    model: &impl SphericalBubbleModel,
) -> Result<BlakeThreshold, InceptionError> {
    let gas = model.gas();
    let sigma = model.liquid().surface_tension_n_m();
    let exponent = 3.0 * gas.polytropic_exponent();
    if gas.reference_pressure_pa() == 0.0 || sigma == 0.0 || exponent <= 1.0 {
        return Err(InceptionError::NoFiniteBlakeMinimum);
    }
    let (critical_radius_m, critical_absolute_pressure_pa, gas_partial_pressure_pa) =
        analytic_turning_point(model)?;
    Ok(BlakeThreshold {
        critical_radius_m,
        critical_absolute_pressure_pa,
        gas_partial_pressure_pa,
    })
}

// Positive gas and surface tension, with 3 kappa != 1. The same stationary
// radius is a minimum for 3 kappa > 1 and a maximum for 3 kappa < 1.
// Return (radius, imposed pressure, gas pressure); only the minimum is Blake.
fn analytic_turning_point(
    model: &impl SphericalBubbleModel,
) -> Result<(f64, f64, f64), CavitationError> {
    let gas = model.gas();
    let sigma = model.liquid().surface_tension_n_m();
    let exponent = 3.0 * gas.polytropic_exponent();
    // Rc/R0 = [3 kappa pg0 R0 / (2 sigma)]^(1/(3 kappa - 1)).
    // Log coordinates avoid forming a potentially overflowing inventory R0^m.
    let log_reference_radius = gas.reference_radius_m().ln();
    let log_radius = log_reference_radius
        + (exponent.ln() + gas.reference_pressure_pa().ln() + log_reference_radius
            - std::f64::consts::LN_2
            - sigma.ln())
            / (exponent - 1.0);
    let critical_radius = log_radius.exp();
    if !positive(critical_radius) {
        return Err(CavitationError::NonFiniteDerived("static turning radius"));
    }
    let capillary = finite(
        2.0 * (sigma / critical_radius),
        "turning capillary pressure",
    )?;
    let gas_pressure = capillary / exponent;
    if !positive(gas_pressure) {
        return Err(CavitationError::NonFiniteDerived("turning gas pressure"));
    }
    Ok((
        critical_radius,
        finite(
            model.liquid().saturation_pressure_pa() - capillary * (1.0 - 1.0 / exponent),
            "static turning pressure",
        )?,
        gas_pressure,
    ))
}

/// Follow either equilibrium branch at a new imposed pressure while preserving
/// gas inventory. No root outside the radius limits is clamped into the domain.
/// The fold itself need not lie inside the admitted radius interval.
pub fn equilibrium_at_pressure(
    model: &impl SphericalBubbleModel,
    absolute_pressure_pa: f64,
    branch: EquilibriumBranch,
    controls: &EquilibriumControls,
) -> Result<BubbleEquilibrium, InceptionError> {
    equilibrium_at_pressure_controlled(model, absolute_pressure_pa, branch, controls, &mut || false)
}

/// Cancellable branch solve with deterministic bounded work and no mutation.
pub fn equilibrium_at_pressure_controlled(
    model: &impl SphericalBubbleModel,
    absolute_pressure_pa: f64,
    branch: EquilibriumBranch,
    controls: &EquilibriumControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<BubbleEquilibrium, InceptionError> {
    poll(cancelled)?;
    controls.validate()?;
    if !absolute_pressure_pa.is_finite() {
        return Err(CavitationError::InvalidInput("static equilibrium pressure").into());
    }
    let limits = model.limits();
    limits.validate()?;
    let gas = model.gas();
    let sigma = model.liquid().surface_tension_n_m();
    let exponent = 3.0 * gas.polytropic_exponent();
    let turning = if gas.reference_pressure_pa() > 0.0 && sigma > 0.0 && exponent != 1.0 {
        Some(analytic_turning_point(model)?)
    } else {
        None
    };
    let mut left = limits.minimum_radius_m;
    let mut right = limits.maximum_radius_m;
    let branch_missing = if let Some((critical_radius, _, _)) = turning {
        // A minimum has the stable branch on its left; a maximum reverses it.
        let lower_radii = (branch == EquilibriumBranch::Stable) == (exponent > 1.0);
        if lower_radii {
            let missing = critical_radius < left;
            if !missing {
                right = right.min(critical_radius);
            }
            missing
        } else {
            let missing = critical_radius > right;
            if !missing {
                left = left.max(critical_radius);
            }
            missing
        }
    } else {
        // With zero sigma, zero gas or exponent one there is no finite fold.
        // The restoring sign is constant. At exponent one compare at R0;
        // avoiding a newly authored equilibrium retains the original inventory.
        let restoring = if sigma == 0.0 {
            gas.reference_pressure_pa()
        } else if gas.reference_pressure_pa() == 0.0 {
            -1.0
        } else {
            finite(
                gas.reference_pressure_pa() - 2.0 * (sigma / gas.reference_radius_m()),
                "monotone restoring pressure",
            )?
        };
        match branch {
            EquilibriumBranch::Stable => restoring <= 0.0,
            EquilibriumBranch::Unstable => restoring >= 0.0,
        }
    };
    let unavailable = InceptionError::EquilibriumOutsideDomain {
        absolute_pressure_pa,
        branch,
        radius_bracket_m: [left, right],
    };
    let tail_branch =
        turning.is_none() || (branch == EquilibriumBranch::Stable) == (exponent < 1.0);
    if branch_missing
        || (tail_branch
            && match branch {
                EquilibriumBranch::Stable => {
                    absolute_pressure_pa <= model.liquid().saturation_pressure_pa()
                }
                EquilibriumBranch::Unstable => {
                    absolute_pressure_pa >= model.liquid().saturation_pressure_pa()
                }
            })
    {
        return Err(unavailable);
    }

    let left_response = static_response(model, left)?;
    let right_response = static_response(model, right)?;
    let at_fold = turning.is_some_and(|(radius, pressure, _)| {
        absolute_pressure_pa == pressure && radius >= left && radius <= right
    });
    let (minimum_pressure, maximum_pressure) = match branch {
        EquilibriumBranch::Stable => (
            right_response.equilibrium_absolute_pressure_pa,
            left_response.equilibrium_absolute_pressure_pa,
        ),
        EquilibriumBranch::Unstable => (
            left_response.equilibrium_absolute_pressure_pa,
            right_response.equilibrium_absolute_pressure_pa,
        ),
    };
    // Both ends matter when the native radius window truncates a branch.
    // The analytic fold has its own residual check below because its pressure
    // and the evaluated static curve can differ by floating-point roundoff.
    if !at_fold
        && (absolute_pressure_pa < minimum_pressure || absolute_pressure_pa > maximum_pressure)
    {
        return Err(unavailable);
    }
    let exact_response = if at_fold {
        Some(
            if (branch == EquilibriumBranch::Stable) == (exponent > 1.0) {
                right_response
            } else {
                left_response
            },
        )
    } else if absolute_pressure_pa == left_response.equilibrium_absolute_pressure_pa {
        Some(left_response)
    } else if absolute_pressure_pa == right_response.equilibrium_absolute_pressure_pa {
        Some(right_response)
    } else {
        None
    };
    if let Some(response) = exact_response {
        let residual = finite(
            response.equilibrium_absolute_pressure_pa - absolute_pressure_pa,
            "equilibrium pressure residual",
        )?;
        if residual.abs() <= controls.absolute_pressure_tolerance_pa {
            poll(cancelled)?;
            return Ok(BubbleEquilibrium {
                response,
                radius_bracket_m: [response.radius_m; 2],
                pressure_residual_pa: residual,
                iterations: 0,
            });
        }
        return Err(InceptionError::ToleranceNotMet { iterations: 0 });
    }

    for iterations in 1..=controls.maximum_iterations {
        poll(cancelled)?;
        let radius = left + 0.5 * (right - left);
        let response = static_response(model, radius)?;
        let residual = finite(
            response.equilibrium_absolute_pressure_pa - absolute_pressure_pa,
            "equilibrium pressure residual",
        )?;
        let radius_tolerance = finite(
            controls.absolute_radius_tolerance_m + controls.relative_radius_tolerance * radius,
            "equilibrium radius tolerance",
        )?;
        if 0.5 * (right - left) <= radius_tolerance
            && residual.abs() <= controls.absolute_pressure_tolerance_pa
        {
            poll(cancelled)?;
            return Ok(BubbleEquilibrium {
                response,
                radius_bracket_m: [left, right],
                pressure_residual_pa: residual,
                iterations,
            });
        }
        if radius == left || radius == right {
            return Err(InceptionError::ToleranceNotMet { iterations });
        }
        let move_left = match branch {
            EquilibriumBranch::Stable => residual > 0.0,
            EquilibriumBranch::Unstable => residual < 0.0,
        };
        if move_left {
            left = radius;
        } else {
            right = radius;
        }
    }
    Err(InceptionError::ToleranceNotMet {
        iterations: controls.maximum_iterations,
    })
}

#[cfg(test)]
#[path = "inception/tests.rs"]
mod tests;
