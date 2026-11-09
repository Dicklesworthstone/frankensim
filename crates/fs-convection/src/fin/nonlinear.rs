//! Quasi-steady straight fins with original local material and film laws.
//!
//! Discretizes `d/dx(k(T) A dT/dx) = P h(x,T) (T-T_fluid)` with nodal
//! control volumes, half volumes at the prescribed base and the exposed tip.
//! Each axial face uses the caller's mean `integral(k(T),T_left,T_right) /
//! (T_right-T_left)`, including its equal-temperature limit. That is the exact
//! Kirchhoff flux for the source-free link; no sampled arithmetic/harmonic
//! average replaces the original conductivity law. Lateral loss uses the
//! local nodal film and trapezoidal quadrature.
//!
//! The conservative element balance follows MIT Unified Thermodynamics 18.2,
//! equation (18.7): https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node128.html.
//! Keeping Fourier flux in divergence form extends that balance to k(T).
//! Positive frozen coefficients give an O(N) admittance sweep; bounded relaxed
//! Picard iteration resolves them at the actual temperatures. Convergence is
//! gated on both thermal correction and independently assembled node/global
//! heat defects. It does not establish uniqueness for an arbitrary callback.
//!
//! This is a one-dimensional, quasi-steady discretization with one uniform
//! fluid temperature. It adds no fin storage, radiation, contact resistance,
//! fluid heating, fin-array interference or correlation applicability claim.

use super::{FinError, StraightFinGeometry, admit_temperatures};

/// The lateral and tip films may be distinct original constitutive laws.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FinSurface {
    Side,
    Tip,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinPoint {
    pub position_m: f64,
    pub temperature_k: f64,
    pub surface: FinSurface,
}

/// Values at one actual local state, with the original owner's source payload.
#[derive(Debug, PartialEq)]
pub struct FinProperties<S> {
    pub conductivity_w_m_k: f64,
    pub coefficient_w_m2_k: f64,
    pub source: S,
}

#[derive(Debug, PartialEq)]
pub struct FinNode<S> {
    pub point: FinPoint,
    pub properties: FinProperties<S>,
}

/// `Convection` queries a separate `FinSurface::Tip` at the actual tip state.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NonlinearFinTip {
    Adiabatic,
    Convection,
}

/// All discretization, nonlinear, source-work and model limits are authored.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonlinearFinControls {
    pub segments: usize,
    pub maximum_segments: usize,
    pub maximum_iterations: usize,
    /// Includes every node, tip and axial-mean query, even discarded iterates.
    pub maximum_source_evaluations: usize,
    pub absolute_heat_tolerance_w: f64,
    pub relative_heat_tolerance: f64,
    pub absolute_temperature_tolerance_k: f64,
    /// In (0,1]; no automatic under-relaxation or hidden nonlinear fallback.
    pub relaxation: f64,
    /// Gates h A/(k P) at every actual local side/tip source evaluation.
    /// This is a sampled local criterion, not a bound between queried states.
    pub maximum_transverse_biot: f64,
}

impl NonlinearFinControls {
    pub fn validate(self) -> Result<(), &'static str> {
        if self.segments == 0
            || self.segments > self.maximum_segments
            || self.segments == usize::MAX
            || self.maximum_iterations == 0
            || self.maximum_source_evaluations == 0
        {
            return Err("positive bounded fin mesh, iterations and source work required");
        }
        if [
            self.absolute_heat_tolerance_w,
            self.absolute_temperature_tolerance_k,
            self.relaxation,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.0)
            || self.relaxation > 1.0
            || !self.relative_heat_tolerance.is_finite()
            || self.relative_heat_tolerance < 0.0
            || !self.maximum_transverse_biot.is_finite()
            || self.maximum_transverse_biot < 0.0
        {
            return Err(
                "finite positive fin tolerances, relaxation in (0,1], and nonnegative relative/Biot limits required",
            );
        }
        Ok(())
    }
}

/// A converged discrete endpoint and its original local source evaluations.
#[derive(Debug, PartialEq)]
pub struct NonlinearFinResponse<S> {
    /// Ordered side nodes, including x=0 and x=L; no source payload is cloned.
    pub nodes: Vec<FinNode<S>>,
    pub tip: Option<FinNode<S>>,
    /// Positive heat leaves the base and enters the fluid.
    pub base_heat_w: f64,
    pub side_heat_w: f64,
    pub tip_heat_w: f64,
    /// Base heat minus independently summed side and tip heat.
    pub energy_residual_w: f64,
    /// Sum of absolute non-Dirichlet control-volume heat defects.
    pub sum_absolute_node_residual_w: f64,
    /// Q_base/(T_base-T_fluid), or the frozen local limiting conductance at
    /// exact thermal equilibrium. Equilibrium still evaluates all sources.
    pub effective_conductance_w_k: f64,
    pub maximum_transverse_biot: f64,
    pub iterations: usize,
    pub source_evaluations: usize,
}

#[derive(Debug, PartialEq)]
pub enum NonlinearFinError<E> {
    Input(FinError),
    Controls(&'static str),
    Allocation,
    Arithmetic,
    Cancelled,
    WorkLimit {
        source_evaluations: usize,
    },
    Source {
        position_m: f64,
        error: E,
    },
    InvalidProperties {
        position_m: f64,
    },
    TransverseBiot {
        position_m: f64,
        biot: f64,
        maximum: f64,
    },
    NotConverged {
        iterations: usize,
        source_evaluations: usize,
        maximum_temperature_correction_k: f64,
        sum_absolute_node_residual_w: f64,
        energy_residual_w: f64,
    },
}

impl<E: core::fmt::Display> core::fmt::Display for NonlinearFinError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source { position_m, error } => {
                write!(f, "nonlinear fin source at x={position_m} m: {error}")
            }
            Self::Controls(reason) => write!(f, "nonlinear fin controls: {reason}"),
            Self::Input(error) => write!(f, "nonlinear fin input: {error}"),
            Self::Allocation => write!(f, "nonlinear fin allocation refused"),
            Self::Arithmetic => write!(f, "nonlinear fin arithmetic is unrepresentable"),
            Self::Cancelled => write!(f, "nonlinear fin cancelled"),
            Self::WorkLimit { source_evaluations } => {
                write!(
                    f,
                    "nonlinear fin source work exhausted after {source_evaluations} evaluations"
                )
            }
            Self::InvalidProperties { position_m } => {
                write!(
                    f,
                    "nonlinear fin needs positive k and nonnegative h at x={position_m} m"
                )
            }
            Self::TransverseBiot {
                position_m,
                biot,
                maximum,
            } => {
                write!(
                    f,
                    "nonlinear fin transverse Bi={biot} exceeds {maximum} at x={position_m} m"
                )
            }
            Self::NotConverged {
                iterations,
                sum_absolute_node_residual_w,
                ..
            } => {
                write!(
                    f,
                    "nonlinear fin did not converge in {iterations} iterations; absolute node defect={sum_absolute_node_residual_w} W"
                )
            }
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for NonlinearFinError<E> {}

/// Resolve original local laws throughout a conservative fin temperature field.
///
/// `axial_mean_k(a,b)` must return the positive integral mean of the SAME
/// temperature-only conductivity law used in `local`, including all interior
/// material knots and its point value when a==b. It is queried once per link
/// and iteration, with no sign assumption on a-b. Local film laws may depend
/// on x and T; their physical applicability remains with the source owner.
/// The tip is a separate source query; no tip coefficient is inferred.
///
/// Accepted temperatures remain between the base and fluid temperatures.
/// Every callback and each bounded linear sweep polls cancellation. Invalid
/// source data and source-domain failures are terminal; no clipped/extrapolated
/// law replaces them. All work and allocation stays within explicit controls.
/// The residual tolerances bound the discrete nonlinear solution, not axial
/// discretization error; use mesh refinement for the latter.
#[allow(clippy::too_many_arguments)]
pub fn solve<S, E>(
    geometry: StraightFinGeometry,
    base_temperature_k: f64,
    fluid_temperature_k: f64,
    tip: NonlinearFinTip,
    controls: NonlinearFinControls,
    mut local: impl FnMut(FinPoint) -> Result<FinProperties<S>, E>,
    mut axial_mean_k: impl FnMut(f64, f64) -> Result<f64, E>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<NonlinearFinResponse<S>, NonlinearFinError<E>> {
    controls.validate().map_err(NonlinearFinError::Controls)?;
    admit_temperatures(base_temperature_k, fluid_temperature_k)
        .map_err(NonlinearFinError::Input)?;
    check_cancel(&mut cancelled)?;
    let n = controls.segments;
    let delta = finite(base_temperature_k - fluid_temperature_k)?;
    let dx = geometry.length_m() / n as f64;
    let axial_scale = geometry.cross_section_m2() / dx;
    let side_area = geometry.exposed_perimeter_m() * dx;
    if [dx, axial_scale, side_area]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
    {
        return Err(NonlinearFinError::Arithmetic);
    }
    let mut ratio = allocated(n + 1, 1.0)?;
    let mut proposed = allocated(n + 1, 1.0)?;
    let mut axial = allocated(n, 0.0)?;
    let mut side = allocated(n + 1, 0.0)?;
    let mut downstream = allocated(n + 1, 0.0)?;
    let mut nodes = Vec::new();
    nodes
        .try_reserve_exact(n + 1)
        .map_err(|_| NonlinearFinError::Allocation)?;
    let mut work = 0;

    for iteration in 1..=controls.maximum_iterations {
        nodes.clear();
        let mut maximum_biot: f64 = 0.0;
        for i in 0..=n {
            let point = FinPoint {
                position_m: if i == n {
                    geometry.length_m()
                } else {
                    i as f64 * dx
                },
                temperature_k: if i == 0 {
                    base_temperature_k
                } else {
                    ratio[i].mul_add(delta, fluid_temperature_k)
                },
                surface: FinSurface::Side,
            };
            meter(&mut work, controls, &mut cancelled)?;
            let properties = local(point).map_err(|error| NonlinearFinError::Source {
                position_m: point.position_m,
                error,
            })?;
            maximum_biot =
                maximum_biot.max(admit_properties(geometry, point, &properties, controls)?);
            side[i] = finite(
                properties.coefficient_w_m2_k
                    * side_area
                    * if i == 0 || i == n { 0.5 } else { 1.0 },
            )?;
            nodes.push(FinNode { point, properties });
        }
        let tip_node = if tip == NonlinearFinTip::Convection {
            let point = FinPoint {
                surface: FinSurface::Tip,
                ..nodes[n].point
            };
            meter(&mut work, controls, &mut cancelled)?;
            let properties = local(point).map_err(|error| NonlinearFinError::Source {
                position_m: point.position_m,
                error,
            })?;
            maximum_biot =
                maximum_biot.max(admit_properties(geometry, point, &properties, controls)?);
            Some(FinNode { point, properties })
        } else {
            None
        };
        let tip_load = finite(tip_node.as_ref().map_or(0.0, |node| {
            node.properties.coefficient_w_m2_k * geometry.cross_section_m2()
        }))?;
        for i in 0..n {
            meter(&mut work, controls, &mut cancelled)?;
            let k = axial_mean_k(
                nodes[i].point.temperature_k,
                nodes[i + 1].point.temperature_k,
            )
            .map_err(|error| NonlinearFinError::Source {
                position_m: (i as f64 + 0.5) * dx,
                error,
            })?;
            if !k.is_finite() || k <= 0.0 {
                return Err(NonlinearFinError::InvalidProperties {
                    position_m: (i as f64 + 0.5) * dx,
                });
            }
            axial[i] = finite(k * axial_scale)?;
            if axial[i] <= 0.0 {
                return Err(NonlinearFinError::Arithmetic);
            }
        }

        // Eliminate from the tip using only positive conductances. Unlike a
        // subtractive Thomas pivot, this preserves the exact insulated limit
        // and keeps every solved ratio in [0,1] without clipping a temperature.
        downstream[n] = finite(side[n] + tip_load)?;
        for i in (1..n).rev() {
            if i % 64 == 0 {
                check_cancel(&mut cancelled)?;
            }
            downstream[i] = finite(side[i] + series(axial[i], downstream[i + 1]))?;
        }
        let limiting_conductance = finite(side[0] + series(axial[0], downstream[1]))?;
        let mut correction_k: f64 = 0.0;
        for i in 1..=n {
            if i % 64 == 0 {
                check_cancel(&mut cancelled)?;
            }
            proposed[i] = proposed[i - 1] * transmission(axial[i - 1], downstream[i]);
            correction_k = correction_k.max(finite((proposed[i] - ratio[i]).abs() * delta.abs())?);
        }

        // Recompute the ORIGINAL nonlinear finite-volume balances at the
        // current source-evaluated temperatures, not at the frozen proposal.
        // Kelvin nodes are the public state and the original source's input.
        // Use their actual differences in the accepted fluxes: a hidden ratio
        // cannot claim heat through a drop too small to represent in that state.
        let base_heat = finite(
            axial[0] * (nodes[0].point.temperature_k - nodes[1].point.temperature_k)
                + side[0] * delta,
        )?;
        let tip_heat = finite(tip_load * (nodes[n].point.temperature_k - fluid_temperature_k))?;
        let mut side_sum = Sum::default();
        let mut defect_sum = Sum::default();
        side_sum.add(finite(side[0] * delta)?);
        for i in 1..=n {
            if i % 64 == 0 {
                check_cancel(&mut cancelled)?;
            }
            let lateral = finite(side[i] * (nodes[i].point.temperature_k - fluid_temperature_k))?;
            side_sum.add(lateral);
            let incoming = finite(
                axial[i - 1] * (nodes[i - 1].point.temperature_k - nodes[i].point.temperature_k),
            )?;
            let outgoing = if i == n {
                tip_heat
            } else {
                finite(
                    axial[i] * (nodes[i].point.temperature_k - nodes[i + 1].point.temperature_k),
                )?
            };
            defect_sum.add(finite(incoming - outgoing - lateral)?.abs());
        }
        let side_heat = finite(side_sum.get())?;
        let residual = finite(base_heat - side_heat - tip_heat)?;
        let defect = finite(defect_sum.get())?;
        let heat_scale = base_heat
            .abs()
            .max(finite(side_heat.abs() + tip_heat.abs())?);
        let tolerance = finite(
            controls.absolute_heat_tolerance_w + controls.relative_heat_tolerance * heat_scale,
        )?;
        if correction_k <= controls.absolute_temperature_tolerance_k
            && defect <= tolerance
            && residual.abs() <= tolerance
        {
            let conductance = if delta == 0.0 {
                limiting_conductance
            } else {
                base_heat / delta
            };
            if !conductance.is_finite() || conductance < 0.0 {
                return Err(NonlinearFinError::Arithmetic);
            }
            check_cancel(&mut cancelled)?;
            return Ok(NonlinearFinResponse {
                nodes,
                tip: tip_node,
                base_heat_w: base_heat,
                side_heat_w: side_heat,
                tip_heat_w: tip_heat,
                energy_residual_w: residual,
                sum_absolute_node_residual_w: defect,
                effective_conductance_w_k: conductance,
                maximum_transverse_biot: maximum_biot,
                iterations: iteration,
                source_evaluations: work,
            });
        }
        if iteration == controls.maximum_iterations {
            return Err(NonlinearFinError::NotConverged {
                iterations: iteration,
                source_evaluations: work,
                maximum_temperature_correction_k: correction_k,
                sum_absolute_node_residual_w: defect,
                energy_residual_w: residual,
            });
        }
        for i in 1..=n {
            if i % 64 == 0 {
                check_cancel(&mut cancelled)?;
            }
            ratio[i] = (1.0 - controls.relaxation) * ratio[i] + controls.relaxation * proposed[i];
        }
    }
    unreachable!("positive bounded iteration count was admitted")
}

fn admit_properties<S, E>(
    geometry: StraightFinGeometry,
    point: FinPoint,
    properties: &FinProperties<S>,
    controls: NonlinearFinControls,
) -> Result<f64, NonlinearFinError<E>> {
    if !point.temperature_k.is_finite()
        || point.temperature_k <= 0.0
        || !properties.conductivity_w_m_k.is_finite()
        || properties.conductivity_w_m_k <= 0.0
        || !properties.coefficient_w_m2_k.is_finite()
        || properties.coefficient_w_m2_k < 0.0
    {
        return Err(NonlinearFinError::InvalidProperties {
            position_m: point.position_m,
        });
    }
    let biot = finite(
        (properties.coefficient_w_m2_k / properties.conductivity_w_m_k)
            * (geometry.cross_section_m2() / geometry.exposed_perimeter_m()),
    )?;
    if biot > controls.maximum_transverse_biot {
        return Err(NonlinearFinError::TransverseBiot {
            position_m: point.position_m,
            biot,
            maximum: controls.maximum_transverse_biot,
        });
    }
    Ok(biot)
}

fn finite<E>(x: f64) -> Result<f64, NonlinearFinError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(NonlinearFinError::Arithmetic)
    }
}
fn allocated<E>(length: usize, value: f64) -> Result<Vec<f64>, NonlinearFinError<E>> {
    let mut out = Vec::new();
    out.try_reserve_exact(length)
        .map_err(|_| NonlinearFinError::Allocation)?;
    out.resize(length, value);
    Ok(out)
}
fn check_cancel<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), NonlinearFinError<E>> {
    if cancelled() {
        Err(NonlinearFinError::Cancelled)
    } else {
        Ok(())
    }
}
fn meter<E>(
    work: &mut usize,
    c: NonlinearFinControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), NonlinearFinError<E>> {
    check_cancel(cancelled)?;
    if *work == c.maximum_source_evaluations {
        return Err(NonlinearFinError::WorkLimit {
            source_evaluations: *work,
        });
    }
    *work += 1;
    Ok(())
}
fn series(a: f64, b: f64) -> f64 {
    if b == 0.0 {
        0.0
    } else {
        a.min(b) / (1.0 + a.min(b) / a.max(b))
    }
}
fn transmission(axial: f64, load: f64) -> f64 {
    if axial >= load {
        1.0 / (1.0 + load / axial)
    } else {
        let ratio = axial / load;
        ratio / (1.0 + ratio)
    }
}
#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let next = self.value + x;
        self.correction += if self.value.abs() >= x.abs() {
            (self.value - next) + x
        } else {
            (x - next) + self.value
        };
        self.value = next;
    }
    fn get(&self) -> f64 {
        self.value + self.correction
    }
}

#[cfg(test)]
mod tests;
