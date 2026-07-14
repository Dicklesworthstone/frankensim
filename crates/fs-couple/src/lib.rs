//! fs-couple — multiphysics composition through port-Hamiltonian Dirac
//! structures. Layer: L3.
//!
//! Ad-hoc FSI staggering suffers added-mass instabilities and energy drift.
//! The implemented DIRAC interconnection is LOSSLESS BY CONSTRUCTION:
//! power-conjugate [`Port`]s use equal effort and opposite flow, so their net
//! interface power is exactly zero. [`EnergyAudit`] records caller-supplied
//! interface balances as a G0 bug alarm. Neither invariant alone proves that
//! the coupled components, discretizations, transfers, iterations, time
//! integrators, sources, or a finite accounting window are passive.
//!
//! For the hard, strongly-coupled cases, [`AitkenRelaxation`] gives dynamic
//! interface relaxation: on the classic ADDED-MASS-INSTABILITY fixture (a light
//! structure in a dense fluid) naive staggering diverges, while Aitken-relaxed
//! coupling converges — demonstrated by [`iterate_fixed_relaxation`] vs
//! [`iterate_aitken`]. Deterministic; no dependencies.

/// The physical type of a power-conjugate port (its effort/flow pair).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortKind {
    /// Mechanical: force (effort) × velocity (flow).
    MechanicalForceVelocity,
    /// Fluid: pressure (effort) × volumetric flux (flow).
    FluidPressureFlux,
    /// Thermal: temperature (effort) × entropy flow.
    ThermalTemperatureEntropy,
}

/// A power port: an effort/flow pair. `power = effort × flow`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Port {
    /// The effort variable (force / pressure / temperature).
    pub effort: f64,
    /// The flow variable (velocity / flux / entropy flow).
    pub flow: f64,
    /// The physical type.
    pub kind: PortKind,
}

impl Port {
    /// A port.
    #[must_use]
    pub fn new(effort: f64, flow: f64, kind: PortKind) -> Port {
        Port { effort, flow, kind }
    }

    /// The instantaneous power `effort × flow`.
    #[must_use]
    pub fn power(&self) -> f64 {
        self.effort * self.flow
    }

    /// Are two ports POWER-CONJUGATE (composable) — the same physical
    /// effort/flow type? (The composition-time type discipline.)
    #[must_use]
    pub fn conjugate_to(&self, other: &Port) -> bool {
        self.kind == other.kind
    }
}

/// A structured coupling failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoupleError {
    /// The ports are not power-conjugate (mismatched physical types).
    IncompatiblePorts {
        /// The first port's kind.
        a: PortKind,
        /// The second port's kind.
        b: PortKind,
    },
}

/// A Dirac interconnection of two ports: shared effort, opposite flow — so the
/// interface power `e·f + e·(−f) = 0` EXACTLY (power-conserving by construction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interconnection {
    /// The first (side A) port.
    pub port_a: Port,
    /// The second (side B) port.
    pub port_b: Port,
    /// The net interface power (`0` by construction).
    pub interface_power: f64,
}

/// Interconnect two subsystems at a shared effort and flow through a Dirac
/// structure (effort continuity + flow balance).
///
/// # Errors
/// [`CoupleError::IncompatiblePorts`] if the ports are not power-conjugate.
pub fn interconnect(
    kind_a: PortKind,
    kind_b: PortKind,
    effort: f64,
    flow: f64,
) -> Result<Interconnection, CoupleError> {
    if kind_a != kind_b {
        return Err(CoupleError::IncompatiblePorts {
            a: kind_a,
            b: kind_b,
        });
    }
    let port_a = Port::new(effort, flow, kind_a);
    let port_b = Port::new(effort, -flow, kind_b);
    Ok(Interconnection {
        interface_power: port_a.power() + port_b.power(),
        port_a,
        port_b,
    })
}

/// The net interface power of a set of ports (`Σ effort·flow`) — `0` for a
/// power-conserving interconnection.
#[must_use]
pub fn interface_power(ports: &[Port]) -> f64 {
    ports.iter().map(Port::power).sum()
}

/// A caller-fed interface-balance audit.
///
/// The legacy [`EnergyAudit::is_passive`] name checks only whether every
/// recorded scalar interface imbalance stays within tolerance. It is not a
/// whole-system passivity certificate.
#[derive(Debug, Clone, Default)]
pub struct EnergyAudit {
    balances: Vec<f64>,
}

impl EnergyAudit {
    /// A fresh audit.
    #[must_use]
    pub fn new() -> EnergyAudit {
        EnergyAudit {
            balances: Vec::new(),
        }
    }

    /// Record one exchange's net interface power.
    pub fn record(&mut self, net_interface_power: f64) {
        self.balances.push(net_interface_power);
    }

    /// The worst interface power generation seen (the bug-alarm metric).
    ///
    /// A recorded NaN interface power means the coupling numerically broke
    /// down — the single worst thing this audit exists to catch. `f64::max`
    /// SILENTLY DROPS NaN (`f64::max(0.0, NaN) == 0.0`), so a plain fold would
    /// report zero imbalance and let the legacy `is_passive` predicate return
    /// true for a blown-up coupling. Poison instead: any NaN balance makes the
    /// metric NaN, and `NaN <= tol` is false, so the audit fails closed.
    /// (`±∞` already survives `f64::max` and alarms correctly.)
    #[must_use]
    pub fn max_generation(&self) -> f64 {
        if self.balances.iter().any(|b| b.is_nan()) {
            return f64::NAN;
        }
        self.balances.iter().map(|b| b.abs()).fold(0.0, f64::max)
    }

    /// Is every recorded interface-power imbalance within `tol`?
    ///
    /// This legacy name does not establish component or closed-window
    /// passivity; callers must audit those obligations separately.
    #[must_use]
    pub fn is_passive(&self, tol: f64) -> bool {
        self.max_generation() <= tol
    }
}

/// Scalar Aitken (Δ²) dynamic relaxation for the strongly-coupled interface
/// fixed point.
#[derive(Debug, Clone)]
pub struct AitkenRelaxation {
    omega: f64,
    omega_max: f64,
    prev_residual: Option<f64>,
}

impl AitkenRelaxation {
    /// A relaxer with an initial ω and a magnitude cap.
    #[must_use]
    pub fn new(omega_init: f64, omega_max: f64) -> AitkenRelaxation {
        AitkenRelaxation {
            omega: omega_init,
            omega_max,
            prev_residual: None,
        }
    }

    /// The Aitken relaxation factor for the current residual:
    /// `ωₖ = −ωₖ₋₁ · rₖ₋₁ / (rₖ − rₖ₋₁)` (scalar), magnitude-capped.
    pub fn next_omega(&mut self, residual: f64) -> f64 {
        if let Some(prev) = self.prev_residual {
            let dr = residual - prev;
            if dr.abs() > 1e-14 {
                let w = -self.omega * prev / dr;
                self.omega = w.clamp(-self.omega_max, self.omega_max);
            }
        }
        self.prev_residual = Some(residual);
        self.omega
    }
}

/// The result of an FSI interface fixed-point iteration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FsiResult {
    /// Did it converge (residual below `tol`) without blowing up?
    pub converged: bool,
    /// Iterations taken (or `max_steps` if it did not converge).
    pub steps: usize,
    /// The final interface value.
    pub solution: f64,
    /// The final residual magnitude.
    pub final_residual: f64,
}

// The classic linearized added-mass interface map: H(x) = −μ·x + c, where μ is
// the added-mass ratio (fluid added mass / structural mass). Naive staggering
// (ω = 1) converges only for μ < 1; a dense fluid on a light structure (μ ≥ 1)
// diverges.
fn interface_map(mu: f64, c: f64, x: f64) -> f64 {
    -mu * x + c
}

const BLOWUP: f64 = 1e12;

/// Iterate the added-mass interface fixed point with FIXED under-relaxation
/// `omega`. Diverges (naive staggering, `omega = 1`) when the added-mass ratio
/// `mu >= 1`.
#[must_use]
pub fn iterate_fixed_relaxation(
    mu: f64,
    c: f64,
    x0: f64,
    omega: f64,
    max_steps: usize,
    tol: f64,
) -> FsiResult {
    let mut x = x0;
    for step in 1..=max_steps {
        let r = interface_map(mu, c, x) - x;
        x += omega * r;
        if !x.is_finite() || x.abs() > BLOWUP {
            return FsiResult {
                converged: false,
                steps: step,
                solution: x,
                final_residual: f64::INFINITY,
            };
        }
        if r.abs() < tol {
            return FsiResult {
                converged: true,
                steps: step,
                solution: x,
                final_residual: r.abs(),
            };
        }
    }
    let r = (interface_map(mu, c, x) - x).abs();
    FsiResult {
        converged: r < tol,
        steps: max_steps,
        solution: x,
        final_residual: r,
    }
}

/// Iterate the same interface fixed point with AITKEN dynamic relaxation, which
/// stabilizes and accelerates it even for `mu >= 1` (the added-mass fix).
#[must_use]
pub fn iterate_aitken(
    mu: f64,
    c: f64,
    x0: f64,
    omega_init: f64,
    omega_max: f64,
    max_steps: usize,
    tol: f64,
) -> FsiResult {
    let mut x = x0;
    let mut aitken = AitkenRelaxation::new(omega_init, omega_max);
    for step in 1..=max_steps {
        let r = interface_map(mu, c, x) - x;
        if r.abs() < tol {
            return FsiResult {
                converged: true,
                steps: step,
                solution: x,
                final_residual: r.abs(),
            };
        }
        let omega = aitken.next_omega(r);
        x += omega * r;
        if !x.is_finite() || x.abs() > BLOWUP {
            return FsiResult {
                converged: false,
                steps: step,
                solution: x,
                final_residual: f64::INFINITY,
            };
        }
    }
    let r = (interface_map(mu, c, x) - x).abs();
    FsiResult {
        converged: r < tol,
        steps: max_steps,
        solution: x,
        final_residual: r,
    }
}
