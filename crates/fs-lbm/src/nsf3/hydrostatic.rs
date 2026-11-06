//! Isothermal hydrostatic initialization of the original compressible solver.
//!
//! At rest the Guo populations are f_q = w_q rho (1 - c_q.g/(2 theta0)).
//! Collision changes the minus to a plus. For axis-aligned g, streaming is
//! therefore stationary when rho[j+1]/rho[j] = (1+b)/(1-b), b=g/(2 theta0).
//! Halfway stationary bounce-back satisfies the same balance at every wall.
//! This is a fixed-point initialization, not a modified time-step operator.

use super::{NsfError, NsfGasMedium, NsfGrid3, NsfSpec, THETA0, WallThermal};
use crate::lattice3::Lattice3;

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Build a resting isothermal gas column with native hydrostatic balance.
    ///
    /// `axis` selects x/y/z (0/1/2), and `acceleration_m_s2` is signed along
    /// its increasing cell index. The reference pressure in `spec` belongs to
    /// the CENTER of `reference_cell` on that axis; every fluid cell in that
    /// plane has lattice density 1. Temperature is `spec.temperature_k`.
    /// Gravity is attached to the returned grid, including its energy work.
    ///
    /// Adjacent layers obey
    /// `p[j+1] - p[j] = g dx (rho_SI[j+1] + rho_SI[j]) / 2`.
    /// This geometric profile is a fixed point of the existing D3Q27 force,
    /// collision and streaming operators, up to floating-point roundoff.
    /// It approaches the ideal-gas barometric exponential with second-order
    /// spatial error. It is deliberately not the sampled continuum
    /// exponential, which would launch a small grid-dependent transient.
    ///
    /// The gravity axis must be nonperiodic for nonzero gravity. Other axes
    /// may be periodic or stationary walls; embedded stationary solids are
    /// allowed. Walls must be adiabatic or held at the reference temperature.
    /// GasSpec, named-mixture and thermochemical inputs retain their original
    /// local material checks and grid limits. Positive populations require
    /// `|g| dx / (2 R T_ref) < 1`; no density or pressure is clipped.
    /// Opening reservoirs, moving walls, changing gravity or adding a thermal
    /// perturbation can subsequently disturb this particular resting state.
    /// No general well-balanced evolution for oblique gravity or nonuniform
    /// temperature is implied. Construction adds one density value per layer
    /// temporarily; it adds no persistent fields or work to `step`.
    ///
    /// # Errors
    /// Refuses an invalid axis/datum, incompatible boundary, unrepresentable
    /// gravity, nonpositive/overflowing profile, or ordinary grid admission.
    pub fn new_isothermal_hydrostatic<G: Into<NsfGasMedium>, S: Into<String>>(
        spec: NsfSpec<G, S>,
        solid: impl Fn([usize; 3]) -> bool,
        axis: usize,
        acceleration_m_s2: f64,
        reference_cell: usize,
    ) -> Result<Self, NsfError> {
        if axis >= 3 || reference_cell >= spec.dims[axis] {
            return Err(NsfError::Input(
                "hydrostatic axis or reference cell is out of range",
            ));
        }
        if !acceleration_m_s2.is_finite() {
            return Err(NsfError::Input("hydrostatic gravity must be finite"));
        }
        if acceleration_m_s2 != 0.0 && spec.periodic[axis] {
            return Err(NsfError::Input(
                "a nonzero hydrostatic gravity axis cannot be periodic",
            ));
        }
        if let WallThermal::Isothermal(t) = spec.wall {
            if t != spec.temperature_k {
                return Err(NsfError::Input(
                    "hydrostatic walls must share the reference temperature",
                ));
            }
        }
        if !(spec.dx_m.is_finite() && spec.dx_m > 0.0) {
            return Err(NsfError::Input(
                "hydrostatic cell size must be positive and finite",
            ));
        }
        // Admit the original cell budget before allocating the layer profile.
        let n = spec
            .dims
            .iter()
            .try_fold(1usize, |n, &d| n.checked_mul(d))
            .ok_or(NsfError::Input("dimension overflow"))?;
        if spec.dims.contains(&0) {
            return Err(NsfError::Input("dimensions must be positive"));
        }
        if n > spec.maximum_cells {
            return Err(NsfError::Cap {
                requested: n,
                maximum: spec.maximum_cells,
            });
        }
        let spec = spec.map_gas(Into::into);
        let temperature_k = spec.temperature_k;
        if acceleration_m_s2 == 0.0 {
            return Self::new(spec, solid, |_| (1.0, [0.0; 3], temperature_k));
        }
        // Use exactly the original constructor and gravity binding's SI scale.
        let state = spec.gas.reference(&spec)?;
        let dt_s = spec.dx_m / (3.0 * state.specific_gas_constant * state.temperature).sqrt();
        let acceleration_lattice = acceleration_m_s2 * dt_s * (dt_s / spec.dx_m);
        let b = acceleration_lattice / (2.0 * THETA0);
        if !(b.is_finite() && b != 0.0 && b.abs() < 1.0) {
            return Err(NsfError::Input(
                "hydrostatic gravity requires 0 < |g_lattice| < 2 theta0",
            ));
        }
        let ratio = (1.0 + b) / (1.0 - b);
        let mut density = vec![1.0; spec.dims[axis]];
        for j in reference_cell + 1..density.len() {
            density[j] = density[j - 1] * ratio;
        }
        for j in (0..reference_cell).rev() {
            density[j] = density[j + 1] / ratio;
        }
        if density.iter().any(|r| !(r.is_finite() && *r > 0.0)) {
            return Err(NsfError::Input(
                "hydrostatic density profile is not representable",
            ));
        }
        let grid = Self::new(spec, solid, |x| (density[x[axis]], [0.0; 3], temperature_k))?;
        let mut gravity = [0.0; 3];
        gravity[axis] = acceleration_m_s2;
        grid.with_gravity_si(gravity)
    }
}
