//! Isothermal hydrostatic initialization of the original compressible solver.
//!
//! At rest the Guo populations are f_q = w_q rho (1 - c_q.g/(2 theta0)).
//! Collision changes the minus to a plus. For axis-aligned g, streaming is
//! therefore stationary when rho[j+1]/rho[j] = (1+b[j])/(1-b[j+1]),
//! b[j]=g[j]/(2 theta0), including a spatially varying acceleration.
//! Halfway stationary bounce-back satisfies the same balance at every wall.
//! This is a fixed-point initialization, not a modified time-step operator.

use super::{NsfError, NsfGasMedium, NsfGrid3, NsfSpec, THETA0, WallThermal};
use crate::lattice3::Lattice3;

fn admit_column<G, S>(
    spec: &NsfSpec<G, S>,
    axis: usize,
    reference_cell: usize,
) -> Result<(), NsfError> {
    if axis >= 3 || reference_cell >= spec.dims[axis] {
        return Err(NsfError::Input(
            "hydrostatic axis or reference cell is out of range",
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
    // Admit the original cell budget before allocating a layer profile.
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
    Ok(())
}

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
        admit_column(&spec, axis, reference_cell)?;
        if !acceleration_m_s2.is_finite() {
            return Err(NsfError::Input("hydrostatic gravity must be finite"));
        }
        if acceleration_m_s2 != 0.0 && spec.periodic[axis] {
            return Err(NsfError::Input(
                "a nonzero hydrostatic gravity axis cannot be periodic",
            ));
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

    /// Build a resting isothermal gas in a varying axial acceleration field.
    ///
    /// Supply one signed acceleration [m/s²] per cell layer along `axis`.
    /// The reference pressure and temperature belong to the CENTER of
    /// `reference_cell`; its lattice density is exactly one. The initialized
    /// pressure balances the actual density-dependent body force on each link:
    /// `p[j+1] - p[j] = dx (rho_SI[j]*a[j] + rho_SI[j+1]*a[j+1]) / 2`.
    /// This admits changing gravity and force reversals without launching
    /// pressure waves from a sampled continuum profile. Constant acceleration
    /// reduces to the same profile as [`Self::new_isothermal_hydrostatic`].
    ///
    /// The field is attached through the original spatial body-force owner,
    /// including its physical half-source populations and energy supply in
    /// `EnergyLedger::body_force_work_j`. Uniform gravity remains zero; callers
    /// must not attach the same acceleration again as a gravity control.
    /// The ordinary force, collision and streaming step preserves this rest
    /// state to roundoff, including stationary walls and embedded solids.
    /// It uses no extra evolution operator or per-step profile evaluation.
    ///
    /// A periodic axis is admitted only when its discrete integrated force
    /// closes the pressure profile around the complete loop. A nonzero uniform
    /// acceleration cannot support periodic hydrostatic pressure. Transverse
    /// axes may be periodic or closed. Walls must be adiabatic or held at the
    /// reference temperature. Nonuniform temperature, oblique forces and
    /// general well-balanced evolution are outside this initialization.
    ///
    /// Construction temporarily stores the one-dimensional density and force
    /// profile. The existing body-force owner retains its bounded native-cell
    /// field; an all-zero profile uses the allocation-free disabled force path.
    ///
    /// # Errors
    /// Refuses invalid shape/datum, incompatible thermal walls or periodic
    /// pressure closure, unrepresentable accelerations, nonpositive/overflowing
    /// density, or any original material, cell-budget or grid admission error.
    /// Every layer requires `|a[j]| dx / (2 R T_ref) < 1` for positive populations.
    pub fn new_isothermal_hydrostatic_profile<G: Into<NsfGasMedium>, S: Into<String>>(
        spec: NsfSpec<G, S>,
        solid: impl Fn([usize; 3]) -> bool,
        axis: usize,
        acceleration_m_s2: &[f64],
        reference_cell: usize,
    ) -> Result<Self, NsfError> {
        admit_column(&spec, axis, reference_cell)?;
        if acceleration_m_s2.len() != spec.dims[axis] {
            return Err(NsfError::Input(
                "hydrostatic acceleration needs one value per layer",
            ));
        }
        if acceleration_m_s2.iter().all(|&a| a == 0.0) {
            return Self::new_isothermal_hydrostatic(spec, solid, axis, 0.0, reference_cell);
        }
        let spec = spec.map_gas(Into::into);
        let state = spec.gas.reference(&spec)?;
        let dt_s = spec.dx_m / (3.0 * state.specific_gas_constant * state.temperature).sqrt();
        let mut b = Vec::with_capacity(acceleration_m_s2.len());
        for &a in acceleration_m_s2 {
            let value = (a * dt_s * (dt_s / spec.dx_m)) / (2.0 * THETA0);
            if !(value.is_finite() && value.abs() < 1.0 && (a == 0.0 || value != 0.0)) {
                return Err(NsfError::Input(
                    "hydrostatic layer acceleration requires |a_lattice| < 2 theta0",
                ));
            }
            b.push(value);
        }
        if spec.periodic[axis] {
            let mut log_cycle = 0.0;
            let mut correction = 0.0;
            let mut log_magnitude = 0.0;
            for &value in &b {
                // The denominators telescope around a periodic loop. log1p
                // retains tiny accelerations without rounding 1+b to one;
                // compensated summation admits a resolved conservative force
                // independently of the arbitrary reference-cell location.
                let term = value.ln_1p() - (-value).ln_1p();
                let corrected = term - correction;
                let sum = log_cycle + corrected;
                correction = (sum - log_cycle) - corrected;
                log_cycle = sum;
                log_magnitude += term.abs();
            }
            if log_cycle.abs() > 64.0 * f64::EPSILON * log_magnitude {
                return Err(NsfError::Input(
                    "hydrostatic acceleration does not close periodic pressure",
                ));
            }
        }
        let mut density = vec![1.0; b.len()];
        for j in reference_cell + 1..density.len() {
            density[j] = density[j - 1] * ((1.0 + b[j - 1]) / (1.0 - b[j]));
        }
        for j in (0..reference_cell).rev() {
            density[j] = density[j + 1] / ((1.0 + b[j]) / (1.0 - b[j + 1]));
        }
        if density.iter().any(|r| !(r.is_finite() && *r > 0.0)) {
            return Err(NsfError::Input(
                "hydrostatic density profile is not representable",
            ));
        }
        let temperature_k = spec.temperature_k;
        let grid = Self::new(spec, solid, |x| (density[x[axis]], [0.0; 3], temperature_k))?;
        let mut field = Vec::new();
        field
            .try_reserve_exact(grid.n)
            .map_err(|_| NsfError::Input("hydrostatic acceleration field allocation failed"))?;
        for i in 0..grid.n {
            let mut acceleration = [0.0; 3];
            acceleration[axis] = acceleration_m_s2[grid.at(i)[axis]];
            field.push(acceleration);
        }
        grid.with_body_acceleration_si(&field, || false)
    }
}
