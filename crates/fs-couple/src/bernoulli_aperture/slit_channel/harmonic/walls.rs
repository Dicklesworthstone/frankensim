//! Resolved frequency-domain pressure and axial moving-wall ports.
//!
//! The original cell equations are (-i*omega*M + G)*v = A*dp + g0*U0 + gn*Un.
//! This module solves that tridiagonal system directly and observes all three
//! physical ports under exp(-i*omega*t). Wall forces act ON the fluid; the
//! mechanical owner receives their opposite. Complex wall power is
//! F_on_fluid*conj(U)/2, so pressure plus wall power is heat-i*2*omega*mean(E).
//! No time step, fitted impedance, frozen profile or numerical heat enters.
//!
//! Applies to the existing fixed-geometry Newtonian laminar reduction. The
//! underlying owner refuses a nonlinear outlet before exposing this operator.

use super::{AcousticRealizeError, ChannelContinuousOperator, finite, invalid};
use fs_math::c64::C64;

/// Conservative arithmetic-work envelope for a solve plus physical observation.
/// All work and shape checks precede allocation; each loop polls cancellation.
pub const WALL_HARMONIC_WORK_PER_CELL: usize = 512;

/// Peak phasors and separate mechanical, pressure, storage and heat accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct ChannelWallHarmonicResponse {
    pub velocities_m_s: Vec<C64>,
    pub pressure_drop_pa: C64,
    /// Lower/upper for a slit, inner/outer for an annulus, axis/outer for a circle.
    pub wall_velocities_m_s: [C64; 2],
    pub flow_m3_s: C64,
    /// Axial forces exerted by the walls ON the fluid [N].
    pub wall_forces_n: [C64; 2],
    /// Cycle-mean fluid kinetic inventory, sum(M_i*|v_i|^2)/4 [J].
    pub mean_stored_energy_j: f64,
    /// Cycle-mean irreversible heat from relative shear [W].
    pub viscous_dissipation_w: f64,
    /// dp*conj(Q)/2 [W]; the real part may be negative during work recovery.
    pub pressure_complex_power_w: C64,
    /// F_on_fluid*conj(U)/2 at each original wall [W], not its conjugate.
    pub wall_complex_power_w: [C64; 2],
    /// Pressure plus both wall complex powers [W].
    pub total_complex_power_w: C64,
    /// Relative-shear heat minus the real total input power [W].
    pub power_residual_w: f64,
    /// Imaginary total input power plus 2*omega*mean(E) [W].
    pub reactive_power_residual_w: f64,
    /// Largest cell-equation residual relative to its uncancelled terms.
    pub solve_relative_residual: f64,
}

impl ChannelContinuousOperator<'_> {
    /// Solve the original continuous momentum equations in O(n) work and space.
    ///
    /// Two bounded vectors hold the upper elimination ratios and the complete
    /// velocity field. There is no dense matrix, transverse-state elimination
    /// or second copy of the solved field in the response. This is a frequency
    /// query, not an allocation-free time-step operation.
    ///
    /// Frequency must be strictly positive because these are peak phasors and
    /// cycle means. The first wall speed must be zero at a circular symmetry
    /// axis. The work budget charges WALL_HARMONIC_WORK_PER_CELL*n before
    /// allocating. No accepted channel history or source is modified.
    pub fn solve_harmonic_with_walls(
        self,
        omega_rad_s: f64,
        pressure_drop_pa: C64,
        walls_m_s: [C64; 2],
        maximum_work: usize,
        maximum_relative_residual: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelWallHarmonicResponse, AcousticRealizeError> {
        self.admit_wall_harmonic(
            omega_rad_s,
            pressure_drop_pa,
            walls_m_s,
            maximum_work,
            maximum_relative_residual,
            cancelled,
        )?;
        let n = self.dimension();
        let mut upper = vec![C64::ZERO; n];
        let mut velocities = vec![C64::ZERO; n];
        for i in 0..n {
            poll(cancelled)?;
            let down = self.conductance[i];
            let up = self.conductance[i + 1];
            let mut pivot = C64::new(down + up, -omega_rad_s * self.mass[i]);
            let mut rhs = pressure_drop_pa.scale(self.area[i]);
            if i == 0 {
                rhs = rhs + walls_m_s[0].scale(down);
            } else {
                pivot = pivot - upper[i - 1].scale(down);
                rhs = rhs + velocities[i - 1].scale(down);
            }
            if i + 1 == n {
                rhs = rhs + walls_m_s[1].scale(up);
            }
            if !finite(pivot) || pivot.re <= 0.0 || pivot.im > 0.0 || !finite(rhs) {
                return Err(invalid(
                    "harmonic channel factorization is not finite passive",
                ));
            }
            velocities[i] = rhs / pivot;
            if i + 1 < n {
                upper[i] = C64::from_re(up) / pivot;
            }
            if !finite(velocities[i]) || !finite(upper[i]) {
                return Err(invalid(
                    "harmonic channel elimination cannot be represented",
                ));
            }
        }
        for i in (0..n - 1).rev() {
            poll(cancelled)?;
            velocities[i] = velocities[i] + upper[i] * velocities[i + 1];
        }
        let mut response = self.inspect_wall_harmonic(
            omega_rad_s,
            &velocities,
            pressure_drop_pa,
            walls_m_s,
            maximum_relative_residual,
            cancelled,
        )?;
        response.velocities_m_s = velocities;
        Ok(response)
    }

    /// Observe a coupled solver's complete field against the same original
    /// equations and physical wall ports. The returned field allocates only
    /// after all equation, work and cancellation checks succeed.
    ///
    /// Relative-shear heat is independently accumulated. The raw real/reactive
    /// power residuals remain visible; their roundoff allowance includes the
    /// uncancelled cell/wall work at nearly rigid, low-frequency translation.
    pub fn observe_harmonic_with_walls(
        self,
        omega_rad_s: f64,
        velocities_m_s: &[C64],
        pressure_drop_pa: C64,
        walls_m_s: [C64; 2],
        maximum_work: usize,
        maximum_relative_residual: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelWallHarmonicResponse, AcousticRealizeError> {
        self.admit_wall_harmonic(
            omega_rad_s,
            pressure_drop_pa,
            walls_m_s,
            maximum_work,
            maximum_relative_residual,
            cancelled,
        )?;
        let mut response = self.inspect_wall_harmonic(
            omega_rad_s,
            velocities_m_s,
            pressure_drop_pa,
            walls_m_s,
            maximum_relative_residual,
            cancelled,
        )?;
        response.velocities_m_s = velocities_m_s.to_vec();
        Ok(response)
    }

    fn admit_wall_harmonic(
        self,
        omega: f64,
        pressure: C64,
        walls: [C64; 2],
        maximum_work: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AcousticRealizeError> {
        poll(cancelled)?;
        let n = self.dimension();
        if n == 0
            || self.area.len() != n
            || self.conductance.len() != n + 1
            || !omega.is_finite()
            || omega <= 0.0
            || !finite(pressure)
            || walls.iter().any(|&u| !finite(u))
            || !tolerance.is_finite()
            || tolerance <= 0.0
            || tolerance >= 1.0
            || n.checked_mul(WALL_HARMONIC_WORK_PER_CELL)
                .is_none_or(|work| work > maximum_work)
        {
            return Err(invalid(
                "harmonic wall ports need positive frequency, finite phasors, tolerance and bounded work",
            ));
        }
        if self.conductance[0] == 0.0 && walls[0] != C64::ZERO {
            return Err(invalid(
                "a circular symmetry axis is not a moving material wall",
            ));
        }
        Ok(())
    }

    // Return the fully admitted account without copying the field. The solver
    // moves its own velocity vector into it; an external observation copies its
    // read-only input only after this same acceptance succeeds.
    fn inspect_wall_harmonic(
        self,
        omega: f64,
        velocities: &[C64],
        pressure: C64,
        walls: [C64; 2],
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ChannelWallHarmonicResponse, AcousticRealizeError> {
        let n = self.dimension();
        if velocities.len() != n {
            return Err(invalid(
                "harmonic wall observation requires every original velocity cell",
            ));
        }
        let roundoff = 64.0 * n as f64 * f64::EPSILON;
        let equation_tolerance = relative_tolerance + roundoff;
        let mut flow = C64::ZERO;
        let mut stored = 0.0;
        let mut heat = 0.0;
        let mut work_roundoff_scale = 0.0;
        let mut largest_residual = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let v = velocities[i];
            let left = if i == 0 { walls[0] } else { velocities[i - 1] };
            let right = if i + 1 == n {
                walls[1]
            } else {
                velocities[i + 1]
            };
            if !finite(v) || !finite(left) || !finite(right) {
                return Err(invalid(
                    "harmonic wall observation needs finite velocity phasors",
                ));
            }
            let inertia = C64::new(0.0, -omega * self.mass[i]) * v;
            let down = (v - left).scale(self.conductance[i]);
            let up = (v - right).scale(self.conductance[i + 1]);
            let forcing = pressure.scale(self.area[i]);
            let residual = inertia + down + up - forcing;
            let scale = inertia.abs()
                + forcing.abs()
                + self.conductance[i] * (v.abs() + left.abs())
                + self.conductance[i + 1] * (v.abs() + right.abs());
            if !finite(residual)
                || !scale.is_finite()
                || residual.abs() > equation_tolerance * scale.max(1e-300)
            {
                return Err(invalid(
                    "harmonic wall field violates original cell momentum",
                ));
            }
            largest_residual = largest_residual.max(residual.abs() / scale.max(1e-300));
            flow = flow + v.scale(self.area[i]);
            stored += 0.25 * self.mass[i] * v.norm_sq();
            heat += 0.5 * self.conductance[i] * (v - left).norm_sq();
            work_roundoff_scale += 0.5 * v.abs() * scale;
        }
        heat += 0.5 * self.conductance[n] * (velocities[n - 1] - walls[1]).norm_sq();
        let forces = [
            (walls[0] - velocities[0]).scale(self.conductance[0]),
            (walls[1] - velocities[n - 1]).scale(self.conductance[n]),
        ];
        let pressure_power = (pressure * flow.conj()).scale(0.5);
        let wall_power = [
            (forces[0] * walls[0].conj()).scale(0.5),
            (forces[1] * walls[1].conj()).scale(0.5),
        ];
        let total_power = pressure_power + wall_power[0] + wall_power[1];
        let real_residual = heat - total_power.re;
        let reactive = 2.0 * omega * stored;
        let reactive_residual = total_power.im + reactive;
        let power_scale =
            heat + pressure_power.abs() + wall_power[0].abs() + wall_power[1].abs() + reactive;
        work_roundoff_scale +=
            0.5 * self.conductance[0] * (walls[0].abs() + velocities[0].abs()) * walls[0].abs()
                + 0.5
                    * self.conductance[n]
                    * (walls[1].abs() + velocities[n - 1].abs())
                    * walls[1].abs();
        // Physical acceptance plus an absolute floating-point allowance from
        // the actual uncancelled work. Finite inventory must not be rejected
        // when nearly common wall/fluid translation makes exchanges vanish.
        let allowance = relative_tolerance * power_scale + roundoff * work_roundoff_scale;
        if !finite(flow)
            || forces.iter().any(|&f| !finite(f))
            || !finite(pressure_power)
            || wall_power.iter().any(|&p| !finite(p))
            || !finite(total_power)
            || ![
                stored,
                heat,
                real_residual,
                reactive,
                reactive_residual,
                power_scale,
                work_roundoff_scale,
                allowance,
            ]
            .iter()
            .all(|x| x.is_finite())
            || real_residual.abs() > allowance
            || reactive_residual.abs() > allowance
        {
            return Err(invalid(
                "harmonic pressure and wall work disagree with fluid storage and relative-shear heat",
            ));
        }
        poll(cancelled)?;
        Ok(ChannelWallHarmonicResponse {
            velocities_m_s: Vec::new(),
            pressure_drop_pa: pressure,
            wall_velocities_m_s: walls,
            flow_m3_s: flow,
            wall_forces_n: forces,
            mean_stored_energy_j: stored,
            viscous_dissipation_w: heat,
            pressure_complex_power_w: pressure_power,
            wall_complex_power_w: wall_power,
            total_complex_power_w: total_power,
            power_residual_w: real_residual,
            reactive_power_residual_w: reactive_residual,
            solve_relative_residual: largest_residual,
        })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticRealizeError> {
    if cancelled() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}
