//! Distributed cavity storage around an existing mass-normalized mechanical bank.
//!
//! Write y_j = sqrt(A_j) C_j.q + omega_j z_j, with A_j = rho*c²/Lambda_j.
//! Then H_air = (|y|²+|p|²)/2, y' = sqrt(A) C v + omega p,
//! p' = -omega y - d p, and the mechanical reaction is -C^T sqrt(A) y.
//! The uniform mode retains y but has NO fictitious momentum coordinate.
//!
//! The shared collective power exchange acts on mechanical v and pressure y;
//! the existing exact oscillator owner supplies free-air transitions. Calling
//! before / mechanical step / after is a symmetric second-order composition,
//! not an exact full-system propagator. Momentum drag remains the caller's
//! explicit causal loss, including its corresponding frequency-domain numerator.
use super::{CavityCoupling, ImpactError, invalid};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticState,
    ModalAcousticTimeBudget, ModalAcousticTimeModel};
use fs_math::c64::C64;
use fs_phs::{PortExchangeBudget, PreparedPortExchange, PreparedStepError};

const MAX_STRUCTURAL_MODES: usize = 128;

#[derive(Clone, Copy, Debug, Default)]
struct Free {
    cavity: usize,
    yy: f64,
    yp: f64,
    py: f64,
    pp: f64,
}

/// Prepared cavity memory in the supplied mass-normalized structural basis.
///
/// Construction starts at zero compression and zero acoustic momentum, so it
/// must precede structural excitation. Frame checkpoint/restore owns only the
/// cavity history; the enclosing mechanical transaction restores its own state.
/// Each before/after call is itself transactional and allocates no storage.
#[derive(Debug)]
pub struct PreparedCavityExchange {
    root_a: Vec<f64>,
    y: Vec<f64>,
    p: Vec<f64>,
    saved_y: Vec<f64>,
    saved_p: Vec<f64>,
    trial_y: Vec<f64>,
    trial_p: Vec<f64>,
    trial_v: Vec<f64>,
    free: Vec<Free>,
    exchange: PreparedPortExchange,
}

impl CavityCoupling {
    fn validate_exchange_basis(&self) -> Result<(), ImpactError> {
        if self.structural > MAX_STRUCTURAL_MODES || !self.necks.is_empty() {
            return Err(invalid("prepared cavity exchange admits at most 128 structural modes and no necks"));
        }
        Ok(())
    }

    // The same signed energy normalization serves time and harmonic responses.
    fn exchange_columns(&self) -> Result<(Vec<f64>, Vec<f64>), ImpactError> {
        self.validate_exchange_basis()?;
        let mut root_a = Vec::with_capacity(self.cavity_modes());
        let mut columns = Vec::with_capacity(self.structural * self.cavity_modes());
        for spring in &self.springs {
            let root = (spring.bulk_modulus_pa / spring.volume_m3).sqrt();
            root_a.push(root);
            for &c in &spring.areas[..self.structural] { columns.push(root*c); }
        }
        if columns.iter().any(|x| !x.is_finite()) {
            return Err(invalid("cavity pressure/velocity coupling overflow"));
        }
        Ok((root_a, columns))
    }

    /// Prepare two reusable half-step flows at the mechanical sample rate.
    ///
    /// The caller declares collective exchange work/rate limits. The existing
    /// 0.9 Nyquist guard applies at the FULL mechanical rate; preparation at
    /// twice that rate only obtains exact free-air half-step coefficients.
    /// All supplied modes are retained, including uniform compression. Neck
    /// states require the existing joint cavity solver and refuse here.
    pub fn prepare_exchange(&self, sample_rate_hz: u32, budget: PortExchangeBudget)
        -> Result<PreparedCavityExchange, ImpactError>
    {
        let half_rate = sample_rate_hz.checked_mul(2)
            .filter(|&r| r > 0).ok_or_else(|| invalid("invalid cavity half-step rate"))?;
        if self.omegas.iter().any(|&w| w >= 0.9*core::f64::consts::PI*f64::from(sample_rate_hz)) {
            return Err(invalid("cavity standing wave exceeds the mechanical Nyquist guard"));
        }
        let (root_a, columns) = self.exchange_columns()?;
        let exchange = PreparedPortExchange::new(self.structural, self.cavity_modes(),
            &columns, 0.5/f64::from(sample_rate_hz), budget)
            .map_err(ImpactError::PreparedSolve)?;
        let dynamic: Vec<_> = self.omegas.iter().enumerate()
            .filter_map(|(j,&w)| (w > 0.0).then_some(j)).collect();
        let mut free: Vec<_> = dynamic.iter().map(|&j| Free { cavity:j, ..Free::default() }).collect();
        if !dynamic.is_empty() {
            let modes = dynamic.iter().map(|&j| ModalAcousticMode {
                angular_frequency_rad_s: self.omegas[j],
                damping_ratio: 0.5*self.damping[j]/self.omegas[j],
                pressure_per_modal_velocity: C64::ZERO,
            }).collect();
            let mut oscillator = ModalAcousticTimeModel::try_new(half_rate, modes,
                ModalAcousticTimeBudget::audible_reference())
                .map_err(|e| ImpactError::Owner(e.to_string()))?;
            // Bounded probing also admits sub-radian frequencies without an
            // arbitrary 1/omega displacement exceeding the preparation budget.
            let mut states: Vec<_> = dynamic.iter().map(|&j| ModalAcousticState {
                displacement_m_sqrt_kg: self.omegas[j].max(1.0).recip(),
                velocity_m_sqrt_kg_per_s: 0.0,
            }).collect();
            let zeros = vec![0.0; dynamic.len()];
            oscillator.restore_states(&states).map_err(|e| ImpactError::Owner(e.to_string()))?;
            oscillator.step(&zeros).map_err(|e| ImpactError::Owner(e.to_string()))?;
            for ((t,s), initial) in free.iter_mut().zip(oscillator.states()).zip(&states) {
                let w = self.omegas[t.cavity];
                let initial_y = w*initial.displacement_m_sqrt_kg;
                t.yy = w*s.displacement_m_sqrt_kg/initial_y;
                t.py = s.velocity_m_sqrt_kg_per_s/initial_y;
            }
            states.fill(ModalAcousticState {
                displacement_m_sqrt_kg: 0.0, velocity_m_sqrt_kg_per_s: 1.0,
            });
            oscillator.restore_states(&states).map_err(|e| ImpactError::Owner(e.to_string()))?;
            oscillator.step(&zeros).map_err(|e| ImpactError::Owner(e.to_string()))?;
            for (t,s) in free.iter_mut().zip(oscillator.states()) {
                t.yp = self.omegas[t.cavity]*s.displacement_m_sqrt_kg;
                t.pp = s.velocity_m_sqrt_kg_per_s;
            }
            if free.iter().any(|t| [t.yy,t.yp,t.py,t.pp].iter().any(|x| !x.is_finite())) {
                return Err(invalid("cavity free-air transition is unrepresentable"));
            }
        }
        let count = self.cavity_modes(); let positive = free.len();
        Ok(PreparedCavityExchange {
            root_a, y:vec![0.0;count], p:vec![0.0;positive],
            saved_y:vec![0.0;count], saved_p:vec![0.0;positive],
            trial_y:vec![0.0;count], trial_p:vec![0.0;positive],
            trial_v:vec![0.0;self.structural], free, exchange,
        })
    }

    /// Reciprocal impedance opposing structural velocity in exp(-i omega t).
    ///
    /// Each pressure column contributes A C C^T (s+d)/(s²+d*s+omega_j²).
    /// In particular the uniform mode contributes the sealed compliance A C C^T/s;
    /// it is neither omitted nor replaced with a low-frequency oscillator.
    /// This is the continuous constitutive response, not a discrete-step fit.
    /// Exact undamped poles and nonrepresentable responses refuse.
    pub fn impedance(&self, omega_rad_s: f64) -> Result<Vec<C64>, ImpactError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("cavity impedance requires positive finite frequency"));
        }
        let (_, columns) = self.exchange_columns()?;
        let n = self.structural; let mut z = vec![C64::ZERO;n*n];
        for (j,column) in columns.chunks_exact(n).enumerate() {
            let den = C64::new(self.omegas[j]*self.omegas[j]-omega_rad_s*omega_rad_s,
                -self.damping[j]*omega_rad_s);
            if !den.re.is_finite() || !den.im.is_finite() || den.abs() == 0.0 {
                return Err(invalid("unresolved lossless cavity pole or impedance denominator"));
            }
            let h = C64::new(self.damping[j], -omega_rad_s)/den;
            for r in 0..n { for c in 0..n {
                z[r*n+c] = z[r*n+c] + h.scale(column[r]*column[c]);
            } }
        }
        if z.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(invalid("cavity impedance overflow"));
        }
        Ok(z)
    }

    /// Full cavity contribution to harmonic dynamic stiffness, without pole
    /// elimination. Row-major coordinates are the original structural prefix
    /// followed by this coupling's actual standing-wave inertias, exactly as in
    /// `total_modes()`. The uniform mode adds its spring without an inertia.
    ///
    /// The caller adds structural stiffness, damping, mass and external loads
    /// to the prefix. This contribution adds A*a*a^T for each actual volume
    /// spring and -omega²-i*omega*d ONLY on appended acoustic diagonals. It
    /// remains finite at an undamped fixed-wall cavity frequency, allowing a
    /// joint bordered solve instead of division by that acoustic pole.
    pub fn dynamic_stiffness(&self, omega_rad_s: f64) -> Result<Vec<C64>, ImpactError> {
        self.validate_exchange_basis()?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("cavity dynamic stiffness requires positive finite frequency"));
        }
        let n = self.total; let mut k = vec![C64::ZERO;n*n];
        for spring in &self.springs {
            let root = (spring.bulk_modulus_pa/spring.volume_m3).sqrt();
            for r in 0..n { for c in 0..n {
                k[r*n+c].re += (root*spring.areas[r])*(root*spring.areas[c]);
            } }
        }
        for (&coordinate,&drag) in self.dynamic.iter().zip(&self.damping) {
            if let Some(j) = coordinate {
                k[j*n+j] = k[j*n+j] + C64::new(-omega_rad_s*omega_rad_s,-omega_rad_s*drag);
            }
        }
        if k.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(invalid("cavity dynamic stiffness overflow"));
        }
        Ok(k)
    }
}

fn energy(y: &[f64], p: &[f64]) -> f64 {
    y.iter().chain(p).map(|x| 0.5*x*x).sum()
}

impl PreparedCavityExchange {
    /// Actual compression and acoustic kinetic energy [J], including uniform air.
    #[must_use]
    pub fn energy(&self) -> f64 { energy(&self.y, &self.p) }

    /// Save acoustic history before an enclosing complete-frame transaction.
    pub fn checkpoint(&mut self) {
        self.saved_y.copy_from_slice(&self.y); self.saved_p.copy_from_slice(&self.p);
    }
    /// Restore exactly the last checkpoint; the host restores its own mechanics.
    pub fn restore(&mut self) {
        self.y.copy_from_slice(&self.saved_y); self.p.copy_from_slice(&self.saved_p);
    }

    /// Pressure coefficients [Pa] in the originally supplied cavity basis.
    /// A dimension or finite-range failure leaves caller storage unchanged.
    pub fn pressures_into(&self, output: &mut [f64]) -> Result<(), ImpactError> {
        if output.len() != self.y.len() {
            return Err(invalid("cavity pressure output must cover the complete pressure basis"));
        }
        let mut pressure = [0.0;8];
        for (j,(&root,&y)) in self.root_a.iter().zip(&self.y).enumerate() {
            pressure[j] = -root*y;
            if !pressure[j].is_finite() { return Err(invalid("cavity pressure overflow")); }
        }
        output.copy_from_slice(&pressure[..self.y.len()]); Ok(())
    }

    fn advance(&mut self, board: &mut [f64], exchange_first: bool) -> Result<f64, ImpactError> {
        if board.len() != self.trial_v.len() || board.iter().any(|v| !v.is_finite()) {
            return Err(invalid("cavity exchange requires the complete finite structural velocity basis"));
        }
        self.trial_y.copy_from_slice(&self.y); self.trial_p.copy_from_slice(&self.p);
        self.trial_v.copy_from_slice(board);
        if exchange_first { self.exchange_half()?; }
        let mut loss = 0.0;
        for (p,t) in self.trial_p.iter_mut().zip(&self.free) {
            let y = self.trial_y[t.cavity]; let momentum = *p;
            let before = 0.5*y*y + 0.5*momentum*momentum;
            self.trial_y[t.cavity] = t.yy*y+t.yp*momentum;
            *p = t.py*y+t.pp*momentum;
            let after = 0.5*self.trial_y[t.cavity]*self.trial_y[t.cavity] + 0.5*(*p)*(*p);
            let tolerance = 128.0*f64::EPSILON*(before+after);
            if !after.is_finite() || !before.is_finite() || !tolerance.is_finite() {
                return Err(invalid("cavity free-air energy overflow"));
            }
            if after-before > tolerance {
                return Err(ImpactError::Energy {residual_j:after-before,tolerance_j:tolerance});
            }
            loss += before-after;
        }
        if !exchange_first { self.exchange_half()?; }
        if !(energy(&self.trial_y,&self.trial_p) + energy(&self.trial_v,&[])).is_finite() {
            return Err(invalid("cavity/structure combined energy overflow"));
        }
        self.y.copy_from_slice(&self.trial_y); self.p.copy_from_slice(&self.trial_p);
        board.copy_from_slice(&self.trial_v); Ok(loss)
    }

    fn exchange_half(&mut self) -> Result<(), ImpactError> {
        self.exchange.apply(&mut self.trial_v, &mut self.trial_y).map_err(|e| match e {
            PreparedStepError::Cancelled => ImpactError::Cancelled,
            PreparedStepError::Solver(e) => ImpactError::PreparedSolve(e),
        })?;
        // Its checked signed floating-point defect is NOT damping loss.
        Ok(())
    }

    /// Exchange half a step, then evolve free air half a step. Returns signed
    /// free-air energy loss [J], with unmodified floating-point roundoff.
    pub fn before(&mut self, board: &mut [f64]) -> Result<f64, ImpactError> {
        self.advance(board, true)
    }
    /// Evolve free air half a step, then exchange half a step. Call after the
    /// existing mechanical full step; this completes the symmetric composition.
    pub fn after(&mut self, board: &mut [f64]) -> Result<f64, ImpactError> {
        self.advance(board, false)
    }
}

