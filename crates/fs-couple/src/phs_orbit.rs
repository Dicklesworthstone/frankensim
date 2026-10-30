//! Finite-amplitude periodic regimes of complete nonlinear pHS ODEs.
//!
//! This is an adapter, not a second orbit solver. fs-orbit owns harmonic balance,
//! shooting and continuation; fs-phs owns energy, interconnection, dissipation
//! and material-memory coordinates. The adapter evaluates the ORIGINAL continuous
//! vector field at every AFT sample and exposes only the default sI frequency port.
//! Descriptor systems are intentionally excluded because their algebraic rows need
//! an explicit descriptor OrbitProblem rather than an inverted mass matrix.
//!
//! The path is useful after small-signal pole analysis: linear poles seed the
//! first periodic branch, while finite amplitude keeps nonlinear storage and every
//! explicit Maxwell coordinate. Harmonic truncation remains disclosed structure;
//! time-domain confirmation is still required for attacks and aperiodic regimes.

use core::cell::{Cell, RefCell};
use fs_math::{c64::C64, det};
use fs_orbit::{HbAnchor, HbBudget, HbOrbit, OrbitError, OrbitProblem, solve_hb};
use fs_phs::PortHamiltonian;

const MAX_INPUT_HARMONICS: usize = 256;

#[derive(Debug, Clone, PartialEq)]
pub enum PhsOrbitError {
    Invalid { what: &'static str },
    /// The original storage callback produced a non-finite effort at an AFT sample.
    SourceRefused,
    Orbit(OrbitError),
}
impl core::fmt::Display for PhsOrbitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "pHS periodic regime refused: {self:?}")
    }
}
impl std::error::Error for PhsOrbitError {}
impl From<OrbitError> for PhsOrbitError {
    fn from(value: OrbitError) -> Self { Self::Orbit(value) }
}
fn bad(what: &'static str) -> PhsOrbitError { PhsOrbitError::Invalid { what } }

/// Peak-periodic physical port input in fs-orbit's exp(+i*n*omega*t) chart.
/// coefficients[port][0] is the real DC value; n>=1 stores the positive-frequency
/// complex Fourier coefficient, so a real cosine of peak amplitude A uses A/2.
#[derive(Debug, Clone)]
pub struct PeriodicPortInput {
    omega_rad_s: f64,
    coefficients: Vec<Vec<C64>>,
}
impl PeriodicPortInput {
    pub fn new(omega_rad_s: f64, coefficients: Vec<Vec<C64>>) -> Result<Self, PhsOrbitError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 || coefficients.len() > 256 {
            return Err(bad("positive finite frequency and at most 256 ports required"));
        }
        let width = coefficients.first().map_or(1, Vec::len);
        if width == 0 || width > MAX_INPUT_HARMONICS + 1
            || coefficients.iter().any(|row| row.len() != width
                || row.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
                || row[0].im != 0.0)
        {
            return Err(bad("complete finite Fourier rows with real DC terms required"));
        }
        Ok(Self { omega_rad_s, coefficients })
    }
    pub fn omega_rad_s(&self) -> f64 { self.omega_rad_s }
    pub fn port_count(&self) -> usize { self.coefficients.len() }
    pub fn harmonic_count(&self) -> usize {
        self.coefficients.first().map_or(0, |row| row.len().saturating_sub(1))
    }
    fn sample_into(&self, time_s: f64, out: &mut [f64]) -> bool {
        if out.len() != self.coefficients.len() || !time_s.is_finite() { return false; }
        for (port, row) in self.coefficients.iter().enumerate() {
            let mut value = row[0].re;
            let mut correction = 0.0;
            for (n, &coefficient) in row.iter().enumerate().skip(1) {
                let phase = n as f64 * self.omega_rad_s * time_s;
                let term = 2.0 * (coefficient.re * det::cos(phase)
                    - coefficient.im * det::sin(phase));
                let adjusted = term - correction;
                let next = value + adjusted;
                correction = (next - value) - adjusted;
                value = next;
            }
            if !value.is_finite() { return false; }
            out[port] = value;
        }
        true
    }
}

/// Private adapter so source-domain refusals can be promoted by the public solve
/// wrappers instead of escaping as an undifferentiated Newton stall.
struct PhsOrbitProblem<'a> {
    system: &'a PortHamiltonian,
    input: InputKind<'a>,
    effort: RefCell<Vec<f64>>,
    ports: RefCell<Vec<f64>>,
    source_refused: Cell<bool>,
}
enum InputKind<'a> {
    Forced(&'a PeriodicPortInput),
    Constant(&'a [f64]),
}
impl PhsOrbitProblem<'_> {
    fn source_refused(&self) -> bool { self.source_refused.get() }
    fn fail(&self, out: &mut [f64]) {
        self.source_refused.set(true);
        out.fill(f64::NAN);
    }
}
impl OrbitProblem for PhsOrbitProblem<'_> {
    fn dim(&self) -> usize { self.system.state_dim() }
    fn autonomous(&self) -> bool { matches!(self.input, InputKind::Constant(_)) }

    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let n = self.system.state_dim();
        let m = self.system.port_dim();
        if x.len() != n || out.len() != n || x.iter().any(|v| !v.is_finite()) {
            self.fail(out); return;
        }
        let mut effort = self.effort.borrow_mut();
        effort.fill(f64::NAN);
        self.system.storage().gradient(x, &mut effort);
        if effort.iter().any(|v| !v.is_finite()) { self.fail(out); return; }

        let mut input = self.ports.borrow_mut();
        let ok = match self.input {
            InputKind::Forced(periodic) => periodic.sample_into(t, &mut input),
            InputKind::Constant(bias) => {
                if bias.len() != m || bias.iter().any(|v| !v.is_finite()) { false }
                else { input.copy_from_slice(bias); true }
            }
        };
        if !ok { self.fail(out); return; }

        let (j, r, g) = self.system.structure();
        for i in 0..n {
            let mut value = 0.0;
            let mut correction = 0.0;
            for k in 0..n {
                let term = (j[i*n+k] - r[i*n+k]) * effort[k];
                let adjusted = term - correction;
                let next = value + adjusted;
                correction = (next - value) - adjusted;
                value = next;
            }
            for p in 0..m {
                let term = g[i*m+p] * input[p];
                let adjusted = term - correction;
                let next = value + adjusted;
                correction = (next - value) - adjusted;
                value = next;
            }
            if !value.is_finite() { self.fail(out); return; }
            out[i] = value;
        }
    }
}

fn forced_problem<'a>(
    system: &'a PortHamiltonian, input: &'a PeriodicPortInput,
) -> Result<PhsOrbitProblem<'a>, PhsOrbitError> {
    if system.state_dim() == 0 || system.port_dim() != input.port_count() {
        return Err(bad("periodic input must cover every original physical port"));
    }
    Ok(PhsOrbitProblem {
        system, input: InputKind::Forced(input),
        effort: RefCell::new(vec![0.0; system.state_dim()]),
        ports: RefCell::new(vec![0.0; system.port_dim()]),
        source_refused: Cell::new(false),
    })
}
fn autonomous_problem<'a>(
    system: &'a PortHamiltonian, bias: &'a [f64],
) -> Result<PhsOrbitProblem<'a>, PhsOrbitError> {
    if system.state_dim() == 0 || bias.len() != system.port_dim()
        || bias.iter().any(|v| !v.is_finite())
    {
        return Err(bad("constant bias must cover every original physical port"));
    }
    Ok(PhsOrbitProblem {
        system, input: InputKind::Constant(bias),
        effort: RefCell::new(vec![0.0; system.state_dim()]),
        ports: RefCell::new(bias.to_vec()),
        source_refused: Cell::new(false),
    })
}

/// Forced finite-amplitude harmonic balance using the original nonlinear pHS.
/// The forcing frequency is taken from input, preventing an inconsistent anchor.
pub fn solve_forced_harmonic_balance(
    system: &PortHamiltonian, input: &PeriodicPortInput,
    guess_amplitude: f64, budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    if !guess_amplitude.is_finite() || guess_amplitude <= 0.0 {
        return Err(bad("positive finite first-state seed amplitude required"));
    }
    let problem = forced_problem(system, input)?;
    let result = solve_hb(&problem, HbAnchor::Forced { omega: input.omega_rad_s() },
        guess_amplitude, budget);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Autonomous finite-amplitude orbit about a constant physical port bias.
pub fn solve_autonomous_harmonic_balance(
    system: &PortHamiltonian, bias: &[f64], omega_guess: f64,
    guess_amplitude: f64, budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    if !omega_guess.is_finite() || omega_guess <= 0.0
        || !guess_amplitude.is_finite() || guess_amplitude <= 0.0
    {
        return Err(bad("positive finite frequency and seed amplitude required"));
    }
    let problem = autonomous_problem(system, bias)?;
    let result = solve_hb(&problem, HbAnchor::Autonomous { omega_guess },
        guess_amplitude, budget);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

/// Conservative/backbone point with first-state amplitude pinned by fs-orbit.
/// The source system is not required to be lossless; if dissipation prevents a
/// periodic backbone the original balance refuses rather than deleting the loss.
pub fn solve_backbone_harmonic_balance(
    system: &PortHamiltonian, bias: &[f64], amplitude: f64, omega_guess: f64,
    budget: &HbBudget,
) -> Result<HbOrbit, PhsOrbitError> {
    if !amplitude.is_finite() || amplitude <= 0.0
        || !omega_guess.is_finite() || omega_guess <= 0.0
    {
        return Err(bad("positive finite backbone amplitude and frequency required"));
    }
    let problem = autonomous_problem(system, bias)?;
    let result = solve_hb(&problem, HbAnchor::Backbone { amplitude, omega_guess },
        amplitude, budget);
    if problem.source_refused() { return Err(PhsOrbitError::SourceRefused); }
    result.map_err(PhsOrbitError::Orbit)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_phs::{PortHamiltonian, QuadraticStorage, RelaxationBranch, Storage};

    fn hb(n: usize) -> HbBudget {
        HbBudget { harmonics:n, max_newton:60, tolerance:1e-11, deflate_equilibrium:false }
    }

    #[test]
    fn forced_linear_oscillator_matches_closed_form_amplitude() {
        let (m,k,d,w,f) = (2.0_f64, 1000.0_f64, 4.0_f64, 15.0_f64, 3.0_f64);
        let system = fs_phs::mass_spring_damper(m,k,d).unwrap();
        let input = PeriodicPortInput::new(w, vec![vec![C64::ZERO,C64::from_re(0.5*f)]]).unwrap();
        let expected = f / det::sqrt((k-m*w*w).powi(2)+(d*w).powi(2));
        let orbit = solve_forced_harmonic_balance(&system,&input,expected,&hb(3)).unwrap();
        assert!((orbit.first_harmonic_amplitude(0)-expected).abs() < 2e-10*expected);
        assert!((orbit.first_harmonic_amplitude(1)-w*expected*m).abs() < 2e-10*w*expected*m);
        assert!(orbit.coeffs.iter().all(|row| row.iter().skip(2).all(|z|z.abs()<1e-10)));
    }

    #[test]
    fn maxwell_memory_is_retained_in_the_periodic_balance() {
        let (k,kb,tau,w,f) = (60_000.0_f64,50_000.0_f64,1.0/300.0,240.0_f64,2.0_f64);
        let system = fs_phs::mass_spring_damper(1.0,k,0.0).unwrap()
            .with_relaxation_branches(vec![RelaxationBranch {
                projection:vec![1.0,0.0], stiffness:kb, relaxation_time_s:tau,
            }]).unwrap();
        let input = PeriodicPortInput::new(w,vec![vec![C64::ZERO,C64::from_re(0.5*f)]]).unwrap();
        let iw = C64::new(0.0,w);
        let modulus = (iw.scale(tau)/(C64::ONE+iw.scale(tau))).scale(kb);
        let denominator = C64::from_re(k-w*w)+modulus;
        let expected = f/denominator.abs();
        let orbit = solve_forced_harmonic_balance(&system,&input,expected,&hb(3)).unwrap();
        assert_eq!(orbit.coeffs.len(),3);
        assert!((orbit.first_harmonic_amplitude(0)-expected).abs()<5e-9*expected);
        assert!(orbit.first_harmonic_amplitude(2)>0.0);
    }

    struct Duffing { k:f64, alpha:f64, inv_mass:f64 }
    impl Storage for Duffing {
        fn hamiltonian(&self,x:&[f64])->f64 {
            0.5*self.k*x[0]*x[0]+0.25*self.alpha*x[0].powi(4)+0.5*self.inv_mass*x[1]*x[1]
        }
        fn gradient(&self,x:&[f64],out:&mut[f64]) {
            out[0]=self.k*x[0]+self.alpha*x[0].powi(3);
            out[1]=self.inv_mass*x[1];
        }
    }
    #[test]
    fn nonlinear_backbone_keeps_the_original_cubic_storage() {
        let (k,alpha,m,a)=(400.0_f64,100.0_f64,1.0_f64,0.1_f64);
        let system=PortHamiltonian::new(2,0,vec![0.,1.,-1.,0.],vec![0.;4],vec![],
            Box::new(Duffing{k,alpha,inv_mass:1.0/m})).unwrap();
        let first_order=det::sqrt((k+0.75*alpha*a*a)/m);
        let orbit=solve_backbone_harmonic_balance(&system,&[],a,first_order,&hb(7)).unwrap();
        assert!((orbit.first_harmonic_amplitude(0)-a).abs()<1e-10);
        assert!((orbit.omega-first_order).abs()/first_order<2e-4);
        assert!(orbit.coeffs[0][3].abs()>0.0, "cubic storage must generate a third harmonic");
    }

    #[test]
    fn input_shape_and_source_failure_are_refusals_not_partial_orbits() {
        assert!(PeriodicPortInput::new(10.0,vec![vec![C64::new(0.0,1.0)]]).is_err());
        let system=fs_phs::mass_spring_damper(1.0,100.0,0.2).unwrap();
        let wrong=PeriodicPortInput::new(10.0,vec![]).unwrap();
        assert!(solve_forced_harmonic_balance(&system,&wrong,0.1,&hb(3)).is_err());
        let ok=PeriodicPortInput::new(10.0,vec![vec![C64::ZERO,C64::from_re(0.5)]]).unwrap();
        let a=solve_forced_harmonic_balance(&system,&ok,0.1,&hb(3)).unwrap();
        let b=solve_forced_harmonic_balance(&system,&ok,0.1,&hb(3)).unwrap();
        assert_eq!(a.omega.to_bits(),b.omega.to_bits());
        assert_eq!(a.coeffs,b.coeffs);
    }
}
