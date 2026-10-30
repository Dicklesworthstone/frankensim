//! Physical port participation of complete continuous pHS descriptor modes.
//! The existing fs-la prepared-pole owner does every factorization and solve.
//! Constraint efforts and material memory remain in their original coordinates.
use super::*;
use fs_la::pole_response::{PoleResponseError, PoleResponseOptions, PreparedPole};

mod derivatives;
pub use derivatives::{HarmonicParameterDirection, HarmonicPoleDerivative};

/// Numerical-owner limits plus physical mapping/output limits. Each limit
/// applies to preparation and independently to each complete batch. The
/// immutable source model has its own HarmonicBudget; no state is truncated.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicPoleResponseOptions {
    pub numerical: PoleResponseOptions,
    /// Conservative live physical payload, including the evaluated pencil,
    /// mapped RHS, cached mode observations and largest physical result batch.
    /// Excludes the separately budgeted fs-la factor, source model and allocator.
    pub max_physical_bytes: usize,
    /// Admission units: 32*(max_vectors+1)*(n*n+n*ports), not measured FLOPs.
    pub max_physical_work: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum HarmonicPoleResponseError {
    Invalid { what: &'static str },
    Physical(HarmonicPoleError),
    Numerical(PoleResponseError),
    Cancelled,
}
impl core::fmt::Display for HarmonicPoleResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physical pole response refused: {self:?}")
    }
}
impl std::error::Error for HarmonicPoleResponseError {}
impl From<PoleResponseError> for HarmonicPoleResponseError {
    fn from(e: PoleResponseError) -> Self { Self::Numerical(e) }
}
impl From<HarmonicPoleError> for HarmonicPoleResponseError {
    fn from(e: HarmonicPoleError) -> Self { Self::Physical(e) }
}
fn invalid(what: &'static str) -> HarmonicPoleResponseError {
    HarmonicPoleResponseError::Invalid { what }
}
fn check_cancel(c: &mut impl FnMut() -> bool) -> Result<(), HarmonicPoleResponseError> {
    if c() { Err(HarmonicPoleResponseError::Cancelled) } else { Ok(()) }
}
fn accumulate(sum: &mut C64, correction: &mut C64, term: C64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

/// Laurent coefficients for the supplied analytic port excitation evaluated
/// at s0: X(s)=state/(s-s0)+regular, E(s)=effort/(s-s0)+regular,
/// Y(s)=output/(s-s0)+regular. These are NOT harmonic amplitudes or powers.
/// In the exp(-i*omega*t) chart, omega0=i*s0 and the corresponding residues
/// are i times these coefficients, since s-s0=-i*(omega-omega0).
#[derive(Clone, Debug)]
pub struct HarmonicPortResidue {
    pub laplace_pole_per_s: C64,
    /// Differential state AND algebraic constraint multipliers, in source order.
    pub state: Vec<C64>,
    pub effort: Vec<C64>,
    pub output: Vec<C64>,
    /// state = prepared.right_state()*participation, in the fixed pivot gauge.
    pub participation: C64,
    pub linear_residual: f64,
}

/// One immutable factor bound to the original model, not a fitted oscillator.
/// It retains unobservable internal modes and supports singular descriptor E.
/// It cannot represent a repeated/defective pole, the polynomial/algebraic
/// feedthrough, or the complete transfer function by one mode's residue.
pub struct PreparedHarmonicPole<'a> {
    model: &'a HarmonicModel,
    numerical: PreparedPole,
    row_scales: Vec<f64>,
    effort: Vec<C64>,
    output: Vec<C64>,
    options: HarmonicPoleResponseOptions,
    physical_residual: f64,
}

fn admit(model: &HarmonicModel, o: &HarmonicPoleResponseOptions)
    -> Result<(), HarmonicPoleResponseError>
{
    let (n, m, k) = (model.n as u128, model.m as u128, o.numerical.max_vectors as u128);
    if !(1..=256).contains(&o.numerical.max_vectors) || model.n > o.numerical.max_dimension
        || !o.numerical.spectral_scale.is_finite() || o.numerical.spectral_scale <= 0.0 {
        return Err(invalid("complete model and bounded batch must fit the numerical limits"));
    }
    // Inputs and observations are counted even when their physical value is zero.
    let bytes = 16*(2*n*n + k*(5*n+2*m) + 8*n+2*m) + 256*k;
    let work = 32*(k+1)*(n*n+n*m);
    if bytes > o.max_physical_bytes as u128 || work > o.max_physical_work as u128 {
        return Err(invalid("complete physical mapping exceeds its payload or work budget"));
    }
    Ok(())
}

impl HarmonicModel {
    /// Prepare port residues at a previously resolved finite simple pole.
    /// The supplied state and pole are rechecked against THIS model, not trusted
    /// because their cached effort, energy, or residual says they are valid.
    /// Row conditioning is shared with local and contour mode discovery; source
    /// forcing is scaled by the same fixed rows. E is never inverted.
    pub fn prepare_pole_response(
        &self, pole: &HarmonicPole, options: HarmonicPoleResponseOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicPole<'_>, HarmonicPoleResponseError> {
        check_cancel(&mut cancelled)?;
        admit(self, &options)?;
        if pole.state.len() != self.n || pole.state.iter().any(|&z| !finite(z))
            || !finite(pole.laplace_pole_per_s) {
            return Err(invalid("complete finite pole state in the original model coordinates required"));
        }
        let mut pencil = PolePencil::new(self, options.numerical.spectral_scale);
        let evaluated = pencil.evaluate(pole.laplace_pole_per_s, &mut cancelled)
            .map_err(HarmonicPoleError::from)?;
        let rows = pencil.conditioning.as_ref().expect("evaluated pencil").0.clone();
        let numerical = PreparedPole::new(pole.laplace_pole_per_s, &pole.state,
            &evaluated, options.numerical, &mut cancelled)?;
        // The shared physical observer checks every original split descriptor
        // row, independently of fs-la's normwise acceptance and the fixed rows.
        let physical = self.observe_pole(NonlinearPole {
            value: numerical.value(), vector: numerical.right_vector().to_vec(),
            residual: numerical.pole_residual(), relative_step: numerical.relative_correction(),
            residual_scale: evaluated.residual_scale, iterations: 0,
        }, &mut cancelled)?;
        check_cancel(&mut cancelled)?;
        Ok(PreparedHarmonicPole {
            model: self, numerical, row_scales: rows, effort: physical.effort,
            output: physical.port_output, options, physical_residual: physical.relative_residual,
        })
    }
}

impl PreparedHarmonicPole<'_> {
    pub fn laplace_pole_per_s(&self) -> C64 { self.numerical.value() }
    /// Fixed normalization v[pivot]=1, not unit energy or a moving norm/phase.
    pub fn right_state(&self) -> &[C64] { self.numerical.right_vector() }
    pub fn mode_effort(&self) -> &[C64] { &self.effort }
    pub fn mode_output(&self) -> &[C64] { &self.output }
    pub fn normalization_pivot(&self) -> usize { self.numerical.normalization_pivot() }
    pub fn physical_residual(&self) -> f64 { self.physical_residual }

    /// Complete port-excitation batches, in original source port order. One
    /// cached factor serves all calls. Complex phase and nonnormal participation
    /// are retained. Zero participation is allowed; no observability cutoff.
    /// No finite-frequency response, mode truncation or artificial damping is
    /// inferred. Failed/cancelled batches publish nothing and leave the factor
    /// reusable. Mapping costs O(n*ports), then the shared O(n^2) bordered solve.
    pub fn port_residues(
        &self, inputs: &[&[C64]], mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<HarmonicPortResidue>, HarmonicPoleResponseError> {
        check_cancel(&mut cancelled)?;
        if inputs.is_empty() || inputs.len() > self.options.numerical.max_vectors
            || inputs.iter().any(|u| u.len() != self.model.m || u.iter().any(|&z| !finite(z))) {
            return Err(invalid("finite complete port batch within the declared ceiling required"));
        }
        let (n, m) = (self.model.n, self.model.m);
        let mut mapped = Vec::with_capacity(inputs.len());
        for &u in inputs {
            let mut rhs = vec![C64::ZERO; n];
            for i in 0..n {
                check_cancel(&mut cancelled)?;
                let mut correction = C64::ZERO;
                for p in 0..m {
                    accumulate(&mut rhs[i], &mut correction, u[p].scale(self.model.g[i*m+p]));
                }
                rhs[i] = divide(rhs[i], self.row_scales[i]);
            }
            if rhs.iter().any(|&z| !finite(z)) { return Err(invalid("port forcing overflow")); }
            mapped.push(rhs);
        }
        let refs: Vec<_> = mapped.iter().map(Vec::as_slice).collect();
        let residues = self.numerical.forcing_residues(&refs, &mut cancelled)?;
        let mut results = Vec::with_capacity(residues.len());
        for r in residues {
            check_cancel(&mut cancelled)?;
            // These observations are linear and were prepared in exactly the
            // numerical owner's gauge; no per-forcing storage projection is needed.
            let effort: Vec<_> = self.effort.iter().map(|&z| z*r.coefficient).collect();
            let output: Vec<_> = self.output.iter().map(|&z| z*r.coefficient).collect();
            if effort.iter().chain(&output).any(|&z| !finite(z)) {
                return Err(invalid("physical residue observation overflow"));
            }
            results.push(HarmonicPortResidue {
                laplace_pole_per_s: self.numerical.value(), state: r.vector, effort, output,
                participation: r.coefficient, linear_residual: r.linear_residual,
            });
        }
        check_cancel(&mut cancelled)?;
        Ok(results)
    }
}
