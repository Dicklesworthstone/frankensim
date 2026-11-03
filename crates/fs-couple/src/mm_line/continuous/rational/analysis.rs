//! Original duct, rational memory and pressure continuity in one pencil.
//!
//! Unknowns are [duct q; duct M*P; original real load states; outward mouth U].
//! The last coordinates are zero-storage flow multipliers. With s=-i*omega,
//! the equations are x'=A_d*x+G_in*(2a-P_in)/Z-G_out*U,
//! z'=A_l*z+B_l*U, and P_out=C_l*z+D_l*U. Neither D nor a bare component's
//! impedance is inverted. This includes zero and mixed-rank direct terms.
//!
//! Generic rational states have no declared storage matrix. Their port work
//! stays signed exchange, never inferred heat. Poles describe the selected
//! rational continuation and finite spatial image, not exact exterior physics.

use super::{
    GasState, ModalGalerkinSource, ModalRationalLoad, RationalModalError, RationalModalLine,
};
use fs_duct::RadiusProfile;
use fs_duct::modal::time_domain::{ModalEnd, assemble_lossless};
use fs_exec::CancelGate;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_la::nonlinear_eig::beyn::{BeynError, BeynReport, beyn_with_cancel};
use fs_la::nonlinear_eig::{
    NonlinearEigError, NonlinearPole, PoleEvaluation, newton_pole_with_cancel,
};
use fs_math::{c64::C64, det};
use fs_phs::PortHamiltonian;
use fs_vfit::mimo::MimoStateSpace;

pub use fs_la::nonlinear_eig::PoleOptions;
pub use fs_la::nonlinear_eig::beyn::{BeynOptions, Ellipse as ResonanceRegion};

/// Resources for the complete physical pencil and one reusable harmonic LU.
/// Pole/contour numerical work has its own explicit shared-solver budget.
#[derive(Clone, Copy, Debug)]
pub struct RationalModalAnalysisBudget {
    pub maximum_dimension: usize,
    /// Conservative simultaneous numerical payload, excluding allocator/RSS.
    pub maximum_payload_bytes: usize,
    /// Bounds preparation, one factorization and original-row observations.
    pub maximum_work: u64,
    /// Fixed positive conditioning scale [1/s]; never a physical damping law.
    pub reference_rate_per_s: f64,
    /// Original, uneliminated rowwise backward-error tolerance in (0, 1e-8].
    /// Uses the fixed numerical coordinate scales, including for dark modes
    /// whose exactly zero physical components acquire roundoff in discovery.
    pub relative_residual_tolerance: f64,
}

#[derive(Debug)]
pub enum RationalModalAnalysisError {
    Source(RationalModalError),
    Invalid(&'static str),
    Linear(LuComplexError),
    Pole(NonlinearEigError),
    Contour(BeynError),
    Residual { relative: f64, tolerance: f64 },
    Cancelled,
}
impl core::fmt::Display for RationalModalAnalysisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint rational modal analysis: {self:?}")
    }
}
impl std::error::Error for RationalModalAnalysisError {}
impl From<RationalModalError> for RationalModalAnalysisError {
    fn from(value: RationalModalError) -> Self {
        Self::Source(value)
    }
}
fn invalid(what: &'static str) -> RationalModalAnalysisError {
    RationalModalAnalysisError::Invalid(what)
}
fn poll(gate: &CancelGate) -> Result<(), RationalModalAnalysisError> {
    if gate.is_requested() {
        Err(RationalModalAnalysisError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}
fn dot(row: &[f64], x: &[C64]) -> C64 {
    row.iter()
        .zip(x)
        .fold(C64::ZERO, |a, (r, x)| a + x.scale(*r))
}

/// Independently observed original fields and boundary exchange.
/// Harmonic queries use peak phasors and cycle means; pole queries use complex
/// envelopes, with twice the quadratic factor. No load storage is inferred.
#[derive(Clone, Debug)]
pub struct RationalModalObservation {
    /// Complete original state, then all zero-storage outward-flow multipliers.
    pub state: Vec<C64>,
    pub inlet_pressure_pa: Vec<C64>,
    pub inlet_flow_m3_s: Vec<C64>,
    pub reflected_pressure_pa: Vec<C64>,
    pub duct_mouth_pressure_pa: Vec<C64>,
    pub load_pressure_pa: Vec<C64>,
    pub outward_mouth_flow_m3_s: Vec<C64>,
    /// Duct storage only, not energy of the complete rational interconnection.
    pub duct_energy_j: f64,
    pub inlet_power_w: f64,
    pub incident_power_w: f64,
    pub reflected_power_w: f64,
    pub duct_mouth_power_w: f64,
    /// Signed exchange into the original rational realization, not heat.
    pub load_power_w: f64,
    /// 2*Re(s)*H_duct - inlet power + duct mouth power.
    pub duct_balance_residual_w: f64,
    /// Independently evaluated duct mouth power minus original load port work.
    pub interface_power_residual_w: f64,
    pub wave_power_residual_w: f64,
    /// Maximum original duct, memory and pressure-continuity rowwise backward
    /// error. Denominator uses |rhs| plus the SPLIT coefficient row envelope
    /// times ||C x||_inf, in the fixed numerical chart C. Opposing body and
    /// interface coefficients are not cancelled in this envelope; exactly zero
    /// physical blocks of dark modes therefore remain meaningful.
    pub relative_equation_residual: f64,
}

#[derive(Clone, Debug)]
pub struct RationalModalPole {
    pub laplace_pole_per_s: C64,
    pub angular_frequency_rad_s: C64,
    /// Whether |Im(s)| lies in the original real-frequency source band.
    /// This is not validation of the continuation away from the frequency axis.
    pub oscillation_within_source_band: bool,
    pub physical: RationalModalObservation,
    pub iterations: usize,
}

#[derive(Clone, Debug)]
pub struct RationalModalSpectrum {
    pub poles: Vec<RationalModalPole>,
    pub unrefined_laplace_poles_per_s: Vec<C64>,
    pub argument_count: usize,
    pub resolved_rank: usize,
    pub quadrature_nodes: usize,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
}

/// Borrowed source identity and immutable full original descriptor coefficients.
/// No accepted history, clock or fit is read into the continuous pencil.
pub struct PreparedRationalModalAnalysis<'a> {
    source: ModalGalerkinSource,
    boundary: &'a ModalRationalLoad,
    nd: usize,
    n: usize,
    q: Vec<f64>,
    duct_a: Vec<f64>,
    g: Vec<f64>,
    pressure: Vec<Vec<f64>>,
    load: MimoStateSpace,
    generator: Vec<f64>,
    rows: Vec<f64>,
    columns: Vec<f64>,
    original_row_envelopes: Vec<f64>,
    a_norm: f64,
    e_norm: f64,
    budget: RationalModalAnalysisBudget,
}

pub struct RationalModalFactor<'m, 'a> {
    model: &'m PreparedRationalModalAnalysis<'a>,
    omega: f64,
    lu: LuComplex,
}

impl RationalModalLine {
    /// Analyze the existing original sources without consulting or changing
    /// accepted state or held-wave coordinates. No time-step warping is used.
    pub fn prepare_analysis(
        &self,
        budget: RationalModalAnalysisBudget,
        gate: &CancelGate,
    ) -> Result<PreparedRationalModalAnalysis<'_>, RationalModalAnalysisError> {
        PreparedRationalModalAnalysis::prepare(
            self.source.clone(),
            &self.duct,
            &self.boundary,
            budget,
            gate,
        )
    }
}

impl<'a> PreparedRationalModalAnalysis<'a> {
    /// Prepare directly from physical sources, including mixed-rank D models
    /// which the exact time-domain elimination route cannot currently evolve.
    /// Every original state and every mouth multiplier is retained.
    pub fn new(
        profile: &RadiusProfile,
        gas: &GasState,
        elements: usize,
        boundary: &'a ModalRationalLoad,
        budget: RationalModalAnalysisBudget,
        gate: &CancelGate,
    ) -> Result<Self, RationalModalAnalysisError> {
        poll(gate)?;
        let p = boundary.load.ports();
        let nd = elements
            .checked_add(1)
            .and_then(|v| v.checked_mul(2 * p))
            .ok_or_else(|| invalid("original modal spatial extent overflow"))?;
        Self::admit(nd, boundary, budget)?;
        let radius = profile.outlet_radius();
        if (boundary.radius_m - radius).abs() > 128.0 * f64::EPSILON * radius
            || boundary.exterior.pressure.to_bits() != gas.pressure.to_bits()
        {
            return Err(invalid("original mouth radius or static pressure mismatch"));
        }
        let (source, duct) = assemble_lossless(
            profile,
            gas,
            boundary.modes(),
            elements,
            ModalEnd::Ports,
            budget.maximum_dimension,
            || gate.is_requested(),
        )
        .map_err(RationalModalError::from)?
        .into_parts();
        Self::prepare(source, &duct, boundary, budget, gate)
    }

    fn admit(
        nd: usize,
        boundary: &ModalRationalLoad,
        budget: RationalModalAnalysisBudget,
    ) -> Result<usize, RationalModalAnalysisError> {
        let p = boundary.load.ports();
        let n = nd
            .checked_add(boundary.load.continuous_operator().dimension())
            .and_then(|v| v.checked_add(p))
            .ok_or_else(|| invalid("joint extent overflow"))?;
        let d = n as u128;
        let bytes = 160 * d * d + 128 * d * p as u128 + 512 * d;
        let work = 4 * d * d * d + 16 * d * d * (p as u128 + 1);
        if n == 0
            || n > 1024
            || n > budget.maximum_dimension
            || bytes > budget.maximum_payload_bytes as u128
            || work > budget.maximum_work as u128
            || !budget.reference_rate_per_s.is_finite()
            || budget.reference_rate_per_s <= 0.0
            || !budget.relative_residual_tolerance.is_finite()
            || budget.relative_residual_tolerance <= 0.0
            || budget.relative_residual_tolerance > 1e-8
        {
            return Err(invalid(
                "complete joint descriptor exceeds declared resources or accuracy",
            ));
        }
        Ok(n)
    }

    fn prepare(
        source: ModalGalerkinSource,
        duct: &PortHamiltonian,
        boundary: &'a ModalRationalLoad,
        budget: RationalModalAnalysisBudget,
        gate: &CancelGate,
    ) -> Result<Self, RationalModalAnalysisError> {
        poll(gate)?;
        let (nd, p) = (duct.state_dim(), boundary.load.ports());
        let n = Self::admit(nd, boundary, budget)?;
        let load = boundary
            .load
            .continuous_operator()
            .state_space(budget.maximum_payload_bytes, &mut || gate.is_requested())
            .map_err(RationalModalError::from)?;
        let nt = load.n;
        let diff = nd + nt;
        let part = duct
            .storage()
            .quadratic_part()
            .ok_or_else(|| invalid("original duct has no complete quadratic storage"))?;
        if !part.exact || part.matrix.len() != nd * nd || duct.port_dim() != 2 * p {
            return Err(invalid(
                "complete original duct storage and both modal faces required",
            ));
        }
        let q = part.matrix;
        let (j, r, g) = duct.structure();
        if r.iter().any(|&v| v != 0.0) {
            return Err(invalid(
                "this original modal duct image requires its lossless source",
            ));
        }
        let mut duct_a = vec![0.0; nd * nd];
        let mut pressure = vec![vec![0.0; nd]; 2 * p];
        for i in 0..nd {
            poll(gate)?;
            for c in 0..nd {
                duct_a[i * nd + c] = (0..nd).map(|k| j[i * nd + k] * q[k * nd + c]).sum();
            }
            for k in 0..2 * p {
                pressure[k][i] = (0..nd).map(|a| g[a * 2 * p + k] * q[a * nd + i]).sum();
            }
        }
        let z = source.inlet_reference_impedance();
        let mut generator = vec![0.0; n * n];
        for i in 0..nd {
            poll(gate)?;
            for c in 0..nd {
                generator[i * n + c] = duct_a[i * nd + c]
                    - (0..p)
                        .map(|k| g[i * 2 * p + k] * pressure[k][c] / z)
                        .sum::<f64>();
            }
            for k in 0..p {
                generator[i * n + diff + k] = -g[i * 2 * p + p + k];
            }
        }
        for i in 0..nt {
            poll(gate)?;
            generator[(nd + i) * n + nd..(nd + i) * n + diff]
                .copy_from_slice(&load.a[i * nt..(i + 1) * nt]);
            generator[(nd + i) * n + diff..(nd + i + 1) * n]
                .copy_from_slice(&load.b[i * p..(i + 1) * p]);
        }
        for k in 0..p {
            let row = (diff + k) * n;
            generator[row..row + nd].copy_from_slice(&pressure[p + k]);
            for c in 0..nt {
                generator[row + nd + c] = -load.c[k * nt + c];
            }
            for c in 0..p {
                generator[row + diff + c] = -load.d[k * p + c];
            }
        }
        if generator
            .iter()
            .chain(&q)
            .chain(&duct_a)
            .chain(pressure.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(invalid("original joint coefficients are not representable"));
        }
        // Constant two-sided equilibration keeps the pencil analytic. Columns
        // are only a reversible numerical chart; returned states stay original.
        let mut rows = vec![0.0; n];
        for i in 0..n {
            rows[i] = generator[i * n..(i + 1) * n].iter().map(|v| v.abs()).fold(
                if i < diff {
                    budget.reference_rate_per_s
                } else {
                    0.0
                },
                f64::max,
            );
            if rows[i] <= 0.0 {
                return Err(invalid("joint descriptor has an unresolved zero row"));
            }
        }
        let mut columns = vec![0.0; n];
        for c in 0..n {
            poll(gate)?;
            columns[c] = (0..n).map(|i| generator[i * n + c].abs() / rows[i]).fold(
                if c < diff {
                    budget.reference_rate_per_s / rows[c]
                } else {
                    0.0
                },
                f64::max,
            );
            if columns[c] <= 0.0 || !columns[c].is_finite() {
                return Err(invalid("joint descriptor has an unresolved coordinate"));
            }
        }
        let mut a_norm = 0.0;
        let mut e_norm = 0.0;
        let pressure_envelopes: Vec<_> = pressure
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&columns)
                    .map(|(v, c)| v.abs() / c)
                    .sum::<f64>()
            })
            .collect();
        let mut original_row_envelopes = vec![0.0; n];
        for i in 0..n {
            poll(gate)?;
            for c in 0..n {
                a_norm = det::hypot(a_norm, generator[i * n + c] / rows[i] / columns[c]);
            }
            if i < diff {
                e_norm = det::hypot(e_norm, (1.0 / rows[i]) / columns[i]);
            }
            original_row_envelopes[i] = if i < nd {
                (0..nd)
                    .map(|c| duct_a[i * nd + c].abs() / columns[c])
                    .sum::<f64>()
                    + (0..p)
                        .map(|k| {
                            g[i * 2 * p + k].abs() / z * pressure_envelopes[k]
                                + g[i * 2 * p + p + k].abs() / columns[diff + k]
                        })
                        .sum::<f64>()
            } else if i < diff {
                let row = i - nd;
                (0..nt)
                    .map(|c| load.a[row * nt + c].abs() / columns[nd + c])
                    .sum::<f64>()
                    + (0..p)
                        .map(|k| load.b[row * p + k].abs() / columns[diff + k])
                        .sum::<f64>()
            } else {
                let k = i - diff;
                pressure_envelopes[p + k]
                    + (0..nt)
                        .map(|c| load.c[k * nt + c].abs() / columns[nd + c])
                        .sum::<f64>()
                    + (0..p)
                        .map(|c| load.d[k * p + c].abs() / columns[diff + c])
                        .sum::<f64>()
            };
        }
        if !a_norm.is_finite()
            || !e_norm.is_finite()
            || e_norm <= 0.0
            || original_row_envelopes.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("joint conditioning is not representable"));
        }
        Ok(Self {
            source,
            boundary,
            nd,
            n,
            q,
            duct_a,
            g: g.to_vec(),
            pressure,
            load,
            generator,
            rows,
            columns,
            original_row_envelopes,
            a_norm,
            e_norm,
            budget,
        })
    }

    pub fn source(&self) -> &ModalGalerkinSource {
        &self.source
    }
    pub fn boundary(&self) -> &ModalRationalLoad {
        self.boundary
    }
    pub fn dimension(&self) -> usize {
        self.n
    }
    pub fn duct_dimension(&self) -> usize {
        self.nd
    }
    pub fn load_dimension(&self) -> usize {
        self.load.n
    }
    pub fn port_count(&self) -> usize {
        self.load.ports
    }

    fn evaluate(
        &self,
        s: C64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<PoleEvaluation, NonlinearEigError> {
        let n = self.n;
        let mut matrix = vec![C64::ZERO; n * n];
        let mut derivative = vec![C64::ZERO; n * n];
        for i in 0..n {
            if cancelled() {
                return Err(NonlinearEigError::Cancelled);
            }
            for c in 0..n {
                matrix[i * n + c] =
                    C64::from_re(-self.generator[i * n + c] / self.rows[i] / self.columns[c]);
            }
            if i < self.nd + self.load.n {
                let e = (1.0 / self.rows[i]) / self.columns[i];
                matrix[i * n + i] = matrix[i * n + i] + s.scale(e);
                derivative[i * n + i] = C64::from_re(e);
            }
        }
        let scale = self.a_norm + s.abs().max(self.budget.reference_rate_per_s) * self.e_norm;
        if !scale.is_finite() || matrix.iter().any(|&v| !finite(v)) {
            return Err(NonlinearEigError::FunctionRefused { at: s });
        }
        Ok(PoleEvaluation {
            matrix,
            derivative,
            residual_scale: scale,
        })
    }

    /// One shared LU at an original continuous frequency, reusable for every
    /// incident modal excitation. The original source frequency band applies.
    pub fn factor(
        &self,
        omega_rad_s: f64,
        gate: &CancelGate,
    ) -> Result<RationalModalFactor<'_, 'a>, RationalModalAnalysisError> {
        poll(gate)?;
        self.boundary
            .load
            .continuous_operator()
            .admit_frequency(omega_rad_s)
            .map_err(RationalModalError::from)?;
        let image = self
            .evaluate(C64::new(0.0, -omega_rad_s), &mut || gate.is_requested())
            .map_err(RationalModalAnalysisError::Pole)?;
        let lu = lu_complex_with_cancel(&image.matrix, self.n, || gate.is_requested())
            .map_err(RationalModalAnalysisError::Linear)?;
        Ok(RationalModalFactor {
            model: self,
            omega: omega_rad_s,
            lu,
        })
    }

    /// Refine an isolated finite mode from a complete ORIGINAL-coordinate seed.
    /// Zero incident inlet waves impose the same original matched reference.
    /// The global rational continuation is selected explicitly by its source;
    /// a band's in-range flag does not establish a complex-frequency fit error.
    pub fn resonance_pole(
        &self,
        laplace_seed_per_s: C64,
        state_seed: &[C64],
        options: &PoleOptions,
        gate: &CancelGate,
    ) -> Result<RationalModalPole, RationalModalAnalysisError> {
        poll(gate)?;
        if state_seed.len() != self.n || state_seed.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "one finite seed for every original joint coordinate required",
            ));
        }
        let seed: Vec<_> = state_seed
            .iter()
            .zip(&self.columns)
            .map(|(x, c)| x.scale(*c))
            .collect();
        let pole = newton_pole_with_cancel(
            self.n,
            laplace_seed_per_s,
            &seed,
            options,
            |s, c| self.evaluate(s, c),
            || gate.is_requested(),
        )
        .map_err(RationalModalAnalysisError::Pole)?;
        self.observe_pole(pole, gate)
    }

    /// Shared Beyn discovery/refinement of the same entire descriptor pencil.
    /// Finite counts are numerical evidence; infinite roots, unresolved
    /// clusters and source approximation errors are not certified by this API.
    pub fn resonance_spectrum(
        &self,
        region: ResonanceRegion,
        options: &BeynOptions,
        gate: &CancelGate,
    ) -> Result<RationalModalSpectrum, RationalModalAnalysisError> {
        poll(gate)?;
        let n = self.n as u128;
        let modes = options.probe_columns as u128;
        let bytes = 64 * modes * n + 512 * (n + modes);
        let work = 16 * modes * n * n;
        if modes == 0
            || modes > n
            || bytes > options.max_dense_bytes as u128
            || work > options.max_work as u128
        {
            return Err(invalid(
                "original joint spectrum observations exceed their resources",
            ));
        }
        let mut numerical = *options;
        numerical.max_dense_bytes -= bytes as usize;
        numerical.max_work -= work as u64;
        let BeynReport {
            poles,
            unrefined_values,
            argument_integral,
            coarse_argument_integral,
            argument_count,
            rank,
            nodes,
            moment_relative_change,
            minimum_reciprocal_condition,
        } = beyn_with_cancel(
            self.n,
            region,
            &numerical,
            |s, c| self.evaluate(s, c),
            || gate.is_requested(),
        )
        .map_err(RationalModalAnalysisError::Contour)?;
        let mut observed = Vec::with_capacity(poles.len());
        for pole in poles {
            observed.push(self.observe_pole(pole, gate)?);
        }
        poll(gate)?;
        Ok(RationalModalSpectrum {
            poles: observed,
            unrefined_laplace_poles_per_s: unrefined_values,
            argument_count,
            resolved_rank: rank,
            quadrature_nodes: nodes,
            argument_integral,
            coarse_argument_integral,
            moment_relative_change,
            minimum_reciprocal_condition,
        })
    }

    fn observe_pole(
        &self,
        pole: NonlinearPole,
        gate: &CancelGate,
    ) -> Result<RationalModalPole, RationalModalAnalysisError> {
        let mut state: Vec<_> = pole
            .vector
            .iter()
            .zip(&self.columns)
            .map(|(x, c)| divide(*x, *c))
            .collect();
        let norm = state.iter().fold(0.0, |a, x| det::hypot(a, x.abs()));
        if !norm.is_finite() || norm <= 0.0 {
            return Err(invalid("original pole normalization overflow"));
        }
        for x in &mut state {
            *x = divide(*x, norm);
        }
        let s = pole.value;
        let band = self.boundary.load.limits().band_rad_s;
        let physical = self.observe(s, state, &vec![C64::ZERO; self.port_count()], 1.0, gate)?;
        Ok(RationalModalPole {
            laplace_pole_per_s: s,
            angular_frequency_rad_s: C64::new(-s.im, s.re),
            oscillation_within_source_band: s.im.abs() >= band[0] && s.im.abs() <= band[1],
            physical,
            iterations: pole.iterations,
        })
    }

    fn observe(
        &self,
        s: C64,
        state: Vec<C64>,
        incident: &[C64],
        factor: f64,
        gate: &CancelGate,
    ) -> Result<RationalModalObservation, RationalModalAnalysisError> {
        let (nd, nt, p) = (self.nd, self.load.n, self.port_count());
        let diff = nd + nt;
        let z = self.source.inlet_reference_impedance();
        let inlet: Vec<_> = self.pressure[..p]
            .iter()
            .map(|row| dot(row, &state[..nd]))
            .collect();
        let mouth: Vec<_> = self.pressure[p..]
            .iter()
            .map(|row| dot(row, &state[..nd]))
            .collect();
        let flow: Vec<_> = inlet
            .iter()
            .zip(incident)
            .map(|(y, a)| (a.scale(2.0) - *y).scale(1.0 / z))
            .collect();
        let reflected: Vec<_> = inlet.iter().zip(incident).map(|(y, a)| *y - *a).collect();
        let outward = state[diff..].to_vec();
        let load: Vec<_> = (0..p)
            .map(|k| {
                dot(&self.load.c[k * nt..(k + 1) * nt], &state[nd..diff])
                    + dot(&self.load.d[k * p..(k + 1) * p], &outward)
            })
            .collect();
        let mut relative = 0.0_f64;
        let state_norm = state
            .iter()
            .zip(&self.columns)
            .map(|(x, c)| x.abs() * c)
            .fold(0.0, f64::max);
        if !state_norm.is_finite() {
            return Err(invalid("original joint norm overflow"));
        }
        // Original split equations, independently of the closed/scaled pencil.
        for i in 0..self.n {
            poll(gate)?;
            let mut residual = C64::ZERO;
            let mut denominator = 0.0;
            let mut add = |v: C64| {
                residual = residual + v;
                denominator += v.abs();
            };
            if i < nd {
                add(s * state[i]);
                for (c, x) in state[..nd].iter().enumerate() {
                    add(x.scale(-self.duct_a[i * nd + c]));
                }
                for k in 0..p {
                    add(flow[k].scale(-self.g[i * 2 * p + k]));
                    add(outward[k].scale(self.g[i * 2 * p + p + k]));
                }
            } else if i < diff {
                let row = i - nd;
                add(s * state[i]);
                for c in 0..nt {
                    add(state[nd + c].scale(-self.load.a[row * nt + c]));
                }
                for k in 0..p {
                    add(outward[k].scale(-self.load.b[row * p + k]));
                }
            } else {
                let k = i - diff;
                // Keep the actual original pressure terms before cancellation.
                for c in 0..nd {
                    add(state[c].scale(-self.pressure[p + k][c]));
                }
                for c in 0..nt {
                    add(state[nd + c].scale(self.load.c[k * nt + c]));
                }
                for c in 0..p {
                    add(outward[c].scale(self.load.d[k * p + c]));
                }
            }
            if !denominator.is_finite() || !finite(residual) {
                return Err(invalid("original joint row accumulation overflow"));
            }
            // A componentwise denominator alone is undefined in an exactly
            // zero block of a dark mode: discovery leaves harmless roundoff
            // there. Use the standard original rowwise backward error in the
            // fixed reversible chart, without projecting any coordinate.
            let row_norm = self.original_row_envelopes[i]
                + if i < diff {
                    s.abs() / self.columns[i]
                } else {
                    0.0
                };
            let rhs = if i < nd {
                (0..p)
                    .fold(C64::ZERO, |a, k| {
                        a + incident[k].scale(2.0 * self.g[i * 2 * p + k] / z)
                    })
                    .abs()
            } else {
                0.0
            };
            let norm_denominator = rhs + row_norm * state_norm;
            if !norm_denominator.is_finite() {
                return Err(invalid("original joint backward-error scale overflow"));
            }
            denominator = denominator.max(norm_denominator);
            let row = if denominator == 0.0 {
                if residual == C64::ZERO {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                residual.abs() / denominator
            };
            if !row.is_finite() {
                return Err(invalid("original joint residual is not representable"));
            }
            relative = relative.max(row);
        }
        if relative > self.budget.relative_residual_tolerance {
            return Err(RationalModalAnalysisError::Residual {
                relative,
                tolerance: self.budget.relative_residual_tolerance,
            });
        }
        let mut energy = 0.0;
        for i in 0..nd {
            poll(gate)?;
            energy += (state[i].conj() * dot(&self.q[i * nd..(i + 1) * nd], &state[..nd])).re
                * 0.5
                * factor;
        }
        let power = |a: &[C64], b: &[C64]| {
            a.iter()
                .zip(b)
                .map(|(a, b)| (a.conj() * *b).re * factor)
                .sum::<f64>()
        };
        let inlet_power = power(&inlet, &flow);
        let mouth_power = power(&mouth, &outward);
        let load_power = power(&load, &outward);
        let incident_power = incident
            .iter()
            .map(|a| a.norm_sq() * factor / z)
            .sum::<f64>();
        let reflected_power = reflected
            .iter()
            .map(|b| b.norm_sq() * factor / z)
            .sum::<f64>();
        let balance = 2.0 * s.re * energy - inlet_power + mouth_power;
        let interface = mouth_power - load_power;
        let wave = incident_power - reflected_power - inlet_power;
        if state
            .iter()
            .chain(&inlet)
            .chain(&flow)
            .chain(&reflected)
            .chain(&mouth)
            .chain(&load)
            .any(|&v| !finite(v))
            || [
                energy,
                inlet_power,
                mouth_power,
                load_power,
                incident_power,
                reflected_power,
                balance,
                interface,
                wave,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "original joint field or power is not representable",
            ));
        }
        Ok(RationalModalObservation {
            state,
            inlet_pressure_pa: inlet,
            inlet_flow_m3_s: flow,
            reflected_pressure_pa: reflected,
            duct_mouth_pressure_pa: mouth,
            load_pressure_pa: load,
            outward_mouth_flow_m3_s: outward,
            duct_energy_j: energy,
            inlet_power_w: inlet_power,
            incident_power_w: incident_power,
            reflected_power_w: reflected_power,
            duct_mouth_power_w: mouth_power,
            load_power_w: load_power,
            duct_balance_residual_w: balance,
            interface_power_residual_w: interface,
            wave_power_residual_w: wave,
            relative_equation_residual: relative,
        })
    }
}

impl RationalModalFactor<'_, '_> {
    /// Peak exp(-i*omega*t) incident pressure in every original radial mode.
    /// This factor can serve many drives; every refusal leaves it reusable.
    pub fn response(
        &self,
        incident: &[C64],
        gate: &CancelGate,
    ) -> Result<RationalModalObservation, RationalModalAnalysisError> {
        poll(gate)?;
        let model = self.model;
        let p = model.port_count();
        if incident.len() != p || incident.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "one finite incident pressure per original inlet mode required",
            ));
        }
        let z = model.source.inlet_reference_impedance();
        let mut state = vec![C64::ZERO; model.n];
        for i in 0..model.nd {
            poll(gate)?;
            state[i] = (0..p).fold(C64::ZERO, |a, k| {
                a + incident[k].scale(2.0 * model.g[i * 2 * p + k] / z / model.rows[i])
            });
        }
        if state.iter().any(|&v| !finite(v)) {
            return Err(invalid("joint incident drive overflow"));
        }
        self.lu
            .solve_with_cancel(&mut state, || gate.is_requested())
            .map_err(RationalModalAnalysisError::Linear)?;
        for (x, c) in state.iter_mut().zip(&model.columns) {
            *x = divide(*x, *c);
        }
        model.observe(C64::new(0.0, -self.omega), state, incident, 0.5, gate)
    }
}
