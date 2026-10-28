//! Transpose of the existing gray-diffuse radiosity equations, at fixed F.
//! The primal remains GrayDiffuseEnclosure::solve. The small nonsymmetric
//! transpose uses fs-solver GMRES and checks the original equation residual.
//! No finite differences, inverse-emissivity formula or black-surface singularity.
use super::{GrayDiffuseEnclosure, STEFAN_BOLTZMANN_W_M2_K4};
use crate::ConductionError;
use fs_exec::Cx;
use fs_solver::{GmresState, LinearOp, norm2};

type Result<T> = std::result::Result<T, ConductionError>;

#[derive(Clone, Debug)]
struct Matrix(Vec<Vec<f64>>);
impl LinearOp for Matrix {
    fn n(&self) -> usize { self.0.len() }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (row, value) in self.0.iter().zip(y) {
            *value = row.iter().zip(x).fold(0.0, |s, (&a,&b)| a.mul_add(b,s));
        }
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.0);
        for (row, &weight) in self.0.iter().zip(x) {
            for (value, &a) in y.iter_mut().zip(row) { *value = a.mul_add(weight,*value); }
        }
    }
}

/// Fixed-geometry gray-diffuse radiosity linearization at the admitted primal state.
///
/// This owns only the surface exchange equations. A coupled thermal derivative
/// must also account for the feedback from these temperatures into the solid.
#[derive(Clone, Debug)]
pub struct EnclosureRadiosityLinearization {
    matrix: Matrix,
    factors: Vec<Vec<f64>>,
    temperature_slopes: Vec<f64>,
    log_emissivity_slopes: Vec<f64>,
}

/// Transpose contraction of net outward surface fluxes at fixed view factors.
#[derive(Clone, Debug)]
pub struct EnclosureRadiosityPullback {
    /// Derivatives with respect to surface mean temperatures [per K].
    pub temperatures: Vec<f64>,
    /// Derivatives with respect to log(emissivity), including reflection.
    pub log_emissivities: Vec<f64>,
    /// Original radiosity-transpose relative residual after physical rescaling.
    pub relative_residual: f64,
    /// Krylov iterations spent on the admitted transpose solve.
    pub iterations: usize,
}

impl EnclosureRadiosityLinearization {
    /// Bind the exact radiosity model and verify its equation and closed-energy
    /// residuals against the supplied positive patch heat allowance [W].
    ///
    /// # Errors
    /// Refuses invalid inputs, failed primal heat gates or cancellation.
    pub fn new(cx: &Cx<'_>, enclosure: &GrayDiffuseEnclosure,
        temperatures: &[f64], heat_tolerance_w: f64) -> Result<Self> {
        poll(cx)?;
        if !heat_tolerance_w.is_finite() || heat_tolerance_w <= 0.0 {
            return Err(error("positive finite radiosity heat tolerance required"));
        }
        let primal = enclosure.solve(cx,temperatures)?;
        let n = enclosure.surfaces().len();
        let area = enclosure.surfaces().iter().map(|s|s.area_m2()).fold(0.0_f64,f64::max);
        if checked(primal.linear_residual_max_w_m2*area)? > heat_tolerance_w
            || checked(primal.enclosure_energy_closure_w)?.abs() > heat_tolerance_w {
            return Err(error("radiosity derivative base failed the primal equation or heat gate"));
        }
        let factors = enclosure.view_factors().factors().to_vec();
        let mut matrix = vec![vec![0.0;n];n];
        let mut temperature_slopes = Vec::with_capacity(n);
        let mut log_emissivity_slopes = Vec::with_capacity(n);
        for (i,surface) in enclosure.surfaces().iter().enumerate() {
            poll(cx)?;
            let epsilon = surface.emissivity().value();
            let t = temperatures[i];
            for (j,&factor) in factors[i].iter().enumerate() {
                matrix[i][j] = (if i==j {1.0} else {0.0}) - (1.0-epsilon)*factor;
            }
            temperature_slopes.push(checked(4.0*epsilon*STEFAN_BOLTZMANN_W_M2_K4*t*t*t)?);
            // d(rhs-MJ)/dln(eps) = eps*(sigma*T^4 - F*J) = q.
            // The primal computes q with a reference-shifted solve, retaining
            // small exchange without subtracting large absolute radiances.
            log_emissivity_slopes.push(checked(primal.net_outward_flux_w_m2[i])?);
        }
        Ok(Self {matrix:Matrix(matrix),factors,temperature_slopes,log_emissivity_slopes})
    }

    /// Weights multiply net outward FLUXES (W/m2), not patch watts. A heat
    /// functional must first multiply its weights by the corresponding areas.
    ///
    /// # Errors
    /// Refuses invalid weights or budgets, cancellation, nonfinite arithmetic,
    /// and any solve that fails its original transpose residual gate.
    pub fn pullback(&self,cx:&Cx<'_>,flux_weights:&[f64],relative:f64,
        max_iterations:usize) -> Result<EnclosureRadiosityPullback> {
        poll(cx)?;
        let n = self.matrix.n();
        if flux_weights.len()!=n || !(relative.is_finite() && relative>0.0 && relative<1.0)
            || max_iterations==0 {
            return Err(error("radiosity pullback needs complete flux weights and explicit solver limits"));
        }
        let mut rhs = flux_weights.to_vec();
        for (i,&weight) in flux_weights.iter().enumerate() {
            checked(weight)?;
            for (j,&factor) in self.factors[i].iter().enumerate() {
                rhs[j] = checked(rhs[j]-factor*weight)?;
            }
        }
        // q=(I-F)J, hence M^T z=(I-F)^T w. Normalize before the norm so a
        // small physical gradient cannot underflow into a false zero solve.
        let scale = rhs.iter().map(|v|v.abs()).fold(0.0_f64,f64::max);
        if scale==0.0 {
            return Ok(EnclosureRadiosityPullback {temperatures:vec![0.0;n],log_emissivities:vec![0.0;n],
                relative_residual:0.0,iterations:0});
        }
        let normalized:Vec<_> = rhs.iter().map(|v|v/scale).collect();
        let denominator = norm2(&normalized);
        let mut state = GmresState::new(&normalized,n.min(60).min(max_iterations));
        let mut product = vec![0.0;n];
        let mut last = f64::INFINITY;
        while state.iters<max_iterations {
            poll(cx)?;
            let before = state.iters;
            state.restart = n.min(60).min(max_iterations-before);
            state.run(&self.matrix,&normalized,relative,1,true);
            poll(cx)?;
            for &value in &state.x { checked(value)?; }
            self.matrix.apply_transpose(&state.x,&mut product);
            let residual = normalized.iter().zip(&product)
                .map(|(&b,&a)|checked(b-a)).collect::<Result<Vec<_>>>()?;
            last = checked(norm2(&residual)/denominator)?;
            if last<=relative {
                let z = state.x.iter().map(|v|checked(v*scale)).collect::<Result<Vec<_>>>()?;
                // Physical rescaling can round subnormal adjoints. Verify the
                // actual returned vector, not only the normalized iterate.
                let rounded = z.iter().map(|v|checked(v/scale)).collect::<Result<Vec<_>>>()?;
                self.matrix.apply_transpose(&rounded,&mut product);
                let rounded_residual = normalized.iter().zip(&product)
                    .map(|(&b,&a)|checked(b-a)).collect::<Result<Vec<_>>>()?;
                let original_residual = checked(norm2(&rounded_residual)/denominator)?;
                if original_residual > relative {
                    return Err(error(format!("radiosity transpose physical rescaling failed its original relative residual gate: {original_residual} > {relative}; no derivative returned")));
                }
                let temperatures = z.iter().zip(&self.temperature_slopes)
                    .map(|(&a,&b)|checked(a*b)).collect::<Result<Vec<_>>>()?;
                let log_emissivities = z.iter().zip(&self.log_emissivity_slopes)
                    .map(|(&a,&b)|checked(a*b)).collect::<Result<Vec<_>>>()?;
                return Ok(EnclosureRadiosityPullback {temperatures,log_emissivities,relative_residual:original_residual,iterations:state.iters});
            }
            if state.iters==before { break; }
        }
        Err(error(format!("radiosity transpose exhausted {max_iterations} Krylov iterations; original relative residual {last}; no derivative returned")))
    }
}

fn checked(value:f64)->Result<f64> {
    if value.is_finite() {Ok(value)} else {Err(error("nonfinite radiosity derivative arithmetic"))}
}
fn poll(cx:&Cx<'_>)->Result<()> {
    cx.checkpoint().map_err(|_|ConductionError::Cancelled {stage:"enclosure-radiosity-adjoint",at:0})
}
fn error(what:impl Into<String>)->ConductionError {
    ConductionError::Radiation {surface:"enclosure-adjoint".into(),what:what.into(),
        fix:"retain the admitted view factors and provide finite objectives with adequate explicit linear budgets".into()}
}
