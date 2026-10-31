//! Mouth flow from the SAME staircase solve as the input impedance.
//!
//! Retain the map from the local forward pressure wave to every mouth-mode
//! volume flow while the original reflection recursion runs backward. Uniform
//! spans multiply by the same decaying exp(i k L); junctions solve pressure and
//! flow continuity in forward-wave coordinates. No growing inverse propagation
//! or separately solved scalar bore is introduced.
use super::*;

/// Loaded modal input response and its complete original mouth-flow map.
#[derive(Debug, Clone)]
pub struct MouthResponse {
    /// Input response, bit-identical to the corresponding impedance-only call.
    pub input: ModalResponse,
    /// Actual first staircase-station radius, not an extrapolated cone endpoint.
    pub inlet_radius_m: f64,
    /// Actual last staircase-station radius used by this discrete load.
    pub mouth_radius_m: f64,
    /// Row-major N x N map: U_mouth = flow_map U_inlet. U_0 is total
    /// volume flow; higher entries use the original wall-normalized modes.
    pub flow_map: Vec<C64>,
}
impl MouthResponse {
    /// Mouth modal flows per unit incident PLANE pressure at a positive real
    /// inlet reference impedance. Higher inlet modal flows are zero: this is
    /// the Z00 source closure used by the scalar multimodal line bank, not a
    /// matched termination of every higher inlet mode.
    pub fn per_incident_plane_pressure(&self, reference: f64) -> Result<Vec<C64>, DuctError> {
        let n = self.input.n_modes;
        if !reference.is_finite() || reference <= 0.0 || n == 0 || n > MAX_MODES
            || self.flow_map.len() != n * n
        {
            return Err(bad("mouth transfer needs a positive reference and a complete modal map"));
        }
        let denominator = self.input.plane_impedance + C64::from_re(reference);
        if !finite(denominator) || denominator.abs() == 0.0 {
            return Err(DuctError::Singular);
        }
        let scale = C64::from_re(2.0) / denominator;
        let out: Vec<_> = (0..n).map(|row| self.flow_map[row * n] * scale).collect();
        if out.iter().any(|&value| !finite(value)) { return Err(DuctError::Singular); }
        Ok(out)
    }
}
fn bad(what: &'static str) -> DuctError { DuctError::BadParameter { what } }
fn finite(z: C64) -> bool { z.re.is_finite() && z.im.is_finite() }

/// Analytic mouth load: closed/pressure-release conditions act on every mode;
/// scalar radiation retains the disclosed matched higher-mode approximation.
/// The existing geometry, mode, frequency and radiation refusals apply.
pub fn mouth_response(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel,
    termination: Termination, modes: usize, extra_slices: usize)
    -> Result<MouthResponse, DuctError>
{
    super::mm_evaluate(duct, gas, omega, loss, MmLoad::Analytic(termination), modes, extra_slices, true)
}

/// Strict table at its declared original mouth plane, radius and gas state.
/// No extrapolation, fallback, or static value is inferred from positive bins.
pub fn mouth_response_tabulated(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel,
    table: &crate::TabulatedLoad, modes: usize, extra_slices: usize)
    -> Result<MouthResponse, DuctError>
{
    let radius = duct.segments.last().ok_or(DuctError::EmptyDuct)?.outlet_radius();
    table.validate_mouth(duct.segments.iter().map(Segment::axial_length).sum(), radius, gas)?;
    super::mm_evaluate(duct, gas, omega, loss, MmLoad::Tabulated(table), modes, extra_slices, true)
}

/// Explicit finite passive plane-mouth impedance; higher modes remain matched.
/// A caller can supply a separately declared static load at its quasi-static
/// probe without extending a table below its positive-frequency support.
pub fn mouth_response_load(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel,
    impedance: C64, modes: usize, extra_slices: usize) -> Result<MouthResponse, DuctError>
{
    if !finite(impedance) || impedance.re < 0.0 {
        return Err(bad("modal mouth impedance must be finite with nonnegative resistance"));
    }
    super::mm_evaluate(duct, gas, omega, loss, MmLoad::Impedance(impedance), modes, extra_slices, true)
}

/// Complete flanged radiation matrix, with the original interior/exterior gases.
/// The input impedance and all observed modal flows use that same coupled load.
pub fn mouth_response_flanged(duct: &Duct, gas: &GasState, exterior: &GasState,
    omega: f64, loss: LossModel, modes: usize, extra_slices: usize)
    -> Result<MouthResponse, DuctError>
{
    super::mm_evaluate(duct, gas, omega, loss, MmLoad::Flanged(exterior), modes, extra_slices, true)
}

/// U = Zc^-1 (I-R) a at the terminal plane.
pub(super) fn terminal_flow(r: &[C64], zc: &[C64], n: usize) -> Vec<C64> {
    let mut map = r.iter().map(|r| -*r).collect::<Vec<_>>();
    for row in 0..n {
        map[row*n+row] = map[row*n+row] + C64::ONE;
        for col in 0..n { map[row*n+col] = map[row*n+col] / zc[row]; }
    }
    map
}

/// H_left = H_right E. Only bounded forward propagation is evaluated.
pub(super) fn propagate(map: &mut Option<Vec<C64>>, e: &[C64], n: usize) {
    if let Some(map) = map {
        for row in 0..n {
            for col in 0..n { map[row*n+col] = map[row*n+col] * e[col]; }
        }
    }
}

/// Transfer forward amplitudes across a load-dependent original junction.
/// Expansion: p_u=F p_d, U_d=F^T U_u. Contraction: p_d=F p_u,
/// U_u=F^T U_d. Eliminate incoming waves using b_d=R_d a_d, but do
/// not invert F or I-R_d. Isolated flow nodes therefore add no new pole.
pub(super) fn junction(map: &mut Option<Vec<C64>>, r: &[C64], z_down: &[C64],
    z_up: &[C64], f: &[C64], expansion: bool, n: usize) -> Result<(), DuctError>
{
    let Some(map) = map else { return Ok(()); };
    let mut p = r.to_vec();
    for i in 0..n { p[i*n+i] = p[i*n+i] + C64::ONE; }
    let y = terminal_flow(r, z_down, n);
    let ft = transpose(f, n);
    let (a, b) = if expansion {
        let mut weighted_f = f.to_vec();
        for row in 0..n {
            for col in 0..n { weighted_f[row*n+col] = weighted_f[row*n+col] / z_up[row]; }
        }
        let mut a = mat_mul(&mat_mul(&ft, &weighted_f, n), &p, n);
        for (a, y) in a.iter_mut().zip(&y) { *a = *a + *y; }
        let mut b = ft;
        for row in 0..n {
            for col in 0..n { b[row*n+col] = b[row*n+col].scale(2.0) / z_up[col]; }
        }
        (a, b)
    } else {
        let mut weighted_ft = ft;
        for row in 0..n {
            for col in 0..n { weighted_ft[row*n+col] = weighted_ft[row*n+col] * z_up[row]; }
        }
        let mut a = mat_mul(&mat_mul(f, &weighted_ft, n), &y, n);
        for (a, p) in a.iter_mut().zip(&p) { *a = *a + *p; }
        (a, f.iter().map(|v| v.scale(2.0)).collect())
    };
    // Shared small complex LU, fixed column order, no second physical solver.
    let lu = lu_complex(&a, n).map_err(|_| DuctError::Singular)?;
    let mut j = vec![C64::ZERO; n*n];
    for col in 0..n {
        let mut rhs: Vec<_> = (0..n).map(|row| b[row*n+col]).collect();
        lu.solve(&mut rhs);
        for row in 0..n { j[row*n+col] = rhs[row]; }
    }
    *map = mat_mul(map, &j, n);
    if map.iter().any(|&v| !finite(v)) { return Err(DuctError::Singular); }
    Ok(())
}

pub(super) fn finish(map: Option<Vec<C64>>, impedance: &[C64], zc: &[C64], n: usize)
    -> Result<Vec<C64>, DuctError>
{
    let Some(map) = map else { return Ok(Vec::new()); };
    // a_in = (Z_in + Zc) U_in / 2, in pressure-wave coordinates.
    let mut input = impedance.to_vec();
    for i in 0..n { input[i*n+i] = input[i*n+i] + zc[i]; }
    for v in &mut input { *v = v.scale(0.5); }
    let flow = mat_mul(&map, &input, n);
    if flow.iter().any(|&v| !finite(v)) { return Err(DuctError::Singular); }
    Ok(flow)
}

#[cfg(test)]
mod tests;
