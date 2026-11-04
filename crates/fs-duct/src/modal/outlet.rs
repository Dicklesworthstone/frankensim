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

/// Input reflection and mouth flow in the original native incident waves.
///
/// The mouth map is retained before input-impedance conversion. Its unit
/// incident-pressure responses remain finite at inlet flow nodes even when
/// `U_mouth / U_inlet` and the input impedance have genuine poles.
#[derive(Debug, Clone)]
pub struct MouthWaveResponse {
    /// Same native input reflection as the corresponding prepared query.
    pub input: ModalReflectionResponse,
    /// Actual first staircase-station radius [m].
    pub inlet_radius_m: f64,
    /// Actual last staircase-station radius [m].
    pub mouth_radius_m: f64,
    /// Row-major N x N map `U_mouth = H a_native`, in m^3/(Pa s).
    /// Positive mouth flow points out of the duct. Modes use the original
    /// signed wall-normalized radial basis, including evanescent modes.
    pub native_flow_map: Vec<C64>,
}

/// One coherent modal source and its observed physical mouth flows.
#[derive(Debug, Clone)]
pub struct MouthExcitation {
    /// Reflected inlet pressure waves [Pa] at the supplied real references.
    pub reflected_pressure: Vec<C64>,
    /// Outward mouth modal volume flows [m^3/s] from that same source solve.
    pub mouth_flow: Vec<C64>,
}

/// Coherent physical source boundaries and their original mouth observation.
#[derive(Debug, Clone)]
pub struct MouthSourceResponse {
    /// Inlet pressures, flows, generator work and source resistance losses.
    pub inlet: ModalSourceResponse,
    /// Outward mouth modal volume flows [m³/s] from the same source solve.
    pub mouth_flow: Vec<C64>,
}

impl MouthWaveResponse {
    /// Mouth flows per unit incident plane pressure, with all higher inlet
    /// volume flows zero. This is the original `Z00` source closure; no input
    /// impedance or mouth/inlet-flow ratio is formed.
    ///
    /// # Errors
    /// Invalid wave references/maps or a singular plane-source boundary.
    pub fn per_incident_plane_pressure(&self, reference: f64) -> Result<Vec<C64>, DuctError> {
        let (_, native) = self.input.plane_excitation(reference)?;
        self.observe(&native)
    }

    /// Excite every retained incident port coherently at positive real
    /// per-mode reference impedances [Pa s/m^3]. A zero incident higher-mode
    /// entry matches that source port, rather than forcing its flow to zero.
    ///
    /// The one native-wave solve supplies both reflection and mouth flow;
    /// interference and mutual radiation terms retain their complex phases.
    /// Input mean power is `sum((|a|^2-|b|^2)/reference)/2` for peak phasors.
    ///
    /// # Errors
    /// Invalid reference/incident/map dimensions, nonfinite data or a singular
    /// source boundary. The original load and propagation limits still apply.
    pub fn excite(&self, references: &[f64], incident_pressure: &[C64])
        -> Result<MouthExcitation, DuctError>
    {
        let (reflected_pressure, native) = self.input.excitation(references, incident_pressure)?;
        Ok(MouthExcitation { reflected_pressure, mouth_flow: self.observe(&native)? })
    }

    /// Apply one physical pressure/flow source law to each retained inlet
    /// mode and observe every mouth mode. The same native-wave solution
    /// supplies inlet work and mouth motion, including reactive source
    /// impedance, cross-mode reflection and coherent source interference.
    ///
    /// See [`ModalReflectionResponse::drive_sources`] for the source law,
    /// zero-flow/zero-pressure distinction and disconnected-resonance rule.
    ///
    /// # Errors
    /// The original source-boundary refusals or an invalid mouth observation.
    pub fn drive_sources(&self, sources: &[crate::network::PortExcitation])
        -> Result<MouthSourceResponse, DuctError>
    {
        let (inlet, native) = self.input.source_excitation(sources)?;
        Ok(MouthSourceResponse { inlet, mouth_flow: self.observe(&native)? })
    }

    fn observe(&self, native: &[C64]) -> Result<Vec<C64>, DuctError> {
        // The shared input solve already admitted n in 1..=MAX_MODES.
        let n = native.len();
        if self.native_flow_map.len() != n*n || self.native_flow_map.iter().any(|&v| !finite(v)) {
            return Err(bad("native mouth observation needs a finite complete modal map"));
        }
        let flow: Vec<_> = self.native_flow_map.chunks_exact(n).map(|row|
            row.iter().zip(native).fold(C64::ZERO, |sum, (&h, &a)| sum+h*a)).collect();
        if flow.iter().any(|&v| !finite(v)) { return Err(DuctError::Singular); }
        Ok(flow)
    }
}

pub(super) fn wave_response(response: super::ModalWaveResponse) -> Result<MouthWaveResponse, DuctError> {
    Ok(MouthWaveResponse {
        input: response.input,
        inlet_radius_m: response.inlet_radius_m,
        mouth_radius_m: response.mouth_radius_m,
        native_flow_map: response.mouth_flow.ok_or_else(|| bad("native mouth map was not requested"))?,
    })
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

/// Complete physical mouth matrix with every original signed modal flow.
/// The input impedance and mouth map use the same coupled boundary. The load
/// must match the exact frequency, retained modes and actual station radius.
pub fn mouth_response_matrix(duct: &Duct, gas: &GasState, omega: f64, loss: LossModel,
    load: &ModalMouthImpedance, modes: usize, extra_slices: usize)
    -> Result<MouthResponse, DuctError>
{
    super::mm_evaluate(duct, gas, omega, loss, MmLoad::Matrix(load), modes, extra_slices, true)
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

/// Exact finite-section transfer at a transverse cutoff. Its local
/// reference is Zc_0, so P_left=P_right-i k_0 Zc_0 L U_right and U_left=U_right
/// give r=b/(2+b), t=2/(2+b), b=-i k_0 L. Other modes retain their native
/// decaying propagation. The finite reference stays private; endpoint
/// observations recover physical pressure/flow, never native cutoff waves.
pub(super) fn propagate_cutoff(
    map: &mut Option<Vec<C64>>, reflection: &mut Vec<C64>, wave: &ModalWave, length: f64,
) -> Result<(), DuctError> {
    let n = wave.k.len();
    let b = C64::new(0.0, -1.0) * wave.k[0].scale(length);
    let denominator = C64::from_re(2.0) + b;
    let mut r = vec![C64::ZERO; n];
    let mut t = Vec::with_capacity(n);
    for (mode, &k) in wave.k.iter().enumerate() {
        t.push(if k == C64::ZERO {
            r[mode] = b / denominator;
            C64::from_re(2.0) / denominator
        } else {
            crate::exp_i(k.scale(length))
        });
    }
    if r.iter().chain(&t).any(|&z| !finite(z)) { return Err(DuctError::Singular); }
    // a_right = J a_left, J=(I-diag(r)R_right)^-1 diag(t).
    // The same feedback solve advances reflection and every mouth observer.
    let mut feedback = reflection.clone();
    for row in 0..n {
        for col in 0..n { feedback[row*n+col] = -r[row] * feedback[row*n+col]; }
        feedback[row*n+row] = feedback[row*n+row] + C64::ONE;
    }
    let mut j = mat_inv(&feedback, n)?;
    for row in 0..n {
        for col in 0..n { j[row*n+col] = j[row*n+col] * t[col]; }
    }
    let mut upstream = mat_mul(reflection, &j, n);
    for row in 0..n {
        if wave.k[row] == C64::ZERO {
            // U_left=U_right, hence I-R_left=(I-R_right)J on this row.
            // Form the conserved flow directly: an exactly closed cutoff
            // channel stays exactly closed, including its genuine impedance
            // pole, without cancellation between r and t R J.
            for col in 0..n {
                let mut flow = C64::ZERO;
                for mode in 0..n {
                    let y = if mode == row { C64::ONE - reflection[row*n+mode] }
                        else { -reflection[row*n+mode] };
                    flow = flow + y * j[mode*n+col];
                }
                upstream[row*n+col] = if col == row { C64::ONE - flow } else { -flow };
            }
        } else {
            for col in 0..n { upstream[row*n+col] = t[row] * upstream[row*n+col]; }
            upstream[row*n+row] = upstream[row*n+row] + r[row];
        }
    }
    if upstream.iter().any(|&z| !finite(z)) { return Err(DuctError::Singular); }
    if let Some(map) = map {
        *map = mat_mul(map, &j, n);
        if map.iter().any(|&z| !finite(z)) { return Err(DuctError::Singular); }
    }
    *reflection = upstream;
    Ok(())
}

/// Reflect and transfer forward amplitudes across an original junction.
/// Expansion: p_u=F p_d, U_d=F^T U_u. Contraction: p_d=F p_u,
/// U_u=F^T U_d. Eliminate incoming waves using b_d=R_d a_d, but do
/// not invert F, I-R_d or I+R_d. Interior pressure/flow nodes therefore
/// add no pole, even when the impedance-only route is requested. The same
/// one LU yields both a_d=J a_u and the upstream reflected wave.
pub(super) fn junction(map: &mut Option<Vec<C64>>, r: &[C64], z_down: &[C64],
    z_up: &[C64], junction: &prepared::PreparedJunction, n: usize)
    -> Result<Vec<C64>, DuctError>
{
    let f = &junction.projection;
    let ft = &junction.transpose;
    let mut p = r.to_vec();
    for i in 0..n { p[i*n+i] = p[i*n+i] + C64::ONE; }
    let y = terminal_flow(r, z_down, n);
    let (a, b) = if junction.expansion {
        let mut weighted_f = f.to_vec();
        for row in 0..n {
            for col in 0..n { weighted_f[row*n+col] = weighted_f[row*n+col] / z_up[row]; }
        }
        let mut a = mat_mul(&mat_mul(ft, &weighted_f, n), &p, n);
        for (a, y) in a.iter_mut().zip(&y) { *a = *a + *y; }
        let mut b = ft.to_vec();
        for row in 0..n {
            for col in 0..n { b[row*n+col] = b[row*n+col].scale(2.0) / z_up[col]; }
        }
        (a, b)
    } else {
        let mut weighted_ft = ft.to_vec();
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
    let reflected = if junction.expansion {
        // p_u = F (I+R_d) J a_u = (I+R_u) a_u.
        let mut value = mat_mul(&mat_mul(f, &p, n), &j, n);
        for k in 0..n { value[k*n+k] = value[k*n+k] - C64::ONE; }
        value
    } else {
        // U_u = F^T Z_d^-1 (I-R_d) J a_u = Z_u^-1 (I-R_u) a_u.
        let mut value = mat_mul(&mat_mul(ft, &y, n), &j, n);
        for row in 0..n {
            for col in 0..n { value[row*n+col] = -z_up[row] * value[row*n+col]; }
            value[row*n+row] = value[row*n+row] + C64::ONE;
        }
        value
    };
    if reflected.iter().any(|&v| !finite(v)) {
        return Err(DuctError::Singular);
    }
    if let Some(map) = map {
        *map = mat_mul(map, &j, n);
        if map.iter().any(|&v| !finite(v)) { return Err(DuctError::Singular); }
    }
    Ok(reflected)
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
