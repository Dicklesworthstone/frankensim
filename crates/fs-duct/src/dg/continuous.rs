//! Continuous access to the original DG spatial operator and its forced response.
//!
//! All gas, boundary, viscothermal and compliant-wall coordinates are retained.
//! Peak phasors use exp(-i omega t). This is the frozen linear spatial model;
//! its spatial and constitutive-band errors remain separate from solve residuals.
use super::*;
/// Transactional source-bound adoption after a complete coupled trajectory.
pub mod state;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;

/// Borrowed physical owner: x' = L x + sum_j B_j U_j, p_j = B_j^T x,
/// with storage x^T x / 2. Every U is volume flow INTO this duct.
/// Keeping L and B permits joint retained-state assembly even where an isolated
/// duct impedance has a pole. No state, basis, material or boundary is reduced.
#[derive(Clone, Copy)]
pub struct DgContinuous<'a> {
    duct: &'a DgDuct,
}

/// Original declared port, in boundary-then-side-port order.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicPort {
    pub id: PortId,
    /// Peak flow INTO the duct [m3/s]; an omitted input is exactly zero.
    pub inflow_m3_s: C64,
    /// Peak conjugate pressure from the original loading column [Pa].
    pub pressure_pa: C64,
    pub mean_input_power_w: f64,
}

/// AC particular solution of the original frozen DG operator.
#[derive(Clone, Debug, PartialEq)]
pub struct DgHarmonicResponse {
    pub omega_rad_s: f64,
    /// Complete original energy-coordinate ordering, including every wall arm.
    pub energy_coordinates: Vec<C64>,
    /// Physical pressure at original DG nodes [Pa].
    pub pressure_pa: Vec<C64>,
    /// Physical axial volume flow at original DG nodes [m3/s]. This is not the
    /// imposed boundary numerical flux; at the outlet an inflow has opposite sign.
    pub volume_flow_m3_s: Vec<C64>,
    pub ports: Vec<DgHarmonicPort>,
    /// Mean AC storage, including gas and all wall coordinates [J].
    pub mean_stored_energy_j: f64,
    /// Included auxiliary storage [J], also containing molecular gas memory
    /// when attached. The legacy field name predates internal gas reservoirs.
    pub mean_wall_energy_j: f64,
    /// Original positive drag/relaxation squares, including compliant walls [W].
    pub wall_dissipation_w: f64,
    /// Internal molecular energy exchange, distinct from thermal wall loss [W].
    pub molecular_dissipation_w: f64,
    /// Original matched characteristic boundaries' mean absorption [W].
    pub matched_absorption_w: f64,
    pub input_power_w: f64,
    /// Independent wall, molecular and boundary loss minus input power [W].
    pub power_residual_w: f64,
    /// Infinity-norm backward residual of (-i omega I - L)x = BU.
    /// This is not a forward-error or conditioning certificate near a pole.
    pub solve_relative_residual: f64,
}

fn invalid(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps: 0 })
    } else {
        Ok(())
    }
}
fn lu_error(error: LuComplexError) -> DuctError {
    match error {
        LuComplexError::Cancelled => DuctError::Cancelled { completed_steps: 0 },
        LuComplexError::Singular { .. } => DuctError::Singular,
    }
}

impl DgDuct {
    /// Immutable, allocation-free view of the exact source-bound spatial model.
    #[must_use]
    pub const fn continuous_operator(&self) -> DgContinuous<'_> {
        DgContinuous { duct: self }
    }

    /// Solve continuous peak-phasor flow forcing with exp(-i omega t).
    /// Missing ports receive zero flow; unknown or repeated IDs refuse.
    /// All memories remain explicit, using the existing cancellable complex LU.
    /// A singular isolated lossless response refuses; coupled consumers retain
    /// these coordinates through continuous_operator instead of inverting it.
    ///
    /// Work is charged before allocation as n^2*(n+ports+8), under the existing
    /// complete-state DG cap. Cancellation is checked during assembly, every LU
    /// row, residual evaluation and observations. Initial transients, finite-step
    /// held-flow response, nonlinear flow and unmodeled external loads are absent.
    /// The admitted viscothermal band is enforced. No accepted state is touched.
    pub fn harmonic_response(
        &self,
        omega_rad_s: f64,
        inflows: &[(PortId, C64)],
        maximum_work: usize,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgHarmonicResponse, DuctError> {
        let view = self.continuous_operator();
        poll(cancelled)?;
        view.admit_frequency(omega_rad_s)?;
        view.admit_work(maximum_work, true)?;
        admit_tolerance(relative_tolerance)?;
        let flows = view.resolve_flows(inflows)?;
        let n = view.dimension();
        let scale = self.l_tilde.iter().fold(omega_rad_s, |a, b| a.max(b.abs()));
        if !scale.is_finite() {
            return Err(invalid("DG harmonic operator is not finite"));
        }
        let mut matrix = vec![C64::ZERO; n * n];
        let mut state = vec![C64::ZERO; n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                matrix[i * n + j] = C64::new(
                    -self.l_tilde[i * n + j] / scale,
                    if i == j { -omega_rad_s / scale } else { 0.0 },
                );
            }
            let rhs = view.forcing(i, &flows);
            state[i] = C64::new(rhs.re / scale, rhs.im / scale);
        }
        if state.iter().any(|&x| !finite(x)) {
            return Err(invalid("DG harmonic forcing overflow"));
        }
        let factor = lu_complex_with_cancel(&matrix, n, &mut *cancelled).map_err(lu_error)?;
        factor
            .solve_with_cancel(&mut state, &mut *cancelled)
            .map_err(lu_error)?;
        view.observe_harmonic(
            omega_rad_s,
            &state,
            inflows,
            maximum_work,
            relative_tolerance,
            cancelled,
        )
    }
}

fn admit_tolerance(tolerance: f64) -> Result<(), DuctError> {
    if !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0 {
        Err(invalid(
            "DG harmonic backward-residual tolerance must be in (0,1)",
        ))
    } else {
        Ok(())
    }
}

impl DgContinuous<'_> {
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.duct.total_dofs()
    }
    /// Original row-major n by n generator [1/s], borrowed without reassembly.
    #[must_use]
    pub fn generator_per_s(&self) -> &[f64] {
        &self.duct.l_tilde
    }
    /// Original inflow columns; their transpose observes conjugate pressure.
    #[must_use]
    pub fn ports(&self) -> &[(PortId, Vec<f64>)] {
        &self.duct.ports
    }
    #[must_use]
    pub fn operator_identity(&self) -> ContentHash {
        self.duct.operator_identity
    }

    /// Positive finite frequency inside any declared viscothermal validity band.
    pub fn admit_frequency(&self, omega_rad_s: f64) -> Result<(), DuctError> {
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("DG harmonic frequency must be positive and finite"));
        }
        if let Some([lo, hi]) = self.duct.viscothermal_band_rad_s()
            && (omega_rad_s < lo || omega_rad_s > hi)
        {
            return Err(invalid(
                "DG harmonic frequency exceeds the retained wall-law band",
            ));
        }
        Ok(())
    }

    fn admit_work(&self, maximum: usize, factor: bool) -> Result<(), DuctError> {
        let n = self.dimension();
        let terms = self.ports().len() + 8 + if factor { n } else { 0 };
        let work = n.checked_mul(n).and_then(|v| v.checked_mul(terms));
        if n == 0 || n > MAX_DG_PROPAGATOR_DOFS || work.is_none_or(|v| v > maximum) {
            Err(invalid(
                "complete DG harmonic operation exceeds its explicit work budget",
            ))
        } else {
            Ok(())
        }
    }
    fn resolve_flows(&self, inputs: &[(PortId, C64)]) -> Result<Vec<C64>, DuctError> {
        if inputs.len() > self.ports().len() {
            return Err(invalid("too many DG harmonic port inputs"));
        }
        let mut flows = vec![C64::ZERO; self.ports().len()];
        for (i, &(id, flow)) in inputs.iter().enumerate() {
            if !finite(flow) || inputs[..i].iter().any(|(old, _)| *old == id) {
                return Err(invalid(
                    "DG harmonic inputs need finite flows and distinct port IDs",
                ));
            }
            let index = self
                .ports()
                .iter()
                .position(|(port, _)| *port == id)
                .ok_or_else(|| invalid("DG harmonic input names an undeclared port"))?;
            flows[index] = flow;
        }
        Ok(flows)
    }
    fn forcing(&self, row: usize, flows: &[C64]) -> C64 {
        self.ports()
            .iter()
            .zip(flows)
            .fold(C64::ZERO, |sum, ((_, column), &flow)| {
                sum + flow.scale(column[row])
            })
    }

    /// Observe a supplied COMPLETE original energy phasor, including one solved
    /// as part of a larger coupled system. This does not invert the isolated
    /// duct, so valid joint responses at a bare duct pole remain representable.
    /// Original equations are checked using the actual supplied port inflows.
    /// Work is charged as n^2*(ports+8), before allocating physical observations.
    /// The field map and positive wall-loss owner are reused for real/imaginary
    /// quadratures; matched-boundary loss uses its actual endpoint gas and area.
    pub fn observe_harmonic(
        &self,
        omega_rad_s: f64,
        energy_coordinates: &[C64],
        inflows: &[(PortId, C64)],
        maximum_work: usize,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgHarmonicResponse, DuctError> {
        poll(cancelled)?;
        self.admit_frequency(omega_rad_s)?;
        self.admit_work(maximum_work, false)?;
        admit_tolerance(relative_tolerance)?;
        let n = self.dimension();
        if energy_coordinates.len() != n || energy_coordinates.iter().any(|&z| !finite(z)) {
            return Err(invalid(
                "DG harmonic observation needs every finite original coordinate",
            ));
        }
        let flows = self.resolve_flows(inflows)?;
        let s = C64::new(0.0, -omega_rad_s);
        let norm_x = energy_coordinates
            .iter()
            .map(|z| z.abs())
            .fold(0.0_f64, f64::max);
        let (mut norm_a, mut norm_b, mut norm_r) = (0.0_f64, 0.0_f64, 0.0_f64);
        for i in 0..n {
            poll(cancelled)?;
            let rhs = self.forcing(i, &flows);
            let mut residual = -rhs;
            let mut row_norm = 0.0;
            for (j, &x) in energy_coordinates.iter().enumerate() {
                let a = C64::from_re(-self.duct.l_tilde[i * n + j])
                    + if i == j { s } else { C64::ZERO };
                residual = residual + a * x;
                row_norm += a.abs();
            }
            if !finite(residual) || !finite(rhs) || !row_norm.is_finite() {
                return Err(invalid("DG harmonic original-equation residual overflow"));
            }
            norm_a = norm_a.max(row_norm);
            norm_b = norm_b.max(rhs.abs());
            norm_r = norm_r.max(residual.abs());
        }
        let denominator = norm_a * norm_x + norm_b;
        if !denominator.is_finite() {
            return Err(invalid("DG harmonic residual scale overflow"));
        }
        let residual = if denominator == 0.0 {
            0.0
        } else {
            norm_r / denominator
        };
        if residual > relative_tolerance {
            return Err(invalid(
                "DG harmonic original-equation residual exceeds tolerance",
            ));
        }
        let mut native = self.duct.zero_state();
        for (x, z) in native.tilde.iter_mut().zip(energy_coordinates) {
            *x = z.re;
        }
        let (real_p, real_q) = self.duct.from_energy_coords(&native.tilde);
        let real_loss = self.duct.wall_dissipated_power_w(&native)?;
        let real_molecular = self.duct.molecular_dissipated_power_w(&native)?;
        poll(cancelled)?;
        for (x, z) in native.tilde.iter_mut().zip(energy_coordinates) {
            *x = z.im;
        }
        let (imag_p, imag_q) = self.duct.from_energy_coords(&native.tilde);
        let wall_dissipation_w = 0.5 * (real_loss + self.duct.wall_dissipated_power_w(&native)?);
        let molecular_dissipation_w = 0.5 * (real_molecular
            + self.duct.molecular_dissipated_power_w(&native)?);
        let pressure_pa: Vec<_> = real_p
            .into_iter()
            .zip(imag_p)
            .map(|(r, i)| C64::new(r, i))
            .collect();
        let volume_flow_m3_s: Vec<_> = real_q
            .into_iter()
            .zip(imag_q)
            .map(|(r, i)| C64::new(r, i))
            .collect();
        let mut matched_absorption_w = 0.0;
        for (boundary, node, gas, area) in [
            (
                self.duct.inlet,
                0,
                &self.duct.gas,
                self.duct.inlet_area_m2(),
            ),
            (
                self.duct.outlet,
                self.duct.field_dofs - 1,
                &self.duct.outlet_gas,
                self.duct.outlet_area_m2(),
            ),
        ] {
            if boundary == DgBoundary::Matched {
                let zc = gas.density * gas.sound_speed / area;
                matched_absorption_w += 0.25
                    * (pressure_pa[node].norm_sq() / zc + zc * volume_flow_m3_s[node].norm_sq());
            }
        }
        let mut ports = Vec::with_capacity(self.ports().len());
        let mut input_power_w = 0.0;
        for ((id, column), &flow) in self.ports().iter().zip(&flows) {
            poll(cancelled)?;
            let pressure = column
                .iter()
                .zip(energy_coordinates)
                .fold(C64::ZERO, |sum, (&b, &x)| sum + x.scale(b));
            let power = 0.5 * (pressure * flow.conj()).re;
            input_power_w += power;
            ports.push(DgHarmonicPort {
                id: *id,
                inflow_m3_s: flow,
                pressure_pa: pressure,
                mean_input_power_w: power,
            });
        }
        let storage = 0.25 * energy_coordinates.iter().map(|x| x.norm_sq()).sum::<f64>();
        let wall_storage = 0.25
            * energy_coordinates[2 * self.duct.field_dofs..]
                .iter()
                .map(|x| x.norm_sq())
                .sum::<f64>();
        let power_residual_w = wall_dissipation_w + molecular_dissipation_w
            + matched_absorption_w - input_power_w;
        if pressure_pa
            .iter()
            .chain(&volume_flow_m3_s)
            .any(|&z| !finite(z))
            || ports
                .iter()
                .any(|p| !finite(p.pressure_pa) || !p.mean_input_power_w.is_finite())
            || [
                storage,
                wall_storage,
                wall_dissipation_w,
                molecular_dissipation_w,
                matched_absorption_w,
                input_power_w,
                power_residual_w,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "DG harmonic physical field, energy or power overflow",
            ));
        }
        poll(cancelled)?;
        Ok(DgHarmonicResponse {
            omega_rad_s,
            energy_coordinates: energy_coordinates.to_vec(),
            pressure_pa,
            volume_flow_m3_s,
            ports,
            mean_stored_energy_j: storage,
            mean_wall_energy_j: wall_storage,
            wall_dissipation_w,
            molecular_dissipation_w,
            matched_absorption_w,
            input_power_w,
            power_residual_w,
            solve_relative_residual: residual,
        })
    }
}

#[cfg(test)]
mod oracle_tests;
#[cfg(test)]
mod tests;
