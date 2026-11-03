//! Forced periodic response of the same complete solid/material/load pencil.
use super::{PreparedViscoPoles, ViscoPoleError, cap, finite, invalid, poll};
use fs_exec::Cx;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::c64::C64;

/// Additional resources for one reusable factor and one complete response.
/// The original prepared source/load caches retain their own admitted budgets.
#[derive(Clone, Copy, Debug)]
pub struct ViscoJointHarmonicBudget {
    pub maximum_coordinates: usize,
    /// Numerical payload, including LU and physical output; excludes RSS.
    pub maximum_payload_bytes: usize,
    /// Conservative factor/solve/refinement and original observation work.
    pub maximum_work: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ViscoJointHarmonicError {
    Pencil(ViscoPoleError),
    /// A true joint singularity or cooperative LU cancellation, never a bare
    /// component impedance failure.
    Linear(LuComplexError),
    PowerBalance {
        relative: f64,
        tolerance: f64,
    },
}
impl From<ViscoPoleError> for ViscoJointHarmonicError {
    fn from(value: ViscoPoleError) -> Self {
        Self::Pencil(value)
    }
}
impl From<LuComplexError> for ViscoJointHarmonicError {
    fn from(value: LuComplexError) -> Self {
        Self::Linear(value)
    }
}
impl core::fmt::Display for ViscoJointHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "joint hereditary harmonic response refused: {self:?}")
    }
}
impl std::error::Error for ViscoJointHarmonicError {}

/// Original peak load phasors and independently observed cycle means.
#[derive(Clone, Debug)]
pub struct ViscoJointHarmonicLoad {
    /// Every original differential load state and algebraic reaction; empty for
    /// an attached K/C/M law. Prepared source ranges identify the two blocks.
    pub state: Vec<C64>,
    pub displacement: Vec<C64>,
    pub flow: Vec<C64>,
    pub effort: Vec<C64>,
    pub average_storage_j: f64,
    pub average_dissipation_w: f64,
    /// Re(flow^H effort)/2, positive into the load; kept distinct from heat.
    pub average_absorbed_power_w: f64,
    /// Absorbed power minus the load's independently calculated dissipation.
    pub power_balance_residual_w: f64,
}

/// Original linear boundary amplitudes and cycle-mean signed port exchange.
/// Storage and heat cannot be inferred from the rational coefficients alone.
#[derive(Clone, Debug)]
pub struct ViscoJointHarmonicStateSpaceLoad {
    pub state: Vec<C64>,
    pub displacement: Vec<C64>,
    pub flow: Vec<C64>,
    pub effort: Vec<C64>,
    /// Re(conj(flow_p)*effort_p)/2; negative mutual exchanges remain signed.
    pub average_absorbed_port_power_w: Vec<f64>,
    pub average_absorbed_power_w: f64,
}

/// Peak exp(-i*omega*t) state and cycle-averaged physical accounts.
#[derive(Clone, Debug)]
pub struct ViscoJointHarmonicResponse {
    pub omega_rad_s: f64,
    /// Free displacement, every original Mandel Maxwell strain, then every
    /// original dynamic load coordinate. No homogeneous state is added.
    pub physical_vector: Vec<C64>,
    pub displacement_m: Vec<C64>,
    pub velocity_m_s: Vec<C64>,
    /// Separate original coefficient-scaled residuals of the three families.
    pub mechanical_residual: f64,
    pub memory_residual: f64,
    pub dynamic_load_residual: f64,
    /// Original equilibrium, consistent inertia and Maxwell strain storage.
    pub average_body_storage_j: f64,
    /// Complete physical body+load mean storage, when every boundary supplies
    /// a physical storage witness. None for a generic state-space boundary.
    pub average_total_storage_j: Option<f64>,
    /// Original Kelvin viscosity and every Maxwell branch, counted once.
    pub average_material_heat_w: f64,
    pub average_input_power_w: f64,
    pub port_load: Option<ViscoJointHarmonicLoad>,
    pub dynamic_load: Option<ViscoJointHarmonicLoad>,
    pub state_space_load: Option<ViscoJointHarmonicStateSpaceLoad>,
    /// Applied power - material heat - actual power into the attached loads.
    pub body_power_balance_residual_w: f64,
    /// Applied power - material heat - each load's independent dissipation.
    /// None when a boundary has no physical storage/dissipation witness;
    /// its actual exchange is still included in the body balance above.
    pub total_power_balance_residual_w: Option<f64>,
}

/// One LU at one physical frequency, reusable across forcing columns/phases.
pub struct ViscoJointHarmonicFactor<'m, 's> {
    model: &'m PreparedViscoPoles<'s>,
    omega_rad_s: f64,
    lu: LuComplex,
}

impl<'s> PreparedViscoPoles<'s> {
    /// Factor the complete original hereditary polynomial at s=-i*omega.
    /// The same assembled coefficients and fixed conditioning used for natural
    /// modes retain all independent load and material states. A bare load or
    /// body resonance is regular whenever the complete forced problem is.
    ///
    /// Positive frequency and the common original source/load band are required.
    /// Preparation is O(N^3); repeated right-hand sides reuse this LU. No time
    /// step, accepted history, condensed impedance or new material fit is used.
    pub fn harmonic_factor(
        &self,
        omega_rad_s: f64,
        budget: ViscoJointHarmonicBudget,
        cx: &Cx<'_>,
    ) -> Result<ViscoJointHarmonicFactor<'_, 's>, ViscoJointHarmonicError> {
        let mut cancelled = || cx.checkpoint().is_err();
        poll(&mut cancelled)?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid("joint harmonic response needs positive finite frequency").into());
        }
        let n = self.dimension;
        let d = n as u128;
        let p = (self.port_load.as_ref().map_or(0, |l| l.port_count())
            + self.dynamic_load.as_ref().map_or(0, |l| l.port_count())
            + self.state_space_load.as_ref().map_or(0, |l| l.port_count())) as u128;
        let quadrature_work: u128 = self
            .source
            .elements
            .iter()
            .map(|e| e.points.len() as u128 * (e.dofs.len() as u128 * 128 + 512))
            .sum();
        cap(
            "joint harmonic coordinates",
            d,
            budget.maximum_coordinates as u128,
        )?;
        cap(
            "joint harmonic numerical payload bytes",
            96 * d * d + 512 * d + 256 * p,
            budget.maximum_payload_bytes as u128,
        )?;
        cap(
            "joint harmonic work",
            16 * d * d * d + 512 * (d * d + d * p + p * p) + quadrature_work,
            budget.maximum_work as u128,
        )?;
        let s = C64::new(0.0, -omega_rad_s);
        self.argument(s, &mut 0, &mut cancelled)?;
        let mut matrix = vec![C64::ZERO; n * n];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                let at = i * n + j;
                matrix[at] = (C64::from_re(self.coefficients[0][at])
                    + s.scale(self.coefficients[1][at])
                    + (s * s).scale(self.coefficients[2][at]))
                .scale(self.column_scales[j] / self.row_scales[i]);
            }
        }
        if matrix.iter().any(|&v| !finite(v)) {
            return Err(invalid("joint harmonic pencil is not representable").into());
        }
        let lu = lu_complex_with_cancel(&matrix, n, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(ViscoJointHarmonicFactor {
            model: self,
            omega_rad_s,
            lu,
        })
    }
}

impl ViscoJointHarmonicFactor<'_, '_> {
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub fn free_dofs(&self) -> &[usize] {
        self.model.free_dofs()
    }
    pub fn dimension(&self) -> usize {
        self.model.dimension
    }

    /// Apply one peak force [N] per original free displacement coordinate.
    /// Every quadrature point passes the original full-cycle material strain
    /// admission before a response is published. Sources and this LU remain
    /// unchanged after a source, numerical, resource or cancellation refusal.
    ///
    /// Original mechanical, Maxwell and independent load rows are checked
    /// separately. Up to three corrective solves reuse this same LU; no
    /// residual is turned into artificial damping or assigned to a loss channel.
    pub fn solve(
        &self,
        force_n: &[C64],
        cx: &Cx<'_>,
    ) -> Result<ViscoJointHarmonicResponse, ViscoJointHarmonicError> {
        let model = self.model;
        let (n, m) = (model.dimension, model.free_dofs().len());
        let mut cancelled = || cx.checkpoint().is_err();
        poll(&mut cancelled)?;
        if force_n.len() != m || force_n.iter().any(|&v| !finite(v)) {
            return Err(
                invalid("one finite peak force per original free coordinate required").into(),
            );
        }
        let s = C64::new(0.0, -self.omega_rad_s);
        let mut state = vec![C64::ZERO; n];
        for i in 0..m {
            state[i] = divide(force_n[i], model.row_scales[i]);
        }
        self.lu.solve_with_cancel(&mut state, &mut cancelled)?;
        for (i, x) in state.iter_mut().enumerate() {
            *x = x.scale(model.column_scales[i]);
        }
        let mut evaluations = 0;
        let mut accepted = None;
        for attempt in 0..4 {
            model.argument(s, &mut evaluations, &mut cancelled)?;
            let (residual, mut correction) =
                self.original_residual(s, &state, force_n, &mut cancelled)?;
            if residual
                .iter()
                .all(|&r| r <= model.budget.physical_residual_tolerance)
            {
                accepted = Some(residual);
                break;
            }
            if attempt == 3 {
                break;
            }
            self.lu.solve_with_cancel(&mut correction, &mut cancelled)?;
            for i in 0..n {
                state[i] = state[i] + correction[i].scale(model.column_scales[i]);
            }
        }
        let residual = if let Some(r) = accepted {
            r
        } else {
            let (r, _) = self.original_residual(s, &state, force_n, &mut cancelled)?;
            return Err(ViscoPoleError::PhysicalResidual {
                mechanical: r[0],
                memory: r[1],
                dynamic_load: r[2],
                tolerance: model.budget.physical_residual_tolerance,
            }
            .into());
        };
        model
            .source
            .validate_harmonic_displacement(&state[..m], cx)
            .map_err(ViscoPoleError::Source)?;
        let port_load = model
            .port_load
            .as_ref()
            .map(|load| {
                let p = load.observe(s, &state[..m], &mut cancelled)?;
                Ok::<_, ViscoPoleError>(ViscoJointHarmonicLoad {
                    state: Vec::new(),
                    displacement: p.displacement,
                    flow: p.flow,
                    effort: p.effort,
                    average_storage_j: 0.5 * p.storage_j,
                    average_dissipation_w: 0.5 * p.dissipation_w,
                    average_absorbed_power_w: 0.5 * p.absorbed_power_w,
                    power_balance_residual_w: 0.5 * p.signed_power_defect_w,
                })
            })
            .transpose()?;
        let mut body = None;
        let dynamic_load = model
            .dynamic_load
            .as_ref()
            .map(|load| {
                let p = load.observe(model, s, &state, &mut cancelled)?;
                body = Some((p.body_storage_j, p.body_dissipation_w));
                Ok::<_, ViscoPoleError>(ViscoJointHarmonicLoad {
                    state: p.state,
                    displacement: p.displacement,
                    flow: p.flow,
                    effort: p.effort,
                    average_storage_j: 0.5 * p.storage_j,
                    average_dissipation_w: 0.5 * p.dissipation_w,
                    average_absorbed_power_w: 0.5 * p.absorbed_power_w,
                    power_balance_residual_w: 0.5 * p.signed_power_defect_w,
                })
            })
            .transpose()?;
        let state_space_load = model
            .state_space_load
            .as_ref()
            .map(|load| {
                let p = load.observe(model, s, &state, &mut cancelled)?;
                body = Some((p.body_storage_j, p.body_dissipation_w));
                Ok::<_, ViscoPoleError>(ViscoJointHarmonicStateSpaceLoad {
                    state: p.state,
                    displacement: p.displacement,
                    flow: p.flow,
                    effort: p.effort,
                    average_absorbed_port_power_w: p
                        .absorbed_port_power_w
                        .into_iter()
                        .map(|x| 0.5 * x)
                        .collect(),
                    average_absorbed_power_w: 0.5 * p.absorbed_power_w,
                })
            })
            .transpose()?;
        let (body_storage, body_heat) = match body {
            Some(x) => x,
            None => super::dynamic::body_energy(model, s, &state, &mut cancelled)?,
        };
        let displacement_m = state[..m].to_vec();
        let velocity_m_s: Vec<_> = displacement_m.iter().map(|&u| s * u).collect();
        let input = 0.5
            * force_n
                .iter()
                .zip(&velocity_m_s)
                .map(|(f, v)| (f.conj() * *v).re)
                .sum::<f64>();
        let mut scale = 0.5
            * force_n
                .iter()
                .zip(&velocity_m_s)
                .map(|(f, v)| f.abs() * v.abs())
                .sum::<f64>();
        let (mut absorbed, mut load_heat, mut load_storage) = (0.0, 0.0, 0.0);
        let mut load_error = 0.0_f64;
        for load in port_load.iter().chain(&dynamic_load) {
            poll(&mut cancelled)?;
            absorbed += load.average_absorbed_power_w;
            load_heat += load.average_dissipation_w;
            load_storage += load.average_storage_j;
            let apparent = 0.5
                * load
                    .flow
                    .iter()
                    .zip(&load.effort)
                    .map(|(u, e)| u.abs() * e.abs())
                    .sum::<f64>();
            // At an antiresonance the interface flow and its real power can
            // vanish while an independent lossless load still oscillates.
            // Its original oscillatory storage sets a physical work scale;
            // normalizing only by roundoff-sized interface power would reject
            // that valid solution. No storage term is counted as mean heat.
            let denom = apparent
                + load.average_dissipation_w.abs()
                + 2.0 * self.omega_rad_s * load.average_storage_j.abs();
            load_error = load_error.max(ratio(load.power_balance_residual_w.abs(), denom));
            scale += denom;
        }
        if let Some(load) = &state_space_load {
            absorbed += load.average_absorbed_power_w;
            scale += 0.5
                * load
                    .flow
                    .iter()
                    .zip(&load.effort)
                    .map(|(v, e)| v.abs() * e.abs())
                    .sum::<f64>();
        }
        let average_body_storage_j = 0.5 * body_storage;
        let average_total_storage_j = state_space_load
            .is_none()
            .then_some(average_body_storage_j + load_storage);
        let average_material_heat_w = 0.5 * body_heat;
        let body_defect = input - average_material_heat_w - absorbed;
        let total_defect = state_space_load
            .is_none()
            .then_some(input - average_material_heat_w - load_heat);
        scale +=
            average_material_heat_w.abs() + 2.0 * self.omega_rad_s * average_body_storage_j.abs();
        let relative = load_error
            .max(ratio(body_defect.abs(), scale))
            .max(total_defect.map_or(0.0, |d| ratio(d.abs(), scale)));
        if [
            input,
            average_body_storage_j,
            average_total_storage_j.unwrap_or(0.0),
            average_material_heat_w,
            body_defect,
            total_defect.unwrap_or(0.0),
            scale,
            relative,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("joint harmonic physical account is not representable").into());
        }
        if relative > model.budget.physical_residual_tolerance {
            return Err(ViscoJointHarmonicError::PowerBalance {
                relative,
                tolerance: model.budget.physical_residual_tolerance,
            });
        }
        poll(&mut cancelled)?;
        Ok(ViscoJointHarmonicResponse {
            omega_rad_s: self.omega_rad_s,
            physical_vector: state,
            displacement_m,
            velocity_m_s,
            mechanical_residual: residual[0],
            memory_residual: residual[1],
            dynamic_load_residual: residual[2],
            average_body_storage_j,
            average_total_storage_j,
            average_material_heat_w,
            average_input_power_w: input,
            port_load,
            dynamic_load,
            state_space_load,
            body_power_balance_residual_w: body_defect,
            total_power_balance_residual_w: total_defect,
        })
    }

    fn original_residual(
        &self,
        s: C64,
        state: &[C64],
        force: &[C64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<([f64; 3], Vec<C64>), ViscoPoleError> {
        let model = self.model;
        let (n, m) = (model.dimension, model.free_dofs().len());
        if state.iter().any(|&v| !finite(v)) {
            return Err(invalid("joint harmonic state overflow"));
        }
        let mut numerators = [0.0_f64; 3];
        let mut denominators = [0.0_f64; 3];
        let algebraic_start = model.individually_checked_load_rows();
        let mut algebraic_residual = 0.0_f64;
        let mut correction = vec![C64::ZERO; n];
        for i in 0..n {
            poll(cancelled)?;
            let rhs = if i < m { force[i] } else { C64::ZERO };
            let (mut sum, mut compensation) = (C64::ZERO, C64::ZERO);
            let mut envelope = rhs.abs();
            for (j, &x) in state.iter().enumerate() {
                let at = i * n + j;
                let a = C64::from_re(model.coefficients[0][at])
                    + s.scale(model.coefficients[1][at])
                    + (s * s).scale(model.coefficients[2][at]);
                let term = a * x - compensation;
                let next = sum + term;
                compensation = (next - sum) - term;
                sum = next;
                envelope += model.envelope(at, s.abs()) * x.abs();
            }
            let defect = rhs - sum;
            correction[i] = divide(defect, model.row_scales[i]);
            if !finite(defect) || !finite(correction[i]) || !envelope.is_finite() {
                return Err(invalid("joint harmonic original equation overflow"));
            }
            let group = if i < m {
                0
            } else if i < model.body_dimension() {
                1
            } else {
                2
            };
            numerators[group] = numerators[group].max(defect.abs());
            denominators[group] = denominators[group].max(envelope);
            if i >= algebraic_start {
                algebraic_residual = algebraic_residual.max(ratio(defect.abs(), envelope));
            }
        }
        let mut residual: [f64; 3] =
            core::array::from_fn(|i| ratio(numerators[i], denominators[i]));
        // Include every original constraint in the same-factor refinement
        // gate, even beside much larger differential rows in different units.
        residual[2] = residual[2].max(algebraic_residual);
        Ok((residual, correction))
    }
}

fn divide(z: C64, scale: f64) -> C64 {
    C64::new(z.re / scale, z.im / scale)
}
fn ratio(numerator: f64, denominator: f64) -> f64 {
    if denominator == 0.0 {
        if numerator == 0.0 { 0.0 } else { f64::INFINITY }
    } else {
        numerator / denominator
    }
}
