//! A retained translation/twist loading spring around the native shear owner.
//!
//! Coordinates are q=[u_x,u_y,l*theta], conjugate to w=[F_x,F_y,M/l].
//! The spring's midpoint equation is
//! w=K*(e_old+(d_stage-d_q)/2). Eliminating it gives
//! d_q=d_stage+2*e_old-2*K^-1*w. With fixed cell kinematics d_u=A*d_q,
//! the existing contact response therefore receives the exact reciprocal
//! addition 2*A*K^-1*A^T. All nonlinear layer trials still start from the
//! same original checkpoint and consume one physical interval.
//!
//! This is an infinitesimal-rotation, fixed-registration loading fixture.
//! The caller's original body response still owns body motion and its work;
//! neither that response nor the loading spring is an extra shear material.

use std::sync::Arc;

use fs_exec::Cx;
use fs_la::factor::{FactorRunError, cholesky_cancellable};
use fs_tribo::TangentialSolveConfig;
use fs_tribo::layer_contact::ThirdBodyLayerDissipation;
use fs_tribo::tangential::multiple::{
    MultipleTangentialError, MultipleTangentialStep, PreparedTangentialContacts,
};

use super::{
    FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt, FiniteGapPlasticLayerContact,
    FiniteGapPlasticLayerError, FiniteGapPlasticLayerNormalAdvance, FiniteGapPlasticLayerState,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt, FiniteGapPlasticReceipt,
    FiniteGapPlasticState,
};

/// Fixed physical spring and explicit whole-fixture admission tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasticLayerShearDriveConfig {
    /// Symmetric positive-definite stiffness [N/m] in coordinates
    /// [translation_x, translation_y, spin_length*rotation].
    /// Cross terms couple translation and twist. No diagonalization is inferred.
    pub stiffness_n_per_m: [[f64; 3]; 3],
    /// Authored positive reference length [m]; not a fitted material parameter.
    /// Changing it requires the conjugate stiffness congruence.
    pub spin_length_m: f64,
    /// Inclusive absolute accumulated carriage-rotation bound [rad].
    /// Positive and at most 0.1 for this infinitesimal-rotation fixture.
    pub maximum_rotation_rad: f64,
    /// Absolute residual tolerance [N] in all three scaled effort coordinates.
    pub generalized_force_tolerance_n: f64,
    /// Absolute tolerance for the complete spring/contact/body work balance [J].
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
struct DriveOwner {
    contact: FiniteGapPlasticLayerContact,
    config: PlasticLayerShearDriveConfig,
    lower: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
    rows: Vec<[f64; 3]>,
    response: PreparedTangentialContacts,
}

/// Full passive loading spring, compiled once around one original contact owner.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticLayerShearDrive {
    owner: Arc<DriveOwner>,
}

/// Indivisible original normal/layer histories and retained loading spring.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticLayerShearDriveState {
    owner: Arc<DriveOwner>,
    contact: FiniteGapPlasticLayerState,
    extension_m: [f64; 3],
    carriage_m: [f64; 3],
}

impl PlasticLayerShearDriveState {
    /// Native contact checkpoint; normal indentation and every shear branch survive.
    pub fn contact_state(&self) -> &FiniteGapPlasticLayerState {
        &self.contact
    }

    /// Spring extension [m,m,rad] in the original fixed tangent frame.
    pub fn spring_extension_m_rad(&self) -> [f64; 3] {
        unscale(self.extension_m, self.owner.config.spin_length_m)
    }

    /// Accumulated carriage translation and infinitesimal twist [m,m,rad].
    pub fn carriage_position_m_rad(&self) -> [f64; 3] {
        unscale(self.carriage_m, self.owner.config.spin_length_m)
    }

    /// Original loading-spring potential [J], separate from layer storage.
    pub fn spring_energy_j(&self) -> f64 {
        let value = transpose_apply(self.owner.lower, self.extension_m);
        0.5 * dot(value, value)
    }
}

/// One simultaneous physical shear interval; no candidate publishes partially.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticLayerShearDriveAdvance {
    /// Adopt this complete state only after the enclosing body/thermal owners accept.
    pub candidate: PlasticLayerShearDriveState,
    /// Original full-matrix local results and independent contact work checks.
    pub shear: MultipleTangentialStep<()>,
    /// Original normal pressures [Pa], in fixed site order.
    pub pressure_pa: Vec<f64>,
    /// Original physical plastic/Maxwell/Kelvin losses [J], in fixed site order.
    pub dissipation: Vec<ThirdBodyLayerDissipation>,
    /// Actual carriage translation and infinitesimal twist increment [m,m,rad].
    pub carriage_increment_m_rad: [f64; 3],
    /// Mean spring effort entering the interface [N,N,N m].
    pub mean_wrench_n_n_m: [f64; 3],
    /// Midpoint spring-equilibrium residual [N,N,N] in scaled coordinates.
    pub spring_equilibrium_residual_n: [f64; 3],
    /// Actual mean effort times prescribed stage increment [J].
    pub stage_work_j: f64,
    /// Signed original spring-potential change [J], evaluated by polarization.
    pub spring_stored_change_j: f64,
    /// Signed work delivered to the caller's original body response [J].
    /// This is sum F*(C_body*F-du_body_free), not heat or invented body storage.
    pub body_response_work_j: f64,
    /// stage_work - spring_change - layer_change - layer_heat - body_work [J].
    pub balance_residual_j: f64,
}

/// A normal substep retains the carriage and loading spring without shear time.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticLayerShearDriveNormalAdvance<B = ()> {
    /// Complete retained fixture state.
    pub candidate: PlasticLayerShearDriveState,
    /// Original normal storage/plastic-work receipt.
    pub normal: FiniteGapPlasticReceipt,
    /// Original normal load/spring receipt, or unit for prescribed approach.
    pub boundary: B,
}

/// No refusal changes original contact or apparatus history.
#[derive(Clone, Debug, PartialEq)]
pub enum PlasticLayerShearDriveError {
    /// Invalid spring, shapes, arithmetic, or a foreign fixture checkpoint.
    Invalid(&'static str),
    /// Original normal, constitutive or simultaneous-contact failure.
    Contact(FiniteGapPlasticLayerError),
    /// The candidate leaves the explicitly bounded small-rotation domain.
    Rotation { value_rad: f64, maximum_rad: f64 },
    /// Original spring equilibrium does not close in scaled force coordinates.
    Equilibrium {
        residual_n: [f64; 3],
        tolerance_n: f64,
    },
    /// The complete physical work budget does not close.
    Energy { residual_j: f64, tolerance_j: f64 },
}

impl From<FiniteGapPlasticLayerError> for PlasticLayerShearDriveError {
    fn from(value: FiniteGapPlasticLayerError) -> Self {
        Self::Contact(value)
    }
}
impl core::fmt::Display for PlasticLayerShearDriveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "plastic layer loading spring: {self:?}")
    }
}
impl std::error::Error for PlasticLayerShearDriveError {}

impl PlasticLayerShearDrive {
    /// Compile the exact spring elimination into the original full response.
    ///
    /// The original matrix and this augmented matrix each retain the existing
    /// contact/term bounds. No iteration factors, contact laws or histories are
    /// replaced. Preparation is O(n^2); every physical step uses the existing
    /// bounded simultaneous shear solve, without an outer nonlinear solve.
    pub fn new(
        contact: FiniteGapPlasticLayerContact,
        config: PlasticLayerShearDriveConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticLayerShearDriveError> {
        poll(cx)?;
        if !config.spin_length_m.is_finite()
            || config.spin_length_m <= 0.0
            || !config.maximum_rotation_rad.is_finite()
            || config.maximum_rotation_rad <= 0.0
            || config.maximum_rotation_rad > 0.1
            || !config.generalized_force_tolerance_n.is_finite()
            || config.generalized_force_tolerance_n <= 0.0
            || !config.energy_tolerance_j.is_finite()
            || config.energy_tolerance_j <= 0.0
        {
            return Err(PlasticLayerShearDriveError::Invalid(
                "spring domain/tolerances",
            ));
        }
        let k = config.stiffness_n_per_m;
        let mut scale = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                if !k[i][j].is_finite() || k[i][j] != k[j][i] {
                    return Err(PlasticLayerShearDriveError::Invalid(
                        "finite exactly reciprocal spring stiffness",
                    ));
                }
                scale = scale.max(k[i][j].abs());
            }
        }
        if scale == 0.0 {
            return Err(PlasticLayerShearDriveError::Invalid(
                "positive spring stiffness",
            ));
        }
        let normalized: Vec<_> = k.iter().flatten().map(|v| v / scale).collect();
        let factor = cholesky_cancellable(&normalized, 3, || cx.checkpoint().is_err())
            .map_err(factor_error)?;
        let root = scale.sqrt();
        let mut lower = [[0.0; 3]; 3];
        let mut inverse_lower = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..=i {
                lower[i][j] = finite(factor.l(i, j) * root)?;
            }
            for j in 0..=i {
                let mut rhs = if i == j { 1.0 } else { 0.0 };
                for m in 0..i {
                    rhs -= lower[i][m] * inverse_lower[m][j];
                }
                inverse_lower[i][j] = finite(rhs / lower[i][i])?;
            }
        }
        let mut inverse = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..=i {
                let value = finite(
                    (0..3)
                        .map(|m| inverse_lower[m][i] * inverse_lower[m][j])
                        .sum(),
                )?;
                inverse[i][j] = value;
                inverse[j][i] = value;
            }
        }
        let grid = contact.normal_model().config().grid;
        let mut rows = Vec::with_capacity(2 * contact.sites().len());
        for site in contact.sites() {
            poll(cx)?;
            let (x, y) = grid.center(site.cell_index);
            rows.push([1.0, 0.0, finite(-y / config.spin_length_m)?]);
            rows.push([0.0, 1.0, finite(x / config.spin_length_m)?]);
        }
        let original = contact.mechanical_response();
        let n = rows.len();
        let projected: Vec<_> = rows.iter().map(|&a| apply(inverse_lower, a)).collect();
        let mut compliance = original.compliance_m_per_n().to_vec();
        for i in 0..n {
            poll(cx)?;
            for j in 0..=i {
                // Gram form preserves exact numerical reciprocity of the
                // added physical spring compliance, including its cross terms.
                let value = finite(2.0 * dot(projected[i], projected[j]))?;
                compliance[i * n + j] = finite(compliance[i * n + j] + value)?;
                if i != j {
                    compliance[j * n + i] = finite(compliance[j * n + i] + value)?;
                }
            }
        }
        let response = PreparedTangentialContacts::new(
            original.contact_count(),
            compliance,
            original.config(),
            || cx.checkpoint().is_err(),
        )
        .map_err(|e| PlasticLayerShearDriveError::Contact(FiniteGapPlasticLayerError::Shear(e)))?;
        poll(cx)?;
        Ok(Self {
            owner: Arc::new(DriveOwner {
                contact,
                config,
                lower,
                inverse,
                rows,
                response,
            }),
        })
    }

    /// Original normal, layers, cell registration and body response.
    pub fn contact(&self) -> &FiniteGapPlasticLayerContact {
        &self.owner.contact
    }

    /// Original physical spring and explicit tolerances.
    pub fn config(&self) -> PlasticLayerShearDriveConfig {
        self.owner.config
    }

    /// Virgin spring and layers, retaining the native unloaded normal state.
    /// Load the mapped cells before any positive-time shear interval.
    pub fn initial_state(&self) -> PlasticLayerShearDriveState {
        self.bind(self.contact().initial_state())
    }

    /// Start unstrained layers and a zero-extension spring on accepted dents.
    /// This is not a transfer/reset of an already sheared composite checkpoint.
    pub fn state_from_normal(
        &self,
        normal: FiniteGapPlasticState,
        cx: &Cx<'_>,
    ) -> Result<PlasticLayerShearDriveState, PlasticLayerShearDriveError> {
        Ok(self.bind(self.contact().state_from_normal(normal, cx)?))
    }

    fn bind(&self, contact: FiniteGapPlasticLayerState) -> PlasticLayerShearDriveState {
        PlasticLayerShearDriveState {
            owner: Arc::clone(&self.owner),
            contact,
            extension_m: [0.0; 3],
            carriage_m: [0.0; 3],
        }
    }

    /// Prescribed normal approach, retaining every spring and shear coordinate.
    pub fn advance_normal(
        &self,
        before: &PlasticLayerShearDriveState,
        approach_m: f64,
        cx: &Cx<'_>,
    ) -> Result<PlasticLayerShearDriveNormalAdvance, PlasticLayerShearDriveError> {
        self.normal_with(before, cx, |state| {
            self.contact().advance_normal(state, approach_m, cx)
        })
    }

    /// Native normal force control, without an intermediate shear/load step.
    pub fn advance_normal_at_load(
        &self,
        before: &PlasticLayerShearDriveState,
        target_force_n: f64,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        PlasticLayerShearDriveNormalAdvance<FiniteGapPlasticLoadReceipt>,
        PlasticLayerShearDriveError,
    > {
        self.normal_with(before, cx, |state| {
            self.contact()
                .advance_normal_at_load(state, target_force_n, controls, cx)
        })
    }

    /// Native normal spring control, independent of the retained shear spring.
    pub fn advance_normal_with_drive(
        &self,
        before: &PlasticLayerShearDriveState,
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
        cx: &Cx<'_>,
    ) -> Result<
        PlasticLayerShearDriveNormalAdvance<FiniteGapPlasticDriveReceipt>,
        PlasticLayerShearDriveError,
    > {
        self.normal_with(before, cx, |state| {
            self.contact()
                .advance_normal_with_drive(state, drive, controls, cx)
        })
    }

    fn normal_with<B>(
        &self,
        before: &PlasticLayerShearDriveState,
        cx: &Cx<'_>,
        step: impl FnOnce(
            &FiniteGapPlasticLayerState,
        )
            -> Result<FiniteGapPlasticLayerNormalAdvance<B>, FiniteGapPlasticLayerError>,
    ) -> Result<PlasticLayerShearDriveNormalAdvance<B>, PlasticLayerShearDriveError> {
        poll(cx)?;
        self.admit(before)?;
        let next = step(&before.contact)?;
        poll(cx)?;
        Ok(PlasticLayerShearDriveNormalAdvance {
            candidate: PlasticLayerShearDriveState {
                owner: Arc::clone(&self.owner),
                contact: next.candidate,
                extension_m: before.extension_m,
                carriage_m: before.carriage_m,
            },
            normal: next.normal,
            boundary: next.boundary,
        })
    }

    /// Advance the held-normal fixture under a prescribed stage increment.
    ///
    /// Stage components are [m,m,rad]. The separate per-site body input is the
    /// original body's released motion with no contact traction [m]; it excludes
    /// this carriage and spring. Zero is appropriate for an otherwise stationary
    /// response. All cells, material trials and the spring share the supplied
    /// positive interval. Body and thermal owners must jointly accept the result.
    pub fn advance(
        &self,
        before: &PlasticLayerShearDriveState,
        stage_increment_m_rad: [f64; 3],
        released_body_increment_m: &[[f64; 2]],
        duration_s: f64,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<PlasticLayerShearDriveAdvance, PlasticLayerShearDriveError> {
        poll(cx)?;
        self.admit(before)?;
        let count = self.contact().sites().len();
        if released_body_increment_m.len() != count {
            return Err(PlasticLayerShearDriveError::Invalid(
                "one original body increment per site",
            ));
        }
        let mut stage = stage_increment_m_rad;
        for value in stage {
            finite(value)?;
        }
        stage[2] = finite(stage[2] * self.owner.config.spin_length_m)?;
        if stage_increment_m_rad[2] != 0.0 && stage[2] == 0.0 {
            return Err(PlasticLayerShearDriveError::Invalid(
                "representable stage twist",
            ));
        }
        let mut released = Vec::with_capacity(count);
        let mut free = [0.0; 3];
        for j in 0..3 {
            free[j] = finite(stage[j] + 2.0 * before.extension_m[j])?;
        }
        for (i, body) in released_body_increment_m.iter().enumerate() {
            poll(cx)?;
            let mut local = [0.0; 2];
            for j in 0..2 {
                finite(body[j])?;
                local[j] = finite(body[j] + dot(self.owner.rows[2 * i + j], free))?;
            }
            released.push(local);
        }
        let next = self.contact().advance_shear_with_response(
            &before.contact,
            &released,
            duration_s,
            controls,
            &self.owner.response,
            cx,
        )?;
        let wrench = [
            next.shear_force_n[0],
            next.shear_force_n[1],
            finite(next.spin_moment_n_m / self.owner.config.spin_length_m)?,
        ];
        let half_extension = apply(self.owner.inverse, wrench);
        let mut extension = [0.0; 3];
        let mut delta_extension = [0.0; 3];
        let mut carriage_increment = [0.0; 3];
        let mut carriage = [0.0; 3];
        let mut midpoint = [0.0; 3];
        for j in 0..3 {
            extension[j] = finite(2.0 * half_extension[j] - before.extension_m[j])?;
            delta_extension[j] = finite(extension[j] - before.extension_m[j])?;
            carriage_increment[j] = finite(stage[j] - delta_extension[j])?;
            carriage[j] = finite(before.carriage_m[j] + carriage_increment[j])?;
            midpoint[j] = 0.5 * extension[j] + 0.5 * before.extension_m[j];
        }
        finite(extension[2] / self.owner.config.spin_length_m)?;
        finite(carriage_increment[2] / self.owner.config.spin_length_m)?;
        let angle = finite(carriage[2] / self.owner.config.spin_length_m)?;
        if angle.abs() > self.owner.config.maximum_rotation_rad {
            return Err(PlasticLayerShearDriveError::Rotation {
                value_rad: angle,
                maximum_rad: self.owner.config.maximum_rotation_rad,
            });
        }
        let spring_force = apply(self.owner.config.stiffness_n_per_m, midpoint);
        let mut force_residual = [0.0; 3];
        for j in 0..3 {
            force_residual[j] = finite(spring_force[j] - wrench[j])?;
        }
        if force_residual
            .iter()
            .any(|r| r.abs() > self.owner.config.generalized_force_tolerance_n)
        {
            return Err(PlasticLayerShearDriveError::Equilibrium {
                residual_n: force_residual,
                tolerance_n: self.owner.config.generalized_force_tolerance_n,
            });
        }
        let stored_change = finite(dot(
            transpose_apply(self.owner.lower, midpoint),
            transpose_apply(self.owner.lower, delta_extension),
        ))?;
        let stage_work = finite(dot(wrench, stage))?;
        let forces: Vec<_> = next
            .shear
            .contacts
            .iter()
            .flat_map(|c| c.traction_n)
            .collect();
        let n = forces.len();
        let matrix = self.contact().mechanical_response().compliance_m_per_n();
        let mut body_work = 0.0;
        let mut correction = 0.0;
        for i in 0..n {
            poll(cx)?;
            let recoil: f64 = matrix[i * n..(i + 1) * n]
                .iter()
                .zip(&forces)
                .map(|(a, b)| a * b)
                .sum();
            let term = finite(forces[i] * (recoil - released_body_increment_m[i / 2][i % 2]))?;
            let value = term - correction;
            let sum = body_work + value;
            correction = (sum - body_work) - value;
            body_work = finite(sum)?;
        }
        let residual = finite(
            stage_work
                - stored_change
                - body_work
                - next.shear.stored_change_j
                - next.shear.dissipated_j,
        )?;
        if residual.abs() > self.owner.config.energy_tolerance_j {
            return Err(PlasticLayerShearDriveError::Energy {
                residual_j: residual,
                tolerance_j: self.owner.config.energy_tolerance_j,
            });
        }
        let candidate = PlasticLayerShearDriveState {
            owner: Arc::clone(&self.owner),
            contact: next.candidate,
            extension_m: extension,
            carriage_m: carriage,
        };
        if !candidate.spring_energy_j().is_finite() {
            return Err(PlasticLayerShearDriveError::Invalid("finite spring energy"));
        }
        poll(cx)?;
        Ok(PlasticLayerShearDriveAdvance {
            candidate,
            shear: next.shear,
            pressure_pa: next.pressure_pa,
            dissipation: next.dissipation,
            carriage_increment_m_rad: unscale(carriage_increment, self.owner.config.spin_length_m),
            mean_wrench_n_n_m: [wrench[0], wrench[1], next.spin_moment_n_m],
            spring_equilibrium_residual_n: force_residual,
            stage_work_j: stage_work,
            spring_stored_change_j: stored_change,
            body_response_work_j: body_work,
            balance_residual_j: residual,
        })
    }

    fn admit(
        &self,
        before: &PlasticLayerShearDriveState,
    ) -> Result<(), PlasticLayerShearDriveError> {
        if !Arc::ptr_eq(&self.owner, &before.owner) && self.owner != before.owner {
            return Err(PlasticLayerShearDriveError::Invalid(
                "original loading spring owner",
            ));
        }
        self.contact().admit(&before.contact)?;
        Ok(())
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn apply(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    a.map(|row| dot(row, x))
}
fn transpose_apply(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|j| (0..3).map(|i| a[i][j] * x[i]).sum())
}
fn unscale(mut x: [f64; 3], length: f64) -> [f64; 3] {
    x[2] /= length;
    x
}
fn finite(x: f64) -> Result<f64, PlasticLayerShearDriveError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(PlasticLayerShearDriveError::Invalid(
            "finite spring arithmetic",
        ))
    }
}
fn poll(cx: &Cx<'_>) -> Result<(), PlasticLayerShearDriveError> {
    super::checkpoint(cx).map_err(|e| PlasticLayerShearDriveError::Contact(e.into()))
}
fn factor_error(error: FactorRunError) -> PlasticLayerShearDriveError {
    match error {
        FactorRunError::Cancelled => PlasticLayerShearDriveError::Contact(
            FiniteGapPlasticLayerError::Shear(MultipleTangentialError::Cancelled),
        ),
        FactorRunError::Factor(_) => {
            PlasticLayerShearDriveError::Invalid("positive-definite loading spring")
        }
    }
}
