//! Low-Tabor DMT adhesion with the original evolving-area Mindlin microslip.
//!
//! The original DMT-M normal owner retains Hertz repulsion on the real contact
//! disk and adhesion outside it. Only that REPULSIVE pressure enters the
//! explicitly sourced constant Coulomb law. Its existing analytic MDR owner
//! evolves elastic shear, nested reversals, contact growth and unloading heat.
//! Signed net normal force never replaces repulsion in the shear capacity.
//!
//! This selects decoupled, isotropic, proportional scalar shear, with constant
//! material/interface data at fixed external conditions. It is not JKR partial
//! slip, a constant-junction-strength model, torsion or shear-induced peeling.
//! Signed affine paths retain exact contact-entry/exit timing. Zero-area gap
//! motion carries no shear, releases no additional heat, and remains in the
//! accumulated displacement before subsequent stress-free recontact.
//!
//! The normal DMT approximation follows Maugis (1992),
//! DOI 10.1016/0021-9797(92)90285-T; the exact retained Hertz--Mindlin/MDR
//! approximation is documented by `fs_tribo::mindlin::loading`.

use super::{DmtError, DmtResponse, SourceDmtSphereContact};
use fs_blake3::{ContentHash, DomainHasher};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedInterfaceStatePoint,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_interface_state_point,
};
use fs_qty::semantic::QuantityKind;
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::mindlin::loading::{
    HertzMindlinContact, HertzMindlinParameters, HertzMindlinState,
    path::{MindlinPathError, MindlinPathTrial},
};
use fs_tribo::{TangentialLaw, TangentialStep, TangentialStepInput, TriboError};

pub use crate::interface_binding::mindlin::loading::path::SourceMindlinPathInput;

/// Original adhesion, elasticity and ordered friction evidence on one owner.
#[derive(Clone, Debug, PartialEq)]
pub struct SourceDmtMindlinContact {
    normal: SourceDmtSphereContact,
    friction: ResolvedInterfaceStatePoint,
    law: HertzMindlinContact,
    temperature_kelvin: f64,
    identity: ContentHash,
}

/// Complete contact/gap history; fields are sealed against source swaps.
#[derive(Clone, Debug, PartialEq)]
pub struct DmtMindlinState {
    identity: ContentHash,
    normal: DmtResponse,
    constitutive: HertzMindlinState,
    tangential_energy_j: f64,
}
impl DmtMindlinState {
    pub const fn indentation_m(&self) -> f64 {
        self.normal.indentation_m
    }
    pub const fn normal_response(&self) -> DmtResponse {
        self.normal
    }
    pub const fn contact_radius_m(&self) -> f64 {
        self.normal.contact_radius_m
    }
    pub const fn tangential_energy_j(&self) -> f64 {
        self.tangential_energy_j
    }
    pub const fn slip_m(&self) -> f64 {
        self.constitutive.slip_m()
    }
    pub fn endpoint_force_n(&self) -> f64 {
        self.constitutive.endpoint_force_n()
    }
    pub fn history_intervals(&self) -> usize {
        self.constitutive.history_intervals()
    }
    /// Read-only original analytic memory; no numerical field substitutes for it.
    pub const fn constitutive_state(&self) -> &HertzMindlinState {
        &self.constitutive
    }
}

/// Allocation-free mechanical probe; final history admission remains mandatory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DmtMindlinPathTrial {
    pub tangential: MindlinPathTrial,
    /// Fraction of the original affine path with positive real contact area.
    /// Retain directly: subtracting event endpoints can lose short intervals.
    pub contact_fraction: f64,
    /// Normalized contact interval, [0, 0] when inactive. Boundary instants
    /// may have zero area; the independent duration above retains small spans.
    pub contact_interval_fraction: [f64; 2],
    pub normal_response: DmtResponse,
    pub mean_normal_force_n: f64,
    pub normal_work_j: f64,
    pub normal_stored_change_j: f64,
}

/// Immutable complete proposal, with separate conservative normal work and heat.
#[derive(Clone, Debug, PartialEq)]
pub struct DmtMindlinPathStep {
    pub tangential: TangentialStep<DmtMindlinState>,
    pub contact_fraction: f64,
    pub contact_interval_fraction: [f64; 2],
    pub normal_response: DmtResponse,
    pub mean_normal_force_n: f64,
    pub normal_work_j: f64,
    pub normal_stored_change_j: f64,
}

/// Refusal preserves the original source and accepted history.
#[derive(Clone, Debug, PartialEq)]
pub enum DmtMindlinError {
    Normal(DmtError),
    Material(MaterialStatePointError),
    Path(MindlinPathError),
    Invalid(&'static str),
    StateMismatch,
}
impl From<DmtError> for DmtMindlinError {
    fn from(e: DmtError) -> Self {
        Self::Normal(e)
    }
}
impl From<MaterialStatePointError> for DmtMindlinError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<MindlinPathError> for DmtMindlinError {
    fn from(e: MindlinPathError) -> Self {
        Self::Path(e)
    }
}
impl From<TriboError> for DmtMindlinError {
    fn from(e: TriboError) -> Self {
        Self::Path(MindlinPathError::Constitutive(e))
    }
}
impl core::fmt::Display for DmtMindlinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bound DMT--Mindlin path refused: {self:?}")
    }
}
impl std::error::Error for DmtMindlinError {}

impl SourceDmtMindlinContact {
    /// Resolve one explicit local Coulomb coefficient from the SAME ordered
    /// interface and query as the original DMT source. It governs both micro-
    /// and gross slip; unequal static/kinetic coefficients are not combined.
    ///
    /// This reuses the native elastic/history law directly. The nonadhesive
    /// SourceHertzMindlinContact cannot admit W>0 and is never bypassed by
    /// fabricating a zero-adhesion material or replacing the source identity.
    pub fn resolve(
        normal: SourceDmtSphereContact,
        friction_requirement: ScalarPropertyRequirement,
        selection: MaterialPropertySelection,
        maximum_history_intervals: usize,
    ) -> Result<Self, DmtMindlinError> {
        if friction_requirement.quantity() != QuantitySpec::dimensional(Dims::NONE) {
            return Err(DmtMindlinError::Invalid(
                "local Coulomb coefficient must be dimensionless",
            ));
        }
        let original = normal.normal_state();
        let elastic = original.elastic();
        let friction = resolve_interface_state_point(
            elastic.interface().card(),
            &original.interface_state().query(),
            &[friction_requirement],
            selection,
        )?;
        for property in elastic
            .surface_a_state()
            .resolved()
            .properties()
            .iter()
            .chain(elastic.surface_b_state().resolved().properties())
            .chain(original.interface_state().properties())
            .chain(friction.properties())
        {
            if !property.is_static_at_fixed_external_conditions() {
                return Err(DmtMindlinError::Invalid(
                    "DMT--Mindlin requires fixed external source conditions",
                ));
            }
            let receipt = &property.answer().receipt;
            for (axis, value) in &receipt.query_point {
                let quantity = receipt.axis_quantities.get(axis).copied();
                let pressure = quantity
                    .and_then(|q| q.semantic_type())
                    .is_some_and(|ty| ty.kind() == QuantityKind::Pressure);
                if pressure || !is_fixed_external_coordinate(axis, quantity, Some(*value)) {
                    return Err(DmtMindlinError::Invalid(
                        "DMT--Mindlin cannot freeze pressure, rate, load, footprint or history queries",
                    ));
                }
            }
        }
        let temperature_kelvin = elastic
            .state_coordinate("T")
            .ok_or(DmtMindlinError::Invalid(
                "original DMT source has no contact temperature",
            ))?;
        let solids = [elastic.surface_a_state(), elastic.surface_b_state()];
        let law = HertzMindlinContact::new(HertzMindlinParameters {
            effective_radius_m: normal.effective_radius_m(),
            shear_moduli_pa: solids.map(|s| 0.5 * s.young_modulus_pa() / (1.0 + s.poisson_ratio())),
            poisson_ratios: solids.map(|s| s.poisson_ratio()),
            friction_coefficient: friction.properties()[0].value_si(),
            maximum_history_intervals,
        })?;
        same(law.normal_modulus_pa(), elastic.reduced_modulus_pa())?;
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.dmt-mindlin.v1");
        hash.update(normal.identity().as_bytes());
        hash.update(friction.identity().as_bytes());
        hash.update(&(maximum_history_intervals as u64).to_le_bytes());
        Ok(Self {
            normal,
            friction,
            law,
            temperature_kelvin,
            identity: hash.finalize(),
        })
    }

    pub const fn normal_source(&self) -> &SourceDmtSphereContact {
        &self.normal
    }
    pub const fn friction_source(&self) -> &ResolvedInterfaceStatePoint {
        &self.friction
    }
    pub const fn law(&self) -> &HertzMindlinContact {
        &self.law
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }

    /// Initialize an admitted signed normal state with zero shear.
    /// Initial normal storage is present, not newly supplied mechanical work.
    /// This cannot consume or erase an existing loaded shear checkpoint.
    pub fn initialize(&self, indentation_m: f64) -> Result<DmtMindlinState, DmtMindlinError> {
        let normal = self.normal.evaluate(indentation_m)?;
        let constitutive = self
            .law
            .normal_step(&self.law.initial_state(), normal.repulsive_normal_force_n)?
            .state;
        match_normal(&constitutive, &normal)?;
        Ok(DmtMindlinState {
            identity: self.identity,
            normal,
            constitutive,
            tangential_energy_j: 0.0,
        })
    }

    /// Trial the original source and scalar recoil solve without allocating
    /// history. Outer mechanics may use this in a bounded normal/shear solve.
    /// Normal work includes adhesion once; shear heat is native microslip only.
    pub fn trial_coupled_path(
        &self,
        state: &DmtMindlinState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DmtMindlinPathTrial, DmtMindlinError> {
        let (normal, native) = self.admit(state, input, cancelled)?;
        let path = self.law.trial_coupled_gap_path(
            &state.constitutive,
            &native,
            [state.indentation_m(), input.indentation_m],
            cancelled,
        )?;
        let (mean, work) = self.normal_work(state, input)?;
        poll(cancelled)?;
        Ok(DmtMindlinPathTrial {
            tangential: path.tangential,
            contact_fraction: path.contact_fraction,
            contact_interval_fraction: path.contact_interval_fraction,
            normal_response: normal,
            mean_normal_force_n: mean,
            normal_work_j: work,
            normal_stored_change_j: work,
        })
    }

    /// Advance simultaneous affine indentation and scalar shear, including
    /// unloading/reversal memory. Optional scalar body recoil uses the SAME
    /// exact-work average force. Positive or tensile net normal effort is
    /// admitted; only real repulsive contact carries shear. Signed paths use
    /// the exact positive-area interval of the original affine indentation.
    /// Gap drift remains in accumulated slip, and its attraction contributes
    /// only the original conservative DMT work. Heat and shear storage change
    /// are never multiplied by the fraction of the step spent in contact.
    ///
    /// The immutable result is a candidate only. Complete native mechanics
    /// owns normal motion, frame transport, step/rate limits and publication.
    pub fn advance_coupled_path(
        &self,
        state: &DmtMindlinState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DmtMindlinPathStep, DmtMindlinError> {
        let (normal, native) = self.admit(state, input, cancelled)?;
        let path = self.law.advance_coupled_gap_path(
            &state.constitutive,
            &native,
            [state.indentation_m(), input.indentation_m],
            cancelled,
        )?;
        let (mean, work) = self.normal_work(state, input)?;
        let candidate = path.path.tangential;
        match_normal(&candidate.state, &normal)?;
        let energy = self.law.stored_energy_j(&candidate.state);
        if !energy.is_finite() || energy < 0.0 {
            return Err(DmtMindlinError::Invalid(
                "unrepresentable native shear storage",
            ));
        }
        poll(cancelled)?;
        Ok(DmtMindlinPathStep {
            tangential: TangentialStep {
                traction_n: candidate.traction_n,
                slip_increment_m: candidate.slip_increment_m,
                state: DmtMindlinState {
                    identity: self.identity,
                    normal,
                    constitutive: candidate.state,
                    tangential_energy_j: energy,
                },
                regime: candidate.regime,
                stored_change_j: candidate.stored_change_j,
                dissipated_j: candidate.dissipated_j,
                closure_residual_j: candidate.closure_residual_j,
                compatibility_residual_m: candidate.compatibility_residual_m,
                iterations: candidate.iterations,
            },
            contact_fraction: path.contact_fraction,
            contact_interval_fraction: path.contact_interval_fraction,
            normal_response: normal,
            mean_normal_force_n: mean,
            normal_work_j: work,
            normal_stored_change_j: work,
        })
    }

    /// Radial Hertz repulsive pressure, original scalar shear traction and
    /// relative elastic displacement; exterior DMT attraction is not pressure
    /// in the solid-contact friction law.
    pub fn surface_response(
        &self,
        state: &DmtMindlinState,
        radius_m: f64,
    ) -> Result<fs_tribo::mindlin::surface::MindlinSurfaceResponse, DmtMindlinError> {
        self.owns(state)?;
        Ok(self.law.surface_response(&state.constitutive, radius_m)?)
    }

    fn owns(&self, state: &DmtMindlinState) -> Result<(), DmtMindlinError> {
        if state.identity != self.identity
            || !self.law.stored_energy_j(&state.constitutive).is_finite()
        {
            Err(DmtMindlinError::StateMismatch)
        } else {
            Ok(())
        }
    }

    fn admit(
        &self,
        state: &DmtMindlinState,
        input: &SourceMindlinPathInput,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(DmtResponse, TangentialStepInput), DmtMindlinError> {
        poll(cancelled)?;
        self.owns(state)?;
        let normal = self.normal.evaluate(input.indentation_m)?;
        if !input.dt_s.is_finite()
            || input.dt_s <= 0.0
            || !((input.indentation_m - state.indentation_m()) / input.dt_s).is_finite()
        {
            return Err(DmtMindlinError::Invalid(
                "positive duration and finite actual normal rate required",
            ));
        }
        Ok((
            normal,
            TangentialStepInput {
                normal_n: normal.repulsive_normal_force_n,
                slip_increment_m: [input.free_slip_increment_m, 0.0],
                dt_s: input.dt_s,
                contact_temperature_k: Some(self.temperature_kelvin),
                nominal_pressure_pa: None,
                compliance_m_per_n: input.compliance_m_per_n,
                solve: input.solve,
            },
        ))
    }

    fn normal_work(
        &self,
        state: &DmtMindlinState,
        input: &SourceMindlinPathInput,
    ) -> Result<(f64, f64), DmtMindlinError> {
        Ok((
            self.normal
                .mean_force_between(state.indentation_m(), input.indentation_m)?,
            self.normal
                .reversible_work_between(state.indentation_m(), input.indentation_m)?,
        ))
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DmtMindlinError> {
    if cancelled() {
        Err(MindlinPathError::Cancelled.into())
    } else {
        Ok(())
    }
}
fn same(actual: f64, expected: f64) -> Result<(), DmtMindlinError> {
    let scale = actual.abs().max(expected.abs()).max(f64::MIN_POSITIVE);
    if !actual.is_finite()
        || !expected.is_finite()
        || (actual - expected).abs() > 512.0 * f64::EPSILON * scale
    {
        Err(DmtMindlinError::Invalid(
            "native Mindlin must retain original DMT repulsive geometry and elasticity",
        ))
    } else {
        Ok(())
    }
}
fn match_normal(native: &HertzMindlinState, normal: &DmtResponse) -> Result<(), DmtMindlinError> {
    for (actual, expected) in [
        (native.indentation_m(), normal.indentation_m.max(0.0)),
        (native.normal_force_n(), normal.repulsive_normal_force_n),
        (native.contact_radius_m(), normal.contact_radius_m),
        (
            native.normal_energy_j(),
            0.4 * normal.repulsive_normal_force_n * normal.indentation_m.max(0.0),
        ),
    ] {
        same(actual, expected)?;
    }
    Ok(())
}
