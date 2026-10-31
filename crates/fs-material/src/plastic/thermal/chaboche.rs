//! Sourced nonisothermal Chaboche history on the existing thermal energy owner.
//!
//! A branch retains the strain-like coordinate `z = 3 X/(2 C)` when temperature
//! changes: `X* = C(T1) X(T0)/C(T0)`. This is the integrated `(C'/C) X dT`
//! contribution of the nonisothermal hardening law; keeping stress-like X fixed
//! would silently change the internal coordinate. The subsequent mechanical
//! increment uses the existing implicit Chaboche return and its exact tangent.
//!
//! Branch free energy is `C z:z/3`, mechanical entropy is `-C' z:z/3`, and
//! frozen-history coefficient work is `(C1-C0) z0:z0/3`. Physical recovery and
//! numerical loss use `X1-X*`, never the thermal rescaling. C is globally affine
//! and strictly positive, or identically zero; gamma is fixed. These restrictions
//! preserve the native saturation bound and add no frozen-history second-T
//! derivative to the existing caloric capacity admission. An explicitly sourced
//! linear or power-law overstress uses the actual physical interval in the same
//! implicit branch return. Kinetic coefficients store no energy. Static recovery
//! and finite-deformation plasticity are not represented.
//!
//! Nonisothermal model family and the positive temperature-rate term (Eq. 9):
//! <https://neml.readthedocs.io/en/stable/hardening/non/chaboche.html>.
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-combinedhardening.htm>.
//! The definitions above use FrankenSim's positive backstress convention.

use super::expansion::ThermalHistoryContribution;
use super::*;
use crate::plastic::chaboche::rate::ChabocheOverstressJ2;
use crate::plastic::chaboche::source::{BackstressProperties, chaboche_branch_requirements};
use crate::plastic::chaboche::{
    BackstressBranch, ChabocheError, ChabocheJ2, ChabocheResponse, ChabocheState,
};
use crate::tensor::contract;
use crate::{MaterialAdmissibility, Tangent6};
use fs_evidence::{Ambition, ModelCard};

/// Complete immutable thermal source plus independently named branch histories.
/// Array order is history order. No branch aggregation or padding occurs.
#[derive(Clone, Debug)]
pub struct ThermalChabocheSource<'a, const N: usize = 1> {
    base: ThermoplasticSource<'a>,
    properties: [BackstressProperties; N],
    reference_moduli_pa: [f64; N],
    modulus_slopes_pa_k: [f64; N],
    recovery_per_strain: [f64; N],
}

/// One temperature transition and its endpoint mechanical law. Incoming
/// ChabocheState belongs to the OLD temperature; returned history belongs to
/// the NEW temperature. Repeated Newton trials must reuse the same old history.
/// Use one accepted physical increment per transition. Resolve a fresh law
/// before advancing the returned state; reusing this old-to-new transport would
/// apply the coefficient change twice. A rate-dependent response also advances
/// plastic time even at unchanged strain: inspect its retained accepted stress,
/// not another return-map evaluation. A new-to-new transition is appropriate for
/// the next physical hold interval, or rate-independent endpoint inspection.
#[derive(Clone, Debug)]
pub struct ThermalChaboche<const N: usize = 1> {
    source: ResolvedJ2Plasticity,
    previous_law: ChabocheJ2<N>,
    law: ChabocheJ2<N>,
    rate: Option<ChabocheOverstressJ2<N>>,
    expansion_strain: f64,
}

/// Audited native candidate and the original thermal energy accounts.
#[derive(Clone, Debug)]
pub struct ThermalChabocheIncrement<const N: usize = 1> {
    /// Complete endpoint history, including each independently retained branch.
    pub state: ChabocheState<N>,
    /// The native mechanical accounts include branch storage and recovery;
    /// `mechanical.state` is exactly `state.plastic`. Thermal coefficient work
    /// and entropy include both the backbone and every active branch.
    pub thermal: ThermoplasticIncrement,
    /// Already included in physical dissipation; do not add it to heat again.
    pub recovery_dissipation_j_m3: f64,
}

/// Preserve the original thermal and native nonlinear return diagnostics.
#[derive(Debug)]
pub enum ThermalChabocheError {
    Thermal(FeedbackError),
    Mechanical(ChabocheError),
}
impl fmt::Display for ThermalChabocheError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Thermal(error) => write!(f, "thermal Chaboche: {error}"),
            Self::Mechanical(error) => write!(f, "thermal Chaboche: {error}"),
        }
    }
}
impl std::error::Error for ThermalChabocheError {}
impl From<FeedbackError> for ThermalChabocheError {
    fn from(error: FeedbackError) -> Self {
        Self::Thermal(error)
    }
}
impl From<ChabocheError> for ThermalChabocheError {
    fn from(error: ChabocheError) -> Self {
        Self::Mechanical(error)
    }
}

impl<'a> ThermoplasticSource<'a> {
    /// Add one through eight independently selected Chaboche branches to an
    /// explicitly selected Voce backbone. Full source lines, including interior
    /// knots and surface crossings, must admit affine C and constant gamma.
    pub fn with_chaboche_branches<const N: usize>(
        self,
        selection: MaterialPropertySelection,
        properties: [BackstressProperties; N],
    ) -> Result<ThermalChabocheSource<'a, N>, FeedbackError> {
        ThermalChabocheSource::try_new(self, selection, properties)
    }
}

impl<'a, const N: usize> ThermalChabocheSource<'a, N> {
    /// Bind branch source roles to the same immutable material and caloric span.
    /// A branch may be explicitly absent (C=0 everywhere). An active branch may
    /// not reach zero, because stress-like history would then lose z. No branch
    /// coefficient may alias an already selected mechanical or thermal role.
    pub fn try_new(
        mut base: ThermoplasticSource<'a>,
        selection: MaterialPropertySelection,
        properties: [BackstressProperties; N],
    ) -> Result<Self, FeedbackError> {
        if !base.has_voce_hardening() {
            return Err(FeedbackError::Invalid(
                "thermal Chaboche requires an explicitly selected Voce backbone",
            ));
        }
        let schema = chaboche_branch_requirements(base.source.hardening_basis(), &properties)
            .map_err(FeedbackError::Source)?;
        let requirements: Vec<_> = schema
            .into_iter()
            .filter(|requirement| {
                properties.iter().any(|names| {
                    requirement.name() == names.modulus || requirement.name() == names.recovery
                })
            })
            .collect();
        for requirement in &requirements {
            if base
                .requirements
                .iter()
                .any(|existing| existing.name() == requirement.name())
            {
                return Err(FeedbackError::Invalid(
                    "Chaboche branch roles must not alias selected thermal or mechanical quantities",
                ));
            }
        }
        let resolved =
            resolve_material_state_point(base.card, &base.context, &requirements, selection)
                .map_err(|e| FeedbackError::Source(e.into()))?;
        let knots = base.caloric.curve().knots();
        let first_temperature = knots[0].temperature_k;
        let last_temperature = knots[knots.len() - 1].temperature_k;
        let lower = point_at(&base.context, &base.temperature_axis, first_temperature)?;
        let upper = point_at(&base.context, &base.temperature_axis, last_temperature)?;
        let mut reference_moduli_pa = [0.0; N];
        let mut modulus_slopes_pa_k = [0.0; N];
        let mut recovery_per_strain = [0.0; N];
        let mut additional_pins = Vec::with_capacity(2 * N);
        for (i, names) in properties.iter().enumerate() {
            for (name, is_modulus) in [(&names.modulus, true), (&names.recovery, false)] {
                let property = resolved
                    .property(name)
                    .expect("resolved branch requirement");
                admit_rate_independent_property(property, &base.temperature_axis)?;
                let reference = property.value_si();
                let selected = property.answer().receipt.selected;
                let claim = base
                    .card
                    .claims()
                    .claim(selected)
                    .ok_or(FeedbackError::Invalid(
                        "Chaboche branch claim left its original material card",
                    ))?;
                let support = base
                    .card
                    .claims()
                    .query_line_typed(
                        &claim.key,
                        &base.temperature_axis,
                        &lower,
                        &upper,
                        ClaimSelection::Pinned(selected),
                    )
                    .map_err(|e| FeedbackError::Support {
                        property: name.clone(),
                        reason: e.to_string(),
                    })?;
                let mut samples = Vec::with_capacity(support.knots.len());
                for (temperature, sample) in support.knots {
                    let value = sample.evidence.value.value;
                    if !value.is_finite() || value < 0.0 {
                        return Err(FeedbackError::Invalid(
                            "Chaboche coefficients must remain finite and nonnegative throughout the source span",
                        ));
                    }
                    if !is_modulus && value != reference {
                        return Err(FeedbackError::Invalid(
                            "thermal Chaboche recovery gamma must remain constant throughout the source span",
                        ));
                    }
                    samples.push((temperature, value));
                }
                if is_modulus {
                    let absent = samples.iter().all(|(_, value)| *value == 0.0);
                    if !absent && samples.iter().any(|(_, value)| *value <= 0.0) {
                        return Err(FeedbackError::Invalid(
                            "an active Chaboche modulus must remain strictly positive throughout the thermal span",
                        ));
                    }
                    reference_moduli_pa[i] = reference;
                    modulus_slopes_pa_k[i] =
                        affine_temperature_slope(&samples, base.reference_temperature_k, reference)
                            .map_err(FeedbackError::Invalid)?;
                } else {
                    recovery_per_strain[i] = reference;
                }
                additional_pins.push((name.clone(), selected));
            }
        }
        let MaterialPropertySelection::PinnedByProperty(pins) = &mut base.selection else {
            return Err(FeedbackError::Invalid(
                "thermal source lost its original claim pins",
            ));
        };
        pins.extend(additional_pins);
        base.requirements.extend(requirements);
        base.source = base.resolve_temperature(base.reference_temperature_k)?;
        let binding = Self {
            base,
            properties,
            reference_moduli_pa,
            modulus_slopes_pa_k,
            recovery_per_strain,
        };
        // Admit the sum of every native return modulus at both affine extrema.
        for temperature in [first_temperature, last_temperature] {
            let source = binding.base.resolve_temperature(temperature)?;
            binding.law_at(&source)?;
        }
        Ok(binding)
    }

    /// Original caloric/expansion owner with the complete pinned branch bundle.
    pub const fn base(&self) -> &ThermoplasticSource<'a> {
        &self.base
    }

    /// Source order is exactly native history order.
    pub const fn properties(&self) -> &[BackstressProperties; N] {
        &self.properties
    }

    fn law_at(&self, source: &ResolvedJ2Plasticity) -> Result<ChabocheJ2<N>, FeedbackError> {
        let mut branches =
            [BackstressBranch::new(0.0, 0.0).map_err(|e| FeedbackError::Source(e.into()))?; N];
        for (i, names) in self.properties.iter().enumerate() {
            let modulus = source
                .material()
                .property(&names.modulus)
                .ok_or(FeedbackError::Invalid(
                    "resolved source omitted a branch modulus",
                ))?
                .value_si();
            let recovery = source
                .material()
                .property(&names.recovery)
                .ok_or(FeedbackError::Invalid(
                    "resolved source omitted branch recovery",
                ))?
                .value_si();
            if recovery != self.recovery_per_strain[i]
                || self.modulus_slopes_pa_k[i] == 0.0 && modulus != self.reference_moduli_pa[i]
                || self.reference_moduli_pa[i] > 0.0 && modulus <= 0.0
            {
                return Err(FeedbackError::Invalid(
                    "branch coefficients left their admitted thermal model",
                ));
            }
            branches[i] = BackstressBranch::new(modulus, recovery)
                .map_err(|e| FeedbackError::Source(e.into()))?;
        }
        let voce = self.base.voce_law(source)?.ok_or(FeedbackError::Invalid(
            "thermal Chaboche lost its selected Voce backbone",
        ))?;
        ChabocheJ2::new(voce, branches).map_err(|e| FeedbackError::Source(e.into()))
    }

    /// Resolve one thermal transition for native continuum assembly. The same
    /// old enthalpy and history must be retained during all endpoint iterations.
    pub fn resolve_enthalpy(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalChaboche<N>), FeedbackError> {
        self.resolve_for_interval(old_h_j_kg, new_h_j_kg, None)
    }

    /// Resolve the sourced joint branch/overstress return for exactly this
    /// physical interval. Newton iterations reuse its duration and old history;
    /// thermal proposals do not create additional plastic time.
    pub fn resolve_enthalpy_for_step(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        duration_s: f64,
    ) -> Result<(EquilibriumPhaseState, ThermalChaboche<N>), FeedbackError> {
        if !positive(duration_s) {
            return Err(FeedbackError::Invalid(
                "positive physical constitutive interval required",
            ));
        }
        self.resolve_for_interval(old_h_j_kg, new_h_j_kg, Some(duration_s))
    }

    fn resolve_for_interval(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        duration_s: Option<f64>,
    ) -> Result<(EquilibriumPhaseState, ThermalChaboche<N>), FeedbackError> {
        let (_, previous) = self.base.resolve_enthalpy(old_h_j_kg)?;
        let (phase, source) = self.base.resolve_enthalpy(new_h_j_kg)?;
        let previous_law = self.law_at(&previous)?;
        let law = self.law_at(&source)?;
        let rate = self.overstress_law(&source, &law, duration_s)?;
        let expansion_strain = self
            .base
            .expansion_strain_between(self.base.reference_temperature_k, phase.temperature_k())?;
        Ok((
            phase,
            ThermalChaboche {
                source,
                previous_law,
                law,
                rate,
                expansion_strain,
            },
        ))
    }

    fn overstress_law(
        &self,
        source: &ResolvedJ2Plasticity,
        law: &ChabocheJ2<N>,
        duration_s: Option<f64>,
    ) -> Result<Option<ChabocheOverstressJ2<N>>, FeedbackError> {
        use crate::plastic::rate::{
            OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY,
            OVERSTRESS_VISCOSITY_PROPERTY,
        };
        let Some(model) = self.base.overstress else {
            return Ok(None);
        };
        let duration_s = duration_s.ok_or(FeedbackError::Invalid(
            "source-defined overstress requires the physical step duration",
        ))?;
        let value = |name| {
            source
                .material()
                .property(name)
                .map(|p| p.value_si())
                .ok_or(FeedbackError::Invalid(
                    "resolved source omitted an overstress coefficient",
                ))
        };
        let viscosity = value(OVERSTRESS_VISCOSITY_PROPERTY)?;
        let rate = match model {
            OverstressSelection::Linear => {
                ChabocheOverstressJ2::new_linear(law.clone(), viscosity, duration_s)
            }
            OverstressSelection::PowerLaw => ChabocheOverstressJ2::new_power_law(
                law.clone(),
                viscosity,
                value(OVERSTRESS_STRESS_SCALE_PROPERTY)?,
                value(OVERSTRESS_EXPONENT_PROPERTY)?,
                duration_s,
            ),
        }
        .map_err(|error| FeedbackError::Source(error.into()))?;
        Ok(Some(rate))
    }

    /// Propose a total-strain/temperature increment without modifying accepted
    /// history. The original native mechanical and thermal energy gates apply.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_energy(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &ChabocheState<N>,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<ThermalChabocheIncrement<N>, ThermalChabocheError> {
        self.advance_for_interval(
            old_h_j_kg, new_h_j_kg, old_strain, old_state, strain, tolerances, None,
        )
    }

    /// Audit the same physical interval used by the spatial constitutive law.
    /// Yield and viscous flow plus dynamic recovery are physical dissipation;
    /// all energetic hardening and numerical integration losses stay separate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_energy_for_step(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &ChabocheState<N>,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: f64,
    ) -> Result<ThermalChabocheIncrement<N>, ThermalChabocheError> {
        if !positive(duration_s) {
            return Err(
                FeedbackError::Invalid("positive physical constitutive interval required").into(),
            );
        }
        self.advance_for_interval(
            old_h_j_kg,
            new_h_j_kg,
            old_strain,
            old_state,
            strain,
            tolerances,
            Some(duration_s),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_for_interval(
        &self,
        old_h_j_kg: f64,
        new_h_j_kg: f64,
        old_strain: &Voigt,
        old_state: &ChabocheState<N>,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        duration_s: Option<f64>,
    ) -> Result<ThermalChabocheIncrement<N>, ThermalChabocheError> {
        let old_phase = self
            .base
            .caloric
            .curve()
            .state_at_specific_enthalpy(old_h_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let (new_phase, resolved) =
            self.resolve_for_interval(old_h_j_kg, new_h_j_kg, duration_s)?;
        let old_mechanical = self
            .base
            .mechanical_strain(old_phase.temperature_k(), old_strain)?;
        let mechanical_strain = self
            .base
            .mechanical_strain(new_phase.temperature_k(), strain)?;
        resolved
            .previous_law
            .validate_history(&old_mechanical, old_state)?;
        let transported = resolved.transport_history(old_state)?;
        let mechanical = match &resolved.rate {
            Some(rate) => rate.advance_with_energy(
                &old_mechanical,
                &transported,
                &mechanical_strain,
                tolerances,
            ),
            None => resolved.law.advance_with_energy(
                &old_mechanical,
                &transported,
                &mechanical_strain,
                tolerances,
            ),
        }?;
        let mut contribution = ThermalHistoryContribution::default();
        for i in 0..N {
            let slope = self.modulus_slopes_pa_k[i];
            if slope == 0.0 {
                continue;
            }
            contribution.variable_storage = true;
            let c0 = resolved.previous_law.branches()[i].modulus_pa();
            let c1 = resolved.law.branches()[i].modulus_pa();
            // X/C is bounded by the admitted plastic history, avoiding C^2
            // overflow/underflow in the otherwise equivalent 3 X:X/(4 C^2).
            let y0 = old_state.backstress_pa[i].map(|x| x / c0);
            let y1 = mechanical.state.backstress_pa[i].map(|x| x / c1);
            let shape0 = 0.75 * contract(&y0, &y0);
            let shape1 = 0.75 * contract(&y1, &y1);
            contribution.old_entropy_j_m3_k -= slope * shape0;
            contribution.new_entropy_j_m3_k -= slope * shape1;
            contribution.coefficient_storage_change_j_m3 += (c1 - c0) * shape0;
        }
        let native = EnergyIncrement {
            strain: mechanical.strain,
            stress_pa: mechanical.stress_pa,
            state: mechanical.state.plastic.clone(),
            energy: mechanical.energy,
        };
        let thermal = self.base.finish_increment_at_temperature(
            old_phase.temperature_k(),
            new_phase.temperature_k(),
            old_strain,
            &old_state.plastic,
            strain,
            &resolved.source,
            native,
            contribution,
            tolerances,
        )?;
        Ok(ThermalChabocheIncrement {
            state: mechanical.state,
            thermal,
            recovery_dissipation_j_m3: mechanical.recovery_dissipation_j_m3,
        })
    }
}

impl<const N: usize> ThermalChaboche<N> {
    /// Complete immutable endpoint bundle, including all thermal/branch receipts.
    pub const fn resolved(&self) -> &ResolvedJ2Plasticity {
        &self.source
    }

    /// Original isothermal endpoint return map. Its history argument must first
    /// be transported with `transport_history` when temperatures differ.
    pub const fn law(&self) -> &ChabocheJ2<N> {
        &self.law
    }

    /// Joint endpoint branch/overstress law at its actual physical duration.
    /// History must first be thermally transported before direct native use.
    pub const fn rate_law(&self) -> Option<&ChabocheOverstressJ2<N>> {
        self.rate.as_ref()
    }

    /// Change only branch stress coordinates; native plastic strain and alpha
    /// remain unchanged. Constant coefficients preserve exact original bytes.
    pub fn transport_history(
        &self,
        state: &ChabocheState<N>,
    ) -> Result<ChabocheState<N>, ChabocheError> {
        self.previous_law.validate_history(&[0.0; 6], state)?;
        let mut transported = state.clone();
        for i in 0..N {
            let c0 = self.previous_law.branches()[i].modulus_pa();
            let c1 = self.law.branches()[i].modulus_pa();
            if c0 != c1 {
                for j in 0..6 {
                    transported.backstress_pa[i][j] = (state.backstress_pa[i][j] / c0) * c1;
                }
            }
        }
        self.law.validate_history(&[0.0; 6], &transported)?;
        Ok(transported)
    }

    fn mechanical_strain(&self, strain: &Voigt) -> Result<Voigt, ChabocheError> {
        self.law
            .backbone()
            .backbone()
            .validate_overstress_history(strain, &self.law.initial_state().plastic)?;
        Ok(core::array::from_fn(|i| {
            strain[i] - if i < 3 { self.expansion_strain } else { 0.0 }
        }))
    }

    /// Exact native tangent at frozen incoming history, including the change in
    /// flow direction. Temperature and its branch transport remain fixed here.
    pub fn response(
        &self,
        strain: &Voigt,
        state: &ChabocheState<N>,
    ) -> Result<ChabocheResponse<N>, ChabocheError> {
        let mechanical = self.mechanical_strain(strain)?;
        let transported = self.transport_history(state)?;
        let mut response = match &self.rate {
            Some(rate) => rate.response(&mechanical, &transported),
            None => self.law.response(&mechanical, &transported),
        }?;
        response.strain = *strain;
        Ok(response)
    }
}

impl<const N: usize> SmallStrainLaw for ThermalChaboche<N> {
    type State = ChabocheState<N>;
    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        match (
            self.mechanical_strain(strain),
            self.transport_history(state),
        ) {
            (Ok(mechanical), Ok(transported)) => self.rate.as_ref().map_or_else(
                || self.law.stress(&mechanical, &transported),
                |rate| rate.stress(&mechanical, &transported),
            ),
            _ => [f64::NAN; 6],
        }
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        match (
            self.mechanical_strain(strain),
            self.transport_history(state),
        ) {
            (Ok(mechanical), Ok(transported)) => self.rate.as_ref().map_or_else(
                || self.law.update_state(&mechanical, &transported),
                |rate| rate.update_state(&mechanical, &transported),
            ),
            _ => ChabocheState {
                plastic: J2State {
                    plastic_strain: [f64::NAN; 6],
                    alpha: f64::NAN,
                },
                backstress_pa: [[f64::NAN; 6]; N],
            },
        }
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            failure_envelope: "small-strain thermal Chaboche with affine positive C, fixed gamma and optional explicit-time overstress",
            ..self
                .rate
                .as_ref()
                .map_or_else(|| self.law.admissibility(), SmallStrainLaw::admissibility)
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.j2-thermal-chaboche",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain; source-pinned affine branch moduli and fixed recovery".into(),
                "fixed strain-like branch history under temperature change".into(),
                "native implicit return and generally nonsymmetric consistent tangent".into(),
                "branch entropy, coefficient storage and physical recovery retained separately"
                    .into(),
                if self.rate.is_some() {
                    "joint source-defined overstress at explicit physical dt"
                } else {
                    "rate-independent response"
                }
                .into(),
            ],
            self.law.card().validity,
            vec![
                "no static recovery, sub-yield creep or finite-deformation plasticity".into(),
                "no experimental thermal-cyclic calibration or fatigue-life claim".into(),
            ],
            0.05,
        )
    }
}
