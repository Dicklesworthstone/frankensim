//! Source-bound thermal histories on the original layered plane-stress section.
//!
//! One prepared interval freezes each layer's source, old/new caloric enthalpy
//! and physical duration. Its endpoint laws use the original return maps and
//! transverse reduction. Auditing replays that same interval from the accepted
//! OLD histories, including every native Chaboche branch or creep memory.
//!
//! Temperature is uniform within each declared layer. The original thickness
//! quadrature integrates separate physical heat, reversible caloric exchange,
//! coefficient-dependent free energy, mechanical entropy and numerical loss.
//! This owner neither advances heat transport nor infers heat from a residual.
//! Reference density/geometry and the source family's solid, small-strain
//! applicability remain unchanged; refining layers resolves thermal gradients.

use super::{
    LayeredPlaneStressSection, LayeredPlaneStressState, PlaneStressLayer, SectionBudget,
    SectionError, ThicknessQuadrature,
};
use crate::SmallStrainLaw;
use crate::plastic::energy::EnergyTolerances;
use crate::plastic::reduction::ReductionConfig;
use core::marker::PhantomData;

/// Native plasticity and creep sources sharing the original section owner.
pub mod source;
use source::{ThermalSectionSourceError, ThermalSectionSourceFamily};

/// Exact integrated native account for one physical layer, or the complete
/// section. All energies are per reference midsurface area, in J/m².
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ThermoplasticLayerEnergy {
    /// Native yield, creep, viscous and recovery heat, with no numerical loss.
    pub physical_heat_j_m2: f64,
    /// Physical creep work, already included once in `physical_heat_j_m2`.
    pub creep_heat_j_m2: f64,
    /// Signed reversible exchange into the original stress-free caloric chart.
    pub reversible_heat_j_m2: f64,
    /// Actual two-temperature elastic/hardening free-energy change, INCLUDING
    /// `coefficient_storage_change_j_m2`. It excludes mechanical entropy.
    pub stored_change_j_m2: f64,
    /// Frozen-history free-energy change caused by temperature-dependent moduli.
    pub coefficient_storage_change_j_m2: f64,
    /// Change in T times mechanical entropy; add to free-energy change for U.
    pub thermoelastic_storage_change_j_m2: f64,
    /// Original constitutive integration loss; excluded from caloric heating.
    pub numerical_loss_j_m2: f64,
    /// Native endpoint stress contracted with the TOTAL tensor strain change.
    pub mechanical_work_j_m2: f64,
    /// Work-conjugate membrane/bending resultants times generalized strain change.
    pub generalized_work_j_m2: f64,
    /// Signed transverse constraint-residual work, never another heat source.
    pub transverse_residual_work_j_m2: f64,
    /// Original membrane forces [N/m] and bending moments [N].
    pub resultants: [f64; 6],
    /// Generalized plus transverse work minus ΔPsi, Δ(Ts), caloric heat and
    /// numerical loss. A tolerance bounds algebraic error, not model accuracy.
    pub balance_residual_j_m2: f64,
}

impl ThermoplasticLayerEnergy {
    /// The signed exchange supplied once to the original caloric owner.
    pub fn caloric_heat_j_m2(&self) -> f64 {
        self.physical_heat_j_m2 + self.reversible_heat_j_m2
    }

    fn finish(
        &mut self,
        strain_change: &[f64; 6],
        layer: Option<usize>,
        tolerance_j_m2: f64,
    ) -> Result<(), ThermoplasticSectionError> {
        self.generalized_work_j_m2 = self
            .resultants
            .iter()
            .zip(strain_change)
            .map(|(stress, strain)| stress * strain)
            .sum();
        let work = self.generalized_work_j_m2 + self.transverse_residual_work_j_m2;
        self.balance_residual_j_m2 = work
            - self.stored_change_j_m2
            - self.thermoelastic_storage_change_j_m2
            - self.caloric_heat_j_m2()
            - self.numerical_loss_j_m2;
        if [
            self.physical_heat_j_m2,
            self.creep_heat_j_m2,
            self.numerical_loss_j_m2,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
            || self.creep_heat_j_m2 > self.physical_heat_j_m2
            || [
                self.reversible_heat_j_m2,
                self.stored_change_j_m2,
                self.coefficient_storage_change_j_m2,
                self.thermoelastic_storage_change_j_m2,
            ]
            .iter()
            .chain(self.resultants.iter())
            .any(|x| !x.is_finite())
        {
            return Err(ThermoplasticSectionError::Invalid(
                "thermal section accounts must be finite and physical/numerical losses nonnegative",
            ));
        }
        for residual in [work - self.mechanical_work_j_m2, self.balance_residual_j_m2] {
            if !residual.is_finite() || residual.abs() > tolerance_j_m2 {
                return Err(ThermoplasticSectionError::Work {
                    layer,
                    residual_j_m2: residual,
                    tolerance_j_m2,
                });
            }
        }
        Ok(())
    }
}

/// Original layer order is preserved for conservative transfer to caloric cells.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermoplasticSectionAudit {
    /// One account per original physical layer, in the original order.
    pub layers: Vec<ThermoplasticLayerEnergy>,
    /// Accumulated in the original point order, including reconstructed stress
    /// resultants for comparison with the native equilibrium endpoint.
    pub total: ThermoplasticLayerEnergy,
}

/// Native source/state refusal, geometric mismatch or independent work gate.
#[derive(Debug)]
pub enum ThermoplasticSectionError {
    /// Invalid coupling controls or unrepresentable arithmetic.
    Invalid(&'static str),
    /// Original geometry, state, budget or transverse-reduction refusal.
    Section(SectionError),
    /// Original source-family refusal, with its physical layer and point.
    Source {
        /// Original physical layer index.
        layer: usize,
        /// Thickness point, or none during initial source resolution.
        point: Option<usize>,
        /// Unmodified native source or energy error.
        source: ThermalSectionSourceError,
    },
    /// An alleged endpoint does not equal the original native update.
    State {
        /// Original thickness point index.
        point: usize,
        /// Exact history or transverse-stress mismatch.
        what: &'static str,
    },
    /// Independently bounded work-conjugacy or thermodynamic closure failure.
    Work {
        /// Original layer, or none for the integrated section.
        layer: Option<usize>,
        /// Signed defect per original reference area.
        residual_j_m2: f64,
        /// Caller-declared positive absolute allowance.
        tolerance_j_m2: f64,
    },
}
impl core::fmt::Display for ThermoplasticSectionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermoplastic section refused: {self:?}")
    }
}
impl std::error::Error for ThermoplasticSectionError {}
impl From<SectionError> for ThermoplasticSectionError {
    fn from(value: SectionError) -> Self {
        Self::Section(value)
    }
}

/// Immutable layer sources and exactly one physical constitutive interval.
/// The source-bearing `PlaneStressLayer::law` is a native source FAMILY here;
/// `section()` contains its resolved original mechanical laws.
#[derive(Debug)]
pub struct ThermoplasticSectionStep<'s, 'a, S: ThermalSectionSourceFamily<'a>> {
    layers: &'s [PlaneStressLayer<S>],
    old_h_j_kg: &'s [f64],
    new_h_j_kg: &'s [f64],
    duration_s: f64,
    section: LayeredPlaneStressSection<S::Material>,
    temperatures_k: Vec<f64>,
    source_lifetime: PhantomData<&'a ()>,
}

impl<'s, 'a, S: ThermalSectionSourceFamily<'a>> ThermoplasticSectionStep<'s, 'a, S> {
    /// Resolve once per physical layer. Geometry, quadrature and state ordering
    /// are admitted by the original section owner. No state is advanced here.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        layers: &'s [PlaneStressLayer<S>],
        old_h_j_kg: &'s [f64],
        new_h_j_kg: &'s [f64],
        duration_s: f64,
        quadrature: ThicknessQuadrature,
        reduction: ReductionConfig,
        budget: SectionBudget,
    ) -> Result<Self, ThermoplasticSectionError> {
        if layers.is_empty()
            || old_h_j_kg.len() != layers.len()
            || new_h_j_kg.len() != layers.len()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || old_h_j_kg.iter().chain(new_h_j_kg).any(|x| !x.is_finite())
        {
            return Err(ThermoplasticSectionError::Invalid(
                "one finite old/new enthalpy per layer and positive physical duration required",
            ));
        }
        // Bound all derived allocations before material resolution. The native
        // section subsequently owns exact geometry and quadrature admission.
        if layers.len() > budget.maximum_layers {
            return Err(SectionError::Budget {
                quantity: "layers",
                requested: layers.len(),
                limit: budget.maximum_layers,
            }
            .into());
        }
        let mut points = 0usize;
        for layer in layers {
            points = layer
                .subdivisions
                .checked_mul(quadrature.points_per_subdivision())
                .and_then(|n| points.checked_add(n))
                .ok_or(ThermoplasticSectionError::Invalid(
                    "thickness point count overflow",
                ))?;
            if points > budget.maximum_quadrature_points {
                return Err(SectionError::Budget {
                    quantity: "quadrature points",
                    requested: points,
                    limit: budget.maximum_quadrature_points,
                }
                .into());
            }
        }
        let mut resolved = Vec::with_capacity(layers.len());
        let mut temperatures_k = Vec::with_capacity(layers.len());
        for (index, layer) in layers.iter().enumerate() {
            // Each source admits both original caloric endpoints before
            // preparing this physical interval's native constitutive law.
            let step = layer
                .law
                .resolve_section_step(old_h_j_kg[index], new_h_j_kg[index], duration_s)
                .map_err(|source| ThermoplasticSectionError::Source {
                    layer: index,
                    point: None,
                    source,
                })?;
            temperatures_k.push(step.temperature_k);
            resolved.push(PlaneStressLayer {
                bottom_z_m: layer.bottom_z_m,
                top_z_m: layer.top_z_m,
                subdivisions: layer.subdivisions,
                law: step.law,
            });
        }
        let section = LayeredPlaneStressSection::new(resolved, quadrature, reduction, budget)?;
        Ok(Self {
            layers,
            old_h_j_kg,
            new_h_j_kg,
            duration_s,
            section,
            temperatures_k,
            source_lifetime: PhantomData,
        })
    }

    /// Original native plane-stress section resolved for this exact interval.
    pub const fn section(&self) -> &LayeredPlaneStressSection<S::Material> {
        &self.section
    }
    /// Actual endpoint temperatures in the original physical layer order.
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    /// The single physical interval used by every native timed return.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Validate and audit an exact accepted section increment. Both generalized
    /// and full tensor strains must belong to the original section geometry.
    /// The returned native histories must match in full, including all branches.
    /// Consumers compare `total.resultants` with their solved stress resultants.
    pub fn audit_energy(
        &self,
        old: &LayeredPlaneStressState<<S::Material as SmallStrainLaw>::State>,
        new: &LayeredPlaneStressState<<S::Material as SmallStrainLaw>::State>,
        material: EnergyTolerances,
        work_tolerance_j_m2: f64,
    ) -> Result<ThermoplasticSectionAudit, ThermoplasticSectionError> {
        self.section.validate_state(old)?;
        self.section.validate_state(new)?;
        if [material.energy_j_m3, material.yield_pa, work_tolerance_j_m2]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(ThermoplasticSectionError::Invalid(
                "positive finite native material and section work allowances required",
            ));
        }
        let mut audit = ThermoplasticSectionAudit {
            layers: vec![ThermoplasticLayerEnergy::default(); self.layers.len()],
            total: ThermoplasticLayerEnergy::default(),
        };
        for (index, point) in self.section.points().iter().enumerate() {
            let before = &old.points[index];
            let after = &new.points[index];
            let layer = point.layer;
            let source = &self.layers[layer].law;
            let native_error = |source| ThermoplasticSectionError::Source {
                layer,
                point: Some(index),
                source,
            };
            let increment = source
                .section_increment(
                    self.old_h_j_kg[layer],
                    self.new_h_j_kg[layer],
                    &before.strain,
                    &before.history,
                    &after.strain,
                    material,
                    self.duration_s,
                )
                .map_err(native_error)?;
            if increment.state != after.history {
                return Err(ThermoplasticSectionError::State {
                    point: index,
                    what: "accepted complete history differs from the original thermal increment",
                });
            }
            let stress = increment.stress_pa;
            if stress.iter().any(|x| !x.is_finite())
                || [stress[2], stress[4], stress[5]]
                    .iter()
                    .any(|x| x.abs() > self.section.reduction_config().stress_tolerance_pa)
            {
                return Err(ThermoplasticSectionError::State {
                    point: index,
                    what: "accepted native stress exceeds the transverse equilibrium allowance",
                });
            }
            let transverse = stress[2] * (after.strain[2] - before.strain[2])
                + 2.0 * stress[4] * (after.strain[4] - before.strain[4])
                + 2.0 * stress[5] * (after.strain[5] - before.strain[5]);
            for account in [&mut audit.layers[layer], &mut audit.total] {
                let w = point.weight_m;
                account.physical_heat_j_m2 += w * increment.physical_heat_j_m3;
                account.creep_heat_j_m2 += w * increment.creep_heat_j_m3;
                account.reversible_heat_j_m2 += w * increment.reversible_heat_j_m3;
                account.stored_change_j_m2 += w * increment.stored_change_j_m3;
                account.coefficient_storage_change_j_m2 +=
                    w * increment.coefficient_storage_change_j_m3;
                account.thermoelastic_storage_change_j_m2 +=
                    w * increment.thermoelastic_storage_change_j_m3;
                account.numerical_loss_j_m2 += w * increment.numerical_loss_j_m3;
                account.mechanical_work_j_m2 += w * increment.mechanical_work_j_m3;
                account.transverse_residual_work_j_m2 += w * transverse;
                for (i, component) in [0, 1, 3].into_iter().enumerate() {
                    account.resultants[i] += w * stress[component];
                    account.resultants[i + 3] += (w * point.z_m) * stress[component];
                }
            }
        }
        let strain_change =
            core::array::from_fn(|i| new.generalized_strain[i] - old.generalized_strain[i]);
        for (layer, account) in audit.layers.iter_mut().enumerate() {
            account.finish(&strain_change, Some(layer), work_tolerance_j_m2)?;
        }
        audit
            .total
            .finish(&strain_change, None, work_tolerance_j_m2)?;
        Ok(audit)
    }
}
