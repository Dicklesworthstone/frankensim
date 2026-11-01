//! Source-bound sharp melting in the original spatial enthalpy transport.
//!
//! Each phase's k(T) is queried only on its own source interval. Coexistence
//! uses a caller-selected parallel or series volume-fraction closure at the
//! original melting temperature. These are explicit effective-medium choices,
//! not a microstructure inferred from Cp, mass fraction or a material name.
//! Nodal k(h) is interpolated with P1 basis functions; its element mean and
//! exact enthalpy tangent feed the existing conductivity assembly and Newton
//! solver. No apparent heat capacity or separate latent-energy account exists.
//!
//! Reference density is resolved once at an explicit single-phase state. The
//! changing equilibrium density is used only to convert phase mass fractions
//! to volume fractions; it never changes reference masses or mesh geometry.
//! This stationary isotropic model does not include shrinkage, convection,
//! pressure work, phase kinetics or source-uncertainty propagation.

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_matdb::{MaterialAnswer, MaterialCard, PropertyKey, QueryPoint};
use fs_material::phase::SolidLiquidPhase;
use fs_material::phase::melting::{MeltingSourceQuery, ResolvedMeltingCurve};

use super::super::{EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError};
use super::MaterialEnthalpyError;
use crate::{
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, ScalarField,
    ThermalBoundary,
};

/// Declared isotropic effective conductivity of coexisting phases.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoexistenceConductivity {
    /// Arithmetic volume mixture, k=(1-phi)*k_s+phi*k_l.
    Parallel,
    /// Harmonic volume mixture, 1/k=(1-phi)/k_s+phi/k_l.
    Series,
}

/// Source selection for caloric storage, phase transport and reference mass.
#[derive(Debug, Clone, Copy)]
pub struct MaterialMeltingSource<'a> {
    /// Original solid/liquid caloric sources and one sharp fusion plateau.
    pub caloric: MeltingSourceQuery<'a>,
    /// Solid then liquid isotropic conductivity keys from the SAME card.
    /// Source conditions are the corresponding caloric branch's exact queries.
    pub conductivity: [&'a PropertyKey; 2],
    /// Explicit solid or liquid branch for the reference-density query.
    /// Coexistence has no unique density at a prescribed temperature and refuses.
    pub reference_phase: SolidLiquidPhase,
    /// All source coordinates match the selected branch; T lies in its support.
    pub reference: &'a QueryPoint,
    /// Specific enthalpy at the solid branch's lower endpoint, J/kg.
    pub lower_specific_enthalpy_j_kg: f64,
    /// Sensible-chart interpolation allowance, not physical source uncertainty.
    pub maximum_interpolation_error_j_kg: f64,
    pub coexistence: CoexistenceConductivity,
}

/// One conductivity value and its exact one-sided admitted-chart tangent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeltingConductivitySample {
    pub conductivity_w_m_k: f64,
    /// d(k)/d(h), W/(m K) per J/kg, including the latent plateau.
    pub derivative_per_j_kg: f64,
    /// Liquid volume fraction, distinct from the chart's liquid mass fraction.
    pub liquid_volume_fraction: f64,
}

/// One retained material source used by the existing spatial transport owner.
#[derive(Debug, Clone)]
pub struct MaterialMeltingTransport {
    card_identity: ContentHash,
    caloric: ResolvedMeltingCurve,
    conductivity: [ConductivityTable; 2],
    transition_conductivity: [f64; 2],
    // Densities divided by their maximum keep coexistence arithmetic scaled.
    transition_density_ratios: [f64; 2],
    reference_query: QueryPoint,
    reference_density: MaterialAnswer,
    coexistence: CoexistenceConductivity,
    numerical_basis: ConductivityModel,
    identity: ContentHash,
}

impl MaterialMeltingTransport {
    /// Resolve both source branches once, without extending either through the
    /// other phase. Cp, rho and k require explicit source temperature support.
    /// The original caloric owner admits the transition and all fixed context.
    pub fn from_isotropic_card(
        card: &MaterialCard,
        source: MaterialMeltingSource<'_>,
    ) -> Result<Self, MaterialEnthalpyError> {
        let caloric = ResolvedMeltingCurve::try_new(
            card,
            source.caloric,
            source.lower_specific_enthalpy_j_kg,
            source.maximum_interpolation_error_j_kg,
        )
        .map_err(MaterialEnthalpyError::Caloric)?;
        let branches = [source.caloric.solid, source.caloric.liquid];
        let table = |index: usize| {
            let branch = branches[index];
            ConductivityTable::from_claims_at_query_points(
                card.claims(),
                source.conductivity[index],
                branch.temperature_axis,
                &[branch.lower.clone(), branch.upper.clone()],
                branch.selection,
            )
            .map_err(MaterialEnthalpyError::Conductivity)
        };
        let conductivity = [table(0)?, table(1)?];
        let axis = source.caloric.solid.temperature_axis;
        let caloric_receipts = caloric.receipts().iter().filter(|receipt| {
            receipt.property != source.caloric.melting_point.name()
                && receipt.property != source.caloric.latent_heat_fusion.name()
        });
        for receipt in
            caloric_receipts.chain(conductivity.iter().flat_map(|table| table.receipts()))
        {
            if !card
                .claims()
                .claim(receipt.selected)
                .is_some_and(|claim| claim.validity.bound(axis).is_some())
            {
                return Err(invalid(
                    "each phase's Cp, density and conductivity require explicit temperature support",
                ));
            }
        }
        let reference_branch = match source.reference_phase {
            SolidLiquidPhase::Solid => branches[0],
            SolidLiquidPhase::Liquid => branches[1],
            SolidLiquidPhase::SolidLiquid => {
                return Err(invalid(
                    "reference density requires an explicit single phase",
                ));
            }
        };
        let reference_k =
            source.reference.axes().get(axis).copied().ok_or_else(|| {
                invalid("reference query lacks the source temperature coordinate")
            })?;
        let expected = reference_branch
            .lower
            .clone()
            .with(axis, reference_k)
            .map_err(MaterialEnthalpyError::Reference)?;
        let lower_k = reference_branch.lower.axes()[axis];
        let upper_k = reference_branch.upper.axes()[axis];
        if !reference_k.is_finite()
            || reference_k < lower_k
            || reference_k > upper_k
            || source.reference.axes() != expected.axes()
            || source.reference.axis_quantities() != expected.axis_quantities()
        {
            return Err(invalid(
                "reference density must remain in its phase's source interval and exact fixed context",
            ));
        }
        let density_at = |point| {
            card.claims()
                .query_typed(reference_branch.density, point, reference_branch.selection)
                .map_err(MaterialEnthalpyError::Reference)
        };
        let reference_density = density_at(source.reference)?;
        let branch_density = density_at(reference_branch.lower)?;
        if !reference_density.evidence.value.value.is_finite()
            || reference_density.evidence.value.value <= 0.0
            || reference_density.receipt.selected != branch_density.receipt.selected
        {
            return Err(invalid(
                "reference density must use the original branch's positive density source",
            ));
        }
        let transition_temperature = caloric.melting_temperature_k();
        let transition_conductivity = [
            conductivity[0].eval(transition_temperature),
            conductivity[1].eval(transition_temperature),
        ];
        let [solid_k, liquid_k] = transition_conductivity;
        let transition_conductivity = [
            solid_k.map_err(MaterialEnthalpyError::Conductivity)?,
            liquid_k.map_err(MaterialEnthalpyError::Conductivity)?,
        ];
        let [solid_h, liquid_h] = caloric.plateau_enthalpy_bounds_j_kg();
        let density = |h| {
            caloric
                .curve()
                .state_at_specific_enthalpy(h)
                .map(|state| state.bulk_density_kg_m3())
                .map_err(MaterialEnthalpyError::Caloric)
        };
        let density = [density(solid_h)?, density(liquid_h)?];
        let scale = density[0].max(density[1]);
        let transition_density_ratios = [density[0] / scale, density[1] / scale];
        if transition_density_ratios
            .iter()
            .any(|&ratio| !ratio.is_finite() || ratio <= 0.0)
        {
            return Err(invalid("phase-density ratio is not representable"));
        }
        let mut identity =
            DomainHasher::new("org.frankensim.fs-conduction.sourced-melting-transport.v1");
        identity.update(card.content_hash().as_bytes());
        identity.update(caloric.curve().identity().as_bytes());
        identity.update(reference_density.receipt.content_hash().as_bytes());
        identity.update(&[match source.coexistence {
            CoexistenceConductivity::Parallel => 0,
            CoexistenceConductivity::Series => 1,
        }]);
        for table in &conductivity {
            for receipt in table.receipts() {
                identity.update(receipt.content_hash().as_bytes());
            }
        }
        Ok(Self {
            card_identity: card.content_hash(),
            caloric,
            conductivity,
            transition_conductivity,
            transition_density_ratios,
            reference_query: source.reference.clone(),
            reference_density,
            coexistence: source.coexistence,
            // Only a numerical basis: every physical W/(m K) comes from the
            // retained nodal source law above. It is not a material claim.
            numerical_basis: ConductivityModel::isotropic_declared(1.0)
                .map_err(MaterialEnthalpyError::Conductivity)?,
            identity: identity.finalize(),
        })
    }

    #[must_use]
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }

    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Original source caloric chart, including every Cp/rho/fusion receipt.
    #[must_use]
    pub const fn caloric(&self) -> &ResolvedMeltingCurve {
        &self.caloric
    }

    /// Solid then liquid k(T) tables, retaining the original source receipts.
    #[must_use]
    pub const fn conductivity_sources(&self) -> &[ConductivityTable; 2] {
        &self.conductivity
    }

    #[must_use]
    pub const fn coexistence(&self) -> CoexistenceConductivity {
        self.coexistence
    }

    #[must_use]
    pub const fn reference_query(&self) -> &QueryPoint {
        &self.reference_query
    }

    #[must_use]
    pub const fn reference_density(&self) -> &MaterialAnswer {
        &self.reference_density
    }

    #[must_use]
    pub fn reference_density_kg_m3(&self) -> f64 {
        self.reference_density.evidence.value.value
    }

    /// Evaluate the original source branches and explicit coexistence law.
    /// At the two fusion endpoints the right enthalpy derivative is selected,
    /// matching the original caloric chart's semismooth Newton convention.
    /// A classical two-sided derivative is not asserted at a source/chart kink.
    pub fn conductivity_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<MeltingConductivitySample, MaterialEnthalpyError> {
        let curve = self.caloric.curve();
        let state = curve
            .state_at_specific_enthalpy(h)
            .map_err(MaterialEnthalpyError::Caloric)?;
        let [solid_h, liquid_h] = self.caloric.plateau_enthalpy_bounds_j_kg();
        let sample = if h < solid_h || h >= liquid_h {
            let phase = usize::from(h >= liquid_h);
            let table = &self.conductivity[phase];
            let temperature = state.temperature_k();
            let slope = curve
                .temperature_derivative_at_specific_enthalpy(h)
                .map_err(MaterialEnthalpyError::Caloric)?;
            let k_slope = table
                .derivative(temperature)
                .map_err(MaterialEnthalpyError::Conductivity)?;
            let derivative = k_slope * slope;
            if k_slope != 0.0 && slope != 0.0 && derivative == 0.0 {
                return Err(invalid("nonzero sensible conductivity tangent underflowed"));
            }
            MeltingConductivitySample {
                conductivity_w_m_k: table
                    .eval(temperature)
                    .map_err(MaterialEnthalpyError::Conductivity)?,
                derivative_per_j_kg: derivative,
                liquid_volume_fraction: phase as f64,
            }
        } else {
            let f = state.liquid_mass_fraction();
            let [solid_rho, liquid_rho] = self.transition_density_ratios;
            let denominator = (1.0 - f) * liquid_rho + f * solid_rho;
            let phi = f * solid_rho / denominator;
            let mass_slope = curve
                .liquid_mass_fraction_derivative_at_specific_enthalpy(h)
                .map_err(MaterialEnthalpyError::Caloric)?;
            let volume_slope = (solid_rho / denominator) * (liquid_rho / denominator) * mass_slope;
            let [ks, kl] = self.transition_conductivity;
            let (k, derivative_per_volume) = match self.coexistence {
                CoexistenceConductivity::Parallel => ((1.0 - phi) * ks + phi * kl, kl - ks),
                CoexistenceConductivity::Series => {
                    let scale = ks.min(kl);
                    let rs = scale / ks;
                    let rl = scale / kl;
                    let denominator = (1.0 - phi) * rs + phi * rl;
                    let k = scale / denominator;
                    (k, k * ((rs - rl) / denominator))
                }
            };
            let derivative = derivative_per_volume * volume_slope;
            if !volume_slope.is_finite()
                || volume_slope <= 0.0
                || (f > 0.0 && f < 1.0 && (phi <= 0.0 || phi >= 1.0))
                || (ks != kl && derivative == 0.0)
            {
                return Err(invalid(
                    "coexistence volume fraction or conductivity tangent is not representable",
                ));
            }
            MeltingConductivitySample {
                conductivity_w_m_k: k,
                derivative_per_j_kg: derivative,
                liquid_volume_fraction: phi,
            }
        };
        if !sample.conductivity_w_m_k.is_finite()
            || sample.conductivity_w_m_k <= 0.0
            || !sample.derivative_per_j_kg.is_finite()
            || !sample.liquid_volume_fraction.is_finite()
            || !(0.0..=1.0).contains(&sample.liquid_volume_fraction)
        {
            return Err(invalid(
                "source conductivity value or tangent is not finite and admissible",
            ));
        }
        Ok(sample)
    }

    /// Bind the retained caloric and conductivity law to original P1 storage.
    /// The ordinary forward, adaptive, radiation and endpoint-adjoint APIs then
    /// use this same law. Spatial budgets precede reference-mass allocation.
    pub fn prepare_uniform<'m, 's>(
        &'s self,
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        budget: EnthalpyBudget,
    ) -> Result<EnthalpyBackwardEuler<'m, 's>, EnthalpyError> {
        let mut storage = EnthalpyBackwardEuler::uniform(
            cx,
            mesh,
            self.caloric.curve(),
            self.reference_density_kg_m3(),
            budget,
        )?;
        storage.sourced_melting = Some(self);
        Ok(storage)
    }

    /// Construct the existing problem with this law's retained numerical basis.
    /// Changing its material tensor or adding element overrides refuses when
    /// advancing: neither may silently multiply the physical source k twice.
    #[must_use]
    pub fn problem<'a>(
        &'a self,
        mesh: &'a ConductionMesh,
        boundary: &'a ThermalBoundary,
        source: &'a ScalarField,
    ) -> ConductionProblem<'a> {
        ConductionProblem {
            mesh,
            boundary,
            material: &self.numerical_basis,
            element_materials: None,
            source,
        }
    }

    pub(in crate::transient::enthalpy) fn matches_numerical_basis(
        &self,
        problem: ConductionProblem<'_>,
    ) -> bool {
        std::ptr::eq(problem.material, &self.numerical_basis) && problem.element_materials.is_none()
    }

    pub(in crate::transient::enthalpy) fn sample_for_vertex(
        &self,
        h: f64,
        vertex: usize,
    ) -> Result<MeltingConductivitySample, EnthalpyError> {
        self.conductivity_at_specific_enthalpy(h)
            .map_err(|error| match error {
                MaterialEnthalpyError::Caloric(source) => EnthalpyError::Phase { vertex, source },
                MaterialEnthalpyError::Conductivity(source) => EnthalpyError::Conduction(source),
                MaterialEnthalpyError::Invalid(what) => EnthalpyError::InvalidInput(what),
                MaterialEnthalpyError::Reference(_) => {
                    EnthalpyError::InvalidInput("unexpected live reference-density query")
                }
            })
    }

    // The original chart guard already refuses its endpoints and slope kinks.
    // Independently sampled k(T) can introduce additional corners between Cp
    // knots, so those must also refuse a classical endpoint adjoint.
    pub(in crate::transient::enthalpy) fn conductivity_kink_at(
        &self,
        h: f64,
        vertex: usize,
    ) -> Result<bool, EnthalpyError> {
        self.sample_for_vertex(h, vertex)?;
        let [solid_h, liquid_h] = self.caloric.plateau_enthalpy_bounds_j_kg();
        if h >= solid_h && h <= liquid_h {
            return Ok(false);
        }
        let phase = usize::from(h > liquid_h);
        let t = self
            .caloric
            .curve()
            .temperature_at_specific_enthalpy(h)
            .map_err(|source| EnthalpyError::Phase { vertex, source })?;
        let knots = self.conductivity[phase].knots();
        if let Ok(index) = knots.binary_search_by(|knot| knot.0.total_cmp(&t)) {
            if index > 0 && index + 1 < knots.len() {
                let left =
                    (knots[index].1 - knots[index - 1].1) / (knots[index].0 - knots[index - 1].0);
                let right =
                    (knots[index + 1].1 - knots[index].1) / (knots[index + 1].0 - knots[index].0);
                return Ok(left != right);
            }
        }
        Ok(false)
    }
}

fn invalid(what: &'static str) -> MaterialEnthalpyError {
    MaterialEnthalpyError::Invalid(what)
}
