//! Source-bound anisotropic Maxwell branches in the material's original axes.
//!
//! Equilibrium elasticity uses the canonical elastic properties. Each branch
//! explicitly names all 21 entries of a symmetric Mandel matrix, including
//! zeros; no direction, shear convention, rate or missing coefficient is
//! inferred from a material name. Plane-stress elimination belongs to the
//! consuming structural owner and must act on the complete hereditary law.

use super::{TensorProny, TensorPronyBranch};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ORTHOTROPIC_POISSON_RATIO_PROPERTIES,
    ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES, ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES,
    ResolvedMaterialStatePoint, ScalarAdmissibility as Admit,
    ScalarPropertyRequirement as Requirement, resolve_material_state_point,
};
use crate::{IsotropicElastic, OrthotropicElastic};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

/// Canonical equilibrium elastic data, before any spatial rotation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TensorElasticity {
    /// `young_modulus` and `poisson_ratio`.
    Isotropic,
    /// The complete three Young moduli, three Poisson ratios and three shears.
    Orthotropic,
}

/// One source branch in Mandel order (11,22,33,sqrt(2)12,sqrt(2)23,sqrt(2)31).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TensorBranchProperties {
    /// Upper triangle in row order: C11,C12,...,C16,C22,...,C66, all in Pa.
    pub stiffness_upper: [String; 21],
    /// Original positive relaxation time in seconds.
    pub relaxation_time: String,
}

/// Reusable selection of immutable source claims for one complete 3-D law.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TensorPronyProperties {
    /// Equilibrium elastic property family, shared with the structural chart.
    pub elasticity: TensorElasticity,
    /// Positive dimensionless Frobenius small-strain limit from the source.
    pub strain_limit: String,
    /// Explicit branch matrices and times; at most 64 source branches.
    pub branches: Vec<TensorBranchProperties>,
}

/// A missing, incompatible or inadmissible source cannot become a material law.
#[derive(Debug)]
pub enum TensorSourceError {
    /// Original material requirement or validity refusal.
    StatePoint(MaterialStatePointError),
    /// Complete source/property or constitutive admission failure.
    Invalid(String),
}
impl core::fmt::Display for TensorSourceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::StatePoint(error) => write!(f, "tensor material source: {error}"),
            Self::Invalid(reason) => write!(f, "tensor material source: {reason}"),
        }
    }
}
impl std::error::Error for TensorSourceError {}
impl From<MaterialStatePointError> for TensorSourceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::StatePoint(error)
    }
}

type Result<T> = core::result::Result<T, TensorSourceError>;

/// Reference inertia and the complete hereditary law from one immutable card.
///
/// Geometry and material orientation belong to the consuming specimen. Density,
/// equilibrium elasticity, every Maxwell tensor and every relaxation time share
/// one physical query and source selection. No equivalent isotropic material,
/// fitted modal damping or density-derived modulus is introduced.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedTensorMaterial {
    state: ResolvedMaterialStatePoint,
    properties: TensorPronyProperties,
    law: TensorProny,
    density_kg_m3: f64,
}

impl ResolvedTensorMaterial {
    /// Atomically resolve the source law and reference density. Frequency
    /// applicability of density also restricts the law's original source band;
    /// a sampled frequency-dependent density cannot become constant inertia.
    /// Temperature and geometry stay fixed for the resulting dynamics. This
    /// binding does not implement thermal expansion or a material transition.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        properties: TensorPronyProperties,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        let mut requirements = properties.requirements()?;
        if requirements.iter().any(|r| r.name() == "density") {
            return Err(TensorSourceError::Invalid(
                "density cannot also name an elastic, strain or relaxation property".into(),
            ));
        }
        requirements.push(Requirement::try_new(
            "density",
            Density::DIMS,
            Admit::StrictlyPositive,
        )?);
        let state = resolve_material_state_point(card, point, &requirements, selection)?;
        Self::from_state(state, properties)
    }

    fn from_state(
        state: ResolvedMaterialStatePoint,
        properties: TensorPronyProperties,
    ) -> Result<Self> {
        let density = state.property("density").expect("resolved source density");
        if density.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
            || !density.is_constant_at_fixed_temperature()
            || !density.value_si().is_finite()
            || density.value_si() <= 0.0
        {
            return Err(TensorSourceError::Invalid(
                "reference density must be a positive SI coefficient constant at the declared external state".into(),
            ));
        }
        let mut law = properties.bind(&state)?;
        if let Some((low, high)) = density.angular_frequency_band_rad_s()? {
            let band = law.band_rad_s();
            law = law
                .with_band((band.0.max(low), band.1.min(high)))
                .map_err(|e| TensorSourceError::Invalid(e.to_string()))?;
        }
        let density_kg_m3 = density.value_si();
        Ok(Self {
            state,
            properties,
            law,
            density_kg_m3,
        })
    }

    /// Re-evaluate the SAME complete source selection at another physical
    /// endpoint. This rebuilds tensors and relaxation times; it never transfers
    /// constitutive history, changes specimen mass, or supplies thermal work.
    /// Those operations belong to the consuming evolution owner.
    pub fn at_query_point(&self, card: &MaterialCard, point: &QueryPoint) -> Result<Self> {
        Self::from_state(
            self.state.at_query_point(card, point)?,
            self.properties.clone(),
        )
    }

    /// Admit the bounded thermal model in which ONLY Maxwell relaxation times
    /// depend on temperature. Stored elastic tensors, strain bounds and reference
    /// density must be validity-wide scalar constants, not two equal samples of
    /// a varying law. There is then no mechanical entropy or coefficient-work
    /// contribution: temperature changes the dissipative kinetics only.
    ///
    /// Every selected source retains its original claim and external context.
    /// One-dimensional linearly interpolated relaxation curves must cover the
    /// complete interval, and every interior knot must remain positive. This
    /// admission supplies no thermal integrator, expansion or phase transition.
    pub fn admit_thermal_relaxation(
        &self,
        card: &MaterialCard,
        temperature_axis: &str,
        interval_k: [f64; 2],
    ) -> Result<()> {
        use fs_matdb::{InterpolationPolicy, PropertyValue};
        use fs_qty::{
            Temperature,
            semantic::{QuantityKind, SemanticType, ValueForm},
        };

        let invalid = |s: &str| TensorSourceError::Invalid(s.into());
        let [low, high] = interval_k;
        let temperature = self
            .state
            .state_coordinate(temperature_axis)
            .ok_or_else(|| {
                invalid("thermal relaxation requires its original temperature coordinate")
            })?;
        let quantity = self.state.axis_quantities().get(temperature_axis);
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || !(low..=high).contains(&temperature)
            || quantity.is_some_and(|q| {
                *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS)
            })
            || quantity.is_none() && temperature_axis != "T"
        {
            return Err(invalid(
                "thermal relaxation needs a bounded absolute-temperature interval containing its original state",
            ));
        }
        for property in self.state.properties() {
            if !property.has_only_fixed_external_conditions(Some(temperature_axis), true) {
                return Err(invalid(
                    "thermal relaxation cannot freeze unowned material response or validity coordinates",
                ));
            }
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or_else(|| invalid("original thermal relaxation claim is absent"))?;
            let relaxation = self
                .properties
                .branches
                .iter()
                .any(|branch| branch.relaxation_time == property.requirement().name());
            match (&claim.value, claim.interpolation) {
                (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {}
                (
                    PropertyValue::Curve {
                        abscissa,
                        abscissa_dims,
                        knots,
                        ..
                    },
                    InterpolationPolicy::LinearInside,
                ) if relaxation
                    && abscissa == temperature_axis
                    && *abscissa_dims == Temperature::DIMS =>
                {
                    if knots.len() < 2
                        || knots.first().unwrap().0 > low
                        || knots.last().unwrap().0 < high
                        || knots.iter().any(|&(t, tau)| {
                            (low..=high).contains(&t) && (!tau.is_finite() || tau <= 0.0)
                        })
                    {
                        return Err(invalid(
                            "relaxation curve must cover the interval with positive interior knots",
                        ));
                    }
                }
                _ => {
                    return Err(invalid(
                        "thermal feedback permits temperature-dependent relaxation only; stored tensors and density require constant source laws",
                    ));
                }
            }
        }
        // Fixed-context interpolation and rectangular validity have no hidden
        // gap between these endpoints. The canonical owner also checks card
        // identity, exact selected claims, units and positive endpoint rates.
        for temperature in interval_k {
            let point = self
                .state
                .query()
                .with(temperature_axis, temperature)
                .map_err(|e| invalid(&e.to_string()))?;
            self.at_query_point(card, &point)?;
        }
        Ok(())
    }

    /// Original same-card query, property selections and receipts.
    pub const fn state(&self) -> &ResolvedMaterialStatePoint {
        &self.state
    }
    /// Original equilibrium family and explicit Maxwell property selectors.
    pub const fn properties(&self) -> &TensorPronyProperties {
        &self.properties
    }
    /// Complete law in its original material frame, including its source band.
    pub const fn law(&self) -> &TensorProny {
        &self.law
    }
    /// Fixed density in the reference geometry [kg/m³].
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Express every original tensor in another frame with the SAME proper
    /// material-to-target rotation. No material direction is inferred from
    /// element ordering. Branch rates, strain limits and source band survive.
    pub fn rotated_law(&self, material_to_target: [[f64; 3]; 3]) -> Result<TensorProny> {
        self.law
            .rotated(material_to_target)
            .map_err(|e| TensorSourceError::Invalid(e.to_string()))
    }
}

impl TensorPronyProperties {
    /// Requirements to add to the original card query. Repeated names resolve
    /// once, while incompatible reuse (for example a modulus as a time) refuses.
    pub fn requirements(&self) -> Result<Vec<Requirement>> {
        if self.branches.len() > 64 {
            return Err(TensorSourceError::Invalid(
                "more than 64 source branches".into(),
            ));
        }
        let mut requirements = Vec::<Requirement>::new();
        let mut add = |name: &str, dims, admit| -> Result<()> {
            if name.len() > 128 {
                return Err(TensorSourceError::Invalid(
                    "property name exceeds 128 bytes".into(),
                ));
            }
            let requirement = Requirement::try_new(name, dims, admit)?;
            if let Some(existing) = requirements.iter().find(|r| r.name() == name) {
                if existing.quantity() != requirement.quantity() {
                    return Err(TensorSourceError::Invalid(format!(
                        "{name}: incompatible property roles"
                    )));
                }
            } else {
                requirements.push(requirement);
            }
            Ok(())
        };
        match self.elasticity {
            TensorElasticity::Isotropic => {
                add("young_modulus", Pressure::DIMS, Admit::StrictlyPositive)?;
                add(
                    "poisson_ratio",
                    Dims::NONE,
                    Admit::OpenInterval {
                        lower: -1.0,
                        upper: 0.5,
                    },
                )?;
            }
            TensorElasticity::Orthotropic => {
                for key in ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES {
                    add(key, Pressure::DIMS, Admit::StrictlyPositive)?;
                }
                for key in ORTHOTROPIC_POISSON_RATIO_PROPERTIES {
                    add(key, Dims::NONE, Admit::Finite)?;
                }
                for key in ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES {
                    add(key, Pressure::DIMS, Admit::StrictlyPositive)?;
                }
            }
        }
        add(&self.strain_limit, Dims::NONE, Admit::StrictlyPositive)?;
        for branch in &self.branches {
            add(
                &branch.relaxation_time,
                Dims([0, 0, 1, 0, 0, 0]),
                Admit::StrictlyPositive,
            )?;
            for key in &branch.stiffness_upper {
                add(key, Pressure::DIMS, Admit::Finite)?;
            }
        }
        Ok(requirements)
    }

    /// Bind every selected coefficient at one original material state.
    /// The admitted frequency band is the intersection of ALL selected source
    /// domains, including equilibrium stiffness and the strain bound. Constant
    /// coefficients at fixed temperature are required: a sampled modulus at one
    /// frequency is not silently treated as an equilibrium spring.
    pub fn bind(&self, state: &ResolvedMaterialStatePoint) -> Result<TensorProny> {
        let requirements = self.requirements()?;
        let mut band = (0.0_f64, f64::MAX);
        for requirement in &requirements {
            let name = requirement.name();
            let property = state.property(name).ok_or_else(|| {
                TensorSourceError::Invalid(format!("{name}: selected source property is missing"))
            })?;
            if property.requirement().quantity() != requirement.quantity()
                || !property.value_si().is_finite()
                || !property.is_constant_at_fixed_temperature()
            {
                return Err(TensorSourceError::Invalid(format!(
                    "{name}: requires dimensioned finite coefficients constant at the original state"
                )));
            }
            if let Some((low, high)) = property.angular_frequency_band_rad_s()? {
                band.0 = band.0.max(low);
                band.1 = band.1.min(high);
            }
        }
        let value = |name: &str| -> f64 {
            // All selected property names and dimensions were checked above.
            state
                .property(name)
                .expect("checked tensor source property")
                .value_si()
        };
        let strain_limit = value(&self.strain_limit);
        if !strain_limit.is_finite() || strain_limit <= 0.0 {
            return Err(TensorSourceError::Invalid(
                "small-strain source bound must be positive".into(),
            ));
        }
        let equilibrium = match self.elasticity {
            TensorElasticity::Isotropic => {
                IsotropicElastic::new(value("young_modulus"), value("poisson_ratio"), strain_limit)
                    .map_err(|e| TensorSourceError::Invalid(e.to_string()))?
                    .stiffness()
            }
            TensorElasticity::Orthotropic => OrthotropicElastic::new(
                ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES.map(value),
                ORTHOTROPIC_POISSON_RATIO_PROPERTIES.map(value),
                ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES.map(value),
                strain_limit,
            )
            .map_err(|e| TensorSourceError::Invalid(e.to_string()))?
            .stiffness(),
        };
        let mut branches = Vec::with_capacity(self.branches.len());
        for source in &self.branches {
            let mut matrix = [[0.0; 6]; 6];
            let mut entry = 0;
            for row in 0..6 {
                for column in row..6 {
                    let coefficient = value(&source.stiffness_upper[entry]);
                    matrix[row][column] = coefficient;
                    matrix[column][row] = coefficient;
                    entry += 1;
                }
            }
            // The existing constitutive owner admits PSD branches, including
            // rank-deficient ones, with its original roundoff policy. This
            // binding adds no regularization or positive loss floor.
            branches.push(
                TensorPronyBranch::new(matrix, value(&source.relaxation_time))
                    .map_err(|e| TensorSourceError::Invalid(e.to_string()))?,
            );
        }
        TensorProny::from_tensor_shear(equilibrium, branches, strain_limit)
            .map_err(|e| TensorSourceError::Invalid(e.to_string()))?
            .with_band(band)
            .map_err(|e| TensorSourceError::Invalid(e.to_string()))
    }
}

#[cfg(test)]
mod tests;
