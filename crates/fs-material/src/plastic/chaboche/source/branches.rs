//! Named source roles for the existing multi-branch Chaboche constitutive law.
//!
//! Every branch has an independently selected modulus and recovery coefficient.
//! Array order is history order, even for an explicitly inactive branch. The
//! complete original Voce bundle and all extra physics receipts are retained.
//! These are fixed-condition coefficients, not a temperature/history transfer.

use super::*;
use crate::plastic::chaboche::MAX_BACKSTRESS_BRANCHES;
use crate::plastic::source::{
    KINEMATIC_HARDENING_MODULUS_PROPERTY, POST_YIELD_TANGENT_MODULUS_PROPERTY,
};
use crate::state_point::DENSITY_PROPERTY;

/// Two explicit property names for one independently calibrated backstress.
/// A single total modulus, a monotonic tangent or a material name cannot define
/// this decomposition. Names must be distinct across every branch and backbone.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BackstressProperties {
    /// C [Pa] in `dX = (2/3) C d(eps_p) - gamma X d(alpha)`.
    pub modulus: String,
    /// Gamma [per unit equivalent plastic strain], a dimensionless coefficient.
    pub recovery: String,
}

impl BackstressProperties {
    /// Name one ordered source pair. Requirements/resolution validate the full
    /// schema together so aliases between separate branches cannot slip through.
    pub fn new(modulus: impl Into<String>, recovery: impl Into<String>) -> Self {
        Self {
            modulus: modulus.into(),
            recovery: recovery.into(),
        }
    }
}

/// Exact Voce/backbone requirements plus every independently named C/gamma pair.
/// The schema admits one through eight branches, preserves their order and
/// refuses aliases with other branches, backbone roles or reference density.
/// Explicit C=0 and gamma=0 are admitted; missing data never select these limits.
pub fn chaboche_branch_requirements<const N: usize>(
    basis: HardeningPropertyBasis,
    properties: &[BackstressProperties; N],
) -> Result<Vec<ScalarPropertyRequirement>, J2SourceError> {
    if N == 0 || N > MAX_BACKSTRESS_BRANCHES {
        return Err(J2SourceError::Invalid(
            "source-selected Chaboche requires one to eight explicit branches",
        ));
    }
    admit_basis(basis)?;
    let mut requirements = voce_requirements(basis)?;
    for branch in properties {
        for (name, dims) in [
            (&branch.modulus, Pressure::DIMS),
            (&branch.recovery, Dims([0; 6])),
        ] {
            if name.trim() != name
                || requirements
                    .iter()
                    .any(|requirement| requirement.name() == name)
                || [
                    KINEMATIC_HARDENING_MODULUS_PROPERTY,
                    POST_YIELD_TANGENT_MODULUS_PROPERTY,
                    DENSITY_PROPERTY,
                ]
                .contains(&name.as_str())
            {
                return Err(J2SourceError::Invalid(
                    "Chaboche source roles must have distinct exact names without branch, backbone or density aliases",
                ));
            }
            requirements.push(ScalarPropertyRequirement::try_new(
                name.clone(),
                dims,
                ScalarAdmissibility::NonNegative,
            )?);
        }
    }
    Ok(requirements)
}

/// Native N-branch cyclic hardening and its full immutable source snapshot.
/// The ordered source schema is retained beside the same native branch array;
/// no sorting, aggregation or inferred allocation of a total modulus occurs.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedChabocheBranches<const N: usize> {
    backbone: ResolvedVoceJ2,
    law: ChabocheJ2<N>,
    properties: [BackstressProperties; N],
}

impl<const N: usize> ResolvedChabocheBranches<N> {
    /// Adopt the complete original bundle without dropping density, transport,
    /// or other consumers' receipts. Every branch uses the fixed-condition and
    /// exact-quantity admission of the original single-branch source owner.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: HardeningPropertyBasis,
        strain_limit: f64,
        properties: [BackstressProperties; N],
    ) -> Result<Self, J2SourceError> {
        chaboche_branch_requirements(basis, &properties)?;
        let mut branches = [BackstressBranch::new(0.0, 0.0)?; N];
        for (index, names) in properties.iter().enumerate() {
            branches[index] = BackstressBranch::new(
                coefficient(&material, &names.modulus, Pressure::DIMS)?.value_si(),
                coefficient(&material, &names.recovery, Dims([0; 6]))?.value_si(),
            )?;
        }
        let backbone = ResolvedVoceJ2::try_from_material(material, basis, strain_limit)?;
        let law = ChabocheJ2::new(backbone.law().clone(), branches)?;
        Ok(Self {
            backbone,
            law,
            properties,
        })
    }

    /// Full shared material bundle, including all extra roles and their receipts.
    pub fn material(&self) -> &ResolvedMaterialStatePoint {
        self.backbone.material()
    }

    /// Original source-bound Voce, isotropic and optional linear Prager owner.
    pub const fn backbone(&self) -> &ResolvedVoceJ2 {
        &self.backbone
    }

    /// Existing native multi-backstress return, energy and analytic tangent.
    pub const fn law(&self) -> &ChabocheJ2<N> {
        &self.law
    }

    /// Explicit independent isotropic or mixed-hardening basis.
    pub const fn hardening_basis(&self) -> HardeningPropertyBasis {
        self.backbone.hardening_basis()
    }

    /// Source-pair order is exactly native backstress/history order.
    pub const fn properties(&self) -> &[BackstressProperties; N] {
        &self.properties
    }

    /// Each original C claim, uncertainty and receipt in native branch order.
    pub fn modulus_sources(&self) -> [&ResolvedScalarProperty; N] {
        core::array::from_fn(|i| {
            self.material()
                .property(&self.properties[i].modulus)
                .expect("the immutable branch bundle retains every admitted modulus")
        })
    }

    /// Each original gamma claim, uncertainty and receipt in native branch order.
    pub fn recovery_sources(&self) -> [&ResolvedScalarProperty; N] {
        core::array::from_fn(|i| {
            self.material()
                .property(&self.properties[i].recovery)
                .expect("the immutable branch bundle retains every admitted recovery coefficient")
        })
    }
}

impl<const N: usize> SmallStrainLaw for ResolvedChabocheBranches<N> {
    type State = ChabocheState<N>;
    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.law.stress(strain, state)
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.law.tangent(strain, state)
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.law.update_state(strain, state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.law.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.law.card()
    }
}

/// Resolve one fixed-condition source snapshot for an ordered branch schema.
/// Independent claims and explicit selection are required for every pair.
pub fn resolve_chaboche_branches<const N: usize>(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    basis: HardeningPropertyBasis,
    strain_limit: f64,
    properties: [BackstressProperties; N],
) -> Result<ResolvedChabocheBranches<N>, J2SourceError> {
    let material = resolve_material_state_point(
        card,
        point,
        &chaboche_branch_requirements(basis, &properties)?,
        selection,
    )?;
    ResolvedChabocheBranches::try_from_material(material, basis, strain_limit, properties)
}
