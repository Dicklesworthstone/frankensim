//! Source-card viscosity and conductivity for the existing gas transport owners.
//!
//! `ThermochemicalGasModel` owns pure-species material resolution, Wilke/WMS
//! mixing and the original caloric identities. A separately supplied, exactly
//! source-matched `TransportLaw` owns only binary diffusivities here; its scalar
//! viscosity/conductivity outputs are explicitly unused. Bulk viscosity is an
//! explicit nonnegative declaration, never inferred from either source.

use crate::multicomponent::dissipation::{
    DissipationError, TransportCoefficients, TransportLaw, validate_binary,
};
use fs_matdb::QueryPoint;
use fs_material::MaterialError;
use fs_material::thermochemical_gas::{ResolvedThermochemicalGasState, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, DynViscosity, Pressure, QuantitySpec, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{Composition, CompositionBasis};

/// Material-card transport composed with an explicit binary-diffusion owner.
///
/// All original species, material cards, claim-selection policies, pack receipts
/// and external query coordinates remain available from this object. Runtime
/// face T,p replace only those two coordinates, with their canonical physical
/// quantity meanings. Composition belongs to the original canonical species
/// axis; exact-zero species are omitted from the material query, without adding
/// trace concentrations. The diffusion owner still supplies its full pair axis.
/// Native face fractions within their roundoff-only unit-sum envelope use the
/// original material owner's generated-simplex closure. Both coefficient owners
/// receive that same canonical composition; genuinely nonunit inputs refuse.
///
/// This adapter performs original material queries and their allocations on
/// each evaluation. It does not claim the allocation-free behavior of the
/// separate collision-data `KineticTransport` adapter. It contains no accepted
/// fluid state or mutable constitutive history.
pub struct MaterialTransport<L> {
    material: ThermochemicalGasModel,
    context: QueryPoint,
    diffusion: L,
    bulk_viscosity: DynViscosity,
    binary_workspace: Vec<f64>,
}
impl<L: TransportLaw> MaterialTransport<L> {
    /// Bind complete original caloric identity before any material or diffusion query.
    /// The material bank must contain every species on the diffusion owner's axis,
    /// including species that may be absent in a particular cell.
    pub fn new(
        material: ThermochemicalGasModel,
        context: QueryPoint,
        diffusion: L,
        bulk_viscosity: DynViscosity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if material.thermodynamics().components() != diffusion.source().models() {
            return Err(DissipationError::Invalid(
                "material and binary diffusion caloric sources differ",
            ));
        }
        if !bulk_viscosity.value().is_finite() || bulk_viscosity.value() < 0.0 {
            return Err(DissipationError::Invalid(
                "material transport bulk viscosity",
            ));
        }
        let n = diffusion.source().models().len();
        let count = n
            .checked_mul(n)
            .ok_or(DissipationError::Invalid("binary diffusion size overflow"))?;
        let mut binary_workspace = Vec::new();
        binary_workspace
            .try_reserve_exact(count)
            .map_err(|_| DissipationError::Flow(crate::multicomponent::Error::Allocation))?;
        binary_workspace.resize(count, 0.0);
        poll(&mut cancelled)?;
        Ok(Self {
            material,
            context,
            diffusion,
            bulk_viscosity,
            binary_workspace,
        })
    }

    /// Original complete material bank, retaining cards, selection and source packs.
    pub const fn material(&self) -> &ThermochemicalGasModel {
        &self.material
    }
    /// Authored external query context, unchanged by all face evaluations.
    pub const fn context(&self) -> &QueryPoint {
        &self.context
    }
    /// Original pair-transport law; only its binary diffusivities are consumed.
    pub const fn diffusion_law(&self) -> &L {
        &self.diffusion
    }
    /// Explicit bulk viscosity, including exact zero for the Stokes assumption.
    pub const fn bulk_viscosity(&self) -> DynViscosity {
        self.bulk_viscosity
    }

    /// Resolve the actual common material state and its original query receipts.
    ///
    /// This does not query binary diffusion or mutate any source. It provides
    /// callers the same source-bearing material state used by `evaluate`, not a
    /// second formula for mixture properties. A refusal returns no partial state.
    pub fn resolve(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalGasState, DissipationError> {
        self.resolve_face(temperature, pressure, mole_fractions, &mut cancelled)
            .map(|(resolved, _)| resolved)
    }

    fn resolve_face(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        cancelled: &mut (impl FnMut() -> bool + ?Sized),
    ) -> Result<(ResolvedThermochemicalGasState, Vec<f64>), DissipationError> {
        poll(cancelled)?;
        admit_state(
            self.diffusion.source().models().len(),
            temperature,
            pressure,
            mole_fractions,
        )?;
        let active = if mole_fractions.iter().sum::<f64>() == 1.0 {
            let composition =
                Composition::new(CompositionBasis::MoleFraction, mole_fractions.to_vec())
                    .map_err(|error| material_input("composition", error))?;
            self.material.with_composition(composition)
        } else {
            // The conservative face owner closes its compensated sum; the
            // material owner requires the canonical ordered sum to be exact.
            // Only already admitted roundoff is reconciled, by the original
            // generated-composition owner, not a second normalization formula.
            let amounts: Vec<_> = mole_fractions.iter().map(|&x| Amount::new(x)).collect();
            self.material.with_mole_amounts(&amounts)
        }
        .map_err(DissipationError::Material)?;
        // Restore the complete diffusion species axis from the actual active
        // material model. A removed species remains exactly zero in both laws.
        let mut canonical = vec![0.0; mole_fractions.len()];
        let mut active_species = active
            .thermodynamics()
            .components()
            .iter()
            .zip(active.thermodynamics().mole_composition().fractions());
        for ((original, &fraction), output) in self
            .material
            .thermodynamics()
            .components()
            .iter()
            .zip(mole_fractions)
            .zip(&mut canonical)
        {
            if fraction > 0.0 {
                let Some((component, &resolved)) = active_species.next() else {
                    return Err(DissipationError::Invalid(
                        "material transport lost an active species",
                    ));
                };
                if original != component || resolved <= 0.0 {
                    return Err(DissipationError::Invalid(
                        "material transport changed an active species",
                    ));
                }
                *output = resolved;
            }
        }
        if active_species.next().is_some() {
            return Err(DissipationError::Invalid(
                "material transport introduced an absent species",
            ));
        }
        let point = self
            .context
            .clone()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                temperature.value(),
            )
            .and_then(|p| {
                p.with_quantity("pressure", typed(QuantityKind::Pressure), pressure.value())
            })
            .map_err(|error| material_input("query coordinates", error))?;
        poll(cancelled)?;
        let state = active
            .evaluate(&point)
            .map_err(DissipationError::Material)?;
        poll(cancelled)?;
        Ok((state, canonical))
    }
}
impl<L: TransportLaw> TransportLaw for MaterialTransport<L> {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.diffusion.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        mole_fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        poll(cancelled)?;
        if binary_diffusivities.len() != self.binary_workspace.len() {
            return Err(DissipationError::Invalid(
                "material transport binary output axis",
            ));
        }
        // Admission precedes BOTH original-source queries; only the original
        // material owner may reconcile native unit-sum floating-point roundoff.
        let (resolved, canonical) =
            self.resolve_face(temperature, pressure, mole_fractions, cancelled)?;
        // Deliberate ownership: this law supplies only D_ij. Its mu/k/bulk are
        // not fallbacks for missing material cards or the explicit bulk input.
        let _unused_scalars = self.diffusion.evaluate(
            temperature,
            pressure,
            &canonical,
            &mut self.binary_workspace,
            cancelled,
        )?;
        validate_binary(&self.binary_workspace, self.source().models().len())?;
        let gas = resolved.state();
        let coefficients = TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(gas.dynamic_viscosity),
            bulk_viscosity: self.bulk_viscosity,
            thermal_conductivity: ThermalConductivityQuantityV1::new(gas.thermal_conductivity),
        };
        poll(cancelled)?;
        binary_diffusivities.copy_from_slice(&self.binary_workspace);
        Ok(coefficients)
    }
}

fn admit_state(n: usize, t: Temperature, p: Pressure, x: &[f64]) -> Result<(), DissipationError> {
    if !t.value().is_finite()
        || t.value() <= 0.0
        || !p.value().is_finite()
        || p.value() <= 0.0
        || x.len() != n
        || x.iter().any(|v| !v.is_finite() || *v < 0.0)
        || (x.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON * n as f64
    {
        return Err(DissipationError::Invalid(
            "material transport T,p or canonical mole fractions outside native roundoff",
        ));
    }
    Ok(())
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn material_input(stage: &str, error: impl core::fmt::Display) -> DissipationError {
    DissipationError::Material(MaterialError::Parameters {
        what: format!("gas transport {stage}: {error}"),
    })
}
fn poll(cancelled: &mut (impl FnMut() -> bool + ?Sized)) -> Result<(), DissipationError> {
    if cancelled() {
        Err(DissipationError::Cancelled)
    } else {
        Ok(())
    }
}
