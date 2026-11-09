//! Charged-nucleus inventory products of the original evaluated decay channels.
//!
//! ENDF-102 (2023), §8.4 and §8.4.4.1, defines RTYP 4 as alpha emission,
//! RTYP 7 as proton emission, and successive digits as ordered transformations:
//! https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf . The same single pass
//! resolves the residual and multiplicities; neither spectra nor Q values create
//! extra products. Modes without one represented residual remain unsupported.

use std::{collections::BTreeMap, sync::Arc};

use super::{DecayAssemblyError, EvaluatedDecayAssembly};
use crate::FissionError;
use crate::burnup::BurnupProduct;
use crate::endf::decay::DecayEvaluation;
use crate::irradiation::NuclideId;

/// Explicit retention of emitted charged nuclei in the original burnup inventory.
/// Default retains only the residual daughter. Enabled species require supplied
/// ground-state He4/H1 evaluations with explicitly stable laws. This policy does
/// not infer kinetic stopping, energy deposition, electron capture, bulk material
/// composition, neutron transport, or changes to the residual-only DecayModel.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DecayProductRetention {
    pub alpha_particles: bool,
    pub protons: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct ChargedMultiplicity {
    alpha_particles: u8,
    protons: u8,
}

pub(super) struct ResolvedMode {
    pub residual: NuclideId,
    pub charged: ChargedMultiplicity,
}

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct RetentionBinding {
    policy: DecayProductRetention,
    alpha: Option<usize>,
    proton: Option<usize>,
}

impl RetentionBinding {
    pub(super) fn new(
        policy: DecayProductRetention,
        indices: &BTreeMap<NuclideId, usize>,
        evaluations: &[Arc<DecayEvaluation>],
        checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DecayAssemblyError> {
        let mut bind =
            |enabled: bool, z: u16, a: u16| -> Result<Option<usize>, DecayAssemblyError> {
                checkpoint()?;
                if !enabled {
                    return Ok(None);
                }
                let id = NuclideId::new(z, a, 0)?;
                let index = *indices
                    .get(&id)
                    .ok_or(DecayAssemblyError::MissingRetainedProduct(id))?;
                if evaluations[index].law().decay_constant().value() != 0.0 {
                    return Err(DecayAssemblyError::UnstableRetainedProduct(id));
                }
                Ok(Some(index))
            };
        let alpha = bind(policy.alpha_particles, 2, 4)?;
        let proton = bind(policy.protons, 1, 1)?;
        checkpoint()?;
        Ok(Self {
            policy,
            alpha,
            proton,
        })
    }

    pub(super) fn policy(&self) -> DecayProductRetention {
        self.policy
    }

    pub(super) fn products(
        &self,
        residual: usize,
        charged: ChargedMultiplicity,
    ) -> RetainedProducts {
        let mut products = RetainedProducts {
            entries: [BurnupProduct {
                nuclide: residual,
                yield_per_event: 1.0,
            }; 3],
            len: 1,
        };
        if let Some(index) = self.alpha {
            products.add(index, charged.alpha_particles);
        }
        if let Some(index) = self.proton {
            products.add(index, charged.protons);
        }
        products
    }
}

impl EvaluatedDecayAssembly {
    /// Select retained charged coproducts for subsequent burnup lowering.
    /// The original residual-lineage model, half-lives and event energies stay
    /// shared unchanged. Every retained coproduct uses an already admitted
    /// evaluation; no isotope, lifetime, or calorimetric partition is synthesized.
    /// The consuming builder publishes no partially changed assembly on refusal.
    pub fn with_burnup_product_retention(
        mut self,
        retention: DecayProductRetention,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, DecayAssemblyError> {
        checkpoint()?;
        self.product_retention = RetentionBinding::new(
            retention,
            &self.indices,
            &self.evaluations,
            &mut checkpoint,
        )?;
        Ok(self)
    }

    pub fn burnup_product_retention(&self) -> DecayProductRetention {
        self.product_retention.policy()
    }

    /// A fixed-capacity list is shared by lowering and exact spectral binding.
    /// At most the residual, He4 and H1 can occur, with coincident identities
    /// merged before the original BurnupModel admits its complete product budget.
    pub(super) fn burnup_products(&self, branch: usize) -> RetainedProducts {
        self.product_retention.products(
            self.model.branches()[branch].daughter,
            self.charged_multiplicities[branch],
        )
    }
}

pub(super) struct RetainedProducts {
    entries: [BurnupProduct; 3],
    len: usize,
}

impl RetainedProducts {
    fn add(&mut self, nuclide: usize, count: u8) {
        if count == 0 {
            return;
        }
        let yield_per_event = f64::from(count);
        if let Some(product) = self.entries[..self.len]
            .iter_mut()
            .find(|product| product.nuclide == nuclide)
        {
            product.yield_per_event += yield_per_event;
        } else {
            self.entries[self.len] = BurnupProduct {
                nuclide,
                yield_per_event,
            };
            self.len += 1;
        }
    }

    pub(super) fn as_slice(&self) -> &[BurnupProduct] {
        &self.entries[..self.len]
    }

    pub(super) fn matches(&self, actual: &[BurnupProduct]) -> bool {
        actual.len() == self.len
            && self
                .as_slice()
                .iter()
                .all(|product| actual.contains(product))
    }
}

pub(super) fn resolve_mode(
    parent: NuclideId,
    rtyp: f64,
    isomer: u16,
) -> Result<ResolvedMode, DecayAssemblyError> {
    let unsupported = || DecayAssemblyError::UnsupportedDecayMode { parent, rtyp };
    let invalid = || DecayAssemblyError::InvalidDaughter { parent, rtyp };
    if !rtyp.is_finite() || !(1.0..8.0).contains(&rtyp) {
        return Err(unsupported());
    }
    // The shortest decimal rendering recovers the ENDF particle-code sequence
    // without truncating binary roundoff in e.g. (1.44 - 1) * 100. The source's
    // eleven-column numeric field bounds the available significant digits.
    let code = rtyp.to_string();
    let mut z = i32::from(parent.atomic_number());
    let mut a = i32::from(parent.mass_number());
    let mut charged = ChargedMultiplicity::default();
    for mode in code.bytes().filter(|&byte| byte != b'.') {
        match mode {
            b'1' => z += 1,
            b'2' => z -= 1,
            b'3' => (),
            b'4' => {
                z -= 2;
                a -= 4;
                charged.alpha_particles =
                    charged.alpha_particles.checked_add(1).ok_or_else(invalid)?;
            }
            b'5' => a -= 1,
            b'7' => {
                z -= 1;
                a -= 1;
                charged.protons = charged.protons.checked_add(1).ok_or_else(invalid)?;
            }
            // SF (6) has multiple residuals; zero and unknown codes are not
            // silently treated as no-op isomeric transitions.
            _ => return Err(unsupported()),
        }
        if z <= 0 || a < z || a > i32::from(u16::MAX) {
            return Err(invalid());
        }
    }
    let residual = NuclideId::new(
        u16::try_from(z).map_err(|_| invalid())?,
        u16::try_from(a).map_err(|_| invalid())?,
        isomer,
    )
    .map_err(|_| invalid())?;
    Ok(ResolvedMode { residual, charged })
}
