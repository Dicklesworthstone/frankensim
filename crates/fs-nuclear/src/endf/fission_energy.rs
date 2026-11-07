//! Evaluated fission-energy components from ENDF-6 MF1/MT458.
//!
//! The nine source components, their uncertainties and any component-specific
//! TAB1 overrides remain distinct. Fragment, neutron and photon kinetic energy
//! are not automatically local heat. Delayed components do not become prompt
//! deposition, and an evaluated total does not replace its partial components.
//!
//! `NPLY=0` supplies thermal/threshold reference values. An untabulated component
//! therefore refuses an incident-energy query until a separate caller-visible
//! continuation/systematics policy is supplied. Polynomial coefficients retain
//! their native ENDF electronvolt convention; no historical unit-error repair,
//! spectrum, transport, lifetime or uncertainty propagation is inferred.
//! Format: ENDF-102 (2023), section 1.5,
//! <https://www.nndc.bnl.gov/endfdocs/ENDF-102-2023.pdf>.

use std::sync::Arc;

use super::reaction_weighted::{TabulatedYield, TabulatedYieldData, TabulatedYieldLimits};
use super::{
    EV_J, EndfError, EndfMetadata, EndfReadLimits, EndfSelection, Reader, Record, budget,
    convert_si, read_metadata, read_nonnegative_tab1,
};
use crate::cross_sections::NeutronReaction;
use crate::irradiation::NuclideId;
use crate::{FissionError, nonnegative};

/// Explicit energy dependence and conservative reaction-weighted group means.
pub mod model;

/// Prompt outgoing products on the original fission-product burnup owner.
pub mod burnup;

/// Evaluated prompt products across processed temperature knots.
pub mod temperature;

/// Original MF1/MT458 component identity, in source IFC order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FissionEnergyComponent {
    /// EFR: fragment kinetic energy after prompt neutron emission.
    Fragments,
    /// ENP: prompt-neutron kinetic energy.
    PromptNeutrons,
    /// END: delayed-neutron kinetic energy.
    DelayedNeutrons,
    /// EGP: prompt-photon energy.
    PromptPhotons,
    /// EGD: delayed-photon energy.
    DelayedPhotons,
    /// EB: delayed beta energy.
    Betas,
    /// ENU: energy carried by neutrinos.
    Neutrinos,
    /// ER: source-declared total excluding neutrinos; not an additional channel.
    Recoverable,
    /// ET: source-declared total; not an additional channel.
    Total,
}

impl FissionEnergyComponent {
    pub const ALL: [Self; 9] = [
        Self::Fragments,
        Self::PromptNeutrons,
        Self::DelayedNeutrons,
        Self::PromptPhotons,
        Self::DelayedPhotons,
        Self::Betas,
        Self::Neutrinos,
        Self::Recoverable,
        Self::Total,
    ];

    /// Zero-based index; the corresponding ENDF IFC is `index() + 1`.
    #[must_use]
    pub const fn index(self) -> usize {
        match self {
            Self::Fragments => 0,
            Self::PromptNeutrons => 1,
            Self::DelayedNeutrons => 2,
            Self::PromptPhotons => 3,
            Self::DelayedPhotons => 4,
            Self::Betas => 5,
            Self::Neutrinos => 6,
            Self::Recoverable => 7,
            Self::Total => 8,
        }
    }
}

/// Exact MF1 identity selection. Zero-temperature original evaluations are
/// admitted independently of the processed MF3 cross-section temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EndfFissionEnergySelection {
    pub material_number: u16,
    pub nuclide: NuclideId,
    pub temperature_k: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EndfFissionEnergyReadLimits {
    pub records: EndfReadLimits,
    /// Highest retained polynomial order, checked before coefficient allocation.
    pub max_polynomial_order: usize,
    /// Total energy/ordinate scalar values across all IFC tables, two per point.
    pub max_table_values: usize,
}

impl Default for EndfFissionEnergyReadLimits {
    fn default() -> Self {
        Self {
            records: EndfReadLimits::default(),
            max_polynomial_order: 16,
            max_table_values: 2_000_000,
        }
    }
}

/// Original coefficient and its declared nonnegative uncertainty. At order k,
/// BOTH values use the source convention eV^(1-k): the polynomial accepts the
/// incident energy numerically in eV and produces energy numerically in eV.
/// Only the zeroth coefficient is itself an energy. Higher coefficients may
/// have either sign. These source uncertainties are not numerical allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FissionEnergyCoefficient {
    pub value: f64,
    pub uncertainty: f64,
}

/// Component-specific energy dependence in the original histogram/linear law.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionEnergyTable {
    derivation: u8,
    // The existing dimensionless interpolation/integration owner represents
    // component energy divided by one eV (EV_J joules), not neutron multiplicity.
    curve: TabulatedYield,
}

impl FissionEnergyTable {
    /// Original LDRV: 1 means derived data; 2 means a primary evaluation.
    #[must_use]
    pub const fn derivation(&self) -> u8 {
        self.derivation
    }

    #[must_use]
    pub fn energies_j(&self) -> &[f64] {
        self.curve.energies_j()
    }

    /// Component energy in J/fission, without incident-energy extrapolation.
    /// An exact jump knot uses the original right-hand value.
    pub fn value_at_energy_j(&self, incident_energy_j: f64) -> Result<f64, EndfError> {
        convert_si(self.curve.value_at_energy(incident_energy_j)?, EV_J)
    }
}

/// Immutable original energy-release evaluation. Clones share coefficient and
/// table payloads. Reported totals are retained even when rounded or evaluated
/// independently from the partial components; no normalization is performed.
#[derive(Clone, Debug, PartialEq)]
pub struct FissionEnergyRelease {
    data_id: Arc<str>,
    metadata: EndfMetadata,
    polynomial_order: usize,
    coefficients: [Arc<[FissionEnergyCoefficient]>; 9],
    reference_values_j: [f64; 9],
    tables: [Option<FissionEnergyTable>; 9],
}

impl FissionEnergyRelease {
    pub fn parse(
        data_id: String,
        selection: EndfFissionEnergySelection,
        text: &str,
        limits: EndfFissionEnergyReadLimits,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, EndfError> {
        checkpoint()?;
        if data_id.trim().is_empty()
            || selection.material_number == 0
            || selection.material_number > 9999
            || !nonnegative(selection.temperature_k)
        {
            return Err(invalid(
                "identified fission energy, valid MAT and nonnegative MF1 temperature required",
            ));
        }
        if text.len() > limits.records.max_bytes {
            return Err(budget("ENDF fission-energy bytes"));
        }
        let mut reader = Reader {
            lines: text.lines(),
            count: 0,
            limits: limits.records,
        };
        let data_id: Arc<str> = data_id.into();
        let mat = i32::from(selection.material_number);
        let mut metadata = None;
        let mut candidate = None;
        while let Some(record) = reader.next(&mut checkpoint)? {
            if record.mat != mat {
                metadata = None;
                continue;
            }
            if record.mf == 1 && record.mt == 451 {
                metadata = read_metadata(
                    &mut reader,
                    record,
                    EndfSelection {
                        material_number: selection.material_number,
                        nuclide: selection.nuclide,
                        reaction: NeutronReaction::Fission,
                        temperature_k: selection.temperature_k,
                    },
                    false,
                    &mut checkpoint,
                )?;
            } else if record.mf == 1 && record.mt == 458 {
                if let Some(identity) = &metadata {
                    if candidate.is_some() {
                        return Err(EndfError::AmbiguousSelection);
                    }
                    candidate = Some(read_release(
                        &mut reader,
                        record,
                        Arc::clone(&data_id),
                        identity.clone(),
                        limits,
                        &mut checkpoint,
                    )?);
                }
            }
        }
        checkpoint()?;
        candidate.ok_or(EndfError::MissingSelection)
    }

    #[must_use]
    pub fn data_id(&self) -> &str {
        &self.data_id
    }

    #[must_use]
    pub const fn metadata(&self) -> &EndfMetadata {
        &self.metadata
    }

    #[must_use]
    pub const fn polynomial_order(&self) -> usize {
        self.polynomial_order
    }

    /// Original native-unit coefficient pairs in increasing polynomial order.
    #[must_use]
    pub fn coefficients(&self, component: FissionEnergyComponent) -> &[FissionEnergyCoefficient] {
        &self.coefficients[component.index()]
    }

    #[must_use]
    pub fn tabulation(&self, component: FissionEnergyComponent) -> Option<&FissionEnergyTable> {
        self.tables[component.index()].as_ref()
    }

    /// Original zeroth coefficient in J/fission, including its reference-point
    /// meaning when NPLY=0. This does not evaluate a component-specific table.
    #[must_use]
    pub fn reference_value_j(&self, component: FissionEnergyComponent) -> f64 {
        self.reference_values_j[component.index()]
    }

    /// Evaluate an explicitly supplied table or a positive-order polynomial.
    /// Untabulated NPLY=0 components require a separate explicit physical
    /// continuation policy, even for a query at zero incident energy.
    ///
    /// Signed coefficients are evaluated by Horner/FMA in their source eV
    /// convention. Nonfinite or negative results refuse; no coefficient repair,
    /// clipping, inferred validity domain or error certificate is supplied.
    pub fn value_at_energy_j(
        &self,
        component: FissionEnergyComponent,
        incident_energy_j: f64,
    ) -> Result<f64, EndfError> {
        if !nonnegative(incident_energy_j) {
            return Err(invalid("finite nonnegative incident energy required"));
        }
        if let Some(table) = self.tabulation(component) {
            return table.value_at_energy_j(incident_energy_j);
        }
        if self.polynomial_order == 0 {
            return Err(EndfError::Unsupported(
                "reference-only fission energy requires an explicit energy-dependence policy",
            ));
        }
        let energy_ev = incident_energy_j / EV_J;
        if !nonnegative(energy_ev) {
            return Err(FissionError::Arithmetic.into());
        }
        let coefficients = self.coefficients(component);
        let mut value_ev = coefficients[self.polynomial_order].value;
        for coefficient in coefficients[..self.polynomial_order].iter().rev() {
            value_ev = value_ev.mul_add(energy_ev, coefficient.value);
            if !value_ev.is_finite() {
                return Err(FissionError::Arithmetic.into());
            }
        }
        if value_ev < 0.0 {
            return Err(invalid(
                "fission-energy polynomial is negative at this energy",
            ));
        }
        convert_si(value_ev, EV_J)
    }
}

fn read_release(
    reader: &mut Reader<'_>,
    head: Record<'_>,
    data_id: Arc<str>,
    metadata: EndfMetadata,
    limits: EndfFissionEnergyReadLimits,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<FissionEnergyRelease, EndfError> {
    let za = u32::from(metadata.nuclide.atomic_number()) * 1000
        + u32::from(metadata.nuclide.mass_number());
    if head.number(0)? != f64::from(za) || head.number(1)? != metadata.atomic_weight_ratio {
        return Err(head.error("fission-energy ZA/AWR disagrees with MF1 identity"));
    }
    head.zero_integer(2)?;
    head.zero_integer(4)?;
    let lfc = head.integer(3)?;
    let table_count = head.count(5)?;
    if !matches!(lfc, 0 | 1)
        || table_count > FissionEnergyComponent::ALL.len()
        || (lfc == 0 && table_count != 0)
        || (lfc == 1 && table_count == 0)
    {
        return Err(head.error("fission-energy LFC/NFC must declare zero or 1..9 tables"));
    }
    let list = reader.section(head.mat, 1, 458, checkpoint)?;
    if list.number(0)? != 0.0 || list.number(1)? != 0.0 {
        return Err(list.error("nonzero fission-energy LIST reserved real field"));
    }
    list.zero_integer(2)?;
    let polynomial_order = list.count(3)?;
    if polynomial_order > limits.max_polynomial_order {
        return Err(budget("ENDF fission-energy polynomial order"));
    }
    if lfc == 1 && polynomial_order != 0 {
        return Err(EndfError::Unsupported(
            "tabulated fission energy requires NPLY=0",
        ));
    }
    let orders = polynomial_order
        .checked_add(1)
        .ok_or_else(|| budget("ENDF fission-energy coefficients"))?;
    let pairs = orders
        .checked_mul(9)
        .ok_or_else(|| budget("ENDF fission-energy coefficients"))?;
    let values = pairs
        .checked_mul(2)
        .ok_or_else(|| budget("ENDF fission-energy coefficients"))?;
    if list.count(4)? != values || list.count(5)? != pairs {
        return Err(list.error("fission-energy LIST shape must be 18 and 9 times NPLY+1"));
    }
    if values / 6 + 1 > reader.limits.max_records.saturating_sub(reader.count) {
        return Err(budget("ENDF fission-energy coefficient records"));
    }
    let mut coefficients: [Vec<FissionEnergyCoefficient>; 9] = std::array::from_fn(|_| Vec::new());
    for component in &mut coefficients {
        component
            .try_reserve_exact(orders)
            .map_err(|_| budget("ENDF fission-energy coefficient memory"))?;
    }
    let mut reference_values_j = [0.0; 9];
    for order in 0..orders {
        for row_index in 0..3 {
            let row = reader.section(head.mat, 1, 458, checkpoint)?;
            for pair in 0..3 {
                let value = row.number(pair * 2)?;
                let uncertainty = row.number(pair * 2 + 1)?;
                if (order == 0 && !nonnegative(value)) || !nonnegative(uncertainty) {
                    return Err(row.error(
                        "nonnegative reference energy and coefficient uncertainty required",
                    ));
                }
                let component = row_index * 3 + pair;
                if order == 0 {
                    reference_values_j[component] = convert_si(value, EV_J)?;
                    // Zeroth-order uncertainties are also physical energies.
                    convert_si(uncertainty, EV_J)?;
                }
                coefficients[component].push(FissionEnergyCoefficient { value, uncertainty });
            }
        }
    }
    let mut tables: [Option<FissionEnergyTable>; 9] = std::array::from_fn(|_| None);
    let mut table_values = 0_usize;
    for _ in 0..table_count {
        let tab = reader.section(head.mat, 1, 458, checkpoint)?;
        if tab.number(0)? != 0.0 || tab.number(1)? != 0.0 {
            return Err(tab.error("nonzero fission-energy TAB1 reserved real field"));
        }
        let derivation = tab.integer(2)?;
        let ifc = tab.count(3)?;
        if !matches!(derivation, 1 | 2) || !(1..=9).contains(&ifc) {
            return Err(tab.error("fission-energy TAB1 requires LDRV=1 or 2 and IFC=1..9"));
        }
        if tables[ifc - 1].is_some() {
            return Err(tab.error("duplicate fission-energy IFC table"));
        }
        table_values = tab
            .count(5)?
            .checked_mul(2)
            .and_then(|count| table_values.checked_add(count))
            .ok_or_else(|| budget("ENDF fission-energy table values"))?;
        if table_values > limits.max_table_values {
            return Err(budget("ENDF fission-energy table values"));
        }
        // Retain E_component / (1 eV) in the original dimensionless owner;
        // converting it to joules below cannot silently lose a positive value.
        let source = read_nonnegative_tab1(reader, tab, 1.0, checkpoint)?;
        for (index, &value) in source.values.iter().enumerate() {
            if index % 64 == 0 {
                checkpoint()?;
            }
            convert_si(value, EV_J)?;
        }
        let curve = TabulatedYield::new_with_checkpoint(
            TabulatedYieldData {
                data_id: format!("{data_id} | MF1/MT458 IFC{ifc}"),
                energies_j: source.energies_j,
                values: source.values,
                interpolation_regions: source.interpolation_regions,
            },
            TabulatedYieldLimits {
                max_points: limits.records.max_points,
                max_interpolation_regions: limits.records.max_interpolation_regions,
            },
            &mut *checkpoint,
        )?;
        tables[ifc - 1] = Some(FissionEnergyTable {
            derivation: derivation as u8,
            curve,
        });
    }
    reader.send(head.mat, 1, checkpoint)?;
    checkpoint()?;
    Ok(FissionEnergyRelease {
        data_id,
        metadata,
        polynomial_order,
        coefficients: coefficients.map(Arc::from),
        reference_values_j,
        tables,
    })
}

fn invalid(reason: &'static str) -> EndfError {
    FissionError::InvalidInput(reason).into()
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod model_tests;
