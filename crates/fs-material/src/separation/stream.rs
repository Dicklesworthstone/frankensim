//! Dry-solid classes, explicit solid-phase element inventories and carrier water.

use super::{Partition, SeparationError, nonnegative, positive, require};
use fs_qty::{Dimensionless, Length, MassFlowRate, chemistry::ElementId};
use std::collections::{BTreeMap, BTreeSet};

/// A caller-defined homogeneous particle class; its name is not a mineral assay.
#[derive(Clone, Debug, PartialEq)]
pub struct ParticleClass {
    id: String,
    diameter: Length,
    dry_solids: MassFlowRate,
    elements: BTreeMap<ElementId, MassFlowRate>,
}

impl ParticleClass {
    /// Admit a representative diameter and solid-phase inventories in kg/s.
    /// Element identities are caller-declared measurands. Partial elemental
    /// coverage is allowed, but no missing element becomes a measured zero.
    pub fn new(
        id: impl Into<String>,
        diameter: Length,
        dry_solids: MassFlowRate,
        elements: BTreeMap<ElementId, MassFlowRate>,
    ) -> Result<Self, SeparationError> {
        let id = id.into();
        require(
            !id.trim().is_empty() && !id.chars().any(char::is_control),
            "particle class needs a usable identity",
        )?;
        positive(
            diameter.value(),
            "representative particle size must be positive",
        )?;
        let mass = nonnegative(
            dry_solids.value(),
            "dry-solids flow must be finite and nonnegative",
        )?;
        for amount in elements.values() {
            nonnegative(
                amount.value(),
                "element flow must be finite and nonnegative",
            )?;
            require(
                amount.value() <= mass,
                "element flow exceeds its class dry-solids flow",
            )?;
        }
        let known = sum(elements.values().map(|m| m.value()))?;
        // Independent rounded elemental entries can differ in their sum by a
        // few ULPs. Retain every value; never normalize or fabricate gangue.
        require(
            known <= mass || (mass > 0.0 && (known - mass) / mass <= 16.0 * f64::EPSILON),
            "known elements exceed dry-solids flow",
        )?;
        Ok(Self {
            id,
            diameter,
            dry_solids,
            elements,
        })
    }

    /// Original particle/mineral/liberation-class identity supplied by the caller.
    pub fn id(&self) -> &str {
        &self.id
    }
    /// Representative particle diameter in metres; no within-bin PSD is inferred.
    pub const fn diameter(&self) -> Length {
        self.diameter
    }
    /// Total dry-solids flow, including unassayed constituents, in kg/s.
    pub const fn dry_solids(&self) -> MassFlowRate {
        self.dry_solids
    }
    /// Explicit solid-phase elemental mass flows, including declared zeroes.
    pub fn elements(&self) -> &BTreeMap<ElementId, MassFlowRate> {
        &self.elements
    }

    fn scaled(&self, fraction: f64) -> Result<Self, SeparationError> {
        let elements = self
            .elements
            .iter()
            .map(|(element, flow)| {
                Ok((
                    element.clone(),
                    MassFlowRate::new(product(flow.value(), fraction)?),
                ))
            })
            .collect::<Result<_, SeparationError>>()?;
        Self::new(
            self.id.clone(),
            self.diameter,
            MassFlowRate::new(product(self.dry_solids.value(), fraction)?),
            elements,
        )
    }
}

/// Immutable size-resolved solids and carrier water, each in kg/s.
/// Dissolved solutes, solution volume and aqueous reactions are not represented.
#[derive(Clone, Debug, PartialEq)]
pub struct ParticleStream {
    classes: Vec<ParticleClass>,
    water: MassFlowRate,
}

impl ParticleStream {
    /// Admit distinct class identities with the SAME explicit elemental axes.
    /// Classes are stored in canonical identity order. Zero flows retain their
    /// declared layout; the solids fraction of an entirely empty flow refuses.
    pub fn new(
        mut classes: Vec<ParticleClass>,
        water: MassFlowRate,
    ) -> Result<Self, SeparationError> {
        require(
            !classes.is_empty(),
            "particle stream requires its declared class layout",
        )?;
        nonnegative(water.value(), "water flow must be finite and nonnegative")?;
        classes.sort_by(|left, right| left.id.cmp(&right.id));
        require(
            classes.windows(2).all(|pair| pair[0].id != pair[1].id),
            "duplicate particle class identity",
        )?;
        let axes: BTreeSet<_> = classes[0].elements.keys().collect();
        for class in &classes {
            require(
                class.elements.keys().collect::<BTreeSet<_>>() == axes,
                "every class must explicitly declare the same elemental axes",
            )?;
        }
        let stream = Self { classes, water };
        stream.dry_solids()?;
        Ok(stream)
    }

    /// Classes in canonical order, with original composition and size identities.
    pub fn classes(&self) -> &[ParticleClass] {
        &self.classes
    }
    /// Carrier-water mass flow in kg/s, not a slurry volume.
    pub const fn water(&self) -> MassFlowRate {
        self.water
    }
    /// Sum of all dry-solid classes, in kg/s.
    pub fn dry_solids(&self) -> Result<MassFlowRate, SeparationError> {
        Ok(MassFlowRate::new(sum(self
            .classes
            .iter()
            .map(|c| c.dry_solids.value()))?))
    }
    /// Actual solids mass fraction; no mineral or water density is assumed.
    pub fn solids_mass_fraction(&self) -> Result<Dimensionless, SeparationError> {
        let solids = self.dry_solids()?.value();
        let scale = solids.max(self.water.value());
        require(
            scale > 0.0,
            "empty flow has no defined solids mass fraction",
        )?;
        let fraction = (solids / scale) / (solids / scale + self.water.value() / scale);
        Ok(Dimensionless::new(fraction))
    }
    /// Total flow of an explicitly declared solid-phase element, in kg/s.
    /// An unassayed element refuses instead of returning a fabricated zero.
    pub fn element_flow(&self, element: &ElementId) -> Result<MassFlowRate, SeparationError> {
        require(
            self.classes[0].elements.contains_key(element),
            "element is not assayed in this stream",
        )?;
        Ok(MassFlowRate::new(sum(self
            .classes
            .iter()
            .map(|class| class.elements[element].value()))?))
    }
    /// Element / total dry solids by mass. `None` denotes a zero-solids stream.
    pub fn element_grade(
        &self,
        element: &ElementId,
    ) -> Result<Option<Dimensionless>, SeparationError> {
        let element = self.element_flow(element)?.value();
        let total = self.dry_solids()?.value();
        if total == 0.0 {
            return Ok(None);
        }
        Ok(Some(Dimensionless::new(super::fraction(element / total)?)))
    }

    /// Combine streams with exactly the same class identities, representative
    /// diameters and declared elemental axes. Missing assays never become zeroes.
    /// Each extensive rate is summed with compensation in magnitude order;
    /// inlet permutation cannot change the result. No composition is normalized.
    pub fn mix(streams: &[&Self]) -> Result<Self, SeparationError> {
        let first = streams.first().ok_or(SeparationError::Invalid(
            "mixing requires at least one particle stream",
        ))?;
        for stream in &streams[1..] {
            require(
                stream.classes.len() == first.classes.len(),
                "mixed particle streams require the same class layout",
            )?;
            for (left, right) in first.classes.iter().zip(&stream.classes) {
                require(
                    left.id == right.id && left.diameter == right.diameter,
                    "mixed class identities and representative diameters must match",
                )?;
                require(
                    left.elements.keys().eq(right.elements.keys()),
                    "mixed classes require the same explicit elemental axes",
                )?;
            }
        }
        let classes = first
            .classes
            .iter()
            .enumerate()
            .map(|(index, class)| {
                let dry_solids = mix_sum(
                    streams
                        .iter()
                        .map(|stream| stream.classes[index].dry_solids.value()),
                )?;
                let elements =
                    class
                        .elements
                        .keys()
                        .map(|element| {
                            Ok((
                                element.clone(),
                                MassFlowRate::new(mix_sum(streams.iter().map(|stream| {
                                    stream.classes[index].elements[element].value()
                                }))?),
                            ))
                        })
                        .collect::<Result<_, SeparationError>>()?;
                ParticleClass::new(
                    class.id.clone(),
                    class.diameter,
                    MassFlowRate::new(dry_solids),
                    elements,
                )
            })
            .collect::<Result<Vec<_>, SeparationError>>()?;
        Self::new(
            classes,
            MassFlowRate::new(mix_sum(streams.iter().map(|stream| stream.water.value()))?),
        )
    }

    /// Multiply every dry-solid, elemental and water flow by a finite factor
    /// greater than or equal to zero. Factors above one are admitted. A zero
    /// factor retains every class and measured-zero axis; positive underflow
    /// and overflow refuse without modifying the original stream.
    pub fn scaled(&self, factor: Dimensionless) -> Result<Self, SeparationError> {
        let factor = nonnegative(
            factor.value(),
            "stream scale must be finite and nonnegative",
        )?;
        Self::new(
            self.classes
                .iter()
                .map(|class| class.scaled(factor))
                .collect::<Result<_, _>>()?,
            MassFlowRate::new(product(self.water.value(), factor)?),
        )
    }

    /// Least representative diameter whose dry-mass cumulative fraction reaches
    /// the requested quantile in `(0, 1]`. Equal diameters form one mass group;
    /// zero-mass classes do not set a passing size. No within-bin interpolation,
    /// continuous PSD, breakage or liberation is inferred. Carrier water is excluded.
    pub fn discrete_passing_size(
        &self,
        fraction: Dimensionless,
    ) -> Result<Length, SeparationError> {
        let fraction = positive(
            fraction.value(),
            "passing fraction must be finite and positive",
        )?;
        require(fraction <= 1.0, "passing fraction cannot exceed one")?;
        let total = positive(
            self.dry_solids()?.value(),
            "passing size requires positive dry solids",
        )?;
        let mut classes: Vec<_> = self
            .classes
            .iter()
            .filter(|class| class.dry_solids.value() > 0.0)
            .collect();
        classes.sort_by(|left, right| {
            left.diameter
                .value()
                .total_cmp(&right.diameter.value())
                .then_with(|| left.id.cmp(&right.id))
        });
        let largest = classes
            .last()
            .ok_or(SeparationError::Invalid(
                "passing size requires a positive class",
            ))?
            .diameter;
        // A rounded total must never erase a positive coarse tail at the 100th percentile.
        if fraction == 1.0 {
            return Ok(largest);
        }
        let mut start = 0;
        let mut cumulative = 0.0_f64;
        let mut correction = 0.0_f64;
        while start < classes.len() {
            let diameter = classes[start].diameter;
            let mut end = start + 1;
            while end < classes.len() && classes[end].diameter == diameter {
                end += 1;
            }
            let group_mass = sum(classes[start..end]
                .iter()
                .map(|class| class.dry_solids.value()))?;
            let next = cumulative + group_mass;
            correction += if cumulative >= group_mass {
                (cumulative - next) + group_mass
            } else {
                (group_mass - next) + cumulative
            };
            cumulative = next;
            let prefix = cumulative + correction;
            if !prefix.is_finite() {
                return Err(SeparationError::Arithmetic(
                    "cumulative dry-solids flow is unrepresentable",
                ));
            }
            // Dividing avoids underflow of total*fraction for tiny positive quantiles.
            if prefix / total >= fraction {
                return Ok(diameter);
            }
            start = end;
        }
        // The complete discrete distribution has CDF=1 by definition, including
        // when its size-ordered sum differs by rounding from the ID-ordered total.
        Ok(largest)
    }

    /// Split every homogeneous solid class with its explicit partition and water
    /// with an independent partition. Input order is arbitrary; exactly one
    /// supplied split per class is required. Both small tails are retained, and
    /// positive underflow refuses. Inputs are unchanged on success and failure.
    pub fn split(
        &self,
        partitions: &[(String, Partition)],
        water: Partition,
    ) -> Result<SeparatedStreams, SeparationError> {
        require(
            partitions.len() == self.classes.len(),
            "supply one partition for every particle class",
        )?;
        let mut by_id = BTreeMap::new();
        for (id, partition) in partitions {
            require(
                by_id.insert(id.as_str(), *partition).is_none(),
                "duplicate class partition",
            )?;
        }
        let mut products = Vec::with_capacity(self.classes.len());
        let mut rejects = Vec::with_capacity(self.classes.len());
        for class in &self.classes {
            let partition = by_id
                .get(class.id.as_str())
                .ok_or(SeparationError::Invalid(
                    "missing or unknown class partition",
                ))?;
            products.push(class.scaled(partition.product().value())?);
            rejects.push(class.scaled(partition.reject().value())?);
        }
        Ok(SeparatedStreams {
            product: Self::new(
                products,
                MassFlowRate::new(product(self.water.value(), water.product().value())?),
            )?,
            reject: Self::new(
                rejects,
                MassFlowRate::new(product(self.water.value(), water.reject().value())?),
            )?,
        })
    }
}

/// Fully formed product/reject streams; no partially applied separation is exposed.
#[derive(Clone, Debug, PartialEq)]
pub struct SeparatedStreams {
    /// Concentrate/underflow according to the caller's explicit partition meaning.
    pub product: ParticleStream,
    /// Tailings/overflow with the complementary solids, element and water flows.
    pub reject: ParticleStream,
}

fn product(flow: f64, fraction: f64) -> Result<f64, SeparationError> {
    let result = flow * fraction;
    if !result.is_finite() || (flow > 0.0 && fraction > 0.0 && result == 0.0) {
        Err(SeparationError::Arithmetic(
            "separated flow is unrepresentable",
        ))
    } else {
        Ok(result)
    }
}

fn sum(values: impl Iterator<Item = f64>) -> Result<f64, SeparationError> {
    let mut total = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = total + value;
        correction += if total >= value {
            (total - next) + value
        } else {
            (value - next) + total
        };
        total = next;
    }
    let result = total + correction;
    if result.is_finite() {
        Ok(result)
    } else {
        Err(SeparationError::Arithmetic(
            "stream total is unrepresentable",
        ))
    }
}

fn mix_sum(values: impl Iterator<Item = f64>) -> Result<f64, SeparationError> {
    let mut values: Vec<_> = values.collect();
    values.sort_by(f64::total_cmp);
    sum(values.into_iter())
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
