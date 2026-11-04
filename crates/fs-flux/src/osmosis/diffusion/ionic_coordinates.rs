//! Feasible original-inventory coordinates for closed ionic molecular graphs.
//!
//! One charged species per connected component and one anchor cell eliminate
//! charge and species-total constraints. Every direction therefore lies in
//! their common nullspace. A single radial factor `1/(1+rho)` per component,
//! with rho the directed box gauge, maps that direction into its physical
//! inventory box. This piecewise differentiable chart covers the full relative
//! interior; it neither projects a constitutive flux nor supplies a solver.
//!
//! Floating-point charge and original concentration caps are checked on every
//! actual decoded state. Exact absent species remain absent. Unrepresentable
//! interiors refuse instead of inserting concentration floors. The caller's
//! original ALL-species backward-Euler equations and balance gates still own
//! physical acceptance, including every eliminated ion and anchor inventory.

use super::{
    MAX_COMPARTMENTS, MAX_OSMOTIC_SPECIES, OsmoticNetworkError, ZeroCurrentDiffusion,
    compensated_add, diffusion_source, finite, nonnegative, poll, positive,
};

#[derive(Debug)]
struct Component {
    nodes: Vec<usize>,
    free_species: Vec<usize>,
    pivot: Option<usize>,
    charged_free: usize,
    start: usize,
    end: usize,
}

/// Stage-local immutable coordinate chart. Its original cells/links are not
/// borrowed; the physical stage may reuse separate candidate solution buffers.
#[derive(Debug)]
pub(super) struct IonicCoordinates {
    species: usize,
    volumes: Vec<f64>,
    concentration_caps: Vec<f64>,
    center: Vec<f64>,
    upper: Vec<f64>,
    components: Vec<Component>,
    independent: Vec<usize>,
    scales: Vec<f64>,
    initial: Vec<f64>,
    ionic: ZeroCurrentDiffusion,
}

impl IonicCoordinates {
    /// Conservative numerical payload admission, including constructor grouping
    /// scratch, owned arrays and policy. This is not an allocator/RSS bound.
    pub(super) fn payload_bound(
        nodes: usize,
        species: usize,
        maximum_unknowns: usize,
    ) -> Result<usize, OsmoticNetworkError> {
        if nodes == 0 || nodes > MAX_COMPARTMENTS || species == 0 || species > MAX_OSMOTIC_SPECIES {
            return Err(OsmoticNetworkError::Invalid(
                "bounded ionic coordinate layout required",
            ));
        }
        nodes
            .checked_mul(species)
            .and_then(|entries| entries.checked_mul(12 * core::mem::size_of::<f64>()))
            .and_then(|bytes| {
                bytes.checked_add(nodes.checked_mul(
                    core::mem::size_of::<Component>()
                        + 12 * core::mem::size_of::<usize>()
                        + core::mem::size_of::<Vec<usize>>(),
                )?)
            })
            .and_then(|bytes| {
                bytes.checked_add(maximum_unknowns.checked_mul(16 * core::mem::size_of::<f64>())?)
            })
            .and_then(|bytes| bytes.checked_add(species.checked_mul(core::mem::size_of::<i32>())?))
            .and_then(|bytes| bytes.checked_add(core::mem::size_of::<Self>()))
            .ok_or(OsmoticNetworkError::Invalid(
                "ionic coordinate payload overflow",
            ))
    }

    /// Flat arrays use original cell-major species order. Component IDs name
    /// ALL positive-area molecular components, including carrier-free paths.
    /// Their numeric labels need only be valid node indices; the lowest member
    /// is selected deterministically as the conserved-inventory anchor.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        volumes: &[f64],
        amounts: &[f64],
        concentration_caps: &[f64],
        ionic: &ZeroCurrentDiffusion,
        component_ids: &[usize],
        maximum_unknowns: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, OsmoticNetworkError> {
        poll(cancelled)?;
        let n = ionic.valences().len();
        let nodes = volumes.len();
        Self::payload_bound(nodes, n, maximum_unknowns)?;
        let entries = nodes.checked_mul(n).ok_or(OsmoticNetworkError::Invalid(
            "ionic coordinate inventory layout overflow",
        ))?;
        if amounts.len() != entries
            || concentration_caps.len() != entries
            || component_ids.len() != nodes
            || volumes.iter().any(|&value| !positive(value))
            || amounts.iter().any(|&value| !nonnegative(value))
            || concentration_caps.iter().any(|&value| !positive(value))
            || component_ids.iter().any(|&id| id >= nodes)
        {
            return Err(OsmoticNetworkError::Invalid(
                "complete finite ionic coordinate layout required",
            ));
        }
        for cell in 0..nodes {
            poll(cancelled)?;
            ionic
                .admit_inventories(&amounts[cell * n..(cell + 1) * n])
                .map_err(diffusion_source)?;
            for i in 0..n {
                if amounts[cell * n + i] / volumes[cell] > concentration_caps[cell * n + i] {
                    return Err(OsmoticNetworkError::Invalid(
                        "ionic coordinate inventory exceeds its original concentration cap",
                    ));
                }
            }
        }
        let mut groups: Vec<Vec<usize>> = filled(nodes, Vec::new())?;
        for (cell, &id) in component_ids.iter().enumerate() {
            groups[id].try_reserve(1).map_err(|_| allocation())?;
            groups[id].push(cell);
        }
        let mut chart = Self {
            species: n,
            volumes: copied(volumes)?,
            concentration_caps: copied(concentration_caps)?,
            center: copied(amounts)?,
            upper: filled(entries, 0.0)?,
            components: Vec::new(),
            independent: Vec::new(),
            scales: Vec::new(),
            initial: Vec::new(),
            ionic: ionic.clone(),
        };
        chart
            .components
            .try_reserve_exact(nodes)
            .map_err(|_| allocation())?;
        for members in groups.into_iter().filter(|members| !members.is_empty()) {
            poll(cancelled)?;
            let mut totals = [0.0; MAX_OSMOTIC_SPECIES];
            let mut corrections = [0.0; MAX_OSMOTIC_SPECIES];
            let (mut volume, mut volume_correction) = (0.0, 0.0);
            for &cell in &members {
                compensated_add(&mut volume, &mut volume_correction, volumes[cell]);
                for i in 0..n {
                    compensated_add(&mut totals[i], &mut corrections[i], amounts[cell * n + i]);
                }
            }
            finite(volume)?;
            for &total in &totals[..n] {
                finite(total)?;
            }
            for &cell in &members {
                for i in 0..n {
                    let k = cell * n + i;
                    // An overflowing cap product cannot constrain a finite
                    // conserved pool. Preserve already-admitted q=0 against
                    // a one-ULP multiplication/division discrepancy; decoding
                    // still rechecks the ORIGINAL concentration by division.
                    chart.upper[k] = totals[i]
                        .min(concentration_caps[k] * volumes[cell])
                        .max(amounts[k]);
                }
            }
            if members.len() == 1 {
                continue;
            }
            let mut free_species = Vec::new();
            free_species
                .try_reserve_exact(n)
                .map_err(|_| allocation())?;
            for i in 0..n {
                if totals[i] != 0.0
                    && !members
                        .iter()
                        .all(|&cell| amounts[cell * n + i] == chart.upper[cell * n + i])
                {
                    free_species.push(i);
                }
            }
            let maximum_total = totals[..n].iter().copied().fold(0.0_f64, f64::max);
            let mut pivot = None;
            let mut pivot_weight = -1.0;
            for &i in &free_species {
                if ionic.valences()[i] != 0 {
                    let weight = (totals[i] / maximum_total) * f64::from(ionic.valences()[i]).abs();
                    if weight > pivot_weight {
                        pivot = Some(i);
                        pivot_weight = weight;
                    }
                }
            }
            free_species.retain(|i| Some(*i) != pivot);
            let charged_free = free_species
                .iter()
                .filter(|&&i| ionic.valences()[i] != 0)
                .count();
            let dimension = (members.len() - 1).checked_mul(free_species.len()).ok_or(
                OsmoticNetworkError::Invalid("ionic coordinate dimension overflow"),
            )?;
            let start = chart.independent.len();
            let end = start
                .checked_add(dimension)
                .ok_or(OsmoticNetworkError::Invalid(
                    "ionic coordinate dimension overflow",
                ))?;
            if end > maximum_unknowns {
                return Err(OsmoticNetworkError::Invalid(
                    "ionic coordinate dimension exceeds the declared unknown budget",
                ));
            }
            if dimension == 0 {
                continue;
            }
            let component = Component {
                nodes: members,
                free_species,
                pivot,
                charged_free,
                start,
                end,
            };
            if !chart.interior(&component) {
                chart.mixed_center(&component, amounts, &totals, volume, cancelled)?;
                if !chart.interior(&component) {
                    return Err(OsmoticNetworkError::Invalid(
                        "no representable interior ionic reference inside the original charge and capacity constraints",
                    ));
                }
            }
            chart
                .independent
                .try_reserve_exact(dimension)
                .map_err(|_| allocation())?;
            chart
                .scales
                .try_reserve_exact(dimension)
                .map_err(|_| allocation())?;
            let anchor = component.nodes[0];
            let peers = (component.nodes.len() - 1) as f64;
            for &cell in &component.nodes[1..] {
                for &i in &component.free_species {
                    let index = cell * n + i;
                    let mut scale = chart
                        .radius(index)
                        .min(chart.radius(anchor * n + i) / peers);
                    if ionic.valences()[i] != 0 {
                        let pivot = component
                            .pivot
                            .expect("an active ion has an eliminated charged species");
                        let ratio = (f64::from(ionic.valences()[i])
                            / f64::from(ionic.valences()[pivot]))
                        .abs();
                        let count = component.charged_free as f64;
                        scale = scale
                            .min(chart.radius(cell * n + pivot) / count / ratio)
                            .min(chart.radius(anchor * n + pivot) / peers / count / ratio);
                    }
                    if !positive(scale) {
                        return Err(OsmoticNetworkError::StepUnderflow);
                    }
                    chart.independent.push(index);
                    chart.scales.push(scale);
                }
            }
            chart.components.push(component);
        }
        chart.initial = filled(chart.independent.len(), 0.0)?;
        chart.admit(&chart.center, cancelled)?;
        Ok(chart)
    }

    pub(super) fn dimension(&self) -> usize {
        self.independent.len()
    }
    pub(super) fn independent_indices(&self) -> &[usize] {
        &self.independent
    }
    pub(super) fn initial_guess(&self) -> &[f64] {
        &self.initial
    }

    pub(super) fn heap_bytes(&self) -> usize {
        let vectors = self.volumes.capacity()
            + self.concentration_caps.capacity()
            + self.center.capacity()
            + self.upper.capacity()
            + self.scales.capacity()
            + self.initial.capacity();
        vectors * core::mem::size_of::<f64>()
            + self.independent.capacity() * core::mem::size_of::<usize>()
            + self.components.capacity() * core::mem::size_of::<Component>()
            + self
                .components
                .iter()
                .map(|component| {
                    (component.nodes.capacity() + component.free_species.capacity())
                        * core::mem::size_of::<usize>()
                })
                .sum::<usize>()
            + self.ionic.heap_bytes()
    }

    fn radius(&self, index: usize) -> f64 {
        self.center[index].min(self.upper[index] - self.center[index])
    }

    fn interior(&self, component: &Component) -> bool {
        component.nodes.iter().all(|&cell| {
            component
                .free_species
                .iter()
                .all(|&i| positive(self.radius(cell * self.species + i)))
                && (component.charged_free == 0
                    || positive(
                        self.radius(cell * self.species + component.pivot.expect("charged pivot")),
                    ))
        })
    }

    fn mixed_center(
        &mut self,
        component: &Component,
        original: &[f64],
        totals: &[f64],
        volume: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        let n = self.species;
        for &cell in &component.nodes[1..] {
            poll(cancelled)?;
            let (mut charge, mut correction) = (0.0, 0.0);
            for &i in &component.free_species {
                let k = cell * n + i;
                self.center[k] = finite((self.volumes[cell] / volume) * totals[i])?;
                if let Some(pivot) = component.pivot {
                    compensated_add(
                        &mut charge,
                        &mut correction,
                        (f64::from(self.ionic.valences()[i])
                            / f64::from(self.ionic.valences()[pivot]))
                            * (self.center[k] - original[k]),
                    );
                }
            }
            if component.charged_free != 0 {
                let k = cell * n + component.pivot.expect("charged pivot");
                self.center[k] = finite(original[k] - charge)?;
            }
        }
        let anchor = component.nodes[0];
        for i in 0..n {
            let (mut change, mut correction) = (0.0, 0.0);
            for &cell in &component.nodes[1..] {
                compensated_add(
                    &mut change,
                    &mut correction,
                    self.center[cell * n + i] - original[cell * n + i],
                );
            }
            self.center[anchor * n + i] = finite(original[anchor * n + i] - change)?;
        }
        self.admit_component(component, &self.center, cancelled)
    }

    /// Decode into reusable candidate SCRATCH. On any error/cancellation its
    /// contents are invalid; the enclosing physical stage must not publish it.
    /// This path allocates no vectors and never queries a substitute material.
    pub(super) fn decode_into(
        &self,
        x: &[f64],
        output: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        poll(cancelled)?;
        if x.len() != self.dimension()
            || output.len() != self.center.len()
            || x.iter().any(|value| !value.is_finite())
        {
            return Err(OsmoticNetworkError::Invalid(
                "finite complete ionic coordinates and candidate storage required",
            ));
        }
        // Frozen components retain their original state without graph searches.
        output.copy_from_slice(&self.center);
        let n = self.species;
        for component in &self.components {
            poll(cancelled)?;
            for &cell in &component.nodes {
                output[cell * n..(cell + 1) * n].fill(0.0);
            }
            let scale = x[component.start..component.end]
                .iter()
                .map(|value| value.abs())
                .fold(1.0_f64, f64::max);
            for j in component.start..component.end {
                output[self.independent[j]] = finite((x[j] / scale) * self.scales[j])?;
            }
            if component.charged_free != 0 {
                let pivot = component.pivot.expect("charged pivot");
                for &cell in &component.nodes[1..] {
                    let (mut charge, mut correction) = (0.0, 0.0);
                    for &i in &component.free_species {
                        compensated_add(
                            &mut charge,
                            &mut correction,
                            (f64::from(self.ionic.valences()[i])
                                / f64::from(self.ionic.valences()[pivot]))
                                * output[cell * n + i],
                        );
                    }
                    output[cell * n + pivot] = finite(-charge)?;
                }
            }
            let anchor = component.nodes[0];
            for i in 0..n {
                let (mut sum, mut correction) = (0.0, 0.0);
                for &cell in &component.nodes[1..] {
                    compensated_add(&mut sum, &mut correction, output[cell * n + i]);
                }
                output[anchor * n + i] = finite(-sum)?;
            }
            let mut gauge: f64 = 0.0;
            for &cell in &component.nodes {
                for i in 0..n {
                    let k = cell * n + i;
                    let direction = output[k];
                    if direction == 0.0 {
                        continue;
                    }
                    let radius = if direction < 0.0 {
                        self.center[k]
                    } else {
                        self.upper[k] - self.center[k]
                    };
                    if !positive(radius) {
                        return Err(OsmoticNetworkError::StepUnderflow);
                    }
                    gauge = gauge.max(finite(direction.abs() / radius)?);
                }
            }
            // Work with the normalized direction to avoid x*scale overflow.
            // Both weights are evaluated independently. At a depletion boundary
            // the small retained weight is not erased by subtracting from one.
            let inverse = 1.0 / scale;
            let denominator = inverse + gauge;
            let keep = inverse / denominator;
            let advance = gauge / denominator;
            for &cell in &component.nodes {
                for i in 0..n {
                    let k = cell * n + i;
                    let direction = output[k];
                    if direction == 0.0 || gauge == 0.0 {
                        output[k] = self.center[k];
                        continue;
                    }
                    let radius = if direction < 0.0 {
                        self.center[k]
                    } else {
                        self.upper[k] - self.center[k]
                    };
                    let fraction = (direction.abs() / radius) / gauge;
                    let boundary = if direction < 0.0 {
                        self.center[k] * (1.0 - fraction)
                    } else if fraction == 1.0 {
                        self.upper[k]
                    } else {
                        self.center[k] + radius * fraction
                    };
                    output[k] = finite(keep * self.center[k] + advance * boundary)?;
                    if output[k] == 0.0 && self.center[k] != 0.0 {
                        return Err(OsmoticNetworkError::StepUnderflow);
                    }
                }
            }
        }
        self.admit(output, cancelled)
    }

    fn admit_component(
        &self,
        component: &Component,
        amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        for &cell in &component.nodes {
            self.admit_cell(cell, amounts, cancelled)?;
        }
        Ok(())
    }

    fn admit(
        &self,
        amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        for cell in 0..self.volumes.len() {
            self.admit_cell(cell, amounts, cancelled)?;
        }
        Ok(())
    }

    fn admit_cell(
        &self,
        cell: usize,
        amounts: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), OsmoticNetworkError> {
        poll(cancelled)?;
        let n = self.species;
        for i in 0..n {
            let k = cell * n + i;
            if !nonnegative(amounts[k])
                || amounts[k] > self.upper[k]
                || amounts[k] / self.volumes[cell] > self.concentration_caps[k]
            {
                return Err(OsmoticNetworkError::Invalid(
                    "decoded ionic inventory leaves an original finite concentration box",
                ));
            }
        }
        self.ionic
            .admit_inventories(&amounts[cell * n..(cell + 1) * n])
            .map_err(diffusion_source)?;
        Ok(())
    }
}

fn allocation() -> OsmoticNetworkError {
    OsmoticNetworkError::Invalid("ionic coordinate payload allocation failed")
}
fn copied<T: Copy>(input: &[T]) -> Result<Vec<T>, OsmoticNetworkError> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(input.len())
        .map_err(|_| allocation())?;
    output.extend_from_slice(input);
    Ok(output)
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, OsmoticNetworkError> {
    let mut output = Vec::new();
    output.try_reserve_exact(count).map_err(|_| allocation())?;
    output.resize(count, value);
    Ok(output)
}

#[cfg(test)]
mod tests;
