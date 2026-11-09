//! Exact original caloric/material/cut-cell binding; no geometry reconstruction.
use std::collections::{BTreeMap, BTreeSet};

use fs_cutfem::HexCell;
use fs_material::elastic::thermal::{
    AffineThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_thermal::enthalpy::{DepositedCell, ThermalLaw};

use super::*;

/// One original native active cell and its actual caloric/material source.
/// The order is exactly the native operator's active-cell order. The law owns
/// source selection, reference expansion and its complete temperature interval.
pub struct CutThermoelasticCell3<'a> {
    pub cell: CellId,
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub point: ThermalMaterialPoint<'a>,
    pub law: AffineThermoelastic,
}

pub(super) struct BoundCell<'a> {
    pub cell: CellId,
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub point: ThermalMaterialPoint<'a>,
    pub reference: DepositedCell,
    input: IsothermalThermoelasticInput,
    bounds: HexCell,
    volume_m3: f64,
}
pub(super) struct Bindings<'a> {
    pub cells: Vec<BoundCell<'a>>,
    pub laws: Vec<AffineThermoelastic>,
    by_bounds: BTreeMap<[u64; 6], usize>,
    thermal_ids: Vec<CellId>,
    thermal_links: Vec<u64>,
}

fn key(bounds: HexCell) -> [u64; 6] {
    [bounds.lo()[0].to_bits(), bounds.lo()[1].to_bits(), bounds.lo()[2].to_bits(),
        bounds.hi()[0].to_bits(), bounds.hi()[1].to_bits(), bounds.hi()[2].to_bits()]
}
fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a - b).abs() <= 512.0 * f64::EPSILON * a.abs().max(b.abs())
}

impl<'a> Bindings<'a> {
    pub fn new<E>(
        native: CutElasticityRef3<'_>, thermal: &DepositionSimulation<'_>,
        sources: Vec<CutThermoelasticCell3<'a>>, maximum_cells: usize, cx: &Cx<'_>,
    ) -> Result<Self, CutThermoelasticError3<E>> {
        poll(cx)?;
        if sources.is_empty() || sources.len() != native::cells(native)
            || sources.len() > maximum_cells || native::has_embedded_supports(native)
        {
            return Err(CutThermoelasticError3::Invalid("one bounded original source per active cell and homogeneous strong supports required"));
        }
        let masses = native::masses(native).ok_or(CutThermoelasticError3::Invalid(
            "original source-density mass inventory required"))?;
        let mut cells = Vec::with_capacity(sources.len());
        let mut laws = Vec::with_capacity(sources.len());
        let mut by_bounds = BTreeMap::new();
        let mut used = BTreeSet::new();
        for (index, source) in sources.into_iter().enumerate() {
            poll(cx)?;
            if !used.insert(source.cell) || source.point.temperature_axis != "T" {
                return Err(CutThermoelasticError3::Invalid("unique active caloric cells and original canonical T source axis required"));
            }
            let reference = thermal.cells().get(&source.cell)
                .ok_or(CutThermoelasticError3::Invalid("inactive cut-solid caloric cell"))?;
            if !thermal.uses_material_law(reference.material_id, source.curve as &dyn ThermalLaw) {
                return Err(CutThermoelasticError3::Invalid("cut solid must retain its actual borrowed caloric chart"));
            }
            let phase = source.curve.state_at_specific_enthalpy(reference.specific_h_j_kg)?;
            let held = native::thermal(native, index)
                .ok_or(CutThermoelasticError3::Invalid("complete native thermoelastic material required in every cell"))?;
            let input = *held.input();
            let bounds = native::bounds(native, index)
                .ok_or(CutThermoelasticError3::Invalid("missing native reference-cell bounds"))?;
            let volume = native::volume(native, index)
                .ok_or(CutThermoelasticError3::Invalid("missing native cut reference volume"))?;
            if !phase.is_fully_solid() || !same_measure(reference.geometry.volume_m3, volume)
                || !same_measure(reference.mass_kg, masses[index])
                || !same_measure(reference.mass_kg, input.reference_density_kg_m3 * volume)
            {
                return Err(CutThermoelasticError3::Invalid("fully solid caloric inventory must equal original cut volume and source-density mass"));
            }
            let resolved = source.law.at_temperature(phase.temperature_k())?;
            let query = source.point.at(phase)?;
            if resolved.elastic.card_identity() != source.point.material.content_hash()
                || resolved.expansion.card_identity() != source.point.material.content_hash()
                || resolved.elastic.query() != query || resolved.expansion.query() != query
                || resolved.isothermal_elasticity(input)? != *held
            {
                return Err(CutThermoelasticError3::Invalid("native stiffness, expansion and thermal source must be the same original material experiment"));
            }
            for property in held.source().elastic.properties().iter()
                .chain(held.source().expansion.properties())
            {
                source.point.admit_static(property)?;
            }
            admit_capacity(&source, reference.mass_kg, volume,
                input.reference_density_kg_m3, native::scales(native)[index], cx)?;
            if by_bounds.insert(key(bounds), index).is_some() {
                return Err(CutThermoelasticError3::Invalid("duplicate native reference-cell bounds"));
            }
            cells.push(BoundCell { cell: source.cell, curve: source.curve, point: source.point,
                reference: reference.clone(), input, bounds, volume_m3: volume });
            laws.push(source.law);
        }
        let thermal_ids = thermal.cells().keys().copied().collect();
        let thermal_links = thermal.links().keys().copied().collect();
        poll(cx)?;
        Ok(Self { cells, laws, by_bounds, thermal_ids, thermal_links })
    }

    pub fn admit_native<E>(&self, native: CutElasticityRef3<'_>, cx: &Cx<'_>)
        -> Result<(), CutThermoelasticError3<E>>
    {
        if native::cells(native) != self.cells.len() || native::has_embedded_supports(native) {
            return Err(CutThermoelasticError3::Invalid("native cell topology or supports changed"));
        }
        let masses = native::masses(native).ok_or(CutThermoelasticError3::Invalid(
            "native reference density inventory disappeared"))?;
        for (i, cell) in self.cells.iter().enumerate() {
            if i % 64 == 0 { poll(cx)?; }
            if native::bounds(native, i).map(key) != Some(key(cell.bounds))
                || native::volume(native, i).map(f64::to_bits) != Some(cell.volume_m3.to_bits())
                || !same_measure(masses[i], cell.reference.mass_kg)
                || native::thermal(native, i).map(|s| *s.input()) != Some(cell.input)
            {
                return Err(CutThermoelasticError3::Invalid("native reference geometry, original mass or material inventory changed"));
            }
        }
        Ok(())
    }

    fn inventory<E>(&self, thermal: &DepositionSimulation<'_>)
        -> Result<(), CutThermoelasticError3<E>>
    {
        if thermal.cells().len() != self.thermal_ids.len()
            || !thermal.cells().keys().copied().eq(self.thermal_ids.iter().copied())
            || thermal.links().len() != self.thermal_links.len()
            || !thermal.links().keys().copied().eq(self.thermal_links.iter().copied())
        {
            return Err(CutThermoelasticError3::Invalid("fixed thermal inventory and contact topology required"));
        }
        Ok(())
    }

    fn phase<E>(&self, index: usize, thermal: &DepositionSimulation<'_>)
        -> Result<(EquilibriumPhaseState, f64), CutThermoelasticError3<E>>
    {
        let bound = &self.cells[index];
        let cell = thermal.cells().get(&bound.cell)
            .ok_or(CutThermoelasticError3::Invalid("inactive original cut-solid caloric cell"))?;
        let old = &bound.reference;
        if cell.id != old.id || cell.material_id != old.material_id
            || cell.segment_id != old.segment_id || cell.body_id != old.body_id
            || cell.geometry != old.geometry || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
            || !thermal.uses_material_law(cell.material_id, bound.curve as &dyn ThermalLaw)
        {
            return Err(CutThermoelasticError3::Invalid("cut-solid reference cell, material or caloric owner changed"));
        }
        let phase = bound.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
        let range = self.laws[index].temperature_range_k();
        if !phase.is_fully_solid() || phase.temperature_k() < range[0]
            || phase.temperature_k() > range[1]
        {
            return Err(CutThermoelasticError3::Invalid("solid thermal state left its source interval or requires a phase-support handoff"));
        }
        Ok((phase, cell.specific_h_j_kg))
    }

    /// Cheap substep admission; no repeated stiffness assembly or source solve.
    pub fn admit_heat<E>(&self, thermal: &DepositionSimulation<'_>, cx: &Cx<'_>)
        -> Result<(), CutThermoelasticError3<E>>
    {
        self.inventory(thermal)?;
        for i in 0..self.cells.len() { poll(cx)?; self.phase::<E>(i, thermal)?; }
        Ok(())
    }

    pub fn admit<E>(&self, thermal: &DepositionSimulation<'_>, cx: &Cx<'_>)
        -> Result<(Vec<f64>, Vec<f64>), CutThermoelasticError3<E>>
    {
        self.inventory(thermal)?;
        let mut temperatures = Vec::with_capacity(self.cells.len());
        let mut enthalpies = Vec::with_capacity(self.cells.len());
        for i in 0..self.cells.len() {
            poll(cx)?;
            let (phase, h) = self.phase(i, thermal)?;
            temperatures.push(phase.temperature_k());
            enthalpies.push(h);
        }
        Ok((temperatures, enthalpies))
    }

    pub fn materials<E>(&self, temperatures: &[f64], cx: &Cx<'_>)
        -> Result<Vec<IsothermalThermoelasticState>, CutThermoelasticError3<E>>
    {
        let mut result = Vec::with_capacity(self.cells.len());
        for (i, bound) in self.cells.iter().enumerate() {
            poll(cx)?;
            result.push(self.laws[i].at_temperature(temperatures[i])?
                .isothermal_elasticity(bound.input)?);
        }
        Ok(result)
    }

    pub fn material_at(
        &self, bounds: HexCell, materials: &[IsothermalThermoelasticState],
    ) -> Result<IsothermalThermoelasticState, ElasticityError3> {
        self.by_bounds.get(&key(bounds)).and_then(|i| materials.get(*i)).cloned()
            .ok_or(ElasticityError3::Invalid("candidate material lacks its original native cell"))
    }
}

fn admit_capacity<E>(
    source: &CutThermoelasticCell3<'_>, mass: f64, volume: f64, density: f64,
    scale: f64, cx: &Cx<'_>,
) -> Result<(), CutThermoelasticError3<E>> {
    let range = source.law.temperature_range_k();
    let knots = source.curve.knots();
    if knots.is_empty() || range[0] < knots[0].temperature_k
        || range[1] > knots[knots.len()-1].temperature_k
    {
        return Err(CutThermoelasticError3::Invalid("original caloric chart must cover the complete thermoelastic interval"));
    }
    let mut minimum = f64::INFINITY;
    for (interval, pair) in knots.windows(2).enumerate() {
        poll(cx)?;
        let dt = pair[1].temperature_k - pair[0].temperature_k;
        if dt <= 0.0 {
            if pair[0].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                return Err(CutThermoelasticError3::Invalid("solid thermal interval intersects an isothermal phase segment"));
            }
            continue;
        }
        if pair[0].temperature_k >= range[1] || pair[1].temperature_k <= range[0] { continue; }
        if pair.iter().any(|k| k.liquid_mass_fraction != 0.0
            || !same_measure(k.bulk_density_kg_m3, density))
        {
            return Err(CutThermoelasticError3::Invalid("fully solid constant reference-density caloric interval required"));
        }
        minimum = minimum.min(finite(mass
            * source.curve.sensible_heat_capacity_bounds_j_kg_k(interval)?[0])?);
    }
    if !minimum.is_finite()
        || finite(minimum + scale * source.law.minimum_heat_capacity_j_k(volume)?)? <= 0.0
    {
        return Err(CutThermoelasticError3::Invalid("caloric plus elastic-field capacity must stay positive over the full original source envelope"));
    }
    Ok(())
}
