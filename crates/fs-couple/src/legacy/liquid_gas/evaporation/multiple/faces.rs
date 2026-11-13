//! Distinct native faces share one original liquid inventory and caloric state.

use super::{
    BTreeMap, CancelGate, Error, GasFilmHeatSolve, LiquidGasEvaporationSurface,
    MaterialFreeSurfaceEnthalpy3, Total, poll,
};
use fs_thermal::evaporation::EvaporationPatch;

pub(super) struct CellFaces<'a> {
    /// Canonical first face retains the original material and binding identity.
    pub(super) surface: LiquidGasEvaporationSurface<'a>,
    pub(super) patches: Vec<EvaporationPatch>,
    heat: Option<(f64, GasFilmHeatSolve)>,
    heat_conductance: Total,
}

impl CellFaces<'_> {
    /// The parallel phase owner supplies the authoritative summed physical area.
    /// All faces see the same original liquid and finite gas temperatures, so
    /// sensible heat adds as conductance while nonlinear phase films stay apart.
    pub(super) fn heat(&self, area_m2: f64) -> Result<Option<(f64, GasFilmHeatSolve)>, Error> {
        let Some((first, controls)) = self.heat else {
            return Ok(None);
        };
        if self.patches.len() == 1 {
            return Ok(Some((first, controls)));
        }
        let conductance = self.heat_conductance.value();
        let coefficient = if area_m2 == 0.0 {
            0.0
        } else {
            conductance / area_m2
        };
        if !coefficient.is_finite()
            || coefficient < 0.0
            || (conductance > 0.0 && coefficient == 0.0)
        {
            return Err(Error::Invalid(
                "combined native sensible-film coefficient is unrepresentable",
            ));
        }
        Ok(Some((coefficient, controls)))
    }
}

pub(super) fn prepare<'a>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    surfaces: impl ExactSizeIterator<
        Item = (
            LiquidGasEvaporationSurface<'a>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    with_heat: bool,
    gate: &CancelGate,
) -> Result<BTreeMap<[usize; 3], CellFaces<'a>>, Error> {
    let mut faces = BTreeMap::new();
    let mut binding = None;
    for (surface, heat) in surfaces {
        poll(gate)?;
        let common = *binding.get_or_insert(surface.binding);
        if surface.binding != common || heat.is_some() != with_heat {
            return Err(Error::Invalid(
                "native faces need one volatile/datum binding and a common heat mode",
            ));
        }
        if !surface.mass_transfer_coefficient_m_s.is_finite()
            || surface.mass_transfer_coefficient_m_s < 0.0
            || heat.is_some_and(|(coefficient, _)| !coefficient.is_finite() || coefficient < 0.0)
        {
            return Err(Error::Invalid(
                "native face transfer coefficients must be finite and nonnegative",
            ));
        }
        // Admit each measured face before any reduction. A corner's total
        // exposure may exceed one dx^2 face, but no individual patch may do so.
        liquid.admit_interface_gas_patch(surface.liquid, surface.gas, surface.area_m2)?;
        if faces
            .insert((surface.liquid, surface.gas), (surface, heat))
            .is_some()
        {
            return Err(Error::Invalid(
                "native phase exposure repeats a liquid/gas face",
            ));
        }
    }
    let mut cells = BTreeMap::<[usize; 3], CellFaces<'a>>::new();
    for ((coordinates, _), (surface, heat)) in faces {
        poll(gate)?;
        let cell = cells.entry(coordinates).or_insert_with(|| CellFaces {
            surface,
            patches: Vec::new(),
            heat,
            heat_conductance: Total::default(),
        });
        if !core::ptr::eq(cell.surface.law, surface.law)
            || cell.heat.map(|(_, controls)| controls) != heat.map(|(_, controls)| controls)
        {
            return Err(Error::Invalid(
                "faces of one native cell require the same original phase law and heat controls",
            ));
        }
        cell.patches.push(EvaporationPatch {
            area_m2: surface.area_m2,
            mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
        });
        if let Some((coefficient, _)) = heat {
            let conductance = surface.area_m2 * coefficient;
            if surface.area_m2 > 0.0 && coefficient > 0.0 && conductance == 0.0 {
                return Err(Error::Invalid(
                    "positive native sensible-film conductance underflows",
                ));
            }
            cell.heat_conductance.add(conductance)?;
        }
    }
    poll(gate)?;
    Ok(cells)
}
