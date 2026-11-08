//! One source-temperature and local Joule-heat path for settled AC consumers.

use super::{SpatialHarmonicHeatingCell, SpatialHarmonicHeatingStep};
use crate::circuit_enthalpy::{
    CircuitEnthalpyError, finite,
    magnetic::{
        Error, MagneticRoutes, MagneticThermalError, MagneticThermalNetwork, MagneticThermalSample,
        check, invalid, material,
    },
    material::{MaterialNetworkError, MaterialNetworkSample},
};
use fs_em::spatial::{
    SpatialMagneticMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, induction::harmonic::SpatialHarmonicInductionResponse,
    },
};
use fs_material::{
    conductor::spatial::SpatialConductorMaterial, field_storage::admit_source_temperature_refresh,
    state_point::ResolvedMaterialStatePoint,
};
use fs_thermal::enthalpy::DepositionSimulation;
use std::{cell::RefCell, collections::BTreeMap};

pub(super) struct HarmonicFieldMaterials<'p, C: SpatialConductorMaterial> {
    problem: &'p SpatialVectorMagnetostatics,
    originals: Vec<&'p ResolvedMaterialStatePoint>,
    conductors: BTreeMap<usize, &'p C>,
    volumes_m3: Vec<f64>,
}

impl<'p, C: SpatialConductorMaterial> HarmonicFieldMaterials<'p, C> {
    pub(super) fn new(
        problem: &'p SpatialVectorMagnetostatics,
        conductors: &'p [Option<C>],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        check(cancel)?;
        let originals = problem
            .materials()
            .iter()
            .map(|m| match m {
                SpatialMagneticMaterial::Linear(m) => Ok(m),
                SpatialMagneticMaterial::Recoil(_) => Err(invalid(
                    "settled AC heating requires original induced-linear magnetic materials",
                )),
            })
            .collect::<Result<Vec<_>, _>>()?;
        let volumes_m3 = (0..problem.mesh().tets.len())
            .map(|cell| {
                check(cancel)?;
                problem
                    .cell_volume_m3(cell)
                    .ok_or_else(|| invalid("original AC tetrahedron has no physical volume"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let conductors = conductors
            .iter()
            .enumerate()
            .filter_map(|(cell, source)| source.as_ref().map(|source| (cell, source)))
            .collect();
        Ok(Self {
            problem,
            originals,
            conductors,
            volumes_m3,
        })
    }

    pub(super) fn cell_count(&self) -> usize {
        self.originals.len()
    }

    pub(super) fn routes(
        &self,
        network: &MagneticThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MagneticRoutes, Error> {
        let conductors = self.conductors.keys().copied().collect::<Vec<_>>();
        network.routes_for_volumes(thermal, &self.volumes_m3, &conductors, cancel)
    }

    pub(super) fn resolve(
        &self,
        network: &MagneticThermalNetwork<'_>,
        state: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        routes: &MagneticRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(Vec<(usize, C)>, Vec<MagneticThermalSample>), Error> {
        network.materials_at_using(
            state,
            &self.originals,
            &self.conductors,
            sample,
            routes,
            cancel,
            |cell, before, after, binding| {
                let field_error =
                    |e| CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(e));
                admit_source_temperature_refresh(before, after).map_err(field_error)?;
                let law = self.problem.magnetic_laws()[cell];
                if law.resolve(before).map_err(field_error)?
                    != law.resolve(after).map_err(field_error)?
                {
                    return Err(invalid(
                        "settled AC heating cannot change magnetic storage laws",
                    ));
                }
                for name in law.source_property_names() {
                    let property = after.property(name).ok_or_else(|| {
                        invalid("settled AC heating lost its original magnetic source")
                    })?;
                    if property
                        .source_axes()
                        .iter()
                        .any(|(axis, _)| axis == binding.point.temperature_axis)
                    {
                        return Err(invalid(
                            "settled AC permeability must be independent of temperature",
                        ));
                    }
                    binding
                        .point
                        .admit_static(property)
                        .map_err(|e| material(MaterialNetworkError::Material(e)))?;
                }
                Ok(())
            },
        )
    }
}

pub(super) fn recover_field_heat<C: SpatialConductorMaterial>(
    response: &SpatialHarmonicInductionResponse<'_, '_, C>,
    dt: f64,
    routes: &MagneticRoutes,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<SpatialHarmonicHeatingStep, Error> {
    if response.cells().len() != routes.conductors.len() {
        return Err(invalid("AC heat omitted original conducting cells"));
    }
    let mut cells = Vec::with_capacity(response.cells().len());
    for loss in response.cells() {
        check(cancel)?;
        let cell = *routes
            .conductors
            .get(&loss.cell)
            .ok_or_else(|| invalid("AC heat lost its original caloric receiver"))?;
        let heat_j = finite(dt * loss.mean_ohmic_power_w)?;
        if heat_j < 0.0 || loss.mean_ohmic_power_w != 0.0 && heat_j == 0.0 {
            return Err(invalid("AC cell heat is negative or unrepresentable"));
        }
        cells.push(SpatialHarmonicHeatingCell {
            tetrahedron: loss.cell,
            cell,
            volume_m3: loss.volume_m3,
            mean_ohmic_power_w: loss.mean_ohmic_power_w,
            heat_j,
        });
    }
    let balance = *response.balance();
    Ok(SpatialHarmonicHeatingStep {
        omega_rad_s: response.omega_rad_s(),
        duration_s: dt,
        cells,
        terminals: response.terminals().to_vec(),
        balance,
        field: response.report().clone(),
        maximum_ampere_residual_a: response.maximum_ampere_residual_a(),
        maximum_charge_residual_a: response.maximum_charge_residual_a(),
        maximum_cell_skin_ratio: response.maximum_cell_skin_ratio(),
        heat_j: finite(dt * balance.mean_ohmic_power_w)?,
        supplied_work_j: finite(dt * balance.mean_supplied_power_w)?,
    })
}
