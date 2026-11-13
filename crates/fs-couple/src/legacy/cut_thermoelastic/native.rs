//! Borrowed dispatch to the two original native displacement spaces.
use std::ops::ControlFlow;
use fs_cutfem::HexCell;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::{CutElasticity3, MaterialRefreshLimits3};
use fs_material::elastic::thermal::{AffineThermoelastic, IsothermalThermoelasticState};
use fs_solver::LinearOp;
use super::*;

/// Exclusive view of an existing native operator; no geometry is cloned.
pub enum CutElasticityMut3<'a> {
    Cartesian(&'a mut CutElasticity3),
    Adaptive(&'a mut AdaptiveElasticity3),
}
impl CutElasticityMut3<'_> {
    pub fn as_ref(&self) -> CutElasticityRef3<'_> {
        match self {
            Self::Cartesian(op) => CutElasticityRef3::Cartesian(op),
            Self::Adaptive(op) => CutElasticityRef3::Adaptive(op),
        }
    }
    pub(super) fn refresh<R, E: From<ElasticityError3>>(
        &mut self, field: ElasticMaterialField3<'_>, limits: MaterialRefreshLimits3,
        attempt: impl FnOnce(CutElasticityRef3<'_>) -> Result<MaterialRefreshDecision3<R>, E>,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<MaterialRefreshDecision3<R>, E> {
        match self {
            Self::Cartesian(op) => op.refresh_material_field_with_decision(field, limits,
                |candidate| attempt(CutElasticityRef3::Cartesian(candidate)), checkpoint),
            Self::Adaptive(op) => op.refresh_material_field_with_decision(field, limits,
                |candidate| attempt(CutElasticityRef3::Adaptive(candidate)), checkpoint),
        }
    }
}
macro_rules! native {
    ($owner:expr, $method:ident $(, $arg:expr)* $(,)?) => {
        match $owner {
            CutElasticityRef3::Cartesian(op) => op.$method($($arg),*),
            CutElasticityRef3::Adaptive(op) => op.$method($($arg),*),
        }
    };
}
pub(super) fn n(owner: CutElasticityRef3<'_>) -> usize { native!(owner, n) }
pub(super) fn cells(owner: CutElasticityRef3<'_>) -> usize { native!(owner, cells) }
pub(super) fn fixed(owner: CutElasticityRef3<'_>) -> &[bool] { native!(owner, fixed) }
pub(super) fn scales(owner: CutElasticityRef3<'_>) -> &[f64] { native!(owner, scales) }
pub(super) fn masses(owner: CutElasticityRef3<'_>) -> Option<&[f64]> {
    native!(owner, reference_cell_masses_kg)
}
pub(super) fn thermal(owner: CutElasticityRef3<'_>, cell: usize)
    -> Option<&IsothermalThermoelasticState>
{
    native!(owner, cell_thermoelastic_state, cell)
}
pub(super) fn bounds(owner: CutElasticityRef3<'_>, cell: usize) -> Option<HexCell> {
    native!(owner, reference_cell_bounds, cell)
}
pub(super) fn volume(owner: CutElasticityRef3<'_>, cell: usize) -> Option<f64> {
    native!(owner, reference_cell_volume_m3, cell)
}
pub(super) fn has_embedded_supports(owner: CutElasticityRef3<'_>) -> bool {
    native!(owner, has_embedded_supports)
}
pub(super) fn validate(owner: CutElasticityRef3<'_>, state: &ThermoelasticEquilibrium3)
    -> Result<(), ElasticityError3>
{
    native!(owner, validate_eigenstrain, &state.eigenstrain)?;
    native!(owner, validate_material_mass, &state.mass)
}
pub(super) fn solve(
    owner: CutElasticityRef3<'_>, forces: &[f64], settings: ThermoelasticEquilibriumSettings3,
    cx: &Cx<'_>,
) -> Result<ThermoelasticEquilibrium3, ElasticityError3> {
    native!(owner, solve_thermoelastic_equilibrium,
        ReferenceLoad3::default(), Some(forces), None, settings, |_| checkpoint(cx))
}
#[allow(clippy::too_many_arguments)]
pub(super) fn exchange(
    owner: CutElasticityRef3<'_>, laws: &[AffineThermoelastic], before_temperatures: &[f64],
    before_u: &[f64], after_u: &[f64], maximum_kinematic_strain: f64,
    limits: ThermoelasticExchangeLimits3, cx: &Cx<'_>,
) -> Result<Vec<ThermoelasticCellExchange3>, ElasticityError3> {
    native!(owner, thermoelastic_exchange, laws, before_temperatures, before_u, after_u,
        maximum_kinematic_strain, limits, || checkpoint(cx))
}

pub(super) fn thermal_parts(
    owner: CutElasticityRef3<'_>, settings: ThermoelasticEquilibriumSettings3, cx: &Cx<'_>,
) -> Result<(fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3,
    fs_cutfem::elastic3::mass::MaterialMass3), ElasticityError3> {
    let eigenstrain = native!(owner, prepare_thermoelastic,
        settings.maximum_points, || checkpoint(cx))?;
    let mass = native!(owner, prepare_material_mass,
        settings.maximum_mass_terms, || checkpoint(cx))?;
    Ok((eigenstrain, mass))
}

pub(super) fn recover_state(
    owner: CutElasticityRef3<'_>, displacement: &[f64],
    eigenstrain: fs_cutfem::elastic3::eigenstrain::PreparedEigenstrain3,
    mass: fs_cutfem::elastic3::mass::MaterialMass3,
    settings: ThermoelasticEquilibriumSettings3, cx: &Cx<'_>,
) -> Result<fs_cutfem::elastic3::thermal::ThermoelasticState3, ElasticityError3> {
    native!(owner, recover_thermoelastic_state, displacement, eigenstrain, mass,
        settings.into(), || checkpoint(cx))
}

pub(super) fn validate_state(
    owner: CutElasticityRef3<'_>,
    state: &fs_cutfem::elastic3::thermal::ThermoelasticState3,
) -> Result<(), ElasticityError3> {
    native!(owner, validate_eigenstrain, &state.eigenstrain)?;
    native!(owner, validate_material_mass, &state.mass)
}
