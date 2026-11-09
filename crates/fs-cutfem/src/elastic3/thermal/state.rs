//! Physical instantaneous thermoelastic states on original native coordinates.
//!
//! A dynamic endpoint need not be a static equilibrium. Recover the source-owned
//! resistance and preparations without solving away its inertial motion.

use std::ops::ControlFlow;

use fs_material::tensor::from_mandel;
use fs_solver::op::LinearOp;

use super::equilibrium::{EquilibriumOwner3, admit_selected_supports};
use super::{ThermalCoordinateMap3, ThermoelasticEquilibriumSettings3, checkpoint};
use crate::elastic3::eigenstrain::PreparedEigenstrain3;
use crate::elastic3::mass::MaterialMass3;
use crate::elastic3::{CutElasticity3, ElasticityError3, PreparedPrestress3};

/// Explicit bounds for an instantaneous physical state; no static solve is run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoelasticStateSettings3 {
    /// Maximum both original physical and independent displacement dimensions.
    pub maximum_native_dofs: usize,
    /// Actual and constitutive TOTAL-strain Frobenius ceiling in (0,1).
    /// Original source total/elastic component domains remain independent.
    pub maximum_kinematic_strain: f64,
    /// Original bulk, selected support and ghost point allowance.
    pub maximum_points: usize,
    /// Original upper-triangular scalar consistent-mass entries, before xyz.
    pub maximum_mass_terms: usize,
    /// Optional signed geometric stiffness; None skips its assembly.
    pub maximum_prestress_entries: Option<usize>,
}

impl From<ThermoelasticEquilibriumSettings3> for ThermoelasticStateSettings3 {
    fn from(settings: ThermoelasticEquilibriumSettings3) -> Self {
        Self {
            maximum_native_dofs: settings.maximum_native_dofs,
            maximum_kinematic_strain: settings.maximum_kinematic_strain,
            maximum_points: settings.maximum_points,
            maximum_mass_terms: settings.maximum_mass_terms,
            maximum_prestress_entries: settings.maximum_prestress_entries,
        }
    }
}

/// An admitted physical endpoint, which may have nonzero free resistance.
///
/// Displacement uses the original independent coordinates with exact zero strong
/// clamps. Resistance is full `K(T)u-b_free(T)`, including the original ghost and
/// homogeneous embedded-support terms and strongly supported rows. It excludes
/// external loads, inertia, damping and geometric stiffness. Supported entries
/// are NOT dynamic reactions: computing those also requires supported inertia.
///
/// Eigenstrain and mass retain the original geometry/material binding; optional
/// prestress additionally freezes the exact stiffness scales and displacement.
/// No equilibrium, stability or total energy balance is asserted by this packet.
#[derive(Debug)]
pub struct ThermoelasticState3 {
    pub displacement: Vec<f64>,
    pub resistance: Vec<f64>,
    pub eigenstrain: PreparedEigenstrain3,
    pub mass: MaterialMass3,
    pub prestress: Option<PreparedPrestress3>,
}

macro_rules! native {
    ($owner:expr, $method:ident $(, $argument:expr)* $(,)?) => {
        match $owner {
            EquilibriumOwner3::Cartesian(op) => op.$method($($argument),*),
            EquilibriumOwner3::Adaptive(op) => op.$method($($argument),*),
        }
    };
}

fn invalid(message: &'static str) -> ElasticityError3 {
    ElasticityError3::Invalid(message)
}

// An arbitrary same-owner eigenstrain has a valid geometry token but is not
// necessarily the thermal source's free strain. Check retained original samples,
// without material callbacks, quadrature reconstruction or force reassembly.
fn admit_source_eigenstrain(
    raw: &CutElasticity3,
    prepared: &PreparedEigenstrain3,
    poll: &mut impl FnMut() -> ControlFlow<()>,
) -> Result<(), ElasticityError3> {
    raw.validate_eigenstrain(prepared)?;
    if prepared.cells.len() != raw.cells.len() {
        return Err(invalid("instantaneous thermoelastic source/cell alignment mismatch"));
    }
    for (id, (cell, free)) in raw.cells.iter().zip(&prepared.cells).enumerate() {
        checkpoint(poll)?;
        let source = cell.material.thermal_state().ok_or(invalid(
            "instantaneous thermoelastic state requires a native thermal source",
        ))?.source();
        let expected = from_mandel(&source.free_strain_mandel());
        let surface_points = raw.embedded.as_ref()
            .map_or(0, |support| support.selected_points(id).len());
        if free.bulk_strains.len() != cell.rules.bulk().len()
            || free.surface_strains.len() != surface_points
            || free.mean_trace != expected[..3].iter().sum::<f64>()
        {
            return Err(invalid("retained eigenstrain is not the current native thermal source"));
        }
        for strain in free.bulk_strains.iter().chain(&free.surface_strains) {
            checkpoint(poll)?;
            if *strain != expected {
                return Err(invalid("retained eigenstrain is not the current native thermal source"));
            }
        }
    }
    checkpoint(poll)
}

/// Shared native implementation. Full scatter retains support rows; the input
/// map alone masks homogeneous clamps. Neither geometry nor constraints is copied.
pub(in crate::elastic3) fn recover_state(
    owner: EquilibriumOwner3<'_>,
    map: ThermalCoordinateMap3<'_>,
    u: &[f64],
    retained: Option<(PreparedEigenstrain3, MaterialMass3)>,
    settings: ThermoelasticStateSettings3,
    scatter: impl Fn(usize, usize, f64, &mut [f64]),
    mut poll: impl FnMut() -> ControlFlow<()>,
) -> Result<ThermoelasticState3, ElasticityError3> {
    checkpoint(&mut poll)?;
    let raw = match owner {
        EquilibriumOwner3::Cartesian(op) => op,
        EquilibriumOwner3::Adaptive(op) => &op.raw,
    };
    let n = map.n();
    if n == 0 || n > settings.maximum_native_dofs
        || raw.n() > settings.maximum_native_dofs || u.len() != n
        || !settings.maximum_kinematic_strain.is_finite()
        || settings.maximum_kinematic_strain <= 0.0
        || settings.maximum_kinematic_strain >= 1.0
    {
        return Err(invalid("invalid instantaneous thermoelastic state bounds or displacement"));
    }
    for chunk in u.chunks(192) {
        checkpoint(&mut poll)?;
        if chunk.iter().any(|value| !value.is_finite()) {
            return Err(invalid("nonfinite instantaneous thermoelastic displacement"));
        }
    }
    // Admit source/point work and actual/assumed physical strains before any
    // state output or new eigenstrain/mass/prestress allocation.
    raw.admit_thermal_temperature_points(settings.maximum_points, &mut poll)?;
    native!(owner, admit_thermoelastic, u, settings.maximum_kinematic_strain,
        settings.maximum_points, &mut poll)?;
    let (eigenstrain, mass) = match retained {
        Some((eigenstrain, mass)) => {
            admit_source_eigenstrain(raw, &eigenstrain, &mut poll)?;
            native!(owner, validate_material_mass, &mass)?;
            if mass.scalar_terms() > settings.maximum_mass_terms {
                return Err(invalid("retained thermoelastic mass exceeds the scalar allowance"));
            }
            (eigenstrain, mass)
        }
        None => (
            native!(owner, prepare_thermoelastic, settings.maximum_points, &mut poll)?,
            native!(owner, prepare_material_mass, settings.maximum_mass_terms, &mut poll)?,
        ),
    };
    // Surface admission is required even when geometric stiffness is omitted.
    let prestress = match settings.maximum_prestress_entries {
        Some(entries) => Some(native!(owner, prepare_eigenstrain_prestress,
            &eigenstrain, u, settings.maximum_kinematic_strain,
            settings.maximum_points, entries, &mut poll)?),
        None => {
            admit_selected_supports(raw, map, &eigenstrain, u,
                settings.maximum_kinematic_strain, 0, &mut |_| poll())?;
            None
        }
    };
    let mut displacement = Vec::with_capacity(n);
    for (node, values) in u.chunks_exact(3).enumerate() {
        if node % 64 == 0 { checkpoint(&mut poll)?; }
        if map.fixed[node] {
            displacement.extend_from_slice(&[0.0; 3]);
        } else {
            displacement.extend_from_slice(values);
        }
    }
    // Reuse original stiffness/ghost matrices and the existing retained thermal
    // force columns. Subtract in-place; no second global load vector is needed.
    let mut resistance = raw.uneliminated_force(n,
        |node, component| map.value(&displacement, node, component),
        &scatter, &mut poll)?;
    raw.add_eigenstrain_full_load(&eigenstrain, &mut resistance,
        |node, component, value, output| scatter(node, component, -value, output),
        &mut poll)?;
    checkpoint(&mut poll)?;
    Ok(ThermoelasticState3 { displacement, resistance, eigenstrain, mass, prestress })
}

impl CutElasticity3 {
    /// Admit an instantaneous native thermoelastic displacement, retaining the
    /// original source free strain, consistent reference mass and optional G.
    /// A nonzero force imbalance is valid; there is no static solve or force gate.
    pub fn prepare_thermoelastic_state(
        &self,
        u: &[f64],
        settings: ThermoelasticStateSettings3,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<ThermoelasticState3, ElasticityError3> {
        recover_state(EquilibriumOwner3::Cartesian(self),
            ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            u, None, settings,
            |node, component, value, output| output[3 * node + component] += value,
            poll)
    }

    /// Recover a new displacement using consumed current native preparations.
    ///
    /// No mass integration or eigenstrain reassembly occurs. Original source
    /// samples, material/geometry tokens and retained mass limits are rechecked;
    /// arbitrary same-owner free strain is refused. G is freshly recovered at u.
    /// Errors/cancellation publish nothing and never mutate this native operator.
    pub fn recover_thermoelastic_state(
        &self,
        u: &[f64],
        eigenstrain: PreparedEigenstrain3,
        mass: MaterialMass3,
        settings: ThermoelasticStateSettings3,
        poll: impl FnMut() -> ControlFlow<()>,
    ) -> Result<ThermoelasticState3, ElasticityError3> {
        recover_state(EquilibriumOwner3::Cartesian(self),
            ThermalCoordinateMap3 { rows: None, fixed: &self.fixed },
            u, Some((eigenstrain, mass)), settings,
            |node, component, value, output| output[3 * node + component] += value,
            poll)
    }
}
