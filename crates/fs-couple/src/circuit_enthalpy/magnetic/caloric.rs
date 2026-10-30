//! Reversible magnetic heat shared by retained-field and circuit-driven owners.
//!
//! Native field solvers supply their own Helmholtz-energy and parameter-work
//! receipts. The original material law supplies entropy and internal energy;
//! this owner routes Q_rev = -P - delta(T*S) into the original caloric cells.
//! It owns no field step, circuit step, constitutive law or thermal integrator.

use super::{
    Error, MagneticRoutes, MagneticThermalError, MagneticThermalNetwork, MagneticThermalSample,
    check, invalid, material, same_material,
};
use crate::circuit_enthalpy::{CircuitEnthalpyError, finite, material::MaterialNetworkError};
use fs_material::{
    field_storage::{
        FieldStorageError, MAGNETIC_PERMEABILITY_PROPERTY,
        nonlinear::{
            REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            thermal::{
                ResolvedThermalReversibleMagneticMaterial, ThermalReversibleMagneticMaterial,
            },
        },
        thermal::{
            FieldThermalState,
            magnetic::{ResolvedThermalMagneticMaterial, ThermalMagneticMaterial},
        },
    },
    state_point::ResolvedMaterialStatePoint,
};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, ValueForm},
};
use fs_rep_mesh::TriComplex2;
use fs_thermal::enthalpy::{CellId, DepositionSimulation};
use std::{cell::RefCell, collections::BTreeMap};

/// One original magnetic triangle's complete source and admitted field envelope.
/// Its receiving cell comes from the shared magnetic thermal-region map.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticThermalLaw<L = ThermalMagneticMaterial> {
    pub triangle: usize,
    pub law: L,
    pub maximum_flux_density_t: f64,
}

/// Supported source laws can coexist in one nonlinear field/heat assembly:
/// e.g. a source-separable magnetic core and ordinary linear winding/air cells.
/// Both variants use the same finite caloric receiver and original source pins.
#[derive(Clone, Debug, PartialEq)]
pub enum MagneticConstitutiveThermalLaw {
    Linear(ThermalMagneticMaterial),
    Separable(ThermalReversibleMagneticMaterial),
}

// Both original constitutive owners provide their own source/energy functions.
// This private dispatch keeps one capacity admission, entropy transfer and
// source-continuation owner; it cannot manufacture an independent thermal law.
pub(in crate::circuit_enthalpy) trait MagneticCaloricLaw {
    type State: MagneticCaloricState;
    fn card(&self) -> &fs_matdb::MaterialCard;
    fn temperature_axis(&self) -> &str;
    fn temperature_range_k(&self) -> [f64; 2];
    fn at_temperature(&self, temperature_k: f64) -> Result<Self::State, FieldStorageError>;
    fn minimum_heat_capacity_j_k(&self, b_t: f64, volume_m3: f64)
    -> Result<f64, FieldStorageError>;
}
pub(in crate::circuit_enthalpy) trait MagneticCaloricState {
    fn material(&self) -> &ResolvedMaterialStatePoint;
    fn temperature_k(&self) -> f64;
    fn state(&self, b_t: f64, volume_m3: f64) -> Result<FieldThermalState, FieldStorageError>;
}
macro_rules! original_thermal_law {
    ($law:ty, $state:ty) => {
        impl MagneticCaloricLaw for $law {
            type State = $state;
            fn card(&self) -> &fs_matdb::MaterialCard {
                <$law>::card(self)
            }
            fn temperature_axis(&self) -> &str {
                <$law>::temperature_axis(self)
            }
            fn temperature_range_k(&self) -> [f64; 2] {
                <$law>::temperature_range_k(self)
            }
            fn at_temperature(&self, t: f64) -> Result<Self::State, FieldStorageError> {
                <$law>::at_temperature(self, t)
            }
            fn minimum_heat_capacity_j_k(&self, b: f64, v: f64) -> Result<f64, FieldStorageError> {
                <$law>::minimum_heat_capacity_j_k(self, b, v)
            }
        }
        impl MagneticCaloricState for $state {
            fn material(&self) -> &ResolvedMaterialStatePoint {
                <$state>::material(self)
            }
            fn temperature_k(&self) -> f64 {
                <$state>::temperature_k(self)
            }
            fn state(&self, b: f64, v: f64) -> Result<FieldThermalState, FieldStorageError> {
                <$state>::state(self, b, v)
            }
        }
    };
}
original_thermal_law!(ThermalMagneticMaterial, ResolvedThermalMagneticMaterial);
original_thermal_law!(
    ThermalReversibleMagneticMaterial,
    ResolvedThermalReversibleMagneticMaterial
);

pub(in crate::circuit_enthalpy) enum ResolvedMagneticConstitutiveThermalLaw {
    Linear(ResolvedThermalMagneticMaterial),
    Separable(ResolvedThermalReversibleMagneticMaterial),
}
impl MagneticCaloricLaw for MagneticConstitutiveThermalLaw {
    type State = ResolvedMagneticConstitutiveThermalLaw;
    fn card(&self) -> &fs_matdb::MaterialCard {
        match self {
            Self::Linear(law) => law.card(),
            Self::Separable(law) => law.card(),
        }
    }
    fn temperature_axis(&self) -> &str {
        match self {
            Self::Linear(law) => law.temperature_axis(),
            Self::Separable(law) => law.temperature_axis(),
        }
    }
    fn temperature_range_k(&self) -> [f64; 2] {
        match self {
            Self::Linear(law) => law.temperature_range_k(),
            Self::Separable(law) => law.temperature_range_k(),
        }
    }
    fn at_temperature(&self, t: f64) -> Result<Self::State, FieldStorageError> {
        match self {
            Self::Linear(law) => law
                .at_temperature(t)
                .map(ResolvedMagneticConstitutiveThermalLaw::Linear),
            Self::Separable(law) => law
                .at_temperature(t)
                .map(ResolvedMagneticConstitutiveThermalLaw::Separable),
        }
    }
    fn minimum_heat_capacity_j_k(&self, b: f64, v: f64) -> Result<f64, FieldStorageError> {
        match self {
            Self::Linear(law) => law.minimum_heat_capacity_j_k(b, v),
            Self::Separable(law) => law.minimum_heat_capacity_j_k(b, v),
        }
    }
}
impl MagneticCaloricState for ResolvedMagneticConstitutiveThermalLaw {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Linear(at) => at.material(),
            Self::Separable(at) => at.material(),
        }
    }
    fn temperature_k(&self) -> f64 {
        match self {
            Self::Linear(at) => at.temperature_k(),
            Self::Separable(at) => at.temperature_k(),
        }
    }
    fn state(&self, b: f64, v: f64) -> Result<FieldThermalState, FieldStorageError> {
        match self {
            Self::Linear(at) => at.state(b, v),
            Self::Separable(at) => at.state(b, v),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MagneticCaloricSample {
    pub triangle: usize,
    pub cell: CellId,
    pub temperature_k: [f64; 2],
    pub flux_density_t: [f64; 2],
    pub state: [FieldThermalState; 2],
    pub parameter_work_j: f64,
    pub reversible_heat_j: f64,
}

/// Original source evaluations and reversible heat over one proposed interval.
#[derive(Clone, Debug, PartialEq)]
pub struct MagneticCaloricExchange {
    pub fields: Vec<MagneticCaloricSample>,
    /// Includes temperature-independent fixed regions with U=psi.
    pub field_internal_energy_j: [f64; 2],
    pub reversible_heat_j: f64,
}

/// A view of the ORIGINAL field owner's per-triangle energy receipt.
/// The ordered endpoints and signed P retain that owner's time discretization.
pub(in crate::circuit_enthalpy) struct MagneticEnergyCell<const D: usize = 2> {
    pub triangle: usize,
    pub flux_density_t: [[f64; D]; 2],
    pub free_energy_j: [f64; 2],
    pub parameter_work_j: f64,
}

struct Initial {
    cell: CellId,
    temperature_k: f64,
    b_t: f64,
    volume_m3: f64,
    state: FieldThermalState,
}

pub(in crate::circuit_enthalpy) struct PreparedMagneticCaloric<'a, L = ThermalMagneticMaterial> {
    laws: BTreeMap<usize, &'a MagneticThermalLaw<L>>,
    initial: BTreeMap<usize, Initial>,
    fixed: BTreeMap<usize, ResolvedMaterialStatePoint>,
    triangle_count: usize,
}

fn field_error(error: FieldStorageError) -> Error {
    CircuitEnthalpyError::Source(MagneticThermalError::Magnetic(error))
}
fn magnitude<const D: usize>(b: [f64; D]) -> f64 {
    b.into_iter().fold(0.0_f64, f64::hypot)
}

impl<'a, L: MagneticCaloricLaw> PreparedMagneticCaloric<'a, L> {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn new(
        network: &MagneticThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        mesh: &TriComplex2,
        materials: &[&ResolvedMaterialStatePoint],
        flux_density_t: &[[f64; 2]],
        bindings: &'a [MagneticThermalLaw<L>],
        routes: &MagneticRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        check(cancel)?;
        let volumes = (0..mesh.faces().len())
            .map(|cell| {
                check(cancel)?;
                mesh.face_measure(cell)
                    .ok_or_else(|| invalid("missing original magnetic volume"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::new_for_volumes(
            network,
            thermal,
            &volumes,
            materials,
            flux_density_t,
            bindings,
            routes,
            cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn new_for_volumes<const D: usize>(
        network: &MagneticThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        volumes_m3: &[f64],
        materials: &[&ResolvedMaterialStatePoint],
        flux_density_t: &[[f64; D]],
        bindings: &'a [MagneticThermalLaw<L>],
        routes: &MagneticRoutes,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        check(cancel)?;
        if D != 2 && D != 3 {
            return Err(invalid(
                "magnetic caloric owner requires planar or spatial vectors",
            ));
        }
        let triangle_count = volumes_m3.len();
        let expected = routes.magnetic.values().filter(|c| c.is_some()).count();
        if bindings.len() != expected
            || materials.len() != triangle_count
            || flux_density_t.len() != triangle_count
            || routes.magnetic.len() != triangle_count
            || routes.magnetic.keys().copied().ne(0..triangle_count)
        {
            return Err(invalid(
                "complete original field triangles and one derivative law per thermal triangle required",
            ));
        }
        let mut initial = BTreeMap::new();
        let mut capacity = BTreeMap::<CellId, f64>::new();
        let mut ranges = BTreeMap::<CellId, [f64; 2]>::new();
        for binding in bindings {
            check(cancel)?;
            let triangle = binding.triangle;
            let id = routes
                .magnetic
                .get(&triangle)
                .copied()
                .flatten()
                .ok_or_else(|| {
                    invalid("magnetic thermal law requires an original receiving cell")
                })?;
            if initial.contains_key(&triangle)
                || !binding.maximum_flux_density_t.is_finite()
                || binding.maximum_flux_density_t < 0.0
            {
                return Err(invalid(
                    "unique original magnetic triangles and finite nonnegative B envelopes required",
                ));
            }
            let source = network
                .cells
                .iter()
                .find(|s| s.cell == id)
                .ok_or_else(|| invalid("magnetic law lacks its caloric material binding"))?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or_else(|| invalid("magnetic caloric receiver is inactive"))?;
            let phase = source
                .curve
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|e| material(MaterialNetworkError::Phase(e)))?;
            let query = source
                .point
                .at(phase)
                .map_err(|e| material(MaterialNetworkError::Material(e)))?;
            let at = binding
                .law
                .at_temperature(phase.temperature_k())
                .map_err(field_error)?;
            if source.point.material.content_hash() != binding.law.card().content_hash()
                || source.point.temperature_axis != binding.law.temperature_axis()
                || query != at.material().query()
                || !same_material(at.material(), materials[triangle])
            {
                return Err(invalid(
                    "magnetic field, derivative law and actual caloric context must use the same original source bundle",
                ));
            }
            let volume = volumes_m3[triangle];
            if !volume.is_finite() || volume <= 0.0 {
                return Err(invalid("original magnetic volume must be finite positive"));
            }
            let b = magnitude(flux_density_t[triangle]);
            if !b.is_finite() || b > binding.maximum_flux_density_t {
                return Err(invalid(
                    "accepted magnetic field exceeds its thermal capacity envelope",
                ));
            }
            let range = binding.law.temperature_range_k();
            let knots = source.curve.knots();
            if range[0] < knots[0].temperature_k || range[1] > knots[knots.len() - 1].temperature_k
            {
                return Err(invalid(
                    "the complete magnetic temperature range requires caloric support",
                ));
            }
            let span = ranges.entry(id).or_insert(range);
            span[0] = span[0].max(range[0]);
            span[1] = span[1].min(range[1]);
            let bound = binding
                .law
                .minimum_heat_capacity_j_k(binding.maximum_flux_density_t, volume)
                .map_err(field_error)?;
            let total = capacity.entry(id).or_insert(0.0);
            *total = finite(*total + bound)?;
            initial.insert(
                triangle,
                Initial {
                    cell: id,
                    temperature_k: phase.temperature_k(),
                    b_t: b,
                    volume_m3: volume,
                    state: at.state(b, volume).map_err(field_error)?,
                },
            );
        }
        for (id, range) in ranges {
            check(cancel)?;
            if range[0] >= range[1] {
                return Err(invalid(
                    "a magnetic cell requires a common nonempty temperature interval",
                ));
            }
            let source = network
                .cells
                .iter()
                .find(|s| s.cell == id)
                .expect("admitted binding");
            let mut minimum = f64::INFINITY;
            for pair in source.curve.knots().windows(2) {
                check(cancel)?;
                let dt = pair[1].temperature_k - pair[0].temperature_k;
                if dt <= 0.0 {
                    return Err(invalid(
                        "magnetic caloric chart requires strictly increasing temperature",
                    ));
                }
                if pair[1].temperature_k >= range[0] && pair[0].temperature_k <= range[1] {
                    let value = finite(
                        thermal.cells()[&id].mass_kg
                            * (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
                            / dt,
                    )?;
                    minimum = minimum.min(value);
                }
            }
            if !minimum.is_finite() || finite(minimum + capacity[&id])? <= 0.0 {
                return Err(invalid(
                    "caloric plus ALL magnetic field heat capacities must stay positive over the declared envelope",
                ));
            }
        }
        // A fixed external region has no caloric receiver. Its changing field
        // may not carry an unowned temperature-dependent entropy contribution.
        let mut fixed = BTreeMap::new();
        for (&triangle, &cell) in &routes.magnetic {
            check(cancel)?;
            if cell.is_some() {
                continue;
            }
            let m = materials[triangle];
            let properties = [
                MAGNETIC_PERMEABILITY_PROPERTY,
                REVERSIBLE_FIELD_STRENGTH_PROPERTY,
            ]
            .map(|name| m.property(name));
            if properties.iter().all(Option::is_none) {
                return Err(invalid(
                    "fixed magnetic region lacks an original field source",
                ));
            }
            if properties.into_iter().flatten().any(|property| {
                property.source_axes().iter().any(|(name, _)| {
                    name == "T"
                        || m.axis_quantities().get(name)
                            == Some(&QuantitySpec::dimensional(Temperature::DIMS))
                        || m.axis_quantities()
                            .get(name)
                            .copied()
                            .and_then(QuantitySpec::semantic_type)
                            .is_some_and(|q| {
                                q.kind() == QuantityKind::AbsoluteTemperature
                                    && q.form() == ValueForm::Static
                            })
                })
            }) {
                return Err(invalid(
                    "temperature-dependent magnetic regions require a thermal receiver and derivative law",
                ));
            }
            fixed.insert(triangle, m.clone());
        }
        Ok(Self {
            laws: bindings.iter().map(|b| (b.triangle, b)).collect(),
            initial,
            fixed,
            triangle_count,
        })
    }

    /// Match the full refreshed source bundles before the native field solve.
    pub(in crate::circuit_enthalpy) fn validate_sources(
        &self,
        thermal: &DepositionSimulation<'_>,
        samples: &[MagneticThermalSample],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error> {
        if samples.len() != self.triangle_count {
            return Err(invalid(
                "magnetic caloric refresh omitted original triangles",
            ));
        }
        for (triangle, sample) in samples.iter().enumerate() {
            check(cancel)?;
            if sample.triangle != triangle {
                return Err(invalid("magnetic caloric refresh changed triangle order"));
            }
            if let Some(before) = self.initial.get(&triangle) {
                let at = self.laws[&triangle]
                    .law
                    .at_temperature(thermal.temperature_k(before.cell)?)
                    .map_err(field_error)?;
                if sample.cell != Some(before.cell)
                    || !same_material(&sample.material, at.material())
                {
                    return Err(invalid(
                        "native magnetic refresh and derivative law selected different source bundles or receivers",
                    ));
                }
            } else if sample.cell.is_some()
                || !same_material(&sample.material, &self.fixed[&triangle])
            {
                return Err(invalid("fixed external magnetic source changed"));
            }
        }
        Ok(())
    }

    /// Add only signed reversible heat; callers separately route native Joule
    /// heat. `sources_w` belongs to the private trial and is discarded on error.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::circuit_enthalpy) fn deposit<const D: usize>(
        &self,
        thermal: &DepositionSimulation<'_>,
        cells: impl IntoIterator<Item = MagneticEnergyCell<D>>,
        sources_w: &mut BTreeMap<CellId, f64>,
        dt: f64,
        energy_tolerance_j: f64,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MagneticCaloricExchange, Error> {
        if !dt.is_finite()
            || dt <= 0.0
            || !energy_tolerance_j.is_finite()
            || energy_tolerance_j < 0.0
        {
            return Err(invalid(
                "positive interval and finite magnetic energy allowance required",
            ));
        }
        let mut fields = Vec::with_capacity(self.initial.len());
        let mut internal = [0.0; 2];
        let mut reversible_heat_j = 0.0;
        let mut count = 0;
        for cell in cells {
            check(cancel)?;
            if cell.triangle != count || count >= self.triangle_count {
                return Err(invalid(
                    "native magnetic receipts changed original triangle order",
                ));
            }
            count += 1;
            finite(cell.parameter_work_j)?;
            for value in cell.free_energy_j {
                finite(value)?;
            }
            if let Some(before) = self.initial.get(&cell.triangle) {
                let binding = self.laws[&cell.triangle];
                let b = magnitude(cell.flux_density_t[1]);
                if !b.is_finite() || b > binding.maximum_flux_density_t {
                    return Err(invalid(
                        "proposed magnetic field exceeds its thermal capacity envelope",
                    ));
                }
                let at = binding
                    .law
                    .at_temperature(thermal.temperature_k(before.cell)?)
                    .map_err(field_error)?;
                let after = at.state(b, before.volume_m3).map_err(field_error)?;
                if finite(after.free_energy_j - cell.free_energy_j[1])?.abs() > energy_tolerance_j
                    || finite(before.state.free_energy_j - cell.free_energy_j[0])?.abs()
                        > energy_tolerance_j
                {
                    return Err(invalid(
                        "native local magnetic energy and sourced thermodynamic law disagree",
                    ));
                }
                let delta_ts = finite(
                    at.temperature_k() * after.entropy_j_k
                        - before.temperature_k * before.state.entropy_j_k,
                )?;
                let heat = finite(-cell.parameter_work_j - delta_ts)?;
                let power = finite(heat / dt)?;
                if heat != 0.0 && power == 0.0 {
                    return Err(invalid(
                        "magnetic reversible heat disappeared during power conversion",
                    ));
                }
                let destination = sources_w.entry(before.cell).or_insert(0.0);
                *destination = finite(*destination + power)?;
                reversible_heat_j = finite(reversible_heat_j + power * dt)?;
                internal[0] = finite(internal[0] + before.state.internal_energy_j)?;
                internal[1] = finite(internal[1] + after.internal_energy_j)?;
                fields.push(MagneticCaloricSample {
                    triangle: cell.triangle,
                    cell: before.cell,
                    temperature_k: [before.temperature_k, at.temperature_k()],
                    flux_density_t: [before.b_t, b],
                    state: [before.state, after],
                    parameter_work_j: cell.parameter_work_j,
                    reversible_heat_j: heat,
                });
            } else {
                if cell.parameter_work_j != 0.0 {
                    return Err(invalid(
                        "fixed region acquired unowned magnetic parameter work",
                    ));
                }
                for (target, value) in internal.iter_mut().zip(cell.free_energy_j) {
                    *target = finite(*target + value)?;
                }
            }
        }
        if count != self.triangle_count {
            return Err(invalid(
                "native magnetic receipts omitted original triangles",
            ));
        }
        Ok(MagneticCaloricExchange {
            fields,
            field_internal_energy_j: internal,
            reversible_heat_j,
        })
    }
}
