//! Retained finite thermal, plastic and creep histories on native P1/Q1 sheets.
//!
//! The original temperature-dependent material family is retained per element,
//! with one temperature and complete inelastic history per quadrature point.
//! Each physical proposal freezes the actual endpoint temperature and duration
//! once, then uses the shared sheet assembly and Newton solver. Every local
//! traction-free return starts from the same accepted history. Old states are
//! inspected at their old temperatures; softening may cause flow during a hold.
//!
//! Native free/internal energy, physical plastic/creep loss, reversible caloric
//! exchange and signed endpoint-work defect remain distinct. Only physical
//! loss plus reversible exchange belongs in an external caloric solve. The
//! incremental potential is the native variational potential, never replaced
//! by full viscous heat. Reference area times original thickness weights all
//! accounts. This owner does not infer caloric storage, conductivity or density.
//!
//! This is an in-plane membrane with solved transverse stretch and shear, not
//! a bending shell or a model of wrinkling, necking, fracture or contact.

mod transfer;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::FiniteJ2State;
use fs_material::plastic::finite::plane_stress::FinitePlaneStressResponse;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteThermalJ2Material, FrozenThermalFiniteJ2,
};

use super::{
    Element, FinitePlasticSheetError, FinitePlasticSheetEvaluation, SheetAssemblyBudget,
    SheetOperator, checkpoint, endpoint_forces, energy_gate, invalid, prepare_geometry,
};
use crate::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use crate::mesh2::Mesh2;
use crate::plastic3::increment::{self, IncrementControls};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

pub use super::FinitePlasticSheetEvaluation as FiniteThermoplasticSheetEvaluation;
/// The same original-volume thermodynamic account used by finite solid points.
pub use crate::hyper3::plastic::thermal::FiniteThermoplasticElementEnergy;
/// Point-density account; this is the original type, not a second stored report.
pub type FiniteThermoplasticSheetPointEnergy = FiniteThermoplasticElementEnergy;
/// Shared whole-solid accounts; no inertial terms are added by this sheet owner.
pub use crate::hyper3::plastic::thermal::FiniteThermoplasticTetEnergy as FiniteThermoplasticSheetEnergy;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Original geometry and source families with an immutable accepted endpoint.
#[derive(Debug)]
pub struct FiniteThermoplasticSheet<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    mesh: &'a Mesh2,
    families: FinitePlasticLawField<'a, L>,
    thickness_m: &'a [f64],
    elements: Arc<Vec<Element>>,
    local: PlaneStressSettings,
    budget: SheetAssemblyBudget,
    displacement_m: Vec<[f64; 2]>,
    states: Vec<FiniteJ2State>,
    deformation: Vec<[f64; 9]>,
    piola_pa: Vec<[f64; 9]>,
    temperatures_k: Vec<f64>,
    quadrature_temperatures_k: Vec<f64>,
    free_energy_j: f64,
    internal_energy_j: f64,
    physical_time_s: f64,
}

/// Complete accepted proposal; the previous sheet remains unchanged.
#[derive(Debug)]
pub struct FiniteThermoplasticSheetStep<'a, L: FiniteThermalJ2Material = AffineThermalFiniteJ2> {
    pub sheet: FiniteThermoplasticSheet<'a, L>,
    pub equilibrium: IncrementRecord,
    /// Native accounts per original reference volume, in element/point order.
    pub quadrature_energy: Vec<FiniteThermoplasticElementEnergy>,
    pub energy: FiniteThermoplasticSheetEnergy,
}

/// Preserve material-point, shared-solver and physical-account diagnoses.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteThermoplasticSheetError {
    Solid(FinitePlasticSheetError),
    Material {
        element: usize,
        point: usize,
        source: FiniteJ2ThermalError,
    },
}

impl From<SmallStrainTetError> for FiniteThermoplasticSheetError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error.into())
    }
}
impl From<FinitePlasticSheetError> for FiniteThermoplasticSheetError {
    fn from(error: FinitePlasticSheetError) -> Self {
        Self::Solid(error)
    }
}
impl core::fmt::Display for FiniteThermoplasticSheetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Material {
                element,
                point,
                source,
            } => {
                write!(
                    f,
                    "finite thermal sheet element {element} point {point}: {source}"
                )
            }
        }
    }
}
impl std::error::Error for FiniteThermoplasticSheetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Material { source, .. } => Some(source),
        }
    }
}

impl<'a, L: FiniteThermalJ2Material> FiniteThermoplasticSheet<'a, L> {
    /// Admit virgin geometry at one original reference temperature per element.
    /// Original source families remain borrowed, including a uniform binding.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        mesh: &'a Mesh2,
        families: FinitePlasticLawField<'a, L>,
        thickness_m: &'a [f64],
        initial_temperatures_k: &[f64],
        local: PlaneStressSettings,
        budget: SheetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticSheetError> {
        Self::build(
            mesh,
            families,
            thickness_m,
            initial_temperatures_k,
            false,
            local,
            budget,
            cx,
        )
    }

    /// Admit one original family per element and one initial temperature per
    /// native quadrature point (one for P1, four for Q1). Each temperature must
    /// equal its family's stress-free reference. Initializing a loaded or heated
    /// sheet requires a real accepted path, not an invented creep interval.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new_with_quadrature_temperatures(
        mesh: &'a Mesh2,
        families: FinitePlasticLawField<'a, L>,
        thickness_m: &'a [f64],
        initial_temperatures_k: &[f64],
        local: PlaneStressSettings,
        budget: SheetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticSheetError> {
        Self::build(
            mesh,
            families,
            thickness_m,
            initial_temperatures_k,
            true,
            local,
            budget,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        mesh: &'a Mesh2,
        families: FinitePlasticLawField<'a, L>,
        thickness_m: &'a [f64],
        initial_temperatures_k: &[f64],
        per_quadrature: bool,
        local: PlaneStressSettings,
        budget: SheetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        if let FinitePlasticLawField::PerElement(values) = families
            && values.len() != mesh.elems.len()
        {
            return Err(SmallStrainTetError::LawCountMismatch {
                expected: mesh.elems.len(),
                got: values.len(),
            }
            .into());
        }
        if !local.relative_tolerance.is_finite()
            || local.relative_tolerance <= 0.0
            || local.relative_tolerance >= 1.0
            || local.max_iterations == 0
            || local.max_backtracks == 0
        {
            return Err(invalid("thermal sheet needs a positive finite local tolerance below one and nonzero local iteration/backtracking budgets").into());
        }
        let elements = prepare_geometry(mesh, thickness_m, budget, cx)?;
        let count = elements.iter().map(|element| element.points.len()).sum();
        let expected = if per_quadrature {
            count
        } else {
            elements.len()
        };
        if initial_temperatures_k.len() != expected {
            return Err(invalid(
                "initial temperatures must match the selected element or quadrature-point layout",
            )
            .into());
        }
        let mut states = Vec::with_capacity(count);
        let mut piola_pa = Vec::with_capacity(count);
        let mut temperatures_k = Vec::with_capacity(elements.len());
        let mut quadrature_temperatures_k = Vec::with_capacity(count);
        let mut free_energy_j = 0.0;
        let mut internal_energy_j = 0.0;
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let family = match families {
                FinitePlasticLawField::Uniform(family) => family,
                FinitePlasticLawField::PerElement(families) => &families[e],
            }
            .thermal_finite_j2();
            let temperature = family.parameters().reference_temperature_k;
            temperatures_k.push(temperature);
            let map = |source| FiniteThermoplasticSheetError::Material {
                element: e,
                point: 0,
                source,
            };
            let law = family
                .at_temperature_for_inspection(temperature)
                .map_err(map)?;
            let initial = family.initial_state();
            law.validate_committed(&IDENTITY, &initial, 0.0)
                .map_err(map)?;
            let value = law.thermodynamics(&IDENTITY, &initial).map_err(map)?;
            for (q, point) in element.points.iter().enumerate() {
                let source_index = if per_quadrature {
                    element.offset + q
                } else {
                    e
                };
                if initial_temperatures_k[source_index] != temperature {
                    return Err(invalid("virgin finite thermal sheet temperatures must equal each original family's reference temperature").into());
                }
                states.push(initial.clone());
                quadrature_temperatures_k.push(temperature);
                piola_pa.push(value.first_piola_pa);
                free_energy_j = point
                    .volume_m3
                    .mul_add(value.free_energy_j_m3, free_energy_j);
                internal_energy_j = point
                    .volume_m3
                    .mul_add(value.internal_energy_j_m3, internal_energy_j);
            }
        }
        if !free_energy_j.is_finite() || !internal_energy_j.is_finite() {
            return Err(invalid("thermal sheet initial storage is not representable").into());
        }
        checkpoint(cx)?;
        Ok(Self {
            mesh,
            families,
            thickness_m,
            elements,
            local,
            budget,
            displacement_m: vec![[0.0; 2]; mesh.nodes.len()],
            states,
            deformation: vec![IDENTITY; count],
            piola_pa,
            temperatures_k,
            quadrature_temperatures_k,
            free_energy_j,
            internal_energy_j,
            physical_time_s: 0.0,
        })
    }

    pub fn families(&self) -> FinitePlasticLawField<'a, L> {
        self.families
    }
    /// Element temperatures for the original element-temperature API. After a
    /// nonuniform quadrature-temperature step these are ORIGINAL-volume weighted
    /// diagnostic means only; mechanics and history always use the retained
    /// [`Self::quadrature_temperatures_k`], never these means.
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    /// Actual accepted constitutive temperatures in element/point order.
    pub fn quadrature_temperatures_k(&self) -> &[f64] {
        &self.quadrature_temperatures_k
    }
    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.displacement_m
    }
    pub fn quadrature_states(&self) -> &[FiniteJ2State] {
        &self.states
    }
    pub fn quadrature_deformation(&self) -> &[[f64; 9]] {
        &self.deformation
    }
    pub fn quadrature_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.piola_pa
    }
    pub fn free_energy_j(&self) -> f64 {
        self.free_energy_j
    }
    /// Mechanical internal energy only; the caloric owner retains its own store.
    pub fn internal_energy_j(&self) -> f64 {
        self.internal_energy_j
    }
    pub fn physical_time_s(&self) -> f64 {
        self.physical_time_s
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.elements
            .get(element)
            .map(|element| element.points.iter().map(|point| point.volume_m3).sum())
    }
    pub fn quadrature_volume_m3(&self, element: usize, point: usize) -> Option<f64> {
        self.elements
            .get(element)?
            .points
            .get(point)
            .map(|point| point.volume_m3)
    }
    pub fn quadrature_thickness_m(&self) -> Vec<f64> {
        self.elements
            .iter()
            .enumerate()
            .flat_map(|(e, element)| {
                (0..element.points.len())
                    .map(move |q| self.thickness_m[e] * self.deformation[element.offset + q][8])
            })
            .collect()
    }

    fn family(&self, element: usize) -> &AffineThermalFiniteJ2 {
        match self.families {
            FinitePlasticLawField::Uniform(family) => family,
            FinitePlasticLawField::PerElement(families) => &families[element],
        }
        .thermal_finite_j2()
    }

    fn expand_element_temperatures(
        &self,
        temperatures_k: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        if temperatures_k.len() != self.elements.len() {
            return Err(invalid("one temperature per original sheet element required").into());
        }
        let mut points = Vec::with_capacity(self.states.len());
        for (element, &temperature) in self.elements.iter().zip(temperatures_k) {
            checkpoint(cx)?;
            points.extend(core::iter::repeat_n(temperature, element.points.len()));
        }
        Ok(points)
    }

    // Preserve the published element-sized diagnostic getter without using a
    // mean to drive mechanics. Constant point temperatures remain exact.
    fn element_temperature_means(
        &self,
        temperatures_k: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticSheetError> {
        let mut means = Vec::with_capacity(self.elements.len());
        for element in self.elements.iter() {
            checkpoint(cx)?;
            let values = &temperatures_k[element.offset..element.offset + element.points.len()];
            let first = values[0];
            let mean = if values.iter().all(|&temperature| temperature == first) {
                first
            } else {
                let minimum = values.iter().copied().fold(f64::INFINITY, f64::min);
                let volume: f64 = element.points.iter().map(|point| point.volume_m3).sum();
                element
                    .points
                    .iter()
                    .zip(values)
                    .fold(minimum, |sum, (point, &temperature)| {
                        (point.volume_m3 / volume).mul_add(temperature - minimum, sum)
                    })
            };
            if !mean.is_finite() || mean <= 0.0 {
                return Err(invalid(
                    "sheet diagnostic element temperature is not positive and representable",
                )
                .into());
            }
            means.push(mean);
        }
        Ok(means)
    }

    fn freeze(
        &self,
        temperatures_k: &[f64],
        duration_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<FrozenThermalFiniteJ2>, FiniteThermoplasticSheetError> {
        if temperatures_k.len() != self.states.len() || !duration_s.is_finite() || duration_s <= 0.0
        {
            return Err(invalid("thermal sheet requires one temperature per quadrature point and a positive finite physical duration").into());
        }
        let mut laws: Vec<FrozenThermalFiniteJ2> = Vec::with_capacity(self.states.len());
        for (e, element) in self.elements.iter().enumerate() {
            for q in 0..element.points.len() {
                checkpoint(cx)?;
                let index = element.offset + q;
                if q > 0 && temperatures_k[index] == temperatures_k[element.offset] {
                    laws.push(laws[element.offset].clone());
                } else {
                    laws.push(
                        self.family(e)
                            .at_temperature_for_step(temperatures_k[index], duration_s)
                            .map_err(|source| FiniteThermoplasticSheetError::Material {
                                element: e,
                                point: q,
                                source,
                            })?,
                    );
                }
            }
        }
        Ok(laws)
    }

    fn response(
        &self,
        laws: &[FrozenThermalFiniteJ2],
        e: usize,
        q: usize,
        u: &[f64],
        states: &[FiniteJ2State],
    ) -> Result<FinitePlaneStressResponse, FiniteJ2ThermalError> {
        let element = &self.elements[e];
        let index = element.offset + q;
        laws[index].plane_stress_response(
            &element.in_plane(&element.points[q], u),
            &states[index],
            [2, 5, 8].map(|i| self.deformation[index][i]),
            self.local,
        )
    }

    /// Equilibrium displacement of the accepted histories at intermediate
    /// point temperatures over the step's own duration. It is only a Newton
    /// start for the endpoint solve; nothing here is committed or published.
    fn stage_equilibrium(
        &self,
        temperatures_k: &[f64],
        start: &[f64],
        controls: IncrementControls<'_>,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticSheetError> {
        let laws = self.freeze(temperatures_k, dt_s, cx)?;
        let operator = SheetOperator {
            mesh: self.mesh,
            elements: &self.elements,
            states: &self.states,
            response: |e, q, u: &[f64], states: &[FiniteJ2State]| {
                self.response(&laws, e, q, u, states).map_err(|error| {
                    invalid(format!("thermal sheet element {e} point {q}: {error}"))
                })
            },
        };
        Ok(
            increment::solve(&operator, controls, Some((start, &self.states)), None, cx)?
                .displacement_m,
        )
    }

    /// Frozen-history evaluation at one endpoint temperature per element and
    /// one actual duration. Accepted temperature, history and time are unchanged.
    pub fn evaluate_trial(
        &self,
        temperatures_k: &[f64],
        displacement_m: &[[f64; 2]],
        dt_s: f64,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetEvaluation, FiniteThermoplasticSheetError> {
        let points = self.expand_element_temperatures(temperatures_k, cx)?;
        self.evaluate_trial_for_step(displacement_m, &points, dt_s, need_tangent, cx)
    }

    /// Evaluate a proposal at explicit quadrature-point temperatures and duration from
    /// the same accepted histories, without publishing a physical step.
    pub fn evaluate_trial_for_step(
        &self,
        displacement_m: &[[f64; 2]],
        temperatures_k: &[f64],
        duration_s: f64,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticSheetEvaluation, FiniteThermoplasticSheetError> {
        let laws = self.freeze(temperatures_k, duration_s, cx)?;
        let operator = SheetOperator {
            mesh: self.mesh,
            elements: &self.elements,
            states: &self.states,
            response: |e, q, u: &[f64], states: &[FiniteJ2State]| {
                self.response(&laws, e, q, u, states).map_err(|error| {
                    invalid(format!("thermal sheet element {e} point {q}: {error}"))
                })
            },
        };
        let result = operator
            .assemble(
                &self.states,
                displacement_m.as_flattened(),
                None,
                need_tangent,
                cx,
            )?
            .0;
        Ok(FinitePlasticSheetEvaluation {
            internal_force_n: result.force.chunks_exact(2).map(|u| [u[0], u[1]]).collect(),
            tangent_n_m: result.tangent,
            incremental_potential_j: result.potential,
        })
    }

    /// Solve one actual interval at one new temperature per element. Old point
    /// temperatures and histories remain canonical; this only expands the NEW
    /// supplied element temperatures to their original material points.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetStep<'a, L>, FiniteThermoplasticSheetError> {
        let points = self.expand_element_temperatures(new_temperatures_k, cx)?;
        self.advance_with_quadrature_temperatures(
            &points,
            prescribed_m,
            nodal_forces_n,
            dt_s,
            settings,
            controls,
            cx,
        )
    }

    /// Solve one physical endpoint using the original accepted point temperatures,
    /// deformation and complete histories. Temperatures are prescribed here;
    /// the returned signed caloric transfer can feed a separate thermal owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_quadrature_temperatures(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        dt_s: f64,
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermoplasticSheetStep<'a, L>, FiniteThermoplasticSheetError> {
        checkpoint(cx)?;
        let external = endpoint_forces(
            self.mesh,
            self.budget,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
        )?;
        let laws = self.freeze(new_temperatures_k, dt_s, cx)?;
        let physical_time_s = self.physical_time_s + dt_s;
        if !physical_time_s.is_finite() || physical_time_s <= self.physical_time_s {
            return Err(invalid("thermal sheet physical time cannot advance representably").into());
        }
        for (e, element) in self.elements.iter().enumerate() {
            for q in 0..element.points.len() {
                checkpoint(cx)?;
                let index = element.offset + q;
                self.family(e)
                    .validate_committed(
                        &self.deformation[index],
                        &self.states[index],
                        self.quadrature_temperatures_k[index],
                        controls.yield_tolerance_pa,
                    )
                    .map_err(|source| FiniteThermoplasticSheetError::Material {
                        element: e,
                        point: q,
                        source,
                    })?;
            }
        }
        let operator = SheetOperator {
            mesh: self.mesh,
            elements: &self.elements,
            states: &self.states,
            response: |e, q, u: &[f64], states: &[FiniteJ2State]| {
                self.response(&laws, e, q, u, states).map_err(|error| {
                    invalid(format!("thermal sheet element {e} point {q}: {error}"))
                })
            },
        };
        let controls_from = || IncrementControls {
            dofs: 2 * self.mesh.nodes.len(),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            settings,
        };
        let accepted_u = self.displacement_m.as_flattened();
        let solved = match increment::solve(
            &operator,
            controls_from(),
            Some((accepted_u, &self.states)),
            None,
            cx,
        ) {
            Ok(solved) => solved,
            // A large temperature change can leave no admissible frozen-history
            // plane-stress return at the accepted in-plane geometry (the
            // transverse strain needs more plastic flow than the card admits),
            // so the solve is refused at its start, before any Newton step.
            // Only then: reach intermediate point temperatures through
            // equilibria of the same accepted histories, and start this same
            // endpoint solve from the first admitted one. A step whose start
            // is admitted keeps its original path and bits.
            Err(error)
                if error != SmallStrainTetError::Cancelled
                    && operator
                        .assemble(&self.states, accepted_u, None, false, cx)
                        .is_err() =>
            {
                let mut start = accepted_u.to_vec();
                let mut reached = 0.0_f64;
                let mut fraction = 0.5_f64;
                let mut stages = 0;
                loop {
                    checkpoint(cx)?;
                    if reached > 0.0
                        && operator
                            .assemble(&self.states, &start, None, false, cx)
                            .is_ok()
                    {
                        break increment::solve(
                            &operator,
                            controls_from(),
                            Some((&start, &self.states)),
                            None,
                            cx,
                        )?;
                    }
                    if stages == 16 || fraction < 1.0 / 64.0 {
                        return Err(error.into());
                    }
                    stages += 1;
                    let s = fraction.mul_add(1.0 - reached, reached);
                    let temperatures: Vec<f64> = self
                        .quadrature_temperatures_k
                        .iter()
                        .zip(new_temperatures_k)
                        .map(|(&old, &new)| (new - old).mul_add(s, old))
                        .collect();
                    match self.stage_equilibrium(&temperatures, &start, controls_from(), dt_s, cx) {
                        Ok(displacement) => {
                            start = displacement;
                            reached = s;
                            fraction = 0.5;
                        }
                        Err(_) => fraction *= 0.5,
                    }
                }
            }
            Err(error) => return Err(error.into()),
        };
        let equilibrium = solved
            .increments
            .into_iter()
            .next()
            .expect("one physical endpoint");
        // The solver has committed its private candidate. Inspect the assembly
        // again from ORIGINAL histories so physical time is never applied twice.
        let assembled = operator
            .assemble(&self.states, &solved.displacement_m, None, false, cx)?
            .0;
        let mut energy = FiniteThermoplasticSheetEnergy::default();
        let mut quadrature_energy = Vec::with_capacity(self.states.len());
        let mut deformation = Vec::with_capacity(self.states.len());
        let mut piola_pa = Vec::with_capacity(self.states.len());
        let (mut old_free, mut new_free, mut old_internal, mut new_internal) = (0.0, 0.0, 0.0, 0.0);
        let tolerance = controls.equilibrium_work_tolerance_j;
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let map = |source| FiniteThermoplasticSheetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                let response = self
                    .response(&laws, e, q, &solved.displacement_m, &self.states)
                    .map_err(map)?;
                let value = self
                    .family(e)
                    .checked_increment_for_step(
                        &self.deformation[index],
                        &self.states[index],
                        self.quadrature_temperatures_k[index],
                        &response.deformation_gradient,
                        new_temperatures_k[index],
                        controls.yield_tolerance_pa,
                        dt_s,
                    )
                    .map_err(map)?;
                if value.mechanical.state != solved.point_states[index]
                    || value.mechanical != response.material
                {
                    return Err(invalid("thermal sheet endpoint audit differs from its accepted frozen-history return").into());
                }
                let account = FiniteThermoplasticSheetPointEnergy::from_increment(&value);
                let v = point.volume_m3;
                if [
                    account.plastic_dissipation_j_m3,
                    account.creep_dissipation_j_m3,
                ]
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
                {
                    return Err(invalid(format!(
                        "thermal sheet element {e} point {q} has invalid physical dissipation"
                    ))
                    .into());
                }
                energy_gate(
                    "point thermal identity",
                    v * account.balance_error_j_m3,
                    tolerance,
                )?;
                energy_gate(
                    "point free energy",
                    v * (account.free_after_j_m3 - value.mechanical.stored_energy_j_m3),
                    tolerance,
                )?;
                old_free = v.mul_add(account.free_before_j_m3, old_free);
                new_free = v.mul_add(account.free_after_j_m3, new_free);
                old_internal = v.mul_add(account.internal_before_j_m3, old_internal);
                new_internal = v.mul_add(account.internal_after_j_m3, new_internal);
                energy.free_energy_change_j = v.mul_add(
                    account.free_after_j_m3 - account.free_before_j_m3,
                    energy.free_energy_change_j,
                );
                energy.internal_energy_change_j = v.mul_add(
                    account.internal_energy_change_j_m3(),
                    energy.internal_energy_change_j,
                );
                energy.plastic_dissipation_j = v.mul_add(
                    account.plastic_dissipation_j_m3,
                    energy.plastic_dissipation_j,
                );
                energy.creep_dissipation_j =
                    v.mul_add(account.creep_dissipation_j_m3, energy.creep_dissipation_j);
                energy.reversible_heat_j =
                    v.mul_add(account.reversible_heat_j_m3, energy.reversible_heat_j);
                energy.caloric_heat_j =
                    v.mul_add(account.caloric_heat_j_m3(), energy.caloric_heat_j);
                energy.coefficient_storage_change_j = v.mul_add(
                    account.coefficient_storage_change_j_m3,
                    energy.coefficient_storage_change_j,
                );
                energy.thermoelastic_storage_change_j = v.mul_add(
                    account.thermoelastic_storage_change_j_m3,
                    energy.thermoelastic_storage_change_j,
                );
                energy.endpoint_work_defect_j = v.mul_add(
                    account.endpoint_work_defect_j_m3,
                    energy.endpoint_work_defect_j,
                );
                energy.incremental_potential_j = v.mul_add(
                    account.incremental_potential_j_m3,
                    energy.incremental_potential_j,
                );
                energy.internal_endpoint_work_j = v.mul_add(
                    account.internal_endpoint_work_j_m3,
                    energy.internal_endpoint_work_j,
                );
                deformation.push(response.deformation_gradient);
                piola_pa.push(value.mechanical.first_piola_pa);
                quadrature_energy.push(account);
            }
        }
        let old_u = self.displacement_m.as_flattened();
        let delta_u: Vec<_> = solved
            .displacement_m
            .iter()
            .zip(old_u)
            .map(|(new, old)| new - old)
            .collect();
        let assembled_work: f64 = assembled
            .force
            .iter()
            .zip(&delta_u)
            .map(|(force, change)| force * change)
            .sum();
        energy.dead_load_endpoint_work_j = external
            .iter()
            .zip(&delta_u)
            .map(|(force, change)| force * change)
            .sum();
        energy.support_endpoint_work_j = prescribed_m
            .iter()
            .zip(&equilibrium.reactions_n)
            .map(|(&(dof, _), reaction)| reaction * delta_u[dof])
            .sum();
        let external_work = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j;
        energy.equilibrium_work_residual_j = external_work - energy.internal_endpoint_work_j;
        energy.balance_error_j = external_work
            - energy.internal_energy_change_j
            - energy.caloric_heat_j
            - energy.endpoint_work_defect_j;
        for (quantity, error) in [
            ("accepted free energy", old_free - self.free_energy_j),
            (
                "accepted internal energy",
                old_internal - self.internal_energy_j,
            ),
            (
                "assembled native potential",
                assembled.potential - energy.incremental_potential_j,
            ),
            (
                "free-energy increment",
                new_free - old_free - energy.free_energy_change_j,
            ),
            (
                "internal-energy increment",
                new_internal - old_internal - energy.internal_energy_change_j,
            ),
            (
                "traction-free point/assembled work",
                assembled_work - energy.internal_endpoint_work_j,
            ),
            (
                "caloric decomposition",
                energy.caloric_heat_j - energy.inelastic_dissipation_j() - energy.reversible_heat_j,
            ),
            (
                "support/dead-load endpoint work",
                energy.equilibrium_work_residual_j,
            ),
            ("total thermal endpoint identity", energy.balance_error_j),
        ] {
            energy_gate(quantity, error, tolerance)?;
        }
        energy_gate(
            "signed mechanical endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if !new_free.is_finite()
            || new_free < 0.0
            || !new_internal.is_finite()
            || !energy.plastic_dissipation_j.is_finite()
            || energy.plastic_dissipation_j < 0.0
            || !energy.creep_dissipation_j.is_finite()
            || energy.creep_dissipation_j < 0.0
            || [
                energy.free_energy_change_j,
                energy.internal_energy_change_j,
                energy.coefficient_storage_change_j,
                energy.thermoelastic_storage_change_j,
                energy.reversible_heat_j,
                energy.caloric_heat_j,
                energy.incremental_potential_j,
                energy.internal_endpoint_work_j,
                energy.dead_load_endpoint_work_j,
                energy.support_endpoint_work_j,
            ]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid("finite thermal sheet storage/accounts and nonnegative physical plastic/creep loss required").into());
        }
        let temperatures_k = self.element_temperature_means(new_temperatures_k, cx)?;
        checkpoint(cx)?;
        Ok(FiniteThermoplasticSheetStep {
            sheet: Self {
                mesh: self.mesh,
                families: self.families,
                thickness_m: self.thickness_m,
                elements: Arc::clone(&self.elements),
                local: self.local,
                budget: self.budget,
                displacement_m: solved
                    .displacement_m
                    .chunks_exact(2)
                    .map(|u| [u[0], u[1]])
                    .collect(),
                states: solved.point_states,
                deformation,
                piola_pa,
                temperatures_k,
                quadrature_temperatures_k: new_temperatures_k.to_vec(),
                free_energy_j: new_free,
                internal_energy_j: new_internal,
                physical_time_s,
            },
            equilibrium,
            quadrature_energy,
            energy,
        })
    }
}
