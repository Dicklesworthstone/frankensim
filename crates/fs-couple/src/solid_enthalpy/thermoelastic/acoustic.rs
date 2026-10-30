//! The accepted solid/heat history supplies local isentropic vibration data.
//!
//! This snapshot eliminates the local reversible temperature increment; it
//! does not evolve a second thermal history. Reference mass, original source
//! pins, actual total strain and the caloric chart's local derivative all come
//! from the existing coupled owner. The resulting materials feed the native
//! tetrahedral operator with incremental displacements about the accepted state.
//!
//! Do not reapply thermal eigenstrain to this perturbation problem or add this
//! adiabatic correction to a simulation already solving the same reversible
//! heat exchange. It is the no-heat-flow incremental limit, not a model of
//! thermoelastic attenuation, finite-frequency diffusion, or prestress geometry.

use super::*;
use fs_material::elastic::thermal::acoustic::{
    AdiabaticThermoelasticInput, AdiabaticThermoelasticState,
};
use fs_solid::linear3::TetThermalStressReport;

/// Frozen local entropy-constrained tangent on the original reference mesh.
pub struct AdiabaticSolidState<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    budget: TetAssemblyBudget,
    time_s: f64,
    mesh_frame: ContentHash,
    mechanical_response_identity: ContentHash,
    cells: Vec<CellId>,
    responses: Vec<AdiabaticThermoelasticState>,
    materials: Vec<TetElasticMaterial>,
}

impl AdiabaticSolidState<'_> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }
    pub const fn mechanical_response_identity(&self) -> ContentHash {
        self.mechanical_response_identity
    }
    /// Persistent thermal cell IDs in the original tetrahedron order.
    pub fn cells(&self) -> &[CellId] {
        &self.cells
    }
    /// Exact source, caloric and total-strain inputs to each local tangent.
    pub fn responses(&self) -> &[AdiabaticThermoelasticState] {
        &self.responses
    }
    /// Ordinary native materials with unchanged reference inertial densities.
    pub fn materials(&self) -> &[TetElasticMaterial] {
        &self.materials
    }

    /// Native operator for INCREMENTAL motion about this accepted state.
    ///
    /// The tangent contains local thermal feedback but no geometric initial-
    /// stress stiffness. Thermal eigenstrain belongs to the reference solution
    /// and must not be applied again to these incremental displacements.
    pub fn problem<'s>(&'s self, fixed_dofs: &'s [usize]) -> TetLinearElasticProblem<'s> {
        TetLinearElasticProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            materials: TetMaterialField::PerElement(&self.materials),
            fixed_dofs,
            budget: self.budget,
        }
    }
}

impl<'s, 'm> ThermoelasticSolid<'s, 'm> {
    /// Derive native adiabatic wave/modal data from the accepted equilibrium.
    /// No material constant, temperature, density or heat capacity is resupplied.
    pub fn adiabatic_state(
        &self,
        cx: &Cx<'_>,
    ) -> Result<AdiabaticSolidState<'m>, ThermoelasticFeedbackError> {
        resolve(
            self.source,
            &self.laws,
            &self.state,
            &self.equilibrium.solution.stress,
            cx,
        )
    }
}

pub(super) fn resolve<'m>(
    source: &SolidEnthalpy<'m>,
    laws: &[AffineIsotropicThermoelastic],
    state: &SolidEnthalpyState<'m>,
    stress: &TetThermalStressReport,
    cx: &Cx<'_>,
) -> Result<AdiabaticSolidState<'m>, ThermoelasticFeedbackError> {
    check(cx)?;
    if source.cells.len() != laws.len()
        || source.cells.len() != state.samples.len()
        || stress.elements().len() != laws.len()
        || stress.mesh_frame() != state.mesh_frame
    {
        return Err(ThermoelasticFeedbackError::Invalid(
            "adiabatic response requires the complete accepted material and stress fields",
        ));
    }
    let mut cells = Vec::with_capacity(laws.len());
    let mut responses = Vec::with_capacity(laws.len());
    let mut materials = Vec::with_capacity(laws.len());
    for (i, (bound, law)) in source.cells.iter().zip(laws).enumerate() {
        check(cx)?;
        let sample = &state.samples[i];
        let at = law.at_materials(&sample.elastic, sample.expansion.current())?;
        let range = law.temperature_range_k();
        if at.temperature_k() <= range[0] || at.temperature_k() >= range[1] {
            return Err(ThermoelasticFeedbackError::Invalid(
                "adiabatic waves require a two-sided thermoelastic temperature neighborhood",
            ));
        }
        let curve = bound.source.curve;
        let h = sample.phase.specific_enthalpy_j_kg();
        let phase = curve
            .state_at_specific_enthalpy(h)
            .map_err(SolidEnthalpyError::from)?;
        if phase != sample.phase
            || !phase.is_fully_solid()
            || phase.material_card_identity() != at.elastic.card_identity()
            || phase.temperature_k().to_bits() != at.temperature_k().to_bits()
        {
            return Err(ThermoelasticFeedbackError::Invalid(
                "adiabatic heat and elasticity must retain the same solid material and temperature",
            ));
        }
        // A fixed caloric chart cannot silently change its original pressure,
        // composition, processing condition or other external coordinates.
        let query = at.elastic.query();
        for receipt in curve.source_receipts() {
            let expected: Vec<_> = query
                .axes()
                .iter()
                .filter(|(axis, _)| axis.as_str() != "T")
                .collect();
            let actual: Vec<_> = receipt
                .query_point
                .iter()
                .filter(|(axis, _)| axis != "T")
                .collect();
            if expected.len() != actual.len()
                || expected.into_iter().zip(actual).any(|((a, x), (b, y))| {
                    a != b
                        || x.to_bits() != y.to_bits()
                        || query.axis_quantities().get(a) != receipt.axis_quantities.get(b)
                })
            {
                return Err(ThermoelasticFeedbackError::Invalid(
                    "adiabatic caloric chart must preserve the original fixed external source context",
                ));
            }
        }
        let density = bound.reference_cell.mass_kg / bound.reference_volume_m3;
        let capacity = caloric_capacity(curve, h, density)?;
        let free = state.strains[i].free_strain_mandel();
        let elastic = stress.elements()[i].elastic_strain_mandel();
        let total_strain_mandel = core::array::from_fn(|k| elastic[k] + free[k]);
        let response = at.adiabatic_acoustics(AdiabaticThermoelasticInput {
            total_strain_mandel,
            reference_density_kg_m3: density,
            caloric_heat_capacity_j_m3_k: capacity,
            mechanical_material_identity: ContentHash(state.materials[i].material_state_identity),
            caloric_state_identity: phase.identity(),
        })?;
        if response.isothermal_stiffness_mandel_pa() != *state.materials[i].stiffness_mandel_pa() {
            return Err(ThermoelasticFeedbackError::Invalid(
                "adiabatic thermal potential must use the actual native elastic stiffness",
            ));
        }
        let acoustic = response.acoustic();
        let material = TetElasticMaterial::try_new_mandel(
            acoustic.density_kg_m3(),
            *acoustic.stiffness_mandel_pa(),
            acoustic.material_state_identity().0,
        )
        .map_err(SolidEnthalpyError::from)?;
        cells.push(sample.cell);
        responses.push(response);
        materials.push(material);
    }
    check(cx)?;
    Ok(AdiabaticSolidState {
        nodes_m: state.nodes_m,
        tetrahedra: state.tetrahedra,
        budget: state.budget,
        time_s: state.time_s,
        mesh_frame: state.mesh_frame,
        mechanical_response_identity: stress.identity(),
        cells,
        responses,
        materials,
    })
}

fn caloric_capacity(
    curve: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    density: f64,
) -> Result<f64, ThermoelasticFeedbackError> {
    let knots = curve.knots();
    if h <= knots[0].specific_enthalpy_j_kg || h >= knots[knots.len() - 1].specific_enthalpy_j_kg {
        return Err(ThermoelasticFeedbackError::Invalid(
            "adiabatic waves require a two-sided caloric neighborhood",
        ));
    }
    let upper = knots.partition_point(|knot| knot.specific_enthalpy_j_kg <= h);
    let exact_knot = knots[upper - 1].specific_enthalpy_j_kg == h;
    let first = if exact_knot { upper - 2 } else { upper - 1 };
    let mut previous = None;
    for pair in knots[first..=upper].windows(2) {
        if pair.iter().any(|knot| {
            knot.liquid_mass_fraction != 0.0 || !same_measure(knot.bulk_density_kg_m3, density)
        }) || pair[1].temperature_k <= pair[0].temperature_k
        {
            return Err(ThermoelasticFeedbackError::Invalid(
                "adiabatic solid response requires fully solid constant-reference-density caloric support",
            ));
        }
        let capacity = (pair[1].specific_enthalpy_j_kg - pair[0].specific_enthalpy_j_kg)
            / (pair[1].temperature_k - pair[0].temperature_k);
        if !capacity.is_finite()
            || capacity <= 0.0
            || previous.is_some_and(|left| !same_measure(left, capacity))
        {
            return Err(ThermoelasticFeedbackError::Invalid(
                "a slope-changing caloric knot has no unique adiabatic linearization",
            ));
        }
        previous = Some(capacity);
    }
    let dt_dh = curve
        .temperature_derivative_at_specific_enthalpy(h)
        .map_err(SolidEnthalpyError::from)?;
    let capacity = finite(density / dt_dh)?;
    if capacity <= 0.0 {
        return Err(ThermoelasticFeedbackError::Invalid(
            "caloric capacity per reference volume must be positive",
        ));
    }
    Ok(capacity)
}
