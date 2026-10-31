//! Source-bound tensor shear projected into the original finite layer geometry.
//!
//! The physical columns are epsilon_i = sym(t_i tensor n)/h and volume A*h.
//! Every original 3-D equilibrium/Maxwell cross term is projected in its declared
//! orientation. This is constrained shear kinematics, not plane-stress/normal
//! contact elimination. Only the already-existing shear wrenches do work here.
//! Shared source thermal laws contribute their complete mechanical entropy and
//! internal energy at each actual endpoint temperature. Temperature may also change relaxation times, Kelvin
//! viscosity and plastic kinetics. The original joint rigid/enthalpy residual
//! evaluates each law at its endpoint trial temperature and accepts all
//! bodies/history/cells once; dissipative heat excludes reversible storage work.

use super::thermal::RigidLayerThermalTransport;
use super::*;
use fs_matdb::{InterpolationPolicy, MaterialCard, PropertyValue};
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, PLASTIC_SHEAR_VISCOSITY_PROPERTY,
    ResolvedMaterialStatePoint, SHEAR_LAYER_YIELD_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::visco_tensor::{
    ProjectedTensorProny, TensorProjectionBudget,
    source::{ResolvedTensorMaterial, TensorSourceError},
    source::thermal::{TensorTemperatureDerivative, TensorThermalLaw},
};
use fs_qty::{DynViscosity, Pressure, Temperature};
use fs_thermal::enthalpy::ThermalLaw;
use fs_tribo::layer_contact::tensor::{ShearStiffness, ShearViscosity, TensorLayerRelaxationArm};
use std::convert::Infallible;

#[derive(Debug)]
pub enum TensorLayerSourceError<E = Infallible> {
    Invalid(&'static str),
    Tensor(TensorSourceError),
    Material(MaterialStatePointError),
    Constitutive(TriboError),
    Transport(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for TensorLayerSourceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source tensor shear layer: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for TensorLayerSourceError<E> {}

/// Complete source receipts at one actual layer temperature, in original order.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorLayerSourceState {
    pub material: ResolvedTensorMaterial,
    pub plastic: ResolvedMaterialStatePoint,
}

/// Immutable material and projection; the runtime history remains with the
/// existing rigid layer owner. There are no authored spring or damping values.
pub struct TensorLayerMaterial<'a> {
    source: TensorLayerSourceState,
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    temperature_axis: String,
    material_to_world: [[f64; 3]; 3],
    geometry: RigidLayerGeometry,
    maximum_angular_frequency_rad_s: f64,
    projection_budget: TensorProjectionBudget,
    thermal_slopes: TensorLayerThermalSlopes,
    thermal_derivative: TensorTemperatureDerivative,
    thermal_projection: ProjectedTensorProny,
}
impl<'a> TensorLayerMaterial<'a> {
    /// Source-bound fully solid, constant-density tensor potential, including
    /// independent affine E/G with exact rational bulk thermodynamics. Every
    /// Maxwell coefficient may have its own affine slope. Retained physical
    /// spring strains preserve even rank-deficient branch history.
    ///
    /// Positive total heat capacity is established over the original full 3-D
    /// strain domain, not just this shear subspace. This sufficient bound can
    /// refuse a safe pure-shear reduction of a strongly curved bulk law.
    /// Changing geometry and thermal expansion remain outside this fixed law.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        material: ResolvedTensorMaterial,
        card: &'a MaterialCard,
        caloric: &'a ResolvedHeatCapacityCurve,
        temperature_axis: &str,
        material_to_world: [[f64; 3]; 3],
        geometry: RigidLayerGeometry,
        plastic_selection: MaterialPropertySelection,
        maximum_angular_frequency_rad_s: f64,
        projection_budget: TensorProjectionBudget,
    ) -> Result<Self, TensorLayerSourceError> {
        geometry
            .validate()
            .map_err(|_| TensorLayerSourceError::Invalid("invalid physical shear geometry"))?;
        let knots = caloric.curve().knots();
        if knots.len() < 2
            || caloric.curve().material_card_identity() != card.content_hash()
            || material.state().card_identity() != card.content_hash()
            || knots.iter().any(|k| {
                k.liquid_mass_fraction != 0.0 || k.bulk_density_kg_m3 != material.density_kg_m3()
            })
            || knots
                .windows(2)
                .any(|p| p[0].temperature_k >= p[1].temperature_k)
            || !maximum_angular_frequency_rad_s.is_finite()
            || maximum_angular_frequency_rad_s <= 0.0
        {
            return Err(TensorLayerSourceError::Invalid(
                "one fully solid source and reference density required",
            ));
        }
        let interval = [knots[0].temperature_k, knots[knots.len() - 1].temperature_k];
        let thermal = TensorThermalLaw::try_new(&material, card, temperature_axis, interval)
            .map_err(TensorLayerSourceError::Tensor)?;
        let derivative = thermal
            .in_frame(material_to_world)
            .map_err(TensorLayerSourceError::Tensor)?;
        let projected = project_layer(&material, material_to_world, geometry, projection_budget)?;
        let temperature = material.state().state_coordinate(temperature_axis).ok_or(
            TensorLayerSourceError::Invalid("missing actual layer temperature"),
        )?;
        let capacity_bound = projected
            .thermal_capacity_bound_j_k2(&derivative, temperature, [0.0; 6])
            .map_err(|_| {
                TensorLayerSourceError::Invalid("tensor shear thermal capacity bound refused")
            })?;
        let mass = projected.volume_m3() * material.density_kg_m3();
        for segment in knots.windows(2) {
            let capacity = mass
                * (segment[1].specific_enthalpy_j_kg - segment[0].specific_enthalpy_j_kg)
                / (segment[1].temperature_k - segment[0].temperature_k)
                - segment[1].temperature_k * capacity_bound;
            if !capacity.is_finite() || capacity <= 0.0 {
                return Err(TensorLayerSourceError::Invalid(
                    "tensor layer cannot establish positive total fixed-strain heat capacity over its original domain",
                ));
            }
        }
        let thermal_slopes = project_thermal_slopes(&projected, &derivative, projection_budget)?;
        let requirements = [
            ScalarPropertyRequirement::try_new(
                SHEAR_LAYER_YIELD_PROPERTY,
                Pressure::DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .map_err(TensorLayerSourceError::Material)?,
            ScalarPropertyRequirement::try_new(
                PLASTIC_SHEAR_VISCOSITY_PROPERTY,
                DynViscosity::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .map_err(TensorLayerSourceError::Material)?,
        ];
        let plastic = resolve_material_state_point(
            card,
            &material.state().query(),
            &requirements,
            plastic_selection,
        )
        .map_err(TensorLayerSourceError::Material)?;
        admit_plastic(&plastic, card, temperature_axis, interval)?;
        let value = Self {
            source: TensorLayerSourceState { material, plastic },
            card,
            caloric,
            temperature_axis: temperature_axis.into(),
            material_to_world,
            geometry,
            maximum_angular_frequency_rad_s,
            projection_budget,
            thermal_slopes,
            thermal_derivative: derivative,
            thermal_projection: projected,
        };
        value.at_temperature(temperature)?;
        Ok(value)
    }
    pub const fn source(&self) -> &TensorLayerSourceState {
        &self.source
    }
    pub const fn geometry(&self) -> RigidLayerGeometry {
        self.geometry
    }
    /// Source derivatives at the ORIGINAL resolved temperature, after the same
    /// material-frame and physical shear projection, in N/(m K). For a nonlinear
    /// source these are reference metadata; use `thermal_slopes_at_temperature`
    /// for another temperature.
    pub const fn thermal_slopes(&self) -> &TensorLayerThermalSlopes {
        &self.thermal_slopes
    }

    /// Complete derivatives at an actual admitted endpoint temperature. The
    /// original source paths and orientation are already admitted; refresh only
    /// the shared thermal snapshot and project it through the fixed geometry.
    pub fn thermal_slopes_at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<TensorLayerThermalSlopes, TensorLayerSourceError> {
        self.resolve_thermal_slopes(temperature_k)
    }

    fn resolve_thermal_slopes<E>(
        &self,
        temperature_k: f64,
    ) -> Result<TensorLayerThermalSlopes, TensorLayerSourceError<E>> {
        let derivative = self
            .thermal_derivative
            .at_temperature(temperature_k)
            .map_err(TensorLayerSourceError::Tensor)?;
        project_thermal_slopes(
            &self.thermal_projection,
            &derivative,
            self.projection_budget,
        )
    }

    /// Requery original claims, rotate the full 3-D law, and project every arm.
    /// The resulting law initializes/continues the existing retained history.
    pub fn at_temperature(
        &self,
        temperature_k: f64,
    ) -> Result<(TensorLayerSourceState, ThirdBodyLayerContact), TensorLayerSourceError> {
        self.resolve(temperature_k)
    }
    fn resolve<E>(
        &self,
        temperature_k: f64,
    ) -> Result<(TensorLayerSourceState, ThirdBodyLayerContact), TensorLayerSourceError<E>> {
        let knots = self.caloric.curve().knots();
        if !temperature_k.is_finite()
            || temperature_k < knots[0].temperature_k
            || temperature_k > knots[knots.len() - 1].temperature_k
        {
            return Err(TensorLayerSourceError::Invalid(
                "layer temperature outside original caloric interval",
            ));
        }
        let query = self
            .source
            .material
            .state()
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| {
                TensorLayerSourceError::Invalid("layer temperature violates original query schema")
            })?;
        let material = self
            .source
            .material
            .at_query_point(self.card, &query)
            .map_err(TensorLayerSourceError::Tensor)?;
        let plastic = self
            .source
            .plastic
            .at_query_point(self.card, &query)
            .map_err(TensorLayerSourceError::Material)?;
        let band = material.law().band_rad_s();
        if band.0 != 0.0 || band.1 < self.maximum_angular_frequency_rad_s {
            return Err(TensorLayerSourceError::Invalid(
                "full tensor source must cover static storage and actual band",
            ));
        }
        for property in plastic.properties() {
            let band = property
                .angular_frequency_band_rad_s()
                .map_err(TensorLayerSourceError::Material)?;
            if !property.is_constant_at_fixed_temperature()
                || band.is_none_or(|b| b.0 != 0.0 || b.1 < self.maximum_angular_frequency_rad_s)
            {
                return Err(TensorLayerSourceError::Invalid(
                    "plastic source cannot freeze evolving response or omit actual band",
                ));
            }
        }
        let h = self.geometry.thickness_m;
        let projected = project_layer(&material, self.material_to_world,
            self.geometry, self.projection_budget)?;
        let k = projected.equilibrium_stiffness();
        let equilibrium = ShearStiffness::positive_definite([[k[0], k[1]], [k[2], k[3]]])
            .map_err(TensorLayerSourceError::Constitutive)?;
        let arms = projected
            .branches()
            .iter()
            .map(|branch| {
                Ok(TensorLayerRelaxationArm {
                    stiffness: ShearStiffness::from_gram_rows(
                        branch.rows().iter().map(|row| [row[0], row[1]]).collect(),
                    )
                    .map_err(TensorLayerSourceError::Constitutive)?,
                    relaxation_time_s: branch.relaxation_time_s(),
                })
            })
            .collect::<Result<Vec<_>, TensorLayerSourceError<E>>>()?;
        let yield_force = self.geometry.area_m2
            * plastic
                .property(SHEAR_LAYER_YIELD_PROPERTY)
                .expect("admitted plastic source")
                .value_si();
        let resistance = self.geometry.area_m2 / h
            * plastic
                .property(PLASTIC_SHEAR_VISCOSITY_PROPERTY)
                .expect("admitted plastic source")
                .value_si();
        // Frobenius norm of this pure shear is |e|/(sqrt(2)*h), independent
        // of orientation. Plastic slip is not bounded by elastic small strain.
        let slip_limit = core::f64::consts::SQRT_2 * h * material.law().maximum_strain_norm();
        let mut law = ThirdBodyLayerContact::new_tensor(
            equilibrium,
            arms,
            yield_force,
            resistance,
            temperature_k,
            slip_limit,
        )
        .map_err(TensorLayerSourceError::Constitutive)?;
        if let Some(rows) = projected.viscosity_gram_rows() {
            let viscosity =
                ShearViscosity::from_gram_rows(rows.iter().map(|row| [row[0], row[1]]).collect())
                    .map_err(TensorLayerSourceError::Constitutive)?;
            law = law
                .with_kelvin_viscosity(viscosity)
                .map_err(TensorLayerSourceError::Constitutive)?;
        }
        Ok((TensorLayerSourceState { material, plastic }, law))
    }
    fn admit<E>(
        &self,
        history: &RigidLayerHistory,
        thermal: &DepositionSimulation<'_>,
        route: &ThirdBodyHeatRoute,
    ) -> Result<(), TensorLayerSourceError<E>> {
        let cell =
            thermal
                .cells()
                .get(&route.layer_cell())
                .ok_or(TensorLayerSourceError::Invalid(
                    "inactive original layer cell",
                ))?;
        let volume = self.geometry.area_m2 * self.geometry.thickness_m;
        let same = |a: f64, b: f64| {
            a.is_finite()
                && b.is_finite()
                && a > 0.0
                && b > 0.0
                && (a - b).abs() <= 512.0 * f64::EPSILON * a.max(b)
        };
        if history.geometry != self.geometry
            || !thermal.uses_material_law(cell.material_id, self.caloric.curve() as &dyn ThermalLaw)
            || !same(cell.geometry.volume_m3, volume)
            || !same(cell.mass_kg, volume * self.source.material.density_kg_m3())
        {
            return Err(TensorLayerSourceError::Invalid(
                "original tensor geometry, caloric source and material mass must agree",
            ));
        }
        let (_, original) = self.resolve(history.law.temperature_kelvin())?;
        if history.law != original {
            return Err(TensorLayerSourceError::Invalid(
                "retained layer belongs to different source or orientation",
            ));
        }
        Ok(())
    }
}

fn project_thermal_slopes<E>(
    projected: &ProjectedTensorProny,
    derivative: &TensorTemperatureDerivative,
    budget: TensorProjectionBudget,
) -> Result<TensorLayerThermalSlopes, TensorLayerSourceError<E>> {
    let derivative = projected
        .thermal_stiffness_derivatives(derivative, budget, || false)
        .map_err(|_| TensorLayerSourceError::Invalid("tensor shear thermal projection refused"))?;
    let matrix = |values: &[f64]| [[values[0], values[1]], [values[2], values[3]]];
    Ok(TensorLayerThermalSlopes {
        equilibrium_n_m_k: matrix(derivative.equilibrium()),
        relaxation_n_m_k: derivative
            .branches()
            .iter()
            .map(|values| matrix(values))
            .collect(),
    })
}

fn project_layer<E>(
    material: &ResolvedTensorMaterial,
    material_to_world: [[f64; 3]; 3],
    geometry: RigidLayerGeometry,
    budget: TensorProjectionBudget,
) -> Result<ProjectedTensorProny, TensorLayerSourceError<E>> {
    let rotated = material.rotated_law(material_to_world)
        .map_err(TensorLayerSourceError::Tensor)?;
    let normal = geometry.frame.normal();
    let tangent = [geometry.frame.longitudinal(), geometry.frame.lateral()];
    let h = geometry.thickness_m;
    let columns = tangent.map(|t| {
        let strain = |i: usize, j: usize| (t[i]*normal[j] + normal[i]*t[j])/(2.0*h);
        [strain(0,0), strain(1,1), strain(2,2), strain(0,1), strain(1,2), strain(2,0)]
    });
    rotated.project_strain(&columns, geometry.area_m2*h, budget, || false)
        .map_err(|_| TensorLayerSourceError::Invalid("tensor shear projection refused"))
}

fn admit_plastic(
    state: &ResolvedMaterialStatePoint,
    card: &MaterialCard,
    axis: &str,
    interval: [f64; 2],
) -> Result<(), TensorLayerSourceError> {
    for property in state.properties() {
        if !property.has_only_fixed_external_conditions(Some(axis), true)
            || property
                .answer()
                .evidence
                .model
                .validity
                .bound(axis)
                .is_none()
        {
            return Err(TensorLayerSourceError::Invalid(
                "plastic law requires explicit temperature and fixed external source conditions",
            ));
        }
        let claim = card
            .claims()
            .claim(property.answer().receipt.selected)
            .ok_or(TensorLayerSourceError::Invalid(
                "missing original plastic claim",
            ))?;
        match (&claim.value, claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {}
            (
                PropertyValue::Curve {
                    abscissa,
                    abscissa_dims,
                    knots,
                    ..
                },
                InterpolationPolicy::LinearInside,
            ) if abscissa == axis && *abscissa_dims == Temperature::DIMS => {
                let yield_law = property.requirement().name() == SHEAR_LAYER_YIELD_PROPERTY;
                if knots.len() < 2
                    || knots[0].0 > interval[0]
                    || knots[knots.len() - 1].0 < interval[1]
                    || knots.iter().any(|&(t, v)| {
                        (interval[0]..=interval[1]).contains(&t)
                            && (!v.is_finite() || v < 0.0 || (!yield_law && v == 0.0))
                    })
                {
                    return Err(TensorLayerSourceError::Invalid(
                        "plastic curve must cover whole interval with admissible knots",
                    ));
                }
            }
            _ => {
                return Err(TensorLayerSourceError::Invalid(
                    "unsupported plastic source temperature law",
                ));
            }
        }
    }
    for temperature in interval {
        let query = state.query().with(axis, temperature).map_err(|_| {
            TensorLayerSourceError::Invalid("plastic source temperature schema refused")
        })?;
        state
            .at_query_point(card, &query)
            .map_err(TensorLayerSourceError::Material)?;
    }
    Ok(())
}

/// Actual tensor material drives the same finite-face force/couple, inertia,
/// retained strain, and heat solve. All source trials retain original receipts.
/// The original retained strains survive tensor coefficient changes. Caloric heat
/// is plastic + Maxwell + Kelvin loss minus coefficient free-energy change
/// minus delta(T*S); only the actual dissipative channels produce friction entropy.
#[allow(clippy::too_many_arguments)]
pub fn step_tensor_heated_rigid_layer_candidate<'t, P, E>(
    bodies: &[RigidBodyState; 2],
    properties: [MassProperties; 2],
    history: &RigidLayerHistory,
    thermal: &DepositionSimulation<'t>,
    route: &ThirdBodyHeatRoute,
    material: &TensorLayerMaterial<'_>,
    target_s: f64,
    normal_n: f64,
    law_solve: TangentialSolveConfig,
    controls: HeatedRigidLayerSolve,
    force_guess_n: [f64; 2],
    energy_guess_j: &[f64],
    resolve_transport: impl Fn(
        &DepositionSimulation<'t>,
        &TensorLayerSourceState,
    ) -> Result<RigidLayerThermalTransport<P>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<HeatedRigidLayerCandidate<'t, P>, HeatedRigidLayerError<TensorLayerSourceError<E>>> {
    material
        .admit(history, thermal, route)
        .map_err(HeatedRigidLayerError::Source)?;
    let nonlinear = !material.thermal_derivative.is_affine();
    let has_thermal_storage = nonlinear || material.thermal_slopes.equilibrium_n_m_k.iter().flatten()
        .chain(material.thermal_slopes.relaxation_n_m_k.iter().flatten().flatten())
        .any(|slope| *slope != 0.0);
    let endpoint_slopes = |temperature| material.resolve_thermal_slopes(temperature);
    let slopes = has_thermal_storage.then(|| {
        if nonlinear {
            LayerThermalSlopes::TensorAt(&endpoint_slopes)
        } else {
            LayerThermalSlopes::Tensor(&material.thermal_slopes)
        }
    });
    step_heated_impl(
        bodies,
        properties,
        history,
        thermal,
        route,
        target_s,
        normal_n,
        law_solve,
        controls,
        force_guess_n,
        energy_guess_j,
        slopes,
        |endpoint, temperature| {
            let (source, law) = material.resolve(temperature)?;
            let transport =
                resolve_transport(endpoint, &source).map_err(TensorLayerSourceError::Transport)?;
            Ok(RigidLayerHeatSource {
                law,
                conductivity_w_m_k: transport.conductivity_w_m_k,
                area_resistance_m2_k_w: transport.area_resistance_m2_k_w,
                transport: transport.transport,
            })
        },
        cancelled,
    )
}
