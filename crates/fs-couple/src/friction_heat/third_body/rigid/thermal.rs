//! Reciprocal source-defined shear softening in the native rigid/heat owner.
//!
//! The original material supplies affine k(T)=A*G(T)/h. Its retained potential
//! is psi=k(T)|e|^2/2 plus every affine Maxwell spring storage term;
//! s=-k'|e|^2/2-sum(k_j'|r_j|^2/2) and U=psi+T*s. The joint residual deposits
//! D-delta_psi_coeff-delta(T*s) into the existing layer caloric cell. D alone
//! produces friction entropy. The source, geometry, strain history, physical
//! clock, face wrenches and thermal cells retain their original owners.
//!
//! This is a fully solid, fixed-geometry, small elastic shear reduction.
//! It does not infer phase change, expansion, nonlinear modulus laws, or any
//! measured coefficient. The existing material source rejects unsupported laws.

use super::*;
use fs_material::phase::ResolvedHeatCapacityCurve;
use fs_material::state_point::{MaterialStatePointError, ResolvedShearLayer, ShearLayerThermalSlopes};
use fs_matdb::MaterialCard;
use fs_thermal::enthalpy::ThermalLaw;
use fs_tribo::layer_contact::ThirdBodyLayerRelaxationArm;
use std::convert::Infallible;

#[derive(Debug)]
pub enum ThermoelasticLayerError<E = Infallible> {
    Material(MaterialStatePointError),
    Constitutive(TriboError),
    Transport(E),
}
impl<E: core::fmt::Debug> core::fmt::Display for ThermoelasticLayerError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-defined thermoelastic layer: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ThermoelasticLayerError<E> {}

fn invalid<E>(quantity: &'static str) -> ThermoelasticLayerError<E> {
    ThermoelasticLayerError::Material(MaterialStatePointError::InvalidDerived { quantity })
}
fn same_measure(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a-b).abs() <= 512.0*f64::EPSILON*a.abs().max(b.abs())
}

/// Immutable original-card binding. Temperature changes re-query the complete
/// pinned material bundle; an independent authored spring slope is not accepted.
pub struct ThermoelasticLayerMaterial<'a> {
    source: ResolvedShearLayer,
    card: &'a MaterialCard,
    caloric: &'a ResolvedHeatCapacityCurve,
    temperature_axis: String,
    thermal_slopes: ShearLayerThermalSlopes,
    maximum_angular_frequency_rad_s: f64,
    density_kg_m3: f64,
}
impl<'a> ThermoelasticLayerMaterial<'a> {
    /// Admit the complete original single-solid caloric interval and constant
    /// reference density. The source owner admits the complete affine equilibrium
    /// and Maxwell potential with original claim identities; kinetics may vary.
    pub fn try_new(source: ResolvedShearLayer, card: &'a MaterialCard,
        caloric: &'a ResolvedHeatCapacityCurve, temperature_axis: &str,
        maximum_angular_frequency_rad_s: f64) -> Result<Self, ThermoelasticLayerError> {
        let curve = caloric.curve();
        let knots = curve.knots();
        if curve.material_card_identity() != source.resolved().card_identity()
            || curve.material_card_identity() != card.content_hash()
            || knots.len() < 2 || knots.iter().any(|k| k.liquid_mass_fraction != 0.0)
            || knots.windows(2).any(|p| p[0].temperature_k >= p[1].temperature_k)
            || knots.iter().any(|k| k.bulk_density_kg_m3 != knots[0].bulk_density_kg_m3) {
            return Err(invalid("layer mechanics and caloric storage need one fully solid, fixed-density original card"));
        }
        source.require_angular_frequency_rad_s(maximum_angular_frequency_rad_s)
            .map_err(ThermoelasticLayerError::Material)?;
        let thermal_slopes = source.thermal_stiffness_slopes_n_m_k(card,
            temperature_axis, (knots[0].temperature_k, knots[knots.len()-1].temperature_k))
            .map_err(ThermoelasticLayerError::Material)?;
        let value = Self { source, card, caloric, temperature_axis: temperature_axis.into(),
            thermal_slopes, maximum_angular_frequency_rad_s,
            density_kg_m3: knots[0].bulk_density_kg_m3 };
        value.at_temperature(value.source.temperature_kelvin())?;
        Ok(value)
    }
    pub const fn source(&self) -> &ResolvedShearLayer { &self.source }
    pub const fn thermal_slopes(&self) -> &ShearLayerThermalSlopes { &self.thermal_slopes }

    /// Resolve the same original claims at an actual caloric temperature. This
    /// also supplies the native law used to initialize the retained history.
    pub fn at_temperature(&self, temperature_k: f64)
        -> Result<(ResolvedShearLayer, ThirdBodyLayerContact), ThermoelasticLayerError> {
        self.resolve(temperature_k)
    }
    fn resolve<E>(&self, temperature_k: f64)
        -> Result<(ResolvedShearLayer, ThirdBodyLayerContact), ThermoelasticLayerError<E>> {
        let knots = self.caloric.curve().knots();
        if !temperature_k.is_finite() || temperature_k < knots[0].temperature_k
            || temperature_k > knots[knots.len()-1].temperature_k {
            return Err(invalid("layer temperature left its original caloric source interval"));
        }
        let query = self.source.resolved().query().with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("layer temperature violates its original query schema"))?;
        let source = self.source.at_query_point(self.card, &query)
            .map_err(ThermoelasticLayerError::Material)?;
        source.require_angular_frequency_rad_s(self.maximum_angular_frequency_rad_s)
            .map_err(ThermoelasticLayerError::Material)?;
        let law = ThirdBodyLayerContact::new(source.stiffness_n_m(), source.yield_force_n(),
            source.plastic_resistance_n_s_m(), source.temperature_kelvin(), source.maximum_elastic_slip_m())
            .and_then(|law| law.with_relaxation(source.relaxation().iter().map(|&(stiffness_n_m, relaxation_time_s)|
                ThirdBodyLayerRelaxationArm { stiffness_n_m, relaxation_time_s }).collect()))
            .map_err(ThermoelasticLayerError::Constitutive)?;
        Ok((source, law))
    }

    fn admit<E>(&self, history: &RigidLayerHistory, thermal: &DepositionSimulation<'_>,
        route: &ThirdBodyHeatRoute) -> Result<(), ThermoelasticLayerError<E>> {
        let cell = thermal.cells().get(&route.layer_cell())
            .ok_or_else(|| invalid("original finite layer cell is inactive"))?;
        if !thermal.uses_material_law(cell.material_id, self.caloric.curve() as &dyn ThermalLaw)
            || !same_measure(history.geometry.area_m2, self.source.area_m2())
            || !same_measure(history.geometry.thickness_m, self.source.thickness_m())
            || !same_measure(cell.geometry.volume_m3, self.source.area_m2()*self.source.thickness_m())
            || !same_measure(cell.mass_kg, self.density_kg_m3*cell.geometry.volume_m3) {
            return Err(invalid("layer mechanics must retain the original caloric law, material mass and source geometry"));
        }
        let (_, original) = self.resolve(history.law.temperature_kelvin())?;
        if original != history.law {
            return Err(invalid("accepted layer history does not belong to these original material claims"));
        }
        Ok(())
    }
}

/// Other heat and the two original finite-volume transport closures, resolved
/// alongside the complete current shear material. This input supplies no second
/// mechanical law, coefficient slope or authored share of frictional heat.
pub struct RigidLayerThermalTransport<P> {
    pub conductivity_w_m_k: [f64; 3],
    pub area_resistance_m2_k_w: [f64; 2],
    pub transport: StepCoefficients<P>,
}

/// Evolve actual body momenta, retained shear/plastic/memory history and finite
/// caloric cells in the original joint solve. All trials start at the same
/// accepted state; failed source, nonlinear, work or heat checks publish nothing.
#[allow(clippy::too_many_arguments)]
pub fn step_thermoelastic_heated_rigid_layer_candidate<'t, P, E>(
    bodies: &[RigidBodyState; 2], properties: [MassProperties; 2],
    history: &RigidLayerHistory, thermal: &DepositionSimulation<'t>,
    route: &ThirdBodyHeatRoute, material: &ThermoelasticLayerMaterial<'_>,
    target_s: f64, normal_n: f64, law_solve: TangentialSolveConfig,
    controls: HeatedRigidLayerSolve, force_guess_n: [f64; 2], energy_guess_j: &[f64],
    resolve_transport: impl Fn(&DepositionSimulation<'t>, &ResolvedShearLayer)
        -> Result<RigidLayerThermalTransport<P>, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<HeatedRigidLayerCandidate<'t, P>, HeatedRigidLayerError<ThermoelasticLayerError<E>>> {
    material.admit(history, thermal, route).map_err(HeatedRigidLayerError::Source)?;
    step_heated_impl(bodies, properties, history, thermal, route, target_s, normal_n,
        law_solve, controls, force_guess_n, energy_guess_j,
        Some(LayerThermalSlopes::Scalar(&material.thermal_slopes)),
        |endpoint, temperature| {
            let (source, law) = material.resolve(temperature)?;
            let resolved = resolve_transport(endpoint, &source).map_err(ThermoelasticLayerError::Transport)?;
            Ok(RigidLayerHeatSource { law, conductivity_w_m_k: resolved.conductivity_w_m_k,
                area_resistance_m2_k_w: resolved.area_resistance_m2_k_w,
                transport: resolved.transport })
        }, cancelled)
}
