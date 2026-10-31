//! Reversible plane-stress potentials of the original layer sources.
//!
//! In the fixed principal frame, engineering strain is `[11,22,gamma12,
//! gamma13,gamma23]` and its work-conjugate stress is `[sigma11,sigma22,
//! tau12,tau13,tau23]`. The original integrated CTE supplies the first two
//! free strains. No thickness-normal elastic constant or strain is invented.
//! Independent affine source constants generally give rational Q(T); exact
//! second-order arithmetic differentiates that reduction before computing
//! `S=-Psi_T`, `U=Psi+T*S`, and `C_epsilon=-T*Psi_TT` at fixed TOTAL strain.
//!
//! Transverse shear here is the UNCORRECTED local G13/G23 potential. A native
//! laminate's Whitney correction depends on the entire stack and belongs to
//! its section owner. `in_plane_state` excludes that energy. A consumer must
//! account for the original corrected shear separately; replacing it by a
//! constant 5/6 per layer is incorrect for a general heterogeneous stack.
//!
//! All quantities use fixed reference volume/geometry/mass. These elastic
//! contributions supplement the original caloric chart; they are not total
//! heat capacity, a thermal integrator, or viscoelastic-history evolution.

use super::source::{
    BoundLaminaSource, IntegratedLaminaExpansion, LaminaElasticSource, LaminaExpansionSource,
    ResolvedLaminaThermal,
};
use crate::elastic::thermal::ThermoelasticError;
use crate::field_storage::{
    FieldStorageError, admit_source_temperature_refresh, thermal::affine_source_path,
};
use crate::state_point::{
    INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
    ORTHOTROPIC_POISSON_RATIO_PROPERTIES as NU, ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES as G,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE3,
    ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES as E, ResolvedIsotropicElasticModuli,
    ResolvedMaterialStatePoint, integrate_isotropic_thermal_expansion,
    integrate_orthotropic_thermal_expansion, pinned_temperature_line,
};
use fs_qty::{Dims, Pressure, QuantitySpec};
use std::ops::{Add, Div, Mul, Neg, Sub};

fn invalid(message: &'static str) -> ThermoelasticError {
    ThermoelasticError::Invalid(message)
}
fn field(error: FieldStorageError) -> ThermoelasticError {
    match error {
        FieldStorageError::Material(e) => ThermoelasticError::Source(e),
        FieldStorageError::Invalid(s) => invalid(s),
    }
}
fn finite(x: f64) -> Result<f64, ThermoelasticError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("plane thermoelastic arithmetic is unrepresentable"))
    }
}
fn value(
    state: &ResolvedMaterialStatePoint,
    name: &str,
    dims: Dims,
) -> Result<f64, ThermoelasticError> {
    let p = state.property(name).ok_or(invalid(
        "original layer coefficient is missing or ambiguous",
    ))?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims) {
        return Err(invalid(
            "plane thermoelasticity retains dimension-only original scalar roles",
        ));
    }
    finite(p.value_si())
}

#[derive(Clone, Debug)]
struct Coefficient {
    name: String,
    dims: Dims,
    reference: f64,
    slope: f64,
    endpoints: [f64; 2],
    constant_over_interval: bool,
}
impl Coefficient {
    fn new(
        source: &BoundLaminaSource,
        state: &ResolvedMaterialStatePoint,
        name: &str,
        dims: Dims,
        positive: bool,
        range: [f64; 2],
    ) -> Result<Self, ThermoelasticError> {
        let reference = value(state, name, dims)?;
        let (_, slope, endpoints) =
            affine_source_path(&source.input().card, state, name, "T", range, positive)
                .map_err(field)?;
        if slope == 0.0 && endpoints[0] != endpoints[1] {
            return Err(invalid(
                "nonconstant layer coefficient derivative underflows",
            ));
        }
        // Affine admission allows representation-scale residuals. A fixed
        // native operator needs the stronger statement that every original
        // support-knot value is exactly constant, not just a zero fitted slope.
        let constant_over_interval = if slope == 0.0 {
            let property = state.property(name).expect("coefficient was resolved");
            let selected = property.answer().receipt.selected;
            let claim = source
                .input()
                .card
                .claims()
                .claim(selected)
                .ok_or(invalid("original layer coefficient claim is absent"))?;
            let point = |t| {
                state
                    .query()
                    .with("T", t)
                    .map_err(|_| invalid("layer temperature convention is invalid"))
            };
            pinned_temperature_line(
                &source.input().card,
                &claim.key,
                "T",
                &point(range[0])?,
                &point(range[1])?,
                selected,
            )?
            .iter()
            .all(|(_, value)| *value == reference)
        } else {
            false
        };
        Ok(Self {
            name: name.to_owned(),
            dims,
            reference,
            slope,
            endpoints,
            constant_over_interval,
        })
    }
    fn at(
        &self,
        state: &ResolvedMaterialStatePoint,
        delta: f64,
    ) -> Result<Jet<f64>, ThermoelasticError> {
        let actual = value(state, &self.name, self.dims)?;
        let expected = finite(self.slope.mul_add(delta, self.reference))?;
        let scale = self
            .endpoints
            .iter()
            .fold(actual.abs().max(expected.abs()), |m, x| m.max(x.abs()));
        if (actual - expected).abs() > 512.0 * f64::EPSILON * scale {
            return Err(invalid(
                "actual layer coefficient leaves the admitted affine source path",
            ));
        }
        Ok(Jet {
            v: actual,
            d: self.slope,
            dd: 0.0,
        })
    }
    fn bound(&self) -> Jet<Interval> {
        Jet {
            v: Interval::range(
                self.endpoints[0].min(self.endpoints[1]),
                self.endpoints[0].max(self.endpoints[1]),
            ),
            d: Interval::from(self.slope),
            dd: Interval::from(0.0),
        }
    }
}

/// Smooth potential on a complete original-source temperature interval.
/// This shares the original card through the bound source's Arc; no source
/// card, measured coefficient or derived three-dimensional law is fabricated.
#[derive(Clone, Debug)]
pub struct AffineLaminaThermoelastic {
    source: BoundLaminaSource,
    range: [f64; 2],
    reference_temperature: f64,
    limit: f64,
    elastic: Vec<Coefficient>,
    expansion: Vec<Coefficient>,
    curvature_bound: f64,
}
impl AffineLaminaThermoelastic {
    /// The explicit component envelope applies to BOTH total and elastic
    /// principal engineering strains. Temperature must be typed absolute `T`.
    /// Each selected coefficient must be affine on the whole interval; source
    /// gaps and slope changes refuse, without an invented transition entropy.
    /// Interval arithmetic may conservatively refuse a poorly bounded rational
    /// denominator; callers may supply a narrower supported interval.
    pub fn new(
        source: &BoundLaminaSource,
        temperature_range_k: [f64; 2],
        maximum_absolute_engineering_strain: f64,
    ) -> Result<Self, ThermoelasticError> {
        let limit = maximum_absolute_engineering_strain;
        if !limit.is_finite() || limit <= 0.0 || limit >= 1.0 {
            return Err(invalid(
                "a principal engineering-strain component bound in (0,1) is required",
            ));
        }
        let schema: Vec<_> = match &source.input().elasticity {
            LaminaElasticSource::Isotropic(p) => p
                .requirements()
                .into_iter()
                .map(|r| (r.name(), r.dims(), r.dims() == Pressure::DIMS))
                .collect(),
            LaminaElasticSource::OrthotropicPlane12 => vec![
                (E[0], Pressure::DIMS, true),
                (E[1], Pressure::DIMS, true),
                (NU[0], Dims::NONE, false),
                (G[0], Pressure::DIMS, true),
                (G[2], Pressure::DIMS, true),
                (G[1], Pressure::DIMS, true),
            ],
        };
        let elastic = schema
            .into_iter()
            .map(|(name, dims, positive)| {
                Coefficient::new(
                    source,
                    source.reference_elasticity(),
                    name,
                    dims,
                    positive,
                    temperature_range_k,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        // alpha3 remains in the original integrated path but supplies no work
        // in this reduced plane. Its unmeasured elastic counterpart is absent.
        let names: &[&str] = match source.input().expansion {
            LaminaExpansionSource::Isotropic => &[CTE],
            LaminaExpansionSource::Orthotropic => &CTE3[..2],
        };
        let expansion = names
            .iter()
            .map(|name| {
                Coefficient::new(
                    source,
                    source.reference_expansion().reference(),
                    name,
                    INVERSE_TEMPERATURE_DIMS,
                    false,
                    temperature_range_k,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let reference_temperature = source
            .reference_elasticity()
            .state_coordinate("T")
            .ok_or(invalid("original layer temperature is absent"))?;
        let mut law = Self {
            source: source.clone(),
            range: temperature_range_k,
            reference_temperature,
            limit,
            elastic,
            expansion,
            curvature_bound: 0.0,
        };
        law.admit_isotropic_interval()?;
        law.curvature_bound = law.bound_curvature()?;
        // Original endpoint/path owners retain full support, positivity and
        // consumer admissibility; density is never part of these refreshes.
        for t in temperature_range_k {
            law.at_temperature(t)?;
        }
        Ok(law)
    }
    /// Original layer card, reference geometry/mass, selections and receipts.
    pub const fn source(&self) -> &BoundLaminaSource {
        &self.source
    }
    /// Complete admitted smooth temperature interval [K].
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.range
    }
    /// Original stress-free expansion reference temperature [K].
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature
    }
    /// Bound on each total AND elastic principal engineering-strain component.
    pub const fn maximum_absolute_engineering_strain(&self) -> f64 {
        self.limit
    }
    /// Certifies fixed Q AND both transverse moduli over the whole admitted
    /// source interval, using exact zero source slopes rather than samples.
    pub fn is_temperature_independent_elasticity(&self) -> bool {
        self.elastic.iter().all(|c| c.constant_over_interval)
    }
    /// Conservative elastic contribution to fixed-strain heat capacity over
    /// the complete interval/envelope. Add to the actual caloric capacity and
    /// require a positive total before coupling. Includes raw local shear;
    /// it does not bound a temperature-dependent Whitney section correction.
    pub fn minimum_heat_capacity_j_k(
        &self,
        reference_volume_m3: f64,
    ) -> Result<f64, ThermoelasticError> {
        if !reference_volume_m3.is_finite() || reference_volume_m3 <= 0.0 {
            return Err(invalid(
                "positive finite reference layer volume is required",
            ));
        }
        Ok(-upper_mul(
            upper_mul(self.range[1], reference_volume_m3)?,
            self.curvature_bound,
        )?)
    }
    /// Resolve the original claims at one actual temperature. Runtime source
    /// data and the signed expansion integral are retained unchanged.
    pub fn at_temperature(
        &self,
        t: f64,
    ) -> Result<ResolvedLaminaThermoelastic, ThermoelasticError> {
        self.admit_temperature(t)?;
        let input = self.source.input();
        let point = input
            .reference
            .clone()
            .with("T", t)
            .map_err(|_| invalid("layer temperature convention is invalid"))?;
        let elastic = self
            .source
            .reference_elasticity()
            .at_query_point(&input.card, &point)?;
        let selected = MaterialPropertySelection::PinnedByProperty(
            self.source
                .reference_expansion()
                .reference()
                .properties()
                .iter()
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let expansion = match input.expansion {
            LaminaExpansionSource::Isotropic => {
                IntegratedLaminaExpansion::Isotropic(integrate_isotropic_thermal_expansion(
                    &input.card,
                    &input.reference,
                    &point,
                    selected,
                )?)
            }
            LaminaExpansionSource::Orthotropic => {
                IntegratedLaminaExpansion::Orthotropic(integrate_orthotropic_thermal_expansion(
                    &input.card,
                    &input.reference,
                    &point,
                    selected,
                )?)
            }
        };
        self.resolve(elastic, expansion)
    }
    /// Adopt the accepted native laminate snapshot, checking exact original
    /// claims, reference expansion path and current context. This performs no
    /// material query or source reselection; actual resolved Q values govern.
    pub fn at_materials(
        &self,
        current: &ResolvedLaminaThermal,
    ) -> Result<ResolvedLaminaThermoelastic, ThermoelasticError> {
        self.resolve(current.elasticity().clone(), current.expansion().clone())
    }
    fn admit_temperature(&self, t: f64) -> Result<(), ThermoelasticError> {
        if !t.is_finite() || !(self.range[0]..=self.range[1]).contains(&t) {
            Err(invalid(
                "temperature leaves the admitted plane thermoelastic interval",
            ))
        } else {
            Ok(())
        }
    }
    fn resolve(
        &self,
        elastic: ResolvedMaterialStatePoint,
        expansion: IntegratedLaminaExpansion,
    ) -> Result<ResolvedLaminaThermoelastic, ThermoelasticError> {
        admit_source_temperature_refresh(self.source.reference_elasticity(), &elastic)
            .map_err(field)?;
        admit_source_temperature_refresh(
            self.source.reference_expansion().reference(),
            expansion.current(),
        )
        .map_err(field)?;
        if expansion.reference() != self.source.reference_expansion().reference()
            || elastic.query() != expansion.current().query()
        {
            return Err(invalid(
                "native layer potential must retain its actual reference path and common current context",
            ));
        }
        let t = elastic
            .state_coordinate("T")
            .ok_or(invalid("actual layer temperature is absent"))?;
        self.admit_temperature(t)?;
        let delta = t - self.reference_temperature;
        let values = self
            .elastic
            .iter()
            .map(|c| c.at(&elastic, delta))
            .collect::<Result<Vec<_>, _>>()?;
        let (mut q, mut g) = plane(&self.source.input().elasticity, &values);
        let (native_q, native_g) =
            native_values(&self.source.input().elasticity, &elastic, &values)?;
        for (jet, value) in q.iter_mut().zip(native_q) {
            jet.v = value;
        }
        for (jet, value) in g.iter_mut().zip(native_g) {
            jet.v = value;
        }
        if q.iter()
            .chain(g.iter())
            .any(|j| ![j.v, j.d, j.dd].iter().all(|v| v.is_finite()))
        {
            return Err(invalid(
                "plane stiffness temperature derivatives are unrepresentable",
            ));
        }
        let alpha = self
            .expansion
            .iter()
            .map(|c| c.at(expansion.current(), delta))
            .collect::<Result<Vec<_>, _>>()?;
        let free = match &expansion {
            IntegratedLaminaExpansion::Isotropic(p) => [p.free_linear_strain(); 2],
            IntegratedLaminaExpansion::Orthotropic(p) => {
                [p.free_principal_strains()[0], p.free_principal_strains()[1]]
            }
        };
        let alpha = core::array::from_fn(|i| alpha[if alpha.len() == 1 { 0 } else { i }]);
        Ok(ResolvedLaminaThermoelastic {
            elastic,
            expansion,
            temperature: t,
            limit: self.limit,
            q,
            g,
            free,
            alpha,
        })
    }
    fn admit_isotropic_interval(&self) -> Result<(), ThermoelasticError> {
        if let LaminaElasticSource::Isotropic(parameters) = &self.source.input().elasticity {
            match parameters {
                IsotropicElasticParameters::YoungPoisson { .. } => {
                    if self.elastic[1]
                        .endpoints
                        .iter()
                        .any(|nu| *nu <= -1.0 || *nu >= 0.5)
                    {
                        return Err(invalid(
                            "isotropic Poisson source must stay in (-1,1/2) throughout the interval",
                        ));
                    }
                }
                IsotropicElasticParameters::YoungShear { .. } => {
                    for i in 0..2 {
                        if 3.0_f64
                            .mul_add(self.elastic[1].endpoints[i], -self.elastic[0].endpoints[i])
                            .next_down()
                            <= 0.0
                        {
                            return Err(invalid(
                                "isotropic E/G interval must retain positive 3G-E",
                            ));
                        }
                    }
                }
                IsotropicElasticParameters::BulkShear { .. } => {}
            }
        }
        Ok(())
    }
    fn bound_curvature(&self) -> Result<f64, ThermoelasticError> {
        let coefficients: Vec<_> = self.elastic.iter().map(Coefficient::bound).collect();
        let (q, g) = plane(&self.source.input().elasticity, &coefficients);
        if q.iter()
            .chain(g.iter())
            .any(|j| ![j.v, j.d, j.dd].iter().all(|v| v.finite()))
        {
            return Err(invalid(
                "whole interval cannot bound a positive plane-stress denominator",
            ));
        }
        let alpha: [f64; 5] = core::array::from_fn(|i| {
            if i >= 2 {
                0.0
            } else {
                let c = &self.expansion[if self.expansion.len() == 1 { 0 } else { i }];
                c.bound().v.abs_max()
            }
        });
        let alpha_first: [f64; 5] = core::array::from_fn(|i| {
            if i >= 2 {
                0.0
            } else {
                self.expansion[if self.expansion.len() == 1 { 0 } else { i }]
                    .slope
                    .abs()
            }
        });
        let mut bound = 0.0;
        for i in 0..5 {
            for j in 0..5 {
                let c = if i < 3 && j < 3 {
                    q[3 * i + j]
                } else if i == j {
                    g[i - 3]
                } else {
                    Jet::constant(0.0)
                };
                let curvature = upper_mul(
                    0.5,
                    upper_mul(c.dd.abs_max(), upper_mul(self.limit, self.limit)?)?,
                )?;
                let mixed = upper_mul(
                    2.0,
                    upper_mul(alpha[i], upper_mul(c.d.abs_max(), self.limit)?)?,
                )?;
                let expansion = upper_mul(alpha[i], upper_mul(c.v.abs_max(), alpha[j])?)?;
                let acceleration =
                    upper_mul(alpha_first[i], upper_mul(c.v.abs_max(), self.limit)?)?;
                bound = upper_add(
                    bound,
                    upper_add(
                        upper_add(curvature, mixed)?,
                        upper_add(expansion, acceleration)?,
                    )?,
                )?;
            }
        }
        // Retain a rounding allowance for equivalent native Q arithmetic and
        // source-affinity admission, in addition to outward interval operations.
        upper_mul(bound, 1.0 + 4096.0 * f64::EPSILON)
    }
}

/// Actual source state and exact temperature derivatives in material axes.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLaminaThermoelastic {
    elastic: ResolvedMaterialStatePoint,
    expansion: IntegratedLaminaExpansion,
    temperature: f64,
    limit: f64,
    q: [Jet<f64>; 9],
    g: [Jet<f64>; 2],
    free: [f64; 2],
    alpha: [Jet<f64>; 2],
}
impl ResolvedLaminaThermoelastic {
    /// Original actual elastic receipts, without refreshed reference density.
    pub const fn elasticity(&self) -> &ResolvedMaterialStatePoint {
        &self.elastic
    }
    /// Original full integrated CTE path, including its reference receipts.
    pub const fn expansion(&self) -> &IntegratedLaminaExpansion {
        &self.expansion
    }
    /// Actual held temperature [K].
    pub const fn temperature_k(&self) -> f64 {
        self.temperature
    }
    /// Row-major principal engineering matrices `[Q,Q_T,Q_TT]`.
    pub fn stiffness_derivatives_pa(&self) -> [[f64; 9]; 3] {
        [
            self.q.map(|j| j.v),
            self.q.map(|j| j.d),
            self.q.map(|j| j.dd),
        ]
    }
    /// Uncorrected principal transverse `[G13,G23]` and its first/second T derivatives.
    pub fn transverse_shear_derivatives_pa(&self) -> [[f64; 2]; 3] {
        [
            self.g.map(|j| j.v),
            self.g.map(|j| j.d),
            self.g.map(|j| j.dd),
        ]
    }
    /// Integrated principal free engineering strain, with no in-plane shear.
    pub fn free_strain(&self) -> [f64; 3] {
        [self.free[0], self.free[1], 0.0]
    }
    /// In-plane-only potential for native laminate thickness integration.
    pub fn in_plane_state(
        &self,
        total: [f64; 3],
        reference_volume_m3: f64,
    ) -> Result<LaminaThermoelasticState, ThermoelasticError> {
        self.state(
            [total[0], total[1], total[2], 0.0, 0.0],
            reference_volume_m3,
        )
    }
    /// Evaluate the local potential at fixed total principal engineering strain.
    /// Last two entries use raw G13/G23, WITHOUT the laminate Whitney factor.
    /// Evaluation allocates nothing and performs no source lookup or rotation.
    pub fn state(
        &self,
        total: [f64; 5],
        reference_volume_m3: f64,
    ) -> Result<LaminaThermoelasticState, ThermoelasticError> {
        let v = reference_volume_m3;
        if !v.is_finite()
            || v <= 0.0
            || total.iter().any(|e| !e.is_finite() || e.abs() > self.limit)
        {
            return Err(invalid(
                "total strain or reference volume leaves the admitted layer envelope",
            ));
        }
        let mut e = total;
        for i in 0..2 {
            e[i] -= self.free[i];
        }
        if e.iter().any(|e| !e.is_finite() || e.abs() > self.limit) {
            return Err(invalid("elastic strain leaves the admitted layer envelope"));
        }
        let alpha = [self.alpha[0].v, self.alpha[1].v, 0.0, 0.0, 0.0];
        let alpha_first = [self.alpha[0].d, self.alpha[1].d, 0.0, 0.0, 0.0];
        let apply = |x: [f64; 5], order: usize| -> [f64; 5] {
            let get = |j: Jet<f64>| [j.v, j.d, j.dd][order];
            core::array::from_fn(|i| {
                if i < 3 {
                    (0..3).map(|j| get(self.q[3 * i + j]) * x[j]).sum()
                } else {
                    get(self.g[i - 3]) * x[i]
                }
            })
        };
        let stress = apply(e, 0);
        let first = apply(e, 1);
        let thermal = apply(alpha, 0);
        let mixed = core::array::from_fn(|i| first[i] - thermal[i]);
        let psi = finite(0.5 * v * dot(e, stress))?;
        let entropy = finite(v * (dot(alpha, stress) - 0.5 * dot(e, first)))?;
        let internal = finite(self.temperature.mul_add(entropy, psi))?;
        let curvature = 0.5 * dot(e, apply(e, 2)) - 2.0 * dot(alpha, first) + dot(alpha, thermal)
            - dot(alpha_first, stress);
        let capacity = finite(-self.temperature * v * curvature)?;
        if psi < 0.0 || stress.iter().chain(mixed.iter()).any(|v| !v.is_finite()) {
            return Err(invalid("layer thermoelastic response is unrepresentable"));
        }
        Ok(LaminaThermoelasticState {
            elastic_strain_engineering: e,
            stress_pa: stress,
            stress_temperature_derivative_pa_k: mixed,
            free_energy_j: psi,
            entropy_j_k: entropy,
            internal_energy_j: internal,
            heat_capacity_j_k: capacity,
        })
    }
}

/// Elastic thermodynamic contributions, separate from the original caloric law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminaThermoelasticState {
    pub elastic_strain_engineering: [f64; 5],
    pub stress_pa: [f64; 5],
    /// Work-conjugate d(stress)/dT at fixed TOTAL engineering strain [Pa/K].
    pub stress_temperature_derivative_pa_k: [f64; 5],
    pub free_energy_j: f64,
    pub entropy_j_k: f64,
    pub internal_energy_j: f64,
    /// Signed fixed-strain elastic contribution, not total heat capacity.
    pub heat_capacity_j_k: f64,
}
fn dot(a: [f64; 5], b: [f64; 5]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

// The same rational expressions evaluate exact second derivatives and their
// outward interval enclosures. No finite-difference source queries enter.
trait Number:
    Copy
    + From<f64>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
}
impl Number for f64 {}
impl Number for Interval {}
#[derive(Clone, Copy, Debug, PartialEq)]
struct Jet<T> {
    v: T,
    d: T,
    dd: T,
}
impl<T: Number> Jet<T> {
    fn constant(v: f64) -> Self {
        Self {
            v: T::from(v),
            d: T::from(0.0),
            dd: T::from(0.0),
        }
    }
}
impl<T: Number> Add for Jet<T> {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        Self {
            v: self.v + b.v,
            d: self.d + b.d,
            dd: self.dd + b.dd,
        }
    }
}
impl<T: Number> Neg for Jet<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: -self.d,
            dd: -self.dd,
        }
    }
}
impl<T: Number> Sub for Jet<T> {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl<T: Number> Mul for Jet<T> {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            d: self.d * b.v + self.v * b.d,
            dd: self.dd * b.v + T::from(2.0) * self.d * b.d + self.v * b.dd,
        }
    }
}
impl<T: Number> Div for Jet<T> {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        let inverse = T::from(1.0) / b.v;
        let relative = b.d * inverse;
        self * Self {
            v: inverse,
            d: -relative * inverse,
            dd: (T::from(2.0) * relative * relative - b.dd * inverse) * inverse,
        }
    }
}
fn plane<T: Number>(schema: &LaminaElasticSource, c: &[Jet<T>]) -> ([Jet<T>; 9], [Jet<T>; 2]) {
    let one = Jet::constant(1.0);
    let two = Jet::constant(2.0);
    let three = Jet::constant(3.0);
    let four = Jet::constant(4.0);
    let (q11, q22, q12, g12, transverse) = match schema {
        LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson { .. }) => {
            let q = c[0] / (one - c[1] * c[1]);
            let g = c[0] / (two * (one + c[1]));
            (q, q, c[1] * q, g, [g; 2])
        }
        LaminaElasticSource::Isotropic(IsotropicElasticParameters::BulkShear { .. }) => {
            let denominator = three * c[0] + four * c[1];
            let q = four * c[1] * ((three * c[0] + c[1]) / denominator);
            let cross = two * c[1] * ((three * c[0] - two * c[1]) / denominator);
            (q, q, cross, c[1], [c[1]; 2])
        }
        LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungShear { .. }) => {
            let denominator = four * c[1] - c[0];
            let q = four * c[1] * (c[1] / denominator);
            let cross = two * c[1] * ((c[0] - two * c[1]) / denominator);
            (q, q, cross, c[1], [c[1]; 2])
        }
        LaminaElasticSource::OrthotropicPlane12 => {
            let denominator = one - c[2] * c[2] * (c[1] / c[0]);
            (
                c[0] / denominator,
                c[1] / denominator,
                c[2] * (c[1] / denominator),
                c[3],
                [c[4], c[5]],
            )
        }
    };
    let z = Jet::constant(0.0);
    ([q11, q12, z, q12, q22, z, z, z, g12], transverse)
}
fn native_values(
    schema: &LaminaElasticSource,
    state: &ResolvedMaterialStatePoint,
    c: &[Jet<f64>],
) -> Result<([f64; 9], [f64; 2]), ThermoelasticError> {
    let (e, nu, g, transverse) = match schema {
        LaminaElasticSource::Isotropic(p) => {
            let law = ResolvedIsotropicElasticModuli::try_from_resolved(state.clone(), p)?;
            (
                [law.young_modulus_pa(); 2],
                law.poisson_ratio(),
                law.shear_modulus_pa(),
                [law.shear_modulus_pa(); 2],
            )
        }
        LaminaElasticSource::OrthotropicPlane12 => {
            ([c[0].v, c[1].v], c[2].v, c[3].v, [c[4].v, c[5].v])
        }
    };
    let denominator = 1.0 - nu * (nu * e[1] / e[0]);
    let cross = nu * e[1] / denominator;
    let q = [
        e[0] / denominator,
        cross,
        0.0,
        cross,
        e[1] / denominator,
        0.0,
        0.0,
        0.0,
        g,
    ];
    crate::sheet::SheetMatrix::from_coefficients(q)
        .map_err(|_| invalid("native plane-stress stiffness is inadmissible"))?;
    if transverse.iter().any(|v| !v.is_finite() || *v <= 0.0) {
        return Err(invalid("positive native transverse shear required"));
    }
    Ok((q, transverse))
}

// Nonfinite interval arithmetic remains an explicit failed bound. Never use
// min/max to discard a NaN corner or invert an interval crossing a pole.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Interval {
    lo: f64,
    hi: f64,
}
impl Interval {
    const WHOLE: Self = Self {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };
    fn range(lo: f64, hi: f64) -> Self {
        if !lo.is_finite() || !hi.is_finite() {
            Self::WHOLE
        } else if lo == 0.0 && hi == 0.0 {
            Self::from(0.0)
        } else {
            Self {
                lo: lo.next_down(),
                hi: hi.next_up(),
            }
        }
    }
    fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite() && self.lo <= self.hi
    }
    fn zero(self) -> bool {
        self.lo == 0.0 && self.hi == 0.0
    }
    fn abs_max(self) -> f64 {
        self.lo.abs().max(self.hi.abs())
    }
}
impl From<f64> for Interval {
    fn from(x: f64) -> Self {
        Self { lo: x, hi: x }
    }
}
impl Add for Interval {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        if self.zero() {
            b
        } else if b.zero() {
            self
        } else {
            Self::range(self.lo + b.lo, self.hi + b.hi)
        }
    }
}
impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
}
impl Sub for Interval {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + (-b)
    }
}
impl Mul for Interval {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        if !self.finite() || !b.finite() {
            return Self::WHOLE;
        }
        if self.zero() || b.zero() {
            return Self::from(0.0);
        }
        let products = [
            self.lo * b.lo,
            self.lo * b.hi,
            self.hi * b.lo,
            self.hi * b.hi,
        ];
        if products.iter().any(|x| !x.is_finite()) {
            return Self::WHOLE;
        }
        Self::range(
            products.iter().copied().fold(f64::INFINITY, f64::min),
            products.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }
}
impl Div for Interval {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        if !b.finite() || b.lo <= 0.0 {
            Self::WHOLE
        } else {
            self * Self::range(1.0 / b.hi, 1.0 / b.lo)
        }
    }
}
fn upper_add(a: f64, b: f64) -> Result<f64, ThermoelasticError> {
    if a == 0.0 && b == 0.0 {
        Ok(0.0)
    } else {
        finite((a + b).next_up())
    }
}
fn upper_mul(a: f64, b: f64) -> Result<f64, ThermoelasticError> {
    if a == 0.0 || b == 0.0 {
        Ok(0.0)
    } else {
        finite((a * b).next_up())
    }
}
