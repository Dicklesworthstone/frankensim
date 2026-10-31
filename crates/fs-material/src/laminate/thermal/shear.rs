//! Temperature derivatives of the ORIGINAL whole-stack shear approximation.
//!
//! The native laminate uses directional Whitney neutral axes and energy
//! integrals, followed by a positive diagonal congruence for off-axis shear.
//! Every layer affects those integrals. Differentiating raw local G, or putting
//! a separate 5/6 on each layer, therefore gives the wrong thermal potential.
//! This module differentiates that same native expression with fs-ad duals.
//!
//! For one section shear strain gamma and reference area A, Psi = A gamma^T H
//! gamma / 2. Its layer-temperature conjugates are S_i = -Psi_{T_i}; the one
//! stack internal energy is U = Psi + sum_i T_i S_i. There is deliberately no
//! arbitrary allocation of shear free/internal energy to individual layers.
//! The entropy Hessian action exposes the nonlocal thermal-capacity response.
//!
//! This retains the native directional Whitney approximation; it does not
//! promote it to a general coupled three-dimensional warping solution. The
//! native numerical 5/6 fallback is refused here, even for a zero derivative
//! direction. No derivative is claimed through that fallback or a source kink.
//! Geometry, orientation, reference density and mass stay fixed.

use super::source::ResolvedLaminateThermal;
use super::thermoelastic::AffineLaminaThermoelastic;
use crate::elastic::thermal::ThermoelasticError;
use crate::laminate::{rotate_reduced, rotate_transverse};
use fs_ad::{Dual, Real};
use std::sync::Arc;

/// Refusal before exposing a complete source-pinned derivative or potential.
#[derive(Debug)]
pub enum LaminateShearThermoelasticError {
    /// One original layer source failed its complete potential admission.
    Source {
        /// Original spatial layer index.
        layer: usize,
        /// Original source refusal.
        source: ThermoelasticError,
    },
    /// Shape, source identity, strain, numerical branch or arithmetic refusal.
    Invalid(&'static str),
    /// Cancellation occurred before a complete result was available.
    Cancelled,
}
impl core::fmt::Display for LaminateShearThermoelasticError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MATERIAL-LAMINATE-SHEAR-THERMOELASTIC: {self:?}")
    }
}
impl std::error::Error for LaminateShearThermoelasticError {}

fn invalid(message: &'static str) -> LaminateShearThermoelasticError {
    LaminateShearThermoelasticError::Invalid(message)
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), LaminateShearThermoelasticError> {
    if cancelled() {
        Err(LaminateShearThermoelasticError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(x: f64) -> Result<f64, LaminateShearThermoelasticError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(invalid("whole-stack shear arithmetic is unrepresentable"))
    }
}
fn positive(x: f64) -> Result<(), LaminateShearThermoelasticError> {
    if x.is_finite() && x > 0.0 {
        Ok(())
    } else {
        Err(invalid(
            "native Whitney expression needs a numerical fallback; no smooth derivative is admitted",
        ))
    }
}

#[derive(Clone, Debug)]
struct Layer {
    bottom: f64,
    top: f64,
    /// Rotated Q and its exact two temperature derivatives.
    q: [[f64; 9]; 3],
    /// Rotated uncorrected G and its exact two temperature derivatives.
    g: [[f64; 4]; 3],
    sin: f64,
    cos: f64,
    strain_limit: f64,
}

/// Source-pinned current section with temperature gradients cached once.
///
/// Preparation is O(n²) scalar work and O(n) storage for n original layers.
/// A directional derivative is O(n), and a Hessian-vector action is O(n²)
/// without allocating a dense n-by-n Hessian. Point-state evaluation uses only
/// the cached native operator and gradients, with no source query or AD.
#[derive(Clone, Debug)]
pub struct ResolvedLaminateShearThermoelastic {
    current: ResolvedLaminateThermal,
    layers: Vec<Layer>,
    temperatures: Vec<f64>,
    stiffness: [f64; 4],
    gradient: Vec<[f64; 4]>,
}

/// One whole-stack elastic shear potential, integrated over reference area.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateShearThermoelasticState {
    /// Section shear force resultants H gamma [N/m], in section xz/yz axes.
    pub force_resultant_n_m: [f64; 2],
    /// Whole-stack Helmholtz free energy [J].
    pub free_energy_j: f64,
    /// Original-cell entropy conjugates -dPsi/dT_i [J/K], in spatial order.
    pub entropy_j_k: Vec<f64>,
    /// Whole-stack U = Psi + sum_i T_i S_i [J].
    pub internal_energy_j: f64,
}

/// Fixed-total-shear thermal response along T_i(s) = T_i + direction_i s.
/// The path parameter s is in kelvin and `direction` is dimensionless.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateShearThermalResponse {
    /// dS_i/ds [J/K²]. This is minus the free-energy Hessian times direction.
    pub entropy_derivative_j_k2: Vec<f64>,
    /// T_i dS_i/ds [J/K] for each thermal coordinate. These contributions can
    /// be negative and must be combined with the original caloric operator.
    pub capacity_action_j_k: Vec<f64>,
    /// dU/ds = sum_i T_i dS_i/ds [J/K] at fixed total shear strain.
    pub internal_energy_derivative_j_k: f64,
}

impl ResolvedLaminateShearThermoelastic {
    /// Adopt exactly the current native source stack. Each law must belong to
    /// its original card, geometry and reference receipts, in spatial order.
    /// No new material query or source reselection occurs here.
    pub fn try_new(
        laws: &[AffineLaminaThermoelastic],
        current: &ResolvedLaminateThermal,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, LaminateShearThermoelasticError> {
        check(&mut cancelled)?;
        let original = current.source().layers();
        if laws.len() != original.len() || current.layers().len() != original.len() {
            return Err(invalid(
                "one exact potential per original source layer is required",
            ));
        }
        let native = current.prepared().section();
        let mut layers = Vec::with_capacity(laws.len());
        let mut temperatures = Vec::with_capacity(laws.len());
        for (i, ((law, source), material)) in
            laws.iter().zip(original).zip(current.layers()).enumerate()
        {
            check(&mut cancelled)?;
            let input = source.input();
            let bound = law.source();
            if !Arc::ptr_eq(&bound.input().card, &input.card)
                || bound.input().thickness_m.to_bits() != input.thickness_m.to_bits()
                || bound.input().orientation_rad.to_bits() != input.orientation_rad.to_bits()
                || bound.reference_material() != source.reference_material()
                || bound.reference_expansion().reference()
                    != source.reference_expansion().reference()
            {
                return Err(invalid(
                    "shear potential must retain each exact original source and layer geometry",
                ));
            }
            let resolved = law
                .at_materials(material)
                .map_err(|source| LaminateShearThermoelasticError::Source { layer: i, source })?;
            let q_principal = resolved.stiffness_derivatives_pa();
            let g_principal = resolved.transverse_shear_derivatives_pa();
            let mut q = [[0.0; 9]; 3];
            let mut g = [[0.0; 4]; 3];
            for derivative in 0..3 {
                q[derivative] = rotate_reduced(&q_principal[derivative], input.orientation_rad)
                    .ok_or(invalid(
                        "rotated layer stiffness derivative is unrepresentable",
                    ))?;
                g[derivative] = rotate_transverse(g_principal[derivative], input.orientation_rad);
                if g[derivative].iter().any(|v| !v.is_finite()) {
                    return Err(invalid("rotated shear derivative is unrepresentable"));
                }
            }
            let geometry = &native.layer_integrals()[i];
            let (sin, cos) = input.orientation_rad.sin_cos();
            layers.push(Layer {
                bottom: geometry.z_bottom_m,
                top: geometry.z_top_m,
                q,
                g,
                sin,
                cos,
                strain_limit: law.maximum_absolute_engineering_strain(),
            });
            temperatures.push(resolved.temperature_k());
        }
        // Evaluate the unmodified smooth Whitney branch before differentiating.
        // A constant direction must not hide the native 5/6 fallback.
        let value = evaluate(&layers, &|_, c| c[0], &mut cancelled)?;
        let native_value = *native.shear_pa_m();
        for (&actual, &expected) in value.iter().zip(&native_value) {
            let scale = actual.abs().max(expected.abs());
            if (actual - expected).abs() > 4096.0 * f64::EPSILON * scale {
                return Err(invalid(
                    "source shear potential disagrees with the native section operator",
                ));
            }
        }
        let mut gradient = Vec::with_capacity(laws.len());
        for coordinate in 0..laws.len() {
            check(&mut cancelled)?;
            let derivative = evaluate(
                &layers,
                &|layer, c| Dual::<f64, 1> {
                    re: c[0],
                    eps: [if layer == coordinate { c[1] } else { 0.0 }],
                },
                &mut cancelled,
            )?;
            let row = derivative.map(|d| d.eps[0]);
            for &x in &row {
                finite(x)?;
            }
            gradient.push(row);
        }
        check(&mut cancelled)?;
        Ok(Self {
            current: current.clone(),
            layers,
            temperatures,
            stiffness: native_value,
            gradient,
        })
    }

    /// Native current section, integrated expansion and original source receipts.
    pub const fn current(&self) -> &ResolvedLaminateThermal {
        &self.current
    }
    /// Actual original layer temperatures [K], in spatial order.
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures
    }
    /// Original native H [Pa m], row-major 2-by-2 in section xz/yz axes.
    pub const fn stiffness_pa_m(&self) -> &[f64; 4] {
        &self.stiffness
    }
    /// dH/dT_i [Pa m/K] for each original layer, including all neutral-axis,
    /// Whitney energy-integral and cross-shear dependencies.
    pub fn stiffness_temperature_derivatives_pa_m_k(&self) -> &[[f64; 4]] {
        &self.gradient
    }

    /// `[H, dH/ds, d²H/ds²]` along the straight layer-temperature direction.
    /// This is a local derivative on the admitted source charts, not permission
    /// to step beyond their intervals or a whole-interval capacity enclosure.
    pub fn directional_derivatives(
        &self,
        direction: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[[f64; 4]; 3], LaminateShearThermoelasticError> {
        check(&mut cancelled)?;
        self.admit_direction(direction)?;
        let values = evaluate(
            &self.layers,
            &|layer, c| seed(c, direction[layer], direction[layer]),
            &mut cancelled,
        )?;
        let answer = [
            self.stiffness,
            values.map(|d| d.re.eps[0]),
            values.map(|d| d.eps[0].eps[0]),
        ];
        for row in &answer {
            for &v in row {
                finite(v)?;
            }
        }
        check(&mut cancelled)?;
        Ok(answer)
    }

    /// Evaluate the one stack potential and its per-cell entropy conjugates.
    /// Both principal transverse engineering-shear components must remain in
    /// every layer's original declared strain envelope.
    pub fn state(
        &self,
        shear_strain: [f64; 2],
        reference_area_m2: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShearThermoelasticState, LaminateShearThermoelasticError> {
        check(&mut cancelled)?;
        self.admit_state(shear_strain, reference_area_m2, &mut cancelled)?;
        let free_energy_j = quadratic(self.stiffness, shear_strain, reference_area_m2)?;
        let force_resultant_n_m = [
            finite(self.stiffness[0] * shear_strain[0] + self.stiffness[1] * shear_strain[1])?,
            finite(self.stiffness[2] * shear_strain[0] + self.stiffness[3] * shear_strain[1])?,
        ];
        let mut entropy_j_k = Vec::with_capacity(self.layers.len());
        let mut internal_energy_j = free_energy_j;
        for (&gradient, &temperature) in self.gradient.iter().zip(&self.temperatures) {
            check(&mut cancelled)?;
            let entropy = -quadratic(gradient, shear_strain, reference_area_m2)?;
            entropy_j_k.push(entropy);
            internal_energy_j = finite(internal_energy_j + temperature * entropy)?;
        }
        check(&mut cancelled)?;
        Ok(LaminateShearThermoelasticState {
            force_resultant_n_m,
            free_energy_j,
            entropy_j_k,
            internal_energy_j,
        })
    }

    /// The complete entropy Hessian action at fixed total shear. All mixed
    /// layer-temperature derivatives are retained; no dense Hessian is stored.
    /// A consumer must add its original caloric and in-plane contributions and
    /// admit the resulting coupled thermal operator before advancing heat.
    pub fn thermal_response(
        &self,
        shear_strain: [f64; 2],
        reference_area_m2: f64,
        direction: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShearThermalResponse, LaminateShearThermoelasticError> {
        check(&mut cancelled)?;
        self.admit_direction(direction)?;
        self.admit_state(shear_strain, reference_area_m2, &mut cancelled)?;
        let mut entropy_derivative_j_k2 = Vec::with_capacity(self.layers.len());
        let mut capacity_action_j_k = Vec::with_capacity(self.layers.len());
        let mut internal_energy_derivative_j_k = 0.0;
        for coordinate in 0..self.layers.len() {
            check(&mut cancelled)?;
            let response = evaluate(
                &self.layers,
                &|layer, c| {
                    seed(
                        c,
                        if layer == coordinate { 1.0 } else { 0.0 },
                        direction[layer],
                    )
                },
                &mut cancelled,
            )?;
            let mixed = response.map(|d| d.eps[0].eps[0]);
            let entropy = -quadratic(mixed, shear_strain, reference_area_m2)?;
            let capacity = finite(self.temperatures[coordinate] * entropy)?;
            entropy_derivative_j_k2.push(entropy);
            capacity_action_j_k.push(capacity);
            internal_energy_derivative_j_k = finite(internal_energy_derivative_j_k + capacity)?;
        }
        check(&mut cancelled)?;
        Ok(LaminateShearThermalResponse {
            entropy_derivative_j_k2,
            capacity_action_j_k,
            internal_energy_derivative_j_k,
        })
    }

    fn admit_direction(&self, direction: &[f64]) -> Result<(), LaminateShearThermoelasticError> {
        if direction.len() != self.layers.len() || direction.iter().any(|v| !v.is_finite()) {
            Err(invalid(
                "one finite temperature direction per original layer is required",
            ))
        } else {
            Ok(())
        }
    }
    fn admit_state(
        &self,
        gamma: [f64; 2],
        area: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), LaminateShearThermoelasticError> {
        if !area.is_finite() || area <= 0.0 || gamma.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "positive reference area and finite total shear strain are required",
            ));
        }
        for layer in &self.layers {
            check(cancelled)?;
            let principal = [
                layer.cos * gamma[0] + layer.sin * gamma[1],
                -layer.sin * gamma[0] + layer.cos * gamma[1],
            ];
            if principal
                .iter()
                .any(|v| !v.is_finite() || v.abs() > layer.strain_limit)
            {
                return Err(invalid(
                    "total principal shear leaves an original layer strain envelope",
                ));
            }
        }
        Ok(())
    }
}

type Mixed = Dual<Dual<f64, 1>, 1>;
fn seed(c: [f64; 3], a: f64, b: f64) -> Mixed {
    Dual {
        re: Dual {
            re: c[0],
            eps: [c[1] * a],
        },
        eps: [Dual {
            re: c[1] * b,
            eps: [c[2] * a * b],
        }],
    }
}
fn quadratic(
    h: [f64; 4],
    gamma: [f64; 2],
    area: f64,
) -> Result<f64, LaminateShearThermoelasticError> {
    finite(
        0.5 * area
            * (gamma[0] * (h[0] * gamma[0] + h[1] * gamma[1])
                + gamma[1] * (h[2] * gamma[0] + h[3] * gamma[1])),
    )
}

// Same neutral-axis and polynomial energy integrals as LaminateSection. The
// only different branch policy is a typed refusal in place of its numerical
// fallback. fs-ad supplies the derivative arithmetic; no local AD engine or
// finite-difference approximation is introduced here.
fn evaluate<T: Real>(
    layers: &[Layer],
    seed: &impl Fn(usize, [f64; 3]) -> T,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[T; 4], LaminateShearThermoelasticError> {
    let mut shear = [T::zero(); 4];
    for (i, layer) in layers.iter().enumerate() {
        check(cancelled)?;
        for (k, value) in shear.iter_mut().enumerate() {
            *value =
                *value + seed(i, layer.g.map(|row| row[k])) * T::from_f64(layer.top - layer.bottom);
        }
    }
    let x = factor(layers, 0, 0, shear[0], seed, cancelled)?;
    let y = factor(layers, 4, 3, shear[3], seed, cancelled)?;
    shear[0] = shear[0] * x;
    shear[3] = shear[3] * y;
    shear[1] = shear[1] * (x.sqrt() * y.sqrt());
    shear[2] = shear[1];
    for value in shear {
        finite(value.value())?;
    }
    Ok(shear)
}

fn factor<T: Real>(
    layers: &[Layer],
    in_plane: usize,
    transverse: usize,
    shear: T,
    seed: &impl Fn(usize, [f64; 3]) -> T,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<T, LaminateShearThermoelasticError> {
    positive(shear.value())?;
    let mut a = T::zero();
    let mut first = T::zero();
    for (i, layer) in layers.iter().enumerate() {
        check(cancelled)?;
        let q = seed(i, layer.q.map(|row| row[in_plane]));
        positive(q.value())?;
        let h = T::from_f64(layer.top - layer.bottom);
        a = a + q * h;
        first = first + q * T::from_f64(0.5) * h * T::from_f64(layer.top + layer.bottom);
    }
    positive(a.value())?;
    let neutral = first / a;
    let mut bending = T::zero();
    let mut energy = T::zero();
    let mut running = T::zero();
    for (i, layer) in layers.iter().enumerate() {
        check(cancelled)?;
        let q = seed(i, layer.q.map(|row| row[in_plane]));
        let g = seed(i, layer.g.map(|row| row[transverse]));
        positive(g.value())?;
        let u0 = T::from_f64(layer.bottom) - neutral;
        let u1 = T::from_f64(layer.top) - neutral;
        bending = bending + q * (u1 * u1 * u1 - u0 * u0 * u0) / T::from_f64(3.0);
        let half_q = T::from_f64(0.5) * q;
        let c0 = running - half_q * u0 * u0;
        let int_1 = u1 - u0;
        let int_u2 = (u1.powi(3) - u0.powi(3)) / T::from_f64(3.0);
        let int_u4 = (u1.powi(5) - u0.powi(5)) / T::from_f64(5.0);
        let r2 =
            c0 * c0 * int_1 + T::from_f64(2.0) * c0 * half_q * int_u2 + half_q * half_q * int_u4;
        energy = energy + r2 / g;
        running = running + half_q * (u1 * u1 - u0 * u0);
    }
    positive(bending.value())?;
    positive(energy.value())?;
    let correction = bending * bending / (shear * energy);
    positive(correction.value())?;
    Ok(correction)
}
