//! Authored Landau–Khalatnikov longitudinal dissipation at one frozen state.
//!
//! Define `q0 = rho_s (v_s-v_n)`. The longitudinal viscous stress is
//! `(4 eta/3 + zeta2) div(v_n) + zeta1 div(q0)` and the dissipative chemical
//! potential is `-zeta1 div(v_n) - zeta3 div(q0)`. Onsager reciprocity fixes
//! `zeta4=zeta1`; nonnegative entropy production requires `zeta2,zeta3>=0`
//! and `zeta1²<=zeta2*zeta3`. Shear viscosity cannot rescue an invalid bulk
//! matrix. The cross coefficient may have either sign.
//!
//! References: Nikuni & Griffin, Phys. Rev. A 63, 033608 (2001), equations
//! (86)–(87), <https://arxiv.org/abs/cond-mat/0009333>; Escobedo, Mannarelli &
//! Manuel, Phys. Rev. A 79, 063623 (2009), equations (8)–(15),
//! <https://arxiv.org/abs/0904.3023>. Only their constitutive form is used;
//! no transport coefficients or helium calibration are inferred.
//!
//! Viscous kinetic loss becomes heat. Fourier conduction instead redistributes
//! internal energy and destroys wave availability: `kappa/T0 * |grad(T')|²`
//! is NOT an additional caloric heat source. This local closure has no spatial
//! boundary conditions, vortex drag, transverse wall drag, thermal feedback,
//! finite-frequency relaxation, or phase-transition model.

use super::{TwoFluidState, invalid, nonnegative};
use crate::MaterialError;

/// Immutable reciprocal second viscosities on the TOTAL mixture-volume basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SecondViscosities {
    zeta1_m2_s: f64,
    zeta2_pa_s: f64,
    zeta3_m5_kg_s: f64,
    // Bulk loss = (factor[0]*div(v_n)+factor[1]*div(q0))²
    //           + (factor[2]*div(q0))².
    factor: [f64; 3],
}

impl SecondViscosities {
    /// Explicit absence of bulk dissipation, not a material-name default.
    pub const ZERO: Self = Self {
        zeta1_m2_s: 0.0,
        zeta2_pa_s: 0.0,
        zeta3_m5_kg_s: 0.0,
        factor: [0.0; 3],
    };

    /// Admit signed `zeta1 [m²/s]`, `zeta2 [Pa s]`, `zeta3 [m⁵/(kg s)]`.
    /// The positive-semidefinite inequality is checked on exact products of
    /// the supplied binary64 values, without overflow or acceptance tolerance.
    pub fn try_new(
        zeta1_m2_s: f64,
        zeta2_pa_s: f64,
        zeta3_m5_kg_s: f64,
    ) -> Result<Self, MaterialError> {
        if !zeta1_m2_s.is_finite() || !nonnegative(zeta2_pa_s) || !nonnegative(zeta3_m5_kg_s) {
            return Err(invalid(
                "second viscosities need finite signed zeta1 and finite nonnegative zeta2/zeta3",
            ));
        }
        let cross = exact_product(zeta1_m2_s.abs(), zeta1_m2_s.abs());
        let diagonal = exact_product(zeta2_pa_s, zeta3_m5_kg_s);
        if cross > diagonal {
            return Err(invalid(
                "second viscosities require zeta1 squared <= zeta2 times zeta3",
            ));
        }
        let a = zeta2_pa_s.sqrt();
        let b = if a == 0.0 { 0.0 } else { zeta1_m2_s / a };
        // Form the Schur complement from the exact product difference. A
        // rounded 1-correlation² can erase a small positive dissipation mode.
        let fraction = product_difference_fraction(diagonal, cross);
        let c = zeta3_m5_kg_s.sqrt() * fraction.sqrt();
        Ok(Self {
            zeta1_m2_s,
            zeta2_pa_s,
            zeta3_m5_kg_s,
            factor: [a, b, c],
        })
    }

    pub const fn zeta1_m2_s(&self) -> f64 {
        self.zeta1_m2_s
    }
    pub const fn zeta2_pa_s(&self) -> f64 {
        self.zeta2_pa_s
    }
    pub const fn zeta3_m5_kg_s(&self) -> f64 {
        self.zeta3_m5_kg_s
    }
}

/// Constitutive fluxes and distinct local accounts for supplied SI gradients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongitudinalResponse {
    /// Positive tensile viscous Cauchy stress [Pa], excluding equilibrium pressure.
    pub longitudinal_viscous_stress_pa: f64,
    /// Additive dissipative specific chemical potential [J/kg].
    pub chemical_potential_correction_j_kg: f64,
    /// Fourier heat flux `-kappa grad(T')` [W/m²]; excludes counterflow heat.
    pub conductive_heat_flux_w_m2: f64,
    /// Fourier entropy flux at the frozen reference temperature [W/(m² K)].
    pub conductive_entropy_flux_w_m2_k: f64,
    /// Kinetic energy converted into caloric heat [W/m³].
    pub viscous_heat_w_m3: f64,
    /// Wave-availability destruction by conduction [W/m³], NOT extra heat.
    pub conductive_availability_loss_w_m3: f64,
    /// Sum of both irreversible entropy-production rates [W/(m³ K)].
    pub entropy_production_w_m3_k: f64,
}

/// One original material and explicit bulk viscosities/thermal conductivity.
/// Preparation caches the uniform operators; evaluation allocates nothing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongitudinalTransport {
    material: TwoFluidState,
    second: SecondViscosities,
    conductivity: f64,
    shear_root: f64,
    longitudinal_viscosity: f64,
    conductive_root: f64,
    resistance: [[f64; 2]; 2],
    kinetic_diffusivity: [[f64; 2]; 2],
    thermal_diffusivity: f64,
}

impl LongitudinalTransport {
    /// Normal shear viscosity comes from `material`; it is never entered twice.
    /// Coefficients are real, frequency independent and frozen at this state.
    pub fn try_new(
        material: TwoFluidState,
        second_viscosities: SecondViscosities,
        thermal_conductivity_w_m_k: f64,
    ) -> Result<Self, MaterialError> {
        if !nonnegative(thermal_conductivity_w_m_k) {
            return Err(invalid(
                "superfluid Fourier conductivity must be finite and nonnegative",
            ));
        }
        let shear = (4.0 / 3.0) * material.normal_viscosity_pa_s();
        let shear_root = shear.sqrt();
        let longitudinal_viscosity = shear + second_viscosities.zeta2_pa_s;
        let conductive_root = thermal_conductivity_w_m_k.sqrt() / material.temperature_k().sqrt();
        let fs = material.superfluid_fraction();
        let rho_s = material.superfluid_density_kg_m3();
        let [a, b, c] = second_viscosities.factor;
        // q0_x=-rho_s*w_x, v_n,x=v_x+f_s*w_x. A Gram representation
        // preserves the semidefinite closure, including its null directions.
        let factors = [
            [shear_root, shear_root * fs],
            [a, a * fs - b * rho_s],
            [0.0, -c * rho_s],
        ];
        let resistance = gram(factors)?;
        let rho_root = material.density_kg_m3().sqrt();
        let mu_root = material.reduced_density_kg_m3().sqrt();
        let normalized = factors.map(|row| [row[0] / rho_root, row[1] / mu_root]);
        let kinetic_diffusivity = gram(normalized)?;
        let thermal_diffusivity = thermal_conductivity_w_m_k
            / material.density_kg_m3()
            / material.specific_heat_cv_j_kg_k();
        if ![thermal_diffusivity, longitudinal_viscosity, conductive_root]
            .into_iter()
            .all(nonnegative)
            || (thermal_conductivity_w_m_k > 0.0 && thermal_diffusivity == 0.0)
        {
            return Err(invalid(
                "superfluid longitudinal transport coefficient is unrepresentable",
            ));
        }
        Ok(Self {
            material,
            second: second_viscosities,
            conductivity: thermal_conductivity_w_m_k,
            shear_root,
            longitudinal_viscosity,
            conductive_root,
            resistance,
            kinetic_diffusivity,
            thermal_diffusivity,
        })
    }

    pub const fn material(&self) -> &TwoFluidState {
        &self.material
    }
    pub const fn second_viscosities(&self) -> SecondViscosities {
        self.second
    }
    pub const fn thermal_conductivity_w_m_k(&self) -> f64 {
        self.conductivity
    }
    /// `R` in `K u_t = (R u_x)_x`, with `u=(v,w=v_n-v_s)` and
    /// `K=diag(rho,rho_n*rho_s/rho)`. Every entry has units [Pa s].
    pub const fn longitudinal_resistance(&self) -> [[f64; 2]; 2] {
        self.resistance
    }
    /// Symmetric diffusion matrix in kinetic energy coordinates
    /// `(sqrt(rho)*v, sqrt(rho_n*rho_s/rho)*w)`, all entries [m²/s].
    pub const fn kinetic_diffusivity_m2_s(&self) -> [[f64; 2]; 2] {
        self.kinetic_diffusivity
    }
    /// The nonzero conductive eigenvalue, `kappa/(rho Cv)` [m²/s].
    /// Local conduction holds density fixed, so Cp is not the denominator.
    pub const fn thermal_diffusivity_m2_s(&self) -> f64 {
        self.thermal_diffusivity
    }

    /// Evaluate `div(v_n) [1/s]`, `div(q0) [kg/(m³ s)]`, `grad(T') [K/m]`.
    /// Here `q0=rho_s(v_s-v_n)`, opposite to positive counterflow velocity.
    pub fn response(
        &self,
        normal_velocity_gradient_s_inv: f64,
        relative_mass_current_gradient_kg_m3_s: f64,
        temperature_gradient_k_m: f64,
    ) -> Result<LongitudinalResponse, MaterialError> {
        let n = normal_velocity_gradient_s_inv;
        let q = relative_mass_current_gradient_kg_m3_s;
        let t = temperature_gradient_k_m;
        if ![n, q, t].into_iter().all(f64::is_finite) {
            return Err(invalid("superfluid transport gradients must be finite"));
        }
        let [a, b, c] = self.second.factor;
        let shear = self.shear_root * n;
        let bulk = a * n + b * q;
        let relative = c * q;
        let viscous_heat = shear * shear + bulk * bulk + relative * relative;
        let conductive = self.conductive_root * t;
        let conductive_loss = conductive * conductive;
        let heat_flux = -self.conductivity * t;
        let stress = self.longitudinal_viscosity * n + self.second.zeta1_m2_s * q;
        let chemical = -self.second.zeta1_m2_s * n - self.second.zeta3_m5_kg_s * q;
        let entropy_flux = heat_flux / self.material.temperature_k();
        let production = (viscous_heat + conductive_loss) / self.material.temperature_k();
        if ![
            stress,
            chemical,
            heat_flux,
            entropy_flux,
            viscous_heat,
            conductive_loss,
            production,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(invalid("superfluid transport response is unrepresentable"));
        }
        Ok(LongitudinalResponse {
            longitudinal_viscous_stress_pa: stress,
            chemical_potential_correction_j_kg: chemical,
            conductive_heat_flux_w_m2: heat_flux,
            conductive_entropy_flux_w_m2_k: entropy_flux,
            viscous_heat_w_m3: viscous_heat,
            conductive_availability_loss_w_m3: conductive_loss,
            entropy_production_w_m3_k: production,
        })
    }
}

fn gram(rows: [[f64; 2]; 3]) -> Result<[[f64; 2]; 2], MaterialError> {
    let a = rows.iter().map(|r| r[0] * r[0]).sum::<f64>();
    let b = rows.iter().map(|r| r[0] * r[1]).sum::<f64>();
    let c = rows.iter().map(|r| r[1] * r[1]).sum::<f64>();
    if ![a, b, c].into_iter().all(f64::is_finite)
        || (a == 0.0 && rows.iter().any(|r| r[0] != 0.0))
        || (c == 0.0 && rows.iter().any(|r| r[1] != 0.0))
    {
        return Err(invalid(
            "superfluid longitudinal transport matrix is unrepresentable",
        ));
    }
    Ok([[a, b], [b, c]])
}

// Exact normalized product: (binary exponent of its leading bit, significand
// shifted to bit 127). Binary64 products have at most 106 significant bits,
// so u128 retains every bit. Tuple ordering compares products without overflow.
fn exact_product(a: f64, b: f64) -> (i32, u128) {
    if a == 0.0 || b == 0.0 {
        return (i32::MIN, 0);
    }
    fn parts(value: f64) -> (u128, i32) {
        let bits = value.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        if exponent == 0 {
            (u128::from(fraction), -1074)
        } else {
            (u128::from(fraction | (1_u64 << 52)), exponent - 1075)
        }
    }
    let (ma, ea) = parts(a);
    let (mb, eb) = parts(b);
    let product = ma * mb;
    let shift = product.leading_zeros();
    (ea + eb + 127 - shift as i32, product << shift)
}

fn product_difference_fraction(large: (i32, u128), small: (i32, u128)) -> f64 {
    if small.1 == 0 {
        return 1.0;
    }
    let shift = (large.0 - small.0) as u32;
    if shift >= 128 {
        return 1.0;
    }
    // Near cancellation the alignment is exact. For a widely separated small
    // product any discarded bits are far below the binary64 rounding scale.
    (large.1 - (small.1 >> shift)) as f64 / large.1 as f64
}
