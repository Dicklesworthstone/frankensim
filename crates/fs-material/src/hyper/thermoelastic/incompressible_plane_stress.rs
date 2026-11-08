//! Exact elastic-volume plane stress with a moving thermal natural map.
//!
//! Write g(T)=det(G(T)), d=det(a), and F22=g/d. Only the two transverse
//! shears are relaxed. The native remainder supplies Psi_r, Pr and Ar;
//! no bulk penalty is evaluated. The compression-positive multiplier uses
//! P=Pr-p*cof(F), with p=Pr22/d, retaining the original remainder's gauge.
//!
//! The moving constraint changes thermal work. For the reduced physical
//! Helmholtz potential, S=Sr-p*g_T and
//! c_a=c_r-T*(beta_r:F_T+p_T*g_T+p*g_TT). Omitting the pressure terms
//! would make entropy depend on the chosen off-constraint energy extension.
//! This is the entropy reaction of a temperature-dependent material constraint
//! (Negahban, 2006, doi:10.1177/1081286506067092, section 7).
//!
//! Densities use the ORIGINAL reference volume. The law imposes an isothermal
//! elastic-volume constraint; it does not certify unrestricted thermo-acoustic
//! stability, provide a caloric law, or infer irreversible heat/history rates.

use super::{
    ElasticTemperatureDerivatives, FiniteThermoelasticRate, FiniteThermoelasticState, IDENTITY,
    ReferencePotentialDerivatives, ThermalReferenceJet, dot, elastic_deformation_rate, invalid,
    product, trace, unrepresentable,
};
use crate::MaterialError;
use crate::hyper::Hyperelastic;
use crate::hyper::incompressible_plane_stress::{
    IncompressiblePlaneStressKinematicResponse, IncompressiblePlaneStressResponse, PositiveShear,
    generic_state, kinematic_response,
};
use crate::hyper::mixed::cofactor;
use crate::hyper::plane_stress::PlaneStressSettings;
use crate::hyper::reference::ReferencedHyperelastic;

const IN_PLANE: [usize; 4] = [0, 1, 3, 4];
const SHEAR: [usize; 2] = [2, 5];
const TRANSVERSE: [usize; 3] = [2, 5, 8];

// The exact bilinear derivative of every 2x2 cofactor. No inverse, finite
// difference, or subtraction of two neighboring determinants is needed.
fn cofactor_rate(f: &[f64; 9], rate: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        let (i, j) = (ij / 3, ij % 3);
        let a = 3 * ((i + 1) % 3) + (j + 1) % 3;
        let b = 3 * ((i + 2) % 3) + (j + 2) % 3;
        let c = 3 * ((i + 1) % 3) + (j + 2) % 3;
        let d = 3 * ((i + 2) % 3) + (j + 1) % 3;
        f[a].mul_add(rate[b], rate[a] * f[b]) - f[c].mul_add(rate[d], rate[c] * f[d])
    })
}

/// Physical thermodynamics on a stable, exactly constrained plane branch.
/// The private remainder is only a derivative carrier: its unconstrained
/// entropy and stress are never exposed as the physical constrained response.
#[derive(Clone, Debug)]
pub struct IncompressibleThermoelasticPlaneStressState {
    remainder: FiniteThermoelasticState,
    kinematics: IncompressiblePlaneStressKinematicResponse,
    natural_volume_derivative_per_k: f64,
    natural_volume_second_derivative_per_k2: f64,
    entropy_j_m3_k: f64,
    internal_energy_j_m3: f64,
    mechanical_heat_capacity_j_m3_k: f64,
    stress_temperature_derivative_pa_k: [f64; 4],
    transverse_temperature_derivative_per_k: [f64; 3],
    pressure_temperature_derivative_pa_k: f64,
}

/// Reversible rates on the same exact-volume thermal branch.
/// Pressure retains the native multiplier convention, while material stress,
/// energy and entropy rates all use the physical constrained response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncompressibleThermoelasticPlaneStressRate {
    material: FiniteThermoelasticRate,
    pressure_rate_pa_s: f64,
    natural_volume_ratio_rate_per_s: f64,
}

impl IncompressibleThermoelasticPlaneStressRate {
    pub const fn material(&self) -> &FiniteThermoelasticRate {
        &self.material
    }
    pub const fn pressure_rate_pa_s(&self) -> f64 {
        self.pressure_rate_pa_s
    }
    /// d det(G)/dt, equal to d det(F)/dt on the constrained trajectory [1/s].
    pub const fn natural_volume_ratio_rate_per_s(&self) -> f64 {
        self.natural_volume_ratio_rate_per_s
    }
}

impl IncompressibleThermoelasticPlaneStressState {
    /// Solve the original exact-volume shear equilibrium, then retain its
    /// constrained mechanical and thermal derivatives. Elastic coefficients
    /// are held; G, G_T and G_TT come from the supplied natural-map owner.
    ///
    /// The optional seed is [F02,F12]; F22=det(G)/det(a) is imposed exactly.
    /// Native stretch, orientation, shear-stability and iteration gates apply.
    pub fn try_new(
        law: &Hyperelastic,
        in_plane: &[f64; 4],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        initial_shear: Option<[f64; 2]>,
        settings: PlaneStressSettings,
    ) -> Result<Self, MaterialError> {
        Self::with_elastic_temperature_derivatives(
            law,
            in_plane,
            temperature_k,
            reference_jet,
            initial_shear,
            settings,
            |_| Ok(ElasticTemperatureDerivatives::default()),
        )
    }

    // Only admitted native source adapters supply explicit coefficient jets.
    // This is not a public callback for an arbitrary incremental potential.
    pub(super) fn with_elastic_temperature_derivatives(
        law: &Hyperelastic,
        in_plane: &[f64; 4],
        temperature_k: f64,
        reference_jet: ThermalReferenceJet,
        initial_shear: Option<[f64; 2]>,
        settings: PlaneStressSettings,
        derivatives: impl FnOnce(&[f64; 9]) -> Result<ElasticTemperatureDerivatives, MaterialError>,
    ) -> Result<Self, MaterialError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(invalid("positive finite absolute temperature required"));
        }
        // Re-admit public parameters without evaluating their bulk penalty.
        let law = Hyperelastic::new(law.model.clone(), law.stretch_limit)?;
        let material = ReferencedHyperelastic::new(&law, *reference_jet.reference());
        let (plane, accepted) = generic_state(&material, in_plane, initial_shear, settings)?;
        let accepted = accepted.ok_or_else(|| invalid("missing accepted remainder tangent"))?;
        let elastic_f = reference_jet
            .reference()
            .elastic_deformation(&plane.deformation_gradient);
        let partial = derivatives(&elastic_f)?;
        let remainder = FiniteThermoelasticState::with_reference_potential(
            &plane.deformation_gradient,
            temperature_k,
            reference_jet,
            ReferencePotentialDerivatives {
                free_energy_j_m3: accepted.potential_pa,
                first_piola_pa: accepted.first_piola_pa,
                isothermal_tangent_pa: accepted.tangent_pa,
            },
            partial,
        )?;
        let kinematics = kinematic_response(plane, accepted, in_plane)?;
        Self::from_prepared(remainder, kinematics)
    }

    fn from_prepared(
        remainder: FiniteThermoelasticState,
        kinematics: IncompressiblePlaneStressKinematicResponse,
    ) -> Result<Self, MaterialError> {
        let jet = remainder.reference_jet();
        let h = jet.reference().elastic_deformation(&IDENTITY);
        let d = product(&h, jet.first_derivative_per_k());
        let d2 = product(&d, &d);
        let trace_d = trace(&d);
        let natural_volume_derivative_per_k = jet.reference().volume_ratio() * trace_d;
        let natural_volume_second_derivative_per_k2 = jet.reference().volume_ratio()
            * (trace_d * trace_d - trace(&d2)
                + trace(&product(&h, jet.second_derivative_per_k2())));
        let f = kinematics.response.deformation_gradient;
        let determinant = f[0] * f[4] - f[1] * f[3];
        let thickness_t = natural_volume_derivative_per_k / determinant;
        let a = remainder.isothermal_tangent_pa();
        let beta = remainder.stress_temperature_derivative_pa_k();
        let block = PositiveShear::new(SHEAR.map(|i| SHEAR.map(|j| a[i][j])))?;
        let forcing = SHEAR.map(|i| a[i][8].mul_add(thickness_t, beta[i]));
        let shear_t = block.solve(forcing.map(|x| -x))?;
        let transverse_temperature_derivative_per_k = [shear_t[0], shear_t[1], thickness_t];
        let mut f_t = [0.0; 9];
        for (i, component) in TRANSVERSE.into_iter().enumerate() {
            f_t[component] = transverse_temperature_derivative_per_k[i];
        }
        let pressure_temperature_derivative_pa_k = (dot(&a[8], &f_t) + beta[8]) / determinant;
        let pressure = kinematics.response.pressure_pa;
        let cof = cofactor(&f);
        let cof_t = cofactor_rate(&f, &f_t);
        let stress_temperature_derivative_pa_k = IN_PLANE.map(|i| {
            (-pressure).mul_add(
                cof_t[i],
                (-pressure_temperature_derivative_pa_k).mul_add(cof[i], dot(&a[i], &f_t) + beta[i]),
            )
        });
        let entropy_j_m3_k =
            (-pressure).mul_add(natural_volume_derivative_per_k, remainder.entropy_j_m3_k());
        let thermal_curvature = pressure.mul_add(
            natural_volume_second_derivative_per_k2,
            pressure_temperature_derivative_pa_k
                .mul_add(natural_volume_derivative_per_k, dot(beta, &f_t)),
        );
        let mechanical_heat_capacity_j_m3_k = (-remainder.temperature_k()).mul_add(
            thermal_curvature,
            remainder.mechanical_heat_capacity_j_m3_k(),
        );
        let internal_energy_j_m3 = remainder
            .temperature_k()
            .mul_add(entropy_j_m3_k, remainder.free_energy_j_m3());
        if [
            natural_volume_derivative_per_k,
            natural_volume_second_derivative_per_k2,
            pressure_temperature_derivative_pa_k,
            entropy_j_m3_k,
            internal_energy_j_m3,
            mechanical_heat_capacity_j_m3_k,
        ]
        .iter()
        .chain(&transverse_temperature_derivative_per_k)
        .chain(&stress_temperature_derivative_pa_k)
        .any(|value| !value.is_finite())
        {
            return Err(unrepresentable());
        }
        Ok(Self {
            remainder,
            kinematics,
            natural_volume_derivative_per_k,
            natural_volume_second_derivative_per_k2,
            entropy_j_m3_k,
            internal_energy_j_m3,
            mechanical_heat_capacity_j_m3_k,
            stress_temperature_derivative_pa_k,
            transverse_temperature_derivative_per_k,
            pressure_temperature_derivative_pa_k,
        })
    }

    /// Original exact-volume response, including physical stress and pressure.
    pub const fn response(&self) -> &IncompressiblePlaneStressResponse {
        &self.kinematics.response
    }
    /// Mechanical derivatives holding temperature and natural volume fixed.
    pub const fn kinematics(&self) -> &IncompressiblePlaneStressKinematicResponse {
        &self.kinematics
    }
    pub const fn reference_jet(&self) -> &ThermalReferenceJet {
        self.remainder.reference_jet()
    }
    pub const fn temperature_k(&self) -> f64 {
        self.remainder.temperature_k()
    }
    pub const fn elastic_deformation_gradient(&self) -> &[f64; 9] {
        self.remainder.elastic_deformation_gradient()
    }
    pub const fn free_energy_j_m3(&self) -> f64 {
        self.remainder.free_energy_j_m3()
    }
    /// -dPsi_reduced/dT, including -p*d det(G)/dT [J/(m³ K)].
    pub const fn entropy_j_m3_k(&self) -> f64 {
        self.entropy_j_m3_k
    }
    pub const fn internal_energy_j_m3(&self) -> f64 {
        self.internal_energy_j_m3
    }
    /// Mechanical -T*d²Psi_reduced/dT² [J/(m³ K)]. This may be negative;
    /// a separate caloric owner must admit the complete heat capacity.
    pub const fn mechanical_heat_capacity_j_m3_k(&self) -> f64 {
        self.mechanical_heat_capacity_j_m3_k
    }
    /// dP_in_plane/dT along relaxed exact-volume equilibrium [Pa/K].
    pub const fn stress_temperature_derivative_pa_k(&self) -> &[f64; 4] {
        &self.stress_temperature_derivative_pa_k
    }
    /// d[F02,F12,F22]/dT at fixed in-plane total deformation [1/K].
    pub const fn transverse_temperature_derivative_per_k(&self) -> &[f64; 3] {
        &self.transverse_temperature_derivative_per_k
    }
    /// Native pressure-multiplier derivative at fixed in-plane F [Pa/K].
    pub const fn pressure_temperature_derivative_pa_k(&self) -> f64 {
        self.pressure_temperature_derivative_pa_k
    }
    pub const fn natural_volume_derivative_per_k(&self) -> f64 {
        self.natural_volume_derivative_per_k
    }
    pub const fn natural_volume_second_derivative_per_k2(&self) -> f64 {
        self.natural_volume_second_derivative_per_k2
    }

    /// Physical rates along prescribed in-plane Fdot and Tdot.
    ///
    /// Thickness follows det(F)=det(G(T)), and both shears follow native
    /// stationarity. Pdot includes -pdot*cof(F)-p*dcof(F)/dt. Mechanical
    /// entropy, Helmholtz and internal-energy rates use the CONSTRAINED
    /// thermal derivatives; raw remainder entropy is not physical entropy.
    /// No irreversible heat, accepted history update or caloric rate is added.
    pub fn rate(
        &self,
        in_plane_rate_per_s: [f64; 4],
        temperature_rate_k_s: f64,
    ) -> Result<IncompressibleThermoelasticPlaneStressRate, MaterialError> {
        if !temperature_rate_k_s.is_finite() {
            return Err(invalid("finite temperature rate required"));
        }
        let mut fdot = self
            .kinematics
            .deformation_gradient_rate(in_plane_rate_per_s)?;
        for (i, component) in TRANSVERSE.into_iter().enumerate() {
            fdot[component] = self.transverse_temperature_derivative_per_k[i]
                .mul_add(temperature_rate_k_s, fdot[component]);
        }
        let pressure_rate_pa_s = self.pressure_temperature_derivative_pa_k.mul_add(
            temperature_rate_k_s,
            self.kinematics.pressure_rate(in_plane_rate_per_s)?,
        );
        let natural_volume_ratio_rate_per_s =
            self.natural_volume_derivative_per_k * temperature_rate_k_s;
        let response = self.response();
        let cof = cofactor(&response.deformation_gradient);
        let cof_dot = cofactor_rate(&response.deformation_gradient, &fdot);
        let a = self.remainder.isothermal_tangent_pa();
        let beta = self.remainder.stress_temperature_derivative_pa_k();
        let first_piola_rate_pa_s = core::array::from_fn(|i| {
            (-response.pressure_pa).mul_add(
                cof_dot[i],
                (-pressure_rate_pa_s).mul_add(
                    cof[i],
                    beta[i].mul_add(temperature_rate_k_s, dot(&a[i], &fdot)),
                ),
            )
        });
        let stress_power = (0..4).fold(0.0, |sum, i| {
            response.piola[i].mul_add(in_plane_rate_per_s[i], sum)
        });
        let entropy_rate_w_m3_k = (self.mechanical_heat_capacity_j_m3_k / self.temperature_k())
            .mul_add(
                temperature_rate_k_s,
                -(0..4).fold(0.0, |sum, i| {
                    self.stress_temperature_derivative_pa_k[i].mul_add(in_plane_rate_per_s[i], sum)
                }),
            );
        let free_energy_rate_w_m3 =
            (-self.entropy_j_m3_k).mul_add(temperature_rate_k_s, stress_power);
        let internal_energy_rate_w_m3 = self
            .temperature_k()
            .mul_add(entropy_rate_w_m3_k, stress_power);
        let elastic_deformation_gradient_rate_per_s = elastic_deformation_rate(
            self.elastic_deformation_gradient(),
            self.reference_jet(),
            &fdot,
            temperature_rate_k_s,
        );
        if [
            pressure_rate_pa_s,
            natural_volume_ratio_rate_per_s,
            entropy_rate_w_m3_k,
            free_energy_rate_w_m3,
            internal_energy_rate_w_m3,
        ]
        .iter()
        .chain(&fdot)
        .chain(&elastic_deformation_gradient_rate_per_s)
        .chain(&first_piola_rate_pa_s)
        .any(|value| !value.is_finite())
        {
            return Err(unrepresentable());
        }
        Ok(IncompressibleThermoelasticPlaneStressRate {
            material: FiniteThermoelasticRate {
                deformation_gradient_rate_per_s: fdot,
                elastic_deformation_gradient_rate_per_s,
                temperature_rate_k_s,
                first_piola_rate_pa_s,
                free_energy_rate_w_m3,
                entropy_rate_w_m3_k,
                internal_energy_rate_w_m3,
            },
            pressure_rate_pa_s,
            natural_volume_ratio_rate_per_s,
        })
    }
}
