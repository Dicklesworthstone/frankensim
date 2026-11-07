//! Physical ply shear on the existing condensed MITC trial space.
//!
//! The material owner prepares the equilibrium profile once per source stack.
//! Each facet caches only its original strain observation and frame rotation.
//! A real or coherent harmonic displacement is observed once; subsequent
//! layer/point queries allocate nothing and retain the complete native field.
//! Product Bernstein bounds cover the triangle and the whole quadratic ply
//! profile, including shear maxima between the physical layer faces.
//! The result inherits the source's declared two-probe bending approximation;
//! it is not a free-edge, delamination or transverse-normal-stress solution.

use super::strain::{
    HarmonicElementStrain, HarmonicShellStrainPoint, PreparedElementStrain, RealElementStrain,
    ShellStrainPoint, hypot, plane_stress_tensor, tensor_cycle_norm, validate_layer,
};
use super::{MitcShellModel, PlateError, ShearShellSection};
use fs_material::laminate::{LaminateShearRecovery, LayerShearResponse};
use fs_math::c64::C64;
use fs_math::det;

fn invalid(what: &'static str) -> PlateError {
    PlateError::BadSection { what }
}

/// Original facet observation bound to the exact oriented laminate operators.
pub struct PreparedLaminateShearFacet<'a> {
    model: &'a MitcShellModel,
    element: usize,
    source: &'a LaminateShearRecovery,
    observation: PreparedElementStrain,
    sin: f64,
    cos: f64,
}

/// Actual native strain at one instant, with the original layer recovery.
/// It can be queried repeatedly without re-observing the nodal displacement.
pub struct LaminateShearField<'a> {
    source: &'a LaminateShearRecovery,
    native: RealElementStrain,
    sin: f64,
    cos: f64,
    // The thermal owner supplies total in-plane strain and its director offset,
    // but this field does not retain the physical normal free-strain closure.
    thermal_origin: bool,
}

/// Coherent peak complex strain on the native trial space, using the original
/// real elastic ply recovery at every common harmonic phase.
pub struct HarmonicLaminateShearField<'a> {
    source: &'a LaminateShearRecovery,
    native: HarmonicElementStrain,
    sin: f64,
    cos: f64,
}

/// The same local response in original ply/section axes and facet ex/ey axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FacetLayerShearResponse {
    /// Original material response; its xz/yz components use SECTION axes.
    pub layer: LayerShearResponse,
    /// Physical [tau_ex,n, tau_ey,n] in facet axes [Pa].
    pub facet_stress_pa: [f64; 2],
    /// Physical engineering [gamma_ex,n, gamma_ey,n] in facet axes.
    pub facet_engineering_strain: [f64; 2],
}

/// Peak complex ply response with the exp(i omega t) convention. Every tensor
/// component retains its original phase; no componentwise amplitude is used
/// to recover another component or to form the mean stored energy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicFacetLayerShearResponse {
    pub z_m: f64,
    pub section_stress_pa: [C64; 2],
    pub section_engineering_strain: [C64; 2],
    pub material_stress_pa: [C64; 2],
    pub material_engineering_strain: [C64; 2],
    pub facet_stress_pa: [C64; 2],
    pub facet_engineering_strain: [C64; 2],
    /// Cycle-mean recoverable shear energy, Re(tau dot gamma*) / 4 [J/m^3].
    /// This elastic recovery does not infer a viscoelastic loss.
    pub average_energy_density_j_m3: f64,
}

/// Whole-triangle, whole-ply strain bound. It combines the original membrane
/// and bending field with the physical recovered transverse shear. Real
/// queries describe one instant; harmonic queries additionally cover all phases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LaminateLayerStrainBound {
    pub frobenius_upper_bound: f64,
    /// Independent component envelopes [xx, yy, zz, xy, yz, zx] in facet axes.
    pub component_amplitudes: [f64; 6],
    /// Native triangle Bernstein control, ordered as `TensorStrainBound`.
    pub ruling_control_point: usize,
    /// 0/2 are the ply's bottom/top coefficients; 1 is its middle degree-two
    /// Bernstein coefficient, NOT the strain at the physical midpoint.
    /// Neither ruling control is a claim about the actual maximum location.
    pub ruling_thickness_control: usize,
}

impl MitcShellModel {
    /// Prepare local ply shear for one original facet and laminate source.
    /// `section_angle_rad` has exactly the meaning used by
    /// [`ShearShellSection::from_laminate`]. All elastic blocks and inertia
    /// moments must match the retained assembly; a stiffness-only lookalike or
    /// changed reference surface refuses. No solve or motion update occurs.
    pub fn prepare_laminate_shear<'a>(
        &'a self,
        element: usize,
        source: &'a LaminateShearRecovery,
        section_angle_rad: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedLaminateShearFacet<'a>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        let actual = self
            .sections
            .get(element)
            .ok_or(invalid("laminate shear facet out of range"))?;
        admit(actual, source, section_angle_rad)?;
        let observation = PreparedElementStrain::new(
            &self.mesh,
            element,
            &self.directors,
            actual,
            &mut cancelled,
        )?;
        Ok(PreparedLaminateShearFacet {
            model: self,
            element,
            source,
            observation,
            sin: det::sin(section_angle_rad),
            cos: det::cos(section_angle_rad),
        })
    }
}

fn admit(
    actual: &ShearShellSection,
    source: &LaminateShearRecovery,
    angle: f64,
) -> Result<(), PlateError> {
    let expected = ShearShellSection::from_laminate(source.section(), angle)?;
    if actual.coupled_matrices() != expected.coupled_matrices()
        || actual.areal_mass_kg_m2() != expected.areal_mass_kg_m2()
        || actual.mass_first_moment_kg_m() != expected.mass_first_moment_kg_m()
        || actual.rotary_inertia_kg() != expected.rotary_inertia_kg()
    {
        return Err(invalid(
            "shear recovery must match the original laminate, reference and facet orientation",
        ));
    }
    Ok(())
}

impl<'a> PreparedLaminateShearFacet<'a> {
    /// Observe a full physical world-coordinate displacement/rotation snapshot.
    /// Validation reads only this facet's nodes, so E facet observations cost
    /// O(E), not O(E times all nodes). Support coordinates must be expanded by
    /// the original model before this call. Thermal director offsets require
    /// [`super::thermal::PreparedThermalShell::laminate_shear_field`] instead.
    pub fn evaluate(
        &self,
        displacement_6n: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShearField<'a>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if displacement_6n.len() != 6 * self.model.mesh.nodes.len() {
            return Err(invalid(
                "one physical displacement/rotation per original shell DOF required",
            ));
        }
        let tri = self.model.mesh.tris[self.element];
        let nodal = core::array::from_fn(|i| displacement_6n[6 * tri[i / 6] + i % 6]);
        let native = self.observation.real_field(&nodal, &mut cancelled)?;
        Ok(LaminateShearField {
            source: self.source,
            native,
            sin: self.sin,
            cos: self.cos,
            thermal_origin: false,
        })
    }

    /// Observe one coherent peak complex physical displacement/rotation field.
    /// As for real snapshots, only this facet's nodes are inspected and the
    /// original elastic bubble condensation is retained. No frequency-specific
    /// constitutive reduction or harmonic loss is introduced by observation.
    pub fn evaluate_harmonic(
        &self,
        displacement_6n: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<HarmonicLaminateShearField<'a>, PlateError> {
        if cancelled() {
            return Err(PlateError::Cancelled);
        }
        if displacement_6n.len() != 6 * self.model.mesh.nodes.len() {
            return Err(invalid(
                "one complex physical displacement/rotation per original shell DOF required",
            ));
        }
        let tri = self.model.mesh.tris[self.element];
        let real = core::array::from_fn(|i| displacement_6n[6 * tri[i / 6] + i % 6].re);
        let imaginary = core::array::from_fn(|i| displacement_6n[6 * tri[i / 6] + i % 6].im);
        let native = self
            .observation
            .harmonic_field(&real, &imaginary, &mut cancelled)?;
        Ok(HarmonicLaminateShearField {
            source: self.source,
            native,
            sin: self.sin,
            cos: self.cos,
        })
    }
}

impl<'a> LaminateShearField<'a> {
    pub(super) fn from_thermal(
        section: &ShearShellSection,
        source: &'a LaminateShearRecovery,
        angle: f64,
        native: RealElementStrain,
    ) -> Result<Self, PlateError> {
        admit(section, source, angle)?;
        Ok(Self {
            source,
            native,
            sin: det::sin(angle),
            cos: det::cos(angle),
            thermal_origin: true,
        })
    }

    /// Complete original native strain polynomial, including any thermal
    /// director offset supplied by the native thermal owner.
    pub const fn native_strain(&self) -> &RealElementStrain {
        &self.native
    }

    /// Physical shear at an interior barycentric location and a fraction
    /// through an original ply. Stress is continuous across bonded interfaces;
    /// strains follow each ply's original transverse compliance.
    pub fn at_barycentric(
        &self,
        barycentric: [f64; 3],
        layer: usize,
        fraction: f64,
    ) -> Result<FacetLayerShearResponse, PlateError> {
        let gamma = self.native.at_barycentric(barycentric)?.transverse_shear;
        recover(self.source, self.sin, self.cos, gamma, layer, fraction)
    }

    /// Full tensor bound with an explicitly selected zero-eigenstrain isotropic
    /// plane-stress closure: epsilon_zz = -nu/(1-nu)*(epsilon_xx+epsilon_yy).
    /// The caller supplies that physical isotropic ply's Poisson ratio; this
    /// closure is not inferred from a general orthotropic reduced source.
    ///
    /// Thermal-derived fields refuse: their total normal strain also requires
    /// the actual thermal free-strain components, which this observation does
    /// not retain. Their native in-plane and recovered shear queries remain
    /// available. Refusal prevents free thermal expansion from being mistaken
    /// for a mechanical Poisson contraction.
    pub fn isotropic_layer_bound(
        &self,
        layer: usize,
        poisson_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LaminateLayerStrainBound, PlateError> {
        if self.thermal_origin {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            return Err(invalid(
                "thermal full-tensor strain requires its original normal free-strain closure",
            ));
        }
        let controls = self
            .native
            .bernstein_controls()
            .map(|real| HarmonicShellStrainPoint {
                real,
                imaginary: ShellStrainPoint {
                    membrane: [0.0; 3],
                    gradient_per_m: [0.0; 3],
                    transverse_shear: [0.0; 2],
                },
            });
        layer_bound(
            self.source,
            self.sin,
            self.cos,
            &controls,
            layer,
            poisson_ratio,
            cancelled,
        )
    }
}

impl HarmonicLaminateShearField<'_> {
    pub const fn native_strain(&self) -> &HarmonicElementStrain {
        &self.native
    }

    /// Coherent complex stress and strain at one original ply point.
    pub fn at_barycentric(
        &self,
        barycentric: [f64; 3],
        layer: usize,
        fraction: f64,
    ) -> Result<HarmonicFacetLayerShearResponse, PlateError> {
        let point = self.native.at_barycentric(barycentric)?;
        let real = recover(
            self.source,
            self.sin,
            self.cos,
            point.real.transverse_shear,
            layer,
            fraction,
        )?;
        let imaginary = recover(
            self.source,
            self.sin,
            self.cos,
            point.imaginary.transverse_shear,
            layer,
            fraction,
        )?;
        let combine = |r: [f64; 2], i: [f64; 2]| core::array::from_fn(|j| C64::new(r[j], i[j]));
        Ok(HarmonicFacetLayerShearResponse {
            z_m: real.layer.z_m,
            section_stress_pa: combine(real.layer.stress_pa, imaginary.layer.stress_pa),
            section_engineering_strain: combine(
                real.layer.engineering_strain,
                imaginary.layer.engineering_strain,
            ),
            material_stress_pa: combine(
                real.layer.material_stress_pa,
                imaginary.layer.material_stress_pa,
            ),
            material_engineering_strain: combine(
                real.layer.material_engineering_strain,
                imaginary.layer.material_engineering_strain,
            ),
            facet_stress_pa: combine(real.facet_stress_pa, imaginary.facet_stress_pa),
            facet_engineering_strain: combine(
                real.facet_engineering_strain,
                imaginary.facet_engineering_strain,
            ),
            average_energy_density_j_m3: f64::midpoint(
                real.layer.energy_density_j_m3,
                imaginary.layer.energy_density_j_m3,
            ),
        })
    }

    /// Cover the complete native triangle, physical ply thickness and common
    /// harmonic phase. The full tensor uses the same explicitly authored
    /// zero-eigenstrain isotropic plane-stress closure as `LaminateShearField`;
    /// it is not a thermal or general orthotropic normal-strain recovery.
    /// Product Bernstein convexity bounds the represented field, not continuum
    /// discretization error, floating-point roundoff or material failure.
    pub fn isotropic_layer_bound(
        &self,
        layer: usize,
        poisson_ratio: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<LaminateLayerStrainBound, PlateError> {
        layer_bound(
            self.source,
            self.sin,
            self.cos,
            &self.native.bernstein_controls(),
            layer,
            poisson_ratio,
            cancelled,
        )
    }
}

fn recover(
    source: &LaminateShearRecovery,
    sin: f64,
    cos: f64,
    gamma: [f64; 2],
    layer: usize,
    fraction: f64,
) -> Result<FacetLayerShearResponse, PlateError> {
    let section_gamma = [
        cos * gamma[0] + sin * gamma[1],
        -sin * gamma[0] + cos * gamma[1],
    ];
    let layer = source
        .from_section_strain(layer, fraction, section_gamma)
        .map_err(|_| invalid("original laminate shear recovery refused"))?;
    let to_facet = |v: [f64; 2]| [cos * v[0] - sin * v[1], sin * v[0] + cos * v[1]];
    let facet_stress_pa = to_facet(layer.stress_pa);
    let facet_engineering_strain = to_facet(layer.engineering_strain);
    if facet_stress_pa
        .iter()
        .chain(&facet_engineering_strain)
        .any(|x| !x.is_finite())
    {
        return Err(invalid("unrepresentable rotated laminate shear response"));
    }
    Ok(FacetLayerShearResponse {
        layer,
        facet_stress_pa,
        facet_engineering_strain,
    })
}

fn layer_bound(
    source: &LaminateShearRecovery,
    sin: f64,
    cos: f64,
    controls: &[HarmonicShellStrainPoint; 6],
    layer: usize,
    poisson_ratio: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<LaminateLayerStrainBound, PlateError> {
    if cancelled() {
        return Err(PlateError::Cancelled);
    }
    let geometry = source
        .section()
        .layer_integrals()
        .get(layer)
        .ok_or(invalid("laminate strain layer out of range"))?;
    validate_layer([geometry.z_bottom_m, geometry.z_top_m], poisson_ratio)?;
    let maps = source
        .strain_bernstein_controls(layer)
        .map_err(|_| invalid("original laminate shear strain profile refused"))?;
    let z = [
        geometry.z_bottom_m,
        f64::midpoint(geometry.z_bottom_m, geometry.z_top_m),
        geometry.z_top_m,
    ];
    let mut bound = LaminateLayerStrainBound {
        frobenius_upper_bound: 0.0,
        component_amplitudes: [0.0; 6],
        ruling_control_point: 0,
        ruling_thickness_control: 0,
    };
    for (control, point) in controls.iter().enumerate() {
        for (thickness_control, map) in maps.iter().enumerate() {
            if cancelled() {
                return Err(PlateError::Cancelled);
            }
            let tensor = |point: &ShellStrainPoint| {
                let gamma = point.transverse_shear;
                let section = [
                    cos * gamma[0] + sin * gamma[1],
                    -sin * gamma[0] + cos * gamma[1],
                ];
                let physical = [
                    map[0] * section[0] + map[1] * section[1],
                    map[2] * section[0] + map[3] * section[1],
                ];
                let coefficients = [
                    point.membrane[0],
                    point.membrane[1],
                    point.membrane[2],
                    point.gradient_per_m[0],
                    point.gradient_per_m[1],
                    point.gradient_per_m[2],
                    cos * physical[0] - sin * physical[1],
                    sin * physical[0] + cos * physical[1],
                ];
                plane_stress_tensor(&coefficients, z[thickness_control], poisson_ratio)
            };
            let real = tensor(&point.real)?;
            let imaginary = tensor(&point.imaginary)?;
            let peak = tensor_cycle_norm(&real, &imaginary)?;
            if peak > bound.frobenius_upper_bound {
                bound.frobenius_upper_bound = peak;
                bound.ruling_control_point = control;
                bound.ruling_thickness_control = thickness_control;
            }
            for i in 0..6 {
                bound.component_amplitudes[i] =
                    bound.component_amplitudes[i].max(hypot(real[i], imaginary[i])?);
            }
        }
    }
    Ok(bound)
}
