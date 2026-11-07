//! Complete harmonic section work in an explicitly fixed physical shear trial field.
//!
//! The eight work coordinates are membrane strain, curvature and the two
//! transverse-shear energy coordinates of an ORIGINAL elastic shear recovery.
//! Actual layer stresses use complete full-tensor Elastic/Prony sources, with
//! only the transverse normal stress eliminated. No in-plane/shear coupling
//! is discarded to fit an A/B/D plus independent-H representation.

use super::{LaminateError, LaminateSection, LaminateShearRecovery, LayerLaw};
use crate::plane_stress_viscosity::MaterialPlane;
use crate::visco::tensor::{HarmonicWorkSum as Sum, NormalStressHarmonic, TensorProny};
use fs_math::c64::C64;
use std::sync::Arc;

/// Source preparation or harmonic-section refusal. No partial result is returned.
#[derive(Clone, Debug, PartialEq)]
pub enum HarmonicLaminateError {
    /// Original laminate/material source refusal, preserving its layer index.
    Source(LaminateError),
    /// Incompatible sources, geometry or nonrepresentable physical work.
    Invalid { what: &'static str },
    /// Explicit caller allowance was exceeded before the relevant allocation.
    Budget { quantity: &'static str, requested: u128, maximum: u128 },
    /// Caller cancellation checkpoint.
    Cancelled,
}

impl core::fmt::Display for HarmonicLaminateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(source) => write!(f, "FS-MATERIAL-HARMONIC-LAMINATE-SOURCE: {source}"),
            Self::Invalid { what } => write!(f, "FS-MATERIAL-HARMONIC-LAMINATE-INVALID: {what}"),
            Self::Budget { quantity, requested, maximum } => write!(
                f, "FS-MATERIAL-HARMONIC-LAMINATE-BUDGET: {quantity}: {requested} > {maximum}"
            ),
            Self::Cancelled => write!(f, "FS-MATERIAL-HARMONIC-LAMINATE-CANCELLED"),
        }
    }
}
impl std::error::Error for HarmonicLaminateError {}
impl From<LaminateError> for HarmonicLaminateError {
    fn from(value: LaminateError) -> Self { Self::Source(value) }
}

/// Allowance for one source preparation, prepared frequency or tangent query.
/// Already owned original laminate/shear sources and fixed stack are excluded.
#[derive(Clone, Copy, Debug)]
pub struct HarmonicLaminateBudget {
    pub maximum_layers: usize,
    /// Additional retained numerical heap payload; allocator metadata excluded.
    pub maximum_payload_bytes: usize,
    /// Conservative operation allowance for this call, including every branch.
    pub maximum_work: u64,
}

type Result<T> = core::result::Result<T, HarmonicLaminateError>;
type StrainMap = [[f64; 8]; 5];

#[derive(Clone, Debug)]
struct LayerTrial {
    // Only an ORIGINAL elastic law needs an exact lossless tensor interface.
    // A Prony law remains untouched in the retained laminate.
    elastic: Option<Box<TensorProny>>,
    controls: [StrainMap; 3],
    quadrature: [(StrainMap, f64); 3],
}

/// Immutable actual layer sources and a separately declared fixed shear trial.
///
/// The caller supplies an existing `LaminateShearRecovery`; no reference
/// material, surrogate, shear correction or profile is inferred from a Prony
/// law. Layer thicknesses, reference faces and material angles must agree
/// exactly. Actual density moments come solely from the actual laminate.
///
/// The reference profile is a FIXED Galerkin strain approximation at every
/// queried frequency. Its elastic two-probe traction equilibration belongs
/// only to its original reference source. General actual anisotropic coupling
/// can create exterior/interface shear tractions even where the reference
/// traction vanishes; this API does not erase those stresses or claim a new
/// frequency-dependent warping/equilibrium solution.
///
/// This is a complete constitutive section, not a finite-element assembler or
/// time integrator. Native enriched MITC coupling must retain its original
/// coordinates, physical mass and numerical gauge when consuming this operator.
#[derive(Clone, Debug)]
pub struct HarmonicLaminate {
    section: Arc<LaminateSection>,
    shear_trial: Arc<LaminateShearRecovery>,
    layers: Vec<LayerTrial>,
}

fn invalid(what: &'static str) -> HarmonicLaminateError {
    HarmonicLaminateError::Invalid { what }
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() { Err(HarmonicLaminateError::Cancelled) } else { Ok(()) }
}
fn cap(quantity: &'static str, requested: u128, maximum: u128) -> Result<()> {
    if requested > maximum {
        Err(HarmonicLaminateError::Budget { quantity, requested, maximum })
    } else {
        Ok(())
    }
}
fn material_error(layer: usize, source: crate::visco::ViscoError) -> HarmonicLaminateError {
    LaminateError::Layer { layer, source }.into()
}

impl HarmonicLaminate {
    /// Bind original material sources to the explicitly supplied fixed profile.
    ///
    /// Reduced-only actual layers refuse: they cannot recover the full original
    /// normal strain or declare its physical strain domain. A full elastic
    /// layer uses its exact lossless tensor interface with the same strain cap.
    /// A Prony layer is borrowed unchanged, with every original branch and cap.
    pub fn try_new(
        section: Arc<LaminateSection>,
        shear_trial: Arc<LaminateShearRecovery>,
        budget: HarmonicLaminateBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        check(&mut cancelled)?;
        let count = section.layers().len();
        cap("layers", count as u128, budget.maximum_layers as u128)?;
        let reference = shear_trial.section();
        if count != reference.layers().len() {
            return Err(invalid("actual and fixed-profile layer counts must agree"));
        }
        let elastic = section.layers().iter()
            .filter(|layer| matches!(&layer.law, LayerLaw::Elastic(_))).count();
        let bytes = (count as u128) * core::mem::size_of::<LayerTrial>() as u128
            + (elastic as u128) * core::mem::size_of::<TensorProny>() as u128;
        cap("source payload bytes", bytes, budget.maximum_payload_bytes as u128)?;
        cap("source work", 20_000 * count as u128, budget.maximum_work as u128)?;
        let mut layers = Vec::with_capacity(count);
        for (layer, actual) in section.layers().iter().enumerate() {
            check(&mut cancelled)?;
            let original = &reference.layers()[layer];
            let geometry = &section.layer_integrals()[layer];
            let reference_geometry = &reference.layer_integrals()[layer];
            if actual.thickness_m != original.thickness_m
                || actual.orientation_rad != original.orientation_rad
                || geometry.z_bottom_m != reference_geometry.z_bottom_m
                || geometry.z_top_m != reference_geometry.z_top_m
            {
                return Err(invalid("actual and fixed-profile layer geometry and axes must agree"));
            }
            let elastic = match &actual.law {
                LayerLaw::Prony(_) => None,
                LayerLaw::Elastic(law) => Some(Box::new(
                    TensorProny::from_tensor_shear(law.stiffness(), Vec::new(), law.strain_limit)
                        .map_err(|_| invalid("original elastic layer has no full tensor strain interface"))?,
                )),
                LayerLaw::PlaneStress { .. } => {
                    return Err(invalid("actual harmonic layers require full original tensor sources"));
                }
            };
            let in_plane = section.principal_strain_transform(layer)?;
            // Both the profile and its material rotation belong to the
            // original shear owner; there is no second shear convention here.
            let shear = shear_trial.material_strain_bernstein_controls(layer)?;
            let bottom = geometry.z_bottom_m;
            let thickness = geometry.z_top_m - bottom;
            if !thickness.is_finite() || thickness <= 0.0 {
                return Err(invalid("original layer thickness is unrepresentable"));
            }
            let z = [bottom, bottom + 0.5 * thickness, geometry.z_top_m];
            let controls: [StrainMap; 3] = core::array::from_fn(|control| {
                let mut map = [[0.0; 8]; 5];
                for i in 0..3 {
                    for j in 0..3 {
                        map[i][j] = in_plane[i][j];
                        map[i][j + 3] = in_plane[i][j] * z[control];
                    }
                }
                for i in 0..2 {
                    for j in 0..2 {
                        map[i + 3][j + 6] = shear[control][2 * i + j];
                    }
                }
                map
            });
            if controls.iter().flatten().flatten().any(|x| !x.is_finite()) {
                return Err(invalid("original layer strain trial is unrepresentable"));
            }
            // Positive three-point Gauss integrates all products of the
            // quadratic Bernstein strain map exactly in real arithmetic.
            const ABSCISSA: f64 = 0.7745966692414834;
            let fractions = [0.5 * (1.0 - ABSCISSA), 0.5, 0.5 * (1.0 + ABSCISSA)];
            let weights = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0];
            let mut quadrature = [([[0.0; 8]; 5], 0.0); 3];
            for point in 0..3 {
                let weight = thickness * weights[point];
                if !weight.is_finite() || weight <= 0.0 {
                    return Err(invalid("positive original thickness weight is unrepresentable"));
                }
                quadrature[point] = (blend(&controls, fractions[point])?, weight);
            }
            layers.push(LayerTrial { elastic, controls, quadrature });
        }
        check(&mut cancelled)?;
        Ok(Self { section, shear_trial, layers })
    }

    /// Actual immutable stack, including every source receipt and reference face.
    pub fn section(&self) -> &LaminateSection { &self.section }

    /// Explicit original elastic shear profile used as the fixed trial field.
    pub fn shear_trial(&self) -> &LaminateShearRecovery { &self.shear_trial }

    /// Original actual [I0,I1,I2] in [kg/m2,kg/m,kg]; no homogeneous averaging.
    pub fn inertia_moments(&self) -> [f64; 3] {
        [self.section.i0_kg_m2(), self.section.i1_kg_m(), self.section.i2_kg()]
    }

    fn layer_law(&self, layer: usize) -> Result<&TensorProny> {
        let actual = self.section.layers().get(layer)
            .ok_or_else(|| invalid("original material layer index is outside the section"))?;
        match &actual.law {
            LayerLaw::Prony(law) => Ok(law),
            LayerLaw::Elastic(_) => self.layers.get(layer).and_then(|l| l.elastic.as_deref())
                .ok_or_else(|| invalid("original elastic layer interface is missing")),
            LayerLaw::PlaneStress { .. } => Err(invalid("unadmitted reduced-only source")),
        }
    }

    fn frequency_budget(&self, budget: HarmonicLaminateBudget, payload: u128) -> Result<()> {
        cap("layers", self.layers.len() as u128, budget.maximum_layers as u128)?;
        cap("frequency payload bytes", payload, budget.maximum_payload_bytes as u128)?;
        let mut work = 0_u128;
        for layer in 0..self.layers.len() {
            work += 500_000 + 4096 * self.layer_law(layer)?.branches().len() as u128;
        }
        cap("frequency work", work, budget.maximum_work as u128)
    }

    /// Prepare the full 8x8 harmonic section at one original source frequency.
    ///
    /// Generalized strain order is `[ex,ey,gxy,kx,ky,kxy,gxz,gyz]`, with
    /// physical in-plane strain `epsilon(z)=epsilon0+z*kappa`. Resultant order
    /// is `[Nx,Ny,Nxy,Mx,My,Mxy,Qs_x,Qs_y]`. The last two entries are
    /// GENERALIZED shear forces conjugate to the declared profile coordinates:
    /// `Qs = integral(P_s^T tau_actual dz)`, where `P_s` maps those coordinates
    /// to physical layer shear strain. They need not equal `integral(tau dz)`
    /// for an actual material different from the reference equilibrated source.
    /// The full complex matrix retains membrane/bending/shear cross blocks,
    /// including signed offset coupling. Existing cached laminate A/B/D/H are
    /// not substituted for this reduction.
    pub fn prepare_harmonic(
        &self,
        omega: f64,
        budget: HarmonicLaminateBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicLaminate<'_>> {
        check(&mut cancelled)?;
        let bytes = self.layers.len() as u128 * core::mem::size_of::<NormalStressHarmonic<'_>>() as u128;
        self.frequency_budget(budget, bytes)?;
        let mut layers = Vec::with_capacity(self.layers.len());
        let mut sums = [[Sum::default(); 64]; 2];
        for (index, trial) in self.layers.iter().enumerate() {
            check(&mut cancelled)?;
            let material = self.layer_law(index)?
                .prepare_normal_stress_harmonic(MaterialPlane::Axes12, omega)
                .map_err(|source| material_error(index, source))?;
            integrate(
                &mut sums, trial, [material.storage_pa(), material.loss_pa()], &mut cancelled,
            )?;
            layers.push(material);
        }
        let [storage, loss] = finish(sums)?;
        check(&mut cancelled)?;
        Ok(PreparedHarmonicLaminate { source: self, layers, storage, loss, omega })
    }
}

/// One frequency, sharing original sources, fixed profiles and normal constraints.
/// Repeated resultant/strain observations allocate no heap storage.
pub struct PreparedHarmonicLaminate<'a> {
    source: &'a HarmonicLaminate,
    layers: Vec<NormalStressHarmonic<'a>>,
    storage: [f64; 64],
    loss: [f64; 64],
    omega: f64,
}

impl<'a> PreparedHarmonicLaminate<'a> {
    pub const fn source(&self) -> &'a HarmonicLaminate { self.source }
    pub const fn omega_rad_s(&self) -> f64 { self.omega }

    /// Row-major complete storage operator. Units follow the work-conjugate
    /// strain/resultant coordinates: A/H [N/m], B [N], D [N m].
    pub const fn storage(&self) -> &[f64; 64] { &self.storage }

    /// Row-major complete loss operator, with passive sign. Physical stress
    /// resultants use `storage - i*loss` for exp(-i omega t).
    pub const fn loss(&self) -> &[f64; 64] { &self.loss }

    /// Original prepared five-component material in the named physical layer.
    pub fn layer_material(&self, layer: usize) -> Option<&NormalStressHarmonic<'a>> {
        self.layers.get(layer)
    }

    /// Admit the entire thickness and entire harmonic cycle at one in-plane
    /// location. The original three quadratic Bernstein controls per layer
    /// bound every physical strain by convexity; they are not sampled maxima.
    /// This sufficient bound may conservatively refuse a field within the cap.
    ///
    /// A native triangle consumer can apply this to each of its ORIGINAL
    /// in-plane Bernstein controls to certify the tensor-product spatial field.
    pub fn admit_generalized_strain(
        &self,
        strain: &[C64; 8],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<()> {
        check(&mut cancelled)?;
        self.admit_input(strain)?;
        for (index, (material, trial)) in self.layers.iter().zip(&self.source.layers).enumerate() {
            check(&mut cancelled)?;
            for map in &trial.controls {
                let retained = apply(map, strain)?;
                material.recover_strain(&retained).map_err(|source| material_error(index, source))?;
            }
        }
        check(&mut cancelled)
    }

    fn admit_input(&self, strain: &[C64; 8]) -> Result<()> {
        if strain.iter().any(|z| !z.re.is_finite() || !z.im.is_finite())
            || (self.omega == 0.0 && strain.iter().any(|z| z.im != 0.0))
        {
            Err(invalid("generalized harmonic strains must be finite and real at DC"))
        } else {
            Ok(())
        }
    }

    /// Full six-component Mandel strain in the ORIGINAL layer material frame
    /// at a physical fraction through that layer, including its actual normal
    /// relaxation and physical transverse shear. This query admits this point;
    /// use `admit_generalized_strain` for the complete thickness.
    pub fn recover_strain(
        &self,
        layer: usize,
        fraction: f64,
        generalized: &[C64; 8],
    ) -> Result<[C64; 6]> {
        self.admit_input(generalized)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(invalid("observation fraction must be inside the original layer"));
        }
        let material = self.layers.get(layer)
            .ok_or_else(|| invalid("original material layer index is outside the section"))?;
        let map = blend(&self.source.layers[layer].controls, fraction)?;
        material.recover_strain(&apply(&map, generalized)?)
            .map_err(|source| material_error(layer, source))
    }

    /// Complete work-conjugate resultants after full-thickness source admission.
    /// The last two entries are the fixed-trial projection `integral(P_s^T tau dz)`,
    /// not independently equilibrated actual shear-force integrals. Signed
    /// membrane/bending/shear interference remains in this same work operator.
    pub fn resultants(
        &self,
        generalized: &[C64; 8],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<[C64; 8]> {
        self.admit_generalized_strain(generalized, &mut cancelled)?;
        let mut result = [C64::ZERO; 8];
        for i in 0..8 {
            check(&mut cancelled)?;
            let (mut re, mut im) = (Sum::default(), Sum::default());
            for j in 0..8 {
                re.add([self.storage[8 * i + j], generalized[j].re, 1.0, 1.0]);
                re.add([self.loss[8 * i + j], generalized[j].im, 1.0, 1.0]);
                im.add([self.storage[8 * i + j], generalized[j].im, 1.0, 1.0]);
                im.add([-self.loss[8 * i + j], generalized[j].re, 1.0, 1.0]);
            }
            result[i] = C64::new(re.value(), im.value());
            if !result[i].re.is_finite() || !result[i].im.is_finite() {
                return Err(invalid("original harmonic resultants are unrepresentable"));
            }
        }
        check(&mut cancelled)?;
        Ok(result)
    }

    /// Exact storage/loss slopes per rad/s with the DECLARED profile held fixed.
    /// Every layer reuses its prepared normal constraint and original material
    /// tangent. Geometry, source rotations, physical mass and profile have zero
    /// derivative; no perturbed frequency or re-equilibrated profile is used.
    /// This bounded query allocates no heap storage.
    pub fn frequency_derivative(
        &self,
        budget: HarmonicLaminateBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<([f64; 64], [f64; 64])> {
        check(&mut cancelled)?;
        self.source.frequency_budget(budget, 0)?;
        let mut sums = [[Sum::default(); 64]; 2];
        for (index, (material, trial)) in self.layers.iter().zip(&self.source.layers).enumerate() {
            check(&mut cancelled)?;
            let (storage, loss) = material.frequency_derivative()
                .map_err(|source| material_error(index, source))?;
            integrate(&mut sums, trial, [&storage, &loss], &mut cancelled)?;
        }
        let [storage, loss] = finish(sums)?;
        check(&mut cancelled)?;
        Ok((storage, loss))
    }
}

fn blend(controls: &[StrainMap; 3], t: f64) -> Result<StrainMap> {
    let weights = [(1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t), t * t];
    let mut out = [[0.0; 8]; 5];
    for i in 0..5 {
        for j in 0..8 {
            let mut sum = Sum::default();
            for k in 0..3 {
                sum.add([controls[k][i][j], weights[k], 1.0, 1.0]);
            }
            out[i][j] = sum.value();
            if !out[i][j].is_finite() {
                return Err(invalid("original strain trial interpolation is unrepresentable"));
            }
        }
    }
    Ok(out)
}

fn apply(map: &StrainMap, strain: &[C64; 8]) -> Result<[C64; 5]> {
    let mut out = [C64::ZERO; 5];
    for i in 0..5 {
        let (mut re, mut im) = (Sum::default(), Sum::default());
        for j in 0..8 {
            re.add([map[i][j], strain[j].re, 1.0, 1.0]);
            im.add([map[i][j], strain[j].im, 1.0, 1.0]);
        }
        out[i] = C64::new(re.value(), im.value());
        if !out[i].re.is_finite() || !out[i].im.is_finite() {
            return Err(invalid("original layer strain is unrepresentable"));
        }
    }
    Ok(out)
}

fn integrate(
    sums: &mut [[Sum; 64]; 2],
    trial: &LayerTrial,
    matrices: [&[f64; 25]; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    for (map, weight) in &trial.quadrature {
        check(cancelled)?;
        for i in 0..8 {
            for j in 0..=i {
                for a in 0..5 {
                    for b in 0..5 {
                        for channel in 0..2 {
                            sums[channel][8 * i + j].add([
                                map[a][i], matrices[channel][5 * a + b], map[b][j], *weight,
                            ]);
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

fn finish(sums: [[Sum; 64]; 2]) -> Result<[[f64; 64]; 2]> {
    let mut out = [[0.0; 64]; 2];
    for channel in 0..2 {
        for i in 0..8 {
            for j in 0..=i {
                let value = sums[channel][8 * i + j].value();
                if !value.is_finite() {
                    return Err(invalid("original section work is unrepresentable"));
                }
                out[channel][8 * i + j] = value;
                out[channel][8 * j + i] = value;
            }
        }
    }
    Ok(out)
}
