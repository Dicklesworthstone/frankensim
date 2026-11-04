//! Reciprocal, three-dimensional, small-signal piezoelectric material.
//!
//! Source convention: `T=C^E S-e^T E`, `D=e S+epsilon^S E`. Mechanical
//! components use the orthonormal Mandel order `[xx,yy,zz,sqrt(2)xy,
//! sqrt(2)yz,sqrt(2)zx]`; electrical components use `[x,y,z]` in the SAME
//! explicitly named frame. `epsilon^S` is absolute fully clamped permittivity,
//! not relative epsilon or the stress-free `epsilon^T` of a datasheet.
//!
//! The independent strain/charge storage is
//! `W(S,D)=S.C^E.S/2+(D-eS).(epsilon^S)^-1.(D-eS)/2`.
//! Its derivatives supply both stress and field. Positive definite `C^E` and
//! `epsilon^S` make this complete potential positive definite for any finite
//! reciprocal e; requiring e itself to be positive would be incorrect.
//! These are the linear stress-charge equations in COMSOL's *Piezoelectricity*
//! theory: <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.043.html>.
//!
//! No diagonal, zero, symmetry mate, poling direction or material coefficient
//! is guessed. All 64 source roles and original query receipts are retained.
//! This is a fixed-condition reversible law. Consumers own small-signal strain
//! and field restrictions; hysteresis, depoling, dielectric/mechanical loss,
//! temperature work, finite deformation and measured accuracy are not inferred.

use fs_blake3::{ContentHash, DomainHasher};
use fs_la::factor::{Cholesky, cholesky};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

use super::{
    COMPLIANCE_DIMS, LongitudinalPiezoelectric, PiezoelectricError, STRAIN_CHARGE_DIMS,
    TensorLongitudinalReduction,
};
use crate::field_storage::ELECTRIC_PERMITTIVITY_DIMS;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use crate::tensor::{mandel_basis, mandel_rotation, matmul3, rotate_mandel};

/// Directional bulk acoustic modes with the electrical wave constraint explicit.
pub mod waves;

/// Stress-charge tensor dimensions [C/m²], including Mandel shear scaling.
pub const STRESS_CHARGE_DIMS: Dims = Dims([-2, 0, 1, 0, 1, 0]);
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
type Result<T> = core::result::Result<T, PiezoelectricError>;
fn invalid(message: &'static str) -> PiezoelectricError {
    PiezoelectricError::Invalid(message)
}

fn name(prefix: &str, row: usize, column: usize) -> String {
    format!("{prefix}_{}{}", row + 1, column + 1)
}
const STIFFNESS: &str = "piezoelectric_stiffness_e_mandel";
const COUPLING: &str = "piezoelectric_stress_charge_mandel";
const PERMITTIVITY: &str = "piezoelectric_permittivity_s";
const COMPLIANCE: &str = "piezoelectric_compliance_e_mandel";
const STRAIN_CHARGE: &str = "piezoelectric_strain_charge_mandel";
const PERMITTIVITY_T: &str = "piezoelectric_permittivity_t";

fn admit_roles(
    material: &ResolvedMaterialStatePoint,
    requirements: &[ScalarPropertyRequirement],
) -> Result<()> {
    for requirement in requirements {
        let property = material
            .property(requirement.name())
            .ok_or_else(|| invalid("missing explicit piezoelectric tensor coefficient"))?;
        if property.requirement().quantity() != QuantitySpec::dimensional(requirement.dims())
            || !property.is_static_at_fixed_external_conditions()
        {
            return Err(invalid(
                "piezoelectric tensor roles require their exact SI units and static source conditions",
            ));
        }
    }
    Ok(())
}
fn source_matrix<const R: usize, const C: usize>(
    material: &ResolvedMaterialStatePoint,
    prefix: &str,
) -> [[f64; C]; R] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            material
                .property(&name(prefix, i, j))
                .expect("all tensor source roles were checked")
                .value_si()
        })
    })
}

fn dot<const N: usize>(a: &[f64; N], b: &[f64; N]) -> f64 {
    a.iter().zip(b).fold(0.0, |sum, (a, b)| a.mul_add(*b, sum))
}
fn apply<const R: usize, const C: usize>(a: &[[f64; C]; R], x: &[f64; C]) -> [f64; R] {
    core::array::from_fn(|i| dot(&a[i], x))
}

/// Native factorization is prepared once after symmetric diagonal scaling.
/// Response calls use only fixed-size stack vectors and triangular solves.
#[derive(Clone, Debug)]
struct ScaledSpd<const N: usize> {
    scale: [f64; N],
    factor: Cholesky,
}
impl<const N: usize> ScaledSpd<N> {
    fn new(matrix: &[[f64; N]; N]) -> Result<Self> {
        for i in 0..N {
            for j in 0..N {
                if !matrix[i][j].is_finite() || matrix[i][j] != matrix[j][i] {
                    return Err(invalid(
                        "piezoelectric diagonal blocks must be finite and exactly symmetric",
                    ));
                }
            }
            if matrix[i][i] <= 0.0 {
                return Err(invalid(
                    "piezoelectric diagonal blocks must be positive definite",
                ));
            }
        }
        let scale = core::array::from_fn(|i| matrix[i][i].sqrt());
        let normalized: Vec<_> = (0..N)
            .flat_map(|i| (0..N).map(move |j| (matrix[i][j] / scale[i]) / scale[j]))
            .collect();
        if normalized.iter().any(|x| !x.is_finite()) {
            return Err(invalid("scaled piezoelectric block is unrepresentable"));
        }
        let factor = cholesky(&normalized, N)
            .map_err(|_| invalid("piezoelectric diagonal blocks must be positive definite"))?;
        Ok(Self { scale, factor })
    }
    fn solve(&self, rhs: [f64; N]) -> Result<[f64; N]> {
        let mut result = core::array::from_fn(|i| rhs[i] / self.scale[i]);
        self.factor.solve(&mut result);
        for (value, scale) in result.iter_mut().zip(self.scale) {
            *value /= scale;
        }
        if result.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "piezoelectric constitutive solve is unrepresentable",
            ));
        }
        Ok(result)
    }
    fn energy(&self, x: &[f64; N]) -> Result<f64> {
        let scaled: [f64; N] = core::array::from_fn(|i| self.scale[i] * x[i]);
        let mut energy = 0.0;
        for i in 0..N {
            let value = (i..N).fold(0.0, |sum, j| self.factor.l(j, i).mul_add(scaled[j], sum));
            energy = (0.5 * value).mul_add(value, energy);
        }
        if !energy.is_finite() || x.iter().any(|value| *value != 0.0) && energy == 0.0 {
            return Err(invalid("piezoelectric storage is unrepresentable"));
        }
        Ok(energy)
    }
}

/// Complete source-resolved stress-charge tensor, prepared for repeated queries.
/// The caller explicitly binds the named source coefficients to `source_frame`;
/// a frame transform does not relabel or replace any original material receipt.
#[derive(Clone, Debug)]
pub struct TensorPiezoelectric {
    material: ResolvedMaterialStatePoint,
    source_frame: ContentHash,
    frame: ContentHash,
    source_to_frame: [[f64; 3]; 3],
    stiffness: [[f64; 6]; 6],
    coupling: [[f64; 6]; 3],
    permittivity: [[f64; 3]; 3],
    open_stiffness: [[f64; 6]; 6],
    elastic_factor: ScaledSpd<6>,
    dielectric_factor: ScaledSpd<3>,
    density: f64,
    identity: ContentHash,
}

/// Work-conjugate response in the material's current frame. Energy is the
/// strain/charge potential W, even when a field-controlled query supplies E.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorPiezoelectricResponse {
    /// Tensile stress in Mandel coordinates [Pa].
    pub stress_mandel_pa: [f64; 6],
    /// Electric displacement in current-frame Cartesian coordinates [C/m²].
    pub displacement_c_m2: [f64; 3],
    /// Electric field in current-frame Cartesian coordinates [V/m].
    pub electric_field_v_m: [f64; 3],
    /// Positive reversible storage per reference volume [J/m³].
    pub stored_energy_j_m3: f64,
}

impl TensorPiezoelectric {
    /// Exact source schema: all C^E(6×6), e(3×6), epsilon^S(3×3), and density.
    /// Matrix suffixes are one-based. Mandel shear values must already include
    /// their sqrt(2) factors; engineering-Voigt or d-coefficients cannot be used
    /// under these names. Zero entries and symmetric mates remain explicit data.
    pub fn requirements() -> Result<Vec<ScalarPropertyRequirement>> {
        Self::source_requirements(false)
    }

    /// Complete strain-charge source schema: s^E(6x6) [Pa^-1], d(3x6)
    /// [C/N=m/V], absolute epsilon^T(3x3) [F/m], and density. Names are
    /// `piezoelectric_compliance_e_mandel_ij`,
    /// `piezoelectric_strain_charge_mandel_ij` and
    /// `piezoelectric_permittivity_t_ij`, with one-based suffixes.
    /// All mechanical indices use the same Mandel order as [`Self::requirements`].
    /// Converting engineering data requires dividing each d shear column and
    /// each compliance shear index by sqrt(2); a shear/shear entry divides by 2.
    /// Neither a convention nor omitted zero entries is inferred.
    pub fn strain_charge_requirements() -> Result<Vec<ScalarPropertyRequirement>> {
        Self::source_requirements(true)
    }

    fn source_requirements(strain_charge: bool) -> Result<Vec<ScalarPropertyRequirement>> {
        let mut requirements = Vec::with_capacity(64);
        let roles = if strain_charge {
            [
                (COMPLIANCE, 6, 6, COMPLIANCE_DIMS),
                (STRAIN_CHARGE, 3, 6, STRAIN_CHARGE_DIMS),
                (PERMITTIVITY_T, 3, 3, ELECTRIC_PERMITTIVITY_DIMS),
            ]
        } else {
            [
                (STIFFNESS, 6, 6, Pressure::DIMS),
                (COUPLING, 3, 6, STRESS_CHARGE_DIMS),
                (PERMITTIVITY, 3, 3, ELECTRIC_PERMITTIVITY_DIMS),
            ]
        };
        for (prefix, rows, columns, dims) in roles {
            for i in 0..rows {
                for j in 0..columns {
                    requirements.push(ScalarPropertyRequirement::try_new(
                        name(prefix, i, j),
                        dims,
                        ScalarAdmissibility::Finite,
                    )?);
                }
            }
        }
        requirements.push(ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?);
        Ok(requirements)
    }

    /// Resolve every coefficient from one immutable card at one physical query.
    /// `source_frame` is the explicit shared Cartesian frame declared for the
    /// source schema; it must not be zero. No off-axis response is synthesized
    /// from the separate one-dimensional longitudinal property family.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        source_frame: ContentHash,
    ) -> Result<Self> {
        Self::from_material_state(
            resolve_material_state_point(card, point, &Self::requirements()?, selection)?,
            source_frame,
        )
    }

    /// Retain an already-resolved complete source bundle. Frequency, strain,
    /// stress and field dependencies AND applicability limits may not be frozen
    /// into this constant-coefficient model. Fixed external conditions survive.
    pub fn from_material_state(
        material: ResolvedMaterialStatePoint,
        source_frame: ContentHash,
    ) -> Result<Self> {
        admit_roles(&material, &Self::requirements()?)?;
        let stiffness = source_matrix(&material, STIFFNESS);
        let coupling = source_matrix(&material, COUPLING);
        let permittivity = source_matrix(&material, PERMITTIVITY);
        let density = material
            .property(DENSITY_PROPERTY)
            .expect("density role was checked")
            .value_si();
        Self::prepare(
            material,
            source_frame,
            source_frame,
            IDENTITY,
            stiffness,
            coupling,
            permittivity,
            density,
        )
    }

    /// Resolve original datasheet-form strain-charge coefficients at one common
    /// source query. The complete coupled Schur complement must be stable.
    pub fn resolve_strain_charge(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        source_frame: ContentHash,
    ) -> Result<Self> {
        Self::from_strain_charge_material_state(
            resolve_material_state_point(
                card,
                point,
                &Self::strain_charge_requirements()?,
                selection,
            )?,
            source_frame,
        )
    }

    /// Convert `S=s^E T+d^T E`, `D=d T+epsilon^T E` exactly once, retaining
    /// the original source bundle. Native scaled solves produce C=(s^E)^-1,
    /// e=d C and epsilon^S=epsilon^T-d C d^T. Positive separate compliance and
    /// stress-free permittivity alone do not establish coupled stability.
    pub fn from_strain_charge_material_state(
        material: ResolvedMaterialStatePoint,
        source_frame: ContentHash,
    ) -> Result<Self> {
        admit_roles(&material, &Self::strain_charge_requirements()?)?;
        let compliance: [[f64; 6]; 6] = source_matrix(&material, COMPLIANCE);
        let d: [[f64; 6]; 3] = source_matrix(&material, STRAIN_CHARGE);
        let epsilon_t: [[f64; 3]; 3] = source_matrix(&material, PERMITTIVITY_T);
        let inverse = ScaledSpd::new(&compliance)?;
        ScaledSpd::new(&epsilon_t)?;
        let mut stiffness = [[0.0; 6]; 6];
        for j in 0..6 {
            let column = inverse.solve(core::array::from_fn(|i| if i == j { 1.0 } else { 0.0 }))?;
            for i in 0..=j {
                stiffness[i][j] = column[i];
                stiffness[j][i] = column[i];
            }
        }
        let coupling = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..6).fold(0.0, |sum, k| d[i][k].mul_add(stiffness[k][j], sum))
            })
        });
        let mut permittivity = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let value = (0..6).fold(epsilon_t[i][j], |sum, k| {
                    (-coupling[i][k]).mul_add(d[j][k], sum)
                });
                permittivity[i][j] = value;
                permittivity[j][i] = value;
            }
        }
        let density = material
            .property(DENSITY_PROPERTY)
            .expect("density role was checked")
            .value_si();
        Self::prepare(
            material,
            source_frame,
            source_frame,
            IDENTITY,
            stiffness,
            coupling,
            permittivity,
            density,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        material: ResolvedMaterialStatePoint,
        source_frame: ContentHash,
        frame: ContentHash,
        source_to_frame: [[f64; 3]; 3],
        stiffness: [[f64; 6]; 6],
        coupling: [[f64; 6]; 3],
        permittivity: [[f64; 3]; 3],
        density: f64,
    ) -> Result<Self> {
        if source_frame == ContentHash([0; 32])
            || frame == ContentHash([0; 32])
            || !density.is_finite()
            || density <= 0.0
            || coupling.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "named frames, finite piezoelectric coupling and positive reference density are required",
            ));
        }
        let elastic_factor = ScaledSpd::new(&stiffness)?;
        let dielectric_factor = ScaledSpd::new(&permittivity)?;
        let mut columns = [[0.0; 3]; 6];
        for (j, column) in columns.iter_mut().enumerate() {
            *column = dielectric_factor.solve(core::array::from_fn(|i| coupling[i][j]))?;
        }
        let mut open_stiffness = stiffness;
        for i in 0..6 {
            for j in i..6 {
                let value = (0..3).fold(stiffness[i][j], |sum, k| {
                    coupling[k][i].mul_add(columns[j][k], sum)
                });
                if !value.is_finite() {
                    return Err(invalid(
                        "open-circuit piezoelectric stiffness is unrepresentable",
                    ));
                }
                open_stiffness[i][j] = value;
                open_stiffness[j][i] = value;
            }
        }
        // A failed numerical Schur complement must not become a mechanical law.
        ScaledSpd::new(&open_stiffness)?;
        let mut hash = DomainHasher::new("org.frankensim.fs-material.piezoelectric-tensor.v1");
        hash.update(material.identity().as_bytes());
        hash.update(source_frame.as_bytes());
        hash.update(frame.as_bytes());
        for value in source_to_frame
            .iter()
            .flatten()
            .chain(stiffness.iter().flatten())
            .chain(coupling.iter().flatten())
            .chain(permittivity.iter().flatten())
        {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.update(&density.to_bits().to_le_bytes());
        Ok(Self {
            material,
            source_frame,
            frame,
            source_to_frame,
            stiffness,
            coupling,
            permittivity,
            open_stiffness,
            elastic_factor,
            dielectric_factor,
            density,
            identity: hash.finalize(),
        })
    }

    /// Complete original same-card/source/query bundle and all selected receipts.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    /// Original common frame of every source coefficient.
    pub const fn source_frame_identity(&self) -> ContentHash {
        self.source_frame
    }
    /// Current frame in which all fields and tensors must be expressed.
    pub const fn frame_identity(&self) -> ContentHash {
        self.frame
    }
    /// Proper rotation carrying original source components into the current frame.
    pub const fn source_to_frame_rotation(&self) -> &[[f64; 3]; 3] {
        &self.source_to_frame
    }
    /// Identity binds the original source, frame mapping and actual tensor values.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    /// Fixed-field stiffness C^E in Mandel coordinates [Pa].
    pub const fn stiffness_mandel_pa(&self) -> &[[f64; 6]; 6] {
        &self.stiffness
    }
    /// Reciprocal stress/charge matrix e in Cartesian-by-Mandel coordinates [C/m²].
    pub const fn coupling_c_m2(&self) -> &[[f64; 6]; 3] {
        &self.coupling
    }
    /// Absolute fully clamped permittivity epsilon^S [F/m].
    pub const fn permittivity_f_m(&self) -> &[[f64; 3]; 3] {
        &self.permittivity
    }
    /// Local fixed-D stiffness C^D=C^E+e^T(epsilon^S)^-1 e [Pa]. Boundary
    /// electrodes in a nonuniform body still require the electrostatic solve.
    pub const fn open_circuit_stiffness_mandel_pa(&self) -> &[[f64; 6]; 6] {
        &self.open_stiffness
    }
    /// Original fixed-reference density at the same material state [kg/m³].
    pub const fn density_kg_m3(&self) -> f64 {
        self.density
    }

    /// Reduce the SAME tensor to a homogeneous longitudinal sensor/actuator law.
    /// With unit direction n in the current frame, stress is `T=t n⊗n` and
    /// electric field is `E=E_n n`: transverse stresses and transverse E vanish,
    /// while transverse strain and transverse D may relax. These are explicit
    /// constitutive constraints, not a finite specimen's lateral boundary solve.
    /// In particular, they are not the charge-free bulk-wave condition n.D=0.
    ///
    /// Writing r=mandel(n⊗n), two prepared elastic solves supply
    /// `s=r.C^-1.r`, `d=n.e.C^-1.r`, and
    /// `epsilon_T=n.(epsilon_S+e.C^-1.e^T).n`. Taking 1/C33 or e33/C33
    /// would miss transverse and shear relaxation in an anisotropic solid.
    /// The original source bundle is retained unchanged; direction/frame and
    /// tensor identity are recorded separately rather than inventing claims.
    /// Geometry, small-signal envelopes, lateral electrical constraints, finite
    /// wavelength effects and thermal evolution remain the consumer's concern.
    pub fn longitudinal(&self, direction: [f64; 3]) -> Result<LongitudinalPiezoelectric> {
        let scale = direction.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        if direction.iter().any(|x| !x.is_finite()) || scale == 0.0 {
            return Err(invalid(
                "longitudinal tensor direction must be finite and nonzero",
            ));
        }
        let scaled = direction.map(|x| x / scale);
        let length = fs_math::det::sqrt(dot(&scaled, &scaled));
        let n = scaled.map(|x| if x == 0.0 { 0.0 } else { x / length });
        let basis = mandel_basis();
        let r: [f64; 6] = core::array::from_fn(|a| {
            (0..3).fold(0.0, |sum, i| {
                (0..3).fold(sum, |sum, j| (basis[a][i][j] * n[i]).mul_add(n[j], sum))
            })
        });
        let electric_stress: [f64; 6] = core::array::from_fn(|a| {
            (0..3).fold(0.0, |sum, i| self.coupling[i][a].mul_add(n[i], sum))
        });
        let mechanical_strain = self.elastic_factor.solve(r)?;
        let electric_strain = self.elastic_factor.solve(electric_stress)?;
        let compliance_e = dot(&r, &mechanical_strain);
        let d = dot(&electric_stress, &mechanical_strain);
        let epsilon_t =
            dot(&n, &apply(&self.permittivity, &n)) + dot(&electric_stress, &electric_strain);
        LongitudinalPiezoelectric::prepare(
            self.material.clone(),
            compliance_e,
            d,
            epsilon_t,
            self.density,
            Some(TensorLongitudinalReduction {
                tensor_identity: self.identity,
                frame_identity: self.frame,
                direction: n,
            }),
        )
    }

    /// Express all three blocks in one new, explicitly named frame. `rotation`
    /// maps CURRENT components to target components. Strain/stress use the
    /// common Mandel rotation M: C'=M C M^T, e'=Q e M^T, epsilon'=Q epsilon Q^T.
    /// Reflections, nonorthogonal matrices and same-frame relabelling refuse.
    pub fn rotated(&self, target_frame: ContentHash, rotation: [[f64; 3]; 3]) -> Result<Self> {
        if target_frame == ContentHash([0; 32])
            || rotation.iter().flatten().any(|x| !x.is_finite())
            || target_frame == self.frame && rotation != IDENTITY
        {
            return Err(invalid(
                "a nonzero target frame and consistent proper rotation are required",
            ));
        }
        for i in 0..3 {
            for j in 0..3 {
                if (dot(&rotation[i], &rotation[j]) - if i == j { 1.0 } else { 0.0 }).abs() > 2e-12
                {
                    return Err(invalid("piezoelectric frame map must be orthonormal"));
                }
            }
        }
        let q = rotation;
        let determinant = q[0][0] * (q[1][1] * q[2][2] - q[1][2] * q[2][1])
            - q[0][1] * (q[1][0] * q[2][2] - q[1][2] * q[2][0])
            + q[0][2] * (q[1][0] * q[2][1] - q[1][1] * q[2][0]);
        if (determinant - 1.0).abs() > 4e-12 {
            return Err(invalid("piezoelectric frame map must preserve handedness"));
        }
        let m = mandel_rotation(rotation);
        let coupling = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let mut value = 0.0;
                for p in 0..3 {
                    for k in 0..6 {
                        value = (rotation[i][p] * self.coupling[p][k]).mul_add(m[j][k], value);
                    }
                }
                value
            })
        });
        let mut permittivity = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let mut value = 0.0;
                for p in 0..3 {
                    for k in 0..3 {
                        value = (rotation[i][p] * self.permittivity[p][k])
                            .mul_add(rotation[j][k], value);
                    }
                }
                // Source symmetry was admitted before rotation; preserve it exactly.
                permittivity[i][j] = value;
                permittivity[j][i] = value;
            }
        }
        Self::prepare(
            self.material.clone(),
            self.source_frame,
            target_frame,
            matmul3(&rotation, &self.source_to_frame),
            rotate_mandel(self.stiffness, rotation),
            coupling,
            permittivity,
            self.density,
        )
    }

    /// Evaluate the stress-charge equations at specified strain and field. The
    /// reported storage is W(S,D(S,E)), NOT electric enthalpy; its derivative
    /// with E held fixed is therefore not asserted to equal stress.
    pub fn field_response(
        &self,
        strain_mandel: [f64; 6],
        electric_field_v_m: [f64; 3],
    ) -> Result<TensorPiezoelectricResponse> {
        let polarization = apply(&self.coupling, &strain_mandel);
        let dielectric = apply(&self.permittivity, &electric_field_v_m);
        self.response(
            strain_mandel,
            electric_field_v_m,
            core::array::from_fn(|i| polarization[i] + dielectric[i]),
        )
    }

    /// Evaluate W(S,D) and its exact constitutive gradients (T,E). Independent
    /// coordinates are strain and electric displacement in the current frame.
    pub fn charge_response(
        &self,
        strain_mandel: [f64; 6],
        displacement_c_m2: [f64; 3],
    ) -> Result<TensorPiezoelectricResponse> {
        let polarization = apply(&self.coupling, &strain_mandel);
        let field = self.dielectric_factor.solve(core::array::from_fn(|i| {
            displacement_c_m2[i] - polarization[i]
        }))?;
        self.response(strain_mandel, field, displacement_c_m2)
    }

    fn response(
        &self,
        strain: [f64; 6],
        field: [f64; 3],
        displacement: [f64; 3],
    ) -> Result<TensorPiezoelectricResponse> {
        if strain
            .iter()
            .chain(&field)
            .chain(&displacement)
            .any(|x| !x.is_finite())
        {
            return Err(invalid("piezoelectric coordinates must be finite"));
        }
        let elastic = apply(&self.stiffness, &strain);
        let stress = core::array::from_fn(|i| {
            (0..3).fold(elastic[i], |sum, j| {
                (-self.coupling[j][i]).mul_add(field[j], sum)
            })
        });
        let stored_energy_j_m3 =
            self.elastic_factor.energy(&strain)? + self.dielectric_factor.energy(&field)?;
        if !stored_energy_j_m3.is_finite() || stress.iter().any(|x| !x.is_finite()) {
            return Err(invalid("piezoelectric tensor response is unrepresentable"));
        }
        Ok(TensorPiezoelectricResponse {
            stress_mandel_pa: stress,
            displacement_c_m2: displacement,
            electric_field_v_m: field,
            stored_energy_j_m3,
        })
    }
}
