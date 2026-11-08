//! Physical section inventory about the original beam reference axis.

use super::{BeamMass, BeamSection, PlateError, bad};

/// Six raw density moments over the original cross-section in local `(y,z)`.
///
/// These are `integral rho*[1,y,z,y²,y*z,z²] dA`, with units respectively
/// kg/m, kg, kg, kg·m, kg·m, kg·m. The elastic reference axis and section
/// orientation are retained; the mass centroid need not lie on that axis.
/// A positive moment matrix is required. Admission does not establish that
/// the inventory fits inside a geometric section whose boundary is not given.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamSectionInertia {
    moments: [f64; 6],
    // Retain the scalar source for exact dispatch when its geometric moments
    // still match. It never replaces or changes the six physical moments.
    uniform_density_kg_m3: Option<f64>,
}

impl BeamSectionInertia {
    /// Admit raw reference moments, preserving signed first/product moments.
    /// # Errors
    /// Nonfinite entries, nonpositive mass, or a nonpositive moment matrix.
    /// No conditioning threshold, diagonal floor or centroid repair is applied.
    pub fn new(
        line_mass_kg_per_m: f64,
        first_y_kg: f64,
        first_z_kg: f64,
        second_yy_kg_m: f64,
        second_yz_kg_m: f64,
        second_zz_kg_m: f64,
    ) -> Result<Self, PlateError> {
        let moments = [
            line_mass_kg_per_m,
            first_y_kg,
            first_z_kg,
            second_yy_kg_m,
            second_yz_kg_m,
            second_zz_kg_m,
        ];
        let [a, b, c, d, e, f] = moments;
        if moments.iter().any(|v| !v.is_finite())
            || a <= 0.0
            || d <= 0.0
            || f <= 0.0
            || !positive_sum([([a, d, 1.0], true, 0), ([b, b, 1.0], false, 0)])
            || !positive_sum([
                ([a, d, f], true, 0),
                ([b, c, e], true, 1),
                ([a, e, e], false, 0),
                ([d, c, c], false, 0),
                ([f, b, b], false, 0),
            ])
        {
            return Err(bad(
                "beam section density moments must have a positive-definite physical moment matrix",
            ));
        }
        Ok(Self {
            moments,
            uniform_density_kg_m3: None,
        })
    }

    /// Retain the original centered, uniform cross-sectional inventory.
    /// `Iy` is the geometric z² integral and `Iz` the y² integral; the torsional
    /// stiffness constant `J` is never a mass moment.
    /// # Errors
    /// Invalid original section or unrepresentable physical moments.
    pub fn from_uniform(section: &BeamSection) -> Result<Self, PlateError> {
        section.validate()?;
        let rho = section.density_kg_m3;
        let mut out = Self::new(
            rho * section.area_m2,
            0.0,
            0.0,
            rho * section.iz_m4,
            0.0,
            rho * section.iy_m4,
        )?;
        out.uniform_density_kg_m3 = Some(rho);
        Ok(out)
    }

    /// `integral rho dA` [kg/m].
    #[must_use]
    pub const fn line_mass_kg_per_m(self) -> f64 {
        self.moments[0]
    }
    /// `integral rho*y dA` [kg], with signed local y.
    #[must_use]
    pub const fn first_y_kg(self) -> f64 {
        self.moments[1]
    }
    /// `integral rho*z dA` [kg], with signed local z.
    #[must_use]
    pub const fn first_z_kg(self) -> f64 {
        self.moments[2]
    }
    /// `integral rho*y² dA` [kg·m].
    #[must_use]
    pub const fn second_yy_kg_m(self) -> f64 {
        self.moments[3]
    }
    /// `integral rho*y*z dA` [kg·m], preserving its sign.
    #[must_use]
    pub const fn second_yz_kg_m(self) -> f64 {
        self.moments[4]
    }
    /// `integral rho*z² dA` [kg·m].
    #[must_use]
    pub const fn second_zz_kg_m(self) -> f64 {
        self.moments[5]
    }

    fn uniform_density_for(self, section: BeamSection) -> Option<f64> {
        let rho = self.uniform_density_kg_m3.unwrap_or(section.density_kg_m3);
        (self.moments
            == [
                rho * section.area_m2,
                0.0,
                0.0,
                rho * section.iz_m4,
                0.0,
                rho * section.iy_m4,
            ])
        .then_some(rho)
    }
}

// Exact signs of sums of finite binary64 triple products. Significands need
// at most 159 bits; aligning their exponents in 100 u64 words covers the full
// binary64 range, including subnormals and the factor two in a determinant.
// Thus disparate physical units/scales cannot overflow, underflow or turn a
// singular inventory into a positive one through a rounded square root.
fn positive_sum<const N: usize>(terms: [([f64; 3], bool, usize); N]) -> bool {
    let mut positive = [0_u64; 100];
    let mut negative = [0_u64; 100];
    for (values, mut plus, doubling) in terms {
        if values.contains(&0.0) {
            continue;
        }
        let mut shift = 3222 + doubling;
        let mut significands = [0_u64; 3];
        for (i, value) in values.into_iter().enumerate() {
            let bits = value.to_bits();
            plus ^= bits >> 63 != 0;
            let exponent = ((bits >> 52) & 0x7ff) as usize;
            significands[i] = bits & ((1_u64 << 52) - 1);
            if exponent == 0 {
                shift -= 1074;
            } else {
                significands[i] |= 1_u64 << 52;
                shift = shift + exponent - 1075;
            }
        }
        let ab = u128::from(significands[0]) * u128::from(significands[1]);
        let low = u128::from(ab as u64) * u128::from(significands[2]);
        let high = (ab >> 64) * u128::from(significands[2]) + (low >> 64);
        let words = [low as u64, high as u64, (high >> 64) as u64];
        let sum = if plus { &mut positive } else { &mut negative };
        let offset = shift / 64;
        let bits = shift % 64;
        for (i, word) in words.into_iter().enumerate() {
            add_word(sum, offset + i, word << bits);
            if bits != 0 {
                add_word(sum, offset + i + 1, word >> (64 - bits));
            }
        }
    }
    positive.iter().rev().cmp(negative.iter().rev()).is_gt()
}

fn add_word(sum: &mut [u64; 100], mut i: usize, mut value: u64) {
    while value != 0 {
        let (next, carry) = sum[i].overflowing_add(value);
        sum[i] = next;
        value = u64::from(carry);
        i += 1;
    }
}

/// Full `integral rho*(v + omega cross r)^2 dA`, before axial quadrature.
/// `r=[0,y,z]`; the six rows are conjugate to `[vx,vy,vz,wx,wy,wz]`.
fn action(q: [f64; 6], u: [f64; 6]) -> [f64; 6] {
    let [m, y, z, yy, yz, zz] = q;
    let [vx, vy, vz, wx, wy, wz] = u;
    [
        m * vx + z * wy - y * wz,
        m * vy - z * wx,
        m * vz + y * wx,
        -z * vy + y * vz + (yy + zz) * wx,
        z * vx + zz * wy - yz * wz,
        -y * vx - yz * wy + yy * wz,
    ]
}

pub(super) fn local_mass(
    section: BeamSection,
    length: f64,
    phi: [f64; 2],
    mass: BeamMass,
    endpoints: [BeamSectionInertia; 2],
) -> Result<[f64; 144], PlateError> {
    if mass != BeamMass::Consistent {
        return Err(bad(
            "full beam section moments require consistent inertia; diagonal lumping loses the signed moment couplings",
        ));
    }
    if let (Some(a), Some(b)) = (
        endpoints[0].uniform_density_for(section),
        endpoints[1].uniform_density_for(section),
    ) {
        return super::local_mass_with_reference_density(section, length, phi, mass, [a, b]);
    }
    let mut out = [0.0; 144];
    if phi.iter().any(|p| !p.is_finite() || *p < 0.0) {
        return Err(bad("unresolved beam shear-to-bending ratio"));
    }
    for (t, weight) in super::GAUSS {
        let q = core::array::from_fn(|k| {
            length * weight * ((1.0 - t) * endpoints[0].moments[k] + t * endpoints[1].moments[k])
        });
        let mut fields = [[0.0; 6]; 12];
        for (dof, value) in [(0, 1.0 - t), (6, t), (3, 1.0 - t), (9, t)] {
            fields[dof][dof % 6] = value;
        }
        for (slots, signs, v, w, rotation_sign, p) in [
            ([1, 5, 7, 11], [1.0; 4], 1, 5, 1.0, phi[0]),
            ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0], 2, 4, -1.0, phi[1]),
        ] {
            let (n, h) = super::bending_fields(length, p, t);
            for i in 0..4 {
                fields[slots[i]][v] = signs[i] * n[i];
                fields[slots[i]][w] = rotation_sign * signs[i] * h[i];
            }
        }
        let weighted = fields.map(|u| action(q, u));
        for i in 0..12 {
            for j in i..12 {
                let value = fields[i]
                    .iter()
                    .zip(weighted[j])
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                out[12 * i + j] += value;
                if i != j {
                    out[12 * j + i] += value;
                }
            }
        }
    }
    if out.iter().any(|v| !v.is_finite()) || (0..12).any(|i| out[13 * i] <= 0.0) {
        return Err(bad("beam section inertia is not representable"));
    }
    Ok(out)
}
