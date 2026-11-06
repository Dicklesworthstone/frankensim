//! Pointwise Saint-Venant torsion in a circular, perfectly bonded annulus.
//! The axial direction is local +x; the sampled coordinates are local [y, z].

/// Elastic torsional contribution at a material point, in the section frame.
/// This excludes axial/bending normal stress and transverse-force shear.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TorsionalStressPoint {
    /// Physical [tau_xy, tau_xz] on the cross-section normal to local +x [Pa].
    pub shear_stress_pa: [f64; 2],
    /// Engineering [gamma_xy, gamma_xz], twice the corresponding tensor strain.
    pub engineering_shear_strain: [f64; 2],
    /// Reversible elastic torsional storage per reference volume [J/m3].
    pub stored_energy_density_j_m3: f64,
}

impl TorsionalStressPoint {
    /// Tensor-Voigt [xx, yy, zz, xy, yz, zx] torsional Cauchy stress [Pa].
    /// Stress shear is physical shear, with no engineering factor of two.
    #[must_use]
    pub fn stress_tensor_pa(&self) -> [f64; 6] {
        [
            0.0,
            0.0,
            0.0,
            self.shear_stress_pa[0],
            0.0,
            self.shear_stress_pa[1],
        ]
    }

    /// The same torsional strain in tensor-Voigt [xx, yy, zz, xy, yz, zx].
    #[must_use]
    pub fn strain_tensor(&self) -> [f64; 6] {
        [
            0.0,
            0.0,
            0.0,
            0.5 * self.engineering_shear_strain[0],
            0.0,
            0.5 * self.engineering_shear_strain[1],
        ]
    }
}

// Source selection belongs to ConcentricShearSection. Keeping the point law
// separate makes its force/energy identities independently testable without a
// second implementation of the constitutive response.
pub(crate) fn sample(
    inner_radius_m: f64,
    outer_radius_m: f64,
    shear_modulus_pa: f64,
    position_yz_m: [f64; 2],
    twist_per_length_rad_m: f64,
) -> Result<TorsionalStressPoint, &'static str> {
    let [y, z] = position_yz_m;
    if ![
        inner_radius_m,
        outer_radius_m,
        shear_modulus_pa,
        y,
        z,
        twist_per_length_rad_m,
    ]
    .into_iter()
    .all(f64::is_finite)
        || inner_radius_m < 0.0
        || outer_radius_m <= inner_radius_m
        || shear_modulus_pa <= 0.0
    {
        return Err(
            "torsional stress needs a finite annulus, positive shear modulus and finite point/twist",
        );
    }
    let radius = fs_math::det::hypot(y, z);
    if radius < inner_radius_m || radius > outer_radius_m {
        return Err("torsional stress point lies outside the selected material layer");
    }
    let engineering_shear_strain = [-twist_per_length_rad_m * z, twist_per_length_rad_m * y];
    let shear_stress_pa = engineering_shear_strain.map(|strain| shear_modulus_pa * strain);
    let stored_energy_density_j_m3 = shear_stress_pa
        .into_iter()
        .zip(engineering_shear_strain)
        .map(|(stress, strain)| 0.5 * stress * strain)
        .sum();
    if !engineering_shear_strain
        .into_iter()
        .chain(shear_stress_pa)
        .chain([stored_energy_density_j_m3])
        .all(f64::is_finite)
        || engineering_shear_strain
            .into_iter()
            .zip([z, y])
            .any(|(strain, coordinate)| twist_per_length_rad_m != 0.0 && coordinate != 0.0 && strain == 0.0)
        || shear_stress_pa
            .into_iter()
            .zip(engineering_shear_strain)
            .any(|(stress, strain)| strain != 0.0 && stress == 0.0)
        || stored_energy_density_j_m3 == 0.0
            && engineering_shear_strain.into_iter().any(|strain| strain != 0.0)
    {
        return Err("torsional point strain, stress or storage is unrepresentable");
    }
    Ok(TorsionalStressPoint {
        shear_stress_pa,
        engineering_shear_strain,
        stored_energy_density_j_m3,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: f64, b: f64) {
        assert!((a - b).abs() <= 2e-14 * a.abs().max(b.abs()), "{a} != {b}");
    }

    #[test]
    fn g1_signed_shear_tensor_and_work_are_work_conjugate() {
        let p = sample(0.0, 0.01, 80e9, [0.006, 0.008], 0.02).unwrap();
        assert_eq!(p.engineering_shear_strain, [-0.00016, 0.00012]);
        close(p.shear_stress_pa[0], -12.8e6);
        close(p.shear_stress_pa[1], 9.6e6);
        close(p.shear_stress_pa[0].hypot(p.shear_stress_pa[1]), 16e6);
        let stress = p.stress_tensor_pa();
        let strain = p.strain_tensor();
        close(
            p.stored_energy_density_j_m3,
            stress[3] * strain[3] + stress[5] * strain[5],
        );
        let reverse = sample(0.0, 0.01, 80e9, [0.006, 0.008], -0.02).unwrap();
        assert_eq!(reverse.shear_stress_pa, p.shear_stress_pa.map(|x| -x));
        assert_eq!(
            reverse.stored_energy_density_j_m3,
            p.stored_energy_density_j_m3
        );
    }

    #[test]
    fn g1_layered_field_integrates_to_original_torque_and_energy() {
        let twist = 0.07;
        let mut torque = 0.0;
        let mut storage = 0.0;
        let mut gj = 0.0;
        for (inner, outer, modulus) in [(0.003_f64, 0.01_f64, 40e9), (0.01, 0.02, 80e9)] {
            // One midpoint in r^2 integrates these radial quadratic fields exactly.
            let radius = (0.5 * (inner * inner + outer * outer)).sqrt();
            let area = core::f64::consts::PI * (outer * outer - inner * inner);
            let p = sample(inner, outer, modulus, [radius, 0.0], twist).unwrap();
            torque += radius * p.shear_stress_pa[1] * area;
            storage += p.stored_energy_density_j_m3 * area;
            gj += modulus * core::f64::consts::FRAC_PI_2 * (outer.powi(4) - inner.powi(4));
        }
        close(torque, gj * twist);
        close(storage, 0.5 * gj * twist * twist);
        let a = sample(0.003, 0.01, 40e9, [0.01, 0.0], twist).unwrap();
        let b = sample(0.01, 0.02, 80e9, [0.01, 0.0], twist).unwrap();
        assert_eq!(a.engineering_shear_strain, b.engineering_shear_strain);
        assert_eq!(b.shear_stress_pa[1], 2.0 * a.shear_stress_pa[1]);
    }

    #[test]
    fn g3_material_subdivision_preserves_the_point_field() {
        for point in [[0.01, 0.0], [0.0, -0.015]] {
            assert_eq!(
                sample(0.0, 0.02, 40e9, point, 0.07).unwrap(),
                sample(0.01, 0.02, 40e9, point, 0.07).unwrap()
            );
        }
    }

    #[test]
    fn g0_bore_exterior_nonfinite_and_overflow_refuse() {
        for point in [[0.0; 2], [0.03, 0.0], [f64::NAN, 0.0], [f64::INFINITY; 2]] {
            assert!(sample(0.01, 0.02, 40e9, point, 0.07).is_err());
        }
        assert!(sample(0.0, 0.02, 40e9, [0.01, 0.0], f64::NAN).is_err());
        assert!(sample(0.0, 0.02, f64::MAX, [0.01, 0.0], 1000.0).is_err());
        let zero = sample(0.0, 0.02, 40e9, [0.01, 0.0], 0.0).unwrap();
        assert_eq!(zero.shear_stress_pa, [0.0; 2]);
        assert_eq!(zero.stored_energy_density_j_m3, 0.0);
    }
}
