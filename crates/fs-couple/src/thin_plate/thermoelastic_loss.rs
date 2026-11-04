//! Diagonal weak-loss thermoelastic image on the original DKT modes.
//!
//! Thermal diffusion is driven by mean curvature, so each quadrature term is
//! a positive coefficient times `(kappa_xx+kappa_yy)^2`. Conservative membrane
//! prestress belongs only in the total modal energy denominator. No elastic
//! matrix surrogate, new eigensolve, or damping of twisting strain is used.
//! This retains the thin-plate adapter's diagonal modal approximation; it does
//! not claim a transient thermal history or strongly coupled complex modes.

use super::{AcousticRealizeError, PlateMesh, PlateModel, ThermoelasticZener, map_plate};
use fs_material::visco::slab_ted::slab_relaxation_shape;

struct ThermalPoint {
    dofs: [usize; 9],
    trace: [f64; 9],
    weight: f64,
}

pub(super) struct ThermoelasticBendingProjection {
    first_relaxation_time_s: f64,
    points: Vec<ThermalPoint>,
}

impl ThermoelasticBendingProjection {
    pub(super) fn new(
        mesh: &PlateMesh,
        law: ThermoelasticZener,
        thickness: f64,
        poisson: f64,
    ) -> Result<Self, AcousticRealizeError> {
        let relaxation = law
            .plate_bending_relaxation_matrix(thickness, poisson)
            .map_err(|_| AcousticRealizeError::InvalidDescription {
                what: "plate thermoelastic bending tensor is invalid or unrepresentable",
            })?;
        let mut points = Vec::with_capacity(3 * mesh.tris.len());
        for (element, tri) in mesh.tris.iter().enumerate() {
            let x = tri.map(|node| mesh.nodes[node].0);
            let y = tri.map(|node| mesh.nodes[node].1);
            let (curvature, twice_area) =
                fs_plate::dkt_curvature_quadrature(&x, &y, element).map_err(map_plate)?;
            let weight = relaxation[0] * (twice_area / 6.0);
            if !(twice_area > 0.0 && weight.is_finite() && weight >= 0.0) {
                return Err(AcousticRealizeError::InvalidDescription {
                    what: "thermoelastic DKT projection requires positive area and finite nonnegative weight",
                });
            }
            for b in curvature {
                points.push(ThermalPoint {
                    dofs: std::array::from_fn(|j| 3 * tri[j / 3] + j % 3),
                    trace: std::array::from_fn(|j| b[0][j] + b[1][j]),
                    weight,
                });
            }
        }
        Ok(Self {
            first_relaxation_time_s: law.relaxation_time(thickness),
            points,
        })
    }

    pub(super) fn modal_zeta(
        &self,
        total: &PlateModel,
        phi: &[f64],
        omega: f64,
        total_energy: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let mut relaxation_energy = 0.0;
        for point in &self.points {
            let mut trace = 0.0;
            for (&dof, &coefficient) in point.dofs.iter().zip(&point.trace) {
                if let Some(reduced) = total.dof_map[dof] {
                    trace += coefficient * phi[reduced];
                }
            }
            relaxation_energy += point.weight * trace * trace;
        }
        let fraction = super::bending_energy_fraction(relaxation_energy, total_energy)?;
        let zeta = 0.5 * fraction * slab_relaxation_shape(omega, self.first_relaxation_time_s);
        if !(zeta.is_finite() && zeta >= 0.0) {
            return Err(AcousticRealizeError::InvalidDescription {
                what: "plate thermoelastic modal loss is unrepresentable",
            });
        }
        Ok(zeta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_sparse::Coo;

    fn thermal() -> ThermoelasticZener {
        ThermoelasticZener {
            e: 70e9,
            alpha_t: 20e-6,
            t0: 300.0,
            rho: 5000.0,
            cp: 600.0,
            conductivity: 100.0,
        }
    }

    #[test]
    fn dkt_saddle_and_twist_do_not_create_thermoelastic_loss() {
        let law = thermal();
        let thickness: f64 = 0.002;
        let poisson: f64 = 0.3;
        let mesh = PlateMesh {
            nodes: vec![(0.1, -0.2), (1.3, 0.1), (0.4, 0.9)],
            tris: vec![[0, 1, 2]],
        };
        let scale = law.e * thickness.powi(3) / (12.0 * (1.0 - poisson * poisson));
        let d = [
            scale,
            poisson * scale,
            0.0,
            poisson * scale,
            scale,
            0.0,
            0.0,
            0.0,
            0.5 * (1.0 - poisson) * scale,
        ];
        let (stiffness, _) =
            fs_plate::dkt_stiffness(&[0.1, 1.3, 0.4], &[-0.2, 0.1, 0.9], &d, 0).unwrap();
        let mut coo = Coo::new(9, 9);
        for i in 0..9 {
            for j in 0..9 {
                coo.push(i, j, stiffness[9 * i + j]);
            }
        }
        let model = PlateModel {
            k: coo.assemble(),
            m: Coo::new(9, 9).assemble(),
            dof_map: (0..9).map(Some).collect(),
            free: 9,
        };
        let projection =
            ThermoelasticBendingProjection::new(&mesh, law, thickness, poisson).unwrap();
        let omega = 1.0 / law.relaxation_time(thickness);
        for [kx, ky, kxy] in [
            [1.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [1.0, -poisson, 0.0],
        ] {
            let phi: Vec<f64> = mesh
                .nodes
                .iter()
                .flat_map(|&(x, y)| {
                    [
                        0.5 * (kx * x * x + ky * y * y) + kxy * x * y,
                        kx * x + kxy * y,
                        ky * y + kxy * x,
                    ]
                })
                .collect();
            let mut kp = vec![0.0; phi.len()];
            model.k.spmv(&phi, &mut kp);
            let energy: f64 = phi.iter().zip(&kp).map(|(p, k)| p * k).sum();
            assert!(energy > 0.0);
            let zeta = projection.modal_zeta(&model, &phi, omega, energy).unwrap();
            if kx + ky == 0.0 {
                assert!(zeta < 1e-28, "trace-free bending created loss {zeta}");
            } else {
                let participation = if ky == 0.0 {
                    (1.0 + poisson) / (1.0 - poisson)
                } else {
                    1.0
                };
                let expected = 0.5 * law.slab_loss_factor(omega, thickness) * participation;
                assert!((zeta / expected - 1.0).abs() < 1e-12);
            }
            let doubled_total = projection
                .modal_zeta(&model, &phi, omega, 2.0 * energy)
                .unwrap();
            assert_eq!(
                doubled_total,
                0.5 * zeta,
                "conservative storage must dilute only the loss ratio"
            );
        }
    }
}
