//! Shared physical planar stress interpretation for original material owners.
//!
//! D_zz=0 gives zeta_2=zeta_3+mu/3 and
//! gamma_3²=(A/tau_s)²+(B/(sqrt(2)*tau_trace))², where the
//! Guo-corrected numerators are A=gamma_dev2*tau_s and
//! B=sqrt(2/3)*|div u|*tau_trace. Units are admitted by each original owner.

pub(super) fn trace_tau(
    density_kg_m3: f64,
    velocity_scale: f64,
    cell_size_m: f64,
    shear_tau: f64,
    bulk_viscosity_pa_s: f64,
) -> Result<f64, &'static str> {
    let bulk_lattice = ((bulk_viscosity_pa_s / density_kg_m3) / velocity_scale) / cell_size_m;
    let bulk_excess = 3.0 * bulk_lattice;
    let tau = 0.5 + ((shear_tau - 0.5) / 3.0 + bulk_excess);
    if !bulk_viscosity_pa_s.is_finite()
        || bulk_viscosity_pa_s < 0.0
        || !bulk_excess.is_finite()
        || (bulk_viscosity_pa_s > 0.0 && 0.5 + bulk_excess <= 0.5)
        || !shear_tau.is_finite()
        || shear_tau <= 0.5
        || !tau.is_finite()
        || tau <= 0.5
    {
        return Err("physical planar bulk viscosity or trace relaxation is unrepresentable");
    }
    Ok(tau)
}

fn norm(a: f64, b: f64) -> f64 {
    let scale = a.max(b);
    if scale == 0.0 {
        0.0
    } else {
        let x = a / scale;
        let y = b / scale;
        scale * (x * x + y * y).sqrt()
    }
}

/// One scalar constitutive root, with normalized original stress weights.
/// Its equivalent relaxation is auxiliary and must never enter collision.
#[derive(Clone, Copy)]
pub(super) struct PlanarStrain {
    pub(super) numerator: f64,
    weights: [f64; 2],
}

impl PlanarStrain {
    pub(super) fn new([a, b]: [f64; 2]) -> Result<Self, &'static str> {
        if !a.is_finite() || a < 0.0 || !b.is_finite() || b < 0.0 {
            return Err("unrepresentable planar shear/bulk stress");
        }
        // Normalize before forming N or inverse relaxation, avoiding N/g
        // overflow and squared-component underflow at tiny physical stress.
        let scale = a.max(b);
        let (numerator, weights) = if scale == 0.0 {
            (0.0, [1.0, 0.0])
        } else {
            let x = a / scale;
            let y = (b / scale) * core::f64::consts::FRAC_1_SQRT_2;
            let n = norm(x, y);
            (scale * n, [x / n, y / n])
        };
        if !numerator.is_finite() || (scale != 0.0 && numerator == 0.0) {
            return Err("unrepresentable planar strain invariant");
        }
        Ok(Self { numerator, weights })
    }

    pub(super) fn equivalent_tau(
        self,
        shear_tau: f64,
        trace_tau: f64,
    ) -> Result<f64, &'static str> {
        // gamma*dt=N/lambda. Both physical taus come from the SAME query.
        let equivalent = if self.numerator == 0.0 {
            shear_tau
        } else {
            1.0 / norm(self.weights[0] / shear_tau, self.weights[1] / trace_tau)
        };
        if !equivalent.is_finite() || equivalent <= 0.5 {
            return Err("unrepresentable equivalent planar relaxation");
        }
        Ok(equivalent)
    }
}
