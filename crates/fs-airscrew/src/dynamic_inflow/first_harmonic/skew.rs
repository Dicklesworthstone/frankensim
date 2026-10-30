//! Pitt--Peters first-harmonic transport for an oblique incident velocity.
//!
//! Chen, NASA TM-102219, equations (13b), (14), (17), (18), and the wind-axis
//! qualification after (16), is the coefficient source:
//! <https://ntrs.nasa.gov/api/citations/19900006622/downloads/19900006622.pdf>.
//! Use total, not perturbation, thrust and the uncorrected mean apparent mass.
//!
//! Let `a=15 pi/64`, `U=hypot(V_t,V_a+u0)`, `c=(V_a+u0)/U`,
//! `Vm=U+c*u0`, and `X=tan(chi/2)=V_t/(U+V_a+u0)`. In downstream/crosswind
//! coordinates the dimensional gain, after converting source moment signs to
//! the native power-conjugate `[T,-My/R,Mx/R]`, is
//! `H=[[1/(2U),-aX/Vm,0],[aX/U,4c/(Vm(1+c)),0],[0,0,4/(Vm(1+c))]]`.
//! The transport operator is `rho A H^-1`; its unequal off-diagonals must
//! not be symmetrized or transposed. A downstream harmonic means greater
//! induction at the trailing edge, opposite rotor velocity relative to air.
//!
//! Only coefficients rotate into the fixed disk reference. Stored induction
//! never reorients with a changing wind direction, including passage through
//! zero transverse speed. The original tangent-frame rate remains separate.
//! The normal-throughflow branch, nonrigid-wake choice, finite-state spatial
//! truncation and all original section limitations remain in force. This is
//! not a far-wake geometry or full fluid/body momentum/passivity theorem.

use super::{
    FirstHarmonicInflow, HarmonicInduction, HarmonicInflowRate, HarmonicInflowStep, Refusal, dot,
};

/// Maximum coefficient iterations in the small physical midpoint problem.
pub const MAX_SKEW_ITERATIONS: usize = 64;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SkewInflowControls {
    pub maximum_iterations: usize,
    /// Absolute midpoint-equation residual, in m/s.
    pub absolute_tolerance_mps: f64,
    pub relative_tolerance: f64,
}

impl Default for SkewInflowControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 32,
            absolute_tolerance_mps: 1e-12,
            relative_tolerance: 1e-11,
        }
    }
}

impl SkewInflowControls {
    /// Refuse unusable or unbounded controls before any iteration.
    pub fn admit(self) -> Result<(), Refusal> {
        if !(1..=MAX_SKEW_ITERATIONS).contains(&self.maximum_iterations)
            || !self.absolute_tolerance_mps.is_finite()
            || self.absolute_tolerance_mps <= 0.0
            || !self.relative_tolerance.is_finite()
            || !(0.0..=1e-3).contains(&self.relative_tolerance)
        {
            return Err(invalid("invalid skew midpoint controls"));
        }
        Ok(())
    }
}

/// Evaluated transport in the ORIGINAL fixed disk modal coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SkewTransport {
    /// `rho A H^-1`, units N s/m; generally nonsymmetric.
    pub operator_ns_per_m: [[f64; 3]; 3],
    /// Positive downstream direction in the disk's two reference coordinates.
    /// No preferred direction exists at exactly zero transverse incident speed.
    pub downstream: Option<[f64; 2]>,
    pub total_speed_mps: f64,
    pub mass_flow_speed_mps: f64,
    pub cosine_skew: f64,
    pub tangent_half_skew: f64,
    wind_operator: [[f64; 3]; 3],
}

impl SkewTransport {
    /// Generalized transport force in N, without moving-reference terms.
    #[must_use]
    pub fn force_n(self, state: HarmonicInduction) -> [f64; 3] {
        self.operator_ns_per_m.map(|row| dot(row, state.modes()))
    }

    fn power_terms(self, state: HarmonicInduction) -> [f64; 2] {
        let d = self.downstream.unwrap_or([1.0, 0.0]);
        let h = state.first_mps;
        let wind = [
            state.mean_mps,
            d[0] * h[0] + d[1] * h[1],
            -d[1] * h[0] + d[0] * h[1],
        ];
        [
            (self.wind_operator[0][0] * wind[0] * wind[0]
                + self.wind_operator[1][1] * wind[1] * wind[1])
                + self.wind_operator[2][2] * wind[2] * wind[2],
            (self.wind_operator[0][1] + self.wind_operator[1][0]) * wind[0] * wind[1],
        ]
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SkewInflowRate {
    pub modal: HarmonicInflowRate,
    pub transport: SkewTransport,
}

/// Accepted nonlinear midpoint solve and its actual modal balance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SkewInflowStep {
    pub modal: HarmonicInflowStep,
    pub transport: SkewTransport,
    /// Zero for the original analytic normal-flow path. Otherwise at most
    /// the admitted cap; transport is evaluated once initially and once per iteration.
    pub iterations: usize,
    pub midpoint_residual_mps: [f64; 3],
    /// Diagonal transport work evaluated in physical wind coordinates, J.
    pub diagonal_transport_work_j: f64,
    /// Signed mean/harmonic cross work, J. This is not friction or heat.
    pub cross_transport_work_j: f64,
}

fn invalid(message: &str) -> Refusal {
    crate::refuse(
        "skew-inflow-invalid",
        message.into(),
        "supply finite normal-throughflow modes, incident velocity and bounded midpoint controls",
    )
}

impl FirstHarmonicInflow {
    /// Evaluate the full source transport at the CURRENT stored mean induction.
    /// `incident_mps` is `[axial, radial-reference, azimuth-reference]` material
    /// velocity relative to incident air, excluding self-induced flow.
    pub fn skew_transport(
        self,
        mean_mps: f64,
        incident_mps: [f64; 3],
    ) -> Result<SkewTransport, Refusal> {
        self.uniform.admit_flow(incident_mps[0], mean_mps)?;
        if !incident_mps.into_iter().all(f64::is_finite) {
            return Err(invalid("nonfinite incident velocity"));
        }
        let transverse = incident_mps[1].hypot(incident_mps[2]);
        let normal = incident_mps[0] + mean_mps;
        let total = transverse.hypot(normal);
        let rho_a = self.uniform.density_kg_m3 * self.uniform.area_m2;
        let (wind, downstream, cosine, half, mass_flow) = if transverse == 0.0 {
            let harmonic = self.damping(incident_mps[0], mean_mps);
            (
                [
                    [2.0 * rho_a * normal, 0.0, 0.0],
                    [0.0, harmonic, 0.0],
                    [0.0, 0.0, harmonic],
                ],
                None,
                1.0,
                0.0,
                incident_mps[0] + 2.0 * mean_mps,
            )
        } else {
            let cosine = normal / total;
            let half = (transverse / total) / (1.0 + cosine);
            let mass_flow = total + cosine * mean_mps;
            let a = 15.0 * core::f64::consts::PI / 64.0;
            let determinant = 2.0 * cosine / (1.0 + cosine) + a * a * half * half;
            if ![transverse, total, cosine, half, mass_flow, determinant]
                .into_iter()
                .all(f64::is_finite)
                || total <= 0.0
                || mass_flow <= 0.0
                || determinant <= 0.0
            {
                return Err(invalid("unrepresentable or singular skew gain"));
            }
            let wind = [
                [
                    rho_a * total * 4.0 * cosine / ((1.0 + cosine) * determinant),
                    rho_a * total * a * half / determinant,
                    0.0,
                ],
                [
                    -rho_a * mass_flow * a * half / determinant,
                    rho_a * mass_flow / (2.0 * determinant),
                    0.0,
                ],
                [0.0, 0.0, 0.25 * rho_a * mass_flow * (1.0 + cosine)],
            ];
            (
                wind,
                Some([-incident_mps[1] / transverse, -incident_mps[2] / transverse]),
                cosine,
                half,
                mass_flow,
            )
        };
        let d = downstream.unwrap_or([1.0, 0.0]);
        let mut operator = [[0.0; 3]; 3];
        operator[0][0] = wind[0][0];
        for i in 0..2 {
            operator[0][i + 1] = wind[0][1] * d[i];
            operator[i + 1][0] = wind[1][0] * d[i];
            for j in 0..2 {
                operator[i + 1][j + 1] =
                    (wind[1][1] - wind[2][2]) * d[i] * d[j] + if i == j { wind[2][2] } else { 0.0 };
            }
        }
        if !operator
            .iter()
            .flatten()
            .chain(&[total, mass_flow])
            .all(|x| x.is_finite())
        {
            return Err(invalid("skew transport operator overflow"));
        }
        Ok(SkewTransport {
            operator_ns_per_m: operator,
            downstream,
            total_speed_mps: total,
            mass_flow_speed_mps: mass_flow,
            cosine_skew: cosine,
            tangent_half_skew: half,
            wind_operator: wind,
        })
    }

    pub fn rate_skewed(
        self,
        state: HarmonicInduction,
        incident_mps: [f64; 3],
        force_n: [f64; 3],
        frame_rate_rad_s: f64,
    ) -> Result<SkewInflowRate, Refusal> {
        let normal = self.rate(state, incident_mps[0], force_n, frame_rate_rad_s)?;
        let transport = self.skew_transport(state.mean_mps, incident_mps)?;
        let modal = if transport.downstream.is_none() {
            normal
        } else {
            self.rate_with_transport(state, force_n, normal.frame_force_n, transport)?
        };
        Ok(SkewInflowRate { modal, transport })
    }

    fn rate_with_transport(
        self,
        state: HarmonicInduction,
        force_n: [f64; 3],
        frame_force_n: [f64; 3],
        transport: SkewTransport,
    ) -> Result<HarmonicInflowRate, Refusal> {
        let transport_force_n = transport.force_n(state);
        let acceleration_mps2 = core::array::from_fn(|j| {
            (force_n[j] - transport_force_n[j] - frame_force_n[j]) / self.mass_kg[j]
        });
        if !transport_force_n
            .iter()
            .chain(&acceleration_mps2)
            .all(|x| x.is_finite())
        {
            return Err(invalid("skew force or acceleration overflow"));
        }
        Ok(HarmonicInflowRate {
            acceleration_mps2,
            transport_force_n,
            frame_force_n,
        })
    }

    /// Solve the actual nonlinear implicit-midpoint modal equation. An analytic
    /// two-mode Schur solve is repeated only for its mean-dependent physical
    /// coefficients; work is bounded by the caller's at-most-64 iteration cap.
    /// Convergence is judged on the original three equations, not change of a
    /// guess. Unresolved roots or inadmissible midpoint/endpoints refuse.
    /// This local iteration is not a global uniqueness or trajectory certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_skewed_held(
        self,
        before: HarmonicInduction,
        incident_mps: [f64; 3],
        force_n: [f64; 3],
        frame_rate_rad_s: f64,
        dt_s: f64,
        controls: SkewInflowControls,
    ) -> Result<SkewInflowStep, Refusal> {
        controls.admit()?;
        let initial = self.rate_skewed(before, incident_mps, force_n, frame_rate_rad_s)?;
        if !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(invalid("invalid skew step duration"));
        }
        if initial.transport.downstream.is_none() {
            let modal = self.step_held(before, incident_mps[0], force_n, frame_rate_rad_s, dt_s)?;
            let transport = self.skew_transport(modal.midpoint.mean_mps, incident_mps)?;
            let rate = self.rate(modal.midpoint, incident_mps[0], force_n, frame_rate_rad_s)?;
            let residual = core::array::from_fn(|j| {
                modal.midpoint.modes()[j]
                    - before.modes()[j]
                    - 0.5 * dt_s * rate.acceleration_mps2[j]
            });
            return self.skew_receipt(modal, transport, 0, residual, dt_s);
        }
        let half = 0.5 * dt_s;
        let mut transport = initial.transport;
        for iteration in 1..=controls.maximum_iterations {
            let midpoint =
                self.solve_frozen_midpoint(before, force_n, frame_rate_rad_s, half, transport)?;
            let actual = self.rate_skewed(midpoint, incident_mps, force_n, frame_rate_rad_s)?;
            let residual = core::array::from_fn::<_, 3, _>(|j| {
                midpoint.modes()[j] - before.modes()[j] - half * actual.modal.acceleration_mps2[j]
            });
            let scale = before
                .modes()
                .into_iter()
                .chain(midpoint.modes())
                .map(f64::abs)
                .fold(0.0, f64::max);
            let tolerance = controls.absolute_tolerance_mps + controls.relative_tolerance * scale;
            if !tolerance.is_finite() {
                return Err(invalid("unrepresentable midpoint residual tolerance"));
            }
            if residual
                .into_iter()
                .all(|x| x.is_finite() && x.abs() <= tolerance)
            {
                let modes = core::array::from_fn::<_, 3, _>(|j| {
                    2.0 * midpoint.modes()[j] - before.modes()[j]
                });
                let after = HarmonicInduction {
                    mean_mps: modes[0],
                    first_mps: [modes[1], modes[2]],
                };
                self.admit(after, incident_mps[0])?;
                let modal =
                    self.account_step(before, after, midpoint, force_n, actual.modal, dt_s)?;
                return self.skew_receipt(modal, actual.transport, iteration, residual, dt_s);
            }
            transport = actual.transport;
        }
        Err(crate::refuse(
            "skew-inflow-unresolved",
            format!(
                "midpoint equations did not converge in {} coefficient iterations",
                controls.maximum_iterations
            ),
            "reduce the timestep or increase the admitted iteration budget; no partial wake state was returned",
        ))
    }

    fn solve_frozen_midpoint(
        self,
        before: HarmonicInduction,
        force_n: [f64; 3],
        frame_rate: f64,
        half: f64,
        transport: SkewTransport,
    ) -> Result<HarmonicInduction, Refusal> {
        let mut matrix = transport.operator_ns_per_m;
        matrix[1][2] -= self.mass_kg[1] * frame_rate;
        matrix[2][1] += self.mass_kg[1] * frame_rate;
        for i in 0..3 {
            for j in 0..3 {
                matrix[i][j] *= half / self.mass_kg[i];
            }
            matrix[i][i] += 1.0;
        }
        let rhs = core::array::from_fn::<_, 3, _>(|i| {
            before.modes()[i] + half * force_n[i] / self.mass_kg[i]
        });
        let determinant = matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];
        if !matrix
            .iter()
            .flatten()
            .chain(&rhs)
            .chain(&[determinant])
            .all(|x| x.is_finite())
            || determinant <= 0.0
        {
            return Err(invalid("unrepresentable harmonic midpoint block"));
        }
        let inverse = |v: [f64; 2]| {
            [
                (matrix[2][2] * v[0] - matrix[1][2] * v[1]) / determinant,
                (-matrix[2][1] * v[0] + matrix[1][1] * v[1]) / determinant,
            ]
        };
        let response = inverse([matrix[1][0], matrix[2][0]]);
        let unforced = inverse([rhs[1], rhs[2]]);
        let schur = matrix[0][0] - matrix[0][1] * response[0] - matrix[0][2] * response[1];
        let mean = (rhs[0] - matrix[0][1] * unforced[0] - matrix[0][2] * unforced[1]) / schur;
        let first = core::array::from_fn(|i| unforced[i] - response[i] * mean);
        if ![schur, mean, first[0], first[1]]
            .into_iter()
            .all(f64::is_finite)
            || schur <= 0.0
        {
            return Err(invalid("unrepresentable mean midpoint Schur equation"));
        }
        Ok(HarmonicInduction {
            mean_mps: mean,
            first_mps: first,
        })
    }

    fn skew_receipt(
        self,
        modal: HarmonicInflowStep,
        transport: SkewTransport,
        iterations: usize,
        residual: [f64; 3],
        dt_s: f64,
    ) -> Result<SkewInflowStep, Refusal> {
        let power = transport.power_terms(modal.midpoint);
        let diagonal_transport_work_j = power[0] * dt_s;
        let cross_transport_work_j = power[1] * dt_s;
        if ![diagonal_transport_work_j, cross_transport_work_j]
            .into_iter()
            .chain(residual)
            .all(f64::is_finite)
        {
            return Err(invalid("skew work accounting overflow"));
        }
        Ok(SkewInflowStep {
            modal,
            transport,
            iterations,
            midpoint_residual_mps: residual,
            diagonal_transport_work_j,
            cross_transport_work_j,
        })
    }
}
