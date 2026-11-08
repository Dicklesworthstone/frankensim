//! Fluid-loaded finite viscoelastic plates, with coherent P/SV conversion.
//!
//! The original solid occupies 0<=z<=h. Both lossless exterior fluids retain
//! their original source. Peak phasors use exp(+i omega t). Incoming pressure
//! amplitudes are referenced to their own physical face and share the same
//! conserved real tangential slowness. Normal velocity and traction agree at
//! both faces; shear traction vanishes, without enforcing tangential no-slip.
//!
//! A six-by-six boundary solve uses inward-decaying waves anchored at opposite
//! faces. For small normal phase, an algebraically equivalent cosine/sinc basis
//! retains the second solution even at EXACT internal critical incidence.
//! Neither growing transfer matrices nor an effective thin-plate material is
//! introduced. Interfering evanescent waves can transmit through a finite plate;
//! see Almeida et al., <https://arxiv.org/abs/2302.08826>, for the lossless limit.
//!
//! Preparation allocates only a bounded LU. Drives, local fields and optional
//! positive constitutive-work quadrature allocate nothing. Quadrature has an
//! explicit work/tolerance allowance and refuses unresolved intervals. A flux
//! deficit is reported separately from this independent bulk/shear work: no
//! heat is assigned to either mathematical interface.
//! Exact undriven guided-mode singularities may refuse. No anisotropy, rough
//! surface, compliant contact, finite aperture or thermal evolution is inferred.

use super::*;

/// Fixed upper allowance for one positive constitutive-work quadrature query.
pub const MAX_VISCOELASTIC_WORK_SEGMENTS: usize = 16_384;

/// Optional depth integration accuracy and work allowance. Each segment uses
/// positive 8- and 16-point Gauss rules; refinement compares their results.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticLayerWorkBudget {
    pub maximum_segments: usize,
    /// Positive relative tolerance, applied separately to bulk and shear work.
    pub relative_tolerance: f64,
    /// Nonnegative absolute tolerance [W/m²], separately for each loss channel.
    pub absolute_tolerance_w_m2: f64,
}

/// Positive original constitutive work integrated on an admitted depth range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticLayerWork {
    pub bulk_absorbed_power_w_m2: f64,
    pub shear_absorbed_power_w_m2: f64,
    /// Absolute 8/16-point differences, not a rigorous interval certificate.
    pub estimated_absolute_error_w_m2: [f64; 2],
    pub segments_used: usize,
}

/// Independent traveling-wave powers and actual coherent face fluxes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticLayerPower {
    /// Incoming powers at z=0 and z=h, respectively.
    pub incident_power_w_m2: [f64; 2],
    pub outgoing_power_w_m2: [f64; 2],
    /// Coherent +z power in the exterior at each physical face.
    pub fluid_face_power_w_m2: [f64; 2],
    /// Coherent +z power in the original solid; absent for a zero-thickness layer.
    pub solid_face_power_w_m2: Option<[f64; 2]>,
    /// Entering minus leaving flux. Its roundoff is retained; use
    /// `absorption_between` for independently integrated weak bulk/shear heat.
    pub net_absorbed_power_w_m2: f64,
    pub boundary_balance_residual_w_m2: f64,
    pub solid_strain_norm_bound: f64,
    pub fluid_relative_density_bounds: [f64; 2],
    pub fluid_particle_mach_bounds: [f64; 2],
}

/// One original homogeneous viscoelastic solid between two original fluids.
pub struct FluidViscoelasticLayer<'a> {
    source: ViscoelasticSlabSource<'a>,
    exterior_sources: [&'a dyn FluidAcoustics; 2],
    fluids: [FluidAcousticState; 2],
    omega: f64,
    thickness: f64,
    tangential: f64,
    cosines: [C64; 2],
    admittance: [C64; 2],
    moduli: [C64; 2],
    basis: Option<[Basis; 2]>,
    scattering: [C64; 4],
    // Four interior coefficients for each of the two incident pressure ports.
    interior: [[C64; 2]; 4],
}

/// An admitted coherent two-face excitation of the same frozen original solid.
pub struct DrivenFluidViscoelasticLayer<'s, 'a> {
    layer: &'s FluidViscoelasticLayer<'a>,
    coefficients: [C64; 4],
    outgoing: [C64; 2],
    fluid_pressure: [C64; 2],
    fluid_velocity: [[C64; 3]; 2],
    faces: Option<[ViscoelasticHalfSpaceField; 2]>,
    power: ViscoelasticLayerPower,
}

#[derive(Clone, Copy)]
struct Basis {
    waves: [Mode; 2],
    even: Mode,
    odd: Mode,
    cosine: C64,
    k0: C64,
    omega: f64,
    thickness: f64,
    regular: bool,
}

fn combine(first: Mode, a: C64, second: Mode, b: C64) -> Mode {
    Mode {
        q: first.q,
        velocity: core::array::from_fn(|i| first.velocity[i] * a + second.velocity[i] * b),
        strain: core::array::from_fn(|i| first.strain[i] * a + second.strain[i] * b),
        trace: first.trace * a + second.trace * b,
        stress: core::array::from_fn(|i| first.stress[i] * a + second.stress[i] * b),
    }
}

fn strain_norm(wave: Mode) -> f64 {
    wave.strain.iter().fold(0.0, |s, e| det::hypot(s, e.abs()))
}

fn minus_i(z: C64) -> C64 {
    C64::new(z.im, -z.re)
}

/// Evaluated only for |z|<=1/4: analytic, bounded and nonsingular at z=0.
fn cosine_sinc(z: C64) -> [C64; 2] {
    let negative_square = -(z * z);
    let mut cosine = C64::ONE;
    let mut sinc = C64::ONE;
    let mut ct = C64::ONE;
    let mut st = C64::ONE;
    for n in 1..=12 {
        ct = (ct * negative_square).scale(1.0 / ((2 * n - 1) * (2 * n)) as f64);
        st = (st * negative_square).scale(1.0 / ((2 * n) * (2 * n + 1)) as f64);
        cosine = cosine + ct;
        sinc = sinc + st;
    }
    [cosine, sinc]
}

impl Basis {
    fn new(
        density: f64,
        moduli: [C64; 2],
        p: f64,
        shear: bool,
        omega: f64,
        thickness: f64,
    ) -> Result<Self, Error> {
        let forward = mode(density, moduli, p, shear)?;
        let mut backward = forward;
        backward.q = -forward.q;
        if shear {
            backward.velocity[0] = -forward.velocity[0];
            for i in 0..3 {
                backward.strain[i] = -forward.strain[i];
                backward.stress[i] = -forward.stress[i];
            }
        } else {
            backward.velocity[2] = -forward.velocity[2];
            backward.strain[5] = -forward.strain[5];
            backward.stress[5] = -forward.stress[5];
        }
        let even = combine(forward, C64::from_re(0.5), backward, C64::from_re(0.5));
        // (forward-backward)/(2*c*q), formed analytically at q=0 as well.
        let mut odd = Mode {
            q: C64::ZERO,
            velocity: [C64::ZERO; 3],
            strain: [C64::ZERO; 6],
            trace: C64::ZERO,
            stress: [C64::ZERO; 6],
        };
        if shear {
            odd.velocity[0] = C64::ONE;
            odd.strain[0] = C64::from_re(-p);
            odd.strain[2] = C64::from_re(p);
        } else {
            odd.velocity[2] = C64::ONE;
            odd.strain[5] = C64::from_re(-core::f64::consts::SQRT_2 * p);
        }
        odd.stress = odd.strain.map(|e| moduli[1].scale(2.0) * e);
        let modulus = if shear {
            moduli[1]
        } else {
            moduli[0] + moduli[1].scale(4.0 / 3.0)
        };
        let speed = divide(modulus.sqrt(), det::sqrt(density));
        let cosine = speed * forward.q;
        let k0 = C64::from_re(omega) / speed;
        let phase = forward.q.scale(omega).scale(thickness);
        if !finite(cosine) || !finite(k0) || !finite(phase) {
            return Err(Error::Invalid(
                "viscoelastic layer phase is unrepresentable",
            ));
        }
        Ok(Self {
            waves: [forward, backward],
            even,
            odd,
            cosine,
            k0,
            omega,
            thickness,
            regular: phase.abs() <= 0.5,
        })
    }

    fn at(self, depth: f64) -> Result<[Mode; 2], Error> {
        if self.regular {
            let offset = depth - 0.5 * self.thickness;
            let phase = self.waves[0].q.scale(self.omega).scale(offset);
            let [cosine, sinc] = cosine_sinc(phase);
            let sine = phase * sinc;
            Ok([
                combine(self.even, cosine, self.odd, minus_i(self.cosine * sine)),
                combine(
                    self.odd,
                    cosine,
                    self.even,
                    minus_i(self.k0.scale(offset) * sinc),
                ),
            ])
        } else {
            let first = propagation(self.waves[0].q, self.omega, depth)?;
            let second = propagation(self.waves[0].q, self.omega, self.thickness - depth)?;
            Ok([
                combine(self.waves[0], first, self.waves[0], C64::ZERO),
                combine(self.waves[1], second, self.waves[1], C64::ZERO),
            ])
        }
    }

    fn strain_bounds(self) -> [f64; 2] {
        if !self.regular {
            return self.waves.map(strain_norm);
        }
        let phase = self.waves[0]
            .q
            .scale(self.omega)
            .scale(0.5 * self.thickness);
        let growth = det::exp(phase.im.abs());
        let even = strain_norm(self.even);
        let odd = strain_norm(self.odd);
        [
            growth * (even + odd * self.cosine.abs() * phase.abs()),
            growth * (odd + even * self.k0.scale(0.5 * self.thickness).abs()),
        ]
    }
}

impl<'a> FluidViscoelasticLayer<'a> {
    /// Prepare a homogeneous solid, retaining the original source and both
    /// exterior fluids. `left_incidence_cosine` is in (0,1]. The transmitted
    /// right fluid may be propagating, at cutoff, or evanescent. A nonpropagating
    /// right port accepts zero incoming excitation only.
    pub fn try_new(
        source: ViscoelasticSlabSource<'a>,
        exterior_sources: [&'a dyn FluidAcoustics; 2],
        angular_frequency_rad_s: f64,
        left_incidence_cosine: f64,
        thickness_m: f64,
        maximum_branches: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        let moduli = source_moduli(source, angular_frequency_rad_s, maximum_branches, cancelled)?;
        if !source.density_kg_m3.is_finite()
            || source.density_kg_m3 <= 0.0
            || moduli
                .iter()
                .any(|m| !finite(*m) || m.re <= 0.0 || m.im < 0.0)
            || !left_incidence_cosine.is_finite()
            || left_incidence_cosine <= 0.0
            || left_incidence_cosine > 1.0
            || !thickness_m.is_finite()
            || thickness_m < 0.0
        {
            return Err(Error::Invalid(
                "invalid viscoelastic layer source, angle or thickness",
            ));
        }
        let left = exterior_sources[0].acoustic_state()?;
        let right = if core::ptr::eq(exterior_sources[0], exterior_sources[1]) {
            left
        } else {
            exterior_sources[1].acoustic_state()?
        };
        let sine = det::sqrt((1.0 - left_incidence_cosine) * (1.0 + left_incidence_cosine));
        let tangential = sine / left.sound_speed_m_s();
        let right_sine = tangential * right.sound_speed_m_s();
        let square = (Dd::ONE - Dd::from_f64(right_sine) * Dd::from_f64(right_sine)).to_f64();
        let right_cosine = if square >= 0.0 {
            C64::from_re(det::sqrt(square))
        } else {
            C64::new(0.0, -det::sqrt(-square))
        };
        if !tangential.is_finite() || (sine > 0.0 && tangential == 0.0) || !finite(right_cosine) {
            return Err(Error::Invalid(
                "viscoelastic layer exterior slowness is unrepresentable",
            ));
        }
        let cosines = [C64::from_re(left_incidence_cosine), right_cosine];
        let fluids = [left, right];
        let admittance = core::array::from_fn(|i| {
            divide(cosines[i], fluids[i].characteristic_impedance_pa_s_m())
        });
        if admittance.iter().any(|y| !finite(*y)) || admittance[0] == C64::ZERO {
            return Err(Error::Invalid(
                "viscoelastic layer exterior admittance is unrepresentable",
            ));
        }
        let basis = if thickness_m == 0.0 {
            None
        } else {
            Some([
                Basis::new(
                    source.density_kg_m3,
                    moduli,
                    tangential,
                    false,
                    angular_frequency_rad_s,
                    thickness_m,
                )?,
                Basis::new(
                    source.density_kg_m3,
                    moduli,
                    tangential,
                    true,
                    angular_frequency_rad_s,
                    thickness_m,
                )?,
            ])
        };
        let mut prepared = Self {
            source,
            exterior_sources,
            fluids,
            omega: angular_frequency_rad_s,
            thickness: thickness_m,
            tangential,
            cosines,
            admittance,
            moduli,
            basis,
            scattering: [C64::ZERO; 4],
            interior: [[C64::ZERO; 2]; 4],
        };
        if basis.is_some() {
            prepared.prepare_boundary(cancelled)?;
        } else {
            let sum = admittance[0] + admittance[1];
            let reflection = (admittance[0] - admittance[1]) / sum;
            prepared.scattering = [
                reflection,
                admittance[1].scale(2.0) / sum,
                admittance[0].scale(2.0) / sum,
                -reflection,
            ];
        }
        if prepared.scattering.iter().any(|s| !finite(*s)) {
            return Err(Error::Invalid(
                "viscoelastic layer scattering is unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(prepared)
    }

    fn traces(&self, depth: f64) -> Result<[Mode; 4], Error> {
        let basis = self
            .basis
            .ok_or(Error::Invalid("an absent layer has no solid field"))?;
        let p = basis[0].at(depth)?;
        let s = basis[1].at(depth)?;
        Ok([p[0], p[1], s[0], s[1]])
    }

    fn prepare_boundary(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
        let faces = [self.traces(0.0)?, self.traces(self.thickness)?];
        let mut matrix = [C64::ZERO; 36];
        let mut inputs = [[C64::ZERO; 6]; 2];
        for side in 0..2 {
            let row = side * 3;
            let sign = if side == 0 { 1.0 } else { -1.0 };
            matrix[6 * row + side] = self.admittance[side].scale(sign);
            matrix[6 * (row + 1) + side] = C64::ONE;
            for j in 0..4 {
                matrix[6 * row + j + 2] = faces[side][j].velocity[2];
                matrix[6 * (row + 1) + j + 2] = faces[side][j].stress[2];
                matrix[6 * (row + 2) + j + 2] = faces[side][j].stress[5];
            }
            inputs[side][row] = self.admittance[side].scale(sign);
            inputs[side][row + 1] = -C64::ONE;
        }
        for row in 0..6 {
            let scale = matrix[6 * row..6 * row + 6]
                .iter()
                .fold(0.0_f64, |a, b| a.max(b.abs()));
            if !scale.is_finite() || scale == 0.0 {
                return Err(Error::Invalid(
                    "viscoelastic layer boundary row is unrepresentable",
                ));
            }
            for value in &mut matrix[6 * row..6 * row + 6] {
                *value = divide(*value, scale);
            }
            for rhs in &mut inputs {
                rhs[row] = divide(rhs[row], scale);
            }
        }
        let factor = lu_complex_with_cancel(&matrix, 6, &mut *cancelled).map_err(solve_error)?;
        for side in 0..2 {
            poll(cancelled)?;
            let mut solved = inputs[side];
            factor
                .solve_with_cancel(&mut solved, &mut *cancelled)
                .map_err(solve_error)?;
            if solved.iter().any(|a| !finite(*a)) {
                return Err(Error::Invalid(
                    "viscoelastic layer solution is unrepresentable",
                ));
            }
            for row in 0..6 {
                let terms: [C64; 6] = core::array::from_fn(|j| matrix[6 * row + j] * solved[j]);
                let residual =
                    (terms.iter().copied().fold(C64::ZERO, |a, b| a + b) - inputs[side][row]).abs();
                let scale = terms.iter().map(|a| a.abs()).sum::<f64>() + inputs[side][row].abs();
                if !residual.is_finite() || (scale > 0.0 && residual / scale > 1e-11) {
                    return Err(Error::Invalid(
                        "viscoelastic layer boundary continuity did not converge",
                    ));
                }
            }
            self.scattering[side] = solved[0];
            self.scattering[2 + side] = solved[1];
            for j in 0..4 {
                self.interior[j][side] = solved[j + 2];
            }
        }
        Ok(())
    }

    pub const fn source(&self) -> ViscoelasticSlabSource<'a> {
        self.source
    }
    pub const fn exterior_sources(&self) -> [&'a dyn FluidAcoustics; 2] {
        self.exterior_sources
    }
    pub const fn exterior_states(&self) -> [FluidAcousticState; 2] {
        self.fluids
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn thickness_m(&self) -> f64 {
        self.thickness
    }
    pub const fn tangential_slowness_s_m(&self) -> f64 {
        self.tangential
    }
    pub const fn exterior_normal_cosines(&self) -> [C64; 2] {
        self.cosines
    }
    pub const fn moduli_pa(&self) -> [C64; 2] {
        self.moduli
    }
    /// Row-major outgoing [left,right] pressure from incoming [left,right].
    /// All amplitudes are referenced to their respective physical face.
    pub const fn pressure_scattering(&self) -> [C64; 4] {
        self.scattering
    }

    /// Reciprocal power-wave normalization when BOTH ports propagate. Cutoff
    /// and evanescent right ports have no positive traveling-power metric.
    pub fn power_normalized_scattering(&self) -> Option<[C64; 4]> {
        if self.admittance.iter().any(|y| y.re <= 0.0 || y.im != 0.0) {
            return None;
        }
        let root = self.admittance.map(|y| det::sqrt(y.re));
        let normalized =
            core::array::from_fn(|i| self.scattering[i].scale(root[i / 2] / root[i % 2]));
        normalized.iter().all(|s| finite(*s)).then_some(normalized)
    }

    pub fn normal_slownesses_s_m(&self) -> Option<[C64; 2]> {
        self.basis.map(|b| b.map(|m| m.waves[0].q))
    }

    /// Admit coherent peak pressures and the two original fluids' small-signal
    /// allowances. The solid strain triangle bound holds over the entire slab.
    pub fn drive<'s>(
        &'s self,
        incident_pressure_pa: [C64; 2],
        limits: [FluidDriveLimits; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenFluidViscoelasticLayer<'s, 'a>, Error> {
        poll(cancelled)?;
        if incident_pressure_pa.iter().any(|p| !finite(*p))
            || (self.admittance[1].re == 0.0 && incident_pressure_pa[1] != C64::ZERO)
        {
            return Err(Error::Invalid(
                "finite incoming propagating fluid pressures are required",
            ));
        }
        let incoming = incident_pressure_pa;
        let outgoing = core::array::from_fn(|i| {
            self.scattering[2 * i] * incoming[0] + self.scattering[2 * i + 1] * incoming[1]
        });
        let coefficients = core::array::from_fn(|i| {
            self.interior[i][0] * incoming[0] + self.interior[i][1] * incoming[1]
        });
        let fluid_pressure = core::array::from_fn(|i| incoming[i] + outgoing[i]);
        let fluid_velocity = core::array::from_fn(|i| {
            [
                divide(
                    fluid_pressure[i].scale(self.tangential),
                    self.fluids[i].density_kg_m3(),
                ),
                C64::ZERO,
                self.admittance[i]
                    * (incoming[i] - outgoing[i]).scale(if i == 0 { 1.0 } else { -1.0 }),
            ]
        });
        let density_bound = core::array::from_fn(|i| {
            (incoming[i].abs() + outgoing[i].abs()) / self.fluids[i].adiabatic_bulk_modulus_pa()
        });
        let mach_bound = core::array::from_fn(|i| {
            let velocity_per_pressure = det::hypot(
                self.tangential / self.fluids[i].density_kg_m3(),
                self.admittance[i].abs(),
            );
            (incoming[i].abs() + outgoing[i].abs()) * velocity_per_pressure
                / self.fluids[i].sound_speed_m_s()
        });
        let strain_bound = self.basis.map_or(0.0, |b| {
            let p = b[0].strain_bounds();
            let s = b[1].strain_bounds();
            [p[0], p[1], s[0], s[1]]
                .iter()
                .zip(coefficients)
                .map(|(n, a)| n * a.abs())
                .sum::<f64>()
        });
        if !strain_bound.is_finite()
            || strain_bound > self.source.maximum_absolute_strain
            || (0..2).any(|i| {
                !density_bound[i].is_finite()
                    || !mach_bound[i].is_finite()
                    || density_bound[i] > limits[i].maximum_relative_density()
                    || mach_bound[i] > limits[i].maximum_particle_mach()
            })
        {
            return Err(Error::OutsideSourceSupport);
        }
        let faces = if self.basis.is_some() {
            Some([
                self.field(coefficients, 0.0)?,
                self.field(coefficients, self.thickness)?,
            ])
        } else {
            None
        };
        let incident_power =
            core::array::from_fn(|i| quadratic(0.5 * self.admittance[i].re, incoming[i]).to_f64());
        let outgoing_power =
            core::array::from_fn(|i| quadratic(0.5 * self.admittance[i].re, outgoing[i]).to_f64());
        let fluid_flux = core::array::from_fn(|i| {
            (product(fluid_pressure[i], fluid_velocity[i][2]) * Dd::from_f64(0.5)).to_f64()
        });
        let solid_flux = faces.map(|f| f.map(|v| v.normal_power_w_m2));
        let flux = solid_flux.unwrap_or(fluid_flux);
        let absorbed = (Dd::from_f64(flux[0]) - Dd::from_f64(flux[1])).to_f64();
        let residual = (0..2)
            .fold(-Dd::from_f64(absorbed), |s, i| {
                s + Dd::from_f64(incident_power[i]) - Dd::from_f64(outgoing_power[i])
            })
            .to_f64();
        let incident_total = incident_power[0] + incident_power[1];
        if outgoing
            .iter()
            .chain(&fluid_pressure)
            .chain(fluid_velocity.iter().flatten())
            .any(|p| !finite(*p))
            || incident_power
                .iter()
                .chain(&outgoing_power)
                .chain(&fluid_flux)
                .chain([absorbed, residual].iter())
                .any(|x| !x.is_finite())
            || (incoming.iter().any(|p| *p != C64::ZERO) && incident_total == 0.0)
            || (incident_total > 0.0
                && (residual.abs() > 1e-9 * incident_total
                    || absorbed < -1e-9 * incident_total
                    || (0..2).any(|i| (flux[i] - fluid_flux[i]).abs() > 1e-9 * incident_total)))
        {
            return Err(Error::Invalid(
                "viscoelastic layer boundary power is unrepresentable or unbalanced",
            ));
        }
        poll(cancelled)?;
        Ok(DrivenFluidViscoelasticLayer {
            layer: self,
            coefficients,
            outgoing,
            fluid_pressure,
            fluid_velocity,
            faces,
            power: ViscoelasticLayerPower {
                incident_power_w_m2: incident_power,
                outgoing_power_w_m2: outgoing_power,
                fluid_face_power_w_m2: fluid_flux,
                solid_face_power_w_m2: solid_flux,
                net_absorbed_power_w_m2: absorbed,
                boundary_balance_residual_w_m2: residual,
                solid_strain_norm_bound: strain_bound,
                fluid_relative_density_bounds: density_bound,
                fluid_particle_mach_bounds: mach_bound,
            },
        })
    }

    fn field(
        &self,
        coefficients: [C64; 4],
        depth: f64,
    ) -> Result<ViscoelasticHalfSpaceField, Error> {
        coherent_field(&self.traces(depth)?, coefficients, self.moduli, self.omega)
    }
}

impl DrivenFluidViscoelasticLayer<'_, '_> {
    pub const fn power(&self) -> ViscoelasticLayerPower {
        self.power
    }
    pub const fn outgoing_pressure_pa(&self) -> [C64; 2] {
        self.outgoing
    }
    pub const fn fluid_face_pressure_pa(&self) -> [C64; 2] {
        self.fluid_pressure
    }
    pub const fn fluid_face_velocity_m_s(&self) -> [[C64; 3]; 2] {
        self.fluid_velocity
    }
    pub const fn solid_face_fields(&self) -> Option<[ViscoelasticHalfSpaceField; 2]> {
        self.faces
    }

    /// Full coherent physical field at x=0 in the original nonempty solid.
    pub fn field(
        &self,
        depth_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticHalfSpaceField, Error> {
        poll(cancelled)?;
        if !depth_m.is_finite() || depth_m < 0.0 || depth_m > self.layer.thickness {
            return Err(Error::Invalid(
                "solid depth must lie inside the original layer",
            ));
        }
        let field = self.layer.field(self.coefficients, depth_m)?;
        poll(cancelled)?;
        Ok(field)
    }

    /// Independently integrate positive original bulk/shear work, including all
    /// standing-wave and P/S interference. Initial segments resolve every mode's
    /// complex phase to at most four radians; 8/16-point differences then drive
    /// bounded refinement. Unresolved oscillatory or opaque ranges refuse with
    /// `Budget`. Neither boundary-power subtraction nor heat clamping is used.
    pub fn absorption_between(
        &self,
        start_m: f64,
        end_m: f64,
        budget: ViscoelasticLayerWorkBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticLayerWork, Error> {
        poll(cancelled)?;
        if !start_m.is_finite()
            || !end_m.is_finite()
            || start_m < 0.0
            || end_m < start_m
            || end_m > self.layer.thickness
            || !budget.relative_tolerance.is_finite()
            || budget.relative_tolerance <= 0.0
            || budget.relative_tolerance >= 1.0
            || !budget.absolute_tolerance_w_m2.is_finite()
            || budget.absolute_tolerance_w_m2 < 0.0
        {
            return Err(Error::Invalid(
                "invalid viscoelastic work interval or tolerance",
            ));
        }
        if budget.maximum_segments == 0 || budget.maximum_segments > MAX_VISCOELASTIC_WORK_SEGMENTS
        {
            return Err(Error::Budget);
        }
        if start_m == end_m
            || self.layer.moduli.iter().all(|m| m.im == 0.0)
            || self.coefficients.iter().all(|a| *a == C64::ZERO)
        {
            return Ok(ViscoelasticLayerWork {
                bulk_absorbed_power_w_m2: 0.0,
                shear_absorbed_power_w_m2: 0.0,
                estimated_absolute_error_w_m2: [0.0; 2],
                segments_used: 0,
            });
        }
        let q = self
            .layer
            .normal_slownesses_s_m()
            .ok_or(Error::Invalid("an absent layer has no interior work"))?;
        let rate = q.iter().fold(0.0_f64, |a, v| a.max(v.abs())) * self.layer.omega;
        let needed = ((end_m - start_m) * rate / 4.0).ceil().max(1.0);
        if !needed.is_finite() || needed > budget.maximum_segments as f64 {
            return Err(Error::Budget);
        }
        let mut segments = needed as usize;
        loop {
            poll(cancelled)?;
            let low = self.integrate(start_m, end_m, segments, &GAUSS_8, cancelled)?;
            let high = self.integrate(start_m, end_m, segments, &GAUSS_16, cancelled)?;
            let errors = core::array::from_fn(|i| (high[i] - low[i]).abs());
            if (0..2).all(|i| {
                errors[i] <= budget.absolute_tolerance_w_m2 + budget.relative_tolerance * high[i]
            }) {
                poll(cancelled)?;
                return Ok(ViscoelasticLayerWork {
                    bulk_absorbed_power_w_m2: high[0],
                    shear_absorbed_power_w_m2: high[1],
                    estimated_absolute_error_w_m2: errors,
                    segments_used: segments,
                });
            }
            segments = segments
                .checked_mul(2)
                .filter(|s| *s <= budget.maximum_segments)
                .ok_or(Error::Budget)?;
        }
    }

    fn integrate(
        &self,
        start: f64,
        end: f64,
        segments: usize,
        rule: &[(f64, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[f64; 2], Error> {
        let mut sum = [Dd::ZERO; 2];
        for segment in 0..segments {
            poll(cancelled)?;
            let left = start + (end - start) * (segment as f64 / segments as f64);
            let right = if segment + 1 == segments {
                end
            } else {
                start + (end - start) * ((segment + 1) as f64 / segments as f64)
            };
            let half = 0.5 * (right - left);
            let middle = left + half;
            if half == 0.0 {
                return Err(Error::Invalid(
                    "viscoelastic work depth resolution is exhausted",
                ));
            }
            for &(node, weight) in rule {
                for sign in [-1.0, 1.0] {
                    let field = self
                        .layer
                        .field(self.coefficients, middle + sign * half * node)?;
                    let heat = [field.bulk_heat_density_w_m3, field.shear_heat_density_w_m3];
                    for i in 0..2 {
                        sum[i] = sum[i]
                            + (Dd::from_f64(heat[i]) * Dd::from_f64(weight)) * Dd::from_f64(half);
                    }
                }
            }
        }
        let result = sum.map(Dd::to_f64);
        if result.iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err(Error::Invalid("viscoelastic layer work is unrepresentable"));
        }
        Ok(result)
    }
}

const GAUSS_8: [(f64, f64); 4] = [
    (0.1834346424956498, 0.3626837833783620),
    (0.5255324099163290, 0.3137066458778873),
    (0.7966664774136267, 0.2223810344533745),
    (0.9602898564975363, 0.1012285362903763),
];
const GAUSS_16: [(f64, f64); 8] = [
    (0.09501250983763744, 0.1894506104550685),
    (0.2816035507792589, 0.1826034150449236),
    (0.4580167776572274, 0.1691565193950025),
    (0.6178762444026438, 0.1495959888165767),
    (0.7554044083550030, 0.1246289712555339),
    (0.8656312023878318, 0.09515851168249278),
    (0.9445750230732326, 0.06225352393864789),
    (0.9894009349916499, 0.02715245941175409),
];
