//! Two coherent physical ports sharing one positive film-loss basis.
use super::*;
use crate::plane_wave::layers::illumination::dot;

/// Cycle-mean powers [W/m²] at one physical face. The inward normal points
/// along +z at the first face and -z at the last face. Incoming/outgoing powers
/// alone do not give the net power when the original exterior medium absorbs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingPortPower {
    pub incoming_w_m2: f64,
    pub outgoing_w_m2: f64,
    /// Signed interference between the actual incoming and outgoing waves.
    /// It may have either sign and is never treated as film absorption.
    pub interference_into_stack_w_m2: f64,
    /// Incoming minus outgoing plus interference. Either face may emit power.
    pub net_into_stack_w_m2: f64,
}

/// Simultaneous actual fields at the first and last faces, in global Cartesian
/// coordinates. Full boundary E/H preserve independent stable field numerators.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AbsorbingTwoSidedReflection {
    pub incoming: [InterfaceBoundaryField; 2],
    pub outgoing: [InterfaceBoundaryField; 2],
    pub boundary: [InterfaceBoundaryField; 2],
    pub power: [AbsorbingPortPower; 2],
    /// Independently evaluated nonnegative loss in the finite films only.
    pub absorbed_power_w_m2: f64,
}

struct RightFilm {
    // Right-launched unit fields, oriented on the ORIGINAL global +z axis.
    fields: [LayerFieldKernel; 2],
    // Shared reference and parity with the original left-launched column.
    loss_modes: [[C64; 2]; 2],
}

pub(in crate::plane_wave::layers::absorbing_incidence) struct PreparedTwoSided<'s> {
    source: &'s Prepared,
    admittances: [[C64; 2]; 2],
    // Polarization -> outgoing face -> incoming face, fixed tangential E.
    scattering: [[[C64; 2]; 2]; 2],
    // Face -> polarization -> E/H -> incoming source.
    boundary: [[[[C64; 2]; 2]; 2]; 2],
    right_films: Vec<RightFilm>,
}

pub(in crate::plane_wave::layers::absorbing_incidence) struct TwoSidedIllumination {
    pub(in crate::plane_wave::layers::absorbing_incidence) response: AbsorbingTwoSidedReflection,
    pub(in crate::plane_wave::layers::absorbing_incidence) absorbed: [f64; MAX_LAYERS],
}

impl<'s> PreparedTwoSided<'s> {
    pub(in crate::plane_wave::layers::absorbing_incidence) fn new(
        source: &'s Prepared,
    ) -> Result<Self> {
        // Keep the original resolved constitutive values; no material law is
        // reevaluated and no epsilon is reconstructed from a rounded index.
        let input = Input {
            epsilon_i: source.input.epsilon_t,
            mu_i: source.input.mu_t,
            epsilon_t: source.input.epsilon_i,
            mu_t: source.input.mu_i,
            ..source.input
        };
        let films: Vec<_> = source.films.iter().rev().map(|f| f.source).collect();
        let reverse = Prepared::with_reference(input, &films, source.reference)?;
        let scattering = core::array::from_fn(|c| {
            [
                [source.reflected[c], reverse.transmitted_e[c]],
                [source.transmitted_e[c], reverse.reflected[c]],
            ]
        });
        let entrance = |prepared: &Prepared, c: usize| -> Result<[C64; 2]> {
            if let Some(first) = prepared.films.first() {
                first.fields[c]
                    .field(0.0)
                    .ok_or_else(|| bad("first-face field is unrepresentable"))
            } else {
                Ok([prepared.transmitted_e[c], prepared.transmitted_h[c]])
            }
        };
        let mut boundary = [[[[C64::ZERO; 2]; 2]; 2]; 2];
        for c in 0..2 {
            let left = entrance(source, c)?;
            let right = entrance(&reverse, c)?;
            boundary[0][c] = [
                [left[0], reverse.transmitted_e[c]],
                [left[1], -reverse.transmitted_h[c]],
            ];
            boundary[1][c] = [
                [source.transmitted_e[c], right[0]],
                [source.transmitted_h[c], -right[1]],
            ];
        }
        let mut right_films = Vec::with_capacity(films.len());
        for reversed in reverse.films.iter().rev() {
            let medium = reversed.medium;
            let reorient = |c: usize| -> Result<LayerFieldKernel> {
                // Rebuild from exact stored face traces instead of evaluating
                // at thickness-z, whose rounded subtraction can lose phase.
                let left = reversed.fields[c]
                    .field(reversed.source.thickness)
                    .ok_or_else(|| bad("reverse film face is unrepresentable"))?;
                let right = reversed.fields[c]
                    .field(0.0)
                    .ok_or_else(|| bad("reverse film face is unrepresentable"))?;
                let (a, b) = if c == 0 {
                    (medium.normal * medium.tm_impedance, medium.omega_epsilon)
                } else {
                    (
                        C64::from_re(medium.omega_mu),
                        medium.normal * medium.te_admittance,
                    )
                };
                LayerFieldKernel::new(
                    reversed.source.thickness,
                    medium.normal,
                    Dd::from_f64(medium.normal.re),
                    a,
                    b,
                    [left[0], -left[1]],
                    [right[0], -right[1]],
                )
                .ok_or_else(|| bad("right-launched film field is unrepresentable"))
            };
            right_films.push(RightFilm {
                fields: [reorient(0)?, reorient(1)?],
                // Swapping the two slab faces leaves the even mode unchanged
                // and reverses the odd mode. Keep that phase before squaring.
                loss_modes: reversed.loss_modes.map(|m| [m[0], -m[1]]),
            });
        }
        Ok(Self {
            source,
            admittances: [source.yi, reverse.yi],
            scattering,
            boundary,
            right_films,
        })
    }

    pub(in crate::plane_wave::layers::absorbing_incidence) fn scattering(
        &self,
    ) -> [[[C64; 2]; 2]; 2] {
        self.scattering
    }
    pub(in crate::plane_wave::layers::absorbing_incidence) fn admittances(&self) -> [[C64; 2]; 2] {
        self.admittances
    }

    pub(in crate::plane_wave::layers::absorbing_incidence) fn illuminate(
        &self,
        left: [C64; 2],
        right: [C64; 2],
    ) -> Result<TwoSidedIllumination> {
        let electric = [left, right];
        if electric.iter().flatten().any(|&e| !finite(e)) {
            return Err(bad("finite incident fields at both stack faces required"));
        }
        let outgoing_e: [[C64; 2]; 2] = core::array::from_fn(|face| {
            core::array::from_fn(|c| dot(self.scattering[c][face], [left[c], right[c]]))
        });
        let media = [self.source.incident, self.source.exit];
        let zero = InterfaceBoundaryField {
            electric_v_m: [C64::ZERO; 3],
            magnetic_a_m: [C64::ZERO; 3],
        };
        let mut incoming = [zero; 2];
        let mut outgoing = incoming;
        let mut boundary = incoming;
        let mut powers = [[Dd::ZERO; 3]; 2];
        let mut net = [Dd::ZERO; 2];
        let mut port_power = [AbsorbingPortPower {
            incoming_w_m2: 0.0,
            outgoing_w_m2: 0.0,
            interference_into_stack_w_m2: 0.0,
            net_into_stack_w_m2: 0.0,
        }; 2];
        for face in 0..2 {
            let sign = if face == 0 { 1.0 } else { -1.0 };
            let yi = self.admittances[face];
            let hi = core::array::from_fn(|c| (yi[c] * electric[face][c]).scale(sign));
            let ho = core::array::from_fn(|c| (yi[c] * outgoing_e[face][c]).scale(-sign));
            incoming[face] = media[face].field(self.source.input.kx, electric[face], hi)?;
            outgoing[face] = media[face].field(self.source.input.kx, outgoing_e[face], ho)?;
            let e = core::array::from_fn(|c| dot(self.boundary[face][c][0], [left[c], right[c]]));
            let h = core::array::from_fn(|c| dot(self.boundary[face][c][1], [left[c], right[c]]));
            boundary[face] = media[face].field(self.source.input.kx, e, h)?;
            powers[face] = port_terms(yi, electric[face], outgoing_e[face]);
            let p = powers[face].map(Dd::to_f64);
            net[face] = powers[face][0] - powers[face][1] + powers[face][2];
            let net_power = net[face].to_f64();
            if p.iter().any(|v| !v.is_finite())
                || !net_power.is_finite()
                || p[0] < 0.0
                || p[1] < 0.0
            {
                return Err(bad("two-sided exterior power is unrepresentable"));
            }
            if electric[face].iter().any(|&e| e != C64::ZERO) && p[0] == 0.0 {
                return Err(bad("two-sided incident power underflows"));
            }
            port_power[face] = AbsorbingPortPower {
                incoming_w_m2: p[0],
                outgoing_w_m2: p[1],
                interference_into_stack_w_m2: p[2],
                net_into_stack_w_m2: net_power,
            };
        }
        let weight = Dd::from_f64(self.source.reference[3]) * Dd::from_f64(0.5);
        let mut absorbed = [0.0; MAX_LAYERS];
        let mut heat = Dd::ZERO;
        for (i, (lf, rf)) in self.source.films.iter().zip(&self.right_films).enumerate() {
            let mut loss = Dd::ZERO;
            for c in 0..2 {
                for mode in 0..2 {
                    let amplitude = dot(
                        [lf.loss_modes[c][mode], rf.loss_modes[c][mode]],
                        [left[c], right[c]],
                    );
                    loss = loss + quadratic(weight, amplitude);
                }
            }
            absorbed[i] = loss.to_f64();
            if !absorbed[i].is_finite() || absorbed[i] < 0.0 {
                return Err(bad("two-sided film heat is unrepresentable"));
            }
            heat = heat + loss;
        }
        let total_heat = heat.to_f64();
        if !total_heat.is_finite() || total_heat < 0.0 {
            return Err(bad("two-sided total heat is unrepresentable"));
        }
        let scale = powers
            .iter()
            .flatten()
            .fold(total_heat, |s, &p| s.max(p.to_f64().abs()));
        if scale > 0.0 {
            let scale = Dd::from_f64(scale);
            let defect = (net[0] / scale + net[1] / scale - heat / scale).to_f64();
            if !defect.is_finite()
                || defect.abs() > 1024.0 * f64::EPSILON * (self.source.films.len() + 1) as f64
            {
                return Err(bad(
                    "two-sided Poynting balance with exterior interference failed",
                ));
            }
        }
        Ok(TwoSidedIllumination {
            response: AbsorbingTwoSidedReflection {
                incoming,
                outgoing,
                boundary,
                power: port_power,
                absorbed_power_w_m2: total_heat,
            },
            absorbed,
        })
    }

    pub(in crate::plane_wave::layers::absorbing_incidence) fn field(
        &self,
        left: [C64; 2],
        right: [C64; 2],
        layer: usize,
        depth: f64,
    ) -> Result<AbsorbingLayerPoint> {
        if left.iter().chain(&right).any(|&e| !finite(e)) {
            return Err(bad("finite incident fields at both stack faces required"));
        }
        let lf = self
            .source
            .films
            .get(layer)
            .ok_or_else(|| bad("field layer is outside the original stack"))?;
        let rf = &self.right_films[layer];
        let mut e = [C64::ZERO; 2];
        let mut h = e;
        for c in 0..2 {
            let l = lf.fields[c]
                .field(depth)
                .ok_or_else(|| bad("field depth is outside the film or is unrepresentable"))?;
            let r = rf.fields[c]
                .field(depth)
                .ok_or_else(|| bad("field depth is outside the film or is unrepresentable"))?;
            e[c] = dot([l[0], r[0]], [left[c], right[c]]);
            h[c] = dot([l[1], r[1]], [left[c], right[c]]);
        }
        field_point(lf.medium, self.source.input.kx, e, h)
    }
}

// For an inward coordinate, E=a+b and H=Y(a-b). The cross flux is
// Im(Y) Im(b conj(a)), which cannot be reconstructed from a reflection ratio
// when a=0 or when b also contains transmission from the opposite source.
fn port_terms(y: [C64; 2], a: [C64; 2], b: [C64; 2]) -> [Dd; 3] {
    let mut p = [Dd::ZERO; 3];
    for c in 0..2 {
        let weight = Dd::from_f64(y[c].re) * Dd::from_f64(0.5);
        p[0] = p[0] + quadratic(weight, a[c]);
        p[1] = p[1] + quadratic(weight, b[c]);
        let yi = Dd::from_f64(y[c].im);
        p[2] = p[2] + (yi * Dd::from_f64(b[c].im)) * Dd::from_f64(a[c].re)
            - (yi * Dd::from_f64(b[c].re)) * Dd::from_f64(a[c].im);
    }
    p
}

#[cfg(test)]
mod tests;
