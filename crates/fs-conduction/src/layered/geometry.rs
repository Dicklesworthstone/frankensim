//! Resistance and source moments for one-dimensional Fourier conduction.
//! All heat rates are divided by one fixed reference area. For A(r)/Aref=a,
//! D=int dr/a, V=int a dr and B=int V(r)/a(r) dr. The exact Kirchhoff drop
//! is q_in*D+s*B; its source weight B/(D*V) is not generally one half.
//! Geometry contains no material law, root finder, storage or source query.

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum RadialKind {
    Cylindrical,
    Spherical,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum ShapeKind {
    Planar,
    Radial {
        kind: RadialKind,
        inner: f64,
        reference: f64,
    },
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct LayerShape {
    kind: ShapeKind,
    thickness: f64,
    volume: f64,
    resistance: f64,
    source_weight: f64,
}

impl LayerShape {
    pub(super) fn planar(thickness: f64) -> Self {
        Self {
            kind: ShapeKind::Planar,
            thickness,
            volume: thickness,
            resistance: thickness,
            source_weight: 0.5,
        }
    }

    pub(super) fn radial(
        kind: RadialKind,
        inner: f64,
        reference: f64,
        thickness: f64,
    ) -> Option<Self> {
        let outer = inner + thickness;
        if ![inner, reference, thickness, outer]
            .into_iter()
            .all(positive)
            || reference > inner
            || outer <= inner
            || !(thickness / inner).is_finite()
        {
            return None;
        }
        let mut shape = Self {
            kind: ShapeKind::Radial {
                kind,
                inner,
                reference,
            },
            thickness,
            volume: 0.0,
            resistance: 0.0,
            source_weight: 0.0,
        };
        shape.volume = shape.volume_length(thickness);
        shape.resistance = shape.resistance_length(thickness);
        shape.source_weight = shape.generation_weight(thickness);
        [
            shape.volume,
            shape.resistance,
            shape.source_weight,
            shape.area_ratio(thickness),
        ]
        .into_iter()
        .all(positive)
        .then_some(shape)
    }

    pub(super) fn thickness(self) -> f64 {
        self.thickness
    }
    pub(super) fn volume(self) -> f64 {
        self.volume
    }
    pub(super) fn resistance(self) -> f64 {
        self.resistance
    }

    /// A(x)/Aref, so q_local=q_reference/area_ratio.
    pub(super) fn area_ratio(self, x: f64) -> f64 {
        match self.kind {
            ShapeKind::Planar => 1.0,
            ShapeKind::Radial {
                kind,
                inner,
                reference,
            } => {
                let ratio = (inner + x) / reference;
                match kind {
                    RadialKind::Cylindrical => ratio,
                    RadialKind::Spherical => ratio * ratio,
                }
            }
        }
    }

    /// Swept physical volume divided by Aref; difference-free thin-shell form.
    pub(super) fn volume_length(self, x: f64) -> f64 {
        match self.kind {
            ShapeKind::Planar => x,
            ShapeKind::Radial {
                kind,
                inner,
                reference,
            } => {
                let t = x / inner;
                let ratio = inner / reference;
                match kind {
                    RadialKind::Cylindrical => ratio * x * (1.0 + 0.5 * t),
                    RadialKind::Spherical => ratio * ratio * x * (1.0 + t + t * t / 3.0),
                }
            }
        }
    }

    /// Aref times the integrated reciprocal area, in metres.
    pub(super) fn resistance_length(self, x: f64) -> f64 {
        match self.kind {
            ShapeKind::Planar => x,
            ShapeKind::Radial {
                kind,
                inner,
                reference,
            } => match kind {
                RadialKind::Cylindrical => reference * (x / inner).ln_1p(),
                RadialKind::Spherical => (reference / inner) * ((reference / (inner + x)) * x),
            },
        }
    }

    fn generation_weight(self, x: f64) -> f64 {
        let ShapeKind::Radial { kind, inner, .. } = self.kind else {
            return 0.5;
        };
        let t = x / inner;
        if t == 0.0 {
            return 0.5;
        }
        match kind {
            RadialKind::Cylindrical => {
                // B=x²*(1/4 + b/2), b=(t-ln(1+t))/t².
                // Twenty bounded series terms avoid cancellation in coatings.
                let log = t.ln_1p();
                let b = if t < 0.125 {
                    let mut b = 1.0 / 21.0;
                    for n in (0..19).rev() {
                        b = 1.0 / (n + 2) as f64 - t * b;
                    }
                    b
                } else {
                    (1.0 - log / t) / t
                };
                (t / log) * ((0.25 + 0.5 * b) / (1.0 + 0.5 * t))
            }
            RadialKind::Spherical => {
                // B=x²*(3+t)/(6*(1+t)); scale the weight for thick shells.
                if t < 1.0 {
                    (3.0 + t) / (2.0 * (t * t + 3.0 * t + 3.0))
                } else {
                    let u = 1.0 / t;
                    u * (1.0 + 3.0 * u) / (2.0 * (1.0 + 3.0 * u + 3.0 * u * u))
                }
            }
        }
    }

    pub(super) fn weighted_flux(self, inlet: f64, generated: f64) -> f64 {
        inlet + generated * self.source_weight
    }

    pub(super) fn potential_drop(self, inlet: f64, generated: f64, x: f64) -> f64 {
        if x == self.thickness {
            return self.resistance * self.weighted_flux(inlet, generated);
        }
        let local_generation = generated * (self.volume_length(x) / self.volume);
        self.resistance_length(x) * (inlet + local_generation * self.generation_weight(x))
    }

    /// Invert the swept-volume fraction without subtracting nearly equal radii.
    pub(super) fn stationary_position(self, fraction: f64) -> f64 {
        let ShapeKind::Radial { kind, inner, .. } = self.kind else {
            return fraction * self.thickness;
        };
        let a = inner / (inner + self.thickness);
        let scaled = match kind {
            RadialKind::Cylindrical => {
                fraction * (1.0 + a) / (a + ((1.0 - fraction) * a * a + fraction).sqrt())
            }
            RadialKind::Spherical => {
                let u = ((1.0 - fraction) * a * a * a + fraction).cbrt();
                fraction * (1.0 + a + a * a) / (u * u + u * a + a * a)
            }
        };
        scaled * self.thickness
    }
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// Regular solid centre. The reference area is the core's outer face.
/// Symmetry fixes Q(0)=0; its singular hollow reciprocal-area integral is
/// never evaluated. Uniform generation gives q(r)=s*r/m, m=2 or3.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct SolidCore {
    kind: RadialKind,
    radius: f64,
    volume: f64,
}

impl SolidCore {
    pub(super) fn new(kind: RadialKind, radius: f64) -> Option<Self> {
        let dimension = match kind {
            RadialKind::Cylindrical => 2.0,
            RadialKind::Spherical => 3.0,
        };
        let volume = radius / dimension;
        (positive(radius) && positive(volume)).then_some(Self {
            kind,
            radius,
            volume,
        })
    }

    pub(super) fn radius(self) -> f64 {
        self.radius
    }
    pub(super) fn volume(self) -> f64 {
        self.volume
    }

    pub(super) fn area_ratio(self, position: f64) -> f64 {
        let ratio = position / self.radius;
        match self.kind {
            RadialKind::Cylindrical => ratio,
            RadialKind::Spherical => ratio * ratio,
        }
    }

    pub(super) fn local_flux(self, generated: f64, position: f64) -> f64 {
        generated * (position / self.radius)
    }

    /// K(Tcentre)-K(T(r))=s*r²/(2*m), in W/m.
    pub(super) fn potential_drop(self, generated: f64, position: f64) -> f64 {
        (0.5 * position) * self.local_flux(generated, position)
    }
}

#[cfg(test)]
#[path = "geometry/solid_tests.rs"]
mod solid_tests;
