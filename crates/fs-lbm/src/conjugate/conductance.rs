//! Conservative nonlinear conduction through two solid half cells and a contact.
//!
//! The caller supplies already-admitted positive conductivity tables. Each half
//! cell obeys its own Kirchhoff integral, and the contact obeys `delta_T = R*q`.
//! Only the common heat flux is solved; no interface temperature is clamped to a
//! cell-centred conductivity estimate. Constant extensions are part of the law.

/// Borrowed native conductivity law, including a grid-axis scale.
/// Table temperatures must be strictly increasing and all coefficients positive.
#[derive(Clone, Copy)]
pub(crate) struct ConductivityCurve<'a> {
    reference: f64,
    axis: f64,
    table: &'a [(f64, f64)],
}

impl<'a> ConductivityCurve<'a> {
    pub(crate) const fn constant(conductivity: f64) -> Self {
        Self {
            reference: 1.0,
            axis: conductivity,
            table: &[],
        }
    }

    pub(crate) const fn tabulated(reference: f64, axis: f64, table: &'a [(f64, f64)]) -> Self {
        Self {
            reference,
            axis,
            table,
        }
    }

    fn raw(self, temperature: f64) -> f64 {
        if self.table.is_empty() {
            return 1.0;
        }
        let upper = self.table.partition_point(|&(t, _)| t < temperature);
        if upper == 0 {
            return self.table[0].1;
        }
        if upper == self.table.len() {
            return self.table[upper - 1].1;
        }
        let (ta, ka) = self.table[upper - 1];
        let (tb, kb) = self.table[upper];
        let fraction = (temperature - ta) / (tb - ta);
        (1.0 - fraction) * ka + fraction * kb
    }

    fn at(self, temperature: f64) -> f64 {
        self.axis * (self.raw(temperature) / self.reference)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// A cursor from one cell toward the other cell through the relevant table
/// pieces. Integrals use fractions of the complete temperature difference and
/// conductivity divided by its interval maximum, avoiding large potentials.
struct Leg<'a> {
    curve: ConductivityCurve<'a>,
    target: f64,
    span: f64,
    maximum: f64,
    resistance: f64,
    first: usize,
    end: usize,
    cursor: usize,
    increasing: bool,
    temperature: f64,
    fraction: f64,
    integral: f64,
    start_k: f64,
    finish: f64,
    finish_k: f64,
    width: f64,
    final_piece: bool,
}

impl<'a> Leg<'a> {
    fn new(curve: ConductivityCurve<'a>, cell: f64, target: f64, distance: f64) -> Option<Self> {
        let lo = cell.min(target);
        let hi = cell.max(target);
        let first = curve.table.partition_point(|&(t, _)| t <= lo);
        let end = curve.table.partition_point(|&(t, _)| t < hi);
        let start = curve.raw(cell);
        let mut maximum = start.max(curve.raw(target));
        for &(_, k) in &curve.table[first..end] {
            maximum = maximum.max(k);
        }
        let maximum_k = curve.axis * (maximum / curve.reference);
        if !positive(maximum_k) {
            return None;
        }
        let resistance = distance / maximum_k;
        if !resistance.is_finite() {
            return None;
        }
        let increasing = cell < target;
        let mut leg = Self {
            curve,
            target,
            span: hi - lo,
            maximum,
            resistance,
            first,
            end,
            cursor: if increasing { first } else { end },
            increasing,
            temperature: cell,
            fraction: 0.0,
            integral: 0.0,
            start_k: start / maximum,
            finish: target,
            finish_k: 0.0,
            width: 0.0,
            final_piece: true,
        };
        leg.load_piece()?;
        Some(leg)
    }

    fn load_piece(&mut self) -> Option<()> {
        let knot = if self.increasing {
            (self.cursor < self.end).then_some(self.cursor)
        } else {
            (self.cursor > self.first).then(|| self.cursor - 1)
        };
        let (finish, k) = match knot {
            Some(i) => self.curve.table[i],
            None => (self.target, self.curve.raw(self.target)),
        };
        self.finish = finish;
        self.finish_k = k / self.maximum;
        self.width = (finish - self.temperature).abs() / self.span;
        self.final_piece = knot.is_none();
        if !positive(self.start_k) || !positive(self.finish_k) || !positive(self.width) {
            return None;
        }
        Some(())
    }

    fn piece_integral(&self) -> f64 {
        self.width * (0.5 * self.start_k + 0.5 * self.finish_k)
    }

    fn boundary(&self) -> f64 {
        if self.resistance == 0.0 {
            f64::INFINITY
        } else {
            (self.integral + self.piece_integral()) / self.resistance
        }
    }

    /// Fractional temperature drop and derivative with respect to scaled flux.
    fn evaluate(&self, flux: f64) -> Option<(f64, f64)> {
        if self.resistance == 0.0 {
            return Some((0.0, 0.0));
        }
        let remaining = (flux * self.resistance - self.integral).max(0.0);
        let scale = self.start_k.max(self.finish_k);
        let a = self.start_k / scale;
        let c = self.finish_k / scale;
        let b = 0.5 * (c - a);
        let average = 0.5 * a + 0.5 * c;
        // The cursor brackets the exact root in this piece. Clamping here
        // removes only accumulated endpoint roundoff, not a physical limit.
        let y = ((remaining / scale) / self.width).min(average);
        let discriminant = if b < 0.0 {
            c * c + (-4.0 * b) * (average - y).max(0.0)
        } else {
            a * a + 4.0 * b * y
        };
        let position = if y == 0.0 {
            0.0
        } else {
            (2.0 * y / (a + discriminant.sqrt())).clamp(0.0, 1.0)
        };
        if !position.is_finite() {
            return None;
        }
        let drop = if position == 1.0 && self.final_piece {
            1.0
        } else {
            (self.fraction + self.width * position).min(1.0)
        };
        let k = (1.0 - position) * self.start_k + position * self.finish_k;
        if !positive(k) {
            return None;
        }
        Some((drop, self.resistance / k))
    }

    fn advance(&mut self) -> Option<()> {
        if self.final_piece {
            return None;
        }
        self.fraction += self.width;
        self.integral += self.piece_integral();
        self.temperature = self.finish;
        self.start_k = self.finish_k;
        if self.increasing {
            self.cursor += 1;
        } else {
            self.cursor -= 1;
        }
        self.load_piece()
    }
}

fn residual(a: &Leg<'_>, b: &Leg<'_>, contact: f64, flux: f64) -> Option<(f64, f64)> {
    let (fa, da) = a.evaluate(flux)?;
    let (fb, db) = b.evaluate(flux)?;
    Some((fa + fb + contact * flux - 1.0, da + db + contact))
}

fn root(a: &Leg<'_>, b: &Leg<'_>, contact: f64, mut lo: f64, mut hi: f64) -> Option<f64> {
    let mut x = lo + 0.5 * (hi - lo);
    for iteration in 0..128 {
        let (f, derivative) = residual(a, b, contact, x)?;
        if f == 0.0 {
            return Some(x);
        }
        if f < 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        let middle = lo + 0.5 * (hi - lo);
        if middle == lo || middle == hi || hi - lo <= 8.0 * f64::EPSILON * hi {
            return Some(middle);
        }
        let newton = x - f / derivative;
        // Every other iteration bisects, bounding the entire inner solve even
        // for an ill-conditioned conductivity contrast or a zero Newton step.
        x = if iteration % 2 == 0 && newton > lo && newton < hi && newton != x {
            newton
        } else {
            middle
        };
    }
    None
}

/// Conductance per area between two cell temperatures, W/(m² K).
/// Distances are half-cell lengths and `contact` is area-specific resistance,
/// m² K/W. Zero-length legs are admitted; two zero lengths require a contact.
/// Returns `None` on nonphysical or unrepresentable coefficient arithmetic.
///
/// Knot events from both laws are merged in heat-flux order. Thus each table
/// piece is visited once, with no allocations; only the final piece pair needs
/// a bounded scalar solve. Reversing the entire face returns identical bits.
pub(crate) fn interface_conductance(
    a: ConductivityCurve<'_>,
    ta: f64,
    da: f64,
    b: ConductivityCurve<'_>,
    tb: f64,
    db: f64,
    contact: f64,
) -> Option<f64> {
    if !ta.is_finite()
        || !tb.is_finite()
        || !da.is_finite()
        || da < 0.0
        || !db.is_finite()
        || db < 0.0
        || !contact.is_finite()
        || contact < 0.0
        || !positive(a.axis)
        || !positive(a.reference)
        || !positive(b.axis)
        || !positive(b.reference)
    {
        return None;
    }
    if ta < tb {
        return interface_conductance(b, tb, db, a, ta, da, contact);
    }
    if ta == tb || (a.table.len() <= 1 && b.table.len() <= 1) || (da == 0.0 && db == 0.0) {
        let ka = a.at(ta);
        let kb = b.at(tb);
        if !positive(ka) || !positive(kb) {
            return None;
        }
        let ra = da / ka;
        let rb = db / kb;
        let scale = ra.max(rb).max(contact);
        if !positive(scale) {
            return None;
        }
        let g = (1.0 / (ra / scale + rb / scale + contact / scale)) / scale;
        return positive(g).then_some(g);
    }
    if !(ta - tb).is_finite() {
        return None;
    }
    let mut hot = Leg::new(a, ta, tb, da)?;
    let mut cold = Leg::new(b, tb, ta, db)?;
    let scale = hot.resistance.max(cold.resistance).max(contact);
    if !positive(scale) {
        return None;
    }
    hot.resistance /= scale;
    cold.resistance /= scale;
    let contact = contact / scale;
    let contact_limit = if contact == 0.0 {
        f64::INFINITY
    } else {
        1.0 / contact
    };
    let events = (hot.end - hot.first)
        .checked_add(cold.end - cold.first)?
        .checked_add(2)?;
    let mut low = 0.0;
    for _ in 0..events {
        let hot_limit = hot.boundary();
        let cold_limit = cold.boundary();
        let high = hot_limit.min(cold_limit).min(contact_limit);
        if !positive(high) {
            return None;
        }
        let (f, _) = residual(&hot, &cold, contact, high)?;
        if f >= 0.0 {
            let flux = root(&hot, &cold, contact, low, high)?;
            let g = flux / scale;
            return positive(g).then_some(g);
        }
        if high == contact_limit
            || (high == hot_limit && hot.final_piece)
            || (high == cold_limit && cold.final_piece)
        {
            // At a complete half-cell drop the true residual is nonnegative.
            // Permit only its final floating-point summation uncertainty.
            let g = high / scale;
            return (f >= -16.0 * f64::EPSILON && positive(g)).then_some(g);
        }
        low = high;
        if hot_limit <= high {
            hot.advance()?;
        }
        if cold_limit <= high {
            cold.advance()?;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-100),
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    #[test]
    fn two_material_interface_and_contact_match_independent_integrals() {
        let at = [(300.0, 1.0), (500.0, 3.0)];
        let bt = [(300.0, 4.0), (500.0, 2.0)];
        let a = ConductivityCurve::tabulated(5.0, 2.0, &at);
        let b = ConductivityCurve::tabulated(2.0, 7.0, &bt);
        // Prescribe the common physical heat flux and surface temperatures.
        // Integrals here are independent trapezoid formulas for affine laws.
        for (sa, sb) in [(400.0, 400.0), (430.0, 390.0), (490.0, 310.0)] {
            let q = 75.0;
            let ka = 1.0 + (sa - 300.0) / 100.0;
            let kb = 4.0 - (sb - 300.0) / 100.0;
            let da = (500.0 - sa) * 0.5 * (3.0 + ka) * (2.0 / 5.0) / q;
            let db = (sb - 300.0) * 0.5 * (4.0 + kb) * (7.0 / 2.0) / q;
            let r = (sa - sb) / q;
            let g = interface_conductance(a, 500.0, da, b, 300.0, db, r).unwrap();
            close(g, q / 200.0);
            assert_eq!(
                g.to_bits(),
                interface_conductance(b, 300.0, db, a, 500.0, da, r)
                    .unwrap()
                    .to_bits()
            );
        }
    }

    #[test]
    fn independently_integrated_crossed_knots_and_endpoint_extensions() {
        let at = [(300.0, 1.0), (400.0, 3.0), (500.0, 2.0)];
        let bt = [(280.0, 4.0), (320.0, 2.0), (360.0, 5.0), (420.0, 1.0)];
        let a = ConductivityCurve::tabulated(1.0, 1.0, &at);
        let b = ConductivityCurve::tabulated(1.0, 1.0, &bt);
        // Integral375..650 of A:68.75+250+300. Integral250..340 of B:120+120+55.
        let q = 125.0;
        let g = interface_conductance(a, 650.0, 618.75 / q, b, 250.0, 295.0 / q, 35.0 / q).unwrap();
        close(g, q / 400.0);
    }

    #[test]
    fn zero_length_legs_and_near_equilibrium_retain_the_series_limit() {
        let table = [(300.0, 1.0), (400.0, 3.0), (500.0, 2.0)];
        let a = ConductivityCurve::tabulated(1.0, 1.0, &table);
        let point = ConductivityCurve::constant(1.0);
        // Cell500, surface350, point300: integral350..500=125+250=375.
        let g =
            interface_conductance(a, 500.0, 375.0 / 20.0, point, 300.0, 0.0, 50.0 / 20.0).unwrap();
        close(g, 0.1);
        close(
            interface_conductance(a, 500.0, 0.0, point, 300.0, 0.0, 2.0).unwrap(),
            0.5,
        );
        assert!(interface_conductance(a, 500.0, 0.0, point, 300.0, 0.0, 0.0).is_none());
        close(
            interface_conductance(a, 350.0, 0.3, point, 350.0, 0.7, 0.2).unwrap(),
            1.0 / 1.05,
        );
        let near = f64::from_bits(350.0f64.to_bits() + 8);
        close(
            interface_conductance(a, near, 0.3, point, 350.0, 0.7, 0.2).unwrap(),
            1.0 / 1.05,
        );
    }

    #[test]
    fn constant_laws_and_refusals_do_not_create_heat_paths() {
        let a = ConductivityCurve::constant(2.0);
        let b = ConductivityCurve::constant(7.0);
        close(
            interface_conductance(a, 700.0, 0.2, b, 200.0, 0.7, 0.4).unwrap(),
            1.0 / 0.6,
        );
        assert!(interface_conductance(a, f64::NAN, 0.2, b, 200.0, 0.7, 0.4).is_none());
        assert!(interface_conductance(a, 700.0, -0.2, b, 200.0, 0.7, 0.4).is_none());
        assert!(interface_conductance(a, 700.0, 0.2, b, 200.0, 0.7, -0.4).is_none());
        let invalid = ConductivityCurve::constant(0.0);
        assert!(interface_conductance(invalid, 700.0, 0.2, b, 200.0, 0.7, 0.4).is_none());
    }
}
