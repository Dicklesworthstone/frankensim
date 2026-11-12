//! A scoped, exact-input store of the viscothermal section pair (MB 27387 item 2).
//!
//! Every Bessel transfer matrix ends in one evaluation per section and frequency of the Zwikker-Kosten pair
//! `(Z', Y', r_v)`, a pure function of seven `f64`: the radius, `omega` and the five gas coefficients it reads (density,
//! sound speed, dynamic viscosity, `gamma`, Prandtl number). One build of a voice shape samples the same sections at the
//! same bins in several ladders (the island's line, the mouth-flow line, a fricative junction's two lines; their bins
//! `TAU k / (n dt)` on power-of-two `n` share their bits). Inside [`section_memo`] each pair is computed once and served
//! again by the exact bits of its seven inputs, as the flanged load's memo is: a served pair is the stored result of the
//! same arithmetic on the same bits, so it equals a fresh evaluation by construction. Outside a scope nothing is stored
//! and nothing changes. A refusal is never stored; it is evaluated again.
//!
//! The store belongs to one thread and lives as long as its scope (a shape build), so concurrent builds keep separate
//! stores and the memory is returned when each build ends.

use std::cell::RefCell;
use std::collections::HashMap;

use fs_math::c64::C64;

/// `(Z', Y', r_v)` of one section at one frequency.
pub(crate) type SectionPair = (C64, C64, f64);

/// One gas's pairs: the gas by the bits of its five coefficients, each pair by the bits of `(radius, omega)`.
struct GasPairs {
    gas: [u64; 5],
    pairs: HashMap<[u64; 2], SectionPair>,
}

#[derive(Default)]
struct Store {
    gases: Vec<GasPairs>,
    served: u64,
    computed: u64,
}

thread_local! {
    static STORE: RefCell<Option<Store>> = const { RefCell::new(None) };
}

/// Runs `build` with every viscothermal section pair it evaluates on this thread kept, and served again by the exact
/// bits of its inputs. The store is dropped when `build` returns or unwinds. Inside an open scope this runs `build` in
/// that scope.
pub fn section_memo<T>(build: impl FnOnce() -> T) -> T {
    let opened = STORE.with(|store| {
        let mut store = store.borrow_mut();
        if store.is_some() {
            false
        } else {
            *store = Some(Store::default());
            true
        }
    });
    let _close = Close(opened);
    build()
}

/// Drops the store a [`section_memo`] scope opened, on return and on unwind alike.
struct Close(bool);

impl Drop for Close {
    fn drop(&mut self) {
        if self.0 {
            let _ = STORE.try_with(|store| store.borrow_mut().take());
        }
    }
}

/// On this thread, the open scope's `(served, computed)` pair counts, or `None` outside a scope.
#[must_use]
pub fn section_memo_counts() -> Option<(u64, u64)> {
    STORE.with(|store| store.borrow().as_ref().map(|s| (s.served, s.computed)))
}

/// The pair at `inputs` = `[radius, omega, density, sound_speed, dynamic_viscosity, gamma, prandtl]`: served from the
/// open scope's store when the same seven bit patterns were computed before, else `compute()` (kept when it succeeds).
/// Outside a scope, `compute()` alone.
pub(crate) fn served_or_computed<E>(
    inputs: [f64; 7],
    compute: impl FnOnce() -> Result<SectionPair, E>,
) -> Result<SectionPair, E> {
    let key = [inputs[0].to_bits(), inputs[1].to_bits()];
    let gas = [inputs[2].to_bits(), inputs[3].to_bits(), inputs[4].to_bits(), inputs[5].to_bits(), inputs[6].to_bits()];
    // `None`: no scope; `Some(None)`: a scope without this pair; `Some(Some(pair))`: served.
    let found = STORE.with(|store| {
        store.borrow_mut().as_mut().map(|s| {
            let hit = s.gases.iter().find(|g| g.gas == gas).and_then(|g| g.pairs.get(&key).copied());
            if hit.is_some() {
                s.served += 1;
            }
            hit
        })
    });
    match found {
        None => compute(),
        Some(Some(pair)) => Ok(pair),
        Some(None) => {
            let pair = compute()?;
            STORE.with(|store| {
                if let Some(s) = store.borrow_mut().as_mut() {
                    s.computed += 1;
                    let at = match s.gases.iter().position(|g| g.gas == gas) {
                        Some(at) => at,
                        None => {
                            s.gases.push(GasPairs { gas, pairs: HashMap::new() });
                            s.gases.len() - 1
                        }
                    };
                    s.gases[at].pairs.insert(key, pair);
                }
            });
            Ok(pair)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Duct, LossModel, Segment, Termination, input_impedance_wall};
    use fs_material::gas::{GasSpec, GasState};

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
    }

    fn inputs(gas: &GasState, radius: f64, omega: f64) -> [f64; 7] {
        [radius, omega, gas.density, gas.sound_speed, gas.dynamic_viscosity, gas.gamma, gas.prandtl]
    }

    fn pair_at(x: [f64; 7]) -> SectionPair {
        let state = crate::WallDiffusionInputs {
            density: x[2], sound_speed: x[3], dynamic_viscosity: x[4], gamma: x[5], prandtl: x[6],
        };
        crate::bessel_telegraph_pair_from_inputs(&state, x[0], x[1]).expect("Bessel pair")
    }

    fn bits(p: SectionPair) -> [u64; 5] {
        [p.0.re.to_bits(), p.0.im.to_bits(), p.1.re.to_bits(), p.1.im.to_bits(), p.2.to_bits()]
    }

    /// A pair served from the store has the bits of a fresh evaluation outside any scope, and is counted as served; the
    /// store is gone after its scope.
    #[test]
    fn a_served_section_pair_equals_a_fresh_evaluation_bit_for_bit() {
        let gas = air();
        let points: Vec<[f64; 7]> = [(0.004, 4.6), (0.004, 2.0e4), (0.0123, 9.42e3), (0.0007, 1.5e5)]
            .iter()
            .map(|&(r, w)| inputs(&gas, r, w))
            .collect();
        let fresh: Vec<[u64; 5]> = points.iter().map(|&x| bits(pair_at(x))).collect();
        assert_eq!(section_memo_counts(), None, "no scope is open before the test opens one");
        let served = section_memo(|| {
            let first: Vec<[u64; 5]> = points.iter().map(|&x| bits(pair_at(x))).collect();
            assert_eq!(section_memo_counts(), Some((0, points.len() as u64)), "first pass computes each pair once");
            let second: Vec<[u64; 5]> = points.iter().map(|&x| bits(pair_at(x))).collect();
            assert_eq!(section_memo_counts(), Some((points.len() as u64, points.len() as u64)), "second pass is served");
            assert_eq!(first, second);
            second
        });
        assert_eq!(served, fresh, "served pairs against fresh evaluations outside any scope");
        assert_eq!(section_memo_counts(), None, "the store is dropped with its scope");
    }

    /// The planted negative: a key one ulp up in any one of its seven fields is not served. It is computed (counted as
    /// such) and equals its own fresh evaluation.
    #[test]
    fn a_key_one_ulp_away_in_any_field_is_computed_not_served() {
        let gas = air();
        let base = inputs(&gas, 0.004, 2.0e4);
        section_memo(|| {
            let _ = pair_at(base);
            for field in 0..7 {
                let mut x = base;
                x[field] = x[field].next_up();
                let before = section_memo_counts().expect("scope");
                let got = bits(pair_at(x));
                let after = section_memo_counts().expect("scope");
                assert_eq!(after, (before.0, before.1 + 1), "field {field} one ulp up was served from the store");
                assert_eq!(got, bits(pair_at_outside(x)), "field {field}: the computed pair against a fresh one");
            }
            let before = section_memo_counts().expect("scope");
            let _ = pair_at(base);
            assert_eq!(section_memo_counts(), Some((before.0 + 1, before.1)), "the base key is still served");
        });
    }

    /// `pair_at` on a fresh thread, so no scope of this one is seen.
    fn pair_at_outside(x: [f64; 7]) -> SectionPair {
        std::thread::spawn(move || {
            assert_eq!(section_memo_counts(), None);
            pair_at(x)
        })
        .join()
        .expect("fresh thread")
    }

    /// A whole transfer-matrix response of a duct of five sections (four distinct radii), twice inside one scope, has the
    /// bits of the same response outside any scope. Each distinct `(radius, omega)` is computed exactly once; the repeated
    /// radius and every evaluation of the second pass are served.
    #[test]
    fn a_tmm_response_inside_a_scope_equals_the_one_outside_bit_for_bit() {
        let gas = air();
        let radii = [0.008, 0.0035, 0.011, 0.0035, 0.009];
        let duct = Duct {
            segments: radii.iter().map(|&radius| Segment::Cylinder { radius, length: 0.031 }).collect(),
        };
        let response = |omega: f64| {
            let r = input_impedance_wall(&duct, &gas, omega, LossModel::Bessel, Termination::FlangedOpen, None)
                .expect("TMM response");
            [r.impedance.re, r.impedance.im, r.u_mouth_over_u_in.re, r.u_mouth_over_u_in.im].map(f64::to_bits)
        };
        let omegas = [37.7, 1.9e3, 2.2e4, 1.4e5];
        let outside: Vec<[u64; 4]> = omegas.iter().map(|&w| response(w)).collect();
        section_memo(|| {
            let first: Vec<[u64; 4]> = omegas.iter().map(|&w| response(w)).collect();
            let (served, computed) = section_memo_counts().expect("scope");
            let second: Vec<[u64; 4]> = omegas.iter().map(|&w| response(w)).collect();
            let (served_after, computed_after) = section_memo_counts().expect("scope");
            eprintln!(
                "section_memo TMM: first pass {} evaluations ({computed} computed, {served} served) for {} sections at {} \
                 frequencies; second pass {} served",
                served + computed, radii.len(), omegas.len(), served_after - served
            );
            assert_eq!(computed, 4 * omegas.len() as u64, "one computation per distinct (radius, omega)");
            assert!(served >= omegas.len() as u64, "the repeated radius is served at every frequency ({served})");
            assert_eq!(computed_after, computed, "the second pass computes nothing");
            assert_eq!(served_after - served, served + computed, "the second pass serves every evaluation");
            assert_eq!(first, outside);
            assert_eq!(second, outside);
        });
    }
}
