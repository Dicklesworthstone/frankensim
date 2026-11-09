//! Reachability of stationary storage to original physical sinks.
//! This is topology admission only; native factorization and physical residuals
//! still decide whether a particular nonlinear operating point is resolved.

#[derive(Debug, PartialEq)]
pub(crate) enum SupportError<E> {
    Invalid,
    Checkpoint(E),
}

/// Reverse edges point from a sink toward storage whose energy can reach it.
/// Conductive edges appear in both directions; advection keeps its direction.
pub(crate) fn first_unanchored<E>(
    reverse: &[Vec<usize>],
    sinks: &[bool],
    checkpoint: &mut impl FnMut() -> Result<(), E>,
) -> Result<Option<usize>, SupportError<E>> {
    checkpoint().map_err(SupportError::Checkpoint)?;
    let n = reverse.len();
    if sinks.len() != n {
        return Err(SupportError::Invalid);
    }
    let mut reached = sinks.to_vec();
    let mut queue = Vec::with_capacity(n);
    for (i, adjacent) in reverse.iter().enumerate() {
        checkpoint().map_err(SupportError::Checkpoint)?;
        for &j in adjacent {
            checkpoint().map_err(SupportError::Checkpoint)?;
            if j >= n {
                return Err(SupportError::Invalid);
            }
        }
        if sinks[i] {
            queue.push(i);
        }
    }
    while let Some(i) = queue.pop() {
        checkpoint().map_err(SupportError::Checkpoint)?;
        for &j in &reverse[i] {
            checkpoint().map_err(SupportError::Checkpoint)?;
            if !reached[j] {
                reached[j] = true;
                queue.push(j);
            }
        }
    }
    checkpoint().map_err(SupportError::Checkpoint)?;
    Ok(reached.iter().position(|seen| !seen))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn check(reverse: &[Vec<usize>], sinks: &[bool]) -> Option<usize> {
        first_unanchored(reverse, sinks, &mut || Ok::<_, ()>(())).unwrap()
    }
    #[test]
    fn directed_channel_reaches_downstream_outlet_only() {
        let reverse = [vec![], vec![0], vec![1]];
        assert_eq!(check(&reverse, &[false, false, true]), None);
        assert_eq!(check(&reverse, &[true, false, false]), Some(1));
    }
    #[test]
    fn closed_circulation_has_no_external_energy_reference() {
        let reverse = [vec![2], vec![0], vec![1]];
        assert_eq!(check(&reverse, &[false; 3]), Some(0));
        assert_eq!(check(&reverse, &[true, false, false]), None);
    }
    #[test]
    fn mixed_conduction_and_branching_throughflow_share_one_support_check() {
        let reverse = [vec![], vec![0, 2], vec![1], vec![1], vec![1]];
        assert_eq!(check(&reverse, &[false, false, false, true, true]), None);
        assert_eq!(
            check(&reverse, &[false, false, false, true, false]),
            Some(4)
        );
    }
    #[test]
    fn disconnected_components_cannot_borrow_another_bath() {
        assert_eq!(
            check(
                &[vec![1], vec![0], vec![3], vec![2]],
                &[true, false, false, false]
            ),
            Some(2)
        );
        assert_eq!(check(&[], &[]), None);
        assert_eq!(check(&[vec![]], &[true]), None);
    }
    #[test]
    fn invalid_indices_and_dimensions_refuse_before_traversal() {
        assert_eq!(
            first_unanchored(&[vec![]], &[], &mut || Ok::<_, ()>(())),
            Err(SupportError::Invalid)
        );
        assert_eq!(
            first_unanchored(&[vec![1]], &[true], &mut || Ok::<_, ()>(())),
            Err(SupportError::Invalid)
        );
    }
    #[test]
    fn cancellation_at_every_checkpoint_refuses_then_retry_is_repeatable() {
        let reverse = [vec![1], vec![0], vec![1]];
        let sinks = [false, false, true];
        let mut count = 0;
        first_unanchored(&reverse, &sinks, &mut || {
            count += 1;
            Ok::<_, ()>(())
        })
        .unwrap();
        for stop in 1..=count {
            let mut calls = 0;
            assert_eq!(
                first_unanchored(&reverse, &sinks, &mut || {
                    calls += 1;
                    if calls == stop { Err(()) } else { Ok(()) }
                }),
                Err(SupportError::Checkpoint(()))
            );
        }
        assert_eq!(check(&reverse, &sinks), None);
    }
}
