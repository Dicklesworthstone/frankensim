use super::*;
use std::sync::{Mutex, atomic::{AtomicBool, AtomicUsize, Ordering}};

#[test]
fn explicit_thread_budget_defaults_to_serial_and_refuses_bad_values() {
    assert_eq!(parse_threads(None).unwrap(), 1);
    for n in [1, 2, 8, MAX_THREADS] {
        assert_eq!(parse_threads(Some(&n.to_string())).unwrap(), n);
    }
    for invalid in ["", "0", "-1", "65", "1.5", "all", " 8", "999999999999999999999999"] {
        assert!(parse_threads(Some(invalid)).unwrap_err().to_string().contains(THREAD_ENV));
    }
}

#[test]
fn each_frequency_runs_once_and_merges_in_identical_serial_order() {
    for samples in [1, 2, 3, 17, 97] {
        for threads in [1, 2, 3, 8, MAX_THREADS] {
            let seen = Mutex::new(vec![0; samples]);
            let first_finished = AtomicBool::new(false);
            let mut actual = Vec::new();
            run(samples, threads, |index| {
                if index == samples - 1 {
                    assert!(!first_finished.swap(true, Ordering::SeqCst));
                } else {
                    assert!(first_finished.load(Ordering::SeqCst));
                }
                seen.lock().unwrap()[index] += 1;
                // An order-sensitive, signed bit pattern including negative zero.
                Ok(-(index as f64))
            }, |index, value| actual.push((index, value.to_bits()))).unwrap();
            let expected: Vec<_> = (0..samples).rev().map(|i| (i, (-(i as f64)).to_bits())).collect();
            assert_eq!(actual, expected);
            assert!(seen.into_inner().unwrap().iter().all(|n| *n == 1));
        }
    }
}

#[test]
fn highest_frequency_refuses_before_any_lower_sample_or_acceptance() {
    for threads in [1, 2, 8] {
        let seen = Mutex::new(Vec::new());
        let error = run(17, threads, |index| -> Result<(), Error> {
            seen.lock().unwrap().push(index);
            Err("underresolved highest frequency".into())
        }, |_, _| panic!("a refused sample cannot be accepted")).unwrap_err();
        assert_eq!(error.to_string(), "underresolved highest frequency");
        assert_eq!(*seen.lock().unwrap(), [16]);
    }
    assert!(run(0, 1, |_| Ok(()), |_, _| {}).is_err());
    assert!(run(2, 0, |_| Ok(()), |_, _| {}).is_err());
}

#[derive(Debug)]
struct LocalError(std::rc::Rc<str>);
impl std::fmt::Display for LocalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { write!(f, "{}", self.0) }
}
impl std::error::Error for LocalError {}

#[test]
fn first_serial_refusal_wins_even_when_a_later_chunk_fails_first() {
    let later_failed = AtomicBool::new(false);
    let visited = Mutex::new(Vec::new());
    // Remaining chunks are [7,6,5,4] and [3,2,1,0]. The first waits for
    // the second's refusal before reaching its own earlier-in-serial refusal.
    let error = run(9, 2, |index| -> Result<usize, Error> {
        visited.lock().unwrap().push(index);
        if index == 7 {
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(10);
            while !later_failed.load(Ordering::SeqCst) {
                assert!(std::time::Instant::now() < deadline, "later worker did not run");
                std::thread::yield_now();
            }
        }
        if index == 3 {
            later_failed.store(true, Ordering::SeqCst);
            return Err(Box::new(LocalError("later frequency".into())));
        }
        if index == 6 { return Err(Box::new(LocalError("earlier frequency".into()))); }
        Ok(index)
    }, |_, _| {}).unwrap_err();
    assert_eq!(error.to_string(), "earlier frequency");
    let mut seen = visited.into_inner().unwrap();
    seen.sort_unstable();
    assert_eq!(seen, [3, 6, 7, 8]);
}

#[test]
fn worker_panics_are_reported_only_after_all_started_workers_are_joined() {
    let completed = AtomicUsize::new(0);
    let error = run(4, 3, |index| {
        if index == 2 { panic!("test worker panic"); }
        if index < 2 { completed.fetch_add(1, Ordering::SeqCst); }
        Ok(index)
    }, |_, _| {}).unwrap_err();
    assert_eq!(error.to_string(), "radiation frequency worker panicked");
    assert_eq!(completed.load(Ordering::SeqCst), 2);
}

#[test]
fn cancellation_after_the_preflight_prevents_lower_frequency_work() {
    for threads in [1, 2, 8] {
        let gate = fs_exec::CancelGate::new_clock_free();
        let work = AtomicUsize::new(0);
        let error = run(17, threads, |_| -> Result<(), Error> {
            if gate.is_requested() { return Err("radiation preparation cancelled".into()); }
            work.fetch_add(1, Ordering::SeqCst);
            gate.request();
            Ok(())
        }, |_, _| {}).unwrap_err();
        assert_eq!(error.to_string(), "radiation preparation cancelled");
        assert_eq!(work.load(Ordering::SeqCst), 1);
    }
}
