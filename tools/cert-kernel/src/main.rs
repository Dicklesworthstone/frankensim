//! Deterministic JSONL driver for the independent certified-arithmetic kernel.

#![deny(unsafe_code)]

use frankensim_cert_kernel::crosscheck::run_audit;

fn main() {
    let samples = match parse_samples(std::env::args().skip(1)) {
        Ok(samples) => samples,
        Err(message) => {
            eprintln!("{message}");
            std::process::exit(2);
        }
    };
    let report = match run_audit(samples) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    print!("{}", report.json_lines());
    if !report.is_green() {
        std::process::exit(1);
    }
}

fn parse_samples(mut args: impl Iterator<Item = String>) -> Result<usize, String> {
    let mut samples = 4096_usize;
    while let Some(argument) = args.next() {
        if argument != "--samples" {
            return Err(format!("unknown argument: {argument}"));
        }
        let value = args
            .next()
            .ok_or_else(|| "--samples requires a positive integer".to_owned())?;
        samples = value
            .parse()
            .map_err(|_| "--samples requires a positive integer".to_owned())?;
        if samples == 0 {
            return Err("--samples requires a positive integer".to_owned());
        }
    }
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::parse_samples;

    #[test]
    fn argument_parser_is_fail_closed() {
        assert_eq!(parse_samples([].into_iter()), Ok(4096));
        assert_eq!(
            parse_samples(["--samples".to_owned(), "17".to_owned()].into_iter()),
            Ok(17)
        );
        assert!(parse_samples(["--samples".to_owned()].into_iter()).is_err());
        assert!(parse_samples(["--samples".to_owned(), "0".to_owned()].into_iter()).is_err());
        assert!(parse_samples(["--unknown".to_owned()].into_iter()).is_err());
    }
}
