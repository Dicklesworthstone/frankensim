# Fleet reference baseline stores (bead c40j)

Governed machine-axis baselines for citable roofline gates (dfh3
design): each JSONL row is a `promote_baseline` record — ≥3 mutually
agreeing quiet probes, named operator, justification, promotion day,
age policy. Gates consume via `FRANKENSIM_BASELINE_STORE=<file>` +
`FRANKENSIM_FIRMWARE_ID=<id>`, or `roofline --baseline <file>`.

Bootstrap/update a machine (re-promotion is the only update path):

    roofline promote --store perf-baselines/<machine>.jsonl \
      --firmware "<os-kernel-id>" --operator "<who>" \
      --justification "<why>" [--probes 3] [--age-days 90]

A loaded host REFUSES promotion (drift bands) — measure quiet.

TRUST BOUNDARY: these stores are operator-trusted and tamper-evident
(content-hashed), NOT independently verified — promotion-authority
signatures are bead fz2.7's layer. Do not present gate results as
third-party-verifiable until that lands.
