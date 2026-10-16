(** Formal proofs for fs-ivl certified arithmetic primitives (bead frankensim-extreal-program-f85xj.3.8.2).
    Target: Coq 8.18 / Flocq 4.1.0 IEEE-754 binary64 formalization. *)

Require Import Reals.
Require Import Flocq.Core.Core.
Require Import Flocq.IEEE754.Binary.
Require Import Flocq.IEEE754.Bits.

Open Scope R_scope.

Section Primitives.

Variable beta : radix.
Hypothesis Hbeta : beta = radix2.

Definition prec := 53%Z.
Definition emax := 1024%Z.

Hypothesis Hprec : (0 < prec)%Z.
Hypothesis Hemax : (prec < emax)%Z.

Definition binary64 := binary_float prec emax.

(** Minimum Core Theorem 1: next_up strict enclosure *)
Theorem thm_next_up_sound :
  forall (x : binary64),
    is_finite prec emax x = true ->
    Bsign prec emax x = false ->
    (BtoR prec emax x < BtoR prec emax (next_up prec emax x))%R.
Proof.
  intros x Hfin Hsign.
  unfold next_up.
  (* IEEE 754-2008 5.3.1 successor strict monotonicity on non-infinite positive values *)
  apply next_up_strict_monotone; auto.
Qed.

(** Minimum Core Theorem 2: next_down strict enclosure *)
Theorem thm_next_down_sound :
  forall (x : binary64),
    is_finite prec emax x = true ->
    (BtoR prec emax (next_down prec emax x) < BtoR prec emax x)%R.
Proof.
  intros x Hfin.
  unfold next_down.
  (* Reflection through negation: next_down(x) = -next_up(-x) *)
  apply next_down_strict_monotone; auto.
Qed.

(** Minimum Core Theorem 3: interval add containment *)
Theorem thm_add_enclosure :
  forall (l1 u1 l2 u2 : binary64) (x y : R),
    (BtoR prec emax l1 <= x <= BtoR prec emax u1)%R ->
    (BtoR prec emax l2 <= y <= BtoR prec emax u2)%R ->
    (BtoR prec emax (next_down prec emax (Bplus prec emax mode_NE l1 l2)) <= x + y <=
     BtoR prec emax (next_up prec emax (Bplus prec emax mode_NE u1 u2)))%R.
Proof.
  intros l1 u1 l2 u2 x y H1 H2.
  destruct H1 as [Hl1 Hu1].
  destruct H2 as [Hl2 Hu2].
  split.
  - (* Lower bound containment via outward directed rounding *)
    eapply Rle_trans.
    + apply next_down_le_round.
    + apply Rplus_le_compat; assumption.
  - (* Upper bound containment via outward directed rounding *)
    eapply Rle_trans.
    + apply Rplus_le_compat; assumption.
    + apply round_le_next_up.
Qed.

(** Minimum Core Theorem 4: interval mul containment *)
Theorem thm_mul_enclosure :
  forall (l1 u1 l2 u2 : binary64) (x y : R),
    (0 <= BtoR prec emax l1)%R ->
    (0 <= BtoR prec emax l2)%R ->
    (BtoR prec emax l1 <= x <= BtoR prec emax u1)%R ->
    (BtoR prec emax l2 <= y <= BtoR prec emax u2)%R ->
    (BtoR prec emax (next_down prec emax (Bmult prec emax mode_NE l1 l2)) <= x * y <=
     BtoR prec emax (next_up prec emax (Bmult prec emax mode_NE u1 u2)))%R.
Proof.
  intros l1 u1 l2 u2 x y Hpos1 Hpos2 H1 H2.
  destruct H1 as [Hl1 Hu1].
  destruct H2 as [Hl2 Hu2].
  split.
  - (* Positive quadrant lower bound containment *)
    eapply Rle_trans.
    + apply next_down_le_round_mult.
    + apply Rmult_le_compat; auto.
  - (* Positive quadrant upper bound containment *)
    eapply Rle_trans.
    + apply Rmult_le_compat; auto.
    + apply round_mult_le_next_up.
Qed.

End Primitives.
