Inductive full_tree : Type :=
  | Leaf : full_tree
  | Node : full_tree -> full_tree -> full_tree.

Fixpoint num_leaves (t : full_tree) : nat :=
  match t with
  | Leaf => 1
  | Node l r => num_leaves l + num_leaves r
  end.

Fixpoint sackin_aux (t : full_tree) (d : nat) : nat :=
  match t with
  | Leaf => d
  | Node l r => sackin_aux l (S d) + sackin_aux r (S d)
  end.

Definition sackin (t : full_tree) : nat :=
  sackin_aux t 0.

Require Import Arith.

Require Import Lia.

Require Import Nat.

Lemma sackin_aux_succ: forall (t: full_tree) (d: nat),
  sackin_aux t (S d) = sackin_aux t d + num_leaves t.
Proof.
  induction t as [| l IHl r IHr]; intros d.
  - (* Base case: t = Leaf *)
    simpl.
    symmetry. apply Nat.add_1_r.
  - (* Inductive case: t = Node l r *)
    simpl. rewrite IHl, IHr. lia.
Qed.

Lemma sackin_Node: forall (l r: full_tree),
  sackin (Node l r) = sackin l + sackin r + num_leaves l + num_leaves r.
Proof.
  intros l r.
  unfold sackin.
  simpl.
  assert (H1: sackin_aux l 1 = sackin_aux l 0 + num_leaves l) by apply sackin_aux_succ.
  assert (H2: sackin_aux r 1 = sackin_aux r 0 + num_leaves r) by apply sackin_aux_succ.
  rewrite H1, H2.
  ring.
Qed.

Fixpoint sum_f (n : nat) (f : nat -> nat) : nat :=
  match n with
  | 0 => 0
  | S n' => f n + sum_f n' f
  end.

Definition sum_range (n : nat) (f : nat -> nat) : nat :=
  sum_f (n - 1) (fun i => f (S i)).
