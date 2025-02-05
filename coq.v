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
