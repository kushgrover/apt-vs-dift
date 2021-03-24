package explicit;

import com.google.common.collect.Sets;
import common.Time;
import de.tum.in.naturals.set.BoundedNatBitSet;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntStack;
import prism.PrismComponent;

public final class PredecessorRelations
{
	private PredecessorRelations()
	{
	}

	private static boolean isPre(PredecessorRelation predecessors, NatBitSet states) {
		Model model = predecessors.getModel();
		if (states == null) {
			states = NatBitSets.fullSet(model.getNumStates());
		}
		IntIterator stateIterator = states.iterator();
		while (stateIterator.hasNext()) {
			int state = stateIterator.nextInt();
			IntIterator successors = model.getSuccessorsIterator(state);
			while (successors.hasNext()) {
				int successor = successors.nextInt();
				IntIterator preIterator = predecessors.getPredecessorsIterator(successor);
				boolean found = false;
				while (preIterator.hasNext()) {
					int predecessor = preIterator.nextInt();
					if (predecessor == state) {
						found = true;
						break;
					}
				}
				if (!found) {
					return false;
				}
			}
		}
		return true;
	}

	private static boolean isPreStar(Model model, NatBitSet preStar, NatBitSet remain, NatBitSet target, NatBitSet absorbing)
	{
		// Check absorbing
		if (absorbing != null) {
			IntIterator absorbingIterator = absorbing.iterator();
			while (absorbingIterator.hasNext()) {
				int absorbingState = absorbingIterator.nextInt();
				if (preStar.contains(absorbingState)) {
					if (!target.contains(absorbingState)) {
						return false;
					}
				}
			}
		}

		// Check reaching
		NatBitSet preStarComputed = NatBitSets.modifiableCopyOf(target);
		BoundedNatBitSet todo = NatBitSets.boundedSet(model.getNumStates());
		if (remain == null) {
			todo = todo.complement();
		} else {
			todo.or(remain);
		}
		todo.andNot(target);
		if (absorbing != null) {
			todo.andNot(absorbing);
		}

		IntSet removeSet = new IntOpenHashSet();
		while (true) {
			IntIterator todoIterator = todo.iterator();
			while (todoIterator.hasNext()) {
				int todoState = todoIterator.nextInt();
				IntIterator successors = model.getSuccessorsIterator(todoState);
				while (successors.hasNext()) {
					int successorState = successors.nextInt();
					if (preStarComputed.contains(successorState)) {
						removeSet.add(todoState);
						preStarComputed.set(todoState);
						break;
					}
				}
			}
			if (removeSet.isEmpty()) {
				break;
			}
			todo.removeAll(removeSet);
			removeSet.clear();
		}

		if (!preStar.equals(preStarComputed)) {
			System.err.println(String.format("States not occurring:%n%s%nUnexpected states:%n%s", Sets.difference(preStarComputed, preStar),
					Sets.difference(preStar, preStarComputed)));
			return false;
		}

		return true;
	}

	/**
	 * Computes the set Pre*(target) via a DFS, i.e., all states that
	 * are in {@code target} or can reach {@code target} via one or more transitions
	 * from states contained in {@code remain}.
	 * <br/>
	 * If the parameter {@code remain} is {@code null}, then
	 * {@code remain} is considered to include all states in the model.
	 * <br/>
	 * If the parameter {@code absorbing} is not {@code null},
	 * then the states in {@code absorbing} are considered to be absorbing,
	 * i.e., to have a single self-loop, disregarding other outgoing edges.
	 *
	 * @param remain    restriction on the states that may occur
	 *                  on the path to target, {@code null} = all states
	 * @param target    The set of target states
	 * @param absorbing (optional) set of states that should be considered to be absorbing,
	 *                  i.e., their outgoing edges are ignored, {@code null} = no states
	 * @return the set of states Pre*(target)
	 */
	public static NatBitSet calculatePreStar(PredecessorRelation predecessors, NatBitSet remain, NatBitSet target, NatBitSet absorbing)
	{
		assert isPre(predecessors, remain);
		// all target states are in Pre*
		NatBitSet result = NatBitSets.modifiableCopyOf(target);
		// the stack of states whose predecessors have to be considered
		IntStack todo = new IntArrayList(target);

		while (!todo.isEmpty()) {
			int s = todo.popInt();
			assert result.contains(s);

			// for each predecessor in the graph
			predecessors.getPredecessorsIterator(s).forEachRemaining((int p) -> {
				if (absorbing != null && absorbing.contains(p)) {
					// predecessor is absorbing, thus the edge is considered to not exist
					return;
				}

				if (remain == null || remain.contains(p)) {
					// can reach result (and is in remain)
					if (result.add(p)) {
						// add to stack
						todo.push(p);
					}
				}
			});
		}
		return result;
	}

	public static PredecessorRelation forModel(PrismComponent parent, Model model, NatBitSet subset)
	{
		Time time = new Time();

		parent.getLog().print(String.format("Calculating sparse predecessor relation for %s on %d states... ",
				model.getModelType().fullName(), subset.size()));
		parent.getLog().flush();
		PredecessorRelation pre = new PredecessorRelationSparse(model, subset);
		parent.getLog().println(String.format("done (%g seconds)", time.elapsedSeconds()));

		return pre;
	}

	/**
	 * Static constructor to compute the predecessor relation for the given model.
	 * Logs diagnostic information to the log of the given PrismComponent.
	 *
	 * @param parent a PrismComponent (for obtaining the log and settings)
	 * @param model  the model for which the predecessor relation should be computed
	 * @return the predecessor relation
	 **/
	public static PredecessorRelation forModel(PrismComponent parent, Model model)
	{
		Time time = new Time();

		parent.getLog().print(String.format("Calculating predecessor relation for %s ... ", model.getModelType().fullName()));
		parent.getLog().flush();
		PredecessorRelation pre = new PredecessorRelationDense(model);
		parent.getLog().println(String.format("done (%g seconds)", time.elapsedSeconds()));

		return pre;
	}
}
