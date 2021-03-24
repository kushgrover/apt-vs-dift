package explicit;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntCollection;
import it.unimi.dsi.fastutil.ints.IntIterable;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntConsumer;

/**
 * A class for storing and accessing the predecessor relation of an explicit Model.
 * <p>
 * As Model only provide easy access to successors of states,
 * the predecessor relation is computed and stored for subsequent efficient access.
 * <p>
 * Note: Naturally, if the model changes, the predecessor relation
 * has to be recomputed to remain accurate.
 */
public class PredecessorRelationDense implements PredecessorRelation
{
	private final Model model;
	/**
	 * pre[i] provides the list of predecessors of state with index i.
	 */
	List<IntList> pre;

	/**
	 * Constructor. Computes the predecessor relation for the given model
	 * by considering the successors of each state.
	 *
	 * @param model the Model
	 */
	public PredecessorRelationDense(Model model)
	{
		pre = new ArrayList<>(model.getNumStates());
		this.model = model;
		// construct the (empty) array list for all states
		int n = model.getNumStates();
		for (int state = 0; state < n; state++) {
			pre.add(state, new IntArrayList());
		}
		for (int s = 0; s < n; s++) {
			int state = s;
			model.getSuccessorsIterator(s).forEachRemaining((IntConsumer) successor -> {
				IntCollection predecessors = pre.get(successor);
				assert !predecessors.contains(state);
				predecessors.add(state);
			});
		}
	}

	@Override public Model getModel()
	{
		return model;
	}

	/**
	 * Get an Iterable over the predecessor states of {@code s}.
	 */
	@Override public IntIterable getPre(int s)
	{
		return pre.get(s);
	}

	/**
	 * Get an Iterator over the predecessor states of {@code s}.
	 */
	@Override public IntIterator getPredecessorsIterator(int s)
	{
		return getPre(s).iterator();
	}

}
